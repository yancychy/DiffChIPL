# library(MASS)
# library(tidyverse)
# library(broom)
# library(glmnet)
# library(limma)
# library(edgeR)
# library("devtools")
# library("roxygen2")
# library(GenomicRanges)
# library(Rsamtools)
# library(bamsignals)
# library(readr)
# library(SGSeq)
# devtools::document()
#' Deduce Pearson Correlation Coefficient between M & A Values
#'
#' @param x,y Two numeric vectors representing the signal intensities of two
#'     samples.
#' @return Safely deduced PCC between \code{(x + y)} and \code{(y - x)}.
#' @examples
#' MA.pcc(1:4, 1:4 + c(1, 2, 4, 9))
#'
#' # The robustness.
#' MA.pcc(1, 0)
#' MA.pcc(1:4, 2:5)
MA.pcc <- function(x, y) {
  if (length(x) < 2) return(NA_real_)
  a <- x + y
  b <- y - x
  if (sd(a) == 0 || sd(b) == 0) return(0)
  cor(a, b, method = "pearson")
}

#' Do log2CPM normalization with or without full library size
#'
#' @param rawcount Matrix representing the raw read counts for samples. Row are the peaks, column are the samples
#' @param libsize Numeric vector representing the full library size for all samples.
#'        If it is null, use the sum of read count in each sample as the library size.
#'        By default, user should use full library size to scale the read counts
#'
#' @return normalized read counts \code{(log2(read counts / sum of read counts + 1))} .
#' @examples
#' cpmD = cpmNorm(countAll)
#'
#' # The robustness.
#' cpmD = cpmNorm(countAll, libsize=fd$lsIP)
#' cpmD = cpmNorm(countAll)
cpmNorm <- function(rawcount, libsize=NULL){
  cX = rawcount
  for(i in 1:ncol(cX)){
    if(is.null(libsize)){
      cX[,i] = log2(10^6 * cX[,i] / sum(cX[,i]) +1)
    }else{
      cX[,i] = log2(10^6 * cX[,i] / libsize[i] + 1)
    }
  }
  cX
}

#' Use edgeR to do TMM normalization with or without full library size
#'
#' @param rawcount Matrix representing the raw read counts for samples. Row are the peaks, column are the samples
#' @param lib.size Numeric vector representing the full library size for all samples.
#'        If it is null, use the sum of read count in each sample as the library size.
#' @param group vector or factor giving the experimental group/condition for each sample/library.
#'        If it is null,no control/condition are used.
#'
#' @return normalized read counts \code{(log2(read counts / sum of read counts + 1))} .
#' @examples
#' cpmD = cpmNorm(countAll)
#'
#' # The robustness.
#' cpmD = cpmNorm(countAll, libsize=fd$lsIP)
#' cpmD = cpmNorm(countAll)
tmm_matrix <-  function(counts, lib.size=NULL, group=c(rep(1,3), rep(2,3))) {
  #cli_alert("Applying trimmed mean of M-values (TMM) normalization.")
  if(is.null(group)){
    cds <- edgeR::DGEList(counts, lib.size=lib.size)
  }else{
    cds <- edgeR::DGEList(counts, lib.size=lib.size, group = group)
  }

  cds <- edgeR::calcNormFactors(cds, method = "TMM")
  tmmD <- edgeR::cpm(cds, normalized.lib.sizes = TRUE, log = TRUE)
  tmmD
}


# From edegR source code
#	TMM between two libraries
#	Mark Robinson
edgeRcalcFactorTMM <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3,
                               sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10){
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
  if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if(max(abs(logR)) < 1e-6) return(1)

  #	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS

  #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #	a fix from leonardo ivan almonacid cardenas, since rank() can return
  #	non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

  if(doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  else
    f <- mean(logR[keep], na.rm=TRUE)

  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}

edgeRcalcFactorQuantile <- function (data, lib.size, p=0.75){
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson and Gordon Smyth
  #	Created 16 Aug 2010. Last modified 12 Sep 2020.
  f <- rep_len(1,ncol(data))
  for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
  if(min(f)==0) warning("One or more quantiles are zero")
  f / lib.size
}

# counts = rawCount
# lib.size = fd$lsIP
TMMnormalizaiton <- function(counts, lib.size=NULL ){

  #	Remove all zero rows
  allzero <- rowSums(counts) == 0L
  if(any(allzero)) counts <- counts[!allzero,,drop=FALSE]
  #	Degenerate cases
  if(nrow(counts)==0 || ncol(counts)==1) method="none"

  if(is.null(lib.size)){
    lib.size = colSums(counts)
  }

  f75 <- suppressWarnings(edgeRcalcFactorQuantile(data=counts, lib.size=lib.size, p=0.75))
  if(median(f75) < 1e-20) {
    refColumn <- which.max(colSums(sqrt(x)))
  } else {
    refColumn <- which.min(abs(f75-mean(f75)))
  }
  id = 1:ncol(counts)
  for(j in id){## Q-Q plot for left sample and ref sample
    ref = counts[,refColumn]
    obs = counts[,j]
    #qqplot(ref,obs, main=j) ##Q-Q plot
    #abline(v=quantile(ref, probs = 0.9), col="red")
  }

  sf = rep(1, ncol(counts))
  for(j in id){
    obs <- as.numeric(counts[,j])
    ref <- as.numeric(counts[,refColumn])
    #qqplot(ref,obs) ##Q-Q plot
    sf[j] =  edgeRcalcFactorTMM(obs, ref, libsize.obs=lib.size[j], libsize.ref=lib.size[refColumn],
                                logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  }
  #	Factors should multiple to one
  scaleF <- sf/exp(mean(log(sf)))
  scaleF
  # tmmFactors

  ef =  lib.size * scaleF
  normC = log2(t(t(counts)/ef) *1e6)

  list(data= normC, scaleFactor = scaleF)
}

# geneLen = df1$len
GeTMMnorm <- function(counts, lib.size = NULL, group=group,
                      geneLen){

  # calculate reads per Kbp of gene length (corrected for gene length)
  # gene length is in bp in exppression dataset and converted to Kbp
  rpk <- ((counts*10^3 )/geneLen)
  # comparing groups
  y <- edgeR::DGEList(counts=rpk, group=group, lib.size=NULL)
  # normalize for library size by cacluating scaling factor using TMM (default method)
  y <- edgeR::calcNormFactors(y)
  # normalization factors for each library
  samples = y$samples
  norm_counts <- edgeR::cpm(y,normalized.lib.sizes = TRUE, log = TRUE)
  list(samples=samples, norm_counts=norm_counts)
}



quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#-------------------------------------------
# Median Ratio Normalization function (MRN)
#-------------------------------------------
#' Use MRN with or without full library size
#'
#' @param rawcount Matrix representing the raw read counts for samples. Row are the peaks, column are the samples
#' @param conditions vector or factor giving the experimental group/condition for each sample/library.
#'        If it is null,no control/condition are used.
#' @param lib.size Numeric vector representing the full library size for all samples.
#'        If it is null, use the sum of read count in each sample as the library size.
#' @return normalized read counts \code{(log2(read counts / sum of read counts + 1))} .
#' @examples
#' cpmD = mrnFactors(countAll,conditions=c(rep(1,2), rep(2,2)), libsize=fd$lsIP)
#'
mrnFactors <- function(rawCounts,conditions, libsize=NULL) {
  rawCounts <- as.matrix(rawCounts)
  if(is.null(libsize)){
    totalCounts <- colSums(rawCounts)
  }else{
    totalCounts = libsize
    names(totalCounts) =  colnames(rawCounts)
  }
  normFactors <- totalCounts
  medianRatios <- rep(1,length(conditions))
  names(medianRatios) <- names(normFactors)

  if (sum(conditions==1)>1){
    meanA <- apply(rawCounts[,conditions==1]%*%diag(1/totalCounts[conditions==1]),1,mean)
  }else{
    meanA <- rawCounts[,conditions==1]/totalCounts[conditions==1]
  }

  for (i in 2:max(conditions)) {
    if (sum(conditions==i)>1){
      meanB <- apply(rawCounts[,conditions==i]%*%diag(1/totalCounts[conditions==i]),1,mean)
    }else{
      meanB <- rawCounts[,conditions==i]/totalCounts[conditions==i]
    }
    id = which(meanA > 0 & meanB>0)
    meanANot0 <- meanA[id]
    meanBNot0 <- meanB[id]
    ratios <- meanBNot0/meanANot0
    medianRatios[conditions==i] <- median(ratios)
    normFactors[conditions==i] <- medianRatios[conditions==i]*totalCounts[conditions==i]
  }
  medianRatios <- medianRatios/exp(mean(log(medianRatios)))
  normFactors <- normFactors/exp(mean(log(normFactors)))
  return(list(medianRatios=medianRatios,normFactors=normFactors))
}

#' Get read counts for each sample from the bam files.
#'
#' @param SampleID Matrix representing the raw read counts for samples. Row are the peaks, column are the samples
#' @param lib.size Numeric vector representing the full library size for all samples.
#'        If it is null, use the sum of read count in each sample as the library size.
#' @param group vector or factor giving the experimental group/condition for each sample/library.
#'        If it is null,no control/condition are used.
#'
#' @return normalized read counts \code{(log2(read counts / sum of read counts + 1))} .
#' @examples
#' cpmD = cpmNorm(countAll)
#'
#' # The robustness.
#' cpmD = cpmNorm(countAll, libsize=fd$lsIP)
#' cpmD = cpmNorm(countAll)
getLibsize <- function(SampleID, bamFils, nCore=1){
  DF <- data.frame(sample_name= SampleID, file_bam = bamFils)
  DF$file_bam <- Rsamtools::BamFileList(DF$file_bam)
  DF_complete <- SGSeq::getBamInfo(DF, cores = nCore)
  #DF_complete$lib_size
  gc()
  DF_complete
}


#' Do merge or overlap over input peaks files based their coordinates
#'
#' @param peakFiles path vector of input peak files
#' @param overlapOrMerge logical If TURE, do overlap; else do merge
#' @param len numeric, peak length to filter peaks
#' @return GRanges objects, merged or ovelapped peak coordinates.
#' @examples
#' peaksO = omPeaks(peakID, overlapOrMerge=TRUE)
#'
omPeaks <- function(peakFiles, overlapOrMerge=TRUE, len=2){
  #pth = -log10(padj)
  daL = list()
  for(i in 1:length(peakFiles)){
    df1 <- readr::read_delim(peakFiles[i], delim="\t", col_names=FALSE)
    peak1 <- GenomicRanges::GRanges(df1$X1, IRanges(df1$X2, df1$X3))
    #speak1 <- GenomeInfoDb::keepStandardChromosomes(peak1, pruning.mode="coarse")
    speak1 <- GenomeInfoDb::sortSeqlevels(peak1)
    speak1 <- BiocGenerics::sort(speak1)
    daL[[i]] = speak1
  }
  if(overlapOrMerge){## Compute the overlapped peaks among the replicates
    if(length(daL)==1){
      po = daL[[1]]
    }else{
      po = IRanges::subsetByOverlaps(daL[[1]], daL[[2]])
      if(length(daL)>len){
        for(i in 3:length(peakFiles) ){
          po = IRanges::subsetByOverlaps(po, daL[[i]])
        }
      }
    }
  }else{## Merge the peaks among different conditions
    if(length(daL)==1){
      po = daL[[1]]
    }else{
      po = GenomicRanges::union(daL[[1]], daL[[2]])
      if(length(daL) > 2){
        for(i in 3:length(peakFiles) ){
          po = GenomicRanges::union(po, daL[[i]])
        }
      }
    }
  }
  po
}

#' Do overlap over input peaks files in same condition, then merge the control and treatment peaks
#'
#' @param fd dataframe the file format is  input file to show the sample name,
#'        TF, Control, Treatment, bam file path, and peak files path.
#'
#' @return list
#'         peakL list of GRanges objects, peak coordinates in each peak file
#'         po GRanges objects, which contains overlapped peak coordinate in each peak file
#'         pu GRanges objects, which contains union peak coordinate in each peak file
#' @examples
#' peaksAll = getPeakRegion(fd)
#'
getPeakRegion <- function(fd, removeUN_Random_Chr = FALSE){
  peakL = list()
  fd$SampleID
  for(i in 1:nrow(fd)){
    df1 <- readr::read_delim(fd$Peaks[i], delim="\t", col_names=FALSE)
    peak1 <- GenomicRanges::GRanges(df1$X1, IRanges::IRanges(df1$X2, df1$X3))
    if(removeUN_Random_Chr){
      peak1 <- GenomeInfoDb::keepStandardChromosomes(peak1, pruning.mode="coarse")
    }
    speak1 <- GenomeInfoDb::sortSeqlevels(peak1)
    speak1 <- BiocGenerics::sort(speak1)
    peakL[[i]] = speak1
  }
  names(peakL) = paste0(fd$SampleID, "_", fd$Condition)
  len = length(peakL)

  if(len==1){
    po = peakL[[1]]
    pu = peakL[[1]]
  }else{
    po = IRanges::subsetByOverlaps(peakL[[1]], peakL[[2]])
    pu = GenomicRanges::union(peakL[[1]], peakL[[2]])
    if(len > 2){
      for(i in 3:len){
        po = IRanges::subsetByOverlaps(po, peakL[[i]])
        pu = GenomicRanges::union(pu, peakL[[i]])
      }
    }
  }
  print(length(po))
  print(length(pu))
  list(peakL = peakL, po = po, pu = pu)
}

#' Do overlap over input peaks files in same condition, then merge the control and treatment peaks
#'
#' @param inputF file path, the file format is same to DiffBind input file to show the sample name,
#'        TF, Control, Treatment, bam file path, and peak files path.
#' @param controlName string The control name
#' @param treatmentName string The treatment name
#' @param peakLen integer The length of peak threshold (>2)
#'
#' @return list
#'         pr GRanges objects, common peak coordinates between two conditions
#'         peakL list of GRanges objects, which contains the peak coordinate in each peak file
#' @examples
#' peaksAll = getCommPeaks(inputF, controlName="", treatmentName="",)
#'
getCommPeaks <- function(inputF, controlName="", treatmentName="", peakLen=2){
  fd = read.table(inputF, sep=",", header = TRUE, stringsAsFactors = FALSE)
  fd$Condition

  peakL = list()
  tc = unique(fd$Condition)
  for(i in 1:length(tc)){
    id1 = which(fd$Condition==tc[i])
    peakID = fd$Peaks[id1]
    peskO = omPeaks(peakID)
    peakL[[i]] = peskO
  }
  names(peakL) = tc
  print(tc)

  contrast = c(treatmentName, controlName)
  tr = contrast[1]
  ct = contrast[2]
  po = GenomicRanges::union(peakL[[tr]] ,peakL[[ct]])
  id = which(po$width > peakLen)
  pr = po[id]
  list(pr = pr, peakL=peakL)
}


#' Get read counts from IP and Input bam file with given peak regions. To remove the backgroup noise,
#' compute the scale factor between the IP and input read counts from the non-peak regions
#' by by (All IP reads - reads of IP peaks regions ) / (All Input reads - reads of Input peaks regions)
#' The final read count of peak region are
#' reads of IP peak regions - scale factor * reads of input peak region
#' @param pr GRanges objects, common or merged peak coordinates between two conditions
#' @param bamIP string IP bam file path
#' @param bamInput string Input bam file path, if it is null, do scale based on IP only
#' @param removedup logical default is true to remove duplicates reads.
#' @param fragment integer default is NULL, the fragment length will be determined by SGSeq::getBamInfo().
#'                 User can set it as 0, to avoid shift the reads.
#' @param scaleControl logical default is TRUE, do scale for read counts
#' @param libszie integer vector corresponding to IP and Input libsize.
#' @param removeBackground logical default is TRUE, remove background read count from each sample;
#'        For IP and Input, the IP will
#' @param paired.end  c("ignore", "filter"), refer to bamsignals-methods {bamsignals}
#'        ignore is used for single end reads, filter is used for paired end reads
#'        If paired.end="filter" then only first reads in proper mapped pairs will be
#'        considered (SAMFLAG 66, i.e. in the flag of the read, the bits in the mask 66 must be all ones).
#'
#' @return list
#'         - count integer vector for each peak region
#'         - peakPos GRanges objects, common peak coordinates between two conditions
#'         - cIP integer read counts of IP sample
#'         - cInput integer read counts of input sample
#' @examples
#' readL = getPeakCounts(pr, bamIP, bamInput,removedup=TRUE, fragment = 300,
#'                          scaleControl=TRUE, libsize=c(1e6,1e6))
#'
getPeakCounts <- function(pr, bamIP, bamInput=NULL,removedup=TRUE, fragment = 300,
                          scaleControl=TRUE, libsize=c(1e6,1e6), 
                          removeBackground=TRUE,  removeUN_Random_Chr = FALSE,
                          paired.end= c("ignore", "filter")){

  if(removedup){#remove duplicate reads
    Flag=1024
  }else{
    Flag=-1
  }
  if(scaleControl){
    if(length(libsize)==1){
      scaleFactor=1
    }else{
      scaleFactor = libsize[1]/libsize[2]
    }
  }else{
    scaleFactor = 1
  }
  ### store the middle point of the peak and the end of the peak
  if(1==0){
    pr <- GenomeInfoDb::sortSeqlevels(pr)
    pr <- BiocGenerics::sort(pr)
    
    ps <- GenomeInfoDb::sortSeqlevels(ps)
    ps <- BiocGenerics::sort(ps)
  }

  ### common peak counts
  cIP <- bamsignals::bamCount(bampath = bamIP, gr = pr, ss=FALSE, filteredFlag=Flag,
                              shift= ceiling(fragment/2) )
  gc()
  message(Sys.time(), ": 1. Read count of IP ")
  if(!is.na(bamInput)){

    ### common peak counts
    cCtr = bamsignals::bamCount(bampath = bamInput, gr =pr, ss=FALSE, filteredFlag=Flag,
                                shift= ceiling(fragment/2) )
    gc()
    message(Sys.time(), ": 2. Read count of Input")
    if(removeBackground){
      #ratio = (libsize[1]- sum(cIPRaw) )/(libsize[2] - sum(cCtrRaw)) # scale factor
      ### Here we used NCIS for calculate the scaling factor
      res <- NCIS(chip.data = bamIP, input.data = bamInput, data.type="BAM", 
                  frag.len=fragment, binsize.shift = ceiling(fragment/2),
                  removeUN_Random_Chr = removeUN_Random_Chr)
      
      ratio = res$est
      count = ceiling(cIP - cCtr*ratio)
      message(Sys.time(), ": 3. Get scaling factor by NCIS")
    }else{
      count = cIP
      ratio = rep(0, length(cIP))
      cCtr = rep(0, length(cIP))
    }
    return(list(count=count, peakPos = pr, cIP=cIP, cInput=cCtr, ratio = ratio))
  }else{### Remove background in IP
    if(removeBackground){
      bamFile <- Rsamtools::BamFile(bamIP)
      chrInfo = GenomeInfoDb::seqinfo(bamFile)
      lens = GenomeInfoDb::seqlengths(bamFile)
      lenA = sum(lens)
      peakLen = pr@ranges@width # common peak length
      ratio = peakLen / (lenA - peakLen)#ratio
      backCount = (libsize[1] - sum(cIP)) ## remove reads in raw peaks
      peakBack =  ratio * backCount ## background noise
      ratio = peakBack/sum(cIP)
      count = ceiling(cIP * (1 - ratio))
      #hist(count, breaks = 100)
    }else{
      count = cIP
      ratio = 0
    }
    return(list(count=count, peakPos = pr, cIP=cIP, cInput=NULL, ratio = ratio))
  }
  # summary(count)
  # hist(count,breaks=512)
}


#' Get read counts from IP and Input bam file with given peak regions. To remove the backgroup noise,
#' compute the scale factor between the IP and input read counts from the non-peak regions
#' by by (All IP reads - reads of IP peaks regions ) / (All Input reads - reads of Input peaks regions)
#' The final read count of peak region are
#' reads of IP peak regions - scale factor * reads of input peak region
#' @param pr GRanges objects, common or merged peak coordinates between two conditions
#' @param ps GRanges objects, raw peak coordinates for the IP sample, if it is null
#' @param bamIP string IP bam file path
#' @param bamInput string Input bam file path, if it is null, do scale based on IP only
#' @param removedup logical default is true to remove duplicates reads.
#' @param fragment integer default is NULL, the fragment length will be determined by SGSeq::getBamInfo().
#'                 User can set it as 0, to avoid shift the reads.
#' @param scaleControl logical default is TRUE, do scale for read counts
#' @param libszie integer vector corresponding to IP and Input libsize.
#' @param removeBackground logical default is TRUE, remove background read count from each sample;
#'        For IP and Input, the IP will
#' @param paired.end  c("ignore", "filter"), refer to bamsignals-methods {bamsignals}
#'        ignore is used for single end reads, filter is used for paired end reads
#'        If paired.end="filter" then only first reads in proper mapped pairs will be
#'        considered (SAMFLAG 66, i.e. in the flag of the read, the bits in the mask 66 must be all ones).
#'
#' @return list
#'         - count integer vector for each peak region
#'         - peakPos GRanges objects, common peak coordinates between two conditions
#'         - cIP integer read counts of IP sample
#'         - cInput integer read counts of input sample
#' @examples
#' readL = getPeakCounts(pr, bamIP, bamInput,removedup=TRUE, fragment = 300,
#'                          scaleControl=TRUE, libsize=c(1e6,1e6))
#'
getPeakCounts0 <- function(pr, ps, bamIP, bamInput=NULL,removedup=TRUE, fragment = 300,
                          scaleControl=TRUE, libsize=c(1e6,1e6), removeBackground=TRUE,
                          paired.end= c("ignore", "filter")){

  if(removedup){#remove duplicate reads
    Flag=1024
  }else{
    Flag=-1
  }
  if(scaleControl){
    if(length(libsize)==1){
      scaleFactor=1
    }else{
      scaleFactor = libsize[1]/libsize[2]
    }
  }else{
    scaleFactor = 1
  }
  ### store the middle point of the peak and the end of the peak
  pr <- GenomeInfoDb::sortSeqlevels(pr)
  pr <- BiocGenerics::sort(pr)

  ps <- GenomeInfoDb::sortSeqlevels(ps)
  ps <- BiocGenerics::sort(ps)

  ### raw peak counts
  cIPRaw <- bamsignals::bamCount(bampath = bamIP, gr = ps, ss=FALSE, filteredFlag=Flag,
                                 shift= ceiling(fragment/2) )
  gc()
  ### common peak counts
  cIP <- bamsignals::bamCount(bampath = bamIP, gr = pr, ss=FALSE, filteredFlag=Flag,
                              shift= ceiling(fragment/2) )
  gc()

  if(!is.null(bamInput)){
    if(removeBackground){
      ### raw peak counts
      cCtrRaw = bamsignals::bamCount(bampath = bamInput, gr =ps, ss=FALSE, filteredFlag=Flag,
                                     shift= ceiling(fragment/2) )
      gc()
      ### common peak counts
      cCtr = bamsignals::bamCount(bampath = bamInput, gr =pr, ss=FALSE, filteredFlag=Flag,
                                  shift= ceiling(fragment/2) )
      gc()
      #ratio = (libsize[1]- sum(cIPRaw) )/(libsize[2] - sum(cCtrRaw)) # scale factor
      ### Here we used NCIS for calculate the scaling factor
      res <- NCIS(chip.data = bamIP, input.data = bamInput, data.type="BAM", frag.len=fragment)
      ratio = res$est
      count = ceiling(cIP - cCtr*ratio)

    }else{
      count = cIP
      ratio = rep(0, length(cIP))
      cInput = rep(0, length(cIP))
    }
    return(list(count=count, peakPos = pr, cIP=cIP, cInput=cCtr, ratio = ratio))
  }else{
    if(removeBackground){
      bamFile <- Rsamtools::BamFile(bamIP)
      chrInfo = GenomeInfoDb::seqinfo(bamFile)
      lens = GenomeInfoDb::seqlengths(bamFile)
      lenA = sum(lens)
      peakLen = pr@ranges@width # common peak length
      ratio = peakLen / (lenA - peakLen)#ratio
      backCount = (libsize[1] - sum(cIPRaw)) ## remove reads in raw peaks
      peakBack =  ratio * backCount ## background noise
      count = ceiling(cIP - peakBack)
      #hist(count, breaks = 100)
    }else{
      count = cIP
    }
    return(list(count=count, peakPos = pr, cIP=cIP, cInput=NULL, ratio = ratio))
  }
  # summary(count)
  # hist(count,breaks=512)
}

#' Get read counts from IP and Input bam file with given peak regions. To remove the backgroup noise,
#' compute the scale factor between the IP and input read counts from the non-peak regions
#' by by (All IP reads - reads of IP peaks regions ) / (All Input reads - reads of Input peaks regions)
#' The bam file must be indexed.
#' The final read count of peak region are
#' reads of IP peak regions - scale factor * reads of input peak region
#' User can also ignore input samples in the configuration files.
#' @param inputF file path, the file format is same to DiffBind input file to show the sample name,
#'        TF, Control, Treatment, bam file path, and peak files path.
#'        #contrast string vector, c(treatmentName, controlName)
#' @param overlap logic default false, use union of the peaks, TRUE, use overlapped peaks
#' @param comPeakFile string default NULL, the path of specified peak bed files. If it is not null,
#'                    @overlap parameters will no be used.
#' @param fragmentLen integer default is NULL, the fragment length will be determined by SGSeq::getBamInfo().
#'                 User can set it as 0, to avoid shift the reads.
#' @param removeBackground logical default is TRUE, remove background read count from each sample;
#'        For IP and Input, the IP will
#' @param removedup logic default true, remove duplicates
#'
#' @return list
#'         - countAll integer matrix, Read count of each sample after removing the background noise
#'         - peakPos GRanges objects, common peak coordinates between two conditions
#'         - fd data frame sample information
#'         - rawcount integer matrix  read counts of input and IP sample
#' @examples
#' readL = getReadCounts(inputF, contrast = c(treatmentName, controlName), peakLen=2,fragment=300)
#'
getReadCount <- function(inputF, overlap=FALSE, comPeakFile=NULL, fragmentLen=0,
                         removeBackground=TRUE, removedup=TRUE, 
                         removeUN_Random_Chr = FALSE,
                         nCore = 1){
  fd = read.table(inputF, sep=",", header = TRUE, stringsAsFactors = FALSE)
  #fd$Condition
  #fd$Peaks
  message(Sys.time(), ": 1. Get overlapped and merged peaks")
  pinfo = getPeakRegion(fd, removeUN_Random_Chr)# get overlapped and merged peaks
  
  if(overlap){
    po = pinfo$po
  }else{
    po = pinfo$pu
  }
  #0. Get specified peaks
  if(!is.null(comPeakFile)){# user specified peak
    df1 <- readr::read_delim(comPeakFile, delim="\t", col_names=FALSE)
    po <- GRanges(df1$X1, IRanges(df1$X2, df1$X3))
    if(removeUN_Random_Chr){
      po <- GenomeInfoDb::keepStandardChromosomes(po, pruning.mode="coarse")
    }
    po <- GenomeInfoDb::sortSeqlevels(po)
    po <- BiocGenerics::sort(po)
  }

  message(Sys.time(), ": 2. Get library size and fragment size")
  ipInfo = getLibsize(fd$SampleID, fd$bamReads, nCore=nCore)
  lsIP = ipInfo$lib_size
  frLenIP = ipInfo$frag_length
  read_lengthIP = ipInfo$read_length
  paired_endIP = ipInfo$paired_end
  if(!is.na(fd$bamControl[1])){# Get library size for Input
     inputInfo = getLibsize(fd$SampleID, fd$bamControl, nCore=nCore)
     lsCtr = inputInfo$lib_size
     frLenCtr = inputInfo$frag_length
     read_lengthInput = inputInfo$read_length
     paired_endInput = inputInfo$paired_end
  }else{## set 0 for Input if not existing
     lsCtr = rep(0, length(lsIP))
     frLenCtr = rep(0, length(lsIP))
     read_lengthInput = rep(0, length(lsIP))
     paired_endInput = rep(NA, length(lsIP))
  }

  fd$lsIP = lsIP
  fd$lsCtr = lsCtr
  fd$frLenIP = frLenIP
  fd$frLenCtr = frLenCtr
  fd$read_lengthIP = read_lengthIP
  fd$read_lengthInput = read_lengthInput
  fd$paired_endIP = paired_endIP
  fd$paired_endInput = paired_endInput

  if(is.null(fragmentLen)){## 1. Use the fragment detected by
    if(is.na(fd$frLenIP[1])){
      fragment =  0
    }else{
      fragment =  min(as.integer(fd$frLenIP))
    }
  }else if(length(fragmentLen)==1){# 2. Use a common fragment length for all samples
    #fragment = fragment #rep(fragment, length(lsIP))
    fragment = fragmentLen
  }else if(length(fragmentLen) == length(lsIP)){# 3. Same fragment number
    fragment = min(as.integer(fragmentLen))
  }else{# 4. Error
    stop("Fragment number is not equal to IP samples", call. = TRUE, domain = NULL)
    return()
  }

  print(length(po))

  message(Sys.time(), ": 3. Get read counts from IP samples")

  countAll = NULL
  countIP = NULL
  countCtr = matrix(0, nrow = length(po), ncol = length(lsCtr))
  countRatio = NULL # IP/Input
  sid= NULL
  for(i in 1:nrow(fd)){
    if(paired_endIP[i] == TRUE){
      pstr = "filter"
    }else{
      pstr = "ignore"
    }
    if(fd$lsCtr[i]==0){
      v1 = getPeakCounts(pr=po, 
                         bamIP = fd$bamReads[i],  bamInput = NA,
                         removedup = removedup, scaleControl = FALSE,
                         fragment=fragment, libsize=c(fd$lsIP[i], 0),
                         removeBackground=removeBackground, 
                         paired.end = pstr, removeUN_Random_Chr = removeUN_Random_Chr)
    }else{
      v1 = getPeakCounts(pr=po,  
                         bamIP = fd$bamReads[i],  bamInput = fd$bamControl[i],
                         removedup = removedup, scaleControl = TRUE,
                         fragment=fragment, libsize=c(fd$lsIP[i], fd$lsCtr[i]),
                         removeBackground=removeBackground, paired.end = pstr,
                         removeUN_Random_Chr = removeUN_Random_Chr
                         )
    }
    countAll = cbind(countAll, v1$count)
    countIP = cbind(countIP, v1$cIP)
    if(is.null(v1$cInput)){

    }else{
      countCtr[,i] =  v1$cInput
    }
    countRatio = cbind(countRatio, v1$ratio)
    sid = c(sid, paste0(fd$SampleID[i], "_", fd$Condition[i]) )
  }
  
  colnames(countIP) = paste0(sid, "_", "IP")
  colnames(countCtr) = paste0(sid, "_", "Input")
  colnames(countAll) = sid
  colnames(countRatio) = paste0(sid, "_", "Ratio")

  rn = paste0(as.character(po@seqnames), "_" , po@ranges@start, 
              "_" , (po@ranges@start + po@ranges@width-1) )
  
  rownames(countIP) = rn
  rownames(countCtr) = rn
  rownames(countAll) = rn
  
  countAll = as.matrix(countAll)
  countIP = as.matrix(countIP)
  countCtr = as.matrix(countCtr)
  peakPos = po
  peakUnion = pinfo$pu
  peakOverlap =  pinfo$po

  message(Sys.time(), ": 4. Finished")
  list(countAll=countAll, fd=fd, peakPos=peakPos, peakOverlap =peakOverlap,
       peakUnion=peakUnion,countIP=countIP, countCtr = countCtr , countRatio = countRatio)
}
 

#' Use ridge regression to shrink the variance and coefficients
#' @param cpmD matrix, normalized read count, row-peaks, col-samples
#' @param meanC vector from lmFit
#' @param coefAll matrix from lmFit
#' @param sigmaAll vector from lmFit
#' @param design0 integer The length of peak threshold (>2)
#' @param ps integer Select differential lambada value for ridge regression
#'
#' @return list
#'         - newCoef matrix, coefficient of ridge regression
#'         - newSigma2 vector, variance of coefficient
#'         - newResidualVar vector, residual variance fitting
#'         - newMSE vector, mean squared error
#' @examples
#' bv = getVB(cpmD, meanC,coefAll, sigmaAll, design0, ps=1)
#'
getVB <- function(cpmD, meanC,coefAll, sigmaAll, design0, ps=1){
  ##1. Given the design matrix Xd, the raw value oy, and the lambda
  ##   The fitted coefficients, sigma
  ##2. Calculate the new coefficients and sigma
  ##   new coefficients =  (1 + lambda)^-1 * coefficients
  ##   new lambda =  sigma^2 * (X'X + lambda)^-1 * X'*X* (X'X + lambda)^-1
  ##1.
  Xd = as.matrix(Xd)
  if(ps==1){
    lambdaAll = (1 - abs(meanC)/(max(meanC)))^2/10
  }else{
    logFC =coefAll[,2]
    svar = sigmaAll/max(sigmaAll)
    smean = abs(meanC)/max(abs(meanC))
    lambdaAll = (smean*svar)/(1+abs(logFC))
  }
  
  # plot(meanC, lambdaAll)
  nv = ncol(Xd)
  X2 = t(Xd) %*% Xd
  Xr = solve(X2)
  
  XL = sapply(1:nrow(cpmD), function(ix) {
    lx=lambdaAll[ix];
    v1=solve(X2 + lx*diag(nv));
    W_lambda = solve(X2 + lambdaAll[ix]*diag(nv)) %*% X2;
    #beta = v1 %*% X2 %*% t(v1);
    beta = v1 %*% t(Xd) %*% cpmD[ix,];
    #residual variance
    residualVar = sum((cpmD[ix,] - Xd %*% beta )^2)/(nrow(Xd) - ncol(Xd))
    ##variance of beta
    v2 = (W_lambda %*% Xr %*% t(W_lambda))
    varB = sigmaAll[ix]^2 * v2
    ##MSE
    m1 = sigmaAll[ix]^2 * sum(diag(v2))
    w1 = W_lambda - diag(ncol(Xd))
    mse0 = m1 + t(beta) %*% t(w1) %*% w1 %*% beta
    
    c(beta, residualVar, diag(varB), mse0[1,1])
  })
  
  newCoef = matrix(0, nrow = nrow(cpmD), ncol = ncol(Xd))
  newSigma2 = matrix(0, nrow = nrow(cpmD), ncol = ncol(Xd))
  newResidualVar = rep(0, nrow(cpmD))
  newMSE = rep(0, nrow(cpmD))
  rownames(newSigma2) = rownames(cpmD)
  rownames(newCoef) = rownames(cpmD)
  names(newResidualVar) = rownames(cpmD)
  names(newMSE) = rownames(cpmD)
  colnames(newCoef) = colnames(coefAll)
  colnames(newSigma2) = colnames(coefAll)
  
  newCoef[,1] = XL[1,]
  newCoef[,2] = XL[2,]
  newResidualVar = XL[3,]
  newSigma2[,1] = XL[4,]
  newSigma2[,2] = XL[5,]
  newMSE = XL[6,]
  return(list(newCoef=newCoef, newSigma2=newSigma2,
              newResidualVar=newResidualVar, newMSE=newMSE))
  
}


#' Use Cui's shrinkage estimator based James-Stein estimator to shrink the variance
#' Cui, X., Hwang, J.G., Qiu, J., Blades, N.J. and Churchill, G.A., 2005.
#' Improved statistical tests for differential gene expression by shrinking variance components estimates. Biostatistics, 6(1), pp.59-75.
#' @param fitq object from lmFit
#'
#' @return sigma2
#' @examples
#' bv = steinShrinkSigma(fitq)
#'
steinShrinkSigma <- function(fitq) {
  vv=c(1:20,25,30,40,50)
  Bv=c(3.53, 1.77, 1.44, 1.31, 1.24, 1.19, 1.16, 1.14, 1.12, 1.11,
       1.10, 1.09, 1.08, 1.08, 1.07, 1.07, 1.06, 1.06, 1.06, 1.05,
       1.04, 1.04, 1.03, 1.02)#B
  V_2_v = c(2.45, 1.64, 1.39, 1.27, 1.22, 1.18, 1.15, 1.13, 1.12, 1.11,
            1.10, 1.09, 1.08, 1.08, 1.07, 1.06, 1.06, 1.06, 1.05, 1.05,
            1.04, 1.03, 1.03, 1.02)# V/(2/v)

  sse = fitq$df.residual[1]*fitq$sigma^2 # Xg be the residual sum of squared errors (denoted by SSE)
  zid = which(sse > 1e-9) #
  zse = sse[zid]

  #2. Calculate (Xg/v)^1/G
  v0 = fitq$df.residual[1]
  id = which(vv==v0)
  gB = Bv[id]
  gV = 2/v0 * V_2_v[id]

  gN = length(zse)
  Xg_v = exp(1/gN * sum(log(zse/v0)))
  mean_ln_Xg = 1/gN * sum(log(zse))

  psum = 1 - ((gN - 3)*gV)/sum((log(zse) - mean_ln_Xg)^2)
  print(psum)
  if(psum < 0){
    psum = 1 # avoid overshrinkage
  }

  e2 =   exp(psum * (log(zse) - mean_ln_Xg) )
  summary(e2 + mean_ln_Xg)

  e3 = e2 + mean_ln_Xg
  summary(exp(e3))

  sigma0 = Xg_v * gB * e2
  sigma2 = fitq$sigma
  sigma2[zid] = sigma0

  plot(fitq$sigma, sqrt(sigma2),
       xlab="raw sigma", ylab="Shrunk sigma")
  sigma2
}


#' Use LOESS regression to correct the fold change
#' @param d0 matrix, normalized read count, row-peaks, col-samples
#' @param group vector the value of conditions
#' @param smean vector The mean value of each peak
#' @param sfold vector The logfold change
#' @param span numeric the parameter which controls the degree of smoothing.
#' @param plotT logical FALSE (default)
#' @param offSets logical TRUE (defulat). User can set it as false when the up and down number are similar.
#' @param sel integer Default is 2. User can set it as 1 because 1 is faster.
#' 
#' @return list
#'         - dnormV vector, shrunk logFold change
#'         - smoothV vector, fitted line
#'         - offSetValue numeric, the value to offset the logFold change
#' @examples
#' loessNormOffSet(d0, group, smean=NULL, sfold=NULL,span=0.6, offSets= TRUE, sel = 2, plotT = FALSE)
#'
loessNormOffSet <- function(d0, group, smean=NULL, sfold=NULL,
                            span1=0.6, offSets= TRUE,
                            plotT = FALSE, titleName=NULL){
  id1 = which(group==1)
  id2 = which(group==0)
  x1 = rowMeans(d0[,id1])
  x2 = rowMeans(d0[,id2])
  
  if(is.null(smean)){
    smean = (x1+x2)/2
  }
  if(is.null(sfold)){
    sfold = x1 - x2
  }
  ###1. correct fold change
  wt = smean/max(smean)
  nf = data.frame(x=smean)
  da = data.frame(x=smean, y=sfold)
  loessMod <- loess(y ~ x, data=da, span= span1) #
  smoothV <- predict(loessMod, newdata = nf) *wt
  summary(smoothV)
  
  if(plotT){
    df0 = data.frame(x = smean, sfold = sfold, smoothV = smoothV)
    p1 = ggplot(data=df0, aes(x=x, y=sfold)) +
      geom_point(size = 0.2, colour="#b3cccc") +
      geom_point(aes(x=x,y=smoothV), size = 0.1, colour="#ff6600") +
      geom_hline(yintercept=0, linetype="dashed") + ##, color = "red"
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="", x ="Mean (log2 CPM)", y = "Log2 Fold change")
    p1
  }
  
  ###3.Offset the fold change
  snormV = sfold - smoothV
 
  if(plotT){
    df0 = data.frame(x = smean, sfold = sfold, snormV = snormV)
    
    p1 = ggplot(data=df0, aes(x=x, y=sfold)) +
      geom_point(size = 0.2, colour="#b3cccc") +
      geom_point(aes(x=x,y=snormV), size = 0.1, colour="#ff6600") +
      geom_hline(yintercept=0, linetype="dashed") + ##, color = "red"
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="", x ="Mean (log2 CPM)", y = "Log2 Fold change")
    p1
  }
  
  if(max(smoothV) >= 0 & min(smoothV) <= 0){
    ths = c(seq(min(smoothV), 0, 0.002), seq(0, max(smoothV), 0.002))
  }else if( max(smoothV) <= 0 ){
    ths = seq(min(smoothV), 0, 0.002)
  }else if(min(smoothV) >= 0){
    ths = seq(0, max(smoothV), 0.002)
  }
  
  wt2 = 1 - abs(sfold)/max(abs(sfold))
  msv = ths
  csv = ths
  for(i in 1:length(ths)){
    f1 = snormV + ths[i]
    wt1 = 1 - abs(f1)/max(abs(f1))
    msv[i] = dist(rbind(f1*wt, sfold*wt))
    dm = data.frame(x = f1*wt, y=sfold*wt)
    ff = lm(y~x, data=dm )
    csv[i] = mean(ff$residuals^2)
  }
  
  if(1==0){
    plot(ths, csv, main="Number difference logFC and raw logFC")
    plot(ths, msv, main="Mean distance between \nfitted logFC and raw logFC")
  }
  
  if(offSets){
    dnormV = snormV +  ths[which.min(csv)]*wt
    
    if(plotT){
      df0 = data.frame(x = smean, sfold = sfold, smoothV = smoothV,
                       snormV = snormV, dnormV=dnormV)
      lineD = data.frame(x=ths, y = msv)
      
      p1 = ggplot(data=df0, aes(x=x, y=sfold)) +
        geom_point(size = 0.1, colour="#999999") + #999999 b3cccc
        geom_point(aes(x=x,y=smoothV), size = 0.1, colour="#ff6600") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
              axis.text = element_text(size = 15,face="bold"),
              axis.title=element_text(size=15,face="bold")) +
        labs(title="A", x ="Mean (log2 CPM)", y = "Log2 fold change")
      
      p2 = ggplot(data=df0, aes(x=x, y=sfold)) +
        geom_point(size = 0.1, colour="#999999") +
        geom_point(aes(x=x,y=snormV), size = 0.1, colour="#ff6600") + #0066ff
        theme_classic() +
        theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
              axis.text = element_text(size = 15,face="bold"),
              axis.title=element_text(size=15,face="bold")) +
        labs(title="B", x ="Mean (log2 CPM)", y = "Log2 fold change")
      
      p3 = ggplot(data=lineD, aes(x=x, y=y)) +
        geom_point(size = 1, colour="#999999") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
              axis.text = element_text(size = 15,face="bold"),
              axis.title=element_text(size=15,face="bold")) +
        labs(title="C", x ="Offset value", y = "Mean residual variance")
      
      p4 = ggplot(data=df0, aes(x=x, y=sfold)) +
        geom_point(size = 0.1, colour="#999999") +
        geom_point(aes(x=x,y=dnormV), size = 0.1, colour="#ff6600")+
        theme_classic() +
        theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
              axis.text = element_text(size = 15,face="bold"),
              axis.title=element_text(size=15,face="bold")) +
        labs(title="D", x ="Mean (log2 CPM)", y = "Log2 fold change")
      
      library(ggplot2)
      library(gridExtra)
      library(grid)
      g1 <- ggplotGrob(p1)
      g2 <- ggplotGrob(p2)
      g3 <- ggplotGrob(p3)
      g4 <- ggplotGrob(p4)
      gList1 = cbind(g1, g2)
      gList2 = cbind(g3, g4)
      gList = rbind(gList1, gList2)
      if(is.null(titleName)){
        titleName = as.character(proc.time()[1]*1000)
      }
      fn = paste0(titleName, "_LOESScorrect.png")
      #pdf(fn, width = 8, height = 6) # Open a new pdf file
      png(fn, width = 2000, height = 1600, res = 300)
      grid::grid.newpage()
      grid::grid.draw(gList)
      dev.off() # Close the file
    }
    
  }else{
    dnormV = snormV
    if(plotT){
      plot(smean, sfold,  main="Loess normalized without offset", cex=0.5,
           xlab="average", ylab="Diff")
      points(x=smean, dnormV, col="red")
      abline(h=0)
    }
  }
  list(dnormV=dnormV, snormV = snormV)
}



#' Fit the residual variance for low read counts
#' @param cpmD matrix, normalized read count, row-peaks, col-samples
#' @param design0 matrix the value of conditions
#' @param sx vector The mean value of each peak
#' @param sy vector The residual variance from lmFit()
#' @param winS numeric the parameter which controls window size, winS=0.02 by default.
#' @param minTh integer Default is 1. 
#'            If the change point is smaller than minTh, it will set change point as minTh
#' @param thN integer Default is 2. If minimum sx is larger than thN, it will return raw input.
#'            If the change point is larger than thN, it will set change point as thN.
#' @param plotT logical FALSE (default)
#' @param shrinkSigma vector The JS sigma for plot
#' 
#' @return list
#'         - sigma vector, fitted new residual variance 
#'         - fid vector, the id of the low count peaks
#' @examples
#' filtLowCounts(cpmD, design0, sx = fit1$Amean, sy = fit1$sigma,winS=0.02, thN=2, plotF=FALSE)
#'
filtLowCounts <- function(cpmDm, design0,
                          sx = fit1$Amean, sy = fit1$sigma,
                          minTh = 1, thN=2, winS=0.02, 
                          plotF=FALSE, titleName="", shrinkSigma) {
  cpmDm = as.matrix(cpmDm)
  # gid1 = which(design0[,2] == 0)
  # gid2 = which(design0[,2] == 1)
  
  if( min(sx) >= thN){
    return(L = list(sigma=sy, fid=NULL ) )
  }
  
  allS = seq(min(sx), max(sx), winS)
  if(allS[length(allS)] < max(sx) ){
    allS = c(allS, max(sx))
  }
  weightS = 0.8 + (1 - (allS - min(allS))/(max(allS) - min(allS)))*0.2
  winSigma = sy
  
  for(i in 2:length(allS)){
    id1 = which(sx >= allS[i-1] & sx < allS[i])
    if(length(id1) > 1){
      newWindowSigma = quantile(sy[id1], probs =weightS[i] )  #summary(sy[id1])
      winSigma[id1] = newWindowSigma 
    }
  }

  ##2.LOESS fit
  l0 <- lowess(sx,winSigma, f = 0.5)
  od = order(sx)
  summary(sx[od] - l0$x)
  smoothY = sx[od]
  smoothY[1:length(smoothY)] = l0$y
  smoothY = smoothY[names(sx)]
  
  ##2.Piecewise
  modelP <- piecewise.linear(sx,smoothY,middle = 1)
  id1 = which(modelP$x <= modelP$change.point)
  da1 = data.frame(x=modelP$x[id1], y=modelP$y[id1])
  da2 = data.frame(x=modelP$x[-id1], y=modelP$y[-id1])
  lf1 <- lm(y ~ x, data=da1)
  lf2 <- lm(y ~ x, data=da2)
  
  if(plotF){
    plot(modelP, cex = 0.1)
    abline(v= sx[which.max(smoothY)])
  }
  
  if(lf1$coefficients[2] < 0 & lf2$coefficients[2] < 0 ){ # both decrease, no change point
    if(min(sx) < minTh){
      midX = minTh
    }else{
      return(L = list(sigma=sy, fid=NULL ) )
    }
  }else if(lf1$coefficients[2] >= 0 & lf2$coefficients[2] <= 0){#increase then decrease
    midX = modelP$change.point
  }else if(lf1$coefficients[2] >= 0 & lf2$coefficients[2] >= 0){#increase only
    midX = modelP$change.point
  }else{
    midX = modelP$change.point
  }
    
  if(midX > thN){
    midX = thN
  }
  
  if(midX < minTh){
    midX = minTh
  }
  
  library(aomisc)
  ##4.Fit a line as variance of low expression data
  mid = which(sx >= midX)
  nd = data.frame(x = sx[mid], y = smoothY[mid])
  model <- drm(y ~ x, fct = DRC.expoDecay(),  data = nd) #
  # DRC.expoDecay EXD.2() DRC.powerCurve
  # Y = a * exp(k*X) NLS.expoDecay DRC.asymReg DRC.logCurve
  #model
  #model$predres
  if(1==0){
    plot(nd$x, nd$y, cex = 0.01 )
    points(nd$x, model$predres[,1], cex = 0.01, col="red")
  }
  rg = PR(model, sx)
  if(1==0){
    plot(sx,  y = sy, ylim=c(0,2), cex = 0.1)
    points(sx, rg, cex = 0.01, col="red")
  }
  fid = which(sx < midX )# & sy < rg
  sigma = sy
  sigma[fid] = rg[fid]
  
  if(plotF){ 
    if(length(fid) >0){
      shrinkSigma[fid] = rg[fid]
    }
    
    df0 = data.frame(x = sx, rawSigma = sy, winSigma = winSigma,
                     smoothSigma = smoothY, expY=rg, sigma=sigma, 
                     JSsigma=shrinkSigma)
    
    lineD = data.frame(x=modelP$x, y = modelP$y)
    p1 = ggplot(data=df0, aes(x=sx, y=rawSigma)) +
      geom_point(size = 0.1, colour="#999999") + #999999 b3cccc
      geom_point(aes(x=x,y=winSigma), size = 0.1, colour="#ff6600") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="A", x ="Mean (log2 CPM)", y = "Residual variance")
    
    p2 = ggplot(data=df0, aes(x=sx, y=rawSigma)) +
      geom_point(size = 0.1, colour="#999999") +
      geom_point(aes(x=x,y=smoothSigma), size = 0.1, colour="#0066ff") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="B", x ="Mean (log2 CPM)", y = "Residual variance")
    
    p3 = ggplot(data=df0, aes(x=sx, y=rawSigma)) +
      geom_point(size = 0.1, colour="#999999") +
      geom_point(aes(x=x,y=smoothSigma), size = 0.1, colour="#0066ff") +
      geom_line(data =lineD, aes(x=modelP$x,y=modelP$y), size = 0.5,colour="#ff0066") +
      geom_vline(xintercept = midX, linetype="dotted", color = "blue", size=0.4) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="C", x ="Mean (log2 CPM)", y = "Residual variance")
    
    p4 = ggplot(data=df0, aes(x=sx, y=rawSigma)) +
      geom_point(size = 0.1, colour="#999999") +
      geom_point(aes(x=x,y=expY), size = 0.1, colour="#ff6600")+
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="D", x ="Mean (log2 CPM)", y = "Residual variance")
    
    p5 = ggplot(data=df0, aes(x=sx, y=rawSigma)) +
      geom_point(size = 0.1, colour="#999999") +
      geom_point(aes(x=x,y=sigma), size = 0.1, colour="#ff6600")+
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="E", x ="Mean (log2 CPM)", y = "Residual variance")
    
    p6 = ggplot(data=df0, aes(x=sx, y=rawSigma)) +
      geom_point(size = 0.1, colour="#999999") +
      geom_point(aes(x=x,y=JSsigma), size = 0.1, colour="#ff6600")+
      theme_classic() +
      theme(plot.title = element_text(hjust = 0,size=20,face="bold"),
            axis.text = element_text(size = 15,face="bold"),
            axis.title=element_text(size=15,face="bold")) +
      labs(title="F", x ="Mean (log2 CPM)", y = "Residual variance")
    
    
    library(ggplot2)
    library(gridExtra)
    library(grid)
    g1 <- ggplotGrob(p1)
    g2 <- ggplotGrob(p2)
    g3 <- ggplotGrob(p3)
    g4 <- ggplotGrob(p4)
    g5 <- ggplotGrob(p5)
    g6 <- ggplotGrob(p6)
    gList1 = cbind(g1, g2, g3)
    gList2 = cbind(g4, g5, g6)
    gList = rbind(gList1, gList2)
    fn = paste0(titleName, "_fitLow.png")
    #pdf(fn, width = 8, height = 6) # Open a new pdf file
    png(fn, width = 2800, height = 1600, res = 300)
    grid::grid.newpage()
    grid::grid.draw(gList)
    dev.off() # Close the file
  }
  
  return(L = list(sigma=sigma, fid=fid, midX=midX) )
}


#' Do differential binding analysis by DiffChIPL
#' @param cpmD matrix, normalized CPM
#' @param design0 matrix The length of peak threshold (>2)
#' @param group0 vector the group information for two groups
#'
#' @return list,
#'         - fitlimma the fitting results of limma, refer the lmFit() in limma
#'         - fitDiffL the fitting results of DiffChIPL, refer the lmFit() in limma
#'         - resDE the DE results of DiffChIPL, refer the topTable() in limma
#'
#' @examples
#' resDE = DiffChIPL(cpmD, design0, group0 = group )
#'
DiffChIPL <- function(cpmD, design0, group0 = group, titleName="DiffChIPL"){
  #1.limma  linear regression
  fit1 <- lmFit(cpmD, design0, weights=NULL, method = "ls")
  rownames(fit1$coefficients) = rownames(cpmD)
  rownames(fit1$stdev.unscaled) = rownames(cpmD)
  names(fit1$Amean) = rownames(cpmD)
  names(fit1$sigma) = rownames(cpmD)

  #2.
  predY = t(as.matrix(design0) %*% t(fit1$coefficients))
  residY = (cpmD - predY)
  colnames(residY) = colnames(cpmD)
  library(corpcor)
  varD = var.shrink(t(residY))
  JS_Sigma = as.vector(sqrt(varD))
  
  #2.Fit the residual varaince for low read counts
  fitRlimmN = fit1
  filtL = filtLowCounts(cpmD, design0,
                        sx = fit1$Amean, sy = fit1$sigma, winS = 0.02,
                        plotF = FALSE, titleName = titleName, 
                        shrinkSigma =  JS_Sigma)
  ## update sigma
  fitRlimmN$sigma = JS_Sigma
  if(!is.null(filtL$fid)){
    fitRlimmN$sigma[filtL$fid] =  filtL$sigma[filtL$fid]
  }
   
  #3.Remove bias by LOESS regression
  resV = loessNormOffSet(d0 = cpmD, group= design0[,2], smean=fitRlimmN$Amean,
                         sfold= fitRlimmN$coefficients[,2], offSets = T)
  fitRlimmN$coefficients[,2] = resV$dnormV

  # 4. limma-trend
  fitRlimmR = fitRlimmN
  fitRlimmR <- eBayes(fitRlimmR, trend = TRUE, robust=TRUE)
  rtRlimmR = topTable(fitRlimmR, n = nrow(cpmD), coef=2)
  rtRlimmR = rtRlimmR[rownames(cpmD),]

  list(fitlimma = fit1, fitDiffL = fitRlimmR,  resDE = rtRlimmR)
}


#' Compare two vector distribution by histogram
#' @param v1 vector the value of v1
#' @param v2 vector the value of v2
#' @param name1 string name of the v1
#' @param name2 string name of the v2
#' @param xname string The axis label
#' @param binN integer 100 by default, bin number
#' @param tN numeric  0.6 by default, [0,1] for transparent
#'
#' @return ggplot2 object
#'
#' @examples
#' histOverlay(d0, group, xs=NULL, xlogfc=NULL,span=0.6, plotT = FALSE)
#'
histOverlay <- function(v1, v2, name1="v1", name2="v2",
                        xname="t-test value", binN = 100, tN = 0.6){
  library("ggplot2")
  c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
  c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
  da = data.frame(value = c(v1, v2),
                  group = c(rep(name1, length(v1)), rep(name2, length(v2))))

  ggplot(da, aes(x = value, fill = group)) +  # Draw overlaying histogram
    geom_histogram(position = "identity", alpha = tN, bins = binN) +
    #scale_fill_manual(values=c("#999999", "#E69F00")) +
    #scale_fill_grey() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15),
          axis.text.y = element_text( size=15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          strip.text = element_text(size = 15),
          plot.title = element_blank(),
          legend.text = element_text(size=18, face="bold"),
          legend.position="top") +
    labs(title="", x =xname, y = "Counts") +
    theme_classic()

}

#' Plot MA plot or volcano plot
#' @param mean vector
#' @param logfc vector
#' @param adj.P.Val vector
#' @param refDE vector known real DE annotation: raw, Up, and Down
#' @param FC numeric, draw logFC line with given FC value
#' @param padj numeric, filter the DE peaks with padj vlaue
#' @param psize numeric, adjust the point size
#' @param title string, the title of plot
#' @param MA logical, if is true, MA plot. else, volcano plot
#'
#' @return p2 ggplot
#'
#' @examples
#' bv = getVB(count, fit0, design)
#'
plotMAVoc <- function(mean, logfc, adj.P.Val, refDE = NULL, FC=0, padjT=0.05,psize=1, title="", MA=FALSE){
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

  res = data.frame(mean=mean,logFC=logfc,adj.P.Val=adj.P.Val)

  if(!is.null(refDE)){#
    shapeDE = rep("", length(refDE))
    shapeDE[which(refDE=="raw")] = paste0("No change:", length(which(refDE=="raw")))
    shapeDE[which(refDE=="Up")] = paste0("Real Up:", length(which(refDE=="Up")))
    shapeDE[which(refDE=="Down")] = paste0("Real Down:" , length(which(refDE=="Down")))
    res$shapeDE = shapeDE

    sizeDE = rep(2, length(refDE))
    sizeDE[which(refDE=="raw")] = 1
    res$sizeDE = sizeDE

    id1 = which(res$logFC > 0 & res$adj.P.Val < padjT)
    id2 = which(res$logFC < 0 & res$adj.P.Val < padjT)

    upLen = length(id1)
    dnLen = length(id2)
    NoLen = nrow(res) - upLen - dnLen

    realUp = length(id3)
    realDn = length(id4)
    realNo = nrow(res) - realUp - realDn

    # add a column of NAs
    res$diffexpressed <- paste0("NO:",NoLen)
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    res$diffexpressed[id1] <- paste0("UP:",upLen)
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    res$diffexpressed[id2] <- paste0("DOWN:",dnLen)

  }else{
    id1 = which(res$logFC > 0 & res$adj.P.Val < padjT)
    id2 = which(res$logFC < 0 & res$adj.P.Val < padjT)

    upLen = length(id1)
    dnLen = length(id2)
    NoLen = nrow(res) - upLen - dnLen

    # add a column of NAs
    res$diffexpressed <- paste0("NO:",NoLen)
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    res$diffexpressed[id1] <- paste0("UP:",upLen)
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    res$diffexpressed[id2] <- paste0("DOWN:",dnLen)

    res$shapeDE = as.factor(rep(1, nrow(res)))
    res$sizeDE = as.factor(rep(psize, nrow(res)))
  }

  if(MA){
    # Re-plot but this time color the points with "diffexpressed"
    p <- ggplot2::ggplot(data=res, aes(x=mean, y=logFC, col=diffexpressed, shape= shapeDE )) +
      geom_point(size = psize) + scale_shape_manual(values = c(1, 2, 6)) +
      theme_minimal()
    #geom_point( size = sizeDE, shape=shapeDE ) + theme_minimal()
    # Add lines as before...
    p2 <- p + geom_hline(yintercept=c(-FC, FC), col="red") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))

  }else{
    # Re-plot but this time color the points with "diffexpressed"
    p <- ggplot2::ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, shape= shapeDE)) +
      geom_point(size = psize) + scale_shape_manual(values = c(1, 2, 6)) + theme_minimal()
    # Add lines as before...
    p2 <- p + geom_vline(xintercept=c(-FC, FC), col="red") +
      geom_hline(yintercept=-log10(padjT), col="red") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  p2
}


#' Plot MA plot or volcano plot
#' @param mean vector
#' @param logfc vector
#' @param adj.P.Val vector
#' @param refDE vector known real DE annotation: raw, Up, and Down
#' @param FC numeric, draw logFC line with given FC value
#' @param padj numeric, filter the DE peaks with padj vlaue
#' @param psize numeric, adjust the point size
#' @param title string, the title of plot
#' @param MA logical, if is true, MA plot. else, volcano plot
#'
#' @return p2 ggplot
#'
#' @examples
#' bv = getVB(count, fit0, design)
#'
plotMAVoc2 <- function(mean, logfc, adj.P.Val, refDE = NULL, FC=0, padjT=0.05,psize=1, title="", MA=FALSE){
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

  res = data.frame(mean=mean,logFC=logfc,adj.P.Val=adj.P.Val)

  if(!is.null(refDE)){#
    sizeDE = rep(2, length(refDE))
    sizeDE[which(refDE=="raw")] = 1
    res$sizeDE = sizeDE

    id1 = which(res$logFC > 0 & res$adj.P.Val <= padjT & refDE == "up") # padj < TH and real Up
    id2 = which(res$logFC < 0 & res$adj.P.Val <= padjT & refDE == "down") # padj < TH and real Down
    id3 = which(res$adj.P.Val <= padjT & refDE == "raw") # padj < TH and real raw
    id4 = which( (res$logFC <= 0 & res$adj.P.Val <= padjT & refDE == "up")  |
                   (res$logFC >= 0 & res$adj.P.Val <= padjT & refDE == "down")  ) # padj < TH but wrong UP/Down

    id5 = which( res$adj.P.Val > padjT & refDE == "raw") # padj < TH and real no change
    id6 = which( res$adj.P.Val > padjT & refDE == "up") # padj < TH and real up
    id7 = which( res$adj.P.Val > padjT & refDE == "down") # padj < TH and real down

    length(id1) + length(id2) + length(id3) + length(id4) + length(id5) + length(id6) + length(id7)

    # add a column of NAs
    res$diffexpressed <- paste0("", nrow(res))
    if(length(id1) >0){
      res$diffexpressed[id1] <- paste0("padj <= ", padjT," and real Up:", length(id1))
    }
    if(length(id2) >0){
      res$diffexpressed[id2] <- paste0("padj <= ", padjT," and real Down:", length(id2))
    }
    if(length(id3) >0){
      res$diffexpressed[id3] <- paste0("padj <= ", padjT," and real no change:", length(id3))
    }
    if(length(id4) >0){
      res$diffexpressed[id4] <- paste0("padj <= ", padjT," but reverse logFC to real Up/Down:", length(id4))
    }
    if(length(id5) >0){
      res$diffexpressed[id5] <- paste0("padj > ", padjT," and real no change:", length(id5))
    }
    if(length(id6) >0){
      res$diffexpressed[id6] <- paste0("padj > ", padjT," and real Up:", length(id6))
    }
    if(length(id7) >0){
      res$diffexpressed[id7] <- paste0("padj > ", padjT," and real Down:", length(id7))
    }
  }else{
    id1 = which(res$logFC > 0 & res$adj.P.Val < padjT)
    id2 = which(res$logFC < 0 & res$adj.P.Val < padjT)

    upLen = length(id1)
    dnLen = length(id2)
    NoLen = nrow(res) - upLen - dnLen

    # add a column of NAs
    res$diffexpressed <- paste0("NO:",NoLen)
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    res$diffexpressed[id1] <- paste0("UP:",upLen)
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    res$diffexpressed[id2] <- paste0("DOWN:",dnLen)

    res$shapeDE = as.factor(rep(1, nrow(res)))
    res$sizeDE = as.factor(rep(psize, nrow(res)))
  }


  if(MA){
    # Re-plot but this time color the points with "diffexpressed"
    p <- ggplot2::ggplot(data=res, aes(x=mean, y=logFC, col=diffexpressed )) +
      geom_point(size = psize)  +
      scale_color_brewer(palette="Paired")  +
      #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
      theme_minimal()
    #geom_point( size = sizeDE, shape=shapeDE ) + theme_minimal()
    # Add lines as before...
    p2 <- p + geom_hline(yintercept=c(-FC, FC), col="red") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=14))

  }else{
    # Re-plot but this time color the points with "diffexpressed"
    p <- ggplot2::ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) +
      geom_point(size = psize) +
      scale_color_brewer(palette="Paired")  +
      theme_minimal()
    # Add lines as before...
    p2 <- p + geom_vline(xintercept=c(-FC, FC), col="red") +
      geom_hline(yintercept=-log10(padjT), col="red") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=14))
  }
  p2
}

#' Plot Box plot of normalized samples
#' @param df input matrix with colnames as samples
#' @param title string
#' @param xn name of x-axis
#' @param yn name of y-axis
#' @return p2 ggplot
#'
#' @examples
#' plotCorBoxs(count)
#'
plotCorBoxs <- function(df, title="", xn="Samples", yn="log2(CPM+1)"){
  comCorData <- data.frame(cbind(1:nrow(df), df))
  colnames(comCorData) = c("Peaks", colnames(df))
  corMelt <- reshape2::melt(comCorData, id.var='Peaks')
  colnames(corMelt) <- c("Peaks", "Samples", "Normalization")
  # 基函数
  ggplot2::ggplot(corMelt, aes(x = factor(Samples), y = Normalization, fill = factor(Samples))) +
    # 箱线图函数
    geom_boxplot(notch = TRUE) +
    theme(legend.text=element_text(size=15),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) +
    ggtitle(title) + xlab(xn) +
    ylab(yn) +  theme(plot.title = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size=10, angle = 90, hjust = 1))

}
