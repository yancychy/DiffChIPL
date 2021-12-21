library(DESeq2)
library(edgeR)
library(GenomicRanges)
library(Rsamtools)
library(bamsignals)
library(readr)
library(SGSeq)
#BiocManager::install("SGSeq")

getCP190Suhw <- function(){
  # CP190
  library(DiffBind)
  load("~/Documents/Projects/ChiPseqdata/comb/cp190_suhw_dba_count.RData") #
  #
  #dba.plotHeatmap(cp190_suhw)
  #r1 = dba.overlap(cp190_suhw,mode=DBA_OLAP_ALL)
  #md = dba.overlap(cp190_suhw, mask=cp190_suhw$masks$Consensus, mode=DBA_OLAP_PEAKS)
  # dba
  # dba.count
  rawcount = matrix(0, nrow = nrow(cp190_suhw$peaks[[1]]), ncol = length(cp190_suhw$peaks))
  dim(rawcount)
  for(i in 1:length(cp190_suhw$peaks)){
    ca = cp190_suhw$peaks[[i]]
    rawcount[,i] = as.vector(ca$cReads)
  }
  rawcount[1,]
  cp190_suhw$binding[1:3,]
  colnames(rawcount) = cp190_suhw$samples$SampleID
  rownames(rawcount) = paste(cp190_suhw$binding[,1],cp190_suhw$binding[,2], cp190_suhw$binding[,3],sep="_")

  lbsize = as.integer(cp190_suhw$class[8,])
  ### Based the edgeR PCA, only a small samples are used for differential anlaysis
  return(list(rawcount=rawcount[,c(5:7, 10:12)],libsize = lbsize[c(5:7, 10:12)]))
}


getCP190Mod <- function(){
  # CP190
  load("~/Documents/Projects/ChiPseqdata/comb/cp190_mod_dba_count.RData") #
  # dba
  # dba.count
  rawcount = matrix(0, nrow = nrow(cp190_mod$peaks[[1]]), ncol = length(cp190_mod$peaks))
  dim(rawcount)
  for(i in 1:length(cp190_mod$peaks)){
    ca = cp190_mod$peaks[[i]]
    rawcount[,i] = as.vector(ca$cReads)
  }
  colnames(rawcount) = cp190_mod$samples$SampleID
  rownames(rawcount) = paste(cp190_mod$binding[,1],cp190_mod$binding[,2], cp190_mod$binding[,3],sep="_")

  lbsize = as.integer(cp190_mod$class[8,])

  return(list(rawcount=rawcount[,c(5:7, 10:13)],libsize = lbsize[c(5:7, 10:13)]))
}


### Get read counts from bam file with given peak regions
## peakfile: The peak regions: chrname, start, end, in 0-base,For example,
## the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
## bamIP: IP bam files
## bamControl: Control bam files, if it is null, only counts the reads in IP, otherwise, reads = IP - control
## window=500, will count the reads in region [peakcenter - window , peakcenter + window], length is 2*window
## strand default is FALSE, it will count the reads whose 5' end map in it.
## If it is TRUE, count the positive strand reads in  [peakcenter - window, peakcenter-1],
## count the negative strand reads in [peakcenter, peakcenter + window]
getPeakCounts <- function(peakfile, bamIP, bamContorl, window=500, strand=FALSE){
  df1 <- read_delim(peakfile, delim="\t", col_names=FALSE)
  dim(df1)

  ### store the middle point of the peak and the end of the peak
  pr <- GRanges(df1$X1, IRanges(round((df1$X3+df1$X2)/2), df1$X3))
  pr <- sortSeqlevels(pr)
  pr <- sort(pr)

  if(!strand){ # counts all reads no matter the strand
    spr = promoters(pr, upstream=window, downstream=window, use.names=TRUE)
    cIP <- bamCount(bampath = bamIP, gr =spr, ss=FALSE)
    if(!is.null(bamContorl)){
      cCtr = bamCount(bampath = bamContorl, gr =spr, ss=FALSE)
      count = cIP - cCtr
    }else{
      count = cIP
    }
    # summary(count)
    # hist(count,breaks=512)
  }else{#
    sprUp = promoters(pr, upstream=window, downstream=0, use.names=TRUE)
    sprDn = promoters(pr, upstream=0, downstream=window, use.names=TRUE)
    cIP_UP <- bamCount(bampath = bamIP, gr =sprUp, ss=TRUE)
    cIP_DN <- bamCount(bampath = bamIP, gr =sprDn, ss=TRUE)

    if(!is.null(bamContorl)){
      cCtr_UP <- bamCount(bampath = bamContorl, gr =sprUp, ss=TRUE)
      cCtr_DN <- bamCount(bampath = bamContorl, gr =sprDn, ss=TRUE)
      count = cIP_UP[1,] + cIP_DN[2,] - (cCtr_UP[1,] + cCtr_DN[2,])
    }else{
      count = cIP_UP[1,] + cIP_DN[2,]
    }
  }
  count
}


getPeakCounts2 <- function(pr, bamIP, bamContorl, window=500, strand=FALSE){

  ### store the middle point of the peak and the end of the peak
  pr <- sortSeqlevels(pr)
  pr <- sort(pr)

  if(!strand){ # counts all reads no matter the strand
    spr = promoters(pr, upstream=window, downstream=window, use.names=TRUE)
    cIP <- bamCount(bampath = bamIP, gr =spr, ss=FALSE)
    if(!is.null(bamContorl)){
      cCtr = bamCount(bampath = bamContorl, gr =spr, ss=FALSE)
      count = cIP - cCtr
    }else{
      count = cIP
    }
    # summary(count)
    # hist(count,breaks=512)
    chr = as.vector( spr@seqnames)
    start = spr@ranges@start
    end = spr@ranges@width +  spr@ranges@start
    peakPos  =  cbind(chr, start, end)

  }else{#
    sprUp = promoters(pr, upstream=window, downstream=0, use.names=TRUE)
    sprDn = promoters(pr, upstream=0, downstream=window, use.names=TRUE)
    cIP_UP <- bamCount(bampath = bamIP, gr =sprUp, ss=TRUE)
    cIP_DN <- bamCount(bampath = bamIP, gr =sprDn, ss=TRUE)

    if(!is.null(bamContorl)){
      cCtr_UP <- bamCount(bampath = bamContorl, gr =sprUp, ss=TRUE)
      cCtr_DN <- bamCount(bampath = bamContorl, gr =sprDn, ss=TRUE)
      count = cIP_UP[1,] + cIP_DN[2,] - (cCtr_UP[1,] + cCtr_DN[2,])
    }else{
      count = cIP_UP[1,] + cIP_DN[2,]
    }

    chr = as.vector( sprUp@seqnames)
    start = sprUp@ranges@start
    end = 2*sprUp@ranges@width  +  sprUp@ranges@start
    peakPos  =  cbind(chr, start, end)
  }
  list(count=count, peakPos = peakPos)
}
##### Do overlap or merge the peaks
##### Do overlap or merge the peaks
##bedType - bed format, narrowPeak format,
##padj - 0.05, filter the peaks in narrowPeak files based the padjust value
getPeaks <- function(peakFiles, overlapOrMerge=TRUE, colNum=9, padj=0.05, bedType="bed"){

  daL = list()
  for(i in 1:length(peakFiles)){
    df1 <- read_delim(peakFiles[i], delim="\t", col_names=FALSE)
    if(bedType=="narrowPeak"){
      #9th: -log10qvalue at peak summit
      logTH = -log10(padj)
      df1 <- df1[which(df1[,9] > logTH),]
    }
    peak1 <- GRanges(df1$X1, IRanges(df1$X2, df1$X3))
    speak1 <- keepStandardChromosomes(peak1, pruning.mode="coarse")
    daL[[i]] = speak1
  }
  if(overlapOrMerge){## Compute the overlapped peaks among the replicates
    if(length(daL)==1){
      po = daL[[1]]
    }else{
      po = subsetByOverlaps(daL[[1]], daL[[2]])
      if(length(daL)>2){
        for(i in 3:length(peakFiles) ){
          po = subsetByOverlaps(po, daL[[i]])
        }
      }
    }
  }else{## Merge the peaks among different conditions
    if(length(daL)==1){
      po = daL[[1]]
    }else{
      po = union(daL[[1]], daL[[2]])
      if(length(daL)>2){
        for(i in 3:length(peakFiles) ){
          po = union(po, daL[[i]])
        }
      }
    }
  }
  po
}

getLibsize <- function(SampleID, bamFils, nCore=1){
  DF <- DataFrame( sample_name= SampleID,
                   file_bam = bamFils)
  #DF$file_bam <- BamFileList(DF$file_bam)
  DF_complete <- getBamInfo(DF, cores = nCore)
  DF_complete$lib_size
}


getPeakCounts3 <- function(pr, bamIP, bamContorl, window=500, strand=FALSE,
                           removedup=FALSE, scaleControl=TRUE, libsize=c(1e6,1e6)){
  if(removedup){
    Flag=1024
  }else{
    Flag=-1
  }
  if(scaleControl){
    scaleFactor = libsize[1]/libsize[2]
  }else{
    scaleFactor = 1
  }
  ### store the middle point of the peak and the end of the peak
  #pr <- GRanges(df1$X1, IRanges(round((df1$X3+df1$X2)/2), df1$X3))
  pr <- sortSeqlevels(pr)
  pr <- sort(pr)

  if(!strand){ # counts all reads no matter the strand
    spr = promoters(pr, upstream=window, downstream=window, use.names=TRUE)
    cIP <- bamCount(bampath = bamIP, gr =spr, ss=FALSE, filteredFlag=Flag)


    if(!is.null(bamContorl)){
      cCtr = bamCount(bampath = bamContorl, gr =spr, ss=FALSE)
      count = ceiling(cIP - cCtr*scaleFactor)
    }else{
      count = cIP
    }
    # summary(count)
    # hist(count,breaks=512)
    chr = as.vector( spr@seqnames)
    start = spr@ranges@start
    end = spr@ranges@width +  spr@ranges@start
    peakPos  =  cbind(chr, start, end)

  }else{#
    sprUp = promoters(pr, upstream=window, downstream=0, use.names=TRUE)
    sprDn = promoters(pr, upstream=0, downstream=window, use.names=TRUE)
    cIP_UP <- bamCount(bampath = bamIP, gr =sprUp, ss=TRUE)
    cIP_DN <- bamCount(bampath = bamIP, gr =sprDn, ss=TRUE)

    if(!is.null(bamContorl)){
      cCtr_UP <- bamCount(bampath = bamContorl, gr =sprUp, ss=TRUE)
      cCtr_DN <- bamCount(bampath = bamContorl, gr =sprDn, ss=TRUE)
      count = ceiling(cIP_UP[1,] + cIP_DN[2,] - (cCtr_UP[1,] + cCtr_DN[2,])*scaleFactor)
    }else{
      count = cIP_UP[1,] + cIP_DN[2,]
    }

    chr = as.vector( sprUp@seqnames)
    start = sprUp@ranges@start
    end = 2*sprUp@ranges@width + sprUp@ranges@start
    peakPos  =  cbind(chr, start, end)
  }
  list(count=count, peakPos = peakPos)
}

getPeakCounts4 <- function(pr, bamIP, bamContorl, window=500, strand=FALSE,
                           removedup=TRUE,fragment = 300, scaleControl=TRUE, libsize=c(1e6,1e6)){

  if(removedup){
    Flag=1024
  }else{
    Flag=-1
  }
  if(scaleControl){
    scaleFactor = libsize[1]/libsize[2]
  }else{
    scaleFactor = 1
  }
  ### store the middle point of the peak and the end of the peak
  pr <- sortSeqlevels(pr)
  pr <- sort(pr)


  cIP <- bamCount(bampath = bamIP, gr =pr, ss=FALSE, filteredFlag=Flag,
                    shift= ceiling(fragment/2) )

  if(!is.null(bamContorl)){
      cCtr = bamCount(bampath = bamContorl, gr =pr, ss=FALSE, filteredFlag=Flag,
                      shift= ceiling(fragment/2) )
      ratio = (libsize[1]- sum(cIP))/(libsize[2] - sum(cCtr)) # scale factor
      count = ceiling(cIP - cCtr*ratio)
      return(list(count=count, peakPos = pr, cIP=cIP, cCtr=cCtr))
  }else{
      ratio = (libsize[1]- sum(cIP))/libsize[1] # background noise ratio
      count = ceiling(cIP*(1-ratio))
      return(list(count=count, peakPos = pr, cIP=cIP, cCtr=NULL))
  }
  # summary(count)
  # hist(count,breaks=512)

}


