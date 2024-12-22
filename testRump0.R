library(GenomicRanges)
########## rename the chromosome in bed file
f1 = "Pc_all_sort.bed"
f2 = "Pc_all_sortNochr.bed"

d1 = read.csv(f1, header = FALSE, sep="\t")

r1 <- GenomicRanges::GRanges(d1[,1], IRanges(d1[,2], d1[,3]) )
r1 <- GenomeInfoDb::keepStandardChromosomes(r1, pruning.mode="coarse")
r1 = GenomeInfoDb::sortSeqlevels(r1)
r1 = BiocGenerics::sort(r1, ignore.strand=TRUE)

chrN = as.character(r1@seqnames)
chrS = as.character(sapply(chrN, function(x){strsplit(x, split = "chr")[[1]][2] } )  )
newD = data.frame(chr=chrS, start= r1@ranges@start, 
                  end = r1@ranges@start + r1@ranges@width -1)
write.table(newD, file = f2, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
