<<<<<<< HEAD
setwd("/Users/cheny48/Documents/Projects/diffA/DiffChIPL/")
=======
setwd("/Users/cheny48/Documents/Projects/diffA/Rlimma0/")
>>>>>>> c212010d18abf5b14c81cfb5bd5e628efd14a151
pth1 = getwd()
setwd(paste0(pth1, "/example/simHist/"))
getwd()

library(DiffChIPL)
flib="sim_hist_1.csv"
countL = getReadCount(inputF=flib,overlap=FALSE, fragment=0,
                      removeBackground=TRUE)
<<<<<<< HEAD

save(a = countL,file = "simHistCount.Data")

=======
>>>>>>> c212010d18abf5b14c81cfb5bd5e628efd14a151
fd = countL$fd
libsize = fd$lsIP

peakPos = countL$peakPos
peakAll = cbind(as.character(peakPos@seqnames), peakPos@ranges@start,
                peakPos@ranges@start+peakPos@ranges@width-1)
rawid = paste0(peakAll[,1],"_" ,peakAll[,2])
countAll = countL$countAll
rownames(countAll) = rawid

str1 = "sim-hist1-sim1.vs.sim2-Ridge"
group= c(1,1,0,0)
ctrName = "sim2"
treatName = "sim1"
<<<<<<< HEAD
groupName = c(treatName, treatName,ctrName,ctrName) 
=======
groupName = c(treatName, treatName,ctrName,ctrName)

deseqName = factor(c( "sim1", "sim1","sim2","sim2"), levels=c("sim2", "sim1"))
coldata= data.frame(group=deseqName)
rownames(coldata) = colnames(countAll)
coldata

>>>>>>> c212010d18abf5b14c81cfb5bd5e628efd14a151
design0 <- cbind(rep(1, 4), c(1,1,0,0))
colnames(design0) <- c(ctrName, treatName)
design0

for(i in 1:ncol(countAll)){
  id = which(countAll[,i] < 1)
  countAll[id,i] = 0
}
cpmD = cpmNorm(countAll, libsize = fd$lsIP)

###limma
fit3 <- lmFit(cpmD, design0, weights=NULL, method = "ls")
head(coef(fit3))
fit3 <- eBayes(fit3,trend=TRUE, robust = TRUE)
rt3 = topTable(fit3, n = nrow(cpmD), coef=2)

###DiffChIPL
resA = DiffChIPL(cpmD, design0, group0 = group )
fitRlimm3 = resA$fitDiffL
rtRlimm3 = resA$resDE

histOverlay(v1=fit3$t[,2], v2=fitRlimm3$t[,2],
            name1="limma", name2="DiffBindL", xname="Moderated t-test value")

histOverlay(v1=rt3$P.Value, v2=rtRlimm3$P.Value,binN = 20, tN = 0.6,
            name1="limma", name2="DiffBindL", xname="P-value")


refF = "hist1_refDE.bed"
d1 = read.table(refF, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
refDE = d1$V4

id_Rlimma_CPM = rownames(rtRlimm3[which(rtRlimm3$adj.P.Val < 0.05),])
rtRlimm3 = rtRlimm3[rownames(cpmD),]
aveE = rtRlimm3$AveExpr
logFC = rtRlimm3$logFC
padj = rtRlimm3$adj.P.Val
p = plotMAVoc2(mean=aveE, logfc=logFC, adj.P.Val=padj, FC=1,
               refDE= refDE, padj=0.05, MA=TRUE,
               title=paste0("DiffChIPL \n", str1,"(padj<0.05)\n",
                            length(id_Rlimma_CPM), " of ", nrow(rtRlimm3) ))
p + theme_classic()

<<<<<<< HEAD

=======
save(a = refDE, b =fit3, c= rtRlimm3, d=fitRlimm3,
     e = countL,
     file = "simHistRes.Data")
>>>>>>> c212010d18abf5b14c81cfb5bd5e628efd14a151


