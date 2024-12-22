library(DiffChIPL)
flib="PC_rumpkd_noInput.csv"

countL = getReadCount(inputF=flib,overlap=FALSE, fragmentLen =250,
                      removeBackground=TRUE, removeUN_Random_Chr = TRUE)
countL$fd
countL$countAll[1:3,]
countL$countIP[1:3,]
countL$countCtr[1:3,]
countL$peakPos[1:3]

save(a = countL,file = "PC_rumpkd_noInput.Data")

load("PC_rumpkd_noInput.Data")

fd = countL$fd
libsize = fd$lsIP

peakPos = countL$peakPos
peakAll = cbind(as.character(peakPos@seqnames), peakPos@ranges@start,
                peakPos@ranges@start+peakPos@ranges@width-1)
rawid = paste0(peakAll[,1],"_" ,peakAll[,2])
countAll = countL$countAll
rownames(countAll) = rawid

str1 = "Pc_mock.vs.rumpkd"
group= c(0,0,0,1, 1,1)
ctrName = "mock"
treatName = "rumpkd"

design0 <- cbind(rep(1, 6), c(0,0,0,1, 1,1))
colnames(design0) <- c(ctrName, treatName)
design0

for(i in 1:ncol(countAll)){
  id = which(countAll[,i] < 1)
  countAll[id,i] = 0
}
cpmD = cpmNorm(countAll, libsize = fd$lsIP)

###DiffChIPL
resA = DiffChIPL(cpmD, design0, group0 = group )
fitRlimm3 = resA$fitDiffL
rtRlimm3 = resA$resDE

