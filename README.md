# DiffChIPL
DiffChIPL: A Differential binding analysis method based on limma for ChIP-seq data featuring biological replicates


## Install
1. User can install source package by Rstudio
2. User can also install the package by 

```R
library(devtools)
install_github("yancychy/DiffChIPL")
```
## Worfflow
![workflow](https://github.com/yancychy/DiffChIPL/blob/main/example/workflow1.jpg)

## Example

### 1. Make a configuration file
The example configuration files are shown in [example/simHist/sim_hist_1.csv](https://github.com/yancychy/DiffChIPL/blob/main/example/simHist/sim_hist_1.csv).

Commonly, the configuration file of DiffBind can be used directly.

### 2.Get read count from simulated histone datasets
We used the simulation code in [csaw](http://bioinf.wehi.edu.au/csaw/) to simulate the histone which is a mixture of complex peak structures.
We simulated 10000 peaks which have 1000 increased peaks  and 1000 decreased peaks.
The simulation histone has two groups. In each group, there has two replicates for each histone condition. 

```R
library(DiffChIPL)
flib="sim_hist_1.csv"
countL = getReadCount(inputF=flib,overlap=FALSE, fragment=0,
                      removeBackground=TRUE)
```

### 3. Make the design matrix and normalization

```R
str1 = "sim-hist1-sim1.vs.sim2-Ridge"
group= c(1,1,0,0)
ctrName = "sim2"
treatName = "sim1"
groupName = c(treatName, treatName,ctrName,ctrName) 
design0 <- cbind(rep(1, 4), c(1,1,0,0))
colnames(design0) <- c(ctrName, treatName)
design0
 
for(i in 1:ncol(countAll)){
  id = which(countAll[,i] < 1)
  countAll[id,i] = 0
}
cpmD = cpmNorm(countAll, libsize = fd$lsIP)

```

### 4. Do differential analysis with DiffChIPL

```R
resA = DiffChIPL(cpmD, design0, group0 = group )
fitRlimm3 = resA$fitDiffL
rtRlimm3 = resA$resDE
```


### 5. Check the differential results

```R
id_Rlimma_CPM = rownames(rtRlimm3[which(rtRlimm3$adj.P.Val < 0.05),])
 
rtRlimm3 = rtRlimm3[rownames(cpmD),]
aveE = rtRlimm3$AveExpr
logFC = rtRlimm3$logFC
padj = rtRlimm3$adj.P.Val
plotMAVoc2(mean=aveE, logfc=logFC, adj.P.Val=padj, FC=1, 
           refDE= refDE, padj=0.05, MA=TRUE,
          title=paste0("Rlimma-CPM \n", str1,"(padj<0.05)\n", 
                       length(id_Rlimma_CPM), " of ", nrow(rtRlimm3) ))

```  
![MA plot](https://github.com/yancychy/DiffChIPL/blob/main/example/simHist/MA_DiffChIPL.png)



