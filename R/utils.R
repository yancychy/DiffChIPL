##
getMapID <- function(rawid, cid){
  sapply(cid, function(x){ which(rawid==x)  })
}

#The Stirling Formula is used to compute log n!
stirling <- function(x){
  1/2*log(2*pi)+(1/2+x)*log(x)-x
}

outlier = function(x, s){
  # x = data; s = number of outliers
  n = length(x) - s
  sigma = sqrt(sum((x-mean(x))^2)/n)
  ###### see footnote 4
  ut = n*log(sigma)+(sqrt(2)*stirling(n))/n*s
  ut

  if(1==0){
    ss = seq(0,100,1)
    aicV = ss
    for(i in 1:length(ss)){
      aicV[i] = outlier(scale(x),ss[i])
    }
    plot(ss,aicV)
  }
}



plotG3D <- function(da, cs, gid){

  q2 <-  data.frame(x = da[,1], y = da[,2], z= da[,3])
  plot_ly(q2, x=~x, y=~y,  z=~z,  mode="markers", color = as.character(gid),
          colors=cs, marker = list(size = 6,opacity=1), scene=gene)

  library(RColorBrewer)
  library(colorRamps)
  library(plotly)
  #coid = which(rownames(brewer.pal.info)==colorname)
  #colors = colorRampPalette(brewer.pal(brewer.pal.info[coid,1], colorname))(length(unique(cluidA)))
  pointsize = 7
  fontsize = 10
  labelsize = 20
  da = data.frame(x = da[,1], y = da[,2], z= da[,3])

  p <- plot_ly() %>%
    add_trace(data = da, x= ~x, y= ~y, z= ~z,
              color =  gid,
              type="scatter3d", colors = cs,
              marker = list(opacity= 0.9, size=fontsize),showlegend = TRUE)

  p <- plotly_build(p)
  markersize <- pointsize
  markerlegendsize <- labelsize
  for(i in seq(1, length(sort(unique(gid) )))) {
    length.group <- nrow(da[which(gid == sort(unique(gid))[i]), ])
    p$x$data[[i]]$marker$size <- c(rep(markersize,length.group), rep(c(-markersize+2*markerlegendsize), length.group))
  }
  p %>%
    plotly::layout(title=list(font=list(family = "sans serif", size = 20,
                                        color = 'Black'), text=paste0("Gene: ", gene) ))
}


