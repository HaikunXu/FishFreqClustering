library(FishFreqClustering)
library(FishFreqTree)
library(tidyverse)

LF1 <- LF1[,c(1,2,3,4,6:20)]
#analysis
bins <- seq(30, 170, 10) # data length bins
nbins <- length(bins)
fcol = 5
lcol = 19
save_dir="D:/OneDrive - IATTC/Git/FishFreqClustering/manual/"

make.meanl.map(LF1,fcol,lcol,bins,save_dir,width=10,height=10)
make.lf.map(LF1,fcol,lcol,bins,save_dir)

# divide by the mean for the year-quarter
LF2 <- lf.demean(LF1, fcol, lcol, bins)

mmd <- LF2[,c(2,4:20)]
min_samplesize <- 5

# setting up input data frames for clustering algorithm
temp = packbylatlon(mmd, 5, 5, nbins)  # lat, lon  
packedmmd3 = temp$table1
packedpdf3 = topdf(packedmmd3, 4, 3 + nbins)
packedcdf3 = tocdf(packedpdf3, 4, 3 + nbins)
mmdt = packedmmd3[packedmmd3[,4 + nbins] >= min_samplesize, ] # select the grids with more than 9 samples
rrs = mmdt[, 4 + nbins] 
mmdtpdf = packedpdf3[packedmmd3[,4 + nbins] >= min_samplesize, ] # remove the grids with no LF data
mmdtpdf[,4 + nbins] = mmdt[,4 + nbins]
mmdtcdf = packedcdf3[packedmmd3[,4 + nbins] >= min_samplesize, ] # remove the grids with no LF data
mmdtcdf[,4 + nbins] = mmdt[,4 + nbins]

densmatx = matrix(0, nrow(mmdt), nbins)
densmaty = matrix(0, nrow(mmdt), nbins)
for(i in 1:nrow(mmdt)){
  weightvec = t(mmdt[i,4:(3 + nbins)])
  weightvec = weightvec/sum(weightvec)
  # tempmmd = density(seq(0.00,2.01, 0.01), weights =weightvec,bw=0.10)
  densmatx[i,] = bins
  densmaty[i,] = t(weightvec)
}

# run distributional clustering with adjacency criterion
adjmat <- adjinf(mmdtpdf[,2],mmdtpdf[,3])
alydens.spatial23 <- hclust.regionsmm(as.matrix(densmaty),adj=TRUE,adjmat=adjmat,rr=sign(rrs))

# cplotu(alydens.spatial23$merges,alydens.spatial23$distseq,hopt='dist')

# making maps of the clusters and corresponding L-F density curves
#     kk is the number of clusters to use
kk<-3
par(mfrow=c(1,2),mar=c(4,4,1,1))

# map of clusters
temp <- putcolor(alydens.spatial23$merges, kk)
drawcells2(mmdtpdf[,2],mmdtpdf[,3],colseq=temp)

# draw density curves by cluster
colcol = rep(c(2,3, 4, 5, 6,7,8),3)
atitle.cl<-"2000-2022"
clusthistd3(kk,colseq=temp,atitle.cl,colcol, ylims = c(0, 0.3))

# save clustering results
cluster <- cbind(mmdt, temp)
names(cluster) <- c("Year", "Lat", "Lon", seq(30,170,10), "Nsamps", "Cell")

write.csv(cluster, file = "cluster_YFT.csv", row.names = FALSE)
