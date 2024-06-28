library(FishFreqClustering)
library(FishFreqTree)
library(tidyverse)

Raw <- read.csv("D:/OneDrive - IATTC/Git/FishFreqClustering/manual/yft_lf_2000_to_2022.csv")
Raw$quarter = ceiling(Raw$month / 3)
Raw$lat = Raw$lat.5deg
Raw$lon = Raw$lon.5deg

LF.OBJ <- Raw %>% filter(class == 6, setype == 1, lat > -25) # 1=DEL; 4=NOA; 5=OBJ

LF <- LF.OBJ[, c("year", "quarter", "lat", "lon", paste0("X", 1:201))] %>%
  group_by(lat, lon) %>% mutate(N = length(unique(paste0(year,"-",quarter))))

plot(LF$lon, LF$lat)

#analysis
bins <- seq(1, 201, 1) # data length bins
nbins <- length(bins)

LF1 <-
  lf.aggregate(
    LF,
    fcol = 5,
    lcol = 205,
    bins,
    bins,
    LengthOnly = FALSE
  )

# divide by the mean for the year-quarter
LF2 <- lf.demean(LF1, 5, 205, bins = bins)

mmd <- LF2[,c(2,4:206)]
min_samplesize <- 50

# setting up input data frames for clustering algorithm
temp = packbylatlon(mmd, 5, 5, nbins)  # lat, lon  
packedmmd3 = temp$table1
packedpdf3 = topdf(packedmmd3, 3, 3 + nbins)
packedcdf3 = tocdf(packedpdf3, 3, 3 + nbins)
mmdt = packedmmd3[packedmmd3[,4 + nbins] >= min_samplesize, ] # select the grids with more than 9 samples
rrs = mmdt[, 4 + nbins] 
mmdtpdf = packedpdf3[packedmmd3[,4 + nbins] >= min_samplesize, ] # remove the grids with no LF data
mmdtpdf[,4 + nbins] = mmdt[,4 + nbins]
mmdtcdf = packedcdf3[packedmmd3[,4 + nbins] >= min_samplesize, ] # remove the grids with no LF data
mmdtcdf[,4 + nbins] = mmdt[,4 + nbins]

densmatx = matrix(0, nrow(mmdt), 512)
densmaty = matrix(0, nrow(mmdt), 512)
for(i in 1:nrow(mmdt)){
  weightvec = t(mmdt[i,3:(3 + nbins)])
  weightvec = weightvec/sum(weightvec)
  tempmmd = density(seq(0.00,2.01, 0.01), weights =weightvec,bw=0.10)
  densmatx[i,] = tempmmd$x
  densmaty[i,] = tempmmd$y
}

# run distributional clustering with adjacency criterion
adjmat <- adjinf(mmdtpdf[,1],mmdtpdf[,2])
alydens.spatial23 <- hclust.regionsmm(as.matrix(densmaty),adj=TRUE,adjmat=adjmat,rr=(rrs))

cplotu(alydens.spatial23$merges,alydens.spatial23$distseq,hopt='dist')

# making maps of the clusters and corresponding L-F density curves
#     kk is the number of clusters to use
kk<-3
par(mfrow=c(1,2),mar=c(4,4,1,1))

# map of clusters
temp <- putcolor(alydens.spatial23$merges, kk)
drawcells2(mmdtpdf[,1],mmdtpdf[,2],colseq=temp)

# draw density curves by cluster
colcol = rep(c(2,3, 4, 5, 6,7,8),3)
atitle.cl<-"2000-2022"
clusthistd3(kk,colseq=temp,atitle.cl,colcol)

