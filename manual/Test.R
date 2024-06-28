library(FishFreqClustering)
library(FishFreqTree)
library(tidyverse)

Raw <- read.csv("D:/OneDrive - IATTC/IATTC/2024/SAC15/Assessment/Data/PS tree/bet lf_2000 to 2022.csv")
Raw$quarter = ceiling(Raw$month / 3)
Raw$lat = Raw$lat.5deg
Raw$lon = Raw$lon.5deg

LF.OBJ <- Raw %>% filter(class == 6, setype == 5, lat > -25) # 1=DEL; 4=NOA; 5=OBJ

LF <- LF.OBJ[, c("year", "quarter", "lat", "lon", paste0("X", 1:201))] %>%
  group_by(lat, lon) %>% mutate(N = length(unique(paste0(year,"-",quarter)))) %>%
  filter(N > 9)

plot(LF$lon, LF$lat)

#analysis
bins <- seq(1, 201, 1) # data length bins

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

# setting up input data frames for clustering algorithm
temp = packbylatlon(mmd, 5, 5, 201)  # lat, lon  
packedmmd3 = temp$table1
packedpdf3 = topdf(packedmmd3, 5, 206)
packedcdf3 = tocdf(packedpdf3, 5, 206)
tempmm = apply(packedmmd3[,5:206],1,sum)
mmdt = packedmmd3[tempmm > 0.1, ]
rrs = mmdt[, 207] 
mmdtpdf = packedpdf3[tempmm > 0.1, ] # remove the grids with no LF data
dimnames(mmdtpdf)[[2]] = dimnames(mmd)[[2]]
dimnames(mmdtpdf)[[2]][207] ="samplesize" 
mmdtpdf[,207] = mmdt[,207]
mmdtcdf = packedcdf3[tempmm > 0.1, ] # remove the grids with no LF data
dimnames(mmdtcdf)[[2]] = dimnames(mmdtpdf)[[2]]
mmdtcdf[,207] = mmdt[,207]

densmatx = matrix(0, nrow(mmdt), 512)
densmaty = matrix(0, nrow(mmdt), 512)
for(i in 1:nrow(mmdt)){
  weightvec = t(mmdt[i,5:206])
  weightvec = weightvec/sum(weightvec)
  tempmmd = density(seq(0.00,2.01, 0.01), weights =weightvec,bw=0.10)
  densmatx[i,] = tempmmd$x
  densmaty[i,] = tempmmd$y
}

# run distributional clustering with adjacency criterion
adjmat <- adjinf(mmdtpdf[,3],mmdtpdf[,4])
alydens.spatial23 <- hclust.regionsmm(as.matrix(densmaty),adj=TRUE,adjmat=adjmat,rr=(rrs))

cplotu(alydens.spatial23$merges,alydens.spatial23$distseq,hopt='dist')

# making maps of the clusters and corresponding L-F density curves
#     kk is the number of clusters to use
kk<-5
par(mfrow=c(1,2),mar=c(4,4,1,1))

# map of clusters
temp <- putcolor(alydens.spatial23$merges, kk)
drawcells2(mmdtpdf[,3],mmdtpdf[,4],colseq=temp)

# draw density curves by cluster
colcol = rep(c(2,3, 4, 5, 6,7,8),3)
atitle.cl<-"2000-2022"
clusthistd3(kk,colseq=temp,atitle.cl,colcol)

