# January 27 2024

# load L-F data (file below is an example; input must match this format)
load("D:/OneDrive - IATTC/IATTC/2024/Irregular clustering/From Cleridy/ps_yft_dph sets_2000-2017_for Mihoko.RData")

# attach function workspace from Mihoko
attach("D:/OneDrive - IATTC/IATTC/2024/Irregular clustering/From Cleridy/20230817_CL.Rdata")


# trim the data in space to match what Mihoko used originally (a limitation of this beta version)
#        5 degree squares, lower right corner: 15S to 25N and 135W to 75W
yftlengths.delsets.20002017.trmed<-yftlengths.delsets.20002017[yftlengths.delsets.20002017$lon.5deg>(-140) & yftlengths.delsets.20002017$lat.5deg>(-20),]


# trim the data to a block of years (note: in this beta version, some blocks of years may not run)
mmd<-yftlengths.delsets.20002017.trmed[yftlengths.delsets.20002017.trmed$year.firstset>=2013 & yftlengths.delsets.20002017.trmed$year.firstset<=2017,]


# setting up input data frames for clustering algorithm
temp = topmf(mmd, 5, 206)  
temp2 = tocdf(temp, 5, 206)
temp = packbylatlonnew(mmd)  # lat, lon  
packedmmd3 = temp$table1
packedpmf3 = topmf(packedmmd3, 5, 206)
packedcdf3 = tocdf(packedpmf3, 5, 206)
tempmm = apply(packedmmd3[,5:206],1,sum)
mmdt = packedmmd3[tempmm > 0.1, ]
rrs = mmdt[, 207] 
mmdtpmf = packedpmf3[tempmm > 0.1, ]
dimnames(mmdtpmf)[[2]] = dimnames(mmd)[[2]]
dimnames(mmdtpmf)[[2]][207] ="samplesize" 
mmdtpmf[,207] = mmdt[,207]
mmdtcdf = packedcdf3[tempmm > 0.1, ]
dimnames(mmdtcdf)[[2]] = dimnames(mmdtpmf)[[2]]
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
adjmat <- adjinf(mmdtpmf[,3],mmdtpmf[,4])
alydens.spatial23 <- hclust.regionsmm(as.matrix(densmaty),adj=TRUE,adjmat=adjmat,rr=rrs)


# making maps of the clusters and corresponding L-F density curves
#     kk is the number of clusters to use
kk<-4
par(mfrow=c(2,2),mar=c(4,4,1,1))

# map of clusters
temp <- putcolor(alydens.spatial23$merges, kk)
drawcells2(mmdtpmf[,3],mmdtpmf[,4],colseq=temp)

# draw density curves by cluster
colcol = rep(c(2,3, 4, 5, 6,7,8),3)
atitle.cl<-"2013-2017"
clusthistd3(4,colseq=temp,atitle.cl,colcol)


teststat <- heterodist(alydens.spatial23$merges,alydens.spatial23$distseq,tempmmd$x,densmaty,
                       sign(rrs),doko=c(1,5),BB=50)

teststats <- heterodists(alydens.spatial23$merges,alydens.spatial23$distseq,tempmmd$x,densmaty,
                         rrs,BB=c(100, 1000), bounds = c(0.005, 0.1))



