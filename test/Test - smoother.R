library(FishFreqClustering)
library(FishFreqTree)
library(tidyverse)

directory <- "D:/OneDrive - IATTC/IATTC/2024/Irregular clustering/YFT DEL/"
setwd(directory)

Raw <- read.csv("yft_lf_2000_to_2022.csv")
Raw$quarter = ceiling(Raw$month / 3)
Raw$lat = Raw$lat.5deg + 2.5
Raw$lon = Raw$lon.5deg - 2.5

LF.DEL <- Raw %>% filter(class == 6, setype == 1) # 1=DEL; 4=NOA; 5=OBJ
LF <- LF.DEL[, c("year", "quarter", "lat", "lon", paste0("X", 1:201))] %>%
  group_by(lat, lon) %>% mutate(N = length(unique(paste0(year, "-", quarter)))) %>%
  filter(N > 4)

bins <- seq(1, 201, 1) # data length bins
new_bins <- seq(1, 201, 1) # bins to  be used in the analysis

LF1 <-
  lf.aggregate(
    LF,
    fcol = 5,
    lcol = 205,
    bins,
    new_bins,
    LengthOnly = FALSE
  )

# LF1 <- LF1[,c(1,2,3,4,6:20)]
#analysis
nbins <- length(bins)
fcol = 5
lcol = 205
save_dir=directory

make.meanl.map(LF1,fcol,lcol,bins,save_dir,width=10,height=10)
# make.lf.map(LF1,fcol,lcol,bins,save_dir)

# divide by the mean for the year-quarter
LF2 <- lf.demean(LF1, fcol, lcol, bins)

mmd <- LF2[,c(2,4:206)]
# mmd <- LF1[,c(1,3:20)]
min_samplesize <- 1

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

densmatx = matrix(0, nrow(mmdt), 512)
densmaty = matrix(0, nrow(mmdt), 512)
for(i in 1:nrow(mmdt)){
  weightvec = t(mmdt[i,4:(3 + nbins)])
  weightvec = weightvec/sum(weightvec)
  tempmmd = density(seq(0.01,2.01, 0.01), weights =weightvec,bw=0.10)
  densmatx[i,] = tempmmd$x
  densmaty[i,] = tempmmd$y
}

# run distributional clustering with adjacency criterion
adjmat <- adjinf(mmdtpdf[,2],mmdtpdf[,3])
alydens.spatial23 <- hclust.regionsmm(as.matrix(densmaty),adj=TRUE,adjmat=adjmat,rr=(rrs))

teststat <- heterodist(alydens.spatial23$merges,alydens.spatial23$distseq,tempmmd$x,densmaty,
                       rrs,doko=c(1,5),BB=50)
write.csv(teststat$statmat, file = "MJS distance.csv", row.names = FALSE)

teststats <- heterodists(alydens.spatial23$merges,alydens.spatial23$distseq,tempmmd$x,densmaty,
                         rrs,BB=c(10, 100), bounds = c(0.005, 0.1))

# cplotu(alydens.spatial23$merges,alydens.spatial23$distseq,hopt='dist')

# making maps of the clusters and corresponding L-F density curves
#     kk is the number of clusters to use
kk<-5
# par(mfrow=c(1,2),mar=c(4,4,1,1))

# # map of clusters
temp <- putcolor(alydens.spatial23$merges, kk)
# drawcells2(mmdtpdf[,2],mmdtpdf[,3],5,5,colseq=temp)
# 
# # draw density curves by cluster
# colcol = rep(c(2,3, 4, 5, 6,7,8),3)
# atitle.cl<-"2000-2022"
# clusthistd3(kk,colseq=temp,atitle.cl,colcol, ylims = c(0, 0.3))

# save clustering results
cluster <- cbind(mmdt[,2:3], factor(temp-1))
names(cluster) <- c("lat", "lon", "cell")

write.csv(cluster, file = "cluster_YFT.csv", row.names = FALSE)

wmap <- map_data("world")
ggplot(data=cluster) +
  geom_tile(aes(x=lon, y=lat, fill=cell), color = "black") +
  geom_polygon(data=wmap,aes(long, lat, group = group),fill = "black",colour = "white",lwd=0.5) +
  coord_quickmap(ylim = c(min(cluster$lat),max(cluster$lat)),xlim = c(min(cluster$lon),max(cluster$lon))) +
  theme_bw()

LF1_cluster <- left_join(LF1, cluster) %>%
  gather(5:20, key = "Length", value = LF) %>%
  mutate(Length = as.numeric(Length)) %>%
  group_by(cell, Length) %>%
  summarise(Mean_LF = mean(LF)) %>%
  group_by(cell) %>%
  mutate(Mean_LF = Mean_LF / sum(Mean_LF))

ggplot(data = LF1_cluster) +
  geom_line(aes(x = Length, y = Mean_LF, color = cell)) +
  geom_point(aes(x = Length, y = Mean_LF, color = cell)) +
  theme_bw()

LF1_cluster <- left_join(LF1, cluster) %>%
  rename(Flag = cell)
make.lf.cell(LF1_cluster,fcol,lcol,bins,save_dir,plot_name = "NewLF")

# compare it with the regression tree result
save_dir <- directory

loop_dir <-paste0(save_dir,"loop/")
dir.create(loop_dir)

my_select_matrix <-
  data.matrix(expand.grid(
    split1 = 1:2,
    split2 = 1:2,
    split3 = 1:1,
    split4 = 1:1
  ))

LF_Tree_Loop <-
  loop_regression_tree(LF2,
                       fcol+1,
                       lcol+1,
                       bins,
                       Nsplit=4,
                       save_dir = loop_dir,
                       select_matrix = my_select_matrix,
                       # year = TRUE,
                       lat.min = 2, # a fishery area should span at least 10 degree
                       lon.min = 2) # a fishery area should span at least 10 degree

# best one
LF_Tree <-
  run_regression_tree(
    LF1,
    fcol,
    lcol,
    bins,
    Nsplit = 4,
    save_dir,
    lat.min = 2,
    lon.min = 2,
    manual = TRUE,
    select = c(1, 1, 1, 1)
  )
make.split.map(
  LF_Tree$LF,
  Nsplit = 4,
  save_dir,
  width = 10,
  height = 10,
  plot_name = "Old"
)

LF1$Flag <- LF_Tree$LF$Flag3
make.lf.cell(LF1,fcol,lcol,bins,save_dir,plot_name = "OldLF")

LF_test <- left_join(LF_Tree$LF, cluster) %>%
  mutate(Flag5 = as.numeric(cell))
make.split.map(
  LF_test,
  Nsplit = 5,
  save_dir,
  width = 10,
  height = 10,
  plot_name = "New"
)

evaluate_regression_tree(LF_test, 5, 20, "Flag4",
                         bins, directory, folder_name = "New")
