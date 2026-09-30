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
  group_by(lat, lon) %>%
  mutate(N = length(unique(paste0(year, "-", quarter)))) %>%
  filter(N > 3, lat > -10) # remove the cells with less than 4 quarters of data since 2000

bins <- seq(1, 201, 1) # data length bins
new_bins <- seq(61, 200, 1) # bins to be used in the clustering analysis

# first aggregate the raw LF to the new bins by quarter
LF1 <- lf.aggregate(LF, fcol = 5, lcol = 205, bins, new_bins, LengthOnly = FALSE)

# Chekcing the data by making two plots
bins <- new_bins # use the new bins
nbins <- length(bins)
fcol = 5
lcol = 4 + length(bins)
save_dir=directory

mmd <- LF1[,c(1,3:lcol)] # mmd is the input data for the clustering analysis - it should have year, lat, lon, and bin numbers

# setting up input data frames for clustering algorithm
temp = packbylatlon(mmd, 5, 5, nbins)  # aggregate the input LF across time for each grid cell

packedmmd3 = temp$table1
packedpdf3 = topdf(packedmmd3, 4, 3 + nbins)
# packedcdf3 = tocdf(packedpdf3, 4, 3 + nbins)
mmdt = packedmmd3[packedmmd3[, 4 + nbins] > 0,]
mmdtpdf = packedpdf3[packedmmd3[, 4 + nbins] > 0,] # PDF sums to 1 for each grid
mmdtpdf[, 4 + nbins] = mmdt[, 4 + nbins]
# mmdtcdf = packedcdf3[packedmmd3[, 4 + nbins] > 0,]
# mmdtcdf[, 4 + nbins] = mmdt[, 4 + nbins]

names(mmdtpdf)[2:3] <- c("lat", "lon")
mmdtpdf$Number <- 1:nrow(mmdtpdf)

rrs = sign(mmdt[, 4 + nbins]) # sample size

densmatx = matrix(0, nrow(mmdt), 420)
densmaty = matrix(0, nrow(mmdt), 420)
for(i in 1:nrow(mmdt)) {
  weightvec = t(mmdt[i, 4:(3 + nbins)])
  weightvec = weightvec / sum(weightvec)
  tempmmd = density(seq(0.61, 2.00, 0.01), weights = weightvec, bw = 0.05)
  ii <- which(tempmmd$x>=0.61 & tempmmd$x<=2.00)
  densmatx[i, ] = tempmmd$x[ii]
  densmaty[i, ] = tempmmd$y[ii]
}

cluster <- read.csv("D:/OneDrive - IATTC/IATTC/2024/Irregular clustering/YFT DEL/cluster_YFT.csv")

colcol = rep(c(2, 3, 4, 5, 6, 7, 8), 3)
densy_df <- clusthistd3(kk = 3, colseq = cluster$area + 1, colcol, plot = FALSE)

f1 <- ggplot(data = densy_df) +
  geom_line(aes(x = Length, y = Density, color = Cell), linewidth = 2) +
  theme_bw() +
  ggtitle("Offshore")

wmap <- map_data("world")
f2 <- ggplot(data = cluster) +
  geom_tile(aes(x = lon, y = lat, fill = factor(area)), color = "black") +
  geom_polygon(
    data = wmap,
    aes(long, lat, group = group),
    fill = "black",
    colour = "white",
    lwd = 0.5
  ) +
  coord_quickmap(ylim = c(min(cluster$lat), max(cluster$lat)), xlim = c(min(cluster$lon), max(cluster$lon))) +
  theme_bw()

cluster <- read.csv("D:/OneDrive - IATTC/IATTC/2024/Irregular clustering/YFT DEL/cluster_YFT2 - Copy.csv")

colcol = rep(c(2, 3, 4, 5, 6, 7, 8), 3)
densy_df <- clusthistd3(kk = 4, colseq = cluster$area, colcol, plot = FALSE)

ggplot(data = densy_df) +
  geom_line(aes(x = Length, y = Density, color = Cell), linewidth = 2) +
  theme_bw()

wmap <- map_data("world")
ggplot(data = cluster) +
  geom_tile(aes(x = lon, y = lat, fill = factor(area-1)), color = "black") +
  geom_polygon(
    data = wmap,
    aes(long, lat, group = group),
    fill = "black",
    colour = "white",
    lwd = 0.5
  ) +
  coord_quickmap(ylim = c(min(cluster$lat), max(cluster$lat)), xlim = c(min(cluster$lon), max(cluster$lon))) +
  theme_bw()
