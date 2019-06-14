library(here)
library(geosphere)
library(readr)
library(reshape2)
library(ggplot2)
library(pbmcapply)

spat_ER <- readr::read_csv(here::here("GR_2017_updates_EA_joined_bw.csv"))

pops <- reshape2::melt(spat_ER[,c("hh_pop","ea_no")],id.vars="ea_no",measure.vars="hh_pop")

# ggplot(data = pops,aes(value,color=as.factor(ea_no))) +
#   geom_freqpoly(alpha=0.35) +
#   guides(color=FALSE) +
#   theme_bw()

ggplot(data = pops,aes(value,fill=as.factor(ea_no))) +
  geom_histogram(alpha=0.5) +
  guides(fill=FALSE) +
  theme_bw()

library(spatstat)
library(spatstat.local)
library(splines)

ea_no <- unique(spat_ER$ea_no)
distances <- pbmclapply(ea_no,function(ea){
  coords <- as.matrix(spat_ER[spat_ER$ea_no==ea,c("longitude","latitude")])
  dmat <- outer(X = 1:nrow(coords),Y = 1:nrow(coords),FUN = function(x,y){
    geosphere::distVincentyEllipsoid(coords[x,],coords[y,])
  })
  dmat
})

sites <- solapply(X = ea_no,function(ea){
  coords <- spat_ER[spat_ER$ea_no==ea,c("longitude","latitude")]
  coordinates(coords) <- c("longitude", "latitude")
  proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
  proj_coords <- spTransform(coords, CRS("+proj=utm +zone=30 ellps=WGS84"))
  xrange <- c(
    min(proj_coords@coords[,"latitude"])-10,
    max(proj_coords@coords[,"latitude"])+10
  )
  yrange <- c(
    min(proj_coords@coords[,"longitude"])-10,
    max(proj_coords@coords[,"longitude"])+10
  )
  ppp(x = proj_coords@coords[,"latitude"],y = proj_coords@coords[,"longitude"],xrange=xrange,yrange=yrange)
})

# try and estimate LGCP on multiple replicates
sites_dat <- hyperframe(Points = sites)
fit_trend <- mppm(Points ~ 1,data = sites_dat)
