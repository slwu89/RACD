rm(list=ls());gc()

library(RACD)
library(tidyverse)
library(spatstat)
library(viridis)

xy_h <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
xy_b <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
theta <- RACD_Parameters()
init <- RACD_Setup(as.matrix.ppx(xy_h),as.matrix.ppx(xy_b),theta)

x=sapply(init$houses,function(x){x$x})
y=sapply(init$houses,function(x){x$y})
z=sapply(init$houses,function(x){x$psi})
x_a = sapply(init$breedingSites,function(x){x$x})
y_a = sapply(init$breedingSites,function(x){x$y})
z_a = sapply(init$breedingSites,function(x){x$sigma})

habitats <- data.frame(x=x_a,y=y_a,sigma=z_a)

library(foreach)
library(doParallel)

grid <- expand.grid(x=seq(-0.2,1.2,by=0.015),y=seq(-0.2,1.2,by=0.015))

cl <- makeCluster(4)
registerDoParallel(cl)

surface <- foreach(xy = iter(grid,by="row"),.combine = "rbind",.inorder = TRUE) %:%
            foreach(hab = iter(habitats,by = "row"),.combine = "+") %dopar% {
              dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
              psi <- dnorm(dist,mean=0,sd=hab$sigma)
              psi
            }


stopCluster(cl)

surface <- as.data.frame(surface)

grid$psi <- surface$V1

houses <- data.frame(x=x,y=y,psi=z)
habitats <- data.frame(x=x_a,y=y_a,sigma=z_a)

# houses
ggplot() +
  geom_raster(aes(x=x,y=y,fill=psi),data=grid) +
  stat_contour(aes(x=x,y=y,z=psi),colour=grey(0.75,0.75),size=0.25,data = grid,geom = "contour") +
  geom_point(aes(x=x,y=y,colour=psi),data=houses) +
  scale_color_viridis(begin=0.5,end=1) +
  geom_contour() +
  theme_bw()

# habitats
ggplot() +
  geom_raster(aes(x=x,y=y,fill=psi),data=grid) +
  stat_contour(aes(x=x,y=y,z=psi),colour=grey(0.75,0.75),size=0.25,data = grid,geom = "contour") +
  geom_point(aes(x=x,y=y,colour=sigma),data=habitats) +
  scale_color_viridis(begin=1,end=0.5) +
  geom_contour() +
  theme_bw()

houseXhabitat <- rbind(data.frame(unname(habitats),type="habitat"),
                       data.frame(unname(houses),type="house"))
names(houseXhabitat)[1:3] <- c("x","y","z")

ggplot() +
  geom_raster(aes(x=x,y=y,fill=psi),data=grid) +
  stat_contour(aes(x=x,y=y,z=psi),colour=grey(0.75,0.75),size=0.25,data = grid,geom = "contour") +
  geom_point(aes(x=x,y=y,shape=type),colour=grey(0.75,0.75),data=houseXhabitat) +
  scale_fill_viridis() +
  geom_contour() +
  theme_bw()
