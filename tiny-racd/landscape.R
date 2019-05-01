rm(list=ls());gc()

library(here)
library(tidyverse)
library(spatstat)
library(viridis)

library(foreach)
library(doSNOW)

nd <- 20 # num dwellings
na <- 40 # num aquatic habitats

xy_d <- rpoispp(lambda = nd,win = owin(c(0,1),c(0,1)))
xy_a <- rpoispp(lambda = na,win = owin(c(0,1),c(0,1)))

# # diagnostic plot
# plot(xy_d,pch=17,cols = adjustcolor("firebrick3",alpha.f = 0.8),main = "Dwelling/Habitat Surface")
# plot(xy_a,pch=16,cols = adjustcolor("steelblue",alpha.f = 0.8),add = T)
# legend(x = "topleft",bg = "transparent",
#        col = c("firebrick3","steelblue"),pch = 17:16,legend = c("Dwelling","Habitat"))

dist_xy <- crossdist(X = xy_a,Y = xy_d)
sigma_a <- apply(dist_xy,1,min)

aqua_df <- data.frame(x=xy_a$x,y=xy_a$y,sigma=sigma_a)
dwell_df <- data.frame(x=xy_d$x,y=xy_d$y)

# the risk on houses
psi_d <- foreach(xy = iter(dwell_df,by="row"),.combine = "rbind",.inorder = TRUE) %:%
  foreach(hab = iter(aqua_df,by = "row"),.combine = "+") %do% {
    dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
    psi <- dnorm(dist,mean=0,sd=hab$sigma)
    psi
  }

dwell_df$psi <- psi_d[,1]/sum(psi_d[,1])

# the risk surface
grid <- expand.grid(x=seq(-0.2,1.2,by=0.015),y=seq(-0.2,1.2,by=0.015))

cl <- makeCluster(4)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(grid)*nrow(aqua_df), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

surface <- foreach(xy = iter(grid,by="row"),.combine = "rbind",.inorder = TRUE) %:%
            foreach(hab = iter(aqua_df,by = "row"),.combine = "+", .options.snow = opts) %dopar% {
              dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
              psi <- dnorm(dist,mean=0,sd=hab$sigma)
              psi
            }


stopCluster(cl)
rm(cl);gc()

surface <- as.data.frame(surface)

grid$psi <- surface$V1

# the surface
ggplot() +
  geom_raster(aes(x=x,y=y,fill=psi),data=grid) +
  stat_contour(aes(x=x,y=y,z=psi),colour=grey(0.5,0.5),size=0.25,data = grid,geom = "contour") +
  geom_point(aes(x=x,y=y),shape=17,colour=grey(0.75,0.75),data=dwell_df) +
  geom_point(aes(x=x,y=y),shape=16,colour=grey(0.75,0.75),data=aqua_df) +
  scale_fill_viridis() +
  geom_contour() +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
