rm(list=ls());gc()

library(here)
library(tidyverse)
library(spatstat)
library(viridis)

library(foreach)
library(doSNOW)

nd <- 100 # num dwellings
na <- 20 # num aquatic habitats

# xy_d <- rpoispp(lambda = nd,win = owin(c(0,1),c(0,1)))
# xy_a <- rpoispp(lambda = na,win = owin(c(0,1),c(0,1)))

rmh_window <- owin(c(0,1),c(0,1))
rmh_mod <- rmhmodel(cif="hardcore",par = list(beta=nd*na,hc=0.05),
                    w = rmh_window,types = c("dwelling","aquatic"))
rmh_con <- rmhcontrol(p = 1,nverb = 1e4,fixall = TRUE)
rmh_start <- rmhstart(n.start = c(dwelling=nd,aquatic=na))
xy_rmh <- rmh(model = rmh_mod,start = rmh_start,control = rmh_con)

xy_a <- data.frame(x=xy_rmh$x[xy_rmh$marks=="aquatic"],y=xy_rmh$y[xy_rmh$marks=="aquatic"])
xy_d <- data.frame(x=xy_rmh$x[xy_rmh$marks=="dwelling"],y=xy_rmh$y[xy_rmh$marks=="dwelling"])

# # diagnostic plot
# plot(xy_d,pch=17,cols = adjustcolor("firebrick3",alpha.f = 0.8),main = "Dwelling/Habitat Surface")
# plot(xy_a,pch=16,cols = adjustcolor("steelblue",alpha.f = 0.8),add = T)
# legend(x = "topleft",bg = "transparent",
#        col = c("firebrick3","steelblue"),pch = 17:16,legend = c("Dwelling","Habitat"))

dist_xy <- crossdist(X = as.ppp(xy_a,W = rmh_window),Y = as.ppp(xy_d,W = rmh_window))
sigma_a <- apply(dist_xy,1,min)

dmat_dwell <- crossdist(X = as.ppp(xy_d,W = rmh_window),Y = as.ppp(xy_d,W = rmh_window))

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
grid_res <- 0.015
grid <- expand.grid(x=seq(0-grid_res,1+grid_res,by=grid_res),y=seq(0-grid_res,1+grid_res,by=grid_res))

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
  geom_point(aes(x=x,y=y,size=cut(psi,breaks=quantile(psi),include.lowest = T)),shape=17,colour=grey(0.9,0.9),data=dwell_df) +
  geom_point(aes(x=x,y=y),shape=16,colour=grey(0.75,0.75),size=1.85,data=aqua_df) +
  scale_fill_viridis() +
  geom_contour() +
  theme_bw() +
  guides(size = FALSE) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
