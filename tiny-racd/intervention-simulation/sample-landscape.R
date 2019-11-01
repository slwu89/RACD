###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Marshall Lab (https://www.marshalllab.com)
#   Sean Wu (slwu89@berkeley.edu)
#   October 2019
#   Sample landscapes under 3 regimes for simulation
#
###############################################################################

rm(list=ls());gc()
library(spatstat)

# simulate landscapes
dwellhab <- 5 # number of houses/habitat
n <- 1e3
n_dwelling <- 1e2
n_hab <- n_dwelling/dwellhab

# clustered (attraction)

samp_clust <- spatstat::rMatClust(kappa = 10,scale = 50/1e3,mu = 10,win = spatstat::square(r = 1))
while(samp_clust$n != 100L){
  samp_clust <- spatstat::rMatClust(kappa = 10,scale = 50/1e3,mu = 10,win = spatstat::square(r = 1))
}

samp_habitats <- spatstat::rpoint(n = n_hab,win = spatstat::square(r = 1))

plot(clusterfield(kppm(samp_clust,~x,"MatClust"),locations = samp_clust),main=paste0(samp_clust$n," points in cluster pattern"))
points(samp_clust,col="white")
points(samp_habitats,pch=17,col="white")

# ddist_lst <- vector("list",1e3)
# for(i in 1:1e3){
#   samp_clust <- spatstat::rMatClust(kappa = 10,scale = 50/1e3,mu = 10,win = spatstat::square(r = 1))
#   while(samp_clust$n != 100L){
#     samp_clust <- spatstat::rMatClust(kappa = 10,scale = 50/1e3,mu = 10,win = spatstat::square(r = 1))
#   }
#   
#   ddist <- crossdist(samp_clust,samp_clust)
#   ddist_lst[[i]] <- ddist[upper.tri(ddist)]
# }
# 
# ddist_flat <- unlist(ddist_lst)
# ddist_ecdf <- ecdf(ddist_flat)
# hist(ddist_flat)
# plot(density(ddist_flat,from = 0))
# x <- seq(0,1.5,by=0.005)
# plot(x,ddist_ecdf(x),type="l")

# random (CSR)

samp_unif <- spatstat::rpoint(n = 1e2,win = spatstat::square(r = 1))

# regular (repulsion)
# we are going to simulate from the hard-core process

mod <- spatstat::rmhmodel(cif="hardcore",par=list(beta=n_dwelling/1e3,hc=50),w=spatstat::square(r = 1e3))
cont <- spatstat::rmhcontrol(p=1,nrep=1e6)

mod <- spatstat::rmhmodel(cif="hardcore",par=list(beta=n_dwelling,hc=50/1e3),w=spatstat::square(r = 1))
cont <- spatstat::rmhcontrol(p=1,nrep=1e6)

samp_hc <- spatstat::rmh(model=mod,start=list(n.start=n_dwelling), control=cont)
plot(samp_hc)


