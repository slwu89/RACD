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

library(here)
library(spatstat)
library(parallel)
library(foreach)
library(iterators)
library(Rcpp)
library(RcppProgress)

set.seed(50)


###############################################################################
# generate landscapes
###############################################################################

landscape <- vector("list",3)
names(landscape) <- c("clustered","CSR","regular")

# bounding box is 2000 meters; standardize to this
bbox <- 2e3

# simulate landscapes
dwellhab <- 5 # number of houses/habitat
n_dwelling <- 1e2
n_hab <- n_dwelling/dwellhab
N <- n_dwelling * 5 # number of people

# clustered (attraction)
lscape_clust_d <- spatstat::rMatClust(kappa = 15,scale = 50/bbox,mu = n_dwelling/15,win = spatstat::square(r = 1))
while(lscape_clust_d$n != n_dwelling){
  lscape_clust_d <- spatstat::rMatClust(kappa = 15,scale = 50/bbox,mu = n_dwelling/15,win = spatstat::square(r = 1))
}
lscape_clust_h <- spatstat::rpoint(n = n_hab,win = spatstat::square(r = 1))

dist_xy <- spatstat::crossdist(X = lscape_clust_h,Y = lscape_clust_d)
sigma_a <- apply(dist_xy,1,min)

landscape$clustered$dwellings <- data.frame(x=lscape_clust_d$x,y=lscape_clust_d$y)
landscape$clustered$habitats <- data.frame(x=lscape_clust_h$x,y=lscape_clust_h$y,sigma=sigma_a)

psi_d <- foreach(xy = iter(landscape$clustered$dwellings,by="row"),.combine = "rbind",.inorder = TRUE) %:%
  foreach(hab = iter(landscape$clustered$habitats,by = "row"),.combine = "+") %do% {
    dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
    psi <- dnorm(dist,mean=0,sd=hab$sigma)
    psi
  }

landscape$clustered$dwellings$psi <- psi_d[,1]/sum(psi_d[,1])

# plot(clusterfield(kppm(lscape_clust_d,~x,"MatClust"),locations = lscape_clust_d),main=paste0(lscape_clust_d$n," dwellings: clustered process"))
# points(lscape_clust_d,col="white")
# points(lscape_clust_h,pch=17,col="white")

# CSR (Poisson)
lscape_csr_d <- spatstat::rpoint(n = n_dwelling,win = spatstat::square(r = 1))
lscape_csr_h <- spatstat::rpoint(n = n_hab,win = spatstat::square(r = 1))

dist_xy <- spatstat::crossdist(X = lscape_csr_h,Y = lscape_csr_d)
sigma_a <- apply(dist_xy,1,min)

landscape$CSR$dwellings <- data.frame(x=lscape_csr_d$x,y=lscape_csr_d$y)
landscape$CSR$habitats <- data.frame(x=lscape_csr_h$x,y=lscape_csr_h$y,sigma=sigma_a)

psi_d <- foreach(xy = iter(landscape$CSR$dwellings,by="row"),.combine = "rbind",.inorder = TRUE) %:%
  foreach(hab = iter(landscape$CSR$habitats,by = "row"),.combine = "+") %do% {
    dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
    psi <- dnorm(dist,mean=0,sd=hab$sigma)
    psi
  }

landscape$CSR$dwellings$psi <- psi_d[,1]/sum(psi_d[,1])

# regular (repulsion)
mod <- spatstat::rmhmodel(cif="hardcore",par=list(beta=n_dwelling,hc=100/bbox),w=spatstat::square(r = 1))
cont <- spatstat::rmhcontrol(p=1,nrep=1e6)

lscape_reg_d <- spatstat::rmh(model=mod,start=list(n.start=n_dwelling), control=cont)
lscape_reg_h <- spatstat::rpoint(n = n_hab,win = spatstat::square(r = 1))

dist_xy <- spatstat::crossdist(X = lscape_reg_h,Y = lscape_reg_d)
sigma_a <- apply(dist_xy,1,min)

landscape$regular$dwellings <- data.frame(x=lscape_reg_d$x,y=lscape_reg_d$y)
landscape$regular$habitats <- data.frame(x=lscape_reg_h$x,y=lscape_reg_h$y,sigma=sigma_a)

psi_d <- foreach(xy = iter(landscape$regular$dwellings,by="row"),.combine = "rbind",.inorder = TRUE) %:%
  foreach(hab = iter(landscape$regular$habitats,by = "row"),.combine = "+") %do% {
    dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
    psi <- dnorm(dist,mean=0,sd=hab$sigma)
    psi
  }

landscape$regular$dwellings$psi <- psi_d[,1]/sum(psi_d[,1])


###############################################################################
# ranges for other parameters
###############################################################################

EIR_vals <- c(0.001,0.003,0.005)
interventions <- c(control=-1,racd=2,rfmda=0,rfmda_rfvc=6,rfvc=1)
p_range <- seq(0.2,0.8,by=0.1)
import_range <- (1/365)*c(1,2,5,10,20)


###############################################################################
# load files we need
###############################################################################

source(here::here("racd-setup.R"))
sourceCpp(here::here("intervention-src/main.cpp"))


###############################################################################
# run the simulations
###############################################################################


run_names <- NULL
i <- 1

for(l in 1:length(landscape)){

  for(e in 1:length(EIR_vals)){

    EIR <- EIR_vals[e]
    lscape <- landscape[[l]]

    # calculate equilibrium
    RACD_init <- RACD_Setup(N = N,EIR_mean = EIR,xy_d = lscape$dwellings,xy_a = lscape$habitats,theta = RACD_theta)

    for(i in 1:length(interventions)){

      for(im in 1:length(import_range)){

        for(p_ix in 1:length(p_range)){

          for(p_n in 1:length(p_range)){

            cat(" --- factorial cell ",i," of 11025 --- \n")

            this_run <- paste0("landscape-",names(landscape)[l],"_EIR-",EIR_vals[e],"_intervention-",names(interventions)[i],
                               "_import-",round(import_range[im],3),"_pi-",p_range[p_ix],"_pn-",p_range[p_n]
                               )
            run_names <- c(run_names,this_run)

            i <- i + 1

          }

        }

      }

    }

  }
}

# how to do parallel foreach

cl <- makeCluster(4)
registerDoSNOW(cl)

clusterEvalQ(cl,{
  Rcpp::sourceCpp(here::here("gillespie-SIR-CXX/gillespie-sim.cpp"))
})

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
