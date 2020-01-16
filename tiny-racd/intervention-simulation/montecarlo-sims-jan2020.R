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
library(tidyverse)
library(reshape2)

library(spatstat)
library(viridis)

library(foreach)
library(doSNOW)

library(Rcpp)
library(RcppProgress)

set.seed(50)


###############################################################################
# generate landscape
###############################################################################

# bounding box is 2000 meters; standardize to this
bbox <- 2e3

# simulate landscapes
dwellhab <- 5 # number of houses/habitat
n_dwelling <- 2e2
n_hab <- n_dwelling/dwellhab
N <- n_dwelling * 5 # number of people

# CSR (Poisson)
lscape_csr_d <- spatstat::rpoint(n = n_dwelling,win = spatstat::square(r = 1))
lscape_csr_h <- spatstat::rpoint(n = n_hab,win = spatstat::square(r = 1))

dist_xy <- spatstat::crossdist(X = lscape_csr_h,Y = lscape_csr_d)
sigma_a <- apply(dist_xy,1,min)

landscape <- list()
landscape$dwellings <- data.frame(x=lscape_csr_d$x,y=lscape_csr_d$y)
landscape$habitats <- data.frame(x=lscape_csr_h$x,y=lscape_csr_h$y,sigma=sigma_a)

psi_d <- foreach(xy = iter(landscape$dwellings,by="row"),.combine = "rbind",.inorder = TRUE) %:%
  foreach(hab = iter(landscape$habitats,by = "row"),.combine = "+") %do% {
    dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
    psi <- dnorm(dist,mean=0,sd=hab$sigma)
    psi
  }

landscape$dwellings$psi <- psi_d[,1]/sum(psi_d[,1])

# # compute risk surface
# grid_res <- 1/50
# grid <- expand.grid(x=seq(0-grid_res,1+grid_res,by=grid_res),y=seq(0-grid_res,1+grid_res,by=grid_res))
# 
# cl <- parallel::makeCluster(4)
# doSNOW::registerDoSNOW(cl)
# 
# pb <- txtProgressBar(max = nrow(grid)*nrow(landscape$habitats), style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# landscape$psi_surface <- foreach(xy = iter(grid,by="row"),.combine = "rbind",.inorder = TRUE) %:%
#   foreach(hab = iter(landscape$habitats,by = "row"),.combine = "+", .options.snow = opts) %dopar% {
#     dist <- as.matrix(dist(x = rbind(as.vector(xy),c(hab$x,hab$y))))[1,2]
#     psi <- dnorm(dist,mean=0,sd=hab$sigma)
#     psi
#   }
# 
# stopCluster(cl)
# rm(cl);gc()
# 
# # plot the risk surface
# surface <- as.data.frame(landscape$psi_surface)
# 
# grid$psi <- surface$V1
# 
# # the surface
# gg_psi_CSR <- ggplot() +
#   geom_raster(aes(x=x,y=y,fill=psi),data=grid) +
#   stat_contour(aes(x=x,y=y,z=psi),colour=grey(0.5,0.5),size=0.25,data = grid,geom = "contour") +
#   geom_point(aes(x=x,y=y,size=cut(psi,breaks=quantile(psi),include.lowest = T)),shape=17,colour=grey(0.9,0.9),data=landscape$dwellings) +
#   geom_point(aes(x=x,y=y),shape=16,colour=grey(0.75,0.75),size=1.85,data=landscape$habitats) +
#   scale_fill_viridis() +
#   geom_contour() +
#   theme_bw() +
#   guides(size = FALSE) +
#   theme(axis.title.x=element_blank(),axis.title.y=element_blank())
# 
# ggsave(filename = here::here("graphics/psi_surface_CSR_jan2020.pdf"),plot = gg_psi_CSR,dpi = 320,width = 10,height = 8)



###############################################################################
#
# MC SIMS
#
###############################################################################

out_folder <- "/Users/slwu89/Dropbox/racd_out/jan2020/"
if(!dir.exists(out_folder)){
  dir.create(out_folder)
}


interventions <- c(control=-1,racd=2,rfmda=0,rfmda_rfvc=6,rfvc=1)

# the range of transmission intensity in which we'd consider RACD
EIR_range <- 1:10/365

# importation from 0:100 cases/1e3 person-yr
import_range <- seq(0,100,by=10)/365

# coverage is fixed for now.
p_ix <- p_n <- 0.75

# run the sims
source(here::here("racd-setup.R"))
dmat_dwell <- crossdist(X = as.ppp(landscape$dwellings[,1:2],square(r=1)),Y = as.ppp(landscape$dwellings[,1:2],square(r=1)))

# setup the cluster
cl <- parallel::makeCluster(5)
parallel::clusterSetRNGStream(cl = cl,iseed = 95812142L)

parallel::clusterEvalQ(cl,{
  Rcpp::sourceCpp(here::here("intervention-src/main.cpp"),rebuild = TRUE)
})

# run sims
i=1
niter <- length(EIR_range)*length(import_range)*length(interventions)
for(eir_i in seq_along(EIR_range)){
  for(import_i in seq_along(import_range)){
    
    eir <- EIR_range[eir_i]
    import <- import_range[import_i]
    
    # calculate equilibrium
    RACD_init <- RACD_Setup(N = N,EIR_mean = eir,xy_d = landscape$dwellings,xy_a = landscape$habitats,theta = RACD_theta)
    
    # calculate immunity parameters for imported cases
    imm_import <- imported_immune(EIR = eir,theta = RACD_theta)
    RACD_theta <- c(RACD_theta,imm_import,import_rate=import)
    
    for(int_i in seq_along(interventions)){

      int <- interventions[int_i]
      mc_iter <- paste0("SIM_intervention-",names(interventions)[int_i],"_eir-",round(eir,3),"_import-",round(import,3))

      parallel::clusterExport(cl = cl,varlist = c("RACD_init","RACD_theta","dmat_dwell","int","p_ix","p_n"))

      # run simulations
      # t_start <- Sys.time()
      mc_reps_cores <- parallel::clusterEvalQ(cl,{
        mc_reps <- vector("list",20)
        for(i in 1:20){
          mc_reps[[i]] <- tiny_racd(humans_param = RACD_init$humans,
                                    house_param = RACD_init$houses,
                                    mosy_param = RACD_init$mosy,
                                    theta = RACD_theta,
                                    tmax = (365*5)+100,
                                    int_type = int,
                                    tstart = 101,
                                    tend = (365*5)+100,
                                    tdelay = 30,
                                    dmat = dmat_dwell,
                                    radius = 500/2e3,
                                    p_index = p_ix,
                                    p_neighbor = p_n,
                                    prog_bar = FALSE)
        }
        return(mc_reps)
      })
      # t_end <- Sys.time()

      saveRDS(object = mc_reps_cores,file = paste0(out_folder,mc_iter,".rds"))
      rm(mc_reps_cores);gc()

      cat(" --- finished intervention cell ",i," of ",length(interventions)," --- \n")
      i = i + 1

    }
  }
}













###############################################################################
#
# ANALYSIS OF MC SIMULATION
#
###############################################################################

library(data.table)

###############################################################################
# clustering
###############################################################################

# read in the files
files <- list.files(path = out_folder)
files_2read <- which(files %in% paste0("landscape-clust_intervention-",names(interventions),".rds"))

out_clust <- vector("list",length(files_2read))
names(out_clust) <- sub(pattern = ".rds",replacement = "",x = files[files_2read])

for(i in 1:length(files_2read)){
  out_clust[[i]] <- read_rds(paste0(out_folder,files[files_2read[i]]))
  out_clust[[i]] <- do.call(c,out_clust[[i]])
}


stateage_clust <- vector("list",length(out_clust))
names(stateage_clust) <- sapply(strsplit(names(out_clust),"_"),function(x){
  if(length(x)==2){
    return(
      sapply(strsplit(x[2],"-"),"[[",2)
    )
  } else {
    return(
      paste0(sapply(strsplit(x[2],"-"),"[[",2),"-",x[3])
    )
  }
})

for(i in 1:length(out_clust)){
  stateage_clust[[i]] <- parallel::mclapply(X = out_clust[[i]],FUN = function(x){
    x_sa <- reshape2::melt(x$state_hist$state_age)
    colnames(x_sa) <- c("state","age","count","day")
    return(x_sa)
  })
  stateage_clust[[i]] <- reshape2::melt(stateage_clust[[i]],id.vars=c("state","age","count","day"))
  colnames(stateage_clust[[i]])[5] <- "run"
  stateage_clust[[i]] <- as.data.table(stateage_clust[[i]])
}

stateage_clust <- parallel::mclapply(X = stateage_clust,FUN = function(dt){
  dt[,
     j=.(hi = quantile(count,0.975),lo = quantile(count,0.025),mean = mean(count)),
     by = .(state,age,day)]
})

stateage_clust_mt <- rbindlist(stateage_clust,idcol="intervention")

gg_stateage_intervention_clust <- ggplot(data = stateage_clust_mt) +
  geom_line(aes(x=day,y=mean,color=state),alpha=0.75) +
  geom_ribbon(aes(x=day,ymin=lo,ymax=hi,fill=state),alpha=0.5) +
  facet_grid(age ~ intervention,scales = "free_y") +
  theme_bw()

ggsave(filename = here::here("graphics/state-age-int-cluster.pdf"),plot = gg_stateage_intervention_clust,dpi = 320,width = 12,height = 10)


###############################################################################
# CSR
###############################################################################

# read in the files
files <- list.files(path = out_folder)
files_2read <- which(files %in% paste0("landscape-CSR_intervention-",names(interventions),".rds"))

out_CSR <- vector("list",length(files_2read))
names(out_CSR) <- sub(pattern = ".rds",replacement = "",x = files[files_2read])

for(i in 1:length(files_2read)){
  out_CSR[[i]] <- read_rds(paste0(out_folder,files[files_2read[i]]))
  out_CSR[[i]] <- do.call(c,out_CSR[[i]])
}


stateage_CSR <- vector("list",length(out_CSR))
names(stateage_CSR) <- sapply(strsplit(names(out_CSR),"_"),function(x){
  if(length(x)==2){
    return(
      sapply(strsplit(x[2],"-"),"[[",2)
    )
  } else {
    return(
      paste0(sapply(strsplit(x[2],"-"),"[[",2),"-",x[3])
    )
  }
})

for(i in 1:length(out_CSR)){
  stateage_CSR[[i]] <- parallel::mclapply(X = out_CSR[[i]],FUN = function(x){
    x_sa <- reshape2::melt(x$state_hist$state_age)
    colnames(x_sa) <- c("state","age","count","day")
    return(x_sa)
  })
  stateage_CSR[[i]] <- reshape2::melt(stateage_CSR[[i]],id.vars=c("state","age","count","day"))
  colnames(stateage_CSR[[i]])[5] <- "run"
  stateage_CSR[[i]] <- as.data.table(stateage_CSR[[i]])
}

stateage_CSR <- parallel::mclapply(X = stateage_CSR,FUN = function(dt){
  dt[,
     j=.(hi = quantile(count,0.975),lo = quantile(count,0.025),mean = mean(count)),
     by = .(state,age,day)]
})

stateage_CSR_mt <- rbindlist(stateage_CSR,idcol="intervention")

gg_stateage_intervention_CSR <- ggplot(data = stateage_CSR_mt) +
  geom_line(aes(x=day,y=mean,color=state),alpha=0.75) +
  geom_ribbon(aes(x=day,ymin=lo,ymax=hi,fill=state),alpha=0.5) +
  facet_grid(age ~ intervention,scales = "free_y") +
  theme_bw()

ggsave(filename = here::here("graphics/state-age-int-CSR.pdf"),plot = gg_stateage_intervention_CSR,dpi = 320,width = 12,height = 10)


###############################################################################
# regular
###############################################################################

# read in the files
files <- list.files(path = out_folder)
files_2read <- which(files %in% paste0("landscape-regular_intervention-",names(interventions),".rds"))

out_regular <- vector("list",length(files_2read))
names(out_regular) <- sub(pattern = ".rds",replacement = "",x = files[files_2read])

for(i in 1:length(files_2read)){
  out_regular[[i]] <- read_rds(paste0(out_folder,files[files_2read[i]]))
  out_regular[[i]] <- do.call(c,out_regular[[i]])
}

# time series
stateage_regular <- vector("list",length(out_regular))
names(stateage_regular) <- sapply(strsplit(names(out_regular),"_"),function(x){
  if(length(x)==2){
    return(
      sapply(strsplit(x[2],"-"),"[[",2)
    )
  } else {
    return(
      paste0(sapply(strsplit(x[2],"-"),"[[",2),"-",x[3])
    )
  }
})

for(i in 1:length(out_regular)){
  stateage_regular[[i]] <- parallel::mclapply(X = out_regular[[i]],FUN = function(x){
    x_sa <- reshape2::melt(x$state_hist$state_age)
    colnames(x_sa) <- c("state","age","count","day")
    return(x_sa)
  })
  stateage_regular[[i]] <- reshape2::melt(stateage_regular[[i]],id.vars=c("state","age","count","day"))
  colnames(stateage_regular[[i]])[5] <- "run"
  stateage_regular[[i]] <- as.data.table(stateage_regular[[i]])
}

stateage_regular <- parallel::mclapply(X = stateage_regular,FUN = function(dt){
  dt[,
     j=.(hi = quantile(count,0.975),lo = quantile(count,0.025),mean = mean(count)),
     by = .(state,age,day)]
})

stateage_regular_mt <- rbindlist(stateage_regular,idcol="intervention")

gg_stateage_intervention_regular <- ggplot(data = stateage_regular_mt) +
  geom_line(aes(x=day,y=mean,color=state),alpha=0.75) +
  geom_ribbon(aes(x=day,ymin=lo,ymax=hi,fill=state),alpha=0.5) +
  facet_grid(age ~ intervention,scales = "free_y") +
  theme_bw()

ggsave(filename = here::here("graphics/state-age-int-regular.pdf"),plot = gg_stateage_intervention_regular,dpi = 320,width = 12,height = 10)


###############################################################################
#
# ANALYSIS OF MC SIMULATION: CUMULATIVE INCIDENCE
#
###############################################################################

# clustered
cum_inc_clust <- lapply(X = out_clust,FUN = function(int_lst){
  sapply(int_lst,function(df){
    colSums(df$state_hist$clinical_incidence)
  })
})
names(cum_inc_clust) <- names(stateage_regular)

cum_inc_clust <- lapply(cum_inc_clust,function(mat){
  mat_sum <- apply(X = mat,FUN = function(x){
    quantile(x = x,probs = c(0.025,0.975))
  },MARGIN = 1)
  mat_sum <- rbind(mat_sum,rowMeans(mat))
  rownames(mat_sum)[3] <- "mean"
  as.data.table(reshape2::melt(t(mat_sum)))
})

cum_inc_clust <- rbindlist(cum_inc_clust,idcol = "intervention")

# CSR
cum_inc_CSR <- lapply(X = out_CSR,FUN = function(int_lst){
  sapply(int_lst,function(df){
    colSums(df$state_hist$clinical_incidence)
  })
})
names(cum_inc_CSR) <- names(stateage_regular)

cum_inc_CSR <- lapply(cum_inc_CSR,function(mat){
  mat_sum <- apply(X = mat,FUN = function(x){
    quantile(x = x,probs = c(0.025,0.975))
  },MARGIN = 1)
  mat_sum <- rbind(mat_sum,rowMeans(mat))
  rownames(mat_sum)[3] <- "mean"
  as.data.table(reshape2::melt(t(mat_sum)))
})

cum_inc_CSR <- rbindlist(cum_inc_CSR,idcol = "intervention")

# regular
cum_inc_regular <- lapply(X = out_regular,FUN = function(int_lst){
  sapply(int_lst,function(df){
    colSums(df$state_hist$clinical_incidence)
  })
})
names(cum_inc_regular) <- names(stateage_regular)

cum_inc_regular <- lapply(cum_inc_regular,function(mat){
  mat_sum <- apply(X = mat,FUN = function(x){
    quantile(x = x,probs = c(0.025,0.975))
  },MARGIN = 1)
  mat_sum <- rbind(mat_sum,rowMeans(mat))
  rownames(mat_sum)[3] <- "mean"
  as.data.table(reshape2::melt(t(mat_sum)))
})

cum_inc_regular <- rbindlist(cum_inc_regular,idcol = "intervention")

cum_inc_all <- rbindlist(list(cluster=cum_inc_clust,
                              CSR=cum_inc_CSR,
                              regular=cum_inc_regular),idcol = "landscape")

cum_inc_all <- dcast(cum_inc_all[Var1 == "all" & Var2 %in% c("mean","2.5%","97.5%"),],landscape + intervention ~ Var2,value.var = "value")
colnames(cum_inc_all)[3:4] <- c("lo",'hi')

gg_cum_inc_allages <- ggplot(data = cum_inc_all) +
  geom_point(aes(x=landscape,y=mean,color=landscape)) +
  geom_errorbar(aes(x=landscape,ymin=lo,ymax=hi,color=landscape)) +
  facet_grid(. ~ intervention) +
  theme_bw() +
  guides(color = FALSE)

ggsave(filename = here::here("graphics/cuminc-allage.pdf"),plot = gg_cum_inc_allages,dpi = 320,width = 12,height = 6)



# ###############################################################################
# # ranges for other parameters
# ###############################################################################
#
# EIR_vals <- c(0.001,0.003,0.005)
# interventions <- c(racd=2,rfmda=0,rfmda_rfvc=6,rfvc=1)
# p_range <- seq(0.2,0.8,by=0.1)
# import_range <- (1/365)*c(1,2,5,10,20)
#
#
# ###############################################################################
# # load files we need
# ###############################################################################
#
# source(here::here("racd-setup.R"))
# # Rcpp::sourceCpp(here::here("intervention-src/main.cpp"),rebuild = TRUE)
#
# ###############################################################################
# # run the simulations
# ###############################################################################
#
# ii <- 1
#
# n_iter <- length(landscape)*length(EIR_vals)*length(interventions)*(length(p_range)^2)*length(import_range)
# run_names <- character(n_iter)
#
# cl <- makeCluster(4)
# clusterSetRNGStream(cl = cl,iseed = 42L)
#
# clusterEvalQ(cl,{
#   Rcpp::sourceCpp(here::here("intervention-src/main.cpp"),rebuild = TRUE)
# })
#
# # # begin factorial sweep
# # l=1
# # e=1
# # i=1
# # im=1
# # p_ix=1
# # p_n=1
#
# for(l in 1:length(landscape)){
#
#   for(e in 1:length(EIR_vals)){
#
#     EIR <- EIR_vals[e]
#     lscape <- landscape[[l]]
#
#     dmat_dwell <- crossdist(X = as.ppp(lscape$dwellings[,1:2],square(r=1)),Y = as.ppp(lscape$dwellings[,1:2],square(r=1)))
#
#     # calculate equilibrium
#     RACD_init <- RACD_Setup(N = N,EIR_mean = EIR,xy_d = lscape$dwellings,xy_a = lscape$habitats,theta = RACD_theta)
#
#     # calculate immunity parameters for imported cases
#     imm_import <- imported_immune(EIR = EIR,theta = RACD_theta)
#     RACD_theta <- c(RACD_theta,imm_import,import_rate=NaN)
#
#     for(i in 1:length(interventions)){
#
#       # for(im in 1:length(import_range)){
#       for(im in 1:length(import_range)){
#
#         for(p_ix in 1:length(p_range)){
#         # for(p_ix in 1:length(p_range)){
#
#           for(p_n in 1:length(p_range)){
#           # for(p_n in 1:length(p_range)){
#
#             this_run <- paste0("landscape-",names(landscape)[l],"_EIR-",EIR_vals[e],"_intervention-",names(interventions)[i],
#                                "_import-",round(import_range[im],3),"_pi-",p_range[p_ix],"_pn-",p_range[p_n]
#                                )
#
#             # do the MC iterations (4*6=24)
#             int_type <- interventions[i]
#             RACD_theta[["import_rate"]] <- import_range[im]
#             p_index <- p_range[p_ix]
#             p_neighbor <- p_range[p_n]
#             clusterExport(cl = cl,varlist = c("RACD_init","RACD_theta","dmat_dwell","int_type","p_index","p_neighbor"))
#
#             # run simulations
#             mc_reps_cores <- clusterEvalQ(cl,{
#               mc_reps <- vector("list",4)
#               for(i in 1:4){
#                 mc_reps[[i]] <- tiny_racd(humans_param = RACD_init$humans,
#                                           house_param = RACD_init$houses,
#                                           mosy_param = RACD_init$mosy,
#                                           theta = RACD_theta,
#                                           tmax = (365*2)+100,
#                                           int_type = int_type,
#                                           tstart = 101,
#                                           tend = 100+365,
#                                           tdelay = 60,
#                                           dmat = dmat_dwell,
#                                           radius = 500/2e3,
#                                           p_index = p_index,
#                                           p_neighbor = p_neighbor,
#                                           prog_bar = FALSE)
#               }
#               return(mc_reps)
#             })
#
#             saveRDS(object = mc_reps_cores,file = paste0("/Users/slwu89/Dropbox/racd_out/",this_run,".rds"))
#             rm(mc_reps_cores);gc()
#
#             cat(" --- finished factorial cell ",ii," of ",n_iter," --- \n")
#             ii <- ii + 1
#
#           }
#
#         }
#
#       }
#
#     }
#
#   }
# } # end factorial sweep
#
# stopCluster(cl)
# rm(cl);gc()
