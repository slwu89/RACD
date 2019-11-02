# source("/Users/slwu89/Desktop/git/RACD/tiny-racd/intervention-simulation/debug-script.R")
rm(list=ls());gc()

library(here)
library(spatstat)
library(foreach)
library(iterators)
library(Rcpp)
library(RcppProgress)


###############################################################################
# SETUP
###############################################################################


library(deSolve)
library(here)

Rcpp::sourceCpp("/Users/slwu89/Desktop/git/RACD/tiny-racd/racd-setup.cpp")
source("/Users/slwu89/Desktop/git/RACD/tiny-racd/mosy-equilibrium.R")
source("/Users/slwu89/Desktop/git/RACD/tiny-racd/racd-parameters.R")

# helper functions
divmod <- function(a,b){
	a <- as.integer(a)
	b <- as.integer(b)
	c(
		quo = a %/% b,
		rem = a %% b
	)
}

# n: number of "things"
# p: set of "places" they can go
distribute <- function(n,p){
	n <- as.integer(n)
	p <- as.integer(p)
	distn <- rep(0L,p)
	div <- divmod(n,p)
	for(i in 0:(p-1)){
		distn[i+1] <- div[["quo"]] + (i < div[["rem"]])
	}
	return(distn)
}

# set up immune (IB,ID,ICA)
# a: my age
# EIR: my EIR
# theta: vector of parameters
immune_initial <- function(a, EIR, theta){

  state <- c("IB"=0,"ID"=0,"ICA"=0)
  time <- seq(from=0,to=a,by=a/1e3)
	immune <- deSolve::lsoda(y = state,times = time,func = immmune_ode,parms = theta, EIR_h = EIR)

  return(c(
    IB = unname(immune[nrow(immune),"IB"]),
    ID = unname(immune[nrow(immune),"ID"]),
    ICA = unname(immune[nrow(immune),"ICA"])
  ))
}

# set up P(state)
# a: my age
# EIR: my EIR
# theta: vector of parameters
state_initial <- function(a, EIR, theta){

	state <- c(S=1, T=0, D=0, A=0, U=0, P=0, IB=0, ICA=0)
	time <- seq(from=0,to=a,by=a/1e3)
	states <- deSolve::lsoda(y = state,times = time,func = state_ode,parms = theta, EIR_h = EIR)

	return(c(
    S = unname(states[nrow(states),"S"]),
		T = unname(states[nrow(states),"T"]),
		D = unname(states[nrow(states),"D"]),
		A = unname(states[nrow(states),"A"]),
		U = unname(states[nrow(states),"U"]),
		P = unname(states[nrow(states),"P"])
  ))
}

# initial infectivity to mosquitos
c_initial <- function(human,theta){
	with(as.list(theta),{
		switch(human$state,
			S = {return(0)},
			E = {return(0)},
			T = {return(cT)},
			D = {return(cD)},
			A = {
				c = cU + (cD - cU)*(human$prDetectAMic^gammaI)
				return(c)
			},
			U = {return(cU)},
			P = {return(0)}
			)
	})
}


# calculate the IB,ID,ICA for imported cases
imported_immune <- function(EIR,theta){
  init <- immune_initial(a = 21, EIR = EIR, theta = theta)
  names(init) <- paste0(names(init),"_imp")
  return(init)
}


RACD_Setup <- function(N, EIR_mean, xy_d, xy_a, theta){


	# theta[["mu"]] <- mu <- 1/(meanAge*365) # Daily death rate as a function of mean age in years

  with(as.list(theta),{

    # calculate EIR and EIR on houses. To get EIR on people, we need pi, we'll do that later
    cat(" --- begin calculating landscape-dependent parameters --- \n")
    numA <- nrow(xy_a)
    numD <- nrow(xy_d)

    psi_house <- xy_d$psi

    EIR_tot <- EIR_mean * N
    EIR_houses <- EIR_tot * psi_house

    household_sizes <- distribute(n=N,p=numD)
    household_assignment <- rep(x = 1:numD,times = household_sizes)

    pi_vec <- rep(0,N) # pi for each person
    pi_house <- rep(0,numD) # summed pi for each house

    houses <- mapply(FUN = function(psii,nn){
      list(psi=psii,n=nn)
    },psii=psi_house,nn=household_sizes,SIMPLIFY = FALSE)


    humans <- vector("list",N)
    for (j in 1:N) {

      # basic parameters
      humans[[j]]$age <- rexp(n = 1, rate = 1/meanAge)
      humans[[j]]$zeta <- rlnorm(n = 1, meanlog = -sigma2/2, sdlog = sqrt(sigma2))
      humans[[j]]$house <- household_assignment[j]

      # interventions (none for now)
      humans[[j]]$ITN <- FALSE

      humans[[j]]$w <- 1
      humans[[j]]$y <- 1
      humans[[j]]$z <- 0
    }

    # calculate pi
    for(i in 1:numD){
      zetas <- sapply(humans[which(household_assignment == i)],function(x){x$zeta})
      ages <- sapply(humans[which(household_assignment == i)],function(x){x$age})
      pi_house[i] <- sum(zetas * (1 - rho * exp(-ages/a0)))
    }
    for(j in 1:N){
      pi_vec[j] <-  psi_house[humans[[j]]$house] * (humans[[j]]$zeta * (1 - rho * exp(-humans[[j]]$age/a0)) / pi_house[humans[[j]]$house])
    }
    cat(" --- done calculating landscape-dependent parameters --- \n")

    # if we wanted the conditional probs of landing on someone at house 1, knowing it lands at 1
    # pi_vec[which(household_assignment==1)] / psi_house[1]

    # calculate immunity and state probabilities for population
    cat("\n --- begin calculating immunity and state probabilities for population --- \n")

    # mean value of ICA for 20-year olds
    theta[["initICA20"]] <- immune_initial(a = 20,EIR = EIR_mean,theta = theta)[["ICA"]]

    pb <- txtProgressBar(min = 1, max = N, initial = 1, style=3)
    for(j in 1:N){
      # run the immune ODEs
      a <- humans[[j]]$age
      epsilon <- EIR_tot * pi_vec[j]
      init_i <- immune_initial(a = a, EIR = epsilon, theta = theta)

      # my initial immune values
      IB <- init_i[["IB"]]
      ID <- init_i[["ID"]]
      ICA <- init_i[["ICA"]]
      ICM <- theta[["initICA20"]] * exp(-a/dM)

      # transmission parameters
      b <- b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
      lambda <- epsilon*b

      # assign to the human
      humans[[j]]$IB <- IB
      humans[[j]]$ID <- ID
      humans[[j]]$ICA <- ICA
      humans[[j]]$ICM <- ICM
      humans[[j]]$epsilon <- epsilon
      humans[[j]]$lambda <- lambda

      # P(clinical disease | infection)
      humans[[j]]$phi <- phi0 * (phi1 + ((1 - phi1)/(1 + ((ICA+ICM)/IC0)^kappaC)))

      # probability of detection
      fD <- 1 - ((1 - fD0)/(1 + (a/aD)^gammaD))
      q <- d1 + ((1 - d1)/(1 + (fD*(ID/ID0)^kappaD)*fD))
      humans[[j]]$prDetectAMic <- q
      humans[[j]]$prDetectAPCR <- q^alphaA
      humans[[j]]$prDetectUPCR <- q^alphaU

      # P(state | immune history, EIR, etc)
      init_s <- state_initial(a = a, EIR = epsilon, theta = theta)
      humans[[j]]$state <- sample(x=names(init_s), size=1, prob = init_s)

      # initial infectivity to mosquito
      humans[[j]]$c <- c_initial(humans[[j]],theta)

      setTxtProgressBar(pb, j)
    }
    close(pb);rm(pb)
    cat(" --- done calculating immunity and state probabilities for population --- \n")

    # set up the mosquitos
    cat("\n --- begin calculating mosquito gonotrophic cycle parameters --- \n")

    delta <- 1.0/(tau1+tau2) # Inverse of gonotrophic cycle without ITNs/IRS
    eggOV <- beta*(exp(muV/delta)-1)/muV # Number of eggs per oviposition per mosquito
    theta["eggOV"] <- eggOV

    w_vec <- sapply(humans,function(x){x$w})
    z_vec <- sapply(humans,function(x){x$z})
    c_vec <- sapply(humans,function(x){x$c})

    # P(successful feed)
    W <- (1 - Q0) + Q0*sum(pi_vec * w_vec)

    # P(repelled w/out feed)
    Z = Q0 * sum(pi_vec * z_vec)

    # feeding rate
    f = 1 / (tau1/(1 - Z) + tau2)

    # survival
    p10 <- exp(-muV*tau1)
    p1 <- p10*W/(1 - Z*p10)
    p2 <- exp(-muV*tau2)
    p12 <- (p1*p2)^f
    mu <- -log(p12)

    # proportion of successful bites on humans & HBR (EIR_tot = a * Iv)
    Q <- 1 - ((1 - Q0)/W)
    a <- f*Q

    # calculate FOI on mosquitos
    C <- sum(pi_vec * c_vec * w_vec)
    lambda_v <- a * C

    # egg laying rate
    betaC <- eggOV*mu/(exp(mu/f) - 1)

    # equilibrum number of infectious vectors
    Iv_eq <- EIR_tot / a
    cat(" --- done calculating mosquito gonotrophic cycle parameters --- \n")

    # solve the mosquitos at equilibrium
    cat("\n --- begin calculating equilibrium values for mosquito population --- \n")

    mosy_eq <- RACD_mosq_equilibrium(theta = theta,dt = 1,IV = Iv_eq,lambdaV = lambda_v,cores=max(2,parallel::detectCores()-2))
    mosy_eq$lambda_v <- lambda_v

    # assign some things to return
    mosy_eq$W <- W
    mosy_eq$Z <- Z
    mosy_eq$mu <- mu
    mosy_eq$Q0 <- Q0
    mosy_eq$C <- C

    cat(" --- done calculating equilibrium values for mosquito population --- \n")

    return(
      list(
        humans=humans,
        mosy=mosy_eq,
        houses=houses,
        theta=theta
      )
    )


  })

}

# RACD_init <- RACD_Setup(N = 200,EIR_mean = 0.05,xy_d = dwell_df,xy_a = aqua_df,theta = RACD_theta)

###############################################################################
# RUN
###############################################################################


# the landscapes


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
p_range <- seq(0.2,0.8,by=0.2)
import_range <- (1/365)*c(1,2,5,10,20)

Rcpp::sourceCpp("/Users/slwu89/Desktop/git/RACD/tiny-racd/intervention-src/main.cpp",rebuild = F)


# begin factorial sweep
l=1
e=1
i=2
im=1
p_ix=1
p_n=2


EIR <- EIR_vals[e]
lscape <- landscape[[l]]

dmat_dwell <- crossdist(X = as.ppp(lscape$dwellings[,1:2],square(r=1)),Y = as.ppp(lscape$dwellings[,1:2],square(r=1)))

# calculate equilibrium
RACD_init <- RACD_Setup(N = N,EIR_mean = EIR,xy_d = lscape$dwellings,xy_a = lscape$habitats,theta = RACD_theta)

# calculate immunity parameters for imported cases
imm_import <- imported_immune(EIR = EIR,theta = RACD_theta)
RACD_theta <- c(RACD_theta,imm_import,import_rate=NaN)

# do the MC iterations (4*6=24)
int_type <- interventions[i]
RACD_theta[["import_rate"]] <- import_range[im]
p_index <- p_range[p_ix]
p_neighbor <- p_range[p_n]

mc_reps <- vector("list",10)
for(iii in 1:30){

  mc_reps[[iii]] <- tiny_racd(humans_param = RACD_init$humans,
                              house_param = RACD_init$houses,
                              mosy_param = RACD_init$mosy,
                              theta = RACD_theta,
                              tmax = (365*2)+100,
                              int_type = int_type,
                              tstart = 101,
                              tend = 100+365,
                              tdelay = 60,
                              dmat = dmat_dwell,
                              radius = 500/2e3,
                              p_index = p_index,
                              p_neighbor = p_neighbor,
                              prog_bar = FALSE)

}
