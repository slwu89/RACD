###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   April 2019
#   Setup the RACD model (at equilibrium)
#
###############################################################################

# load things we need
library(Rcpp)
library(deSolve)
library(here)

Rcpp::sourceCpp(here::here("racd-setup.cpp"))
source(here::here("mosy-equilibrium.R"))
source(here::here("racd-parameters.R"))

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
c_initial <- function(human){
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
}


RACD_Setup <- function(N, EIR_mean, xy_d, xy_a, theta){

	# extract variables that we need and assign here
	invisible(mapply(FUN = function(val,name){
    assign(x = name,value = val,pos = 1)
	},val = unname(theta),name = names(theta)))

	theta[["mu"]] <- mu <- 1/(meanAge*365) # Daily death rate as a function of mean age in years

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
		humans[[j]]$c <- c_initial(humans[[j]])

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
  mu <- -f*log(p1*p1)

	# proportion of successful bites on humans & HBR (EIR_tot = a * Iv)
	Q <- 1 - ((1 - Q0)/W)
	a <- f*Q

	# calculate FOI on mosquitos
	lambda_v <- a * sum(pi_vec * c_vec * w_vec)

	# egg laying rate
	betaC <- eggOV*mu/(exp(mu/f) - 1)

	# equilibrum number of infectious vectors
	Iv_eq <- EIR_tot / a
	cat(" --- done calculating mosquito gonotrophic cycle parameters --- \n")

	# solve the mosquitos at equilibrium
	cat("\n --- begin calculating equilibrium values for mosquito population --- \n")

	mosy_eq <- RACD_mosq_equilibrium(theta = theta,dt = 1,IV = Iv_eq,lambdaV = lambda_v,cores=max(2,parallel::detectCores()-2))

	cat(" --- done calculating equilibrium values for mosquito population --- \n")

	return(
		list(
			humans=humans,
			mosy=mosy_eq
		)
	)
}

RACD_init <- RACD_Setup(N = 200,EIR_mean = 0.05,xy_d = dwell_df,xy_a = aqua_df,theta = RACD_theta)
