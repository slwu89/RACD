# -------------------------------------------------------------------------------- #
#
#   P. falciparum simulation (based on Imperial Model)
#   Equilibrium solutions for joint human-mosquito model
#   January 2020
#   Sean Wu (slwu89@berkeley.edu)
#
# -------------------------------------------------------------------------------- #

# load things we need
library(Rcpp)
library(deSolve)
library(here)
library(doParallel)

# Rcpp::sourceCpp(here::here("racd-setup.cpp"))
# source(here::here("mosy-equilibrium.R"))
# source(here::here("racd-parameters.R"))

path2ode <- here::here("setup/humans-ode.cpp")
Rcpp::sourceCpp(path2ode)

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


pfsim_setup <- function(N, EIR_mean, xy_d, xy_a, theta, cores = 4, parseed = NaN){

	if(is.nan(parseed)){
		stop("please enter a seed for 'parseed'")
	}

	# -------------------------------------------------------------------------------- #
	#		Extract parameters
	# -------------------------------------------------------------------------------- #

	beta  <- theta[["beta"]]
	muEL  <- theta[["muEL"]]
	muLL  <- theta[["muLL"]]
	muPL  <- theta[["muPL"]]
	durEL  <- theta[["durEL"]]
	durLL  <- theta[["durLL"]]
	durPL  <- theta[["durPL"]]
	durEV  <- theta[["durEV"]]
	gamma  <- theta[["gamma"]]
	tau1  <- theta[["tau1"]]
	tau2  <- theta[["tau2"]]

	muV  <- theta[["muV"]]
	Q0  <- theta[["Q0"]]
	phiB  <- theta[["phiB"]]
	phiI  <- theta[["phiI"]]
	rITN  <- theta[["rITN"]]
	sITN  <- theta[["sITN"]]
	rIRS  <- theta[["rIRS"]]
	sIRS  <- theta[["sIRS"]]
	delta <- theta[["delta"]]
	eggOV <- theta[["eggOV"]]

	fT  <- theta[["fT"]]
	dE  <- theta[["dE"]]
	dT  <- theta[["dT"]]
	dD  <- theta[["dD"]]
	dA  <- theta[["dA"]]
	dU  <- theta[["dU"]]
	dP  <- theta[["dP"]]

	cD  <- theta[["cD"]]
	cT  <- theta[["cT"]]
	cU  <- theta[["cU"]]
	gammaI  <- theta[["gammaI"]]

	rho  <- theta[["rho"]]
	a0  <- theta[["a0"]]
	sigma2  <- theta[["sigma2"]]

	d1  <- theta[["d1"]]
	dID  <- theta[["dID"]]
	ID0  <- theta[["ID0"]]
	kappaD  <- theta[["kappaD"]]
	uD  <- theta[["uD"]]
	aD  <- theta[["aD"]]
	fD0  <- theta[["fD0"]]
	gammaD  <- theta[["gammaD"]]
	alphaA  <- theta[["alphaA"]]
	alphaU  <- theta[["alphaU"]]

	b0  <- theta[["b0"]]
	b1  <- theta[["b1"]]
	dB  <- theta[["dB"]]
	IB0  <- theta[["IB0"]]
	kappaB  <- theta[["kappaB"]]
	uB  <- theta[["uB"]]

	phi0  <- theta[["phi0"]]
	phi1  <- theta[["phi1"]]
	dC  <- theta[["dC"]]
	IC0  <- theta[["IC0"]]
	kappaC  <- theta[["kappaC"]]
	uC  <- theta[["uC"]]
	PM  <- theta[["PM"]]
	dM  <- theta[["dM"]]

	rW  <- theta[["rW"]]
	rP  <- theta[["rP"]]

	meanAge  <- theta[["meanAge"]]
	mu <- theta[["mu"]]

	# -------------------------------------------------------------------------------- #
	#		calculate EIR and EIR on houses; then sample biting heterogeneity and age to get epsilon (personal EIR)
	# -------------------------------------------------------------------------------- #

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
		humans[[j]]$pi_h <- pi_vec[j]
	}
	cat(" --- done calculating landscape-dependent parameters --- \n")


	# -------------------------------------------------------------------------------- #
	# 	begin calculating immunity and state probabilities for population
	# -------------------------------------------------------------------------------- #

	cat("\n --- begin calculating immunity and state probabilities for population --- \n")
	# mean value of ICA for 20-year olds
	theta[["initICA20"]] <- immune_initial(a = 20,EIR = EIR_mean,theta = theta)[["ICA"]]

	# set up parallel core; compile necessary C++ on them
	cl <- makeCluster(cores)
	registerDoParallel(cl)
	clusterSetRNGStream(cl = cl, iseed = parseed)
	clusterExport(cl = cl,varlist = c("path2ode"))
	clusterEvalQ(cl,{
  	Rcpp::sourceCpp(path2ode)
	})

	h_imm <- 	foreach(h = iter(humans), .export = c("immune_initial","state_initial","c_initial","EIR_tot","theta")) %dopar% {

		# output
		h_out <- list()

		# immune state ODEs
		a <- h$age
		epsilon <- EIR_tot * h$pi_h
		init_i <- immune_initial(a = a, EIR = epsilon, theta = theta)

		# my initial immune values
		IB <- init_i[["IB"]]
		ID <- init_i[["ID"]]
		ICA <- init_i[["ICA"]]
		ICM <- theta[["initICA20"]] * exp(-a/theta[["dM"]])

		# transmission parameters
		b <- theta[["b0"]]*(theta[["b1"]] + ((1-theta[["b1"]])/(1 + (IB/theta[["IB0"]])^theta[["kappaB"]])))
		lambda <- epsilon*b

		# assign to the human
		h_out$IB <- IB
		h_out$ID <- ID
		h_out$ICA <- ICA
		h_out$ICM <- ICM
		h_out$epsilon <- epsilon
		h_out$lambda <- lambda

		# P(clinical disease | infection)
		h_out$phi <- theta[["phi0"]] * (theta[["phi1"]] + ((1 - theta[["phi1"]])/(1 + ((ICA+ICM)/theta[["IC0"]])^theta[["kappaC"]])))

		# probability of detection
		fD <- 1 - ((1 - theta[["fD0"]])/(1 + (a/theta[["aD"]])^theta[["gammaD"]]))
		q <- theta[["d1"]] + ((1 - theta[["d1"]])/(1 + (fD*(ID/theta[["ID0"]])^theta[["kappaD"]])*fD))
		h_out$prDetectAMic <- q
		h_out$prDetectAPCR <- q^theta[["alphaA"]]
		h_out$prDetectUPCR <- q^theta[["alphaU"]]

		# P(state | immune history, EIR, etc)
		init_s <- state_initial(a = a, EIR = epsilon, theta = theta)
		h_out$state <- sample(x=names(init_s), size=1, prob = init_s)

		# initial infectivity to mosquito
		h_out$c <- c_initial(h_out,theta)

		# what is returned
		h_out
	}

	# merge the results
	humans <- mapply(function(h0,h1){
    c(h0,h1)
	},h0=humans,h1=h_imm,SIMPLIFY = FALSE)

	stopCluster(cl);rm(cl);gc()
	cat(" --- done calculating immunity and state probabilities for population --- \n")


	# -------------------------------------------------------------------------------- #
	# 	begin calculating mosquito gonotrophic cycle parameters
	# -------------------------------------------------------------------------------- #

	cat("\n --- begin calculating mosquito gonotrophic cycle parameters --- \n")

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

	mosy_eq <- RACD_mosq_equilibrium(theta = theta,dt = 1,IV = Iv_eq,lambdaV = lambda_v)
	mosy_eq$lambda_v <- lambda_v

	# assign some things to return
	mosy_eq$W <- W
	mosy_eq$Z <- Z
	mosy_eq$mu <- mu
	mosy_eq$Q0 <- Q0
	mosy_eq$C <- C

	cat(" --- done calculating equilibrium values for mosquito population --- \n")

	# RETURN TO R
	return(
		list(
			humans=humans,
			mosy=mosy_eq,
			houses=houses,
			theta=theta
		)
	)

}
