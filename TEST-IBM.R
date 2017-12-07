#########################################################################
## Malaria village-scale individual-based model                        ##
## John Marshall (john.marshall@berkeley.edu)                          ##
## 11/Jun/2015                                                         ##
#########################################################################

# Clear all stored parameters:
rm(list=ls())

# Set working directory:
setwd("C://Users/John/My Documents/Berkeley/MEI/MicroEpiModeling/")

# Load required libraries:
library(deSolve)
library(ggplot2)
library(RColorBrewer)

theta <- c(
  ## Variable model parameters:
  epsilon0 = 5/365, # Mean EIR for adults (per day)
  fT = 0.4, # Proportion of clinical disease cases successfully treated

  ## Model parameters taken from Griffin et al. (2014):
  ## Human infection durations:
  dE = 12, # Duration of latent period (days)
  dT = 5, # Duration of treated clinical disease (days)
  dD = 5, # Duration of untreated clinical disease (days)
  dA = 200, # Duration of patent infection (days)
  dU = 110, # Duration of sub-patent infection (days) (fitted)
  dP = 25, # Duration of prophylactic protection following treatment (days)

  ## Infectiousness of humans to mosquitoes:
  cD = 0.068, # Infectiousness with untreated disease & no immunity (fitted)
  cT = 0.32 * 0.068, # Infectiousness after treatment
  cU = 0.0062, # Infectiousness with sub-patent infection (fitted)
  gammaI = 1.82, # Relates infectiousness to probability of detection (fitted)

  ## Age and heterogeneity parameters:
  rho = 0.85, # Age-dependent biting parameter
  a0 = 8, # Age-dependent biting parameter (years)
  sigma2 = 1.67, # Variance of log of heterogeneity in biting rates

  ## Effect of immunity on reducing probability of detection:
  d1 = 0.161, # Probability of detection with maximum immunity (fitted)
  dID = 10*365, # Inverse of decay rate (days)
  ID0 = 1.58, # Immunity scale parameter (fitted)
  kappaD = 0.477, # Immunity shape parameter (fitted)
  uD = 9.45, # Duration in which immunity is not boosted (days) (fitted)
  aD = 21.9, # Scale parameter relating age to immunity (years) (fitted)
  fD0 = 0.0071, # Parameter relating age to immunity (fitted)
  gammaD = 4.81, # Shape parameter relating age to immunity (fitted)
  alphaA = 0.757, # PCR prevalence parameter (fitted)
  alphaU = 0.186, # PCR prevalence parameter (fitted)

  ## Immunity reducing probability of infection:
  b0 = 0.590, # Probabiliy with no immunity (fitted)
  b1 = 0.5, # Maximum relative reduction
  dB = 10*365, # Inverse of decay rate (days)
  IB0 = 43.9, # Scale parameter (fitted)
  kappaB = 2.16, # Shape parameter (fitted)
  uB = 7.20, # Duration in which immunity is not boosted (days) (fitted)

  ## Immunity reducing probability of clinical disease:
  phi0 = 0.792, # Probability with no immunity
  phi1 = 0.00074, # Maximum relative reduction
  dC = 30*365, # Inverse decay rate (days)
  IC0 = 18.0, # Scale parameter
  kappaC = 2.37, # Shape parameter
  uC = 6.06, # Duration in which immunity is not boosted (days)
  PM = 0.774, # New-born immunity relative to mother's immunity
  dM = 67.7, # Inverse decay rate of maternal immunity (days)

  ## Case detection (recorded incidence relative to daily active case
  ## detection):
  rW = 0.723, # Weekly active case detection
  rP = 0.342, # Weekly passive case detection

  ## Demographic parameters:
  meanAge = 17.4, # Mean age in Tanzania (males and females, years)
  N = 200, # Village population size

  ## Geographic parameters:
  meanNumPeoplePerHouse = 6.5, # Mean number of people per house (from Misungu data set)
  numHousesPerBreedingSite = 5) # Number of houses per breeding site

##############################################################################
## Functions for calculating levels of immunity when the simulation starts: ##
##############################################################################

## Differential equations for:
## 1. Pre-erythrocytic immunity (IB, reduces the probability of infection
##    following an infectious challenge).
## 2. Detection immunity (ID, a.k.a. blood-stage immunity, reduces the
##    probability of detection and reduces infectiousness to mosquitoes).
## 3. Acquired clinical immunity (ICA, reduces the probability of clinical
##    disease, acquired from previous exposure).

I_ode <- function(time, state, theta) {
  ## States:
  IB <- state["IB"]
  ID <- state["ID"]
  ICA <- state["ICA"]

  ## Parameters:
  a0 <- theta[["a0"]]
  rho <- theta[["rho"]]
  durB <- theta[["dB"]] / 365 # Inverse decay rate (years)
  uB <- theta[["uB"]] / 365 # Duration in which immunity is not boosted (years)
  durD <- theta[["dID"]] / 365 # Inverse decay rate (years)
  uD <- theta[["uD"]] / 365 # Duration in which immunity is not boosted (years)
  durC <- theta[["dC"]] / 365 # Inverse decay rate (years)
  uC <- theta[["uC"]] / 365 # Duration in which immunity is not boosted (years)
  zita <- theta[["zita"]]
  psi <- theta[["psi"]]
  b0 <- theta[["b0"]]
  b1 <- theta[["b1"]]
  IB0 <- theta[["IB0"]]
  kappaB <- theta[["kappaB"]]
  epsilon0 <- theta[["epsilon0"]] * 365 # Entomological innoculation rate (annual)
  epsilon <- epsilon0*zita*(1 - rho*exp(-time/a0))*psi
  # b <- b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
  lambda <- epsilon*b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))

  # cat("power in lambda: ",(IB/IB0)^kappaB,"\n")
  # cat("exp in epsilon: ",exp(-time/a0),"\n")

  # cat("time: ",time," epsilon: ",epsilon," b: ",b," lambda: ",lambda,"\n")

  ## ODEs:
  dIB <- epsilon/(epsilon*uB + 1) - IB/durB
  dID <- lambda/(lambda*uD + 1) - ID/durD
  dICA <- lambda/(lambda*uC + 1) - ICA/durC

  return(list(c(dIB, dID, dICA)))
}

## Wrapper function to calculate pre-erythrocytic immunity at age a for
## given EIR heterogeneities (biting rate & geographic location):
initialI <- function(a, zita, psi, theta) {
  initState <- c(IB=0, ID=0, ICA=0)
  thetaI <- c(a0=theta[["a0"]], rho=theta[["rho"]], dB=theta[["dB"]],
              uB=theta[["uB"]], epsilon0=theta[["epsilon0"]],
              dID=theta[["dID"]], uD=theta[["uD"]], dC=theta[["dC"]],
              uC=theta[["uC"]], b0=theta[["b0"]], b1=theta[["b1"]],
              IB0=theta[["IB0"]], kappaB=theta[["kappaB"]], psi=psi, zita=zita)
  IvsAge <- data.frame(ode(y=initState, times=seq(0, a, by=a/1000), func=I_ode,
                           parms=thetaI, method="ode45"))
  c(IB=IvsAge$IB[length(IvsAge$IB)], ID=IvsAge$ID[length(IvsAge$ID)],
    ICA=IvsAge$ICA[length(IvsAge$ICA)])
}

a=20
zita=1
psi=1
time = seq(from=0,to=a,by=a/1e3)
thetaODE = c(a0=theta[["a0"]], rho=theta[["rho"]], dB=theta[["dB"]],
             uB=theta[["uB"]], epsilon0=theta[["epsilon0"]],
             dID=theta[["dID"]], uD=theta[["uD"]], dC=theta[["dC"]],
             uC=theta[["uC"]], b0=theta[["b0"]], b1=theta[["b1"]],
             IB0=theta[["IB0"]], kappaB=theta[["kappaB"]], psi=psi, zita=zita)
testOut = ode(y=c(IB=0, ID=0, ICA=0), times=time, func=I_ode,
              parms=thetaODE, method="ode45")


thetaC = c(a0=8,rho=0.85,zeta=1,psi=1,
           durB=(3650/365),uB=(7.2/365),b0=0.59,b1=0.5,IB0=43.9,kappaB=2.16,
           durD=(3650/365),uD=(9.45/365),
           durC=(10950/365),uC=(6.06/365),epsilon0=(0.01369863*365))
testC = RACD::immune_ode(time = time,theta = thetaC,state = c(IB=0, ID=0, ICA=0))


# microbenchmark::microbenchmark(
#   {ode(y=c(IB=0, ID=0, ICA=0), times=time, func=I_ode,parms=thetaODE)},
#   {RACD::immune_ode(time = time,theta = thetaC,state = c(IB=0, ID=0, ICA=0))}
# )
##############################################################################
## Functions for calculating distribution of states when simulation starts: ##
##############################################################################

## Differential equations for calculating the probability that an individual
## is in each state given their age and EIR heterogeneity attributes:

malaria_ode <- function(time, state, theta) {
  ## States:
  prS <- state["prS"]
  prT <- state["prT"]
  prD <- state["prD"]
  prA <- state["prA"]
  prU <- state["prU"]
  prP <- state["prP"]
  IB <- state["IB"]
  ICA <- state["ICA"]

  # Parameters:
  ## Variable model parameters:
  epsilon0 <- theta[["epsilon0"]] * 365 # Mean EIR for adults (per year)
  fT <- theta[["fT"]] # Proportion of clinical disease cases successfully treated

  ## Model parameters taken from Griffin et al. (2014):
  ## Human infection durations:
  dE <- theta[["dE"]] / 365 # Duration of latent period (years)
  dT <- theta[["dT"]] / 365 # Duration of treated clinical disease (years)
  dD <- theta[["dD"]] / 365 # Duration of untreated clinical disease (years)
  dA <- theta[["dA"]] / 365 # Duration of patent infection (years)
  dU <- theta[["dU"]] / 365 # Duration of sub-patent infection (years) (fitted)
  dP <- theta[["dP"]] / 365 # Duration of prophylactic protection following treatment (years)

  ## Age parameters:
  rho <- theta[["rho"]] # Age-dependent biting parameter
  a0 <- theta[["a0"]] # Age-dependent biting parameter (years)

  ## Immunity reducing probability of infection:
  b0 <- theta[["b0"]] # Probabiliy with no immunity (fitted)
  b1 <- theta[["b1"]] # Maximum relative reduction
  dB <- theta[["dB"]] / 365 # Inverse of decay rate (years)
  IB0 <- theta[["IB0"]] # Scale parameter (fitted)
  kappaB <- theta[["kappaB"]] # Shape parameter (fitted)
  uB <- theta[["uB"]] / 365 # Duration in which immunity is not boosted (years) (fitted)

  ## Immunity reducing probability of clinical disease:
  phi0 <- theta[["phi0"]] # Probability with no immunity
  phi1 <- theta[["phi1"]] # Maximum relative reduction
  dC <- theta[["dC"]] / 365 # Inverse decay rate (years)
  IC0 <- theta[["IC0"]] # Scale parameter
  kappaC <- theta[["kappaC"]] # Shape parameter
  uC <- theta[["uC"]] / 365 # Duration in which immunity is not boosted (years)
  PM <- theta[["PM"]] # New-born immunity relative to mother's immunity
  dM <- theta[["dM"]] / 365 # Inverse decay rate of maternal immunity (years)
  initICA20 <- theta["initICA20"]
  ICM <- initICA20*exp(-time/dM)

  ## Individual heterogeneity parameters:
  zita <- theta[["zita"]]
  psi <- theta[["psi"]]
  epsilon <- epsilon0*zita*(1 - rho*exp(-time/a0))*psi
  b <- b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
  lambda <- epsilon*b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
  phi <- phi0*(phi1 + ((1 - phi1)/(1 + ((ICA + ICM)/IC0)^kappaC)))
  cat("epsilon0 ",epsilon0," zeta ",zita," rho ",rho," time ",time," a0 ",a0," psi ",psi,"\n")
  cat("time: ",time," ICM: ",ICM," epsilon: ",epsilon," b: ",b," lambda: ",lambda," phi: ",phi,"\n")

  ## ODEs:
  dprS <- - lambda*prS + prP/dP + prU/dU
  dprT <- phi*fT*lambda*(prS + prA + prU) - prT/dT
  dprD <- phi*(1 - fT)*lambda*(prS + prA + prU) - prD/dD
  dprA <- (1 - phi)*lambda*(prS + prA + prU) + prD/dD - lambda*prA - prA/dA
  dprU <- prA/dA - prU/dU - lambda*prU
  dprP <- prT/dT - prP/dP
  dIB <- epsilon/(epsilon*uB + 1) - IB/dB
  dICA <- lambda/(lambda*uC + 1) - ICA/dC

  return(list(c(dprS, dprT, dprD, dprA, dprU, dprP, dIB, dICA)))
}

## Wrapper function for proportion belonging to each state for each age
## at the beginning of the simulation:
initialStates <- function(a, zita, psi, theta) {
  ## All individuals begin as susceptible when born:
  initState <- c(prS=1, prT=0, prD=0, prA=0, prU=0, prP=0,
                 IB=0, ICA=0)

  ## The vector of parameters, with additional heterogeneity parameters:
  thetaI <- theta
  thetaI["zita"] <- zita
  thetaI["psi"] <- psi
  initI20 <- initialI(a=20, zita=1, psi=1, theta=theta)
  thetaI["initICA20"] <- initI20[["ICA"]]

  ## Running the system of ODEs:
  prStates <- data.frame(ode(y=initState, times=seq(0, a, by=a/1000),
                             func=malaria_ode, parms=thetaI, method="ode45"))

  c(prS=prStates$prS[length(prStates$prS)],
    prT=prStates$prT[length(prStates$prT)],
    prD=prStates$prD[length(prStates$prD)],
    prA=prStates$prA[length(prStates$prA)],
    prU=prStates$prU[length(prStates$prU)],
    prP=prStates$prP[length(prStates$prP)])
}

a=1
zita=1
psi=1
thetaI <- theta
thetaI["zita"] <- zita
thetaI["psi"] <- psi
initI20 <- initialI(a=20, zita=1, psi=1, theta=theta)
thetaI["initICA20"] <- initI20[["ICA"]]

times = seq(0, a, by=a/4)
rout = ode(y=c(prS=1, prT=0, prD=0, prA=0, prU=0, prP=0,IB=0, ICA=0), times=times,func=malaria_ode, parms=thetaI)

thetaC = c(thetaI["epsilon0"],thetaI["fT"],thetaI["dE"],thetaI["dT"],thetaI["dD"],thetaI["dA"],thetaI["dU"],thetaI["dP"],thetaI["rho"],thetaI["a0"],thetaI["b0"],thetaI["b1"],thetaI["dB"],thetaI["IB0"],thetaI["kappaB"],thetaI["uB"],thetaI["phi0"],thetaI["phi1"],thetaI["dC"],thetaI["IC0"],thetaI["kappaC"],thetaI["uC"],thetaI["PM"],thetaI["dM"])
thetaC = c(thetaC,initICA20=initI20[["ICA"]],zeta=1,psi=1)
cout = infection_ode(time = times,theta = thetaC,state = c(prS=1, prT=0, prD=0, prA=0, prU=0, prP=0,IB=0, ICA=0))

numIter <- 10*365 # Simulation over 10 years


rInitialStates = initialStates(a = 20,zita = 1,psi = 1,theta = theta)
cInitialStates = infection_initial(a = 20,zeta = 1,psi = 1,method="ode45")

#########################################################################
## Coding the Imperial College model as an IBM:                        ##
## Part 1. Setting up the IBM framework                                ##
#########################################################################

malaria_ibm <- function(theta, numIter) {
	# Parameters:
	## Variable model parameters:
	epsilon0 <- theta[["epsilon0"]] # Mean EIR for adults (per day)
	fT <- theta[["fT"]] # Proportion of clinical disease cases successfully treated

	## Model parameters taken from Griffin et al. (2014):
	## Human infection durations:
	dE <- theta[["dE"]] # Duration of latent period (days)
	dT <- theta[["dT"]] # Duration of treated clinical disease (days)
	dD <- theta[["dD"]] # Duration of untreated clinical disease (days)
	dA <- theta[["dA"]] # Duration of patent infection (days)
	dU <- theta[["dU"]] # Duration of sub-patent infection (days) (fitted)
	dP <- theta[["dP"]] # Duration of prophylactic protection following treatment (days)

	## Infectiousness of humans to mosquitoes:
	cD <- theta[["cD"]] # Infectiousness with untreated disease & no immunity (fitted)
	cT <- theta[["cT"]] # Infectiousness after treatment
	cU <- theta[["cU"]] # Infectiousness with sub-patent infection (fitted)
	gammaI <- theta[["gammaI"]] # Relates infectiousness to probability of detection (fitted)

	## Age and heterogeneity parameters:
	rho <- theta[["rho"]] # Age-dependent biting parameter
	a0 <- theta[["a0"]] # Age-dependent biting parameter (years)
	sigma2 <- theta[["sigma2"]] # Variance of log of heterogeneity in biting rates

	## Effect of immunity on reducing probability of detection:
	d1 <- theta[["d1"]] # Probability of detection with maximum immunity (fitted)
	dID <- theta[["dID"]] # Inverse of decay rate (days)
	ID0 <- theta[["ID0"]] # Immunity scale parameter (fitted)
	kappaD <- theta[["kappaD"]] # Immunity shape parameter (fitted)
	uD <- theta[["uD"]] # Duration in which immunity is not boosted (fitted)
	aD <- theta[["aD"]] # Scale parameter relating age to immunity (years) (fitted)
	fD0 <- theta[["fD0"]] # Parameter relating age to immunity (fitted)
	gammaD <- theta[["gammaD"]] # Shape parameter relating age to immunity (fitted)
	alphaA <- theta[["alphaA"]] # PCR prevalence parameter (fitted)
	alphaU <- theta[["alphaU"]] # PCR prevalence parameter (fitted)

	## Immunity reducing probability of infection:
	b0 <- theta[["b0"]] # Probabiliy with no immunity (fitted)
	b1 <- theta[["b1"]] # Maximum relative reduction
	dB <- theta[["dB"]] # Inverse of decay rate (days)
	IB0 <- theta[["IB0"]] # Scale parameter (fitted)
	kappaB <- theta[["kappaB"]] # Shape parameter (fitted)
	uB <- theta[["uB"]] # Duration in which immunity is not boosted (days) (fitted)

	## Immunity reducing probability of clinical disease:
	phi0 <- theta[["phi0"]] # Probability with no immunity
	phi1 <- theta[["phi1"]] # Maximum relative reduction
	dC <- theta[["dC"]] # Inverse decay rate (days)
	IC0 <- theta[["IC0"]] # Scale parameter
	kappaC <- theta[["kappaC"]] # Shape parameter
	uC <- theta[["uC"]] # Duration in which immunity is not boosted
	PM <- theta[["PM"]] # New-born immunity relative to mother's immunity
	dM <- theta[["dM"]] # Inverse decay rate of maternal immunity

	## Case detection (recorded incidence relative to daily active case
	## detection):
	rW <- theta[["rW"]] # Weekly active case detection
	rP <- theta[["rP"]] # Weekly passive case detection

	## Demographic parameters:
	N <- theta[["N"]] # Village population size
	meanAge <- theta[["meanAge"]] # Mean age in Tanzania (males and females, years)
	mu <- 1/(meanAge*365) # Daily death rate as a function of mean age in years

	## Geographic parameters:
	meanNumPeoplePerHouse <- theta[["meanNumPeoplePerHouse"]] # Mean number of people per house (from Misungu data set)
	numHouses <- round(N/meanNumPeoplePerHouse)
	numHousesPerBreedingSite <- theta[["numHousesPerBreedingSite"]] # Number of houses per breeding site
	numBreedingSites <- round(numHouses/numHousesPerBreedingSite) # Number of breeding sites

	# Use a multi-level list where elements represent individuals and
	# sub-levels contain attribute information about the individuals:
	indiv <- vector(mode = "list", N)

	# Randomly assign age attributes:
	# (Age is sampled from an exponential distribution with mean equal to
	# the mean age in the country being considered)
	for (j in 1:N) {
		indiv[[j]]$age <- rexp(n = 1, rate = 1/meanAge)
	}
	# Checking initial age attributes:
	# ageVector <- sapply(indiv, function(x) x$age)
	# hist(ageVector)

	# At the beginning of the simulation, all individuals are alive:
	for (j in 1:N) {
		indiv[[j]]$alive <- 1
	}

	# Randomly assign house coordinates:
	longHouse <- rep(0, numHouses) # Vector of house longitudes
	latHouse <- rep(0, numHouses) # Vector of house latitudes
	longHouse[1] <- runif(1) # First house longitude coordinate
	latHouse[1] <- runif(1) # First house latitude coordinate
	for (j in 2:numHouses) {
		distOtherJHouses <- rep(0, j-1)
		# Ensure that the remainder of the houses are separated
		# by a minimum distance of 0.05
		while(min(distOtherJHouses) < 0.05) {
			longHouse[j] <- runif(1)
			latHouse[j] <- runif(1)
			distOtherJHouses <- rep(0, j-1)
			for (k in 1:(j-1)) {
				distOtherJHouses[k] <- sqrt((longHouse[j]-longHouse[k])^2
								  + (latHouse[j]-latHouse[k])^2)
			}
		}
	}
	# Checking house distribution:
	# plot(longHouse, latHouse, pch=16, col='blue')

	# Randomly assign individuals to houses:
	householdSize <- rep(0, numHouses) # Vector of household sizes
	for (j in 1:N) {
		# Assign individuals to one of the samllest houses:
		smallestHouse <- which(householdSize==min(householdSize))
		if (length(smallestHouse)==1) { indiv[[j]]$house <- smallestHouse }
		else { indiv[[j]]$house <- sample(smallestHouse, 1) }
		householdSize[indiv[[j]]$house] <- householdSize[indiv[[j]]$house] + 1
	}
	# Checking household size distribution:
	# householdSize

	# Randomly assign breeding site coordinates:
	longBreedingSite <- rep(0, numBreedingSites) # Vector of breeding site longitudes
	latBreedingSite <- rep(0, numBreedingSites) # Vector of breeding site latitudes
	for (j in 1:numBreedingSites) {
		distHousesBreedingSites <- rep(0, numHouses+j-1)
		# Ensure that houses and remainder of breeding sites are separated
		# by a minimum distance of 0.05
		while(min(distHousesBreedingSites) < 0.05) {
			longBreedingSite[j] <- runif(1)
			latBreedingSite[j] <- runif(1)
			for (k in 1:numHouses) {
				distHousesBreedingSites[k] <- sqrt((longBreedingSite[j]-longHouse[k])^2
							               + (latBreedingSite[j]-latHouse[k])^2)
			}
			if (j>1) {
				for (k in 1:(j-1)) {
					distHousesBreedingSites[numHouses+k] <- sqrt((longBreedingSite[j]-longBreedingSite[k])^2
								               		 + (latBreedingSite[j]-latBreedingSite[k])^2)
				}
			}
		}
	}
	# Checking breeding site distribution:
	# plot(longBreedingSite, latBreedingSite, type='p', pch=16, col='red')
	# points(longHouse, latHouse, pch=16, col='blue')

	# Define geographical risk surface due to breeding sites:
	# We model the risk due to breeding sites at each household as the sum
	# of multivariate normal distributions centered at each breeding site
	# with the standard deviation of the normal distribution being equal to
	# the distance between the breeding site and the nearest household.
	# 1. First, we calculate the standard deviation of the risk surface for
	#    each breeding site.
	sigmaBreedingSite <- rep(0, numBreedingSites)
	for (j in 1:numBreedingSites) {
		distHousesBreedingSites <- rep(0, numHouses)
		for (k in 1:numHouses) {
			distHousesBreedingSites[k] <- sqrt((longBreedingSite[j]-longHouse[k])^2
						               + (latBreedingSite[j]-latHouse[k])^2)
		}
		sigmaBreedingSite[j] <- min(distHousesBreedingSites)
	}
	# 2. Next, calculate the contribution of each breeding site to the
	#    relative risk value for each house.
	psiHouse <- rep(0, numHouses)
	for (j in 1:numHouses) {
		for (k in 1:numBreedingSites) {
			psiHouse[j] <- ( psiHouse[j] + (1/(sigmaBreedingSite[k]*sqrt(2*pi)))
					     * exp(-((longHouse[j]-longBreedingSite[k])^2 +
                                   (latHouse[j]-latBreedingSite[k])^2)
                                   / (2*sigmaBreedingSite[k]^2)))
		}
	}
	# 3. Finally, we normalize the psiHouse values so they have a mean of 1.
	psiHouse <- psiHouse*numHouses/sum(psiHouse)
	# Check distribution of psiHouse:
	# mean(psiHouse)
	# psiHouse

	# Randomly assign biting heterogeneity attributes:
	# (Biting heterogeneity is sampled from a normal distribution with mean
	# 1 and standard deviation sigma)
	# (Referred to as zita in Griffin et al. (2014))
	for (j in 1:N) {
		indiv[[j]]$bitingHet <- rlnorm(n = 1, meanlog = -sigma2/2, sdlog = sqrt(sigma2))
	}
	# Checking biting heterogeneity attributes:
	# bitingHetVector <- sapply(indiv, function(x) x$bitingHet)
	# hist(bitingHetVector)

	# Lambda (the force of infection) is calculated for each individual. It
	# varies according to age and biting heterogeneity group.
	# Immunity values are also calculated here since these depend on the
	# the values of epsilon (the entomological inoculation rate) and lambda:
	# 1. Pre-erythrocytic immunity (IB, reduces the probability of infection
	#    following an infectious challenge)
	# 2. Acquired clinical immunity (ICA, reduces the probability of clinical
	#    disease, acquired from previous exposure)
	# 3. Maternal clinical immunity (ICM, reduces the probability of clinical
	#    disease, acquired maternally)
	# 4. Detection immunity (ID, a.k.a. blood-stage immunity, reduces the
	#    probability of detection and reduces infectiousness to mosquitoes)
	for (j in 1:N) {
		zita <- indiv[[j]]$bitingHet
		psi <- psiHouse[indiv[[j]]$house]
		a <- indiv[[j]]$age
		# Calculate initial immunity levels from their differential equations:
		initI <- initialI(a=a, zita=zita, psi=psi, theta=theta)
		IB <- initI[["IB"]]; ID <- initI[["ID"]]; ICA <- initI[["ICA"]]
		initI20 <- initialI(a=20, zita=1, psi=1, theta=theta)
		ICM <- initI[["ICA"]] * exp(-a/dM)
		epsilon <- epsilon0 * zita * (1 - rho * exp(-a/a0)) * psi
		b <- b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
		lambda <- epsilon*b
		indiv[[j]]$IB <- IB
		indiv[[j]]$ID <- ID
		indiv[[j]]$ICA <- ICA
		indiv[[j]]$ICM <- ICM
		indiv[[j]]$epsilon <- epsilon
		indiv[[j]]$lambda <- lambda
	}
	# Checking EIR and FOI attributes:
	# EIRVector <- sapply(indiv, function(x) x$epsilon)
	# hist(EIRVector)
	# FOIVector <- sapply(indiv, function(x) x$lambda)
	# hist(FOIVector)

	# Phi (the probability of acquiring clinical disease upon infection) is
	# also calculated for each individual. It varies according to immune
	# status:
	for (j in 1:N) {
		ICA <- indiv[[j]]$ICA
		ICM <- indiv[[j]]$ICM
		indiv[[j]]$phi <- phi0 * (phi1 + ((1 - phi1)/(1 + ((ICA+ICM)/IC0)^kappaC)))
	}
	# Checking phi (probability of acquiring clinical disease) attributes:
	# phiVector <- sapply(indiv, function(x) x$phi)
	# hist(phiVector)

	# q (the probability that an asymptomatic infection is detected by
	# microscopy) is also calculated for each individual, as well as the
	# probability of detection by PCR for asymptomatic infections in states
	# A (patent) and U (subpatent). This also varies according to immune
	# status:
	for (j in 1:N) {
		a <- indiv[[j]]$age
		ID <- indiv[[j]]$ID
		fD <- 1 - ((1 - fD0)/(1 + (a/aD)^gammaD))
		q <- d1 + ((1 - d1)/(1 + (fD*(ID/ID0)^kappaD)*fD))
		indiv[[j]]$prDetectAMic <- q
		indiv[[j]]$prDetectAPCR <- q^alphaA
		indiv[[j]]$prDetectUPCR <- q^alphaU
	}
	# Checking phi (probability of acquiring clinical disease) attributes:
	# prDetectAMicVector <- sapply(indiv, function(x) x$prDetectAMic)
	# hist(prDetectAMicVector)

	# For each invididual, use the ODE transmission model to determine the
	# probability that they are in each state given their age and EIR
	# heterogeneity attributes:
	cat("Initializing population...\n")
	pb <- txtProgressBar(min = 1, max = N, initial = 1, style=3)
	for (j in 1:N) {
		a <- indiv[[j]]$age
		zita <- indiv[[j]]$bitingHet
		psi <- psiHouse[indiv[[j]]$house]
		prStateVector <- initialStates(a=a, zita=zita, psi=psi, theta=theta)
		prS <- prStateVector[["prS"]]
		prT <- prStateVector[["prT"]]
		prD <- prStateVector[["prD"]]
		prA <- prStateVector[["prA"]]
		prU <- prStateVector[["prU"]]
		prP <- prStateVector[["prP"]]
		randNum <- runif(1)
		if (randNum <= prS) {
			indiv[[j]]$state <- "S"
		} else if ((randNum > prS) & (randNum <= (prS+prT))) {
			indiv[[j]]$state <- "T"
		} else if ((randNum > (prS+prT)) & (randNum <= (prS+prT+prD))) {
			indiv[[j]]$state <- "D"
		} else if ((randNum > (prS+prT+prD)) & (randNum <= (prS+prT+prD+prA))) {
			indiv[[j]]$state <- "A"
		} else if ((randNum > (prS+prT+prD+prA)) & (randNum <= (prS+prT+prD+prA+prU))) {
			indiv[[j]]$state <- "U"
		} else {
			indiv[[j]]$state <- "P"
		}
		setTxtProgressBar(pb, j)
	}
	close(pb)
	# Checking state attributes:
	# stateVector <- sapply(indiv, function(x) x$state)
	# stateVector

	# Add an attribute to keep track of number of days of latent infection:
	for (j in 1:N) {
		indiv[[j]]$daysLatent <- 0
	}

	## Make empty vectors to record population statistics over time:
	time <- seq(numIter)

	# Number of individuals belonging to each state:
	numS <- NA * time # number of individuals in S state
	numE <- NA * time # number of individuals in E state
	numT <- NA * time # number of individuals in T state
	numD <- NA * time # number of individuals in D state
	numA <- NA * time # number of individuals in A state
	numU <- NA * time # number of individuals in U state
	numP <- NA * time # number of individuals in P state

	NAll <- NA * time # total number of individuals
	N2_10 <- NA * time # total number of individuals (ages 2-10)
	N0_5 <- NA * time # total number of individuals (ages 0-5)
	N5_10 <- NA * time # total number of individuals (ages 5-10)
	N10_15 <- NA * time # total number of individuals (ages 10-15)
	N15Plus <- NA * time # total number of individuals (ages 15+)

	clinIncAll <- 0 * time # clinical incidence (all ages)
	clinInc2_10 <- 0 * time # clinical incidence (ages 2-10)
	clinInc0_5 <- 0 * time # clinical incidence (ages 0-5)
	clinInc5_10 <- 0 * time # clinical incidence (ages 5-10)
	clinInc10_15 <- 0 * time # clinical incidence (ages 10-15)
	clinInc15Plus <- 0 * time # clinical incidence (ages 15+)

	slidPosAll <- 0 * time # clinical incidence (all ages)
	slidPos2_10 <- 0 * time # clinical incidence (ages 2-10)
	slidPos0_5 <- 0 * time # clinical incidence (ages 0-5)
	slidPos5_10 <- 0 * time # clinical incidence (ages 5-10)
	slidPos10_15 <- 0 * time # clinical incidence (ages 10-15)
	slidPos15Plus <- 0 * time # clinical incidence (ages 15+)

	cat("Initialization complete!\n")

	## Run the simulation:
	for (i in 1:numIter) {

		# Determine which individuals are currently in each state:
		isS <- which( sapply(indiv, function(x) x$state) == "S")
		isE <- which( sapply(indiv, function(x) x$state) == "E")
		isT <- which( sapply(indiv, function(x) x$state) == "T")
		isD <- which( sapply(indiv, function(x) x$state) == "D")
		isA <- which( sapply(indiv, function(x) x$state) == "A")
		isU <- which( sapply(indiv, function(x) x$state) == "U")
		isP <- which( sapply(indiv, function(x) x$state) == "P")

		## Update number of individuals belonging to each state:
		numS[i] <- length(isS)
		numE[i] <- length(isE)
		numT[i] <- length(isT)
		numD[i] <- length(isD)
		numA[i] <- length(isA)
		numU[i] <- length(isU)
		numP[i] <- length(isP)
		NAll[i] <- length(indiv)
		N <- NAll[i]

		## Determine which individuals are in each age group:
		is2_10 <- which( ( sapply(indiv, function(x) x$age) >= 2 ) &
			           ( sapply(indiv, function(x) x$age) < 10 ) )
		is0_5 <- which( sapply(indiv, function(x) x$age) < 5 )
		is5_10 <- which( ( sapply(indiv, function(x) x$age) >= 5 ) &
			           ( sapply(indiv, function(x) x$age) < 10 ) )
		is10_15 <- which( ( sapply(indiv, function(x) x$age) >= 10 ) &
			            ( sapply(indiv, function(x) x$age) < 15 ) )
		is15Plus <- which( sapply(indiv, function(x) x$age) >= 15 )

		## Update number of individuals belonging to each age group:
		N2_10[i] <- length(is2_10)
		N0_5[i] <- length(is0_5)
		N5_10[i] <- length(is5_10)
		N10_15[i] <- length(is10_15)
		N15Plus[i] <- length(is15Plus)

		## Of asymptomatic diseased individuals, determine which are in each
		## age group:
		isA2_10 <- which( ( sapply(indiv, function(x) x$age) >= 2 ) &
			            ( sapply(indiv, function(x) x$age) < 10 ) &
					( sapply(indiv, function(x) x$state) == "A" ) )
		isA0_5 <- which( ( sapply(indiv, function(x) x$age) < 5 ) &
				     ( sapply(indiv, function(x) x$state) == "A" ) )
		isA5_10 <- which( ( sapply(indiv, function(x) x$age) >= 5 ) &
			            ( sapply(indiv, function(x) x$age) < 10 ) &
					( sapply(indiv, function(x) x$state) == "A" ) )
		isA10_15 <- which( ( sapply(indiv, function(x) x$age) >= 10 ) &
			             ( sapply(indiv, function(x) x$age) < 15 ) &
					 ( sapply(indiv, function(x) x$state) == "A" ) )
		isA15Plus <- which( ( sapply(indiv, function(x) x$age) >= 15 ) &
					  ( sapply(indiv, function(x) x$state) == "A" ) )

		## Calculate slide positivity by age group:
		stateVector <- sapply(indiv, function(x) x$state)
		prDetectAMicVector <- sapply(indiv, function(x) x$prDetectAMic)

		slidPos2_10[i] <- sum(stateVector[is2_10] == "T") +
					sum(stateVector[is2_10] == "D") +
					if (length(isA2_10)>0) { sum(sapply(prDetectAMicVector[isA2_10], function(x) rbinom(n=1, size=1, prob=x))) } else { 0 }
		slidPos0_5[i] <- sum(stateVector[is0_5] == "T") +
				     sum(stateVector[is0_5] == "D") +
				     if (length(isA0_5)>0) { sum(sapply(prDetectAMicVector[isA0_5], function(x) rbinom(n=1, size=1, prob=x))) } else { 0 }
		slidPos5_10[i] <- sum(stateVector[is5_10] == "T") +
					sum(stateVector[is5_10] == "D") +
					if (length(isA5_10)>0) { sum(sapply(prDetectAMicVector[isA5_10], function(x) rbinom(n=1, size=1, prob=x))) } else { 0 }
		slidPos10_15[i] <- sum(stateVector[is10_15] == "T") +
					 sum(stateVector[is10_15] == "D") +
					 if (length(isA10_15)>0) { sum(sapply(prDetectAMicVector[isA10_15], function(x) rbinom(n=1, size=1, prob=x))) } else { 0 }
		slidPos15Plus[i] <- sum(stateVector[is15Plus] == "T") +
					  sum(stateVector[is15Plus] == "D") +
					  if (length(isA15Plus)>0) { sum(sapply(prDetectAMicVector[isA15Plus], function(x) rbinom(n=1, size=1, prob=x))) } else { 0 }
		slidPosAll[i] <- slidPos0_5[i] + slidPos5_10[i] + slidPos10_15[i] + slidPos15Plus[i]

		# Determine which individuals die during the current time step:
		for (j in 1:N) {
			randNum <- runif(1)
			if(randNum <= mu) {
				indiv[[j]]$alive <- 0
				householdSize[indiv[[j]]$house] <- householdSize[indiv[[j]]$house] - 1
			}
		}

		## Transitions from S compartment:
		for (j in isS) {
			if (indiv[[j]]$alive == 1) {
				# Extract information for jth individual.
				lambda <- indiv[[j]]$lambda

				# Draw a random number between 1 and 0 for each susceptible
				# individual.
				randNum <- runif(1)

				## Latent infection (S -> E):
				# If the random number is less than lambda, that individual
				# develops a latent infection (E) in the next time step.
				if (randNum <= lambda) {
					indiv[[j]]$state <- "E"
					indiv[[j]]$daysLatent <- 1
				}
			}
		}

		## Transitions from E compartment:
		for (j in isE) {
			if (indiv[[j]]$alive == 1) {
				if (indiv[[j]]$daysLatent < dE) {
					indiv[[j]]$daysLatent <- indiv[[j]]$daysLatent + 1
				} else {
					# Extract information for jth individual.
					phi <- indiv[[j]]$phi

					# Draw a random number between 1 and 0 for each
					# latently-infected individual.
					randNum <- runif(1)

					## Treated clinical infection (E -> T):
					# If the random number is less than phi*fT, that
					# individual develops a treated clinical infection (T)
					# in the next time step.
					if (randNum <= phi*fT) {
						indiv[[j]]$state <- "T"
						indiv[[j]]$daysLatent <- 0

						# Keep track of clinical incidence in different age groups:
						clinIncAll[i] <- clinIncAll[i] + 1
						if ((indiv[[j]]$age >= 2) & (indiv[[j]]$age < 10)) {
							clinInc2_10[i] <- clinInc2_10[i] + 1
						}
						if (indiv[[j]]$age < 5) {
							clinInc0_5[i] <- clinInc0_5[i] + 1
						} else if ((indiv[[j]]$age >= 5) & (indiv[[j]]$age < 10)) {
							clinInc5_10[i] <- clinInc5_10[i] + 1
						} else if ((indiv[[j]]$age >= 10) & (indiv[[j]]$age < 15)) {
							clinInc10_15[i] <- clinInc10_15[i] + 1
						} else if (indiv[[j]]$age >= 15) {
							clinInc15Plus[i] <- clinInc15Plus[i] + 1
						}
					}

					## Untreated clinical infection (E -> D):
					# If the random number is greater than phi*fT and less
					# than phi, that individual develops an untreated
					# clinical infection (D) in the next time step.
					if ((randNum > phi*fT) & (randNum <= phi)) {
						indiv[[j]]$state <- "D"
						indiv[[j]]$daysLatent <- 0

						# Keep track of clinical incidence in different age groups:
						clinIncAll[i] <- clinIncAll[i] + 1
						if ((indiv[[j]]$age >= 2) & (indiv[[j]]$age < 10)) {
							clinInc2_10[i] <- clinInc2_10[i] + 1
						}
						if (indiv[[j]]$age < 5) {
							clinInc0_5[i] <- clinInc0_5[i] + 1
						} else if ((indiv[[j]]$age >= 5) & (indiv[[j]]$age < 10)) {
							clinInc5_10[i] <- clinInc5_10[i] + 1
						} else if ((indiv[[j]]$age >= 10) & (indiv[[j]]$age < 15)) {
							clinInc10_15[i] <- clinInc10_15[i] + 1
						} else if (indiv[[j]]$age >= 15) {
							clinInc15Plus[i] <- clinInc15Plus[i] + 1
						}
					}

					## Asymptomatic infection (E -> A):
					# If the random number is greater than phi, that
					# individual develops an asymptomatic infection (A) in
					# the next time step.
					if (randNum > phi) {
						indiv[[j]]$state <- "A"
						indiv[[j]]$daysLatent <- 0
					}
				}
			}
		}

		## Transitions from T compartment:
		for (j in isT) {
			if (indiv[[j]]$alive == 1) {
				# Draw a random number between 1 and 0 for each treated individual.
				randNum <- runif(1)

				## Prophylactic protection (T -> P):
				# If the random number is less than 1/dT, that individual enters the
				# phase of prophylactic protection (P) in the next time step.
				if (randNum <= (1/dT)) {
					indiv[[j]]$state <- "P"
				}
			}
		}

		## Transitions from D compartment:
		for (j in isD) {
			if (indiv[[j]]$alive == 1) {
				# Draw a random number between 1 and 0 for each treated individual.
				randNum <- runif(1)

				## Progression from diseased to asymptomatic (D -> A):
				# If the random number is less than 1/dD, that individual enters the
				# phase of asymptomatic patent infection (A) in the next time step.
				if (randNum <= (1/dD)) {
					indiv[[j]]$state <- "A"
				}
			}
		}

		## Transitions from A compartment:
		for (j in isA) {
			if (indiv[[j]]$alive == 1) {
				# Extract information for jth individual.
				lambda <- indiv[[j]]$lambda
				phi <- indiv[[j]]$phi

				# Draw a random number between 1 and 0 for each susceptible
				# individual.
				randNum <- runif(1)

				## Treated clinical infection (A -> T):
				# If the random number is less than phi*fT*lambda, that
				# individual develops a treated clinical infection (T) in
				# the next time step.
				if (randNum <= phi*fT*lambda) {
					indiv[[j]]$state <- "T"

					# Keep track of clinical incidence in different age groups:
					clinIncAll[i] <- clinIncAll[i] + 1
					if ((indiv[[j]]$age >= 2) & (indiv[[j]]$age < 10)) {
						clinInc2_10[i] <- clinInc2_10[i] + 1
					}
					if (indiv[[j]]$age < 5) {
						clinInc0_5[i] <- clinInc0_5[i] + 1
					} else if ((indiv[[j]]$age >= 5) & (indiv[[j]]$age < 10)) {
						clinInc5_10[i] <- clinInc5_10[i] + 1
					} else if ((indiv[[j]]$age >= 10) & (indiv[[j]]$age < 15)) {
						clinInc10_15[i] <- clinInc10_15[i] + 1
					} else if (indiv[[j]]$age >= 15) {
						clinInc15Plus[i] <- clinInc15Plus[i] + 1
					}
				}

				## Untreated clinical infection (A -> D):
				# If the random number is greater than phi*fT*lambda and
				# less than phi*lambda, that individual develops an
				# untreated clinical infection (D) in the next time step.
				if ((randNum > phi*fT*lambda) & (randNum <= phi*lambda)) {
					indiv[[j]]$state <- "D"

					# Keep track of clinical incidence in different age groups:
					clinIncAll[i] <- clinIncAll[i] + 1
					if ((indiv[[j]]$age >= 2) & (indiv[[j]]$age < 10)) {
						clinInc2_10[i] <- clinInc2_10[i] + 1
					}
					if (indiv[[j]]$age < 5) {
						clinInc0_5[i] <- clinInc0_5[i] + 1
					} else if ((indiv[[j]]$age >= 5) & (indiv[[j]]$age < 10)) {
						clinInc5_10[i] <- clinInc5_10[i] + 1
					} else if ((indiv[[j]]$age >= 10) & (indiv[[j]]$age < 15)) {
						clinInc10_15[i] <- clinInc10_15[i] + 1
					} else if (indiv[[j]]$age >= 15) {
						clinInc15Plus[i] <- clinInc15Plus[i] + 1
					}
				}

				## Progression to asymptomatic sub-patent infection (A -> U):
				# If the random number is greater than phi*lambda and less
				# than (phi*lambda + 1/dA), that individual develops an asymptomatic
				# infection (A) in the next time step.
				if ((randNum > phi*lambda) & (randNum <= (phi*lambda + (1/dA)))) {
					indiv[[j]]$state <- "U"
				}
			}
		}

		## Transitions from U compartment:
		for (j in isU) {
			if (indiv[[j]]$alive == 1) {
				# Extract information for jth individual.
				lambda <- indiv[[j]]$lambda
				phi <- indiv[[j]]$phi

				# Draw a random number between 1 and 0 for each susceptible
				# individual.
				randNum <- runif(1)

				## Treated clinical infection (U -> T):
				# If the random number is less than phi*fT*lambda, that
				# individual develops a treated clinical infection (T) in
				# the next time step.
				if (randNum <= phi*fT*lambda) {
					indiv[[j]]$state <- "T"

					# Keep track of clinical incidence in different age groups:
					clinIncAll[i] <- clinIncAll[i] + 1
					if ((indiv[[j]]$age >= 2) & (indiv[[j]]$age < 10)) {
						clinInc2_10[i] <- clinInc2_10[i] + 1
					}
					if (indiv[[j]]$age < 5) {
						clinInc0_5[i] <- clinInc0_5[i] + 1
					} else if ((indiv[[j]]$age >= 5) & (indiv[[j]]$age < 10)) {
						clinInc5_10[i] <- clinInc5_10[i] + 1
					} else if ((indiv[[j]]$age >= 10) & (indiv[[j]]$age < 15)) {
						clinInc10_15[i] <- clinInc10_15[i] + 1
					} else if (indiv[[j]]$age >= 15) {
						clinInc15Plus[i] <- clinInc15Plus[i] + 1
					}
				}

				## Untreated clinical infection (U -> D):
				# If the random number is greater than phi*fT*lambda and
				# less than phi*lambda, that individual develops an
				# untreated clinical infection (D) in the next time step.
				if ((randNum > phi*fT*lambda) & (randNum <= phi*lambda)) {
					indiv[[j]]$state <- "D"

					# Keep track of clinical incidence in different age groups:
					clinIncAll[i] <- clinIncAll[i] + 1
					if ((indiv[[j]]$age >= 2) & (indiv[[j]]$age < 10)) {
						clinInc2_10[i] <- clinInc2_10[i] + 1
					}
					if (indiv[[j]]$age < 5) {
						clinInc0_5[i] <- clinInc0_5[i] + 1
					} else if ((indiv[[j]]$age >= 5) & (indiv[[j]]$age < 10)) {
						clinInc5_10[i] <- clinInc5_10[i] + 1
					} else if ((indiv[[j]]$age >= 10) & (indiv[[j]]$age < 15)) {
						clinInc10_15[i] <- clinInc10_15[i] + 1
					} else if (indiv[[j]]$age >= 15) {
						clinInc15Plus[i] <- clinInc15Plus[i] + 1
					}
				}

				## Asymptomatic infection (U -> A):
				# If the random number is greater than phi*lambda and
				# less than lambda, that individual develops a patent
				# asymptomatic infection (A) in the next time step.
				if ((randNum > phi*lambda) & (randNum <= lambda)) {
					indiv[[j]]$state <- "A"
				}

				## Progression to asymptomatic sub-patent infection (U -> S):
				# If the random number is greater than lambda and less
				# than (lambda + 1/dU), that individual returns to the susceptible
				# state (S) in the next time step.
				if ((randNum > lambda) & (randNum <= (lambda + (1/dU)))) {
					indiv[[j]]$state <- "S"
				}
			}
		}

		## Transitions from P compartment:
		for (j in isP) {
			if (indiv[[j]]$alive == 1) {
				# Draw a random number between 1 and 0 for each treated individual.
				randNum <- runif(1)

				## Prophylactic protection (P -> S):
				# If the random number is less than 1/dP, that individual returns to
				# the susceptible state (S) in the next time step.
				if (randNum <= (1/dP)) {
					indiv[[j]]$state <- "S"
				}
			}
		}

		## Update ages by time step (one day):
		for (j in 1:N) {
			if (indiv[[j]]$alive == 1) {
				indiv[[j]]$age <- indiv[[j]]$age + (1/365)
			}
		}

		## Immunity:
		# Update values for:
		# 1. Pre-erythrocytic immunity (IB, reduces the probability of infection
		#    following an infectious challenge)
		# 2. Acquired clinical immunity (ICA, reduces the probability of clinical
		#    disease, acquired from previous exposure)
		# 3. Maternal clinical immunity (ICM, reduces the probability of clinical
		#    disease, acquired maternally)
		# 4. Detection immunity (ID, a.k.a. blood-stage immunity, reduces the
		#    probability of detection and reduces infectiousness to mosquitoes)
		for (j in 1:N) {
			if (indiv[[j]]$alive == 1) {
				epsilon <- indiv[[j]]$epsilon
				lambda <- indiv[[j]]$lambda
				indiv[[j]]$IB <- indiv[[j]]$IB + (epsilon/(epsilon*uB + 1)) - (indiv[[j]]$IB)*(1/dB)
				indiv[[j]]$ICA <- indiv[[j]]$ICA + (lambda/(lambda*uC + 1)) - (indiv[[j]]$ICA)*(1/dC)
				indiv[[j]]$ICM <- indiv[[j]]$ICM - (indiv[[j]]$ICM)*(1/dM)
				indiv[[j]]$ID <- indiv[[j]]$ID + (lambda/(lambda*uD + 1)) - (indiv[[j]]$ID)*(1/dID)
			}
		}

		# Lambda (the force of infection) is calculated for each individual. It
		# varies according to age and biting heterogeneity group:
		for (j in 1:N) {
			if (indiv[[j]]$alive == 1) {
				zita <- indiv[[j]]$bitingHet
				psi <- psiHouse[indiv[[j]]$house]
				a <- indiv[[j]]$age
				IB <- indiv[[j]]$IB
				epsilon <- epsilon0 * zita * (1 - rho * exp(-a/a0)) * psi
				b <- b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
				indiv[[j]]$epsilon <- epsilon
				indiv[[j]]$lambda <- epsilon*b
			}
		}

		# Phi (the probability of acquiring clinical disease upon infection) is
		# also calculated for each individual. It varies according to immune
		# status:
		for (j in 1:N) {
			if (indiv[[j]]$alive == 1) {
				ICA <- indiv[[j]]$ICA
				ICM <- indiv[[j]]$ICM
				indiv[[j]]$phi <- phi0 * (phi1 + ((1 - phi1)/(1 + ((ICA+ICM)/IC0)^kappaC)))
			}
		}

		# q (the probability that an asymptomatic infection is detected by
		# microscopy) is also calculated for each individual, as well as the
		# probability of detection by PCR for asymptomatic infections in states
		# A (patent) and U (subpatent). This also varies according to immune
		# status:
		for (j in 1:N) {
			if (indiv[[j]]$alive == 1) {
				a <- indiv[[j]]$age
				ID <- indiv[[j]]$ID
				fD <- 1 - ((1 - fD0)/(1 + (a/aD)^gammaD))
				q <- d1 + ((1 - d1)/(1 + (fD*(ID/ID0)^kappaD)*fD))
				indiv[[j]]$prDetectAMic <- q
				indiv[[j]]$prDetectAPCR <- q^alphaA
				indiv[[j]]$prDetectUPCR <- q^alphaU
			}
		}

		# Incorporate new births at the end of the time step:
		numNewBirths <- rbinom(n = 1, size = N, prob = mu)
		if (numNewBirths > 0) {
			is18_22 <- which( ( sapply(indiv, function(x) x$age) >= 18 ) &
			                  ( sapply(indiv, function(x) x$age) < 22 ) )
			meanICA18_22 <- mean(sapply(indiv[is18_22], function(x) x$ICA))
			for (j in 1:numNewBirths) {
				# Assign new birth to one of the samllest houses:
				smallestHouse <- which(householdSize==min(householdSize))
				house <- sample(smallestHouse, 1)
				householdSize[house] <- householdSize[house] + 1
				psi <- psiHouse[house]
				bitingHet <- rlnorm(n = 1, meanlog = -sigma2/2, sdlog = sqrt(sigma2))
				epsilon <- epsilon0*bitingHet*(1-rho)*psi
				phi <- phi0 * (phi1 + ((1 - phi1)/(1 + (PM*meanICA18_22/IC0)^kappaC)))
				indiv[[N+j]] <- list(state="S", age=0, alive=1, bitingHet=bitingHet,
							IB=0, ICA=0, ICM=PM*meanICA18_22, ID=0, daysLatent=0,
							epsilon=epsilon, lambda=epsilon*b0, phi=phi,
							house=house,
							prDetectAMic=1, prDetectAPCR=1, prDetectUPCR=1)
			}
		}

		# Remove dead individuals:
		isDead <- which( sapply(indiv, function(x) x$alive) == 0)
		if (length(isDead > 0)) {
			indiv <- indiv[-isDead]
		}

		# Print year annually:
		if (i%%365==0) { cat("Year", i/365, "of", numIter/365, "\n") }

		############################
		## Visualize cases:       ##
		############################

		# if (i%%7==0) { # Visualize cases weekly
		if (i%%1==0) { # Visualize cases daily
			# 1. Plot risk surface on a 1x1 grid:
			# First, calculate risk surface for the map grid:
			riskMap <- matrix(rep(0,101*101),101)
			for (j in 1:101) {
				longMap <- (j-1)/100 # Convert index to 0-1 map scale
				for (k in 1:101) {
					latMap <- (k-1)/100 # Convert index to 0-1 map scale
					for (l in 1:numBreedingSites) {
						# Compute cumulative risk due to all of the breeding sites:
						riskMap[j,k] <- ( riskMap[j,k] + (1/(sigmaBreedingSite[l]*sqrt(2*pi)))
								    * exp(-((longMap-longBreedingSite[l])^2 +
		            	                       (latMap-latBreedingSite[l])^2)
		                  	                / (2*sigmaBreedingSite[l]^2)))
					}
				}
			}
			numContours <- 9
			contour(x = seq(0, 1, length.out = nrow(riskMap)),
				  y = seq(0, 1, length.out = ncol(riskMap)),
				  z= riskMap, drawlabels=FALSE, nlevels=numContours, lwd=2,
				  col=brewer.pal(numContours, "Blues"), axes=FALSE)

			# 2. Plot breeding sites
			# for (j in 1:numBreedingSites) {
			#	points(longBreedingSite[j], latBreedingSite[j], col="red", lwd=3)
			# }

			# 3. Plot houses
			houseEdgeSize <- 0.015
			for (j in 1:numHouses) {
				polygon(x=c(longHouse[j]+houseEdgeSize/2, longHouse[j]+houseEdgeSize/2, longHouse[j]-houseEdgeSize/2, longHouse[j]-houseEdgeSize/2, longHouse[j]),
				        y=c(latHouse[j]+houseEdgeSize/2, latHouse[j]-houseEdgeSize/2, latHouse[j]-houseEdgeSize/2, latHouse[j]+houseEdgeSize/2, latHouse[j]+houseEdgeSize),
					  col="brown", border="black", lwd=2)
			}

			# 3. Plot T (green), D (red), A (orange), U (yellow) as they come
			#    up (around house):
			radCases <- 0.025
			for (j in 1:numHouses) {
				isTDAUHouseJ <- which( ( sapply(indiv, function(x) x$house) == j) &
							     ( ( sapply(indiv, function(x) x$state) == "T") |
								 ( sapply(indiv, function(x) x$state) == "D") |
								 ( sapply(indiv, function(x) x$state) == "A") |
								 ( sapply(indiv, function(x) x$state) == "U") ) )
				numCasesHouseJ <- length(isTDAUHouseJ)
				if (numCasesHouseJ > 0) {
					for (k in 1:numCasesHouseJ) {
						if (indiv[[isTDAUHouseJ[k]]]$state == "T") {
							points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="green", lwd=2, pch=16)
						} else if (indiv[[isTDAUHouseJ[k]]]$state == "D") {
							points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="red", lwd=2, pch=16)
						} else if (indiv[[isTDAUHouseJ[k]]]$state == "A") {
							points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="orange", lwd=2, pch=16)
						} else if (indiv[[isTDAUHouseJ[k]]]$state == "U") {
							points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="yellow", lwd=2, pch=16)
						}
					}
				}
			}
		}
	}

	## Return the simulation results as a data frame:
	simResults <- data.frame(time, numS, numE, numT, numD, numA, numU, numP,
	  clinIncAll, clinInc2_10, clinInc0_5, clinInc5_10, clinInc10_15, clinInc15Plus,
	  slidPosAll, slidPos2_10, slidPos0_5, slidPos5_10, slidPos10_15, slidPos15Plus,
	  NAll, N2_10, N0_5, N5_10, N10_15, N15Plus)
	return(simResults)
}



#########################
## Run the simulation: ##
#########################


# numIter <- 365 # Simulation over 1 year
sim1 <- malaria_ibm(theta=theta, numIter=numIter) # Runs the simulation

##########################################################
## Plot number of individuals in each state over time:  ##
##########################################################
ggplot(sim1, aes(x = time, y = sim1, color = State)) +
	geom_line(aes(y = numS, col = "S"), size = 1.2) +
	geom_line(aes(y = numE, col = "E"), size = 1.2) +
	geom_line(aes(y = numT, col = "T"), size = 1.2) +
	geom_line(aes(y = numD, col = "D"), size = 1.2) +
	geom_line(aes(y = numA, col = "A"), size = 1.2) +
	geom_line(aes(y = numU, col = "U"), size = 1.2) +
	geom_line(aes(y = numP, col = "P"), size = 1.2) +
	labs(x = "Time (days)", y = "Prevalence")

## Plot clinical incidence by age group over time:
ggplot(sim1, aes(x = time, y = sim1, color = "Age group")) +
	geom_line(aes(y = clinIncAll/NAll, col = "All"), size = 1.2) +
	geom_line(aes(y = clinInc2_10/N2_10, col = "2-10"), size = 1.2) +
	geom_line(aes(y = clinInc0_5/N0_5, col = "0-5"), size = 1.2) +
	geom_line(aes(y = clinInc5_10/N5_10, col = "5-10"), size = 1.2) +
	geom_line(aes(y = clinInc10_15/N10_15, col = "10-15"), size = 1.2) +
	geom_line(aes(y = clinInc15Plus/N15Plus, col = "15+"), size = 1.2) +
	labs(x = "Time (days)", y = "Clinical incidence")

## Plot slide positivity by age group over time:
ggplot(sim1, aes(x = time, y = sim1, color = "Age group")) +
	geom_line(aes(y = slidPosAll, col = "All"), size = 1.2) +
	geom_line(aes(y = slidPos2_10, col = "2-10"), size = 1.2) +
	geom_line(aes(y = slidPos0_5, col = "0-5"), size = 1.2) +
	geom_line(aes(y = slidPos5_10, col = "5-10"), size = 1.2) +
	geom_line(aes(y = slidPos10_15, col = "10-15"), size = 1.2) +
	geom_line(aes(y = slidPos15Plus, col = "15+"), size = 1.2) +
	labs(x = "Time (days)", y = "Slide positive")

## Plot population size by age group over time:
ggplot(sim1, aes(x = time, y = sim1, color = "Age group")) +
	geom_line(aes(y = NAll, col = "All"), size = 1.2) +
	geom_line(aes(y = N2_10, col = "2-10"), size = 1.2) +
	geom_line(aes(y = N0_5, col = "0-5"), size = 1.2) +
	geom_line(aes(y = N5_10, col = "5-10"), size = 1.2) +
	geom_line(aes(y = N10_15, col = "10-15"), size = 1.2) +
	geom_line(aes(y = N15Plus, col = "15+"), size = 1.2) +
	labs(x = "Time (days)", y = "Number of individuals")

#################################
## Save simulation output:     ##
#################################

head(sim1)
write.table(sim1,file="MalariaIBMsim.csv",sep=",",row.names=F,col.names=T)
