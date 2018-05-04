#' Initialize State for RACD Model
#'
#' do something
#'
#' @param theta named vector of parameters (see \code{\link{RACD_Parameters}})
#'
#' @export
RACD_Setup1 <- function(theta){

  cat("begin initializing RACD model state\n")

  # return from within the 'with' statement; constructs local scope
  with(as.list(theta),{
    browser()
    mu = 1/(meanAge*365) # Daily death rate as a function of mean age in years

    # geographic parameters
    numHouses <- round(N/meanNumPeoplePerHouse)
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

  	# At the beginning of the simulation, all individuals are alive:
  	for (j in 1:N) {
  		indiv[[j]]$alive <- TRUE
  	}

    cat("begin assigning house coordinates\n")
  	# Randomly assign house coordinates:
  	longHouse <- rep(0, numHouses) # Vector of house longitudes
  	latHouse <- rep(0, numHouses) # Vector of house latitudes
  	longHouse[1] <- runif(1) # First house longitude coordinate
  	latHouse[1] <- runif(1) # First house latitude coordinate

    # progress bar
    progress_bar = utils::txtProgressBar(min=1,max=numHouses,initial = 1, style=3)

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
      utils::setTxtProgressBar(progress_bar,j)
  	}
    close(progress_bar)
    cat("\n")
    cat("done assigning house coordinates\n")

    cat("begin assigning individuals to houses\n")
  	# Randomly assign individuals to houses:
  	householdSize <- rep(0, numHouses) # Vector of household sizes
  	for (j in 1:N) {
  		# Assign individuals to one of the samllest houses:
  		smallestHouse <- which(householdSize==min(householdSize))
  		if (length(smallestHouse)==1) { indiv[[j]]$house <- smallestHouse }
  		else { indiv[[j]]$house <- sample(smallestHouse, 1) }
  		householdSize[indiv[[j]]$house] <- householdSize[indiv[[j]]$house] + 1
  	}
    cat("end assigning individuals to houses\n")

    cat("begin assigning breeding site coordinates\n")
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
    cat("done assigning breeding site coordinates\n")

    cat("begin calculating risk surface psi\n")
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
    cat("done calculating risk surface psi\n")

    cat("begin calculating individual immune status\n")
  	# Randomly assign biting heterogeneity attributes:
  	# (Biting heterogeneity is sampled from a normal distribution with mean
  	# 1 and standard deviation sigma)
  	# (Referred to as zeta in Griffin et al. (2014))
  	for (j in 1:N) {
  		indiv[[j]]$bitingHet <- rlnorm(n = 1, meanlog = -sigma2/2, sdlog = sqrt(sigma2))
  	}

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

    # progress bar
    progress_bar = utils::txtProgressBar(min=1,max=N,initial = 1, style=3)

  	for (j in 1:N) {
  		zeta <- indiv[[j]]$bitingHet
  		psi <- psiHouse[indiv[[j]]$house]
  		a <- indiv[[j]]$age
  		# Calculate initial immunity levels from their differential equations:
      initI <- immune_initial(a=a,zeta=zeta,psi=psi,theta=theta)
  		IB <- initI[["IB"]]; ID <- initI[["ID"]]; ICA <- initI[["ICA"]]
  		initI20 <- immune_initial(a=20, zeta=1, psi=1, theta=theta)
  		ICM <- initI[["ICA"]] * exp(-a/dM)
  		epsilon <- epsilon0 * zeta * (1 - rho * exp(-a/a0)) * psi
  		b <- b0*(b1 + ((1-b1)/(1 + (IB/IB0)^kappaB)))
  		lambda <- epsilon*b
  		indiv[[j]]$IB <- IB
  		indiv[[j]]$ID <- ID
  		indiv[[j]]$ICA <- ICA
  		indiv[[j]]$ICM <- ICM
  		indiv[[j]]$epsilon <- epsilon
  		indiv[[j]]$lambda <- lambda
      utils::setTxtProgressBar(progress_bar,j)
  	}
    close(progress_bar)
    cat("\n")
    cat("done calculating individual immune status\n")

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

  	# For each invididual, use the ODE transmission model to determine the
  	# probability that they are in each state given their age and EIR
  	# heterogeneity attributes:
  	cat("Initializing population...\n")
  	pb <- txtProgressBar(min = 1, max = N, initial = 1, style=3)
  	for (j in 1:N) {
  		a <- indiv[[j]]$age
  		zeta <- indiv[[j]]$bitingHet
  		psi <- psiHouse[indiv[[j]]$house]
  		prStateVector <- infection_initial(a=a, zeta=zeta, psi=psi, theta=theta)
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
    cat("\n")
    cat("Done initializing population...\n")

  	# Add an attribute to keep track of number of days of latent infection:
  	for (j in 1:N) {
  		indiv[[j]]$daysLatent <- 0
  	}

    return(list(
      humans = indiv,
      breedingSites = mapply(longBreedingSite,latBreedingSite,sigmaBreedingSite,FUN = function(x,y,z){list(x=x,y=y,sigma=z)},SIMPLIFY = FALSE),
      houses = mapply(longHouse,latHouse,householdSize,psiHouse,FUN = function(x,y,z,w){list(x=x,y=y,size=z,psi=w)},SIMPLIFY = FALSE)
    ))
  })
}
