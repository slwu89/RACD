###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   December 2017
#   Setup Simulation
#
###############################################################################

#' Initialize State for RACD Model
#'
#' Initialize initial conditions for the RACD simulation model
#'
#' @param xy_h two column matrix of house coordinates
#' @param xy_b two column matrix of breeding site coordinates
#' @param theta named vector of parameters (see \code{\link{RACD_Parameters}})
#'
#' @examples
#' \dontrun{
#' library(RACD)
#' library(RACDaux)
#' library(tidyverse)
#' library(spatstat)
#' xy_h <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
#' xy_b <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
#' theta <- RACD_Parameters()
#' init <- RACD_Setup(as.matrix.ppx(xy_h),as.matrix.ppx(xy_b),theta)
#' outfile = "/Users/slwu89/Desktop/log_trans.csv"
#' RACD_Simulation(365,theta,init$humans,init$houses,123,outfile)
#' state = RACDaux::RACD_StateVector(outfile)
#' state %>% as.tibble %>% gather(state,value,-time) %>% ggplot(aes(x=time,y=value,color=state)) + geom_line() + theme_bw()
#' }
#'
#' @export
RACD_Setup <- function(xy_h, xy_b, theta){

  with(as.list(theta),{

    mu <- 1/(meanAge*365) # daily mortality st mean age does not change

    numHouses <- nrow(xy_h)
    numBreedingSites <- nrow(xy_b)

    indiv <- vector(mode="list",length=N)

    # Randomly assign age attributes:
    # (Age is sampled from an exponential distribution with mean equal to
    # the mean age in the country being considered)
    for (j in 1:N) {
      indiv[[j]]$age <- rexp(n = 1, rate = 1/meanAge)
    }

    # At the beginning of the simulation, all humans are alive:
    for (j in 1:N) {
      indiv[[j]]$alive <- TRUE
    }

    # Randomly assign individuals to houses:
    householdSize <- rep(0, numHouses) # Vector of household sizes
    for (j in 1:N) {
      # Assign individuals to one of the samllest houses:
      smallestHouse <- which(householdSize==min(householdSize))
      if(length(smallestHouse)==1){
        indiv[[j]]$house <- smallestHouse
      } else {
        indiv[[j]]$house <- sample(smallestHouse, 1)
      }
      householdSize[indiv[[j]]$house] <- householdSize[indiv[[j]]$house] + 1
    }

    cat("begin calculating risk surface psi\n")
    # Define geographical risk surface due to breeding sites:
    # We model the risk due to breeding sites at each household as the sum
    # of multivariate normal distributions centered at each breeding site
    # with the standard deviation of the normal distribution being equal to
    # the distance between the breeding site and the nearest household.
    # 1. First, we calculate the standard deviation of the risk surface for
    #    each breeding site.
    sigmaBreedingSite <- rep(0,numBreedingSites)
    pb <- txtProgressBar(min=1,max=(numBreedingSites+numHouses),initial = 1, style=3)
    for(j in 1:numBreedingSites){
      distHousesBreedingSites <- as.matrix(dist(x=rbind(xy_b[j,],xy_h)))[1,2:(nrow(xy_h)+1)] # vector of distance from breeding site j to all houses
      sigmaBreedingSite[j] <- min(distHousesBreedingSites) # sigma(sd) of breeding site k is nearest distance to house
      setTxtProgressBar(pb,j)
    }

    # 2. Next, calculate the contribution of each breeding site to the
    #    relative risk value for each house.
    psiHouse <- rep(0,numHouses)
    for(j in 1:numHouses){
      for(k in 1:numBreedingSites){
        d_jk <- as.matrix(dist(x = rbind(xy_h[j,],xy_b[k,])))[1,2] # distance between house j and breeding site k
        psiHouse[j] <- psiHouse[j] + dnorm(x = d_jk,mean = 0,sd = sigmaBreedingSite[k])
      }
      setTxtProgressBar(pb,(j+numBreedingSites))
    }
    close(pb)
    # 3. Finally, we normalize the psiHouse values so they have a mean of 1.
    psiHouse <- psiHouse*numHouses/sum(psiHouse) # normalize to have mean=1, sum=numHouses
    cat("done calculating risk surface psi\n")
    cat("\n")

    cat("begin calculating individual immune status\n")
    pb = txtProgressBar(min=1,max=N,initial = 1, style=3)

    for (j in 1:N) {

      # Randomly assign biting heterogeneity attributes:
      # (Biting heterogeneity is sampled from a log-normal distribution with mean
      # 1 and standard deviation sigma)
      # (Referred to as zeta in Griffin et al. (2014))
      indiv[[j]]$bitingHet <- zeta <- rlnorm(n = 1, meanlog = -sigma2/2, sdlog = sqrt(sigma2))

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

      # Phi (the probability of acquiring clinical disease upon infection) is
      # also calculated for each individual. It varies according to immune
      # status:
      indiv[[j]]$phi <- phi0 * (phi1 + ((1 - phi1)/(1 + ((ICA+ICM)/IC0)^kappaC)))

      # q (the probability that an asymptomatic infection is detected by
      # microscopy) is also calculated for each individual, as well as the
      # probability of detection by PCR for asymptomatic infections in states
      # A (patent) and U (subpatent). This also varies according to immune
      # status:
      fD <- 1 - ((1 - fD0)/(1 + (a/aD)^gammaD))
      q <- d1 + ((1 - d1)/(1 + (fD*(ID/ID0)^kappaD)*fD))
      indiv[[j]]$prDetectAMic <- q
      indiv[[j]]$prDetectAPCR <- q^alphaA
      indiv[[j]]$prDetectUPCR <- q^alphaU

      setTxtProgressBar(pb,j)
    }
    close(pb)
    cat("done calculating individual immune status\n")
    cat("\n")

    # For each invididual, use the ODE transmission model to determine the
    # probability that they are in each state given their age and EIR
    # heterogeneity attributes:
    cat("Initializing individual state probabilities...\n")
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
      indiv[[j]]$state <- sample(x=c("S","T","D","A","U","P"), size=1, prob = c(prS,prT,prD,prA,prU,prP))

      # Add an attribute to keep track of number of days of latent infection:
      indiv[[j]]$daysLatent <- 0

      setTxtProgressBar(pb, j)
    }
    close(pb)
    cat("Done initializing individual state probabilities...\n")
    cat("\n")

    return(list(
      humans = indiv,
      breedingSites = mapply(xy_b[,1],xy_b[,2],sigmaBreedingSite,FUN = function(x,y,z){list(x=x,y=y,sigma=z)},SIMPLIFY = FALSE),
      houses = mapply(xy_h[,1],xy_h[,2],householdSize,psiHouse,FUN = function(x,y,z,w){list(x=x,y=y,size=z,psi=w)},SIMPLIFY = FALSE)
    ))
  })
}
