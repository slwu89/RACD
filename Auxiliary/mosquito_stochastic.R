################################################################################
#   Exponential distribution of EIP
#   Deterministic Model
################################################################################

rm(list=ls());gc()

# make a single node
make_node <- function(){

  node <- list()

  # new eggs are generated from a conditionally independent Poisson process
  node$EL_new <- 0
  storage.mode(node$EL_new) <- "integer"

  # probabilities & transitions for early-stage instars (names are destinations of jumps)
  node$EL_probs <- rep(0,3)
  node$EL_probs <- setNames(node$EL_probs,c("EL","D","LL"))
  node$EL_transitions <- rep(0,3)
  node$EL_transitions <- setNames(node$EL_transitions,c("EL","D","LL"))
  storage.mode(node$EL_transitions) <- "integer"

  # probabilities & transitions for late-stage instars
  node$LL_probs <- rep(0,3)
  node$LL_probs <- setNames(node$LL_probs,c("LL","D","PL"))
  node$LL_transitions <- rep(0,3)
  node$LL_transitions <- setNames(node$LL_transitions,c("LL","D","PL"))
  storage.mode(node$LL_transitions) <- "integer"

  # probabilities & transitions for pupae
  node$PL_probs <- rep(0,4)
  node$PL_probs <- setNames(node$PL_probs,c("PL","D","SV_F","SV_M"))
  node$PL_transitions <- rep(0,4)
  node$PL_transitions <- setNames(node$PL_transitions,c("PL","D","SV_F","SV_M"))
  storage.mode(node$PL_transitions) <- "integer"

  # probabilities & transitions for susceptible vectors
  node$SV_probs <- rep(0,3)
  node$SV_probs <- setNames(node$SV_probs,c("SV","D","EV"))
  node$SV_transitions <- rep(0,3)
  node$SV_transitions <- setNames(node$SV_transitions,c("SV","D","EV"))
  storage.mode(node$SV_transitions) <- "integer"

  # probabilities & transitions for incubating vectors
  node$EV_probs <- rep(0,3)
  node$EV_probs <- setNames(node$EV_probs,c("EV","D","IV"))
  node$EV_transitions <- rep(0,3)
  node$EV_transitions <- setNames(node$EV_transitions,c("EV","D","IV"))
  storage.mode(node$EV_transitions) <- "integer"

  # probabilities & transitions for infectious vectors
  node$IV_probs <- rep(0,2)
  node$IV_probs <- setNames(node$IV_probs,c("IV","D"))
  node$IV_transitions <- rep(0,2)
  node$IV_transitions <- setNames(node$IV_transitions,c("IV","D"))
  storage.mode(node$IV_transitions) <- "integer"

  # state space
  node$EL <- 0L
  node$LL <- 0L
  node$PL <- 0L
  node$SV <- 0L
  node$EV <- 0L
  node$IV <- 0L

  list2env(node,hash=TRUE)
}

# execute an Euler-step for the model
# rates -> probabilities are defined by discretisation of the underlying
# CTMC system (competing hazards)
euler_step <- function(node,pars,tnow,dt){
  with(pars,{

    ########################################
    # INTERVENTION-DEPENDENT PARAMETERS
    ########################################

    delta <- 1/(tau1+tau2) # Inverse of gonotrophic cycle without ITNs/IRS
    e_ov <- beta*(exp(muV/delta)-1)/muV # Number of eggs per oviposition per mosquito

    ## Derived parameters which depend on intervention status:
    if(tnow > pars$time_ITN_on){
      ITNcov_t <- pars$ITNcov
    } else {
      ITNcov_t <- 0
    }
    if(tnow > pars$time_IRS_on){
      IRScov_t <- pars$IRScov
    } else {
      IRScov_t <- 0
    }

    # zCom: Probability of a mosquito being repelled from an ITN or IRS-treated house:
    c0 <- 1 - ITNcov_t - IRScov_t + ITNcov_t*IRScov_t
    cITN <- ITNcov_t - ITNcov_t*IRScov_t
    cIRS <- IRScov_t - ITNcov_t*IRScov_t
    cCom <- ITNcov_t*IRScov_t
    rCom <- rIRS + (1-rIRS)*rITN
    sCom  <- (1-rIRS)*sITN*sIRS
    zCom <- Q0*cITN*phiB*rITN + Q0*cIRS*phiI*rIRS + Q0*cCom*(phiI-phiB)*rIRS + Q0*cCom*phiB*rCom

    # deltaCom: Inverse of gonotrophic cycle length with ITNs & IRS:
    deltaCom <- 1/(tau1/(1-zCom) + tau2)

    # wCom: Probability that a surviving mosquito succeeds in feeding during a single attempt:
    wCom <- 1 - Q0 + Q0*c0 + Q0*cITN*(1-phiB+phiB*sITN) + Q0*cIRS*(1-phiI+phiI*sIRS) + Q0*cCom*((phiI-phiB)*sIRS + 1-phiI + phiB*sCom)

    # muVCom: Female mosquito death rate in presence of ITNs & IRS:
    p10 <- exp(-muV*tau1)
    p1Com <- p10*wCom/(1 - zCom*p10)
    p2 <- exp(-muV*tau2)
    pCom <- (p1Com*p2)^deltaCom
    muVCom <- -log(pCom)

    # betaCom: Eggs laid per day by female mosquitoes in presence of ITNs & IRS:
    betaCom <- e_ov*muVCom/(exp(muVCom/deltaCom) - 1)

    ########################################
    # EARLY-STAGE LARVAL INSTARS (EL)
    ########################################

    # inbound oviposition to EL
    NV <- sum(node$SV,node$EV,node$IV)
    node$EL_new <- rpois(n = 1,lambda = betaCom*NV*dt)

    # instantaneous hazards for EL
    haz_EL_mort <- muEL*(1 + ((node$EL+node$LL)/K))
    haz_EL_2LL <- 1/durEL

    # jump probabilities
    node$EL_probs[["EL"]] <- exp(-(haz_EL_mort + haz_EL_2LL)*dt)
    node$EL_probs[["D"]] <- (1 - node$EL_probs[["EL"]])*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)) # death
    node$EL_probs[["LL"]] <- (1 - node$EL_probs[["EL"]])*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)) # to late-instar

    # jump sizes
    node$EL_transitions[c("EL","D","LL")] <- as.vector(rmultinom(n=1,size=node$EL,prob=node$EL_probs))

    ########################################
    # LATE-STAGE LARVAL INSTARS (LL)
    ########################################

    # instantaneous hazards for LL
    haz_LL_mort <- muLL*(1 + gamma*((node$EL+node$LL)/K))
    haz_LL_2PL <- 1/durLL

    # jump probabilities
    node$LL_probs[["LL"]] <- exp(-(haz_LL_mort + haz_LL_2PL)*dt)
    node$LL_probs[["D"]] <- (1 - node$LL_probs[["LL"]])*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)) # death
    node$LL_probs[["PL"]] <- (1 - node$LL_probs[["LL"]])*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)) # to pupae

    # jump sizes
    node$LL_transitions[c("LL","D","PL")] <- as.vector(rmultinom(n=1,size=node$LL,prob=node$LL_probs))

    ########################################
    # PUPAE (PL)
    ########################################

    # instantaneous hazards for PL
    haz_PL_mort <- muPL
    haz_PL_2SV_F <- (1/durPL)*0.5
    haz_PL_2SV_M <- (1/durPL)*0.5

    # jump probabilities
    node$PL_probs[["PL"]] <- exp(-(haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)*dt)
    node$PL_probs[["D"]] <- (1 - node$PL_probs[["PL"]])*(haz_PL_mort / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)) # death
    node$PL_probs[["SV_F"]] <- (1 - node$PL_probs[["PL"]])*(haz_PL_2SV_F / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)) # to susceptible female
    node$PL_probs[["SV_M"]] <- (1 - node$PL_probs[["PL"]])*(haz_PL_2SV_M / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)) # to susceptible males

    # jump sizes
    node$PL_transitions[c("PL","D","SV_F","SV_M")] <- as.vector(rmultinom(n=1,size=node$PL,prob=node$PL_probs))

    ########################################
    # SUSCEPTIBLE VECTORS (SV)
    ########################################

    # instantaneous hazards for SV
    haz_SV_mort <-  muVCom
    haz_SV_inf <- lambdaV

    # jump probabilities
    node$SV_probs[["SV"]] <- exp(-(haz_SV_mort + haz_SV_inf)*dt)
    node$SV_probs[["D"]] <- (1 - node$SV_probs[["SV"]])*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)) # death
    node$SV_probs[["EV"]] <- (1 - node$SV_probs[["SV"]])*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)) # to incubating

    # jump sizes
    node$SV_transitions[c("SV","D","EV")] <- as.vector(rmultinom(n=1,size=node$SV,prob=node$SV_probs))

    ########################################
    # INCUBATING VECTORS (EV)
    ########################################

    # instantaneous hazards for EV
    haz_EV_mort <-  muVCom
    haz_EV_inc <- 1/durEV

    # jump probabilities
    node$EV_probs[["EV"]] <- exp(-(haz_EV_mort + haz_EV_inc)*dt)
    node$EV_probs[["D"]] <- (1 - node$EV_probs[["EV"]])*(haz_EV_mort / (haz_EV_mort + haz_EV_inc)) # death
    node$EV_probs[["IV"]] <- (1 - node$EV_probs[["EV"]])*(haz_EV_inc / (haz_EV_mort + haz_EV_inc)) # to infectious

    # jump sizes
    node$EV_transitions[c("EV","D","IV")] <- as.vector(rmultinom(n=1,size=node$EV,prob=node$EV_probs))

    ########################################
    # INFECTIOUS VECTORS (IV)
    ########################################

    # instantaneous hazards for IV
    haz_IV_mort <- muVCom

    # jump probabilities
    node$IV_probs[["IV"]] <- exp(-haz_IV_mort*dt)
    node$IV_probs[["D"]] <- (1 - node$IV_probs[["IV"]])

    # jump sizes
    node$IV_transitions[c("IV","D")] <- as.vector(rmultinom(n=1,size=node$IV,prob=node$IV_probs))

    ########################################
    # UPDATE POPULATION
    ########################################

    node$EL <- node$EL_transitions[["EL"]] + node$EL_new
    node$LL <- node$LL_transitions[["LL"]] + node$EL_transitions[["LL"]]
    node$PL <- node$PL_transitions[["PL"]] + node$LL_transitions[["PL"]]
    node$SV <- node$SV_transitions[["SV"]] + node$PL_transitions[["SV_F"]]
    node$EV <- node$EV_transitions[["EV"]] + node$SV_transitions[["EV"]]
    node$IV <- node$IV_transitions[["IV"]] + node$EV_transitions[["IV"]]

  })
}

## Model parameters:
theta <- list(
	## Mosquito life cycle parameters:
	beta = 21.19, # Number of eggs laid per day by female mosquito
	muEL = 0.034, # Early larval instar daily mortality
	muLL = 0.035, # Late larval instar daily mortality
	muPL = 0.25, # Pupal daily mortality
	durEL = 6.64, # Duration of early instar stage
	durLL = 3.72, # Duration of late instar stage
	durPL = 0.64, # Duration of pupal stage
	durEV = 10, # Duration of latent period in mosquito (days)
	gamma = 13.25, # Effect of density-dependence on late instars relative to early instars
	tau1 = 0.68, # Time spent foraging for a blood meal at 0% ITN coverage
	tau2 = 2.32, # Time spent resting and ovipositing by a mosquito

	## Intervention parameters (variable):
	ITNcov = 0.5, # ITN coverage
	IRScov = 0.5, # IRS coverave
	time_ITN_on = 1e3, # When ITNs are applied (days)
	time_IRS_on = 1e3, # When IRS is applied (days)

	## Species-specific parameters:
	## An. gambiae:
	muV = 1/7.6, # Adult mosquito daily mortality
	Q0 = 0.92, # Human blood index
	phiB = 0.89, # Proportion of bites on a person while they are in bed
	phiI = 0.97, # Proportion of bites on a person while they are indoors
	rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
	sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
	rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
	sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

	## An. arabiensis:
	# muV = 1/7.6, # Adult mosquito daily mortality
	# Q0 = 0.71, # Human blood index
	# phiB = 0.90, # Proportion of bites on a person while they are in bed
	# phiI = 0.96, # Proportion of bites on a person while they are indoors
	# rITN = 0.48, # Probability of mosquito repeating a feeding attempt due to IRS
	# sITN = 0.39, # Probability of mosquito feeding and surviving in presence of ITNs
	# rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
	# sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

	## An. funestus:
	# muV = 1/8.9, # Adult mosquito daily mortality
	# Q0 = 0.94, # Human blood index
	# phiB = 0.90, # Proportion of bites on a person while they are in bed
	# phiI = 0.98, # Proportion of bites on a person while they are indoors
	# rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
	# sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
	# rIRS = 0.63, # Probability of mosquito repeating a feeding attempt due to IRS
	# sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

	## Additional transmission parameters:
	f0 = 1/3, # Daily biting rate by mosquitoes on animals and humans
	epsilon0 = 10/365, # Daily entomological inolculation rate
	iH_eq = 0.45, # Equilibrium malaria prevalence in humans
	NH_eq = 200, # Equilibrium human population size
	bV = 0.15 # Probability of transmission from human to vector per infectious bite
)

with(theta,{

  a0 <<- Q0*f0 # Human biting rate at equilibrium
  lambdaV <<- a0*iH_eq*bV # Force of infection in mosquitoes at equilibrium

  IV_eq <<- 100
  EV_eq <<- durEV*IV_eq*muV
  SV_eq <<- ((durEV*IV_eq*muV) + ((durEV^2)*IV_eq*(muV^2))) / (durEV*lambdaV)
  NV_eq <<- IV_eq + EV_eq + SV_eq

  omega_1 <<- -0.5 * ( (gamma*(muLL/muEL)) - (durEL/durLL) + ((gamma-1)*muLL*durEL) )
  omega_2 <<- sqrt(
    (0.25 * (( (gamma*(muLL/muEL)) - (durEL/durLL) + ((gamma-1)*muLL*durEL) )^2)) + (gamma * ((beta*muLL*durEL) / (2*muEL*muV*durLL*(1 + (durPL*muPL)))) )
  )
  omega <<- omega_1 + omega_2

  PL_eq <<- 2*durPL*muV*NV_eq
  LL_eq <<- 2*muV*durLL*(1 + (durPL*muPL))*NV_eq
  EL_eq <<- 2*omega*muV*durLL*(1 + (durPL*muPL))*NV_eq

  K <<- (NV_eq*2*durLL*muV*(1 + (durPL*muPL))*gamma*(omega+1)) / ((omega/(muLL*durEL)) - (1/(muLL*durLL)) - 1)
})


# # equilibrium solutions
# with(theta,{
#
#   ## Derived parameters:
#   # NV <<- SV + EV + IV # Total mosquito population size
#   delta <<- 1/(tau1+tau2) # Inverse of gonotrophic cycle without ITNs/IRS
#   e_ov <<- beta*(exp(muV/delta)-1)/muV # Number of eggs per oviposition per mosquito
#   b_omega <<- gamma*muLL/muEL - durEL/durLL + (gamma-1)*muLL*durEL
#   omega <<- -0.5*b_omega + sqrt(0.25*b_omega^2 + gamma*beta*muLL*durEL/(2*muEL*muV*durLL*(1+durPL*muPL)))
#
#   a0 <<- Q0*f0 # Human biting rate at equilibrium
#
#   lambdaV <<- a0*iH_eq*bV # Force of infection in mosquitoes at equilibrium
#
#
#   iV_eq <<- lambdaV*exp(-muV*durEV)/(lambdaV + muV)
#   sV_eq <<- iV_eq*muV/(lambdaV*exp(-muV*durEV))
#   eV_eq <<- 1 - sV_eq - iV_eq
#
#   NV_eq <<- epsilon0*NH_eq/(iV_eq*a0)
#
#   EL_eq <<- 2*omega*muV*durLL*(1 + muPL*durPL)*NV_eq
#   LL_eq <<- 2*muV*durLL*(1 + muPL*durPL)*NV_eq
#   PL_eq <<- 2*muV*durPL*NV_eq
#
#   SV_eq <<- sV_eq*NV_eq
#   EV_eq <<- eV_eq*NV_eq
#   IV_eq <<- iV_eq*NV_eq
#
#   K <<- 2*NV_eq*muV*durLL*(1 + muPL*durPL)*gamma*(omega+1)/(omega/(muLL*durEL) - 1/(muLL*durLL) - 1) # Larval carrying capacity
#
# })

theta1 <- c(theta,K=K,lambdaV=lambdaV)

node <- make_node()
node$EL <- as.integer(EL_eq)
node$LL <- as.integer(LL_eq)
node$PL <- as.integer(PL_eq)
node$SV <- as.integer(SV_eq)
node$EV <- as.integer(EV_eq)
node$IV <- as.integer(IV_eq)

tmax <- 500
dt <- 0.05
time <- seq(from=1,to=tmax,by=dt)

# sampling grid
sample_grid <- tsamp <- c(0,seq(from=10,to = tmax,by = 5))
sample_pop <- matrix(0,nrow=6,ncol=length(sample_grid),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(sample_grid)))

sample_pop["EL",1] <- node$EL
sample_pop["LL",1] <- node$LL
sample_pop["PL",1] <- node$PL
sample_pop["SV",1] <- node$SV
sample_pop["EV",1] <- node$EV
sample_pop["IV",1] <- node$IV
sample_grid <- sample_grid[-1]

# run simulation
pb <- txtProgressBar(min = 1,max = length(time))
for(t in 1:length(time)){

  # euler step
  euler_step(node = node,pars = theta1,tnow = time[t],dt = dt)

  # sample the population (done at the very end of the time-step, because its not part of the dynamics)
  if(time[t] == sample_grid[1]){

    sample_pop["EL",as.character(sample_grid[1])] <- node$EL
    sample_pop["LL",as.character(sample_grid[1])] <- node$LL
    sample_pop["PL",as.character(sample_grid[1])] <- node$PL
    sample_pop["SV",as.character(sample_grid[1])] <- node$SV
    sample_pop["EV",as.character(sample_grid[1])] <- node$EV
    sample_pop["IV",as.character(sample_grid[1])] <- node$IV
    sample_grid <- sample_grid[-1]

  }
  setTxtProgressBar(pb = pb,value = t)
}

sample_pop_t <- t(sample_pop)
par(mfrow=c(1,2))
ylim <- max(sample_pop_t[,4:6])
plot(x = tsamp,y = sample_pop_t[,"SV"],col="blue",lwd=2,ylim=c(0,ylim),type="l")
lines(x = tsamp,y = sample_pop_t[,"EV"],col="green",lwd=2)
lines(x = tsamp,y = sample_pop_t[,"IV"],col="red",lwd=2)

ylim <- max(sample_pop_t[,5:6])
plot(x = tsamp,y = sample_pop_t[,"EV"],col="green",lwd=2,ylim=c(0,ylim),type="l")
lines(x = tsamp,y = sample_pop_t[,"IV"],col="red",lwd=2)
par(mfrow=c(1,1))

tail(sample_pop_t)
mean(sample_pop_t[,"IV"])



################################################################################
# run ensemble of simulations
################################################################################

nruns <- 100
tmax <- 250
dt <- 0.5
time <- seq(from=1,to=tmax,by=dt)

# sampling grid
tsamp <- c(0,seq(from=10,to = tmax,by = 1))
sample_pop <- array(0,dim=c(6,length(tsamp),nruns),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(tsamp),paste0(1:nruns)))


pb <- txtProgressBar(min = 1,max = nruns)
for(i in 1:nruns){

  sample_grid <- tsamp

  # make the node
  node <- make_node()
  node$EL <- as.integer(EL_eq)
  node$LL <- as.integer(LL_eq)
  node$PL <- as.integer(PL_eq)
  node$SV <- as.integer(SV_eq)
  node$EV <- as.integer(EV_eq)
  node$IV <- as.integer(IV_eq)

  # record output
  sample_pop["EL",1,i] <- node$EL
  sample_pop["LL",1,i] <- node$LL
  sample_pop["PL",1,i] <- node$PL
  sample_pop["SV",1,i] <- node$SV
  sample_pop["EV",1,i] <- node$EV
  sample_pop["IV",1,i] <- node$IV
  sample_grid <- sample_grid[-1]

  for(t in 1:length(time)){

    # euler step
    euler_step(node = node,pars = theta1,tnow = time[t],dt = dt)

    # sample the population (done at the very end of the time-step, because its not part of the dynamics)
    if(time[t] == sample_grid[1]){

      sample_pop["EL",as.character(sample_grid[1]),i] <- node$EL
      sample_pop["LL",as.character(sample_grid[1]),i] <- node$LL
      sample_pop["PL",as.character(sample_grid[1]),i] <- node$PL
      sample_pop["SV",as.character(sample_grid[1]),i] <- node$SV
      sample_pop["EV",as.character(sample_grid[1]),i] <- node$EV
      sample_pop["IV",as.character(sample_grid[1]),i] <- node$IV
      sample_grid <- sample_grid[-1]

    }
  }

  setTxtProgressBar(pb,i)
}

# mean_E <- colMeans(sample_pop["EL",,])
# mean_L <- colMeans(sample_pop["LL",,])
# mean_P <- colMeans(sample_pop["PL",,])

mean_pops <- apply(X = sample_pop,MARGIN = c(1,2),FUN = mean)

mean_SV <- rowMeans(sample_pop["SV",,])
mean_EV <- rowMeans(sample_pop["EV",,])
mean_IV <- rowMeans(sample_pop["IV",,])

quant_SV <- apply(X = sample_pop["SV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_EV <- apply(X = sample_pop["EV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_IV <- apply(X = sample_pop["IV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})

plot_datEV <- data.frame(time=tsamp,EV=mean_EV,EV_l=quant_EV[1,],EV_h=quant_EV[2,])
plot_datIV <- data.frame(time=tsamp,IV=mean_IV,IV_l=quant_IV[1,],IV_h=quant_IV[2,])

library(ggplot2)

ggplot() +
  geom_line(data=plot_datEV,aes(x=time,y=EV),color="steelblue") +
  geom_ribbon(data=plot_datEV,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="steelblue") +
  geom_line(data=plot_datIV,aes(x=time,y=IV),color="firebrick3") +
  geom_ribbon(data=plot_datIV,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick3") +
  theme_bw()


################################################################################
# run ensemble of simulations (with better equilibrium)
################################################################################

# test our numerical/analytic "exact" equilibria
dt <- 1
eq <- calc_eq(theta = theta1,dt = dt,IV = IV_eq,lambdaV = lambdaV)
theta2 <- theta1
theta2$K <- eq$K_eq

# ensemble parameters
nruns <- 100
tmax <- 1e3
dt <- 1
time <- seq(from=1,to=tmax,by=dt)

# sampling grid
tsamp <- c(0,seq(from=10,to = tmax,by = 1))
sample_pop <- array(0,dim=c(6,length(tsamp),nruns),dimnames=list(c("EL","LL","PL","SV","EV","IV"),paste0(tsamp),paste0(1:nruns)))


pb <- txtProgressBar(min = 1,max = nruns)
for(i in 1:nruns){

  sample_grid <- tsamp

  # make the node
  node <- make_node()
  node$EL <- as.integer(eq$EL_eq)
  node$LL <- as.integer(eq$LL_eq)
  node$PL <- as.integer(eq$PL_eq)
  node$SV <- as.integer(eq$SV_eq)
  node$EV <- as.integer(eq$EV_eq)
  node$IV <- as.integer(IV_eq)

  # record output
  sample_pop["EL",1,i] <- node$EL
  sample_pop["LL",1,i] <- node$LL
  sample_pop["PL",1,i] <- node$PL
  sample_pop["SV",1,i] <- node$SV
  sample_pop["EV",1,i] <- node$EV
  sample_pop["IV",1,i] <- node$IV
  sample_grid <- sample_grid[-1]

  for(t in 1:length(time)){

    # euler step
    euler_step(node = node,pars = theta2,tnow = time[t],dt = dt)

    # sample the population (done at the very end of the time-step, because its not part of the dynamics)
    if(time[t] == sample_grid[1]){

      sample_pop["EL",as.character(sample_grid[1]),i] <- node$EL
      sample_pop["LL",as.character(sample_grid[1]),i] <- node$LL
      sample_pop["PL",as.character(sample_grid[1]),i] <- node$PL
      sample_pop["SV",as.character(sample_grid[1]),i] <- node$SV
      sample_pop["EV",as.character(sample_grid[1]),i] <- node$EV
      sample_pop["IV",as.character(sample_grid[1]),i] <- node$IV
      sample_grid <- sample_grid[-1]

    }
  }

  setTxtProgressBar(pb,i)
}

# plot the output
mean_SV <- rowMeans(sample_pop["SV",,])
mean_EV <- rowMeans(sample_pop["EV",,])
mean_IV <- rowMeans(sample_pop["IV",,])

quant_SV <- apply(X = sample_pop["SV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_EV <- apply(X = sample_pop["EV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})
quant_IV <- apply(X = sample_pop["IV",,],MARGIN = 1,FUN = function(x){
  quantile(x,probs = c(0.05,0.95))
})

plot_datEV <- data.frame(time=tsamp,EV=mean_EV,EV_l=quant_EV[1,],EV_h=quant_EV[2,])
plot_datIV <- data.frame(time=tsamp,IV=mean_IV,IV_l=quant_IV[1,],IV_h=quant_IV[2,])

library(ggplot2)

ggplot() +
  geom_line(data=plot_datEV,aes(x=time,y=EV),color="steelblue") +
  geom_ribbon(data=plot_datEV,aes(x=time,ymin=EV_l,ymax=EV_h),alpha=0.35,fill="steelblue") +
  geom_line(data=plot_datIV,aes(x=time,y=IV),color="firebrick3") +
  geom_ribbon(data=plot_datIV,aes(x=time,ymin=IV_l,ymax=IV_h),alpha=0.35,fill="firebrick3") +
  theme_bw()
