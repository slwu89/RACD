# IF THE DELAYS DONT WORK, JUST USE EXPONENTIAL EIP

# inbound jump size to EL

# instantaneous hazards for EL
haz_EL_mort <- muEL*(1 + ((EL+LL)/K))
haz_EL_2LL <- 1/durEL

# jump probabilities
p_EL_0 <- exp(-(haz_EL_mort + haz_EL_2LL)*dt)
p_EL_mort <- (1 - p_EL_0)*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)) # death
p_EL_2LL <- (1 - p_EL_0)*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)) # to late-instar

# instantaneous hazards for LL
haz_LL_mort <- muLL*(1 + gamma*((EL+LL)/K))
haz_LL_2PL <- 1/durLL

# jump probabilities
p_LL_0 <- exp(-(haz_LL_mort + haz_LL_2PL)*dt)
p_LL_mort <- (1 - p_LL_0)*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)) # death
p_LL_2PL <- (1 - p_LL_0)*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)) # to pupae

# instantaneous hazards for PL
haz_PL_mort <- muPL
haz_PL_2PL <- 1/durPL

# jump probabilities
p_PL_0 <- exp(-(haz_PL_mort + haz_PL_2PL)*dt)
p_PL_mort <- (1 - p_PL_0)*(haz_PL_mort / (haz_PL_mort + haz_PL_2PL)) # death
p_PL_2PL <- (1 - p_PL_0)*(haz_PL_2PL / (haz_PL_mort + haz_PL_2PL)) # to susceptible

# instantaneous hazards for SV
haz_SV_mort <-  muVCom
haz_SV_inf <- lambdaV

# jump probabilities
p_SV_0 <- exp(-(haz_SV_mort + haz_SV_inf)*dt)
p_SV_mort <- (1 - p_SV_0)*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)) # death
p_SV_2EV <- (1 - p_SV_0)*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)) # to incubating

# instantaneous hazards for EV
haz_EV_mort <-  muVCom
haz_EV_inc <- lambdaV

# jump probabilities
p_SV_0 <- exp(-(haz_SV_mort + haz_SV_inf)*dt)
p_SV_mort <- (1 - p_SV_0)*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)) # death
p_SV_2EV <- (1 - p_SV_0)*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)) # to incubating

# dEL <- betaCom*NV - muEL*(1 + ((EL+LL)/K))*EL - EL/durEL
# dLL <- EL/durEL - muLL*(1 + gamma*((EL+LL)/K))*LL - LL/durLL
# dPL <- LL/durLL - muPL*PL - PL/durPL
# dSV <- 0.5*PL/durPL - lambdaV*SV - muVCom*SV
# dEV <- lambdaV*SV - lambdaV*SVLag*exp(-muVCom*durEV) - muVCom*EV
# dIV <- lambdaV*SVLag*exp(-muVCom*durEV) - muVCom*IV