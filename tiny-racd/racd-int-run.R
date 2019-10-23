###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Marshall Lab (https://www.marshalllab.com)
#   Sean Wu (slwu89@berkeley.edu)
#   June 2019
#   Run RACD model w/interventions
#
###############################################################################

rm(list=ls());gc()

library(ggplot2)
library(gridExtra)
library(reshape2)
library(here)
source(here::here("racd-setup.R"))

# assume a landscape has been made with appropriate data frames; see landscape.R for how to do this
radius <- quantile(dmat_dwell,probs = (1:10)*1e-1)[[1]]

meanEIR <- 0.05
N <- 500

RACD_init <- RACD_Setup(N = N,EIR_mean = meanEIR,xy_d = dwell_df,xy_a = aqua_df,theta = RACD_theta)

# calculate immunity parameters for imported cases
imm_import <- imported_immune(EIR = meanEIR,theta = RACD_theta)

RACD_theta <- c(RACD_theta,imm_import,import_rate=0.01) 

# compile the simulation
sourceCpp(here::here("intervention-src/main.cpp"))

# run the simulation
RACD_out <- tiny_racd(humans_param = RACD_init$humans,
                      house_param = RACD_init$houses,
                      mosy_param = RACD_init$mosy,
                      theta = RACD_theta,
                      tmax = 365*30,
                      int_type = 1,
                      tstart = 365*10,
                      tend = 365*20,
                      tdelay = 5,
                      dmat = dmat_dwell,
                      radius = radius,
                      p_index = 0.95,
                      p_neighbor = 0.9,
                      prog_bar = TRUE)

out_stateage <- melt(RACD_out$state_hist$state_age)
colnames(out_stateage) <- c("state","age","count","day")

out_cinc <- melt(RACD_out$state_hist$clinical_incidence)
colnames(out_cinc) <- c("day","age","count")

out_mosy <- melt(RACD_out$state_hist$mosquito)
colnames(out_mosy) <- c("day","state","count")

out_eir_b <- rbind(
  data.frame(par="b",day=seq_along(RACD_out$state_hist$b),RACD_out$state_hist$b),
  data.frame(par="EIR",day=seq_along(RACD_out$state_hist$b),RACD_out$state_hist$EIR)
)

ggplot(data = out_stateage) +
  geom_line(aes(x=day,y=count,color=state),alpha=0.75) +
  facet_wrap(. ~ age) +
  theme_bw()

ggplot(data = out_cinc) +
  geom_line(aes(x=day,y=count,color=age),alpha=0.75) +
  theme_bw()

ggplot(data = out_mosy) +
  geom_line(aes(x=day,y=count,color=state),alpha=0.75) +
  theme_bw()

ggplot(data = out_eir_b) +
  geom_line(aes(x=day,y=mean,color=par)) +
  geom_linerange(aes(x=day,ymin=pmax(mean-sqrt(var),0),ymax=mean+sqrt(var),color=par),alpha=0.005) +
  facet_wrap(. ~ par) +
  theme_bw()
