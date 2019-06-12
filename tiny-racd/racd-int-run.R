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
radius <- quantile(dmat_dwell)[[3]]

meanEIR <- 0.01
N <- 500

RACD_init <- RACD_Setup(N = N,EIR_mean = meanEIR,xy_d = dwell_df,xy_a = aqua_df,theta = RACD_theta)

# compile the simulation
sourceCpp(here::here("intervention-src/main.cpp"))

# run the simulation
RACD_out <- tiny_racd(humans_param = RACD_init$humans,
                      house_param = RACD_init$houses,
                      mosy_param = RACD_init$mosy,
                      theta = RACD_theta,
                      tmax = 365*20,
                      int_type = 1,
                      dmat = dmat_dwell,
                      radius = radius)


# state variables
plot_outS <- ggplot(data = melt(RACD_out$state,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.75) +
  theme_bw()

plot_outA <- ggplot(data = melt(RACD_out$age,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

plot_outFOI <- ggplot(data = melt(RACD_out$trans[,1:2],id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

plot_outM <- ggplot(data = melt(RACD_out$mosy[,!names(RACD_out$mosy) == "S"],id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.65) +
  theme_bw()

grid.arrange(plot_outS,plot_outA,plot_outFOI,plot_outM)


# table(RACD_out$intervention)