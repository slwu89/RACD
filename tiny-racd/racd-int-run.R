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
                      tend = 365*12,
                      dmat = dmat_dwell,
                      radius = radius,
                      prog_bar = TRUE)


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

# grid.arrange(plot_outS,plot_outA,plot_outFOI,plot_outM)
invisible(gc())

# transmission to humans

plot_EIR <- ggplot(data = melt(RACD_out$trans[,c("time","EIR_mean","EIR_var")],id.vars = c("time","EIR_mean","EIR_var"))) +
  geom_ribbon(aes(x=time,ymin=pmax(EIR_mean - sqrt(EIR_var),0),ymax=EIR_mean + sqrt(EIR_var)),alpha=0.75) +
  geom_line(aes(x=time,y=EIR_mean)) +
  xlab("Time (days)") + ylab("EIR (entomological inoculation rate)") +
  theme_bw()

# plot_FOI_h <- ggplot(data = melt(RACD_out$trans[,c("time","lambda_h_mean","lambda_h_var")],id.vars = c("time","lambda_h_mean","lambda_h_var"))) +
#   geom_ribbon(aes(x=time,ymin=pmax(lambda_h_mean - sqrt(lambda_h_var),0),ymax=lambda_h_mean + sqrt(lambda_h_var)),alpha=0.75) +
#   geom_line(aes(x=time,y=lambda_h_mean)) +
#   xlab("Time (days)") + ylab("FOI on humans") +
#   theme_bw()

plot_b <- ggplot(data = melt(RACD_out$trans[,c("time","b_mean","b_var")],id.vars = c("time","b_mean","b_var"))) +
  geom_ribbon(aes(x=time,ymin=pmax(b_mean - sqrt(b_var),0),ymax=b_mean + sqrt(b_var)),alpha=0.75) +
  geom_line(aes(x=time,y=b_mean)) +
  xlab("Time (days)") + ylab("b (mosquito to human transmission efficiency)") +
  theme_bw()

# grid.arrange(plot_EIR,plot_FOI_h,nrow=1)
invisible(gc())
# table(RACD_out$intervention)

grid.arrange(plot_outS,plot_outA,plot_outFOI,plot_outM,plot_EIR,plot_b)
