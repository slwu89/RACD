###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Marshall Lab (https://www.marshalllab.com)
#   Sean Wu (slwu89@berkeley.edu)
#   April 2019
#   Run RACD model
#
###############################################################################

rm(list=ls());gc()

library(here)
source(here("racd-setup.R"))

# just a test; dead simple model
dwell_df <- data.frame(x=c(1,2),y=c(1,2),psi=rep(0.5,2))
aqua_df <- data.frame(x=1,y=1)
RACD_init <- RACD_Setup(N = 2000,EIR_mean = 0.05,xy_d = dwell_df,xy_a = aqua_df,theta = RACD_theta)

# model diagnostics
library(ggplot2)
library(gridExtra)

plot_c <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$c})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("P(infect | blood feed on host)") +
  ggtitle("Population Infectivity to Mosquitos")

plot_phi <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$phi})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("P(clinical malaria | effective colonization by merozoites)") +
  ggtitle("Clinical incidence probability")

plot_age <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$age})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("age") +
  ggtitle("Age distribution of population")

plot_state <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$state})))+
  geom_histogram(aes(x=x),stat="count") +
  theme_bw() +
  xlab("State") +
  ggtitle("Initial Distribution over States")

grid.arrange(plot_c,plot_phi,plot_age,plot_state)

# compile the simulation
sourceCpp(here("racd-src/main.cpp"))

# run the simulation
RACD_out <- tiny_racd(humans_param = RACD_init$humans,
                      house_param = RACD_init$houses,
                      mosy_param = RACD_init$mosy,
                      theta = RACD_theta,
                      tmax = 365*3)

library(reshape2)

ggplot(data = melt(RACD_out$state,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

ggplot(data = melt(RACD_out$age,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

ggplot(data = melt(RACD_out$clinical_incidence,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

ggplot(data = melt(RACD_out$mosy,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()
