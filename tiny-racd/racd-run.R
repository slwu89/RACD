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

meanEIR <- 0.01
N <- 500

RACD_init <- RACD_Setup(N = N,EIR_mean = meanEIR,xy_d = dwell_df,xy_a = aqua_df,theta = RACD_theta)

# model diagnostics
library(ggplot2)
library(gridExtra)

# basic diagnostic plots
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

# immunity
plot_ib <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$IB})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("IB") +
  ggtitle("Pre-erythrocytic immunity",subtitle = "reduces the probability of infection following an infectious challenge")

plot_id <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$ID})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("ID") +
  ggtitle("Blood-stage immunity",subtitle = "reduces the probability of detection and reduces infectiousness to mosquitoes")

plot_icm <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$ICM})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("ICM") +
  ggtitle("Maternal clinical immunity",subtitle = "reduces the probability of clinical disease, acquired maternally")

plot_ica <- ggplot(data=data.frame(x=sapply(RACD_init$humans,function(x){x$ICA})))+
  geom_histogram(aes(x=x)) +
  theme_bw() +
  xlab("ICA") +
  ggtitle("Acquired clinical immunity", subtitle = "reduces the probability of clinical disease, acquired from previous exposure")

grid.arrange(plot_ib,plot_id,plot_icm,plot_ica)

# compile the simulation
sourceCpp(here("intervention-src/main.cpp"))

# run the simulation
RACD_out <- tiny_racd(humans_param = RACD_init$humans,
                      house_param = RACD_init$houses,
                      mosy_param = RACD_init$mosy,
                      theta = RACD_theta,
                      tmax = 365*20)

library(reshape2)

plot_outS <- ggplot(data = melt(RACD_out$state,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.75) +
  theme_bw()

plot_outA <- ggplot(data = melt(RACD_out$age,id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

# ggplot(data = melt(RACD_out$clinical_incidence,id.vars="time")) +
#   geom_line(aes(x=time,y=value,color=variable)) +
#   theme_bw()

ggplot(data = melt(RACD_out$trans[,-2],id.vars = c("time","EIR_mean","EIR_var"))) +
  geom_ribbon(aes(x=time,ymin=pmax(EIR_mean - sqrt(EIR_var),0),ymax=EIR_mean + sqrt(EIR_var)),alpha=0.85) +
  geom_line(aes(x=time,y=EIR_mean)) +
  xlab("Time (days)") + ylab("EIR (entomological inoculation rate)") +
  theme_bw()

plot_outFOI <- ggplot(data = melt(RACD_out$trans[,1:2],id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

plot_outM <- ggplot(data = melt(RACD_out$mosy[,!names(RACD_out$mosy) == "S"],id.vars="time")) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.65) +
  theme_bw()

grid.arrange(plot_outS,plot_outA,plot_outFOI,plot_outM)
