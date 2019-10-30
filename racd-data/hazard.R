rm(list=ls());gc()
library(readstata13)
library(here)
library(readr)
library(tidyverse)


house_data <- read.dta13(here("prev/household_appended.dta"))
index_data <- read.dta13(here("prev/index_cases_all_appended_ea_and_intcode_joined_clean_August2018_2.dta"))
inc_data <- read.dta13(here("prev/locality_incidence_cases.dta"))
census_data <- read.dta13(here("prev/swazi census data.dta"))
table1_data <- read.dta13(here("prev/table1_dataset_final_Apr11.dta"))
weights_data <- read.dta13(here("prev/weights2.dta"))
xsect_data <- read.dta13(here("prev/crosssect.dta"))

caseInc <- readr::read_csv(here("prev/weeklyincidence.csv"))

caseIncM <- reshape2::melt(caseInc[,1:53],id.vars="ea_no")
caseIncM$variable <- as.integer(caseIncM$variable)

ggplot(data=caseIncM) +
  geom_line(aes(x=variable,y=value,group=ea_no,color=as.factor(ea_no))) +
  geom_smooth(aes(x=variable,y=value)) +
  guides(color=FALSE) +
  theme_bw()

# attack rate
AR <- (caseInc$allweeks/caseInc$ea_pop_update)
FOI <- -log(1-AR)/365

# calculate daily EIR with standard b parameter
b <- 0.5
EIR <- FOI/b

EIR_dens <- density(EIR,kernel = "gaussian",from=0)
EIR_draws <- rnorm(1e6, sample(EIR, size = 1e6, replace = TRUE), EIR_dens$bw)
EIR_draws <- EIR_draws[EIR_draws>0]
