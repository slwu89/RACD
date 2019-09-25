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

ggplot(data=caseIncM) +
  geom_line(aes(x=variable,y=value,group=ea_no,color=as.factor(ea_no))) +
  guides(color=FALSE)

# attack rate
AR <- (caseInc$allweeks/caseInc$ea_pop_update)
FOI <- -log(1-AR)/365

# calculate daily EIR with standard b parameter
b <- 0.55
EIR <- FOI/b
