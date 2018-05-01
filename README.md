# RACD
stochastic simulation of reactive case detection for p. falciparum

* RACDaux is a directory containing an R package that includes C code for within-host immune dynamics compatible with `deSolve` ODE solvers as well as auxiliary setup and graphics R functions. It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACDaux")`
* RACD is the main R package that implements the individual based model in C++11/14 (requires C++11/14 compatible compiler). It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACD")`

## RACD Work Queue
1. Calculate mosquito population size required to generate required EIR
2. Link vector model to human model in malaria IBM (see MalariaODEVectorsOnlyMay30.R, reference https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)
3. Incorporate targeted interventions (all 3 interventions: RACD, etc.)
4. Visualize targeted interventions
5. Incorporate possibility for seasonality in mosquito numbers
6. Validate human model
7. Make sure that functions all make sense
8. Validate modified treatment node in human model

## Assumptions to check and change:
* Lifespan follows geometric distribution, will be right as the ensemble average but very very wrong for individuals and small samples. Should probably allow for discrete versions of Weibull and Gamma age distributions.
* "superinfection" is just going back to T or D from other compartments, barely more sophisticated than the RM model. Changing this would be a substantial undertaking.

## RACD References
* Human Model:
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923296/
  * https://www.ncbi.nlm.nih.gov/pubmed/20711482/
* Vector Model:
  * https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153
