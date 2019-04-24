# RACD
stochastic simulation of reactive case detection for p. falciparum

* RACD is an R package that implements the individual based model in C++11/14 (requires C++11/14 compatible compiler). It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACD")`

We are currently working to extend the RACD package to include mosquito dynamics (a standalone mosquito simulator can be found in Auxiliary folder).

## RACD Work Queue (from JM)
1. Calculate mosquito population size required to generate required EIR - DONE
2. Link vector model to human model in malaria IBM (see MalariaODEVectorsOnlyMay30.R, reference https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153) - IN PROGRESS
3. Incorporate targeted interventions (all 3 interventions: RACD, etc.) - NEED TO HEAR FROM MICHELLE
4. Visualize targeted interventions - NOT YET
5. Incorporate possibility for seasonality in mosquito numbers - NOT YET (easy though, just implement forcing on K)
6. Validate human model - HOW TO VALIDATE?
7. Make sure that functions all make sense - NOTHING YET
8. Validate modified treatment node in human model - NOT SURE WHAT THIS MEANS

## RACD References
* Human Model:
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923296/
  * https://www.ncbi.nlm.nih.gov/pubmed/20711482/
* Vector Model:
  * https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153

## RACD Parameters

New and changed parameters from initial version of model

### Main simulation

  * IRSduration: IRS works until it doesn't (has eff=1 until an exp RV w/lambda=1/IRSduration turns it off)
  * ITNduration: same

### Entomological

## Possible Optimizations

for humans: compartment_funs, infectiousness_funs should be static class members initialized once for all people


## to-do

  * all simulation pars need to be in a single hash-table in the village class
