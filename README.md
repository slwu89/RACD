# RACD

## to-do
  * make sure people aren't getting ITNs when we do IRS
  * RfVC happens once a season; need to have a lag window where we don't re-spray immediately
  * need to add importation at some rate
  * need to add more failure probabilities for the main interventions



<!-- # RACD
stochastic simulation of reactive case detection for p. falciparum

* RACD is an R package that implements the individual based model in C++11/14 (requires C++11/14 compatible compiler). It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACD")`

We are currently working to extend the RACD package to include mosquito dynamics (a standalone mosquito simulator can be found in Auxiliary folder).

<!-- ## RACD Work Queue (from JM)
1. Calculate mosquito population size required to generate required EIR - DONE
2. Link vector model to human model in malaria IBM (see MalariaODEVectorsOnlyMay30.R, reference https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153) - IN PROGRESS
3. Incorporate targeted interventions (all 3 interventions: RACD, etc.) - NEED TO HEAR FROM MICHELLE
4. Visualize targeted interventions - NOT YET
5. Incorporate possibility for seasonality in mosquito numbers - NOT YET (easy though, just implement forcing on K)
6. Validate human model - HOW TO VALIDATE?
7. Make sure that functions all make sense - NOTHING YET
8. Validate modified treatment node in human model - NOT SURE WHAT THIS MEANS -->

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

  * eggOV: Number of eggs per oviposition per mosquito (w/0 interventions)

## To-do

  - [x] all simulation pars need to be in a single hash-table in the village class
  - [x] if its possible to have empty houses, we need to renormalize the Psi vector to give them probabiltiy 0, since all bites are conditioned on the bite going at least _somewhere_
  - [x] in main sim loop, don't simulate empty houses
  - [ ] develop a central "intervention manager" which houses have pointer to (so humans can tell the manager if they got sick that day, for example). This object at the end of the day can fire up interventions (RACD, MSAT, Rf VC)
  - [ ] very low priority: compartment/infectiousness fn's could be `static`, and also use lambdas for callbacks.

## Classes

### Human -->
