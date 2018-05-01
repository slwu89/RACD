# RACD
stochastic simulation of reactive case detection for p. falciparum

* RACDaux is a directory containing an R package that includes C code for within-host immune dynamics compatible with `deSolve` ODE solvers as well as auxiliary setup and graphics R functions. It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACDaux")`
* RACD is the main R package that implements the individual based model in C++11/14 (requires C++11/14 compatible compiler). It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACD")`

## RACD Work Queue (from JM)
1. Calculate mosquito population size required to generate required EIR
2. Link vector model to human model in malaria IBM (see MalariaODEVectorsOnlyMay30.R, reference https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)
3. Incorporate targeted interventions (all 3 interventions: RACD, etc.)
4. Visualize targeted interventions
5. Incorporate possibility for seasonality in mosquito numbers
6. Validate human model
7. Make sure that functions all make sense
8. Validate modified treatment node in human model

## Stuff to investigate:
* Lifespan follows geometric distribution, will be right as the ensemble average but very very wrong for individuals and small samples. Should probably allow for discrete versions of Weibull and Gamma age distributions.
* "superinfection" is just going back to T or D from other compartments, barely more sophisticated than the RM model. Changing this would be a substantial undertaking.
* look into using Gibbs processes to sample houses and breeding sites (can use `spatstat`), maybe also bivariate marked point processes for between type covariance.
  * Cox & Cluster processes:
    * Crucial assumptions: offspring mutually independent within a cluster (only depend on position of parents), parents are independent thus implying clusters are independent (can only depend on first-order intensity terms; in general will follow a Poisson or Cox process).
    * Remember: all Neyman-Scott cluster processes are just Cox processes with different modulating intensity. In fact they all belong to the class of "shot-noise field" models because they are Cox processes with driving intensity Lambda equal to the superposition of the offspring kernels centered on parent points. The intensity surface is a random field because the distribution of the parent generator points is random.
    * Matérn cluster process when represented as a Cox process has driving intensity: $\Lambda(u) = \sum_{i} h(u-y_{i})$, where $h(z) = \frac{\mu}{\pi R^{2}}$ if $\left \| z \right \|<R$, $0$ otherwise.
    * See similar results for Thomas process (isotropic Gaussian density of offspring points), Cauchy, variance-Gamma.
    * For non Neyman-Scott clustering process see `rGaussPoisson` for Gauss-Poisson clustering process (violates following common assumptions: offspring independent within a cluster, Poisson number of offspring, isotropic clusters).
  * Gibbs processes:
    * Pairwise interaction:
      * Hard-core process (inhibitory):
      * Strauss process (inhibitory):
      * Strauss/Hard-core process (inhibitory):

## RACD References
* Human Model:
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923296/
  * https://www.ncbi.nlm.nih.gov/pubmed/20711482/
* Vector Model:
  * https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153
