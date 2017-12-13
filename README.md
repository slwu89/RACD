# RACD
stochastic simulation of reactive case detection for p. falciparum

* RACDaux is a directory containing an R package that includes C code for within-host immune dynamics compatible with `deSolve` ODE solvers as well as auxiliary setup and graphics R functions. It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACDaux")`
* RACD is the main R package that implements the individual based model in C++11/14 (requires C++11/14 compatible compiler). It can be installed from github with `devtools::install_github(repo = "slwu89/RACD",subdir = "RACD")`
