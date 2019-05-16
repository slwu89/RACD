# tiny-racd: a small, portable version of the RACD model

## tiny-human:
  * tiny-human.cpp: functions to operate on humans (as a list)
  * tiny-human.R: run the RACD (human) model with the C++ functions (just to check it is working)

## tiny-population:
  * tiny-population.cpp: run the entire RACD human simulation in C++
  
## racd-src
contains the source code of the full RACD dynamic model (humans + mosquitos) in C++

## setup files
  * landscape.R: make a random landscape, calculate sigma for aquatic habitats and psi for dwellings
  * racd-setup.R: using a random landscape, set up the model
  * racd-setup.cpp: helper functions (immunity initialization) for above file

# scratchpad

Obvious code optimizations:
  * fewer loops over people when updating C,W,Z global variables
  * only draw from RNG once to choose when each person dies (at birth)
