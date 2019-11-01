# Simulation Design

The Monte Carlo experiments will need to sweep over a few dimensions:

  * intervention type: RACD, RfMDA, RfMDA + RfVC, RfVC, control
  * landscape type: clustered, CSR, hardcore
  * parameters: EIR, importation, p_index, p_neighbor

There will be some assumptions baked in, but most of them will be related to model paramters, which were fitted
by the Imperial College/MRC modelling centre.

Once we've sampled locations of dwellings and habitats, we need to calculate the sigma value (SD) for each aquatic habitat.
Following JM, it's going to be the distance from the habitat to the nearest dwelling; this risk dispersion is an assumption, but
generally encodes the idea that in regions where habitats are dense, competition between mosquitoes means that each habitat does not produce many
mosquitoes; therefore dispersion of bites is more limited. Productive habitats in low density areas can produce mosquito populations able
to spread bites to a large geographic region.

The radius for spatial targeting is 500m. There are about 4 people per house, and about 100 hourses per EA. Standardize the bounding box to 2km.
The consecutive lag for IRS is 60 days.
