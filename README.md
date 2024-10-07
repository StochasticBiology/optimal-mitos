# optimal-mitos

`simulate.c` runs quasi-physical simulation of mitochondria in the cell, constructing trajectories and adjacency matrices. It does a parameter sweep and several specific parameterisation cases.

`plot-simulate.R` plots the output from these simulations and also a comparison with bio data. For this comparison with experiments, pull https://github.com/StochasticBiology/plant-mito-dynamics into `./plant-mito-dynamics-main/` and run `wrapper.sh` from there. This analyses the trajectories in that repo and summarises the adjacency matrices and coordinates for inclusion here.

`generate-video.R` produces a visualisation of the social network construction process.
