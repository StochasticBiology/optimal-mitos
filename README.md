# optimal-mitos

Investigation of the collective behaviour of plant mitochondria, through the lens of multi-objective optimisation

<img width="2692" height="1566" alt="image" src="https://github.com/user-attachments/assets/48f0984f-2651-449d-b542-fa5158e05a1c" />

Overview
----

This repo contains agent-based simulation code, data analysis code, and summaries of new experiments. This new data is analysed alongside existing data from https://github.com/StochasticBiology/plant-mito-dynamics . `pipeline.sh` should wrap this whole process.

Details
----
`simulate.c` runs quasi-physical simulation of mitochondria in the cell, constructing trajectories and adjacency matrices. It does a parameter sweep and several specific parameterisation cases.

`plot-simulate-cipro.R` plots the output from these simulations and also a comparison with bio data. For this comparison with experiments, pull https://github.com/StochasticBiology/plant-mito-dynamics into `./plant-mito-dynamics-main/` and run (the first part of) `wrapper.sh` from there. This analyses the trajectories in that repo and summarises the adjacency matrices and coordinates for inclusion here.

`cipro-rawtrajectories/` contains XML files describing mitochondrial trajectories for the new ciprofloxacin experiments. The first part of the `wrapper.sh` script from the repo above can also be used to summarise these dynamics.

`generate-video.R` produces a visualisation of the social network construction process. This also calls the Bash scripts `concatenate-video.sh` and `stacker.sh` which use `ffmpeg` to process and compile the videos.
