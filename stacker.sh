ffmpeg -i Simulation-11-simulation.mp4 -i Simulation-8-simulation.mp4 -i Simulation-9-simulation.mp4 \
-filter_complex "vstack=inputs=3" \
several-simulations.mp4

ffmpeg -i mtGFP-1-combined.mp4 -i mtGFP-2-combined.mp4 \
-filter_complex "vstack=inputs=2" \
several-experiments.mp4

