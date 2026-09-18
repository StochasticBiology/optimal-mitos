# pull existing experimental data and analysis script 
git clone https://github.com/StochasticBiology/plant-mito-dynamics

# run analysis script on existing experimental data
cd plant-mito-dynamics
find . -name "*[0-9].xml" -print0 | while read -d $'\0' file
do  
    Rscript trajectory-analysis.R "$file" 1.6 10 50 100
done

# put analysis script in folder with new data
cp trajectory-analysis.R ../cipro-rawtrajectories

# run analysis script on new data
cd ../cipro-rawtrajectories
find . -name "*[0-9].xml" -print0 | while read -d $'\0' file
do  
    Rscript trajectory-analysis.R "$file" 1.6 10 50 100
done

# run simulations
cd ..
gcc -o3 simulate.c -lm -o simulate.ce
# given parameter sets
./simulate.ce 0 > tmp0 &
# parameter scans (take ~hours)
./simulate.ce 1 > tmp1 &
./simulate.ce 2 > tmp2 &
./simulate.ce 3 > tmp3 &
./simulate.ce 4 > tmp4 &
./simulate.ce 5 > tmp5 &

# run analysis and visualisation code
Rscript plot-simulate-cipro.R

# generate example videos (takes ~dozens of minutes)
Rscript generate-video.R
