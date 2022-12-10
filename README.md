# CautiousLearning

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CautiousLearning

To (locally) reproduce this project, do the following:

1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths for the `Julia` scripts.

2. Run the required scripts for reproducing the part of the paper you are interested in. Below is a description of the relevant script and their functionality.

### scripts/admissionsICU.jl

This `julia` script reproduces the application to the ICU data in Section 5 of the main paper.
To reproduce the figures and table, run the script as
```
   julia scripts/admissionsICU.jl
```
The script will produce the required figures and table in the `plots/ICUadmissions` folder.

### scripts/plot-path-comparison-cl.jl

This `julia` script reproduces Figure 1 in the main paper, with each plot in a separate PDF file.
To reproduce the Figure, run the script as
```
   julia scripts/plot-path.comparison.jl
```
The script will produce the required plots `plots/sims/window-of-opportuniy` folder.

### scripts/simulate-all.sh
The script reproduces all figures and tables in the simulation study of Section  4 of the main paper, as well as the additional simulations contained in the supplemental material, in separate PDF and TeX files.

To run the script, make sure all relevant scripts are executable by running
```
   chmod +x scripts/simulate-all.sh
   chmod +x scripts/collect-and-plot.sh
   chmod +x scripts/trim.sh
```
and then run the full simulation study using `./scripts/simulate-all.sh`.

The script requires the choice of the number of cores to run the process in parallel.
For instance, by setting all occurrences of `<ncores>` in the script to 1 the computation is not run in parallel.
The default number of cores is the one used in our local machine, 16.
If the number of cores are not a perfect divisor of the number of simulations, the script automatically checks for exceeding simulations and removes them.

```
totsim=200
ncores=<ncores>
add_existing=false
echo "ncores: $ncores"
echo "add_existing: $add_existing"
neach=$((($totsim+$ncores-1)/$ncores))
echo "neach: $neach"

parallel --ungroup julia simulate/run-GICP.jl -n $neach -t {1} -a $add_existing ::: {1..<ncores>}
wait

parallel --ungroup julia simulate/run-ARL.jl -n $neach -t {1} -a $add_existing ::: {1..<ncores>}
wait

for ((i=1; i<$ncores; i++))
do
    nsim=$(($totsim+$i))
    echo $nsim
    simstring="simulation=$nsim"
    echo $simstring
    find ../data/sims/ -name "Simula*_$simstring*" -exec rm '{}' \;
done

./collect-and-plot.sh
```
**Note**: to produce plots and tables using the `scripts/collect-and-plot.sh` script, as part of the `scripts/simulate-all.sh` procedure, please change to your local directory path at the first line of the `scripts/summary/plots-tables.R` script.

All figures and tables are saved in subdirectories of the `plots/sims` directory.
