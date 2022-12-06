# CautiousLearning

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CautiousLearning

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

To run the simulations and produce the plots, a `bash` script that performs all the necessary computations can be launched by
```
./scripts/simulate-all.sh
```

The script requires the choice of the number of cores to run the process in parallel.
For instance, by setting `<number-of-cores>` to be 1 the computation is not run in parallel.
If the number of cores are not a perfect divisor of the number of simulations, the script automatically checks for exceeding simulations and removes them.

**simulate-all.sh**
```
totsim=200
ncores=<number-of-cores>
add_existing=false
echo "ncores: $ncores"
echo "add_existing: $add_existing"
neach=$((($totsim+$ncores-1)/$ncores))
echo "neach: $neach"

parallel --ungroup julia simulate/run-GICP.jl -n $neach -t {1} -a $add_existing ::: {1..<number-of-cores>}
wait

parallel --ungroup julia simulate/run-ARL.jl -n $neach -t {1} -a $add_existing ::: {1..<number-of-cores>}
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
