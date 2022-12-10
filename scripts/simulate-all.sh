totsim=200
ncores=16
remove_excessive=true
add_existing=false

echo "ncores: $ncores"
echo "add_existing: $add_existing"
neach=$((($totsim+$ncores-1)/$ncores))
echo "neach: $neach"

parallel --ungroup julia simulate/run-GICP.jl -n $neach -t {1} -a $add_existing ::: {1..16}
wait

parallel --ungroup julia simulate/run-ARL.jl -n $neach -t {1} -a $add_existing ::: {1..16}
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
