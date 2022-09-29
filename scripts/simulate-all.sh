totsim=600
ncores=16
add_existing=true
# ncores=$(grep -c ^processor /proc/cpuinfo)
echo "ncores: $ncores"
echo "add_existing: $add_existing"
neach=$((($totsim+$ncores-1)/$ncores))
echo "neach: $neach"

parallel --ungroup julia simulate/run-GICP.jl -n $neach -t {1} -a $add_existing ::: {1..16}
wait

parallel --ungroup julia simulate/run-ARL.jl -n $neach -t {1} -a $add_existing ::: {1..16}
wait

./collect-and-plot.sh
