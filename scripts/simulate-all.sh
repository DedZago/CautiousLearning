totsim=200
# ncores=16
ncores=$(grep -c ^processor /proc/cpuinfo)
echo "ncores: $ncores"

neach=$((($totsim+$ncores-1)/$ncores))
echo "neach: $neach"

parallel --ungroup julia simulate/run-GICP.jl -n $neach -t {1} -a false ::: {1..$ncores}
wait

parallel --ungroup julia simulate/run-ARL.jl -n $neach -t {1} -a false ::: {1..$ncores}
wait

./collect-and-plot.sh
