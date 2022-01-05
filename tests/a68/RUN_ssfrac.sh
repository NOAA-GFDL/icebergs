#!/bin/bash

#ocean drag values
ss=(200)
# 400 500 1000)
gcoef=(1.e6)
nstress=(9000)
tstress=(2500)
name=(lsf_nsf_200_gc1e5_t2.5k_n9k.nc)
# sstest_400_ssf.nc sstest_500_ssf.nc sstest_1000_ssf.nc)

for ((i = 0; i < ${#ss[@]}; ++i)); do
    substeps=(${ss[$i]})
    gc=(${gcoef[$i]})
    fname2=(${name[$i]})
    ns=(${nstress[$i]})
    ts=(${tstress[$i]})

    echo $substeps
    echo $gc
    echo $fname2
    echo $ns
    echo $ts

    sed  "s/<ss>/$substeps/; s/<gc>/$gc/; s/<ns>/$ns/; s/<ts>/$ts/; s/<fname>/$fname2/g" \
	 check_ss3.nml > input.nml
    rm $fname2
    ./RUN1p
    # &>/dev/null &
done
