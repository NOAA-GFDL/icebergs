#!/bin/bash

#ocean drag values
ss=(200)
#(100 300 750)
gcoef=(1.e6)
nstress=(20000)
tstress=(8000)
name=(ssf_nstressfix_200_gc1e6_t8k_n20k.nc)
#(sstest_100_ssf.nc sstest_300_ssf.nc sstest_750_ssf.nc)

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

    sed  "s/<ss>/$substeps/;  s/<gc>/$gc/; s/<ns>/$ns/; s/<ts>/$ts/; s/<fname>/$fname2/g" \
	 check_ss3_ssf.nml > input.nml
    rm $fname2
    ./RUN1p
    # &>/dev/null &
done
