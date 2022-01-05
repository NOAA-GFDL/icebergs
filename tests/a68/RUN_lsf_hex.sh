#!/bin/bash

#ocean drag values
ss=(10000)
gcoef=(1.e20)
nstress=(68000)
tstress=(50000)
damp=(0.1)
name=(zh_lsf_10kss_el_gl1e20_t50k_n68k_od7_gt_realloc_dc.1.nc)
#name=(zhex_lsf_el_g1e60_t50k_n400k_od7_2kss.nc)
radius=6363827
#6.371229e6

for ((i = 0; i < ${#ss[@]}; ++i)); do
    substeps=(${ss[$i]})
    gc=(${gcoef[$i]})
    fname2=(${name[$i]})
    ns=(${nstress[$i]})
    ts=(${tstress[$i]})
    dc=(${damp[$i]})

    echo $substeps
    echo $gc
    echo $fname2
    echo $ns
    echo $ts

    sed  "s/<ss>/$substeps/; s/<gc>/$gc/; s/<dc>/$dc/; s/<rad>/$radius/; s/<ns>/$ns/; s/<ts>/$ts/; s/<fname>/$fname2/g" \
	 hex_r_lsf.nml > input.nml
    rm $fname2
    ./RUN1p
    # &>/dev/null &
done
