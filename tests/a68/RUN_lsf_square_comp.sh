#!/bin/bash

#ocean drag values
ss=(200)
gcoef=(1.e8)
odrag=(15)
nstress=(8000)
tstress=(8000)
A68dy=0.
A68dx=0.
name=(sq_lsf_gc1e8_od15_c8k.nc)
Eradius=6.378e6
#6363827

for ((i = 0; i < ${#ss[@]}; ++i)); do
    substeps=(${ss[$i]})
    gc=(${gcoef[$i]})
    fname2=(${name[$i]})
    #bname2=(${bond_name[$i]})
    ns=(${nstress[$i]})
    ts=(${tstress[$i]})
    od=(${odrag[$i]})

    echo $substeps
    echo $gc
    echo $fname2
    echo $ns
    echo $ts
    echo $od

    sed  "s/<ss>/$substeps/;  s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<gc>/$gc/; s/<rad>/$Eradius/; s/<od>/$od/; s/<ns>/$ns/; s/<ts>/$ts/; s/<name>/$fname2/g" \
	 square_r_lsf_comp.nml > input.nml
    rm $fname2
    #rm $bname2
    ./RUN1p
done
