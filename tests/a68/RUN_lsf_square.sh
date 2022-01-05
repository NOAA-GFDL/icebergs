#!/bin/bash

#ocean drag values
ss=(200)
gcoef=(1.e6)
odrag=(15)
nstress=(18)
tstress=(2.75)
# A68dy=0.
# A68dx=0.
a68=.false.
#.true.
A68dx=0.
#-0.4017
A68dy=0.
#-0.05
dxw=0.
#n.4017
dyw=0.
#n.05
#name=(sq_lsf_gc1e6_od15_n18k_t2.75k_cd4e3.nc)
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
    echo $a68

    sed  "s/<ss>/$substeps/;  s/<a68>/$a68/; s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<xdw>/$dxw/; s/<ydw>/$dyw/; s/<gc>/$gc/; s/<rad>/$Eradius/; s/<od>/$od/; s/<ns>/$ns/; s/<ts>/$ts/; s/<name>/$fname2/g" \
	 square_r_lsf.nml > input.nml
    rm $fname2
    #rm $bname2
    ./RUN1p
done
