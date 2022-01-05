#!/bin/bash

# odrag=(7.2)
# nstress=(82000)
# tstress=(82000)
# A68dy=0.1
# A68dx=-0.08
# name=(zh_ssf_el_gl1e30_od7.2_compd82.k_tripy.1xn.08_cspr1en7_10e3cd_powgnd.nc)
# radius=6363827

#ocean drag values
ss=(200)
gcoef=(1.e4)
pg=.false.
#.true.
odrag=(6)
#5
nstress=105
#(85)
tstress=55
#(55)
a68=.true.
#.false.
ssg=.true.
#.true.
A68dy=0.06792000000000087
#-0.02 #for 5.4
#0.06
A68dx=-0.27071999999997587
#-0.323
#-0.15 #for 5.4
#-0.2

#name=(zzh_ssf_yn.075_gl1e10_od6_nt80k_a68_y.06xn.2.nc)
#name=(zzh_ssf_yn.1_gl1e10_od5.4_nt82k_a68_yn.02xn.15.nc)
#name=(hz_ssf_yn.1_gl1e6_od5_n90_t60_a68.nc)
radius=6363827

for ((i = 0; i < ${#ss[@]}; ++i)); do
    substeps=(${ss[$i]})
    gc=(${gcoef[$i]})
    #fname2=(${name[$i]})
    #bname2=(${bond_name[$i]})
    ns=(${nstress[$i]})
    ts=(${tstress[$i]})
    od=(${odrag[$i]})

    echo $substeps
    echo $gc
    #echo $fname2
    echo $ns
    echo $ts
    echo $od
    echo $A68dx
    echo $A68dy
    echo $pg
    echo $ssg

    sed  "s/<ss>/$substeps/; s/<pg>/$pg/; s/<a68>/$a68/; s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<gc>/$gc/; s/<rad>/$radius/; s/<ssg>/$ssg/; s/<od>/$od/; s/<ns>/$ns/; s/<ts>/$ts/g" \
	 hex_r_ssf.nml > input.nml
    rm $fname2
    rm $bname2
    ./RUN1p
done
