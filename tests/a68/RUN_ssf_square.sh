#!/bin/bash

#ocean drag values
ss=200
gcoef=(1.e4 5.e3 1.e1)
#gcoef=(1.e4 5.e3 1.e1)
rdir=S_17.8_ydn.1_r
odrag=17.8
#19.5
#21.25
#15
nstress=(24 23 23)
#25
#(25 30 35 30 25)
tstress=(6.5 6.5 6.5)
#tstress=(7 7 7)
#7
#(6 5 5 4 4)
#nstress=(50000)
#tstress=(3000)
a68=.true.
#.false

# A68dx=-0.01872
# A68dy=-0.09338

#A68dx=-37.2954
#A68dy=-55.2975

A68dx=-37.5285
A68dy=-55.2166


#name=(sq_ssf_gc1e6_od15_n20k_t5k.nc)
fname2=loop
#Eradius=6.378e6
Eradius=6363827
pg=.false.
#.true.
#6363827

for ((i = 0; i < ${#nstress[@]}; ++i)); do
    #substeps=(${ss[$i]})
    gc=(${gcoef[$i]})
    #fname2=(${name[$i]})
    #bname2=(${bond_name[$i]})
    ns=(${nstress[$i]})
    ts=(${tstress[$i]})
    #od=(${odrag[$i]})

    #echo $substeps
    echo $gc
    #echo $fname2
    echo $a68
    echo $ns
    echo $ts
    echo $pg
    echo $rdir
    #echo $od

    sed  "s/<ss>/$ss/; s/<a68>/$a68/; s/<pg>/$pg/; s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<rdir>/$rdir/; s/<gc>/$gc/; s/<rad>/$Eradius/; s/<od>/$odrag/; s/<ns>/$ns/; s/<ts>/$ts/; s/<name>/$fname2/g" \
	 square_r_ssf.nml > input.nml
    #rm $fname2
    #rm $bname2
    ./RUN1p
done
