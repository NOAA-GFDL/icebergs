#!/bin/bash

#grounding coef
gcoef=(1.e4)
#number of substeps
ss=90
#restart directory
rdir=S_17.8_ydn.1_ss90_r
#scale factor for ocean drag coefficients
odrag=17.8
#normal and shear stress thresholds (kPa)
nstress=(18)
tstress=(10)
#Use A68a grounding zone?
a68=.true.
#lower-left coords of grounding zone
A68dx=-37.51
A68dy=-55.2166
#Earth radius
Eradius=6363827
#Use power-law grounding drag? If false, uses linear grounding drag.
pg=.false.

for ((i = 0; i < ${#nstress[@]}; ++i)); do
    gc=(${gcoef[$i]})
    ns=(${nstress[$i]})
    ts=(${tstress[$i]})

    echo $gc
    echo $a68
    echo $ns
    echo $ts
    echo $pg
    echo $rdir
    echo $ss
    echo $A68dx
    echo $A68dy

    sed  "s/<ss>/$ss/; s/<a68>/$a68/; s/<pg>/$pg/; s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<rdir>/$rdir/; s/<gc>/$gc/; s/<rad>/$Eradius/; s/<od>/$odrag/; s/<ns>/$ns/; s/<ts>/$ts/g" \
	 square_r_ssf.nml > input.nml
    ./RUN1p
done

##Another good run, with more sub-steps:
#gcoef=(1.e4)
#ss=200
#rdir=S_17.8_ydn.1_r
#A68dx=-37.5285
#A68dy=-55.2166
#Eradius=6363827
#pg=.false.
