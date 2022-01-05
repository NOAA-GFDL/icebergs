#!/bin/bash

#ocean drag values
ss=200
gcoef=(1.e4)
rdir=S_17.8_ydn.1_r
odrag=17.8
nstress=(23)
tstress=(10)
a68=.true.

A68dx=-37.5285
A68dy=-55.2166

fname2=sq_ssf_ssg1.e4_od17.8_addsm.true._rectx-37.5285_y-55.2166_n23_t10.nc
Eradius=6363827
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

    sed  "s/<ss>/$ss/; s/<a68>/$a68/; s/<pg>/$pg/; s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<rdir>/$rdir/; s/<gc>/$gc/; s/<rad>/$Eradius/; s/<od>/$odrag/; s/<ns>/$ns/; s/<ts>/$ts/; s/<name>/$fname2/g" \
	 square_r_ssf.nml > input.nml
    ./RUN1p
done
