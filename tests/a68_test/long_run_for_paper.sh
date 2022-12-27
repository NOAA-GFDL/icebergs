#!/bin/bash

Rearth=6363827
odrag=(17.8)
name=(longrun_od_17.8_tn_ng_ydn.1_ss90_ns18_rx-37.51_ry-55.2166)
ss=90
gcoef=(1.e4)
nstress=(18)
A68dx=-37.51
A68dy=-55.2166

for ((i = 0; i < ${#odrag[@]}; ++i)); do
    drag=(${odrag[$i]})
    fname=(${name[$i]})
    gc=(${gcoef[$i]})
    ns=(${nstress[$i]})

    echo $drag
    echo $fname
    echo $Rearth
    echo $ss
    echo $A68dx
    echo $A68dy
    echo $ns

    sed  "s/<od>/$drag/; s/<re>/$Rearth/; s/<ss>/$ss/; s/<ns>/$ns/; s/<gc>/$gc/; s/<xd>/$A68dx/; s/<yd>/$A68dy/; s/<name>/$fname/g" \
	 long_run.nml > input.nml
    ./RUN1p
    # &>/dev/null &
done
