#!/bin/bash

Rearth=6363827
odrag=(17.8)
name=(od_17.8_tn_ng_ydn.1.nc)

for ((i = 0; i < ${#odrag[@]}; ++i)); do
    drag=(${odrag[$i]})
    fname=(${name[$i]})

    echo $drag
    echo $fname
    echo $Rearth

    sed  "s/<od>/$drag/; s/<re>/$Rearth/; s/<name>/$fname/g" \
	 input_fill2.nml > input.nml
    ./RUN1p
    # &>/dev/null &
done
