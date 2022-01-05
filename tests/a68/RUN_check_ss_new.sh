#!/bin/bash

#ocean drag values
ss=(200)
name=(sstest_200.nc)

for ((i = 0; i < ${#ss[@]}; ++i)); do
    substeps=(${ss[$i]})
    fname2=(${name[$i]})

    echo $substeps
    echo $fname2

    sed  "s/<ss>/$substeps/; s/<fname>/$fname2/g" \
	 check_ss.nml > input.nml
    ./RUN1p
    # &>/dev/null &
done
