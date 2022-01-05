#!/bin/bash

#ocean drag values
ss=(200 300 400 750 500 100 1000)
name=(sstest2_200.nc sstest2_300.nc sstest2_400.nc sstest2_750.nc sstest2_500.nc sstest2_100.nc sstest2_1000.nc)

for ((i = 0; i < ${#ss[@]}; ++i)); do
    substeps=(${ss[$i]})
    fname2=(${name[$i]})

    echo $substeps
    echo $fname2

    sed  "s/<ss>/$substeps/; s/<fname>/$fname2/g" \
	 check_ss2.nml > input.nml
    ./RUN1p
    # &>/dev/null &
done
