#!/bin/bash

#ocean drag values
ss=(200 300 400 750 500 100 1000)
name=(sstest_200.nc sstest_300.nc sstest_400.nc sstest_750.nc sstest_500.nc sstest_100.nc sstest_1000.nc)

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
