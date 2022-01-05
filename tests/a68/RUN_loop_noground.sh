#!/bin/bash

#ocean drag values
odrag=(7)
#(12 15 20 1)
name=(test_5e6_od7_hex_check.nc)
#(h_yn.05_ss100e1e6_nf_od12.nc h_yn.05ss100e1e6_nf_od15.nc h_yn.05ss100e1e6_nf_od20.nc h_yn.05ss100e1e6_nf_od1.nc)
youngs=5e6
#5e6
#odrag=(7.2)
#name=(hex_ydn.05_wgs84_gc0_nofrac_od7.2.nc)

radius=6363827

echo $radius

for ((i = 0; i < ${#odrag[@]}; ++i)); do
    drag=(${odrag[$i]})
    fname=(${name[$i]})

    echo $drag
    echo $fname

    sed  "s/<od>/$drag/; s/<rad>/$radius/; s/<sc>/$youngs/; s/<name>/$fname/g" \
	 input_fill_hex_noground.nml > input.nml
    ./RUN1p
    # &>/dev/null &
done
