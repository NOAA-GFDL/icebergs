#!/bin/bash

#ocean drag values
#odrag=(6.5 6 5 7 4 8 3 1)
#name=(hex_wgs84_nofrac_od6.5.nc hex_wgs84_nofrac_od6.nc hex_wgs84_nofrac_od5.nc hex_wgs84_nofrac_od7.nc hex_wgs84_nofrac_od4.nc hex_wgs84_nofrac_od8.nc hex_wgs84_nofrac_od3.nc hex_wgs84_nofrac_od1.nc)

#odrag=(6 6.5 7)
#name=(hex_ydn.5_wgs84_nofrac_od6.nc hex_ydn.5_wgs84_nofrac_od6.5.nc hex_ydn.5_wgs84_nofrac_od7.nc)

#this was to test having rev_mind=.true. and rev_mind=.false., to determine where it places the seamount.
#odrag=(6.)
#name=(test_rm2_false.nc)

#odrag=(6.075 6.025)
#name=(hex_yn.075_nf_od6.075.nc hex_yn.075_nf_od6.025.nc)

#odrag=(6. 6.5 5.5)
#name=(hex_yn.1_nf_od6.nc hex_yn.1_nf_od6.5.nc hex_yn.1_nf_od5.5.nc)

#odrag=(4.75 4.85)
#name=(hex_yn.125_nf_od4.75.nc hex_yn.125_nf_od4.85.nc)

#odrag=(2.9 2.95)
#name=(hex_yn.175_nf_od2.9.nc hex_yn.175_nf_od2.95.nc)

odrag=(9.5 10.5)
name=(hex_yn.1_nf_od9.5_rad_drag_E1e6.nc hex_yn.1_nf_od10.5_rad_drag_E1e6.nc)


#odrag=(5.5 5)
#name=(hex_yn.1_nf_od5.5_radius_drag.nc hex_yn.1_nf_od5_radius_drag.nc)


radius=6363827
damp=1
gc=0

echo $radius
echo $damp
echo $gc

for ((i = 0; i < ${#odrag[@]}; ++i)); do
    drag=(${odrag[$i]})
    fname=(${name[$i]})

    echo $drag
    echo $fname

    sed  "s/<od>/$drag/; s/<gc>/$gc/; s/<rad>/$radius/; s/<damp>/$damp/; s/<name>/$fname/g" \
	 input_fill_hex.nml > input.nml
    ./RUN1p
    # &>/dev/null &
done
