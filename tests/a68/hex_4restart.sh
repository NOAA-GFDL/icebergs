#!/bin/bash

#ocean drag values
#odrag=(6.5 6 5 7 4 8 3 1)
#name=(hex_wgs84_nofrac_od6.5.nc hex_wgs84_nofrac_od6.nc hex_wgs84_nofrac_od5.nc hex_wgs84_nofrac_od7.nc hex_wgs84_nofrac_od4.nc hex_wgs84_nofrac_od8.nc hex_wgs84_nofrac_od3.nc hex_wgs84_nofrac_od1.nc)

#odrag=(7.2)
#name=(hex_ydn.05_wgs84_nf_od7.2_4restart.nc)

# odrag=(6)
# name=(hex_yn.075_nf_od6_4restart.nc)
# bond_name=(hex_yn.075_nf_od6_4restart_bonds.nc)
# radius=6363827

# odrag=(5.4)
# name=(hex_yn.1_nf_od5.4_4restart.nc)
# bond_name=(hex_yn.1_nf_od5.4_4restart_bonds.nc)
# radius=6363827

odrag=(5)
name=(hex_yn.1_nf_od5_4restart.nc)
bond_name=(hex_yn.1_nf_od5_4restart_bonds.nc)
radius=6363827

echo $radius

for ((i = 0; i < ${#odrag[@]}; ++i)); do
    drag=(${odrag[$i]})
    fname=(${name[$i]})
    bname=(${bond_name[$i]})
    echo $drag
    echo $fname

    sed  "s/<od>/$drag/; s/<rad>/$radius/; s/<bname>/$bname/; s/<name>/$fname/g" \
	 input_fill_hex_4restart.nml > input.nml
    rm $fname
    rm $bname
    ./RUN1p
    # &>/dev/null &
done
