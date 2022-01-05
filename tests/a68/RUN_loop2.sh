#!/bin/bash

#ocean drag values
#odrag=(20 22 24 26 28)
#name=(od_20_tn_ng_nw.nc od_22_tn_ng_nw.nc od_24_tn_ng_nw.nc od_26_tn_ng_nw.nc od_28_tn_ng_nw.nc)

# odrag=(1 5 10 15 20 22 24 26)
# name=(od_1_tn_ng_ellr.nc od_5_tn_ng_ellr.nc od_10_tn_ng_ellr.nc od_15_tn_ng_ellr.nc od_20_tn_ng_ellr.nc od_22_tn_ng_ellr.nc od_24_tn_ng_ellr.nc od_26_tn_ng_ellr.nc)

# odrag=(15 17.5 20)
# name=(od_15_tn_ng_noyd.nc od_17.5_tn_ng_noyd.nc od_20_tn_ng_noyd.nc)

#odrag=(19.7)
#name=(od_19.7_tn_ng_yd.05.nc)

#odrag=(13 14 15)
#name=(ydp.5_od_13.nc ydp.5_od_14.nc ydp.5_od_15.nc)

Rearth=6363827
odrag=(17.8 17.9)
name=(od_17.8_tn_ng_ydn.1.nc od_17.9_tn_ng_ydn.1.nc)

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
