ncdump -v  VI ../ISOMIP_ocean_ice/Shelf/rho/00010101.ice_day.nc > a 
ncdump -v  VI ../ISOMIP_ocean_ice/Bergs/rho/00010101.ice_day.nc > b
diff a b > c
more c

