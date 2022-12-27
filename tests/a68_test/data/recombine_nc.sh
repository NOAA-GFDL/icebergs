
#using cdo
cdo mergetime ./split_files/ocean*.nc a68_experiment_ocean_surf_vel_oscar_dec2020_HOURLY_ll_p125.nc
cdo mergetime ./split_files/ssh*.nc a68_experiment_ssh_duacs_dec2020_HOURLY_ll_p125.nc
cdo mergetime ./split_files/wind*.nc a68_experiment_wind_vel_ncep_10m_dec2020_HOURLY_ll_p125.nc

#using ncrcat (nco)
#ncrcat ./split_files/ocean00000?.nc a68_experiment_ocean_surf_vel_oscar_dec2020_HOURLY_ll_p125.nc
#ncrcat ./split_files/ssh00000?.nc a68_experiment_ssh_duacs_dec2020_HOURLY_ll_p125.nc
#ncrcat ./split_files/wind00000?.nc a68_experiment_wind_vel_ncep_10m_dec2020_HOURLY_ll_p125.nc

ln -s a68_experiment_ocean_surf_vel_oscar_dec2020_HOURLY_ll_p125.nc ../../../postprocessing/data/a68_experiment_ocean_surf_vel_oscar_dec2020_HOURLY_ll_p125.nc
