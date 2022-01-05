#!/usr/bin/env python

import xarray as xr
fname='./square_for_restart_od17.8_ydn.1.nc'
test=xr.load_dataset(fname,decode_timedelta=False,decode_times=False)
test2=test.copy()
test2=test2.where(test2.day<test2.day.max(),drop=True)
test2.to_netcdf('./square_for_restart_od17.8_ydn.1_without_last_ts.nc', format='NETCDF3_CLASSIC')


# import xarray as xr
# fname='./sq_ssf_ssg1.e4_od17.8_addsm.true._rectx-37.5285_y-55.2166_n23_t10.nc'
# test=xr.load_dataset(fname,decode_timedelta=False,decode_times=False)
# test2=test.copy()
# test2=test2.where(test2.day<test2.day.max(),drop=True)
# test2.to_netcdf('./sq_ssf_ssg1.e4_od17.8_addsm.true._rectx-37.5285_y-55.2166_n23_t10_without_last_ts.nc', format='NETCDF3_CLASSIC')


# ncrcat square_for_restart_od17.8_ydn.1_without_last_ts.nc sq_ssf_ssg1.e4_od17.8_addsm.true._rectx-37.5285_y-55.2166_n23_t10_without_last_ts.nc sq_ssf_ssg1.e4_od17.8_addsm.true._rectx-37.5285_y-55.2166_n23_t10_full2.nc
