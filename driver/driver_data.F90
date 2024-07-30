!> Handles reading/writing of restart files and trajectory-based diagnostic files
module driver_data
! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.

#ifndef USE_FMS2_IO
  use driver_data_fms
#else
  use driver_data_fms2
#endif

end module driver_data
