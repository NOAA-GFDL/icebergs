!> Handles reading/writing of restart files and trajectory-based diagnostic files
module ice_bergs_io
! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.

#ifndef USE_FMS2_IO
  use ice_bergs_fmsio
#else
  use ice_bergs_fms2io
#endif

end module
