!> Handles reading/writing of restart files and trajectory-based diagnostic files
module ice_bergs_io

! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.
#ifdef use_fms_io

use icebergs_fmsio
#else
use icebergs_fms2io

#endif

end module
