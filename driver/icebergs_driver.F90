program icebergs_driver

use fms_mod, only : fms_init
use fms_mod, only : fms_end
use fms_mod, only : open_namelist_file
use fms_mod, only : close_file
use fms_mod, only : error_mesg
use fms_mod, only : FATAL
use mpp_domains_mod, only : mpp_define_layout
use mpp_domains_mod, only : mpp_define_domains
use mpp_domains_mod, only : mpp_define_io_domain
use mpp_domains_mod, only : mpp_get_compute_domain
use mpp_domains_mod, only : mpp_get_data_domain
use mpp_domains_mod, only : domain2D
use mpp_mod, only : mpp_npes
use mpp_mod, only : mpp_pe
use mpp_mod, only : mpp_root_pe
use time_manager_mod, only : time_type
use time_manager_mod, only : set_date
use time_manager_mod, only : set_calendar_type
use time_manager_mod, only : THIRTY_DAY_MONTHS
use diag_manager_mod, only : diag_manager_init
use diag_axis_mod, only : diag_axis_init

use ice_bergs, only : icebergs
use ice_bergs, only : icebergs_init
use ice_bergs, only : icebergs_run
use ice_bergs, only : icebergs_end

! For parallelism
integer :: nprocs !< Number of processes
integer :: layout(2) !< Layout of MPI processes
integer :: io_layout(2) !< Layout of I/O processes
integer :: iflags=0 !< For i-direction flags
integer :: jflags=0 !< For j-direction flags
type(domain2D) :: mpp_domain !< Parallel decomposition
! For declarations
integer :: isd !< Start of i-index for data domain (used in declarations)
integer :: ied !< End of i-index for data domain (used in declarations)
integer :: jsd !< Start of j-index for data domain (used in declarations)
integer :: jed !< End of j-index for data domain (used in declarations)
! I/O
integer :: nmunit !< Namelist file unit
integer :: axes(2) !< Diagnostic axes
! Run-time parameters for this driver
integer :: ni=0 !< Global number of grid cells in i-direction (mandatory input)
integer :: nj=0 !< Global number of grid cells in j-direction (mandatory input)
logical :: debug=.true. !< If true, do some debugging
integer :: halo=1 !< Width of halo in parallel decomposition
namelist /icebergs_driver_nml/ debug, ni, nj, halo
! For loops
integer :: isc !< Start of i-index for computational domain (used for loops)
integer :: iec !< End of i-index for computational domain (used for loops)
integer :: jsc !< Start of j-index for computational domain (used for loops)
integer :: jec !< End of j-index for computational domain (used for loops)
! For icebergs
type(icebergs), pointer :: bergs !< Container for icebergs
type (time_type) :: time !< Model time
real, allocatable :: lon(:,:) !< Longitude or position in x
real, allocatable :: lat(:,:) !< Latitude or position in y
real, allocatable :: wet(:,:) !< Wet-dry mask to indicate ocean/land
real, allocatable :: dx(:,:) !< Length of northern edge of cell (m)
real, allocatable :: dy(:,:) !< Length of eatern edge of cell (m)
real, allocatable :: area(:,:) !< Area of cell (m2)
real, allocatable :: cos_rot(:,:) !< Cosine of angle of grid orientation
real, allocatable :: sin_rot(:,:) !< Sine of angle of grid orientation
real, allocatable :: depth(:,:) !< Depth of ocean (m)
! Work variables
integer :: i,j
real, allocatable :: coord(:) !< One dimensional coordinate for initializing diagnostics

! Boot FMS
call fms_init()

! Read run-time parameters
nmunit = open_namelist_file(file='input.nml')
read(nmunit, nml=icebergs_driver_nml)
call close_file(nmunit)
if (min(ni,nj)<1) call error_mesg('icebergs_driver:', &
  'Both ni,nj must be positive in icebergs_driver_nml namelist', FATAL)

! Ask FMS to figure our how to decompose the global ni x nj grid
layout = (/0,0/) ! 0,0 lets FMS find an optimal layout
nprocs = mpp_npes() ! How many processes are we using?
call mpp_define_layout((/1, ni, 1, nj/), nprocs, layout)
if (debug .and. mpp_root_pe()==0) print *,'layout=',layout,'nprocs=',nprocs,'PE=',mpp_pe()

! Create a decomposed domain with the given layout from above
call mpp_define_domains((/1,ni,1,nj/), layout, mpp_domain, &
         xflags=iflags, yflags=jflags, &
         xhalo=halo, yhalo=halo, &
         symmetry=.false., name='iceberg_driver')

! Set up I/O
io_layout = (/1,1/) ! 1,1 tells FMS to use just one processor for I/O
call mpp_define_io_domain(mpp_domain, io_layout)

! Diagnostics manager
call diag_manager_init()
! Axes for diagnostics
allocate( coord(ni) )
do i = 1, ni
  coord(i) = real(i)
enddo
axes(1) = diag_axis_init('i', coord, 'index', 'X', 'cell index i', Domain2=mpp_domain)
deallocate( coord )
allocate( coord(nj) )
do j = 1, nj
  coord(j) = real(j)
enddo
axes(2) = diag_axis_init('j', coord, 'index', 'Y', 'cell index j', Domain2=mpp_domain)
deallocate( coord )

! Query for declaration ranges
call mpp_get_data_domain(mpp_domain, isd, ied, jsd, jed)
if (debug .and. mpp_root_pe()==0) print *,'isd,ied,jsc,jed=',isd,ied,jsc,jed

! Query for loop ranges
call mpp_get_compute_domain(mpp_domain, isc, iec, jsc, jec)
if (debug .and. mpp_root_pe()==0) print *,'isd,isc,iec,ied=',isd,isc,iec,ied
if (debug .and. mpp_root_pe()==0) print *,'jsd,jsc,jec,jed=',jsd,jsc,jec,jed

allocate( lon(isd:ied,jsd:jed) )
allocate( lat(isd:ied,jsd:jed) )
allocate( wet(isd:ied,jsd:jed) )
allocate( dx(isd:ied,jsd:jed) )
allocate( dy(isd:ied,jsd:jed) )
allocate( area(isd:ied,jsd:jed) )
allocate( cos_rot(isd:ied,jsd:jed) )
allocate( sin_rot(isd:ied,jsd:jed) )
allocate( depth(isd:ied,jsd:jed) )

do j = jsc, jec
  do i = isc, iec
    lon(i,j) = 10.e3*real(i)
    lat(i,j) = 10.e3*real(j)
    dx(i,j) = 10.e3
    dy(i,j) = 10.e3
    area(i,j) = 1.e8
    wet(i,j) = 1.0
    cos_rot(i,j) = 1.0
    sin_rot(i,j) = 0.0
    depth(i,j) = 1000.
  enddo
enddo

call set_calendar_type(THIRTY_DAY_MONTHS)
time = set_date(1,1,1,0,0,0)

call icebergs_init(bergs, &
                   ni, nj, layout, io_layout, axes, iflags, jflags, dt, time, &
                   lon(isc:iec,jsc:jec), lat(isc:iec,jsc:jec), wet, dx, dy, area(isc:iec,jsc:jec), &
                   cos_rot, sin_rot, depth(isc:iec,jsc:jec), fractional_area=.false.)

end program icebergs_driver
