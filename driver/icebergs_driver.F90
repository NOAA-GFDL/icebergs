program icebergs_driver

use constants_mod, only: pi
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
use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN,GLOBAL_DATA_DOMAIN
use mpp_mod, only : mpp_npes
use mpp_mod, only : mpp_pe
use mpp_mod, only : mpp_root_pe
use mpp_mod, only: mpp_sync
use mpp_mod, only: mpp_exit
use time_manager_mod, only : time_type
use time_manager_mod, only : set_date
use time_manager_mod, only : set_calendar_type
use time_manager_mod, only : THIRTY_DAY_MONTHS
use time_manager_mod, only : increment_date
use time_manager_mod, only : real_to_time_type,time_type_to_real
use time_manager_mod, only : operator(+), operator(-), operator(*), operator(/)
use time_manager_mod, only : operator(>), operator(<), operator(>=)
use diag_manager_mod, only : diag_manager_init
use time_manager_mod, only: get_date
use diag_axis_mod, only : diag_axis_init
use ice_bergs, only : icebergs
use ice_bergs, only : icebergs_init
use ice_bergs, only : icebergs_run
use ice_bergs, only : icebergs_end
use ice_bergs, only : icebergs_save_restart
use driver_data

! For parallelism
integer :: nprocs !< Number of processes
integer :: layout(2) !< Layout of MPI processes
integer :: io_layout(2) !< Layout of I/O processes
integer :: iflags=CYCLIC_GLOBAL_DOMAIN !< For i-direction flags
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
logical :: saverestart=.false. !< If true, save a berg restart file at the end
logical :: collision_test=.false.
logical :: chaotic_test=.false.
logical :: grounding_test=.false.
logical :: big_grounding_test=.false.
logical :: dem_test_contact=.false.
logical :: a68_test=.false.
logical :: vert_int_ocean_vel=.false.
logical :: tau_is_velocity=.true.
logical :: reverse_a68_forcings=.false.
character(len=1000) :: data_dir = 'data/'
integer :: transient_a68_data_start_ind=0 !<set to .gt. 0 to use the hourly JRA-55 wind data for the A68
!! test. The values indicate the time index to start from in the hourly data
integer :: halo=1 !< Width of halo in parallel decomposition
real :: ibdt=3600.0
real :: ibuo=0.0
real :: ibvo=0.0
real :: ibui=0.0
real :: ibvi=0.0
real :: gridres=1.e3
real :: bump_depth=0
real :: sst=-2
logical :: Old_wrong_Rearth_and_SSH=.false. !< use Rearth=6.36e6 and SSH on B grid
real :: REarth=6.378e6
!logical :: mom_Rearth=.false. !< Use Rearth=6.378e6
logical :: no_wind=.false.
integer :: ibhrs=2
integer :: nmax = 2000000000 !<max number of iteration
integer :: write_time_inc=1
namelist /icebergs_driver_nml/ debug, ni, nj, halo, ibhrs, ibdt, ibuo, ibvo, nmax, &
  saverestart,ibui,ibvi,collision_test,chaotic_test,grounding_test,&
  gridres,write_time_inc,bump_depth, sst, a68_test, data_dir, vert_int_ocean_vel,&
  reverse_a68_forcings,tau_is_velocity,transient_a68_data_start_ind,&
  Old_wrong_Rearth_and_SSH,no_wind,REarth,big_grounding_test,dem_test_contact !mom_Rearth
! For loops
integer :: isc !< Start of i-index for computational domain (used for loops)
integer :: iec !< End of i-index for computational domain (used for loops)
integer :: jsc !< Start of j-index for computational domain (used for loops)
integer :: jec !< End of j-index for computational domain (used for loops)
! Time
type(time_type) :: Time !< Model time
type(time_type) :: Time_end  !<end time for the experiment (s)
real :: dt !<time increment
integer :: ns !<timestep counter
real :: ns2
integer :: iyr, imon, iday, ihr, imin, isec
real :: time_begin, time_finish
! For icebergs
type(icebergs), pointer :: bergs !< Container for icebergs
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
integer :: i,j, dit
real, allocatable :: coord(:) !< One dimensional coordinate for initializing diagnostics

! For icebergs_run
real, dimension(:,:), allocatable :: calving !<Calving (kg/s)
real, dimension(:,:), allocatable :: uo,vo !<zonal and meridonal ocean velocities (m/s)
real, dimension(:,:), allocatable :: ui,vi !<zonal and meridonal ice velocities (m/s)
real, dimension(:,:), allocatable :: tauxa,tauya !< Zonal and meridonal wind stress (Pa)
real, dimension(:,:,:), allocatable :: tauxa_hr,tauya_hr !< Hourly jra-55 zonal and meridonal wind velocity (m/s) for A68 test
real, dimension(:,:,:), allocatable :: uo_hr,vo_hr !<Hourly OSCAR zonal and meridonal ocean surface velocities (m/s) for A68 test
real, dimension(:,:), allocatable :: ssh !< Effective sea-surface height (m)
real, dimension(:,:,:), allocatable :: ssh_hr !< Effective sea-surface height (m, hourly)
real, dimension(:,:), allocatable :: sstemp !< Sea-surface temperature (C or K)
real, dimension(:,:), allocatable :: cn !< Sea-ice concentration (nondim)
real, dimension(:,:), allocatable :: hi !< Sea-ice thickness (m)
real, dimension(:,:), allocatable :: calving_hflx !< Calving heat flux (W/m2)
integer :: stagger, stress_stagger
real, dimension(:,:), allocatable :: sss !< Sea-surface salinity (1e-3)
real, dimension(:,:), pointer :: mass_berg !< Mass of bergs (kg)
real, dimension(:,:), pointer  :: ustar_berg !< Friction velocity on base of bergs (m/s)
real, dimension(:,:), pointer :: area_berg !< Area of bergs (m2)
!grounding_test
real :: a,bx,by,c,xc,yc !gaussian params
!grounding/collision tests
real :: mid

! Boot FMS
call fms_init()

! Read run-time parameters
nmunit = open_namelist_file(file='input.nml')
read(nmunit, nml=icebergs_driver_nml)
call close_file(nmunit)

if (a68_test) then
  call a68_dims(data_dir,ni,nj)
  if (mpp_pe()==0) print *,'ni,nj',ni,nj
  iflags=0
endif

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
if (debug .and. mpp_root_pe()==0) print *,'isd,ied,jsd,jed=',isd,ied,jsd,jed

! Query for loop ranges
call mpp_get_compute_domain(mpp_domain, isc, iec, jsc, jec)
if (debug .and. mpp_root_pe()==0) print *,'isd,isc,iec,ied=',isd,isc,iec,ied
if (debug .and. mpp_root_pe()==0) print *,'jsd,jsc,jec,jed=',jsd,jsc,jec,jed

! for icebergs_init
allocate( lon(isd:ied,jsd:jed) )
allocate( lat(isd:ied,jsd:jed) )
allocate( wet(isd:ied,jsd:jed) )
allocate( dx(isd:ied,jsd:jed) )
allocate( dy(isd:ied,jsd:jed) )
allocate( area(isd:ied,jsd:jed) )
allocate( cos_rot(isd:ied,jsd:jed) )
allocate( sin_rot(isd:ied,jsd:jed) )
allocate( depth(isd:ied,jsd:jed) )

! for icebergs_run:
allocate( calving(isd:ied,jsd:jed) )
allocate( uo(isd:ied,jsd:jed) )
allocate( vo(isd:ied,jsd:jed) )
allocate( ui(isd:ied,jsd:jed) )
allocate( vi(isd:ied,jsd:jed) )
allocate( tauxa(isd:ied,jsd:jed) )
allocate( tauya(isd:ied,jsd:jed) )
allocate( ssh(isd:ied,jsd:jed) )
allocate( sstemp(isd:ied,jsd:jed) )
allocate( cn(isd:ied,jsd:jed) )
allocate( hi(isd:ied,jsd:jed) )
allocate( calving_hflx(isd:ied,jsd:jed) )

calving = 0.0 !"calving" flux (kg/s) due to "melt" of bergs
uo = ibuo !zonal ocean velocities (m/s)
vo = ibvo !meridonal ocean velocities (m/s)
ui = ibui !zonal ice velocities (m/s)
vi = ibvi !meridonal ice velocities (m/s)
tauxa = 0.0 !zonal wind stress (Pa)
tauya = 0.0 !meridonal wind stress (Pa)
ssh = 0.0 !eff sea-surf height
sstemp = sst !sea surface temperature (C or K; if K, will automatically adjust to C)
cn = 0.0 !sea-ice concentration (nondim)
hi = 0.0 !sea-ice thickness (m)
calving_hflx = 0.0 !calving heat flux (W/m2)

! optional:
! allocate( sss(isd:ied,jsd:jed) )
! allocate( mass_berg(isd:ied,jsd:jed) )
! allocate( ustar_berg(isd:ied,jsd:jed) )
! allocate( area_berg(isd:ied,jsd:jed) )
! stagger =  !scalar integer
! stress_stagger = !scalar integer
! sss = !sea surface salinity (1e-3)
! mass_berg = !mass of bergs (kg)
! ustar_berg = !Friction velocity on base of bergs (m/s)
! area_berg = !Area of bergs (m2)

if (a68_test) then
  if (Old_wrong_Rearth_and_SSH) REarth=6360000.
  call mpp_sync()
  call a68_prep(data_dir,mpp_domain,lon,lat,dx,dy,area,depth,uo,vo,tauxa,tauya,ssh,vert_int_ocean_vel,tau_is_velocity,Old_wrong_Rearth_and_SSH,REarth)
  !mom_Rearth)
  call mpp_sync()
  if (no_wind) then
    tauxa=0.0; tauya=0.0
  endif
  !where (lon<(-37.6+360.))
  !  depth=depth+1000.
  !end where
  wet=1.0
  cos_rot=1.0
  sin_rot=0.0
  if (reverse_a68_forcings) then
    uo=-uo; vo=-vo
    tauxa=-tauxa; tauya=-tauya
  endif

  if (transient_a68_data_start_ind>0) then
    if (.not. (ibdt==3600.0 .or. ibdt==1800.0)) then
      call error_mesg('icebergs_driver:', &
      'in order to use transient_a68_data_start_ind>0, timestep must be 30 min or 1 hr', FATAL)
    endif
    if (.not. tau_is_velocity) call error_mesg('icebergs_driver:', &
      'tau_is_velocity must be true if transient_a68_data_start_ind>-1', FATAL)
    allocate( tauxa_hr(isd:ied,jsd:jed,739) )
    allocate( tauya_hr(isd:ied,jsd:jed,739) )
    allocate( uo_hr(isd:ied,jsd:jed,744) )
    allocate( vo_hr(isd:ied,jsd:jed,744) )
    allocate( ssh_hr(isd:ied,jsd:jed,744) )
    call a68_prep_3d(data_dir,mpp_domain,tauxa_hr,tauya_hr,uo_hr,vo_hr,ssh_hr)
    if (reverse_a68_forcings) then
      tauxa_hr=-tauxa_hr; tauya_hr=-tauya_hr
      uo_hr=-uo_hr; vo_hr=-vo_hr
      ssh_hr=-ssh_hr
      dit=-1
    else
      dit=1
    endif
    if (no_wind) then
      tauxa_hr=0.0; tauya_hr=0.0
    endif
  endif
else
  do j = jsd, jed
    do i = isd, ied
      lon(i,j) = gridres*real(i)
      lat(i,j) = gridres*real(j)
      dx(i,j) = gridres
      dy(i,j) = gridres
      area(i,j) = gridres*gridres
      wet(i,j) = 1.0
      cos_rot(i,j) = 1.0
      sin_rot(i,j) = 0.0
      depth(i,j) = 1000.
    enddo
  enddo

  if (chaotic_test) then
    !assign land cells at the N and S portions of the domain.
    !input.nml: set coastal drift parameter to prevent grounding
    do j=jsd,jed; do i=isd,ied
      if (lat(i,j)<=0.0+2*gridres .or. lat(i,j)>=1000.e3-2*gridres) wet(i,j)=0.0
    enddo;enddo
  endif

  if (big_grounding_test) then
    lat(:,:)=lat(:,:)-0.45 !-25000.1
    lon(:,:)=lon(:,:)-0.45
    !assign land cells at the N and S portions of the domain.
    !input.nml: set coastal drift parameter to prevent grounding
    do j=jsd,jed; do i=isd,ied
      if (lat(i,j)<=(-5.e3) .or. lat(i,j)>=(220.e3)) wet(i,j)=0.0
    enddo;enddo

    !introduce gaussian bump to bathymetry
    a = 1000.0-bump_depth !max height(m)
    c = 5e3; !controls width
    bx = 63.e3; by = 60.e3 !50.e3
    do j=jsd,jed; do i=isd,ied
      xc=lon(i,j)-(gridres/2.); yc=lat(i,j)-(gridres/2.) !define at cell centers
      depth(i,j)=a*exp(-((xc-bx)*(xc-bx)/(2.*c*c)+(yc-by)*(yc-by)/(2.*c*c)))
    enddo; enddo
    print *,'min, max depth',minval(depth(:,:)),maxval(depth(:,:))
    depth = 1000.0-depth
  endif

   if (dem_test_contact) then
    lat(:,:)=lat(:,:)-0.45 !-25000.1
    lon(:,:)=lon(:,:)-0.45
    !assign land cells at the N and S portions of the domain.
    !input.nml: set coastal drift parameter to prevent grounding
    do j=jsd,jed; do i=isd,ied
      if (lat(i,j)<=(-5.e3) .or. lat(i,j)>=(220.e3)) wet(i,j)=0.0
    enddo;enddo

    mid=50.e3

    do j=jsd,jed; do i=isd,ied
      if (lon(i,j)>65 .or. lat(i,j)==mid) then
        vo(i,j)=0.0
      else
        if (lat(i,j)>mid) then
          vo(i,j)=-ibvo
        else
          vo(i,j)=ibvo
        endif
      endif
  enddo;enddo


  endif

  if (grounding_test) then
    lat(:,:)=lat(:,:)-15000.0
    !assign land cells at the N and S portions of the domain.
    !input.nml: set coastal drift parameter to prevent grounding
    do j=jsd,jed; do i=isd,ied
      if (lat(i,j)<=0.0+2*gridres .or. lat(i,j)>=30.e3-2*gridres) wet(i,j)=0.0
    enddo;enddo

    !introduce gaussian bump to bathymetry
    a = 1000.0-bump_depth !max height(m)
    c = 2.5e3; !controls width
    bx = 15.e3; by = 15.e3
    do j=jsd,jed; do i=isd,ied
      xc=lon(i,j)-(gridres/2.); yc=lat(i,j)-(gridres/2.) !define at cell centers
      depth(i,j)=a*exp(-((xc-bx)*(xc-bx)/(2.*c*c)+(yc-by)*(yc-by)/(2.*c*c)))
    enddo; enddo
    print *,'min, max depth',minval(depth(:,:)),maxval(depth(:,:))
    depth = 1000.0-depth
  endif
endif


call set_calendar_type(THIRTY_DAY_MONTHS)

if (transient_a68_data_start_ind==355) then
  time = set_date(1,1,7,5,0,0)
else
  time = set_date(1,1,1,0,0,0)
endif
dt = ibdt

call icebergs_init(bergs, &
                   ni, nj, layout, io_layout, axes, iflags, jflags, dt, time, &
                   lon(isc:iec,jsc:jec), lat(isc:iec,jsc:jec), wet(isc-1:iec+1,jsc-1:jec+1), &
                   dx(isc-1:iec+1,jsc-1:jec+1), dy(isc-1:iec+1,jsc-1:jec+1), area(isc:iec,jsc:jec), &
                   cos_rot(isc-1:iec+1,jsc-1:jec+1), sin_rot(isc-1:iec+1,jsc-1:jec+1), &
                   depth(isc:iec,jsc:jec), fractional_area=.false.)

if (collision_test .or. grounding_test) then
  if (grounding_test) then
    mid=15.e3
  else
    mid=10.e3
  endif

  do j=jsd,jed; do i=isd,ied
    if (lon(i,j)>mid .or. lon(i,j)<=0.0 .or. lat(i,j)==mid) then
      vo(i,j)=0.0
      ! if (grounding_test .and. lon(i,j)>mid+5.e3) then
      !   if (lat(i,j)>mid) then
      !     vo(i,j)=+ibvo
      !   else
      !     vo(i,j)=-ibvo
      !   endif
      ! endif

    else
      if (lat(i,j)>mid) then
        vo(i,j)=-ibvo
      else
        vo(i,j)=ibvo
      endif
    endif
  enddo;enddo
endif

if (chaotic_test) then
  do j=jsd,jed; do i=isd,ied
    if (lat(i,j)==500.e3 .or. lon(i,j)==500.e3) then
      vo(i,j)=0.0
    elseif (lon(i,j)>500.e3) then
      if (lat(i,j)>500.e3) then
        vo(i,j)=0.5*ibvo
      else
        vo(i,j)=-0.5*ibvo
      endif
    else
      if (lat(i,j)>500.e3) then
        vo(i,j)=-ibvo
      else
        vo(i,j)=ibvo
      endif
    endif
  enddo;enddo
endif

ns = 1 !timestep counter
if (transient_a68_data_start_ind>0) ns2 = 1

! Time_end = increment_date(Time, years, months, days, hours, minutes, seconds)
Time_end = increment_date(Time,0,0,0,ibhrs,0,0)

call cpu_time(time_begin)

do while ((ns <= nmax) .and. (Time < Time_end))
  if (mpp_pe()==0 .and. debug .and. mod(ns-1,write_time_inc)==0) then
    call get_date(Time, iyr, imon, iday, ihr, imin, isec)
    write(*,*) ''
    write(*,'(a)') '-------------------------------------------'
    write(*,'(a,i5)')' Timestep   ',ns
    write(*,'(a,i5,a,i5,a,i5)')' year',iyr,'  month ',imon,'  day   ',iday
    write(*,'(a,i5,a,i5,a,i5)')' hour',ihr,'  minute',imin,'  second',isec
    call cpu_time(time_finish)
    write(*,*) 'clock-time elapsed: ', (time_finish - time_begin)/60., 'minutes'
    write(*,'(a)') '-------------------------------------------'
    write(*,*) ''
  end if

  if (transient_a68_data_start_ind>0) then

    if ((ibdt==3600.0) .or. (mod(ns2,1.)==0)) then
      tauxa=tauxa_hr(:,:,transient_a68_data_start_ind+dit*(int(ns2)-1))
      tauya=tauya_hr(:,:,transient_a68_data_start_ind+dit*(int(ns2)-1))
      uo=uo_hr(:,:,transient_a68_data_start_ind+dit*(int(ns2)-1))
      vo=vo_hr(:,:,transient_a68_data_start_ind+dit*(int(ns2)-1))
      ssh=ssh_hr(:,:,transient_a68_data_start_ind+dit*(int(ns2)-1))
    else
      tauxa=0.5*(tauxa+tauxa_hr(:,:,transient_a68_data_start_ind+dit*(int(ceiling(ns2))-1)))
      tauya=0.5*(tauya+tauya_hr(:,:,transient_a68_data_start_ind+dit*(int(ceiling(ns2))-1)))
      uo=0.5*(uo+uo_hr(:,:,transient_a68_data_start_ind+dit*(int(ceiling(ns2))-1)))
      vo=0.5*(vo+vo_hr(:,:,transient_a68_data_start_ind+dit*(int(ceiling(ns2))-1)))
      ssh=ssh_hr(:,:,transient_a68_data_start_ind+dit*(int(ns2)-1))
    endif
  endif

  ! The main driver the steps updates icebergs
  call icebergs_run(bergs, time, calving(isc:iec,jsc:jec), &
                    uo(isc-1:iec+1,jsc-1:jec+1), vo(isc-1:iec+1,jsc-1:jec+1), &
                    ui(isc-1:iec+1,jsc-1:jec+1), vi(isc-1:iec+1,jsc-1:jec+1), &
                    tauxa(isc:iec,jsc:jec), tauya(isc:iec,jsc:jec), &
                    ssh(isc-1:iec+1,jsc-1:jec+1), sstemp(isc:iec,jsc:jec), &
                    calving_hflx(isc:iec,jsc:jec), &
                    cn(isc-1:iec+1,jsc-1:jec+1), hi(isc-1:iec+1,jsc-1:jec+1))

  ! full form, with optional vars:
  ! call icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sstemp, &
  ! calving_hflx, cn, hi,stagger, stress_stagger, sss, mass_berg, ustar_berg, area_berg)

  Time = Time + real_to_time_type(dt)
  ns = ns+1

  if (transient_a68_data_start_ind>0) then
    if (ibdt==3600.0) then
      ns2=ns2+1
    elseif (ibdt==1800.0) then
      ns2=ns2+0.5
    endif
  endif
enddo


if (a68_test) then
  call cpu_time(time_finish)
  call get_date(Time, iyr, imon, iday, ihr, imin, isec)
  write(*,*) ''
  write(*,'(a)') '-------------------------------------------'
  write(*,'(a,i5)')' Saving after timestep   ',ns-1
  print *,' ns2   ',ns2
  if (mod(ns2,1.)==0) then
    print *,'transient_a68_data_start_ind+dit*(int(ns2)-1)',&
      transient_a68_data_start_ind+dit*(int(ns2)-1)
  else
    print *,'transient_a68_data_start_ind+dit*(int(ceiling(ns2))-1))',&
      transient_a68_data_start_ind+dit*(int(ceiling(ns2))-1)
  endif
  write(*,'(a)')' Restart time is:'
  write(*,'(a,i5,a,i5,a,i5)')' year',iyr,'  month ',imon,'  day   ',iday
  write(*,'(a,i5,a,i5,a,i5)')' hour',ihr,'  minute',imin,'  second',isec
  write(*,*) ''
  write(*,*) 'clock-time elapsed: ', (time_finish - time_begin)/60., 'minutes'
  write(*,*) 'clock-time per day: ', (time_finish - time_begin)/(ns-1)*(60*60*24)/ibdt
  write(*,'(a)') '-------------------------------------------'
  write(*,*) ''
endif

if (saverestart) call icebergs_save_restart(bergs)

call cpu_time(time_finish)

if (mpp_pe()==0) then
  write(*,*) 'Simulation finished in ', (time_finish - time_begin)/60., 'minutes'
end if

call mpp_exit()

!Deallocate all memory and disassociated pointer
call icebergs_end(bergs)

end program icebergs_driver
