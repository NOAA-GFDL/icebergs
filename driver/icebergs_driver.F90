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
  use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN,GLOBAL_DATA_DOMAIN
  use mpp_mod, only : mpp_npes
  use mpp_mod, only : mpp_pe
  use mpp_mod, only : mpp_root_pe
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


  ! For parallelism
  integer :: nprocs !< Number of processes
  integer :: layout(2) !< Layout of MPI processes
  integer :: io_layout(2) !< Layout of I/O processes

  integer :: iflags= CYCLIC_GLOBAL_DOMAIN !< For i-direction flags
  integer :: jflags= 0  !< For j-direction flags
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
  integer :: halo=1 !< Width of halo in parallel decomposition
  real :: ibdt=3600.0
  real :: ibuo=0.1
  real :: ibvo=0.0
  real :: ibui=0.0
  real :: ibvi=0.0
  real :: gridres
  integer :: ibhrs=2
  integer :: nmax = 2000000000 !<max number of iterations    
  namelist /icebergs_driver_nml/ debug, ni, nj, halo, ibhrs, ibdt, ibuo, ibvo, nmax, &
       saverestart,ibui,ibvi,collision_test
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
  integer :: iyr, imon, iday, ihr, imin, isec
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
  integer :: i,j
  real, allocatable :: coord(:) !< One dimensional coordinate for initializing diagnostics

  ! For icebergs_run
  real, dimension(:,:), allocatable :: calving !<Calving (kg/s)
  real, dimension(:,:), allocatable :: uo,vo !<zonal and meridonal ocean velocities (m/s)
  real, dimension(:,:), allocatable :: ui,vi !<zonal and meridonal ice velocities (m/s)
  real, dimension(:,:), allocatable :: tauxa,tauya !< Zonal and meridonal wind stress (Pa)
  real, dimension(:,:), allocatable :: ssh !< Effective sea-surface height (m)
  real, dimension(:,:), allocatable :: sst !< Sea-surface temperature (C or K)
  real, dimension(:,:), allocatable :: cn !< Sea-ice concentration (nondim)
  real, dimension(:,:), allocatable :: hi !< Sea-ice thickness (m)
  real, dimension(:,:), allocatable :: calving_hflx !< Calving heat flux (W/m2)
  integer :: stagger, stress_stagger
  real, dimension(:,:), allocatable :: sss !< Sea-surface salinity (1e-3)
  real, dimension(:,:), pointer :: mass_berg !< Mass of bergs (kg)
  real, dimension(:,:), pointer  :: ustar_berg !< Friction velocity on base of bergs (m/s)
  real, dimension(:,:), pointer :: area_berg !< Area of bergs (m2)

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
     coord(i) =  real(i)
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
  if (debug .and. mpp_root_pe()==0) print *,''

  allocate( lon(isd:ied,jsd:jed) )
  allocate( lat(isd:ied,jsd:jed) )
  allocate( wet(isd:ied,jsd:jed) )
  allocate( dx(isd:ied,jsd:jed) )
  allocate( dy(isd:ied,jsd:jed) )
  allocate( area(isd:ied,jsd:jed) )
  allocate( cos_rot(isd:ied,jsd:jed) )
  allocate( sin_rot(isd:ied,jsd:jed) )
  allocate( depth(isd:ied,jsd:jed) )

  gridres = 1.e3


  do j = jsd, jed
     do i = isd, ied
        lon(i,j) = gridres*real(i)
        lat(i,j) = gridres*real(j)
        dx(i,j) = gridres
        dy(i,j) = gridres
        area(i,j) = gridres*gridres
        ! if (j>18.0 .or. j < 2.0) then
        !    wet(i,j)=0.0
        ! else
        wet(i,j) = 1.0
        ! end if
        cos_rot(i,j) = 1.0 
        sin_rot(i,j) = 0.0
        depth(i,j) = 1000.
     enddo
  enddo


  call set_calendar_type(THIRTY_DAY_MONTHS)
  time = set_date(1,1,1,0,0,0)
  dt = ibdt !3600.0

  call icebergs_init(bergs, &
       ni, nj, layout, io_layout, axes, iflags, jflags, dt, time, &
       lon(isc:iec,jsc:jec), lat(isc:iec,jsc:jec), wet(isc-1:iec+1,jsc-1:jec+1), &
       dx(isc-1:iec+1,jsc-1:jec+1), dy(isc-1:iec+1,jsc-1:jec+1), area(isc:iec,jsc:jec), &
       cos_rot(isc-1:iec+1,jsc-1:jec+1), sin_rot(isc-1:iec+1,jsc-1:jec+1), &
       depth(isc:iec,jsc:jec), fractional_area=.false.)

  ! for icebergs_run: 
  allocate( calving(isd:ied,jsd:jed) )
  allocate( uo(isd:ied,jsd:jed) )
  allocate( vo(isd:ied,jsd:jed) )
  allocate( ui(isd:ied,jsd:jed) )
  allocate( vi(isd:ied,jsd:jed) )
  allocate( tauxa(isd:ied,jsd:jed) )
  allocate( tauya(isd:ied,jsd:jed) )
  allocate( ssh(isd:ied,jsd:jed) )
  allocate( sst(isd:ied,jsd:jed) ) 
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
  sst = -2.0 !sea surface temperature (C or K; if K, will automatically adjust to C)
  cn = 0.0 !sea-ice concentration (nondim)
  hi = 0.0 !sea-ice thickness (m)
  calving_hflx = 0.0 !calving heat flux (W/m2)

  if (collision_test) then
    do j = jsc,jec
       if (j>10) vo = -ibvo
    enddo
    do i = isc,iec
       if (i>10) vo = 0.0
    enddo
  end if

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

  ns = 1 !timestep counter

  ! Time_end = increment_date(Time, years, months, days, hours, minutes, seconds)
  Time_end = increment_date(Time,0,0,0,ibhrs,0,0) 

  do while ((ns < nmax) .and. (Time < Time_end))

     if (mpp_pe()==0 .and. debug) then
        call get_date(Time, iyr, imon, iday, ihr, imin, isec)
        write(*,*),''
        write(*,'(a)'),'-------------------------------------------'
        write(*,'(a,i5)')' Timestep   ',ns
        write(*,'(a,i5,a,i5,a,i5)')' year',iyr,'  month ',imon,'  day   ',iday
        write(*,'(a,i5,a,i5,a,i5)')' hour',ihr,'  minute',imin,'  second',isec       
        write(*,'(a)'),'-------------------------------------------'
        write(*,*),''
     end if

     ! The main driver the steps updates icebergs
     call icebergs_run(bergs, time, calving(isc:iec,jsc:jec), &
          uo(isc-1:iec+1,jsc-1:jec+1), vo(isc-1:iec+1,jsc-1:jec+1), &
          ui(isc-1:iec+1,jsc-1:jec+1), vi(isc-1:iec+1,jsc-1:jec+1), &
          tauxa(isc:iec,jsc:jec), tauya(isc:iec,jsc:jec), &
          ssh(isc-1:iec+1,jsc-1:jec+1), sst(isc:iec,jsc:jec), &
          calving_hflx(isc:iec,jsc:jec), &
          cn(isc-1:iec+1,jsc-1:jec+1), hi(isc-1:iec+1,jsc-1:jec+1))

     ! full form, with optional vars:
     ! call icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, &
     ! calving_hflx, cn, hi,stagger, stress_stagger, sss, mass_berg, ustar_berg, area_berg)       

     ! if (ns == 1) then
     !    calving = 0.0
     ! endif


     Time = Time + real_to_time_type(dt)
     ns = ns+1
  enddo

  if (saverestart) call icebergs_save_restart(bergs)

  !Deallocate all memory and disassociated pointer
  call icebergs_end(bergs)

end program icebergs_driver
