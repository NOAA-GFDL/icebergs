module ice_bergs_io

use mpp_domains_mod, only: domain2D
use mpp_domains_mod, only: mpp_domain_is_tile_root_pe,mpp_get_domain_tile_root_pe
use mpp_domains_mod, only: mpp_get_tile_pelist,mpp_get_tile_npes,mpp_get_io_domain,mpp_get_tile_id

use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_gather, mpp_chksum
use mpp_mod, only: COMM_TAG_11, COMM_TAG_12, COMM_TAG_13, COMM_TAG_14

use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING, NOTE
use fms_mod, only: field_exist, file_exist, read_data, write_data

use fms_io_mod, only: get_instance_filename
use fms_io_mod, only : save_restart, restart_file_type, free_restart_type, set_meta_global
use fms_io_mod, only : register_restart_axis, register_restart_field, set_domain, nullify_domain
use fms_io_mod, only : read_unlimited_axis =>read_compressed, field_exist, get_field_size

use mpp_mod,    only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod,    only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP
use fms_mod,    only : clock_flag_default

use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)

use ice_bergs_framework, only: icebergs_gridded, xyt, iceberg, icebergs, buffer, bond
use ice_bergs_framework, only: pack_berg_into_buffer2,unpack_berg_from_buffer2
use ice_bergs_framework, only: pack_traj_into_buffer2,unpack_traj_from_buffer2
use ice_bergs_framework, only: find_cell,find_cell_by_search,count_bergs,is_point_in_cell,pos_within_cell,append_posn
use ice_bergs_framework, only: count_bonds, form_a_bond, find_individual_iceberg
use ice_bergs_framework, only: push_posn
use ice_bergs_framework, only: add_new_berg_to_list,destroy_iceberg
use ice_bergs_framework, only: increase_ibuffer,grd_chksum2,grd_chksum3
use ice_bergs_framework, only: sum_mass,sum_heat,bilin
!params !Niki: write a subroutine to get these
use ice_bergs_framework, only: nclasses, buffer_width, buffer_width_traj
use ice_bergs_framework, only: verbose, really_debug, debug, restart_input_dir,make_calving_reproduce
use ice_bergs_framework, only: ignore_ij_restart, use_slow_find,generate_test_icebergs,print_berg
use ice_bergs_framework, only: force_all_pes_traj
use ice_bergs_framework, only: check_for_duplicates_in_parallel

implicit none ; private

include 'netcdf.inc'

public ice_bergs_io_init
public read_restart_bergs,read_restart_bergs_orig,write_restart,write_trajectory
public read_restart_calving, read_restart_bonds
public read_ocean_depth

!Local Vars
integer, parameter :: file_format_major_version=0
integer, parameter :: file_format_minor_version=1
!I/O vars
type(domain2d), pointer, save :: io_domain=>NULL()
integer, save :: io_tile_id(1), io_tile_root_pe, io_npes
integer, allocatable,save :: io_tile_pelist(:)
logical :: is_io_tile_root_pe = .true.

integer :: clock_trw,clock_trp

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

subroutine ice_bergs_io_init(bergs, io_layout)
type(icebergs), pointer :: bergs
integer, intent(in) :: io_layout(2)

integer :: np
integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "ice_bergs_framework: "//trim(version)

  !I/O layout init 
  io_tile_id=-1
  io_domain => mpp_get_io_domain(bergs%grd%domain)
  if(associated(io_domain)) then
     io_tile_id = mpp_get_tile_id(io_domain)
     is_io_tile_root_pe = mpp_domain_is_tile_root_pe(io_domain)
     io_tile_root_pe = mpp_get_domain_tile_root_pe(io_domain)
     np=mpp_get_tile_npes(io_domain)
     allocate(io_tile_pelist(np))
     call mpp_get_tile_pelist(io_domain,io_tile_pelist)
     io_npes = io_layout(1)*io_layout(2)
  endif

  clock_trw=mpp_clock_id( 'Icebergs-traj write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  clock_trp=mpp_clock_id( 'Icebergs-traj prepare', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

end subroutine ice_bergs_io_init

! ##############################################################################

subroutine write_restart(bergs)
! Arguments
type(icebergs), pointer :: bergs
type(bond), pointer :: current_bond
! Local variables
integer :: i,j,id
character(len=35) :: filename
character(len=35) :: filename_bonds
type(iceberg), pointer :: this=>NULL()
integer :: stderrunit
!I/O vars
type(restart_file_type) :: bergs_restart
type(restart_file_type) :: bergs_bond_restart
integer :: nbergs, nbonds
integer :: n_static_bergs
logical :: check_bond_quality 
type(icebergs_gridded), pointer :: grd
real, allocatable, dimension(:) :: lon,          &
                                   lat,          &
                                   uvel,         &
                                   vvel,         &
                                   mass,         &
                                   axn,          &
                                   ayn,          &
                                   bxn,          &
                                   byn,          &
                                   thickness,    &
                                   width,        &
                                   length,       &
                                   start_lon,    &
                                   start_lat,    &
                                   start_day,    &
                                   start_mass,   &
                                   mass_scaling, &
                                   mass_of_bits, &
                                   static_berg,  &
                                   heat_density

integer, allocatable, dimension(:) :: ine,              &
                                      jne,              &
                                      iceberg_num,      &
                                      start_year,       &
                                      first_berg_num,   &
                                      other_berg_num,   &
                                      first_berg_jne,         &
                                      first_berg_ine,         &
                                      other_berg_jne,         &
                                      other_berg_ine


integer :: grdi, grdj
  
! Get the stderr unit number
 stderrunit=stderr()


  ! For convenience
  grd=>bergs%grd
 
  !First add the bergs on the io_tile_root_pe (if any) to the I/O list
  nbergs = 0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      nbergs = nbergs +1
      this=>this%next
    enddo
  enddo ; enddo

   allocate(lon(nbergs))
   allocate(lat(nbergs))
   allocate(uvel(nbergs))
   allocate(vvel(nbergs))
   allocate(mass(nbergs))
   allocate(axn(nbergs))    !Alon
   allocate(ayn(nbergs))    !Alon
   allocate(bxn(nbergs)) !Alon
   allocate(byn(nbergs)) !Alon
   allocate(thickness(nbergs))
   allocate(width(nbergs))
   allocate(length(nbergs))
   allocate(start_lon(nbergs))
   allocate(start_lat(nbergs))
   allocate(start_day(nbergs))
   allocate(start_mass(nbergs))
   allocate(mass_scaling(nbergs))
   allocate(mass_of_bits(nbergs))
   allocate(heat_density(nbergs))
   allocate(static_berg(nbergs))

   allocate(ine(nbergs))
   allocate(jne(nbergs))
   allocate(start_year(nbergs))
   allocate(iceberg_num(nbergs))


  call get_instance_filename("icebergs.res.nc", filename)
  call set_domain(bergs%grd%domain)
  call register_restart_axis(bergs_restart,filename,'i',nbergs)
  call set_meta_global(bergs_restart,'file_format_major_version',ival=(/file_format_major_version/))
  call set_meta_global(bergs_restart,'file_format_minor_version',ival=(/file_format_minor_version/))
  call set_meta_global(bergs_restart,'time_axis',ival=(/0/))

  !Now start writing in the io_tile_root_pe if there are any bergs in the I/O list

  ! Define Variables
  id = register_restart_field(bergs_restart,filename,'lon',lon,longname='longitude',units='degrees_E')
  id = register_restart_field(bergs_restart,filename,'lat',lat,longname='latitude',units='degrees_N')
  id = register_restart_field(bergs_restart,filename,'uvel',uvel,longname='zonal velocity',units='m/s')
  id = register_restart_field(bergs_restart,filename,'vvel',vvel,longname='meridional velocity',units='m/s')
  id = register_restart_field(bergs_restart,filename,'mass',mass,longname='mass',units='kg')
  if (.not. bergs%Runge_not_Verlet) then
    id = register_restart_field(bergs_restart,filename,'axn',axn,longname='explicit zonal acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'ayn',ayn,longname='explicit meridional acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'bxn',bxn,longname='inplicit zonal acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'byn',byn,longname='implicit meridional acceleration',units='m/s^2')
  endif
  id = register_restart_field(bergs_restart,filename,'ine',ine,longname='i index',units='none')
  id = register_restart_field(bergs_restart,filename,'jne',jne,longname='j index',units='none')
  id = register_restart_field(bergs_restart,filename,'thickness',thickness,longname='thickness',units='m')
  id = register_restart_field(bergs_restart,filename,'width',width,longname='width',units='m')
  id = register_restart_field(bergs_restart,filename,'length',length,longname='length',units='m')
  id = register_restart_field(bergs_restart,filename,'start_lon',start_lon, &
                                            longname='longitude of calving location',units='degrees_E')
  id = register_restart_field(bergs_restart,filename,'start_lat',start_lat, &
                                            longname='latitude of calving location',units='degrees_N')
  id = register_restart_field(bergs_restart,filename,'start_year',start_year, &
                                            longname='calendar year of calving event', units='years')
  id = register_restart_field(bergs_restart,filename,'iceberg_num',iceberg_num, &
                                            longname='identification of the iceberg', units='dimensionless')
  id = register_restart_field(bergs_restart,filename,'start_day',start_day, &
                                            longname='year day of calving event',units='days')
  id = register_restart_field(bergs_restart,filename,'start_mass',start_mass, &
                                            longname='initial mass of calving berg',units='kg')
  id = register_restart_field(bergs_restart,filename,'mass_scaling',mass_scaling, &
                                            longname='scaling factor for mass of calving berg',units='none')
  id = register_restart_field(bergs_restart,filename,'mass_of_bits',mass_of_bits, &
                                            longname='mass of bergy bits',units='kg')
  id = register_restart_field(bergs_restart,filename,'heat_density',heat_density, &
                                            longname='heat density',units='J/kg')

  !Checking if any icebergs are static in order to decide whether to save static_berg
  n_static_bergs = 0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      n_static_bergs=n_static_bergs+this%static_berg
      this=>this%next
    enddo
  enddo ; enddo
  call mpp_sum(n_static_bergs) 
  if (n_static_bergs .gt. 0) &
    id = register_restart_field(bergs_restart,filename,'static_berg',static_berg, &
                                              longname='static_berg',units='dimensionless')

  ! Write variables
   
  i = 0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))
      i = i + 1
      lon(i) = this%lon; lat(i) = this%lat
      uvel(i) = this%uvel; vvel(i) = this%vvel
      ine(i) = this%ine; jne(i) = this%jne
      mass(i) = this%mass; thickness(i) = this%thickness
      axn(i) = this%axn; ayn(i) = this%ayn !Added by Alon
      bxn(i) = this%bxn; byn(i) = this%byn !Added by Alon
      width(i) = this%width; length(i) = this%length
      start_lon(i) = this%start_lon; start_lat(i) = this%start_lat
      start_year(i) = this%start_year; start_day(i) = this%start_day
      start_mass(i) = this%start_mass; mass_scaling(i) = this%mass_scaling
      static_berg(i) = this%static_berg 
      iceberg_num(i) = this%iceberg_num; 
      mass_of_bits(i) = this%mass_of_bits; heat_density(i) = this%heat_density
      this=>this%next
    enddo
  enddo ; enddo


  call save_restart(bergs_restart)
  call free_restart_type(bergs_restart)

  deallocate(              &
             lon,          &
             lat,          &
             uvel,         &
             vvel,         &
             mass,         &
             axn,          &
             ayn,          &
             bxn,          &
             byn,          &
             thickness,    &
             width,        &
             length,       &
             start_lon,    &
             start_lat,    &
             start_day,    &
             start_mass,   &
             mass_scaling, &
             mass_of_bits, &
             static_berg,    &
             heat_density )

  deallocate(           &
             ine,       &
             jne,       &
             iceberg_num,       &
             start_year )

  call nullify_domain()

!########## Creating bond restart file ######################

   !Allocating restart memory for bond related variables.
   nbonds=0
   if (bergs%iceberg_bonds_on) then
     check_bond_quality=.true.
     call count_bonds(bergs, nbonds,check_bond_quality)

   allocate(first_berg_num(nbonds))
   allocate(other_berg_num(nbonds))
   allocate(first_berg_ine(nbonds))
   allocate(first_berg_jne(nbonds))
   allocate(other_berg_ine(nbonds))
   allocate(other_berg_jne(nbonds))

  call get_instance_filename("bonds_iceberg.res.nc", filename_bonds)
  call set_domain(bergs%grd%domain)
  call register_restart_axis(bergs_bond_restart,filename,'i',nbonds)
  call set_meta_global(bergs_bond_restart,'file_format_major_version',ival=(/file_format_major_version/))
  call set_meta_global(bergs_bond_restart,'file_format_minor_version',ival=(/file_format_minor_version/))
  call set_meta_global(bergs_bond_restart,'time_axis',ival=(/0/))

  !Now start writing in the io_tile_root_pe if there are any bergs in the I/O list

  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_berg_ine',first_berg_ine,longname='iceberg ine of first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_berg_jne',first_berg_jne,longname='iceberg jne of first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_berg_num',first_berg_num,longname='iceberg id first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_berg_ine',other_berg_ine,longname='iceberg ine of second berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_berg_jne',other_berg_jne,longname='iceberg jne of second berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_berg_num',other_berg_num,longname='iceberg id second berg in bond',units='dimensionless')
  
  
  ! Write variables
   
  i = 0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this)) !Loops over all bergs
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        i = i + 1
        first_berg_ine(i)=this%ine
        first_berg_jne(i)=this%jne
        first_berg_num(i)= this%iceberg_num
        other_berg_num(i)=current_bond%other_berg_num
        other_berg_ine(i)=current_bond%other_berg%ine
        other_berg_jne(i)=current_bond%other_berg%jne

        current_bond=>current_bond%next_bond
      enddo !End of loop over bonds
      this=>this%next
    enddo!End of loop over bergs
  enddo; enddo !End of loop over grid

  call save_restart(bergs_bond_restart)
  call free_restart_type(bergs_bond_restart)


  deallocate(                 &
             first_berg_num,  &
             other_berg_num, &
             first_berg_ine,        &
             first_berg_jne,        &
             other_berg_ine,        &
             other_berg_jne )


  call nullify_domain()
   endif
!#############################################################################################

  ! Write stored ice
  filename='RESTART/calving.res.nc'
  if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(stderrunit,'(2a)') 'diamonds, write_restart: writing ',filename

  call grd_chksum3(bergs%grd, bergs%grd%stored_ice, 'write stored_ice')
  call write_data(filename, 'stored_ice', bergs%grd%stored_ice, bergs%grd%domain)
  call grd_chksum2(bergs%grd, bergs%grd%stored_heat, 'write stored_heat')
  call write_data(filename, 'stored_heat', bergs%grd%stored_heat, bergs%grd%domain)
  !call grd_chksum2(bergs%grd, bergs%grd%iceberg_counter_grd, 'write iceberg_counter_grd')
  call write_data(filename, 'iceberg_counter_grd', bergs%grd%iceberg_counter_grd, bergs%grd%domain)
  if (bergs%tau_calving>0.) then
    call grd_chksum2(bergs%grd, bergs%grd%rmean_calving, 'write mean calving')
    call write_data(filename, 'rmean_calving', bergs%grd%rmean_calving, bergs%grd%domain)
    call grd_chksum2(bergs%grd, bergs%grd%rmean_calving_hflx, 'write mean calving_hflx')
    call write_data(filename, 'rmean_calving_hflx', bergs%grd%rmean_calving_hflx, bergs%grd%domain)
  endif
  contains

  function last_berg(berg)
  ! Arguments
  type(iceberg), pointer :: last_berg, berg
  ! Local variables

    last_berg=>berg
    do while (associated(last_berg%next))
      last_berg=>last_berg%next
    enddo

  end function last_berg

end subroutine write_restart

! ##############################################################################

subroutine read_restart_bergs_orig(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: Time
! Local variables
integer, dimension(:), allocatable :: found_restart_int
integer :: k, ierr, ncid, dimid, nbergs_in_file
integer :: lonid, latid,  uvelid, vvelid, ineid, jneid
integer :: axnid, aynid, uvel_oldid, vvel_oldid, bxnid, bynid
integer :: massid, thicknessid, widthid, lengthid
integer :: start_lonid, start_latid, start_yearid, iceberg_numid, start_dayid, start_massid
integer :: scaling_id, mass_of_bits_id, heat_density_id, static_bergid
logical :: lres, found_restart, multiPErestart
real :: lon0, lon1, lat0, lat1
character(len=33) :: filename, filename_base
type(icebergs_gridded), pointer :: grd
type(iceberg) :: localberg ! NOT a pointer but an actual local variable
integer :: stderrunit, iNg, jNg, i, j

  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>bergs%grd
  iNg=(grd%ieg-grd%isg+1) ! Total number of points globally in i direction
  jNg=(grd%jeg-grd%jsg+1) ! Total number of points globally in j direction

  ! Find a restart file
  multiPErestart=.false.

  ! Zero out nbergs_in_file
  nbergs_in_file = 0

  filename_base=trim(restart_input_dir)//'icebergs.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart, io_tile_id(1))

  ! Check if no restart found on any pe
  allocate(found_restart_int(mpp_npes()))
  if (found_restart .eqv. .true.) then
     k=1
  else
     k=0
  endif
  call mpp_gather((/k/),found_restart_int)
  if (sum(found_restart_int)==0.and.mpp_pe()==mpp_root_pe())&
       & write(*,'(a)') 'diamonds, read_restart_bergs: no restart file found'
  deallocate(found_restart_int)

  if (.not.found_restart) then

  multiPErestart=.true. ! This is to force sanity checking in a mulit-PE mode if no file was found on this PE

  elseif (found_restart) then ! if (.not.found_restart)
  ! only do the following if a file was found
  
  if (verbose.and.mpp_pe()==mpp_root_pe()) write(*,'(2a)') 'diamonds, read_restart_bergs: found restart file = ',filename

  ierr=nf_open(filename, NF_NOWRITE, ncid)
  if (ierr .ne. NF_NOERR) write(stderrunit,*) 'diamonds, read_restart_bergs: nf_open failed'

  ierr=nf_inq_unlimdim(ncid, dimid)
  if (ierr .ne. NF_NOERR) write(stderrunit,*) 'diamonds, read_restart_bergs: nf_inq_unlimdim failed'

  ierr=nf_inq_dimlen(ncid, dimid, nbergs_in_file)
  if (ierr .ne. NF_NOERR) write(stderrunit,*) 'diamonds, read_restart_bergs: nf_inq_dimlen failed'
 !write(stderrunit,*) 'diamonds, read_restart_bergs: nbergs in file =', nbergs_in_file

  lonid=inq_var(ncid, 'lon')
  latid=inq_var(ncid, 'lat')
  uvelid=inq_var(ncid, 'uvel')
  vvelid=inq_var(ncid, 'vvel')
  massid=inq_var(ncid, 'mass')
  axnid=inq_var(ncid, 'axn', unsafe=.true.)
  aynid=inq_var(ncid, 'ayn', unsafe=.true.)
  bxnid=inq_var(ncid, 'bxn', unsafe=.true.)
  bynid=inq_var(ncid, 'byn', unsafe=.true.)
  thicknessid=inq_var(ncid, 'thickness')
  widthid=inq_var(ncid, 'width')
  lengthid=inq_var(ncid, 'length')
  start_lonid=inq_var(ncid, 'start_lon')
  start_latid=inq_var(ncid, 'start_lat')
  start_yearid=inq_var(ncid, 'start_year')
  iceberg_numid=inq_var(ncid, 'icberg_num', unsafe=.true.)
  start_dayid=inq_var(ncid, 'start_day')
  start_massid=inq_var(ncid, 'start_mass')
  scaling_id=inq_var(ncid, 'mass_scaling')
  static_bergid=inq_var(ncid, 'static_berg', unsafe=.true.)
  mass_of_bits_id=inq_var(ncid, 'mass_of_bits', unsafe=.true.)
  heat_density_id=inq_var(ncid, 'heat_density', unsafe=.true.)
  ineid=inq_var(ncid, 'ine',unsafe=.true.)
  jneid=inq_var(ncid, 'jne',unsafe=.true.)

  ! Find approx outer bounds for tile
  lon0=minval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lon1=maxval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat0=minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat1=maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
 !write(stderrunit,'(a,i3,a,4f10.3)') 'diamonds, read_restart_bergs: (',mpp_pe(),') ',lon0,lon1,lat0,lat1
  do k=1, nbergs_in_file
   !write(stderrunit,*) 'diamonds, read_restart_bergs: reading berg ',k
    localberg%lon=get_double(ncid, lonid, k)
    localberg%lat=get_double(ncid, latid, k)
    if (ineid>0 .and. jneid>0 .and. .not. ignore_ij_restart) then ! read i,j position and avoid the "find" step
      localberg%ine=get_int(ncid, ineid, k)
      localberg%jne=get_int(ncid, jneid, k)
      if ( localberg%ine>=grd%isc .and. localberg%ine<=grd%iec .and. &
           localberg%jne>=grd%jsc .and.localberg%jne<=grd%jec ) then
        lres=.true.
      else
        lres=.false.
      endif
    else ! i,j are not available from the file so we search the grid to find out if we reside on this PE
      if (use_slow_find) then
        lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      else
        lres=find_cell_by_search(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      endif
    endif
    if (really_debug) then
      write(stderrunit,'(a,i8,a,2f9.4,a,i8)') 'diamonds, read_restart_bergs: berg ',k,' is at ',localberg%lon,localberg%lat,&
           & ' on PE ',mpp_pe()
      write(stderrunit,*) 'diamonds, read_restart_bergs: lres = ',lres
    endif
    if (lres) then ! true if we reside on this PE grid
      localberg%uvel=get_double(ncid, uvelid, k)
      localberg%vvel=get_double(ncid, vvelid, k)
      localberg%mass=get_double(ncid, massid, k)
      localberg%axn=get_real_from_file(ncid, axnid, k, value_if_not_in_file=0.)
      localberg%ayn=get_real_from_file(ncid, aynid, k, value_if_not_in_file=0.)
      localberg%bxn=get_real_from_file(ncid, bxnid, k, value_if_not_in_file=0.)
      localberg%byn=get_real_from_file(ncid, bynid, k, value_if_not_in_file=0.)
      localberg%thickness=get_double(ncid, thicknessid, k)
      localberg%width=get_double(ncid, widthid, k)
      localberg%length=get_double(ncid, lengthid, k)
      localberg%start_lon=get_double(ncid, start_lonid, k)
      localberg%start_lat=get_double(ncid, start_latid, k)
      localberg%start_year=get_int(ncid, start_yearid, k)
      if (iceberg_numid>0) then
        localberg%iceberg_num=get_int(ncid, iceberg_numid, k)
      else
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
      endif
      localberg%start_day=get_double(ncid, start_dayid, k)
      localberg%start_mass=get_double(ncid, start_massid, k)
      localberg%mass_scaling=get_double(ncid, scaling_id, k)
      localberg%halo_berg=0.
      localberg%static_berg=get_real_from_file(ncid, static_bergid, k, value_if_not_in_file=0.)
      localberg%mass_of_bits=get_real_from_file(ncid, mass_of_bits_id, k, value_if_not_in_file=0.)
      localberg%heat_density=get_real_from_file(ncid, heat_density_id, k, value_if_not_in_file=0.)
      if (really_debug) lres=is_point_in_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, explain=.true.)
      lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
     !call add_new_berg_to_list(bergs%first, localberg, quick=.true.)
      call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg)
      if (really_debug) call print_berg(stderrunit, bergs%list(localberg%ine,localberg%jne)%first, 'read_restart_bergs, add_new_berg_to_list')
    elseif (multiPErestart .and. io_tile_id(1) .lt. 0) then
      call error_mesg('diamonds, read_restart_bergs', 'berg in PE file was not on PE!', FATAL)
    endif
  enddo

  else ! if no restart file was read on this PE
    nbergs_in_file=0
  endif ! if (.not.found_restart)
  
  ! Sanity check
  k=count_bergs(bergs)
  if (verbose) write(*,'(2(a,i8))') 'diamonds, read_restart_bergs: # bergs =',k,' on PE',mpp_pe()
  if (multiPErestart)  then
      if (.NOT. is_io_tile_root_pe) nbergs_in_file=0 !If io_layout specified only tile root pes should count bergs
      call mpp_sum(nbergs_in_file) ! In case PE 0 didn't open a file
  endif
  call mpp_sum(k)
  bergs%nbergs_start=k
  if (mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,i8,a,i8,a)') 'diamonds, read_restart_bergs: there were',nbergs_in_file,' bergs in the restart file and', &
     k,' bergs have been read'
  endif

  if (k.ne.nbergs_in_file) call error_mesg('diamonds, read_restart_bergs', 'wrong number of bergs read!', FATAL) 

  if (.not. found_restart .and. bergs%nbergs_start==0 .and. generate_test_icebergs) call generate_bergs(bergs,Time)

  bergs%floating_mass_start=sum_mass(bergs)
  call mpp_sum( bergs%floating_mass_start )
  bergs%icebergs_mass_start=sum_mass(bergs,justbergs=.true.)
  call mpp_sum( bergs%icebergs_mass_start )
  bergs%bergy_mass_start=sum_mass(bergs,justbits=.true.)
  call mpp_sum( bergs%bergy_mass_start )
  if (mpp_pe().eq.mpp_root_pe().and.verbose) write(*,'(a)') 'diamonds, read_restart_bergs: completed'

contains

  real function get_real_from_file(ncid, varid, k, value_if_not_in_file)
  integer, intent(in) :: ncid, varid, k
  real, optional :: value_if_not_in_file

  if (varid<1.and.present(value_if_not_in_file)) then
    get_real_from_file=value_if_not_in_file
  else
    get_real_from_file=get_double(ncid, ncid, k)
  endif
  end function get_real_from_file
  
  subroutine generate_bergs(bergs,Time)
  ! Arguments
  type(icebergs), pointer :: bergs
  type(time_type), intent(in) :: Time
  ! Local variables
  integer :: i,j
  integer :: iNg, jNg  !Total number of points gloablly in i and j direction
  type(iceberg) :: localberg ! NOT a pointer but an actual local variable
  integer :: iyr, imon, iday, ihr, imin, isec

    ! For convenience
    grd=>bergs%grd
    iNg=(grd%ieg-grd%isg+1) ! Total number of points globally in i direction
    jNg=(grd%jeg-grd%jsg+1) ! Total number of points globally in j direction

    call get_date(Time, iyr, imon, iday, ihr, imin, isec)

    do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (grd%msk(i,j)>0. .and. abs(grd%latc(i,j))>60.) then
        localberg%xi=0.5
        localberg%yj=0.5
        localberg%ine=i
        localberg%jne=j
        localberg%lon=bilin(grd, grd%lon, i, j, localberg%xi, localberg%yj)
        localberg%lat=bilin(grd, grd%lat, i, j, localberg%xi, localberg%yj)
        localberg%lon_old=bilin(grd, grd%lon, i, j, localberg%xi, localberg%yj) !Alon
        localberg%lat_old=bilin(grd, grd%lat, i, j, localberg%xi, localberg%yj) !Alon
        localberg%mass=bergs%initial_mass(1)
        localberg%thickness=bergs%initial_thickness(1)
        localberg%width=bergs%initial_width(1)
        localberg%length=bergs%initial_length(1)
        localberg%start_lon=localberg%lon
        localberg%start_lat=localberg%lat
        localberg%start_year=iyr
        localberg%start_day=float(iday)+(float(ihr)+float(imin)/60.)/24.
        localberg%start_mass=localberg%mass
        localberg%mass_scaling=bergs%mass_scaling(1)
        localberg%mass_of_bits=0.
        localberg%halo_berg=0.
        localberg%static_berg=0.
        localberg%heat_density=0.
        localberg%axn=0. !Alon
        localberg%ayn=0. !Alon
        localberg%uvel_old=0. !Alon
        localberg%vvel_old=0. !Alon
        localberg%bxn=0. !Alon
        localberg%byn=0. !Alon
        
        !Berg A
        localberg%uvel=1.
        localberg%vvel=0.
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
        !Berg B
        localberg%uvel=-1.
        localberg%vvel=0.
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
        !Berg C
        localberg%uvel=0.
        localberg%vvel=1.
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
        !Berg D
        localberg%uvel=0.
        localberg%vvel=-1.
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
      endif
    enddo; enddo

    bergs%nbergs_start=count_bergs(bergs)
    call mpp_sum(bergs%nbergs_start)
    if (mpp_pe().eq.mpp_root_pe()) &
      write(*,'(a,i8,a)') 'diamonds, generate_bergs: ',bergs%nbergs_start,' were generated'

  end subroutine generate_bergs
  
end subroutine read_restart_bergs_orig

! ##############################################################################

subroutine read_restart_bergs(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: Time
! Local variables
integer :: k, siz(4), nbergs_in_file, nbergs_read
logical :: lres, found_restart, found
logical :: explain
logical :: multiPErestart  ! Not needed with new restart read; currently kept for compatibility
real :: lon0, lon1, lat0, lat1
real :: pos_is_good, pos_is_good_all_pe
character(len=33) :: filename, filename_base
type(icebergs_gridded), pointer :: grd
type(iceberg) :: localberg ! NOT a pointer but an actual local variable
integer :: stderrunit, iNg, jNg, i, j

real, allocatable, dimension(:) :: lon,          &
                                   lat,          &
                                   uvel,         &
                                   vvel,         &
                                   mass,         &
                                   axn,          &
                                   ayn,          &
                                   bxn,          &
                                   byn,          &
                                   thickness,    &
                                   width,        &
                                   length,       &
                                   start_lon,    &
                                   start_lat,    &
                                   start_day,    &
                                   start_mass,   &
                                   mass_scaling, &
                                   mass_of_bits, &
                                   static_berg,    &
                                   heat_density
integer, allocatable, dimension(:) :: ine,       &
                                      jne,       &
                                      iceberg_num,       &
                                      start_year

!integer, allocatable, dimension(:,:) :: iceberg_counter_grd

  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>bergs%grd
  iNg=(grd%ieg-grd%isg+1) ! Total number of points globally in i direction
  jNg=(grd%jeg-grd%jsg+1) ! Total number of points globally in j direction

  ! Zero out nbergs_in_file
  nbergs_in_file = 0

  filename_base=trim(restart_input_dir)//'icebergs.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart, io_tile_id(1))
  call error_mesg('read_restart_bergs_new', 'Using new icebergs restart read', NOTE)

  if (found_restart) then
     filename = filename_base
     call get_field_size(filename,'i',siz, field_found=found, domain=bergs%grd%domain)
     nbergs_in_file = siz(1)
  endif

  if(nbergs_in_file > 0) then
     allocate(lon(nbergs_in_file))
     allocate(lat(nbergs_in_file))
     allocate(uvel(nbergs_in_file))
     allocate(vvel(nbergs_in_file))
     allocate(mass(nbergs_in_file))
     allocate(axn(nbergs_in_file)) !Alon
     allocate(ayn(nbergs_in_file)) !Alon
     allocate(bxn(nbergs_in_file)) !Alon
     allocate(byn(nbergs_in_file)) !Alon
     allocate(thickness(nbergs_in_file))
     allocate(width(nbergs_in_file))
     allocate(length(nbergs_in_file))
     allocate(start_lon(nbergs_in_file))
     allocate(start_lat(nbergs_in_file))
     allocate(start_day(nbergs_in_file))
     allocate(start_mass(nbergs_in_file))
     allocate(mass_scaling(nbergs_in_file))
     allocate(mass_of_bits(nbergs_in_file))
     allocate(static_berg(nbergs_in_file))
     allocate(heat_density(nbergs_in_file))
     allocate(ine(nbergs_in_file))
     allocate(jne(nbergs_in_file))
     allocate(start_year(nbergs_in_file))
     allocate(iceberg_num(nbergs_in_file))
  endif

  if (found_restart) then
     call read_unlimited_axis(filename,'lon',lon,domain=grd%domain)
     call read_unlimited_axis(filename,'lat',lat,domain=grd%domain)
     call read_unlimited_axis(filename,'uvel',uvel,domain=grd%domain)
     call read_unlimited_axis(filename,'vvel',vvel,domain=grd%domain)
     call read_unlimited_axis(filename,'mass',mass,domain=grd%domain)
     call read_real_vector(filename,'axn',axn,grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'ayn',ayn,grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'bxn',bxn,grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'byn',byn,grd%domain,value_if_not_in_file=0.)
     call read_unlimited_axis(filename,'thickness',thickness,domain=grd%domain)
     call read_unlimited_axis(filename,'width',width,domain=grd%domain)
     call read_unlimited_axis(filename,'length',length,domain=grd%domain)
     call read_unlimited_axis(filename,'start_lon',start_lon,domain=grd%domain)
     call read_unlimited_axis(filename,'start_lat',start_lat,domain=grd%domain)
     call read_unlimited_axis(filename,'start_day',start_day,domain=grd%domain)
     call read_unlimited_axis(filename,'start_mass',start_mass,domain=grd%domain)
     call read_unlimited_axis(filename,'mass_scaling',mass_scaling,domain=grd%domain)
     call read_real_vector(filename,'mass_of_bits',mass_of_bits,domain=grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'heat_density',heat_density,domain=grd%domain,value_if_not_in_file=0.)
     call read_unlimited_axis(filename,'ine',ine,domain=grd%domain)
     call read_unlimited_axis(filename,'jne',jne,domain=grd%domain)
     call read_unlimited_axis(filename,'start_year',start_year,domain=grd%domain)
     call read_int_vector(filename,'iceberg_num',iceberg_num,grd%domain,value_if_not_in_file=-1)
     call read_real_vector(filename,'static_berg',static_berg,grd%domain,value_if_not_in_file=0.)
  endif

  ! Find approx outer bounds for tile
  lon0=minval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lon1=maxval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat0=minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat1=maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
     
  do k=1, nbergs_in_file
    localberg%lon=lon(k)
    localberg%lat=lat(k)
    if (.not. ignore_ij_restart) then ! read i,j position and avoid the "find" step
      localberg%ine=ine(k)
      localberg%jne=jne(k)
      if ( localberg%ine>=grd%isc .and. localberg%ine<=grd%iec .and. &
         localberg%jne>=grd%jsc .and.localberg%jne<=grd%jec ) then
        lres=.true.
      else
        lres=.false.
      endif
    else ! i,j are not available from the file so we search the grid to find out if we reside on this PE
      if (use_slow_find) then
        lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      else
        lres=find_cell_by_search(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      endif
    endif
    if (really_debug) then
      write(stderrunit,'(a,i8,a,2f9.4,a,i8)') 'diamonds, read_restart_bergs: berg ',k,' is at ',localberg%lon,localberg%lat,&
           & ' on PE ',mpp_pe()
      write(stderrunit,*) 'diamonds, read_restart_bergs: lres = ',lres
    endif
    if (lres) then ! true if we reside on this PE grid
      localberg%uvel=uvel(k)
      localberg%vvel=vvel(k)
      localberg%mass=mass(k)
      localberg%axn=axn(k) !Alon
      localberg%ayn=ayn(k) !Alon
      localberg%uvel_old=uvel(k) !Alon
      localberg%vvel_old=vvel(k) !Alon
      localberg%lon_old=lon(k) !Alon
      localberg%lat_old=lat(k) !Alon
      localberg%bxn=bxn(k) !Alon
      localberg%byn=byn(k) !Alon
      localberg%thickness=thickness(k)
      localberg%width=width(k)
      localberg%length=length(k)
      localberg%start_lon=start_lon(k)
      localberg%start_lat=start_lat(k)
      localberg%start_year=start_year(k)
      localberg%iceberg_num=iceberg_num(k)
      localberg%start_day=start_day(k)
      localberg%start_mass=start_mass(k)
      localberg%mass_scaling=mass_scaling(k)
      localberg%mass_of_bits=mass_of_bits(k)
      localberg%halo_berg=0.
      localberg%static_berg=static_berg(k)
      localberg%heat_density=heat_density(k)
      localberg%first_bond=>null()

      if (really_debug) lres=is_point_in_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, explain=.true.)
      lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
      !call add_new_berg_to_list(bergs%first, localberg, quick=.true.)

      if (bergs%grd%area(localberg%ine,localberg%jne) .ne. 0)  then
        if (iceberg_num(k)==-1) then ! If using an old_restart then iceberg_num needs to be generated
          localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(localberg%ine,localberg%jne))+(localberg%ine+(iNg*(localberg%jne-1)))
          grd%iceberg_counter_grd(localberg%ine,localberg%jne)=grd%iceberg_counter_grd(localberg%ine,localberg%jne)+1
        endif
        call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg)
      else
        if (mpp_pe().eq.mpp_root_pe()) then
          print * , 'Grounded iceberg: ', lat(k),lon(k), iceberg_num(k)
          call error_mesg('diamonds, read_restart_bergs', 'Iceberg not added because it is grounded', WARNING)
        endif
       endif

      if (really_debug) call print_berg(stderrunit, bergs%list(localberg%ine,localberg%jne)%first, 'read_restart_bergs, add_new_berg_to_list')
    elseif (multiPErestart .and. io_tile_id(1) .lt. 0) then
      call error_mesg('diamonds, read_restart_bergs', 'berg in PE file was not on PE!', FATAL)
    endif
  enddo
  
  if(nbergs_in_file > 0) then
    deallocate(              &
               lon,          &
               lat,          &
               uvel,         &
               vvel,         &
               mass,         &
               axn,          &
               ayn,          &
               bxn,          &
               byn,          &
               thickness,    &
               width,        &
               length,       &
               start_lon,    &
               start_lat,    &
               start_day,    &
               start_mass,   &
               mass_scaling, &
               mass_of_bits, &
               static_berg,  &
               heat_density )
    deallocate(              &
               ine,          &
               jne,          &
               iceberg_num,  &
               start_year )

    ! This block only works for IO_LAYOUT=1,1 or 0,0 but not for arbitrary layouts.
    ! I'm commenting this out until we find a way to implement the same sorts of checks
    ! that work for all i/o layouts. -AJA
    !Checking the total number of icebergs read from the restart file.
    !nbergs_read=count_bergs(bergs)
    !call mpp_sum(nbergs_read)
    !if (mpp_pe().eq.mpp_root_pe()) then
    !  write(*,'(a,i8,a,i8,a)') 'diamonds, read_restart_bergs: Number of Icebergs in restart file=',nbergs_in_file,' Number of Icebergs read=', nbergs_read
    !  if (nbergs_read .gt. nbergs_in_file) then
    !    call error_mesg('diamonds, read_restart_bergs', 'More icebergs read than exist in restart file.', FATAL)
    !  elseif (nbergs_read .lt. nbergs_in_file) then
    !    if (bergs%ignore_missing_restart_bergs) then
    !      call error_mesg('diamonds, read_restart_bergs', 'Some Icebergs from restart file were not found (ignore_missing flag is on)', WARNING)
    !    else
    !      call error_mesg('diamonds, read_restart_bergs', 'Some Icebergs from restart file were not found', FATAL)
    !    endif
    !  elseif (nbergs_read .eq. nbergs_in_file) then
    !    write(*,'(a,i8,a,i8,a)') 'diamonds, read_restart_bergs: Number of icebergs read (#',nbergs_read,') matches the number of icebergs in the file'
    !  endif
    !endif

  elseif(.not. found_restart .and. bergs%nbergs_start==0 .and. generate_test_icebergs) then
    call generate_bergs(bergs,Time)
  endif

  call check_for_duplicates_in_parallel(bergs)

  bergs%floating_mass_start=sum_mass(bergs)
  call mpp_sum( bergs%floating_mass_start )
  bergs%icebergs_mass_start=sum_mass(bergs,justbergs=.true.)
  call mpp_sum( bergs%icebergs_mass_start )
  bergs%bergy_mass_start=sum_mass(bergs,justbits=.true.)
  call mpp_sum( bergs%bergy_mass_start )
  if (mpp_pe().eq.mpp_root_pe().and.verbose) write(*,'(a)') 'diamonds, read_restart_bergs: completed'

contains

  subroutine read_real_vector(filename, varname, values, domain, value_if_not_in_file)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: varname
    real,             intent(out) :: values(:)
    type(domain2D),   intent(in)  :: domain
    real, optional,   intent(in)  :: value_if_not_in_file

    if (present(value_if_not_in_file).and..not.field_exist(filename, varname)) then
      values(:)=value_if_not_in_file
    else
      call read_unlimited_axis(filename,varname,values,domain=domain)
    endif
  end subroutine read_real_vector

  subroutine read_int_vector(filename, varname, values, domain, value_if_not_in_file)
    character(len=*),  intent(in)  :: filename
    character(len=*),  intent(in)  :: varname
    integer,           intent(out) :: values(:)
    type(domain2D),    intent(in)  :: domain
    integer, optional, intent(in)  :: value_if_not_in_file

    if (present(value_if_not_in_file).and..not.field_exist(filename, varname)) then
      values(:)=value_if_not_in_file
    else
      call read_unlimited_axis(filename,varname,values,domain=domain)
    endif
  end subroutine read_int_vector
  
  subroutine generate_bergs(bergs,Time)
  ! Arguments
  type(icebergs), pointer :: bergs
  type(time_type), intent(in) :: Time
  ! Local variables
  integer :: i,j
  integer :: iNg, jNg  !Total number of points gloablly in i and j direction
  type(iceberg) :: localberg ! NOT a pointer but an actual local variable
  integer :: iyr, imon, iday, ihr, imin, isec
  logical :: lres

    ! For convenience
    grd=>bergs%grd
    iNg=(grd%ieg-grd%isg+1) ! Total number of points globally in i direction
    jNg=(grd%jeg-grd%jsg+1) ! Total number of points globally in j direction

    call get_date(Time, iyr, imon, iday, ihr, imin, isec)

    do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (grd%msk(i,j)>0. .and. abs(grd%latc(i,j))>80.0) then
        if (max(grd%lat(i,j),grd%lat(i-1,j),grd%lat(i,j-1),grd%lat(i-1,j-1))>89.999) cycle ! Cannot use this at Pole cells
        localberg%xi=-999.
        localberg%yj=-999.
        localberg%ine=i
        localberg%jne=j
        localberg%lon=grd%lonc(i,j)
        localberg%lat=grd%latc(i,j)
        localberg%lon_old=localberg%lon
        localberg%lat_old=localberg%lat
        localberg%mass=bergs%initial_mass(1)
        localberg%thickness=bergs%initial_thickness(1)
        localberg%width=bergs%initial_width(1)
        localberg%length=bergs%initial_length(1)
        localberg%start_lon=localberg%lon
        localberg%start_lat=localberg%lat
        localberg%start_year=iyr
        localberg%start_day=float(iday)+(float(ihr)+float(imin)/60.)/24.
        localberg%start_mass=localberg%mass
        localberg%mass_scaling=bergs%mass_scaling(1)
        localberg%mass_of_bits=0.
        localberg%halo_berg=0.
        localberg%static_berg=0.
        localberg%heat_density=0.
        localberg%axn=0. !Alon
        localberg%ayn=0. !Alon
        localberg%uvel_old=0. !Alon
        localberg%vvel_old=0. !Alon
        localberg%bxn=0. !Alon
        localberg%byn=0. !Alon

        !Berg A
        call loc_set_berg_pos(grd, 0.9, 0.5, 1., 0., localberg)
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
        !Berg B
        call loc_set_berg_pos(grd, 0.1, 0.5, -1., 0., localberg)
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
        !Berg C
        call loc_set_berg_pos(grd, 0.5, 0.9, 0., 1., localberg)
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
        !Berg D
        call loc_set_berg_pos(grd, 0.5, 0.1, 0., -1., localberg)
        localberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i +(iNg*(j-1)))  ! unique number for each iceberg
        grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
      endif
    enddo; enddo

    bergs%nbergs_start=count_bergs(bergs)
    call mpp_sum(bergs%nbergs_start)
    if (mpp_pe().eq.mpp_root_pe()) &
      write(*,'(a,i8,a)') 'diamonds, generate_bergs: ',bergs%nbergs_start,' were generated'

  end subroutine generate_bergs

  subroutine loc_set_berg_pos(grd, xi, yj, uvel, vvel, berg)
    type(icebergs_gridded), pointer :: grd
    real, intent(in) :: xi, yj, uvel, vvel
    type(iceberg), intent(inout) :: berg
    integer :: i, j
    logical :: lres
    i = berg%ine ; j = berg%jne
    if (max(grd%lat(i,j),grd%lat(i-1,j),grd%lat(i,j-1),grd%lat(i-1,j-1))>89.999) then
      berg%lon=grd%lonc(i,j)
      berg%lat=grd%latc(i,j)
      berg%xi=0.5 ; berg%yj=0.5
    else
      berg%lon=bilin(grd, grd%lon, i, j, xi, yj)
      berg%lat=bilin(grd, grd%lat, i, j, xi, yj)
      berg%xi=xi ; berg%yj=yj
    endif
    berg%uvel=uvel ; berg%vvel=vvel
    berg%lon_old=berg%lon ; berg%lat_old=berg%lat
    berg%start_lon=berg%lon ; berg%start_lat=berg%lat
    lres=pos_within_cell(grd, berg%lon, berg%lat, berg%ine, berg%jne, berg%xi, berg%yj)
    if (.not. lres) then
      lres=pos_within_cell(grd, berg%lon, berg%lat, berg%ine, berg%jne, berg%xi, berg%yj, explain=.true.)
      write(0,*) lres, i, j, xi, yj, uvel, vvel
      write(0,*) lres, berg%ine, berg%jne, berg%xi, berg%yj
      write(0,*) 'bx=',berg%lon, 'gx=',grd%lon(i-1,j-1), grd%lon(i,j-1), grd%lon(i,j), grd%lon(i-1,j),'cx=', grd%lonc(i,j)
      write(0,*) 'by=',berg%lat, 'gy=',grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i,j), grd%lat(i-1,j),'cy=', grd%latc(i,j)
      stop 'generate_bergs, loc_set_berg_pos(): VERY FATAL!'
    endif
  end subroutine loc_set_berg_pos
  
end subroutine read_restart_bergs


! ##############################################################################
subroutine read_restart_bonds(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: Time
! Local variables
integer :: k, siz(4), nbonds_in_file
logical :: lres, found_restart, found
logical :: first_berg_found, second_berg_found
logical :: multiPErestart  ! Not needed with new restart read; currently kept for compatibility
real :: lon0, lon1, lat0, lat1
character(len=33) :: filename, filename_base
type(icebergs_gridded), pointer :: grd
type(iceberg) :: localberg ! NOT a pointer but an actual local variable
type(iceberg) , pointer :: this, first_berg, second_berg
type(bond) , pointer :: current_bond
integer :: stderrunit
integer :: number_first_bonds_matched !How many first bond bergs found on pe
integer :: number_second_bonds_matched !How many second bond bergs found on pe
integer :: number_perfect_bonds ! How many complete bonds formed
integer :: number_partial_bonds ! How many either complete/partial bonds formed.
integer :: number_perfect_bonds_with_first_on_pe ! How many bonds with first bond on the compuational domain
integer :: all_pe_number_perfect_bonds, all_pe_number_partial_bonds
integer :: all_pe_number_first_bonds_matched, all_pe_number_second_bonds_matched
integer :: all_pe_number_perfect_bonds_with_first_on_pe
integer :: ine, jne
logical :: search_data_domain
real :: berg_found, berg_found_all_pe
integer, allocatable, dimension(:) :: first_berg_num,   &
                                      other_berg_num,   &
                                      first_berg_jne,   &
                                      first_berg_ine,   &
                                      other_berg_jne,   &
                                      other_berg_ine
!integer, allocatable, dimension(:,:) :: iceberg_counter_grd

  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>bergs%grd

  ! Zero out nbergs_in_file
  nbonds_in_file = 0
  all_pe_number_perfect_bonds=0

  filename_base=trim(restart_input_dir)//'bonds_iceberg.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart, io_tile_id(1))
  call error_mesg('read_restart_bonds_bergs_new', 'Using icebergs bond restart read', NOTE)

  filename = filename_base
  call get_field_size(filename,'i',siz, field_found=found, domain=bergs%grd%domain)
  nbonds_in_file = siz(1)
  
    if (mpp_pe() .eq. mpp_root_pe()) then
      write(stderrunit,*)  'diamonds, bond read restart : ','Number of bonds in file',  nbonds_in_file
    endif

  if (nbonds_in_file .gt. 0) then

    allocate(first_berg_num(nbonds_in_file))
    allocate(other_berg_num(nbonds_in_file))
    allocate(first_berg_jne(nbonds_in_file))
    allocate(first_berg_ine(nbonds_in_file))
    allocate(other_berg_ine(nbonds_in_file))
    allocate(other_berg_jne(nbonds_in_file))


    call read_unlimited_axis(filename,'first_berg_num',first_berg_num,domain=grd%domain)
    call read_unlimited_axis(filename,'other_berg_num',other_berg_num,domain=grd%domain)
    call read_unlimited_axis(filename,'first_berg_jne',first_berg_jne,domain=grd%domain)
    call read_unlimited_axis(filename,'first_berg_ine',first_berg_ine,domain=grd%domain)
    call read_unlimited_axis(filename,'other_berg_jne',other_berg_jne,domain=grd%domain)
    call read_unlimited_axis(filename,'other_berg_ine',other_berg_ine,domain=grd%domain)

    number_first_bonds_matched=0
    number_second_bonds_matched=0
    number_perfect_bonds=0
    number_partial_bonds=0
    number_perfect_bonds_with_first_on_pe=0

    do k=1, nbonds_in_file
      

       ! If i,j in restart files are not good, then we find the berg position of the bond addresses manually:
       if (ignore_ij_restart) then 
         !Finding first iceberg in bond
         ine=999 ; jne=999 ; berg_found=0.0 ; search_data_domain=.true.
         call find_individual_iceberg(bergs,first_berg_num(k), ine, jne,berg_found,search_data_domain)
         berg_found_all_pe=berg_found
         call mpp_sum(berg_found_all_pe)
         if (berg_found_all_pe .gt. 0.5) then
             first_berg_ine(k)=ine
             first_berg_jne(k)=jne
         else
           print * , 'First bond berg not located: ', first_berg_num(k),berg_found, mpp_pe(),ine, jne
           call error_mesg('read_restart_bonds_bergs_new', 'First iceberg in bond not found on any pe', FATAL)
         endif
         !else

         !Finding other iceberg other iceberg
         ine=999 ; jne=999 ; berg_found=0.0 ; search_data_domain =.true.
         call find_individual_iceberg(bergs,other_berg_num(k), ine, jne, berg_found,search_data_domain)
         berg_found_all_pe=berg_found
         call mpp_sum(berg_found_all_pe)
         if (berg_found_all_pe .gt. 0.5) then
           !if (berg_found_all_pe .gt. 1.5) then
           !  call error_mesg('read_restart_bonds_bergs_new', 'Other iceberg bond found on more than one pe', FATAL)
           !else
             other_berg_ine(k)=ine
             other_berg_jne(k)=jne
           !endif
         else
          call error_mesg('read_restart_bonds_bergs_new', 'Other iceberg in bond not found on any pe', FATAL)
         endif
         if (berg_found_all_pe .lt. 0.5) then
                 print * , 'First bond berg not located: ', other_berg_num(k),berg_found, mpp_pe(),ine, jne
             call error_mesg('read_restart_bonds_bergs_new', 'First bond iceberg not located', FATAL)
         endif
       endif

      ! Decide whether the first iceberg is on the processeor
      if ( (first_berg_ine(k)>=grd%isd) .and. (first_berg_ine(k)<=grd%ied) .and. &
        (first_berg_jne(k)>=grd%jsd) .and. (first_berg_jne(k)<=grd%jed) ) then
        number_first_bonds_matched=number_first_bonds_matched+1
        
        ! Search for the first berg, which the bond belongs to
        first_berg_found=.false.
        first_berg=>null()
        this=>bergs%list(first_berg_ine(k),first_berg_jne(k))%first
        do while(associated(this))
          if (this%iceberg_num == first_berg_num(k)) then
            first_berg_found=.true.
            first_berg=>this
            !if (first_berg%halo_berg.gt.0.5) print *, 'bonding halo berg:', first_berg_num(k),  first_berg_ine(k),first_berg_jne(k) ,grd%isc, grd%iec, mpp_pe()
            this=>null()
          else  
            this=>this%next
          endif
        enddo
     

        ! Decide whether the second iceberg is on the processeor (data domain)
        second_berg_found=.false.
        !if ( other_berg_ine(k)>=grd%isc-1 .and. other_berg_ine(k)<=grd%iec+1 .and. &
        !  other_berg_jne(k)>=grd%jsc-1 .and.other_berg_jne(k)<=grd%jec+1 ) then
        if ( (other_berg_ine(k)>=grd%isd) .and. (other_berg_ine(k)<=grd%ied) .and. &
          (other_berg_jne(k)>=grd%jsd) .and.(other_berg_jne(k)<=grd%jed) ) then
          number_second_bonds_matched=number_second_bonds_matched+1

          ! Search for the second berg, which the bond belongs to
          second_berg=>null()
          this=>bergs%list(other_berg_ine(k),other_berg_jne(k))%first
          do while(associated(this))
            if (this%iceberg_num == other_berg_num(k)) then
              second_berg_found=.true.
              second_berg=>this
              this=>null()
            else  
              this=>this%next
            endif
          enddo
        endif
         
        if (first_berg_found) then
          number_partial_bonds=number_partial_bonds+1
          if (second_berg_found) then
            call form_a_bond(first_berg, other_berg_num(k), other_berg_ine(k), other_berg_jne(k),  second_berg)
            number_perfect_bonds=number_perfect_bonds+1
    
            !Counting number of bonds where the first bond is in the computational domain
            if ( (first_berg_ine(k)>=grd%isc) .and. (first_berg_ine(k)<=grd%iec) .and. &
              (first_berg_jne(k)>=grd%jsc) .and. (first_berg_jne(k)<=grd%jec) ) then
               number_perfect_bonds_with_first_on_pe=number_perfect_bonds_with_first_on_pe+1
        endif
     
          else
            !print *, 'Forming a bond of the second type', mpp_pe(), first_berg_num(k),  other_berg_num(k)
            !call form_a_bond(first_berg, other_berg_num(k),other_berg_ine(k),other_berg_jne(k))
          endif
        else
          write(stderrunit,*) 'diamonds, bond read restart : ','Not enough partial bonds formed', k, mpp_pe(), nbonds_in_file
          call error_mesg('read_restart_bonds_bergs_new', 'Failure with reading bonds: First bond not found on pe', FATAL)
        endif
      endif
    enddo

    !Analyse how many bonds were created and take appropriate action
    all_pe_number_perfect_bonds=number_perfect_bonds
    all_pe_number_perfect_bonds_with_first_on_pe=number_perfect_bonds_with_first_on_pe
    all_pe_number_partial_bonds=number_partial_bonds
    all_pe_number_first_bonds_matched=number_first_bonds_matched
    all_pe_number_second_bonds_matched=number_second_bonds_matched
    call mpp_sum(all_pe_number_perfect_bonds)
    call mpp_sum(all_pe_number_partial_bonds)
    call mpp_sum(all_pe_number_perfect_bonds_with_first_on_pe)

    if (all_pe_number_partial_bonds .lt. nbonds_in_file) then
      write(stderrunit,*) 'diamonds, bond read restart : ','Not enough partial bonds formed', all_pe_number_partial_bonds , nbonds_in_file
      call error_mesg('read_restart_bonds_bergs_new', 'Not enough partial bonds formed', FATAL)
    endif
    
    if (all_pe_number_perfect_bonds .lt. nbonds_in_file) then
      call mpp_sum(all_pe_number_first_bonds_matched)
      call mpp_sum(all_pe_number_second_bonds_matched)
      write(stderrunit,*)  'diamonds, bond read restart : ','Warning, some bonds are not fully formed',  all_pe_number_first_bonds_matched , nbonds_in_file
      write(stderrunit,*)  'diamonds, bond read restart : ','Number of first and second bonds matched:', all_pe_number_second_bonds_matched , nbonds_in_file
      call error_mesg('read_restart_bonds_bergs_new', 'Not enough perfect bonds formed', NOTE)
    endif

    if (all_pe_number_perfect_bonds_with_first_on_pe .ne. nbonds_in_file) then
      call mpp_sum(all_pe_number_first_bonds_matched)
      call mpp_sum(all_pe_number_second_bonds_matched)
      write(stderrunit,*)  'diamonds, bond read restart : ','Warning, # bonds with first bond on computational domain, does not match file',  all_pe_number_first_bonds_matched , nbonds_in_file
      write(stderrunit,*)  'diamonds, bond read restart : ','Computational bond, first second:', all_pe_number_second_bonds_matched , nbonds_in_file
      call error_mesg('read_restart_bonds_bergs_new', 'Computational perfect bonds do not match those in file', NOTE)
    endif

    deallocate(               &
            first_berg_num,   &
            other_berg_num,  &
            first_berg_ine,   &
            first_berg_jne,   &
            other_berg_ine,  &
            other_berg_jne )
  endif
    
  if (mpp_pe() .eq. mpp_root_pe()) then
    write(stderrunit,*)  'diamonds, bond read restart : ','Number of bonds (including halos)',  all_pe_number_perfect_bonds
    write(stderrunit,*)  'diamonds, bond read restart : ','Number of true bonds created',  all_pe_number_perfect_bonds_with_first_on_pe
  endif

end subroutine read_restart_bonds

! ##############################################################################

subroutine read_restart_calving(bergs)
use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream
! Arguments
type(icebergs), pointer :: bergs
! Local variables
integer :: k,i,j
character(len=37) :: filename, actual_filename
type(icebergs_gridded), pointer :: grd
real, allocatable, dimension(:,:) :: randnum
type(randomNumberStream) :: rns

  ! For convenience
  grd=>bergs%grd

  ! Read stored ice
  filename=trim(restart_input_dir)//'calving.res.nc'
  if (file_exist(filename)) then
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') &
     'diamonds, read_restart_calving: reading ',filename
    call read_data(filename, 'stored_ice', grd%stored_ice, grd%domain)
    if (field_exist(filename, 'stored_heat')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'diamonds, read_restart_calving: reading stored_heat from restart file.'
      call read_data(filename, 'stored_heat', grd%stored_heat, grd%domain)
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: stored_heat WAS NOT FOUND in the file. Setting to 0.'
      grd%stored_heat(:,:)=0.
    endif
    if (field_exist(filename, 'rmean_calving')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'diamonds, read_restart_calving: reading rmean_calving from restart file.'
      call read_data(filename, 'rmean_calving', grd%rmean_calving, grd%domain)
      grd%rmean_calving_initialized=.true.
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: rmean_calving WAS NOT FOUND in the file. Setting to 0.'
    endif
    if (field_exist(filename, 'rmean_calving_hflx')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'diamonds, read_restart_calving: reading rmean_calving_hflx from restart file.'
      call read_data(filename, 'rmean_calving_hflx', grd%rmean_calving_hflx, grd%domain)
      grd%rmean_calving_hflx_initialized=.true.
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: rmean_calving_hflx WAS NOT FOUND in the file. Setting to 0.'
    endif
    if (field_exist(filename, 'iceberg_counter_grd')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'diamonds, read_restart_calving: reading iceberg_counter_grd from restart file.'
      call read_data(filename, 'iceberg_counter_grd', grd%iceberg_counter_grd, grd%domain) 
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: iceberg_counter_grd WAS NOT FOUND in the file. Setting to 0.'
      grd%iceberg_counter_grd(:,:)=1
    endif
    bergs%restarted=.true.
  else
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: initializing stored ice to random numbers'
    if ( make_calving_reproduce ) then
       allocate(randnum(1,nclasses))
       do j=grd%jsc, grd%jec
          do i=grd%isc, grd%iec
             rns=initializeRandomNumberStream(i+10000*j)
             call getRandomNumbers(rns,randnum(1,:))
             do k=1, nclasses
                grd%stored_ice(i,j,k)=randnum(1,k) * grd%msk(i,j) * bergs%initial_mass(k) * bergs%mass_scaling(k)
             end do
          end do
       end do
    else
       allocate(randnum(grd%jsc:grd%jec,nclasses))
       do i=grd%isc, grd%iec
          rns = initializeRandomNumberStream(i)
          call getRandomNumbers(rns,randnum)
          do k=1, nclasses
             grd%stored_ice(i,grd%jsc:grd%jec,k) = randnum(:,k) * grd%msk(i,grd%jsc:grd%jec) * &
                  & bergs%initial_mass(k) * bergs%mass_scaling(k)
          end do
       end do
    end if
    deallocate(randnum)
  endif

  call grd_chksum3(bergs%grd, bergs%grd%stored_ice, 'read_restart_calving, stored_ice')
  call grd_chksum2(bergs%grd, bergs%grd%stored_heat, 'read_restart_calving, stored_heat')
  call grd_chksum2(bergs%grd, bergs%grd%rmean_calving, 'read_restart_calving, rmean_calving')
  call grd_chksum2(bergs%grd, bergs%grd%rmean_calving_hflx, 'read_restart_calving, rmean_calving_hflx')

  bergs%stored_start=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
  bergs%rmean_calving_start=sum( grd%rmean_calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%rmean_calving_hflx_start=sum( grd%rmean_calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( bergs%stored_start )
  bergs%stored_heat_start=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( bergs%stored_heat_start )
  bergs%floating_heat_start=sum_heat(bergs)
  call mpp_sum( bergs%floating_heat_start )

end subroutine read_restart_calving

! ##############################################################################

subroutine read_ocean_depth(grd)
! Arguments
! Local variables
character(len=37) :: filename 
type(icebergs_gridded), pointer :: grd

  ! Read stored ice
  filename=trim(restart_input_dir)//'topog.nc'
  if (file_exist(filename)) then
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') &
     'diamonds, read_ocean_depth: reading ',filename
    if (field_exist(filename, 'depth')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'diamonds, read_restart_calving: reading stored_heat from restart file.'
      call read_data(filename, 'depth', grd%ocean_depth, grd%domain)
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: stored_heat WAS NOT FOUND in the file. Setting to 0.'
      !grd%ocean_depth(:,:)=0.
    endif
  else
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_ocean_depth: Ocean depth file (topog.nc) not present)'
  endif

  !call grd_chksum2(bergs%grd, bergs%grd%ocean_depth, 'read_ocean_depth, ocean_depth')
end subroutine read_ocean_depth

! ##############################################################################

subroutine write_trajectory(trajectory, save_short_traj)
! Arguments
type(xyt), pointer :: trajectory
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, uvelid, vvelid, iceberg_numid
integer :: uoid, void, uiid, viid, uaid, vaid, sshxid, sshyid, sstid, sssid
integer :: cnid, hiid
integer :: mid, did, wid, lid, mbid, hdid
character(len=37) :: filename
character(len=7) :: pe_name
type(xyt), pointer :: this, next
integer :: stderrunit
logical, intent(in) :: save_short_traj
!I/O vars
type(xyt), pointer :: traj4io=>null()
integer :: ntrajs_sent_io,ntrajs_rcvd_io
integer :: from_pe,np
type(buffer), pointer :: obuffer_io=>null(), ibuffer_io=>null()
logical :: io_is_in_append_mode

  ! Get the stderr unit number
  stderrunit=stderr()
  traj4io=>null()
  obuffer_io=>null()
  ibuffer_io=>null()

  !Assemble the list of trajectories from all pes in this I/O tile
  call mpp_clock_begin(clock_trp)

  !First add the trajs on the io_tile_root_pe (if any) to the I/O list
  if(is_io_tile_root_pe .OR. force_all_pes_traj ) then
     if(associated(trajectory)) then
        this=>trajectory
        do while (associated(this))
           call append_posn(traj4io, this)
           this=>this%next
        enddo
        trajectory => null()
     endif
  endif

  if(.NOT. force_all_pes_traj ) then

  !Now gather and append the bergs from all pes in the io_tile to the list on corresponding io_tile_root_pe
  ntrajs_sent_io =0
  ntrajs_rcvd_io =0 

  if(is_io_tile_root_pe) then
     !Receive trajs from all pes in this I/O tile !FRAGILE!SCARY!
     do np=2,size(io_tile_pelist) ! Note: np starts from 2 to exclude self
        from_pe=io_tile_pelist(np)
        call mpp_recv(ntrajs_rcvd_io, glen=1, from_pe=from_pe, tag=COMM_TAG_11)
        if (ntrajs_rcvd_io .gt. 0) then
           call increase_ibuffer(ibuffer_io, ntrajs_rcvd_io,buffer_width_traj)
           call mpp_recv(ibuffer_io%data, ntrajs_rcvd_io*buffer_width_traj,from_pe=from_pe, tag=COMM_TAG_12)
           do i=1, ntrajs_rcvd_io
              call unpack_traj_from_buffer2(traj4io, ibuffer_io, i, save_short_traj)
           enddo
       endif
     enddo
  else
     ! Pack and send trajectories to the root PE for this I/O tile
     do while (associated(trajectory))
       ntrajs_sent_io = ntrajs_sent_io +1
       call pack_traj_into_buffer2(trajectory, obuffer_io, ntrajs_sent_io, save_short_traj)
       this => trajectory ! Need to keep pointer in order to free up the links memory
       trajectory => trajectory%next ! This will eventually result in trajectory => null()
       deallocate(this) ! Delete the link from memory
     enddo
        
     call mpp_send(ntrajs_sent_io, plen=1, to_pe=io_tile_root_pe, tag=COMM_TAG_11)
     if (ntrajs_sent_io .gt. 0) then
        call mpp_send(obuffer_io%data, ntrajs_sent_io*buffer_width_traj, to_pe=io_tile_root_pe, tag=COMM_TAG_12)
     endif
  endif

  endif !.NOT. force_all_pes_traj

  call mpp_clock_end(clock_trp)

  !Now start writing in the io_tile_root_pe if there are any bergs in the I/O list
  call mpp_clock_begin(clock_trw)

  if((force_all_pes_traj .OR. is_io_tile_root_pe) .AND. associated(traj4io)) then
 
    call get_instance_filename("iceberg_trajectories.nc", filename)
    if(io_tile_id(1) .ge. 0 .AND. .NOT. force_all_pes_traj) then !io_tile_root_pes write
       if(io_npes .gt. 1) then !attach tile_id  to filename only if there is more than one I/O pe
          if (io_tile_id(1)<10000) then
             write(filename,'(A,".",I4.4)') trim(filename), io_tile_id(1) 
          else
             write(filename,'(A,".",I6.6)') trim(filename), io_tile_id(1) 
          endif
       endif
    else !All pes write, attach pe# to filename
       if (mpp_npes()<10000) then
          write(filename,'(A,".",I4.4)') trim(filename), mpp_pe() 
       else
          write(filename,'(A,".",I6.6)') trim(filename), mpp_pe() 
       endif
    endif

    io_is_in_append_mode = .false.
    iret = nf_create(filename, NF_NOCLOBBER, ncid)
    if (iret .ne. NF_NOERR) then
      iret = nf_open(filename, NF_WRITE, ncid)
      io_is_in_append_mode = .true.
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_open failed'
    endif
    if (verbose) then
      if (io_is_in_append_mode) then
        write(*,'(2a)') 'diamonds, write_trajectory: appending to ',filename
      else
        write(*,'(2a)') 'diamonds, write_trajectory: creating ',filename
      endif
    endif

    if (io_is_in_append_mode) then
      iret = nf_inq_dimid(ncid, 'i', i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_inq_dimid i failed'
      lonid = inq_varid(ncid, 'lon')
      latid = inq_varid(ncid, 'lat')
      yearid = inq_varid(ncid, 'year')
      dayid = inq_varid(ncid, 'day')
      iceberg_numid = inq_varid(ncid, 'iceberg_num')
      if (.not.save_short_traj) then
        uvelid = inq_varid(ncid, 'uvel')
        vvelid = inq_varid(ncid, 'vvel')
        uoid = inq_varid(ncid, 'uo')
        void = inq_varid(ncid, 'vo')
        uiid = inq_varid(ncid, 'ui')
        viid = inq_varid(ncid, 'vi')
        uaid = inq_varid(ncid, 'ua')
        vaid = inq_varid(ncid, 'va')
        mid = inq_varid(ncid, 'mass')
        mbid = inq_varid(ncid, 'mass_of_bits')
        hdid = inq_varid(ncid, 'heat_density')
        did = inq_varid(ncid, 'thickness')
        wid = inq_varid(ncid, 'width')
        lid = inq_varid(ncid, 'length')
        sshxid = inq_varid(ncid, 'ssh_x')
        sshyid = inq_varid(ncid, 'ssh_y')
        sstid = inq_varid(ncid, 'sst')
        sssid = inq_varid(ncid, 'sss')
        cnid = inq_varid(ncid, 'cn')
        hiid = inq_varid(ncid, 'hi')
      endif
    else
      ! Dimensions
      iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_def_dim i failed'

      ! Variables
      lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
      latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
      yearid = def_var(ncid, 'year', NF_INT, i_dim)
      dayid = def_var(ncid, 'day', NF_DOUBLE, i_dim)
      iceberg_numid = def_var(ncid, 'iceberg_num', NF_INT, i_dim)
      if (.not. save_short_traj) then
        uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
        vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
        uoid = def_var(ncid, 'uo', NF_DOUBLE, i_dim)
        void = def_var(ncid, 'vo', NF_DOUBLE, i_dim)
        uiid = def_var(ncid, 'ui', NF_DOUBLE, i_dim)
        viid = def_var(ncid, 'vi', NF_DOUBLE, i_dim)
        uaid = def_var(ncid, 'ua', NF_DOUBLE, i_dim)
        vaid = def_var(ncid, 'va', NF_DOUBLE, i_dim)
        mid = def_var(ncid, 'mass', NF_DOUBLE, i_dim)
        mbid = def_var(ncid, 'mass_of_bits', NF_DOUBLE, i_dim)
        hdid = def_var(ncid, 'heat_density', NF_DOUBLE, i_dim)
        did = def_var(ncid, 'thickness', NF_DOUBLE, i_dim)
        wid = def_var(ncid, 'width', NF_DOUBLE, i_dim)
        lid = def_var(ncid, 'length', NF_DOUBLE, i_dim)
        sshxid = def_var(ncid, 'ssh_x', NF_DOUBLE, i_dim)
        sshyid = def_var(ncid, 'ssh_y', NF_DOUBLE, i_dim)
        sstid = def_var(ncid, 'sst', NF_DOUBLE, i_dim)
        sssid = def_var(ncid, 'sss', NF_DOUBLE, i_dim)
        cnid = def_var(ncid, 'cn', NF_DOUBLE, i_dim)
        hiid = def_var(ncid, 'hi', NF_DOUBLE, i_dim)
      endif

      ! Attributes
      iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_major_version', NF_INT, 1, 0)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_minor_version', NF_INT, 1, 1)
      call put_att(ncid, lonid, 'long_name', 'longitude')
      call put_att(ncid, lonid, 'units', 'degrees_E')
      call put_att(ncid, latid, 'long_name', 'latitude')
      call put_att(ncid, latid, 'units', 'degrees_N')
      call put_att(ncid, yearid, 'long_name', 'year')
      call put_att(ncid, yearid, 'units', 'years')
      call put_att(ncid, dayid, 'long_name', 'year day')
      call put_att(ncid, dayid, 'units', 'days')
      call put_att(ncid, iceberg_numid, 'long_name', 'iceberg id number')
      call put_att(ncid, iceberg_numid, 'units', 'dimensionless')
      
      if (.not. save_short_traj) then
        call put_att(ncid, uvelid, 'long_name', 'zonal spped')
        call put_att(ncid, uvelid, 'units', 'm/s')
        call put_att(ncid, vvelid, 'long_name', 'meridional spped')
        call put_att(ncid, vvelid, 'units', 'm/s')
        call put_att(ncid, uoid, 'long_name', 'ocean zonal spped')
        call put_att(ncid, uoid, 'units', 'm/s')
        call put_att(ncid, void, 'long_name', 'ocean meridional spped')
        call put_att(ncid, void, 'units', 'm/s')
        call put_att(ncid, uiid, 'long_name', 'ice zonal spped')
        call put_att(ncid, uiid, 'units', 'm/s')
        call put_att(ncid, viid, 'long_name', 'ice meridional spped')
        call put_att(ncid, viid, 'units', 'm/s')
        call put_att(ncid, uaid, 'long_name', 'atmos zonal spped')
        call put_att(ncid, uaid, 'units', 'm/s')
        call put_att(ncid, vaid, 'long_name', 'atmos meridional spped')
        call put_att(ncid, vaid, 'units', 'm/s')
        call put_att(ncid, mid, 'long_name', 'mass')
        call put_att(ncid, mid, 'units', 'kg')
        call put_att(ncid, mbid, 'long_name', 'mass_of_bits')
        call put_att(ncid, mbid, 'units', 'kg')
        call put_att(ncid, hdid, 'long_name', 'heat_density')
        call put_att(ncid, hdid, 'units', 'J/kg')
        call put_att(ncid, did, 'long_name', 'thickness')
        call put_att(ncid, did, 'units', 'm')
        call put_att(ncid, wid, 'long_name', 'width')
        call put_att(ncid, wid, 'units', 'm')
        call put_att(ncid, lid, 'long_name', 'length')
        call put_att(ncid, lid, 'units', 'm')
        call put_att(ncid, sshxid, 'long_name', 'sea surface height gradient_x')
        call put_att(ncid, sshxid, 'units', 'non-dim')
        call put_att(ncid, sshyid, 'long_name', 'sea surface height gradient_y')
        call put_att(ncid, sshyid, 'units', 'non-dim')
        call put_att(ncid, sstid, 'long_name', 'sea surface temperature')
        call put_att(ncid, sstid, 'units', 'degrees_C')
        call put_att(ncid, sssid, 'long_name', 'sea surface salinity')
        call put_att(ncid, sssid, 'units', 'psu')
        call put_att(ncid, cnid, 'long_name', 'sea ice concentration')
        call put_att(ncid, cnid, 'units', 'none')
        call put_att(ncid, hiid, 'long_name', 'sea ice thickness')
        call put_att(ncid, hiid, 'units', 'm')
      endif
    endif

    ! End define mode
    iret = nf_enddef(ncid)
         
    ! Write variables
    this=>traj4io
    if (io_is_in_append_mode) then
      iret = nf_inq_dimlen(ncid, i_dim, i)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_inq_dimlen i failed'
    else
      i = 0
    endif
    do while (associated(this))
      i=i+1
      call put_double(ncid, lonid, i, this%lon)
      call put_double(ncid, latid, i, this%lat)
      call put_int(ncid, yearid, i, this%year)
      call put_double(ncid, dayid, i, this%day)
      call put_int(ncid, iceberg_numid, i, this%iceberg_num)
      if (.not. save_short_traj) then
        call put_double(ncid, uvelid, i, this%uvel)
        call put_double(ncid, vvelid, i, this%vvel)
        call put_double(ncid, uoid, i, this%uo)
        call put_double(ncid, void, i, this%vo)
        call put_double(ncid, uiid, i, this%ui)
        call put_double(ncid, viid, i, this%vi)
        call put_double(ncid, uaid, i, this%ua)
        call put_double(ncid, vaid, i, this%va)
        call put_double(ncid, mid, i, this%mass)
        call put_double(ncid, hdid, i, this%heat_density)
        call put_double(ncid, did, i, this%thickness)
        call put_double(ncid, wid, i, this%width)
        call put_double(ncid, lid, i, this%length)
        call put_double(ncid, sshxid, i, this%ssh_x)
        call put_double(ncid, sshyid, i, this%ssh_y)
        call put_double(ncid, sstid, i, this%sst)
        call put_double(ncid, sssid, i, this%sss)
        call put_double(ncid, cnid, i, this%cn)
        call put_double(ncid, hiid, i, this%hi)
      endif
      next=>this%next
      deallocate(this)
      this=>next
    enddo

    ! Finish up
    iret = nf_close(ncid)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_close failed',mpp_pe(),filename

  endif !(is_io_tile_root_pe .AND. associated(traj4io))
  call mpp_clock_end(clock_trw)

end subroutine write_trajectory


! ##############################################################################

integer function inq_var(ncid, var, unsafe)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
logical, optional, intent(in) :: unsafe
! Local variables
integer :: iret
integer :: stderrunit
logical :: unsafely=.false.

if(present(unsafe)) unsafely=unsafe
  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_inq_varid(ncid, var, inq_var)
  if (iret .ne. NF_NOERR) then
    if (.not. unsafely) then
      write(stderrunit,*) 'diamonds, inq_var: nf_inq_varid ',var,' failed'
      call error_mesg('diamonds, inq_var', 'netcdf function returned a failure!', FATAL)
    else
      inq_var=-1
    endif
  endif

end function inq_var

! ##############################################################################

integer function def_var(ncid, var, ntype, idim)
! Arguments
integer, intent(in) :: ncid, ntype, idim
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_def_var(ncid, var, ntype, 1, idim, def_var)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, def_var: nf_def_var failed for ',trim(var)
    call error_mesg('diamonds, def_var', 'netcdf function returned a failure!', FATAL)
  endif

end function def_var

! ##############################################################################

integer function inq_varid(ncid, var)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_inq_varid(ncid, var, inq_varid)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, inq_varid: nf_inq_varid failed for ',trim(var)
    call error_mesg('diamonds, inq_varid', 'netcdf function returned a failure!', FATAL)
  endif

end function inq_varid

! ##############################################################################

subroutine put_att(ncid, id, att, attval)
! Arguments
integer, intent(in) :: ncid, id
character(len=*), intent(in) :: att, attval
! Local variables
integer :: vallen, iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  vallen=len_trim(attval)
  iret = nf_put_att_text(ncid, id, att, vallen, attval)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_att: nf_put_att_text failed adding', &
      trim(att),' = ',trim(attval)
    call error_mesg('diamonds, put_att', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_att

! ##############################################################################

real function get_double(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_double(ncid, id, i, get_double)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, get_double: nf_get_var1_double failed reading'
    call error_mesg('diamonds, get_double', 'netcdf function returned a failure!', FATAL)
  endif

end function get_double

! ##############################################################################

integer function get_int(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_int(ncid, id, i, get_int)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, get_int: nf_get_var1_int failed reading'
    call error_mesg('diamonds, get_int', 'netcdf function returned a failure!', FATAL)
  endif

end function get_int

! ##############################################################################

subroutine put_double(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i
real, intent(in) :: val
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_double(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_double: nf_put_vara_double failed writing'
    call error_mesg('diamonds, put_double', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_double

! ##############################################################################

subroutine put_int(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i, val
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_int(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_int: nf_put_vara_int failed writing'
    call error_mesg('diamonds, put_int', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_int


! ##############################################################################

logical function find_restart_file(filename, actual_file, multiPErestart, tile_id)
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: actual_file
  logical, intent(out) :: multiPErestart
  integer, intent(in) :: tile_id

  character(len=6) :: pe_name

  find_restart_file = .false.

  ! If running as ensemble, add the ensemble id string to the filename
  call get_instance_filename(filename, actual_file)
    
  ! Prefer combined restart files.
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) return
    
  ! Uncombined restart
  if(tile_id .ge. 0) then
    write(actual_file,'(A,".",I4.4)') trim(actual_file), tile_id
  else
  if (mpp_npes()>10000) then
     write(pe_name,'(a,i6.6)' )'.', mpp_pe()    
  else
     write(pe_name,'(a,i4.4)' )'.', mpp_pe()    
  endif
  actual_file=trim(actual_file)//trim(pe_name)
  endif
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) then
     multiPErestart=.true.
     return
  endif

  ! No file found, Reset all return parameters
  find_restart_file=.false.
  actual_file = ''
  multiPErestart=.false.

end function find_restart_file


end module
