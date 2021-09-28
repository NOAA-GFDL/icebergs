!> Handles reading/writing of restart files and trajectory-based diagnostic files
module ice_bergs_io

! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.

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

use ice_bergs_framework, only: icebergs_gridded, xyt, bond_xyt, iceberg, icebergs, buffer, bond
use ice_bergs_framework, only: pack_traj_into_buffer2,unpack_traj_from_buffer2
use ice_bergs_framework, only: find_cell,find_cell_by_search,count_bergs,is_point_in_cell,pos_within_cell,append_posn
use ice_bergs_framework, only: count_bonds, form_a_bond, find_individual_iceberg
use ice_bergs_framework, only: push_posn
use ice_bergs_framework, only: add_new_berg_to_list,destroy_iceberg
use ice_bergs_framework, only: increase_ibuffer,grd_chksum2,grd_chksum3
use ice_bergs_framework, only: sum_mass,sum_heat,bilin
!params !Niki: write a subroutine to get these
use ice_bergs_framework, only: nclasses, buffer_width, buffer_width_traj, buffer_width_bond_traj
use ice_bergs_framework, only: verbose, really_debug, debug, restart_input_dir,make_calving_reproduce
use ice_bergs_framework, only: ignore_ij_restart, use_slow_find,generate_test_icebergs,print_berg
use ice_bergs_framework, only: force_all_pes_traj
use ice_bergs_framework, only: check_for_duplicates_in_parallel
use ice_bergs_framework, only: split_id, id_from_2_ints, generate_id
! for MTS/DEM/fracture:
use ice_bergs_framework, only: mts,save_bond_traj
use ice_bergs_framework, only: push_bond_posn, append_bond_posn
use ice_bergs_framework, only: pack_bond_traj_into_buffer2,unpack_bond_traj_from_buffer2
use ice_bergs_framework, only: monitor_energy, use_damage, dem, fracture_criterion, iceberg_bonds_on

implicit none ; private

include 'netcdf.inc'

public ice_bergs_io_init
public read_restart_bergs, write_restart, write_trajectory, write_bond_trajectory
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
integer :: clock_btrw,clock_btrp !bond trajectories

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

!> Initialize parallel i/o
subroutine ice_bergs_io_init(bergs, io_layout)
type(icebergs), pointer :: bergs !< Icebergs container
integer, intent(in) :: io_layout(2) !< Decomposition of i/o processors

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
  clock_btrw=mpp_clock_id( 'Iceberg-bonds-traj write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  clock_btrp=mpp_clock_id( 'Iceberg-bonds-traj prepare', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

end subroutine ice_bergs_io_init

!> Write an iceberg restart file
subroutine write_restart(bergs, time_stamp)
! Arguments
type(icebergs), pointer :: bergs !< Icebergs container
character(len=*), intent(in), optional :: time_stamp !< Timestamp for restart file

! Local variables
type(bond), pointer :: current_bond
integer :: i,j,id
character(len=35) :: filename
character(len=35) :: filename_bonds
type(iceberg), pointer :: this=>NULL()
integer :: stderrunit
!I/O vars
type(restart_file_type) :: bergs_restart
type(restart_file_type) :: bergs_bond_restart
type(restart_file_type) :: calving_restart
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
                                   fl_k,         &
                                   start_lon,    &
                                   start_lat,    &
                                   start_day,    &
                                   start_mass,   &
                                   mass_scaling, &
                                   mass_of_bits, &
                                   mass_of_fl_bits, &
                                   mass_of_fl_bergy_bits, &
                                   static_berg,  &
                                   heat_density, &
                                   axn_fast,     &
                                   ayn_fast,     &
                                   bxn_fast,     &
                                   byn_fast,     &
                                   damage,       &
                                   ang_vel,      &
                                   ang_accel,    &
                                   rot,          &
                                   lon_prev,     &
                                   lat_prev,     &
                                   accum_bond_rotation, &
                                   tangd1,       &
                                   tangd2,       &
                                   nstress,      &
                                   sstress,      &
                                   rotation,     &
                                   rel_rotation, &
                                   n_frac_var,   &
                                   n_strain_rate,&
                                   spring_pe

integer, allocatable, dimension(:) :: ine,              &
                                      jne,              &
                                      id_cnt,           &
                                      id_ij,            &
                                      start_year,       &
                                      first_id_cnt,     &
                                      other_id_cnt,     &
                                      first_id_ij,      &
                                      other_id_ij,      &
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
   allocate(fl_k(nbergs))
   allocate(start_lon(nbergs))
   allocate(start_lat(nbergs))
   allocate(start_day(nbergs))
   allocate(start_mass(nbergs))
   allocate(mass_scaling(nbergs))
   allocate(mass_of_bits(nbergs))
   allocate(mass_of_fl_bits(nbergs))
   allocate(mass_of_fl_bergy_bits(nbergs))
   allocate(heat_density(nbergs))
   allocate(static_berg(nbergs))

   allocate(ine(nbergs))
   allocate(jne(nbergs))
   allocate(start_year(nbergs))
   allocate(id_cnt(nbergs))
   allocate(id_ij(nbergs))

   if (mts) then
     allocate(axn_fast(nbergs))
     allocate(ayn_fast(nbergs))
     allocate(bxn_fast(nbergs))
     allocate(byn_fast(nbergs))
   endif
   if (dem) then
     allocate(ang_vel(nbergs))
     allocate(ang_accel(nbergs))
     allocate(rot(nbergs))
   elseif (fracture_criterion.ne.'none') then
     allocate(lon_prev(nbergs))
     allocate(lat_prev(nbergs))
     allocate(accum_bond_rotation(nbergs))
   endif

  filename = trim("icebergs.res.nc")
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
    id = register_restart_field(bergs_restart,filename,'bxn',bxn,longname='implicit zonal acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'byn',byn,longname='implicit meridional acceleration',units='m/s^2')
  endif
  id = register_restart_field(bergs_restart,filename,'ine',ine,longname='i index',units='none')
  id = register_restart_field(bergs_restart,filename,'jne',jne,longname='j index',units='none')
  id = register_restart_field(bergs_restart,filename,'thickness',thickness,longname='thickness',units='m')
  id = register_restart_field(bergs_restart,filename,'width',width,longname='width',units='m')
  id = register_restart_field(bergs_restart,filename,'length',length,longname='length',units='m')
  id = register_restart_field(bergs_restart,filename,'fl_k',fl_k,longname='footloose calving k',units='m')
  id = register_restart_field(bergs_restart,filename,'start_lon',start_lon, &
                                            longname='longitude of calving location',units='degrees_E')
  id = register_restart_field(bergs_restart,filename,'start_lat',start_lat, &
                                            longname='latitude of calving location',units='degrees_N')
  id = register_restart_field(bergs_restart,filename,'start_year',start_year, &
                                            longname='calendar year of calving event', units='years')
  id = register_restart_field(bergs_restart,filename,'id_cnt',id_cnt, &
                                            longname='counter component of iceberg id', units='dimensionless')
  id = register_restart_field(bergs_restart,filename,'id_ij',id_ij, &
                                            longname='position component of iceberg id', units='dimensionless')
  id = register_restart_field(bergs_restart,filename,'start_day',start_day, &
                                            longname='year day of calving event',units='days')
  id = register_restart_field(bergs_restart,filename,'start_mass',start_mass, &
                                            longname='initial mass of calving berg',units='kg')
  id = register_restart_field(bergs_restart,filename,'mass_scaling',mass_scaling, &
                                            longname='scaling factor for mass of calving berg',units='none')
  id = register_restart_field(bergs_restart,filename,'mass_of_bits',mass_of_bits, &
                                            longname='mass of bergy bits',units='kg')
  id = register_restart_field(bergs_restart,filename,'mass_of_fl_bits',mass_of_fl_bits, &
                                            longname='mass of footloose bits',units='kg')
  id = register_restart_field(bergs_restart,filename,'mass_of_fl_bergy_bits',mass_of_fl_bergy_bits, &
                                            longname='mass of bergy bits associated with footloose bits',units='kg')
  id = register_restart_field(bergs_restart,filename,'heat_density',heat_density, &
                                            longname='heat density',units='J/kg')

  if (mts) then
    id = register_restart_field(bergs_restart,filename,'axn_fast',axn_fast,&
      longname='sub-step explicit zonal acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'ayn_fast',ayn_fast,&
      longname='sub-step explicit meridional acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'bxn_fast',bxn_fast,&
      longname='sub-step implicit zonal acceleration',units='m/s^2')
    id = register_restart_field(bergs_restart,filename,'byn_fast',byn_fast,&
      longname='sub-step implicit meridional acceleration',units='m/s^2')
  endif

  if (dem) then
    id = register_restart_field(bergs_restart,filename,'ang_vel',ang_vel,&
      longname='dem angular velocity',units='rad/s')
    id = register_restart_field(bergs_restart,filename,'ang_accel',ang_accel,&
      longname='dem angular acceleration',units='rad/s^2')
    id = register_restart_field(bergs_restart,filename,'rot',rot,&
      longname='dem accumulated rotation',units='rad')
  elseif (fracture_criterion.ne.'none') then
    id = register_restart_field(bergs_restart,filename,'lon_prev',lon_prev,longname='previous longitude',units='degrees_E')
    id = register_restart_field(bergs_restart,filename,'lat_prev',lat_prev,longname='previous latitude',units='degrees_N')
    id = register_restart_field(bergs_restart,filename,'accum_bond_rotation',accum_bond_rotation,&
      longname='accumulated bond rotation',units='rad')
  endif

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
      fl_k(i) = this%fl_k
      start_lon(i) = this%start_lon; start_lat(i) = this%start_lat
      start_year(i) = this%start_year; start_day(i) = this%start_day
      start_mass(i) = this%start_mass; mass_scaling(i) = this%mass_scaling
      static_berg(i) = this%static_berg
      call split_id(this%id, id_cnt(i), id_ij(i))
      mass_of_bits(i) = this%mass_of_bits; mass_of_fl_bits(i) = this%mass_of_fl_bits
      mass_of_fl_bergy_bits(i) = this%mass_of_fl_bergy_bits; heat_density(i) = this%heat_density
      if (mts) then
        axn_fast(i) = this%axn_fast
        ayn_fast(i) = this%ayn_fast
        bxn_fast(i) = this%bxn_fast
        byn_fast(i) = this%byn_fast
      endif
      if (dem) then
        ang_vel(i) = this%ang_vel
        ang_accel(i) = this%ang_accel
        rot(i) = this%rot
      elseif (fracture_criterion.ne.'none') then
        lon_prev(i) = this%lon_prev
        lat_prev(i) = this%lat_prev
        accum_bond_rotation(i) = this%accum_bond_rotation
      endif
      this=>this%next
    enddo
  enddo ; enddo


  call save_restart(bergs_restart, time_stamp)
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
             fl_k,         &
             start_lon,    &
             start_lat,    &
             start_day,    &
             start_mass,   &
             mass_scaling, &
             mass_of_bits, &
             mass_of_fl_bits, &
             mass_of_fl_bergy_bits, &
             static_berg,  &
             heat_density)
  if (mts) then
    deallocate(            &
             axn_fast,     &
             ayn_fast,     &
             bxn_fast,     &
             byn_fast)
  endif
  if (dem) then
    deallocate(            &
             ang_vel,      &
             ang_accel,    &
             rot)
  elseif (fracture_criterion.ne.'none') then
    deallocate(            &
             lon_prev,     &
             lat_prev,     &
             accum_bond_rotation)
  endif

  deallocate(           &
             ine,       &
             jne,       &
             id_cnt,    &
             id_ij,     &
             start_year )

  call nullify_domain()

!########## Creating bond restart file ######################

  !Allocating restart memory for bond related variables.
  nbonds=0
  if (bergs%iceberg_bonds_on) then
    check_bond_quality=.true.
    call count_bonds(bergs, nbonds,check_bond_quality)

  allocate(first_id_cnt(nbonds))
  allocate(other_id_cnt(nbonds))
  allocate(first_id_ij(nbonds))
  allocate(other_id_ij(nbonds))
  allocate(first_berg_ine(nbonds))
  allocate(first_berg_jne(nbonds))
  allocate(other_berg_ine(nbonds))
  allocate(other_berg_jne(nbonds))

  if (use_damage) allocate(damage(nbonds))
  if (dem) then
    allocate(tangd1(nbonds))
    allocate(tangd2(nbonds))
    allocate(nstress(nbonds))
    allocate(sstress(nbonds))
  elseif (fracture_criterion.ne.'none') then
    allocate(rotation(nbonds))
    allocate(rel_rotation(nbonds))
    allocate(n_frac_var(nbonds))
    if (fracture_criterion.eq.'strain_rate') then
      allocate(n_strain_rate(nbonds))
    elseif (fracture_criterion.eq.'energy') then
      allocate(spring_pe(nbonds))
    endif
  endif

  call get_instance_filename("bonds_iceberg.res.nc", filename_bonds)
  call set_domain(bergs%grd%domain)
  call register_restart_axis(bergs_bond_restart,filename,'i',nbonds)
  call set_meta_global(bergs_bond_restart,'file_format_major_version',ival=(/file_format_major_version/))
  call set_meta_global(bergs_bond_restart,'file_format_minor_version',ival=(/file_format_minor_version/))
  call set_meta_global(bergs_bond_restart,'time_axis',ival=(/0/))

  !Now start writing in the io_tile_root_pe if there are any bergs in the I/O list

  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_berg_ine',first_berg_ine,longname='iceberg ine of first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_berg_jne',first_berg_jne,longname='iceberg jne of first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_id_cnt',first_id_cnt,longname='counter component of iceberg id first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'first_id_ij',first_id_ij,longname='position component of iceberg id first berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_berg_ine',other_berg_ine,longname='iceberg ine of second berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_berg_jne',other_berg_jne,longname='iceberg jne of second berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_id_cnt',other_id_cnt,longname='counter component of iceberg id second berg in bond',units='dimensionless')
  id = register_restart_field(bergs_bond_restart,filename_bonds,'other_id_ij',other_id_ij,longname='position component of iceberg id second berg in bond',units='dimensionless')

  if (use_damage) then
    id = register_restart_field(bergs_bond_restart,filename_bonds,'damage',damage,&
      longname='damage',units='dimensionless')
  endif
  if (dem) then
    id = register_restart_field(bergs_bond_restart,filename_bonds,'tangd1',tangd1,&
      longname='zonal tangential displacement',units='m')
    id = register_restart_field(bergs_bond_restart,filename_bonds,'tangd2',tangd2,&
      longname='meridional tangential displacement',units='m')
    id = register_restart_field(bergs_bond_restart,filename_bonds,'nstress',nstress,&
      longname='normal stress',units='Pa')
    id = register_restart_field(bergs_bond_restart,filename_bonds,'sstress',sstress,&
      longname='shear stress',units='Pa')
  elseif (fracture_criterion.ne.'none') then
    id = register_restart_field(bergs_bond_restart,filename_bonds,'rotation',rotation,&
      longname='rotation',units='rad')
    id = register_restart_field(bergs_bond_restart,filename_bonds,'rel_rotation',rel_rotation,&
      longname='relative rotation',units='rad')
    id = register_restart_field(bergs_bond_restart,filename_bonds,'n_frac_var',n_frac_var,&
      longname='normal fracture variable',units='dimensionless/s^-1/Pa/J')
    if (fracture_criterion.eq.'strain_rate') then
      id = register_restart_field(bergs_bond_restart,filename_bonds,'n_strain_rate',n_strain_rate,&
        longname='normal strain-rate',units='1/seconds')
    elseif (fracture_criterion.eq.'energy') then
      id = register_restart_field(bergs_bond_restart,filename_bonds,'spring_pe',spring_pe,&
        longname='spring potential energy',units='J')
    endif
  endif

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
        call split_id( this%id, first_id_cnt(i), first_id_ij(i) )
        call split_id( current_bond%other_id, other_id_cnt(i), other_id_ij(i) )
        other_berg_ine(i)=current_bond%other_berg%ine
        other_berg_jne(i)=current_bond%other_berg%jne

        if (use_damage) damage(i)=current_bond%damage
        if (dem) then
          tangd1(i)=current_bond%tangd1
          tangd2(i)=current_bond%tangd2
          nstress(i)=current_bond%nstress
          sstress(i)=current_bond%sstress
        elseif (fracture_criterion.ne.'none') then
          rotation(i)=current_bond%rotation
          rel_rotation(i)=current_bond%rel_rotation
          n_frac_var(i)=current_bond%n_frac_var
          if (fracture_criterion.eq.'strain_rate') then
            n_strain_rate(i)=current_bond%n_strain_rate
          elseif (fracture_criterion.eq.'energy') then
            spring_pe(i)=current_bond%spring_pe
          endif
        endif

        current_bond=>current_bond%next_bond
      enddo !End of loop over bonds
      this=>this%next
    enddo!End of loop over bergs
  enddo; enddo !End of loop over grid

  call save_restart(bergs_bond_restart, time_stamp)
  call free_restart_type(bergs_bond_restart)


  deallocate(first_id_cnt,          &
             other_id_cnt,          &
             first_id_ij,           &
             other_id_ij,           &
             first_berg_ine,        &
             first_berg_jne,        &
             other_berg_ine,        &
             other_berg_jne)
  if (dem) then
    deallocate(                     &
             tangd1,                &
             tangd2,                &
             nstress,               &
             sstress)
  elseif (fracture_criterion.ne.'none') then
    deallocate(                     &
             rotation,              &
             rel_rotation,          &
             n_frac_var)
    if (fracture_criterion.eq.'strain_rate') then
      deallocate(n_strain_rate)
    elseif (fracture_criterion.eq.'energy') then
      deallocate(spring_pe)
    endif
  endif

  call nullify_domain()
  endif

  ! Write stored ice
  filename='calving.res.nc'
  if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(stderrunit,'(2a)') 'KID, write_restart: writing ',filename
  call grd_chksum3(bergs%grd, bergs%grd%stored_ice, 'write stored_ice')
  call grd_chksum2(bergs%grd, bergs%grd%stored_heat, 'write stored_heat')
  if (bergs%tau_calving>0.) then
    call grd_chksum2(bergs%grd, bergs%grd%rmean_calving, 'write mean calving')
    call grd_chksum2(bergs%grd, bergs%grd%rmean_calving_hflx, 'write mean calving_hflx')
  endif

  call set_domain(bergs%grd%domain)
  id = register_restart_field(calving_restart,filename,'stored_ice',bergs%grd%stored_ice,longname='STORED_ICE',units='none')
  id = register_restart_field(calving_restart,filename,'stored_heat',bergs%grd%stored_heat,longname='STORED_HEAT',units='none')
  id = register_restart_field(calving_restart,filename,'iceberg_counter_grd',bergs%grd%iceberg_counter_grd,longname='ICEBERG_COUNTER_GRD',units='none')
  if (bergs%tau_calving>0.) then
    id = register_restart_field(calving_restart,filename,'rmean_calving',bergs%grd%rmean_calving,longname='RMEAN_CALVING',units='none')
    id = register_restart_field(calving_restart,filename,'rmean_calving_hflx',bergs%grd%rmean_calving_hflx,longname='RMEAN_CALVING_HFLX',units='none')
  endif

  call save_restart(calving_restart, time_stamp)
  call free_restart_type(calving_restart)

end subroutine write_restart

!> Find the last berg in a linked list.
function last_berg(berg)
! Arguments
type(iceberg), pointer :: berg !< Pointer to an iceberg
type(iceberg), pointer :: last_berg
! Local variables

  last_berg=>berg
  do while (associated(last_berg%next))
    last_berg=>last_berg%next
  enddo

end function last_berg

!> Read a real value from a file and optionally return a default value if variable is missing
real function get_real_from_file(ncid, varid, k, value_if_not_in_file)
integer, intent(in) :: ncid !< file id
integer, intent(in) :: varid !< variable id
integer, intent(in) :: k !< index
real, optional :: value_if_not_in_file !< default value, if variable missing in file

if (varid<1.and.present(value_if_not_in_file)) then
  get_real_from_file=value_if_not_in_file
else
  get_real_from_file=get_double(ncid, ncid, k)
endif
end function get_real_from_file

!> Read an iceberg restart file
subroutine read_restart_bergs(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs !< Icebergs container
type(time_type), intent(in) :: Time !< Model time
! Local variables
integer :: k, siz(4), nbergs_in_file, nbergs_read
logical :: lres, found_restart, found, replace_iceberg_num
logical :: explain
logical :: multiPErestart  ! Not needed with new restart read; currently kept for compatibility
real :: lon0, lon1, lat0, lat1
real :: pos_is_good, pos_is_good_all_pe
character(len=33) :: filename, filename_base
type(icebergs_gridded), pointer :: grd
type(iceberg) :: localberg ! NOT a pointer but an actual local variable
integer :: stderrunit, i, j, cnt, ij

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
                                   fl_k,         &
                                   start_lon,    &
                                   start_lat,    &
                                   start_day,    &
                                   start_mass,   &
                                   mass_scaling, &
                                   mass_of_bits, &
                                   mass_of_fl_bits, &
                                   mass_of_fl_bergy_bits, &
                                   static_berg,  &
                                   heat_density, &
                                   axn_fast,     &
                                   ayn_fast,     &
                                   bxn_fast,     &
                                   byn_fast,     &
                                   ang_vel,      &
                                   ang_accel,    &
                                   rot,          &
                                   lon_prev,     &
                                   lat_prev,     &
                                   accum_bond_rotation

integer, allocatable, dimension(:) :: ine,        &
                                      jne,        &
                                      iceberg_num,&
                                      id_cnt,     &
                                      id_ij,      &
                                      start_year

!integer, allocatable, dimension(:,:) :: iceberg_counter_grd

  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>bergs%grd

  ! Zero out nbergs_in_file
  nbergs_in_file = 0

  filename_base=trim(restart_input_dir)//'icebergs.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart, io_tile_id(1))

  if (found_restart) then
     filename = filename_base
     call get_field_size(filename,'i',siz, field_found=found, domain=bergs%grd%domain)
     nbergs_in_file = siz(1)

     replace_iceberg_num = field_exist(filename, 'iceberg_num') ! True if using 32-bit iceberg_num in restart file
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
     allocate(fl_k(nbergs_in_file))
     allocate(start_lon(nbergs_in_file))
     allocate(start_lat(nbergs_in_file))
     allocate(start_day(nbergs_in_file))
     allocate(start_mass(nbergs_in_file))
     allocate(mass_scaling(nbergs_in_file))
     allocate(mass_of_bits(nbergs_in_file))
     allocate(mass_of_fl_bits(nbergs_in_file))
     allocate(mass_of_fl_bergy_bits(nbergs_in_file))
     allocate(static_berg(nbergs_in_file))
     allocate(heat_density(nbergs_in_file))
     allocate(ine(nbergs_in_file))
     allocate(jne(nbergs_in_file))
     allocate(start_year(nbergs_in_file))
     if (replace_iceberg_num) then
       call error_mesg('read_restart_bergs', "Calculating new iceberg ID's", WARNING)
       allocate(iceberg_num(nbergs_in_file))
     else
       allocate(id_cnt(nbergs_in_file))
       allocate(id_ij(nbergs_in_file))
     endif

      if (iceberg_bonds_on) then
        allocate(localberg%n_bonds)
        localberg%n_bonds=0.
      end if

      if (monitor_energy) then
        allocate(localberg%Ee,localberg%Ed,localberg%Eext,localberg%Ee_contact,localberg%Ed_contact,localberg%Efrac,&
          localberg%Ee_temp,localberg%Ed_temp,localberg%Eext_temp,localberg%Ee_contact_temp,localberg%Ed_contact_temp)
        localberg%Ee=0.0; localberg%Ed=0.0; localberg%Eext=0.0; localberg%Ee_contact=0.0; localberg%Ed_contact=0.0
        localberg%Efrac=0.0; localberg%Ee_temp=0.0; localberg%Ed_temp=0.0; localberg%Eext_temp=0.0
        localberg%Ee_contact_temp=0.0; localberg%Ed_contact_temp=0.0
      endif

      if (mts) then
        allocate(localberg%axn_fast)
        allocate(localberg%ayn_fast)
        allocate(localberg%bxn_fast)
        allocate(localberg%byn_fast)
        allocate(localberg%conglom_id)
        localberg%conglom_id=0.

       allocate(axn_fast(nbergs_in_file))
       allocate(ayn_fast(nbergs_in_file))
       allocate(bxn_fast(nbergs_in_file))
       allocate(byn_fast(nbergs_in_file))
     endif
     if (dem) then
       allocate(localberg%ang_vel)
       allocate(localberg%ang_accel)
       allocate(localberg%rot)
       allocate(ang_vel(nbergs_in_file))
       allocate(ang_accel(nbergs_in_file))
       allocate(rot(nbergs_in_file))
     elseif (fracture_criterion.ne.'none') then
       allocate(lon_prev(nbergs_in_file))
       allocate(lat_prev(nbergs_in_file))
       allocate(accum_bond_rotation(nbergs_in_file))
     endif

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
     call read_real_vector(filename,'fl_k',fl_k,grd%domain,value_if_not_in_file=0.)
     call read_unlimited_axis(filename,'start_lon',start_lon,domain=grd%domain)
     call read_unlimited_axis(filename,'start_lat',start_lat,domain=grd%domain)
     call read_unlimited_axis(filename,'start_day',start_day,domain=grd%domain)
     call read_unlimited_axis(filename,'start_mass',start_mass,domain=grd%domain)
     call read_unlimited_axis(filename,'mass_scaling',mass_scaling,domain=grd%domain)
     call read_real_vector(filename,'mass_of_bits',mass_of_bits,domain=grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'mass_of_fl_bits',mass_of_fl_bits,domain=grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'mass_of_fl_bergy_bits',mass_of_fl_bergy_bits,domain=grd%domain,value_if_not_in_file=0.)
     call read_real_vector(filename,'heat_density',heat_density,domain=grd%domain,value_if_not_in_file=0.)
     call read_unlimited_axis(filename,'ine',ine,domain=grd%domain)
     call read_unlimited_axis(filename,'jne',jne,domain=grd%domain)
     call read_unlimited_axis(filename,'start_year',start_year,domain=grd%domain)
     if (replace_iceberg_num) then
       call read_int_vector(filename,'iceberg_num',iceberg_num,grd%domain,value_if_not_in_file=-1)
     else
       call read_int_vector(filename,'id_cnt',id_cnt,grd%domain)
       call read_int_vector(filename,'id_ij',id_ij,grd%domain)
     endif
     call read_real_vector(filename,'static_berg',static_berg,grd%domain,value_if_not_in_file=0.)

     if (mts) then
       call read_real_vector(filename,'axn_fast',axn_fast,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'ayn_fast',ayn_fast,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'bxn_fast',bxn_fast,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'byn_fast',byn_fast,grd%domain,value_if_not_in_file=0.)
     endif

     if (dem) then
       call read_real_vector(filename,'ang_vel'  ,ang_vel  ,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'ang_accel',ang_accel,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'rot'      ,rot      ,grd%domain,value_if_not_in_file=0.)
     elseif (fracture_criterion.ne.'none') then
       call read_real_vector(filename,'lon_prev',lon_prev,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'lat_prev',lat_prev,grd%domain,value_if_not_in_file=0.)
       call read_real_vector(filename,'accum_bond_rotation',accum_bond_rotation,grd%domain,value_if_not_in_file=0.)
     endif
  elseif (bergs%require_restart) then
     stop 'read_restart_bergs, RESTART NOT FOUND!'
  endif

   if (all(lon_prev.eq.0)) lon_prev=lon
   if (all(lat_prev.eq.0)) lat_prev=lat

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
      write(stderrunit,'(a,i8,a,2f9.4,a,i8)') 'KID, read_restart_bergs: berg ',k,' is at ',localberg%lon,localberg%lat,&
           & ' on PE ',mpp_pe()
      write(stderrunit,*) 'KID, read_restart_bergs: lres = ',lres
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
      localberg%lon_prev=lon(k)
      localberg%lat_prev=lat(k)
      localberg%bxn=bxn(k) !Alon
      localberg%byn=byn(k) !Alon
      localberg%thickness=thickness(k)
      localberg%width=width(k)
      localberg%length=length(k)
      localberg%fl_k=fl_k(k)
      localberg%start_lon=start_lon(k)
      localberg%start_lat=start_lat(k)
      localberg%start_year=start_year(k)
      if (replace_iceberg_num) then
        localberg%id = generate_id(grd, localberg%ine, localberg%jne)
      else
        localberg%id=id_from_2_ints(id_cnt(k), id_ij(k))
      endif
      localberg%start_day=start_day(k)
      localberg%start_mass=start_mass(k)
      localberg%mass_scaling=mass_scaling(k)
      localberg%mass_of_bits=mass_of_bits(k)
      localberg%mass_of_fl_bits=mass_of_fl_bits(k)
      localberg%mass_of_fl_bergy_bits=mass_of_fl_bergy_bits(k)
      localberg%halo_berg=0.
      localberg%static_berg=static_berg(k)
      localberg%heat_density=heat_density(k)
      localberg%first_bond=>null()

     if (mts) then
       localberg%axn_fast=axn_fast(k)
       localberg%ayn_fast=ayn_fast(k)
       localberg%bxn_fast=bxn_fast(k)
       localberg%byn_fast=byn_fast(k)
     endif

     if (dem) then
       localberg%ang_vel  =ang_vel(k)
       localberg%ang_accel=ang_accel(k)
       localberg%rot      =rot(k)
     elseif (fracture_criterion.ne.'none') then
       localberg%lon_prev=lon_prev(k)
       localberg%lat_prev=lat_prev(k)
       localberg%accum_bond_rotation=accum_bond_rotation(k)
     endif

      if (really_debug) lres=is_point_in_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, explain=.true.)
      lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
      !call add_new_berg_to_list(bergs%first, localberg)

      if (bergs%grd%area(localberg%ine,localberg%jne) .ne. 0)  then
        call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg)
      else
        if (mpp_pe().eq.mpp_root_pe()) then
          print * , 'Grounded iceberg: ', lat(k),lon(k), localberg%id
          call error_mesg('KID, read_restart_bergs', 'Iceberg not added because it is grounded', WARNING)
        endif
       endif

      if (really_debug) call print_berg(stderrunit, bergs%list(localberg%ine,localberg%jne)%first, 'read_restart_bergs, add_new_berg_to_list')
    elseif (multiPErestart .and. io_tile_id(1) .lt. 0) then
      call error_mesg('KID, read_restart_bergs', 'berg in PE file was not on PE!', FATAL)
    endif
  enddo

  if (found_restart) then
    deallocate(lon,          &
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
               fl_k,         &
               start_lon,    &
               start_lat,    &
               start_day,    &
               start_mass,   &
               mass_scaling, &
               mass_of_bits, &
               mass_of_fl_bits, &
               mass_of_fl_bergy_bits, &
               static_berg,  &
               heat_density, &
               ine,          &
               jne,          &
               start_year )

    if (mts) then
      deallocate(            &
               axn_fast,     &
               ayn_fast,     &
               bxn_fast,     &
               byn_fast)
    endif
    if (dem) then
      deallocate(            &
               ang_vel,      &
               ang_accel,    &
               rot)
    elseif (fracture_criterion.ne.'none') then
      deallocate(            &
               lon_prev,     &
               lat_prev,     &
               accum_bond_rotation)
    endif

    if (replace_iceberg_num) then
      deallocate(iceberg_num)
    else
      deallocate(id_cnt)
      deallocate(id_ij)
    endif

    ! This block only works for IO_LAYOUT=1,1 or 0,0 but not for arbitrary layouts.
    ! I'm commenting this out until we find a way to implement the same sorts of checks
    ! that work for all i/o layouts. -AJA
    !Checking the total number of icebergs read from the restart file.
    !nbergs_read=count_bergs(bergs)
    !call mpp_sum(nbergs_read)
    !if (mpp_pe().eq.mpp_root_pe()) then
    !  write(*,'(a,i8,a,i8,a)') 'KID, read_restart_bergs: Number of Icebergs in restart file=',nbergs_in_file,' Number of Icebergs read=', nbergs_read
    !  if (nbergs_read .gt. nbergs_in_file) then
    !    call error_mesg('KID, read_restart_bergs', 'More icebergs read than exist in restart file.', FATAL)
    !  elseif (nbergs_read .lt. nbergs_in_file) then
    !    if (bergs%ignore_missing_restart_bergs) then
    !      call error_mesg('KID, read_restart_bergs', 'Some Icebergs from restart file were not found (ignore_missing flag is on)', WARNING)
    !    else
    !      call error_mesg('KID, read_restart_bergs', 'Some Icebergs from restart file were not found', FATAL)
    !    endif
    !  elseif (nbergs_read .eq. nbergs_in_file) then
    !    write(*,'(a,i8,a,i8,a)') 'KID, read_restart_bergs: Number of icebergs read (#',nbergs_read,') matches the number of icebergs in the file'
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
  bergs%fl_bits_mass_start=sum_mass(bergs,justflbits=.true.)
  call mpp_sum( bergs%fl_bits_mass_start )
  if (mpp_pe().eq.mpp_root_pe().and.verbose) write(*,'(a)') 'KID, read_restart_bergs: completed'

end subroutine read_restart_bergs

!> Read a vector of reals from file and use a default value if variable is missing
subroutine read_real_vector(filename, varname, values, domain, value_if_not_in_file)
  character(len=*), intent(in)  :: filename !< Name of file to read from
  character(len=*), intent(in)  :: varname !< Name of variable to read
  real,             intent(out) :: values(:) !< Returned vector of reals
  type(domain2D),   intent(in)  :: domain !< Parallel decomposition
  real, optional,   intent(in)  :: value_if_not_in_file !< Value to use if variable is not in file

  if (present(value_if_not_in_file).and..not.field_exist(filename, varname)) then
    values(:)=value_if_not_in_file
  else
    call read_unlimited_axis(filename,varname,values,domain=domain)
  endif
end subroutine read_real_vector

!> Read a vector of integers from file and use a default value if variable is missing
subroutine read_int_vector(filename, varname, values, domain, value_if_not_in_file)
  character(len=*),  intent(in)  :: filename !< Name of file to read from
  character(len=*),  intent(in)  :: varname !< Name of variable to read
  integer,           intent(out) :: values(:) !< Returned vector of integers
  type(domain2D),    intent(in)  :: domain !< Parallel decomposition
  integer, optional, intent(in)  :: value_if_not_in_file !< Value to use if variable is not in file

  if (present(value_if_not_in_file).and..not.field_exist(filename, varname)) then
    values(:)=value_if_not_in_file
  else
    call read_unlimited_axis(filename,varname,values,domain=domain)
  endif
end subroutine read_int_vector

!> Generate bergs for the purpose of debugging
subroutine generate_bergs(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs !< Icebergs container
type(time_type), intent(in) :: Time !< Model time
! Local variables
type(icebergs_gridded), pointer :: grd
integer :: i,j
type(iceberg) :: localberg ! NOT a pointer but an actual local variable
integer :: iyr, imon, iday, ihr, imin, isec
logical :: lres

  ! For convenience
  grd=>bergs%grd

  call get_date(Time, iyr, imon, iday, ihr, imin, isec)

  if (mts) then
    allocate(localberg%axn_fast)
    allocate(localberg%ayn_fast)
    allocate(localberg%bxn_fast)
    allocate(localberg%byn_fast)
  endif
  if (dem) then
    allocate(localberg%ang_vel)
    allocate(localberg%ang_accel)
    allocate(localberg%rot)
  elseif (fracture_criterion.ne.'none') then
    allocate(localberg%accum_bond_rotation)
  endif

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
      localberg%mass=bergs%initial_mass_s(1)
      localberg%thickness=bergs%initial_thickness_s(1)
      localberg%width=bergs%initial_width_s(1)
      localberg%length=bergs%initial_length_s(1)
      localberg%fl_k=0.
      localberg%start_lon=localberg%lon
      localberg%start_lat=localberg%lat
      localberg%start_year=iyr
      localberg%start_day=float(iday)+(float(ihr)+float(imin)/60.)/24.
      localberg%start_mass=localberg%mass
      localberg%mass_scaling=bergs%mass_scaling_s(1)
      localberg%mass_of_bits=0.
      localberg%mass_of_fl_bits=0.
      localberg%mass_of_fl_bergy_bits=0.
      localberg%halo_berg=0.
      localberg%static_berg=0.
      localberg%heat_density=0.
      localberg%axn=0. !Alon
      localberg%ayn=0. !Alon
      localberg%uvel_old=0. !Alon
      localberg%vvel_old=0. !Alon
      localberg%bxn=0. !Alon
      localberg%byn=0. !Alon
      localberg%uvel_prev=0. !Alon
      localberg%vvel_prev=0. !Alon
      localberg%lon_prev=0.
      localberg%lat_prev=0.

      if (mts) then
        localberg%axn_fast=0.
        localberg%ayn_fast=0.
        localberg%bxn_fast=0.
        localberg%byn_fast=0.
      endif
      if (dem) then
        localberg%ang_vel=0.
        localberg%ang_accel=0.
        localberg%rot=0.
      elseif (fracture_criterion.ne.'none') then
        localberg%accum_bond_rotation=0.
      endif

      !Berg A
      call loc_set_berg_pos(grd, 0.9, 0.5, 1., 0., localberg)
      localberg%id = generate_id(grd, i, j)
      call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
      !Berg B
      call loc_set_berg_pos(grd, 0.1, 0.5, -1., 0., localberg)
      localberg%id = generate_id(grd, i, j)
      call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
      !Berg C
      call loc_set_berg_pos(grd, 0.5, 0.9, 0., 1., localberg)
      localberg%id = generate_id(grd, i, j)
      call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
      !Berg D
      call loc_set_berg_pos(grd, 0.5, 0.1, 0., -1., localberg)
      localberg%id = generate_id(grd, i, j)
      call add_new_berg_to_list(bergs%list(i,j)%first, localberg)
    endif
  enddo; enddo

  bergs%nbergs_start=count_bergs(bergs)
  call mpp_sum(bergs%nbergs_start)
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'(a,i8,a)') 'KID, generate_bergs: ',bergs%nbergs_start,' were generated'

end subroutine generate_bergs

!> For generating bergs for debugging, assign initial berg positions and velocities
subroutine loc_set_berg_pos(grd, xi, yj, uvel, vvel, berg)
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  real, intent(in) :: xi !< Non-dimensional x-position within cell to give berg
  real, intent(in) :: yj !< Non-dimensional y-position within cell to give berg
  real, intent(in) :: uvel !< Zonal velocity to give berg
  real, intent(in) :: vvel !< Meridional velocity to give berg
  type(iceberg), intent(inout) :: berg !< An iceberg
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

!> Read bond restart file
subroutine read_restart_bonds(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs !< Icebergs container
type(time_type), intent(in) :: Time !< Model time
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
integer, allocatable, dimension(:) :: id_cnt, id_ij,    &
                                      first_berg_jne,   &
                                      first_berg_ine,   &
                                      other_berg_jne,   &
                                      other_berg_ine
real, allocatable, dimension(:) ::    tangd1,           &
                                      tangd2,           &
                                      nstress,          &
                                      sstress,          &
                                      rotation,         &
                                      rel_rotation,     &
                                      n_frac_var,       &
                                      n_strain_rate,    &
                                      spring_pe
integer(kind=8), allocatable, dimension(:) :: first_id,   &
                                              other_id
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
      write(stderrunit,*)  'KID, bond read restart : ','Number of bonds in file',  nbonds_in_file
    endif

  if (nbonds_in_file .gt. 0) then

    allocate(first_id(nbonds_in_file))
    allocate(other_id(nbonds_in_file))
    allocate(id_cnt(nbonds_in_file))
    allocate(id_ij(nbonds_in_file))
    allocate(first_berg_jne(nbonds_in_file))
    allocate(first_berg_ine(nbonds_in_file))
    allocate(other_berg_ine(nbonds_in_file))
    allocate(other_berg_jne(nbonds_in_file))

    if (dem) then
      allocate(tangd1(nbonds_in_file))
      allocate(tangd2(nbonds_in_file))
      allocate(nstress(nbonds_in_file))
      allocate(sstress(nbonds_in_file))
    elseif (fracture_criterion.ne.'none') then
      allocate(rotation(nbonds_in_file))
      allocate(rel_rotation(nbonds_in_file))
      allocate(n_frac_var(nbonds_in_file))
      if (fracture_criterion.eq.'strain_rate') then
        allocate(n_strain_rate(nbonds_in_file))
      elseif (fracture_criterion.eq.'energy') then
        allocate(spring_pe(nbonds_in_file))
      endif
    endif

    call read_unlimited_axis(filename,'first_id_cnt',id_cnt,domain=grd%domain)
    call read_unlimited_axis(filename,'first_id_ij',id_ij,domain=grd%domain)
    do k=1, nbonds_in_file
      first_id(k) = id_from_2_ints( id_cnt(k), id_ij(k) )
    enddo
    call read_unlimited_axis(filename,'other_id_cnt',id_cnt,domain=grd%domain)
    call read_unlimited_axis(filename,'other_id_ij',id_ij,domain=grd%domain)
    do k=1, nbonds_in_file
      other_id(k) = id_from_2_ints( id_cnt(k), id_ij(k) )
    enddo
    deallocate(id_cnt, id_ij)
    call read_unlimited_axis(filename,'first_berg_jne',first_berg_jne,domain=grd%domain)
    call read_unlimited_axis(filename,'first_berg_ine',first_berg_ine,domain=grd%domain)
    call read_unlimited_axis(filename,'other_berg_jne',other_berg_jne,domain=grd%domain)
    call read_unlimited_axis(filename,'other_berg_ine',other_berg_ine,domain=grd%domain)

    if (dem) then
      call read_real_vector(filename,'tangd1',tangd1,grd%domain,value_if_not_in_file=0.)
      call read_real_vector(filename,'tangd2',tangd2,grd%domain,value_if_not_in_file=0.)
      call read_real_vector(filename,'nstress',nstress,grd%domain,value_if_not_in_file=0.)
      call read_real_vector(filename,'sstress',sstress,grd%domain,value_if_not_in_file=0.)
    elseif (fracture_criterion.ne.'none') then
      call read_real_vector(filename,'rotation',rotation,grd%domain,value_if_not_in_file=0.)
      call read_real_vector(filename,'rel_rotation',rel_rotation,grd%domain,value_if_not_in_file=0.)
      call read_real_vector(filename,'n_frac_var',n_frac_var,grd%domain,value_if_not_in_file=0.)
      if (fracture_criterion.eq.'strain_rate') then
        call read_real_vector(filename,'n_strain_rate',n_strain_rate,grd%domain,value_if_not_in_file=0.)
      elseif (fracture_criterion.eq.'energy') then
        call read_real_vector(filename,'spring_pe',spring_pe,grd%domain,value_if_not_in_file=0.)
      endif
    endif

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
         call find_individual_iceberg(bergs,first_id(k), ine, jne,berg_found,search_data_domain)
         berg_found_all_pe=berg_found
         call mpp_sum(berg_found_all_pe)
         if (berg_found_all_pe .gt. 0.5) then
             first_berg_ine(k)=ine
             first_berg_jne(k)=jne
         else
           print * , 'First bond berg not located: ', first_id(k),berg_found, mpp_pe(),ine, jne
           call error_mesg('read_restart_bonds_bergs_new', 'First iceberg in bond not found on any pe', FATAL)
         endif
         !else

         !Finding other iceberg other iceberg
         ine=999 ; jne=999 ; berg_found=0.0 ; search_data_domain =.true.
         call find_individual_iceberg(bergs,other_id(k), ine, jne, berg_found,search_data_domain)
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
                 print * , 'First bond berg not located: ', other_id(k),berg_found, mpp_pe(),ine, jne
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
          if (this%id == first_id(k)) then
            first_berg_found=.true.
            first_berg=>this
            !if (first_berg%halo_berg.gt.0.5) print *, 'bonding halo berg:', first_id(k),  first_berg_ine(k),first_berg_jne(k) ,grd%isc, grd%iec, mpp_pe()
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
            if (this%id == other_id(k)) then
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
            call form_a_bond(first_berg, other_id(k), other_berg_ine(k), other_berg_jne(k),  second_berg)
            number_perfect_bonds=number_perfect_bonds+1

            current_bond=>first_berg%first_bond

            if (dem) then
              current_bond%tangd1=tangd1(k)
              current_bond%tangd2=tangd2(k)
              current_bond%nstress=nstress(k)
              current_bond%sstress=sstress(k)
            elseif (fracture_criterion.ne.'none') then
              current_bond%tangd1=rotation(k)
              current_bond%tangd2=rel_rotation(k)
              current_bond%nstress=n_frac_var(k)
              if (fracture_criterion.eq.'strain_rate') then
                current_bond%n_strain_rate=n_strain_rate(k)
              elseif (fracture_criterion.eq.'energy') then
                current_bond%spring_pe=spring_pe(k)
              endif
            endif

            !Counting number of bonds where the first bond is in the computational domain
            if ( (first_berg_ine(k)>=grd%isc) .and. (first_berg_ine(k)<=grd%iec) .and. &
              (first_berg_jne(k)>=grd%jsc) .and. (first_berg_jne(k)<=grd%jec) ) then
               number_perfect_bonds_with_first_on_pe=number_perfect_bonds_with_first_on_pe+1
            endif

          else
            !print *, 'Forming a bond of the second type', mpp_pe(), first_id(k),  other_id(k)
            !call form_a_bond(first_berg, other_id(k), other_berg_ine(k),other_berg_jne(k))
          endif
        else
          write(stderrunit,*) 'KID, bond read restart : ','Not enough partial bonds formed', k, mpp_pe(), nbonds_in_file
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
      write(stderrunit,*) 'KID, bond read restart : ','Not enough partial bonds formed', all_pe_number_partial_bonds , nbonds_in_file
      call error_mesg('read_restart_bonds_bergs_new', 'Not enough partial bonds formed', FATAL)
    endif

    if (all_pe_number_perfect_bonds .lt. nbonds_in_file) then
      call mpp_sum(all_pe_number_first_bonds_matched)
      call mpp_sum(all_pe_number_second_bonds_matched)
      write(stderrunit,*)  'KID, bond read restart : ','Warning, some bonds are not fully formed',  all_pe_number_first_bonds_matched , nbonds_in_file
      write(stderrunit,*)  'KID, bond read restart : ','Number of first and second bonds matched:', all_pe_number_second_bonds_matched , nbonds_in_file
      call error_mesg('read_restart_bonds_bergs_new', 'Not enough perfect bonds formed', NOTE)
    endif

    if (all_pe_number_perfect_bonds_with_first_on_pe .ne. nbonds_in_file) then
      call mpp_sum(all_pe_number_first_bonds_matched)
      call mpp_sum(all_pe_number_second_bonds_matched)
      write(stderrunit,*)  'KID, bond read restart : ','Warning, # bonds with first bond on computational domain, does not match file',  all_pe_number_first_bonds_matched , nbonds_in_file
      write(stderrunit,*)  'KID, bond read restart : ','Computational bond, first second:', all_pe_number_second_bonds_matched , nbonds_in_file
      call error_mesg('read_restart_bonds_bergs_new', 'Computational perfect bonds do not match those in file', NOTE)
    endif

    deallocate(               &
            first_id,         &
            other_id,         &
            first_berg_ine,   &
            first_berg_jne,   &
            other_berg_ine,   &
            other_berg_jne )

    if (dem) then
      deallocate(                     &
               tangd1,                &
               tangd2,                &
               nstress,               &
               sstress)
    elseif (fracture_criterion.ne.'none') then
      deallocate(                     &
               rotation,              &
               rel_rotation,          &
               n_frac_var)
      if (fracture_criterion.eq.'strain_rate') then
        deallocate(n_strain_rate)
      elseif (fracture_criterion.eq.'energy') then
        deallocate(spring_pe)
      endif
    endif
  endif

  if (mpp_pe() .eq. mpp_root_pe()) then
    write(stderrunit,*)  'KID, bond read restart : ','Number of bonds (including halos)',  all_pe_number_perfect_bonds
    write(stderrunit,*)  'KID, bond read restart : ','Number of true bonds created',  all_pe_number_perfect_bonds_with_first_on_pe
  endif

end subroutine read_restart_bonds

!> Reading calving and gridded restart data
subroutine read_restart_calving(bergs)
use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream
! Arguments
type(icebergs), pointer :: bergs !< Icebergs container
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
     'KID, read_restart_calving: reading ',filename
    call read_data(filename, 'stored_ice', grd%stored_ice, grd%domain)
    if (field_exist(filename, 'stored_heat')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'KID, read_restart_calving: reading stored_heat from restart file.'
      call read_data(filename, 'stored_heat', grd%stored_heat, grd%domain)
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_restart_calving: stored_heat WAS NOT FOUND in the file. Setting to 0.'
      grd%stored_heat(:,:)=0.
    endif
    if (field_exist(filename, 'rmean_calving')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'KID, read_restart_calving: reading rmean_calving from restart file.'
      call read_data(filename, 'rmean_calving', grd%rmean_calving, grd%domain)
      grd%rmean_calving_initialized=.true.
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_restart_calving: rmean_calving WAS NOT FOUND in the file. Setting to 0.'
    endif
    if (field_exist(filename, 'rmean_calving_hflx')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'KID, read_restart_calving: reading rmean_calving_hflx from restart file.'
      call read_data(filename, 'rmean_calving_hflx', grd%rmean_calving_hflx, grd%domain)
      grd%rmean_calving_hflx_initialized=.true.
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_restart_calving: rmean_calving_hflx WAS NOT FOUND in the file. Setting to 0.'
    endif
    if (field_exist(filename, 'iceberg_counter_grd')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'KID, read_restart_calving: reading iceberg_counter_grd from restart file.'
      call read_data(filename, 'iceberg_counter_grd', grd%iceberg_counter_grd, grd%domain)
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_restart_calving: iceberg_counter_grd WAS NOT FOUND in the file. Setting to 0.'
      grd%iceberg_counter_grd(:,:) = 0
    endif
    bergs%restarted=.true.
  else
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_restart_calving: initializing stored ice to random numbers'
    if ( make_calving_reproduce ) then
       allocate(randnum(1,nclasses))
       do j=grd%jsc, grd%jec
          do i=grd%isc, grd%iec
             rns=initializeRandomNumberStream(i+10000*j)
             call getRandomNumbers(rns,randnum(1,:))
             do k=1, nclasses
               if (grd%lat(i,j)<0.) then
                 grd%stored_ice(i,j,k)=randnum(1,k) * grd%msk(i,j) * bergs%initial_mass_s(k) * bergs%mass_scaling_s(k)
               else
                 grd%stored_ice(i,j,k)=randnum(1,k) * grd%msk(i,j) * bergs%initial_mass_n(k) * bergs%mass_scaling_n(k)
               endif
             end do
          end do
       end do
    else
       allocate(randnum(grd%jsc:grd%jec,nclasses))
       do i=grd%isc, grd%iec
          rns = initializeRandomNumberStream(i)
          call getRandomNumbers(rns,randnum)
          do k=1, nclasses
            where (grd%lat(i,grd%jsc:grd%jec)<0.)
             grd%stored_ice(i,grd%jsc:grd%jec,k) = randnum(:,k) * grd%msk(i,grd%jsc:grd%jec) * &
               & bergs%initial_mass_s(k) * bergs%mass_scaling_s(k)
           elsewhere
             grd%stored_ice(i,grd%jsc:grd%jec,k) = randnum(:,k) * grd%msk(i,grd%jsc:grd%jec) * &
               & bergs%initial_mass_n(k) * bergs%mass_scaling_n(k)
           end where
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

!> Read ocean depth from file
subroutine read_ocean_depth(grd)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
! Local variables
character(len=37) :: filename

  ! Read stored ice
  filename=trim(restart_input_dir)//'topog.nc'
  if (file_exist(filename)) then
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') &
     'KID, read_ocean_depth: reading ',filename
    if (field_exist(filename, 'depth')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'KID, read_ocean_depth: reading depth from topog file.'
      call read_data(filename, 'depth', grd%ocean_depth, grd%domain)
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_ocean_depth: depth WAS NOT FOUND in the file. Setting to 0.'
      !grd%ocean_depth(:,:)=0.
    endif
  else
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'KID, read_ocean_depth: Ocean depth file (topog.nc) not present)'
  endif

  !call grd_chksum2(bergs%grd, bergs%grd%ocean_depth, 'read_ocean_depth, ocean_depth')
end subroutine read_ocean_depth

!> Write a trajectory-based diagnostics file
subroutine write_trajectory(trajectory, save_short_traj, save_fl_traj, fl_r)
! Arguments
type(xyt), pointer :: trajectory !< An iceberg trajectory
logical, intent(in) :: save_short_traj !< If true, record less data
logical, intent(in) :: save_fl_traj !< If true, save masses and footloose data
real, intent(in) :: fl_r !< If >0 and save_fl_traj, save footloose params
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, uvelid, vvelid, idcntid, idijid
integer :: uvelpid,vvelpid
integer :: uoid, void, uiid, viid, uaid, vaid, sshxid, sshyid, sstid, sssid
integer :: cnid, hiid, hsid
integer :: mid, smid, did, wid, lid, mbid, mflbid, mflbbid, hdid, nbid, odid, flkid
integer :: axnid,aynid,bxnid,bynid,axnfid,aynfid,bxnfid,bynfid, msid
integer :: eecid,edcid,eeid,edid,aeid,efid
integer :: eectid, edctid, eetid, edtid, aetid
integer :: avid, aaid, rid, abrid
character(len=37) :: filename
character(len=7) :: pe_name
type(xyt), pointer :: this, next
integer :: stderrunit, cnt, ij
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
              call unpack_traj_from_buffer2(traj4io, ibuffer_io, i, save_short_traj, save_fl_traj, fl_r)
           enddo
       endif
     enddo
  else
     ! Pack and send trajectories to the root PE for this I/O tile
     do while (associated(trajectory))
       ntrajs_sent_io = ntrajs_sent_io +1
       call pack_traj_into_buffer2(trajectory, obuffer_io, ntrajs_sent_io, save_short_traj, save_fl_traj, fl_r)
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
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_trajectory: nf_open failed'
    endif
    if (verbose) then
      if (io_is_in_append_mode) then
        write(*,'(2a)') 'KID, write_trajectory: appending to ',filename
      else
        write(*,'(2a)') 'KID, write_trajectory: creating ',filename
      endif
    endif

    if (io_is_in_append_mode) then
      iret = nf_inq_dimid(ncid, 'i', i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_trajectory: nf_inq_dimid i failed'
      lonid = inq_varid(ncid, 'lon')
      latid = inq_varid(ncid, 'lat')
      yearid = inq_varid(ncid, 'year')
      dayid = inq_varid(ncid, 'day')
      idcntid = inq_varid(ncid, 'id_cnt')
      idijid = inq_varid(ncid, 'id_ij')
      if (save_fl_traj) then
        mid = inq_varid(ncid, 'mass')
        smid = inq_varid(ncid, 'start_mass')
        did = inq_varid(ncid, 'thickness')
        mbid = inq_varid(ncid, 'mass_of_bits')
        uvelid = inq_varid(ncid, 'uvel')
        vvelid = inq_varid(ncid, 'vvel')
        if (fl_r>0) then
          msid = inq_varid(ncid, 'mass_scaling')
          mflbid = inq_varid(ncid, 'mass_of_fl_bits')
          mflbbid = inq_varid(ncid, 'mass_of_fl_bergy_bits')
          flkid = inq_varid(ncid, 'fl_k')
        endif
      endif
      if (.not.save_short_traj) then
        !uvelid = inq_varid(ncid, 'uvel')
        !vvelid = inq_varid(ncid, 'vvel')
        uvelpid = inq_varid(ncid, 'uvel_prev')
        vvelpid = inq_varid(ncid, 'vvel_prev')
        uoid = inq_varid(ncid, 'uo')
        void = inq_varid(ncid, 'vo')
        uiid = inq_varid(ncid, 'ui')
        viid = inq_varid(ncid, 'vi')
        uaid = inq_varid(ncid, 'ua')
        vaid = inq_varid(ncid, 'va')
        hdid = inq_varid(ncid, 'heat_density')
        wid = inq_varid(ncid, 'width')
        lid = inq_varid(ncid, 'length')
        sshxid = inq_varid(ncid, 'ssh_x')
        sshyid = inq_varid(ncid, 'ssh_y')
        sstid = inq_varid(ncid, 'sst')
        sssid = inq_varid(ncid, 'sss')
        cnid = inq_varid(ncid, 'cn')
        hiid = inq_varid(ncid, 'hi')
        axnid = inq_varid(ncid, 'axn')
        aynid = inq_varid(ncid, 'ayn')
        bxnid = inq_varid(ncid, 'bxn')
        bynid = inq_varid(ncid, 'byn')
        hsid = inq_varid(ncid, 'halo_berg')
        odid = inq_varid(ncid, 'od')

        if (mts) then
          axnfid = inq_varid(ncid, 'axn_fast')
          aynfid = inq_varid(ncid, 'ayn_fast')
          bxnfid = inq_varid(ncid, 'bxn_fast')
          bynfid = inq_varid(ncid, 'byn_fast')
        endif

        if (iceberg_bonds_on) then
          nbid = inq_varid(ncid, 'n_bonds')
        end if

        if (monitor_energy) then
          eeid = inq_varid(ncid, 'Ee')
          edid = inq_varid(ncid, 'Ed')
          aeid = inq_varid(ncid, 'Eext')
          eecid = inq_varid(ncid, 'Ee_contact')
          edcid = inq_varid(ncid, 'Ed_contact')
          efid =  inq_varid(ncid, 'Efrac')
          ! eectid = inq_varid(ncid, 'Ee_contact_temp')
          ! edctid = inq_varid(ncid, 'Ed_contact_temp')
          ! eetid= inq_varid(ncid, 'Ee_temp')
          ! edtid= inq_varid(ncid, 'Ed_temp')
          ! aetid= inq_varid(ncid, 'Eext_temp')
        endif

        if (dem) then
          avid = inq_varid(ncid, 'ang_vel')
          aaid = inq_varid(ncid, 'ang_accel')
          rid  = inq_varid(ncid, 'rot')
        elseif (fracture_criterion .ne. 'none') then
          abrid = inq_varid(ncid, 'accum_bond_rotation')
        endif

      endif
    else
      ! Dimensions
      iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_trajectory: nf_def_dim i failed'

      ! Variables
      lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
      latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
      yearid = def_var(ncid, 'year', NF_INT, i_dim)
      dayid = def_var(ncid, 'day', NF_DOUBLE, i_dim)
      idcntid = def_var(ncid, 'id_cnt', NF_INT, i_dim)
      idijid = def_var(ncid, 'id_ij', NF_INT, i_dim)
      if (save_fl_traj) then
        mid = def_var(ncid, 'mass', NF_DOUBLE, i_dim)
        smid = def_var(ncid, 'start_mass', NF_DOUBLE, i_dim)
        did = def_var(ncid, 'thickness', NF_DOUBLE, i_dim)
        mbid = def_var(ncid, 'mass_of_bits', NF_DOUBLE, i_dim)
        uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
        vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
        if (fl_r>0) then
          msid = def_var(ncid, 'mass_scaling', NF_DOUBLE, i_dim)
          mflbid = def_var(ncid, 'mass_of_fl_bits', NF_DOUBLE, i_dim)
          mflbbid = def_var(ncid, 'mass_of_fl_bergy_bits', NF_DOUBLE, i_dim)
          flkid = def_var(ncid, 'fl_k', NF_DOUBLE, i_dim)
        endif
      endif
      if (.not. save_short_traj) then
        !uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
        !vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
        uvelpid = def_var(ncid, 'uvel_prev', NF_DOUBLE, i_dim)
        vvelpid = def_var(ncid, 'vvel_prev', NF_DOUBLE, i_dim)
        uoid = def_var(ncid, 'uo', NF_DOUBLE, i_dim)
        void = def_var(ncid, 'vo', NF_DOUBLE, i_dim)
        uiid = def_var(ncid, 'ui', NF_DOUBLE, i_dim)
        viid = def_var(ncid, 'vi', NF_DOUBLE, i_dim)
        uaid = def_var(ncid, 'ua', NF_DOUBLE, i_dim)
        vaid = def_var(ncid, 'va', NF_DOUBLE, i_dim)
        hdid = def_var(ncid, 'heat_density', NF_DOUBLE, i_dim)
        wid = def_var(ncid, 'width', NF_DOUBLE, i_dim)
        lid = def_var(ncid, 'length', NF_DOUBLE, i_dim)
        sshxid = def_var(ncid, 'ssh_x', NF_DOUBLE, i_dim)
        sshyid = def_var(ncid, 'ssh_y', NF_DOUBLE, i_dim)
        sstid = def_var(ncid, 'sst', NF_DOUBLE, i_dim)
        sssid = def_var(ncid, 'sss', NF_DOUBLE, i_dim)
        cnid = def_var(ncid, 'cn', NF_DOUBLE, i_dim)
        hiid = def_var(ncid, 'hi', NF_DOUBLE, i_dim)
        axnid = def_var(ncid, 'axn', NF_DOUBLE, i_dim)
        aynid = def_var(ncid, 'ayn', NF_DOUBLE, i_dim)
        bxnid = def_var(ncid, 'bxn', NF_DOUBLE, i_dim)
        bynid = def_var(ncid, 'byn', NF_DOUBLE, i_dim)
        hsid = def_var(ncid, 'halo_berg', NF_DOUBLE, i_dim)
        odid = def_var(ncid, 'od', NF_DOUBLE, i_dim)

        if (mts) then
          axnfid = def_var(ncid, 'axn_fast', NF_DOUBLE, i_dim)
          aynfid = def_var(ncid, 'ayn_fast', NF_DOUBLE, i_dim)
          bxnfid = def_var(ncid, 'bxn_fast', NF_DOUBLE, i_dim)
          bynfid = def_var(ncid, 'byn_fast', NF_DOUBLE, i_dim)
        endif

        if (iceberg_bonds_on) then
          nbid = def_var(ncid, 'n_bonds', NF_INT, i_dim)
        end if

        if (monitor_energy) then
          eeid = def_var(ncid, 'Ee', NF_DOUBLE, i_dim)
          edid = def_var(ncid, 'Ed', NF_DOUBLE, i_dim)
          aeid = def_var(ncid, 'Eext', NF_DOUBLE, i_dim)
          eecid = def_var(ncid, 'Ee_contact', NF_DOUBLE, i_dim)
          edcid = def_var(ncid, 'Ed_contact', NF_DOUBLE, i_dim)
          efid = def_var(ncid, 'Efrac', NF_DOUBLE, i_dim)
          ! eectid = def_var(ncid, 'Ee_contact_temp', NF_DOUBLE, i_dim)
          ! edctid = def_var(ncid, 'Ed_contact_temp', NF_DOUBLE, i_dim)
          ! eetid= def_var(ncid, 'Ee_temp', NF_DOUBLE, i_dim)
          ! edtid= def_var(ncid, 'Ed_temp', NF_DOUBLE, i_dim)
          ! aetid= def_var(ncid, 'Eext_temp', NF_DOUBLE, i_dim)
        endif

        if (dem) then
          avid = def_var(ncid, 'ang_vel', NF_DOUBLE, i_dim)
          aaid = def_var(ncid, 'ang_accel', NF_DOUBLE, i_dim)
          rid  = def_var(ncid, 'rot', NF_DOUBLE, i_dim)
        elseif (fracture_criterion .ne. 'none') then
          abrid = def_var(ncid, 'accum_bond_rotation', NF_DOUBLE, i_dim)
        endif
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
      call put_att(ncid, idcntid, 'long_name', 'counter component of iceberg id')
      call put_att(ncid, idcntid, 'units', 'dimensionless')
      call put_att(ncid, idijid, 'long_name', 'position component of iceberg id')
      call put_att(ncid, idijid, 'units', 'dimensionless')
      if (save_fl_traj) then
        call put_att(ncid, mid, 'long_name', 'mass')
        call put_att(ncid, mid, 'units', 'kg')
        call put_att(ncid, smid, 'long_name', 'start_mass')
        call put_att(ncid, smid, 'units', 'kg')
        call put_att(ncid, did, 'long_name', 'thickness')
        call put_att(ncid, did, 'units', 'm')
        call put_att(ncid, mbid, 'long_name', 'mass_of_bits')
        call put_att(ncid, mbid, 'units', 'kg')
        call put_att(ncid, uvelid, 'long_name', 'zonal spped')
        call put_att(ncid, uvelid, 'units', 'm/s')
        call put_att(ncid, vvelid, 'long_name', 'meridional spped')
        call put_att(ncid, vvelid, 'units', 'm/s')
        if (fl_r>0) then
          call put_att(ncid, msid, 'long_name', 'mass_scaling')
          call put_att(ncid, msid, 'units', 'dimensionless')
          call put_att(ncid, mflbid, 'long_name', 'mass_of_fl_bits')
          call put_att(ncid, mflbid, 'units', 'kg')
          call put_att(ncid, mflbbid, 'long_name', 'mass_of_fl_bergy_bits')
          call put_att(ncid, mflbbid, 'units', 'kg')
          call put_att(ncid, flkid, 'long_name', 'footloose calving k')
          call put_att(ncid, flkid, 'units', 'none')
        endif
      endif
      if (.not. save_short_traj) then
        ! call put_att(ncid, uvelid, 'long_name', 'zonal spped')
        ! call put_att(ncid, uvelid, 'units', 'm/s')
        ! call put_att(ncid, vvelid, 'long_name', 'meridional spped')
        ! call put_att(ncid, vvelid, 'units', 'm/s')
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
        call put_att(ncid, hdid, 'long_name', 'heat_density')
        call put_att(ncid, hdid, 'units', 'J/kg')
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
        call put_att(ncid, axnid, 'long_name', 'explicit zonal acceleration')
        call put_att(ncid, axnid, 'units', 'm')
        call put_att(ncid, aynid, 'long_name', 'explicit meridional acceleration')
        call put_att(ncid, aynid, 'units', 'm')
        call put_att(ncid, bxnid, 'long_name', 'implicit zonal acceleration')
        call put_att(ncid, bxnid, 'units', 'm')
        call put_att(ncid, bynid, 'long_name', 'implicit meridional acceleration')
        call put_att(ncid, bynid, 'units', 'm')
        call put_att(ncid, hsid, 'long_name', 'halo status')
        call put_att(ncid, hsid, 'units', 'non-dim')
        call put_att(ncid, odid, 'long_name', 'ocean_depth')
        call put_att(ncid, odid, 'units', 'm')

        if (mts) then
          call put_att(ncid, uvelpid, 'long_name', 'zonal speed mts')
          call put_att(ncid, uvelpid, 'units', 'm/s')
          call put_att(ncid, vvelpid, 'long_name', 'meridional speed mts')
          call put_att(ncid, vvelpid, 'units', 'm/s')
          call put_att(ncid, axnfid, 'long_name', 'explicit fast step zonal acceleration')
          call put_att(ncid, axnfid, 'units', 'm')
          call put_att(ncid, aynfid, 'long_name', 'explicit fast step meridional acceleration')
          call put_att(ncid, aynfid, 'units', 'm')
          call put_att(ncid, bxnfid, 'long_name', 'implicit fast step zonal acceleration')
          call put_att(ncid, bxnfid, 'units', 'm')
          call put_att(ncid, bynfid, 'long_name', 'implicit fast step meridional acceleration')
          call put_att(ncid, bynfid, 'units', 'm')
        endif

        if (iceberg_bonds_on) then
          call put_att(ncid, nbid, 'long_name', 'number of bonds')
          call put_att(ncid, nbid, 'units', 'dimensionless')
        end if

        if (monitor_energy) then
          call put_att(ncid, eeid, 'long_name', 'bonded elastic energy')
          call put_att(ncid, eeid, 'units', 'J')
          call put_att(ncid, edid, 'long_name', 'bonded damping energy')
          call put_att(ncid, edid, 'units', 'J')
          call put_att(ncid, aeid, 'long_name', 'external energy')
          call put_att(ncid, aeid, 'units', 'J')
          call put_att(ncid, eecid, 'long_name', 'contact elastic energy')
          call put_att(ncid, eecid, 'units', 'J')
          call put_att(ncid, edcid, 'long_name', 'contact damping energy')
          call put_att(ncid, edcid, 'units', 'J')
          call put_att(ncid, efid, 'long_name', 'dissipated fracture energy')
          call put_att(ncid, efid, 'units', 'J')
          ! call put_att(ncid, eectid, 'long_name', 'contact elastic energy temp')
          ! call put_att(ncid, eectid, 'units', 'J')
          ! call put_att(ncid, edctid, 'long_name', 'contact damping energy temp')
          ! call put_att(ncid, edctid, 'units', 'J')
          ! call put_att(ncid, eetid, 'long_name', 'bonded elastic energy temp')
          ! call put_att(ncid, eetid, 'units', 'J')
          ! call put_att(ncid, edtid, 'long_name', 'bonded damping energy temp')
          ! call put_att(ncid, edtid, 'units', 'J')
          ! call put_att(ncid, aetid, 'long_name', 'external energy temp')
          ! call put_att(ncid, aetid, 'units', 'J')
        endif

        if (dem) then
          call put_att(ncid, avid, 'long_name', 'angular velocity')
          call put_att(ncid, avid, 'units', 'rad/s')
          call put_att(ncid, aaid, 'long_name', 'angular acceleration')
          call put_att(ncid, aaid, 'units', 'rad/s^2')
          call put_att(ncid, rid, 'long_name', 'accumulated rotation')
          call put_att(ncid, rid, 'units', 'rad')
        elseif (fracture_criterion .ne. 'none') then
          call put_att(ncid, abrid, 'long_name', 'accumulated bond rotation')
          call put_att(ncid, abrid, 'units', 'radians')
        endif
      endif
    endif

    ! End define mode
    iret = nf_enddef(ncid)

    ! Write variables
    this=>traj4io
    if (io_is_in_append_mode) then
      iret = nf_inq_dimlen(ncid, i_dim, i)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_trajectory: nf_inq_dimlen i failed'
    else
      i = 0
    endif
    do while (associated(this))
      i=i+1
      call put_double(ncid, lonid, i, this%lon)
      call put_double(ncid, latid, i, this%lat)
      call put_int(ncid, yearid, i, this%year)
      call put_double(ncid, dayid, i, this%day)
      call split_id(this%id, cnt, ij)
      call put_int(ncid, idcntid, i, cnt)
      call put_int(ncid, idijid, i, ij)
      if (save_fl_traj) then
        call put_double(ncid, mid, i, this%mass)
        call put_double(ncid, smid, i, this%start_mass)
        call put_double(ncid, did, i, this%thickness)
        call put_double(ncid, mbid, i, this%mass_of_bits)
        call put_double(ncid, uvelid, i, this%uvel)
        call put_double(ncid, vvelid, i, this%vvel)
        if (fl_r>0) then
          call put_double(ncid, msid, i, this%mass_scaling)
          call put_double(ncid, mflbid, i, this%mass_of_fl_bits)
          call put_double(ncid, mflbbid, i, this%mass_of_fl_bergy_bits)
          call put_double(ncid, flkid, i, this%fl_k)
        endif
      endif
      if (.not. save_short_traj) then
        ! call put_double(ncid, uvelid, i, this%uvel)
        ! call put_double(ncid, vvelid, i, this%vvel)
        call put_double(ncid, uvelpid, i, this%uvel_prev)
        call put_double(ncid, vvelpid, i, this%vvel_prev)
        call put_double(ncid, uoid, i, this%uo)
        call put_double(ncid, void, i, this%vo)
        call put_double(ncid, uiid, i, this%ui)
        call put_double(ncid, viid, i, this%vi)
        call put_double(ncid, uaid, i, this%ua)
        call put_double(ncid, vaid, i, this%va)
        call put_double(ncid, hdid, i, this%heat_density)
        call put_double(ncid, wid, i, this%width)
        call put_double(ncid, lid, i, this%length)
        call put_double(ncid, sshxid, i, this%ssh_x)
        call put_double(ncid, sshyid, i, this%ssh_y)
        call put_double(ncid, sstid, i, this%sst)
        call put_double(ncid, sssid, i, this%sss)
        call put_double(ncid, cnid, i, this%cn)
        call put_double(ncid, hiid, i, this%hi)
        call put_double(ncid, axnid, i, this%axn)
        call put_double(ncid, aynid, i, this%ayn)
        call put_double(ncid, bxnid, i, this%bxn)
        call put_double(ncid, bynid, i, this%byn)
        call put_double(ncid, hsid, i, this%halo_berg)
        call put_double(ncid, odid, i, this%od)

        if (mts) then
          call put_double(ncid, axnfid, i, this%axn_fast)
          call put_double(ncid, aynfid, i, this%ayn_fast)
          call put_double(ncid, bxnfid, i, this%bxn_fast)
          call put_double(ncid, bynfid, i, this%byn_fast)
        endif

        if (iceberg_bonds_on) then
          call put_int(ncid, nbid, i, this%n_bonds)
        end if

        if (monitor_energy) then
          call put_double(ncid, eeid,   i, this%Ee)
          call put_double(ncid, edid,   i, this%Ed)
          call put_double(ncid, aeid,   i, this%Eext)
          call put_double(ncid, eecid,  i, this%Ee_contact)
          call put_double(ncid, edcid,  i, this%Ed_contact)
          call put_double(ncid, efid,   i, this%Efrac)
          ! call put_double(ncid, eectid, i, this%Ee_contact_temp)
          ! call put_double(ncid, edctid, i, this%Ed_contact_temp)
          ! call put_double(ncid, eetid,  i, this%Ee_temp)
          ! call put_double(ncid, edtid,  i, this%Ed_temp)
          ! call put_double(ncid, aetid,  i, this%Eext_temp)
        endif

        if (dem) then
          call put_double(ncid, avid, i, this%ang_vel)
          call put_double(ncid, aaid, i, this%ang_accel)
          call put_double(ncid, rid,  i, this%rot)
        elseif (fracture_criterion .ne. 'none') then
          call put_double(ncid, abrid, i, this%accum_bond_rotation)
        endif

      endif
      next=>this%next
      deallocate(this)
      this=>next
    enddo

    ! Finish up
    iret = nf_close(ncid)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_trajectory: nf_close failed',mpp_pe(),filename

  endif !(is_io_tile_root_pe .AND. associated(traj4io))
  call mpp_clock_end(clock_trw)

end subroutine write_trajectory

!> Write a trajectory-based diagnostics file for bonds
subroutine write_bond_trajectory(trajectory)
! Arguments
type(bond_xyt), pointer :: trajectory !< An iceberg bond trajectory
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, lenid,n1id, n2id, peid
integer :: rotid,rrotid,nsid,nsrid, damid
integer :: idcnt1_id, idcnt2_id, idij1_id, idij2_id, eeid, edid
integer :: axid,ayid,bxid,byid
integer :: td1id,td2id,dnsid,dssid
character(len=34) :: filename
character(len=7) :: pe_name
type(bond_xyt), pointer :: this, next
integer :: stderrunit, cnt, ij
!I/O vars
type(bond_xyt), pointer :: traj4io=>null()
integer :: ntrajs_sent_io,ntrajs_rcvd_io
integer :: from_pe,np
type(buffer), pointer :: obuffer_io=>null(), ibuffer_io=>null()
logical :: io_is_in_append_mode

  stderrunit=stderr()  ! Get the stderr unit number
  traj4io=>null()
  ibuffer_io=>null(); obuffer_io=>null()

  !Assemble the list of trajectories from all pes in this I/O tile
  call mpp_clock_begin(clock_btrp)

  !First add the trajs on the io_tile_root_pe (if any) to the I/O list
  if(is_io_tile_root_pe .OR. force_all_pes_traj ) then
    if(associated(trajectory)) then
      this=>trajectory
      do while (associated(this))
        call append_bond_posn(traj4io, this)
        this=>this%next
      enddo
      trajectory => null()
    endif
  endif

  if(.NOT. force_all_pes_traj ) then
    !Now gather and append the bergs from all pes in the io_tile to the list on corresponding io_tile_root_pe
    ntrajs_sent_io =0; ntrajs_rcvd_io =0

    if(is_io_tile_root_pe) then
      !Receive trajs from all pes in this I/O tile !FRAGILE!SCARY!
      do np=2,size(io_tile_pelist) ! Note: np starts from 2 to exclude self
        from_pe=io_tile_pelist(np)
        call mpp_recv(ntrajs_rcvd_io, glen=1, from_pe=from_pe, tag=COMM_TAG_11)
        if (ntrajs_rcvd_io .gt. 0) then
          call increase_ibuffer(ibuffer_io, ntrajs_rcvd_io,buffer_width_bond_traj)
          call mpp_recv(ibuffer_io%data, ntrajs_rcvd_io*buffer_width_bond_traj,from_pe=from_pe, tag=COMM_TAG_12)
          do i=1, ntrajs_rcvd_io
            call unpack_bond_traj_from_buffer2(traj4io, ibuffer_io, i)
          enddo
        endif
      enddo
    else
      ! Pack and send trajectories to the root PE for this I/O tile
      do while (associated(trajectory))
        ntrajs_sent_io = ntrajs_sent_io +1
        call pack_bond_traj_into_buffer2(trajectory, obuffer_io, ntrajs_sent_io)
        this => trajectory ! Need to keep pointer in order to free up the links memory
        trajectory => trajectory%next ! This will eventually result in trajectory => null()
        deallocate(this) ! Delete the link from memory
      enddo

      call mpp_send(ntrajs_sent_io, plen=1, to_pe=io_tile_root_pe, tag=COMM_TAG_11)
      if (ntrajs_sent_io .gt. 0) then
        call mpp_send(obuffer_io%data, ntrajs_sent_io*buffer_width_bond_traj, to_pe=io_tile_root_pe, tag=COMM_TAG_12)
      endif
    endif
  endif !.NOT. force_all_pes_traj

  call mpp_clock_end(clock_btrp)

  !Now start writing in the io_tile_root_pe if there are any bergs in the I/O list
  call mpp_clock_begin(clock_btrw)

  if((force_all_pes_traj .OR. is_io_tile_root_pe) .AND. associated(traj4io)) then

    call get_instance_filename("bond_trajectories.nc", filename)
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
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_bond_trajectory: nf_open failed'
    endif
    if (verbose) then
      if (io_is_in_append_mode) then
        write(*,'(2a)') 'KID, write_bond_trajectory: appending to ',filename
      else
        write(*,'(2a)') 'KID, write_bond_trajectory: creating ',filename
      endif
    endif

    if (io_is_in_append_mode) then
      iret = nf_inq_dimid(ncid, 'i', i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_bond_trajectory: nf_inq_dimid i failed'
      lonid = inq_varid(ncid, 'lon');         latid = inq_varid(ncid, 'lat')
      yearid = inq_varid(ncid, 'year');       dayid = inq_varid(ncid, 'day')
      lenid = inq_varid(ncid, 'length');
      n1id = inq_varid(ncid, 'n1');           n2id = inq_varid(ncid, 'n2')

      idcnt1_id = inq_varid(ncid, 'id_cnt1'); idij1_id = inq_varid(ncid, 'id_ij1')
      idcnt2_id = inq_varid(ncid, 'id_cnt2'); idij2_id = inq_varid(ncid, 'id_ij2')

      if (use_damage) then
        damid = inq_varid(ncid, 'damage')
      endif

      if (monitor_energy) then
        eeid = inq_varid(ncid, 'Ee');       edid = inq_varid(ncid, 'Ed')
        axid = inq_varid(ncid, 'axn_fast'); ayid = inq_varid(ncid, 'ayn_fast')
        bxid = inq_varid(ncid, 'bxn_fast'); byid = inq_varid(ncid, 'byn_fast')
      endif

      if (dem) then
        td1id = inq_varid(ncid, 'tangd1');  td2id = inq_varid(ncid, 'tangd2')
        dnsid = inq_varid(ncid, 'nstress'); dssid = inq_varid(ncid, 'sstress')
      elseif (fracture_criterion.ne.'none') then
        rotid = inq_varid(ncid, 'rotation');    rrotid = inq_varid(ncid, 'rel_rotation')
        nsid = inq_varid(ncid, 'n_frac_var')
        if (fracture_criterion.eq.'strain_rate') then
          nsrid = inq_varid(ncid, 'n_strain_rate')
        elseif (fracture_criterion.eq.'energy') then
          peid = inq_varid(ncid, 'spring_pe')
        endif
      endif

    else
      ! Dimensions
      iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_bond_trajectory: nf_def_dim i failed'

      ! Variables
      lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim);      latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
      yearid = def_var(ncid, 'year', NF_INT, i_dim);       dayid = def_var(ncid, 'day', NF_DOUBLE, i_dim)
      lenid = def_var(ncid, 'length', NF_DOUBLE, i_dim);
      n1id = def_var(ncid, 'n1', NF_DOUBLE, i_dim);        n2id = def_var(ncid, 'n2', NF_DOUBLE, i_dim)

      idcnt1_id = def_var(ncid, 'id_cnt1', NF_INT, i_dim); idij1_id = def_var(ncid, 'id_ij1', NF_INT, i_dim)
      idcnt2_id = def_var(ncid, 'id_cnt2', NF_INT, i_dim); idij2_id = def_var(ncid, 'id_ij2', NF_INT, i_dim)

      if (use_damage) then
        damid = def_var(ncid, 'damage', NF_DOUBLE, i_dim)
      endif

      if (monitor_energy) then
        eeid = def_var(ncid, 'Ee', NF_DOUBLE, i_dim);       edid = def_var(ncid, 'Ed', NF_DOUBLE, i_dim)
        axid = def_var(ncid, 'axn_fast', NF_DOUBLE, i_dim); ayid = def_var(ncid, 'ayn_fast', NF_DOUBLE, i_dim)
        bxid = def_var(ncid, 'bxn_fast', NF_DOUBLE, i_dim); byid = def_var(ncid, 'byn_fast', NF_DOUBLE, i_dim)
      endif

      if (dem) then
        td1id = def_var(ncid, 'tangd1',  NF_DOUBLE, i_dim); td2id = def_var(ncid, 'tangd2',  NF_DOUBLE, i_dim)
        dnsid = def_var(ncid, 'nstress', NF_DOUBLE, i_dim); dssid = def_var(ncid, 'sstress', NF_DOUBLE, i_dim)
      elseif (fracture_criterion.ne.'none') then
        rotid = def_var(ncid, 'rotation', NF_DOUBLE, i_dim); rrotid = def_var(ncid, 'rel_rotation', NF_DOUBLE, i_dim)
        nsid = def_var(ncid, 'n_frac_var', NF_DOUBLE, i_dim)
        if (fracture_criterion.eq.'strain_rate') then
          nsrid = def_var(ncid, 'n_strain_rate', NF_DOUBLE, i_dim)
        elseif (fracture_criterion.eq.'energy') then
          peid = def_var(ncid, 'spring_pe', NF_DOUBLE, i_dim)
        endif
      endif

      ! Attributes
      iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_major_version', NF_INT, 1, 0)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_minor_version', NF_INT, 1, 1)
      call put_att(ncid, lonid, 'long_name', 'longitude');    call put_att(ncid, lonid, 'units', 'degrees_E')
      call put_att(ncid, latid, 'long_name', 'latitude');     call put_att(ncid, latid, 'units', 'degrees_N')
      call put_att(ncid, yearid, 'long_name', 'year');        call put_att(ncid, yearid, 'units', 'years')
      call put_att(ncid, dayid, 'long_name', 'year day');     call put_att(ncid, dayid, 'units', 'days')
      call put_att(ncid, lenid, 'long_name', 'length');       call put_att(ncid, lenid, 'units', 'm')
      call put_att(ncid, n1id, 'long_name', 'unit vector 1'); call put_att(ncid, n1id, 'units', 'non-dim')
      call put_att(ncid, n2id, 'long_name', 'unit vector 2'); call put_att(ncid, n2id, 'units', 'non-dim')

      call put_att(ncid, idcnt1_id, 'long_name', 'counter component of first connected iceberg id')
      call put_att(ncid, idcnt1_id, 'units', 'dimensionless')
      call put_att(ncid, idij1_id, 'long_name', 'position component of first connected iceberg id')
      call put_att(ncid, idij1_id, 'units', 'dimensionless')
      call put_att(ncid, idcnt2_id, 'long_name', 'counter component of second connected iceberg id')
      call put_att(ncid, idcnt2_id, 'units', 'dimensionless')
      call put_att(ncid, idij2_id, 'long_name', 'position component of second connected iceberg id')
      call put_att(ncid, idij2_id, 'units', 'dimensionless')

      if (use_damage) then
        call put_att(ncid, damid, 'long_name', 'damage')
        call put_att(ncid, damid, 'units', 'none')
      endif

      if (monitor_energy) then
        call put_att(ncid, eeid, 'long_name', 'elastic energy'); call put_att(ncid, eeid, 'units', 'J')
        call put_att(ncid, edid, 'long_name', 'damping energy'); call put_att(ncid, edid, 'units', 'J')
        call put_att(ncid, axid, 'long_name', 'explicit short-step zonal acceleration')
        call put_att(ncid, axid, 'units', 'm/s^2')
        call put_att(ncid, ayid, 'long_name', 'explicit short-step meridional acceleration')
        call put_att(ncid, ayid, 'units', 'm/s^2')
        call put_att(ncid, bxid, 'long_name', 'implicit short-step zonal acceleration')
        call put_att(ncid, bxid, 'units', 'm/s^2')
        call put_att(ncid, byid, 'long_name', 'implicit short-step meridional acceleration')
        call put_att(ncid, byid, 'units', 'm/s^2')
      endif

      if (dem) then
        call put_att(ncid, td1id, 'long_name', 'zonal tangential displacement')
        call put_att(ncid, td1id, 'units', 'm')
        call put_att(ncid, td2id, 'long_name', 'meridional tangential displacement')
        call put_att(ncid, td2id, 'units', 'm')
        call put_att(ncid, dnsid, 'long_name', 'normal stress'); call put_att(ncid, dnsid, 'units', 'Pa')
        call put_att(ncid, dssid, 'long_name', 'shear stress');  call put_att(ncid, dssid, 'units', 'Pa')
      elseif (fracture_criterion.ne.'none') then
        call put_att(ncid, rotid, 'long_name', 'rotation');                call put_att(ncid, rotid, 'units', 'rad')
        call put_att(ncid, rrotid, 'long_name', 'relative rotation');      call put_att(ncid, rrotid, 'units', 'rad')
        call put_att(ncid, nsid, 'long_name', 'normal fracture variable'); call put_att(ncid, nsid, 'units', 'varies')
        if (fracture_criterion.eq.'strain_rate') then
          call put_att(ncid, nsrid, 'long_name', 'normal strain-rate'); call put_att(ncid, nsrid, 'units', '1/seconds')
        elseif (fracture_criterion.eq.'energy') then
          call put_att(ncid, peid, 'long_name', 'spring potential energy'); call put_att(ncid, peid, 'units', 'J')
        endif
      endif
    endif

    ! End define mode
    iret = nf_enddef(ncid)

    ! Write variables
    this=>traj4io
    if (io_is_in_append_mode) then
      iret = nf_inq_dimlen(ncid, i_dim, i)
      if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_bond_trajectory: nf_inq_dimlen i failed'
    else
      i = 0
    endif
    do while (associated(this))
      i=i+1
      call put_double(ncid, lonid, i, this%lon);    call put_double(ncid, latid, i, this%lat)
      call put_int(ncid, yearid, i, this%year);     call put_double(ncid, dayid, i, this%day)
      call put_double(ncid, lenid, i, this%length);
      call put_double(ncid, n1id, i, this%n1);      call put_double(ncid, n2id, i, this%n2)

      call split_id(this%id1, cnt, ij)
      call put_int(ncid, idcnt1_id, i, cnt);        call put_int(ncid, idij1_id, i, ij)
      call split_id(this%id2, cnt, ij)
      call put_int(ncid, idcnt2_id, i, cnt);        call put_int(ncid, idij2_id, i, ij)

      if (use_damage) then
        call put_double(ncid,damid,i,this%damage)
      endif

      if (monitor_energy) then
        call put_double(ncid,eeid,i,this%Ee);       call put_double(ncid,edid,i,this%Ed)
        call put_double(ncid,axid,i,this%axn_fast); call put_double(ncid,ayid,i,this%ayn_fast)
        call put_double(ncid,bxid,i,this%bxn_fast); call put_double(ncid,byid,i,this%byn_fast)
      endif

      if (dem) then
        call put_double(ncid,td1id,i,this%tangd1);  call put_double(ncid,td2id,i,this%tangd2)
        call put_double(ncid,dnsid,i,this%nstress); call put_double(ncid,dssid,i,this%sstress)
      elseif (fracture_criterion.ne.'none') then
        call put_double(ncid,rotid,i,this%rotation);  call put_double(ncid,rrotid,i,this%rel_rotation)
        call put_double(ncid,nsid,i,this%n_frac_var)
        if (fracture_criterion.eq.'strain_rate') then
          call put_double(ncid,nsrid,i,this%n_strain_rate)
        elseif (fracture_criterion.eq.'energy') then
          call put_double(ncid,peid,i,this%spring_pe)
        endif
      endif

      next=>this%next
      deallocate(this)
      this=>next
    enddo

    ! Finish up
    iret = nf_close(ncid)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'KID, write_bond_trajectory: nf_close failed',mpp_pe(),filename

  endif !(is_io_tile_root_pe .AND. associated(traj4io))
  call mpp_clock_end(clock_btrw)

end subroutine write_bond_trajectory

!> Returns netcdf id of variable
integer function inq_var(ncid, var, unsafe)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
character(len=*), intent(in) :: var !< Name of variable
logical, optional, intent(in) :: unsafe !< If present and true, do not fail if variable is not in file
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
      write(stderrunit,*) 'KID, inq_var: nf_inq_varid ',var,' failed'
      call error_mesg('KID, inq_var', 'netcdf function returned a failure!', FATAL)
    else
      inq_var=-1
    endif
  endif

end function inq_var

!> Define a netcdf variable
integer function def_var(ncid, var, ntype, idim)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
character(len=*), intent(in) :: var !< Name of variable
integer, intent(in) :: ntype !< Netcdf type of variable
integer, intent(in) :: idim !< Length of vector
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_def_var(ncid, var, ntype, 1, idim, def_var)
  if (iret .ne. NF_NOERR) then
    call error_mesg('KID, def_var', nf_strerror(iret), FATAL)
  endif

end function def_var

!> Returns id of variable
integer function inq_varid(ncid, var)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
character(len=*), intent(in) :: var !< Name of variable
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_inq_varid(ncid, var, inq_varid)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'KID, inq_varid: nf_inq_varid failed for ',trim(var)
    call error_mesg('KID, inq_varid', 'netcdf function returned a failure!', FATAL)
  endif

end function inq_varid

!> Add a string attribute to a netcdf variable
subroutine put_att(ncid, id, att, attval)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
character(len=*), intent(in) :: att !< Name of attribute
character(len=*), intent(in) :: attval !< Value of attribute
! Local variables
integer :: vallen, iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  vallen=len_trim(attval)
  iret = nf_put_att_text(ncid, id, att, vallen, attval)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'KID, put_att: nf_put_att_text failed adding', &
      trim(att),' = ',trim(attval)
    call error_mesg('KID, put_att', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_att

!> Read a real from a netcdf file
real function get_double(ncid, id, i)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
integer, intent(in) :: i !< Index to read
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_double(ncid, id, i, get_double)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'KID, get_double: nf_get_var1_double failed reading'
    call error_mesg('KID, get_double', 'netcdf function returned a failure!', FATAL)
  endif

end function get_double

!> Read an integer from a netcdf file
integer function get_int(ncid, id, i)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
integer, intent(in) :: i !< Index to read
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_int(ncid, id, i, get_int)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'KID, get_int: nf_get_var1_int failed reading'
    call error_mesg('KID, get_int', 'netcdf function returned a failure!', FATAL)
  endif

end function get_int

!> Write a real to a netcdf file
subroutine put_double(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
integer, intent(in) :: i !< Index of position to write
real, intent(in) :: val !< Value to write
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_double(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'KID, put_double: nf_put_vara_double failed writing'
    call error_mesg('KID, put_double', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_double

!> Write an integer to a netcdf file
subroutine put_int(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid !< Handle to netcdf file
integer, intent(in) :: id !< Netcdf id of variable
integer, intent(in) :: i !< Index of position to write
integer, intent(in) :: val !< Value to write
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_int(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'KID, put_int: nf_put_vara_int failed writing'
    call error_mesg('KID, put_int', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_int

!> True is a restart file can be found
logical function find_restart_file(filename, actual_file, multiPErestart, tile_id)
  character(len=*), intent(in) :: filename !< Base-name of restart file
  character(len=*), intent(out) :: actual_file !< Actual name of file, if found
  logical, intent(out) :: multiPErestart !< True if found, false otherwise
  integer, intent(in) :: tile_id !< Parallel tile number of file

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
