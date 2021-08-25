!> Top-level/entry functions that step forward the governing equations
module ice_bergs

! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.

use constants_mod, only: pi, omega, HLF
use fms_mod, only: open_namelist_file, check_nml_error, close_file
use fms_mod, only: field_exist, get_global_att_value
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: write_version_number, read_data, write_data, file_exist
use mosaic_mod, only: get_mosaic_ntiles, get_mosaic_ncontacts
use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_chksum, mpp_sync
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP
use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream
use random_numbers_mod, only: constructSeed
use mpp_mod, only: mpp_gather
use fms_mod, only: clock_flag_default
use fms_io_mod, only: get_instance_filename
use mpp_domains_mod, only: domain2D, mpp_update_domains, mpp_define_domains
use mpp_parameter_mod, only: CGRID_NE, BGRID_NE, CORNER, AGRID
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use diag_manager_mod, only: diag_axis_init

use ice_bergs_framework, only: ice_bergs_framework_init
use ice_bergs_framework, only: icebergs_gridded, xyt, iceberg, icebergs, buffer, bond
use ice_bergs_framework, only: verbose, really_debug,debug,old_bug_rotated_weights,budget,use_roundoff_fix
use ice_bergs_framework, only: find_cell,find_cell_by_search,find_cell_wide,count_bergs,is_point_in_cell,pos_within_cell
use ice_bergs_framework, only: count_bonds, form_a_bond,connect_all_bonds,show_all_bonds, bond_address_update
use ice_bergs_framework, only: nclasses,old_bug_bilin
use ice_bergs_framework, only: sum_mass,sum_heat,bilin,yearday,count_bergs,bergs_chksum,count_bergs_in_list
use ice_bergs_framework, only: checksum_gridded,add_new_berg_to_list,list_chksum
use ice_bergs_framework, only: send_bergs_to_other_pes,move_trajectory,move_all_trajectories
use ice_bergs_framework, only: move_berg_between_cells, update_halo_icebergs
use ice_bergs_framework, only: record_posn,check_position,print_berg,print_bergs,print_fld
use ice_bergs_framework, only: add_new_berg_to_list,delete_iceberg_from_list,destroy_iceberg
use ice_bergs_framework, only: grd_chksum2,grd_chksum3
use ice_bergs_framework, only: fix_restart_dates, offset_berg_dates
use ice_bergs_framework, only: orig_read  ! Remove when backward compatibility no longer needed
use ice_bergs_framework, only: monitor_a_berg
use ice_bergs_framework, only: is_point_within_xi_yj_bounds
use ice_bergs_framework, only: test_check_for_duplicate_ids_in_list
use ice_bergs_framework, only: generate_id
use ice_bergs_framework, only: update_latlon,set_conglom_ids,transfer_mts_bergs,quad_interp_from_agrid
use ice_bergs_framework, only: save_bond_traj
use ice_bergs_framework, only: fracture_testing_initialization, orig_bond_length
use ice_bergs_framework, only: update_and_break_bonds,break_bonds_dem,assign_n_bonds,reset_bond_rotation
use ice_bergs_framework, only: update_bond_angles
use ice_bergs_framework, only: monitor_energy, energy_tests_init, mts, new_mts
use ice_bergs_framework, only: dem_tests_init
use ice_bergs_framework, only: dem, init_dem_params
use ice_bergs_framework, only: set_constant_interaction_length_and_width
use ice_bergs_framework, only: use_damage, damage_test_1_init

use ice_bergs_io,        only: ice_bergs_io_init, write_restart, write_trajectory, write_bond_trajectory
use ice_bergs_io,        only: read_restart_bergs, read_restart_calving
use ice_bergs_io,        only: read_restart_bonds
use ice_bergs_io,        only: read_ocean_depth

implicit none ; private

public icebergs_init, icebergs_end, icebergs_run, icebergs_stock_pe, icebergs
public icebergs_incr_mass, icebergs_save_restart

real, parameter :: pi_180=pi/180.  !< Converts degrees to radians
real, parameter :: r180_pi=180./pi !< Converts radians to degrees
real, parameter :: Rearth=6360000. !< Radius of earth (m)
real, parameter :: rho_ice=916.7 !< Density of fresh ice @ 0oC (kg/m^3)
real, parameter :: rho_water=999.8 !< Density of fresh water @ 0oC (kg/m^3)
real, parameter :: rho_air=1.1 !< Density of air @ 0oC (kg/m^3)
real, parameter :: rho_seawater=1025. !< Approx. density of surface sea water @ 0oC (kg/m^3)
real, parameter :: gravity=9.8 !< Gravitational acceleratio (m/s^2)
real, parameter :: Cd_av=1.3 !< (Vertical) Drag coefficient between bergs and atmos
real, parameter :: Cd_ah=0.0055 !< (Horizontal) Drag coefficient between bergs and atmos
real, parameter :: Cd_wv=0.9 !< (Vertical) Drag coefficient between bergs and ocean
real, parameter :: Cd_wh=0.0012 !< (Horizontal) Drag coefficient between bergs and ocean
real, parameter :: Cd_iv=0.9 !< (Vertical) Drag coefficient between bergs and sea-ice
!TOM> no horizontal drag for sea ice! real, parameter :: Cd_ih=0.0012 !< (Horizontal) Drag coefficient between bergs and sea-ice

#ifdef _FILE_VERSION
character(len=128) :: version = _FILE_VERSION !< Version of file
#else
character(len=128) :: version = 'unknown' !< Version of file
#endif

contains

!> Initializes icebergs container "bergs"
subroutine icebergs_init(bergs, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth, maskmap, fractional_area)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer, intent(in) :: gni !< Number of global points in i-direction
  integer, intent(in) :: gnj !< Number of global points in j-direction
  integer, intent(in) :: layout(2) !< Parallel decomposition of computational processors in i/j direction
  integer, intent(in) :: io_layout(2) !< Parallel decomposition of i/o processors in i/j direction
  integer, intent(in) :: axes(2) !< Diagnostic axes
  integer, intent(in) :: dom_x_flags !< Decomposition flags for i-direction
  integer, intent(in) :: dom_y_flags !< Decomposition flags for j-direction
  real, intent(in) :: dt !< Time step (s)
  type (time_type), intent(in) :: Time !< Model time
  real, dimension(:,:), intent(in) :: ice_lon !< Longitude of cell corners using NE convention (degree E)
  real, dimension(:,:), intent(in) :: ice_lat !< Latitude of cell corners using NE conventino (degree N)
  real, dimension(:,:), intent(in) :: ice_wet !< Wet/dry mask (1 is wet, 0 is dry) of cell centers
  real, dimension(:,:), intent(in) :: ice_dx !< Zonal length of cell on northern side (m)
  real, dimension(:,:), intent(in) :: ice_dy !< Meridional length of cell on eastern side (m)
  real, dimension(:,:), intent(in) :: ice_area !< Area of cells (m^2, or non-dim is fractional_area=True)
  real, dimension(:,:), intent(in) :: cos_rot !< Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), intent(in) :: sin_rot !< Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), intent(in),optional :: ocean_depth !< Depth of ocean bottom (m)
  logical, intent(in), optional :: maskmap(:,:) !< Masks out parallel cores
  logical, intent(in), optional :: fractional_area !< If true, ice_area contains cell area as fraction of entire spherical surface
  ! Local variables
  type(icebergs_gridded), pointer :: grd => null()
  integer :: nbonds
  logical :: check_bond_quality
  integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "ice_bergs: "//trim(version)

  call ice_bergs_framework_init(bergs, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth=ocean_depth, maskmap=maskmap, fractional_area=fractional_area)

  call unit_testing(bergs)

  call mpp_clock_begin(bergs%clock_ior)
  call ice_bergs_io_init(bergs,io_layout)
  call read_restart_calving(bergs)
  if (orig_read) then
    call error_mesg('KID, icebergs_init: ', 'Parameter "orig_read" is no longer supported!', FATAL)
  else
    call read_restart_bergs(bergs,Time)
  endif
  call bergs_chksum(bergs, 'read_restart bergs')
  if (fix_restart_dates) call offset_berg_dates(bergs,Time)
  call mpp_clock_end(bergs%clock_ior)

  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_init, initial status')

  !Reading ocean depth from a file
  if (bergs%read_ocean_depth_from_file) call read_ocean_depth(bergs%grd)

  if (monitor_energy) call energy_tests_init(bergs)
  if (bergs%dem) call init_dem_params(bergs)

  if (bergs%iceberg_bonds_on) then
    if (bergs%manually_initialize_bonds) then
      call initialize_iceberg_bonds(bergs)
    else
      call read_restart_bonds(bergs,Time)
    endif
    call update_halo_icebergs(bergs)
    if (bergs%manually_initialize_bonds) call initialize_iceberg_bonds(bergs)
    if (bergs%mts) then
      call transfer_mts_bergs(bergs)
    else
      call update_halo_icebergs(bergs)
      call connect_all_bonds(bergs)
    endif
    nbonds=0
    check_bond_quality=.True.
    call count_bonds(bergs, nbonds,check_bond_quality)
    call assign_n_bonds(bergs)
    !call fracture_testing_initialization(bergs)
  endif

  if (bergs%uniaxial_test .or. bergs%dem_beam_test>0) call dem_tests_init(bergs)
  if (bergs%damage_test_1) call damage_test_1_init(bergs)

  if (bergs%constant_interaction_LW .and. (bergs%constant_length==0. .or. bergs%constant_width==0.)) then
    call set_constant_interaction_length_and_width(bergs)
  endif
end subroutine icebergs_init

subroutine debugwriteandstop(bergs)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  bergs%debug_write=.true.
  call record_posn(bergs)
  call move_all_trajectories(bergs)
  call write_trajectory(bergs%trajectories, bergs%save_short_traj, bergs%save_fl_traj, bergs%fl_r)
  if (save_bond_traj) call write_bond_trajectory(bergs%bond_trajectories)
  call mpp_sync()
  call error_mesg('KID', 'WRITE AND STOP!!!', FATAL)
end subroutine debugwriteandstop


!> Invoke some unit testing
subroutine unit_testing(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory

  call hexagon_test()
  call point_in_triangle_test()
  call basal_melt_test(bergs)
  call test_check_for_duplicate_ids_in_list()

end subroutine unit_testing

!> Test find_basal_melt()
subroutine basal_melt_test(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  real :: dvo,lat,salt,temp, basal_melt, thickness
  integer(kind=8) :: id
  logical :: Use_three_equation_model

  if (mpp_pe() .eq. mpp_root_pe() ) print *, 'Begining Basal Melting Unit Test'
  dvo=0.2 ;lat=0.0 ; salt=35.0 ; temp=2.0 ;thickness=100.; id=0
  Use_three_equation_model=.False.
  call find_basal_melt(bergs,dvo,lat,salt,temp,Use_three_equation_model,thickness,basal_melt,id)
  if (mpp_pe() .eq. mpp_root_pe()) print *, 'Two equation model basal_melt =',basal_melt

  Use_three_equation_model=.True.
  call find_basal_melt(bergs,dvo,lat,salt,temp,Use_three_equation_model,thickness,basal_melt,id)
   if (mpp_pe() .eq. mpp_root_pe()) print *, 'Three equation model basal_melt =',basal_melt

end subroutine basal_melt_test

!> Test point_in_triangle()
subroutine point_in_triangle_test()
  ! Local variables
  real :: Ax,Ay,Bx,By,Cx,Cy  !Position of icebergs
  logical :: fail_unit_test
  integer :: stderrunit

  ! Get the stderr unit number.
  stderrunit = stderr()
  Ax= -2.695732526092343E-012
  Ay=0.204344508198090
  Bx=-2.695750202346321E-012
  By= -8.433062639672301E-002
  Cx=0.249999999997304
  Cy=6.000694090068343E-002

  fail_unit_test=(.not. point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))
  if (fail_unit_test) call error_mesg('KID, hexagon unit testing:', 'Point in triangle test does not pass!', FATAL)

end subroutine point_in_triangle_test

!> Test Hexagon_into_quadrants_using_triangles()
subroutine hexagon_test()
  ! Local variables
  real :: x0,y0  !Position of icebergs
  real :: H,theta,S !Apothen of iceberg and angle.
  real :: Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4 ! Areas of icebergs
  real :: tol
  logical :: fail_unit_test
  integer :: stderrunit

  ! Get the stderr unit number.
  stderrunit = stderr()

  fail_unit_test=.False.

  tol=1.e-10
  theta=0.0
  H=1.
  S=2.*H/sqrt(3.)

  !Test 1: center at origin: Areas should be equal
  x0=0.  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (abs(Area_hex - ((3.*sqrt(3.)/2.)*(S*S)))>tol) then
    call error_mesg('KID, hexagon unit testing:', 'Hexagon at origin has the wrong area!', WARNING)
    if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    fail_unit_test=.True.
  endif
  if (((abs((Area_hex/4)-Area_Q1 )>tol) .or.  (abs((Area_hex/4)-Area_Q2 )>tol)) .or. ((abs((Area_hex/4)-Area_Q3 )>tol) .or. (abs((Area_hex/4)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon at origin divides into unqual parts!', WARNING)
    fail_unit_test=.True.
  endif

  ! Test 2:  Hexagon split into two quadrants
  !Test 2a: center on x>0 axis
  x0=S  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q1 )>tol) .or.  (abs(0-Area_Q2 )>tol)) .or. ((abs(0-Area_Q3 )>tol) .or. (abs((Area_hex/2)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon split btw 1 and 4!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 2b: center on x<0 axis
  x0=-S  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q2 )>tol) .or.  (abs(0-Area_Q1 )>tol)) .or. ((abs(0-Area_Q4 )>tol) .or. (abs((Area_hex/2)-Area_Q3 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon split btw 2 and 3!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 2c: center on y>0 axis
  x0=0.  ;  y0=H
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q1 )>tol) .or.  (abs(0-Area_Q3 )>tol)) .or. ((abs(0-Area_Q4 )>tol) .or. (abs((Area_hex/2)-Area_Q2 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon split btw 1 and 2!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 3d: center on y<0 axis
  x0=0.  ;  y0=-H
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q3 )>tol) .or.  (abs(0-Area_Q1 )>tol)) .or. ((abs(0-Area_Q2 )>tol) .or. (abs((Area_hex/2)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon split btw 3 and 4!', WARNING)
    fail_unit_test=.True.
  endif

  ! Test 3:  Two corners of hex on the axis
  !Test 3a: center on x>0 axis
  x0=S/2.  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((2.5*Area_hex/6.)-Area_Q1 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q2 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q3 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon split two coners of hex (x>0)!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 3b: center on x<0 axis
  x0=-S/2.  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((2.5*Area_hex/6.)-Area_Q2 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q1 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q4 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q3 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('KID, hexagon unit testing:', 'Hexagon split two coners of hex (x<0)!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 3c: center on y>0 axis
  !x0=0.  ;  y0=H/2.
  !call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  !if (((abs((2.5*Area_hex/6.)-Area_Q1 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q3 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q4 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q2 )>tol))) then
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon errors =', (abs((2.5*Area_hex/6.)-Area_Q1 )), (abs((0.5*Area_hex/6.)-Area_Q3 )),&
  !  call error_mesg('KID, hexagon unit testing:', 'Hexagon split two coners of hex (y>0)!', WARNING)
  !  fail_unit_test=.True.
  !endif
  !!Test 3d: center on y<0 axis
  !x0=0.  ;  y0=-H/2.
  !call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  !if (((abs((2.5*Area_hex/6.)-Area_Q3 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q2 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q1 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q4 )>tol))) then
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'KID, hexagon errots =', (abs((2.5*Area_hex/6.)-Area_Q3 )), (abs((0.5*Area_hex/6.)-Area_Q2 )),&
  !  call error_mesg('KID, hexagon unit testing:', 'Hexagon split two coners of hex (y<0)!', WARNING)
  !  fail_unit_test=.True.
  !endif


  if (fail_unit_test) call error_mesg('KID, hexagon unit testing:', 'Hexagon unit testing does not pass!', FATAL)

end subroutine hexagon_test

!> Initializes bonds
subroutine initialize_iceberg_bonds(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(iceberg), pointer :: berg
  type(iceberg), pointer :: other_berg
  type(icebergs_gridded), pointer :: grd
  type(bond) , pointer :: current_bond
  logical :: already_bonded
  real :: T1, L1, W1, lon1, lat1, x1, y1, R1, A1   !Current iceberg
  real :: T2, L2, W2, lon2, lat2, x2, y2, R2, A2   !Other iceberg
  real :: dlon,dlat
  real :: dx_dlon,dy_dlat, lat_ref
  real :: r_dist_x, r_dist_y, r_dist
  real :: radius1,radius2,rdenom
  integer :: grdi_outer, grdj_outer
  integer :: grdi_inner, grdj_inner

  if (bergs%manually_initialize_bonds_from_radii) then
    if (bergs%hexagonal_icebergs) then
      rdenom=1./(2.*sqrt(3.))
    else
      !rdenom=1./pi
      rdenom=1./4.
    endif
  endif

  ! For convenience
  grd=>bergs%grd
  !Should update halos before doing this
  ! do grdj_outer = grd%jsc,grd%jec ; do grdi_outer = grd%isc,grd%iec  !Should you be on the data domain??
  do grdj_outer = grd%jsd,grd%jed ; do grdi_outer = grd%isd,grd%ied  !using data domain -Alex
    berg=>bergs%list(grdi_outer,grdj_outer)%first
    do while (associated(berg)) ! loop over all bergs

      lon1=berg%lon; lat1=berg%lat
      !call rotpos_to_tang(lon1,lat1,x1,y1)  !Is this correct? Shouldn't it only be on tangent plane?

      ! do grdj_inner = grd%jsc,grd%jec ; do grdi_inner = grd%isc,grd%iec  !This line uses n^2 steps
      do grdj_inner = grd%jsd,grd%jed ; do grdi_inner = grd%isd,grd%ied !Uses n^2 steps. Change to data domain-Alex
!     do grdj_inner = berg%jne-1,berg%jne+1 ; do grdi_inner = berg%ine-1,berg%ine+1   !Only looping through adjacent cells.
        other_berg=>bergs%list(grdi_inner,grdj_inner)%first
        do while (associated(other_berg)) ! loop over all other bergs

          if (berg%id .ne. other_berg%id) then
            !first, make sure the bergs are not bonded already
            already_bonded=.false.
            current_bond=>berg%first_bond
            do while (associated(current_bond))
              if (current_bond%other_id .ne. other_berg%id) then
                current_bond=>current_bond%next_bond
              else
                current_bond=>null()
                already_bonded=.true.
              endif
            enddo

            if (.not. already_bonded) then
              lon2=other_berg%lon; lat2=other_berg%lat
              dlon=lon1-lon2;      dlat=lat1-lat2
              lat_ref=0.5*(lat1+lat2)
              call convert_from_grid_to_meters(lat_ref,grd%grid_is_latlon,dx_dlon,dy_dlat)
              r_dist_x=dlon*dx_dlon
              r_dist_y=dlat*dy_dlat
              r_dist=sqrt( (r_dist_x**2) + (r_dist_y**2) )

              if (bergs%manually_initialize_bonds_from_radii) then
                radius1=sqrt(berg%length*berg%width*rdenom)
                radius2=sqrt(other_berg%length*other_berg%width*rdenom)
                !radius=sqrt(min(berg%length*berg%width,other_berg%length*other_berg%width))
                if (r_dist.lt.1.1*(radius1+radius2)) &
                  call form_a_bond(berg, other_berg%id, other_berg%ine, other_berg%jne, other_berg)
              elseif (r_dist.lt.bergs%length_for_manually_initialize_bonds) then
                ! If the bergs are closer than bergs%length_for_manually_initialize_bonds, then form a bond -Alex
                call form_a_bond(berg, other_berg%id, other_berg%ine, other_berg%jne, other_berg)
              endif
            endif
          endif
          other_berg=>other_berg%next
        enddo  ! End of looping through all other bergs in the inner list
      enddo ; enddo;  !End of inner loop
      berg=>berg%next
    enddo ! End of looping through all bergs in the outer list
  enddo ; enddo; !End of outer loop.

end subroutine initialize_iceberg_bonds

!> Returns metric converting grid distances to meters
subroutine convert_from_grid_to_meters(lat_ref, grid_is_latlon, dx_dlon, dy_dlat)
  ! Arguments
  real, intent(in) :: lat_ref !< Latitude at which to make metric conversion (degree N)
  logical, intent(in) :: grid_is_latlon !< True if grid model grid is in lat-lon coordinates
  real, intent(out) :: dx_dlon !< Metric dx/dlon
  real, intent(out) :: dy_dlat !< Metric dy/dlat

  if (grid_is_latlon) then
    dx_dlon=(pi/180.)*Rearth*cos((lat_ref)*(pi/180.))
    dy_dlat=(pi/180.)*Rearth
  else
    dx_dlon=1.
    dy_dlat=1.
  endif

end subroutine convert_from_grid_to_meters

!> Returns metric converting distance in meters to grid distance
subroutine  convert_from_meters_to_grid(lat_ref,grid_is_latlon ,dlon_dx,dlat_dy)
  ! Arguments
  real, intent(in) :: lat_ref !< Latitude at which to make metric conversion (degree N)
  logical, intent(in) :: grid_is_latlon !< True if grid model grid is in lat-lon coordinates
  real, intent(out) :: dlon_dx !< Metric dlon/dx
  real, intent(out) :: dlat_dy !< Metric dlat/dy

  if (grid_is_latlon) then
    dlon_dx=(180./pi)/(Rearth*cos((lat_ref)*(pi/180.)))
    dlat_dy=(180./pi)/Rearth
  else
    dlon_dx=1.
    dlat_dy=1.
  endif

end subroutine convert_from_meters_to_grid

!> Calculates interactions between a berg and all bergs in range
subroutine interactive_force(bergs, berg, IA_x, IA_y, u0, v0, u1, v1,&
                             P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Primary iceberg
  type(iceberg), pointer :: other_berg !< Berg that primary is interacting with
  real, intent(in) :: u0 !< Zonal velocity of primary berg (m/s)
  real, intent(in) :: v0 !< Meridional velocity of primary berg (m/s)
  real, intent(in) :: u1 !< Zonal velocity of other berg (m/s)
  real, intent(in) :: v1 !< Meridional velocity of other berg (m/s)
  real, intent(out) :: IA_x !< Net zonal acceleration of berg due to interactions (m/s2)
  real, intent(out) :: IA_y !< Net meridional acceleration of berg due to interactions (m/s2)
  real, intent(out) :: P_ia_11 !< Damping projection matrix, xx component (kg/s)
  real, intent(out) :: P_ia_12 !< Damping projection matrix, xy component (kg/s)
  real, intent(out) :: P_ia_22 !< Damping projection matrix, yy component (kg/s)
  real, intent(out) :: P_ia_21 !< Damping projection matrix, yx component (kg/s)
  real, intent(out) :: P_ia_times_u_x !< Zonal damping projection (without new berg velocity) (kg m/s)
  real, intent(out) :: P_ia_times_u_y !< Meridional damping projection (without new berg velocity) (kg m/s)
  ! Local variables
  type(bond), pointer :: current_bond
  type(icebergs_gridded), pointer :: grd
  real :: u2, v2
  logical :: critical_interaction_damping_on
  integer :: grdi, grdj
  logical :: iceberg_bonds_on
  logical :: bonded
  integer :: nc_x,nc_y

  IA_x=0.; IA_y=0.
  P_ia_11=0. ; P_ia_12=0. ;  P_ia_21=0.;  P_ia_22=0.
  P_ia_times_u_x=0. ; P_ia_times_u_y=0.

  !no interactive force allowed for footloose child bergs until they are out of contact range of any of other berg
  !at least once.
  if (berg%fl_k.eq.-1) return

  nc_x=bergs%contact_cells_lon; nc_y=bergs%contact_cells_lat
  iceberg_bonds_on=bergs%iceberg_bonds_on

  if (bergs%mts .or. (bergs%contact_distance>0.) .or. (bergs%contact_spring_coef .ne. bergs%spring_coef) ) then

    ! Interactions for bonded berg elements only. For the mts scheme, this is only called if mts_part==3
    if ( (.not. bergs%mts) .or. (bergs%mts .and. bergs%mts_part .eq. 3)) then
      bonded=.true.
      if (iceberg_bonds_on) then
        current_bond=>berg%first_bond
        do while (associated(current_bond)) ! loop over all bonds
          other_berg=>current_bond%other_berg
          if (.not. associated(other_berg)) then
            call error_mesg('KID,bond interactions', 'Trying to do Bond interactions with unassosiated berg!' ,FATAL)
          else
            call calculate_force(bergs, berg, other_berg, IA_x, IA_y, u0, v0, u1, v1,  &
              P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded)
            other_berg%id=-other_berg%id !mark, as to not repeat in contact search below
          endif
          current_bond=>current_bond%next_bond
        enddo
        !contact search: find any surrounding non-bonded bergs in the same conglom, and processes
        !them for contact using the regular (non-contact) spring constant and crit_dist based on radii
        bonded=.false.
        do grdj = max(berg%jne-2,bergs%grd%jsd+1),min(berg%jne+2,bergs%grd%jed);&
          do grdi = max(berg%ine-2,bergs%grd%isd+1),min(berg%ine+2,bergs%grd%ied)
          other_berg=>bergs%list(grdi,grdj)%first
          do while (associated(other_berg))
            if (other_berg%id>0 .and. other_berg%conglom_id.eq.berg%conglom_id) then
              call calculate_force(bergs, berg, other_berg, IA_x, IA_y, u0, v0, u1, v1,  &
                P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded,c_crit_dist=.true.)
            endif
            other_berg=>other_berg%next
          enddo
        enddo;enddo
        !unmark the bonded bergs
        current_bond=>berg%first_bond
        do while (associated(current_bond))
          other_berg=>current_bond%other_berg
          other_berg%id=abs(other_berg%id)
          current_bond=>current_bond%next_bond
        enddo
      endif
    endif

    ! Interactions for non-bonded berg elements only. For the mts scheme, this is only called if mts_part==1
    if (.not. (bergs%mts .and. bergs%mts_part .eq. 3)) then
      bonded=.false. ! Interactions for non-bonded berg elements
      do grdj = max(berg%jne-nc_y,bergs%grd%jsd),min(berg%jne+nc_y,bergs%grd%jed);&
        do grdi = max(berg%ine-nc_x,bergs%grd%isd),min(berg%ine+nc_x,bergs%grd%ied)

        other_berg=>bergs%list(grdi,grdj)%first
        do while (associated(other_berg)) ! loop over all other bergs
          if (other_berg%conglom_id.ne.berg%conglom_id) then
            call calculate_force(bergs, berg, other_berg, IA_x, IA_y, u0, v0, u1, v1,  &
              P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded)
          endif
          other_berg=>other_berg%next
        enddo
      enddo; enddo
      if (bergs%use_spring_for_land_contact) then
        grd=>bergs%grd
        do grdj = max(berg%jne-nc_y,bergs%grd%jsd),min(berg%jne+nc_y,bergs%grd%jed);&
          do grdi = max(berg%ine-nc_x,bergs%grd%isd),min(berg%ine+nc_x,bergs%grd%ied)
          if (grd%msk(grdi,grdj)==0) call calculate_force_land_contact(bergs, berg, grd, grdi, grdj, &
            IA_x, IA_y, P_ia_11, P_ia_12, P_ia_21, P_ia_22) !, P_ia_times_u_x, P_ia_times_u_y)
        enddo; enddo
      endif
    endif

  else
    !Alon's original code, only usuable for simulations WITHOUT MTS Verlet, a specified contact_distance,
    !or a separate spring constant for collision.
    bonded=.false. ! 'Unbonded' iceberg interactions (here, not to be taken literally)
    do grdj = berg%jne-1,berg%jne+1 ; do grdi = berg%ine-1,berg%ine+1
      other_berg=>bergs%list(grdi,grdj)%first
      do while (associated(other_berg)) ! loop over all other bergs
        call calculate_force(bergs, berg, other_berg, IA_x, IA_y, u0, v0, u1, v1, &
          P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y, bonded)
        other_berg=>other_berg%next
      enddo ! loop over all bergs
    enddo; enddo

    bonded=.true. ! Interactions due to iceberg bonds
    if (iceberg_bonds_on) then ! MP1
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        other_berg=>current_bond%other_berg
        if (.not. associated(other_berg)) then
          call error_mesg('KID,bond interactions', 'Trying to do Bond interactions with unassosiated berg!' ,FATAL)
        else
          call calculate_force(bergs, berg, other_berg, IA_x, IA_y, u0, v0, u1, v1,  &
            P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded)
        endif
        current_bond=>current_bond%next_bond
      enddo
    endif
    if (bergs%use_spring_for_land_contact) then
      grd=>bergs%grd
      do grdj = max(berg%jne-nc_y,bergs%grd%jsd),min(berg%jne+nc_y,bergs%grd%jed);&
        do grdi = max(berg%ine-nc_x,bergs%grd%isd),min(berg%ine+nc_x,bergs%grd%ied)
        if (grd%msk(grdi,grdj)==0) call calculate_force_land_contact(bergs, berg, grd, grdi, grdj, &
          IA_x, IA_y, P_ia_11, P_ia_12, P_ia_21, P_ia_22) !, P_ia_times_u_x, P_ia_times_u_y)
      enddo; enddo
    endif
  endif
end subroutine interactive_force


!> Calculate interactive forces between two bergs
subroutine calculate_force(bergs, berg, other_berg, IA_x, IA_y, u0, v0, u1, v1, &
  P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y, bonded, c_crit_dist)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Primary berg
  type(iceberg), pointer :: other_berg !< Berg that primary is interacting with
  real, intent(inout) :: IA_x !< Net zonal acceleration of berg due to interactions (m/s2)
  real, intent(inout) :: IA_y !< Net meridional acceleration of berg due to interactions (m/s2)
  real, intent(in) :: u0 !< Zonal velocity of primary berg (m/s)
  real, intent(in) :: v0 !< Meridional velocity of primary berg (m/s)
  real, intent(in) :: u1 !< Zonal velocity of other berg (m/s)
  real, intent(in) :: v1 !< Meridional velocity of other berg (m/s)
  real, intent(inout) :: P_ia_11 !< Damping projection matrix, xx component (kg/s)
  real, intent(inout) :: P_ia_12 !< Damping projection matrix, xy component (kg/s)
  real, intent(inout) :: P_ia_22 !< Damping projection matrix, yy component (kg/s)
  real, intent(inout) :: P_ia_21 !< Damping projection matrix, yx component (kg/s)
  real, intent(inout) :: P_ia_times_u_x !< Zonal damping projection (without new berg velocity) (kg m/s)
  real, intent(inout) :: P_ia_times_u_y !< Meridional damping projection (without new berg velocity) (kg m/s)
  logical ,intent(in) :: bonded !< If T, berg and other_berg are bonded
  logical ,intent(in), optional :: c_crit_dist !< Modified critical interaction distance for contact (m)
  ! Local variables
  real :: T1, L1, W1, lon1, lat1, x1, y1, R1, A1 ! Current iceberg
  real :: T2, L2, W2, lon2, lat2, x2, y2, R2, A2 ! Other iceberg
  real :: dlon, dlat
  real :: r_dist_x, r_dist_y, r_dist, A_o, A_min, trapped, T_min
  real :: P_11, P_12, P_21, P_22
  real :: M1, M2, M_min
  real :: u2, v2
  real :: lat_ref, dx_dlon, dy_dlat
  logical :: critical_interaction_damping_on,tbonded,crit_dist_from_radius_only
  real :: spring_coef, accel_spring, radial_damping_coef, p_ia_coef, tangental_damping_coef, bond_coef
  real :: crit_dist, vr1,vr2,vn,vn1,vn2,vt1,vt2

  if ((berg%id .ne. other_berg%id) .and. (berg%fl_k.ne.-1) .and. (other_berg%fl_k.ne.-1)) then
    ! From Berg 1
    T1=berg%thickness
    lon1=berg%lon_old; lat1=berg%lat_old
    !call rotpos_to_tang(lon1,lat1,x1,y1)

    ! From Berg 2
    T2=other_berg%thickness
    lon2=other_berg%lon_old; lat2=other_berg%lat_old !Old values are used to make it order invariant

    !call rotpos_to_tang(lon2,lat2,x2,y2)
    u2=other_berg%uvel_old;  v2=other_berg%vvel_old !Old values are used to make it order invariant

    if (bergs%constant_interaction_LW .and. bergs%mts .and. bonded) then
      ! use constant length and width here, just for interactions
      A1=bergs%constant_length*bergs%constant_width !Berg 1 area
      M1=A1*T1*bergs%rho_bergs                      !Berg 1 mass
      A2=A1                                         !Berg 2 area
      M2=A2*T2*bergs%rho_bergs                      !Berg 2 mass
    else
      ! use actual length and width of bergs
      L1=berg%length; W1=berg%width; M1=berg%mass; A1=L1*W1                   !From Berg 1
      L2=other_berg%length; W2=other_berg%width; M2=other_berg%mass; A2=L2*W2 !From Berg 2
    endif

    dlon=lon1-lon2
    dlat=lat1-lat2

    ! Note that this is not the exact distance along a great circle.
    ! Approximation for small distances. Should be fine.
    !r_dist_x=x1-x2 ; r_dist_y=y1-y2
    !r_dist=sqrt( ((x1-x2)**2) + ((y1-y2)**2) )
    lat_ref=0.5*(lat1+lat2)
    call convert_from_grid_to_meters(lat_ref,bergs%grd%grid_is_latlon,dx_dlon,dy_dlat)

    r_dist_x=dlon*dx_dlon
    r_dist_y=dlat*dy_dlat
    r_dist=sqrt( (r_dist_x**2) + (r_dist_y**2) ) !(Stern et al 2017, Eqn 3)

    !Stern et al 2017, Eqn 4: radius of circle inscribed within hexagon or square with area A1 or A2
    if (bergs%hexagonal_icebergs) then
      R1=sqrt(A1/(2.*sqrt(3.)))
      R2=sqrt(A2/(2.*sqrt(3.)))
    else !square packing
      if (bergs%iceberg_bonds_on) then
        R1=0.5*sqrt(A1)
        R2=0.5*sqrt(A2)
      else
        R1=sqrt(A1/pi) ! Interaction radius of the iceberg (assuming circular icebergs)
        R2=sqrt(A2/pi) ! Interaction radius of the other iceberg
      endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!debugging!!!!!!!!!!!!!!!!!!!!!!!!!!MP1
    ! if (berg%id .eq. 1) then
    !   print *, 'Comparing longitudes: ', lon1, lon2, r_dist_x, dlon
    !   print *, 'Comparing latitudes: ', lat1, lat2, r_dist_y, dlat
    !   print *, 'Outside, id, r_dist', berg%id, r_dist,bonded
    !   print *, 'Halo_status', berg%halo_berg,other_berg%halo_berg
    ! endif
    ! print *, 'outside the loop',R1, R2,r_dist, bonded
    !!!!!!!!!!!!!!!!!!!!!!!!!!!debugging!!!!!!!!!!!!!!!!!!!!!!!!!!


    !call overlap_area(R1,R2,r_dist,A_o,trapped)
    !T_min=min(T1,T2)
    !A_min = min((pi*R1**R1),(pi*R2*R2))
    M_min=min(M1,M2)
    !Calculating spring force (Stern et al 2017, Eqn 6):
    if (bonded) then
      crit_dist=R1+R2 !Stern et al 2017, Eqn 5
      spring_coef=bergs%spring_coef
    else
      spring_coef=bergs%contact_spring_coef
      if (present(c_crit_dist)) then
        if (c_crit_dist) then
          crit_dist=R1+R2
          spring_coef=bergs%spring_coef
        else
          crit_dist=max(R1+R2,bergs%contact_distance)
        endif
      else
        crit_dist=max(R1+R2,bergs%contact_distance)
      endif
    end if

    radial_damping_coef=bergs%radial_damping_coef
    tangental_damping_coef=bergs%tangental_damping_coef
    critical_interaction_damping_on=bergs%critical_interaction_damping_on

    ! Using critical values for damping rather than manually setting the damping.
    if (critical_interaction_damping_on) then
      radial_damping_coef=2.*sqrt(spring_coef) ! Critical damping
      if (bergs%tang_crit_int_damp_on) then
        tangental_damping_coef=(2.*sqrt(spring_coef))/4 ! Critical damping (just a guess)
      endif
    endif

    tbonded=bonded
    !if bonded and  STS,contact_dist=0, and  contact_spring coeff==spring coeff
    if (bonded .and. .not. (bergs%mts .or. (bergs%contact_distance>0.) .or. &
      (bergs%contact_spring_coef .ne. bergs%spring_coef) )) then
      if (.not. (r_dist>crit_dist)) tbonded=.false.
    endif

    if (bergs%uniaxial_test) then
      if (berg%start_lon==bergs%dem_tests_start_lon) IA_x=IA_x-1.e-5
      if (berg%start_lon==bergs%dem_tests_end_lon)   IA_x=IA_x+1.e-5
    endif

    if  ((r_dist>0.) .and. ( tbonded .or. (r_dist<crit_dist .and. .not. bonded) )) then

      !Spring force (Stern et al 2017, Eqn 7):
      !accel_spring=spring_coef*(T_min/T1)*(A_o/A1) ! Old version dependent on area
      accel_spring=spring_coef*(M_min/M1)*(crit_dist-r_dist)
      IA_x=IA_x+(accel_spring*(r_dist_x/r_dist))
      IA_y=IA_y+(accel_spring*(r_dist_y/r_dist))

      ! if (r_dist < 5*(R1+R2)) then
        ! Damping force (Stern et al 2017, Eqn 8):
        ! Paralel velocity
        !projection matrix in Stern et al 2017, Eqn 8:
        P_11=(r_dist_x*r_dist_x)/(r_dist**2)
        P_12=(r_dist_x*r_dist_y)/(r_dist**2)
        P_21=(r_dist_x*r_dist_y)/(r_dist**2)
        P_22=(r_dist_y*r_dist_y)/(r_dist**2)
        !p_ia_coef=radial_damping_coef*(T_min/T1)*(A_min/A1)
        p_ia_coef=radial_damping_coef*(M_min/M1)
        !the following was done originally, but it is wrong:
        if (bergs%scale_damping_by_pmag) then
          p_ia_coef=p_ia_coef*(0.5*(sqrt((((P_11*(u2-u1))+(P_12*(v2-v1)))**2)+ (((P_12*(u2-u1))+(P_22*(v2-v1)))**2)) &
            + sqrt((((P_11*(u2-u0))+(P_12*(v2-v0)))**2)+(((P_12*(u2-u0)) +(P_22*(v2-v0)))**2))))
        endif

        P_ia_11=P_ia_11+p_ia_coef*P_11
        P_ia_12=P_ia_12+p_ia_coef*P_12
        P_ia_21=P_ia_21+p_ia_coef*P_21
        P_ia_22=P_ia_22+p_ia_coef*P_22
        P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
        P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))
        !print *, 'Paralel: ',berg%id,  p_ia_coef, IA_x, P_ia_11, P_ia_21,P_ia_12, P_ia_22

        ! Normal velocities
        P_11=1-P_11  ;  P_12=-P_12 ; P_21= -P_21 ;    P_22=1-P_22
        !p_ia_coef=tangental_damping_coef*(T_min/T1)*(A_min/A1)
        p_ia_coef=tangental_damping_coef*(M_min/M1)
        !the following was done originally, but it is wrong:
        if (bergs%scale_damping_by_pmag) then
          p_ia_coef=p_ia_coef*(0.5*(sqrt((((P_11*(u2-u1))+(P_12*(v2-v1)))**2)+ (((P_12*(u2-u1))+(P_22*(v2-v1)))**2))  &
            + sqrt((((P_11*(u2-u0))+(P_12*(v2-v0)))**2)+(((P_12*(u2-u0)) +(P_22*(v2-v0)))**2))))
        endif
        P_ia_11=P_ia_11+p_ia_coef*P_11
        P_ia_12=P_ia_12+p_ia_coef*P_12
        P_ia_21=P_ia_21+p_ia_coef*P_21
        P_ia_22=P_ia_22+p_ia_coef*P_22
        P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
        P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))
        !print *, 'Perp: ',berg%id,  p_ia_coef, IA_x, P_ia_11, P_ia_21,P_ia_12, P_ia_22
        !print *, 'P_11',P_11
        !print *, 'P_21',P_21
        !print *, 'P_12',P_12
        !print *, 'P_22',P_22
      !endif
    endif
  endif

end subroutine calculate_force

!> Calculate interactive forces between two bergs using a discrete element method (DEM)
!> that includes additional terms for tangential forces and torque.
subroutine calculate_force_dem(bergs, berg, other_berg, current_bond, &
  dt, F_x, F_y, T, Fd_x, Fd_y, T_d, savestress)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Primary berg
  type(iceberg), pointer :: other_berg !< Berg that primary is interacting with
  type(bond), pointer :: current_bond !< Bond being processed
  real, intent(in) :: dt !< Time step size
  real, intent(inout) :: F_x !< Net x-force without damping (N)
  real, intent(inout) :: F_y !< Net y-force without damping (N)
  real, intent(inout) :: T !< Net torque without damping (Nm)
  real, intent(inout) :: Fd_x !< Net x-force from damping (N)
  real, intent(inout) :: Fd_y !< Net y-force from damping (N)
  real, intent(inout) :: T_d !< Net torque from damping (Nm)
  logical, intent(in) :: savestress !< Save stress/tangential displacement on bond or not
  ! Local variables
  real :: A1,T1,M1,u1,v1,lon1,lat1,R1 ! Current iceberg
  real :: A2,T2,M2,u2,v2,lon2,lat2,R2 ! Other iceberg
  real :: dlon, dlat, lat_ref, dx_dlon, dy_dlat
  real :: r_dist_x, r_dist_y, r_dist
  real :: Rmin, delta, L, l0, youngs_dem, Thick, T_Rmin, n1, n2
  real :: ur, vr, up, vp, rotu, rotv, rtu, rtv, RR1, RR2
  real :: ur2, vr2
  real :: ss_factor, ss1, ss2
  real :: K_damp, Meff, damping_coef
  real :: Fn_x, Fn_y, Fs_x, Fs_y, Ts, Tr
  real :: tangdotnt,tangd1p,tangd2p, tmag, tmagp
  real :: RR1x,RR1y,RR2x,RR2y

  !This DEM formulation follows Wang 2020 "A scale-invariant bonded particle model for
  !simulating large deformation and failure of continua", which differs from "classic"
  !DEM (e.g. Potyondy & Cundall, 2004) in its definition of bond width and torque from
  !relative particle rotation. Forces and torques were further modified here to depend
  !on ice thickness.

  if ((berg%id .ne. other_berg%id) .and. (berg%fl_k.ne.-1) .and. (other_berg%fl_k.ne.-1)) then

    ! From Berg 1
    T1=berg%thickness
    lon1=berg%lon_old; lat1=berg%lat_old

    ! From Berg 2
    T2=other_berg%thickness
    lon2=other_berg%lon_old; lat2=other_berg%lat_old

    if (bergs%constant_interaction_LW) then
      ! use constant length and width here, just for interactions
      A1=bergs%constant_length*bergs%constant_width !Berg 1 area
      M1=A1*T1*bergs%rho_bergs                      !Berg 1 mass
      A2=A1                                         !Berg 2 area
      M2=A2*T2*bergs%rho_bergs                      !Berg 2 mass
    else
      ! use actual length and width of bergs
      M1=berg%mass; A1=berg%length*berg%width                   !From Berg 1
      M2=other_berg%mass; A2=other_berg%length*other_berg%width !From Berg 2
    endif

    dlon=lon1-lon2; dlat=lat1-lat2

    ! Note that this is not the exact distance along a great circle.
    ! Approximation for small distances. Should be fine.
    lat_ref=0.5*(lat1+lat2)
    call convert_from_grid_to_meters(lat_ref,bergs%grd%grid_is_latlon,dx_dlon,dy_dlat)

    r_dist_x=dlon*dx_dlon
    r_dist_y=dlat*dy_dlat
    r_dist=sqrt( (r_dist_x**2) + (r_dist_y**2) ) !(Stern et al 2017, Eqn 3)
    current_bond%length=r_dist

    !unit vec
    n1=r_dist_x/r_dist; n2=r_dist_y/r_dist

    !Stern et al 2017, Eqn 4: radius of circle inscribed within hexagon or square with area A1 or A2
    if (bergs%hexagonal_icebergs) then
      R1=sqrt(A1/(2.*sqrt(3.)))
      R2=sqrt(A2/(2.*sqrt(3.)))
    else !square packing
      if (bergs%iceberg_bonds_on) then
        R1=0.5*sqrt(A1)
        R2=0.5*sqrt(A2)
      else
        R1=sqrt(A1/pi) ! Interaction radius of the iceberg (assuming circular icebergs)
        R2=sqrt(A2/pi) ! Interaction radius of the other iceberg
      endif
    endif

    !Rmin=min(R1,R2) !old
    if (R1<R2) then
      Rmin=R1
      T_Rmin=T1 !thickness of element w/ minimum radius
    else
      Rmin=R2
      T_Rmin=T2
    endif
    l0=R1+R2
    delta=l0-r_dist

    !element distances to contact point
    RR1=R1-0.5*delta; RR2=R2-0.5*delta !magnitude
    RR1x=RR1*n1;  RR1y=RR1*n2 !element 1 components
    RR2x=RR2*n1;  RR2y=RR2*n2 !element 2 components

    !bond width (determined at contact point)
    L=2.0*(Rmin+(Rmin-0.5*delta)*abs(R1-R2)/r_dist)

    Thick=T_Rmin+(Rmin-0.5*delta)*abs(T1-T2)/r_dist !thickness as determined at contact point
    !Thick=min(T1,T2) !minimum thickness
    !Thick=0.5*(T1+T2) !average thickness

    youngs_dem=bergs%dem_spring_coef

    !normal force:
    !Fn_x=youngs_dem*Thick*delta*n1*L/l0 !Fn_y=youngs_dem*Thick*delta*n2*L/l0
    !normal stiffness is youngs_dem*Thick*L/l0
    if (savestress) current_bond%nstress=-youngs_dem*delta/l0 !normal stress
    Fn_x=youngs_dem*Thick*delta*L/l0; Fn_y=Fn_x*n2; Fn_x=Fn_x*n1

    !relative translational velocity
    ur=berg%uvel_old-other_berg%uvel_old; vr=berg%vvel_old-other_berg%vvel_old

    if (savestress) then
      !Rotation of contact plane:
      tmag=sqrt(current_bond%tangd1**2+current_bond%tangd2**2) !old magnitude
      !project old tangential displacement to current tangent plane
      tangdotnt=current_bond%tangd1*n1+current_bond%tangd2*n2
      tangd1p=current_bond%tangd1-tangdotnt*n1; tangd2p=current_bond%tangd2-tangdotnt*n2
      !rescale to the old magnitude
      tmagp=sqrt(tangd1p**2+tangd2p**2)
      if (tmagp>0.) then
        tangd1p=(tmag/tmagp)*tangd1p; tangd2p=(tmag/tmagp)*tangd2p
      else
        tangd1p=0.; tangd2p=0.
      endif
      current_bond%tangd1=tangd1p; current_bond%tangd2=tangd2p

      !relative tangential velocities from rotation
      rotu=RR1y*berg%ang_vel + RR2y*other_berg%ang_vel
      rotv=-(RR1x*berg%ang_vel + RR2x*other_berg%ang_vel)
      ur2=ur+rotu; vr2=vr+rotv

      !'parallel' velocity = dot(relative velocity,unit vec)*unit vec
      up=ur2*n1+vr2*n2; vp=up*n2; up=up*n1
      rtu=ur2-up; rtv=vr2-vp !the relative tangential velocities

      !add the new tangential displacement
      current_bond%tangd1=current_bond%tangd1+rtu*dt; current_bond%tangd2=current_bond%tangd2+rtv*dt
    endif

    !shear terms
    ss_factor=-L*Thick*youngs_dem/(l0*2.0*(1.0+bergs%poisson)) !shear stiffness
    if (bergs%ignore_tangential_force) ss_factor=0.
    Fs_x=ss_factor*current_bond%tangd1; Fs_y=ss_factor*current_bond%tangd2 !shear forces
    if (savestress) current_bond%sstress=sqrt(Fs_x**2+Fs_y**2)/(L*Thick) !shear stress

    !Torque from shearing force: Ts = R .cross. Fs
    Ts = -(RR1x*Fs_y-RR1y*Fs_x)

    !Torque from relative particle rotation
    Tr = -youngs_dem*Thick*(L*L*L)*sin(berg%rot-other_berg%rot)/(12.*l0)

    !damping
    K_damp = 2.*youngs_dem/(3.*(1.-bergs%poisson**2)) !damping factor
    Meff = M1*M2/(M1+M2) !effective mass

    damping_coef = bergs%dem_damping_coef*sqrt(K_damp*Meff)

    Fd_x = Fd_x-damping_coef*ur; Fd_y = Fd_y-damping_coef*vr  !linear damping force
    T_d  = T_d-damping_coef*(berg%ang_vel-other_berg%ang_vel) !damping torque

    if (bergs%dem_shear_for_frac_only) then
      Fs_x=0.; Fs_y=0.
    endif

    !final forces
    T = T + (Ts + Tr) !torque
    F_x = F_x + Fn_x + Fs_x; F_y = F_y + Fn_y + Fs_y !linear forces
  endif

end subroutine calculate_force_dem

!> Experimental subroutine to calculate interactive force between a berg and a land cell
!> Should prevent grounding.
subroutine calculate_force_land_contact(bergs, berg, grd, i, j, IA_x, IA_y, &
  P_ia_11, P_ia_12, P_ia_21, P_ia_22) !, P_ia_times_u_x, P_ia_times_u_y)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Primary berg
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  integer, intent(in) :: i !< Grid cell i to contact
  integer, intent(in) :: j !< Grid cell j to contact
  real, intent(inout) :: IA_x !< Net zonal acceleration of berg due to interactions (m/s2)
  real, intent(inout) :: IA_y !< Net meridional acceleration of berg due to interactions (m/s2)
  real, intent(inout) :: P_ia_11 !< Damping projection matrix, xx component (kg/s)
  real, intent(inout) :: P_ia_12 !< Damping projection matrix, xy component (kg/s)
  real, intent(inout) :: P_ia_22 !< Damping projection matrix, yy component (kg/s)
  real, intent(inout) :: P_ia_21 !< Damping projection matrix, yx component (kg/s)
  !not needed because "other" iceberg is assumed to have velocity = 0
  !real, intent(inout) :: P_ia_times_u_x !< Zonal damping projection (without new berg velocity) (kg m/s)
  !real, intent(inout) :: P_ia_times_u_y !< Meridional damping projection (without new berg velocity) (kg m/s)

  ! Local variables
  real :: T1, L1, W1, lon1, lat1, R1, A1 ! Current iceberg
  real :: lon2, lat2 ! "Other" iceberg
  real :: dlon, dlat
  real :: r_dist_x, r_dist_y, r_dist
  real :: P_11, P_12, P_21, P_22
  !real :: u2, v2
  real :: lat_ref, dx_dlon, dy_dlat, crit_dist
  logical :: critical_interaction_damping_on
  real :: spring_coef, accel_spring, radial_damping_coef, p_ia_coef, tangental_damping_coef

  if (i == berg%ine .and. j == berg%jne) return

  lon1=berg%lon_old; lat1=berg%lat_old  ! From Berg 1

  ! Find the location on the edge/corner of the land cell that is closest to the berg. The force is
  ! calculated  as if there is a second berg at this location, with identical size to the input berg,
  ! but with velocities set to zero
  !u2=0;  v2=0 !"other berg" velocities are assumed to be zero
  if (berg%ine==i) then
    lon2=lon1
    if (berg%jne>j) then !N of cell
      lat2=grd%lat(i,j)
    else !S of cell
      lat2=grd%lat(i,j-1)
    endif
  elseif (berg%ine>i) then
    if (berg%jne==j) then !E of cell
      lon2=grd%lon(i,j)
      lat2=lat1
    elseif (berg%jne>j) then !NE of cell
      lon2=grd%lon(i,j)
      lat2=grd%lat(i,j)
    else !SE of cell
      lon2=grd%lon(i,j-1)
      lat2=grd%lat(i,j-1)
    endif
  else
    if (berg%jne==j) then !W of cell
      lon2=grd%lon(i-1,j)
      lat2=lat1
    elseif (berg%jne>j) then !NW of cell
      lon2=grd%lon(i-1,j)
      lat2=grd%lat(i-1,j)
    else !SW of cell
      lon2=grd%lon(i-1,j-1)
      lat2=grd%lat(i-1,j-1)
    endif
  endif

  L1=berg%length; W1=berg%width; A1=L1*W1
  dlon=lon1-lon2; dlat=lat1-lat2

  ! Note that this is not the exact distance along a great circle.
  ! Approximation for small distances. Should be fine.
  lat_ref=0.5*(lat1+lat2)
  call convert_from_grid_to_meters(lat_ref,bergs%grd%grid_is_latlon,dx_dlon,dy_dlat)

  r_dist_x=dlon*dx_dlon; r_dist_y=dlat*dy_dlat
  r_dist=sqrt( (r_dist_x**2) + (r_dist_y**2) ) !(Stern et al 2017, Eqn 3)

  !Stern et al 2017, Eqn 4: radius of circle inscribed within hexagon or square with area A1 or A2
  if (bergs%hexagonal_icebergs) then
    R1=sqrt(A1/(2.*sqrt(3.)))
  else !square packing
    if (bergs%iceberg_bonds_on) then
      R1=0.5*sqrt(A1)
    else
      R1=sqrt(A1/pi) ! Interaction radius of the iceberg (assuming circular icebergs)
    endif
  endif

  !Calculating spring force (Stern et al 2017, Eqn 6):
  spring_coef=bergs%contact_spring_coef
  crit_dist=max(2*R1,bergs%contact_distance)

  radial_damping_coef=bergs%radial_damping_coef
  tangental_damping_coef=bergs%tangental_damping_coef
  critical_interaction_damping_on=bergs%critical_interaction_damping_on

  ! Using critical values for damping rather than manually setting the damping.
  if (critical_interaction_damping_on) then
    radial_damping_coef=2.*sqrt(spring_coef) ! Critical damping
    if (bergs%tang_crit_int_damp_on) then
      tangental_damping_coef=(2.*sqrt(spring_coef))/4 ! Critical damping (just a guess)
    endif
  endif

  if  (r_dist<crit_dist) then
    !Spring force (Stern et al 2017, Eqn 7):
    accel_spring=spring_coef*(crit_dist-r_dist)
    IA_x=IA_x+(accel_spring*(r_dist_x/r_dist))
    IA_y=IA_y+(accel_spring*(r_dist_y/r_dist))

    ! Damping force (Stern et al 2017, Eqn 8):
    ! Parallel velocity
    !projection matrix in Stern et al 2017, Eqn 8:
    P_11=(r_dist_x*r_dist_x)/(r_dist**2)
    P_12=(r_dist_x*r_dist_y)/(r_dist**2)
    P_21=(r_dist_x*r_dist_y)/(r_dist**2)
    P_22=(r_dist_y*r_dist_y)/(r_dist**2)
    p_ia_coef=radial_damping_coef
    P_ia_11=P_ia_11+p_ia_coef*P_11
    P_ia_12=P_ia_12+p_ia_coef*P_12
    P_ia_21=P_ia_21+p_ia_coef*P_21
    P_ia_22=P_ia_22+p_ia_coef*P_22
    !not needed because u2 and v2 are zero
    !P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
    !P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))

    ! Normal velocities
    P_11=1-P_11  ;  P_12=-P_12 ; P_21= -P_21 ;    P_22=1-P_22
    p_ia_coef=tangental_damping_coef
    P_ia_11=P_ia_11+p_ia_coef*P_11
    P_ia_12=P_ia_12+p_ia_coef*P_12
    P_ia_21=P_ia_21+p_ia_coef*P_21
    P_ia_22=P_ia_22+p_ia_coef*P_22
    !not needed because u2 and v2 are zero
    !P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
    !P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))
  endif
end subroutine calculate_force_land_contact

!> Calculates area of overlap between two circular bergs
subroutine overlap_area(R1, R2, d, A, trapped)
  ! Arguments
  real, intent(in) :: R1 !< Radius of berg 1 (m)
  real, intent(in) :: R2 !< Radius of berg 2 (m)
  real, intent(in) :: d !< Separation of berg centers (m)
  real, intent(out) :: A !< Overlap area (m2)
  real, intent(out) :: Trapped !< =1. if one berg is completely inside the other, =0. otherwise
  ! Local variables
  real :: R1_sq, R2_sq, d_sq

  R1_sq=R1**2
  R2_sq=R2**2
  d_sq=d**2
  Trapped=0.

  if (d>0.) then
    if (d<(R1+R2)) then
      if (d>abs(R1-R2)) then
        A=(R1_sq*acos((d_sq+R1_sq-R2_sq)/(2.*d*R1))) + (R2_sq*acos((d_sq+R2_sq-R1_sq)/(2.*d*R2))) - (0.5*sqrt((-d+R1+R2)*(d+R1-R2)*(d-R1+R2)*(d+R1+R2)))
      else
        A=min(pi*R1_sq,pi*R2_sq)
        Trapped=1.
      endif
    else
      A=0.
    endif
  else
    A=0. ! No area of perfectly overlapping bergs (ie: a berg interacting with itself)
  endif

end subroutine overlap_area

!> For the MTS scheme, calculates the instantaneous acceleration of an iceberg
subroutine accel_mts(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, rx, ry, &
  ax, ay, axn, ayn, bxn, byn, save_bond_energy, Fec_x, Fec_y, Fdc_x, Fdc_y, debug_flag)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< An iceberg
  integer, intent(in) :: i !< i-index of cell berg is in
  integer, intent(in) :: j !< j-index of cell berg is in
  real, intent(in) :: xi !< Non-dimensional x-position within cell of berg
  real, intent(in) :: yj !< Non-dimensional y-position within cell of berg
  real, intent(in) :: lat !< Latitude of berg (degree N)
  real, intent(inout) :: uvel !< Zonal velocity of berg (m/s)
  real, intent(inout) :: vvel !< Meridional velocity of berg (m/s)
  real, intent(inout) :: uvel0 !< Zonal velocity of berg at beginning of time-step (m/s)
  real, intent(inout) :: vvel0 !< Meridional velocity of berg at beginning of time-step (m/s)
  real, intent(in) :: dt !< Time step (s)
  real, intent(in) :: rx !< Random number between -1 and 1 for use in x-component of stochastic tidal parameterization
  real, intent(in) :: ry !< Random number between -1 and 1 for use in y-component of stochastic tidal parameterization
  real, intent(out) :: ax !< Zonal acceleration (m/s2)
  real, intent(out) :: ay !< Meridional acceleration (m/s2)
  real, intent(inout) :: axn !< Explicit estimate of zonal acceleration (m/s2)
  real, intent(inout) :: ayn !< Explicit estimate of meridional acceleration (m/s2)
  real, intent(inout) :: bxn !< Implicit component of zonal acceleration (m/s2)
  real, intent(inout) :: byn !< Implicit component of meridional acceleration (m/s2)
  logical, intent(in) :: save_bond_energy !< Track energy terms for bonds (when monitor_energy==T)
  real, optional, intent(out) :: Fec_x !< If bergs%mts_part==1: Zonal collisional elastic force (N)
  real, optional, intent(out) :: Fec_y !< If bergs%mts_part==1: Meridional collisional elastic force (N)
  real, optional, intent(out) :: Fdc_x !< If bergs%mts_part==1: Zonal collisional damping force (N)
  real, optional, intent(out) :: Fdc_y !< If bergs%mts_part==1: Meridional collisional damping force (N)
  logical, optional :: debug_flag !< If true, print debugging
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: uo, vo, ui, vi, ua, va, uwave, vwave, ssh_x, ssh_y, sst, sss, cn, hi, od
  real :: f_cori, T, D, W, L, M, F
  real :: drag_ocn, drag_atm, drag_ice, wave_rad
  real :: c_ocn, c_atm, c_ice
  real :: ampl, wmod, Cr, Lwavelength, Lcutoff, Ltop
  real, parameter :: accel_lim=1.e-2, Cr0=0.06, vel_lim=15.
  real :: lambda, detA, A11, A12, A21, A22, RHS_x, RHS_y, D_hi
  real :: uveln, vveln, us, vs, speed, loc_dx, new_speed
  real :: u_star, v_star    !Added by Alon
  real :: IA_x, IA_y    !Added by Alon
  real :: P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y    !Added by Alon
  logical :: dumpit
  logical :: interactive_icebergs_on  ! Flag to decide whether to use forces between icebergs.
  logical :: Runge_not_Verlet  ! Flag to specify whether it is Runge-Kutta or Verlet
  logical :: use_new_predictive_corrective !Flad to use Bob's predictive corrective scheme. (default off)
  integer :: itloop
  integer :: stderrunit
  real :: dragfrac, N_bonds, N_max, groundfrac, c_gnd, drag_gnd, scaling
  real :: minfricvel=3.17e-12 !m/s (~0.0001 m/a)
  type(bond), pointer :: current_bond
  real :: IAd_x, IAd_y, sum_bond_Ee, sum_bond_Ed, oneoverkn, M2
  type(iceberg), pointer :: other_berg

  Runge_not_Verlet=bergs%Runge_not_Verlet  ! Loading directly from namelist/default , Alon
  interactive_icebergs_on=bergs%interactive_icebergs_on  ! Loading directly from namelist/default , Alon
  use_new_predictive_corrective=bergs%use_new_predictive_corrective  ! Loading directly from namelist/default , Alon

  u_star=uvel0+(axn*(dt/2.)); v_star=vvel0+(ayn*(dt/2.))  !Alon
  stderrunit = stderr()   ! Get the stderr unit number.
  grd=>bergs%grd   ! For convenience

  if (bergs%mts_part == 1) then
    !in this case, u_star and v_star equal u_k and v_k, respectively, the
    !the velocity components from the previous cycle final MTS 'fast' force iteration, k.
    u_star=berg%uvel; v_star=berg%vvel
    !might as well set it for the berg velocity, too.
    uvel0=berg%uvel;  vvel0=berg%vvel
    uvel=berg%uvel;   vvel=berg%vvel
  endif

  if (new_mts) then
    scaling=0.5
  else
    scaling=1.0
  endif

  ! Initializing accelerations
  axn=0.; ayn=0.; bxn=0.; byn=0.

  if (.not. bergs%only_interactive_forces) then
    !for mts, gridded fields already saved on berg
    uo=berg%uo; vo=berg%vo; ua=berg%ua; va=berg%va; ui=berg%ui; vi=berg%vi;
    ssh_x=berg%ssh_x; ssh_y=berg%ssh_y; sst=berg%sst; sss=berg%sss;  cn=berg%cn; hi=berg%hi; od=berg%od

    if ((grd%grid_is_latlon) .and. (.not. bergs%use_f_plane)) then
      f_cori=(2.*omega)*sin(pi_180*lat)
    else
      f_cori=(2.*omega)*sin(pi_180*bergs%lat_ref)
    endif
    !f_cori=0.

    M=berg%mass
    T=berg%thickness ! total thickness
    D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
    F=T-D ! freeboard
    W=berg%width; L=berg%length

    hi=min(hi,D)
    D_hi=max(0.,D-hi)


    !start to feel grounded with draught of bergs%h_to_init_grounding meters of the sea floor topography
    !at groundfrac=0, apply no grounding force. at groundfrac=1, apply max grounding force
    if (bergs%h_to_init_grounding>0.0) then
      groundfrac=1.0-(od-D)/bergs%h_to_init_grounding
      groundfrac=max(groundfrac,0.0); groundfrac=min(groundfrac,1.0)
    else
      if (D>od) then
        groundfrac=1.0
      else
        groundfrac=0.0
      endif
    endif
    if (groundfrac>0.0) then
      c_gnd=(bergs%cdrag_grounding*W*L*groundfrac)/M
    else
      c_gnd=0.0
    endif

    ! Wave radiation (Stern et al 2017, Eqs A4-A5)
    uwave=ua-uo; vwave=va-vo  ! Use wind speed rel. to ocean for wave model (aja)
    wmod=uwave*uwave+vwave*vwave ! The wave amplitude and length depend on the wind speed relative to the ocean current
    ! actually wmod is wmod**2 here.
    ampl=0.5*0.02025*wmod ! This is "a", the wave amplitude
    Lwavelength=0.32*wmod ! Surface wave length fitted to data in table at
    ! http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
    Lcutoff=0.125*Lwavelength
    Ltop=0.25*Lwavelength
    Cr=Cr0*min(max(0.,(L-Lcutoff)/((Ltop-Lcutoff)+1.e-30)),1.) ! Wave radiation coefficient fitted to
    ! graph from Carrieres et al.,  POAC Drift Model.
    wave_rad=0.5*rho_seawater/M*Cr*gravity*ampl*min(ampl,F)*(2.*W*L)/(W+L)
    wmod = sqrt(ua*ua+va*va) ! Wind speed
    if (wmod.ne.0.) then
      uwave=ua/wmod ! Wave radiation force acts in wind direction ...
      vwave=va/wmod
    else
      uwave=0.; vwave=0.; wave_rad=0. ! ... and only when wind is present.
    endif

    dragfrac = 1.0
    if ((bergs%iceberg_bonds_on) .and. (bergs%internal_bergs_for_drag)) then
      N_bonds=0.
      N_max=4.0  !Maximum number of bonds that element can form based on shape
      if (bergs%hexagonal_icebergs) N_max=6.0
      ! Determining number of bonds
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        N_bonds=N_bonds+1.0
        current_bond=>current_bond%next_bond
      enddo
      dragfrac = ((N_max-N_bonds)/N_max)
    endif

    ! Weighted drag coefficients (Stern et al 2017, Eqs A1-A3)
    c_ocn=rho_seawater/M*(0.5*Cd_wv*dragfrac*W*(D_hi)+Cd_wh*W*L)
    c_atm=rho_air     /M*(0.5*Cd_av*dragfrac*W*F     +Cd_ah*W*L)
    if (abs(hi).eq.0.) then
      c_ice=0.
    else
      c_ice=rho_ice   /M*(0.5*Cd_iv*dragfrac*W*hi              )
    endif
    if (abs(ui)+abs(vi).eq.0.) c_ice=0.

    !Turning drag off for testing - Alon
    !c_ocn=0.; c_atm=0.; c_ice=0.
    ! Stern et al 2017, explicit accel due to sea surface slope and wave radiation force
    ! (forces F_R and F_SS, respectively, in Eq. 1. Also see Eqn (A4+A6)) :
    ! Half half accelerations  - axn, ayn
    axn=-gravity*ssh_x +wave_rad*uwave; ayn=-gravity*ssh_y +wave_rad*vwave

    ! Interactive spring acceleration - (Does the spring part need to be called twice?)
    if (interactive_icebergs_on) then
      call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0, &
        P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
      axn=axn + IA_x; ayn=ayn + IA_y
    endif
    axn=axn+f_cori*v_star; ayn=ayn-f_cori*u_star

  else
    !only_interactive_forces
    if (interactive_icebergs_on) then
      call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0, &
        P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces
       !axn=axn + IA_x; ayn=ayn + IA_y
    endif
  endif !if (.not. bergs%only_interactive_forces)

  uveln=uvel0;           vveln=vvel0
  us=uvel0   ;           vs=vvel0
  do itloop=1,2 ! Iterate on drag coefficients
    if (itloop .eq. 2) then
      us=uveln ; vs=vveln
    endif

    if (bergs%only_interactive_forces) then
      if (interactive_icebergs_on) then
        if (itloop>1) then
          call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, us,vs, &
            P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
        endif

        ! Solve for implicit accelerations
        RHS_x=(IA_x/2) -scaling*(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x)
        RHS_y=(IA_y/2) -scaling*(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y)

        A11=1+(scaling*dt*P_ia_11); A22=1+(scaling*dt*P_ia_22)
        A12=(scaling*dt*P_ia_12);   A21=(scaling*dt*P_ia_21)
      endif

    else
      !Stern et al 2017, Eqn A1-A3 !Alon's proposed change - using Bob's improved scheme.
      drag_ocn=c_ocn*0.5*(sqrt( (uveln-uo)**2+(vveln-vo)**2 )+sqrt( (uvel0-uo)**2+(vvel0-vo)**2 ))
      drag_atm=c_atm*0.5*(sqrt( (uveln-ua)**2+(vveln-va)**2 )+sqrt( (uvel0-ua)**2+(vvel0-va)**2 ))
      drag_ice=c_ice*0.5*(sqrt( (uveln-ui)**2+(vveln-vi)**2 )+sqrt( (uvel0-ui)**2+(vvel0-vi)**2 ))
      drag_gnd=c_gnd*max(0.5*(sqrt(uveln**2+vveln**2)+sqrt(uvel0**2+vvel0**2)),minfricvel)**(1.0/3.0 - 1.0)
      !RHS ~ similar to accel terms in Stern et al 2017, Eqn B5
      RHS_x=(axn/2) + scaling*(-drag_ocn*(u_star-uo) -drag_atm*(u_star-ua) -drag_ice*(u_star-ui) -drag_gnd*u_star)
      RHS_y=(ayn/2) + scaling*(-drag_ocn*(v_star-vo) -drag_atm*(v_star-va) -drag_ice*(v_star-vi) -drag_gnd*v_star)

      if (interactive_icebergs_on) then
        if (itloop>1) then
          call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, us,vs, &
            P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
        endif
        RHS_x=RHS_x -scaling*(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x)
        RHS_y=RHS_y -scaling*(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y)
      endif

      ! Solve for implicit accelerations
      lambda=drag_ocn+drag_atm+drag_ice+drag_gnd
      !Matrix A = [A11 A12; A21 A22] is the denominator on the rhs of Stern et al 2017, Eqn B7
      A11=1.+scaling*dt*lambda; A22=1.+scaling*dt*lambda
      A12=-scaling*dt*f_cori; A21=scaling*dt*f_cori
      !A12=dt*f_cori  !Removed by ALon (in order to have the entire matrix. I hope the sign is correct)

      ! For Crank-Nicolson Coriolis term.
      A12=A12/2.; A21=A21/2.

      if (interactive_icebergs_on) then
        A11=A11+(scaling*dt*P_ia_11); A22=A22+(scaling*dt*P_ia_22)
        A12=A12+(scaling*dt*P_ia_12); A21=A21+(scaling*dt*P_ia_21)
      endif
    endif !if (bergs%only_interactive_forces)

    !inverse of matrix A * RHS
    detA=1./((A11*A22)-(A12*A21))
    ax=detA*(A22*RHS_x-A12*RHS_y); ay=detA*(A11*RHS_y-A21*RHS_x)

    ! Stern et al 2017, Eqn B5 and B6. Here, ax accounts for both acceleration terms (separated below)
    uveln=u_star+dt*ax; vveln=v_star+dt*ay        ! Alon
  enddo ! itloop

  !(axn & ayn) / (bxn & byn) are the new (explicit) / (implicit) acceleration terms in Stern et al 2017, Eqn B5
  !Saving the totally explicit part of the acceleration to use in finding the next position and u_star -Alon
  if (bergs%only_interactive_forces) then
    axn=IA_x; ayn=IA_y
  else
    axn=-gravity*ssh_x +wave_rad*uwave; ayn=-gravity*ssh_y +wave_rad*vwave
    if (interactive_icebergs_on) then
      axn=axn + IA_x; ayn=ayn + IA_y
    endif
    axn=axn+f_cori*vveln; ayn=ayn-f_cori*uveln    !for Crank Nicolson Coriolis
  endif

  if (new_mts) then
    !ax = 0.5*(axn + bxn)
    !where uvel_new = uvel_old + dt*ax
    bxn = 2*ax-axn; byn = 2*ay-ayn
  else
    !ax = axn/2 + bxn
    !where uvel_new = uvel_old + dt*ax
    bxn= ax-(axn/2); byn= ay-(ayn/2) !Alon
  endif

  if (bergs%mts_part==1) then
    !collisional elastic forces
    if (present(Fec_x) .and. present(Fec_y)) then
      Fec_x=berg%mass*IA_x
      Fec_y=berg%mass*IA_y
    endif
    !collisional damping forces
    if (present(Fdc_x) .and. present(Fdc_y)) then
      Fdc_x=berg%mass*( P_ia_times_u_x - (P_ia_11*uveln+P_ia_12*vveln))
      Fdc_y=berg%mass*( P_ia_times_u_y - (P_ia_21*uveln+P_ia_22*vveln))
    endif
  elseif (save_bond_energy) then
    !save bond energies
    if (monitor_energy .or. bergs%fracture_criterion=='energy') then
      sum_bond_Ee=0.0
      sum_bond_Ed=0.0

      if (bergs%constant_interaction_LW) then
        ! mass according to constant length and width here, just for interactions
        M=bergs%constant_length*bergs%constant_width*berg%thickness*bergs%rho_bergs
      else
        M=berg%mass
      endif

      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        other_berg=>current_bond%other_berg

        if (monitor_energy) then
          !first, the half-step:
          current_bond%Ee=current_bond%Ee - 0.25*dt*(&
            (uvel0 + u_star)*M*current_bond%axn_fast +&
            (vvel0 + v_star)*M*current_bond%ayn_fast)

          if (new_mts) then
            current_bond%Ed=current_bond%Ed - 0.25*dt*(&
              (uvel0 + u_star)*M*current_bond%bxn_fast +&
              (vvel0 + v_star)*M*current_bond%byn_fast)
          endif
        endif

        !the rest:
        IA_x=0.0; IA_y=0.0 !explicit elastic force
        P_ia_11=0.0; P_ia_12=0.0; P_ia_21=0.0; P_ia_22=0.0
        call calculate_force(bergs, berg, other_berg, IA_x, IA_y, uveln, vveln, uveln, vveln,  &
          P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded=.true.)
        !implicit damping force
        IAd_x=P_ia_11*(other_berg%uvel_old-uveln)+P_ia_12*(other_berg%vvel_old-vveln)
        IAd_y=P_ia_12*(other_berg%uvel_old-uveln)+P_ia_22*(other_berg%vvel_old-vveln)

        if (bergs%fracture_criterion=='energy') then
          oneoverkn = 1./bergs%spring_coef
          if (bergs%constant_interaction_LW) then
            ! other_berg mass according to constant length and width here, just for interactions
            M2=bergs%constant_length*bergs%constant_width*other_berg%thickness*bergs%rho_bergs
          else
            M2=other_berg%mass
          endif

          oneoverkn=oneoverkn*(M/min(M,other_berg%mass))
          current_bond%spring_pe = 0.5*((M*IA_x)**2+(M*IA_y)**2)*oneoverkn
        endif

        if (monitor_energy) then
          !work
          current_bond%Ee=current_bond%Ee - 0.25*dt*(&
            (u_star + uveln)*M*IA_x + &
            (v_star + vveln)*M*IA_y)


          if (new_mts) then
            !V{i+0.5} to V{i+1}
            current_bond%Ed=current_bond%Ed- 0.25*dt*(&
              (u_star + uveln)*M*IAd_x + &
              (v_star + vveln)*M*IAd_y)
          else
            current_bond%Ed=current_bond%Ed - 0.5*dt*(&
              (u_star + uveln)*M*IAd_x + &
              (v_star + vveln)*M*IAd_y)
          endif

          sum_bond_Ee = sum_bond_Ee + current_bond%Ee
          sum_bond_Ed = sum_bond_Ed + current_bond%Ed

          current_bond%axn_fast=IA_x;  current_bond%ayn_fast=IA_y
          current_bond%bxn_fast=IAd_x; current_bond%byn_fast=IAd_y
        endif

        current_bond=>current_bond%next_bond
      enddo
    endif
  endif

  ! Limit speed of bergs based on a CFL criteria
  if ((bergs%speed_limit>0.) .or. (bergs%speed_limit .eq.-1.)) then
    speed=sqrt(uveln*uveln+vveln*vveln) ! Speed of berg
    if (speed>0.) then
      loc_dx=min(0.5*(grd%dx(i,j)+grd%dx(i,j-1)),0.5*(grd%dy(i,j)+grd%dy(i-1,j))) ! min(dx,dy)
      !new_speed=min(loc_dx/dt*bergs%speed_limit,speed) ! Restrict speed to dx/dt x factor
      new_speed=loc_dx/dt*bergs%speed_limit ! Speed limit as a factor of dx / dt
      if (new_speed<speed) then
        if (bergs%speed_limit>0.) then
          uveln=uveln*(new_speed/speed) ! Scale velocity to reduce speed
          vveln=vveln*(new_speed/speed) ! without changing the direction
          bergs%nspeeding_tickets=bergs%nspeeding_tickets+1
        else
          call error_mesg('KID, Speeding icebergs', 'Faster than the CFL!', WARNING)
          write(stderrunit,*) 'KID, Speeding berg1! =',mpp_pe(), berg%id
          write(stderrunit,*) 'KID, Speeding berg2, speed =',speed, loc_dx/dt
          write(stderrunit,*) 'KID, Speeding berg3, lat, lon =',lat,xi,yj
        endif
      endif
    endif
  endif

  dumpit=.false.
  if (abs(uveln)>vel_lim.or.abs(vveln)>vel_lim) then
    if (debug) then
      dumpit=.true.
      write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive velocity'
    else
      !write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Excessive velocity detected'
    endif
  endif
  if (abs(ax)>accel_lim.or.abs(ay)>accel_lim) then
    if (debug) then
      dumpit=.true.
      write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive acceleration'
    else
      !write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Excessive acceleration detected'
    endif
  endif
  if (present(debug_flag)) then
    if (debug_flag) dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Debug dump flagged by arguments'
  endif
  if (dumpit .and. (.not. bergs%only_interactive_forces)) then
100 format('pe=',i3,a15,9(x,a8,es12.3))
200 format('pe=',i3,a5,i14,(a8,i5),(a7,es12.3))
    write(stderrunit,200) mpp_pe(),'id0=',berg%id, &
      'yr0=',berg%start_year, 'day0=',berg%start_day, &
      'lon0=',berg%start_lon, 'lat0=',berg%start_lat, 'mass0=',berg%start_mass, &
      'sclng=',berg%mass_scaling
    write(stderrunit,100) mpp_pe(),'Geometry:', &
      'M=',M, 'T=',T, 'D=',D, 'F=',F, 'W=',W, 'L=',L
    write(stderrunit,100) mpp_pe(),'delta U:', &
      'u(n)=',uvel0, 'u(*)=', uvel, 'u(n+1)=',uvel+dt*ax, 'del u=',dt*ax
    write(stderrunit,100) mpp_pe(),'U terms', &
      'f*v=',f_cori*vvel, &
      'g*H_x=',-gravity*ssh_x, &
      'wave*ua=',wave_rad*uwave, &
      'd*(u-uo)=',-drag_ocn*(uvel-uo), &
      'd*(u-ua)=',-drag_atm*(uvel-ua), &
      'd*(u-ui)=',-drag_ice*(uvel-ui)
    write(stderrunit,100) mpp_pe(),'U accel.', &
      'RHS_x=',RHS_x, &
      'ax=',ax, &
      'ax(cori)=',detA*(A11*(f_cori*vvel)+A12*(-f_cori*uvel)), &
      'ax(grav)=',detA*(A11*(-gravity*ssh_x)+A12*(-gravity*ssh_y)), &
      'ax(wave)=',detA*(A11*(wave_rad*uwave)+A12*(wave_rad*vwave)), &
      'ax(ocn)=',detA*(A11*(-drag_ocn*(uvel-uo))+A12*(-drag_ocn*(vvel-vo))), &
      'ax(atm)=',detA*(A11*(-drag_atm*(uvel-ua))+A12*(-drag_atm*(vvel-va))), &
      'ax(ice)=',detA*(A11*(-drag_ice*(uvel-ui))+A12*(-drag_ice*(vvel-vi)))
    write(stderrunit,100) mpp_pe(),'delta V:', &
      'v(n)=',vvel0, 'v(*)=', vvel, 'v(n+1)=',vvel+dt*ay, 'del v=',dt*ay
    write(stderrunit,100) mpp_pe(),'V terms', &
      'f*u=',-f_cori*uvel, &
      'g*H_y=',-gravity*ssh_y, &
      'wave*va=',wave_rad*vwave, &
      'd*(v-vo)=',-drag_ocn*(vvel-vo), &
      'd*(v-va)=',-drag_atm*(vvel-va), &
      'd*(v-vi)=',-drag_ice*(vvel-vi)
    write(stderrunit,100) mpp_pe(),'V accel. pe=', &
      'RHS_y=',RHS_y, &
      'ay=',ay, &
      'ay(cori)=',detA*(-A12*(f_cori*vvel)+A11*(-f_cori*uvel)), &
      'ay(grav)=',detA*(-A12*(-gravity*ssh_x)+A11*(-gravity*ssh_y)), &
      'ay(wave)=',detA*(-A12*(wave_rad*uwave)+A11*(wave_rad*vwave)), &
      'ay(ocn)=',detA*(-A12*(-drag_ocn*(uvel-uo))+A11*(-drag_ocn*(vvel-vo))), &
      'ay(atm)=',detA*(-A12*(-drag_atm*(uvel-ua))+A11*(-drag_atm*(vvel-va))), &
      'ay(ice)=',detA*(-A12*(-drag_ice*(uvel-ui))+A11*(-drag_ice*(vvel-vi)))
    write(stderrunit,100) mpp_pe(),'Vel scales', &
      '|va-vo|=',sqrt((ua-uo)**2+(va-vo)**2), &
      '|vo-vb|=',sqrt((uvel-uo)**2+(vvel-vo)**2), &
      '|va-vb|=',sqrt((uvel-ua)**2+(vvel-va)**2), &
      '|vi-vb|=',sqrt((uvel-ui)**2+(vvel-vi)**2), &
      '|vb|=',sqrt((uvel)**2+(vvel)**2), &
      '|va|=',sqrt((ua)**2+(va)**2), &
      '|vo|=',sqrt((uo)**2+(vo)**2), &
      '|vi|=',sqrt((ui)**2+(vi)**2)
    write(stderrunit,100) mpp_pe(),'Time scales', &
      'f=',f_cori, 'wave_rad=',wave_rad, 'do=',drag_ocn, 'da=',drag_atm, 'di=',drag_ice
    write(stderrunit,100) mpp_pe(),'u*', &
      'd*=',lambda, &
      'u*=',(drag_ocn*uo+drag_atm*ua+drag_ice*ui)/lambda, &
      'uo*=',(drag_ocn*uo)/lambda, &
      'ua*=',(drag_atm*ua)/lambda, &
      'ui*=',(drag_ice*ui)/lambda
    write(stderrunit,100) mpp_pe(),'v*', &
      'd*=',lambda, &
      'v*=',(drag_ocn*vo+drag_atm*va+drag_ice*vi)/lambda, &
      'vo*=',(drag_ocn*vo)/lambda, &
      'va*=',(drag_atm*va)/lambda, &
      'vi*=',(drag_ice*vi)/lambda
    write(stderrunit,100) mpp_pe(),'params', &
      'a=',ampl, 'Lwl=',Lwavelength, 'Lcut=',Lcutoff, 'Ltop=',Ltop, 'hi=',hi, 'Cr=',Cr
    write(stderrunit,100) mpp_pe(),'Position', &
      'xi=',xi, 'yj=',yj, 'lat=',lat
    call dump_locfld(grd,i,j,grd%msk,'MSK')
    call dump_locfld(grd,i,j,grd%ssh,'SSH')
    call dump_locfld(grd,i,j,grd%sst,'SST')
    call dump_locfld(grd,i,j,grd%sss,'SSS')
    call dump_locvel(grd,i,j,grd%uo,'Uo')
    call dump_locvel(grd,i,j,grd%vo,'Vo')
    call dump_locvel(grd,i,j,grd%ua,'Ua')
    call dump_locvel(grd,i,j,grd%va,'Va')
    call dump_locvel(grd,i,j,grd%ui,'Ui')
    call dump_locvel(grd,i,j,grd%vi,'Vi')
    call dump_locfld(grd,i,j,grd%hi,'HI')
    call dump_locfld(grd,i,j,grd%cn,'CN')
    call dump_locvel(grd,i,j,grd%lon,'Lon')
    call dump_locvel(grd,i,j,grd%lat,'Lat')
    call print_berg(stderrunit,berg,'KID, accel, large accel')
  endif

  !Used for testing the ocean response to fixed iceberg motion.
  if (bergs%override_iceberg_velocities) then
    ax  = 0.0;  ay  = 0.0
    axn = 0.0;  ayn = 0.0
    bxn = 0.0;  byn = 0.0
  endif

end subroutine accel_mts

!> Calculates the instantaneous acceleration of an iceberg due to berg interactions
!> For the inner steps of the mts scheme only. Explicit velocity verlet (no implicit terms)
subroutine accel_explicit_inner_mts(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt,  &
  ax, ay, axn, ayn, save_bond_energy, debug_flag)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< An iceberg
  integer, intent(in) :: i !< i-index of cell berg is in
  integer, intent(in) :: j !< j-index of cell berg is in
  real, intent(in) :: xi !< Non-dimensional x-position within cell of berg
  real, intent(in) :: yj !< Non-dimensional y-position within cell of berg
  real, intent(in) :: lat !< Latitude of berg (degree N)
  real, intent(inout) :: uvel !< Zonal velocity of berg (m/s)
  real, intent(inout) :: vvel !< Meridional velocity of berg (m/s)
  real, intent(inout) :: uvel0 !< Zonal velocity of berg at beginning of time-step (m/s)
  real, intent(inout) :: vvel0 !< Meridional velocity of berg at beginning of time-step (m/s)
  real, intent(in) :: dt !< Time step (s)
  real, intent(out) :: ax !< Zonal acceleration (m/s2)
  real, intent(out) :: ay !< Meridional acceleration (m/s2)
  real, intent(inout) :: axn !< Explicit estimate of zonal acceleration (m/s2)
  real, intent(inout) :: ayn !< Explicit estimate of meridional acceleration (m/s2)
  logical, intent(in) :: save_bond_energy !< Track energy terms for bonds (when monitor_energy==T)
  logical, optional :: debug_flag !< If true, print debugging
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: bxn, byn
  real :: M
  real :: uveln, vveln, speed, loc_dx, new_speed
  real :: u_star, v_star
  real :: IA_x, IA_y
  real :: P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y    !Added by Alon
  integer :: stderrunit
  type(bond), pointer :: current_bond
  real :: IAd_x, IAd_y, sum_bond_Ee, sum_bond_Ed, oneoverkn, M2
  type(iceberg), pointer :: other_berg
  integer :: grdi,grdj
  logical :: bonded
  ! DEM-mode local variables
  real :: F_x, F_y, T, Fd_x, Fd_y, T_d, R1, at

  u_star=uvel0+(axn*(dt/2.)); v_star=vvel0+(ayn*(dt/2.))  !Alon
  stderrunit = stderr()   ! Get the stderr unit number.
  grd=>bergs%grd   ! For convenience

  ! Initializing accelerations
  axn=0.; ayn=0.; bxn=0.; byn=0.
  IA_x=0.; IA_y=0.; IAd_x=0.; IAd_y=0.

  if (bergs%dem) then
    F_x=0.; F_y=0.; T=0.; Fd_x=0.; Fd_y=0.; T_d=0.
  endif

  bonded=.true.
  if (bergs%iceberg_bonds_on) then
    current_bond=>berg%first_bond
    do while (associated(current_bond)) ! loop over all bonds
      other_berg=>current_bond%other_berg
      if (.not. associated(other_berg)) then
        call error_mesg('KID,bond interactions', 'Trying to do Bond interactions with unassosiated berg!' ,FATAL)
      else
        if (bergs%dem) then
          call calculate_force_dem(bergs,berg,other_berg,current_bond,&
            dt,F_x,F_y,T,Fd_x,Fd_y,T_d,savestress=.true.)
          other_berg%id=-other_berg%id !mark, as to not repeat in contact search below
        else
          P_ia_11=0. ; P_ia_12=0. ;  P_ia_21=0.;  P_ia_22=0.
          P_ia_times_u_x=0. ; P_ia_times_u_y=0.
          call calculate_force(bergs, berg, other_berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0,  &
            P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded)
          other_berg%id=-other_berg%id !mark, as to not repeat in contact search below
          IAd_x=IAd_x+P_ia_11*(other_berg%uvel_old-berg%uvel_old)+P_ia_12*(other_berg%vvel_old-berg%vvel_old)
          IAd_y=IAd_y+P_ia_12*(other_berg%uvel_old-berg%uvel_old)+P_ia_22*(other_berg%vvel_old-berg%vvel_old)
        endif
      endif
      current_bond=>current_bond%next_bond
    enddo
    !contact search: find any surrounding non-bonded bergs in the same conglom, and processes
    !them for contact using the regular (non-contact) spring constant and crit_dist based on radii
    bonded=.false.
    do grdj = max(berg%jne-2,bergs%grd%jsd+1),min(berg%jne+2,bergs%grd%jed);&
      do grdi = max(berg%ine-2,bergs%grd%isd+1),min(berg%ine+2,bergs%grd%ied)
      other_berg=>bergs%list(grdi,grdj)%first
      do while (associated(other_berg))
        if (other_berg%id>0 .and. other_berg%conglom_id.eq.berg%conglom_id) then
          P_ia_11=0. ; P_ia_12=0. ;  P_ia_21=0.;  P_ia_22=0.
          P_ia_times_u_x=0. ; P_ia_times_u_y=0.
          call calculate_force(bergs, berg, other_berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0,  &
            P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded,c_crit_dist=.true.)
          IAd_x=IAd_x+P_ia_11*(other_berg%uvel_old-berg%uvel_old)+P_ia_12*(other_berg%vvel_old-berg%vvel_old)
          IAd_y=IAd_y+P_ia_12*(other_berg%uvel_old-berg%uvel_old)+P_ia_22*(other_berg%vvel_old-berg%vvel_old)
        endif
        other_berg=>other_berg%next
      enddo
    enddo;enddo

    !unmark the bonded bergs
    current_bond=>berg%first_bond
    do while (associated(current_bond))
      other_berg=>current_bond%other_berg
      other_berg%id=abs(other_berg%id)
      current_bond=>current_bond%next_bond
    enddo
  endif

  if (bergs%dem) then

    if (bergs%hexagonal_icebergs) then
      R1=sqrt(berg%length*berg%width/(2.*sqrt(3.)))
    else !square packing
      if (bergs%iceberg_bonds_on) then
        R1=0.5*sqrt(berg%length*berg%width)
      else
        R1=sqrt(berg%length*berg%width/pi) ! Interaction radius of the iceberg (assuming circular icebergs)
      endif
    endif

    if (bergs%dem_beam_test>0) then
      if (bergs%dem_beam_test==1) then
        !3.1 (Wang,2020): Simply supported beam test
        !no vertical load on ends of beam, add vertical load to center of beam
        if (berg%start_lon==bergs%dem_tests_start_lon .or. berg%start_lon==bergs%dem_tests_end_lon) then
          F_y=0.0; Fd_y=0.0
        elseif (berg%start_lon==0.5*(bergs%dem_tests_start_lon+bergs%dem_tests_end_lon)) then
          F_y = F_y - 1.5e5
        endif
      elseif (bergs%dem_beam_test==2) then
        !3.2 (Wang,2020): Cantilever beam test
        !add vertical load to end of the cantilever beam
        if (berg%start_lon==bergs%dem_tests_end_lon) then
          F_y = F_y - 1.5e10/3.
        endif
      elseif (bergs%dem_beam_test==3) then
        !Angular velocity test (starting element has constant rot=-pi/2)
        if (berg%start_lon==bergs%dem_tests_start_lon) then
          T=0.0; T_d=0.0
        endif
      elseif (bergs%dem_beam_test==4) then
        if (berg%start_lon==bergs%dem_tests_start_lon) then
          F_y=0.0; Fd_y=0.0
        endif
      endif
    endif

    if (bergs%constant_interaction_LW) then
      ! mass according to constant length and width here, just for interactions
      M=bergs%constant_length*bergs%constant_width*berg%thickness*bergs%rho_bergs
    else
      M=berg%mass
    endif

    IA_x=IA_x+F_x/M; IA_y=IA_y+F_y/M
    IAd_x=IAd_x+Fd_x/M; IAd_y=IAd_y+Fd_y/M

    at=(T+T_d)/(0.5*M*R1**2) !torque acceleration
    berg%ang_accel = at
  endif

  axn=IA_x+IAd_x; ayn=IA_y+IAd_y
  bxn=0.0; byn=0.0

  if (new_mts) then
    ax=0.5*(axn+bxn); ay=0.5*(ayn+byn)
  else
    ax=0.5*axn+bxn; ay=0.5*ayn+byn
  endif

  uveln=u_star+dt*ax; vveln=v_star+dt*ay

  if (save_bond_energy) then

    !save bond energies
    if (monitor_energy .or. bergs%fracture_criterion=='energy') then
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        other_berg=>current_bond%other_berg

        if (monitor_energy) then
          !first, the half-step:
          current_bond%Ee=current_bond%Ee - 0.25*dt*(&
            (uvel0 + u_star)*M*current_bond%axn_fast +&
            (vvel0 + v_star)*M*current_bond%ayn_fast)

          if (new_mts) then
            !here, bxn_fast and byn_fast are explicit damping accel (not implicit!)
            current_bond%Ed=current_bond%Ed - 0.25*dt*(&
              (uvel0 + u_star)*M*current_bond%bxn_fast +&
              (vvel0 + v_star)*M*current_bond%byn_fast)
          endif
        endif

        !the rest:

        if (bergs%dem) then
          F_x=0.; F_y=0.; T=0.; Fd_x=0.; Fd_y=0.; T_d=0.
          call calculate_force_dem(bergs,berg,other_berg,current_bond,dt,&
            F_x,F_y,T,Fd_x,Fd_y,T_d,savestress=.false.)
          if (bergs%dem_beam_test>0) then
            if (bergs%dem_beam_test==1) then
              !3.1 (Wang,2020): Simply supported beam test
              !no vertical load on ends of beam, add vertical load to center of beam
              if (berg%start_lon==bergs%dem_tests_start_lon .or. berg%start_lon==bergs%dem_tests_end_lon) then
                F_y=0.0; Fd_y=0.0
              elseif (berg%start_lon==0.5*(bergs%dem_tests_start_lon+bergs%dem_tests_end_lon)) then
                F_y = F_y - 1.5e5
              endif
            elseif (bergs%dem_beam_test==2) then
              !3.2 (Wang,2020): Cantilever beam test
              !add vertical load to end of the cantilever beam
              if (berg%start_lon==bergs%dem_tests_end_lon) then
                F_y = F_y - 1.5e10/3.
              endif
            elseif (bergs%dem_beam_test==3) then
              !Angular velocity test (starting element has constant rot=-pi/2)
              if (berg%start_lon==bergs%dem_tests_start_lon) then
                T=0.0; T_d=0.0
              endif
            elseif (bergs%dem_beam_test==4) then
              if (berg%start_lon==bergs%dem_tests_start_lon) then
                F_y=0.0; Fd_y=0.0
              endif
            endif
          endif
          IA_x=F_x/M; IA_y=F_y/M; IAd_x=Fd_x/M; IAd_y=Fd_y/M
        else
          IA_x=0.0; IA_y=0.0 !explicit elastic force
          P_ia_11=0.0; P_ia_12=0.0; P_ia_21=0.0; P_ia_22=0.0
          call calculate_force(bergs, berg, other_berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0,  &
            P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded=.true.)
          !damping force
          IAd_x=P_ia_11*(other_berg%uvel_old-berg%uvel_old)+P_ia_12*(other_berg%vvel_old-berg%vvel_old)
          IAd_y=P_ia_12*(other_berg%uvel_old-berg%uvel_old)+P_ia_22*(other_berg%vvel_old-berg%vvel_old)
        endif

        if (bergs%fracture_criterion=='energy') then
          if (bergs%constant_interaction_LW) then
            ! other_berg mass according to constant length and width here, just for interactions
            M2=bergs%constant_length*bergs%constant_width*other_berg%thickness*bergs%rho_bergs
          else
            M2=other_berg%mass
          endif

          oneoverkn = 1./bergs%spring_coef
          oneoverkn=oneoverkn*(M/min(M,M2))
          current_bond%spring_pe = 0.5*((M*IA_x)**2+(M*IA_y)**2)*oneoverkn
        endif

        if (monitor_energy) then
          !work
          current_bond%Ee=current_bond%Ee - 0.25*dt*(&
            (u_star + uveln)*M*IA_x + &
            (v_star + vveln)*M*IA_y)


          if (new_mts) then
            !V{i+0.5} to V{i+1}
            current_bond%Ed=current_bond%Ed- 0.25*dt*(&
              (u_star + uveln)*M*IAd_x + &
              (v_star + vveln)*M*IAd_y)
          else
            current_bond%Ed=current_bond%Ed - 0.5*dt*(&
              (u_star + uveln)*M*IAd_x + &
              (v_star + vveln)*M*IAd_y)
          endif

          current_bond%axn_fast=IA_x;  current_bond%ayn_fast=IA_y

          !here, bxn_fast and byn_fast are explicit damping accel (not implicit!)
          current_bond%bxn_fast=IAd_x; current_bond%byn_fast=IAd_y
        endif

        current_bond=>current_bond%next_bond
      enddo
    endif
  endif

  ! Limit speed of bergs based on a CFL criteria
  if ((bergs%speed_limit>0.) .or. (bergs%speed_limit .eq.-1.)) then
    speed=sqrt(uveln*uveln+vveln*vveln) ! Speed of berg
    if (speed>0.) then
      loc_dx=min(0.5*(grd%dx(i,j)+grd%dx(i,j-1)),0.5*(grd%dy(i,j)+grd%dy(i-1,j))) ! min(dx,dy)
      !new_speed=min(loc_dx/dt*bergs%speed_limit,speed) ! Restrict speed to dx/dt x factor
      new_speed=loc_dx/dt*bergs%speed_limit ! Speed limit as a factor of dx / dt
      if (new_speed<speed) then
        if (bergs%speed_limit>0.) then
          uveln=uveln*(new_speed/speed) ! Scale velocity to reduce speed
          vveln=vveln*(new_speed/speed) ! without changing the direction
          bergs%nspeeding_tickets=bergs%nspeeding_tickets+1
        else
          call error_mesg('KID, Speeding icebergs', 'Faster than the CFL!', WARNING)
          write(stderrunit,*) 'KID, Speeding berg1! =',mpp_pe(), berg%id
          write(stderrunit,*) 'KID, Speeding berg2, speed =',speed, loc_dx/dt
          write(stderrunit,*) 'KID, Speeding berg3, lat, lon =',lat,xi,yj
        endif
      endif
    endif
  endif

  !Used for testing the ocean response to fixed iceberg motion.
  if (bergs%override_iceberg_velocities) then
    ax  = 0.0;  ay  = 0.0
    axn = 0.0;  ayn = 0.0
    bxn = 0.0;  byn = 0.0
  endif

end subroutine accel_explicit_inner_mts

!> Calculates the instantaneous acceleration of an iceberg
subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, rx, ry, ax, ay, axn, ayn, bxn, byn, debug_flag)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< An iceberg
  integer, intent(in) :: i !< i-index of cell berg is in
  integer, intent(in) :: j !< j-index of cell berg is in
  real, intent(in) :: xi !< Non-dimensional x-position within cell of berg
  real, intent(in) :: yj !< Non-dimensional y-position within cell of berg
  real, intent(in) :: lat !< Latitude of berg (degree N)
  real, intent(in) :: uvel !< Zonal velocity of berg (m/s)
  real, intent(in) :: vvel !< Meridional velocity of berg (m/s)
  real, intent(in) :: uvel0 !< Zonal velocity of berg at beginning of time-step (m/s)
  real, intent(in) :: vvel0 !< Meridional velocity of berg at beginning of time-step (m/s)
  real, intent(in) :: dt !< Time step (s)
  real, intent(in) :: rx !< Random number between -1 and 1 for use in x-component of stochastic tidal parameterization
  real, intent(in) :: ry !< Random number between -1 and 1 for use in y-component of stochastic tidal parameterization
  real, intent(out) :: ax !< Zonal acceleration (m/s2)
  real, intent(out) :: ay !< Meridional acceleration (m/s2)
  real, intent(inout) :: axn !< Explicit estimate of zonal acceleration (m/s2)
  real, intent(inout) :: ayn !< Explicit estimate of meridional acceleration (m/s2)
  real, intent(inout) :: bxn !< Implicit component of zonal acceleration (m/s2)
  real, intent(inout) :: byn !< Implicit component of meridional acceleration (m/s2)
  logical, optional :: debug_flag !< If true, print debugging
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: uo, vo, ui, vi, ua, va, uwave, vwave, ssh_x, ssh_y, sst, sss, cn, hi, od
  real :: f_cori, T, D, W, L, M, F
  real :: drag_ocn, drag_atm, drag_ice, wave_rad
  real :: c_ocn, c_atm, c_ice
  real :: ampl, wmod, Cr, Lwavelength, Lcutoff, Ltop
  real, parameter :: accel_lim=1.e-2, Cr0=0.06, vel_lim=15.
  real :: alpha, beta, C_N
  real :: lambda, detA, A11, A12, A21, A22, RHS_x, RHS_y, D_hi
  real :: uveln, vveln, us, vs, speed, loc_dx, new_speed
  real :: u_star, v_star    !Added by Alon
  real :: IA_x, IA_y    !Added by Alon
  real :: P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y    !Added by Alon
  logical :: dumpit
  logical :: interactive_icebergs_on  ! Flag to decide whether to use forces between icebergs.
  logical :: Runge_not_Verlet  ! Flag to specify whether it is Runge-Kutta or Verlet
  logical :: use_new_predictive_corrective !Flad to use Bob's predictive corrective scheme. (default off)
  integer :: itloop
  integer :: stderrunit
  real :: dragfrac, N_bonds, N_max, groundfrac, c_gnd, drag_gnd
  real :: minfricvel=3.17e-12 !m/s (~0.0001 m/a)
  type(bond), pointer :: current_bond

  Runge_not_Verlet=bergs%Runge_not_Verlet  ! Loading directly from namelist/default , Alon
  interactive_icebergs_on=bergs%interactive_icebergs_on  ! Loading directly from namelist/default , Alon
  use_new_predictive_corrective=bergs%use_new_predictive_corrective  ! Loading directly from namelist/default , Alon

  !These values are no longer set as parameters, but rather can be changed as variables.
  alpha=0.0
  beta=1.0
  C_N=0.0

  !Alon: Verlet requires implicit Coriolis and implicit drag.
  !Alon: Also, I think that the implicit Coriolis with RK gives icebergs which do not complete inertial circles.
  if (.not.Runge_not_Verlet) then
    alpha=1.0
    C_N=1.0
    beta=1.0
    use_new_predictive_corrective=.True.
  endif


  u_star=uvel0+(axn*(dt/2.))  !Alon
  v_star=vvel0+(ayn*(dt/2.))  !Alon

  ! Get the stderr unit number.
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  ! Interpolate gridded fields to berg
  !if (bergs%mts) then
    !gridded fields already saved on berg
    uo=berg%uo; vo=berg%vo; ua=berg%ua; va=berg%va; ui=berg%ui; vi=berg%vi;
    ssh_x=berg%ssh_x; ssh_y=berg%ssh_y; sst=berg%sst; sss=berg%sss;  cn=berg%cn; hi=berg%hi; od=berg%od
  !else
  !  call interp_flds(grd, berg%lon, berg%lat, i, j, xi, yj, rx, ry, uo, vo, ui, vi, ua, va, ssh_x, &
  !    ssh_y, sst, sss, cn, hi, od)
  !end if

  if ((grd%grid_is_latlon) .and. (.not. bergs%use_f_plane)) then
     f_cori=(2.*omega)*sin(pi_180*lat)
  else
     f_cori=(2.*omega)*sin(pi_180*bergs%lat_ref)
  endif
  !f_cori=0.

  M=berg%mass
  T=berg%thickness ! total thickness
  D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
  F=T-D ! freeboard
  W=berg%width
  L=berg%length

  ! Initializing accelerations - Alon. (I am not 100% sure this is needed). I'm not sure what is output if variable is not defined in the subroutine.
  axn=0.
  ayn=0.
  bxn=0.
  byn=0.

  hi=min(hi,D)
  D_hi=max(0.,D-hi)

  !start to feel grounded with draught of bergs%h_to_init_grounding meters of the sea floor topography
  !at groundfrac=0, apply no grounding force. at groundfrac=1, apply max grounding force
  if (bergs%h_to_init_grounding>0.0) then
    groundfrac=1.0-(od-D)/bergs%h_to_init_grounding
    groundfrac=max(groundfrac,0.0); groundfrac=min(groundfrac,1.0)
  else
    if (D>od) then
      groundfrac=1.0
    else
      groundfrac=0.0
    endif
  endif
  if (groundfrac>0.0) then
    c_gnd=(bergs%cdrag_grounding*W*L*groundfrac)/M
  else
    c_gnd=0.0
  endif

  ! Wave radiation (Stern et al 2017, Eqs A4-A5)
  uwave=ua-uo; vwave=va-vo  ! Use wind speed rel. to ocean for wave model (aja)
  wmod=uwave*uwave+vwave*vwave ! The wave amplitude and length depend on the wind speed relative to the ocean current
                               ! actually wmod is wmod**2 here.
  ampl=0.5*0.02025*wmod ! This is "a", the wave amplitude
  Lwavelength=0.32*wmod ! Surface wave length fitted to data in table at
                        ! http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
  Lcutoff=0.125*Lwavelength
  Ltop=0.25*Lwavelength
  Cr=Cr0*min(max(0.,(L-Lcutoff)/((Ltop-Lcutoff)+1.e-30)),1.) ! Wave radiation coefficient fitted to
                                                             ! graph from Carrieres et al.,  POAC Drift Model.
  wave_rad=0.5*rho_seawater/M*Cr*gravity*ampl*min(ampl,F)*(2.*W*L)/(W+L)
  wmod = sqrt(ua*ua+va*va) ! Wind speed
  if (wmod.ne.0.) then
    uwave=ua/wmod ! Wave radiation force acts in wind direction ...
    vwave=va/wmod
  else
    uwave=0.; vwave=0.; wave_rad=0. ! ... and only when wind is present.
  endif

  dragfrac = 1.0
  if ((bergs%iceberg_bonds_on) .and. (bergs%internal_bergs_for_drag)) then
    N_bonds=0.
    N_max=4.0  !Maximum number of bonds that element can form based on shape
    if (bergs%hexagonal_icebergs) N_max=6.0
    ! Determining number of bonds
    current_bond=>berg%first_bond
    do while (associated(current_bond)) ! loop over all bonds
      N_bonds=N_bonds+1.0
      current_bond=>current_bond%next_bond
    enddo
    dragfrac = ((N_max-N_bonds)/N_max)
  endif

  ! Weighted drag coefficients (Stern et al 2017, Eqs A1-A3)
  c_ocn=rho_seawater/M*(0.5*Cd_wv*dragfrac*W*(D_hi)+Cd_wh*W*L)
  c_atm=rho_air     /M*(0.5*Cd_av*dragfrac*W*F     +Cd_ah*W*L)
  if (abs(hi).eq.0.) then
    c_ice=0.
  else
    c_ice=rho_ice     /M*(0.5*Cd_iv*dragfrac*W*hi              )
  endif
  if (abs(ui)+abs(vi).eq.0.) c_ice=0.

!Turning drag off for testing - Alon
!c_ocn=0.
!c_atm=0.
!c_ice=0.

  ! Stern et al 2017, explicit accel due to
  ! sea surface slope and wave radiation force
  ! (forces F_R and F_SS, respectively, in Eq. 1.
  ! Also see Eqn (A4+A6)) :
    ! Half half accelerations  - axn, ayn
    if (.not.Runge_not_Verlet) then
      axn=-gravity*ssh_x +wave_rad*uwave
      ayn=-gravity*ssh_y +wave_rad*vwave
    else
      ! Not half half accelerations  - for RK
      bxn=-gravity*ssh_x +wave_rad*uwave
      byn=-gravity*ssh_y +wave_rad*vwave
    endif

  ! Interactive spring acceleration - (Does the spring part need to be called twice?)
  if (interactive_icebergs_on) then
    call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
    if (.not.Runge_not_Verlet) then
      axn=axn + IA_x
      ayn=ayn + IA_y
    else
      bxn=bxn + IA_x
      byn=byn + IA_y
    endif
  endif

  if (alpha>0.) then ! If implicit Coriolis, use u_star rather than RK4 latest  !Alon
    if (C_N>0.) then !  C_N=1 for Crank Nicolson Coriolis, C_N=0 for full implicit Coriolis !Alon
      axn=axn+f_cori*v_star
      ayn=ayn-f_cori*u_star
    else
      bxn=bxn+f_cori*v_star
      byn=byn-f_cori*u_star
    endif
  else
    bxn=bxn+f_cori*vvel
    byn=byn-f_cori*uvel
  endif

  if (use_new_predictive_corrective) then
    uveln=uvel0; vveln=vvel0 ! Discuss this change with Alistair. Alon thinks that it is needed.
  else
    uveln=uvel; vveln=vvel
  endif

  us=uvel0   ; vs=vvel0
  do itloop=1,2 ! Iterate on drag coefficients
    if (itloop .eq. 2) then
      us=uveln ; vs=vveln
    endif
   !Stern et al 2017, Eqn A1-A3
  if (use_new_predictive_corrective) then
    !Alon's proposed change - using Bob's improved scheme.
    drag_ocn=c_ocn*0.5*(sqrt( (uveln-uo)**2+(vveln-vo)**2 )+sqrt( (uvel0-uo)**2+(vvel0-vo)**2 ))
    drag_atm=c_atm*0.5*(sqrt( (uveln-ua)**2+(vveln-va)**2 )+sqrt( (uvel0-ua)**2+(vvel0-va)**2 ))
    drag_ice=c_ice*0.5*(sqrt( (uveln-ui)**2+(vveln-vi)**2 )+sqrt( (uvel0-ui)**2+(vvel0-vi)**2 ))
    drag_gnd=c_gnd*max(0.5*(sqrt(uveln**2+vveln**2)+sqrt(uvel0**2+vvel0**2)),minfricvel)**(1.0/3.0 - 1.0)
  else
    !Original Scheme
    us=0.5*(uveln+uvel); vs=0.5*(vveln+vvel)
    drag_ocn=c_ocn*sqrt( (us-uo)**2+(vs-vo)**2 )
    drag_atm=c_atm*sqrt( (us-ua)**2+(vs-va)**2 )
    drag_ice=c_ice*sqrt( (us-ui)**2+(vs-vi)**2 )
    drag_gnd=c_gnd*max(sqrt(us**2+vs**2),minfricvel)**(1.0/3.0 - 1.0)
  endif

  !RHS ~ similar to accel terms in Stern et al 2017, Eqn B5
  RHS_x=(axn/2) + bxn
  RHS_y=(ayn/2) + byn

  if (beta>0.) then ! If implicit, use u_star, v_star rather than RK4 latest
    RHS_x=RHS_x - drag_ocn*(u_star-uo) -drag_atm*(u_star-ua) -drag_ice*(u_star-ui) -drag_gnd*u_star
    RHS_y=RHS_y - drag_ocn*(v_star-vo) -drag_atm*(v_star-va) -drag_ice*(v_star-vi) -drag_gnd*v_star
  else
    RHS_x=RHS_x - drag_ocn*(uvel-uo) -drag_atm*(uvel-ua) -drag_ice*(uvel-ui) -drag_gnd*uvel
    RHS_y=RHS_y - drag_ocn*(vvel-vo) -drag_atm*(vvel-va) -drag_ice*(vvel-vi) -drag_gnd*vvel
  endif

  if (interactive_icebergs_on) then
    if (itloop>1) then
      call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, us,vs, &
        P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
    endif
    if (beta>0.) then ! If implicit, use u_star, v_star rather than RK4 latest
      RHS_x=RHS_x -(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x)
      RHS_y=RHS_y -(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y)
    else
      RHS_x=RHS_x - (((P_ia_11*uvel)+(P_ia_12*vvel))-P_ia_times_u_x)
      RHS_y=RHS_y - (((P_ia_21*uvel)+(P_ia_22*vvel))-P_ia_times_u_y)
    endif
    !print *,'Before calculation:', berg%id, IA_x, IA_y, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y
  endif

  ! Solve for implicit accelerations
  if (alpha+beta.gt.0.) then
    if (bergs%only_interactive_forces) then
      RHS_x=(IA_x/2) -(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x)
      RHS_y=(IA_y/2) -(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y)
      A11=1+(dt*P_ia_11)
      A12=(dt*P_ia_12)
      A21=(dt*P_ia_21)
      A22=1+(dt*P_ia_22)
    else
      lambda=drag_ocn+drag_atm+drag_ice+drag_gnd
      !Matrix A = [A11 A12; A21 A22] is the denominator on the rhs of Stern et al 2017, Eqn B7
      A11=1.+beta*dt*lambda
      A22=1.+beta*dt*lambda
      A12=-alpha*dt*f_cori
      A21=alpha*dt*f_cori
      !A12=dt*f_cori  !Removed by ALon (in order to have the entire matrix. I hope the sign is correct)

      if (C_N>0.) then ! For Crank-Nicolson Coriolis term.
        A12=A12/2.
        A21=A21/2.
      endif

      if (interactive_icebergs_on) then
        A11=A11+(dt*P_ia_11)
        A12=A12+(dt*P_ia_12)
        A21=A21+(dt*P_ia_21)
        A22=A22+(dt*P_ia_22)
      endif
    endif

    !inverse of matrix A * RHS
    detA=1./((A11*A22)-(A12*A21))
    ax=detA*(A22*RHS_x-A12*RHS_y)
    ay=detA*(A11*RHS_y-A21*RHS_x)

!Alistair's version removed by Alon
!      detA=1./(A11**2+A12**2)
!      ax=detA*(A11*RHS_x+A12*RHS_y)
!      ay=detA*(A11*RHS_y-A12*RHS_x)
    else
      ax=RHS_x; ay=RHS_y
    endif

    ! Stern et al 2017, Eqn B5 and B6. Here, ax accounts for both acceleration terms (separated below)
    uveln=u_star+dt*ax        ! Alon
    vveln=v_star+dt*ay        ! Alon
  enddo ! itloop

  !Below: axn and ayn are the new explicit acceleration terms in Stern et al 2017, Eqn B5
  !       bxn and byn are the new implicit accerelation terms
  !Saving the totally explicit part of the acceleration to use in finding the next position and u_star -Alon
  if (bergs%only_interactive_forces) then
    axn=IA_x
    ayn=IA_y
  else
    axn=0.; ayn=0.
    if (.not.Runge_not_Verlet) then
      axn=-gravity*ssh_x +wave_rad*uwave
      ayn=-gravity*ssh_y +wave_rad*vwave
      if (interactive_icebergs_on) then
        axn=axn + IA_x
        ayn=ayn + IA_y
      endif
    endif
    if (C_N>0.) then !  C_N=1 for Crank Nicolson Coriolis, C_N=0 for full implicit Coriolis !Alon
      axn=axn+f_cori*vveln
      ayn=ayn-f_cori*uveln
    endif
  endif

  bxn= ax-(axn/2); byn= ay-(ayn/2) !Alon

  ! Limit speed of bergs based on a CFL criteria
  if ((bergs%speed_limit>0.) .or. (bergs%speed_limit .eq.-1.)) then
    speed=sqrt(uveln*uveln+vveln*vveln) ! Speed of berg
    if (speed>0.) then
      loc_dx=min(0.5*(grd%dx(i,j)+grd%dx(i,j-1)),0.5*(grd%dy(i,j)+grd%dy(i-1,j))) ! min(dx,dy)
      !new_speed=min(loc_dx/dt*bergs%speed_limit,speed) ! Restrict speed to dx/dt x factor
      new_speed=loc_dx/dt*bergs%speed_limit ! Speed limit as a factor of dx / dt
      if (new_speed<speed) then
        if (bergs%speed_limit>0.) then
          uveln=uveln*(new_speed/speed) ! Scale velocity to reduce speed
          vveln=vveln*(new_speed/speed) ! without changing the direction
          bergs%nspeeding_tickets=bergs%nspeeding_tickets+1
        else
          call error_mesg('KID, Speeding icebergs', 'Faster than the CFL!', WARNING)
          write(stderrunit,*) 'KID, Speeding berg1! =',mpp_pe(), berg%id
          write(stderrunit,*) 'KID, Speeding berg2, speed =',speed, loc_dx/dt
          write(stderrunit,*) 'KID, Speeding berg3, lat, lon =',lat,xi,yj
        endif
      endif
    endif
  endif

  dumpit=.false.
  if (abs(uveln)>vel_lim.or.abs(vveln)>vel_lim) then
    if (debug) then
      dumpit=.true.
      write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive velocity'
    else
      !write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Excessive velocity detected'
    endif
  endif
  if (abs(ax)>accel_lim.or.abs(ay)>accel_lim) then
    if (debug) then
      dumpit=.true.
      write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive acceleration'
    else
      !write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Excessive acceleration detected'
    endif
  endif
  if (present(debug_flag)) then
    if (debug_flag) dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Debug dump flagged by arguments'
  endif
  if (dumpit) then
100 format('pe=',i3,a15,9(x,a8,es12.3))
200 format('pe=',i3,a5,i14,(a8,i5),(a7,es12.3))
    write(stderrunit,200) mpp_pe(),'id0=',berg%id, &
      'yr0=',berg%start_year, 'day0=',berg%start_day, &
      'lon0=',berg%start_lon, 'lat0=',berg%start_lat, 'mass0=',berg%start_mass, &
      'sclng=',berg%mass_scaling
    write(stderrunit,100) mpp_pe(),'Geometry:', &
      'M=',M, 'T=',T, 'D=',D, 'F=',F, 'W=',W, 'L=',L
    write(stderrunit,100) mpp_pe(),'delta U:', &
      'u(n)=',uvel0, 'u(*)=', uvel, 'u(n+1)=',uvel+dt*ax, 'del u=',dt*ax
    write(stderrunit,100) mpp_pe(),'U terms', &
      'f*v=',f_cori*vvel, &
      'g*H_x=',-gravity*ssh_x, &
      'wave*ua=',wave_rad*uwave, &
      'd*(u-uo)=',-drag_ocn*(uvel-uo), &
      'd*(u-ua)=',-drag_atm*(uvel-ua), &
      'd*(u-ui)=',-drag_ice*(uvel-ui)
    write(stderrunit,100) mpp_pe(),'U accel.', &
      'RHS_x=',RHS_x, &
      'ax=',ax, &
      'ax(cori)=',detA*(A11*(f_cori*vvel)+A12*(-f_cori*uvel)), &
      'ax(grav)=',detA*(A11*(-gravity*ssh_x)+A12*(-gravity*ssh_y)), &
      'ax(wave)=',detA*(A11*(wave_rad*uwave)+A12*(wave_rad*vwave)), &
      'ax(ocn)=',detA*(A11*(-drag_ocn*(uvel-uo))+A12*(-drag_ocn*(vvel-vo))), &
      'ax(atm)=',detA*(A11*(-drag_atm*(uvel-ua))+A12*(-drag_atm*(vvel-va))), &
      'ax(ice)=',detA*(A11*(-drag_ice*(uvel-ui))+A12*(-drag_ice*(vvel-vi)))
    write(stderrunit,100) mpp_pe(),'delta V:', &
      'v(n)=',vvel0, 'v(*)=', vvel, 'v(n+1)=',vvel+dt*ay, 'del v=',dt*ay
    write(stderrunit,100) mpp_pe(),'V terms', &
      'f*u=',-f_cori*uvel, &
      'g*H_y=',-gravity*ssh_y, &
      'wave*va=',wave_rad*vwave, &
      'd*(v-vo)=',-drag_ocn*(vvel-vo), &
      'd*(v-va)=',-drag_atm*(vvel-va), &
      'd*(v-vi)=',-drag_ice*(vvel-vi)
    write(stderrunit,100) mpp_pe(),'V accel. pe=', &
      'RHS_y=',RHS_y, &
      'ay=',ay, &
      'ay(cori)=',detA*(-A12*(f_cori*vvel)+A11*(-f_cori*uvel)), &
      'ay(grav)=',detA*(-A12*(-gravity*ssh_x)+A11*(-gravity*ssh_y)), &
      'ay(wave)=',detA*(-A12*(wave_rad*uwave)+A11*(wave_rad*vwave)), &
      'ay(ocn)=',detA*(-A12*(-drag_ocn*(uvel-uo))+A11*(-drag_ocn*(vvel-vo))), &
      'ay(atm)=',detA*(-A12*(-drag_atm*(uvel-ua))+A11*(-drag_atm*(vvel-va))), &
      'ay(ice)=',detA*(-A12*(-drag_ice*(uvel-ui))+A11*(-drag_ice*(vvel-vi)))
    write(stderrunit,100) mpp_pe(),'Vel scales', &
      '|va-vo|=',sqrt((ua-uo)**2+(va-vo)**2), &
      '|vo-vb|=',sqrt((uvel-uo)**2+(vvel-vo)**2), &
      '|va-vb|=',sqrt((uvel-ua)**2+(vvel-va)**2), &
      '|vi-vb|=',sqrt((uvel-ui)**2+(vvel-vi)**2), &
      '|vb|=',sqrt((uvel)**2+(vvel)**2), &
      '|va|=',sqrt((ua)**2+(va)**2), &
      '|vo|=',sqrt((uo)**2+(vo)**2), &
      '|vi|=',sqrt((ui)**2+(vi)**2)
    write(stderrunit,100) mpp_pe(),'Time scales', &
      'f=',f_cori, 'wave_rad=',wave_rad, 'do=',drag_ocn, 'da=',drag_atm, 'di=',drag_ice
    write(stderrunit,100) mpp_pe(),'u*', &
      'd*=',lambda, &
      'u*=',(drag_ocn*uo+drag_atm*ua+drag_ice*ui)/lambda, &
      'uo*=',(drag_ocn*uo)/lambda, &
      'ua*=',(drag_atm*ua)/lambda, &
      'ui*=',(drag_ice*ui)/lambda
    write(stderrunit,100) mpp_pe(),'v*', &
      'd*=',lambda, &
      'v*=',(drag_ocn*vo+drag_atm*va+drag_ice*vi)/lambda, &
      'vo*=',(drag_ocn*vo)/lambda, &
      'va*=',(drag_atm*va)/lambda, &
      'vi*=',(drag_ice*vi)/lambda
    write(stderrunit,100) mpp_pe(),'params', &
      'a=',ampl, 'Lwl=',Lwavelength, 'Lcut=',Lcutoff, 'Ltop=',Ltop, 'hi=',hi, 'Cr=',Cr
    write(stderrunit,100) mpp_pe(),'Position', &
      'xi=',xi, 'yj=',yj, 'lat=',lat
    call dump_locfld(grd,i,j,grd%msk,'MSK')
    call dump_locfld(grd,i,j,grd%ssh,'SSH')
    call dump_locfld(grd,i,j,grd%sst,'SST')
    call dump_locfld(grd,i,j,grd%sss,'SSS')
    call dump_locvel(grd,i,j,grd%uo,'Uo')
    call dump_locvel(grd,i,j,grd%vo,'Vo')
    call dump_locvel(grd,i,j,grd%ua,'Ua')
    call dump_locvel(grd,i,j,grd%va,'Va')
    call dump_locvel(grd,i,j,grd%ui,'Ui')
    call dump_locvel(grd,i,j,grd%vi,'Vi')
    call dump_locfld(grd,i,j,grd%hi,'HI')
    call dump_locfld(grd,i,j,grd%cn,'CN')
    call dump_locvel(grd,i,j,grd%lon,'Lon')
    call dump_locvel(grd,i,j,grd%lat,'Lat')
    call print_berg(stderrunit,berg,'KID, accel, large accel')
  endif

  !Used for testing the ocean response to fixed iceberg motion.
  if (bergs%override_iceberg_velocities) then
    ax  = 0.0;  ay  = 0.0
    axn = 0.0;  ayn = 0.0
    bxn = 0.0;  byn = 0.0
  endif

end subroutine accel

!> Print 3x3 cells from 2d array A
subroutine dump_locfld(grd, i0, j0, A, lbl)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  integer, intent(in) :: i0 !< i-index of center of 3x3 patch to print
  integer, intent(in) :: j0 !< j-index of center of 3x3 patch to print
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: A !< Field to print
  character(len=*) :: lbl !< Label to add to messages
  ! Local variables
  integer :: i, j, ii, jj
  real :: B(-1:1,-1:1), fac
  integer :: stderrunit
  stderrunit = stderr()

  do jj=-1,1
    j=max(grd%jsd,min(grd%jed,jj+j0))
    do ii=-1,1
      i=max(grd%isd,min(grd%ied,ii+i0))
      B(ii,jj)=A(i,j)
      if ((i.ne.ii+i0).or.(j.ne.jj+j0)) B(ii,jj)=-9.999999e-99
    enddo
  enddo
  write(stderrunit,'("pe=",i3,x,a8,3i12)') mpp_pe(),lbl,(i0+ii,ii=-1,1)
  do jj=1,-1,-1
    write(stderrunit,'("pe=",i3,x,i8,3es12.4)') mpp_pe(),j0+jj,(B(ii,jj),ii=-1,1)
  enddo
end subroutine dump_locfld

!> Print 2x2 cells from 2d array A
subroutine dump_locvel(grd, i0, j0, A, lbl)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  integer, intent(in) :: i0 !< i-index of NE-cell of 2x2 patch to print
  integer, intent(in) :: j0 !< j-index of NE-cell of 2x2 patch to print
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: A !< Field to print
  character(len=*) :: lbl !< Label to add to messages
  ! Local variables
  integer :: i, j, ii, jj
  real :: B(-1:0,-1:0), fac
  integer :: stderrunit
  stderrunit = stderr()

  do jj=-1,0
    j=max(grd%jsd,min(grd%jed,jj+j0))
    do ii=-1,0
      i=max(grd%isd,min(grd%ied,ii+i0))
      B(ii,jj)=A(i,j)
      if ((i.ne.ii+i0).or.(j.ne.jj+j0)) B(ii,jj)=-9.999999e-99
    enddo
  enddo
  write(stderrunit,'("pe=",i3,x,a8,3i12)') mpp_pe(),lbl,(i0+ii,ii=-1,0)
  do jj=0,-1,-1
    write(stderrunit,'("pe=",i3,x,i8,3es12.4)') mpp_pe(),j0+jj,(B(ii,jj),ii=-1,0)
  enddo
end subroutine dump_locvel


!> Footloose (FL) mechanism
!> Based on England et al. (2020) Modeling the breakup of tabular icebergs. Sci. Adv.
subroutine footloose_calving(bergs, time)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  ! Local Variables
  type(iceberg), pointer :: this
  type(icebergs_gridded), pointer :: grd
  type(bond), pointer :: current_bond
  real, save :: N_max, youngs, poisson, l_c, lw_c, B_c, exp_nlambda, rn
  real, save :: e1, drho, sigmay, lfootparam
  real :: foot_l,foot_mass, foot_area
  type(randomNumberStream),save :: rns ! Random numbers for stochastic tidal parameterization
  logical, save :: Visited=.false.
  integer :: grdi, grdj
  integer, dimension(8) :: seed
  real :: T, W, L, N_bonds, Ln
  real :: IC, max_k, k, l_w, l_b, pu, c, l_b3, Lmin, Wmin, ds, Wn, dA
  real :: fl_disp_x, fl_disp_y, interp_loc, dM_fl_bits

  ! For convenience
  grd=>bergs%grd

  if (.not. Visited) then
    !define some constants:
    !Maximum number of bonds that element can form based on shape
    if (bergs%iceberg_bonds_on) then
      if (bergs%hexagonal_icebergs) then
        N_max=6.0
      else
        N_max=4.0
      endif
    else
      N_max=0.0
    endif

    !for fl_use_l_scale
    e1=exp(0.25*pi)
    drho=rho_seawater-bergs%rho_bergs
    sigmay=bergs%fl_strength*1000 !strength (Pa)
    !to get foot length for calving, multiply lfootparam by thickness/l_w
    lfootparam=e1*rho_seawater*sigmay/(6*bergs%rho_bergs*gravity*drho)

    poisson=0.3; youngs=bergs%fl_youngs
    l_c  = pi/(2.*sqrt(2.)) !for length-scale of child iceberg
    lw_c = 1./(gravity*rho_seawater) !for buoyancy length
    B_c  = youngs/(12.*(1.-poisson**2.)) !for bending stiffness
    exp_nlambda=exp(-bergs%fl_r) !e^(-r*dt)
    seed = constructSeed(mpp_pe(),mpp_pe(),time) !Seed random numbers for Poisson distribution
    rns = initializeRandomNumberStream(seed)
    call getRandomNumbers(rns, rn)
    Visited=.true.
  endif

  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec !computational domain only
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))

      !only non-static, non-footloose-child bergs are eligible for footloose calving:
      if (.not. (this%static_berg.eq.1 .or. this%fl_k.lt.0)) then

        T=this%thickness; W=this%width; L=this%length

        N_bonds=0
        if (bergs%iceberg_bonds_on) then
          N_bonds=this%n_bonds
          if (N_bonds>0) then
            call error_mesg('KID,footloose_calving', &
              'Bonded footloose calving not yet fully implemented!', FATAL)
          endif
        endif

        l_w  = (lw_c*B_c*(T**3.))**0.25  !buoyancy length
        l_b  = l_c*l_w                   !length of child icebergs
        l_b3 = 3*l_b

        !---Determine k, the number of child bergs to calve via the footloose mechanism---!
        !Bergs%fl_r is the average number of bergs to calve over the time increment
        !given in seconds by bergs%fl_r_s. K is set according to a Poisson distribution
        !or may be set directly equal to bergs%fl_r.
        !England et al., 2020 used fl_r=4, fl_r_s= 86400 seconds (1 day), and daily
        !timesteps, so that k (a whole number) gave the number of bergs to calve each
        !timestep. Here, smaller timesteps are typically used, so k is scaled accordingly.
        !Therefore, k is no longer a whole number, and is instead accumulated over time
        !on this%fl_k. Then, floor(this%fl_k) child bergs are calved each timestep, after
        !which, this number of child bergs is subtracted from this%fl_k

        !A max possible number of child bergs to calve, max_k, is set to prevent footloose
        !calving on very small icebergs or interior elements of bonded conglomerates:
        if (bergs%iceberg_bonds_on .and. N_bonds>0.) then
          if (N_bonds==N_max) then
            max_k=0          !no FL calving allowed for interior bergs
          else
            call max_k_for_edge_elements
          endif
        else
          !non-bonded berg length must be greater than 3*(l_b) for FL calving. Set max_k accordingly.
          if (bergs%fl_k_scale_by_perimeter>0) then
            !the minimum length and width the parent berg can have after footloose
            c = ceiling((L-l_b3)/l_b3); Lmin = L - c*l_b3;
            c = ceiling((W-l_b3)/l_b3); Wmin = W - c*l_b3;

            max_k = max(floor((L*W - Lmin*Wmin)/(l_b3*l_b)),0)
          else
            max_k = max(ceiling((L-l_b3)*W/(3*l_b**2.)),0)
            if (L-max_k*3.*(l_b**2.)/W<=0.) max_k = max_k-1
          endif
        endif

        if (bergs%fl_use_l_scale) then
          if (max_k==0) then
            k=0
          else
            ! !here, fl_k is the accumulated mass from side melt/erosion
            ! foot_l = lfootparam*T/l_w
            ! !this is the mass of side melt/erosion needed to calve a single footloose berg
            ! foot_mass = foot_l * T*(1.-bergs%rho_bergs/rho_seawater) * l_b3 * bergs%rho_bergs
            ! !If bergs%fl_l_scale=1, all side mass loss is above sea level, and goes into creating the
            ! !foot. Adjust to >1 if some mass loss evenly melts the sides at all depths...
            ! foot_mass = bergs%fl_l_scale * foot_mass/T
            ! k=floor(this%fl_k/foot_mass)
            ! if (k>max_k) k=max_k
            ! this%fl_k=this%fl_k-k*foot_mass

            !here, fl_k is the accumulated mass from side melt/erosion
            foot_l = lfootparam*T/l_w
            !this is the area of side melt/erosion needed to calve a single footloose berg
            foot_area = foot_l * l_b3
            foot_area = bergs%fl_l_scale * foot_area
            k=floor(this%fl_k/foot_area)
            if (k>max_k) k=max_k
            this%fl_k=this%fl_k-k*foot_area
          endif
        else
          if (max_k==0) then
            k=0
            this%fl_k=0
          else
            ! If sea ice concentration is <=50%, generate footloose child bergs
            IC=min(1.,grd%cn(grdi,grdj)+bergs%sicn_shift) !sea ice only known on grid when this routine is called
            if (.not. IC>0.5) then

              if (bergs%fl_use_poisson_distribution) then
                k=0; pu=1.0
                do while (pu >= exp_nlambda)
                  call getRandomNumbers(rns, rn)
                  pu=pu*rn
                  k=k+1
                enddo
                k=k-1.0 !FL calving rate (per day or year)
              else
                k=bergs%fl_r
              endif

              if (bergs%iceberg_bonds_on) k=k*((N_max-N_bonds)/N_max) !reduce FL calve rate by number of bonds
              k=k*(bergs%dt/bergs%fl_r_s) !scale k to dt
              if (bergs%fl_k_scale_by_perimeter>0) then
                if (c == 0) then
                  !W not eligible for footloose, so only scale by length
                  k = k * 2.*L/bergs%fl_k_scale_by_perimeter
                else
                  k = k * 2.*(L+W)/bergs%fl_k_scale_by_perimeter
                endif
              endif
              this%fl_k=this%fl_k+k       !update fl_k
            end if

            k=floor(this%fl_k)
          endif

          if (k>max_k) then
            k=max_k
            this%fl_k=0.
          else
            this%fl_k=this%fl_k-k !save the remainder for the next timestep
          endif
        endif

        !if there is footloose calving, create the new bergs and reduce length of parent berg
        if (k>0) then

          !new parent berg length (Ln) and width (Wn), and change in area (dA):
          if (bergs%fl_k_scale_by_perimeter>0 .and. c>0) then
            !for scale by perimeter, both the length and width are reduced by ds
            ds = 0.5*( (L+W) - sqrt( (L+W)**2. - 4. * (l_b3*l_b*k) ) )
            Ln = L-ds; Wn=W-ds
            !Corrections to keep Wn >= Wmin. Ln is adjusted accordingly to retain the same dA.
            if (Wn<Wmin) then
              Ln=Ln*(1-(Wmin-Wn)/Wmin)
              Wn=Wmin
            endif
          else
            ds = k*3.*(l_b**2.)/W !FL mechanism length reduction for parent berg
            Ln=L-ds
            Wn=W
          endif
          dA = L*W - Ln*Wn !change in area. Should equal k*3*l_b**2

          if (bergs%fl_style.eq.'new_bergs') then
            !calve and track footloose bergs
            if (bergs%displace_fl_bergs .and. .not. bergs%fl_use_poisson_distribution) call getRandomNumbers(rns, rn)
            call get_footloose_displacement
            call calve_fl_icebergs(bergs,this,k,l_b,fl_disp_x,fl_disp_y)
            bergs%nbergs_calved_fl=bergs%nbergs_calved_fl+1
          else
            !put FL berg mass in footloose bergs bits.
            dM_fl_bits=bergs%rho_bergs*T*dA
            this%mass_of_fl_bits=this%mass_of_fl_bits+dM_fl_bits
            ! mass flux into fl bits (kg/m2/s)
            if (grd%area(grdi,grdj).ne.0.) then
              grd%fl_bits_src(grdi,grdj)=grd%fl_bits_src(grdi,grdj)+&
                dM_fl_bits/(bergs%dt*grd%area(grdi,grdj))*this%mass_scaling
            endif
          endif

          if (Ln .le. 0 .or. Wn .le. 0) then
            if (N_bonds==0) then
              print *,'l_b,L,W,k,max_k',l_b,L,W,k,max_k
              call error_mesg('KID,footloose_calving', &
                'non-edge element has fully calved from footloose mechanism', FATAL)
            endif
            this%fl_k=-3 !marks this berg to be deleted later
          else
            !update length, width, thickness,and mass of parent berg
            if (bergs%allow_bergs_to_roll .and. N_bonds.eq.0.) call rolling(bergs,T,Wn,Ln)
            this%thickness=T; this%width=Wn; this%length=Ln
            this%mass=this%length*this%width*this%thickness*bergs%rho_bergs
          endif
        endif
      endif

      !Optionally, create a new berg from the FL bits if their mass exceeds a threshold
      if (this%mass_of_fl_bits*this%mass_scaling > bergs%new_berg_from_fl_bits_mass_thres) then
        if (bergs%displace_fl_bergs .and. .not. bergs%fl_use_poisson_distribution) call getRandomNumbers(rns, rn)
        call get_footloose_displacement
        k=floor(this%mass_of_fl_bits*this%mass_scaling/bergs%new_berg_from_fl_bits_mass_thres)
        call calve_fl_icebergs(bergs,this,k,l_b,fl_disp_x,fl_disp_y,berg_from_bits=.true.)
        bergs%nbergs_calved_fl=bergs%nbergs_calved_fl+1
        if (grd%area(grdi,grdj).ne.0.) then
          grd%fl_bits_src(grdi,grdj)=grd%fl_bits_src(grdi,grdj)-&
            k*bergs%new_berg_from_fl_bits_mass_thres/(bergs%dt*grd%area(grdi,grdj))
        endif
      endif

      this=>this%next
    enddo
  enddo; enddo

contains

  !> Get maximum number of footloose bergs that can calve over a timestep for bonded edge elements.
  !> \todo finish this
  subroutine max_k_for_edge_elements
    !for most edge elements, maxk=huge(1.0)
    max_k=huge(1.0)
  end subroutine max_k_for_edge_elements

  !> Determine random displacement of footloose child bergs from parent berg
  subroutine get_footloose_displacement
    real :: lon1,lat1,x1,y1,dxdl1,dydl,xdot2,ydot2
    logical :: on_tangential_plane

    if (.not. bergs%displace_fl_bergs) then
      fl_disp_x=0.0; fl_disp_y=0.0
      return
    endif

    !displace child berg to a random location along one of the sides of the rectangular berg
    if (rn<0.25) then !north side
      interp_loc=4.*rn
      fl_disp_x=this%length*(interp_loc-0.5)
      fl_disp_y=0.5*this%width
    elseif (rn<0.5) then !east side
      interp_loc=4.*(rn-0.25)
      fl_disp_x=0.5*this%length
      fl_disp_y=this%width*(interp_loc-0.5)
    elseif (rn<0.75) then !south side
      interp_loc=4.*(rn-0.5)
      fl_disp_x=this%length*(interp_loc-0.5)
      fl_disp_y=-0.5*this%width
    else !west side
      interp_loc=4.*(rn-0.75)
      fl_disp_x=-0.5*this%length
      fl_disp_y=0.5*this%width*(interp_loc-0.5)
    endif

    if (grd%grid_is_latlon) then
      on_tangential_plane=.false.
      if ((this%lat>89.) .and. (grd%grid_is_latlon)) on_tangential_plane=.true.
      lon1=this%lon; lat1=this%lat
      if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1,this%id)
      call convert_from_meters_to_grid(lat1,grd%grid_is_latlon,dxdl1,dydl)
      if (on_tangential_plane) then
        call rotvec_to_tang(lon1,fl_disp_x,fl_disp_y,xdot2,ydot2)
        x1=x1+xdot2; y1=y1+ydot2
        call rotpos_from_tang(x1,y1,lon1,lat1) !lon1 & lat1 = new FL berg position
      else
        lon1=lon1+fl_disp_x*dxdl1; lat1=lat1+fl_disp_y*dydl !new FL berg position
      endif
      fl_disp_x=lon1-this%lon; fl_disp_y=lat1-this%lat !convert back to displacement from parent berg
    endif

  end subroutine get_footloose_displacement
end subroutine footloose_calving

!> Delete any edge elements that fully calved from the footloose mechanism
subroutine delete_fully_fl_calved_edge_elements(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local Variables
  type(iceberg), pointer :: this, kick_the_bucket
  type(icebergs_gridded), pointer :: grd
  integer :: grdi,grdj

  ! For convenience
  grd=>bergs%grd

  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs
      if (this%fl_k==-3) then
        !berg is marked to delete
        kick_the_bucket=>this
        this=>this%next
        print *,'deleting berg',this%id
        call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
      else
        this=>this%next
      endif
    enddo
  enddo; enddo
end subroutine delete_fully_fl_calved_edge_elements

!> Adjust footloose child bergs to be eligible for contact interactions once they are
!> out of contact range of any other berg for the first time
subroutine adjust_fl_berg_interactivity(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local Variables
  type(iceberg), pointer :: this, other_berg
  type(icebergs_gridded), pointer :: grd
  integer :: grdi, grdj, grdi2, grdj2, nc_x, nc_y
  logical :: contact,radial_contact
  real :: lat1,lat2,lon1,lon2,dlon,dlat,dx_dlon,dy_dlat,lat_ref
  real :: r_dist,crit_dist,R1,R2,rdenom

  ! For convenience
  grd=>bergs%grd

  nc_x=bergs%contact_cells_lon; nc_y=bergs%contact_cells_lat

  if (nc_x.eq.1 .and. nc_y.eq.1) then
    radial_contact=.true. !contact based on berg radii can occur
    if (bergs%hexagonal_icebergs) then
      rdenom=1./(2.*sqrt(3.))
    else
      if (bergs%iceberg_bonds_on) then
        rdenom=1./pi
      else
        rdenom=1./4.
      endif
    endif
  else
    radial_contact=.false.
    crit_dist=bergs%contact_distance**2
  endif

  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%fl_k.eq.-1) then
        !It is a non-interactive, footloose child berg. Test if within contact range of another berg
        contact=.false.
        lat1=this%lat; lon1=this%lon
        if (radial_contact) R1=sqrt(this%length*this%width*rdenom) !radius of current berg
        do grdj2=max(grdj-nc_y,grd%jsd+1),min(grdj+nc_y,grd%jed)
          do grdi2=max(grdi-nc_x,grd%isd+1),min(grdi+nc_x,grd%ied)
            other_berg=>bergs%list(grdi2,grdj2)%first
            do while (associated(other_berg))
              if (other_berg%id.ne.this%id) then
                lat2=other_berg%lat; lon2=other_berg%lon
                dlon=lon2-lon1; dlat=lat2-lat1
                if (radial_contact) then
                  R2=sqrt(other_berg%length*other_berg%width*rdenom)
                  crit_dist=max(R1+R2,bergs%contact_distance)**2
                endif
                if (bergs%grd%grid_is_latlon) then
                  lat_ref=0.5*(lat1+lat2)
                  dx_dlon=pi_180*Rearth*cos(lat_ref*pi_180); dy_dlat=pi_180*Rearth
                  r_dist=(dlon*dx_dlon)**2 + (dlat*dy_dlat)**2
                else
                  r_dist=dlon**2 + dlat**2
                endif
                if (r_dist<crit_dist) then
                  contact=.true.
                  exit
                endif
              endif
              other_berg=>other_berg%next
            enddo
            if (contact) exit
          enddo
          if (contact) exit
        enddo
        !if outside of contact range, mark the berg eligible for future interaction/contact
        if (.not. contact) this%fl_k=-2
      endif
      this=>this%next
    enddo
  enddo; enddo

end subroutine adjust_fl_berg_interactivity

!> Steps forward thermodynamic state of all bergs
subroutine thermodynamics(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(bond), pointer :: current_bond
  real :: M, T, W, L, SST, Vol, Ln, Ln1, Wn, Wn1, Tn, nVol, IC, Dn
  real :: Mv, Me, Mb, melt, dvo, dva, dM, Ss, dMe, dMb, dMv
  real :: Mnew, Mnew1, Mnew2, Hocean
  real :: Mbits, nMbits, dMbitsE, dMbitsM, Lbits, Abits, Mbb
  real :: Ms, N_bonds, N_max !Ice shelf melt, Number of bonds, Max_number of bonds
  integer :: i,j, stderrunit
  type(iceberg), pointer :: this, next
  real, parameter :: perday=1./86400.
  integer :: grdi, grdj, k
  real :: SSS !Temporarily here
  ! Footloose bits stuff
  real :: Lfl, Wfl, Tfl, Mfl, Afl, Me_fl, Mv_fl, Mb_fl, dMfl
  real :: Volfl, Lnfl, Wnfl, Tnfl, nVolfl, Mnew1_fl, Mnew2_fl, Mnew_fl
  real :: dMb_fl, dMv_fl, dMe_fl, dMe_l, dMv_l
  real :: Mbits_fl, dMbitsE_fl, nMbits_fl, Lbits_fl, Abits_fl, Mbb_fl, dMbitsM_fl
  real :: prev_mass_of_fl_bits, prev_mass_of_fl_bergy_bits, M_edit, Mscale_edit,l_b3,fb, kd
  real, parameter :: l_c=pi/(2.*sqrt(2.)),lw_c = 1./(gravity*rho_seawater),B_c=1./(12.*(1.-0.3**2.))

  ! For convenience
  grd=>bergs%grd

  !Initializing
  grd%Uvel_on_ocean(:,:,:)=0.
  grd%Vvel_on_ocean(:,:,:)=0.

  if (bergs%use_mixed_melting .or. bergs%allow_bergs_to_roll) then
    !Maximum number of bonds that element can form based on shape
    if (bergs%hexagonal_icebergs) then
      N_max=6.0
    else
      N_max=4.0
    endif
  endif

  ! Thermodynamics of first halo row is calculated, so that spread mass to ocean works correctly
  do grdj = grd%jsc-1,grd%jec+1 ; do grdi = grd%isc-1,grd%iec+1
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))
      if (debug) call check_position(grd, this, 'thermodynamics (top)')
      ! call interp_flds(grd, this%lon, this%lat, this%ine, this%jne, this%xi, this%yj, 0., 0., &
      !                  this%uo, this%vo, this%ui, this%vi, this%ua, this%va, this%ssh_x, &
      !                  this%ssh_y, this%sst, this%sss,this%cn, this%hi)
      SST=this%sst
      SSS=this%sss
      IC=min(1.,this%cn+bergs%sicn_shift) ! Shift sea-ice concentration
      M=this%mass
      T=this%thickness ! total thickness
      !D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
      !F=T-D ! freeboard
      W=this%width
      L=this%length
      i=this%ine
      j=this%jne
      Vol=T*W*L

      ! Environment
      dvo=sqrt((this%uvel-this%uo)**2+(this%vvel-this%vo)**2)
      dva=sqrt((this%ua-this%uo)**2+(this%va-this%vo)**2)
      Ss=1.5*(dva**0.5)+0.1*dva ! Sea state

      ! Melt rates in m/s
      Mv=max( 7.62e-3*SST+1.29e-3*(SST**2), 0.) &! Buoyant convection at sides
          *perday ! convert to m/s
      Mb=max( 0.58*(dvo**0.8)*(SST+4.0)/(L**0.2), 0.) &! Basal turbulent melting
          *perday ! convert to m/s
      Me=max( 1./12.*(SST+2.)*Ss*(1+cos(pi*(IC**3))) ,0.) &! Wave erosion
          *perday ! convert to m/s

      !same these values Mv and Me separately for footloose bits,
      !in case Mv and Me change later due to number of bonds:
      if (this%mass_of_fl_bits>0.) then
        Mv_fl=Mv; Me_fl=Me !Note: Mb_fl is calculated below
      endif

        if (bergs%use_mixed_melting .or. bergs%allow_bergs_to_roll) then
          N_bonds=0.
          if (bergs%iceberg_bonds_on) then
            ! Determining number of bonds
            current_bond=>this%first_bond
            do while (associated(current_bond)) ! loop over all bonds
              N_bonds=N_bonds+1.0
              current_bond=>current_bond%next_bond
            enddo
          endif
          if  (this%static_berg .eq. 1)  N_bonds=N_max  !Static icebergs melt like ice shelves
        endif

      !For icebergs acting as ice shelves
      if ((bergs%melt_icebergs_as_ice_shelf) .or.(bergs%use_mixed_melting))  then
        if (.not. bergs%use_mixed_layer_salinity_for_thermo)  SSS=35.0
        call find_basal_melt(bergs,dvo,this%lat,SSS,SST,bergs%Use_three_equation_model,T,Ms,this%id)
        Ms=max(Ms,0.) !No refreezing allowed for now
        !Set melt to zero if ocean is too thin.
        if ((bergs%melt_cutoff >=0.) .and. (bergs%apply_thickness_cutoff_to_bergs_melt)) then
          Dn=(bergs%rho_bergs/rho_seawater)*this%thickness ! draught (keel depth)
          if ((grd%ocean_depth(i,j)-Dn) < bergs%melt_cutoff) then
            Ms=0.
          endif
        endif

        if (bergs%use_mixed_melting) then
          Me=((N_max-N_bonds)/N_max)*(Mv+Me)
          Mv=0.0
          Mb=(((N_max-N_bonds)/N_max)*(Mb)) + (N_bonds/N_max)*Ms
        else  !Using Three equation model only.
          Mv=0.0
          Me=0.0
          Mb=Ms
        endif
      endif

      if (bergs%set_melt_rates_to_zero) then
        Mv=0.0
        Mb=0.0
        Me=0.0
      endif

      if (bergs%use_operator_splitting) then
        ! Operator split update of volume/mass
        Tn=max(T-Mb*bergs%dt,0.) ! new total thickness (m)
        nVol=Tn*W*L ! new volume (m^3)
        Mnew1=(nVol/Vol)*M ! new mass (kg)
        dMb=M-Mnew1 ! mass lost to basal melting (>0) (kg)

        Ln1=max(L-Mv*bergs%dt,0.) ! new length (m)
        Wn1=max(W-Mv*bergs%dt,0.) ! new width (m)
        nVol=Tn*Wn1*Ln1 ! new volume (m^3)
        Mnew2=(nVol/Vol)*M ! new mass (kg)
        dMv=Mnew1-Mnew2 ! mass lost to buoyant convection (>0) (kg)

        Ln=max(Ln1-Me*bergs%dt,0.) ! new length (m)
        Wn=max(Wn1-Me*bergs%dt,0.) ! new width (m)
        nVol=Tn*Wn*Ln ! new volume (m^3)
        Mnew=(nVol/Vol)*M ! new mass (kg)
        dMe=Mnew2-Mnew ! mass lost to erosion (>0) (kg)
        dM=M-Mnew ! mass lost to all erosion and melting (>0) (kg)
      else
        ! Update dimensions of berg
        Ln=max(L-(Mv+Me)*(bergs%dt),0.) ! (m)
        Wn=max(W-(Mv+Me)*(bergs%dt),0.) ! (m)
        Tn=max(T-Mb*(bergs%dt),0.) ! (m)
        ! Update volume and mass of berg
        nVol=Tn*Wn*Ln ! (m^3)
        Mnew=(nVol/Vol)*M ! (kg)
        dM=M-Mnew ! (kg)
        dMb=(M/Vol)*(W*L)*Mb*bergs%dt ! approx. mass loss to basal melting (kg)
        dMe=(M/Vol)*(T*(W+L))*Me*bergs%dt ! approx. mass lost to erosion (kg)
        dMv=(M/Vol)*(T*(W+L))*Mv*bergs%dt ! approx. mass loss to buoyant convection (kg)
      endif

      !if footloose is based on length of foot, accumulate side mass loss on fl_k
      !note: this only works when bergs%use_operator_splitting=.true.
      if (bergs%fl_use_l_scale .and. this%fl_k>=0) then
        l_b3 = 3.*l_c*(lw_c*bergs%fl_youngs*B_c*(Tn**3.))**0.25 !child berg length x 3
        if (L>l_b3) then !do not accumulate side made loss for sides < l_b3
          fb = Tn*(1.-bergs%rho_bergs/rho_seawater) !freeboard
          kd = Tn-fb !keel depth
          if (W>l_b3) then
            if (bergs%fl_l_scale_erosion_only) then
              this%fl_k=this%fl_k + (dMe/fb - dMv/kd)/bergs%rho_bergs !horiz area from erosion - from buoy conv
              if (this%fl_k<0) this%fl_k=0
            else
              this%fl_k= this%fl_k + (dMe + dMv)/Tn
            endif
          else
            dMv_l=dMv*(Wn1 + W)/(2.*(Ln1 + W)) !mass loss from length from buoyant convection
            dMe_l=dMe*(Wn+ Wn1)/(2.*(Ln+ Wn1)) !mass loss from length from erosion
            if (bergs%fl_l_scale_erosion_only) then
              this%fl_k=this%fl_k + (dMe_l/fb - dMv_l/kd)/bergs%rho_bergs
              if (this%fl_k<0) this%fl_k=0
            else
              this%fl_k= this%fl_k + (dMe_l + dMv_l)/Tn
            endif
          endif
        endif
      endif

      ! Footloose bits (FL bits).
      if (this%mass_of_fl_bits>0.) then
        call fl_bits_dimensions(bergs,this,Lfl,Wfl,Tfl)
        Mfl=this%mass_of_fl_bits
        Volfl=Lfl*Wfl*Tfl
        !basal turbulent melting
        Mb_fl=max( 0.58*(dvo**0.8)*(SST+4.0)/(Lfl**0.2), 0.) *perday
        Tnfl=max(Tfl-  Mb_fl*bergs%dt,0.) ! new FL thickness (m)
        if (bergs%use_operator_splitting) then
          !basal turbulent melting, continued
          nVolfl=Tnfl*Wfl*Lfl ! new volume (m^3)
          Mnew1_fl=(nVolfl/Volfl)*Mfl ! new mass (kg)
          dMb_fl=Mfl-Mnew1_fl ! mass lost to basal melting (>0) (kg)
          !buoyant convection
          Lnfl=max(Lfl-Mv_fl*bergs%dt,0.) ! new FL length (m)
          Wnfl=max(Wfl-Mv_fl*bergs%dt,0.) ! new FL width (m)
          nVolfl=Tnfl*Wnfl*Lnfl ! new footloose volume (m^3)
          Mnew2_fl=(nVolfl/Volfl)*Mfl !new FL mass (kg), after buoyant convection
          dMv_fl=Mnew1_fl-Mnew2_fl
          !erosion
          Lnfl=max(Lnfl-Me_fl*bergs%dt,0.) ! new FL length (m)
          Wnfl=max(Wnfl-Me_fl*bergs%dt,0.) ! new FL width (m)
          nVolfl=Tnfl*Wnfl*Lnfl ! new footloose volume (m^3)
          Mnew_fl=(nVolfl/Volfl)*Mfl !new FL mass (kg), after all melt and erosion
          dMe_fl=Mnew2_fl-Mnew_fl !FL mass lost to erosion (>0) (kg)
        else
          Lnfl=max(Lfl-(Mv_fl+Me_fl)*bergs%dt,0.) ! (m)
          Wnfl=max(Wfl-(Mv_fl+Me_fl)*bergs%dt,0.) ! (m)
          nVolfl=Tnfl*Wnfl*Lnfl ! new footloose volume (m^3)
          Mnew_fl=(nVolfl/Volfl)*Mfl ! new footloose mass (kg)
          dMb_fl=(Mfl/Volfl)*(Wfl*Lfl)*Mb_fl*bergs%dt       !approx. FL mass loss to basal melting (kg)
          dMe_fl=(Mfl/Volfl)*(Tfl*(Wfl+Lfl))*Me_fl*bergs%dt !approx. FL mass loss to erosion (kg)
          dMv_fl=(Mfl/Volfl)*(Tfl*(Wfl+Lfl))*Mv_fl*bergs%dt !approx. FL mass loss to buoyant convection (kg)
        endif
        dMfl=Mfl-Mnew_fl ! total footloose mass lost to all erosion and melting (>0) (kg)
      else
        dMfl=0.; dMb_fl=0.; dMv_fl=0.; dMe_fl=0.
        Mnew_fl=this%mass_of_fl_bits ! retain previous value incase non-zero
      endif

      ! Bergy bits
      if (bergs%bergy_bit_erosion_fraction>0.) then
        !Parent Berg bergy bits
        Mbits=this%mass_of_bits ! mass of bergy bits (kg)
        dMbitsE=bergs%bergy_bit_erosion_fraction*dMe ! change in mass of bits (kg)
        nMbits=Mbits+dMbitsE ! add new bergy bits to mass (kg)
        Lbits=min(L,W,T,40.) ! assume bergy bits are smallest dimension or 40 meters
        Abits=(Mbits/bergs%rho_bergs)/Lbits ! Effective bottom area (assuming T=Lbits)
        Mbb=max( 0.58*(dvo**0.8)*(SST+2.0)/(Lbits**0.2), 0.) &! Basal turbulent melting (for bits)
             *perday ! convert to m/s
        Mbb=bergs%rho_bergs*Abits*Mbb ! in kg/s
        dMbitsM=min(Mbb*bergs%dt,nMbits) ! bergy bits mass lost to melting (kg)
        nMbits=nMbits-dMbitsM ! remove mass lost to bergy bits melt
        if (Mnew==0.) then ! if parent berg has completely melted then
          dMbitsM=dMbitsM+nMbits ! instantly melt all the bergy bits
          nMbits=0.
        endif

        !Footloose bits bergy bits
        if (this%mass_of_fl_bits>0.) then
          Mbits_fl=this%mass_of_fl_bergy_bits ! mass of footloose bergy bits (kg)
          dMbitsE_fl=bergs%bergy_bit_erosion_fraction*dMe_fl ! change in mass of bits (kg)
          nMbits_fl=Mbits_fl+dMbitsE_fl ! add new bergy bits to mass (kg)
          Lbits_fl=min(Lfl,Wfl,Tfl,40.) ! assume bergy bits are smallest dimension or 40 meters
          Abits_fl=(Mbits_fl/bergs%rho_bergs)/Lbits_fl ! Effective bottom area (assuming T=Lbits)
          Mbb_fl=max( 0.58*(dvo**0.8)*(SST+2.0)/(Lbits_fl**0.2), 0.) &! Basal turbulent melting (for bits)
            *perday ! convert to m/s
          Mbb_fl=bergs%rho_bergs*Abits_fl*Mbb_fl ! in kg/s
          dMbitsM_fl=min(Mbb_fl*bergs%dt,nMbits_fl) ! bergy bits mass lost to melting (kg)
          nMbits_fl=nMbits_fl-dMbitsM_fl ! remove mass lost to bergy bits melt
          if (Mnew_fl==0.) then ! if FL berg associated with these bits has completely melted then
            dMbitsM_fl=dMbitsM_fl+nMbits_fl ! instantly melt all the bergy bits
            nMbits_fl=0.
          endif
        else
          dMbitsE_fl=0.; dMbitsM_fl=0.; Abits_fl=0.; nMbits_fl=0.
        endif
      else
        ! retain previous value incase non-zero
        Abits=0.;    dMbitsE=0.;    dMbitsM=0.;    nMbits=this%mass_of_bits
        Abits_fl=0.; dMbitsE_fl=0.; dMbitsM_fl=0.; nMbits_fl=this%mass_of_fl_bergy_bits
      endif

      ! Add melting to the grid and field diagnostics
      if (grd%area(i,j).ne.0.) then

        melt=(dM-(dMbitsE-dMbitsM)+dMfl-(dMbitsE_fl-dMbitsM_fl))/bergs%dt ! kg/s
        grd%floating_melt(i,j)=grd%floating_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s

        if (grd%id_melt_by_class>0) then
          if (this%lat<0.) then
            k=minloc(abs(bergs%initial_mass_s-this%start_mass),1)
          else
            k=minloc(abs(bergs%initial_mass_n-this%start_mass),1)
          endif
          grd%melt_by_class(i,j,k)=grd%melt_by_class(i,j,k)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
        endif

        melt=melt*this%heat_density ! kg/s x J/kg = J/s
        grd%calving_hflx(i,j)=grd%calving_hflx(i,j)+melt/grd%area(i,j)*this%mass_scaling ! W/m2
        bergs%net_heat_to_ocean=bergs%net_heat_to_ocean+melt*this%mass_scaling*bergs%dt ! J

        melt=dM/bergs%dt ! kg/s
        grd%berg_melt(i,j)=grd%berg_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s

        melt=(dMbitsE+dMbitsE_fl)/bergs%dt ! mass flux into bergy bits in kg/s
        grd%bergy_src(i,j)=grd%bergy_src(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s

        melt=(dMbitsM+dMbitsM_fl)/bergs%dt ! melt rate of bergy bits in kg/s
        grd%bergy_melt(i,j)=grd%bergy_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s

        melt=dMfl/bergs%dt ! melt+erosion rate of fl bits in kg/s
        grd%fl_bits_melt(i,j)=grd%fl_bits_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s

        if (this%fl_k>=0) then
          !Non-footloose berg:
          if(grd%id_fl_parent_melt>0) then
            melt=(dM-(dMbitsE-dMbitsM))/bergs%dt !total melt of the "parent" berg
            grd%fl_parent_melt(i,j)=grd%fl_parent_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling !kg/m2/s
          endif
          if(grd%id_fl_child_melt>0) then
            melt=(dMfl-(dMbitsE_fl-dMbitsM_fl))/bergs%dt !total melt of the "child" FL bits
            grd%fl_child_melt(i,j)=grd%fl_child_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling !kg/m2/s
          endif
          if(grd%id_melt_buoy>0) then
            melt=dMb/bergs%dt ! melt rate due to buoyancy term in kg/s
            grd%melt_buoy(i,j)=grd%melt_buoy(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
          endif
          if(grd%id_melt_eros>0) then
            melt=dMe/bergs%dt ! erosion rate in kg/s
            grd%melt_eros(i,j)=grd%melt_eros(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
          endif
          if(grd%id_melt_conv>0) then
            melt=dMv/bergs%dt ! melt rate due to convection term in kg/s
            grd%melt_conv(i,j)=grd%melt_conv(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
          endif

          if (dMfl>0) then
            if(grd%id_melt_buoy_fl>0) then
              melt=dMb_fl/bergs%dt ! melt rate due to buoyancy term in kg/s
              grd%melt_buoy_fl(i,j)=grd%melt_buoy_fl(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
            endif
            if(grd%id_melt_eros_fl>0) then
              melt=dMe_fl/bergs%dt ! erosion rate in kg/s
              grd%melt_eros_fl(i,j)=grd%melt_eros_fl(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
            endif
            if(grd%id_melt_conv_fl>0) then
              melt=dMv_fl/bergs%dt ! melt rate due to convection term in kg/s
              grd%melt_conv_fl(i,j)=grd%melt_conv_fl(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
            endif
          endif
        else
          !Independently-tracked footloose "child" berg:
          if(grd%id_fl_child_melt>0) then
            melt=(dM-(dMbitsE-dMbitsM))/bergs%dt !total melt of the "child" berg (here, melts like a parent berg)
            grd%fl_child_melt(i,j)=grd%fl_child_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling !kg/m2/s
          endif
          if(grd%id_melt_buoy_fl>0) then
            melt=dMb/bergs%dt ! melt rate due to buoyancy term in kg/s
            grd%melt_buoy_fl(i,j)=grd%melt_buoy_fl(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
          endif
          if(grd%id_melt_eros_fl>0) then
            melt=dMe/bergs%dt ! erosion rate in kg/s
            grd%melt_eros_fl(i,j)=grd%melt_eros_fl(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
          endif
          if(grd%id_melt_conv_fl>0) then
            melt=dMv/bergs%dt ! melt rate due to convection term in kg/s
            grd%melt_conv_fl(i,j)=grd%melt_conv_fl(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
          endif
        endif
      else
        stderrunit = stderr()
        write(stderrunit,*) 'KID, thermodynamics: berg appears to have grounded!!!! PE=',mpp_pe(),i,j
        call print_berg(stderrunit,this,'thermodynamics, grounded')
        if (associated(this%trajectory)) &
          write(stderrunit,*) 'traj=',this%trajectory%lon,this%trajectory%lat
        write(stderrunit,*) 'msk=',grd%msk(i,j),grd%area(i,j)
        call error_mesg('KID, thermodynamics', 'berg appears to have grounded!', FATAL)
      endif

      if (bergs%allow_bergs_to_roll .and. N_bonds.eq.0.) call rolling(bergs,Tn,Wn,Ln)
      !Dn=(bergs%rho_bergs/rho_seawater)*Tn ! re-calculate draught (keel depth) for grounding

      !This option allows iceberg melt fluxes to enter the ocean without the icebergs changing shape
      if (bergs%Iceberg_melt_without_decay) then
        !In this case, the iceberg dimension are reset to their values before
        !the thermodynamics are applied.
        !If the spread_mass is being used to calculate melt, we calculate this
        !before reseting

        !temporarily update the iceberg's footloose and footloose bergy bits masses, which are added
        !to the total mass within subroutine spread_mass_across_ocean_cells:
        prev_mass_of_fl_bits = this%mass_of_fl_bits; prev_mass_of_fl_bergy_bits    = this%mass_of_fl_bergy_bits
        this%mass_of_fl_bits = Mnew_fl;              this%mass_of_fl_bergy_bits = nMbits_fl

        if (bergs%find_melt_using_spread_mass) then
          if (Mnew>0.) then !If the berg still exists
            call spread_mass_across_ocean_cells(bergs,this, i, j, this%xi, this%yj, Mnew, &
              nMbits, this%mass_scaling, Ln*Wn,  Tn, addfootloose=.true.)
          elseif (Mnew_fl>0.) then
            !In the unlikely case that the parent berg has melted completely, but footloose bits still
            !exist, spread mass using footloose in place of the parent berg.
            M_edit=Lnfl*Wnfl*Tnfl*bergs%rho_bergs
            Mscale_edit = Mnew_fl*this%mass_scaling/M_edit
            call spread_mass_across_ocean_cells(bergs,this, i, j, this%xi, this%yj, M_edit, &
              nMbits, Mscale_edit, Lnfl*Wnfl,  Tnfl, addfootloose=.false.)
          endif
        endif
        !Reset all the values
        this%mass_of_fl_bits=prev_mass_of_fl_bits
        this%mass_of_fl_bergy_bits=prev_mass_of_fl_bergy_bits
        Mnew=this%mass
        nMbits=this%mass_of_bits
        Mnew_fl=this%mass_of_fl_bits
        nMbits_fl=this%mass_of_fl_bergy_bits
        Tn=this%thickness
        Wn=this%width
        Ln=this%length
        if (bergs%bergy_bit_erosion_fraction>0.) then
          Mbits=this%mass_of_bits ! mass of bergy bits (kg)
          Lbits=min(L,W,T,40.) ! assume bergy bits are smallest dimension or 40 meters
          Abits=(Mbits/bergs%rho_bergs)/Lbits ! Effective bottom area (assuming T=Lbits)
          if (this%mass_of_fl_bits>0.) then
            Mbits_fl=this%mass_of_fl_bergy_bits ! mass of footloose bergy bits (kg)
            Lbits_fl=min(Lfl,Wfl,Tfl,40.) ! assume bergy bits are smallest dimension or 40 meters
            Abits_fl=(Mbits_fl/bergs%rho_bergs)/Lbits_fl ! Effective bottom area (assuming T=Lbits)
          endif
        endif
      else
        ! Store the new state of iceberg (with L>W)
        this%mass=Mnew
        this%mass_of_bits=nMbits
        this%mass_of_fl_bits=Mnew_fl
        this%mass_of_fl_bergy_bits=nMbits_fl
        this%thickness=Tn
        this%width=min(Wn,Ln)
        this%length=max(Wn,Ln)
      endif
      next=>this%next

      ! Did berg completely melt?
      if (Mnew<=0.) then
        if (Mnew_fl>0) then
          ! The parent berg is melted, but the associated footloose bergs are not.
          ! Replace the parent with the footloose
          bergs%nbergs_calved_fl=bergs%nbergs_calved_fl+1
          this%mass=Lnfl*Wnfl*Tnfl*bergs%rho_bergs
          this%length=Lnfl; this%width=Wnfl; this%thickness=Tnfl
          nMbits_fl=nMbits_fl*this%mass_scaling
          this%mass_scaling = Mnew_fl*this%mass_scaling/this%mass
          this%mass_of_bits=nMbits_fl/this%mass_scaling
          this%mass_of_fl_bits=0.
          this%mass_of_fl_bergy_bits=0.
          this%fl_k=-1.
          this%start_year=bergs%current_year
          this%start_day=bergs%current_yearday
          if (grd%area(i,j).ne.0.) then
            grd%fl_bits_src(i,j)=grd%fl_bits_src(i,j)-&
              this%mass*this%mass_scaling/(bergs%dt*grd%area(i,j))
          endif
        else
          ! Delete the berg
          if (.not. bergs%debug_write) call move_trajectory(bergs, this)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first, this)
        endif
        bergs%nbergs_melted=bergs%nbergs_melted+1
      endif
      this=>next
    enddo
  enddo ; enddo
end subroutine thermodynamics

!> Rolling
!> There are now 3 iceberg rolling schemes:
!> 1) Rolling based on aspect ratio threshold (iceberg of constant density)
!> 2) Rolling based on corrected Weeks and Mellor scheme
!> 3) Rolling based on incorrect Weeks and Mellor scheme - kept for legacy reasons
subroutine rolling(bergs,Tn,Wn,Ln)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  real :: Tn !< Berg thickness
  real :: Wn !< Berg width
  real :: Ln !< Berg length
  ! Local variables
  real, parameter :: Delta=6.0
  real :: Dn, q, tip_parameter

  Dn=(bergs%rho_bergs/rho_seawater)*Tn ! draught (keel depth)
  if ( Dn>0. ) then
    if ( (.not.bergs%use_updated_rolling_scheme) .and. (bergs%tip_parameter<999.) ) then    !Use Rolling Scheme 3
      if ( max(Wn,Ln)<sqrt(0.92*(Dn**2)+58.32*Dn) ) then
        call swap_variables(Tn,Wn)
        if (Wn>Ln) call swap_variables(Wn,Ln)
      endif
    else
      if (Wn>Ln) call swap_variables(Ln,Wn)  !Make sure that Wn is the smaller dimension

      if ( (.not.bergs%use_updated_rolling_scheme) .and. (bergs%tip_parameter>=999.) ) then    !Use Rolling Scheme 2
        q=bergs%rho_bergs/rho_seawater
        if (Wn<sqrt((6.0*q*(1-q)*(Tn**2))-(12*Delta*q*Tn))) then
          call swap_variables(Tn,Wn)
          if (Wn>Ln) call swap_variables(Wn,Ln)
        endif
      endif

      if (bergs%use_updated_rolling_scheme) then    !Use Rolling Scheme 1
        if (bergs%tip_parameter>0.) then
          tip_parameter=bergs%tip_parameter
        else
          ! Equation 27 from Burton et al 2012, or equivalently, Weeks and Mellor 1979 with constant density
          ! Tip_parameter=0.92 if using default values
          tip_parameter=sqrt(6*(bergs%rho_bergs/rho_seawater)*(1-(bergs%rho_bergs/rho_seawater)))
        endif
        if ((tip_parameter*Tn)>Wn)  then     !note that we use the Thickness instead of the Draft
          call swap_variables(Tn,Wn)
          if (Wn>Ln) call swap_variables(Wn,Ln)
        endif
      endif
    endif
    !Dn=(bergs%rho_bergs/rho_seawater)*Tn ! re-calculate draught (keel depth) for grounding
  endif

contains

  !> swap values between two parameters
  subroutine swap_variables(x,y)
    ! Arguments
    real, intent(inout) :: x !< first parameter
    real, intent(inout) :: y !< second parameter
    real :: temp
    temp=x
    x=y
    y=temp
  end subroutine swap_variables
end subroutine rolling

!> Estimate the dimensions of a footloose berg that represents the footloose bits.
!! The goal is that the meltwater flux from this representative berg, multiplied by
!! (the mass of all footloose bits)/(mass of the representative berg), matches the meltwater
!! flux of all footloose bits (which if they were individually tracked, could vary in size).
subroutine fl_bits_dimensions(bergs,this,L_fl,W_fl,T_fl)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: this !< Parent berg to the footloose bits
  real :: L_fl !< Footloose berg bits length
  real :: W_fl !< Footloose berg bits width
  real :: T_fl !< Footloose berg bits thickness
  ! Local variables
  real,parameter :: l_c=pi/(2.*sqrt(2.)), lw_c = 1./(gravity*rho_seawater)
  real,parameter :: B_c=1./(12.*(1.-0.3**2.)) !poisson=0.3
  real :: l_w,l_b

  l_w  = (lw_c*bergs%fl_youngs*B_c*(this%thickness**3.))**0.25  !buoyancy length
  l_b  = l_c*l_w !length of a freshly calved footloose child berg
  L_fl = bergs%fl_bits_scale_l*3.*l_b; W_fl=bergs%fl_bits_scale_w*l_b
  T_fl=bergs%fl_bits_scale_t*this%thickness
  call rolling(bergs,T_fl,W_fl,L_fl)
end subroutine fl_bits_dimensions

!> Generate gridded fields from icebergs
subroutine create_gridded_icebergs_fields(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
! Local variables
type(icebergs_gridded), pointer :: grd
type(iceberg), pointer :: this
integer i,j
integer :: grdi, grdj
real :: Hocean, Dn,Tn,dvo, mass_tmp
real :: ustar_h, ustar
real :: orientation
real :: ave_thickness, ave_draft
real, dimension(bergs%grd%isd:bergs%grd%ied,bergs%grd%jsd:bergs%grd%jed)  :: spread_mass_tmp
real :: tmp

  ! For convenience
  grd=>bergs%grd

  spread_mass_tmp(:,:)=0. !Initializing temporary variable to use in iceberg melt calculation

  !Special case for icebergs not decaying, but mass diffence being used for melt rates
  if ((bergs%find_melt_using_spread_mass) .and.  (bergs%Iceberg_melt_without_decay)) then
    call sum_up_spread_fields(bergs, spread_mass_tmp(grd%isc:grd%iec,grd%jsc:grd%jec),'mass')
  endif

  !Loop through icebergs and spread mass on ocean
  call calculate_mass_on_ocean(bergs, with_diagnostics=.true.)

  !Finding the spread fields
  if ((grd%id_spread_uvel>0)  .or. (bergs%pass_fields_to_ocean_model)) then
    grd%spread_uvel(:,:)=0.
    call sum_up_spread_fields(bergs, grd%spread_uvel(grd%isc:grd%iec,grd%jsc:grd%jec), 'Uvel')
  endif
  if ( (grd%id_spread_vvel>0)  .or. (bergs%pass_fields_to_ocean_model)) then
    grd%spread_vvel(:,:)=0.
    call sum_up_spread_fields(bergs, grd%spread_vvel(grd%isc:grd%iec,grd%jsc:grd%jec), 'Vvel')
  endif
  if ( (grd%id_spread_area>0)  .or. (bergs%pass_fields_to_ocean_model)) then
    grd%spread_area(:,:)=0.
    call sum_up_spread_fields(bergs, grd%spread_area(grd%isc:grd%iec,grd%jsc:grd%jec), 'area')
  endif
  !Always find spread_mass since it is used for so many things.
  grd%spread_mass(:,:)=0.
    call sum_up_spread_fields(bergs, grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec),'mass')

  !Using spread_mass_to_ocean to calculate melt rates (if this option is chosen)
  if (bergs%find_melt_using_spread_mass) then
    if (.not. bergs%Iceberg_melt_without_decay) &
         spread_mass_tmp(grd%isc:grd%iec,grd%jsc:grd%jec)= grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec)
    do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
      if (grd%area(i,j)>0.0) then
        grd%floating_melt(i,j)=max((grd%spread_mass_old(i,j) - spread_mass_tmp(i,j))/(bergs%dt),0.0)
      else
        grd%floating_melt(i,j)=0.0
      endif
    enddo ;enddo
    grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)=grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*HLF !Not 100% sure this is correct.
  endif

  ! Dividing the gridded iceberg momentum diagnostic by the iceberg mass to get velocities
  if ((grd%id_u_iceberg>0) .or. (grd%id_v_iceberg>0)) then
    do j = grd%jsc,grd%jec ; do i = grd%isc,grd%iec
      if (grd%mass(i,j)>0.) then
        if (grd%id_u_iceberg>0) &
          grd%u_iceberg(i,j)=grd%u_iceberg(i,j)/grd%mass(i,j)
        if (grd%id_v_iceberg>0) &
          grd%v_iceberg(i,j)=grd%v_iceberg(i,j)/grd%mass(i,j)
      else
        if (grd%id_u_iceberg>0)  grd%u_iceberg(i,j)=0.
        if (grd%id_v_iceberg>0)  grd%v_iceberg(i,j)=0.
      endif
    enddo; enddo
  endif

  !Calculating ustar_iceberg (gridded)
  grd%ustar_iceberg(:,:)=0.
  if  ((grd%id_ustar_iceberg>0) .or. (bergs%pass_fields_to_ocean_model)) then   !Update diagnostic of iceberg mass spread on ocean
    do j = grd%jsc,grd%jec ; do i = grd%isc,grd%iec
      dvo=sqrt((grd%spread_uvel(i,j)-grd%uo(i,j))**2+(grd%spread_vvel(i,j)-grd%vo(i,j))**2)
      ustar = sqrt(bergs%cdrag_icebergs*(dvo**2  + bergs%utide_icebergs**2))
      ustar_h = max(bergs%ustar_icebergs_bg, ustar)
      if (grd%spread_area(i,j) ==0.0) ustar_h=0.
        grd%ustar_iceberg(i,j)=ustar_h
    enddo; enddo
  endif

  !Only allowing melt in ocean above a minimum cutoff thickness
  if (bergs%apply_thickness_cutoff_to_gridded_melt) then
    do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
      if ((bergs%melt_cutoff >=0.) .and. (grd%spread_area(i,j)>0.)) then
        ave_thickness=grd%spread_mass(i,j)/(grd%spread_area(i,j)*bergs%rho_bergs)
        ave_draft=ave_thickness*(bergs%rho_bergs/rho_seawater)
        if ((grd%ocean_depth(i,j)-ave_draft) < bergs%melt_cutoff) then
          grd%floating_melt(i,j)=0.0
          grd%calving_hflx(i,j)=0.0
        endif
      endif
    enddo ;enddo
  endif
end subroutine create_gridded_icebergs_fields

!> Calculates basal melt for given thermodynamic properties
subroutine find_basal_melt(bergs, dvo, lat, salt, temp, Use_three_equation_model, thickness, basal_melt, id)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  real, intent(in) :: dvo !< Speed of iceberg relative to ocean mixed layer
  real, intent(in) :: lat !< Latitude (for boundary layer calculation)
  real, intent(in) :: salt !< Salinity of mixed layer
  real, intent(in) :: temp !< Temperature of mixed layer
  logical, intent(in) :: Use_three_equation_model !< True uses the 3 equation model, False uses the 2 equation model.
  real, intent(in) :: thickness !< Ice thickness - needed to work out the pressure below the ice
  real, intent(out) :: basal_melt !< Melt rate underneath the icebergs
  integer(kind=8), intent(in) :: id !< Iceberg id, used for debugging (error messages)
  ! Local variables
  real :: ustar, f_cori, absf,tfreeze
  real :: Hml !Mixed layer depth

  !These could also be useful output variables if needed.
  real :: t_flux, exch_vel_t, exch_vel_s,tflux_shelf,lprec

  real ::  Rhoml   ! Ocean mixed layer density in kg m-3.
  real ::  p_int   ! The pressure at the ice-ocean interface, in Pa.

  real, parameter :: VK    = 0.40     ! Von Karman's constant - dimensionless
  real :: ZETA_N = 0.052   ! The fraction of the boundary layer over which the
                           ! viscosity is linearly increasing. (Was 1/8. Why?)
  real, parameter :: RC    = 0.20     ! critical flux Richardson number.
  real :: I_ZETA_N  ! The inverse of ZETA_N.
  real :: I_LF  ! Inverse of Latent Heat of fusion (J kg-1)
  real :: I_VK      ! The inverse of VK.
  real :: PR, SC    ! The Prandtl number and Schmidt number, nondim.

  ! 3 equation formulation variables
  real :: Sbdry     !   Salinities in the ocean at the interface with the
  real :: Sbdry_it  ! the ice shelf, in PSU.
  real :: dS_it     ! The interface salinity change during an iteration, in PSU.
  real :: hBL_neut  ! The neutral boundary layer thickness, in m.
  real :: hBL_neut_h_molec ! The ratio of the neutral boundary layer thickness
                           ! to the molecular boundary layer thickness, ND.
  real :: wT_flux ! The vertical fluxes of heat and buoyancy just inside the
  real :: wB_flux ! ocean, in C m s-1 and m2 s-3, ###CURRENTLY POSITIVE UPWARD.
  real :: dB_dS  ! The derivative of buoyancy with salinity, in m s-2 PSU-1.
  real :: dB_dT  ! The derivative of buoyancy with temperature, in m s-2 C-1.
  real :: I_n_star, n_star_term
  real :: dIns_dwB  ! The partial derivative of I_n_star with wB_flux, in ???.
  real :: dT_ustar, dS_ustar
  real :: ustar_h
  real :: Gam_turb
  real :: Gam_mol_t, Gam_mol_s
  real :: RhoCp
  real :: I_RhoLF
  real :: Rho0
  real :: ln_neut
  real :: mass_exch
  real :: Sb_min, Sb_max
  real :: dS_min, dS_max
  real :: density_ice

  ! Variables used in iterating for wB_flux.
  real :: wB_flux_new, DwB, dDwB_dwB_in
  real :: I_Gam_T, I_Gam_S
  real :: dG_dwB, iDens
  logical :: Sb_min_set, Sb_max_set
  logical :: out_of_bounds

  real, parameter :: c2_3 = 2.0/3.0
  integer ::  it1, it3

  !Parameters copied ice shelf module defaults (could be entered in the namelist later)
  real, parameter :: dR0_dT = -0.038357 ! Partial derivative of the mixed layer density with temperature, in units of kg m-3 K-1.
  real, parameter :: dR0_dS = 0.805876 ! Partial derivative of the mixed layer density with salinity, in units of kg m-3 psu-1.
  real, parameter :: RHO_T0_S0 = 999.910681 ! Density of water with T=0, S=0 for linear EOS
  real, parameter :: Salin_Ice =0.0 !Salinity of ice
  real, parameter :: Temp_Ice = -15.0 !Salinity of ice
  real, parameter :: kd_molec_salt=  8.02e-10 !The molecular diffusivity of salt in sea water at the freezing point
  real, parameter :: kd_molec_temp=  1.41e-7 !The molecular diffusivity of heat in sea water at the freezing point
  real, parameter :: kv_molec=  1.95e-6 !The molecular molecular kinematic viscosity of sea water at the freezing point
  real, parameter :: Cp_Ice =  2009.0 !Specific heat capacity of ice, taking from HJ99 (Holland and Jenkins 1999)
  real, parameter :: Cp_ml =  3974.0 !Specific heat capacity of mixed layer, taking from HJ99 (Holland and Jenkins 1999)
  real, parameter :: LF =  3.335e5 !Latent heat of fusion, taken from HJ99 (Holland and Jenkins 1999)
  real, parameter :: gamma_t =  0.0 ! Exchange velocity used in 2 equation model. Whn gamma_t is >0, the exchange velocity is independent of u_star.
                                  ! When gamma_t=0.0, then gamma_t is not used, and the exchange velocity is found using u_star.
  real, parameter :: p_atm =  101325 ! Average atmospheric pressure (Pa) - from Google.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  density_ice = bergs%rho_bergs
  Rho0=rho_seawater  !Note that the ice shelf code has a default of Rho0=1035
  Hml =10.      !Mixed layer depth. This is an approximate value. It looks like the code is not sensitive to it (since it enters in log(Hml)
  p_int= p_atm+(gravity*thickness*density_ice)    ! The pressure at the ice-ocean interface, in Pa.

  ! Find the ocean mixed layer density in kg m-3.
  call calculate_density(temp, salt, p_int, Rhoml, Rho_T0_S0, dR0_dT, dR0_dS)

  ! This routine finds the melt at the base of the icebergs using the 2 equation
  ! model or 3 equation model. This code is adapted from the ice shelf code. Once
  ! the iceberg model is inside the ocean model, we should use the same code.

  I_ZETA_N = 1.0 / ZETA_N
  I_RhoLF = 1.0/(Rho0*LF)
  I_LF = 1.0 / LF
  SC = kv_molec/kd_molec_salt
  PR = kv_molec/kd_molec_temp
  I_VK = 1.0/VK
  RhoCp = Rho0 * Cp_ml

  !first calculate molecular component
  Gam_mol_t = 12.5 * (PR**c2_3) - 6
  Gam_mol_s = 12.5 * (SC**c2_3) - 6

  iDens = 1.0/Rho0

  !Preparing the mixed layer properties for use in both 2 and 3 equation version
  ustar = sqrt(bergs%cdrag_icebergs*(dvo**2  + bergs%utide_icebergs**2))
  ustar_h = max(bergs%ustar_icebergs_bg, ustar)

  ! Estimate the neutral ocean boundary layer thickness as the minimum of the
  ! reported ocean mixed layer thickness and the neutral Ekman depth.
  !(Note that in Dan's code, f is spread over adjacent grid cells)
  if ((bergs%grd%grid_is_latlon) .and. (.not. bergs%use_f_plane)) then
     f_cori=(2.*omega)*sin(pi_180*lat)
  else
     f_cori=(2.*omega)*sin(pi_180*bergs%lat_ref)
  endif
  absf = abs(f_cori)  !Absolute value of the Coriolis parameter
  if ((absf*Hml <= VK*ustar_h) .or. (absf.eq.0.))  then
    hBL_neut = Hml
  else
    hBL_neut = (VK*ustar_h) / absf
  endif
  hBL_neut_h_molec = ZETA_N * ((hBL_neut * ustar_h) / (5.0 * Kv_molec))
  ln_neut = 0.0 ; if (hBL_neut_h_molec > 1.0) ln_neut = log(hBL_neut_h_molec)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (Use_three_equation_model) then   ! Use 3 equation model
    ! 3 equation model solves for the melt rates iteratively. This is not working right now, because we don't have access to the mixed layer
    ! temperature and salinity gradients

    ! Guess sss as the iteration starting point for the boundary salinity.
    Sbdry = salt ; Sb_max_set = .false. ; Sb_min_set = .false.
    out_of_bounds=.false.

    ! Determine the mixed layer buoyancy flux, wB_flux.
    dB_dS = (gravity / Rhoml) * dR0_dS
    dB_dT = (gravity / Rhoml) * dR0_dT

    do it1 = 1,20
      ! Determine the potential temperature at the ice-ocean interface.
      call calculate_TFreeze(Sbdry, p_int, tfreeze)

      dT_ustar = (temp - tfreeze) * ustar_h
      dS_ustar = (salt - Sbdry) * ustar_h

      ! First, determine the buoyancy flux assuming no effects of stability
      ! on the turbulence.  Following H & J '99, this limit also applies
      ! when the buoyancy flux is destabilizing.

      if (bergs%const_gamma) then ! if using a constant gamma_T
        I_Gam_T = bergs%Gamma_T_3EQ
        I_Gam_S = bergs%Gamma_T_3EQ/35.
      else
        Gam_turb = I_VK * (ln_neut + (0.5 * I_ZETA_N - 1.0))
        I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
        I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
      endif
      wT_flux = dT_ustar * I_Gam_T
      wB_flux = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * wT_flux

      if (wB_flux > 0.0) then
        ! The buoyancy flux is stabilizing and will reduce the tubulent
        ! fluxes, and iteration is required.
        n_star_term = (ZETA_N/RC) * (hBL_neut * VK) / ustar_h**3
        do it3 = 1,30
          ! n_star <= 1.0 is the ratio of working boundary layer thickness
          ! to the neutral thickness.
          ! hBL = n_star*hBL_neut ; hSub = 1/8*n_star*hBL
          I_n_star = sqrt(1.0 + n_star_term * wB_flux)
          dIns_dwB = 0.5 * n_star_term / I_n_star
          if (hBL_neut_h_molec > I_n_star**2) then
            Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + &
                       (0.5*I_ZETA_N*I_n_star - 1.0))
            dG_dwB =  I_VK * ( -2.0 / I_n_star + (0.5 * I_ZETA_N)) * dIns_dwB
          else
            !   The layer dominated by molecular viscosity is smaller than
            ! the assumed boundary layer.  This should be rare!
            Gam_turb = I_VK * (0.5 * I_ZETA_N*I_n_star - 1.0)
            dG_dwB = I_VK * (0.5 * I_ZETA_N) * dIns_dwB
          endif

          if (bergs%const_gamma) then ! if using a constant gamma_T
            I_Gam_T = bergs%Gamma_T_3EQ
            I_Gam_S = bergs%Gamma_T_3EQ/35.
          else
            I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
            I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
          endif

          wT_flux = dT_ustar * I_Gam_T
          wB_flux_new = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * wT_flux

          ! Find the root where dwB = 0.0
          DwB = wB_flux_new - wB_flux
          if (abs(wB_flux_new - wB_flux) < &
            1e-4*(abs(wB_flux_new) + abs(wB_flux))) exit

          dDwB_dwB_in = -dG_dwB * (dB_dS * (dS_ustar * I_Gam_S**2) + &
                                         dB_dT * (dT_ustar * I_Gam_T**2)) - 1.0
          ! This is Newton's method without any bounds. ( ### SHOULD BOUNDS BE NEEDED?)
          wB_flux_new = wB_flux - DwB / dDwB_dwB_in
        enddo !it3
      endif

      t_flux  = RhoCp * wT_flux
      exch_vel_t = ustar_h * I_Gam_T
      exch_vel_s = ustar_h * I_Gam_S

      if (t_flux <= 0.0) then  ! Freezing occurs, so zero ice heat flux.
        lprec = I_LF * t_flux
        tflux_shelf = 0.0
      else
      !no conduction/perfect insulator
      tflux_shelf = 0.0
      lprec = I_LF * t_flux
      ! With melting, from H&J 1999, eqs (31) & (26)...
      !   Q_ice ~= cp_ice * (Temp_Ice-T_freeze) * lprec
      !   RhoLF*lprec = Q_ice + t_flux
      !   lprec = (t_flux) / (LF + cp_ice * (T_freeze-Temp_Ice))
      !   lprec = t_flux /  (LF + Cp_ice * (tfreeze - Temp_Ice))
      !   tflux_shelf = t_flux - LF*lprec
      !other options: dTi/dz linear through shelf
      !            dTi_dz = (Temp_Ice - tfreeze)/draft
      !            tflux_shelf = - Rho_Ice * Cp_ice * KTI * dTi_dz
      endif

      mass_exch = exch_vel_s * Rho0
      Sbdry_it = (salt * mass_exch + Salin_Ice * lprec) / (mass_exch + lprec)
      dS_it = Sbdry_it - Sbdry
      if (abs(dS_it) < 1e-4*(0.5*(salt + Sbdry + 1.e-10))) exit

      if (dS_it < 0.0) then ! Sbdry is now the upper bound.
        if (Sb_max_set .and. (Sbdry > Sb_max)) then
          if (debug) then
            call error_mesg('KID,Find basal melt', 'shelf_calc_flux: Irregular iteration for Sbdry (max).' ,WARNING)
            print *, 'Sbdry error: id,dvo,temp,salt,lat,thickness :',id,dvo,temp,salt,lat,thickness
          endif
          out_of_bounds=.true.
          exit
        endif
        Sb_max = Sbdry ; dS_max = dS_it ; Sb_max_set = .true.
      else ! Sbdry is now the lower bound.
        if (Sb_min_set .and. (Sbdry < Sb_min)) then
          if (debug) then
            call error_mesg('KID,Find basal melt', 'shelf_calc_flux: Irregular iteration for Sbdry (min).' ,WARNING)
            print *, 'Sbdry error: id,dvo,temp,salt,lat,thickness :',id,dvo,temp,salt,lat,thickness
          endif
          out_of_bounds=.true.
          exit
        endif
        Sb_min = Sbdry ; dS_min = dS_it ; Sb_min_set = .true.
      endif
      if (Sb_min_set .and. Sb_max_set) then
        ! Use the false position method for the next iteration.
        Sbdry = Sb_min + (Sb_max-Sb_min) * (dS_min / (dS_min - dS_max))
      else
        Sbdry = Sbdry_it
      endif
      Sbdry = Sbdry_it
    enddo !it1
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((.not. Use_three_equation_model) .or. (out_of_bounds)) then
    ! In the 2-equation form, the mixed layer turbulent exchange velocity
    ! is specified and large enough that the ocean salinity at the interface
    ! is about the same as the boundary layer salinity.
    ! Alon: I have adapted the code so that the turbulent exchange velocoty is not constant, but rather proportional to the frictional velocity.
    ! This should give you the same answers as the 3 equation model when salinity gradients in the mixed layer are zero (I think/hope)
    ! Use 2-equation model when 3 equation version fails.

    call calculate_TFreeze(salt, p_int, tfreeze)

    Gam_turb = I_VK * (ln_neut + (0.5 * I_ZETA_N - 1.0))
    I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)

    exch_vel_t= ustar_h * I_Gam_T
    if (gamma_t>0.0) exch_vel_t = gamma_t  !Option to set the exchange to a constant, independent of the frictional velocity (as was previously coded)
    wT_flux = exch_vel_t *(temp - tfreeze)

    t_flux  = RhoCp * wT_flux
    tflux_shelf = 0.0
    lprec = I_LF * t_flux
  endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! melt in m/s (melts of ice melted per second)
  basal_melt = lprec /density_ice

end subroutine find_basal_melt

!> Calculates freezing point potential temperature of seawater using a linear relation
!!
!! This subroutine computes the freezing point potential temperature
!! (in deg C) from salinity (in psu), and pressure (in Pa) using a simple
!! linear expression, with coefficients passed in as arguments.
!!
!! Copied from subroutine calculate_TFreeze_linear_scalar (in MOM/equation_of_state)
subroutine calculate_TFreeze(S, pres, T_Fr)
  ! Arguments
  real, intent(in) :: S !< Salinity (1e-3)
  real, intent(in) :: pres !< Presure (Pa)
  real, intent(out) :: T_Fr !< Freezing point (C)
  ! Local variables
  real, parameter :: dTFr_dp    = -7.53E-08    !DTFREEZE_DP in MOM_input
  real, parameter :: dTFr_dS    = -0.0573      !DTFREEZE_DS in MOM_input
  real, parameter :: TFr_S0_P0  =0.0832        !TFREEZE_S0_P0 in MOM_input
  ! TFr_S0_P0 - The freezing point at S=0, p=0, in deg C.
  ! dTFr_dS - The derivatives of freezing point with salinity, in deg C PSU-1.
  ! dTFr_dp - The derivatives of freezing point with pressure, in deg C Pa-1.
  T_Fr = (TFr_S0_P0 + dTFr_dS*S) + dTFr_dp*pres
end subroutine calculate_TFreeze

!> Calculates density of seawater using a linear equation of state
!!
!! This subroutine computes the density of sea water with a trivial
!! linear equation of state (in kg/m^3) from salinity (sal in psu),
!! potential temperature (T in deg C), and pressure in Pa.
!!
!! Copied from subroutine calculate_density_scalar_linear (in MOM/equation_of_state)
subroutine calculate_density(T, S, pressure, rho, Rho_T0_S0, dRho_dT, dRho_dS)
  !Arguments
  real, intent(in)  :: T !< Potential temperature (C)
  real, intent(in)  :: S !< Salinity (1e-3)
  real, intent(in)  :: pressure !< Pressure (Pa)
  real, intent(out) :: rho !< In situ density (kg/3)
  real, intent(in)  :: Rho_T0_S0 !< Density at T=0, S=0 (kg/m3)
  real, intent(in)  :: dRho_dT !< Derivative of density w.r.t. potential temperature (kg/m3/C)
  real, intent(in)  :: dRho_dS !< Derivative of density w.r.t. salinity (1e3 kg/m3)
  rho = Rho_T0_S0 + dRho_dT*T + dRho_dS*S
end subroutine calculate_density

!> Returns orientation of a berg determined by its bonds
subroutine find_orientation_using_iceberg_bonds(grd, berg, orientation)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  type(iceberg), pointer :: berg !< Berg for which orientation is needed
  real, intent(inout) :: orientation !< Angle of orientation (radians)
  ! Local variables
  type(iceberg), pointer :: other_berg
  type(bond), pointer :: current_bond
  real :: angle, lat1,lat2,lon1,lon2,dlat,dlon
  real :: r_dist_x, r_dist_y
  real :: lat_ref, dx_dlon, dy_dlat
  real :: theta, bond_count, Average_angle

  bond_count=0.
  Average_angle=0.
  !Don't check orientation of the edges of halo,  since they can contain unassosiated bonds  (this is why halo width must be larger >= 2 to use bonds)
  if  (  ((berg%ine .gt.  grd%isd) .and. (berg%ine .lt. grd%ied)) .and. ((berg%jne .ge.  grd%jsd) .and. (berg%jne .le. grd%jed) ) ) then
    current_bond=>berg%first_bond
    lat1=berg%lat
    lon1=berg%lon
    do while (associated(current_bond)) ! loop over all bonds
      other_berg=>current_bond%other_berg
      if (.not. associated(other_berg)) then !good place for debugging
        !One valid option: current iceberg is on the edge of halo, with other berg on the next pe (not influencing mass spreading)
        !print *, 'Iceberg bond details:',berg%id, current_bond%other_id,berg%halo_berg, mpp_pe()
        !print *, 'Iceberg bond details2:',berg%ine, berg%jne, current_bond%other_berg_ine, current_bond%other_berg_jne
        !print *, 'Iceberg isd,ied,jsd,jed:',grd%isd, grd%ied, grd%jsd, grd%jed
        !print *, 'Iceberg isc,iec,jsc,jec:',grd%isc, grd%iec, grd%jsc, grd%jec
        !call error_mesg('KID,calculating orientation', 'Looking at bond interactions of unassosiated berg!' ,FATAL)
        !endif
      else
        lat2=other_berg%lat
        lon2=other_berg%lon

        dlat=lat2-lat1
        dlon=lon2-lon1

        lat_ref=0.5*(lat1+lat2)
        call convert_from_grid_to_meters(lat_ref,grd%grid_is_latlon,dx_dlon,dy_dlat)
        r_dist_x=dlon*dx_dlon
        r_dist_y=dlat*dy_dlat

        if (r_dist_x .eq. 0.) then
          angle=pi/2.
        else
          angle=atan(r_dist_y/r_dist_x)
          angle= ((pi/2.)  - (orientation*(pi/180.)))  - angle
          !print *, 'angle: ', angle*(180/pi), initial_orientation
          angle=modulo(angle ,pi/3.)
        endif
        bond_count=bond_count+1.
        Average_angle=Average_angle+angle
      endif
      current_bond=>current_bond%next_bond
    enddo  !End loop over bonds
    if (bond_count.gt.0) then
      Average_angle =Average_angle/bond_count
    else
      Average_angle =0.
    endif
    orientation=modulo(Average_angle ,pi/3.)
  endif

end subroutine find_orientation_using_iceberg_bonds

!> Spread mass of a berg around cells centered on i,j
subroutine spread_mass_across_ocean_cells(bergs, berg, i, j, x, y, Mberg, Mbits, scaling, Area, Tn, addfootloose)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Berg whose mass is being considered
  integer, intent(in) :: i !< i-index of cell contained center of berg
  integer, intent(in) :: j !< j-index of cell contained center of berg
  real, intent(in) :: x !< Longitude of berg (degree E)
  real, intent(in) :: y !< Latitude of berg (degree N)
  real, intent(in) :: Mberg !< Mass of berg (kg)
  real, intent(in) :: Mbits !< Combined mass of bergy and footloose bergy bits (kg)
  real, intent(in) :: scaling !< Multiplier to scale mass (nondim)
  real, intent(in) :: Area !< Area of berg (m2)
  real, intent(in) :: Tn !< Thickness of berg (m)
  logical, intent(in) :: addfootloose !< Add contributions of footloose bergs
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: xL, xC, xR, yD, yC, yU, Mass, L
  real :: yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR
  real :: S, H, origin_x, origin_y, x0, y0
  real :: Area_Q1,Area_Q2 , Area_Q3,Area_Q4, Area_hex
  real :: fraction_used !fraction of iceberg mass included (part of the mass near the boundary is discarded sometimes)
  real :: I_fraction_used !Inverse of fraction used
  real :: tol
  real :: Dn, Hocean
  real, parameter :: rho_seawater=1035.
  integer :: stderrunit
  logical :: debug
  real :: orientation, Mass_berg
  real :: Mfl,Lfl,Wfl,Tfl,Mbits_fl

  ! Get the stderr unit number
  stderrunit = stderr()

  tol=1.e-10
  grd=>bergs%grd
  Mass_berg=Mberg

  if (addfootloose) then
    Mfl=berg%mass_of_fl_bits
    Mbits_fl=berg%mass_of_fl_bergy_bits
  else
    Mfl=0.
    Mbits_fl=0.
  endif

  ! Trimming icebergs to account for grounded fraction.
  if (bergs%grounding_fraction>0.) then
    Hocean=bergs%grounding_fraction*(grd%ocean_depth(i,j)+grd%ssh(i,j))
    Dn=(bergs%rho_bergs/rho_seawater)*Tn ! re-calculate draught (keel depth)
    if (Dn>Hocean) Mass_berg=Mass_berg*min(1.,Hocean/Dn)
    if (Mfl>0. .and. addfootloose) then
      call fl_bits_dimensions(bergs,berg,Lfl,Wfl,Tfl)
      Dn=(bergs%rho_bergs/rho_seawater)*Tfl ! re-calculate draught (keel depth) for FL bits
      if (Dn>Hocean) Mfl=Mfl*min(1.,Hocean/Dn)
    endif
  endif

  Mass_berg=Mass_berg+Mfl

  Mass=(Mass_berg+Mbits+Mbits_fl)*scaling
  ! This line attempts to "clip" the weight felt by the ocean. The concept of
  ! clipping is non-physical and this step should be replaced by grounding.
  if (grd%clipping_depth>0.) Mass=min(Mass,grd%clipping_depth*grd%area(i,j)*rho_seawater)

  !Initialize weights for each cell
  yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0.
  yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

  if (.not. bergs%hexagonal_icebergs) then ! Treat icebergs as rectangles of size L: (this is the default)

    ! L is the non dimensional length of the iceberg [ L=(Area of berg/ Area of grid cell)^0.5 ] or something like that.
    if (grd%area(i,j)>0) then
      L=min( sqrt(Area / grd%area(i,j)),1.0)
    else
      L=1.
    endif

    if (bergs%use_old_spreading) then
      ! Old version before icebergs were given size L
      xL=min(0.5, max(0., 0.5-x))
      xR=min(0.5, max(0., x-0.5))
      xC=max(0., 1.-(xL+xR))
      yD=min(0.5, max(0., 0.5-y))
      yU=min(0.5, max(0., y-0.5))
      yC=max(0., 1.-(yD+yU))
    else
      xL=min(0.5, max(0., 0.5-(x/L)))
      xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
      xC=max(0., 1.-(xL+xR))
      yD=min(0.5, max(0., 0.5-(y/L)))
      yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
      yC=max(0., 1.-(yD+yU))
    endif

    yDxL=yD*xL*grd%msk(i-1,j-1)
    yDxC=yD*xC*grd%msk(i  ,j-1)
    yDxR=yD*xR*grd%msk(i+1,j-1)
    yCxL=yC*xL*grd%msk(i-1,j  )
    yCxR=yC*xR*grd%msk(i+1,j  )
    yUxL=yU*xL*grd%msk(i-1,j+1)
    yUxC=yU*xC*grd%msk(i  ,j+1)
    yUxR=yU*xR*grd%msk(i+1,j+1)
    yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

    fraction_used=1. ! rectangular bergs do share mass with boundaries (all mass is included in cells)

  else ! Spread mass as if elements area hexagonal

    orientation=bergs%initial_orientation
    if ((bergs%iceberg_bonds_on) .and. (bergs%rotate_icebergs_for_mass_spreading)) call find_orientation_using_iceberg_bonds(grd,berg,orientation)

    if (grd%area(i,j)>0) then
      H=min(( (sqrt(Area/(2.*sqrt(3.))) / sqrt(grd%area(i,j)))),1.) ! Non-dimensionalize element length by grid area. (This gives the non-dim Apothem of the hexagon)
    else
      H=(sqrt(3.)/2)*(0.49) ! Largest allowable H, since this makes S=0.49, and S has to be less than 0.5 (Not sure what the implications of this are)
    endif
    S=(2/sqrt(3.))*H !Side of the hexagon

    if (S>0.5) then
      ! The width of an iceberg should not be greater than half the grid cell, or else it can spread over 3 cells  (i.e. S must be less than 0.5 non-dimensionally)
      !print 'Elements must be smaller than a whole grid cell', 'i.e.: S= ' , S , '>=0.5'
      call error_mesg('KID, hexagonal spreading', 'Diameter of the iceberg is larger than a grid cell. Use smaller icebergs', WARNING)
    endif

    !Subtracting the position of the nearest corner from x,y  (The mass will then be spread over the 4 cells connected to that corner)
    origin_x=1. ; origin_y=1.
    if (x<0.5) origin_x=0.
    if (y<0.5) origin_y=0.

    !Position of the hexagon center, relative to origin at the nearest vertex
    x0=(x-origin_x)
    y0=(y-origin_y)

    call Hexagon_into_quadrants_using_triangles(x0,y0,H,orientation,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)

    if (min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <-tol) then
      call error_mesg('KID, hexagonal spreading', 'Intersection with hexagons should not be negative!!!', WARNING)
      write(stderrunit,*) 'KID, yU,yC,yD', Area_Q1, Area_Q2, Area_Q3, Area_Q4
    endif

    Area_Q1=Area_Q1/Area_hex
    Area_Q2=Area_Q2/Area_hex
    Area_Q3=Area_Q3/Area_hex
    Area_Q4=Area_Q4/Area_hex

    !Now, you decide which quadrant belongs to which mass on ocean cell.
    if ((x.ge. 0.5) .and. (y.ge. 0.5)) then !Top right vertex
      yUxR=Area_Q1
      yUxC=Area_Q2
      yCxC=Area_Q3
      yCxR=Area_Q4
    elseif ((x .lt. 0.5) .and. (y.ge. 0.5)) then  !Top left vertex
      yUxC=Area_Q1
      yUxL=Area_Q2
      yCxL=Area_Q3
      yCxC=Area_Q4
    elseif ((x.lt.0.5) .and. (y.lt. 0.5)) then !Bottom left vertex
      yCxC=Area_Q1
      yCxL=Area_Q2
      yDxL=Area_Q3
      yDxC=Area_Q4
    elseif ((x.ge.0.5) .and. (y.lt. 0.5)) then!Bottom right vertex
      yCxR=Area_Q1
      yCxC=Area_Q2
      yDxC=Area_Q3
      yDxR=Area_Q4
    endif

    !Temporary for debugging reasons.
    if (mpp_pe()==mpp_root_pe()) then
      !write(stderrunit,*) 'KID, You are in the hexagonal domain now!!!'
    endif

    !Double check that all the mass is being used.
    if ((abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>tol) .and. (mpp_pe().eq. mpp_root_pe())) then
      !call error_mesg('KID, hexagonal spreading', 'All the mass is not being used!!!', WARNING)
      write(stderrunit,*) 'KID, hexagonal, H,x0,y0', H, x0 , y0
      write(stderrunit,*) 'KID, hexagonal, Areas',(Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      debug=.True.
      !call Hexagon_into_quadrants_using_triangles(x0,y0,H,orientation,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4, debug)
      call error_mesg('KID, hexagonal spreading', 'All the mass is not being used!!!', FATAL)
    endif

    !Scale each cell by (1/fraction_used) in order to redisribute ice mass which landed up on the land, back into the ocean
    !Note that for the square elements, the mass has already been reassigned, so fraction_used shoule be equal to 1 aready
    fraction_used= ((yDxL*grd%msk(i-1,j-1)) + (yDxC*grd%msk(i  ,j-1))  +(yDxR*grd%msk(i+1,j-1)) +(yCxL*grd%msk(i-1,j  )) +  (yCxR*grd%msk(i+1,j  ))&
                   +(yUxL*grd%msk(i-1,j+1)) +(yUxC*grd%msk(i  ,j+1))   +(yUxR*grd%msk(i+1,j+1)) + (yCxC**grd%msk(i,j)))
    if  (berg%static_berg .eq. 1)  fraction_used=1.  !Static icebergs do not share their mass with the boundary
                                                ! (this allows us to easily  initialize hexagonal icebergs in regular arrangements against boundaries)
  endif
  I_fraction_used=1./fraction_used !Invert this so that the arithmatec reprocudes

  !Spreading the iceberg mass onto the ocean
  call spread_variable_across_cells(grd, grd%mass_on_ocean, Mass, i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)
  !Spreading the iceberg area onto the ocean
  call spread_variable_across_cells(grd, grd%area_on_ocean, Area*scaling , i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR,I_fraction_used)
  !Spreading the iceberg x momentum onto the ocean
  call spread_variable_across_cells(grd,grd%Uvel_on_ocean, berg%uvel*Area*scaling , i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)
  !Spreading the iceberg y momentum onto the ocean
  call spread_variable_across_cells(grd,grd%Vvel_on_ocean, berg%vvel*Area*scaling , i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)

end subroutine spread_mass_across_ocean_cells

!> Distribute a quantity among nine cells on a grid centered at cell i,j
subroutine spread_variable_across_cells(grd, variable_on_ocean, Var, i, j, &
           yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)
  ! Arguments
  type(icebergs_gridded), pointer, intent(in) :: grd !< Container for gridded fields
  real, dimension(grd%isd:grd%ied, grd%jsd:grd%jed, 9), intent(inout) :: variable_on_ocean !< Gridded field to augment
  real, intent(in) :: Var !< Variable to be spread accross cell
  real, intent(in) :: yDxL !< Weight for the cell at i-1,j-1
  real, intent(in) :: yDxC !< Weight for the cell at i-1,j
  real, intent(in) :: yDxR !< Weight for the cell at i-1,j+1
  real, intent(in) :: yCxL !< Weight for the cell at i,j-1
  real, intent(in) :: yCxC !< Weight for the cell at i,j
  real, intent(in) :: yCxR !< Weight for the cell at i,j-1
  real, intent(in) :: yUxL !< Weight for the cell at i+1,j-1
  real, intent(in) :: yUxC !< Weight for the cell at i+1,j
  real, intent(in) :: yUxR !< Weight for the cell at i+1,j+1
  real, intent(in) :: I_fraction_used !< Amount of iceberg used (inverse)
  integer, intent(in) :: i !< i-index of cell containing center of berg
  integer, intent(in) :: j !< j-index of cell containing center of berg

  !Spreading the iceberg mass onto the ocean
  variable_on_ocean(i,j,1)=variable_on_ocean(i,j,1)+(yDxL*Var*I_fraction_used)
  variable_on_ocean(i,j,2)=variable_on_ocean(i,j,2)+(yDxC*Var*I_fraction_used)
  variable_on_ocean(i,j,3)=variable_on_ocean(i,j,3)+(yDxR*Var*I_fraction_used)
  variable_on_ocean(i,j,4)=variable_on_ocean(i,j,4)+(yCxL*Var*I_fraction_used)
  variable_on_ocean(i,j,5)=variable_on_ocean(i,j,5)+(yCxC*Var*I_fraction_used)
  variable_on_ocean(i,j,6)=variable_on_ocean(i,j,6)+(yCxR*Var*I_fraction_used)
  variable_on_ocean(i,j,7)=variable_on_ocean(i,j,7)+(yUxL*Var*I_fraction_used)
  variable_on_ocean(i,j,8)=variable_on_ocean(i,j,8)+(yUxC*Var*I_fraction_used)
  variable_on_ocean(i,j,9)=variable_on_ocean(i,j,9)+(yUxR*Var*I_fraction_used)

end subroutine spread_variable_across_cells

!> Returns area of a triangle
real function Area_of_triangle(Ax, Ay, Bx, By, Cx, Cy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  Area_of_triangle    =   abs(    0.5*((Ax*(By-Cy))+(Bx*(Cy-Ay))+(Cx*(Ay-By))) )
end function Area_of_triangle

!> Returns x rounded of to sig_fig
!! \todo What the heck is this for? -AJA
real function roundoff(x,sig_fig)
  ! Arguments
  real, intent(in) :: x !< A quantity with 15 significant figures of useful information
  integer, intent(in) :: sig_fig !< Number of significant figures to keep
  !roundoff=round(x*(10**(sig_fig))
  roundoff=(FLOAT(INT(x * (10.**sig_fig) + 0.5)) / (10.**sig_fig))
end function roundoff

!> Returns true of a point is in or on the rectangle with opposite corners A and B
logical function point_in_interval(Ax, Ay, Bx, By, px, py)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: px !< x-position of point
  real, intent(in) :: py !< y-position of point
  point_in_interval=.False.
  if ((px <= max(Ax,Bx)) .and. (px >= min(Ax,Bx))) then
    if ((py <= max(Ay,By)) .and. (py >= min(Ay,By))) then
      point_in_interval=.True.
    endif
  endif
end function point_in_interval

!> Returns true if point q is on a line through points A and B
logical function point_is_on_the_line(Ax, Ay, Bx, By, qx, qy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of point A
  real, intent(in) :: Ay !< y-position of point A
  real, intent(in) :: Bx !< x-position of point B
  real, intent(in) :: By !< y-position of point B
  real, intent(in) :: qx !< x-position of point q
  real, intent(in) :: qy !< y-position of point q
  ! Local variables
  real :: tol, dxc,dyc,dxl,dyl,cross
  !tol=1.e-12
  tol=0.0
  dxc = qx - Ax
  dyc = qy - Ay
  dxl = Bx - Ax
  dyl = By - Ay
  cross = dxc * dyl - dyc * dxl
  if (abs(cross)<=tol) then
    point_is_on_the_line=.True.
  else
   point_is_on_the_line=.False.
  endif
end function point_is_on_the_line

!> Returns True if a point q is inside a triangle ABC
!!
!! This function decides whether a point (qx,qy) is inside the triangle ABC.
!! There is also the option to include the boundary of the triangle.
logical function point_in_triangle(Ax, Ay, Bx, By, Cx, Cy, qx, qy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  real, intent(in) :: qx !< x-position of point q
  real, intent(in) :: qy !< y-position of point q
  ! Local variables
  real :: l0,l1,l2,p0,p1,p2
  real :: v0x,v1x,v2x,v0y,v1y,v2y,dot00,dot01,dot02,dot11,dot12

  point_in_triangle = .False.
  if ((Ax==qx .and. Ay==qy) .or. (Bx==qx .and. By==qy) .or. (Cx==qx .and. Cy==qy)) then ! Exclude the pathelogical case
      point_in_triangle = .False.
  else
    if (((point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) .or. (point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy))) .or. (point_is_on_the_line(Bx,By,Cx,Cy,qx,qy)))) then
      point_in_triangle = .False.
    else
      ! Compute point in triangle using Barycentric coordinates (the same as sum_sign_dot_prod routines)
      l0=(qx-Ax)*(By-Ay)-(qy-Ay)*(Bx-Ax)
      l1=(qx-Bx)*(Cy-By)-(qy-By)*(Cx-Bx)
      l2=(qx-Cx)*(Ay-Cy)-(qy-Cy)*(Ax-Cx)

      p0=sign(1., l0); if (l0==0.)  p0=0.
      p1=sign(1., l1); if (l1==0.)  p1=0.
      p2=sign(1., l2); if (l2==0.)  p2=0.

      if ( (abs(p0)+abs(p2))+(abs(p1)) == abs((p0+p2)+(p1)) )  point_in_triangle = .True.
    endif
  endif
end function point_in_triangle

!> Calculates the two areas of a triangle divided by an axis line
!!
!! This function calculates the area of a triangle on opposite sides of an axis when the
!! triangle is split with two points on one side, and one point on the other.
!! In this function, A is the point on one side of the axis, and B,C are on the opposite sides.
!! \todo You should change this name a little, so that it not similar the other routine.
subroutine Area_of_triangle_across_axes(Ax, Ay, Bx, By, Cx, Cy, axis1, Area_positive, Area_negative)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  character, intent(in) :: axis1 !< Either 'x' or 'y'
  real, intent(out) :: Area_positive !< Area on negative side of axis line
  real, intent(out) :: Area_negative !< Area on positive side of axis line
  ! Local variables
  real :: pABx, pABy, pACx, pACy, A0
  real :: A_half_triangle, A_triangle

  A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)

  call intercept_of_a_line(Ax,Ay,Bx,By,axis1,pABx, pABy)
  call intercept_of_a_line(Ax,Ay,Cx,Cy,axis1,pACx, pACy)

  if (axis1=='x')  A0=Ay; !Value used for if statements (deciding up/down vs left/right)
  if (axis1=='y')  A0=Ax; !Value used for if statements (deciding up/down vs left/right)

  A_half_triangle=Area_of_triangle(Ax,Ay,pABx,pABy,pACx,pACy)
  if (A0>=0.) then
    Area_positive= A_half_triangle
    Area_negative= A_triangle-A_half_triangle
  else
    Area_positive= A_triangle-A_half_triangle
    Area_negative= A_half_triangle
  endif

end subroutine Area_of_triangle_across_axes

!> Returns the axis intercept of a line AB
!!
!! This routine returns the position (x0,y0) at which a line AB intercepts the x or y axis.
!! The value No_intercept_val is returned when the line does not intercept the axis.
subroutine intercept_of_a_line(Ax, Ay, Bx, By, axes1, x0, y0)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  character, intent(in) :: axes1 !< Either 'x' or 'y'
  real, intent(out) :: x0 !< x-position of intercept
  real, intent(out) :: y0 !< y-position of intercept
  ! Local variables
  real :: No_intercept_val ! Huge value used to make sure that the intercept is outside the triangle in the parallel case.

  No_intercept_val=100000000000. ! Huge value used to make sure that the intercept is outside the triangle in the parallel case.
  x0=No_intercept_val
  y0=No_intercept_val

  if (axes1=='x') then ! x intercept
    if (Ay.ne.By) then
      x0=Ax -(((Ax-Bx)/(Ay-By))*Ay)
      y0=0.
    endif
  endif

  if (axes1=='y') then ! y intercept
    if (Ax.ne.Bx) then
      x0=0.
      y0=-(((Ay-By)/(Ax-Bx))*Ax)+Ay
    endif
  endif
end subroutine intercept_of_a_line

!> Calculates the area of a triangle on either side of an axis, if any.
!!
!! This routine gives you the area of a triangle on opposite sides of the axis specified.
!! It also takes care of the special case where the triangle is totally on one side.
!! This routine calls Area_of_triangle_across_axes to calculate the areas when the triangles are split.
subroutine divding_triangle_across_axes(Ax, Ay, Bx, By, Cx, Cy, axes1, Area_positive, Area_negative)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  character, intent(in) :: axes1 !< Either 'x' or 'y'
  real, intent(out) :: Area_positive !< Area on negative side of axis line
  real, intent(out) :: Area_negative !< Area on positive side of axis line
  ! Local variables
  real :: A0,B0,C0
  real A_triangle

  if (axes1=='x') then ! Use the y-coordinates for if statements to see which side of the line you are on
    A0=Ay
    B0=By
    C0=Cy
  endif
  if (axes1=='y') then ! Use the y-coordinates for if statements to see which side of the line you are on
    A0=Ax
    B0=Bx
    C0=Cx
  endif

  A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)
  if ((B0*C0)>0.) then ! B and C are on the same side  (and non-zero)
    if ((A0*B0).ge.0.) then ! all three on the same side (if it equals zero, then A0=0 and the others are not)
      if ((A0>0.)  .or.  ((A0==0.) .and.  (B0>0.))) then
        Area_positive= A_triangle
        Area_negative= 0.
      else
        Area_positive= 0.
        Area_negative= A_triangle
      endif
    else  !A is on the opposite side to B and C
      call Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative)
    endif

  elseif ((B0*C0)<0.) then !B and C are on the opposite sides
    if ((A0*B0).ge. 0.) then !C is all alone
      call Area_of_triangle_across_axes(Cx,Cy,Bx,By,Ax,Ay,axes1,Area_positive, Area_negative)
    else !B is all alone
      call Area_of_triangle_across_axes(Bx,By,Cx,Cy,Ax,Ay,axes1,Area_positive, Area_negative)
    endif

  else  !This is the case when either B or C is equal to zero (or both), A0 could be zero too.
    if (((A0.eq.0.) .and. (B0.eq.0.)) .and. (C0.eq.0.)) then
      Area_positive= 0.
      Area_negative= 0.
    elseif ((A0*B0<0.)  .or.  (A0*C0<0.)) then    !A, B are on opposite sides, and C is zero.  OR  A, C are on opposite sides, and B is zero.
      call Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative)
    elseif (((A0*B0>0.) .or. (A0*C0>0.)) .or. (((abs(A0)>0.) .and. (B0==0.)) .and. (C0==0.))) then
      if (A0>0.) then
        Area_positive= A_triangle
        Area_negative= 0.
      else
        Area_positive= 0.
        Area_negative= A_triangle
      endif

    elseif (A0.eq. 0.) then   !(one of B,C is zero too)
      if ((B0>0.) .or. (C0>0.)) then
        Area_positive= A_triangle
        Area_negative= 0.
      elseif ((B0<0.) .or. (C0<0.)) then
        Area_positive= 0.
        Area_negative= A_triangle
      else
        call error_mesg('KID, iceberg_run', 'Logical error inside triangle dividing routine', FATAL)
      endif
    else
      call error_mesg('KID, iceberg_run', 'Another logical error inside triangle dividing routine', FATAL)
    endif
  endif
end subroutine divding_triangle_across_axes

!> Areas of a triangle divided into quadrants
!!
!! This routine takes a triangle, and finds the intersection with the four quadrants.
subroutine Triangle_divided_into_four_quadrants(Ax, Ay, Bx, By, Cx, Cy, Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  real, intent(out) :: Area_triangle !< Are of triangle
  real, intent(out) :: Area_Q1 !< Are in quadrant 1
  real, intent(out) :: Area_Q2 !< Are in quadrant 2
  real, intent(out) :: Area_Q3 !< Are in quadrant 2
  real, intent(out) :: Area_Q4 !< Are in quadrant 4
  ! Local variables
  real :: Area_Upper, Area_Lower, Area_Right, Area_Left
  real :: px, py , qx , qy
  real :: Area_key_quadrant,Error
  real :: tol
  integer :: Key_quadrant
  integer ::sig_fig
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()
  tol=1.e-10

  Area_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)

  ! Calculating area across axes
  call divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'x',Area_Upper ,Area_Lower)
  call divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'y',Area_Right ,Area_Left)

  ! Decide if the origin is in the triangle. If so, then you have to divide the area 4 ways
  ! This is done by finding a quadrant where the intersection between the triangle and quadrant forms a new triangle
  ! (This occurs when on of the sides of the triangle  intersects both the x and y axis)
  if (point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.)) then
    ! Find a line in the triangle that cuts both axes in/on the triangle
    call intercept_of_a_line(Ax,Ay,Bx,By,'x',px,py); !x_intercept
    call intercept_of_a_line(Ax,Ay,Bx,By,'y',qx,qy); !y_intercept
    ! Note that the 1. here means that we include points on the boundary of the triangle.
    if (.not.((point_in_interval(Ax,Ay,Bx,By,px,py)) .and. (point_in_interval(Ax,Ay,Bx,By,qx,qy)))) then
      call intercept_of_a_line(Ax,Ay,Cx,Cy,'x',px,py); !x_intercept
      call intercept_of_a_line(Ax,Ay,Cx,Cy,'y',qx,qy); !y_intercept
      if (.not.((point_in_interval(Ax,Ay,Cx,Cy,px,py)) .and. (point_in_interval(Ax,Ay,Cx,Cy,qx,qy)))) then
        call intercept_of_a_line(Bx,By,Cx,Cy,'x',px,py); !x_intercept
        call intercept_of_a_line(Bx,By,Cx,Cy,'y',qx,qy); !y_intercept
        if (.not.((point_in_interval(Bx,By,Cx,Cy,px,py)) .and. (point_in_interval(Bx,By,Cx,Cy,qx,qy)))) then
          ! You should not get here, but there might be some bugs in the code to do with points exactly falling on axes.
          !if (mpp_pe().eq.12) then
            write(stderrunit,*) 'KID,corners', Ax,Ay,Bx,By,Cx,Cy
          !endif
          call error_mesg('KID, iceberg_run', 'Something went wrong with Triangle_divide_into_four_quadrants', FATAL)
        endif
      endif
    endif

    ! Assigning quadrants. Key_quadrant is the quadrant with the baby triangle in it.
    Area_key_quadrant=Area_of_triangle(px,py,qx,qy,0.,0.)
    if ((px.ge. 0.) .and. (qy.ge. 0.)) then  !First quadrant
      Key_quadrant=1
    elseif ((px.lt.0.) .and. (qy.ge. 0.)) then  !Second quadrant
      Key_quadrant=2
    elseif ((px.lt. 0.) .and. (qy.lt. 0.)) then !Third quadrant
      Key_quadrant=3
    elseif ((px.ge. 0.) .and. (qy.lt. 0.)) then !Forth quadrant
      Key_quadrant=4
    else  !
      call error_mesg('KID, iceberg_run', 'None of the quadrants are Key', WARNING)
      write(stderrunit,*) 'KID, Triangle, px,qy', px,qy
    endif

  else ! At least one quadrant is empty, and this can be used to find the areas in the other quadrant.  Assigning quadrants. Key_quadrant is the empty quadrant.
    Area_key_quadrant=0
    if      ( (.not. ((((Ax>0.) .and. (Ay>0.)) .or. ((Bx>0.) .and. (By> 0.))) .or. ((Cx>0.) .and. (Cy> 0.)))) .and. ((Area_Upper+Area_Right).le.Area_triangle) ) then
      ! No points land in this quadrant and triangle does not cross the quadrant
      Key_quadrant=1
    elseif  ( (.not. ((((Ax<0.) .and. (Ay>0)) .or. ((Bx<0.) .and. (By>0.))) .or. ((Cx<0.) .and. (Cy>0.)))) .and. ((Area_Upper+Area_Left).le. Area_triangle) ) then
      Key_quadrant=2
    elseif  ( (.not. ((((Ax<0.) .and. (Ay<0.)) .or. ((Bx<0.) .and. (By< 0.))) .or. ((Cx<0.) .and. (Cy< 0.)))) .and. ((Area_Lower+Area_Left) .le.Area_triangle) ) then
      Key_quadrant=3
    else
      Key_quadrant=4
    endif
  endif

  ! Assign values to quadrants
  if (Key_quadrant .eq. 1) then
    Area_Q1=Area_key_quadrant
    Area_Q2=Area_Upper-Area_Q1
    Area_Q4=Area_Right-Area_Q1
    !Area_Q3=Area_Left-Area_Q2 ! These lines have been changes so that the sum of the 4 quadrants exactly matches the triangle area.
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4)
  elseif (Key_quadrant .eq. 2) then
    Area_Q2=Area_key_quadrant
    Area_Q1=Area_Upper-Area_Q2
    Area_Q4=Area_Right-Area_Q1
    !Area_Q3=Area_Left-Area_Q2
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4)
  elseif (Key_quadrant==3) then
    Area_Q3=Area_key_quadrant
    Area_Q2=Area_Left-Area_Q3
    Area_Q1=Area_Upper-Area_Q2
    !Area_Q4=Area_Right-Area_Q1
    Area_Q4=Area_triangle-(Area_Q1+Area_Q2+Area_Q3)
  elseif (Key_quadrant==4) then
    Area_Q4=Area_key_quadrant
    Area_Q1=Area_Right-Area_Q4
    Area_Q2=Area_Upper-Area_Q1
    !Area_Q3=Area_Left-Area_Q2
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4)
  else
    call error_mesg('KID, iceberg_run', 'Logical error inside triangle into four quadrants. Should not get here.', FATAL)
  endif

  Area_Q1=max(Area_Q1,0.)
  Area_Q2=max(Area_Q2,0.)
  Area_Q3=max(Area_Q3,0.)
  Area_Q4=max(Area_Q4,0.)


  Error=abs(Area_Q1+Area_Q2+Area_Q3+Area_Q4-Area_triangle)
  if (Error>tol) then
    call error_mesg('KID, triangle spreading', 'Triangle not evaluated accurately!!', WARNING)
    !if (mpp_pe().eq.mpp_root_pe()) then
    if (mpp_pe().eq. 20) then
      write(stderrunit,*) 'KID, Triangle corners:',Ax,Ay,Bx,By,Cx,Cy
      write(stderrunit,*) 'KID, Triangle, Full Area', Area_Q1+ Area_Q2+ Area_Q3+ Area_Q4
      write(stderrunit,*) 'KID, Triangle, Areas', Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      write(stderrunit,*) 'KID, Triangle, Areas', Error
      write(stderrunit,*) 'KID, Key quadrant',Key_quadrant,Area_key_quadrant
      write(stderrunit,*) 'KID, point in triangle',(point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))
      write(stderrunit,*) 'KID, halves',Area_Upper,Area_Lower,Area_Right,Area_Left
    endif
  endif

end subroutine Triangle_divided_into_four_quadrants

!> Rotates a point clockwise about origin and then translates by x0,y0
subroutine rotate_and_translate(px, py, theta, x0, y0)
  ! Arguments
  real, intent(in) :: x0 !< x-direction shift
  real, intent(in) :: y0 !< y-direction shift
  real, intent(in) :: theta !< Angle to rotate (degrees)
  real, intent(inout) :: px !< x-coordinate of point
  real, intent(inout) :: py !< y-coordinate of point
  ! Local variables
  real :: px_temp,py_temp

  ! Rotation
  px_temp = ( cos(theta*pi/180)*px) + (sin(theta*pi/180)*py)
  py_temp = (-sin(theta*pi/180)*px) + (cos(theta*pi/180)*py)

  ! Translation
  px= px_temp + x0
  py= py_temp + y0
end subroutine rotate_and_translate

!> Areas of a hexagon divided into quadrants
!!
!! This subroutine divides a regular hexagon centered at x0,y0 with apothem H, and orientation theta into its intersection with the 4 quadrants.
!! Theta=0 assumes that the apothem points upwards.
!! Routine works by finding the corners of the 6 triangles, and then finding the intersection of each of these with each quadrant.
!! \todo (also the rotation is not working yet)
subroutine Hexagon_into_quadrants_using_triangles(x0, y0, H, theta, Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  ! Arguments
  real, intent(in) :: x0 !< x-coordinate of center of hexagon
  real, intent(in) :: y0 !< y-coordinate of center of hexagon
  real, intent(in) :: H !< Apothem (inner radius of hexagon)
  real, intent(in) :: theta !< Orientation angle of hexagon
  real, intent(out) :: Area_hex !< Area of hexagon
  real, intent(out) :: Area_Q1 !< Are in quadrant 1
  real, intent(out) :: Area_Q2 !< Are in quadrant 2
  real, intent(out) :: Area_Q3 !< Are in quadrant 2
  real, intent(out) :: Area_Q4 !< Are in quadrant 4
  ! Local variables
  real :: C1x, C2x, C3x, C4x, C5x, C6x
  real :: C1y, C2y, C3y, C4y, C5y, C6y
  real :: T12_Area, T12_Q1, T12_Q2, T12_Q3, T12_Q4
  real :: T23_Area, T23_Q1, T23_Q2, T23_Q3, T23_Q4
  real :: T34_Area, T34_Q1, T34_Q2, T34_Q3, T34_Q4
  real :: T45_Area, T45_Q1, T45_Q2, T45_Q3, T45_Q4
  real :: T56_Area, T56_Q1, T56_Q2, T56_Q3, T56_Q4
  real :: T61_Area, T61_Q1, T61_Q2, T61_Q3, T61_Q4
  real :: S, exact_hex_area, Error
  real :: tol
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()
  tol=1.e-10

  ! Length of side of Hexagon
  S=(2/sqrt(3.))*H

  ! Finding positions of corners
  C1x=S           ; C1y=0.  !Corner 1 (right)
  C2x=H/sqrt(3.)  ; C2y=H;  !Corner 2 (top right)
  C3x=-H/sqrt(3.) ; C3y=H;  !Corner 3 (top left)
  C4x=-S          ; C4y=0.; !Corner 4 (left)
  C5x=-H/sqrt(3.) ; C5y=-H; !Corner 5 (bottom left)
  C6x=H/sqrt(3.)  ; C6y=-H; !Corner 6 (bottom right)

  ! Finding positions of corners
  call rotate_and_translate(C1x,C1y,theta,x0,y0)
  call rotate_and_translate(C2x,C2y,theta,x0,y0)
  call rotate_and_translate(C3x,C3y,theta,x0,y0)
  call rotate_and_translate(C4x,C4y,theta,x0,y0)
  call rotate_and_translate(C5x,C5y,theta,x0,y0)
  call rotate_and_translate(C6x,C6y,theta,x0,y0)

  ! Area of Hexagon is the sum of the triangles
  call Triangle_divided_into_four_quadrants(x0,y0,C1x,C1y,C2x,C2y,T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4); !Triangle 012
  call Triangle_divided_into_four_quadrants(x0,y0,C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4); !Triangle 023
  call Triangle_divided_into_four_quadrants(x0,y0,C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4); !Triangle 034
  call Triangle_divided_into_four_quadrants(x0,y0,C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4); !Triangle 045
  call Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4); !Triangle 056
  call Triangle_divided_into_four_quadrants(x0,y0,C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4); !Triangle 061

  ! Summing up the triangles
  Area_hex=T12_Area+T23_Area+T34_Area+T45_Area+T56_Area+T61_Area
  Area_Q1=T12_Q1+T23_Q1+T34_Q1+T45_Q1+T56_Q1+T61_Q1
  Area_Q2=T12_Q2+T23_Q2+T34_Q2+T45_Q2+T56_Q2+T61_Q2
  Area_Q3=T12_Q3+T23_Q3+T34_Q3+T45_Q3+T56_Q3+T61_Q3
  Area_Q4=T12_Q4+T23_Q4+T34_Q4+T45_Q4+T56_Q4+T61_Q4

  Area_Q1=max(Area_Q1,0.)
  Area_Q2=max(Area_Q2,0.)
  Area_Q3=max(Area_Q3,0.)
  Area_Q4=max(Area_Q4,0.)

  Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
  if ((abs(Error)>tol))then
    if (mpp_pe().eq.mpp_root_pe()) then
      call error_mesg('KID, hexagonal spreading', 'Hexagon error is large!!', WARNING)
      write(stderrunit,*) 'KID, hex error, H,x0,y0, Error', H, x0 , y0, Error
      write(stderrunit,*) 'KID, hex error, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      write(stderrunit,*) 'KID, Triangle1',C1x,C1y,C2x,C2y,T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4,(T12_Q1+T12_Q2+T12_Q3+T12_Q4-T12_Area)
      write(stderrunit,*) 'KID, Triangle2',C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4,(T23_Q1+T23_Q2+T23_Q3+T23_Q4-T23_Area)
      write(stderrunit,*) 'KID, Triangle3',C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4,(T34_Q1+T34_Q2+T34_Q3+T34_Q4-T34_Area)
      write(stderrunit,*) 'KID, Triangle4',C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4,(T45_Q1+T45_Q2+T45_Q3+T45_Q4-T45_Area)
      write(stderrunit,*) 'KID, Triangle5',C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4,(T56_Q1+T56_Q2+T56_Q3+T56_Q4-T56_Area)
      write(stderrunit,*) 'KID, Triangle6',C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4,(T61_Q1+T61_Q2+T61_Q3+T61_Q4-T61_Area)
    endif
  endif

  exact_hex_area=((3.*sqrt(3.)/2)*(S*S))
  if (abs(Area_hex-exact_hex_area)>tol) then
    call error_mesg('KID, hexagonal spreading', 'Hexagon not evaluated accurately!!', WARNING)
    if (mpp_pe().eq.mpp_root_pe()) then
      write(stderrunit,*) 'KID, hex calculations, H,x0,y0', H, x0 , y0
      write(stderrunit,*) 'KID, hex calculations, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
    endif
  endif

  ! Adjust Areas so that the error is zero by subtracting the error from the largest sector.
   if  (((Area_Q1>=Area_Q2) .and. (Area_Q1>=Area_Q3)) .and. (Area_Q1>=Area_Q4)) then
     Area_Q1=Area_Q1+Error
   elseif  (((Area_Q2>=Area_Q1) .and. (Area_Q2>=Area_Q3)) .and. (Area_Q2>=Area_Q4)) then
     Area_Q2=Area_Q2+Error
   elseif  (((Area_Q3>=Area_Q1) .and. (Area_Q3>=Area_Q2)) .and. (Area_Q3>=Area_Q4)) then
     Area_Q3=Area_Q3+Error
   elseif  (((Area_Q4>=Area_Q1) .and. (Area_Q4>=Area_Q2)) .and. (Area_Q4>=Area_Q3)) then
     Area_Q4=Area_Q4+Error
   else
     call error_mesg('KID, hexagonal spreading', 'Error in hexagon is larger than any quadrant!!', WARNING)
     if (mpp_pe().eq.mpp_root_pe()) then
      write(stderrunit,*) 'KID, hex quadrants, H,x0,y0', H, x0 , y0, Error
      write(stderrunit,*) 'KID, hex quadrants, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
     endif
   endif

 end subroutine Hexagon_into_quadrants_using_triangles

!> Loop through all bergs and call interp_flds
subroutine interp_gridded_fields_to_bergs(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg
  integer :: grdj,grdi,i,is,ie,js,je
  type(randomNumberStream) :: rns ! Random numbers for stochastic tidal parameterization
  real :: rx,ry

  ! For convenience
  grd=>bergs%grd
  rx=0.0; ry=0.0;

  if (bergs%mts) then
    is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
  else
    is=grd%isc-1; ie=grd%iec+1; js=grd%jsc-1; je=grd%jec+1
  endif

  do grdj = js,je ; do grdi = is,ie ! just for computational domain
    berg=>bergs%list(grdi,grdj)%first
    if (associated(berg) .and. grd%tidal_drift>0.) then
      ! Seed random numbers based on space and "time"
      rns = initializeRandomNumberStream( grdi + 10000*grdj + &
                                          int( 16384.*abs( sin(262144.*grd%ssh(grdi,grdj)) ) ) )
    endif
    do while (associated(berg)) ! loop over all bergs
      if (berg%halo_berg.lt.0.5) then
        if (grd%tidal_drift>0.) then
          call getRandomNumbers(rns, rx)
          rx = 2.*rx - 1.
          call getRandomNumbers(rns, ry)
          ry = 2.*ry - 1.
        endif
        call interp_flds(grd, berg%lon, berg%lat, berg%ine, berg%jne, berg%xi, berg%yj, rx, ry, berg%uo, berg%vo, &
          berg%ui, berg%vi, berg%ua, berg%va, berg%ssh_x, berg%ssh_y, berg%sst, berg%sss, berg%cn, berg%hi, berg%od)
      endif
      berg=>berg%next
    enddo
  enddo;enddo

end subroutine interp_gridded_fields_to_bergs

!> Interpolate ocean, sea ice, and atmosphere fields from grid to iceberg
subroutine interp_flds(grd, x, y, i, j, xi, yj, rx, ry, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, sss, cn, hi, od)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  integer, intent(in) :: i !< i-index of cell in which to interpolate
  integer, intent(in) :: j !< j-index of cell in which to interpolate
  real, intent(in) :: x !<Longitude of position
  real, intent(in) :: y !<Latitude of position
  real, intent(in) :: xi !< Non-dimensional x-position within cell to interpolate to
  real, intent(in) :: yj !< Non-dimensional y-position within cell to interpolate to
  real, intent(in) :: rx !< Random number between -1 and 1 for use in x-component of stochastic tidal parameterization
  real, intent(in) :: ry !< Random number between -1 and 1 for use in y-component of stochastic tidal parameterization
  real, intent(out) :: uo !< Ocean zonal velocity at point xi,yj (m/s)
  real, intent(out) :: vo !< Ocean meridional velocity at point xi,yj (m/s)
  real, intent(out) :: ui !< Ice zonal velocity at point xi,yj (m/s)
  real, intent(out) :: vi !< Ice meridional velocity at point xi,yj (m/s)
  real, intent(out) :: ua !< Atmospheric zonal velocity at point xi,yj (m/s)
  real, intent(out) :: va !< Atmospheric meridional velocity at point xi,yj (m/s)
  real, intent(out) :: ssh_x !< Zonal slope of sea-surface height (nondim)
  real, intent(out) :: ssh_y !< Meridional slope of sea-surface height (nondim)
  real, intent(out) :: sst !< Sea-surface temperature (C)
  real, intent(out) :: sss !< Sea-surface salinity (1e-3)
  real, intent(out) :: cn !< Sea-ice concentration (nondim)
  real, intent(out) :: hi !< Sea-ice thickness (m)
  real, optional, intent(out) :: od !< Ocean depth (m)
  ! Local variables
  real :: cos_rot, sin_rot, du, dv
#ifdef USE_OLD_SSH_GRADIENT
  real :: dxm, dx0, dxp
  real, parameter :: ssh_coast=0.00
#endif
  real :: hxm, hxp, ssh
  integer :: stderrunit
  integer :: ii, jj

  ! Get the stderr unit number
  stderrunit = stderr()

  cos_rot=bilin(grd, grd%cos, i, j, xi, yj) ! If true, uses the inverted bilin function
  sin_rot=bilin(grd, grd%sin, i, j, xi, yj)

  uo=bilin(grd, grd%uo, i, j, xi, yj)
  vo=bilin(grd, grd%vo, i, j, xi, yj)
  ui=bilin(grd, grd%ui, i, j, xi, yj)
  vi=bilin(grd, grd%vi, i, j, xi, yj)
  ua=bilin(grd, grd%ua, i, j, xi, yj)
  va=bilin(grd, grd%va, i, j, xi, yj)

  ! The following block accelerates a berg away from coastlines towards open water
  ! which is a bias needed to avoid piling up in coarse resolution models
  if (grd%coastal_drift > 0.) then
    ! If a cell to the west is land (msk=0), accelerate to the east, and vice versa.
    uo = uo + grd%coastal_drift * ( grd%msk(i+1,j) - grd%msk(i-1,j) ) * grd%msk(i,j)
    ui = ui + grd%coastal_drift * ( grd%msk(i+1,j) - grd%msk(i-1,j) ) * grd%msk(i,j)
    ! If a cell to the south is land (msk=0), accelerate to the north, and vice versa.
    vo = vo + grd%coastal_drift * ( grd%msk(i,j+1) - grd%msk(i,j-1) ) * grd%msk(i,j)
    vi = vi + grd%coastal_drift * ( grd%msk(i,j+1) - grd%msk(i,j-1) ) * grd%msk(i,j)
  endif

  ! Stochastic acceleration to represent unresolved tidal and wave motions.
  ! Note: this scheme should account for the time-step and have memory!
  if (grd%tidal_drift > 0.) then
    ! The acceleration is modulated to not move particles towards land cells.
    du = ( min(0., rx) * grd%msk(i-1,j) + max(0., rx) * grd%msk(i+1,j) ) &
         * ( 1. - grd%msk(i,j-1) * grd%msk(i,j+1) ) ! Do not apply in open ocean
    dv = ( min(0., ry) * grd%msk(i,j-1) + max(0., ry) * grd%msk(i,j+1) ) &
         * ( 1. - grd%msk(i-1,j) * grd%msk(i+1,j) ) ! Do not apply in open ocean
    du = du * grd%tidal_drift * grd%msk(i,j)
    dv = dv * grd%tidal_drift * grd%msk(i,j)
    uo = uo + du
    ui = ui + du
    vo = vo + dv
    vi = vi + dv
  endif

  if (ua.ne.ua) then
    if (mpp_pe().eq.9) then
      write(stderrunit,'(a3,32i7)') 'ua',(ii,ii=grd%isd,grd%ied)
      do jj=grd%jed,grd%jsd,-1
        write(stderrunit,'(i3,32f7.1)') jj,(grd%ua(ii,jj),ii=grd%isd,grd%ied)
      enddo
  !   write(stderrunit,'(a3,32i7)') 'Lat',(i,i=grd%isd,grd%ied)
  !   do j=grd%jed,grd%jsd,-1
  !     write(stderrunit,'(i3,32f7.1)') j,(grd%lat(i,j),i=grd%isd,grd%ied)
  !   enddo
  !   write(stderrunit,'(a3,32i7)') 'Msk',(i,i=grd%isd,grd%ied)
  !   do j=grd%jed,grd%jsd,-1
  !     write(stderrunit,'(i3,32f7.1)') j,(grd%msk(i,j),i=grd%isd,grd%ied)
  !   enddo
  ! endif
      call error_mesg('KID, interp fields', 'ua is NaNs', FATAL)
    endif
  endif

  ! These fields are cell centered (A-grid) and would
  ! best be interpolated using PLM. For now we use PCM!
  sst=grd%sst(i,j) ! A-grid
  sss=grd%sss(i,j) ! A-grid
  cn=grd%cn(i,j) ! A-grid
  hi=grd%hi(i,j) ! A-grid

  ! Estimate SSH gradient in X direction
#ifdef USE_OLD_SSH_GRADIENT
  dxp=0.5*(grd%dx(i+1,j)+grd%dx(i+1,j-1))
  dx0=0.5*(grd%dx(i,j)+grd%dx(i,j-1))
  dxm=0.5*(grd%dx(i-1,j)+grd%dx(i-1,j-1))
  hxm=2.*(grd%ssh(i,j)-grd%ssh(i-1,j))/(dx0+dxm)*grd%msk(i-1,j) &
        +(-ssh_coast)/(dx0+dxm)*(1.-grd%msk(i-1,j)) ! force to drive bergs away from coasts
  hxp=2.*(grd%ssh(i+1,j)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i+1,j) &
        +(+ssh_coast)/(dx0+dxp)*(1.-grd%msk(i+1,j)) ! force to drive bergs away from coasts
#else
  if (yj>=0.5) then
    hxp=(yj-0.5)*ddx_ssh(grd,i  ,j+1)+(1.5-yj)*ddx_ssh(grd,i  ,j  )
    hxm=(yj-0.5)*ddx_ssh(grd,i-1,j+1)+(1.5-yj)*ddx_ssh(grd,i-1,j  )
  else
    hxp=(yj+0.5)*ddx_ssh(grd,i  ,j  )+(0.5-yj)*ddx_ssh(grd,i  ,j-1)
    hxm=(yj+0.5)*ddx_ssh(grd,i-1,j  )+(0.5-yj)*ddx_ssh(grd,i-1,j-1)
  endif
#endif
  ! ssh_x is at the u-point on a C-grid
  ssh_x=xi*hxp+(1.-xi)*hxm

  ! Estimate SSH gradient in Y direction
#ifdef USE_OLD_SSH_GRADIENT
  dxp=0.5*(grd%dy(i,j+1)+grd%dy(i-1,j+1))
  dx0=0.5*(grd%dy(i,j)+grd%dy(i-1,j))
  dxm=0.5*(grd%dy(i,j-1)+grd%dy(i-1,j-1))
  hxm=2.*(grd%ssh(i,j)-grd%ssh(i,j-1))/(dx0+dxm)*grd%msk(i,j-1) &
        +(-ssh_coast)/(dx0+dxm)*(1.-grd%msk(i,j-1)) ! force to drive bergs away from coasts
  hxp=2.*(grd%ssh(i,j+1)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i,j+1) &
        +(+ssh_coast)/(dx0+dxp)*(1.-grd%msk(i,j+1)) ! force to drive bergs away from coasts
#else
  if (xi>=0.5) then
    hxp=(xi-0.5)*ddy_ssh(grd,i+1,j  )+(1.5-xi)*ddy_ssh(grd,i  ,j  )
    hxm=(xi-0.5)*ddy_ssh(grd,i+1,j-1)+(1.5-xi)*ddy_ssh(grd,i  ,j-1)
  else
    hxp=(xi+0.5)*ddy_ssh(grd,i  ,j  )+(0.5-xi)*ddy_ssh(grd,i-1,j  )
    hxm=(xi+0.5)*ddy_ssh(grd,i  ,j-1)+(0.5-xi)*ddy_ssh(grd,i-1,j-1)
  endif
#endif
  ! ssh_y is at the v-point on a C-grid
  ssh_y=yj*hxp+(1.-yj)*hxm

  ! Rotate vectors from local grid to lat/lon coordinates
  call rotate(uo, vo, cos_rot, sin_rot)
  call rotate(ui, vi, cos_rot, sin_rot)
  call rotate(ua, va, cos_rot, sin_rot)
  call rotate(ssh_x, ssh_y, cos_rot, sin_rot)

  !There are some issues with the boundaries ssh gradient calculation in a finite domain. This is a temporary fix
  if (ssh_x.ne.ssh_x) ssh_x=0.
  if (ssh_y.ne.ssh_y) ssh_y=0.

  if (((((uo.ne.uo) .or. (vo.ne.vo)) .or. ((ui.ne.ui) .or. (vi.ne.vi))) .or. &
       (((ua.ne.ua) .or. (va.ne.va)) .or. ((ssh_x.ne.ssh_x) .or. (ssh_y.ne.ssh_y)))) .or. &
       (((sst.ne. sst) .or. (sss.ne. sss) .or. (cn.ne.cn)) .or. (hi.ne. hi))) then
    write(stderrunit,*) 'KID, Error in interpolate: uo,vo,ui,vi',uo, vo, ui, vi
    write(stderrunit,*) 'KID, Error in interpolate: ua,va,ssh_x,ssh_y', ua, va, ssh_x, ssh_y
    write(stderrunit,*) 'KID, Error in interpolate: sst,cn,hi', sst, sss, cn, hi, mpp_pe()
    call error_mesg('KID, interp fields', 'field interpaolations has NaNs', FATAL)

  endif

  ! Quadratic interpolation of ocean depth+sea surface height (A-grid)
  ! Ocean depth is only needed for grounding friction, used with the MTS model only?
  if (present(od)) then
    if (mts) then
      od=quad_interp_from_agrid(grd,grd%ocean_depth,x,y,i,j,xi,yj)+quad_interp_from_agrid(grd,grd%ssh,x,y,i,j,xi,yj)
    else
      od=grd%ocean_depth(i,j)+grd%ssh(i,j)
    endif
  endif
end subroutine interp_flds

!> Returns zonal slope of sea-surface height across the east face of cell i,j
real function ddx_ssh(grd,i,j)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  integer, intent(in) :: i !< i-index of cell
  integer, intent(in) :: j !< j-index of cell
  ! Local variables
  real :: dxp,dx0
  dxp=0.5*(grd%dx(i+1,j)+grd%dx(i+1,j-1))
  dx0=0.5*(grd%dx(i,j)+grd%dx(i,j-1))
  ddx_ssh=2.*(grd%ssh(i+1,j)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i+1,j)*grd%msk(i,j)
end function ddx_ssh

!> Returns meridional slope of sea-surface height across the northern face of cell i,j
real function ddy_ssh(grd,i,j)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  integer, intent(in) :: i !< i-index of cell
  integer, intent(in) :: j !< j-index of cell
  ! Local variables
  real :: dyp,dy0
  dyp=0.5*(grd%dy(i,j+1)+grd%dy(i-1,j+1))
  dy0=0.5*(grd%dy(i,j)+grd%dy(i-1,j))
  ddy_ssh=2.*(grd%ssh(i,j+1)-grd%ssh(i,j))/(dy0+dyp)*grd%msk(i,j+1)*grd%msk(i,j)
end function ddy_ssh

! Rotates vector (u,v) using rotation matrix with elements cos_rot and sin_rot
subroutine rotate(u, v, cos_rot, sin_rot)
  ! Arguments
  real, intent(inout) :: u !< x-component of vector
  real, intent(inout) :: v !< y-component of vector
  real, intent(in) :: cos_rot !< Cosine of rotation angle
  real, intent(in) :: sin_rot !< Sine of rotation angle
  ! Local variables
  real :: u_old, v_old

  u_old=u
  v_old=v
  u=cos_rot*u_old+sin_rot*v_old
  v=cos_rot*v_old-sin_rot*u_old

end subroutine rotate

!> Calculates bergs%grd%mass_on_ocean
subroutine calculate_mass_on_ocean(bergs, with_diagnostics)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  logical, intent(in) :: with_diagnostics !< If true, calculate diagnostics
  ! Local variables
  type(iceberg), pointer :: berg
  type(icebergs_gridded), pointer :: grd
  integer :: grdj, grdi
  integer :: j, i

  ! For convenience
  grd=>bergs%grd

  !Initialize fields
  grd%mass_on_ocean(:,:,:)=0.
  grd%area_on_ocean(:,:,:)=0.
  grd%Uvel_on_ocean(:,:,:)=0.
  grd%Vvel_on_ocean(:,:,:)=0.

  do grdj = grd%jsc-1,grd%jec+1 ; do grdi = grd%isc-1,grd%iec+1
    berg=>bergs%list(grdi,grdj)%first
    do while(associated(berg))
      if (berg%halo_berg<2 .or. .not. bergs%mts) then
        i=berg%ine  ;     j=berg%jne
        if (grd%area(i,j) > 0.) then

          !Increasing Mass on ocean
          if ((bergs%add_weight_to_ocean .and. .not. bergs%time_average_weight) .or.(bergs%find_melt_using_spread_mass)) then
            call spread_mass_across_ocean_cells(bergs, berg, berg%ine, berg%jne, berg%xi, berg%yj, berg%mass,&
              berg%mass_of_bits, berg%mass_scaling, berg%length*berg%width, berg%thickness,addfootloose=.true.)
          endif

          !Calculated some iceberg diagnositcs
          if (with_diagnostics) call calculate_sum_over_bergs_diagnositcs(bergs,grd,berg,i,j)

        endif
      endif
      berg=>berg%next
    enddo
  enddo ;enddo

end subroutine calculate_mass_on_ocean

!> Projects additional diagnostics of bergs on to the grid
subroutine calculate_sum_over_bergs_diagnositcs(bergs, grd, berg, i, j)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  type(iceberg), pointer :: berg !< An iceberg
  integer, intent(in) :: i !< i-index of cell containing berg
  integer, intent(in) :: j !< j-index of cell containing berg
  ! Local variables
  real ::  Abits, Abits_fl, Abits_fl_bergy, Lbits, Mbits
  real :: L_fl,W_fl,T_fl

  !Virtual area diagnostic
  if (grd%id_virtual_area>0) then
    if (bergs%bergy_bit_erosion_fraction>0.) then
      Lbits=min(berg%length,berg%width,berg%thickness,40.) ! assume bergy bits are smallest dimension or 40 meters
      Abits=(berg%mass_of_bits/bergs%rho_bergs)/Lbits ! Effective bottom area (assuming T=Lbits)
    else
      Abits=0.0
    endif
    if (bergs%fl_style.eq.'fl_bits') then
      call fl_bits_dimensions(bergs,berg,L_fl,W_fl,T_fl)
      Abits_fl=(berg%mass_of_fl_bits/bergs%rho_bergs)/T_fl ! Effective bottom area
      if (bergs%bergy_bit_erosion_fraction>0.) then
        Lbits=min(L_fl,W_fl,T_fl,40.)
        Abits_fl_bergy=(berg%mass_of_fl_bergy_bits/bergs%rho_bergs)/Lbits
      else
        Abits_fl_bergy=0.0
      endif
    else
      Abits_fl=0.0
      Abits_fl_bergy=0.0
    endif
    grd%virtual_area(i,j)=grd%virtual_area(i,j)+&
      (berg%width*berg%length+Abits+Abits_fl+Abits_fl_bergy)*berg%mass_scaling ! m^2
  endif

  !Mass diagnostic (also used in u_iceberg, v_iceberg
  if ((grd%id_mass>0 ) .or. ((grd%id_u_iceberg>0) .or. (grd%id_v_iceberg>0)))   &
       & grd%mass(i,j)=grd%mass(i,j)+berg%mass/grd%area(i,j)*berg%mass_scaling ! kg/m2

  !Finding the average iceberg velocity in a grid cell (mass weighted)
  if (grd%id_u_iceberg>0) &
  grd%u_iceberg(i,j)=grd%u_iceberg(i,j)+((berg%mass/grd%area(i,j)*berg%mass_scaling)*berg%uvel) ! kg/m2
  if (grd%id_v_iceberg>0) &
  grd%v_iceberg(i,j)=grd%v_iceberg(i,j)+((berg%mass/grd%area(i,j)*berg%mass_scaling)*berg%vvel) ! kg/m2

  !Mass of bergy bits
  if (grd%id_bergy_mass>0 .or. bergs%add_weight_to_ocean)&
    & grd%bergy_mass(i,j)=grd%bergy_mass(i,j)+(berg%mass_of_bits+berg%mass_of_fl_bergy_bits)/grd%area(i,j)*berg%mass_scaling ! kg/m2

  !Mass of footloose bits
  if (grd%id_fl_bits_mass>0 .or. bergs%add_weight_to_ocean)&
    & grd%fl_bits_mass(i,j)=grd%fl_bits_mass(i,j)+berg%mass_of_fl_bits/grd%area(i,j)*berg%mass_scaling ! kg/m2

  !Mass of footloose bits
  if (grd%id_fl_bergy_bits_mass>0 .or. bergs%add_weight_to_ocean)&
    & grd%fl_bergy_bits_mass(i,j)=grd%fl_bergy_bits_mass(i,j)+berg%mass_of_fl_bergy_bits/grd%area(i,j)*berg%mass_scaling ! kg/m2
end subroutine calculate_sum_over_bergs_diagnositcs

!> The main driver the steps updates icebergs
subroutine icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, calving_hflx, cn, hi, &
                        stagger, stress_stagger, sss, mass_berg, ustar_berg, area_berg)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  real, dimension(:,:), intent(inout) :: calving !< Calving (kg/s). This field is updated with melt by bergs.
  real, dimension(:,:), intent(inout) :: calving_hflx !< Calving heat flux (W/m2)
  real, dimension(:,:), intent(in) :: uo !< Ocean zonal velocity (m/s)
  real, dimension(:,:), intent(in) :: vo !< Ocean meridional velocity (m/s)
  real, dimension(:,:), intent(in) :: ui !< Ice zonal velocity (m/s)
  real, dimension(:,:), intent(in) :: vi !< Ice meridional velocity (m/s)
  real, dimension(:,:), intent(in) :: tauxa !< Zonal wind stress (Pa)
  real, dimension(:,:), intent(in) :: tauya !< Meridional wind stress (Pa)
  real, dimension(:,:), intent(in) :: ssh !< Effective sea-surface height (m)
  real, dimension(:,:), intent(in) :: sst !< Sea-surface temperature (C or K)
  real, dimension(:,:), intent(in) :: cn !< Sea-ice concentration (nondim)
  real, dimension(:,:), intent(in) :: hi !< Sea-ice thickness (m)
  integer, optional, intent(in) :: stagger !< Enumerated value indicating staggering of ocean/ice u,v variables
  integer, optional, intent(in) :: stress_stagger !< Enumerated value indicating staggering of stress variables
  real, dimension(:,:), optional, intent(in) :: sss !< Sea-surface salinity (1e-3)
  real, dimension(:,:), optional, pointer :: mass_berg !< Mass of bergs (kg)
  real, dimension(:,:), optional, pointer :: ustar_berg !< Friction velocity on base of bergs (m/s)
  real, dimension(:,:), optional, pointer :: area_berg !< Area of bergs (m2)
  ! Local variables
  integer :: iyr, imon, iday, ihr, imin, isec, k
  type(icebergs_gridded), pointer :: grd
  logical :: lerr, sample_traj, write_traj, lbudget, lverbose, check_bond_quality
  real :: unused_calving, tmpsum, grdd_berg_mass, grdd_bergy_mass, grdd_fl_bits_mass, grdd_spread_mass, grdd_spread_area
  real :: grdd_u_iceberg, grdd_v_iceberg, grdd_ustar_iceberg, grdd_spread_uvel, grdd_spread_vvel
  integer :: i, j, Iu, ju, iv, Jv, Iu_off, ju_off, iv_off, Jv_off
  real :: mask, max_SST
  real, dimension(:,:), allocatable :: uC_tmp, vC_tmp, uA_tmp, vA_tmp
  integer :: vel_stagger, str_stagger
  real, dimension(:,:), allocatable :: iCount
  integer :: nbonds
  integer :: stderrunit
  logical :: Visited=.false.
  save :: Visited

  ! Get the stderr unit number
  stderrunit = stderr()

  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_int)

  vel_stagger = BGRID_NE ; if (present(stagger)) vel_stagger = stagger
  str_stagger = vel_stagger ; if (present(stress_stagger)) str_stagger = stress_stagger

  ! For convenience
  grd=>bergs%grd

  grd%floating_melt(:,:)=0.
  grd%berg_melt(:,:)=0.
  grd%melt_buoy(:,:)=0.
  grd%melt_eros(:,:)=0.
  grd%melt_conv(:,:)=0.
  grd%bergy_src(:,:)=0.
  grd%bergy_melt(:,:)=0.
  grd%bergy_mass(:,:)=0.
  grd%fl_bits_src(:,:)=0.
  grd%fl_bits_melt(:,:)=0.
  grd%fl_bits_mass(:,:)=0.
  grd%fl_bergy_bits_mass(:,:)=0.
  grd%spread_mass_old(:,:)=0.
  !grd%spread_mass(:,:)=0.  !Don't zero this out yet, because we can first use this an add it onto the SSH
  grd%spread_area(:,:)=0.
  grd%u_iceberg(:,:)=0.
  grd%v_iceberg(:,:)=0.
  grd%spread_uvel(:,:)=0.
  grd%spread_vvel(:,:)=0.
  grd%ustar_iceberg(:,:)=0.
  grd%mass(:,:)=0.
  grd%virtual_area(:,:)=0.
  grd%melt_by_class(:,:,:)=0.
  grd%melt_buoy_fl(:,:)=0.
  grd%melt_eros_fl(:,:)=0.
  grd%melt_conv_fl(:,:)=0.
  grd%fl_parent_melt(:,:)=0.
  grd%fl_child_melt(:,:)=0.

  !Initializing _on_ocean_fields
  grd%mass_on_ocean(:,:,:)=0. ;   grd%area_on_ocean(:,:,:)=0.
  grd%Uvel_on_ocean(:,:,:)=0. ;   grd%Vvel_on_ocean(:,:,:)=0.

  if (present(mass_berg)) then ;  if (associated(mass_berg)) then
    mass_berg(:,:)=0.0
  endif ;  endif
  if (present(ustar_berg)) then ; if (associated(ustar_berg)) then
    ustar_berg(:,:)=0.0
  endif ;  endif
  if (present(area_berg)) then ;  if (associated(area_berg)) then
    area_berg(:,:)=0.0
  endif ;  endif

  ! Manage time
  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  bergs%current_year=iyr
  bergs%current_yearday=yearday(imon, iday, ihr, imin, isec)
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if ( (bergs%traj_sample_hrs>0)  .and. (.not. bergs%ignore_traj) ) then
    if (mod(60*60*24*iday+ 60*60*ihr + 60.*imin + isec ,60*60*bergs%traj_sample_hrs).eq.0) &
      sample_traj=.true.
  elseif (bergs%traj_sample_hrs==-1) then
    sample_traj=.true.
  endif
  write_traj=.false.
  if ((bergs%traj_write_hrs>0) .and. (.not. bergs%ignore_traj))  then
     if (mod(60*60*24*iday+ 60*60*ihr + 60.*imin + isec ,60*60*bergs%traj_write_hrs).eq.0) &
       write_traj=.true.
   elseif (bergs%traj_write_hrs==-1) then
     write_traj=.true.
  endif
  lverbose=.false.
  if (bergs%verbose_hrs>0) then
    !if (mod(24*iday+ihr+(imin/60.),float(bergs%verbose_hrs)).eq.0) lverbose=verbose
    if (mod(24*iday+ihr+(imin/60.),bergs%verbose_hrs).eq.0) lverbose=verbose
  endif
  lbudget=.false.
  if (bergs%verbose_hrs>0) then
     !if (mod(24*iday+ihr+(imin/60.),float(bergs%verbose_hrs)).eq.0) lbudget=budget  !Added minutes, so that it does not repeat when smaller time steps are used.
     if (mod(24*iday+ihr+(imin/60.),bergs%verbose_hrs).eq.0) lbudget=budget
  endif
  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'KID: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', bergs%current_year, bergs%current_yearday


 !call sanitize_field(grd%calving,1.e20)
  tmpsum=sum( calving(:,:)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_calving_received=bergs%net_calving_received+tmpsum*bergs%dt

  ! Adapt calving heat flux from coupler
  grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)=calving_hflx(:,:) & ! Units of W/m2
       *grd%msk(grd%isc:grd%iec,grd%jsc:grd%jec)

  ! Adapt calving flux from coupler for use here
  grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:) & ! Units of kg/m2/s
       *grd%msk(grd%isc:grd%iec,grd%jsc:grd%jec)

  ! Running means of calving and calving_hflx
  if (bergs%tau_calving>0.) then
    call get_running_mean_calving(bergs, grd%calving, grd%calving_hflx)
    calving(:,:)=grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) ! Copy back from grd%calving if using running-mean
    calving_hflx(:,:)=grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)
  endif

  grd%calving(:,:)=grd%calving(:,:)*grd%msk(:,:)*grd%area(:,:) ! Convert to kg/s
  tmpsum=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_incoming_calving=bergs%net_incoming_calving+tmpsum*bergs%dt
  if (grd%id_calving>0) &
    lerr=send_data(grd%id_calving, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  grd%calving_hflx(:,:)=grd%calving_hflx(:,:)*grd%msk(:,:) ! Mask (just in case)
  if (grd%id_calving_hflx_in>0) &
    lerr=send_data(grd%id_calving_hflx_in, grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  tmpsum=sum( grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_incoming_calving_heat=bergs%net_incoming_calving_heat+tmpsum*bergs%dt ! Units of J

  if (grd%id_ocean_depth>0) &
    lerr=send_data(grd%id_ocean_depth, grd%ocean_depth(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  if (vel_stagger == BGRID_NE) then
    ! Copy ocean and ice velocities. They are already on B-grid u-points.
    grd%uo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1) = uo(:,:)
    grd%vo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1) = vo(:,:)
    call mpp_update_domains(grd%uo, grd%vo, grd%domain, gridtype=BGRID_NE)
    grd%ui(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1) = ui(:,:)
    grd%vi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1) = vi(:,:)
    call mpp_update_domains(grd%ui, grd%vi, grd%domain, gridtype=BGRID_NE)
  elseif (vel_stagger == CGRID_NE) then
    ! The u- and v- points will have different offsets with symmetric memory.
    Iu_off = (size(uo,1) - (grd%iec - grd%isc))/2 - grd%isc + 1
    ju_off = (size(uo,2) - (grd%jec - grd%jsc))/2 - grd%jsc + 1
    iv_off = (size(vo,1) - (grd%iec - grd%isc))/2 - grd%isc + 1
    Jv_off = (size(vo,2) - (grd%jec - grd%jsc))/2 - grd%jsc + 1
    do I=grd%isc-1,grd%iec ; do J=grd%jsc-1,grd%jec
      ! Interpolate ocean and ice velocities from C-grid velocity points.
      Iu = i + Iu_off ; ju = j + ju_off ; iv = i + iv_off ; Jv = j + Jv_off
      ! This masking is needed for now to prevent icebergs from running up on to land.
      mask = min(grd%msk(i,j), grd%msk(i+1,j), grd%msk(i,j+1), grd%msk(i+1,j+1))
      grd%uo(I,J) = mask * 0.5*(uo(Iu,ju)+uo(Iu,ju+1))
      grd%ui(I,J) = mask * 0.5*(ui(Iu,ju)+ui(Iu,ju+1))
      grd%vo(I,J) = mask * 0.5*(vo(iv,Jv)+vo(iv+1,Jv))
      grd%vi(I,J) = mask * 0.5*(vi(iv,Jv)+vi(iv+1,Jv))
    enddo ; enddo
  else
    call error_mesg('KID, iceberg_run', 'Unrecognized value of stagger!', FATAL)
  endif

  if (str_stagger == BGRID_NE) then
    ! Copy wind stress components on B-grid u-points.
    grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec) = tauxa(:,:)
    grd%va(grd%isc:grd%iec,grd%jsc:grd%jec) = tauya(:,:)
  elseif (str_stagger == CGRID_NE) then
    ! The u- and v- points will have different offsets with symmetric memory.
    Iu_off = (size(tauxa,1) - (grd%iec - grd%isc))/2 - grd%isc + 1
    ju_off = (size(tauxa,2) - (grd%jec - grd%jsc))/2 - grd%jsc + 1
    iv_off = (size(tauya,1) - (grd%iec - grd%isc))/2 - grd%isc + 1
    Jv_off = (size(tauya,2) - (grd%jec - grd%jsc))/2 - grd%jsc + 1
    allocate(uC_tmp(grd%isd:grd%ied,grd%jsd:grd%jed), &
             vC_tmp(grd%isd:grd%ied,grd%jsd:grd%jed))
    uC_tmp(:,:) = 0. ! This avoids uninitialized values that might remain in halo
    vC_tmp(:,:) = 0. ! regions after the call to mpp_update_domains() below.
    !   If the iceberg model used symmetric memory, the starting value of these
    ! copies would need to be decremented by 1.
    do i=grd%isc,grd%iec ; do j=grd%jsc,grd%jec
      uC_tmp(i,j) = tauxa(i+Iu_off, j+ju_off)
      vC_tmp(i,J) = tauya(i+iv_off, J+Jv_off)
    enddo ; enddo

    call mpp_update_domains(uC_tmp, vC_tmp, grd%domain, gridtype=CGRID_NE)
    do I=grd%isc-1,grd%iec ; do J=grd%jsc-1,grd%jec
      ! Interpolate wind stresses from C-grid velocity-points.
      ! This masking is needed for now to prevent icebergs from running up on to land.
      mask = min(grd%msk(i,j), grd%msk(i+1,j), grd%msk(i,j+1), grd%msk(i+1,j+1))
       grd%ua(I,J) = mask * 0.5*(uC_tmp(I,j)+uC_tmp(I,j+1))
       grd%va(I,J) = mask * 0.5*(vC_tmp(i,J)+vC_tmp(i+1,J))
    enddo ; enddo
    deallocate(uC_tmp, vC_tmp)
  elseif (str_stagger == AGRID) then
    ! Copy into arrays with local index conventions and halos.
    allocate(uA_tmp(grd%isd:grd%ied,grd%jsd:grd%jed), &
             vA_tmp(grd%isd:grd%ied,grd%jsd:grd%jed))
    uA_tmp(:,:) = 0.0 ! This avoids uninitialized values that might remain in halo
    vA_tmp(:,:) = 0.0 ! regions after the call to mpp_update_domains() below.
    uA_tmp(grd%isc:grd%iec,grd%jsc:grd%jec) = tauxa(:,:)
    vA_tmp(grd%isc:grd%iec,grd%jsc:grd%jec) = tauya(:,:)
    call mpp_update_domains(uA_tmp, vA_tmp, grd%domain, gridtype=AGRID)
    do I=grd%isc-1,grd%iec ; do J=grd%jsc-1,grd%jec
      ! Interpolate wind stresses from A-grid tracer points to the corner B-grid points.
      ! This masking is needed for now to prevent icebergs from running up on to land.
      mask = min(grd%msk(i,j), grd%msk(i+1,j), grd%msk(i,j+1), grd%msk(i+1,j+1))
      grd%ua(I,J) = mask * 0.25*((uA_tmp(i,j) + uA_tmp(i+1,j+1)) + &
                                 (uA_tmp(i+1,j) + uA_tmp(i,j+1)))
      grd%va(I,J) = mask * 0.25*((vA_tmp(i,j) + vA_tmp(i+1,j+1)) + &
                                 (vA_tmp(i+1,j) + vA_tmp(i,j+1)))
    enddo ; enddo

    deallocate(uA_tmp, vA_tmp)
  else
    call error_mesg('KID, iceberg_run', 'Unrecognized value of stress_stagger!', FATAL)
  endif

  call mpp_update_domains(grd%uo, grd%vo, grd%domain, gridtype=BGRID_NE)
  call mpp_update_domains(grd%ui, grd%vi, grd%domain, gridtype=BGRID_NE)

  call invert_tau_for_du(grd%ua, grd%va) ! Note rough conversion from stress to speed
 !grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec)=sign(sqrt(abs(tauxa(:,:))/0.01),tauxa(:,:))  ! Note rough conversion from stress to speed
 !grd%va(grd%isc:grd%iec,grd%jsc:grd%jec)=sign(sqrt(abs(tauya(:,:))/0.01),tauya(:,:))  ! Note rough conversion from stress to speed
  call mpp_update_domains(grd%ua, grd%va, grd%domain, gridtype=BGRID_NE)

  ! Copy sea surface height and temperature(resides on A grid)
  grd%ssh(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=ssh(:,:)
  if (bergs%add_iceberg_thickness_to_SSH) then
    !We might need to make sure spread_mass is defined on halos (or this might be done automatically. I need to look into this)
    do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
      if (grd%area(i,j)>0) then
        grd%ssh(i,j) =   ((grd%spread_mass(i,j)/grd%area(i,j))*(bergs%rho_bergs/rho_seawater))  !Is this an appropriate sea water density to use? Should be freezing point.
      endif
    enddo ;enddo
  endif

  call mpp_update_domains(grd%ssh, grd%domain)
  max_SST = maxval(sst(:,:))
  if (max_SST > 120.0) then ! The input sst is in degrees Kelvin, otherwise the water would be boiling.
    grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec) = sst(:,:)-273.15 ! Note convert from Kelvin to Celsius
  else  ! The input sst is already in degrees Celsius.
    grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec) = sst(:,:) ! Note no conversion necessary.
  endif
  call mpp_update_domains(grd%sst, grd%domain)
  ! Copy sea-ice concentration and thickness (resides on A grid)
  grd%cn(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=cn(:,:)
  call mpp_update_domains(grd%cn, grd%domain)
  grd%hi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=hi(:,:)
  call mpp_update_domains(grd%hi, grd%domain)

  ! Adding gridded salinity.
  if (present(sss)) then
    grd%sss(grd%isc:grd%iec,grd%jsc:grd%jec)=sss(:,:)
  else
    grd%sss(grd%isc:grd%iec,grd%jsc:grd%jec)=-1.0
    if ((bergs%use_mixed_layer_salinity_for_thermo) .and. (bergs%melt_icebergs_as_ice_shelf))  then
      call error_mesg('KID, icebergs_run', 'Can not use salinity for thermo. Ocean ML salinity not present!', FATAL)
    endif
  endif

  ! Make sure that gridded values agree with mask  (to get rid of NaN values)
  do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
    ! Initializing all gridded values to zero
    if (grd%msk(i,j).lt. 0.5) then
      grd%ua(i,j) = 0.0 ;  grd%va(i,j) = 0.0
      grd%uo(i,j) = 0.0 ;  grd%vo(i,j) = 0.0
      grd%ui(i,j) = 0.0 ;  grd%vi(i,j) = 0.0
      grd%sst(i,j)= 0.0;  grd%sss(i,j)= 0.0
      grd%cn(i,j) = 0.0 ;  grd%hi(i,j) = 0.0
    endif
    if (grd%ua(i,j) .ne. grd%ua(i,j)) grd%ua(i,j)=0.
    if (grd%va(i,j) .ne. grd%va(i,j)) grd%va(i,j)=0.
    if (grd%uo(i,j) .ne. grd%uo(i,j)) grd%uo(i,j)=0.
    if (grd%vo(i,j) .ne. grd%vo(i,j)) grd%vo(i,j)=0.
    if (grd%ui(i,j) .ne. grd%ui(i,j)) grd%ui(i,j)=0.
    if (grd%vi(i,j) .ne. grd%vi(i,j)) grd%vi(i,j)=0.
    if (grd%sst(i,j) .ne. grd%sst(i,j)) grd%sst(i,j)=0.
    if (grd%sss(i,j) .ne. grd%sss(i,j)) grd%sss(i,j)=0.
    if (grd%cn(i,j) .ne. grd%cn(i,j)) grd%cn(i,j)=0.
    if (grd%hi(i,j) .ne. grd%hi(i,j)) grd%hi(i,j)=0.
  enddo; enddo

  if (debug) call bergs_chksum(bergs, 'run bergs (top)')
  if (debug) call checksum_gridded(bergs%grd, 'top of s/r run')

  ! Accumulate ice from calving
  call accumulate_calving(bergs)
  if (grd%id_accum>0) then
    grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:)
    grd%tmp(:,:)=grd%tmp(:,:)*grd%msk(:,:)*grd%area(:,:)-grd%calving(:,:)
    lerr=send_data(grd%id_accum, grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  endif
  if (grd%id_unused>0) &
    lerr=send_data(grd%id_unused, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  unused_calving=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_clock_end(bergs%clock_int)

  call mpp_clock_begin(bergs%clock_cal)
  ! Calve excess stored ice into icebergs
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, before calving()   ')
  call calve_icebergs(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (calved)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after calving')
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after calving()    ')
  call mpp_clock_end(bergs%clock_cal)

  !need to redo MTS transfers after a restart
  if (.not. Visited) then
    Visited=.true.
    if (bergs%mts) then
      call interp_gridded_fields_to_bergs(bergs)
      call transfer_mts_bergs(bergs)
    elseif ((bergs%contact_distance>0.) .or. (bergs%contact_spring_coef .ne. bergs%spring_coef)) then
      call set_conglom_ids(bergs)
    endif
    if (bergs%iceberg_bonds_on) call orig_bond_length(bergs)
    if (monitor_energy) call energy_tests_init(bergs)
    if (mpp_pe()==0) write(*,'(a)') 'KID, iceberg_run: completed first visit initialization'
  endif

  if (.not. bergs%mts) call interp_gridded_fields_to_bergs(bergs)

  ! For each berg, evolve
  call mpp_clock_begin(bergs%clock_mom)

  if (.not.bergs%Static_icebergs) then
    if (bergs%iceberg_bonds_on .and. (.not. bergs%dem)) call reset_bond_rotation(bergs)

    if (bergs%mts) then
      call evolve_icebergs_mts(bergs)
    else
      call evolve_icebergs(bergs) !The single time step (STS) scheme
    end if
    if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after evolve()     ')
  endif
  call move_berg_between_cells(bergs)  !Markpoint6
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after move_lists() ')
  if (debug) call bergs_chksum(bergs, 'run bergs (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(bergs%grd, 's/r run after evolve')
  call mpp_clock_end(bergs%clock_mom)

  ! Send bergs to other PEs
  call mpp_clock_begin(bergs%clock_com1)
  if (bergs%iceberg_bonds_on)  call  bond_address_update(bergs)

  call send_bergs_to_other_pes(bergs)
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after send_bergs() ')
  call mpp_clock_end(bergs%clock_com1)

  call mpp_clock_begin(bergs%clock_fl1)
  ! Footloose mechanism part 1: calve the child icebergs
  if (bergs%fl_r>0.) call footloose_calving(bergs, time)
  call mpp_clock_end(bergs%clock_fl1)

  call mpp_clock_begin(bergs%clock_com2)
  if (bergs%mts) then
    call interp_gridded_fields_to_bergs(bergs)
    call transfer_mts_bergs(bergs)
  else
    if ((bergs%interactive_icebergs_on) .or. (bergs%iceberg_bonds_on)) then
      call update_halo_icebergs(bergs)
      if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after update_halo()')
      if (bergs%iceberg_bonds_on) then
        call connect_all_bonds(bergs)
      else
        if (grd%Lx>0.) call update_latlon(bergs)
      endif
      if ((bergs%contact_distance>0.) .or. (bergs%contact_spring_coef .ne. bergs%spring_coef)) then
        call set_conglom_ids(bergs)
      endif
    endif
    call interp_gridded_fields_to_bergs(bergs)
  endif
  if (debug) call bergs_chksum(bergs, 'run bergs (exchanged)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after exchange')
  call mpp_clock_end(bergs%clock_com2)

  call mpp_clock_begin(bergs%clock_fl2)
  ! Footloose mechanism part 2:
  if (bergs%fl_r>0.) then
    !delete any edge elements that have fully calved from the footloose mechanism
    if (bergs%iceberg_bonds_on) call delete_fully_fl_calved_edge_elements(bergs)
    !child bergs become interactive once they are out of contact range of any other
    !berg for the first time
    if (bergs%interactive_icebergs_on) call adjust_fl_berg_interactivity(bergs)
  endif
  call mpp_clock_end(bergs%clock_fl2)

  if (bergs%find_melt_using_spread_mass .or. bergs%mts) then
    call calculate_mass_on_ocean(bergs, with_diagnostics=.false.)
  end if

  ! Calculate mass on ocean before thermodynamics, to use in melt rate calculation
  if (bergs%find_melt_using_spread_mass) then
    grd%spread_mass_old(:,:)=0.
    call sum_up_spread_fields(bergs,grd%spread_mass_old(grd%isc:grd%iec,grd%jsc:grd%jec), 'mass')
    !Reset fields
    grd%mass_on_ocean(:,:,:)=0.  ;  grd%area_on_ocean(:,:,:)=0.
    grd%Uvel_on_ocean(:,:,:)=0.  ;  grd%Vvel_on_ocean(:,:,:)=0.
  endif

  ! Iceberg thermodynamics (melting) + rolling
  call mpp_clock_begin(bergs%clock_the)
  call thermodynamics(bergs)
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after thermodyn()  ')
  if (debug) call bergs_chksum(bergs, 'run bergs (thermo)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after thermodynamics')
  call mpp_clock_end(bergs%clock_the)

  !Creating gridded fields from new icebergs
  call create_gridded_icebergs_fields(bergs)

  ! For each berg, record
  call mpp_clock_begin(bergs%clock_dia)
  if (sample_traj .or. bergs%writeandstop) call record_posn(bergs)
  if (write_traj .or. bergs%writeandstop) then
    call move_all_trajectories(bergs)
    call write_trajectory(bergs%trajectories, bergs%save_short_traj, bergs%save_fl_traj, bergs%fl_r)
    if (save_bond_traj) call write_bond_trajectory(bergs%bond_trajectories)
  endif

  if (bergs%writeandstop) then
    call mpp_sync()
    call error_mesg('KID', 'WRITE AND STOP!!!', FATAL)
  endif

  ! Gridded diagnostics
  if (grd%id_uo>0) &
    lerr=send_data(grd%id_uo, grd%uo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_vo>0) &
    lerr=send_data(grd%id_vo, grd%vo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_ui>0) &
    lerr=send_data(grd%id_ui, grd%ui(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_vi>0) &
    lerr=send_data(grd%id_vi, grd%vi(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_ua>0) &
    lerr=send_data(grd%id_ua, grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_va>0) &
    lerr=send_data(grd%id_va, grd%va(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_sst>0) &
    lerr=send_data(grd%id_sst, grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_sss>0) &
    lerr=send_data(grd%id_sss, grd%sss(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_cn>0) &
    lerr=send_data(grd%id_cn, grd%cn(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_hi>0) &
    lerr=send_data(grd%id_hi, grd%hi(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_floating_melt>0) &
    lerr=send_data(grd%id_floating_melt, grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_m_per_year>0) &
    lerr=send_data(grd%id_melt_m_per_year, grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec)* (86400.0*365.0/bergs%rho_bergs), Time)
  if (grd%id_berg_melt>0) &
    lerr=send_data(grd%id_berg_melt, grd%berg_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_buoy>0) &
    lerr=send_data(grd%id_melt_buoy, grd%melt_buoy(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_eros>0) &
    lerr=send_data(grd%id_melt_eros, grd%melt_eros(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_conv>0) &
    lerr=send_data(grd%id_melt_conv, grd%melt_conv(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_virtual_area>0) &
    lerr=send_data(grd%id_virtual_area, grd%virtual_area(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_bergy_src>0) &
    lerr=send_data(grd%id_bergy_src, grd%bergy_src(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_bergy_melt>0) &
    lerr=send_data(grd%id_bergy_melt, grd%bergy_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_bergy_mass>0) &
    lerr=send_data(grd%id_bergy_mass, grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fl_bits_src>0) &
    lerr=send_data(grd%id_fl_bits_src, grd%fl_bits_src(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fl_bits_melt>0) &
    lerr=send_data(grd%id_fl_bits_melt, grd%fl_bits_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fl_bits_mass>0) &
    lerr=send_data(grd%id_fl_bits_mass, grd%fl_bits_mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fl_bergy_bits_mass>0) &
    lerr=send_data(grd%id_fl_bergy_bits_mass, grd%fl_bergy_bits_mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_spread_mass>0) &
    lerr=send_data(grd%id_spread_mass, grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_spread_area>0) &
    lerr=send_data(grd%id_spread_area, grd%spread_area(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_u_iceberg>0) &
    lerr=send_data(grd%id_u_iceberg, grd%u_iceberg(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_v_iceberg>0) &
    lerr=send_data(grd%id_v_iceberg, grd%v_iceberg(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_spread_uvel>0) &
    lerr=send_data(grd%id_spread_uvel, grd%spread_uvel(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_spread_vvel>0) &
    lerr=send_data(grd%id_spread_vvel, grd%spread_vvel(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_ustar_iceberg>0) &
    lerr=send_data(grd%id_ustar_iceberg, grd%ustar_iceberg(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_mass>0) &
    lerr=send_data(grd%id_mass, grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_stored_ice>0) &
    lerr=send_data(grd%id_stored_ice, grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)
  if (grd%id_rmean_calving>0) &
    lerr=send_data(grd%id_rmean_calving, grd%rmean_calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_rmean_calving_hflx>0) &
    lerr=send_data(grd%id_rmean_calving_hflx, grd%rmean_calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_real_calving>0) &
    lerr=send_data(grd%id_real_calving, grd%real_calving(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)
  if (grd%id_ssh>0) &
    lerr=send_data(grd%id_ssh, grd%ssh(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fax>0) &
    lerr=send_data(grd%id_fax, tauxa(:,:), Time)
  if (grd%id_fay>0) &
    lerr=send_data(grd%id_fay, tauya(:,:), Time)
  if (grd%id_melt_by_class>0) &
    lerr=send_data(grd%id_melt_by_class, grd%melt_by_class(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)
  if (grd%id_melt_buoy_fl>0) &
    lerr=send_data(grd%id_melt_buoy_fl, grd%melt_buoy_fl(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_eros_fl>0) &
    lerr=send_data(grd%id_melt_eros_fl, grd%melt_eros_fl(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_conv_fl>0) &
    lerr=send_data(grd%id_melt_conv_fl, grd%melt_conv_fl(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fl_parent_melt>0) &
    lerr=send_data(grd%id_fl_parent_melt, grd%fl_parent_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fl_child_melt>0) &
    lerr=send_data(grd%id_fl_child_melt, grd%fl_child_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_count>0) then
    allocate( iCount(grd%isc:grd%iec,grd%jsc:grd%jec) ); iCount(:,:)=0
    do j = grd%jsc, grd%jec ; do i = grd%isc, grd%iec
      iCount(i,j) = count_bergs_in_list(bergs%list(i,j)%first)
    enddo ; enddo
    lerr=send_data(grd%id_count, iCount(:,:), Time)
    deallocate( iCount )
  endif
  if (grd%id_chksum>0) then
    allocate( iCount(grd%isc:grd%iec,grd%jsc:grd%jec) ); iCount(:,:)=0
    do j = grd%jsc, grd%jec ; do i = grd%isc, grd%iec
      iCount(i,j) = list_chksum(bergs%list(i,j)%first)
    enddo ; enddo
    lerr=send_data(grd%id_chksum, iCount(:,:), Time)
    deallocate( iCount )
  endif

  ! Dump icebergs to screen
  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_run, status')

  ! Dump icebergs bonds to screen
  if (really_debug)  call show_all_bonds(bergs)

  call mpp_clock_end(bergs%clock_dia)

  ! This is the point in the algorithm which determines which fields get passed to the ice model
  ! Return what ever calving we did not use and additional icebergs melt

  ! Making sure that spread_mass has the correct mass
  !grd%spread_mass(:,:)=0.0
  !call icebergs_incr_mass(bergs, grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec), within_iceberg_model=.True.)


  ! Return what ever calving we did not use and additional icebergs melt
  call mpp_clock_begin(bergs%clock_int)
  if (.not. bergs%passive_mode) then
    where (grd%area(grd%isc:grd%iec,grd%jsc:grd%jec)>0.)
      calving(:,:)=grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)/grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) &
                  +grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec)
    elsewhere
      calving(:,:)=0.
    end where
    calving_hflx(:,:)=grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)
    !Return iceberg mass, area and ustar to pass on to ocean model
    if (present(mass_berg)) then
      if (associated(mass_berg)) then
        if (bergs%add_weight_to_ocean) &
          mass_berg(:,:)=grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec)
      endif
    endif
    if (present(ustar_berg)) then
      if (associated(ustar_berg)) then
        ustar_berg(:,:)=grd%ustar_iceberg(grd%isc:grd%iec,grd%jsc:grd%jec)
      endif
    endif
    if (present(area_berg)) then
      if (associated(area_berg)) then
        area_berg(:,:)=grd%spread_area(grd%isc:grd%iec,grd%jsc:grd%jec)
      endif
    endif
  endif

  call mpp_clock_end(bergs%clock_int)

  ! Diagnose budgets
  call mpp_clock_begin(bergs%clock_dia)
  tmpsum=sum( grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_melt=bergs%net_melt+tmpsum*bergs%dt
  bergs%net_outgoing_calving=bergs%net_outgoing_calving+(unused_calving+tmpsum)*bergs%dt
  tmpsum=sum( grd%berg_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%berg_melt=bergs%berg_melt+tmpsum*bergs%dt
  tmpsum=sum( grd%bergy_src(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%bergy_src=bergs%bergy_src+tmpsum*bergs%dt
  tmpsum=sum( grd%bergy_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%bergy_melt=bergs%bergy_melt+tmpsum*bergs%dt
  tmpsum=sum( grd%fl_bits_src(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%fl_bits_src=bergs%fl_bits_src+tmpsum*bergs%dt
  tmpsum=sum( grd%fl_bits_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%fl_bits_melt=bergs%fl_bits_melt+tmpsum*bergs%dt
  tmpsum=sum( calving(:,:)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_calving_returned=bergs%net_calving_returned+tmpsum*bergs%dt
  tmpsum=sum( grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_outgoing_calving_heat=bergs%net_outgoing_calving_heat+tmpsum*bergs%dt ! Units of J
  if (lbudget) then
    bergs%stored_end=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    bergs%rmean_calving_end=sum( grd%rmean_calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%rmean_calving_hflx_end=sum( grd%rmean_calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%stored_heat_end=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%floating_mass_end=sum_mass(bergs)
    bergs%icebergs_mass_end=sum_mass(bergs,justbergs=.true.)
    bergs%bergy_mass_end=sum_mass(bergs,justbits=.true.)
    bergs%fl_bits_mass_end=sum_mass(bergs,justflbits=.true.)
    !bergs%spread_mass_end=sum_mass(bergs) !Not sure what this is
    !bergs%spread_area_end=sum_mass(bergs) !Not sure what this is
    !bergs%u_iceberg_end=sum_mass(bergs) !Not sure what this is
    !bergs%v_iceberg_end=sum_mass(bergs) !Not sure what this is
    bergs%floating_heat_end=sum_heat(bergs)
    grd%tmpc(:,:)=0.
    !Finding spread mass
    call mpp_clock_end(bergs%clock); call mpp_clock_end(bergs%clock_dia) ! To enable calling of public s/r
    call sum_up_spread_fields(bergs, grd%tmpc, 'mass')
    call mpp_clock_begin(bergs%clock_dia); call mpp_clock_begin(bergs%clock) ! To enable calling of public s/r
    bergs%returned_mass_on_ocean=sum( grd%tmpc(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    !Finding spread area
    call mpp_clock_end(bergs%clock); call mpp_clock_end(bergs%clock_dia) ! To enable calling of public s/r
    call sum_up_spread_fields(bergs, grd%tmpc, 'area')
    call mpp_clock_begin(bergs%clock_dia); call mpp_clock_begin(bergs%clock) ! To enable calling of public s/r
    bergs%returned_area_on_ocean=sum( grd%tmpc(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%nbergs_end=count_bergs(bergs)
    call mpp_sum(bergs%stored_end)
    call mpp_sum(bergs%stored_heat_end)
    call mpp_sum(bergs%floating_mass_end)
    call mpp_sum(bergs%icebergs_mass_end)
    call mpp_sum(bergs%bergy_mass_end)
    call mpp_sum(bergs%fl_bits_mass_end)
    call mpp_sum(bergs%spread_mass_end)
    call mpp_sum(bergs%spread_area_end)
    call mpp_sum(bergs%u_iceberg_end)
    call mpp_sum(bergs%v_iceberg_end)
    call mpp_sum(bergs%spread_uvel_end)
    call mpp_sum(bergs%spread_vvel_end)
    call mpp_sum(bergs%ustar_iceberg_end)
    call mpp_sum(bergs%floating_heat_end)
    call mpp_sum(bergs%returned_mass_on_ocean)
    call mpp_sum(bergs%nbergs_end)
    call mpp_sum(bergs%nbergs_calved)
    call mpp_sum(bergs%nbergs_calved_fl)
    do k=1,nclasses; call mpp_sum(bergs%nbergs_calved_by_class_s(k)); enddo
    do k=1,nclasses; call mpp_sum(bergs%nbergs_calved_by_class_n(k)); enddo
    call mpp_sum(bergs%nbergs_melted)
    call mpp_sum(bergs%nspeeding_tickets)
    call mpp_sum(bergs%net_calving_returned)
    call mpp_sum(bergs%net_outgoing_calving)
    call mpp_sum(bergs%net_calving_received)
    call mpp_sum(bergs%net_incoming_calving)
    call mpp_sum(bergs%net_incoming_calving_heat)
    call mpp_sum(bergs%net_incoming_calving_heat_used)
    call mpp_sum(bergs%net_outgoing_calving_heat)
    call mpp_sum(bergs%net_calving_used)
    call mpp_sum(bergs%net_calving_to_bergs)
    call mpp_sum(bergs%net_heat_to_bergs)
    call mpp_sum(bergs%net_heat_to_ocean)
    call mpp_sum(bergs%net_melt)
    call mpp_sum(bergs%berg_melt)
    call mpp_sum(bergs%bergy_src)
    call mpp_sum(bergs%bergy_melt)
    call mpp_sum(bergs%fl_bits_src)
    call mpp_sum(bergs%fl_bits_melt)
    grdd_berg_mass=sum( grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_berg_mass)
    grdd_bergy_mass=sum( grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_bergy_mass)
    grdd_fl_bits_mass=sum( grd%fl_bits_mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_fl_bits_mass)
    grdd_spread_mass=sum( grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_spread_mass)
    grdd_spread_area=sum( grd%spread_area(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_spread_area)
    if (mpp_pe().eq.mpp_root_pe()) then
 100 format("KID: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
 200 format("KID: ",a19,10(a18,"=",es14.7,x,a2,:,","))
      call report_state('stored ice','kg','',bergs%stored_start,'',bergs%stored_end,'')
      call report_state('floating','kg','',bergs%floating_mass_start,'',bergs%floating_mass_end,'',bergs%nbergs_end)
      call report_state('icebergs','kg','',bergs%icebergs_mass_start,'',bergs%icebergs_mass_end,'')
      call report_state('bits','kg','',bergs%bergy_mass_start,'',bergs%bergy_mass_end,'')
      call report_state('fl_bits','kg','',bergs%fl_bits_mass_start,'',bergs%fl_bits_mass_end,'')
      call report_state('spread icebergs','kg','',bergs%spread_mass_start,'',bergs%spread_mass_end,'')
      call report_state('spread icebergs','m^2','',bergs%spread_area_start,'',bergs%spread_area_end,'')
      call report_istate('berg #','',bergs%nbergs_start,'',bergs%nbergs_end,'')
      call report_ibudget('berg #','calved',bergs%nbergs_calved, &
                                   'FL calved',bergs%nbergs_calved_fl, &
                                   'melted',bergs%nbergs_melted, &
                                   '#',bergs%nbergs_start,bergs%nbergs_end)
      call report_budget('stored mass','kg','calving used',bergs%net_calving_used, &
                                            'bergs',bergs%net_calving_to_bergs, &
                                            'stored mass',bergs%stored_start,bergs%stored_end)
      call report_budget('floating mass','kg','calving used',bergs%net_calving_to_bergs, &
                                              'bergs',bergs%net_melt, &
                                              'stored mass',bergs%floating_mass_start,bergs%floating_mass_end)
      call report_budget('berg mass','kg','calving',bergs%net_calving_to_bergs, &
                                          'melt+eros+fl',bergs%berg_melt+bergs%fl_bits_src, &
                                          'berg mass',bergs%icebergs_mass_start,bergs%icebergs_mass_end)
      call report_budget('bits mass','kg','eros used',bergs%bergy_src, &
                                          'bergs',bergs%bergy_melt, &
                                          'stored mass',bergs%bergy_mass_start,bergs%bergy_mass_end)
      call report_budget('fl bits mass','kg','fl calving',bergs%fl_bits_src, &
                                          'fl melt+eros',bergs%fl_bits_melt, &
                                          'fl mass',bergs%fl_bits_mass_start,bergs%fl_bits_mass_end)
      call report_budget('net mass','kg','recvd',bergs%net_calving_received, &
                                         'rtrnd',bergs%net_calving_returned, &
                                         'net mass',bergs%stored_start+bergs%floating_mass_start, &
                                                    bergs%stored_end+bergs%floating_mass_end)
      call report_consistant('iceberg mass','kg','gridded',grdd_berg_mass,'bergs',bergs%icebergs_mass_end)
      call report_consistant('spread mass','kg','gridded',grdd_spread_mass,'bergs',bergs%spread_mass_end)
      call report_consistant('spread area','kg','gridded',grdd_spread_area,'bergs',bergs%spread_area_end)
      call report_consistant('bits mass','kg','gridded',grdd_bergy_mass,'bits',bergs%bergy_mass_end)
      call report_consistant('fl bits mass','kg','gridded',grdd_fl_bits_mass,'fl bits',bergs%fl_bits_mass_end)
      call report_consistant('wieght','kg','returned',bergs%returned_mass_on_ocean,'floating',bergs%floating_mass_end)
      call report_state('net heat','J','',bergs%stored_heat_start+bergs%floating_heat_start,'',&
           & bergs%stored_heat_end+bergs%floating_heat_end,'')
      call report_state('stored heat','J','',bergs%stored_heat_start,'',bergs%stored_heat_end,'')
      call report_state('floating heat','J','',bergs%floating_heat_start,'',bergs%floating_heat_end,'')
      call report_budget('net heat','J','net heat',bergs%net_incoming_calving_heat, &
                                        'net heat',bergs%net_outgoing_calving_heat, &
                                        'net heat',bergs%stored_heat_start+bergs%floating_heat_start, &
                                                   bergs%stored_heat_end+bergs%floating_heat_end)
      call report_budget('stored heat','J','calving used',bergs%net_incoming_calving_heat_used, &
                                           'bergs',bergs%net_heat_to_bergs, &
                                           'net heat',bergs%stored_heat_start,bergs%stored_heat_end)
      call report_budget('flting heat','J','calved',bergs%net_heat_to_bergs, &
                                           'melt',bergs%net_heat_to_ocean, &
                                           'net heat',bergs%floating_heat_start,bergs%floating_heat_end)
      if (debug) then
        call report_consistant('top interface','kg','from SIS',bergs%net_incoming_calving,'seen by KID',&
             & bergs%net_calving_received)
        call report_consistant('bot interface','kg','sent',bergs%net_outgoing_calving,'seen by SIS',bergs%net_calving_returned)
      endif
      write(*,'("KID: calved by class S. hemisphere = ",i4,20(",",i4))') (bergs%nbergs_calved_by_class_s(k),k=1,nclasses)
      write(*,'("KID: calved by class N. hemisphere = ",i4,20(",",i4))') (bergs%nbergs_calved_by_class_n(k),k=1,nclasses)
      if (bergs%nspeeding_tickets>0) write(*,'("KID: speeding tickets issued = ",i4)') bergs%nspeeding_tickets
    endif
    bergs%nbergs_start=bergs%nbergs_end
    bergs%stored_start=bergs%stored_end
    bergs%nbergs_melted=0
    bergs%nbergs_calved=0
    bergs%nbergs_calved_fl=0
    bergs%nbergs_calved_by_class_s(:)=0
    bergs%nbergs_calved_by_class_n(:)=0
    bergs%nspeeding_tickets=0
    bergs%stored_heat_start=bergs%stored_heat_end
    bergs%floating_heat_start=bergs%floating_heat_end
    bergs%floating_mass_start=bergs%floating_mass_end
    bergs%icebergs_mass_start=bergs%icebergs_mass_end
    bergs%bergy_mass_start=bergs%bergy_mass_end
    bergs%fl_bits_mass_start=bergs%fl_bits_mass_end
    bergs%spread_mass_start=bergs%spread_mass_end
    bergs%spread_area_start=bergs%spread_area_end
    bergs%net_calving_used=0.
    bergs%net_calving_to_bergs=0.
    bergs%net_heat_to_bergs=0.
    bergs%net_heat_to_ocean=0.
    bergs%net_calving_received=0.
    bergs%net_calving_returned=0.
    bergs%net_incoming_calving=0.
    bergs%net_outgoing_calving=0.
    bergs%net_incoming_calving_heat=0.
    bergs%net_incoming_calving_heat_used=0.
    bergs%net_outgoing_calving_heat=0.
    bergs%net_melt=0.
    bergs%berg_melt=0.
    bergs%bergy_melt=0.
    bergs%bergy_src=0.
    bergs%fl_bits_melt=0.
    bergs%fl_bits_src=0.

    if (bergs%iceberg_bonds_on) then
      check_bond_quality=.true.
      nbonds=0
      call count_bonds(bergs, nbonds,check_bond_quality)
    endif
  endif

  if (debug) call bergs_chksum(bergs, 'run bergs (bot)')
  if (debug) call checksum_gridded(bergs%grd, 'end of s/r run')
  call mpp_clock_end(bergs%clock_dia)

  call mpp_clock_end(bergs%clock)

end subroutine icebergs_run

!> Prints summary of start and end states
subroutine report_state(budgetstr, budgetunits, startstr, startval, endstr, endval, delstr, nbergs)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: budgetunits !< Units of budgeted quantity
  character*(*), intent(in) :: startstr !< Start label
  real, intent(in) :: startval !< Start value for budget
  character*(*), intent(in) :: endstr !< End label
  real, intent(in) :: endval !< End value for budget
  character*(*), intent(in) :: delstr !< Delta label
  integer, intent(in), optional :: nbergs !< Number of bergs
  ! Local variables
  if (present(nbergs)) then
    write(*,100) budgetstr//' state:', &
                        startstr//' start',startval,budgetunits, &
                        endstr//' end',endval,budgetunits, &
                        'Delta '//delstr,endval-startval,budgetunits, &
                        '# of bergs',nbergs
  else
    write(*,100) budgetstr//' state:', &
                        startstr//' start',startval,budgetunits, &
                        endstr//' end',endval,budgetunits, &
                        delstr//'Delta',endval-startval,budgetunits
  endif
  100 format("KID: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
end subroutine report_state

!> Prints consistency summary of start and end states
subroutine report_consistant(budgetstr, budgetunits, startstr, startval, endstr, endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: budgetunits !< Units of budgeted quantity
  character*(*), intent(in) :: startstr !< Start label
  real, intent(in) :: startval !< Start value for budget
  character*(*), intent(in) :: endstr !< End label
  real, intent(in) :: endval !< End value for budget
  ! Local variables
  write(*,200) budgetstr//' check:', &
                      startstr,startval,budgetunits, &
                      endstr,endval,budgetunits, &
                      'error',(endval-startval)/((endval+startval)+1e-30),'nd'
  200 format("KID: ",a19,10(a18,"=",es14.7,x,a2,:,","))
end subroutine report_consistant

!> Prints a budget
subroutine report_budget(budgetstr, budgetunits, instr, inval, outstr, outval, delstr, startval, endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: budgetunits !< Units of budgeted quantity
  character*(*), intent(in) :: instr !< Incoming label
  real, intent(in) :: inval !< Incoming value
  character*(*), intent(in) :: outstr !< Outgoing label
  real, intent(in) :: outval !< Outgoing value
  character*(*), intent(in) :: delstr !< Delta label
  real, intent(in) :: startval !< Start value for budget
  real, intent(in) :: endval !< End value for budget
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval,budgetunits, &
                      outstr//' out',outval,budgetunits, &
                      'Delta '//delstr,inval-outval,budgetunits, &
                      'error',((endval-startval)-(inval-outval))/max(1.e-30,max(abs(endval-startval),abs(inval-outval))),'nd'
  200 format("KID: ",a19,3(a18,"=",es14.7,x,a2,:,","),a8,"=",es10.3,x,a2)
end subroutine report_budget

!> Prints summary of start and end states
subroutine report_istate(budgetstr, startstr, startval, endstr, endval, delstr)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: startstr !< Start label
  integer, intent(in) :: startval !< Start value for budget
  character*(*), intent(in) :: endstr !< End label
  integer, intent(in) :: endval !< End value for budget
  character*(*), intent(in) :: delstr !< Delta label
  ! Local variables
  write(*,100) budgetstr//' state:', &
                        startstr//' start',startval, &
                        endstr//' end',endval, &
                        delstr//'Delta',endval-startval
  100 format("KID: ",a19,3(a18,"=",i14,x,:,","))
end subroutine report_istate

!> Prints a budget
subroutine report_ibudget(budgetstr,instr1,inval1,instr2,inval2,outstr,outval,delstr,startval,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: instr1 !< Incoming label 1
  integer, intent(in) :: inval1 !< Incoming value 1
  character*(*), intent(in) :: instr2 !< Incoming label 2
  integer, intent(in) :: inval2 !< Incoming value 2
  character*(*), intent(in) :: outstr !< Outgoing label
  integer, intent(in) :: outval !< Outgoing value
  character*(*), intent(in) :: delstr !< Delta label
  integer, intent(in) :: startval !< Start value for budget
  integer, intent(in) :: endval !< End value for budget
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr1//' in',inval1, &
                      instr2//' in',inval2, &
                      outstr//' out',outval, &
                      'Delta '//delstr,inval1+inval2-outval, &
                      'error',((endval-startval)-(inval1+inval2-outval))
  200 format("KID: ",a19,10(a18,"=",i14,x,:,","))
end subroutine report_ibudget

!> Time-filter calving and calving_hflx with a running mean.
!!
!! This subroutine takes in the new calving and calving_hflx, and uses them to time step a running-mean_calving value.
!! The time stepping uses a time scale tau. When tau is equal to zero, the
!! running mean is exactly equal to the new calving value.
subroutine get_running_mean_calving(bergs, calving, calving_hflx)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  real, dimension(:,:), intent(inout) :: calving !< Calving (kg/s)
  real, dimension(:,:), intent(inout) :: calving_hflx !< Calving heat flux (W/m2)
  ! Local variables
  real :: alpha  !Parameter used for calving relaxation time stepping.  (0<=alpha<1)
  real :: tau  !Relaxation timescale in seconds
  real :: beta  ! = 1-alpha (0<=beta<1)

  ! For the first time-step, initialize the running mean with the current data
  if (.not. bergs%grd%rmean_calving_initialized) then
    bergs%grd%rmean_calving(:,:)=calving(:,:)
    bergs%grd%rmean_calving_initialized=.true.
  endif
  if (.not. bergs%grd%rmean_calving_hflx_initialized) then
    bergs%grd%rmean_calving_hflx(:,:)=calving_hflx(:,:)
    bergs%grd%rmean_calving_hflx_initialized=.true.
  endif

  !Applying "Newton cooling" with timescale tau, to smooth out the calving field.
  tau=bergs%tau_calving/(365.*24*60*60) !Converting time scale from years to seconds
  alpha=tau/(tau+bergs%dt)
  if (alpha==0.) return ! Avoids unnecessary copying of arrays
  if (alpha>0.5) then ! beta is small
    beta=bergs%dt/(tau+bergs%dt)
    alpha=1.-beta
  else ! alpha is small
    beta=1.-alpha
  endif

  ! For non-negative alpha and beta, these expressions for the running means are sign preserving
  bergs%grd%rmean_calving(:,:)=beta*calving(:,:) + alpha*bergs%grd%rmean_calving(:,:)
  bergs%grd%rmean_calving_hflx(:,:)=beta*calving_hflx(:,:) + alpha*bergs%grd%rmean_calving_hflx(:,:)

  !Setting calving used by the iceberg model equal to the running mean
  calving(:,:)=bergs%grd%rmean_calving(:,:)
  calving_hflx(:,:)=bergs%grd%rmean_calving_hflx(:,:)

end subroutine get_running_mean_calving

!> Increments a gridded mass field with the mass of bergs (called from outside icebergs_run)
!!
!! This routine is called from SIS, (and older versions of SIS2), but not within
!! the iceberg model. The routine adds the spread iceberg mass to mass provided
!! the add weight to ocean flag is on, and passive mode is off. It also appears to
!! play some role in diagnostics
subroutine icebergs_incr_mass(bergs, mass, Time)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  real, dimension(bergs%grd%isc:bergs%grd%iec,bergs%grd%jsc:bergs%grd%jec), intent(inout) :: mass !< Mass field to increment
  type(time_type), intent(in), optional :: Time !< Model time
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  integer :: i, j
  logical :: lerr

  if (.not. associated(bergs)) return
  if (.not. bergs%add_weight_to_ocean) return

  ! For convenience
  grd=>bergs%grd

  !Start the clocks
  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_int)

  do j=grd%jsc, grd%jec; do i=grd%isc, grd%iec
    if (.not. bergs%passive_mode)   mass(i,j)=mass(i,j) + grd%spread_mass(i,j)
  enddo  ;enddo

  !Stop the clocks
  call mpp_clock_end(bergs%clock_int)
  call mpp_clock_end(bergs%clock)

end subroutine icebergs_incr_mass

!> Sums up a berg property into a gridded field
subroutine sum_up_spread_fields(bergs, field, field_name)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  real, dimension(bergs%grd%isc:bergs%grd%iec,bergs%grd%jsc:bergs%grd%jec), intent(out) :: field !< Gridded field
  character(len=4), intent(in) :: field_name !< Name of field to grid
  ! Local variables
  integer :: i, j
  type(icebergs_gridded), pointer :: grd
  real :: dmda
  logical :: lerr
  real, dimension(bergs%grd%isd:bergs%grd%ied, bergs%grd%jsd:bergs%grd%jed,9) :: var_on_ocean   !Variable being spread onto the ocean  (mass, area, Uvel, Vvel)
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()
  ! For convenience
  grd=>bergs%grd

  field(:,:)=0.

  !Deciding which varibale to spread across cells across grid cells
  if (field_name=='mass') var_on_ocean(:,:,:)=grd%mass_on_ocean(:,:,:)
  if (field_name=='area') var_on_ocean(:,:,:)=grd%area_on_ocean(:,:,:)
  if (field_name=='Uvel') var_on_ocean(:,:,:)=grd%Uvel_on_ocean(:,:,:)
  if (field_name=='Vvel') var_on_ocean(:,:,:)=grd%Vvel_on_ocean(:,:,:)

  !This line has been removed, for that routine can be used for other fields
  !if (.not. bergs%add_weight_to_ocean) return

  !Update the halos of the var_on_ocean
  call mpp_update_domains(var_on_ocean, grd%domain)

  !Rotatine when old_bug_rotated_weights is on  - we should remove this.
  if (.not. old_bug_rotated_weights) then
    do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
      if (grd%parity_x(i,j)<0.) then
        ! This block assumes both parity_x and parity_y are negative
        ! (i.e. a 180 degree rotation). In general, we should handle
        ! +/- 90 degree rotations as well but in CM2*-class models
        ! this is not necessary. -aja
        dmda=var_on_ocean(i,j,9); var_on_ocean(i,j,9)=var_on_ocean(i,j,1); var_on_ocean(i,j,1)=dmda
        dmda=var_on_ocean(i,j,8); var_on_ocean(i,j,8)=var_on_ocean(i,j,2); var_on_ocean(i,j,2)=dmda
        dmda=var_on_ocean(i,j,7); var_on_ocean(i,j,7)=var_on_ocean(i,j,3); var_on_ocean(i,j,3)=dmda
        dmda=var_on_ocean(i,j,6); var_on_ocean(i,j,6)=var_on_ocean(i,j,4); var_on_ocean(i,j,4)=dmda
      endif
    enddo; enddo
  endif

  !Here we add the contribution of the 9 cells. This is the heart of the routine.
  do j=grd%jsc, grd%jec; do i=grd%isc, grd%iec
    dmda=var_on_ocean(i,j,5) &
         + ( ( (var_on_ocean(i-1,j-1,9)+var_on_ocean(i+1,j+1,1))   &
         +     (var_on_ocean(i+1,j-1,7)+var_on_ocean(i-1,j+1,3)) ) &
         +   ( (var_on_ocean(i-1,j  ,6)+var_on_ocean(i+1,j  ,4))   &
         +     (var_on_ocean(i  ,j-1,8)+var_on_ocean(i  ,j+1,2)) ) )
    if (grd%area(i,j)>0) dmda=dmda/grd%area(i,j)*grd%msk(i,j)

    !Make sure that area <=1.0
    if (field_name=='area') dmda=min(dmda,1.0)

    field(i,j)=dmda
  enddo; enddo

  if (debug) then
    grd%tmp(:,:)=0.; grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=field
    if (field_name=='mass') then
      call grd_chksum3(grd, grd%mass_on_ocean, 'mass bergs (incr)')
      call grd_chksum2(grd, grd%tmp, 'mass out (incr)')
    elseif (field_name=='area') then
      call grd_chksum3(grd, grd%area_on_ocean, 'area bergs (incr)')
      call grd_chksum2(grd, grd%tmp, 'area out (incr)')
    endif
 endif
end subroutine sum_up_spread_fields

!> Adds calving (from driver) to the coastal "buckets"
subroutine accumulate_calving(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: remaining_dist_s, remaining_dist_n, net_calving_used
  integer :: k, i, j
  logical, save :: first_call=.true.
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  ! This is a hack to simplify initialization
  if (first_call.and..not.bergs%restarted) then
    first_call=.false.
   !do k=1, nclasses
   !  where (grd%calving==0.) grd%stored_ice(:,:,k)=0.
   !enddo
    bergs%stored_start=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    call mpp_sum( bergs%stored_start )
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,es13.6,a)') &
        'KID, accumulate_calving: initial stored mass=',bergs%stored_start,' kg'
    do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (grd%calving(i,j).ne.0.) grd%stored_heat(i,j)= & ! Need units of J
            sum(grd%stored_ice(i,j,:)) & ! initial stored ice in kg
           *grd%calving_hflx(i,j)*grd%area(i,j) & ! J/s/m2 x m^2 = J/s
           /grd%calving(i,j) ! /calving in kg/s
    enddo; enddo
    bergs%stored_heat_start=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum( bergs%stored_heat_start )
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,es13.6,a)') &
        'KID, accumulate_calving: initial stored heat=',bergs%stored_heat_start,' J'
   endif

  remaining_dist_s=1.; remaining_dist_n=1.
  do k=1, nclasses
    where (grd%lat<0.)
      grd%stored_ice(:,:,k)=grd%stored_ice(:,:,k)+bergs%dt*grd%calving(:,:)*bergs%distribution_s(k)
    elsewhere
      grd%stored_ice(:,:,k)=grd%stored_ice(:,:,k)+bergs%dt*grd%calving(:,:)*bergs%distribution_n(k)
    end where
    remaining_dist_s=remaining_dist_s-bergs%distribution_s(k)
    remaining_dist_n=remaining_dist_n-bergs%distribution_n(k)
  enddo
  if (remaining_dist_s.lt.0. .or. remaining_dist_n.lt.0.) then
    write(stderrunit,*) 'KID, accumulate_calving: sum(distribution)>1!!!',remaining_dist_s, remaining_dist_n
    call error_mesg('KID, accumulate_calving', 'calving is OVER distributed!', WARNING)
  endif
  where (grd%lat<0.)
    grd%tmp=remaining_dist_s
  elsewhere
    grd%tmp=remaining_dist_n
  end where
  net_calving_used=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) *(1.-grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)))
  bergs%net_calving_used=bergs%net_calving_used+net_calving_used*bergs%dt
  ! Remove the calving accounted for by accumulation
  grd%calving(:,:)=grd%calving(:,:)*grd%tmp(:,:)

  ! Do the same for heat (no separate classes needed)
  grd%calving_hflx(:,:)=grd%calving_hflx(:,:)*grd%tmp(:,:)
  grd%tmp(:,:)=bergs%dt*grd%calving_hflx(:,:)*grd%area(:,:)*(1.-grd%tmp(:,:))
  bergs%net_incoming_calving_heat_used=bergs%net_incoming_calving_heat_used+sum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  grd%stored_heat(:,:)=grd%stored_heat(:,:)+grd%tmp(:,:) ! +=bergs%dt*grd%calving_hflx(:,:)*grd%area(:,:)*(1.-remaining_dist)

end subroutine accumulate_calving

!> Generate icebergs from overflowing "buckets"
subroutine calve_icebergs(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  integer :: i,j,k,icnt,icntmax
  integer :: iNg, jNg ! Total number of points globally in i and j direction
  type(iceberg) :: newberg
  logical :: lret
  real :: xi, yj, ddt, calving_to_bergs, calved_to_berg, heat_to_bergs, heat_to_berg
  integer :: stderrunit
  real, pointer :: initial_mass, mass_scaling, initial_thickness, initial_width, initial_length
  logical :: allocations_done
  type(randomNumberStream) :: rns ! Random numbers for stochastic tidal parameterization
  real :: rx,ry

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  iNg=(grd%ieg-grd%isg+1) ! Total number of points globally in i direction
  jNg=(grd%jeg-grd%jsg+1) ! Total number of points globally in j direction

  grd%real_calving(:,:,:)=0.
  calving_to_bergs=0.
  heat_to_bergs=0.
  icntmax=0

  allocations_done=.false.

  do k=1, nclasses
    do j=grd%jsc, grd%jec
      do i=grd%isc, grd%iec
        ddt=0.; icnt=0
        if (grd%lat(i,j)<0.) then
          initial_mass=>bergs%initial_mass_s(k); mass_scaling=>bergs%mass_scaling_s(k)
          initial_thickness=>bergs%initial_thickness_s(k)
          initial_width=>bergs%initial_width_s(k); initial_length=>bergs%initial_length_s(k)
        else
          initial_mass=>bergs%initial_mass_n(k); mass_scaling=>bergs%mass_scaling_n(k)
          initial_thickness=>bergs%initial_thickness_n(k)
          initial_width=>bergs%initial_width_n(k); initial_length=>bergs%initial_length_n(k)
        endif
        do while (grd%stored_ice(i,j,k).ge.initial_mass*mass_scaling)
          newberg%lon=0.25*((grd%lon(i,j)+grd%lon(i-1,j-1))+(grd%lon(i-1,j)+grd%lon(i,j-1)))
          newberg%lat=0.25*((grd%lat(i,j)+grd%lat(i-1,j-1))+(grd%lat(i-1,j)+grd%lat(i,j-1)))
         !write(stderr(),*) 'KID, calve_icebergs: creating new iceberg at ',newberg%lon,newberg%lat
          lret=pos_within_cell(grd, newberg%lon, newberg%lat, i, j, xi, yj)
          if (.not.lret) then
            write(stderrunit,*) 'KID, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('KID, calve_icebergs', 'berg is not in the correct cell!', FATAL)
          endif
          if (debug.and.(xi<0..or.xi>1..or.yj<0..or.yj>1.)) then
            write(stderrunit,*) 'KID, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('KID, calve_icebergs', 'berg xi,yj is not correct!', FATAL)
          endif
          if (grd%msk(i,j)<0.5) then
            write(stderrunit,*) 'KID, calve_icebergs: WARNING!!! Iceberg born in land cell',i,j,newberg%lon,newberg%lat
            if (debug) call error_mesg('KID, calve_icebergs', 'Iceberg born in Land Cell!', FATAL)
          endif
          newberg%ine=i
          newberg%jne=j
          newberg%xi=xi
          newberg%yj=yj
          newberg%uvel=0.
          newberg%vvel=0.
          !--added by Alex:
          newberg%uvel_prev=0.
          newberg%vvel_prev=0.
          newberg%uvel_old=0.
          newberg%vvel_old=0.
          newberg%lon_prev=newberg%lon
          newberg%lat_prev=newberg%lat
          newberg%lon_old=newberg%lon
          newberg%lat_old=newberg%lat
          newberg%fl_k=0.
          !--added by Alon
          newberg%axn=0.
          newberg%ayn=0.
          newberg%bxn=0.
          newberg%byn=0.
          !---
          newberg%mass=initial_mass
          newberg%thickness=initial_thickness
          newberg%width=initial_width
          newberg%length=initial_length
          newberg%start_lon=newberg%lon
          newberg%start_lat=newberg%lat
          newberg%start_year=bergs%current_year
          newberg%id = generate_id(grd, i, j)
          newberg%start_day=bergs%current_yearday+ddt/86400.
          newberg%start_mass=initial_mass
          newberg%mass_scaling=mass_scaling
          newberg%mass_of_bits=0.
          newberg%mass_of_fl_bits=0.
          newberg%mass_of_fl_bergy_bits=0.
          newberg%halo_berg=0.
          newberg%static_berg=0.
          newberg%heat_density=grd%stored_heat(i,j)/grd%stored_ice(i,j,k) ! This is in J/kg

          if (bergs%mts) then
            if (.not. allocations_done) then
              if (.not. allocated(newberg%axn_fast)) allocate(newberg%axn_fast)
              if (.not. allocated(newberg%ayn_fast)) allocate(newberg%ayn_fast)
              if (.not. allocated(newberg%bxn_fast)) allocate(newberg%bxn_fast)
              if (.not. allocated(newberg%byn_fast)) allocate(newberg%byn_fast)
              if (.not. allocated(newberg%conglom_id)) allocate(newberg%conglom_id)
            endif
            newberg%axn_fast=0.; newberg%ayn_fast=0.; newberg%bxn_fast=0.; newberg%byn_fast=0.; newberg%conglom_id=0
          endif

          if (bergs%iceberg_bonds_on) then
            if (.not. allocations_done) then
              if (.not. allocated(newberg%n_bonds))  allocate(newberg%n_bonds)
            endif
            newberg%n_bonds=0
          endif

          if (bergs%dem) then
            if (.not. allocations_done) then
              if (.not. allocated(newberg%ang_vel)) allocate(newberg%ang_vel)
              if (.not. allocated(newberg%ang_accel)) allocate(newberg%ang_accel)
              if (.not. allocated(newberg%rot)) allocate(newberg%rot)
            endif
            newberg%ang_vel=0.; newberg%ang_accel=0.; newberg%rot=0.
          endif

          if (monitor_energy) then
            if (.not. allocations_done) then
              if (.not. allocated(newberg%Ee)) allocate(newberg%Ee)
              if (.not. allocated(newberg%Ed)) allocate(newberg%Ed)
              if (.not. allocated(newberg%Eext)) allocate(newberg%Eext)
              if (.not. allocated(newberg%Ee_contact)) allocate(newberg%Ee_contact)
              if (.not. allocated(newberg%Ed_contact)) allocate(newberg%Ed_contact)
              if (.not. allocated(newberg%Efrac)) allocate(newberg%Efrac)
              if (.not. allocated(newberg%Ee_temp)) allocate(newberg%Ee_temp)
              if (.not. allocated(newberg%Ed_temp)) allocate(newberg%Ed_temp)
              if (.not. allocated(newberg%Eext_temp)) allocate(newberg%Eext_temp)
              if (.not. allocated(newberg%Ee_contact_temp)) allocate(newberg%Ee_contact_temp)
              if (.not. allocated(newberg%Ed_contact_temp)) allocate(newberg%Ed_contact_temp)
            endif

            newberg%Ee=0.; newberg%Ed=0.; newberg%Eext=0.; newberg%Ee_contact=0.; newberg%Ed_contact=0.
            newberg%Efrac=0.; newberg%Ee_temp=0.; newberg%Ed_temp=0.;newberg%Eext_temp=0.;
            newberg%Ee_contact_temp=0.; newberg%Ed_contact_temp =0.
          endif

          if (bergs%fracture_criterion .ne. 'none' .and. (.not. bergs%dem)) then
            if (.not. allocations_done) then
              if (.not. allocated(newberg%accum_bond_rotation)) allocate(newberg%accum_bond_rotation)
            endif
            newberg%accum_bond_rotation=0.
          endif

          !interpolate gridded variables to new iceberg
          if (grd%tidal_drift>0.) then
            call getRandomNumbers(rns, rx)
            call getRandomNumbers(rns, ry)
            rx = 2.*rx - 1.; ry = 2.*ry - 1.
          endif
          call interp_flds(grd, newberg%lon, newberg%lat, i, j, xi, yj, rx, ry, newberg%uo, newberg%vo, newberg%ui, &
            newberg%vi, newberg%ua, newberg%va, newberg%ssh_x, newberg%ssh_y, newberg%sst, newberg%sss, newberg%cn, &
            newberg%hi, newberg%od)

          call add_new_berg_to_list(bergs%list(i,j)%first, newberg)
          calved_to_berg=initial_mass*mass_scaling ! Units of kg
          ! Heat content
          heat_to_berg=calved_to_berg*newberg%heat_density ! Units of J
          grd%stored_heat(i,j)=grd%stored_heat(i,j)-heat_to_berg
          heat_to_bergs=heat_to_bergs+heat_to_berg
          ! Stored mass
          grd%stored_ice(i,j,k)=grd%stored_ice(i,j,k)-calved_to_berg
          calving_to_bergs=calving_to_bergs+calved_to_berg
          grd%real_calving(i,j,k)=grd%real_calving(i,j,k)+calved_to_berg/bergs%dt
          ddt=ddt-bergs%dt*2./17. ! Minor offset to start day (negative offsets)
          icnt=icnt+1
          bergs%nbergs_calved=bergs%nbergs_calved+1
          if (grd%lat(i,j)<0.) then
            bergs%nbergs_calved_by_class_s(k)=bergs%nbergs_calved_by_class_s(k)+1
          else
            bergs%nbergs_calved_by_class_n(k)=bergs%nbergs_calved_by_class_n(k)+1
          endif

          allocations_done=.true.
        enddo
        icntmax=max(icntmax,icnt)
      enddo
    enddo
  enddo

  if (debug.and.icntmax>1) write(stderrunit,*) 'calve_icebergs: icnt=',icnt,' on',mpp_pe()

  bergs%net_calving_to_bergs=bergs%net_calving_to_bergs+calving_to_bergs
  bergs%net_heat_to_bergs=bergs%net_heat_to_bergs+heat_to_bergs

end subroutine calve_icebergs

  !> Calve footloose icebergs from parent berg
subroutine calve_fl_icebergs(bergs,pberg,k,l_b,fl_disp_x,fl_disp_y,berg_from_bits)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: pberg !< Parent berg
  type(real),intent(in) :: k !< Number of child bergs to calve
  type(real),intent(in) :: l_b !< The width, and 1/3 the length, of a child berg
  type(real) :: fl_disp_x !< Child berg x-displacement from parent berg
  type(real) :: fl_disp_y !< Child berg x-displacement from parent berg
  logical, optional, intent(in) :: berg_from_bits !< True to create new berg from footloose bits
  ! Local variables
  type(iceberg) :: cberg ! The new child berg
  type(icebergs_gridded), pointer :: grd
  logical :: lres, displace
  real :: Lfl, Wfl, Tfl, percent_fl

  ! For convenience
  grd=>bergs%grd

  !Initially, cberg%fl_k is set to -1 to mark the berg as a child berg.
  !If fl_k<0, the berg is not eligible to decay via footloose calving (see subroutine footloose_calving)
  !If fl_k==-1, the berg will also not interact with other berg elements. This is necessary because the child berg
  !may initially overlap its parent berg, which would otherwise result in unrealistic contact forces.
  !However, if (bergs%interactive_icebergs_on), a berg with fl_k==-1 will be set to fl_k=-2 once it is out of
  !contact range of any other berg. With fl_k=-2, interactive contact forces may now be applied

  !Set variable values for the new child berg:

  displace=bergs%displace_fl_bergs
  !Lat/lon and cell position for new berg:
  if (displace) then
    cberg%lon = pberg%lon + fl_disp_x
    cberg%lat = pberg%lat + fl_disp_y
    lres= find_cell(grd, cberg%lon, cberg%lat, cberg%ine, cberg%jne)

    !If new berg is not on current PE (computational domain), correct it so that it is.
    if (.not. lres) then
      !try the corners
      cberg%lon = pberg%lon - 0.5*pberg%length
      cberg%lat = pberg%lat - 0.5*pberg%width
      lres= find_cell(grd, cberg%lon, cberg%lat, cberg%ine, cberg%jne)
      if (.not. lres) then
        cberg%lon = pberg%lon - 0.5*pberg%length
        cberg%lat = pberg%lat + 0.5*pberg%width
        lres= find_cell(grd, cberg%lon, cberg%lat, cberg%ine, cberg%jne)
        if (.not. lres) then
          cberg%lon = pberg%lon + 0.5*pberg%length
          cberg%lat = pberg%lat + 0.5*pberg%width
          lres= find_cell(grd, cberg%lon, cberg%lat, cberg%ine, cberg%jne)
          if (.not. lres) then
            cberg%lon = pberg%lon + 0.5*pberg%length
            cberg%lat = pberg%lat - 0.5*pberg%width
            lres= find_cell(grd, cberg%lon, cberg%lat, cberg%ine, cberg%jne)
          endif
        endif
      endif
      if (.not. lres) then
        !if all else fails, give the child berg the same coords as the parent
        fl_disp_x=0.0; fl_disp_y=0.0; displace=.false.
        !call error_mesg('KID, calve_fl_icebergs', &
        !'corrected new berg position still not on current PE!', FATAL)
      else
        fl_disp_x=pberg%lon-cberg%lon; fl_disp_y=pberg%lat-cberg%lat
      endif
    endif
    if (displace) then
      !if new berg is located within a grounded cell, change its position to the same as the parent berg
      if (grd%area(cberg%ine,cberg%jne).eq.0.) then
        fl_disp_x=0.0; fl_disp_y=0.0; displace=.false.
      else
        lres=pos_within_cell(grd, cberg%lon, cberg%lat, cberg%ine, cberg%jne, cberg%xi, cberg%yj)
      endif
    endif
  endif

  if (.not. displace) then
    cberg%lon = pberg%lon
    cberg%lat = pberg%lat
    cberg%ine = pberg%ine
    cberg%jne = pberg%jne
    cberg%xi  = pberg%xi
    cberg%yj  = pberg%yj
  endif

  if (present(berg_from_bits)) then
    !use scaling for fl_bits to calculate new berg T,L,W, and M. Also affects mass scaling.
    call fl_bits_dimensions(bergs,pberg,Lfl,Wfl,Tfl)
    cberg%length = Lfl; cberg%width = Wfl; cberg%thickness = Tfl
    cberg%mass = Tfl * Lfl * Wfl * bergs%rho_bergs
    cberg%mass_scaling = k*bergs%new_berg_from_fl_bits_mass_thres/cberg%mass
    !the new berg will take a fraction of the parent berg footloose bergy bits mass as bergy bits mass
    percent_fl = (cberg%mass*cberg%mass_scaling)/(pberg%mass_of_fl_bits*pberg%mass_scaling)
    cberg%mass_of_bits = (percent_fl * pberg%mass_of_fl_bergy_bits*pberg%mass_scaling)/cberg%mass_scaling
    pberg%mass_of_fl_bergy_bits = (1-percent_fl)*pberg%mass_of_fl_bergy_bits
    pberg%mass_of_fl_bits=pberg%mass_of_fl_bits-k*bergs%new_berg_from_fl_bits_mass_thres/pberg%mass_scaling
  else
    cberg%length       = l_b*3.
    cberg%width        = l_b
    cberg%thickness    = pberg%thickness
    cberg%mass         = cberg%width * cberg%length * pberg%thickness * bergs%rho_bergs
    cberg%mass_scaling = pberg%mass_scaling * k !k is the number of icebergs cberg represents
    cberg%mass_of_bits = 0.0
  endif

  cberg%start_lon    = cberg%lon
  cberg%start_lat    = cberg%lat
  cberg%lon_prev     = pberg%lon_prev + fl_disp_x
  cberg%lat_prev     = pberg%lat_prev + fl_disp_y
  cberg%lon_old      = pberg%lon_old  + fl_disp_x
  cberg%lat_old      = pberg%lat_old  + fl_disp_y
  cberg%start_day    = bergs%current_yearday

  cberg%mass_of_fl_bits = 0.0
  cberg%mass_of_fl_bergy_bits = 0.0
  cberg%fl_k         = -1.0
  cberg%start_year   = bergs%current_year
  cberg%id           = generate_id(grd, pberg%ine, pberg%jne)
  cberg%halo_berg    = 0.0

  !always same values as parent:
  cberg%start_mass   = pberg%start_mass !used for tracking the different size classes calved from grounded ice, so set to pberg.
  cberg%uvel         = pberg%uvel
  cberg%vvel         = pberg%vvel
  cberg%axn          = pberg%axn
  cberg%ayn          = pberg%ayn
  cberg%bxn          = pberg%bxn
  cberg%byn          = pberg%byn
  cberg%uvel_prev    = pberg%uvel_prev
  cberg%vvel_prev    = pberg%vvel_prev
  cberg%uvel_old     = pberg%uvel_old
  cberg%vvel_old     = pberg%vvel_old
  cberg%heat_density = pberg%heat_density
  cberg%static_berg  = pberg%static_berg
  cberg%uo           = pberg%uo
  cberg%vo           = pberg%vo
  cberg%ui           = pberg%ui
  cberg%vi           = pberg%vi
  cberg%ua           = pberg%ua
  cberg%va           = pberg%va
  cberg%ssh_x        = pberg%ssh_x
  cberg%ssh_y        = pberg%ssh_y
  cberg%sst          = pberg%sst
  cberg%sss          = pberg%sss
  cberg%cn           = pberg%cn
  cberg%hi           = pberg%hi
  cberg%od           = pberg%od

  !optional variables
  if (bergs%mts) then
    allocate(cberg%axn_fast,cberg%ayn_fast,cberg%bxn_fast,cberg%byn_fast,cberg%conglom_id)
    cberg%axn_fast=pberg%axn_fast
    cberg%ayn_fast=pberg%ayn_fast
    cberg%bxn_fast=pberg%bxn_fast
    cberg%byn_fast=pberg%byn_fast
    cberg%conglom_id=pberg%conglom_id
  endif

  if (bergs%iceberg_bonds_on) then
    allocate(cberg%n_bonds)
    cberg%n_bonds=0
  endif

  if (bergs%dem) then
    allocate(cberg%ang_vel,cberg%ang_accel,cberg%rot)
    cberg%ang_vel=0.; cberg%ang_accel=0.; cberg%rot=0.
  endif

  if (monitor_energy) then
    allocate(cberg%Ee,cberg%Ed,cberg%Eext,cberg%Ee_contact,cberg%Ed_contact,cberg%Efrac,&
      cberg%Ee_temp,cberg%Ed_temp,cberg%Eext_temp,cberg%Ee_contact_temp,cberg%Ed_contact_temp)
    cberg%Ee=0.; cberg%Ed=0.; cberg%Eext=0.; cberg%Ee_contact=0.; cberg%Ed_contact=0.
    cberg%Efrac=0.; cberg%Ee_temp=0.; cberg%Ed_temp=0.;cberg%Eext_temp=0.;
    cberg%Ee_contact_temp=0.; cberg%Ed_contact_temp =0.
  endif

  if (bergs%fracture_criterion .ne. 'none' .and. (.not. bergs%dem)) then
    allocate(cberg%accum_bond_rotation)
    cberg%accum_bond_rotation=0.
  endif

  call add_new_berg_to_list(bergs%list(cberg%ine,cberg%jne)%first, cberg)
end subroutine calve_fl_icebergs

!> Evolves icebergs forward by updating velocity and position with a multiple-time-step Velocity Verlet
!! scheme. Experimental option: each short/long step can be iterated until a given convergence tolerance
!! is met, which better enforces momentum conservation during collision...this iterative scheme may or may
!! not be needed to yield consistent fracture behavior. The short steps can also be evaluated explicitly,
!! which is the default for dem mode.
subroutine evolve_icebergs_mts(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg
  real :: uveln, vveln, lonn, latn
  real :: axn, ayn, bxn, byn
  real :: xi, yj, rx, ry
  integer :: i, j, k
  integer :: grdi, grdj
  integer :: stderrunit
  logical :: only_interactive_forces
  real :: lon1, lat1, dxdl1, dydl
  real :: uvel1, vvel1, uvel2, vvel2, uvel3, vvel3
  real :: xddot1, yddot1, xdot2, ydot2, xdot3, ydot3, xdotn, ydotn
  real :: u2, v2, x1, y1, xn, yn
  real :: ax1,ay1
  real :: dt, dt_2
  logical :: on_tangential_plane, error_flag, bounced
  !energy stuff
  real :: dx,dy,M,Fex,Fey,Fex_prev,Fey_prev,Fdx_prev,Fdy_prev,Fdx,Fdy,Fx_tot,Fy_tot,Fx_tot_prev,Fy_tot_prev
  logical :: save_bond_energy
  real :: Fex_i, Fey_i, Fec_x, Fec_y, Fdc_x, Fdc_y, ustar, vstar, uveln_star, vveln_star
  integer :: ii,jj,maxii,minii,maxjj,minjj
  real :: usum,usum1,usum2,normchange,denom
  logical :: finished,last_iter

  ! Multiple Time Step Velocity Verlet:
  ! NOTE: To avoid the need for computationally expensive transfers between processors
  ! at the end of each MTS sub-step, each PE updates a complete copy of any
  ! conglomerate of bonded iceberg elements that spans over the computational domain of
  ! several processors (see subroutine transfer_mts_bergs). Berg elements that lie within
  ! the contact distance of the conglomerate are also included, as to account for collision
  ! forces. The conglomerate iceberg elements with coords that lie outside of a PE's domain
  ! are assigned to the nearest halo cell, but the coords are not changed.

  ! Get the stderr unit number
  stderrunit = stderr()
  only_interactive_forces=bergs%only_interactive_forces

  ! For convenience
  grd=>bergs%grd

  !Checking if everything is ok:
  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec ! just for computational domain
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      if (berg%static_berg .lt. 0.5 .and. berg%halo_berg .lt. 0.5) then
        if (.not. is_point_in_cell(bergs%grd, berg%lon, berg%lat, berg%ine, berg%jne) ) then
          write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
          do j=grd%jed,grd%jsd,-1
            write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
          enddo
          write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
          do j=grd%jed,grd%jsd,-1
            write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
          enddo
          call print_berg(stderrunit, berg, 'evolve_icebergs_mts, berg is not in proper starting cell')
          write(stderrunit,'(a,i3,2(i4,3f8.2))') 'evolve_icebergs_mts: pe,lon/lat(i,j)=', mpp_pe(), &
            berg%ine,berg%lon,grd%lon(berg%ine-1,berg%jne-1),grd%lon(berg%ine,berg%jne), &
            berg%jne,berg%lat,grd%lat(berg%ine-1,berg%jne-1),grd%lat(berg%ine,berg%jne)
          if (debug) call error_mesg('KID, evolve_icebergs_mts','berg is in wrong starting cell!',FATAL)
        endif
        if (debug) call check_position(grd, berg, 'evolve_icebergs_mts (top)')
      end if
      berg=>berg%next
    enddo
  enddo; enddo

  ! PART 1: solve for V_n+1 (from the previous timestep)
  bergs%mts_part=1
  dt = bergs%dt
  dt_2=0.5*dt
  rx = 0.0; ry = 0.0 !not needed, because interp_gridded_fields_to_bergs was already called
  ii = 0
  save_bond_energy=.false.
  usum=0.0; usum1=0.0; usum2=0.0
  finished=.false.

  if (bergs%force_convergence) then
    last_iter=.false.
  else
    last_iter=.true.
  endif

  do while (.not. finished)

    ii=ii+1

    do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg)) ! loop over all bergs
        !only evolve non-static bergs that overlap, or are part of a conglom that overlaps, the computational domain:
        if (berg%static_berg .lt. 0.5 .and. (berg%conglom_id.ne.0 .or. bergs%force_convergence)) then
          latn = berg%lat ;   lonn = berg%lon
          uvel1=berg%uvel ;   vvel1=berg%vvel !V_n (previous cycle)
          axn  = 0.0      ;   ayn  = 0.0 !note these are redefined at the start of subroutine accel_mts
          bxn  = 0.0      ;   byn  = 0.0
          i=berg%ine      ;   j=berg%jne
          xi=berg%xi      ;   yj=berg%yj

          call accel_mts(bergs, berg, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, rx, ry, &
            ax1, ay1, axn, ayn, bxn, byn, save_bond_energy,Fec_x, Fec_y, Fdc_x, Fdc_y)

          if (monitor_energy .and. last_iter) &
            call mts_energy_part_1(bergs,berg,ax1,ay1,axn,ayn,bxn,byn,Fec_x,Fec_y,Fdc_x,Fdc_y)

          ! Saving all the iceberg variables
          berg%axn=axn; berg%ayn=ayn
          berg%bxn=bxn; berg%byn=byn

          if (bergs%force_convergence) then
            berg%uvel_prev=berg%uvel+(dt*ax1); berg%vvel_prev=berg%vvel+(dt*ay1) !the new velocity
            if (ii==1) usum=usum+berg%uvel_old**2 + berg%vvel_old**2
            usum1=usum1+berg%uvel_prev**2+berg%vvel_prev**2
            usum2=usum2+(berg%uvel_prev-berg%uvel_old)**2+(berg%vvel_prev-berg%vvel_old)**2
          else
            berg%uvel=berg%uvel+(dt*ax1); berg%vvel=berg%vvel+(dt*ay1)

            !The final velocity from the previous cycle, which can be used
            !during post-processing to calculate kinetic energy and momentum
            berg%uvel_prev=berg%uvel    ; berg%vvel_prev=berg%vvel
          endif
        endif
        berg=>berg%next
      enddo ! loop over all bergs
    enddo; enddo

    !update uvel_old and vvel_old for damping interaction with other elements
    if (bergs%force_convergence) then
      do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
        berg=>bergs%list(grdi,grdj)%first
        do while (associated(berg)) ! loop over all bergs
          if (berg%static_berg .lt. 0.5) then
            berg%uvel_old=berg%uvel_prev; berg%vvel_old=berg%vvel_prev
          endif
          berg=>berg%next
        enddo
      enddo; enddo
    endif

    if (bergs%force_convergence .and. .not. last_iter) then
      if (ii>1) then
        denom=sqrt(usum)+sqrt(usum1)
        if (denom>0) then
          normchange=2.0*sqrt(usum2)/denom
        else
          normchange=0.0
        endif
        if (normchange<bergs%convergence_tolerance) last_iter=.true.
      endif
      usum=usum1 !previous norm (squared)
    else
      finished=.true.
    endif

    if (last_iter .and. .not. monitor_energy) finished=.true.

    usum1=0.0; usum2=0.0

  enddo


  if (.not. bergs%dem) then
    call update_bond_angles(bergs)
    call update_and_break_bonds(bergs)
  else
    call break_bonds_dem(bergs)
  endif

  !PART 2: X_0 and V_0, before fast sub-steps
  !X_0=X_n (no change)
  !update V_0 on uvel_old and vvel_old
  !update lat and lon with implicit slow force component
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      if (berg%static_berg .lt. 0.5 .and. berg%conglom_id.ne.0) then
        berg%lat_prev = berg%lat; berg%lon_prev = berg%lon
        berg%uvel=berg%uvel_prev; berg%vvel=berg%vvel_prev

        if (new_mts) then
          berg%uvel = berg%uvel + dt_2*(berg%axn+berg%bxn)
          berg%vvel = berg%vvel + dt_2*(berg%ayn+berg%byn)
        else
          berg%uvel = berg%uvel + dt_2*berg%axn
          berg%vvel = berg%vvel + dt_2*berg%ayn
          berg%lat = berg%lat + 0.5*(dt**2)*berg%byn
          berg%lon = berg%lon + 0.5*(dt**2)*berg%bxn
        endif

        berg%uvel_old=berg%uvel; berg%vvel_old=berg%vvel

        if (bergs%force_convergence) then
          berg%axn=berg%axn_fast; berg%ayn=berg%ayn_fast
          berg%bxn=berg%bxn_fast; berg%byn=berg%byn_fast
        endif

        !if (monitor energy): already did energy update in part I
      endif
      berg=>berg%next
    enddo
  enddo; enddo


  ! PART 3: MTS fast sub-steps (bonded interactions only):
  ! Position and velocity are updated by:
  ! X2 = X1+dt*V1+((dt^2)/2)*a_k +((dt^2)/2)*b_k
  ! V2 = V1+dt/2*a_k +dt/2*a_kp1 +dt*b_k+1

  bergs%only_interactive_forces=.true. !so that external forcings are ignored in subroutine accel
  bergs%mts_part=3 !so that only bonded interactions are considered in subroutine interactive_force
  dt = bergs%mts_fast_dt
  dt_2=0.5*dt
  maxjj=0; minjj=huge(1)

  do k = 1,bergs%mts_sub_steps ! loop over sub-steps

    do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied ! update positions
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg)) ! loop over all bergs
        if (berg%static_berg .lt. 0.5 .and. berg%conglom_id.ne.0) then

          on_tangential_plane=.false.
          if ((berg%lat>89.) .and. (grd%grid_is_latlon)) on_tangential_plane=.true.

          lon1=berg%lon; lat1=berg%lat
          if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1,berg%id)
          !dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
          !dydl=r180_pi/Rearth
          call convert_from_meters_to_grid(lat1,grd%grid_is_latlon ,dxdl1,dydl)
          uvel1=berg%uvel; vvel1=berg%vvel

          ! Loading past accelerations
          axn=berg%axn_fast; ayn=berg%ayn_fast; bxn=berg%bxn_fast; byn=berg%byn_fast

          ! Velocities used to update the position
          uvel2=uvel1+(dt_2*axn)+(dt_2*bxn); vvel2=vvel1+(dt_2*ayn)+(dt_2*byn)

          if (on_tangential_plane) call rotvec_to_tang(lon1,uvel2,vvel2,xdot2,ydot2)
          u2=uvel2*dxdl1; v2=vvel2*dydl

          ! Solving for new position
          if (on_tangential_plane) then
            xn=x1+(dt*xdot2) ; yn=y1+(dt*ydot2)
            call rotpos_from_tang(xn,yn,lonn,latn)
          else
            lonn=lon1+(dt*u2) ; latn=lat1+(dt*v2)
          endif

          berg%lon=lonn      ;  berg%lat=latn
          berg%lon_old=lonn  ;  berg%lat_old=latn

          !Set uvel_old to ustar for berg interactions?
          if (new_mts) then
            berg%uvel_old=berg%uvel+dt_2*(berg%axn_fast+berg%bxn_fast)
            berg%vvel_old=berg%vvel+dt_2*(berg%ayn_fast+berg%bxn_fast)
          else
            berg%uvel_old=berg%uvel+dt_2*(berg%axn_fast)
            berg%vvel_old=berg%vvel+dt_2*(berg%ayn_fast)
          endif
        end if
        berg=>berg%next
      enddo !loop over all bergs
    enddo;enddo !update positions

    jj=0
    usum=0.0; usum1=0.0; usum2=0.0
    finished=.false.
    save_bond_energy=.false.
    if (bergs%force_convergence .and. (.not. bergs%explicit_inner_mts)) then
      last_iter=.false.
    else
      last_iter=.true.
    endif

    do while (.not. finished)

      jj=jj+1

      if (monitor_energy .and. last_iter) save_bond_energy=.true.

    ! update velocities
    do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg)) ! loop over all bergs
        if (berg%static_berg .lt. 0.5 .and. berg%conglom_id.ne.0) then

          latn = berg%lat ;   lonn = berg%lon
          axn  = berg%axn_fast ;   ayn  = berg%ayn_fast
          bxn  = berg%bxn_fast ;   byn  = berg%byn_fast
          uvel1=berg%uvel ;   vvel1=berg%vvel
          i=berg%ine      ;   j=berg%jne
          xi=berg%xi      ;   yj=berg%yj

          if (new_mts) then
            axn = axn+bxn; ayn=ayn+byn
          endif

          ! Turn the velocities into u_star, v_star.(uvel3 is v_star) - (possible issues with tangent plane?)
          uvel3=uvel1+(dt_2*axn)                  !(Stern et al 2017, Eq B4)
          vvel3=vvel1+(dt_2*ayn)

          if (monitor_energy .and. last_iter) call mts_energy_part_3(bergs,berg,uvel3,vvel3,wt=1)

          if (bergs%explicit_inner_mts) then
            call accel_explicit_inner_mts(bergs, berg, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, &
              ax1, ay1, axn, ayn, save_bond_energy)
              bxn=0.; byn=0.
          else
            call accel_mts(bergs, berg, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, rx, ry, &
              ax1, ay1, axn, ayn, bxn, byn, save_bond_energy)
          endif


          ! Solving for the new velocity (Stern et al 2017, Eqn B5)
          on_tangential_plane=.false.
          if ((berg%lat>89.) .and. (bergs%grd%grid_is_latlon)) on_tangential_plane=.true.
          if (on_tangential_plane) then
            call rotvec_to_tang(lonn,uvel3,vvel3,xdot3,ydot3)
            call rotvec_to_tang(lonn,ax1,ay1,xddot1,yddot1)
            xdotn=xdot3+(dt*xddot1); ydotn=ydot3+(dt*yddot1)
            call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
          else
            !uvel3, vvel3 become uveln and vveln after they are put into lat/lon co-ordinates
            uveln=uvel3+(dt*ax1); vveln=vvel3+(dt*ay1)
          endif

          if (bergs%force_convergence .and. (.not. bergs%explicit_inner_mts)) then
            if (jj==1) usum=usum+berg%uvel_old**2 + berg%vvel_old**2
            usum1=usum1 + uveln**2 + vveln**2
            usum2=usum2+(uveln-berg%uvel_old)**2+(vveln-berg%vvel_old)**2
          endif

          ! Saving all the iceberg variables.
          berg%axn_fast=axn; berg%ayn_fast=ayn
          berg%bxn_fast=bxn; berg%byn_fast=byn
          berg%uvel=uveln  ; berg%vvel=vveln

          if (monitor_energy .and. last_iter) call mts_energy_part_3(bergs,berg,uvel3,vvel3,wt=2)

        endif
        berg=>berg%next
      enddo ! loop over all bergs
    enddo; enddo ! update velocities

    if (bergs%force_convergence .and. .not. last_iter) then
      if (jj>1) then
        denom=sqrt(usum)+sqrt(usum1)
        if (denom>0) then
          normchange=2.0*sqrt(usum2)/denom
        else
          normchange=0.0
        endif
        if (normchange<bergs%convergence_tolerance) last_iter=.true.
      endif
      usum=usum1 !previous norm (squared)
    else
      finished=.true.
    endif

    if (last_iter .and. .not. monitor_energy) finished=.true.

    if (bergs%force_convergence .and. (.not. finished)) then
      do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
        berg=>bergs%list(grdi,grdj)%first
        do while (associated(berg)) ! loop over all bergs
          if (berg%static_berg .lt. 0.5 .and. berg%conglom_id.ne.0) then

            !update old velocities for interactions
            berg%uvel_old=berg%uvel; berg%vvel_old=berg%vvel

            if (new_mts) then
              berg%uvel=berg%uvel-dt_2*(berg%axn_fast+berg%bxn_fast)-dt_2*(berg%axn+berg%bxn)
              berg%vvel=berg%vvel-dt_2*(berg%ayn_fast+berg%byn_fast)-dt_2*(berg%ayn+berg%byn)
            else
              berg%uvel=berg%uvel-dt_2*(berg%axn_fast+2.0*berg%bxn_fast)-dt_2*berg%axn
              berg%vvel=berg%vvel-dt_2*(berg%ayn_fast+2.0*berg%byn_fast)-dt_2*berg%ayn
            endif
            berg%axn_fast=berg%axn; berg%ayn_fast=berg%ayn
            berg%bxn_fast=berg%bxn; berg%byn_fast=berg%byn
          endif
          berg=>berg%next
        enddo
      enddo; enddo
    endif

    usum1=0.0; usum2=0.0

    if (jj<minjj .and. finished) minjj=jj
    if (jj>maxjj .and. finished) maxjj=jj
  enddo

    !update 'old' velocities used for interactions
    do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg)) ! loop over all bergs
        if (berg%static_berg .lt. 0.5 .and. berg%conglom_id.ne.0) then
          berg%uvel_old=berg%uvel; berg%vvel_old=berg%vvel

          if (bergs%dem) then
            berg%ang_vel=berg%ang_vel+dt*berg%ang_accel
            berg%rot=berg%rot+dt*berg%ang_vel
          endif

          if (bergs%force_convergence) then
            berg%axn=berg%axn_fast; berg%ayn=berg%ayn_fast
            berg%bxn=berg%bxn_fast; berg%byn=berg%byn_fast
          endif
        endif
        berg=>berg%next
      enddo ! loop over all bergs
    enddo; enddo ! update 'old' velocities
  enddo ! loop over sub-steps

  !reset interactive_forces
  bergs%only_interactive_forces=only_interactive_forces

  !call update_and_break_bonds(bergs)

  !update indices
  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    !no reason to worry about halo/conglomerate bergs here
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      if (berg%static_berg .lt. 0.5 .and. berg%halo_berg<1) then
        berg%uvel_old=berg%uvel; berg%vvel_old=berg%vvel
        uveln=berg%uvel; vveln=berg%vvel
        lonn=berg%lon  ; latn=berg%lat
        i=berg%ine     ; j=berg%jne
        xi=berg%xi     ; yj=berg%yj
        ! finalize new iceberg positions and index
        call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag, berg%id)
        berg%lon=lonn      ;  berg%lat=latn
        berg%lon_old=lonn  ;  berg%lat_old=latn
        berg%ine=i    ;  berg%jne=j
        berg%xi=xi    ;  berg%yj=yj
      endif
      berg=>berg%next
    enddo ! loop over all bergs
  enddo; enddo ! update 'old' velocities


  ! for testing:
  ! maxii=ii; minii=ii
  ! call mpp_max(maxii); call mpp_min(minii); call mpp_max(maxjj); call mpp_min(minjj)
  !if (maxii>2 .or. maxjj>2) print *,'outer max iters, inner max iters',maxii,maxjj

end subroutine evolve_icebergs_mts


!> Energy calculations for MTS_part==1
subroutine mts_energy_part_1(bergs, berg, ax1, ay1, axn, ayn, bxn, byn, Fec_x, Fec_y, Fdc_x, Fdc_y)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Berg for which to do calculations
  real, intent(in) :: ax1 !< Zonal acceleration [m s-2]
  real, intent(in) :: ay1 !< Meridional acceleration [m s-2]
  real, intent(in) :: axn !< Explicit estimate of zonal acceleration [m s-2]
  real, intent(in) :: ayn !< Explicit estimate of meridional acceleration [m s-2]
  real, intent(in) :: bxn !< Implicit estimate of zonal acceleration [m s-2]
  real, intent(in) :: byn !< Implicit estimate of meridional acceleration [m s-2]
  real, intent(in) :: Fec_x !< Zonal elastic force from collision with other conglomerate [N]
  real, intent(in) :: Fec_y !< Meridional elastic force from collision with other conglomerate [N]
  real, intent(in) :: Fdc_x !< !< Zonal damping force from collision with other conglomerate [N]
  real, intent(in) :: Fdc_y !< !< Meridional damping force from collision with other conglomerate [N]
  ! Local
  real :: dt, M, Fex, Fey, Fex_i, Fey_i
  real :: uveln, vveln, uveln_star, vveln_star, ustar, vstar

  dt=bergs%dt

  !Again, MTS part 1 solves for V_n+1 from the previous timestep
  !"Temp" energy variables are from the previous cycle
  !At the end of MTS part 1, the non-"temp" energy variables are saved for post-processing/analysis
  !see plot_energy.py

  !derived from Mayrhofer,Hager, and Kloss 2017 (MHK17) and Asmar et al 2003
  M = berg%mass
  if (new_mts) then
    Fex = M * 2*ax1; Fey = M * 2*ay1 !external forces
    Fex = Fex - Fec_x - Fdc_x; Fey = Fey - Fec_y - Fdc_y !external forces minus collisional

    uveln=berg%uvel+(dt*ax1); vveln=berg%vvel+(dt*ay1)
    ustar=berg%uvel; vstar=berg%vvel

    !dEext=0.5*(V_{i+1} + V_{i+0.5})*F_{ext,i+1}*dt/2
    berg%Eext=berg%Eext_temp - 0.25*dt*(&
      (uveln + ustar)*Fex + (vveln + vstar)*Fey)

    berg%Ee_contact=berg%Ee_contact_temp - 0.25*dt*(&
      (uveln + ustar)*Fec_x + (vveln + vstar)*Fec_y)

    berg%Ed_contact=berg%Ed_contact_temp - 0.25*dt*(&
      (uveln + ustar)*Fdc_x + (vveln + vstar)*Fdc_y)

    !The following is for part II
    !V{i+1.5}
    uveln_star = uveln + dt*ax1; vveln_star = vveln + dt*ay1

    berg%Eext_temp=berg%Eext - 0.25*dt*(&
      (uveln + uveln_star)*Fex + (vveln + vveln_star)*Fey)

    berg%Ee_contact_temp=berg%Ee_contact - 0.25*dt*(&
      (uveln + uveln_star)*Fec_x + (vveln + vveln_star)*Fec_y)

    berg%Ed_contact_temp=berg%Ed_contact - 0.25*dt*(&
      (uveln + uveln_star)*Fdc_x + (vveln + vveln_star)*Fdc_y)

  else
    Fex = M * axn; Fey = M * ayn !explicit forces
    Fex = Fex - Fec_x; Fey = Fey - Fec_y !explicit forces minus collisional elastic
    Fex_i = M * bxn;     Fey_i= M * byn !implicit forces
    Fex_i = Fex_i - Fdc_x; Fey_i = Fey_i - Fdc_y !implicit forces minus collisional damping

    uveln=berg%uvel+(dt*ax1); vveln=berg%vvel+(dt*ay1)
    ustar=berg%uvel; vstar=berg%vvel

    !explicit:
    !dEext=0.5*(V_{i+1} + V_{i+0.5})*F_{ext,i+1}*dt/2
    berg%Eext=berg%Eext_temp - 0.25*dt*(&
      (uveln + ustar)*Fex + (vveln + vstar)*Fey)

    berg%Ee_contact=berg%Ee_contact_temp - 0.25*dt*(&
      (uveln + ustar)*Fec_x + (vveln + vstar)*Fec_y)

    berg%Eext=berg%Eext - 0.5*dt*(&
      (uveln + ustar)*Fex_i + (vveln + vstar)*Fey_i)

    berg%Ed_contact=berg%Ed_contact - 0.5*dt*(&
      (uveln + ustar)*Fdc_x + (vveln + vstar)*Fdc_y)

    !The following is for part II
    !V{i+1.5}
    uveln_star = uveln + 0.5*dt*axn; vveln_star = vveln + 0.5*dt*ayn

    berg%Eext_temp=berg%Eext - 0.25*dt*(&
      (uveln + uveln_star)*Fex + (vveln + vveln_star)*Fey)

    berg%Ee_contact_temp=berg%Ee_contact - 0.25*dt*(&
      (uveln + uveln_star)*Fec_x + (vveln + vveln_star)*Fec_y)
  endif

  berg%Ee=berg%Ee_temp; berg%Ed=berg%Ed_temp
end subroutine mts_energy_part_1

!> Energy calculations for MTS_part==3
subroutine mts_energy_part_3(bergs,berg,ustar,vstar,wt)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Berg to do calculations for
  real,intent(in) :: ustar !< u-component of half-step velocity [m s-1]
  real,intent(in) :: vstar !< v-component of half-step velocity [m s-1]
  integer,intent(in) :: wt !< integrate from time i to i+1/2 (wt==1) or i+1/2 to i+1 (wt==2)
  ! Local
  real :: dt, M, Fex, Fey, Fdx, Fdy


  !MTS sub-step forces
  dt = bergs%mts_fast_dt
  M = berg%mass
  Fex = M * berg%axn_fast; Fey = M * berg%ayn_fast !elastic forces

  !for i to i+1/2 or i+1/2 to i+1
  berg%Ee_temp=berg%Ee_temp - 0.25*dt*(&
    (ustar + berg%uvel)*Fex + &
    (vstar + berg%vvel)*Fey)

  if (new_mts) then   !for i to i+1/2 or i+1/2 to i+1
    Fdx = M * berg%bxn_fast; Fdy = M * berg%byn_fast !damping forces
    berg%Ed_temp=berg%Ed_temp - 0.25*dt*(&
      (ustar + berg%uvel)*Fdx + &
      (vstar + berg%vvel)*Fdy)
  elseif (wt==2) then !for i+1/2 to i+1
    Fdx = M * berg%bxn_fast; Fdy = M * berg%byn_fast !damping forces
    berg%Ed_temp=berg%Ed_temp - 0.5*dt*(&
      (ustar + berg%uvel)*Fdx + &
      (vstar + berg%vvel)*Fdy)
  endif
end subroutine mts_energy_part_3

!> Evolves icebergs forward by updating velocity and position with a time-stepping scheme
subroutine evolve_icebergs(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg
  real :: uveln, vveln, lonn, latn
  real :: axn, ayn, bxn, byn          ! Added by Alon - explicit and implicit accelerations from the previous step
  real :: xi, yj, rx, ry
  integer :: i, j
  integer :: grdi, grdj
  integer :: stderrunit
  logical :: bounced, interactive_icebergs_on, Runge_not_Verlet
  type(randomNumberStream) :: rns ! Random numbers for stochastic tidal parameterization
  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd
  rx = 0.
  ry = 0.

  interactive_icebergs_on=bergs%interactive_icebergs_on
  Runge_not_Verlet=bergs%Runge_not_Verlet

  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    berg=>bergs%list(grdi,grdj)%first
    if (associated(berg) .and. grd%tidal_drift>0.) then
      ! Seed random numbers based on space and "time"
      rns = initializeRandomNumberStream( grdi + 10000*grdj + &
                                          int( 16384.*abs( sin(262144.*grd%ssh(grdi,grdj)) ) ) )
    endif
    do while (associated(berg)) ! loop over all bergs
      if (berg%static_berg .lt. 0.5) then  !Only allow non-static icebergs to evolve

        !Checking it everything is ok:
        if (.not. is_point_in_cell(bergs%grd, berg%lon, berg%lat, berg%ine, berg%jne) ) then
          write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
          do j=grd%jed,grd%jsd,-1
            write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
          enddo
          write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
          do j=grd%jed,grd%jsd,-1
            write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
          enddo
          call print_berg(stderrunit, berg, 'evolve_iceberg, berg is not in proper starting cell')
          write(stderrunit,'(a,i3,2(i4,3f8.2))') 'evolve_iceberg: pe,lon/lat(i,j)=', mpp_pe(), &
                   berg%ine,berg%lon,grd%lon(berg%ine-1,berg%jne-1),grd%lon(berg%ine,berg%jne), &
                   berg%jne,berg%lat,grd%lat(berg%ine-1,berg%jne-1),grd%lat(berg%ine,berg%jne)
          if (debug) call error_mesg('KID, evolve_iceberg','berg is in wrong starting cell!',FATAL)
        endif
        if (debug) call check_position(grd, berg, 'evolve_iceberg (top)')

        if (grd%tidal_drift>0.) then
          call getRandomNumbers(rns, rx)
          rx = 2.*rx - 1.
          call getRandomNumbers(rns, ry)
          ry = 2.*ry - 1.
        endif

          !Time stepping schemes:
          if (Runge_not_Verlet) then
            call Runge_Kutta_stepping(bergs, berg, axn, ayn, bxn, byn, uveln, vveln, &
                                      lonn, latn, i, j, xi, yj, rx, ry)
          endif
          if (.not.Runge_not_Verlet) then
            call verlet_stepping(bergs, berg, axn, ayn, bxn, byn, uveln, vveln, rx, ry)
          endif

          !Used for testing the ocean response to fixed iceberg motion.
          if (bergs%override_iceberg_velocities) then
            uveln  = bergs%u_override
            vveln  = bergs%v_override
          endif

        ! Saving all the iceberg variables.
        berg%axn=axn
        berg%ayn=ayn
        berg%bxn=bxn
        berg%byn=byn
        berg%uvel=uveln
        berg%vvel=vveln

        if (Runge_not_Verlet) then
          berg%lon=lonn  ;   berg%lat=latn
          berg%ine=i     ;   berg%jne=j
          berg%xi=xi     ;   berg%yj=yj
        else
          if (.not. interactive_icebergs_on) call update_verlet_position(bergs,berg)
        endif

        !call interp_flds(grd, berg%lon, berg%lat, i, j, xi, yj, berg%uo, berg%vo, berg%ui, &
        !berg%vi, berg%ua, berg%va, berg%ssh_x, berg%ssh_y, berg%sst,berg%od)
        !if (debug) call print_berg(stderr(), berg, 'evolve_iceberg, final posn.')
        if (debug) call check_position(grd, berg, 'evolve_iceberg (bot)')
      endif
      berg=>berg%next
    enddo ! loop over all bergs
  enddo ; enddo

  ! When we are using interactive icebergs, we update the (old) iceberg positions and velocities in a second loop, all together (to make code order invarient)
  if (interactive_icebergs_on) then
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg)) ! loop over all bergs
        if (berg%static_berg .lt. 0.5) then  !Only allow non-static icebergs to evolve
         if (.not. Runge_not_Verlet) call update_verlet_position(bergs,berg)

         !Updating old velocities (for use in iceberg interactions)
          berg%uvel_old=berg%uvel
          berg%vvel_old=berg%vvel
          berg%lon_old=berg%lon  ! lon_old, lat_old are not really needed for Verlet. But are needed for RK
          berg%lat_old=berg%lat
        endif
        berg=>berg%next
      enddo ! loop over all bergs
    enddo ; enddo
  endif

end subroutine evolve_icebergs

!> Calculate explicit and implicit accelerations, and new velocity, using the Verlet method
subroutine verlet_stepping(bergs,berg, axn, ayn, bxn, byn, uveln, vveln, rx, ry)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer, intent(inout) :: berg !< Iceberg
  real, intent(out) :: axn !< Explicit zonal acceleration (m/s2)
  real, intent(out) :: ayn !< Explicit meridional acceleration (m/s2)
  real, intent(out) :: bxn !< Implicit zonal acceleration (m/s2)
  real, intent(out) :: byn !< Implicit meridional acceleration (m/s2)
  real, intent(out) :: uveln !< New zonal velocity (m/s)
  real, intent(out) :: vveln !< New meridional velocity (m/s)
  real, intent(in) :: rx !< Random number between -1 and 1 for use in x-component of stochastic tidal parameterization
  real, intent(in) :: ry !< Random number between -1 and 1 for use in y-component of stochastic tidal parameterization
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: lonn, latn
  real :: uvel1, vvel1,  uvel2, vvel2, uvel3, vvel3
  real :: ax1, ay1
  real :: x1,  y1, xddot1, yddot1, xi, yj
  real :: xdot3, ydot3
  real :: xdotn, ydotn
  real :: dt, dt_2, dt_6, dydl
  real :: orientation
  logical :: bounced, on_tangential_plane, error_flag
  integer :: i, j
  integer :: stderrunit

  ! Initialize variables

  ! In this scheme a_n and b_n are saved from the previous timestep, giving the explicit and implicit parts of the acceleration, and a_np1, b_np1 are for the next time step
  ! Note that ax1=a_np1/2 +b_np1, as calculated by the acceleration subrouting
  ! Positions and velocity is updated by
  ! X2 = X1+dt*V1+((dt^2)/2)*a_n +((dt^2)/2)*b_n = X1+dt*u_star +((dt^2)/2)*b_n
  ! V2 = V1+dt/2*a_n +dt/2*a_np1 +dt*b_n = u_star + dt/2*a_np1 + dt*b_np1 = u_star +dt*ax

  !*************************************************************************************************

  ! Get the stderr unit number
  stderrunit = stderr()
  ! For convenience
  grd=>bergs%grd
  ! Common constants
  dt=bergs%dt
  dt_2=0.5*dt

  orientation=bergs%initial_orientation
  if ((bergs%iceberg_bonds_on) .and. (bergs%rotate_icebergs_for_mass_spreading)) call find_orientation_using_iceberg_bonds(grd,berg,orientation)

  lonn = berg%lon ;   latn = berg%lat
  axn  = berg%axn ;   ayn  = berg%ayn
  bxn=   berg%bxn ;   byn  = berg%byn
  uvel1=berg%uvel ;   vvel1=berg%vvel
  i=berg%ine      ;   j=berg%jne
  xi=berg%xi      ;   yj=berg%yj

  berg%uvel_prev=berg%uvel-dt_2*berg%bxn; berg%vvel_prev=berg%vvel-dt_2*berg%byn

  ! Turn the velocities into u_star, v_star.(uvel3 is v_star) - Alon (not sure how this works with tangent plane)
  uvel3=uvel1+(dt_2*axn)                  !Alon (Stern et al 2017, Eq B4)
  vvel3=vvel1+(dt_2*ayn)                  !Alon

  ! Note, the mass scaling is equal to 1 (rather than 0.25 as in RK), since
  ! this is only called once in Verlet stepping.
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 1.0*berg%mass_scaling,&
    berg%length*berg%width, berg%thickness, addfootloose=.true.)

  ! Calling the acceleration   (note that the velocity is converted to u_star inside the accel script)
  call accel(bergs, berg, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, rx, ry, ax1, ay1, axn, ayn, bxn, byn) !axn, ayn, bxn, byn - Added by Alon

  ! Solving for the new velocity (Stern et al 2017, Eqn B5)
  on_tangential_plane=.false.
  if ((berg%lat>89.) .and. (bergs%grd%grid_is_latlon)) on_tangential_plane=.true.
  if (on_tangential_plane) then
    call rotvec_to_tang(lonn,uvel3,vvel3,xdot3,ydot3)
    call rotvec_to_tang(lonn,ax1,ay1,xddot1,yddot1)
    xdotn=xdot3+(dt*xddot1); ydotn=ydot3+(dt*yddot1)                                    !Alon
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
  else
    uveln=uvel3+(dt*ax1); vveln=vvel3+(dt*ay1)    !Alon , we call it uvel3, vvel3 until it is put into lat/long co-ordinates, where it becomes uveln, vveln
  endif

  !if (berg%id .eq. 1) print *, 'New velocity: ', uveln, vveln


  !!!!!!!!!!!!!!! Debugging  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  error_flag=.false.
  if (.not.error_flag) then
    if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lonn=',lonn,berg%lon
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: latn=',latn,berg%lat
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: u3,un,u0=',uvel3,uveln,berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: v3,vn,v0=',vvel3,vveln,berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: id=',berg%id
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ax1=',&
         & dt*ax1
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ay1=',&
         & dt*ay1
    write(stderrunit,*) 'KID, evolve_iceberg: on_tangential_plane=',on_tangential_plane
    write(stderrunit,*) 'Acceleration terms for position 1'
    error_flag=pos_within_cell(grd, lonn, latn, i, j, xi,  yj)
    call accel(bergs, berg, i, j, xi, yj, latn, uvel3, vvel3, uvel1, vvel1, dt_2, rx, ry, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn, bxn, byn - Added by Alon

    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of cell at end!')
    bounced=is_point_in_cell(bergs%grd, lonn, latn, i, j,explain=.true.)
    if (debug) call error_mesg('KID, evolve_iceberg','berg is out of posn at end!',FATAL)
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
  endif

end subroutine verlet_stepping

!> Calculate explicit and implicit accelerations, new velocity, and new position, using the fourth order Runge-Kutta  method
subroutine Runge_Kutta_stepping(bergs, berg, axn, ayn, bxn, byn, uveln, vveln, lonn, latn, i, j, xi, yj, rx, ry)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer, intent(inout) :: berg !< Iceberg
  real, intent(out) :: axn !< Explicit zonal acceleration (m/s2)
  real, intent(out) :: ayn !< Explicit meridional acceleration (m/s2)
  real, intent(out) :: bxn !< Implicit zonal acceleration (m/s2)
  real, intent(out) :: byn !< Implicit meridional acceleration (m/s2)
  real, intent(out) :: uveln !< New zonal velocity (m/s)
  real, intent(out) :: vveln !< New meridional velocity (m/s)
  real, intent(out) :: lonn !< New longitude (degree E)
  real, intent(out) :: latn !< New latitude (degree N)
  integer, intent(out) :: i !< New i-index of containing cell
  integer, intent(out) :: j !< New i-index of containing cell
  real, intent(out) :: xi !< New non-dimensional x-position
  real, intent(out) :: yj !< New non-dimensional y-position
  real, intent(in) :: rx !< Random number between -1 and 1 for use in x-component of stochastic tidal parameterization
  real, intent(in) :: ry !< Random number between -1 and 1 for use in y-component of stochastic tidal parameterization
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1, ax1, ay1, axn1, ayn1
  real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2, ax2, ay2, axn2, ayn2
  real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3, ax3, ay3, axn3, ayn3
  real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4, ax4, ay4, axn4, ayn4
  real :: x1, xdot1, xddot1, y1, ydot1, yddot1, xddot1n, yddot1n
  real :: x2, xdot2, xddot2, y2, ydot2, yddot2, xddot2n, yddot2n
  real :: x3, xdot3, xddot3, y3, ydot3, yddot3, xddot3n, yddot3n
  real :: x4, xdot4, xddot4, y4, ydot4, yddot4, xddot4n, yddot4n
  real :: xn, xdotn, xddotn, yn, ydotn, yddotn, xddotnn, yddotnn
  real :: dt, dt_2, dt_6, dydl
  integer :: i1,j1,i2,j2,i3,j3,i4,j4
  integer :: stderrunit
  logical :: bounced, on_tangential_plane, error_flag
  ! 4th order Runge-Kutta to solve:
  !    d/dt X = V,  d/dt V = A
  ! with I.C.'s:
  !    X=X1 and V=V1
  !
  !  A1 = A(X1)
  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  !  X4 = X1+  dt*V3 ; V4 = V1+  dt*A3; A4=A(X4)
  !
  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6


  ! Get the stderr unit number
  stderrunit = stderr()
  ! For convenience
  grd=>bergs%grd
  ! Common constants
  dt=bergs%dt
  dt_2=0.5*dt
  dt_6=dt/6.

  i=berg%ine
  j=berg%jne
  xi=berg%xi
  yj=berg%yj
  bounced=.false.
  on_tangential_plane=.false.
  if ((berg%lat>89.) .and. (bergs%grd%grid_is_latlon)) on_tangential_plane=.true.
  i1=i;j1=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness, addfootloose=.true.)

  ! Loading past accelerations - Alon
  axn=berg%axn; ayn=berg%ayn !Alon
  axn1=axn; axn2=axn; axn3=axn; axn4=axn
  ayn1=ayn; ayn2=ayn; ayn3=ayn; ayn4=ayn

  ! A1 = A(X1)
  lon1=berg%lon; lat1=berg%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)

  call  convert_from_meters_to_grid(lat1,bergs%grd%grid_is_latlon ,dxdl1,dydl)
  !dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  !dydl=r180_pi/Rearth
  uvel1=berg%uvel; vvel1=berg%vvel
  if (on_tangential_plane) call rotvec_to_tang(lon1,uvel1,vvel1,xdot1,ydot1)
  u1=uvel1*dxdl1; v1=vvel1*dydl

  call accel(bergs, berg, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, rx, ry, ax1, ay1, axn1, ayn1, bxn, byn) !axn,ayn, bxn, byn  - Added by Alon
  !call accel(bergs, berg, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt, rx, ry, ax1, ay1, axn1, ayn1, bxn, byn) !Note change to dt. Markpoint_1
  if (on_tangential_plane) call rotvec_to_tang(lon1,ax1,ay1,xddot1,yddot1)
  if (on_tangential_plane) call rotvec_to_tang(lon1,axn1,ayn1,xddot1n,yddot1n) !Alon

  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  !if (debug) write(stderr(),*) 'KID, evolve: x2=...'
  if (on_tangential_plane) then
    x2=x1+dt_2*xdot1; y2=y1+dt_2*ydot1
    xdot2=xdot1+dt_2*xddot1; ydot2=ydot1+dt_2*yddot1
    call rotpos_from_tang(x2,y2,lon2,lat2)
    call rotvec_from_tang(lon2,xdot2,ydot2,uvel2,vvel2)
  else
    lon2=lon1+dt_2*u1; lat2=lat1+dt_2*v1
    uvel2=uvel1+dt_2*ax1; vvel2=vvel1+dt_2*ay1
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag, berg%id)
  i2=i; j2=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness, addfootloose=.true.)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon2,lat2,x2,y2)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon2, lat2, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: i1,i2,i=',i1,i2,i
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: j1,j2,j=',j1,j2,j
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lon1,lon2=',lon1,lon2,berg%lon
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lat1,lat2=',lat1,lat2,berg%lat
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: u1,u2,u0=',uvel1,uvel2,berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: v1,v2,v0=',vvel1,vvel2,berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ax1,ax2=',dt*ax1,dt*ax2
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ay1,ay2=',dt*ay1,dt*ay2
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u0=',dt*uvel1,dt*uvel2,dt*berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v0=',dt*vvel1,dt*vvel2,dt*berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2 (deg)=',dt*u1,dt*u2
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2 (deg)=',dt*v1,dt*v2
    write(stderrunit,*) 'Acceleration terms for position 1'
    error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
    call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, rx, ry, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn,- Added by Alon
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 2')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos2 i,j,lon,lat,xi,yj=',i,j,lon2,lat2,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos2 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j,explain=.true.)
    call error_mesg('KID, evolve_iceberg','berg is out of posn at 2!',FATAL)
  endif
  call  convert_from_meters_to_grid(lat2,bergs%grd%grid_is_latlon ,dxdl2,dydl)
  !dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
  u2=uvel2*dxdl2; v2=vvel2*dydl
  call accel(bergs, berg, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, rx, ry, ax2, ay2, axn2, ayn2, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  !call accel(bergs, berg, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt, rx, ry, ax2, ay2, axn2, ayn2, bxn, byn) !Note change to dt. Markpoint_1
  if (on_tangential_plane) call rotvec_to_tang(lon2,ax2,ay2,xddot2,yddot2)
  if (on_tangential_plane) call rotvec_to_tang(lon2,axn2,ayn2,xddot2n,yddot2n) !Alon

  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  !if (debug) write(stderr(),*) 'KID, evolve: x3=...'
  if (on_tangential_plane) then
    x3=x1+dt_2*xdot2; y3=y1+dt_2*ydot2
    xdot3=xdot1+dt_2*xddot2; ydot3=ydot1+dt_2*yddot2
    call rotpos_from_tang(x3,y3,lon3,lat3)
    call rotvec_from_tang(lon3,xdot3,ydot3,uvel3,vvel3)
  else
    lon3=lon1+dt_2*u2; lat3=lat1+dt_2*v2
    uvel3=uvel1+dt_2*ax2; vvel3=vvel1+dt_2*ay2
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag, berg%id)
  i3=i; j3=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness,addfootloose=.true.)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon3,lat3,x3,y3)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon3, lat3, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: i1,i2,i3,i=',i1,i2,i3,i
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: j1,j2,j3,j=',j1,j2,j3,j
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lon1,lon2,lon3=',lon1,lon2,lon3,berg%lon
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lat1,lat2,lat3=',lat1,lat2,lat3,berg%lat
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: u1,u2,u3,u0=',uvel1,uvel2,uvel3,berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: v1,v2,v3,v0=',vvel1,vvel2,vvel3,berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ax1,ax2,ax3=',dt*ax1,dt*ax2,dt*ax3
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ay1,ay2,ay3=',dt*ay1,dt*ay2,dt*ay3
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u3,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v3,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u3 (deg)=',dt*u1,dt*u2,dt*u3
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v3 (deg)=',dt*v1,dt*v2,dt*v3
    write(stderrunit,*) 'Acceleration terms for position 1'
    error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
    call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, rx, ry, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn - Added by Alon
    write(stderrunit,*) 'Acceleration terms for position 2'
    error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
    call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, rx, ry, ax2, ay2, axn2, ayn2, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 3')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos3 i,j,lon,lat,xi,yj=',i,j,lon3,lat3,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos3 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j,explain=.true.)
    call error_mesg('KID, evolve_iceberg','berg is out of posn at 3!',FATAL)
  endif
  call  convert_from_meters_to_grid(lat3,bergs%grd%grid_is_latlon ,dxdl3,dydl)
  !dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
  u3=uvel3*dxdl3; v3=vvel3*dydl
  call accel(bergs, berg, i, j, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, rx, ry, ax3, ay3, axn3, ayn3, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  if (on_tangential_plane) call rotvec_to_tang(lon3,ax3,ay3,xddot3,yddot3)
  if (on_tangential_plane) call rotvec_to_tang(lon3,axn3,ayn3,xddot3n,yddot3n) !Alon

  !  X4 = X1+dt*V3 ; V4 = V1+dt*A3; A4=A(X4)
  !if (debug) write(stderr(),*) 'KID, evolve: x4=...'
  if (on_tangential_plane) then
    x4=x1+dt*xdot3; y4=y1+dt*ydot3
    xdot4=xdot1+dt*xddot3; ydot4=ydot1+dt*yddot3
    call rotpos_from_tang(x4,y4,lon4,lat4)
    call rotvec_from_tang(lon4,xdot4,ydot4,uvel4,vvel4)
  else
    lon4=lon1+dt*u3; lat4=lat1+dt*v3
    uvel4=uvel1+dt*ax3; vvel4=vvel1+dt*ay3
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag, berg%id)
  i4=i; j4=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon4,lat4,x4,y4)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon4, lat4, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lon1,lon2,lon3,lon4=',lon1,lon2,lon3,lon4,berg%lon
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lat1,lat2,lat3,lat4=',lat1,lat2,lat3,lat4,berg%lat
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: u1,u2,u3,u4,u0=',uvel1,uvel2,uvel3,uvel4,berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: v1,v2,v3,v4,v0=',vvel1,vvel2,vvel3,vvel4,berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ax1,ax2,ax3,ax4=',dt*ax1,dt*ax2,dt*ax3,dt*ax4
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ay1,ay2,ay3,ay4=',dt*ay1,dt*ay2,dt*ay3,dt*ay4
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u3,u4,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v3,v4,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u3,u4 (deg)=',dt*u1,dt*u2,dt*u3,dt*u4
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v3,v4 (deg)=',dt*v1,dt*v2,dt*v3,dt*v4
    write(stderrunit,*) 'Acceleration terms for position 1'
    error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
    call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, rx, ry, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
    write(stderrunit,*) 'Acceleration terms for position 2'
    error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
    call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, rx, ry, ax2, ay2, axn2, ayn2, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
    write(stderrunit,*) 'Acceleration terms for position 3'
    error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
    call accel(bergs, berg, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, rx, ry, ax3, ay3, axn3, ayn3, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 4')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos4 i,j,lon,lat,xi,yj=',i,j,lon4,lat4,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos4 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('KID, evolve_iceberg','berg is out of posn at 4!',FATAL)
  endif
  call  convert_from_meters_to_grid(lat4,bergs%grd%grid_is_latlon ,dxdl4,dydl)
  !dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
  u4=uvel4*dxdl4; v4=vvel4*dydl
  call accel(bergs, berg, i, j, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, rx, ry, ax4, ay4, axn4, ayn4, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  if (on_tangential_plane) call rotvec_to_tang(lon4,ax4,ay4,xddot4,yddot4)
  if (on_tangential_plane) call rotvec_to_tang(lon4,axn4,ayn4,xddot4n,yddot4n)

  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
  if (on_tangential_plane) then
    xn=x1+dt_6*( (xdot1+xdot4)+2.*(xdot2+xdot3) )
    yn=y1+dt_6*( (ydot1+ydot4)+2.*(ydot2+ydot3) )
    xdotn=xdot1+dt_6*( (xddot1+xddot4)+2.*(xddot2+xddot3) )
    ydotn=ydot1+dt_6*( (yddot1+yddot4)+2.*(yddot2+yddot3) )
    xddotn=( (xddot1n+xddot4n)+2.*(xddot2n+xddot3n) )/6.  !Alon
    yddotn=( (yddot1n+yddot4n)+2.*(yddot2n+yddot3n) )/6.  !Alon
    call rotpos_from_tang(xn,yn,lonn,latn)
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
    call rotvec_from_tang(lonn,xddotn,yddotn,axn,ayn) !Alon
  else
    lonn=berg%lon+dt_6*( (u1+u4)+2.*(u2+u3) )
    latn=berg%lat+dt_6*( (v1+v4)+2.*(v2+v3) )
    uveln=berg%uvel+dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
    vveln=berg%vvel+dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
    axn=( (axn1+axn4)+2.*(axn2+axn3) )/6. !Alon
    ayn=( (ayn1+ayn4)+2.*(ayn2+ayn3) )/6. !Alon
    bxn=(((ax1+ax4)+2.*(ax2+ax3) )/6)  - (axn/2)
    byn=(((ay1+ay4)+2.*(ay2+ay3) )/6)  - (ayn/2)
  endif

  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag, berg%id)
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness, addfootloose=.true.)

  if (.not.error_flag) then
    if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
    write(stderrunit,'(a,6i5)') 'KID, evolve_iceberg: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lon1,lon2,lon3,lon4,lonn=',lon1,lon2,lon3,lon4,lonn,berg%lon
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: lat1,lat2,lat3,lat4,latn=',lat1,lat2,lat3,lat4,latn,berg%lat
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: u1,u2,u3,u4,un,u0=',uvel1,uvel2,uvel3,uvel4,uveln,berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: v1,v2,v3,v4,vn,v0=',vvel1,vvel2,vvel3,vvel4,vveln,berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ax1,ax2,ax3,ax4,axn=',&
         & dt*ax1,dt*ax2,dt*ax3,dt*ax4,dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* ay1,ay2,ay3,ay4,ayn=',&
         & dt*ay1,dt*ay2,dt*ay3,dt*ay4,dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u3,u4,un,u0=',&
         & dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*uveln,dt*berg%uvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v3,v4,vn,v0=',&
         & dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*vveln,dt*berg%vvel
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* u1,u2,u3,u4,u_rk (deg)=',&
         & dt*u1,dt*u2,dt*u3,dt*u4,dt_6*( (u1+u4)+2.*(u2+u3) )
    write(stderrunit,'(a,6es9.3)') 'KID, evolve_iceberg: dt* v1,v2,v3,v4,v_rk (deg)=',&
         & dt*v1,dt*v2,dt*v3,dt*v4,dt_6*( (v1+v4)+2.*(v2+v3) )
    write(stderrunit,*) 'KID, evolve_iceberg: on_tangential_plane=',on_tangential_plane
    write(stderrunit,*) 'Acceleration terms for position 1'
    error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
    call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, rx, ry, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
    write(stderrunit,*) 'Acceleration terms for position 2'
    error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
    call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, rx, ry, ax2, ay2, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
    write(stderrunit,*) 'Acceleration terms for position 3'
    error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
    call accel(bergs, berg, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, rx, ry, ax3, ay3, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
    write(stderrunit,*) 'Acceleration terms for position 4'
    error_flag=pos_within_cell(grd, lon4, lat4, i4, j4, xi, yj)
    call accel(bergs, berg, i4, j4, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, rx, ry, ax4, ay4, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of cell at end!')
    bounced=is_point_in_cell(bergs%grd, lonn, latn, i, j, explain=.true.)
    if (debug) call error_mesg('KID, evolve_iceberg','berg is out of posn at end!',FATAL)
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
  endif
end subroutine Runge_Kutta_stepping

!> Updates a bergs position using the Verlet algorithm
!!
!! \todo The intent(in) are not consistent with usage, or are not even needed.
subroutine update_verlet_position(bergs, berg)
  type(icebergs), intent(in), pointer :: bergs !< Container for all types and memory
  type(iceberg), intent(in), pointer :: berg !< Iceberg
  !Local variable
  type(icebergs_gridded), pointer :: grd
  real :: lonn, latn
  real :: xi, yj
  real :: uvel3, vvel3
  real :: lon1, lat1, dxdl1, dydl
  real :: uvel1, vvel1, uvel2, vvel2
  real :: axn, ayn, bxn, byn
  real :: xdot2, ydot2
  real :: u2, v2, x1, y1, xn, yn
  real :: dx, dt, dt_2
  integer :: i, j
  logical :: on_tangential_plane, error_flag, bounced
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd
  ! Common constants
  dt=bergs%dt
  dt_2=0.5*dt


  on_tangential_plane=.false.
  if ((berg%lat>89.) .and. (grd%grid_is_latlon)) on_tangential_plane=.true.

  lon1=berg%lon; lat1=berg%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1,berg%id)
  !dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  !dydl=r180_pi/Rearth
  call  convert_from_meters_to_grid(lat1,grd%grid_is_latlon ,dxdl1,dydl)
  uvel1=berg%uvel; vvel1=berg%vvel

  ! Loading past acceleartions - Alon
  axn=berg%axn; ayn=berg%ayn !Alon
  bxn=berg%bxn; byn=berg%byn !Alon

  ! Velocities used to update the position
  uvel2=uvel1+(dt_2*axn)+(dt_2*bxn)                    !Alon
  vvel2=vvel1+(dt_2*ayn)+(dt_2*byn)                    !Alon

  !dx=(dt*(uvel1+(dt_2*axn)+(dt_2*bxn)))

  if (on_tangential_plane) call rotvec_to_tang(lon1,uvel2,vvel2,xdot2,ydot2)
  u2=uvel2*dxdl1; v2=vvel2*dydl

  ! Solving for new position
  if (on_tangential_plane) then
    xn=x1+(dt*xdot2) ; yn=y1+(dt*ydot2)             !Alon
    call rotpos_from_tang(xn,yn,lonn,latn)
  else
    lonn=lon1+(dt*u2) ; latn=lat1+(dt*v2)  !Alon
  endif

  ! Turn the velocities into u_star, v_star.(uvel3 is v_star) - Alon (not sure how this works with tangent plane)
  uvel3=uvel1+(dt_2*axn)                  !Alon
  vvel3=vvel1+(dt_2*ayn)                  !Alon

  ! Adjusting mass...
  !MP3
  i=berg%ine;  j=berg%jne;  xi = berg%xi;  yj = berg%yj
  call adjust_index_and_ground(grd, lonn, latn, uvel3, vvel3, i, j, xi, yj, bounced, error_flag, berg%id)  !Alon:"unclear which velocity to use here?"

  !if (bounced) then
  !  print *, 'you have been bounce: big time!',mpp_pe(),berg%id,lonn, latn, uvel3, vvel3, i, j, xi, yj, bounced, error_flag
  !  berg%axn=0.0  ;  berg%ayn=0.0
  !  berg%bxn=0.0  ;  berg%byn=0.0
  !  berg%uvel=0.0 ;  berg%vvel=0.0
  !endif

  !Updating positions and index
  berg%lon=lonn      ;  berg%lat=latn
  berg%ine=i    ;  berg%jne=j
  berg%xi=xi    ;  berg%yj=yj

end subroutine update_verlet_position

!> Calculate longitude-latitude from tangent plane coordinates
subroutine rotpos_from_tang(x, y, lon, lat)
  ! Arguments
  real, intent(in) :: x !< x-coordinate in tangent plane
  real, intent(in) :: y !< y-coordinate in tangent plane
  real, intent(out) :: lon !< Longitude (degree E)
  real, intent(out) :: lat !< Latitude (degree N)
  ! Local variables
  real :: r

  r=sqrt(x**2+y**2)
  lat=90.-(r180_pi*r/Rearth)
  lon=r180_pi*acos(x/r)*sign(1.,y)

end subroutine rotpos_from_tang

!> Calculates tangent plane velocity from velocity in velocity oriented in geographic coordinates
subroutine rotvec_to_tang(lon, uvel, vvel, xdot, ydot)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: uvel !< Zonal velocity (m/s)
  real, intent(in) :: vvel !< Meridional velocity (m/s)
  real, intent(out) :: xdot !< x-component of velocity in tangent plane (m/s)
  real, intent(out) :: ydot !< y-component of velocity in tangent plane (m/s)
  ! Local variables
  real :: clon,slon

  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  xdot=-slon*uvel-clon*vvel
  ydot=clon*uvel-slon*vvel

end subroutine rotvec_to_tang

!> Calculate velocity oriented in geographic coordinates from tangent plane velocity
subroutine rotvec_from_tang(lon, xdot, ydot, uvel, vvel)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: xdot !< x-component of velocity in tangent plane (m/s)
  real, intent(in) :: ydot !< y-component of velocity in tangent plane (m/s)
  real, intent(out) :: uvel !< Zonal velocity (m/s)
  real, intent(out) :: vvel !< Meridional velocity (m/s)
  ! Local variables
  real :: clon,slon

  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  uvel=-slon*xdot+clon*ydot
  vvel=-clon*xdot-slon*ydot

end subroutine rotvec_from_tang

!> Moves berg's cell indexes,(i,j), checking for collisional with coasts
subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error, id)
  ! Arguments
  type(icebergs_gridded), pointer :: grd !< Container for gridded fields
  real, intent(inout) :: lon !< Longitude (degree E)
  real, intent(inout) :: lat !< Latitude (degree N)
  real, intent(inout) :: uvel !< Zonal velocity (m/s)
  real, intent(inout) :: vvel !< Meridional velocity (m/s)
  real, intent(inout) :: xi !< Non-dimension x-position within cell
  real, intent(inout) :: yj !< Non-dimension y-position within cell
  integer, intent(inout) :: i !< i-index of cell
  integer, intent(inout) :: j !< j-index of cell
  logical, intent(out) :: bounced !< True if berg collided with coast
  logical, intent(out) :: error !< True if adjustments could not be made consistently
  integer(kind=8), intent(in) :: id !< Berg identifier
  ! Local variables
  logical lret, lpos
  real, parameter :: posn_eps=0.05
  integer :: icount, i0, j0, inm, jnm
  real :: xi0, yj0, lon0, lat0
  integer :: stderrunit
  logical :: point_in_cell_using_xi_yj

  ! Get the stderr unit number
  stderrunit = stderr()

  bounced=.false.
  error=.false.
  lon0=lon; lat0=lat ! original position
  i0=i; j0=j ! original i,j
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
!  print *, 'Alon:', lon, lat, i, j, xi, yj, lret
  xi0=xi; yj0=yj ! original xi,yj


  !Removing this while debuggin
  if (debug) then
    !Sanity check lret, xi and yj
    lret=is_point_in_cell(grd, lon, lat, i, j)
    point_in_cell_using_xi_yj=is_point_within_xi_yj_bounds(xi,yj)
    if (.not. point_in_cell_using_xi_yj) then

      if (lret) then
        write(stderrunit,*) 'KID, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'KID, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'KID, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'KID, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'KID, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'KID, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'KID, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,explain=.true.)
        write(stderrunit,*) 'KID, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
        write(stderrunit,*) 'KID, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','Iceberg is_point_in_cell=True but xi, yi are out of cell',FATAL)
        error=.true.; return
     endif
    else
      if (.not.lret) then
        write(stderrunit,*) 'KID, adjust: WARNING!!! lret=F but |xi,yj|<1',mpp_pe()
        write(stderrunit,*) 'KID, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'KID, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'KID, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'KID, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'KID, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'KID, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,  explain=.true.)
        write(stderrunit,*) 'KID, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'KID, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','Iceberg is_point_in_cell=False but xi, yi are out of cell',FATAL)
        error=.true.; return
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  endif ! debug

  if (lret) return ! Berg was already in cell

  ! Find inm, jnm (as if adjusting i,j) based on xi,yj
  ! ignoring the mand mask.
  ! NOTE:  This search appears to have *NO* active role
  ! in the algorithm other than to flag a warning.
  icount=0
  inm=i0; jnm=j0 ! original i,j
  do while (debug .and. .not.lret .and. icount<4)
    icount=icount+1
    if (xi.lt.0.) then
      if (inm>grd%isd) then
        inm=inm-1
      endif
    elseif (xi.gt.1.) then
!    elseif (xi.ge.1.) then   !Alon: maybe it should be .ge.
      if (inm<grd%ied) then
        inm=inm+1
      endif
    endif
    if (yj.lt.0.) then
      if (jnm>grd%jsd) then
        jnm=jnm-1
      endif
    elseif (yj.gt.1.) then
!    elseif (yj.ge.1.) then   !Alon:maybe it should be .ge.
      if (jnm<grd%jed) then
        jnm=jnm+1
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj) ! Update xi and yj
  enddo
  if (abs(inm-i0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'KID, adjust: inm,i0,inm-i0=',inm,i0,inm-i0
   !stop 'Moved too far in i without mask!'
  endif
  if (abs(jnm-j0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'KID, adjust: jnm,i0,jnm-j0=',jnm,j0,inm-j0
   !stop 'Moved too far in j without mask!'
  endif

  ! Adjust i,j based on xi,yj while bouncing off of masked land cells
  icount=0
  lret=pos_within_cell(grd, lon, lat, i0, j0, xi, yj)
  do while ( .not.lret.and. icount<4 )
    icount=icount+1
    if (xi.lt.0.) then
      if (i>grd%isd) then
        if (grd%msk(i-1,j)>0.) then
          if (i>grd%isd+1) i=i-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'KID, adjust: bouncing berg from west',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (xi.ge.1.) then    !Alon!!!!
!    elseif (xi.gt.1.) then
      if (i<grd%ied) then
        if (grd%msk(i+1,j)>0.) then
          if (i<grd%ied) i=i+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'KID, adjust: bouncing berg from east',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (yj.lt.0.) then
      if (j>grd%jsd) then
        if (grd%msk(i,j-1)>0.) then
          if (j>grd%jsd+1) j=j-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'KID, adjust: bouncing berg from south',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (yj.ge.1.) then     !Alon.
!    elseif (yj.gt.1.) then
      if (j<grd%jed) then
        if (grd%msk(i,j+1)>0.) then
          if (j<grd%jed) j=j+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'KID, adjust: bouncing berg from north',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (bounced) then
      if (xi>=1.) xi=1.-posn_eps   !Alon.
!      if (xi>1.) xi=1.-posn_eps   !
      if (xi<0.) xi=posn_eps
      if (yj>=1.) yj=1.-posn_eps  !Alon.
!      if (yj>1.) yj=1.-posn_eps
      if (yj<0.) yj=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
    endif
    if (debug) then
      if (grd%msk(i,j)==0.) stop 'KID, adjust: Berg is in land! This should not happen...'
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  enddo
 !if (debug) then
 !  if (abs(i-i0)>2) then
 !    stop 'KID, adjust: Moved too far in i!'
 !  endif
 !  if (abs(j-j0)>2) then
 !    stop 'KID, adjust: Moved too far in j!'
 !  endif
 !endif

  if (.not.bounced.and.lret.and.grd%msk(i,j)>0.) return ! Landed in ocean without bouncing so all is well

  if (.not.bounced.and..not.lret) then ! This implies the berg traveled many cells without getting far enough
    if (debug) then
      write(stderrunit,*) 'KID, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'KID, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'KID, adjust: i0,j0=',i0,j0
      write(stderrunit,*) 'KID, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'KID, adjust: xi,yj=',xi,yj
      write(stderrunit,*) 'KID, adjust: i,j=',i,j
      write(stderrunit,*) 'KID, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'KID, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
      write(stderrunit,*) 'KID, adjust: lret=',lret
    endif

    if (abs(i-i0)+abs(j-j0)==0) then
      if (use_roundoff_fix) then
        ! This is a special case due to round off where is_point_in_cell()
        ! returns false but xi and yj are between 0 and 1.
        ! It occurs very rarely but often enough to have brought down
        ! ESM2G four times since the spin-up began. (as of 8/10/2010)
        ! This temporary fix arbitrarily moves the berg toward the
        ! center of the current cell.
        xi=(xi-0.5)*(1.-posn_eps)+0.5
        yj=(yj-0.5)*(1.-posn_eps)+0.5
      endif
      call error_mesg('KID, adjust', 'Berg did not move or bounce during iterations AND was not in cell. Adjusting!', WARNING)
      write(stderrunit,*) 'KID, adjust: The adjusting iceberg is: ', id,  mpp_pe()
      write(stderrunit,*) 'KID, adjust: The adjusting lon,lat,u,v: ', lon, lat, uvel, vvel
      write(stderrunit,*) 'KID, adjust: The adjusting xi,ji: ', xi, yj
      lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj,explain=.true.)
    else
      call error_mesg('KID, adjust', 'Berg iterated many times without bouncing!', WARNING)
    endif
  endif
!  if (xi>1.) xi=1.-posn_eps    !Alon
  if (xi>=1.) xi=1.-posn_eps
  if (xi<0.) xi=posn_eps
  if (yj>1.) yj=1.-posn_eps
!  if (yj>1.) yj=1.-posn_eps
  if (yj<=0.) yj=posn_eps        !Alon
  lon=bilin(grd, grd%lon, i, j, xi, yj)
  lat=bilin(grd, grd%lat, i, j, xi, yj)
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  if (.not. lret) then
    write(0,*) 'i0, j0,=', i0,j0
    write(0,*) 'xi0, yj0,=', xi0,yj0
    write(0,*) 'grd%msk(i0, j0)=', grd%msk(i0,j0)
    write(0,*) 'lon0, lat0,=', lon0,lat0
    write(0,*) 'i,j,lon, lat,grd%msk(i,j)=', i,j,lon,lat,grd%msk(i,j)
    write(stderrunit,*) 'KID, adjust: Should not get here! Berg is not in cell after adjustment', id, mpp_pe()
    if (debug) error=.true.
  endif
 end subroutine adjust_index_and_ground

!> Calculate longitude-latitude from tangent plane coordinates
subroutine rotpos_to_tang(lon, lat, x, y, id_in)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: lat !< Latitude (degree N)
  real, intent(out) :: x !< x-coordinate in tangent plane
  real, intent(out) :: y !< y-coordinate in tangent plane
  integer(kind=8), intent(in), optional :: id_in !< Berg identifier
  ! Local variables
  real :: r,colat,clon,slon
  integer :: stderrunit, id

  stderrunit = stderr()
  id=0
  if (present(id_in)) then
        id=id_in
  endif

  if (lat>90.) then
      write(stderrunit,*) 'KID, rotpos_to_tang: lat>90 already!',lat, lon, id
      call error_mesg('KID, rotpos_to_tang','Something went very wrong!',FATAL)
  endif
  if (lat==90.) then
    write(stderrunit,*) 'KID, rotpos_to_tang: lat==90 already!',lat, lon
    call error_mesg('KID, rotpos_to_tang','Something went wrong!',FATAL)
  endif

  colat=90.-lat
  r=Rearth*(colat*pi_180)
  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  x=r*clon
  y=r*slon

end subroutine rotpos_to_tang

!> Calculate stocks of water and heat
subroutine icebergs_stock_pe(bergs, index, value)
  ! Modules
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer, intent(in) :: index !< =ISTOCK_WATER or ISTOCK_HEAT
  real, intent(out) :: value !< Amount of ice or water
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: berg_mass, stored_mass

  ! For convenience
  grd=>bergs%grd

  select case (index)

  case (ISTOCK_WATER)
    berg_mass=sum_mass(bergs)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=stored_mass+berg_mass

  case (ISTOCK_HEAT)
    berg_mass=sum_mass(bergs)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=-(stored_mass+berg_mass)*HLF ! HLF is in (J/kg) from constants_mod

  case default
    value = 0.0

  end select

end subroutine icebergs_stock_pe

!> Write restart files
subroutine icebergs_save_restart(bergs, time_stamp)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  character(len=*),    intent(in), optional :: time_stamp !< Timestamp for restart file
  ! Local variables

  if (.not.associated(bergs)) return

  call mpp_clock_begin(bergs%clock_iow)
  call bergs_chksum(bergs, 'write_restart bergs')
  call write_restart(bergs, time_stamp)
  call mpp_clock_end(bergs%clock_iow)

end subroutine icebergs_save_restart

!> Deallocate all memory and disassociated pointer
subroutine icebergs_end(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(iceberg), pointer :: this, next

  if (.not.associated(bergs)) return

  ! icebergs_save_restart() is called directly by SIS1 and SIS2 so
  ! we do not need to call it a second time. If icebergs were controlled
  ! by the coupler then the icebergs would need to take responsibility for
  ! the restarts at the end of the run.
  !call icebergs_save_restart(bergs)

  call mpp_clock_begin(bergs%clock_ini)
  ! Delete bergs and structures
  call move_all_trajectories(bergs, delete_bergs=.true.)

  if (.not. bergs%ignore_traj) then
    call write_trajectory(bergs%trajectories, bergs%save_short_traj, bergs%save_fl_traj, bergs%fl_r)
    if (save_bond_traj) call write_bond_trajectory(bergs%bond_trajectories)
  endif

  deallocate(bergs%grd%lon)
  deallocate(bergs%grd%lat)
  deallocate(bergs%grd%lonc)
  deallocate(bergs%grd%latc)
  deallocate(bergs%grd%dx)
  deallocate(bergs%grd%dy)
  deallocate(bergs%grd%area)
  deallocate(bergs%grd%msk)
  deallocate(bergs%grd%cos)
  deallocate(bergs%grd%sin)
  deallocate(bergs%grd%ocean_depth)
  deallocate(bergs%grd%calving)
  deallocate(bergs%grd%calving_hflx)
  deallocate(bergs%grd%stored_heat)
  deallocate(bergs%grd%floating_melt)
  deallocate(bergs%grd%berg_melt)
  deallocate(bergs%grd%melt_buoy)
  deallocate(bergs%grd%melt_eros)
  deallocate(bergs%grd%melt_conv)
  deallocate(bergs%grd%bergy_src)
  deallocate(bergs%grd%bergy_melt)
  deallocate(bergs%grd%bergy_mass)
  deallocate(bergs%grd%fl_bits_src)
  deallocate(bergs%grd%fl_bits_melt)
  deallocate(bergs%grd%fl_bits_mass)
  deallocate(bergs%grd%fl_bergy_bits_mass)
  deallocate(bergs%grd%spread_mass)
  deallocate(bergs%grd%spread_mass_old)
  deallocate(bergs%grd%spread_area)
  deallocate(bergs%grd%virtual_area)
  deallocate(bergs%grd%mass)
  deallocate(bergs%grd%mass_on_ocean)
  deallocate(bergs%grd%area_on_ocean)
  deallocate(bergs%grd%tmp)
  deallocate(bergs%grd%tmpc)
  deallocate(bergs%grd%stored_ice)
  deallocate(bergs%grd%rmean_calving)
  deallocate(bergs%grd%rmean_calving_hflx)
  deallocate(bergs%grd%real_calving)
  deallocate(bergs%grd%uo)
  deallocate(bergs%grd%vo)
  deallocate(bergs%grd%ui)
  deallocate(bergs%grd%vi)
  deallocate(bergs%grd%ua)
  deallocate(bergs%grd%va)
  deallocate(bergs%grd%ssh)
  deallocate(bergs%grd%sst)
  deallocate(bergs%grd%sss)
  deallocate(bergs%grd%cn)
  deallocate(bergs%grd%hi)
  deallocate(bergs%grd%melt_by_class)
  deallocate(bergs%grd%melt_buoy_fl)
  deallocate(bergs%grd%melt_eros_fl)
  deallocate(bergs%grd%melt_conv_fl)
  deallocate(bergs%grd%fl_parent_melt)
  deallocate(bergs%grd%fl_child_melt)
  deallocate(bergs%grd%domain)
  deallocate(bergs%grd)
  deallocate(bergs%initial_mass_s)
  deallocate(bergs%distribution_s)
  deallocate(bergs%initial_thickness_s)
  deallocate(bergs%initial_width_s)
  deallocate(bergs%initial_length_s)
  deallocate(bergs%initial_mass_n)
  deallocate(bergs%distribution_n)
  deallocate(bergs%initial_thickness_n)
  deallocate(bergs%initial_width_n)
  deallocate(bergs%initial_length_n)
  call dealloc_buffer(bergs%obuffer_n)
  call dealloc_buffer(bergs%obuffer_s)
  call dealloc_buffer(bergs%obuffer_e)
  call dealloc_buffer(bergs%obuffer_w)
  call dealloc_buffer(bergs%ibuffer_n)
  call dealloc_buffer(bergs%ibuffer_s)
  call dealloc_buffer(bergs%ibuffer_e)
  call dealloc_buffer(bergs%ibuffer_w)
  call dealloc_buffer(bergs%ibuffer_io)
  call dealloc_buffer(bergs%ibuffer_io)
  call mpp_clock_end(bergs%clock_ini)
  deallocate(bergs)

  if (mpp_pe()==mpp_root_pe()) write(*,'(a,i8)') 'KID: icebergs_end complete',mpp_pe()

end subroutine icebergs_end

!> deallocate a buffer
subroutine dealloc_buffer(buff)
  ! Arguments
  type(buffer), pointer :: buff !< the buffer
  ! Local variables
    if (associated(buff)) then
      if (associated(buff%data)) deallocate(buff%data)
      deallocate(buff)
    endif
end subroutine dealloc_buffer

!> Approximately convert a wind-stress into a velocity difference
subroutine invert_tau_for_du(u, v)
  ! Arguments
  real, dimension(:,:), intent(inout) :: u !< On entry, zonal wind stress (Pa). On exit, zonal velocity difference (m/s).
  real, dimension(:,:), intent(inout) :: v !< On entry, meridional wind stress (Pa). On exit, meridional velocity difference (m/s).
  ! Local variables
  integer :: i, j
  real :: cd, cddvmod, tau2

  cd=0.0015

  do j=lbound(u,2), ubound(u,2)
    do i=lbound(u,1), ubound(u,1)
      tau2=u(i,j)*u(i,j)+v(i,j)*v(i,j)
      cddvmod=sqrt(cd*sqrt(tau2))
      if (cddvmod.ne.0.) then
        u(i,j)=u(i,j)/cddvmod
        v(i,j)=v(i,j)/cddvmod
      else
        u(i,j)=0.
        v(i,j)=0.
      endif
    enddo
  enddo

end subroutine invert_tau_for_du

end module
