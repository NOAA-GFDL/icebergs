module ice_bergs

use constants_mod, only: pi, omega, HLF
use fms_mod, only: open_namelist_file, check_nml_error, close_file
use fms_mod, only: field_exist, get_global_att_value
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: write_version_number, read_data, write_data, file_exist
use mosaic_mod, only: get_mosaic_ntiles, get_mosaic_ncontacts
use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_chksum
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP

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
use ice_bergs_framework, only: find_cell,find_cell_by_search,count_bergs,is_point_in_cell,pos_within_cell
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

use ice_bergs_io,        only: ice_bergs_io_init,write_restart,write_trajectory
use ice_bergs_io,        only: read_restart_bergs,read_restart_bergs_orig,read_restart_calving
use ice_bergs_io,        only: read_restart_bonds
use ice_bergs_io,        only: read_ocean_depth

implicit none ; private

public icebergs_init, icebergs_end, icebergs_run, icebergs_stock_pe, icebergs
public icebergs_incr_mass, icebergs_save_restart

real, parameter :: pi_180=pi/180.  ! Converts degrees to radians
real, parameter :: r180_pi=180./pi ! Converts radians to degrees
real, parameter :: Rearth=6360000. ! Radius of earth (m)
real, parameter :: rho_ice=916.7 ! Density of fresh ice @ 0oC (kg/m^3)
real, parameter :: rho_water=999.8 ! Density of fresh water @ 0oC (kg/m^3)
real, parameter :: rho_air=1.1 ! Density of air @ 0oC (kg/m^3) ???
real, parameter :: rho_seawater=1025. ! Approx. density of surface sea water @ 0oC (kg/m^3)
real, parameter :: gravity=9.8 ! Gravitational acceleratio (m/s^2)
real, parameter :: Cd_av=1.3 ! (Vertical) Drag coefficient between bergs and atmos (?)
real, parameter :: Cd_ah=0.0055 ! (Horizontal) Drag coefficient between bergs and atmos (?)
real, parameter :: Cd_wv=0.9 ! (Vertical) Drag coefficient between bergs and ocean (?)
real, parameter :: Cd_wh=0.0012 ! (Horizontal) Drag coefficient between bergs and ocean (?)
real, parameter :: Cd_iv=0.9 ! (Vertical) Drag coefficient between bergs and sea-ice (?)
!TOM> no horizontal drag for sea ice! real, parameter :: Cd_ih=0.0012 ! (Horizontal) Drag coefficient between bergs and sea-ice (?)

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

! ##############################################################################
subroutine icebergs_init(bergs, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth, maskmap, fractional_area)
! Arguments
type(icebergs), pointer :: bergs
type(icebergs_gridded), pointer :: grd => null()
integer, intent(in) :: gni, gnj, layout(2), io_layout(2), axes(2)
integer, intent(in) :: dom_x_flags, dom_y_flags
real, intent(in) :: dt
type (time_type), intent(in) :: Time ! current time
real, dimension(:,:), intent(in) :: ice_lon, ice_lat, ice_wet
real, dimension(:,:), intent(in) :: ice_dx, ice_dy, ice_area
real, dimension(:,:), intent(in) :: cos_rot, sin_rot
real, dimension(:,:), intent(in), optional :: ocean_depth
logical, intent(in), optional :: maskmap(:,:)
logical, intent(in), optional :: fractional_area
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
  call read_restart_calving(bergs)  !This is moved to before restart_bergs (by Alon) so that generate icebergs can have the correct counter
  if(orig_read) then
     call read_restart_bergs_orig(bergs,Time)
  else
     call read_restart_bergs(bergs,Time)
  endif
  call bergs_chksum(bergs, 'read_restart bergs')
  if (fix_restart_dates) call offset_berg_dates(bergs,Time)
  !call read_restart_calving(bergs)
  call mpp_clock_end(bergs%clock_ior)

  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_init, initial status')
  
  !Reading ocean depth from a file
  if (bergs%read_ocean_depth_from_file) call read_ocean_depth(bergs%grd)

  if (bergs%iceberg_bonds_on) then
    call update_halo_icebergs(bergs)
    if (bergs%manually_initialize_bonds) then
      call initialize_iceberg_bonds(bergs)
    else
      call read_restart_bonds(bergs,Time)
    endif
    call update_halo_icebergs(bergs)
    call connect_all_bonds(bergs)
    nbonds=0
    check_bond_quality=.True.
    call count_bonds(bergs, nbonds,check_bond_quality)
  endif

end subroutine icebergs_init


! ##############################################################################
subroutine unit_testing(bergs)
! Arguments
type(icebergs), pointer :: bergs

call hexagon_test()
call point_in_triangle_test()
call basal_melt_test(bergs)
call test_check_for_duplicate_ids_in_list()

end subroutine unit_testing

subroutine basal_melt_test(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs
  real :: dvo,lat,salt,temp, basal_melt, thickness
  integer :: iceberg_num
  logical :: Use_three_equation_model

  if (mpp_pe() .eq. mpp_root_pe() ) print *, 'Begining Basal Melting Unit Test'
  dvo=0.2 ;lat=0.0 ; salt=35.0 ; temp=2.0 ;thickness=100.; iceberg_num=0
  Use_three_equation_model=.False.
  call find_basal_melt(bergs,dvo,lat,salt,temp,Use_three_equation_model,thickness,basal_melt,iceberg_num)
  if (mpp_pe() .eq. mpp_root_pe()) print *, 'Two equation model basal_melt =',basal_melt

  Use_three_equation_model=.True.
  call find_basal_melt(bergs,dvo,lat,salt,temp,Use_three_equation_model,thickness,basal_melt,iceberg_num)
   if (mpp_pe() .eq. mpp_root_pe()) print *, 'Three equation model basal_melt =',basal_melt

end subroutine basal_melt_test

subroutine point_in_triangle_test()
! Arguments
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
  if (fail_unit_test) call error_mesg('diamonds, hexagon unit testing:', 'Point in triangle test does not pass!', FATAL)

end subroutine point_in_triangle_test

subroutine hexagon_test()
! Arguments
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
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon at origin has the wrong area!', WARNING)
    if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    fail_unit_test=.True.
  endif
  if (((abs((Area_hex/4)-Area_Q1 )>tol) .or.  (abs((Area_hex/4)-Area_Q2 )>tol)) .or. ((abs((Area_hex/4)-Area_Q3 )>tol) .or. (abs((Area_hex/4)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon at origin divides into unqual parts!', WARNING)
    fail_unit_test=.True.
  endif

  ! Test 2:  Hexagon split into two quadrants
  !Test 2a: center on x>0 axis
  x0=S  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q1 )>tol) .or.  (abs(0-Area_Q2 )>tol)) .or. ((abs(0-Area_Q3 )>tol) .or. (abs((Area_hex/2)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split btw 1 and 4!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 2b: center on x<0 axis
  x0=-S  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q2 )>tol) .or.  (abs(0-Area_Q1 )>tol)) .or. ((abs(0-Area_Q4 )>tol) .or. (abs((Area_hex/2)-Area_Q3 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split btw 2 and 3!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 2c: center on y>0 axis
  x0=0.  ;  y0=H
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q1 )>tol) .or.  (abs(0-Area_Q3 )>tol)) .or. ((abs(0-Area_Q4 )>tol) .or. (abs((Area_hex/2)-Area_Q2 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split btw 1 and 2!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 3d: center on y<0 axis
  x0=0.  ;  y0=-H
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((Area_hex/2)-Area_Q3 )>tol) .or.  (abs(0-Area_Q1 )>tol)) .or. ((abs(0-Area_Q2 )>tol) .or. (abs((Area_hex/2)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split btw 3 and 4!', WARNING)
    fail_unit_test=.True.
  endif
  
  ! Test 3:  Two corners of hex on the axis
  !Test 3a: center on x>0 axis
  x0=S/2.  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((2.5*Area_hex/6.)-Area_Q1 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q2 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q3 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q4 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split two coners of hex (x>0)!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 3b: center on x<0 axis
  x0=-S/2.  ;  y0=0.
  call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  if (((abs((2.5*Area_hex/6.)-Area_Q2 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q1 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q4 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q3 )>tol))) then
  if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
    call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split two coners of hex (x<0)!', WARNING)
    fail_unit_test=.True.
  endif
  !Test 3c: center on y>0 axis
  !x0=0.  ;  y0=H/2.
  !call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  !if (((abs((2.5*Area_hex/6.)-Area_Q1 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q3 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q4 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q2 )>tol))) then
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon errors =', (abs((2.5*Area_hex/6.)-Area_Q1 )), (abs((0.5*Area_hex/6.)-Area_Q3 )),&
  !  call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split two coners of hex (y>0)!', WARNING)
  !  fail_unit_test=.True.
  !endif
  !!Test 3d: center on y<0 axis
  !x0=0.  ;  y0=-H/2.
  !call Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  !if (((abs((2.5*Area_hex/6.)-Area_Q3 )>tol) .or.  (abs((0.5*Area_hex/6.)-Area_Q2 )>tol)) .or. ((abs((0.5*Area_hex/6.)-Area_Q1 )>tol) .or. (abs((2.5*Area_hex/6.)-Area_Q4 )>tol))) then
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon areas =', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
  !if (mpp_pe() .eq. mpp_root_pe()) write(stderrunit,*) 'diamonds, hexagon errots =', (abs((2.5*Area_hex/6.)-Area_Q3 )), (abs((0.5*Area_hex/6.)-Area_Q2 )),&
  !  call error_mesg('diamonds, hexagon unit testing:', 'Hexagon split two coners of hex (y<0)!', WARNING)
  !  fail_unit_test=.True.
  !endif


  if (fail_unit_test) call error_mesg('diamonds, hexagon unit testing:', 'Hexagon unit testing does not pass!', FATAL)


end subroutine hexagon_test

! ##############################################################################

subroutine initialize_iceberg_bonds(bergs)

type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
type(iceberg), pointer :: other_berg
type(icebergs_gridded), pointer :: grd
real :: T1, L1, W1, lon1, lat1, x1, y1, R1, A1   !Current iceberg
real :: T2, L2, W2, lon2, lat2, x2, y2, R2, A2   !Other iceberg
real :: dlon,dlat
real :: dx_dlon,dy_dlat, lat_ref
real :: r_dist_x, r_dist_y, r_dist
integer :: grdi_outer, grdj_outer
integer :: grdi_inner, grdj_inner


  ! For convenience
    grd=>bergs%grd
    !Should update halos before doing this 
  do grdj_outer = grd%jsc,grd%jec ; do grdi_outer = grd%isc,grd%iec  !Should you be on the data domain??
    berg=>bergs%list(grdi_outer,grdj_outer)%first
    do while (associated(berg)) ! loop over all bergs
    
      lon1=berg%lon; lat1=berg%lat
      !call rotpos_to_tang(lon1,lat1,x1,y1)  !Is this correct? Shouldn't it only be on tangent plane?

      do grdj_inner = grd%jsc,grd%jec ; do grdi_inner = grd%isc,grd%iec  !This line uses n^2 steps
!     do grdj_inner = berg%jne-1,berg%jne+1 ; do grdi_inner = berg%ine-1,berg%ine+1   !Only looping through adjacent cells.
        other_berg=>bergs%list(grdi_inner,grdj_inner)%first
        do while (associated(other_berg)) ! loop over all other bergs  
          
          if (berg%iceberg_num .ne. other_berg%iceberg_num) then
            lon2=other_berg%lon; lat2=other_berg%lat
            !call rotpos_to_tang(lon2,lat2,x2,y2) !Is this correct? Shouldn't it only be on tangent plane?
            !r_dist_x=x1-x2 ; r_dist_y=y1-y2
            !r_dist=sqrt( ((x1-x2)**2) + ((y1-y2)**2) )

            dlon=lon1-lon2
            dlat=lat1-lat2
            lat_ref=0.5*(lat1+lat2)
            call convert_from_grid_to_meters(lat_ref,grd%grid_is_latlon,dx_dlon,dy_dlat)
            r_dist_x=dlon*dx_dlon
            r_dist_y=dlat*dy_dlat
            r_dist=sqrt( (r_dist_x**2) + (r_dist_y**2) )
        
            !if (r_dist.gt.1000.) then  ! If the bergs are close together, then form a bond
              call form_a_bond(berg, other_berg%iceberg_num, other_berg%ine, other_berg%jne, other_berg)
            !endif       
          endif
          other_berg=>other_berg%next
        enddo  ! End of looping through all other bergs in the inner list
      enddo ; enddo;  !End of inner loop
      berg=>berg%next
    enddo ! End of looping through all bergs in the outer list
  enddo ; enddo; !End of outer loop.


end subroutine initialize_iceberg_bonds

subroutine  convert_from_grid_to_meters(lat_ref,grid_is_latlon ,dx_dlon,dy_dlat)
  ! Arguments
  real, intent(in) :: lat_ref
  logical, intent(in) :: grid_is_latlon
  real, intent(out) :: dx_dlon,dy_dlat

  if (grid_is_latlon) then
    dx_dlon=(pi/180.)*Rearth*cos((lat_ref)*(pi/180.))
    dy_dlat=(pi/180.)*Rearth

  else
    dx_dlon=1.
    dy_dlat=1.

  endif
end subroutine  convert_from_grid_to_meters

subroutine  convert_from_meters_to_grid(lat_ref,grid_is_latlon ,dlon_dx,dlat_dy)
  ! Arguments
  real, intent(in) :: lat_ref
  logical, intent(in) :: grid_is_latlon
  real, intent(out) :: dlon_dx,dlat_dy

  if (grid_is_latlon) then
    dlon_dx=(180./pi)/(Rearth*cos((lat_ref)*(pi/180.)))
    dlat_dy=(180./pi)/Rearth

  else
    dlon_dx=1.
    dlat_dy=1.

  endif
end subroutine  convert_from_meters_to_grid
! ##############################################################################

subroutine interactive_force(bergs,berg,IA_x, IA_y, u0, v0, u1, v1, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) !Calculating interactive force between icebergs. Alon,  Markpoint_4
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
type(iceberg), pointer :: other_berg
type(bond), pointer :: current_bond
real, intent(in) :: u0,v0, u1, v1
real :: u2, v2
logical :: critical_interaction_damping_on
real, intent(out) :: IA_x, IA_y 
real, intent(out) :: P_ia_11, P_ia_12, P_ia_22, P_ia_21, P_ia_times_u_x, P_ia_times_u_y
integer :: grdi, grdj
logical :: iceberg_bonds_on
logical :: bonded
iceberg_bonds_on=bergs%iceberg_bonds_on

  IA_x=0.
  IA_y=0.
  P_ia_11=0. ; P_ia_12=0. ;  P_ia_21=0.;  P_ia_22=0. 
  P_ia_times_u_x=0. ; P_ia_times_u_y=0.

  bonded=.false. !Unbonded iceberg interactions
  do grdj = berg%jne-1,berg%jne+1 ; do grdi = berg%ine-1,berg%ine+1  !Note: need  to make sure this is wide enough, but less than the halo width
    other_berg=>bergs%list(grdi,grdj)%first
    do while (associated(other_berg)) ! loop over all other bergs  
      call calculate_force(bergs,berg,other_berg,IA_x, IA_y, u0, v0, u1, v1, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y, bonded) 
      other_berg=>other_berg%next
    enddo ! loop over all bergs
  enddo ; enddo

  bonded=.true.  !Interactions due to iceberg bonds
  if (iceberg_bonds_on) then ! MP1  
    current_bond=>berg%first_bond
    do while (associated(current_bond)) ! loop over all bonds
      other_berg=>current_bond%other_berg
      if (.not. associated(other_berg)) then
        call error_mesg('diamonds,bond interactions', 'Trying to do Bond interactions with unassosiated berg!' ,FATAL)
      else
        call calculate_force(bergs,berg,other_berg,IA_x, IA_y, u0, v0, u1, v1, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded) 
      endif
      current_bond=>current_bond%next_bond
    enddo
  endif

  !print *,'IA_x=',IA_x,'IA_y',IA_y, berg%iceberg_num
  !print *,'P_ia_11',P_ia_11,'P_ia_12',P_ia_12, 'P_ia_21',P_ia_21,'P_ia_22', P_ia_22
  !print *, 'P_ia_times_u_x', P_ia_times_u_x, 'P_ia_times_u_y', P_ia_times_u_y
  contains


   subroutine calculate_force(bergs,berg,other_berg,IA_x, IA_y, u0, v0, u1, v1, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y,bonded) 
   !Arguments
     type(icebergs), pointer :: bergs
     type(iceberg), pointer :: berg
     type(iceberg), pointer :: other_berg
     real :: T1, L1, W1, lon1, lat1, x1, y1, R1, A1   !Current iceberg
     real :: T2, L2, W2, lon2, lat2, x2, y2, R2, A2   !Other iceberg
     real :: dlon, dlat
     real :: r_dist_x, r_dist_y, r_dist, A_o, A_min, trapped, T_min
     real, intent(in) :: u0,v0, u1, v1
     real :: P_11, P_12, P_21, P_22
     real :: M1, M2, M_min
     real :: u2, v2
     real :: lat_ref, dx_dlon, dy_dlat
     logical :: critical_interaction_damping_on
     real :: spring_coef, accel_spring, radial_damping_coef, p_ia_coef, tangental_damping_coef, bond_coef
     real, intent(inout) :: IA_x, IA_y 
     real, intent(inout) :: P_ia_11, P_ia_12, P_ia_22, P_ia_21, P_ia_times_u_x, P_ia_times_u_y
     logical ,intent(in) :: bonded

      spring_coef=bergs%spring_coef
      !bond_coef=bergs%bond_coef 
      radial_damping_coef=bergs%radial_damping_coef
      tangental_damping_coef=bergs%tangental_damping_coef
      critical_interaction_damping_on=bergs%critical_interaction_damping_on

      !Using critical values for damping rather than manually setting the damping.
      if (critical_interaction_damping_on) then
        radial_damping_coef=2.*sqrt(spring_coef)  ! Critical damping  
        tangental_damping_coef=(2.*sqrt(spring_coef))/4  ! Critical damping   (just a guess)
      endif

      if (berg%iceberg_num .ne. other_berg%iceberg_num) then
        !From Berg 1
        L1=berg%length
        W1=berg%width
        T1=berg%thickness 
        M1=berg%mass 
        A1=L1*W1 
        lon1=berg%lon_old; lat1=berg%lat_old
        !call rotpos_to_tang(lon1,lat1,x1,y1)

        !From Berg 1
        L2=other_berg%length
        W2=other_berg%width
        T2=other_berg%thickness 
        M2=other_berg%mass 
        u2=other_berg%uvel_old !Old values are used to make it order invariant 
        v2=other_berg%vvel_old !Old values are used to make it order invariant 
        A2=L2*W2
        lon2=other_berg%lon_old; lat2=other_berg%lat_old !Old values are used to make it order invariant

        !call rotpos_to_tang(lon2,lat2,x2,y2)

        dlon=lon1-lon2
        dlat=lat1-lat2
        
        !Note that this is not the exact distance along a great circle.
        !Approximation for small distances. Should be fine.
        !r_dist_x=x1-x2 ; r_dist_y=y1-y2
        !r_dist=sqrt( ((x1-x2)**2) + ((y1-y2)**2) )
        lat_ref=0.5*(lat1+lat2)
        call convert_from_grid_to_meters(lat_ref,bergs%grd%grid_is_latlon,dx_dlon,dy_dlat)

        r_dist_x=dlon*dx_dlon
        r_dist_y=dlat*dy_dlat
        r_dist=sqrt( (r_dist_x**2) + (r_dist_y**2) )
       
        if (bergs%hexagonal_icebergs) then 
          R1=sqrt(A1/(2.*sqrt(3.)))
          R2=sqrt(A2/(2.*sqrt(3.)))
        else !square packing
          R1=sqrt(A1/pi) ! Interaction radius of the iceberg (assuming circular icebergs)
          R2=sqrt(A2/pi) ! Interaction radius of the other iceberg
        endif
       !!!!!!!!!!!!!!!!!!!!!!!!!!!debugging!!!!!!!!!!!!!!!!!!!!!!!!!!MP1 
       ! if (berg%iceberg_num .eq. 1) then
       !   print *, 'Comparing longitudes: ', lon1, lon2, r_dist_x, dlon
       !   print *, 'Comparing latitudes: ', lat1, lat2, r_dist_y, dlat
       !   print *, 'Outside, iceberg_num, r_dist', berg%iceberg_num, r_dist,bonded
       !   print *, 'Halo_status', berg%halo_berg,other_berg%halo_berg 
       ! endif
       ! print *, 'outside the loop',R1, R2,r_dist, bonded
       !!!!!!!!!!!!!!!!!!!!!!!!!!!debugging!!!!!!!!!!!!!!!!!!!!!!!!!! 


       !call overlap_area(R1,R2,r_dist,A_o,trapped)
       !T_min=min(T1,T2)
       !A_min = min((pi*R1**R1),(pi*R2*R2)) 
       M_min=min(M1,M2)
       !Calculating spring force  (later this should only be done on the first time around)
       if ((r_dist>0.) .AND. ((r_dist< (R1+R2).AND. (.not. bonded)) .OR. ( (r_dist> (R1+R2)) .AND. (bonded) ) )) then
         !Spring force
         !accel_spring=spring_coef*(T_min/T1)*(A_o/A1) ! Old version dependent on area
         accel_spring=spring_coef*(M_min/M1)*(R1+R2-r_dist)
         IA_x=IA_x+(accel_spring*(r_dist_x/r_dist))
         IA_y=IA_y+(accel_spring*(r_dist_y/r_dist))
         

        if (r_dist < 5*(R1+R2)) then
          
          !MP1
          !if (berg%iceberg_num .eq. 1) then
          !  !print *,  '************************************************************'
          !  print *, 'INSIDE, r_dist', berg%iceberg_num, other_berg%iceberg_num, r_dist, bonded
          !endif
          !print *, 'in the loop1', spring_coef, (M_min/M1), accel_spring,(R1+R2-r_dist) 
          !print *, 'in the loop2', IA_x, IA_y, R1, R2,r_dist, berg%iceberg_num,other_berg%iceberg_num
          !Damping force:
          !Paralel velocity
           P_11=(r_dist_x*r_dist_x)/(r_dist**2)
           P_12=(r_dist_x*r_dist_y)/(r_dist**2)
           P_21=(r_dist_x*r_dist_y)/(r_dist**2)
           P_22=(r_dist_y*r_dist_y)/(r_dist**2)
           !p_ia_coef=radial_damping_coef*(T_min/T1)*(A_min/A1)
           p_ia_coef=radial_damping_coef*(M_min/M1)
           p_ia_coef=p_ia_coef*(0.5*(sqrt((((P_11*(u2-u1))+(P_12*(v2-v1)))**2)+ (((P_12*(u2-u1))+(P_22*(v2-v1)))**2)) &
           + sqrt((((P_11*(u2-u0))+(P_12*(v2-v0)))**2)+(((P_12*(u2-u0)) +(P_22*(v2-v0)))**2))))
          
           P_ia_11=P_ia_11+p_ia_coef*P_11
           P_ia_12=P_ia_12+p_ia_coef*P_12
           P_ia_21=P_ia_21+p_ia_coef*P_21
           P_ia_22=P_ia_22+p_ia_coef*P_22
           P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
           P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))
           !print *, 'Paralel: ',berg%iceberg_num,  p_ia_coef, IA_x, P_ia_11, P_ia_21,P_ia_12, P_ia_22

           !Normal velocities
           P_11=1-P_11  ;  P_12=-P_12 ; P_21= -P_21 ;    P_22=1-P_22
           !p_ia_coef=tangental_damping_coef*(T_min/T1)*(A_min/A1)
           p_ia_coef=tangental_damping_coef*(M_min/M1)
           p_ia_coef=p_ia_coef*(0.5*(sqrt((((P_11*(u2-u1))+(P_12*(v2-v1)))**2)+ (((P_12*(u2-u1))+(P_22*(v2-v1)))**2))  &
           + sqrt((((P_11*(u2-u0))+(P_12*(v2-v0)))**2)+(((P_12*(u2-u0)) +(P_22*(v2-v0)))**2))))
           P_ia_11=P_ia_11+p_ia_coef*P_11
           P_ia_12=P_ia_12+p_ia_coef*P_12
           P_ia_21=P_ia_21+p_ia_coef*P_21
           P_ia_22=P_ia_22+p_ia_coef*P_22
           P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
           P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))
           !print *, 'Perp: ',berg%iceberg_num,  p_ia_coef, IA_x, P_ia_11, P_ia_21,P_ia_12, P_ia_22
           !print *, 'P_11',P_11
           !print *, 'P_21',P_21
           !print *, 'P_12',P_12
           !print *, 'P_22',P_22
          endif
        endif
      endif


      end subroutine calculate_force


  subroutine overlap_area(R1,R2,d,A,trapped)
    real, intent(in) :: R1, R2, d
    real, intent(out) :: A, Trapped
    real :: R1_sq, R2_sq, d_sq  
    R1_sq=R1**2
    R2_sq=R2**2
    d_sq=d**2
    Trapped=0.

    if (d>0.) then
      if (d<(R1+R2)) then
        if (d>abs(R1-R2)) then
          A= (R1_sq*acos((d_sq+R1_sq-R2_sq)/(2.*d*R1)))  +  (R2_sq*acos((d_sq+R2_sq-R1_sq)/(2.*d*R2)))  - (0.5*sqrt((-d+R1+R2)*(d+R1-R2)*(d-R1+R2)*(d+R1+R2)))
        else
          A=min(pi*R1_sq,pi*R2_sq) 
          Trapped=1.
        endif
      else
        A=0.
      endif
    else
      A=0.     ! No area of perfectly overlapping bergs (ie: a berg interacting with itself)
    endif

   end subroutine overlap_area




end subroutine interactive_force


! ##############################################################################


subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, axn, ayn, bxn, byn, debug_flag) !Saving  acceleration for Verlet, Adding Verlet flag - Alon  MP1
!subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, debug_flag) !old version commmented out by Alon
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
integer, intent(in) :: i, j
real, intent(in) :: xi, yj, lat, uvel, vvel, uvel0, vvel0, dt
real, intent(out) :: ax, ay
real, intent(inout) :: axn, ayn, bxn, byn ! Added implicit and explicit accelerations to output -Alon
logical, optional :: debug_flag
! Local variables
type(icebergs_gridded), pointer :: grd
real :: uo, vo, ui, vi, ua, va, uwave, vwave, ssh_x, ssh_y, sst, sss, cn, hi
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


!print *, 'axn=',axn,'ayn=',ayn
  u_star=uvel0+(axn*(dt/2.))  !Alon
  v_star=vvel0+(ayn*(dt/2.))  !Alon

  ! Get the stderr unit number.
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  ! Interpolate gridded fields to berg     - Note: It should be possible to move this to evolve, so that it only needs to be called once. !!!!
  call interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, sss, cn, hi)

  if ((grd%grid_is_latlon) .and. (.not. bergs%use_f_plane)) then
     f_cori=(2.*omega)*sin(pi_180*lat)
  else
     f_cori=(2.*omega)*sin(pi_180*bergs%lat_ref)
  endif
!  f_cori=0.

  M=berg%mass
  T=berg%thickness ! total thickness
  D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
  F=T-D ! freeboard
  W=berg%width
  L=berg%length

!Initializing accelerations - Alon. (I am not 100% sure this is needed). I'm not sure what is output if variable is not defined in the subroutine.
  axn=0.
  ayn=0.
  bxn=0.
  byn=0.

  hi=min(hi,D)
  D_hi=max(0.,D-hi)

  ! Wave radiation
  uwave=ua-uo; vwave=va-vo  ! Use wind speed rel. to ocean for wave model (aja)?
  wmod=uwave*uwave+vwave*vwave ! The wave amplitude and length depend on the wind speed relative to the ocean current;
                               !  actually wmod is wmod**2 here.
  ampl=0.5*0.02025*wmod ! This is "a", the wave amplitude
  Lwavelength=0.32*wmod ! Surface wave length fitted to data in table at
  !      http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
  Lcutoff=0.125*Lwavelength
  Ltop=0.25*Lwavelength
  Cr=Cr0*min(max(0.,(L-Lcutoff)/((Ltop-Lcutoff)+1.e-30)),1.) ! Wave radiation coefficient
  !     fitted to graph from Carrieres et al.,  POAC Drift Model.
  wave_rad=0.5*rho_seawater/M*Cr*gravity*ampl*min(ampl,F)*(2.*W*L)/(W+L)
  wmod = sqrt(ua*ua+va*va) ! Wind speed
  if (wmod.ne.0.) then
    uwave=ua/wmod ! Wave radiation force acts in wind direction ...
    vwave=va/wmod
  else
    uwave=0.; vwave=0.; wave_rad=0. ! ... and only when wind is present.
  endif

  ! Weighted drag coefficients
  c_ocn=rho_seawater/M*(0.5*Cd_wv*W*(D_hi)+Cd_wh*W*L)
  c_atm=rho_air     /M*(0.5*Cd_av*W*F     +Cd_ah*W*L)
  if (abs(hi).eq.0.) then
    c_ice=0.
  else
    c_ice=rho_ice     /M*(0.5*Cd_iv*W*hi              )
  endif
  if (abs(ui)+abs(vi).eq.0.) c_ice=0.

!Turning drag off for testing - Alon
!c_ocn=0.
!c_atm=0.
!c_ice=0.

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
      us=uveln   ; vs=vveln
    endif
  if (use_new_predictive_corrective) then
    !Alon's proposed change - using Bob's improved scheme.
    drag_ocn=c_ocn*0.5*(sqrt( (uveln-uo)**2+(vveln-vo)**2 )+sqrt( (uvel0-uo)**2+(vvel0-vo)**2 ))
    drag_atm=c_atm*0.5*(sqrt( (uveln-ua)**2+(vveln-va)**2 )+sqrt( (uvel0-ua)**2+(vvel0-va)**2 ))
    drag_ice=c_ice*0.5*(sqrt( (uveln-ui)**2+(vveln-vi)**2 )+sqrt( (uvel0-ui)**2+(vvel0-vi)**2 ))
  else
    !Original Scheme
    us=0.5*(uveln+uvel); vs=0.5*(vveln+vvel)
    drag_ocn=c_ocn*sqrt( (us-uo)**2+(vs-vo)**2 )
    drag_atm=c_atm*sqrt( (us-ua)**2+(vs-va)**2 )
    drag_ice=c_ice*sqrt( (us-ui)**2+(vs-vi)**2 )
  endif

     RHS_x=(axn/2) + bxn
     RHS_y=(ayn/2) + byn

     if (beta>0.) then ! If implicit, use u_star, v_star rather than RK4 latest
         RHS_x=RHS_x - drag_ocn*(u_star-uo) -drag_atm*(u_star-ua) -drag_ice*(u_star-ui) 
         RHS_y=RHS_y - drag_ocn*(v_star-vo) -drag_atm*(v_star-va) -drag_ice*(v_star-vi)  
    else
         RHS_x=RHS_x - drag_ocn*(uvel-uo) -drag_atm*(uvel-ua) -drag_ice*(uvel-ui)  
         RHS_y=RHS_y - drag_ocn*(vvel-vo) -drag_atm*(vvel-va) -drag_ice*(vvel-vi)  
    endif

if (interactive_icebergs_on) then
   if (itloop>1) then
        call interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, us,vs, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
    endif
     if (beta>0.) then ! If implicit, use u_star, v_star rather than RK4 latest
         RHS_x=RHS_x -(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x) 
         RHS_y=RHS_y -(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y) 
    else
         RHS_x=RHS_x - (((P_ia_11*uvel)+(P_ia_12*vvel))-P_ia_times_u_x)
         RHS_y=RHS_y - (((P_ia_21*uvel)+(P_ia_22*vvel))-P_ia_times_u_y)       
    endif
    !print *,'Before calculation:', berg%iceberg_num, IA_x, IA_y, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y
endif


  ! Solve for implicit accelerations
    if (alpha+beta.gt.0.) then
      lambda=drag_ocn+drag_atm+drag_ice
      A11=1.+beta*dt*lambda
      A22=1.+beta*dt*lambda
      A12=-alpha*dt*f_cori  
      A21=alpha*dt*f_cori  
      !A12=dt*f_cori  !Removed by ALon (in order to have the entire matrix. I hope the sign is correct)

      if (C_N>0.) then   !For Crank-Nicolson Coriolis term.   
          A12=A12/2.
          A21=A21/2.
      endif

      if (interactive_icebergs_on) then
            A11=A11+(dt*P_ia_11)
            A12=A12+(dt*P_ia_12)
            A21=A21+(dt*P_ia_21)
            A22=A22+(dt*P_ia_22)
     endif

     !This is for testing the code using only interactive forces
     if (bergs%only_interactive_forces) then
       RHS_x=(IA_x/2) -(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x) 
       RHS_y=(IA_y/2) -(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y) 
       A11=1+(dt*P_ia_11)
       A12=(dt*P_ia_12)
       A21=(dt*P_ia_21)
       A22=1+(dt*P_ia_22)
     endif



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

      uveln=u_star+dt*ax        ! Alon
      vveln=v_star+dt*ay        ! Alon
!MP4
  !    if (berg%iceberg_num .eq. 1) then
  !      print *, '***************************************************'
  !      print *,'Iceberg_num, itloop', berg%iceberg_num, itloop
  !      print *, 'P matrix:', P_ia_11, P_ia_12,P_ia_21,P_ia_22,P_ia_times_u_x, P_ia_times_u_x 
  !      print *,'A_matrix', A11, A12, A21, A22 
  !      print *,'IA_x IA_y', IA_x, IA_y 
  !      print *, 'RHS, ustar, uvel,ax: ', RHS_x, u_star,uveln, ax
  !    endif
  enddo ! itloop

!Saving the totally explicit part of the acceleration to use in finding the next position and u_star -Alon
    axn=0.
    ayn=0.
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

     !This is for testing the code using only interactive forces
     if (bergs%only_interactive_forces) then
       axn=IA_x
       ayn=IA_y
     endif

    bxn= ax-(axn/2) !Alon
    byn= ay-(ayn/2) !Alon

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
          call error_mesg('diamonds, Speeding icebergs', 'Faster than the CFL!', WARNING)
          write(stderrunit,*) 'diamonds, Speeding berg1! =',mpp_pe(), berg%iceberg_num
          write(stderrunit,*) 'diamonds, Speeding berg2, speed =',speed, loc_dx/dt
          write(stderrunit,*) 'diamonds, Speeding berg3, lat, lon =',lat,xi,yj
        endif
      endif
    endif
  endif

  dumpit=.false.
  if (abs(uveln)>vel_lim.or.abs(vveln)>vel_lim) then
    dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive velocity'
  endif
  if (abs(ax)>accel_lim.or.abs(ay)>accel_lim) then
    dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive acceleration'
  endif
  if (present(debug_flag)) then
    if (debug_flag) dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Debug dump flagged by arguments'
  endif
  if (dumpit) then
 100 format('pe=',i3,a15,9(x,a8,es12.3))
 200 format('pe=',i3,a15,(x,a8,i12),9(x,a8,es12.3))
    write(stderrunit,200) mpp_pe(),'Starting pars:', &
      'yr0=',berg%start_year, 'day0=',berg%start_day, &
      'lon0=',berg%start_lon, 'lat0=',berg%start_lat, 'mass0=',berg%start_mass, &
      'sclng=',berg%mass_scaling, 'num0=',berg%iceberg_num
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
    call print_berg(stderrunit,berg,'diamonds, accel, large accel')
  endif

  !Used for testing the ocean response to fixed iceberg motion.
  if (bergs%override_iceberg_velocities) then
    ax  = 0.0;  ay  = 0.0;
    axn = 0.0;  ayn = 0.0;
    bxn = 0.0;  byn = 0.0;
  endif

  contains

  subroutine dump_locfld(grd,i0,j0,A,lbl)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i0, j0
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: A
  character(len=*) :: lbl
! Local variables
  integer :: i, j, ii, jj
  real :: B(-1:1,-1:1), fac

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
  
  subroutine dump_locvel(grd,i0,j0,A,lbl)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i0, j0
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: A
  character(len=*) :: lbl
! Local variables
  integer :: i, j, ii, jj
  real :: B(-1:0,-1:0), fac

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
  
end subroutine accel

! ##############################################################################

subroutine thermodynamics(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: M, T, W, L, SST, Vol, Ln, Wn, Tn, nVol, IC, Dn
real :: Mv, Me, Mb, melt, dvo, dva, dM, Ss, dMe, dMb, dMv
real :: Mnew, Mnew1, Mnew2, Hocean
real :: Mbits, nMbits, dMbitsE, dMbitsM, Lbits, Abits, Mbb
real :: tip_parameter
real :: Delta, q
integer :: i,j, stderrunit
type(iceberg), pointer :: this, next
real, parameter :: perday=1./86400.
integer :: grdi, grdj
real :: SSS !Temporarily here

  ! For convenience
  grd=>bergs%grd
  
  !Initializing 
  grd%Uvel_on_ocean(:,:,:)=0.
  grd%Vvel_on_ocean(:,:,:)=0.

  ! Thermodynamics of first halo row is calculated, so that spread mass to ocean works correctly
  do grdj = grd%jsc-1,grd%jec+1 ; do grdi = grd%isc-1,grd%iec+1  
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))
      if (debug) call check_position(grd, this, 'thermodynamics (top)')
      call interp_flds(grd, this%ine, this%jne, this%xi, this%yj, this%uo, this%vo, &
              this%ui, this%vi, this%ua, this%va, this%ssh_x, this%ssh_y, this%sst, &
              this%sss,this%cn, this%hi)
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

      !For icebergs acting as ice shelves 
      if (bergs%melt_icebergs_as_ice_shelf) then
        Mv=0.0
        Mb=0.0
        Me=0.0
        if (.not. bergs%use_mixed_layer_salinity_for_thermo)  SSS=35.0  
        call find_basal_melt(bergs,dvo,this%lat,SSS,SST,bergs%Use_three_equation_model,T,Mb,this%iceberg_num)
        Mb=max(Mb,0.) !No refreezing allowed for now
        !Set melt to zero if ocean is too thin.
        if ((bergs%melt_cutoff >=0.) .and. (bergs%apply_thickness_cutoff_to_bergs_melt)) then
          Dn=(bergs%rho_bergs/rho_seawater)*this%thickness ! draught (keel depth)
          if ((grd%ocean_depth(i,j)-Dn) < bergs%melt_cutoff) then
            Mb=0.
          endif
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
  
        Ln=max(L-Mv*bergs%dt,0.) ! new length (m)
        Wn=max(W-Mv*bergs%dt,0.) ! new width (m)
        nVol=Tn*Wn*Ln ! new volume (m^3)
        Mnew2=(nVol/Vol)*M ! new mass (kg)
        dMv=Mnew1-Mnew2 ! mass lost to buoyant convection (>0) (kg)
  
        Ln=max(Ln-Me*bergs%dt,0.) ! new length (m)
        Wn=max(Wn-Me*bergs%dt,0.) ! new width (m)
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
  
      ! Bergy bits
      if (bergs%bergy_bit_erosion_fraction>0.) then
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
      else
        Abits=0.
        dMbitsE=0.
        dMbitsM=0.
        nMbits=this%mass_of_bits ! retain previous value incase non-zero
      endif
  
      ! Add melting to the grid and field diagnostics
      if (grd%area(i,j).ne.0.) then
        melt=(dM-(dMbitsE-dMbitsM))/bergs%dt ! kg/s
        grd%floating_melt(i,j)=grd%floating_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
        melt=melt*this%heat_density ! kg/s x J/kg = J/s
        grd%calving_hflx(i,j)=grd%calving_hflx(i,j)+melt/grd%area(i,j)*this%mass_scaling ! W/m2
        bergs%net_heat_to_ocean=bergs%net_heat_to_ocean+melt*this%mass_scaling*bergs%dt ! J
        melt=dM/bergs%dt ! kg/s
        grd%berg_melt(i,j)=grd%berg_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
        melt=dMbitsE/bergs%dt ! mass flux into bergy bits in kg/s
        grd%bergy_src(i,j)=grd%bergy_src(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
        melt=dMbitsM/bergs%dt ! melt rate of bergy bits in kg/s
        grd%bergy_melt(i,j)=grd%bergy_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
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
      else
        stderrunit = stderr()
        write(stderrunit,*) 'diamonds, thermodynamics: berg appears to have grounded!!!! PE=',mpp_pe(),i,j
        call print_berg(stderrunit,this,'thermodynamics, grounded')
        if (associated(this%trajectory)) &
          write(stderrunit,*) 'traj=',this%trajectory%lon,this%trajectory%lat
        write(stderrunit,*) 'msk=',grd%msk(i,j),grd%area(i,j)
        call error_mesg('diamonds, thermodynamics', 'berg appears to have grounded!', FATAL)
      endif
  
 
   ! Rolling 
   !There are now 3 iceberg rolling schemes:
   !1) Rolling based on aspect ratio threshold (iceberg of constant density)
   !2) Rolling based on corrected Weeks and Mellor scheme
   !3) Rolling based on incorrect Weeks and Mellor scheme - kept for legacy reasons
    if (bergs%allow_bergs_to_roll) then
      Dn=(bergs%rho_bergs/rho_seawater)*Tn ! draught (keel depth)
      if ( Dn>0. ) then
        if ( (.not.bergs%use_updated_rolling_scheme) .and. (bergs%tip_parameter<999.) ) then    !Use Rolling Scheme 3
            if ( max(Wn,Ln)<sqrt(0.92*(Dn**2)+58.32*Dn) ) then
              call swap_variables(Tn,Wn)
            endif
        else
          if (Wn>Ln) call swap_variables(Ln,Wn)  !Make sure that Wn is the smaller dimension
        
          if ( (.not.bergs%use_updated_rolling_scheme) .and. (bergs%tip_parameter>=999.) ) then    !Use Rolling Scheme 2
            q=bergs%rho_bergs/rho_seawater
            Delta=6.0
            if (Wn<sqrt((6.0*q*(1-q)*(Tn**2))-(12*Delta*q*Tn))) then
                call swap_variables(Tn,Wn)
            endif
          endif

          if (bergs%use_updated_rolling_scheme) then    !Use Rolling Scheme 1
            if (bergs%tip_parameter>0.) then
              tip_parameter=bergs%tip_parameter
            else
              ! Equation 27 from Burton et al 2012, or equivolently, Weeks and Mellor 1979 with constant density
              tip_parameter=sqrt(6*(bergs%rho_bergs/rho_seawater)*(1-(bergs%rho_bergs/rho_seawater)))   !using default values gives 0.92
            endif
            if ((tip_parameter*Tn)>Wn)  then     !note that we use the Thickness instead of the Draft
                call swap_variables(Tn,Wn)
            endif
          endif
        endif
        Dn=(bergs%rho_bergs/rho_seawater)*Tn ! re-calculate draught (keel depth) for grounding
      endif
    endif

       !This option allows iceberg melt fluxes to enter the ocean without the icebergs changing shape
       if (bergs%Iceberg_melt_without_decay) then
         !In this case, the iceberg dimension are reset to their values before
         !the thermodynamics are applied. 
         !If the spread_mass is being used to calculate melt, we calculate this
         !before reseting 
         if (bergs%find_melt_using_spread_mass) then
           if (Mnew>0.) then !If the berg still exists
             call spread_mass_across_ocean_cells(bergs,this, i, j, this%xi, this%yj,Mnew , nMbits, this%mass_scaling, Ln*Wn,  Tn)
           endif
         endif
         !Reset all the values
         Mnew=this%mass
         nMbits=this%mass_of_bits
         Tn=this%thickness
         Wn=this%width
         Ln=this%length
         if (bergs%bergy_bit_erosion_fraction>0.) then
           Mbits=this%mass_of_bits ! mass of bergy bits (kg)
           Lbits=min(L,W,T,40.) ! assume bergy bits are smallest dimension or 40 meters
           Abits=(Mbits/bergs%rho_bergs)/Lbits ! Effective bottom area (assuming T=Lbits)
         endif
      else
        ! Store the new state of iceberg (with L>W)
        this%mass=Mnew
        this%mass_of_bits=nMbits
        this%thickness=Tn
        this%width=min(Wn,Ln)
        this%length=max(Wn,Ln)
      endif
      next=>this%next
  
      ! Did berg completely melt?
      if (Mnew<=0.) then ! Delete the berg
        call move_trajectory(bergs, this)
        call delete_iceberg_from_list(bergs%list(grdi,grdj)%first, this)
        bergs%nbergs_melted=bergs%nbergs_melted+1
      endif
      this=>next
    enddo
  enddo ; enddo


  contains

  subroutine swap_variables(x,y)
  ! Arguments
  real, intent(inout) :: x, y
  real :: temp
  temp=x
  x=y
  y=temp
  end subroutine swap_variables

end subroutine thermodynamics


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_gridded_icebergs_fields(bergs)
! Arguments
type(icebergs), pointer :: bergs
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

  !Dividng the gridded iceberg momentum diagnostic by the iceberg mass to get velocities
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_basal_melt(bergs,dvo,lat,salt,temp,Use_three_equation_model,thickness,basal_melt,iceberg_num)
  ! Arguments
  type(icebergs), pointer :: bergs
  ! Local variables
  real , intent(out) :: basal_melt !Melt rate underneath the icebergs
  real , intent(in) :: dvo !Speed of iceberg relative to ocean mixed layer
  real , intent(in) :: salt !Salinity of mixed layer
  real , intent(in) :: temp !Temperature of mixed layer
  real , intent(in) :: lat !Latitude (for boundary layer calculation)
  integer , intent(in) :: iceberg_num !Iceberg number, used for debugging (error messages)
  real , intent(in) :: thickness !Ice thickness - needed to work out the pressure below the ice
  logical , intent(in) :: Use_three_equation_model !True uses the 3 equation model, False uses the 2 equation model.
  
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
!
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
!
  ! Variables used in iterating for wB_flux.
  real :: wB_flux_new, DwB, dDwB_dwB_in
  real :: I_Gam_T, I_Gam_S
  real :: dG_dwB, iDens
  logical :: Sb_min_set, Sb_max_set
  logical :: out_of_bounds

  real, parameter :: c2_3 = 2.0/3.0
  integer ::  it1, it3
 
  !Parameters copied ice shelf module defaults (could be entered in the namelist later)
  real, parameter ::  dR0_dT = -0.038357 ! Partial derivative of the mixed layer density with temperature, in units of kg m-3 K-1. 
  real, parameter ::  dR0_dS = 0.805876 ! Partial derivative of the mixed layer density with salinity, in units of kg m-3 psu-1.
  real, parameter ::  RHO_T0_S0 = 999.910681 ! Density of water with T=0, S=0 for linear EOS
  real, parameter :: Salin_Ice =0.0 !Salinity of ice
  real, parameter :: Temp_Ice = -15.0 !Salinity of ice
  real, parameter :: kd_molec_salt=  8.02e-10 !The molecular diffusivity of salt in sea water at the freezing point
  real, parameter :: kd_molec_temp=  1.41e-7 !The molecular diffusivity of heat in sea water at the freezing point
  real, parameter :: kv_molec=  1.95e-6 !The molecular molecular kinimatic viscosity of sea water at the freezing point
  real, parameter :: Cp_Ice =  2009.0 !Specific heat capacity of ice, taking from HJ99 (Holland and Jenkins 1999)
  real, parameter :: Cp_ml =  3974.0 !Specific heat capacity of mixed layer, taking from HJ99 (Holland and Jenkins 1999)
  real, parameter :: LF =  3.335e5 !Latent heat of fusion, taken from HJ99 (Holland and Jenkins 1999)
  real, parameter :: gamma_t =  0.0 ! Exchange velcoity used in 2 equation model. Whn gamma_t is >0, the exchange velocity is independ of u_star. 
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
    ! temperature and salinty gradients


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
            call error_mesg('diamonds,Find basal melt', 'shelf_calc_flux: Irregular iteration for Sbdry (max).' ,WARNING)
            print *, 'Sbdry error: iceberg_num,dvo,temp,salt,lat,thickness :',iceberg_num,dvo,temp,salt,lat,thickness
          endif
          out_of_bounds=.true.
          exit
        endif
        Sb_max = Sbdry ; dS_max = dS_it ; Sb_max_set = .true.
      else ! Sbdry is now the lower bound.
        if (Sb_min_set .and. (Sbdry < Sb_min)) then
          if (debug) then
            call error_mesg('diamonds,Find basal melt', 'shelf_calc_flux: Irregular iteration for Sbdry (min).' ,WARNING)
            print *, 'Sbdry error: iceberg_num,dvo,temp,salt,lat,thickness :',iceberg_num,dvo,temp,salt,lat,thickness
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

contains

  subroutine calculate_TFreeze(S, pres, T_Fr)
    !Arguments
    real,    intent(in)  :: S, pres
    real,    intent(out) :: T_Fr
    real, parameter :: dTFr_dp    = -7.53E-08    !DTFREEZE_DP in MOM_input
    real, parameter :: dTFr_dS    = -0.0573      !DTFREEZE_DS in MOM_input
    real, parameter :: TFr_S0_P0  =0.0832        !TFREEZE_S0_P0 in MOM_input
    !    This subroutine computes the freezing point potential temparature
    !  (in deg C) from salinity (in psu), and pressure (in Pa) using a simple
    !  linear expression, with coefficients passed in as arguments.
    !  Copied from subroutine calculate_TFreeze_linear_scalar (in MOM/equation_of_state)
    !
    ! Arguments: S - salinity in PSU.
    !  (in)      pres - pressure in Pa.
    !  (out)     T_Fr - Freezing point potential temperature in deg C.
    !  (in)      TFr_S0_P0 - The freezing point at S=0, p=0, in deg C.
    !  (in)      dTFr_dS - The derivatives of freezing point with salinity, in
    !                      deg C PSU-1.
    !  (in)      dTFr_dp - The derivatives of freezing point with pressure, in
    !                      deg C Pa-1.
    T_Fr = (TFr_S0_P0 + dTFr_dS*S) + dTFr_dp*pres
  end subroutine calculate_TFreeze

  subroutine calculate_density(T, S, pressure, rho, Rho_T0_S0, dRho_dT, dRho_dS)
    !Arguments
    real,    intent(in)  :: T, S, pressure
    real,    intent(out) :: rho
    real,    intent(in)  :: Rho_T0_S0, dRho_dT, dRho_dS
    ! *  This subroutine computes the density of sea water with a trivial  *
    ! *  linear equation of state (in kg/m^3) from salinity (sal in psu),  *
    ! *  potential temperature (T in deg C), and pressure in Pa.           *
    !    Copied from subroutine calculate_density_scalar_linear (in MOM/equation_of_state)
    ! *                                                                    *
    ! * Arguments: T - potential temperature relative to the surface in C. *
    ! *  (in)      S - salinity in PSU.                                    *
    ! *  (in)      pressure - pressure in Pa.                              *
    ! *  (out)     rho - in situ density in kg m-3.                        *
    ! *  (in)      start - the starting point in the arrays.               *
    ! *  (in)      npts - the number of values to calculate.               *
    ! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
    ! *  (in)      dRho_dT - The derivatives of density with temperature   *
    ! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *
    rho = Rho_T0_S0 + dRho_dT*T + dRho_dS*S
  end subroutine calculate_density

end subroutine find_basal_melt
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine find_orientation_using_iceberg_bonds(grd,berg,orientation)
  ! Arguments
  type(iceberg), pointer :: berg
  real, intent(inout) :: orientation 
  type(icebergs_gridded), pointer :: grd
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
          !print *, 'Iceberg bond details:',berg%iceberg_num, current_bond%other_berg_num,berg%halo_berg, mpp_pe()
          !print *, 'Iceberg bond details2:',berg%ine, berg%jne, current_bond%other_berg_ine, current_bond%other_berg_jne
          !print *, 'Iceberg isd,ied,jsd,jed:',grd%isd, grd%ied, grd%jsd, grd%jed
          !print *, 'Iceberg isc,iec,jsc,jec:',grd%isc, grd%iec, grd%jsc, grd%jec
          !call error_mesg('diamonds,calculating orientation', 'Looking at bond interactions of unassosiated berg!' ,FATAL)
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

subroutine spread_mass_across_ocean_cells(bergs, berg, i, j, x, y, Mberg, Mbits, scaling, Area, Tn)
  ! Arguments
  type(icebergs), pointer :: bergs
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg
  integer, intent(in) :: i, j
  real, intent(in) :: x, y, Mberg, Mbits, scaling, Area
  real, intent(in) :: Tn
  ! Local variables
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

  ! Get the stderr unit number
  stderrunit = stderr()

  tol=1.e-10
  grd=>bergs%grd
  Mass_berg=Mberg

  !Trimming icebergs to account for grounded fraction.
  if (bergs%grounding_fraction>0.) then
    Hocean=bergs%grounding_fraction*(grd%ocean_depth(i,j)+grd%ssh(i,j))
    Dn=(bergs%rho_bergs/rho_seawater)*Tn ! re-calculate draught (keel depth)
    if (Dn>Hocean) Mass_berg=Mberg*min(1.,Hocean/Dn)
  endif

  Mass=(Mass_berg+Mbits)*scaling
  ! This line attempts to "clip" the weight felt by the ocean. The concept of
  ! clipping is non-physical and this step should be replaced by grounding.
  if (grd%clipping_depth>0.) Mass=min(Mass,grd%clipping_depth*grd%area(i,j)*rho_seawater)
 
  !Initialize weights for each cell       
  yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0.
  yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

  if (.not. bergs%hexagonal_icebergs) then  !Treat icebergs as rectangles of size L:  (this is the default)
    
    !L is the non dimensional length of the iceberg [  L=(Area of berg/ Area of grid cell)^0.5  ] or something like that.
    if (grd%area(i,j)>0) then
      L=min( sqrt(Area / grd%area(i,j)),1.0)
    else 
      L=1.
    endif
  
    if (bergs%use_old_spreading) then
      !Old version before icebergs were given size L
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
    
    fraction_used=1.  !rectangular bergs do share mass with boundaries (all mass is included in cells)

  else !Spread mass as if elements area hexagonal
    
    orientation=bergs%initial_orientation
    if ((bergs%iceberg_bonds_on) .and. (bergs%rotate_icebergs_for_mass_spreading)) call find_orientation_using_iceberg_bonds(grd,berg,orientation)

    if (grd%area(i,j)>0) then
      H = min(( (sqrt(Area/(2.*sqrt(3.)))  / sqrt(grd%area(i,j)))),1.) ;  !Non dimensionalize element length by grid area. (This gives the non-dim Apothen of the hexagon)
    else 
      H= (sqrt(3.)/2)*(0.49)  !Largest allowable H, since this makes S=0.49, and S has to be less than 0.5  (Not sure what the implications of this are)
    endif
    S=(2/sqrt(3.))*H !Side of the hexagon

    if (S>0.5) then
      !The width of an iceberg should not be greater than half the gridcell, or else it can spread over 3 cells  (i.e. S must be less than 0.5 nondimensionally)
      !print 'Elements must be smaller than a whole gridcell', 'i.e.: S= ' , S , '>=0.5'
      call error_mesg('diamonds, hexagonal spreading', 'Diameter of the iceberg is larger than a grid cell. Use smaller icebergs', WARNING)
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
      call error_mesg('diamonds, hexagonal spreading', 'Intersection with hexagons should not be negative!!!', WARNING)
      write(stderrunit,*) 'diamonds, yU,yC,yD', Area_Q1, Area_Q2, Area_Q3, Area_Q4
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
            !write(stderrunit,*) 'diamonds, You are in the hexagonal domain now!!!'
    endif

      !Double check that all the mass is being used.
      if ((abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>tol) .and. (mpp_pe().eq. mpp_root_pe())) then
        !call error_mesg('diamonds, hexagonal spreading', 'All the mass is not being used!!!', WARNING)
        write(stderrunit,*) 'diamonds, hexagonal, H,x0,y0', H, x0 , y0
        write(stderrunit,*) 'diamonds, hexagonal, Areas',(Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
        debug=.True.
        !call Hexagon_into_quadrants_using_triangles(x0,y0,H,orientation,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4, debug)
        call error_mesg('diamonds, hexagonal spreading', 'All the mass is not being used!!!', FATAL)
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

subroutine spread_variable_across_cells(grd, variable_on_ocean, Var,i,j, &
           yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR,I_fraction_used)
  ! Arguments
  type(icebergs_gridded), pointer, intent(in) :: grd
  real, dimension(grd%isd:grd%ied, grd%jsd:grd%jed, 9), intent(inout) :: variable_on_ocean 
  real, intent(in) :: Var !Variable to be spread accross cell
  real, intent(in) :: yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR !Weights
  real, intent(in) :: I_fraction_used !Amount of iceberg used (inverse)
  integer, intent(in) :: i, j

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


! ##############################################################################

real function Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)
  ! Arguments
  real, intent(in) :: Ax,Ay,Bx,By,Cx,Cy
  Area_of_triangle    =   abs(    0.5*((Ax*(By-Cy))+(Bx*(Cy-Ay))+(Cx*(Ay-By))) );
end function Area_of_triangle

real function roundoff(x,sig_fig)
  ! Arguments
  real, intent(in) :: x
  integer, intent(in) :: sig_fig
    !roundoff=round(x*(10**(sig_fig))
    roundoff=(FLOAT (INT(x * (10.**sig_fig) + 0.5)) / (10.**sig_fig))
end function roundoff

logical function point_in_interval(Ax,Ay,Bx,By,px,py)
  ! Arguments
  real, intent(in) :: Ax,Ay,Bx,By,px,py
  point_in_interval=.False.
  if ((px <= max(Ax,Bx)) .and. (px >= min(Ax,Bx))) then
    if ((py <= max(Ay,By)) .and. (py >= min(Ay,By))) then
      point_in_interval=.True.
    endif
  endif
end function point_in_interval


logical function point_is_on_the_line(Ax,Ay,Bx,By,qx,qy)
  ! Arguments
  real, intent(in) :: Ax,Ay,Bx,By,qx,qy
  real :: tol, dxc,dyc,dxl,dyl,cross
    !tol=1.e-12;
    tol=0.0;
    dxc = qx - Ax;
    dyc = qy - Ay;
    dxl = Bx - Ax;
    dyl = By - Ay;
    cross = dxc * dyl - dyc * dxl;
    if (abs(cross)<=tol) then
      point_is_on_the_line=.True.
    else
     point_is_on_the_line=.False.
    endif
end function point_is_on_the_line

logical function point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,qx,qy)
  !This function decides whether a point (qx,qy) is inside the triangle ABC.
  !There is also the option to include the boundary of the triangle.
  ! Arguments
  real, intent(in) :: Ax,Ay,Bx,By,Cx,Cy,qx,qy
  real :: l0,l1,l2,p0,p1,p2
  real :: v0x,v1x,v2x,v0y,v1y,v2y,dot00,dot01,dot02,dot11,dot12
  
  point_in_triangle = .False.
  if ((Ax==qx .and. Ay==qy) .or. (Bx==qx .and. By==qy) .or. (Cx==qx .and. Cy==qy)) then  !Exclude the pathelogical case
      point_in_triangle = .False.
  else
    if (((point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) .or. (point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy))) .or. (point_is_on_the_line(Bx,By,Cx,Cy,qx,qy)))) then
      point_in_triangle = .False.
    else
      !Compute point in triangle using Barycentric coordinates (the same as sum_sign_dot_prod routines)
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


subroutine Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axis1,Area_positive, Area_negative)  !You should change this name a little, so that it not similar the other routine.
!This function calculates the area of a triangle on opposited sides of an axis when the triangle is split with two points on one side, and one point on the other.
!In this fuction, A is the point on one side of the axis, and B,C are on the opposite sides
  ! Arguments
  real , intent(in) :: Ax,Ay,Bx,By,Cx,Cy
  character , intent(in) :: axis1
  real, intent(out) :: Area_positive, Area_negative
  real :: pABx, pABy, pACx, pACy, A0
  real :: A_half_triangle, A_triangle

  A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);

  call intercept_of_a_line(Ax,Ay,Bx,By,axis1,pABx, pABy);
  call intercept_of_a_line(Ax,Ay,Cx,Cy,axis1,pACx, pACy);

  if (axis1=='x')  A0=Ay; !Value used for if statements (deciding up/down vs left/right)
  if (axis1=='y')  A0=Ax; !Value used for if statements (deciding up/down vs left/right)

  A_half_triangle=Area_of_triangle(Ax,Ay,pABx,pABy,pACx,pACy);
  if (A0>=0.) then
    Area_positive= A_half_triangle;
    Area_negative= A_triangle-A_half_triangle
  else
    Area_positive= A_triangle-A_half_triangle;
    Area_negative= A_half_triangle;
  endif

end subroutine Area_of_triangle_across_axes

subroutine intercept_of_a_line(Ax,Ay,Bx,By,axes1,x0,y0)
!This routine returns the position (x0,y0) at which a line AB intercepts the x or y axis
!The value No_intercept_val is returned when the line does not intercept the axis
  !Arguments
  real, intent(in) :: Ax,Ay,Bx,By
  character, intent(in) ::axes1
  real, intent(out) :: x0,y0
  real :: No_intercept_val !Huge value used to make sure that the intercept is outside the triange in the parralel case.
  

  No_intercept_val=100000000000.; !Huge value used to make sure that the intercept is outside the triange in the parralel case. 
  x0=No_intercept_val
  y0=No_intercept_val
  
  if (axes1=='x') then  !x intercept
    if (Ay.ne.By) then
      x0=Ax -(((Ax-Bx)/(Ay-By))*Ay)
      y0=0.
    endif
  endif

  if (axes1=='y') then !y intercept
    if (Ax.ne.Bx) then
      x0=0.
      y0=-(((Ay-By)/(Ax-Bx))*Ax)+Ay
    endif
  endif
end subroutine intercept_of_a_line


subroutine divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative)
!This routine gives you the area of a triangle on opposite sides of the axis specified. 
!It also takes care of the special case where the triangle is totally on one side
!This routine calls Area_of_triangle_across_axes to calculate the areas when the triangles are split.

  !Arguments
  real, intent(in) :: Ax,Ay,Bx,By,Cx,Cy
  character, intent(in) ::axes1
  real, intent(out) :: Area_positive, Area_negative
  real :: A0,B0,C0
  real A_triangle 

  if (axes1=='x') then  !Use the y-coordinates for if statements to see which side of the line you are on
    A0=Ay
    B0=By
    C0=Cy
  endif
  if (axes1=='y') then  !Use the y-coordinates for if statements to see which side of the line you are on
    A0=Ax
    B0=Bx
    C0=Cx
  endif
  
  A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
  if ((B0*C0)>0.) then !B and C are on the same side  (and non-zero)
    if ((A0*B0).ge.0.) then !all three on the the same side (if it equals zero, then A0=0 and the otehrs are not)
      if ((A0>0.)  .or.  ((A0==0.) .and.  (B0>0.))) then
        Area_positive= A_triangle;
        Area_negative= 0.;
      else
        Area_positive= 0.;
        Area_negative= A_triangle;
      endif
    else  !A is on the opposite side to B and C
      call Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative);
    endif

  elseif ((B0*C0)<0.) then !B and C are on the opposite sides
    if ((A0*B0).ge. 0.) then !C is all alone
      call Area_of_triangle_across_axes(Cx,Cy,Bx,By,Ax,Ay,axes1,Area_positive, Area_negative);
    else !B is all alone
      call Area_of_triangle_across_axes(Bx,By,Cx,Cy,Ax,Ay,axes1,Area_positive, Area_negative);
    endif

  else  !This is the case when either B or C is equal to zero (or both), A0 could be zero too.
    if (((A0.eq.0.) .and. (B0.eq.0.)) .and. (C0.eq.0.)) then
      Area_positive= 0.;
      Area_negative= 0.;
    elseif ((A0*B0<0.)  .or.  (A0*C0<0.)) then    !A, B are on opposite sides, and C is zero.  OR  A, C are on opposite sides, and B is zero.
      call Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative);
    elseif (((A0*B0>0.) .or. (A0*C0>0.)) .or. (((abs(A0)>0.) .and. (B0==0.)) .and. (C0==0.))) then
      if (A0>0.) then
        Area_positive= A_triangle;
        Area_negative= 0.;
      else
        Area_positive= 0.;
        Area_negative= A_triangle;
      endif

    elseif (A0.eq. 0.) then   !(one of B,C is zero too)
      if ((B0>0.) .or. (C0>0.)) then
        Area_positive= A_triangle;
        Area_negative= 0.;
      elseif ((B0<0.) .or. (C0<0.)) then
        Area_positive= 0.;
        Area_negative= A_triangle;
      else
        call error_mesg('diamonds, iceberg_run', 'Logical error inside triangle dividing routine', FATAL)
      endif
    else
      call error_mesg('diamonds, iceberg_run', 'Another logical error inside triangle dividing routine', FATAL)
    endif
  endif
end subroutine divding_triangle_across_axes


subroutine Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy,Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4)
!This routine takes a triangle, and finds the intersection with the four quadrants
  !Arguments
  real, intent(in) :: Ax,Ay,Bx,By,Cx,Cy
  real, intent(out) :: Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4
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

  Area_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
  
  !Calculating area across axes
  call divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'x',Area_Upper ,Area_Lower);
  call divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'y',Area_Right ,Area_Left);

  !Decide if the origin is in the triangle. If so, then you have to divide the area 4 ways
  !This is done by finding a quadrant where the intersection between the triangle and quadrant forms a new triangle
  !(This occurs when on of the sides of the triangle  intersects both the x and y axis)
  if (point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.)) then
    !Find a line in the triangle that cuts both axes in/on the trianlge
    call intercept_of_a_line(Ax,Ay,Bx,By,'x',px,py); !x_intercept
    call intercept_of_a_line(Ax,Ay,Bx,By,'y',qx,qy); !y_intercept
    !Note that the 1. here means that we include points on the boundary of the triange.
    if (.not.((point_in_interval(Ax,Ay,Bx,By,px,py)) .and. (point_in_interval(Ax,Ay,Bx,By,qx,qy)))) then
      call intercept_of_a_line(Ax,Ay,Cx,Cy,'x',px,py); !x_intercept
      call intercept_of_a_line(Ax,Ay,Cx,Cy,'y',qx,qy); !y_intercept
      if (.not.((point_in_interval(Ax,Ay,Cx,Cy,px,py)) .and. (point_in_interval(Ax,Ay,Cx,Cy,qx,qy)))) then
        call intercept_of_a_line(Bx,By,Cx,Cy,'x',px,py); !x_intercept
        call intercept_of_a_line(Bx,By,Cx,Cy,'y',qx,qy); !y_intercept
        if (.not.((point_in_interval(Bx,By,Cx,Cy,px,py)) .and. (point_in_interval(Bx,By,Cx,Cy,qx,qy)))) then
          !You should not get here, but there might be some bugs in the code to do with points exactly falling on axes.
          !if (mpp_pe().eq.12) then
            write(stderrunit,*) 'diamonds,corners', Ax,Ay,Bx,By,Cx,Cy
          !endif
          call error_mesg('diamonds, iceberg_run', 'Something went wrong with Triangle_divide_into_four_quadrants', FATAL)
        endif
      endif
    endif

    !Assigning quadrants. Key_quadrant is the quadrant with the baby triangle in it.
    Area_key_quadrant=Area_of_triangle(px,py,qx,qy,0.,0.)
    if ((px.ge. 0.) .and. (qy.ge. 0.)) then  !First quadrant
      Key_quadrant=1;
    elseif ((px.lt.0.) .and. (qy.ge. 0.)) then  !Second quadrant
      Key_quadrant=2
    elseif ((px.lt. 0.) .and. (qy.lt. 0.)) then !Third quadrant
      Key_quadrant=3;  
    elseif ((px.ge. 0.) .and. (qy.lt. 0.)) then !Forth quadrant
      Key_quadrant=4 
    else  !
      call error_mesg('diamonds, iceberg_run', 'None of the quadrants are Key', WARNING)
      write(stderrunit,*) 'diamonds, Triangle, px,qy', px,qy
    endif
         
  else  !At least one quadrant is empty, and this can be used to find the areas in the other quadrant.  Assigning quadrants. Key_quadrant is the empty quadrant.
    Area_key_quadrant=0;
    if      ( (.not. ((((Ax>0.) .and. (Ay>0.)) .or. ((Bx>0.) .and. (By> 0.))) .or. ((Cx>0.) .and. (Cy> 0.)))) .and. ((Area_Upper+Area_Right).le.Area_triangle) ) then
      !No points land in this quadrant and triangle does not cross the quadrant
      Key_quadrant=1;
    elseif  ( (.not. ((((Ax<0.) .and. (Ay>0)) .or. ((Bx<0.) .and. (By>0.))) .or. ((Cx<0.) .and. (Cy>0.)))) .and. ((Area_Upper+Area_Left).le. Area_triangle) ) then
      Key_quadrant=2
    elseif  ( (.not. ((((Ax<0.) .and. (Ay<0.)) .or. ((Bx<0.) .and. (By< 0.))) .or. ((Cx<0.) .and. (Cy< 0.)))) .and. ((Area_Lower+Area_Left) .le.Area_triangle) ) then
      Key_quadrant=3;
    else
      Key_quadrant=4
    endif
  endif


  !Assign values to quadrants
  if (Key_quadrant .eq. 1) then
    Area_Q1=Area_key_quadrant;
    Area_Q2=Area_Upper-Area_Q1;
    Area_Q4=Area_Right-Area_Q1;
    !Area_Q3=Area_Left-Area_Q2;   !These lines have been changes so that the sum of the 4 quadrants exactly matches the triangle area.
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4);  
  elseif (Key_quadrant .eq. 2) then
    Area_Q2=Area_key_quadrant;
    Area_Q1=Area_Upper-Area_Q2;
    Area_Q4=Area_Right-Area_Q1;
    !Area_Q3=Area_Left-Area_Q2;
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4);
  elseif (Key_quadrant==3) then
    Area_Q3=Area_key_quadrant;
    Area_Q2=Area_Left-Area_Q3;
    Area_Q1=Area_Upper-Area_Q2;
    !Area_Q4=Area_Right-Area_Q1;
    Area_Q4=Area_triangle-(Area_Q1+Area_Q2+Area_Q3);
  elseif (Key_quadrant==4) then
    Area_Q4=Area_key_quadrant;
    Area_Q1=Area_Right-Area_Q4;
    Area_Q2=Area_Upper-Area_Q1;
    !Area_Q3=Area_Left-Area_Q2;
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4);
  else
    call error_mesg('diamonds, iceberg_run', 'Logical error inside triangle into four quadrants. Should not get here.', FATAL)
  endif

  Area_Q1=max(Area_Q1,0.);
  Area_Q2=max(Area_Q2,0.);
  Area_Q3=max(Area_Q3,0.);
  Area_Q4=max(Area_Q4,0.);


  Error=abs(Area_Q1+Area_Q2+Area_Q3+Area_Q4-Area_triangle)
  if (Error>tol) then
    call error_mesg('diamonds, triangle spreading', 'Triangle not evaluated accurately!!', WARNING)
    !if (mpp_pe().eq.mpp_root_pe()) then
    if (mpp_pe().eq. 20) then
      write(stderrunit,*) 'diamonds, Triangle corners:',Ax,Ay,Bx,By,Cx,Cy
      write(stderrunit,*) 'diamonds, Triangle, Full Area', Area_Q1+ Area_Q2+ Area_Q3+ Area_Q4
      write(stderrunit,*) 'diamonds, Triangle, Areas', Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      write(stderrunit,*) 'diamonds, Triangle, Areas', Error
      write(stderrunit,*) 'diamonds, Key quadrant',Key_quadrant,Area_key_quadrant
      write(stderrunit,*) 'diamonds, point in triangle',(point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))
      write(stderrunit,*) 'diamonds, halves',Area_Upper,Area_Lower,Area_Right,Area_Left
    endif
  endif


end subroutine Triangle_divided_into_four_quadrants

subroutine rotate_and_translate(px,py,theta,x0,y0)
  !This function takes a point px,py, and rotates it clockwise around the origin by theta degrees, and then translates by (x0,y0)
  ! Arguments
  real, intent(in) :: x0,y0,theta
  real, intent(inout) :: px,py
  real :: px_temp,py_temp

  !Rotation
  px_temp = ( cos(theta*pi/180)*px) + (sin(theta*pi/180)*py)
  py_temp = (-sin(theta*pi/180)*px) + (cos(theta*pi/180)*py)
 
  !Translation
  px= px_temp + x0
  py= py_temp + y0
end subroutine rotate_and_translate

subroutine Hexagon_into_quadrants_using_triangles(x0,y0,H,theta,Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  !This subroutine divides a regular hexagon centered at x0,y0 with apothen H, and orientation theta into its intersection with the 4 quadrants
  !Theta=0 assumes that the apothen points upwards. (also the rotation is not working yet)
  !Script works by finding the corners of the 6 triangles, and then finding the intersection of each of these with each quadrant.
  !Arguments
  real, intent(in) :: x0,y0,H,theta
  real, intent(out) :: Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4
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

  !Length of side of Hexagon
  S=(2/sqrt(3.))*H
  
  !Finding positions of corners
  C1x=S           ; C1y=0.  !Corner 1 (right)
  C2x=H/sqrt(3.)  ; C2y=H;  !Corner 2 (top right)
  C3x=-H/sqrt(3.) ; C3y=H;  !Corner 3 (top left)
  C4x=-S          ; C4y=0.; !Corner 4 (left)
  C5x=-H/sqrt(3.) ; C5y=-H; !Corner 5 (bottom left)
  C6x=H/sqrt(3.)  ; C6y=-H; !Corner 6 (bottom right)

  !Finding positions of corners
  call rotate_and_translate(C1x,C1y,theta,x0,y0)
  call rotate_and_translate(C2x,C2y,theta,x0,y0)
  call rotate_and_translate(C3x,C3y,theta,x0,y0)
  call rotate_and_translate(C4x,C4y,theta,x0,y0)
  call rotate_and_translate(C5x,C5y,theta,x0,y0)
  call rotate_and_translate(C6x,C6y,theta,x0,y0)

  !Area of Hexagon is the sum of the triangles
  call Triangle_divided_into_four_quadrants(x0,y0,C1x,C1y,C2x,C2y,T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4); !Triangle 012
  call Triangle_divided_into_four_quadrants(x0,y0,C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4); !Triangle 023
  call Triangle_divided_into_four_quadrants(x0,y0,C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4); !Triangle 034
  call Triangle_divided_into_four_quadrants(x0,y0,C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4); !Triangle 045
  call Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4); !Triangle 056
  call Triangle_divided_into_four_quadrants(x0,y0,C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4); !Triangle 061

  !Summing up the triangles
  Area_hex=T12_Area+T23_Area+T34_Area+T45_Area+T56_Area+T61_Area;
  Area_Q1=T12_Q1+T23_Q1+T34_Q1+T45_Q1+T56_Q1+T61_Q1;
  Area_Q2=T12_Q2+T23_Q2+T34_Q2+T45_Q2+T56_Q2+T61_Q2;
  Area_Q3=T12_Q3+T23_Q3+T34_Q3+T45_Q3+T56_Q3+T61_Q3;
  Area_Q4=T12_Q4+T23_Q4+T34_Q4+T45_Q4+T56_Q4+T61_Q4; 
  
  Area_Q1=max(Area_Q1,0.);
  Area_Q2=max(Area_Q2,0.);
  Area_Q3=max(Area_Q3,0.);
  Area_Q4=max(Area_Q4,0.);

  Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
  if ((abs(Error)>tol))then
    if (mpp_pe().eq.mpp_root_pe()) then
      call error_mesg('diamonds, hexagonal spreading', 'Hexagon error is large!!', WARNING)
      write(stderrunit,*) 'diamonds, hex error, H,x0,y0, Error', H, x0 , y0, Error
      write(stderrunit,*) 'diamonds, hex error, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      write(stderrunit,*) 'diamonds, Triangle1',C1x,C1y,C2x,C2y,T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4,(T12_Q1+T12_Q2+T12_Q3+T12_Q4-T12_Area)
      write(stderrunit,*) 'diamonds, Triangle2',C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4,(T23_Q1+T23_Q2+T23_Q3+T23_Q4-T23_Area)
      write(stderrunit,*) 'diamonds, Triangle3',C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4,(T34_Q1+T34_Q2+T34_Q3+T34_Q4-T34_Area)
      write(stderrunit,*) 'diamonds, Triangle4',C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4,(T45_Q1+T45_Q2+T45_Q3+T45_Q4-T45_Area)
      write(stderrunit,*) 'diamonds, Triangle5',C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4,(T56_Q1+T56_Q2+T56_Q3+T56_Q4-T56_Area)
      write(stderrunit,*) 'diamonds, Triangle6',C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4,(T61_Q1+T61_Q2+T61_Q3+T61_Q4-T61_Area)
    endif
  endif

  exact_hex_area=((3.*sqrt(3.)/2)*(S*S))
  if (abs(Area_hex-exact_hex_area)>tol) then
    call error_mesg('diamonds, hexagonal spreading', 'Hexagon not evaluated accurately!!', WARNING)
    if (mpp_pe().eq.mpp_root_pe()) then
      write(stderrunit,*) 'diamonds, hex calculations, H,x0,y0', H, x0 , y0
      write(stderrunit,*) 'diamonds, hex calculations, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
    endif
  endif

  !Adjust Areas so that the error is zero by subtracting the error from the largest sector.
   if  (((Area_Q1>=Area_Q2) .and. (Area_Q1>=Area_Q3)) .and. (Area_Q1>=Area_Q4)) then
     Area_Q1=Area_Q1+Error
   elseif  (((Area_Q2>=Area_Q1) .and. (Area_Q2>=Area_Q3)) .and. (Area_Q2>=Area_Q4)) then
     Area_Q2=Area_Q2+Error
   elseif  (((Area_Q3>=Area_Q1) .and. (Area_Q3>=Area_Q2)) .and. (Area_Q3>=Area_Q4)) then
     Area_Q3=Area_Q3+Error
   elseif  (((Area_Q4>=Area_Q1) .and. (Area_Q4>=Area_Q2)) .and. (Area_Q4>=Area_Q3)) then
     Area_Q4=Area_Q4+Error
   else
     call error_mesg('diamonds, hexagonal spreading', 'Error in hexagon is larger than any quadrant!!', WARNING)
     if (mpp_pe().eq.mpp_root_pe()) then
      write(stderrunit,*) 'diamonds, hex quadrants, H,x0,y0', H, x0 , y0, Error
      write(stderrunit,*) 'diamonds, hex quadrants, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
     endif
   endif

end subroutine Hexagon_into_quadrants_using_triangles


subroutine interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, sss, cn, hi)
! Arguments
type(icebergs_gridded), pointer :: grd
integer, intent(in) :: i, j
real, intent(in) :: xi, yj
real, intent(out) :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, sss, cn, hi
! Local variables
real :: cos_rot, sin_rot
#ifdef USE_OLD_SSH_GRADIENT
real :: dxm, dx0, dxp
#endif
real :: hxm, hxp
real, parameter :: ssh_coast=0.00
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
      call error_mesg('diamonds, interp fields', 'ua is NaNs', FATAL)
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
    write(stderrunit,*) 'diamonds, Error in interpolate: uo,vo,ui,vi',uo, vo, ui, vi
    write(stderrunit,*) 'diamonds, Error in interpolate: ua,va,ssh_x,ssh_y', ua, va, ssh_x, ssh_y
    write(stderrunit,*) 'diamonds, Error in interpolate: sst,cn,hi', sst, sss, cn, hi, mpp_pe()
    call error_mesg('diamonds, interp fields', 'field interpaolations has NaNs', FATAL)

  endif 
  contains

  real function ddx_ssh(grd,i,j)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  ! Local variables
  real :: dxp,dx0
    dxp=0.5*(grd%dx(i+1,j)+grd%dx(i+1,j-1))
    dx0=0.5*(grd%dx(i,j)+grd%dx(i,j-1))
    ddx_ssh=2.*(grd%ssh(i+1,j)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i+1,j)*grd%msk(i,j)
  end function ddx_ssh

  real function ddy_ssh(grd,i,j)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  ! Local variables
  real :: dyp,dy0
    dyp=0.5*(grd%dy(i,j+1)+grd%dy(i-1,j+1))
    dy0=0.5*(grd%dy(i,j)+grd%dy(i-1,j))
    ddy_ssh=2.*(grd%ssh(i,j+1)-grd%ssh(i,j))/(dy0+dyp)*grd%msk(i,j+1)*grd%msk(i,j)
  end function ddy_ssh

  subroutine rotate(u, v, cos_rot, sin_rot)
  ! Arguments
  real, intent(inout) :: u, v
  real, intent(in) :: cos_rot, sin_rot
  ! Local variables
  real :: u_old, v_old

    u_old=u
    v_old=v
    u=cos_rot*u_old+sin_rot*v_old
    v=cos_rot*v_old-sin_rot*u_old

  end subroutine rotate

end subroutine interp_flds


subroutine calculate_mass_on_ocean(bergs, with_diagnostics)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
type(icebergs_gridded), pointer :: grd
logical, intent(in) :: with_diagnostics
! Local variables
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
      i=berg%ine  ;     j=berg%jne
      if (grd%area(i,j) > 0.) then
        
        !Increasing Mass on ocean
        if ((bergs%add_weight_to_ocean .and. .not. bergs%time_average_weight) .or.(bergs%find_melt_using_spread_mass)) then
          call spread_mass_across_ocean_cells(bergs, berg, berg%ine, berg%jne, berg%xi, berg%yj, berg%mass,berg%mass_of_bits, berg%mass_scaling, &
                berg%length*berg%width, berg%thickness)
        endif

        !Calculated some iceberg diagnositcs  
        if (with_diagnostics) call calculate_sum_over_bergs_diagnositcs(bergs,grd,berg,i,j)

      endif
      berg=>berg%next
    enddo
  enddo ;enddo

  contains

  subroutine calculate_sum_over_bergs_diagnositcs(bergs,grd,berg,i,j)
  ! Arguments
  type(icebergs), pointer :: bergs
  type(iceberg), pointer :: berg
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  ! Local variables
  real ::  Abits, Lbits, Mbits

           !Virtual area diagnostic
           if (grd%id_virtual_area>0) then
             if (bergs%bergy_bit_erosion_fraction>0.) then
               Lbits=min(berg%length,berg%width,berg%thickness,40.) ! assume bergy bits are smallest dimension or 40 meters
               Abits=(berg%mass_of_bits/bergs%rho_bergs)/Lbits ! Effective bottom area (assuming T=Lbits)
             else
               Abits=0.0
             endif
             grd%virtual_area(i,j)=grd%virtual_area(i,j)+(berg%width*berg%length+Abits)*berg%mass_scaling ! m^2
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
                & grd%bergy_mass(i,j)=grd%bergy_mass(i,j)+berg%mass_of_bits/grd%area(i,j)*berg%mass_scaling ! kg/m2
 end subroutine calculate_sum_over_bergs_diagnositcs

end subroutine calculate_mass_on_ocean

! ##############################################################################

subroutine icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, calving_hflx, cn, hi, &
                        stagger, stress_stagger, sss, mass_berg, ustar_berg, area_berg)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: time
real, dimension(:,:), intent(inout) :: calving, calving_hflx
real, dimension(:,:), intent(in) :: uo, vo, ui, vi, tauxa, tauya, ssh, sst, cn, hi
integer,    optional, intent(in) :: stagger, stress_stagger
real, dimension(:,:), optional, intent(in) :: sss
real, dimension(:,:), optional, pointer ::  mass_berg, ustar_berg, area_berg
! Local variables
integer :: iyr, imon, iday, ihr, imin, isec, k
type(icebergs_gridded), pointer :: grd
logical :: lerr, sample_traj, write_traj, lbudget, lverbose, check_bond_quality 
real :: unused_calving, tmpsum, grdd_berg_mass, grdd_bergy_mass,grdd_spread_mass, grdd_spread_area
real :: grdd_u_iceberg, grdd_v_iceberg, grdd_ustar_iceberg, grdd_spread_uvel, grdd_spread_vvel
integer :: i, j, Iu, ju, iv, Jv, Iu_off, ju_off, iv_off, Jv_off
real :: mask, max_SST
real, dimension(:,:), allocatable :: uC_tmp, vC_tmp, uA_tmp, vA_tmp
integer :: vel_stagger, str_stagger
real, dimension(:,:), allocatable :: iCount
integer :: nbonds
integer :: stderrunit

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
    if (mod(60*60*24*iday+ 60*60*ihr + 60*imin + isec ,60*60*bergs%traj_sample_hrs).eq.0) &
        sample_traj=.true.
  endif
  write_traj=.false.
  if ((bergs%traj_write_hrs>0) .and. (.not. bergs%ignore_traj))  then
     if (mod(60*60*24*iday+ 60*60*ihr + 60*imin + isec ,60*60*bergs%traj_write_hrs).eq.0) &
         write_traj=.true.
  endif
  lverbose=.false.
  if (bergs%verbose_hrs>0) then
     if (mod(24*iday+ihr+(imin/60.),float(bergs%verbose_hrs)).eq.0) lverbose=verbose
  endif
  lbudget=.false.
  if (bergs%verbose_hrs>0) then
     if (mod(24*iday+ihr+(imin/60.),float(bergs%verbose_hrs)).eq.0) lbudget=budget  !Added minutes, so that it does not repeat when smaller time steps are used.
  endif
  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'diamonds: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
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
    call error_mesg('diamonds, iceberg_run', 'Unrecognized value of stagger!', FATAL)
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
    call error_mesg('diamonds, iceberg_run', 'Unrecognized value of stress_stagger!', FATAL)
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
  
  !Adding gridded salinity.
  if (present(sss)) then
    grd%sss(grd%isc:grd%iec,grd%jsc:grd%jec)=sss(:,:)
  else
    grd%sss(grd%isc:grd%iec,grd%jsc:grd%jec)=-1.0
    if ((bergs%use_mixed_layer_salinity_for_thermo) .and. (bergs%melt_icebergs_as_ice_shelf))  then
      call error_mesg('diamonds, icebergs_run', 'Can not use salinity for thermo. Ocean ML salinity not present!', FATAL)
    endif
  endif

 !Make sure that gridded values agree with mask  (to get ride of NaN values)
  do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
  !Initializing all gridded values to zero
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

  ! For each berg, evolve
  call mpp_clock_begin(bergs%clock_mom)

  if (.not.bergs%Static_icebergs) then
    call evolve_icebergs(bergs)
    if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after evolve()     ')
  endif
  call move_berg_between_cells(bergs)  !Markpoint6
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after move_lists() ')
  if (debug) call bergs_chksum(bergs, 'run bergs (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(bergs%grd, 's/r run after evolve')
  call mpp_clock_end(bergs%clock_mom)

  ! Send bergs to other PEs
  call mpp_clock_begin(bergs%clock_com)
  if (bergs%iceberg_bonds_on)  call  bond_address_update(bergs)

  call send_bergs_to_other_pes(bergs)
  if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after send_bergs() ')
  if ((bergs%interactive_icebergs_on) .or. (bergs%iceberg_bonds_on)) then
    call update_halo_icebergs(bergs)
    if (bergs%debug_iceberg_with_id>0) call monitor_a_berg(bergs, 'icebergs_run, after update_halo()')
    if (bergs%iceberg_bonds_on)  call connect_all_bonds(bergs)
  endif
  if (debug) call bergs_chksum(bergs, 'run bergs (exchanged)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after exchange')
  call mpp_clock_end(bergs%clock_com)

  !Caculate mass on ocean before thermodynamics, to use in melt rate calculation
  if (bergs%find_melt_using_spread_mass) then 
    call calculate_mass_on_ocean(bergs, with_diagnostics=.false.)
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
  if (sample_traj) call record_posn(bergs)
  if (write_traj) then
    call move_all_trajectories(bergs)
    call write_trajectory(bergs%trajectories, bergs%save_short_traj)
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

  !This is the point in the algorithem which determines which fields get passed to the ice model
  !Return what ever calving we did not use and additional icebergs melt
  
  !Making sure that spread_mass has the correct mass
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
    !bergs%spread_mass_end=sum_mass(bergs) !Not sure what this is
    !bergs%spread_area_end=sum_mass(bergs) !Not sure what this is
    !bergs%u_iceberg_end=sum_mass(bergs) !Not sure what this is
    !bergs%v_iceberg_end=sum_mass(bergs) !Not sure what this is
    bergs%floating_heat_end=sum_heat(bergs)
    grd%tmpc(:,:)=0.;
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
    do k=1,nclasses; call mpp_sum(bergs%nbergs_calved_by_class(k)); enddo
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
    grdd_berg_mass=sum( grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_berg_mass)
    grdd_bergy_mass=sum( grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_bergy_mass)
    grdd_spread_mass=sum( grd%spread_mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_spread_mass)
    grdd_spread_area=sum( grd%spread_area(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_spread_area)
    if (mpp_pe().eq.mpp_root_pe()) then
 100 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
 200 format("diamonds: ",a19,10(a18,"=",es14.7,x,a2,:,","))
      call report_state('stored ice','kg','',bergs%stored_start,'',bergs%stored_end,'')
      call report_state('floating','kg','',bergs%floating_mass_start,'',bergs%floating_mass_end,'',bergs%nbergs_end)
      call report_state('icebergs','kg','',bergs%icebergs_mass_start,'',bergs%icebergs_mass_end,'')
      call report_state('bits','kg','',bergs%bergy_mass_start,'',bergs%bergy_mass_end,'')
      call report_state('spread icebergs','kg','',bergs%spread_mass_start,'',bergs%spread_mass_end,'')
      call report_state('spread icebergs','m^2','',bergs%spread_area_start,'',bergs%spread_area_end,'')
      call report_istate('berg #','',bergs%nbergs_start,'',bergs%nbergs_end,'')
      call report_ibudget('berg #','calved',bergs%nbergs_calved, &
                                   'melted',bergs%nbergs_melted, &
                                   '#',bergs%nbergs_start,bergs%nbergs_end)
      call report_budget('stored mass','kg','calving used',bergs%net_calving_used, &
                                            'bergs',bergs%net_calving_to_bergs, &
                                            'stored mass',bergs%stored_start,bergs%stored_end)
      call report_budget('floating mass','kg','calving used',bergs%net_calving_to_bergs, &
                                              'bergs',bergs%net_melt, &
                                              'stored mass',bergs%floating_mass_start,bergs%floating_mass_end)
      call report_budget('berg mass','kg','calving',bergs%net_calving_to_bergs, &
                                          'melt+eros',bergs%berg_melt, &
                                          'berg mass',bergs%icebergs_mass_start,bergs%icebergs_mass_end)
      call report_budget('bits mass','kg','eros used',bergs%bergy_src, &
                                          'bergs',bergs%bergy_melt, &
                                          'stored mass',bergs%bergy_mass_start,bergs%bergy_mass_end)
      call report_budget('net mass','kg','recvd',bergs%net_calving_received, &
                                         'rtrnd',bergs%net_calving_returned, &
                                         'net mass',bergs%stored_start+bergs%floating_mass_start, &
                                                    bergs%stored_end+bergs%floating_mass_end)
      call report_consistant('iceberg mass','kg','gridded',grdd_berg_mass,'bergs',bergs%icebergs_mass_end)
      call report_consistant('spread mass','kg','gridded',grdd_spread_mass,'bergs',bergs%spread_mass_end)
      call report_consistant('spread area','kg','gridded',grdd_spread_area,'bergs',bergs%spread_area_end)
      call report_consistant('bits mass','kg','gridded',grdd_bergy_mass,'bits',bergs%bergy_mass_end)
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
        call report_consistant('top interface','kg','from SIS',bergs%net_incoming_calving,'seen by diamonds',&
             & bergs%net_calving_received)
        call report_consistant('bot interface','kg','sent',bergs%net_outgoing_calving,'seen by SIS',bergs%net_calving_returned)
      endif
      write(*,'("diamonds: calved by class = ",i4,20(",",i4))') (bergs%nbergs_calved_by_class(k),k=1,nclasses)
      if (bergs%nspeeding_tickets>0) write(*,'("diamonds: speeding tickets issued = ",i4)') bergs%nspeeding_tickets
    endif
    bergs%nbergs_start=bergs%nbergs_end
    bergs%stored_start=bergs%stored_end
    bergs%nbergs_melted=0
    bergs%nbergs_calved=0
    bergs%nbergs_calved_by_class(:)=0
    bergs%nspeeding_tickets=0
    bergs%stored_heat_start=bergs%stored_heat_end
    bergs%floating_heat_start=bergs%floating_heat_end
    bergs%floating_mass_start=bergs%floating_mass_end
    bergs%icebergs_mass_start=bergs%icebergs_mass_end
    bergs%bergy_mass_start=bergs%bergy_mass_end
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

  contains

  subroutine report_state(budgetstr,budgetunits,startstr,startval,endstr,endval,delstr,nbergs)
  ! Arguments
  character*(*), intent(in) :: budgetstr, budgetunits, startstr, endstr, delstr
  real, intent(in) :: startval, endval
  integer, intent(in), optional :: nbergs
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
  100 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
  end subroutine report_state

  subroutine report_consistant(budgetstr,budgetunits,startstr,startval,endstr,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr, budgetunits, startstr, endstr
  real, intent(in) :: startval, endval
  ! Local variables
  write(*,200) budgetstr//' check:', &
                      startstr,startval,budgetunits, &
                      endstr,endval,budgetunits, &
                      'error',(endval-startval)/((endval+startval)+1e-30),'nd'
  200 format("diamonds: ",a19,10(a18,"=",es14.7,x,a2,:,","))
  end subroutine report_consistant

  subroutine report_budget(budgetstr,budgetunits,instr,inval,outstr,outval,delstr,startval,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr, budgetunits, instr, outstr, delstr
  real, intent(in) :: inval, outval, startval, endval
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval,budgetunits, &
                      outstr//' out',outval,budgetunits, &
                      'Delta '//delstr,inval-outval,budgetunits, &
                      'error',((endval-startval)-(inval-outval))/max(1.e-30,max(abs(endval-startval),abs(inval-outval))),'nd'
  200 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a8,"=",es10.3,x,a2)
  end subroutine report_budget

  subroutine report_istate(budgetstr,startstr,startval,endstr,endval,delstr)
  ! Arguments
  character*(*), intent(in) :: budgetstr, startstr, endstr, delstr
  integer, intent(in) :: startval, endval
  ! Local variables
  write(*,100) budgetstr//' state:', &
                        startstr//' start',startval, &
                        endstr//' end',endval, &
                        delstr//'Delta',endval-startval
  100 format("diamonds: ",a19,3(a18,"=",i14,x,:,","))
  end subroutine report_istate

  subroutine report_ibudget(budgetstr,instr,inval,outstr,outval,delstr,startval,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr, instr, outstr, delstr
  integer, intent(in) :: inval, outval, startval, endval
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval, &
                      outstr//' out',outval, &
                      'Delta '//delstr,inval-outval, &
                      'error',((endval-startval)-(inval-outval))
  200 format("diamonds: ",a19,10(a18,"=",i14,x,:,","))
  end subroutine report_ibudget

  subroutine get_running_mean_calving(bergs,calving,calving_hflx)
  ! Arguments
  type(icebergs), pointer :: bergs
  real, dimension(:,:), intent(inout) :: calving, calving_hflx
  ! Local variables
  real :: alpha  !Parameter used for calving relaxation time stepping.  (0<=alpha<1)
  real :: tau  !Relaxation timescale in seconds
  real :: beta  ! = 1-alpha (0<=beta<1)
  !This subroutine takes in the new calving and calving_hflx, and uses them to time step a running-mean_calving value
  !The time stepping uses a time scale tau. When tau is equal to zero, the
  !running mean is exactly equal to the new calving value.

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

end subroutine icebergs_run

! ##############################################################################

subroutine icebergs_incr_mass(bergs, mass, Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in), optional :: Time
type(icebergs_gridded), pointer :: grd
integer :: i, j
logical :: lerr
real, dimension(bergs%grd%isc:bergs%grd%iec,bergs%grd%jsc:bergs%grd%jec), intent(inout) :: mass

!This routine is called from SIS, (and older versions of SIS2), but not within
!the iceberg model. The routine adds the spread iceberg mass to mass provided
!the add weight to ocean flag is on, and passive mode is off. It also appears to
!play some role in diagnostics

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


subroutine sum_up_spread_fields(bergs, field, field_name)
! Arguments
type(icebergs), pointer :: bergs
real, dimension(bergs%grd%isc:bergs%grd%iec,bergs%grd%jsc:bergs%grd%jec), intent(out) :: field
character(len=4), intent(in) :: field_name
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

! ##############################################################################

subroutine accumulate_calving(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: remaining_dist, net_calving_used
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
        'diamonds, accumulate_calving: initial stored mass=',bergs%stored_start,' kg'
    do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (grd%calving(i,j).ne.0.) grd%stored_heat(i,j)= & ! Need units of J
            sum(grd%stored_ice(i,j,:)) & ! initial stored ice in kg
           *grd%calving_hflx(i,j)*grd%area(i,j) & ! J/s/m2 x m^2 = J/s
           /grd%calving(i,j) ! /calving in kg/s
    enddo; enddo
    bergs%stored_heat_start=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum( bergs%stored_heat_start )
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,es13.6,a)') &
        'diamonds, accumulate_calving: initial stored heat=',bergs%stored_heat_start,' J'
   endif

  remaining_dist=1.
  do k=1, nclasses
    grd%stored_ice(:,:,k)=grd%stored_ice(:,:,k)+bergs%dt*grd%calving(:,:)*bergs%distribution(k)
    remaining_dist=remaining_dist-bergs%distribution(k)
  enddo
  if (remaining_dist.lt.0.) then
    write(stderrunit,*) 'diamonds, accumulate_calving: sum(distribution)>1!!!',remaining_dist
    call error_mesg('diamonds, accumulate_calving', 'calving is OVER distributed!', WARNING)
  endif
  net_calving_used=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )*(1.-remaining_dist)
  bergs%net_calving_used=bergs%net_calving_used+net_calving_used*bergs%dt
  ! Remove the calving accounted for by accumulation
  grd%calving(:,:)=grd%calving(:,:)*remaining_dist

  ! Do the same for heat (no separate classes needed)
  grd%tmp(:,:)=bergs%dt*grd%calving_hflx(:,:)*grd%area(:,:)*(1.-remaining_dist)
  bergs%net_incoming_calving_heat_used=bergs%net_incoming_calving_heat_used+sum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  grd%stored_heat(:,:)=grd%stored_heat(:,:)+grd%tmp(:,:) ! +=bergs%dt*grd%calving_hflx(:,:)*grd%area(:,:)*(1.-remaining_dist)
  grd%calving_hflx(:,:)=grd%calving_hflx(:,:)*remaining_dist

end subroutine accumulate_calving

! ##############################################################################

subroutine calve_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
integer :: i,j,k,icnt,icntmax
integer :: iNg, jNg  !Total number of points gloablly in i and j direction
type(iceberg) :: newberg
logical :: lret
real :: xi, yj, ddt, calving_to_bergs, calved_to_berg, heat_to_bergs, heat_to_berg
integer :: stderrunit

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

  do k=1, nclasses
    do j=grd%jsc, grd%jec
      do i=grd%isc, grd%iec
        ddt=0.; icnt=0
        do while (grd%stored_ice(i,j,k).ge.bergs%initial_mass(k)*bergs%mass_scaling(k))
          newberg%lon=0.25*((grd%lon(i,j)+grd%lon(i-1,j-1))+(grd%lon(i-1,j)+grd%lon(i,j-1)))
          newberg%lat=0.25*((grd%lat(i,j)+grd%lat(i-1,j-1))+(grd%lat(i-1,j)+grd%lat(i,j-1)))
         !write(stderr(),*) 'diamonds, calve_icebergs: creating new iceberg at ',newberg%lon,newberg%lat
          lret=pos_within_cell(grd, newberg%lon, newberg%lat, i, j, xi, yj)
          if (.not.lret) then
            write(stderrunit,*) 'diamonds, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('diamonds, calve_icebergs', 'berg is not in the correct cell!', FATAL)
          endif
          if (debug.and.(xi<0..or.xi>1..or.yj<0..or.yj>1.)) then
            write(stderrunit,*) 'diamonds, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('diamonds, calve_icebergs', 'berg xi,yj is not correct!', FATAL)
          endif
          if (grd%msk(i,j)<0.5) then
            write(stderrunit,*) 'diamonds, calve_icebergs: WARNING!!! Iceberg born in land cell',i,j,newberg%lon,newberg%lat
            if (debug) call error_mesg('diamonds, calve_icebergs', 'Iceberg born in Land Cell!', FATAL)
          endif
          newberg%ine=i
          newberg%jne=j
          newberg%xi=xi
          newberg%yj=yj
          newberg%uvel=0.
          newberg%vvel=0.
          newberg%axn=0. !Added by Alon
          newberg%ayn=0. !Added by Alon
          newberg%bxn=0. !Added by Alon
          newberg%byn=0. !Added by Alon
          newberg%mass=bergs%initial_mass(k)
          newberg%thickness=bergs%initial_thickness(k)
          newberg%width=bergs%initial_width(k)
          newberg%length=bergs%initial_length(k)
          newberg%start_lon=newberg%lon
          newberg%start_lat=newberg%lat
          newberg%start_year=bergs%current_year
          newberg%iceberg_num=((iNg*jNg)*grd%iceberg_counter_grd(i,j))+(i+(iNg*(j-1)))  ! unique number for each iceberg
          newberg%start_day=bergs%current_yearday+ddt/86400.
          newberg%start_mass=bergs%initial_mass(k)
          newberg%mass_scaling=bergs%mass_scaling(k)
          newberg%mass_of_bits=0.
          newberg%halo_berg=0.
          newberg%static_berg=0.
          newberg%heat_density=grd%stored_heat(i,j)/grd%stored_ice(i,j,k) ! This is in J/kg
          call add_new_berg_to_list(bergs%list(i,j)%first, newberg)
          calved_to_berg=bergs%initial_mass(k)*bergs%mass_scaling(k) ! Units of kg
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
          bergs%nbergs_calved_by_class(k)=bergs%nbergs_calved_by_class(k)+1
          grd%iceberg_counter_grd(i,j)=grd%iceberg_counter_grd(i,j)+1
        enddo
        icntmax=max(icntmax,icnt)
      enddo
    enddo
  enddo

  if (debug.and.icntmax>1) write(stderrunit,*) 'calve_icebergs: icnt=',icnt,' on',mpp_pe()

  bergs%net_calving_to_bergs=bergs%net_calving_to_bergs+calving_to_bergs
  bergs%net_heat_to_bergs=bergs%net_heat_to_bergs+heat_to_bergs

end subroutine calve_icebergs

! ##############################################################################

subroutine evolve_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
type(iceberg), pointer :: berg
real :: uveln, vveln, lonn, latn  
real :: axn, ayn, bxn, byn                                           ! Added by Alon - explicit and implicit accelations from the previous step
real :: xi, yj
integer :: i, j
integer :: grdi, grdj
integer :: stderrunit
logical :: bounced, interactive_icebergs_on, Runge_not_Verlet 
  
  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  interactive_icebergs_on=bergs%interactive_icebergs_on
  Runge_not_Verlet=bergs%Runge_not_Verlet


  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    berg=>bergs%list(grdi,grdj)%first
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
          if (debug) call error_mesg('diamonds, evolve_iceberg','berg is in wrong starting cell!',FATAL)
        endif
        if (debug) call check_position(grd, berg, 'evolve_iceberg (top)')

          !Time stepping schemes:
          if (Runge_not_Verlet) then 
            call Runge_Kutta_stepping(bergs,berg, axn, ayn, bxn, byn, uveln, vveln,lonn, latn, i, j, xi, yj)
          endif 
          if (.not.Runge_not_Verlet) then 
            call verlet_stepping(bergs,berg, axn, ayn, bxn, byn, uveln, vveln)
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
          if (.not. interactive_icebergs_on)  call update_verlet_position(bergs,berg) 
        endif

        !call interp_flds(grd, i, j, xi, yj, berg%uo, berg%vo, berg%ui, berg%vi, berg%ua, berg%va, berg%ssh_x, berg%ssh_y, berg%sst)
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
         if (.not. Runge_not_Verlet)  call update_verlet_position(bergs,berg) 
          
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

!contains
end subroutine evolve_icebergs

!######################################################################

subroutine verlet_stepping(bergs,berg, axn, ayn, bxn, byn, uveln, vveln)
type(icebergs), pointer :: bergs
type(iceberg), pointer, intent(inout) :: berg
type(icebergs_gridded), pointer :: grd
! Local variables
real, intent(out) :: axn, ayn, bxn, byn, uveln, vveln
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

  !Initialize variables

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
        i=berg%ine     ;   j=berg%jne
        xi=berg%xi      ;   yj=berg%yj

        ! Turn the velocities into u_star, v_star.(uvel3 is v_star) - Alon (not sure how this works with tangent plane)
        uvel3=uvel1+(dt_2*axn)                  !Alon
        vvel3=vvel1+(dt_2*ayn)                  !Alon

        !Note, the mass scaling is equal to 1 (rather than 0.25 as in RK), since
        !this is only called once in Verlet stepping.
        if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
          call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 1.0*berg%mass_scaling,berg%length*berg%width, berg%thickness)

        ! Calling the acceleration   (note that the velocity is converted to u_star inside the accel script)
        call accel(bergs, berg, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, ax1, ay1, axn, ayn, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  
        !Solving for the new velocity
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

        !if (berg%iceberg_num .eq. 1) print *, 'New velocity: ', uveln, vveln

       
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
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lonn=',lonn,berg%lon
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: latn=',latn,berg%lat
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u3,un,u0=',uvel3,uveln,berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v3,vn,v0=',vvel3,vveln,berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: iceberg_num=',berg%iceberg_num
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1=',&
              & dt*ax1
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1=',&
              & dt*ay1
         write(stderrunit,*) 'diamonds, evolve_iceberg: on_tangential_plane=',on_tangential_plane
         write(stderrunit,*) 'Acceleration terms for position 1'
         error_flag=pos_within_cell(grd, lonn, latn, i, j, xi,  yj)
         call accel(bergs, berg, i, j, xi, yj, latn, uvel3, vvel3, uvel1, vvel1, dt_2, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn, bxn, byn - Added by Alon
         
          write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
          write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
          call print_berg(stderrunit, berg, 'evolve_iceberg, out of cell at end!')
          bounced=is_point_in_cell(bergs%grd, lonn, latn, i, j,explain=.true.)
          if (debug) call error_mesg('diamonds, evolve_iceberg','berg is out of posn at end!',FATAL)
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

!######################################################################

subroutine Runge_Kutta_stepping(bergs, berg, axn, ayn, bxn, byn, uveln, vveln,lonn, latn, i, j, xi, yj)
type(icebergs), pointer :: bergs
type(iceberg), pointer, intent(inout) :: berg
type(icebergs_gridded), pointer :: grd
real , intent(out) :: axn, ayn, bxn, byn, uveln, vveln,lonn, latn, xi, yj
integer, intent(out) :: i, j
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
          call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness)

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

        call accel(bergs, berg, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn) !axn,ayn, bxn, byn  - Added by Alon
        !call accel(bergs, berg, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt, ax1, ay1, axn1, ayn1, bxn, byn) !Note change to dt. Markpoint_1
        if (on_tangential_plane) call rotvec_to_tang(lon1,ax1,ay1,xddot1,yddot1)
        if (on_tangential_plane) call rotvec_to_tang(lon1,axn1,ayn1,xddot1n,yddot1n) !Alon
        
        !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
       !if (debug) write(stderr(),*) 'diamonds, evolve: x2=...'
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
        call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag, berg%iceberg_num)
        i2=i; j2=j
        if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
          call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness)
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
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i=',i1,i2,i
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j=',j1,j2,j
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2=',lon1,lon2,berg%lon
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2=',lat1,lat2,berg%lat
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u0=',uvel1,uvel2,berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v0=',vvel1,vvel2,berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2=',dt*ax1,dt*ax2
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2=',dt*ay1,dt*ay2
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u0=',dt*uvel1,dt*uvel2,dt*berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v0=',dt*vvel1,dt*vvel2,dt*berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2 (deg)=',dt*u1,dt*u2
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2 (deg)=',dt*v1,dt*v2
         write(stderrunit,*) 'Acceleration terms for position 1'
         error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
         call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn,- Added by Alon
          call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 2')
          write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos2 i,j,lon,lat,xi,yj=',i,j,lon2,lat2,xi,yj
          write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos2 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
          bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j,explain=.true.)
          call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 2!',FATAL)
        endif
        call  convert_from_meters_to_grid(lat2,bergs%grd%grid_is_latlon ,dxdl2,dydl)
        !dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
        u2=uvel2*dxdl2; v2=vvel2*dydl
        call accel(bergs, berg, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn2, ayn2, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
        !call accel(bergs, berg, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt, ax2, ay2, axn2, ayn2, bxn, byn) !Note change to dt. Markpoint_1
        if (on_tangential_plane) call rotvec_to_tang(lon2,ax2,ay2,xddot2,yddot2)
        if (on_tangential_plane) call rotvec_to_tang(lon2,axn2,ayn2,xddot2n,yddot2n) !Alon
        
        !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
       !if (debug) write(stderr(),*) 'diamonds, evolve: x3=...'
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
        call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag, berg%iceberg_num)
        i3=i; j3=j
        if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
          call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness)
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
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i3,i=',i1,i2,i3,i
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j3,j=',j1,j2,j3,j
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2,lon3=',lon1,lon2,lon3,berg%lon
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2,lat3=',lat1,lat2,lat3,berg%lat
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u3,u0=',uvel1,uvel2,uvel3,berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v3,v0=',vvel1,vvel2,vvel3,berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2,ax3=',dt*ax1,dt*ax2,dt*ax3
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2,ay3=',dt*ay1,dt*ay2,dt*ay3
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3 (deg)=',dt*u1,dt*u2,dt*u3
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3 (deg)=',dt*v1,dt*v2,dt*v3
         write(stderrunit,*) 'Acceleration terms for position 1'
         error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
         call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn - Added by Alon
         write(stderrunit,*) 'Acceleration terms for position 2'
         error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
         call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn2, ayn2, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
          call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 3')
          write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos3 i,j,lon,lat,xi,yj=',i,j,lon3,lat3,xi,yj
          write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos3 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
          bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j,explain=.true.)
          call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 3!',FATAL)
        endif
        call  convert_from_meters_to_grid(lat3,bergs%grd%grid_is_latlon ,dxdl3,dydl)
        !dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
        u3=uvel3*dxdl3; v3=vvel3*dydl
        call accel(bergs, berg, i, j, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, axn3, ayn3, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
        if (on_tangential_plane) call rotvec_to_tang(lon3,ax3,ay3,xddot3,yddot3)
        if (on_tangential_plane) call rotvec_to_tang(lon3,axn3,ayn3,xddot3n,yddot3n) !Alon
        
        !  X4 = X1+dt*V3 ; V4 = V1+dt*A3; A4=A(X4)
       !if (debug) write(stderr(),*) 'diamonds, evolve: x4=...'
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
        call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag, berg%iceberg_num)
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
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2,lon3,lon4=',lon1,lon2,lon3,lon4,berg%lon
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2,lat3,lat4=',lat1,lat2,lat3,lat4,berg%lat
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u3,u4,u0=',uvel1,uvel2,uvel3,uvel4,berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v3,v4,v0=',vvel1,vvel2,vvel3,vvel4,berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2,ax3,ax4=',dt*ax1,dt*ax2,dt*ax3,dt*ax4
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2,ay3,ay4=',dt*ay1,dt*ay2,dt*ay3,dt*ay4
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4 (deg)=',dt*u1,dt*u2,dt*u3,dt*u4
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4 (deg)=',dt*v1,dt*v2,dt*v3,dt*v4
         write(stderrunit,*) 'Acceleration terms for position 1'
         error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
         call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
         write(stderrunit,*) 'Acceleration terms for position 2'
         error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
         call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn2, ayn2, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
         write(stderrunit,*) 'Acceleration terms for position 3'
         error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
         call accel(bergs, berg, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, axn3, ayn3, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
          call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 4')
          write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos4 i,j,lon,lat,xi,yj=',i,j,lon4,lat4,xi,yj
          write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos4 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
          bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
          call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 4!',FATAL)
        endif
        call  convert_from_meters_to_grid(lat4,bergs%grd%grid_is_latlon ,dxdl4,dydl)
        !dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
        u4=uvel4*dxdl4; v4=vvel4*dydl
        call accel(bergs, berg, i, j, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, ax4, ay4, axn4, ayn4, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
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
        call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag, berg%iceberg_num)
        if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
          call spread_mass_across_ocean_cells(bergs, berg, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling,berg%length*berg%width, berg%thickness)
  
        if (.not.error_flag) then
          if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j)) error_flag=.true.
        endif
        if (error_flag) then
         call print_fld(grd, grd%msk, 'msk')
         call print_fld(grd, grd%ssh, 'ssh')
         call print_fld(grd, grd%sst, 'sst')
         call print_fld(grd, grd%sss, 'sss')
         call print_fld(grd, grd%hi, 'hi')
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
         write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2,lon3,lon4,lonn=',lon1,lon2,lon3,lon4,lonn,berg%lon
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2,lat3,lat4,latn=',lat1,lat2,lat3,lat4,latn,berg%lat
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u3,u4,un,u0=',uvel1,uvel2,uvel3,uvel4,uveln,berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v3,v4,vn,v0=',vvel1,vvel2,vvel3,vvel4,vveln,berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2,ax3,ax4,axn=',&
              & dt*ax1,dt*ax2,dt*ax3,dt*ax4,dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2,ay3,ay4,ayn=',&
              & dt*ay1,dt*ay2,dt*ay3,dt*ay4,dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4,un,u0=',&
              & dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*uveln,dt*berg%uvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4,vn,v0=',&
              & dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*vveln,dt*berg%vvel
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4,u_rk (deg)=',&
              & dt*u1,dt*u2,dt*u3,dt*u4,dt_6*( (u1+u4)+2.*(u2+u3) )
         write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4,v_rk (deg)=',&
              & dt*v1,dt*v2,dt*v3,dt*v4,dt_6*( (v1+v4)+2.*(v2+v3) )
         write(stderrunit,*) 'diamonds, evolve_iceberg: on_tangential_plane=',on_tangential_plane
         write(stderrunit,*) 'Acceleration terms for position 1'
         error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
         call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
         write(stderrunit,*) 'Acceleration terms for position 2'
         error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
         call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
         write(stderrunit,*) 'Acceleration terms for position 3'
         error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
         call accel(bergs, berg, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
         write(stderrunit,*) 'Acceleration terms for position 4'
         error_flag=pos_within_cell(grd, lon4, lat4, i4, j4, xi, yj)
         call accel(bergs, berg, i4, j4, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, ax4, ay4, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
          write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
          write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
          call print_berg(stderrunit, berg, 'evolve_iceberg, out of cell at end!')
          bounced=is_point_in_cell(bergs%grd, lonn, latn, i, j, explain=.true.)
          if (debug) call error_mesg('diamonds, evolve_iceberg','berg is out of posn at end!',FATAL)
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


!#######################################################################
!MP6
subroutine update_verlet_position(bergs,berg)
type(icebergs), intent(in), pointer :: bergs 
type(iceberg), intent(in), pointer :: berg 
type(icebergs_gridded), pointer :: grd
!Local variable
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
        if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1,berg%iceberg_num)
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
        call adjust_index_and_ground(grd, lonn, latn, uvel3, vvel3, i, j, xi, yj, bounced, error_flag, berg%iceberg_num)  !Alon:"unclear which velocity to use here?"

        !if (bounced) then
        !  print *, 'you have been bounce: big time!',mpp_pe(),berg%iceberg_num,lonn, latn, uvel3, vvel3, i, j, xi, yj, bounced, error_flag 
        !  berg%axn=0.0  ;  berg%ayn=0.0
        !  berg%bxn=0.0  ;  berg%byn=0.0
        !  berg%uvel=0.0 ;  berg%vvel=0.0
        !endif

        !Updating positions and index
        berg%lon=lonn      ;  berg%lat=latn
        berg%ine=i    ;  berg%jne=j
        berg%xi=xi    ;  berg%yj=yj

end subroutine update_verlet_position

!#######################################################################

  subroutine rotpos_from_tang(x, y, lon, lat)
  ! Arguments
  real, intent(in) :: x, y
  real, intent(out) :: lon, lat
  ! Local variables
  real :: r

    r=sqrt(x**2+y**2)
    lat=90.-(r180_pi*r/Rearth)
    lon=r180_pi*acos(x/r)*sign(1.,y)

  end subroutine rotpos_from_tang

  subroutine rotvec_to_tang(lon, uvel, vvel, xdot, ydot)
  ! Arguments
  real, intent(in) :: lon, uvel, vvel
  real, intent(out) :: xdot, ydot
  ! Local variables
  real :: clon,slon

    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    xdot=-slon*uvel-clon*vvel
    ydot=clon*uvel-slon*vvel

  end subroutine rotvec_to_tang

  subroutine rotvec_from_tang(lon, xdot, ydot, uvel, vvel)
  ! Arguments
  real, intent(in) :: lon, xdot, ydot
  real, intent(out) :: uvel, vvel
  ! Local variables
  real :: clon,slon

    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    uvel=-slon*xdot+clon*ydot
    vvel=-clon*xdot-slon*ydot

  end subroutine rotvec_from_tang

! ##############################################################################

subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error, iceberg_num)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(inout) :: lon, lat, uvel, vvel, xi, yj
integer, intent(inout) :: i,j
integer, intent(in) :: iceberg_num
logical, intent(out) :: bounced, error
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
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','Iceberg is_point_in_cell=True but xi, yi are out of cell',FATAL)
        error=.true.; return
     endif
    else
      if (.not.lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=F but |xi,yj|<1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,  explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
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
    write(stderrunit,*) 'pe=',mpp_pe(),'diamonds, adjust: inm,i0,inm-i0=',inm,i0,inm-i0
   !stop 'Moved too far in i without mask!'
  endif
  if (abs(jnm-j0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'diamonds, adjust: jnm,i0,jnm-j0=',jnm,j0,inm-j0
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
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from west',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (xi.ge.1.) then    !Alon!!!!
!    elseif (xi.gt.1.) then
      if (i<grd%ied) then
        if (grd%msk(i+1,j)>0.) then
          if (i<grd%ied) i=i+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from east',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (yj.lt.0.) then
      if (j>grd%jsd) then
        if (grd%msk(i,j-1)>0.) then
          if (j>grd%jsd+1) j=j-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from south',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (yj.ge.1.) then     !Alon.
!    elseif (yj.gt.1.) then
      if (j<grd%jed) then
        if (grd%msk(i,j+1)>0.) then
          if (j<grd%jed) j=j+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from north',lon,lat,xi,yj,uvel,vvel,mpp_pe()
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
      if (grd%msk(i,j)==0.) stop 'diamonds, adjust: Berg is in land! This should not happen...'
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  enddo
 !if (debug) then
 !  if (abs(i-i0)>2) then
 !    stop 'diamonds, adjust: Moved too far in i!'
 !  endif
 !  if (abs(j-j0)>2) then
 !    stop 'diamonds, adjust: Moved too far in j!'
 !  endif
 !endif

  if (.not.bounced.and.lret.and.grd%msk(i,j)>0.) return ! Landed in ocean without bouncing so all is well

  if (.not.bounced.and..not.lret) then ! This implies the berg traveled many cells without getting far enough
    if (debug) then
      write(stderrunit,*) 'diamonds, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'diamonds, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'diamonds, adjust: i0,j0=',i0,j0 
      write(stderrunit,*) 'diamonds, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'diamonds, adjust: xi,yj=',xi,yj 
      write(stderrunit,*) 'diamonds, adjust: i,j=',i,j
      write(stderrunit,*) 'diamonds, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'diamonds, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
      write(stderrunit,*) 'diamonds, adjust: lret=',lret
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
      call error_mesg('diamonds, adjust', 'Berg did not move or bounce during iterations AND was not in cell. Adjusting!', WARNING)
      write(stderrunit,*) 'diamonds, adjust: The adjusting iceberg is: ', iceberg_num,  mpp_pe()
      write(stderrunit,*) 'diamonds, adjust: The adjusting lon,lat,u,v: ', lon, lat, uvel, vvel
      write(stderrunit,*) 'diamonds, adjust: The adjusting xi,ji: ', xi, yj
      lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj,explain=.true.)
    else
      call error_mesg('diamonds, adjust', 'Berg iterated many times without bouncing!', WARNING)
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
    write(stderrunit,*) 'diamonds, adjust: Should not get here! Berg is not in cell after adjustment', iceberg_num, mpp_pe()
    if (debug) error=.true.
  endif
 end subroutine adjust_index_and_ground

!end subroutine evolve_icebergs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rotpos_to_tang(lon, lat, x, y, iceberg_num_in)
  ! Arguments
  real, intent(in) :: lon, lat
  real, intent(out) :: x, y
  integer, intent(in) , optional :: iceberg_num_in
  ! Local variables
  real :: r,colat,clon,slon
  integer :: stderrunit, iceberg_num

  stderrunit = stderr()
  iceberg_num=000
  if (present(iceberg_num_in)) then
        iceberg_num=iceberg_num_in
  endif
        
  if (lat>90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat>90 already!',lat, lon, iceberg_num
      call error_mesg('diamonds, rotpos_to_tang','Something went very wrong!',FATAL)
    endif
    if (lat==90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat==90 already!',lat, lon
      call error_mesg('diamonds, rotpos_to_tang','Something went wrong!',FATAL)
    endif

    colat=90.-lat
    r=Rearth*(colat*pi_180)
    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    x=r*clon
    y=r*slon

  end subroutine rotpos_to_tang

! ##############################################################################

subroutine icebergs_stock_pe(bergs, index, value)
! Modules
use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
! Arguments
type(icebergs), pointer :: bergs
integer, intent(in) :: index
real, intent(out) :: value
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

! ##############################################################################

subroutine icebergs_save_restart(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables

  if (.not.associated(bergs)) return

  call mpp_clock_begin(bergs%clock_iow)
  call bergs_chksum(bergs, 'write_restart bergs')
  call write_restart(bergs)
  call mpp_clock_end(bergs%clock_iow)

end subroutine icebergs_save_restart

! ##############################################################################

subroutine icebergs_end(bergs)
! Arguments
type(icebergs), pointer :: bergs
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

  if (.not. bergs%ignore_traj) &
  call write_trajectory(bergs%trajectories, bergs%save_short_traj)

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
  deallocate(bergs%grd%domain)
  deallocate(bergs%grd)
  deallocate(bergs%initial_mass)
  deallocate(bergs%distribution)
  deallocate(bergs%initial_thickness)
  deallocate(bergs%initial_width)
  deallocate(bergs%initial_length)
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

  if (mpp_pe()==mpp_root_pe()) write(*,'(a,i8)') 'diamonds: icebergs_end complete',mpp_pe()

  contains

  subroutine dealloc_buffer(buff)
  ! Arguments
  type(buffer), pointer :: buff
  ! Local variables
    if (associated(buff)) then
      if (associated(buff%data)) deallocate(buff%data)
      deallocate(buff)
    endif
  end subroutine dealloc_buffer

end subroutine icebergs_end

! ##############################################################################

subroutine invert_tau_for_du(u,v)
! Arguments
real, dimension(:,:),intent(inout) :: u, v
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

! ##############################################################################


end module
