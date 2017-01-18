module ice_bergs_framework

use constants_mod, only: radius, pi, omega, HLF

use mpp_domains_mod, only: domain2D

use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_pe, mpp_root_pe, mpp_chksum
use mpp_mod, only: COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
use mpp_mod, only: COMM_TAG_5, COMM_TAG_6, COMM_TAG_7, COMM_TAG_8
use mpp_mod, only: COMM_TAG_9, COMM_TAG_10
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING

use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)

implicit none ; private

integer :: buffer_width=27 !Changed from 20 to 28 by Alon 
integer :: buffer_width_traj=31  !Changed from 23 by Alon
!integer, parameter :: buffer_width=26 !Changed from 20 to 26 by Alon 
!integer, parameter :: buffer_width_traj=29  !Changed from 23 by Alon
integer, parameter :: nclasses=10 ! Number of ice bergs classes

!Local Vars
! Global data (minimal for debugging)
logical :: folded_north_on_pe = .false.
logical :: verbose=.false. ! Be verbose to stderr
logical :: budget=.true. ! Calculate budgets
logical :: debug=.false. ! Turn on debugging
logical :: really_debug=.false. ! Turn on debugging
logical :: parallel_reprod=.true. ! Reproduce across different PE decompositions
logical :: use_slow_find=.true. ! Use really slow (but robust) find_cell for reading restarts
logical :: ignore_ij_restart=.false. ! Read i,j location from restart if available (needed to use restarts on different grids)
logical :: generate_test_icebergs=.false. ! Create icebergs in absence of a restart file
logical :: use_roundoff_fix=.true. ! Use a "fix" for the round-off discrepancy between is_point_in_cell() and pos_within_cell()
logical :: old_bug_rotated_weights=.false. ! Skip the rotation of off-center weights for rotated halo updates
logical :: make_calving_reproduce=.false. ! Make the calving.res.nc file reproduce across pe count changes.
logical :: old_bug_bilin=.true. ! If true, uses the inverted bilin function (use False to get correct answer)
character(len=10) :: restart_input_dir = 'INPUT/'
integer, parameter :: delta_buf=25 ! Size by which to increment buffers
real, parameter :: pi_180=pi/180. ! Converts degrees to radians
logical :: fix_restart_dates=.true. ! After a restart, check that bergs were created before the current model date
logical :: do_unit_tests=.false. ! Conduct some unit tests
logical :: force_all_pes_traj=.false. ! Force all pes write trajectory files regardless of io_layout

!Public params !Niki: write a subroutine to expose these
public nclasses,buffer_width,buffer_width_traj
public verbose, really_debug, debug, restart_input_dir,make_calving_reproduce,old_bug_bilin,use_roundoff_fix
public ignore_ij_restart, use_slow_find,generate_test_icebergs,old_bug_rotated_weights,budget
public orig_read, force_all_pes_traj

!Public types
public icebergs_gridded, xyt, iceberg, icebergs, buffer, bond 

!Public subs
public ice_bergs_framework_init
public send_bergs_to_other_pes
public update_halo_icebergs
public pack_berg_into_buffer2, unpack_berg_from_buffer2
public pack_traj_into_buffer2, unpack_traj_from_buffer2
public increase_ibuffer 
public add_new_berg_to_list, count_out_of_order, check_for_duplicates
public insert_berg_into_list, create_iceberg, delete_iceberg_from_list, destroy_iceberg
public print_fld,print_berg, print_bergs,record_posn, push_posn, append_posn, check_position
public move_trajectory, move_all_trajectories
public form_a_bond, connect_all_bonds, show_all_bonds, bond_address_update
public find_cell, find_cell_by_search, count_bergs, is_point_in_cell, pos_within_cell, count_bonds
public sum_mass, sum_heat, bilin, yearday, bergs_chksum, list_chksum, count_bergs_in_list
public checksum_gridded
public grd_chksum2,grd_chksum3
public fix_restart_dates, offset_berg_dates
public move_berg_between_cells
public find_individual_iceberg
public monitor_a_berg
public is_point_within_xi_yj_bounds
public test_check_for_duplicate_ids_in_list
public check_for_duplicates_in_parallel

type :: icebergs_gridded
  type(domain2D), pointer :: domain ! MPP domain
  integer :: halo ! Nominal halo width
  integer :: isc, iec, jsc, jec ! Indices of computational domain
  integer :: isd, ied, jsd, jed ! Indices of data domain
  integer :: isg, ieg, jsg, jeg ! Indices of global domain
  integer :: my_pe, pe_N, pe_S, pe_E, pe_W ! MPI PE idenLx ! Length of domain, for periodic boundary condition (Ly to be adde later if needed)
  logical :: grid_is_latlon !Flag to say whether the coordinate is in lat lon degrees, or meters
  logical :: grid_is_regular !Flag to say whether point in cell can be found assuming regular cartesian grid
  real :: Lx !Length of the domain in x direction
  real, dimension(:,:), pointer :: lon=>null() ! Longitude of cell corners
  real, dimension(:,:), pointer :: lat=>null() ! Latitude of cell corners
  real, dimension(:,:), pointer :: lonc=>null() ! Longitude of cell centers
  real, dimension(:,:), pointer :: latc=>null() ! Latitude of cell centers
  real, dimension(:,:), pointer :: dx=>null() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: dy=>null() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: area=>null() ! Area of cell (m^2)
  real, dimension(:,:), pointer :: msk=>null() ! Ocean-land mask (1=ocean)
  real, dimension(:,:), pointer :: cos=>null() ! Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: sin=>null() ! Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: ocean_depth=>NULL() ! Depth of ocean (m)
  real, dimension(:,:), pointer :: uo=>null() ! Ocean zonal flow (m/s)
  real, dimension(:,:), pointer :: vo=>null() ! Ocean meridional flow (m/s)
  real, dimension(:,:), pointer :: ui=>null() ! Ice zonal flow (m/s)
  real, dimension(:,:), pointer :: vi=>null() ! Ice meridional flow (m/s)
  real, dimension(:,:), pointer :: ua=>null() ! Atmosphere zonal flow (m/s)
  real, dimension(:,:), pointer :: va=>null() ! Atmosphere meridional flow (m/s)
  real, dimension(:,:), pointer :: ssh=>null() ! Sea surface height (m)
  real, dimension(:,:), pointer :: sst=>null() ! Sea surface temperature (oC)
  real, dimension(:,:), pointer :: sss=>null() ! Sea surface salinity (psu)
  real, dimension(:,:), pointer :: cn=>null() ! Sea-ice concentration (0 to 1)
  real, dimension(:,:), pointer :: hi=>null() ! Sea-ice thickness (m)
  real, dimension(:,:), pointer :: calving=>null() ! Calving mass rate [frozen runoff] (kg/s) (into stored ice)
  real, dimension(:,:), pointer :: calving_hflx=>null() ! Calving heat flux [heat content of calving] (W/m2) (into stored ice)
  real, dimension(:,:), pointer :: floating_melt=>null() ! Net melting rate to icebergs + bits (kg/s/m^2)
  real, dimension(:,:), pointer :: berg_melt=>null() ! Melting+erosion rate of icebergs (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_buoy=>null() ! Buoyancy componenet of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_eros=>null() ! Erosion component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_conv=>null() ! Convective component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_src=>null() ! Mass flux from berg erosion into bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_melt=>null() ! Melting rate of bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_mass=>null() ! Mass distribution of bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: spread_mass=>null() ! Mass of icebergs after spreading (kg/m^2)
  real, dimension(:,:), pointer :: spread_mass_old=>null() ! Mass of icebergs after spreading old (kg/m^2)
  real, dimension(:,:), pointer :: spread_area=>null() ! Area of icebergs after spreading (m^2/m^2)
  real, dimension(:,:), pointer :: u_iceberg=>null() ! Average iceberg velocity in grid cell (mass weighted - but not spread mass weighted)
  real, dimension(:,:), pointer :: v_iceberg=>null() ! Average iceberg velocity in grid cell (mass weighted - but not spread mass weighted)
  real, dimension(:,:), pointer :: spread_uvel=>null() ! Average iceberg velocity in grid cell (spread area weighted)
  real, dimension(:,:), pointer :: spread_vvel=>null() ! Average iceberg velocity in grid cell (spread area weighted)
  real, dimension(:,:), pointer :: ustar_iceberg=>null() ! Frictional velocity below icebergs to be passed to ocean
  real, dimension(:,:), pointer :: virtual_area=>null() ! Virtual surface coverage by icebergs (m^2)
  real, dimension(:,:), pointer :: mass=>null() ! Mass distribution (kg/m^2)
  real, dimension(:,:,:), pointer :: mass_on_ocean=>null() ! Mass distribution partitioned by neighbor (kg)  
  real, dimension(:,:,:), pointer :: area_on_ocean=>null() ! Area distribution partitioned by neighbor (m^2)  
  real, dimension(:,:,:), pointer :: Uvel_on_ocean=>null() ! zonal velocity distribution partitioned by neighbor (m^2* m/s)  
  real, dimension(:,:,:), pointer :: Vvel_on_ocean=>null() ! meridional momentum distribution partitioned by neighbor (m^2 m/s)  
  real, dimension(:,:), pointer :: tmp=>null() ! Temporary work space
  real, dimension(:,:), pointer :: tmpc=>null() ! Temporary work space
  real, dimension(:,:,:), pointer :: stored_ice=>null() ! Accumulated ice mass flux at calving locations (kg)
  real, dimension(:,:), pointer :: rmean_calving=>null() ! Running mean for ice calving
  real, dimension(:,:), pointer :: rmean_calving_hflx=>null() ! Running mean for ice calving
  real, dimension(:,:), pointer :: stored_heat=>null() ! Heat content of stored ice (J)
  real, dimension(:,:,:), pointer :: real_calving=>null() ! Calving rate into iceberg class at calving locations (kg/s)
  real, dimension(:,:), pointer :: iceberg_heat_content=>null() ! Distributed heat content of bergs (J/m^2)
  real, dimension(:,:), pointer :: parity_x=>null() ! X component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  real, dimension(:,:), pointer :: parity_y=>null() ! Y component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  integer, dimension(:,:), pointer :: iceberg_counter_grd=>null() ! Counts icebergs created for naming purposes
  logical :: rmean_calving_initialized = .false. ! True if rmean_calving(:,:) has been filled with meaningful data
  logical :: rmean_calving_hflx_initialized = .false. ! True if rmean_calving_hflx(:,:) has been filled with meaningful data
  ! Diagnostics handles
  integer :: id_uo=-1, id_vo=-1, id_calving=-1, id_stored_ice=-1, id_accum=-1, id_unused=-1, id_floating_melt=-1
  integer :: id_melt_buoy=-1, id_melt_eros=-1, id_melt_conv=-1, id_virtual_area=-1, id_real_calving=-1
  integer :: id_calving_hflx_in=-1, id_stored_heat=-1, id_melt_hflx=-1, id_heat_content=-1
  integer :: id_mass=-1, id_ui=-1, id_vi=-1, id_ua=-1, id_va=-1, id_sst=-1, id_cn=-1, id_hi=-1
  integer :: id_bergy_src=-1, id_bergy_melt=-1, id_bergy_mass=-1, id_berg_melt=-1
  integer :: id_rmean_calving=-1, id_rmean_calving_hflx=-1
  integer :: id_spread_mass=-1, id_spread_area=-1
  integer :: id_ssh=-1, id_fax=-1, id_fay=-1
  integer :: id_count=-1, id_chksum=-1, id_u_iceberg=-1, id_v_iceberg=-1, id_sss=-1, id_ustar_iceberg
  integer :: id_spread_uvel=-1, id_spread_vvel=-1
  integer :: id_melt_m_per_year=-1
  integer :: id_ocean_depth=-1

  real :: clipping_depth=0. ! The effective depth at which to clip the weight felt by the ocean [m].

end type icebergs_gridded

type :: xyt
  real :: lon, lat, day
  real :: mass, thickness, width, length, uvel, vvel
  real :: axn, ayn, bxn, byn, uvel_old, vvel_old, lat_old, lon_old  !Explicit and implicit accelerations !Alon 
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, sss, cn, hi, halo_berg, static_berg
  real :: mass_of_bits, heat_density
  integer :: year, iceberg_num
  type(xyt), pointer :: next=>null()
end type xyt

type :: iceberg
  type(iceberg), pointer :: prev=>null(), next=>null()
  ! State variables (specific to the iceberg, needed for restarts)
  real :: lon, lat, uvel, vvel, mass, thickness, width, length
  real :: axn, ayn, bxn, byn, uvel_old, vvel_old, lon_old, lat_old !Explicit and implicit accelerations !Alon 
  real :: start_lon, start_lat, start_day, start_mass, mass_scaling
  real :: mass_of_bits, heat_density
  real :: halo_berg  ! Equal to zero for bergs on computational domain, and =1 for bergs on the halo
  real :: static_berg  ! Equal to 1 for icebergs which are static (not allowed to move). Might be extended to grounding later.
  integer :: start_year
  integer :: iceberg_num
  integer :: ine, jne ! nearest index in NE direction (for convenience)
  real :: xi, yj ! Non-dimensional coords within current cell (0..1)
  ! Environment variables (as seen by the iceberg)
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, sss, cn, hi
  type(xyt), pointer :: trajectory=>null()
  type(bond), pointer :: first_bond=>null()  !First element of bond list.
end type iceberg

type :: bond
  type(bond), pointer :: prev_bond=>null(), next_bond=>null()
  type(iceberg), pointer :: other_berg=>null()
  integer :: other_berg_num, other_berg_ine, other_berg_jne 
end type bond

type :: buffer
  integer :: size=0
  real, dimension(:,:), pointer :: data
end type buffer

type :: linked_list
  type(iceberg), pointer :: first=>null()
end type linked_list

type :: icebergs !; private!Niki: Ask Alistair why this is private. ice_bergs_io cannot compile if this is private!
  type(icebergs_gridded), pointer :: grd
  type(linked_list), dimension(:,:), allocatable :: list
  type(xyt), pointer :: trajectories=>null()
  real :: dt           ! Time-step between iceberg calls (should make adaptive?)
  integer :: current_year
  real :: current_yearday ! 1.00-365.99
  integer :: traj_sample_hrs, traj_write_hrs
  integer :: verbose_hrs
  integer :: max_bonds
  integer :: clock, clock_mom, clock_the, clock_int, clock_cal, clock_com, clock_ini, clock_ior, clock_iow, clock_dia ! ids for fms timers
  integer :: clock_trw, clock_trp
  real :: rho_bergs ! Density of icebergs [kg/m^3]
  real :: spring_coef  ! Spring contant for iceberg interactions 
  real :: bond_coef  ! Spring contant for iceberg bonds 
  real :: radial_damping_coef     ! Coef for relative iceberg motion damping (radial component) -Alon
  real :: tangental_damping_coef     ! Coef for relative iceberg motion damping (tangental component) -Alon
  real :: LoW_ratio ! Initial ratio L/W for newly calved icebergs
  real :: bergy_bit_erosion_fraction ! Fraction of erosion melt flux to divert to bergy bits
  real :: sicn_shift ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
  real :: lat_ref=0. ! Reference latitude for f-plane (when this option is on)
  real :: u_override=0.0 ! Overrides the u velocity of icebergs (for ocean testing)
  real :: v_override=0.0 ! Overrides the v velocity of icebergs (for ocean testing)
  real :: utide_icebergs= 0.      ! Tidal speeds, set to zero for now.
  real :: ustar_icebergs_bg=0.001 ! Background u_star under icebergs. This should be linked to a value felt by the ocean boundary layer
  real :: cdrag_icebergs =  1.5e-3 !Momentum Drag coef, taken from HJ99  (Holland and Jenkins 1999)
  real :: initial_orientation=0. ! Iceberg orientaion relative to this angle (in degrees). Used for hexagonal mass spreading.
  real :: Gamma_T_3EQ=0.022 ! Nondimensional heat-transfer coefficient
  real :: melt_cutoff=-1.0 !Minimum ocean thickness for melting to occur (is not applied for values < 0)
  logical :: const_gamma=.True. !If true uses a constant heat tranfer coefficient, from which the salt transfer is calculated
  real, dimension(:), pointer :: initial_mass, distribution, mass_scaling
  real, dimension(:), pointer :: initial_thickness, initial_width, initial_length
  logical :: restarted=.false. ! Indicate whether we read state from a restart or not
  logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
  logical :: add_weight_to_ocean=.true. ! Add weight of bergs to ocean
  logical :: passive_mode=.false. ! Add weight of icebergs + bits to ocean
  logical :: time_average_weight=.false. ! Time average the weight on the ocean
  logical :: Runge_not_Verlet=.True.  !True=Runge Kuttai, False=Verlet.  - Added by Alon 
  logical :: apply_thickness_cutoff_to_gridded_melt=.False.  !Prevents melt for ocean thickness below melt_cuttoff (applied to gridded melt fields)
  logical :: apply_thickness_cutoff_to_bergs_melt=.False.  !Prevents melt for ocean thickness below melt_cuttoff (applied to bergs)
  logical :: use_updated_rolling_scheme=.false. ! True to use the aspect ratio based rolling scheme rather than incorrect version of WM scheme   (set tip_parameter=1000. for correct WM scheme)
  logical :: pass_fields_to_ocean_model=.False. !Iceberg area, mass and ustar fields are prepared to pass to ocean model
  logical :: use_mixed_layer_salinity_for_thermo=.False.  !If true, then model uses ocean salinity for 3 and 2 equation melt model.
  logical :: find_melt_using_spread_mass=.False.  !If true, then the model calculates ice loss by looping at the spread_mass before and after.
  logical :: Use_three_equation_model=.True.  !Uses 3 equation model for melt when ice shelf type thermodynamics are used.
  logical :: melt_icebergs_as_ice_shelf=.False.  !Uses iceshelf type thermodynamics
  logical :: Iceberg_melt_without_decay=.False.  !Allows icebergs meltwater fluxes to enter the ocean, without the iceberg decaying or changing shape.
  logical :: add_iceberg_thickness_to_SSH=.False.  !Adds the iceberg contribution to SSH.   
  logical :: override_iceberg_velocities=.False.  !Allows you to set a fixed iceberg velocity for all non-static icebergs.
  logical :: use_f_plane=.False.  !Flag to use a f-plane for the rotation
  logical :: rotate_icebergs_for_mass_spreading=.True.  !Flag allows icebergs to rotate for spreading their mass (in hexagonal spreading mode)
  logical :: set_melt_rates_to_zero=.False.  !Sets all melt rates to zero, for testing purposes (thermodynamics routine is still run)
  logical :: hexagonal_icebergs=.False. !True treats icebergs as rectangles, False as hexagonal elements (for the purpose of mass spreading)
  logical :: allow_bergs_to_roll=.True. !Allows icebergs to roll over when rolling conditions are met
  logical :: ignore_missing_restart_bergs=.False.  !True Allows the model to ignorm icebergs missing in the restart. 
  logical :: Static_icebergs=.False.  !True= icebergs do no move
  logical :: only_interactive_forces=.False.  !Icebergs only feel interactive forces, and not ocean, wind... 
  logical :: halo_debugging=.False.  !Use for debugging halos (remove when its working) 
  logical :: save_short_traj=.True.  !True saves only lon,lat,time,iceberg_num in iceberg_trajectory.nc 
  logical :: ignore_traj=.False.  !If true, then model does not traj trajectory data at all 
  logical :: iceberg_bonds_on=.False.  !True=Allow icebergs to have bonds, False=don't allow. 
  logical :: manually_initialize_bonds=.False.  !True= Bonds are initialize manually. 
  logical :: use_new_predictive_corrective =.False.  !Flag to use Bob's predictive corrective iceberg scheme- Added by Alon 
  logical :: interactive_icebergs_on=.false.  !Turn on/off interactions between icebergs  - Added by Alon 
  logical :: critical_interaction_damping_on=.true.  !Sets the damping on relative iceberg velocity to critical value - Added by Alon 
  logical :: use_old_spreading=.true. ! If true, spreads iceberg mass as if the berg is one grid cell wide
  logical :: read_ocean_depth_from_file=.false. ! If true, ocean depth is read from a file.
  integer :: debug_iceberg_with_id = -1 ! If positive, monitors a berg with this id

  real :: speed_limit=0. ! CFL speed limit for a berg [m/s]
  real :: tau_calving=0. ! Time scale for smoothing out calving field (years)
  real :: tip_parameter=0. ! parameter to override iceberg rollilng critica ratio (use zero to get parameter directly from ice and seawater densities) 
  real :: grounding_fraction=0. ! Fraction of water column depth at which grounding occurs
  type(buffer), pointer :: obuffer_n=>null(), ibuffer_n=>null()
  type(buffer), pointer :: obuffer_s=>null(), ibuffer_s=>null()
  type(buffer), pointer :: obuffer_e=>null(), ibuffer_e=>null()
  type(buffer), pointer :: obuffer_w=>null(), ibuffer_w=>null()
  type(buffer), pointer :: obuffer_io=>null(), ibuffer_io=>null()
  ! Budgets
  real :: net_calving_received=0., net_calving_returned=0.
  real :: net_incoming_calving=0., net_outgoing_calving=0.
  real :: net_incoming_calving_heat=0., net_outgoing_calving_heat=0.
  real :: net_incoming_calving_heat_used=0., net_heat_to_bergs=0.
  real :: stored_start=0., stored_end=0.
  real :: rmean_calving_start=0., rmean_calving_end=0.
  real :: rmean_calving_hflx_start=0., rmean_calving_hflx_end=0.
  real :: stored_heat_start=0., stored_heat_end=0., net_heat_to_ocean=0.
  real :: net_calving_used=0., net_calving_to_bergs=0.
  real :: floating_mass_start=0., floating_mass_end=0.
  real :: floating_heat_start=0., floating_heat_end=0.
  real :: icebergs_mass_start=0., icebergs_mass_end=0.
  real :: bergy_mass_start=0., bergy_mass_end=0.
  real :: spread_mass_start=0., spread_mass_end=0.
  real :: spread_area_start=0., spread_area_end=0.
  real :: u_iceberg_start=0., u_iceberg_end=0.
  real :: v_iceberg_start=0., v_iceberg_end=0.
  real :: spread_uvel_start=0., spread_uvel_end=0.
  real :: spread_vvel_start=0., spread_vvel_end=0.
  real :: ustar_iceberg_start=0., ustar_iceberg_end=0.
  real :: returned_mass_on_ocean=0.
  real :: returned_area_on_ocean=0.
  real :: net_melt=0., berg_melt=0., bergy_src=0., bergy_melt=0.
  integer :: nbergs_calved=0, nbergs_melted=0, nbergs_start=0, nbergs_end=0
  integer :: nspeeding_tickets=0
  integer :: nbonds=0
  integer, dimension(:), pointer :: nbergs_calved_by_class=>null()
end type icebergs

! Needs to be module global so can be public to icebergs_mod.
! Remove when backward compatibility no longer needed
logical :: orig_read=.false.

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains


! ##############################################################################

subroutine ice_bergs_framework_init(bergs, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth, maskmap, fractional_area)

use mpp_parameter_mod, only: SCALAR_PAIR, CGRID_NE, BGRID_NE, CORNER, AGRID
use mpp_domains_mod, only: mpp_update_domains, mpp_define_domains
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
use mpp_domains_mod, only: CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use mpp_domains_mod, only: mpp_define_io_domain

use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id, input_nml_file
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP

use fms_mod, only: open_namelist_file, check_nml_error, close_file
use fms_mod, only: clock_flag_default

use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use diag_manager_mod, only: diag_axis_init

! Arguments
type(icebergs), pointer :: bergs
integer, intent(in) :: gni, gnj, layout(2), io_layout(2), axes(2)
integer, intent(in) :: dom_x_flags, dom_y_flags
real, intent(in) :: dt
type (time_type), intent(in) :: Time ! current time
real, dimension(:,:), intent(in) :: ice_lon, ice_lat, ice_wet
real, dimension(:,:), intent(in) :: ice_dx, ice_dy, ice_area
real, dimension(:,:), intent(in) :: cos_rot, sin_rot
real, dimension(:,:), intent(in),optional :: ocean_depth
logical, intent(in), optional :: maskmap(:,:)
logical, intent(in), optional :: fractional_area

! Namelist parameters (and defaults)
integer :: halo=4 ! Width of halo region
integer :: traj_sample_hrs=24 ! Period between sampling of position for trajectory storage
integer :: traj_write_hrs=480 ! Period between writing sampled trajectories to disk
integer :: verbose_hrs=24 ! Period between verbose messages
integer :: max_bonds=6 ! Maximum number of iceberg bond passed between processors
real :: rho_bergs=850. ! Density of icebergs
real :: spring_coef=1.e-8  ! Spring contant for iceberg interactions (this seems to be the highest stable value)
real :: bond_coef=1.e-8 ! Spring contant for iceberg bonds - not being used right now
real :: radial_damping_coef=1.e-4     ! Coef for relative iceberg motion damping (radial component) -Alon
real :: tangental_damping_coef=2.e-5     ! Coef for relative iceberg motion damping (tangental component) -Alon
real :: LoW_ratio=1.5 ! Initial ratio L/W for newly calved icebergs
real :: bergy_bit_erosion_fraction=0. ! Fraction of erosion melt flux to divert to bergy bits
real :: sicn_shift=0. ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
real :: lat_ref=0. ! Reference latitude for f-plane (when this option is on)
real :: u_override=0.0 ! Overrides the u velocity of icebergs (for ocean testing)
real :: v_override=0.0 ! Overrides the v velocity of icebergs (for ocean testing)
real :: Lx=360. ! Length of domain in x direction, used for periodicity (use a huge number for non-periodic)
real :: initial_orientation=0. ! Iceberg orientaion relative to this angle (in degrees). Used for hexagonal mass spreading.
real :: utide_icebergs= 0.      ! Tidal speeds, set to zero for now.
real :: ustar_icebergs_bg=0.001 ! Background u_star under icebergs. This should be linked to a value felt by the ocean boundary layer
real :: cdrag_icebergs =  1.5e-3 !Momentum Drag coef, taken from HJ99  (Holland and Jenkins 1999)
real :: Gamma_T_3EQ=0.022 ! Nondimensional heat-transfer coefficient
real :: melt_cutoff=-1.0 !Minimum ocean thickness for melting to occur (is not applied for values < 0)
logical :: const_gamma=.True. !If true uses a constant heat tranfer coefficient, from which the salt transfer is calculated
logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
logical :: add_weight_to_ocean=.true. ! Add weight of icebergs + bits to ocean
logical :: passive_mode=.false. ! Add weight of icebergs + bits to ocean
logical :: time_average_weight=.false. ! Time average the weight on the ocean
real :: speed_limit=0. ! CFL speed limit for a berg
real :: tau_calving=0. ! Time scale for smoothing out calving field (years)
real :: tip_parameter=0. ! parameter to override iceberg rollilng critica ratio (use zero to get parameter directly from ice and seawater densities
real :: grounding_fraction=0. ! Fraction of water column depth at which grounding occurs
logical :: Runge_not_Verlet=.True.  !True=Runge Kutta, False=Verlet.  - Added by Alon 
logical :: apply_thickness_cutoff_to_gridded_melt=.False.  !Prevents melt for ocean thickness below melt_cuttoff (applied to gridded melt fields)
logical :: apply_thickness_cutoff_to_bergs_melt=.False.  !Prevents melt for ocean thickness below melt_cuttoff (applied to bergs)
logical :: use_updated_rolling_scheme=.false. ! Use the corrected Rolling Scheme rather than the erronios one
logical :: pass_fields_to_ocean_model=.False. !Iceberg area, mass and ustar fields are prepared to pass to ocean model
logical :: use_mixed_layer_salinity_for_thermo=.False.  !If true, then model uses ocean salinity for 3 and 2 equation melt model.
logical :: find_melt_using_spread_mass=.False.  !If true, then the model calculates ice loss by looping at the spread_mass before and after.
logical :: Use_three_equation_model=.True.  !Uses 3 equation model for melt when ice shelf type thermodynamics are used.
logical :: melt_icebergs_as_ice_shelf=.False.  !Uses iceshelf type thermodynamics
logical :: Iceberg_melt_without_decay=.False.  !Allows icebergs meltwater fluxes to enter the ocean, without the iceberg decaying or changing shape.
logical :: add_iceberg_thickness_to_SSH=.False.  !Adds the iceberg contribution to SSH.   
logical :: override_iceberg_velocities=.False.  !Allows you to set a fixed iceberg velocity for all non-static icebergs.
logical :: use_f_plane=.False.  !Flag to use a f-plane for the rotation
logical :: grid_is_latlon=.True.  !True means that the grid is specified in lat lon, and uses to radius of the earth to convert to distance
logical :: grid_is_regular=.True. !Flag to say whether point in cell can be found assuming regular cartesian grid
logical :: rotate_icebergs_for_mass_spreading=.True.  !Flag allows icebergs to rotate for spreading their mass (in hexagonal spreading mode)
logical :: set_melt_rates_to_zero=.False.  !Sets all melt rates to zero, for testing purposes (thermodynamics routine is still run)
logical :: allow_bergs_to_roll=.True. !Allows icebergs to roll over when rolling conditions are met
logical :: hexagonal_icebergs=.False. !True treats icebergs as rectangles, False as hexagonal elements (for the purpose of mass spreading)
logical :: ignore_missing_restart_bergs=.False.  !True Allows the model to ignorm icebergs missing in the restart. 
logical :: Static_icebergs=.False.  !True= icebergs do no move
logical :: only_interactive_forces=.False.  !Icebergs only feel interactive forces, and not ocean, wind... 
logical :: halo_debugging=.False.  !Use for debugging halos (remove when its working) 
logical :: save_short_traj=.True.  !True saves only lon,lat,time,iceberg_num in iceberg_trajectory.nc 
logical :: ignore_traj=.False.  !If true, then model does not traj trajectory data at all 
logical :: iceberg_bonds_on=.False.  !True=Allow icebergs to have bonds, False=don't allow. 
logical :: manually_initialize_bonds=.False.  !True= Bonds are initialize manually. 
logical :: use_new_predictive_corrective =.False.  !Flag to use Bob's predictive corrective iceberg scheme- Added by Alon 
logical :: interactive_icebergs_on=.false.  !Turn on/off interactions between icebergs  - Added by Alon 
logical :: critical_interaction_damping_on=.true.  !Sets the damping on relative iceberg velocity to critical value - Added by Alon 
logical :: do_unit_tests=.false. ! Conduct some unit tests
logical :: input_freq_distribution=.false. ! Flag to show if input distribution is freq or mass dist (=1 if input is a freq dist, =0 to use an input mass dist)
logical :: read_old_restarts=.false. ! Legacy option that does nothing
logical :: use_old_spreading=.true. ! If true, spreads iceberg mass as if the berg is one grid cell wide
logical :: read_ocean_depth_from_file=.false. ! If true, ocean depth is read from a file.
real, dimension(nclasses) :: initial_mass=(/8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11/) ! Mass thresholds between iceberg classes (kg)
real, dimension(nclasses) :: distribution=(/0.24, 0.12, 0.15, 0.18, 0.12, 0.07, 0.03, 0.03, 0.03, 0.02/) ! Fraction of calving to apply to this class (non-dim) , 
real, dimension(nclasses) :: mass_scaling=(/2000, 200, 50, 20, 10, 5, 2, 1, 1, 1/) ! Ratio between effective and real iceberg mass (non-dim)
real, dimension(nclasses) :: initial_thickness=(/40., 67., 133., 175., 250., 250., 250., 250., 250., 250./) ! Total thickness of newly calved bergs (m)
integer :: debug_iceberg_with_id = -1 ! If positive, monitors a berg with this id


namelist /icebergs_nml/ verbose, budget, halo,  traj_sample_hrs, initial_mass, traj_write_hrs, max_bonds, save_short_traj,Static_icebergs,  &
         distribution, mass_scaling, initial_thickness, verbose_hrs, spring_coef,bond_coef, radial_damping_coef, tangental_damping_coef, only_interactive_forces, &
         rho_bergs, LoW_ratio, debug, really_debug, use_operator_splitting, bergy_bit_erosion_fraction, iceberg_bonds_on, manually_initialize_bonds, ignore_missing_restart_bergs, &
         parallel_reprod, use_slow_find, sicn_shift, add_weight_to_ocean, passive_mode, ignore_ij_restart, use_new_predictive_corrective, halo_debugging, hexagonal_icebergs, &
         time_average_weight, generate_test_icebergs, speed_limit, fix_restart_dates, use_roundoff_fix, Runge_not_Verlet, interactive_icebergs_on, critical_interaction_damping_on, &
         old_bug_rotated_weights, make_calving_reproduce,restart_input_dir, orig_read, old_bug_bilin,do_unit_tests,grounding_fraction, input_freq_distribution, force_all_pes_traj, &
         allow_bergs_to_roll,set_melt_rates_to_zero,lat_ref,initial_orientation,rotate_icebergs_for_mass_spreading,grid_is_latlon,Lx,use_f_plane,use_old_spreading, &
         grid_is_regular,override_iceberg_velocities,u_override,v_override,add_iceberg_thickness_to_SSH,Iceberg_melt_without_decay,melt_icebergs_as_ice_shelf, &
         Use_three_equation_model,find_melt_using_spread_mass,use_mixed_layer_salinity_for_thermo,utide_icebergs,ustar_icebergs_bg,cdrag_icebergs, pass_fields_to_ocean_model, &
         const_gamma, Gamma_T_3EQ, ignore_traj, debug_iceberg_with_id,use_updated_rolling_scheme, tip_parameter, read_old_restarts, tau_calving, read_ocean_depth_from_file, melt_cutoff,&
         apply_thickness_cutoff_to_gridded_melt, apply_thickness_cutoff_to_bergs_melt 

! Local variables
integer :: ierr, iunit, i, j, id_class, axes3d(3), is,ie,js,je,np
type(icebergs_gridded), pointer :: grd
real :: lon_mod, big_number
logical :: lerr
integer :: stdlogunit, stderrunit
real :: Total_mass  !Added by Alon 

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "ice_bergs_framework: "//trim(version)

! Read namelist parameters
 !write(stderrunit,*) 'diamonds: reading namelist'
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=icebergs_nml, iostat=ierr)
#else
  iunit = open_namelist_file()
  read  (iunit, icebergs_nml,iostat=ierr)
  call close_file (iunit)
#endif
  ierr = check_nml_error(ierr,'icebergs_nml')

  if (really_debug) debug=.true. ! One implies the other...

  write (stdlogunit, icebergs_nml)


! Allocate overall structure
 !write(stderrunit,*) 'diamonds: allocating bergs'
  allocate(bergs)
  allocate(bergs%grd)
  grd=>bergs%grd ! For convenience to avoid bergs%grd%X
 !write(stderrunit,*) 'diamonds: allocating domain'
  allocate(grd%domain)

! Clocks
  bergs%clock=mpp_clock_id( 'Icebergs', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  bergs%clock_mom=mpp_clock_id( 'Icebergs-momentum', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_the=mpp_clock_id( 'Icebergs-thermodyn', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_int=mpp_clock_id( 'Icebergs-interface', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_cal=mpp_clock_id( 'Icebergs-calving', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_com=mpp_clock_id( 'Icebergs-communication', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_ini=mpp_clock_id( 'Icebergs-initialization', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_ior=mpp_clock_id( 'Icebergs-I/O read', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_iow=mpp_clock_id( 'Icebergs-I/O write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_dia=mpp_clock_id( 'Icebergs-diagnostics', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_ini)

! Set up iceberg domain
 !write(stderrunit,*) 'diamonds: defining domain'
  call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
                           maskmap=maskmap, &
                           xflags=dom_x_flags, xhalo=halo,  &
                           yflags=dom_y_flags, yhalo=halo, name='diamond')

  call mpp_define_io_domain(grd%domain, io_layout)

 !write(stderrunit,*) 'diamond: get compute domain'
  call mpp_get_compute_domain( grd%domain, grd%isc, grd%iec, grd%jsc, grd%jec )
  call mpp_get_data_domain( grd%domain, grd%isd, grd%ied, grd%jsd, grd%jed )
  call mpp_get_global_domain( grd%domain, grd%isg, grd%ieg, grd%jsg, grd%jeg )

  call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
  call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
  call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
  call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)


  folded_north_on_pe = ((dom_y_flags == FOLD_NORTH_EDGE) .and. (grd%jec == gnj)) 
 !write(stderrunit,'(a,6i4)') 'diamonds, icebergs_init: pe,n,s,e,w =',mpp_pe(),grd%pe_N,grd%pe_S,grd%pe_E,grd%pe_W, NULL_PE

 !if (verbose) &
 !write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'diamonds, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
 !     grd%isc,grd%iec,grd%jsc,grd%jec, &
 !     ' [lon|lat][min|max]=', minval(ice_lon),maxval(ice_lon),minval(ice_lat),maxval(ice_lat)
 !write(stderrunit,*) 'diamonds, int args = ', mpp_pe(),gni, gnj, layout, axes

 ! Allocate grid of pointers
  allocate( bergs%list(grd%isd:grd%ied, grd%jsd:grd%jed) )
  do j = grd%jsd,grd%jed ; do i = grd%isd,grd%ied
    bergs%list(i,j)%first => null()
  enddo ; enddo

  big_number=1.0E15
 !write(stderrunit,*) 'diamonds: allocating grid'
  allocate( grd%lon(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=big_number
  allocate( grd%lat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=big_number
  allocate( grd%lonc(grd%isd:grd%ied, grd%jsd:grd%jed) );grd%lon(:,:)=big_number
  allocate( grd%latc(grd%isd:grd%ied, grd%jsd:grd%jed) );grd%lat(:,:)=big_number
  allocate( grd%dx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dx(:,:)=0.
  allocate( grd%dy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dy(:,:)=0.
  allocate( grd%area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%area(:,:)=0.
  allocate( grd%msk(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%msk(:,:)=0.
  allocate( grd%cos(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cos(:,:)=1.
  allocate( grd%sin(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sin(:,:)=0.
  allocate( grd%ocean_depth(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ocean_depth(:,:)=0.
  allocate( grd%calving(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%calving(:,:)=0.
  allocate( grd%calving_hflx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%calving_hflx(:,:)=0.
  allocate( grd%stored_heat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%stored_heat(:,:)=0.
  allocate( grd%floating_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%floating_melt(:,:)=0.
  allocate( grd%berg_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%berg_melt(:,:)=0.
  allocate( grd%melt_buoy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_buoy(:,:)=0.
  allocate( grd%melt_eros(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_eros(:,:)=0.
  allocate( grd%melt_conv(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_conv(:,:)=0.
  allocate( grd%bergy_src(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%bergy_src(:,:)=0.
  allocate( grd%bergy_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%bergy_melt(:,:)=0.
  allocate( grd%bergy_mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%bergy_mass(:,:)=0.
  allocate( grd%spread_mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%spread_mass(:,:)=0.
  allocate( grd%spread_mass_old(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%spread_mass_old(:,:)=0.
  allocate( grd%spread_area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%spread_area(:,:)=0.
  allocate( grd%u_iceberg(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%u_iceberg(:,:)=0.
  allocate( grd%v_iceberg(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%v_iceberg(:,:)=0.
  allocate( grd%spread_uvel(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%spread_uvel(:,:)=0.
  allocate( grd%spread_vvel(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%spread_vvel(:,:)=0.
  allocate( grd%ustar_iceberg(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ustar_iceberg(:,:)=0.
  allocate( grd%virtual_area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%virtual_area(:,:)=0.
  allocate( grd%mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%mass(:,:)=0.
  allocate( grd%mass_on_ocean(grd%isd:grd%ied, grd%jsd:grd%jed, 9) ); grd%mass_on_ocean(:,:,:)=0.
  allocate( grd%area_on_ocean(grd%isd:grd%ied, grd%jsd:grd%jed, 9) ); grd%area_on_ocean(:,:,:)=0.
  allocate( grd%Uvel_on_ocean(grd%isd:grd%ied, grd%jsd:grd%jed, 9) ); grd%Uvel_on_ocean(:,:,:)=0.
  allocate( grd%Vvel_on_ocean(grd%isd:grd%ied, grd%jsd:grd%jed, 9) ); grd%Vvel_on_ocean(:,:,:)=0.
  allocate( grd%stored_ice(grd%isd:grd%ied, grd%jsd:grd%jed, nclasses) ); grd%stored_ice(:,:,:)=0.
  allocate( grd%rmean_calving(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%rmean_calving(:,:)=0.
  allocate( grd%rmean_calving_hflx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%rmean_calving_hflx(:,:)=0.
  allocate( grd%real_calving(grd%isd:grd%ied, grd%jsd:grd%jed, nclasses) ); grd%real_calving(:,:,:)=0.
  allocate( grd%uo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%uo(:,:)=0.
  allocate( grd%vo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vo(:,:)=0.
  allocate( grd%ui(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ui(:,:)=0.
  allocate( grd%vi(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vi(:,:)=0.
  allocate( grd%ua(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ua(:,:)=0.
  allocate( grd%va(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%va(:,:)=0.
  allocate( grd%ssh(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ssh(:,:)=0.
  allocate( grd%sst(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sst(:,:)=0.
  allocate( grd%sss(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sss(:,:)=0.
  allocate( grd%cn(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cn(:,:)=0.
  allocate( grd%hi(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%hi(:,:)=0.
  allocate( grd%tmp(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%tmp(:,:)=0.
  allocate( grd%tmpc(grd%isc:grd%iec, grd%jsc:grd%jec) ); grd%tmpc(:,:)=0.
  allocate( bergs%nbergs_calved_by_class(nclasses) ); bergs%nbergs_calved_by_class(:)=0
  allocate( grd%parity_x(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_x(:,:)=1.
  allocate( grd%parity_y(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_y(:,:)=1.
  allocate( grd%iceberg_counter_grd(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%iceberg_counter_grd(:,:)=1

 !write(stderrunit,*) 'diamonds: copying grid'
  ! Copy data declared on ice model computational domain
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
  grd%lon(is:ie,js:je)=ice_lon(:,:)
  grd%lat(is:ie,js:je)=ice_lat(:,:)
  grd%area(is:ie,js:je)=ice_area(:,:) !sis2 has *(4.*pi*radius*radius)

  !!!!!!!!!!!!!!!debugging!!!!!!!!!!!!!!!!!!
  !if (mpp_pe().eq.5) then
  ! write(stderrunit,'(a3,32i7)') 'LB',(i,i=grd%isd,grd%ied)
  ! do j=grd%jed,grd%jsd,-1
  !   write(stderrunit,'(i3,32f7.1)') j,(grd%lon(i,j),i=grd%isd,grd%ied)
  ! enddo
  ! write(stderrunit,'(a3,32i7)') 'Ice lon',(i,i=grd%isd,grd%ied)
  ! do j=grd%jed,grd%jsd,-1
  !   write(stderrunit,'(i3,32f7.1)') j,(ice_lon(i,j),i=grd%isd,grd%ied)
  ! enddo
  ! write(stderrunit,'(a3,32i7)') 'LA',(i,i=grd%isd,grd%ied)
  ! do j=grd%jed,grd%jsd,-1
  !   write(stderrunit,'(i3,32f7.1)') j,(grd%lon(i,j),i=grd%isd,grd%ied)
  ! enddo
  !endif
  !!!!!!!!!!!!!!!debugging!!!!!!!!!!!!!!!!!!

  !For SIS not to change answers
  if(present(fractional_area)) then
    if(fractional_area) grd%area(is:ie,js:je)=ice_area(:,:) *(4.*pi*radius*radius)
  endif
  if(present(ocean_depth)) grd%ocean_depth(is:ie,js:je)=ocean_depth(:,:)
  
  ! Copy data declared on ice model data domain
  is=grd%isc-1; ie=grd%iec+1; js=grd%jsc-1; je=grd%jec+1
  grd%dx(is:ie,js:je)=ice_dx(:,:)
  grd%dy(is:ie,js:je)=ice_dy(:,:)
  grd%msk(is:ie,js:je)=ice_wet(:,:)
  grd%cos(is:ie,js:je)=cos_rot(:,:)
  grd%sin(is:ie,js:je)=sin_rot(:,:)

  call mpp_update_domains(grd%lon, grd%domain, position=CORNER)
  call mpp_update_domains(grd%lat, grd%domain, position=CORNER)
  call mpp_update_domains(grd%dy, grd%dx, grd%domain, gridtype=CGRID_NE, flags=SCALAR_PAIR)
  call mpp_update_domains(grd%area, grd%domain)
  call mpp_update_domains(grd%msk, grd%domain)
  call mpp_update_domains(grd%cos, grd%domain, position=CORNER)
  call mpp_update_domains(grd%sin, grd%domain, position=CORNER)
  call mpp_update_domains(grd%ocean_depth, grd%domain)
  call mpp_update_domains(grd%parity_x, grd%parity_y, grd%domain, gridtype=AGRID) ! If either parity_x/y is -ve, we need rotation of vectors

  ! Sanitize lon and lat in the southern halo
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=grd%lon(i,j+1)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
  enddo; enddo

  ! fix halos on edge of the domain
  !1) South
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i,j+1)-grd%lon(i,j+2)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
  enddo; enddo
  !2) North
  do j=grd%jec+1,grd%jed; do i=grd%isd,grd%ied
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i,j-1)-grd%lon(i,j-2)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i,j-1)-grd%lat(i,j-2)
  enddo; enddo
  !3) West
  do i=grd%isc-1,grd%isd,-1; do j=grd%jsd,grd%jed
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i+1,j)-grd%lon(i+2,j)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i+1,j)-grd%lat(i+2,j)
  enddo; enddo
  !4) East
  do i=grd%iec+1,grd%ied; do j=grd%jsd,grd%jed
      if (grd%lon(i,j).ge.big_number) grd%lon(i,j)=2.*grd%lon(i-1,j)-grd%lon(i-2,j)
      if (grd%lat(i,j).ge.big_number) grd%lat(i,j)=2.*grd%lat(i-1,j)-grd%lat(i-2,j)
  enddo; enddo

  if (.not. present(maskmap)) then ! Using a maskmap causes tickles this sanity check
    do j=grd%jsd,grd%jed; do i=grd%isd,grd%ied
      !if (grd%lon(i,j).ge.big_number) write(stderrunit,*) 'bad lon: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1,grd%lon(i,j)
      !if (grd%lat(i,j).ge.big_number) write(stderrunit,*) 'bad lat: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1,grd%lat(i,j)
    enddo; enddo
  endif

  if ((Lx.gt.1E15 ) .and. (mpp_pe().eq.mpp_root_pe())) then
          call error_mesg('diamonds, framework', 'Model does not enjoy the domain being larger than 1E15. Not sure why. Probably to do with floating point precision.', WARNING) 
  endif
  if ((.not. grid_is_latlon) .and. (Lx.eq.360.)) then
    if (mpp_pe().eq.mpp_root_pe())  then
            call error_mesg('diamonds, framework', 'Since the lat/lon grid is off, the x-direction is being set as non-periodic. Set Lx not equal to 360 override.', WARNING) 
    endif
    Lx=-1.
  endif


 !The fix to reproduce across PE layout change, from AJA
  if (Lx>0.) then
    j=grd%jsc; do i=grd%isc+1,grd%ied
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i-1,j),Lx)
      if (abs(grd%lon(i,j)-lon_mod)>(Lx/2.)) &
        grd%lon(i,j)= lon_mod
    enddo
    j=grd%jsc; do i=grd%isc-1,grd%isd,-1
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i+1,j) ,Lx)
      if (abs(grd%lon(i,j)-  lon_mod )>(Lx/2.)) &
        grd%lon(i,j)= lon_mod
    enddo
    do j=grd%jsc+1,grd%jed; do i=grd%isd,grd%ied
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i,j-1) ,Lx)
      if (abs(grd%lon(i,j)-(lon_mod ))>(Lx/2.)) &
        grd%lon(i,j)= lon_mod
    enddo; enddo
    do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      lon_mod = apply_modulo_around_point(grd%lon(i,j),grd%lon(i,j+1) ,Lx)
      if (abs(grd%lon(i,j)- lon_mod )>(Lx/2.)) &
        grd%lon(i,j)=  lon_mod
    enddo; enddo
  endif



  ! lonc, latc used for searches
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
    grd%lonc(i,j)=0.25*( (grd%lon(i,j)+grd%lon(i-1,j-1)) &
                        +(grd%lon(i-1,j)+grd%lon(i,j-1)) )
    grd%latc(i,j)=0.25*( (grd%lat(i,j)+grd%lat(i-1,j-1)) &
                        +(grd%lat(i-1,j)+grd%lat(i,j-1)) )
  enddo; enddo

  if (debug) then
    write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'diamonds, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
         grd%isc,grd%iec,grd%jsc,grd%jec, &
         ' [lon|lat][min|max]=', minval(grd%lon),maxval(grd%lon),minval(grd%lat),maxval(grd%lat)
  endif

 !if (mpp_pe().eq.5) then
 !  write(stderrunit,'(a3,32i7)') 'Lon',(i,i=grd%isd,grd%ied)
 !  do j=grd%jed,grd%jsd,-1
 !    write(stderrunit,'(i3,32f7.1)') j,(grd%lon(i,j),i=grd%isd,grd%ied)
 !  enddo
 !  write(stderrunit,'(a3,32i7)') 'Lat',(i,i=grd%isd,grd%ied)
 !  do j=grd%jed,grd%jsd,-1
 !    write(stderrunit,'(i3,32f7.1)') j,(grd%lat(i,j),i=grd%isd,grd%ied)
 !  enddo
 !  write(stderrunit,'(a3,32i7)') 'Msk',(i,i=grd%isd,grd%ied)
 !  do j=grd%jed,grd%jsd,-1
 !    write(stderrunit,'(i3,32f7.1)') j,(grd%msk(i,j),i=grd%isd,grd%ied)
 !  enddo
 !endif

! Final check for NaN's in the latlon grid:
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
    if (grd%lat(i,j) .ne. grd%lat(i,j)) then
      write(stderrunit,*) 'Lat not defined properly', mpp_pe(),i,j,grd%lat(i,j)
      call error_mesg('diamonds,grid defining', 'Latitude contains NaNs', FATAL)
    endif
    if (grd%lon(i,j) .ne. grd%lon(i,j)) then
      write(stderrunit,*) 'Lon not defined properly', mpp_pe(),i,j,grd%lon(i,j)
      call error_mesg('diamonds, grid defining', 'Longatudes contains NaNs', FATAL)
    endif
  enddo; enddo


!Added by Alon  - If a freq distribution is input, we have to convert the freq distribution to a mass flux distribution)
if (input_freq_distribution) then
     Total_mass=0.
     do j=1,nclasses
          Total_mass=Total_mass+(distribution(j)*initial_mass(j))
     enddo
     do j=1,nclasses
           distribution(j)=(distribution(j)*initial_mass(j))/Total_mass
     enddo
endif 

if ((halo .lt. 3) .and. (rotate_icebergs_for_mass_spreading .and. iceberg_bonds_on) )   then
    halo=3
    call error_mesg('diamonds, framework', 'Setting iceberg halos =3, since halos must be >= 3 for rotating icebergs for mass spreading', WARNING) 
elseif  ((halo .lt. 2) .and. (interactive_icebergs_on .or. iceberg_bonds_on) )   then
    halo=2
    call error_mesg('diamonds, framework', 'Setting iceberg halos =2, since halos must be >= 2 for interactions', WARNING) 
endif

if (interactive_icebergs_on) then
  if (Runge_not_Verlet) then
    !Runge_not_Verlet=.false.  ! Iceberg interactions only with Verlet
    call error_mesg('diamonds, framework', 'It is unlcear whther interactive icebergs work with Runge Kutta stepping.', WARNING) 
  endif
endif
if (.not.interactive_icebergs_on) then
  if (iceberg_bonds_on) then  
    !iceberg_bonds_on=.false.  
    call error_mesg('diamonds, framework', 'Interactive icebergs off requires iceberg bonds off (turning bonds off).', WARNING) 
  endif
endif
if (.not. iceberg_bonds_on) then
   max_bonds=0
else
  buffer_width=buffer_width+(max_bonds*3) ! Increase buffer width to include bonds being passed between processors 
endif
if (save_short_traj) buffer_width_traj=5 ! This is the length of the short buffer used for abrevated traj
if (ignore_traj) buffer_width_traj=0 ! If this is true, then all traj files should be ignored


 ! Parameters
  bergs%dt=dt
  bergs%traj_sample_hrs=traj_sample_hrs
  bergs%traj_write_hrs=traj_write_hrs
  bergs%save_short_traj=save_short_traj
  bergs%ignore_traj=ignore_traj
  bergs%verbose_hrs=verbose_hrs
  bergs%grd%halo=halo
  bergs%grd%Lx=Lx
  bergs%grd%grid_is_latlon=grid_is_latlon  
  bergs%grd%grid_is_regular=grid_is_regular 
  bergs%max_bonds=max_bonds
  bergs%rho_bergs=rho_bergs
  bergs%spring_coef=spring_coef
  bergs%bond_coef=bond_coef
  bergs%radial_damping_coef=radial_damping_coef
  bergs%tangental_damping_coef=tangental_damping_coef
  bergs%LoW_ratio=LoW_ratio
  bergs%use_operator_splitting=use_operator_splitting
  bergs%bergy_bit_erosion_fraction=bergy_bit_erosion_fraction
  bergs%sicn_shift=sicn_shift
  bergs%passive_mode=passive_mode
  bergs%time_average_weight=time_average_weight
  bergs%speed_limit=speed_limit
  bergs%tau_calving=tau_calving
  bergs%tip_parameter=tip_parameter
  bergs%use_updated_rolling_scheme=use_updated_rolling_scheme  !Alon
  bergs%Runge_not_Verlet=Runge_not_Verlet   
  bergs%apply_thickness_cutoff_to_bergs_melt=apply_thickness_cutoff_to_bergs_melt
  bergs%apply_thickness_cutoff_to_gridded_melt=apply_thickness_cutoff_to_gridded_melt
  bergs%melt_cutoff=melt_cutoff 
  bergs%read_ocean_depth_from_file=read_ocean_depth_from_file
  bergs%const_gamma=const_gamma 
  bergs%Gamma_T_3EQ=Gamma_T_3EQ
  bergs%pass_fields_to_ocean_model=pass_fields_to_ocean_model 
  bergs%ustar_icebergs_bg=ustar_icebergs_bg   
  bergs%utide_icebergs=utide_icebergs  
  bergs%cdrag_icebergs=cdrag_icebergs  
  bergs%use_mixed_layer_salinity_for_thermo=use_mixed_layer_salinity_for_thermo 
  bergs%find_melt_using_spread_mass=find_melt_using_spread_mass 
  bergs%Use_three_equation_model=Use_three_equation_model 
  bergs%melt_icebergs_as_ice_shelf=melt_icebergs_as_ice_shelf 
  bergs%Iceberg_melt_without_decay=Iceberg_melt_without_decay 
  bergs%add_iceberg_thickness_to_SSH=add_iceberg_thickness_to_SSH  
  bergs%override_iceberg_velocities=override_iceberg_velocities 
  bergs%use_f_plane=use_f_plane 
  bergs%rotate_icebergs_for_mass_spreading=rotate_icebergs_for_mass_spreading 
  bergs%lat_ref=lat_ref
  bergs%u_override=u_override
  bergs%v_override=v_override
  bergs%initial_orientation=initial_orientation
  bergs%set_melt_rates_to_zero=set_melt_rates_to_zero 
  bergs%allow_bergs_to_roll=allow_bergs_to_roll 
  bergs%hexagonal_icebergs=hexagonal_icebergs 
  bergs%ignore_missing_restart_bergs=ignore_missing_restart_bergs
  bergs%Static_icebergs=Static_icebergs 
  bergs%only_interactive_forces=only_interactive_forces
  bergs%halo_debugging=halo_debugging
  bergs%iceberg_bonds_on=iceberg_bonds_on   !Alon
  bergs%manually_initialize_bonds=manually_initialize_bonds   !Alon
  bergs%critical_interaction_damping_on=critical_interaction_damping_on   !Alon
  bergs%interactive_icebergs_on=interactive_icebergs_on   !Alon
  bergs%use_new_predictive_corrective=use_new_predictive_corrective  !Alon
  bergs%grounding_fraction=grounding_fraction
  bergs%add_weight_to_ocean=add_weight_to_ocean
  bergs%use_old_spreading=use_old_spreading
  bergs%debug_iceberg_with_id=debug_iceberg_with_id
  allocate( bergs%initial_mass(nclasses) ); bergs%initial_mass(:)=initial_mass(:)
  allocate( bergs%distribution(nclasses) ); bergs%distribution(:)=distribution(:)
  allocate( bergs%mass_scaling(nclasses) ); bergs%mass_scaling(:)=mass_scaling(:)
  allocate( bergs%initial_thickness(nclasses) ); bergs%initial_thickness(:)=initial_thickness(:)
  allocate( bergs%initial_width(nclasses) )
  allocate( bergs%initial_length(nclasses) )
  bergs%initial_width(:)=sqrt(initial_mass(:)/(LoW_ratio*rho_bergs*initial_thickness(:)))
  bergs%initial_length(:)=LoW_ratio*bergs%initial_width(:)

  if (read_old_restarts) call error_mesg('diamonds, ice_bergs_framework_init', 'Setting "read_old_restarts=.true." is obsolete and does nothing!', WARNING)

  ! Diagnostics
  id_class = diag_axis_init('mass_class', initial_mass, 'kg','Z', 'iceberg mass')
  axes3d(1:2)=axes
  axes3d(3)=id_class
  grd%id_calving=register_diag_field('icebergs', 'calving', axes, Time, &
     'Incoming Calving mass rate', 'kg/s')
  grd%id_calving_hflx_in=register_diag_field('icebergs', 'calving_hflx_in', axes, Time, &
     'Incoming Calving heat flux', 'J/s')
  grd%id_accum=register_diag_field('icebergs', 'accum_calving', axes, Time, &
     'Accumulated calving mass rate', 'kg/s')
  grd%id_unused=register_diag_field('icebergs', 'unused_calving', axes, Time, &
     'Unused calving mass rate', 'kg/s')
  grd%id_floating_melt=register_diag_field('icebergs', 'melt', axes, Time, &
     'Melt rate of icebergs + bits', 'kg/(m^2*s)')
  grd%id_melt_m_per_year=register_diag_field('icebergs', 'melt_m_per_year', axes, Time, &
     'Melt rate of icebergs + bits (m/yr)', 'm/yr')
  grd%id_berg_melt=register_diag_field('icebergs', 'berg_melt', axes, Time, &
     'Melt rate of icebergs', 'kg/(m^2*s)')
  grd%id_melt_buoy=register_diag_field('icebergs', 'melt_buoy', axes, Time, &
     'Buoyancy component of iceberg melt rate', 'kg/(m^2*s)')
  grd%id_melt_eros=register_diag_field('icebergs', 'melt_eros', axes, Time, &
     'Erosion component of iceberg melt rate', 'kg/(m^2*s)')
  grd%id_melt_conv=register_diag_field('icebergs', 'melt_conv', axes, Time, &
     'Convective component of iceberg melt rate', 'kg/(m^2*s)')
  grd%id_bergy_src=register_diag_field('icebergs', 'bergy_src', axes, Time, &
     'Mass source of bergy bits', 'kg/(m^2*s)')
  grd%id_bergy_melt=register_diag_field('icebergs', 'bergy_melt', axes, Time, &
     'Melt rate of bergy bits', 'kg/(m^2*s)')
  grd%id_bergy_mass=register_diag_field('icebergs', 'bergy_mass', axes, Time, &
     'Bergy bit density field', 'kg/(m^2)')
  grd%id_spread_mass=register_diag_field('icebergs', 'spread_mass', axes, Time, &
     'Iceberg mass after spreading', 'kg/(m^2)')
  grd%id_spread_area=register_diag_field('icebergs', 'spread_area', axes, Time, &
     'Iceberg area after spreading', 'm^2/(m^2)')
  grd%id_u_iceberg=register_diag_field('icebergs', 'u_iceberg', axes, Time, &
     'Iceberg u velocity (m/s)')
  grd%id_v_iceberg=register_diag_field('icebergs', 'v_iceberg', axes, Time, &
     'Iceberg v velocity (m/s)')
  grd%id_spread_uvel=register_diag_field('icebergs', 'spread_uvel', axes, Time, &
     'Iceberg u velocity spread (m/s)')
  grd%id_spread_vvel=register_diag_field('icebergs', 'spread_vvel', axes, Time, &
     'Iceberg v velocity spread (m/s)')
  grd%id_ustar_iceberg=register_diag_field('icebergs', 'ustar_iceberg', axes, Time, &
     'Iceberg frictional velocity (m/s)')
  grd%id_virtual_area=register_diag_field('icebergs', 'virtual_area', axes, Time, &
     'Virtual coverage by icebergs', 'm^2')
  grd%id_mass=register_diag_field('icebergs', 'mass', axes, Time, &
     'Iceberg density field', 'kg/(m^2)')
  grd%id_stored_ice=register_diag_field('icebergs', 'stored_ice', axes3d, Time, &
     'Accumulated ice mass by class', 'kg')
  grd%id_real_calving=register_diag_field('icebergs', 'real_calving', axes3d, Time, &
     'Calving into iceberg class', 'kg/s')
  grd%id_rmean_calving=register_diag_field('icebergs', 'running_mean_calving', axes, Time, &
     'Running mean of calving', 'kg/s')
  grd%id_rmean_calving_hflx=register_diag_field('icebergs', 'running_mean_calving_hflx', axes, Time, &
     'Running mean of calving heat flux', 'J/s')
  grd%id_count=register_diag_field('icebergs', 'bergs_per_cell', axes, Time, &
     'Number of bergs per cell', '#')
  grd%id_chksum=register_diag_field('icebergs', 'list_chksum', axes, Time, &
     'mpp_chksum on bergs in each cell', '#')
  grd%id_uo=register_diag_field('icebergs', 'uo', axes, Time, &
     'Ocean zonal component of velocity', 'm s^-1')
  grd%id_vo=register_diag_field('icebergs', 'vo', axes, Time, &
     'Ocean meridional component of velocity', 'm s^-1')
  grd%id_ui=register_diag_field('icebergs', 'ui', axes, Time, &
     'Ice zonal component of velocity', 'm s^-1')
  grd%id_vi=register_diag_field('icebergs', 'vi', axes, Time, &
     'Ice meridional component of velocity', 'm s^-1')
  grd%id_ua=register_diag_field('icebergs', 'ua', axes, Time, &
     'Atmos zonal component of velocity', 'm s^-1')
  grd%id_va=register_diag_field('icebergs', 'va', axes, Time, &
     'Atmos meridional component of velocity', 'm s^-1')
  grd%id_sst=register_diag_field('icebergs', 'sst', axes, Time, &
     'Sea surface temperature', 'degrees_C')
  grd%id_sss=register_diag_field('icebergs', 'sss', axes, Time, &
     'Sea surface salinity', 'psu')
  grd%id_cn=register_diag_field('icebergs', 'cn', axes, Time, &
     'Sea ice concentration', '(fraction)')
  grd%id_hi=register_diag_field('icebergs', 'hi', axes, Time, &
     'Sea ice thickness', 'm')
  grd%id_ssh=register_diag_field('icebergs', 'ssh', axes, Time, &
     'Sea surface hieght', 'm')
  grd%id_fax=register_diag_field('icebergs', 'taux', axes, Time, &
     'X-stress on ice from atmosphere', 'N m^-2')
  grd%id_fay=register_diag_field('icebergs', 'tauy', axes, Time, &
     'Y-stress on ice from atmosphere', 'N m^-2')
  grd%id_ocean_depth=register_diag_field('icebergs', 'Depth', axes, Time, &
     'Ocean Depth', 'm')

  ! Static fields
  id_class=register_static_field('icebergs', 'lon', axes, &
               'longitude (corners)', 'degrees_E')
  if (id_class>0) lerr=send_data(id_class, grd%lon(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'lat', axes, &
               'latitude (corners)', 'degrees_N')
  if (id_class>0) lerr=send_data(id_class, grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'area', axes, &
               'cell area', 'm^2')
  if (id_class>0) lerr=send_data(id_class, grd%area(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'mask', axes, &
               'wet point mask', 'none')
  if (id_class>0) lerr=send_data(id_class, grd%msk(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'ocean_depth_static', axes, &
               'ocean depth static', 'm')
  if (id_class>0) lerr=send_data(id_class, grd%ocean_depth(grd%isc:grd%iec,grd%jsc:grd%jec))

  if (debug) then
    call grd_chksum2(grd, grd%lon, 'init lon')
    call grd_chksum2(grd, grd%lat, 'init lat')
    call grd_chksum2(grd, grd%lonc, 'init lonc')
    call grd_chksum2(grd, grd%latc, 'init latc')
    call grd_chksum2(grd, grd%area, 'init area')
    call grd_chksum2(grd, grd%msk, 'init msk')
    call grd_chksum2(grd, grd%cos, 'init cos')
    call grd_chksum2(grd, grd%sin, 'init sin')
    call grd_chksum2(grd, grd%ocean_depth, 'init ocean_depth')
  endif

  if (do_unit_tests) then
   if (unitTests(bergs)) call error_mesg('diamonds, icebergs_init', 'Unit tests failed!', FATAL)
  endif

 !write(stderrunit,*) 'diamonds: done'
  call mpp_clock_end(bergs%clock_ini)
  call mpp_clock_end(bergs%clock)

end subroutine ice_bergs_framework_init
! ##############################################################################

subroutine offset_berg_dates(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: Time
! Local variables
type(iceberg), pointer :: this
integer :: iyr, imon, iday, ihr, imin, isec, yr_offset
real :: latest_start_year, berg_start_year
real :: current_time_val
integer :: grdi, grdj

  call get_date(Time, iyr, imon, iday, ihr, imin, isec)
  latest_start_year=iyr-999999.

  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      berg_start_year=float(this%start_year)+this%start_day/367.
      if (berg_start_year>latest_start_year) latest_start_year=berg_start_year
      this=>this%next
    enddo
  enddo ; enddo
  call mpp_max(latest_start_year)

  current_time_val=float(iyr)+yearday(imon, iday, ihr, imin, isec)/367.
  if (latest_start_year<=current_time_val) return ! No conflicts!

  yr_offset=int(latest_start_year+1.)-iyr
  if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,i8,a)') &
    'diamonds: Bergs found with creation dates after model date! Adjusting berg dates by ',yr_offset,' years'
  call bergs_chksum(bergs, 'before adjusting start dates')
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      this%start_year=this%start_year-yr_offset
      this=>this%next
    enddo
  enddo ; enddo
  call bergs_chksum(bergs, 'after adjusting start dates')

end subroutine offset_berg_dates

! #############################################################################

subroutine move_berg_between_cells(bergs)  !Move icebergs onto the correct lists if they have moved from cell to cell.
! Arguments
type(icebergs), pointer :: bergs
type(icebergs_gridded), pointer :: grd => null()
type(iceberg), pointer :: moving_berg => null(), this => null()
integer :: grdi, grdj
logical :: quick
! For convenience
grd=>bergs%grd

do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))

      if ((this%ine.ne.grdi) .or. (this%jne.ne.grdj))  then
        moving_berg=>this
        this=>this%next
        
        !Removing the iceberg from the old list
        if (associated(moving_berg%prev)) then
          moving_berg%prev%next=>moving_berg%next
        else
          bergs%list(grdi,grdj)%first=>moving_berg%next
        endif
        if (associated(moving_berg%next)) moving_berg%next%prev=>moving_berg%prev

        !Inserting the iceberg into the new list 
        call insert_berg_into_list(bergs%list(moving_berg%ine,moving_berg%jne)%first,moving_berg)

        !Clear moving_berg
        moving_berg=>null()

      else
        this=>this%next
      endif
    enddo
enddo ; enddo

end subroutine move_berg_between_cells


! #############################################################################

subroutine update_halo_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: kick_the_bucket, this
integer :: nbergs_to_send_e, nbergs_to_send_w
integer :: nbergs_to_send_n, nbergs_to_send_s
integer :: nbergs_rcvd_from_e, nbergs_rcvd_from_w
integer :: nbergs_rcvd_from_n, nbergs_rcvd_from_s
type(icebergs_gridded), pointer :: grd
integer :: i, nbergs_start, nbergs_end
integer :: stderrunit
integer :: grdi, grdj
integer :: halo_width
integer :: temp1, temp2
real :: current_halo_status
logical :: halo_debugging

halo_width=bergs%grd%halo  
halo_debugging=bergs%halo_debugging  

 ! Get the stderr unit number
   stderrunit = stderr()

 ! For convenience
   grd=>bergs%grd

!For debugging, MP1
if (halo_debugging) then
  do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
        write(stderrunit,*) 'A', this%iceberg_num, mpp_pe(), this%halo_berg, grdi, grdj
      this=>this%next
    enddo
  enddo; enddo
  ! Use when debugging:
  call show_all_bonds(bergs)
endif


! Step 1: Clear the current halos

  call mpp_sync_self()
  do grdj = grd%jsd,grd%jsc-1 ;  do grdi = grd%isd,grd%ied
    call delete_all_bergs_in_list(bergs, grdj, grdi)
  enddo ; enddo

  do grdj = grd%jec+1,grd%jed ;  do grdi = grd%isd,grd%ied
    call delete_all_bergs_in_list(bergs,grdj,grdi)
  enddo ; enddo

  do grdj = grd%jsd,grd%jed ;    do grdi = grd%isd,grd%isc-1
    call delete_all_bergs_in_list(bergs,grdj,grdi)
  enddo ; enddo

  do grdj = grd%jsd,grd%jed ;    do grdi = grd%iec+1,grd%ied
    call delete_all_bergs_in_list(bergs,grdj,grdi)
  enddo ; enddo

  call mpp_sync_self()
!##############################

!For debugging
if (halo_debugging) then
  do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
      write(stderrunit,*) 'B', this%iceberg_num, mpp_pe(), this%halo_berg, grdi, grdj
    this=>this%next
    enddo
  enddo; enddo
endif
  if (debug) then
    nbergs_start=count_bergs(bergs)
  endif

  call mpp_sync_self()
!#######################################################

! Step 2: Updating the halos  - This code is mostly copied from send_to_other_pes


  ! Find number of bergs that headed east/west
  nbergs_to_send_e=0
  nbergs_to_send_w=0
  !Bergs on eastern side of the processor
  do grdj = grd%jsc,grd%jec ; do grdi = grd%iec-halo_width+2,grd%iec  
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
    !write(stderrunit,*)  'sending east', this%iceberg_num, this%ine, this%jne, mpp_pe()
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_e=nbergs_to_send_e+1
        current_halo_status=kick_the_bucket%halo_berg
        kick_the_bucket%halo_berg=1.
        call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_e, nbergs_to_send_e, bergs%max_bonds)
        kick_the_bucket%halo_berg=current_halo_status
    enddo
  enddo; enddo


  !Bergs on the western side of the processor
  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%isc+halo_width-1 
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nbergs_to_send_w=nbergs_to_send_w+1
       current_halo_status=kick_the_bucket%halo_berg
       kick_the_bucket%halo_berg=1.
       call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_w, nbergs_to_send_w, bergs%max_bonds)
       kick_the_bucket%halo_berg=current_halo_status
    enddo 
  enddo; enddo


  ! Send bergs east
  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_e, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
    if (nbergs_to_send_e.gt.0) then
      call mpp_send(bergs%obuffer_e%data, nbergs_to_send_e*buffer_width, grd%pe_E, tag=COMM_TAG_2)
    endif
  endif

  ! Send bergs west
  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_w, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
    if (nbergs_to_send_w.gt.0) then
      call mpp_send(bergs%obuffer_w%data, nbergs_to_send_w*buffer_width, grd%pe_W, tag=COMM_TAG_4)
    endif
  endif

  ! Receive bergs from west
  if (grd%pe_W.ne.NULL_PE) then
    nbergs_rcvd_from_w=-999
    call mpp_recv(nbergs_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
    if (nbergs_rcvd_from_w.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_w.gt.0) then
      call increase_ibuffer(bergs%ibuffer_w, nbergs_rcvd_from_w,buffer_width)
      call mpp_recv(bergs%ibuffer_w%data, nbergs_rcvd_from_w*buffer_width, grd%pe_W, tag=COMM_TAG_2)
      do i=1, nbergs_rcvd_from_w
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_w, i, grd, max_bonds_in=bergs%max_bonds )
      enddo
    endif
  else
    nbergs_rcvd_from_w=0
  endif

  ! Receive bergs from east
  if (grd%pe_E.ne.NULL_PE) then
    nbergs_rcvd_from_e=-999
    call mpp_recv(nbergs_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
    if (nbergs_rcvd_from_e.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_e.gt.0) then
      call increase_ibuffer(bergs%ibuffer_e, nbergs_rcvd_from_e,buffer_width)
      call mpp_recv(bergs%ibuffer_e%data, nbergs_rcvd_from_e*buffer_width, grd%pe_E, tag=COMM_TAG_4)
      do i=1, nbergs_rcvd_from_e
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_e, i, grd, max_bonds_in=bergs%max_bonds )
      enddo
    endif
  else
    nbergs_rcvd_from_e=0
  endif


 ! Find number of bergs that headed north/south
  nbergs_to_send_n=0
  nbergs_to_send_s=0
  

  !Bergs on north side of the processor
  do grdj = grd%jec-halo_width+2,grd%jec ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nbergs_to_send_n=nbergs_to_send_n+1
      current_halo_status=kick_the_bucket%halo_berg
      kick_the_bucket%halo_berg=1.
      call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_n, nbergs_to_send_n, bergs%max_bonds )
      kick_the_bucket%halo_berg=current_halo_status
    enddo
  enddo; enddo


  !Bergs on south side of the processor
  do grdj = grd%jsc,grd%jsc+halo_width-1 ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      kick_the_bucket=>this
      this=>this%next
      nbergs_to_send_s=nbergs_to_send_s+1
       current_halo_status=kick_the_bucket%halo_berg
       kick_the_bucket%halo_berg=1.
       call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_s, nbergs_to_send_s,bergs%max_bonds )
       kick_the_bucket%halo_berg=current_halo_status
    enddo
  enddo; enddo


 ! Send bergs north
  if (grd%pe_N.ne.NULL_PE) then
    if(folded_north_on_pe) then
       call mpp_send(nbergs_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_send(nbergs_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
    endif
    if (nbergs_to_send_n.gt.0) then
       if(folded_north_on_pe) then
          call mpp_send(bergs%obuffer_n%data, nbergs_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
       else
          call mpp_send(bergs%obuffer_n%data, nbergs_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_6)
       endif
    endif
  endif

  ! Send bergs south
  if (grd%pe_S.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_s, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
    if (nbergs_to_send_s.gt.0) then
      call mpp_send(bergs%obuffer_s%data, nbergs_to_send_s*buffer_width, grd%pe_S, tag=COMM_TAG_8)
    endif
  endif


  ! Receive bergs from south
  if (grd%pe_S.ne.NULL_PE) then
    nbergs_rcvd_from_s=-999
    call mpp_recv(nbergs_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
    if (nbergs_rcvd_from_s.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_s.gt.0) then
      call increase_ibuffer(bergs%ibuffer_s, nbergs_rcvd_from_s,buffer_width)
      call mpp_recv(bergs%ibuffer_s%data, nbergs_rcvd_from_s*buffer_width, grd%pe_S, tag=COMM_TAG_6)
      do i=1, nbergs_rcvd_from_s
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_s, i, grd, max_bonds_in=bergs%max_bonds  )
      enddo
    endif
  else
    nbergs_rcvd_from_s=0
  endif

  ! Receive bergs from north
  if (grd%pe_N.ne.NULL_PE) then
    nbergs_rcvd_from_n=-999
    if(folded_north_on_pe) then
       call mpp_recv(nbergs_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_recv(nbergs_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
    endif
    if (nbergs_rcvd_from_n.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_n.gt.0) then
      call increase_ibuffer(bergs%ibuffer_n, nbergs_rcvd_from_n,buffer_width)
      if(folded_north_on_pe) then
         call mpp_recv(bergs%ibuffer_n%data, nbergs_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
         call mpp_recv(bergs%ibuffer_n%data, nbergs_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_8)
      endif
      do i=1, nbergs_rcvd_from_n
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_n, i, grd, max_bonds_in=bergs%max_bonds )
      enddo
    endif
  else
    nbergs_rcvd_from_n=0
  endif



!For debugging
if (halo_debugging) then
  call mpp_sync_self()
  do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      write(stderrunit,*)  'C', this%iceberg_num, mpp_pe(), this%halo_berg,  grdi, grdj
      this=>this%next
    enddo
  enddo; enddo
  call show_all_bonds(bergs)
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
  if (debug) then
    nbergs_end=count_bergs(bergs)
    i=nbergs_rcvd_from_n+nbergs_rcvd_from_s+nbergs_rcvd_from_e+nbergs_rcvd_from_w &
     -nbergs_to_send_n-nbergs_to_send_s-nbergs_to_send_e-nbergs_to_send_w
    if (nbergs_end-(nbergs_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_end=',nbergs_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_start=',nbergs_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: error=',nbergs_end-(nbergs_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_to_send_n=',nbergs_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_to_send_s=',nbergs_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_to_send_e=',nbergs_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_to_send_w=',nbergs_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_rcvd_from_n=',nbergs_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_rcvd_from_s=',nbergs_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_rcvd_from_e=',nbergs_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, update_halos: nbergs_rcvd_from_w=',nbergs_rcvd_from_w,' on PE',mpp_pe()
      !call error_mesg('diamonds, update_halos:', 'We lost some bergs!', FATAL)
    endif
  endif
  if (debug) then
    i=0
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        call check_position(grd, this, 'exchange (bot)')
        if (this%ine.lt.bergs%grd%isc .or. &
            this%ine.gt.bergs%grd%iec .or. &
            this%jne.lt.bergs%grd%jsc .or. &
            this%jne.gt.bergs%grd%jec) i=i+1
        this=>this%next
      enddo ! while
    enddo ; enddo
    call mpp_sum(i)
    if (i>0 .and. mpp_pe()==mpp_root_pe()) then
      write(stderrunit,'(a,i4)') 'diamonds, update_halos: # of bergs outside computational domain = ',i
      call error_mesg('diamonds, update_halos:', 'there are bergs still in halos!', FATAL)
    endif ! root_pe
  endif ! debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

end subroutine update_halo_icebergs



subroutine delete_all_bergs_in_list(bergs,grdj,grdi)
  type(icebergs), pointer :: bergs
  ! Local variables
  type(iceberg), pointer :: kick_the_bucket, this
  integer :: grdi, grdj
  this=>bergs%list(grdi,grdj)%first
  do while (associated(this))
    kick_the_bucket=>this
    this=>this%next
    call destroy_iceberg(kick_the_bucket)
!    call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
  enddo
  bergs%list(grdi,grdj)%first=>null()
end  subroutine delete_all_bergs_in_list


! #############################################################################

subroutine send_bergs_to_other_pes(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: kick_the_bucket, this
integer :: nbergs_to_send_e, nbergs_to_send_w
integer :: nbergs_to_send_n, nbergs_to_send_s
integer :: nbergs_rcvd_from_e, nbergs_rcvd_from_w
integer :: nbergs_rcvd_from_n, nbergs_rcvd_from_s
type(icebergs_gridded), pointer :: grd
integer :: i, nbergs_start, nbergs_end
integer :: stderrunit
integer :: grdi, grdj

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  if (debug) then
    nbergs_start=count_bergs(bergs, with_halos=.true.)
  endif

  ! Find number of bergs that headed east/west
  nbergs_to_send_e=0
  nbergs_to_send_w=0
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%halo_berg .lt. 0.5) then
        if (this%ine.gt.bergs%grd%iec) then
          kick_the_bucket=>this
          this=>this%next
          nbergs_to_send_e=nbergs_to_send_e+1
          call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_e, nbergs_to_send_e, bergs%max_bonds  )
          call move_trajectory(bergs, kick_the_bucket)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
        elseif (this%ine.lt.bergs%grd%isc) then
          kick_the_bucket=>this
          this=>this%next
          nbergs_to_send_w=nbergs_to_send_w+1
          call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_w, nbergs_to_send_w, bergs%max_bonds  )
          call move_trajectory(bergs, kick_the_bucket)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
        else
          this=>this%next
        endif
      else
         this=>this%next
      endif
    enddo
  enddo ; enddo

  ! Send bergs east
  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_e, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
    if (nbergs_to_send_e.gt.0) then
      call mpp_send(bergs%obuffer_e%data, nbergs_to_send_e*buffer_width, grd%pe_E, tag=COMM_TAG_2)
    endif
  endif

  ! Send bergs west
  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_w, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
    if (nbergs_to_send_w.gt.0) then
      call mpp_send(bergs%obuffer_w%data, nbergs_to_send_w*buffer_width, grd%pe_W, tag=COMM_TAG_4)
    endif
  endif

  ! Receive bergs from west
  if (grd%pe_W.ne.NULL_PE) then
    nbergs_rcvd_from_w=-999
    call mpp_recv(nbergs_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
    if (nbergs_rcvd_from_w.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_w.gt.0) then
      call increase_ibuffer(bergs%ibuffer_w, nbergs_rcvd_from_w,buffer_width)
      call mpp_recv(bergs%ibuffer_w%data, nbergs_rcvd_from_w*buffer_width, grd%pe_W, tag=COMM_TAG_2)
      do i=1, nbergs_rcvd_from_w
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_w, i, grd, max_bonds_in=bergs%max_bonds  )
      enddo
    endif
  else
    nbergs_rcvd_from_w=0
  endif

  ! Receive bergs from east
  if (grd%pe_E.ne.NULL_PE) then
    nbergs_rcvd_from_e=-999
    call mpp_recv(nbergs_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
    if (nbergs_rcvd_from_e.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_e.gt.0) then
      call increase_ibuffer(bergs%ibuffer_e, nbergs_rcvd_from_e,buffer_width)
      call mpp_recv(bergs%ibuffer_e%data, nbergs_rcvd_from_e*buffer_width, grd%pe_E, tag=COMM_TAG_4)
      do i=1, nbergs_rcvd_from_e
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_e, i, grd, max_bonds_in=bergs%max_bonds)
      enddo
    endif
  else
    nbergs_rcvd_from_e=0
  endif

  ! Find number of bergs that headed north/south
  ! (note: this block should technically go ahead of the E/W recv block above
  !  to handle arbitrary orientation of PEs. But for simplicity, it is
  !  here to accomodate diagonal transfer of bergs between PEs -AJA)
  nbergs_to_send_n=0
  nbergs_to_send_s=0
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%halo_berg .lt. 0.5) then
        if (this%jne.gt.bergs%grd%jec) then
          kick_the_bucket=>this
          this=>this%next
          nbergs_to_send_n=nbergs_to_send_n+1
          call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_n, nbergs_to_send_n,bergs%max_bonds)
          call move_trajectory(bergs, kick_the_bucket)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
        elseif (this%jne.lt.bergs%grd%jsc) then
          kick_the_bucket=>this
          this=>this%next
          nbergs_to_send_s=nbergs_to_send_s+1
          call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_s, nbergs_to_send_s,bergs%max_bonds)
          call move_trajectory(bergs, kick_the_bucket)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
        else
          this=>this%next
        endif
      else
          this=>this%next
      endif
    enddo
  enddo ; enddo

  ! Send bergs north
  if (grd%pe_N.ne.NULL_PE) then
    if(folded_north_on_pe) then
       call mpp_send(nbergs_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_9)
    else 
       call mpp_send(nbergs_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
    endif
    if (nbergs_to_send_n.gt.0) then
       if(folded_north_on_pe) then
          call mpp_send(bergs%obuffer_n%data, nbergs_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
       else
          call mpp_send(bergs%obuffer_n%data, nbergs_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_6)
       endif
    endif
  endif

  ! Send bergs south
  if (grd%pe_S.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_s, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
    if (nbergs_to_send_s.gt.0) then
      call mpp_send(bergs%obuffer_s%data, nbergs_to_send_s*buffer_width, grd%pe_S, tag=COMM_TAG_8)
    endif
  endif

  ! Receive bergs from south
  if (grd%pe_S.ne.NULL_PE) then
    nbergs_rcvd_from_s=-999
    call mpp_recv(nbergs_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
    if (nbergs_rcvd_from_s.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_s.gt.0) then
      call increase_ibuffer(bergs%ibuffer_s, nbergs_rcvd_from_s,buffer_width)
      call mpp_recv(bergs%ibuffer_s%data, nbergs_rcvd_from_s*buffer_width, grd%pe_S, tag=COMM_TAG_6)
      do i=1, nbergs_rcvd_from_s
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_s, i, grd, max_bonds_in=bergs%max_bonds )
      enddo
    endif
  else
    nbergs_rcvd_from_s=0
  endif

  ! Receive bergs from north
  if (grd%pe_N.ne.NULL_PE) then
    nbergs_rcvd_from_n=-999
    if(folded_north_on_pe) then
       call mpp_recv(nbergs_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_recv(nbergs_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
    endif
    if (nbergs_rcvd_from_n.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_n.gt.0) then
      call increase_ibuffer(bergs%ibuffer_n, nbergs_rcvd_from_n,buffer_width)
      if(folded_north_on_pe) then
         call mpp_recv(bergs%ibuffer_n%data, nbergs_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
         call mpp_recv(bergs%ibuffer_n%data, nbergs_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_8)
      endif
      do i=1, nbergs_rcvd_from_n
        call unpack_berg_from_buffer2(bergs, bergs%ibuffer_n, i, grd, max_bonds_in=bergs%max_bonds)
      enddo
    endif
  else
    nbergs_rcvd_from_n=0
  endif

  if (debug) then
    nbergs_end=count_bergs(bergs, with_halos=.true.)
    i=nbergs_rcvd_from_n+nbergs_rcvd_from_s+nbergs_rcvd_from_e+nbergs_rcvd_from_w &
     -nbergs_to_send_n-nbergs_to_send_s-nbergs_to_send_e-nbergs_to_send_w
    if (nbergs_end-(nbergs_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_end=',nbergs_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_start=',nbergs_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: error=',nbergs_end-(nbergs_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_n=',nbergs_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_s=',nbergs_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_e=',nbergs_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_w=',nbergs_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_n=',nbergs_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_s=',nbergs_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_e=',nbergs_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_w=',nbergs_rcvd_from_w,' on PE',mpp_pe()
      call error_mesg('diamonds, send_bergs_to_other_pes:', 'We lost some bergs!', FATAL)
    endif
  endif

  if (debug) then
    i=0
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        call check_position(grd, this, 'exchange (bot)', grdi, grdj)
        if (this%ine.lt.bergs%grd%isc .or. &
            this%ine.gt.bergs%grd%iec .or. &
            this%jne.lt.bergs%grd%jsc .or. &
            this%jne.gt.bergs%grd%jec) i=i+1
        this=>this%next
      enddo ! while
    enddo ; enddo
    call mpp_sum(i)
    if (i>0 .and. mpp_pe()==mpp_root_pe()) then
      write(stderrunit,'(a,i4)') 'diamonds, send_bergs_to_other_pes: # of bergs outside computational domain = ',i
      call error_mesg('diamonds, send_bergs_to_other_pes:', 'there are bergs still in halos!', FATAL)
    endif ! root_pe
  endif ! debug

  call mpp_sync_self()

end subroutine send_bergs_to_other_pes

  subroutine pack_berg_into_buffer2(berg, buff, n, max_bonds_in)
  ! Arguments
  type(iceberg), pointer :: berg
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  integer, optional :: max_bonds_in
  !integer, intent(in) :: max_bonds  ! Change this later
  ! Local variables
  integer :: counter, k, max_bonds
  type(bond), pointer :: current_bond
 
    max_bonds=0
    if (present(max_bonds_in)) max_bonds=max_bonds_in


    if (.not.associated(buff)) call increase_ibuffer(buff,n,buffer_width)
    if (n>buff%size) call increase_ibuffer(buff,n,buffer_width)

    buff%data(1,n)=berg%lon
    buff%data(2,n)=berg%lat
    buff%data(3,n)=berg%uvel
    buff%data(4,n)=berg%vvel
    buff%data(5,n)=berg%xi
    buff%data(6,n)=berg%yj
    buff%data(7,n)=berg%start_lon
    buff%data(8,n)=berg%start_lat
    buff%data(9,n)=float(berg%start_year)
    buff%data(10,n)=berg%start_day
    buff%data(11,n)=berg%start_mass
    buff%data(12,n)=berg%mass
    buff%data(13,n)=berg%thickness
    buff%data(14,n)=berg%width
    buff%data(15,n)=berg%length
    buff%data(16,n)=berg%mass_scaling
    buff%data(17,n)=berg%mass_of_bits
    buff%data(18,n)=berg%heat_density
    buff%data(19,n)=berg%ine
    buff%data(20,n)=berg%jne
    buff%data(21,n)=berg%axn  !Alon
    buff%data(22,n)=berg%ayn  !Alon
    buff%data(23,n)=berg%bxn  !Alon
    buff%data(24,n)=berg%byn  !Alon
    buff%data(25,n)=float(berg%iceberg_num)
    buff%data(26,n)=berg%halo_berg 
    buff%data(27,n)=berg%static_berg 

    if (max_bonds .gt. 0) then
      counter=27 !how many data points being passed so far (must match above)
      current_bond=>berg%first_bond
      do k = 1,max_bonds
        if (associated(current_bond)) then
          buff%data(counter+(3*(k-1)+1),n)=float(current_bond%other_berg_num) 
          buff%data(counter+(3*(k-1)+2),n)=float(current_bond%other_berg_ine)
          buff%data(counter+(3*(k-1)+3),n)=float(current_bond%other_berg_jne)
          current_bond=>current_bond%next_bond
        else
          buff%data(counter+(3*(k-1)+1),n)=0. 
          buff%data(counter+(3*(k-1)+2),n)=0.
          buff%data(counter+(3*(k-1)+3),n)=0.
        endif
      enddo
    endif
   
    ! Clearing berg pointer from partner bonds
    !if (berg%halo_berg .lt. 0.5) then
    !  call clear_berg_from_partners_bonds(berg)
    !endif

  end subroutine pack_berg_into_buffer2


!###########################################################################3

  subroutine clear_berg_from_partners_bonds(berg) 
  !Arguments
  type(iceberg), intent(in), pointer :: berg
  type(iceberg), pointer :: other_berg
  type(bond), pointer :: current_bond, matching_bond
  integer ::  stderrunit
  ! Get the stderr unit number
  stderrunit = stderr()

    current_bond=>berg%first_bond
    do while (associated(current_bond)) !Looping over bonds
      other_berg=>current_bond%other_berg
      if (associated(other_berg)) then
        !write(stderrunit,*) , 'Other berg', berg%iceberg_num, other_berg%iceberg_num, mpp_pe()
        matching_bond=>other_berg%first_bond
        do while (associated(matching_bond))  ! Looping over possible matching bonds in other_berg
          if (matching_bond%other_berg_num .eq. berg%iceberg_num) then
            !write(stderrunit,*) , 'Clearing', berg%iceberg_num, matching_bond%other_berg_num,other_berg%iceberg_num, mpp_pe()
            matching_bond%other_berg=>null()
            matching_bond=>null()
          else
            matching_bond=>matching_bond%next_bond
          endif
        enddo
      else
       ! Note: This is meant to be unmatched after you have cleared the first berg       
       ! call error_mesg('diamonds, clear berg from partners', 'The bond you are trying to clear is unmatched!', WARNING) 
      endif
      current_bond=>current_bond%next_bond
    enddo !End loop over bonds

  end subroutine clear_berg_from_partners_bonds


  subroutine unpack_berg_from_buffer2(bergs, buff, n,grd, force_append, max_bonds_in)
  ! Arguments
  type(icebergs), pointer :: bergs
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  type(icebergs_gridded), pointer :: grd  
  logical, optional :: force_append
  integer, optional :: max_bonds_in
 ! Local variables
 !real :: lon, lat, uvel, vvel, xi, yj

 !real :: start_lon, start_lat, start_day, start_mass
 !integer :: ine, jne, start_year
  logical :: lres
  type(iceberg) :: localberg
  type(iceberg), pointer :: this
  integer :: other_berg_num, other_berg_ine, other_berg_jne
  integer :: counter, k, max_bonds
  integer :: stderrunit
  logical :: force_app
  logical :: quick

  ! Get the stderr unit number
  stderrunit = stderr()
 
  quick=.false.
  max_bonds=0
  if (present(max_bonds_in)) max_bonds=max_bonds_in

  force_app = .false.
  if(present(force_append)) force_app = force_append
     
    localberg%lon=buff%data(1,n)
    localberg%lat=buff%data(2,n)
    localberg%uvel=buff%data(3,n)
    localberg%vvel=buff%data(4,n)
    localberg%xi=buff%data(5,n)
    localberg%yj=buff%data(6,n)
    localberg%start_lon=buff%data(7,n)
    localberg%start_lat=buff%data(8,n)
    localberg%start_year=nint(buff%data(9,n))
    localberg%start_day=buff%data(10,n)
    localberg%start_mass=buff%data(11,n)
    localberg%mass=buff%data(12,n)
    localberg%thickness=buff%data(13,n)
    localberg%width=buff%data(14,n)
    localberg%length=buff%data(15,n)
    localberg%mass_scaling=buff%data(16,n)
    localberg%mass_of_bits=buff%data(17,n)
    localberg%heat_density=buff%data(18,n)

    localberg%axn=buff%data(21,n) 
    localberg%ayn=buff%data(22,n) 
    localberg%bxn=buff%data(23,n) 
    localberg%byn=buff%data(24,n) 
    localberg%iceberg_num=nint(buff%data(25,n))
    localberg%halo_berg=buff%data(26,n) 
    localberg%static_berg=buff%data(27,n) 
    counter=27 !how many data points being passed so far (must match largest number directly above)

    !These quantities no longer need to be passed between processors
    localberg%uvel_old=localberg%uvel
    localberg%vvel_old=localberg%vvel
    localberg%lon_old=localberg%lon 
    localberg%lat_old=localberg%lat

    ! force_app=.true.
    if(force_app) then !force append with origin ine,jne (for I/O)

      localberg%ine=buff%data(19,n) 
      localberg%jne=buff%data(20,n) 
      call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg,quick,this) 
    else
      lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      if (lres) then
        lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
        call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg,quick,this)
      else
        lres=find_cell_wide(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
        if (lres) then
          lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
          call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg,quick,this)
        else
          write(stderrunit,'("diamonds, unpack_berg_from_buffer pe=(",i3,a,2i4,a,2f8.2)')&
           & mpp_pe(),') Failed to find i,j=',localberg%ine,localberg%jne,' for lon,lat=',localberg%lon,localberg%lat
          write(stderrunit,*) localberg%lon,localberg%lat
          write(stderrunit,*) localberg%uvel,localberg%vvel
          write(stderrunit,*) localberg%axn,localberg%ayn !Alon
          write(stderrunit,*) localberg%bxn,localberg%byn !Alon
          write(stderrunit,*) localberg%uvel_old,localberg%vvel_old 
          write(stderrunit,*) localberg%lon_old,localberg%lat_old 
          write(stderrunit,*) grd%isc,grd%iec,grd%jsc,grd%jec
          write(stderrunit,*) grd%isd,grd%ied,grd%jsd,grd%jed
          write(stderrunit,*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
          write(stderrunit,*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
          write(stderrunit,*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
          write(stderrunit,*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
          write(stderrunit,*) lres
          call error_mesg('diamonds, unpack_berg_from_buffer', 'can not find a cell to place berg in!', FATAL)
        endif
      endif
    endif

    !#  Do stuff to do with bonds here MP1

    this%first_bond=>null()
    if (max_bonds .gt. 0) then
      do k = 1,max_bonds
        other_berg_num=nint(buff%data(counter+(3*(k-1)+1),n))
        other_berg_ine=nint(buff%data(counter+(3*(k-1)+2),n))
        other_berg_jne=nint(buff%data(counter+(3*(k-1)+3),n))
        if (other_berg_num .gt. 0.5) then
          call form_a_bond(this, other_berg_num, other_berg_ine, other_berg_jne)
        endif
      enddo
    endif
    this=>null()

    !##############################

  end subroutine unpack_berg_from_buffer2

  subroutine increase_ibuffer(old,num_bergs,width)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: num_bergs,width
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size, old_size
  !This routine checks if the buffer size is smaller than nbergs
  !If it is, the buffer size is increased by delta_buf
  !The buffer increases by more than 1 so that the buffer does not have to increase every time

    if (.not.associated(old)) then
      new_size=num_bergs+delta_buf
      old_size=0
    else
      old_size=old%size
      if (num_bergs<old%size) then
        new_size=old%size
      else
        new_size=num_bergs+delta_buf
      endif
    endif

    if (old_size.ne.new_size) then
      allocate(new)
      !allocate(new%data(buffer_width,new_size))
      allocate(new%data(width,new_size))
      new%size=new_size
      if (associated(old)) then
        new%data(:,1:old%size)=old%data(:,1:old%size)
        deallocate(old%data)
        deallocate(old)
      endif
      old=>new
     !write(stderr(),*) 'diamonds, increase_ibuffer',mpp_pe(),' increased to',new_size
    endif

  end subroutine increase_ibuffer

  subroutine pack_traj_into_buffer2(traj, buff, n, save_short_traj)
  ! Arguments
  type(xyt), pointer :: traj
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  logical, intent(in) :: save_short_traj
  ! Local variables

    if (.not.associated(buff)) call increase_ibuffer(buff,n,buffer_width_traj)
    if (n>buff%size) call increase_ibuffer(buff,n,buffer_width_traj)

    buff%data(1,n)=traj%lon
    buff%data(2,n)=traj%lat
    buff%data(3,n)=float(traj%year)
    buff%data(4,n)=traj%day
    buff%data(5,n)=float(traj%iceberg_num)
    if (.not. save_short_traj) then
      buff%data(6,n)=traj%uvel
      buff%data(7,n)=traj%vvel
      buff%data(8,n)=traj%mass
      buff%data(9,n)=traj%mass_of_bits
      buff%data(10,n)=traj%heat_density
      buff%data(11,n)=traj%thickness
      buff%data(12,n)=traj%width
      buff%data(13,n)=traj%length
      buff%data(14,n)=traj%uo
      buff%data(15,n)=traj%vo
      buff%data(16,n)=traj%ui
      buff%data(17,n)=traj%vi
      buff%data(18,n)=traj%ua
      buff%data(19,n)=traj%va
      buff%data(20,n)=traj%ssh_x
      buff%data(21,n)=traj%ssh_y
      buff%data(22,n)=traj%sst
      buff%data(23,n)=traj%cn
      buff%data(24,n)=traj%hi
      buff%data(25,n)=traj%axn !Alon
      buff%data(26,n)=traj%ayn !Alon
      buff%data(27,n)=traj%bxn !Alon
      buff%data(28,n)=traj%byn !Alon
      buff%data(29,n)=traj%halo_berg !Alon
      buff%data(30,n)=traj%static_berg !Alon
      buff%data(31,n)=traj%sss
    endif

  end subroutine pack_traj_into_buffer2

  subroutine unpack_traj_from_buffer2(first, buff, n, save_short_traj)
  ! Arguments
  type(xyt), pointer :: first
  type(buffer), pointer :: buff
  integer, intent(in) :: n
 ! Local variables
  type(xyt) :: traj
  integer :: stderrunit
  logical, intent(in) :: save_short_traj 
  ! Get the stderr unit number
  stderrunit = stderr()

    traj%lon=buff%data(1,n)
    traj%lat=buff%data(2,n)
    traj%year=nint(buff%data(3,n))
    traj%day=buff%data(4,n)
    traj%iceberg_num=nint(buff%data(5,n))
    if (.not. save_short_traj) then
      traj%uvel=buff%data(6,n)
      traj%vvel=buff%data(7,n)
      traj%mass=buff%data(8,n)
      traj%mass_of_bits=buff%data(9,n)
      traj%heat_density=buff%data(10,n)
      traj%thickness=buff%data(11,n)
      traj%width=buff%data(12,n)
      traj%length=buff%data(13,n)
      traj%uo=buff%data(14,n)
      traj%vo=buff%data(15,n)
      traj%ui=buff%data(16,n)
      traj%vi=buff%data(17,n)
      traj%ua=buff%data(18,n)
      traj%va=buff%data(19,n)
      traj%ssh_x=buff%data(20,n)
      traj%ssh_y=buff%data(21,n)
      traj%sst=buff%data(22,n)
      traj%cn=buff%data(23,n)
      traj%hi=buff%data(24,n)
      traj%axn=buff%data(25,n) !Alon
      traj%ayn=buff%data(26,n) !Alon
      traj%bxn=buff%data(27,n) !Alon
      traj%byn=buff%data(28,n) !Alon
      traj%halo_berg=buff%data(29,n) !Alon
      traj%static_berg=buff%data(30,n) !Alon
      traj%sss=buff%data(31,n)
    endif
    call append_posn(first, traj)

  end subroutine unpack_traj_from_buffer2


! ##############################################################################

subroutine add_new_berg_to_list(first, bergvals, quick, newberg_return)
! Arguments
type(iceberg), pointer :: first
type(iceberg), intent(in) :: bergvals
type(iceberg), intent(out), pointer, optional :: newberg_return
logical, intent(in), optional :: quick
! Local variables
type(iceberg), pointer :: new=>null()

  new=>null()
  call create_iceberg(new, bergvals)

  if (present(newberg_return)) then
    newberg_return=>new
    !newberg_return=>null()
  endif

  if (present(quick)) then
    if(quick) then
      call insert_berg_into_list(first, new, quick=.true.)
    else
      call insert_berg_into_list(first, new)
    endif
  else
    call insert_berg_into_list(first, new)
  endif

  !Clear new
  new=>null()

end subroutine add_new_berg_to_list

! ##############################################################################

subroutine count_out_of_order(bergs,label)
! Arguments
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
type(iceberg), pointer :: this, next
integer :: i, icnt1, icnt2, icnt3
integer :: grdi, grdj

  icnt1=0; icnt3=0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    next=>null()
    if (associated(this)) then
      if (associated(this%next)) next=>this%next
    endif
    do while (associated(next))
      if (.not. inorder(this,next)) icnt1=icnt1+1
      if (inorder(this,next).and.inorder(next,this)) icnt3=icnt3+1
      this=>next
      next=>next%next
    enddo
  enddo;enddo
  call mpp_sum(icnt1)

  i=0;
  icnt2=0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      i=1
      if (this%ine<bergs%grd%isc .or. &
          this%ine>bergs%grd%iec .or. &
          this%jne<bergs%grd%jsc .or. &
          this%jne>bergs%grd%jec) icnt2=icnt2+1
      this=>this%next
      if (i>1.and..not.associated(this%prev)) then
        call error_mesg('diamonds, count_out_of_order', 'Pointer %prev is unassociated. This should not happen!', FATAL)
      endif
    enddo
  enddo; enddo
  call mpp_sum(icnt2)

  if ((debug.or.icnt1.ne.0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,3(x,a,i6),x,a)') 'diamonds, count_out_of_order:', &
      '# out of order=', icnt1,'# in halo=',icnt2,'# identicals=',icnt3,label
  endif

  call check_for_duplicates(bergs,label)

end subroutine count_out_of_order

! ##############################################################################

subroutine check_for_duplicates(bergs,label)
! Arguments
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
type(iceberg), pointer :: this1, next1, this2, next2
integer :: icnt_id, icnt_same
integer :: grdi, grdj
integer :: grdi_inner, grdj_inner

  icnt_id=0
  icnt_same=0
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this1=>bergs%list(grdi,grdj)%first
    do while (associated(this1))
      do grdj_inner = grdj,bergs%grd%jec ; do grdi_inner = bergs%grd%isc,bergs%grd%iec
        if ( .not.(grdj_inner==grdj .and. grdi_inner<grdi) ) then
          if (grdj_inner==grdj .and. grdi_inner==grdi ) then
            this2=>this1%next
          else
            this2=>bergs%list(grdi_inner,grdj_inner)%first
          endif
          do while (associated(this2))
            if (sameid(this1,this2)) icnt_id=icnt_id+1
            if (sameberg(this1,this2)) icnt_same=icnt_same+1
            this2=>this2%next
          enddo
        endif
      enddo ; enddo
      this1=>this1%next
    enddo
  enddo ; enddo
  call mpp_sum(icnt_id)
  call mpp_sum(icnt_same)

  if ((debug.or.icnt_id>0.or.icnt_same>0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,2(x,a,i9),x,a)') 'diamonds, check_for_duplicates:', &
      '# with same id=', icnt_id,'# identical bergs=',icnt_same,label
  endif

end subroutine check_for_duplicates

! ##############################################################################

subroutine monitor_a_berg(bergs, label)
! Arguments
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
type(iceberg), pointer :: this
integer :: grdi, grdj
integer :: stderrunit

  if (bergs%debug_iceberg_with_id<0) return
  stderrunit=stderr() ! Get the stderr unit number

  do grdj = bergs%grd%jsd,bergs%grd%jed ; do grdi = bergs%grd%isd,bergs%grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%iceberg_num == bergs%debug_iceberg_with_id) then
        call print_berg(stderrunit, this, 'MONITOR: '//label, grdi, grdj)
      endif
      this=>this%next
    enddo
  enddo ; enddo

end subroutine monitor_a_berg

! ##############################################################################

subroutine insert_berg_into_list(first, newberg, quick)
! Arguments
type(iceberg), pointer :: first, newberg
logical, intent(in), optional :: quick
! Local variables
type(iceberg), pointer :: this, prev
logical :: quickly

  quickly = .false.

  if (associated(first)) then
    if (.not. parallel_reprod .or. quickly) then
      newberg%next=>first
      newberg%prev=>null()
      first%prev=>newberg
      first=>newberg
    else
      if (inorder(newberg,first)) then
        ! Insert at front of list
        newberg%next=>first
        newberg%prev=>null()
        first%prev=>newberg
        first=>newberg
      else
        this=>first
        prev=>null()
        do while( associated(this) )
          if (inorder(newberg,this) ) then
            exit
          endif
          prev=>this
          this=>this%next
        enddo
        prev%next=>newberg
        newberg%prev=>prev
        if (associated(this)) this%prev=>newberg
        newberg%next=>this
      endif
    endif
  else
    ! list is empty so create it
    first=>newberg
    first%next=>null()
    first%prev=>null()
  endif

end subroutine insert_berg_into_list

! ##############################################################################

logical function inorder(berg1, berg2)  !MP Alon - Change to include iceberg_num
! Arguments
type(iceberg), pointer :: berg1, berg2
! Local variables
  if (berg1%start_year<berg2%start_year) then ! want newer first
    inorder=.true.
    return
  else if (berg1%start_year>berg2%start_year) then
    inorder=.false.
    return
  endif
  if (berg1%start_day<berg2%start_day) then ! want newer first
    inorder=.true.
    return
  else if (berg1%start_day>berg2%start_day) then
    inorder=.false.
    return
  endif
  if (berg1%start_mass<berg2%start_mass) then ! want lightest first
    inorder=.true.
    return
  else if (berg1%start_mass>berg2%start_mass) then
    inorder=.false.
    return
  endif
  if (berg1%start_lon<berg2%start_lon) then ! want eastward first
    inorder=.true.
    return
  else if (berg1%start_lon>berg2%start_lon) then
    inorder=.false.
    return
  endif
  if (berg1%start_lat<berg2%start_lat) then ! want southern first
    inorder=.true.
    return
  else if (berg1%start_lat>berg2%start_lat) then
    inorder=.false.
    return
  endif
  inorder=.true. ! passing the above tests mean the bergs 1 and 2 are identical?
end function inorder

! ##############################################################################

  real function time_hash(berg)!  Alon: Think about removing this.
  ! Arguments
  type(iceberg), pointer :: berg
    time_hash=berg%start_day+366.*float(berg%start_year)
  end function time_hash

! ##############################################################################

  real function pos_hash(berg)
  ! Arguments
  type(iceberg), pointer :: berg
    pos_hash=berg%start_lon+360.*(berg%start_lat+90.)
  end function pos_hash

! ##############################################################################

logical function sameid(berg1, berg2) !  Alon: MP updat this.
! Arguments
type(iceberg), pointer :: berg1, berg2
! Local variables
  sameid=.false.
  if (berg1%start_year.ne.berg2%start_year) return
  if (berg1%start_day.ne.berg2%start_day) return
  if (berg1%start_mass.ne.berg2%start_mass) return
  if (berg1%start_lon.ne.berg2%start_lon) return
  if (berg1%start_lat.ne.berg2%start_lat) return
  sameid=.true. ! passing the above tests means that bergs 1 and 2 have the same id
end function sameid

! ##############################################################################

logical function sameberg(berg1, berg2)
! Arguments
type(iceberg), pointer :: berg1, berg2
! Local variables
  sameberg=.false.
  if (.not. sameid(berg1, berg2)) return
  if (berg1%lon.ne.berg2%lon) return
  if (berg1%lat.ne.berg2%lat) return
  if (berg1%mass.ne.berg2%mass) return
  if (berg1%uvel.ne.berg2%uvel) return
  if (berg1%vvel.ne.berg2%vvel) return
  if (berg1%thickness.ne.berg2%thickness) return
  if (berg1%width.ne.berg2%width) return
  if (berg1%length.ne.berg2%length) return
  if (berg1%axn.ne.berg2%axn) return  !Alon
  if (berg1%ayn.ne.berg2%ayn) return  !Alon
  if (berg1%bxn.ne.berg2%bxn) return  !Alon
  if (berg1%byn.ne.berg2%byn) return  !Alon
  if (berg1%uvel_old.ne.berg2%uvel_old) return  !Alon
  if (berg1%vvel_old.ne.berg2%vvel_old) return  !Alon
  if (berg1%lon_old .ne.berg2%lon_old) return  !Alon
  if (berg1%lat_old.ne.berg2%lat_old) return  !Alon
  sameberg=.true. ! passing the above tests mean that bergs 1 and 2 are identical
end function sameberg

! ##############################################################################

real function yearday(imon, iday, ihr, imin, isec)
! Arguments
integer, intent(in) :: imon, iday, ihr, imin, isec

  yearday=float(imon-1)*31.+float(iday-1)+(float(ihr)+(float(imin)+float(isec)/60.)/60.)/24.

end function yearday

! ##############################################################################

subroutine create_iceberg(berg, bergvals)
! Arguments
type(iceberg), pointer :: berg
type(iceberg), intent(in) :: bergvals
! Local variables
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  if (associated(berg)) then
    write(stderrunit,*) 'diamonds, create_iceberg: berg already associated!!!!',mpp_pe()
    call error_mesg('diamonds, create_iceberg', 'berg already associated. This should not happen!', FATAL)
  endif
  allocate(berg)
  berg=bergvals
  berg%prev=>null()
  berg%next=>null()

end subroutine create_iceberg


! ##############################################################################

subroutine delete_iceberg_from_list(first, berg)
! Arguments
type(iceberg), pointer :: first, berg
! Local variables

  ! Connect neighbors to each other
  if (associated(berg%prev)) then
    berg%prev%next=>berg%next
  else
    first=>berg%next
  endif
  if (associated(berg%next)) berg%next%prev=>berg%prev

  ! Bye-bye berg
  call destroy_iceberg(berg)

end subroutine delete_iceberg_from_list

! ##############################################################################

subroutine destroy_iceberg(berg)
! Arguments
type(iceberg), pointer :: berg
! Local variables

  ! Clears all matching bonds before deallocint memory
  call clear_berg_from_partners_bonds(berg)

  ! Bye-bye berg
  deallocate(berg)

end subroutine destroy_iceberg

! ##############################################################################

subroutine print_berg(iochan, berg, label, il, jl)
! Arguments
integer, intent(in) :: iochan
type(iceberg), pointer :: berg
character(len=*) :: label
integer, optional, intent(in) :: il, jl !< Indices of cell berg should be in
! Local variables

  write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,a,2f10.4,i5,f7.2,es12.4,f5.1)') &
    label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, ' start lon,lat,yr,day,mass,hb=', &
    berg%start_lon, berg%start_lat, berg%start_year, berg%start_day, berg%start_mass, berg%halo_berg
  if (present(il).and.present(jl)) then
    write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,a,2i5)') &
      label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, ' List i,j=',il,jl
  endif
  write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,a,2i5,a,2l2)') &
    label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, &
    ' i,j=', berg%ine, berg%jne, &
    ' p,n=', associated(berg%prev), associated(berg%next)
  write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,3(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, &
    ' xi,yj=', berg%xi, berg%yj, &
    ' lon,lat=', berg%lon, berg%lat, &
    ' lon_old,lat_old=', berg%lon_old, berg%lat_old
  write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,2(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, &
    ' u,v=', berg%uvel, berg%vvel, &
    ' uvel_old,vvel_old=', berg%uvel_old, berg%vvel_old
  write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,2(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, &
    ' axn,ayn=', berg%axn, berg%ayn, &
    ' bxn,byn=', berg%bxn, berg%byn
  write(iochan,'("diamonds, print_berg: ",2a,i5,a,i12,3(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%iceberg_num, &
    ' uo,vo=', berg%uo, berg%vo, &
    ' ua,va=', berg%ua, berg%va, &
    ' ui,vi=', berg%ui, berg%vi
end subroutine print_berg

! ##############################################################################

subroutine print_bergs(iochan, bergs, label)
! Arguments
integer, intent(in) :: iochan
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
integer :: nbergs, nnbergs
type(iceberg), pointer :: this
integer :: grdi, grdj

  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))
      call print_berg(iochan, this, label)
      this=>this%next
    enddo
  enddo ; enddo
  nbergs=count_bergs(bergs)
  nnbergs=nbergs
  call mpp_sum(nnbergs)
  if (nbergs.gt.0) write(iochan,'("diamonds, ",a," there are",i5," bergs out of",i6," on PE ",i4)') label, nbergs, nnbergs, mpp_pe()

end subroutine print_bergs


! ##############################################################################

subroutine form_a_bond(berg, other_berg_num, other_berg_ine, other_berg_jne, other_berg)

type(iceberg), pointer :: berg
type(iceberg), optional,  pointer :: other_berg
type(bond) , pointer :: new_bond, first_bond
integer, intent(in) :: other_berg_num
integer, optional  :: other_berg_ine, other_berg_jne
integer :: stderrunit

 stderrunit = stderr()
    
if (berg%iceberg_num .ne. other_berg_num) then
              
 !write (stderrunit,*) , 'Forming a bond!!!', mpp_pe(), berg%iceberg_num, other_berg_num, berg%halo_berg, berg%ine, berg%jne

  ! Step 1: Create a new bond
  allocate(new_bond)
  new_bond%other_berg_num=other_berg_num
  if(present(other_berg)) then
    new_bond%other_berg=>other_berg
    new_bond%other_berg_ine=other_berg%ine
    new_bond%other_berg_jne=other_berg%jne
  else
    new_bond%other_berg=>null()
    if (present(other_berg_ine)) then
      new_bond%other_berg_ine=other_berg_ine
      new_bond%other_berg_jne=other_berg_jne
    endif
  endif

  ! Step 2: Put this new bond at the start of the bond list
   first_bond=>berg%first_bond
   if (associated(first_bond)) then
        new_bond%next_bond=>first_bond
        new_bond%prev_bond=>null()  !This should not be needed
        first_bond%prev_bond=>new_bond
        berg%first_bond=>new_bond
   else
       new_bond%next_bond=>null()  !This should not be needed
       new_bond%prev_bond=>null()  !This should not be needed
       berg%first_bond=>new_bond
   endif
   new_bond=>null()
 else
   call error_mesg('diamonds, bonds', 'An iceberg is trying to bond with itself!!!', FATAL)
 endif

end subroutine form_a_bond

! #############################################################################

subroutine bond_address_update(bergs)
type(icebergs), pointer :: bergs
type(iceberg), pointer :: other_berg, berg
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj, nbonds
type(bond) , pointer :: current_bond

 ! For convenience
  grd=>bergs%grd

  ! This could be done for only the bergs near the halos
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        if  (associated(current_bond%other_berg)) then
          current_bond%other_berg_ine=current_bond%other_berg%ine
          current_bond%other_berg_jne=current_bond%other_berg%jne
        else
          if (berg%halo_berg .lt. 0.5) then     
            call error_mesg('diamonds, bond address update', 'other berg in bond not assosiated!', FATAL)
          endif
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo; enddo

  call mpp_sync_self()

end subroutine bond_address_update

!###################################################################################################

subroutine show_all_bonds(bergs)
type(icebergs), pointer :: bergs
type(iceberg), pointer :: other_berg, berg
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj, nbonds
type(bond) , pointer :: current_bond

 ! For convenience
  grd=>bergs%grd

  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        print *, 'Show Bond1 :', berg%iceberg_num, current_bond%other_berg_num, current_bond%other_berg_ine, current_bond%other_berg_jne,  mpp_pe()
        !print *, 'Current:', berg%iceberg_num, berg%ine, berg%jne,berg%halo_berg, mpp_pe()
        if  (associated(current_bond%other_berg)) then
          if (current_bond%other_berg%iceberg_num .ne. current_bond%other_berg_num) then
            print *, 'Bond matching', berg%iceberg_num,current_bond%other_berg%iceberg_num, current_bond%other_berg_num,&
            berg%halo_berg,current_bond%other_berg%halo_berg ,mpp_pe()
            call error_mesg('diamonds, show all bonds:', 'The bonds are not matching properly!', FATAL)
          endif
        else
            print *, 'This bond has an non-assosiated other berg :', berg%iceberg_num, current_bond%other_berg_num,&
            current_bond%other_berg_ine, current_bond%other_berg_jne, berg%halo_berg,  mpp_pe()
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo; enddo

end subroutine show_all_bonds


subroutine connect_all_bonds(bergs)
type(icebergs), pointer :: bergs
type(iceberg), pointer :: other_berg, berg
type(icebergs_gridded), pointer :: grd
integer :: i, j
integer :: grdi, grdj
integer :: grdi_inner, grdj_inner
type(bond) , pointer :: current_bond, other_berg_bond
logical :: bond_matched, missing_bond, check_bond_quality
integer nbonds

missing_bond=.false.
bond_matched=.false.

! For convenience
  grd=>bergs%grd

  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
! do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec  ! Don't connect halo bergs
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        !code to find parter bond goes here
        if (.not.associated(current_bond%other_berg)) then
          bond_matched=.false.
          i = current_bond%other_berg_ine ; j = current_bond%other_berg_jne
          if ( (i.gt. grd%isd-1) .and. (i .lt. grd%ied+1) .and. (j .gt. grd%jsd-1) .and. (j .lt. grd%jed+1)) then
            other_berg=>bergs%list(i,j)%first
            do while (associated(other_berg)) ! loop over all other bergs
              if (other_berg%iceberg_num .eq. current_bond%other_berg_num) then
                current_bond%other_berg=>other_berg
                other_berg=>null()
                bond_matched=.true. 
              else
                other_berg=>other_berg%next
              endif 
            enddo
          endif
          if (.not.bond_matched) then
            ! If you are stil not matched, then search adjacent cells
            do grdj_inner = j-1,j+1 ; do grdi_inner = i-1,i+1
              if (.not. bond_matched) then
                if    ((grdj_inner .gt. grd%jsd-1) .and. (grdj_inner .lt. grd%jed+1)        &
                .and.  (grdi_inner .gt. grd%isd-1) .and. (grdi_inner .lt. grd%ied+1)      &
                .and. ((grdi_inner .ne. i) .or. (grdj_inner .ne. j)) ) then
                  other_berg=>bergs%list(grdi_inner,grdj_inner)%first
                  do while (associated(other_berg)) ! loop over all other bergs
                    if (other_berg%iceberg_num .eq. current_bond%other_berg_num) then
                      current_bond%other_berg=>other_berg
                      other_berg=>null()
                      bond_matched=.true.  
                    else
                      other_berg=>other_berg%next
                    endif 
                  enddo
                endif
              endif
            enddo;enddo
          endif
          if (.not.bond_matched) then      
            if (berg%halo_berg .lt. 0.5) then
              missing_bond=.true.    
              print * ,'non-halo berg unmatched: ', berg%iceberg_num, mpp_pe(), current_bond%other_berg_num, current_bond%other_berg_ine 
              call error_mesg('diamonds, connect_all_bonds', 'A non-halo bond is missing!!!', FATAL)
            else  ! This is not a problem if the partner berg is not yet in the halo
              !if (  (current_bond%other_berg_ine .gt.grd%isd-1) .and. (current_bond%other_berg_ine .lt.grd%ied+1) &
                !.and.  (current_bond%other_berg_jne .gt.grd%jsd-1) .and. (current_bond%other_berg_jne .lt.grd%jed+1) ) then      
                !print * ,'halo berg unmatched: ',mpp_pe(),  berg%iceberg_num, current_bond%other_berg_num, current_bond%other_berg_ine,current_bond%other_berg_jne 
                !call error_mesg('diamonds, connect_all_bonds', 'A halo bond is missing!!!', WARNING)
              !endif
            endif
          endif
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo;enddo

  if (debug) then
    check_bond_quality=.true.
    nbonds=0
    call count_bonds(bergs, nbonds,check_bond_quality)
  endif
end subroutine connect_all_bonds


! #############################################################################
subroutine count_bonds(bergs, number_of_bonds, check_bond_quality)

type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
type(iceberg), pointer :: other_berg
type(icebergs_gridded), pointer :: grd
type(bond) , pointer :: current_bond, other_berg_bond
integer, intent(out) :: number_of_bonds
integer :: number_of_bonds_all_pe
integer :: grdi, grdj
logical :: bond_is_good
logical, intent(inout), optional :: check_bond_quality
logical :: quality_check
integer :: num_unmatched_bonds,num_unmatched_bonds_all_pe
integer :: num_unassosiated_bond_pairs, num_unassosiated_bond_pairs_all_pe 
integer :: stderrunit

!  print *, "starting bond_check"
 stderrunit = stderr()
 quality_check=.false.
 if(present(check_bond_quality)) quality_check = check_bond_quality
 check_bond_quality=.false.
 num_unmatched_bonds=0
 num_unassosiated_bond_pairs=0

 ! For convenience
  grd=>bergs%grd

  number_of_bonds=0  ! This is a bond counter.
  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs

      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        number_of_bonds=number_of_bonds+1

        ! ##### Beginning Quality Check on Bonds ######
        !      print *, 'Quality check', mpp_pe(), berg%iceberg_num
        if (quality_check) then
          num_unmatched_bonds=0
          num_unassosiated_bond_pairs=0
          bond_is_good=.False.
          other_berg=>current_bond%other_berg
          if (associated(other_berg)) then
            other_berg_bond=>other_berg%first_bond
            do while (associated(other_berg_bond))  !loops over the icebergs in the other icebergs bond list           
              if (associated(other_berg_bond%other_berg)) then
                if (other_berg_bond%other_berg%iceberg_num .eq.berg%iceberg_num) then
                  bond_is_good=.True.  !Bond_is_good becomes true when the corresponding bond is found
                endif
              endif                 
              if (bond_is_good) then
                other_berg_bond=>null()
              else
                other_berg_bond=>other_berg_bond%next_bond
              endif
            enddo  ! End of loop over the other berg's bonds.

            if (bond_is_good) then
              if (debug) write(stderrunit,*) 'Perfect quality Bond:', berg%iceberg_num, current_bond%other_berg_num
            else  
              if (debug) write(stderrunit,*) 'Non-matching bond...:', berg%iceberg_num, current_bond%other_berg_num
              num_unmatched_bonds=num_unmatched_bonds+1
            endif
          else
            if (debug) write(stderrunit,*) 'Opposite berg is not assosiated:', berg%iceberg_num, current_bond%other_berg%iceberg_num    
            num_unassosiated_bond_pairs=0
          endif
        endif
        ! ##### Ending Quality Check on Bonds ######
        current_bond=>current_bond%next_bond
      enddo !End of loop over current bonds
      berg=>berg%next
    enddo ! End of loop over all bergs
  enddo; enddo !End of loop over all grid cells

    number_of_bonds_all_pe=number_of_bonds
    call mpp_sum(number_of_bonds_all_pe)

    bergs%nbonds=number_of_bonds_all_pe !Total number of bonds across all pe's
    if (debug) then 
      if (number_of_bonds .gt. 0) then
        write(stderrunit,*) "Bonds on PE:",number_of_bonds, "Total bonds", number_of_bonds_all_PE, "on PE number:",  mpp_pe()
      endif
    endif

    if (quality_check) then
      num_unmatched_bonds_all_pe=num_unmatched_bonds
      num_unassosiated_bond_pairs_all_pe = num_unassosiated_bond_pairs
      call mpp_sum(num_unmatched_bonds_all_pe)
      call mpp_sum(num_unassosiated_bond_pairs_all_pe)

      if (num_unmatched_bonds_all_pe .gt. 0) then
        call error_mesg('diamonds, bonds', 'Bonds are not matching!', FATAL)
      endif
      if (num_unassosiated_bond_pairs_all_pe .ne. 0) then
        call error_mesg('diamonds, bonds', 'Bonds partners not located!', Warning)
        if (num_unassosiated_bond_pairs .ne. 0) then
          write(*,'(2a)') 'diamonds, Bonds parnters not located!!!! PE=', mpp_pe()
        endif
      endif
      if ((num_unmatched_bonds_all_pe .eq. 0)  .and. (num_unassosiated_bond_pairs_all_pe .eq. 0)) then
        if (mpp_pe().eq.mpp_root_pe()) then
                write(stderrunit,*)  "Total number of bonds is: ", number_of_bonds_all_PE, "All iceberg bonds are connected and working well"
        endif
        check_bond_quality=.true.
      else
        if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'diamonds: Warning, Broken Bonds! '
      endif
    endif

end subroutine count_bonds



! ##############################################################################

integer function count_bergs(bergs, with_halos)
! Arguments
type(icebergs), pointer :: bergs
logical, optional :: with_halos
! Local variables
integer :: grdi, grdj, is, ie, js, je
logical :: include_halos

  include_halos = .false.
  if (present(with_halos)) include_halos = with_halos
  if (include_halos) then
   is = bergs%grd%isd ; ie = bergs%grd%ied ; js = bergs%grd%jsd ; je = bergs%grd%jed
  else
   is = bergs%grd%isc ; ie = bergs%grd%iec ; js = bergs%grd%jsc ; je = bergs%grd%jec
  endif

  count_bergs=0
  do grdj = js,je ; do grdi = is,ie
    count_bergs=count_bergs+count_bergs_in_list(bergs%list(grdi,grdj)%first)
  enddo ; enddo

end function count_bergs

! ##############################################################################

integer function count_bergs_in_list(first)
! Arguments
type(iceberg), pointer :: first
! Local variables
type(iceberg), pointer :: this

  count_bergs_in_list=0
  this=>first
  do while(associated(this))
    count_bergs_in_list=count_bergs_in_list+1
    this=>this%next
  enddo

end function count_bergs_in_list

! ##############################################################################

subroutine record_posn(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(xyt) :: posn
type(iceberg), pointer :: this
integer :: grdi, grdj

  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      posn%lon=this%lon
      posn%lat=this%lat
      posn%year=bergs%current_year
      posn%day=bergs%current_yearday
      posn%iceberg_num=posn%iceberg_num
      if (.not. bergs%save_short_traj) then !Not totally sure that this is correct
        posn%uvel=this%uvel
        posn%vvel=this%vvel
        posn%mass=this%mass
        posn%mass_of_bits=this%mass_of_bits
        posn%heat_density=this%heat_density
        posn%thickness=this%thickness
        posn%width=this%width
        posn%length=this%length
        posn%uo=this%uo
        posn%vo=this%vo
        posn%ui=this%ui
        posn%vi=this%vi
        posn%ua=this%ua
        posn%va=this%va
        posn%ssh_x=this%ssh_x
        posn%ssh_y=this%ssh_y
        posn%sst=this%sst
        posn%sss=this%sss
        posn%cn=this%cn
        posn%hi=this%hi
        posn%axn=this%axn
        posn%ayn=this%ayn
        posn%bxn=this%bxn
        posn%byn=this%byn
        posn%uvel_old=this%uvel_old
        posn%vvel_old=this%vvel_old
        posn%lon_old=this%lon_old
        posn%lat_old=this%lat_old
        posn%halo_berg=this%halo_berg
        posn%static_berg=this%static_berg
      endif
  
      call push_posn(this%trajectory, posn)
  
      this=>this%next
    enddo
  enddo ; enddo

end subroutine record_posn

! ##############################################################################

subroutine push_posn(trajectory, posn_vals)
! Arguments
type(xyt), pointer :: trajectory
type(xyt) :: posn_vals
! Local variables
type(xyt), pointer :: new_posn

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>trajectory
  trajectory=>new_posn

end subroutine push_posn

subroutine append_posn(trajectory, posn_vals)
! This routine appends a new position leaf to the end of the given trajectory 
! Arguments
type(xyt), pointer :: trajectory
type(xyt) :: posn_vals
! Local variables
type(xyt), pointer :: new_posn,next,last

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>null()
  if(.NOT. associated(trajectory)) then
     trajectory=>new_posn
  else
     ! Find end of the trajectory and point it to the  new leaf
     next=>trajectory
     do while (associated(next))
        last=>next
        next=>next%next
     enddo
     last%next=>new_posn
  endif
end subroutine append_posn

! ##############################################################################

subroutine move_trajectory(bergs, berg)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
! Local variables
type(xyt), pointer :: next, last
type(xyt) :: vals

  if (bergs%ignore_traj) return

  ! If the trajectory is empty, ignore it
  if (.not.associated(berg%trajectory)) return

  ! Push identifying info into first posn (note reverse order due to stack)
  vals%lon=berg%start_lon
  vals%lat=berg%start_lat
  vals%year=berg%start_year
  vals%iceberg_num=berg%iceberg_num
  vals%day=berg%start_day
  vals%mass=berg%start_mass
  call push_posn(berg%trajectory, vals)
  vals%lon=0.
  vals%lat=99.
  vals%year=0
  vals%day=0.
  vals%mass=berg%mass_scaling
  call push_posn(berg%trajectory, vals)

  ! Find end of berg trajectory and point it to start of existing trajectories
  next=>berg%trajectory
  do while (associated(next))
    last=>next
    next=>next%next
  enddo
  last%next=>bergs%trajectories

  bergs%trajectories=>berg%trajectory
  berg%trajectory=>null()

end subroutine move_trajectory

! ##############################################################################

subroutine move_all_trajectories(bergs, delete_bergs)
! Arguments
type(icebergs),    pointer    :: bergs
logical, optional, intent(in) :: delete_bergs
! Local variables
type(iceberg), pointer :: this, next
logical :: delete_bergs_after_moving_traj
integer :: grdi, grdj
  
  if (bergs%ignore_traj) return

  delete_bergs_after_moving_traj = .false.
  if (present(delete_bergs)) delete_bergs_after_moving_traj = delete_bergs
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      next=>this%next
      call move_trajectory(bergs, this)
   !  if (delete_bergs_after_moving_traj) call destroy_iceberg(this)
      this=>next
    enddo
  enddo ; enddo

end subroutine move_all_trajectories

! ##############################################################################

logical function find_cell_by_search(grd, x, y, i, j)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: x, y
integer, intent(inout) :: i, j
! Local variables
integer :: is,ie,js,je,di,dj,io,jo,icnt
real :: d0,d1,d2,d3,d4,d5,d6,d7,d8,dmin
logical :: explain=.false.
real :: Lx
 
911 continue

  Lx=grd%Lx  
  find_cell_by_search=.false.
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec

  ! Start at nearest corner
  d1=dcost(x,y,grd%lonc(is+1,js+1),grd%latc(is+1,js+1),Lx)
  d2=dcost(x,y,grd%lonc(ie-1,js+1),grd%latc(ie-1,js+1),Lx)
  d3=dcost(x,y,grd%lonc(ie-1,je-1),grd%latc(ie-1,je-1),Lx)
  d4=dcost(x,y,grd%lonc(is+1,je-1),grd%latc(is+1,je-1),Lx)
  dmin=min(d1,d2,d3,d4)
  if (d1==dmin) then; i=is+1; j=js+1
  elseif (d2==dmin) then; i=ie-1; j=js+1
  elseif (d3==dmin) then; i=ie-1; j=je-1
  elseif (d4==dmin) then; i=is+1; j=je-1
  else
    call error_mesg('diamonds, find_cell_by_search:', 'This should never EVER happen! (1)', FATAL)
  endif

  if (explain) then
    write(0,'(i3,a,2i4,f9.4)') mpp_pe(),'Initial corner i-is,j-js,cost=',i-is,j-js,dmin
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d4,d3
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d1,d2
  endif

  if (is_point_in_cell(grd, x, y, i, j)) then
    find_cell_by_search=.true.
    return
  endif
    
  do icnt=1, 1*(ie-is+je-js)
    io=i; jo=j

    d0=dcost(x,y,grd%lonc(io,jo),grd%latc(io,jo),Lx)
    d1=dcost(x,y,grd%lonc(io,jo+1),grd%latc(io,jo+1),Lx)
    d2=dcost(x,y,grd%lonc(io-1,jo+1),grd%latc(io-1,jo+1),Lx)
    d3=dcost(x,y,grd%lonc(io-1,jo),grd%latc(io-1,jo),Lx)
    d4=dcost(x,y,grd%lonc(io-1,jo-1),grd%latc(io-1,jo-1),Lx)
    d5=dcost(x,y,grd%lonc(io,jo-1),grd%latc(io,jo-1),Lx)
    d6=dcost(x,y,grd%lonc(io+1,jo-1),grd%latc(io+1,jo-1),Lx)
    d7=dcost(x,y,grd%lonc(io+1,jo),grd%latc(io+1,jo),Lx)
    d8=dcost(x,y,grd%lonc(io+1,jo+1),grd%latc(io+1,jo+1),Lx)

  ! dmin=min(d0,d1,d3,d5,d7)
    dmin=min(d0,d1,d2,d3,d4,d5,d6,d7,d8)
    if (d0==dmin) then; di=0; dj=0
    elseif (d2==dmin) then; di=-1; dj=1
    elseif (d4==dmin) then; di=-1; dj=-1
    elseif (d6==dmin) then; di=1; dj=-1
    elseif (d8==dmin) then; di=1; dj=1
    elseif (d1==dmin) then; di=0; dj=1
    elseif (d3==dmin) then; di=-1; dj=0
    elseif (d5==dmin) then; di=0; dj=-1
    elseif (d7==dmin) then; di=1; dj=0
    else
      call error_mesg('diamonds, find_cell_by_search:', 'This should never EVER happen!', FATAL)
    endif

    i=min(ie, max(is, io+di))
    j=min(je, max(js, jo+dj))

    if (explain) then
      write(0,'(i3,a,2i4,f9.5,a,2i4,a,2i4)') mpp_pe(),'Current position i,j,cost=',i,j,dmin,' di,dj=',di,dj,' old io,jo=',io,jo
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    endif

    if (is_point_in_cell(grd, x, y, i, j)) then
      find_cell_by_search=.true.
      return
    endif
    
    if ((i==io.and.j==jo) &
        .and. .not.find_better_min(grd, x, y, 3, i, j) &
       ) then
      ! Stagnated
      find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 1, i, j)
      if (.not. find_cell_by_search) find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 3, i, j)
      i=min(ie, max(is, i))
      j=min(je, max(js, j))
      if (is_point_in_cell(grd, x, y, i, j)) then
        find_cell_by_search=.true.
      else
  !     find_cell_by_search=find_cell(grd, x, y, io, jo)
  !     if (find_cell_by_search) then
  !       if (explain) then
  !         call print_fld(grd, grd%lat, 'Lat')
  !         call print_fld(grd, grd%lon, 'Lat')
  !         do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
  !           grd%tmp(i,j) = dcost(x,y,grd%lonc(i,j),grd%latc(i,j))
  !         enddo; enddo
  !         call print_fld(grd, grd%tmp, 'Cost')
  !         stop 'Avoid recursing'
  !       endif
  !       write(0,'(i3,a,2i5,a,2i3,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative io,jo=',io,jo,' di,dj=',di,dj,' targ=',x,y
  !       explain=.true.; goto 911
  !     endif
      endif
      return
    endif

  enddo

  find_cell_by_search=find_cell(grd, x, y, i, j)
  if (find_cell_by_search) then
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    write(0,'(i3,a,2f9.5)') mpp_pe(),'x,y ',x,y
    write(0,'(i3,a,4i5)') mpp_pe(),'io,jo ',io,jo,di,dj
    write(0,'(i3,a,2i5,a,2i3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 i,j=',i-is,j-js,' di,dj=',di,dj
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 io,jo=',io,jo
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 i,j=',i,j,' targ=',x,y
    return
  endif
  find_cell_by_search=.false.

  contains

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  real function dcost(x1, y1, x2, y2,Lx)
  ! Arguments
  real, intent(in) :: x1, x2, y1, y2,Lx
  ! Local variables
  real :: x1m

    x1m=apply_modulo_around_point(x1,x2,Lx)
  ! dcost=(x2-x1)**2+(y2-y1)**2
    dcost=(x2-x1m)**2+(y2-y1)**2
  end function dcost

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  logical function find_better_min(grd, x, y, w, oi, oj)
  ! Arguments
  type(icebergs_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye
  real :: dmin, dcst
  real :: Lx

  Lx=grd%Lx

  xs=max(grd%isc, oi-w)
  xe=min(grd%iec, oi+w)
  ys=max(grd%jsc, oj-w)
  ye=min(grd%jec, oj+w)

  find_better_min=.false.
  dmin=dcost(x,y,grd%lonc(oi,oj),grd%latc(oi,oj),Lx)
  do j=ys,ye; do i=xs,xe
      dcst=dcost(x,y,grd%lonc(i,j),grd%latc(i,j),Lx)
      if (dcst<dmin) then
        find_better_min=.true.
        dmin=dcst
        oi=i;oj=j
      endif
  enddo; enddo

  end function find_better_min

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  logical function find_cell_loc(grd, x, y, is, ie, js, je, w, oi, oj)
  ! Arguments
  type(icebergs_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: is, ie, js, je, w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye

    xs=max(is, oi-w)
    xe=min(ie, oi+w)
    ys=max(js, oj-w)
    ye=min(je, oj+w)

    find_cell_loc=.false.
    do j=ys,ye; do i=xs,xe
        if (is_point_in_cell(grd, x, y, i, j)) then
          oi=i; oj=j; find_cell_loc=.true.
          return
        endif
    enddo; enddo

  end function find_cell_loc

end function find_cell_by_search


! ##############################################################################

subroutine find_individual_iceberg(bergs,iceberg_num, ine, jne, berg_found, search_data_domain)
type(icebergs), pointer :: bergs
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
integer, intent(in) :: iceberg_num
logical, intent(in) :: search_data_domain
integer, intent(out) :: ine, jne
real, intent(out) :: berg_found
integer :: ilim1, ilim2, jlim1, jlim2

berg_found=0.0
ine=999
jne=999
  ! For convenience
    grd=>bergs%grd
    
    if (search_data_domain) then
        ilim1 = grd%isd  ; ilim2=grd%ied  ; jlim1 = grd%jsd  ; jlim2=grd%jed
    else
        ilim1 = grd%isc  ; ilim2=grd%iec  ; jlim1 = grd%jsc  ; jlim2=grd%jec
    endif


    do grdj = jlim1, jlim2 ; do grdi = ilim1, ilim2
    !do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    !do grdj = bergs%grd%jsd,bergs%grd%jed ; do grdi = bergs%grd%isd,bergs%grd%ied
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        if (iceberg_num .eq. this%iceberg_num) then
          ine=this%ine
          jne=this%jne
          berg_found=1.0
          !print *, 'found this one'
          return
        endif
        this=>this%next
      enddo
    enddo ; enddo                                                                             
end subroutine  find_individual_iceberg 


! ##############################################################################

logical function find_cell(grd, x, y, oi, oj) 
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell=.false.; oi=-999; oj=-999

  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell=.true.
        return
      endif
  enddo; enddo

end function find_cell

! ##############################################################################

logical function find_cell_wide(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell_wide=.false.; oi=-999; oj=-999

  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell_wide=.true.
        return
      endif
  enddo; enddo

end function find_cell_wide

! ##############################################################################

logical function is_point_in_cell(grd, x, y, i, j, explain)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
logical, intent(in), optional :: explain
! Local variables
real :: xlo, xhi, ylo, yhi
integer :: stderrunit
real :: Lx
real :: tol

  ! Get the stderr unit number
  stderrunit=stderr()
  Lx=grd%Lx

  ! Safety check index bounds
  if (i-1.lt.grd%isd.or.i.gt.grd%ied.or.j-1.lt.grd%jsd.or.j.gt.grd%jed) then
    write(stderrunit,'(a,i3,(a,3i4))') &
                     'diamonds, is_point_in_cell: pe=(',mpp_pe(),') i,s,e=', &
                     i,grd%isd,grd%ied,' j,s,e=', j,grd%jsd,grd%jed
    call error_mesg('diamonds, is_point_in_cell', 'test is off the PE!', FATAL)
  endif

  is_point_in_cell=.false.

  ! Test crude bounds
  xlo=min( apply_modulo_around_point(grd%lon(i-1,j-1)   ,x,  Lx), & 
           apply_modulo_around_point(grd%lon(i  ,j-1)   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i-1,j  )   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i  ,j  )   ,x,  Lx) )
  xhi=max( apply_modulo_around_point(grd%lon(i-1,j-1)   ,x,  Lx), & 
           apply_modulo_around_point(grd%lon(i  ,j-1)   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i-1,j  )   ,x,  Lx), &
           apply_modulo_around_point(grd%lon(i  ,j  )   ,x,  Lx) )

  ! The modolo function inside sum_sign_dot_prod leads to a roundoff.
  !Adding adding a tolorance to the crude bounds avoids excluding the cell which
  !would be correct after roundoff. This is a bit of a hack.
  tol=0.1 
  if (x.lt.(xlo-tol) .or. x.gt.(xhi+tol)) return

  ylo=min( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  yhi=max( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  if (y.lt.ylo .or. y.gt.yhi) return
  
  if ((grd%lat(i,j).gt.89.999).and. (grd%grid_is_latlon))   then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, Lx,explain=explain) 
  elseif ((grd%lat(i-1,j).gt.89.999) .and. (grd%grid_is_latlon))  then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i  ,j  ),grd%lat(i-1,j  ), &
                                        grd%lon(i-1,j-1),grd%lat(i-1,j  ), &
                                        x, y,Lx, explain=explain) 
  elseif ((grd%lat(i-1,j-1).gt.89.999) .and. (grd%grid_is_latlon))  then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j  ),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y,Lx, explain=explain) 
  elseif ((grd%lat(i,j-1).gt.89.999) .and. (grd%grid_is_latlon)) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i-1,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, Lx,explain=explain) 
  else
  is_point_in_cell=sum_sign_dot_prod4(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y,Lx, explain=explain) 
  endif

end function is_point_in_cell

! ##############################################################################

logical function sum_sign_dot_prod4(x0, y0, x1, y1, x2, y2, x3, y3, x, y,Lx, explain)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x, y, Lx
logical, intent(in), optional :: explain
! Local variables
real :: p0,p1,p2,p3,xx
real :: l0,l1,l2,l3
real :: xx0,xx1,xx2,xx3
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  sum_sign_dot_prod4=.false.
  xx= apply_modulo_around_point(x,x0,Lx)
  xx0= apply_modulo_around_point(x0,x0,Lx)
  xx1= apply_modulo_around_point(x1,x0,Lx)
  xx2= apply_modulo_around_point(x2,x0,Lx)
  xx3= apply_modulo_around_point(x3,x0,Lx)


  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y0-y3)-(y-y3)*(xx0-xx3)

  !We use an assymerty between South and East line boundaries and North and East
  !to avoid icebergs appearing to two cells (half values used for debugging)
  !This is intended to make the South and East boundaries be part of the
  !cell, while the North and West are not part of the cell.
  p0=sign(1., l0); if (l0.eq.0.) p0=-0.5
  p1=sign(1., l1); if (l1.eq.0.) p1=0.5
  p2=sign(1., l2); if (l2.eq.0.) p2=0.5
  p3=sign(1., l3); if (l3.eq.0.) p3=-0.5

  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    sum_sign_dot_prod4=.true.
  endif


  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: x=',mpp_pe(),':', &
                           x0,x1,x2,x3, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: y=',mpp_pe(),':', &
                           y0,y1,y2,y3, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: l=',mpp_pe(),':', &
                           l0,l1,l2,l3
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: p=',mpp_pe(),':', &
                           p0,p1,p2,p3, abs( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) - abs((p0+p2)+(p1+p3)) )
   endif
  endif

end function sum_sign_dot_prod4

! ##############################################################################

logical function sum_sign_dot_prod5(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y, Lx, explain)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y, Lx
logical, intent(in), optional :: explain
! Local variables
real :: p0,p1,p2,p3,p4,xx
real :: l0,l1,l2,l3,l4
real :: xx0,xx1,xx2,xx3,xx4
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  sum_sign_dot_prod5=.false.
  xx= apply_modulo_around_point(x,x0,Lx)
  xx0= apply_modulo_around_point(x0,x0,Lx)
  xx1= apply_modulo_around_point(x1,x0,Lx)
  xx2= apply_modulo_around_point(x2,x0,Lx)
  xx3= apply_modulo_around_point(x3,x0,Lx)
  xx4= apply_modulo_around_point(x4,x0,Lx)


  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y4-y3)-(y-y3)*(xx4-xx3)
  l4=(xx-xx4)*(y0-y4)-(y-y4)*(xx0-xx4)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.
  p4=sign(1., l4); if (l4.eq.0.) p4=0.

  if ( ((abs(p0)+abs(p2))+(abs(p1)+abs(p3)))+abs(p4) - abs(((p0+p2)+(p1+p3))+p4) .lt. 0.5 ) then
    sum_sign_dot_prod5=.true.
  endif

  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: x=',mpp_pe(),':', &
                           x0,x1,x2,x3,x4, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3,xx4, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: y=',mpp_pe(),':', &
                           y0,y1,y2,y3,y4, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod5: l=',mpp_pe(),':', &
                           l0,l1,l2,l3,l4
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: p=',mpp_pe(),':', &
                           p0,p1,p2,p3,p4
   endif
  endif

end function sum_sign_dot_prod5

! ##############################################################################

logical function pos_within_cell(grd, x, y, i, j, xi, yj, explain)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
real, intent(out) :: xi, yj
logical, intent(in), optional :: explain
! Local variables
real :: x1,y1,x2,y2,x3,y3,x4,y4,xx,yy,fac
integer :: stderrunit
real :: Lx, dx,dy
real :: Delta_x
logical :: is_point_in_cell_using_xi_yj

  ! Get the stderr unit number
  stderrunit=stderr()
  Lx=grd%Lx
  pos_within_cell=.false.; xi=-999.; yj=-999.
  if (i-1<grd%isd) return
  if (j-1<grd%jsd) return
  if (i>grd%ied) return
  if (j>grd%jed) return

  x1=grd%lon(i-1,j-1)
  y1=grd%lat(i-1,j-1)
  x2=grd%lon(i  ,j-1)
  y2=grd%lat(i  ,j-1)
  x3=grd%lon(i  ,j  )
  y3=grd%lat(i  ,j  )
  x4=grd%lon(i-1,j  )
  y4=grd%lat(i-1,j  )

  if (present(explain)) then
    if(explain) then
    write(stderrunit,'(a,4(f12.6,a))') 'pos_within_cell: lon=[',x1,',',x2,',',x3,',',x4,']'
    write(stderrunit,'(a,4(f12.6,a))') 'pos_within_cell: lat=[',y1,',',y2,',',y3,',',y4,']'
    write(stderrunit,'(2(a,f12.6))') 'pos_within_cell: x,y=',x,',',y
    endif
  endif

  !This part only works for a regular cartesian grid. For more complex grids, we
  !should use calc_xiyj
  if ((.not. grd%grid_is_latlon) .and. (grd%grid_is_regular))  then
    dx=abs((grd%lon(i  ,j  )-grd%lon(i-1  ,j  )))
    dy=abs((grd%lat(i  ,j  )-grd%lat(i  ,j-1  )))
    x1=grd%lon(i  ,j  )-(dx/2)
    y1=grd%lat(i  ,j  )-(dy/2)

    Delta_x= apply_modulo_around_point(x,x1,Lx)-x1

    xi=((Delta_x)/dx)+0.5
    !xi=((x-x1)/dx)+0.5
    yj=((y-y1)/dy)+0.5

  elseif ((max(y1,y2,y3,y4)<89.999) .or.(.not. grd%grid_is_latlon)) then
    ! This returns non-dimensional position xi,yj for quad cells (not at a pole)
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, Lx, explain=explain)
  else
    ! One of the cell corners is at the north pole so we switch to a tangent plane with
    ! co-latitude as a radial coordinate.
    xx=(90.-y)*cos(x*pi_180)
    yy=(90.-y)*sin(x*pi_180)
    x1=(90.-y1)*cos(grd%lon(i-1,j-1)*pi_180)
    y1=(90.-y1)*sin(grd%lon(i-1,j-1)*pi_180)
    x2=(90.-y2)*cos(grd%lon(i  ,j-1)*pi_180)
    y2=(90.-y2)*sin(grd%lon(i  ,j-1)*pi_180)
    x3=(90.-y3)*cos(grd%lon(i  ,j  )*pi_180)
    y3=(90.-y3)*sin(grd%lon(i  ,j  )*pi_180)
    x4=(90.-y4)*cos(grd%lon(i-1,j  )*pi_180)
    y4=(90.-y4)*sin(grd%lon(i-1,j  )*pi_180)
    if (present(explain)) then
      if(explain) then
      write(stderrunit,'(a,4(f12.6,a))') 'pos_within_cell: lon=[',x1,',',x2,',',x3,',',x4,']'
      write(stderrunit,'(a,4(f12.6,a))') 'pos_within_cell: lat=[',y1,',',y2,',',y3,',',y4,']'
      write(stderrunit,'(2(a,f12.6))') 'pos_within_cell: x,y=',xx,',',yy
      endif
    endif
    ! Calculate non-dimensional position xi,yj within a quad in the tangent plane.
    ! This quad has straight sides in the plane and so is not the same on the
    ! projection of the spherical quad.
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, xx, yy, xi, yj, Lx,explain=explain)
    if (is_point_in_cell(grd, x, y, i, j)) then
      ! The point is within the spherical quad
      is_point_in_cell_using_xi_yj=is_point_within_xi_yj_bounds(xi,yj)
      if (.not. is_point_in_cell_using_xi_yj) then
        ! If the non-dimensional position is found to be outside (possible because the
        ! projection of the spherical quad and the quad in the tangent plane are different)
        ! then scale non-dimensional coordinates to be consistent with is_point_in_cell()
        ! Note: this is intended to fix the inconsistency between the tangent plane
        ! and lat-lon calculations and is a work around only for the four polar cells.
        fac=2.1*max( abs(xi-0.5), abs(yj-0.5) ); fac=max(1., fac)
        xi=0.5+(xi-0.5)/fac
        yj=0.5+(yj-0.5)/fac
        if (debug) call error_mesg('diamonds, pos_within_cell', 'in cell but scaling internal coordinates!', WARNING)
      endif
    else
      ! The point is not inside the spherical quad
      if (abs(xi-0.5)<0.5.and.abs(yj-0.5)<0.5) then
        ! The projection of the spherical quad onto the tangent plane should be larger than
        ! quad in the tangent plane so we should never be able to get here.
        call error_mesg('diamonds, pos_within_cell', 'not in cell but coordinates <0.5!', FATAL)
      endif
    endif
  endif

  if (present(explain)) then
    if(explain) write(stderrunit,'(a,2f12.6)') 'pos_within_cell: xi,yj=',xi,yj
  endif

  ! Check for consistency with test for whether point is inside a polygon
  pos_within_cell=is_point_in_cell(grd, x, y, i, j,explain=explain)
  is_point_in_cell_using_xi_yj=is_point_within_xi_yj_bounds(xi,yj)
  if (is_point_in_cell_using_xi_yj) then
    ! Based on coordinate, the point is out of cell
    if (.not. pos_within_cell) then
      ! Based on is_point_in_cell() the point is within cell so we have an inconsistency
      if (debug) then
        write(stderrunit,'(a,1p6e12.4)') 'values of xi, yj ',xi, yj
        pos_within_cell=is_point_in_cell(grd, x, y, i, j,explain=.True.)
        call error_mesg('diamonds, pos_within_cell', 'pos_within_cell is False BUT is_point_in_cell disagrees!', FATAL)
      endif
    endif
  else
    ! Based on coordinate, the point is within cell
    if (pos_within_cell) then
      ! Based on is_point_in_cell() the point is out of cell so we have an inconsistency
      if (debug) then
        write(stderrunit,'(a,1p6e12.4)') 'values of xi, yj ',xi, yj
        pos_within_cell=is_point_in_cell(grd, x, y, i, j,explain=.True.)
        call error_mesg('diamonds, pos_within_cell', 'pos_within_cell is True BUT is_point_in_cell disagrees!', FATAL)
      endif
    endif
  endif

  contains

  subroutine calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj,Lx, explain)
  ! Arguments
  real,  intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4, x, y, Lx
  real, intent(out) :: xi, yj
  logical, intent(in), optional :: explain
  ! Local variables
  real :: alpha, beta, gamma, delta, epsilon, kappa, a, b, c, d, dx, dy, yy1, yy2
  logical :: expl=.false.
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  expl=.false.
  if (present(explain)) then
     if(explain) expl=.true.
  endif

  alpha=x2-x1
  delta=y2-y1
  beta=x4-x1
  epsilon=y4-y1
  gamma=(x3-x1)-(alpha+beta)
  kappa=(y3-y1)-(delta+epsilon)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs alpha,beta,gamma',alpha,beta,gamma,delta,epsilon,kappa
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs delta,epsilon,kappa',alpha,beta,gamma,delta,epsilon,kappa

  a=(kappa*beta-gamma*epsilon)
  dx= apply_modulo_around_point(x,x1,Lx)-x1

  dy=y-y1
  b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
  c=(alpha*dy-delta*dx)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs dx,dy=',dx,dy
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs A,B,C=',a,b,c
  if (abs(a)>1.e-12) then
    d=0.25*(b**2)-a*c
    if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs D=',d
    if (d.ge.0.) then
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for b/2a, sqrt(d) = ',-0.5*b/a,sqrt(d)/a
      yy1=-(0.5*b+sqrt(d))/a
      yy2=-(0.5*b-sqrt(d))/a
      if (abs(yy1-0.5).lt.abs(yy2-0.5)) then; yj=yy1; else; yj=yy2; endif
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for y = ',yy1,yy2,yj
    else
      write(stderrunit,'(a,i3,a,4(f8.2,a))') 'calc_xiyj: ',mpp_pe(),'lon=[',x1,',',x2,',',x3,',',x4,']'
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,a,4(f8.2,a))') 'calc_xiyj: ',mpp_pe(),'lat=[',y1,',',y2,',',y3,',',y4,']'
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3)') 'calc_xiyj: b<0 in quadratic root solver!!!!',mpp_pe()
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs a,b,c,d,dx,dy',mpp_pe(),a,b,c,d,dx,dy
      call error_mesg('diamonds, calc_xiyj', 'We have complex roots. The grid must be very distorted!', FATAL)
    endif
  else
    if (b.ne.0.) then
      yj=-c/b
    else
      yj=0.
    endif
  endif

  a=(alpha+gamma*yj)
  b=(delta+kappa*yj)
  if (a.ne.0.) then
    xi=(dx-beta*yj)/a
  elseif (b.ne.0.) then
    xi=(dy-epsilon*yj)/b
  else
    c=(epsilon*alpha-beta*delta)+(epsilon*gamma-beta*kappa)*yj
    if (c.ne.0.) then
      xi=(epsilon*dx-beta*dy)/c
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3,1p2e12.4)') 'calc_xiyj: coeffs a,b',mpp_pe(),a,b
      call error_mesg('diamonds, calc_xiyj', 'Can not invert either linear equaton for xi! This should not happen!', FATAL)
    endif
  endif
  if (expl) write(stderrunit,'(a,2e12.4)') 'calc_xiyj: xi,yj=',xi,yj

  end subroutine calc_xiyj

end function pos_within_cell

! ##############################################################################

logical function is_point_within_xi_yj_bounds(xi,yj)
! Arguments
real, intent(in) :: xi, yj
! Local variables
!Includes South and East boundaries, and excludes North and West  (double check this is the way that is needed)
  is_point_within_xi_yj_bounds=.False.
  if ((xi .ge. 0 )  .and.  (xi .lt. 1)) then
    if ((yj .ge. 0 )  .and.  (yj .lt. 1)) then
      is_point_within_xi_yj_bounds=.True.
    endif
  endif
end function is_point_within_xi_yj_bounds

real function apply_modulo_around_point(x,y,Lx)
! Arguments
real, intent(in) :: x ,y ,Lx
!Local_variables
real ::Lx_2
!Gives the modula value of x in an interval [y-(Lx/2)  y+(Lx/2)]  , modulo Lx
!If Lx<=0, then it returns x without applying modulo arithmetic.

  if (Lx>0.) then
    Lx_2=Lx/2.
    apply_modulo_around_point=modulo(x-(y-Lx_2),Lx)+(y-Lx_2)
  else
    apply_modulo_around_point=x
  endif

end function apply_modulo_around_point

subroutine check_position(grd, berg, label, il, jl)
! Arguments
type(icebergs_gridded), pointer :: grd
type(iceberg), pointer :: berg
character(len=*) :: label
integer, optional, intent(in) :: il, jl !< Indices of cell berg should be in
! Local variables
real :: xi, yj
logical :: lret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  lret=pos_within_cell(grd, berg%lon, berg%lat, berg%ine, berg%jne, xi, yj)
  if (xi.ne.berg%xi.or.yj.ne.berg%yj) then
    write(stderrunit,'("diamonds: check_position (",i4,") b%x,x,-=",3(es12.4,x),a)') mpp_pe(),berg%xi,xi,berg%xi-xi,label
    write(stderrunit,'("diamonds: check_position (",i4,") b%y,y,-=",3(es12.4,x),a)') mpp_pe(),berg%yj,yj,berg%yj-yj,label
    call print_berg(stderrunit, berg, 'check_position', il, jl)
    call error_mesg('diamonds, check_position, '//trim(label),'berg has inconsistent xi,yj!',FATAL)
  endif
  if (grd%msk(berg%ine, berg%jne)==0.) then
    call print_berg(stderrunit, berg, 'check_position, '//trim(label), il, jl)
    call error_mesg('diamonds, check_position, '//trim(label),'berg is in a land cell!',FATAL)
  endif

end subroutine check_position

! ##############################################################################

real function sum_mass(bergs,justbits,justbergs)
! Arguments
type(icebergs), pointer :: bergs
logical, intent(in), optional :: justbits, justbergs
! Local variables
type(iceberg), pointer :: this
integer :: grdi, grdj

  sum_mass=0.
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))
      if (present(justbergs)) then
        sum_mass=sum_mass+this%mass*this%mass_scaling
      elseif (present(justbits)) then
        sum_mass=sum_mass+this%mass_of_bits*this%mass_scaling
      else
        sum_mass=sum_mass+(this%mass+this%mass_of_bits)*this%mass_scaling
      endif
      this=>this%next
    enddo
  enddo ; enddo

end function sum_mass

! ##############################################################################

real function sum_heat(bergs,justbits,justbergs)
! Arguments
type(icebergs), pointer :: bergs
logical, intent(in), optional :: justbits, justbergs
! Local variables
type(iceberg), pointer :: this
real :: dm
integer :: grdi, grdj

  sum_heat=0.
  do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while(associated(this))
      dm=0.
      if (present(justbergs)) then
        dm=this%mass*this%mass_scaling
      elseif (present(justbits)) then
        dm=this%mass_of_bits*this%mass_scaling
      else
        dm=(this%mass+this%mass_of_bits)*this%mass_scaling
      endif
      sum_heat=sum_heat+dm*this%heat_density
      this=>this%next
    enddo
  enddo ; enddo

end function sum_heat


subroutine sanitize_field(arr,val)
! Arguments
real, dimension(:,:),intent(inout) :: arr
real, intent(in) :: val
! Local variables
integer :: i, j

  do j=lbound(arr,2), ubound(arr,2)
    do i=lbound(arr,1), ubound(arr,1)
      if (abs(arr(i,j)).ge.val) arr(i,j)=0.
    enddo
  enddo

end subroutine sanitize_field

! ##############################################################################




! ##############################################################################

subroutine checksum_gridded(grd, label)
! Arguments
type(icebergs_gridded), pointer :: grd
character(len=*) :: label
! Local variables

  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'diamonds: checksumming gridded data @ ',trim(label)

  ! external forcing
  call grd_chksum2(grd, grd%uo, 'uo')
  call grd_chksum2(grd, grd%vo, 'vo')
  call grd_chksum2(grd, grd%ua, 'ua')
  call grd_chksum2(grd, grd%va, 'va')
  call grd_chksum2(grd, grd%ui, 'ui')
  call grd_chksum2(grd, grd%vi, 'vi')
  call grd_chksum2(grd, grd%ssh, 'ssh')
  call grd_chksum2(grd, grd%sst, 'sst')
  call grd_chksum2(grd, grd%sss, 'sss')
  call grd_chksum2(grd, grd%hi, 'hi')
  call grd_chksum2(grd, grd%cn, 'cn')
  call grd_chksum2(grd, grd%calving, 'calving')
  call grd_chksum2(grd, grd%calving_hflx, 'calving_hflx')

  ! state
  call grd_chksum2(grd, grd%mass, 'mass')
  call grd_chksum3(grd, grd%mass_on_ocean, 'mass_on_ocean')
  call grd_chksum3(grd, grd%area_on_ocean, 'area_on_ocean')
  call grd_chksum3(grd, grd%Uvel_on_ocean, 'Uvel_on_ocean')
  call grd_chksum3(grd, grd%Vvel_on_ocean, 'Vvel_on_ocean')
  call grd_chksum3(grd, grd%stored_ice, 'stored_ice')
  call grd_chksum2(grd, grd%rmean_calving, 'rmean_calving')
  call grd_chksum2(grd, grd%rmean_calving_hflx, 'rmean_calving_hflx')
  call grd_chksum2(grd, grd%stored_heat, 'stored_heat')
  call grd_chksum2(grd, grd%melt_buoy, 'melt_b')
  call grd_chksum2(grd, grd%melt_eros, 'melt_e')
  call grd_chksum2(grd, grd%melt_conv, 'melt_v')
  call grd_chksum2(grd, grd%bergy_src, 'bergy_src')
  call grd_chksum2(grd, grd%bergy_melt, 'bergy_melt')
  call grd_chksum2(grd, grd%bergy_mass, 'bergy_mass')
  call grd_chksum2(grd, grd%spread_mass, 'spread_mass')
  call grd_chksum2(grd, grd%spread_area, 'spread_area')
  call grd_chksum2(grd, grd%u_iceberg, 'u_iceberg')
  call grd_chksum2(grd, grd%v_iceberg, 'v_iceberg')
  call grd_chksum2(grd, grd%spread_uvel, 'spread_uvel')
  call grd_chksum2(grd, grd%spread_vvel, 'spread_vvel')
  call grd_chksum2(grd, grd%ustar_iceberg, 'ustar_iceberg')
  call grd_chksum2(grd, grd%virtual_area, 'varea')
  call grd_chksum2(grd, grd%floating_melt, 'floating_melt')
  call grd_chksum2(grd, grd%berg_melt, 'berg_melt')

  ! static
  call grd_chksum2(grd, grd%lon, 'lon')
  call grd_chksum2(grd, grd%lat, 'lat')
  call grd_chksum2(grd, grd%lonc, 'lonc')
  call grd_chksum2(grd, grd%latc, 'latc')
  call grd_chksum2(grd, grd%dx, 'dx')
  call grd_chksum2(grd, grd%dy, 'dy')
  call grd_chksum2(grd, grd%area, 'area')
  call grd_chksum2(grd, grd%msk, 'msk')
  call grd_chksum2(grd, grd%cos, 'cos')
  call grd_chksum2(grd, grd%sin, 'sin')
  call grd_chksum2(grd, grd%ocean_depth, 'depth')

end subroutine checksum_gridded

! ##############################################################################

subroutine grd_chksum3(grd, fld, txt)
! Arguments
type(icebergs_gridded), pointer :: grd
real, dimension(:,:,:), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, k, halo, icount, io, jo
real :: mean, rms, SD, minv, maxv
real, dimension(lbound(fld,1):ubound(fld,1), lbound(fld,2):ubound(fld,2), lbound(fld,3):ubound(fld,3)) :: tmp

  halo=grd%halo
  mean=0.
  rms=0.
  sd=0.
  icount=0
  i=lbound(fld,1)+halo
  j=lbound(fld,2)+halo
  k=lbound(fld,3)
  minv=fld(i,j,k)
  maxv=fld(i,j,k)
  tmp(:,:,:)=0.
  io=grd%isd-lbound(fld,1)
  jo=grd%jsd-lbound(fld,2)
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        icount=icount+1
        mean=mean+fld(i,j,k)
        rms=rms+fld(i,j,k)**2
        minv=min(minv,fld(i,j,k))
        maxv=max(maxv,fld(i,j,k))
        tmp(i,j,k)=fld(i,j,k)*float(i+io+2*(j+jo)+3*(k-1))
      enddo
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        sd=sd+(fld(i,j,k)-mean)**2
      enddo
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  j=mpp_chksum( tmp(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum3: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  j=mpp_chksum( tmp(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum3* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#endif

end subroutine grd_chksum3

! ##############################################################################

subroutine grd_chksum2(grd, fld, txt)
! Arguments
type(icebergs_gridded), pointer :: grd
real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, icount
real :: mean, rms, SD, minv, maxv

  grd%tmp(:,:)=0.

  mean=0.
  rms=0.
  sd=0.
  icount=0
  minv=fld(grd%isc,grd%jsc)
  maxv=fld(grd%isc,grd%jsc)
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      icount=icount+1
      mean=mean+fld(i,j)
      rms=rms+fld(i,j)**2
      minv=min(minv,fld(i,j))
      maxv=max(maxv,fld(i,j))
      grd%tmp(i,j)=fld(i,j)*float(i+2*j)
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      sd=sd+(fld(i,j)-mean)**2
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(grd%isc:grd%iec,grd%jsc:grd%jec) )
  j=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum2: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i8)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(grd%isd:grd%ied,grd%jsd:grd%jed) )
  j=mpp_chksum( grd%tmp(grd%isd:grd%ied,grd%jsd:grd%jed) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum2* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#endif

end subroutine grd_chksum2

! ##############################################################################

subroutine bergs_chksum(bergs, txt, ignore_halo_violation)
! Arguments
type(icebergs), pointer :: bergs
character(len=*), intent(in) :: txt
logical, optional :: ignore_halo_violation
! Local variables
integer :: i, nbergs, ichk1, ichk2, ichk3, ichk4, ichk5, iberg
real, allocatable :: fld(:,:), fld2(:,:)
integer, allocatable :: icnt(:,:)
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
logical :: check_halo
integer :: grdi, grdj

! For convenience
  grd=>bergs%grd

  nbergs=count_bergs(bergs)
  call mpp_max(nbergs)
  nbergs = max(nbergs, 1)
  allocate( fld( nbergs, 19 ) ) !Changed from 11 to 19 by Alon
  allocate( fld2( nbergs, 19 ) ) !Changed from 11 to 19 by Alon
  allocate( icnt( grd%isd:grd%ied, grd%jsd:grd%jed ) )
  fld(:,:)=0.
  fld2(:,:)=0.
  icnt(:,:)=0
  grd%tmp(:,:)=0.

  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    this=>bergs%list(grdi,grdj)%first
    i=0; ichk5=0
    do while(associated(this))
      i=i+1
      iberg=berg_chksum(this)
      fld(i,1) = this%lon
      fld(i,2) = this%lat
      fld(i,3) = this%uvel
      fld(i,4) = this%vvel
      fld(i,5) = this%mass
      fld(i,6) = this%thickness
      fld(i,7) = this%width
      fld(i,8) = this%length
      fld(i,9) = this%axn !added by Alon
      fld(i,10) = this%ayn !added by Alon
      fld(i,11) = this%bxn !added by Alon
      fld(i,12) = this%byn !added by Alon
      fld(i,13) = this%uvel_old !added by Alon
      fld(i,14) = this%vvel_old !added by Alon
      fld(i,15) = this%lon_old !added by Alon
      fld(i,16) = this%lat_old !added by Alon
      fld(i,17) = time_hash(this) !Changed from 9 to 17 by Alon
      fld(i,18) = pos_hash(this) !Changed from 10 to 18 by Alon
      fld(i,19) = float(iberg) !Changed from 11 to 19 by Alon
      icnt(this%ine,this%jne)=icnt(this%ine,this%jne)+1
      fld2(i,:) = fld(i,:)*float( icnt(this%ine,this%jne) ) !*float( i )
      grd%tmp(this%ine,this%jne)=grd%tmp(this%ine,this%jne)+time_hash(this)*pos_hash(this)+log(this%mass)
      ichk5=ichk5+iberg
      this=>this%next
    enddo
  enddo ; enddo

  ichk1=mpp_chksum( fld )
  ichk2=mpp_chksum( fld2 )
  ichk3=mpp_chksum( grd%tmp )
  ichk4=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( ichk5 )
  nbergs=count_bergs(bergs)

  if (nbergs.ne.sum(icnt(:,:))) then
    write(*,'("diamonds, bergs_chksum: ",2(a,i8))') &
      '# bergs =', nbergs, ' sum(icnt) =',sum(icnt(:,:))
    call error_mesg('diamonds, bergs_chksum:', 'mismatch in berg count!', FATAL)
  endif

  check_halo=.true.
  if (present(ignore_halo_violation)) then
    if (ignore_halo_violation) check_halo=.false.
  endif
  if (check_halo.and.nbergs.ne.sum(icnt(grd%isc:grd%iec, grd%jsc:grd%jec))) then
    write(*,'("diamonds, bergs_chksum: ",2(a,i8))') &
      '# bergs =', nbergs, ' sum(icnt(comp_dom)) =',sum(icnt(:,:))
    call error_mesg('diamonds, bergs_chksum:', 'mismatch in berg count on computational domain!', FATAL)
  endif

  call mpp_sum(nbergs)
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, bergs_chksum: ",a18,6(x,a,"=",i22))') &
      txt, 'chksum', ichk1, 'chksum2', ichk2, 'chksum3', ichk3, 'chksum4', ichk4, 'chksum5', ichk5, '#', nbergs

  grd%tmp(:,:)=real(icnt(:,:))
  call grd_chksum2(grd,grd%tmp,'# of bergs/cell')

  deallocate( fld )
  deallocate( fld2 )
  deallocate( icnt )

  if (debug) call count_out_of_order(bergs,txt)

end subroutine bergs_chksum

! ##############################################################################

integer function list_chksum(first)
! Arguments
type(iceberg), pointer :: first
! Local variables
integer :: i
type(iceberg), pointer :: this

  this=>first
  i=0; list_chksum=0
  do while(associated(this))
    i=i+1
    list_chksum=list_chksum+berg_chksum(this)*i
    this=>this%next
  enddo

end function list_chksum

! ##############################################################################

integer function berg_chksum(berg)
! Arguments
type(iceberg), pointer :: berg
! Local variables
real :: rtmp(38) !Changed from 28 to 34 by Alon
integer :: itmp(38+4), i8=0, ichk1, ichk2, ichk3 !Changed from 28 to 34 by Alon
integer :: i

  rtmp(:)=0.
  rtmp(1)=berg%lon
  rtmp(2)=berg%lat
  rtmp(3)=berg%uvel
  rtmp(4)=berg%vvel
  rtmp(5)=berg%mass
  rtmp(6)=berg%thickness
  rtmp(7)=berg%width
  rtmp(8)=berg%length
  rtmp(9)=berg%start_lon
  rtmp(10)=berg%start_lat
  rtmp(11)=berg%start_day
  rtmp(12)=berg%start_mass
  rtmp(13)=berg%mass_scaling
  rtmp(14)=berg%mass_of_bits
  rtmp(15)=berg%heat_density
  rtmp(16)=berg%xi
  rtmp(17)=berg%yj
  rtmp(19)=berg%uo
  rtmp(20)=berg%vo
  rtmp(21)=berg%ui
  rtmp(22)=berg%vi
  rtmp(23)=berg%ua
  rtmp(24)=berg%va
  rtmp(25)=berg%ssh_x
  rtmp(26)=berg%ssh_y
  rtmp(27)=berg%cn
  rtmp(28)=berg%hi
  rtmp(29)=berg%axn 
  rtmp(30)=berg%ayn 
  rtmp(31)=berg%bxn 
  rtmp(32)=berg%byn 
  rtmp(33)=berg%uvel_old 
  rtmp(34)=berg%vvel_old 
  rtmp(35)=berg%lat_old 
  rtmp(36)=berg%lon_old 
  itmp(37)=berg%halo_berg 
  itmp(38)=berg%static_berg 
  itmp(1:38)=transfer(rtmp,i8) 
  itmp(39)=berg%start_year 
  itmp(40)=berg%ine 
  itmp(41)=berg%jne 
  itmp(42)=berg%iceberg_num 

  ichk1=0; ichk2=0; ichk3=0
  do i=1,38+4 
   ichk1=ichk1+itmp(i)
   ichk2=ichk2+itmp(i)*i
   ichk3=ichk3+itmp(i)*i*i
  enddo
  berg_chksum=ichk1+ichk2+ichk3

end function berg_chksum

! ##############################################################################

real function bilin(grd, fld, i, j, xi, yj)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
! Local variables

  if (old_bug_bilin) then
    bilin=(fld(i,j  )*(1.-xi)+fld(i-1,j  )*xi)*(1.-yj) &
         +(fld(i,j-1)*(1.-xi)+fld(i-1,j-1)*xi)*yj
  else
    bilin=(fld(i,j  )*xi+fld(i-1,j  )*(1.-xi))*yj &
         +(fld(i,j-1)*xi+fld(i-1,j-1)*(1.-xi))*(1.-yj)
  endif
end function bilin

! ##############################################################################

subroutine print_fld(grd, fld, label)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed)
character(len=*) :: label
! Local variables
integer :: i, j
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  write(stderrunit,'("pe=",i3,x,a8,32i10)') mpp_pe(),label,(i,i=grd%isd,grd%ied)
  do j=grd%jed,grd%jsd,-1
    write(stderrunit,'("pe=",i3,x,i8,32es10.2)') mpp_pe(),j,(fld(i,j),i=grd%isd,grd%ied)
  enddo

end subroutine print_fld

! ##############################################################################

logical function unitTests(bergs)
  type(icebergs), pointer :: bergs
  type(icebergs_gridded), pointer :: grd
  ! Local variables
  integer :: stderrunit,i,j

  ! This function returns True is a unit test fails
  unitTests=.false.
  ! For convenience
  grd=>bergs%grd
  stderrunit=stderr()
  
  i=grd%isc; j=grd%jsc
  call localTest( bilin(grd, grd%lon, i, j, 0., 1.), grd%lon(i-1,j) )
  call localTest( bilin(grd, grd%lon, i, j, 1., 1.), grd%lon(i,j) )
  call localTest( bilin(grd, grd%lat, i, j, 1., 0.), grd%lat(i,j-1) )
  call localTest( bilin(grd, grd%lat, i, j, 1., 1.), grd%lat(i,j) )

  contains
  subroutine localTest(answer, rightAnswer)
  real, intent(in) :: answer, rightAnswer
  if (answer==rightAnswer) return
  unitTests=.true.
  write(stderrunit,*) 'a=',answer,'b=',rightAnswer
  end subroutine localTest
end function unitTests

! ##############################################################################

!> Check for duplicates of icebergs on and across processors and issue an error
!! if any are detected
subroutine check_for_duplicates_in_parallel(bergs)
  type(icebergs), pointer :: bergs !< Icebergs
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: this
  integer :: stderrunit, i, j, k, nbergs, nbergs_total
  integer, dimension(:), allocatable :: ids ! List of ids of all bergs on this processor

  stderrunit=stderr()
  grd=>bergs%grd
  nbergs = count_bergs(bergs)
  nbergs_total = nbergs
  call mpp_sum(nbergs_total) ! Total number of bergs
  if (nbergs_total==0) return ! Skip the rest of the test

  k = 0
  if (nbergs>0) then
    allocate(ids(nbergs))
    do j = grd%jsc,grd%jec ; do i = grd%isc,grd%iec
      this=>bergs%list(i,j)%first
      do while (associated(this))
        k = k + 1
        ids(k) = this%iceberg_num
        this=>this%next
      enddo
    enddo ; enddo
  endif
  if (k /= nbergs) then
    write(stderrunit,*) 'counted bergs=',k,'count_bergs()=',nbergs
    call error_mesg('diamonds, check_for_duplicates:', 'Mismatch between concatenation of lists and count_bergs()!', FATAL)
  endif
  k = check_for_duplicate_ids_in_list(nbergs, ids, verbose=.true.)
  if (k /= 0) call error_mesg('diamonds, check_for_duplicates:', 'Duplicate berg detected across PEs!', FATAL)
  if (nbergs>0) deallocate(ids)
end subroutine check_for_duplicates_in_parallel

!> Returns error count of duplicates of integer values in a distributed list
integer function check_for_duplicate_ids_in_list(nbergs, ids, verbose)
  integer,               intent(in)    :: nbergs !< Length of ids
  integer, dimension(:), intent(inout) :: ids !< List of ids
  logical,               intent(in)    :: verbose !< True if messages should be written
  integer :: stderrunit, i, j, k, l, nbergs_total, ii, lowest_id, nonexistent_id
  logical :: have_berg

  stderrunit=stderr()
  nbergs_total = nbergs
  call mpp_sum(nbergs_total) ! Total number of bergs

  ! Establish lowest id or 0 across all PEs
  lowest_id = 0
  if (nbergs>0) lowest_id = minval(ids)
  call mpp_min(lowest_id)
  i = lowest_id
  nonexistent_id = lowest_id - 1
  if (nonexistent_id >= lowest_id) then
    write(stderrunit,*) 'Underflow in iceberg ids!',nonexistent_id,lowest_id,mpp_pe()
  endif
  ! Sort the list "ids" (largest first)
  do j = 1, nbergs-1
    do i = j+1, nbergs
      if (ids(j) < ids(i)) then
        ! Swap
        k = ids(i)
        ids(i) = ids(j)
        ids(j) = k
      endif
    enddo
  enddo
  ! Check for duplicates on processor
  check_for_duplicate_ids_in_list = 0
  do k = 1, nbergs-1
    if (ids(k) == ids(k+1)) then
      if (verbose) write(stderrunit,*) 'Duplicated berg on PE with id=',ids(k),'pe=',mpp_pe()
      check_for_duplicate_ids_in_list = check_for_duplicate_ids_in_list + 1
    endif
  enddo
  ! Check for duplicates across processor
  j = 1 ! Pointer to first berg in my list
  do k = 1, nbergs_total
    ! Set i to first id in my list
    if (j <= nbergs) then
      i = ids(j)
      have_berg = .true.
    else
      i = nonexistent_id
      have_berg = .false.
    endif
    l = i
    call mpp_max(l)
    if (have_berg .and. i == l) then
      ii = 1 ! This berg is mine
      j = j + 1
    else
      ii = 0 ! This berg is not mine
    endif
    call mpp_sum(ii)
    if (ii > 1) then
      if (verbose) write(stderrunit,*) 'Duplicated berg across PEs with id=',i,l,' seen',ii,' times pe=',mpp_pe(),k,j,nbergs
      check_for_duplicate_ids_in_list = check_for_duplicate_ids_in_list + 1
    elseif (ii == 0) then
      if (verbose) write(stderrunit,*) 'Berg not accounted for on all PEs with id=',i,l,' seen',ii,' times pe=',mpp_pe(),k,j,nbergs
    endif
  enddo

end function check_for_duplicate_ids_in_list

subroutine test_check_for_duplicate_ids_in_list()
  integer :: k
  integer, dimension(:), allocatable :: ids
  integer :: error_count

  allocate(ids(5))
  do k = 1,5
    ids(k) = k + 5*mpp_pe()
  enddo
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count /= 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('diamonds, test_check_for_duplicate_ids_in_list:', 'Unit test for clean list failed!', FATAL)
  endif
  if (mpp_pe() == mpp_root_pe()) ids(5) = ids(4)
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count == 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('diamonds, test_check_for_duplicate_ids_in_list:', 'Unit test for dirty list failed!', FATAL)
  endif
  if (mpp_pe() == mpp_root_pe()) ids(5) = 7 + 5*mpp_pe()
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count == 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('diamonds, test_check_for_duplicate_ids_in_list:', 'Unit test for a really dirty list failed!', FATAL)
  endif
  deallocate(ids)
end subroutine test_check_for_duplicate_ids_in_list

end module
