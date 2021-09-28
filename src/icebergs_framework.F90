!> Provides utilites for managing bergs in linked lists, and bonds between bergs
module ice_bergs_framework

! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.

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

integer :: buffer_width=36 ! This should be a parameter
integer :: buffer_width_traj=32 ! This should be a parameter
integer :: buffer_width_bond_traj=11 !This should be a parameter
integer, parameter :: nclasses=10 ! Number of ice bergs classes

!Local Vars
! Global data (minimal for debugging)
logical :: folded_north_on_pe = .false. !< If true, indicates the presence of the tri-polar grid
logical :: verbose=.false. !< Be verbose to stderr
logical :: budget=.true. !< Calculate budgets
logical :: debug=.false. !< Turn on debugging
logical :: really_debug=.false. !< Turn on debugging
logical :: parallel_reprod=.true. !< Reproduce across different PE decompositions
logical :: use_slow_find=.true. !< Use really slow (but robust) find_cell for reading restarts
logical :: ignore_ij_restart=.false. !< Read i,j location from restart if available (needed to use restarts on different grids)
logical :: generate_test_icebergs=.false. !< Create icebergs in absence of a restart file
logical :: use_roundoff_fix=.true. !< Use a "fix" for the round-off discrepancy between is_point_in_cell() and pos_within_cell()
logical :: old_bug_rotated_weights=.false. !< Skip the rotation of off-center weights for rotated halo updates
logical :: make_calving_reproduce=.false. !< Make the calving.res.nc file reproduce across pe count changes.
logical :: old_bug_bilin=.true. !< If true, uses the inverted bilinear function (use False to get correct answer)
character(len=10) :: restart_input_dir = 'INPUT/' !< Directory to look for restart files
integer, parameter :: delta_buf=25 !< Size by which to increment buffers
real, parameter :: pi_180=pi/180. !< Converts degrees to radians
real, parameter :: Rearth=6360000. !< Radius of earth (m)
logical :: fix_restart_dates=.true. !< After a restart, check that bergs were created before the current model date
logical :: do_unit_tests=.false. !< Conduct some unit tests
logical :: force_all_pes_traj=.false. !< Force all pes write trajectory files regardless of io_layout
logical :: mts=.false. !< Use multiple time stepping scheme
logical :: new_mts=.false. !If T, implicit accel is added 50% at the current time step and 50% at the next time step
logical :: save_bond_traj=.false. !<Save trajectory files for bonds
logical :: ewsame=.false. !<(F) set T if periodic and 2 PEs along the x direction (zonal) (i.e. E/W PEs are the same)
logical :: monitor_energy=.false. !<monitors energies: elastic (spring+collision), external, dissipated, fracture
logical :: iceberg_bonds_on=.False. ! True=Allow icebergs to have bonds, False=don't allow.
logical :: dem=.false. !< If T, run in DEM-mode with angular terms, variable stiffness, etc
character(len=11) :: fracture_criterion='none' !<'energy','stress','strain_rate','strain',or 'none'
logical :: use_damage=.false. !< Damage on bonds. Can evolve and serve as fracture criterion, or just to represent weaker bond

!Public params !Niki: write a subroutine to expose these
public nclasses,buffer_width,buffer_width_traj,buffer_width_bond_traj
public verbose, really_debug, debug, restart_input_dir,make_calving_reproduce,old_bug_bilin,use_roundoff_fix
public ignore_ij_restart, use_slow_find,generate_test_icebergs,old_bug_rotated_weights,budget
public orig_read, force_all_pes_traj
public mts,new_mts,save_bond_traj,ewsame,monitor_energy,iceberg_bonds_on
public dem, fracture_criterion

!Public types
public icebergs_gridded, xyt, iceberg, icebergs, buffer, bond, bond_xyt

!Public subs
public ice_bergs_framework_init
public send_bergs_to_other_pes
public update_halo_icebergs
public pack_traj_into_buffer2, unpack_traj_from_buffer2
public increase_ibuffer
public add_new_berg_to_list, count_out_of_order, check_for_duplicates
public insert_berg_into_list, create_iceberg, delete_iceberg_from_list, destroy_iceberg
public print_fld,print_berg, print_bergs,record_posn, push_posn, append_posn, check_position
public move_trajectory, move_all_trajectories
public form_a_bond, connect_all_bonds, show_all_bonds, bond_address_update
public find_cell, find_cell_by_search, find_cell_wide, count_bergs, is_point_in_cell, pos_within_cell, count_bonds
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
public split_id, id_from_2_ints, generate_id, cij_from_old_id, convert_old_id
public update_latlon,set_conglom_ids,transfer_mts_bergs,quad_interp_from_agrid
public pack_bond_traj_into_buffer2, unpack_bond_traj_from_buffer2, push_bond_posn, append_bond_posn
public update_and_break_bonds,break_bonds_dem,assign_n_bonds,reset_bond_rotation,update_bond_angles
public fracture_testing_initialization, orig_bond_length
public energy_tests_init,dem_tests_init
public init_dem_params
public set_constant_interaction_length_and_width,use_damage,damage_test_1_init

!> Container for gridded fields
type :: icebergs_gridded
  type(domain2D), pointer :: domain !< MPP parallel domain
  integer :: halo !< Nominal halo width
  integer :: isc !< Start i-index of computational domain
  integer :: iec !< End i-index of computational domain
  integer :: jsc !< Start j-index of computational domain
  integer :: jec !< End j-index of computational domain
  integer :: isd !< Start i-index of data domain
  integer :: ied !< End i-index of data domain
  integer :: jsd !< Start j-index of data domain
  integer :: jed !< End j-index of data domain
  integer :: isg !< Start i-index of global domain
  integer :: ieg !< End i-index of global domain
  integer :: jsg !< Start j-index of global domain
  integer :: jeg !< End j-index of global domain
  integer :: my_pe !< MPI PE index
  integer :: pe_N !< MPI PE index of PE to the north
  integer :: pe_S !< MPI PE index of PE to the south
  integer :: pe_E !< MPI PE index of PE to the east
  integer :: pe_W !< MPI PE index of PE to the west
  logical :: grid_is_latlon !< Flag to say whether the coordinate is in lat-lon degrees, or meters
  logical :: grid_is_regular !< Flag to say whether point in cell can be found assuming regular Cartesian grid
  real :: Lx !< Length of the domain in x direction
  real :: maxlon_c !< maximum x coord of computational domain
  real :: minlon_c !< minimum x coord of computational domain
  real, dimension(:,:), pointer :: lon=>null() !< Longitude of cell corners (degree E)
  real, dimension(:,:), pointer :: lat=>null() !< Latitude of cell corners (degree N)
  real, dimension(:,:), pointer :: lonc=>null() !< Longitude of cell centers (degree E)
  real, dimension(:,:), pointer :: latc=>null() !< Latitude of cell centers (degree N)
  real, dimension(:,:), pointer :: dx=>null() !< Length of cell edge (m)
  real, dimension(:,:), pointer :: dy=>null() !< Length of cell edge (m)
  real, dimension(:,:), pointer :: area=>null() !< Area of cell (m^2)
  real, dimension(:,:), pointer :: msk=>null() !< Ocean-land mask (1=ocean)
  real, dimension(:,:), pointer :: cos=>null() !< Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: sin=>null() !< Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: ocean_depth=>NULL() !< Depth of ocean (m)
  real, dimension(:,:), pointer :: uo=>null() !< Ocean zonal flow (m/s)
  real, dimension(:,:), pointer :: vo=>null() !< Ocean meridional flow (m/s)
  real, dimension(:,:), pointer :: ui=>null() !< Ice zonal flow (m/s)
  real, dimension(:,:), pointer :: vi=>null() !< Ice meridional flow (m/s)
  real, dimension(:,:), pointer :: ua=>null() !< Atmosphere zonal flow (m/s)
  real, dimension(:,:), pointer :: va=>null() !< Atmosphere meridional flow (m/s)
  real, dimension(:,:), pointer :: ssh=>null() !< Sea surface height (m)
  real, dimension(:,:), pointer :: sst=>null() !< Sea surface temperature (oC)
  real, dimension(:,:), pointer :: sss=>null() !< Sea surface salinity (psu)
  real, dimension(:,:), pointer :: cn=>null() !< Sea-ice concentration (0 to 1)
  real, dimension(:,:), pointer :: hi=>null() !< Sea-ice thickness (m)
  real, dimension(:,:,:), pointer :: melt_by_class=>null() !< Total icebergs melt rate by mass class (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_buoy_fl=>null() !< Footloose bergs buoyancy component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_eros_fl=>null() !< Footloose bergs erosion component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_conv_fl=>null() !< Footloose bergs convective component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: fl_parent_melt=>null() !< Footloose parent bergs melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: fl_child_melt=>null() !< Footloose child bergs melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: calving=>null() !< Calving mass rate [frozen runoff] (kg/s) (into stored ice)
  real, dimension(:,:), pointer :: calving_hflx=>null() !< Calving heat flux [heat content of calving] (W/m2) (into stored ice)
  real, dimension(:,:), pointer :: floating_melt=>null() !< Net melting rate to icebergs + bits (kg/s/m^2)
  real, dimension(:,:), pointer :: berg_melt=>null() !< Melting+erosion rate of parent icebergs (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_buoy=>null() !< Buoyancy component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_eros=>null() !< Erosion component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_conv=>null() !< Convective component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_src=>null() !< Mass flux from berg erosion into bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_melt=>null() !< Melting rate of bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_mass=>null() !< Mass distribution of bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: fl_bits_src=>null() !< Mass flux from berg into footloose bits (kg/s/m^2)
  real, dimension(:,:), pointer :: fl_bits_melt=>null() !< Melting rate of footloose bits (kg/s/m^2)
  real, dimension(:,:), pointer :: fl_bits_mass=>null() !< Mass distribution of footloose bits (kg/s/m^2)
  real, dimension(:,:), pointer :: fl_bergy_bits_mass=>null() !< Mass distribution of footloose bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: spread_mass=>null() !< Mass of icebergs after spreading (kg/m^2)
  real, dimension(:,:), pointer :: spread_mass_old=>null() !< Mass of icebergs after spreading old (kg/m^2)
  real, dimension(:,:), pointer :: spread_area=>null() !< Area of icebergs after spreading (m^2/m^2)
  real, dimension(:,:), pointer :: u_iceberg=>null() !< Average iceberg velocity in grid cell (mass weighted - but not spread mass weighted)
  real, dimension(:,:), pointer :: v_iceberg=>null() !< Average iceberg velocity in grid cell (mass weighted - but not spread mass weighted)
  real, dimension(:,:), pointer :: spread_uvel=>null() !< Average iceberg velocity in grid cell (spread area weighted)
  real, dimension(:,:), pointer :: spread_vvel=>null() !< Average iceberg velocity in grid cell (spread area weighted)
  real, dimension(:,:), pointer :: ustar_iceberg=>null() !< Frictional velocity below icebergs to be passed to ocean
  real, dimension(:,:), pointer :: virtual_area=>null() !< Virtual surface coverage by icebergs (m^2)
  real, dimension(:,:), pointer :: mass=>null() !< Mass distribution (kg/m^2)
  real, dimension(:,:,:), pointer :: mass_on_ocean=>null() !< Mass distribution partitioned by neighbor (kg)
  real, dimension(:,:,:), pointer :: area_on_ocean=>null() !< Area distribution partitioned by neighbor (m^2)
  real, dimension(:,:,:), pointer :: Uvel_on_ocean=>null() !< zonal velocity distribution partitioned by neighbor (m^2* m/s)
  real, dimension(:,:,:), pointer :: Vvel_on_ocean=>null() !< meridional momentum distribution partitioned by neighbor (m^2 m/s)
  real, dimension(:,:), pointer :: tmp=>null() !< Temporary work space
  real, dimension(:,:), pointer :: tmpc=>null() !< Temporary work space
  real, dimension(:,:,:), pointer :: stored_ice=>null() !< Accumulated ice mass flux at calving locations (kg)
  real, dimension(:,:), pointer :: rmean_calving=>null() !< Running mean for ice calving
  real, dimension(:,:), pointer :: rmean_calving_hflx=>null() !< Running mean for ice calving
  real, dimension(:,:), pointer :: stored_heat=>null() !< Heat content of stored ice (J)
  real, dimension(:,:,:), pointer :: real_calving=>null() !< Calving rate into iceberg class at calving locations (kg/s)
  real, dimension(:,:), pointer :: iceberg_heat_content=>null() !< Distributed heat content of bergs (J/m^2)
  real, dimension(:,:), pointer :: parity_x=>null() !< X component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  real, dimension(:,:), pointer :: parity_y=>null() !< Y component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  integer, dimension(:,:), pointer :: iceberg_counter_grd=>null() !< Counts icebergs created for naming purposes
  logical :: rmean_calving_initialized = .false. !< True if rmean_calving(:,:) has been filled with meaningful data
  logical :: rmean_calving_hflx_initialized = .false. !< True if rmean_calving_hflx(:,:) has been filled with meaningful data
  real :: coastal_drift=0. ! A velocity added to ocean currents to cause bergs to drift away from land cells
  real :: tidal_drift=0. ! Amplitude of a stochastic tidal velocity added to ocean currents to cause bergs to drift randomly
  !>@{
  !! Diagnostic handle
  integer :: id_uo=-1, id_vo=-1, id_calving=-1, id_stored_ice=-1, id_accum=-1, id_unused=-1, id_floating_melt=-1
  integer :: id_melt_buoy=-1, id_melt_eros=-1, id_melt_conv=-1, id_virtual_area=-1, id_real_calving=-1
  integer :: id_calving_hflx_in=-1, id_stored_heat=-1, id_melt_hflx=-1, id_heat_content=-1
  integer :: id_mass=-1, id_ui=-1, id_vi=-1, id_ua=-1, id_va=-1, id_sst=-1, id_cn=-1, id_hi=-1
  integer :: id_bergy_src=-1, id_bergy_melt=-1, id_bergy_mass=-1, id_berg_melt=-1
  integer :: id_fl_bits_src=-1, id_fl_bits_melt=-1, id_fl_bits_mass=-1, id_fl_bergy_bits_mass=-1
  integer :: id_rmean_calving=-1, id_rmean_calving_hflx=-1
  integer :: id_spread_mass=-1, id_spread_area=-1
  integer :: id_ssh=-1, id_fax=-1, id_fay=-1
  integer :: id_count=-1, id_chksum=-1, id_u_iceberg=-1, id_v_iceberg=-1, id_sss=-1, id_ustar_iceberg
  integer :: id_spread_uvel=-1, id_spread_vvel=-1
  integer :: id_melt_m_per_year=-1
  integer :: id_ocean_depth=-1
  integer :: id_melt_by_class=-1, id_melt_buoy_fl=-1, id_melt_eros_fl=-1, id_melt_conv_fl=-1
  integer :: id_fl_parent_melt=-1, id_fl_child_melt=-1
  !>@}

  real :: clipping_depth=0. !< The effective depth at which to clip the weight felt by the ocean [m].

end type icebergs_gridded

!> A link in the trajectory record (diagnostic)
type :: xyt
  real :: lon !< Longitude of berg (degree N or unit of grid coordinate)
  real :: lat !< Latitude of berg (degree N or unit of grid coordinate)
  real :: day !< Day of this record (days)
  real :: mass !< Mass of berg (kg)
  real :: thickness !< Thickness of berg (m)
  real :: width !< Width of berg (m)
  real :: length !< Length of berg (m)
  real :: uvel !< Zonal velocity of berg (m/s)
  real :: vvel !< Meridional velocity of berg (m/s)
  real :: axn
  real :: ayn
  real :: bxn
  real :: byn
  real :: uvel_prev
  real :: vvel_prev
  real :: fl_k
  real :: uo !< Zonal velocity of ocean (m/s)
  real :: vo !< Meridional velocity of ocean (m/s)
  real :: ui !< Zonal velocity of ice (m/s)
  real :: vi !< Meridional velocity of ice (m/s)
  real :: ua !< Zonal velocity of atmosphere (m/s)
  real :: va !< Meridional velocity of atmosphere (m/s)
  real :: ssh_x !< Zonal gradient of sea-surface height (nondim)
  real :: ssh_y !< Meridional gradient of sea-surface height (nondim)
  real :: sst !< Sea-surface temperature (Celsius)
  real :: sss !< Sea-surface salinity (1e-3)
  real :: cn !< Sea-ice concentration (nondim)
  real :: hi !< Sea-ice thickness (m)
  real :: halo_berg
  real :: static_berg
  real :: mass_of_bits !< Mass of bergy bits (kg)
  real :: mass_of_fl_bits !< Mass of footloose bits (kg)
  real :: mass_of_fl_bergy_bits !< Mass of bergy bits associated with the footloose bits (kg)
  real :: heat_density !< Heat density of berg (J/kg)
  real :: od !< Ocean depth
  real :: start_mass
  real :: mass_scaling
  integer :: year !< Year of this record (years)
  integer(kind=8) :: id = -1 !< Iceberg identifier
  type(xyt), pointer :: next=>null() !< Next link in list

    ! For MTS
  real, allocatable :: axn_fast
  real, allocatable :: ayn_fast
  real, allocatable :: bxn_fast
  real, allocatable :: byn_fast
  ! If iceberg_bonds_on
  integer, allocatable :: n_bonds
  ! If frac
  real, allocatable :: accum_bond_rotation
  ! For optional energy monitoring
  real, allocatable :: Ee !< Energy from bonds, elastic (J)
  real, allocatable :: Ed !< Energy from bonds, damping (J)
  real, allocatable :: Eext !< Energy from external forces (J)
  real, allocatable :: Ee_contact !< Energy from contact, elastic (J)
  real, allocatable :: Ed_contact !< Energy from contanct, damping (J)
  real, allocatable :: Efrac !< Energy from breaking bonds (J)
  real, allocatable :: Ee_temp !< Temporary storage: Energy from bonds, elastic (J)
  real, allocatable :: Ed_temp !< Temporary storage: Energy from bonds, damping (J)
  real, allocatable :: Eext_temp !< Temporary storage: Energy from external forces (J)
  real, allocatable :: Ee_contact_temp !< Temporary storage: Energy from contact, elastic (J)
  real, allocatable :: Ed_contact_temp !< Temporary storage: Energy from damping, elastic (J)
  ! For DEM-mode
  real, allocatable :: ang_vel !< Angular velocity
  real, allocatable :: ang_accel !< Angular acceleration
  real, allocatable :: rot !< Accumulated rotation
end type xyt

!> An iceberg object, used as a link in a linked list
type :: iceberg
  type(iceberg), pointer :: prev=>null() !< Previous link in list
  type(iceberg), pointer :: next=>null() !< Next link in list
  ! State variables (specific to the iceberg, needed for restarts)
  real :: lon !< Longitude of berg (degree N or unit of grid coordinate)
  real :: lat !< Latitude of berg (degree E or unit of grid coordinate)
  real :: lon_prev !< Previous timestep longitude of berg
  real :: lat_prev !< Previous timestep latitude of berg
  real :: uvel !< Zonal velocity of berg (m/s)
  real :: vvel !< Meridional velocity of berg (m/s)
  real :: mass !< Mass of berg (kg)
  real :: thickness !< Thickness of berg (m)
  real :: width !< Width of berg (m)
  real :: length !< Length of berg (m)
  real :: axn
  real :: ayn
  real :: bxn
  real :: byn
  real :: uvel_prev
  real :: vvel_prev
  real :: uvel_old
  real :: vvel_old
  real :: lon_old
  real :: lat_old
  real :: start_lon !< Longitude where berg was created (degree N or unit of grid coordinate)
  real :: start_lat !< Latitude where berg was created (degree E or unit of grid coordinate)
  real :: start_day !< Day that berg was created (days)
  real :: start_mass !< Mass berg had when created (kg)
  real :: mass_scaling !< Multiplier to scale mass when interpreting berg as a cloud of bergs (nondim)
  real :: mass_of_bits !< Mass of bergy bits following berg (kg)
  real :: mass_of_fl_bits !< Mass of footloose bergy bits following berg (kg)
  real :: mass_of_fl_bergy_bits !< Mass of bergy bits associated with the footloose bits (kg)
  real :: fl_k !< Cumulative number of footloose bergs to calve
  real :: heat_density !< Heat density of berg (J/kg)
  real :: halo_berg  ! Equal to zero for bergs on computational domain, and =1 for bergs on the halo
  real :: static_berg  ! Equal to 1 for icebergs which are static (not allowed to move). Might be extended to grounding later.
  integer :: start_year !< Year that berg was created (years)
  integer(kind=8) :: id !< Iceberg identifier
  integer :: ine !< Nearest i-index in NE direction (for convenience)
  integer :: jne !< Nearest j-index in NE direction (for convenience)
  real :: xi !< Non-dimensional x-coordinate within current cell (0..1)
  real :: yj !< Non-dimensional y-coordinate within current cell (0..1)
  ! Environment variables (as seen by the iceberg)
  real :: uo !< Zonal velocity of ocean (m/s)
  real :: vo !< Meridional velocity of ocean (m/s)
  real :: ui !< Zonal velocity of ice (m/s)
  real :: vi !< Meridional velocity of ice (m/s)
  real :: ua !< Zonal velocity of atmosphere (m/s)
  real :: va !< Meridional velocity of atmosphere (m/s)
  real :: ssh_x !< Zonal gradient of sea-surface height (nondim)
  real :: ssh_y !< Meridional gradient of sea-surface height (nondim)
  real :: sst !< Sea-surface temperature (Celsius)
  real :: sss !< Sea-surface salinity (1e-3)
  real :: cn !< Sea-ice concentration (nondim)
  real :: hi !< Sea-ice thickness (m)
  real :: od !< Ocean depth
  type(xyt), pointer :: trajectory=>null() !< Trajectory for this berg
  type(bond), pointer :: first_bond=>null() !< First element of bond list.

  ! For MTS
  real, allocatable :: axn_fast
  real, allocatable :: ayn_fast
  real, allocatable :: bxn_fast
  real, allocatable :: byn_fast
  integer, allocatable :: conglom_id !<ID for a conglomerate
  ! If iceberg_bonds_on
  integer, allocatable :: n_bonds
  ! If frac
  real, allocatable :: accum_bond_rotation
  ! For optional energy monitoring
  real, allocatable :: Ee !< Energy from bonds, elastic (J)
  real, allocatable :: Ed !< Energy from bonds, damping (J)
  real, allocatable :: Eext !< Energy from external forces (J)
  real, allocatable :: Ee_contact !< Energy from contact, elastic (J)
  real, allocatable :: Ed_contact !< Energy from contanct, damping (J)
  real, allocatable :: Efrac !< Energy from breaking bonds (J)
  real, allocatable :: Ee_temp !< Temporary storage: Energy from bonds, elastic (J)
  real, allocatable :: Ed_temp !< Temporary storage: Energy from bonds, damping (J)
  real, allocatable :: Eext_temp !< Temporary storage: Energy from external forces (J)
  real, allocatable :: Ee_contact_temp !< Temporary storage: Energy from contact, elastic (J)
  real, allocatable :: Ed_contact_temp !< Temporary storage: Energy from damping, elastic (J)
  ! For DEM-mode
  real, allocatable :: ang_vel !< Angular velocity
  real, allocatable :: ang_accel !< Angular acceleration
  real, allocatable :: rot !< Accumulated rotation
end type iceberg

!> A bond object connecting two bergs, used as a link in a linked list
type :: bond
  type(bond), pointer :: prev_bond=>null() !< Previous link in list
  type(bond), pointer :: next_bond=>null() !< Next link in list
  type(iceberg), pointer :: other_berg=>null()
  integer(kind=8) :: other_id !< ID of other berg
  integer :: other_berg_ine
  integer :: other_berg_jne
  real :: length
  type(bond_xyt), pointer :: bond_trajectory=>null()
  ! Fracture parameters:
  !The following are accumulated over all timesteps, unless fracture_criterion=='strain_rate',
  !in which case the following correspond to only the most recent timestep and units are s^(-1)
  real, allocatable :: damage
  real, allocatable :: rotation !radians
  real, allocatable :: rel_rotation !<rotation relative to the mean of all bonds for a berg
  real, allocatable :: n_frac_var !<normal strain, stress, or spring energy
  real, allocatable :: n_strain_rate !<for fracture_criterion=='strain_rate' only
  real, allocatable :: spring_pe
  ! For optional energy monitoring
  real, allocatable :: Ee !< Bond energy, elastic (J)
  real, allocatable :: Ed !< Bond energy, damping (J)
  real, allocatable :: axn_fast !< Zonal acceleration contribution to parent berg, explicit (m s^-2)
  real, allocatable :: ayn_fast !< Meridional  acceleration contribution to parent berg, explicit (m s^-2)
  real, allocatable :: bxn_fast !< Zonal acceleration contribution to parent berg, implicit (m s^-2)
  real, allocatable :: byn_fast !< Meridional acceleration contribution to parent berg, implicit (m s^-2)
  ! For DEM-mode
  real, allocatable :: tangd1 !< Accumulated tangential displacement, x-component (m)
  real, allocatable :: tangd2 !< Accumulated tangential displacement, y-component (m)
  real, allocatable :: nstress !< Normal stress (Pa)
  real, allocatable :: sstress !< Shear stress (Pa)
end type bond

!> A link in the bond trajectory record (diagnostic)
type :: bond_xyt
  real :: lon !< Longitude of bond (degree E or unit of grid coordinate)
  real :: lat !< Latitude of bond (degree N or unit of grid coordinate)
  integer :: year !< Year of this record (years)
  real :: day !< Day of this record (days)
  real :: length !<Length of bond
  real :: n1 !<Unit vector, x-component
  real :: n2 !<Unit vector, y-component
  integer(kind=8) :: id1 !<ID of first berg
  integer(kind=8) :: id2 !<ID of second berg
  type(bond_xyt), pointer :: next=>null() !< Next link in list
  ! Fracture parameters:
  real, allocatable :: damage
  real, allocatable :: rotation !radians
  real, allocatable :: rel_rotation !<rotation relative to the mean of all bonds for a berg
  real, allocatable :: n_frac_var !<normal strain, stress, or spring energy
  real, allocatable :: n_strain_rate !<for fracture_criterion=='strain_rate' only
  real, allocatable :: spring_pe
  ! For optional energy monitoring
  real, allocatable :: Ee !< Bond energy (elastic)
  real, allocatable :: Ed !< Bond energy (damping)
  real, allocatable :: axn_fast !< Zonal acceleration contribution to parent berg, explicit (m s^-2)
  real, allocatable :: ayn_fast !< Meridional  acceleration contribution to parent berg, explicit (m s^-2)
  real, allocatable :: bxn_fast !< Zonal acceleration contribution to parent berg, implicit (m s^-2)
  real, allocatable :: byn_fast !< Meridional acceleration contribution to parent berg, implicit (m s^-2)
  ! For DEM-mode
  real, allocatable :: tangd1 !< Accumulated tangential displacement, x-component (m)
  real, allocatable :: tangd2 !< Accumulated tangential displacement, y-component (m)
  real, allocatable :: nstress
  real, allocatable :: sstress
end type bond_xyt

! A dynamic buffer, used for communication, that packs types into rectangular memory
type :: buffer
  integer :: size=0 !< Size of buffer
  real, dimension(:,:), pointer :: data !< Buffer memory
end type buffer

!> A wrapper for the iceberg linked list (since an array of pointers is not allowed)
type :: linked_list
  type(iceberg), pointer :: first=>null() !< Pointer to the beginning of a linked list of bergs
end type linked_list

!> Container for all types and memory
type :: icebergs !; private !Niki: Ask Alistair why this is private. ice_bergs_io cannot compile if this is private!
  type(icebergs_gridded), pointer :: grd !< Container with all gridded data
  type(linked_list), dimension(:,:), allocatable :: list !< Linked list of icebergs
  type(xyt), pointer :: trajectories=>null() !< A linked list for detached segments of trajectories
  type(bond_xyt), pointer :: bond_trajectories=>null() !< A linked list for detached segments of bond trajectories
  real :: dt !< Time-step between iceberg calls
             !! \todo Should make dt adaptive?
  integer :: current_year !< Current year (years)
  real :: current_yearday !< Current year-day, 1.00-365.99, (days)
  real :: traj_area_thres !< Threshold for berg area (km^2) that must be exceeded to save a non-bonded berg trajectory
  real :: traj_area_thres_sntbc !< Threshold for berg area (km^2) for saving trajectory when using save_nonfl_traj_by_class
  real :: traj_area_thres_fl !< Threshold for berg area (km^2) that must be exceeded to save a non-bonded FL berg trajectory
  real :: traj_sample_hrs !< Period between sampling for trajectories (hours)
  real :: traj_write_hrs !< Period between writing of trajectories (hours)
  real :: verbose_hrs !< Period between terminal status reports (hours)
  integer :: max_bonds
  !>@{
  !! Handles for clocks
  integer :: clock, clock_mom, clock_the, clock_int, clock_cal, clock_com1, clock_fl1, clock_com2, clock_fl2
  integer :: clock_ini, clock_ior, clock_iow, clock_dia
  integer :: clock_trw, clock_trp
  integer :: clock_btrw, clock_btrp
  !>@}
  real :: rho_bergs !< Density of icebergs [kg/m^3]
  real :: spring_coef !< Spring constant for bonded iceberg interactions
  real :: contact_spring_coef !<Spring coefficient for berg collisions -Alex
  real :: cdrag_grounding=0.0 !< Drag coefficient against ocean bottom
  real :: h_to_init_grounding
  character(len=11) :: fracture_criterion='none' !<'damage','energy','stress','strain_rate','strain',or 'none'
  real :: frac_thres_n=0.0 !normal fracture strain threshold
  real :: frac_thres_t=0.0 !tangential fracture strain threshold
  logical :: damage_test_1=.false. !< Sets initial damage for bonds that overlap y=17000.0
  logical :: debug_write=.false. !< Sets traj_sample_hours=traj_write_hours & includes halo bergs in write
  real :: bond_coef !< Spring constant for iceberg bonds
  real :: radial_damping_coef !< Coefficient for relative iceberg motion damping (radial component) -Alon
  real :: tangental_damping_coef !< Coefficient for relative iceberg motion damping (tangential component) -Alon
  real :: LoW_ratio !< Initial ratio L/W for newly calved icebergs
  real :: bergy_bit_erosion_fraction !< Fraction of erosion melt flux to divert to bergy bits
  real :: sicn_shift !< Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
  real :: lat_ref=0. !< Reference latitude for f-plane (when this option is on)
  real :: u_override=0.0 !< Overrides the u velocity of icebergs (for ocean testing)
  real :: v_override=0.0 !< Overrides the v velocity of icebergs (for ocean testing)
  real :: utide_icebergs= 0. !< Tidal speeds, set to zero for now.
  real :: ustar_icebergs_bg=0.001 !< Background u_star under icebergs. This should be linked to a value felt by the ocean boundary layer
  real :: cdrag_icebergs =  1.5e-3 !< Momentum Drag coef, taken from HJ99 (Holland and Jenkins 1999)
  real :: initial_orientation=0. !< Iceberg orientation relative to this angle (in degrees). Used for hexagonal mass spreading.
  real :: Gamma_T_3EQ=0.022 !< Non-dimensional heat-transfer coefficient
  real :: melt_cutoff=-1.0 !< Minimum ocean thickness for melting to occur (is not applied for values < 0)
  logical :: const_gamma=.True. !< If true uses a constant heat transfer coefficient, from which the salt transfer is calculated
  real, dimension(:), pointer :: initial_mass_s, distribution_s, mass_scaling_s !< Southern hemisphere
  real, dimension(:), pointer :: initial_thickness_s, initial_width_s, initial_length_s !< Southern hemisphere
  logical :: separate_distrib_for_n_hemisphere=.False. ! Flag to use a separate berg distribution/mass/mass scaling/init thickness for N hemisphere
  real, dimension(:), pointer :: initial_mass_n, distribution_n, mass_scaling_n !< Northern hemisphere
  real, dimension(:), pointer :: initial_thickness_n, initial_width_n, initial_length_n !< Northern hemisphere
  logical :: restarted=.false. !< Indicate whether we read state from a restart or not
  logical :: use_operator_splitting=.true. !< Use first order operator splitting for thermodynamics
  logical :: add_weight_to_ocean=.true. !< Add weight of bergs to ocean
  logical :: passive_mode=.false. !< Add weight of icebergs + bits to ocean
  logical :: time_average_weight=.false. !< Time average the weight on the ocean
  logical :: Runge_not_Verlet=.True. !< True=Runge-Kutta, False=Verlet.
  logical :: use_mixed_melting=.False. !< If true, then the melt is determined partly using 3 eq model partly using iceberg parameterizations (according to iceberg bond number)
  logical :: internal_bergs_for_drag=.False. !< True=reduces side drag for bonded elements in momentum equation.
  logical :: apply_thickness_cutoff_to_gridded_melt=.False. !< Prevents melt for ocean thickness below melt_cuttoff (applied to gridded melt fields)
  logical :: apply_thickness_cutoff_to_bergs_melt=.False. !< Prevents melt for ocean thickness below melt_cuttoff (applied to bergs)
  logical :: use_updated_rolling_scheme=.false. !< True to use the aspect ratio based rolling scheme rather than incorrect version of WM scheme (set tip_parameter=1000. for correct WM scheme)
  logical :: pass_fields_to_ocean_model=.False. !< Iceberg area, mass and ustar fields are prepared to pass to ocean model
  logical :: use_mixed_layer_salinity_for_thermo=.False. !< If true, then model uses ocean salinity for 3 and 2 equation melt model.
  logical :: find_melt_using_spread_mass=.False. !< If true, then the model calculates ice loss by looping at the spread_mass before and after.
  logical :: Use_three_equation_model=.True. !< Uses 3 equation model for melt when ice shelf type thermodynamics are used.
  logical :: melt_icebergs_as_ice_shelf=.False. !< Uses iceshelf type thermodynamics
  logical :: Iceberg_melt_without_decay=.False. !< Allows icebergs meltwater fluxes to enter the ocean, without the iceberg decaying or changing shape.
  logical :: add_iceberg_thickness_to_SSH=.False. !< Adds the iceberg contribution to SSH.
  logical :: override_iceberg_velocities=.False. !< Allows you to set a fixed iceberg velocity for all non-static icebergs.
  logical :: use_f_plane=.False. !< Flag to use a f-plane for the rotation
  logical :: rotate_icebergs_for_mass_spreading=.True. !< Flag allows icebergs to rotate for spreading their mass (in hexagonal spreading mode)
  logical :: set_melt_rates_to_zero=.False. !< Sets all melt rates to zero, for testing purposes (thermodynamics routine is still run)
  logical :: hexagonal_icebergs=.False. !< True treats icebergs as rectangles, False as hexagonal elements (for the purpose of mass spreading)
  logical :: allow_bergs_to_roll=.True. !< Allows icebergs to roll over when rolling conditions are met
  logical :: ignore_missing_restart_bergs=.False. !< True Allows the model to ignore icebergs missing in the restart.
  logical :: require_restart=.false. !< If true, requires a restart file to be present when starting the model.
  logical :: Static_icebergs=.False. !< True= icebergs do no move
  logical :: only_interactive_forces=.False. !< Icebergs only feel interactive forces, and not ocean, wind...
  logical :: halo_debugging=.False. !< Use for debugging halos (remove when its working)
  logical :: save_short_traj=.True. !< True saves only lon,lat,time,id in iceberg_trajectory.nc
  logical :: save_fl_traj=.True. ! True saves short traj, plus masses and footloose parameters in iceberg_trajectory.nc
  real :: save_all_traj_year=huge(0.0) ! Year at which all trajectories (for all berg areas) are saved
  logical :: save_nonfl_traj_by_class=.false. ! Save non-footloose trajectories based on initial mass class
  real :: save_traj_by_class_start_mass_thres_n=0.0 ! The northern hemisphere mass thres for save_nonfl_traj_by_class
  real :: save_traj_by_class_start_mass_thres_s=0.0 ! The southern hemisphere mass thres for save_nonfl_traj_by_class
  logical :: ignore_traj=.False. !< If true, then model does not write trajectory data at all
  logical :: iceberg_bonds_on=.False. !< True=Allow icebergs to have bonds, False=don't allow.
  logical :: manually_initialize_bonds=.False. !< True= Bonds are initialize manually.
  logical :: use_new_predictive_corrective =.False. !< Flag to use Bob's predictive corrective iceberg scheme- Added by Alon
  logical :: interactive_icebergs_on=.false. !< Turn on/off interactions between icebergs  - Added by Alon
  logical :: scale_damping_by_pmag=.true. !< Scales damping by magnitude of (projection matrix dot relative velocity)
  logical :: critical_interaction_damping_on=.true. !< Sets the damping on relative iceberg velocity to critical value - Added by Alon
  logical :: tang_crit_int_damp_on=.true. !<crit interaction damping for tangential component?
  logical :: use_old_spreading=.true. !< If true, spreads iceberg mass as if the berg is one grid cell wide
  logical :: read_ocean_depth_from_file=.false. !< If true, ocean depth is read from a file.
  integer(kind=8) :: debug_iceberg_with_id = -1 !< If positive, monitors a berg with this id

  real :: length_for_manually_initialize_bonds=1000.0 !< If manually initialing bonds, only  bond if dist between particles is < this length -Added by Alex
  logical :: manually_initialize_bonds_from_radii=.false. !< If manually initialing bonds, form bonds if dist between particles is < 1.25x the smaller radii-Alex
  real :: speed_limit=0. !< CFL speed limit for a berg [m/s]
  real :: tau_calving=0. !< Time scale for smoothing out calving field (years)
  real :: tip_parameter=0. !< If non-zero, overrides iceberg rolling critical ratio (use zero to get parameter directly from ice and seawater densities)
  real :: grounding_fraction=0. !< Fraction of water column depth at which grounding occurs
  type(buffer), pointer :: obuffer_n=>null() !< Buffer for outgoing bergs to the north
  type(buffer), pointer :: ibuffer_n=>null() !< Buffer for incoming bergs from the north
  type(buffer), pointer :: obuffer_s=>null() !< Buffer for outgoing bergs to the south
  type(buffer), pointer :: ibuffer_s=>null() !< Buffer for incoming bergs from the south
  type(buffer), pointer :: obuffer_e=>null() !< Buffer for outgoing bergs to the east
  type(buffer), pointer :: ibuffer_e=>null() !< Buffer for incoming bergs from the east
  type(buffer), pointer :: obuffer_w=>null() !< Buffer for outgoing bergs to the west
  type(buffer), pointer :: ibuffer_w=>null() !< Buffer for incoming bergs from the west
  type(buffer), pointer :: obuffer_io=>null() !< Buffer for outgoing bergs during i/o
  type(buffer), pointer :: ibuffer_io=>null() !< Buffer for incoming bergs during i/o
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
  real :: fl_bits_mass_start=0., fl_bits_mass_end=0.
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
  real :: fl_bits_src=0., fl_bits_melt=0.
  integer :: nbergs_calved=0, nbergs_calved_fl=0, nbergs_melted=0, nbergs_start=0, nbergs_end=0
  integer :: nspeeding_tickets=0
  integer :: nbonds=0
  integer, dimension(:), pointer :: nbergs_calved_by_class_s=>null()
  integer, dimension(:), pointer :: nbergs_calved_by_class_n=>null()
  ! mts parameters - added by Alex
  logical :: mts=.false. !< Use multiple timestepping scheme (substep size is automatically determined)
  integer :: mts_sub_steps
  real :: mts_fast_dt
  integer :: mts_part !for turning on/off berg interactions/collisions during different mts scheme
  logical :: remove_unused_bergs=.true. !remove unneeded bergs after PEs transfers
  real :: contact_distance
  integer :: contact_cells_lon=1,contact_cells_lat=1 !how many cells to search to find contact pairs
  logical :: writeandstop=.false. !for debugging
  logical :: force_convergence !experimental MTS convergence scheme to better conserve momentum during collisions
  logical :: explicit_inner_mts=.false. !If T, inner MTS iterations are treated explicitly
  real :: convergence_tolerance=1.e-8 !tolerance for the MTS convergence scheme
  ! DEM-mode parameters
  logical :: dem=.false. !if T, run in DEM-mode with angular terms, variable stiffness, etc
  logical :: ignore_tangential_force=.false.
  real :: poisson=0.3 ! Poisson's ratio
  real :: dem_spring_coef=0.
  real :: dem_damping_coef=0.1
  logical :: dem_shear_for_frac_only=.false. !< If true, DEM shear is calculated for fracture, but zeroed for berg interactions
  integer :: dem_beam_test=0 !1=Simply supported beam,2=cantilever beam,3=angular vel tes
  logical :: uniaxial_test=.false. !adds a tensile stress to the east-most berg element
  real :: dem_tests_start_lon !starting lon of west-most berg element for uniaxial/dem tests
  real :: dem_tests_end_lon !starting lon of east-most berg element for uniaxial/dem tests
  ! Element interactions
  logical :: constant_interaction_LW=.false. !< Use a constant element length & width during berg interactions
  real :: constant_length=0.!< If constant_interaction_LW, the constant length. If 0, will be set to max initial length
  real :: constant_width=0. !< If constant_interaction_LW, the constant width.  If 0, will be set to max initial width
  logical :: use_spring_for_land_contact=.false. !< Treat contact with masked (land) cells like contact with a static berg
  ! Footloose calving parameters [England et al (2020) Modeling the breakup of tabular icebergs. Sci. Adv.]
  logical :: fl_use_poisson_distribution=.true. !< fl_r is (T) mean of Poisson distribution to determine k, or (F) k=fl_r
  logical :: fl_use_perimeter=.false. !< scale number of footloose bergs to calve by perimeter of the parent berg
  real :: fl_youngs=1.e8 !< Young's modulus for footloose calculations (Pa)
  real :: fl_strength=500. !< yield stress for footloose calculations (kPa)
  logical :: fl_use_l_scale=.false. !< footloose based on side melt and/or erosion
  real :: fl_l_scale=1. !< scaling factor for footloose based on side melt
  logical :: fl_l_scale_erosion_only=.true. !< fl according to erosion - buoyant convection
  real :: fl_r=0. !< footloose average number of bergs calved per fl_r_s
  real :: fl_r_s=0. !< seconds over which fl_r footloose bergs calve
  logical :: displace_fl_bergs=.true. !< footloose berg positions are randomly assigned along edges of parent berg
  character(len=11) :: fl_style='new_bergs' !< Evolve footloose bergs individually as 'new_bergs', or as a group with size 'fl_bits'
  logical :: fl_bits_erosion_to_bergy_bits=.true. !< Erosion from footloose bits becomes bergy bits
  real :: fl_k_scale_by_perimeter=0 !< If greater than 0, scales FL k by (berg perimeter (m))/(this scaling (m))
  real :: new_berg_from_fl_bits_mass_thres=huge(0.) ! Create a new berg from FL bits when mass_of_fl_bits exceeds this value
  real :: fl_bits_scale_l=0.9 !< For determining dimensions of FL bits berg; FL_bits length    = fl_bits_scale_l * 3*l_b
  real :: fl_bits_scale_w=0.9 !< For determining dimensions of FL bits berg; FL_bits width     = fl_bits_scale_w * 3*l_b
  real :: fl_bits_scale_t=0.9 !< For determining dimensions of FL bits berg; FL_bits thickness = fl_bits_scale_t * T
end type icebergs

!> Read original restarts. Needs to be module global so can be public to icebergs_mod.
!! \todo Remove when backward compatibility no longer needed
logical :: orig_read=.false.

!> Version of file provided by CPP macro (usually set to git hash)
#ifdef _FILE_VERSION
character(len=128) :: version = _FILE_VERSION !< Version of file
#else
character(len=128) :: version = 'unknown' !< Version of file
#endif

!> Set a value in the buffer at position (counter,n) after incrementing counter
interface push_buffer_value
  module procedure push_buffer_rvalue, push_buffer_ivalue
end interface

!> Get a value in the buffer at position (counter,n) after incrementing counter
interface pull_buffer_value
  module procedure pull_buffer_rvalue, pull_buffer_ivalue
end interface

contains

!> Initializes parallel framework
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
type(icebergs), pointer :: bergs !< Container for all types and memory
integer, intent(in) :: gni !< Number grid cells in i-direction
integer, intent(in) :: gnj !< Number grid cells in j-direction
integer, intent(in) :: layout(2) !< Number of processing cores in i,j direction
integer, intent(in) :: io_layout(2) !< Number of i/o cores in i,j direction
integer, intent(in) :: axes(2) !< Diagnostic axes
integer, intent(in) :: dom_x_flags !< Domain flags in i-direction
integer, intent(in) :: dom_y_flags !< Domain flags in j-direction
real, intent(in) :: dt !< Time-step (s)
type (time_type), intent(in) :: Time !< Current model time
real, dimension(:,:), intent(in) :: ice_lon !< Longitude of cell corners using NE convention (degree E)
real, dimension(:,:), intent(in) :: ice_lat !< Latitude of cell corners using NE convention (degree N)
real, dimension(:,:), intent(in) :: ice_wet !< Wet/dry mask (1 is wet, 0 is dry) of cell centers
real, dimension(:,:), intent(in) :: ice_dx !< Zonal length of cell on northern side (m)
real, dimension(:,:), intent(in) :: ice_dy !< Meridional length of cell on eastern side (m)
real, dimension(:,:), intent(in) :: ice_area !< Area of cells (m^2, or non-dim is fractional_area=True)
real, dimension(:,:), intent(in) :: cos_rot !< Cosine from rotation matrix to lat-lon coords
real, dimension(:,:), intent(in) :: sin_rot !< Sine from rotation matrix to lat-lon coords
real, dimension(:,:), intent(in),optional :: ocean_depth !< Depth of ocean bottom (m)
logical, intent(in), optional :: maskmap(:,:) !< Masks out parallel cores
logical, intent(in), optional :: fractional_area !< If true, ice_area contains cell area as fraction of entire spherical surface

! Namelist parameters (and defaults)
integer :: halo=4 ! Width of halo region
real :: traj_area_thres=0. ! Threshold for berg area (km^2) that must be exceeded to save a non-bonded berg trajectory
real :: traj_area_thres_sntbc=0. !< Threshold for berg area (km^2) for saving trajectory when using save_nonfl_traj_by_class
real :: traj_area_thres_fl=huge(0.0) ! Threshold for berg area (km^2) that must be exceeded to save a non-bonded FL berg trajectory
real :: traj_sample_hrs=24. ! Period between sampling of position for trajectory storage
real :: traj_write_hrs=480. ! Period between writing sampled trajectories to disk
real :: verbose_hrs=24. ! Period between verbose messages
integer :: max_bonds=6 ! Maximum number of iceberg bond passed between processors
real :: rho_bergs=850. ! Density of icebergs
real :: spring_coef=1.e-8 ! Spring constant for iceberg interactions (this seems to be the highest stable value)
real :: contact_spring_coef=0. !Spring coef for berg collisions (is set to spring_coef if not specified)
logical :: uniaxial_test=.false. !adds a tensile stress to the east-most berg element
real :: cdrag_grounding=0.0 ! Drag coefficient against ocean bottom
real :: h_to_init_grounding=100.0
!character(len=11) :: fracture_criterion ! 'energy','stress','strain_rate','strain',or 'none'
real :: frac_thres_n=0.0 !normal fracture strain threshold
real :: frac_thres_t=0.0 !tangential fracture strain threshold
logical :: damage_test_1=.false. ! Sets initial damage for bonds that overlap y=17000.0
real :: frac_thres_scaling=1.0 !scaling factor for frac_thres_n and frac_thres_t, useful for tuning
logical :: debug_write=.false. !Sets traj_sample_hours=traj_write_hours & includes halo bergs in write
real :: bond_coef=1.e-8 ! Spring constant for iceberg bonds - not being used right now
real :: radial_damping_coef=1.e-4 ! Coefficient for relative iceberg motion damping (radial component) -Alon
real :: tangental_damping_coef=2.e-5 ! Coefficient for relative iceberg motion damping (tangential component) -Alon
real :: LoW_ratio=1.5 ! Initial ratio L/W for newly calved icebergs
real :: bergy_bit_erosion_fraction=0. ! Fraction of erosion melt flux to divert to bergy bits
real :: sicn_shift=0. ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
real :: lat_ref=0. ! Reference latitude for f-plane (when this option is on)
real :: u_override=0.0 ! Overrides the u velocity of icebergs (for ocean testing)
real :: v_override=0.0 ! Overrides the v velocity of icebergs (for ocean testing)
real :: Lx=360. ! Length of domain in x direction, used for periodicity (use a huge number for non-periodic)
real :: initial_orientation=0. ! Iceberg orientation relative to this angle (in degrees). Used for hexagonal mass spreading.
real :: utide_icebergs= 0. ! Tidal speeds, set to zero for now.
real :: ustar_icebergs_bg=0.001 ! Background u_star under icebergs. This should be linked to a value felt by the ocean boundary layer
real :: cdrag_icebergs =  1.5e-3 ! Momentum Drag coef, taken from HJ99  (Holland and Jenkins 1999)
real :: Gamma_T_3EQ=0.022 ! Non-dimensional heat-transfer coefficient
real :: melt_cutoff=-1.0 ! Minimum ocean thickness for melting to occur (is not applied for values < 0)
logical :: const_gamma=.True. ! If true uses a constant heat transfer coefficient, from which the salt transfer is calculated
logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
logical :: add_weight_to_ocean=.true. ! Add weight of icebergs + bits to ocean
logical :: passive_mode=.false. ! Add weight of icebergs + bits to ocean
logical :: time_average_weight=.false. ! Time average the weight on the ocean
real :: length_for_manually_initialize_bonds=1000.0 ! if manually init bonds, only  bond if dist between particles is .lt. this length - Alex
logical :: manually_initialize_bonds_from_radii=.false. ! if manually init bonds, form bonds if dist between particles is < 1.25x the smaller radii - Alex
real :: speed_limit=0. ! CFL speed limit for a berg
real :: tau_calving=0. ! Time scale for smoothing out calving field (years)
real :: tip_parameter=0. ! Parameter to override iceberg rolling critical ratio (use zero to get parameter directly from ice and seawater densities
real :: grounding_fraction=0. ! Fraction of water column depth at which grounding occurs
real :: coastal_drift=0. ! A velocity added to ocean currents to cause bergs to drift away from land cells
real :: tidal_drift=0. ! Amplitude of a stochastic tidal velocity added to ocean currents to cause bergs to drift randomly
logical :: Runge_not_Verlet=.True. ! True=Runge Kutta, False=Verlet.
logical :: use_mixed_melting=.False. ! If true, then the melt is determined partly using 3 eq model partly using iceberg parameterizations (according to iceberg bond number)
logical :: internal_bergs_for_drag=.False. ! True=reduces side drag for bonded elements in momentum equation.
logical :: apply_thickness_cutoff_to_gridded_melt=.False. ! Prevents melt for ocean thickness below melt_cuttoff (applied to gridded melt fields)
logical :: apply_thickness_cutoff_to_bergs_melt=.False. ! Prevents melt for ocean thickness below melt_cuttoff (applied to bergs)
logical :: use_updated_rolling_scheme=.false. ! Use the corrected Rolling Scheme rather than the erroneous one
logical :: pass_fields_to_ocean_model=.False. ! Iceberg area, mass and ustar fields are prepared to pass to ocean model
logical :: use_mixed_layer_salinity_for_thermo=.False. ! If true, then model uses ocean salinity for 3 and 2 equation melt model.
logical :: find_melt_using_spread_mass=.False. ! If true, then the model calculates ice loss by looping at the spread_mass before and after.
logical :: Use_three_equation_model=.True. ! Uses 3 equation model for melt when ice shelf type thermodynamics are used.
logical :: melt_icebergs_as_ice_shelf=.False. ! Uses iceshelf type thermodynamics
logical :: Iceberg_melt_without_decay=.False. ! Allows icebergs meltwater fluxes to enter the ocean, without the iceberg decaying or changing shape.
logical :: add_iceberg_thickness_to_SSH=.False. ! Adds the iceberg contribution to SSH.
logical :: override_iceberg_velocities=.False. ! Allows you to set a fixed iceberg velocity for all non-static icebergs.
logical :: use_f_plane=.False. ! Flag to use a f-plane for the rotation
logical :: grid_is_latlon=.True. ! True means that the grid is specified in lat lon, and uses to radius of the earth to convert to distance
logical :: grid_is_regular=.True. ! Flag to say whether point in cell can be found assuming regular Cartesian grid
logical :: rotate_icebergs_for_mass_spreading=.True. ! Flag allows icebergs to rotate for spreading their mass (in hexagonal spreading mode)
logical :: set_melt_rates_to_zero=.False. ! Sets all melt rates to zero, for testing purposes (thermodynamics routine is still run)
logical :: allow_bergs_to_roll=.True. ! Allows icebergs to roll over when rolling conditions are met
logical :: hexagonal_icebergs=.False. ! True treats icebergs as rectangles, False as hexagonal elements (for the purpose of mass spreading)
logical :: ignore_missing_restart_bergs=.False. ! True Allows the model to ignore icebergs missing in the restart.
logical :: require_restart=.false. ! If true, requires a restart file to be present when starting the model.
logical :: Static_icebergs=.False. ! True= icebergs do no move
logical :: only_interactive_forces=.False. ! Icebergs only feel interactive forces, and not ocean, wind...
logical :: halo_debugging=.False. ! Use for debugging halos (remove when its working)
logical :: save_short_traj=.True. ! True saves only lon,lat,time,id in iceberg_trajectory.nc
logical :: save_fl_traj=.True. ! True saves short traj, plus masses and footloose parameters in iceberg_trajectory.nc
real :: save_all_traj_year=huge(0.0) ! Year at which all trajectories (for all berg areas) are saved
logical :: save_nonfl_traj_by_class=.false. ! Save non-footloose trajectories based on initial mass class
real :: save_traj_by_class_start_mass_thres_n=0.0 ! The northern hemisphere mass thres for save_nonfl_traj_by_class
real :: save_traj_by_class_start_mass_thres_s=0.0 ! The southern hemisphere mass thres for save_nonfl_traj_by_class
logical :: ignore_traj=.False. ! If true, then model does not traj trajectory data at all
!logical :: iceberg_bonds_on=.False. ! True=Allow icebergs to have bonds, False=don't allow.
logical :: manually_initialize_bonds=.False. ! True= Bonds are initialize manually.
logical :: use_new_predictive_corrective =.False. ! Flag to use Bob's predictive corrective iceberg scheme- Added by Alon
logical :: interactive_icebergs_on=.false. ! Turn on/off interactions between icebergs  - Added by Alon
logical :: scale_damping_by_pmag=.true. ! Scales damping by magnitude of (projection matrix dot relative velocity)
logical :: critical_interaction_damping_on=.true. ! Sets the damping on relative iceberg velocity to critical value - Added by Alon
logical :: tang_crit_int_damp_on=.true. ! Critical interaction damping for tangential component?
logical :: do_unit_tests=.false. ! Conduct some unit tests
logical :: input_freq_distribution=.true. ! Flag to show if input distribution is freq or mass dist (=1 if input is a freq dist, =0 to use an input mass dist)
logical :: read_old_restarts=.false. ! Legacy option that does nothing
logical :: use_old_spreading=.true. ! If true, spreads iceberg mass as if the berg is one grid cell wide
logical :: read_ocean_depth_from_file=.false. ! If true, ocean depth is read from a file.
integer :: mts_sub_steps=-1 ! If -1, the number of mts sub-steps will be automatically determined
logical :: remove_unused_bergs=.true. ! Remove unneeded bergs after PEs transfers
real :: contact_distance=0.0 ! For unbonded berg interactions, collision is assumed at max(contact_distance,sum of the 2 bergs radii)
logical :: force_convergence=.false. ! Experimental MTS convergence scheme that better preserves momentum during collisions
logical :: explicit_inner_mts=.false. !If T, inner MTS iterations are treated explicitly
real :: convergence_tolerance=1.e-8 ! Tolerance for the MTS force_convergence scheme
! Initial mass, distribution, scaling, thickness at calving. Default is gladstone et al 2001 for Southern hemisphere:
real, dimension(nclasses) :: initial_mass=(/8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11/) ! Mass thresholds between iceberg classes (kg)
real, dimension(nclasses) :: distribution=(/0.24, 0.12, 0.15, 0.18, 0.12, 0.07, 0.03, 0.03, 0.03, 0.02/) ! Fraction of calving to apply to this class (non-dim) ,
real, dimension(nclasses) :: mass_scaling=(/2000, 200, 50, 20, 10, 5, 2, 1, 1, 1/) ! Ratio between effective and real iceberg mass (non-dim)
real, dimension(nclasses) :: initial_thickness=(/40., 67., 133., 175., 250., 250., 250., 250., 250., 250./) ! Total thickness of newly calved bergs (m)
!Default is Bigg et al 1997 for Northern hemisphere:
logical :: separate_distrib_for_n_hemisphere=.False. ! Flag to use a separate berg distribution/mass/mass scaling/init thickness for N hemisphere
real, dimension(nclasses) :: initial_mass_n=(/4.58e8, 3.61e9, 1.22e10, 2.91e10, 5.09e10, 7.34e10, 1.15e11, 1.65e11, 2.94e11, 5.59e11/) !for N hemisphere
real, dimension(nclasses) :: distribution_n=(/0.14, 0.15, 0.20, 0.15, 0.08, 0.07, 0.05, 0.05, 0.05, 0.05/) ! for N hemisphere
real, dimension(nclasses) :: mass_scaling_n=(/200, 50, 25, 13, 8, 5, 2, 1, 1, 1/) ! for N hemisphere
real, dimension(nclasses) :: initial_thickness_n=(/80.4, 159.5, 240., 320., 360., 360., 360., 360., 360., 360./) ! for N hemisphere
integer(kind=8) :: debug_iceberg_with_id = -1 ! If positive, monitors a berg with this id
! DEM-mode parameters
!logical :: dem=.false. !if T, run in DEM-mode with angular terms, variable stiffness, etc
logical :: ignore_tangential_force=.false.
real :: poisson=0.3 ! Poisson's ratio
real :: dem_spring_coef=0.
real :: dem_damping_coef=0.1
logical :: dem_shear_for_frac_only=.false. ! If true, DEM shear is calculated for fracture, but zeroed for berg interactions
integer :: dem_beam_test=0 !1=Simply supported beam,2=cantilever beam,3=angular vel test
! Element Interactions
logical :: constant_interaction_LW=.false. ! Always use the initial, globally constant, element length & width during berg interactions
real :: constant_length=0. ! If constant_interaction_LW, the constant length used. If zero in the nml, will be set to max initial L
real :: constant_width=0. ! If constant_interaction_LW, the constant width used. If zero in the nml, will be set to max initial W
logical :: use_spring_for_land_contact=.false. ! Treat contact with masked (land) cells like contact with a static berg
! Footloose calving parameters [England et al (2020) Modeling the breakup of tabular icebergs. Sci. Adv.]
logical :: fl_use_poisson_distribution=.true. ! fl_r is (T) mean of Poisson distribution to determine k, or (F) k=fl_r
logical :: fl_use_perimeter=.false. ! scale number of footloose bergs to calve by perimeter of the parent berg
real :: fl_youngs=1.e8 !< Young's modulus for footloose calculations (Pa)
real :: fl_strength=500. !< yield stress for footloose calculations (kPa)
logical :: fl_use_l_scale=.false. !< footloose based on side melt
real :: fl_l_scale=1. !< scaling factor for footloose based on side melt and/or erosion
logical :: fl_l_scale_erosion_only=.true. !< fl according to erosion - buoyant convection
real :: fl_r=0. ! footloose average number of bergs calved per fl_r_s
real :: fl_r_s=0. ! seconds over which fl_r footloose bergs calve
logical :: displace_fl_bergs=.true. ! footloose berg positions are randomly assigned along edges of parent berg
character(len=11) :: fl_style='new_bergs' ! Evolve footloose bergs individually as 'fl_bits', or as a group with size 'bergy_bits' or 'mean_size'
logical :: fl_bits_erosion_to_bergy_bits=.true. ! Erosion from footloose bits becomes bergy bits
real :: fl_k_scale_by_perimeter=0 ! If greater than 0, scales FL k by (berg perimeter (m))/(this scaling (m))
real :: new_berg_from_fl_bits_mass_thres=huge(0.) ! Create a new berg from FL bits when mass_of_fl_bits exceeds this value
real :: fl_bits_scale_l=0.9 ! For determining dimensions of FL bits berg; FL_bits length    = fl_bits_scale_l * 3*l_b
real :: fl_bits_scale_w=0.9 ! For determining dimensions of FL bits berg; FL_bits width     = fl_bits_scale_w * l_b
real :: fl_bits_scale_t=0.9 ! For determining dimensions of FL bits berg; FL_bits thickness = fl_bits_scale_t * T

namelist /icebergs_nml/ verbose, budget, halo,  traj_sample_hrs, initial_mass, traj_write_hrs, max_bonds, save_short_traj,&
         traj_area_thres, Static_icebergs,distribution, mass_scaling, initial_thickness, verbose_hrs, spring_coef,bond_coef,&
         radial_damping_coef, tangental_damping_coef, only_interactive_forces, rho_bergs, LoW_ratio, debug, really_debug,&
         use_operator_splitting, bergy_bit_erosion_fraction, iceberg_bonds_on, manually_initialize_bonds,&
         ignore_missing_restart_bergs,  parallel_reprod, use_slow_find, sicn_shift, add_weight_to_ocean, passive_mode,&
         ignore_ij_restart, use_new_predictive_corrective, halo_debugging, hexagonal_icebergs, time_average_weight,&
         generate_test_icebergs, speed_limit, fix_restart_dates, use_roundoff_fix, Runge_not_Verlet, interactive_icebergs_on,&
         scale_damping_by_pmag, critical_interaction_damping_on, tang_crit_int_damp_on, require_restart,&
         old_bug_rotated_weights, make_calving_reproduce,restart_input_dir, orig_read, old_bug_bilin,do_unit_tests,&
         grounding_fraction, input_freq_distribution, force_all_pes_traj,allow_bergs_to_roll,set_melt_rates_to_zero,lat_ref,&
         initial_orientation,rotate_icebergs_for_mass_spreading,grid_is_latlon,Lx,use_f_plane,use_old_spreading,&
         grid_is_regular,override_iceberg_velocities,u_override,v_override,add_iceberg_thickness_to_SSH,&
         Iceberg_melt_without_decay,melt_icebergs_as_ice_shelf, Use_three_equation_model,find_melt_using_spread_mass,&
         use_mixed_layer_salinity_for_thermo,utide_icebergs,ustar_icebergs_bg,cdrag_icebergs, pass_fields_to_ocean_model,&
         const_gamma, Gamma_T_3EQ, ignore_traj, debug_iceberg_with_id,use_updated_rolling_scheme, tip_parameter, &
         read_old_restarts, tau_calving, read_ocean_depth_from_file, melt_cutoff,apply_thickness_cutoff_to_gridded_melt,&
         apply_thickness_cutoff_to_bergs_melt, use_mixed_melting, internal_bergs_for_drag, coastal_drift, tidal_drift,&
         mts,new_mts,ewsame,monitor_energy,mts_sub_steps,contact_distance,length_for_manually_initialize_bonds,&
         manually_initialize_bonds_from_radii,contact_spring_coef,fracture_criterion, damage_test_1, uniaxial_test, &
         debug_write,cdrag_grounding,h_to_init_grounding,frac_thres_scaling,frac_thres_n,frac_thres_t,save_bond_traj,&
         remove_unused_bergs,force_convergence,explicit_inner_mts,convergence_tolerance,dem,ignore_tangential_force,poisson,&
         dem_spring_coef,dem_damping_coef,dem_beam_test,constant_interaction_LW,constant_length,constant_width,&
         dem_shear_for_frac_only,use_damage,fl_use_poisson_distribution, fl_use_perimeter, fl_r, fl_r_s,displace_fl_bergs,&
         fl_style,fl_bits_erosion_to_bergy_bits,fl_k_scale_by_perimeter,use_spring_for_land_contact, save_fl_traj,&
         new_berg_from_fl_bits_mass_thres,fl_bits_scale_l,fl_bits_scale_w,fl_bits_scale_t,separate_distrib_for_n_hemisphere,&
         initial_mass_n, distribution_n, mass_scaling_n, initial_thickness_n, fl_use_l_scale, fl_l_scale,&
         fl_l_scale_erosion_only, fl_youngs, fl_strength,  save_all_traj_year, save_nonfl_traj_by_class,&
         save_traj_by_class_start_mass_thres_n, save_traj_by_class_start_mass_thres_s,traj_area_thres_sntbc,&
         traj_area_thres_fl

! Local variables
integer :: ierr, iunit, i, j, id_class, axes3d(3), is,ie,js,je,np
type(icebergs_gridded), pointer :: grd
real :: lon_mod, big_number
logical :: lerr
integer :: stdlogunit, stderrunit
real :: Total_mass_s, Total_mass_n
real :: remaining_dist_s,remaining_dist_n
integer :: last_dist_j, last_dist_j_n
real :: mts_fast_dt=0.0 !Added by Alex
real :: maxlon_c,minlon_c !Added by Alex, for mts verlet periodicity
integer :: k,maxk
real :: dx,dy,dx_dlon,dy_dlat,lat_ref2,lon_ref

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "ice_bergs_framework: "//trim(version)

! Read namelist parameters
 !write(stderrunit,*) 'KID: reading namelist'
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
 !write(stderrunit,*) 'KID: allocating bergs'
  allocate(bergs)
  allocate(bergs%grd)
  grd=>bergs%grd ! For convenience to avoid bergs%grd%X
 !write(stderrunit,*) 'KID: allocating domain'
  allocate(grd%domain)

! Clocks
  bergs%clock=mpp_clock_id( 'Icebergs', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  bergs%clock_mom=mpp_clock_id( 'Icebergs-momentum', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_the=mpp_clock_id( 'Icebergs-thermodyn', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_int=mpp_clock_id( 'Icebergs-interface', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_cal=mpp_clock_id( 'Icebergs-calving', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_com1=mpp_clock_id( 'Icebergs-communication1', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_fl1=mpp_clock_id( 'Icebergs-footloose1', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_com2=mpp_clock_id( 'Icebergs-communication2', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_fl2=mpp_clock_id( 'Icebergs-footloose2', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_ini=mpp_clock_id( 'Icebergs-initialization', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_ior=mpp_clock_id( 'Icebergs-I/O read', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_iow=mpp_clock_id( 'Icebergs-I/O write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_dia=mpp_clock_id( 'Icebergs-diagnostics', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_ini)

! Set up iceberg domain
 !write(stderrunit,*) 'KID: defining domain'
  call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
                           maskmap=maskmap, &
                           xflags=dom_x_flags, xhalo=halo,  &
                           yflags=dom_y_flags, yhalo=halo, name='KID')

  call mpp_define_io_domain(grd%domain, io_layout)

 !write(stderrunit,*) 'KID: get compute domain'
  call mpp_get_compute_domain( grd%domain, grd%isc, grd%iec, grd%jsc, grd%jec )
  call mpp_get_data_domain( grd%domain, grd%isd, grd%ied, grd%jsd, grd%jed )
  call mpp_get_global_domain( grd%domain, grd%isg, grd%ieg, grd%jsg, grd%jeg )

  call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
  call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
  call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
  call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)


  folded_north_on_pe = ((dom_y_flags == FOLD_NORTH_EDGE) .and. (grd%jec == gnj))
 !write(stderrunit,'(a,6i4)') 'KID, icebergs_init: pe,n,s,e,w =',mpp_pe(),grd%pe_N,grd%pe_S,grd%pe_E,grd%pe_W, NULL_PE

 !if (verbose) &
 !write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'KID, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
 !     grd%isc,grd%iec,grd%jsc,grd%jec, &
 !     ' [lon|lat][min|max]=', minval(ice_lon),maxval(ice_lon),minval(ice_lat),maxval(ice_lat)
 !write(stderrunit,*) 'KID, int args = ', mpp_pe(),gni, gnj, layout, axes

 ! Allocate grid of pointers
  allocate( bergs%list(grd%isd:grd%ied, grd%jsd:grd%jed) )
  do j = grd%jsd,grd%jed ; do i = grd%isd,grd%ied
    bergs%list(i,j)%first => null()
  enddo ; enddo

  big_number=1.0E15
 !write(stderrunit,*) 'KID: allocating grid'
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
  allocate( grd%fl_bits_src(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%fl_bits_src(:,:)=0.
  allocate( grd%fl_bits_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%fl_bits_melt(:,:)=0.
  allocate( grd%fl_bits_mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%fl_bits_mass(:,:)=0.
  allocate( grd%fl_bergy_bits_mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%fl_bergy_bits_mass(:,:)=0.
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
  allocate( grd%melt_by_class(grd%isd:grd%ied, grd%jsd:grd%jed, nclasses) ); grd%melt_by_class(:,:,:)=0.
  allocate( grd%melt_buoy_fl(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_buoy_fl(:,:)=0.
  allocate( grd%melt_eros_fl(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_eros_fl(:,:)=0.
  allocate( grd%melt_conv_fl(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_conv_fl(:,:)=0.
  allocate( grd%fl_parent_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%fl_parent_melt(:,:)=0.
  allocate( grd%fl_child_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%fl_child_melt(:,:)=0.
  allocate( grd%tmp(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%tmp(:,:)=0.
  allocate( grd%tmpc(grd%isc:grd%iec, grd%jsc:grd%jec) ); grd%tmpc(:,:)=0.
  allocate( bergs%nbergs_calved_by_class_s(nclasses) ); bergs%nbergs_calved_by_class_s(:)=0
  allocate( bergs%nbergs_calved_by_class_n(nclasses) ); bergs%nbergs_calved_by_class_n(:)=0
  allocate( grd%parity_x(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_x(:,:)=1.
  allocate( grd%parity_y(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_y(:,:)=1.
  allocate( grd%iceberg_counter_grd(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%iceberg_counter_grd(:,:)=0

 !write(stderrunit,*) 'KID: copying grid'
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
  !if(present(ocean_depth)) grd%ocean_depth(grd%isd:grd%ied,grd%jsd:grd%jed)=ocean_depth(:,:)

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


  if ( (mpp_pe() == 0) .and. (Lx .ne. 360.) .and. .not. (grid_is_latlon) ) then
    print *,''
    print *,'pe0 x-domain before periodicity fix'
    write(*,'(f12.2)') (grd%lon(i,grd%jsd), i=grd%isd,grd%ied)
  endif


  if (.not. present(maskmap)) then ! Using a maskmap causes tickles this sanity check
    do j=grd%jsd,grd%jed; do i=grd%isd,grd%ied
      !if (grd%lon(i,j).ge.big_number) write(stderrunit,*) 'bad lon: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1,grd%lon(i,j)
      !if (grd%lat(i,j).ge.big_number) write(stderrunit,*) 'bad lat: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1,grd%lat(i,j)
    enddo; enddo
  endif

  if ((Lx.gt.1E15 ) .and. (mpp_pe().eq.mpp_root_pe())) then
          call error_mesg('KID, framework', 'Model does not enjoy the domain being larger than 1E15. Not sure why. Probably to do with floating point precision.', WARNING)
  endif
  if ((.not. grid_is_latlon) .and. (Lx.eq.360.)) then
    if (mpp_pe().eq.mpp_root_pe())  then
            call error_mesg('KID, framework', 'Since the lat/lon grid is off, the x-direction is being set as non-periodic. Set Lx not equal to 360 override.', WARNING)
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
   ! write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'KID, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
   !      grd%isc,grd%iec,grd%jsc,grd%jec, &
   !      ' [lon|lat][min|max]=', minval(grd%lon),maxval(grd%lon),minval(grd%lat),maxval(grd%lat)

    write(stderrunit,'(a,i3,a,4i4)') 'KID, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
          grd%isc,grd%iec,grd%jsc,grd%jec

    !print *,'minval(grd%lon)',minval(grd%lon)
    !print *,'maxval(grd%lon)',maxval(grd%lon)
    !print *,'minval(grd%lat)',minval(grd%lat)
    !print *,'maxval(grd%lat)',maxval(grd%lat)

    write(stderrunit,'(a,4f10.2)') '[Lon|lat][min|max]=', &
          minval(grd%lon),maxval(grd%lon),minval(grd%lat),maxval(grd%lat)
  endif

 ! if ( (mpp_pe() == 0) .and. (Lx .ne. -1.) .and. .not. (grid_is_latlon) ) then
 !   print *,''
 !   print *,'pe 0 x-domain after periodicity fix'
 !   write(*,'(f12.2)') (grd%lon(i,grd%jsd), i=grd%isd,grd%ied)
 ! end if

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
      call error_mesg('KID,grid defining', 'Latitude contains NaNs', FATAL)
    endif
    if (grd%lon(i,j) .ne. grd%lon(i,j)) then
      write(stderrunit,*) 'Lon not defined properly', mpp_pe(),i,j,grd%lon(i,j)
      call error_mesg('KID, grid defining', 'Longatudes contains NaNs', FATAL)
    endif
  enddo; enddo


  ! If a separate calving distribution for the northern hemisphere is not used, the following northern
  ! hemisphere parameters inherit the values of their respective southern hemisphere parameters
  if (.not. separate_distrib_for_n_hemisphere) then
    initial_mass_n=initial_mass; distribution_n=distribution
    mass_scaling_n=mass_scaling; initial_thickness_n=initial_thickness
  endif


!Added by Alon  - If a freq distribution is input, we have to convert the freq distribution to a mass flux distribution)
if (input_freq_distribution) then
     Total_mass_s=0.; Total_mass_n=0.
     do j=1,nclasses
          Total_mass_s=Total_mass_s+(distribution(j)*initial_mass(j))
          Total_mass_n=Total_mass_n+(distribution_n(j)*initial_mass_n(j))
     enddo
     do j=1,nclasses
           distribution(j)=(distribution(j)*initial_mass(j))/Total_mass_s
           distribution_n(j)=(distribution_n(j)*initial_mass_n(j))/Total_mass_n
     enddo

     if (mpp_pe()==0) then
       print *,'uncorrected distribution_s',distribution
       print *,'uncorrected distribution_n',distribution_n
     endif
     !floating point fix: correct the distributions so that
     !remaining_dist_s and remaining_dist_n in icebergs.F90/accumulate_calving are >= 0
     do j=1,nclasses
       if (distribution(j)>0.) last_dist_j=j
       if (distribution_n(j)>0.) last_dist_j_n=j
     enddo
     remaining_dist_s=1.; remaining_dist_n=1.
     do j=1,last_dist_j-1
       remaining_dist_s=remaining_dist_s-distribution(j)
     enddo
     distribution(last_dist_j)=remaining_dist_s
     do j=1,last_dist_j_n-1
       remaining_dist_n=remaining_dist_n-distribution_n(j)
     enddo
     distribution_n(last_dist_j_n)=remaining_dist_n
     if (mpp_pe()==0) then
       print *,'remaining_dist_s, last_dist_j_s, remaining_dist_n','last_dist_j_n',&
         remaining_dist_s-distribution(last_dist_j),last_dist_j,&
         remaining_dist_n-distribution_n(last_dist_j_n),last_dist_j_n
       print *,'corrected distribution_s',distribution
       print *,'corrected distribution_n',distribution_n
     endif
endif

if ((halo .lt. 3) .and. (rotate_icebergs_for_mass_spreading .and. iceberg_bonds_on) )   then
    halo=3
    call error_mesg('KID, framework', 'Setting iceberg halos =3, since halos must be >= 3 for rotating icebergs for mass spreading', WARNING)
elseif  ((halo .lt. 2) .and. (interactive_icebergs_on .or. iceberg_bonds_on) )   then
    halo=2
    call error_mesg('KID, framework', 'Setting iceberg halos =2, since halos must be >= 2 for interactions', WARNING)
endif

if (interactive_icebergs_on) then
  if (Runge_not_Verlet) then
    !Runge_not_Verlet=.false.  ! Iceberg interactions only with Verlet
    call error_mesg('KID, framework', 'It is unlcear whther interactive icebergs work with Runge Kutta stepping.', WARNING)
  endif
endif
if (.not.interactive_icebergs_on) then
  if (iceberg_bonds_on) then
    !iceberg_bonds_on=.false.
    call error_mesg('KID, framework', 'Interactive icebergs off requires iceberg bonds off (turning bonds off).', WARNING)
  endif
endif
if (.not. iceberg_bonds_on) then
  max_bonds=0
else
  buffer_width=buffer_width+(max_bonds*5) ! Increase buffer width to include bonds being passed between processors
endif
if (save_short_traj) buffer_width_traj=6 ! This is the length of the short buffer used for abrevated traj
if (save_fl_traj) then
  if (fl_r>0) then
    buffer_width_traj=buffer_width_traj+8+2
  else
    buffer_width_traj=buffer_width_traj+4+2
  endif
endif
if (ignore_traj) buffer_width_traj=0 ! If this is true, then all traj files should be ignored

if (use_damage) then
  buffer_width=buffer_width+(max_bonds*1)
  buffer_width_bond_traj=buffer_width_bond_traj+1
endif

if (monitor_energy) then
  buffer_width=buffer_width+11+(max_bonds*6)
  buffer_width_traj=buffer_width_traj+6
  !buffer_width_traj=buffer_width_traj+5 !if include temp terms
  buffer_width_bond_traj=buffer_width_bond_traj+6
endif

if (dem) then
  buffer_width=buffer_width+3+(max_bonds*4)
  buffer_width_traj=buffer_width_traj+3
  buffer_width_bond_traj=buffer_width_bond_traj+4
elseif (fracture_criterion .ne. 'none') then
  buffer_width=buffer_width+1+(max_bonds*3)
  buffer_width_traj=buffer_width_traj+1
  buffer_width_bond_traj=buffer_width_bond_traj+3
  if (fracture_criterion .eq. 'strain_rate') then
    buffer_width=buffer_width+(max_bonds*1)
    buffer_width_bond_traj=buffer_width_bond_traj+1
  elseif (fracture_criterion .eq. 'energy') then
    buffer_width=buffer_width+(max_bonds*1)
    buffer_width_bond_traj=buffer_width_bond_traj+1
  endif
endif

if (iceberg_bonds_on) then
  buffer_width=buffer_width+1
  buffer_width_traj=buffer_width_traj+1
endif

!must use verlet with mts - Alex
if (mts) then
  buffer_width_traj=buffer_width_traj+4 !to accomodate n_bonds,fast-step accel terms, overall accel terms
  buffer_width=buffer_width+17 !to accomodate external forcing (conglom_id,uo,vo,ua,va,ui,vi,ssh_x,ssh_y,sst,sss,cn,hi,accel)

  !if mts_sub_steps is not given in the nml already,
  !then detemine according to the max dt according to the spring coef
  if (mts_sub_steps .eq. -1) then
    mts_fast_dt = 0.3/sqrt(spring_coef)        !critical dt for mts scheme w/ safety multiplier of 0.75
    mts_sub_steps = ceiling(dt/mts_fast_dt)    !the number of substeps
  end if

  mts_fast_dt = dt/mts_sub_steps !adjust mts_fast_dt so that dt is an integer multiple of mts_fast_dt

  if (Runge_not_Verlet) then
    call error_mesg('KID, framework', &
      'Multiple time stepping does not work with Runge Kutta stepping, switching to Verlet.', WARNING)
    Runge_not_Verlet=.false.
  end if
end if

!by default, set contact_spring_coef equal to the spring_coef for bonded interactions, unless
!contact_spring_coef is specified as > 0.
if (contact_spring_coef.le.0.) contact_spring_coef=spring_coef
if (debug_write) then
  if (traj_sample_hrs.ne.traj_write_hrs) then
    call error_mesg('KID, framework', &
      'Debug_write activated: setting traj_sample_hrs=traj_write_hrs', WARNING)
    traj_sample_hrs=traj_write_hrs
  endif
  if (.not. force_all_pes_traj) then
    call error_mesg('KID, framework', &
      'Debug_write activated: setting force_all_pes_traj=.true.', WARNING)
    force_all_pes_traj=.true.
  endif
endif


 ! Parameters
  bergs%dt=dt
  bergs%traj_area_thres=traj_area_thres
  bergs%traj_area_thres_sntbc=traj_area_thres_sntbc
  bergs%traj_area_thres_fl=traj_area_thres_fl
  bergs%traj_sample_hrs=traj_sample_hrs
  bergs%traj_write_hrs=traj_write_hrs
  bergs%save_short_traj=save_short_traj
  bergs%save_fl_traj=save_fl_traj
  bergs%save_all_traj_year=save_all_traj_year
  bergs%save_nonfl_traj_by_class=save_nonfl_traj_by_class
  bergs%save_traj_by_class_start_mass_thres_n=save_traj_by_class_start_mass_thres_n
  bergs%save_traj_by_class_start_mass_thres_s=save_traj_by_class_start_mass_thres_s
  bergs%ignore_traj=ignore_traj
  bergs%verbose_hrs=verbose_hrs
  bergs%grd%halo=halo
  bergs%grd%Lx=Lx
  bergs%grd%grid_is_latlon=grid_is_latlon
  bergs%grd%grid_is_regular=grid_is_regular
  bergs%max_bonds=max_bonds
  bergs%rho_bergs=rho_bergs
  bergs%spring_coef=spring_coef
  bergs%contact_spring_coef=contact_spring_coef !Alex
  bergs%fracture_criterion=fracture_criterion
  bergs%damage_test_1=damage_test_1
  bergs%uniaxial_test=uniaxial_test
  bergs%cdrag_grounding=cdrag_grounding
  bergs%h_to_init_grounding=h_to_init_grounding
  bergs%frac_thres_n=frac_thres_n*frac_thres_scaling
  bergs%frac_thres_t=frac_thres_t*frac_thres_scaling
  bergs%debug_write=debug_write
  bergs%bond_coef=bond_coef
  bergs%radial_damping_coef=radial_damping_coef
  bergs%tangental_damping_coef=tangental_damping_coef
  bergs%LoW_ratio=LoW_ratio
  bergs%separate_distrib_for_n_hemisphere=separate_distrib_for_n_hemisphere
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
  bergs%internal_bergs_for_drag=internal_bergs_for_drag
  bergs%use_mixed_melting=use_mixed_melting
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
  bergs%require_restart=require_restart
  bergs%Static_icebergs=Static_icebergs
  bergs%only_interactive_forces=only_interactive_forces
  bergs%halo_debugging=halo_debugging
  bergs%iceberg_bonds_on=iceberg_bonds_on   !Alon
  bergs%manually_initialize_bonds=manually_initialize_bonds   !Alon
  bergs%length_for_manually_initialize_bonds=length_for_manually_initialize_bonds !Alex
  bergs%manually_initialize_bonds_from_radii=manually_initialize_bonds_from_radii !Alex
  bergs%scale_damping_by_pmag=scale_damping_by_pmag
  bergs%critical_interaction_damping_on=critical_interaction_damping_on   !Alon
  bergs%tang_crit_int_damp_on=tang_crit_int_damp_on
  bergs%interactive_icebergs_on=interactive_icebergs_on   !Alon
  bergs%use_new_predictive_corrective=use_new_predictive_corrective  !Alon
  bergs%grounding_fraction=grounding_fraction
  bergs%add_weight_to_ocean=add_weight_to_ocean
  bergs%use_old_spreading=use_old_spreading
  bergs%debug_iceberg_with_id=debug_iceberg_with_id
  bergs%use_spring_for_land_contact=use_spring_for_land_contact
  bergs%mts=mts
  bergs%mts_fast_dt = mts_fast_dt
  bergs%mts_sub_steps = mts_sub_steps
  bergs%mts_part = 1
  bergs%remove_unused_bergs = remove_unused_bergs
  bergs%contact_distance=contact_distance
  bergs%force_convergence=force_convergence
  bergs%explicit_inner_mts=explicit_inner_mts
  bergs%convergence_tolerance=convergence_tolerance
  ! DEM-mode parameters
  bergs%dem=dem
  bergs%ignore_tangential_force=ignore_tangential_force
  if (bergs%dem) bergs%explicit_inner_mts=.true.
  bergs%poisson=poisson
  bergs%dem_spring_coef=dem_spring_coef
  bergs%dem_damping_coef=dem_damping_coef
  bergs%dem_beam_test=dem_beam_test
  bergs%constant_interaction_LW=constant_interaction_LW
  bergs%constant_length=constant_length
  bergs%constant_width=constant_width
  bergs%dem_shear_for_frac_only=dem_shear_for_frac_only
  ! Footloose calving parameters
  bergs%fl_use_poisson_distribution=fl_use_poisson_distribution
  bergs%fl_use_perimeter=fl_use_perimeter
  bergs%fl_youngs=fl_youngs
  bergs%fl_strength=fl_strength
  bergs%fl_use_l_scale=fl_use_l_scale
  bergs%fl_l_scale=fl_l_scale
  bergs%fl_l_scale_erosion_only=fl_l_scale_erosion_only
  bergs%new_berg_from_fl_bits_mass_thres=new_berg_from_fl_bits_mass_thres
  bergs%fl_r=fl_r
  bergs%fl_r_s=fl_r_s
  bergs%displace_fl_bergs=displace_fl_bergs
  bergs%fl_style=fl_style
  bergs%fl_bits_erosion_to_bergy_bits=fl_bits_erosion_to_bergy_bits
  bergs%fl_k_scale_by_perimeter=fl_k_scale_by_perimeter
  bergs%fl_bits_scale_l=fl_bits_scale_l
  bergs%fl_bits_scale_w=fl_bits_scale_w
  bergs%fl_bits_scale_t=fl_bits_scale_t
  if (bergs%fl_use_l_scale) then
    bergs%fl_use_perimeter=.true.
    bergs%fl_k_scale_by_perimeter=1
    if (.not. bergs%use_operator_splitting) then
      call error_mesg('KID, ice_bergs_framework_init', &
        'use_operator_splitting must be true when fl_use_l_scale=.true.', FATAL)
    endif
  endif

  if (monitor_energy) then
    if (bergs%constant_interaction_LW) then
      call error_mesg('KID, ice_bergs_framework_init', &
        'cannot use monitor_energy with constant_interaction_LW yet', FATAL)
    endif
    if (bergs%dem .and. .not. bergs%ignore_tangential_force) then
      call error_mesg('KID, ice_bergs_framework_init', &
        'monitor_energy does not yet account for dem tangential force terms', FATAL)
    endif
  endif

  if (bergs%contact_distance>0) then
    dx_dlon=1; dy_dlat=1
    if (grd%grid_is_latlon) dy_dlat=pi_180*Rearth
    maxk=0
    do j=grd%jsd,grd%jed; do i=grd%isd,grd%ied
      if (grd%grid_is_latlon) then
        lat_ref2=grd%lat(i,j)
        dx_dlon=pi_180*Rearth*cos(lat_ref2*pi_180)
      endif
      lon_ref=grd%lon(i,j)
      k=0
      do while((k+i)<grd%ied)
        k=k+1
        dx=(grd%lon(k+i,j)-lon_ref)*dx_dlon
        if (k>maxk) maxk=k
        if (dx>=contact_distance) exit
      enddo
    enddo;enddo

    if (grd%grid_is_latlon) dy_dlat=pi_180*Rearth
    dy=(grd%lat(grd%isc,grd%jsc+1)-grd%lat(grd%isc,grd%jsc))*dy_dlat

    bergs%contact_cells_lon = max(maxk,1)
    bergs%contact_cells_lat = max(int(ceiling(contact_distance/dy)),1)
  else
    bergs%contact_cells_lon = 1
    bergs%contact_cells_lat = 1
  endif

  !print *,'# contact cells lon/lat',bergs%contact_cells_lon,bergs%contact_cells_lat

  !necessary?
  if (.not. mts) then
    if ((halo-1)<bergs%contact_cells_lon .or. (halo-1)<bergs%contact_cells_lat) then
      if (mpp_pe()==0) then
        write(stderrunit,'(a,i3)') 'KID, icebergs_init: halo width minus 1 ',halo - 1
        write(stderrunit,'(a,f8.2,a,2i3)') 'KID, icebergs_init:  contact distance ',bergs%contact_distance,&
          ' num contact cells lon and lat',bergs%contact_cells_lon,bergs%contact_cells_lat

        call error_mesg('KID, ice_bergs_framework_init', &
          'halo width must be increased to accomodate specified contact distance!!!', FATAL)
      endif
    endif
  endif
  allocate( bergs%initial_mass_s(nclasses) ); bergs%initial_mass_s(:)=initial_mass(:)
  allocate( bergs%distribution_s(nclasses) ); bergs%distribution_s(:)=distribution(:)
  allocate( bergs%mass_scaling_s(nclasses) ); bergs%mass_scaling_s(:)=mass_scaling(:)
  allocate( bergs%initial_thickness_s(nclasses) ); bergs%initial_thickness_s(:)=initial_thickness(:)
  allocate( bergs%initial_width_s(nclasses) )
  allocate( bergs%initial_length_s(nclasses) )
  bergs%initial_width_s(:)=sqrt(initial_mass(:)/(LoW_ratio*rho_bergs*initial_thickness(:)))
  bergs%initial_length_s(:)=LoW_ratio*bergs%initial_width_s(:)

  allocate( bergs%initial_mass_n(nclasses) ); bergs%initial_mass_n(:)=initial_mass_n(:)
  allocate( bergs%distribution_n(nclasses) ); bergs%distribution_n(:)=distribution_n(:)
  allocate( bergs%mass_scaling_n(nclasses) ); bergs%mass_scaling_n(:)=mass_scaling_n(:)
  allocate( bergs%initial_thickness_n(nclasses) ); bergs%initial_thickness_n(:)=initial_thickness_n(:)
  allocate( bergs%initial_width_n(nclasses) )
  allocate( bergs%initial_length_n(nclasses) )
  bergs%initial_width_n(:)=sqrt(initial_mass_n(:)/(LoW_ratio*rho_bergs*initial_thickness_n(:)))
  bergs%initial_length_n(:)=LoW_ratio*bergs%initial_width_n(:)

  bergs%writeandstop=.false.
  grd%coastal_drift = coastal_drift
  grd%tidal_drift = tidal_drift

  !needed for periodicity
  maxlon_c = grd%lon(grd%iec,grd%jec); minlon_c = grd%lon(grd%isc-1,grd%jsc-1)
  call mpp_max(maxlon_c); call mpp_min(minlon_c)
  grd%maxlon_c=maxlon_c; grd%minlon_c=minlon_c

  if (read_old_restarts) call error_mesg('KID, ice_bergs_framework_init', 'Setting "read_old_restarts=.true." is obsolete and does nothing!', WARNING)

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
  grd%id_fl_bits_src=register_diag_field('icebergs', 'fl_bits_src', axes, Time, &
     'Mass source of footloose bits', 'kg/(m^2*s)')
  grd%id_fl_bits_melt=register_diag_field('icebergs', 'fl_bits_melt', axes, Time, &
     'Melt rate of footloose bits', 'kg/(m^2*s)')
  grd%id_fl_bits_mass=register_diag_field('icebergs', 'fl_bits_mass', axes, Time, &
    'Footloose bits density field', 'kg/(m^2)')
  grd%id_fl_bergy_bits_mass=register_diag_field('icebergs', 'fl_bergy_bits_mass', axes, Time, &
    'Footloose bergy bits density field', 'kg/(m^2)')
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
     'Accumulated ice mass by class (z-axis labels correspond to Southern hemisphere classes)', &
     'kg')
  grd%id_real_calving=register_diag_field('icebergs', 'real_calving', axes3d, Time, &
     'Calving into iceberg class (z-axis labels correspond to Southern hemisphere classes)', &
     'kg/s')
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
  grd%id_melt_by_class=register_diag_field('icebergs', 'melt_by_class', axes3d, Time, &
     'Total ice melt (bergs+bits+FL_bits) by class (z-axis labels correspond to Southern hemisphere classes)', &
     'kg/(m^2*s)')
  grd%id_melt_buoy_fl=register_diag_field('icebergs', 'melt_buoy_fl', axes, Time, &
     'Buoyancy component of footloose iceberg melt rate', 'kg/(m^2*s)')
  grd%id_melt_eros_fl=register_diag_field('icebergs', 'melt_eros_fl', axes, Time, &
     'Erosion component of footloose iceberg melt rate', 'kg/(m^2*s)')
  grd%id_melt_conv_fl=register_diag_field('icebergs', 'melt_conv_fl', axes, Time, &
     'Convective component of footloose iceberg melt rate', 'kg/(m^2*s)')
  grd%id_fl_parent_melt=register_diag_field('icebergs', 'fl_parent_melt', axes, Time, &
     'Melt rate of footloose parent bergs', 'kg/(m^2*s)')
  grd%id_fl_child_melt=register_diag_field('icebergs', 'fl_child_melt', axes, Time, &
     'Melt rate of footloose child bergs', 'kg/(m^2*s)')

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
   if (unit_tests(bergs)) call error_mesg('KID, icebergs_init', 'Unit tests failed!', FATAL)
  endif

 !write(stderrunit,*) 'KID: done'
  call mpp_clock_end(bergs%clock_ini)
  call mpp_clock_end(bergs%clock)

end subroutine ice_bergs_framework_init

!> Adjust berg dates to allow use of restarts from later dates
subroutine offset_berg_dates(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
type(time_type), intent(in) :: Time !< Model time
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
    'KID: Bergs found with creation dates after model date! Adjusting berg dates by ',yr_offset,' years'
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

!> Moves icebergs between lists if they have moved from cell to cell
subroutine move_berg_between_cells(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
! Local variables
type(icebergs_gridded), pointer :: grd => null()
type(iceberg), pointer :: moving_berg => null(), this => null()
integer :: grdi, grdj

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

!> Populates the halo lists with bergs from neighbor processers
subroutine update_halo_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
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

  ! For debugging, MP1
  if (halo_debugging) then
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        write(stderrunit,'(a,5i8)') 'A', this%id, mpp_pe(), int(this%halo_berg), grdi, grdj
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

  ! if (bergs%mts) then
  !   !remove conglomerate bergs outside halo
  !   do grdj = grd%jsc,grd%jec ;    do grdi = grd%isc,grd%iec
  !     this=>bergs%list(grdi,grdj)%first
  !     do while (associated(this))
  !       kick_the_bucket=>this
  !       this=>this%next
  !       if (kick_the_bucket%halo_berg .ne. 0) then
  !         call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
  !       end if
  !     enddo
  !   enddo;enddo
  ! endif

  call mpp_sync_self()

  ! For debugging
  if (halo_debugging) then
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        write(stderrunit,'(a,5i8)') 'B', this%id, mpp_pe(), int(this%halo_berg), grdi, grdj
        this=>this%next
      enddo
    enddo; enddo
  endif
  if (debug) then
    nbergs_start=count_bergs(bergs, with_halos=.true.)
  endif

  call mpp_sync_self()

  ! Step 2: Updating the halos  - This code is mostly copied from send_to_other_pes

  ! Find number of bergs that headed east/west
  nbergs_to_send_e=0
  nbergs_to_send_w=0
  ! Bergs on eastern side of the processor
  do grdj = grd%jsc,grd%jec ; do grdi = grd%iec-halo_width+2,grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
    !write(stderrunit,*)  'sending east', this%id, this%ine, this%jne, mpp_pe()
      kick_the_bucket=>this
      this=>this%next
      nbergs_to_send_e=nbergs_to_send_e+1
      current_halo_status=kick_the_bucket%halo_berg
      kick_the_bucket%halo_berg=1.
      call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_e, nbergs_to_send_e, bergs%max_bonds)
      kick_the_bucket%halo_berg=current_halo_status
    enddo
  enddo; enddo

  ! Bergs on the western side of the processor
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

  ! Bergs on north side of the processor
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

  ! Bergs on south side of the processor
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

  ! For debugging
  if (halo_debugging) then
    call mpp_sync_self()
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        write(stderrunit,'(a,5i8)') 'C', this%id, mpp_pe(), int(this%halo_berg),  grdi, grdj
        this=>this%next
      enddo
    enddo; enddo
    call show_all_bonds(bergs)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (debug) then
    nbergs_end=count_bergs(bergs, with_halos=.true.)
    i=nbergs_rcvd_from_n+nbergs_rcvd_from_s+nbergs_rcvd_from_e+nbergs_rcvd_from_w &
     -nbergs_to_send_n-nbergs_to_send_s-nbergs_to_send_e-nbergs_to_send_w
    if (nbergs_end-(nbergs_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_end=',nbergs_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_start=',nbergs_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: error=',nbergs_end-(nbergs_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_to_send_n=',nbergs_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_to_send_s=',nbergs_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_to_send_e=',nbergs_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_to_send_w=',nbergs_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_rcvd_from_n=',nbergs_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_rcvd_from_s=',nbergs_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_rcvd_from_e=',nbergs_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, update_halos: nbergs_rcvd_from_w=',nbergs_rcvd_from_w,' on PE',mpp_pe()
      !call error_mesg('KID, update_halos:', 'We lost some bergs!', FATAL)
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
      write(stderrunit,'(a,i4)') 'KID, update_halos: # of bergs outside computational domain = ',i
      call error_mesg('KID, update_halos:', 'there are bergs still in halos!', FATAL)
    endif ! root_pe
  endif ! debug
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

end subroutine update_halo_icebergs

!> For the multiple-timestepping velocity verlet scheme, populates the current PE with the following
!! bergs from neighboring PEs: halo bergs, bergs that comprise any conglomerate that overlaps both
!! PEs, and bergs within the contact distance of these halo and conglomerate bergs.
Subroutine transfer_mts_bergs(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
! Local variables
type(icebergs_gridded), pointer :: grd
integer :: stderrunit
integer :: grdi, grdj,i
type(iceberg), pointer :: this,kick_the_bucket
integer :: nbergs_to_send_e, nbergs_to_send_w
integer :: nbergs_to_send_n, nbergs_to_send_s

  grd=>bergs%grd ! for convenience
  stderrunit = stderr() ! Get the stderr unit number

  ! Step 1: Clear the current halos
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

  !Remove contact-only/conglom bergs outside halo. Reset conglom IDs to zero for all other bergs
  do grdj = grd%jsc,grd%jec ;    do grdi = grd%isc,grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%halo_berg .ne. 0) then
        kick_the_bucket=>this
        this=>this%next
        call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
      else
        this%conglom_id=0
        this=>this%next
      end if
    enddo
  enddo;enddo

  !Copy bergs between PEs
  do i = 1,2 !run twice to account for diagonal transfers and guarantee robust transfers of conglomerates
    call connect_all_bonds(bergs,ignore_unmatched=.true.)
    nbergs_to_send_e=0; nbergs_to_send_w=0; nbergs_to_send_n=0; nbergs_to_send_s=0

    call mts_pack_in_dir(bergs,nbergs_to_send_e,"e")
    call mts_pack_in_dir(bergs,nbergs_to_send_w,"w")
    call mts_pack_in_dir(bergs,nbergs_to_send_n,"n")
    call mts_pack_in_dir(bergs,nbergs_to_send_s,"s")

    call mts_send_and_receive(bergs,nbergs_to_send_e,nbergs_to_send_w,nbergs_to_send_n,nbergs_to_send_s)
  enddo

  if (debug) then
    call connect_all_bonds(bergs,ignore_unmatched=.false.)
  else
    call connect_all_bonds(bergs,ignore_unmatched=.true.)
  endif

  call set_conglom_ids(bergs)
  if (bergs%remove_unused_bergs) call mts_remove_unused_bergs(bergs)

  ! For debugging
  if (bergs%halo_debugging) then
    do grdj = grd%jsd,grd%jed ;  do grdi = grd%isd,grd%ied
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        write(stderrunit,'(a,5i8)') 'transfer_mts_bergs', this%id, mpp_pe(), int(this%halo_berg),  grdi, grdj
        this=>this%next
      enddo
    enddo; enddo
    call show_all_bonds(bergs)
  endif

end subroutine transfer_mts_bergs

!> For the MTS scheme, packs bergs for transfers between PEs to the N,S,E, and W
subroutine mts_pack_in_dir(bergs, nbergs_to_send, dir)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer :: nbergs_to_send !< Counter for number of bergs to send
  character(len=1) :: dir !< Indicates PE to the north,south,east,or west
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg
  integer :: halo_width,nc_x,nc_y
  real :: cds !<contact distance squared
  logical :: edgecontact
  integer :: is,ie,js,je,grdi,grdj
  integer :: ics,ice,jcs,jce
  integer :: inhs,inhe,jnhs,jnhe
  real :: pfix !<used to adjust for periodicity
  real :: rhc(4)
  real :: clat,clon,dlat,dlon,r_dist
  integer :: current_conglom_id
  real :: current_halo_id

  !---INITIALIZATION---!
  grd=>bergs%grd
  halo_width=bergs%grd%halo
  nc_x=bergs%contact_cells_lon; nc_y=bergs%contact_cells_lat
  cds=bergs%contact_distance**2
  edgecontact=.false.
  pfix=0.0

  select case (dir)
  case ("e")
    !congloms overlapping this range of cells are always sent:
    is=grd%iec-halo_width+2; ie=grd%ied; js=grd%jsd; je=grd%jed
    if (grd%Lx>0.0 .and. grd%lon(grd%iec,grd%jec).eq.grd%maxlon_c) pfix=grd%Lx
    !lat/lon range of halo cells
    rhc(1)=grd%lon(is-1,js);    rhc(2)=grd%lat(is-1,js)    !min corner
    rhc(3)=grd%lon(grd%iec,je); rhc(4)=grd%lat(grd%iec,je) !max corner
    !range for non-halo cells (for contact with conglom)
    inhs=grd%isd; inhe=is-1; jnhs=js; jnhe=je
    !range for cells w/in contact dist of edge of comp domain:
    ics=max(grd%iec-nc_x,grd%isd+1)
    if (ics<is) then
      edgecontact=.true.; ice=grd%iec; jcs=js; jce=je
      clon=grd%lon(grd%iec,grd%jec) !closest lon from the east PE computational domain
    endif
  case ("w")
    is=grd%isd; ie=grd%isc+halo_width-1; js=grd%jsd; je=grd%jed
    if (grd%Lx>0.0 .and. grd%lon(grd%isc-1,grd%jsc-1).eq.grd%minlon_c) pfix=-grd%Lx
    rhc(1)=grd%lon(grd%isc-1,js); rhc(2)=grd%lat(grd%isc-1,js)
    rhc(3)=grd%lon(ie,je);        rhc(4)=grd%lat(ie,je)
    inhs=ie+1; inhe=grd%ied; jnhs=js; jnhe=je
    ice=min(grd%isc+nc_x,grd%ied)
    if (ice>ie) then
      edgecontact=.true.; ics=grd%isc; jcs=js; jce=je
      clon=grd%lon(grd%isc-1,grd%jsc-1)
    endif
  case ("n")
    is=grd%isd; ie=grd%ied; js=grd%jec-halo_width+2; je=grd%jed
    rhc(1)=grd%lon(is,js-1);    rhc(2)=grd%lat(is,js-1)
    rhc(3)=grd%lon(ie,grd%jec); rhc(4)=grd%lat(ie,grd%jec)
    inhs=is; inhe=ie; jnhs=grd%jsd; jnhe=js-1
    jcs=max(grd%jec-nc_y,grd%jsd+1)
    if (jcs<js) then
      edgecontact=.true.; ics=is; ice=ie; jce=grd%jec
      clat=grd%lat(grd%iec,grd%jec)
    endif
  case ("s")
    is=grd%isd; ie=grd%ied; js=grd%jsd; je=grd%jsc+halo_width-1
    rhc(1)=grd%lon(is,grd%jsc-1); rhc(2)=grd%lat(is,grd%jsc-1)
    rhc(3)=grd%lon(ie,je);        rhc(4)=grd%lat(ie,je)
    inhs=is; inhe=ie; jnhs=je+1; jnhe=grd%jed
    jce=min(grd%jsc+nc_y,grd%jed)
    if (jce>je) then
      edgecontact=.true.; ics=is; ice=ie; jcs=grd%jsc
      clat=grd%lat(grd%isc-1,grd%jsc-1)
    endif
  case default
    call error_mesg('mts_pack_in_dir', 'dir not specified correctly!', FATAL)
  end select

  !---MAIN ALGORITHM---:
  !
  !The packing procedure is split into the three routines below. A note on bookkeeping:
  !
  !-In routines 1 & 2, the sign of berg%id is used to mark if a berg has already been processed
  !
  !-The packing routines rely on knowing which neighboring PEs already contain (or will contain,
  ! after sending) copies of a berg. The locations of copies are tracked using bergs%conglom_id, which
  ! initially equals zero at the start of subroutine transfer_mts_bergs. Whenever a copy of a berg is
  ! packed to send to a neighboring PE, the original berg’s conglom_id=conglom_id+a, where a is an integer
  ! associated with the direction of the PE that the copy will reside (East=4, West=8, North=2, South=1).
  ! The copy’s conglom_id=conglom_id+b, where b is the integer for the direction of the original
  ! berg’s PE relative to the PE to which the copy will be sent. From a berg’s growing conglom_id, it is
  ! possible to backtrack which neighboring PEs contain (or will contain after sending) a copy of the
  ! berg by using the function mts_berg_sent, which is called in all routines below. Note: during the MTS
  ! scheme, conglom_id is reset and repurposed for marking bergs that belong to the same conglomerate.

  !1. Pack bergs to send to the neighboring PE in the given dir if they are in the
  !   corresponding halo cells, or connected to bergs in the halo cells as a conglomerate.
  do grdj=js,je; do grdi=is,ie
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg))
      if (berg%id>0) then !bergs that have already been processed are marked with negative berg%id
        call mts_mark_and_pack_halo_and_congloms(bergs,berg,dir,nbergs_to_send,pfix,rhc)
      endif
      berg=>berg%next
    enddo
  enddo;enddo

  !2. Find and pack bergs within contact distance of the bergs on the edge of the conglom
  do grdj=js,je; do grdi=is,ie
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg))
      if (berg%id<0) then !berg was processed in 1. Will be reset to positive berg%id here.
        call mts_pack_contact_bergs(bergs,berg,dir,pfix,nbergs_to_send,inhs,inhe,jnhs,jnhe)
      endif
      berg=>berg%next
    enddo
  enddo;enddo

  !3. Find and pack bergs within contact distance to the closest possible coords on the neighbor
  !   PE, to account for the possibility of a berg being located there.
  if (edgecontact) then
    do grdj=jcs,jce; do grdi=ics,ice
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg))
        !has the berg already been sent/packed in this dir?
        if (.not. mts_berg_sent(berg%conglom_id,dir)) then
          if (dir.eq."e" .or. dir.eq."w") then
            dlon=berg%lon-clon; dlat=0.0 !lat_ref=berg%lat
            if (grd%grid_is_latlon) dlon=dlon*pi_180*Rearth*cos(berg%lat*pi_180)
          else
            dlat=berg%lat-clat; dlon=0.0
            if (grd%grid_is_latlon) dlat=dlat*pi_180*Rearth
          endif
          r_dist=dlat**2+dlon**2
          if (r_dist<cds) then
            current_conglom_id=berg%conglom_id; current_halo_id=berg%halo_berg
            berg%halo_berg=10
            nbergs_to_send=nbergs_to_send+1
            select case (dir)
            case ("s")
              berg%conglom_id=2
              call pack_berg_into_buffer2(berg, bergs%obuffer_s, nbergs_to_send, bergs%max_bonds)
              berg%conglom_id=current_conglom_id+1
            case ("n")
              berg%conglom_id=1
              call pack_berg_into_buffer2(berg, bergs%obuffer_n, nbergs_to_send, bergs%max_bonds)
              berg%conglom_id=current_conglom_id+2
            case ("w")
              berg%lon=berg%lon-pfix; berg%conglom_id=4
              berg%lon_prev=berg%lon_prev-pfix
              call pack_berg_into_buffer2(berg, bergs%obuffer_w, nbergs_to_send, bergs%max_bonds)
              berg%lon=berg%lon+pfix; berg%conglom_id=current_conglom_id+8
              berg%lon_prev=berg%lon_prev+pfix
            case ("e")
              berg%lon=berg%lon-pfix; berg%conglom_id=8
              berg%lon_prev=berg%lon_prev-pfix
              call pack_berg_into_buffer2(berg, bergs%obuffer_e, nbergs_to_send, bergs%max_bonds)
              berg%lon=berg%lon+pfix; berg%conglom_id=current_conglom_id+4
              berg%lon_prev=berg%lon_prev+pfix
            end select
            berg%halo_berg=current_halo_id
          endif
        endif
        berg=>berg%next
      enddo
    enddo;enddo
  endif
end subroutine mts_pack_in_dir

!> Find, mark, and pack the halo and conglomerate bergs
recursive subroutine mts_mark_and_pack_halo_and_congloms(bergs, berg, dir, nbergs_to_send,pfix, rhc)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Berg to pack
  character(len=1) :: dir !< north,south,east,or west
  integer :: nbergs_to_send !< Number of bergs to send
  real :: pfix !<for periodicity
  real :: rhc(4) !<lat/lon bounds of halo region for receiving cell
  ! Local variables
  type(iceberg), pointer :: other_berg
  type(bond) , pointer :: current_bond
  integer :: k !<bond counter
  integer :: current_conglom_id
  real :: current_halo_id

  !pack the berg for transfer if it has not been packed already
  if (.not. mts_berg_sent(berg%conglom_id,dir)) then
    current_halo_id=berg%halo_berg; current_conglom_id=berg%conglom_id
    nbergs_to_send=nbergs_to_send+1

    if (berg%lon<rhc(1).or.berg%lat<rhc(2).or.berg%lon>rhc(3).or.berg%lat>rhc(4)) then
      berg%halo_berg=max(2.,current_halo_id) !outside neighboring PE halo
    else
      berg%halo_berg=1 !inside neighboring PE halo
    endif

    select case (dir)
    case ("e")
      berg%conglom_id=8;                    berg%lon=berg%lon-pfix; berg%lon_prev=berg%lon_prev-pfix
      call pack_berg_into_buffer2(berg,bergs%obuffer_e, nbergs_to_send, bergs%max_bonds)
      berg%conglom_id=current_conglom_id+4; berg%lon=berg%lon+pfix; berg%lon_prev=berg%lon_prev+pfix
    case ("w")
      berg%conglom_id=4;                    berg%lon=berg%lon-pfix; berg%lon_prev=berg%lon_prev-pfix
      call pack_berg_into_buffer2(berg,bergs%obuffer_w, nbergs_to_send, bergs%max_bonds)
      berg%conglom_id=current_conglom_id+8; berg%lon=berg%lon+pfix; berg%lon_prev=berg%lon_prev+pfix
    case ("n")
      berg%conglom_id=1
      call pack_berg_into_buffer2(berg,bergs%obuffer_n, nbergs_to_send, bergs%max_bonds)
      berg%conglom_id=current_conglom_id+2
    case ("s")
      berg%conglom_id=2
      call pack_berg_into_buffer2(berg,bergs%obuffer_s, nbergs_to_send, bergs%max_bonds)
      berg%conglom_id=current_conglom_id+1
    end select

    berg%halo_berg=current_halo_id
  endif

  berg%id=-berg%id !marks that the berg has been processed

  !process connected bergs:
  k=0 !bond counter
  current_bond=>berg%first_bond
  do while (associated(current_bond))
    k=k+1
    if  (associated(current_bond%other_berg)) then
      other_berg=>current_bond%other_berg
      if (other_berg%id>0) then
        call mts_mark_and_pack_halo_and_congloms(bergs,other_berg,dir,nbergs_to_send,pfix,rhc)
      endif
    endif
    current_bond=>current_bond%next_bond
  enddo

  if (k.ne.bergs%max_bonds) then !mark the berg element as belonging to the outer layer of a conglomerate
    berg%conglom_id=-berg%conglom_id
  endif

end subroutine mts_mark_and_pack_halo_and_congloms

!> Pack the bergs within contact distance of conglom/halo bergs
recursive subroutine mts_pack_contact_bergs(bergs, berg, dir, pfix, nbergs_to_send, inhs, inhe, jnhs, jnhe)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: berg !< Berg to find contact with
  character(len=1) :: dir !< north,south,east,or west
  real :: pfix !< for periodicity
  integer :: nbergs_to_send !< Number of bergs to send
  integer :: inhs !< Start i bounds of non-halo domain, i.e. where contact bergs may exist
  integer :: inhe !< End i bounds of non-halo domain, i.e. where contact bergs may exist
  integer :: jnhs !< Start j bounds of non-halo domain, i.e. where contact bergs may exist
  integer :: jnhe !< End j bounds of non-halo domain, i.e. where contact bergs may exist
  ! Local variables
  type(bond) , pointer :: current_bond
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: other_berg
  integer :: nc_x,nc_y,grdj,grdi,js,je,is,ie
  integer :: current_conglom_id
  logical :: radial_contact
  real :: R1,R2,dlon,dlat,lat_ref,dx_dlon,dy_dlat,r_dist,crit_dist,current_halo_id
  real :: rdenom
  logical :: pack_contacts

  pack_contacts=.false.
  radial_contact=.false.
  berg%id=-berg%id !unmark

  if (berg%conglom_id<0) then
    berg%conglom_id=-berg%conglom_id !unmark
    grd=>bergs%grd
    nc_x=bergs%contact_cells_lon; nc_y=bergs%contact_cells_lat
    grdj=berg%jne; grdi=berg%ine
    js=max(grdj-nc_y,jnhs); je=min(grdj+nc_y,jnhe)
    is=max(grdi-nc_x,inhs); ie=min(grdi+nc_x,inhe)
    !only pack contacs on computational domain (bergs on data domain already sent)
    if (.not. (je<js .or. ie<is)) pack_contacts=.true.
  endif

  if (pack_contacts) then
    if (bergs%hexagonal_icebergs) then
      rdenom=1./(2.*sqrt(3.))
    else
      if (bergs%iceberg_bonds_on) then
        rdenom=1./4.
      else
        rdenom=1./pi
      endif
    endif
    if (nc_x.eq.1 .and. nc_y.eq.1) then
      radial_contact=.true. !contact based on berg radii can occur
      R1=sqrt(berg%length*berg%width*rdenom) !radius of current berg
    else
      crit_dist=bergs%contact_distance**2
    endif

    do grdj=js,je; do grdi=is,ie
      other_berg=>bergs%list(grdi,grdj)%first
      do while (associated(other_berg)) ! inner berg loop
        if (other_berg%id>0) then
          !has the berg already been sent/packed in this dir?
          if (.not. mts_berg_sent(other_berg%conglom_id,dir)) then
            if (radial_contact) then
              R2=sqrt(other_berg%length*other_berg%width*rdenom)
              crit_dist=max(R1+R2,bergs%contact_distance)**2
            endif

            if (grd%grid_is_latlon) then
              lat_ref=0.5*(berg%lat+other_berg%lat)
              dx_dlon=pi_180*Rearth*cos(lat_ref*pi_180); dy_dlat=pi_180*Rearth
              dlon=berg%lon-other_berg%lon; dlat=berg%lat-other_berg%lat
              r_dist=(dlon*dx_dlon)**2 + (dlat*dy_dlat)**2
            else
              r_dist=(berg%lon-other_berg%lon)**2 + (berg%lat-other_berg%lat)**2
            endif

            if (r_dist<crit_dist) then       !other_berg is within the contact distance
              current_halo_id=other_berg%halo_berg; current_conglom_id=other_berg%conglom_id
              other_berg%halo_berg=10
              nbergs_to_send=nbergs_to_send+1
              select case (dir)
              case ("s")
                other_berg%conglom_id=2
                call pack_berg_into_buffer2(other_berg,bergs%obuffer_s, nbergs_to_send, bergs%max_bonds)
                other_berg%conglom_id=current_conglom_id+1
              case ("n")
                other_berg%conglom_id=1
                call pack_berg_into_buffer2(other_berg,bergs%obuffer_n, nbergs_to_send, bergs%max_bonds)
                other_berg%conglom_id=current_conglom_id+2
              case ("w")
                other_berg%lon=other_berg%lon-pfix; other_berg%conglom_id=4
                other_berg%lon_prev=other_berg%lon_prev-pfix
                call pack_berg_into_buffer2(other_berg,bergs%obuffer_w, nbergs_to_send, bergs%max_bonds)
                other_berg%lon=other_berg%lon+pfix; other_berg%conglom_id=current_conglom_id+8
                other_berg%lon_prev=other_berg%lon_prev+pfix
              case ("e")
                other_berg%lon=other_berg%lon-pfix; other_berg%conglom_id=8
                other_berg%lon_prev=other_berg%lon_prev-pfix
                call pack_berg_into_buffer2(other_berg,bergs%obuffer_e, nbergs_to_send, bergs%max_bonds)
                other_berg%lon=other_berg%lon+pfix; other_berg%conglom_id=current_conglom_id+4
                other_berg%lon_prev=other_berg%lon_prev+pfix
              end select
              other_berg%halo_berg=current_halo_id
            endif
          endif
        endif
        other_berg=>other_berg%next
      enddo !inner berg loop
    enddo;enddo
  endif

  !process connected bergs:
  current_bond=>berg%first_bond
  do while (associated(current_bond))
    if  (associated(current_bond%other_berg)) then
      other_berg=>current_bond%other_berg
      if (other_berg%id<0) then
        call mts_pack_contact_bergs(bergs,other_berg,dir,pfix,nbergs_to_send,inhs,inhe,jnhs,jnhe)
      endif
    endif
    current_bond=>current_bond%next_bond
  enddo
end subroutine mts_pack_contact_bergs

!> For MTS scheme, determines if a berg has been packed to send to the neighboring PE in a given direction
!! based on its evolving conglom_id from subroutines transfer_mts_bergs and mts_pack_in_dir
logical function mts_berg_sent(conglom_id, dir)
  integer :: conglom_id !< Id of conglomerate
  character(len=1) :: dir !< Direction for sending

  mts_berg_sent=.false.
  select case (dir)
  case ("e")
    if (mod(conglom_id,8)>3) mts_berg_sent=.true.
    if (ewsame) then
      if (conglom_id>7) mts_berg_sent=.true.
    endif
  case ("w")
    if (conglom_id>7) mts_berg_sent=.true.
    if (ewsame) then
      if (mod(conglom_id,8)>3) mts_berg_sent=.true.
    endif
  case ("n")
    if (mod(conglom_id,4)>1) mts_berg_sent=.true.
  case ("s")
    if (mod(conglom_id,2).eq.1) mts_berg_sent=.true.
  end select
end function mts_berg_sent

!> Clear and reassign conglom ids
subroutine set_conglom_ids(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: this
  integer :: grdi,grdj,new_conglom_id

  grd=>bergs%grd

  !Reset all conglom ids to zero
  do grdj = grd%jsd,grd%jed; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      this%conglom_id=0
      this=>this%next
    enddo
  enddo;enddo

  !Set the conglom_id for conglomerate bergs that overlap the computational domain
  new_conglom_id=0
  do grdj = grd%jsc,grd%jec; do grdi=grd%isc,grd%iec
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%conglom_id.eq.0) then
        new_conglom_id=new_conglom_id+1
        this%conglom_id=new_conglom_id
        call label_conglomerates(this,new_conglom_id)
      endif
      this=>this%next
    enddo
  enddo;enddo

end subroutine set_conglom_ids

!> Assign initial velocities and energies for energy tracking tests
subroutine energy_tests_init(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: this
  integer :: grdi,grdj
  real :: vol
  type(bond), pointer :: current_bond

  grd=>bergs%grd

  !Reset all conglom ids to zero
  do grdj = grd%jsd,grd%jed; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))

      this%uvel_old=this%uvel; this%uvel_prev=this%uvel
      this%vvel_old=this%vvel; this%vvel_prev=this%vvel

      this%Ee=0.0; this%Ed=0.0; this%Eext=0.0; this%Ee_contact=0.0; this%Ed_contact=0.0; this%Efrac=0.0
      this%Ee_temp=0.0; this%Ed_temp=0.0; this%Eext_temp=0.0; this%Ee_contact_temp=0.0; this%Ed_contact_temp=0.0

      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        current_bond%Ee=0.0; current_bond%Ed=0.0
        current_bond%axn_fast=0.0; current_bond%ayn_fast=0.0
        current_bond%bxn_fast=0.0; current_bond%byn_fast=0.0
        current_bond=>current_bond%next_bond
      enddo

      this=>this%next
    enddo
  enddo;enddo
end subroutine energy_tests_init

!> Identify all bergs in the same conglomerate as the given berg, and assign them the given conglom_id
recursive subroutine label_conglomerates(berg, new_conglom_id)
  type(iceberg), pointer :: berg !< Berg to start search from
  integer :: new_conglom_id !< Conglomerate id to assign
  ! Local variables
  type(iceberg), pointer :: other_berg
  type(bond) , pointer :: current_bond

  current_bond=>berg%first_bond
  do while (associated(current_bond))
    if  (associated(current_bond%other_berg)) then
      other_berg=>current_bond%other_berg
      if (other_berg%conglom_id.ne.new_conglom_id) then
        other_berg%conglom_id=new_conglom_id
        if (other_berg%halo_berg==10) other_berg%halo_berg=2
        call label_conglomerates(other_berg,new_conglom_id)
      endif
    endif
    current_bond=>current_bond%next_bond
  enddo
end subroutine label_conglomerates

!> Keep only the bergs which are part of, or within contact distance of, a conglomerate that overlaps
!! the computational domain, or which are halo_berg=1
subroutine mts_remove_unused_bergs(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: this,other_berg,kick_the_bucket
  integer :: i,nc_x,nc_y,is,ie,js,je
  integer :: grdi, grdj, grdi2, grdj2
  logical :: contact,radial_contact
  real :: lat1,lat2,lon1,lon2,dlon,dlat,dx_dlon,dy_dlat,lat_ref
  real :: r_dist,crit_dist,R1,R2,rdenom

  grd=>bergs%grd
  nc_x=bergs%contact_cells_lon; nc_y=bergs%contact_cells_lat
  if (nc_x.eq.1 .and. nc_y.eq.1) then
    radial_contact=.true. !contact based on berg radii can occur
    if (bergs%hexagonal_icebergs) then
      rdenom=1./(2.*sqrt(3.))
    else
      if (bergs%iceberg_bonds_on) then
        rdenom=1./4.
      else
        rdenom=1./pi
      endif
    endif
  else
    radial_contact=.false.
    crit_dist=bergs%contact_distance**2
  endif

  do i =1,4
    select case (i) !halo cells on each side of the current PE
    case (1)
      is=grd%isd+1; ie=grd%isc-1; js=grd%jsd+1; je=grd%jed
    case (2)
      is=grd%iec+1; ie=grd%ied;   js=grd%jsd+1; je=grd%jed
    case (3)
      is=grd%isc;   ie=grd%iec;   js=grd%jsd+1; je=grd%jsc-1
    case (4)
      is=grd%isc;   ie=grd%iec;   js=grd%jec+1; je=grd%jed
    end select

    ! Loop over the current PE's halo cells
    do grdj = js,je ;  do grdi = is,ie
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this))
        if (this%halo_berg<2 .or. this%conglom_id>0) then !keep true halo and congloms bergs
          this=>this%next
        else !if not a halo/conglom berg, then test if w/in contact dist to a conglom/halo berg
          contact=.false.
          lat1=this%lat; lon1=this%lon
          if (radial_contact) R1=sqrt(this%length*this%width*rdenom) !radius of current berg
          do grdj2=max(grdj-nc_y,grd%jsd+1),min(grdj+nc_y,grd%jed)
            do grdi2=max(grdi-nc_x,grd%isd+1),min(grdi+nc_x,grd%ied)
              other_berg=>bergs%list(grdi2,grdj2)%first
              do while (associated(other_berg))
                if (other_berg%conglom_id>0) then !test contact w/ conglom berg
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

          if (.not. contact) then  !remove if not w/in contact dist to a true halo/conglom berg
            kick_the_bucket=>this
            this=>this%next
            call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
          else
            this%halo_berg=10 !contact berg
            this=>this%next
          endif
        endif
      enddo
    enddo; enddo
  enddo
end subroutine mts_remove_unused_bergs

!> For MTS velocity verlet, the send and receive call for for transfer_mts_bergs
subroutine mts_send_and_receive(bergs, nbergs_to_send_e, nbergs_to_send_w, nbergs_to_send_n, nbergs_to_send_s)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer :: nbergs_to_send_e !< Number of bergs to send to the east
  integer :: nbergs_to_send_w !< Number of bergs to send to the west
  integer :: nbergs_to_send_n !< Number of bergs to send to the north
  integer :: nbergs_to_send_s !< Number of bergs to send to the south
  ! Local vars
  type(icebergs_gridded), pointer :: grd
  integer :: nbergs_rcvd_from_e, nbergs_rcvd_from_w, nbergs_rcvd_from_n, nbergs_rcvd_from_s
  integer :: i,stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

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
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!'
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
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!'
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

  call mpp_sync_self()

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
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!'
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
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!'
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

  call mpp_sync_self()
end subroutine mts_send_and_receive


!> Destroys all bergs in a list
subroutine delete_all_bergs_in_list(bergs, grdj, grdi)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer :: grdi !< i-index of list
  integer :: grdj !< j-index of list
  ! Local variables
  type(iceberg), pointer :: kick_the_bucket, this
  this=>bergs%list(grdi,grdj)%first
  do while (associated(this))
    kick_the_bucket=>this
    this=>this%next
    call destroy_iceberg(kick_the_bucket)
   !call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
  enddo
  bergs%list(grdi,grdj)%first=>null()
end  subroutine delete_all_bergs_in_list


!> Send bergs in halo lists to other processors
subroutine send_bergs_to_other_pes(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
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
          if (.not. bergs%debug_write) call move_trajectory(bergs, kick_the_bucket)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
        elseif (this%ine.lt.bergs%grd%isc) then
          kick_the_bucket=>this
          this=>this%next
          nbergs_to_send_w=nbergs_to_send_w+1
          call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_w, nbergs_to_send_w, bergs%max_bonds  )
          if (.not. bergs%debug_write) call move_trajectory(bergs, kick_the_bucket)
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
          if (.not. bergs%debug_write) call move_trajectory(bergs, kick_the_bucket)
          call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,kick_the_bucket)
        elseif (this%jne.lt.bergs%grd%jsc) then
          kick_the_bucket=>this
          this=>this%next
          nbergs_to_send_s=nbergs_to_send_s+1
          call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_s, nbergs_to_send_s,bergs%max_bonds)
          if (.not. bergs%debug_write) call move_trajectory(bergs, kick_the_bucket)
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
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_end=',nbergs_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_start=',nbergs_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: error=',nbergs_end-(nbergs_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_to_send_n=',nbergs_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_to_send_s=',nbergs_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_to_send_e=',nbergs_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_to_send_w=',nbergs_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_rcvd_from_n=',nbergs_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_rcvd_from_s=',nbergs_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_rcvd_from_e=',nbergs_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'KID, send_bergs_to_other_pes: nbergs_rcvd_from_w=',nbergs_rcvd_from_w,' on PE',mpp_pe()
      call error_mesg('KID, send_bergs_to_other_pes:', 'We lost some bergs!', FATAL)
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
      write(stderrunit,'(a,i4)') 'KID, send_bergs_to_other_pes: # of bergs outside computational domain = ',i
      call error_mesg('KID, send_bergs_to_other_pes:', 'there are bergs still in halos!', FATAL)
    endif ! root_pe
  endif ! debug

  call mpp_sync_self()

end subroutine send_bergs_to_other_pes

!> Pack a berg into a buffer
subroutine pack_berg_into_buffer2(berg, buff, n, max_bonds_in)
! Arguments
type(iceberg), pointer :: berg !< Iceberg to pack into buffer
type(buffer), pointer :: buff !< Buffer to pack berg into
integer, intent(in) :: n !< Position in buffer to place berg
integer, optional :: max_bonds_in !< undocumented
!integer, intent(in) :: max_bonds  ! Change this later
! Local variables
integer :: counter, k, max_bonds, id_cnt, id_ij
type(bond), pointer :: current_bond

  max_bonds=0
  if (present(max_bonds_in)) max_bonds=max_bonds_in

  if (.not.associated(buff)) call increase_ibuffer(buff,n,buffer_width)
  if (n>buff%size) call increase_ibuffer(buff,n,buffer_width)

  counter = 0
  call push_buffer_value(buff%data(:,n), counter, berg%lon)
  call push_buffer_value(buff%data(:,n), counter, berg%lat)
  call push_buffer_value(buff%data(:,n), counter, berg%lon_prev)
  call push_buffer_value(buff%data(:,n), counter, berg%lat_prev)
  call push_buffer_value(buff%data(:,n), counter, berg%uvel)
  call push_buffer_value(buff%data(:,n), counter, berg%vvel)
  call push_buffer_value(buff%data(:,n), counter, berg%uvel_prev)
  call push_buffer_value(buff%data(:,n), counter, berg%vvel_prev)
  call push_buffer_value(buff%data(:,n), counter, berg%xi)
  call push_buffer_value(buff%data(:,n), counter, berg%yj)
  call push_buffer_value(buff%data(:,n), counter, berg%start_lon)
  call push_buffer_value(buff%data(:,n), counter, berg%start_lat)
  call push_buffer_value(buff%data(:,n), counter, berg%start_year)
  call push_buffer_value(buff%data(:,n), counter, berg%start_day)
  call push_buffer_value(buff%data(:,n), counter, berg%start_mass)
  call push_buffer_value(buff%data(:,n), counter, berg%mass)
  call push_buffer_value(buff%data(:,n), counter, berg%thickness)
  call push_buffer_value(buff%data(:,n), counter, berg%width)
  call push_buffer_value(buff%data(:,n), counter, berg%length)
  call push_buffer_value(buff%data(:,n), counter, berg%fl_k)
  call push_buffer_value(buff%data(:,n), counter, berg%mass_scaling)
  call push_buffer_value(buff%data(:,n), counter, berg%mass_of_bits)
  call push_buffer_value(buff%data(:,n), counter, berg%mass_of_fl_bits)
  call push_buffer_value(buff%data(:,n), counter, berg%mass_of_fl_bergy_bits)
  call push_buffer_value(buff%data(:,n), counter, berg%heat_density)
  call push_buffer_value(buff%data(:,n), counter, berg%ine)
  call push_buffer_value(buff%data(:,n), counter, berg%jne)
  call push_buffer_value(buff%data(:,n), counter, berg%axn)
  call push_buffer_value(buff%data(:,n), counter, berg%ayn)
  call push_buffer_value(buff%data(:,n), counter, berg%bxn)
  call push_buffer_value(buff%data(:,n), counter, berg%byn)
  call push_buffer_value(buff%data(:,n), counter, berg%halo_berg)
  call push_buffer_value(buff%data(:,n), counter, berg%static_berg)
  call split_id(berg%id, id_cnt, id_ij)
  call push_buffer_value(buff%data(:,n), counter, id_cnt)
  call push_buffer_value(buff%data(:,n), counter, id_ij)
  call push_buffer_value(buff%data(:,n), counter, berg%od)

  if (mts) then
    call push_buffer_value(buff%data(:,n), counter, berg%conglom_id)
    call push_buffer_value(buff%data(:,n), counter, berg%axn_fast)
    call push_buffer_value(buff%data(:,n), counter, berg%ayn_fast)
    call push_buffer_value(buff%data(:,n), counter, berg%bxn_fast)
    call push_buffer_value(buff%data(:,n), counter, berg%byn_fast)
    call push_buffer_value(buff%data(:,n), counter, berg%uo)
    call push_buffer_value(buff%data(:,n), counter, berg%vo)
    call push_buffer_value(buff%data(:,n), counter, berg%ua)
    call push_buffer_value(buff%data(:,n), counter, berg%va)
    call push_buffer_value(buff%data(:,n), counter, berg%ui)
    call push_buffer_value(buff%data(:,n), counter, berg%vi)
    call push_buffer_value(buff%data(:,n), counter, berg%ssh_x)
    call push_buffer_value(buff%data(:,n), counter, berg%ssh_y)
    call push_buffer_value(buff%data(:,n), counter, berg%sst)
    call push_buffer_value(buff%data(:,n), counter, berg%sss)
    call push_buffer_value(buff%data(:,n), counter, berg%cn)
    call push_buffer_value(buff%data(:,n), counter, berg%hi)
  endif

  if (iceberg_bonds_on) then
    call push_buffer_value(buff%data(:,n), counter, berg%n_bonds)
  end if

  if (monitor_energy) then
    call push_buffer_value(buff%data(:,n), counter, berg%Ee)
    call push_buffer_value(buff%data(:,n), counter, berg%Ed)
    call push_buffer_value(buff%data(:,n), counter, berg%Eext)
    call push_buffer_value(buff%data(:,n), counter, berg%Ee_contact)
    call push_buffer_value(buff%data(:,n), counter, berg%Ed_contact)
    call push_buffer_value(buff%data(:,n), counter, berg%Efrac)
    call push_buffer_value(buff%data(:,n), counter, berg%Ee_contact_temp)
    call push_buffer_value(buff%data(:,n), counter, berg%Ed_contact_temp)
    call push_buffer_value(buff%data(:,n), counter, berg%Ee_temp)
    call push_buffer_value(buff%data(:,n), counter, berg%Ed_temp)
    call push_buffer_value(buff%data(:,n), counter, berg%Eext_temp)
  endif

  if (dem) then
    call push_buffer_value(buff%data(:,n), counter, berg%ang_vel)
    call push_buffer_value(buff%data(:,n), counter, berg%ang_accel)
    call push_buffer_value(buff%data(:,n), counter, berg%rot)
  elseif (fracture_criterion .ne. 'none') then
    call push_buffer_value(buff%data(:,n), counter, berg%accum_bond_rotation)
  endif

  if (max_bonds .gt. 0) then
    current_bond=>berg%first_bond
    do k = 1,max_bonds
      if (associated(current_bond)) then
        call split_id(current_bond%other_id, id_cnt, id_ij)
        call push_buffer_value(buff%data(:,n), counter, id_cnt)
        call push_buffer_value(buff%data(:,n), counter, id_ij)
        call push_buffer_value(buff%data(:,n), counter, current_bond%other_berg_ine)
        call push_buffer_value(buff%data(:,n), counter, current_bond%other_berg_jne)
        call push_buffer_value(buff%data(:,n), counter, current_bond%length)

        if (use_damage) then
          call push_buffer_value(buff%data(:,n), counter, current_bond%damage)
        endif

        if (monitor_energy) then
          call push_buffer_value(buff%data(:,n), counter, current_bond%Ee)
          call push_buffer_value(buff%data(:,n), counter, current_bond%Ed)
          call push_buffer_value(buff%data(:,n), counter, current_bond%axn_fast)
          call push_buffer_value(buff%data(:,n), counter, current_bond%ayn_fast)
          call push_buffer_value(buff%data(:,n), counter, current_bond%bxn_fast)
          call push_buffer_value(buff%data(:,n), counter, current_bond%byn_fast)
        endif

        if (dem) then
          call push_buffer_value(buff%data(:,n), counter, current_bond%tangd1)
          call push_buffer_value(buff%data(:,n), counter, current_bond%tangd2)
          call push_buffer_value(buff%data(:,n), counter, current_bond%nstress)
          call push_buffer_value(buff%data(:,n), counter, current_bond%sstress)
        elseif (fracture_criterion .ne. 'none') then
          call push_buffer_value(buff%data(:,n), counter, current_bond%rotation)
          call push_buffer_value(buff%data(:,n), counter, current_bond%rel_rotation)
          call push_buffer_value(buff%data(:,n), counter, current_bond%n_frac_var)
          if (fracture_criterion .eq. 'strain_rate') then
            call push_buffer_value(buff%data(:,n), counter, current_bond%n_strain_rate)
          endif
          if (fracture_criterion .eq. 'energy') then
            call push_buffer_value(buff%data(:,n), counter, current_bond%spring_pe)
          endif
        endif

        current_bond=>current_bond%next_bond
      else
        call push_buffer_value(buff%data(:,n), counter, 0)
        call push_buffer_value(buff%data(:,n), counter, 0)
        call push_buffer_value(buff%data(:,n), counter, 0)
        call push_buffer_value(buff%data(:,n), counter, 0)
        call push_buffer_value(buff%data(:,n), counter, 0)

        if (use_damage) then
          call push_buffer_value(buff%data(:,n), counter, 0)
        endif

        if (monitor_energy) then
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
        endif

        if (dem) then
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
        elseif (fracture_criterion .ne. 'none') then
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          call push_buffer_value(buff%data(:,n), counter, 0)
          if (fracture_criterion .eq. 'strain_rate') then
            call push_buffer_value(buff%data(:,n), counter, 0)
          endif
          if (fracture_criterion .eq. 'energy') then
            call push_buffer_value(buff%data(:,n), counter, 0)
          endif
        endif

      endif
    enddo
  endif

  ! Clearing berg pointer from partner bonds
  !if (berg%halo_berg .lt. 0.5) then
  !  call clear_berg_from_partners_bonds(berg)
  !endif

end subroutine pack_berg_into_buffer2

!> Set a real value in the buffer at position (counter) after incrementing counter
subroutine push_buffer_rvalue(vbuf, counter, val)
  real, dimension(:), intent(inout) :: vbuf    !< Buffer vector to push value into
  integer,            intent(inout) :: counter !< Position to increment
  real,               intent(in)    :: val     !< Value to place in buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in push_buffer_rvalue'
  vbuf(counter) = val

end subroutine push_buffer_rvalue

!> Set an integer value in the buffer at position (counter) after incrementing counter
subroutine push_buffer_ivalue(vbuf, counter, val)
  real, dimension(:), intent(inout) :: vbuf    !< Buffer vector to push value into
  integer,            intent(inout) :: counter !< Position to increment
  integer,            intent(in)    :: val     !< Value to place in buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in push_buffer_ivalue'
  vbuf(counter) = float(val)

end subroutine push_buffer_ivalue

!> Get a real value from the buffer at position (counter) after incrementing counter
subroutine pull_buffer_rvalue(vbuf, counter, val)
  real, dimension(:), intent(in)    :: vbuf    !< Buffer vector to pull value from
  integer,            intent(inout) :: counter !< Position to increment
  real,               intent(out)   :: val     !< Value to get from buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in pull_buffer_rvalue'
  val = vbuf(counter)

end subroutine pull_buffer_rvalue

!> Get an integer value from the buffer at position (counter) after incrementing counter
subroutine pull_buffer_ivalue(vbuf, counter, val)
  real, dimension(:), intent(in)    :: vbuf    !< Buffer vector to pull value from
  integer,            intent(inout) :: counter !< Position to increment
  integer,            intent(out)   :: val     !< Value to get from buffer

  counter = counter + 1
  if (counter > size(vbuf)) stop 'OOB in pull_buffer_ivalue'
  val = nint(vbuf(counter))

end subroutine pull_buffer_ivalue

!> Nullify the other_berg pointer for a bond pointing to the current berg
subroutine clear_berg_from_partners_bonds(berg)
! Arguments
type(iceberg), intent(in), pointer :: berg !< Iceberg to clear from partner bonds
! Local variables
type(iceberg), pointer :: other_berg
type(bond), pointer :: current_bond, matching_bond
integer ::  stderrunit
! Get the stderr unit number
stderrunit = stderr()

  current_bond=>berg%first_bond
  do while (associated(current_bond)) !Looping over bonds
    other_berg=>current_bond%other_berg
    if (associated(other_berg)) then
      !write(stderrunit,*) , 'Other berg', berg%id, other_berg%id, mpp_pe()
      matching_bond=>other_berg%first_bond
      do while (associated(matching_bond))  ! Looping over possible matching bonds in other_berg
        if (matching_bond%other_id .eq. berg%id) then
          !write(stderrunit,*) , 'Clearing', berg%id, matching_bond%other_id,other_berg%id, mpp_pe()
          matching_bond%other_berg=>null()
          matching_bond=>null()
          if (iceberg_bonds_on) then !should not find bonds unless they are on anyway, but just to be safe...
            if (other_berg%n_bonds>0) other_berg%n_bonds=other_berg%n_bonds-1
          endif
        else
          matching_bond=>matching_bond%next_bond
        endif
      enddo
    else
     ! Note: This is meant to be unmatched after you have cleared the first berg
     ! call error_mesg('KID, clear berg from partners', 'The bond you are trying to clear is unmatched!', WARNING)
    endif
    current_bond=>current_bond%next_bond
  enddo !End loop over bonds

end subroutine clear_berg_from_partners_bonds

!> Unpacks a berg entry from a buffer to a new berg
subroutine unpack_berg_from_buffer2(bergs, buff, n, grd, force_append, max_bonds_in)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
type(buffer), pointer :: buff !< Buffer from which to unpack berg
integer, intent(in) :: n !< Position in buffer to unpack
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
logical, optional :: force_append !< undocumented
integer, optional :: max_bonds_in !< undocumented
! Local variables
!real :: lon, lat, uvel, vvel, xi, yj
!real :: start_lon, start_lat, start_day, start_mass
!integer :: ine, jne, start_year
logical :: lres
type(iceberg) :: localberg
type(iceberg), pointer :: this,kick_the_bucket
type(bond) , pointer :: current_bond
integer :: other_berg_ine, other_berg_jne
integer :: counter, k, max_bonds, id_cnt, id_ij
integer(kind=8) :: id
integer :: stderrunit
logical :: force_app,duplicate
real :: temp_lon,temp_lat,length

  ! Get the stderr unit number
  stderrunit = stderr()

  max_bonds=0
  if (present(max_bonds_in)) max_bonds=max_bonds_in

  force_app = .false.
  if(present(force_append)) force_app = force_append

  if (mts) then
    allocate(localberg%axn_fast,localberg%ayn_fast,localberg%bxn_fast,localberg%byn_fast,localberg%conglom_id)
  endif

  if (iceberg_bonds_on) then
    allocate(localberg%n_bonds)
  endif

  if (monitor_energy) then
    allocate(localberg%Ee,localberg%Ed,localberg%Eext,localberg%Ee_contact,localberg%Ed_contact,localberg%Efrac,&
      localberg%Ee_temp,localberg%Ed_temp,localberg%Eext_temp,localberg%Ee_contact_temp,localberg%Ed_contact_temp)
  endif

  if (dem) then
    allocate(localberg%ang_vel,localberg%ang_accel,localberg%rot)
  elseif (fracture_criterion .ne. 'none') then
    allocate(localberg%accum_bond_rotation)
  endif

  counter = 0
  call pull_buffer_value(buff%data(:,n), counter, localberg%lon)
  call pull_buffer_value(buff%data(:,n), counter, localberg%lat)
  call pull_buffer_value(buff%data(:,n), counter, localberg%lon_prev)
  call pull_buffer_value(buff%data(:,n), counter, localberg%lat_prev)
  call pull_buffer_value(buff%data(:,n), counter, localberg%uvel)
  call pull_buffer_value(buff%data(:,n), counter, localberg%vvel)
  call pull_buffer_value(buff%data(:,n), counter, localberg%uvel_prev)
  call pull_buffer_value(buff%data(:,n), counter, localberg%vvel_prev)
  call pull_buffer_value(buff%data(:,n), counter, localberg%xi)
  call pull_buffer_value(buff%data(:,n), counter, localberg%yj)
  call pull_buffer_value(buff%data(:,n), counter, localberg%start_lon)
  call pull_buffer_value(buff%data(:,n), counter, localberg%start_lat)
  call pull_buffer_value(buff%data(:,n), counter, localberg%start_year)
  call pull_buffer_value(buff%data(:,n), counter, localberg%start_day)
  call pull_buffer_value(buff%data(:,n), counter, localberg%start_mass)
  call pull_buffer_value(buff%data(:,n), counter, localberg%mass)
  call pull_buffer_value(buff%data(:,n), counter, localberg%thickness)
  call pull_buffer_value(buff%data(:,n), counter, localberg%width)
  call pull_buffer_value(buff%data(:,n), counter, localberg%length)
  call pull_buffer_value(buff%data(:,n), counter, localberg%fl_k)
  call pull_buffer_value(buff%data(:,n), counter, localberg%mass_scaling)
  call pull_buffer_value(buff%data(:,n), counter, localberg%mass_of_bits)
  call pull_buffer_value(buff%data(:,n), counter, localberg%mass_of_fl_bits)
  call pull_buffer_value(buff%data(:,n), counter, localberg%mass_of_fl_bergy_bits)
  call pull_buffer_value(buff%data(:,n), counter, localberg%heat_density)
  call pull_buffer_value(buff%data(:,n), counter, localberg%ine)
  call pull_buffer_value(buff%data(:,n), counter, localberg%jne)
  call pull_buffer_value(buff%data(:,n), counter, localberg%axn)
  call pull_buffer_value(buff%data(:,n), counter, localberg%ayn)
  call pull_buffer_value(buff%data(:,n), counter, localberg%bxn)
  call pull_buffer_value(buff%data(:,n), counter, localberg%byn)
  call pull_buffer_value(buff%data(:,n), counter, localberg%halo_berg)
  call pull_buffer_value(buff%data(:,n), counter, localberg%static_berg)
  call pull_buffer_value(buff%data(:,n), counter, id_cnt)
  call pull_buffer_value(buff%data(:,n), counter, id_ij)
  localberg%id = id_from_2_ints(id_cnt, id_ij)
  call pull_buffer_value(buff%data(:,n), counter, localberg%od)

  if (bergs%mts) then
    call pull_buffer_value(buff%data(:,n), counter, localberg%conglom_id)
    call pull_buffer_value(buff%data(:,n), counter, localberg%axn_fast)
    call pull_buffer_value(buff%data(:,n), counter, localberg%ayn_fast)
    call pull_buffer_value(buff%data(:,n), counter, localberg%bxn_fast)
    call pull_buffer_value(buff%data(:,n), counter, localberg%byn_fast)
    call pull_buffer_value(buff%data(:,n), counter, localberg%uo)
    call pull_buffer_value(buff%data(:,n), counter, localberg%vo)
    call pull_buffer_value(buff%data(:,n), counter, localberg%ua)
    call pull_buffer_value(buff%data(:,n), counter, localberg%va)
    call pull_buffer_value(buff%data(:,n), counter, localberg%ui)
    call pull_buffer_value(buff%data(:,n), counter, localberg%vi)
    call pull_buffer_value(buff%data(:,n), counter, localberg%ssh_x)
    call pull_buffer_value(buff%data(:,n), counter, localberg%ssh_y)
    call pull_buffer_value(buff%data(:,n), counter, localberg%sst)
    call pull_buffer_value(buff%data(:,n), counter, localberg%sss)
    call pull_buffer_value(buff%data(:,n), counter, localberg%cn)
    call pull_buffer_value(buff%data(:,n), counter, localberg%hi)
  endif

  if (iceberg_bonds_on) then
    call pull_buffer_value(buff%data(:,n), counter, localberg%n_bonds)
  endif

  if (monitor_energy) then
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ee)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ed)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Eext)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ee_contact)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ed_contact)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Efrac)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ee_contact_temp)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ed_contact_temp)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ee_temp)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Ed_temp)
    call pull_buffer_value(buff%data(:,n), counter, localberg%Eext_temp)
  endif

  if (dem) then
    call pull_buffer_value(buff%data(:,n), counter, localberg%ang_vel)
    call pull_buffer_value(buff%data(:,n), counter, localberg%ang_accel)
    call pull_buffer_value(buff%data(:,n), counter, localberg%rot)
  elseif (fracture_criterion .ne. 'none') then
    call pull_buffer_value(buff%data(:,n), counter, localberg%accum_bond_rotation)
  endif

  !These quantities no longer need to be passed between processors
  localberg%uvel_old=localberg%uvel
  localberg%vvel_old=localberg%vvel
  localberg%lon_old=localberg%lon
  localberg%lat_old=localberg%lat

  !when using MTS Verlet, iceberg elements that are part of a conglomerate may lie outside
  !of the PE's grid. They are assigned to the berg list of the nearest halo cell.
  !The global (lon/lat) and local (xi,yi) coords of these elements are not adjusted.
  if (bergs%mts .and. abs(localberg%halo_berg)>1) then
    temp_lon = localberg%lon; temp_lat = localberg%lat

    if (temp_lon > grd%lon(grd%ied,grd%jed)) then
      temp_lon = grd%lonc(grd%ied,grd%jed)
    elseif (temp_lon < grd%lon(grd%isd,grd%jsd)) then
      temp_lon = grd%lonc(grd%isd+1,grd%jsd+1)
    endif
    if (temp_lat > grd%lat(grd%ied,grd%jed)) then
      temp_lat = grd%latc(grd%ied,grd%jed)
    elseif (temp_lat < grd%lat(grd%isd,grd%jsd)) then
      temp_lat = grd%latc(grd%isd+1,grd%jsd+1)
    endif

    lres=find_cell_wide(grd, temp_lon, temp_lat, localberg%ine, localberg%jne)
    if (lres) then
      call unpack_duplicate_check(bergs,localberg,duplicate)
      if (duplicate) return
      call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg,this)
    else
      write(stderrunit,'("KID, unpack_berg_from_buffer pe=(",i3,a,2i4,a,2f8.2)')&
        & mpp_pe(),') Failed to find i,j=',localberg%ine,localberg%jne,&
        ' for mts conglom berg lon,lat=',temp_lon, temp_lat
      write(stderrunit,*) temp_lon, temp_lat
      write(stderrunit,*) localberg%uvel,localberg%vvel
      write(stderrunit,*) localberg%axn,localberg%ayn !Alon
      write(stderrunit,*) localberg%bxn,localberg%byn !Alon
      if (bergs%mts) then
        write(stderrunit,*) localberg%axn_fast,localberg%ayn_fast !Alex
        write(stderrunit,*) localberg%bxn_fast,localberg%byn_fast !Alex
      endif
      write(stderrunit,*) localberg%uvel_old,localberg%vvel_old
      write(stderrunit,*) localberg%lon_old,localberg%lat_old
      write(stderrunit,*) grd%isc,grd%iec,grd%jsc,grd%jec
      write(stderrunit,*) grd%isd,grd%ied,grd%jsd,grd%jed
      write(stderrunit,*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
      write(stderrunit,*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
      write(stderrunit,*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
      write(stderrunit,*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
      write(stderrunit,*) temp_lat,temp_lon
      write(stderrunit,*) lres
      call error_mesg('KID, unpack_berg_from_buffer', 'can not find a cell to place berg in!', FATAL)
    endif
  else
    ! force_app=.true.
    if(force_app) then !force append with origin ine,jne (for I/O)
      call unpack_duplicate_check(bergs,localberg,duplicate)
      if (duplicate) return
      call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg, this)
    else
      ! lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      !this version is faster, by checking if current ine and jne are appropriate before searching rest of cells.
      lres=check_and_find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      if (lres) then
        call unpack_duplicate_check(bergs,localberg,duplicate)
        if (duplicate) return
        lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
        call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg,this)
      else
        lres=find_cell_wide(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
        if (lres) then
          call unpack_duplicate_check(bergs,localberg,duplicate)
          if (duplicate) return
          lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
          call add_new_berg_to_list(bergs%list(localberg%ine,localberg%jne)%first, localberg,this)
        else
          write(stderrunit,'("KID, unpack_berg_from_buffer pe=(",i3,a,2i4,a,2f8.2)')&
            & mpp_pe(),') Failed to find i,j=',localberg%ine,localberg%jne,' for lon,lat=',localberg%lon,localberg%lat
          write(stderrunit,*) localberg%halo_berg
          write(stderrunit,*) localberg%lon,localberg%lat
          write(stderrunit,*) localberg%uvel,localberg%vvel
          write(stderrunit,*) localberg%axn,localberg%ayn !Alon
          write(stderrunit,*) localberg%bxn,localberg%byn !Alon
          if (bergs%mts) then
            write(stderrunit,*) localberg%axn_fast,localberg%ayn_fast !Alex
            write(stderrunit,*) localberg%bxn_fast,localberg%byn_fast !Alex
          endif
          write(stderrunit,*) localberg%uvel_old,localberg%vvel_old
          write(stderrunit,*) localberg%lon_old,localberg%lat_old
          write(stderrunit,*) grd%isc,grd%iec,grd%jsc,grd%jec
          write(stderrunit,*) grd%isd,grd%ied,grd%jsd,grd%jed
          write(stderrunit,*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
          write(stderrunit,*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
          write(stderrunit,*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
          write(stderrunit,*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
          write(stderrunit,*) lres
          call error_mesg('KID, unpack_berg_from_buffer', 'can not find a cell to place berg in!', FATAL)
        endif
      endif
    endif
  endif

  !#  Do stuff to do with bonds here MP1
  this%first_bond=>null()
  if (max_bonds .gt. 0) then
    do k = 1,max_bonds
      call pull_buffer_value(buff%data(:,n), counter, id_cnt)
      call pull_buffer_value(buff%data(:,n), counter, id_ij)
      id = id_from_2_ints(id_cnt, id_ij)
      call pull_buffer_value(buff%data(:,n), counter, other_berg_ine)
      call pull_buffer_value(buff%data(:,n), counter, other_berg_jne)
      call pull_buffer_value(buff%data(:,n), counter, length)

      if (id_ij > 0) then
        call form_a_bond(this, id, other_berg_ine, other_berg_jne)
        current_bond=>this%first_bond
        current_bond%length=length

        if (use_damage) then
          call pull_buffer_value(buff%data(:,n), counter, current_bond%damage)
        endif

        if (monitor_energy) then
          call pull_buffer_value(buff%data(:,n), counter, current_bond%Ee)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%Ed)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%axn_fast)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%ayn_fast)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%bxn_fast)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%byn_fast)
        endif

        if (dem) then
          call pull_buffer_value(buff%data(:,n), counter, current_bond%tangd1)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%tangd2)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%nstress)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%sstress)
        elseif (fracture_criterion .ne. 'none') then
          call pull_buffer_value(buff%data(:,n), counter, current_bond%rotation)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%rel_rotation)
          call pull_buffer_value(buff%data(:,n), counter, current_bond%n_frac_var)
          if (fracture_criterion .eq. 'strain_rate') then
            call pull_buffer_value(buff%data(:,n), counter, current_bond%n_strain_rate)
          endif
          if (fracture_criterion .eq. 'energy') then
            call pull_buffer_value(buff%data(:,n), counter, current_bond%spring_pe)
          endif
        endif
      endif
    enddo
  endif
  this=>null()

end subroutine unpack_berg_from_buffer2

!> Increase size of buffer
!!
!! This routine checks if the buffer size is smaller than nbergs
!! If it is, the buffer size is increased by delta_buf.
!! The buffer increases by more than 1 so that the buffer does not have to increase every time.
subroutine increase_ibuffer(old, num_bergs, width)
! Arguments
type(buffer), pointer :: old !< Buffer to expand
integer, intent(in) :: num_bergs !< Number of bergs
integer, intent(in) :: width !< Width of buffer (first dimension)
! Local variables
type(buffer), pointer :: new
integer :: new_size, old_size

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
   !write(stderr(),*) 'KID, increase_ibuffer',mpp_pe(),' increased to',new_size
  endif

end subroutine increase_ibuffer

!> Packs a trajectory entry into a buffer
subroutine pack_traj_into_buffer2(traj, buff, n, save_short_traj, save_fl_traj, fl_r)
  ! Arguments
  type(xyt), pointer :: traj !< Trajectory entry to pack
  type(buffer), pointer :: buff !< Buffer to pack entry into
  integer, intent(in) :: n !< Position in buffer to place entry
  logical, intent(in) :: save_short_traj !< If true, only use a subset of trajectory data
  logical, intent(in) :: save_fl_traj !< If true, save masses and footloose parameters
  real, intent(in) :: fl_r !< If >0 and save_fl_traj, save footloose parameters
  ! Local variables
  integer :: counter ! Position in stack
  integer :: cnt, ij

  if (.not.associated(buff)) call increase_ibuffer(buff,n,buffer_width_traj)
  if (n>buff%size) call increase_ibuffer(buff,n,buffer_width_traj)

  counter = 0
  call push_buffer_value(buff%data(:,n),counter,traj%lon)
  call push_buffer_value(buff%data(:,n),counter,traj%lat)
  call push_buffer_value(buff%data(:,n),counter,traj%year)
  call push_buffer_value(buff%data(:,n),counter,traj%day)
  call split_id(traj%id, cnt, ij)
  call push_buffer_value(buff%data(:,n),counter,cnt)
  call push_buffer_value(buff%data(:,n),counter,ij)
  if (save_fl_traj) then
    call push_buffer_value(buff%data(:,n),counter,traj%mass)
    call push_buffer_value(buff%data(:,n),counter,traj%start_mass)
    call push_buffer_value(buff%data(:,n),counter,traj%thickness)
    call push_buffer_value(buff%data(:,n),counter,traj%mass_of_bits)
    call push_buffer_value(buff%data(:,n),counter,traj%uvel)
    call push_buffer_value(buff%data(:,n),counter,traj%vvel)
    if (fl_r>0) then
      call push_buffer_value(buff%data(:,n),counter,traj%mass_scaling)
      call push_buffer_value(buff%data(:,n),counter,traj%mass_of_fl_bits)
      call push_buffer_value(buff%data(:,n),counter,traj%mass_of_fl_bergy_bits)
      call push_buffer_value(buff%data(:,n),counter,traj%fl_k)
    endif
  endif
  if (.not. save_short_traj) then
    !call push_buffer_value(buff%data(:,n),counter,traj%uvel)
    !call push_buffer_value(buff%data(:,n),counter,traj%vvel)
    call push_buffer_value(buff%data(:,n),counter,traj%uvel_prev)
    call push_buffer_value(buff%data(:,n),counter,traj%vvel_prev)
    call push_buffer_value(buff%data(:,n),counter,traj%heat_density)
    call push_buffer_value(buff%data(:,n),counter,traj%width)
    call push_buffer_value(buff%data(:,n),counter,traj%length)
    call push_buffer_value(buff%data(:,n),counter,traj%uo)
    call push_buffer_value(buff%data(:,n),counter,traj%vo)
    call push_buffer_value(buff%data(:,n),counter,traj%ui)
    call push_buffer_value(buff%data(:,n),counter,traj%vi)
    call push_buffer_value(buff%data(:,n),counter,traj%ua)
    call push_buffer_value(buff%data(:,n),counter,traj%va)
    call push_buffer_value(buff%data(:,n),counter,traj%ssh_x)
    call push_buffer_value(buff%data(:,n),counter,traj%ssh_y)
    call push_buffer_value(buff%data(:,n),counter,traj%sst)
    call push_buffer_value(buff%data(:,n),counter,traj%sss)
    call push_buffer_value(buff%data(:,n),counter,traj%cn)
    call push_buffer_value(buff%data(:,n),counter,traj%hi)
    call push_buffer_value(buff%data(:,n),counter,traj%axn)
    call push_buffer_value(buff%data(:,n),counter,traj%ayn)
    call push_buffer_value(buff%data(:,n),counter,traj%bxn)
    call push_buffer_value(buff%data(:,n),counter,traj%byn)
    call push_buffer_value(buff%data(:,n),counter,traj%halo_berg)
    call push_buffer_value(buff%data(:,n),counter,traj%static_berg)
    call push_buffer_value(buff%data(:,n),counter,traj%od)

    if (mts) then
      call push_buffer_value(buff%data(:,n), counter, traj%axn_fast)
      call push_buffer_value(buff%data(:,n), counter, traj%ayn_fast)
      call push_buffer_value(buff%data(:,n), counter, traj%bxn_fast)
      call push_buffer_value(buff%data(:,n), counter, traj%byn_fast)
    endif

    if (iceberg_bonds_on) then
      call push_buffer_value(buff%data(:,n), counter, traj%n_bonds)
    end if

    if (monitor_energy) then
      call push_buffer_value(buff%data(:,n), counter, traj%Ee)
      call push_buffer_value(buff%data(:,n), counter, traj%Ed)
      call push_buffer_value(buff%data(:,n), counter, traj%Eext)
      call push_buffer_value(buff%data(:,n), counter, traj%Ee_contact)
      call push_buffer_value(buff%data(:,n), counter, traj%Ed_contact)
      call push_buffer_value(buff%data(:,n), counter, traj%Efrac)
      ! call push_buffer_value(buff%data(:,n), counter, traj%Ee_contact_temp)
      ! call push_buffer_value(buff%data(:,n), counter, traj%Ed_contact_temp)
      ! call push_buffer_value(buff%data(:,n), counter, traj%Ee_temp)
      ! call push_buffer_value(buff%data(:,n), counter, traj%Ed_temp)
      ! call push_buffer_value(buff%data(:,n), counter, traj%Eext_temp)
    endif

    if (dem) then
      call push_buffer_value(buff%data(:,n), counter, traj%ang_vel)
      call push_buffer_value(buff%data(:,n), counter, traj%ang_accel)
      call push_buffer_value(buff%data(:,n), counter, traj%rot)
    elseif (fracture_criterion .ne. 'none') then
      call push_buffer_value(buff%data(:,n), counter, traj%accum_bond_rotation)
    endif

  endif

end subroutine pack_traj_into_buffer2

!> Unpacks a trajectory entry from a buffer
subroutine unpack_traj_from_buffer2(first, buff, n, save_short_traj, save_fl_traj, fl_r)
  ! Arguments
  type(xyt), pointer :: first !< Trajectory list
  type(buffer), pointer :: buff !< Buffer from which to unpack
  integer, intent(in) :: n !< Position in buffer to unpack
  logical, intent(in) :: save_short_traj !< If true, only use a subset of trajectory data
  logical, intent(in) :: save_fl_traj !< If true, save masses and footloose parameters
  real, intent(in) :: fl_r !< If >0 and save_fl_traj, save footloose parameters
  ! Local variables
  type(xyt) :: traj
  integer :: counter ! Position in stack
  integer :: cnt, ij

  if (mts) then
    allocate(traj%axn_fast,traj%ayn_fast,traj%bxn_fast,traj%byn_fast)
  endif

  if (iceberg_bonds_on) then
    allocate(traj%n_bonds)
  endif

  if (monitor_energy) then
    allocate(traj%Ee,traj%Ed,traj%Eext,traj%Ee_contact,traj%Ed_contact,traj%Efrac)
    !allocate(traj%Ee_temp,traj%Ed_temp,traj%Eext_temp,traj%Ee_contact_temp,traj%Ed_contact_temp)
  endif

  if (dem) then
    allocate(traj%ang_vel,traj%ang_accel,traj%rot)
  elseif (fracture_criterion .ne. 'none') then
    allocate(traj%accum_bond_rotation)
  endif

  counter = 0
  call pull_buffer_value(buff%data(:,n),counter,traj%lon)
  call pull_buffer_value(buff%data(:,n),counter,traj%lat)
  call pull_buffer_value(buff%data(:,n),counter,traj%year)
  call pull_buffer_value(buff%data(:,n),counter,traj%day)
  call pull_buffer_value(buff%data(:,n),counter,cnt)
  call pull_buffer_value(buff%data(:,n),counter,ij)
  traj%id = id_from_2_ints(cnt, ij)
  if (save_fl_traj) then
    call pull_buffer_value(buff%data(:,n),counter,traj%mass)
    call pull_buffer_value(buff%data(:,n),counter,traj%start_mass)
    call pull_buffer_value(buff%data(:,n),counter,traj%thickness)
    call pull_buffer_value(buff%data(:,n),counter,traj%mass_of_bits)
    call pull_buffer_value(buff%data(:,n),counter,traj%uvel)
    call pull_buffer_value(buff%data(:,n),counter,traj%vvel)
    if (fl_r>0) then
      call pull_buffer_value(buff%data(:,n),counter,traj%mass_scaling)
      call pull_buffer_value(buff%data(:,n),counter,traj%mass_of_fl_bits)
      call pull_buffer_value(buff%data(:,n),counter,traj%mass_of_fl_bergy_bits)
      call pull_buffer_value(buff%data(:,n),counter,traj%fl_k)
    endif
  endif
  if (.not. save_short_traj) then
    !call pull_buffer_value(buff%data(:,n),counter,traj%uvel)
    !call pull_buffer_value(buff%data(:,n),counter,traj%vvel)
    call pull_buffer_value(buff%data(:,n),counter,traj%uvel_prev)
    call pull_buffer_value(buff%data(:,n),counter,traj%vvel_prev)
    call pull_buffer_value(buff%data(:,n),counter,traj%heat_density)
    call pull_buffer_value(buff%data(:,n),counter,traj%width)
    call pull_buffer_value(buff%data(:,n),counter,traj%length)
    call pull_buffer_value(buff%data(:,n),counter,traj%uo)
    call pull_buffer_value(buff%data(:,n),counter,traj%vo)
    call pull_buffer_value(buff%data(:,n),counter,traj%ui)
    call pull_buffer_value(buff%data(:,n),counter,traj%vi)
    call pull_buffer_value(buff%data(:,n),counter,traj%ua)
    call pull_buffer_value(buff%data(:,n),counter,traj%va)
    call pull_buffer_value(buff%data(:,n),counter,traj%ssh_x)
    call pull_buffer_value(buff%data(:,n),counter,traj%ssh_y)
    call pull_buffer_value(buff%data(:,n),counter,traj%sst)
    call pull_buffer_value(buff%data(:,n),counter,traj%sss)
    call pull_buffer_value(buff%data(:,n),counter,traj%cn)
    call pull_buffer_value(buff%data(:,n),counter,traj%hi)
    call pull_buffer_value(buff%data(:,n),counter,traj%axn)
    call pull_buffer_value(buff%data(:,n),counter,traj%ayn)
    call pull_buffer_value(buff%data(:,n),counter,traj%bxn)
    call pull_buffer_value(buff%data(:,n),counter,traj%byn)
    call pull_buffer_value(buff%data(:,n),counter,traj%halo_berg)
    call pull_buffer_value(buff%data(:,n),counter,traj%static_berg)
    call pull_buffer_value(buff%data(:,n),counter,traj%od)

    if (mts) then
      call pull_buffer_value(buff%data(:,n), counter, traj%axn_fast)
      call pull_buffer_value(buff%data(:,n), counter, traj%ayn_fast)
      call pull_buffer_value(buff%data(:,n), counter, traj%bxn_fast)
      call pull_buffer_value(buff%data(:,n), counter, traj%byn_fast)
    endif

    if (iceberg_bonds_on) then
      call pull_buffer_value(buff%data(:,n), counter, traj%n_bonds)
    end if

    if (monitor_energy) then
      call pull_buffer_value(buff%data(:,n), counter, traj%Ee)
      call pull_buffer_value(buff%data(:,n), counter, traj%Ed)
      call pull_buffer_value(buff%data(:,n), counter, traj%Eext)
      call pull_buffer_value(buff%data(:,n), counter, traj%Ee_contact)
      call pull_buffer_value(buff%data(:,n), counter, traj%Ed_contact)
      call pull_buffer_value(buff%data(:,n), counter, traj%Efrac)
      ! call pull_buffer_value(buff%data(:,n), counter, traj%Ee_contact_temp)
      ! call pull_buffer_value(buff%data(:,n), counter, traj%Ed_contact_temp)
      ! call pull_buffer_value(buff%data(:,n), counter, traj%Ee_temp)
      ! call pull_buffer_value(buff%data(:,n), counter, traj%Ed_temp)
      ! call pull_buffer_value(buff%data(:,n), counter, traj%Eext_temp)
    endif

    if (dem) then
      call pull_buffer_value(buff%data(:,n), counter, traj%ang_vel)
      call pull_buffer_value(buff%data(:,n), counter, traj%ang_accel)
      call pull_buffer_value(buff%data(:,n), counter, traj%rot)
    elseif (fracture_criterion .ne. 'none') then
      call pull_buffer_value(buff%data(:,n), counter, traj%accum_bond_rotation)
    endif
  endif
  call append_posn(first, traj)

end subroutine unpack_traj_from_buffer2

!> Packs a trajectory entry into a buffer
subroutine pack_bond_traj_into_buffer2(bond_traj, buff, n)
  ! Arguments
  type(bond_xyt), pointer :: bond_traj !< Trajectory entry to pack
  type(buffer), pointer :: buff !< Buffer to pack entry into
  integer, intent(in) :: n !< Position in buffer to place entry
  ! Local variables
  integer :: counter ! Position in stack
  integer :: cnt1, ij1, cnt2, ij2

  if (.not.associated(buff)) call increase_ibuffer(buff,n,buffer_width_bond_traj)
  if (n>buff%size) call increase_ibuffer(buff,n,buffer_width_bond_traj)

  counter = 0
  call push_buffer_value(buff%data(:,n),counter,bond_traj%lon)
  call push_buffer_value(buff%data(:,n),counter,bond_traj%lat)
  call push_buffer_value(buff%data(:,n),counter,bond_traj%year)
  call push_buffer_value(buff%data(:,n),counter,bond_traj%day)
  call push_buffer_value(buff%data(:,n),counter,bond_traj%length)
  call push_buffer_value(buff%data(:,n),counter,bond_traj%n1)
  call push_buffer_value(buff%data(:,n),counter,bond_traj%n2)

  call split_id(bond_traj%id1, cnt1, ij1)
  call push_buffer_value(buff%data(:,n),counter,cnt1)
  call push_buffer_value(buff%data(:,n),counter,ij1)
  call split_id(bond_traj%id2, cnt2, ij2)
  call push_buffer_value(buff%data(:,n),counter,cnt2)
  call push_buffer_value(buff%data(:,n),counter,ij2)

  if (use_damage) then
    call push_buffer_value(buff%data(:,n),counter,bond_traj%damage)
  endif

  if (monitor_energy) then
    call push_buffer_value(buff%data(:,n),counter,bond_traj%Ee)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%Ed)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%axn_fast)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%ayn_fast)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%bxn_fast)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%byn_fast)
  endif

  if (dem) then
    call push_buffer_value(buff%data(:,n),counter,bond_traj%tangd1)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%tangd2)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%nstress)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%sstress)
  elseif (fracture_criterion.ne.'none') then
    call push_buffer_value(buff%data(:,n),counter,bond_traj%rotation)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%rel_rotation)
    call push_buffer_value(buff%data(:,n),counter,bond_traj%n_frac_var)
    if (fracture_criterion.eq.'strain_rate') then
      call push_buffer_value(buff%data(:,n),counter,bond_traj%n_strain_rate)
    elseif (fracture_criterion.eq.'energy') then
      call push_buffer_value(buff%data(:,n),counter,bond_traj%spring_pe)
    endif
  endif

end subroutine pack_bond_traj_into_buffer2

!> Unpacks a trajectory entry from a buffer
subroutine unpack_bond_traj_from_buffer2(first, buff, n)
  ! Arguments
  type(bond_xyt), pointer :: first !< Bond trajectory list
  type(buffer), pointer :: buff !< Buffer from which to unpack
  integer, intent(in) :: n !< Position in buffer to unpack
  ! Local variables
  type(bond_xyt) :: bond_traj
  integer :: counter ! Position in stack
  integer :: cnt1, ij1, cnt2, ij2

  if (use_damage) then
    allocate(bond_traj%damage)
  endif
  if (monitor_energy) then
    allocate(bond_traj%Ee,bond_traj%Ed,bond_traj%axn_fast,bond_traj%ayn_fast,bond_traj%bxn_fast,bond_traj%byn_fast)
  endif
  if (dem) then
    allocate(bond_traj%tangd1,bond_traj%tangd2,bond_traj%nstress,bond_traj%sstress)
  elseif (fracture_criterion.ne.'none') then
    allocate(bond_traj%rotation,bond_traj%rel_rotation,bond_traj%n_frac_var)
    if (fracture_criterion.eq.'strain_rate') then
      allocate(bond_traj%n_strain_rate)
    elseif (fracture_criterion.eq.'energy') then
      allocate(bond_traj%spring_pe)
    endif
  endif

  counter = 0
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%lon)
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%lat)
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%year)
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%day)
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%length)
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%n1)
  call pull_buffer_value(buff%data(:,n),counter,bond_traj%n2)
  call pull_buffer_value(buff%data(:,n),counter,cnt1)
  call pull_buffer_value(buff%data(:,n),counter,ij1)
  bond_traj%id1 = id_from_2_ints(cnt1, ij1)
  call pull_buffer_value(buff%data(:,n),counter,cnt2)
  call pull_buffer_value(buff%data(:,n),counter,ij2)
  bond_traj%id2 = id_from_2_ints(cnt2, ij2)

  if (use_damage) then
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%damage)
  endif

  if (monitor_energy) then
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%Ee)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%Ed)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%axn_fast)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%ayn_fast)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%bxn_fast)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%byn_fast)
  endif

  if (dem) then
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%tangd1)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%tangd2)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%nstress)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%sstress)
  elseif (fracture_criterion.ne.'none') then
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%rotation)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%rel_rotation)
    call pull_buffer_value(buff%data(:,n),counter,bond_traj%n_frac_var)
    if (fracture_criterion.eq.'strain_rate') then
      call pull_buffer_value(buff%data(:,n),counter,bond_traj%n_strain_rate)
    elseif (fracture_criterion.eq.'energy') then
      call pull_buffer_value(buff%data(:,n),counter,bond_traj%spring_pe)
    endif
  endif

  call append_bond_posn(first, bond_traj)

end subroutine unpack_bond_traj_from_buffer2

!> Add a new berg to a list by copying values
!!
!! The input berg are a berg with set values whose memory is assumed to be
!! temporary. This routine allocates memory for a new berg and copies the
!! the input values into it. The memory for the new berg is pointed to
!! by newberg_return (if present).
subroutine add_new_berg_to_list(first, bergvals, newberg_return)
! Arguments
type(iceberg), pointer :: first !< List of icebergs
type(iceberg), intent(in) :: bergvals !< Berg values to copy
type(iceberg), intent(out), pointer, optional :: newberg_return !< New berg
! Local variables
type(iceberg), pointer :: new=>null()

  new=>null()
  call create_iceberg(new, bergvals)

  if (present(newberg_return)) then
    newberg_return=>new
    !newberg_return=>null()
  endif

  call insert_berg_into_list(first, new)

  !Clear new
  new=>null()

end subroutine add_new_berg_to_list


!> Scans all lists and checks that bergs are in a sorted order in each list
subroutine count_out_of_order(bergs,label)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
character(len=*) :: label !< Label to add to messages
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
        call error_mesg('KID, count_out_of_order', 'Pointer %prev is unassociated. This should not happen!', FATAL)
      endif
    enddo
  enddo; enddo
  call mpp_sum(icnt2)

  if ((debug.or.icnt1.ne.0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,3(x,a,i6),x,a)') 'KID, count_out_of_order:', &
      '# out of order=', icnt1,'# in halo=',icnt2,'# identicals=',icnt3,label
  endif

  call check_for_duplicates(bergs,label)

end subroutine count_out_of_order

!> Scans a single cell for duplicate bergs
subroutine unpack_duplicate_check(bergs,localberg,duplicate)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg) :: localberg !< Input iceberg
  logical :: duplicate !< True if a duplicate is detected
  ! Local variables
  type(iceberg), pointer :: this ! bergs for comparison

  duplicate=.false.
  if (.not. bergs%mts) return
  this=>bergs%list(localberg%ine,localberg%jne)%first
  do while (associated(this) .and. .not. duplicate)
    if (this%id.eq.localberg%id) then
      if (localberg%halo_berg==10) then
        this%halo_berg=10
      elseif (localberg%halo_berg<this%halo_berg) then
        this%halo_berg=localberg%halo_berg
      endif
      duplicate=.true.
    endif
    this=>this%next
  end do
  this=>null()
end subroutine unpack_duplicate_check

!> Scans all lists and checks for duplicate identifiers between lists
subroutine check_for_duplicates(bergs,label)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
character(len=*) :: label !< Label to add to message
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
    write(*,'(a,2(x,a,i9),x,a)') 'KID, check_for_duplicates:', &
      '# with same id=', icnt_id,'# identical bergs=',icnt_same,label
  endif

end subroutine check_for_duplicates

!> Generate an iceberg id from a counter at calving-location and the calving location itself.
!! Note that this updates grd%iceberg_counter_grd.
!!
!! \todo If we initialized grd%iceberg_counter_grd to 0 when the model is first run we could move
!! this increment line to before the id generation and then the counter would be an actual count.
integer(kind=8) function generate_id(grd, i, j)
  type(icebergs_gridded), pointer    :: grd !< Container for gridded fields
  integer,                intent(in) :: i   !< i-index of calving location
  integer,                intent(in) :: j   !< j-index of calving location
  ! Local variables
  integer :: ij ! Hash of i,j

  ! Increment counter in calving cell
  grd%iceberg_counter_grd(i,j) = grd%iceberg_counter_grd(i,j) + 1
  ! ij is unique number for each grid cell (32-bit integers allow for ~1/100th degree global resolution)
  ij = ij_component_of_id(grd, i, j)
  ! Generate a 64-bit id
  generate_id = id_from_2_ints( grd%iceberg_counter_grd(i,j), ij )

end function generate_id

!> Convert an old 32-bit id to a 64-bit id
integer(kind=8) function convert_old_id(grd, old_id)
  type(icebergs_gridded), pointer    :: grd    !< Container for gridded fields
  integer,                intent(in) :: old_id !< 32-bit iceberg id
  ! Local variables
  integer :: cnt ! Counter component
  integer :: ij ! Hash of i,j
  integer :: i,j ! Cell indexes

  call cij_from_old_id(grd, old_id, cnt, i, j)
  ij = ij_component_of_id(grd, i, j)
  convert_old_id = id_from_2_ints( cnt, ij )

end function convert_old_id

!> Recover i,j an old 32-bit id
subroutine cij_from_old_id(grd, old_id, cnt, i, j)
  type(icebergs_gridded), pointer     :: grd    !< Container for gridded fields
  integer,                intent(in)  :: old_id !< 32-bit iceberg id
  integer,                intent(out) :: cnt    !< Counter component of old id
  integer,                intent(out) :: i      !< i-index of calving cell
  integer,                intent(out) :: j      !< j-index of calving cell
  ! Local variables
  integer :: ij ! Hash of i,j
  integer :: iNg, jNg, ncells ! Shape and size of the global grid

  ! Number cells in the grid
  iNg = grd%ieg - grd%isg + 1
  jNg = grd%jeg - grd%jsg + 1
  ncells = iNg * jNg

  ! cnt is the cell-based counter
  cnt = old_id / ncells
  ! ij is the has of i,j
  ij = mod( old_id, ncells )

  ! i,j should be the cells the berg was calved in
  j = ij / iNg
  i = mod( ij, iNg )

end subroutine cij_from_old_id

!> Calculate the location-derived component of an iceberg id which is a hash of the i,j-indexes for the cell
integer function ij_component_of_id(grd, i, j)
  type(icebergs_gridded), pointer    :: grd !< Container for gridded fields
  integer,                intent(in) :: i   !< i-index of calving location
  integer,                intent(in) :: j   !< j-index of calving location
  ! Local variables
  integer :: ij ! Hash of i,j
  integer :: iNg ! Zonal size of the global grid

  ! Using the current grid shape maximizes the numbers of IDs that can be represented
  ! allowing up to 30-minute uniform global resolution, or potentially finer if non-uniform.
  iNg = grd%ieg - grd%isg + 1

  ! ij_component_of_id is unique number for each grid cell (32-bit integers allow for ~1/100th degree global resolution)
  ij_component_of_id = i + ( iNg * ( j - 1 ) )

end function ij_component_of_id

!> Prints a particular berg's vitals
!!
!! All lists are scanned and if a berg has the identifier equal to
!! debug_iceberg_with_id then the state of that berg is printed.
subroutine monitor_a_berg(bergs, label)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
character(len=*) :: label !< Label to add to message
! Local variables
type(iceberg), pointer :: this
integer :: grdi, grdj
integer :: stderrunit

  if (bergs%debug_iceberg_with_id<0) return
  stderrunit=stderr() ! Get the stderr unit number

  do grdj = bergs%grd%jsd,bergs%grd%jed ; do grdi = bergs%grd%isd,bergs%grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      if (this%id == bergs%debug_iceberg_with_id) then
        call print_berg(stderrunit, this, 'MONITOR: '//label, grdi, grdj)
      endif
      this=>this%next
    enddo
  enddo ; enddo

end subroutine monitor_a_berg

!> Inserts a berg into a list
subroutine insert_berg_into_list(first, newberg)
! Arguments
type(iceberg), pointer :: first !< List of bergs
type(iceberg), pointer :: newberg !< New berg to insert
! Local variables
type(iceberg), pointer :: this, prev

  if (associated(first)) then
    if (.not. parallel_reprod) then
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

!> Returns True when berg1 and berg2 are in sorted order
!! \todo inorder() should use the iceberg identifier for efficiency and simplicity
!! instead of dates and properties
logical function inorder(berg1, berg2)  !MP Alon - Change to use iceberg id
! Arguments
type(iceberg), pointer :: berg1 !< An iceberg
type(iceberg), pointer :: berg2 !< An iceberg
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


!> Returns a hash of a berg's start year and day
!! \todo Should be able to remove this function if using identifiers properly
real function time_hash(berg)!  Alon: Think about removing this.
! Arguments
type(iceberg), pointer :: berg !< An iceberg
  time_hash=berg%start_day+366.*float(berg%start_year)
end function time_hash

!> Returns a hash of a berg's start position
!! \todo Should be able to remove this function if using identifiers properly
real function pos_hash(berg)
! Arguments
type(iceberg), pointer :: berg !< An iceberg
  pos_hash=berg%start_lon+360.*(berg%start_lat+90.)
end function pos_hash

!> Returns True if berg1 and berg2 have the identifying properties
!!
!! This function compares the start year, day, mass and position of the bergs.
logical function sameid(berg1, berg2) !  Alon: MP updat this.
! Arguments
type(iceberg), pointer :: berg1 !< An iceberg
type(iceberg), pointer :: berg2 !< An iceberg
! Local variables
  sameid=.false.
  if (berg1%start_year.ne.berg2%start_year) return
  if (berg1%start_day.ne.berg2%start_day) return
  if (berg1%start_mass.ne.berg2%start_mass) return
  if (berg1%start_lon.ne.berg2%start_lon) return
  if (berg1%start_lat.ne.berg2%start_lat) return
  sameid=.true. ! passing the above tests means that bergs 1 and 2 have the same id
end function sameid

!> Returns True if berg1 and berg2 are identical in both identifying properties
!! and dynamic properties
logical function sameberg(berg1, berg2)
! Arguments
type(iceberg), pointer :: berg1 !< An iceberg
type(iceberg), pointer :: berg2 !< An iceberg
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

  if (mts) then
    if (berg1%axn_fast.ne.berg2%axn_fast) return  !Alex
    if (berg1%ayn_fast.ne.berg2%ayn_fast) return  !Alex
    if (berg1%bxn_fast.ne.berg2%bxn_fast) return  !Alex
    if (berg1%byn_fast.ne.berg2%byn_fast) return  !Alex
  endif
  sameberg=.true. ! passing the above tests mean that bergs 1 and 2 are identical
end function sameberg

!> Returns the year day (a single float for the day of the year, range 0-365.999...)
real function yearday(imon, iday, ihr, imin, isec)
! Arguments
integer, intent(in) :: imon !< Month of year (1-12)
integer, intent(in) :: iday !< Day of month (1-31)
integer, intent(in) :: ihr !< Hour of day (0-23)
integer, intent(in) :: imin !< Minute of hour (0-59)
integer, intent(in) :: isec !< Second of minute (0-59)

  yearday=float(imon-1)*31.+float(iday-1)+(float(ihr)+(float(imin)+float(isec)/60.)/60.)/24.

end function yearday

!> Create a new berg with given values
subroutine create_iceberg(berg, bergvals)
! Arguments
type(iceberg), pointer :: berg !< Berg to be created
type(iceberg), intent(in) :: bergvals !< Values to assign
! Local variables
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  if (associated(berg)) then
    write(stderrunit,*) 'KID, create_iceberg: berg already associated!!!!',mpp_pe()
    call error_mesg('KID, create_iceberg', 'berg already associated. This should not happen!', FATAL)
  endif
  allocate(berg)

  if (mts) then
    allocate(berg%axn_fast,berg%ayn_fast,berg%bxn_fast,berg%byn_fast,berg%conglom_id)
  endif

  if (iceberg_bonds_on) then
    allocate(berg%n_bonds)
    berg%n_bonds=0
  endif

  if (monitor_energy) then
    allocate(berg%Ee,berg%Ed,berg%Eext,berg%Ee_contact,berg%Ed_contact,berg%Efrac,&
      berg%Ee_temp,berg%Ed_temp,berg%Eext_temp,berg%Ee_contact_temp,berg%Ed_contact_temp)
  endif
  if (dem) then
    allocate(berg%ang_vel,berg%ang_accel,berg%rot)
  elseif (fracture_criterion .ne. 'none') then
    allocate(berg%accum_bond_rotation)
  endif
  berg=bergvals
  berg%prev=>null()
  berg%next=>null()

end subroutine create_iceberg


!> Delete a berg from a list and destroy the memory for the berg
!!
!! first is needed when berg is the first in the list
subroutine delete_iceberg_from_list(first, berg)
! Arguments
type(iceberg), pointer :: first !< List of bergs
type(iceberg), pointer :: berg !< Berg to be deleted
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

!> Destroy a berg
!!
!! Deallocates memory for a berg after deleting links to the berg recorded in bonds
subroutine destroy_iceberg(berg)
! Arguments
type(iceberg), pointer :: berg
! Local variables

  ! Clears all matching bonds before deallocating memory
  call clear_berg_from_partners_bonds(berg)

  ! Bye-bye berg
  deallocate(berg)

end subroutine destroy_iceberg

!> Print the state of a particular berg
subroutine print_berg(iochan, berg, label, il, jl)
! Arguments
integer, intent(in) :: iochan !< Standard channel to use (usually stdout or stderr)
type(iceberg), pointer :: berg !< Berg to print
character(len=*) :: label !< Label to use in messages
integer, optional, intent(in) :: il !< i-index of cell berg should be in
integer, optional, intent(in) :: jl !< j-index of cell berg should be in
! Local variables

  write(iochan,'("KID, print_berg: ",2a,i5,a,i12,a,2f16.4,i5,f7.2,es12.4,f5.1)') &
    label, 'pe=(', mpp_pe(), ') #=', berg%id, ' start lon,lat,yr,day,mass,hb=', &
    berg%start_lon, berg%start_lat, berg%start_year, berg%start_day, berg%start_mass, berg%halo_berg
  if (present(il).and.present(jl)) then
    write(iochan,'("KID, print_berg: ",2a,i5,a,i12,a,2i5)') &
      label, 'pe=(', mpp_pe(), ') #=', berg%id, ' List i,j=',il,jl
  endif
  write(iochan,'("KID, print_berg: ",2a,i5,a,i12,a,2i5,a,2l2)') &
    label, 'pe=(', mpp_pe(), ') #=', berg%id, &
    ' i,j=', berg%ine, berg%jne, &
    ' p,n=', associated(berg%prev), associated(berg%next)
  write(iochan,'("KID, print_berg: ",2a,i5,a,i12,3(a,2f16.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%id, &
    ' xi,yj=', berg%xi, berg%yj, &
    ' lon,lat=', berg%lon, berg%lat, &
    ' lon_old,lat_old=', berg%lon_old, berg%lat_old
  write(iochan,'("KID, print_berg: ",2a,i5,a,i12,2(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%id, &
    ' u,v=', berg%uvel, berg%vvel, &
    ' uvel_old,vvel_old=', berg%uvel_old, berg%vvel_old
  write(iochan,'("KID, print_berg: ",2a,i5,a,i12,2(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%id, &
    ' axn,ayn=', berg%axn, berg%ayn, &
    ' bxn,byn=', berg%bxn, berg%byn
  if (mts) then
    write(iochan,'("KID, print_berg: ",2a,i5,a,i12,2(a,2f14.8))') &
      label, 'pe=(', mpp_pe(), ') #=', berg%id, &
      ' axn_fast,ayn_fast=', berg%axn_fast, berg%ayn_fast, &
      ' bxn_fast,byn_fast=', berg%bxn_fast, berg%byn_fast
  endif
  write(iochan,'("KID, print_berg: ",2a,i5,a,i12,3(a,2f14.8))') &
    label, 'pe=(', mpp_pe(), ') #=', berg%id, &
    ' uo,vo=', berg%uo, berg%vo, &
    ' ua,va=', berg%ua, berg%va, &
    ' ui,vi=', berg%ui, berg%vi
end subroutine print_berg

!> Print the state of all bergs
subroutine print_bergs(iochan, bergs, label)
! Arguments
integer, intent(in) :: iochan !< Standard channel to use (usually stdout or stderr)
type(icebergs), pointer :: bergs !< Container for all types and memory
character(len=*) :: label !< Label to use in messages
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
  if (nbergs.gt.0) write(iochan,'("KID, ",a," there are",i5," bergs out of",i6," on PE ",i4)') label, nbergs, nnbergs, mpp_pe()

end subroutine print_bergs

!> Set initial length of the bond
subroutine orig_bond_length(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: other_berg, berg
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
type(bond) , pointer :: current_bond,prev,next
real :: dist

  grd=>bergs%grd

  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        other_berg=>current_bond%other_berg
        if (associated(current_bond%other_berg)) then
          dist=(berg%lon-other_berg%lon)**2+(berg%lat-other_berg%lat)**2
          current_bond%length=sqrt(dist)
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo; enddo
end subroutine orig_bond_length

!> Initialize fracture parameters for non-DEM style fracture
subroutine fracture_testing_initialization(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this
type(bond) , pointer :: current_bond
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj

  grd=>bergs%grd
  do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      this%accum_bond_rotation=0.0
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        current_bond%rotation=0.0
        current_bond%rel_rotation=0.0
        current_bond%n_frac_var=0.0
        current_bond%n_strain_rate=0.0
        current_bond=>current_bond%next_bond
      enddo
      this=>this%next
    enddo
  enddo;enddo
end subroutine fracture_testing_initialization

!< Sets initial damage for bonds that overlap y=17000.0
subroutine damage_test_1_init(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: berg, other_berg
type(bond) , pointer :: current_bond
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj

  print *,'DAMAGE TEST 1 INITIALIZATION'
  grd=>bergs%grd
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      current_bond=>berg%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        other_berg=>current_bond%other_berg
        if (associated(current_bond%other_berg)) then
          if ((other_berg%lat>17000.0 .and. berg%lat<=17000.0) .or. &
            (berg%lat>17000.0 .and. other_berg%lat<=17000.0)) then
            current_bond%damage=0.75
          endif
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo; enddo
end subroutine damage_test_1_init

!> Save number of bonds on element
subroutine assign_n_bonds(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this, other_berg
type(bond) , pointer :: current_bond
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj

  grd=>bergs%grd
  do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      this%n_bonds=0
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        this%n_bonds=this%n_bonds+1
        current_bond=>current_bond%next_bond
      enddo
      this=>this%next
    enddo
  enddo;enddo
end subroutine assign_n_bonds

!> For if using constant_LW_for_interactions: if constant_length and constant_width
!! are not given in the nml or are both set to zero, determine their values here
subroutine set_constant_interaction_length_and_width(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
real :: elem_sum, l_sum, w_sum

  ! constant_length/constant_width are set to the mean values of the lengths/widths
  ! that are currently initialized on the domain. Note that using the mean should be
  ! overkill because ideally, the initialized elements should all have the same
  ! length and width anyway...

  if (.not. bergs%constant_interaction_LW) return

  elem_sum = 0. !counter for number of elements
  l_sum = 0.    !sum of lengths
  w_sum = 0.    !sum of widths
  grd=>bergs%grd
  do grdj=grd%jsc,grd%jec ; do grdi=grd%isc,grd%iec!loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      elem_sum=elem_sum+1.
      l_sum=l_sum+this%length
      w_sum=w_sum+this%width
      this=>this%next
    enddo
  enddo;enddo

  call mpp_sum(elem_sum); call mpp_sum(l_sum); call mpp_sum(w_sum)
  bergs%constant_length = l_sum/elem_sum
  bergs%constant_width  = w_sum/elem_sum

end subroutine set_constant_interaction_length_and_width

!> Initialization for the DEM beam tests
subroutine dem_tests_init(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
real :: minlon,maxlon

  maxlon=-huge(1.0); minlon=huge(1.0)
  grd=>bergs%grd
  do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      this%start_lon=this%lon; this%start_lat=this%lat
      if (this%lon>maxlon) maxlon=this%lon
      if (this%lon<minlon) minlon=this%lon
      this=>this%next
    enddo
  enddo;enddo

  call mpp_max(maxlon); call mpp_min(minlon)
  bergs%dem_tests_start_lon=minlon; bergs%dem_tests_end_lon=maxlon

  print *,'pe,dem tests start lon,dem tests end lon',mpp_pe(),minlon,maxlon

  if (bergs%dem_beam_test>=3) then
    do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
      this=>bergs%list(grdi,grdj)%first
      do while (associated(this)) ! loop over all bergs in cell
        if (bergs%dem_beam_test==3) then
          if (this%start_lon==minlon) this%rot=-pi/2.
        elseif (bergs%dem_beam_test==4) then
          if (this%start_lon==maxlon) then
            this%vvel=-40.0
            this%vvel_old=-40.0
            this%vvel_prev=-40.0
          endif
        endif
        this=>this%next
      enddo
    enddo;enddo
  endif

end subroutine dem_tests_init

!> Initialize berg element parameters for DEM-style simulations
subroutine init_dem_params(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
real :: maxlon

  grd=>bergs%grd
  do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      this%ang_vel=0.
      this%ang_accel=0.
      this%rot=0.
      this=>this%next
    enddo
  enddo;enddo
end subroutine init_dem_params

!> Set the accumulated bond rotation to zero for all bergs
subroutine reset_bond_rotation(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this, other_berg
type(bond) , pointer :: current_bond
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj

if (.not. bergs%fracture_criterion=='strain_rate') return

  grd=>bergs%grd
  do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      this%accum_bond_rotation=0.0
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        current_bond%rotation=0.0
        current_bond=>current_bond%next_bond
      enddo
      this=>this%next
    enddo
  enddo;enddo
end subroutine reset_bond_rotation

!>Saves rotation of bond on each bond, and accumulates the rotation of all bonds connected to a particle
!!onto the particle.
subroutine update_bond_angles(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: this, other_berg
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
type(bond) , pointer :: current_bond
real :: dx_dlon1,dx_dlon2,dx_dlon,dy_dlat,lat_ref,b1_u,b1_v,b2_u,b2_v
real :: dxa,dya,dxb,dyb,vr1,vr2
real :: angular_momentum,moment_of_inertia,tangent_velocity
real :: cross_b1b2,dot_b1b2,angle

  if (bergs%fracture_criterion=='none') return

  grd=>bergs%grd

  if (.not. bergs%grd%grid_is_latlon) then
    dx_dlon=1.0; dy_dlat=1.0
  else
    dy_dlat=pi_180*Rearth
  endif

  do grdj=grd%jsd,grd%jed ; do grdi=grd%isd,grd%ied !loop over all cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs in cell
      if (this%conglom_id>0) then
        current_bond=>this%first_bond
        if (associated(current_bond)) then

          do while (associated(current_bond)) ! loop over all bonds
            other_berg=>current_bond%other_berg
            if (associated(other_berg)) then

              !relative coords from previous cycle
              if (bergs%grd%grid_is_latlon) then
                lat_ref=0.5*(this%lat_prev+other_berg%lat_prev); dx_dlon=pi_180*Rearth*cos(lat_ref*pi_180)
              endif
              dxa=(other_berg%lon_prev-this%lon_prev)*dx_dlon;   dya=(other_berg%lat_prev-this%lat_prev)*dy_dlat

              !relative coords from current cycle
              if (bergs%grd%grid_is_latlon) then
                lat_ref=0.5*(this%lat+other_berg%lat);         dx_dlon=pi_180*Rearth*cos(lat_ref*pi_180)
              endif
              dxb=(other_berg%lon-this%lon)*dx_dlon;           dyb=(other_berg%lat-this%lat)*dy_dlat

              current_bond%length=sqrt(dxb*dxb+dyb*dyb) !update bond length

              cross_b1b2=dxa*dyb-dya*dxb !([dxa dya 0] x [dxb dyb 0]) . [0 0 1]
              dot_b1b2 = dxa*dxb+dya*dyb ![dxa dya] . [dxa dxb]
              angle = atan2(cross_b1b2,dot_b1b2)

              current_bond%rotation=current_bond%rotation+angle
            endif
            current_bond=>current_bond%next_bond
          enddo !loop over all bonds
        endif
      endif
      this=>this%next
    enddo ! loop over all bergs in cell
  enddo; enddo !loop over all cells
end subroutine update_bond_angles


!> For each berg element, calculates the average rotation rate of all bonds connected to a particle. Bonds are
!! then broken based on their individual rotation rate compared to the average rotation rates accumulated on
!! the berg elements that the bonds connect, and based on the rate of bond stretching (normal strain-rate).
!! Alternatively, bonds can be broken based on the total accumulated relative rotation and normal strain (not rates).
subroutine update_and_break_bonds(bergs)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: other_berg, this
  type(icebergs_gridded), pointer :: grd
  integer :: grdi, grdj
  type(bond) , pointer :: current_bond,other_bond,kick_the_bucket
  real :: frac_thres_n,frac_thres_t,spring_coef
  real :: inv_dt, rdenom, mean_denom, R1, R2, n_frac_var, BondEetot

  frac_thres_n=bergs%frac_thres_n; frac_thres_t=bergs%frac_thres_t
  if (frac_thres_n<=0.0 .and. frac_thres_t<=0.0) return !all fracture ignored
  if (frac_thres_n<=0.0) frac_thres_n=huge(1.0) !no fracture from normal strain
  if (frac_thres_t<=0.0) frac_thres_t=huge(1.0) !no fracture from bond rotation

  if (bergs%fracture_criterion=='strain_rate') then
    inv_dt=1.0/bergs%dt
  elseif (bergs%fracture_criterion=='strain' .or. bergs%fracture_criterion=='stress') then
    inv_dt=1.0
  else
    return
  endif

  if (bergs%hexagonal_icebergs) then
    rdenom=1./(2.*sqrt(3.))
  else
    if (bergs%iceberg_bonds_on) then
      rdenom=1./4.
    else
      rdenom=1./pi
    endif
  endif

  grd=>bergs%grd

  !1. accumulate bond rotation/rotation-rate to the berg element
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs
      this%accum_bond_rotation=0.0
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        current_bond%rotation=current_bond%rotation*inv_dt
        this%accum_bond_rotation=this%accum_bond_rotation+current_bond%rotation
        current_bond=>current_bond%next_bond
      enddo
      this=>this%next
    enddo
  enddo;enddo

  !2. Compare rotation/rotation-rate of each bond to the berg-based mean w/o the bond it in
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied !loop over cells
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs
      if (this%conglom_id>0 .and. this%n_bonds>1) then
        R1=sqrt(this%length*this%width*rdenom)
        if (this%n_bonds>1) mean_denom=1.0/(this%n_bonds-1)

        !2. Bond updates: normal strain/strain-rate, and relative rotation/rotation-rate
        current_bond=>this%first_bond
        do while (associated(current_bond)) ! loop over all bonds
          other_berg=>current_bond%other_berg
          if (.not. associated(other_berg)) &
            call error_mesg('KID,update_and_break_bonds','other_berg not associated!' ,FATAL)
          R2=sqrt(other_berg%length*other_berg%width*rdenom)

          !the relative rotation (or rotation-rate) of the current bond relative to the mean
          !rotation (or rotation-rate) of bonds on the parent berg, where the current bond is
          !excluded from the mean.
          if (this%n_bonds>1) then
            current_bond%rel_rotation = current_bond%rotation-&
              (this%accum_bond_rotation-current_bond%rotation)*mean_denom
          else
            current_bond%rel_rotation = 0.0
          endif

          !normal strain or strain-rate: note original length is taken as the sum of the
          !radii (calculated according to area) of the two bergs the bond connects.
          if (bergs%fracture_criterion=='strain_rate') then
            n_frac_var=current_bond%length/(R1+R2) - 1.0
            current_bond%n_strain_rate=(n_frac_var-current_bond%n_frac_var)*inv_dt
            current_bond%n_frac_var=n_frac_var
            n_frac_var=current_bond%n_strain_rate
          elseif (bergs%fracture_criterion=='strain') then
            n_frac_var=current_bond%length/(R1+R2) - 1.0
            current_bond%n_frac_var=n_frac_var
          elseif (bergs%fracture_criterion=='stress') then
            n_frac_var=(current_bond%length-(R1+R2))*min(this%mass,other_berg%mass)*&
              (bergs%spring_coef)/(min(this%thickness,other_berg%thickness)*(R1+R2))
            current_bond%rel_rotation = current_bond%rel_rotation*R1*min(this%mass,other_berg%mass)*&
              bergs%spring_coef/(min(this%thickness,other_berg%thickness)*(R1+R2))
            current_bond%n_frac_var=n_frac_var
          elseif (bergs%fracture_criterion=='energy') then
            current_bond%n_frac_var=current_bond%spring_pe
          endif

        !mark bond and matching bond on other_berg as fractured
        if (abs(current_bond%rel_rotation)>frac_thres_t.or.n_frac_var>frac_thres_n) then
          current_bond%other_id=-1

          !mark matching bond on other_berg
          other_bond=>other_berg%first_bond
          do while (associated(other_bond)) ! loop over all bonds
            if (other_bond%other_id.eq.this%id) other_bond%other_id=-1
            other_bond=>other_bond%next_bond
          enddo
        endif
        current_bond=>current_bond%next_bond
      enddo
    endif
    this=>this%next
  enddo
enddo; enddo

  !4. final cleanup
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        if (current_bond%other_id==-1) then
          if (grdj>=grd%jsc .and. grdj<=grd%jec .and. grdi>=grd%isc .and. grdi<=grd%iec) then
            if (bergs%fracture_criterion=='strain_rate') then
              n_frac_var=current_bond%n_strain_rate
            else
              n_frac_var=current_bond%n_frac_var
            endif
            if (abs(current_bond%rel_rotation)>frac_thres_t.and.n_frac_var>frac_thres_n) then
              print *,'angle and n_frac_var break',current_bond%rel_rotation,n_frac_var,this%id
            elseif (abs(current_bond%rel_rotation)>frac_thres_t) then
              print *,'angle break',current_bond%rel_rotation,n_frac_var,this%id
            elseif (n_frac_var>frac_thres_n) then
              print *,'n_frac_var break',current_bond%rel_rotation,n_frac_var,this%id
            else
              print *,'break due to bond match',current_bond%rel_rotation,n_frac_var,this%id
            endif
          endif
          this%accum_bond_rotation=this%accum_bond_rotation-current_bond%rotation
          this%n_bonds=this%n_bonds-1
          if (monitor_energy) then
            !convert the elastic spring energy associated with the bond to dissipated energy
            !on the parent element.
            this%Ee=this%Ee-current_bond%Ee
            this%Ee_temp=this%Ee_temp-current_bond%Ee
            this%Efrac=this%Efrac+current_bond%Ee
          endif
          kick_the_bucket=>current_bond
          current_bond=>current_bond%next_bond
          call delete_bond_from_list(this,kick_the_bucket)
        else
          current_bond=>current_bond%next_bond
        endif
      enddo
      this=>this%next
    enddo
  enddo;enddo

end subroutine update_and_break_bonds

!> DEM-style fracture
subroutine break_bonds_dem(bergs)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(iceberg), pointer :: other_berg, this
  type(icebergs_gridded), pointer :: grd
  integer :: grdi, grdj
  type(bond) , pointer :: current_bond,other_bond,kick_the_bucket
  real :: frac_thres_n,frac_thres_t

  frac_thres_n=bergs%frac_thres_n; frac_thres_t=bergs%frac_thres_t

  if (frac_thres_n<=0.0 .and. frac_thres_t<=0.0) return !all fracture ignored
  if (frac_thres_n<=0.0) frac_thres_n=huge(1.0) !no fracture from normal strain
  if (frac_thres_t<=0.0) frac_thres_t=huge(1.0) !no fracture from bond rotation
  grd=>bergs%grd
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds

        if (use_damage) then
          frac_thres_n=bergs%frac_thres_n*(1.-current_bond%damage)
          frac_thres_t=bergs%frac_thres_t*(1.-current_bond%damage)
        endif

        if (current_bond%other_id.ne.-1) then

          if (current_bond%nstress>frac_thres_n .or. current_bond%sstress>frac_thres_t) then
            if (current_bond%nstress>frac_thres_n .and. current_bond%sstress>frac_thres_t) then
              print *,'normal and shear stress break',this%id,current_bond%nstress,current_bond%sstress
            elseif (current_bond%nstress>frac_thres_n) then
              print *,'normal stress break',this%id,current_bond%nstress,current_bond%sstress
            else
              print *,'shear stress break',this%id,current_bond%nstress,current_bond%sstress
            endif
            current_bond%other_id=-1
          endif

          if (current_bond%other_id.eq.-1) then
            other_berg=>current_bond%other_berg
            if (associated(other_berg)) then
              !mark matching bond on other_berg
              other_bond=>other_berg%first_bond
              do while (associated(other_bond)) ! loop over all bonds
                if (other_bond%other_id.eq.this%id) other_bond%other_id=-1
                other_bond=>other_bond%next_bond
              enddo
            endif
          endif
        endif
        current_bond=>current_bond%next_bond
      enddo
      this=>this%next
    enddo
  enddo;enddo


  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this)) ! loop over all bergs
      current_bond=>this%first_bond
      do while (associated(current_bond)) ! loop over all bonds
        if (current_bond%other_id.eq.-1) then
          this%n_bonds=this%n_bonds-1
          if (monitor_energy) then
            !convert the elastic spring energy associated with the bond to dissipated energy
            !on the parent element.
            this%Ee=this%Ee-current_bond%Ee
            this%Ee_temp=this%Ee_temp-current_bond%Ee
            this%Efrac=this%Efrac+current_bond%Ee
          endif
          kick_the_bucket=>current_bond
          current_bond=>current_bond%next_bond
          call delete_bond_from_list(this,kick_the_bucket)
        else
          current_bond=>current_bond%next_bond
        endif
      enddo
      this=>this%next
    enddo
  enddo;enddo
end subroutine break_bonds_dem

!> Delete current_bond from the list of its parent berg
subroutine delete_bond_from_list(berg,bond_to_delete)
type(iceberg), intent(in), pointer :: berg !<parent berg to bond_to_delete
type(bond), pointer :: bond_to_delete !<deleting this bond
type(bond), pointer :: prev,next,current_bond
  prev=>bond_to_delete%prev_bond
  next=>bond_to_delete%next_bond
  if (associated(prev)) then
    prev%next_bond=>next
  else
    berg%first_bond=>next
  endif
  if (associated(next)) next%prev_bond=>prev
  deallocate(bond_to_delete)
end subroutine delete_bond_from_list

!> Bond two bergs together
subroutine form_a_bond(berg, other_id, other_berg_ine, other_berg_jne, other_berg)
! Arguments
type(iceberg), pointer :: berg !< first berg
integer(kind=8), intent(in) :: other_id !< ID of the other other berg
integer, optional  :: other_berg_ine !< second berg zonal cell index
integer, optional :: other_berg_jne !< second berg meridional cell index
type(iceberg), optional,  pointer :: other_berg !< second berg
! Local variables
type(bond) , pointer :: new_bond, first_bond
integer :: stderrunit

 stderrunit = stderr()

if (berg%id .ne. other_id) then

 !write (stderrunit,*) , 'Forming a bond!!!', mpp_pe(), berg%id, other_id, berg%halo_berg, berg%ine, berg%jne

  ! Step 1: Create a new bond
  allocate(new_bond)
  if (use_damage) then
    allocate(new_bond%damage); new_bond%damage=0.
  endif
  if (monitor_energy) then
    allocate(new_bond%Ee,new_bond%Ed,new_bond%axn_fast,new_bond%ayn_fast,new_bond%bxn_fast,new_bond%byn_fast)
    new_bond%Ee=0.; new_bond%Ed=0.; new_bond%axn_fast=0.; new_bond%ayn_fast=0.; new_bond%bxn_fast=0.; new_bond%byn_fast=0.
  endif

  if (dem) then
    allocate(new_bond%tangd1,new_bond%tangd2)
    new_bond%tangd1=0.; new_bond%tangd2=0.
    allocate(new_bond%nstress,new_bond%sstress)
    new_bond%nstress=0.; new_bond%sstress=0.
  elseif (fracture_criterion.ne.'none') then
    allocate(new_bond%rotation,new_bond%rel_rotation,new_bond%n_frac_var)
    new_bond%rotation=0.; new_bond%rel_rotation=0.; new_bond%n_frac_var=0.
    if (fracture_criterion.eq.'strain_rate') then
      allocate(new_bond%n_strain_rate)
      new_bond%n_strain_rate=0.
    endif
    if (fracture_criterion.eq.'energy') then
      allocate(new_bond%spring_pe)
      new_bond%spring_pe=0.
    endif
  endif

  new_bond%other_id=other_id
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
   call error_mesg('KID, bonds', 'An iceberg is trying to bond with itself!!!', FATAL)
 endif

end subroutine form_a_bond

! #############################################################################

subroutine bond_address_update(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
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
            call error_mesg('KID, bond address update', 'other berg in bond not assosiated!', FATAL)
          endif
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo; enddo

  call mpp_sync_self()

end subroutine bond_address_update

subroutine show_all_bonds(bergs)
type(icebergs), pointer :: bergs !< Container for all types and memory
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
        !print *, 'Show Bond1 :', berg%id, current_bond%other_id, current_bond%other_berg_ine, current_bond%other_berg_jne,  mpp_pe()
        !write(*,'(a,5i)')'Show Bond1 :',berg%id, current_bond%other_id, current_bond%other_berg_ine, current_bond%other_berg_jne, mpp_pe()
        !print *, 'Current:', berg%id, berg%ine, berg%jne,berg%halo_berg, mpp_pe()
        if  (associated(current_bond%other_berg)) then
          if (current_bond%other_berg%id .ne. current_bond%other_id) then
            !print *, 'Bond matching', berg%id,current_bond%other_berg%id, current_bond%other_id,&
            !  berg%halo_berg,current_bond%other_berg%halo_berg ,mpp_pe()
            write(*,'(a,3i8,2i3,i4)')'Bond matching :',berg%id,current_bond%other_berg%id, current_bond%other_id,&
              int(berg%halo_berg),current_bond%other_berg%halo_berg ,mpp_pe()
            call error_mesg('KID, show all bonds:', 'The bonds are not matching properly!', FATAL)
          endif
        else
          ! print *, 'This bond has an non-assosiated other berg :', berg%id, current_bond%other_id,&
          !   current_bond%other_berg_ine, current_bond%other_berg_jne, berg%halo_berg,  mpp_pe()
          ! write(*,'(a,4i,i3,i4)')'This bond has an non-assosiated other berg :', berg%id, current_bond%other_id,&
          !   current_bond%other_berg_ine, current_bond%other_berg_jne, int(berg%halo_berg),  mpp_pe()
        endif
        current_bond=>current_bond%next_bond
      enddo
      berg=>berg%next
    enddo
  enddo; enddo

end subroutine show_all_bonds

!> Sweep across all bergs filling in bond data
subroutine connect_all_bonds(bergs, ignore_unmatched)
type(icebergs), pointer :: bergs !< Container for all types and memory
logical,optional :: ignore_unmatched !< If true, do not call error if the other berg in a bond is not found
! Local variables
type(iceberg), pointer :: other_berg, berg
type(icebergs_gridded), pointer :: grd
integer :: i, j
integer :: grdi, grdj
integer :: grdi_inner, grdj_inner
type(bond) , pointer :: current_bond, other_berg_bond
logical :: bond_matched, missing_bond, check_bond_quality,check_match
integer nbonds

  check_match = .true.
  if (present(ignore_unmatched)) then
    if (ignore_unmatched) check_match = .false.
  end if

missing_bond=.false.
bond_matched=.false.

! For convenience
  grd=>bergs%grd

  !If interacting bergs & domain is periodic, then correct lat & lon for halo bergs - Alex
  call update_latlon(bergs)

  do grdj = grd%jsd+1,grd%jed ; do grdi = grd%isd+1,grd%ied
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
              if (other_berg%id .eq. current_bond%other_id) then
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
                    if (other_berg%id .eq. current_bond%other_id) then
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
            ! Finally, if still not matched, search adjacent cells to current berg -Alex (periodicity fix)
            i = berg%ine; j = berg%jne
            do grdj_inner = j-2,j+2 ; do grdi_inner = i-2,i+2 !probably wider search than necessary...
              if (.not. bond_matched) then
                if ((grdj_inner .gt. grd%jsd-1) .and. (grdj_inner .lt. grd%jed+1) &
                    .and.  (grdi_inner .gt. grd%isd-1) .and. (grdi_inner .lt. grd%ied+1) ) then
                  other_berg=>bergs%list(grdi_inner,grdj_inner)%first
                  do while (associated(other_berg)) ! loop over all other bergs
                    if (other_berg%id .eq. current_bond%other_id) then
                      current_bond%other_berg=>other_berg
                      current_bond%other_berg_ine = grdi_inner
                      current_bond%other_berg_jne = grdj_inner
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

          if (check_match .and. .not.bond_matched .and. berg%halo_berg.ne.10) then
            if (berg%halo_berg .lt. 0.5) then
              missing_bond=.true.
              print * ,'non-halo berg unmatched: ', berg%id, berg%ine, berg%jne, mpp_pe(), &
                current_bond%other_id, current_bond%other_berg_ine, current_bond%other_berg_jne
              call error_mesg('KID, connect_all_bonds', 'A non-halo bond is missing!!!', FATAL)
            else  ! This is not a problem if the partner berg is not yet in the halo
              if (bergs%mts) then
                call error_mesg('KID, connect_all_bonds', 'A halo bond is missing!!!', WARNING)
                print *,'mpp_pe(),berg,target,ine,jne,halo_stat,lon,lat',&
                  mpp_pe(),berg%id,current_bond%other_id,i,j,berg%halo_berg,berg%lon,berg%lat
                print *,'mpp_pe(),latbounds:',mpp_pe(),grd%lat(grd%isd,grd%jsd),grd%lat(grd%isc-1,grd%jsc-1),&
                  grd%lat(grd%iec,grd%jec),grd%lat(grd%ied,grd%jed)
                print *,'mpp_pe(),lonbounds:',mpp_pe(),grd%lon(grd%isd,grd%jsd),grd%lon(grd%isc-1,grd%jsc-1),&
                  grd%lon(grd%iec,grd%jec),grd%lon(grd%ied,grd%jed)
                !bergs%writeandstop=.true.
              endif
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

!> If periodic domain, update lat & lon for bergs located beyond Lx - Alex
subroutine update_latlon(bergs)
  !Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg
  real :: dlon, dlat,dlon2,dlat2,xi,yj,Lx
  integer :: grdi,grdj
  logical :: lret

  ! For convenience
  grd=>bergs%grd
  Lx=grd%Lx

  if (Lx > 0.) then
    if ((grd%lon(grd%isc-1,grd%jsc-1) == grd%minlon_c) .or.  (grd%lon(grd%iec,grd%jec) == grd%maxlon_c)) then
      do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
        berg=>bergs%list(grdi,grdj)%first
        do while (associated(berg))
          if ((.not. bergs%mts) .or. (berg%halo_berg .le. 1)) then
            !for bergs%mts skip halo elements that only exist due to their
            !inclusion in a conglomerate (berg%halo_berg > 1). Their coords have
            !already been corrected (if needed), as the scheme below would fail.
            dlon = berg%lon - berg%lon_old
            dlat = berg%lat - berg%lat_old
            dlon2 = berg%lon - berg%lon_prev
            dlat2 = berg%lat - berg%lat_prev

            berg%lon=bilin(grd, grd%lon, berg%ine, berg%jne, berg%xi, berg%yj)
            berg%lat=bilin(grd, grd%lat, berg%ine, berg%jne, berg%xi, berg%yj)

            berg%lon_old = berg%lon-dlon
            berg%lat_old = berg%lat-dlat
            berg%lon_prev = berg%lon-dlon2
            berg%lat_prev = berg%lat-dlat2

            !need to update local positon in cell from updated global coords, due to roundoff
            lret=pos_within_cell(grd, berg%lon, berg%lat, berg%ine, berg%jne, berg%xi, berg%yj)
          endif
          berg=>berg%next
        enddo !loop over all bergs
      enddo; enddo
    endif
  endif

end subroutine update_latlon

!> Counts (and error checks) bonds
subroutine count_bonds(bergs, number_of_bonds, check_bond_quality)
type(icebergs), pointer :: bergs !< Container for all types and memory
integer, intent(out) :: number_of_bonds !< Number of bonds
logical, intent(inout), optional :: check_bond_quality !< If true, check bond quality
! Local variables
type(iceberg), pointer :: berg
type(iceberg), pointer :: other_berg
type(icebergs_gridded), pointer :: grd
type(bond) , pointer :: current_bond, other_berg_bond
integer :: number_of_bonds_all_pe
integer :: grdi, grdj
logical :: bond_is_good
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
        !      print *, 'Quality check', mpp_pe(), berg%id
        if (quality_check) then
          !num_unmatched_bonds=0
          !num_unassosiated_bond_pairs=0
          bond_is_good=.False.
          other_berg=>current_bond%other_berg
          if (associated(other_berg)) then
            other_berg_bond=>other_berg%first_bond
            do while (associated(other_berg_bond))  !loops over the icebergs in the other icebergs bond list
              if (associated(other_berg_bond%other_berg)) then
                if (other_berg_bond%other_berg%id .eq.berg%id) then
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
              !if (debug) write(stderrunit,*) 'Perfect quality Bond:', berg%id, current_bond%other_id
            else
              if (debug) write(stderrunit,*) 'Non-matching bond...:', berg%id, current_bond%other_id
              num_unmatched_bonds=num_unmatched_bonds+1
            endif
          else
            if (debug) write(stderrunit,*) 'Opposite berg is not assosiated:', berg%id, current_bond%other_id, mpp_pe()
            !num_unassosiated_bond_pairs=0
            num_unassosiated_bond_pairs=num_unassosiated_bond_pairs+1
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
        call error_mesg('KID, bonds', 'Bonds are not matching!', FATAL)
      endif
      if (num_unassosiated_bond_pairs_all_pe .ne. 0) then
        call error_mesg('KID, bonds', 'Bonds partners not located!', Warning)
        if (num_unassosiated_bond_pairs .ne. 0) then
          write(*,'(2a)') 'KID, Bonds parnters not located!!!! PE=', mpp_pe()
        endif
      endif
      if ((num_unmatched_bonds_all_pe .eq. 0)  .and. (num_unassosiated_bond_pairs_all_pe .eq. 0)) then
        if (mpp_pe().eq.mpp_root_pe()) then
                write(stderrunit,*)  "Total number of bonds is: ", number_of_bonds_all_PE, "All iceberg bonds are connected and working well"
        endif
        check_bond_quality=.true.
      else
        if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'KID: Warning, Broken Bonds! '
      endif
    endif

end subroutine count_bonds

!> Returns number of bergs across all lists
integer function count_bergs(bergs, with_halos)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
logical, optional :: with_halos !< If true, include halo lists
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

!> Returns number of bergs in a list
integer function count_bergs_in_list(first)
! Arguments
type(iceberg), pointer :: first !< List of bergs
! Local variables
type(iceberg), pointer :: this

  count_bergs_in_list=0
  this=>first
  do while(associated(this))
    count_bergs_in_list=count_bergs_in_list+1
    this=>this%next
  enddo

end function count_bergs_in_list

!> Add a record to the trajectory of each berg
subroutine record_posn(bergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
! Local variables
type(xyt) :: posn
type(iceberg), pointer :: this, other_berg
integer :: grdi, grdj
integer :: js,je,is,ie
type(bond) , pointer :: current_bond
real :: area_thres,area_thres2,area_thres3
! Local bond variables
type(bond_xyt) :: bond_posn
real :: dx_dlon,dy_dlat,lat_ref,dx,dy,n1,n2,id1,id2,berg_area
logical :: loc_save_nonfl_traj_by_class,save_fl_berg

if (bergs%debug_write) then
  js=bergs%grd%jsd+1; je=bergs%grd%jed; is=bergs%grd%isd+1; ie=bergs%grd%ied
else
  js=bergs%grd%jsc;   je=bergs%grd%jec; is=bergs%grd%isc;   ie=bergs%grd%iec
endif

  if (mts) then
    allocate(posn%axn_fast,posn%ayn_fast,posn%bxn_fast,posn%byn_fast)
  endif

  if (iceberg_bonds_on) then
    allocate(posn%n_bonds)
  endif

  if (monitor_energy) then
    allocate(posn%Ee,posn%Ed,posn%Eext,posn%Ee_contact,posn%Ed_contact,posn%Efrac)
    !allocate(posn%Ee_temp,posn%Ed_temp,posn%Eext_temp,posn%Ee_contact_temp,posn%Ed_contact_temp)
  endif

  if (dem) then
    allocate(posn%ang_vel,posn%ang_accel,posn%rot)
  elseif (fracture_criterion .ne. 'none') then
    allocate(posn%accum_bond_rotation)
  endif

  if (save_bond_traj) then
    if (use_damage) allocate(bond_posn%damage)
    if (monitor_energy) then
      allocate(bond_posn%Ee,bond_posn%Ed,bond_posn%axn_fast,bond_posn%ayn_fast,bond_posn%bxn_fast,bond_posn%byn_fast)
    endif
    if (dem) then
      allocate(bond_posn%tangd1,bond_posn%tangd2,bond_posn%nstress,bond_posn%sstress)
    elseif (fracture_criterion.ne.'none') then
      allocate(bond_posn%rotation,bond_posn%rel_rotation,bond_posn%n_frac_var)
      if (fracture_criterion.eq.'strain_rate') allocate(bond_posn%n_strain_rate)
      if (fracture_criterion.eq.'energy') allocate(bond_posn%spring_pe)
    endif
  endif

  area_thres=bergs%traj_area_thres*1.e6 ! convert from km^2 to m^2
  area_thres2=bergs%traj_area_thres_sntbc*1.e6
  area_thres3=bergs%traj_area_thres_fl*1.e6
  loc_save_nonfl_traj_by_class=.false.

  do grdj = js,je ; do grdi = is,ie
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      berg_area=this%mass/(bergs%rho_bergs*this%thickness)
      if (bergs%save_nonfl_traj_by_class) then
        loc_save_nonfl_traj_by_class=.false.
        if (this%fl_k>=0. .and. berg_area>area_thres2) then
          if (this%lat<0.) then
            if (this%start_mass>=bergs%save_traj_by_class_start_mass_thres_s) loc_save_nonfl_traj_by_class=.true.
          else
            if (this%start_mass>=bergs%save_traj_by_class_start_mass_thres_n) loc_save_nonfl_traj_by_class=.true.
          endif
        endif
      endif
      if (this%fl_k<0 .and. berg_area>area_thres3) then
        save_fl_berg=.true.
      else
        save_fl_berg=.false.
      endif
      current_bond=>this%first_bond
      if ( bergs%current_year>bergs%save_all_traj_year .or. loc_save_nonfl_traj_by_class .or. &
        berg_area >= area_thres .or. associated(current_bond) .or. save_fl_berg) then
        posn%lon=this%lon
        posn%lat=this%lat
        posn%year=bergs%current_year
        posn%day=bergs%current_yearday
        posn%id=this%id
        if (bergs%save_fl_traj) then
          posn%mass=this%mass
          posn%start_mass=this%start_mass
          posn%thickness=this%thickness
          posn%mass_of_bits=this%mass_of_bits
          posn%uvel=this%uvel
          posn%vvel=this%vvel
          if (bergs%fl_r>0) then
            posn%mass_scaling=this%mass_scaling
            posn%mass_of_fl_bits=this%mass_of_fl_bits
            posn%mass_of_fl_bergy_bits=this%mass_of_fl_bergy_bits
            posn%fl_k=this%fl_k
          endif
        endif
        if (.not. bergs%save_short_traj) then !Not totally sure that this is correct
          !posn%uvel=this%uvel
          !posn%vvel=this%vvel
          posn%uvel_prev=this%uvel_prev
          posn%vvel_prev=this%vvel_prev
          posn%heat_density=this%heat_density
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
          posn%halo_berg=this%halo_berg
          posn%static_berg=this%static_berg
          posn%od=this%od

          if (bergs%mts) then
            posn%axn_fast=this%axn_fast
            posn%ayn_fast=this%ayn_fast
            posn%bxn_fast=this%bxn_fast
            posn%byn_fast=this%byn_fast
          endif

          if (iceberg_bonds_on) then
            posn%n_bonds=this%n_bonds
          end if

          if (monitor_energy) then
            posn%Ee=this%Ee
            posn%Ed=this%Ed
            posn%Eext=this%Eext
            posn%Ee_contact=this%Ee_contact
            posn%Ed_contact=this%Ed_contact
            posn%Efrac=this%Efrac
            ! posn%Ee_contact_temp=this%Ee_contact_temp
            ! posn%Ed_contact_temp=this%Ed_contact_temp
            ! posn%Ee_temp=this%Ee_temp
            ! posn%Ed_temp=this%Ed_temp
            ! posn%Eext_temp=this%Eext_temp
          endif

          if (dem) then
            posn%ang_vel=this%ang_vel
            posn%ang_accel=this%ang_accel
            posn%rot=this%rot
          elseif (fracture_criterion .ne. 'none') then
            posn%accum_bond_rotation=this%accum_bond_rotation
          endif
        endif

        call push_posn(this%trajectory, posn)

        if (save_bond_traj) then
          do while (associated(current_bond)) ! loop over all bonds
            if  (associated(current_bond%other_berg)) then
              other_berg=>current_bond%other_berg

              ! From icebergs.F90 --> subroutine convert_from_grid_to_meters
              if (bergs%grd%grid_is_latlon) then
                lat_ref=0.5*(this%lat+other_berg%lat)
                dx_dlon=pi_180*Rearth*cos(lat_ref*pi_180);  dy_dlat=pi_180*Rearth
              else
                dx_dlon=1.0;                                dy_dlat=1.0
              endif

              dx=(this%lon-other_berg%lon)*dx_dlon; dy=(this%lat-other_berg%lat)*dy_dlat
              n1=dx/this%length; n2=dy/this%length !unit vector
              bond_posn%lon=this%lon-dx/2;          bond_posn%lat=this%lat-dy/2
              bond_posn%year=bergs%current_year;    bond_posn%day=bergs%current_yearday
              bond_posn%length=current_bond%length;
              bond_posn%n1=n1;                      bond_posn%n2=n2
              bond_posn%id1=this%id;                bond_posn%id2=other_berg%id

              if (use_damage) then
                bond_posn%damage=current_bond%damage
              endif

              if (monitor_energy) then
                bond_posn%Ee=current_bond%Ee
                bond_posn%Ed=current_bond%Ed
                bond_posn%axn_fast=current_bond%axn_fast
                bond_posn%ayn_fast=current_bond%ayn_fast
                bond_posn%bxn_fast=current_bond%bxn_fast
                bond_posn%byn_fast=current_bond%byn_fast
              endif

              if (dem) then
                bond_posn%tangd1=current_bond%tangd1
                bond_posn%tangd2=current_bond%tangd2
                bond_posn%nstress=current_bond%nstress
                bond_posn%sstress=current_bond%sstress
              elseif (fracture_criterion.ne.'none') then
                bond_posn%rotation=current_bond%rotation
                bond_posn%rel_rotation=current_bond%rel_rotation
                bond_posn%n_frac_var=current_bond%n_frac_var
                if (fracture_criterion.eq.'strain_rate') bond_posn%n_strain_rate=current_bond%n_strain_rate
                if (fracture_criterion.eq.'energy') bond_posn%spring_pe=current_bond%spring_pe
              endif

              call push_bond_posn(current_bond%bond_trajectory, bond_posn)
            endif
            current_bond=>current_bond%next_bond
          enddo
        endif
      endif

      this=>this%next
    enddo
  enddo ; enddo

end subroutine record_posn

!> Add trajectory values as a new record in a trajectory
subroutine push_posn(trajectory, posn_vals)
! Arguments
type(xyt), pointer :: trajectory !< Trajectory list
type(xyt) :: posn_vals !< Values to add
! Local variables
type(xyt), pointer :: new_posn

  allocate(new_posn)

  if (mts) then
    allocate(new_posn%axn_fast,new_posn%ayn_fast,new_posn%bxn_fast,new_posn%byn_fast)
  endif

  if (iceberg_bonds_on) then
    allocate(new_posn%n_bonds)
  endif

  if (monitor_energy) then
    allocate(new_posn%Ee,new_posn%Ed,new_posn%Eext,new_posn%Ee_contact,new_posn%Ed_contact,new_posn%Efrac)
    !allocate(new_posn%Ee_temp,new_posn%Ed_temp,new_posn%Eext_temp,new_posn%Ee_contact_temp,new_posn%Ed_contact_temp)
  endif

  if (dem) then
    allocate(new_posn%ang_vel,new_posn%ang_accel,new_posn%rot)
  elseif (fracture_criterion .ne. 'none') then
    allocate(new_posn%accum_bond_rotation)
  endif

  new_posn=posn_vals
  new_posn%next=>trajectory
  trajectory=>new_posn

end subroutine push_posn

!> Add bond trajectory values as a new record in a bond trajectory
subroutine push_bond_posn(bond_trajectory, bond_posn_vals)
! Arguments
type(bond_xyt), pointer :: bond_trajectory !< Trajectory list
type(bond_xyt) :: bond_posn_vals !< Values to add
! Local variables
type(bond_xyt), pointer :: new_bond_posn

  allocate(new_bond_posn)

  if (use_damage) allocate(new_bond_posn%damage)
  if (monitor_energy) then
    allocate(new_bond_posn%Ee,new_bond_posn%Ed,new_bond_posn%axn_fast,new_bond_posn%ayn_fast,&
      new_bond_posn%bxn_fast,new_bond_posn%byn_fast)
  endif
  if (dem) then
    allocate(new_bond_posn%tangd1,new_bond_posn%tangd2,new_bond_posn%nstress,new_bond_posn%sstress)
  elseif (fracture_criterion.ne.'none') then
    allocate(new_bond_posn%rotation,new_bond_posn%rel_rotation,new_bond_posn%n_frac_var)
    if (fracture_criterion.eq.'strain_rate') allocate(new_bond_posn%n_strain_rate)
    if (fracture_criterion.eq.'energy') allocate(new_bond_posn%spring_pe)
  endif

  new_bond_posn=bond_posn_vals
  new_bond_posn%next=>bond_trajectory
  bond_trajectory=>new_bond_posn

end subroutine push_bond_posn

!> Appends trajectory values to the end of the trajectory list (slow)
!! \todo append_posn() is very slow and should be removed a.s.a.p.
subroutine append_posn(trajectory, posn_vals)
! Arguments
type(xyt), pointer :: trajectory !< Trajectory list
type(xyt) :: posn_vals !< Values to add
! Local variables
type(xyt), pointer :: new_posn,next,last

  allocate(new_posn)

  if (mts) then
    allocate(new_posn%axn_fast,new_posn%ayn_fast,new_posn%bxn_fast,new_posn%byn_fast)
  endif

  if (iceberg_bonds_on) then
    allocate(new_posn%n_bonds)
  endif

  if (monitor_energy) then
    allocate(new_posn%Ee,new_posn%Ed,new_posn%Eext,new_posn%Ee_contact,new_posn%Ed_contact,new_posn%Efrac)
    !allocate(new_posn%Ee_temp,new_posn%Ed_temp,new_posn%Eext_temp,new_posn%Ee_contact_temp,new_posn%Ed_contact_temp)
  endif

  if (dem) then
    allocate(new_posn%ang_vel,new_posn%ang_accel,new_posn%rot)
  elseif (fracture_criterion .ne. 'none' .and. (.not. dem)) then
    allocate(new_posn%accum_bond_rotation)
  endif

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

!> Appends bond trajectory values to the end of the bond trajectory list (slow)
!! \todo append_bond_posn() is very slow and should be removed a.s.a.p.
subroutine append_bond_posn(bond_trajectory, bond_posn_vals)
! Arguments
type(bond_xyt), pointer :: bond_trajectory !< Trajectory list
type(bond_xyt) :: bond_posn_vals !< Values to add
! Local variables
type(bond_xyt), pointer :: new_bond_posn,next,last

  allocate(new_bond_posn)

  if (use_damage) allocate(new_bond_posn%damage)
  if (monitor_energy) then
    allocate(new_bond_posn%Ee,new_bond_posn%Ed,new_bond_posn%axn_fast,new_bond_posn%ayn_fast,&
      new_bond_posn%bxn_fast,new_bond_posn%byn_fast)
  endif
  if (dem) then
    allocate(new_bond_posn%tangd1,new_bond_posn%tangd2,new_bond_posn%nstress,new_bond_posn%sstress)
  elseif (fracture_criterion.ne.'none') then
    allocate(new_bond_posn%rotation,new_bond_posn%rel_rotation,new_bond_posn%n_frac_var)
    if (fracture_criterion.eq.'strain_rate') allocate(new_bond_posn%n_strain_rate)
    if (fracture_criterion.eq.'energy') allocate(new_bond_posn%spring_pe)
  endif

  new_bond_posn=bond_posn_vals
  new_bond_posn%next=>null()
  if(.NOT. associated(bond_trajectory)) then
     bond_trajectory=>new_bond_posn
  else
     ! Find end of the trajectory and point it to the  new leaf
     next=>bond_trajectory
     do while (associated(next))
        last=>next
        next=>next%next
     enddo
     last%next=>new_bond_posn
  endif
end subroutine append_bond_posn

!> Disconnect a trajectory from a berg and add it to a list of trajectory segments
subroutine move_trajectory(bergs, berg)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
type(iceberg), pointer :: berg !< Berg containing trajectory
! Local variables
type(xyt), pointer :: next, last
type(xyt) :: vals
type(bond) , pointer :: current_bond

  if (bergs%ignore_traj) return

  ! If the trajectory is empty, ignore it
  if (.not.associated(berg%trajectory)) return

  ! Find end of berg trajectory and point it to start of existing trajectories
  next=>berg%trajectory
  do while (associated(next))
    last=>next
    next=>next%next
  enddo
  last%next=>bergs%trajectories

  bergs%trajectories=>berg%trajectory

  if (save_bond_traj) then
    current_bond=>berg%first_bond
    do while (associated(current_bond)) ! loop over all bonds
      call move_bond_trajectory(bergs, current_bond)
      current_bond=>current_bond%next_bond
    enddo
  endif

  berg%trajectory=>null()

end subroutine move_trajectory

!> Disconnect a trajectory from a bond and add it to a list of trajectory segments
subroutine move_bond_trajectory(bergs, current_bond)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
type(bond) , pointer :: current_bond !< Bond containing trajectory
! Local variables
type(bond_xyt), pointer :: next, last
type(bond_xyt) :: vals

  if (bergs%ignore_traj) return

  ! If the trajectory is empty, ignore it
  if (.not.associated(current_bond%bond_trajectory)) return

  ! Find end of berg trajectory and point it to start of existing trajectories
  next=>current_bond%bond_trajectory
  do while (associated(next))
    last=>next
    next=>next%next
  enddo
  last%next=>bergs%bond_trajectories

  bergs%bond_trajectories=>current_bond%bond_trajectory
  current_bond%bond_trajectory=>null()

end subroutine move_bond_trajectory

!> Scan all bergs in a list and disconnect trajectories and more to the list of trajectory segments
!! \todo The argument delete_bergs should be removed.
subroutine move_all_trajectories(bergs, delete_bergs)
! Arguments
type(icebergs),    pointer    :: bergs !< Container for all types and memory
logical, optional, intent(in) :: delete_bergs !< If true, delete bergs after disconnecting its trajectory
! Local variables
type(iceberg), pointer :: this, next
logical :: delete_bergs_after_moving_traj
integer :: grdi, grdj, is, ie, js, je

  if (bergs%ignore_traj) return

  delete_bergs_after_moving_traj = .false.
  if (present(delete_bergs)) delete_bergs_after_moving_traj = delete_bergs

  if (bergs%debug_write) then
    js=bergs%grd%jsd+1; je=bergs%grd%jed; is=bergs%grd%isd+1; ie=bergs%grd%ied
  else
    js=bergs%grd%jsc;   je=bergs%grd%jec; is=bergs%grd%isc;   ie=bergs%grd%iec
  endif

  !do grdj = bergs%grd%jsc,bergs%grd%jec ; do grdi = bergs%grd%isc,bergs%grd%iec
  do grdj = js,je ; do grdi = is,ie
    this=>bergs%list(grdi,grdj)%first
    do while (associated(this))
      next=>this%next
      call move_trajectory(bergs, this)
   !  if (delete_bergs_after_moving_traj) call destroy_iceberg(this)
      this=>next
    enddo
  enddo ; enddo

end subroutine move_all_trajectories

!> Search the grid for a cell containing position x,y
logical function find_cell_by_search(grd, x, y, i, j)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: x !< Longitude of position
real, intent(in) :: y !< Latitude of position
integer, intent(inout) :: i !< i-index of cell containing x,y
integer, intent(inout) :: j !< j-index of cell containing x,y
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
    call error_mesg('KID, find_cell_by_search:', 'This should never EVER happen! (1)', FATAL)
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
      call error_mesg('KID, find_cell_by_search:', 'This should never EVER happen!', FATAL)
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
  !       write(0,'(i3,a,2i5,a,2i3,a,2f8.3)') mpp_pe(),'KID, find_cell_by_search: false negative io,jo=',io,jo,' di,dj=',di,dj,' targ=',x,y
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
    write(0,'(i3,a,2i5,a,2i3)') mpp_pe(),'KID, find_cell_by_search: false negative 2 i,j=',i-is,j-js,' di,dj=',di,dj
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'KID, find_cell_by_search: false negative 2 io,jo=',io,jo
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'KID, find_cell_by_search: false negative 2 i,j=',i,j,' targ=',x,y
    return
  endif
  find_cell_by_search=.false.

end function find_cell_by_search

!< A cost function with minimum when x2==x1 and y2==y1
real function dcost(x1, y1, x2, y2, Lx)
  ! Arguments
  real, intent(in) :: x1 !< x-coordinate or point 1
  real, intent(in) :: x2 !< x-coordinate or point 2
  real, intent(in) :: y1 !< y-coordinate or point 1
  real, intent(in) :: y2 !< y-coordinate or point 2
  real, intent(in) :: Lx !< Periodicity of x-coordinate
  ! Local variables
  real :: x1m

    x1m=apply_modulo_around_point(x1,x2,Lx)
  ! dcost=(x2-x1)**2+(y2-y1)**2
    dcost=(x2-x1m)**2+(y2-y1)**2
end function dcost

!> A better way to find the indices (oi,oj) or the cell containing position (x,y).
!! Returns True if a cell is found and updates oi,oj.
logical function find_better_min(grd, x, y, w, oi, oj)
  ! Arguments
  type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
  real, intent(in) :: x !< x-coordinate being searched for
  real, intent(in) :: y !< y-coordinate being searched for
  integer, intent(in) :: w !< Number of surrounding cells to the N,S,E, and W to search
  integer, intent(inout) :: oi !< i-index of cell containing x,y
  integer, intent(inout) :: oj !< j-index of cell containing x,y
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

! Returns true and updates oi and oj, if a cell containing x,y is found within w cells of the input oi, oj
logical function find_cell_loc(grd, x, y, is, ie, js, je, w, oi, oj)
  ! Arguments
  type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
  real, intent(in) :: x !< Longitude
  real, intent(in) :: y !< Latitude
  integer, intent(in) :: is !< Start i-index of computational grid
  integer, intent(in) :: ie !< End i-index of computational grid
  integer, intent(in) :: js !< Start j-index of computational grid
  integer, intent(in) :: je !< End j-index of computational grid
  integer, intent(in) :: w !< Width to search around
  integer, intent(inout) :: oi !< Starting and updated i-index
  integer, intent(inout) :: oj !< Starting and updated j-index
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

!> Returns the i,j of cell containing an iceberg with the given identifier
subroutine find_individual_iceberg(bergs, id, ine, jne, berg_found, search_data_domain)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
integer(kind=8), intent(in) :: id !< Berg identifier
integer, intent(out) :: ine !< i-index of cell containing berg
integer, intent(out) :: jne !< j-index of cell containing berg
logical, intent(in) :: search_data_domain !< If true, search halos too
real, intent(out) :: berg_found !< Returns 1.0 if berg is found, 0. otherwise
! Local variables
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
integer :: grdi, grdj
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
        if (id .eq. this%id) then
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

!> Scans each computational grid cell until is_point_in_cell() is true,
!! checking the current cell first.
logical function check_and_find_cell(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
real, intent(in) :: x !< Longitude of position
real, intent(in) :: y !< Latitude of position
integer, intent(out) :: oi !< i-index of cell containing position or -999
integer, intent(out) :: oj !< j-index of cell containing position or -999
! Local variables
integer :: i,j

check_and_find_cell=.false.

if (.not. (oi-1.lt.grd%isd.or.oi.gt.grd%ied.or.oj-1.lt.grd%jsd.or.oj.gt.grd%jed)) then
  if (is_point_in_cell(grd, x, y, oi, oj)) then
    check_and_find_cell=.true.
        return
      endif
endif

!find cell for structured grid without looping:
oi=floor((x-grd%lon(grd%isd,grd%jsd))/(grd%lon(grd%isd+1,grd%jsd+1)-grd%lon(grd%isd,grd%jsd)))+grd%isd+1
oj=floor((y-grd%lat(grd%isd,grd%jsd))/(grd%lat(grd%isd+1,grd%jsd+1)-grd%lat(grd%isd,grd%jsd)))+grd%jsd+1
if (.not. (oi-1.lt.grd%isd.or.oi.gt.grd%ied.or.oj-1.lt.grd%jsd.or.oj.gt.grd%jed)) then
  if (is_point_in_cell(grd, x, y, oi, oj)) then
    check_and_find_cell=.true.; return
  endif
endif

  oi=-999; oj=-999
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; check_and_find_cell=.true.
        return
      endif
    enddo; enddo
end function check_and_find_cell

!> Scans each computational grid cell until is_point_in_cell() is true
logical function find_cell(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
real, intent(in) :: x !< Longitude of position
real, intent(in) :: y !< Latitude of position
integer, intent(out) :: oi !< i-index of cell containing position or -999
integer, intent(out) :: oj !< j-index of cell containing position or -999
! Local variables
integer :: i,j

find_cell=.false.


!find cell for structured grid without looping:
oi=floor((x-grd%lon(grd%isd,grd%jsd))/(grd%lon(grd%isd+1,grd%jsd+1)-grd%lon(grd%isd,grd%jsd)))+grd%isd+1
oj=floor((y-grd%lat(grd%isd,grd%jsd))/(grd%lat(grd%isd+1,grd%jsd+1)-grd%lat(grd%isd,grd%jsd)))+grd%jsd+1
if (oi>grd%isc-1 .and. oi<grd%iec+1 .and. oj>grd%jsc-1 .and. oj<grd%jec+1) then
  if (is_point_in_cell(grd, x, y, oi, oj)) then
    find_cell=.true.; return
  endif
endif

  oi=-999; oj=-999
  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell=.true.
        return
      endif
  enddo; enddo

end function find_cell

!> Scans each all grid cells until is_point_in_cell() is true (includes halos)
logical function find_cell_wide(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
real, intent(in) :: x !< Longitude of position
real, intent(in) :: y !< Latitude of position
integer, intent(out) :: oi !< i-index of cell containing position or -999
integer, intent(out) :: oj !< j-index of cell containing position or -999
! Local variables
integer :: i,j

find_cell_wide=.false.
!find cell for structured grid without looping:
oi=floor((x-grd%lon(grd%isd,grd%jsd))/(grd%lon(grd%isd+1,grd%jsd+1)-grd%lon(grd%isd,grd%jsd)))+grd%isd+1
oj=floor((y-grd%lat(grd%isd,grd%jsd))/(grd%lat(grd%isd+1,grd%jsd+1)-grd%lat(grd%isd,grd%jsd)))+grd%jsd+1
if (.not. (oi-1.lt.grd%isd.or.oi.gt.grd%ied.or.oj-1.lt.grd%jsd.or.oj.gt.grd%jed)) then
  if (is_point_in_cell(grd, x, y, oi, oj)) then
    find_cell_wide=.true.; return
  endif
endif

oi=-999; oj=-999

  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell_wide=.true.
        return
      endif
  enddo; enddo

end function find_cell_wide

!> Returns True if x,y is in cell i,j
logical function is_point_in_cell(grd, x, y, i, j, explain)
! Arguments
type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
real, intent(in) :: x !< Longitude of position
real, intent(in) :: y !< Latitude of position
integer, intent(in) :: i !< i-index of cell
integer, intent(in) :: j !< j-index of cell
logical, intent(in), optional :: explain !< If true, print debugging
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
                     'KID, is_point_in_cell: pe=(',mpp_pe(),') i,s,e=', &
                     i,grd%isd,grd%ied,' j,s,e=', j,grd%jsd,grd%jed
    call error_mesg('KID, is_point_in_cell', 'test is off the PE!', FATAL)
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

!> Returns true if point x,y is inside polygon with four corners
logical function sum_sign_dot_prod4(x0, y0, x1, y1, x2, y2, x3, y3, x, y, Lx, explain)
! Arguments
real, intent(in) :: x0 !< Longitude of first corner
real, intent(in) :: y0 !< Latitude of first corner
real, intent(in) :: x1 !< Longitude of second corner
real, intent(in) :: y1 !< Latitude of second corner
real, intent(in) :: x2 !< Longitude of third corner
real, intent(in) :: y2 !< Latitude of third corner
real, intent(in) :: x3 !< Longitude of fourth corner
real, intent(in) :: y3 !< Latitude of fourth corner
real, intent(in) :: x !< Longitude of point
real, intent(in) :: y !< Latitude of point
real, intent(in) :: Lx !< Length of domain in zonal direction
logical, intent(in), optional :: explain !< If true, print debugging
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

  !We use an asymmetry between South and East line boundaries and North and East
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

!> Returns true if point x,y is inside polygon with five corners
logical function sum_sign_dot_prod5(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y, Lx, explain)
! Arguments
real, intent(in) :: x0 !< Longitude of first corner
real, intent(in) :: y0 !< Latitude of first corner
real, intent(in) :: x1 !< Longitude of second corner
real, intent(in) :: y1 !< Latitude of second corner
real, intent(in) :: x2 !< Longitude of third corner
real, intent(in) :: y2 !< Latitude of third corner
real, intent(in) :: x3 !< Longitude of fourth corner
real, intent(in) :: y3 !< Latitude of fourth corner
real, intent(in) :: x4 !< Longitude of fifth corner
real, intent(in) :: y4 !< Latitude of fifth corner
real, intent(in) :: x !< Longitude of point
real, intent(in) :: y !< Latitude of point
real, intent(in) :: Lx !< Length of domain in zonal direction
logical, intent(in), optional :: explain !< If true, print debugging
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

!> Calculates non-dimensional position with cell i,j and returns true if point is in cell
logical function pos_within_cell(grd, x, y, i, j, xi, yj, explain)
! Arguments
type(icebergs_gridded), intent(in) :: grd !< Container for gridded fields
real, intent(in) :: x !< Longitude of position
real, intent(in) :: y !< Latitude of position
integer, intent(in) :: i !< i-index of cell
integer, intent(in) :: j !< j-index of cell
real, intent(out) :: xi !< Non-dimensional x-position within cell
real, intent(out) :: yj !< Non-dimensional y-position within cell
logical, intent(in), optional :: explain !< If true, print debugging
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
        if (debug) call error_mesg('KID, pos_within_cell', 'in cell but scaling internal coordinates!', WARNING)
      endif
    else
      ! The point is not inside the spherical quad
      if (abs(xi-0.5)<0.5.and.abs(yj-0.5)<0.5) then
        ! The projection of the spherical quad onto the tangent plane should be larger than
        ! quad in the tangent plane so we should never be able to get here.
        call error_mesg('KID, pos_within_cell', 'not in cell but coordinates <0.5!', FATAL)
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
        call error_mesg('KID, pos_within_cell', 'pos_within_cell is False BUT is_point_in_cell disagrees!', FATAL)
      endif
    endif
  else
    ! Based on coordinate, the point is within cell
    if (pos_within_cell) then
      ! Based on is_point_in_cell() the point is out of cell so we have an inconsistency
      if (debug) then
        write(stderrunit,'(a,1p6e12.4)') 'values of xi, yj ',xi, yj
        pos_within_cell=is_point_in_cell(grd, x, y, i, j,explain=.True.)
        call error_mesg('KID, pos_within_cell', 'pos_within_cell is True BUT is_point_in_cell disagrees!', FATAL)
      endif
    endif
  endif

end function pos_within_cell

!> Calculates non-dimension position of x,y within a polygon with four corners
subroutine calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, Lx, explain)
! Arguments
real, intent(in) :: x1 !< Longitude of first corner
real, intent(in) :: y1 !< Latitude of first corner
real, intent(in) :: x2 !< Longitude of second corner
real, intent(in) :: y2 !< Latitude of second corner
real, intent(in) :: x3 !< Longitude of third corner
real, intent(in) :: y3 !< Latitude of third corner
real, intent(in) :: x4 !< Longitude of fourth corner
real, intent(in) :: y4 !< Latitude of fourth corner
real, intent(in) :: x !< Longitude of point
real, intent(in) :: y !< Latitude of point
real, intent(out) :: xi !< Non-dimensional x-position within cell
real, intent(out) :: yj !< Non-dimensional y-position within cell
real, intent(in) :: Lx !< Length of domain in zonal direction
logical, intent(in), optional :: explain !< If true, print debugging
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
    call error_mesg('KID, calc_xiyj', 'We have complex roots. The grid must be very distorted!', FATAL)
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
    call error_mesg('KID, calc_xiyj', 'Can not invert either linear equaton for xi! This should not happen!', FATAL)
  endif
endif
if (expl) write(stderrunit,'(a,2e12.4)') 'calc_xiyj: xi,yj=',xi,yj

end subroutine calc_xiyj

!> Returns true if non-dimensional position xi,yj is in unit interval
!!
!! Includes South and East boundaries, and excludes North and West.
!! \todo Double check definition of is_point_within_xi_yj_bounds()
logical function is_point_within_xi_yj_bounds(xi,yj)
! Arguments
real, intent(in) :: xi !< Non-dimensional x-position
real, intent(in) :: yj !< Non-dimensional y-position
! Local variables
!Includes South and East boundaries, and excludes North and West  (double check this is the way that is needed)
  is_point_within_xi_yj_bounds=.False.
  if ((xi .ge. 0 )  .and.  (xi .lt. 1)) then
    if ((yj .ge. 0 )  .and.  (yj .lt. 1)) then
      is_point_within_xi_yj_bounds=.True.
    endif
  endif
end function is_point_within_xi_yj_bounds

!> Modulo value of x in an interval [y-(Lx/2)  y+(Lx/2)]
!!
!! Gives the modulo value of x in an interval [y-(Lx/2)  y+(Lx/2)]  , modulo Lx
!! If Lx<=0, then it returns x without applying modulo arithmetic.
real function apply_modulo_around_point(x, y, Lx)
! Arguments
real, intent(in) :: x !< Value to apply modulo arithmetic to
real, intent(in) :: y !< Center of modulo range
real, intent(in) :: Lx !< Modulo width
!Local_variables
real ::Lx_2

  if (Lx>0.) then
    Lx_2=Lx/2.
    apply_modulo_around_point=modulo(x-(y-Lx_2),Lx)+(y-Lx_2)
  else
    apply_modulo_around_point=x
  endif

end function apply_modulo_around_point

!> Checks that a berg's position metrics are consistent
subroutine check_position(grd, berg, label, il, jl)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
type(iceberg), pointer :: berg !< Berg to check
character(len=*) :: label !< Label to add to messages
integer, optional, intent(in) :: il !< i-index of cell berg should be in
integer, optional, intent(in) :: jl !< j-index of cell berg should be in
! Local variables
real :: xi, yj
logical :: lret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  lret=pos_within_cell(grd, berg%lon, berg%lat, berg%ine, berg%jne, xi, yj)
  if (xi.ne.berg%xi.or.yj.ne.berg%yj) then
    write(stderrunit,'("KID: check_position (",i4,") b%x,x,-=",3(es12.4,x),a)') mpp_pe(),berg%xi,xi,berg%xi-xi,label
    write(stderrunit,'("KID: check_position (",i4,") b%y,y,-=",3(es12.4,x),a)') mpp_pe(),berg%yj,yj,berg%yj-yj,label
    call print_berg(stderrunit, berg, 'check_position', il, jl)
    call error_mesg('KID, check_position, '//trim(label),'berg has inconsistent xi,yj!',FATAL)
  endif
  if (grd%msk(berg%ine, berg%jne)==0.) then
    call print_berg(stderrunit, berg, 'check_position, '//trim(label), il, jl)
    call error_mesg('KID, check_position, '//trim(label),'berg is in a land cell!',FATAL)
  endif

end subroutine check_position

!> Add up the mass of bergs and/or bergy bits
real function sum_mass(bergs, justbits, justflbits, justbergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
logical, intent(in), optional :: justbits !< If present, add up mass of just bergy bits
logical, intent(in), optional :: justflbits !< If present, add up mass of just footloose bits
logical, intent(in), optional :: justbergs !< If present, add up mass of just bergs
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
        sum_mass=sum_mass+(this%mass_of_bits+this%mass_of_fl_bergy_bits)*this%mass_scaling
      elseif (present(justflbits)) then
        sum_mass=sum_mass+this%mass_of_fl_bits*this%mass_scaling
      else
        sum_mass=sum_mass+(this%mass+this%mass_of_bits+this%mass_of_fl_bits+this%mass_of_fl_bergy_bits)*this%mass_scaling
      endif
      this=>this%next
    enddo
  enddo ; enddo

end function sum_mass

!> Add up the heat content of bergs and/or bergy bits
real function sum_heat(bergs,justbits,justflbits, justbergs)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
logical, intent(in), optional :: justbits !< If present, add up heat content of just bergy bits
logical, intent(in), optional :: justflbits !< If present, add up mass of just footloose bits
logical, intent(in), optional :: justbergs !< If present, add up heat content of just bergs
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
        dm=(this%mass_of_bits+this%mass_of_fl_bergy_bits)*this%mass_scaling
      elseif (present(justflbits)) then
        dm=this%mass_of_fl_bits*this%mass_scaling
      else
        dm=(this%mass+this%mass_of_bits+this%mass_of_fl_bits+this%mass_of_fl_bergy_bits)*this%mass_scaling
      endif
      sum_heat=sum_heat+dm*this%heat_density
      this=>this%next
    enddo
  enddo ; enddo

end function sum_heat

!> Set elements of an array to 0. if the absolute value is larger than a given value
subroutine sanitize_field(arr, val)
! Arguments
real, dimension(:,:),intent(inout) :: arr !< Array to sanitize
real, intent(in) :: val !< Threshold value to use for sanitizing
! Local variables
integer :: i, j

  do j=lbound(arr,2), ubound(arr,2)
    do i=lbound(arr,1), ubound(arr,1)
      if (abs(arr(i,j)).ge.val) arr(i,j)=0.
    enddo
  enddo

end subroutine sanitize_field

!> Calculates checksums for all gridded fields
subroutine checksum_gridded(grd, label)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
character(len=*) :: label !< Label to use in messages
! Local variables

  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'KID: checksumming gridded data @ ',trim(label)

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
  call grd_chksum2(grd, grd%fl_bits_src, 'fl_bits_src')
  call grd_chksum2(grd, grd%fl_bits_melt, 'fl_bits_melt')
  call grd_chksum2(grd, grd%fl_bits_mass, 'fl_bits_mass')
  call grd_chksum2(grd, grd%fl_bergy_bits_mass, 'fl_bergy_bits_mass')
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
  call grd_chksum3(grd, grd%melt_by_class, 'melt_by_class')
  call grd_chksum2(grd, grd%melt_buoy_fl, 'melt_b_fl')
  call grd_chksum2(grd, grd%melt_eros_fl, 'melt_e_fl')
  call grd_chksum2(grd, grd%melt_conv_fl, 'melt_v_fl')
  call grd_chksum2(grd, grd%fl_parent_melt, 'fl_parent_melt')
  call grd_chksum2(grd, grd%fl_child_melt, 'fl_child_melt')

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

!> Calculates checksum for a 3d field
subroutine grd_chksum3(grd, fld, txt)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
real, dimension(:,:,:), intent(in) :: fld !< Field to checksum
character(len=*), intent(in) :: txt !< Label to use in message
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
    write(*,'("KID, grd_chksum3: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  j=mpp_chksum( tmp(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("KID, grd_chksum3* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#endif

end subroutine grd_chksum3

!> Calculates checksum for a 2d field
subroutine grd_chksum2(grd, fld, txt)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: fld !< Field to checksum
character(len=*), intent(in) :: txt !< Label to use in message
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
    write(*,'("KID, grd_chksum2: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i8)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(grd%isd:grd%ied,grd%jsd:grd%jed) )
  j=mpp_chksum( grd%tmp(grd%isd:grd%ied,grd%jsd:grd%jed) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("KID, grd_chksum2* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#endif

end subroutine grd_chksum2

!> Calculates checksums for all bergs
subroutine bergs_chksum(bergs, txt, ignore_halo_violation)
! Arguments
type(icebergs), pointer :: bergs !< Container for all types and memory
character(len=*), intent(in) :: txt !< Label to use in messages
logical, optional :: ignore_halo_violation !< If true, skip error checking
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
    write(*,'("KID, bergs_chksum: ",2(a,i8))') &
      '# bergs =', nbergs, ' sum(icnt) =',sum(icnt(:,:))
    call error_mesg('KID, bergs_chksum:', 'mismatch in berg count!', FATAL)
  endif

  check_halo=.true.
  if (present(ignore_halo_violation)) then
    if (ignore_halo_violation) check_halo=.false.
  endif
  if (check_halo.and.nbergs.ne.sum(icnt(grd%isc:grd%iec, grd%jsc:grd%jec))) then
    write(*,'("KID, bergs_chksum: ",2(a,i8))') &
      '# bergs =', nbergs, ' sum(icnt(comp_dom)) =',sum(icnt(:,:))
    call error_mesg('KID, bergs_chksum:', 'mismatch in berg count on computational domain!', FATAL)
  endif

  call mpp_sum(nbergs)
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("KID, bergs_chksum: ",a18,6(x,a,"=",i22))') &
      txt, 'chksum', ichk1, 'chksum2', ichk2, 'chksum3', ichk3, 'chksum4', ichk4, 'chksum5', ichk5, '#', nbergs

  grd%tmp(:,:)=real(icnt(:,:))
  call grd_chksum2(grd,grd%tmp,'# of bergs/cell')

  deallocate( fld )
  deallocate( fld2 )
  deallocate( icnt )

  if (debug) call count_out_of_order(bergs,txt)

end subroutine bergs_chksum

!> Checksum a list of bergs
integer function list_chksum(first)
! Arguments
type(iceberg), pointer :: first !< List of bergs
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

!> Checksum a berg
integer function berg_chksum(berg)
! Arguments
type(iceberg), pointer :: berg !< An iceberg
! Local variables
real :: rtmp(36)
integer :: itmp(36+7), i8=0, ichk1, ichk2, ichk3
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
  itmp(1:36)=transfer(rtmp,i8)
  itmp(37)=berg%halo_berg
  itmp(38)=berg%static_berg
  itmp(39)=berg%start_year
  itmp(40)=berg%ine
  itmp(41)=berg%jne
  call split_id( berg%id, itmp(42), itmp(43) )

  ichk1=0; ichk2=0; ichk3=0
  do i=1,36+7
   ichk1=ichk1+itmp(i)
   ichk2=ichk2+itmp(i)*i
   ichk3=ichk3+itmp(i)*i*i
  enddo
  berg_chksum=ichk1+ichk2+ichk3

end function berg_chksum

!> Bi-linear interpolate a field at corners in cell i,j to non-dimensional position xi,yj
real function bilin(grd, fld, i, j, xi, yj)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed) !< Field to interpolate
real, intent(in) :: xi !< Non-dimensional x-position within cell
real, intent(in) :: yj !< Non-dimensional y-position within cell
integer, intent(in) :: i !< i-index of cell
integer, intent(in) :: j !< j-index of cell
! Local variables

  if (old_bug_bilin) then
    bilin=(fld(i,j  )*(1.-xi)+fld(i-1,j  )*xi)*(1.-yj) &
         +(fld(i,j-1)*(1.-xi)+fld(i-1,j-1)*xi)*yj
  else
    bilin=(fld(i,j  )*xi+fld(i-1,j  )*(1.-xi))*yj &
         +(fld(i,j-1)*xi+fld(i-1,j-1)*(1.-xi))*(1.-yj)
  endif
end function bilin

!> Quadratic interpolation of an A-grid field to a location within a cell
real function quad_interp_from_agrid(grd, fld, x, y, i, j, xi, yj)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed) !< Field to interpolate
real, intent(in) :: x !<Longitude of position
real, intent(in) :: y !<Latitude of position
real, intent(in) :: xi !< Non-dimensional x-position within cell
real, intent(in) :: yj !< Non-dimensional y-position within cell
integer, intent(in) :: i !< i-index of cell
integer, intent(in) :: j !< j-index of cell
! Local variables
integer :: is,ie,js,je
real :: x1,x2,x3,x4,y1,y2,y3,y4
real :: dx,dy,Delta_x,xx,yy
real :: xloc,yloc !local coords
real :: xb(3,3), yb(3,3) !basis functions

!Uses bi-quadratic lagrange quadrilateral basis functions (like in FEM), where the 9
!nodes of an "element" are defined by a 3x3 array of neigboring cell-centers

!bounds of "node" array (is:ie,js:je).
if (mod(i,2)==1) then
  if (xi>=0.5) then
    is=i; ie=i+2
  else
    is=i-2; ie=i
  endif
else
  is=i-1; ie=i+1
endif
if (mod(j,2)==1) then
  if (yj>=0.5) then
    js=j; je=j+2
  else
    js=j-2; je=j
  endif
else
  js=j-1; je=j+1
endif

!four corners
x1=grd%lonc(is,js); y1=grd%latc(is,js)
x2=grd%lonc(ie,js); y2=grd%latc(ie,js)
x3=grd%lonc(ie,je); y3=grd%latc(ie,je)
x4=grd%lonc(is,je); y4=grd%latc(is,je)

!get non-dimensional, local (natural) coords xloc and yloc:
!largely copied from to pos_within_cell
if ((.not. grd%grid_is_latlon) .and. (grd%grid_is_regular))  then
  dx=abs(x3-x4); dy=abs(y3-y2)
  x1=x3-(dx/2); y1=y3-(dy/2)
  Delta_x= apply_modulo_around_point(x,x1,grd%Lx)-x1
  xloc=((Delta_x)/dx)+0.5; yloc=((y-y1)/dy)+0.5
elseif ((max(y1,y2,y3,y4)<89.999) .or.(.not. grd%grid_is_latlon)) then
  ! This returns non-dimensional position xi,yj for quad cells (not at a pole)
  call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xloc, yloc, grd%Lx)
else
  ! One of the cell corners is at the north pole so we switch to a tangent plane with
  ! co-latitude as a radial coordinate.
  xx=(90.-y)*cos(x*pi_180); yy=(90.-y)*sin(x*pi_180)
  x1=(90.-y1)*cos(grd%lon(is,js)*pi_180); y1=(90.-y1)*sin(grd%lon(is,js)*pi_180)
  x2=(90.-y2)*cos(grd%lon(ie,je)*pi_180); y2=(90.-y2)*sin(grd%lon(ie,je)*pi_180)
  x3=(90.-y3)*cos(grd%lon(ie,je)*pi_180); y3=(90.-y3)*sin(grd%lon(ie,je)*pi_180)
  x4=(90.-y4)*cos(grd%lon(is,je)*pi_180); y4=(90.-y4)*sin(grd%lon(is,je)*pi_180)
  call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, xx, yy, xloc, yloc, grd%Lx)
endif

!convert local coords to lie within the range [-1,1] rather than [0,1]
xloc=xloc*2-1; yloc=yloc*2-1

!basis functions:
xb(1,:)=0.5*xloc*(xloc-1); yb(:,1)=0.5*yloc*(yloc-1)
xb(2,:)=(1+xloc)*(1-xloc); yb(:,2)=(1+yloc)*(1-yloc)
xb(3,:)=0.5*xloc*(xloc+1); yb(:,3)=0.5*yloc*(yloc+1)

!interpolate:
quad_interp_from_agrid=sum(xb*yb*fld(is:ie,js:je))

end function quad_interp_from_agrid


!> Prints a field
subroutine print_fld(grd, fld, label)
! Arguments
type(icebergs_gridded), pointer :: grd !< Container for gridded fields
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed) !< Field to print
character(len=*) :: label !< Label to use in title
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

!> Combine a counter and i,j hash into an id
integer(kind=8) function id_from_2_ints(counter, ijhash)
  integer, intent(in) :: counter !< The counter value assigned at calving
  integer, intent(in) :: ijhash  !< A hash of i,j calving location

  id_from_2_ints = int(counter,8) * (int(2,8)**32) + int(ijhash,8)

end function id_from_2_ints

!> Split an iceberg ID into two parts
subroutine split_id(id, counter, ijhash)
  integer(kind=8), intent(in)  :: id      !< A unique id assigned when a berg is created
  integer,         intent(out) :: counter !< The counter value assigned at calving
  integer,         intent(out) :: ijhash  !< A hash of i,j calving location
  ! Local variables
  integer(kind=8) :: i8

  counter = ishft(id,-32)
  !counter = i8
  ijhash = int(id,4)

end subroutine split_id

!> Invoke some unit tests
logical function unit_tests(bergs)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  integer :: stderrunit,i,j,c1,c2
  integer(kind=8) :: id

  ! This function returns True is a unit test fails
  unit_tests=.false.
  ! For convenience
  grd=>bergs%grd
  stderrunit=stderr()

  i=grd%isc; j=grd%jsc
  call localTest( unit_tests, bilin(grd, grd%lon, i, j, 0., 1.), grd%lon(i-1,j) )
  call localTest( unit_tests, bilin(grd, grd%lon, i, j, 1., 1.), grd%lon(i,j) )
  call localTest( unit_tests, bilin(grd, grd%lat, i, j, 1., 0.), grd%lat(i,j-1) )
  call localTest( unit_tests, bilin(grd, grd%lat, i, j, 1., 1.), grd%lat(i,j) )

  ! Test 64-bit ID conversion
  i = 1440*1080 ; c1 = 2**30 + 2**4 + 1
  id = id_from_2_ints(c1, i)
  call split_id(id,c2,j)
  if (j /= i .or. c2 /= c1) then
    write(0,*) 'i,c in:',i,c1,' id=',id,' i,c out:',j,c2
    unit_tests=.true.
  endif

end function unit_tests

!> Checks answer to right answer and prints results if different
subroutine localTest(unit_test, answer, rightAnswer)
  logical, intent(inout) :: unit_test !< Set to true answer is wrong
  real, intent(in) :: answer !< Calculated answer
  real, intent(in) :: rightAnswer !< Correct answer
  ! Local variables
  integer :: stderrunit
  stderrunit=stderr()
  if (answer==rightAnswer) return
  unit_test=.true.
  write(stderrunit,*) 'a=',answer,'b=',rightAnswer
end subroutine localTest

!> Check for duplicates of icebergs on and across processors and issue an error
!! if any are detected
subroutine check_for_duplicates_in_parallel(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: this
  integer :: stderrunit, i, j, k, nbergs, nbergs_total
  integer(kind=8), dimension(:), allocatable :: ids ! List of ids of all bergs on this processor

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
        ids(k) = this%id
        this=>this%next
      enddo
    enddo ; enddo
  endif
  if (k /= nbergs) then
    write(stderrunit,*) 'counted bergs=',k,'count_bergs()=',nbergs
    call error_mesg('KID, check_for_duplicates:', 'Mismatch between concatenation of lists and count_bergs()!', FATAL)
  endif
  k = check_for_duplicate_ids_in_list(nbergs, ids, verbose=.true.)
  if (k /= 0) call error_mesg('KID, check_for_duplicates:', 'Duplicate berg detected across PEs!', FATAL)
  if (nbergs>0) deallocate(ids)
end subroutine check_for_duplicates_in_parallel

!> Returns error count of duplicates of integer values in a distributed list
integer function check_for_duplicate_ids_in_list(nbergs, ids, verbose)
  ! Arguments
  integer,                       intent(in)    :: nbergs !< Length of ids
  integer(kind=8), dimension(:), intent(inout) :: ids !< List of ids
  logical,                       intent(in)    :: verbose !< True if messages should be written
  ! Local variables
  integer :: stderrunit, i, j, k, nbergs_total, ii
  integer(kind=8) :: lowest_id, nonexistent_id, id, lid
  logical :: have_berg

  stderrunit=stderr()
  nbergs_total = nbergs
  call mpp_sum(nbergs_total) ! Total number of bergs

  ! Establish lowest id or 0 across all PEs
  lowest_id = 0
  if (nbergs>0) lowest_id = minval(ids)
  call mpp_min(lowest_id)
  id = lowest_id
  nonexistent_id = lowest_id - 1
  if (nonexistent_id >= lowest_id) then
    write(stderrunit,*) 'Underflow in iceberg ids!',nonexistent_id,lowest_id,mpp_pe()
  endif
  ! Sort the list "ids" (largest first)
  do j = 1, nbergs-1
    do i = j+1, nbergs
      if (ids(j) < ids(i)) then
        ! Swap
        id = ids(i)
        ids(i) = ids(j)
        ids(j) = id
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
      id = ids(j)
      have_berg = .true.
    else
      id = nonexistent_id
      have_berg = .false.
    endif
    lid = id
    call mpp_max(lid)
    if (have_berg .and. id == lid) then
      ii = 1 ! This berg is mine
      j = j + 1
    else
      ii = 0 ! This berg is not mine
    endif
    call mpp_sum(ii)
    if (ii > 1) then
      if (verbose) write(stderrunit,*) 'Duplicated berg across PEs with id=',id,lid,' seen',ii,' times pe=',mpp_pe(),k,j,nbergs
      check_for_duplicate_ids_in_list = check_for_duplicate_ids_in_list + 1
    elseif (ii == 0) then
      if (verbose) write(stderrunit,*) 'Berg not accounted for on all PEs with id=',id,lid,' seen',ii,' times pe=',mpp_pe(),k,j,nbergs
    endif
  enddo

end function check_for_duplicate_ids_in_list

!> Unit test for check_for_duplicate_ids_in_list()
subroutine test_check_for_duplicate_ids_in_list()
  ! Local variables
  integer :: k
  integer(kind=8), dimension(:), allocatable :: ids
  integer :: error_count

  allocate(ids(5))
  do k = 1,5
    ids(k) = k + 5*mpp_pe()
  enddo
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count /= 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('KID, test_check_for_duplicate_ids_in_list:', 'Unit test for clean list failed!', FATAL)
  endif
  if (mpp_pe() == mpp_root_pe()) ids(5) = ids(4)
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count == 0) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('KID, test_check_for_duplicate_ids_in_list:', 'Unit test for dirty list failed!', FATAL)
  endif
  if (mpp_pe() == mpp_root_pe()) ids(5) = 7 + 5*mpp_pe()
  error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.false.)
  call mpp_sum(error_count)
  if (error_count == 0 .and. mpp_npes()>1) then
    error_count = check_for_duplicate_ids_in_list(5, ids, verbose=.true.)
    call error_mesg('KID, test_check_for_duplicate_ids_in_list:', 'Unit test for a really dirty list failed!', FATAL)
  endif
  deallocate(ids)
end subroutine test_check_for_duplicate_ids_in_list

end module
