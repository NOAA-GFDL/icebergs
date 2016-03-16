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
use ice_bergs_framework, only: icebergs_gridded, xyt, iceberg, icebergs, buffer
use ice_bergs_framework, only: verbose, really_debug,debug,old_bug_rotated_weights,budget,use_roundoff_fix
use ice_bergs_framework, only: find_cell,find_cell_by_search,count_bergs,is_point_in_cell,pos_within_cell
use ice_bergs_framework, only: nclasses,old_bug_bilin
use ice_bergs_framework, only: sum_mass,sum_heat,bilin,yearday,count_bergs,bergs_chksum
use ice_bergs_framework, only: checksum_gridded,add_new_berg_to_list
use ice_bergs_framework, only: send_bergs_to_other_pes,move_trajectory,move_all_trajectories
use ice_bergs_framework, only: record_posn,check_position,print_berg,print_bergs,print_fld
use ice_bergs_framework, only: add_new_berg_to_list,delete_iceberg_from_list,destroy_iceberg
use ice_bergs_framework, only: grd_chksum2,grd_chksum3
use ice_bergs_framework, only: fix_restart_dates, offset_berg_dates
use ice_bergs_framework, only: orig_read  ! Remove when backward compatibility no longer needed

use ice_bergs_io,        only: ice_bergs_io_init,write_restart,write_trajectory
use ice_bergs_io,        only: read_restart_bergs,read_restart_bergs_orig,read_restart_calving

implicit none ; private

public icebergs_init, icebergs_end, icebergs_run, icebergs_stock_pe, icebergs
public icebergs_incr_mass, icebergs_save_restart

real, parameter :: pi_180=pi/180. ! Converts degrees to radians
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

integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "ice_bergs: "//trim(version)

  call ice_bergs_framework_init(bergs, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth=ocean_depth, maskmap=maskmap, fractional_area=fractional_area)

  call mpp_clock_begin(bergs%clock_ior)
  call ice_bergs_io_init(bergs,io_layout)
  if(orig_read) then
     call read_restart_bergs_orig(bergs,Time)
  else
     call read_restart_bergs(bergs,Time)
  endif
  call bergs_chksum(bergs, 'read_restart bergs')
  if (fix_restart_dates) call offset_berg_dates(bergs,Time)
  call read_restart_calving(bergs)
  call mpp_clock_end(bergs%clock_ior)

  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_init, initial status')


end subroutine icebergs_init


! ##############################################################################

subroutine interactive_force(bergs,berg,IA_x, IA_y, u0, v0, u1, v1, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) !Calculating interactive force between icebergs. Alon,  Markpoint_4
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
type(iceberg), pointer :: other_berg
real :: T1, L1, W1, lon1, lat1, x1, y1, R1, A1   !Current iceberg
real :: T2, L2, W2, lon2, lat2, x2, y2, R2, A2   !Other iceberg
real :: r_dist_x, r_dist_y, r_dist, A_o, trapped, T_min
real, intent(in) :: u0,v0, u1, v1
real :: P_11, P_12, P_21, P_22
real :: u2, v2
real :: Rearth
logical :: critical_interaction_damping_on
real :: spring_coef, accel_spring, radial_damping_coef, p_ia_coef, tangental_damping_coef
real, intent(out) :: IA_x, IA_y 
real, intent(out) :: P_ia_11, P_ia_12, P_ia_22, P_ia_21, P_ia_times_u_x, P_ia_times_u_y
integer :: stderrunit

Rearth=6360.e3
!spring_coef=1.e-4
spring_coef=bergs%spring_coef
radial_damping_coef=bergs%radial_damping_coef
tangental_damping_coef=bergs%tangental_damping_coef
critical_interaction_damping_on=bergs%critical_interaction_damping_on

!Using critical values for damping rather than manually setting the damping.
if (critical_interaction_damping_on) then
     radial_damping_coef=2.*sqrt(spring_coef)  ! Critical damping  
     tangental_damping_coef=(2.*sqrt(spring_coef)/5)  ! Critical damping  /5   (just a guess)
endif


! Get the stderr unit number.  Not sure what this does
  stderrunit = stderr()

IA_x=0.
IA_y=0.
P_ia_11=0. ; P_ia_12=0. ;  P_ia_21=0.;  P_ia_22=0. 
P_ia_times_u_x=0. ; P_ia_times_u_y=0.


L1=berg%length
W1=berg%width
T1=berg%thickness 
A1=L1*W1 
R1=sqrt(A1/pi) ! Interaction radius of the iceberg (assuming circular icebergs)
lon1=berg%lon; lat1=berg%lat
call rotpos_to_tang(lon1,lat1,x1,y1)

  other_berg=>bergs%first

!Note: This summing should be made order invarient. 
!Note: Need to limit how many icebergs we search over
  do while (associated(other_berg)) ! loop over all other bergs  
       L2=other_berg%length
       W2=other_berg%width
       T2=other_berg%thickness 
       u2=other_berg%uvel_old !Old values are used to make it order invariant 
       v2=other_berg%vvel_old !Old values are used to make it order invariant 
       A2=L2*W2
       R2=sqrt(A2/pi) ! Interaction radius of the other iceberg
       lon2=berg%lon_old; lat2=berg%lat_old !Old values are used to make it order invariant
       call rotpos_to_tang(lon2,lat2,x2,y2)

       r_dist_x=x1-x2 ; r_dist_y=y1-y2
       r_dist=sqrt( ((x1-x2)**2) + ((y1-y2)**2) )

      call overlap_area(R1,R2,r_dist,A_o,trapped)
      T_min=min(T1,T2)

      !Calculating spring force  (later this should only be done on the first time around)
      accel_spring=spring_coef*(T_min/T1)*(A_o/A1)
      if ((r_dist>0.) .AND. (r_dist< (R1+R2)) ) then
        IA_x=IA_x+(accel_spring*(r_dist_x/r_dist))
        IA_y=IA_y+(accel_spring*(r_dist_y/r_dist))


       !Working out the damping

       !Paralel velocity
        P_11=(r_dist_x*r_dist_x)/(r_dist**2)
        P_12=(r_dist_x*r_dist_y)/(r_dist**2)
        P_21=(r_dist_x*r_dist_y)/(r_dist**2)
        P_22=(r_dist_y*r_dist_y)/(r_dist**2)
        p_ia_coef=radial_damping_coef*(T_min/T1)*(A_o/A1)
        p_ia_coef=p_ia_coef*(0.5*(sqrt((((P_11*(u2-u1))+(P_12*(v2-v1)))**2)+ (((P_12*(u2-u1))+(P_22*(v2-v1)))**2))+sqrt((((P_11*(u2-u0))+(P_12*(v2-v0)))**2)+(((P_12*(u2-u0)) +(P_22*(v2-v0)))**2))))
        P_ia_11=P_ia_11+p_ia_coef*P_11
        P_ia_12=P_ia_12+p_ia_coef*P_12
        P_ia_21=P_ia_21+p_ia_coef*P_21
        P_ia_22=P_ia_22+p_ia_coef*P_22
        P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
        P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))


        !Normal velocities
        P_11=1-P_11  ;  P_12=-P_12 ; P_22=1-P_22
        p_ia_coef=tangental_damping_coef*(T_min/T1)*(A_o/A1)
        p_ia_coef=p_ia_coef*(0.5*(sqrt((((P_11*(u2-u1))+(P_12*(v2-v1)))**2)+ (((P_12*(u2-u1))+(P_22*(v2-v1)))**2))+sqrt((((P_11*(u2-u0))+(P_12*(v2-v0)))**2)+(((P_12*(u2-u0)) +(P_22*(v2-v0)))**2))))
        P_ia_11=P_ia_11+p_ia_coef*P_11
        P_ia_12=P_ia_12+p_ia_coef*P_12
        P_ia_21=P_ia_21+p_ia_coef*P_21
        P_ia_22=P_ia_22+p_ia_coef*P_22
        P_ia_times_u_x=P_ia_times_u_x+ (p_ia_coef* ((P_11*u2) +(P_12*v2)))
        P_ia_times_u_y=P_ia_times_u_y+ (p_ia_coef* ((P_12*u2) +(P_22*v2)))

!print *, 'P_11',P_11
!print *, 'P_21',P_21
!print *, 'P_12',P_12
!print *, 'P_22',P_22

    endif

       other_berg=>other_berg%next
  enddo ! loop over all bergs

  contains


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


  subroutine rotpos_to_tang(lon, lat, x, y)
  ! Arguments
  real, intent(in) :: lon, lat
  real, intent(out) :: x, y
  ! Local variables
  real :: r,colat,clon,slon

    if (lat>90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat>90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went very wrong!',FATAL)
    endif
    if (lat==90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat==90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went wrong!',FATAL)
    endif

    colat=90.-lat
    r=Rearth*(colat*pi_180)
    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    x=r*clon
    y=r*slon

  end subroutine rotpos_to_tang

end subroutine interactive_force


! ##############################################################################


subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, axn, ayn, bxn, byn, debug_flag) !Saving  acceleration for Verlet, Adding Verlet flag - Alon  MP1
!subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, debug_flag) !old version commmented out by Alon
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
integer, intent(in) :: i, j
real, intent(in) :: xi, yj, lat, uvel, vvel, uvel0, vvel0, dt
real, intent(inout) :: ax, ay
real, intent(inout) :: axn, ayn, bxn, byn ! Added implicit and explicit accelerations to output -Alon
logical, optional :: debug_flag
! Local variables
type(icebergs_gridded), pointer :: grd
real :: uo, vo, ui, vi, ua, va, uwave, vwave, ssh_x, ssh_y, sst, cn, hi
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
  call interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi)

   f_cori=(2.*omega)*sin(pi_180*lat)
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
  c_ice=rho_ice     /M*(0.5*Cd_iv*W*hi              )
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
           call Interactive_force(bergs, berg, IA_x, IA_y, uvel0, vvel0, uvel0, vvel0, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
           if (.not.Runge_not_Verlet) then
               axn=axn + IA_x
               ayn=ayn + IA_y
           else
               bxn=bxn + IA_x
               byn=byn + IA_y
           endif
!print *,'IA_x=',IA_x
!print *,'IA_x=',IA_x,'IA_y',IA_y
!print *,'P_ia_11',P_ia_11,'P_ia_12',P_ia_12, 'P_ia_21',P_ia_21,'P_ia_22', P_ia_22
!print *, 'P_ia_times_u_x', P_ia_times_u_x, 'P_ia_times_u_y', P_ia_times_u_y

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
 
  do itloop=1,2 ! Iterate on drag coefficients

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
        call Interactive_force(bergs, berg, IA_x, IA_y, us, vs, uvel0, vvel0, P_ia_11, P_ia_12, P_ia_21, P_ia_22, P_ia_times_u_x, P_ia_times_u_y) ! Spring forces, Made by Alon.
    endif
     if (beta>0.) then ! If implicit, use u_star, v_star rather than RK4 latest
         RHS_x=RHS_x -(((P_ia_11*u_star)+(P_ia_12*v_star))-P_ia_times_u_x) 
         RHS_y=RHS_y -(((P_ia_21*u_star)+(P_ia_22*v_star))-P_ia_times_u_y) 
    else
         RHS_x=RHS_x - (((P_ia_11*uvel)+(P_ia_12*vvel))-P_ia_times_u_x)
         RHS_y=RHS_y - (((P_ia_21*uvel)+(P_ia_22*vvel))-P_ia_times_u_y)       
    endif
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
            A11=A11+P_ia_11
            A12=A12+P_ia_12
            A21=A21+P_ia_21
            A22=A22+P_ia_22
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

  enddo ! itloop
 

!Saving the totally explicit part of the acceleration to use in finding the next position and u_star -Alon
    axn=0.
    ayn=0.
    if (.not.Runge_not_Verlet) then
        axn=-gravity*ssh_x +wave_rad*uwave + IA_x
        ayn=-gravity*ssh_y +wave_rad*vwave + IA_y
    endif
    if (C_N>0.) then !  C_N=1 for Crank Nicolson Coriolis, C_N=0 for full implicit Coriolis !Alon 
      axn=axn+f_cori*vveln
      ayn=ayn-f_cori*uveln
    endif
    bxn= ax-(axn/2) !Alon
    byn= ay-(ayn/2) !Alon

 
  ! Limit speed of bergs based on a CFL criteria
  if (bergs%speed_limit>0.) then
    speed=sqrt(uveln*uveln+vveln*vveln) ! Speed of berg
    if (speed>0.) then
      loc_dx=min(0.5*(grd%dx(i,j)+grd%dx(i,j-1)),0.5*(grd%dy(i,j)+grd%dy(i-1,j))) ! min(dx,dy)
     !new_speed=min(loc_dx/dt*bergs%speed_limit,speed) ! Restrict speed to dx/dt x factor
      new_speed=loc_dx/dt*bergs%speed_limit ! Speed limit as a factor of dx / dt 
      if (new_speed<speed) then
        uveln=uveln*(new_speed/speed) ! Scale velocity to reduce speed
        vveln=vveln*(new_speed/speed) ! without changing the direction
        bergs%nspeeding_tickets=bergs%nspeeding_tickets+1
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
integer :: i,j, stderrunit
type(iceberg), pointer :: this, next
real, parameter :: perday=1./86400.

  ! For convenience
  grd=>bergs%grd

  this=>bergs%first
  do while(associated(this))
    if (debug) call check_position(grd, this, 'thermodynamics (top)')

    call interp_flds(grd, this%ine, this%jne, this%xi, this%yj, this%uo, this%vo, &
            this%ui, this%vi, this%ua, this%va, this%ssh_x, this%ssh_y, this%sst, &
            this%cn, this%hi)
    SST=this%sst
    IC=min(1.,this%cn+bergs%sicn_shift) ! Shift sea-ice concentration 
    M=this%mass
    T=this%thickness ! total thickness
  ! D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
  ! F=T-D ! freeboard
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

    ! Rolling  - The corrected scheme has been included. The old scheme is here for legacy reasons
    Dn=(bergs%rho_bergs/rho_seawater)*Tn ! draught (keel depth)
    if ( Dn>0. ) then
      if (bergs%use_updated_rolling_scheme) then
        if (bergs%tip_parameter>0.) then
                tip_parameter=bergs%tip_parameter
        else
          ! Equation 27 from Burton et al 2012, or equivolently, Weeks and Mellor 1979 with constant density
          tip_parameter=sqrt(6*(bergs%rho_bergs/rho_seawater)*(1-(bergs%rho_bergs/rho_seawater)))   !using default values gives 0.92
        endif
        !print *, 'tip_parameter',tip_parameter
        if (Th<(tip_parameter* min(Wn,Ln)))  then     !note that we use the Thickness instead of the Draft
          if (Wn<Ln) then
            T=Tn
            Tn=Wn
            Wn=T
          else
            T=Tn
            Tn=Ln
            Ln=T
          endif
        endif
      else     
        !print *, 'using old tipping scheme'
        if ( max(Wn,Ln)<sqrt(0.92*(Dn**2)+58.32*Dn) ) then
          T=Tn
          Tn=Wn
          Wn=T
        end if
      end if
      Dn=(bergs%rho_bergs/rho_seawater)*Tn ! re-calculate draught (keel depth) for grounding
    endif

    ! Store the new state of iceberg (with L>W)
    this%mass=Mnew
    this%mass_of_bits=nMbits
    this%thickness=Tn
    this%width=min(Wn,Ln)
    this%length=max(Wn,Ln)

    next=>this%next

    ! Did berg completely melt?
    if (Mnew<=0.) then ! Delete the berg
      call move_trajectory(bergs, this)
      call delete_iceberg_from_list(bergs%first, this)
      bergs%nbergs_melted=bergs%nbergs_melted+1
    else ! Diagnose mass distribution on grid
      if (grd%id_virtual_area>0)&
           & grd%virtual_area(i,j)=grd%virtual_area(i,j)+(Wn*Ln+Abits)*this%mass_scaling ! m^2
      if (grd%id_mass>0 .or. bergs%add_weight_to_ocean)&
           & grd%mass(i,j)=grd%mass(i,j)+Mnew/grd%area(i,j)*this%mass_scaling ! kg/m2
      if (grd%id_bergy_mass>0 .or. bergs%add_weight_to_ocean)&
           & grd%bergy_mass(i,j)=grd%bergy_mass(i,j)+nMbits/grd%area(i,j)*this%mass_scaling ! kg/m2
      if (bergs%add_weight_to_ocean .and. .not. bergs%time_average_weight) then
        if (bergs%grounding_fraction>0.) then
          Hocean=bergs%grounding_fraction*(grd%ocean_depth(i,j)+grd%ssh(i,j))
          if (Dn>Hocean) Mnew=Mnew*min(1.,Hocean/Dn)
        endif
        call spread_mass_across_ocean_cells(grd, i, j, this%xi, this%yj, Mnew, nMbits, this%mass_scaling)
      endif
    endif
  
    this=>next
  enddo

end subroutine thermodynamics

! ##############################################################################

subroutine spread_mass_across_ocean_cells(grd, i, j, x, y, Mberg, Mbits, scaling)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  real, intent(in) :: x, y, Mberg, Mbits, scaling
  ! Local variables
  real :: xL, xC, xR, yD, yC, yU, Mass
  real :: yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR
  real, parameter :: rho_seawater=1035.
  
  Mass=(Mberg+Mbits)*scaling
  ! This line attempts to "clip" the weight felt by the ocean. The concept of
  ! clipping is non-physical and this step should be replaced by grounding.
  if (grd%clipping_depth>0.) Mass=min(Mass,grd%clipping_depth*grd%area(i,j)*rho_seawater)

  xL=min(0.5, max(0., 0.5-x))
  xR=min(0.5, max(0., x-0.5))
  xC=max(0., 1.-(xL+xR))
  yD=min(0.5, max(0., 0.5-y))
  yU=min(0.5, max(0., y-0.5))
  yC=max(0., 1.-(yD+yU))

  yDxL=yD*xL*grd%msk(i-1,j-1)
  yDxC=yD*xC*grd%msk(i  ,j-1)
  yDxR=yD*xR*grd%msk(i+1,j-1)
  yCxL=yC*xL*grd%msk(i-1,j  )
  yCxR=yC*xR*grd%msk(i+1,j  )
  yUxL=yU*xL*grd%msk(i-1,j+1)
  yUxC=yU*xC*grd%msk(i  ,j+1)
  yUxR=yU*xR*grd%msk(i+1,j+1)
  yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

  grd%mass_on_ocean(i,j,1)=grd%mass_on_ocean(i,j,1)+yDxL*Mass
  grd%mass_on_ocean(i,j,2)=grd%mass_on_ocean(i,j,2)+yDxC*Mass
  grd%mass_on_ocean(i,j,3)=grd%mass_on_ocean(i,j,3)+yDxR*Mass
  grd%mass_on_ocean(i,j,4)=grd%mass_on_ocean(i,j,4)+yCxL*Mass
  grd%mass_on_ocean(i,j,5)=grd%mass_on_ocean(i,j,5)+yCxC*Mass
  grd%mass_on_ocean(i,j,6)=grd%mass_on_ocean(i,j,6)+yCxR*Mass
  grd%mass_on_ocean(i,j,7)=grd%mass_on_ocean(i,j,7)+yUxL*Mass
  grd%mass_on_ocean(i,j,8)=grd%mass_on_ocean(i,j,8)+yUxC*Mass
  grd%mass_on_ocean(i,j,9)=grd%mass_on_ocean(i,j,9)+yUxR*Mass

end subroutine spread_mass_across_ocean_cells

! ##############################################################################

subroutine interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi)
! Arguments
type(icebergs_gridded), pointer :: grd
integer, intent(in) :: i, j
real, intent(in) :: xi, yj
real, intent(out) :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi
! Local variables
real :: cos_rot, sin_rot
#ifdef USE_OLD_SSH_GRADIENT
real :: dxm, dx0, dxp
#endif
real :: hxm, hxp
real, parameter :: ssh_coast=0.00

  cos_rot=bilin(grd, grd%cos, i, j, xi, yj) ! If true, uses the inverted bilin function
  sin_rot=bilin(grd, grd%sin, i, j, xi, yj)

  uo=bilin(grd, grd%uo, i, j, xi, yj)
  vo=bilin(grd, grd%vo, i, j, xi, yj)
  ui=bilin(grd, grd%ui, i, j, xi, yj)
  vi=bilin(grd, grd%vi, i, j, xi, yj)
  ua=bilin(grd, grd%ua, i, j, xi, yj)
  va=bilin(grd, grd%va, i, j, xi, yj)
  ! These fields are cell centered (A-grid) and would
  ! best be interpolated using PLM. For now we use PCM!
  sst=grd%sst(i,j) ! A-grid
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


! ##############################################################################

subroutine icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, calving_hflx, cn, hi, &
                        stagger, stress_stagger)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: time
real, dimension(:,:), intent(inout) :: calving, calving_hflx
real, dimension(:,:), intent(in) :: uo, vo, ui, vi, tauxa, tauya, ssh, sst, cn, hi
integer,    optional, intent(in) :: stagger, stress_stagger
! Local variables
integer :: iyr, imon, iday, ihr, imin, isec, k
type(icebergs_gridded), pointer :: grd
logical :: lerr, sample_traj, write_traj, lbudget, lverbose
real :: unused_calving, tmpsum, grdd_berg_mass, grdd_bergy_mass
integer :: i, j, Iu, ju, iv, Jv, Iu_off, ju_off, iv_off, Jv_off
real :: mask
real, dimension(:,:), allocatable :: uC_tmp, vC_tmp
integer :: vel_stagger, str_stagger

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
  grd%mass(:,:)=0.
  if (bergs%add_weight_to_ocean) grd%mass_on_ocean(:,:,:)=0.
  grd%virtual_area(:,:)=0.

  ! Manage time
  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  bergs%current_year=iyr
  bergs%current_yearday=yearday(imon, iday, ihr, imin, isec)
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if (bergs%traj_sample_hrs>0) then
     if (mod(24*iday+ihr,bergs%traj_sample_hrs).eq.0) sample_traj=.true.
  end if
  write_traj=.false.
  if (bergs%traj_write_hrs>0) then
     if (mod(24*iday+ihr,bergs%traj_write_hrs).eq.0) write_traj=.true.
  end if
  lverbose=.false.
  if (bergs%verbose_hrs>0) then
     if (mod(24*iday+ihr,bergs%verbose_hrs).eq.0) lverbose=verbose
  end if
  lbudget=.false.
  if (bergs%verbose_hrs>0) then
     if (mod(24*iday+ihr,bergs%verbose_hrs).eq.0) lbudget=budget
  end if
  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'diamonds: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', bergs%current_year, bergs%current_yearday

  ! Adapt calving flux from coupler for use here
 !call sanitize_field(grd%calving,1.e20)
  tmpsum=sum( calving(:,:)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_calving_received=bergs%net_calving_received+tmpsum*bergs%dt
  grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:) ! Units of kg/m2/s
  grd%calving(:,:)=grd%calving(:,:)*grd%msk(:,:)*grd%area(:,:) ! Convert to kg/s
  tmpsum=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_incoming_calving=bergs%net_incoming_calving+tmpsum*bergs%dt
  if (grd%id_calving>0) &
    lerr=send_data(grd%id_calving, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  ! Adapt calving heat flux from coupler
  grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)=calving_hflx(:,:) ! Units of W/m2
  grd%calving_hflx(:,:)=grd%calving_hflx(:,:)*grd%msk(:,:) ! Mask (just in case)
  if (grd%id_calving_hflx_in>0) &
    lerr=send_data(grd%id_calving_hflx_in, grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  tmpsum=sum( grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_incoming_calving_heat=bergs%net_incoming_calving_heat+tmpsum*bergs%dt ! Units of J

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
      grd%vi(I,J) = mask * 0.5*(vi(iv,Jv)+vo(iv+1,Jv))
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
    !   If the iceberg model used symmetric memory, the starting value of these
    ! copies would need to be decremented by 1.
    do I=grd%isc,grd%iec ; do j=grd%jsc,grd%jec
      uC_tmp(I,j) = tauxa(I+Iu_off, j+ju_off)
    enddo ; enddo
    do i=grd%isc,grd%iec ; do J=grd%jsc,grd%jec
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
  call mpp_update_domains(grd%ssh, grd%domain)
  grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec)=sst(:,:)-273.15 ! Note convert from Kelvin to Celsius
  call mpp_update_domains(grd%sst, grd%domain)
  ! Copy sea-ice concentration and thickness (resides on A grid)
  grd%cn(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=cn(:,:)
  call mpp_update_domains(grd%cn, grd%domain)
  grd%hi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=hi(:,:)
  call mpp_update_domains(grd%hi, grd%domain)

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
  call calve_icebergs(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (calved)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after calving')
  call mpp_clock_end(bergs%clock_cal)

  ! For each berg, evolve
  call mpp_clock_begin(bergs%clock_mom)
  if (associated(bergs%first)) call evolve_icebergs(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(bergs%grd, 's/r run after evolve')
  call mpp_clock_end(bergs%clock_mom)

  ! Send bergs to other PEs
  call mpp_clock_begin(bergs%clock_com)
  call send_bergs_to_other_pes(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (exchanged)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after exchange')
  call mpp_clock_end(bergs%clock_com)

  ! Iceberg thermodynamics (melting) + rolling
  call mpp_clock_begin(bergs%clock_the)
  if (associated(bergs%first)) call thermodynamics(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (thermo)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after thermodynamics')
  call mpp_clock_end(bergs%clock_the)

  ! For each berg, record
  call mpp_clock_begin(bergs%clock_dia)
  if (sample_traj.and.associated(bergs%first)) call record_posn(bergs)
  if (write_traj) then
    call move_all_trajectories(bergs)
    call write_trajectory(bergs%trajectories)
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
  if (grd%id_cn>0) &
    lerr=send_data(grd%id_cn, grd%cn(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_hi>0) &
    lerr=send_data(grd%id_hi, grd%hi(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_floating_melt>0) &
    lerr=send_data(grd%id_floating_melt, grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
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
  if (grd%id_mass>0) &
    lerr=send_data(grd%id_mass, grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_stored_ice>0) &
    lerr=send_data(grd%id_stored_ice, grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)
  if (grd%id_real_calving>0) &
    lerr=send_data(grd%id_real_calving, grd%real_calving(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)
  if (grd%id_ssh>0) &
    lerr=send_data(grd%id_ssh, grd%ssh(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_fax>0) &
    lerr=send_data(grd%id_fax, tauxa(:,:), Time)
  if (grd%id_fay>0) &
    lerr=send_data(grd%id_fay, tauya(:,:), Time)


  ! Dump icebergs to screen
  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_run, status')
  call mpp_clock_end(bergs%clock_dia)

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
    bergs%stored_heat_end=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%floating_mass_end=sum_mass(bergs%first)
    bergs%icebergs_mass_end=sum_mass(bergs%first,justbergs=.true.)
    bergs%bergy_mass_end=sum_mass(bergs%first,justbits=.true.)
    bergs%floating_heat_end=sum_heat(bergs%first)
    grd%tmpc(:,:)=0.;
    call mpp_clock_end(bergs%clock); call mpp_clock_end(bergs%clock_dia) ! To enable calling of public s/r
    call icebergs_incr_mass(bergs, grd%tmpc)
    call mpp_clock_begin(bergs%clock_dia); call mpp_clock_begin(bergs%clock) ! To enable calling of public s/r
    bergs%returned_mass_on_ocean=sum( grd%tmpc(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%nbergs_end=count_bergs(bergs)
    call mpp_sum(bergs%stored_end)
    call mpp_sum(bergs%stored_heat_end)
    call mpp_sum(bergs%floating_mass_end)
    call mpp_sum(bergs%icebergs_mass_end)
    call mpp_sum(bergs%bergy_mass_end)
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
    if (mpp_pe().eq.mpp_root_pe()) then
 100 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
 200 format("diamonds: ",a19,10(a18,"=",es14.7,x,a2,:,","))
      call report_state('stored ice','kg','',bergs%stored_start,'',bergs%stored_end,'')
      call report_state('floating','kg','',bergs%floating_mass_start,'',bergs%floating_mass_end,'',bergs%nbergs_end)
      call report_state('icebergs','kg','',bergs%icebergs_mass_start,'',bergs%icebergs_mass_end,'')
      call report_state('bits','kg','',bergs%bergy_mass_start,'',bergs%bergy_mass_end,'')
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

end subroutine icebergs_run

! ##############################################################################

subroutine icebergs_incr_mass(bergs, mass, Time)
! Arguments
type(icebergs), pointer :: bergs
real, dimension(bergs%grd%isc:bergs%grd%iec,bergs%grd%jsc:bergs%grd%jec), intent(inout) :: mass
type(time_type), intent(in), optional :: Time
! Local variables
integer :: i, j
type(icebergs_gridded), pointer :: grd
real :: dmda
logical :: lerr

  if (.not. associated(bergs)) return

  if (.not. bergs%add_weight_to_ocean) return

  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_int)

  ! For convenience
  grd=>bergs%grd

  ! Add iceberg+bits mass field to non-haloed SIS field (kg/m^2)
 !mass(:,:)=mass(:,:)+( grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec) &
 !                    + grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec) )


  if (debug) then
    grd%tmp(:,:)=0.; grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=mass
    call grd_chksum2(grd, grd%tmp, 'mass in (incr)')
  endif

  call mpp_update_domains(grd%mass_on_ocean, grd%domain)
  if (.not. old_bug_rotated_weights) then
    do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
      if (grd%parity_x(i,j)<0.) then
        ! This block assumes both parity_x and parity_y are negative
        ! (i.e. a 180 degree rotation). In general, we should handle
        ! +/- 90 degree rotations as well but in CM2*-class models
        ! this is not necessary. -aja
        dmda=grd%mass_on_ocean(i,j,9); grd%mass_on_ocean(i,j,9)=grd%mass_on_ocean(i,j,1); grd%mass_on_ocean(i,j,1)=dmda
        dmda=grd%mass_on_ocean(i,j,8); grd%mass_on_ocean(i,j,8)=grd%mass_on_ocean(i,j,2); grd%mass_on_ocean(i,j,2)=dmda
        dmda=grd%mass_on_ocean(i,j,7); grd%mass_on_ocean(i,j,7)=grd%mass_on_ocean(i,j,3); grd%mass_on_ocean(i,j,3)=dmda
        dmda=grd%mass_on_ocean(i,j,6); grd%mass_on_ocean(i,j,6)=grd%mass_on_ocean(i,j,4); grd%mass_on_ocean(i,j,4)=dmda
      endif
    enddo; enddo
  endif
  do j=grd%jsc, grd%jec; do i=grd%isc, grd%iec
    dmda=grd%mass_on_ocean(i,j,5) &
         + ( ( (grd%mass_on_ocean(i-1,j-1,9)+grd%mass_on_ocean(i+1,j+1,1))   &
         +     (grd%mass_on_ocean(i+1,j-1,7)+grd%mass_on_ocean(i-1,j+1,3)) ) &
         +   ( (grd%mass_on_ocean(i-1,j  ,6)+grd%mass_on_ocean(i+1,j  ,4))   &
         +     (grd%mass_on_ocean(i  ,j-1,8)+grd%mass_on_ocean(i  ,j+1,2)) ) )
    if (grd%area(i,j)>0) dmda=dmda/grd%area(i,j)*grd%msk(i,j)
    if (.not. bergs%passive_mode) mass(i,j)=mass(i,j)+dmda
    if (grd%id_mass_on_ocn>0) grd%tmp(i,j)=dmda
  enddo; enddo
  if (present(Time).and. (grd%id_mass_on_ocn>0)) &
    lerr=send_data(grd%id_mass_on_ocn, grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  if (debug) then
    grd%tmp(:,:)=0.; grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=mass
    call grd_chksum3(grd, grd%mass_on_ocean, 'mass bergs (incr)')
    call grd_chksum2(grd, grd%tmp, 'mass out (incr)')
  endif

  call mpp_clock_end(bergs%clock_int)
  call mpp_clock_end(bergs%clock)

end subroutine icebergs_incr_mass

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
type(iceberg) :: newberg
logical :: lret
real :: xi, yj, ddt, calving_to_bergs, calved_to_berg, heat_to_bergs, heat_to_berg
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

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
          newberg%start_day=bergs%current_yearday+ddt/86400.
          newberg%start_mass=bergs%initial_mass(k)
          newberg%mass_scaling=bergs%mass_scaling(k)
          newberg%mass_of_bits=0.
          newberg%heat_density=grd%stored_heat(i,j)/grd%stored_ice(i,j,k) ! This is in J/kg
          call add_new_berg_to_list(bergs%first, newberg)
          calved_to_berg=bergs%initial_mass(k)*bergs%mass_scaling(k) ! Units of kg
          ! Heat content
          heat_to_berg=calved_to_berg*newberg%heat_density ! Units of J
          grd%stored_heat(i,j)=grd%stored_heat(i,j)-heat_to_berg
          heat_to_bergs=heat_to_bergs+heat_to_berg
          ! Stored mass
          grd%stored_ice(i,j,k)=grd%stored_ice(i,j,k)-calved_to_berg
          calving_to_bergs=calving_to_bergs+calved_to_berg
          grd%real_calving(i,j,k)=grd%real_calving(i,j,k)+calved_to_berg/bergs%dt
          ddt=ddt+bergs%dt*2./17. ! Minor offset to start day
          icnt=icnt+1
          bergs%nbergs_calved=bergs%nbergs_calved+1
          bergs%nbergs_calved_by_class(k)=bergs%nbergs_calved_by_class(k)+1
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
real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1, ax1, ay1, axn1, ayn1 
real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2, ax2, ay2, axn2, ayn2
real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3, ax3, ay3, axn3, ayn3
real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4, ax4, ay4, axn4, ayn4
real :: uveln, vveln, lonn, latn, un, vn, dxdln
real :: x1, xdot1, xddot1, y1, ydot1, yddot1, xddot1n, yddot1n 
real :: x2, xdot2, xddot2, y2, ydot2, yddot2, xddot2n, yddot2n
real :: x3, xdot3, xddot3, y3, ydot3, yddot3, xddot3n, yddot3n
real :: x4, xdot4, xddot4, y4, ydot4, yddot4, xddot4n, yddot4n
real :: xn, xdotn, yn, ydotn, xddotn, yddotn
real :: bxddot, byddot                                               ! Added by Alon
real :: axn, ayn, bxn, byn                                           ! Added by Alon - explicit and implicit accelations from the previous step
real :: r180_pi, dt, dt_2, dt_6, dydl, Rearth
integer :: i, j
integer :: i1,j1,i2,j2,i3,j3,i4,j4
real :: xi, yj
logical :: bounced, on_tangential_plane, error_flag
logical :: Runge_not_Verlet  ! Runge_not_Verlet=1 for Runge Kutta, =0 for Verlet method. Added by Alon
type(iceberg), pointer :: berg
integer :: stderrunit
logical :: interactive_icebergs_on  ! Flag to decide whether to use forces between icebergs.

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

  interactive_icebergs_on=bergs%interactive_icebergs_on  ! Loading directly from namelist/default , Alon

  ! Common constants
  r180_pi=1./pi_180
  dt=bergs%dt
  dt_2=0.5*dt
  dt_6=dt/6.
  Rearth=6360.e3

  !Choosing time stepping scheme - Alon
   !Runge_not_Verlet=.False.    !Loading manually: true=Runge Kutta, False=Verlet   , Alon
   Runge_not_Verlet=bergs%Runge_not_Verlet  ! Loading directly from namelist/default , Alon

  berg=>bergs%first
  do while (associated(berg)) ! loop over all bergs

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

  i=berg%ine
  j=berg%jne
  xi=berg%xi
  yj=berg%yj
  bounced=.false.
  on_tangential_plane=.false.
  if (berg%lat>89.) on_tangential_plane=.true.
  i1=i;j1=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)



  if (Runge_not_Verlet) then !Start of the Runge-Kutta Loop -Added by Alon, MP2

 !Loading past acceleartions - Alon
  axn=berg%axn; ayn=berg%ayn !Alon
  axn1=axn; axn2=axn; axn3=axn; axn4=axn
  ayn1=ayn; ayn2=ayn; ayn3=ayn; ayn4=ayn



 ! A1 = A(X1)
  lon1=berg%lon; lat1=berg%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)
  dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  dydl=r180_pi/Rearth
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
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag)
  i2=i; j2=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon2,lat2,x2,y2)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon2, lat2, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
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
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 2!',FATAL)
  endif
  dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
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
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag)
  i3=i; j3=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon3,lat3,x3,y3)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon3, lat3, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
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
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 3!',FATAL)
  endif
  dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
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
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag)
  i4=i; j4=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon4,lat4,x4,y4)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon4, lat4, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
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
  dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
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
  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag)
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)

  if (.not.error_flag) then
    if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
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
 
  endif ! End of the Runge-Kutta Loop -added by Alon  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (.not.Runge_not_Verlet) then !Start of the Verlet time_stepping -Whole loop added by Alon

 ! In this scheme a_n and b_n are saved from the previous timestep, giving the explicit and implicit parts of the acceleration, and a_np1, b_np1 are for the next time step
 ! Note that ax1=a_np1/2 +b_np1, as calculated by the acceleration subrouting
 ! Positions and velocity is updated by
 ! X2 = X1+dt*V1+((dt^2)/2)*a_n +((dt^2)/2)*b_n = X1+dt*u_star +((dt^2)/2)*b_n 
 ! V2 = V1+dt/2*a_n +dt/2*a_np1 +dt*b_n = u_star + dt/2*a_np1 + dt*b_np1 = u_star +dt*ax

!print *, 'you are here!'

lon1=berg%lon; lat1=berg%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)
  dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  dydl=r180_pi/Rearth
  uvel1=berg%uvel; vvel1=berg%vvel

!Loading past acceleartions - Alon
  axn=berg%axn; ayn=berg%ayn !Alon
  bxn=berg%bxn; byn=berg%byn !Alon

print *, 'first', axn, bxn, lon1, lat1, uvel1, i, j ,xi, yj

! Velocities used to update the position
  uvel2=uvel1+(dt_2*axn)+(dt_2*bxn)                    !Alon
  vvel2=vvel1+(dt_2*ayn)+(dt_2*byn)                    !Alon

if (on_tangential_plane) call rotvec_to_tang(lon1,uvel2,vvel2,xdot2,ydot2)
  u2=uvel2*dxdl1; v2=vvel2*dydl


!Solving for new position
  if (on_tangential_plane) then
    xn=x1+(dt*xdot2) ; yn=y1+(dt*ydot2)             !Alon
    call rotpos_from_tang(xn,yn,lonn,latn)
  else
    lonn=lon1+(dt*u2) ; latn=lat1+(dt*v2)  !Alon
  endif
  dxdln=r180_pi/(Rearth*cos(latn*pi_180))

! Turn the velocities into u_star, v_star.(uvel3 is v_star)
  uvel3=uvel1+(dt_2*axn)                  !Alon
  vvel3=vvel1+(dt_2*ayn)                  !Alon


!Adjusting mass...                      Alon decided to move this before calculating the new velocities (so that acceleration can be a fn(r_np1)
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lonn, latn, uvel3, vvel3, i, j, xi, yj, bounced, error_flag)  !Alon:"unclear which velocity to use here?"
!  call adjust_index_and_ground(grd, lonn, latn, uvel1, vvel1, i, j, xi, yj, bounced, error_flag)  !Alon:"unclear which velocity to use here?"

if (bounced) then  !This is the case when the iceberg changes direction due to  topography
  axn=0.
  ayn=0.
  bxn=0.
  byn=0.
endif


  i2=i; j2=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)


print *, 'second', axn, bxn, lon1, lat1, uvel1, i , j , xi, yj
!Calling the acceleration   (note that the velocity is converted to u_star inside the accel script)
  call accel(bergs, berg, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, ax1, ay1, axn, ayn, bxn, byn) !axn, ayn, bxn, byn - Added by Alon

print *, 'third', axn, bxn, lon1, lat1, uvel1, i, j, xi, yj
!Solving for the new velocity
  if (on_tangential_plane) then
    call rotvec_to_tang(lonn,uvel3,vvel3,xdot3,ydot3)
    call rotvec_to_tang(lon1,ax1,ay1,xddot1,yddot1)
    xdotn=xdot3+(dt*xddot1); ydotn=ydot3+(dt*yddot1)                                    !Alon
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
  else
    uvel4=uvel3+(dt*ax1); vvel4=vvel3+(dt*ay1)    !Alon , we call it uvel3, vvel3 until it is put into lat/long co-ordinates, where it becomes uveln, vveln
  endif
!  uveln=uvel4*dxdln; vveln=vvel4*dydl    !Converted to degrees.  (Perhaps this should not be here)
  uveln=uvel4
  vveln=vvel4 

print *, 'forth', axn, bxn, lon1, lat1, uvel1, i, j, xi, yj, uveln
!Debugging
  if (.not.error_flag) then
    if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
   call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i=',i1,i2,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j=',j1,j2,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lonn=',lon1,lonn,berg%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,latn=',lat1,latn,berg%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u3,un,u0=',uvel3,uveln,berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v3,vn,v0=',vvel3,vveln,berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1=',&
        & dt*ax1
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1=',&
        & dt*ay1
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u3,un,u0=',&
        & dt*uvel3,dt*uveln,dt*berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v3,vn,v0=',&
        & dt*vvel3,dt*vveln,dt*berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u_n (deg)=',&
        & dt*u1,dt*uveln
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v_n (deg)=',&
        & dt*v1,dt*vveln
   write(stderrunit,*) 'diamonds, evolve_iceberg: on_tangential_plane=',on_tangential_plane
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(bergs, berg, i2, j2, xi, yj, latn, uvel3, vvel3, uvel1, vvel1, dt_2, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn, bxn, byn - Added by Alon
   
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


print *, 'fifth', axn, bxn, lon1, lat1, uvel1, i, j, xi, yj, uveln

  endif ! End of the Verlet Stepiing -added by Alon  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Saving all the iceberg variables.
  berg%axn=axn !Alon
  berg%ayn=ayn !Alon
  berg%bxn=bxn !Alon
  berg%byn=byn !Alon
  berg%lon=lonn
  berg%lat=latn
  berg%uvel=uveln
  berg%vvel=vveln
  berg%ine=i
  berg%jne=j
  berg%xi=xi
  berg%yj=yj
  !call interp_flds(grd, i, j, xi, yj, berg%uo, berg%vo, berg%ui, berg%vi, berg%ua, berg%va, berg%ssh_x, berg%ssh_y, berg%sst)
  !if (debug) call print_berg(stderr(), berg, 'evolve_iceberg, final posn.')
  if (debug) call check_position(grd, berg, 'evolve_iceberg (bot)')


  berg=>berg%next
  enddo ! loop over all bergs

! When we are using interactive icebergs, we update the (old) iceberg positions and velocities in a second loop, all together (to make code order invarient)
 if (interactive_icebergs_on) then
  berg=>bergs%first
  do while (associated(berg)) ! loop over all bergs

      !Updating iceberg positions and velocities
      berg%lon_old=berg%lon
      berg%lat_old=berg%lat
      berg%uvel_old=berg%uvel
      berg%vvel_old=berg%vvel

      berg=>berg%next
  enddo ! loop over all bergs


 endif
  contains

  subroutine rotpos_to_tang(lon, lat, x, y)
  ! Arguments
  real, intent(in) :: lon, lat
  real, intent(out) :: x, y
  ! Local variables
  real :: r,colat,clon,slon

    if (lat>90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat>90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went very wrong!',FATAL)
    endif
    if (lat==90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat==90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went wrong!',FATAL)
    endif

    colat=90.-lat
    r=Rearth*(colat*pi_180)
    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    x=r*clon
    y=r*slon

  end subroutine rotpos_to_tang

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

subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(inout) :: lon, lat, uvel, vvel, xi, yj
integer, intent(inout) :: i,j
logical, intent(out) :: bounced, error
! Local variables
logical lret, lpos
real, parameter :: posn_eps=0.05
integer :: icount, i0, j0, inm, jnm
real :: xi0, yj0, lon0, lat0

  bounced=.false.
  error=.false.
  lon0=lon; lat0=lat ! original position
  i0=i; j0=j ! original i,j
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  xi0=xi; yj0=yj ! original xi,yj
  if (debug) then
    !Sanity check lret, xi and yj
    lret=is_point_in_cell(grd, lon, lat, i, j)
    if (xi<0. .or. xi>1. .or. yj<0. .or. yj>1.) then
      if (lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
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
        lret=is_point_in_cell(grd, lon, lat, i, j, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
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
      if (inm<grd%ied) then
        inm=inm+1
      endif
    endif
    if (yj.lt.0.) then
      if (jnm>grd%jsd) then
        jnm=jnm-1
      endif
    elseif (yj.gt.1.) then
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
    elseif (xi.gt.1.) then
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
    elseif (yj.gt.1.) then
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
      if (xi>1.) xi=1.-posn_eps
      if (xi<0.) xi=posn_eps
      if (yj>1.) yj=1.-posn_eps
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
                                       ! OR that it did not move at all (round-off problem)
    if (debug) then
      write(stderrunit,*) 'diamonds, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'diamonds, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'diamonds, adjust: i0,j0=',i0,j0 
      write(stderrunit,*) 'diamonds, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'diamonds, adjust: xi,yj=',xi,yj 
      write(stderrunit,*) 'diamonds, adjust: i,j=',i,j
      write(stderrunit,*) 'diamonds, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'diamonds, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
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
    else
      call error_mesg('diamonds, adjust', 'Berg iterated many times without bouncing!', WARNING)
    endif
  endif
  if (xi>1.) xi=1.-posn_eps
  if (xi<0.) xi=posn_eps
  if (yj>1.) yj=1.-posn_eps
  if (yj<0.) yj=posn_eps
  lon=bilin(grd, grd%lon, i, j, xi, yj)
  lat=bilin(grd, grd%lat, i, j, xi, yj)
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  if (.not. lret) then
    write(stderrunit,*) 'diamonds, adjust: Should not get here! Berg is not in cell after adjustment'
    if (debug) error=.true.
  endif

end subroutine adjust_index_and_ground

end subroutine evolve_icebergs



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
    berg_mass=sum_mass(bergs%first)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=stored_mass+berg_mass

  case (ISTOCK_HEAT)
    berg_mass=sum_mass(bergs%first)
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

  call write_trajectory(bergs%trajectories)

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
  deallocate(bergs%grd%virtual_area)
  deallocate(bergs%grd%mass)
  deallocate(bergs%grd%mass_on_ocean)
  deallocate(bergs%grd%tmp)
  deallocate(bergs%grd%tmpc)
  deallocate(bergs%grd%stored_ice)
  deallocate(bergs%grd%real_calving)
  deallocate(bergs%grd%uo)
  deallocate(bergs%grd%vo)
  deallocate(bergs%grd%ui)
  deallocate(bergs%grd%vi)
  deallocate(bergs%grd%ua)
  deallocate(bergs%grd%va)
  deallocate(bergs%grd%ssh)
  deallocate(bergs%grd%sst)
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
