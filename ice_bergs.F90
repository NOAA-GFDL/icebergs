module ice_bergs

use constants_mod, only: radius, pi, omega, HLF
use fms_mod, only: open_namelist_file, check_nml_error, close_file
use fms_mod, only: field_exist, read_data, get_global_att_value
use fms_mod, only: stdlog, stdout, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: write_version_number, read_data, write_data, file_exist
use mosaic_mod, only: get_mosaic_ntiles, get_mosaic_ncontacts
use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_sum, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id, CLOCK_COMPONENT
use fms_mod, only: clock_flag_default
use mpp_domains_mod, only: domain2D, mpp_update_domains, mpp_define_domains
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use diag_manager_mod, only: diag_axis_init


implicit none ; private

include 'netcdf.inc'

public icebergs_init, icebergs_end, icebergs_run, icebergs_stock_pe

type :: icebergs_gridded
  type(domain2D), pointer :: domain ! MPP domain
  integer :: halo ! Nominal halo width
  integer :: isc, iec, jsc, jec ! Indices of computational domain
  integer :: isd, ied, jsd, jed ! Indices of data domain
  integer :: my_pe, pe_N, pe_S, pe_E, pe_W ! MPI PE identifiers
  real, dimension(:,:), pointer :: lon=>NULL() ! Longitude of cell corners
  real, dimension(:,:), pointer :: lat=>NULL() ! Latitude of cell corners
  real, dimension(:,:), pointer :: dx=>NULL() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: dy=>NULL() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: area=>NULL() ! Area of cell (m^2)
  real, dimension(:,:), pointer :: msk=>NULL() ! Ocean-land mask (1=ocean)
  real, dimension(:,:), pointer :: cos=>NULL() ! Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: sin=>NULL() ! Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: uo=>NULL() ! Ocean zonal flow (m/s)
  real, dimension(:,:), pointer :: vo=>NULL() ! Ocean meridional flow (m/s)
  real, dimension(:,:), pointer :: ui=>NULL() ! Ice zonal flow (m/s)
  real, dimension(:,:), pointer :: vi=>NULL() ! Ice meridional flow (m/s)
  real, dimension(:,:), pointer :: ua=>NULL() ! Atmosphere zonal flow (m/s)
  real, dimension(:,:), pointer :: va=>NULL() ! Atmosphere meridional flow (m/s)
  real, dimension(:,:), pointer :: ssh=>NULL() ! Sea surface height (m)
  real, dimension(:,:), pointer :: sst=>NULL() ! Sea surface temperature (oC)
  real, dimension(:,:), pointer :: cn=>NULL() ! Sea-ice concentration (0 to 1)
  real, dimension(:,:), pointer :: hi=>NULL() ! Sea-ice thickness (m)
  real, dimension(:,:), pointer :: calving=>NULL() ! Calving mass rate [frozen runoff] (kg/s)
  real, dimension(:,:), pointer :: melt=>NULL() ! Iceberg melting mass rate (kg/s/m^2)
  real, dimension(:,:), pointer :: mass=>NULL() ! Mass distribution (kg/m^2)
  real, dimension(:,:), pointer :: tmp=>NULL() ! Temporary work space
  real, dimension(:,:,:), pointer :: stored_ice=>NULL() ! Accumulated ice mass flux at calving locations (kg)
  ! Diagnostics handles
  integer :: id_uo=-1, id_vo=-1, id_calving=-1, id_stored_ice=-1, id_accum=-1, id_unused=-1, id_melt=-1
  integer :: id_mass=-1, id_ui=-1, id_vi=-1, id_ua=-1, id_va=-1, id_sst=-1
end type icebergs_gridded

type :: xyt
  real :: lon, lat, day
  real :: mass, depth, width, length, uvel, vvel
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst
  integer :: year
  type(xyt), pointer :: next=>NULL()
end type xyt

type :: iceberg
  type(iceberg), pointer :: prev=>NULL(), next=>NULL()
  ! State variables (specific to the iceberg, needed for restarts)
  real :: lon, lat, uvel, vvel, mass, depth, width, length
  real :: start_lon, start_lat, start_day, start_mass, mass_scaling
  integer :: start_year
  integer :: ine,jne ! nearest index in NE direction (for convenience)
  real :: xi, yj ! Non-dimensional coords within current cell (0..1)
  ! Environment variables (as seen by the iceberg)
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst
  type(xyt), pointer :: trajectory=>NULL()
end type iceberg

type :: buffer
  integer :: size=0
  real, dimension(:,:), pointer :: data
end type buffer

type, public :: icebergs ; private
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: first=>NULL()
  type(xyt), pointer :: trajectories=>NULL()
  real :: dt           ! Time-step between iceberg calls (should make adaptive?)
  integer :: current_year
  real :: current_yearday ! 1.00-365.99
  integer :: traj_sample_hrs
  integer :: verbose_hrs
  integer :: clock ! id for fms timers
  real, dimension(:), pointer :: initial_mass, distribution, mass_scaling
  real, dimension(:), pointer :: initial_depth, initial_width, initial_length
  logical :: restarted=.false. ! Indicate whether we read state from a restart or not
  type(buffer), pointer :: obuffer_e, ibuffer_e, obuffer_w, ibuffer_w
end type icebergs

! Global constants
integer, parameter :: nclasses=10 ! Number of ice bergs classes
integer, parameter :: file_format_major_version=0
integer, parameter :: file_format_minor_version=1
character(len=*), parameter :: version = '$Id: ice_bergs.F90,v 1.1.2.7 2008/06/05 17:01:56 aja Exp $'
character(len=*), parameter :: tagname = '$Name:  $'
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

! Global data (minimal for debugging)
logical :: verbose=.false. ! Be verbose to stderr
logical :: budget=.true. ! Calculate budgets
logical :: debug=.false. ! Turn on debugging

contains

! ##############################################################################

subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, ax, ay)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
integer, intent(in) :: i, j
real, intent(in) :: xi, yj, lat, uvel, vvel
real, intent(inout) :: ax, ay
! Local variables
type(icebergs_gridded), pointer :: grd
real :: uo, vo, ui, vi, ua, va, uw, vw, ssh_x, ssh_y, sst, hi
real :: f_cori, D, W, L, M, F, drag_ocn, drag_atm, drag_ice, wave_rad
real :: a2, dva, Cr, Lwavelength, Lcutoff, Ltop
real, parameter :: alpha=0.0, beta=1.0, accel_lim=1.e-3, Cr0=0.06, vel_lim=5.
real :: lambda, detA, A11, A12, axe, aye
logical :: dumpit

  ! For convenience
  grd=>bergs%grd

  ! Interpolate gridded fields to berg 
  call interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, hi)

  f_cori=(2.*omega)*sin(pi_180*lat)

  M=berg%mass
  D=berg%depth
  W=berg%width
  L=berg%length
  F=D/5.

  drag_ocn=rho_seawater/M*(0.5*Cd_wv*W*(D-hi)+Cd_wh*W*L)*sqrt( (uvel-uo)**2+(vvel-vo)**2 )
  drag_atm=rho_air     /M*(0.5*Cd_av*W*F     +Cd_ah*W*L)*sqrt( (uvel-ua)**2+(vvel-va)**2 )
  drag_ice=rho_ice     /M*(0.5*Cd_iv*W*hi              )*sqrt( (uvel-ui)**2+(vvel-vi)**2 )
  if (abs(ui)+abs(vi).eq.0.) drag_ice=0.

  uw=ua-uo; vw=va-vo  ! Use wind speed rel. to ocean for wave model (aja)?
  dva=sqrt(uw*uw+vw*vw)
  a2=(0.5*0.02025*(dva**2))**2
  Lwavelength=0.32*(dva**2) ! Surface wave length fitted to data in table at
  !      http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
  Lcutoff=0.125*Lwavelength
  Ltop=0.25*Lwavelength
  Cr=Cr0*min(max(0.,(L-Lcutoff)/(Ltop-Lcutoff)),1.) ! Wave radiation coefficient
  !     fitted to graph from Carrieres et al.,  POAC Drift Model.
  wave_rad=0.5*rho_seawater/M*Cr*gravity*a2*sqrt(W*L)

  ! Explicit accelerations
  dva=sqrt((uw-uvel)**2+(vw-vvel)**2)
  if (dva.ne.0.) then
    wave_rad=wave_rad/dva
  else
    wave_rad=0.
  endif
  axe= f_cori*vvel -gravity*ssh_x +wave_rad*(uw-uvel) &
      -drag_ocn*(uvel-uo) -drag_atm*(uvel-ua) -drag_ice*(uvel-ui)
  aye=-f_cori*uvel -gravity*ssh_y +wave_rad*(vw-vvel) &
      -drag_ocn*(vvel-vo) -drag_atm*(vvel-va) -drag_ice*(vvel-vi)

  ! Solve for implicit accelerations
  if (alpha+beta.gt.0.) then
    lambda=drag_ocn+drag_atm+drag_ice  +wave_rad
    A11=1.+beta*bergs%dt*lambda
    A12=alpha*bergs%dt*f_cori
    detA=1./(A11**2+A12**2)
    ax=detA*(A11*axe+A12*aye)
    ay=detA*(A11*aye-A12*axe)
  else
    ax=axe; ay=aye
  endif
  
  dumpit=.false.; if (abs(uvel)>vel_lim.or.abs(vvel)>vel_lim) dumpit=.true.
  if (dumpit.or.abs(ax)>accel_lim) write(stderr(),'(a,i3,9(1xa,1pe12.3))') 'Large ax: pe=',mpp_pe(), &
      'f*v=',f_cori*vvel, &
      'g*H_x=',gravity*ssh_x, &
      'wave*ua=',wave_rad*(uw-uvel), &
      'd*(u-uo)=',drag_ocn*(uvel-uo), &
      'd*(u-ua)=',drag_atm*(uvel-ua), &
      'd*(u-ui)=',drag_ice*(uvel-ui)
  if (dumpit.or.abs(ay)>accel_lim) write(stderr(),'(a,i3,9(1xa,1pe12.3))') 'Large ay: pe=',mpp_pe(), &
      'f*u=',f_cori*uvel, &
      'g*H_y=',gravity*ssh_y, &
      'wave*va=',wave_rad*(vw-vvel), &
      'd*(v-vo)=',drag_ocn*(vvel-vo), &
      'd*(v-va)=',drag_atm*(vvel-va), &
      'd*(v-vi)=',drag_ice*(vvel-vi)
  if (dumpit.or.abs(ax)>accel_lim.or.abs(ay)>accel_lim) then
    write(stderr(),'(a,i3,9(1xa,1pe12.3))') '          pe=',mpp_pe(), &
      'dva=',dva, &
      'wave=',wave_rad, &
      'xi=',xi, &
      'yj=',yj
    write(stderr(),'(a,i3,9(1xa,1pe12.3))') '          pe=',mpp_pe(), &
      'a=',sqrt(a2), &
      'Lo=',Lwavelength, &
      'Cr=',Cr
    write(stderr(),'(a,i3,9(1xa,1pe12.3))') '          pe=',mpp_pe(), &
      'M=',M, &
      'D=',D, &
      'W=',W, &
      'L=',L
    write(stderr(),'(a,i3,9(1xa,1pe12.3))') '          pe=',mpp_pe(), &
      'ax=',ax, &
      'ay=',ay, &
      'Dvo=',sqrt((uvel-uo)**2+(vvel-vo)**2), &
      'Dva=',sqrt((uvel-ua)**2+(vvel-va)**2), &
      'Dvi=',sqrt((uvel-ui)**2+(vvel-vi)**2), &
      'U=',sqrt(uvel**2+vvel**2)
    write(stderr(),'(a,i3,9(1xa,1pe12.3))') '          pe=',mpp_pe(), &
      'axe=',axe, &
      'aye=',aye
    call print_berg(stderr(),berg,'diamond, accel, large accel')
  endif

 !write(stderr(),'(i3,a,8(1pe10.2,1x,a))') mpp_pe(),'diamond, accel: drag_ocn=',drag_ocn,'dt*drag=',bergs%dt*drag_ocn,"rho/M=",rho_seawater/M*(0.5*Cd_wv*W*D+Cd_wh*W*L),"|u|=",sqrt( (uvel-uo)**2+(vvel-vo)**2 ),'f=',f_cori
 !write(stderr(),'(i3,a,8(1pe10.2,1x,a))') mpp_pe(),'diamond, accel: M=',M,'W=',W,"rho=",rho_seawater,"Cd_wv=",Cd_wv*W*D
end subroutine accel

! ##############################################################################

subroutine thermodynamics(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: M, D, W, L, SST, F, Vol, Ln, Wn, Dn, nVol, IC
real :: Mv, Me, Mb, melt, Mnew, dvo, dva, dM, Ss
integer :: i,j
type(iceberg), pointer :: this, next

  ! For convenience
  grd=>bergs%grd

  this=>bergs%first
  do while(associated(this))

    call interp_flds(grd, this%ine, this%jne, this%xi, this%yj, this%uo, this%vo, &
            this%ui, this%vi, this%ua, this%va, this%ssh_x, this%ssh_y, this%sst, IC)
    SST=this%sst
    M=this%mass
    D=this%depth
    F=0.2*D
    W=this%width
    L=this%length
    i=this%ine
    j=this%jne
    Vol=D*W*L

    ! Melt rates in m/day
    dvo=sqrt((this%uvel-this%uo)**2+(this%vvel-this%vo)**2)
    dva=sqrt((this%uvel-this%ua)**2+(this%vvel-this%va)**2)
    Ss=1.5*(dva**0.5)+0.1*dva ! Sea state
    Mb=max( 0.58*(dvo**0.8)*(SST+4.0)/(L**0.2), 0.) ! Basal turbulent melting
    Mv=max( 7.62e-3*SST+1.29e-3*(SST**2), 0.) ! Buoyant convection at sides
    Me=max( 1./12.*(SST+2.)*Ss*(1+cos(pi*(IC**3))) ,0.)  ! Wave erosion

    ! Sum and convert to kg/s
  ! melt=rho_ice/86400.*( (2.*(L+W)*D*Mv +F*(W+L)*Me) + W*L*Mb )
  ! dM=bergs%dt*melt ! Change in mass per time-step (kg)
  ! dM=min( dM, M ) ! Limit mass change by available mass
  ! Mnew=M-dM ! New mass
  ! melt=dM/bergs%dt ! kg/s
  ! D=D*(Mnew/M) ! Modify volume in proportion to change in mass
  ! this%depth=D

    ! Convert melt rates from m/day to m
    Ln=max(L-(Mv+Me)*(bergs%dt/86400.),0.)
    Wn=max(W-(Mv+Me)*(bergs%dt/86400.),0.)
    Dn=max(D-Mb*(bergs%dt/86400.),0.)
    nVol=Dn*Wn*Ln
    Mnew=(nVol/Vol)*M
    dM=M-Mnew

    this%mass=Mnew
    this%depth=Dn
    this%width=Wn
    this%length=Ln

    ! Rolling
    if ( Dn>0. .and. max(Wn,Ln)<sqrt(0.92*(Dn**2)+58.32*Dn) ) then
      this%depth=Wn
      this%width=Dn
    endif

    ! Add melting to the grid
    if (grd%area(i,j).ne.0.) then
      melt=dM/bergs%dt ! kg/s
      grd%melt(i,j)=grd%melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
    else
      write(stderr(),*) 'diamond, thermodynamics: berg appears to have grounded!!!! PE=',mpp_pe(),i,j
      call print_berg(stderr(),this,'thermodynamics, grounded')
      if (associated(this%trajectory)) &
        write(stderr(),*) 'traj=',this%trajectory%lon,this%trajectory%lat
      write(stderr(),*) 'msk=',grd%msk(i,j),grd%area(i,j)
      call error_mesg('diamond, thermodynamics', 'berg appears to have grounded!', FATAL)
    endif

    next=>this%next

    ! Did berg completely melt?
    if (Mnew<=0.) then ! Delete the berg
      call move_trajectory(bergs, this)
      call delete_iceberg_from_list(bergs%first, this)
    else ! Diagnose mass distribution on grid
      if(grd%id_mass>0) grd%mass(i,j)=grd%mass(i,j)+Mnew/grd%area(i,j)
    endif
  
    this=>next
  enddo

end subroutine thermodynamics

! ##############################################################################

subroutine interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi)
! Arguments
type(icebergs_gridded), pointer :: grd
integer, intent(in) :: i, j
real, intent(in) :: xi, yj
real, intent(out) :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst
real, intent(out), optional :: cn, hi
! Local variables
real :: cos_rot, sin_rot
real :: dxm, dx0, dxp, hxm, hxp
real, parameter :: ssh_coast=0.00

  cos_rot=bilin(grd, grd%cos, i, j, xi, yj)
  sin_rot=bilin(grd, grd%sin, i, j, xi, yj)

  uo=bilin(grd, grd%uo, i, j, xi, yj)
  vo=bilin(grd, grd%vo, i, j, xi, yj)
  ui=bilin(grd, grd%ui, i, j, xi, yj)
  vi=bilin(grd, grd%vi, i, j, xi, yj)
  ua=bilin(grd, grd%ua, i, j, xi, yj)
  va=bilin(grd, grd%va, i, j, xi, yj)
  sst=grd%sst(i,j) ! A-grid, instead of sst=bilin(grd, grd%sst, i, j, xi, yj)
  if (present(cn)) then
     cn=grd%cn(i,j) ! A-grid
  endif
  if (present(hi)) then
     hi=grd%hi(i,j) ! A-grid
  endif

  ! Estimate SSH gradient in X direction
  dxp=0.5*(grd%dx(i+1,j)+grd%dx(i+1,j-1))
  dx0=0.5*(grd%dx(i,j)+grd%dx(i,j-1))
  dxm=0.5*(grd%dx(i-1,j)+grd%dx(i-1,j-1))
  hxm=2.*(grd%ssh(i,j)-grd%ssh(i-1,j))/(dx0+dxm)*grd%msk(i-1,j) &
        +(-ssh_coast)/(dx0+dxm)*(1.-grd%msk(i-1,j)) ! force to drive bergs away from coasts
  hxp=2.*(grd%ssh(i+1,j)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i+1,j) &
        +(+ssh_coast)/(dx0+dxp)*(1.-grd%msk(i+1,j)) ! force to drive bergs away from coasts
  ssh_x=xi*hxp+(1.-xi)*hxm

  ! Estimate SSH gradient in Y direction
  dxp=0.5*(grd%dy(i,j+1)+grd%dy(i-1,j+1))
  dx0=0.5*(grd%dy(i,j)+grd%dy(i-1,j))
  dxm=0.5*(grd%dy(i,j-1)+grd%dy(i-1,j-1))
  hxm=2.*(grd%ssh(i,j)-grd%ssh(i,j-1))/(dx0+dxm)*grd%msk(i,j-1) &
        +(-ssh_coast)/(dx0+dxm)*(1.-grd%msk(i,j-1)) ! force to drive bergs away from coasts
  hxp=2.*(grd%ssh(i,j+1)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i,j+1) &
        +(+ssh_coast)/(dx0+dxp)*(1.-grd%msk(i,j+1)) ! force to drive bergs away from coasts
  ssh_y=yj*hxp+(1.-yj)*hxm
  
  ! Rotate vectors from ocean/ice grid to lat/lon coordinates
  call rotate(uo,vo)
  call rotate(ui,vi)
  call rotate(ua,va)
  call rotate(ssh_x,ssh_y)
 !u=uo; uo=cos_rot*uo+sin_rot*vo; vo=-sin_rot*u+cos_rot*vo
 !u=ui; ui=cos_rot*ui+sin_rot*vi; vi=-sin_rot*u+cos_rot*vi
 !u=ua; ua=cos_rot*ua+sin_rot*va; va=-sin_rot*u+cos_rot*va
 !u=ssh_x; ua=cos_rot*ssh_x+sin_rot*ssh_y; ssh_y=-sin_rot*u+cos_rot*ssh_y

  contains

  subroutine rotate(u,v)
  ! Arguments
  real, intent(inout) :: u,v
  ! Local variables
  real :: u_old, v_old

    u_old=u
    v_old=v
    u=cos_rot*u_old+sin_rot*v_old
    v=cos_rot*v_old-sin_rot*u_old

  end subroutine rotate

end subroutine interp_flds

! ##############################################################################

real function bilin(grd, fld, i, j, xi, yj)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
! Local variables

  bilin=(fld(i,j  )*(1.-xi)+fld(i-1,j  )*xi)*(1.-yj) &
       +(fld(i,j-1)*(1.-xi)+fld(i-1,j-1)*xi)*yj

end function bilin

! ##############################################################################

subroutine icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, cn, hi)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: time
real, dimension(:,:), intent(inout) :: calving
real, dimension(:,:), intent(in) :: uo, vo, ui, vi, tauxa, tauya, ssh, sst, cn, hi
! Local variables
integer :: iyr, imon, iday, ihr, imin, isec, nbergs
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
logical :: lerr, sample_traj, lbudget, lverbose
real :: incoming_calving, unused_calving, stored_mass, total_iceberg_mass, meltmass

  call mpp_clock_begin(bergs%clock)

  ! For convenience
  grd=>bergs%grd

  ! Manage time
  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  bergs%current_year=iyr
  bergs%current_yearday=float(imon-1)*31.+float(iday-1)+(float(ihr)+(float(imin)+float(isec)/60.)/60.)/24.
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if (bergs%traj_sample_hrs>0 .and. mod(24*iday+ihr,bergs%traj_sample_hrs).eq.0) sample_traj=.true.
  lverbose=.false.
  if (bergs%verbose_hrs>0 .and. mod(24*iday+ihr,bergs%verbose_hrs).eq.0) lverbose=verbose
  lbudget=.false.
  if (bergs%verbose_hrs>0 .and. mod(24*iday+ihr,bergs%verbose_hrs).eq.0) lbudget=budget
  if (mpp_pe()==mpp_root_pe().and.lverbose) write(stderr(),'(a,3i5,a,3i5,a,i5,f8.3)') &
       'diamond: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', bergs%current_year, bergs%current_yearday

  ! Adapt calving flux from coupler for use here
  grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:)
 !call sanitize_field(grd%calving,1.e20)
  grd%calving(:,:)=grd%calving(:,:)*grd%msk(:,:)*grd%area(:,:) ! Convert to kg/s from kg/m2/s
  if (lbudget) incoming_calving=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  if (grd%id_calving>0) &
    lerr=send_data(grd%id_calving, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  ! Copy ocean flow (resides on B grid)
  grd%uo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=uo(:,:)
  grd%vo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=vo(:,:)
  call mpp_update_domains(grd%uo, grd%vo, grd%domain)
  ! Copy ice flow (resides on B grid)
  grd%ui(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=ui(:,:)
  grd%vi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=vi(:,:)
  call mpp_update_domains(grd%ui, grd%vi, grd%domain)
  ! Copy atmospheric stress (resides on A grid)
  grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec)=tauxa(:,:) ! Note rough conversion from stress to speed
  grd%va(grd%isc:grd%iec,grd%jsc:grd%jec)=tauya(:,:) ! Note rough conversion from stress to speed
  call invert_tau_for_du(grd%ua, grd%va) ! Note rough conversion from stress to speed
 !grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec)=sign(sqrt(abs(tauxa(:,:))/0.01),tauxa(:,:))  ! Note rough conversion from stress to speed
 !grd%va(grd%isc:grd%iec,grd%jsc:grd%jec)=sign(sqrt(abs(tauya(:,:))/0.01),tauya(:,:))  ! Note rough conversion from stress to speed
  call mpp_update_domains(grd%ua, grd%va, grd%domain)
  ! Copy sea surface height and temperature(resides on A grid)
  grd%ssh(grd%isc:grd%iec,grd%jsc:grd%jec)=ssh(:,:)
  call mpp_update_domains(grd%ssh, grd%domain)
  grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec)=sst(:,:)-273.15 ! Note convert from Kelvin to Celsius
  call mpp_update_domains(grd%sst, grd%domain)
  ! Copy sea-ice concentration and thickness (resides on A grid)
  grd%cn(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=cn(:,:)
  call mpp_update_domains(grd%cn, grd%domain)
  grd%hi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=hi(:,:)
  call mpp_update_domains(grd%hi, grd%domain)

  ! Accumulate ice from calving
  call accumulate_calving(bergs)
  if (grd%id_accum>0) then
    grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:)
    grd%tmp(:,:)=grd%tmp(:,:)*grd%msk(:,:)*grd%area(:,:)-grd%calving(:,:)
    lerr=send_data(grd%id_accum, grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  endif
  if (grd%id_unused>0) &
    lerr=send_data(grd%id_unused, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (lbudget) unused_calving=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )

  ! Calve excess stored ice into icebergs
  call calve_icebergs(bergs)

  ! For each berg, evolve
  if (associated(bergs%first)) then
    this=>bergs%first
    do while (associated(this))
      call evolve_iceberg(bergs,this)
      this=>this%next
    enddo
  endif

  ! Send bergs to other PEs
  call send_bergs_to_other_pes(bergs)

  ! Ice berg thermodynamics (melting) + rolling
  grd%melt(:,:)=0.
  grd%mass(:,:)=0.
  if (associated(bergs%first)) call thermodynamics(bergs)

  ! For each berg, record
  if (sample_traj.and.associated(bergs%first)) then
    this=>bergs%first
    do while (associated(this))
      call record_posn(this, bergs%current_year ,bergs%current_yearday)
      this=>this%next
    enddo
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
  if (grd%id_melt>0) &
    lerr=send_data(grd%id_melt, grd%melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_mass>0) &
    lerr=send_data(grd%id_mass, grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_stored_ice>0) &
    lerr=send_data(grd%id_stored_ice, grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)

  ! Dump icebergs to screen
 !if (lverbose) call print_bergs(stderr(),bergs,'icebergs_run, status')

  ! Diagnose budgets
  if (lbudget) then
    call mpp_sum(incoming_calving)
    call mpp_sum(unused_calving)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    call mpp_sum(stored_mass)
    total_iceberg_mass=sum_icebergs_mass(bergs%first)
    call mpp_sum(total_iceberg_mass)
    nbergs=count_bergs(bergs)
    call mpp_sum(nbergs)
    meltmass=sum( grd%melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(meltmass)
    if (mpp_pe().eq.mpp_root_pe()) write(stderr(),'(a,5(1pe10.4,a),i6)') &
        'diamond, budget: incoming calving=',incoming_calving, &
        ' kg/s, unused calving=',unused_calving, &
        ' kg/s, melt=',meltmass, &
        ' kg/s, stored mass=',stored_mass, &
        ' kg, ice berg mass=',total_iceberg_mass, &
        ' kg, # of bergs=',nbergs
  endif

  ! Return what ever calving we did not use and additional icebergs melt
  where (grd%area(:,:)>0.)
    calving(:,:)=grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)/grd%area(:,:) &
                 +grd%melt(grd%isc:grd%iec,grd%jsc:grd%jec)
  elsewhere
    calving(:,:)=0.
  end where

  call mpp_clock_end(bergs%clock)
end subroutine icebergs_run

! ##############################################################################

subroutine accumulate_calving(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: remaining_dist
integer :: k
logical, save :: first_call=.true.

  ! For convenience
  grd=>bergs%grd

  ! This is a hack to simplify initialization
  if (first_call.and..not.bergs%restarted) then
    first_call=.false.
    do k=1, nclasses
      where (grd%calving==0.) grd%stored_ice(:,:,k)=0.
    enddo
  endif

  remaining_dist=1.
  do k=1, nclasses
    grd%stored_ice(:,:,k)=grd%stored_ice(:,:,k)+bergs%dt*grd%calving(:,:)*bergs%distribution(k)
    remaining_dist=remaining_dist-bergs%distribution(k)
  enddo
  if (remaining_dist.lt.0.) then
    write(stderr(),*) 'diamond, accumulate_calving: sum(distribution)>1!!!',remaining_dist
    call error_mesg('diamond, accumulate_calving', 'calving is OVER distributed!', WARNING)
  endif
  ! Remove the calving accounted for by accumulation
  grd%calving(:,:)=grd%calving(:,:)*remaining_dist

end subroutine accumulate_calving

! ##############################################################################

subroutine calve_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
integer :: i,j,k
type(iceberg) :: newberg
logical :: lret
real :: xi, yj

  ! For convenience
  grd=>bergs%grd

  do k=1, nclasses
    do j=grd%jsc, grd%jec
      do i=grd%isc, grd%iec
        if (grd%stored_ice(i,j,k).ge.bergs%initial_mass(k)*bergs%mass_scaling(k)) then
          newberg%lon=0.25*((grd%lon(i,j)+grd%lon(i-1,j-1))+(grd%lon(i-1,j)+grd%lon(i-1,j)))
          newberg%lat=0.25*((grd%lat(i,j)+grd%lat(i-1,j-1))+(grd%lat(i-1,j)+grd%lat(i-1,j)))
         !write(stderr(),*) 'diamond, calve_icebergs: creating new iceberg at ',newberg%lon,newberg%lat
          lret=pos_within_cell(grd, newberg%lon, newberg%lat, i, j, xi, yj)
          if (.not.lret) then
            write(stderr(),*) 'diamond, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('diamond, calve_icebergs', 'berg is not in the correct cell!', FATAL)
          endif
          newberg%ine=i
          newberg%jne=j
          newberg%xi=xi
          newberg%yj=yj
          newberg%uvel=0.
          newberg%vvel=0.
          newberg%mass=bergs%initial_mass(k)
          newberg%depth=bergs%initial_depth(k)
          newberg%width=bergs%initial_width(k)
          newberg%length=bergs%initial_length(k)
          newberg%start_lon=grd%lon(i,j)
          newberg%start_lat=grd%lat(i,j)
          newberg%start_year=bergs%current_year
          newberg%start_day=bergs%current_yearday
          newberg%start_mass=bergs%initial_mass(k)
          newberg%mass_scaling=bergs%mass_scaling(k)
          call add_new_berg_to_list(bergs%first, newberg)
         !if (verbose) call print_berg(stderr(), bergs%first, 'calve_icebergs, new berg created')
          grd%stored_ice(i,j,k)=grd%stored_ice(i,j,k)-bergs%initial_mass(k)*bergs%mass_scaling(k)
        endif
      enddo
    enddo
  enddo

end subroutine calve_icebergs

! ##############################################################################

subroutine evolve_iceberg(bergs, berg)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
! Local variables
type(icebergs_gridded), pointer :: grd
real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1, ax1, ay1
real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2, ax2, ay2
real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3, ax3, ay3
real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4, ax4, ay4
real :: uveln, vveln, lonn, latn
real :: r180_pi, dt, dt_2, dt_6, dydl, Rearth
integer :: i, j, i0, j0
real :: xi, yj, xi0, yj0
logical :: bounced

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
  
  ! For convenience
  grd=>bergs%grd

  if (.not. is_point_in_cell(bergs%grd, berg%lon, berg%lat, berg%ine, berg%jne) ) then
    call print_berg(stderr(), berg, 'evolve_iceberg, berg is not in proper starting cell')
    write(stderr(),'(a,i3,2(i4,3f8.2))') 'evolve_iceberg: berg is not in proper starting cell!!!!!', &
             mpp_pe(), &
             berg%ine,berg%lon,grd%lon(berg%ine-1,berg%jne-1),grd%lon(berg%ine,berg%jne), &
             berg%jne,berg%lat,grd%lon(berg%ine-1,berg%jne-1),grd%lon(berg%ine,berg%jne)
  endif

  r180_pi=1./pi_180
  dt=bergs%dt
  dt_2=0.5*dt
  dt_6=dt/6.
  Rearth=6360.e3
  i=berg%ine
  j=berg%jne
  xi=berg%xi
  yj=berg%yj

  i0=i; j0=j; xi0=xi; yj0=yj;
  
  ! A1 = A(X1)
  lon1=berg%lon
  lat1=berg%lat
  dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  dydl=r180_pi/Rearth
  uvel1=berg%uvel
  vvel1=berg%vvel
  u1=uvel1*dxdl1
  v1=vvel1*dydl
  call accel(bergs, berg, i, j, xi, yj, lat1, uvel1, vvel1, ax1, ay1)
  
  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  lon2=lon1+dt_2*u1
  lat2=lat1+dt_2*v1
  uvel2=uvel1+dt_2*ax1
  vvel2=vvel1+dt_2*ay1
  bounced=.false.
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced)
  if (.not. is_point_in_cell(bergs%grd, lon2, lat2, i, j) ) then
    call print_berg(stderr(), berg, 'evolve_iceberg, out of position at 2')
    write(stderr(),'(a,i3,a,2i3,4f8.3)') 'pe=',mpp_pe(),'pos2 lon,lat,i,j,xi,yj=',i,j,lon2,lat2,xi,yj
    write(stderr(),'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos2 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
  endif
  dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
  u2=uvel2*dxdl2
  v2=vvel2*dydl
  call accel(bergs, berg, i, j, xi, yj, lat2, uvel2, vvel2, ax2, ay2)
  
  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  lon3=lon1+dt_2*u2
  lat3=lat1+dt_2*v2
  uvel3=uvel1+dt_2*ax2
  vvel3=vvel1+dt_2*ay2
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced)
  if (.not. is_point_in_cell(bergs%grd, lon3, lat3, i, j) ) then
    call print_berg(stderr(), berg, 'evolve_iceberg, out of position at 3')
    write(stderr(),'(a,i3,a,2i3,4f8.3)') 'pe=',mpp_pe(),'pos3 lon,lat,i,j,xi,yj=',i,j,lon3,lat3,xi,yj
    write(stderr(),'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos3 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
  endif
  dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
  u3=uvel3*dxdl3
  v3=vvel3*dydl
  call accel(bergs, berg, i, j, xi, yj, lat3, uvel3, vvel3, ax3, ay3)
  
  !  X4 = X1+dt*V3 ; V4 = V1+dt*A3; A4=A(X4)
  lon4=lon1+dt*u3
  lat4=lat1+dt*v3
  uvel4=uvel1+dt*ax3
  vvel4=vvel1+dt*ay3
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced)
  if (.not. is_point_in_cell(bergs%grd, lon4, lat4, i, j) ) then
    call print_berg(stderr(), berg, 'evolve_iceberg, out of position at 4')
    write(stderr(),'(a,i3,a,2i3,4f8.3)') 'pe=',mpp_pe(),'pos4 lon,lat,i,j,xi,yj=',i,j,lon4,lat4,xi,yj
    write(stderr(),'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos4 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
  endif
  dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
  u4=uvel4*dxdl4
  v4=vvel4*dydl
  call accel(bergs, berg, i, j, xi, yj, lat4, uvel4, vvel4, ax4, ay4)
  
  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
  lonn=berg%lon+dt_6*( (u1+u4)+2.*(u2+u3) )
  latn=berg%lat+dt_6*( (v1+v4)+2.*(v2+v3) )
  uveln=berg%uvel+dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
  vveln=berg%vvel+dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
  
 !if (verbose) call print_berg(stderr(), berg, 'evolve_iceberg, initial pos')


  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced)

  if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j) ) then
    call print_berg(stderr(), berg, 'evolve_iceberg, out of cell at end!')
    write(stderr(),'(a,i3,a,2i3,4f8.3)') 'pe=',mpp_pe(),'posn lon,lat,i,j,xi,yj=',i,j,lonn,latn,xi,yj
    write(stderr(),'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
  endif

  berg%lon=lonn
  berg%lat=latn
  berg%uvel=uveln
  berg%vvel=vveln
  berg%ine=i
  berg%jne=j
  berg%xi=xi
  berg%yj=yj
 !call interp_flds(grd, i, j, xi, yj, berg%uo, berg%vo, berg%ui, berg%vi, berg%ua, berg%va, berg%ssh_x, berg%ssh_y, berg%sst)

 !if (verbose) call print_berg(stderr(), berg, 'evolve_iceberg, final posn.')

 contains

! ##############################################################################

subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(inout) :: lon, lat, uvel, vvel, xi, yj
integer, intent(inout) :: i,j
logical, intent(out) :: bounced
! Local variables
logical lret
real, parameter :: posn_eps=1e-3, sanity_acc=1.e-11
integer :: icount

  bounced=.false.
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  lret=is_point_in_cell(grd, lon, lat, i, j)
  ! Sanity check lret, xi and yj
  if (lret.and. (abs(xi-0.5)-0.5>sanity_acc.or.abs(yj-0.5)-0.5>sanity_acc)) then
    write(stderr(),*) 'diamond, adjust: WARNING!!! lret=T but |xi,yj|>1',lret,xi,yj,lon,lat,i,j
!   stop 'This should not ever happen!'
  elseif ((.not.lret).and. (abs(xi-0.5)-0.5<sanity_acc.and.abs(yj-0.5)-0.5<sanity_acc)) then
    write(stderr(),*) 'diamond, adjust: WARNING!!! lret=F but |xi,yj|<1',lret,xi,yj,lon,lat,i,j
!   stop 'This should not ever happen!'
  endif
  icount=0
 !do while(.not.lret.and.icount<=grd%halo-1) ! Note -1 needed because using non-symmetric arrays
  do while( .not.lret.and. icount<10 .and. &
      (i.gt.grd%isd.and.i.lt.grd%ied.and. &
       j.gt.grd%jsd.and.j.lt.grd%jed) )
  icount=icount+1
  if (xi.lt.0.) then
    if (grd%msk(i-1,j)>0.) then
      i=i-1
    else
     !write(stderr(),'(a,6f8.3,i)') 'diamond, adjust: bouncing berg from west',lon,lat,xi,yj,uvel,vvel,mpp_pe()
      xi=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
      bounced=.true.
    endif
  elseif (xi.gt.1.) then
    if (grd%msk(i+1,j)>0.) then
      i=i+1
    else
     !write(stderr(),'(a,6f8.3,i)') 'diamond, adjust: bouncing berg from east',lon,lat,xi,yj,uvel,vvel,mpp_pe()
      xi=1.-posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
      bounced=.true.
    endif
  endif
  if (yj.lt.0.) then
    if (grd%msk(i,j-1)>0.) then
      j=j-1
    else
     !write(stderr(),'(a,6f8.3,i)') 'diamond, adjust: bouncing berg from south',lon,lat,xi,yj,uvel,vvel,mpp_pe()
      yj=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
      bounced=.true.
    endif
  elseif (yj.gt.1.) then
    if (grd%msk(i,j+1)>0.) then
      j=j+1
    else
     !write(stderr(),'(a,6f8.3,i)') 'diamond, adjust: bouncing berg from north',lon,lat,xi,yj,uvel,vvel,mpp_pe()
      yj=1.-posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
      bounced=.true.
    endif
  endif
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  enddo

  if (.not.lret) &
     write(stderr(),'(a,2i4,6f8.3,2i3)') 'diamond, adjust: FAILED!!!',i,j,lon,lat,xi,yj,uvel,vvel,icount,mpp_pe()
  if (icount>grd%halo) &
     write(stderr(),'(a,2i4,6f8.3,2i3)') 'diamond, adjust: Large icount!!!',i,j,lon,lat,xi,yj,uvel,vvel,icount,mpp_pe()

end subroutine adjust_index_and_ground

end subroutine evolve_iceberg

! ##############################################################################

subroutine send_bergs_to_other_pes(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: kick_the_bucket, this
integer :: nbergs_to_send_e, nbergs_to_send_w, nbergs_incoming, i
integer, parameter :: max_bergs_in_buffer=50, buffer_width=16
real, dimension(buffer_width,max_bergs_in_buffer) :: obuffer_e, ibuffer_e
real, dimension(buffer_width,max_bergs_in_buffer) :: obuffer_w, ibuffer_w
type(icebergs_gridded), pointer :: grd

  ! For convenience
  grd=>bergs%grd

  ! Find number of bergs that headed east/west
  nbergs_to_send_e=0
  nbergs_to_send_w=0
  if (associated(bergs%first)) then
    this=>bergs%first
    do while (associated(this))
      if (this%ine.gt.bergs%grd%iec) then
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_e=nbergs_to_send_e+1
        call pack_berg_into_buffer(kick_the_bucket, obuffer_e, nbergs_to_send_e)
       !call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_e, nbergs_to_send_e)
        call move_trajectory(bergs, kick_the_bucket)
        call delete_iceberg_from_list(bergs%first,kick_the_bucket)
       !if (verbose) write(stderr(),*) 'diamond, send_bergs_to_other_pes: packing berg for E',mpp_pe(),associated(bergs%first)
      elseif (this%ine.lt.bergs%grd%isc) then
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_w=nbergs_to_send_w+1
        call pack_berg_into_buffer(kick_the_bucket, obuffer_w, nbergs_to_send_w)
       !call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_w, nbergs_to_send_w)
        call move_trajectory(bergs, kick_the_bucket)
        call delete_iceberg_from_list(bergs%first,kick_the_bucket)
       !if (verbose) write(stderr(),*) 'diamond, send_bergs_to_other_pes: packing berg for W',mpp_pe(),associated(bergs%first)
      else
        this=>this%next
      endif
    enddo
  endif

  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_e, plen=1, to_pe=grd%pe_E)
   !if (nbergs_to_send_e.gt.0) then
   !  write(stderr(),*) 'pe=',mpp_pe(),' sent the number',nbergs_to_send_e,' to',grd%pe_E,' (E)'
   !endif
    if (nbergs_to_send_e.gt.0) then
      call mpp_send(obuffer_e, nbergs_to_send_e*buffer_width, grd%pe_E)
     !if (verbose) write(stderr(),'(a,i3,a,i3,a)') 'diamond, send_bergs_to_other_pes:',mpp_pe(),' sending',nbergs_to_send_e,' bergs E'
    endif
  endif

  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_w, plen=1, to_pe=grd%pe_W)
   !if (nbergs_to_send_w.gt.0) then
   !  write(stderr(),*) 'pe=',mpp_pe(),' sent the number',nbergs_to_send_w,' to',grd%pe_W,' (W)'
   !endif
    if (nbergs_to_send_w.gt.0) then
      call mpp_send(obuffer_w, nbergs_to_send_w*buffer_width, grd%pe_W)
     !if (verbose) write(stderr(),'(a,i3,a,i3,a)') 'diamond, send_bergs_to_other_pes:',mpp_pe(),' sending',nbergs_to_send_w,' bergs W'
    endif
  endif

  if (grd%pe_W.ne.NULL_PE) then
    nbergs_incoming=-999
    call mpp_recv(nbergs_incoming, glen=1, from_pe=grd%pe_W)
    if (nbergs_incoming.lt.0) then
      write(stderr(),*) 'pe=',mpp_pe(),' received a bad number',nbergs_incoming,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_incoming.gt.0) then
     !write(stderr(),*) 'pe=',mpp_pe(),' received the number',nbergs_incoming,' from',grd%pe_W,' (W)'
      call mpp_recv(ibuffer_w, nbergs_incoming*buffer_width, grd%pe_W)
      do i=1, nbergs_incoming
        call unpack_berg_from_buffer(bergs%first, ibuffer_w, i)
      enddo
     !if (verbose) write(stderr(),*) 'diamond, send_bergs_to_other_pes: receiving bergs from W',mpp_pe()
    endif
  endif

  if (grd%pe_E.ne.NULL_PE) then
    nbergs_incoming=-999
    call mpp_recv(nbergs_incoming, glen=1, from_pe=grd%pe_E)
    if (nbergs_incoming.lt.0) then
      write(stderr(),*) 'pe=',mpp_pe(),' received a bad number',nbergs_incoming,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_incoming.gt.0) then
     !write(stderr(),*) 'pe=',mpp_pe(),' received the number',nbergs_incoming,' from',grd%pe_E,' (E)'
      call mpp_recv(ibuffer_e, nbergs_incoming*buffer_width, grd%pe_E)
      do i=1, nbergs_incoming
        call unpack_berg_from_buffer(bergs%first, ibuffer_e, i)
      enddo
     !if (verbose) write(stderr(),*) 'diamond, send_bergs_to_other_pes: receiving bergs from E',mpp_pe()
    endif
  endif

  call mpp_sync_self()

contains

  subroutine pack_berg_into_buffer2(berg, buff, n)
  ! Arguments
  type(iceberg), pointer :: berg
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  ! Local variables


    if (.not.associated(buff)) call increase_buffer(buff)
    if (n>buff%size) call increase_buffer(buff)

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
    buff%data(13,n)=berg%depth
    buff%data(14,n)=berg%width
    buff%data(15,n)=berg%length
    buff%data(16,n)=berg%mass_scaling

  end subroutine pack_berg_into_buffer2

  subroutine increase_buffer(old)
  ! Arguments
  type(buffer), pointer :: old
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size

    if (.not.associated(old)) then
      new_size=3
    else
      new_size=old%size+3
    endif
    allocate(new)
    allocate(new%data(buffer_width,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new

  end subroutine increase_buffer

  subroutine pack_berg_into_buffer(berg, buffer, n)
  ! Arguments
  type(iceberg), pointer :: berg
  real, dimension(buffer_width,max_bergs_in_buffer), intent(inout) :: buffer
  integer, intent(in) :: n
  ! Local variables

    if (n.gt.max_bergs_in_buffer) then
      write(stderr(),*) 'diamond, pack_berg_into_buffer: too many icebergs being packed. Use a bigger buffer!'
      call error_mesg('diamond, pack_berg_into_buffer', 'communication buffer is too small!', FATAL)
    endif

    buffer(1,n)=berg%lon
    buffer(2,n)=berg%lat
    buffer(3,n)=berg%uvel
    buffer(4,n)=berg%vvel
    buffer(5,n)=berg%xi
    buffer(6,n)=berg%yj
    buffer(7,n)=berg%start_lon
    buffer(8,n)=berg%start_lat
    buffer(9,n)=float(berg%start_year)
    buffer(10,n)=berg%start_day
    buffer(11,n)=berg%start_mass
    buffer(12,n)=berg%mass
    buffer(13,n)=berg%depth
    buffer(14,n)=berg%width
    buffer(15,n)=berg%length
    buffer(16,n)=berg%mass_scaling

  end subroutine pack_berg_into_buffer

  subroutine unpack_berg_from_buffer(first, buffer, n)
  ! Arguments
  type(iceberg), pointer :: first
  real, dimension(buffer_width,max_bergs_in_buffer), intent(inout) :: buffer
  integer, intent(in) :: n
  ! Local variables
 !real :: lon, lat, uvel, vvel, xi, yj
 !real :: start_lon, start_lat, start_day, start_mass
 !integer :: ine, jne, start_year
  logical :: lres
  type(iceberg) :: localberg

    localberg%lon=buffer(1,n)
    localberg%lat=buffer(2,n)
    localberg%uvel=buffer(3,n)
    localberg%vvel=buffer(4,n)
    localberg%xi=buffer(5,n)
    localberg%yj=buffer(6,n)
    localberg%start_lon=buffer(7,n)
    localberg%start_lat=buffer(8,n)
    localberg%start_year=int(buffer(9,n)+0.5)
    localberg%start_day=buffer(10,n)
    localberg%start_mass=buffer(11,n)
    localberg%mass=buffer(12,n)
    localberg%depth=buffer(13,n)
    localberg%width=buffer(14,n)
    localberg%length=buffer(15,n)
    localberg%mass_scaling=buffer(16,n)
    lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
    if (lres) then
      call add_new_berg_to_list(first, localberg)
    else
      lres=find_cell_wide(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      if (lres) then
        call add_new_berg_to_list(first, localberg)
      else
        write(stderr(),'("diamond, unpack_berg_from_buffer pe=(",i3,a,2i4,a,2f8.2)') mpp_pe(),') Failed to find i,j=',localberg%ine,localberg%jne,' for lon,lat=',localberg%lon,localberg%lat
        write(stderr(),*) localberg%uvel,localberg%vvel
        write(stderr(),*) grd%isc,grd%iec,grd%jsc,grd%jec
        write(stderr(),*) grd%isd,grd%ied,grd%jsd,grd%jed
        write(stderr(),*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
        write(stderr(),*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
        write(stderr(),*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
        write(stderr(),*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
        debug=.true.
        lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
        write(stderr(),*) lres
        call error_mesg('diamond, unpack_berg_from_buffer', 'can not find a cell to place berg in!', FATAL)
      endif
    endif

  end subroutine unpack_berg_from_buffer

end subroutine send_bergs_to_other_pes

! ##############################################################################

subroutine icebergs_init(bergs, &
             gni, gnj, layout, axes, maskmap, x_cyclic, tripolar_grid, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot)
! Arguments
type(icebergs), pointer :: bergs
integer, intent(in) :: gni, gnj, layout(2), axes(2)
logical, intent(in), optional :: maskmap(:,:)
logical, intent(in) :: x_cyclic, tripolar_grid
real, intent(in) :: dt
type (time_type), intent(in) :: Time ! current time
real, dimension(:,:), intent(in) :: ice_lon, ice_lat, ice_wet
real, dimension(:,:), intent(in) :: ice_dx, ice_dy, ice_area
real, dimension(:,:), intent(in) :: cos_rot, sin_rot
! Namelist parameters (and defaults)
integer :: halo=4 ! Width of halo region
integer :: traj_sample_hrs=4 ! Period between sampling of position for trajectory storage
integer :: verbose_hrs=12 ! Period between verbose messages
real, dimension(nclasses) :: initial_mass=(/8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11/) ! Mass thresholds between iceberg classes (kg)
real, dimension(nclasses) :: distribution=(/0.24, 0.12, 0.15, 0.18, 0.12, 0.07, 0.03, 0.03, 0.03, 0.02/) ! Fraction of calving to apply to this class (non-dim)
real, dimension(nclasses) :: mass_scaling=100. ! Ratio between effective and real iceberg mass (non-dim)
real, dimension(nclasses) :: initial_depth=(/40., 67., 133., 175., 250., 250., 250., 250., 250., 250./) ! Depth of newly calved bergs (m)
namelist /icebergs_nml/ verbose, budget, halo, traj_sample_hrs, initial_mass, distribution, mass_scaling, initial_depth, verbose_hrs
! Local variables
integer :: ierr, iunit, i, j, id_class, axes3d(3), is,ie,js,je
type(icebergs_gridded), pointer :: grd
real :: minl
logical :: lerr

! Read namelist parameters
 !write(stderr(),*) 'diamond: reading namelist'
  iunit = open_namelist_file()
  read  (iunit, icebergs_nml,iostat=ierr)
  ierr = check_nml_error(ierr, 'icebergs_nml')
  call close_file(iunit)

! Log version and parameters
  call write_version_number(version, tagname)
  write (stdlog(), icebergs_nml)

! Allocate overall structure
 !write(stderr(),*) 'diamond: allocating bergs'
  allocate(bergs)
  allocate(bergs%grd)
  grd=>bergs%grd ! For convenience to avoid bergs%grd%X
 !write(stderr(),*) 'diamond: allocating domain'
  allocate(grd%domain)

! Clocks
  bergs%clock=mpp_clock_id( 'Icebergs', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  call mpp_clock_begin(bergs%clock)

! Set up iceberg domain
 !write(stderr(),*) 'diamond: defining domain'
  if(tripolar_grid) then
    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
!                            maskmap=maskmap, &
                             xflags=CYCLIC_GLOBAL_DOMAIN, xhalo=halo,  &
                             yflags=FOLD_NORTH_EDGE, yhalo=halo, name='diamond')
  else if(x_cyclic) then
    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
!                            maskmap=maskmap, &
                             xflags=CYCLIC_GLOBAL_DOMAIN, &
                             xhalo=halo, yhalo=halo, name='diamond')
  else
    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
!                            maskmap=maskmap, &
                             xhalo=halo, yhalo=halo, name='diamond')
  endif

 !write(stderr(),*) 'diamond: get compute domain'
  call mpp_get_compute_domain( grd%domain, grd%isc, grd%iec, grd%jsc, grd%jec )
  call mpp_get_data_domain( grd%domain, grd%isd, grd%ied, grd%jsd, grd%jed )

  call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
  call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
  call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
  call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)
 !write(stderr(),'(a,6i4)') 'diamond, icebergs_init: pe,n,s,e,w =',mpp_pe(),grd%pe_N,grd%pe_S,grd%pe_E,grd%pe_W, NULL_PE

 !if (verbose) &
 !write(stderr(),'(a,i3,a,4i4,a,4f8.2)') 'diamond, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
 !     grd%isc,grd%iec,grd%jsc,grd%jec, &
 !     ' [lon|lat][min|max]=', minval(ice_lon),maxval(ice_lon),minval(ice_lat),maxval(ice_lat)
 !write(stderr(),*) 'diamond, int args = ', mpp_pe(),gni, gnj, layout, axes


 !write(stderr(),*) 'diamond: allocating grid'
  allocate( grd%lon(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=999.
  allocate( grd%lat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=999.
  allocate( grd%dx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dx(:,:)=0.
  allocate( grd%dy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dy(:,:)=0.
  allocate( grd%area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%area(:,:)=0.
  allocate( grd%msk(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%msk(:,:)=0.
  allocate( grd%cos(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cos(:,:)=1.
  allocate( grd%sin(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sin(:,:)=0.
  allocate( grd%calving(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%calving(:,:)=0.
  allocate( grd%melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt(:,:)=0.
  allocate( grd%mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%mass(:,:)=0.
  allocate( grd%stored_ice(grd%isd:grd%ied, grd%jsd:grd%jed, nclasses) ); grd%stored_ice(:,:,:)=0.
  allocate( grd%uo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%uo(:,:)=0.
  allocate( grd%vo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vo(:,:)=0.
  allocate( grd%ui(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ui(:,:)=0.
  allocate( grd%vi(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vi(:,:)=0.
  allocate( grd%ua(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ua(:,:)=0.
  allocate( grd%va(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%va(:,:)=0.
  allocate( grd%ssh(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ssh(:,:)=0.
  allocate( grd%sst(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sst(:,:)=0.
  allocate( grd%cn(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cn(:,:)=0.
  allocate( grd%hi(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%hi(:,:)=0.
  allocate( grd%tmp(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%tmp(:,:)=0.

 !write(stderr(),*) 'diamond: copying grid'
  ! Copy data declared on ice model computational domain
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
  grd%lon(is:ie,js:je)=ice_lon(:,:)
  grd%lat(is:ie,js:je)=ice_lat(:,:)
  grd%area(is:ie,js:je)=ice_area(:,:)*(4.*pi*radius*radius)
  ! Copy data declared on ice model data domain
  is=grd%isc-1; ie=grd%iec+1; js=grd%jsc-1; je=grd%jec+1
  grd%dx(is:ie,js:je)=ice_dx(:,:)
  grd%dy(is:ie,js:je)=ice_dy(:,:)
  grd%msk(is:ie,js:je)=ice_wet(:,:)
  grd%cos(is:ie,js:je)=cos_rot(:,:)
  grd%sin(is:ie,js:je)=sin_rot(:,:)

  call mpp_update_domains(grd%lon, grd%domain)
  call mpp_update_domains(grd%lat, grd%domain)
  call mpp_update_domains(grd%dy, grd%dx, grd%domain)
  call mpp_update_domains(grd%area, grd%domain)
  call mpp_update_domains(grd%msk, grd%domain)
  call mpp_update_domains(grd%cos, grd%domain)
  call mpp_update_domains(grd%sin, grd%domain)

  ! Sanitize lon and lat at the SW edges
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      if (grd%lon(i,j).gt.900.) grd%lon(i,j)=grd%lon(i,j+1)
      if (grd%lat(i,j).gt.900.) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
  enddo; enddo

  do j=grd%jsd,grd%jed; do i=grd%isd,grd%ied
      if (grd%lon(i,j).gt.900.) write(stderr(),*) 'bad lon: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1
      if (grd%lat(i,j).gt.900.) write(stderr(),*) 'bad lat: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1
  enddo; enddo

  ! Sanitize lon for the tile (need continuous longitudes within one tile)
  do j=grd%jsc+1,grd%jed; i=grd%isc
      minl=grd%lon(i,j-1)-180.
      grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo
  do j=grd%jsc-1,grd%jed,-1; i=grd%isc
      minl=grd%lon(i,j+1)-180.
      grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo
  do j=grd%jsd,grd%jed; do i=grd%isc+1,grd%ied
      minl=grd%lon(i-1,j)-180.
      grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo; enddo
  do j=grd%jsd,grd%jed; do i=grd%isc-1,grd%isd,-1
      minl=grd%lon(i+1,j)-180.
      grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo; enddo

 !if (mpp_pe().eq.0) then
 !  write(stderr(),'(a3,16i7)') ' ',(i,i=grd%isd,grd%ied)
 !  do j=grd%jsd,grd%jed
 !    write(stderr(),'(i3,16f7.1)') j,(grd%lon(i,j),i=grd%isd,grd%ied)
 !  enddo
 !endif

 ! Parameters
  bergs%dt=dt
  bergs%traj_sample_hrs=traj_sample_hrs
  bergs%verbose_hrs=verbose_hrs
  bergs%grd%halo=halo
  allocate( bergs%initial_mass(nclasses) ); bergs%initial_mass(:)=initial_mass(:)
  allocate( bergs%distribution(nclasses) ); bergs%distribution(:)=distribution(:)
  allocate( bergs%mass_scaling(nclasses) ); bergs%mass_scaling(:)=mass_scaling(:)
  allocate( bergs%initial_depth(nclasses) ); bergs%initial_depth(:)=initial_depth(:)
  allocate( bergs%initial_width(nclasses) ); bergs%initial_width(:)=sqrt((2./3.)*initial_mass(:)/(rho_ice*initial_depth(:)))
  allocate( bergs%initial_length(nclasses) ); bergs%initial_length(:)=1.5*bergs%initial_width(:)

  call read_restart_bergs(bergs)
  call read_restart_calving(bergs)

 !if (verbose) call print_bergs(stderr(),bergs,'icebergs_init, initial status')

  ! Diagnostics
  id_class = diag_axis_init('mass', initial_mass, 'kg','Z', 'iceberg mass')
  axes3d(1:2)=axes
  axes3d(3)=id_class
  grd%id_calving=register_diag_field('icebergs', 'calving', axes, Time, &
     'Incoming Calving mass rate', 'kg/s')
  grd%id_accum=register_diag_field('icebergs', 'accum_calving', axes, Time, &
     'Accumulated calving mass rate', 'kg/s')
  grd%id_unused=register_diag_field('icebergs', 'unused_calving', axes, Time, &
     'Unused calving mass rate', 'kg/s')
  grd%id_melt=register_diag_field('icebergs', 'melt', axes, Time, &
     'Iceberg melt mass rate', 'kg/(m^2*s)')
  grd%id_mass=register_diag_field('icebergs', 'mass', axes, Time, &
     'Iceberg density field', 'kg/(m^2)')
  grd%id_stored_ice=register_diag_field('icebergs', 'stored_ice', axes3d, Time, &
     'Accumulated ice mass by class', 'kg')
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

  ! Static fields
  id_class=register_static_field('icebergs', 'lon', axes, &
               'longitude (corners)', 'degrees_E',require=.false.)
  if (id_class>0) lerr=send_data(id_class, grd%lon(grd%isc:grd%iec,grd%jsc:grd%jec), Time);
  id_class=register_static_field('icebergs', 'lat', axes, &
               'latitude (corners)', 'degrees_N',require=.false.)
  if (id_class>0) lerr=send_data(id_class, grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec), Time);
  id_class=register_static_field('icebergs', 'area', axes, &
               'cell area', 'm^2',require=.false.)
  if (id_class>0) lerr=send_data(id_class, grd%area(grd%isc:grd%iec,grd%jsc:grd%jec), Time);
  id_class=register_static_field('icebergs', 'mask', axes, &
               'wet point mask', 'none',require=.false.)
  if (id_class>0) lerr=send_data(id_class, grd%msk(grd%isc:grd%iec,grd%jsc:grd%jec), Time);

 !write(stderr(),*) 'diamond: done'
  call mpp_clock_end(bergs%clock)
end subroutine icebergs_init

! ##############################################################################

subroutine read_restart_bergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
integer :: k, ierr, ncid, dimid, nbergs_in_file
integer :: lonid, latid, uvelid, vvelid
integer :: massid, depthid, widthid, lengthid
integer :: start_lonid, start_latid, start_yearid, start_dayid, start_massid
integer :: scaling_id
logical :: lres, found_restart
real :: lon0, lon1, lat0, lat1
character(len=30) :: filename
type(icebergs_gridded), pointer :: grd
type(iceberg) :: localberg ! NOT a pointer but an actual local variable

  ! For convenience
  grd=>bergs%grd

  ! Find a restart file
  do
    filename='INPUT/icebergs.res.nc'; inquire(file=filename,exist=found_restart)
    if (found_restart) exit
    write(filename(1:27),'("INPUT/icebergs.res.nc.",I4.4)') mpp_pe()
    inquire(file=filename,exist=found_restart)
    if (found_restart) exit
    write(filename(1:21),'("icebergs.res.nc.",I4.4)') mpp_pe()
    inquire(file=filename,exist=found_restart)
    if (found_restart) exit
    filename='icebergs.res.nc'; inquire(file=filename,exist=found_restart)
    if (found_restart) exit
    if (verbose.and.mpp_pe()==mpp_root_pe()) write(stderr(),*) 'diamond, read_restart_bergs: no restart file found'
    return
  enddo 

  if (verbose.and.mpp_pe()==mpp_root_pe()) write(stderr(),*) 'diamond, read_restart_bergs: found restart file = ',filename

  ierr=nf_open(filename, NF_NOWRITE, ncid)
  if (ierr .ne. NF_NOERR) write(stderr(),*) 'diamond, read_restart_bergs: nf_open failed'

  ierr=nf_inq_unlimdim(ncid, dimid)
  if (ierr .ne. NF_NOERR) write(stderr(),*) 'diamond, read_restart_bergs: nf_inq_unlimdim failed'

  ierr=nf_inq_dimlen(ncid, dimid, nbergs_in_file)
  if (ierr .ne. NF_NOERR) write(stderr(),*) 'diamond, read_restart_bergs: nf_inq_dimlen failed'
 !write(stderr(),*) 'diamond, read_restart_bergs: nbergs in file =', nbergs_in_file

  lonid=inq_var(ncid, 'lon')
  latid=inq_var(ncid, 'lat')
  uvelid=inq_var(ncid, 'uvel')
  vvelid=inq_var(ncid, 'vvel')
  massid=inq_var(ncid, 'mass')
  depthid=inq_var(ncid, 'depth')
  widthid=inq_var(ncid, 'width')
  lengthid=inq_var(ncid, 'length')
  start_lonid=inq_var(ncid, 'start_lon')
  start_latid=inq_var(ncid, 'start_lat')
  start_yearid=inq_var(ncid, 'start_year')
  start_dayid=inq_var(ncid, 'start_day')
  start_massid=inq_var(ncid, 'start_mass')
  scaling_id=inq_var(ncid, 'mass_scaling')

  ! Find approx outer bounds for tile
  lon0=minval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lon1=maxval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat0=minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat1=maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
 !write(stderr(),'(a,i3,a,4f10.3)') 'diamond, read_restart_bergs: (',mpp_pe(),') ',lon0,lon1,lat0,lat1
  do k=1, nbergs_in_file
   !write(stderr(),*) 'diamond, read_restart_bergs: reading berg ',k
    localberg%lon=get_double(ncid, lonid, k)
    localberg%lat=get_double(ncid, latid, k)
    ! Test if this berg is within the maximum possible bounds of tile
    if ( sum_sign_dot_prod4(lon0,lat0,lon1,lat0,lon1,lat1,lon0,lat1,localberg%lon,localberg%lat) ) then
     !write(stderr(),*) 'diamond, read_restart_bergs: berg ',k,' is at ',lon,lat
      lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      if (lres) then
        localberg%uvel=get_double(ncid, uvelid, k)
        localberg%vvel=get_double(ncid, vvelid, k)
        localberg%mass=get_double(ncid, massid, k)
        localberg%depth=get_double(ncid, depthid, k)
        localberg%width=get_double(ncid, widthid, k)
        localberg%length=get_double(ncid, lengthid, k)
        localberg%start_lon=get_double(ncid, start_lonid, k)
        localberg%start_lat=get_double(ncid, start_latid, k)
        localberg%start_year=get_int(ncid, start_yearid, k)
        localberg%start_day=get_double(ncid, start_dayid, k)
        localberg%start_mass=get_double(ncid, start_massid, k)
        localberg%mass_scaling=get_double(ncid, scaling_id, k)
       !write(stderr(),'(a,i3,a,2i4,2f8.2)') 'diamond, read_restart_bergs: found cell (',mpp_pe(),') ',i,j,grd%lon(i,j),grd%lat(i,j)
        lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
        call add_new_berg_to_list(bergs%first, localberg)
       !if (verbose) call print_berg(stderr(), bergs%first, 'read_restart_bergs, add_new_berg_to_list')
      endif
    endif
  enddo

end subroutine read_restart_bergs

! ##############################################################################

subroutine read_restart_calving(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
integer :: k
character(len=30) :: filename
type(icebergs_gridded), pointer :: grd

  ! For convenience
  grd=>bergs%grd

  ! Read stored ice
  filename='INPUT/calving.res.nc'
  if (file_exist(filename)) then
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(stderr(),*) 'diamond, read_restart_calving: reading ',filename
    call read_data(filename, 'stored_ice', bergs%grd%stored_ice, bergs%grd%domain)
    bergs%restarted=.true.
  else
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(stderr(),*) 'diamond, read_restart_calving: initialization stored ice to random numbers '
    call random_number(bergs%grd%stored_ice)
    do k=1, nclasses
      bergs%grd%stored_ice(:,:,k)=bergs%grd%stored_ice(:,:,k)*bergs%grd%msk(:,:) &
                                 *bergs%initial_mass(k)*bergs%mass_scaling(k)
    enddo
  endif

end subroutine read_restart_calving

! ##############################################################################

subroutine add_new_berg_to_list(first, bergvals)
! Arguments
type(iceberg), pointer :: first
type(iceberg), intent(in) :: bergvals
! Local variables
type(iceberg), pointer :: new=>NULL()

  new=>NULL()
  call create_iceberg(new, bergvals)

  if (associated(first)) then
   !write(stderr(),*) 'diamond, add_new_berg_to_list: moving old pointers',mpp_pe()
    new%next=>first
    first%prev=>new
    first=>new
  else
   !write(stderr(),*) 'diamond, add_new_berg_to_list: setting new pointer',mpp_pe()
    first=>new
  endif
 !call print_berg(stderr(), first, 'add_new_berg_to_list')

  !Clear new
  new=>NULL()

end subroutine add_new_berg_to_list

! ##############################################################################

subroutine create_iceberg(berg, bergvals)
! Arguments
type(iceberg), pointer :: berg
type(iceberg), intent(in) :: bergvals
! Local variables

  if (associated(berg)) then
    write(stderr(),*) 'diamond, create_iceberg: berg already associated!!!!',mpp_pe()
    call error_mesg('diamond, create_iceberg', 'berg already associated. This should not happen!', FATAL)
  endif
  allocate(berg)
  berg=bergvals
  berg%prev=>NULL()
  berg%next=>NULL()

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

 !call print_berg(stderr(), berg, 'destroy_iceberg')

  ! Bye-bye berg
  deallocate(berg)

end subroutine destroy_iceberg

! ##############################################################################

subroutine print_berg(iochan, berg, label)
! Arguments
integer, intent(in) :: iochan
type(iceberg), pointer :: berg
character(len=*) :: label
! Local variables

  write(iochan,'("diamond, print_berg: ",a," pe=(",i3,a,2i5,3(a,2f10.4),a,2l2)') &
    label, mpp_pe(), ') i,j=',berg%ine, berg%jne, &
    ' xi,yj=', berg%xi, berg%yj, &
    ' lon,lat=', berg%lon, berg%lat, &
    ' u,v=', berg%uvel, berg%vvel, &
    ' p,n=', associated(berg%prev), associated(berg%next)
  write(iochan,'("diamond, print_berg: ",a," pe=(",i3,") start lon,lat,yr,day,mass=",2f10.4,i5,f7.2,1pe12.4)') &
    label, mpp_pe(), berg%start_lon, berg%start_lat, &
    berg%start_year, berg%start_day, berg%start_mass

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

  this=>bergs%first
  do while(associated(this))
    call print_berg(iochan, this, label)
    this=>this%next
  enddo
  nbergs=count_bergs(bergs)
  nnbergs=nbergs
  call mpp_sum(nnbergs)
  if (nbergs.gt.0) write(iochan,'("diamond, ",a," there are",i5," bergs out of",i6," on PE ",i4)') label, nbergs, nnbergs, mpp_pe()

end subroutine print_bergs

! ##############################################################################

integer function count_bergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: this

  count_bergs=0
  this=>bergs%first
  do while(associated(this))
    count_bergs=count_bergs+1
    this=>this%next
  enddo

end function count_bergs

! ##############################################################################

subroutine record_posn(berg, current_year, current_yearday)
! Arguments
type(iceberg), pointer :: berg
integer, intent(in) :: current_year
real, intent(in) :: current_yearday
! Local variables
type(xyt) :: posn

  posn%lon=berg%lon
  posn%lat=berg%lat
  posn%year=current_year
  posn%day=current_yearday
  posn%uvel=berg%uvel
  posn%vvel=berg%vvel
  posn%mass=berg%mass
  posn%depth=berg%depth
  posn%width=berg%width
  posn%length=berg%length
  posn%uo=berg%uo
  posn%vo=berg%vo
  posn%ui=berg%ui
  posn%vi=berg%vi
  posn%ua=berg%ua
  posn%va=berg%va
  posn%ssh_x=berg%ssh_x
  posn%ssh_y=berg%ssh_y
  posn%sst=berg%sst

  call push_posn(berg%trajectory, posn)

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

! ##############################################################################

subroutine move_trajectory(bergs, berg)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
! Local variables
type(xyt), pointer :: next, last
type(xyt) :: vals

  ! If the trajectory is empty, ignore it
  if (.not.associated(berg%trajectory)) return

  ! Push identifying info into first posn (note reverse order due to stack)
  vals%lon=berg%start_lon
  vals%lat=berg%start_lat
  vals%year=berg%start_year
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
  berg%trajectory=>NULL()

end subroutine move_trajectory

! ##############################################################################

logical function find_cell(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell=.false.; oi=-999; oj=-999

  ! This should be a bisection search but ...
  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell=.true.
       !write(stderr(),'(a,i3,a,10f7.1,2i4)') 'diamond, find_cell: (',mpp_pe(),') ', &
       !                   grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
       !                   grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
       !                   grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
       !                   grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
       !                   x, y, i, j
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

  ! This should be a bisection search but ...
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell_wide=.true.
        return
      endif
  enddo; enddo

end function find_cell_wide

! ##############################################################################

logical function is_point_in_cell(grd, x, y, i, j)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
! Local variables
real :: xlo, xhi, ylo, yhi, xx

  ! Safety check index bounds
  if (i-1.lt.grd%isd.or.i.gt.grd%ied.or.j-1.lt.grd%jsd.or.j.gt.grd%jed) then
    write(stderr(),'(a,i3,(a,3i4))') &
                     'diamond, is_point_in_cell: pe=(',mpp_pe(),') i,s,e=', &
                     i,grd%isd,grd%ied,' j,s,e=', j,grd%jsd,grd%jed
    call error_mesg('diamond, is_point_in_cell', 'berg is off the PE!', FATAL)
  endif

  is_point_in_cell=.false.

  ! Test crude bounds
  xlo=min( grd%lon(i-1,j-1), grd%lon(i,j-1), grd%lon(i-1,j), grd%lon(i,j) )
  xhi=max( grd%lon(i-1,j-1), grd%lon(i,j-1), grd%lon(i-1,j), grd%lon(i,j) )
  xx=modulo(x-(xlo-180.),360.)+(xlo-180.)
  if (xx.lt.xlo .or. xx.gt.xhi) return
  ylo=min( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  yhi=max( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  if (y.lt.ylo .or. y.gt.yhi) return
  
  is_point_in_cell=sum_sign_dot_prod4(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y) 
  if (debug) then
    write(stderr(),*) 'diamond, is_point_in_cell: inside crude bounds but no in cell',mpp_pe()
    write(stderr(),*) 'diamond, is_point_in_cell: ',xlo,xhi,x,xx
    write(stderr(),*) 'diamond, is_point_in_cell: ',ylo,yhi,y
  endif

end function is_point_in_cell

! ##############################################################################

logical function sum_sign_dot_prod4(x0, y0, x1, y1, x2, y2, x3, y3, x, y)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x, y
! Local variables
real :: p0,p1,p2,p3,xx
real :: l0,l1,l2,l3

  sum_sign_dot_prod4=.false.
  xx=modulo(x-(x0-180.),360.)+(x0-180.)

  l0=(xx-x0)*(y1-y0)-(y-y0)*(x1-x0)
  l1=(xx-x1)*(y2-y1)-(y-y1)*(x2-x1)
  l2=(xx-x2)*(y3-y2)-(y-y2)*(x3-x2)
  l3=(xx-x3)*(y0-y3)-(y-y3)*(x0-x3)

  p0=sign(1., l0);!if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3);!if (l3.eq.0.) p3=0.

 !if ( abs( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) - abs((p0+p2)+(p1+p3)) ).lt.0.5 ) then
  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    sum_sign_dot_prod4=.true.
  endif

 !if (((x.ge.min(x0,x1,x2,x3)).and.(x.le.max(x0,x1,x2,x3)) &
 !.and.(y.ge.min(y0,y1,y2,y3)).and.(y.le.max(y0,y1,y2,y3))) .or. debug) &
 !then
  if (debug) then
   write(stderr(),'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: x=',mpp_pe(),':', &
                           x0,x1,x2,x3, x
   write(stderr(),'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: y=',mpp_pe(),':', &
                           y0,y1,y2,y3, y
   write(stderr(),'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: l=',mpp_pe(),':', &
                           l0,l1,l2,l3
   write(stderr(),'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: p=',mpp_pe(),':', &
                           p0,p1,p2,p3, abs( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) - abs((p0+p2)+(p1+p3)) )
  endif

end function sum_sign_dot_prod4

! ##############################################################################

logical function pos_within_cell(grd, x, y, i, j, xi, yj)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
real, intent(out) :: xi, yj
! Local variables
real :: x1,y1,x2,y2,x3,y3,x4,y4,dx,dy
real :: alpha, beta, gamma, delta, epsilon, kappa, a, b, c, d

  pos_within_cell=.true.; xi=-999.; yj=-999.

  x1=grd%lon(i-1,j-1)
  y1=grd%lat(i-1,j-1)
  x2=grd%lon(i  ,j-1)
  y2=grd%lat(i  ,j-1)
  x3=grd%lon(i  ,j  )
  y3=grd%lat(i  ,j  )
  x4=grd%lon(i-1,j  )
  y4=grd%lat(i-1,j  )
 !write(stderr(),'(a,i3,4f8.2)') 'pos_within_cell: x1..x4 ',mpp_pe(),x1,x2,x3,x4
 !write(stderr(),'(a,i3,3f8.2)') 'pos_within_cell: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
 !write(stderr(),'(a,i3,4f8.2)') 'pos_within_cell: y1..y4 ',mpp_pe(),y1,y2,y3,y4
 !write(stderr(),'(a,i3,3f8.2)') 'pos_within_cell: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1

 !if (.not. is_point_in_cell(grd, x, y, i, j) ) then
 !   write(stderr(),'(a,i3,a,8f8.2,a)') 'diamond, pos_within_cell: (',mpp_pe(),') ', &
 !                   x1, y1, x2, y2, x3, y3, x4, y4, ' NOT IN CELL!'
 !endif

  alpha=x2-x1
  delta=y2-y1
  beta=x4-x1
  epsilon=y4-y1
  gamma=(x3-x1)-(alpha+beta)
  kappa=(y3-y1)-(delta+epsilon)
 !write(stderr(),'(a,i3,1p6e12.4)') 'pos_within_cell: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa

  a=(kappa*beta-gamma*epsilon)
  dx=modulo(x-(x1-180.),360.)+(x1-180.)-x1
  dy=y-y1
  b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
  c=(alpha*dy-delta*dx)
  if (a.ne.0.) then
    d=0.25*(b**2)-a*c
   !write(stderr(),'(a,i3,1p6e12.4)') 'pos_within_cell: coeffs a,b,c,d,dx,dy',mpp_pe(),a,b,c,dx,dy
    if (d.ge.0.) then
      y1=-0.5*b/a-sqrt(d)/a
      y2=-0.5*b/a+sqrt(d)/a
      if (abs(y1-0.5).lt.abs(y2-0.5)) then; yj=y1; else; yj=y2; endif
     !write(stderr(),'(a,i3,1p3e12.4)') 'Roots for y = ',mpp_pe(),y1,y2,yj
    else
      call dump_pos
      write(stderr(),'(a,i3)') 'pos_within_cell: b<0 in quadratic root solver!!!!',mpp_pe()
     !write(stderr(),'(a,i3,1p6e12.4)') 'pos_within_cell: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderr(),'(a,i3,1p6e12.4)') 'pos_within_cell: coeffs a,b,c,d,dx,dy',mpp_pe(),a,b,c,d,dx,dy
      call error_mesg('diamond, pos_within_cell', 'We have complex roots. The grid must be very distorted!', FATAL)
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
      call dump_pos
     !write(stderr(),'(a,i3,1p6e12.4)') 'pos_within_cell: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderr(),'(a,i3,1p2e12.4)') 'pos_within_cell: coeffs a,b',mpp_pe(),a,b
      call error_mesg('diamond, pos_within_cell', 'Can not invert either linear equaton for xi! This should not happen!', FATAL)
    endif
  endif
 !write(stderr(),'(a,i3,2e12.4)') 'pos_within_cell: xi,yj=',mpp_pe(),xi,yj

  if (abs(xi-0.5).gt.0.5) pos_within_cell=.false.
  if (abs(yj-0.5).gt.0.5) pos_within_cell=.false.

  contains

  subroutine dump_pos
    write(stderr(),'(a,i3,4f8.2)') 'pos_within_cell: x1..x4 ',mpp_pe(),x1,x2,x3,x4
    write(stderr(),'(a,i3,3f8.2)') 'pos_within_cell: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
    write(stderr(),'(a,i3,4f8.2)') 'pos_within_cell: y1..y4 ',mpp_pe(),y1,y2,y3,y4
    write(stderr(),'(a,i3,3f8.2)') 'pos_within_cell: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
    write(stderr(),'(a,i3,1p6e12.4)') 'pos_within_cell: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
  end subroutine dump_pos

end function pos_within_cell

! ##############################################################################

real function sum_icebergs_mass(first)
! Arguments
type(iceberg), pointer :: first
! Local variables
type(iceberg), pointer :: this

  sum_icebergs_mass=0.
  this=>first
  do while(associated(this))
    sum_icebergs_mass=sum_icebergs_mass+this%mass*this%mass_scaling
    this=>this%next
  enddo

end function sum_icebergs_mass

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
    berg_mass=sum_icebergs_mass(bergs%first)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=stored_mass+berg_mass

  case (ISTOCK_HEAT)
    berg_mass=sum_icebergs_mass(bergs%first)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=(stored_mass+berg_mass)*HLF ! HLF is in (J/kg) from constants_mod

  case default
    value = 0.0

end select

end subroutine icebergs_stock_pe

! ##############################################################################

subroutine icebergs_end(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: this, next

  if (.not.associated(bergs)) return

  call write_restart(bergs)

  ! Delete bergs and structures
  this=>bergs%first
  do while (associated(this))
    next=>this%next
    call move_trajectory(bergs, this)
    call destroy_iceberg(this)
    this=>next
  enddo

  if (associated(bergs%trajectories)) then
    call write_trajectory(bergs%trajectories)
  endif

  deallocate(bergs%grd%lon)
  deallocate(bergs%grd%lat)
  deallocate(bergs%grd%dx)
  deallocate(bergs%grd%dy)
  deallocate(bergs%grd%area)
  deallocate(bergs%grd%msk)
  deallocate(bergs%grd%cos)
  deallocate(bergs%grd%sin)
  deallocate(bergs%grd%calving)
  deallocate(bergs%grd%melt)
  deallocate(bergs%grd%mass)
  deallocate(bergs%grd%tmp)
  deallocate(bergs%grd%stored_ice)
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
  deallocate(bergs%initial_depth)
  deallocate(bergs%initial_width)
  deallocate(bergs%initial_length)
  deallocate(bergs)

 !if (verbose) write(stderr(),*) 'diamond: icebergs_end complete',mpp_pe()

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

subroutine write_restart(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, uvelid, vvelid
integer :: massid, depthid, lengthid, widthid
integer :: start_lonid, start_latid, start_yearid, start_dayid, start_massid
integer :: scaling_id
character(len=30) :: filename
type(iceberg), pointer :: this

  ! Only create a restart file for this PE if we have anything to say
  if (associated(bergs%first)) then

    write(filename(1:28),'("RESTART/icebergs.res.nc.",I4.4)') mpp_pe()
    if (verbose) write(stderr(),*) 'diamond, write_restart: creating ',filename

    iret = nf_create(filename, NF_CLOBBER, ncid)
    if (iret .ne. NF_NOERR) write(stderr(),*) 'diamond, write_restart: nf_create failed'

    ! Dimensions
    iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
    if (iret .ne. NF_NOERR) write(stderr(),*) 'diamond, write_restart: nf_def_dim i failed'

    ! Variables
    lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
    latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
    uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
    vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
    massid = def_var(ncid, 'mass', NF_DOUBLE, i_dim)
    depthid = def_var(ncid, 'depth', NF_DOUBLE, i_dim)
    widthid = def_var(ncid, 'width', NF_DOUBLE, i_dim)
    lengthid = def_var(ncid, 'length', NF_DOUBLE, i_dim)
    start_lonid = def_var(ncid, 'start_lon', NF_DOUBLE, i_dim)
    start_latid = def_var(ncid, 'start_lat', NF_DOUBLE, i_dim)
    start_yearid = def_var(ncid, 'start_year', NF_INT, i_dim)
    start_dayid = def_var(ncid, 'start_day', NF_DOUBLE, i_dim)
    start_massid = def_var(ncid, 'start_mass', NF_DOUBLE, i_dim)
    scaling_id = def_var(ncid, 'mass_scaling', NF_DOUBLE, i_dim)

    ! Attributes
    call put_att(ncid, lonid, 'long_name', 'longitude')
    call put_att(ncid, lonid, 'units', 'degrees_E')
    call put_att(ncid, latid, 'long_name', 'latitude')
    call put_att(ncid, latid, 'units', 'degrees_N')
    call put_att(ncid, uvelid, 'long_name', 'zonal velocity')
    call put_att(ncid, uvelid, 'units', 'm/s')
    call put_att(ncid, vvelid, 'long_name', 'meridional velocity')
    call put_att(ncid, vvelid, 'units', 'm/s')
    call put_att(ncid, massid, 'long_name', 'mass')
    call put_att(ncid, massid, 'units', 'kg')
    call put_att(ncid, depthid, 'long_name', 'depth')
    call put_att(ncid, depthid, 'units', 'm')
    call put_att(ncid, widthid, 'long_name', 'width')
    call put_att(ncid, widthid, 'units', 'm')
    call put_att(ncid, lengthid, 'long_name', 'length')
    call put_att(ncid, lengthid, 'units', 'm')
    call put_att(ncid, start_lonid, 'long_name', 'longitude of calving location')
    call put_att(ncid, start_lonid, 'units', 'degrees_E')
    call put_att(ncid, start_latid, 'long_name', 'latitude of calving location')
    call put_att(ncid, start_latid, 'units', 'degrees_N')
    call put_att(ncid, start_yearid, 'long_name', 'calendar year of calving event')
    call put_att(ncid, start_yearid, 'units', 'years')
    call put_att(ncid, start_dayid, 'long_name', 'year day of calving event')
    call put_att(ncid, start_dayid, 'units', 'days')
    call put_att(ncid, start_massid, 'long_name', 'initial mass of calving berg')
    call put_att(ncid, start_massid, 'units', 'kg')
    call put_att(ncid, scaling_id, 'long_name', 'scaling factor for mass of calving berg')
    call put_att(ncid, scaling_id, 'units', 'none')
    iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_major_version', NF_INT, 1, file_format_major_version)
    iret = nf_put_att_int(ncid, NCGLOBAL, 'file_format_minor_version', NF_INT, 1, file_format_minor_version)

    ! End define mode
    iret = nf_enddef(ncid)
         
    ! Write variables
    this=>bergs%first; i=0
    do while (associated(this))
      i=i+1
      call put_double(ncid, lonid, i, this%lon)
      call put_double(ncid, latid, i, this%lat)
      call put_double(ncid, uvelid, i, this%uvel)
      call put_double(ncid, vvelid, i, this%vvel)
      call put_double(ncid, massid, i, this%mass)
      call put_double(ncid, depthid, i, this%depth)
      call put_double(ncid, widthid, i, this%width)
      call put_double(ncid, lengthid, i, this%length)
      call put_double(ncid, start_lonid, i, this%start_lon)
      call put_double(ncid, start_latid, i, this%start_lat)
      call put_int(ncid, start_yearid, i, this%start_year)
      call put_double(ncid, start_dayid, i, this%start_day)
      call put_double(ncid, start_massid, i, this%start_mass)
      call put_double(ncid, scaling_id, i, this%mass_scaling)
      this=>this%next
    enddo
         
    ! Finish up
    iret = nf_close(ncid)
    if (iret .ne. NF_NOERR) write(stderr(),*) 'diamond, write_restart: nf_close failed'

  endif ! associated(bergs%first)

  ! Write stored ice
  filename='RESTART/calving.res.nc'
  if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(stderr(),*) 'diamond, write_restart: writing ',filename
  call write_data(filename, 'stored_ice', bergs%grd%stored_ice, bergs%grd%domain)

end subroutine write_restart

! ##############################################################################

subroutine write_trajectory(trajectory)
! Arguments
type(xyt), pointer :: trajectory
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, uvelid, vvelid
integer :: uoid, void, uiid, viid, uaid, vaid, sshxid, sshyid, sstid
integer :: mid, did, wid, lid
character(len=30) :: filename
type(xyt), pointer :: this, next

  write(filename(1:30),'("iceberg_trajectories.nc.",I4.4)') mpp_pe()
  if (verbose) write(stderr(),*) 'diamond, write_tracjectory: creating ',filename

  iret = nf_create(filename, NF_CLOBBER, ncid)
  if (iret .ne. NF_NOERR) write(stderr(),*) 'diamond, write_restart: nf_create failed'

  ! Dimensions
  iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
  if (iret .ne. NF_NOERR) write(stderr(),*) 'diamond, write_restart: nf_def_dim i failed'

  ! Variables
  lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
  latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
  yearid = def_var(ncid, 'year', NF_INT, i_dim)
  dayid = def_var(ncid, 'day', NF_DOUBLE, i_dim)
  uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
  vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
  uoid = def_var(ncid, 'uo', NF_DOUBLE, i_dim)
  void = def_var(ncid, 'vo', NF_DOUBLE, i_dim)
  uiid = def_var(ncid, 'ui', NF_DOUBLE, i_dim)
  viid = def_var(ncid, 'vi', NF_DOUBLE, i_dim)
  uaid = def_var(ncid, 'ua', NF_DOUBLE, i_dim)
  vaid = def_var(ncid, 'va', NF_DOUBLE, i_dim)
  mid = def_var(ncid, 'mass', NF_DOUBLE, i_dim)
  did = def_var(ncid, 'depth', NF_DOUBLE, i_dim)
  wid = def_var(ncid, 'width', NF_DOUBLE, i_dim)
  lid = def_var(ncid, 'length', NF_DOUBLE, i_dim)
  sshxid = def_var(ncid, 'ssh_x', NF_DOUBLE, i_dim)
  sshyid = def_var(ncid, 'ssh_y', NF_DOUBLE, i_dim)
  sstid = def_var(ncid, 'sst', NF_DOUBLE, i_dim)

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
  call put_att(ncid, did, 'long_name', 'depth')
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

  ! End define mode
  iret = nf_enddef(ncid)
       
  ! Write variables
  this=>trajectory; i=0
  do while (associated(this))
    i=i+1
    call put_double(ncid, lonid, i, this%lon)
    call put_double(ncid, latid, i, this%lat)
    call put_int(ncid, yearid, i, this%year)
    call put_double(ncid, dayid, i, this%day)
    call put_double(ncid, uvelid, i, this%uvel)
    call put_double(ncid, vvelid, i, this%vvel)
    call put_double(ncid, uoid, i, this%uo)
    call put_double(ncid, void, i, this%vo)
    call put_double(ncid, uiid, i, this%ui)
    call put_double(ncid, viid, i, this%vi)
    call put_double(ncid, uaid, i, this%ua)
    call put_double(ncid, vaid, i, this%va)
    call put_double(ncid, mid, i, this%mass)
    call put_double(ncid, did, i, this%depth)
    call put_double(ncid, wid, i, this%width)
    call put_double(ncid, lid, i, this%length)
    call put_double(ncid, sshxid, i, this%ssh_x)
    call put_double(ncid, sshyid, i, this%ssh_y)
    call put_double(ncid, sstid, i, this%sst)
    next=>this%next
    deallocate(this)
    this=>next
  enddo
  trajectory=>NULL()
       
  ! Finish up
  iret = nf_close(ncid)
  if (iret .ne. NF_NOERR) write(stderr(),*) 'diamond, write_trajectory: nf_close failed'

end subroutine write_trajectory


! ##############################################################################

integer function inq_var(ncid, var)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
! Local variables
integer :: iret

  iret=nf_inq_varid(ncid, var, inq_var)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, inq_var: nf_inq_varid ',var,' failed'
    call error_mesg('diamond, inq_var', 'netcdf function returned a failure!', FATAL)
  endif

end function inq_var

! ##############################################################################

integer function def_var(ncid, var, ntype, idim)
! Arguments
integer, intent(in) :: ncid, ntype, idim
character(len=*), intent(in) :: var
! Local variables
integer :: iret

  iret = nf_def_var(ncid, var, ntype, 1, idim, def_var)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, def_var: nf_def_var failed for ',trim(var)
    call error_mesg('diamond, def_var', 'netcdf function returned a failure!', FATAL)
  endif

end function def_var

! ##############################################################################

subroutine put_att(ncid, id, att, attval)
! Arguments
integer, intent(in) :: ncid, id
character(len=*), intent(in) :: att, attval
! Local variables
integer :: vallen, iret

  vallen=LEN_TRIM(attval)
  iret = nf_put_att_text(ncid, id, att, vallen, attval)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, put_att: nf_put_att_text failed adding', &
      trim(att),' = ',trim(attval)
    call error_mesg('diamond, put_att', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_att

! ##############################################################################

real function get_double(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret

  iret=nf_get_var1_double(ncid, id, i, get_double)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, get_double: nf_get_var1_double failed reading'
    call error_mesg('diamond, get_double', 'netcdf function returned a failure!', FATAL)
  endif

end function get_double

! ##############################################################################

integer function get_int(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret

  iret=nf_get_var1_int(ncid, id, i, get_int)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, get_int: nf_get_var1_int failed reading'
    call error_mesg('diamond, get_int', 'netcdf function returned a failure!', FATAL)
  endif

end function get_int

! ##############################################################################

subroutine put_double(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i
real, intent(in) :: val
! Local variables
integer :: iret

  iret = nf_put_vara_double(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, put_double: nf_put_vara_double failed writing'
    call error_mesg('diamond, put_double', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_double

! ##############################################################################

subroutine put_int(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i, val
! Local variables
integer :: iret

  iret = nf_put_vara_int(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderr(),*) 'diamond, put_int: nf_put_vara_int failed writing'
    call error_mesg('diamond, put_int', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_int

! ##############################################################################

end module
