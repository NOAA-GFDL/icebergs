module driver_data

use fms_mod, only : read_data
use fms_mod, only : file_exist
use fms_mod, only : field_exist
use fms_mod, only : error_mesg
use fms_mod, only : FATAL
use mpp_domains_mod, only : domain2D
use mpp_mod, only : mpp_pe
use mpp_mod, only : mpp_root_pe
use mpp_mod, only: mpp_sync
use constants_mod, only: pi
use netcdf
!include 'netcdf.inc'

implicit none

contains

  subroutine a68_dims(data_dir,ni,nj)
    ! Arguments
    character(len=*) :: data_dir
    integer, intent(out):: ni,nj
    ! Local variables
    character(len=1000):: infile
    integer(kind=4) :: ncid
    integer :: status

    infile=trim(trim(data_dir)//'a68_experiment_ll_p125_grid.nc')

    !Open netCDF file
    call handle_err(nf90_open(infile, nf90_nowrite, ncid))
    !Inquire about the dimensions
    call handle_err(nf90_inquire_dimension(ncid,1,len=nj))
    call handle_err(nf90_inquire_dimension(ncid,2,len=ni))
    !Close netCDF file
    call handle_err(nf90_close(ncid))
  end subroutine a68_dims

  subroutine handle_err(status)
    integer, intent (in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine handle_err


  subroutine a68_prep(data_dir,mpp_domain,lon,lat,dx,dy,area,depth,uo,vo,tauxa,tauya,ssh,&
    vert_int_ocean_vel,tau_is_velocity,Old_wrong_Rearth_and_SSH,REarth)
    !mom_Rearth)
  ! Arguments
  character(len=*) :: data_dir
  type(domain2D) :: mpp_domain
  real :: lon(:,:) !< Longitude or position in x
  real :: lat(:,:) !< Latitude or position in y
  real :: dx(:,:) !< Length of northern edge of cell (m)
  real :: dy(:,:) !< Length of eatern edge of cell (m)
  real :: area(:,:) !< Area of cell (m2)
  real :: depth(:,:) !< Depth of ocean (m)
  real :: uo(:,:) !< Zonal ocean velocities (m/s)
  real :: vo(:,:) !< Meridional ocean velocities (m/s)
  real :: tauxa(:,:) !< Zonal wind stress (Pa)
  real :: tauya(:,:) !< Meridional wind stress (Pa)
  real :: ssh(:,:) !< Effective sea surface height (m)
  logical :: Old_wrong_Rearth_and_SSH !< use Rearth=6.36e6 and SSH on B grid
  logical :: mom_Rearth !< Use Rearth=6.378e6
  logical :: vert_int_ocean_vel !< If true, use ocean velocity vertically-integrated over berg draft
  logical :: tau_is_velocity !< If true, read in JRA-55 wind velocity (m/s), not ORAS5 stress (Pa)
  real :: REarth
  ! Local variables
  character(len=1000) :: filename
  real, parameter :: gres=0.125

  !Grid specs: lon,lat,dx,dy,area
  filename=trim(trim(data_dir)//'a68_experiment_ll_p125_grid.nc')
  call read_a68_data(mpp_domain,filename,'longitude',lon)
  lon=lon+360
  call read_a68_data(mpp_domain,filename,'latitude',lat)

  call haversine_dist_and_area(REarth,gres,lon,lat,dx,dy,area)

  ! if (Old_wrong_Rearth_and_SSH) then
  !   call read_a68_data(mpp_domain,filename,'dx_Rold',dx)
  !   call read_a68_data(mpp_domain,filename,'dy_Rold',dy)
  !   call read_a68_data(mpp_domain,filename,'area_Rold',area)
  ! elseif (mom_Rearth) then
  !   call read_a68_data(mpp_domain,filename,'dx_mom',dx)
  !   call read_a68_data(mpp_domain,filename,'dy_mom',dy)
  !   call read_a68_data(mpp_domain,filename,'area_mom',area)
  ! else
  !   call read_a68_data(mpp_domain,filename,'dx',dx)
  !   call read_a68_data(mpp_domain,filename,'dy',dy)
  !   call read_a68_data(mpp_domain,filename,'area',area)
  ! endif

  !Depth
  filename=trim(trim(data_dir)//'a68_experiment_gebco_ll_p125.nc')
  call read_a68_data(mpp_domain,filename,'elevation',depth)
  depth=-depth !depth is negative in the data, but needs to be positive in the model

  !Ocean velocity
  !(vertically-integrated over draft of 200 m berg if vert_int_ocean_vel==true, otherwise, surface vel)
  filename=trim(trim(data_dir)//'a68_experiment_ocean_surf_and_mean_vel_fixed_oras5_dec2020_ll_p125.nc')
  !filename=trim(trim(data_dir)//'a68_experiment_ocean_test_vel.nc')

  if (vert_int_ocean_vel) then
    call read_a68_data(mpp_domain,filename,'uo_mean_over_berg_depth',uo)
    call read_a68_data(mpp_domain,filename,'vo_mean_over_berg_depth',vo)
  else
    call read_a68_data(mpp_domain,filename,'uo',uo)
    call read_a68_data(mpp_domain,filename,'vo',vo)
  end if

  if (tau_is_velocity) then
  !Wind velocities -- note this is saved on tauxa and tauya, the wind stress variables.
  !However, these variables are interpreted in the code as wind velocity if
  !tau_is_velocity=.true. in the nml,
    filename=trim(trim(data_dir)//'a68_experiment_wind_vel_jra55_dec2020_ll_p125.nc')
    call read_a68_data(mpp_domain,filename,'ua',tauxa)
    call read_a68_data(mpp_domain,filename,'va',tauya)
  else
    filename=trim(trim(data_dir)//'a68_experiment_wind_stress_oras5_dec2020_ll_p125.nc')
    call read_a68_data(mpp_domain,filename,'taux',tauxa)
    call read_a68_data(mpp_domain,filename,'tauy',tauya)
  endif

  !Effect sea surface height
  filename=trim(trim(data_dir)//'a68_experiment_ssh_oras5_dec2020_ll_p125.nc')
  if (Old_wrong_Rearth_and_SSH) then
    call read_a68_data(mpp_domain,filename,'ssh_bgrid',ssh)
  else
    call read_a68_data(mpp_domain,filename,'SSH',ssh)
  endif
end subroutine a68_prep

subroutine haversine_dist_and_area(REarth,gres,lon1,lat1,dx,dy,area)
  !Arguments
  real, intent(in) :: REarth    !Earth radius (m)
  real, intent(in) :: gres      !grid resolution (degrees)
  real, intent(in) :: lon1(:,:) !grid longitude
  real, intent(in) :: lat1(:,:) !grid latitude
  real, intent(out) :: dx(:,:)  !grid dx (m)
  real, intent(out) :: dy(:,:)   !grid dy (m)
  real, intent(out) :: area(:,:) !grid area (m^2)
  ! Local variables
  real, parameter :: pi_180=pi/180.
  real, dimension(size(lon1,1),size(lon1,2)) :: lon2,lat2,p1,p2,dp,dm,a,c

  !calculate dx
  lon2=lon1-gres
  lat2=lat1
  p1 = lat1 * pi_180; p2 = lat2 * pi_180; dp = (lat2-lat1) * pi_180; dm = (lon2-lon1) * pi_180
  a = sin(0.5*dp)**2. + cos(p1) * cos(p2) * sin(0.5*dm)**2.; c = 2. * atan2(sqrt(a), sqrt(1-a))
  dx = REarth * c

  !calculate dy
  lat2=lat1-gres
  lon2=lon1
  p1 = lat1 * pi_180; p2 = lat2 * pi_180; dp = (lat2-lat1) * pi_180; dm = (lon2-lon1) * pi_180
  a = sin(0.5*dp)**2. + cos(p1) * cos(p2) * sin(0.5*dm)**2.; c = 2. * atan2(sqrt(a), sqrt(1-a))
  dy = REarth * c

  !calculate area
  area=pi_180*(REarth**2.)*abs(sin(lat1*pi_180)-sin(lat2*pi_180))*abs(gres)
end subroutine haversine_dist_and_area


!> Read ocean and atmospheric data for A68 experiment from file
subroutine read_a68_data(domain,filename,var,field)
  ! Arguments
  type(domain2D) :: domain !< Parallel decomposition
  character(len=*), intent(in) :: filename, var
  real :: field(:,:)

  call mpp_sync()
  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'icebergs_driver: Reading ',var
  if (file_exist(filename)) then
    if (field_exist(filename, var)) then
      call read_data(filename, var, field, domain)
    else
      if (mpp_pe().eq.mpp_root_pe()) then
        call error_mesg('icebergs_driver:','Variable missing in file!', FATAL)
      endif
    endif
  else
    if (mpp_pe().eq.mpp_root_pe()) then
      call error_mesg('icebergs_driver:','File for variable not found!', FATAL)
    endif
  endif
  call mpp_sync()
end subroutine read_a68_data

subroutine a68_prep_3d(data_dir,mpp_domain,tauxa_hr,tauya_hr,uo_hr,vo_hr,ssh_hr)
  ! Arguments
  character(len=*) :: data_dir

  type(domain2D) :: mpp_domain
  real :: tauxa_hr(:,:,:) !< Zonal jra55 wind velocity (m/s, hourly)
  real :: tauya_hr(:,:,:) !< Meridional jra55 wind velocity (m/s, hourly)
  real :: uo_hr(:,:,:) !< Zonal OSCAR ocean surface velocity (m/s, hourly)
  real :: vo_hr(:,:,:) !< Meridional OSCAR ocean surface velocity (m/s, hourly)
  real :: ssh_hr(:,:,:) !< Effective sea surface height (m, hourly)
  ! Local variables
  character(len=1000) :: filename

  !filename=trim(trim(data_dir)//'a68_experiment_wind_vel_jra55_dec2020_HOURLY_ll_p125.nc')
  filename=trim(trim(data_dir)//'a68_experiment_wind_vel_ncep_10m_dec2020_HOURLY_ll_p125.nc')
  call read_a68_data_3d(mpp_domain,filename,'ua',tauxa_hr)
  call read_a68_data_3d(mpp_domain,filename,'va',tauya_hr)

  filename=trim(trim(data_dir)//'a68_experiment_ocean_surf_vel_oscar_dec2020_HOURLY_ll_p125.nc')
  call read_a68_data_3d(mpp_domain,filename,'uo',uo_hr)
  call read_a68_data_3d(mpp_domain,filename,'vo',vo_hr)

  filename=trim(trim(data_dir)//'a68_experiment_ssh_duacs_dec2020_HOURLY_ll_p125.nc')
  call read_a68_data_3d(mpp_domain,filename,'SSH',ssh_hr)
end subroutine a68_prep_3d

subroutine read_a68_data_3d(domain,filename,var,field)
  ! Arguments
  type(domain2D) :: domain !< Parallel decomposition
  character(len=*), intent(in) :: filename, var
  real :: field(:,:,:)

  call mpp_sync()
  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'icebergs_driver: Reading ',var
  if (file_exist(filename)) then
    if (field_exist(filename, var)) then
      call read_data(filename, var, field, domain)
    else
      if (mpp_pe().eq.mpp_root_pe()) then
        call error_mesg('icebergs_driver:','Variable missing in file!', FATAL)
      endif
    endif
  else
    if (mpp_pe().eq.mpp_root_pe()) then
      call error_mesg('icebergs_driver:','File for variable not found!', FATAL)
    endif
  endif
  call mpp_sync()
end subroutine read_a68_data_3d

end module driver_data
