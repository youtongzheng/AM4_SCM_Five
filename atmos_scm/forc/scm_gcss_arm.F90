module scm_gcss_arm_mod 

   use            mpp_mod, only:  mpp_pe, mpp_root_pe, stdlog, mpp_chksum
   use         mpp_io_mod, only:  mpp_open, MPP_RDONLY
   use            fms_mod, only:  write_version_number, open_file,        &
                                check_nml_error,  &
                                FILE_EXIST, ERROR_MESG,     &                               
                                CLOSE_FILE, FATAL,  NOTE,   &
                                read_data, write_data,      &
                                mpp_error
#ifdef INTERNAL_FILE_NML
   USE              mpp_mod, ONLY: input_nml_file
#else
   USE              fms_mod, ONLY: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod
   use   constants_mod, only:  rdgas, cp_air, tfreeze, kappa, &
                               hlv, rvgas, grav, stefan, PI, SECONDS_PER_DAY

   use   surface_flux_mod, only:  surface_flux

   use  scm_utils_mod, only:  us_std_atm, thetal_to_temp

   use vert_advection_mod, only:  vert_advection, SECOND_CENTERED, &
                                  FOURTH_CENTERED, FINITE_VOLUME_LINEAR, &
                                  FINITE_VOLUME_PARABOLIC, &
                                  SECOND_CENTERED_WTS, FOURTH_CENTERED_WTS, &
                                  ADVECTIVE_FORM
   use   field_manager_mod, only: MODEL_ATMOS
   use  tracer_manager_mod, only: get_tracer_index

use             fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                               endlon, rlonb, rlatb,  cold_start, ncnst, &
                               pnats, consv_te, ptop, fv_init, fv_domain, &
                               fv_end, change_time, p_var, restart_format, area, &
                               ak, bk, lon,lat, rlon, rlat, ng_d, f_d, nt_prog, get_eta_level

use  sat_vapor_pres_mod,only : sat_vapor_pres_init   ! ZNT 02/20/2020

   implicit none
!   include 'netcdf.inc'

   public       gcss_arm_data_read,                 &
                gcss_arm_forc_diagnostic_init,      &
                gcss_arm_forc_init,                 &
                update_gcss_arm_forc,               &
                gcss_arm_forc_end,                  &
                get_gcss_arm_flx,                   &
                gcss_arm_sfcflx,                    &
                gcss_arm_snd,                       &
                diag_ustar,                         &
                gcss_arm_surface_flux_loop    ! ZNT 05/20/2020

character(len=12) :: mod_name = 'scm_gcss_arm'
character(len=7) :: mod_name_diag = 'forcing'

real, public, allocatable, dimension(:,:,:)  :: u_geos, v_geos


!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following subroutines:
!
!            arm_data_read     reads forcing data files and initializes
!                              any needed parameters of the run
!            arm_forc_init     initializes prognostic variables:           
!                              pt,u,v,qv,ql,qi
!            arm_forc_end      deallocates allocated spaces
!            update_arm_forc   a call to this routine returns the 
!                              forcing parameters needed.
!            interp_time       finds the index of the TIME_TYPE whose
!                              time is just less than and just greater 
!                              than the time being interpolated to
!            interp_2d_field   interpolates a 2d field given as a 2 1d 
!                              vectors, assumed to be dimensioned over 
!                              pressure, to the desired time and level.
!            rem_miss_var      removes all missing data, indicated by 
!                              being less than a specified input value, 
!                              and replaces it with the value derived 
!                              from linear interpolation in time at a  
!                              given level
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              GLOBAL STORAGE VARIABLES
!
!  FORCING DATA (I.E. OBSERVED FIELDS)
!
!  STANDARD ATMOSPHERE PROFILES
!  
!       plev_stand     pressure levels for standard atmosphere (Pascals)
!       t_stand        temperature for standard  atmosphere (K)
!       q_stand        specific humidity of vapor (kg vapor/kg air)
!       O3_stand       specific humidity of ozone (kg ozone/kg air)
!
!  COORDINATE FIELDS
!
!       plev_forc      pressure levels for observed data (Pascals) 
!       time_layer     time_type variable containing times of layer data
!       time_sfc       time_type variable containing times of sfc data
!       Time_var_start time type variable containing time of Variational
!                      dataset start
!
!  LAYERED FIELDS
!
!       t_forc         mean T of observed data (K)
!       q_forc         mean specific humidity of observed data 
!                      (kg vapor/kg air)
!       u_forc         mean zonal wind (m/s)
!       v_forc         mean meridional wind (m/s)
!       omega_forc     vertical p-velocity for data (Pascals/second)
!       tadv_tot_forc  total temperature advection (K/s)
!       qadv_tot_forc  total specific humidity advection (kg/kg/sec)
!       tadv_hor_forc  horizontal temperature advection (K/s)
!       qadv_hor_forc  horizontal specific humidity advection(kg/kg/sec)
!       tau_inv_forc   timescale for advection (1/sec)
!
!       Q1_obs         Q1 Heating in variational data (K/day)
!       Q2_obs         Q2 Drying in variational data (g/kg/day)
!
!       arscl          ARSCL Cloud Fraction (fraction)
!
!  SURFACE DATA
!
!       ps_forc        surface pressure (Pa)
!       dps_dt_forc    surface pressure tendency (Pascals/second)
!       sen_flux_forc  surface sensible heat flux (W/m2)
!       lat_flux_forc  surface latent heat flux (W/m2)
!       tair_surf_forc surface air temperature (K)
!
!      (more diagnostics from Minghua Zhang's surface dataset)
!
!       precip         precipitation (kg/m2/sec)
!       tground        ground temperature (K)
!       rh             relative humidity (fraction)
!       uwnd           zonal wind (m/s)
!       vwnd           meridional wind (m/s)
!       netraddnsfc    net downward radiation at the surface (W/m2)
!       olr            longwave up at the toa (W/m2)
!       swup_toa       upward shortwave at toa (W/m2)
!       swdn_toa       downward shortwave at toa (W/m2)
!       low_cld_amt    low cloud amount (percent)
!       mid_cld_amt    mid cloud amount (percent)
!       hgh_cld_amt    high cloud amount (percent)   
!       tot_cld_amt    total cloud amount (percent)
!       cld_thk        cloud thickness (km)
!       cld_hgt        cloud top height (km)
!       s_surface      surface value of dry static energy (J/kg)
!       qs_surface     saturated value of surface specific humidity (kg/kg)
!       wvp            water vapor path (kg/m2)
!       lwp            liquid water path (kg/m2)
!       qtend          total tendency in column water vapor (kg/m2/sec)
!       advqtend       total advective tendency in 
!                      column water vapor (kg/m2/sec)
!       evap           evaporative flux (kg/m2/sec)
!       ttend          total tendency in energy (W/m2)
!       advttend       total advective tendency in energy (W/m2)
!       radttend       total radiative tendency in energy (W/m2)
!       latttend       total latent heating  (W/m2)
!
!       swdn_sfc       downward shortwave at the surface (W/m2)
!       swup_sfc       upward shortwave at the surface (W/m2)
!       lwdn_sfc       downward longwave at the surface (W/m2)
!       lwup_sfc       upward longwave at the surface (W/m2)
!
!       tskin_forc     tskin provided by the surface dataset (Kelvin)
!
!       alb_24hr       surface albedo provided by the surface dataset
!                      (dimensionless)
!       alb_day        surface albedo provided by the surface dataset
!                      (dimensionless)
!
!       Note that alb_24hr has values at all times of day even night,
!       but alb_day has values only during the day.
!

!----------Diagnostic data----------------------------------------------
     integer          ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                                                &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,     id_qndt_vadv,     & !h1g
            id_udt_vadv, id_udt_geos, id_udt_lf, id_udt_nudge,          &
            id_vdt_vadv, id_vdt_geos, id_vdt_lf, id_vdt_nudge,          &
            id_pf_forc, id_ph_forc, id_zf_forc, id_zh_forc,                  &
            id_u_geos, id_v_geos
  
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         INITIALIZE PARAMETERS OF THE MODULE
!
!       ivert_stand  # of vertical levels to standard atmosphere
!       ivert        # of vertical levels of VAR data

!
!       p00          reference pressure (Pa)
!       rough_mom    momentum roughness length characteristic for ARM 
!                    site taken from Andy Brown's set up for ARM 
!                    shallow cu intercomparison experiment
!       rough_heat   heat roughness length = rough_mom/7.4
!       TSKIN        skin temperature (K)
!       SENFLUX      sensible heat flux (W/m2)
!       EVAPFLUX     evaporation flux (kg water/m2/sec)
!       ALBEDO_OBS   observed albedo (dimensionless)
!
!       FORCMETHOD  letter indicating which forcing method is being used
!
!       FORCMETHOD            SPECIFIED FORCING
!       ----------     -----------------------------------
!           D          revealed or total advective forcing
!           F               revealed plus nudging
!
!       do_qc_adv_loss            Should condensate be horizontally
!                                 advected out of the grid box without
!                                 no influx?
!
!       tracer_vert_advec_scheme  Which advection scheme should be used?
!
!       temp_vert_advec_scheme    Which advection scheme should be used?
!
!                                 SECOND_CENTERED         1
!                                 FOURTH_CENTERED         2
!                                 FINITE_VOLUME_LINEAR    3 
!                                 FINITE_VOLUME_PARABOLIC 4
!                                 SECOND_CENTERED_WTS     5 
!                                 FOURTH_CENTERED_WTS     6
!
!       vert_advec_cond           Should condensate be vertically
!                                 advected?
!
!       do_iop_mean_albedo        If do_specified_albedo is true, should
!                                 one use the IOP mean surface albedo as
!                                 a constant or use the time varying 
!                                 observed albedo
!
!       wind_relaxation_tau       relaxation time scale for horizontal
!                                 winds (seconds)
!
!       relaxation_tau            relaxation time scale for temperature
!                                 and water vapor (seconds)
!
!       use_given_omega_surface   should you use the provided 
!                                 omega_surface? In many cases you 
!                                 should not as it is zero.

REAL,    PRIVATE               :: p00 = 1000.e2

REAL,    PUBLIC                :: missing_value = -999.
LOGICAL, PUBLIC                :: do_netcdf_restart = .TRUE.

character(len=64)              :: configuration = 'unused'

integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3

logical, public                :: do_geo = .true.

real,    private               :: psfc = 970.e2  
real,    private               :: zsfc = 0.0

logical, private               :: initialized = .false.
integer, dimension(6)          :: init_date = (/ 0, 0, 0, 0, 0, 0 /)
integer, dimension(6)          :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
integer                        :: dt_atmos = 60
integer, save                  :: current_days0,  current_sec0
integer                        :: current_days,  current_sec

type(time_type)                :: Time_init
type(time_type)                :: Time_current
type (time_type)               :: Time_step_atmos

character(len=15)              :: rad_scheme = "simplified"

!--------------------- version number ----------------------------------
!
character(len=128) :: Version = '$Id$'
character(len=128) :: Tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)

NAMELIST /SCM_GCSS_ARM_NML/                                 &
                         tracer_vert_advec_scheme,          &
                         temp_vert_advec_scheme,            &
                         momentum_vert_advec_scheme,        &
                         do_geo,                            &
                         rad_scheme,                        & 
                         init_date, current_date, dt_atmos

LOGICAL :: arm_forc_initialized = .FALSE.

real, parameter   :: d622 = rdgas/rvgas
real, parameter   :: d378 = 1.0-d622
real, parameter   :: d608 = d378/d622

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

! ---> h1g, 2010-10-12
integer ::    id_qvdt_forc_col
integer ::    id_qldt_vadv_col
! <--- h1g, 2010-10-12


CONTAINS


!############################################################
subroutine gcss_arm_data_read(kmax)

implicit none

integer,  intent (in)                  :: kmax
integer                                :: i
integer                                :: unit,ierr,io,logunit

character*64                           :: fname_res='INPUT/gcss_arm.res.nc'

integer :: year, month, day

if (initialized) return
initialized = .true.

! if the constructor is already called, then return, otherwise 
! set arm_forc_initialized .TRUE.
#ifdef INTERNAL_FILE_NML
   READ (input_nml_file, nml=scm_gcss_arm_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_gcss_arm_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1; 
      do while (ierr /= 0)
         read  (unit, nml=scm_gcss_arm_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_gcss_arm_nml')
      enddo
10    call close_file (unit)
   endif
#endif

Time_init =set_date(init_date(1), init_date(2), init_date(3), init_date(4), init_date(5), init_date(6) )
call get_time(Time_init, current_sec0, current_days0)

Time_current=set_date(current_date(1), current_date(2),current_date(3), current_date(4),current_date(5), current_date(6))
call get_time(Time_current, current_sec, current_days)

Time_step_atmos = set_time (dt_atmos,0)
! in get_gcss_arm_flx,   Time_current will added 1 time step
Time_current = Time_current - Time_step_atmos

!--------- write version number and namelist --------

call write_version_number ( version, tagname )
if(mpp_pe() == mpp_root_pe() ) then
  logunit =stdlog()
  write(logunit,nml=scm_gcss_arm_nml)
endif

!--------- allocate memory ---------
if (allocated(u_geos)) deallocate(u_geos);  allocate(u_geos(1,1,kmax))
if (allocated(v_geos)) deallocate(v_geos);  allocate(v_geos(1,1,kmax))

! ---> h1g, 2010-10-12
if (file_exist('INPUT/gcss_arm.res.nc') ) then
    if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_gcss_arm_mod', &
         'Reading netCDF formatted restart file: INPUT/gcss_arm.res.nc', NOTE)
    if( allocated(u_geos) ) call read_data (fname_res, 'u_geos',  u_geos)
    if( allocated(v_geos) ) call read_data (fname_res, 'v_geos',  v_geos)
endif
! <--- h1g, 2010-10-12

end subroutine  gcss_arm_data_read


!#######################################################################
subroutine gcss_arm_forc_init(time_interp, As, Bs)
#include "fv_arrays.h"
!      VARIABLES
!
!      ------
!      INPUT:
!      ------
!
!      time_interp  time_type variable containing time we are
!                   interpolating to.
!      pdamp        maximum pressure for which profile is set to a 
!                   reference midlatitude profile 
!      As, Bs       A's and B's of half levels in hybrid coordinate
!                   phalf(k) = A(k) + B(k) * psurf
!
!      -------------
!      INPUT/OUTPUT:
!      -------------
!
!      pt      temperature in degrees Kelvin 
!      u       zonal wind (m/s)
!      v       meridional wind (m/s)
!      qv      water vapor specific humidity (kg water vapor/kg air)
!      ql      liquid water condensate specific humidity 
!              (kg condensate/ kg air)
!      qi      ice H2O condensate specific humidity 
!              (kg condensate/kg dry air)
!      qa      cloud fraction (fraction)
!      qn          cloud droplet number concentration  !h1g
!


TYPE(TIME_TYPE)                          :: time_interp
REAL,  INTENT (IN)   , DIMENSION(:)      :: As,Bs
REAL,  DIMENSION(size(pt,1),size(pt,2))    :: elev

!  Internal variables
!  ------------------

integer  :: kdim, k

integer, parameter :: itmax = 10
real, dimension(size(pt,3)+1) :: eta, peta
real, dimension(size(pt,3))   :: T_us_std, qv_us_std
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, zh, zhnew
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf, zf, zfnew

! ZNT 02/20/2020: Test
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph0
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf0


real :: u_snd, v_snd, ug_snd, vg_snd, thetal_snd, qt_snd, T_snd, qv_snd, ql_snd

real maxerror

integer i
integer :: j
#include "fv_point.inc"
       nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
       nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
       nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
       nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
       nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   !h1g
       nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )   !h1g

!   Initialize surface pressure and topography

ps   = psfc
elev = zsfc  !phis ! ZNT 03/20/2020

!   ZNT 02/20/2020: initialize vapor pressure table
    call sat_vapor_pres_init

!   Setup temporary vertical grid: this grid is needed in order to use

kdim = size(pt,3)

       ! --- Create delp (from hydro_eq in init_dry_atm.F90)
       do k=1,KDIM
         do j=1,size(ps,2)
           do i=1,size(ps,1)
             delp(i,j,k) = As(k+1)-As(k) + ps(i,j)*(Bs(k+1)-Bs(k))
           enddo
         enddo
       enddo
       call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,     &
              pe, peln,  pk,  pkz,  kappa, q, ng_d, ncnst, .false. )

       ! ZNT 02/20/2020: Note - '.false.' means no dry air mass adj,
       !     thus the calculations are consistent, given Bs(1)=0, and
       !     (ak,bk) from fv_pack = (As,Bs) as input, and
       !     i.e., ph (with get_eta level) = pe = p_half (from compute_p_z), 
       !           pf (with get_eta_level) = p_full (with compute_p_z).
       ! 
       ! p_var: pe(1) = ptop = ak(1); 
       !        pe(k) = pe(k-1) + delp(k-1)
       !              = pe(k-1) + As(k)-As(k-1) + ps*(Bs(k)-Bs(k-1))
       !              = ... = pe(1) + As(k)-As(1) + ps*(Bs(k)-Bs(1))
       !              = As(k) + Bs(k)*ps + (ak(1)-As(1)-Bs(1)*ps)
       ! compute_p_z: p_half = pe; p_full(k) = delp(k)/(peln(k+1)-peln(k))
       ! 
       ! get_eta_level: ph(1) = ak(1); ph(k) = ak(k) + bk(k)*ps
       !                pf(k) = (ph(k+1)-ph(k))/log(ph(k+1)/ph(k))
       !                (ak(1) > ptop_min = 1e-6 is always satisfied).

       do j=1,size(ps,2)
         do i=1,size(ps,1)
           ! call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
           call get_eta_level(nlev, ps(i,j) , pf0(i,j,:), ph0(i,j,:))  ! ZNT
         enddo
       enddo

!   Initialize pt, qv first with US standard atmosphere

do k=1,kdim
  ! call us_std_atm( pf(1,1,k), T_us_std(k), qv_us_std(k) )
  call us_std_atm( pf0(1,1,k), T_us_std(k), qv_us_std(k) )
  pt(:,:,k)  = T_us_std(k)

  q(:,:,k,nsphum) = qv_us_std(k)
end do
q(:,:,:,nql) = 0.0
ua  = 0.0
va  = 0.0

!   Compute heights
    call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zf, zh)
    ! ZNT 02/20/2020: DEBUG output
    write (*,*) 'PF0 - PF', pf0 - pf
    write (*,*) 'PH0 - PH', ph0 - ph
    write (*,*) 'ZF', zf
    write (*,*) 'ZH', zh
    ! ZNT 02/20/2020: DEBUG result: ph0 == ph, pf0 ~= pf (differ by 1e-9)

    !  call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zf, zh )
    !        call error_mesg ('scm_gcss_arm',   &
    !        'need to compute height', FATAL )

!   Iteration to compute ARM sounding
i=0
maxerror = 1.0
do while ( maxerror > 0.001 .and. i < itmax )
  i = i + 1
  do k=1,kdim
    call gcss_arm_snd( zf(1,1,k), u_snd, v_snd, ug_snd, vg_snd, thetal_snd, qt_snd )
    ua(:,:,k) = u_snd
    va(:,:,k) = v_snd
    u_geos(:,:,k) = ug_snd
    v_geos(:,:,k) = vg_snd
    call thetal_to_temp( thetal_snd, qt_snd, pf(1,1,k), T_snd, qv_snd, ql_snd )
    pt(:,:,k)  = max ( T_snd, 200.0 )
    q(:,:,k,nsphum) = qv_snd
    q(:,:,k,nql) = ql_snd
  enddo
     call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zfnew, zhnew)

     ! call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zfnew, zhnew )
     !       call error_mesg ('scm_gcss_arm',   &
     !       'need to compute height', FATAL )

  maxerror = maxval( abs( zfnew - zf ) )
  zf = zfnew
  zh = zhnew
enddo

    ! ZNT 03/31/2020: DEBUG output
    write (*,*) 'ZFnew', zf
    write (*,*) 'ZHnew', zh
    write (*,*) 'THL1', thetal_snd
    write (*,*) 'QT1', qt_snd
    write (*,*) 'TEMP', pt(1,1,:)
    write (*,*) 'QV', q(1,1,:,nsphum)


if ( i >= itmax ) then
    call error_mesg('gcss_arm_forc_init',  &
                    'failed to converge while creating blended sounding', FATAL)
endif

! Initialize cloud amount

where ( q(:,:,:,nql) > 0.0 )
    q(:,:,:,nqa) = 1.0
elsewhere
    q(:,:,:,nqa) = 0.0
endwhere
if( nqn > 0 ) q(:,:,:,nqn) = 0.0
 
    ! ZNT 02/23/2020 - Note: Also need to initialize u_srf and v_srf
    u_srf(:,:)=ua(:,:,KDIM)
    v_srf(:,:)=va(:,:,KDIM)

return 
end subroutine gcss_arm_forc_init

!#######################################################################
subroutine gcss_arm_forc_diagnostic_init(axes, Time)

implicit none

integer, dimension(3) :: half = (/1,2,4/)
integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time
   
! --- initialize axes -------------------------------------------------!
id_pf_forc = register_diag_field (mod_name_diag, 'pf_forc', axes(1:3), Time, &
     'Pressure at full level', 'hPa',  missing_value = missing_value)

id_ph_forc = register_diag_field (mod_name_diag, 'ph_forc', axes(half), Time, &
     'Pressure at half level', 'hPa',  missing_value = missing_value)

id_zf_forc = register_diag_field (mod_name_diag, 'zf_forc', axes(1:3), Time, &
     'Height at full level', 'm',  missing_value = missing_value)

id_zh_forc = register_diag_field (mod_name_diag, 'zh_forc', axes(half), Time, &
     'Height at half level', 'm',  missing_value = missing_value)

id_tdt_radf = register_diag_field (mod_name_diag, 'tdt_radf', axes(1:3), Time, &
     'Temperature tendencies due to radiative forcing', 'K/s',  missing_value = missing_value)

id_tdt_vadv = register_diag_field (mod_name_diag, 'tdt_vadv', axes(1:3), Time, &
     'Temperature tendencies due to vertical advection', 'K/s', missing_value = missing_value)

id_udt_vadv = register_diag_field (mod_name_diag, 'udt_vadv', axes(1:3), Time, &
     'U tendencies due to vertical advection', 'm/s2', missing_value = missing_value)

id_vdt_vadv = register_diag_field (mod_name_diag, 'vdt_vadv', axes(1:3), Time, &
     'V tendencies due to vertical advection', 'm/s2', missing_value = missing_value)

id_udt_geos = register_diag_field (mod_name_diag, 'udt_geos', axes(1:3), Time, &
     'U tendencies due to geostrophic wind', 'm/s2', missing_value = missing_value)

id_vdt_geos = register_diag_field (mod_name_diag, 'vdt_geos', axes(1:3), Time, &
     'V tendencies due to geostrophic wind', 'm/s2', missing_value = missing_value)

id_udt_nudge = register_diag_field (mod_name_diag, 'udt_nudge', axes(1:3), Time, &
     'U tendencies due to nudging', 'm/s2', missing_value = missing_value)

id_vdt_nudge = register_diag_field (mod_name_diag, 'vdt_nudge', axes(1:3), Time, &
     'V tendencies due to nudging', 'm/s2', missing_value = missing_value)

id_qvdt_vadv = register_diag_field (mod_name_diag, 'qvdt_vadv', axes(1:3), Time, &
     'Vapor tendencies due to vertical advection', 'kg/kg/s', missing_value = missing_value)

id_qldt_vadv = register_diag_field (mod_name_diag, 'qldt_vadv', axes(1:3), Time, &
     'liquid tendencies due to vertical advection', 'kg/kg/s', missing_value = missing_value)

id_qidt_vadv = register_diag_field (mod_name_diag, 'qidt_vadv', axes(1:3), Time, &
     'ice tendencies due to vertical advection', 'kg/kg/s', missing_value = missing_value)

id_qadt_vadv = register_diag_field (mod_name_diag, 'qadt_vadv', axes(1:3), Time, &
     'cloud amount tendencies due to vertical advection', '1/s', missing_value = missing_value)

id_qndt_vadv = register_diag_field (mod_name_diag, 'qndt_vadv', axes(1:3), Time, &
     'liquid droplet number concentration tendencies due to vertical advection', '1/cm3/s', missing_value = missing_value)

id_udt_lf = register_diag_field (mod_name_diag, 'udt_lf', axes(1:3), Time, &
     'U tendencies due to large-scale forcing', 'm/s2', missing_value = missing_value)

id_vdt_lf = register_diag_field (mod_name_diag, 'vdt_lf', axes(1:3), Time, &
     'V tendencies due to large-scale forcing', 'm/s2', missing_value = missing_value)

id_tdt_lf = register_diag_field (mod_name_diag, 'tdt_lf', axes(1:3), Time, &
     'Temperature tendencies due to large-scale forcing', 'K/s', missing_value = missing_value)

id_qvdt_lf = register_diag_field (mod_name_diag, 'qvdt_lf', axes(1:3), Time, &
     'Vapor tendencies due to large-scale forcing', 'kg/kg/s', missing_value = missing_value)

id_u_geos = register_diag_field (mod_name_diag, 'u_geos', axes(1:3), Time, &
     'U geostrophic wind', 'm/s',  missing_value = missing_value)

id_v_geos = register_diag_field (mod_name_diag, 'v_geos', axes(1:3), Time, &
     'V geostrophic wind', 'm/s',  missing_value = missing_value)

! ---> h1g, 2010-10-12
 id_qvdt_forc_col =  register_diag_field (mod_name_diag, 'qvdt_forc_col', axes(1:2), Time, &
     'column integrated vapor forcing', 'kg/m2/s',  missing_value = missing_value)
 id_qldt_vadv_col =  register_diag_field (mod_name_diag, 'qldt_vadv_col', axes(1:2), Time, &
     'column integrated cloud water vertical advection', 'kg/m2/s',  missing_value = missing_value)
! <--- h1g, 2010-10-12

end subroutine gcss_arm_forc_diagnostic_init



!#######################################################################
subroutine update_gcss_arm_forc(time_interp,time_diag )
#include "fv_arrays.h"
!       VARIABLES
!
!       ------
!       INPUT:
!      ------
!      time_interp  time for interpolation
!      time_diag    time for diagnostics
!      dt_int       time step
!      Hgrid        horizontal grid
!      Vgrid        vertical grid
!      elev         elevation
!
!       ------------
!       INPUT/OUTPUT
!       ------------
!
!       pfull   pressure in PASCALS at full model levels
!       phalf   pressure in PASCALS at half model levels
!               NOTE that it is assumed that p(j)<p(j+1)
!       pt       temperature in degrees Kelvin 
!       u       zonal wind (m/s)
!       v       meridional wind (m/s)
!       qv      water vapor specific humidity (kg water vapor/kg air)
!       ql      liquid water condensate specific humidity 
!                           (kg condensate/ kg air)
!       qi      ice H2O condensate specific humidity (kg condensate/kg dry air)
!       qa      cloud fraction (fraction)
! 
!
!       ------
!       OUTPUT:
!       ------
!
!       AT      temperature tendency due to large-scale forcing (K/sec)
!       AU      zonal wind  tendency due to large-scale forcing (m/sec/sec)
!       AU meridional wind  tendency due to large-scale forcing (m/sec/sec)
!       AQ      water vapor tendency due to large-scale forcing 
!                           (kg vapor/kg air/sec)
!       AL      liquid water condensate tendency due to large-scale forcing
!                           (kg condensate/ kg air/ sec)
!       AI      H2O ice condensate tendency due to large-scale forcing
!                           (kg condensate/ kg air/ sec)
!       AA      cloud fraction tendency due to large-scale forcing
!                           (fraction / sec)
!       omega_full   omega interpolated to full model levels (Pa/sec)

TYPE(TIME_TYPE), intent(in)              ::  time_interp,  time_diag

integer                                          ::  i, j, k, kdim
logical                                          ::  used
real                                             ::  fcriolis
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf, zf

real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: exner

real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, zh
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dT_adi, dT_vadv, dT_lf, dT_rad, &
                                                       dqv_vadv, dqv_lf, dql_vadv, dqi_vadv, dqa_vadv, dqn_vadv !h1g
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: du_vadv, dv_vadv, du_geos, dv_geos, du_lf, dv_lf, du_nudge, dv_nudge
! ---> h1g, 2010-10-12
real, dimension(size(pt,1),size(pt,2))  ::  qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::  qldt_vadv_col
! <--- h1g, 2010-10-12
real, dimension(size(pt,1),size(pt,2))              :: elev

!  Internal variables
!  ------------------
real       ::    theta_tmp,    rad_tmp,   rt_tmp
real       ::    true_time,      a,   b
integer    ::    current_sec,  current_days
integer    ::   i1, i2

! Constant Parameters
real, parameter, dimension(6) ::  & 
  atheta = (/ 0.000, 0.000,  0.000, -0.080, -0.160, -0.160/), & 
  rtheta = (/-0.125, 0.000,  0.000,  0.000,  0.000, -0.100/), & 
  art    = (/ 0.080, 0.020, -0.040, -0.100, -0.160, -0.300/)

!  ------------------------------------------------------------------------------------------------------------

kdim = size(pt,3)
omga = 0.0

! --- update pf, ph, and zf, zh
elev=zsfc   ! znt 20200226
ps=psfc
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo
       call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zf, zh)

       ! call compute_height    (Vgrid, grav*elev, T, qv, pf, ph, zf, zh)
       !     call error_mesg ('scm_gcss_arm',   &
       !     'need to compute height', FATAL )

exner = ( pf/p00 )**kappa

call get_time( time_interp, current_sec,  current_days)
true_time = real(current_sec)+86400.0*real(current_days-current_days0)

! Interpolate in time to get theta and rt tendency
i1 = floor( ( true_time - 41400. ) / 10800. ) + 1
i1 = min( max( i1, 1 ), 5 )
i2 = i1 + 1

if (i1 < 5) then
   a = ( true_time - (41400. + 10800. * real(i1-1)) ) / 10800.
else
   a = ( true_time - (41400. + 10800. * real(i1-1)) ) / 9000.
end if

if ( trim( rad_scheme ) == "simplified" ) then
    theta_tmp = ( 1. - a ) * ( atheta(i1) ) &
            + a * ( atheta(i2) ) 
    rad_tmp = ( 1. - a ) * ( rtheta(i1) ) & 
            + a * ( rtheta(i2) )

else ! Factor in radiation later
    theta_tmp = ( 1. - a ) * ( atheta(i1) + 0.0 ) & 
            + a * ( atheta(i2) + 0.0 )
    rad_tmp   = 0.0
end if

rt_tmp = ( 1. - a ) * art(i1) + a * art(i2)

! Convert to the right units
theta_tmp = theta_tmp / 3600.
rad_tmp   = rad_tmp / 3600.
rt_tmp    = rt_tmp / ( 3600. * 1000. )

! Interpolate with respect to height

do k = 1, kdim
  select case( int( zf( 1, 1, k ) ) )
  case ( 0:999 )
    dqv_lf( 1, 1, k )  =  rt_tmp  
!  ZNT 05/19/2020: strictly, rt_tmp is the mixing ratio tendency.
!  Using it directly as the specific humidity tendency would result in a small error.

!  ZNT 03/30/2020 - bug corrected
!    dT_lf( 1, 1, k )   =  (theta_tmp + rad_tmp) * exner( 1, 1, k ) 
    dT_lf( 1, 1, k )   =  theta_tmp * exner( 1, 1, k )
    dT_rad( 1, 1, k )  =  rad_tmp * exner( 1, 1, k ) 

!  ZNT 03/30/2020 - bug corrected
!  case ( 1000:2999 )
!    b                 = 1. - ( zf( 1, 1,k) - 1000. ) / 2000.
  case ( 1000:1999 )
    b                 = 1. - ( zf( 1, 1,k) - 1000. ) / 1000.
    dqv_lf( 1, 1, k ) =  b * rt_tmp
!  ZNT 03/30/2020 - bug corrected
!    dT_lf( 1, 1, k )  =  b * ( theta_tmp + rad_tmp ) * exner( 1, 1, k ) 
    dT_lf( 1, 1, k )  =  b * theta_tmp * exner( 1, 1, k ) 
    dT_rad( 1, 1, k ) =  b * rad_tmp * exner( 1, 1, k ) 

  case default
    dqv_lf( 1, 1, k ) = 0.0
    dT_lf( 1, 1, k )  = 0.0
    dT_rad( 1, 1, k ) = 0.0

  end select
end do ! k=1, kdim

! --- large-scale forcing tendencies
du_lf = 0.0
dv_lf = 0.0

! --- large-scale subsidence tendencies
dT_vadv=0.0; dT_adi=0.0; dqv_vadv=0.0; dql_vadv=0.0; dqi_vadv=0.0;dqa_vadv=0.0;dqn_vadv=0.0;!h1g
du_vadv=0.0; dv_vadv=0.0;


! --- nudge winds
du_nudge = 0.0
dv_nudge = 0.0

! --- compute geostrophic tendencies
du_geos=0.; dv_geos=0.;
if (do_geo) then
   ! fcriolis=4.*PI/SECONDS_PER_DAY*sin(rlat(1,1))
   ! ZNT 05/19/2020: fc = 8.55e-5 s-1 for 36 deg lat; Brown 2002 uses fc = 8.5e-5 s-1
   fcriolis=f_d(1,1)
   do k=1, kdim
      du_geos(:,:,k) =   fcriolis * (va(1,1,k)-v_geos(1,1,k))
      dv_geos(:,:,k) = - fcriolis * (ua(1,1,k)-u_geos(1,1,k))
   end do
end if

! --- sum all tendencies
u_dt = du_vadv + du_geos + du_lf + du_nudge
v_dt = dv_vadv + dv_geos + dv_lf + dv_nudge
t_dt = dT_rad + dT_vadv + dT_lf
q_dt(:,:,:,nsphum) = dqv_vadv + dqv_lf
q_dt(:,:,:,nql) = dql_vadv
q_dt(:,:,:,nqi) = dqi_vadv
q_dt(:,:,:,nqa) = dqa_vadv

if( nqn > 0 )  q_dt(:,:,:,nqn) = dqn_vadv

if (id_pf_forc > 0)  used = send_data( id_pf_forc,  pf      (:,:,:), time_diag, 1, 1)

if (id_ph_forc > 0)  used = send_data( id_ph_forc,  ph      (:,:,:), time_diag, 1, 1)

if (id_zf_forc > 0)  used = send_data( id_zf_forc,  zf      (:,:,:), time_diag, 1, 1)

if (id_zh_forc > 0)  used = send_data( id_zh_forc,  zh      (:,:,:), time_diag, 1, 1)

if (id_tdt_vadv > 0) used = send_data( id_tdt_vadv, dT_vadv (:,:,:), time_diag, 1, 1)

if (id_udt_vadv > 0) used = send_data( id_udt_vadv, dU_vadv (:,:,:), time_diag, 1, 1)

if (id_vdt_vadv > 0) used = send_data( id_vdt_vadv, dV_vadv (:,:,:), time_diag, 1, 1)

if (id_udt_geos > 0) used = send_data( id_udt_geos, du_geos (:,:,:), time_diag, 1, 1)

if (id_vdt_geos > 0) used = send_data( id_vdt_geos, dv_geos (:,:,:), time_diag, 1, 1)

if (id_udt_nudge > 0) used = send_data( id_udt_nudge, du_nudge (:,:,:), time_diag, 1, 1)

if (id_vdt_nudge > 0) used = send_data( id_vdt_nudge, dv_nudge (:,:,:), time_diag, 1, 1)

if (id_qvdt_vadv> 0) used = send_data( id_qvdt_vadv,dqv_vadv(:,:,:), time_diag, 1, 1)

if (id_qldt_vadv> 0) used = send_data( id_qldt_vadv,dql_vadv(:,:,:), time_diag, 1, 1)

if (id_qidt_vadv> 0) used = send_data( id_qidt_vadv,dqi_vadv(:,:,:), time_diag, 1, 1)

if (id_qadt_vadv> 0) used = send_data( id_qadt_vadv,dqa_vadv(:,:,:), time_diag, 1, 1)

! ---> h1g, 2010-10-12
   qvdt_forcing_col = 0.0
   do k=1, kdim
         qvdt_forcing_col =  qvdt_forcing_col &
            + ( dqv_vadv( :,:,k )  + dqv_lf( :,:,k ) ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do

   qldt_vadv_col = 0.0
   do k=1, kdim
         qldt_vadv_col =  qldt_vadv_col &
            + (  dql_vadv( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do
! <--- h1g, 2010-10-12

if ( nqn > 0 ) then 
    if (id_qndt_vadv> 0) used = send_data( id_qndt_vadv,dqn_vadv(:,:,:), time_diag, 1, 1)
endif

if (id_tdt_radf > 0) used = send_data( id_tdt_radf, dT_rad  (:,:,:), time_diag, 1, 1)

if (id_udt_lf > 0)   used = send_data( id_udt_lf,   du_lf   (:,:,:), time_diag, 1, 1)

if (id_vdt_lf > 0)   used = send_data( id_vdt_lf,   dv_lf   (:,:,:), time_diag, 1, 1)

if (id_tdt_lf > 0)   used = send_data( id_tdt_lf,   dT_lf   (:,:,:), time_diag, 1, 1)

if (id_qvdt_lf > 0)  used = send_data( id_qvdt_lf,  dqv_lf  (:,:,:), time_diag, 1, 1)

if (id_u_geos > 0)   used = send_data( id_u_geos,   u_geos  (:,:,:), time_diag, 1, 1)

if (id_v_geos > 0)   used = send_data( id_v_geos,   v_geos  (:,:,:), time_diag, 1, 1)

! ---> h1g, 2010-10-12
if ( id_qvdt_forc_col > 0 )  used = send_data(  id_qvdt_forc_col, qvdt_forcing_col(:,:), time_diag, 1, 1 )
if ( id_qldt_vadv_col > 0 )  used = send_data(  id_qldt_vadv_col, qldt_vadv_col(:,:), time_diag, 1, 1 )
! <--- h1g, 2010-10-12
        print*,' forc, ql= ', mpp_chksum( q(:,:,:,nql) ) 
end subroutine update_gcss_arm_forc



!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_gcss_arm_flx( z, dn0, thlm_sfc, tm_sfc, qm_sfc, um_sfc, vm_sfc,  &
                             ustar, flux_t, flux_q, flux_u, flux_v )

implicit none
real, intent(in),  dimension(:)  ::  & 
  z,               &   ! Height at zt(2)       [m]
  dn0,             &   ! Density at zm(1)      [kg/m^3]
  thlm_sfc,        &   ! Theta_l at zt(2)      [K]
  tm_sfc,          &   ! T at zt(2)            [K], ZNT
  qm_sfc,          &   ! qv at zt(2)           [kg/kg], ZNT
  um_sfc,          &   ! um at zt(2)           [m/s]
  vm_sfc               ! vm at zt(2)           [m/s]

real, intent(out), dimension(:) :: ustar, flux_t, flux_q, flux_u, flux_v

! ARM roughness height
real, parameter :: z0 = 0.035  ! momentum roughness height
  
real, parameter :: ubmin = 0.25 ! [m/s]

! Internal variables
real,   dimension( size(dn0, 1)  )   ::  & 
  heat_flx, moisture_flx,       &
  heat_flx2, moisture_flx2,   &
  bflx, ubar, & 
  tvm_sfc  ! ZNT 05/19/2020

real       ::    true_time
integer    ::    current_sec,  current_days
integer    ::    i

! -----------------------------------------------------------------------------------

! Compute heat and moisture fluxes from ARM data in (W/m2)

Time_current = Time_current + Time_step_atmos
call get_time(Time_current, current_sec,  current_days)
true_time = real(current_sec)+86400.0*real(current_days-current_days0)

! write(*,*) z(1), dn0(1), thlm_sfc(1), tm_sfc(1), qm_sfc(1),um_sfc(1),vm_sfc(1)

call gcss_arm_sfcflx( true_time, heat_flx(1), moisture_flx(1) )

flux_t( : )  =  heat_flx(1)
flux_q( : ) =  moisture_flx(1)
 
heat_flx( : )        =  heat_flx(1)
moisture_flx( : ) =  moisture_flx(1)

! Convert heat_flx and moisture_flx to natural units
heat_flx2     = heat_flx / ( cp_air * dn0 )    ! (K m/s)
moisture_flx2 = moisture_flx / ( hlv * dn0 )! (m/s)

! Heat flux in units of (m2/s3) (needed by diag_ustar)
! bflx = grav/thlm_sfc * heat_flx2

! ZNT 05/19/2020: including virtual effect
tvm_sfc = tm_sfc*(1.0+d608*qm_sfc)
bflx = grav/tvm_sfc*(heat_flx2*(1.0+d608*qm_sfc) + &
                 moisture_flx2*d608*tm_sfc)

! Compute ubar
ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )
! Compute ustar

do i= 1, size(dn0, 1)
    ustar( i ) = diag_ustar( z( i ),  bflx( i ),  ubar( i ),  z0 )
enddo

flux_u = -um_sfc * ustar * ustar/ ubar
flux_v = -vm_sfc * ustar * ustar / ubar

! write(*,*) heat_flx2(1), moisture_flx2(1), tvm_sfc(1), bflx(1), ubar(1), ustar(1)

end subroutine get_gcss_arm_flx
!########################################################################



!#######################################################################
!
! To use specified surface fluxes for Single column model
!
! Note that this routine does iteration if do_specified_flux = .true.
! Otherwise this does nothing more than to call surface_flux 
!
! The algorithm is to do newtonian iteration for both the surface 
! temperature (canopy air temperature over land) and the surface
! specific humidity. The iterations are nested with the inner loop
! iterating over temperature.  A fixed number of iterations is 
! permitted. A special condition handle the case of non-convergence
! to the observed surface fluxes. 
!
     
subroutine gcss_arm_surface_flux_loop (                                &
                 t_atm,     q_atm,      u_atm,       v_atm,            &
                 p_atm,     z_atm,      p_surf,      t_surf,           &
                 t_ca,      q_surf,     u_surf,      v_surf,           &
                 rough_mom, rough_heat, rough_moist, rough_scale, gust,&
                 flux_t,    flux_q,     flux_lw,     flux_u,           &
                 flux_v,    cd_m,       cd_t,        cd_q,             &
                 w_atm,     u_star,     b_star,      q_star,           &
                 dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
                 dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
                 dt,        land,       seawater,    avail )

  ! ---- arguments -----------------------------------------------------
  logical, intent(in), dimension(:) :: land,  seawater, avail
  real, intent(in),  dimension(:) :: &
       t_atm,     q_atm,      u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,  &
       rough_mom, rough_heat, rough_moist,  rough_scale, gust
  real, intent(out), dimension(:) :: &
       flux_t,    flux_q,     flux_lw,   flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q
  real, intent(inout), dimension(:) :: q_surf
  real, intent(in) :: dt

  !---------------------------------------------------------------------

  integer :: k, nouter, ninner, npoints
  integer :: niter_max = 20
  logical :: do_outer_loop, do_inner_loop
  real    :: max_flux_error = 0.1
  real    :: bflxs,bfact,chu,cmu
  real, dimension(size(t_atm(:))) :: t_ca_tmp, t_surf_tmp, q_surf_tmp, q_surf_tmp1, flux_lw_tmp
  real, dimension(size(t_atm(:))) :: rhotmp
  logical, dimension(size(t_atm(:))) :: land_flg

  real       ::    true_time
  integer    ::    current_sec,  current_days
  integer    ::    i
  real       ::    shflx, lhflx

   ! -----------------------------------------------------------------------------------

   ! Compute heat and moisture fluxes from ARM data in (W/m2)
   Time_current = Time_current + Time_step_atmos
   call get_time(Time_current, current_sec,  current_days)
   true_time = real(current_sec)+86400.0*real(current_days-current_days0)
   call gcss_arm_sfcflx( true_time, shflx, lhflx )
   lhflx = lhflx/hlv
   land_flg(:) = .true.

   !-------------------------------------------------------!
   ! do first calculation of surface flux

   call surface_flux (                                                 &
                 t_atm,     q_atm,      u_atm,       v_atm,            &
                 p_atm,     z_atm,      p_surf,      t_surf,           &
                 t_ca,      q_surf,     u_surf,      v_surf,           &
                 rough_mom, rough_heat, rough_moist, rough_scale, gust,&
                 flux_t,    flux_q,     flux_lw,     flux_u,           &
                 flux_v,    cd_m,       cd_t,        cd_q,             &
                 w_atm,     u_star,     b_star,      q_star,           &
                 dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
                 dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
                 dt,        land_flg,   seawater,    avail             )  ! ZNT 05/20/2020: force land = .true.


   !----------- initialize iteration variables  -----------------------!
   !-- Note that with specified surface fluxes, q_surf from the land
   !-- model may be garbage; thus initialize with 120% of q_atm
   t_surf_tmp = t_surf
   t_ca_tmp   = t_surf
   q_surf_tmp = 1.2*q_atm  ! q_surf
   flux_lw_tmp = flux_lw


   !----------- outer loop -------------!
   
   do_outer_loop = .true.      
   nouter = 1
   npoints = size(t_atm(:))
      
   do while (do_outer_loop)


     !------ inner loop -------------!

     do_inner_loop = .true.
     ninner = 1

     do while (do_inner_loop)


       !---------------------------
       ! calculate surface flux
       q_surf_tmp1(:) = q_surf_tmp(:)
       call surface_flux (                                       &
           t_atm,     q_atm,      u_atm,       v_atm,            &
           p_atm,     z_atm,      p_surf,      t_surf_tmp,       &
           t_ca_tmp,  q_surf_tmp1, u_surf,      v_surf,           &
           rough_mom, rough_heat, rough_moist, rough_scale, gust,&  
           flux_t,    flux_q,     flux_lw,     flux_u,           &
           flux_v,    cd_m,       cd_t,        cd_q,             &
           w_atm,     u_star,     b_star,      q_star,           &
           dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
           dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
           dt,        land_flg,   seawater,    avail             )  ! ZNT 05/20/2020: force land = .true.      

       !------------------------------
       ! should inner loop be redone? 

       do_inner_loop = .false.
       do k = 1, npoints
         if ( ninner .lt. niter_max .and.                     &
              abs(flux_t(k)-shflx) .gt. max_flux_error ) then

           do_inner_loop = .true.
           t_ca_tmp(k) = max( 200., min( 350.,             &
                t_ca_tmp(k) - ( (flux_t(k)-shflx) /      &
                max(dhdt_surf(k),1.e-10) ) ) )
           t_surf_tmp(k) = max( 200., min( 350.,             &
                t_surf_tmp(k) - ( (flux_t(k)-shflx) /      &
                max(dhdt_surf(k),1.e-10) ) ) )

         end if                            
       enddo        

       ninner = ninner + 1

     enddo  !------------- end of inner loop


     !---------------------------
     ! calculate surface flux
     q_surf_tmp1(:) = q_surf_tmp(:)
     call surface_flux (                                            &
              t_atm,     q_atm,      u_atm,       v_atm,            &
              p_atm,     z_atm,      p_surf,      t_surf_tmp,       &
              t_ca_tmp,  q_surf_tmp1, u_surf,      v_surf,           &
              rough_mom, rough_heat, rough_moist, rough_scale, gust,&
              flux_t,    flux_q,     flux_lw,     flux_u,           &
              flux_v,    cd_m,       cd_t,        cd_q,             &
              w_atm,     u_star,     b_star,      q_star,           &
              dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
              dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
              dt,        land_flg,   seawater,    avail             )  ! ZNT 05/20/2020: force land = .true.  

     !------------------------------
     ! should outer loop be redone?

     do_outer_loop = .false.
     do k = 1, npoints
       if (  nouter .lt. niter_max .and.                           &
             hlv*abs(flux_q(k)- lhflx) .gt. max_flux_error ) then

         do_outer_loop = .true.
         q_surf_tmp(k) = max ( 1.e-08, min(0.05, q_surf_tmp(k)- &
                         ((flux_q(k)-lhflx) / max(-1.0*dedq_atm(k),1.e-10))))
!                         ((flux_q(k)-lhflx) / max(dedq_surf(k),1.e-10))))
         if (mpp_pe() == mpp_root_pe() .and. k == 1) then 
             ! print *, 'q-iter #', nouter, ninner, flux_q(k)*hlv, lhflx*hlv, &
             !          ((flux_q(k)-lhflx) / max(-1.0*dedq_atm(k),1.e-10)), &
             !          q_surf_tmp(k)
         end if

       end if
     enddo

     nouter = nouter + 1

   enddo  !----------- end of outer loop

   !----------------------------
   ! deal with non-convergence
   
   rhotmp = p_atm/(rdgas*t_atm*(1.+d608*q_atm))  

   if ( nouter .gt. niter_max .or. ninner .gt. niter_max ) then

     do k = 1, npoints

       if ( hlv*abs(flux_q(k)- lhflx) .gt. max_flux_error .or.&
                abs(flux_t(k)- shflx) .gt. max_flux_error) then 

         if (mpp_pe() == mpp_root_pe() .and. k == 1) then 
             print *, 'NON-CONVERGENCE ', flux_t(k), shflx, &
                          flux_q(k)*hlv, lhflx*hlv, u_star(k)
         end if

         flux_t(k) = shflx
         flux_q(k) = lhflx

         u_star(k) = sqrt ( sqrt ( flux_u(k)**2. +            &
                                   flux_v(k)**2. + 1.e-20 ) / rhotmp(k))

         q_star(k) = lhflx / rhotmp(k) / u_star(k)

         !-----------------------------------------------------
         ! the buoyancy flux, bflxs, is computed assuming 
         ! unsaturated air at the surface from the formula:
         !
         !  bflxs = (chu*senflux + cmu*evapflux) / rho
         !
         !  units of chu = (1/m)
         !  units of cmu = (kg air/kg water)*(m/s*s)

         bfact     = grav/(t_atm(k)*(1.+d608*q_atm(k)))
         chu       =  (1.+d608*q_atm(k))*bfact/cp_air
         cmu       = d608 * bfact * t_atm(k) 
         bflxs     = (chu*shflx + cmu*lhflx) / rhotmp(k)
         b_star(k) = bflxs / u_star(k)

       end if  ! for outside of convergence bounds
          
     end do ! for flux errors

   end if  ! for lack of convergence

   if (mpp_pe() == mpp_root_pe()) then
   ! print *,  'Input variables for sfc flx calculation: '
   ! print *,  '(t_atm, q_atm, u_atm, v_atm)', &  
   !     t_atm(1),     q_atm(1),      u_atm(1),     v_atm(1)
   ! print *,  '(p_atm, z_atm, t_ca)', & 
   !     p_atm(1),     z_atm(1),      t_ca(1)
   ! print *,  '(p_surf, t_surf, u_surf, v_surf)', & 
   !     p_surf(1),    t_surf(1),     u_surf(1),    v_surf(1)
   ! print *,  '(rough_mom, rough_heat, rough_moist, rough_scale, gust)', & 
   !     rough_mom(1), rough_heat(1), rough_moist(1),  rough_scale(1), gust(1)

   print *, 'desired fluxes  (sensible, latent): ',     &
        shflx, hlv*lhflx
   print *, 'iterated fluxes (sensible, latent): ',     &
        flux_t(1), hlv*flux_q(1)
   print *, 'iterated sfc values (t_ca_tmp, q_surf_tmp): ',     &
        t_ca_tmp(1), q_surf_tmp(1)
   print *, 'iterated fluxes (uflux, vflux): ',     &
        flux_u(1), flux_v(1)
   print *, 'scales (u_star, b_star): ',     &
        u_star(1), b_star(1)
   end if

   ! return original longwave flux
   ! 
   ! Note this is needed so that the upward longwave flux is not affected
   ! by the changed t_surf during iteration.

   flux_lw = flux_lw_tmp

   ! reset q_surf, surface specific humidity
   q_surf = q_atm + flux_q / (rhotmp*cd_q*w_atm)  

!
!-----------------------------------------------------------------------
           
end subroutine gcss_arm_surface_flux_loop


!------------------------------------------------------------------------
subroutine gcss_arm_sfcflx( time, heat_flx, moisture_flx )

!       Description:
!       This subroutine computes surface heat and moisture for a specific time
!       according to GCSS ARM specifications. Flux returned are in (W/m2)
!------------------------------------------------------------------------
implicit none

! Parameter constants
integer, parameter :: ntimes = 7

real, parameter, dimension(ntimes) ::  & 
  times = (/ 41400., 55800., 64800., 68400., & 
             77400., 86400., 93600. /), & 
  ! H and LE specifications
  H  = (/-30.,  90., 140., 140., 100., -10., -10./), & 
  LE = (/  5., 250., 450., 500., 420., 180.,   0./)

! Input variable
real, intent(in) :: time !  Current time [s]

! Output variables
real, intent(out) :: heat_flx, moisture_flx

! Local variables
integer :: i1, i2
real :: a

if ( time <= times(1) ) then
   heat_flx     = H(1)
   moisture_flx = LE(1)
else if ( time >= times(ntimes) ) then
   heat_flx     = H(ntimes)
   moisture_flx = LE(ntimes)
else
   i1 = 1
   do while ( i1 <= ntimes-1 )
      i2 = i1 + 1
      if ( time >= times(i1) .and. time < times(i2) ) then
         a            = (time-times(i1))/(times(i2)-times(i1))
         heat_flx     = ( 1. - a ) * H(i1) + a * H(i2)
         moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
         i1           = ntimes
      end if
      i1 = i2
   end do
end if ! time <= times(1)

write(*,*) 'Time, H, LE', time, heat_flx, moisture_flx

return
end subroutine gcss_arm_sfcflx


!########################################################################
subroutine gcss_arm_snd( z, u, v, ug, vg, thetal, qt )
implicit none

real, intent(in) :: z
real, intent(out) :: u, v, ug, vg, thetal, qt
!------------------------------------------------------------------------

! u, v wind;    ug, vg wind;
u = 10.0
v = 0.0

ug = 10.0
vg = 0.0

! qt (specific total water content), g/kg
if( z <= 50.0 ) then
    qt = 15.2 + (15.17 - 15.2 )/50.0 * z
else if ( z <= 350.0 ) then
          qt = 15.17 + ( 14.98 - 15.17 ) / 300.0 * (z - 50.0)
else if ( z <= 650.0 ) then
          qt = 14.98 + ( 14.8 - 14.98 ) / 300.0 * (z - 350.0)
else if ( z <= 700.0 ) then
          qt = 14.8 + ( 14.7 - 14.8 ) / 50.0 * (z - 650.0)
else if ( z <= 1300.0 ) then
          qt = 14.7 + ( 13.5 - 14.7 ) / 600.0 * (z - 700.0)
else if ( z <= 2500.0 ) then
          qt = 13.5 + ( 3.0 - 13.5 ) / 1200.0 * (z - 1300.0)
else if ( z <= 5500.0 ) then
          qt = 3.0
else
     qt = 3.0 - 8.E-3 * (z - 5500.0)
     qt = max( 0.0, qt )
endif
!      convert units to kg/kg
qt = 0.001 * qt
! ZNT 05/19/2020: Note the moisture profile is already given 
! as mixing ratio in Brown (2002). 
 
! thetal, potential temperature [k]
if( z <= 50.0 ) then
    thetal = 299. + (301.5 - 299. )/50.0 * z
else if ( z <= 350.0 ) then
          thetal = 301.5 + ( 302.5 - 301.5 ) / 300.0 * (z - 50.0)
else if ( z <= 650.0 ) then
          thetal = 302.5 + ( 303.53 - 302.5 ) / 300.0 * (z - 350.0)
else if ( z <= 700.0 ) then
          thetal = 303.53 + ( 303.7 - 303.53 ) / 50.0 * (z - 650.0)
else if ( z <= 1300.0 ) then
          thetal = 303.7 + ( 307.13 - 303.7 ) / 600.0 * (z - 700.0)
else if ( z <= 2500.0 ) then
          thetal = 307.13  + ( 314.0 - 307.13  ) / 1200.0 * (z - 1300.0)
else if ( z <= 5500.0 ) then
         thetal = 314.0 + ( 343.2 - 314.  ) / 3000.0 * (z - 2500.0)
else
         thetal =  343.2 +  0.0E-3 * (z - 5500.)
endif

return
end subroutine gcss_arm_snd
!#######################################################################


real function diag_ustar( z, bflx, wnd, z0 )
implicit none

real, parameter ::   & 
!    grav = 9.81,     & ! Gravitational acceleration     [m/s^2]
    vonk   = 0.4,    & ! Accepted value is 0.40 (+/-) 0.01      [-]
    pi = 3.141592654 ! The ratio of radii to their circumference
!Recommendation replace pi above by value from constants.

real, parameter      :: am   =  4.8   !   "          "         "
real, parameter      :: bm   = 19.3   !   "          "         "

real, intent (in)    :: z             ! height where u locates
real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
real, intent (in)    :: wnd           ! wind speed at z
real, intent (in)    :: z0            ! momentum roughness height

integer :: iterate
real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

lnz   = log( z / z0 )
klnz  = vonk/lnz
c1    = pi / 2.0 - 3.0*log( 2.0 )

ustar =  wnd*klnz
!      if (bflx /= 0.0) then
if (abs(bflx) > 1.e-6) then
!      if (abs(bflx) > 1.e-4) then
  do iterate=1,4
!          lmo   = -bflx * vonk/(ustar**3 + eps)
    lmo   = -ustar**3 / ( vonk * bflx )
    zeta  = z/lmo
    if (zeta > 0.) then
      if ( zeta > 1.e10 ) then ! -dschanen UWM for large zeta
        ustar = 1e-10
        exit
      else
        ustar =  vonk*wnd  /(lnz + am*zeta)
        ! ZNT: this iteration seems problematic for the stable branch. Need to change.
      end if
    else
      x     = sqrt( sqrt( 1.0 - bm*zeta ) )
      psi1  = 2.*log( 1.0+x ) + log( 1.0+x*x ) - 2.*atan( x ) + c1
      ustar = wnd*vonk/(lnz - psi1)
    end if
  end do ! 1..4
end if

diag_ustar = ustar

return
end function diag_ustar


!#######################################################################

subroutine gcss_arm_forc_end ()
character*64                 :: fname_res='RESTART/gcss_arm.res.nc'

if (.not.initialized) return
initialized = .false.

! ---> h1g, 2010-10-12
if( do_netcdf_restart ) then
    if (mpp_pe() == mpp_root_pe()) then
        call mpp_error ('scm_gcss_arm_mod', 'Writing netCDF formatted restart file: RESTART/gcss_arm.res.nc', NOTE)
    endif
    call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)	   
    call write_data (fname_res,  'u_geos',  u_geos)
    call write_data (fname_res,  'v_geos',  v_geos)
endif
! <--- h1g, 2010-10-12

deallocate (u_geos, v_geos)
 
end subroutine gcss_arm_forc_end

!########################################################################

! ZNT 02/20/2020: Note - Existing code for calculating z
!
! ------ atmos_fv_dynamics/tools/fv_diagnostics.F90, Line 653 ------
! Compute height at layer edges
!         do i=1,im
!            wz(i,j,km+1) = phis(i,j) * ginv
!         enddo
!
!         do k=km,1,-1
!            do i=1,im
! #ifdef SW_DYN
!                wz(i,j,k) = wz(i,j,k+1) + gg*(pk(i,j,k+1)-pk(i,j,k))
! #else
!                wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))   &
!                     *(peln(i,k+1,j)-peln(i,k,j))
! #endif
!             enddo
!          enddo
!       end do
!
! ------ atmos_scm/driver/coupled/atmosphere.F90, fv_compute_p_z ------
!  if (hydrostatic ) then
!       do k=npz,1,-1
!         do j=1,size(phis,2)
!           do i=1,size(phis,1)
!             tvm = rrg*pt(i,j,k)*(1.+zvir*q_sph(i,j,k))
!             p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
!             z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
!             z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(i,k+1,j)-peln(i,k,j))
!           enddo
!         enddo
!       enddo
!  endif

subroutine compute_p_z (npz, phis, pe, peln, delp, pt, q_sph, p_full, p_half, z_full, z_half)
  ! ZNT 02/20/2020: adapted from the hydrostatic version of 'fv_compute_p_z' subroutine in atmosphere.F90
     integer, intent(in)  :: npz
     real, dimension(:,:),   intent(in)  :: phis
     real, dimension(:,:,:), intent(in)  :: pe, peln, delp, pt, q_sph
     real, dimension(:,:,:), intent(out) :: p_full, p_half, z_full, z_half
  !--- local variables
     integer i,j,k
     real    tvm
     real    :: zvir, rrg, ginv

     zvir = rvgas/rdgas - 1.
     ginv = 1./ grav
     rrg  = rdgas / grav

  !----------------------------------------------------
  ! Compute pressure and height at full and half levels
  !----------------------------------------------------
     z_half(:,:,npz+1) = phis(:,:) * ginv

     do k=1,npz+1
       do j=1,size(phis,2)
          do i=1,size(phis,1)
            p_half(i,j,k) = pe(i,k,j)
         enddo
       enddo
     enddo

     do k=npz,1,-1
       do j=1,size(phis,2)
         do i=1,size(phis,1)
           tvm = rrg*pt(i,j,k)*(1.+zvir*q_sph(i,j,k))
           p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
           z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
           z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(i,k+1,j)-peln(i,k,j))
         enddo
       enddo
     enddo

end subroutine compute_p_z

!#######################################################################
end module scm_gcss_arm_mod
