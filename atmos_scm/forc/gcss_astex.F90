module gcss_astex_mod
 
   use            mpp_mod, only:  mpp_pe, mpp_root_pe, stdlog
   use         mpp_io_mod, only:  mpp_open, MPP_RDONLY
   use            fms_mod, only: write_version_number, open_file,        &
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

   use sat_vapor_pres_mod, only:  lookup_es,  compute_qs
   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav, pi, SECONDS_PER_DAY

   use      scm_utils_mod, only:  us_std_atm, locate

   use vert_advection_mod, only:  vert_advection, SECOND_CENTERED, &
                                  FOURTH_CENTERED, FINITE_VOLUME_LINEAR, &
                                  FINITE_VOLUME_PARABOLIC, &
                                  SECOND_CENTERED_WTS, FOURTH_CENTERED_WTS, &
                                  ADVECTIVE_FORM
   use   field_manager_mod, only: MODEL_ATMOS
   use  tracer_manager_mod, only: get_tracer_index

use       constants_mod, only: kappa
use             fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                               endlon, rlonb, rlatb,  cold_start, ncnst, &
                               pnats, consv_te, ptop, fv_init, fv_domain, &
                               fv_end, change_time, p_var, restart_format, area, &
                               ak, bk, rlon, rlat, ng_d, f_d, nt_prog, get_eta_level

   implicit none
   private

   public gcss_astex_data_read, gcss_astex_forc_init, gcss_astex_forc_end, update_gcss_astex_forc, &
          gcss_astex_forc_diagnostic_init, add_gcss_astex_tdtlw, add_gcss_astex_tdtsw,       &
          get_gcss_astex_flx

   character(len=10) :: mod_name = 'gcss_astex'
   character(len=7) :: mod_name_diag = 'forcing'

   real, public, allocatable, dimension(:,:,:)  :: u_geos, v_geos, tdt_lw, tdt_sw

   ! arrays to hold initial sounding

   integer, parameter :: ksnd = 42
   real, dimension(ksnd) :: p_snd,  T_snd,  qv_snd, ql_snd, u_snd,  v_snd,  ug_snd,  vg_snd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       p00                      reference pressure (pascals)
!
!       tracer_vert_advec_scheme Which advection scheme should be used?
!
!       temp_vert_advec_scheme   Which advection scheme should be used?
!
!                                SECOND_CENTERED         1
!                                FOURTH_CENTERED         2
!                                FINITE_VOLUME_LINEAR    3 
!                                FINITE_VOLUME_PARABOLIC 4
!                                SECOND_CENTERED_WTS     5 
!                                FOURTH_CENTERED_WTS     6
!
!       vert_advec_cond          Should condensate be vertically 
!                                advected?
!
!       divf                     large scale divergence to compute omega
!
!       p_omega_zero             pressure above which omega is set to
!                                to zero
!
!       do_rad                   include radiation forcings?
!       do_geo                   include geostrophic forcings?
!       do_vadv                  include large scale subsidence?
!
!       psfc                     surface pressure (Pa)
!       zsfc                     topography height (m)
!       ustar_sfc                surface u_star (m/s)
!       flux_tw                  surface sensible flux (K m/s)
!       flux_qw                  surface latent heat flux (m/s)

real,    private               :: p00 = 100000.
character(len=5)              :: configuration = "gcss"

integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
logical, public                :: vert_advec_cond = .true.

real,    public                :: divf = 5.e-6
real,    public                :: p_omega_zero = 85000.

logical, public                :: do_rad = .true.
logical, public                :: do_geo = .true.
logical, public                :: do_vadv = .true.

real,    private               :: psfc = 1029.e2
real,    private               :: zsfc =    0.0
real,    private               :: ustar_sfc =  0.25
real,    private               :: flux_tw   =  1.e-2
real,    private               :: flux_qw   =  1.e-5

real,    private               :: missing_value = -999.

logical, private               :: initialized = .false.
logical                        :: do_netcdf_restart = .true.

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)
        
! namelist
namelist /gcss_astex_nml/ tracer_vert_advec_scheme,                   &
                        temp_vert_advec_scheme,                     &
                        momentum_vert_advec_scheme,                 &
                        vert_advec_cond,                            &
                        p_omega_zero, divf,                         &
                        do_rad, do_geo, do_vadv,                    &
                        ustar_sfc, flux_tw, flux_qw, configuration

! diagnostics

integer ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                           &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,    id_qndt_vadv,         &
            id_udt_vadv, id_udt_geos, id_udt_lf,                &
            id_vdt_vadv, id_vdt_geos, id_vdt_lf,                &
            id_flx_radf, id_zi_forc,                            &
            id_pf_forc, id_ph_forc, id_zf_forc, id_zh_forc,     &
            id_u_geos, id_v_geos

! ---> h1g, 2010-09-27
integer ::    id_qvdt_forc_col
integer ::    id_qldt_vadv_col
! <--- h1g, 2010-09_27

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!#######################################################################
! Subroutine to read case specific namelist

subroutine gcss_astex_data_read(kmax)

implicit none

integer,  intent (in)                     :: kmax
integer                                   :: unit,ierr,io, logunit
character*23                              :: tracer_ascheme,temp_ascheme
character*64                              :: fname

character*64                              :: fname_res='INPUT/gcss_astex.res.nc'
! integer                                   :: vers  !restart file version
integer k

! local variables
character*10                              :: A_tmp 
real                                      :: theta_snd, Nd_snd
integer                                   :: p_snd_int, itype
integer                                   :: kk, k_dry

   if (initialized) return
   initialized = .true.
      
!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   READ (input_nml_file, nml=gcss_astex_nml, iostat=io)
   ierr = check_nml_error(io, 'gcss_astex_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=gcss_astex_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'gcss_astex_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=gcss_astex_nml)
   endif

!--------- allocate memory ---------

   if (allocated(u_geos)) deallocate(u_geos);  allocate(u_geos(1,1,kmax))
   if (allocated(v_geos)) deallocate(v_geos);  allocate(v_geos(1,1,kmax))
   if (allocated(tdt_lw)) deallocate(tdt_lw);  allocate(tdt_lw(1,1,kmax))
   if (allocated(tdt_sw)) deallocate(tdt_sw);  allocate(tdt_sw(1,1,kmax))

! ---> h1g, 2010-10-05
   if (file_exist('INPUT/gcss_astex.res.nc') ) then
      if(mpp_pe() == mpp_root_pe() ) call mpp_error ('gcss_astex_mod', &
         'Reading netCDF formatted restart file: INPUT/gcss_astex.res.nc', NOTE)
      if( allocated(u_geos) ) call read_data (fname_res, 'u_geos',  u_geos)
      if( allocated(v_geos) ) call read_data (fname_res, 'v_geos',  v_geos)
   endif
! <--- h1g, 2010-10-05

       
!--------- read initial sounding from file ---------
! initial condition is the hour-12 of ASTEX Lagrangian 1 (June 13.167) 


   fname='INPUT/gcss_astex_lagr1_h12' 
   call mpp_open(unit,fname,action=MPP_RDONLY)
   do k=1, 7
      read(unit,'(A10)')  A_tmp
   enddo

   k_dry = 0
   do k=1,ksnd
      read( unit,'(i5,  3(2x, f7.2),  5(1x, f7.2),  i5)' )  &
        p_snd_int, &
        theta_snd,qv_snd(k),u_snd(k), &
        v_snd(k),ug_snd(k),vg_snd(k),ql_snd(k),Nd_snd, &
        itype
        p_snd(k) = real(p_snd_int) * 100.0    ! convert from mb to pa
        qv_snd(k) = qv_snd(k) + ql_snd(k)
        ql_snd(k) = 0.0

        if ( trim(configuration)== "dry" ) then
               print*, 'k=',k,    qv_snd(k),  k_dry
           if( qv_snd(k) >= 9.0 .and. k_dry==0 ) then
            k_dry = k-1
            do  kk = 1, k_dry
                qv_snd(kk) = qv_snd(kk) * 0.44
               ! convert from temperature to potential temperature
                T_snd(kk)  =  T_snd(kk) / ( ( p_snd(kk)/p00 )**kappa )
                T_snd(kk)  =  T_snd(kk) + 3.2 * (p_snd(kk) - p_snd(1))/(p_snd(k_dry) - p_snd(1))
               ! convert  back from potential temperature to temperature
                T_snd(kk)  =  T_snd(kk) * ( ( p_snd(kk)/p00 )**kappa )
            enddo 
           endif  ! qv_snd(k+1) > 9.0 .and. k_dry==0

        endif  ! configuration
        T_snd(k) = theta_snd*( p_snd(k)/p00 )**kappa
   end do

   call close_file(unit)

end subroutine gcss_astex_data_read

!#######################################################################
! Subroutine to initialize case forcings

subroutine  gcss_astex_forc_init(time_interp,As,Bs,elev)
#include "fv_arrays.h"

!      VARIABLES
!      ------
!      INPUT:
!      ------
!      time_interp     time
!      As, Bs          A's and B's of half levels in hybrid coordinate
!                      ph(k) = A(k) + B(k) * psurf
!
!      -------
!      OUTPUT:
!      -------
!      elev            topography elevation
!      ps              surface pressure
!
!      -------------
!      INPUT/OUTPUT:
!      -------------
!      pt               temperature in degrees Kelvin 
!      u               zonal wind (m/s)
!      v               meridional wind (m/s)
!      qv              water vapor specific humidity (kg water vapor/     &
!                                                     kg air)
!      ql              liquid water condensate specific humidity 
!                      (kg condensate/ kg air)
!      qi              ice H2O condensate specific humidity 
!                      (kg condensate/kg dry air)
!      qa              cloud fraction (fraction)

type(time_type)                          :: time_interp
real,  intent (in), dimension(:)         :: As,Bs
real,  intent (out)  , dimension(:,:)    :: elev

!  Internal variables
!  ------------------------------------------------------

integer  :: kdim, k

real, dimension(size(pt,3)+1) :: eta, peta
real, dimension(size(pt,3))   :: T_us_std, qv_us_std

real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf
!  ------------------------------------------------------
integer :: i, j
#include "fv_point.inc"
       nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
       nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
       nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
       nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
       nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   !h1g
       nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )   !h1g

!   Initialize surface pressure and topography

    ps = psfc
    elev = zsfc

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
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo
    
!   Initialize pt, qv first with US standard atmosphere

    do k=1,kdim
      call us_std_atm( pf(1,1,k), T_us_std(k), qv_us_std(k) )
      pt(:,:,k)  = T_us_std(k)
      q(:,:,k,nsphum) = qv_us_std(k)
    end do
    q(:,:,:,nql) = 0.0
    ua  = 0.0
    va  = 0.0
    u_geos = 0.0
    v_geos = 0.0

    do k=1,kdim
      call gcss_astex_snd( pf(1,1,k),  ua(1,1,k),  va(1,1,k),  u_geos(1,1,k), v_geos(1,1,k), pt(1,1,k),  q(1,1,k,nsphum), q(1,1,k,nql) )
      ua(:,:,k)       =  ua(1,1,k)
      va(:,:,k)       =  va(1,1,k)
      u_geos(:,:,k)   =  u_geos(1,1,k)
      v_geos(:,:,k)   =  v_geos(1,1,k)
      pt(:,:,k)       =  pt(1,1,k)
      q(:,:,k,nsphum) =  q(1,1,k,nsphum) / ( 1.0 + q(1,1,k,nsphum) )   ! convert from mixing ratio to specific humidity
      q(:,:,k,nql)    =  q(1,1,k,nql)
    enddo

!   Initialize cloud amount
    where ( q(:,:,:,nql) > 0.0 )
      q(:,:,:,nqa) = 1.0
    elsewhere
      q(:,:,:,nqa) = 0.0
    endwhere
    if( nqn > 0) q(:,:,:,nqn) = 0.0

end subroutine  gcss_astex_forc_init

!#######################################################################

subroutine gcss_astex_forc_end ()
character*64                 :: fname_res='RESTART/gcss_astex.res.nc'

  if (.not.initialized) return
  initialized = .false.
  
! ---> h1g, 2010-10-05
     if( do_netcdf_restart ) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('gcss_astex_mod', 'Writing netCDF formatted restart file: RESTART/gcss_astex.res.nc', NOTE)
           endif
           call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)
           call write_data (fname_res,  'u_geos',  u_geos)
           call write_data (fname_res,  'v_geos',  v_geos)
      endif
! <--- h1g, 2010-10-05

  deallocate (  u_geos, v_geos, tdt_lw, tdt_sw )
 
end subroutine  gcss_astex_forc_end

!#######################################################################

subroutine  gcss_astex_forc_diagnostic_init(axes, Time)

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

id_flx_radf = register_diag_field (mod_name_diag, 'flx_radf', axes(half), Time, &
     'Vertical radiative flux', 'W/m2', missing_value = missing_value)

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

id_zi_forc = register_diag_field (mod_name_diag, 'zi_forc', axes(1:2), Time, &
     'inversion height', 'm', missing_value = missing_value)

id_u_geos = register_diag_field (mod_name_diag, 'u_geos', axes(1:3), Time, &
     'U geostrophic wind', 'm/s',  missing_value = missing_value)

id_v_geos = register_diag_field (mod_name_diag, 'v_geos', axes(1:3), Time, &
     'V geostrophic wind', 'm/s',  missing_value = missing_value)


! ---> h1g, 2010-09-27
 id_qvdt_forc_col =  register_diag_field (mod_name_diag, 'qvdt_forc_col', axes(1:2), Time, &
     'column integrated vapor forcing', 'kg/m2/s',  missing_value = missing_value)
 id_qldt_vadv_col =  register_diag_field (mod_name_diag, 'qldt_vadv_col', axes(1:2), Time, &
     'column integrated cloud water vertical advection', 'kg/m2/s',  missing_value = missing_value)
! <--- h1g, 2010-09-27

end subroutine gcss_astex_forc_diagnostic_init

!#######################################################################
! Subroutine to apply gcss_astex forcings

subroutine update_gcss_astex_forc(time_interp,time_diag,dt_int,elev)
#include "fv_arrays.h"

!      ------
!      INPUT:
!      ------
!      time_interp  time for interpolation
!      time_diag    time for diagnostics
!      dt_int       time step
!      Hgrid        horizontal grid
!      Vgrid        vertical grid
!      elev         elevation

!      ------------
!      INPUT/OUTPUT
!      ------------
!
!      pt            temperature in degrees Kelvin 
!      u            zonal wind (m/s)
!      v            meridional wind (m/s)
!      qv           water vapor specific humidity
!                   (kg water vapor/kg air)
!      ql           liquid water condensate specific humidity 
!                   (kg condensate/ kg air)
!      qi           ice H2O condensate specific humidity
!                   (kg condensate/kg dry air)
!      qa           cloud fraction (fraction)
! 
!      ------
!      OUTPUT:
!      ------
!
!      ps               surface pressure
!      AT               temperature tendency due to large-scale forcing
!                       (K/sec)
!      AU               zonal wind tendency due to large-scale forcing
!                       (m/s*s)
!      AU               meridional wind tendency due to large-scale 
!                       forcing (m/s*s)
!      AQ               water vapor tendency due to large-scale forcing 
!                       (kg vapor/kg air/sec)
!      AL               liquid water condensate tendency due to 
!                       large-scale forcing (kg condensate/ kg air/sec)
!      AI               H2O ice condensate tendency due to large-scale
!                       forcing (kg condensate/ kg air/ sec)
!      AA               cloud fraction tendency due to large-scale 
!                       forcing (fraction / sec)
!      omega_f          omega interpolated to full model levels 
!                       (Pa/sec)
!

type(time_type), intent(in)              :: time_interp,time_diag,dt_int
real,  intent (in),    dimension(:,:)    :: elev

integer                                          :: i,j,k,kdim
integer                                          :: dt_seconds,dt_days
logical                                          :: used
real                                             :: fcriolis, term3, dzi, Qz1, Qz2, dts
real, dimension(size(pt,3))                       :: klwp
real, dimension(size(pt,1),size(pt,2))             :: zi
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf,zf,th
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pi_fac, dp
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, omega_h, zh, frad
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dT_rad, dT_adi, dT_vadv, dT_lf, &
                                                    dqv_vadv, dqv_lf, dql_vadv, dqi_vadv, dqa_vadv, dqn_vadv
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: du_vadv, dv_vadv, du_geos, dv_geos, du_lf, dv_lf

! ---> h1g, 2010-09-27
real, dimension(size(pt,1),size(pt,2))  ::  qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::  qldt_vadv_col
! <--- h1g, 2010-09-27
!  ------------------------------------------------------------------------------------------------------------
#include "fv_point.inc"

! --- update pf, ph, and zf, zh
ps=psfc; ph(:,:,1)=0.;
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo
       call error_mesg ('rf02 ',   &
            'need to compute height', FATAL )

! --- compute large-scale subsidence
kdim = size(pt,3)
do k=1,kdim
   if (pf(1,1,k) > p_omega_zero) then
      omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav * zf(1,1,k)*divf
   else
      omga(:,:,k) =0.
   end if
end do

! --- compute dp, pi_fac, theta
do k = 2,kdim+1
   dp(:,:,k-1) = ph(:,:,k) - ph(:,:,k-1)
   pi_fac(:,:,k-1)= (pf(:,:,k-1)/p00)**(rdgas/cp_air)
enddo
th(:,:,:) = pt(:,:,:) / pi_fac(:,:,:)

! --- compute geostrophic tendencies
du_geos = 0.0
dv_geos = 0.0
if (do_geo) then
   fcriolis = f_d(1,1)
   do k=1, kdim
      du_geos(:,:,k) =   fcriolis * (va(1,1,k)-v_geos(1,1,k))
      dv_geos(:,:,k) = - fcriolis * (ua(1,1,k)-u_geos(1,1,k))
   end do
end if

! --- large-scale forcing
du_lf = 0.0
dv_lf = 0.0
dT_lf = 0.0
dqv_lf = 0.0

tdt_lw = 0.0
tdt_sw = 0.0
dT_rad = 0.0

if (do_rad) then
   do k=1, kdim
      klwp(k)=85.* q(1,1,k,nql) * dp(1,1,k)/grav
   end do
   frad(1,1,:)=0.;

      do i = kdim, 2, -1
         if ( ( q(1,1,i,nsphum)/( 1.-q(1,1,i,nsphum) )+q(1,1,i,nql)   )<0.01 ) then 
!---> h1g, get inversion height zi
           zi(1,1) =  zf(1,1,i)*( (  q(1,1,i+1,nsphum)/( 1.- q(1,1,i+1,nsphum) )+q(1,1,i+1,nql))-0.01) &
               +zf(1,1,i+1)*( 0.01-( q(1,1,i,nsphum)/( 1.- q(1,1,i,nsphum) ) +q(1,1,i,nql)) )

           zi(1,1) = zi(1,1) / ( ( q(1,1,i+1,nsphum)/( 1.- q(1,1,i+1,nsphum) )+q(1,1,i+1,nql)) &
                              -( q(1,1,i,nsphum)/(1.- q(1,1,i,nsphum) )+q(1,1,i,nql)) )
           print*, 'zi = ',  zi
          exit 
!<--- h1g, end zi
         endif
      end do

   do k=1, kdim
      Qz1=0.; Qz2=0.;
      do i=1, k
         Qz1 = Qz1 + klwp(i)
      end do
      do i=k+1, kdim
         Qz2 = Qz2 + klwp(i)
      end do

      dzi=zf(1,1,k)-zi(1,1);

      if ( dzi > 0. ) then
         term3= cp_air * divf * (0.25*dzi**(4./3.) + zi(1,1)*dzi**(1./3.))
         term3= term3*pf(1,1,k)/(rdgas*pt(1,1,k))
      else
         term3=0.;
      end if

      frad(1,1,k+1)=70.*exp(-Qz1) + 22.*exp(-Qz2) + term3
      ! if ( k > 31 .and. k < 80 ) print*, ' k=', k, pf(1,1,k+1)*1.e-2, frad(1,1,k+1)
   end do
   frad(1,1,1)=frad(1,1,2);
   do k=1,kdim
      dT_rad(1,1,k)=grav/cp_air*(frad(1,1,k+1)-frad(1,1,k))/dp(1,1,k)  
   end do
end if

tdt_lw = dT_rad

dT_vadv=0.; dT_adi=0.; dqv_vadv=0.; dql_vadv=0.; dqi_vadv=0.;dqa_vadv=0.;dqn_vadv=0.;
du_vadv=0.; dv_vadv=0.;

if (do_vadv) then
   call get_time(dt_int,dt_seconds,dt_days)
   dts = real(dt_seconds + 86400*dt_days)

   dT_adi=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pf(:,:,:)

   omega_h=0.0
   do k=2,kdim
      if (ph(1,1,k) > p_omega_zero) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1))) * grav * zh(1,1,k)*divf
      else
         omega_h(:,:,k) =0.
      end if
   end do
 
   select case (temp_vert_advec_scheme)
            CASE(1)
                 vadvec_scheme = SECOND_CENTERED
            CASE(2)
                 vadvec_scheme = FOURTH_CENTERED
            CASE(3)
                 vadvec_scheme = FINITE_VOLUME_LINEAR
            CASE(4)
                 vadvec_scheme = FINITE_VOLUME_PARABOLIC
            CASE(5)
                 vadvec_scheme = SECOND_CENTERED_WTS
            CASE(6)
                 vadvec_scheme = FOURTH_CENTERED_WTS
   END SELECT          
   call vert_advection(dts,omega_h,delp,pt,dT_vadv,&
                       scheme=vadvec_scheme,form=ADVECTIVE_FORM)

   SELECT CASE (tracer_vert_advec_scheme)
          CASE(1)
               vadvec_scheme = SECOND_CENTERED
          CASE(2)
               vadvec_scheme = FOURTH_CENTERED
          CASE(3)
               vadvec_scheme = FINITE_VOLUME_LINEAR
          CASE(4)
               vadvec_scheme = FINITE_VOLUME_PARABOLIC
          CASE(5)
               vadvec_scheme = SECOND_CENTERED_WTS
          CASE(6)
               vadvec_scheme = FOURTH_CENTERED_WTS
   END SELECT          
   call vert_advection(dts,omega_h,delp,q(:,:,:,nsphum),dqv_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   call vert_advection(dts,omega_h,delp,q(:,:,:,nql),dql_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   call vert_advection(dts,omega_h,delp,q(:,:,:,nqi),dqi_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   call vert_advection(dts,omega_h,delp,q(:,:,:,nqa),dqa_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   if( nqn > 0 ) &
     call vert_advection(dts,omega_h,delp,q(:,:,:,nqn),dqn_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)


   select case (momentum_vert_advec_scheme)
   case(1)
      call vert_advection(dts,omega_h,delp,ua,du_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,delp,va,dv_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
   case(2)
      call vert_advection(dts,omega_h,delp,pt,dT_vadv,scheme=FOURTH_CENTERED,form=ADVECTIVE_FORM)
   case(3)
      call vert_advection(dts,omega_h,delp,pt,dT_vadv,scheme=FINITE_VOLUME_LINEAR,form=ADVECTIVE_FORM)
   case(4)
      call vert_advection(dts,omega_h,delp,pt,dT_vadv,scheme=FINITE_VOLUME_PARABOLIC,form=ADVECTIVE_FORM)
   case(5)
      call vert_advection(dts,omega_h,delp,pt,dT_vadv,scheme=SECOND_CENTERED_WTS,form=ADVECTIVE_FORM)
   case(6)
      call vert_advection(dts,omega_h,delp,pt,dT_vadv,scheme=FOURTH_CENTERED_WTS,form=ADVECTIVE_FORM)
   end select

end if

dT_vadv= dT_vadv + dT_adi
t_dt = t_dt + dT_rad + dT_vadv
q_dt(:,:,:,nsphum) = dqv_vadv
q_dt(:,:,:,nql) = dql_vadv
q_dt(:,:,:,nqi) = dqi_vadv
q_dt(:,:,:,nqa) = dqa_vadv
if(nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv

u_dt = du_vadv + du_geos
v_dt = dv_vadv + dv_geos

if (id_pf_forc > 0)  used = send_data( id_pf_forc,  pf      (:,:,:), time_diag, 1, 1)

if (id_ph_forc > 0)  used = send_data( id_ph_forc,  ph      (:,:,:), time_diag, 1, 1)

if (id_zf_forc > 0)  used = send_data( id_zf_forc,  zf      (:,:,:), time_diag, 1, 1)

if (id_zh_forc > 0)  used = send_data( id_zh_forc,  zh      (:,:,:), time_diag, 1, 1)

if (id_tdt_vadv > 0) used = send_data( id_tdt_vadv, dT_vadv (:,:,:), time_diag, 1, 1)

if (id_udt_vadv > 0) used = send_data( id_udt_vadv, du_vadv (:,:,:), time_diag, 1, 1)

if (id_vdt_vadv > 0) used = send_data( id_vdt_vadv, dv_vadv (:,:,:), time_diag, 1, 1)

if (id_udt_geos > 0) used = send_data( id_udt_geos, du_geos (:,:,:), time_diag, 1, 1)

if (id_vdt_geos > 0) used = send_data( id_vdt_geos, dv_geos (:,:,:), time_diag, 1, 1)

if (id_qvdt_vadv> 0) used = send_data( id_qvdt_vadv,dqv_vadv(:,:,:), time_diag, 1, 1)

if (id_qldt_vadv> 0) used = send_data( id_qldt_vadv,dql_vadv(:,:,:), time_diag, 1, 1)

if (id_qidt_vadv> 0) used = send_data( id_qidt_vadv,dqi_vadv(:,:,:), time_diag, 1, 1)

if (id_qadt_vadv> 0) used = send_data( id_qadt_vadv,dqa_vadv(:,:,:), time_diag, 1, 1)


! ---> h1g, 2010-09-27
   qvdt_forcing_col = 0.0
   do k=1, kdim
         qvdt_forcing_col =  qvdt_forcing_col &
            + ( dqv_vadv( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do

   qldt_vadv_col = 0.0
   do k=1, kdim
         qldt_vadv_col =  qldt_vadv_col &
            + (  dql_vadv( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do
! <--- h1g, 2010-09-27



if( nqn > 0 ) then
  if (id_qndt_vadv> 0) used = send_data( id_qndt_vadv,dqn_vadv(:,:,:), time_diag, 1, 1)
endif

if (id_tdt_radf > 0) used = send_data( id_tdt_radf, dT_rad  (:,:,:), time_diag, 1, 1)

if (id_flx_radf > 0) used = send_data( id_flx_radf, frad    (:,:,:), time_diag, 1, 1)

if (id_udt_lf > 0)   used = send_data( id_udt_lf,   du_lf   (:,:,:), time_diag, 1, 1)

if (id_vdt_lf > 0)   used = send_data( id_vdt_lf,   dv_lf   (:,:,:), time_diag, 1, 1)

if (id_tdt_lf > 0)   used = send_data( id_tdt_lf,   dT_lf   (:,:,:), time_diag, 1, 1)

if (id_qvdt_lf > 0)  used = send_data( id_qvdt_lf,  dqv_lf  (:,:,:), time_diag, 1, 1)

if (id_zi_forc > 0)  used = send_data( id_zi_forc,  zi      (:,:),   time_diag, 1, 1)

if (id_u_geos > 0)   used = send_data( id_u_geos,   u_geos  (:,:,:), time_diag, 1, 1)

if (id_v_geos > 0)   used = send_data( id_v_geos,   v_geos  (:,:,:), time_diag, 1, 1)

! ---> h1g, 2010-09-27
if ( id_qvdt_forc_col > 0 )  used = send_data(  id_qvdt_forc_col, qvdt_forcing_col(:,:), time_diag, 1, 1 )
if ( id_qldt_vadv_col > 0 )  used = send_data(  id_qldt_vadv_col, qldt_vadv_col(:,:), time_diag, 1, 1 )
! <--- h1g, 2010-09-27

end subroutine update_gcss_astex_forc

!########################################################################
! This subroutine adds longwave radiative heating from forcings
! to input arrays

subroutine add_gcss_astex_tdtlw( x )

  implicit none
  real, intent(inout) :: x(:,:,:)

  if (allocated(tdt_lw)) x = x + tdt_lw

end subroutine add_gcss_astex_tdtlw

!########################################################################
! This subroutine adds shortwave radiative heating from forcings
! to input arrays

subroutine add_gcss_astex_tdtsw( x )

  implicit none
  real, intent(inout) :: x(:,:,:)

  if (allocated(tdt_sw)) x = x + tdt_sw

end subroutine add_gcss_astex_tdtsw

!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_gcss_astex_flx( ustar, flux_t, flux_q )

  implicit none
  real, intent(out), dimension(:) :: ustar, flux_t, flux_q

  ustar = ustar_sfc
  flux_t = flux_tw
  flux_q = flux_qw

end subroutine get_gcss_astex_flx
 !########################################################################




!########################################################################
! Subroutine retuns gcss_astex initial sounding
       subroutine gcss_astex_snd( p, u, v, ug, vg, T, qv, ql )
       implicit none

       real, intent(in) :: p
       real, intent(out) :: u, v, ug, vg, T, qv, ql

       integer k
       real a, b

!      Interpolate to get T, qv, ql

       call locate( p_snd, ksnd, p, k )
       if ( k < 1 .or. k >= ksnd ) then
          goto 120
     !    call error_mesg('gcss_astex_snd',  &
     !                    'p value out of sounding bounds', FATAL)
       endif

       b = ( p - p_snd(k) ) / ( p_snd(k+1) - p_snd(k) )
       a = 1.0 - b
       u  = a * u_snd(k)  + b * u_snd(k+1)
       v  = a * v_snd(k)  + b * v_snd(k+1)
       ug = a * ug_snd(k) + b * ug_snd(k+1)
       vg = a * vg_snd(k) + b * vg_snd(k+1)
       T  = a * T_snd(k)  + b * T_snd(k+1)
       qv = ( a * qv_snd(k) + b * qv_snd(k+1) ) * 0.001
       qv = qv / (1. + qv )   ! convert from mixing ratio to specific humidity
       ql = ( a * ql_snd(k) + b * ql_snd(k+1) ) * 0.001

120    continue
       return
       end subroutine gcss_astex_snd
!########################################################################

end module gcss_astex_mod
 
