module scm_rf02_mod
 
   use            mpp_mod, only:  stdlog
   use            fms_mod, only:  check_nml_error,  &
                                  mpp_pe, mpp_root_pe,                  &
                                  write_version_number,                 &
                                  open_file, close_file, file_exist,    &
                                  read_data, write_data, mpp_error,     &
                                  error_mesg, FATAL, NOTE
#ifdef INTERNAL_FILE_NML
   USE              mpp_mod, ONLY: input_nml_file
#else
   USE              fms_mod, ONLY: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod

   use sat_vapor_pres_mod, only:  lookup_es, sat_vapor_pres_init   ! ZNT 02/20/2020
   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav, pi, seconds_per_day

   use      scm_utils_mod, only:  us_std_atm, thetal_to_temp

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

   public rf02_data_read, rf02_forc_init, rf02_forc_end, update_rf02_forc, &
          rf02_forc_diagnostic_init, add_rf02_tdtlw, add_rf02_tdtsw,       &
          get_rf02_flx

   character(len=8) :: mod_name = 'scm_rf02'
   character(len=7) :: mod_name_diag = 'forcing'

   real, public, allocatable, dimension(:,:,:)  :: u_geos, v_geos, tdt_lw, tdt_sw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       p00                      reference pressure (pascals)
!
!       configuration            case configuration (used only if multiple
!                                configurations share the same module)
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
!       rho_sfc                  surface layer air density for flux computations (kg/m3)
!       ustar_sfc                surface u_star (m/s)
!       flux_t_sfc               surface sensible flux (W/m2)
!       flux_q_sfc               surface latent heat flux (W/m2)

real,    private               :: p00 = 100000.

character(len=64)              :: configuration = 'base'

integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
logical, public                :: vert_advec_cond = .true.

real,    public                :: divf = 3.75e-6
real,    public                :: p_omega_zero = 70000.

logical, public                :: do_rad = .true.
logical, public                :: do_geo = .true.
logical, public                :: do_vadv = .true.

real,    private               :: psfc = 1017.8e2
real,    private               :: zsfc =    0.0
real,    private               :: rho_sfc     = 1.21
real,    private               :: ustar_sfc   = 0.25
real,    private               :: flux_t_sfc  = 16.0
real,    private               :: flux_q_sfc  = 93.0

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

namelist /scm_rf02_nml/ tracer_vert_advec_scheme,                   &
                        temp_vert_advec_scheme,                     &
                        momentum_vert_advec_scheme,                 &
                        vert_advec_cond,                            &
                        p_omega_zero, divf,                         &
                        do_rad, do_geo, do_vadv,                    &
                        rho_sfc, ustar_sfc, flux_t_sfc, flux_q_sfc, &
                        configuration

! diagnostics

integer ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                           &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,   id_qndt_vadv,        &
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

subroutine rf02_data_read(kmax)

implicit none

integer,  intent (in)                     :: kmax
integer                                   :: unit,ierr,io, logunit
character*23                              :: tracer_ascheme,temp_ascheme

character*64                              :: fname_res='INPUT/rf02.res.nc'


   if (initialized) return
   initialized = .true.
      
!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   READ (input_nml_file, nml=scm_rf02_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_rf02_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=scm_rf02_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_rf02_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=scm_rf02_nml)
   endif

!--------- allocate memory ---------

   if (allocated(u_geos)) deallocate(u_geos);  allocate(u_geos(1,1,kmax))
   if (allocated(v_geos)) deallocate(v_geos);  allocate(v_geos(1,1,kmax))
   if (allocated(tdt_lw)) deallocate(tdt_lw);  allocate(tdt_lw(1,1,kmax))
   if (allocated(tdt_sw)) deallocate(tdt_sw);  allocate(tdt_sw(1,1,kmax))

! ---> h1g, 2010-10-07
   if (file_exist('INPUT/rf02.res.nc') ) then
      if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_rf02_mod', &
         'Reading netCDF formatted restart file: INPUT/rf02.res.nc', NOTE)
      if( allocated(u_geos) ) call read_data (fname_res, 'u_geos',  u_geos)
      if( allocated(v_geos) ) call read_data (fname_res, 'v_geos',  v_geos)
   endif
! <--- h1g, 2010-10-07
       
end subroutine rf02_data_read

!#######################################################################
! Subroutine to initialize case forcings

subroutine rf02_forc_init(time_interp, As, Bs)
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

!  Internal variables
!  ------------------

integer  :: kdim, k

integer, parameter :: itmax = 10

real, dimension(size(pt,3)+1) :: eta, peta
real, dimension(size(pt,3))   :: T_us_std, qv_us_std
real, dimension(size(pt,3))   :: thetal_rf02, rt_rf02, u_rf02, v_rf02, T_rf02, qv_rf02, ql_rf02 

real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, zh, zhnew
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf, zf, zfnew

! ZNT 02/20/2020: Test
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph0
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf0


real maxerror

integer i
real,  dimension(size(pt,1),size(pt,2))    :: elev
integer :: j
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
    !        call error_mesg ('rf02 ',   &
    !        'need to compute height', FATAL )

!   Iteration to compute blended RF02 sounding

    i=0
    maxerror = 1.0
    do while ( maxerror > 0.001 .and. i < itmax )

     i = i + 1
     do k=1,kdim

      call rf02_snd( zf(1,1,k), u_rf02(k), v_rf02(k), thetal_rf02(k), rt_rf02(k) )
      call thetal_to_temp( thetal_rf02(k), rt_rf02(k), pf(1,1,k),  &
                           T_rf02(k), qv_rf02(k), ql_rf02(k) )
      pt(:,:,k)  = max ( T_rf02(k), 200.0 )
      ua(:,:,k) = u_rf02(k)
      va(:,:,k) = v_rf02(k)
      q(:,:,k,nsphum) = qv_rf02(k)
      q(:,:,k,nql) = ql_rf02(k)
     enddo
     
     print*,'=======interpolate initial sounding ========='
     print*,'rt_rf02(kdim)', rt_rf02(kdim), q(:,:,kdim,nsphum),  q(:,:,kdim,nql)
     print*,'================================'

     call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zfnew, zhnew)

     ! call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zfnew, zhnew )
     !       call error_mesg ('rf02 ',   &
     !       'need to compute height', FATAL )

     maxerror = maxval( abs( zfnew - zf ) )
     zf = zfnew
     zh = zhnew

    enddo

    if ( i >= itmax ) then
      call error_mesg('rf02_forc_init',  &
                      'failed to converge while creating blended sounding', FATAL)
    endif

!   Initialize cloud amount

    where ( q(:,:,:,nql) > 0.0 )
      q(:,:,:,nqa) = 1.0
    elsewhere
      q(:,:,:,nqa) = 0.0
    endwhere
    if( nqn > 0) q(:,:,:,nqn) = 0.0

!   Initialize geostrophic winds

    u_geos = ua
    v_geos = va

    ! ZNT 02/23/2020 - Note: Also need to initialize u_srf and v_srf
    u_srf(:,:)=ua(:,:,KDIM)
    v_srf(:,:)=va(:,:,KDIM)

end subroutine rf02_forc_init

!#######################################################################

subroutine rf02_forc_end ()
character*64                 :: fname_res='RESTART/rf02.res.nc'

  if (.not.initialized) return
  initialized = .false.

! ---> h1g, 2010-10-07
     if( do_netcdf_restart ) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('scm_rf02_mod', 'Writing netCDF formatted restart file: RESTART/rf02.res.nc', NOTE)
           endif
           call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)	   
           call write_data (fname_res,  'u_geos',  u_geos)
           call write_data (fname_res,  'v_geos',  v_geos)
      endif
! <--- h1g, 2010-10-07

  deallocate (  u_geos, v_geos, tdt_lw, tdt_sw )
 
end subroutine rf02_forc_end

!#######################################################################

subroutine rf02_forc_diagnostic_init(axes, Time)

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

end subroutine rf02_forc_diagnostic_init

!#######################################################################
! Subroutine to apply RF02 forcings

subroutine update_rf02_forc(time_interp,time_diag,dt_int)
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
real,  dimension(size(pt,1),size(pt,2))    :: elev  ! znt 20200226

! ---> h1g, 2010-09-27
real, dimension(size(pt,1),size(pt,2))  ::  qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::  qldt_vadv_col
! <--- h1g, 2010-09-27
#include "fv_point.inc"

! --- update pf, ph, and zf, zh
elev=zsfc   ! znt 20200226
ps=psfc; ph(:,:,1)=0.;
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo

       call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zf, zh)

       ! call compute_height    (Vgrid, grav*elev, T, qv, pf, ph, zf, zh)
       ! call error_mesg ('rf02 ',   &
       !      'need to compute height', FATAL )


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
   fcriolis=f_d(1,1)
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

! --- longwave radiative forcing
!     F(z) = F0*exp(-Q(z,inf))
!          + F1*exp(-Q(0,z))
!          + alpha*rho_i*cp*Divf*H(z-zi)*(0.25(z-zi)4/3+zi(z-zi)1/3)
dT_rad=0.
if (do_rad) then
   do k=1, kdim
      klwp(k)=85.* q(1,1,k,nql) * dp(1,1,k)/grav
   end do
   frad(1,1,:)=0.;
   do k=1, kdim
      Qz1=0.; Qz2=0.;
      do i=1, k
         Qz1 = Qz1 + klwp(i)
      end do
      do i=k+1, kdim
         Qz2 = Qz2 + klwp(i)
      end do

      do i = kdim, 2, -1
         if ( (q(1,1,i,nsphum)+q(1,1,i,nql))<0.008 ) then 
	!h1g
          !  zi(1,1) = zf(1,1,i) ! original cjg
           zi(1,1) =  zf(1,1,i)*( (q(1,1,i+1,nsphum)+q(1,1,i+1,nql))-0.008) &
	               +zf(1,1,i+1)*( 0.008-(q(1,1,i,nsphum)+q(1,1,i,nql)) )
           		       
           zi(1,1) = zi(1,1) / ( (q(1,1,i+1,nsphum)+q(1,1,i+1,nql)) &
	                              -(q(1,1,i,nsphum)+q(1,1,i,nql)) )
          exit 
	!h1g
         endif
      end do
      dzi=zf(1,1,k)-zi(1,1);
      !h1g
      ! if ( (dzi > 0.) .and. (ph(1,1,k+1) > p_omega_zero)) then
      !    term3=1. * 1.12 * cp_air * divf * (0.25*dzi**(4./3.) + zi(1,1)*dzi**(1./3.))

      if (  dzi > 0. ) then
         term3= cp_air * divf * (0.25*dzi**(4./3.) + zi(1,1)*dzi**(1./3.))
	 term3= term3*pf(1,1,k)/(rdgas*pt(1,1,k))
      elseif ( dzi == 0. ) then
         term3= cp_air * divf * (0.25*dzi**(4./3.) + zi(1,1)*dzi**(1./3.))
	 term3= 0.5* term3*pf(1,1,k)/(rdgas*pt(1,1,k))
      else
         term3=0.;
      end if
      !h1g

      frad(1,1,k+1)=70.*exp(-Qz1) + 22.*exp(-Qz2) + term3
   end do
   frad(1,1,1)=frad(1,1,2);
   do k=1,kdim
      dT_rad(1,1,k)=grav/cp_air*(frad(1,1,k+1)-frad(1,1,k))/dp(1,1,k)
   !    dT_rad(1,1,k)= 0.0
   end do
end if
tdt_lw = dT_rad
tdt_sw = 0.0

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
   case(1)
      call vert_advection(dts,omega_h,delp,pt,dT_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
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
t_dt = t_dt +  dT_rad + dT_vadv
q_dt(:,:,:,nsphum) = dqv_vadv
q_dt(:,:,:,nql) = dql_vadv
q_dt(:,:,:,nqi) = dqi_vadv
q_dt(:,:,:,nqa) = dqa_vadv
if(nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv

u_dt = du_vadv+du_geos; v_dt = dv_vadv+dv_geos;

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


if ( nqn > 0 ) then
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

end subroutine update_rf02_forc

!########################################################################
! This subroutine adds longwave radiative heating from forcings
! to input arrays

subroutine add_rf02_tdtlw( x )

  implicit none
  real, intent(inout) :: x(:,:,:)

  if (allocated(tdt_lw)) x = x + tdt_lw

end subroutine add_rf02_tdtlw

!########################################################################
! This subroutine adds shortwave radiative heating from forcings
! to input arrays

subroutine add_rf02_tdtsw( x )

  implicit none
  real, intent(inout) :: x(:,:,:)

  if (allocated(tdt_sw)) x = x + tdt_sw

end subroutine add_rf02_tdtsw

!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_rf02_flx( rho, ustar, flux_t, flux_q )

  implicit none
  real, intent(out), dimension(:) :: rho, ustar, flux_t, flux_q

  rho = rho_sfc
  ustar = ustar_sfc
  flux_t = flux_t_sfc
  flux_q = flux_q_sfc

end subroutine get_rf02_flx

!########################################################################
! Subroutine retuns RF02 initial sounding

       subroutine rf02_snd( z, u, v, thetal, qt )
       implicit none

       real, intent(in) :: z
       real, intent(out) :: u, v, thetal, qt

       real, parameter :: zi = 795.0
       real, parameter :: ztop = 2500.0

       if ( configuration.ne.'dry' ) then

!        Original sounding
!        Specified sounding below ztop

         if ( z <= ztop ) then
           u =  3.0 + 0.0043 * z
           v = -9.0 + 0.0056 * z
           if ( z < zi ) then
             thetal = 288.3
             qt = 9.45e-3
           else
             thetal = 295.0 + (z-zi)**(1.0/3.0)
             qt = 5.0e-3 - 3.0e-3*( 1.0 - exp( (zi-z)/500. ) )
           endif

!        Above ztop, we simply extend the sounding using constant
!        dtheta/dz, drt/dz, u, v.

         else
           u =  3.0 + 0.0043 * ztop
           v = -9.0 + 0.0056 * ztop
           thetal = 295.0 + (ztop-zi)**(1.0/3.0)                &
                    + 3.5e-3*(z-ztop)
           qt = max( 5.0e-3 - 3.0e-3*( 1.0 - exp( (zi-ztop)/500. ) ) &
                     - 2.0e-6*(z-ztop)                               &
                    ,0.0 )
         endif

       else

!        Alternate sounding modified with drier air above the inversion
!        Specified sounding below ztop

         if ( z <= ztop ) then
           u =  3.0 + 0.0043 * z
           v = -9.0 + 0.0056 * z
           if ( z < zi ) then
             thetal = 288.3
             qt = 9.45e-3
           else
             thetal = 295.0 + (z-zi)**(1.0/3.0)
             qt = 4.0e-3 - 3.0e-3*( 1.0 - exp( (zi-z)/500. ) )
           endif

!        Above ztop, we simply extend the sounding using constant
!        dtheta/dz, drt/dz, u, v.

         else
           u =  3.0 + 0.0043 * ztop
           v = -9.0 + 0.0056 * ztop
           thetal = 295.0 + (ztop-zi)**(1.0/3.0)                &
                    + 3.5e-3*(z-ztop)
           qt = max( 4.0e-3 - 3.0e-3*( 1.0 - exp( (zi-ztop)/500. ) ) &
                     - 2.0e-6*(z-ztop)                               &
                    ,0.0 )
         endif

       endif

       return
       end subroutine rf02_snd

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


end module scm_rf02_mod
 
