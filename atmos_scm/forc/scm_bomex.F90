module scm_bomex_mod
 
   use            mpp_mod, only:  stdlog
   use            fms_mod, only:  check_nml_error,                      &
                                  mpp_pe, mpp_root_pe,                  &
                                  write_version_number,                 &
                                  open_file, close_file, file_exist,    &
                                  read_data, write_data,                &
                                  mpp_error,                            &
                                  error_mesg, FATAL, NOTE
#ifdef INTERNAL_FILE_NML
use              mpp_mod, only: input_nml_file
#else
use              fms_mod, only: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod

   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav, pi

   use      scm_utils_mod, only:  us_std_atm, thetal_to_temp

   use vert_advection_mod, only:  vert_advection, SECOND_CENTERED, &
                                  FOURTH_CENTERED, FINITE_VOLUME_LINEAR, &
                                  FINITE_VOLUME_PARABOLIC, &
                                  SECOND_CENTERED_WTS, FOURTH_CENTERED_WTS, &
                                  ADVECTIVE_FORM
   use   field_manager_mod, only: MODEL_ATMOS
   use  tracer_manager_mod, only: get_tracer_index

! use       constants_mod, only: kappa
use             fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                               endlon, rlonb, rlatb,  cold_start, ncnst, &
                               pnats, consv_te, ptop, fv_init, fv_domain, &
                               fv_end, change_time, p_var, restart_format, area, &
                               ak, bk, rlon, rlat, ng_d, f_d, nt_prog, get_eta_level

use  sat_vapor_pres_mod,only : sat_vapor_pres_init   ! ZNT 02/20/2020
                                  
   implicit none
   private

   public bomex_data_read, bomex_forc_init, bomex_forc_end, update_bomex_forc, &
          bomex_forc_diagnostic_init, get_bomex_flx

character(len=9) :: mod_name = 'scm_bomex'
character(len=7) :: mod_name_diag = 'forcing'

real, public, allocatable, dimension(:,:,:)  :: u_geos, v_geos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       p00                      reference pressure (pascals)
!
!       configuration            case configuration (unused)
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
!       do_rad                   include radiation forcings?
!       do_geo                   include geostrophic forcings?
!       do_vadv                  include large scale subsidence?
!       do_nudge                 nudge u,v winds to their geostrophic value
!                                above a certain level?
!
!       psfc                     surface pressure (Pa)
!       zsfc                     topography height (m)
!       ustar_sfc                surface u_star (m/s)
!       flux_t_sfc               surface sensible flux (W/m2)
!       flux_q_sfc               surface latent heat flux (W/m2)

real,    private               :: p00 = 100000.

character(len=64)              :: configuration = 'unused'

integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
logical, public                :: vert_advec_cond = .true.

logical, public                :: do_rad = .true.
logical, public                :: do_geo = .true.
logical, public                :: do_vadv = .true.
logical, public                :: do_nudge = .true.

real,    private               :: psfc = 1015.0e2
real,    private               :: zsfc =    0.0
real,    private               :: ustar_sfc   =   0.28
! ZNT 05/15/2020: Minor difference from Siebesma paper 
! in which the non-massweighted values are prescribed (8e-3 K*m/s, 5.2e-5 kg/kg*m/s)
real,    private               :: flux_t_sfc  =   9.52
real,    private               :: flux_q_sfc  = 153.4
! ZNT 05/18/2020: Original Siebesma value
real,    private               :: wpthlp_sfc = 8.e-3
real,    private               :: wpqtp_sfc  = 5.2e-5
real,    private               :: fcor = 0.376e-4
! ZNT 06/14/2021: Stability for profile above 3000m (need to be non-convecting)
real,    private               :: dthldz_up = 3.65E-3

real,    private               :: missing_value = -999.

logical, private               :: initialized = .false.
logical                        :: do_netcdf_restart = .true.
logical                        :: do_siebesma = .true.   ! ZNT 05/18/2020

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)
        
! namelist

namelist /scm_bomex_nml/ tracer_vert_advec_scheme,                   &
                         temp_vert_advec_scheme,                     &
                         momentum_vert_advec_scheme,                 &
                         vert_advec_cond,                            &
                         do_rad, do_geo, do_vadv,                    &
                         do_nudge,                                   &
                         ustar_sfc, flux_t_sfc, flux_q_sfc,          &
                         configuration,                              &
                         wpthlp_sfc, wpqtp_sfc, fcor,  & ! ZNT 05/18/2020
                         do_siebesma, &                  ! ZNT 05/18/2020
                         dthldz_up                       ! ZNT 06/14/2021

! diagnostics

integer ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                           &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,     id_qndt_vadv,     & !h1g
            id_udt_vadv, id_udt_geos, id_udt_lf, id_udt_nudge,  &
            id_vdt_vadv, id_vdt_geos, id_vdt_lf, id_vdt_nudge,  &
            id_pf_forc, id_ph_forc, id_zf_forc, id_zh_forc,     &
            id_u_geos, id_v_geos

! ---> h1g, 2010-09-29
integer ::    id_qvdt_forc_col
integer ::    id_qldt_vadv_col
! <--- h1g, 2010-09-29

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!#######################################################################
! Subroutine to read case specific namelist

subroutine bomex_data_read(kmax)

implicit none

integer,  intent (in)                     :: kmax
integer                                   :: i
integer                                   :: unit,ierr,io, logunit
character*23                              :: tracer_ascheme,temp_ascheme

character*64                              :: fname_res='INPUT/bomex.res.nc'
     
integer :: year, month, day

   if (initialized) return
   initialized = .true.
      
!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=scm_bomex_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_bomex_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=scm_bomex_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_bomex_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=scm_bomex_nml)
   endif

!--------- allocate memory ---------

   if (allocated(u_geos)) deallocate(u_geos);  allocate(u_geos(1,1,kmax))
   if (allocated(v_geos)) deallocate(v_geos);  allocate(v_geos(1,1,kmax))

! ---> h1g, 2010-10-07
   if (file_exist('INPUT/bomex.res.nc') ) then
      if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_bomex_mod', &
         'Reading netCDF formatted restart file: INPUT/bomex.res.nc', NOTE)
      if( allocated(u_geos) ) call read_data (fname_res, 'u_geos',  u_geos)
      if( allocated(v_geos) ) call read_data (fname_res, 'v_geos',  v_geos)
   endif
! <--- h1g, 2010-10-07

end subroutine bomex_data_read

!#######################################################################
! Subroutine to initialize case forcings

subroutine bomex_forc_init(time_interp,As,Bs,elev)
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
!      pt            temperature in degrees Kelvin 
!      u            zonal wind (m/s)
!      v            meridional wind (m/s)
!      qv           water vapor specific humidity (kg water vapor/     &
!                                                  kg air)
!      ql           liquid water condensate specific humidity 
!                   (kg condensate/ kg air)
!      qi           ice H2O condensate specific humidity 
!                   (kg condensate/kg dry air)
!      qa           cloud fraction (fraction)
!      qn          cloud droplet number concentration  !h1g

type(time_type)                          :: time_interp
real,  intent (in), dimension(:)         :: As,Bs
real,  intent (out)  , dimension(:,:)    :: elev

!  Internal variables
!  ------------------

integer  :: kdim, k

integer, parameter :: itmax = 10

real, dimension(size(pt,3)+1) :: eta, peta
real, dimension(size(pt,3))   :: T_us_std, qv_us_std
real :: u_snd, v_snd, ug_snd, vg_snd, thetal_snd, qt_snd, T_snd, qv_snd, ql_snd

real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, zh, zhnew
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf, zf, zfnew

! ZNT 02/20/2020: Test
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph0
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf0


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
    !        call error_mesg ('rf01 ',   &
    !        'need to compute height', FATAL )

!   Iteration to compute BOMEX sounding

    i=0
    maxerror = 1.0
    do while ( maxerror > 0.001 .and. i < itmax )

     i = i + 1
     do k=1,kdim

      call bomex_snd( zf(1,1,k), u_snd, v_snd, ug_snd, vg_snd, thetal_snd, qt_snd )
      u_geos(:,:,k) = ug_snd
      v_geos(:,:,k) = vg_snd
      ! ZNT 05/15/2020: Note the following subroutine takes mixing ratio as input but outputs sphum.
      call thetal_to_temp( thetal_snd, qt_snd, pf(1,1,k), T_snd, qv_snd, ql_snd )
      pt(:,:,k)  = max ( T_snd, 200.0 )
      ua(:,:,k) = u_snd
      va(:,:,k) = v_snd
      q(:,:,k,nsphum) = qv_snd
      q(:,:,k,nql) = ql_snd

     enddo
     call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zfnew, zhnew)

     ! call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zfnew, zhnew )
     !       call error_mesg ('bomex ',   &
     !       'need to compute height', FATAL )
     maxerror = maxval( abs( zfnew - zf ) )
     zf = zfnew
     zh = zhnew

    enddo

    if ( i >= itmax ) then
      call error_mesg('bomex_forc_init',  &
                      'failed to converge while creating blended sounding', FATAL)
    endif

!   Initialize cloud amount

    where ( q(:,:,:,nql) > 0.0 )
      q(:,:,:,nqa) = 1.0
    elsewhere
      q(:,:,:,nqa) = 0.0
    endwhere
    if( nqn > 0) q(:,:,:,nqn) = 0.0

    ! ZNT 02/23/2020 - Note: Also need to initialize u_srf and v_srf
    u_srf(:,:)=ua(:,:,KDIM)
    v_srf(:,:)=va(:,:,KDIM)

end subroutine bomex_forc_init

!#######################################################################

subroutine bomex_forc_end ()
character*64                 :: fname_res='RESTART/bomex.res.nc'

  if (.not.initialized) return
  initialized = .false.

! ---> h1g, 2010-10-07
     if( do_netcdf_restart ) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('scm_bomex_mod', 'Writing netCDF formatted restart file: RESTART/bomex.res.nc', NOTE)
           endif
           call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)	   
           call write_data (fname_res,  'u_geos',  u_geos)
           call write_data (fname_res,  'v_geos',  v_geos)
      endif
! <--- h1g, 2010-10-07

  deallocate ( u_geos, v_geos )
 
end subroutine bomex_forc_end

!#######################################################################

subroutine bomex_forc_diagnostic_init(axes, Time)

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

! ---> h1g, 2010-09-29
 id_qvdt_forc_col =  register_diag_field (mod_name_diag, 'qvdt_forc_col', axes(1:2), Time, &
     'column integrated vapor forcing', 'kg/m2/s',  missing_value = missing_value)
 id_qldt_vadv_col =  register_diag_field (mod_name_diag, 'qldt_vadv_col', axes(1:2), Time, &
     'column integrated cloud water vertical advection', 'kg/m2/s',  missing_value = missing_value)
! <--- h1g, 2010-09-29

end subroutine bomex_forc_diagnostic_init

!#######################################################################
! Subroutine to apply BOMEX forcings

subroutine update_bomex_forc(time_interp,time_diag,dt_int,elev)
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
!
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
!      qn          droplet number concentration (#/cm3)
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
!      omega_f           omega interpolated to full model levels 
!                       (Pa/sec)
!

type(time_type), intent(in)              :: time_interp,time_diag,dt_int
real,  intent (in),    dimension(:,:)    :: elev

integer                                          :: i,j,k,kdim
integer                                          :: dt_seconds,dt_days
logical                                          :: used
real                                             :: fcriolis, fnudge, dts
real, dimension(size(pt,1),size(pt,2))             :: zi
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf,zf,th
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pi_fac, dp
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, omega_h, zh, frad
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dT_adi, dT_vadv, dT_lf, dT_rad, &
                                                    dqv_vadv, dqv_lf, dql_vadv, dqi_vadv, dqa_vadv, dqn_vadv !h1g
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: du_vadv, dv_vadv, du_geos, dv_geos, du_lf, dv_lf, du_nudge, dv_nudge

! ---> h1g, 2010-09-29
real, dimension(size(pt,1),size(pt,2))  ::  qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::  qldt_vadv_col
! <--- h1g, 2010-09-29
!  ------------------------------------------------------------------------------------------------------------
#include "fv_point.inc"


! --- update pf, ph, and zf, zh
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
   if ( zf(1,1,k) < 1500.0 ) then
      omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav  &
                       * (0.0065/1500.0)*zf(1,1,k)
   elseif ( zf(1,1,k) < 2100.0 ) then
      omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav  &
                       * ( 0.0065 - (0.0065/(2100.-1500.))*(zf(1,1,k)-1500.) )
   else
      omga(:,:,k) = 0.
   end if
end do

omega_h = 0.0
do k=2,kdim
   if ( zh(1,1,k) < 1500.0 ) then
      omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                       * (0.0065/1500.0)*zh(1,1,k)                           ! ZNT corrected 05/18/2020
   elseif ( zh(1,1,k) < 2100.0 ) then
      omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                       * (0.0065 - (0.0065/(2100.-1500.))*(zh(1,1,k)-1500.)) ! ZNT corrected 05/18/2020
   else
      omega_h(:,:,k) = 0.
   end if
end do

! --- compute dp, pi_fac, theta
do k = 2,KDIM+1
   dp(:,:,k-1) = ph(:,:,k) - ph(:,:,k-1)
   pi_fac(:,:,k-1)= (pf(:,:,k-1)/p00)**(rdgas/cp_air)
enddo
th(:,:,:) = pt(:,:,:) / pi_fac(:,:,:)

! --- radiative cooling
dT_rad=0.;
if (do_rad) then
   do k=1, kdim
      if ( zf(1,1,k) < 1500.0 ) then
         dT_rad(1,1,k) = -2.315e-5 * pi_fac(1,1,k)
      else if ( zf(1,1,k) < 2500.0 ) then
         dT_rad(1,1,k) = ( -2.315e-5 + 2.315e-5/(2500.0-1500.0)*(zf(1,1,k)-1500.0) ) &
                         * pi_fac(1,1,k)
      else
         dT_rad(1,1,k) = 0.0
      end if
   end do
end if

! --- large-scale forcing tendencies
du_lf = 0.0
dv_lf = 0.0
dT_lf = 0.0

dqv_lf=0.0
do k=1, kdim
   if ( zf(1,1,k) < 300.0 ) then
      dqv_lf(1,1,k) = -1.2e-8
   else if ( zf(1,1,k) < 500.0 ) then
      dqv_lf(1,1,k) = -( 1.2e-8 - 1.2e-8*(zf(1,1,k)-300.0)/(500.0-300.0) )
   else
      dqv_lf(1,1,k) = 0.0
   end if
end do

! --- large-scale subsidence tendencies
dT_vadv=0.0; dT_adi=0.0; dqv_vadv=0.0; dql_vadv=0.0; dqi_vadv=0.0;dqa_vadv=0.0;dqn_vadv=0.0;!h1g
du_vadv=0.0; dv_vadv=0.0;

if (do_vadv) then
   call get_time(dt_int,dt_seconds,dt_days)
   dts = real(dt_seconds + 86400*dt_days)

   dT_adi=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pf(:,:,:)

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

   dT_vadv= dT_vadv + dT_adi

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
      call vert_advection(dts,omega_h,dp,q(:,:,:,nqn),dqn_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)


   select case (momentum_vert_advec_scheme)
   case(1)
      call vert_advection(dts,omega_h,delp,ua,dU_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,delp,va,dV_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
   end select

   end if

! --- nudge winds
du_nudge = 0.0
dv_nudge = 0.0

! --- compute geostrophic tendencies
du_geos=0.; dv_geos=0.;
if (do_geo) then
   ! ZNT 05/15/2020: minor difference from Siebesma (2003) paper (computed: 0.3775e-4; paper: 0.376e-4)
   if (do_siebesma) then
      fcriolis = fcor
   else
      fcriolis = f_d(1,1)
   endif   
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
if(nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv
 
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

! ---> h1g, 2010-09-29
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
! <--- h1g, 2010-09-29

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

! ---> h1g, 2010-09-29
if ( id_qvdt_forc_col > 0 )  used = send_data(  id_qvdt_forc_col, qvdt_forcing_col(:,:), time_diag, 1, 1 )
if ( id_qldt_vadv_col > 0 )  used = send_data(  id_qldt_vadv_col, qldt_vadv_col(:,:), time_diag, 1, 1 )
! <--- h1g, 2010-09-29

end subroutine update_bomex_forc

!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_bomex_flx( rho, ustar, flux_t, flux_q )

  implicit none
  real, intent(in),  dimension(:) :: rho
  real, intent(out), dimension(:) :: ustar, flux_t, flux_q

  ustar = ustar_sfc
  if (do_siebesma) then
      flux_t = rho*cp_air*wpthlp_sfc*(psfc/p00)**(rdgas/cp_air)
      flux_q = rho*hlv*wpqtp_sfc
  else
      flux_t = flux_t_sfc
      flux_q = flux_q_sfc
  endif

end subroutine get_bomex_flx

!########################################################################
! Subroutine retuns BOMEX initial sounding

       subroutine bomex_snd( z, u, v, ug, vg, thetal, qt )
       implicit none

       real, intent(in) :: z
       real, intent(out) :: u, v, ug, vg, thetal, qt

!      u wind

       if ( z < 700.0 ) then
         u = -8.75
       else if ( z < 5500.0 ) then
         u = -8.75 + 1.8e-3 * (z - 700.)
       else
         u = 0.0
       end if

!      v wind

       v = 0.0

!      qt (specific total water content)

       if ( z < 520.0 ) then
         qt = 17.0 + (16.3 - 17.0) / (520.) * z
       else if ( z < 1480.0 ) then
         qt = 16.3 + (10.7 - 16.3) / (1480. - 520.) * (z - 520.)
       else if ( z < 2000.0 ) then
         qt = 10.7 + (4.2 - 10.7) / (2000. - 1480.) * (z - 1480.)
       else if ( z < 5500.0 ) then
         qt = 4.2 - 1.2e-3 * (z - 2000.)
       else
         qt = 0.0
       end if

!      convert units to kg/kg

       qt = 0.001 * qt

!      convert qt from specific humidity to mixing ratio

       qt = qt / ( 1.0 - qt )

!      thetal

       if ( z < 520.0 ) then
         thetal = 298.7
       else if ( z < 1480.0 ) then
         thetal = 298.7 + (302.4 - 298.7)/(1480. - 520.)  * (z -  520.)
       else if ( z < 2000.0 ) then
         thetal = 302.4 + (308.2 - 302.4)/(2000. - 1480.) * (z - 1480.)
       else if ( z < 3000.0 ) then
         thetal = 308.2 +  3.65E-3 * (z - 2000.)
       else   ! ZNT 06/14/2021
         thetal = 311.85 + dthldz_up * (z - 3000.)
       end if

!      u geostrophic wind

       if ( z < 5500.0 ) then
         ug = -10.0 + 1.8e-3 * z
       else
         ug = 0.0
       end if

!      v geostrophic wind

       vg = 0.0

       return
       end subroutine bomex_snd

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


end module scm_bomex_mod
