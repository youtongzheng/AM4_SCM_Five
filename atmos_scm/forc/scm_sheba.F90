module scm_sheba_mod
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       Surface Heat Budget of the Arctic Ocean (SHEBA) FORCING MODULE
!
!       Jan 2011
!       Contact person: Huan Guo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This module provides the large scale forcing necessary
!       to run the single column model to simulate a mixed-phase Arctic stratus cloud 
!       on 07 May 1998 during SHEBA.
!
!       This case is derived from observation on May 7-8, 1998 during the Surface 
!       Heat Budget of the Arctic Ocean (SHEBA). This case was characterized by 
!       a single layer mixed-phase cloud at the top of a well-mixed boundary layer 
!       and minimum temperature of -20C above a sea ice-covered surface with weak 
!       turbulent surface fluxes. 
!
!       References
!       ----------
!
!       Morrison, H., J. O. Pinto, J. A. Curry, and G. M. McFarquhar (2008), Sensitivity 
!       of modeled arctic mixed-phase stratocumulus to cloud condensation and ice 
!       nuclei over regionally varying surface conditions, J. Geophys. Res., 113, D05203, 
!       doi:10. 1029/2007JD008729.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

   use sat_vapor_pres_mod, only:  lookup_es, compute_qs
   use polysvp_mod,        only:  polysvp_l,  polysvp_i
   use            mpp_mod, only:  mpp_pe, mpp_root_pe, stdlog
   use         mpp_io_mod, only:  mpp_open,MPP_RDONLY
   use            fms_mod, only:  write_version_number, open_namelist_file,  &
                                  check_nml_error,      file_exist,          &
                                  error_mesg,           close_file,   FATAL, &
                                  NOTE,                 read_data,    write_data, &
                                  mpp_error

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod
   use      constants_mod, only:  rdgas, cp_air, tfreeze, hlv, hls, rvgas, &
                                  grav, stefan, pi

   use   surface_flux_mod, only:  surface_flux

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

   public sheba_data_read, sheba_forc_init, sheba_forc_end, update_sheba_forc, &
          sheba_forc_diagnostic_init,                                          &
          get_sheba_flx, get_sheba_sst, ice_nucl_sheba

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following subroutines:
!
!            sheba_data_read      reads forcing data files and initializes
!                                 any needed parameters of the run
!            sheba_forc_init      initializes prognostic variables:           
!                                 T,u,v,qv,ql,qi
!            sheba_forc_end       deallocates allocated spaces
!            update_sheba_forc    a call to this routine returns the 
!                                 forcing parameters needed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              GLOBAL STORAGE VARIABLES
!
!  FORCING DATA (I.E. OBSERVED FIELDS)

!  LAYERED FIELDS
!       plev           pressure levels for initial and forcing data
!       temp_init      initial air temperature (K)
!       qv_init        initial vapor specific humidity (kg /kg moist air)) 
!       ql_init        initial liquid water specific humidity
!                       (kg liquid/kg moist air)
!       u_init         initial zonal wind (m/s)
!       v_init         initial meridional wind (m/s)
!       omega_init     vertical pressure velocity (Pa/s, positive indicates downward velocity)
!       t_forc         large-scale horizontal temperature advection (K/s)
!       qv_forc        large-scale horizontal water vapor advection (Kg/Kg/s)
!
!  SURFACE DATA
!       lat_forc       latitude of parcel (radians)
!       lon_forc       longitude of parcel (radians)

character(len=8) :: mod_name = 'sheba'
character(len=7) :: mod_name_diag = 'forcing'


REAL, ALLOCATABLE, DIMENSION(:)       :: plev,  temp_init, qv_init, ql_init, u_init, v_init, omega_init
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: u_init_sv, v_init_sv, omega_f_sv,  omega_h_sv
REAL, ALLOCATABLE, DIMENSION(:,:)     :: t_forc,    qv_forc
REAL, ALLOCATABLE, DIMENSION(:,:)     :: t_forc_sv, qv_forc_sv

!----------Diagnostic data----------------------------------------------
integer ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                           &
            id_qldt_vadv, id_qidt_vadv,                         & !h1g
            id_udt_vadv, id_udt_geos, id_udt_lf, id_udt_nudge,  &
            id_vdt_vadv, id_vdt_geos, id_vdt_lf, id_vdt_nudge,  &
            id_pf_forc, id_ph_forc, id_zf_forc, id_zh_forc,     &
            id_u_geos, id_v_geos

integer ::  id_qvdt_forc_col, id_qldt_vadv_col
integer ::  id_qadt_vadv,     id_qndt_vadv,     id_qnidt_vadv
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       PARAMETERS OF THE MODULE
!
!       ivert                    # of vertical levels of SHEBA initial
!                                  sounding data
!       itime                    # of times of mpace initial sounding
!
!       p00                      reference pressure (pascals)
!
!       tskin                    sea surface temperature (K)
!       SENFLUX                  sensible heat flux (W/m2)
!       EVAPFLUX                 evaporation flux (kg water/m2/sec)
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
!
!       p_omega_zero             pressure above which omega is set to
!                                to zero
!
!       p_cld_zero               pressure above which clouds are forced
!                                to be zero
!

INTEGER, PUBLIC                :: tracer_vert_advec_scheme = 3
INTEGER, PUBLIC                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3

logical, public                :: do_rad = .true.
logical, public                :: do_geo = .false.
logical, public                :: do_vadv = .true.
logical, public                :: do_nudge = .true.

INTEGER, PRIVATE               :: itime = 1
INTEGER, PRIVATE               :: ivert = 106
REAL,    PRIVATE               :: missing_value = -999.
REAL,    PRIVATE               :: p00 = 100000.
REAL,    PUBLIC                :: tskin = 257.4
REAL,    PUBLIC                :: SENFLUX = 7.98  !(w/m2)
REAL,    PUBLIC                :: EVAPFLUX = 2.86 !(w/m2)

real,    public                :: Ni_ini = 1.7e3   ! initial ice nuclei (#/m3)

real,    private               :: psfc = 1017e2
real,    private               :: zsfc = 0.0

real,    private               :: lat_forc = pi*76.0/180.
real,    private               :: lon_forc = pi*(360.-165.0)/180.

LOGICAL, PUBLIC                :: vert_advec_cond = .TRUE.
LOGICAL, PRIVATE               :: sheba_forc_initialized = .FALSE.
LOGICAL, PUBLIC                :: do_netcdf_restart = .true.

real, parameter                :: d622 = rdgas/rvgas

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!

character(len=128) :: Version = '$Id: scm_sheba.F90,v 1.1.2.1 2014/02/10 18:17:23 wfc Exp $'
character(len=128) :: Tagname = '$Name: new_gate_wfc $'
integer, dimension(1) :: restart_versions = (/ 1 /)
        
NAMELIST /scm_sheba_nml/ do_netcdf_restart,                          &
                         tracer_vert_advec_scheme,                   &
                         temp_vert_advec_scheme,                     &
                         momentum_vert_advec_scheme,                 &
                         vert_advec_cond,                            &
                         do_rad, do_geo, do_vadv,                    &
                         do_nudge,     Ni_ini 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CONTAINS


!#######################################################################
!#######################################################################

subroutine sheba_data_read(kmax)

implicit none
                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine reads in the SCM forcing data and puts it into
! allocated global storage variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,l,t            counting integers
!       unit                 unit number for I/O file
!       io,ierr              dummy integer variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Internal variables
!  ------------------
integer,  intent (in)                  :: kmax
INTEGER                                :: i,j,k,l,t,unit,io,ierr
CHARACTER*23                           :: tracer_ascheme,temp_ascheme
real                                   :: rh
real                                   :: esn

character*64                           :: fname_res='INPUT/sheba.res.nc'
     
       ! if the constructor is already called, then return, otherwise 
       ! set sheba_forc_initialized .TRUE.

       if (sheba_forc_initialized) return
       sheba_forc_initialized = .TRUE.

!    --------- read namelist --------
      
       if (file_exist('input.nml')) then
            unit = open_namelist_file()
            ierr=1; do while (ierr /= 0)
            read  (unit, nml=scm_sheba_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'scm_sheba_nml')
            enddo
10          call close_file (unit)
       endif

 !--------- write version number and namelist --------
       call write_version_number ( version, tagname )
       if(mpp_pe() == mpp_root_pe() ) write(stdlog(),nml=scm_sheba_nml)


       if (ALLOCATED(u_init_sv) )   DEALLOCATE (u_init_sv)
          ALLOCATE(u_init_sv(1,1, kmax))
       if (ALLOCATED(v_init_sv) )   DEALLOCATE (v_init_sv)
          ALLOCATE(v_init_sv(1,1, kmax))
       if (ALLOCATED(omega_f_sv) )   DEALLOCATE (omega_f_sv)
          ALLOCATE(omega_f_sv(1,1, kmax))
       if (ALLOCATED(omega_h_sv) )   DEALLOCATE (omega_h_sv)
          ALLOCATE(omega_h_sv(1,1, kmax+1))

       if (ALLOCATED(t_forc_sv) )   DEALLOCATE (t_forc_sv)
          ALLOCATE(t_forc_sv(kmax, itime))
       if (ALLOCATED(qv_forc_sv) )   DEALLOCATE (qv_forc_sv)
          ALLOCATE(qv_forc_sv(kmax, itime))

       u_init_sv  = 0.0
       v_init_sv  = 0.0
       t_forc_sv  = 0.0
       qv_forc_sv = 0.0
       if (file_exist('INPUT/sheba.res.nc') ) then
         if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_sheba_mod', &
         'Reading netCDF formatted restart file: INPUT/sheba.res.nc', NOTE)
         if( allocated( u_init_sv ) ) call read_data (fname_res, 'u_init_sv',  u_init_sv)
         if( allocated( v_init_sv ) ) call read_data (fname_res, 'v_init_sv',  v_init_sv)
         if( allocated( omega_f_sv) ) call read_data (fname_res, 'omega_f_sv', omega_f_sv)
         if( allocated( omega_h_sv) ) call read_data (fname_res, 'omega_h_sv', omega_h_sv)
         if( allocated( t_forc_sv ) ) call read_data (fname_res, 't_forc_sv',  t_forc_sv, no_domain=.true.)
         if( allocated( qv_forc_sv) ) call read_data (fname_res, 'qv_forc_sv', qv_forc_sv, no_domain=.true.)
       else
!-----------------------------------------------------------------------
!        ALLOCATE STORAGE
         if (ALLOCATED(plev) )   DEALLOCATE (plev)
           ALLOCATE(plev(ivert))
         if (ALLOCATED(temp_init) )   DEALLOCATE (temp_init)
           ALLOCATE(temp_init(ivert))
         if (ALLOCATED(qv_init) )   DEALLOCATE (qv_init)
           ALLOCATE(qv_init(ivert))
         if (ALLOCATED(ql_init) )   DEALLOCATE (ql_init)
           ALLOCATE(ql_init(ivert))
         if (ALLOCATED(u_init) )   DEALLOCATE (u_init)
           ALLOCATE(u_init(ivert))
         if (ALLOCATED(v_init) )   DEALLOCATE (v_init)
           ALLOCATE(v_init(ivert))
         if (ALLOCATED(omega_init) )   DEALLOCATE (omega_init)
           ALLOCATE(omega_init(ivert))
         if (ALLOCATED(t_forc) )   DEALLOCATE (t_forc)
           ALLOCATE(t_forc(ivert, itime))
         if (ALLOCATED(qv_forc) )   DEALLOCATE (qv_forc)
           ALLOCATE(qv_forc(ivert, itime))
!----------------------------------------------------------------------- 
!      READ IN SHEBA INITIAL SOUNDING and FORCING DATA
!
         call mpp_open(unit,file='INPUT/formatted.sheba.initial.sounding',&
                          action=MPP_RDONLY)
         do 1000 t = 1, itime
            read(unit,*)
            read(unit,*)
            read(unit,*)
            read(unit,*)
            read(unit,*)
            read(unit,*)

            do k = 1, ivert
                 read(unit,13) plev(k), temp_init(k), rh, ql_init(k), &
                               u_init(k), v_init(k), omega_init(k), t_forc(k,t), qv_forc(k,t)

! convert relative humidity to vapor specific humidity                
! Note: "lookup_es" calculates saturation vapor pressure with respect to mixed liquid and ice if temperature < 0C
!        "polysvp_l" calculates saturation vapor pressure with respect to liquid; 
                     
                 esn =  polysvp_l( temp_init(k) )
                 qv_init(k) = rh * 0.01 * esn
                 qv_init(k) = d622 * qv_init(k)/( plev(k) - (1.0- d622)*qv_init(k) )

! convert cloud liquid mixing ratio to specific content
                 ql_init(k) = ql_init(k) / (1.0 + ql_init(k) )
             enddo
13          format( 10(E15.5) )
1000     continue
         call close_file(unit)

       endif  !file_exist('INPUT/sheba.res.nc' )
        print*, 'after read   sheba.res.nc  '

end subroutine sheba_data_read


!#######################################################################
!#######################################################################


subroutine sheba_forc_init(time_interp, As, Bs)
#include "fv_arrays.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a given time and initializes the model 
!      variables for that time. This involves interpolating the global 
!      storage fields to a given time and pressure.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      VARIABLES
!
!      ------
!      INPUT:
!      ------
!
!      time_interp     time_type variable containing time we are inter-
!                      polating to.
!      As, Bs          A's and B's of half levels in hybrid coordinate
!                      phalf(k) = A(k) + B(k) * psurf
!
!      -------------
!      INPUT/OUTPUT:
!      -------------
!      T            temperature in degrees Kelvin 
!      u            zonal wind (m/s)
!      v            meridional wind (m/s)
!      qv           water vapor specific humidity (kg water vapor/     &
!                                                  kg air)
!      ql           liquid water condensate specific humidity 
!                   (kg condensate/ kg air)
!      qi           ice H2O condensate specific humidity 
!                   (kg condensate/kg dry air)
!      qa           cloud fraction (fraction)
!
!      -------------------
!      INTERNAL VARIABLES:
!      -------------------
!
!      k               counting integers
!      itime_more      indicates indice of TIME_TYPE for which 
!                      TIME_TYPE(itime_more)>= time_interp
!      itime_less      indicates indice of TIME_TYPE for which 
!                      TIME_TYPE(itime_more)<= time_interp
!      KDIM            no. of vertical levels to T array
!      weight_more     real number indicating closeness of interpolated
!                      time to TIME_TYPE(itime_more)
!      weight_less     real number indicating closeness of interpolated
!                      time to TIME_TYPE(itime_less)
!      plev_int        array of pressure coordinates to be passed to 
!                      interp_2d_field (Pa)
!      field_int_less  1D array of field at less time
!      field_int_more  1D array of field at more time 
!      tmp_p           temporary pressure levels for interpolating
!      tmp_ans         temporary array of answer from interpolating
!      ps_init         surface pressure initial value (Pa)
!      qsat            saturation specific humidity (kg vapor/kg air)
!      gamma           temporary variable        
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

TYPE(TIME_TYPE)                          :: time_interp
REAL,  INTENT (IN)   , DIMENSION(:)      :: As,Bs
REAL,  DIMENSION(size(pt,1),size(pt,2))    :: elev

!  Internal variables
!  ------------------

integer  :: kdim, k, klev

real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, lph
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf
real                                             :: weight_k,  weight_km1
integer :: i,j
#include "fv_point.inc"
       nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
       nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
       nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
       nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
       nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   !h1g
       nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )   !h1g
       
!   Initialize surface pressure and topography
    ps   = psfc
    elev = phis

!--- find out # of vertical levels
    KDIM = size(pt,3)

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
       do k = 1,KDIM
         do klev= ivert, 2, -1

          if( pf(1,1,k) .le. plev(klev) .and. pf(1,1,k) .gt. plev(klev-1) ) then
            weight_k   = ( pf(1,1,k) - plev(klev-1) ) / ( plev(klev)-plev(klev-1))
            weight_km1 = ( plev(klev)- pf(1,1,k)  ) / ( plev(klev)-plev(klev-1))
            pt(1,1,k)  = temp_init(klev) * weight_k + temp_init(klev-1) * weight_km1
            q(1,1,k,nsphum) = qv_init(klev)   * weight_k + qv_init(klev-1)   * weight_km1
            q(1,1,k,nql) = ql_init(klev)   * weight_k + ql_init(klev-1)   * weight_km1
            ua(1,1,k)  = u_init(klev)    * weight_k + u_init(klev-1)    * weight_km1
            va(1,1,k)  = v_init(klev)    * weight_k + v_init(klev-1)    * weight_km1

            omega_f_sv(1,1,k) = omega_init(klev) * weight_k + omega_init(klev-1) * weight_km1
            t_forc_sv( k, 1) =t_forc(klev,1) * weight_k + t_forc(klev-1,1) * weight_km1
            qv_forc_sv( k, 1)=qv_forc(klev,1)* weight_k + qv_forc(klev-1,1)* weight_km1
          endif

! above the observation height
          if( pf(1,1,k) .gt. plev(ivert) ) then
            pt(1,1,k)  = temp_init(ivert)
            q(1,1,k,nsphum) = qv_init(ivert)
            q(1,1,k,nql) = ql_init(ivert)
            ua(1,1,k)  = u_init(ivert)
            va(1,1,k)  = v_init(ivert)

            omega_f_sv(1,1,k) = omega_init(ivert)
            t_forc_sv( k, 1) = t_forc( ivert,1)
            qv_forc_sv( k, 1)= qv_forc(ivert,1)
          endif

! below the observation height
          if( pf(1,1,k) .lt. plev(1) ) then
            pt(1,1,k)  = temp_init(1)
            q(1,1,k,nsphum) = qv_init(1)
            q(1,1,k,nql) = ql_init(1)
            ua(1,1,k)  = u_init(1)
            va(1,1,k)  = v_init(1)

            omega_f_sv(1,1,k) = omega_init(1)
            t_forc_sv( k, 1) =t_forc( 1,1)
            qv_forc_sv( k, 1)=qv_forc(1,1)
          endif

         enddo ! klev
         pt(:,:,k)  = pt(1,1,k)
         q(:,:,k,nsphum) = q(1,1,k,nsphum)
         q(:,:,k,nql) = q(1,1,k,nql)
         ua(:,:,k)  = ua(1,1,k)
         va(:,:,k)  = va(1,1,k)

         u_init_sv(1,1,k)  = ua(1,1,k)
         v_init_sv(1,1,k)  = va(1,1,k)
       enddo ! k


! --- interpolation to half levels
      do k = 1,KDIM+1
         do klev= ivert, 2, -1

          if( ph(1,1,k) .le. plev(klev) .and. ph(1,1,k) .gt. plev(klev-1) ) then
            weight_k   = ( ph(1,1,k) - plev(klev-1) ) / ( plev(klev)-plev(klev-1))
            weight_km1 = ( plev(klev)- ph(1,1,k)  ) / ( plev(klev)-plev(klev-1))
            omega_h_sv(1,1,k) = omega_init(klev) * weight_k + omega_init(klev-1) * weight_km1
          endif

! above the observation height
          if( ph(1,1,k) .gt. plev(ivert) ) then
            omega_h_sv(1,1,k) = omega_init(ivert)
           endif

! below the observation height
          if( ph(1,1,k) .lt. plev(1) ) then
            omega_h_sv(1,1,k) = omega_init(1)
           endif

         enddo ! klev
       enddo ! k

!   Initialize cloud amount
    where ( q(:,:,:,nql) > 0.0 )
      q(:,:,:,nqa) = 1.0  
     elsewhere
      q(:,:,:,nqa) = 0.0
    endwhere
    q(:,:,:,nqi) = 0.0
    if( nqn > 0)  q(:,:,:,nqn) = 0.0
    if( nqni > 0) q(:,:,:,nqni)= 0.0

!-----------------------------------------------------------------------
 end subroutine sheba_forc_init



!#######################################################################
subroutine sheba_forc_end()
character*64                 :: fname_res='RESTART/sheba.res.nc'

  if (.NOT.sheba_forc_initialized) return
  sheba_forc_initialized = .FALSE.

  if (ALLOCATED(plev) )         DEALLOCATE (plev)
  if (ALLOCATED(temp_init) )    DEALLOCATE (temp_init)
  if (ALLOCATED(qv_init) )      DEALLOCATE (qv_init)
  if (ALLOCATED(ql_init) )      DEALLOCATE (ql_init)
  if (ALLOCATED(u_init) )       DEALLOCATE (u_init)
  if (ALLOCATED(v_init) )       DEALLOCATE (v_init)
  if (ALLOCATED(omega_init) )   DEALLOCATE (omega_init)
  if (ALLOCATED(t_forc) )       DEALLOCATE (t_forc)
  if (ALLOCATED(qv_forc) )      DEALLOCATE (qv_forc) 

  if( do_netcdf_restart ) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('scm_sheba_mod', 'Writing netCDF formatted restart file: RESTART/sheba.res.nc', NOTE)
           endif
           call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))), no_domain=.true.)   
           call write_data (fname_res,  'u_init_sv',  u_init_sv)
           call write_data (fname_res,  'v_init_sv',  v_init_sv)
           call write_data (fname_res,  'omega_f_sv', omega_f_sv)
           call write_data (fname_res,  'omega_h_sv', omega_h_sv)
           call write_data (fname_res,  't_forc_sv',  t_forc_sv, no_domain=.true.)
           call write_data (fname_res,  'qv_forc_sv', qv_forc_sv, no_domain=.true.)
   endif

  deallocate ( u_init_sv, v_init_sv )
  deallocate ( t_forc_sv, qv_forc_sv )
  
end subroutine sheba_forc_end

!#######################################################################


subroutine sheba_forc_diagnostic_init(axes, Time)

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
     'liquid droplet number concentration tendencies due to vertical advection', '1/kg/s', missing_value = missing_value)

id_qnidt_vadv = register_diag_field (mod_name_diag, 'qnidt_vadv', axes(1:3), Time, &
     'ice number concentration tendencies due to vertical advection', '1/kg/s', missing_value = missing_value)

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

id_qvdt_forc_col =  register_diag_field (mod_name_diag, 'qvdt_forc_col', axes(1:2), Time, &
    'column integrated vapor forcing', 'kg/m2/s',  missing_value = missing_value)

id_qldt_vadv_col =  register_diag_field (mod_name_diag, 'qldt_vadv_col', axes(1:2), Time, &
    'column integrated cloud water vertical advection', 'kg/m2/s',  missing_value = missing_value)
!-----------------------------------------------------------------------
end subroutine sheba_forc_diagnostic_init


!#######################################################################
!#######################################################################
subroutine update_sheba_forc(time_interp,time_diag,dt_int)
#include "fv_arrays.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a given time and returns the 
!      large-scale forcing appropriate for that given time.
!      This involves interpolating the global storage fields
!      to a given time an pressure, to do the resulting
!      calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!      VARIABLES
!
!      ------
!      INPUT:
!      ------
!
!      time_interp  time_type variable containing time we are
!                   interpolating to.
!      time_diag    diagnostic time to use for netcdf output
!      dt_int       time type variable containing the delta-t of the
!                   integration
!      Vgrid        vertical grid
!      ------------
!      INPUT/OUTPUT
!      ------------
!      T            temperature in degrees Kelvin 
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
!
!      ------
!      OUTPUT:
!      ------
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
!      omega_f       omega interpolated to full model levels 
!                       (Pa/sec)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

TYPE(TIME_TYPE), intent(in)              :: time_interp,time_diag,dt_int

!  Internal variables
!  ------------------
integer                                          :: i, j, k, kdim, klev
integer                                          :: dt_seconds,dt_days
logical                                          :: used
real                                             :: dts
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dp
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, omega_h
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dT_adi, dT_vadv, dT_lf, dT_rad, &
                                                    dqv_vadv, dqv_lf, dql_vadv, dqi_vadv, dqa_vadv, &
                                                    dqn_vadv,  dqni_vadv!h1g
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: du_vadv, dv_vadv, du_geos, dv_geos, du_lf, dv_lf, du_nudge, dv_nudge

real, dimension(size(pt,1),size(pt,2))  ::   qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::   qldt_vadv_col
real                                  ::   weight_k,  weight_km1
real                                  ::   tao_nudging

#include "fv_point.inc"

! ------------------------------------------------------------------------------------------------------------
! --- find out # of vertical levels
       KDIM = size(pt,3)
       ps   = psfc
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo

       dT_lf(1,1, :)   = t_forc_sv(:, 1)
       dqv_lf(1,1,:)   = qv_forc_sv(:, 1)
       omega_h(1,1,:)  = omega_h_sv(1,1,:)
       do k = 1, KDIM
          omga(:,:,k)  = omega_f_sv(1,1,k)
       enddo
       du_lf = 0.0
       dv_lf = 0.0

! --- compute dp, pi_fac, theta
do k = 2,KDIM+1
   dp(:,:,k-1) = ph(:,:,k) - ph(:,:,k-1)
enddo

! --- large-scale subsidence tendencies
dT_vadv=0.0; dT_adi=0.0; dqv_vadv=0.0; dql_vadv=0.0; dqi_vadv=0.0;dqa_vadv=0.0;
dqn_vadv=0.0;dqni_vadv=0.0;!h1g
du_vadv=0.0; dv_vadv=0.0;

if (do_vadv) then
   call get_time(dt_int,dt_seconds,dt_days)
   dts = real(dt_seconds + 86400*dt_days)

   dT_adi=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pf(:,:,:)

  select case (temp_vert_advec_scheme)
   case(1)
      vadvec_scheme=SECOND_CENTERED
   case(2)
      vadvec_scheme=FOURTH_CENTERED
   case(3)
      vadvec_scheme=FINITE_VOLUME_LINEAR
   case(4)
      vadvec_scheme=FINITE_VOLUME_PARABOLIC
   case(5)
      vadvec_scheme=SECOND_CENTERED_WTS
   case(6)
      vadvec_scheme=FOURTH_CENTERED_WTS
   end select

   call vert_advection(dts,omega_h,delp,pt,dT_vadv, &
                       scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   dT_vadv= dT_vadv + dT_adi

   select case (tracer_vert_advec_scheme)
   case(1)
      vadvec_scheme=SECOND_CENTERED
   case(2)
      vadvec_scheme=FOURTH_CENTERED
   case(3)
      vadvec_scheme=FINITE_VOLUME_LINEAR
   case(4)
      vadvec_scheme=FINITE_VOLUME_PARABOLIC
   case(5)
      vadvec_scheme=SECOND_CENTERED_WTS
   case(6)
      vadvec_scheme=FOURTH_CENTERED_WTS
   end select

   call vert_advection(dts,omega_h,delp,q(:,:,:,nsphum),dqv_vadv, &
                       scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   call vert_advection(dts,omega_h,delp,q(:,:,:,nql),dql_vadv, &
                       scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   call vert_advection(dts,omega_h,delp,q(:,:,:,nqi),dqi_vadv, &
                       scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   call vert_advection(dts,omega_h,delp,q(:,:,:,nqa),dqa_vadv, &
                       scheme=vadvec_scheme,form=ADVECTIVE_FORM)
   if( nqn > 0 ) &!h1g
     call vert_advection(dts,omega_h,delp,q(:,:,:,nqn),dqn_vadv, &
                         scheme=vadvec_scheme,form=ADVECTIVE_FORM)!h1g
   if( nqni > 0 ) &!h1g
     call vert_advection(dts,omega_h,delp,q(:,:,:,nqni),dqni_vadv, &
                         scheme=vadvec_scheme,form=ADVECTIVE_FORM)!h1g


   select case (momentum_vert_advec_scheme)
   case(1)
      call vert_advection(dts,omega_h,delp,ua,dU_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
      call vert_advection(dts,omega_h,delp,ua,dV_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
   end select

   end if


! --- radiative cooling
   dT_rad=0.;

! --- nudge winds
   if( do_nudge ) then
   ! nudging time scale is 1 hr
     tao_nudging = 3600.0 * 1.0
     do k=1, kdim
       du_nudge(1, 1, k) = ( u_init_sv(1, 1, k) - u(1, 1, k) ) / tao_nudging
       dv_nudge(1, 1, k) = ( v_init_sv(1, 1, k) - v(1, 1, k) ) / tao_nudging

       du_nudge(:, :, k) = du_nudge(1, 1, k)
       dv_nudge(:, :, k) = dv_nudge(1, 1, k)
     enddo
   else
     du_nudge = 0.0
     dv_nudge = 0.0
   endif

! --- compute geostrophic tendencies
   du_geos=0.; dv_geos=0.

! --- sum all tendencies
   u_dt = du_vadv + du_geos + du_lf + du_nudge
   v_dt = dv_vadv + dv_geos + dv_lf + dv_nudge
   t_dt = dT_rad + dT_vadv + dT_lf
   q_dt(:,:,:,nsphum) = dqv_vadv + dqv_lf
   q_dt(:,:,:,nql) = dql_vadv
   q_dt(:,:,:,nqi) = dqi_vadv
   q_dt(:,:,:,nqa) = dqa_vadv
   if( nqn  > 0 )  q_dt(:,:,:,nqn)  = dqn_vadv
   if( nqni > 0 ) q_dt(:,:,:,nqni) = dqni_vadv

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

if( nqn > 0) then
    if (id_qndt_vadv> 0) used = send_data( id_qndt_vadv,dqn_vadv(:,:,:), time_diag, 1, 1)
endif

if( nqni > 0 ) then
    if (id_qnidt_vadv> 0) used = send_data( id_qnidt_vadv,dqni_vadv(:,:,:), time_diag, 1, 1)
endif

if (id_pf_forc > 0)  used = send_data( id_pf_forc,  pf      (:,:,:), time_diag, 1, 1)

if (id_ph_forc > 0)  used = send_data( id_ph_forc,  ph      (:,:,:), time_diag, 1, 1)

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

if (id_tdt_radf > 0) used = send_data( id_tdt_radf, dT_rad  (:,:,:), time_diag, 1, 1)

if (id_udt_lf > 0)   used = send_data( id_udt_lf,   du_lf   (:,:,:), time_diag, 1, 1)

if (id_vdt_lf > 0)   used = send_data( id_vdt_lf,   dv_lf   (:,:,:), time_diag, 1, 1)

if (id_tdt_lf > 0)   used = send_data( id_tdt_lf,   dT_lf   (:,:,:), time_diag, 1, 1)

if (id_qvdt_lf > 0)  used = send_data( id_qvdt_lf,  dqv_lf  (:,:,:), time_diag, 1, 1)

! if (id_u_geos > 0)   used = send_data( id_u_geos,   u_geos  (:,:,:), time_diag, 1, 1)

! if (id_v_geos > 0)   used = send_data( id_v_geos,   v_geos  (:,:,:), time_diag, 1, 1)

if ( id_qvdt_forc_col > 0 )  used = send_data(  id_qvdt_forc_col, qvdt_forcing_col(:,:), time_diag, 1, 1 )

if ( id_qldt_vadv_col > 0 )  used = send_data(  id_qldt_vadv_col, qldt_vadv_col(:,:), time_diag, 1, 1 )

!-----------------------------------------------------------------------
end subroutine update_sheba_forc




!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_sheba_flx( ustar, flux_t, flux_q )

  implicit none
  real, intent(out), dimension(:) :: ustar, flux_t, flux_q
 ! ustar is estimated from Persson et al., JGR, Ocean, Vol 107, C10, 8045, 2002.  

  ustar  = 0.175
  flux_t = SENFLUX
  flux_q = EVAPFLUX

end subroutine get_sheba_flx

!########################################################################
! This subroutine returns imposed SST

subroutine get_sheba_sst( sst )

  implicit none
  real, intent(out) :: sst

  sst = tskin

end subroutine get_sheba_sst
!########################################################################




!########################################################################
subroutine   ice_nucl_sheba( rhi_in, Ni,   rh_crit_1d )
 implicit none
  real, intent(in)  ::   rhi_in
  real, intent(out) ::   Ni,  rh_crit_1d
  
  ! internal
     real, parameter ::    rhi_cri = 1.05
 ! ------------------------------------------------------------------------------
     if( rhi_in .gt.  rhi_cri ) then 
        Ni         = Ni_ini
        rh_crit_1d = rhi_in
     else
        Ni = 0.0
        rh_crit_1d = 1.0
     endif

return
end subroutine  ice_nucl_sheba
!########################################################################
!########################################################################



end module scm_sheba_mod
