module scm_twp_mod 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       scm forcing module --- compatible with all Minghua Zhang
!                              varational analysis data sets
!
!       September 2002
!       Contact person: Steve Klein
!
!       YLL: modified to run TWP-ICE, add the get_twpice_flx, get_twpice_sst
!       YLL: 2009/06/09, modified to run for TWP-ICE SCM intercomparison
!            tskin = 29 + 273.15, albedo = 0.07, sfc_flux, etc.   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This module provides the large scale forcing necessary
!       to run the single column model.  The large-scale forcing
!       read comes from the output of the variational analysis
!       of Zhang and Lin (JAS 1997); Zhang et al. (MWR 2001).  
!       Variational analysis of data has been performed currently 
!       only at the ARM CART Site in Oklahoma.  It is anticipated 
!       that other ARM sites will eventually have forcing data.
!
!       References
!
!       Zhang, M.-H. and J. L. Lin, 1997: Constrained variational anal-
!       ysis of sounding data based on column-integrated budgets of 
!       mass, heat, moisture, and momentum: Approach and application
!       to ARM measurements. J. Atmos. Sci., vol. 54, pp. 1503-1524.
!
!       Zhang, M.-H., J. L. Lin, R. T. Cederwall, J. J. Yio, and S. C.
!       Xie, 2001: Objective analysis of ARM IOP data: Method and
!       sensitivity. Mon. Wea. Rev., vol. 129, pp. 295-311.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

   use            mpp_mod, only:  stdlog
   use         mpp_io_mod, only:  mpp_open, mpp_close, MPP_RDONLY
   use            fms_mod, only:  mpp_pe, FATAL,   &
                                  close_file, write_version_number,    &
                                  file_exist, error_mesg,              &
                                  check_nml_error, open_restart_file,  &
                                  mpp_root_pe,                         &
                                  read_data, write_data,               &
                                  mpp_error, NOTE, field_size

#ifdef INTERNAL_FILE_NML
   use            mpp_mod, only: input_nml_file
#else
   use            fms_mod, only: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod, only:  time_type, get_time, get_date, set_date, &
                                  operator(-),  operator(<=), operator(>=), &
                                  operator(<), operator(>)
   use      constants_mod, only:  rdgas, cp_air, tfreeze, &
                                  hlv, rvgas, grav, stefan

   use   surface_flux_mod, only:  surface_flux

   use vert_advection_mod, only:  vert_advection, SECOND_CENTERED,     &
                                  FOURTH_CENTERED, SECOND_CENTERED_WTS,&
                                  FOURTH_CENTERED_WTS, ADVECTIVE_FORM, &
                                  FINITE_VOLUME_LINEAR,                &
                                  FINITE_VOLUME_PARABOLIC  
   use   field_manager_mod, only: MODEL_ATMOS
   use  tracer_manager_mod, only: get_tracer_index

use       constants_mod, only: kappa
use             fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                               endlon, rlonb, rlatb,  cold_start, ncnst, &
                               pnats, consv_te, ptop, fv_init, fv_domain, &
                               fv_end, change_time, p_var, restart_format, area, &
                               ak, bk, rlon, rlat, ng_d, nt_prog, get_eta_level

   implicit none
   private

   public twp_data_read, twp_forc_init, twp_forc_end, update_twp_forc, &
          twp_forc_diagnostic_init,                      &
          twp_surface_flux_loop, get_twp_sfc, get_twp_sst

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following subroutines:
!
!            twp_data_read     reads forcing data files and initializes
!                              any needed parameters of the run
!            twp_forc_init     initializes prognostic variables:           
!                              T,u,v,qv,ql,qi
!            twp_forc_end      deallocates allocated spaces
!            update_twp_forc   a call to this routine returns the 
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
character(len=8) :: mod_name = 'forcing'


REAL, ALLOCATABLE, DIMENSION(:,:) :: plev_stand,t_stand,q_stand,O3_stand
REAL, ALLOCATABLE, DIMENSION(:)   :: plev_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: t_forc,q_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: u_forc,v_forc,omega_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: tadv_tot_forc,tadv_hor_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: qadv_tot_forc,qadv_hor_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: tau_inv_forc,q1_obs,q2_obs,arscl
REAL, ALLOCATABLE, DIMENSION(:)   :: ps_forc,dps_dt_forc
REAL, ALLOCATABLE, DIMENSION(:)   :: sen_flux_forc,lat_flux_forc
REAL, ALLOCATABLE, DIMENSION(:)   :: tair_surf_forc,tskin_forc
REAL, ALLOCATABLE, DIMENSION(:)   :: alb_24hr, alb_day
REAL, ALLOCATABLE, DIMENSION(:)   :: gFq, q2_integral

REAL, ALLOCATABLE, DIMENSION(:)   :: precip,tground,rh,uwnd,vwnd,      &
                                     netraddnsfc,olr,swup_toa,swdn_toa,&
                                     low_cld_amt,mid_cld_amt,          &
                                     hgh_cld_amt,tot_cld_amt,cld_thk,  &
                                     cld_hgt, s_surface, qs_surface,   &
                                     wvp, lwp, qtend, advqtend,        &
                                     evap, ttend, advttend, radttend,  &
                                     latttend, swdn_sfc, swup_sfc,     &
                                     lwdn_sfc, lwup_sfc, albedo_sfc,   &
                                     avg_albedo_sfc,coltcal,colqcal,   &
                                     phalf_forc

TYPE(TIME_TYPE), ALLOCATABLE, DIMENSION(:)  :: time_layer,time_sfc
TYPE(TIME_TYPE),  PUBLIC                    :: Time_var_start

!----------Diagnostic data----------------------------------------------

     integer          :: id_temp_obs, id_qv_obs, id_uwnd_obs,          &
                         id_vwnd_obs, id_omega, id_q1_obs, id_arscl,   &
                         id_q2_obs, id_tair_obs, id_tground_obs,       &
                         id_rh_obs, id_uref_obs, id_vref_obs,          &
                         id_precip_obs, id_latent_obs, id_sens_obs,    &
                         id_netdnsfc_obs, id_olr_obs,                  &
                         id_swup_toa_obs, id_swdn_toa_obs,             &
                         id_lwup_sfc_obs, id_lwdn_sfc_obs,             &
                         id_swup_sfc_obs, id_swdn_sfc_obs,             &
                         id_tot_cld_amt_obs, id_lwp_obs, id_pw_obs,    &
                         id_low_cld_amt_obs, id_mid_cld_amt_obs,       &
                         id_hgh_cld_amt_obs, id_cld_thk_obs,           &
                         id_cld_hgt_obs, id_qtend_obs, id_advqtend_obs,&
                         id_evap_obs, id_ttend_obs, id_advttend_obs,   &
                         id_radttend_obs,id_latttend_obs,              &
                         id_albedo_sfc_obs, id_nudg_t, id_nudg_q,      &
                         id_nudg_u, id_nudg_v,id_tdt_lf,id_qvdt_lf,    &
                         id_Tadv_hor_forc, id_vadvT, id_adiT, id_vadvQ,&
                         id_Qadv_hor_forc

!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         INITIALIZE PARAMETERS OF THE MODULE
!
!       ivert_stand  # of vertical levels to standard atmosphere
!       ivert        # of vertical levels of VAR data
!       itime        # of times of VAR data
!
!       p00          reference pressure (Pa)
!       rough_mom    momentum roughness length characteristic for TWP-ICE
!                    site taken from Andy Brown's set up for TWP-ICE
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
!           N               nudging only
!           S          semi-prognistic. Tendencies are zero.
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

CHARACTER*1, PUBLIC            :: FORCMETHOD = 'D'
INTEGER, PUBLIC                :: tracer_vert_advec_scheme = 3
INTEGER, PUBLIC                :: temp_vert_advec_scheme = 1
INTEGER, PRIVATE               :: itime, ivert
INTEGER, PRIVATE               :: nfields = 0
INTEGER, PRIVATE               :: nlayer_fields = 0
!INTEGER, PRIVATE               :: max_itime = 792
!!INTEGER, PRIVATE               :: max_ivert = 42
!INTEGER, PRIVATE               :: max_ivert = 107 ! for up to 10hPa forcing data
INTEGER, PRIVATE               :: ivert_stand = 60
REAL,    PRIVATE               :: p00 = 100000.
REAL,    PUBLIC                :: rough_mom = 0.035
REAL,    PUBLIC                :: rough_heat = 0.035 / 7.4
REAL,    PUBLIC                :: TSKIN = 302.15
REAL,    PUBLIC                :: SENFLUX = 10.
REAL,    PUBLIC                :: EVAPFLUX = 100./hlv
REAL,    PUBLIC                :: ALBEDO_OBS = 0.2
REAL,    PUBLIC                :: wind_relaxation_tau = 6.*3600.
REAL,    PUBLIC                :: relaxation_tau = 6.*3600.
REAL,    PUBLIC                :: missing_value = -999.
REAL,    PUBLIC                :: a1_precfrac = 1.0 ! The fraction of precipitation that we want the Donner a1 factor to capture.
LOGICAL, PUBLIC                :: vert_advec_cond = .TRUE.
LOGICAL, PUBLIC                :: do_qc_adv_loss = .FALSE.
LOGICAL, PUBLIC                :: do_iop_mean_albedo = .FALSE.
LOGICAL, PUBLIC                :: use_given_omega_surface = .FALSE.
LOGICAL, PUBLIC                :: do_netcdf_restart = .TRUE.

!--------------------- version number ----------------------------------
!
        
character(len=128) :: Version = '$Id: scm_twpice.F90,v 1.1.2.2.4.1 2014/07/30 18:04:56 wfc Exp $'
character(len=128) :: Tagname = '$Name: a1_closure_wfc $'
        

NAMELIST /SCM_twp_NML/  do_netcdf_restart,                             &
                        FORCMETHOD, tracer_vert_advec_scheme,          &
                        temp_vert_advec_scheme, vert_advec_cond,       &
                        do_qc_adv_loss, do_iop_mean_albedo,            &
                        relaxation_tau, wind_relaxation_tau,           &
                        use_given_omega_surface, a1_precfrac

LOGICAL :: twp_forc_initialized = .FALSE.

real, parameter   :: d622 = rdgas/rvgas
real, parameter   :: d378 = 1.0-d622
real, parameter   :: d608 = d378/d622
real, parameter   :: INIT_DATA = -9.999e+99 !-123456789.0

integer :: nsphum, nql, nqi, nqa, nqn, nqni, nq2_clo, nprec_clo
integer :: vadvec_scheme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS


!#######################################################################
!#######################################################################

subroutine twp_data_read()

implicit none

                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine reads in the SCM forcing data and puts it into
! allocated global storage variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k                counting integers
!       unit                 unit number for I/O file
!       io                   dummy variable
!       adum                 dummy character variable
!       beg_year_date        date of beginning of the year
!       dum_var,dum_mat_var  dummy arrays
!       year_var,month_var,  
!       day_var,hour_var,    year, month,day,hour,minute from VAR files
!       minute_var
!       nvar_var             # of variables in VAR layer file
!       varname              variable name in VAR layer file
!       vardat               layer data, tmp variable
!       fname                file name
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Internal variables
!  ------------------

INTEGER                                   :: i,j,k
INTEGER                                   :: unit,io,ierr,jstart,jend
CHARACTER*20                              :: adum
CHARACTER*23                              :: tracer_ascheme,temp_ascheme
TYPE(TIME_TYPE)                           :: beg_year_date
INTEGER                                   :: nvar_var, kmax
CHARACTER*128                              :: varname,fname
INTEGER                                   :: days, seconds, logunit, t
REAL                                      :: tmp1, tmp2, levthick
REAL, ALLOCATABLE, DIMENSION(:)           :: dum_var,year_var,month_var 
REAL, ALLOCATABLE, DIMENSION(:)           :: day_var,hour_var,minute_var
REAL, ALLOCATABLE, DIMENSION(:,:)         :: vardat
REAL, ALLOCATABLE, DIMENSION(:,:)         :: dum_mat_var


CHARACTER(len=64)                         :: fname_nc='INPUT/twp_forc.res.nc'
CHARACTER(len=64)                         :: forcing_nc='INPUT/varanal_forcing.nc'
CHARACTER(len=64)                         :: layer_name='INPUT/layer.dat'
CHARACTER(len=64)                         :: surf_name='INPUT/surface.dat'

INTEGER, DIMENSION(4)                     :: siz
logical                                   :: ASCII_VARA
!-----------------------------------------------------------------------
!      ----- read namelist -----

      ! if the constructor is already called, then return, otherwise 
      ! set twp_forc_initialized .TRUE.

      if (twp_forc_initialized) return
      twp_forc_initialized = .TRUE.


      if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=SCM_TWP_NML, iostat=io)
        ierr = check_nml_error(io, 'SCM_TWP_NML')
#else
        unit = open_namelist_file ()
        io=1; do while (io .ne. 0)
          read  (unit, nml=SCM_TWP_NML, iostat=io, end=10)
        enddo
  10    call close_file (unit)
#endif

!--------- write version and namelist to standard log ------------

        call write_version_number ( version, tagname )
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() ) &
          write ( logunit, nml=SCM_TWP_NML )

      endif

      !dummy check -------         
       if (temp_vert_advec_scheme .gt. 6 ) then
            write (adum,'(i1)') temp_vert_advec_scheme
            call error_mesg ( 'twp_data_read in scm_twp_mod',         &
                 'Bad temp_vert_advec_scheme : '//adum, FATAL )
       end if

       if (tracer_vert_advec_scheme .gt. 6 ) then
            write (adum,'(i1)') tracer_vert_advec_scheme
            call error_mesg ( 'twp_data_read in scm_twp_mod',         &
                 'Bad tracer_vert_advec_scheme : '//adum, FATAL )
       end if

       SELECT CASE(tracer_vert_advec_scheme)
            CASE(1)
                 write (tracer_ascheme,'(a)') 'SECOND_CENTERED        '
            CASE(2)
                 write (tracer_ascheme,'(a)') 'FOURTH_CENTERED        '
            CASE(3)
                 write (tracer_ascheme,'(a)') 'FINITE_VOLUME_LINEAR   '
            CASE(4)
                 write (tracer_ascheme,'(a)') 'FINITE_VOLUME_PARABOLIC'
            CASE(5)
                 write (tracer_ascheme,'(a)') 'SECOND_CENTERED_WTS    '
            CASE(6)
                 write (tracer_ascheme,'(a)') 'FOURTH_CENTERED_WTS    '
       END SELECT
       
       SELECT CASE(temp_vert_advec_scheme)
            CASE(1)
                 write (temp_ascheme,'(a)') 'SECOND_CENTERED        '
            CASE(2)
                 write (temp_ascheme,'(a)') 'FOURTH_CENTERED        '
            CASE(3)
                 write (temp_ascheme,'(a)') 'FINITE_VOLUME_LINEAR   '
            CASE(4)
                 write (temp_ascheme,'(a)') 'FINITE_VOLUME_PARABOLIC'
            CASE(5)
                 write (temp_ascheme,'(a)') 'SECOND_CENTERED_WTS    '
            CASE(6)
                 write (temp_ascheme,'(a)') 'FOURTH_CENTERED_WTS    '
       END SELECT
      
       unit = stdlog()
       if ( mpp_pe() == 0 ) then
           write (unit,'(/,80("="),/(a))') trim(version), trim(Tagname)
           write (unit,'(a15,a1)')  ' forcmethod =  ',forcmethod
           write (unit,'(a28,a23)') ' tracer_vert_advec_scheme = ',    &
                                      tracer_ascheme
           write (unit,'(a26,a23)') ' temp_vert_advec_scheme = ',      &
                                      temp_ascheme
           write (unit,*) ' vert_advec_cond    = ',vert_advec_cond
           write (unit,*) ' do_qc_adv_loss     = ',do_qc_adv_loss
           write (unit,*) ' do_iop_mean_albedo = ',do_iop_mean_albedo   
           write (unit,*) ' use_given_omega_surface = ',               &
                            use_given_omega_surface
       endif      
       call close_file (unit)
      
       !set up number of levels and number of times and number of surface
       !fields according to the data being used by reading ivert, itime,
       !and nfields from file
       
       ASCII_VARA = .true. ! Default to reading variational forcing data from ascii files.
       
       if( file_exist(trim(forcing_nc)) ) then
         call field_size (forcing_nc, 'time', siz)
         itime= siz(1)
         call field_size (forcing_nc, 'lev', siz)
         ivert= siz(1)
         ! For TWP_ICE, assume the number of fields is 43.
         nfields = 43
         ! For TWP_ICE, we assume the number of layer fields is 19 (needed for arscl variable).
         nlayer_fields = 19
         ASCII_VARA = .FALSE.
       endif
       if( file_exist(trim(layer_name)) .and. ASCII_VARA ) then
         call mpp_open(unit,layer_name,action=MPP_RDONLY)       
         read(unit,*)
         read(unit,111) itime
111      format(9X,i3)
         read(unit,*)
         read(unit,112) ivert
112      format(5X,i3)
         call mpp_close(unit)  
       
         if (mpp_pe() == mpp_root_pe()) then
           print *, 'itime = ', itime, 'ivert= ', ivert
         endif
       else
! ASCII Variational data is missing but we are expecting to be reading it.
         if (ASCII_VARA ) &
           call mpp_error ('scm_twp_mod',                            &
                         trim(layer_name)//' variational forcing file is missing.', &
                          FATAL )
       endif

       if( file_exist(trim(surf_name))  .and. ASCII_VARA ) then
         call mpp_open(unit,surf_name,action=MPP_RDONLY)
         read(unit,*)
         read(unit,211) nfields
211      format(5X,i3)       
         call mpp_close(unit)
       else
! ASCII Variational data is missing but we are expecting to be reading it.
         if (ASCII_VARA ) &
           call mpp_error ('scm_twp_mod',                            &
                           trim(surf_name)//' variational forcing file is missing.', &
                            FATAL )
       endif

!-----------------------------------------------------------------------
!
!      ALLOCATE STORAGE
!


       ALLOCATE(dum_var(itime),year_var(itime),month_var(itime) )
       ALLOCATE(day_var(itime),hour_var(itime),minute_var(itime))
       ALLOCATE(vardat(ivert,itime))
       ALLOCATE(dum_mat_var(ivert,itime))

       if (ALLOCATED(plev_stand)) DEALLOCATE (plev_stand)
           ALLOCATE(plev_stand(ivert_stand,2)) ; plev_stand = INIT_DATA
       if (ALLOCATED(t_stand)) DEALLOCATE (t_stand)
           ALLOCATE(t_stand(ivert_stand,2)) ; t_stand = INIT_DATA
       if (ALLOCATED(q_stand)) DEALLOCATE (q_stand)
           ALLOCATE(q_stand(ivert_stand,2)) ; q_stand = INIT_DATA
       if (ALLOCATED(O3_stand)) DEALLOCATE (O3_stand)
           ALLOCATE(O3_stand(ivert_stand,2)) ; O3_stand = INIT_DATA
       
       if (ALLOCATED(plev_forc)) DEALLOCATE (plev_forc)
           ALLOCATE(plev_forc(ivert)) ; plev_forc = INIT_DATA
       if (ALLOCATED(time_layer)) DEALLOCATE (time_layer)
           ALLOCATE(time_layer(itime))
       if (ALLOCATED(time_sfc)) DEALLOCATE (time_sfc)
           ALLOCATE(time_sfc(itime))
       
       if (ALLOCATED(t_forc)) DEALLOCATE (t_forc)
           ALLOCATE(t_forc(ivert,itime)) ; t_forc = INIT_DATA
       if (ALLOCATED(q_forc)) DEALLOCATE (q_forc)
           ALLOCATE(q_forc(ivert,itime)) ; q_forc = INIT_DATA
       if (ALLOCATED(u_forc)) DEALLOCATE (u_forc)
           ALLOCATE(u_forc(ivert,itime)) ; u_forc = INIT_DATA
       if (ALLOCATED(v_forc)) DEALLOCATE (v_forc)
           ALLOCATE(v_forc(ivert,itime)) ; v_forc = INIT_DATA
       if (ALLOCATED(omega_forc)) DEALLOCATE (omega_forc)
           ALLOCATE(omega_forc(ivert,itime)) ; omega_forc = INIT_DATA
       if (ALLOCATED(tadv_tot_forc)) DEALLOCATE (tadv_tot_forc)
           ALLOCATE(tadv_tot_forc(ivert,itime)) ; tadv_tot_forc = INIT_DATA
       if (ALLOCATED(qadv_tot_forc)) DEALLOCATE (qadv_tot_forc)
           ALLOCATE(qadv_tot_forc(ivert,itime)) ; qadv_tot_forc = INIT_DATA
       if (ALLOCATED(tadv_hor_forc)) DEALLOCATE (tadv_hor_forc)
           ALLOCATE(tadv_hor_forc(ivert,itime)) ; tadv_hor_forc = INIT_DATA
       if (ALLOCATED(qadv_hor_forc)) DEALLOCATE (qadv_hor_forc)
           ALLOCATE(qadv_hor_forc(ivert,itime)) ; qadv_hor_forc = INIT_DATA
       if (ALLOCATED(tau_inv_forc)) DEALLOCATE (tau_inv_forc)
           ALLOCATE(tau_inv_forc(ivert,itime)) ; tau_inv_forc = INIT_DATA

       if (ALLOCATED(ps_forc)) DEALLOCATE (ps_forc)
           ALLOCATE(ps_forc(itime)) ; ps_forc = INIT_DATA
       if (ALLOCATED(dps_dt_forc)) DEALLOCATE (dps_dt_forc)
           ALLOCATE(dps_dt_forc(itime)) ; dps_dt_forc = INIT_DATA
       if (ALLOCATED(sen_flux_forc)) DEALLOCATE (sen_flux_forc)
           ALLOCATE(sen_flux_forc(itime)) ; sen_flux_forc = INIT_DATA
       if (ALLOCATED(lat_flux_forc)) DEALLOCATE (lat_flux_forc)
           ALLOCATE(lat_flux_forc(itime)) ; lat_flux_forc = INIT_DATA
       if (ALLOCATED(tair_surf_forc)) DEALLOCATE (tair_surf_forc)
           ALLOCATE(tair_surf_forc(itime)) ; tair_surf_forc = INIT_DATA
       
       !allocate space for netcdf variables
       if (ALLOCATED(q1_obs)) DEALLOCATE (q1_obs)
           ALLOCATE(q1_obs(ivert,itime)) ; q1_obs = INIT_DATA
       if (ALLOCATED(q2_obs)) DEALLOCATE (q2_obs)
           ALLOCATE(q2_obs(ivert,itime)) ; q2_obs = INIT_DATA
       if (ALLOCATED(arscl)) DEALLOCATE (arscl)
           ALLOCATE(arscl(ivert,itime)) ; arscl = INIT_DATA
       if (ALLOCATED(precip)) DEALLOCATE (precip)
           ALLOCATE(precip(itime)) ; precip = INIT_DATA
       if (ALLOCATED(tground)) DEALLOCATE (tground)
           ALLOCATE(tground(itime)) ; tground = INIT_DATA
       if (ALLOCATED(rh)) DEALLOCATE (rh)
           ALLOCATE(rh(itime)) ; rh = INIT_DATA
       if (ALLOCATED(uwnd)) DEALLOCATE (uwnd)
           ALLOCATE(uwnd(itime)) ; uwnd = INIT_DATA
       if (ALLOCATED(vwnd)) DEALLOCATE (vwnd)
           ALLOCATE(vwnd(itime)) ; vwnd = INIT_DATA
       if (ALLOCATED(netraddnsfc)) DEALLOCATE (netraddnsfc)
           ALLOCATE(netraddnsfc(itime)) ; netraddnsfc = INIT_DATA

       if (ALLOCATED(olr)) DEALLOCATE (olr)
           ALLOCATE(olr(itime)) ; olr = INIT_DATA
       if (ALLOCATED(swup_toa)) DEALLOCATE (swup_toa)
           ALLOCATE(swup_toa(itime)) ; swup_toa = INIT_DATA
       if (ALLOCATED(swdn_toa)) DEALLOCATE (swdn_toa)
           ALLOCATE(swdn_toa(itime)) ; swdn_toa = INIT_DATA
       if (ALLOCATED(low_cld_amt)) DEALLOCATE (low_cld_amt)
           ALLOCATE(low_cld_amt(itime)) ; low_cld_amt = INIT_DATA

       if (ALLOCATED(mid_cld_amt)) DEALLOCATE (mid_cld_amt)
           ALLOCATE(mid_cld_amt(itime)) ; mid_cld_amt = INIT_DATA
       if (ALLOCATED(hgh_cld_amt)) DEALLOCATE (hgh_cld_amt)
           ALLOCATE(hgh_cld_amt(itime)) ; hgh_cld_amt = INIT_DATA
       if (ALLOCATED(tot_cld_amt)) DEALLOCATE (tot_cld_amt)
           ALLOCATE(tot_cld_amt(itime)) ; tot_cld_amt = INIT_DATA

       if (ALLOCATED(cld_thk)) DEALLOCATE (cld_thk)
           ALLOCATE(cld_thk(itime)) ; cld_thk = INIT_DATA
       if (ALLOCATED(cld_hgt)) DEALLOCATE (cld_hgt)
           ALLOCATE(cld_hgt(itime)) ; cld_hgt = INIT_DATA
       if (ALLOCATED(lwp)) DEALLOCATE (lwp)
           ALLOCATE(lwp(itime)) ; lwp = INIT_DATA

       if (ALLOCATED(qtend)) DEALLOCATE (qtend)
           ALLOCATE(qtend(itime)) ; qtend = INIT_DATA
       if (ALLOCATED(advqtend)) DEALLOCATE (advqtend)
           ALLOCATE(advqtend(itime)) ; advqtend = INIT_DATA
       if (ALLOCATED(evap)) DEALLOCATE (evap)
           ALLOCATE(evap(itime)) ; evap = INIT_DATA

       if (ALLOCATED(ttend)) DEALLOCATE (ttend)
           ALLOCATE(ttend(itime)) ; ttend = INIT_DATA
       if (ALLOCATED(advttend)) DEALLOCATE (advttend)
           ALLOCATE(advttend(itime)) ; advttend = INIT_DATA
       if (ALLOCATED(radttend)) DEALLOCATE (radttend)
           ALLOCATE(radttend(itime)) ; radttend = INIT_DATA
       if (ALLOCATED(latttend)) DEALLOCATE (latttend)
           ALLOCATE(latttend(itime)) ; latttend = INIT_DATA
       if (ALLOCATED(coltcal)) DEALLOCATE (coltcal)
           ALLOCATE(coltcal(itime)) ; coltcal = INIT_DATA
       if (ALLOCATED(colqcal)) DEALLOCATE (colqcal)
           ALLOCATE(colqcal(itime)) ; colqcal = INIT_DATA
       
       if (ALLOCATED(gFq)) DEALLOCATE (gFq)
           ALLOCATE(gFq(itime)) ; gFq = INIT_DATA
       if (ALLOCATED(q2_integral)) DEALLOCATE (q2_integral)
           ALLOCATE(q2_integral(itime)) ; q2_integral = INIT_DATA
       if (ALLOCATED(phalf_forc)) DEALLOCATE (phalf_forc)
           ALLOCATE(phalf_forc(ivert+1)) ; phalf_forc = INIT_DATA

       if (nfields.ge.43) then
       if (ALLOCATED(swup_sfc)) DEALLOCATE (swup_sfc)
           ALLOCATE(swup_sfc(itime)) ; swup_sfc = INIT_DATA
       if (ALLOCATED(swdn_sfc)) DEALLOCATE (swdn_sfc)
           ALLOCATE(swdn_sfc(itime)) ; swdn_sfc = INIT_DATA
       if (ALLOCATED(lwdn_sfc)) DEALLOCATE (lwdn_sfc)
           ALLOCATE(lwdn_sfc(itime)) ; lwdn_sfc = INIT_DATA
       if (ALLOCATED(lwup_sfc)) DEALLOCATE (lwup_sfc)
           ALLOCATE(lwup_sfc(itime)) ; lwup_sfc = INIT_DATA
       if (ALLOCATED(s_surface)) DEALLOCATE (s_surface)
          ALLOCATE(s_surface(itime)) ; s_surface = INIT_DATA
       if (ALLOCATED(qs_surface)) DEALLOCATE (qs_surface)
          ALLOCATE(qs_surface(itime))   ; qs_surface = INIT_DATA
       if (ALLOCATED(wvp)) DEALLOCATE (wvp)
           ALLOCATE(wvp(itime)) ; wvp = INIT_DATA
       end if              
   
       if (ALLOCATED(albedo_sfc)) DEALLOCATE (albedo_sfc)
           ALLOCATE(albedo_sfc(itime)) ; albedo_sfc = INIT_DATA
       if (ALLOCATED(avg_albedo_sfc)) DEALLOCATE (avg_albedo_sfc)
           ALLOCATE(avg_albedo_sfc(1)) ; avg_albedo_sfc = INIT_DATA

!       if (nfields.ge.44) then
       if (ALLOCATED(tskin_forc)) DEALLOCATE (tskin_forc)
           ALLOCATE(tskin_forc(itime)) ; tskin_forc = INIT_DATA
!       end if
       
       if (nfields.ge.45) then
       if (ALLOCATED(alb_24hr)) DEALLOCATE (alb_24hr)
           ALLOCATE(alb_24hr(itime)) ; alb_24hr = INIT_DATA
       end if
       
       if (nfields.ge.46) then
       if (ALLOCATED(alb_day)) DEALLOCATE (alb_day)
           ALLOCATE(alb_day(itime)) ; alb_day = INIT_DATA
       end if
       
!-----------------------------------------------------------------------
! 
!     READ IN TROPICAL ATMOSHPHERE DATA
!

       !----------------------------
       ! read summer

       fname='INPUT/tropic.dat'
       call mpp_open(unit,fname,action=MPP_RDONLY)
      
       !skip to beginning of data
       do j=1,4
            read(unit,'(a)') adum
       enddo
            
       !read data
       do j=1,ivert_stand
            read (unit,'(0p,2f10.3,1p,2e12.4)') &
                 plev_stand(ivert_stand+1-j,2), t_stand(ivert_stand+1-j,2),&
                 q_stand(ivert_stand+1-j,2),O3_stand(ivert_stand+1-j,2)  
       end do
       call mpp_close(unit)
      
!----- since only read in one file, let 

       plev_stand(:,1) = plev_stand(:,2)
       t_stand(:,1) = t_stand(:,2)
       q_stand(:,1) = q_stand(:,2)
       O3_stand(:,1) = O3_stand(:,2)
 
       !----record which file was read
       unit= stdlog()
       if (mpp_pe() == 0 ) Write (unit,'(a)') 'used data file: '// fname 
       call close_file(unit)
       
      
       !change units
       
       !(a) mb to Pascals
       plev_stand(:,:)=100.*plev_stand(:,:)
             
       !(b) mass mixing ratio to specific humidity
       q_stand(:,:)=q_stand(:,:)/(1.+q_stand(:,:))
       O3_stand(:,:)=O3_stand(:,:)/(1.+O3_stand(:,:))
           
       
!-----------------------------------------------------------------------
! 
!      READ IN VARIATIONAL FORCING DATA
!
!
!       

!----- read layer files -----
              
       if( ASCII_VARA ) then
         call mpp_open(unit,layer_name,action=MPP_RDONLY)
         read(unit,*)
         read(unit,*)
         read(unit,*)
         read(unit,*)
         read(unit,*)
         read(unit,12) plev_forc
12       format(5e15.7)
         read(unit,*)
         read(unit,12) dum_var(1:itime)
         read(unit,*)
         read(unit,12) year_var(1:itime)
         read(unit,*)
         read(unit,12) month_var(1:itime)
         read(unit,*)
         read(unit,12) day_var(1:itime)
         read(unit,*)
         read(unit,12) hour_var(1:itime)
         read(unit,*)
         read(unit,12) minute_var(1:itime)
         read(unit,*)
         read(unit,*) nvar_var
         nlayer_fields = nvar_var
       
         !read data
         do i=1,nvar_var

            read(unit,11)varname
11          format(a128)
            do j=1,ivert
                 read(unit,12)(vardat(j,k),k=1,itime)
            enddo
            if (mpp_pe() == mpp_root_pe()) then
               print *, 'var = ', varname
               print *, 'sample value = ', vardat(1,1)
            endif
            if (i.eq.1) t_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.2) q_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.3) u_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.4) v_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.5) omega_forc(:,:)=vardat(1:ivert,1:itime)
              
            !i.eq.7 = horizontal temperature advection
            !i.eq.12 = horizontal dry static energy advection
            if (i.eq.7) tadv_hor_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.13) tadv_tot_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.9)  qadv_hor_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.10) qadv_tot_forc(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.17) q1_obs(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.18) q2_obs(:,:)=vardat(1:ivert,1:itime)
            if (i.eq.19) arscl(:,:)=vardat(1:ivert,1:itime)  
         enddo

         call mpp_close(unit)
       
         !----record which file was read
         unit= stdlog()
         if (mpp_pe() == 0 ) Write (unit,'(a)') 'used data file: '//layer_name
         call close_file (unit)

    
!----- read surface file -----  

         call mpp_open(unit,surf_name,action=MPP_RDONLY)
       
         read(unit,*)
         read(unit,*)      
         read(unit,*)
         read(unit,*)
             
         !read data
         do i=1,nfields

            read(unit,11)  varname
            read(unit,12) (vardat(1,k),k=1,itime)

            if (mpp_pe() == mpp_root_pe()) then
              print *, varname
              print *, vardat(1,1)
            endif

            if (i.eq.7)  precip(:)         = vardat(1,1:itime)
            if (i.eq.8)  lat_flux_forc(:)  = vardat(1,1:itime)
            if (i.eq.9)  sen_flux_forc(:)  = vardat(1,1:itime)
            if (i.eq.10) ps_forc(:)        = vardat(1,1:itime)
! T_srf and T_skin data in file is in degrees C. 
! They will be converted to K later.
            if (i.eq.12) tair_surf_forc(:) = vardat(1,1:itime)
            if (i.eq.13) tground(:)        = vardat(1,1:itime)
            if (i.eq.14) rh(:)             = vardat(1,1:itime)
            if (i.eq.16) uwnd(:)           = vardat(1,1:itime)
            if (i.eq.17) vwnd(:)           = vardat(1,1:itime)
            if (i.eq.18) netraddnsfc(:)    = vardat(1,1:itime)
            if (i.eq.19) olr(:)            = vardat(1,1:itime)      
            if (i.eq.20) swup_toa(:)       = vardat(1,1:itime)           
            if (i.eq.21) then
                 swdn_toa(:)               = vardat(1,1:itime)
                 !convert netdn toa to upward reflected
                 swup_toa(:) = swdn_toa(:) - swup_toa(:)
            end if
            if (i.eq.22) low_cld_amt(:)    = vardat(1,1:itime)
            if (i.eq.23) mid_cld_amt(:)    = vardat(1,1:itime)
            if (i.eq.24) hgh_cld_amt(:)    = vardat(1,1:itime)
            if (i.eq.25) tot_cld_amt(:)    = vardat(1,1:itime)
            if (i.eq.26) cld_thk(:)        = vardat(1,1:itime)
            if (i.eq.27) cld_hgt(:)        = vardat(1,1:itime)
            if (i.eq.28) lwp(:)            = vardat(1,1:itime)
            if (i.eq.29) qtend(:)          = vardat(1,1:itime)
            if (i.eq.30) advqtend(:)       = vardat(1,1:itime)
            if (i.eq.31) evap(:)           = vardat(1,1:itime)
            if (i.eq.32) ttend(:)          = vardat(1,1:itime)
            if (i.eq.33) advttend(:)       = vardat(1,1:itime)
            if (i.eq.34) radttend(:)       = vardat(1,1:itime)
            if (i.eq.35) latttend(:)       = vardat(1,1:itime)            
            if (i.eq.36) dps_dt_forc(:)=vardat(1,1:itime)
            if (i.eq.37) qs_surface(:)=vardat(1,1:itime)
            if (i.eq.38) s_surface(:)=vardat(1,1:itime)
            if (i.eq.39) wvp(:)=vardat(1,1:itime)
            if (i.eq.40) lwup_sfc(:)=vardat(1,1:itime)
            if (i.eq.41) lwdn_sfc(:)=vardat(1,1:itime)
            if (i.eq.42) swup_sfc(:)=vardat(1,1:itime)
            if (i.eq.43) swdn_sfc(:)=vardat(1,1:itime)
!--- YLL, let tskin_forc = TSKIN
                         tskin_forc(:)=TSKIN
!wfc
! Should the TSkin forcing be the ground temperature? No. Requested to be 29C = TSKIN = 302.15
!            if (i.eq.13) tskin_forc(:)=vardat(1,1:itime) + Tfreeze
!wfc
!   these three var are not in the var file. so
!            if (i.eq.44) tskin_forc(:)=vardat(1,1:itime)
!            if (i.eq.45) alb_24hr(:)=vardat(1,1:itime)
!            if (i.eq.46) alb_day(:)=vardat(1,1:itime)
         enddo
       
         call mpp_close(unit)  
             
         !----record which file was read
         unit= stdlog()
         if (mpp_pe() == 0 ) Write (unit,'(a)') 'used data file: '// surf_name
         call close_file (unit)

       else   ! Read NetCDF data

         if(mpp_pe() == mpp_root_pe() )                              &
           call mpp_error ('scm_twp_mod',                            &
                           'Reading netCDF formatted forcing file.', &
                            NOTE )

         call read_data(forcing_nc,   'lev',     plev_forc(:)          ,no_domain=.true.)
         do t = 1,itime
           call read_data(forcing_nc, 'year',    year_var(t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'month',   month_var(t)          ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'day',     day_var(t)            ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'hour',    hour_var(t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'minute',  minute_var(t)         ,no_domain=.true., timelevel=t)

           call read_data(forcing_nc, 'Temp', t_forc(1:ivert,t)        ,no_domain=.true., timelevel=t)

           call read_data(forcing_nc, 'q', q_forc(1:ivert,t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'u', u_forc(1:ivert,t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'v', v_forc(1:ivert,t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'omega', omega_forc(1:ivert,t)   ,no_domain=.true., timelevel=t)
             
           !i.eq.7 = horizontal temperature advection
           !i.eq.12 = horizontal dry static energy advection
           call read_data(forcing_nc, 'T_adv_h', tadv_hor_forc(:,t)    ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 's_adv_v', tadv_tot_forc(:,t)    ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'q_adv_h', qadv_hor_forc(:,t)    ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'q_adv_v', qadv_tot_forc(:,t)    ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'q1',      q1_obs(:,t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'q2',      q2_obs(:,t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'cld',     arscl(:,t)            ,no_domain=.true., timelevel=t)

!----- read surface data -----  
           call read_data(forcing_nc, 'prec_srf',     precip(t)        ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'LH',           lat_flux_forc(t) ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'SH',           sen_flux_forc(t) ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'p_srf_aver',   ps_forc(t)       ,no_domain=.true., timelevel=t)
!wfc T_srf and T_skin data in file is in degrees C. 
!They are converted to K later.
           call read_data(forcing_nc, 'T_srf',        tair_surf_forc(t),no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'T_skin',       tground(t)       ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'RH_srf',       rh(t)            ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'u_srf',        uwnd(t)          ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'v_srf',        vwnd(t)          ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'rad_net_srf',  netraddnsfc(t)   ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'lw_net_toa',   olr(t)           ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'sw_net_toa',   swup_toa(t)      ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'sw_dn_toa',    swdn_toa(t)      ,no_domain=.true., timelevel=t)
           !convert netdn toa to upward reflected
           swup_toa(t) = swdn_toa(t) - swup_toa(t)
           call read_data(forcing_nc, 'cld_low',      low_cld_amt(t)  ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'cld_mid',      mid_cld_amt(t)  ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'cld_high',     hgh_cld_amt(t)  ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'cld_tot',      tot_cld_amt(t)  ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'cld_thick',    cld_thk(t)      ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'cld_top',      cld_hgt(t)      ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'LWP',          lwp(t)          ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'dh2odt_col',   qtend(t)        ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'h2o_adv_col',  advqtend(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'evap_srf',     evap(t)         ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'dsdt_col',     ttend(t)        ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 's_adv_col',    advttend(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'rad_heat_col', radttend(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'LH_col',       latttend(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'omega_srf',    dps_dt_forc(t)  ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'q_srf',        qs_surface(t)   ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 's_srf',        s_surface(t)    ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'PW',           wvp(t)          ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'lw_up_srf',    lwup_sfc(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'lw_dn_srf',    lwdn_sfc(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'sw_up_srf',    swup_sfc(t)     ,no_domain=.true., timelevel=t)
           call read_data(forcing_nc, 'sw_dn_srf',    swdn_sfc(t)     ,no_domain=.true., timelevel=t)
         enddo
!--- YLL, let tskin_forc = TSKIN
!         tskin_forc(:)=TSKIN
!wfc
! Should the TSkin forcing be the ground temperature? No. Requested to be 29C = TSKIN = 302.15
         tskin_forc(:)=tground(:) +Tfreeze
!wfc
       endif  

        
       !set time variable for VAR data 
       do j = 1, itime
            time_layer(j) = set_date(int(year_var(j)),                 &
                                     int(month_var(j)),                &
                                     int(day_var(j)),                  &
                                     int(hour_var(j)),                 &
                                     int(minute_var(j)),0)             
       enddo        
       time_sfc=time_layer
       Time_var_start = time_layer(1)

       !remove bad values
       !for temperature, fill missing values in theta coordinates then
       !translate back        
       !---- yll, SOMEHOW, t_forc after the following is wrong
       do i = 1, itime 
            where (t_forc(:,i) .gt. 1.) 
                 dum_mat_var (1:ivert,i) = t_forc(:,i) *              &
                                    (1000./plev_forc(:))**(rdgas/cp_air) 
            elsewhere
                 dum_mat_var (1:ivert,i) = 0.
            end where
       end do
       call rem_miss_var(1.,time_layer,dum_mat_var(1:ivert,1:itime))
       do i = 1, itime
            t_forc(:,i) = dum_mat_var(1:ivert,i) /                    &
                                    (1000./plev_forc(:))**(rdgas/cp_air) 
       end do
       !call rem_miss_var(100.,  time_layer,t_forc)
       call rem_miss_var(1.e-06,time_layer,q_forc)
       call rem_miss_var(-9998.,time_layer,u_forc)
       call rem_miss_var(-9998.,time_layer,v_forc)
       call rem_miss_var(-9998.,time_layer,omega_forc)
       call rem_miss_var(-9998.,time_layer,tadv_tot_forc)
       call rem_miss_var(-9998.,time_layer,tadv_hor_forc)
       call rem_miss_var(-9998.,time_layer,qadv_tot_forc)
       call rem_miss_var(-9998.,time_layer,qadv_hor_forc)
       call rem_miss_var(-9998.,time_layer,q1_obs)
       call rem_miss_var(-9998.,time_layer,q2_obs)
       if (nlayer_fields.ge.19) &
           call rem_miss_var(-9998.,time_layer,arscl)
       dum_mat_var(1,1:itime)=precip
       call rem_miss_var(-0.5,time_sfc,dum_mat_var(:,1:itime),1)
       precip=dum_mat_var(1,1:itime)
       dum_mat_var(1,1:itime)=ps_forc
       call rem_miss_var(-9998.,time_sfc,dum_mat_var(:,1:itime),1)
       ps_forc=dum_mat_var(1,1:itime)
       dum_mat_var(1,1:itime)=sen_flux_forc
       call rem_miss_var(-9998.,time_sfc,dum_mat_var(:,1:itime),1)
       sen_flux_forc=dum_mat_var(1,1:itime)
       dum_mat_var(1,1:itime)=lat_flux_forc
       call rem_miss_var(-9998.,time_sfc,dum_mat_var(:,1:itime),1)
       lat_flux_forc=dum_mat_var(1,1:itime)
       dum_mat_var(1,1:itime)=tair_surf_forc
       call rem_miss_var(-9998.,time_sfc,dum_mat_var(:,1:itime),1)
       tair_surf_forc=dum_mat_var(1,1:itime)
    
       if (ALLOCATED(tskin_forc)) then
            dum_mat_var(1,1:itime)=tskin_forc
            call rem_miss_var(-9998.,time_sfc,dum_mat_var(:,1:itime),1)
            tskin_forc=dum_mat_var(1,1:itime)
       end if
   
       if (ALLOCATED(alb_24hr)) then
            dum_mat_var(1,1:itime)=alb_24hr
            call rem_miss_var(-8887.,time_sfc,dum_mat_var(:,1:itime),1)
            alb_24hr=dum_mat_var(1,1:itime)
       end if
   
       if (ALLOCATED(alb_day)) then
            dum_mat_var(1,1:itime)=alb_day
            call rem_miss_var(-8887.,time_sfc,dum_mat_var(:,1:itime),1)
            alb_day=dum_mat_var(1,1:itime)
       end if
       
       ! units change
       ! ------------
       ! (a) convert mixing ratio (g/kg) to specific humidity (kg/kg)
       q_forc(:,:) = q_forc(:,:)/1000.
       q_forc(:,:) = q_forc(:,:)/(1.+q_forc(:,:))

       ! (b) change omega from mb/hr to Pascals/second,
       !     change pressure variables from mb to Pascals
       omega_forc(:,:) =  100. * omega_forc(:,:)/3600.
       plev_forc(:)    =  100. * plev_forc(:)
       ps_forc(:)      =  100. * ps_forc(:)

       ! (c) convert temperature advections K/hour to K/sec
       !     add vertical s-advection to TOTAL T ADV
       tadv_hor_forc(:,:) = tadv_hor_forc(:,:)/3600.
       tadv_tot_forc(:,:) =(tadv_tot_forc(:,:)/3600.)+tadv_hor_forc(:,:)
        
       ! (d) convert Qadv (g/kg/hour) to (kg/kg/sec)
       !     combine horizontal and vertical advection of Q
       qadv_hor_forc(:,:) = qadv_hor_forc(:,:)/3600./1000.
       qadv_tot_forc(:,:) =(qadv_tot_forc(:,:)/3600./1000.) +          &
                            qadv_hor_forc(:,:)
        
       ! (e) convert mixing ratio advection to specific humidity adv.
       qadv_hor_forc(:,:) = qadv_hor_forc(:,:) *                       &
                            (1.-q_forc(:,:))*(1.-q_forc(:,:))
       qadv_tot_forc(:,:) = qadv_tot_forc(:,:) *                       &
                            (1.-q_forc(:,:))*(1.-q_forc(:,:))
                            
       ! (f) convert surface air and ground temperature from C to K
       tair_surf_forc(:)=tair_surf_forc(:)+Tfreeze
       tground(:) = tground(:) + Tfreeze

       ! (g) convert q1 K/hr to K/day, and q2 K/hr to g/kg/day
       q1_obs(:,:) = q1_obs(:,:)*24.
!       q2_obs(:,:) = q2_obs(:,:)*24.*1000.*cp_air/hlv
       ! Convert q2 K/Hr to kg/kg/s for forcing 
       q2_obs(:,:) = q2_obs(:,:)*cp_air/hlv/3600.
       ! (h) convert rh (%) to rh (fraction)
       rh(:) = rh(:)/100.

       ! (i) convert cm lwp to kg/m2
       lwp(:) = lwp(:)*10.

       ! (j) convert mm/hr to kg/m2/sec
       qtend(:) = qtend(:) / 3600.
       advqtend(:) = advqtend(:) / 3600.        
       evap(:) = evap(:) / 3600.
!       precip(:) = precip(:) /3600.
       ! For a1 forcing with precipitation, precip should be mm/day!!
       precip(:) = precip(:) *24.

       ! (k) handle dps_dt_forc
       !calculate it for those cases for which omega-surface is not 
       !available
       if (use_given_omega_surface .and. nfields .ge. 36) then
            !change units from mb/hr to Pa/s
            dps_dt_forc(:) = dps_dt_forc(:)*100./3600.       
       else
            dps_dt_forc(1)=-9999.
            do i=2,itime-1
                 call get_time((time_sfc(i+1) - time_sfc(i-1)),seconds,&
                      days )
                 dps_dt_forc(i) = (ps_forc(i+1)-ps_forc(i-1)) /        &
                      REAL(seconds+days*86400)
            enddo
            dps_dt_forc(itime)=-9999.        
        
            !remove bad values
            dum_mat_var(1,1:itime)=dps_dt_forc
            call rem_miss_var(-9998.,time_sfc,dum_mat_var(:,1:itime),1)
            dps_dt_forc=dum_mat_var(1,1:itime)
       endif
      
       ! (l) change units on wvp from cm to kg/m2
       if (allocated(wvp)) wvp(:) = wvp(:) * 10.

       ! (m) compute tau_inv_forc
       tau_inv_forc(:,:) = 2.*SQRT(u_forc(:,:)*u_forc(:,:) +           &
            v_forc(:,:)*v_forc(:,:))/(370.*1000.)

       
       ! (n) handle surface albedo
       !
       !     If provided use that, if not calculate it here.
       !
       !     If not, compute surface albedo and its time mean
       !     note that the time evolving surface albedo is calculated 
       !     from obserations only when the swdn_sfc > 50. W/m2 to 
       !     prevent ridiculous values being set at night.  At night
       !     values are interpolated using the rem_miss_var subroutine.
       !     

       if (allocated(alb_24hr)) then
             
             albedo_sfc(:) = alb_24hr(:)
             tmp1 = 0.
             tmp2 = 0.
             do i = 1, itime
                if (swdn_sfc(i) .gt. 50.) then
                    tmp1 = tmp1 + alb_24hr(i)
                    tmp2 = tmp2 + 1.
                end if    
             enddo
             avg_albedo_sfc = tmp1/tmp2
       
       else            
       
       if (allocated(swup_sfc).and.allocated(swdn_sfc)) then
             tmp1 = 0.
             tmp2 = 0.
             albedo_sfc(:) = -0.9
             do i = 1, itime
                    tmp1 = tmp1 + swup_sfc(i)
                    tmp2 = tmp2 + swdn_sfc(i)
                    if (swdn_sfc(i) .gt. 50.) then
                         albedo_sfc(i) = swup_sfc(i)/swdn_sfc(i)
                    end if
             enddo
             avg_albedo_sfc = tmp1/tmp2
        
             dum_mat_var(1,1:itime)=albedo_sfc
             call rem_miss_var(-0.89,time_sfc,dum_mat_var(:,1:itime),1)
             albedo_sfc = dum_mat_var(1,1:itime)            
       else
             albedo_sfc = 0.19
             avg_albedo_sfc = 0.19
       end if
       
       end if
       
       !(o) make sure rough_heat = rough_mom / 7.4
       rough_heat = rough_mom / 7.4

       !(p) change units on qs_surface to specific humidity
       if (allocated(qs_surface)) qs_surface = qs_surface / (1.+qs_surface)
       
       !(q) change s_surface from K to J/kg
       if (allocated(s_surface)) s_surface(:) = s_surface(:) * cp_air

       !(r) change tskin from degC to degK
       !if (ALLOCATED(tskin_forc)) tskin_forc = tskin_forc + Tfreeze


       !(s) check balances
       coltcal(:) = 0.
       colqcal(:) = 0.
       
       !compute model level thickness in Pascals
       levthick =  abs( plev_forc(2) - plev_forc(1) )
       
       fname='budget.out'
       call mpp_open(unit,fname)
       write (unit,*) ' '
       write (unit,*) ' i   coltcal  advttend   tomegas   colqcal  advqtend   qomegas'        
       write (unit,*) '--------------------------------------------------------------'
       

       !There is no time dependence on the forcing plevs so 
       !compute last level which is above 100 hPa once.
       jend = 1
       do while(plev_forc(jend).gt.10000.)
          jend=jend+1
          if(jend .gt. size(plev_forc) ) exit
       enddo
       jend=jend-1
            
       !add in surface level
       do i=1,itime
       
            !compute starting level
            jstart =1
            do while (ps_forc(i).lt.plev_forc(jstart))
                  jstart = jstart + 1
            enddo
            
            !do surface level
            coltcal(i) = coltcal(i) + tadv_tot_forc(jstart,i) * &
                     cp_air*(0.5*levthick+(ps_forc(i)-plev_forc(jstart)))/grav
            colqcal(i) = colqcal(i) + qadv_tot_forc(jstart,i) * &
                        (0.5*levthick+(ps_forc(i)-plev_forc(jstart)))/grav
            
            !do main levels
            do j = jstart+1,jend-1
            coltcal(i) = coltcal(i) + tadv_tot_forc(j,i) * &
                     cp_air*levthick/grav
            colqcal(i) = colqcal(i) + qadv_tot_forc(j,i) * &
                        levthick/grav
            enddo

            !do top level
            coltcal(i) = coltcal(i) + tadv_tot_forc(jend,i) * &
                     cp_air*(0.5*levthick+plev_forc(jend)-10000.)/grav
            colqcal(i) = colqcal(i) + qadv_tot_forc(jend,i) * &
                        (0.5*levthick+plev_forc(jend)-10000.)/grav
            
            !calculate 
            if (allocated(s_surface).and.allocated(qs_surface)) then       
            write(unit,'(i3,6(1X,f9.3))') i,coltcal(i),advttend(i),&
                    dps_dt_forc(i)*s_surface(i)/grav,colqcal(i)*86400.,&
                    advqtend(i)*86400.,&
                    86400.*dps_dt_forc(i)*qs_surface(i)/grav
            else
            write(unit,'(i3,6(1X,E9.3))') i,coltcal(i),advttend(i),&
                    0.,colqcal(i)*86400.,advqtend(i)*86400.,0.
            end if

       end do
       
       write (unit,*) ' '
       write (unit,*) 'NOTE: coltcal will not equal advttend because: '
       write (unit,*) ' '
       write (unit,*) '      coltcal = column integrated temperature advection in W/m2'
       write (unit,*) '     advttend = column integrated dry static energy advection in W/m2'

       call mpp_close(unit)  
        
       !(t) change units on ARSCL cloud fraction from % to fraction       
       if (nlayer_fields.ge.19) arscl = arscl / 100.

       ! -1/Lg * integral (Q2*dp, p=0->Pg) + Fq = 1/g*integral (QR*dp, p=0->Pg)
       ! Donner et al 2001 Eq. 2
       ! QR < 0 => convection
       ! => 1/Lg * integral (Q2*dp, p=0->Pg) - Fq > 0 => convection

       ! Calculate integral of Q2*dp
       ! plev_forc is mid level (pfull) pressure.
       ! Need to get limits. 

       ! plev_forc has no time dependency. Calculate boundary pressures.
       ! NB phalf_forc(1) will be time dependent through ps_forc but ps_forc
       ! is constant in the data set in use.
       ! NB Move this inside the time loop if ps_forc varies with time.
       ! plev_forc is on pfull levels in GCM terms. 
       ! March up from there.
       kmax = size(q2_obs,1)
       
       fname='a1_limit.out'
       call mpp_open(unit,fname)
       write (unit,*) ' '
       write (unit,*) ' i     colq2_obs       SurfaceFlux*grav/HLV   SF + colq2_obs   '
       write (unit,*) '---------------------------------------------------------------'
       do i = 1, itime

          jstart = 1
          do while (ps_forc(i).lt.plev_forc(jstart))
            jstart = jstart + 1
          enddo

          !compute starting level
          do k=jstart,kmax
            phalf_forc(k) = (plev_forc(k) + plev_forc(k-1))/2.0
          enddo


          q2_integral(i) = q2_obs(jstart,i) * &
                            (ps_forc(i) - plev_forc(jstart))  
          
          !March down through the atmosphere

          do k = kmax,jstart+1,-1
!            q2_integral(i) = q2_integral(i) + q2_obs(k,i) * &
!                            (phalf_forc(k) - phalf_forc(k+1))  
            tmp1 = q2_obs(k,i) * &
                            (phalf_forc(k) - phalf_forc(k+1))  
            q2_integral(i) = q2_integral(i) + tmp1
          enddo
          gFq(i) = lat_flux_forc(i) ! W/m^2

!We want to end up with 
! a1= (1+r)^2*[g*Fq - 1/HLV*integral(0:Pg)(Q2*dp)]/integral(0:Pg)(~Qr*dp)
! Q2 above is actually Q2(GATE)*cp/L so need to multiply by HLV below.
! Actually units of ~Qr*dp in donner_deep is kg(H20)/kg/s*Pa 
! Therefore need to divide the equation above by HLV, so q2_integral is left 
! alone and divide Fq by HLV
          gFq(i) = gFq(i) * grav / HLV
!          q2_integral(i) = q2_integral(i) * HLV
          write(unit,'(i3,3E20.12)') i,q2_integral(i), gFq(i), &
                                  gFq(i) + q2_integral(i)
        enddo
        call mpp_close(unit)
       
!-----------------------------------------------------------------------
! 
!      HANDLE RESTART INFO
!
    ! If restart file exists override values with the restart values

    if( file_exist(trim(fname_nc)) ) then

          if(mpp_pe() == mpp_root_pe() )                            &
          call mpp_error ('scm_twp_mod',                            &
                          'Reading netCDF formatted restart file.', &
                           NOTE )

          call read_data(fname_nc, 'SENFLUX',    SENFLUX   )
          call read_data(fname_nc, 'EVAPFLUX',   EVAPFLUX  )
          call read_data(fname_nc, 'TSKIN',      TSKIN     )
          call read_data(fname_nc, 'ALBEDO_OBS', ALBEDO_OBS)

    else if (file_exist('INPUT/twp_forc.res')) then

          unit=open_restart_file('INPUT/twp_forc.res',action='read')
          read(unit)  SENFLUX, EVAPFLUX, TSKIN, ALBEDO_OBS
          call close_file(unit)
     
    endif
! 
!-----------------------------------------------------------------------
       if (mpp_pe() == mpp_root_pe()) then
         print*,'t_forc_at16 = ', t_forc(:,16)
         print*,'t_forc_at17 = ', t_forc(:,17)
       endif
end subroutine twp_data_read


!#######################################################################
!#######################################################################

subroutine twp_forc_init(time_interp, pdamp, As, Bs)
#include "fv_arrays.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a given time and initializes the model 
!      variables for that time. This involves interpolating the global 
!      storage fields to a given time and pressure, to do the resulting
!      calculation.
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
!
!      -------------------
!      INTERNAL VARIABLES:
!      -------------------
!
!
!      k               counting integers
!      itime_more      indicates indice of TIME_TYPE for which 
!                      TIME_TYPE(itime_more)>= time_interp
!      itime_less      indicates indice of TIME_TYPE for which 
!                      TIME_TYPE(itime_more)<= time_interp
!      KDIM            no. of vertical levels to pt array
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
!      pi_fac          (p/pref(1000mb))**(Rd/cp) = T/theta
!
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

TYPE(TIME_TYPE)                          :: time_interp
REAL,  INTENT (IN)                       :: pdamp
REAL,  INTENT (IN)   , DIMENSION(:)      :: As,Bs
!  Internal variables
!  ------------------

INTEGER  :: itime_less, itime_more, KDIM, k, jstart, jstart_stand
INTEGER  :: nlev_interp, days, months, years,seconds, minutes, hours
REAL     :: weight_less, weight_more, summerweight

REAL, DIMENSION(ivert+2+ivert_stand)                :: plev_int
REAL, DIMENSION(ivert+2+ivert_stand)                :: field_int_less, &
                                                       field_int_more
real, dimension(size(q,3))                         :: tmp_p,tmp_ans
real, dimension(size(q,1),size(q,2))               :: ps_init
real, dimension(size(q,1),size(q,2),size(q,3))     :: pfull
real, dimension(size(q,1),size(q,2),size(q,3)+1)   :: phalf

!
!CODE
!
integer :: i,j
#include "fv_point.inc"
       nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
       nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
       nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
       nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
       nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   !h1g
       nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )   !h1g
       nq2_clo = get_tracer_index ( MODEL_ATMOS, 'a1q2_closure' )
       nprec_clo = get_tracer_index ( MODEL_ATMOS, 'a1prec_closure' )
       
       ! --- find out # of vertical levels
       KDIM = SIZE(pt,3)


! --- CREATE PFULL AND PHALF AND PI_FAC  ----------------------------- !

       ! --- find out time in sfc pressure field
       call interp_time(time_interp,time_sfc,itime_less,itime_more,    &
                        weight_less,weight_more)
 
       ! --- interpolate surface pressure field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=ps_forc(itime_less)
       field_int_more(:)=ps_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       ps_init(:,:) = tmp_ans(1)
       ps = ps_init
       

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
           call get_eta_level(nlev, ps(i,j) , pfull(i,j,:), phalf(i,j,:))
         enddo
       enddo

! --- compute starting level for variational dataset ----------------- !

      jstart =1
      do while (ps(1,1).lt.plev_forc(jstart))
          jstart = jstart + 1
      enddo
           
! --- compute starting level for midlatitude reference profile  ------ !

      jstart_stand = 1
      do while (plev_forc(ivert).lt.plev_stand(jstart_stand,1))
          jstart_stand = jstart_stand + 1
      enddo       
      
! --- create average mid-latitude profile based upon month ----------- !       

      call get_date(time_interp,years,months,days,hours,minutes,seconds)
      summerweight = 1. - ((real(abs(months - 7)))/6.)
      print *,'summerweight = ', summerweight
      summerweight = 0.0    !-------- only one tropic data available

! --- CREATE INITIAL pt  ---------------------------------------------- !

       ! --- find nearest time to time_interp in LAYERED fields
       call interp_time(time_interp,time_layer,itime_less,itime_more,  &
                        weight_less,weight_more)
       
       ! --- create field data for interpolator
       ! --- NOTES:
       ! ---
       ! --- level 1 corresponds to the surface and it is assumed that &
       ! ---     theta(surface) = theta(BOT of Layered fields)
            
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=t_forc(jstart,itime_less)*                    &
            ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
       field_int_more(1)=t_forc(jstart,itime_more)*                    &
            ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
       plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
       field_int_less(2:(ivert-jstart+2))=                             &
                                   t_forc(jstart:ivert,itime_less)
       field_int_more(2:(ivert-jstart+2))=                             &
                                   t_forc(jstart:ivert,itime_more)
       
       plev_int((ivert-jstart+3):                                      &
                (ivert-jstart+2+ivert_stand-jstart_stand+1))=                    &
            plev_stand(jstart_stand:ivert_stand,1)
       field_int_less((ivert-jstart+3):                                &
                (ivert-jstart+2+ivert_stand-jstart_stand+1))=          &
            ((1.-summerweight)*t_stand(jstart_stand:ivert_stand,1)+    &
                 summerweight *t_stand(jstart_stand:ivert_stand,2))
       field_int_more((ivert-jstart+3):                                &
                (ivert-jstart+2+ivert_stand-jstart_stand+1))=          &
            ((1.-summerweight)*t_stand(jstart_stand:ivert_stand,1)+    &
                 summerweight *t_stand(jstart_stand:ivert_stand,2))
       nlev_interp = ivert-jstart+2+ivert_stand-jstart_stand+1
       tmp_p(:)=pfull(1,1,:)
        
       ! --- do interpolation

        if (mpp_pe() == mpp_root_pe()) then
         print*, 'itime_less = ', itime_less
         print*, 'itime_more = ', itime_more
         print*,'tmp_p = ', tmp_p
         print*,'plev_int = ', plev_int
         print*,'field_int_less = ', field_int_less
         print*,'field_int_more = ', field_int_more
         print*,'t_forc_less = ', t_forc(:,itime_less)
         print*,'t_forc_more = ', t_forc(:,itime_more)
         print*,'===================================================================================='
        endif 
         
         
         
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)

       ! --- assign pt
       do k = 1, KDIM
            pt(:,:,k)=tmp_ans(k)
       enddo
           
! --- CREATE INITIAL Qv  --------------------------------------------- !

 
       ! --- find nearest time to time_interp in Layered fields ------ !
       ! --- no need since this was done for T 
      
       ! --- create field data for interpolator
       ! --- NOTE:
       ! ---
       ! --- level 1 corresponds to the surface and it is assumed that
       !          qv(surface) = qv(BOT)
           
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=q_forc(jstart,itime_less)
       field_int_more(1)=q_forc(jstart,itime_more)
       plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
       field_int_less(2:(ivert-jstart+2))=                             &
                                   q_forc(jstart:ivert,itime_less)
       field_int_more(2:(ivert-jstart+2))=                             &
                                   q_forc(jstart:ivert,itime_more)
       
       plev_int((ivert-jstart+3):                                      &
                (ivert-jstart+2+ivert_stand-jstart_stand+1))=                    &
            plev_stand(jstart_stand:ivert_stand,1)
       field_int_less((ivert-jstart+3):                                &
                (ivert-jstart+2+ivert_stand-jstart_stand+1))=          &
            ((1.-summerweight)*q_stand(jstart_stand:ivert_stand,1)+    &
                 summerweight *q_stand(jstart_stand:ivert_stand,2))
       field_int_more((ivert-jstart+3):                                &
                (ivert-jstart+2+ivert_stand-jstart_stand+1))=          &
            ((1.-summerweight)*q_stand(jstart_stand:ivert_stand,1)+    &
                 summerweight *q_stand(jstart_stand:ivert_stand,2))
       nlev_interp = ivert-jstart+2+ivert_stand-jstart_stand+1
       tmp_p(:)=pfull(1,1,:)

       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
      
       ! --- assign answer to qv
       do k = 1, KDIM
            q(:,:,k,nsphum)=tmp_ans(k)
       enddo
       
! --- CREATE INITIAL U,V  -------------------------------------------- !
      
       ! --- find nearest time to time_interp in LAYERED FIELDS
       ! --- no need since this was done for T 
      
       ! --- create field data for interpolator
       ! --- NOTE: 
       ! --- 
       ! --- level 1 corresponds to the surface and it is assumed that
       !           u,v(surface) = u,v(BOT)
             
      
       tmp_ans = get_column(u_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                            weight_less,weight_more, pfull(1,1,:))
       ! --- assign answer to u
       do k = 1, KDIM
            ua(:,:,k)=tmp_ans(k)
       enddo
       u_srf(:,:)=ua(:,:,KDIM)

       ! --- create field data for interpolator
                
       tmp_ans = get_column(v_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                            weight_less,weight_more, pfull(1,1,:))
       ! --- assign answer to v
       do k = 1, KDIM
            va(:,:,k)=tmp_ans(k)
       enddo
       ! Initialize surface winds
       v_srf(:,:)=va(:,:,KDIM)

! --- YLL use fixed SST 29 oC
! --- compute ground skin temperature -------------------------------- !
!       TSKIN = 302.15 
! --- YLL no surface flux iteration being called
! --- interpolate sensible and latent heat flux and surface albedo  -- !
     
       !---- interpolate closure field
       if (nq2_clo > 0 ) then
         tmp_p(:)=50000.
         plev_int(:)=p00
         field_int_less(:)= -1.0 * gFq(itime_less) - q2_integral(itime_less)
         field_int_more(:)= -1.0 * gFq(itime_more) - q2_integral(itime_more)
         call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
              field_int_less,field_int_more,tmp_ans)
         q(:,:,:,nq2_clo) = tmp_ans(1)
       endif

       if (nprec_clo > 0 ) then
         field_int_less(:)=precip(itime_less)
         field_int_more(:)=precip(itime_more)
         call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                  field_int_less,field_int_more,tmp_ans)
         q(:,:,:,nprec_clo) = tmp_ans(1)* a1_precfrac
       endif

       !---- interpolate sensible heat flux field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=sen_flux_forc(itime_less)
       field_int_more(:)=sen_flux_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       SENFLUX = tmp_ans(1) 

       !---- interpolate latent heat flux field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=lat_flux_forc(itime_less)
       field_int_more(:)=lat_flux_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       EVAPFLUX = tmp_ans(1)/hlv
             
       !---- interpolate surface albedo field
       if (do_iop_mean_albedo) then
            ALBEDO_OBS = avg_albedo_sfc(1)
       else
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=albedo_sfc(itime_less)
            field_int_more(:)=albedo_sfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,        &
                 plev_int,field_int_less,field_int_more,tmp_ans)
            ALBEDO_OBS = tmp_ans(1) 
      end if
      
! --- CREATE INITIAL ql,qi,qa ---------------------------------------- !

       q(:,:,:,nql)=0.0
       q(:,:,:,nqa)=0.0
       q(:,:,:,nqi)=0.0
       if ( nqn  > 0 ) q(:,:,:,nqn)  = 0.0
       if ( nqni > 0 ) q(:,:,:,nqni) = 0.0

!------ YLL print the init T,qv,U,V profile and compare with input data

    if (mpp_pe() == mpp_root_pe()) then
     print *, 't_stand = ', t_stand(:,1)
     print *, 't_forc  = ', t_forc(:,itime_less)
     print *, 'q_forc  = ', q_forc(:,1)
     print *, 'u_forc  = ', u_forc(:,1)
     print *, 'phalf = ', phalf
     print *,'init T profile = ', pt(1,1,:)
     print *,'init q profile = ', q(1,1,:,nsphum)
     print *,'init U profile = ', ua(1,1,:)
     print *,'init V profile = ', va(1,1,:)
     print *, 'ALBEDO = ', ALBEDO_OBS
     print *, 'SENFLUX = ', SENFLUX
     print *, 'EVAPFLUX = ', EVAPFLUX
    endif
! 
!-----------------------------------------------------------------------
 
 end subroutine twp_forc_init

!#######################################################################
!#######################################################################


subroutine twp_forc_end ()

  integer :: unit
  character(len=64) :: fname_nc='RESTART/twp_forc.res.nc'
  
  if (.NOT.twp_forc_initialized) return
  twp_forc_initialized = .FALSE.

  DEALLOCATE (  plev_stand, t_stand, q_stand, O3_stand, plev_forc,     &
                t_forc, q_forc, u_forc, v_forc, omega_forc,            &
                tadv_tot_forc, tadv_hor_forc, qadv_tot_forc,           &
                qadv_hor_forc, tau_inv_forc, q1_obs, q2_obs,           &
                ps_forc, dps_dt_forc, sen_flux_forc, lat_flux_forc,    &
                tair_surf_forc, precip, tground, rh, uwnd, vwnd,       &
                netraddnsfc, olr, swup_toa, swdn_toa, low_cld_amt,     &
                mid_cld_amt, hgh_cld_amt, tot_cld_amt, cld_thk,        &
                cld_hgt, lwp, qtend, advqtend, evap, ttend, advttend,  &
                radttend, latttend, albedo_sfc, avg_albedo_sfc )
                
  if (allocated(s_surface)) deallocate(s_surface)
  if (allocated(qs_surface)) deallocate(qs_surface)
  if (allocated(wvp)) deallocate(wvp)
  if (allocated(swdn_sfc)) deallocate(swdn_sfc)
  if (allocated(swup_sfc)) deallocate(swup_sfc)
  if (allocated(lwdn_sfc)) deallocate(lwdn_sfc)
  if (allocated(lwup_sfc)) deallocate(lwup_sfc)
  if (allocated(tskin_forc)) deallocate(tskin_forc)
  if (allocated(alb_24hr)) deallocate(alb_24hr)
  if (allocated(alb_day)) deallocate(alb_day)
  if (allocated(arscl)) deallocate(arscl)
  
! write to the restart file

  if( do_netcdf_restart ) then

    if (mpp_pe() == mpp_root_pe()) then
        call write_data (fname_nc, 'SENFLUX',    SENFLUX   )
        call write_data (fname_nc, 'EVAPFLUX',   EVAPFLUX  )
        call write_data (fname_nc, 'TSKIN',      TSKIN     )
        call write_data (fname_nc, 'ALBEDO_OBS', ALBEDO_OBS)
    endif

  else

       unit=open_restart_file('RESTART/twp_forc.res',action='write')
       if(mpp_pe() == mpp_root_pe())  write(unit) SENFLUX, EVAPFLUX, TSKIN, &
                                             ALBEDO_OBS
       call close_file(unit)

  endif

! 
!-----------------------------------------------------------------------
 

end subroutine twp_forc_end


!#######################################################################
!#######################################################################

subroutine twp_forc_diagnostic_init(axes, Time)

implicit none

integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time


   
! --- initialize axes -------------------------------------------------!

       id_temp_obs = register_diag_field (mod_name, 'temp_obs',        &
          axes(1:3), Time,                                             &
          'Observed Temperature', 'K',                                 &
          missing_value = missing_value )
       id_qv_obs   = register_diag_field (mod_name, 'qv_obs',          &
          axes(1:3), Time,                                             &
          'Observed Water Vapor Specific Humidity', 'g/kg',            &
         missing_value = missing_value )
       id_uwnd_obs   = register_diag_field (mod_name, 'uwnd_obs',      &
          axes(1:3), Time, 'Observed Zonal Wind', 'm/s',               &
         missing_value = missing_value )
       id_vwnd_obs   = register_diag_field (mod_name, 'vwnd_obs',      &
          axes(1:3), Time, 'Observed Meridional Wind', 'm/s',          &
         missing_value = missing_value )
       id_omega = register_diag_field (mod_name, 'omega',              &
          axes(1:3),  Time,                                            &
          'Observed Pressure Velocity', 'hPa/day',                     &
          missing_value = missing_value )
       id_q1_obs = register_diag_field (mod_name, 'q1_obs',            &
          axes(1:3), Time,                                             &
          'Apparent Heat Source', 'K/day',                             &
          missing_value = missing_value )
       id_q2_obs = register_diag_field (mod_name, 'q2_obs',            &
          axes(1:3), Time,                                             &
          'Apparent Moisture Sink', 'g/kg/day',                        &
          missing_value = missing_value )

       if (nlayer_fields .ge. 19) then       
       id_arscl = register_diag_field (mod_name, 'arscl_obs',          &
          axes(1:3), Time,                                             &
          'ARSCL Cloud Fraction', 'fraction',                          &
          missing_value = missing_value )
       end if
              
       id_tair_obs = register_diag_field (mod_name, 'tair_obs',        &
          axes(1:2), Time,                                             &
         'Observed Reference Height Temperature', 'K',                 &
         missing_value = missing_value )
       id_tground_obs = register_diag_field (mod_name, 'tground_obs',  &
          axes(1:2), Time,                                             &
         'Observed Ground Temperature', 'K',                           &
         missing_value = missing_value )
       id_rh_obs = register_diag_field (mod_name, 'rh_obs',            &
          axes(1:2), Time,                                             &
         'Observed Reference Height Relative Humidity', 'percent',     &
         missing_value = missing_value )
       id_uref_obs = register_diag_field (mod_name, 'uref_obs',        &
          axes(1:2), Time,                                             &
         'Observed Reference Height Zonal Wind', 'm/s',                &
         missing_value = missing_value )
       id_vref_obs = register_diag_field (mod_name, 'vref_obs',        &
          axes(1:2), Time,                                             &
         'Observed Reference Height Meridional Wind', 'm/s',           &
         missing_value = missing_value )
       id_precip_obs = register_diag_field (mod_name,                  &
         'precip_obs',  axes(1:2), Time,                               &
         'Observed Precipitation', 'kg/m2/sec',                        &
         missing_value = missing_value )
       id_latent_obs = register_diag_field (mod_name,                  &
         'latent_obs',  axes(1:2), Time,                               &
         'Observed Surface Latent Heat Flux', 'W/m2',                  &
         missing_value = missing_value )
       id_sens_obs = register_diag_field (mod_name, 'sens_obs',        &
          axes(1:2), Time,                                             &
         'Observed Surface Sensible Heat Flux', 'W/m2',                &
         missing_value = missing_value )
       id_netdnsfc_obs = register_diag_field (mod_name,                &
         'netdnsfc_obs',                                               &
         axes(1:2), Time,                                              &
         'Observed Net Downward Surface Radiation', 'W/m2',            &
         missing_value = missing_value )
       id_olr_obs = register_diag_field (mod_name, 'olr_obs',          &
          axes(1:2), Time,                                             &
         'Observed Outgoing Longwave Radiation', 'W/m2',               &
         missing_value = missing_value )
       id_swup_toa_obs = register_diag_field (mod_name,                &
         'swup_toa_obs',  axes(1:2), Time,                             &
         'Observed Outgoing Shortwave Radiation', 'W/m2',              &
         missing_value = missing_value )
       id_swdn_toa_obs = register_diag_field (mod_name,                &
         'swdn_toa_obs',  axes(1:2), Time,                             &
         'Observed Incoming Shortwave Radiation', 'W/m2',              &
         missing_value = missing_value )
       id_tot_cld_amt_obs = register_diag_field (mod_name,             &
         'tot_cld_amt_obs',  axes(1:2), Time,                          &
         'Observed Total Cloud Amount', 'percent',                     &
         missing_value = missing_value )
       id_lwp_obs = register_diag_field (mod_name, 'lwp_obs',          &
          axes(1:2), Time,                                             &
          'Observed Liquid Water Path', 'kg/m2',                       &
          missing_value = missing_value )   
       id_low_cld_amt_obs = register_diag_field (mod_name,             &
         'low_cld_amt_obs',  axes(1:2), Time,                          &
         'Observed Low Cloud Amount', 'percent',                       &
         missing_value = missing_value )
       id_mid_cld_amt_obs = register_diag_field (mod_name,             &
         'mid_cld_amt_obs',  axes(1:2), Time,                          &
         'Observed Middle Cloud Amount', 'percent',                    &
         missing_value = missing_value )
       id_hgh_cld_amt_obs = register_diag_field (mod_name,             &
         'hgh_cld_amt_obs',  axes(1:2), Time,                          &
         'Observed High Cloud Amount', 'percent',                      &
         missing_value = missing_value )
       id_cld_thk_obs = register_diag_field (mod_name,                 &
         'cld_thk_obs',  axes(1:2), Time,                              &
         'Observed Cloud Thickness', 'km',                             &
         missing_value = missing_value )
       id_cld_hgt_obs = register_diag_field (mod_name,                 &
         'cld_hgt_obs',  axes(1:2), Time,                              &
         'Observed Cloud Height', 'km',                                &
         missing_value = missing_value )
       id_qtend_obs = register_diag_field (mod_name, 'qtend_obs',      &
          axes(1:2), Time,                                             &
         'Observed Column Moisture Tendency', 'kg/m2/sec',             &
         missing_value = missing_value )
       id_advqtend_obs = register_diag_field (mod_name,                &
         'advqtend_obs',  axes(1:2), Time,                             &
         'Observed Column Moisture Advection Tendency', 'kg/m2/sec',   &
         missing_value = missing_value )
       id_evap_obs = register_diag_field (mod_name, 'evap_obs',        &
          axes(1:2), Time,                                             &
          'Observed Surface Evaporation', 'kg/m2/sec',                 &
          missing_value = missing_value )
       id_ttend_obs = register_diag_field (mod_name, 'ttend_obs',      &
          axes(1:2), Time,                                             &
          'Observed Column Energy tendency', 'W/m2',                   &
          missing_value = missing_value )
       id_advttend_obs = register_diag_field (mod_name,                &
         'advttend_obs',  axes(1:2), Time,                             &
         'Observed Column Energy Advection Tendency', 'W/m2',          &
         missing_value = missing_value )
       id_radttend_obs = register_diag_field (mod_name,                &
         'radttend_obs',  axes(1:2), Time,                             &
         'Observed Column Energy Radiation Tendency', 'W/m2',          &
         missing_value = missing_value )
       id_latttend_obs = register_diag_field (mod_name,                &
         'latttend_obs',  axes(1:2), Time,                             &
         'Observed Column Energy Latent Heating Tendency', 'W/m2',     &
         missing_value = missing_value )       
       

       id_Tadv_hor_forc = register_diag_field (mod_name, 'dT_hor_forc', &
          axes(1:3), Time,                                             &
          'Temperature tendency from SCM horizontal advection', 'K/sec', &
          missing_value = missing_value )
       id_Qadv_hor_forc = register_diag_field (mod_name, 'dQ_hor_forc', &
          axes(1:3), Time,                                             &
          'Specific humidity tendency from SCM horizontal advection', 'kg/kg/sec', &
          missing_value = missing_value )
       id_vadvT = register_diag_field (mod_name, 'dT_vadv',            &
          axes(1:3), Time,                                             &
          'Temperature tendency from SCM vertical advection', 'K/sec', &
          missing_value = missing_value )
       id_adiT = register_diag_field (mod_name, 'dT_adi',            &
          axes(1:3), Time,                                             &
          'Temperature tendency from SCM vertical advection2', 'K/sec', &
          missing_value = missing_value )
       id_vadvQ = register_diag_field (mod_name, 'dq_vadv',            &
          axes(1:3), Time,                                             &
          'Moisture tendency from SCM vertical advection', 'K/sec', &
          missing_value = missing_value )

       id_nudg_t = register_diag_field (mod_name, 'nudg_t',            &
          axes(1:3), Time,                                             &
          'Temperature tendency from Nudging', 'K/sec',                &
          missing_value = missing_value )
       id_nudg_q = register_diag_field (mod_name, 'nudg_q',            &
          axes(1:3), Time,                                             &
          'Water vapor tendency from Nudging', 'kg/kg/sec',            &
          missing_value = missing_value )
       
       id_nudg_u = register_diag_field (mod_name, 'nudg_u',            &
          axes(1:3), Time,                                             &
          'U tendency from Nudging', 'm/s2',                &
          missing_value = missing_value )
       id_nudg_v = register_diag_field (mod_name, 'nudg_v',            &
          axes(1:3), Time,                                             &
          'V tendency from Nudging', 'm/s2',            &
          missing_value = missing_value )

id_tdt_lf = register_diag_field (mod_name, 'tdt_lf', axes(1:3), Time, &
     'Temperature tendencies due to large-scale forcing', 'K/s', missing_value = missing_value)

id_qvdt_lf = register_diag_field (mod_name, 'qvdt_lf', axes(1:3), Time, &
     'Vapor tendencies due to large-scale forcing', 'kg/kg/s', missing_value = missing_value)

       if (nfields .ge. 43) then
       id_lwup_sfc_obs = register_diag_field (mod_name,                &
         'lwup_sfc_obs',  axes(1:2), Time,                             &
         'Observed Surface Upward Longwave Radiation', 'W/m2',         &
         missing_value = missing_value )
       id_lwdn_sfc_obs = register_diag_field (mod_name,                &
         'lwdn_sfc_obs',  axes(1:2), Time,                             &
         'Observed Surface Downward Longwave Radiation', 'W/m2',       &
         missing_value = missing_value )
       id_swup_sfc_obs = register_diag_field (mod_name,                &
         'swup_sfc_obs',  axes(1:2), Time,                             &
         'Observed Surface Upward Shortwave Radiation', 'W/m2',        &
         missing_value = missing_value )
       id_swdn_sfc_obs = register_diag_field (mod_name,                &
         'swdn_sfc_obs',  axes(1:2), Time,                             &
         'Observed Surface Downward Shortwave Radiation', 'W/m2',      &
         missing_value = missing_value )
       id_pw_obs = register_diag_field (mod_name, 'pw_obs',            &
          axes(1:2), Time,                                             &
          'Observed Precipitable Water', 'kg/m2',                      &
           missing_value = missing_value )
       id_albedo_sfc_obs = register_diag_field (mod_name,              &
         'albedo_sfc_obs', axes(1:2), Time,                            &
         'Observed Surface Albedo', 'dimensionless',                   &
         missing_value = missing_value )
      end if


! 
!-----------------------------------------------------------------------
 
end subroutine twp_forc_diagnostic_init


!#######################################################################
!#######################################################################


subroutine update_twp_forc(time_interp,time_diag,dt_int,pdamp )
#include "fv_arrays.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a given time and returns the 
!      large-scale forcing appropriate for that given time.
!      This involves interpolating the global storage fields
!      to a given time and pressure, to do the resulting
!      calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       ------
!       INPUT:
!       ------
!
!       time_interp  time_type variable containing time we are
!                    interpolating to.
!       time_diag    diagnostic time to use for netcdf output
!       dt_int   time type variable containing the delta-t of the integration
!       pdamp        pressure for which all tendencies at p < pdamp will
!                    be set to zero.
!       As, Bs       A's and B's of half levels in hybrid coordinate
!                    phalf(k) = A(k) + B(k) * psurf
!
!       ------------
!       INPUT/OUTPUT
!       ------------
!
!       pfull   pressure in PASCALS at full model levels
!       phalf   pressure in PASCALS at half model levels
!               NOTE that it is assumed that p(j)<p(j+1)
!       pt      temperature in degrees Kelvin 
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
!       AV meridional wind  tendency due to large-scale forcing (m/sec/sec)
!       AQ      water vapor tendency due to large-scale forcing 
!                           (kg vapor/kg air/sec)
!       AL      liquid water condensate tendency due to large-scale forcing
!                           (kg condensate/ kg air/ sec)
!       AI      H2O ice condensate tendency due to large-scale forcing
!                           (kg condensate/ kg air/ sec)
!       AA      cloud fraction tendency due to large-scale forcing
!                           (fraction / sec)
!       omega_full   omega interpolated to full model levels (Pa/sec)
!
!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!       j,k            counting integers
!       itime_more   indicates indice of time_layer for which 
!                         time_layer(itime_more)>= time_interp
!       itime_less   indicates indice of time_layer for which 
!                         time_layer(itime_more)<= time_interp
!       KDIM         # of vertical levels to pt array
!       weight_more  real number indicating closeness of interpolated
!                     time to time_layer(itime_more)
!       weight_less  real number indicating closeness of interpolated
!                     time to time_layer(itime_less)
!       plev_int     array of pressure coordinates to be passed to 
!                               interp_2d_field (Pa)
!       field_int_less      1D array of field at less time
!       field_int_more      1D array of field at more time 
!       tmp_p        temporary pressure levels for interpolating
!       tmp_ans      temporary array of answer from interpolating
!       tmp_phalf    temporary half-pressure levels for interpolating
!
!       dt_seconds   time step in seconds  
!       dt_days      time step in days
!       ps_init      surface pressure initial value (Pa)
!       pi_fac       (p/pref(1000mb))**(Rd/cp) = pt/theta
!       omega_surface surface pressure tendency (Pascals/second)
!       omega_half    pressure velocity at model half levels (Pascals/second)
!       adv_vert      vertical advection of a quantity (units/second)
!       vartmp        observed value of prognostic variable (K or kg/kg)
!       nudg_t        (T_observed - T_model)/relaxation_tau (K/sec)
!       nudg_q        (q_observed - q_model)/relaxation_tau 
!                                              (kg vapor/kg air/sec)
!       dp            pressure thickness of model levels (Pa)
!       airdensity    airdensity (kg air/m3)
!       tairobs       surface skin temperature (Kelvin)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

TYPE(TIME_TYPE), intent(in)              :: time_interp,time_diag,dt_int
REAL,  INTENT (IN)                       :: pdamp


!  Internal variables
!  ------------------
real,  dimension(size(pt,1),size(pt,2))    :: tairobs
real,  dimension(size(pt,1),size(pt,2),size(pt,3)) :: pfull,t_obs,q_obs
real,  dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: phalf,omega_half
real,  dimension(size(pt,3))               :: tmp_p,tmp_ans
real,  dimension(size(pt,3)+1)             :: tmp_phalf, tmp_ans_half

INTEGER                                  :: jstart, nlev_interp
INTEGER                                  :: i,j,k,itime_less,itime_more,KDIM
INTEGER                                  :: dt_seconds,dt_days
REAL                                     :: weight_less,weight_more
REAL,     DIMENSION(ivert+2)             :: plev_int
REAL,     DIMENSION(ivert+2)             :: field_int_less,field_int_more
INTEGER                                  :: days, months, years,seconds,&
                                            minutes, hours
REAL                                     :: omega_surface
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2))  :: ps_init, tmp_res2
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2))  :: SEN_FLUX, LAT_FLUX
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2))  :: ALBEDO_DIAG
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2))  :: tmp_res3, mass_tot
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2),SIZE(pfull,3)) :: adv_vert,dT_adi, dT_vadv,dqv_vadv
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2),SIZE(pfull,3)) :: tmp_res
REAL, DIMENSION(SIZE(pfull,3))           :: vartmp
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2),SIZE(pfull,3)) :: nudg_t,  &
                                                              nudg_q
REAL, DIMENSION(SIZE(pfull,1),SIZE(pfull,2),SIZE(pfull,3)) :: airdensity
CHARACTER*1                              :: adum
LOGICAL                                  :: used

#include "fv_point.inc"

!CODE
!

!check to see if time is outside of range 
if (time_interp < time_layer(1)) then
    call get_date(time_interp,years,months,days,hours,minutes,seconds)
    if (mpp_pe() == mpp_root_pe()) then
    print *, 'DATE REQUESTED: ',years,months,days,hours,minutes,seconds
    call error_mesg ( 'update_twp_forc in scm_twp_mod',                &
         'requested time is before beginning of variational forcing '//&
         'dataset', FATAL )
    end if
end if
if (time_interp > time_layer(itime)) then
    call get_date(time_interp,years,months,days,hours,minutes,seconds)
    if (mpp_pe() == mpp_root_pe()) then
    print *, 'DATE REQUESTED: ',years,months,days,hours,minutes,seconds
    call error_mesg ( 'update_twp_forc in scm_twp_mod',                &
         'requested time is after the end of variational forcing '//   &
         'dataset', FATAL )
    end if       
end if

! --- find out # of vertical levels
KDIM = size(pt,3)
! --- compute timestep interval in seconds and days
call get_time(dt_int,dt_seconds,dt_days)
dt_seconds = dt_seconds + 86400*dt_days

! --- UPDATE PFULL AND PHALF, COMPUTE OMEGA SURFACE, AND PI FACTOR -------- !

      ! --- find out time in sfc pressure field ---- !
      call interp_time(time_interp,time_sfc,&
                        itime_less,itime_more,weight_less,weight_more)
 
      !---- interpolate surface pressure field
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=ps_forc(itime_less)
            field_int_more(:)=ps_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            ps_init(:,:) = tmp_ans(1)
            ps = ps_init
      !---- interpolate surface pressure tendency field
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=dps_dt_forc(itime_less)
            field_int_more(:)=dps_dt_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            omega_surface = tmp_ans(1)

       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pfull(i,j,:), phalf(i,j,:))
         enddo
       enddo

      ! --- compute starting level for variational dataset ----- !
      jstart =1
      do while (ps_init(1,1).lt.plev_forc(jstart))
          jstart = jstart + 1
      enddo


   
! ----------------- COMPUTE OMEGA ON FULL LEVELS --------------------------- !

! --- find out time in layered field ---- !
      call interp_time(time_interp,time_layer,&
                        itime_less,itime_more,weight_less,weight_more)

      plev_int(1)=ps_init(1,1)
      field_int_less(1)=dps_dt_forc(itime_less)
      field_int_more(1)=dps_dt_forc(itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    omega_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    omega_forc(jstart:ivert,itime_more)             
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=omega_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=omega_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3      
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_p(:)=pfull(1,1,:)


      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)

      !---- assign to omega vector
      do k=1,KDIM
           omga(:,:,k) = tmp_ans(k) 
      enddo

! ----------------- COMPUTE OMEGA ON HALF LEVELS --------------------------- !

      plev_int(1)=ps_init(1,1)
      field_int_less(1)=dps_dt_forc(itime_less)
      field_int_more(1)=dps_dt_forc(itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    omega_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    omega_forc(jstart:ivert,itime_more)             
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=omega_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=omega_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3      
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_phalf(:)=phalf(1,1,:)

      !---- do interpolation
      call interp_2d_field(tmp_phalf,weight_less,weight_more,          &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans_half)
            
      !---- assign to omega vectors
      do k=1,KDIM+1
           omega_half(:,:,k) = tmp_ans_half(k)
      enddo


! ----------------- PRODUCE TEMPERATURE FORCING ---------------------------- !
! --- YLL per TWPICE need, only use hor_adv_tendency

! --- find out time in layered field ---- !
      call interp_time(time_interp,time_layer,&
                        itime_less,itime_more,weight_less,weight_more)
 

! --- initialize fields
      t_dt   = 0.
      nudg_t = 0.
     dT_vadv = 0.0
     dT_adi  = 0.0
     dqv_vadv = 0.0
select case(FORCMETHOD)
case("D")


      !interpolate TOTAL ADVECTIVE T TENDENCY
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=0.
      field_int_more(1)=0.
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    tadv_hor_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    tadv_hor_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=tadv_hor_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=tadv_hor_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3      
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_p(:)=pfull(1,1,:)

      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,          &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
            
      !---- assign answer to AT
      do k=1,KDIM
             t_dt(:,:,k)   = tmp_ans(k)
      enddo
      if (id_Tadv_hor_forc > 0) then
        used = send_data (id_Tadv_hor_forc, t_dt(:,:,:), time_diag, 1, 1, 1)
      end if  
 
!140325 wfc and cjs
! Laura Davies TWP-ICE Single Column Model case (May 21, 2009)
! http://users.monash.edu.au/~ladavies/SCM_TWP-ICEdoc.pdf
! Do not use the vertical advection terms. refer to Figure 3.
! 140326 wfc  This was not very good. Temperatures of +60K (similar to what Figure 3 "fixes") occurred!
! yll add T_vadv
      dT_adi = rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pfull(:,:,:)
      select case (temp_vert_advec_scheme)
        case(1)
                  vadvec_scheme = SECOND_CENTERED
        case(2)
                  vadvec_scheme = FOURTH_CENTERED
        case(3)
                  vadvec_scheme = FINITE_VOLUME_LINEAR
        case(4)
                  vadvec_scheme = FINITE_VOLUME_PARABOLIC
        case(5)
                  vadvec_scheme = SECOND_CENTERED_WTS
        case(6)
                  vadvec_scheme = FOURTH_CENTERED_WTS
      end select
      call vert_advection(REAL(dt_seconds),omega_half,delp,pt, &
           dT_vadv,scheme=vadvec_scheme, form=ADVECTIVE_FORM)
      t_dt = t_dt + dT_vadv + dT_adi

      if (id_vadvT > 0) then
        used = send_data (id_vadvT, dT_vadv(:,:,:), time_diag, 1, 1, 1)
      end if
      if (id_adiT > 0) then
        used = send_data (id_adiT, dT_adi(:,:,:), time_diag, 1, 1, 1)
      end if  
!wfc and cjs
case("F")

      !interpolate TOTAL ADVECTIVE T TENDENCY
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=0.
      field_int_more(1)=0.
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    tadv_tot_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    tadv_tot_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=tadv_tot_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=tadv_tot_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3      
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_p(:)=pfull(1,1,:)

      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,          &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
            
      !---- assign answer to AT
      do k=1,KDIM
             t_dt(:,:,k)   = tmp_ans(k)
      enddo
 
      if (id_Tadv_hor_forc > 0) then
        used = send_data (id_Tadv_hor_forc, t_dt(:,:,:), time_diag, 1, 1, 1)
      end if  

  !---- COMPUTE NUDGING TERM ------------------------------------!


      !---- interpolate t_forc
      plev_int(1)=ps_init(1,1)
!The surface field is Theta but the other levels are T ????
      field_int_less(1)=t_forc(jstart,itime_less) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      field_int_more(1)=t_forc(jstart,itime_more) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=t_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=t_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
    
      !---- COMPUTE NUDGING term
      do k=1,KDIM
          nudg_t(:,:,k)=(tmp_ans(k)-pt(:,:,k))/                         &
                     max(relaxation_tau,REAL(dt_seconds))          
      enddo

  !---- COMPUTE FINAL ADVECTIVE T TENDENCY FROM THE SUM OF :  ---------!
  !---- TOTAL ADVECTIVE TENDENCY AND DAMPING TERM             ---------!
      
      do k=1,KDIM
             t_dt(:,:,k)= t_dt(:,:,k) + nudg_t(:,:,k)
      enddo
  
case("N")

  !---- COMPUTE NUDGING TERM ------------------------------------!


      !---- interpolate t_forc
      plev_int(1)=ps_init(1,1)
!The surface field is Theta but the other levels are T ????
      field_int_less(1)=t_forc(jstart,itime_less) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      field_int_more(1)=t_forc(jstart,itime_more) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=t_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=t_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
    
      !---- COMPUTE NUDGING term
      do k=1,KDIM
          nudg_t(:,:,k)=(tmp_ans(k)-pt(:,:,k))/                         &
                     max(relaxation_tau,REAL(dt_seconds))          
      enddo

  !---- COMPUTE FINAL ADVECTIVE T TENDENCY FROM THE SUM OF :  ---------!
  !---- TOTAL ADVECTIVE TENDENCY AND DAMPING TERM             ---------!
      
      do k=1,KDIM
             t_dt(:,:,k)=t_dt(:,:,k)+nudg_t(:,:,k)
      enddo
  
case("S")
      ! Semi prognostic

      !---- interpolate t_forc
      plev_int(1)=ps_init(1,1)
!The surface field is Theta but the other levels are T ????
      field_int_less(1)=t_forc(jstart,itime_less) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      field_int_more(1)=t_forc(jstart,itime_more) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=t_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=t_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)

      !---- No tendency when using fixed forcing.
      do k=1,KDIM
          t_dt(:,:,k) = 0.0
          pt(:,:,k) = tmp_ans(k)
      enddo

case default
        call error_mesg ( 'update_twp_forc in scm_twp_mod',            &
                 'Bad forcing method : '//FORCMETHOD, FATAL )
end select

 
! ----------------- PRODUCE MOISTURE FORCING ------------------------------ !

! --- find out time in layered field ---- !
      call interp_time(time_interp,time_layer,&
                        itime_less,itime_more,weight_less,weight_more)
 
! --- initialize fields
      q_dt(:,:,:,nsphum)     = 0.
      nudg_q = 0.

select case(FORCMETHOD)
case("D")


      !interpolate TOTAL ADVECTIVE Q TENDENCY
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=0.
      field_int_more(1)=0.
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    qadv_hor_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    qadv_hor_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=qadv_hor_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=qadv_hor_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3      
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_p(:)=pfull(1,1,:)

      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
      
      !---- assign answer to AQ
      do k=1,KDIM
             q_dt(:,:,k,nsphum)=tmp_ans(k)
      enddo
      if (id_Qadv_hor_forc > 0) then
        used = send_data (id_Qadv_hor_forc, q_dt(:,:,:,nsphum), time_diag, 1, 1, 1)
      end if  
    
      select case (tracer_vert_advec_scheme)
        case(1)
                  vadvec_scheme = SECOND_CENTERED
        case(2)
                  vadvec_scheme = FOURTH_CENTERED
        case(3)
                  vadvec_scheme = FINITE_VOLUME_LINEAR
        case(4)
                  vadvec_scheme = FINITE_VOLUME_PARABOLIC
        case(5)
                  vadvec_scheme = SECOND_CENTERED_WTS
        case(6)
                  vadvec_scheme = FOURTH_CENTERED_WTS
      end select
      call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nsphum), &
           dqv_vadv,scheme=vadvec_scheme, form=ADVECTIVE_FORM)
      q_dt(:,:,:,nsphum) = q_dt(:,:,:,nsphum) + dqv_vadv
      if (id_vadvQ > 0) then
        used = send_data (id_vadvQ, dqv_vadv(:,:,:), time_diag, 1, 1, 1)
      end if

case("F")

      !interpolate TOTAL ADVECTIVE Q TENDENCY
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=0.
      field_int_more(1)=0.
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    qadv_tot_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    qadv_tot_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=qadv_tot_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=qadv_tot_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3      
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_p(:)=pfull(1,1,:)

      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
      
      !---- assign answer to AQ
      do k=1,KDIM
          q_dt(:,:,k,nsphum)=tmp_ans(k)
      enddo
      if (id_Qadv_hor_forc > 0) then
        used = send_data (id_Qadv_hor_forc, q_dt(:,:,:,nsphum), time_diag, 1, 1, 1)
      end if  

  !---- COMPUTE NUDGING TERM ------------------------------------!

      !---- interpolate q_forc
          
      tmp_ans = get_column(q_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                           weight_less,weight_more, pfull(1,1,:), q_stand)
      !---- COMPUTE NUDGING term
      do k=1,KDIM
          nudg_q(:,:,k)=(tmp_ans(k)-q(:,:,k,nsphum))/                  &
                     max(relaxation_tau,REAL(dt_seconds))          
      enddo
      

  !---- COMPUTE FINAL ADVECTIVE Q TENDENCY FROM THE SUM OF    ---------!
  !---- TOTAL ADVECTIVE TENDENCY AND NUDGING TERM             ---------!
      
      do k=1,KDIM
             q_dt(:,:,k,nsphum)=q_dt(:,:,k,nsphum)+nudg_q(:,:,k)
      enddo
    
case("N")

  !---- COMPUTE NUDGING TERM ------------------------------------!

      !---- interpolate q_forc

      tmp_ans = get_column(q_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                           weight_less,weight_more, pfull(1,1,:),q_stand)
      !---- COMPUTE NUDGING term
      do k=1,KDIM
          nudg_q(:,:,k)=(tmp_ans(k)-q(:,:,k,nsphum))/                  &
                     max(relaxation_tau,REAL(dt_seconds))          
      enddo
      

  !---- COMPUTE FINAL ADVECTIVE Q TENDENCY FROM THE SUM OF    ---------!
  !---- TOTAL ADVECTIVE TENDENCY AND NUDGING TERM             ---------!
      
      do k=1,KDIM
             q_dt(:,:,k,nsphum)=nudg_q(:,:,k)
      enddo
    
case("S")

      ! Semi prognostic 

      !---- interpolate q_forc

      tmp_ans = get_column(q_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                           weight_less,weight_more, pfull(1,1,:),q_stand)

      !---- No tendency when using fixed forcing.
      do k=1,KDIM
          q(:,:,k,nsphum) = tmp_ans(k)
          q_dt(:,:,k,nsphum)=0.0
      enddo

case default
        call error_mesg ( 'update_twp_forc in scm_twp_mod',            &
                 'Bad forcing method : '//FORCMETHOD, FATAL )
end select

 
! ----------------- PRODUCE MOMENTUM FORCING ------------------------------ !

! --- initialize fields
      u_dt = 0.
      v_dt = 0.
      
select case(FORCMETHOD)
case("D","F","N")

!---- interpolate u_forc

      tmp_ans = get_column(u_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                           weight_less,weight_more, pfull(1,1,:))
      !---- COMPUTE NUDGING term
      do k=1,KDIM
          u_dt(:,:,k)=(tmp_ans(k)-ua(:,:,k))/                            &
                     max(wind_relaxation_tau,REAL(dt_seconds))
      enddo
   
      !---- output netcdf result
      tmp_res(1,1,:) = tmp_ans
      if (id_uwnd_obs > 0) then
      used = send_data (id_uwnd_obs, tmp_res(:,:,:), time_diag, 1, 1, 1)
      end if  
   
!---- interpolate v_forc

      tmp_ans = get_column(v_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                           weight_less,weight_more, pfull(1,1,:))
      !---- COMPUTE NUDGING term
      do k=1,KDIM          
          v_dt(:,:,k)=(tmp_ans(k)-va(:,:,k))/                            &
                     max(wind_relaxation_tau,REAL(dt_seconds))
      enddo

      !---- output netcdf result
      tmp_res(1,1,:) = tmp_ans
      if (id_vwnd_obs > 0) then
      used = send_data (id_vwnd_obs, tmp_res(:,:,:), time_diag, 1, 1, 1)
      end if  
    
case("S")

  !---- FORCE SEMI PROGNOSTIC ------------------------------------!
      !Force U and V to calm conditions
      ua = 0.0
      va = 0.0
      
 case default
        call error_mesg ( 'update_twp_forc in scm_gate_mod',            &
                 'Bad forcing method : '//FORCMETHOD, FATAL )
end select


! ----------------- PRODUCE CONDENSATE FORCING ------------------------------ !

      q_dt(:,:,:,nql)=0.
      q_dt(:,:,:,nqi)=0. 
      q_dt(:,:,:,nqa)=0.
      if(nqn > 0 ) q_dt(:,:,:,nqn)=0.  ! h1g
      if(nqni > 0 ) q_dt(:,:,:,nqni)=0.  ! h1g

      if (vert_advec_cond) then
       
          select case (tracer_vert_advec_scheme)
            case(1)
                      vadvec_scheme = SECOND_CENTERED
            case(2)
                      vadvec_scheme = FOURTH_CENTERED
            case(3)
                      vadvec_scheme = FINITE_VOLUME_LINEAR
            case(4)
                      vadvec_scheme = FINITE_VOLUME_PARABOLIC
            case(5)
                      vadvec_scheme = SECOND_CENTERED_WTS
            case(6)
                      vadvec_scheme = FOURTH_CENTERED_WTS
          end select

          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nql),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
          q_dt(:,:,:,nql) = q_dt(:,:,:,nql) + adv_vert    

          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqi),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
          q_dt(:,:,:,nqi) = q_dt(:,:,:,nqi) + adv_vert

          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqa),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
          q_dt(:,:,:,nqa) = q_dt(:,:,:,nqa) + adv_vert

          if (nqn > 0 ) then 
            call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqn),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
            q_dt(:,:,:,nqn) = q_dt(:,:,:,nqn) + adv_vert
          endif

          if ( nqni > 0 ) then
            call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqni),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
            q_dt(:,:,:,nqni) = q_dt(:,:,:,nqni) + adv_vert
         endif
! <--- h1g

      end if

      if (do_qc_adv_loss) then

          !---- interpolate tau_inv_forc

          tmp_ans = get_column(tau_inv_forc,ps_init(1,1),itime_less,itime_more,ivert,jstart,&
                               weight_less,weight_more, pfull(1,1,:))
          !---- compute condensate loss term
          do k=1,KDIM
              q_dt(:,:,k,nql)= q_dt(:,:,k,nql)  - (q(:,:,k,nql)*tmp_ans(k))
              q_dt(:,:,k,nqi)= q_dt(:,:,k,nqi)  - (q(:,:,k,nqi)*tmp_ans(k))
              q_dt(:,:,k,nqa)= q_dt(:,:,k,nqa)  - (q(:,:,k,nqa)*tmp_ans(k))                    
          enddo

      end if

! --- COMPUTE SURFACE DATA: SENSIBLE AND LATENT HEAT FLUX AND -------- !
! ---                       SURFACE ALBEDO AND TSKIN --------- !

      ! --- do fluxes first
      ! --- find out time in sfc pressure field ---- !
      call interp_time(time_interp,time_sfc,&
                        itime_less,itime_more,weight_less,weight_more)


      !---- interpolate sensible heat flux field
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=sen_flux_forc(itime_less)
            field_int_more(:)=sen_flux_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            SEN_FLUX(:,:) = tmp_ans(1)
            SENFLUX = tmp_ans(1)

      !---- interpolate latent heat flux field
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=lat_flux_forc(itime_less)
            field_int_more(:)=lat_flux_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            LAT_FLUX(:,:) = tmp_ans(1)
            EVAPFLUX = tmp_ans(1)/hlv

      !---- interpolate TAIR_SFC field
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=tair_surf_forc(itime_less)
            field_int_more(:)=tair_surf_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tairobs(:,:) = tmp_ans(1)
            
            !---- interpolate tskin_forc (if available) to retrieve tskin
            
            if (allocated(tskin_forc)) THEN
     
            ! --- interpolate tskin_forc
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=tskin_forc(itime_less)
            field_int_more(:)=tskin_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,        &
                 plev_int,field_int_less,field_int_more,tmp_ans)
            TSKIN = tmp_ans(1)
            
            else
            
            !---- interpolate lwup_sfc (if available) to retrieve tskin
            
            if (allocated(lwup_sfc)) THEN
     
            ! --- interpolate lwup_sfc
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=lwup_sfc(itime_less)
            field_int_more(:)=lwup_sfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,        &
                 plev_int,field_int_less,field_int_more,tmp_ans)
            TSKIN = sqrt(sqrt(tmp_ans(1)/stefan))
            
            else
      
            ! --- if lwup_sfc not available use surface air temp   --- !      
     
            ! --- interpolate TAIR_SFC field
            tmp_p(:)=50000.
            plev_int(:)=p00
            field_int_less(:)=tair_surf_forc(itime_less)
            field_int_more(:)=tair_surf_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,        &
                 plev_int,field_int_less,field_int_more,tmp_ans)
            TSKIN = tmp_ans(1)
       
            end if
            
            end if
            
            !---- interpolate surface albedo field
            if (do_iop_mean_albedo) then
                 ALBEDO_OBS = avg_albedo_sfc(1)
                 ALBEDO_DIAG(:,:) = avg_albedo_sfc(1)
            else
                 tmp_p(:)=50000.
                 plev_int(:)=p00
                 field_int_less(:)=albedo_sfc(itime_less)
                 field_int_more(:)=albedo_sfc(itime_more)
                 call interp_2d_field(tmp_p,weight_less,weight_more,   &
                      plev_int,field_int_less,field_int_more,tmp_ans)
                 ALBEDO_DIAG(:,:) = tmp_ans(1)
                 ALBEDO_OBS = tmp_ans(1) 
            end if
      
! ----------------- do NETCDF DIAGNOSTIC OUTPUT ------------------------- !
!



      !------------  do DIAGNOSTICS  ---------------------!
      !------------ ALSO OUTPUT OBSERVED T and Q  ---------------------!

      ! --- find out time in layered field ---- !
      call interp_time(time_interp,time_layer,&
                        itime_less,itime_more,weight_less,weight_more)
 
      !---- interpolate t_forc
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=t_forc(jstart,itime_less) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      field_int_more(1)=t_forc(jstart,itime_more) * &
                  ( (  (ps_init(1,1)/plev_forc(jstart)) )**(rdgas/cp_air) )
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    t_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=t_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=t_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !
      !adjust nearest surface pressure so that plev reflects the
      !mean pressure of the layer
      !
      ! Note that this adjustment is performed only for netcdf diagnostics
      ! so that evaluation of sounding can take place against truer
      ! observations.  Initialization of the SCM as well as forcing
      ! fields do not have this correction so as to be consistent with
      ! the way all the other SCM SCM models are forced
      !
!      plev_int(2)=0.5*(0.5*(plev_forc(jstart)+plev_forc(jstart+1))+plev_int(1))
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
            
      ! Set to missing above top level of variational analysis
      !
      where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
      
      !---- output netcdf result
      t_obs(1,1,:) = tmp_ans
      if (id_temp_obs > 0) then
      used = send_data ( id_temp_obs, &
          t_obs(:,:,:), time_diag, 1, 1, 1)
      end if  
   
      !---- interpolate q_forc
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=q_forc(jstart,itime_less)
      field_int_more(1)=q_forc(jstart,itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    q_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    q_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=q_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=q_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !
      !adjust nearest surface pressure so that plev reflects the
      !mean pressure of the layer
      !
      ! Note that this adjustment is performed only for netcdf diagnostics
      ! so that evaluation of sounding can take place against truer
      ! observations.  Initialization of the SCM as well as forcing
      ! fields do not have this correction so as to be consistent with
      ! the way all the other SCM SCM models are forced
      !
!      plev_int(2)=0.5*(0.5*(plev_forc(jstart)+plev_forc(jstart+1))+plev_int(1))
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
         
      ! Set to missing above top level of variational analysis
      ! also convert kg/kg to g/kg
      
      where (tmp_p .lt. plev_forc(ivert))
           tmp_ans= missing_value
      elsewhere
           tmp_ans= 1000.*tmp_ans
      end where
      
      !---- output netcdf result
      q_obs(1,1,:) = tmp_ans
      if (id_qv_obs > 0) then
      used = send_data( id_qv_obs, &
          q_obs(:,:,:), time_diag, 1, 1, 1)
      end if  
   
      !---- interpolate omega_forc
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=dps_dt_forc(itime_less)
      field_int_more(1)=dps_dt_forc(itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    omega_forc(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    omega_forc(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=omega_forc(ivert,itime_less)
      field_int_more(ivert-jstart+3)=omega_forc(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      where (plev_int .lt. pdamp) 
           field_int_less = 0.
           field_int_more = 0.
      endwhere
      tmp_p(:)=pfull(1,1,:)
      
      !
      !adjust nearest surface pressure so that plev reflects the
      !mean pressure of the layer
      !
      ! Note that this adjustment is performed only for netcdf diagnostics
      ! so that evaluation of sounding can take place against truer
      ! observations.  Initialization of the SCM as well as forcing
      ! fields do not have this correction so as to be consistent with
      ! the way all the other SCM SCM models are forced
      !
      plev_int(2)=0.5*(0.5*(plev_forc(jstart)+plev_forc(jstart+1))+plev_int(1))
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
              
      !---- output netcdf result (convert Pa/s to hPa/day)
      tmp_res(1,1,:) = 86400.*tmp_ans(:)/100.
      if (id_omega > 0) then
      used = send_data( id_omega, &
          tmp_res(:,:,:), time_diag, 1, 1, 1)
      end if  
      
      !---- interpolate q1_obs
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=q1_obs(jstart,itime_less)
      field_int_more(1)=q1_obs(jstart,itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    q1_obs(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    q1_obs(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=q1_obs(ivert,itime_less)
      field_int_more(ivert-jstart+3)=q1_obs(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !
      !adjust nearest surface pressure so that plev reflects the
      !mean pressure of the layer
      !
      ! Note that this adjustment is performed only for netcdf diagnostics
      ! so that evaluation of sounding can take place against truer
      ! observations.  Initialization of the SCM as well as forcing
      ! fields do not have this correction so as to be consistent with
      ! the way all the other SCM SCM models are forced
      !
      plev_int(2)=0.5*(0.5*(plev_forc(jstart)+plev_forc(jstart+1))+plev_int(1))
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
              
      !---- output netcdf result
      tmp_res(1,1,:) = tmp_ans(:)
      if (id_q1_obs > 0) then
      used = send_data( id_q1_obs, &
          tmp_res(:,:,:), time_diag, 1, 1, 1)
      end if  
      
      !---- interpolate q2_obs
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=q2_obs(jstart,itime_less)
      field_int_more(1)=q2_obs(jstart,itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    q2_obs(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    q2_obs(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=q2_obs(ivert,itime_less)
      field_int_more(ivert-jstart+3)=q2_obs(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !
      !adjust nearest surface pressure so that plev reflects the
      !mean pressure of the layer
      !
      ! Note that this adjustment is performed only for netcdf diagnostics
      ! so that evaluation of sounding can take place against truer
      ! observations.  Initialization of the SCM as well as forcing
      ! fields do not have this correction so as to be consistent with
      ! the way all the other SCM SCM models are forced
      !
      plev_int(2)=0.5*(0.5*(plev_forc(jstart)+plev_forc(jstart+1))+plev_int(1))
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
              
      !---- output netcdf result
      tmp_res(1,1,:) = tmp_ans(:)
      if (id_q2_obs > 0) then
      used = send_data( id_q2_obs, &
          tmp_res(:,:,:), time_diag, 1, 1, 1)
      end if  

      if (nlayer_fields.ge.19) then     
      
      !---- interpolate arscl
      plev_int(1)=ps_init(1,1)
      field_int_less(1)=arscl(jstart,itime_less)
      field_int_more(1)=arscl(jstart,itime_more)
      plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
      field_int_less(2:(ivert-jstart+2))=&
                    arscl(jstart:ivert,itime_less)
      field_int_more(2:(ivert-jstart+2))=&
                    arscl(jstart:ivert,itime_more)
      plev_int(ivert-jstart+3)=0.
      field_int_less(ivert-jstart+3)=arscl(ivert,itime_less)
      field_int_more(ivert-jstart+3)=arscl(ivert,itime_more)
      nlev_interp = ivert-jstart+3 
      tmp_p(:)=pfull(1,1,:)
      
      !
      !adjust nearest surface pressure so that plev reflects the
      !mean pressure of the layer
      !
      ! Note that this adjustment is performed only for netcdf diagnostics
      ! so that evaluation of sounding can take place against truer
      ! observations.  Initialization of the SCM as well as forcing
      ! fields do not have this correction so as to be consistent with
      ! the way all the other SCM SCM models are forced
      !
      plev_int(2)=0.5*(0.5*(plev_forc(jstart)+plev_forc(jstart+1))+plev_int(1))
      
      !---- do interpolation
      call interp_2d_field(tmp_p,weight_less,weight_more,              &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
              
      !---- output netcdf result
      tmp_res(1,1,:) = tmp_ans(:)
      if (id_arscl > 0) then
      used = send_data( id_arscl, &
          tmp_res(:,:,:), time_diag, 1, 1, 1)
      end if  

      end if
            
      !---- nudging tendencies
      if (id_nudg_t > 0) then
           used = send_data ( id_nudg_t, nudg_t, time_diag, 1, 1, 1)
      end if
      
      if (id_nudg_q > 0) then
           used = send_data ( id_nudg_q, nudg_q, time_diag, 1, 1, 1)
      end if

! --- YLL output U, and V nudging term
      
      if (id_nudg_u > 0) &
           used = send_data ( id_nudg_u, u_dt, time_diag, 1, 1, 1)
      if (id_nudg_v > 0) &
           used = send_data ( id_nudg_v, v_dt, time_diag, 1, 1, 1)
      if (id_tdt_lf > 0)   &
           used = send_data( id_tdt_lf,   t_dt (:,:,:), time_diag, 1, 1)
      if (id_qvdt_lf > 0)  &
           used = send_data( id_qvdt_lf,  q_dt (:,:,:,nsphum), time_diag, 1, 1)



      !---- output sensible and latent heat fluxes, tairobs, and albedo
      if (id_latent_obs > 0) &
           used = send_data( id_latent_obs, LAT_FLUX, time_diag, 1, 1)
      if (id_sens_obs > 0) &
           used = send_data( id_sens_obs,   SEN_FLUX, time_diag, 1, 1)
      if (id_tair_obs > 0) &
           used = send_data( id_tair_obs,   tairobs,  time_diag, 1, 1)
      if (id_albedo_sfc_obs > 0) &
           used = send_data( id_albedo_sfc_obs, ALBEDO_DIAG, time_diag, 1, 1)
           
      ! --- find out time in sfc pressure field ---- !
      call interp_time(time_interp,time_sfc,&
                        itime_less,itime_more,weight_less,weight_more)
      
      !---- set tmp_p and plev_int for surface fields
      tmp_p(:)=50000.
      plev_int(:)=p00

! --- YLL output tskin_forc instead of tground            
      !---- interpolate tground
            field_int_less(:)=tskin_forc(itime_less)
            field_int_more(:)=tskin_forc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_tground_obs > 0) then
            used = send_data( id_tground_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate rh
            field_int_less(:)=rh(itime_less)
            field_int_more(:)=rh(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_rh_obs > 0) then
            used = send_data( id_rh_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
      !---- interpolate uwnd
            field_int_less(:)=uwnd(itime_less)
            field_int_more(:)=uwnd(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_uref_obs > 0) then
            used = send_data( id_uref_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
      !---- interpolate vwnd
            field_int_less(:)=vwnd(itime_less)
            field_int_more(:)=vwnd(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_vref_obs > 0) then
            used = send_data( id_vref_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      

     
      !---- interpolate precip
            field_int_less(:)=precip(itime_less)
            field_int_more(:)=precip(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_precip_obs > 0) then
            used = send_data( id_precip_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if
            if ( nprec_clo > 0 ) then
              q(:,:,:,nprec_clo) = tmp_res2(1,1) * a1_precfrac
            endif
 
       !---- interpolate closure field
            if ( nq2_clo > 0 ) then
              tmp_p(:)=50000.
              plev_int(:)=p00
              field_int_less(:)= -1.0 * gFq(itime_less) - q2_integral(itime_less)
              field_int_more(:)= -1.0 * gFq(itime_more) - q2_integral(itime_more)
              call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
                   field_int_less,field_int_more,tmp_ans)
              q(:,:,:,nq2_clo) = tmp_ans(1) 
            endif

      !---- interpolate netraddnsfc
            field_int_less(:)=netraddnsfc(itime_less)
            field_int_more(:)=netraddnsfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_netdnsfc_obs > 0) then
            used = send_data( id_netdnsfc_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
      !---- interpolate olr
            field_int_less(:)=olr(itime_less)
            field_int_more(:)=olr(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_olr_obs > 0) then
            used = send_data( id_olr_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate swup_toa
            field_int_less(:)=swup_toa(itime_less)
            field_int_more(:)=swup_toa(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_swup_toa_obs > 0) then
            used = send_data( id_swup_toa_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate swdn_toa
            field_int_less(:)=swdn_toa(itime_less)
            field_int_more(:)=swdn_toa(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_swdn_toa_obs > 0) then
            used = send_data( id_swdn_toa_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate tot_cld_amt
            field_int_less(:)=tot_cld_amt(itime_less)
            field_int_more(:)=tot_cld_amt(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_tot_cld_amt_obs > 0) then
            used = send_data( id_tot_cld_amt_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate lwp
            field_int_less(:)=lwp(itime_less)
            field_int_more(:)=lwp(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_lwp_obs > 0) then
            used = send_data( id_lwp_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate wvp
            if (allocated(wvp)) then
            field_int_less(:)=wvp(itime_less)
            field_int_more(:)=wvp(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_pw_obs > 0) then
            used = send_data( id_pw_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if        
            end if

       
       !---- interpolate low_cld_amt
            field_int_less(:)=low_cld_amt(itime_less)
            field_int_more(:)=low_cld_amt(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_low_cld_amt_obs > 0) then
            used = send_data( id_low_cld_amt_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate mid_cld_amt
            field_int_less(:)=mid_cld_amt(itime_less)
            field_int_more(:)=mid_cld_amt(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_mid_cld_amt_obs > 0) then
            used = send_data( id_mid_cld_amt_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate hgh_cld_amt
            field_int_less(:)=hgh_cld_amt(itime_less)
            field_int_more(:)=hgh_cld_amt(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_hgh_cld_amt_obs > 0) then
            used = send_data( id_hgh_cld_amt_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate cld_thk
            field_int_less(:)=cld_thk(itime_less)
            field_int_more(:)=cld_thk(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_cld_thk_obs > 0) then
            used = send_data( id_cld_thk_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate cld_hgt
            field_int_less(:)=cld_hgt(itime_less)
            field_int_more(:)=cld_hgt(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_cld_hgt_obs > 0) then
            used = send_data( id_cld_hgt_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
      
       !---- interpolate qtend
            field_int_less(:)=qtend(itime_less)
            field_int_more(:)=qtend(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_qtend_obs > 0) then
            used = send_data( id_qtend_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate advqtend
            field_int_less(:)=advqtend(itime_less)
            field_int_more(:)=advqtend(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_advqtend_obs > 0) then
            used = send_data( id_advqtend_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate evap
            field_int_less(:)=evap(itime_less)
            field_int_more(:)=evap(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_evap_obs > 0) then
            used = send_data( id_evap_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate ttend
            field_int_less(:)=ttend(itime_less)
            field_int_more(:)=ttend(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_ttend_obs > 0) then
            used = send_data( id_ttend_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate advttend
            field_int_less(:)=advttend(itime_less)
            field_int_more(:)=advttend(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_advttend_obs > 0) then
            used = send_data( id_advttend_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate radttend
            field_int_less(:)=radttend(itime_less)
            field_int_more(:)=radttend(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_radttend_obs > 0) then
            used = send_data( id_radttend_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate latttend
            field_int_less(:)=latttend(itime_less)
            field_int_more(:)=latttend(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_latttend_obs > 0) then
            used = send_data( id_latttend_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
      
       if (allocated(lwdn_sfc)) then
       !---- interpolate lwdn_sfc
            field_int_less(:)=lwdn_sfc(itime_less)
            field_int_more(:)=lwdn_sfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_lwdn_sfc_obs > 0) then
            used = send_data( id_lwdn_sfc_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate lwup_sfc
            field_int_less(:)=lwup_sfc(itime_less)
            field_int_more(:)=lwup_sfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_lwup_sfc_obs > 0) then
            used = send_data( id_lwup_sfc_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  

       !---- interpolate swdn_sfc
            field_int_less(:)=swdn_sfc(itime_less)
            field_int_more(:)=swdn_sfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_swdn_sfc_obs > 0) then
            used = send_data( id_swdn_sfc_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
       !---- interpolate swup_sfc
            field_int_less(:)=swup_sfc(itime_less)
            field_int_more(:)=swup_sfc(itime_more)
            call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
            tmp_res2(1,1)=tmp_ans(1)
            if (id_swup_sfc_obs > 0) then
            used = send_data( id_swup_sfc_obs, &
                  tmp_res2, time_diag, 1, 1)
            end if  
      
      end if
                 
! 
!-----------------------------------------------------------------------
 
end subroutine update_twp_forc


!########################################################################
!########################################################################


subroutine interp_time (time_interp,TIME_VEC,&
                        itime_less,itime_more,weight_less,weight_more)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine finds the integer corresponding to the
!      nearest time in the TIME_VEC set of times that is less than
!      and more than the given time (time_interp)
!
!      THE subroutine ASSUMES THAT TIME_VEC IS ARRANGED IN ASCendING ORDER
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       ------
!       INPUT:
!       ------
!
!       time_interp  time_type variable containing time we are
!                    interpolating to.
!       TIME_VEC     list of times we are searching through to
!                    find nearest time
!
!       -------------
!       INPUT/OUTPUT:
!       -------------
!     
!       itime_more   indicates indice of TIME_VEC for which 
!                         TIME_VEC(itime_more)>= time_interp
!       itime_less   indicates indice of TIME_VEC for which 
!                         TIME_VEC(itime_more)<= time_interp
!                      
!               NOTE if time_interp is outside of the range of
!                       TIME_VEC   itime_more and itime_less 
!                       equal 1 (for time_interp < TIME_VEC(1))
!                    or equal itime_vec
!                       (for time_interp > TIME_VEC(itime_vec))
!       
!       weight_more   real number indicating closeness of interpolated
!                     time to TIME_VEC(itime_more)
!                     (=1 if time_interp = TIME_VEC(itime_more))
!       weight_less   real number indicating closeness of interpolated
!                     time to TIME_VEC(itime_less)
!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!       itime_vec    # of times in TIME_VEC set
!       j            counting integers
!       seconds_frac, days_frac   # of seconds and days in fraction of 
!                                   time interval
!       seconds_int, days_int     # of seconds and days in time interval
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------


TYPE(TIME_TYPE),INTENT(IN)        :: time_interp
TYPE(TIME_TYPE),INTENT(IN),DIMENSION(:)        :: TIME_VEC
INTEGER, INTENT (INOUT)           :: itime_more,itime_less
REAL, INTENT (INOUT)              :: weight_more,weight_less

!  Internal variables
!  ------------------

INTEGER                           :: itime_vec,j
INTEGER                           :: seconds_frac,seconds_int
INTEGER                           :: days_frac, days_int
!
! Code
! ----

!findout size of TIME_VEC
itime_vec = SIZE(TIME_VEC,1)

! --- find nearest time to time_interp --- !
      if (time_interp <= TIME_VEC(1) .or. &
          time_interp >= TIME_VEC(itime_vec)) THEN
          if (time_interp <= TIME_VEC(1)) then
               itime_more = 1
               itime_less = 1
               weight_more = 1.
               weight_less = 0.
          end if
          if (time_interp >= TIME_VEC(itime_vec)) THEN
               itime_more = itime_vec
               itime_less = itime_vec
               weight_more = 1.
               weight_less = 0.
          end if
      ELSE 
          itime_more = 0
          do j = 2, itime_vec
              if (TIME_VEC(j) >= time_interp) THEN
              if (itime_more .eq. 0) THEN
                  itime_more = j
                  itime_less = j-1
              end if
              end if
          enddo
          call get_time((time_interp - TIME_VEC(itime_less)),&
                        seconds_frac,days_frac)
          call get_time((TIME_VEC(itime_more) - TIME_VEC(itime_less)),&
                        seconds_int,days_int)
          weight_more = real(seconds_frac+86400*days_frac) / &
                        real(seconds_int+86400*days_int)
          weight_less = 1.-weight_more
      end if

! 
!-----------------------------------------------------------------------
 
end subroutine interp_time 


!########################################################################
!########################################################################


subroutine interp_2d_field(p_interp,weight_less,weight_more,plev,&
                     field_less,field_more,field_interp)


implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the 2D time/pressure interpolated
!      value of a field given at the two time (less and more).
!      The field is linearly interpolated in time and in pressure.
!
!      THE subroutine ASSUMES THAT plev(1) > plev(2) > plev(3)...
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       ------
!       INPUT:
!       ------
!
!       p_interp     pressure levels being interpolated to (Pascals)
!       weight_more   real number indicating closeness of interpolated
!                     time to itime_more (range: 0 to 1)
!       weight_less   real number indicating closeness of interpolated
!                     time to itime_less (range: 0 to 1)
!       plev          pressure at levels of field vector (Pascals)
!       field_less    field at itime_less on the pressure levels given (plev)
!       field_more    field at itime_more on the pressure levels given (plev)
!
!       -------------
!       INPUT/OUTPUT:
!       -------------
!
!       field_interp  the field interpolated in pressure and time to the 
!                     desired levels
!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!       field_tmp    field interpolated in time using weights
!       j            counting integers
!       klev_more    integer array containing indexes of plev for which
!                        plev(klev_more) >= p_interp
!       klev_less    integer array containing indexes of plev for which
!                        plev(klev_less) <= p_interp
!       ptmp         temporary pressure
!       FIND_PINTERP logical used to identify klev_more and klev_less
!       pweight_more weighting factor for klev_more level
!       pweight_less weighting factor for klev_less level
!       nfield_lev   # of pressure levels to field data
!       ninterp      # of pressure levels to interpolate data to
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL, INTENT (IN)                 :: weight_more,weight_less
REAL, INTENT (IN), DIMENSION(:)   :: p_interp
REAL, INTENT (INOUT),DIMENSION(:) :: plev,field_less,field_more
REAL, INTENT (INOUT),DIMENSION(:) :: field_interp

!  Internal variables
!  ------------------

REAL,    DIMENSION(SIZE(field_less,1))  :: field_tmp
REAL                                    :: ptmp
INTEGER, DIMENSION(SIZE(p_interp,1))    :: klev_more,klev_less
INTEGER                                 :: j,nfield_lev,ninterp
LOGICAL, DIMENSION(SIZE(p_interp,1))    :: FIND_PINTERP
REAL,    DIMENSION(SIZE(p_interp,1))    :: pweight_more,pweight_less

!
! Code
! ----

! ---  count # of levels to the field and # of interpolating levels
nfield_lev = SIZE(field_less,1) 
ninterp = SIZE(p_interp,1)


! --- resort bottom levels if necessary
      if (plev(1).lt.plev(2)) THEN
             !change pressures first
             ptmp=plev(2)
             plev(2)=plev(1)
             plev(1)=ptmp
             !change field values
             field_more(2)=field_more(1)
             field_less(2)=field_less(1)             
      end if

! --- interpolate field in time
field_tmp = weight_more * field_more + weight_less*field_less    


! --- find pressure levels nearest to p_interp --- !
      !set find flag to false and initialize weights
      FIND_PINTERP(:)=.FALSE.
      pweight_more(:)=0.
      pweight_less(:)=0.

      WHERE(p_interp(:) .ge. plev(1))
            FIND_PINTERP(:)=.TRUE.
            klev_more(:)=1
            klev_less(:)=1
            pweight_more(:)=1.
            pweight_less(:)=0.
      end WHERE

      WHERE(p_interp(:) .le. plev(nfield_lev))
            FIND_PINTERP(:)=.TRUE.
            klev_more(:)=nfield_lev
            klev_less(:)=nfield_lev
            pweight_more(:)=1.
            pweight_less(:)=0.
      end WHERE

      do j = 2, nfield_lev
         WHERE(plev(j) .le. p_interp(:) .and. .not. FIND_PINTERP(:))
            FIND_PINTERP(:)=.TRUE.
            klev_more(:)=j
            klev_less(:)=j-1
            pweight_more(:)=(p_interp(:) - plev(klev_less))/&
                            (plev(klev_more) - plev(klev_less))
            pweight_less(:)=1.-pweight_more(:)
         end WHERE
      enddo  !for j loop  

! --- do interpolation
    field_interp(:) = pweight_more(:)*field_tmp(klev_more(:)) + &
                      pweight_less(:)*field_tmp(klev_less(:))
    

                      
! 
!-----------------------------------------------------------------------
 
end subroutine interp_2d_field 

!########################################################################
!########################################################################


subroutine rem_miss_var(bad_val,TIME_VEC,data_mat,adjlev)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine removes all values of data_mat that are
!      less than or equal to bad val and replaces them with 
!      linear interpolation in time or the vertical depending
!      on the direction in which the number of points to interpol-
!      ate over is less.  Missing values on the
!      end of the data record or at the bottom or top of the 
!      profile are taken as equal to the value of the nearest 
!      record.
!
!      data_mat is assumed to be 2dimension with the
!      first dimension being the physical coordinate and the
!      second coordinate referring to the time variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       ------
!       INPUT:
!       ------
!
!       bad_val      value for which if data_mat <= bad_val
!                    the value is replaced
!       TIME_VEC     time_type variable containing times  of data_mat
!       adjlev       # of levels to be adjusted
!
!       ------------
!       INPUT/OUTPUT:
!       ------------
!
!       data_mat     data matrix (#x values, #times)
!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!       time_less    time of data at less point
!       time_more    time of data at more point
!       itime_less   index of data_mat for less time
!       itime_more   index of data_mat for more time
!       nlev,ntime   #of levels (or times)
!       itime,jlev,icur   counting integers
!       seconds_frac, days_frac   # of seconds and days in fraction of 
!                                   time interval
!       seconds_int, days_int     # of seconds and days in time interval
!       weight_more   real number indicating closeness of interpolated
!                     time to TIME_VEC(itime_more)
!                     (=1 if time_interp = TIME_VEC(itime_more))
!       weight_less   real number indicating closeness of interpolated
!                     time to TIME_VEC(itime_less)
!
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

REAL,     INTENT (IN)                     :: bad_val
TYPE(TIME_TYPE), INTENT (IN), DIMENSION(:):: TIME_VEC
REAL,     INTENT (INOUT), DIMENSION(:,:)  :: data_mat
INTEGER,  OPTIONAL, INTENT (IN)           :: adjlev

!  Internal variables
!  ------------------

TYPE(TIME_TYPE)                          :: time_less,time_more
INTEGER                                  :: itime_less,itime_more
INTEGER                                  :: jlev_less,jlev_more,jcur
INTEGER                                  :: dpointslev,dpointstime
INTEGER                                  :: nlev,ntime,jlev,itime,icur,itime2
REAL                                     :: weight_less,weight_more
INTEGER                                  :: seconds_frac,seconds_int
INTEGER                                  :: days_frac, days_int
REAL, DIMENSION(SIZE(data_mat,1),SIZE(data_mat,2)) :: data_mat_new
!CODE
!


! --- find out # of vertical levels
nlev = SIZE(data_mat,1)
if (present(adjlev)) nlev=adjlev
ntime = SIZE(data_mat,2)

! If nlev .le. 1 then jlev_more and jlev_less are not assigned values but are
! used in
!      dpointslev = abs(jlev_more-jlev) + abs(jlev_less-jlev)
! Therefore assign values here.

if ( nlev .eq. 1 ) then 
  jlev_more = nlev
  jlev_less = nlev
endif

! --- assign data_mat_new with old values
data_mat_new = data_mat

! --------  do interpolation in time -------------------------------- !

! --- start general loop
do jlev=1,nlev
do itime=1,ntime

       !find out if data is bad
       if (data_mat(jlev,itime).le.bad_val) THEN

              !find less time
              icur=itime
              do while (data_mat(jlev,icur).le.bad_val.and.icur.ge.1)
                       icur=icur-1
                   if ( icur == 0 ) go to 10
              enddo
              10 continue
              itime_less=icur
              
              !find more time
              icur=itime
              do while (data_mat(jlev,icur).le.bad_val.and.icur.le.ntime)
                       icur=icur+1
                   if ( icur == ntime + 1 ) go to 11
              enddo
              11 continue
              itime_more=icur

              !correct limits
              if (itime_less.eq.0.and.itime_more.le.ntime) then
                 itime_less=itime_more
              end if
              if (itime_more.eq.(ntime+1).and.itime_less.ge.1) then
                 itime_more=itime_less
              end if
              if (itime_more.eq.(ntime+1).and.itime_less.eq.0) then
                 !whole timeseries consists of missing values
                 return
              end if

              !find out time weighting
              if (itime_more.ne.itime_less) then
              call get_time((TIME_VEC(itime)      - TIME_VEC(itime_less)),&
                        seconds_frac,days_frac)
              call get_time((TIME_VEC(itime_more) - TIME_VEC(itime_less)),&
                        seconds_int ,days_int )
              weight_more = real(seconds_frac+86400*days_frac) / &
                            real(seconds_int+86400*days_int)
              weight_less = 1.-weight_more
              else
              weight_more=1.
              weight_less=0.
              end if

              if (nlev .gt. 1) then


              !find less lev
              jcur=jlev
              do while (data_mat(jcur,itime).le.bad_val.and.jcur.ge.1)
                       jcur=jcur-1
              enddo
              jlev_less=jcur
              
              !find more time
              jcur=jlev
              do while (data_mat(jcur,itime).le.bad_val.and.jcur.le.nlev)
                       jcur=jcur+1
              enddo
              jlev_more=jcur

              !correct limits
              if (jlev_less.eq.0.and.jlev_more.le.nlev) then
                 jlev_less=jlev_more
              end if
              if (jlev_more.eq.(nlev+1).and.jlev_less.ge.1) then
                 jlev_more=jlev_less
              end if
              if (jlev_more.eq.(nlev+1).and.jlev_less.eq.0) then
                 if (mpp_pe() == mpp_root_pe()) then
                 print *, 'ERROR...jlev_more , jlev_less = ',&
                           jlev_more,jlev_less
                 end if
              end if
              
              end if  !for nlev .gt. 1

              !if interpolation in the vertical is over a smaller
              !number of points, then interpolate in the vertical
              !if it is not, then interpolate in time

              dpointslev = abs(jlev_more-jlev) + abs(jlev_less-jlev)
              dpointstime = abs(itime_more-itime) + abs(itime_less-itime)


              !produce answer
              if (dpointstime .le. dpointslev .or. nlev .eq. 1 .or. &
                  dpointslev  .ge. 2*(nlev+2)) then

              data_mat_new(jlev,itime)= &
                                   weight_more*data_mat(jlev,itime_more)+ &
                                   weight_less*data_mat(jlev,itime_less)              
              else

              weight_more = real(abs(jlev_less-jlev)) / real(dpointslev)
              weight_less = real(abs(jlev_more-jlev)) / real(dpointslev)
              
              data_mat_new(jlev,itime)= &
                                   weight_more*data_mat(jlev_more,itime)+ &
                                   weight_less*data_mat(jlev_less,itime)              
              
              
              end if

       end if  !for data being bad

enddo  !end do for itime
enddo  !end do for jlev

!assign to data_mat the new matrix
data_mat = data_mat_new

! 
!-----------------------------------------------------------------------
 
end subroutine rem_miss_var

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
     
subroutine twp_surface_flux_loop (                                     &
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
  real, dimension(size(t_atm(:))) :: t_ca_tmp, t_surf_tmp, q_surf_tmp, flux_lw_tmp
  real, dimension(size(t_atm(:))) :: rhotmp

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
                 dt,        land,       seawater,    avail             )


   !----------- initialize iteration variables  -----------------------!
   !-- Note that with specified surface fluxes, q_surf from the land
   !-- model may be garbage; thus initialize with 120% of q_atm
   t_ca_tmp   = t_surf
   t_surf_tmp   = t_surf
   q_surf_tmp = 1.2*q_atm
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

       call surface_flux (                                       &
           t_atm,     q_atm,      u_atm,       v_atm,            &
           p_atm,     z_atm,      p_surf,      t_surf_tmp,       &
           t_ca_tmp,  q_surf_tmp, u_surf,      v_surf,           &
           rough_mom, rough_heat, rough_moist, rough_scale, gust,&
           flux_t,    flux_q,     flux_lw,     flux_u,           &
           flux_v,    cd_m,       cd_t,        cd_q,             &
           w_atm,     u_star,     b_star,      q_star,           &
           dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
           dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
           dt,        land,       seawater,    avail             )        

       !------------------------------
       ! should inner loop be redone? 

       do_inner_loop = .false.
       do k = 1, npoints
         if ( ninner .lt. niter_max .and.                     &
              abs(flux_t(k)-SENFLUX) .gt. max_flux_error ) then

           do_inner_loop = .true.
           ! surface flux uses t_ca for land points and t_surf for non-land points.
           ! Therefore need to iterate on both.
           t_surf_tmp(k) = max( 200., min( 350.,             &
                t_surf_tmp(k) - ( (flux_t(k)-SENFLUX) /      &
                max(dhdt_surf(k),1.e-10) ) ) )
           t_ca_tmp(k) = max( 200., min( 350.,             &
                t_ca_tmp(k) - ( (flux_t(k)-SENFLUX) /      &
                max(dhdt_surf(k),1.e-10) ) ) )

         end if                            
       enddo        

       ninner = ninner + 1

     enddo  !------------- end of inner loop


     !---------------------------
     ! calculate surface flux

     call surface_flux (                                            &
              t_atm,     q_atm,      u_atm,       v_atm,            &
              p_atm,     z_atm,      p_surf,      t_surf_tmp,       &
              t_ca_tmp,  q_surf_tmp, u_surf,      v_surf,           &
              rough_mom, rough_heat, rough_moist, rough_scale, gust,&
              flux_t,    flux_q,     flux_lw,     flux_u,           &
              flux_v,    cd_m,       cd_t,        cd_q,             &
              w_atm,     u_star,     b_star,      q_star,           &
              dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
              dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
              dt,        land,       seawater,    avail             )     

     !------------------------------
     ! should outer loop be redone?

     do_outer_loop = .false.
     do k = 1, npoints
       if (  nouter .lt. niter_max .and.                           &
             hlv*abs(flux_q(k)- EVAPFLUX) .gt. max_flux_error ) then

         do_outer_loop = .true.
         q_surf_tmp(k) = max ( 1.e-08, min(0.05, q_surf_tmp(k)-&
                         ((flux_q(k)-EVAPFLUX) / max(-1.0*dedq_atm(k),1.e-10))))
!need to use dedq_atm(k) as dedq_surf=0 for non-land points.
!                         ((flux_q(k)-EVAPFLUX) / max(dedq_surf(k),1.e-10))))

       end if
     enddo

     nouter = nouter + 1

   enddo  !----------- end of outer loop

   !----------------------------
   ! deal with non-convergence
   
   rhotmp = p_atm/(rdgas*t_atm*(1.+d608*q_atm))  

   if ( nouter .gt. niter_max .or. ninner .gt. niter_max ) then

     do k = 1, npoints
     
       if ( hlv*abs(flux_q(k)- EVAPFLUX) .gt. max_flux_error .or.&
                abs(flux_t(k)- SENFLUX ) .gt. max_flux_error) then 
       
         if (mpp_pe() == mpp_root_pe()) print *, 'NON-CONVERGENCE '
                                        
         flux_t(k) = SENFLUX
         flux_q(k) = EVAPFLUX

         u_star(k) = sqrt ( sqrt ( flux_u(k)**2. +            &
                                   flux_v(k)**2. + 1.e-20 ) / rhotmp(k))
                                  
         q_star(k) = EVAPFLUX / rhotmp(k) / u_star(k)
         
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
         bflxs     = (chu*SENFLUX + cmu*EVAPFLUX) / rhotmp(k)
         b_star(k) = bflxs / u_star(k)

       end if  ! for outside of convergence bounds
          
     end do ! for flux errors

   end if  ! for lack of convergence

   if (mpp_pe() == mpp_root_pe()) then
   print *, 'desired fluxes  (sensible, latent): ',     &
        SENFLUX, hlv*EVAPFLUX
   print *, 'iterated fluxes (sensible, latent): ',     &
        flux_t, hlv*flux_q
   print *, 'iterated fluxes (t_ca_tmp, q_surf_tmp): ',     &
        t_ca_tmp, q_surf_tmp
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
           
end subroutine twp_surface_flux_loop

!#######################################################################
!########################################################################
! This subroutine returns surface conditions

subroutine get_twp_sfc( albedo_obs_out, tskin_out,     &
                        senflux_out, evapflux_out,     &
                        rough_mom_out, rough_heat_out )

  implicit none
  real, intent(out) :: albedo_obs_out, tskin_out,      &
                       senflux_out, evapflux_out,      &
                       rough_mom_out, rough_heat_out

  albedo_obs_out = ALBEDO_OBS
  tskin_out = TSKIN
  senflux_out = SENFLUX
  evapflux_out = EVAPFLUX
  rough_mom_out = ROUGH_MOM
  rough_heat_out = ROUGH_HEAT

 if (mpp_pe() == mpp_root_pe()) then

 print *, 'albedo_obs_out =', ALBEDO_OBS
 print *, ' tskin_out =', TSKIN
 print *, ' senflux_out =', SENFLUX
 print *, ' evapflux_out =', EVAPFLUX
 print *, ' rough_mom_out =', ROUGH_MOM
 print *, ' rough_heat_out =', ROUGH_HEAT

 endif

end subroutine get_twp_sfc

!########################################################################
! This subroutine returns imposed SST

subroutine get_twp_sst( sst )

  implicit none
  real, intent(out) :: sst

  sst = tskin_forc(1)

end subroutine get_twp_sst

!########################################################################


       function  get_column(field,ps_init,itime_less,itime_more,ivert,jstart,&
                            weight_less,weight_more, pfull, Xstand)
REAL, INTENT(IN), DIMENSION(:,:) :: field
REAL, INTENT(IN)                 :: ps_init
INTEGER, INTENT(IN)              :: itime_less, itime_more, ivert, jstart
REAL, INTENT(IN)                 :: weight_less, weight_more
REAL, INTENT(IN), DIMENSION(:) :: pfull
REAL, INTENT(IN), DIMENSION(:,:), OPTIONAL :: Xstand ! Standard atmosphere data.
REAL, DIMENSION(size(pfull,1)) :: get_column


REAL, DIMENSION(ivert+2+ivert_stand)                :: plev_int
REAL, DIMENSION(ivert+2+ivert_stand)                :: field_int_less, &
                                                       field_int_more
REAL, DIMENSION(SIZE(pfull,1))                      :: tmp_p,tmp_ans
INTEGER            :: nlev_interp, jstart_stand

! Set the lowest pressure level to the surface pressure 
       plev_int(1)=ps_init
! Set the lowest interpolation value to the value of the forcing field at the surface pressure.
       field_int_less(1)=field(jstart,itime_less)
       field_int_more(1)=field(jstart,itime_more)
! Set the interpolation pressure level to the forcing pressure levels
       plev_int(2:(ivert-jstart+2))=plev_forc(jstart:ivert)
! Set the interpolation forcing values to the value of the forcing field at the forcing pressure levels.
       field_int_less(2:(ivert-jstart+2))=field(jstart:ivert,itime_less)
       field_int_more(2:(ivert-jstart+2))=field(jstart:ivert,itime_more)
! Set a upper limit for the pressure levels.
       plev_int(ivert-jstart+3)=0.
! Set the interpolation forcing values to the forcing field at the top of the forcing data.
       field_int_less(ivert-jstart+3)=field(ivert,itime_less)
       field_int_more(ivert-jstart+3)=field(ivert,itime_more)
       nlev_interp = ivert-jstart+3


       if ( present(Xstand) ) then
! --- compute starting level for midlatitude reference profile  ------ !

         jstart_stand = 1
         do while (plev_forc(ivert).lt.plev_stand(jstart_stand,1))
             jstart_stand = jstart_stand + 1
         enddo       
         plev_int((ivert-jstart+3):                                      &
                  (ivert-jstart+2+ivert_stand-jstart_stand+1))=                    &
              plev_stand(jstart_stand:ivert_stand,1)
         field_int_less((ivert-jstart+3):                                &
                  (ivert-jstart+2+ivert_stand-jstart_stand+1))=          &
              Xstand(jstart_stand:ivert_stand,1)
         field_int_more((ivert-jstart+3):                                &
                  (ivert-jstart+2+ivert_stand-jstart_stand+1))=          &
              Xstand(jstart_stand:ivert_stand,1)
         nlev_interp = ivert-jstart+2+ivert_stand-jstart_stand+1
       endif
       
       tmp_p(:)=pfull(:) 
       ! --- do interpolation to the model resolution
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)
      
       get_column = tmp_ans
       end function get_column

!########################################################################
end module scm_twp_mod
