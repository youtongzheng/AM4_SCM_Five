module scm_mpaceb_mod
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MPACE PERIOD B FORCING MODULE
!
!       July 2006
!       Contact person: Steve Klein
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This module provides the large scale forcing necessary
!       to run the single column model to simulate MPACE PERIOD B.
!
!       The ARM Mixed Phase Arctic Cloud Experiment (MPACE) took place
!       in October 2004 at the North Slope of Alaska. Two periods were
!       selected for analysis.  
!
!       Period A :  14Z 5 October - 14Z 8 October, 2004
!       Period B :  17Z 9 October - 5Z 10 October, 2004
!
!
!       References
!       ----------
!
!       Klein, S. A., A. Fridlind, R. McCoy, G. McFarquhar, S. Menon, 
!       H. Morrison, D. Veron, S. Xie, J. J. Yio, and M. Zhang, 2006:
!       ARM Cloud Parameterization and Modeling Working Group - GCSS
!       Polar Cloud Working Group Model Intercomparison. Procedures for
!       ARM CPMWG Case 5/ GCSS Polar Cloud WG SCM/CRM/LES Intercomparison
!       Case f2004: ARM Mxied-Phase Arctic Cloud Experiment (M-PACE):
!       October 5-22, 2004. Available from:
!       http://science.arm.gov/wg/cpm/scimc5/docs/Document_ic5_v2.7.pdf
!      
!       Xie, S., S. A. Klein, J. J. Yio, A. C. M. Beljaars, C. N. Long, 
!       and M. Zhang (2006), An assessment of ECMWF analyses and model 
!       forecasts over the North Slope of Alaska using observations from 
!       the ARM Mixed-Phase Arctic Cloud Experiment, J. Geophys. Res., 
!       111, D05107, doi:10.1029/2005JD006509.        
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

   use sat_vapor_pres_mod, only:  lookup_es
   use            mpp_mod, only:  stdlog
   use         mpp_io_mod, only:  mpp_open,MPP_RDONLY,MPP_APPEND
   use            fms_mod, only:  write_version_number,     &
                                  mpp_pe, close_file, file_exist,   &
                                  check_nml_error, error_mesg, FATAL, &
                                  mpp_root_pe, open_restart_file,  & 
                                  read_data, write_data,           &
                                  nullify_domain, mpp_error, NOTE

#ifdef INTERNAL_FILE_NML
   USE              mpp_mod, ONLY: input_nml_file
#else
   USE              fms_mod, ONLY: open_namelist_file
#endif

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

   public mpaceb_data_read, mpaceb_forc_init, mpaceb_forc_end, update_mpaceb_forc, &
          mpaceb_forc_diagnostic_init,                                             &
          mpaceb_surface_flux_loop, get_mpaceb_flx, get_mpaceb_sst, ice_nucl_mpaceb

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following subroutines:
!
!            mpaceb_data_read   reads forcing data files and initializes
!                               any needed parameters of the run
!            mpaceb_forc_init   initializes prognostic variables:           
!                               T,u,v,qv,ql,qi
!            mpaceb_forc_end    deallocates allocated spaces
!            update_mpaceb_forc a call to this routine returns the 
!                               forcing parameters needed.
!            interp_time        finds the index of the TIME_TYPE whose
!                               time is just less than and just greater 
!                               than the time being interpolated to
!            interp_2d_field    interpolates a 2d field given as a 2 1d 
!                               vectors, assumed to be dimensioned over 
!                               pressure, to the desired time and level.
!            rem_miss_var       removes all missing data, indicated by 
!                               being less than a specified input value, 
!                               and replaces it with the value derived 
!                               from linear interpolation in time at a  
!                               given level
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              GLOBAL STORAGE VARIABLES
!
!  FORCING DATA (I.E. OBSERVED FIELDS)
!
!  ECMWF ATMOSPHERE PROFILES - TO BE USED FOR STRATOSPHERIC PROFILES
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
!
!  LAYERED FIELDS
!
!       t_forc         mean T of observed data (K)
!       q_forc         mean specific humidity of observed data 
!                      (kg vapor/kg air)
!       ql_forc        mean liquid water specific humidity of observed
!                      data (kg liquid/kg air)
!       N_forc         cloud drop number (1/m3)
!       u_forc         mean observed zonal wind (m/s)
!       v_forc         mean observed meridional wind (m/s)
!       ug_forc        geostrophic zonal wind (m/s)
!       vg_forc        geostrophic meridional wind (m/s) 
!
!  SURFACE DATA
!
!       ps_forc        surface pressure (Pa)
!       sst_forc       sea surface temperature (K)
!       div_forc       divergence of MBL air (s-1)
!       lat_forc       latitude of parcel (radians)
!       lon_forc       longitude of parcel (radians)
!       lwdn_forc      downward longwave radiation at 700 hPa (W/m2)
!       lwp_forc       liquid water path (kg/m2)
!       pib_forc       inversion base pressure (Pa)
!       cdn_forc       cloud drop # (1/m3)
!

character(len=8) :: mod_name = 'forcing'


REAL, ALLOCATABLE, DIMENSION(:)   :: plev_stand,t_stand,q_stand,O3_stand
REAL, ALLOCATABLE, DIMENSION(:)   :: plev_forc,ps_forc,sst_forc,div_forc
REAL, ALLOCATABLE, DIMENSION(:)   :: lon_forc,lat_forc,lwdn_forc
REAL, ALLOCATABLE, DIMENSION(:)   :: cdn_forc,lwp_forc,pib_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: t_forc,q_forc,ql_forc,u_forc,v_forc
REAL, ALLOCATABLE, DIMENSION(:,:) :: ug_forc, vg_forc, N_forc

TYPE(TIME_TYPE), ALLOCATABLE, DIMENSION(:)  :: time_layer

!----------Diagnostic data----------------------------------------------

integer :: id_temp_obs, id_qv_obs, id_ql_obs, id_omega_obs,id_uwnd_obs,&
           id_vwnd_obs, id_lwp_obs, id_pib_obs, id_pib_mod, &
           id_vadv_t, id_hadv_t, id_vadv_q, id_hadv_q, &
           id_vadv_ql,  id_vadv_qi
                         
! ---> h1g
integer ::  id_qadt_vadv,  id_qndt_vadv,     id_qnidt_vadv
! <--- h1g

!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       ivert_stand              # of vertical levels to standard atmos.
!       ivert                    # of vertical levels of mpace initial
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
!       wind_relaxation_tau      relaxation time scale for horizontal 
!                                winds (seconds)
!
!       p_omega_zero             pressure above which omega is set to
!                                to zero
!
!       p_cld_zero               pressure above which clouds are forced
!                                to be zero
!

INTEGER, PUBLIC                :: tracer_vert_advec_scheme = 3
INTEGER, PUBLIC                :: temp_vert_advec_scheme = 3
INTEGER, PRIVATE               :: itime = 1
INTEGER, PRIVATE               :: ivert = 550
INTEGER, PRIVATE               :: ivert_stand = 60
REAL,    PRIVATE               :: missing_value = -999.
REAL,    PRIVATE               :: p00 = 100000.
REAL,    PUBLIC                :: tskin = 274.01
REAL,    PUBLIC                :: SENFLUX = 136.5
REAL,    PUBLIC                :: EVAPFLUX = 107.7/hlv

real,    public                :: Ni_ini = 2.0e3   ! initial ice nuclei (#/m3)

REAL,    PUBLIC                :: p_omega_zero = 50000.
REAL,    PUBLIC                :: p_cld_zero = 50000.
REAL,    PUBLIC                :: wind_relaxation_tau = 3.*3600.
LOGICAL, PUBLIC                :: vert_advec_cond = .TRUE.
LOGICAL, PRIVATE               :: mpaceb_forc_initialized = .FALSE.
LOGICAL, PUBLIC                :: do_netcdf_restart = .true.

real, parameter   :: d622 = rdgas/rvgas
real, parameter   :: d378 = 1.0-d622
real, parameter   :: d608 = d378/d622

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!
        
character(len=128) :: Version = '$Id$'
character(len=128) :: Tagname = '$Name$'
        
NAMELIST /SCM_MPACEB_NML/ do_netcdf_restart,                           &
                          tracer_vert_advec_scheme,                    &
                          temp_vert_advec_scheme, vert_advec_cond,     &
                          wind_relaxation_tau,                         &
                          p_omega_zero, p_cld_zero, Ni_ini


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CONTAINS


!#######################################################################
!#######################################################################

subroutine mpaceb_data_read()

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
!       i,j,k,l,t            counting integers
!       unit                 unit number for I/O file
!       io,ierr              dummy integer variable
!       adum                 dummy character variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Internal variables
!  ------------------

INTEGER                                   :: i,j,k,l,t,unit,io,ierr
CHARACTER*23                              :: tracer_ascheme,temp_ascheme
CHARACTER*20                              :: adum
REAL                                      :: tmpr, tmprl
     
       ! if the constructor is already called, then return, otherwise 
       ! set mpaceb_forc_initialized .TRUE.

       if (mpaceb_forc_initialized) return
       mpaceb_forc_initialized = .TRUE.

!-----------------------------------------------------------------------
!      ----- read namelist -----
#ifdef INTERNAL_FILE_NML
       READ (input_nml_file, nml=SCM_MPACEB_NML, iostat=io)
       ierr = check_nml_error(io, 'SCM_MPACEB_NML')
#else
       if (file_exist('input.nml')) then
            unit = open_namelist_file()
            io=1
            do while (io .ne. 0)
                 read  (unit, nml=SCM_MPACEB_NML, iostat=io, end=10)
            enddo
  10        call close_file (unit)
       endif
#endif
 
       !dummy check -------         
       if (temp_vert_advec_scheme .gt. 6 ) then
            write (adum,'(i1)') temp_vert_advec_scheme
            call error_mesg ( 'scm_data_read in scm_mpaceb_mod',         &
                 'Bad temp_vert_advec_scheme : '//adum, FATAL )
       end if

       if (tracer_vert_advec_scheme .gt. 6 ) then
            write (adum,'(i1)') tracer_vert_advec_scheme
            call error_mesg ( 'scm_data_read in scm_mpaceb_mod',         &
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
            write (unit,'(a28,a23)') ' tracer_vert_advec_scheme = ', &
                                       tracer_ascheme
            write (unit,'(a26,a23)') ' temp_vert_advec_scheme = ', &
                                       temp_ascheme
            write (unit,*) ' vert_advec_cond  = ', vert_advec_cond      
            write (unit,*) ' wind relaxation time = ',                 &
                             wind_relaxation_tau
            write (unit,*) ' p_omega_zero = ', p_omega_zero
            write (unit,*) ' p_cld_zero = ',p_cld_zero
       endif
      
       call close_file (unit)
      
       
!-----------------------------------------------------------------------
!
!      ALLOCATE STORAGE
!

       if (ALLOCATED(plev_stand)) DEALLOCATE (plev_stand)
           ALLOCATE(plev_stand(ivert_stand))
       if (ALLOCATED(t_stand)) DEALLOCATE (t_stand)
           ALLOCATE(t_stand(ivert_stand))
       if (ALLOCATED(q_stand)) DEALLOCATE (q_stand)
           ALLOCATE(q_stand(ivert_stand))
       if (ALLOCATED(O3_stand)) DEALLOCATE (O3_stand)
           ALLOCATE(O3_stand(ivert_stand))
       
       if (ALLOCATED(plev_forc)) DEALLOCATE (plev_forc)
           ALLOCATE(plev_forc(ivert))
       if (ALLOCATED(time_layer)) DEALLOCATE (time_layer)
           ALLOCATE(time_layer(itime))
       
       if (ALLOCATED(t_forc)) DEALLOCATE (t_forc)
           ALLOCATE(t_forc(ivert,itime))
       if (ALLOCATED(q_forc)) DEALLOCATE (q_forc)
           ALLOCATE(q_forc(ivert,itime))
       if (ALLOCATED(ql_forc)) DEALLOCATE (ql_forc)
           ALLOCATE(ql_forc(ivert,itime))
       if (ALLOCATED(u_forc)) DEALLOCATE (u_forc)
           ALLOCATE(u_forc(ivert,itime))
       if (ALLOCATED(v_forc)) DEALLOCATE (v_forc)
           ALLOCATE(v_forc(ivert,itime))
       if (ALLOCATED(ug_forc)) DEALLOCATE (ug_forc)
           ALLOCATE(ug_forc(ivert,itime))
       if (ALLOCATED(vg_forc)) DEALLOCATE (vg_forc)
           ALLOCATE(vg_forc(ivert,itime))
       if (ALLOCATED(N_forc)) DEALLOCATE (N_forc)
           ALLOCATE(N_forc(ivert,itime))
       
       if (ALLOCATED(ps_forc)) DEALLOCATE (ps_forc)
           ALLOCATE(ps_forc(itime))
       if (ALLOCATED(sst_forc)) DEALLOCATE (sst_forc)
           ALLOCATE(sst_forc(itime))
       if (ALLOCATED(div_forc)) DEALLOCATE (div_forc)
           ALLOCATE(div_forc(itime))
       if (ALLOCATED(lat_forc)) DEALLOCATE (lat_forc)
           ALLOCATE(lat_forc(itime))
       if (ALLOCATED(lon_forc)) DEALLOCATE (lon_forc)
           ALLOCATE(lon_forc(itime))
       if (ALLOCATED(lwdn_forc)) DEALLOCATE (lwdn_forc)
           ALLOCATE(lwdn_forc(itime))
       if (ALLOCATED(cdn_forc)) DEALLOCATE (cdn_forc)
           ALLOCATE(cdn_forc(itime))
       if (ALLOCATED(lwp_forc)) DEALLOCATE (lwp_forc)
           ALLOCATE(lwp_forc(itime))
       if (ALLOCATED(pib_forc)) DEALLOCATE (pib_forc)
           ALLOCATE(pib_forc(itime))
       
!-----------------------------------------------------------------------
! 
!     READ IN ECMWF STRATOSPHERIC SOUNDING
!

       call mpp_open(unit,file='INPUT/arctic_profile.ascii',action=MPP_RDONLY)
      
       !skip to beginning of data
       do j=1,4
            read(unit,'(a)') adum
       enddo
            
       !read data
       do j=1,ivert_stand
            read (unit,'(f15.7,2X,f15.7,2X,f15.7,2X,f15.7)') &
                 plev_stand(j), t_stand(j),&
                 q_stand(j),O3_stand(j)  
       end do
       call close_file(unit)
       
       !----record which file was read
       unit= stdlog()
       if (mpp_pe() == 0 ) Write (unit,'(a)') 'used data file: INPUT/arctic_profile.ascii' 
       call close_file(unit)
       
       !change units
       
       !(a) mb to Pascals
       plev_stand(:)=100.*plev_stand(:)
             
       !(b) mass mixing ratio to specific humidity
       q_stand(:)=q_stand(:)/1000.
       q_stand(:)=q_stand(:)/(1.+q_stand(:))
       
       !(c) convert ppmv to specific humidity
       O3_stand(:)=(48.00/28.97)*1.e-06*O3_stand(:)
       O3_stand(:)=O3_stand(:)/(1.+O3_stand(:))
           
!-----------------------------------------------------------------------
! 
!      READ IN MPACE INITIAL SOUNDING
!
!
!       

       do 1000 t = 1, itime

            call mpp_open(unit,file='INPUT/formatted.mpaceB.initial.sounding',&
                          action=MPP_RDONLY)

            read(unit,*)
            read(unit,*)
            read(unit,*)
            read(unit,*)
            
            do k = 1, ivert
                 read(unit,13) plev_forc(k),u_forc(k,t),v_forc(k,t), &
                               t_forc(k,t),q_forc(k,t),ql_forc(k,t)
13               format(f8.3,1X,f5.1,1X,f5.1,1X,f7.3,1X,f7.5,1X,f7.5)
            end do
            
            call close_file(unit)  

            !----record which file was read
            call mpp_open(unit,'scm.out',action=MPP_APPEND)
            if (mpp_pe() == 0 ) Write (unit,'(a)')                  &
                 'used data file: INPUT/formatted.mpaceB.initial.sounding'
            call close_file (unit)

            !set time variable for VAR data 
            time_layer(t) = set_date(2004,10,9,17,0,0)
             
1000   continue

            
       ! units change
       ! ------------
       
       ! (a) convert (g/kg) to (kg/kg)
       q_forc(:,:)  =  q_forc(:,:)/1000.
       ql_forc(:,:) = ql_forc(:,:)/1000.

       ! (b) convert mixing ratio to specific humidity
       do t = 1, itime
       do k = 1, ivert
             tmpr  = q_forc(k,t)
             tmprl = ql_forc(k,t)
             q_forc(k,t)  = tmpr / (1. + tmpr + tmprl)
             ql_forc(k,t) = tmprl / (1. + tmpr + tmprl)
       enddo
       enddo
              
       ! (d) change pressure variables from mb to Pascals
       plev_forc(:) =  100.*plev_forc(:)
       ps_forc(:)   =  100.*  1010.00

       ! (f) convert 10-6 sec-1 to sec-1
       div_forc(:) = 5.8*1.E-06

       ! (g) convert latitude and longitude to radians
       lat_forc(:) = pi*71.75/180.
       lon_forc(:) = pi*(360.-151.0)/180.

       ! (h) convert cloud drop number to per cubic meter
       cdn_forc(:) = 20.*1.e+06
       N_forc(:,:) = 20.*1.e+06
       
       ! (i) compute lwp
       lwp_forc=0.
       do t = 1, itime
            do k = 2, ivert-1
                 tmpr=(plev_forc(k-1)+plev_forc(k  )) -          &
                      (plev_forc(k  )+plev_forc(k+1))
                 tmpr=tmpr*0.5
                 lwp_forc(t)  = lwp_forc(t) + ql_forc(k,t)*      &
                      tmpr/grav
            enddo
       enddo                 
            
       ! (j) set inversion pressure     
       pib_forc(:)=850.*100.
            
       ! (k) set SST
       sst_forc(:)=274.01
       
       ! (l) set lwdn_forc
       lwdn_forc(:) =  140.00
           
       !--------------------------------------
       ! reading tskin from the restart file.

       if( file_exist('INPUT/scm_mpaceb.res.nc') ) then

           if(mpp_pe() == mpp_root_pe() )                            &
           call mpp_error ('scm_mpaceb_mod',                           &
                           'Reading netCDF formatted restart file.', &
                            NOTE )
          call read_data('INPUT/scm_mpaceb.res.nc', 'SENFLUX',    SENFLUX   )
          call read_data('INPUT/scm_mpaceb.res.nc', 'EVAPFLUX',   EVAPFLUX  )
          call read_data('INPUT/scm_mpaceb.res.nc', 'Tskin', tskin )

       else if (file_exist('INPUT/scm_mpaceb.res')) then

          unit=open_restart_file('INPUT/scm_mpaceb.res',action='read')
          read(unit)  SENFLUX, EVAPFLUX, tskin
          call close_file(unit)

       endif       
        

!-----------------------------------------------------------------------
        
end subroutine mpaceb_data_read


!#######################################################################
!#######################################################################


subroutine mpaceb_forc_init(time_interp,pdamp, As, Bs)
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
!      pdamp           maximum pressure for which profile is set to the
!                      ECMWF reference profile  <<< NOTE THAT THIS IS IGNORED
!      As, Bs          A's and B's of half levels in hybrid coordinate
!                      phalf(k) = A(k) + B(k) * psurf
!
!      -------------
!      INPUT/OUTPUT:
!      -------------
!
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
REAL,  INTENT (IN)                       :: pdamp
REAL,  INTENT (IN)   , DIMENSION(:)      :: As,Bs

!  Internal variables
!  ------------------

INTEGER  :: itime_less, itime_more, KDIM, k
INTEGER  :: nlev_interp, days, months, years, seconds, minutes, hours
REAL     :: weight_less, weight_more


REAL, DIMENSION(ivert+1)                         :: plev_int
REAL, DIMENSION(ivert+1)                         :: field_int_less, &
                                                    field_int_more
REAL, DIMENSION(size(pt,3))                       :: tmp_p,tmp_ans
REAL, DIMENSION(size(pt,1),size(pt,2))             :: ps_init
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3))   :: pfull,gamma,qsat
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3)+1) :: phalf,lphalf

integer :: i,j
#include "fv_point.inc"
       nsphum = get_tracer_index(MODEL_ATMOS, 'sphum')
       nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
       nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
       nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
       nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )   !h1g
       nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )   !h1g

!
!CODE
!
       
       
       ! --- find out # of vertical levels
       KDIM = size(pt,3)

       ! --- set all itime_less, itime_more, weight_less, weight_more
       
       itime_less = 1
       itime_more = 1
       weight_less = 1.
       weight_more = 0.

! --- CREATE PFULL AND PHALF  ---------------------------------------- !

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
       
! --- CREATE INITIAL pt  ---------------------------------------------- !
            
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=  t_forc(1:ivert,itime_less)
       plev_int(ivert+1)=0.
       field_int_less(ivert+1)=t_forc(ivert,itime_less)
       field_int_more = field_int_less
       tmp_p(:)=pfull(1,1,:)
       
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int,field_int_less,     &
            field_int_more,tmp_ans)

       !---- assign T
       DO k = 1,KDIM
              pt(:,:,k) = tmp_ans(k)
       ENDDO
      
! --- CREATE INITIAL Qv  --------------------------------------------- !

       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=  q_forc(1:ivert,itime_less)
       plev_int(ivert+1)=0.
       field_int_less(ivert+1)=q_forc(ivert,itime_less)
       field_int_more = field_int_less
       tmp_p(:)=pfull(1,1,:)
       
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int,field_int_less,     &
            field_int_more,tmp_ans)
      
       !---- assign answer to qv
       DO k = 1,KDIM
              q(:,:,k,nsphum) = tmp_ans(k)
       ENDDO
       
! --- CREATE INITIAL QL  --------------------------------------------- !
       
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=  ql_forc(1:ivert,itime_less)
       plev_int(ivert+1)=0.
       field_int_less(ivert+1)=ql_forc(ivert,itime_less)
       field_int_more = field_int_less
       tmp_p(:)=pfull(1,1,:)
       
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int,field_int_less,     &
            field_int_more,tmp_ans)
      
       !---- assign answer to ql
       DO k = 1,KDIM
              q(:,:,k,nql) = tmp_ans(k)
       ENDDO

! --- CREATE INITIAL QI,QA  ------------------------------------------ !

       q(:,:,:,nqi)  =  0.
       WHERE (q(:,:,:,nql) .gt. 0.) 
              q(:,:,:,nqa)  =  1.
       END WHERE

! ---> h1g
       if( nqn  > 0)   q(:,:,:,nqn)  = 0.0   !h1g
       if( nqni > 0)   q(:,:,:,nqni) = 0.0   !h1g
! <--- h1g

! --- CREATE INITIAL U,V  -------------------------------------------- !
      
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=  u_forc(1:ivert,itime_less)
       plev_int(ivert+1)=0.
       field_int_less(ivert+1)=u_forc(ivert,itime_less)
       field_int_more = field_int_less
       tmp_p(:)=pfull(1,1,:)
       
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int,field_int_less,     &
            field_int_more,tmp_ans)
      
       !---- assign answer to u
       DO k = 1,KDIM
              ua(:,:,k) = tmp_ans(k)
       ENDDO
       u_srf(:,:)=ua(:,:,KDIM)

       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=  v_forc(1:ivert,itime_less)
       plev_int(ivert+1)=0.
       field_int_less(ivert+1)=v_forc(ivert,itime_less)
       field_int_more = field_int_less
       tmp_p(:)=pfull(1,1,:)
       
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int,field_int_less,     &
            field_int_more,tmp_ans)
      
       !---- assign answer to v
       DO k = 1,KDIM
              va(:,:,k) = tmp_ans(k)
       ENDDO
       v_srf(:,:)=va(:,:,KDIM)
       
! ---- interpolate sst ----------------------------------------------- !
     
       !---- interpolate sst field
       tmp_p(:) = 50000.
       plev_int(:) = p00
       field_int_less(:) = sst_forc(itime_less)
       field_int_more(:) = sst_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
       tskin = tmp_ans(1)
 
! 
!-----------------------------------------------------------------------
 
 end subroutine mpaceb_forc_init

!#######################################################################
!#######################################################################


subroutine mpaceb_forc_end ()

  integer :: unit

  if (.NOT.mpaceb_forc_initialized) return
  mpaceb_forc_initialized = .FALSE.

  DEALLOCATE (  plev_stand, t_stand, q_stand, O3_stand )
  DEALLOCATE (  plev_forc , t_forc, q_forc, ql_forc, u_forc, v_forc )
  DEALLOCATE (  ug_forc, vg_forc, N_forc )
  DEALLOCATE (  lat_forc, lon_forc, ps_forc, sst_forc, div_forc )
  DEALLOCATE (  lwdn_forc, cdn_forc, time_layer )
  
  !--------------------------------------
  ! write tskin to the restart file.

  if( do_netcdf_restart ) then
 
     if (mpp_pe() == mpp_root_pe()) &
        call write_data ('RESTART/scm_mpaceb.res.nc', 'SENFLUX',    SENFLUX   )
        call write_data ('RESTART/scm_mpaceb.res.nc', 'EVAPFLUX',   EVAPFLUX  )
        call write_data( 'RESTART/scm_mpaceb.res.nc', 'Tskin', tskin )
 
  else
        
     unit=open_restart_file('RESTART/scm_mpaceb.res',action='write')
     if(mpp_pe() == mpp_root_pe())  write(unit)  SENFLUX, EVAPFLUX, tskin
     call close_file(unit)       

  endif        
! 
!-----------------------------------------------------------------------
 

end subroutine mpaceb_forc_end


!#######################################################################
!#######################################################################

subroutine mpaceb_forc_diagnostic_init(axes, Time)

implicit none

integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time
   
! --- initialize axes -------------------------------------------------!

       id_temp_obs = register_diag_field (mod_name, 'temp_obs',        &
          axes(1:3), Time, 'Observed Temperature', 'K',                &
          missing_value = missing_value )
       id_qv_obs   = register_diag_field (mod_name, 'qv_obs',          &
          axes(1:3), Time, 'Observed Water Vapor Specific Humidity',   &
          'g/kg', missing_value = missing_value )
       id_ql_obs   = register_diag_field (mod_name, 'ql_obs',          &
          axes(1:3), Time, 'Observed Water Vapor Specific Humidity',   &
          'g/kg', missing_value = missing_value )
       id_omega_obs = register_diag_field (mod_name, 'omega_obs',      &
          axes(1:3),  Time, 'Observed Pressure Velocity', 'Pa/s',      &
          missing_value = missing_value )
       id_uwnd_obs = register_diag_field (mod_name, 'uwnd_obs',        &
          axes(1:3), Time, 'Observed Zonal Wind', 'm/s',               &
          missing_value = missing_value )
       id_vwnd_obs = register_diag_field (mod_name, 'vwnd_obs',        &
          axes(1:3), Time, 'Observed Meridional Wind', 'm/s',          &
         missing_value = missing_value )
       
       id_lwp_obs = register_diag_field (mod_name, 'lwp_obs',          &
          axes(1:2), Time, 'Observed Liquid Water Path', 'kg/m2',      &
          missing_value = missing_value )
       id_pib_obs = register_diag_field (mod_name, 'pib_obs',          &
          axes(1:2), Time, 'Observed Inversion Base', 'Pa',            &
          missing_value = missing_value )
       id_pib_mod = register_diag_field (mod_name, 'pib_mod',          &
          axes(1:2), Time, 'Model Inversion Base', 'Pa',               &
          missing_value = missing_value )
       
       id_vadv_t = register_diag_field (mod_name, 'tdt_vadv',          &
          axes(1:3), Time,                                             &
          'Temperature tendency from Vertical Advection', 'K/sec',     &
          missing_value = missing_value )
          
       id_vadv_q = register_diag_field (mod_name, 'qdt_vadv',          &
          axes(1:3), Time,                                             &
          'Water vapor tendency from Vertical Advection', 'kg/kg/sec', &
          missing_value = missing_value )
       
       id_hadv_t = register_diag_field (mod_name, 'tdt_hadv',          &
          axes(1:3), Time,                                             &
          'Temperature tendency from Horizontal Advection', 'K/sec',   &
          missing_value = missing_value )
          
       id_hadv_q = register_diag_field (mod_name, 'qdt_hadv',          &
          axes(1:3), Time,                                             &
          'Water vapor tendency from Horizontal Advection','kg/kg/sec',&
          missing_value = missing_value )
          
       id_vadv_ql = register_diag_field (mod_name, 'qldt_vadv',        &
          axes(1:3), Time,                                             &
          'Liquid Water tendency from Vertical Advection', 'kg/kg/sec',&
          missing_value = missing_value )
       
       id_vadv_qi = register_diag_field (mod_name, 'qidt_vadv',        &
          axes(1:3), Time,                                             &
          'Ice Water tendency from Vertical Advection', 'kg/kg/sec',   &
          missing_value = missing_value )

       id_qadt_vadv = register_diag_field (mod_name, 'qadt_vadv', axes(1:3), Time, &
          'cloud amount tendencies due to vertical advection', '1/s', missing_value = missing_value)

       id_qndt_vadv = register_diag_field (mod_name, 'qndt_vadv', axes(1:3), Time, &
          'liquid droplet number concentration tendencies due to vertical advection', '1/cm3/s', missing_value = missing_value)

       id_qnidt_vadv = register_diag_field (mod_name, 'qnidt_vadv', axes(1:3), Time, &
          'ice number concentration tendencies due to vertical advection', '1/cm3/s', missing_value = missing_value)
       
       
                                                          
! 
!-----------------------------------------------------------------------
 
end subroutine mpaceb_forc_diagnostic_init


!#######################################################################
!#######################################################################


subroutine update_mpaceb_forc(time_interp,time_diag,dt_int,pdamp )
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
!      pdamp        pressure for which all tendencies at p < pdamp will
!                   be set to zero. 
!      As, Bs       A's and B's of half levels in hybrid coordinate
!                   phalf(k) = A(k) + B(k) * psurf
!
!      ------------
!      INPUT/OUTPUT
!      ------------
!
!      pfull        pressure in PASCALS at full model levels
!      phalf        pressure in PASCALS at half model levels
!                   NOTE that it is assumed that p(j)<p(j+1)
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
!
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
!      omega_full       omega interpolated to full model levels 
!                       (Pa/sec)
!
!
!      -------------------
!      INTERNAL VARIABLES:
!      -------------------
!
!      j,k              counting integers
!      itime_more       indicates indice of time_layer for which 
!                       time_layer(itime_more)>= time_interp
!      itime_less       indicates indice of time_layer for which 
!                       time_layer(itime_more)<= time_interp
!      KDIM             # of vertical levels to T array
!      weight_more      real number indicating closeness of inter-
!                       polated time to time_layer(itime_more)
!      weight_less      real number indicating closeness of inter-
!                       polated time to time_layer(itime_less)
!      plev_int         array of pressure coordinates to be passed to 
!                       interp_2d_field (Pa)
!      field_int_less   1D array of field at less time
!      field_int_more   1D array of field at more time 
!      tmp_p            temporary pressure levels for interpolating
!      tmp_ans          temporary array of answer from interpolating
!      dt_seconds       time step (seconds) 
!      dt_days          time step (days)
!      ps_init          surface pressure initial value (Pa)
!      omega_half       pressure velocity at model half levels (Pa/sec)
!      adv_vert         vertical advection of a quantity (units/second)
!      dp               pressure thickness of model levels (Pa)
!      t_obs            observed temperature (K)
!      q_obs            observed qv (kg water/kg air)
!      ql_obs           observed liquid water specific humidity 
!                       (kg condensate/kg air)
!      u_obs            observed zonal wind (m/s)
!      v_obs            observed meridional wind (m/s)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

TYPE(TIME_TYPE), intent(in)              :: time_interp,time_diag,dt_int
REAL,  INTENT (IN)                       :: pdamp


!  Internal variables
!  ------------------

INTEGER                                          :: days,months,years
INTEGER                                          :: seconds,minutes
INTEGER                                          :: hours
INTEGER                                          :: i,j,k,itime_less
INTEGER                                          :: itime_more,KDIM
INTEGER                                          :: dt_seconds,dt_days
LOGICAL                                          :: used
REAL                                             :: tmplwp, tmplwpforc
REAL                                             :: weight_less
REAL                                             :: weight_more
REAL, DIMENSION(ivert+1)                         :: plev_int
REAL, DIMENSION(ivert+1)                         :: field_int_less
REAL, DIMENSION(ivert+1)                         :: field_int_more
REAL, DIMENSION(size(pt,3))                       :: tmp_p,tmp_ans
REAL, DIMENSION(size(pt,1),size(pt,2))             :: ps_init, tmp_res2
REAL, DIMENSION(size(pt,1),size(pt,2))             :: pinv_mod,pinv_obs
REAL, DIMENSION(size(pt,1),size(pt,2))             :: tmpr
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3))   :: pfull,t_obs,q_obs
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3))   :: ql_obs,u_obs,v_obs
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3))   :: pi_fac, dp
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3))   :: adv_hor
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3))   :: adv_vert,tmp_res
REAL, DIMENSION(size(pt,1),size(pt,2),size(pt,3)+1) :: phalf, omega_half

#include "fv_point.inc"

!CODE
!

! --- find out # of vertical levels
KDIM = size(pt,3)

! --- compute timestep interval in seconds and days
call get_time(dt_int,dt_seconds,dt_days)
dt_seconds = dt_seconds + 86400*dt_days

! --- UPDATE PFULL AND PHALF, COMPUTE OMEGA SURFACE, AND PI FACTOR --- !

       ! --- find out time in sfc pressure field ---- !
       call interp_time(time_interp,time_layer,itime_less,itime_more,  &
            weight_less,weight_more)
 
       itime_less = 1
       itime_more = 1
       weight_less = 1.
       weight_more = 0.

       !---- interpolate surface pressure field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=ps_forc(itime_less)
       field_int_more(:)=field_int_less
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       ps_init(:,:) = tmp_ans(1)
       ps=ps_init
       
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pfull(i,j,:), phalf(i,j,:))
         enddo
       enddo
! --- ZERO OUT UPPER LEVEL CLOUDS ------------------------------------ !

       where (pfull .lt. p_cld_zero)
         pt = pt - hlv*q(:,:,:,nql)/cp_air - hls*q(:,:,:,nqi)/cp_air
         q(:,:,:,nsphum) = q(:,:,:,nsphum) + q(:,:,:,nql) + q(:,:,:,nqi)
         q(:,:,:,nql) = 0.0
         q(:,:,:,nqi) = 0.0
         q(:,:,:,nqa) = 0.0
       endwhere

! ----------------- COMPUTE INVERSION LEVEL -------------------------- !

       !---- interpolate pib field
       tmp_p(:) = 50000.
       plev_int(:) = p00
       field_int_less(:) = pib_forc(itime_less)
       field_int_more(:) = field_int_less
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       pinv_obs(:,:) = tmp_ans(1)
       
       if (id_pib_obs > 0) then
       used = send_data( id_pib_obs, pinv_obs(:,:), time_diag, 1, 1)
       end if  

       !---- find model inversion level as level with temperature
       !---- minimum.
       !---- Note that the model definition uses the half level above the
       !---- temperature minimum.  This is consistent with the sub-grid 
       !---- model for temperature.  
       
       tmpr = 0.
       pinv_mod = pinv_obs
       do k = KDIM-1,2,-1
            where (pt(:,:,k) .lt. pt(:,:,k-1) .and.  &
                   pt(:,:,k) .lt. pt(:,:,k+1) .and.  &
                   abs(phalf(:,:,k)-pinv_obs).lt.10000.)
                  pinv_mod = phalf(:,:,k)
            end where
       enddo

       if (id_pib_mod > 0) then
       used = send_data( id_pib_mod, pinv_mod(:,:), time_diag, 1, 1)
       end if  

       if (minval(pinv_obs) .lt. p_omega_zero)                         &
            call error_mesg ('update_mpaceb_forc in scm_mpaceb_mod',   &
            'pinv_obs < p_omega_zero', FATAL )
            
! ----------------- COMPUTE OMEGA ON HALF AND FULL LEVELS ------------ !

       !---- interpolate surface divergence field
       tmp_p(:) = 50000.
       plev_int(:) = p00
       field_int_less(:) = div_forc(itime_less)
       field_int_more(:) = field_int_more
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)

       !---- compute omega_half
       omega_half = 0.0
       do k = 1, KDIM + 1
            where(phalf(:,:,k) .gt. pinv_obs(:,:))
                 omega_half(:,:,k) = tmp_ans(1)*(ps_init(:,:) -        &
                      phalf(:,:,k))
            endwhere
            where(phalf(:,:,k) .le. pinv_obs(:,:))
                 omega_half(:,:,k) = tmp_ans(1)*( ps_init(:,:) -       &
                                                 pinv_obs(:,:) )
            endwhere 
       enddo
       
       !---- compute omega_full
       omga = 0.0
       do k = 1, KDIM
            where(pfull(:,:,k) .gt. pinv_obs(:,:))
                 omga(:,:,k) = tmp_ans(1)*(ps_init(:,:) -        &
                      pfull(:,:,k))
            endwhere
            where(pfull(:,:,k) .le. pinv_obs(:,:))
                 omga(:,:,k) = tmp_ans(1)*( ps_init(:,:) -       &
                                           pinv_obs(:,:) )
            endwhere 
       enddo        

! ----------------- PRODUCE TEMPERATURE FORCING ---------------------- !

       !---- initialize variables
       t_dt     = 0.
       adv_vert = 0.       
       adv_hor  = 0.
       
       !--------- COMPUTE ADIABATIC COMPRESSION  --------------------- !
      
       !---- compute adiabatic compression
       t_dt(:,:,:)=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pfull(:,:,:)

       SELECT CASE (temp_vert_advec_scheme)
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
       call vert_advection(REAL(dt_seconds),omega_half,delp,pt, &
            adv_vert,scheme = vadvec_scheme, form = ADVECTIVE_FORM)
       
       !---- COMPUTE HORIZONTAL T ADVECTION -------------------------- !
       
       do k = 1, KDIM   
            tmpr(:,:) = 0.0
            tmpr(:,:) = (ps_init(:,:)-pfull(:,:,k))/21818.
            tmpr(:,:) = MIN(-4.,-15.*(1.-tmpr(:,:)))/86400.
            adv_hor(:,:,k) = tmpr(:,:)
       enddo     
       
       !---- set tendencies to zero above pdamp
       do k=1,KDIM          
            where(pfull(:,:,k).le.pdamp)
                 t_dt(:,:,k)       = 0.0
                 adv_vert(:,:,k) = 0.0
                 adv_hor(:,:,k)  = 0.0
            endwhere
       enddo
                       
       !---- interpolate t_forc
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=t_forc(1:ivert,itime_less)
       field_int_more(1:ivert)=field_int_less(1:ivert)
       plev_int(ivert+1) = 0.
       field_int_less(ivert+1) = t_forc(ivert,itime_less)
       field_int_more(ivert+1) = field_int_less(ivert+1)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+1)),field_int_less(1:(ivert+1)),         &
            field_int_more(1:(ivert+1)),tmp_ans)
       
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result
       t_obs(1,1,:) = tmp_ans
       if (id_temp_obs > 0) then
       used = send_data ( id_temp_obs, t_obs(:,:,:), time_diag, 1, 1, 1)
       end if  
       
       !------ tendency diagnostics
       if (id_vadv_t > 0) then
            used = send_data ( id_vadv_t, adv_vert+t_dt, time_diag, 1, 1, 1)
       end if
 
       if (id_hadv_t > 0) then
            used = send_data ( id_hadv_t, adv_hor, time_diag, 1, 1, 1)
       end if

       !---- COMPUTE FINAL ADVECTIVE T TENDENCY FROM THE SUM OF :  --- !
       !---- RELAXATION TENDENCY, VERTICAL ADVECTION AND ADIABATIC --- !
       !---- COMPRESSION                                           --- !
       t_dt(:,:,:) = t_dt(:,:,:) + adv_vert(:,:,:) + adv_hor(:,:,:) 
      
 
! ----------------- PRODUCE MOISTURE FORCING ------------------------- !

       !---- initialize variables
       q_dt(:,:,:,nsphum)       = 0.
       adv_vert = 0.       
       adv_hor  = 0.
              
       !---- COMPUTE VERTICAL Q ADVECTION ---------------------------- !
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

       call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nsphum),&
            adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)

       !---- COMPUTE HORIZONTAL Q ADVECTION ------------------- !
       
       do k = 1, KDIM   
            tmpr(:,:) = 0.0
            tmpr(:,:) = (ps_init(:,:)-pfull(:,:,k))/15171.
            tmpr(:,:) = MIN(0.164,-3.*(1.-tmpr(:,:)))/1000./86400.
            adv_hor(:,:,k) = tmpr(:,:)
       enddo     
       
       !---- Zero out advection terms above pdamp
       do k=1,KDIM          
            where(pfull(:,:,k).le.pdamp)
                 adv_vert(:,:,k) = 0.0
                 adv_hor(:,:,k)  = 0.0
            endwhere
       enddo
                   
       !---- interpolate q_forc
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)= q_forc(1:ivert,itime_less)
       field_int_more(1:ivert)=field_int_less(1:ivert)
       plev_int(ivert+1) = 0.
       field_int_less(ivert+1) = q_forc(ivert,itime_less)
       field_int_more(ivert+1) = field_int_less(ivert+1)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+1)),field_int_less(1:(ivert+1)),         &
            field_int_more(1:(ivert+1)),tmp_ans)
       
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result (convert kg/kg to g/kg)
       q_obs(1,1,:) = tmp_ans(:)
       if (id_qv_obs > 0) then
            used = send_data( id_qv_obs, 1000.*q_obs(:,:,:), time_diag,&
                 1, 1, 1) 
       end if  
      
       !------ tendency diagnostics
       if (id_vadv_q > 0) then
            used = send_data ( id_vadv_q, adv_vert, time_diag, 1, 1, 1)
       end if
 
       if (id_hadv_q > 0) then
            used = send_data ( id_hadv_q, adv_hor, time_diag, 1, 1, 1)
       end if

       !---- COMPUTE FINAL ADVECTIVE Q TENDENCY FROM THE SUM OF :  --- !
       !---- RELAXATION TENDENCY AND VERTICAL ADVECTION            --- !                                      --- !
       q_dt(:,:,:,nsphum) = adv_vert(:,:,:) + adv_hor(:,:,:)
            
! ----------------- PRODUCE MOMENTUM FORCING ------------------------- !

       !---- initialize variables
       u_dt     = 0.
       v_dt     = 0.     
       
       !---- interpolate u_forc
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=u_forc(1:ivert,itime_less)
       field_int_more(1:ivert)=field_int_less(1:ivert)
       plev_int(ivert+1) = 0.
       field_int_less(ivert+1) = u_forc(ivert,itime_less)
       field_int_more(ivert+1) = field_int_less(ivert+1)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+1)),field_int_less(1:(ivert+1)),         &
            field_int_more(1:(ivert+1)),tmp_ans)
       
       !---- compute damping term
       do k=1,KDIM          
            u_dt(:,:,k)=(tmp_ans(k)-ua(:,:,k))/  &
                 max(wind_relaxation_tau,REAL(dt_seconds))
       enddo
      
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result 
       u_obs(1,1,:) = tmp_ans(:)
       if (id_uwnd_obs > 0) then
       used = send_data( id_uwnd_obs, u_obs(:,:,:), time_diag, 1, 1, 1)
       end if  
      
       !---- interpolate u_forc
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=v_forc(1:ivert,itime_less)
       field_int_more(1:ivert)=field_int_less(1:ivert)
       plev_int(ivert+1) = 0.
       field_int_less(ivert+1) = v_forc(ivert,itime_less)
       field_int_more(ivert+1) = field_int_less(ivert+1)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+1)),field_int_less(1:(ivert+1)),         &
            field_int_more(1:(ivert+1)),tmp_ans)
       
       !---- compute damping term
       do k=1,KDIM
            v_dt(:,:,k)=(tmp_ans(k)-va(:,:,k))/ &
                 max(wind_relaxation_tau,REAL(dt_seconds))
       enddo
       
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result 
       v_obs(1,1,:) = tmp_ans(:)
       if (id_vwnd_obs > 0) then
       used = send_data( id_vwnd_obs, v_obs(:,:,:), time_diag, 1, 1, 1)
       end if  
      
       
! ----------------- PRODUCE CONDENSATE FORCING ----------------------- !

       q_dt(:,:,:,nql)=0.
       q_dt(:,:,:,nqi)=0. 
       q_dt(:,:,:,nqa)=0.
       if(nqn > 0) q_dt(:,:,:,nqn) = 0.0

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

          ! Vertical advection of cloud liquid
          call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nql),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
          q_dt(:,:,:,nql) = q_dt(:,:,:,nql) + adv_vert
          
          if ( id_vadv_ql > 0) then
            used = send_data ( id_vadv_ql,  q_dt(:,:,:,nql), time_diag, 1, 1, 1)
          end if

          ! Vertical advection of cloud ice
          call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nqi),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
          q_dt(:,:,:,nqi) = q_dt(:,:,:,nqi) + adv_vert
          
          if ( id_vadv_qi > 0) then
            used = send_data ( id_vadv_qi,  q_dt(:,:,:,nqi), time_diag, 1, 1, 1)
          end if

          ! Vertical advection of cloud area
          call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nqa),&
               adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
          q_dt(:,:,:,nqa) = q_dt(:,:,:,nqa) + adv_vert
          
          if ( id_qadt_vadv > 0) then
            used = send_data ( id_qadt_vadv,  q_dt(:,:,:,nqa), time_diag, 1, 1, 1)
          end if

          ! Vertical advection of cloud droplet number
          if ( nqn > 0 ) then
            call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nqn),&
                 adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
            q_dt(:,:,:,nqn) = q_dt(:,:,:,nqn) + adv_vert
           
            if ( id_qndt_vadv > 0) then
              used = send_data ( id_qndt_vadv,  q_dt(:,:,:,nqn), time_diag, 1, 1, 1)
            end if
          endif

          ! Vertical advection of cloud ice number
          if ( nqni > 0 ) then
            call vert_advection(REAL(dt_seconds),omega_half,delp,q(:,:,:,nqni),&
                 adv_vert, scheme = vadvec_scheme, form = ADVECTIVE_FORM)
            q_dt(:,:,:,nqni) = q_dt(:,:,:,nqni) + adv_vert
          
            if ( id_qnidt_vadv > 0) then
              used = send_data ( id_qnidt_vadv,  q_dt(:,:,:,nqni), time_diag, 1, 1, 1)
            end if
          endif

       end if  !  vert_advec_cond


! ------------------- COMPUTE SURFACE DATA:  sst --------------------- !

       !---- interpolate sst field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=sst_forc(itime_less)
       field_int_more(:)=field_int_less(:)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
                     field_int_less,field_int_more,tmp_ans)
       tskin = tmp_ans(1)
        
! ----------------- do NETCDF DIAGNOSTIC OUTPUT ---------------------- !

       !---- extract omega_forc
       tmp_res(1,1,:) = 86400.*omga(1,1,:)/100.
       if (id_omega_obs > 0)                                           &
       used = send_data(id_omega_obs, tmp_res(:,:,:), time_diag, 1, 1,1)  

       !---- interpolate lwp field
       tmp_p(:) = 50000.
       plev_int(:) = p00
       field_int_less(:) = lwp_forc(itime_less)
       field_int_more(:) = field_int_less(:)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       tmp_res2(:,:) = tmp_ans(1)       
       tmplwpforc = tmp_res2(1,1)
       
       if (id_lwp_obs > 0) then
       used = send_data( id_lwp_obs, tmp_res2(:,:), time_diag, 1, 1)
       end if  
          
       ! ---- interpolate ql ---------------------------------- !
       
       plev_int(1:ivert)=plev_forc(1:ivert)
       field_int_less(1:ivert)=  ql_forc(1:ivert,itime_less)
       field_int_more = field_int_less
       plev_int(ivert+1) = 0.
       field_int_less(ivert+1) = ql_forc(ivert,itime_less)
       field_int_more(ivert+1) = field_int_less(ivert+1)
       tmp_p(:)=pfull(1,1,:)
       
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more, &
            plev_int(1:ivert),field_int_less(1:ivert),     &
            field_int_more(1:ivert),tmp_ans)
                          
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
              
       !---- output netcdf result (convert kg/kg to g/kg)
       ql_obs(1,1,:) = 1000.*tmp_ans(:)
       if (id_ql_obs > 0)                                              &
       used = send_data( id_ql_obs, ql_obs(:,:,:), time_diag, 1, 1, 1)
                                                                           
! 
!-----------------------------------------------------------------------
 
end subroutine update_mpaceb_forc


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
!      THE subroutine ASSUMES THAT TIME_VEC IS ARRANGED IN ASCENDNG ORDER
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
itime_vec = size(TIME_VEC,1)

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

REAL,    DIMENSION(size(field_less,1))  :: field_tmp
REAL                                    :: ptmp
INTEGER, DIMENSION(size(p_interp,1))    :: klev_more,klev_less
INTEGER                                 :: j,nfield_lev,ninterp
LOGICAL, DIMENSION(size(p_interp,1))    :: FIND_PINTERP
REAL,    DIMENSION(size(p_interp,1))    :: pweight_more,pweight_less

!
! Code
! ----

! ---  count # of levels to the field and # of interpolating levels
nfield_lev = size(field_less,1) 
ninterp = size(p_interp,1)


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
REAL, DIMENSION(size(data_mat,1),size(data_mat,2)) :: data_mat_new
!CODE
!


! --- find out # of vertical levels
nlev = size(data_mat,1)
if (present(adjlev)) nlev=adjlev
ntime = size(data_mat,2)


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
                 print *, 'ERROR...jlev_more , jlev_less = ',&
                           jlev_more,jlev_less
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



!########################################################################
!########################################################################


!=======================================================================
!
!      These subroutines contain the follow interfaces:
!
!   vert_advec_2nd   2nd order scheme with Euler-backward time
!                    differencing
!   vert_advec_4th   4th order scheme with Euler-backward time
!                    differencing
!
!-----------------------------------------------------------------------
!
!    The vertical grid indexing
!
!                        ---------------
!               deta(1)      var(1)
!                        ---------------  etadot(1)
!               deta(2)      var(2)
!                        ---------------  etadot(2)
!               deta(3)      var(3)
!                              :
!                              :
!                              :
!          deta(kdim-1)     var(kdim-1)
!                        ---------------  etadot(kdim-1)
!            deta(kdim)     var(kdim)
!                        ---------------
!
!   where  deta   = the eta (or sigma) thickness of model layers
!          etadot = time tendency of eta (or sigma); i.e., d(eta)/dt
!          var    = variable being advected
!          var_dt = the time tendency, d(var)/dt, for vertical advection
!          mask   = masking array of var, used for the eta coordinate,
!                   if var is underground then mask = 0.0, if var is
!                   above ground then mask = 1.0. The mask array
!                   is only needed for the fourth order scheme.
!
!   note: var, var_dt, and mask are located at model levels and all
!         have the same dimensions.
!
!=======================================================================

!#######################################################################
!#######################################################################

   subroutine vert_advec_2nd (dt, deta, var, etadot, var_dt, w)

!-----------------------------------------------------------------------
!
!     Routine for calculating vertical advection of any variable
!     located at full model levels using Euler-backward scheme.
!
!-----------------------------------------------------------------------
!
!   INPUT:   dt     = time step in seconds
!            deta   = eta thickness of model layers
!            var    = Input data at model levels
!            etadot = time tendency of eta (or sigma); i.e., eta-dot,
!                     the vertical velocity on the model's grid;
!                     located at model layer interfaces
!
!   OUTPUT:  var_dt = time tendency for the vertical advection of data
!                     in "var-units" per second.
!
!   OPTIONAL:     w = must be between 0.0 and 1.0, if w=1.0 (default)
!                     the Euler-backward scheme, if w=0.0 then
!                     Euler-forward scheme, if 0.< w < 1. then
!                     modified Euler-backward scheme.
!
!   NOTE:    The horizonal (1st & 2nd dimensions) of all 3d arrays
!            must be equal.  The vertical (or 3rd dimension) of
!            var & var_dt must be equal, while etadot must have
!            one fewer vertical levels.
!
!-----------------------------------------------------------------------

   real, intent(in)  :: dt, deta(:), var(:,:,:), etadot(:,:,:)
   real, intent(out) :: var_dt(:,:,:)
   real, intent(in), optional :: w

!-----------------------------------------------------------------------

   real, dimension(size(var,1),size(var,2),0:size(var,3)) :: a
   real, dimension(size(var,1),size(var,2),size(var,3)) :: vars
   real, dimension(size(deta)) :: rdeta

   integer  k, kdim, unit
   real     dtinv, wt, wdt

      wt = 1.0; if (present(w)) wt = w

      dtinv = 1./dt
      wdt   = - wt * dt * 0.5
      kdim  = size(var,3)
      rdeta = 1./deta

!---- define advective components -----

         a = 0.
      do k = 1, kdim-1
         a(:,:,k) = etadot(:,:,k) * (var(:,:,k+1) - var(:,:,k))
      enddo

!---- first pass (only apply if backward scheme) ----

      if ( wt > 0.0 ) then
         do k = 1, kdim
            vars(:,:,k) = var(:,:,k) + wdt * (a(:,:,k) + a(:,:,k-1))  &
                                     * rdeta(k)
         enddo
      else
         vars = var
      endif

!---- re-define advective components -----

      do k = 1, kdim-1
         a(:,:,k) = etadot(:,:,k) * (vars(:,:,k+1) - vars(:,:,k))
      enddo

!---- second pass ----

      do k = 1, kdim
         var_dt(:,:,k) = -0.5 * (a(:,:,k) + a(:,:,k-1)) * rdeta(k)
      enddo

!-----------------------------------------------------------------------

! 
!-----------------------------------------------------------------------
 
   end subroutine vert_advec_2nd


!#######################################################################


   subroutine vert_advec_4th (dt, deta, var, etadot, var_dt, w, mask)

!-----------------------------------------------------------------------
!
!   Routine for calculating 4th order vertical advection of any
!   variable located at full model levels using Euler-backward scheme.
!
!-----------------------------------------------------------------------

   real, intent(in)  :: dt, deta(:), var(:,:,:), etadot(:,:,:)
   real, intent(out) :: var_dt(:,:,:)
   real, intent(in), optional :: w, mask(:,:,:)

!-----------------------------------------------------------------------

   real, dimension(size(var,1),size(var,2),-1:size(var,3)+1) :: a
   real, dimension(size(var,1),size(var,2),size(var,3)) :: vars
   real, dimension(size(deta)) :: rdeta

   integer  k, kdim, unit
   real     dtinv, wt, wdt, f2, f4

      wt = 1.0; if (present(w)) wt = w

      dtinv = 1./dt
      wdt   = - wt * dt * 0.5
      kdim  = size(var,3)
      rdeta = 1./deta
      f2 = 7./6.; f4 = -1./6.

!---- define advective components -----

         a = 0.
      do k = 1, kdim-1
         a(:,:,k) = etadot(:,:,k) * (var(:,:,k+1) - var(:,:,k))
      enddo

!---- first pass (only apply if backward scheme) ----

      if ( wt > 0.0 ) then
         do k = 1, kdim
            vars(:,:,k) = var(:,:,k) + wdt *  &
                             (f2 * (a(:,:,k)   + a(:,:,k-1)) +  &
                              f4 * (a(:,:,k+1) + a(:,:,k-2))) * rdeta(k)
         enddo
         if (present(mask)) then
             where (mask < 0.1) vars = var
         endif
      else
         vars = var
      endif

!---- re-define advective components -----

      do k = 1, kdim-1
         a(:,:,k) = etadot(:,:,k) * (vars(:,:,k+1) - vars(:,:,k))
      enddo

!---- second pass ----

      do k = 1, kdim
         var_dt(:,:,k) = -0.5 * (f2 * (a(:,:,k)   + a(:,:,k-1)) +  &
                                 f4 * (a(:,:,k+1) + a(:,:,k-2)))   &
                               * rdeta(k)
      enddo

      if (present(mask)) var_dt = var_dt*mask

!-----------------------------------------------------------------------

   end subroutine vert_advec_4th

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
     
subroutine mpaceb_surface_flux_loop (                                  &
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
  real, dimension(size(t_atm(:))) :: t_surf_tmp, q_atm_tmp, flux_lw_tmp
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
   t_surf_tmp   = t_surf
   q_atm_tmp = q_atm
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
                 t_atm,     q_atm_tmp,  u_atm,       v_atm,            &
                 p_atm,     z_atm,      p_surf,      t_surf_tmp,       &
                 t_ca,      q_surf,     u_surf,      v_surf,           &
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
                       t_surf_tmp(k) = max( 200., min( 350.,             &
                            t_surf_tmp(k) - ( (flux_t(k)-SENFLUX) /      &
                            max(dhdt_surf(k),1.e-10) ) ) )
                            
                  end if                            
             enddo        
        
             ninner = ninner + 1
                     
        enddo  !------------- end of inner loop
        
        
        !---------------------------
        ! calculate surface flux

        call surface_flux (                                            &
                 t_atm,     q_atm_tmp,  u_atm,       v_atm,            &
                 p_atm,     z_atm,      p_surf,      t_surf_tmp,       &
                 t_ca,      q_surf,     u_surf,      v_surf,           &
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
                 q_atm_tmp(k) = max ( 1.e-08, min(0.05, q_atm_tmp(k)-&
                                 ((flux_q(k)-EVAPFLUX) /dedq_atm(k))))

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
                                            flux_v(k)**2. ) / rhotmp(k))
                                           
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

! ---> h1g
!                  if (mpp_pe() == mpp_root_pe()) then
!                  print *, 'desired fluxes  (sensible, latent): ',     &
!                       SENFLUX, hlv*EVAPFLUX
!                  print *, 'iterated fluxes (sensible, latent): ',     &
!                       flux_t, hlv*flux_q
!                  end if
! <--- h1g

   ! return original longwave flux
   ! 
   ! Note this is needed so that the upward longwave flux is not affected
   ! by the changed t_surf during iteration.

   flux_lw = flux_lw_tmp

   ! reset q_surf, surface specific humidity
   q_surf = q_atm + flux_q / (rhotmp*cd_q*w_atm)  
                     
!
!-----------------------------------------------------------------------
           
end subroutine mpaceb_surface_flux_loop

!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_mpaceb_flx( flux_t, flux_q )

  implicit none
  real, intent(out) :: flux_t, flux_q

  flux_t = SENFLUX
  flux_q = EVAPFLUX

end subroutine get_mpaceb_flx

!########################################################################
! This subroutine returns imposed SST

subroutine get_mpaceb_sst( sst )

  implicit none
  real, intent(out) :: sst

  sst = tskin

end subroutine get_mpaceb_sst

!########################################################################
!########################################################################


!########################################################################
 subroutine   ice_nucl_mpaceb( rhi_in, Ni,   rh_crit_1d )
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
end subroutine  ice_nucl_mpaceb
!########################################################################
!########################################################################


end module scm_mpaceb_mod
