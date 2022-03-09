module scm_astex_mod
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       ASTEX Lagrangians I and II FORCING MODULE
!
!       February 2002
!       Contact person: Steve Klein
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This module provides the large scale forcing necessary
!       to run the single column model to simulate either ASTEX
!       Lagrangian I or II.  
!
!       The two ASTEX Lagrangian experiments were comprehensive studies
!       of the evolution of a boundary layer airmass over roughly 36 
!       hour periods.  In each Lagrangian, aircraft observations were 
!       taken an continuously as possible following the mean motion of 
!       the air in the marine boundary layer (MBL).
!
!       Periods of the Lagrangians:
!
!       Lagrangian 1: 16 UTC 12 June-10 UTC 14 June, 1992 
!       Lagrangian 2: 22 UTC 18 June-14 UTC 20 June, 1992
!
!       References
!       ----------
!
!       Bretherton, C. S., S. K. Krueger, M. C. Wyant, P. Bechtold,
!       E. van Meijgaard, B. Stevens, and J. Teixeira, 1999: A GCSS
!       boundary-layer cloud model intercomparison study of the first 
!       ASTEX lagrangian experiment. Boundary Layer Meteorology, 
!       vol. 93, pp. 341-380.
!
!       Bretherton, C. S. and R. Pincus 1995; Cloudiness and Marine
!       Boundary Layer Dynamics in the ASTEX Lagrangian Experiments.
!       Part I:  Synoptic setting and vertical structure.
!       J. Atmos. Sci., vol. 52, pp. 2707-2723.
!
!       Bretherton, C. S., P. A. Austin and S. T. Siems, 1995:
!       Cloudiness and Marine Boundary Layer Dynamics in the ASTEX 
!       Lagrangian Experiments.  Part II:  Cloudiness, drizzle, 
!       surface fluxes and entrainment.  J. Atmos. Sci., vol. 52,
!       2724-2735.
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
use              mpp_mod, only: input_nml_file
#else
use              fms_mod, only: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod
   use      constants_mod, only:  rdgas, cp_air, tfreeze, hlv, hls, rvgas, &
                                  grav, stefan, pi

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
                               ak, bk, rlon, rlat, ng_d, nt_prog, get_eta_level
                                  
   implicit none
   private

   public astex_data_read, astex_forc_init, astex_forc_end, update_astex_forc, &
          astex_forc_diagnostic_init, get_astex_sst

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following subroutines:
!
!            astex_data_read   reads forcing data files and initializes
!                              any needed parameters of the run
!            astex_forc_init   initializes prognostic variables:           
!                              T,u,v,qv,ql,qi
!            astex_forc_end    deallocates allocated spaces
!            update_astex_forc a call to this routine returns the 
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
!       lwp_forc       liquid water path (kg/m3)
!       pib_forc       inversion base pressure (Pa)
!       cdn_forc       cloud drop # (1/m3)
!

character(len=8) :: mod_name = 'forcing'


real, allocatable, dimension(:,:) :: plev_stand,t_stand,q_stand,O3_stand
real, allocatable, dimension(:)   :: plev_forc,ps_forc,sst_forc,div_forc
real, allocatable, dimension(:)   :: lon_forc,lat_forc,lwdn_forc
real, allocatable, dimension(:)   :: cdn_forc,lwp_forc,pib_forc
real, allocatable, dimension(:,:) :: t_forc,q_forc,ql_forc,u_forc,v_forc
real, allocatable, dimension(:,:) :: ug_forc, vg_forc, N_forc

type(TIME_TYPE), allocatable, dimension(:)  :: time_layer

!----------Diagnostic data----------------------------------------------

integer :: id_temp_obs, id_qv_obs, id_ql_obs, id_omega_obs,id_uwnd_obs,&
           id_vwnd_obs, id_lwp_obs, id_pib_obs, id_pib_mod, id_nudg_t, &
           id_nudg_q
                         
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       ivert_stand              # of vertical levels to standard atmos.
!       ivert                    # of vertical levels of forcing data
!       itime                    # of times of astex data being used
!       itime_lagr               # of times in each astex lagrangian
!
!       p00                      reference pressure (pascals)
!
!       tskin                    sea surface temperature (K)
!
!       LAGRNUM                  number indicating which lagrangian is 
!                                being run (1 or 2)
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
!       thermo_relaxation_tau    relaxation time scale for temperature
!                                and moisture profiles above the inver-
!                                sion (seconds)
!
!       p_omega_zero             pressure above which omega is set to
!                                to zero
!
!       p_cld_zero               pressure above which clouds are forced
!                                to be zero
!

integer, public                :: LAGRNUM = 1
integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, private               :: itime
integer, private               :: ivert = 41
integer, private               :: ivert_stand = 60
integer, dimension(2), private :: itime_lagr
real,    private               :: missing_value = -999.
real,    private               :: p00 = 100000.
real,    public                :: tskin = 300.
real,    public                :: p_omega_zero = 70000.
real,    public                :: p_cld_zero = 50000.
real,    public                :: wind_relaxation_tau = 3.*3600.
real,    public                :: thermo_relaxation_tau = 3.*3600.
logical, public                :: vert_advec_cond = .true.
logical, private               :: astex_forc_initialized = .false.
logical, public                :: do_netcdf_restart = .true.

data itime_lagr     / 43, 41 /

!--------------------- version number ----------------------------------
!
        
character(len=128) :: Version = '$Id$'
character(len=128) :: Tagname = '$Name$'
        
namelist /SCM_ASTEX_NML/ do_netcdf_restart,                           &
                        LAGRNUM, tracer_vert_advec_scheme,           &
                        temp_vert_advec_scheme, vert_advec_cond,     &
                        wind_relaxation_tau,                         &
                        thermo_relaxation_tau, p_omega_zero, p_cld_zero

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains


!#######################################################################
!#######################################################################

subroutine astex_data_read()

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
!       io,ierr,tmpint       dummy integer variable
!       adum,ahour1,ahour2   dummy character variables
!       fname,prename        file name variables
!       junday               day in June 1992
!       hour                 hours
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Internal variables
!  ------------------

integer                                   :: i,j,k,l,t,unit,io,ierr
integer                                   :: tmpint
character*23                              :: tracer_ascheme,temp_ascheme
character*20                              :: adum, fname
character*13                              :: prename
character*1                               :: ahour1
character*2                               :: ahour2
real                                      :: junday, hour
real                                      :: tmpr, tmprl
     
       ! if the constructor is already called, then return, otherwise 
       ! set astex_forc_initialized .TRUE.

       if (astex_forc_initialized) return
       astex_forc_initialized = .true.

!-----------------------------------------------------------------------
!      ----- read namelist -----
#ifdef INTERNAL_FILE_NML
       read (input_nml_file, nml=SCM_ASTEX_NML, iostat=io)
       ierr = check_nml_error(io, 'SCM_ASTEX_NML')
#else
       if (file_exist('input.nml')) then
            unit = open_namelist_file()
            io=1
            do while (io .ne. 0)
                 read  (unit, nml=SCM_ASTEX_NML, iostat=io, end=10)
            enddo
  10        call close_file (unit)
       endif
#endif 

       !dummy check -------         
       if (temp_vert_advec_scheme .gt. 6 ) then
            write (adum,'(i1)') temp_vert_advec_scheme
            call error_mesg ( 'astex_data_read in scm_astex_mod',     &
                 'Bad temp_vert_advec_scheme : '//adum, FATAL )
       end if

       if (tracer_vert_advec_scheme .gt. 6 ) then
            write (adum,'(i1)') tracer_vert_advec_scheme
            call error_mesg ( 'astex_data_read in scm_astex_mod',     &
                 'Bad tracer_vert_advec_scheme : '//adum, FATAL )
       end if

       select case(tracer_vert_advec_scheme)
            case(1)
                 write (tracer_ascheme,'(a)') 'SECOND_CENTERED        '
            case(2)
                 write (tracer_ascheme,'(a)') 'FOURTH_CENTERED        '
            case(3)
                 write (tracer_ascheme,'(a)') 'FINITE_VOLUME_LINEAR   '
            case(4)
                 write (tracer_ascheme,'(a)') 'FINITE_VOLUME_PARABOLIC'
            case(5)
                 write (tracer_ascheme,'(a)') 'SECOND_CENTERED_WTS    '
            case(6)
                 write (tracer_ascheme,'(a)') 'FOURTH_CENTERED_WTS    '
       end select
       
       select case(temp_vert_advec_scheme)
            case(1)
                 write (temp_ascheme,'(a)') 'SECOND_CENTERED        '
            case(2)
                 write (temp_ascheme,'(a)') 'FOURTH_CENTERED        '
            case(3)
                 write (temp_ascheme,'(a)') 'FINITE_VOLUME_LINEAR   '
            case(4)
                 write (temp_ascheme,'(a)') 'FINITE_VOLUME_PARABOLIC'
            case(5)
                 write (temp_ascheme,'(a)') 'SECOND_CENTERED_WTS    '
            case(6)
                 write (temp_ascheme,'(a)') 'FOURTH_CENTERED_WTS    '
       end select
       
       unit = stdlog()
       if ( mpp_pe() == 0 ) then
            write (unit,'(/,80("="),/(a))') trim(version), trim(Tagname)
            write (unit,'(a24,i1)') ' ASTEX lagrangian #  = ',LAGRNUM
            write (unit,'(a28,a23)') ' tracer_vert_advec_scheme = ', &
                                       tracer_ascheme
            write (unit,'(a26,a23)') ' temp_vert_advec_scheme = ', &
                                       temp_ascheme
            write (unit,*) ' vert_advec_cond  = ', vert_advec_cond
            write (unit,*) ' wind relaxation time = ',                 &
                             wind_relaxation_tau
            write (unit,*) ' thermo relaxation time = ',               &
                             thermo_relaxation_tau
            write (unit,*) ' p_omega_zero = ', p_omega_zero
            write (unit,*) ' p_cld_zero = ',p_cld_zero
       endif
      
       call close_file (unit)
      
       if (LAGRNUM .gt. size(itime_lagr,1)) then
            write (adum,'(i1)') LAGRNUM
            call error_mesg ( 'astex_data_read in scm_astex_mod',    &
                 'Bad case number : '//adum, FATAL )
       end if
        
       !set up number of times in the lagrangian   
       itime = itime_lagr(LAGRNUM)
      
!-----------------------------------------------------------------------
!
!      ALLOCATE STORAGE
!

       if (allocated(plev_stand)) deallocate (plev_stand)
           allocate(plev_stand(ivert_stand,2))
       if (allocated(t_stand)) deallocate (t_stand)
           allocate(t_stand(ivert_stand,2))
       if (allocated(q_stand)) deallocate (q_stand)
           allocate(q_stand(ivert_stand,2))
       if (allocated(O3_stand)) deallocate (O3_stand)
           allocate(O3_stand(ivert_stand,2))
       
       if (allocated(plev_forc)) deallocate (plev_forc)
           allocate(plev_forc(ivert))
       if (allocated(time_layer)) deallocate (time_layer)
           allocate(time_layer(itime))
       
       if (allocated(t_forc)) deallocate (t_forc)
           allocate(t_forc(ivert,itime))
       if (allocated(q_forc)) deallocate (q_forc)
           allocate(q_forc(ivert,itime))
       if (allocated(ql_forc)) deallocate (ql_forc)
           allocate(ql_forc(ivert,itime))
       if (allocated(u_forc)) deallocate (u_forc)
           allocate(u_forc(ivert,itime))
       if (allocated(v_forc)) deallocate (v_forc)
           allocate(v_forc(ivert,itime))
       if (allocated(ug_forc)) deallocate (ug_forc)
           allocate(ug_forc(ivert,itime))
       if (allocated(vg_forc)) deallocate (vg_forc)
           allocate(vg_forc(ivert,itime))
       if (allocated(N_forc)) deallocate (N_forc)
           allocate(N_forc(ivert,itime))
       
       if (allocated(ps_forc)) deallocate (ps_forc)
           allocate(ps_forc(itime))
       if (allocated(sst_forc)) deallocate (sst_forc)
           allocate(sst_forc(itime))
       if (allocated(div_forc)) deallocate (div_forc)
           allocate(div_forc(itime))
       if (allocated(lat_forc)) deallocate (lat_forc)
           allocate(lat_forc(itime))
       if (allocated(lon_forc)) deallocate (lon_forc)
           allocate(lon_forc(itime))
       if (allocated(lwdn_forc)) deallocate (lwdn_forc)
           allocate(lwdn_forc(itime))
       if (allocated(cdn_forc)) deallocate (cdn_forc)
           allocate(cdn_forc(itime))
       if (allocated(lwp_forc)) deallocate (lwp_forc)
           allocate(lwp_forc(itime))
       if (allocated(pib_forc)) deallocate (pib_forc)
           allocate(pib_forc(itime))
       
!-----------------------------------------------------------------------
! 
!     READ IN MIDLATITUDE ATMOSHPHERE DATA
!

       fname='INPUT/mls.dat'
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
       call close_file(unit)
       
       !----record which file was read
       unit= stdlog()
       if (mpp_pe() == 0 ) write (unit,'(a)') 'used data file: '// fname 
       call close_file(unit)
       
       !----------------------------
       ! read winter

       fname='INPUT/mlw.dat'
       call mpp_open(unit,fname,action=MPP_RDONLY)
      
       !skip to beginning of data
       do j=1,4
            read(unit,'(a)') adum
       enddo
            
       !read data
       do j=1,ivert_stand
            read (unit,'(0p,2f10.3,1p,2e12.4)') &
                 plev_stand(ivert_stand+1-j,1), t_stand(ivert_stand+1-j,1),&
                 q_stand(ivert_stand+1-j,1),O3_stand(ivert_stand+1-j,1)  
       end do
       call close_file(unit)
       
       !----record which file was read
       unit= stdlog()
       if (mpp_pe() == 0 ) write (unit,'(a)') 'used data file: '// fname 
       call close_file(unit)
      
       !change units
       
       !(a) mb to Pascals
       plev_stand(:,:)=100.*plev_stand(:,:)
             
       !(b) mass mixing ratio to specific humidity
       q_stand(:,:)=q_stand(:,:)/(1.+q_stand(:,:))
       O3_stand(:,:)=O3_stand(:,:)/(1.+O3_stand(:,:))
           
!-----------------------------------------------------------------------
! 
!      READ IN ASTEX FORCING DATA
!
!
!       

       ! ---> read layer files
       select case(LAGRNUM)
            case(1)
                 prename='INPUT/lagr1_h'
            case(2)
                 prename='INPUT/lagr2_h'
       end select

       do 1000 t = 1, itime

            !create file name
            if (t .le. 10) then
                 write(ahour1,'(i1)') t-1
                 fname = prename//'0'//ahour1
            else
                 write(ahour2,'(i2)') t-1
                 fname = prename//ahour2
            end if
            !write(6,*) fname
            call mpp_open(unit,fname,action=MPP_RDONLY)
             
            read(unit,*)
            read(unit,*)
            read(unit,*)
            read(unit,12) junday, lat_forc(t), lon_forc(t),            &
                 sst_forc(t), ps_forc(t), div_forc(t), lwdn_forc(t),   &
                 tmpint
            cdn_forc(t)=real(tmpint)
12          format(1X,f8.5,1X,f7.3,1X,f7.3,1X,f7.3,1X,f6.1,1X,f4.1,2X, &
                   f5.1,2X,i3)
            !write(6,12) junday, lat_forc(t), lon_forc(t), sst_forc(t), &
            !     ps_forc(t), div_forc(t), lwdn_forc(t), tmpint
            read(unit,*)
            read(unit,*)
            lwp_forc(t) = 0.0
            do k = 1, ivert
                 l=ivert+1-k
                 read(unit,13) tmpint,t_forc(l,t),q_forc(l,t),         &
                               u_forc(l,t),v_forc(l,t),ug_forc(l,t),   &
                               vg_forc(l,t),ql_forc(l,t),N_forc(l,t)     
13               format(1X,i4,3(f9.2),4(f8.2),f8.0)
                 plev_forc(l) = real(tmpint)
                 lwp_forc(t)  = lwp_forc(t) + ql_forc(l,t)*1.e-03*      &
                      10.*100./grav
            end do
            call close_file(unit)  

            !----record which file was read
            call mpp_open(unit,'scm.out',action=MPP_APPEND)
            if (mpp_pe() == 0 ) write (unit,'(a)')                  &
                 'used data file: '//fname
            call close_file (unit)

            !set time variable for VAR data 
            hour = 24. * (junday - real(int(junday)))
            time_layer(t) = set_date(1992,6,int(junday),int(hour),0,0)
             
1000   continue

       ! ---> read scalar file
       select case(LAGRNUM)
            case(1)
                 fname='INPUT/lagr1_scalars'
            case(2)
                 fname='INPUT/lagr2_scalars'
       end select

       !write(6,*) fname
       call mpp_open(unit,fname,action=MPP_RDONLY)
       read(unit,*)
       read(unit,*)
       read(unit,*)
       read(unit,*)
       read(unit,*)
       do t = 1, itime
            read(unit,14)  tmpint
            pib_forc(t)=real(tmpint)*100.
14          format(66X,i3)
            !write(6,'(f12.1)') pib_forc(t)
       enddo
       call close_file(unit)  
            
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
       
       ! (c) impose a minimum boundary on specific humidity
       do t = 1, itime
       do k = 1, ivert
             q_forc(k,t)  = max(4.e-06,q_forc(k,t))
       enddo
       enddo
              
       ! (d) change pressure variables from mb to Pascals
       plev_forc(:) =  100.*plev_forc(:)
       ps_forc(:)   =  100.*  ps_forc(:)

       ! (e) convert theta to Temperature
       do k = 1, ivert
            t_forc(k,:) = t_forc(k,:)*((plev_forc(k)/100000.)**(rdgas/cp_air))
       enddo

       ! (f) convert 10-6 sec-1 to sec-1
       div_forc(:) = div_forc(:)*1.E-06

       ! (g) convert latitude and longitude to radians
       lat_forc(:) = pi*lat_forc(:)/180.
       lon_forc(:) = pi*(360.+lon_forc(:))/180.

       ! (h) convert cloud drop number to per cubic meter
       cdn_forc(:) = cdn_forc(:)*1.e+06
       N_forc(:,:) = N_forc(:,:)*1.e+06
       

       !--------------------------------------
       ! reading tskin from the restart file.

       if( file_exist('INPUT/scm_forc.res.nc') ) then

           if(mpp_pe() == mpp_root_pe() )                            &
           call mpp_error ('scm_forc_mod',                           &
                           'Reading netCDF formatted restart file.', &
                            NOTE )
           call read_data('INPUT/scm_forc.res.nc', 'Tskin', tskin )

       else if (file_exist('INPUT/scm_forc.res')) then

          unit=open_restart_file('INPUT/scm_forc.res',action='read')
          read(unit)  tskin
          call close_file(unit)

       endif       
        

!-----------------------------------------------------------------------
        
end subroutine astex_data_read


!#######################################################################
!#######################################################################


subroutine astex_forc_init(time_interp,pdamp,As,Bs)
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
!      pdamp           maximum pressure for which profile is set to a 
!                      reference midlatitude profile 
!      As, Bs          A's and B's of half levels in hybrid coordinate
!                      phalf(k) = A(k) + B(k) * psurf
!
!      -------------
!      INPUT/OUTPUT:
!      -------------
!
!      pt            temperature in degrees Kelvin 
!      u            zonal wind (m/s)
!      v            meridional wind (m/s)
!      qv           water vapor specific humidity (kg water vapor/ kg air)
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
!      qsat            saturation specific humidity (kg vapor/kg air)
!      gamma           temporary variable        
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

type(TIME_TYPE)                          :: time_interp
real,  intent (IN)                       :: pdamp
real,  intent (IN)   , dimension(:)      :: As,Bs


!  Internal variables
!  ------------------

integer  :: itime_less, itime_more, KDIM, k, jstart, jstart_stand
integer  :: nlev_interp, days, months, years, seconds, minutes, hours
real     :: weight_less, weight_more, summerweight
real     :: tmplwp, tmplwpforc

real, dimension(ivert+2+ivert_stand)             :: plev_int
real, dimension(ivert+2+ivert_stand)             :: field_int_less, &
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
       
       ! --- find out # of vertical levels
       KDIM = size(pt,3)

! --- CREATE PFULL AND PHALF  ---------------------------------------- !

       ! --- find out time in sfc pressure field
       call interp_time(time_interp,time_layer,itime_less,itime_more,  &
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

! --- compute starting level for midlatitude reference profile --- !

      jstart = 1
      jstart_stand = 1
      do while (plev_forc(ivert).lt.plev_stand(jstart_stand,1))
         jstart_stand = jstart_stand + 1
      enddo

! --- create average mid-latitude profile based upon month ----------- !       

      call get_date(time_interp,years,months,days,hours,minutes,seconds)
      summerweight = 1. - ((real(abs(months - 7)))/6.)
      
! --- CREATE INITIAL T  ---------------------------------------------- !

       ! --- find nearest time to time_interp in LAYERED fields
       call interp_time(time_interp,time_layer,itime_less,itime_more,  &
                        weight_less,weight_more)
       
       ! --- create field data for interpolator
       ! --- NOTES:
       ! ---
       ! --- level 1 corresponds to the surface and it is assumed that &
       ! ---     theta(surface) = theta(BOT of Layered fields)
       ! ---
       ! --- pdamp is ignored
            
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
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:nlev_interp),field_int_less(1:nlev_interp),     &
            field_int_more(1:nlev_interp),tmp_ans)

       !---- assign pt
       do k = 1,KDIM
              pt(:,:,k) = tmp_ans(k)
       enddo
      
! --- CREATE INITIAL Qv  --------------------------------------------- !

 
       ! --- find nearest time to time_interp in Layered fields ------ !
       ! --- no need since this was done for T 
      
       !--- create field data for interpolator
       !--- note: *level 1 corresponds to the surface and it is assumed
       !    that   qv(surface) = qv(BOT)
       !           qv(p<ptop)  = qv(TOP)
       !      
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
       do k = 1,KDIM
            q(:,:,k,nsphum)=tmp_ans(k)
       enddo
  
! --- CREATE INITIAL QL  --------------------------------------------- !

 
       ! --- find nearest time to time_interp in Layered fields ------ !
       ! --- no need since this was done for T 
       
       !--- create field data for interpolator
       !--- note: *level 1 corresponds to the surface and it is assumed 
       !    that   ql(surface)  =  ql(BOT)
       !          *ql(p<ptop)   =  ql(TOP)
       !      
       
       plev_int(1) = ps_init(1,1)
       field_int_less(1) = ql_forc(1,itime_less)
       field_int_more(1) = ql_forc(1,itime_more)
       plev_int(2:(ivert+1)) = plev_forc(1:ivert)
       field_int_less(2:(ivert+1)) = ql_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1)) = ql_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = 0.
       field_int_more(ivert+2) = 0.
       tmp_p(:) = pfull(1,1,:)
        
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
               plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),      &
               field_int_more(1:(ivert+2)),tmp_ans)
       
       !---- assign answer to ql
       do k = 1,KDIM
              q(:,:,k,nql) = tmp_ans(k)
       enddo

! --- CREATE INITIAL QI,QA  ------------------------------------------ !

       q(:,:,:,nqi)  =  0.
       where (q(:,:,:,nql) .gt. 0.) 
              q(:,:,:,nqa)  =  1.
       end where
            
      
! --- CREATE INITIAL U,V  -------------------------------------------- !
      
       ! --- find nearest time to time_interp in LAYERED FIELDS
       ! --- no need since this was done for T 
      
       ! --- create field data for interpolator
       ! --- NOTE: 
       ! --- 
       ! --- level 1 corresponds to the surface and it is assumed that
       !           u,v(surface)  =  u,v(BOT)
       !           u,v(p<ptop)   =  u,v(TOP)
             
       plev_int(1) = ps_init(1,1)
       field_int_less(1) = u_forc(1,itime_less)
       field_int_more(1) = u_forc(1,itime_more)
       plev_int(2:(ivert+1)) = plev_forc(1:ivert)
       field_int_less(2:(ivert+1)) = u_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1)) = u_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = u_forc(ivert,itime_less)
       field_int_more(ivert+2) = u_forc(ivert,itime_more)
       tmp_p(:) = pfull(1,1,:) 

       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
      
       ! --- assign answer to u
       do k  =  1, KDIM
            ua(:,:,k) = tmp_ans(k)
       enddo
       u_srf(:,:)=ua(:,:,KDIM)

       ! --- create field data for interpolator
                
       plev_int(1) = ps_init(1,1)
       field_int_less(1) = v_forc(1,itime_less)
       field_int_more(1) = v_forc(1,itime_more)
       plev_int(2:(ivert+1)) = plev_forc(1:ivert)
       field_int_less(2:(ivert+1)) = v_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1)) = v_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = v_forc(ivert,itime_less)
       field_int_more(ivert+2) = v_forc(ivert,itime_more)
       tmp_p(:) = pfull(1,1,:)
      
       ! --- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
      
       ! --- assign answer to v
       do k  =  1, KDIM
            va(:,:,k) = tmp_ans(k)
       enddo
       ! Initialize surface winds
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

! ---- interpolate lwp_forc and scale ql ----------------------------- !
     
       !---- interpolate lwp field
       tmp_p(:) = 50000.
       plev_int(:) = p00
       field_int_less(:) = lwp_forc(itime_less)
       field_int_more(:) = lwp_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,&
                     field_int_less,field_int_more,tmp_ans)
       tmplwpforc = tmp_ans(1)

       tmplwp = 0.
       do k = 1, kdim
               tmplwp = tmplwp + q(1,1,k,nql)*&
                                 (phalf(1,1,k+1)-phalf(1,1,k))/grav
       enddo
       
       if (tmplwp .gt. 0. .and. tmplwpforc .gt. 0.) then
            q(:,:,:,nql) = q(:,:,:,nql) * tmplwpforc / tmplwp
       end if       
       
! 
!-----------------------------------------------------------------------
 
 end subroutine astex_forc_init

!#######################################################################
!#######################################################################


subroutine astex_forc_end ()

  integer :: unit

  if (.not.astex_forc_initialized) return
  astex_forc_initialized = .false.

  deallocate (  plev_stand, t_stand, q_stand, O3_stand )
  deallocate (  plev_forc , t_forc, q_forc, ql_forc, u_forc, v_forc )
  deallocate (  ug_forc, vg_forc, N_forc )
  deallocate (  lat_forc, lon_forc, ps_forc, sst_forc, div_forc )
  deallocate (  lwdn_forc, cdn_forc, time_layer )
  
  !--------------------------------------
  ! write tskin to the restart file.

  if( do_netcdf_restart ) then
 
     if (mpp_pe() == mpp_root_pe()) &
     call write_data( 'RESTART/scm_forc.res.nc', 'Tskin', tskin )
 
  else
        
     unit=open_restart_file('RESTART/scm_forc.res',action='write')
     if(mpp_pe() == mpp_root_pe())  write(unit) tskin
     call close_file(unit)       

  endif        
! 
!-----------------------------------------------------------------------
 

end subroutine astex_forc_end


!#######################################################################
!#######################################################################

subroutine astex_forc_diagnostic_init(axes, Time)

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
       
     
       id_nudg_t = register_diag_field (mod_name, 'nudg_t',            &
          axes(1:3), Time,                                             &
          'Temperature tendency from Nudging', 'K/sec',                &
          missing_value = missing_value )
          
       id_nudg_q = register_diag_field (mod_name, 'nudg_q',            &
          axes(1:3), Time,                                             &
          'Water vapor tendency from Nudging', 'kg/kg/sec',            &
          missing_value = missing_value )
 
                                                                   
! 
!-----------------------------------------------------------------------
 
end subroutine astex_forc_diagnostic_init


!#######################################################################
!#######################################################################


subroutine update_astex_forc(time_interp,time_diag,dt_int,pdamp )
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
!                   be set to zero.  <NOTE THIS IS IGNORED>
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
!      pt           temperature in degrees Kelvin 
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
!      KDIM             # of vertical levels to pt array
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
!      nudg_t           (T_obs - T_model)/thermo_relaxation_tau (K/sec)
!      nudg_q           (q_obs - q_model)/thermo_relaxation_tau 
!                       (kg vapor/kg air/sec)
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

type(TIME_TYPE), intent(in)              :: time_interp,time_diag,dt_int
real,  intent (IN)                       :: pdamp

!  Internal variables
!  ------------------

integer                                          :: days,months,years
integer                                          :: seconds,minutes
integer                                          :: hours
integer                                          :: i,j,k,itime_less
integer                                          :: itime_more,KDIM
integer                                          :: dt_seconds,dt_days
logical                                          :: used
real                                             :: tmplwp, tmplwpforc
real                                             :: weight_less
real                                             :: weight_more
real, dimension(ivert+2)                         :: plev_int
real, dimension(ivert+2)                         :: field_int_less
real, dimension(ivert+2)                         :: field_int_more
real, dimension(size(pt,3))                         :: tmp_p,tmp_ans
real, dimension(size(pt,1),size(pt,2))              :: ps_init, tmp_res2
real, dimension(size(pt,1),size(pt,2))              :: pinv_mod,pinv_obs
real, dimension(size(pt,1),size(pt,2))              :: tmpr
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pfull,t_obs,q_obs
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: ql_obs,u_obs,v_obs
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pi_fac, dp
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: adv_vert,tmp_res
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: nudg_t,nudg_q
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: phalf, omega_half

#include "fv_point.inc"

!CODE
!

!check to see if time is outside of range 
if (time_interp < time_layer(1)) then
    call get_date(time_interp,years,months,days,hours,minutes,seconds)
    print *, 'DATE REQUESTED: ',years,months,days,hours,minutes,seconds
    call error_mesg ( 'update_astex_forc in scm_astex_mod',            &
         'requested time is before beginning of variational forcing '//&
         'dataset', FATAL )
    
end if
if (time_interp > time_layer(itime)) then
    call get_date(time_interp,years,months,days,hours,minutes,seconds)
    print *, 'DATE REQUESTED: ',years,months,days,hours,minutes,seconds
    call error_mesg ( 'update_astex_forc in scm_astex_mod',            &
         'requested time is after the end of variational forcing '//   &
         'dataset', FATAL )
end if

! --- find out # of vertical levels
KDIM = size(pt,3)

! --- initialize nudg tendencies
nudg_t = 0.
nudg_q = 0.

! --- compute timestep interval in seconds and days
call get_time(dt_int,dt_seconds,dt_days)
dt_seconds = dt_seconds + 86400*dt_days

! --- UPDATE PFULL AND PHALF, COMPUTE OMEGA SURFACE, AND PI FACTOR --- !

       ! --- find out time in sfc pressure field ---- !
       call interp_time(time_interp,time_layer,itime_less,itime_more,  &
            weight_less,weight_more)
 
       !---- interpolate surface pressure field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=ps_forc(itime_less)
       field_int_more(:)=ps_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       ps_init(:,:) = tmp_ans(1)
       ps = ps_init

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
       field_int_more(:) = pib_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       pinv_obs(:,:) = tmp_ans(1)
       
       if (id_pib_obs > 0) then
       used = send_data( id_pib_obs, pinv_obs(:,:), time_diag, 1, 1)
       end if  

       !---- find model inversion level as level with temperature
       !---- minimum,  this is the same definition as used in Bretherton
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
            call error_mesg ('update_astex_forc in scm_astex_mod',     &
            'pinv_obs < p_omega_zero', FATAL )
            
! ----------------- COMPUTE OMEGA ON HALF AND FULL LEVELS ------------ !

       !---- interpolate surface divergence field
       tmp_p(:) = 50000.
       plev_int(:) = p00
       field_int_less(:) = div_forc(itime_less)
       field_int_more(:) = div_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)

       !---- compute omega_half
       omega_half = 0.0
       do k = 1, KDIM + 1
            where(phalf(:,:,k) .gt. pinv_obs(:,:))
                 omega_half(:,:,k) = tmp_ans(1)*(ps_init(:,:) -        &
                      phalf(:,:,k))
            endwhere
            where(phalf(:,:,k) .le. pinv_obs(:,:) .and. phalf(:,:,k)   &
                 .gt. p_omega_zero)
                 omega_half(:,:,k) = tmp_ans(1)*(ps_init(:,:) -        &
                      pinv_obs(:,:))*(phalf(:,:,k) - p_omega_zero) /   &
                      (pinv_obs-p_omega_zero)
            endwhere
       enddo
       
       !---- compute omega_full
       omga = 0.0
       do k = 1, KDIM   
            where(pfull(:,:,k) .gt. pinv_obs(:,:))
                 omga(:,:,k) = tmp_ans(1)*(ps_init(:,:) -        &
                      pfull(:,:,k))
            endwhere
            where(pfull(:,:,k) .le. pinv_obs(:,:) .and. pfull(:,:,k)   &
                 .gt. p_omega_zero)
                 omga(:,:,k) = tmp_ans(1)*(ps_init(:,:) -        &
                      pinv_obs(:,:))*(pfull(:,:,k)-p_omega_zero) /     &  
                      (pinv_obs-p_omega_zero)
            endwhere
       enddo        

! ----------------- PRODUCE TEMPERATURE FORCING ---------------------- !

       !---- initialize variables
       t_dt     = 0.
       nudg_t   = 0.
       adv_vert = 0.       
       
       !--------- COMPUTE ADIABATIC COMPRESSION  --------------------- !
      
       !---- compute adiabatic compression
       t_dt(:,:,:)=rdgas*pt(:,:,:)*omga(:,:,:)/cp_air/pfull(:,:,:)
      
       !---- COMPUTE VERTICAL T ADVECTION ---------------------------- !
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
          call vert_advection(real(dt_seconds),omega_half,delp,pt,&
               adv_vert, scheme = vadvec_scheme,          &
               form = ADVECTIVE_FORM)

       !---- interpolate t_forc
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=t_forc(1,itime_less)* &
                           ((plev_int(1)/plev_forc(1))**(rdgas/cp_air))
       field_int_more(1)=t_forc(1,itime_more)* &
                           ((plev_int(1)/plev_forc(1))**(rdgas/cp_air))
       plev_int(2:(ivert+1))=plev_forc(1:ivert)
       field_int_less(2:(ivert+1))=t_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1))=t_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = t_forc(ivert,itime_less)
       field_int_more(ivert+2) = t_forc(ivert,itime_more)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
       
       !---- compute damping term
       !     note that damping to observed profile occurs only at levels
       !     ~150 meters above the nomial inversion. also note assumed 
       !     density of 1 kg/m3
       
       do k=1,KDIM
            where(pfull(:,:,k).le.(pinv_obs(:,:)-grav*150.))
                 nudg_t(:,:,k)=(tmp_ans(k)-pt(:,:,k))/                  &
                      max(thermo_relaxation_tau,real(dt_seconds))
                 t_dt(:,:,k)     = 0.0
                 adv_vert(:,:,k) = 0.0
            endwhere
       enddo
       
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result
       t_obs(1,1,:) = tmp_ans
       if (id_temp_obs > 0) then
       used = send_data ( id_temp_obs, t_obs(:,:,:), time_diag, 1, 1, 1)
       end if  
       
       !---- COMPUTE FINAL ADVECTIVE T TENDENCY FROM THE SUM OF :  --- !
       !---- RELAXATION TENDENCY, VERTICAL ADVECTION AND ADIABATIC --- !
       !---- COMPRESSION                                           --- !
       t_dt(:,:,:) = t_dt(:,:,:) + adv_vert(:,:,:) + nudg_t(:,:,:)
      
! ----------------- PRODUCE MOISTURE FORCING ------------------------- !

       !---- initialize variables
       q_dt(:,:,:,nsphum) = 0.0
       nudg_q   = 0.
       adv_vert = 0.       
       
       !---- COMPUTE VERTICAL Q ADVECTION ---------------------------- !
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
          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nsphum),&
               adv_vert, scheme = vadvec_scheme,          &
               form = ADVECTIVE_FORM)

       !---- interpolate q_forc
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=q_forc(1,itime_less)
       field_int_more(1)=q_forc(1,itime_more)
       plev_int(2:(ivert+1))=plev_forc(1:ivert)
       field_int_less(2:(ivert+1))=q_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1))=q_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = q_forc(ivert,itime_less)
       field_int_more(ivert+2) = q_forc(ivert,itime_more)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
       
       !---- compute damping term
       do k=1,KDIM          
            where(pfull(:,:,k).le.(pinv_obs(:,:)-grav*150.))
                 nudg_q(:,:,k)=(tmp_ans(k)-q(:,:,k,nsphum))/                            &
                      max(thermo_relaxation_tau,real(dt_seconds))
                 adv_vert(:,:,k)       = 0.0
            endwhere
       enddo
       
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result (convert kg/kg to g/kg)
       q_obs(1,1,:) = tmp_ans(:)
       if (id_qv_obs > 0) then
            used = send_data( id_qv_obs, 1000.*q_obs(:,:,:), time_diag,&
                 1, 1, 1)
       end if  
      
       !---- COMPUTE FINAL ADVECTIVE Q TENDENCY FROM THE SUM OF :  --- !
       !---- RELAXATION TENDENCY AND VERTICAL ADVECTION            --- !                                      --- !
       q_dt(:,:,:,nsphum) = adv_vert(:,:,:) + nudg_q(:,:,:)
            
! ----------------- PRODUCE MOMENTUM FORCING ------------------------- !

       !---- initialize variables
       u_dt       = 0.
       v_dt       = 0.
       
       !---- interpolate u_forc
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=u_forc(1,itime_less)
       field_int_more(1)=u_forc(1,itime_more)
       plev_int(2:(ivert+1))=plev_forc(1:ivert)
       field_int_less(2:(ivert+1))=u_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1))=u_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = u_forc(ivert,itime_less)
       field_int_more(ivert+2) = u_forc(ivert,itime_more)
       tmp_p(:)=pfull(1,1,:)
      
       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
       
       !---- compute damping term
       do k=1,KDIM          
            u_dt(:,:,k)=(tmp_ans(k)-ua(:,:,k))/                         &
                 max(wind_relaxation_tau,real(dt_seconds))
       enddo
      
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
       
       !---- output netcdf result 
       u_obs(1,1,:) = tmp_ans(:)
       if (id_uwnd_obs > 0) then
       used = send_data( id_uwnd_obs, u_obs(:,:,:), time_diag, 1, 1, 1)
       end if  
      
       !---- interpolate v_forc
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=v_forc(1,itime_less)
       field_int_more(1)=v_forc(1,itime_more)
       plev_int(2:(ivert+1))=plev_forc(1:ivert)
       field_int_less(2:(ivert+1))=v_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1))=v_forc(1:ivert,itime_more)
       plev_int(ivert+2) = 0.
       field_int_less(ivert+2) = v_forc(ivert,itime_less)
       field_int_more(ivert+2) = v_forc(ivert,itime_more)
       tmp_p(:)=pfull(1,1,:)

       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
      
       !---- compute damping term
       do k=1,KDIM
            v_dt(:,:,k)=(tmp_ans(k)-va(:,:,k))/                        &
                 max(wind_relaxation_tau,real(dt_seconds))
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
          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nql),&
               adv_vert, scheme = vadvec_scheme,          &
               form = ADVECTIVE_FORM)

          q_dt(:,:,:,nql) = q_dt(:,:,:,nql) + adv_vert    
    
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
          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqi),&
               adv_vert, scheme = vadvec_scheme,          &
               form = ADVECTIVE_FORM)

          q_dt(:,:,:,nqi) = q_dt(:,:,:,nqi) + adv_vert    

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
          call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqa),&
               adv_vert, scheme = vadvec_scheme,          &
               form = ADVECTIVE_FORM)

          q_dt(:,:,:,nqa) = q_dt(:,:,:,nqa) + adv_vert

         if (nqn> 0) then 
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
           call vert_advection(real(dt_seconds),omega_half,delp,q(:,:,:,nqn),&
                adv_vert, scheme = vadvec_scheme,          &
                form = ADVECTIVE_FORM)

           q_dt(:,:,:,nqn) = q_dt(:,:,:,nqn) + adv_vert
         endif
                                        
       end if

! ------------------- COMPUTE SURFACE DATA:  sst --------------------- !

       !---- interpolate sst field
       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=sst_forc(itime_less)
       field_int_more(:)=sst_forc(itime_more)
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
       field_int_more(:) = lwp_forc(itime_more)
       call interp_2d_field(tmp_p,weight_less,weight_more,plev_int,    &
            field_int_less,field_int_more,tmp_ans)
       tmp_res2(:,:) = tmp_ans(1)       
       tmplwpforc = tmp_res2(1,1)
       
       if (id_lwp_obs > 0) then
       used = send_data( id_lwp_obs, tmp_res2(:,:), time_diag, 1, 1)
       end if  
          
       ! ---- interpolate and scale ql ----------------------------- !
       plev_int(1)=ps_init(1,1)
       field_int_less(1)=ql_forc(1,itime_less)
       field_int_more(1)=ql_forc(1,itime_more)
       plev_int(2:(ivert+1))=plev_forc(1:ivert)
       field_int_less(2:(ivert+1))=ql_forc(1:ivert,itime_less)
       field_int_more(2:(ivert+1))=ql_forc(1:ivert,itime_more)
       tmp_p(:)=pfull(1,1,:)
       plev_int(ivert+2)=0.
       field_int_less(ivert+2)=ql_forc(ivert,itime_less)
       field_int_more(ivert+2)=ql_forc(ivert,itime_more)
       tmp_p(:)=pfull(1,1,:)

       !---- do interpolation
       call interp_2d_field(tmp_p,weight_less,weight_more,             &
            plev_int(1:(ivert+2)),field_int_less(1:(ivert+2)),         &
            field_int_more(1:(ivert+2)),tmp_ans)
                   
       tmplwp = 0.
       do k = 1, kdim
               tmplwp = tmplwp + tmp_ans(k)*&
                                 (phalf(1,1,k+1)-phalf(1,1,k))/grav
       enddo

       if (tmplwp .gt. 0. .and. tmplwpforc .gt. 0.) then
            tmp_ans = tmp_ans * tmplwpforc / tmplwp
       end if       
     
       !---- Set to missing above top level of variational analysis
       !
       where (tmp_p .lt. plev_forc(ivert)) tmp_ans = missing_value
              
       !---- output netcdf result (convert kg/kg to g/kg)
       ql_obs(1,1,:) = 1000.*tmp_ans(:)
       if (id_ql_obs > 0)                                              &
       used = send_data( id_ql_obs, ql_obs(:,:,:), time_diag, 1, 1, 1)



   !---- nudging tendencies
   if (id_nudg_t > 0) then
        used = send_data ( id_nudg_t, nudg_t, time_diag, 1, 1, 1)
   end if
 
   if (id_nudg_q > 0) then
        used = send_data ( id_nudg_q, nudg_q, time_diag, 1, 1, 1)
   end if
 
                                                                           
! 
!-----------------------------------------------------------------------
 
end subroutine update_astex_forc


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


type(TIME_TYPE),intent(IN)        :: time_interp
type(TIME_TYPE),intent(IN),dimension(:)        :: TIME_VEC
integer, intent (INOUT)           :: itime_more,itime_less
real, intent (INOUT)              :: weight_more,weight_less

!  Internal variables
!  ------------------

integer                           :: itime_vec,j
integer                           :: seconds_frac,seconds_int
integer                           :: days_frac, days_int
!
! Code
! ----

!findout size of TIME_VEC
itime_vec = size(TIME_VEC,1)

! --- find nearest time to time_interp --- !
      if (time_interp <= TIME_VEC(1) .or. &
          time_interp >= TIME_VEC(itime_vec)) then
          if (time_interp <= TIME_VEC(1)) then
               itime_more = 1
               itime_less = 1
               weight_more = 1.
               weight_less = 0.
          end if
          if (time_interp >= TIME_VEC(itime_vec)) then
               itime_more = itime_vec
               itime_less = itime_vec
               weight_more = 1.
               weight_less = 0.
          end if
      else 
          itime_more = 0
          do j = 2, itime_vec
              if (TIME_VEC(j) >= time_interp) then
              if (itime_more .eq. 0) then
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

real, intent (IN)                 :: weight_more,weight_less
real, intent (IN), dimension(:)   :: p_interp
real, intent (INOUT),dimension(:) :: plev,field_less,field_more
real, intent (INOUT),dimension(:) :: field_interp

!  Internal variables
!  ------------------

real,    dimension(size(field_less,1))  :: field_tmp
real                                    :: ptmp
integer, dimension(size(p_interp,1))    :: klev_more,klev_less
integer                                 :: j,nfield_lev,ninterp
logical, dimension(size(p_interp,1))    :: FIND_PINTERP
real,    dimension(size(p_interp,1))    :: pweight_more,pweight_less

!
! Code
! ----

! ---  count # of levels to the field and # of interpolating levels
nfield_lev = size(field_less,1) 
ninterp = size(p_interp,1)


! --- resort bottom levels if necessary
      if (plev(1).lt.plev(2)) then
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
      FIND_PINTERP(:)=.false.
      pweight_more(:)=0.
      pweight_less(:)=0.

      where(p_interp(:) .ge. plev(1))
            FIND_PINTERP(:)=.true.
            klev_more(:)=1
            klev_less(:)=1
            pweight_more(:)=1.
            pweight_less(:)=0.
      end where

      where(p_interp(:) .le. plev(nfield_lev))
            FIND_PINTERP(:)=.true.
            klev_more(:)=nfield_lev
            klev_less(:)=nfield_lev
            pweight_more(:)=1.
            pweight_less(:)=0.
      end where

      do j = 2, nfield_lev
         where(plev(j) .le. p_interp(:) .and. .not. FIND_PINTERP(:))
            FIND_PINTERP(:)=.true.
            klev_more(:)=j
            klev_less(:)=j-1
            pweight_more(:)=(p_interp(:) - plev(klev_less))/&
                            (plev(klev_more) - plev(klev_less))
            pweight_less(:)=1.-pweight_more(:)
         end where
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

real,     intent (IN)                     :: bad_val
type(TIME_TYPE), intent (IN), dimension(:):: TIME_VEC
real,     intent (INOUT), dimension(:,:)  :: data_mat
integer,  optional, intent (IN)           :: adjlev

!  Internal variables
!  ------------------

type(TIME_TYPE)                          :: time_less,time_more
integer                                  :: itime_less,itime_more
integer                                  :: jlev_less,jlev_more,jcur
integer                                  :: dpointslev,dpointstime
integer                                  :: nlev,ntime,jlev,itime,icur,itime2
real                                     :: weight_less,weight_more
integer                                  :: seconds_frac,seconds_int
integer                                  :: days_frac, days_int
real, dimension(size(data_mat,1),size(data_mat,2)) :: data_mat_new
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
       if (data_mat(jlev,itime).le.bad_val) then

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


!########################################################################
! This subroutine returns imposed SST

subroutine get_astex_sst( sst )

  implicit none
  real, intent(out) :: sst

  sst = tskin

end subroutine get_astex_sst

!########################################################################
!########################################################################


end module scm_astex_mod
 

     


