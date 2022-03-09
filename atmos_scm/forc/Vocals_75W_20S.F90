module scm_Vocals_75W_20S_mod
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       Vocals_75W_20S (Oct.-Nov., 2008) FORCING MODULE
!
!       Dec 2011
!       Contact person: Huan Guo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This module read initial conditions and forcing data from CAPT2 simulations
!       References
!       ----------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   use            mpp_mod, only:  mpp_pe, mpp_root_pe, stdlog
   use            fms_mod, only:  write_version_number, open_namelist_file,  &
                                  check_nml_error,      file_exist,          &
                                  error_mesg,           close_file,   FATAL, &
                                  NOTE,                 read_data,           &
                                  mpp_error

#ifdef INTERNAL_FILE_NML
   use              mpp_mod, only: input_nml_file
#else
   use              fms_mod, only: open_namelist_file
#endif

   use   diag_manager_mod, only:  register_diag_field, send_data
   use   time_manager_mod
   use  scm_horiz_grid_mod,only:  horiz_grid_type
   use      constants_mod, only:  kappa
   use   field_manager_mod, only: MODEL_ATMOS
   use  tracer_manager_mod, only: get_tracer_index

use             fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                               endlon, rlonb, rlatb,  cold_start, ncnst, &
                               pnats, consv_te, ptop, fv_init, fv_domain, &
                               fv_end, change_time, p_var, restart_format, area, &
                               ak, bk, rlon, rlat, ng_d, f_d, nt_prog, get_eta_level

   implicit none
   private

   public Vocals_75W_20S_data_read, Vocals_75W_20S_forc_init, Vocals_75W_20S_forc_end, update_Vocals_75W_20S_forc, &
          Vocals_75W_20S_forc_diagnostic_init,        get_Vocals_75W_20S_flx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following subroutines:
!
!            Vocals_75W_20S_data_read      reads forcing data files and initializes
!                                 any needed parameters of the run
!            Vocals_75W_20S_forc_init      initializes prognostic variables:           
!                                 T,u,v,qv,ql,qi
!            Vocals_75W_20S_forc_end       deallocates allocated spaces
!            update_Vocals_75W_20S_forc    a call to this routine returns the 
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

character(len=14) :: mod_name = 'Vocals_75W_20S'
character(len=7)  :: mod_name_diag = 'forcing'

!----------Diagnostic data----------------------------------------------
integer ::  id_tdt_dyn,    id_udt_dyn,  id_vdt_dyn
integer ::  id_qadt_dyn,   id_qdt_dyn,  id_qldt_dyn,  id_qidt_dyn
integer ::  id_qndt_dyn,   id_qnidt_dyn

!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       PARAMETERS OF THE MODULE
!       tskin                    sea surface temperature (K)
!       SENFLUX                  sensible heat flux (W/m2)
!       EVAPFLUX                 evaporation flux (kg water/m2/sec)

REAL,    PRIVATE               :: missing_value = -999.
REAL,    PUBLIC                :: tskin = 257.4
REAL,    PUBLIC                :: SENFLUX = 7.98  !(w/m2)
REAL,    PUBLIC                :: EVAPFLUX = 2.86 !(w/m2)

real,    private               :: zsfc = 0.0
real,    private               :: p00  = 1000.0  !reference pressure: 1000 mb

integer, dimension(6)          :: forc_begin_date=(/ 0, 0, 0, 0, 0, 0 /)
type(TIME_TYPE)                :: Time_step_forc
integer                        :: dt_forc = 1800    !  time interval of forcing data (sec)

integer                        :: n_forc_time   = 17520
integer                        :: n_forc_levels = 48
integer                        :: n_forc_lat    = 1
integer                        :: n_forc_lon    = 1


type(TIME_TYPE), allocatable,  dimension( : )                :: time_forc
real,            allocatable,  dimension( : )                :: plev_forc
real,            allocatable,  dimension( :, :, :, : )       :: t_profile,    &
                                                                qa_profile,   &
                                                                qv_profile,   &
                                                                ql_profile,   &
                                                                qi_profile,   &
                                                                nqn_profile,  &
                                                                nqni_profile, &
                                                                u_profile,    &
                                                                v_profile

real,            allocatable,  dimension( :, :, :, : )       :: t_forc,       &
                                                                qa_forc,      &
                                                                qv_forc,      &
                                                                ql_forc,      &
                                                                qi_forc,      &
                                                                nqn_forc,     &
                                                                nqni_forc,    &
                                                                u_forc,       &
                                                                v_forc

real,            allocatable,  dimension( :, :, : )          :: ps_forc,      &
                                                                u_star_forc,  &
                                                                shflx_forc,   &
                                                                evap_forc

integer, save                                                :: itime_less,   itime_more
real, save                                                   :: weight_less,  weight_more

real,            allocatable,  dimension( : )                :: ph,  lph
real,            allocatable,  dimension( : )                :: pf


logical, public                :: use_const_forcing = .true.
integer, dimension(6)          :: avg_t1 = (/ 0, 0, 0, 0, 0, 0 /)
integer, dimension(6)          :: avg_t2 = (/ 0, 0, 0, 0, 0, 0 /)
TYPE(TIME_TYPE)                :: Time_t1,  Time_t2

logical, public                :: do_nudging = .false.
real,    public                :: tao_nudge = 21600.0

logical, public                :: do_geo = .true.

character*64                   :: fname_forc='INPUT/19810101.atmos_col_004.nc'

LOGICAL, PRIVATE               :: Vocals_75W_20S_forc_initialized = .FALSE.

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!
character(len=128) :: Version = '$Id: Vocals_75W_20S.F90,v 1.1.2.1 2014/02/10 18:17:23 wfc Exp $'
character(len=128) :: Tagname = '$Name: new_gate_wfc $'
integer, dimension(1) :: restart_versions = (/ 1 /)
        
NAMELIST /scm_Vocals_75W_20S_nml/   fname_forc, forc_begin_date,  dt_forc,                      &
                                    n_forc_lon, n_forc_lat,       n_forc_levels, n_forc_time,   &
                                    do_nudging, tao_nudge,                                      &
                                    use_const_forcing,   avg_t1,     avg_t2,                    &
                                    do_geo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!#######################################################################

subroutine Vocals_75W_20S_data_read( )

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
!       unit                 unit number for I/O file
!       io,ierr              dummy integer variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Internal variables
!  ------------------
INTEGER                        :: unit,io,ierr
integer                        :: sec0, day0
integer                        :: i_forc_time

       ! if this subroutine is already called, then return, otherwise 
       ! set Vocals_75W_20S_forc_initialized .TRUE.

       if ( Vocals_75W_20S_forc_initialized ) return
       Vocals_75W_20S_forc_initialized = .TRUE.

!    --------- read namelist --------
       if (file_exist('input.nml')) then
            unit = open_namelist_file()
            ierr=1; do while (ierr /= 0)
            read  (unit, nml=scm_Vocals_75W_20S_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'scm_Vocals_75W_20S_nml')
            enddo
10          call close_file (unit)
       endif

 !--------- write version number and namelist --------
       call write_version_number ( version, tagname )
       if(mpp_pe() == mpp_root_pe() ) write(stdlog(),nml=scm_Vocals_75W_20S_nml)


       if ( file_exist( fname_forc ) ) then
          if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_Vocals_75W_20S_mod', &
         'Reading netCDF forcing data from: ',   NOTE)

!
!      ALLOCATE STORAGE
!      
         if ( allocated( time_forc ) )   deallocate ( time_forc )
         allocate( time_forc(n_forc_time) )

         if ( allocated( plev_forc ) )   deallocate ( plev_forc )
         allocate( plev_forc(n_forc_levels) )

         if ( allocated( t_profile ) )   deallocate ( t_profile )
         allocate( t_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( qa_profile ) )   deallocate ( qa_profile )
         allocate( qa_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( qv_profile ) )   deallocate ( qv_profile )
         allocate( qv_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( ql_profile ) )   deallocate ( ql_profile )
         allocate( ql_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( qi_profile ) )   deallocate ( qi_profile )
         allocate( qi_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( nqn_profile ) )   deallocate ( nqn_profile )
         allocate( nqn_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( nqni_profile ) )   deallocate ( nqni_profile )
         allocate( nqni_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( u_profile ) )   deallocate ( u_profile )
         allocate( u_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( v_profile ) )   deallocate ( v_profile )
         allocate( v_profile( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )

         if ( allocated( t_forc ) )   deallocate ( t_forc )
         allocate( t_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( qa_forc ) )   deallocate ( qa_forc )
         allocate( qa_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( qv_forc ) )   deallocate ( qv_forc )
         allocate( qv_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( ql_forc ) )   deallocate ( ql_forc )
         allocate( ql_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( qi_forc ) )   deallocate ( qi_forc )
         allocate( qi_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( nqn_forc ) )   deallocate ( nqn_forc )
         allocate( nqn_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( nqni_forc ) )   deallocate ( nqni_forc )
         allocate( nqni_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( u_forc ) )   deallocate ( u_forc )
         allocate( u_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )
         if ( allocated( v_forc ) )   deallocate ( v_forc )
         allocate( v_forc( n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time ) )


         if ( allocated( ps_forc ) )   deallocate ( ps_forc )
         allocate( ps_forc( n_forc_lon, n_forc_lat, n_forc_time ) )
         if ( allocated( u_star_forc ) )   deallocate ( u_star_forc )
         allocate( u_star_forc( n_forc_lon, n_forc_lat, n_forc_time ) )
         if ( allocated( shflx_forc ) )   deallocate ( shflx_forc )
         allocate( shflx_forc( n_forc_lon, n_forc_lat, n_forc_time ) )
         if ( allocated( evap_forc ) )   deallocate ( evap_forc )
         allocate( evap_forc( n_forc_lon, n_forc_lat, n_forc_time ) )

         t_profile    = 0.0
         qa_profile   = 0.0
         qv_profile   = 0.0
         ql_profile   = 0.0
         qi_profile   = 0.0
         nqn_profile  = 0.0
         nqni_profile = 0.0
         u_profile    = 0.0
         v_profile    = 0.0

         t_forc       = 0.0
         qa_forc      = 0.0
         qv_forc      = 0.0
         ql_forc      = 0.0
         qi_forc      = 0.0
         nqn_forc     = 0.0
         nqni_forc    = 0.0
         u_forc       = 0.0
         v_forc       = 0.0

         ps_forc      = 0.0
         u_star_forc  = 0.0
         shflx_forc   = 0.0
         evap_forc    = 0.0

         call read_data (fname_forc, 'pfull',    plev_forc( 1:n_forc_levels), no_domain=.true. )

         Time_step_forc  = set_time ( dt_forc, 0 )
         Time_forc( 1 )  = set_date(forc_begin_date(1), forc_begin_date(2), forc_begin_date(3), forc_begin_date(4), forc_begin_date(5), forc_begin_date(6) )
         do  i_forc_time = 1, n_forc_time
           if ( i_forc_time > 1 ) Time_forc( i_forc_time ) = Time_forc( i_forc_time -1 ) + Time_step_forc
           call read_data (fname_forc, 'temp',      t_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),     timelevel=i_forc_time   )
           call read_data (fname_forc, 'cld_amt',   qa_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),    timelevel=i_forc_time   )
           call read_data (fname_forc, 'sphum',     qv_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),    timelevel=i_forc_time   )
           call read_data (fname_forc, 'liq_wat',   ql_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),    timelevel=i_forc_time   )
           call read_data (fname_forc, 'ice_wat',   qi_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),    timelevel=i_forc_time   )
           call read_data (fname_forc, 'liq_drp',   nqn_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),   timelevel=i_forc_time   )
           call read_data (fname_forc, 'ucomp',     u_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),     timelevel=i_forc_time   )
           call read_data (fname_forc, 'vcomp' ,    v_profile(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),     timelevel=i_forc_time   )

           call read_data (fname_forc, 'tdt_dyn',   t_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),        timelevel=i_forc_time   )
           call read_data (fname_forc, 'qdt_dyn',   qv_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),       timelevel=i_forc_time   )
           call read_data (fname_forc, 'qadt_dyn',  qa_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),       timelevel=i_forc_time   )
           call read_data (fname_forc, 'qldt_dyn',  ql_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),       timelevel=i_forc_time   )
           call read_data (fname_forc, 'qidt_dyn',  qi_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),       timelevel=i_forc_time   )
           call read_data (fname_forc, 'qndt_dyn',  nqn_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),      timelevel=i_forc_time   )
           call read_data (fname_forc, 'qnidt_dyn', nqni_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),     timelevel=i_forc_time   )
           call read_data (fname_forc, 'udt_dyn',   u_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),        timelevel=i_forc_time   )
           call read_data (fname_forc, 'vdt_dyn',   v_forc(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time),        timelevel=i_forc_time   )

           call read_data (fname_forc, 'ps',        ps_forc(:,:,i_forc_time),         timelevel=i_forc_time   )

           call read_data (fname_forc, 'u_star',    u_star_forc(:,:,i_forc_time),     timelevel=i_forc_time   )
           call read_data (fname_forc, 'shflx',     shflx_forc(:,:,i_forc_time),      timelevel=i_forc_time   )
           call read_data (fname_forc, 'evap',      evap_forc(:,:,i_forc_time),       timelevel=i_forc_time   )

         enddo

         if ( use_const_forcing ) then
            Time_t1 =set_date(avg_t1(1), avg_t1(2), avg_t1(3), avg_t1(4), avg_t1(5), avg_t1(6) )
            Time_t2 =set_date(avg_t2(1), avg_t2(2), avg_t2(3), avg_t2(4), avg_t2(5), avg_t2(6) )

            call avg_forc ( t_forc )
            call avg_forc ( qv_forc )
            call avg_forc ( qa_forc )
            call avg_forc ( ql_forc )
            call avg_forc ( qi_forc )
            call avg_forc ( nqn_forc )
            call avg_forc ( nqni_forc )
         endif  ! use_const_forcing

       endif     !file_exist('INPUT/*.atmos_col_004.nc' )
!-----------------------------------------------------------------------
                            
end subroutine Vocals_75W_20S_data_read


!#######################################################################
subroutine Vocals_75W_20S_forc_init(time_interp, As, Bs)
#include "fv_arrays.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a given time and initializes the model 
!      variables for that time. This involves interpolating the global 
!      storage fields to a given time and pressure.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!      k               counting integers
!      KDIM            no. of vertical levels to T array
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
integer                                          :: kdim, k
real, dimension( n_forc_levels )                 :: field_int_less, field_int_more, plev_int
real, dimension( size(pt,3) )                     :: field_interp_out, tmp_p

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

       ! --- find out time in sfc pressure field
       call interp_time(time_interp,time_forc,itime_less,itime_more,    &
                        weight_less,weight_more)
       print*, ' init_profile, itime_less = ', itime_less,  weight_less
       print*, ' init_profile, itime_more = ', itime_more,  weight_more

       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=ps_forc(1,1,itime_less)
       field_int_more(:)=ps_forc(1,1,itime_more)
       call interp_2d_field( tmp_p, weight_less,weight_more, plev_int,    &
            field_int_less, field_int_more, field_interp_out(1:KDIM))
       ps(:,:) = field_interp_out(1)
       elev = zsfc
       print*, 'initial ps= ',  ps(:,:)

! --- Create half pressure levels
       if ( allocated( ph ) )    deallocate ( ph )
       allocate( ph( KDIM+1 ) )
       if ( allocated( lph ) )   deallocate ( lph )
       allocate( lph( KDIM+1 ) )
       if ( allocated( pf ) )    deallocate ( pf )
       allocate( pf( KDIM ) )


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
           call get_eta_level(nlev, ps(i,j) , pf(:), ph(:))
         enddo
       enddo

       field_int_less(1:n_forc_levels) = t_profile(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = t_profile(1,1, 1:n_forc_levels, itime_more)

       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       ! --- assign T
       do k = 1, KDIM
            pt(:,:,k)= field_interp_out(k)
       enddo

      field_int_less(1:n_forc_levels) = qv_profile(1,1, 1:n_forc_levels, itime_less) &
                                       +ql_profile(1,1, 1:n_forc_levels, itime_less) &
                                       +qi_profile(1,1, 1:n_forc_levels, itime_less)
      field_int_more(1:n_forc_levels) = qv_profile(1,1, 1:n_forc_levels, itime_more) &
                                       +ql_profile(1,1, 1:n_forc_levels, itime_more) &
                                       +qi_profile(1,1, 1:n_forc_levels, itime_more)

       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )

       do k = 1, KDIM
            q(:,:,k,nsphum)= field_interp_out(k)
       enddo
       q(:, :, :,nqa)  = 0.0
       q(:, :, :,nql)  = 0.0
       q(:, :, :,nqi)  = 0.0

       if(nqn  > 0) q(:, :, :,nqn)  =  0.0
       if(nqni > 0) q(:, :, :,nqni) =  0.0

       field_int_less(1:n_forc_levels) = u_profile(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = u_profile(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
            ua(:,:,k)= field_interp_out(k)
       enddo

       field_int_less(1:n_forc_levels) = v_profile(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = v_profile(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
            va(:,:,k)= field_interp_out(k)
       enddo
!-----------------------------------------------------------------------
 end subroutine Vocals_75W_20S_forc_init



!#######################################################################
subroutine Vocals_75W_20S_forc_end()

  if (.NOT.Vocals_75W_20S_forc_initialized) return
  Vocals_75W_20S_forc_initialized = .FALSE.

end subroutine Vocals_75W_20S_forc_end

!#######################################################################


subroutine Vocals_75W_20S_forc_diagnostic_init(axes, Time)

implicit none

integer, dimension(3) :: half = (/1,2,4/)
integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time
! --- initialize axes -------------------------------------------------!

  id_tdt_dyn =  register_diag_field (mod_name_diag, 'tdt_dyn', axes(1:3), Time, &
    'temperatue tendency from dynamics', 'K/s',  missing_value = missing_value)

  id_udt_dyn =  register_diag_field (mod_name_diag, 'udt_dyn', axes(1:3), Time, &
    'u-wind tendency from dynamics', 'm/s/s',  missing_value = missing_value)

  id_vdt_dyn =  register_diag_field (mod_name_diag, 'vdt_dyn', axes(1:3), Time, &
    'v-wind tendency from dynamics', 'm/s/s',  missing_value = missing_value)

  id_qadt_dyn =  register_diag_field (mod_name_diag, 'qadt_dyn', axes(1:3), Time, &
    'cloud amount tendency from dynamics', '/s',  missing_value = missing_value)

  id_qdt_dyn =  register_diag_field (mod_name_diag, 'qdt_dyn', axes(1:3), Time, &
    'water vapor tendency from dynamics', 'Kg/Kg/s',  missing_value = missing_value)

  id_qldt_dyn =  register_diag_field (mod_name_diag, 'qldt_dyn', axes(1:3), Time, &
    'liquid water tendency from dynamics', 'Kg/Kg/s',  missing_value = missing_value)

  id_qidt_dyn =  register_diag_field (mod_name_diag, 'qidt_dyn', axes(1:3), Time, &
    'ice water tendency from dynamics', 'Kg/Kg/s',  missing_value = missing_value)

  id_qndt_dyn =  register_diag_field (mod_name_diag, 'qndt_dyn', axes(1:3), Time, &
    'droplet number tendency from dynamics', '/Kg/s',  missing_value = missing_value)

  id_qnidt_dyn =  register_diag_field (mod_name_diag, 'qnidt_dyn', axes(1:3), Time, &
    'ice number tendency from dynamics', '/Kg/s',  missing_value = missing_value)
!-----------------------------------------------------------------------
end subroutine Vocals_75W_20S_forc_diagnostic_init


!#######################################################################
subroutine update_Vocals_75W_20S_forc(time_interp,time_diag,dt_int)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------
!

TYPE(TIME_TYPE), intent(in)              :: time_interp,time_diag,dt_int
!  Internal variables
!  ------------------
integer                                            :: kdim, k
real, dimension( n_forc_levels )                   :: field_int_less, field_int_more, plev_int
real, dimension( size(pt,3) )                       :: field_interp_out, tmp_p
real                                               :: fcriolis

real, dimension(size(pt,1),size(pt,2),size(pt,3))     :: du_geos, dv_geos, u_geos, v_geos
logical                                            :: used
! ------------------------------------------------------------------------------------------------------------

       KDIM = size(pt,3)

       call interp_time(time_interp,time_forc,itime_less,itime_more,    &
                        weight_less,weight_more)
       print*, ' update_forc, itime_less = ', itime_less,  weight_less
       print*, ' update_forc, itime_more = ', itime_more,  weight_more

       tmp_p(:)=50000.
       plev_int(:)=p00
       field_int_less(:)=ps_forc(1,1,itime_less)
       field_int_more(:)=ps_forc(1,1,itime_more)

       call interp_2d_field( tmp_p, weight_less,weight_more, plev_int,    &
            field_int_less, field_int_more, field_interp_out(1:KDIM))
       ps(:,:) = field_interp_out(1)

       field_int_less(1:n_forc_levels) = t_forc(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = t_forc(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
            t_dt(:,:,k)= field_interp_out(k)
       enddo

       if ( do_nudging ) then
         field_int_less(1:n_forc_levels) = u_profile(1,1, 1:n_forc_levels, itime_less)
         field_int_more(1:n_forc_levels) = u_profile(1,1, 1:n_forc_levels, itime_more)
         call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
         do k = 1, KDIM
            u_dt(:,:,k) = - ( u(:,:,k) - field_interp_out(k) ) / tao_nudge
         enddo

         field_int_less(1:n_forc_levels) = v_profile(1,1, 1:n_forc_levels, itime_less)
         field_int_more(1:n_forc_levels) = v_profile(1,1, 1:n_forc_levels, itime_more)
         call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
         do k = 1, KDIM
            v_dt(:,:,k) = - ( V(:,:,k) - field_interp_out(k) ) / tao_nudge
         enddo

       else
         field_int_less(1:n_forc_levels) = u_forc(1,1, 1:n_forc_levels, itime_less)
         field_int_more(1:n_forc_levels) = u_forc(1,1, 1:n_forc_levels, itime_more)
         call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
         do k = 1, KDIM
            u_dt(:,:,k)= field_interp_out(k)
         enddo

         field_int_less(1:n_forc_levels) = v_forc(1,1, 1:n_forc_levels, itime_less)
         field_int_more(1:n_forc_levels) = v_forc(1,1, 1:n_forc_levels, itime_more)
         call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
         do k = 1, KDIM
            v_dt(:,:,k)= field_interp_out(k)
         enddo
       endif

! --- compute geostrophic tendencies
       du_geos=0.
       dv_geos=0.
       if (do_geo) then
           fcriolis=f_d(1,1)

           field_int_less(1:n_forc_levels) = u_profile(1,1, 1:n_forc_levels, itime_less)
           field_int_more(1:n_forc_levels) = u_profile(1,1, 1:n_forc_levels, itime_more)
           call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
           do k = 1, KDIM
             u_geos(:,:,k) = field_interp_out(k)
           enddo
           field_int_less(1:n_forc_levels) = v_profile(1,1, 1:n_forc_levels, itime_less)
           field_int_more(1:n_forc_levels) = v_profile(1,1, 1:n_forc_levels, itime_more)
           call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
           do k = 1, KDIM
             v_geos(:,:,k) = field_interp_out(k)
           enddo

           do k=1, KDIM
              du_geos(:,:,k) =   fcriolis * (v(1,1,k)-v_geos(1,1,k))
              dv_geos(:,:,k) = - fcriolis * (u(1,1,k)-u_geos(1,1,k))
              u_dt(:,:,k)      = u_dt(:,:,k) + du_geos(:,:,k)
              v_dt(:,:,k)      = v_dt(:,:,k) + dv_geos(:,:,k)
            end do
       end if

       field_int_less(1:n_forc_levels) = qv_forc(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = qv_forc(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
           q_dt(:,:,k,nsphum)= field_interp_out(k)
       enddo

       field_int_less(1:n_forc_levels) = ql_forc(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = ql_forc(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
           q_dt(:,:,k,nql)= field_interp_out(k)
       enddo

       field_int_less(1:n_forc_levels) = qi_forc(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = qi_forc(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
           q_dt(:,:,k,nqi)= field_interp_out(k)
       enddo

       field_int_less(1:n_forc_levels) = qa_forc(1,1, 1:n_forc_levels, itime_less)
       field_int_more(1:n_forc_levels) = qa_forc(1,1, 1:n_forc_levels, itime_more)
       call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
            field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
       do k = 1, KDIM
           q_dt(:,:,k,nqa)= field_interp_out(k)
       enddo

       if( nqn > 0 ) then
         field_int_less(1:n_forc_levels) = nqn_forc(1,1, 1:n_forc_levels, itime_less)
         field_int_more(1:n_forc_levels) = nqn_forc(1,1, 1:n_forc_levels, itime_more)
         call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
         do k = 1, KDIM
           q_dt(:,:,k,nqn)= field_interp_out(k)
         enddo
       endif

       if( nqni > 0 ) then
         field_int_less(1:n_forc_levels) = nqni_forc(1,1, 1:n_forc_levels, itime_less)
         field_int_more(1:n_forc_levels) = nqni_forc(1,1, 1:n_forc_levels, itime_more)
         call interp_2d_field( pf(1:KDIM), weight_less, weight_more, plev_forc(1:n_forc_levels),    &
              field_int_less(1:n_forc_levels), field_int_more (1:n_forc_levels), field_interp_out(1:KDIM) )
         do k = 1, KDIM
           q_dt(:,:,k,nqni)= field_interp_out(k)
         enddo
       endif

       omga  = 0.0

       if ( id_tdt_dyn > 0 )        used = send_data(  id_tdt_dyn, t_dt(:,:,:), Time_diag, 1, 1 )

       if ( id_udt_dyn > 0 )        used = send_data(  id_udt_dyn, u_dt(:,:,:), time_diag, 1, 1 )

       if ( id_vdt_dyn > 0 )        used = send_data(  id_vdt_dyn, v_dt(:,:,:), time_diag, 1, 1 )

       if ( id_qadt_dyn > 0 )       used = send_data(  id_qadt_dyn, q_dt(:,:,:,nqa), time_diag, 1, 1 )

       if ( id_qdt_dyn > 0 )        used = send_data(  id_qdt_dyn, q_dt(:,:,:,nsphum), time_diag, 1, 1 )

       if ( id_qldt_dyn > 0 )       used = send_data(  id_qldt_dyn, q_dt(:,:,:,nql), time_diag, 1, 1 )

       if ( id_qidt_dyn > 0 )       used = send_data(  id_qidt_dyn, q_dt(:,:,:,nqi), time_diag, 1, 1 )

       if( (nqn > 0) .and. id_qndt_dyn > 0)      used = send_data(  id_qndt_dyn,  q_dt(:,:,:,nqn), time_diag, 1, 1 )

       if( (nqni > 0) .and. id_qnidt_dyn > 0)    used = send_data(  id_qnidt_dyn, q_dt(:,:,:,nqni), time_diag, 1, 1 )

! the surface fluxes are called before update_Vocals_75W_20S_forc, 
! so the next time index is calculated for the surface fluxes at next time step!
       call interp_time(time_interp+dt_int,time_forc,itime_less,itime_more,    &
                        weight_less,weight_more)
       print*, ' update_forc, itime_less = ', itime_less,  weight_less
       print*, ' update_forc, itime_more = ', itime_more,  weight_more
!-----------------------------------------------------------------------
end subroutine update_Vocals_75W_20S_forc



!########################################################################
! This subroutine returns imposed surface fluxes

subroutine get_Vocals_75W_20S_flx( ustar, flux_t, flux_q )

  implicit none
  real, intent(out), dimension(:) :: ustar, flux_t, flux_q
 ! ustar: m/s;    flux_t: w/m2;    flux_q: kg/m2/s
!-----------------------------------------------------------------------

  print*, ' sfc_flux, itime_less = ', itime_less,  weight_less
  print*, ' sfc_flux, itime_more = ', itime_more,  weight_more
  ustar  = u_star_forc(1, 1, itime_less) * weight_less + u_star_forc(1, 1, itime_more) * weight_more
  flux_t = shflx_forc( 1, 1, itime_less) * weight_less + shflx_forc( 1, 1, itime_more) * weight_more
  flux_q = evap_forc(  1, 1, itime_less) * weight_less + evap_forc(  1, 1, itime_more) * weight_more


end subroutine get_Vocals_75W_20S_flx
!-----------------------------------------------------------------------


!########################################################################
! This subroutine calculates average forcing

subroutine  avg_forc( forc_ls )
  real, INTENT (INout),   dimension(n_forc_lon, n_forc_lat, n_forc_levels, n_forc_time)  ::  forc_ls

  integer                        :: sec0, day0
  integer                        :: i_forc_time
  integer                        :: iforc_t1, iforc_t2

  real,    dimension(n_forc_lon, n_forc_lat, n_forc_levels)       ::  forc_tmp
!-----------------------------------------------------------------------
  call get_time( Time_t1 - Time_forc(1), sec0, day0)
  iforc_t1 = (day0*86400 + sec0) / dt_forc + 1
  iforc_t1 = min ( iforc_t1, n_forc_time)

  call get_time( Time_t2 - Time_forc(1), sec0, day0)
  iforc_t2 = (day0*86400 + sec0) / dt_forc + 1
  iforc_t2 = min ( iforc_t2, n_forc_time)

  if ( iforc_t2 >= iforc_t1 .and. iforc_t1 > 0 ) then
   forc_tmp = 0.0
   do  i_forc_time = iforc_t1, iforc_t2
     forc_tmp(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels) =  forc_tmp(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels) &
                                                            + forc_ls(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time)
   enddo
   forc_tmp(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels) = &
           forc_tmp(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels) / real( iforc_t2 - iforc_t1 + 1)
   do  i_forc_time = 1, n_forc_time
       forc_ls(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels, i_forc_time ) = forc_tmp(1:n_forc_lon, 1:n_forc_lat, 1:n_forc_levels)
   enddo
  else
   if ( iforc_t1 <= 0 )      print*, 'averaging forcing t1 < initial forcing data,  use time-dependent forcing'
   if ( iforc_t2 < iforc_t1) print*, 'averaging forcing t2 < t1,  use time-dependent forcing', iforc_t2, iforc_t1
  endif !  iforc_t2 >= iforc_t1 .and. iforc_t1 > 0

end subroutine   avg_forc
!-----------------------------------------------------------------------
!########################################################################





! copy from scm_arm interpolations
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




! copy from scm_arm interpolations. Since the pressure levels in the ARM forcing data 
! increases upwards (i.e., plev(1) > plev(2) > plev(3)), and 
! they are increases downwards in GFDL models, I made slight changes.
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
real, intent (IN),dimension(:)    :: plev,field_less,field_more
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

! --- interpolate field in time
field_tmp = weight_more * field_more + weight_less*field_less


! --- find pressure levels nearest to p_interp --- !
      !set find flag to false and initialize weights
      FIND_PINTERP(:)=.false.
      pweight_more(:)=0.
      pweight_less(:)=0.

      where(p_interp(:) .ge. plev(nfield_lev))
            FIND_PINTERP(:)=.true.
            klev_more(:)= nfield_lev
            klev_less(:)= nfield_lev
            pweight_more(:)=1.
            pweight_less(:)=0.
      end where

      where(p_interp(:) .le. plev(1))
            FIND_PINTERP(:)=.true.
            klev_more(:)= 1
            klev_less(:)= 1
            pweight_more(:)=1.
            pweight_less(:)=0.
      end where

      do j = 2, nfield_lev
         where(plev(j) .ge. p_interp(:) .and. .not. FIND_PINTERP(:))
            FIND_PINTERP(:)=.true.
            klev_more(:)=j
            klev_less(:)=j-1
            pweight_more(:)=(p_interp(:) - plev(klev_less))/&
                            (plev(klev_more) - plev(klev_less))
            pweight_less(:)=1.-pweight_more(:)
         end where
      enddo  !for j loop  
! <-- h1g, 2012-01-12

! --- do interpolation
    field_interp(:) = pweight_more(:)*field_tmp(klev_more(:)) + &
                      pweight_less(:)*field_tmp(klev_less(:))

! 
!-----------------------------------------------------------------------
 
end subroutine interp_2d_field 

!########################################################################
!########################################################################

end module scm_Vocals_75W_20S_mod
