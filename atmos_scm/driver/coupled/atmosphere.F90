module atmosphere_mod
#include <fms_platform.h>
!-----------------------------------------------------------------------
!
!         interface for b-grid dynamical core and physics
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use scm_core_driver_mod, only: scm_core_driver_init, &
                               scm_core_driver,      &
                               scm_core_driver_end

use    scm_prog_var_mod, only: prog_var_type, prog_var_init,  &
                               var_init, prog_var_time_diff,  &
                               prog_var_times_scalar

use scm_diagnostics_mod, only: scm_diagnostics,      &
                               scm_diagnostics_init
use        scm_forc_mod, only: scm_forc_diagnostic_init, scm_data_read, &
                               update_scm_forc
use  scm_grid_masks_mod, only: grid_mask_type, grid_masks_init


use    time_manager_mod, only: time_type, get_time, set_time, &
                               operator(+), operator(-)

use             fms_mod, only: error_mesg, FATAL, stdlog,  &
                               write_version_number,       &
                               mpp_pe, mpp_root_pe,        &
                               mpp_npes,                   &
                               file_exist, field_size,     &
                               open_restart_file,          &
                               check_nml_error,            &
                               close_file, set_domain,     &
                               write_data, field_exist,    &
                               read_data, mpp_error, NOTE,     &
                               mpp_clock_id, mpp_clock_begin,  &
                               mpp_clock_end, CLOCK_SUBCOMPONENT, &
                               clock_flag_default, nullify_domain

#ifdef INTERNAL_FILE_NML
use             mpp_mod, only: input_nml_file, mpp_get_current_pelist
#else
use             fms_mod, only: open_namelist_file
#endif


use     mpp_domains_mod, only: domain2d
use   field_manager_mod, only: MODEL_ATMOS
use diag_manager_mod,   only: register_diag_field, send_data
use       constants_mod, only: rdgas, grav, rvgas, radius, PI
 
use       constants_mod, only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use  tracer_manager_mod, only: get_number_tracers, get_tracer_index, &
                               get_tracer_names, NO_TRACER
 
use      topography_mod, only: get_topog_mean,get_ocean_mask

use           xgrid_mod, only: grid_box_type

use   block_control_mod, only: block_control_type
use  physics_driver_mod, only: surf_diff_type
use   physics_types_mod, only: physics_type, &
                               physics_tendency_type
use radiation_types_mod, only: radiation_type, compute_g_avg
use atmos_cmip_diag_mod, only: atmos_cmip_diag_init
use atmos_global_diag_mod, only: atmos_global_diag_init, &
                                 atmos_global_diag_end

!-----------------
! FV core modules:
!-----------------
use            fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                              endlon, rlonb, rlatb,  cold_start, ncnst, &
                              pnats, consv_te, ptop, fv_init, fv_domain, &
                              fv_end, change_time, p_var, restart_format, area, &
                              ak, bk, rlon, rlat, ng_d, nt_prog, get_eta_level, &
                              esm2_bugs, target_lon, target_lat, SCM_RES_X, SCM_RES_Y
use fv_diagnostics,     only: fv_diag_init, fv_diag, fv_time
use timingModule,       only: timing_on, timing_off
use fv_restart_mod,     only: fv_restart, write_fv_rst
use atmos_nudge_mod,    only: atmos_nudge_init, atmos_nudge_end
use fv_arrays_mod,      only: fv_print_chksums
use update_fv_phys_mod, only: update_fv_phys

!yzheng
use five_mod, only: ua_five, va_five, pt_five, q_five, &
                    u_dt_five, v_dt_five, t_dt_five, q_dt_five

!-----------------------------------------------------------------------

implicit none
private

!--- driver routines
public :: atmosphere_init,       atmosphere_end,  atmosphere_restart,      &
          atmosphere_dynamics, atmosphere_state_update

!--- utility routines
!public :: atmosphere_down,        atmosphere_up,           &
public :: atmosphere_resolution,  atmosphere_boundary,     &
          atmosphere_grid_center, atmosphere_domain,       &
          atmosphere_cell_area,   atmosphere_control_data, &
          atmosphere_pref,                                 &
          get_atmosphere_axes,    get_bottom_mass,         &
          get_bottom_wind,        get_stock_pe,            &
          set_atmosphere_pelist,  reset_atmos_tracers

!--- physics/radiation data exchange routines
public :: atmos_radiation_driver_inputs, atmos_physics_driver_inputs

!yzheng
public :: atmosphere_state_update_five

public  surf_diff_type

integer sec
integer seconds, days
!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! halo            The number of halo rows along all (NEWS) boundaries.
!                 There is currently no namelist option that allows an
!                 unequal halo boundary.  [integer, default, halo = 1]
!
!20150820  
! Removed 
!   do_netcdf_restart from namelist and remove non-netcdf restart output.
!   num_adjust_dt from namelist as it was not used.
!   num_dynam_dt from namelist as it was not used.
!   nt_zero from namelist as it was not used.
!   do_conserve_energy from namelist as it was not used.
!   decomp from namelist as it was not used.
!
!  IMPORTANT NOTE:
!      In versions of the model with an option to generate
!      it's own initial condition the variables ntrace and ntprog
!      do not have to be specified.
!


! namelist scm_atmosphere_nml
!
! Kdim       Number of vertical levels.
!
! As         A's for the vertical coordinate.
!
! Bs         B's for the vertical coordinate.
!
! pdamp      The pressure level above which damping is applied.

   integer, parameter :: maxpes = 32

   integer   :: halo = 1
! physics_window is replaced by nxblocks, nyblocks in atmos_model_nml
!   integer, dimension(2) :: physics_window = (/0,1/)

   logical   :: initialize_noncloud = .false.

   namelist /atmosphere_nml/ halo,        &
                             initialize_noncloud

   integer   :: nx=1, ny=1, Kdim=0
   real      :: pdamp=5000. 
   real, dimension(200):: As, Bs
   real, allocatable, dimension(:,:)   :: ztend  
   real, allocatable, dimension(:,:,:) :: phalf, pfull, zhalf, zfull
   real, allocatable, dimension(:,:,:) :: damp_t, damp_q, fill_q
real :: ps




!----------------namelist variables---------------------------------
   namelist /scm_atmosphere_nml/ Kdim, As, Bs, pdamp


!-----------------------------------------------------------------------
!---- private data ----


real                :: dt_atmos
! public              :: dt_atmos !yzheng
type    (time_type) :: Time_step_atmos

real :: tph0d_in = 0.0, tlm0d_in = 0.0

real, parameter :: rog  = rdgas/grav
real, parameter :: d608 = (rvgas-rdgas)/rdgas

integer           :: irestart_format = 1

integer           :: ntrace, ntprog, ntdiag

!-----------------------------------------------------------------------
!---- private data ----

  integer, dimension(4)              :: atmos_axes
  integer :: id_dynam, id_phys_down, id_phys_up, id_fv_diag
  type    (prog_var_type), save        :: Forc_tend
  type    (prog_var_type), save        :: Phys_tend

!-----------------------------------------------------------------------
  real, allocatable :: pref(:,:), dum1d(:)
  logical :: do_atmos_nudge


! &fv_core_nml
!        layout = 1,1
!        nlon=1
!        mlat=1
!        nlev= 24
!        target_lon  = 262.55
!        target_lat  = 36.65
!/


!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff, Grid_box)
#include "fv_arrays.h"

 type (time_type),     intent(in)    :: Time_init, Time, Time_step
 type(surf_diff_type), intent(inout) :: Surf_diff
 type(grid_box_type),  intent(inout) :: Grid_box

  integer :: unit
  integer :: ix, jx, kx, nt, ntp, ierr, io, n,i,j
  integer :: nql, nqi, nqa, sphum
  integer :: nqn, nqni, n_q2clo, n_prec
  integer :: logunit
  character(len=256) :: grid_file='INPUT/grid_spec.nc'

!-----------------------------------------------------------------------

  real    :: DEG2RAD
  integer :: ss, ds, ks
  integer :: dims(4)
  real,    allocatable, dimension(:)       :: xba,yba,xta,yta
#include "fv_point.inc"

!-----------------------------------------------------------------------
    
!----- read namelist -----
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=atmosphere_nml, iostat=io)
   ierr = check_nml_error(io, 'atmosphere_nml')
#else
   if (file_exist('input.nml')) then
        unit=open_namelist_file ();
        ierr=1
        do while (ierr /= 0)
           read (unit, nml=atmosphere_nml, iostat=io, end=5)
           ierr=check_nml_error (io, 'atmosphere_nml')
        enddo
 5      call close_file (unit);
   endif
#endif

    unit = stdlog()
    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number (version, tagname)
         write (unit, nml=atmosphere_nml)
    endif
    call close_file (unit)

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=scm_atmosphere_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_atmosphere_nml')
#else
    if (file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1 
        do while (ierr /= 0)
           read (unit, nml=scm_atmosphere_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'scm_atmosphere_nml')
        enddo
 10     call close_file (unit)
    endif
#endif

    unit = stdlog()
    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number (version, tagname)
         write (unit, nml=scm_atmosphere_nml)
    endif
    call close_file (unit)



   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)


!   get tracer indices

    sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    if (sphum <= 0) call error_mesg ('atmosphere_mod', &
         'specific humidity tracer not found', FATAL)
  
    nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
    if (nql <= 0) call error_mesg ('atmosphere_mod', &
         'liquid water tracer not found', FATAL)
    nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
    if (nqi <= 0) call error_mesg ('atmosphere_mod', &
         'ice water tracer not found', FATAL)
    nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
    if (nqa <= 0) call error_mesg ('atmosphere_mod', &
         'cloud amout tracer not found', FATAL)
    if (min(nql,nqi,nqa) <= 0) call error_mesg ('atmosphere_mod', &
         'stratiform cloud tracer(s) not found', FATAL)
    if (nql == nqi .or. nqa == nqi .or. nql == nqa) &
         call error_mesg ('atmosphere_mod',  &
         'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)
  
    nqn    = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
    nqni   = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

    ! Diagnostic tracer to allow forcing data be passed to 
    ! convection scheme for closure calculation.
    n_q2clo = get_tracer_index ( MODEL_ATMOS, 'a1q2_closure' )
    n_prec  = get_tracer_index ( MODEL_ATMOS, 'a1prec_closure' )

!===========================================================
!  Initialize tracer

! how many tracers have been registered?
    call get_number_tracers ( MODEL_ATMOS, num_tracers=ntrace, &
         num_prog=ntprog, num_diag=ntdiag )

!  ----- write tracer info to log file -----

    if (mpp_pe() == mpp_root_pe()) then
       logunit = stdlog()
       write (logunit, '(a,i3)') 'Number of tracers =', ntrace
       write (logunit, '(a,i3)') 'Number of prognostic tracers =', ntprog
       write (logunit, '(a,i3)') 'Number of diagnostic tracers =', ntdiag
    endif

!----- initialize FV dynamical core -----
!Read the AREA_ATM variable in the grid_spec file to determine the grid size.
! Should be 1x1.
    nx = 1 ; ny = 1
    if(field_exist(grid_file, 'AREA_ATM')) then
      call field_size(grid_file, "AREA_ATM", dims)
      nx = dims(1)
      ny = dims(2)
    else
      call error_mesg ('atmosphere_mod','cannot find AREA_ATM on INPUT/grid_spec.nc', FATAL)
    endif

    allocate ( xba(1:nx+1), yba(1:ny+1))

! Read the boundaries of the gridbox. 
! These will determine the resolution of the SCM
    call read_data(grid_file, "xba",  xba)
    call read_data(grid_file, "yba",  yba)

    SCM_RES_X = abs(xba(nx+1)-xba(nx))
    SCM_RES_Y = abs(yba(ny+1)-yba(ny))

! Read the cell center of the gridbox. 
! These will determine the location of the SCM experiment.
    allocate(xta(nx), yta(ny))

    call read_data(grid_file, "xta",  xta)
    call read_data(grid_file, "yta",  yta)
    target_lon = xta(1)
    target_lat = yta(1)

    !set up model boundaries
    !make sure requested longitude lies in the 0. to 360. range
    do while (target_lon .lt. 0.) 
      target_lon = target_lon + 360.
    enddo
    do while (target_lon .gt. 360.)
      target_lon = target_lon - 360.
    enddo

!----- make sure that latitude is within 90S - 90N -------
! Check your grid_spec file if this error is invoked.
    if (target_lat < -90. .or. target_lat > 90.) call error_mesg ('SCM_atmosphere',  &
          'requested site latitude outside of bounds', FATAL)

 
    call fv_init( sec )  ! allocates Atm components
    call prog_var_init (ntprog, Forc_tend) !Allocate the Forcing Tendency type
    call prog_var_init (ntprog, Phys_tend) !Allocate the Physics Tendency type
    allocate ( damp_t(beglon:endlon,beglat:endlat,nlev),          &
               damp_q(beglon:endlon,beglat:endlat,nlev),          &
               fill_q(beglon:endlon,beglat:endlat,nlev))


    call set_eta(nlev, ks, ptop, ak, bk)
! Read the observational data.
    call scm_data_read(nlev)

    DEG2RAD = PI/180.
    if (get_topog_mean(xba*DEG2RAD, yba*DEG2RAD, phis)) then
    else
       call error_mesg('atmosphere_init', 'error reading topography data', FATAL)
    endif

!Interpolate the first observed data onto the model grid.
    call scm_core_driver_init(Time, pdamp, phis, ak, bk) 

    !phis is height (m) when coming from grid_spec or initialization routine.
    !Convert to geopotential height.
    phis = phis * grav
    if ( cold_start .or. change_time ) then
        fv_time = time
    else
        fv_time = set_time (seconds, days)
        call get_time (Time, ss,  ds)

        if( seconds /= ss .or. days /= ds )   call  error_mesg         &
            ('FV_init:','Time inconsistent between fv_rst and INPUT/atmos_model.res', FATAL)
    endif

!----- initialize atmos_axes and fv_dynamics diagnostics

    call fv_diag_init( atmos_axes, Time )
    call scm_diagnostics_init ( Time, phis, atmos_axes )

    call scm_forc_diagnostic_init(atmos_axes, Time)

!    call fv_print_chksums( 'after fv_diag_init' )
   

!----- initialize physics interface -----
!----- initialize domains for reading global physics data -----

    call set_domain ( fv_domain )

!  --- initialize clocks for dynamics, physics_down and physics_up

    id_dynam     = mpp_clock_id ('FV dynamical core',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    id_fv_diag   = mpp_clock_id ('FV Diag',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

!--- allocate pref
    allocate(pref(nlev+1,2), dum1d(nlev+1))
!---------- reference profile -----------
    pref(nlev+1,1) = 101325.
    pref(nlev+1,2) = 81060.

    call get_eta_level ( nlev, pref(nlev+1,1), pref(1,1), dum1d )
    call get_eta_level ( nlev, pref(nlev+1,2), pref(1,2), dum1d )

!--- setup Grid_box area for physics
    allocate(Grid_box%area  (beglon:endlon, beglat:endlat))
    Grid_box%area  (beglon:endlon, beglat:endlat) = area (beglon:endlon, beglat:endlat)

!---- initialize cmip diagnostic output ----

    call atmos_cmip_diag_init   ( ak, bk, pref(1,1), atmos_axes, Time )
    call atmos_global_diag_init ( atmos_axes, Grid_box%area )

!--- initialize nudging module ---
    call atmos_nudge_init ( Time, atmos_axes(1:3), flag=do_atmos_nudge )

    call fv_print_chksums( 'Exiting  atmosphere_init' )

    
    ix= nlon ; jx = mlat ; kx = nlev

!-----------------------------------------------------------------------

 end subroutine atmosphere_init

!#######################################################################

 subroutine atmosphere_dynamics (Time, Surf_diff)
!
!        Time = time at the current time level
!
#include "fv_arrays.h"

   type(time_type),intent(in)    :: Time
   type(surf_diff_type), intent(in) :: Surf_diff

!---- dynamics -----


     real, dimension(size(pe,1),size(pe,3),size(pe,2))     :: p_full, z_full, delz
     real, dimension(size(pe,1),size(pe,3),size(pe,2) + 1) :: p_half, z_half
     integer :: sphum
#include "fv_point.inc"

    sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    if (sphum <= 0) call error_mesg ('atmosphere_mod', &
         'specific humidity tracer not found', FATAL)
  

!     call fv_compute_p_z (nlev, phis, pe, peln, delp, delz, &
!                          pt, q(:,:,:,sphum), &
!                          p_full, p_half, z_full, z_half, .true.)
   call scm_core_driver (Time, Time_step_atmos, pdamp, phis, Forc_tend)

 end subroutine atmosphere_dynamics

!#######################################################################

 subroutine atmosphere_end (Time, Grid_box)
 type (time_type),       intent(in) :: Time
 type(grid_box_type), intent(inout) :: Grid_box

 integer :: unit, pe
 character(len=64) :: fname_nc='RESTART/scm_atmosphere.res.nc'
 integer :: sphum, nql, nqi, nqa

   sphum = get_tracer_index ( MODEL_ATMOS, 'sphum'   )
   nql   = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
   nqi   = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
   nqa   = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )

!----- initialize domains for writing global physics data -----

    call set_domain ( fv_domain )
    call get_time (Time, seconds,  days)
    call write_fv_rst( 'RESTART/fv_rst.res', days, seconds, grav, &
         restart_format )

    call atmos_nudge_end

    call atmos_global_diag_end
    call fv_end(days, seconds)      

    deallocate (pref, dum1d)

    call scm_core_driver_end ()


!----- initialize domains for writing global physics data -----

!    call scm_physics_end   (Time)                           

 end subroutine atmosphere_end

  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  dummy routine.
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call error_mesg ('atmosphere_restart in atmosphere_mod', &
                     'intermediate restart capability is not implemented for this model', FATAL)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>

!#######################################################################
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_resolution (nlon, nlat, global)

  integer, intent(out)          :: nlon, nlat
  logical, intent(in), optional :: global

!---- return the size of the grid used for physics computations ----

   logical :: local

   local = .TRUE.
   if( PRESENT(global) )local = .NOT.global

   if( local )then
       nlon = endlon - beglon + 1
       nlat = endlat - beglat + 1
   else
       nlon = nlon
       nlat = mlat
   end if

 end subroutine atmosphere_resolution

 
  subroutine atmosphere_pref (p_ref)
    real, dimension(:,:), intent(inout) :: p_ref
 
    real,                 dimension(2,1)     :: pssl
    real,    allocatable, dimension(:,:,:)   :: dummy_pfull

    integer :: i,k
    
    p_ref = pref

  end subroutine atmosphere_pref
 
  subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro, do_uni_zfull)
    integer, intent(out)           :: i1, i2, j1, j2, kt
    logical, intent(out), optional :: p_hydro, hydro, do_uni_zfull

    i1 = beglon
    i2 = endlon
    j1 = beglat
    j2 = endlat
    kt = nlev
 !---non-hydrostatic is not supported by the SCM
    if (present(p_hydro)) p_hydro = .true.
    if (present(  hydro))   hydro = .true.
 
  end subroutine atmosphere_control_data
 
  subroutine atmosphere_cell_area  (area_out)
     real, dimension(:,:),  intent(out) :: area_out
 
    area_out(1:size(area_out,1), 1:size(area_out,2)) =  &
                                   area (beglon:endlon, beglat:endlat)
 
  end subroutine atmosphere_cell_area
 
 
 
 !---this is a dummy routine needed for compatibility with the 
 !---decoupling of physics and radiation from the dynamic cores
  subroutine atmosphere_grid_center (lon, lat)
 !---------------------------------------------------------------
 !    returns the longitude and latitude cell centers
 !---------------------------------------------------------------
     real,    intent(out) :: lon(:,:), lat(:,:)   ! Unit: radian
 ! Local data:
     real :: lon1(size(lon,1)), lat1(size(lat,2))
     integer i,j
 
    do j = beglat,endlat
      do i = beglon,endlon
        lon(i-beglon+1,j-beglat+1) = rlon(i,j)
        lat(i-beglon+1,j-beglat+1) = rlat(i,j)
      end do
    end do

  end subroutine atmosphere_grid_center
 
!#######################################################################
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_boundary (blon, blat, global)

    real,    intent(out)          :: blon(:,:), blat(:,:)
    logical, intent(in), optional :: global

! Local:
    integer i,j
    logical :: local

    local = .TRUE.
    if( PRESENT(global) )local = .NOT.global

    if( local )then
        do i = beglon,endlon+1
           blon(i-beglon+1,:) = rlonb(i)
        end do
        do j = beglat,endlat+1
           blat(:,j-beglat+1) = rlatb(j)
        end do
    else
        do i=1,nlon+1
           blon(i,:) = rlonb(i)
        end do
        do j=1,mlat+1
           blat(:,j) = rlatb(j)
        end do
    end if

 end subroutine atmosphere_boundary

!#######################################################################

 subroutine set_atmosphere_pelist ()
!--- no-op for SCM
 end subroutine set_atmosphere_pelist


!#######################################################################
!    returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

 subroutine atmosphere_domain (Domain)
 type(domain2d), intent(inout) :: Domain

   Domain = fv_domain

 end subroutine atmosphere_domain

!#######################################################################
!    returns the axis indices associated with the coupling grid

 subroutine get_atmosphere_axes ( axes )

   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----

     if ( size(axes) < 0 .or. size(axes) > 4 ) call error_mesg (    &
                           'get_atmosphere_axes in atmosphere_mod', &
                           'size of argument is incorrect', FATAL   )

     axes (1:size(axes(:))) = atmos_axes (1:size(axes(:)))

 end subroutine get_atmosphere_axes

!#######################################################################
! returns temp, sphum, pres, height at the lowest model level
!         and surface pressure

 subroutine get_bottom_mass (t_bot, tr_bot, p_bot, z_bot, p_surf, slp)

#include "fv_arrays.h"
! returns temp, sphum, pres, height at the lowest model level
!         and surface pressure and sea level pressure

   real, intent(out), dimension(beglon:endlon,beglat:endlat)  &
        :: t_bot, p_bot, z_bot, p_surf, slp
   real, intent(out), dimension(beglon:endlon,beglat:endlat,ncnst-pnats):: tr_bot
   integer :: i, j, k, kr
   real zvir, rrg, sigtop, sigbot
   real, dimension(beglon:endlon,beglat:endlat) :: tref
   real, parameter :: tlaps = 6.5e-3
#include "fv_point.inc"

   rrg  = rdgas / grav
   zvir = rvgas/rdgas - 1.


   ! determine 0.8 sigma reference level
   sigtop = ak(1)/pstd_mks+bk(1)
   do k = 1, nlev 
      sigbot = ak(k+1)/pstd_mks+bk(k+1)
      if (sigbot+sigtop > 1.6) then
         kr = k
         exit
      endif   
      sigtop = sigbot
   enddo

!$omp parallel do default(shared)    &
!$omp private (i, j)
     do j = beglat, endlat
        do i = beglon, endlon
           p_surf(i,j) =  ps(i,j)
           t_bot(i,j) =  pt(i,j,nlev)

           p_bot(i,j) = delp(i,j,nlev)/(peln(i,nlev+1,j)-peln(i,nlev,j))
           z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*q(i,j,nlev,1))*  &
                  (1. - pe(i,nlev,j)/p_bot(i,j))
           ! sea level pressure
           tref(i,j) = pt(i,j,kr)*(delp(i,j,kr)/((peln(i,kr+1,j)-peln(i,kr,j))*ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = ps(i,j)*(1.+tlaps*phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
! Copy tracers
!$omp parallel do default(shared)    &
!$omp private (i, j, k)
     do k = 1,ncnst-pnats
        do j = beglat,endlat
           do i = beglon,endlon
              tr_bot(i,j,k) = q(i,j,nlev,k)
           enddo
        enddo
     enddo


 end subroutine get_bottom_mass

!#######################################################################
! returns u and v on the mass grid at the lowest model level

 subroutine get_bottom_wind (u_bot, v_bot)

#include "fv_arrays.h"
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------

   real, intent(out), dimension(beglon:,beglat:) :: u_bot, v_bot

   integer i, j
#include "fv_point.inc"
!Balaji: this cannot work unless beglon=1: corrected declaration Lbounds
   do j=beglat,endlat
      do i=beglon,endlon
         u_bot(i,j) = u_srf(i,j)
         v_bot(i,j) = v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind

!#######################################################################

 subroutine get_bottom_data ( a, b, a_bot, b_bot, k_bot )

  real   , intent(in) , dimension(:,:,:) :: a    , b
  real   , intent(out), dimension(:,:)   :: a_bot, b_bot
  integer, intent(in) , dimension(:,:), optional :: k_bot

  integer :: i, j, kb, kb_min
    
                      kb_min = size(a,3)
  if (present(k_bot)) kb_min = minval (k_bot)

  if ( kb_min == size(a,3) ) then
          a_bot = a(:,:,kb_min)
          b_bot = b(:,:,kb_min)
  else
       do j = 1, size(a,2)
       do i = 1, size(a,1)
          kb = k_bot(i,j)
          a_bot(i,j) = a(i,j,kb)
          b_bot(i,j) = b(i,j,kb)
       enddo
       enddo
  endif

 end subroutine get_bottom_data

!#######################################################################

 subroutine put_bottom_data ( a_bot, b_bot, a, b, k_bot )

  real   , intent(in)   , dimension(:,:)   :: a_bot, b_bot
  real   , intent(inout), dimension(:,:,:) :: a    , b
  integer, intent(in)   , dimension(:,:), optional :: k_bot

  integer :: i, j, kb, kb_min
    
                      kb_min = size(a,3)
  if (present(k_bot)) kb_min = minval (k_bot)

  if ( kb_min == size(a,3) ) then
          a(:,:,kb_min) = a_bot
          b(:,:,kb_min) = b_bot
  else
       do j = 1, size(a,2)
       do i = 1, size(a,1)
          kb = k_bot(i,j)
          a(i,j,kb) = a_bot(i,j)
          b(i,j,kb) = b_bot(i,j)
       enddo
       enddo
  endif

 end subroutine put_bottom_data

!#######################################################################

! subroutine atmosphere_cell_area  (area_out)
!
!   real, dimension(:,:),  intent(out)          :: area_out       
!
!   area_out = 1.0                      
!
! end subroutine atmosphere_cell_area 

!#######################################################################

 subroutine get_stock_pe(index, value)

 ! This is a dummy routine.
 ! It is neccessary to satisfy revision 13.0.4.3.2.1 of atmos_coupled/atmos_model.f90
 ! Since that revision of atmos_coupled/atmos_model.f90 does nothing with the result,
 ! this routine can be a dummy.
 ! If and when the result is needed, it should be the total water content of the
 ! global atmosphere (Kg), including vapor, liquid and ice.

 integer, intent(in) :: index
 real, intent(out)   :: value

 value = 0.0

 end subroutine get_stock_pe

!#######################################################################

 subroutine atmosphere_state_update (Time, Physics_tendency, Physics, Atm_block)
#include "fv_arrays.h"

    type(time_type),              intent(in) :: Time
    type (physics_tendency_type), intent(in) :: Physics_tendency
    type (physics_type),          intent(in) :: Physics
    type (block_control_type),    intent(in) :: Atm_block
    type(time_type) :: Time_next, Time_step
 !--- local variables ---
    integer:: ibs, ibe, jbs, jbe, nb, sphum
    real:: zvir
#include "fv_point.inc"

    Time_step = set_time(int(dt_atmos), 0)
    Time_next = Time + Time_step
    zvir = rvgas/rdgas - 1.
 
 !--- put u/v tendencies into haloed arrays u_dt and v_dt
 !$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
    do nb = 1,Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)
 

     u_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%u_dt
     v_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%v_dt
     t_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%t_dt
     q_dt(ibs:ibe,jbs:jbe,:,1:ntprog) = Physics_tendency%block(nb)%q_dt

!--- diagnostic tracers are being updated in-place
!--- tracer fields must be returned to the Atm structure
     q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst) = Physics_tendency%block(nb)%qdiag

    enddo

    Phys_tend%T = t_dt - Forc_tend%t
    Phys_tend%u = u_dt - Forc_tend%u
    Phys_tend%v = v_dt - Forc_tend%v
    Phys_tend%r = q_dt(:,:,:,1:nt_prog) - Forc_tend%r(:,:,:,1:nt_prog)

   !---get tracer index for water vapor
   sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
   if (sphum <= 0) call error_mesg ('atmosphere_mod', &
        'specific humidity tracer not found', FATAL)
    
! If damping required then pfull below needs to be allocated and filled in.
    damp_t = 0.
    damp_q = 0.
   !---compute filler for water vapor
   where ( (q(:,:,:,sphum) + q_dt(:,:,:,sphum)*dt_atmos) &
        .lt. 0.5e-15 )             
        fill_q = (0.5e-15 - q(:,:,:,sphum) -        &
             q_dt(:,:,:,sphum)*dt_atmos)/dt_atmos           
   elsewhere
        fill_q = 0.        
   ENDWHERE
         
   q_dt(:,:,:,sphum) = q_dt(:,:,:,sphum) + fill_q

!---- diagnostics for FV dynamics -----
    call timing_on('FV_DIAG')
 
    call fv_diag(Time_next, nlon, mlat, nlev, beglat, endlat, &
         ntrace, zvir, dt_atmos, .false.)
    call timing_off('FV_DIAG')


    call scm_diagnostics ( Physics%block(1)%p_half,  &
                           Forc_tend, Phys_tend, Physics%block(1)%p_full, &
                           Physics%block(1)%z_full, Physics%block(1)%z_half,     &
                           damp_t, damp_q, fill_q, Time_next)

! Emulate the old style updating of the prognostic variables.
! Do this to allow for semi-prognostic ability.
!------ time differencing using physics tendencies -------
   
     call prog_var_time_diff (dt_atmos)

 
  end subroutine atmosphere_state_update

  !yzheng
  subroutine atmosphere_state_update_five (Time, Physics_tendency, Physics, Atm_block)
   type(time_type),              intent(in) :: Time
   type (physics_tendency_type), intent(in) :: Physics_tendency
   type (physics_type),          intent(in) :: Physics
   type (block_control_type),    intent(in) :: Atm_block
   type(time_type) :: Time_next, Time_step
!--- local variables ---
   integer:: ibs, ibe, jbs, jbe, nb, nt_tot, nt_prog
   real:: zvir
#include "fv_point.inc"
   
   Time_step = set_time(int(dt_atmos), 0)
   Time_next = Time + Time_step
   zvir = rvgas/rdgas - 1.
   
   call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)  
!--- put u/v tendencies into haloed arrays u_dt and v_dt
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)


    u_dt_five(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%u_dt
    v_dt_five(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%v_dt
    t_dt_five(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%t_dt
    q_dt_five(ibs:ibe,jbs:jbe,:,1:nt_prog) = Physics_tendency%block(nb)%q_dt

!--- diagnostic tracers are being updated in-place
!--- tracer fields must be returned to the Atm structure
    q_five(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst) = Physics_tendency%block(nb)%qdiag

   enddo
 
   call prog_var_time_diff_five (dt_atmos)
 
 end subroutine atmosphere_state_update_five

 subroutine prog_var_time_diff_five (dt, nt)
      real,                intent(in)    :: dt
      integer, optional,   intent(in)    :: nt
   
      integer :: ntp, n, nt_tot, nt_prog
   #include "fv_point.inc"
     call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)  

      ntp = nt_prog
      if (present(nt)) ntp = min(nt_prog, nt)
   
     ua_five = ua_five + dt * u_dt_five
     va_five = va_five + dt * v_dt_five
     pt_five = pt_five + dt * t_dt_five
     q_five(:,:,:,1:ntp) = q_five(:,:,:,1:ntp) + &
                               dt * q_dt_five(:,:,:,1:ntp)
   
   !----- zero out tendencies -----
      u_dt_five = 0.0
      v_dt_five = 0.0
      t_dt_five = 0.0
      q_dt_five(:,:,:,1:ntp) = 0.0
   
 end subroutine prog_var_time_diff_five


!#######################################################################

 subroutine atmos_physics_driver_inputs (Physics, Atm_block, Physics_tendency)
#include "fv_arrays.h"
    type (physics_type),  intent(inout) :: Physics
    type (block_control_type), intent(in) :: Atm_block
    type (physics_tendency_type), intent(inout), optional :: Physics_tendency
 !--- local variabls
    integer :: nb, ibs, ibe, jbs, jbe
    real,                 dimension(2,1)     :: pssl
#include "fv_point.inc"
 

 !---------------------------------------------------------------------
 ! use most up to date atmospheric properties when running serially
 !---------------------------------------------------------------------
 
 !$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)
 
     Physics%block(nb)%phis = phis(ibs:ibe,jbs:jbe)
     Physics%block(nb)%u    = ua(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%v    = va(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%t    = pt(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%q    = q(ibs:ibe,jbs:jbe,:,1:nt_prog)
     Physics%block(nb)%omega= omga(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%pe   = pe(ibs:ibe,:,jbs:jbe)
     Physics%block(nb)%peln = peln(ibs:ibe,:,jbs:jbe)
     Physics%block(nb)%delp = delp(ibs:ibe,jbs:jbe,:)
     if (.not.Physics%control%phys_hydrostatic) &
        call mpp_error(FATAL,'atmosphere: the non-hydrostatic option is not supported (phys_hydrostatic=.false.)')
     if (_ALLOCATED(Physics%block(nb)%tmp_4d)) &
        Physics%block(nb)%tmp_4d = q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst)

     call fv_compute_p_z (Atm_block%npz, Physics%block(nb)%phis, Physics%block(nb)%pe, &
                          Physics%block(nb)%peln, Physics%block(nb)%delp, Physics%block(nb)%delz, &
                          Physics%block(nb)%t, Physics%block(nb)%q(:,:,:,Physics%control%sphum), &
                          Physics%block(nb)%p_full, Physics%block(nb)%p_half, &
                          Physics%block(nb)%z_full, Physics%block(nb)%z_half, &
                          Physics%control%phys_hydrostatic)

     if (PRESENT(Physics_tendency)) then
!--- copy the dynamics tendencies into the physics tendencies
!--- if one wants to run physics concurrent with dynamics,
!--- these values would be zeroed out and accumulated
!--- in the atmosphere_state_update

       Physics_tendency%block(nb)%u_dt = u_dt(ibs:ibe,jbs:jbe,:)
       Physics_tendency%block(nb)%v_dt = v_dt(ibs:ibe,jbs:jbe,:)
       Physics_tendency%block(nb)%t_dt = t_dt(ibs:ibe,jbs:jbe,:)
       Physics_tendency%block(nb)%q_dt = q_dt(ibs:ibe,jbs:jbe,:,1:ntprog)
       Physics_tendency%block(nb)%qdiag = q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst)
     endif

    enddo
 
  end subroutine atmos_physics_driver_inputs
 
!#######################################################################
 
  subroutine atmos_radiation_driver_inputs (Time, Radiation, Atm_block )
#include "fv_arrays.h"
    type (time_type),      intent(in)    :: Time
    type (radiation_type), intent(inout) :: Radiation
    type (block_control_type), intent(in) :: Atm_block
 !--- local variables
    integer :: nb, ibs, ibe, jbs, jbe
   real,                 dimension(2,1)     :: pssl
#include "fv_point.inc"

 !---------------------------------------------------------------------
 ! use most up to date atmospheric properties when running serially
 !---------------------------------------------------------------------
 !$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
    do nb = 1,Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)
 
     Radiation%block(nb)%phis = phis(ibs:ibe,jbs:jbe)
     Radiation%block(nb)%t    = pt(ibs:ibe,jbs:jbe,:)
     Radiation%block(nb)%q    = q(ibs:ibe,jbs:jbe,:,1:nt_prog)
     Radiation%block(nb)%pe   = pe(ibs:ibe,:,jbs:jbe)
     Radiation%block(nb)%peln = peln(ibs:ibe,:,jbs:jbe)
     Radiation%block(nb)%delp = delp(ibs:ibe,jbs:jbe,:)
     if (.not.Radiation%control%phys_hydrostatic) &
        call mpp_error(FATAL,'atmosphere: the non-hydrostatic option is not supported (phys_hydrostatic=.false.)')

     call fv_compute_p_z (Atm_block%npz, Radiation%block(nb)%phis, Radiation%block(nb)%pe, &
                          Radiation%block(nb)%peln, Radiation%block(nb)%delp, Radiation%block(nb)%delz, &
                          Radiation%block(nb)%t, Radiation%block(nb)%q(:,:,:,Radiation%control%sphum), &
                          Radiation%block(nb)%p_full, Radiation%block(nb)%p_half, &
                          Radiation%block(nb)%z_full, Radiation%block(nb)%z_half, &
                          Radiation%control%phys_hydrostatic)
    enddo
 
 !----------------------------------------------------------------------
 ! obtain pressure-weighted global mean co2 dry volume mixing ratio for
 ! use by radiation package.
 !----------------------------------------------------------------------
 ! compute_g_avg must be called here because it contains
 ! mpp_sums that cannot be called during the concurrent radiation
 ! phase due to the way in which MPI interacts with nested OpenMP
 !----------------------------------------------------------------------
 !--- esm2_bugs is a logical switch (default=.FALSE.) read in via the nml
 !--- in order to allow early esm2 simulations to reproduce a bug
 !--- inside of compute_g_avg
!     call compute_g_avg(Time, 'co2', Radiation, Atm_block, esm2_bugs)
 
  end subroutine atmos_radiation_driver_inputs

!#######################################################################

subroutine fv_compute_p_z (npz, phis, pe, peln, delp, delz, pt, q_sph, p_full, p_half, z_full, z_half, hydrostatic)
     integer, intent(in)  :: npz
     real, dimension(:,:),   intent(in)  :: phis
     real, dimension(:,:,:), intent(in)  :: pe, peln, delp, delz, pt, q_sph
     real, dimension(:,:,:), intent(out) :: p_full, p_half, z_full, z_half
     logical, intent(in)  :: hydrostatic
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
 
 !--------- Hydrostatic option ----------------------------------------------
     if (hydrostatic ) then
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
     else
!--------- Non-Hydrostatic option ------------------------------------------
       do k=npz,1,-1
         do j=1,size(phis,2)
           do i=1,size(phis,1)
             p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
             z_half(i,j,k) = z_half(i,j,k+1) - delz(i,j,k)
             z_full(i,j,k) = 0.5*(z_half(i,j,k) + z_half(i,j,k+1))
           enddo
         enddo
       enddo
     endif
 
  end subroutine fv_compute_p_z

!#######################################################################

  subroutine reset_atmos_tracers (Physics, Physics_tendency, Atm_block)
#include "fv_arrays.h"
    type (physics_type), intent(in) :: Physics
    type (physics_tendency_type), intent(in) :: Physics_tendency
    type (block_control_type), intent(in) :: Atm_block
 !--- local variables
    integer :: nb, ibs, ibe, jbs, jbe
#include "fv_point.inc"
 
 !--- After initialization by the physics, tracer fields must be
 !--- returned to the Atm structure.  This is because tracer driver
 !--- can reset the initial values
 !$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
     do nb = 1, Atm_block%nblks
       ibs = Atm_block%ibs(nb)
       ibe = Atm_block%ibe(nb)
       jbs = Atm_block%jbs(nb)
       jbe = Atm_block%jbe(nb)
 
       q(ibs:ibe,jbs:jbe,:,1:nt_prog)       = Physics%block(nb)%q
       q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst) = Physics_tendency%block(nb)%qdiag
     enddo
 
  end subroutine reset_atmos_tracers
 
!#######################################################################

end module atmosphere_mod
