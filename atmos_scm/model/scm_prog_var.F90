
module scm_prog_var_mod

!-----------------------------------------------------------------------
!
!       allocates storage for the basic dynamical variables
!
!-----------------------------------------------------------------------
!--------- public defined data type prog_var_type ----------------------
!
!     nlon = number of longitude points (first dimension)
!            includes 2 halo points (1 west, 1 east)
!     nlat = number of latitude points (second dimension)
!            includes 3 halo points
!     nlev = number of vertical levels
!
!     ntrace = number of tracers
!
!     u    = zonal wind component
!     v    = meridional wind component
!     t    = temperature
!     r    = arbitrary number of tracers (includes specific humidity)
!
!     ps   = surface pressure
!     pssl = surface pressure adjust to eta=1. (for eta coordinate)
!
!-----------------------------------------------------------------------

use scm_horiz_grid_mod, only: horiz_grid_type
#ifdef INTERNAL_FILE_NML
use              mpp_mod, only: input_nml_file
#else
use              fms_mod, only: open_namelist_file
#endif

use              mpp_mod, only: stdlog  
use              fms_mod, only: file_exist,                 &
                                check_nml_error,            &
                                write_version_number,       &
                                mpp_pe, mpp_npes,           &
                                mpp_root_pe,                &
                                close_file
use    field_manager_mod, only: MODEL_ATMOS
use   tracer_manager_mod, only: get_tracer_index

use            fv_pack, only: nlon, mlat, &
                              pnats, consv_te, ptop, fv_init, fv_domain, &
                              fv_end, change_time, p_var, restart_format, area, &
                              ak, bk, rlon, rlat, ng_d, nt_prog
use    mpp_domains_mod, only: mpp_get_data_domain

implicit none
private

public  prog_var_type, prog_var_init, var_init,  &
        copy_prog_var, assignment (=), prog_var_time_diff,  &
        prog_var_times_scalar

type prog_var_type
     integer       :: nlon, nlat, nlev, ntrace
     integer       :: ilb, iub, jlb, jub, klb, kub
     real, pointer :: ps(:,:)=>NULL(), pssl(:,:)=>NULL()
     real, pointer :: u(:,:,:)=>NULL(), v(:,:,:)=>NULL(), t(:,:,:)=>NULL(), r(:,:,:,:)=>NULL()
end type prog_var_type

interface var_init
    module procedure var_init_type_4d, var_init_bound_4d, &
                     var_init_type_3d, var_init_bound_3d, &
                     var_init_type_2d, var_init_bound_2d
end interface

interface assignment (=)
    module procedure prog_var_equals_scalar
end interface


logical :: do_semi_prog = .false. ! Flag to allow tendency updates to be ignored.

!----------------namelist variables---------------------------------

namelist /scm_prog_var_nml/ do_semi_prog

!-----------------------------------------------------------------------

integer :: nsphum

contains

!#######################################################################

 subroutine prog_var_init (ntrs, Vars)
#include "fv_arrays.h"
    integer,               intent(in)  :: ntrs
    type(prog_var_type), intent(inout) :: Vars
!-----------------------------------------------------------------------
    integer :: n, io, ierr, unit
#include "fv_point.inc"

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=scm_prog_var_nml, iostat=io)
    ierr = check_nml_error(io, 'scm_prog_var_nml')
#else
    if (file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1 
        do while (ierr /= 0)
           read (unit, nml=scm_prog_var_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'scm_atmosphere_nml')
        enddo
 10     call close_file (unit)
    endif
#endif

!    unit = stdlog()
!    if ( mpp_pe() == mpp_root_pe() ) then
!         call write_version_number (version, tagname)
!         write (unit, nml=scm_prog_var_nml)
!    endif
!    call close_file (unit)

    call mpp_get_data_domain   ( fv_domain, Vars % ilb, Vars % iub, &
                                            Vars % jlb, Vars % jub )
    Vars % klb = 1
    Vars % kub = nlev

    Vars % nlon = nlon
    Vars % nlat = mlat
    Vars % nlev = nlev
    Vars % ntrace = ntrs
    
    Vars % ps   => var_init_bound_2d (Vars % ilb, Vars % iub, &
                                      Vars % jlb, Vars % jub)

    Vars % pssl => var_init_bound_2d (Vars % ilb, Vars % iub, &
                                      Vars % jlb, Vars % jub)

    Vars % u => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % v => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % t => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % r => var_init_bound_4d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev, ntrs)
    nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )

 end subroutine prog_var_init

!#######################################################################

 function var_init_bound_2d (ilb, iub, jlb, jub) result (var)

    integer, intent(in)           :: ilb, iub, jlb, jub
    real, dimension(:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub) )
    var = 0.0

 end function var_init_bound_2d

!#######################################################################

 function var_init_type_2d (Hgrid) result (var)

    type(horiz_grid_type), intent(in) :: Hgrid
    real, dimension(:,:), pointer     :: var

    var => var_init_bound_2d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub)

 end function var_init_type_2d

!#######################################################################

 function var_init_bound_3d (ilb, iub, jlb, jub, kdim) result (var)

    integer, intent(in)             :: ilb, iub, jlb, jub, kdim
    real, dimension(:,:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub, 1:kdim) )
    var = 0.0

 end function var_init_bound_3d

!#######################################################################

 function var_init_type_3d (Hgrid, kdim) result (var)

    type(horiz_grid_type), intent(in) :: Hgrid
    integer, intent(in)               :: kdim
    real, dimension(:,:,:), pointer   :: var

    var => var_init_bound_3d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub, kdim)

 end function var_init_type_3d

!#######################################################################

 function var_init_bound_4d (ilb, iub, jlb, jub, kdim, ntrace) &
                     result (var)

    integer, intent(in)               :: ilb, iub, jlb, jub, &
                                         kdim, ntrace
    real, dimension(:,:,:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub, 1:kdim, 1:ntrace) )
    var = 0.0

 end function var_init_bound_4d

!#######################################################################

 function var_init_type_4d (Hgrid, kdim, ntrace) result (var)

    type(horiz_grid_type), intent(in)   :: Hgrid
    integer, intent(in)                 :: kdim, ntrace
    real, dimension(:,:,:,:), pointer   :: var

    var => var_init_bound_4d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub, kdim, ntrace)

 end function var_init_type_4d

!#######################################################################

!subroutine prog_var_equals_prog_var (Left, Right)
 subroutine copy_prog_var (Left, Right)

    type(prog_var_type), intent(inout) :: Left
    type(prog_var_type), intent(in)    :: Right

        Left % u = Right % u
        Left % v = Right % v
        Left % t = Right % t
        Left % r = Right % r

        Left % ps   = Right % ps
        Left % pssl = Right % pssl

!!      Left % nlon = Right % nlon
!!      Left % nlat = Right % nlat
!!      Left % nlev = Right % nlev
!!      Left % ntrace = Right % ntrace
!!      Left % ntprog = Right % ntprog

!!      Left % ilb = Right % ilb
!!      Left % iub = Right % iub
!!      Left % jlb = Right % jlb
!!      Left % jub = Right % jub
!!      Left % klb = Right % klb
!!      Left % kub = Right % kub

 end subroutine copy_prog_var
!end subroutine prog_var_equals_prog_var

!#######################################################################

 subroutine prog_var_equals_scalar (Left, scalar)

    type(prog_var_type), intent(inout) :: Left
    real               , intent(in)    :: scalar

        Left % u = scalar
        Left % v = scalar
        Left % t = scalar
        Left % r = scalar

        Left % ps   = scalar
        Left % pssl = scalar

 end subroutine prog_var_equals_scalar

!#######################################################################

 subroutine prog_var_times_scalar (Var, scalar)

    type(prog_var_type), intent(inout) :: Var
    real               , intent(in)    :: scalar

        Var % u = Var % u * scalar
        Var % v = Var % v * scalar
        Var % t = Var % t * scalar
        Var % r = Var % r * scalar

        Var % ps   = Var % ps   * scalar
        Var % pssl = Var % pssl * scalar

 end subroutine prog_var_times_scalar

!#######################################################################

 subroutine prog_var_time_diff (dt, nt)
#include "fv_arrays.h"

   real,                intent(in)    :: dt
   integer, optional,   intent(in)    :: nt

   integer :: ntp, n
#include "fv_point.inc"

   ntp = nt_prog
   if (present(nt)) ntp = min(nt_prog, nt)

   if ( .not. do_semi_prog) then
!     Var % ps   = Var % ps   + dt * Var_dt % ps
!     Var % pssl = Var % pssl + dt * Var_dt % pssl

     ua = ua + dt * u_dt
     va = va + dt * v_dt
     pt = pt + dt * t_dt
   endif

   if ( .not. do_semi_prog) then
     q(:,:,:,1:ntp) = q(:,:,:,1:ntp) + &
                                dt * q_dt(:,:,:,1:ntp)
   else
!Update tracers other than specific humidity
     do n = 1, ntp
       if ( n /= nsphum ) &
         q(:,:,:,n) = q(:,:,:,n) + &
                                  dt * q_dt(:,:,:,n)
     enddo
   endif

! ZNT 04/10/2020: update Usfc and Vsfc
   u_srf(:,:) = ua(:,:,nlev)
   v_srf(:,:) = va(:,:,nlev)

!----- zero out tendencies -----

!   Var_dt % ps   = 0.0
!   Var_dt % pssl = 0.0

   u_dt = 0.0
   v_dt = 0.0
   t_dt = 0.0

   q_dt(:,:,:,1:ntp) = 0.0


 end subroutine prog_var_time_diff

!#######################################################################

end module scm_prog_var_mod

