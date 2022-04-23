module five_mod
#include <fms_platform.h>

use physics_types_mod,  only: physics_tendency_type, physics_type, &
                              physics_input_block_type, physics_tendency_block_type
use       constants_mod, only: rdgas, grav, rvgas, radius, PI

use       constants_mod, only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks                             
use block_control_mod,  only: block_control_type
use physics_radiation_exch_mod,only: radiation_flux_block_type, &
                                      radiation_flux_type
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index, get_number_tracers
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
!-----------------
! FV core modules:
!-----------------
use            fv_pack, only: ak, bk, nlon, mlat, nlev, ncnst, get_eta_level
 
!-----------------------------------------------------------------
implicit none
private

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

public five_init, atmos_physics_driver_inputs_five, &
                  five_tend_low_to_high, five_tend_high_to_low, &
                  five_var_high_to_low_4d, &
                  update_bomex_forc_five, &
                  atmosphere_pref_five

!------------------------------------------------------------------------
!---prognostic variables and their tendencies in FIVE grid---------------
!------------------------------------------------------------------------
real, allocatable:: ua_five(:,:,:)      ! (ua, va) are mostly used as the A grid winds
real, allocatable:: va_five(:,:,:)
real, allocatable :: pt_five(:,:,:)   ! temperature (K)
real, allocatable :: delp_five(:,:,:) ! pressure thickness (pascal)
real, allocatable :: q_five(:,:,:,:)  ! specific humidity and constituents

real, allocatable :: ps_five (:,:)      ! Surface pressure (pascal)
real, allocatable :: pe_five (:,:,: )   ! edge pressure (pascal)
real, allocatable :: pk_five  (:,:,:)   ! pe**cappa
real, allocatable :: peln_five(:,:,:)   ! ln(pe)
real, allocatable :: pkz_five (:,:,:)   ! finite-volume mean pk
real, allocatable :: phis_five(:,:)     ! Surface geopotential (g*Z_surf)
real, allocatable :: omga_five(:,:,:)   ! Vertical pressure velocity (pa/s)

real, allocatable, dimension(:,:,:) :: u_dt_five, v_dt_five, t_dt_five
real, allocatable :: q_dt_five(:,:,:,:)

public ua_five, va_five, pt_five, q_five
public u_dt_five, v_dt_five, t_dt_five, q_dt_five


!------------------------------------------------------------------------
!---namelist-------------------------------------------------------------
!------------------------------------------------------------------------
logical :: do_five = .FALSE.
public do_five

! This is the number of layers to add between native levels
!  NOTE: This must be an EVEN number, due to limitations
!  in the tendency interpolation scheme
integer :: five_add_nlevels = 2
public five_add_nlevels

! The bottom layer to which we will add layers to (set
!   to a very large value to add all the way to surface, though
!   currently setting this value to surface results in model
!   crashes in SCM, so need to investigate)
real, private :: five_bot_toadd = 100000.

! The top layer to which we will add layers to
real, private :: five_top_toadd = 50000.

NAMELIST / five_nml / do_five, five_add_nlevels, five_bot_toadd, five_top_toadd

!------------------------------------------------------------------------
!---Internal variables used for calculation------------------------------
!------------------------------------------------------------------------
integer :: five_bot_k, five_top_k ! indicees where levels are added
real, allocatable :: pf0_five(:) ! midpoint FIVE pressures (pascals), the 1D version (only used in the initialization)
real, allocatable :: ph0_five(:) ! interface FIVE pressures (pascals), the 1D version (only used in the initialization)
real, allocatable :: ak_five(:)
real, allocatable :: bk_five(:) 

!variables that are used for interpolation/averaging purpose.
real, allocatable :: p_half_five(:,:,:) ! interface FIVE pressures (pascals)
real, allocatable :: p_full_five(:,:,:) ! midpoint FIVE pressures (pascals)
real, allocatable :: z_full_five(:,:,:) ! midpoint FIVE pressures (pascals)
real, allocatable :: delz_five(:,:,:)
real, allocatable :: rho_five(:,:,:)

!variables from the host grid, also used for interpolation/averaging purpose
real, allocatable :: p_half_host(:,:,:)
real, allocatable :: p_full_host(:,:,:)
real, allocatable :: delz_host(:,:,:)
real, allocatable :: rho_host(:,:,:)

integer :: nlev_five !total levels
public :: nlev_five

logical :: nonzero_rad_flux_init = .false.

real, allocatable :: pref_five(:,:), dum1d(:)
real, parameter::ptop_min = 1.E-6  ! minimum pressure (pa) at model top to avoid
! floating point exception; this is not needed
! if model top is not at zero

!------------------------------------------------------------------------
!------ constants-------
!-----------------------------------------------------------------------
real, parameter :: ps0 = 1017.8e2
real, parameter :: zsfc = 0

!------------------------------------------------------------------------
!------ variables for loop-------
!-----------------------------------------------------------------------
integer :: i, j, k, p

contains

subroutine five_init(Physics_five, Physics_tendency_five, Rad_flux_five, &
                    Atm_block, lonb, latb, &
                    p_hydro, hydro, do_uni_zfull)
#include "fv_arrays.h"
      
    type (physics_type), intent(inout)          :: Physics_five
    type (physics_tendency_type), intent(inout) :: Physics_tendency_five
    type (radiation_flux_type), intent(inout) :: Rad_flux_five(:)
    type (block_control_type), intent(in) :: Atm_block
    real,    dimension(:,:),      intent(in)    :: lonb, latb
    logical,               intent(in) :: p_hydro, hydro, do_uni_zfull

    !--- local varialbes
    integer :: n, nb, ix, jx, npz, nt_tot, nt_prog
    integer          ::  id, jd
    integer :: unit, ierr, io
#include "fv_point.inc"

!----- read namelist -----
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=five_nml, iostat=io)
   ierr = check_nml_error(io, 'five_nml')
#else
   if (file_exist('input.nml')) then
        unit=open_namelist_file ();
        ierr=1
        do while (ierr /= 0)
           read (unit, nml=five_nml, iostat=io, end=5)
           ierr=check_nml_error (io, 'five_nml')
        enddo
 5      call close_file (unit);
   endif
#endif

    unit = stdlog()
    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number (version, tagname)
         write (unit, nml=five_nml)
    endif
    call close_file (unit)

    if (.not. do_five) return

    !compute nlev_five and ph
    !should replace ps0 with the ps(1,1)
    call five_pressure_init(pe(1,:,1), ak/(0.01*ps(1,1)),bk/0.01, ps(1,1), nlev_five)

    write (*,*) 'initial ps', ps

    !allocate five variables
    allocate ( ua_five(nlon, mlat, nlev_five) )       ; ua_five    = 0.0
    allocate ( va_five(nlon, mlat, nlev_five) )       ; va_five    = 0.0
    allocate (delp_five(nlon, mlat,   nlev_five))     ; delp_five  = 0.0
    allocate (  pt_five(nlon, mlat, nlev_five))       ; pt_five    = 0.0
    allocate (q_five(nlon, mlat, nlev_five, ncnst))   ; q_five     = 0.0

    allocate (phis_five(nlon, mlat))                  ; phis_five  = 0.0
    allocate (ps_five  (nlon, mlat))                  ; ps_five    = 0.0

    allocate (pkz_five (nlon, mlat, nlev_five))       ; pkz_five   = 0.0
    allocate (pk_five  (nlon, mlat, nlev_five+1))     ; pk_five    = 0.0
    allocate (pe_five  (nlon, nlev_five+1, mlat))     ; pe_five    = 0.0
    allocate (peln_five(nlon, nlev_five+1, mlat))     ; peln_five  = 0.0

    allocate (omga_five(nlon, mlat, nlev_five))       ; omga_five  = 0.0

    allocate( u_dt_five(nlon, mlat, nlev_five) )     ; u_dt_five   = 0.0
    allocate( v_dt_five(nlon, mlat, nlev_five) )     ; v_dt_five   = 0.0
    allocate( t_dt_five(nlon, mlat, nlev_five) )     ; t_dt_five   = 0.0
    allocate( q_dt_five(nlon, mlat, nlev_five, ncnst) ) ; q_dt_five   = 0.0

    allocate ( p_half_five(nlon, mlat, nlev_five + 1) )       ; p_half_five    = 0.0
    allocate ( p_full_five(nlon, mlat, nlev_five) )       ; p_full_five    = 0.0
    allocate ( z_full_five(nlon, mlat, nlev_five) )       ; z_full_five    = 0.0
    allocate ( delz_five(nlon, mlat, nlev_five) )       ; delz_five    = 0.0
    allocate ( rho_five(nlon, mlat, nlev_five) )       ; rho_five    = 0.0

    allocate ( p_half_host(nlon, mlat, nlev + 1) )       ; p_half_host    = 0.0
    allocate ( p_full_host(nlon, mlat, nlev) )       ; p_full_host    = 0.0
    allocate ( delz_host(nlon, mlat, nlev) )       ; delz_host    = 0.0
    allocate ( rho_host(nlon, mlat, nlev) )       ; rho_host    = 0.0

    !allocate Physics_five and Physics_tendency_five variables
    !---control data
    npz = nlev_five
    Physics_five%control%phys_hydrostatic = p_hydro
    Physics_five%control%hydrostatic = hydro
    Physics_five%control%do_uni_zfull = do_uni_zfull !miz

    call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)
    Physics_five%control%sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

    !---allocate global quantities
    !--- set pref
    allocate (Physics_five%glbl_qty%pref(npz+1,2))
    Physics_five%glbl_qty%pref = 0.

    !---allocate input block data
    allocate (Physics_five%block(Atm_block%nblks))
    do n = 1, Atm_block%nblks
        ix = Atm_block%ibe(n)-Atm_block%ibs(n)+1
        jx = Atm_block%jbe(n)-Atm_block%jbs(n)+1
        allocate (Physics_five%block(n)%phis   (ix,jx),             &
                Physics_five%block(n)%tmp_4d (ix,jx,npz,nt_prog+1:nt_tot), &
                Physics_five%block(n)%u      (ix,jx,npz),         &
                Physics_five%block(n)%v      (ix,jx,npz),         &
                Physics_five%block(n)%w      (ix,jx,npz),         &
                Physics_five%block(n)%t      (ix,jx,npz),         &
                Physics_five%block(n)%q      (ix,jx,npz,nt_prog), &
                Physics_five%block(n)%omega  (ix,jx,npz),         &
                Physics_five%block(n)%pe     (ix,npz+1,jx),       &
                Physics_five%block(n)%peln   (ix,npz+1,jx),       &
                Physics_five%block(n)%delp   (ix,jx,npz),         &
                Physics_five%block(n)%delz   (ix,jx,npz),         &
                Physics_five%block(n)%p_full (ix,jx,npz),         &
                Physics_five%block(n)%p_half (ix,jx,npz+1),       &
                Physics_five%block(n)%z_full (ix,jx,npz),         &
                Physics_five%block(n)%z_half (ix,jx,npz+1)        )

        Physics_five%block(n)%phis   = 0.
        Physics_five%block(n)%tmp_4d   = 0.
        Physics_five%block(n)%u      = 0.
        Physics_five%block(n)%v      = 0.
        Physics_five%block(n)%w      = 0.
        Physics_five%block(n)%t      = 0.
        Physics_five%block(n)%q      = 0.
        Physics_five%block(n)%omega  = 0.
        Physics_five%block(n)%pe     = 0.
        Physics_five%block(n)%peln   = 0.
        Physics_five%block(n)%delp   = 0.
        Physics_five%block(n)%delz   = 0.
        Physics_five%block(n)%p_full = 0.
        Physics_five%block(n)%p_half = 0.
        Physics_five%block(n)%z_full = 0.
        Physics_five%block(n)%z_half = 0.
    enddo

    allocate ( Physics_tendency_five%block(Atm_block%nblks) )
    do n = 1,Atm_block%nblks
      allocate ( Physics_tendency_five%block(n)%u_dt(ix,jx,npz), &
                 Physics_tendency_five%block(n)%v_dt(ix,jx,npz), &
                 Physics_tendency_five%block(n)%t_dt(ix,jx,npz), &
                 Physics_tendency_five%block(n)%q_dt(ix,jx,npz,nt_prog), &
                 Physics_tendency_five%block(n)%qdiag(ix,jx,npz,nt_prog+1:nt_tot) )

      Physics_tendency_five%block(n)%u_dt = 0.
      Physics_tendency_five%block(n)%v_dt = 0.
      Physics_tendency_five%block(n)%t_dt = 0.
      Physics_tendency_five%block(n)%q_dt = 0.
      Physics_tendency_five%block(n)%qdiag = 0.
    enddo

    !initialize Rad_flux_five
    do n = 1, size(Rad_flux_five,1)
        allocate (Rad_flux_five(n)%block(Atm_block%nblks))
        do nb = 1, Atm_block%nblks
          ix = Atm_block%ibe(nb) - Atm_block%ibs(nb) + 1
          jx = Atm_block%jbe(nb) - Atm_block%jbs(nb) + 1
          call Rad_flux_five(n)%block(nb)%alloc ( ix, jx, npz, nonzero_rad_flux_init )
        end do
    end do

    !Compute pfull for initial interpolation
    do k=size(pt,3),1,-1
        do j=1,size(pt,2)
          do i=1,size(pt,1)
            p_full_host(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
          enddo
        enddo
    enddo

    !initialize some parameters used to compute p_half, z_half, p_full, z_full
    ! --- Create delp
    ps_five = ps(1,1)

    do k=1,size(pt_five,3)
      do j=1,size(ps_five,2)
        do i=1,size(ps_five,1)
          delp_five(i,j,k) = ak_five(k+1)-ak_five(k) + ps_five(i,j)*(bk_five(k+1)-bk_five(k))
        enddo
      enddo
    enddo
    call p_var(nlon, mlat, nlev_five, 1, mlat, ak_five(1), delp_five, ps_five,     &
           pe_five, peln_five,  pk_five,  pkz_five,  2./7.)

    !Compute pfull for initial interpolation
    do k=size(pt_five,3),1,-1
        do j=1,size(pt_five,2)
            do i=1,size(pt_five,1)
              p_full_five(i,j,k) = delp_five(i,j,k)/(peln_five(i,k+1,j)-peln_five(i,k,j))
            enddo
        enddo
    enddo

    write (*,*) 'pe_five', pe_five
    write (*,*) 'p_full_five', p_full_five

    !initialize state variables: u, v, t, q
    call five_profiles_init(ua, va, pt, q, omga, p_full_host)

end subroutine five_init

subroutine atmos_physics_driver_inputs_five (Physics_five, Atm_block, Physics_tendency_five)

  type (physics_type),  intent(inout) :: Physics_five
  type (block_control_type), intent(in) :: Atm_block
  type (physics_tendency_type), intent(inout), optional :: Physics_tendency_five
  !--- local variabls
  integer :: nb, ibs, ibe, jbs, jbe, nt_tot, nt_prog

  call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)

   !---------------------------------------------------------------------
   ! use most up to date atmospheric properties when running serially
   !---------------------------------------------------------------------
  do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)

      Physics_five%block(nb)%phis = phis_five(ibs:ibe,jbs:jbe)
      Physics_five%block(nb)%u    = ua_five(ibs:ibe,jbs:jbe,:)
      Physics_five%block(nb)%v    = va_five(ibs:ibe,jbs:jbe,:)
      Physics_five%block(nb)%t    = pt_five(ibs:ibe,jbs:jbe,:)
      Physics_five%block(nb)%q    = q_five(ibs:ibe,jbs:jbe,:,1:nt_prog)
      Physics_five%block(nb)%omega= omga_five(ibs:ibe,jbs:jbe,:)
      Physics_five%block(nb)%pe   = pe_five(ibs:ibe,:,jbs:jbe)
      Physics_five%block(nb)%peln = peln_five(ibs:ibe,:,jbs:jbe)
      Physics_five%block(nb)%delp = delp_five(ibs:ibe,jbs:jbe,:)
      if (_ALLOCATED(Physics_five%block(nb)%tmp_4d)) &
          Physics_five%block(nb)%tmp_4d = q_five(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst)

      !   Compute heights
      call fv_compute_p_z (nlev_five, Physics_five%block(nb)%phis, Physics_five%block(nb)%pe, &
      Physics_five%block(nb)%peln, Physics_five%block(nb)%delp, Physics_five%block(nb)%delz, &
      Physics_five%block(nb)%t, Physics_five%block(nb)%q(:,:,:,Physics_five%control%sphum), &
      Physics_five%block(nb)%p_full, Physics_five%block(nb)%p_half, &
      Physics_five%block(nb)%z_full, Physics_five%block(nb)%z_half, &
      Physics_five%control%phys_hydrostatic)
      
      if (PRESENT(Physics_tendency_five)) then
          !--- copy the dynamics tendencies into the physics tendencies
          !--- if one wants to run physics concurrent with dynamics,
          !--- these values would be zeroed out and accumulated
          !--- in the atmosphere_state_update
          Physics_tendency_five%block(nb)%u_dt = u_dt_five(ibs:ibe,jbs:jbe,:)
          Physics_tendency_five%block(nb)%v_dt = v_dt_five(ibs:ibe,jbs:jbe,:)
          Physics_tendency_five%block(nb)%t_dt = t_dt_five(ibs:ibe,jbs:jbe,:)
          Physics_tendency_five%block(nb)%q_dt = q_dt_five(ibs:ibe,jbs:jbe,:,1:nt_prog)
          Physics_tendency_five%block(nb)%qdiag = q_five(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst)
      endif
  enddo

end subroutine atmos_physics_driver_inputs_five

subroutine five_tend_low_to_high (Physics_input_block, Physics_tendency_block, Rad_flux_block, &
    Physics_five_input_block, Physics_tendency_five_block, Rad_flux_five_block)
    type(physics_input_block_type), intent(in)      :: Physics_input_block
    type (physics_tendency_block_type), intent(in)  :: Physics_tendency_block
    type (radiation_flux_block_type), intent(in)  :: Rad_flux_block

    type(physics_input_block_type), intent(in)      :: Physics_five_input_block
    type (physics_tendency_block_type), intent(inout)  :: Physics_tendency_five_block
    type (radiation_flux_block_type), intent(inout)  :: Rad_flux_five_block

    ! Local variables	
    integer :: nt_tot, nt_prog
  
  !---------------------------------------------------------------------
  !    set up local pointers into the physics input and physics tendency
  !    blocks.
  !---------------------------------------------------------------------
    real, dimension(:,:,:), pointer :: z_full_host, z_half_host, delp_host
    real, dimension(:,:,:), pointer :: u_dt_host, v_dt_host, t_dt_host
    real, dimension(:,:,:,:), pointer :: r_dt_host, rdiag_host
  
    real, dimension(:,:,:), pointer :: z_half_five0, delp_five0
    real, dimension(:,:,:), pointer :: u_dt_five0, v_dt_five0, t_dt_five0
    real, dimension(:,:,:,:), pointer :: r_dt_five0, rdiag_five0

    z_full_host => Physics_input_block%z_full
    z_half_host => Physics_input_block%z_half
    delp_host => Physics_input_block%delp
    u_dt_host => Physics_tendency_block%u_dt
    v_dt_host => Physics_tendency_block%v_dt
    t_dt_host => Physics_tendency_block%t_dt
    r_dt_host => Physics_tendency_block%q_dt
    rdiag_host => Physics_tendency_block%qdiag

    z_half_five0 => Physics_five_input_block%z_half
    delp_five0 => Physics_five_input_block%delp
    u_dt_five0 => Physics_tendency_five_block%u_dt
    v_dt_five0 => Physics_tendency_five_block%v_dt
    t_dt_five0 => Physics_tendency_five_block%t_dt
    r_dt_five0 => Physics_tendency_five_block%q_dt
    rdiag_five0 => Physics_tendency_five_block%qdiag

    !pass values to  variables that will not only be used here but also be shared by other modules
    p_full_host = Physics_input_block%p_full
    p_half_host = Physics_input_block%p_half
    p_full_five = Physics_five_input_block%p_full
    p_half_five = Physics_five_input_block%p_half
    z_full_five = Physics_five_input_block%z_full

    do k=1,nlev
      do j=1,nlon
        do i=1,mlat
          delz_host(i,j,k) = z_half_host(i,j,k) - z_half_host(i,j,k+1)
          rho_host(i,j,k) = delp_host(i,j,k)/(delz_host(i,j,k)*grav)
        enddo
      enddo
    enddo
  
    do k=1,nlev_five
      do j=1,nlon
        do i=1,mlat
          delz_five(i,j,k) = z_half_five0(i,j,k) - z_half_five0(i,j,k+1)
          rho_five(i,j,k) = delp_five0(i,j,k)/(delz_five(i,j,k)*grav)
        enddo
      enddo
    enddo

    call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)
    
    !interpolate from five grid to host grid
    do j=1,nlon
      do i=1,mlat
    
        call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
        dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(t_dt_host(i,j,:)),t_dt_five0(i,j,:)) 
  
        call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
        dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(u_dt_host(i,j,:)),u_dt_five0(i,j,:)) 
  
        call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
        dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(v_dt_host(i,j,:)),v_dt_five0(i,j,:)) 
        
        do p = 1, nt_prog
          call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
          dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(r_dt_host(i,j,:, p)),r_dt_five0(i,j,:, p)) 
        enddo

        do p = nt_prog + 1, ncnst
          call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
          dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(rdiag_host(i,j,:, p)),rdiag_five0(i,j,:, p)) 
        enddo

        call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
        dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(Rad_flux_block%tdt_rad(i,j,:)),Rad_flux_five_block%tdt_rad(i,j,:)) 
          
        call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
        dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(Rad_flux_block%tdt_lw(i,j,:)),Rad_flux_five_block%tdt_lw(i,j,:)) 
          
        call tendency_low_to_high(dble(z_full_host(i,j,:)),dble(z_half_host(i,j,:)),dble(z_full_five(i,j,:)),&
        dble(rho_host(i,j,:)),dble(rho_five(i,j,:)),dble(Rad_flux_block%extinction(i,j,:)),Rad_flux_five_block%extinction(i,j,:)) 
      enddo
    enddo

    z_full_host => null()
    z_half_host => null()
    u_dt_host => null()
    v_dt_host => null()
    t_dt_host => null()
    r_dt_host => null()
    rdiag_host => null()

    z_half_five0 => null()
    delp_five0 => null()
    u_dt_five0 => null()
    v_dt_five0 => null()
    t_dt_five0 => null()
    r_dt_five0 => null()
    rdiag_five0 => null()

end subroutine five_tend_low_to_high

subroutine five_tend_high_to_low (Physics_five_input_block, Physics_tendency_five_block, &
  Physics_input_block, Physics_tendency_block)

  type(physics_input_block_type), intent(in)      :: Physics_five_input_block
  type (physics_tendency_block_type), intent(in)  :: Physics_tendency_five_block

  type(physics_input_block_type), intent(inout)      :: Physics_input_block
  type (physics_tendency_block_type), intent(inout)  :: Physics_tendency_block
    ! Local variables	
  integer :: nt_tot, nt_prog

!---------------------------------------------------------------------
!    set up local pointers into the physics input and physics tendency
!    blocks.
!---------------------------------------------------------------------
  real, dimension(:,:,:), pointer :: z_full_host, z_half_host
  real, dimension(:,:,:), pointer :: u_dt_host, v_dt_host, t_dt_host
  real, dimension(:,:,:,:), pointer :: r_dt_host, rdiag_host

  real, dimension(:,:,:), pointer :: z_half_five0
  real, dimension(:,:,:), pointer :: u_dt_five0, v_dt_five0, t_dt_five0
  real, dimension(:,:,:,:), pointer :: r_dt_five0, rdiag_five0

  z_full_host => Physics_input_block%z_full
  z_half_host => Physics_input_block%z_half
  u_dt_host => Physics_tendency_block%u_dt
  v_dt_host => Physics_tendency_block%v_dt
  t_dt_host => Physics_tendency_block%t_dt
  r_dt_host => Physics_tendency_block%q_dt
  rdiag_host => Physics_tendency_block%qdiag

  z_half_five0 => Physics_five_input_block%z_half
  u_dt_five0 => Physics_tendency_five_block%u_dt
  v_dt_five0 => Physics_tendency_five_block%v_dt
  t_dt_five0 => Physics_tendency_five_block%t_dt
  r_dt_five0 => Physics_tendency_five_block%q_dt
  rdiag_five0 => Physics_tendency_five_block%qdiag

  call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)
  !average from host grid to five grid
  do j=1,nlon
    do i=1,mlat
      
      call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
      p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), t_dt_five0(i,j,:),t_dt_host(i,j,:))

      call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
      p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), u_dt_five0(i,j,:),u_dt_host(i,j,:))

      call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
      p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), v_dt_five0(i,j,:),v_dt_host(i,j,:))

      do p = 1, nt_prog
        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
        p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), r_dt_five0(i,j,:,p),r_dt_host(i,j,:, p))
      enddo

      do p = nt_prog + 1, ncnst
        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
        p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), rdiag_five0(i,j,:,p),rdiag_host(i,j,:,p))
      enddo
    enddo
  enddo

  z_full_host => null()
  z_half_host => null()
  u_dt_host => null()
  v_dt_host => null()
  t_dt_host => null()
  r_dt_host => null()
  rdiag_host => null()

  z_half_five0 => null()
  u_dt_five0 => null()
  v_dt_five0 => null()
  t_dt_five0 => null()
  r_dt_five0 => null()
  rdiag_five0 => null()
end subroutine five_tend_high_to_low

subroutine atmosphere_pref_five (p_ref_five)
  real, dimension(:,:), intent(inout) :: p_ref_five

  !--- allocate pref
  allocate(pref_five(nlev_five+1,2), dum1d(nlev_five+1))

  !---------- reference profile -----------
  pref_five(nlev_five+1,1) = 101325.
  pref_five(nlev_five+1,2) = 81060.

  call get_eta_level_five ( nlev_five, pref_five(nlev_five+1,1), pref_five(1,1), dum1d )
  call get_eta_level_five ( nlev_five, pref_five(nlev_five+1,2), pref_five(1,2), dum1d )

  p_ref_five = pref_five

end subroutine atmosphere_pref_five

!-----------------------------------------------------------------------
subroutine get_eta_level_five(km, p_s, pf, ph, pscale)

integer, intent(in) :: km
real, intent(in)  :: p_s            ! unit: pascal
real, intent(out) :: pf(km)
real, intent(out) :: ph(km+1)
real, intent(in), optional :: pscale

integer k
real    ak1

if(nlev_five /= km)  &
call error_mesg('get_eta_level:','dimensionally inconsistent', FATAL)

    ph(1) = ak_five (1)

do k=2,nlev_five+1
    ph(k) = ak_five(k) + bk_five(k)*p_s
enddo

if ( present(pscale) ) then
    do k=1,nlev_five+1
        ph(k) = pscale*ph(k)
    enddo
endif

if( ak_five(1) > ptop_min ) then
    pf(1) = (ph(2) - ph(1)) / log(ph(2)/ph(1))
else
    ak1 = (kappa + 1.) / kappa
    pf(1) = (ph(2) - ph(1)) / ak1
!    if(master) write(*,*) 'Modified p_full(1)=',pf(1), 0.5*(ph(1)+ph(2))
endif

do k=2,nlev_five
    pf(k) = (ph(k+1) - ph(k)) / log(ph(k+1)/ph(k))
enddo

end subroutine get_eta_level_five

subroutine five_var_high_to_low (varin, varout)
  ! #include "fv_arrays.h"

  real, dimension(:,:,:), intent(in)  :: varin
  real, dimension(:,:,:), intent(out) :: varout

  do j=1,nlon
    do i=1,mlat
        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
        p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), varin(i,j,:),varout(i,j,:))
    enddo
  enddo
end subroutine five_var_high_to_low

subroutine five_var_high_to_low_4d (varin, Physics_input_block, Physics_five_input_block, varout)
  ! #include "fv_arrays.h"

  real, dimension(:,:,:,:), intent(in)  :: varin
  real, dimension(:,:,:,:), intent(out) :: varout
  type(physics_input_block_type), intent(in)      :: Physics_five_input_block
  type(physics_input_block_type), intent(in)      :: Physics_input_block
  
  real, dimension(:,:,:), pointer :: z_half_host, z_half_five0, delp_host,delp_five0

  p_full_host = Physics_input_block%p_full
  p_half_host = Physics_input_block%p_half
  p_full_five = Physics_five_input_block%p_full

  z_half_host => Physics_input_block%z_half
  z_half_five0 => Physics_five_input_block%z_half
  delp_host => Physics_input_block%delp
  delp_five0 => Physics_five_input_block%delp

  do k=1,nlev
    do j=1,nlon
      do i=1,mlat
        delz_host(i,j,k) = z_half_host(i,j,k) - z_half_host(i,j,k+1)
        rho_host(i,j,k) = delp_host(i,j,k)/(delz_host(i,j,k)*grav)
      enddo
    enddo
  enddo

  do k=1,nlev_five
    do j=1,nlon
      do i=1,mlat
        delz_five(i,j,k) = z_half_five0(i,j,k) - z_half_five0(i,j,k+1)
        rho_five(i,j,k) = delp_five0(i,j,k)/(delz_five(i,j,k)*grav)
      enddo
    enddo
  enddo

  do j=1,nlon
    do i=1,mlat
      do p = 1, size(varin, 4)
        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),delz_host(i,j,:),delz_five(i,j,:),&
        p_half_host(i,j,:),p_full_five(i,j,:),p_full_host(i,j,:), varin(i,j,:,p),varout(i,j,:,p))

      enddo
    enddo
  enddo
end subroutine five_var_high_to_low_4d

subroutine five_pressure_init(ph0, ak,bk,ps0, nlev_five)

    ! Purpose is to initialize FIVE heights.  This is called 
    !   xxx.  We need to do this here to initialize
    !   height coordinates for output purposes, so output
    !   can be written to CAM history tapes.  
    !   NOTE: that the pressure values determined here
    !   will not necessarily be the values used in the 
    !   simulation.  This is because pressure levels can 
    !   vary from timestep to timestep based on surface pressure.
    !   Here we will find the hybrid points
    !   and the actual pressure values will be computed in
    !   init_five_profiles.
    
    ! Input variables
    real, intent(in) :: ph0(nlev + 1)  ! hybrid points from host
    real, intent(in) :: ak(nlev + 1)  ! hybrid points from host
    real, intent(in) :: bk(nlev + 1)
    
    real, intent(in) :: ps0
    integer, intent(inout) :: nlev_five
    
    ! Local variables
    real, allocatable :: ak_tmp_five(:)
    real, allocatable :: bk_tmp_five(:)

    integer :: low_ind, high_ind
    integer :: kh, ki, k, i
    real :: incr, incr2, incr3
    
    low_ind = nlev
    high_ind = 1   
    
    ! First we need to determine how many nlev_five and
    !   and pverp_five levels there are.  This may seem repetitive 
    !   with the code below it but is needed to appease E3SM 
    !   order of operations and to allow us to add levels
    !   on the fly
    kh=0
    do k=1,nlev + 1
    
    kh=kh+1
    if (k .eq. nlev + 1) goto 50
    ! Test to see if this is within the bounds of the grid
    ! we want to add layers to
    if (ph0(k) .le. five_bot_toadd .and. &
    ph0(k) .ge. five_top_toadd) then
    
    do ki=1,five_add_nlevels
    kh=kh+1
    enddo
    
    endif
    enddo
    
    50  continue
    
    ! Set nlev_five and pverp_five
    nlev_five=kh-1
    
    allocate(pf0_five(nlev_five)) ! midpoint FIVE pressures (pascals)
    allocate(ph0_five(nlev_five + 1)) ! interface FIVE pressures (pascals)
    allocate(ak_tmp_five(nlev_five + 1))
    allocate(bk_tmp_five(nlev_five + 1))  
    allocate(ak_five(nlev_five + 1))
    allocate(bk_five(nlev_five + 1))              
    
    ! Note that here we are actually adding "layers" to the
    !   hybrid coordinates.  The reason is that the 
    !   pressure levels of E3SM can change from time step
    !   to timestep based on the surface pressure.  
    !   The hybrid coordinates are the only thing related to
    !   the height profile that remains constant.    
    kh=1 ! kh is layer index on FIVE grid
    do k=1,nlev + 1 ! k is layer index on E3SM grid
    
    ! Copy preexisting points to FIVE grid
    ak_tmp_five(kh) = ak(k)
    bk_tmp_five(kh) = bk(k)
    
    kh=kh+1
    if (k .eq. nlev + 1) goto 100
    ! Test to see if this is within the bounds of the grid
    ! we want to add layers to
    if (ph0(k) .le. five_bot_toadd .and. &
    ph0(k) .ge. five_top_toadd) then
    
    ! save the indicee
    if (high_ind .eq. 1) high_ind = k
    low_ind = k
    
    ! If we are inside the portion of grid we want to add layers
    ! to then compute the pressure increment (resolution)
    incr2=(ak(k+1)-ak(k))/(five_add_nlevels+1)
    incr3=(bk(k+1)-bk(k))/(five_add_nlevels+1)
    
    ! Define the new layer
    do ki=1,five_add_nlevels
    ak_tmp_five(kh) = ak_tmp_five(kh-1)+incr2
    bk_tmp_five(kh) = bk_tmp_five(kh-1)+incr3
    kh=kh+1
    enddo
    
    endif
    
    enddo
    
    100 continue
    
    ! Save the indicees of the layers
    five_bot_k = low_ind
    five_top_k = high_ind 
    
    ! Now define five interface layers.  At this point
    !  the reference and base pressure will be the same
    !  i.e. ps0 = 1.0e5_r8
    
    ph0_five(:nlev_five + 1) = 0.01*ps0*(ak_tmp_five(:nlev_five + 1) + bk_tmp_five(:nlev_five + 1)) 
    
    ! Now define five_mid layers
    do k=1,nlev_five
    pf0_five(k) = (ph0_five(k)+ph0_five(k+1))/2.0
    enddo      
    
    ! ak_five(:nlev_five + 1)  = 0.01*ps0*ak_tmp_five(:nlev_five + 1) 
    ! bk_five(:nlev_five + 1)  = 0.01*bk_tmp_five(:nlev_five + 1) 
    
    ak_five  = 0.01*ps0*ak_tmp_five
    bk_five  = 0.01*bk_tmp_five    

    write (*,*) 'ak_five', ak_five
    write (*,*) 'bk_five', bk_five

end subroutine five_pressure_init 

subroutine five_profiles_init(u_host, v_host, t_host, q_host, omga_host, pf_host)

    real, intent(in) :: u_host(nlon,mlat,nlev)
    real, intent(in) :: v_host(nlon,mlat,nlev)
    real, intent(in) :: t_host(nlon,mlat,nlev)
    real, intent(in) :: q_host(nlon,mlat,nlev,ncnst)
    real, intent(in) :: omga_host(nlon,mlat,nlev)
    real, intent(in) :: pf_host(nlon,mlat,nlev)
  
      do j=1,size(pt_five,2)
        do i=1,size(pt_five,1)
          ! Now interpolate onto the FIVE grid
          call linear_interp(pf_host(i,j,:),p_full_five(i,j,:),t_host(i,j,:),pt_five(i,j,:),nlev,nlev_five)
          call linear_interp(pf_host(i,j,:),p_full_five(i,j,:),u_host(i,j,:),ua_five(i,j,:),nlev,nlev_five)
          call linear_interp(pf_host(i,j,:),p_full_five(i,j,:),v_host(i,j,:),va_five(i,j,:),nlev,nlev_five)
          call linear_interp(pf_host(i,j,:),p_full_five(i,j,:),omga_host(i,j,:),omga_five(i,j,:),nlev,nlev_five)

          ! For Q constituents 
          do p=1,ncnst
            call linear_interp(pf_host(i,j,:),p_full_five(i,j,:),q_host(i,j,:,p),q_five(i,j,:,p),nlev,nlev_five)
          enddo
        enddo
      enddo
end subroutine five_profiles_init

!#######################################################################
! Subroutine to apply BOMEX forcings
subroutine update_bomex_forc_five()
#include "fv_arrays.h"
#include "fv_point.inc"
  
  write (*,*) 'ps', ps
  
  ! --- update pe_five, peln_five, and delp_five
  ps_five = ps
  do k=1,size(pt_five,3)
    do j=1,size(ps_five,2)
      do i=1,size(ps_five,1)
        delp_five(i,j,k) = ak_five(k+1)-ak_five(k) + ps_five(i,j)*(bk_five(k+1)-bk_five(k))
      enddo
    enddo
  enddo
  
  call p_var(nlon, mlat, nlev_five, 1, mlat, ak_five(1), delp_five, ps_five,     &
         pe_five, peln_five,  pk_five,  pkz_five,  2./7.)
  
    ! --- compute large-scale subsidence
    do k=1,nlev_five
       if ( z_full_five(1,1,k) < 1500.0 ) then
          omga_five(:,:,k) = p_full_five(1,1,k)/(rdgas*pt_five(1,1,k))*grav  &
                           * (0.0065/1500.0)*z_full_five(1,1,k)
       elseif ( z_full_five(1,1,k) < 2100.0 ) then
          omga_five(:,:,k) = p_full_five(1,1,k)/(rdgas*pt_five(1,1,k))*grav  &
                           * ( 0.0065 - (0.0065/(2100.-1500.))*(z_full_five(1,1,k)-1500.) )
       else
          omga_five(:,:,k) = 0.
       end if
    end do
end subroutine update_bomex_forc_five

subroutine tendency_low_to_high(zm_in, zi_in, &
  zm_five_in, &
  rho_low_in, rho_five_in, &
  ten_low,ten_high)

  ! Subroutine to compute the tendency from the low
  !   resolution to high resolution E3SM grid.  This code
  !   was adopted from the FIVE implementation into SAM,
  !   obtained by Peter Bogenschutz from Tak Yamaguchi.
  !   SAM code assumes that the surface is index 1, whereas
  !   E3SM assumes that the surface is index pver(p).  Thus,
  !   for simplicity, we need to "flip" the arrays.  
  ! This interpolation routine is based on that of 
  !   Sheng and Zwiers (1998) and implemented as 
  !   documented in Yamaguchi et al. (2017) appendix B.
  implicit none

  ! Input variables
  real*8, intent(in) :: zm_in(nlev) ! midpoint levels E3SM
  real*8, intent(in) :: zi_in(nlev + 1) ! interface levels E3SM
  real*8, intent(in) :: zm_five_in(nlev_five) ! midpoint levels for FIVE
  real*8, intent(in) :: ten_low(nlev)
  real*8, intent(in) :: rho_low_in(nlev)
  real*8, intent(in) :: rho_five_in(nlev_five)

  ! Output variables
  real, intent(out) :: ten_high(nlev_five)

  ! Local variables
  integer, parameter :: ml=1, mu=1, lda=2*ml+mu+1
  real*8 :: zm(nlev) ! flipped version of zm_in
  real*8 :: zi(nlev + 1) ! flipped version of zi_in
  real*8 :: df_z(nlev) 
  real*8 :: rho_low(nlev)
  real*8 :: rho_zs(nlev_five)
  real*8 :: zm_five(nlev_five) ! flipped
  real*8 :: dz ! distance of lowest level
  real*8 :: alpha
  real*8 :: adz_dn(nlev)
  real*8 :: adz_up(nlev)
  real*8 :: a(lda,nlev)
  real*8 :: adz(nlev) ! ratio of the thickness of scalar levels to dz
  real*8 :: adzw(nlev) ! ratio of the thickness of w levels to dz
  real*8 :: ipiv(nlev)
  integer :: info
  real*8 :: weight(nlev_five)
  real*8 :: rdf_z(nlev)
  real*8 :: rdf_zs_ml(nlev) ! mid-layer target value
  real*8 :: df_zs(nlev_five)
  real*8, dimension(nlev + 1) :: rdf_zm
  real*8, dimension(nlev) :: rdf_zm_dn, rdf_zm_up
  real*8, dimension(nlev) :: c0, c1, c2, c3, ic1, ic2, ic3
  real*8 :: b(nlev)

  logical, dimension(nlev) :: spurious
  logical :: cnd1, cnd2, cnd3, cnd4, cnd5

  logical :: do_limit

  integer :: i, i1, i2, i3, i4
  integer :: k, km3, km2, km1, k00
  integer :: kp1, kp2, ierr
  integer :: zi1, zi2
  integer :: i4zi1

  ! external dgbtrf
  ! external dgbtrs

  ! Construct tendency profile from E3SM grid to FIVE grid

  ! The constructed tendency profile satisfies layer mean average

  ! If no levels are added via FIVE just return a copy
  if (five_add_nlevels .eq. 0) then
  ten_high(1:nlev) = ten_low(1:nlev)
  return
  endif

  ! First let's "flip" things so lowest level index = 1
  do k=1,nlev
  zm(k) = zm_in(nlev-k+1)
  rho_low(k) = rho_low_in(nlev-k+1)
  df_z(k) = ten_low(nlev-k+1)
  enddo

  do k=1,nlev + 1
  zi(k) = zi_in(nlev + 1-k+1)
  enddo

  do k=1,nlev_five
  zm_five(k) = zm_five_in(nlev_five-k+1)
  rho_zs(k) = rho_five_in(nlev_five-k+1)
  enddo  
  
  zi1 = nlev-five_bot_k + 1
  zi2 = nlev-five_top_k + 1

  dz=zi(2)

  ! define adz and adzw
  do k=2,nlev
  adzw(k) = (zm(k)-zm(k-1))/dz
  enddo
  adzw(1) = 1._8

  adz(1) = 1._8
  do k=2,nlev-1
  adz(k) = 0.5*(zm(k+1)-zm(k-1))/dz
  enddo
  adz(nlev) = adzw(nlev)

  ! Prepare coefficients

  ! If the location of the mid-layer point is optionally specified then following variables
  ! are required to be computed every time this function is called.

  ! Distance from the mid-layer level to the host model lower/upper interface value divided
  ! by dz. Here the mid-layer level is z(k).

  do k=1,nlev
  adz_dn(k) = (zm(k) - zi(k))/dz
  enddo

  do k=1,nlev
  adz_up(k) = (zi(k+1) - zm(k))/dz
  enddo

  ! For solving system of equations
  ! The j-th column of the matrix A is stored in the j-th column of the array a as follows:
  !    a(ml+mu+1+i-j,j) = A(i,j) for max(1,j-mu)<=i<=min(nzm,j+ml)
  ! Set up superdiagonal, diagonal, subdiagonal component of the matrix
  a(:,:) = 0._8
  ! Superdiagonal 
  do k=2,nlev
  a(2,k)=adz_up(k-1)**2/(adz(k)*adzw(k-1))*0.5_8
  enddo
  ! Diagonal
  k = 1
  a(3,1) = (adz_dn(1) + adzw(1) + adz_up(1) * adz_dn(2) / adz(2))/adzw(1) * 0.5_8
  do k=2,nlev
  kp1=min(k+1,nlev)
  a(3,k) = (adz_up(k-1) * adz_dn(k) / adz(k) + adzw(k) & 
    + adz_up(k) * adz_dn(kp1) / adz(kp1) ) / adzw(k) * 0.5
  !+ adz_dn(kp1) used to be adz_dn(k)
  enddo
  ! Subdiagonal
  do k=1, nlev
  a(4,k) = adz_dn(k)**2 / (adz(k) * adzw(k))*0.5
  enddo

  ! Factor the matrix with LAPACK, BLAS
  call dgbtrf( nlev, nlev, ml, mu, a, lda, ipiv, info )    

  ! For interpolation
  i4=0
  do k=1,zi1-1
  i4 = i4 + 1
  enddo

  c0(:) = 0._8
  c1(:) = 0._8
  c2(:) = 0._8
  c3(:) = 0._8
  ic1(:) = 0._8
  ic2(:) = 0._8
  ic3(:) = 0._8

  weight(:) = 0._8

  i4zi1 = i4
  do k = zi1, zi2
  i1 = i4 + 1
  if (mod(five_add_nlevels,2) .ne. 0 ) then
  i2 = i1 + (five_add_nlevels+1) / 2 - 1 
  else
  i2 = i1 + five_add_nlevels / 2 - 1 
  end if
  i3 = i2 + 1
  i4 = i1 + five_add_nlevels 
  ! weight for linear interpolation
  weight(i1:i2) = ( zm_five(i1:i2) - zi(k) ) / ( zm(k) - zi(k) )
  weight(i3:i4) = ( zm_five(i3:i4) - zi(k+1) ) / ( zm(k) - zi(k+1) ) 
  ! c1, c2, c3
  c0(k) = 2.0_8 * (zi(k+1) - zi(k))
  c1(k) = zi(k+1) - zi(k)
  c2(k) = zm(k) - zi(k)
  c3(k) = zi(k+1) - zm(k)

  ic1(k) = 1.0_8 / c1(k)
  ic2(k) = 1.0_8 / c2(k) 
  ic3(k) = 1.0_8 / c3(k)            

  enddo

  ! add flag computed?

  ! Mass weight inout value
  rdf_z(:) = rho_low(:) * df_z(:)

  ! Solve system of equations to get mid-layer target value
  ! Solve the linear system with LAPACK, BLAS
  ! input b wil be solution when output
  b(:) = rdf_z(:)
  call dgbtrs('n', nlev, ml, mu, 1, a, lda, ipiv, b, nlev, info)
  rdf_zs_ml(:) = b(:)

  ! Interface target value
  rdf_zm(1) = rdf_zs_ml(1)
  do k = 2, nlev
  rdf_zm(k) = adz_dn(k) / adz(k) * rdf_zs_ml(k-1) & 
  + adz_up(k-1) / adz(k) * rdf_zs_ml(k)
  enddo
  rdf_zm(nlev + 1) = 0.0_8 !domain top tendency

  do_limit = .true.
  if (do_limit) then 

  ! Detection and correction of grid-scale violation for df_zm 
  !  Zerroukat et al. (2005 QJRMS)
  spurious(:) = .false.
  do k = 1, nlev
  km3 = MAX( k-3, 1 )
  km2 = MAX( k-2, 1 )
  km1 = MAX( k-1, 1 )
  k00 = MIN( k, nlev )
  kp1 = MIN( k+1, nlev )
  kp2 = MIN( k+2, nlev )
  cnd1 = ( rdf_zm(k) - rdf_z(km1) ) * ( rdf_z(k00) - rdf_zm(k) ) < 0.0
  cnd2 = ( rdf_z(km1) - rdf_z(km2) ) * ( rdf_z(kp1) - rdf_z(k00) ) >= 0.0
  cnd3 = ( rdf_z(km1) - rdf_z(km2) ) * ( rdf_z(km2) - rdf_z(km3) ) <= 0.0
  cnd4 = ( rdf_z(kp1) - rdf_z(k00) ) * ( rdf_z(kp2) - rdf_z(kp1) ) <= 0.0
  cnd5 = ( rdf_zm(k) - rdf_z(km1) ) * ( rdf_z(km1) - rdf_z(km2) ) <= 0.0
  if ( cnd1 .and. ( cnd2 .or. cnd3 .or. cnd4 .or. cnd5 ) ) then
  spurious(k) = .false.
  alpha = ABS( rdf_zm(k) - rdf_z(k00) ) - ABS( rdf_zm(k) - rdf_z(km1) )
  alpha = SIGN( 0.5_8, alpha ) + 0.5_8 ! Heaviside step function, alpha = 0 or 1
  rdf_zm(k) = alpha * rdf_z(km1) + ( 1.0_8 - alpha ) * rdf_z(k00)
  endif
  enddo

  ! Store rdf_zm into rdf_zm_up and rdf_zm_dn
  rdf_zm_up(:) = rdf_zm(2:nlev + 1)
  rdf_zm_dn(:) = rdf_zm(1:nlev)  

  ! Detection and correction of grid-scale violation for rdf_zs_ml for monotonic layer
  ! - Recompute rdf_zs_zl with updated rdf_zm
  ! - For monotonic layer, check if it is bounded between rdf_zm(k) and rdf_zm(k+1)
  !   If not bounded, assign rdf_zs_ml to closest rdf_zm, and compute other interface value
  !   and store it into either rdf_zm_up or rdf_zm_dn. That level will have two interface
  !   values for upper and lower layer.    
  spurious(:) = .FALSE.
  do k = 1, nlev
  km1 = MAX( k-1, 1 )
  kp1 = MIN( k+1, nlev )
  rdf_zs_ml(k) = ic1(k) * ( c0(k)*rdf_z(k) - c2(k)*rdf_zm(k) - c3(k)*rdf_zm(k+1) ) !+PAB change
  cnd1 = ( rdf_z(k) - rdf_z(km1) ) * ( rdf_z(kp1) - rdf_z(k) ) >= 0.0_8
  cnd2 = ( rdf_zs_ml(k) - rdf_zm(k) ) * ( rdf_zm(k+1) - rdf_zs_ml(k) ) < 0.0_8 !+PAB change
  if ( cnd1 .AND. cnd2 ) then
  ! Inflection within a monotonic layer
  spurious(k) = .TRUE.
  alpha = ABS( rdf_zs_ml(k) - rdf_zm(k) ) - ABS( rdf_zs_ml(k) - rdf_zm(k+1) ) !+PAB change
  alpha = SIGN( 0.5_8, alpha ) + 0.5_8 ! alpha = 0 or 1
  rdf_zs_ml(k) = alpha * rdf_zm(k+1) + ( 1.0_8 - alpha ) * rdf_zm(k) !+PAB change
  if ( alpha < 0.5_8 ) then
  rdf_zm_up(k) = ic3(k) * ( c0(k)*rdf_z(k) - c1(k)*rdf_zs_ml(k) - c2(k)*rdf_zm(k) )
  else
  rdf_zm_dn(k) = ic2(k) * (c0(k)*rdf_z(k) - c1(k)*rdf_zs_ml(k) - c3(k)*rdf_zm(k+1)) !+PAB change
  endif
  endif
  enddo

  ! Remove discountinuity at interface level as many as possible
  ! - For monotonic layer, set rdf_z_dn(k) = rdf_z_up(k-1) and rdf_z_up(k) = rdf_z_dn(k+1),
  !   then re-compute rdf_zs_ml.
  ! - Check if new rdf_zs_ml is bounded between rdf_zm_dn and rdf_zm_up, and if not, set
  !   rdf_zs_ml to the closer interface value and compute the other interfaec value.              
  do k = 1, nlev
  km1 = MAX( k-1, 1 )
  kp1 = MIN( k+1, nlev )
  cnd1 = ( rdf_zs_ml(k) - rdf_zm_dn(k) ) * ( rdf_zm_up(k) - rdf_zs_ml(k) ) >= 0.0_8
  if ( cnd1 ) then
  ! Monotonic layer

  ! Re-set df_zm_dn and df_zm_up
  rdf_zm_dn(k) = rdf_zm_up(km1)
  rdf_zm_up(k) = rdf_zm_dn(kp1)
  ! Re-compute df_zs
  rdf_zs_ml(k) = ic1(k) * ( c0(k)*rdf_z(k) - c2(k)*rdf_zm_dn(k) - c3(k)*rdf_zm_up(k) )
  !
  cnd2 = ( rdf_zs_ml(k) - rdf_zm_dn(k) ) * ( rdf_zm_up(k) - rdf_zs_ml(k) ) < 0.0_8

  if ( cnd2 ) then
  ! Non-monotonic profile
  spurious(k) = .TRUE.
  alpha = ABS( rdf_zs_ml(k) - rdf_zm_dn(k) ) - ABS( rdf_zs_ml(k) - rdf_zm_up(k) )
  alpha = SIGN( 0.5_8, alpha ) + 0.5_8
  rdf_zs_ml(k) = alpha * rdf_zm_up(k) + ( 1.0_8 - alpha ) * rdf_zm_dn(k)
  if ( alpha < 0.5_8 ) then
  rdf_zm_up(k) = ic3(k) * (c0(k)*rdf_z(k)-c1(k)*rdf_zs_ml(k)-c2(k)*rdf_zm_dn(k))
  else
  rdf_zm_dn(k) = ic2(k) * (c0(k)*rdf_z(k)-c1(k)*rdf_zs_ml(k)-c3(k)*rdf_zm_up(k))
  endif
  endif

  endif
  enddo

  else

  rdf_zm_up(:) = rdf_zm(2:nlev + 1)
  rdf_zm_dn(:) = rdf_zm(1:nlev)

  endif

  ! Construct the tendency profile
  i4 = 0
  ! Below zi1
  do k=1,zi1-1
  i4=i4+1
  df_zs(i4) = df_z(k)
  enddo

  ! between zi1 and zi2
  do k = zi1, zi2
  i1 = i4 + 1
  if (mod(five_add_nlevels,2) .ne. 0 ) then
  i2 = i1 + (five_add_nlevels+1) / 2 - 1 
  else
  i2 = i1 + five_add_nlevels / 2 - 1 
  end if
  i3 = i2 + 1
  i4 = i1 + five_add_nlevels 
  ! Compute df_zs for all zs levels
  ! zi(k) <= zs(i1:i2) < z(k)
  df_zs(i1:i2) = ( 1.0_8 - weight(i1:i2) ) * rdf_zm_dn(k) + weight(i1:i2) * rdf_zs_ml(k)
  ! z(k) <= zs(i3:i4) < zi(k+1)
  df_zs(i3:i4) = ( 1.0_8 - weight(i3:i4) ) * rdf_zm_up(k) + weight(i3:i4) * rdf_zs_ml(k)
  ! Non-mass weighted
  df_zs(i1:i4) = df_zs(i1:i4) / rho_zs(i1:i4)
  enddo

  ! Above zi2
  do k = zi2+1,nlev
  i4 = i4 + 1
  df_zs(i4) = df_z(k)
  enddo    

  ! Finally, "unflip" the tendency
  do k = 1,nlev_five
  ten_high(k) = df_zs(nlev_five-k+1)
  enddo
end subroutine tendency_low_to_high 

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

subroutine linear_interp(x1,x2,y1,y2,km1,km2)
    implicit none

    ! Simple linear interpolation routine
    ! WILL LIKELY BE REPLACED
    
    integer :: km1, km2, k1, k2
    real :: x1(km1), y1(km1)
    real :: x2(km2), y2(km2)

    do k2=1,km2

      if( x2(k2) <= x1(1) ) then
        y2(k2) = y1(1) + (y1(2)-y1(1))*(x2(k2)-x1(1))/(x1(2)-x1(1))
      elseif( x2(k2) >= x1(km1) ) then
        y2(k2) = y1(km1) + (y1(km1)-y1(km1-1))*(x2(k2)-x1(km1))/(x1(km1)-x1(km1-1))    
      else
        do k1 = 2,km1
          if( (x2(k2)>=x1(k1-1)).and.(x2(k2)<x1(k1)) ) then
            y2(k2) = y1(k1-1) + (y1(k1)-y1(k1-1))*(x2(k2)-x1(k1-1))/(x1(k1)-x1(k1-1))
          endif
        enddo
      endif
     
    enddo

end subroutine linear_interp  

subroutine masswgt_vert_avg(&
    rho_host,rho_high,dz_host,dz_high,&
  pint_host,pmid_high,pmid_host,&
  var_high,var_host)
  
  ! Purpose is to compute the mass weighted vertical 
  !  average from the FIVE high resolution grid onto the
  !  E3SM grid
  implicit none
  
  real, intent(in) :: rho_high(nlev_five)
  real, intent(in) :: rho_host(nlev)   
  real, intent(in) :: dz_high(nlev_five)
  real, intent(in) :: dz_host(nlev)
  real, intent(in) :: pmid_high(nlev_five)
  real, intent(in) :: pint_host(nlev + 1)
  real, intent(in) :: pmid_host(nlev)
  real, intent(in) :: var_high(nlev_five)
  real, intent(out) :: var_host(nlev)
  
  integer :: i, k, kh
  
  real :: rho_host_avg(nlev)
  
  ! Initialize host variable
  var_host(:) = 0.
  rho_host_avg(:) = 0.
  
  kh=1 ! Vertical index for FIVE 
  do k=1,nlev ! Vertical index for E3SM
  
  ! Check to see how many FIVE layers there are within
  !   an E3SM layer 
  do while ( pmid_high(kh) .lt. pint_host(k+1) .and. &
      pmid_high(kh) .gt. pint_host(k))
  
  var_host(k) = var_host(k) + rho_high(kh) * &
  var_high(kh) * dz_high(kh)
  rho_host_avg(k) = rho_host_avg(k) + rho_high(kh) * &
  dz_high(kh)
  
  kh = kh + 1 ! increase high res model by one layer
  if (kh .gt. nlev_five) goto 10
  
  end do ! end while loop for kh
  10 continue
  
  ! Compute rho on host grid
  rho_host_avg(k) = rho_host_avg(k)/dz_host(k)
  
  var_host(k) = var_host(k)/(rho_host_avg(k)*dz_host(k))
  
  enddo
  
  return
  
end subroutine masswgt_vert_avg    

subroutine p_var(im, jm, km, jfirst, jlast, ptop,    &
  delp,  ps,  pe, peln, pk, pkz,      &
  cappa)

! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
!

! Input:
integer,  intent(in):: im, jm, km               ! Dimensions
integer,  intent(in):: jfirst, jlast            ! Latitude strip

real, intent(in):: ptop
real, intent(in):: cappa
real, intent(inout):: delp(im,jfirst:jlast, km)
! Output:
real, intent(out) ::   ps(im, jfirst:jlast)
real, intent(out) ::   pk(im, jfirst:jlast, km+1)
real, intent(out) ::   pe(im, km+1, jfirst:jlast)    ! Edge pressure
real, intent(out) :: peln(im, km+1, jfirst:jlast)    ! Edge pressure
real, intent(out) ::  pkz(im, jfirst:jlast, km)

! Local
real pek
real lnp
real ak1
integer i, j, k

pek = ptop ** cappa
ak1 = (cappa + 1.) / cappa

!omp parallel do private(i,j,k)
do j=jfirst,jlast !jfirst,jlast
do i=1,im
pe(i,1,j) = ptop
pk(i,j,1) = pek
enddo

do k=2,km+1
do i=1,im
pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
peln(i,k,j) = log(pe(i,k,j))
pk(i,j,k) = pe(i,k,j)**cappa
enddo
enddo

do i=1,im
ps(i,j) = pe(i,km+1,j)
enddo

!---- GFDL modification
if( ptop < 1.E-6 ) then
do i=1,im
peln(i,1,j) = peln(i,2,j) - ak1
enddo
else
lnp = log( ptop )
do i=1,im
peln(i,1,j) = lnp
enddo
endif
!---- GFDL modification

do k=1,km
do i=1,im
pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
enddo
enddo
enddo

end subroutine p_var

end module five_mod