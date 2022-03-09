module scm_rico_mod
 
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

   use   surface_flux_mod, only:  surface_flux

   use sat_vapor_pres_mod, only:  lookup_es, compute_qs, sat_vapor_pres_init   ! ZNT 02/20/2020
   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav, pi, seconds_per_day

   use      scm_utils_mod, only:  us_std_atm

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
   include 'netcdf.inc'
   private

   public rico_data_read, rico_forc_init, rico_forc_end, update_rico_forc, &
          rico_forc_diagnostic_init, rico_sfclyr,   rico_sfcflux_online

character(len=8) :: mod_name = 'scm_rico'
character(len=7) :: mod_name_diag = 'forcing'

real, public                                 :: psfc
real, parameter                         :: zsfc =  0.0
real, public, allocatable, dimension(:,:,:)  :: u_geos, v_geos

! Forcing data

integer :: rcode, ncid, dimid(12), dimid_time, dimid_level_f, varid

integer :: ntime, nlevel
type(time_type), allocatable, dimension(:) :: time_forc
integer,      allocatable, dimension(:)    :: date
real(kind=4), allocatable, dimension(:)    :: tmp1d
real(kind=4), allocatable, dimension(:,:)  :: tmp2d
real, allocatable, dimension(:)            :: time, ps_forc
real, allocatable, dimension(:,:)          :: pres_forc, dtdt_forc, dqvdt_forc, &
                                              dudt_forc, dvdt_forc,             &
                                              u_geos_forc, v_geos_forc

real, allocatable, dimension(:)            :: plev_forc, field_forc_less, field_forc_more

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       PARAMETERS OF THE MODULE
!
!
!       p00                      reference pressure (pascals)
!
!       configuration            case configuration ('short' for composite
!                                case, 'long' for long case)
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

real,    private               :: p00 = 100000.

character(len=64)              :: configuration = 'short'

integer, public                :: tracer_vert_advec_scheme = 3
integer, public                :: temp_vert_advec_scheme = 3
integer, public                :: momentum_vert_advec_scheme = 3
logical, public                :: vert_advec_cond = .true.

logical, public                :: do_rad = .true.
logical, public                :: do_geo = .true.
logical, public                :: do_vadv = .true.
logical, public                :: do_nudge = .true.

real,    private               :: missing_value = -999.

logical, private               :: initialized = .false.
logical                            :: do_netcdf_restart = .true.

integer :: nsphum, nql, nqi, nqa, nqn, nqni
integer :: vadvec_scheme

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
integer, dimension(1) :: restart_versions = (/ 1 /)

! namelist

namelist /scm_rico_nml/ tracer_vert_advec_scheme,   &
                        temp_vert_advec_scheme,     &
                        momentum_vert_advec_scheme, &
                        vert_advec_cond,            &
                        do_rad, do_geo, do_vadv,    &
                        do_nudge, configuration

! diagnostics

integer ::  id_tdt_radf, id_tdt_vadv, id_tdt_lf,                &
            id_qvdt_vadv, id_qvdt_lf,                           &
            id_qldt_vadv, id_qidt_vadv, id_qadt_vadv,     id_qndt_vadv,        &
            id_udt_vadv, id_udt_geos, id_udt_lf, id_udt_nudge,  &
            id_vdt_vadv, id_vdt_geos, id_vdt_lf, id_vdt_nudge,  &
            id_flx_radf, id_zi_forc,                            &
            id_pf_forc, id_ph_forc, id_zf_forc, id_zh_forc,     &
            id_u_geos, id_v_geos
	    
! ---> h1g, 2010-07-19
integer ::    id_qvdt_forc_col
integer ::    id_qldt_vadv_col
! <--- h1g, 2010-07-19

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!#######################################################################
! Subroutine to read case specific namelist

subroutine rico_data_read(kmax)

implicit none

integer,  intent (in)                     :: kmax
integer                                   :: i
integer                                   :: unit,ierr,io, logunit
character*23                              :: tracer_ascheme,temp_ascheme
character*64                              :: fname_res='INPUT/rico.res.nc'

integer :: year, month, day

   if (initialized) return
   initialized = .true.
      
!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   READ (input_nml_file, nml=scm_rico_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_rico_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=scm_rico_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_rico_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=scm_rico_nml)
   endif

!--------- allocate memory ---------

   if (allocated(u_geos)) deallocate(u_geos);  allocate(u_geos(1,1,kmax))
   if (allocated(v_geos)) deallocate(v_geos);  allocate(v_geos(1,1,kmax))

! ---> h1g, 2010-10-07
    if (trim(configuration)=='short') then
      psfc = 1015.4e2
    elseif (trim(configuration)=='long') then
      psfc = 1014.0e2
    else
      call error_mesg('rico_forc_init','invalid configuration', FATAL)
    endif
    
   if (file_exist('INPUT/rico.res.nc') ) then
      if(mpp_pe() == mpp_root_pe() ) call mpp_error ('scm_rico_mod', &
         'Reading netCDF formatted restart file: INPUT/rico.res.nc', NOTE)
      if( allocated(u_geos) ) call read_data (fname_res, 'u_geos',  u_geos)
      if( allocated(v_geos) ) call read_data (fname_res, 'v_geos',  v_geos)
   endif
! <--- h1g, 2010-10-07
       
!-----------------------------------------------------------------------
! Read forcing data if for 'long' configuration

   if (trim(configuration)=='long') then

!    Open netCDF file

     rcode = nf_open('INPUT/lsforcings.nc',NF_SHARE,ncid)
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot open INPUT/lsforcings.nc', FATAL)

!    Read dimension sizes

     rcode = nf_inq_dimid(ncid, 'time', dimid_time)
     rcode = nf_inq_dimlen(ncid, dimid_time, ntime)
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read dimension time', FATAL)

     rcode = nf_inq_dimid(ncid, 'level_f', dimid_level_f)
     rcode = nf_inq_dimlen(ncid, dimid_level_f, nlevel)
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read dimension level_f', FATAL)

!    Allocate scracth arrays of rel(kind=4) type

     if (allocated(tmp1d)) deallocate (tmp1d)
     allocate(tmp1d(ntime))
     if (allocated(tmp2d)) deallocate (tmp2d)
     allocate(tmp2d(nlevel,ntime))
  
!    Read time

     if (allocated(time)) deallocate (time)
     allocate(time(ntime))
     rcode = nf_inq_varid(ncid, 'time', varid)
     rcode = nf_get_var_real(ncid, varid, tmp1d)
     time = tmp1d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable time', FATAL)

!    Read date

     if (allocated(date)) deallocate (date)
     allocate(date(ntime))
     rcode = nf_inq_varid(ncid, 'date', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable date', FATAL)
     endif
     rcode = nf_get_var_int(ncid, varid, date)
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable date', FATAL)

!    Read ps_forc

     if (allocated(ps_forc)) deallocate (ps_forc)
     allocate(ps_forc(ntime))
     rcode = nf_inq_varid(ncid, 'surfpres', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable surfpres', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp1d)
     ps_forc = tmp1d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable ps_forc', FATAL)

!    Read pres_forc

     if (allocated(pres_forc)) deallocate (pres_forc)
     allocate(pres_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'pres', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable pres', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     pres_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable pres', FATAL)

!    Read dtdt_forc

     if (allocated(dtdt_forc)) deallocate (dtdt_forc)
     allocate(dtdt_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'dtdt_ls', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable dtdt_ls', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     dtdt_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable dtdt_ls', FATAL)

!    Read dqvdt_forc

     if (allocated(dqvdt_forc)) deallocate (dqvdt_forc)
     allocate(dqvdt_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'dqvdt_ls', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable dqvdt_ls', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     dqvdt_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable dqvdt_ls', FATAL)

!    Read dudt_forc

     if (allocated(dudt_forc)) deallocate (dudt_forc)
     allocate(dudt_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'dudt_ls', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable dudt_ls', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     dudt_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable dudt_ls', FATAL)

!    Read dvdt_forc

     if (allocated(dvdt_forc)) deallocate (dvdt_forc)
     allocate(dvdt_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'dvdt_ls', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable dvdt_ls', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     dvdt_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable dvdt_ls', FATAL)

!    Read u_geos_forc

     if (allocated(u_geos_forc)) deallocate (u_geos_forc)
     allocate(u_geos_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'ugeo', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable ugeo', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     u_geos_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable ugeo', FATAL)

!    Read v_geos_forc

     if (allocated(v_geos_forc)) deallocate (v_geos_forc)
     allocate(v_geos_forc(nlevel,ntime))
     rcode = nf_inq_varid(ncid, 'vgeo', varid)
     rcode = nf_inq_vardimid(ncid, varid, dimid)
     if (dimid(1)/=dimid_level_f .or. dimid(2)/=dimid_time) then
       call error_mesg ('scm_rico_mod','dimension error for variable vgeo', FATAL)
     endif
     rcode = nf_get_var_real(ncid, varid, tmp2d)
     v_geos_forc = tmp2d
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot read variable vgeo', FATAL)

!    Close netCDF file

     rcode = nf_close(ncid)
     if (rcode/=0) call error_mesg ('scm_rico_mod','cannot close INPUT/lsforcings.nc', FATAL)
     deallocate( tmp1d, tmp2d )

!    Define time_forc array

     if (allocated(time_forc)) deallocate (time_forc)
     allocate(time_forc(ntime))
  
     year = date(1)/10000
     month = (date(1) - 10000*year)/100
     day = date(1) - 10000*year - 100*month
     if (int(time(1))/=0) call error_mesg ('scm_rico_mod','time(1) should be zero', FATAL)
     time_forc(1) = set_date(year, month, day)
     do i=2,ntime
       time_forc(i) = increment_date(time_forc(1), 0, 0, 0, 0, 0, int(time(i)), 0)
     end do
  
!    Allocate scratch arrays needed for interpolation

     if (allocated(plev_forc)) deallocate(plev_forc)
     allocate( plev_forc(nlevel) )
     if (allocated(field_forc_less)) deallocate(field_forc_less)
     allocate( field_forc_less(nlevel) )
     if (allocated(field_forc_more)) deallocate(field_forc_more)
     allocate( field_forc_more(nlevel) )

   endif
        
end subroutine rico_data_read

!#######################################################################
! Subroutine to initialize case forcings

subroutine rico_forc_init(time_interp, As, Bs)
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
!      pt           temperature in degrees Kelvin 
!      u            zonal wind (m/s)
!      v            meridional wind (m/s)
!      qv           water vapor specific humidity (kg water vapor/     &
!                                                  kg air)
!      ql           liquid water condensate specific humidity 
!                   (kg condensate/ kg air)
!      qi           ice H2O condensate specific humidity 
!                   (kg condensate/kg dry air)
!      qa           cloud fraction (fraction)

type(time_type)                          :: time_interp
real,  intent (in), dimension(:)         :: As,Bs


!  Internal variables
!  ------------------

integer  :: kdim, k

integer, parameter :: itmax = 10

real, dimension(size(pt,3)+1) :: eta, peta
real, dimension(size(pt,3))   :: T_us_std, qv_us_std
real, dimension(size(pt,3))   :: u_rico, v_rico, ug_rico, vg_rico, T_rico, qv_rico

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

    if (trim(configuration)=='short') then
      psfc = 1015.4e2
    elseif (trim(configuration)=='long') then
      psfc = 1014.0e2
    else
      call error_mesg('rico_forc_init','invalid configuration', FATAL)
    endif
    ps = psfc
    elev = zsfc
!   ZNT 02/20/2020: initialize vapor pressure table
    call sat_vapor_pres_init

!   Setup temporary vertical grid: this grid is needed in order to use
!   the subroutine 

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
    !        call error_mesg ('rico_forc_init ',   &
    !        'need to compute height', FATAL )

!   Iteration to compute RICO sounding

    i=0
    maxerror = 1.0
    do while ( maxerror > 0.001 .and. i < itmax )

     i = i + 1
     do k=1,kdim

      call rico_snd( configuration, zf(1,1,k),                      & 
                     u_rico(k), v_rico(k), ug_rico(k), vg_rico(k),  &
                     T_rico(k), qv_rico(k) )
      u_geos(:,:,k) = ug_rico(k)
      v_geos(:,:,k) = vg_rico(k)
      pt(:,:,k)  = T_rico(k)
      ua(:,:,k) = u_rico(k)
      va(:,:,k) = v_rico(k)
      q(:,:,k,nsphum) = 0.001 * qv_rico(k)
      q(:,:,k,nql) = 0.0
     enddo

     call compute_p_z (nlev, grav*elev, pe, peln, delp, pt, q(:,:,:,nsphum), pf, ph, zfnew, zhnew)

     ! call compute_height (Vgrid, grav*elev, T, qv, pf, ph, zfnew, zhnew )
     !       call error_mesg ('rico_forc_init ',   &
     !       'need to compute height', FATAL )

     maxerror = maxval( abs( zfnew - zf ) )
     zf = zfnew
     zh = zhnew

    enddo

    if ( i >= itmax ) then
      call error_mesg('rico_forc_init',  &
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

end subroutine rico_forc_init

!#######################################################################

subroutine rico_forc_end ()
character*64                 :: fname_res='RESTART/rico.res.nc'

  if (.not.initialized) return
  initialized = .false.
  
 ! ---> h1g, 2010-10-07
     if( do_netcdf_restart ) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('scm_rico_mod', 'Writing netCDF formatted restart file: RESTART/rico.res.nc', NOTE)
           endif
           call write_data (fname_res, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)	   
           call write_data (fname_res,  'u_geos',  u_geos)
           call write_data (fname_res,  'v_geos',  v_geos)
      endif
! <--- h1g, 2010-10-07

  deallocate ( u_geos, v_geos )
 
  if (trim(configuration)=='long') then

    deallocate( date )
    deallocate( time, ps_forc )
    deallocate( pres_forc, dtdt_forc, dqvdt_forc )
    deallocate( dudt_forc, dvdt_forc )
    deallocate( u_geos_forc, v_geos_forc )
    deallocate( time_forc )
    deallocate( plev_forc, field_forc_less, field_forc_more )

  endif

end subroutine rico_forc_end

!#######################################################################

subroutine rico_forc_diagnostic_init(axes, Time)

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
   'liquid droplet number concentration tendencies due to vertical advection', '1/cm3/s', &
    missing_value = missing_value)

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
     
! ---> h1g, 2010-07-19
 id_qvdt_forc_col =  register_diag_field (mod_name_diag, 'qvdt_forc_col', axes(1:2), Time, &
     'column integrated vapor forcing', 'kg/m2/s',  missing_value = missing_value)
 id_qldt_vadv_col =  register_diag_field (mod_name_diag, 'qldt_vadv_col', axes(1:2), Time, &
     'column integrated cloud water vertical advection', 'kg/m2/s',  missing_value = missing_value)
! <--- h1g, 2010-07-19

end subroutine rico_forc_diagnostic_init

!#######################################################################
! Subroutine to apply RICO forcings

subroutine update_rico_forc(time_interp,time_diag,dt_int)
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
real,  dimension(size(pt,1),size(pt,2))    :: elev

integer                                          :: i,j,k,kdim
integer                                          :: dt_seconds,dt_days
logical                                          :: used
real                                             :: fcriolis, fnudge, dts
real, dimension(size(pt,1),size(pt,2))             :: zi
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: pf,zf
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dp
real, dimension(size(pt,1),size(pt,2),size(pt,3)+1) :: ph, omega_h, zh, frad
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: dT_adi, dT_vadv, dT_lf, dT_rad, &
                                                       dqv_vadv, dqv_lf, dql_vadv, dqi_vadv, dqa_vadv, dqn_vadv
real, dimension(size(pt,1),size(pt,2),size(pt,3))   :: du_vadv, dv_vadv, du_geos, dv_geos, du_lf, dv_lf, du_nudge, dv_nudge

real, dimension(size(pt,3))                       :: tmp_ans

integer :: itime_less,  itime_more
real    :: weight_less, weight_more

! ---> h1g, 2010-07-19
real, dimension(size(pt,1),size(pt,2))  ::  qvdt_forcing_col
real, dimension(size(pt,1),size(pt,2))  ::  qldt_vadv_col
integer      ::   current_sec0,  current_days0

#include "fv_point.inc"
! <--- h1g, 2010-07-19


kdim = size(pt,3)

if (trim(configuration)=='short') then

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
       !     call error_mesg ('rico_forc_init ',   &
       !     'need to compute height', FATAL )
 
   ! --- compute large-scale subsidence
   do k=1,kdim
      if ( zf(1,1,k) < 2260.0 ) then
         omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav * (0.005/2260.0)*zf(1,1,k)
      elseif ( zf(1,1,k) < 4000.0 ) then
         omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav * (0.005)

! ---> h1g
      elseif ( zf(1,1,k) < 5000.0 ) then
         omga(:,:,k) = pf(1,1,k)/(rdgas*pt(1,1,k))*grav * (0.005 - (0.005/(5000.-4000.))*(zf(1,1,k)-4000.))
! <--- h1g

      else
         omga(:,:,k) = 0.
      end if
   end do

   omega_h = 0.0
   do k=2,kdim
      if ( zh(1,1,k) < 2260.0 ) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                          * (0.005/2260.0)*zh(1,1,k)
      elseif ( zh(1,1,k) < 4000.0 ) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                          * (0.005)

! ---> h1g
      elseif ( zh(1,1,k) < 5000.0 ) then
         omega_h(:,:,k) = ph(1,1,k)/(rdgas*0.5*(pt(1,1,k)+pt(1,1,k-1)))*grav &
                          * (0.005 - (0.005/(5000.-4000.))*(zh(1,1,k)-4000.))
! <--- h1g

      else
         omega_h(:,:,k) = 0.
      end if
   end do

   ! --- compute dp
   do k = 2,kdim+1
      dp(:,:,k-1) = ph(:,:,k) - ph(:,:,k-1)
   enddo

   ! --- radiative cooling: included in the large-scale temperature forcing
   dT_rad=0.

   ! --- large-scale forcing tendencies
   du_lf = 0.0
   dv_lf = 0.0

   dT_lf = 0.0
   do k=1, kdim
      if ( zf(1,1,k) < 4000.0 ) then
         dT_lf(1,1,k) = -2.51/86400.0 + (-2.18+2.51)/(86400.0*4000.0) * zf(1,1,k)

! ---> h1g
       else if ( zf(1,1,k) < 5000.0 ) then
         dT_lf(1,1,k) = -2.18/86400.0 + (2.18)/(86400.0*(5000.0-4000.0)) * (zf(1,1,k)-4000.0)
! <--- h1g
      else
         dT_lf(1,1,k) = 0.0
      end if
   end do

   dqv_lf = 0.0
   do k=1, kdim
      if ( zf(1,1,k) < 3000.0 ) then
         dqv_lf(1,1,k) = 0.001*( -1.0/86400.0 + (0.345+1.0)/(86400.0*3000.0) * zf(1,1,k) )
      else if ( zf(1,1,k) < 4000.0 ) then
         dqv_lf(1,1,k) = 0.001*( 0.345/86400.0 )

! ---> h1g
      else if ( zf(1,1,k) < 5000.0 ) then
         dqv_lf(1,1,k) = 0.001*( 0.345/86400.0 + (-0.345)/(86400.0*(5000.0-4000.0)) * (zf(1,1,k)-4000.0) )
! <--- h1g

      else
         dqv_lf(1,1,k) = 0.0
      end if
   end do
  
   ! --- large-scale subsidence tendencies
   dT_vadv=0.0; dT_adi=0.0; dqv_vadv=0.0; dql_vadv=0.0; dqi_vadv=0.0;dqa_vadv=0.0;dqn_vadv=0.0;
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
     if( nqn >0 ) &
       call vert_advection(dts,omega_h,delp,q(:,:,:,nqn),dqn_vadv,scheme=vadvec_scheme,form=ADVECTIVE_FORM)

      select case (momentum_vert_advec_scheme)
      case(1)
         call vert_advection(dts,omega_h,delp,ua,dU_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
         call vert_advection(dts,omega_h,delp,va,dV_vadv,scheme=SECOND_CENTERED,form=ADVECTIVE_FORM)
      end select

   end if

   ! --- nudge winds
   du_nudge = 0.0
   dv_nudge = 0.0

elseif (trim(configuration)=='long') then

   omga=0.0; omega_h=0.0;
   dT_rad=0.0
   du_lf=0.0; dv_lf=0.0; dT_lf=0.0; dqv_lf=0.0;
   dT_vadv=0.0; dT_adi=0.0; dqv_vadv=0.0; dql_vadv=0.0; dqi_vadv=0.0;dqa_vadv=0.0; dqn_vadv=0.0;
   du_vadv=0.0; dv_vadv=0.0;

   ! --- find out time in forcing data ---- !
   call interp_time(time_interp,time_forc,itime_less,itime_more,  &
                    weight_less,weight_more)

   !---- interpolate surface pressure field
   call interpolate_1d_field                                               &
        (weight_less,weight_more,ps_forc(itime_less),ps_forc(itime_more),  &
         tmp_ans(1))
   ps = tmp_ans(1)

   ! --- update pf, ph, and zf, zh
   ph(:,:,1)=0.;
       do j=1,size(ps,2)
         do i=1,size(ps,1)
           call get_eta_level(nlev, ps(i,j) , pf(i,j,:), ph(i,j,:))
         enddo
       enddo
            call error_mesg ('rico_forc_init ',   &
            'need to compute height', FATAL )
 
   !---- interpolate forcing data: dtdt, dqvdt, dudt, dvdt, u_geos, v_geos

   call interpolate_2d_field                                   &
        (weight_less, weight_more,                             &
         pres_forc(:,itime_less), dudt_forc(:,itime_less),     &
         pres_forc(:,itime_more), dudt_forc(:,itime_more),     &
         pf(1,1,:), du_lf(1,1,:))
   
   call interpolate_2d_field                                   &
        (weight_less, weight_more,                             &
         pres_forc(:,itime_less), dvdt_forc(:,itime_less),     &
         pres_forc(:,itime_more), dvdt_forc(:,itime_more),     &
         pf(1,1,:), dv_lf(1,1,:))
 
   call interpolate_2d_field                                   &
        (weight_less, weight_more,                             &
         pres_forc(:,itime_less), dtdt_forc(:,itime_less),     &
         pres_forc(:,itime_more), dtdt_forc(:,itime_more),     &
         pf(1,1,:), dT_lf(1,1,:))
   
   call interpolate_2d_field                                   &
        (weight_less, weight_more,                             &
         pres_forc(:,itime_less), dqvdt_forc(:,itime_less),    &
         pres_forc(:,itime_more), dqvdt_forc(:,itime_more),    &
         pf(1,1,:), dqv_lf(1,1,:))

    call interpolate_2d_field                                  &
        (weight_less, weight_more,                             &
         pres_forc(:,itime_less), u_geos_forc(:,itime_less),   &
         pres_forc(:,itime_more), u_geos_forc(:,itime_more),   &
         pf(1,1,:), u_geos(1,1,:))

     call interpolate_2d_field                                 &
        (weight_less, weight_more,                             &
         pres_forc(:,itime_less), v_geos_forc(:,itime_less),   &
         pres_forc(:,itime_more), v_geos_forc(:,itime_more),   &
         pf(1,1,:), v_geos(1,1,:))
   
     ! --- nudge winds
     du_nudge = 0.0
     dv_nudge = 0.0
     if (do_nudge) then
        do k=1, kdim
           fnudge = 0.5 * ( 1.0 - tanh((pf(1,1,k)-400.0e2)/15.0e2) )
           du_nudge(:,:,k) = fnudge * (u_geos(1,1,k)-ua(1,1,k))/3600.0
           dv_nudge(:,:,k) = fnudge * (v_geos(1,1,k)-va(1,1,k))/3600.0
        end do
     end if

end if   !  configuration

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

! --- sum all tendencies
u_dt = du_vadv + du_geos + du_lf + du_nudge
v_dt = dv_vadv + dv_geos + dv_lf + dv_nudge
t_dt = dT_rad + dT_vadv + dT_lf
q_dt(:,:,:,nsphum) = dqv_vadv + dqv_lf
q_dt(:,:,:,nql) = dql_vadv
q_dt(:,:,:,nqi) = dqi_vadv
q_dt(:,:,:,nqa) = dqa_vadv

if( nqn > 0) q_dt(:,:,:,nqn) = dqn_vadv

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

! ---> h1g, 2010-07-19
   qvdt_forcing_col = 0.0
   do k=1, kdim
         qvdt_forcing_col =  qvdt_forcing_col &
            + ( dqv_vadv( :,:,k ) + dqv_lf( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do

   qldt_vadv_col = 0.0
   do k=1, kdim
         qldt_vadv_col =  qldt_vadv_col &
            + (  dql_vadv( :,:,k )  ) * ( ph( :,:,k+1 ) -  ph( :,:,k ) ) / grav
   end do
! <--- h1g, 2010-07-19

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

! ---> h1g, 2010-07-19
if ( id_qvdt_forc_col > 0 )  used = send_data(  id_qvdt_forc_col, qvdt_forcing_col(:,:), time_diag, 1, 1 )
if ( id_qldt_vadv_col > 0 )  used = send_data(  id_qldt_vadv_col, qldt_vadv_col(:,:), time_diag, 1, 1 )
! <--- h1g, 2010-07-19

end subroutine update_rico_forc

!########################################################################
! Subroutine retuns RICO initial sounding

       subroutine rico_snd( configuration, z, u, v, ug, vg, T, qv )
       implicit none

       character(len=64), intent(in) :: configuration
       real, intent(in) :: z
       real, intent(out) :: u, v, ug, vg, T, qv

!      u wind

       if (trim(configuration)=='short') then
         if ( z < 4000.0 ) then
           u = -9.9 + (-1.9 + 9.9) / 4000.0 * z 
         else if ( z < 12000.0 ) then
           u = -1.9 + (30.0 + 1.9) / (12000.0 - 4000.0) * (z - 4000.0)
         else if ( z < 13000.0 ) then
           u = 30.0
         else if ( z < 20000.0 ) then
           u = 30.0 - (30.0) / (20000.0 - 13000.0) * (z - 13000.0)
         else
           u = 0.0
         end if
       else if (trim(configuration)=='long') then
         if ( z < 570.0 ) then
           u = -7.0 + (-8.85 + 7.0) / 570.0 * z 
         else if ( z < 6800.0 ) then
           u = -8.85 + (4.0 + 8.85) / (6800.0 - 570.0) * (z - 570.0)
         else if ( z < 13400.0 ) then
           u = 4.0 + (26.2 - 4.0) / (13400.0 - 6800.0) * (z - 6800.0)
         else if ( z < 18600.0 ) then
           u = 26.2 + (0.5 - 26.2) / (18600.0 - 13400.0) * (z - 13400.0)
         else if ( z < 30000.0 ) then
           u = 0.5 + (-6.0 - 0.5) / (30000.0 - 18600.0) * (z - 18600.0)
         else
           u = -6.0
         end if
       end if

!      v wind

       if (trim(configuration)=='short') then
         v = -3.8
       else if (trim(configuration)=='long') then
         if ( z < 9500.0 ) then
           v = -1.5 + (-3.5 + 1.5) / 9500.0 * z 
         else if ( z < 13500.0 ) then
           v = -3.5 + (-13.5 + 3.5) / (13500.0 - 9500.0) * (z - 9500.0)
         else if ( z < 16000.0 ) then
           v = -13.5 + (-5.0 + 13.5) / (16000.0 - 13500.0) * (z - 13500.0)
         else if ( z < 30000.0 ) then
           v = -5.0 + (3.0 + 5.0) / (30000.0 - 16000.0) * (z - 16000.0)
         else
           v = 3.0
         end if
       end if

!      qv
!h1g
       if ( z < 740.0 ) then
         qv = 16.0 + (13.8 - 16.0) / (740.0) * z
       else if ( z < 3260.0 ) then
         qv = 13.8 + (2.4 - 13.8) / (3260.0 - 740.0) * (z - 740.0)
       else if ( z < 4000.0 ) then
         qv = 2.4 + (1.8 - 2.4) / (4000.0 - 3260.0) * (z - 3260.0)
       else if ( z < 9000.0 ) then
         qv = 1.8 + (0. - 1.8) / (10000.0 - 4000.0) * (z - 4000.0)
       else
         qv = 0.0
       end if

!copy from CLUBB
     !  if ( z < 740.0 ) then
     !    qv = 16.26 + (13.993 - 16.26) / (740.0) * z
     !  else if ( z < 3260.0 ) then
     !    qv = 13.993 + (2.406 - 13.993) / (3260.0 - 740.0) * (z - 740.0)
     !  else if ( z < 4000.0 ) then
     !    qv = 2.406 + (1.803 - 2.406) / (4000.0 - 3260.0) * (z - 3260.0)
     !  else if ( z < 9000.0 ) then
    !     qv = 1.803 + (0 - 1.803) / (9000.0 - 4000.0) * (z - 4000.0)
     !  else
     !    qv = 0.0
    !   end if
!h1g

!      Temperature

       if ( z < 740.0 ) then
         T = 299.2 + (292.0 - 299.2) / (740.0) * z
       else if ( z < 4000.0 ) then
         T = 292.0 + (278.0 - 292.0) / (4000.0 - 740.0) * (z - 740.0)
       else if ( z < 15000.0 ) then
         T = 278.0 + (203.0 - 278.0) / (15000.0 - 4000.0) * (z - 4000.0)
       else if ( z < 17500.0 ) then
         T = 203.0 + (194.0 - 203.0) / (17500.0 - 15000.0) * (z - 15000.0)
       else if ( z < 20000.0 ) then
         T = 194.0 + (206.0 - 194.0) / (20000.0 - 17500.0) * (z - 17500.0)
       else if ( z < 60000.0 ) then
         T = 206.0 + (270.0 - 206.0) / (60000.0 - 20000.0) * (z - 20000.0)
       else
         T = 270.0
       end if

!      u geostrophic wind

       if ( z < 4000.0 ) then
         ug = -9.9 + (-1.9 + 9.9) / 4000.0 * z 
       else if ( z < 12000.0 ) then
         ug = -1.9 + (30.0 + 1.9) / (12000.0 - 4000.0) * (z - 4000.0)
       else if ( z < 13000.0 ) then
         ug = 30.0
       else if ( z < 20000.0 ) then
         ug = 30.0 - (30.0) / (20000.0 - 13000.0) * (z - 13000.0)
       else
         ug = 0.0
       end if

!      v geostrophic wind

       vg = -3.8

       return
       end subroutine rico_snd

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

end subroutine interp_time 

!########################################################################

subroutine interpolate_1d_field                                       &
            (weight_less, weight_more, field_in_less, field_in_more,  &
             field_out)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the 1D time interpolated
!      value of a field given at the two time (less and more).
!      The field is linearly interpolated in time.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       ------
!       INPUT:
!       ------
!
!       weight_less    real number indicating closeness of interpolated
!                      time to time 'less' (range: 0 to 1)
!       weight_more    real number indicating closeness of interpolated
!                      time to time 'more' (range: 0 to 1)
!       field_in_less  input field at time 'less' on pressure levels p_in_less
!       field_in_more  input field at time 'more' on pressure levels p_in_more
!
!       -------
!       OUTPUT:
!       -------
!
!       field_out     the field interpolated in pressure and time to the 
!                     desired levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

real, intent (in)  :: weight_less, weight_more
real, intent (in)  :: field_in_less, field_in_more
real, intent (out) :: field_out

! --- interpolate field in time

   field_out = weight_more*field_in_more + weight_less*field_in_less    

! 
!-----------------------------------------------------------------------

end subroutine interpolate_1d_field

!########################################################################

subroutine interpolate_2d_field                                   &
            (weight_less, weight_more,                            &
             p_in_less, field_in_less, p_in_more, field_in_more,  &
             p_out, field_out)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the 2D time/pressure interpolated
!      value of a field given at the two time (less and more).
!      The field is linearly interpolated in time and in pressure.
!
!      THE subroutine ASSUMES THAT p_in_(1) < p_in_(2) < p_in_(3)...
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       VARIABLES
!
!       ------
!       INPUT:
!       ------
!
!       weight_less    real number indicating closeness of interpolated
!                      time to time 'less' (range: 0 to 1)
!       weight_more    real number indicating closeness of interpolated
!                      time to time 'more' (range: 0 to 1)
!       p_in_less      pressure levels of input field at time 'less' (Pa)
!       field_in_less  input field at time 'less' on pressure levels p_in_less
!       p_in_more      pressure levels of input field at time 'more' (Pa)
!       field_in_more  input field at time 'more' on pressure levels p_in_more
!       p_out          pressure levels being interpolated to (Pa)
!
!       -------
!       OUTPUT:
!       -------
!
!       field_out     the field interpolated in pressure and time to the 
!                     desired levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

real, intent (in)                :: weight_less, weight_more
real, intent (in),  dimension(:) :: p_in_less, field_in_less, p_in_more, field_in_more
real, intent (in),  dimension(:) :: p_out
real, intent (out), dimension(:) :: field_out

!  Internal variables
!  ------------------

real,     dimension(size(p_out,1))   :: field_out_less, field_out_more
real,     dimension(size(p_out,1))   :: pweight_less, pweight_more
integer,  dimension(size(p_out,1))   :: klev_less, klev_more
logical,  dimension(size(p_out,1))   :: find_p_out

integer n_p_in_less, n_p_in_more, n_p_out

integer j

!
! Code
! ----

! ---  Array sizes

   n_p_in_less = size(p_in_less,1)
   n_p_in_more = size(p_in_more,1)
   n_p_out = size(p_out,1)
      
! ---  Interpolate 'less' field to desired pressure levels

   find_p_out(:) = .false.
   pweight_more(:) = 0.0
   pweight_less(:) = 0.0

   where ( p_out(:) .ge. p_in_less(n_p_in_less) )
      find_p_out(:) = .true.
      klev_more(:) = n_p_in_less
      klev_less(:) = n_p_in_less
      pweight_more(:) = 1.0
      pweight_less(:) = 0.0
   end where

   where ( p_out(:) .le. p_in_less(1) )
      find_p_out(:) = .true.
      klev_more(:) = 1
      klev_less(:) = 1
      pweight_more(:) = 1.0
      pweight_less(:) = 0.0
   end where

   do j = 2, n_p_in_less
      where ( p_in_less(j) .ge. p_out(:) .and. .not. find_p_out(:) )
         find_p_out(:) = .true.
         klev_less(:) = j-1
         klev_more(:) = j
         pweight_more(:) =  ( p_out(:)             - p_in_less(klev_less) )  &
                           /( p_in_less(klev_more) - p_in_less(klev_less) )
         pweight_less(:) = 1.0 - pweight_more(:)
      end where
   enddo

!  interpolate

   field_out_less(:) = pweight_more(:)*field_in_less(klev_more(:))   &
                     + pweight_less(:)*field_in_less(klev_less(:))
 
! ---  Interpolate 'more' field to desired pressure levels

   find_p_out(:) = .false.
   pweight_more(:) = 0.0
   pweight_less(:) = 0.0

   where ( p_out(:) .ge. p_in_more(n_p_in_more) )
      find_p_out(:) = .true.
      klev_more(:) = n_p_in_more
      klev_less(:) = n_p_in_more
      pweight_more(:) = 1.0
      pweight_less(:) = 0.0
   end where

   where ( p_out(:) .le. p_in_more(1) )
      find_p_out(:) = .true.
      klev_more(:) = 1
      klev_less(:) = 1
      pweight_more(:) = 1.0
      pweight_less(:) = 0.0
   end where

   do j = 2, n_p_in_more
      where ( p_in_more(j) .ge. p_out(:) .and. .not. find_p_out(:) )
         find_p_out(:) = .true.
         klev_less(:) = j-1
         klev_more(:) = j
         pweight_more(:) =  ( p_out(:)             - p_in_more(klev_less) )  &
                           /( p_in_more(klev_more) - p_in_more(klev_less) )
         pweight_less(:) = 1.0 - pweight_more(:)
      end where
   enddo

!  interpolate

   field_out_more(:) = pweight_more(:)*field_in_more(klev_more(:))   &
                     + pweight_less(:)*field_in_more(klev_less(:))
 

! --- interpolate field in time

   field_out = weight_more*field_out_more + weight_less*field_out_less    

! 
!-----------------------------------------------------------------------
 
end subroutine interpolate_2d_field 

!#########################################################
! This subroutine returns imposed surface fluxes
! copy from CLUBB
subroutine rico_sfclyr( um_sfc, vm_sfc, thlm_sfc, rtm_sfc, & 
                          lowestlevel, sst, psfc,    & 
                          upwp_sfc,    vpwp_sfc,   &
			  wpthlp_sfc,  wprtp_sfc)			 
			  
!----------------------------------------------------------------------
!        Description:
!          Surface forcing subroutine for RICO case.  Written
!          December 2006 by Michael Falk.
!
!          Updated to use specific formulations for surface fluxes
!          as specified in the RICO 3D LES specification, in hopes that
!          they'll be more accurate.
!
!        References:
!          ATEX: http://www.atmos.ucla.edu/~bstevens/gcss/setup.html
!          RICO: http://www.knmi.nl/samenw/rico/setup3d.html
implicit none

  ! Constants
  real, parameter :: & 
    ubmin   = 0.25,      & ! Minimum value for ubar.
!    ustar   = 0.3,       & ! Defined by ATEX specification
    C_10    = 0.0013,    & ! Drag coefficient, defined by ATEX specification
    C_m_20  = 0.001229,  & ! Drag coefficient, defined by RICO 3D specification
    C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
    C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
    z0      = 0.00015      ! Roughness length, defined by ATEX specification

  real :: & 
    Cz,   & ! This is C_10 scaled to the height of the lowest model level.
    Cm,   & ! This is C_m_20 scaled to the height of the lowest model level.
    Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
    Cq      ! This is C_q_20 scaled to the height of the lowest model level.

  logical :: & 
    l_use_old_atex  ! if true, use ATEX version; if not, use RICO-specific

  ! Input variables
  real, intent(in),  dimension(:)  :: & 
    um_sfc,        & ! This is u at the lowest above-ground model level.  [m/s]
    vm_sfc,        & ! This is v at the lowest above-ground model level.  [m/s]
    thlm_sfc,      & ! This is theta-l at the lowest above-ground model level.  
                     ! (DOES THIS NEED A CORRECTION FOR THETA-L TO THETA?)  [K]
    rtm_sfc            ! This is rt at the lowest above-ground model level.  [kg/kg]
    
 real, intent(in)  :: & 
    lowestlevel,   & ! This is z at the lowest above-ground model level.  [m]
    sst,               & ! This is the sea surface temperature [K].
    psfc             ! This is the surface pressure [Pa].

! ---> h1g, 10-09-07
     real   ::  qsat 
! <--- h1g, 10-09-07 

  ! Output variables
  real, intent(out) ,  dimension(:) ::  & 
    upwp_sfc,   & ! The upward flux of u-momentum         [(m^2 s^-2]
    vpwp_sfc,   & ! The Upward flux of v-momentum         [(m^2 s^-2]
    wpthlp_sfc, & ! The upward flux of theta-l            [K m s^-1]
    wprtp_sfc    ! The upward flux of rtm (total water)  [kg kg^-1 m s^-1]
    
  ! Internal variables
  real, dimension(size(um_sfc))  :: & 
    ubar   ! This is root (u^2 + v^2), per ATEX and RICO spec.

   real  :: ustar         ! surface friction velocity             [m/s]
  ! Declare the value of ustar.
  ustar = 0.3
! ------------------------------------------------------------------------------------------

  ! Choose which scheme to use
!h1g to use CLUBB in AM3p5
  l_use_old_atex = .FALSE.   
!h1g

  ! Define variable values
  ubar = max(ubmin, sqrt(um_sfc*um_sfc + vm_sfc*vm_sfc))
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cz   = C_10 * ((log(10./z0))/(log(lowestlevel/z0))) * & 
         ((log(10./z0))/(log(lowestlevel/z0)))         
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cm   = C_m_20 * ((log(20./z0))/(log(lowestlevel/z0))) * & 
         ((log(20./z0))/(log(lowestlevel/z0)))             
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Ch   = C_h_20 * ((log(20./z0))/(log(lowestlevel/z0))) * & 
         ((log(20./z0))/(log(lowestlevel/z0)))          
         ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cq   = C_q_20 * ((log(20./z0))/(log(lowestlevel/z0))) * & 
         ((log(20./z0))/(log(lowestlevel/z0)))            

! Compute heat and moisture fluxes
! qsat: saturated specific content  (  UNIT="kg(vapor) / kg (moist air)" )

  call compute_qs( sst, psfc, qsat )
  if (l_use_old_atex) then ! Use ATEX version
    wpthlp_sfc = -Cz * ubar * ( thlm_sfc - sst * (p00/psfc)**kappa ) ! [K m s^-1
    wprtp_sfc  = -Cz * ubar * ( rtm_sfc -  qsat ) ! [kg kg^-1  m s^-1]
    upwp_sfc   = -um_sfc * ustar*ustar / ubar                    ! [m^2 s^-2]
    vpwp_sfc   = -vm_sfc * ustar*ustar / ubar                    ! [m^2 s^-2]
 
  else ! Use RICO version
    wpthlp_sfc = -Ch * ubar * ( thlm_sfc - sst * (p00/psfc)**kappa ) ! K m s^-1
    wprtp_sfc  = -Cq * ubar * ( rtm_sfc - qsat ) ! [kg kg^-1  m s^-1]
    upwp_sfc   = -Cm * ubar * um_sfc    ! m^2 s^-2
    vpwp_sfc   = -Cm * ubar * vm_sfc    ! m^2 s^-2
  end if

  return
end subroutine rico_sfclyr
!#####################################



subroutine  rico_sfcflux_online (&
                 t_atm,     q_atm,      u_atm,       v_atm,            &
                 p_atm,     z_atm,      p_surf,      t_surf,           &
                 t_ca,      q_surf,     u_surf,      v_surf,           &
                 rough_mom, rough_heat, rough_moist, rough_scale, gust,&
                 flux_t,    flux_q,     flux_lw,     flux_u,           &
                 flux_v,    cd_m,       cd_t,        cd_q,             &
                 w_atm,     u_star,     b_star,      q_star,           &
                 dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
                 dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
                 dt,        land,       seawater,    avail  )  

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
  
  real, dimension(size(t_atm(:))) :: t_surf_tmp, q_surf_tmp,   p_surf_tmp
  real :: sst,     qsat,  psfc 
  !---------------------------------------------------------------------

          sst   = 299.8
          psfc  =  101540.
          call compute_qs( sst, psfc, qsat )

          t_surf_tmp     =  sst 
          q_surf_tmp  =  qsat
          p_surf_tmp  =   psfc
  
            call surface_flux (                                       &
                 t_atm,      q_atm,          u_atm,            v_atm,            &
                 p_atm,     z_atm,          p_surf_tmp,     t_surf_tmp,           &
                 t_ca,        q_surf_tmp,   u_surf,            v_surf,           &
                 rough_mom, rough_heat, rough_moist, rough_scale, gust,&  
                 flux_t,    flux_q,     flux_lw,     flux_u,           &
                 flux_v,    cd_m,       cd_t,        cd_q,             &
                 w_atm,     u_star,     b_star,      q_star,           &
                 dhdt_surf, dedt_surf,  dedq_surf,   drdt_surf,        &
                 dhdt_atm,  dedq_atm,   dtaudu_atm,  dtaudv_atm,       &
                 dt,        land,       seawater,    avail             )        

  return
end subroutine  rico_sfcflux_online
!#####################################

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

end module scm_rico_mod
