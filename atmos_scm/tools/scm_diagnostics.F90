
module scm_diagnostics_mod

!-----------------------------------------------------------------------

use            fms_mod, only: file_exist, open_namelist_file, error_mesg, &
                              NOTE, close_file, check_nml_error, mpp_pe,  &
                              mpp_root_pe,write_version_number, FATAL
use            mpp_mod, only: stdlog
use   time_manager_mod, only: time_type

use scm_horiz_grid_mod, only: horiz_grid_type
use   scm_prog_var_mod, only: prog_var_type

use   diag_manager_mod, only: diag_axis_init, register_diag_field, &
                              register_static_field, send_data
use      constants_mod, only: grav, cp_air, PI

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_names, get_tracer_index, &
                              get_number_tracers

use mpp_domains_mod, only: mpp_get_compute_domain
use            fv_pack, only: nlon, mlat, nlev,  beglat, endlat, beglon, &
                              endlon, lon, lat, lonb, latb, fv_domain, &
                              ak, bk, get_eta_level, nt_prog
 
implicit none
private

public :: scm_diagnostics_init, &
          scm_diagnostics


!-----------------------------------------------------------------------
!------------------------- axis names ----------------------------------

character(len=8) :: axiset = 'dynamics'
character(len=8) :: mod_name = 'dynamics'

!-----------------------------------------------------------------------

   integer, parameter :: mxtr = 10
!   real,    parameter :: ginv = 1./grav

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

integer :: id_hlonb, id_hlon , id_hlatb, id_hlat , &
           id_vlonb, id_vlon , id_vlatb, id_vlat , &
           id_phalf, id_pfull, id_hlat_wgt, id_vlat_wgt

integer :: id_bk   , id_pk   , id_zsurf, id_res  , id_wspd,          &
           id_ps   , id_ucomp, id_vcomp, id_temp ,                   &
           id_omega, id_div  , id_vor  , id_pgfx , id_pgfy            
 
integer, allocatable :: id_tracer(:), id_tracer_tend(:)         

integer :: id_zfull     , id_zhalf     , id_pfull_tv  ,  id_phalf_tv,   &
           id_dv_phys   , id_du_phys   , &
           id_dv_adv    , id_du_adv    , &
           id_dt_phys   , id_dq_phys, id_dt_damp, id_dq_damp,           &
           id_dt_adv    , id_dq_adv , id_dq_fill, id_dl_adv ,           &
           id_dl_phys   , id_di_adv , id_di_phys,                       &
           id_pw        , id_lwp    , id_iwp    , id_ttend ,  id_ztend, &
           id_advttend  , id_dampttend,  id_physttend,  id_qtend, &
           id_advqtend  , id_dampqtend,  id_fillqtend, &
           id_physqtend , id_advltend,  id_physltend, &
           id_advitend  , id_physitend

 real, save :: ginv
 type(time_type) :: fv_time
 real, allocatable, save :: phalf(:)
 save fv_time
!-----------------------------------------------------------------------

 character(len=128) :: version = '$Id$'
 character(len=128) :: tagname = '$Name$'

!-----------------------------------------------------------------------

contains

!#######################################################################


subroutine scm_diagnostics_init ( Time, fis, axes)
#include "fv_arrays.h"

   type(time_type),       intent(in)  :: Time
   real, dimension(:,:),  intent(in)  :: fis
   integer, dimension(4), intent(out) :: axes

   real, dimension(beglon:endlon+1) :: hlonb
   real, dimension(beglon:endlon+1) :: vlonb
   real, dimension(beglat:endlat+1) :: hlatb
   real, dimension(beglat:endlat+1) :: vlatb

   real, dimension(beglon:endlon) :: hlon
   real, dimension(beglon:endlon) :: vlon
   real, dimension(beglat:endlat) :: hlat
   real, dimension(beglat:endlat) :: vlat

   real, dimension(beglat:endlat) :: hlat_wgt, vlat_wgt

   real, dimension(nlev)   :: pfull

   integer :: id_lonb, id_lon, id_latb, id_lat, id_phalf, id_pfull
   real    :: vrange(2), trange(2), arange(2), slprange(2)
   integer :: i, j, n, ntrace, ntprog, logunit
   integer :: isg, ieg, jsg, jeg
   integer :: is, ie, js, je
   logical :: used
   character(len=128) :: tname
   character(len=256) :: longname, units
   integer, dimension(3) :: half = (/1,2,4/)

      real :: missing_value = -999.

#include "fv_point.inc"
!-----------------------------------------------------------------------

      call write_version_number (version,tagname)


   ! diagnostics for all tracers
   call get_number_tracers ( MODEL_ATMOS, num_tracers=ntrace, &
                                          num_prog=ntprog )
   logunit=stdlog()

!--------------------------- set up axes -------------------------------

100 format ('Diagnostics for the following tracer fields are available for module name = ',a)
110 format (3x,a,' (',a,'; ',a,')')


!-----------------------------------------------------------------------
!         -------- tendencies ---------

    id_du_phys = register_diag_field (mod_name,'udt_phys',    &
         axes(1:3), Time,                                &
         'zonal wind Tendency from Physics', 'm/s2',          &
         missing_value = missing_value )
    id_du_adv = register_diag_field (mod_name,'udt_adv',      &
         axes(1:3), Time,                                &
         'zonal wind Tendency from SCM forcing', 'm/s2',      &
         missing_value = missing_value )

    id_dv_phys = register_diag_field (mod_name,'vdt_phys',    &
         axes(1:3), Time,                                &
         'meridional wind Tendency from Physics', 'm/s2',     &
         missing_value = missing_value )
    id_dv_adv = register_diag_field (mod_name,'vdt_adv',      &
         axes(1:3), Time,                                &
         'meridional wind Tendency from SCM forcing', 'm/s2', &
         missing_value = missing_value )

  
   allocate (id_tracer_tend(ntprog))
 ! need to do this for just prognostic variables
   do n = 1, ntprog
     call get_tracer_names ( MODEL_ATMOS, n, tname, longname, units )
     tname    = trim(tname)   //'_dt_dyn'
     longname = trim(longname)//' tendency for dynamics'
     units    = trim(units)   //'/s'
     if (units == 'none') units = '1/sec'
     if (mpp_pe() == mpp_root_pe()) write(logunit,110) trim(tname),trim(longname),trim(units)
     id_tracer_tend(n) = register_diag_field ( mod_name, trim(tname),  &
                                 axes(1:3), Time, trim(longname), &
                                 trim(units), missing_value=-999.      )
   enddo

!-----------------------------------------------------------------------

!The above register_diag_field calls were present in bgrid_diagnostics_init


    id_zfull= register_diag_field (mod_name, 'zfull',          &
         axes(1:3), Time,'Geopotential Height on Full Levels',  'm', &
         missing_value = missing_value )
   
    id_pfull_tv= register_diag_field (mod_name, 'pfull_tv',       &
         axes(1:3), Time,'Time varying Full Level Pressures',  'Pa', &
         missing_value = missing_value )
   
    id_zhalf= register_diag_field (mod_name, 'zhalf',          &
         axes(half), Time,'Geopotential Height on Half Levels',  'm',&
         missing_value = missing_value )
   
    id_phalf_tv= register_diag_field (mod_name, 'phalf_tv',       &
         axes(half), Time,'Time varying Half Level Pressures',  'Pa', &
         missing_value = missing_value )
   
    id_dt_phys = register_diag_field (mod_name,'dt_phys',      &
         axes(1:3), Time,                                 &
         'Temperature Tendency from Physics', 'K/sec',         &
         missing_value = missing_value )
    id_dq_phys = register_diag_field (mod_name, 'dq_phys',     &
         axes(1:3), Time,                                 &
         'Water Vapor Tendency from Physics', 'kg/kg/sec',     &
         missing_value = missing_value )
    id_dt_damp = register_diag_field (mod_name, 'dt_damp',     &
         axes(1:3), Time,                                 &
         'Temperature Tendency from Damping', 'K/sec',         &
          missing_value = missing_value )
    id_dq_damp = register_diag_field (mod_name, 'dq_damp',     &
         axes(1:3), Time,                                 &
         'Water Vapor Tendency from Damping', 'kg/kg/sec',     &
         missing_value = missing_value )
    id_dt_adv  = register_diag_field (mod_name, 'dt_adv',      &
         axes(1:3), Time,                                 &
         'Temperature Tendency from Advection', 'K/sec',       &
         missing_value = missing_value )
    id_dq_adv  = register_diag_field (mod_name, 'dq_adv',      &
         axes(1:3), Time,                                 &
         'Water Vapor Tendency from Advection', 'kg/kg/sec',   &
         missing_value = missing_value )
    id_dq_fill = register_diag_field (mod_name, 'dq_fill',     &
         axes(1:3), Time,                                 &
         'Water Vapor Tendency from Filling', 'kg/kg/sec',     &
         missing_value = missing_value )
    id_dl_adv  = register_diag_field (mod_name, 'dl_adv',      &
         axes(1:3), Time,                                 &
         'Liquid Water Tendency from Advection', 'kg/kg/sec',  &
         missing_value = missing_value )
    id_dl_phys = register_diag_field (mod_name, 'dl_phys',     &
         axes(1:3), Time,                                 &
         'Liquid Water Tendency from Physics', 'kg/kg/sec',    &
         missing_value = missing_value )
    id_di_adv  = register_diag_field (mod_name, 'di_adv',      &
         axes(1:3), Time,                                 &
         'Ice Water Tendency from Advection', 'kg/kg/sec',     &
         missing_value = missing_value )
    id_di_phys = register_diag_field (mod_name, 'di_phys',     &
         axes(1:3), Time,                                 &
         'Ice Water Tendency from Physics', 'kg/kg/sec',       &
         missing_value = missing_value )
   
    id_pw =    register_diag_field (mod_name, 'PW',            &
         axes(1:2), Time, 'Precipitable Water', 'kg/m2',  &
         missing_value = missing_value )
    id_lwp =   register_diag_field (mod_name, 'LWP',           &
         axes(1:2), Time, 'Liquid Water Path', 'kg/m2',   &
         missing_value = missing_value )
    id_iwp =   register_diag_field (mod_name, 'IWP',           &
         axes(1:2), Time, 'Ice Water Path', 'kg/m2',      &
         missing_value = missing_value )
    id_ttend = register_diag_field (mod_name, 'ttend',         &
         axes(1:2), Time, 'Column Energy Tendency',       &
         'W/m2', missing_value = missing_value )
    id_ztend = register_diag_field (mod_name, 'ztend',         &
         axes(1:2), Time,                                 &
         'Column Geopotential Energy Tendency', 'W/m2',        &
         missing_value = missing_value )
    id_advttend = register_diag_field (mod_name, 'advttend',   &
         axes(1:2), Time,                                 &
         'Column Advective Energy Tendency', 'W/m2',           &
         missing_value = missing_value )
    id_dampttend = register_diag_field (mod_name, 'dampttend', &
         axes(1:2), Time, 'Column Damping Energy Tendency', &
         'W/m2', missing_value = missing_value )
    id_physttend = register_diag_field (mod_name, 'physttend', &
         axes(1:2), Time, 'Column Physics Energy Tendency', &
         'W/m2', missing_value = missing_value )  
    id_qtend = register_diag_field (mod_name, 'qtend',         &
         axes(1:2), Time, 'Column Moisture Tendency',     &
         'kg/m2/sec', missing_value = missing_value )
    id_advqtend = register_diag_field (mod_name, 'advqtend',   &
         axes(1:2), Time,                                 &
         'Column Moisture Advection Tendency', 'kg/m2/sec',    &
         missing_value = missing_value )
    id_dampqtend = register_diag_field (mod_name, 'dampqtend', &
         axes(1:2), Time,                                 &
         'Column Moisture Damping Tendency', 'kg/m2/sec',      &
         missing_value = missing_value )
    id_fillqtend = register_diag_field (mod_name, 'fillqtend', &
         axes(1:2), Time,                                 &
         'Column Moisture Filling Tendency', 'kg/m2/sec',      &
         missing_value = missing_value )
    id_physqtend = register_diag_field (mod_name, 'physqtend', &
         axes(1:2), Time,                                 &
         'Column Moisture Physics Tendency', 'kg/m2/sec',      &
         missing_value = missing_value )
    id_advltend = register_diag_field (mod_name, 'advltend',   &
         axes(1:2), Time,                                 &
         'Column Liquid Water Advection Tendency', 'kg/m2/sec', &
         missing_value = missing_value )
    id_physltend = register_diag_field (mod_name, 'physltend', &
         axes(1:2), Time,                                 &
         'Column Liquid Water Physics Tendency', 'kg/m2/sec',  &
         missing_value = missing_value )
    id_advitend = register_diag_field (mod_name, 'advitend',   &
         axes(1:2), Time,                                 &
         'Column Ice Water Advection Tendency', 'kg/m2/sec',   &
         missing_value = missing_value )
    id_physitend = register_diag_field (mod_name, 'physitend', &
         axes(1:2), Time,                                 &
         'Column Ice Water Physics Tendency', 'kg/m2/sec',     &
         missing_value = missing_value )
   


   end subroutine scm_diagnostics_init

!#######################################################################

 subroutine scm_diagnostics ( phalf, Forc_tend, Phys_tend, pfull,      &
                              zfull, zhalf, damp_t, damp_q,            &
                              fill_q, Time)
#include "fv_arrays.h"

!-----------------------------------------------------------------------
   type (prog_var_type),  intent(in) :: Forc_tend, Phys_tend
   type(time_type),       intent(in) :: Time
   real, intent(in), dimension(beglon:endlon,beglat:endlat,nlev+1) :: phalf, zhalf
   real, intent(in), dimension(beglon:endlon,beglat:endlat,nlev) :: fill_q, damp_t, &
                                                                    damp_q, pfull,  &
                                                                    zfull

                                                       
   real, dimension(beglon:endlon,beglat:endlat,nlev) :: wspd, vor, div, adp, udp, vdp
   real, dimension(beglon:endlon,beglat:endlat) :: ttend,  &
         advttend, dampttend, physttend, qtend, advqtend,   &
         dampqtend, physqtend, physltend, physitend,        &
         fillqtend, advltend, advitend

   real, dimension(size(pt,1),size(pt,2))   :: PW,LWP,IWP
#include "fv_point.inc"


!-----------------------------------------------------------------------
   integer :: is, ie, hs, he, vs, ve, n, k
   integer :: sphum, nql, nqi
   logical :: used
!-----------------------------------------------------------------------

      call mpp_get_compute_domain ( fv_domain, is, ie, hs, he )
      vs = hs;  ve = he
  
  !-----------------------------------------------------------------------
!   get tracer indices

    sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    if (sphum <= 0) call error_mesg ('scm_diagnostics_mod', &
         'specific humidity tracer not found', FATAL)
  
    nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
    if (nql <= 0) call error_mesg ('scm_diagnostics_mod', &
         'liquid water tracer not found', FATAL)
    nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
    if (nqi <= 0) call error_mesg ('scm_diagnostics_mod', &
         'ice water tracer not found', FATAL)
    if (min(nql,nqi) <= 0) call error_mesg ('scm_diagnostics_mod', &
         'stratiform cloud tracer(s) not found', FATAL)
    if (nql == nqi)  call error_mesg ('scm_diagnostics_mod',  &
         'tracers indices cannot be the same (i.e., nql=nqi).', FATAL)
   
!-----------------------------------------------------------------------
!---------------- surface fields ---------------------------------------

! Already passed out in fv_diag
!      if ( id_ps > 0 ) &
!      used = send_data ( id_ps , ps(is:ie,hs:he), Time )

!---------------- 3d momentum fields (u & v) ---------------------------

! Already passed out in fv_diag
!      if ( id_ucomp > 0 ) &
!      used = send_data ( id_ucomp, ua(is:ie,vs:ve,:), Time )
!
!      if ( id_vcomp > 0 ) &
!      used = send_data ( id_vcomp, va(is:ie,vs:ve,:), Time )

      if ( id_temp > 0 ) &
      used = send_data ( id_temp, pt(is:ie,hs:he,:), Time)

!      do n = 1, ncnst
!        if ( id_tracer(n) > 0 ) &
!        used = send_data ( id_tracer(n), q(is:ie,hs:he,:,n), Time )
!      enddo
!
!      if ( id_omega > 0 ) &
!      used = send_data ( id_omega, omga(is:ie,hs:he,:), Time )

!------------ wind speed, divergence, vorticity ------------------------

      if ( id_wspd > 0 ) then
          wspd(is:ie,vs:ve,:) = sqrt &
                      ( ua(is:ie,vs:ve,:)*ua(is:ie,vs:ve,:) + &
                        va(is:ie,vs:ve,:)*va(is:ie,vs:ve,:) )
          used = send_data ( id_wspd, wspd(is:ie,vs:ve,:), Time)
      endif

!------------ pressure and heights of model levels ---------------------

    if (id_zfull > 0) then
        used = send_data( id_zfull, zfull, Time, &
                          1, 1, 1)
    end if 
    
    if (id_pfull_tv > 0) then
        used = send_data( id_pfull_tv, pfull, Time, &
                          1, 1, 1)
    end if 
    
    if (id_zhalf > 0) then
        used = send_data( id_zhalf, zhalf, Time, &
                          1, 1, 1)
    end if 
    
    if (id_phalf_tv > 0) then
        used = send_data( id_phalf_tv, phalf, Time, &
                          1, 1, 1)
    end if 
    
!--------------- optional arguments ------------------------------------
!--------------- pressure gradient components --------------------------



    !-----------------
    ! Zonal wind
    !

    if (id_du_phys > 0) then
         used = send_data( id_du_phys, Phys_tend%u, Time, &
                           1, 1, 1)
    end if 
    
    if (id_du_adv > 0) then
         used = send_data( id_du_adv, Forc_tend%u, Time, &
                           1, 1, 1)
    end if 


    !-----------------
    ! Meridional wind
    !

    if (id_dv_phys > 0) then
         used = send_data( id_dv_phys, Phys_tend%v, Time, &
                           1, 1, 1)
    end if 
    
    if (id_dv_adv > 0) then
         used = send_data( id_dv_adv, Forc_tend%v, Time, &
                           1, 1, 1)
    end if 

    !-----------------
    ! Temperature
    !

    if (id_dt_phys > 0) then
         used = send_data( id_dt_phys, Phys_tend%T, Time, &
                           1, 1, 1)
    end if 
    
    if (id_dt_adv > 0) then
         used = send_data( id_dt_adv, Forc_Tend%T, Time, &
                           1, 1, 1)
    end if 
    
    if (id_dt_damp > 0) then
         used = send_data( id_dt_damp, damp_t, Time, &
                           1, 1, 1)
    end if
    
    !-----------------
    ! Zonal wind
    !

    if (id_dq_phys > 0) then
         used = send_data( id_dq_phys, Phys_tend%R(:,:,:,sphum), Time, &
                           1, 1, 1)
    end if 
   
    if (id_dq_adv > 0) then
         used = send_data( id_dq_adv, Forc_tend%R(:,:,:,sphum), Time, &
                           1, 1, 1)
    end if 
    
    if (id_dq_fill > 0) then
         used = send_data( id_dq_fill, fill_q, Time, &
                           1, 1, 1)
    end if 
    
    if (id_dq_damp > 0) then
         used = send_data( id_dq_damp, damp_q, Time, &
                           1, 1, 1)
    end if 
    
    !----------------
    ! Cloud variables
    !
    
    if (id_dl_phys > 0) then
         used = send_data( id_dl_phys, Phys_tend%R(:,:,:,nql), Time, &
                           1, 1, 1)
    end if
    
    if (id_dl_adv > 0) then
         used = send_data( id_dl_adv, Forc_tend%R(:,:,:,nql), Time, &
                           1, 1, 1)
    end if 
     
    if (id_di_phys > 0) then
         used = send_data( id_di_phys, Phys_tend%R(:,:,:,nqi), Time, &
                           1, 1, 1)
    end if
    
    if (id_di_adv > 0) then
         used = send_data( id_di_adv, Forc_tend%R(:,:,:,nqi), Time, &
                           1, 1, 1)
    end if
 
 
      do n = 1, nt_prog
       if ( id_tracer_tend(n) > 0 ) &
       used = send_data ( id_tracer_tend(n), q(is:ie,hs:he,:,n), Time )
      enddo
    
    !compute column integrated prognostic variables

    PW=0.
    LWP=0.
    IWP=0.
     
    !loop over levels

    do k=1, size(q,3)
          PW =  PW + (q(:,:,k,sphum) *  &  
              (phalf(:,:,k+1)-phalf(:,:,k))/grav)
         LWP = LWP + (q(:,:,k,nql)   *  &
              (phalf(:,:,k+1)-phalf(:,:,k))/grav)
         IWP = IWP + (q(:,:,k,nqi)   *  &
              (phalf(:,:,k+1)-phalf(:,:,k))/grav)
    enddo 

    !update output fields
    if (id_pw > 0) then
         used = send_data( id_pw, PW, Time, 1, 1)
    end if  
    if (id_lwp > 0) then
         used = send_data( id_lwp, LWP, Time, 1, 1)
    end if  
    if (id_iwp > 0) then
         used = send_data( id_iwp, IWP, Time, 1, 1)
    end if 

    !----------------------------------------------
    !compute column integrated tendency budgets
   

    ttend = 0.; advttend = 0.; dampttend = 0.; physttend = 0.
    qtend = 0.; advqtend = 0.; dampqtend = 0.; physqtend = 0.; fillqtend = 0.
    advltend = 0.; physltend = 0.; advitend = 0.; physitend = 0.
   
    !loop over levels
    do k=1, size(pt,3)
        ttend =     ttend  +  (    t_dt(:,:,k) *  cp_air * &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
        qtend =     qtend  +  (    q_dt(:,:,k,sphum) *   &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
     advttend =  advttend  +  ( Forc_tend%T(:,:,k) *  cp_air * &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
    dampttend = dampttend  +  (      damp_t(:,:,k) *  cp_air * &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
    physttend = physttend  +  ( Phys_tend%T(:,:,k) *  cp_air * &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )

     advqtend =  advqtend  +  ( Forc_tend%R(:,:,k,sphum) *   &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
    dampqtend = dampqtend  +  (      damp_q(:,:,k)       *   &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
    fillqtend = fillqtend  +  (      fill_q(:,:,k)       *   &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
    physqtend = physqtend  +  ( Phys_tend%R(:,:,k,sphum) *   &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav )
     
     advltend =  advltend  +  ( Forc_tend%R(:,:,k,nql) *     &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav)
    physltend = physltend  +  ( Phys_tend%R(:,:,k,nql) *     &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav)
     advitend =  advitend  +  ( Forc_tend%R(:,:,k,nqi) *     &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav)
    physitend = physitend  +  ( Phys_tend%R(:,:,k,nqi) *     &
                         (phalf(:,:,k+1)-phalf(:,:,k)) / grav)
                                     
    enddo

    !update output fields   
    if (id_ttend > 0) then
         used = send_data( id_ttend, ttend, Time, 1, 1)
    end if
    if (id_advttend > 0) then
         used = send_data( id_advttend, advttend, Time, 1, 1)
    end if
    if (id_dampttend > 0) then
         used = send_data( id_dampttend, dampttend, Time, 1, 1)
    end if
    if (id_physttend > 0) then
         used = send_data( id_physttend, physttend, Time, 1, 1)
    end if
    if (id_qtend > 0) then
         used = send_data( id_qtend, qtend, Time, 1, 1)
    end if
    if (id_advqtend > 0) then
         used = send_data( id_advqtend, advqtend, Time, 1, 1)
    end if
    if (id_dampqtend > 0) then
         used = send_data( id_dampqtend, dampqtend, Time, 1, 1)
    end if
    if (id_fillqtend > 0) then
         used = send_data( id_fillqtend, fillqtend, Time, 1, 1)
    end if
    if (id_physqtend > 0) then
         used = send_data( id_physqtend, physqtend, Time, 1, 1)
    end if
    if (id_advltend > 0) then
         used = send_data( id_advltend, advltend, Time, 1, 1)
    end if
    if (id_physltend > 0) then
         used = send_data( id_physltend, physltend, Time, 1, 1)
    end if
    if (id_advitend > 0) then
         used = send_data( id_advitend, advitend, Time, 1, 1)
    end if
    if (id_physitend > 0) then
         used = send_data( id_physitend, physitend, Time, 1, 1)
    end if




!-----------------------------------------------------------------------

 end subroutine scm_diagnostics


!#######################################################################

end module scm_diagnostics_mod

