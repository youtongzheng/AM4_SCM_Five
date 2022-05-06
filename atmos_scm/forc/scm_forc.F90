module scm_forc_mod
   
   use            mpp_mod, only:  stdlog
   use            fms_mod, only:  check_nml_error,                      &
                                  mpp_pe, mpp_root_pe,                  &
                                  write_version_number,                 &
                                  open_file, close_file, file_exist,    &
                                  open_restart_file,                    &
                                  read_data, write_data,                &
                                  error_mesg, FATAL,                    &
                                  mpp_error, NOTE
#ifdef INTERNAL_FILE_NML
   use              mpp_mod, only: input_nml_file
#else
   use              fms_mod, only: open_namelist_file
#endif
   use   time_manager_mod

   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, grav

   use scm_horiz_grid_mod, only:  horiz_grid_type

   use   surface_flux_mod, only:  surface_flux

   use       scm_AMIE_mod, only:  AMIE_data_read, AMIE_forc_init, AMIE_forc_end,    &
                                  update_AMIE_forc, AMIE_forc_diagnostic_init,      &
                                  AMIE_surface_flux_loop, get_AMIE_sfc, get_AMIE_sst

   use        scm_arm_mod, only:  arm_data_read, arm_forc_init, arm_forc_end,        &
                                  update_arm_forc, arm_forc_diagnostic_init,         &
                                  arm_surface_flux_loop, get_arm_sfc

   use    scm_arm2003_mod, only:  arm2003_data_read, arm2003_forc_init, arm2003_forc_end, &
                                  update_arm2003_forc, arm2003_forc_diagnostic_init,      &
                                  arm2003_surface_flux_loop, get_arm2003_sfc

   use      scm_astex_mod, only:  astex_data_read, astex_forc_init, astex_forc_end,  &
                                  update_astex_forc, astex_forc_diagnostic_init,     &
                                  get_astex_sst

   use      scm_bomex_mod, only:  bomex_data_read, bomex_forc_init, bomex_forc_end,  &
                                  update_bomex_forc, bomex_forc_diagnostic_init,     &
                                  get_bomex_flx

   use       scm_dcbl_mod, only:  dcbl_data_read, dcbl_forc_init, dcbl_forc_end,  &
                                  update_dcbl_forc, dcbl_forc_diagnostic_init,     &
                                  get_dcbl_flx

   use     scm_mpacea_mod, only:  mpacea_data_read, mpacea_forc_init,                &
                                  mpacea_forc_end,                                   &
                                  update_mpacea_forc, mpacea_forc_diagnostic_init,   &
                                  mpacea_surface_flux_loop, get_mpacea_sfc

   use     scm_mpaceb_mod, only:  mpaceb_data_read, mpaceb_forc_init,                &
                                  mpaceb_forc_end,                                   &
                                  update_mpaceb_forc, mpaceb_forc_diagnostic_init,   &
                                  mpaceb_surface_flux_loop,                          &
                                  get_mpaceb_flx, get_mpaceb_sst

   use       scm_rf01_mod, only:  rf01_data_read, rf01_forc_init, rf01_forc_end,     &
                                  update_rf01_forc, rf01_forc_diagnostic_init,       &
                                  add_rf01_tdtlw, add_rf01_tdtsw, get_rf01_flx, get_rf01_flx_online

   use       scm_rf02_mod, only:  rf02_data_read, rf02_forc_init, rf02_forc_end,     &
                                  update_rf02_forc, rf02_forc_diagnostic_init,       &
                                  add_rf02_tdtlw, add_rf02_tdtsw, get_rf02_flx

   use       scm_rico_mod, only:  rico_data_read, rico_forc_init, rico_forc_end,     &
                                  update_rico_forc, rico_forc_diagnostic_init, rico_sfclyr,      & 
                                  rico_sfcflux_online  !h1g

   use       scm_atex_mod, only:  atex_data_read, atex_forc_init, atex_forc_end, &   !h1g
                                  update_atex_forc,  atex_forc_diagnostic_init,  &   !h1g
                                  add_atex_tdtlw, add_atex_tdtsw, atex_sfclyr        !h1g

   use       scm_gcss_arm_mod, only:  gcss_arm_data_read,     gcss_arm_forc_init,  &            !h1g
                                      update_gcss_arm_forc,   gcss_arm_forc_diagnostic_init,  & !h1g
                                      get_gcss_arm_flx,       gcss_arm_forc_end, &              !h1g
                                      gcss_arm_surface_flux_loop    ! ZNT 05/20/2020

    use      scm_sheba_mod, only:    sheba_data_read,      sheba_forc_init,   &   !h1g
                                     update_sheba_forc,    sheba_forc_diagnostic_init,   &   !h1g
                                     get_sheba_flx,        sheba_forc_end,   get_sheba_sst

   use      gcss_astex_mod, only:   gcss_astex_data_read, gcss_astex_forc_init,   &
                                     gcss_astex_forc_end,  update_gcss_astex_forc, &
                                     gcss_astex_forc_diagnostic_init,              &
                                     add_gcss_astex_tdtlw, add_gcss_astex_tdtsw,   &
                                     get_gcss_astex_flx

    use      scm_vocals_mod, only:   vocals_data_read, vocals_forc_init,           &
                                     vocals_forc_end, update_vocals_forc,          &
                                     vocals_forc_diagnostic_init,                  &
                                     get_vocals_flx


    use      scm_Vocals_75W_20S_mod, only:   Vocals_75W_20S_data_read, Vocals_75W_20S_forc_init,           &
                                     Vocals_75W_20S_forc_end, update_Vocals_75W_20S_forc,          &
                                     Vocals_75W_20S_forc_diagnostic_init,                  &
                                     get_Vocals_75W_20S_flx


    use      scm_Hawaii_143W_8N_mod, only:   Hawaii_143W_8N_data_read, Hawaii_143W_8N_forc_init,           &
                                     Hawaii_143W_8N_forc_end, update_Hawaii_143W_8N_forc,          &
                                     Hawaii_143W_8N_forc_diagnostic_init,                  &
                                     get_Hawaii_143W_8N_flx


! copy from Lin Yanluan
   use       scm_twp_mod, only:  twp_data_read, twp_forc_init, twp_forc_end,        &
                                 update_twp_forc, twp_forc_diagnostic_init,         &
                                 get_twp_sst, twp_surface_flux_loop, get_twp_sfc

   use      scm_gate_mod, only:  gate_data_read, gate_forc_init, gate_forc_end,        &
                                 update_gate_forc, gate_forc_diagnostic_init,         &
                                 get_gate_sst, gate_surface_flux_loop, get_gate_sfc

   use      scm_mc3e_mod, only:  mc3e_data_read, mc3e_forc_init, mc3e_forc_end,        &
                                 update_mc3e_forc, mc3e_forc_diagnostic_init,         &
                                 get_mc3e_sst, mc3e_surface_flux_loop, get_mc3e_sfc
   !yzheng: modules and variables related to FIVE
   use      five_mod, only: update_bomex_forc_five, update_rf01_forc_five

implicit none

public scm_data_read, scm_forc_init, scm_forc_end, update_scm_forc,  &
       scm_forc_diagnostic_init,                                     &
       use_scm_rad, add_scm_tdtlw, add_scm_tdtsw,                    &
       do_specified_flux, scm_surface_flux,                          &
       do_specified_tskin, TSKIN,  QSKIN,                                   &
       do_specified_albedo, ALBEDO_OBS,                              &
       do_specified_rough_leng, ROUGH_MOM, ROUGH_HEAT,               &
       do_specified_land

logical, private  :: initialized = .false.
character(len=8)  :: mod_name = 'forcing'
character(len=64) :: experiment = 'rico'

!
!       use_scm_rad               Use the SCM forcings radiative forcings?
!
!       do_specified_flux         Should surface fluxes be overwritten?
!
!       do_specified_tskin        Should surface skin temperature be 
!                                 specified from observations?
!
!       do_specified_albedo       Should the surface albedo be specified
!                                 from observations?
!
!       do_specified_rough_leng   Should roughness lengths be specified
!                                 from observations? 
!
!       do_specified_land         Should land conditions be specified?
!

logical           :: use_scm_rad             = .false.
logical           :: do_specified_flux       = .false.
logical           :: do_specified_tskin      = .false.
logical           :: do_specified_albedo     = .false.
logical           :: do_specified_rough_leng = .false.
logical           :: do_specified_land       = .false.
logical           :: use_rico_online_sfcflux = .false.
logical           :: use_gcss_arm_sfcflx_old = .false.    ! ZNT 05/20/2020

real              :: TSKIN, QSKIN
real              :: SENFLUX, EVAPFLUX
real              :: ALBEDO_OBS
real              :: ROUGH_MOM, ROUGH_HEAT

real, parameter   :: d622 = rdgas/rvgas
real, parameter   :: d378 = 1.0-d622
real, parameter   :: d608 = d378/d622

!--------------------- version number ----------------------------------
!
character(len=128) :: Version = '$Id$'
character(len=128) :: Tagname = '$Name$'
        
namelist /scm_forc_nml/ experiment, use_scm_rad,                      &
                        do_specified_flux,   do_specified_tskin,      &
                        do_specified_albedo, do_specified_rough_leng, &
                        do_specified_land,   use_rico_online_sfcflux, &
                        use_gcss_arm_sfcflx_old    ! ZNT 05/20/2020

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine scm_data_read(kmax)
implicit none

integer, intent (in) :: kmax

integer              :: unit,ierr,io, logunit
     
character(len=64)    :: fname_nc='INPUT/scm_forc.res.nc'

   if (initialized) return
   initialized = .true.
      
!-------- read namelist --------
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=scm_forc_nml, iostat=io)
   ierr = check_nml_error(io, 'scm_forc_nml')
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=scm_forc_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'scm_forc_nml')
      enddo
10    call close_file (unit)
   endif
#endif

!--------- write version number and namelist --------

   call write_version_number ( version, tagname )
   if(mpp_pe() == mpp_root_pe() ) then
     logunit =stdlog()
     write(logunit,nml=scm_forc_nml)
   endif

!--------- case specific initialization --------

   select case (trim(experiment))

      case ('arm')
         call arm_data_read()

      case ('arm2003')
         call arm2003_data_read()

      case ('AMIE')
         call AMIE_data_read()

! copy from Lin Yanluan, 2011-09-23
      case ('twp')
         call twp_data_read()

      case ('gate')
         call gate_data_read()

      case ('mc3e')
         call mc3e_data_read()

      case ('astex')
         call astex_data_read()

      case ('bomex')
         call bomex_data_read(kmax)

      case ('dcbl')
         call dcbl_data_read(kmax)

      case ('mpacea')
         call mpacea_data_read()

      case ('mpaceb')
          call mpaceb_data_read()
 
      case ('rf01')
         call rf01_data_read(kmax)

      case ('rf02')
         call rf02_data_read(kmax)

      case ('rico')
         call rico_data_read(kmax)

      case ('atex')
         call atex_data_read(kmax)

! ---> h1g, 2010-06-21
       case('gcss_arm')
          call  gcss_arm_data_read(kmax)
! <--- h1g, 2010-06-21

! ---> h1g, 2011-01-18
      case('sheba')
         call  sheba_data_read(kmax)
! <--- h1g, 2011-01-18

! ---> h1g, 2011-02-25
       case('gcss_astex')
          call  gcss_astex_data_read(kmax)
! <--- h1g, 2011-02-25


! ---> h1g, 2011-12-29
      case('vocals')
         call  vocals_data_read( )
! <--- h1g, 2011-12-29

! ---> h1g, 2012-2-27
      case('Vocals_75W_20S')
         call  Vocals_75W_20S_data_read( )
! <--- h1g, 2012-2-27

! ---> h1g, 2012-3-2
      case('Hawaii_143W_8N')
         call  Hawaii_143W_8N_data_read( )
! <--- h1g, 2012-3-2

    case  default
       call mpp_error ('scm_data_read',                          &
                       'Data ingestion for '//trim(experiment)//' not supported.', &
                       FATAL )

   end select

!--------- restart file --------

    ! If restart file exists override values with the restart values

    if( file_exist(trim(fname_nc)) ) then

          if(mpp_pe() == mpp_root_pe() )                            &
          call mpp_error ('scm_forc_mod',                           &
                          'Reading netCDF formatted restart file.', &
                           NOTE )

          call read_data(fname_nc, 'SENFLUX',    SENFLUX   )
          call read_data(fname_nc, 'EVAPFLUX',   EVAPFLUX  )
          call read_data(fname_nc, 'TSKIN',      TSKIN     )
          call read_data(fname_nc, 'ALBEDO_OBS', ALBEDO_OBS)
          call read_data(fname_nc, 'ROUGH_MOM',  ROUGH_MOM )
          call read_data(fname_nc, 'ROUGH_HEAT', ROUGH_HEAT)

    end if

end subroutine scm_data_read

!#######################################################################

subroutine scm_forc_init(time_interp,pdamp,elev,As,Bs)
implicit none

type(time_type)                          :: time_interp
real,  intent (in)                       :: pdamp
real,  intent (in), dimension(:)         :: As,Bs
real,  intent (inout), dimension(:,:)    :: elev

! As and Bs and Vgrid are passed here for historical reasons.
! Vgrid pressure grid is computed differently from the grid with As and Bs.

select case (trim(experiment))

   case ('arm')  !h1g
      call arm_forc_init(time_interp,pdamp,As,Bs)!h1g
      call get_arm_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g


   case ('arm2003')  !h1g
      call arm2003_forc_init(time_interp,pdamp,As,Bs)!h1g
      call get_arm2003_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g



   case ('AMIE')
     call AMIE_forc_init(time_interp,pdamp, As, Bs)
!     if ( do_specified_tskin ) call get_AMIE_sst( TSKIN )
     call get_amie_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g

   case ('twp')
     call twp_forc_init(time_interp, pdamp, As, Bs)
!     if ( do_specified_tskin ) call get_twp_sst( TSKIN )
     call get_twp_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g

   case ('gate')
     call gate_forc_init(time_interp,pdamp, As, Bs)
!     if ( do_specified_tskin ) call get_gate_sst( TSKIN )
     call get_gate_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g

   case ('mc3e')
     call mc3e_forc_init(time_interp,pdamp, As, Bs)
!     if ( do_specified_tskin ) call get_mc3e_sst( TSKIN )
     call get_mc3e_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g


   case ('astex')!h1g
      call astex_forc_init(time_interp,pdamp,As,Bs)!h1g
       if ( do_specified_tskin ) call get_astex_sst( TSKIN )!h1g

   case ('bomex')!h1g
      call bomex_forc_init(time_interp,As,Bs,elev)

   case ('dcbl')!h1g
      call dcbl_forc_init(time_interp,As,Bs,elev)

   case ('mpacea')!h1g
#ifdef BGRID
     if( present(qn) ) then!h1g
       if( present(qni) ) then!h1g
          call mpacea_forc_init(time_interp,pdamp,As,Bs,ps,T,u,v,qv,ql,qi,qa,qn, qni )!h1g
       else
          call mpacea_forc_init(time_interp,pdamp,As,Bs,ps,T,u,v,qv,ql,qi,qa,qn)!h1g
        endif !  present(qni)
     else!h1g
      call mpacea_forc_init(time_interp,pdamp,As,Bs,ps,T,u,v,qv,ql,qi,qa)!h1g
     endif!h1g
#endif
      call mpacea_forc_init(time_interp,pdamp,As,Bs)
      call get_mpacea_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )!h1g

   case ('mpaceb')!h1g
      call mpaceb_forc_init(time_interp,pdamp,As,Bs)
      if ( do_specified_tskin ) call get_mpaceb_sst( TSKIN )!h1g
      QSKIN=  1.95e-3

   case ('rf01')!h1g
      call rf01_forc_init(time_interp,As,Bs)

   case ('rf02')!h1g
       call rf02_forc_init(time_interp,As,Bs)

   case ('rico')!h1g
      call rico_forc_init(time_interp,As,Bs)

   case ('atex')!h1g
#ifdef BGRID
     if( present(qn) ) then!h1g
      call atex_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa,qn)!h1g
     else!h1g
      call atex_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa)!h1g
     endif!h1g
#endif

   case ('gcss_arm')!h1g
     call  gcss_arm_forc_init(time_interp,As,Bs)
     if ( do_specified_albedo )    ALBEDO_OBS  = 0.2 
     if ( do_specified_tskin )     TSKIN       = 300.0
     if ( do_specified_rough_leng ) then
          rough_mom  = 0.035
          rough_heat = 0.035/7.4
    endif

   case ('sheba')!h1g
#ifdef BGRID
     if( present(qn) .and.  present(qni) ) then!h1g
       call  sheba_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa,qn, qni)!h1g
     elseif( present(qn) ) then !h1g
       call  sheba_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa,qn)!h1g
     else
       call  sheba_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa)!h1g
     endif!h1g
#endif

     if ( do_specified_albedo )    ALBEDO_OBS = 0.827
     if ( do_specified_tskin )     TSKIN      = 257.4
     if ( do_specified_rough_leng ) then
          rough_mom  = 4.0e-4
          rough_heat = 0.0
    endif

! ---> h1g, 2011-02-25
   case ('gcss_astex')!h1g
#ifdef BGRID
     if( present(qn) ) then!h1g
      call gcss_astex_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa,qn)!h1g
     else!h1g
      call gcss_astex_forc_init(time_interp,As,Bs,elev,ps,T,u,v,qv,ql,qi,qa)!h1g
    endif!h1g
#endif
! <--- h1g, 2011-02-25


! ---> h1g, 2011-12-29
   case ('vocals')!h1g
      call vocals_forc_init(time_interp,As,Bs)
! <--- h1g, 2011-12-29


! ---> h1g, 2012-2-27
   case ('Vocals_75W_20S')!h1g
      call Vocals_75W_20S_forc_init(time_interp,As,Bs)
! <--- h1g, 2012-2-27

! ---> h1g, 2012-3-2
   case ('Hawaii_143W_8N')!h1g
      call Hawaii_143W_8N_forc_init(time_interp,As,Bs)
! <--- h1g, 2012-3-2

   case default 
       call mpp_error ('scm_forc_init',                          &
                       'Experiment name '//trim(experiment)//' not supported.', &
                       FATAL )

end select

end subroutine scm_forc_init

!#######################################################################

subroutine scm_forc_end ()

character(len=64)    :: fname_nc='RESTART/scm_forc.res.nc'

select case (trim(experiment))

   case ('arm')

      call arm_forc_end()

! ---> h1g
  case ('arm2003')
      call arm2003_forc_end()
! <--- h1g

   case ('AMIE')
      call AMIE_forc_end()

! copy from Lin Yanluan
   case ('twp')
      call twp_forc_end()

   case ('gate')
      call gate_forc_end()

   case ('mc3e')
      call mc3e_forc_end()

   case ('astex')

      call astex_forc_end()

   case ('bomex')

      call bomex_forc_end()

   case ('dcbl')

      call dcbl_forc_end()

   case ('mpacea')

      call mpacea_forc_end()

   case ('mpaceb')

      call mpaceb_forc_end()

   case ('rf01')

      call rf01_forc_end()

   case ('rf02')

      call rf02_forc_end()

   case ('rico')

      call rico_forc_end()

! ---> h1g
   case ('atex')
      call atex_forc_end()
! <--- h1g


! ---> h1g, 2010-06-21
   case ('gcss_arm')
      call  gcss_arm_forc_end ()
! <--- h1g, 2010-06-21


! ---> h1g, 2011-01-18
   case ('sheba')
      call  sheba_forc_end()
! <--- h1g, 2011-01-18


! ---> h1g, 2011-02-25
  case ('gcss_astex')
      call gcss_astex_forc_end()
! <--- h1g, 2011-02-25

! ---> h1g, 2011-12-29
  case ('vocals')
      call vocals_forc_end()
! <--- h1g, 2011-12-29


! ---> h1g, 2012-2-27
  case ('Vocals_75W_20S')
      call Vocals_75W_20S_forc_end()
! <--- h1g, 2012-2-27

! ---> h1g, 2012-3-2
  case ('Hawaii_143W_8N')
      call Hawaii_143W_8N_forc_end()
! <--- h1g, 2012-3-2

    case  default
       call mpp_error ('scm_forc_end',                          &
                       'Forcing destructor for '//trim(experiment)//' not called.', &
                       FATAL )
end select

!--------- restart file --------

if (mpp_pe() == mpp_root_pe()) then
      call write_data(fname_nc, 'SENFLUX',    SENFLUX   )
      call write_data(fname_nc, 'EVAPFLUX',   EVAPFLUX  )
      call write_data(fname_nc, 'TSKIN',      TSKIN     )
      call write_data(fname_nc, 'ALBEDO_OBS', ALBEDO_OBS)
      call write_data(fname_nc, 'ROUGH_MOM',  ROUGH_MOM )
      call write_data(fname_nc, 'ROUGH_HEAT', ROUGH_HEAT)
end if

end subroutine scm_forc_end

!#######################################################################

subroutine scm_forc_diagnostic_init(axes, Time)
implicit none

integer, dimension(:),intent(in) :: axes
type(time_type), intent(in)      :: Time

select case (trim(experiment))

   case ('arm')

      call arm_forc_diagnostic_init(axes, Time)

! ---> h1g
   case ('arm2003')
      call arm2003_forc_diagnostic_init(axes, Time)
! <--- h1g

   case ('AMIE')
      call AMIE_forc_diagnostic_init(axes, Time)

! Copy from Lin Yanluan, 2011-09-23
   case ('twp')
      call twp_forc_diagnostic_init(axes, Time)

   case ('gate')
      call gate_forc_diagnostic_init(axes, Time)

   case ('mc3e')
      call mc3e_forc_diagnostic_init(axes, Time)

   case ('astex')

      call astex_forc_diagnostic_init(axes, Time)

   case ('bomex')

      call bomex_forc_diagnostic_init(axes, Time)

   case ('dcbl')

      call dcbl_forc_diagnostic_init(axes, Time)

   case ('mpacea')

      call mpacea_forc_diagnostic_init(axes, Time)

   case ('mpaceb')

      call mpaceb_forc_diagnostic_init(axes, Time)

   case ('rf01')

      call rf01_forc_diagnostic_init(axes, Time)

   case ('rf02')

      call rf02_forc_diagnostic_init(axes, Time)

   case ('rico')

      call rico_forc_diagnostic_init(axes, Time)

! ---> h1g
   case ('atex')
      call  atex_forc_diagnostic_init(axes, Time)
! <--- h1g

! ---> h1g, 2010-06-21
   case ('gcss_arm')
      call  gcss_arm_forc_diagnostic_init(axes, Time)
! <--- h1g, 2010-06-21

! ---> h1g, 2011-01-18
   case ('sheba')
      call  sheba_forc_diagnostic_init(axes, Time)
! <--- h1g, 2011-01-18

! ---> h1g, 2011-02-25
  case ('gcss_astex')
      call gcss_astex_forc_diagnostic_init(axes, Time)
! <--- h1g, 2011-02-25

! ---> h1g, 2011-12-29
  case ('vocals')
      call vocals_forc_diagnostic_init(axes, Time)
! <--- h1g, 2011-12-29

! ---> h1g, 2012-2-27
  case ('Vocals_75W_20S')
      call Vocals_75W_20S_forc_diagnostic_init(axes, Time)
! <--- h1g, 2012-2-27


! ---> h1g, 2012-3-2
  case ('Hawaii_143W_8N')
      call Hawaii_143W_8N_forc_diagnostic_init(axes, Time)
! <--- h1g, 2012-3-2

  case  default
       call mpp_error ('scm_forc_diagnostics_init',                          &
                       'Diagnostics initialization for '//trim(experiment)//' not called.', &
                       FATAL )

end select
  
end subroutine scm_forc_diagnostic_init

!#######################################################################

subroutine update_scm_forc( time_interp, time_diag, dt_int, pdamp, elev, do_five) !yzheng
implicit none

type(time_type), intent(in)              :: time_interp,time_diag,dt_int
real,  intent (in)                       :: pdamp
real,  intent (in),    dimension(:,:)    :: elev
logical,                      intent(in)  :: do_five !yzheng

! ---> h1g


integer ::  klev
! <--- h1g


select case (trim(experiment))

   case ('arm')
      call update_arm_forc(time_interp,time_diag,dt_int, pdamp)!, omega_f)   !h1g
      call get_arm_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT ) !h1g

! ---> h1g
   case ('arm2003')
      call update_arm2003_forc(time_interp,time_diag,dt_int, pdamp)!, omega_f)   !h1g
      call get_arm2003_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT ) !h1g



   case ('AMIE')
     call update_AMIE_forc(time_interp,time_diag,dt_int,pdamp)
     if ( do_specified_tskin ) call get_AMIE_sst( TSKIN )

! Copy from Lin Yanluan, 2011-09-23
   case ('twp')
     call update_twp_forc(time_interp,time_diag,dt_int,pdamp)
     if ( do_specified_tskin ) call get_twp_sst( TSKIN )
!     call get_twp_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT ) !h1g

   case ('gate')
     call update_gate_forc(time_interp,time_diag,dt_int,pdamp)
     if ( do_specified_tskin ) call get_gate_sst( TSKIN )

   case ('mc3e')
     call update_mc3e_forc(time_interp,time_diag,dt_int,pdamp)
     if ( do_specified_tskin ) call get_mc3e_sst( TSKIN )


   case ('astex') !h1g
    call update_astex_forc(time_interp,time_diag,dt_int,pdamp)
    if ( do_specified_tskin ) call get_astex_sst( TSKIN ) !h1g

   case ('bomex') !h1g
      call update_bomex_forc(time_interp,time_diag,dt_int,elev)
      !yzheng
      if (do_five) call update_bomex_forc_five()

   case ('dcbl') !h1g
      call update_dcbl_forc(time_interp,time_diag,dt_int,elev)

   case ('mpacea') !h1g
      call update_mpacea_forc(time_interp,time_diag,dt_int,pdamp)
      call get_mpacea_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT ) !h1g

   case ('mpaceb') !h1g
     call update_mpaceb_forc(time_interp,time_diag,dt_int,pdamp)
     if ( do_specified_tskin ) call get_mpaceb_sst( TSKIN ) !h1g

   case ('rf01') !h1g
      call update_rf01_forc(time_interp,time_diag,dt_int)
      !yzheng
      if (do_five) call update_rf01_forc_five()
      
   case ('rf02') !h1g
     call update_rf02_forc(time_interp,time_diag,dt_int)

   case ('rico') !h1g
      call update_rico_forc(time_interp,time_diag,dt_int)

   case ('atex') !h1g
      call update_atex_forc(time_interp,time_diag,dt_int,elev)

! ---> h1g, 2010-06-21
   case ('gcss_arm') !h1g
      call update_gcss_arm_forc( time_interp,  time_diag)

    case ('sheba') !h1g
      call update_sheba_forc(time_interp,time_diag,dt_int)
 
! ---> h1g, 2011-02-15
    case ('gcss_astex') !h1g
      call update_gcss_astex_forc(time_interp,time_diag,dt_int,elev)

! ---> h1g, 2011-12-29
    case ('vocals') !h1g
      call update_vocals_forc( time_interp,time_diag,dt_int)

! ---> h1g, 2012-2-27
    case ('Vocals_75W_20S') !h1g
      call update_Vocals_75W_20S_forc( time_interp,time_diag,dt_int)

! ---> h1g, 2012-3-2
    case ('Hawaii_143W_8N') !h1g
      call update_Hawaii_143W_8N_forc( time_interp,time_diag,dt_int)

! <--- h1g, 2012-3-2
    case  default
       call mpp_error ('update_scm_forc',                          &
                       'Forcing update for '//trim(experiment)//' not supported.', &
                       FATAL )

end select

end subroutine update_scm_forc

!########################################################################

subroutine add_scm_tdtlw( tdtlw )

real, dimension(:,:,:) :: tdtlw

select case (trim(experiment))

   case ('bomex')

   case ('dcbl')

   case ('rf01')

     call add_rf01_tdtlw( tdtlw )

   case ('rf02')

     call add_rf02_tdtlw( tdtlw )

   case ('rico')

 ! ---> h1g
   case ('atex')
     call add_atex_tdtlw( tdtlw )!h1g
! <--- h1g

! ---> h1g, 2011-02-15
   case ('gcss_astex')
     call add_gcss_astex_tdtlw( tdtlw )!h1g
! <--- h1g, 2011-02-15
   case  default
     call mpp_error ('add_scm_tdtlw',                          &
                     'tdtlw addition for '//trim(experiment)//' not implemented.', &
                      NOTE )

end select
 
end subroutine add_scm_tdtlw

!########################################################################

subroutine add_scm_tdtsw( tdtsw )

real, dimension(:,:,:) :: tdtsw

select case (trim(experiment))

   case ('bomex')

   case ('dcbl')

   case ('rf01')

     call add_rf01_tdtsw( tdtsw )

   case ('rf02')

     call add_rf02_tdtsw( tdtsw )

   case ('rico')

! ---> h1g
   case ('atex')  !h1g
      call add_atex_tdtsw( tdtsw )!h1g
! <--- h1g

! ---> h1g, 2011-02-15
   case ('gcss_astex')
     call add_gcss_astex_tdtsw( tdtsw )!h1g
! <--- h1g, 2011-02-15

   case  default
     call mpp_error ('add_scm_tdtsw',                          &
                     'tdtsw addition for '//trim(experiment)//' not implemented.', &
                     NOTE )

end select
 
end subroutine add_scm_tdtsw
 
!########################################################################

subroutine scm_surface_flux(                                           &
     t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
     flux_t, flux_q, flux_r, flux_u, flux_v,                           &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dt,        land,      seawater,     avail,                        &
     dhdt_surf_forland,  dedt_surf_forland,  dedq_surf_forland  )

    real, intent(inout),  dimension(:) :: &
     t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
     flux_t, flux_q, flux_r, flux_u, flux_v,                           &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dhdt_surf_forland,  dedt_surf_forland,  dedq_surf_forland
    real, intent(in) :: dt
    logical, intent(in), dimension(:) :: land,  seawater, avail

    real, parameter :: gust_const =  1.0
    real, parameter :: d622   = rdgas/rvgas
    real, parameter :: d378   = 1.-d622
    real            :: d608   = d378/d622

    real, dimension(size(t_atm)) :: rho, tv_atm
    
     real, dimension(size(t_atm)) ::   thlm_atm,  wpthlp_sfc, wprtp_sfc   
     real           ::  sum_flux
     integer      ::  isum, nsum

    real, dimension(size(t_surf)) :: tsurf_mc3e
    real, dimension(size(flux_q)) :: fluxq_mc3e
! ------------------------------------------------------------------------------------------

    nsum =   size(t_atm)
    select case (trim(experiment))

     case ('arm')

      call get_arm_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )
      call arm_surface_flux_loop(                                           &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0

! --->h1g
    case ('arm2003')

      call get_arm2003_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )
      call arm2003_surface_flux_loop(                                           &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0
! <--- h1g

! copy from Lin Yanluan
     case ('twp')
     ! Call get_twp_sfc in order to get the latest SENFLUX and EVAPFLUX. 
     ! The other variables here are not actually used.
      call get_twp_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )
      call twp_surface_flux_loop(                                           &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0
!     call error_mesg('scm_forc_mod','scm_surface_flux not applicable to TWPICE',FATAL)

     case ('gate')
      call get_gate_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )
!      t_surf  = TSKIN
!      t_ca    = TSKIN
      call gate_surface_flux_loop(                                           &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0

     case ('mc3e')
      call get_mc3e_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )
!      tsurf_mc3e(:) = TSKIN
!      fluxq_mc3e(:) = EVAPFLUX
!      t_surf  = TSKIN
!      t_ca    = TSKIN
      call mc3e_surface_flux_loop(                                           &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
!          p_surf,    tsurf_mc3e,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
!          flux_t, fluxq_mc3e, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0



     case ('bomex')
      tv_atm = t_atm  * (1.0 + d608*q_atm)
      rho = p_atm / (rdgas * tv_atm)
      call get_bomex_flx( rho, u_star, flux_t, flux_q )

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      dtaudu_atm = 0.0
      dtaudv_atm = 0.0

      flux_q  = flux_q/hlv  ! latent heat convert from (W/m2) to (kg/m2/s)
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale

     case ('dcbl')
      tv_atm = t_atm  * (1.0 + d608*q_atm)
      rho = p_atm / (rdgas * tv_atm)
      call get_dcbl_flx( u_star, wpthlp_sfc, wprtp_sfc)

      thlm_atm = t_atm/( p_atm/ 1000.0e2 )**(rdgas/cp_air)

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      dtaudu_atm = 0.0
      dtaudv_atm = 0.0

      flux_t  = wpthlp_sfc * rho * cp_air
      flux_q  = wprtp_sfc * rho
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! ZNT 05/19/2020: Assuming no evap -> buoyancy flux ~ Sensible heat flux
      b_star = grav/thlm_atm * flux_t/(cp_air*rho) /u_star !buoyancy scale
      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale


     case ('rf01')
        thlm_atm = t_atm/( p_atm/ 1000.0e2 )**(rdgas/cp_air)
        tv_atm = t_atm  * (1.0 + d608*q_atm)
        rho = p_atm / (rdgas * tv_atm)

       call get_rf01_flx( u_star, flux_t, flux_q )
      !  call  get_rf01_flx_online( u_atm, v_atm,  thlm_atm, q_atm,  rho, & 
      !                                          flux_t, flux_q,  u_star)
          
      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      dtaudu_atm = 0.0
      dtaudv_atm = 0.0

      flux_q  = flux_q/hlv  ! latent heat convert from (W/m2) to (kg/m2/s)
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star

      q_star = flux_q/rho / u_star  ! moisture scale

     case ('rf02')

      tv_atm = t_atm  * (1.0 + d608*q_atm)
      call get_rf02_flx( rho, u_star, flux_t, flux_q )

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      dtaudu_atm = 0.0
      dtaudv_atm = 0.0

      flux_q  = flux_q/hlv  ! latent heat convert from (W/m2) to (kg/m2/s)
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale

! ---> h1g
     case ('rico')

! ---> h1g, 2010-09-14, using online surface fluxes
         if ( use_rico_online_sfcflux ) then
             call rico_sfcflux_online (&
               t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
               p_surf,    t_surf,     t_ca,      q_surf,                         &
               u_surf,    v_surf,                                                &
               rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
               flux_t, flux_q, flux_r, flux_u, flux_v,                           &
               cd_m,      cd_t,       cd_q,                                      &
               w_atm,     u_star,     b_star,     q_star,                        &
               dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
               dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
               dt,        land,      seawater,     avail  )


         else     
              thlm_atm = t_atm/( p_atm/ 1000.0e2 )**(rdgas/cp_air)
              call   rico_sfclyr ( u_atm, v_atm, thlm_atm, q_atm, &
                                          20.08,  299.8,  101540.,           & 
                                          flux_u,         flux_v,                 &
                                          wpthlp_sfc,  wprtp_sfc) 
                       u_star = 0.3
                       w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )

                       u_star =  sqrt( - flux_u *  w_atm / u_atm )
 
                       tv_atm = t_atm  * (1.0 + d608*q_atm)
                       rho = p_atm / (rdgas * tv_atm)
                      ! ZNT 05/07/2020: flux_u and flux_v should be rho-weighted
                      flux_u = flux_u*rho
                      flux_v = flux_v*rho

                      flux_t   = wpthlp_sfc * rho * cp_air
                      flux_q  = wprtp_sfc * rho

                      ! b_star = grav/292. * flux_t/(rho * cp_air) /u_star !buoyancy scale
                      ! ZNT 05/19/2020: including virtual effect
                      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
                      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                                              flux_q/rho*d608*t_atm) /u_star
                      q_star = flux_q/rho / u_star  ! moisture scale

                     dhdt_atm    = 0.0
                     dedq_atm   = 0.0
                     dtaudu_atm = 0.0
                     dtaudv_atm = 0.0

 
                    !  print*, ' clubb  surface flux'
                      sum_flux = 0.0
                      do isum = 1,  nsum 
                          sum_flux = sum_flux +  flux_t (isum )
                      enddo
      
                      sum_flux = 0.0
                      do isum = 1,  nsum 
                          sum_flux = sum_flux +  flux_q (isum )
                      enddo
         endif
! <--- h1g, 2010-09-14, using online surface fluxes


     case ('atex')
       wpthlp_sfc = 0.0
       wprtp_sfc  = 0.0
       thlm_atm = t_atm/( p_atm/ 1000.0e2 )**(rdgas/cp_air) 
       call    atex_sfclyr( u_atm, v_atm, thlm_atm, q_atm, &
                            flux_u, flux_v, wpthlp_sfc, wprtp_sfc)
  
       u_star = 0.3
       w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )
  
       tv_atm = t_atm  * (1.0 + d608*q_atm)
       rho = p_atm / (rdgas * tv_atm)

       ! flux_u  = - u_atm*u_star*u_star/w_atm
       ! flux_v  = - v_atm*u_star*u_star/w_atm 
       flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
       flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor

       flux_t   = wpthlp_sfc * rho * cp_air
       flux_q  = wprtp_sfc * rho

       ! b_star = grav/292. * flux_t/(rho * cp_air) /u_star !buoyancy scale
       ! ZNT 05/19/2020: including virtual effect
       ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
       b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                               flux_q/rho*d608*t_atm) /u_star

       q_star = flux_q/rho / u_star  ! moisture scale

       dhdt_atm    = 0.0
       dedq_atm    = 0.0
       dtaudu_atm  = 0.0
       dtaudv_atm  = 0.0
! <--- h1g

     case ('mpacea')

      call get_mpacea_sfc( ALBEDO_OBS, TSKIN, SENFLUX, EVAPFLUX, ROUGH_MOM, ROUGH_HEAT )
      call mpacea_surface_flux_loop(                                        &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0


     case ('mpaceb')

      call get_mpaceb_flx( SENFLUX, EVAPFLUX )
      call mpaceb_surface_flux_loop(                                        &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )
 
      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0

! ---> h1g, 2010-06-21
    case ('gcss_arm')

   if (use_gcss_arm_sfcflx_old) then
       tv_atm = t_atm  * (1.0 + d608*q_atm)
       rho = p_atm / (rdgas * tv_atm)
       thlm_atm = t_atm/( p_atm/ 1000.0e2 )**(rdgas/cp_air)
 
      call  get_gcss_arm_flx( z_atm, rho,  thlm_atm, t_atm, q_atm, u_atm, v_atm,    &
                              u_star, flux_t, flux_q, flux_u, flux_v )        

       ! ZNT 05/07/2020: flux_u and flux_v should be rho-weighted
       flux_u = flux_u*rho
       flux_v = flux_v*rho

       flux_q = flux_q / hlv  ! latent heat convert from (W/m2) to (kg/m2/s)

       ! b_star = grav/thlm_atm * flux_t/( rho * cp_air ) /u_star !buoyancy scale
       ! ZNT 05/19/2020: including virtual effect
       ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
       b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                               flux_q/rho*d608*t_atm) /u_star
       ! write(*,*) tv_atm, rho, t_atm, q_atm, d608, flux_t, flux_q, u_star, b_star
       q_star = flux_q/ rho / u_star  ! moisture scale

      if (do_specified_land) then
          dhdt_surf_forland = dhdt_surf
          dedt_surf_forland = dedt_surf
          dedq_surf_forland = dedq_surf
       end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0

      if ( do_specified_albedo )    ALBEDO_OBS = 0.2
      if ( do_specified_tskin )     TSKIN      = 300.0
      if ( do_specified_rough_leng ) then
          rough_mom     = 0.035
          rough_heat    = 0.035 / 7.4
      endif

   else
      rough_scale = rough_m
      call gcss_arm_surface_flux_loop(                                      &
          t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
          p_surf,    t_surf,     t_ca,      q_surf,                         &
          u_surf,    v_surf,                                                &
          rough_m,   rough_h,    rough_moist, rough_scale, gust,            &
          flux_t, flux_q, flux_r, flux_u, flux_v,                           &
          cd_m,      cd_t,       cd_q,                                      &
          w_atm,     u_star,     b_star,     q_star,                        &
          dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
          dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
          dt,        land,      seawater,     avail  )

      if (do_specified_land) then

        dhdt_surf_forland = dhdt_surf
        dedt_surf_forland = dedt_surf
        dedq_surf_forland = dedq_surf

      end if

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0

      if ( do_specified_albedo )    ALBEDO_OBS = 0.2
      if ( do_specified_tskin )     TSKIN      = 300.0
      if ( do_specified_rough_leng ) then
          rough_mom     = 0.035
          rough_heat      = 0.035 / 7.4
      endif
   endif

! ---> h1g, 2011-01-18
      case ('sheba')

      tv_atm = t_atm  * (1.0 + d608*q_atm)
      rho = p_atm / (rdgas * tv_atm)
      call get_sheba_flx( u_star, flux_t, flux_q )

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      dtaudu_atm = 0.0
      dtaudv_atm = 0.0

      flux_q  = flux_q / hlv        ! latent heat convert from (W/m2) to (kg/m2/s)
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale

      dhdt_surf = 0.0
      dedt_surf = 0.0
      dedq_surf = 0.0
      dhdt_atm  = 0.0
      dedq_atm  = 0.0

      if ( do_specified_albedo )    ALBEDO_OBS  = 0.827
      if ( do_specified_tskin )     TSKIN       = 257.4
      if ( do_specified_rough_leng ) then
           rough_mom   = 4.0e-4
           rough_heat    = 0.0
     endif

     case ('gcss_astex')
        tv_atm = t_atm  * (1.0 + d608*q_atm)
        rho = p_atm / (rdgas * tv_atm)

        call get_gcss_astex_flx( u_star, flux_t, flux_q )
        w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
        ! flux_u  = - u_atm*u_star*u_star/w_atm
        ! flux_v  = - v_atm*u_star*u_star/w_atm
        flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
        flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
        dtaudu_atm = 0.0
        dtaudv_atm = 0.0

        dhdt_atm = 0.0
        dedq_atm = 0.0
  
        ! b_star = grav/292.* flux_t/ u_star
        ! ZNT 05/19/2020: including virtual effect
        ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
        b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                                flux_q/rho*d608*t_atm) /u_star
        q_star = flux_q / u_star  ! moisture scale


      case ('vocals')
      tv_atm = t_atm  * (1.0 + d608*q_atm)
      rho = p_atm / (rdgas * tv_atm)
      call get_vocals_flx( u_star, flux_t, flux_q )

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor

      dtaudu_atm = 0.0
      dtaudv_atm = 0.0
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale

! --->h1g, 2012-2-27
      case ('Vocals_75W_20S')
      tv_atm = t_atm  * (1.0 + d608*q_atm)
      rho = p_atm / (rdgas * tv_atm)
      call get_Vocals_75W_20S_flx( u_star, flux_t, flux_q )

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor

      dtaudu_atm = 0.0
      dtaudv_atm = 0.0
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale
! <---h1g, 2012-2-27



! --->h1g, 2012-3-2
      case ('Hawaii_143W_8N')
      tv_atm = t_atm  * (1.0 + d608*q_atm)
      rho = p_atm / (rdgas * tv_atm)
      call get_Hawaii_143W_8N_flx( u_star, flux_t, flux_q )

      w_atm = max( sqrt(u_atm*u_atm + v_atm*v_atm), gust_const )  
      ! flux_u  = - u_atm*u_star*u_star/w_atm
      ! flux_v  = - v_atm*u_star*u_star/w_atm
      flux_u  = - rho*u_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor
      flux_v  = - rho*v_atm*u_star*u_star/w_atm   ! ZNT 05/07/2020: rho factor

      dtaudu_atm = 0.0
      dtaudv_atm = 0.0
      dhdt_atm = 0.0
      dedq_atm = 0.0

      ! b_star = grav/292. * flux_t/(cp_air*rho) /u_star
      ! ZNT 05/19/2020: including virtual effect
      ! b_star = grav/Tv*Tv_star = grav/Tv*(T_star*(1+0.608*q)+q_star*0.608*T)
      b_star = grav/tv_atm * (flux_t/(cp_air*rho)*(1.0+d608*q_atm) + &
                              flux_q/rho*d608*t_atm) /u_star
      q_star = flux_q/rho / u_star  ! moisture scale
! <---h1g, 2012-3-2

    case  default
       call mpp_error ('scm_surface_flux',                          &
                       'Surface flux for '//trim(experiment)//' not supported.', &
                       FATAL )

 end select

end subroutine scm_surface_flux

!#######################################################################

end module scm_forc_mod
