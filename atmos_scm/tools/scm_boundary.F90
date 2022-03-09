
module scm_boundary_mod

use  scm_horiz_grid_mod, only: horiz_grid_type, update_np, update_sp, &
                             TGRID, VGRID
use         fms_mod, only: error_mesg, FATAL
use mpp_domains_mod, only: mpp_update_domains, SUPDATE, NUPDATE,  &
                                               WUPDATE, EUPDATE

implicit none
private


!------------ public interfaces ------------

public :: update_halo,  horiz_grid_boundary, vel_flux_boundary

interface update_halo
    module procedure  update_halo_2d, update_halo_3d, update_halo_4d
end interface

interface vel_flux_boundary
    module procedure  vel_flux_boundary_2d, vel_flux_boundary_3d
end interface

!  public parameters

integer, parameter, public :: TEMP = 21, UWND = 22, VWND = 23

integer, parameter, public :: SOUTH = 1, NORTH = 2
integer, parameter, public ::  WEST = 4,  EAST = 8
integer, parameter         ::  ALL = SOUTH+NORTH+WEST+EAST

integer, parameter, public :: NOPOLE   = 16
integer, parameter, public :: POLEONLY = 32


contains

!#######################################################################

 subroutine update_halo_3d (Hgrid, field, data, flags)

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(Hgrid%ilb:,Hgrid%jlb:,:)
   integer, optional,     intent(in)    :: flags

   integer :: is, ie, iflags, n, xygrid
   integer :: domain_flags
   logical :: no_pole_vel, do_pole_only, update_sbnd, update_nbnd
   logical :: decomp_2d

!  ----- check dimensions ------

   if (size(data,2) /= Hgrid % jsize)  call error_mesg  &
             ('update_halo', 'j dimension has wrong size', FATAL)

   if (size(data,1) /= Hgrid % isize)  call error_mesg  &
                ('update_halo', 'i dimension has wrong size', FATAL)

!  ----- check/set optional flag arguments ----

   iflags = ALL
   if (present(flags)) iflags = flags

   if ( iflags <= 0 .or. iflags > ALL+NOPOLE+POLEONLY ) &
   call error_mesg ('update_halo', 'invalid value for flags', FATAL)

   domain_flags = 0

!  ------ need to determine and check grid -------

   select case (field)
     case (TEMP)
        is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
        xygrid = TGRID
     case (UWND:VWND)
        is = Hgrid % Vel % is;  ie = Hgrid % Vel % ie
        xygrid = VGRID
     case default
        call error_mesg ('update_halo', 'invalid field', FATAL)
   end select

!  ------ south and north boundary ------

   if ( btest(iflags,0) ) domain_flags = domain_flags + SUPDATE
   if ( btest(iflags,1) ) domain_flags = domain_flags + NUPDATE

!  ------ west and east boundary ------

   if ( btest(iflags,2) ) then
      if (Hgrid%decompx) domain_flags = domain_flags + WUPDATE
   endif

   if ( btest(iflags,3) ) then
      if (Hgrid%decompx) domain_flags = domain_flags + EUPDATE
   endif

!  ------ flags related to polar boundary ------

   no_pole_vel  = btest(iflags,4)
   do_pole_only = btest(iflags,5)
   update_sbnd  = btest(iflags,0)
   update_nbnd  = btest(iflags,1)


!  ----- update non-polar boundaries -----

   if (.not.do_pole_only) then
     select case (field)
       case (TEMP)
         call mpp_update_domains (data, Hgrid%Tmp%Domain, domain_flags)
       case (UWND:VWND)
         call mpp_update_domains (data, Hgrid%Vel%Domain, domain_flags)
     end select
   endif

!  ----- update east-west cyclic boundaries (for 1-d decomp only) ----

   if (.not.Hgrid%decompx) then
       n = Hgrid%ihalo
     ! update west halo
       if (btest(iflags,2)) data(is-n:is-1,:,:) = data(ie-n+1:ie,:,:)
     ! update east halo
       if (btest(iflags,3)) data(ie+1:ie+n,:,:) = data(is:is+n-1,:,:)
   endif

!  ------ update south pole ------

   if ( (update_sbnd.or.do_pole_only) .and. update_sp (Hgrid,xygrid) ) then
      call south_boundary_3d (Hgrid, field, data(:,:,:), no_pole_vel)
   endif

!  ------ update north pole ------

   if ( (update_nbnd.or.do_pole_only) .and. update_np (Hgrid,xygrid) ) then
      call north_boundary_3d (Hgrid, field, data(:,:,:), no_pole_vel)
   endif


 end subroutine update_halo_3d

!#######################################################################

 subroutine north_boundary_3d (Hgrid, field, data, nopole)

   type(horiz_grid_type), intent(in)    :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(:,Hgrid%jlb:,:)
   logical,               intent(in)    :: nopole

   integer :: js, je, jeg, halo

      halo = Hgrid % jhalo

! --- update north pole boundary ---

      select case (field)
         case (TEMP)
!        ---- mass ----
            js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je
            data (:, je+1:je+halo, :) = data (:, je:je-halo+1:-1, :)
         case (UWND)
!        ---- u comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je;  jeg = Hgrid % Vel % jeg
            if (.not. nopole) then
                if ( jeg+1 <= je+halo ) data (:, jeg+1,:) = 0.0
            endif
            if ( jeg+2 <= je+halo ) &
            data (:, jeg+2:jeg+halo, :) = data (:, je:je-halo+2:-1, :)
         case (VWND)
!        ---- v comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je;  jeg = Hgrid % Vel % jeg
            if (.not. nopole) then
                if ( jeg+1 <= je+halo ) data (:, jeg+1, :) = 0.0
            endif
            if ( jeg+2 <= je+halo ) &
            data (:, jeg+2:jeg+halo, :) = - data (:, je:je-halo+2:-1, :)
      end select

 end subroutine north_boundary_3d

!#######################################################################

 subroutine south_boundary_3d (Hgrid, field, data, nopole)

   type(horiz_grid_type), intent(in)    :: Hgrid
   integer,               intent(in)    ::  field
   real,                  intent(inout) :: data(:,Hgrid%jlb:,:)
   logical,               intent(in)    :: nopole
      

   integer :: js, je, halo

      halo = Hgrid % jhalo

! --- update south pole boundary ---

      select case (field)
         case (TEMP)
!        ---- mass ----
            js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je
            data (:, js-1:js-halo:-1, :) = data (:, js:js+halo-1, :)
         case (UWND)
!        ---- u comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je
            if (.not. nopole) data (:, js-1, :) = 0.0
            data (:, js-2:js-halo:-1, :) = data (:, js:js+halo-2, :)
         case (VWND)
!        ---- v comp ----
            js = Hgrid % Vel % js;  je = Hgrid % Vel % je
            if (.not. nopole) data (:, js-1, :) = 0.0
            data (:, js-2:js-halo:-1, :) = - data (:, js:js+halo-2, :)
      end select

 end subroutine south_boundary_3d

!#######################################################################

 subroutine update_halo_2d (Hgrid, field, data, flags)

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(:,:)
   integer, optional,     intent(in)    :: flags

   real, dimension(size(data,1),size(data,2),1) :: data3
   integer :: n

   data3(:,:,1) = data
   call update_halo_3d (Hgrid, field, data3, flags)
   data = data3(:,:,1)

 end subroutine update_halo_2d

!#######################################################################

 subroutine update_halo_4d (Hgrid, field, data, flags)

   type(horiz_grid_type), intent(inout) :: Hgrid
   integer,               intent(in)    :: field
   real,                  intent(inout) :: data(:,:,:,:)
   integer, optional,     intent(in)    :: flags

   integer :: n

   if (size(data,4) == 0) return

   do n = 1, size(data,4)
      call update_halo_3d (Hgrid, field, data(:,:,:,n), flags)
   enddo

 end subroutine update_halo_4d

!#######################################################################

 subroutine horiz_grid_boundary (Hgrid)

    type(horiz_grid_type), intent(inout) :: Hgrid

!---- updates boundaries for 2d arrays in horiz_grid_type -----

     call update_halo (Hgrid, UWND, Hgrid % Vel % dx, ALL+NOPOLE)
     call update_halo (Hgrid, TEMP, Hgrid % Tmp % dx)

     call update_halo (Hgrid, TEMP, Hgrid % Tmp %  area)
     call update_halo (Hgrid, TEMP, Hgrid % Tmp % rarea)

     call update_halo (Hgrid, UWND, Hgrid % Vel %  area, ALL+NOPOLE)
     call update_halo (Hgrid, UWND, Hgrid % Vel % rarea, ALL+NOPOLE)

 end subroutine horiz_grid_boundary

!#######################################################################

 subroutine vel_flux_boundary_2d (Hgrid, data)

   type(horiz_grid_type), intent(in)    :: Hgrid
   real,                  intent(inout) :: data(:,Hgrid%jlb:)

      if ( update_sp (Hgrid,VGRID) ) then
                                data (:, Hgrid%Vel%js  ) = 0.0
           if (Hgrid%jhalo > 1) data (:, Hgrid%Vel%js-1) = 0.0
      endif

      if ( Hgrid%jhalo == 0 ) return
      if ( update_np (Hgrid,VGRID) ) then
          if (Hgrid%Vel%jeg+1 <= Hgrid%Vel%je+Hgrid%jhalo) &
                                           data (:,Hgrid%Vel%jeg+1) = 0.
          if (Hgrid%Vel%jeg+2 <= Hgrid%Vel%je+Hgrid%jhalo) &
                                           data (:,Hgrid%Vel%jeg+2) = 0.
      endif

 end subroutine vel_flux_boundary_2d

!#######################################################################

 subroutine vel_flux_boundary_3d (Hgrid, data)

   type(horiz_grid_type), intent(in)    :: Hgrid
   real,                  intent(inout) :: data(:,Hgrid%jlb:,:)

      if ( update_sp (Hgrid,VGRID) ) then
                                data (:, Hgrid%Vel%js  , :) = 0.0
           if (Hgrid%jhalo > 1) data (:, Hgrid%Vel%js-1, :) = 0.0
      endif

      if ( Hgrid%jhalo == 0 ) return
      if ( update_np (Hgrid,VGRID) ) then
          if (Hgrid%Vel%jeg+1 <= Hgrid%Vel%je+Hgrid%jhalo) &
                                         data (:,Hgrid%Vel%jeg+1,:) = 0.
          if (Hgrid%Vel%jeg+2 <= Hgrid%Vel%je+Hgrid%jhalo) &
                                         data (:,Hgrid%Vel%jeg+2,:) = 0.
      endif

 end subroutine vel_flux_boundary_3d

!#######################################################################

end module scm_boundary_mod

