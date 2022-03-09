
module scm_grid_masks_mod

!-----------------------------------------------------------------------
!
!    allocates storage and initializes masks and vertical indexing
!    associated with the step mountain vertical coordinate.
!
!-----------------------------------------------------------------------

use scm_horiz_grid_mod, only:  horiz_grid_type
use   scm_boundary_mod, only:  update_halo, NOPOLE, UWND
use  scm_vert_grid_mod, only:  vert_grid_type
use            fms_mod, only:  mpp_pe
use            mpp_mod, only:  mpp_max
implicit none
private

!-----------------------------------------------------------------------
!    ---- public data types ----

   public      mask_type
   public grid_mask_type

   type mask_type
      real,    pointer, dimension(:,:,:) :: mask=>NULL()
      integer, pointer, dimension(:,:)   :: kbot=>NULL()
      integer                            :: kbotmin
   end type mask_type

!  mask  = step-mountain topography mask (0.0 or 1.0) for
!            mass (height) grid points
!  kbot  = lowest model level above ground
!
!   note:  for the sigma coordinate, mask = 1.0 everywhere, and
!          kbot = number of vertical model levels

   type grid_mask_type
      type(mask_type) :: Tmp, Vel
      logical :: sigma
   end type grid_mask_type

!  Tmp = grid masking values for the temperature/mass grid
!  Vel = grid masking values for the velocity/momentum grid
!  sigma = logical flag that specific whether vertical coordinate is
!            the step-mountain (eta) or sigma coordinate
!
!-----------------------------------------------------------------------
!    ---- public interfaces (and some that may become public) ----

   public  grid_masks_init

   private compute_mass_mask, compute_vel_mask, compute_lowest_level

   logical :: sigma  ! local variable
!-----------------------------------------------------------------------

contains

!#######################################################################

   function grid_masks_init (Hgrid, Vgrid, res) result (Mask)

!-----------------------------------------------------------------------
!
!     arguments (intent in)
!
!     res  =  reciprical of eta at the surface
!
!-----------------------------------------------------------------------

   type (horiz_grid_type), intent(inout) :: Hgrid
   type  (vert_grid_type), intent(in)    :: Vgrid
   real,                   intent(in)    :: res(Hgrid%ilb:,Hgrid%jlb:)

   type (grid_mask_type) :: Mask

!-----------------------------------------------------------------------
!   sigma = logical flag for terrian-following versus step-mountain
!           vertical coordinate.

   logical :: sigma
   integer :: kx
   real    :: maxres

   kx = Vgrid % nlev

!--------------allocate global storage----------------------------------

 allocate ( Mask%Tmp%mask (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, kx), &
            Mask%Vel%mask (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, kx), &
            Mask%Tmp%kbot (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub),     &
            Mask%Vel%kbot (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub)      )

!--------------is this an eta of sigma coordinate mountain? ------------

      maxres = maxval(res)
      call mpp_max (maxres)

      if (maxres > 1.0001) then
         Mask % sigma=.false.
         sigma=.false.
         if (mpp_pe() == 0) write (*,100)
      else
         Mask % sigma=.true.
         sigma=.true.
         if (mpp_pe() == 0) write (*,200)
      endif

 100  format (' B-grid dynamical core has been initialized with eta coordinate data.')
 200  format (' B-grid dynamical core has been initialized with sigma coordinate data.')

!--------------topography masks ----------------------------------------

      Mask % Tmp % mask = compute_mass_mask (res, Vgrid % aeta)
      Mask % Vel % mask = compute_vel_mask  (res, Vgrid % aeta)
      call update_halo (Hgrid, UWND, Mask%Vel%mask, flags=NOPOLE)
!!!!! call update_halo (Hgrid, UWND, Mask%Vel%mask)   ! sets mask=0 at poles

!------------- compute the lowest model level --------------------------

      Mask % Tmp % kbot = compute_lowest_level (Mask % Tmp % mask)
      Mask % Vel % kbot = compute_lowest_level (Mask % Vel % mask)

!     ----- global values -----

      Mask % Tmp % kbotmin = minval(Mask % Tmp % kbot)
      Mask % Vel % kbotmin = minval(Mask % Vel % kbot)

!-----------------------------------------------------------------------

end function grid_masks_init

!#######################################################################

   function compute_mass_mask (res, aeta) result (mask)

   real, intent(in) :: res(:,:), aeta(:)
   real, dimension(size(res,1),size(res,2),size(aeta)) :: mask
   integer  i, j, k

      mask = 1.0

      if (.not.sigma) then
         do j=1,size(res,2); do i=1,size(res,1)
         do k=1,size(aeta)
            if (aeta(k) > (1.0/res(i,j))) mask(i,j,k) = 0.0
         enddo; enddo; enddo
      endif

   end function compute_mass_mask

!#######################################################################

   function compute_vel_mask (res, aeta) result (mask)

   real, intent(in) :: res(:,:), aeta(:)
   real, dimension(size(res,1),size(res,2),size(aeta)) :: mask
   integer  i, j, k

      mask = 1.0

      if (.not.sigma) then
         do j=2,size(res,2); do i=2,size(res,1)
         do k=1,size(aeta)
            if (aeta(k) > (1.0/res(i,j))) then
                mask(i-1,j-1,k) = 0.0
                mask(i  ,j-1,k) = 0.0
                mask(i-1,j  ,k) = 0.0
                mask(i  ,j  ,k) = 0.0
            endif
         enddo; enddo; enddo
      endif

   end function compute_vel_mask

!#######################################################################

   function compute_lowest_level (mask) result (kbot)

   real, intent(in) :: mask(:,:,:)
   integer, dimension(size(mask,1),size(mask,2)) :: kbot
   integer   i, j, k, kdim

      kdim = size(mask,3)
      kbot = kdim

      if (.not.sigma) then
         do j=1,size(mask,2); do i=1,size(mask,1)
         do k=1,kdim
            if (mask(i,j,k) < 0.50) then
               kbot(i,j)=k-1; exit
            endif
         enddo; enddo; enddo
       ! must not be zero
         kbot = max(kbot,1)
      endif

   end function compute_lowest_level

!#######################################################################

end module scm_grid_masks_mod

