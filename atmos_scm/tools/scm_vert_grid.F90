
module scm_vert_grid_mod

!-----------------------------------------------------------------------
!
!     allocates storage and initializes vertical grid constants
!
!     contains several interfaces for compute pressure and height
!
!-----------------------------------------------------------------------

use constants_mod, only: grav, rdgas, rvgas
use       fms_mod, only: error_mesg, FATAL, mpp_pe

implicit none
private

!-----------------------------------------------------------------------
!  public defined data type (vert_grid_type)
!  -----------------------------------------
!    nlev  = number of vertical levels (integer)
!    nplev = number of pure pressure levels at the top of the model,
!            will equal zero when not using hybrid coordinate (integer)
!    deta  = vertical eta thickness/depth of model layers
!    aeta  = eta values at model levels (full levels)
!    eta   = eta values at model layer interfaces (half levels)
!    dpeta = vertical pressure thickness/depth of model layers
!    apeta = pressure values at model levels (full levels)
!    peta  = pressure values at model layer interfaces (half levels)
!    dfl   = reference values of geopotental height at half levels
!    
!    psmin  = minimum allowable surface pressure to avoid negative
!             mass in a model layer (important for hybrid coord)
!    hybrid = logical flag (true for hybrid coordinate)
!    pzero  = logical flag (true for pres = 0 at top of model)
!
!    pref,tref,gamma = reference values for computing dfl
!
!
!  public interfaces
!  -----------------
!    vert_grid_init      - initializes the vert_grid_type
!    compute_height_half - computes the geopotential height (in meters)
!                            at half model levels
!    compute_pres_depth  - computes the pressure weight (mass) of a
!                            model layer
!
!-----------------------------------------------------------------------
!------- interfaces -------

public  vert_grid_init, compute_geop_height, compute_height,   &
        compute_pres_depth, &
        compute_pres_weights, compute_pressures, compute_height_bottom

!------- public defined data type -------

public  vert_grid_type

type vert_grid_type
   integer                     :: nlev, nplev
   real, pointer, dimension(:) :: deta=>NULL(),  aeta=>NULL(),  eta=>NULL(),  dfl=>NULL(),  &
                                  dpeta=>NULL(), apeta=>NULL(), peta=>NULL(), wta=>NULL(), wtb=>NULL()
   real                        :: pref, tref, gamma, psmin
   logical                     :: hybrid, pzero
end type vert_grid_type

!-----------------------------------------------------------------------

  real, parameter :: d608 = (rvgas-rdgas)/rdgas
  real, parameter :: ginv = 1./grav

!------ parameters for eta coordinate reference surface heights --------

  real, parameter :: pref = 101325., tref = 288., gamma = 0.0065

!-----------------------------------------------------------------------

contains

!#######################################################################

function vert_grid_init (eta, peta, verbose) result (Vgrid)

!-----------------------------------------------------------------------
!
!    all arguments are input only
!    ----------------------------
!
!    deta_in  =  vertical grid spacing, an array of order one
!

     real, intent (in)           ::  eta(:)
     real, intent (in), optional :: peta(:)
  integer, intent (in), optional :: verbose

     type(vert_grid_type),target :: Vgrid

!-----------------------------------------------------------------------
integer :: k, kx, lverbose
real    :: rgog
real, dimension(size(eta))   :: phalf, lphalf, pres
real, dimension(size(eta)-1) ::        lpfull

!-----------------------------------------------------------------------

   lverbose = 0;  if (present(verbose)) lverbose = verbose

!--------------derived vertical constants-------------------------------

      kx = size(eta) - 1
      allocate (Vgrid% deta(kx), Vgrid% aeta(kx), Vgrid% eta(kx+1),  &
                Vgrid%dpeta(kx), Vgrid%apeta(kx), Vgrid%peta(kx+1),  &
                Vgrid%wta  (kx), Vgrid%wtb  (kx),                    &
                Vgrid% dfl(kx+1))

      Vgrid % nlev = kx

!--------- set-up eta values and hybrid pressure levels ----------
!--------- note: eta(1) and eta(kx+1) have set values -----
!--------- also note: peta(kx+1) = 0.0 -----

        Vgrid % eta(1)    = 0.0
        Vgrid % eta(kx+1) = 1.0
        Vgrid % eta(2:kx) = eta(2:kx)

        Vgrid % peta = 0.0
        if (present(peta)) Vgrid % peta(1:kx) = peta(1:kx)

      do k = 1, kx
        Vgrid %  deta(k) =  Vgrid %  eta(k+1) - Vgrid %  eta(k)
        Vgrid %  aeta(k) = (Vgrid %  eta(k+1) + Vgrid %  eta(k)) * 0.5
        Vgrid % dpeta(k) =  Vgrid % peta(k+1) - Vgrid % peta(k)
        Vgrid % apeta(k) =  0.0
        Vgrid % wta  (k) =  0.0
        Vgrid % wtb  (k) =  0.0
      enddo

!----------- is this a hybrid coordinate ??? -----

      Vgrid % hybrid = .false.

      do k = 1, kx+1
         if ( Vgrid % peta(k) > 0.0 ) then
              Vgrid % hybrid = .true.
              exit
         endif
      enddo

!----------- find lowest pure pressure level --------------

      Vgrid % nplev = 0

      do k = 1, kx
         if ( Vgrid % deta(k) > 0.0 ) exit
         Vgrid % nplev = k
      enddo

!     ---- need average pressure in these layers ----

      Vgrid % pzero = .true.

      if ( Vgrid % nplev >= 1 ) then
         phalf(:) = Vgrid%peta(:) + Vgrid%eta(:)*pref
         if ( phalf(1) <= epsilon(phalf) ) then
            lphalf(1) = 0.0
            lphalf(2:) = log(phalf(2:))
         else
            lphalf(:) = log(phalf(:))
            Vgrid % pzero = .false.
         endif

         do k = 1, kx
            lpfull(k) = (phalf(k+1)*lphalf(k+1) - phalf(k)*lphalf(k)) &
                        / (phalf(k+1)-phalf(k)) - 1.0
            Vgrid % apeta(k) = exp(lpfull(k))
            Vgrid % wtb  (k) = lphalf(k+1) - lpfull(k)
            Vgrid % wta  (k) = lpfull(k)   - lphalf(k)
         enddo
         if (Vgrid % pzero) Vgrid % wta(1) =Vgrid % wtb(1)

      endif

!----------- find the minimum allowable surface pressure ------

      Vgrid % psmin = 0.0

      do k = 1, kx
         if ( Vgrid % deta(k) > 0.0 )  Vgrid % psmin =  &
               max ( Vgrid % psmin, -Vgrid % dpeta(k)/Vgrid % deta(k) )
      enddo


!--- optional output of coordinate values ----

   if (lverbose > 1) then
       print *,  ' eta=',Vgrid %   eta
       print *,  'deta=',Vgrid %  deta
       print *,  'aeta=',Vgrid %  aeta
       print *, ' peta=',Vgrid %  peta
       print *, 'dpeta=',Vgrid % dpeta
       print *, 'apeta=',Vgrid % apeta
       print *, 'nplev=',Vgrid % nplev
       print *, 'minimum allowable surface pressure = ', Vgrid % psmin
   endif

!---------- set-up eta coordinate geopotential heights -------------

      rgog = rdgas*gamma/grav

      do k=1,kx
        pres(k) = Vgrid%peta(k) + Vgrid%eta(k)*pref
        Vgrid % dfl(k) = grav*tref*(1.0-(pres(k)/pref)**rgog)/gamma
      enddo
        Vgrid % dfl(kx+1) = 0.0

        Vgrid % pref  = pref
        Vgrid % tref  = tref
        Vgrid % gamma = gamma

!-----------------------------------------------------------------------

end function vert_grid_init

!#######################################################################

subroutine compute_pres_depth (Vgrid, pssl, pdepth)

  type(vert_grid_type), intent(in)  :: Vgrid
       real           , intent(in)  :: pssl(:,:)
       real           , intent(out) :: pdepth(:,:,:)

       integer :: k, kp, ke

!-----------------------------------------------------------------------
!   compute the pressure weight (depth) for model layers
!-----------------------------------------------------------------------

  kp = Vgrid % nplev
  ke = Vgrid % nlev

  if (size(pdepth,3) /= ke) call error_mesg (  &
                              'compute_pres_depth in vert_grid_mod', &
                              'incorrect dimension 3 for pdepth', FATAL)

! --- check for zero/negative depth layers ---

  if (Vgrid % hybrid) then
      if (minval(pssl) <= Vgrid % psmin) call error_mesg  &
                 ('compute_pres_depth in vert_grid_mod',  &
                  'pressure depth <= 0.0', FATAL)
  endif

! --- compute depth ---

  do k = 1, kp
    pdepth(:,:,k) = Vgrid % dpeta(k)
  enddo
    
  do k = kp+1, ke
    pdepth(:,:,k) = Vgrid % dpeta(k) + Vgrid % deta(k) * pssl(:,:)
  enddo

end subroutine compute_pres_depth

!#######################################################################

subroutine compute_pressures (Vgrid, pssl, phalf, pfull, dpde, wta, wtb)

  type(vert_grid_type), intent(in)  :: Vgrid
  real                , intent(in)  :: pssl(:,:)
  real                , intent(out) :: phalf(:,:,:), pfull(:,:,:)
  real, optional      , intent(out) :: dpde(:,:,:)
  real, optional      , intent(out) :: wta(:,:,:), wtb(:,:,:)

  real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)+1) :: ph,lph
  real, dimension(size(pfull,1),size(pfull,2),size(pfull,3))   :: dp,lpf
  integer :: k, kp, ke
!-----------------------------------------------------------------------
!      compute the pressure at full model levels
!-----------------------------------------------------------------------

  kp = Vgrid % nplev
  ke = Vgrid % nlev

  if (size(pfull,3) /= ke) call error_mesg (  &
                               'compute_pressures in vert_grid_mod', &
                               'incorrect dimension 3 for pfull', FATAL)

!--- set or compute optional arguments ---

  call compute_pres_half (Vgrid, pssl, ph)
  phalf = ph

  call compute_pres_depth (Vgrid, pssl, dp)
  if (present(dpde))  dpde = dp

!--- compute p*logp at half levels ---

  if ( Vgrid % pzero ) then
      lph(:,:,1) = 0.0
       ph(:,:,1) = 0.0
      lph(:,:,2:ke+1) = log(ph(:,:,2:ke+1))
       ph(:,:,2:ke+1) = ph(:,:,2:ke+1) * lph(:,:,2:ke+1)
  else
      lph(:,:,:) = log(ph(:,:,:))
       ph(:,:,:) = ph(:,:,:) * lph(:,:,:)
  endif

!--- compute pressure at full levels ---

  do k = 1, kp
    pfull(:,:,k) = Vgrid % apeta(k)
  enddo

! if (present(wta) .or. present(wtb)) then
!   do k = 1, kp
!     lpf(:,:,k) = log(pfull(:,:,k))
!   enddo
! endif

  do k = kp+1, ke
    lpf(:,:,k) = (ph(:,:,k+1)-ph(:,:,k))/dp(:,:,k) - 1.0
    pfull(:,:,k) = exp( lpf(:,:,k) )
  enddo
    
!--- compute weights at full levels ---

   if (present(wtb)) then
      do k = 1, kp
        wtb(:,:,k) = Vgrid % wtb(k)
      enddo
      do k = kp+1, size(wtb,3)
        wtb(:,:,k) = lph(:,:,k+1) - lpf(:,:,k)
      enddo
   endif

   if (present(wta)) then
      do k = 1, kp
        wta(:,:,k) = Vgrid % wta(k)
      enddo
      do k = kp+1, size(wta,3)
        wta(:,:,k) = lpf(:,:,k) - lph(:,:,k)
      enddo
      if (Vgrid % pzero .and. kp == 0) wta(:,:,1) = wtb(:,:,1)
   endif

end subroutine compute_pressures

!#######################################################################


subroutine compute_pres_weights ( Vgrid, phalf, pfull, wta, wtb )

 type(vert_grid_type), intent(in)  :: Vgrid
 real, intent(in)  :: phalf(:,:,:), pfull(:,:,:)
 real, intent(out) :: wta(:,:,:), wtb(:,:,:)

 real, dimension(size(phalf,1),size(phalf,2),size(phalf,3)) :: logph
 real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: logpf
 integer :: k, kp, kx, ks

   kp = Vgrid % nplev
   kx = size(pfull,3)

   if (Vgrid%pzero) then
       ks = max(2,kp+1)
   else
       ks = kp+1
   endif

   logph(:,:,ks  :kx+1) = log(phalf(:,:,ks  :kx+1))
   logpf(:,:,kp+1:kx  ) = log(pfull(:,:,kp+1:kx  ))

   do k = 1, kp
     wtb(:,:,k) = Vgrid % wtb(k)
   enddo
   do k = kp+1, kx
     wtb(:,:,k) = logph(:,:,k+1) - logpf(:,:,k)
   enddo

   if (Vgrid%pzero .and. kp == 0) wta(:,:,1) = wtb(:,:,1)
   do k = 1, kp
     wta(:,:,k) = Vgrid % wta(k)
   enddo
   do k = ks, kx
     wta(:,:,k) = logpf(:,:,k) - logph(:,:,k)
   enddo


end subroutine compute_pres_weights

!#######################################################################

subroutine compute_geop_height (Vgrid, fssl, vtemp, wta, wtb, &
                                zfull, zhalf, mask)

    type(vert_grid_type), intent(in)   :: Vgrid
    real, intent(in), dimension(:,:)   :: fssl
    real, intent(in), dimension(:,:,:) :: vtemp, wta, wtb
    real, intent(out)                  :: zfull(:,:,:)
    real, intent(out), optional        :: zhalf(:,:,:)
    real, intent(in),  optional        :: mask(:,:,:)

       integer :: k, klev
       real, dimension(size(vtemp,1),size(vtemp,2)) :: zb, zt, rt

!-----------------------------------------------------------------------
!    compute the height (in meters) at the interface between
!              model layers (half model levels)
!-----------------------------------------------------------------------

  klev = Vgrid % nlev

  if (size(zfull,3) /= klev) call error_mesg (  &
                               'compute_geop_height in vert_grid_mod', &
                               'incorrect dimension 3 for zfull', FATAL)

  if (present(zhalf)) then
     if (size(zhalf,3) /= klev+1) call error_mesg (  &
                               'compute_geop_height in vert_grid_mod', &
                               'incorrect dimension 3 for zhalf', FATAL)
  endif

      zb(:,:) = fssl(:,:)
      if (present(zhalf)) zhalf(:,:,klev+1) = zb(:,:)

!------- vertical integration loop (bottom to top) ----------

      do k = klev, 1, -1

         rt(:,:) = rdgas * vtemp(:,:,k)
         zt(:,:) = zb(:,:) + rt(:,:) * (wta(:,:,k)+wtb(:,:,k))
         zfull(:,:,k) = zb(:,:) + rt(:,:) * wtb(:,:,k)
         if (present(mask)) then
            where (mask(:,:,k) < 0.5)
                zt(:,:) = Vgrid % dfl(k)
                zfull(:,:,k) = 0.5*(Vgrid % dfl(k)+Vgrid % dfl(k+1))
            endwhere
         endif
         zb(:,:) = zt(:,:)
         if (present(zhalf)) zhalf(:,:,k) = zb(:,:)

      enddo

end subroutine compute_geop_height

!#######################################################################

subroutine compute_height (Vgrid, fssl, temp, sphum, pfull, phalf, &
                           zfull, zhalf, mask)

   type(vert_grid_type), intent(in)    :: Vgrid
   real, intent(in),  dimension(:,:)   :: fssl
   real, intent(in),  dimension(:,:,:) :: temp, sphum, pfull, phalf
   real, intent(out), dimension(:,:,:) :: zfull, zhalf
   real, intent(in),  optional         :: mask(:,:,:)

   real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  &
                                                   wta, wtb, vtemp

!-----------------------------------------------------------------------
!    compute the height (in meters) at the interface between
!              model layers (half model levels)
!    assumes that specific humidity will be used to compute
!                the virtual temperature
!-----------------------------------------------------------------------

  call compute_pres_weights ( Vgrid, phalf, pfull, wta, wtb )

  vtemp = temp * (1.0+d608*sphum)

  if (present(mask)) then
     call compute_geop_height ( Vgrid, fssl, vtemp, wta, wtb, &
                                zfull, zhalf, mask )
  else
     call compute_geop_height ( Vgrid, fssl, vtemp, wta, wtb, &
                                zfull, zhalf )
  endif

  zfull = zfull * ginv
  zhalf = zhalf * ginv

end subroutine compute_height

!#######################################################################

subroutine compute_height_bottom ( Vgrid, pssl, tbot, qbot,  &
                                   zbot, pbot, kbot )

   type(vert_grid_type), intent(in)  :: Vgrid
   real, intent(in),  dimension(:,:) :: pssl, tbot, qbot
   real, intent(out), dimension(:,:) :: zbot, pbot
   integer, intent(in),  optional    :: kbot(:,:)

   real, dimension(size(pssl,1),size(pssl,2)) :: rt, dp, phb, pht, &
                                                 lphb, lpht, lpf
   integer :: i, j, kb

!  ----- pressure at top and bottom interface of bottom level -----

   if (present(kbot)) then
      do j = 1, size(pssl,2)
      do i = 1, size(pssl,1)
         kb = kbot(i,j)
         pht(i,j) = Vgrid%peta(kb  ) + Vgrid%eta(kb  )*pssl(i,j)
         phb(i,j) = Vgrid%peta(kb+1) + Vgrid%eta(kb+1)*pssl(i,j)
      enddo
      enddo
   else
         kb = Vgrid%nlev
         pht(:,:) = Vgrid%peta(kb+1) + Vgrid%eta(kb+1)*pssl(:,:)
         phb(:,:) = Vgrid%peta(kb  ) + Vgrid%eta(kb  )*pssl(:,:)
   endif


   rt = ginv*rdgas * (tbot * (1.+d608*qbot))
   dp = phb - pht
   lphb = log(phb)
   lpht = log(pht)
   lpf = (phb*lphb-pht*lpht)/dp -1.
   
   zbot = rt * (lphb-lpf)
   pbot = exp(lpf)


end subroutine compute_height_bottom

!#######################################################################

end module scm_vert_grid_mod

