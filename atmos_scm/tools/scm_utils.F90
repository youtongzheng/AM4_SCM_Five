module scm_utils_mod

   use sat_vapor_pres_mod, only:  lookup_es
   use      constants_mod, only:  rdgas, rvgas, cp_air, hlv, kappa, grav
   real,    private  :: p00 = 100000.

   public us_std_atm, thetal_to_temp, locate

contains

!########################################################################
! Subroutine returns US standard atmosphere temperature and specific 
! humidity for a given pressure.
!
! Reference: http://mtp.jpl.nasa.gov/notes/altitude/StdAtmos1976.html
!
       subroutine us_std_atm( p, T, qv )
       implicit none

       real, intent(in)  :: p
       real, intent(out) :: T, qv

       integer i
       real dz

!      US standard atmosphere data:
      
       integer, parameter :: nmax = 8
       real, dimension(nmax) :: Zs, Ts, LRs, ps

!      Geopotential height (km)
       data Zs / 0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852 /
!      Temperature (K)
       data Ts / 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95 /
!      Lapse rate (K/km)
       data LRs / -6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0 /
!      Pressure (Pa)
       data ps / 1013.25e2, 226.3206e2, 54.7489e2, 8.6802e2, 1.1091e2, 0.6694e2, 0.0396e2, 0.0037e2 /

!      Identify reference level

       i = nmax
       do while ( p >= ps(i) .and. i >= 2 )
         i = i - 1
       enddo
        
       if ( i < nmax ) then

!        Compute height above reference level
         if ( LRs(i) == 0.0 ) then
           dz = -  rdgas * Ts(i) / ( 1000.0 * grav ) * log( p / ps(i) )
         else
           dz = - (Ts(i)/LRs(i)) * ( 1. - (p/ps(i))**(-rdgas*LRs(i)/(1000.0*grav)) )
         endif

!        Compute interpolated temperature
         T = Ts(i) + LRs(i) * dz

       else

         T = Ts(i)

       endif

!      Water vapor is always zero

       qv = 0.0

       return
       end subroutine us_std_atm

!########################################################################
! Subroutines to convert thetal to temperature

      subroutine thetal_to_temp(thetal_in,rt_in,pres_in,temp,qv,qc)
      implicit none

      real, intent(in)  :: thetal_in, rt_in, pres_in
      real, intent(out) :: temp, qv, qc

!     Internal variables

      real ex
      real esat, rsat
      real temp_unsat, temp_sat, temp1, temp2

!     Conversion from theta_l to temperature

      ex = ( pres_in / p00 )**(kappa)

!     unsaturated case

      temp_unsat = thetal_in * ex

!     return right away if there is no moisture

      if ( rt_in <= 1.e-10 ) then
        temp = temp_unsat
        qv = rt_in
        qc = 0.0
        return
      endif

!     saturated case

      temp1 = temp_unsat - 30.0
      temp2 = temp_unsat + 30.0
      temp_sat = rtsec(temp1,temp2,thetal_in,rt_in,pres_in,ex,1.e-3)

!     Choose between one of the two choices, the highest value is temp
!     Compute the water vapor and cloud liquid water specific humidities

      if (temp_sat .gt. temp_unsat ) then
        temp = temp_sat
        call lookup_es(temp,esat)
        rsat = rdgas/rvgas * esat / ( pres_in - esat )
        qv   =          rsat / ( 1.0 + rt_in )
        qc   = ( rt_in - rsat ) / ( 1.0 + rt_in )
      else
        temp = temp_unsat
        qv   = rt_in / ( 1.0 + rt_in )
        qc   = 0.0
      endif

      return
      end subroutine thetal_to_temp

!     -----------------------------------------------------------------
!     To diagnose temp from thetal and rt, we need to find the zero
!     of this function

      real function saturated(temp,thetal,rt,pres,ex)
      implicit none

      real, intent(in) :: temp, thetal, rt, pres, ex

      real esat, rsat

      call lookup_es(temp,esat)
      rsat = rdgas/rvgas * esat / ( pres - esat )

      saturated = thetal*ex - temp + hlv/cp_air * ( rt - rsat )

      return
      end function saturated

!     -----------------------------------------------------------------
!     Function to find zero of function 'saturated' using the secant 
!     method

      real function rtsec(x1,x2,ya,yb,yc,yd,xacc)
      implicit none

      integer maxit
      parameter (maxit=30)

      real x1,x2,xacc
      real ya,yb,yc,yd

      integer j
      real dx,f,fl,swap,xl

      fl = saturated(x1,ya,yb,yc,yd)
      f = saturated(x2,ya,yb,yc,yd)

      if (abs(fl).lt.abs(f)) then
         rtsec = x1
         xl = x2
         swap = fl
         fl = f
         f = swap
      else
         xl = x1
         rtsec = x2
      endif
      do j=1,maxit
          dx = (xl-rtsec)*f/(f-fl)
          xl = rtsec
          fl = f
          rtsec = rtsec + dx
          f = saturated(rtsec,ya,yb,yc,yd)
          if (abs(dx).lt.xacc.or.f.eq.0.) return
      enddo

      write(*,*) 'Warning: rtsec exceed maximum iterations'
      write(*,*) 'thetal = ',ya
      write(*,*) 'qt     = ',yb
      write(*,*) 'pres   = ',yc
      write(*,*) 'ex     = ',yd
      write(*,*) 'x1     = ',x1
      write(*,*) 'x2     = ',x2
      write(*,*) 'xacc   = ',xacc
!      call stopcode('rtsec: exceed maximum iterations')
      stop 'rtsec: exceed maximum iterations'

      return
      end function rtsec

!########################################################################
! Subroutines to find j such that x falls between xx(j) and xx(j+1)

      subroutine locate(xx,n,x,j)

      integer, intent(in)  :: n
      real,    intent(in)  :: xx(n), x
      integer, intent(out) :: j

      integer jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl

      return
      end subroutine locate

!########################################################################

end module scm_utils_mod

