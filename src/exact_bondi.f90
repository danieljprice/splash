!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
! Exact solution for Bondi flow
! ----------------------------------------------------------------------
module bondi
 implicit none

 public :: exact_bondi

 real, private, parameter :: pi = 4.*atan(1.)

 ! Constants from user input
 real, private :: rcrit, rhocrit     ! for non-rel
 real, private :: den0, en0          ! for GR geodesic flow
 real, private :: adiabat            ! for GR sonic-point flow (and rcrit as well)

 real, private :: C1,C2,Tc,n         ! for GR sonic-point flow: intermediate constants not input by user
 real, private :: mass1
 logical, private :: iswind

 private

contains

subroutine exact_bondi(iplot,time,gamma,const1,const2,m,relativistic,geodesic_flow,is_wind,xpts,ypts,ierr)
 integer, intent(in)  :: iplot
 integer, intent(out) :: ierr
 real, intent(in)     :: time,gamma,const1,const2,m
 logical, intent(in)  :: relativistic, geodesic_flow,is_wind
 real, dimension(:), intent(in) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 integer :: i,npts
 real    :: r,rhor,vr,ur

 npts = size(xpts)

 print "(a,es10.3)",' Plotting exact Bondi solution at t = ',time

 if (.not.relativistic) then
    rcrit   = const1
    rhocrit = const2
 elseif (relativistic) then
    if (geodesic_flow) then
       den0 = const1
       en0  = const2
    elseif (.not.geodesic_flow) then
       rcrit   = const1
       adiabat = const2
       iswind  = is_wind
    endif
 endif
 !
 !--determine which parameter to plot
 !
 do i=1,npts
    r = xpts(i)

    vr = 0.
    ur = 0.
    rhor = 0.

    ! Note: Lambert functions solutions not great below 0.3 for some rcrit and rhocrit
    if (.not. relativistic .and. r>0.3) then
       call get_bondi_nonrel(rhor,vr,ur,r,m,gamma)
    elseif (relativistic .and. r>2.) then
       if (geodesic_flow) then
          call get_bondi_geodesic(rhor,vr,ur,r,m,gamma)
       elseif (.not. geodesic_flow) then
          call get_bondi_sonicpoint(rhor,vr,ur,r,m,gamma)
       endif
    endif

    select case(iplot)
    case(1)
       ypts(i) = vr
       if (.not.is_wind) ypts(i) = -vr
    case(2)
       ypts(i) = ur
    case(3)
       ypts(i) = rhor
    case default
       ypts(i) = 0.
    end select
 enddo
 ierr = 0

 return
end subroutine exact_bondi


!------------------------------------------------------------------------
!--- Non-relativistic solution ------------------------------------------
!------------------------------------------------------------------------

subroutine get_bondi_nonrel(rho,v,u,r,mass,gamma)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass,gamma
 real :: cs2,vr,mdot

 cs2  = mass/(2.*rcrit)
 mdot = rhocrit*4.*pi*rcrit**2*sqrt(cs2)

 if (r>=rcrit) then
    vr = sqrt(-cs2*lambertw_0(-func(r,rcrit)))
 else
    vr = sqrt(-cs2*lambertw_neg1(-func(r,rcrit)))
 endif

 v   = vr
 rho = mdot/(4.*pi*abs(v)*r**2)
 u   = 0.

end subroutine get_bondi_nonrel

! See Barry, Parlange & Li 2000, for the analytic approximations used for the Lambert W function.
!
!--- Lambert W function for the principal branch (k=0) approaching from the negative
real function lambertw_0(x)
 real, intent(in) :: x
 real, parameter  :: exp1 = exp(1.)
 real :: eta, N1, N2
 eta = 2. + 2.*exp1*x
 N2  = 6. + 3.*Sqrt(2.) - ((-5764. - 4108.*Sqrt(2.) + (2237. + 1457.*Sqrt(2.))*exp1)*eta)/&
      (-796. - 430.*Sqrt(2.) + (215. + 199.*Sqrt(2.))*exp1)
 N1  = (1. - 1./Sqrt(2.))*(Sqrt(2.) + N2)
 lambertw_0 = -1 + Sqrt(eta)/(1 + (Sqrt(eta)*N1)/(Sqrt(eta) + N2))
end function lambertw_0

!--- Lambert W function for the k=-1 branch
real function lambertw_neg1(x)
 real, intent(in) :: x
 real, parameter  :: M1 = 0.3361, M2 = -0.0042, M3 = -0.0201
 real :: sigma
 sigma = -1. - Log(-x)
 lambertw_neg1 = -1. - sigma - (2.*(1. - 1./(1. + (M1*Sqrt(sigma)*(1. + exp(M3*Sqrt(sigma))*M2*sigma))/Sqrt(2.))))/M1
end function lambertw_neg1

!--- Function used in the non-rel solution of velocity [D(r)] (See: Cranmer 2004)
real function func(r,rc)
 real, intent(in) :: r,rc
 func = (rc/r)**4 * exp(4.*(1.-rc/r) - 1.)
end function func


!------------------------------------------------------------------------
!--- GR geodesic flow solution ------------------------------------------
!------------------------------------------------------------------------

subroutine get_bondi_geodesic(rho,v,u,r,m,gamma)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,m,gamma
 real :: sqrtg,alpha,dfunc,efunc

 dfunc = den0/(r**2*sqrt(2.*m/r*(1.- 2.*m/r)))
 efunc = en0/((sqrt(2.*m/r)*r**2)**gamma * sqrt(1.- 2.*m/r))

 sqrtg = 1.
 alpha = sqrt(1. - 2.*m/r)
 rho = sqrtg/alpha*dfunc
 v   = sqrt(2.*m/r)*(1. - 2.*m/r)
 u   = efunc/dfunc

end subroutine get_bondi_geodesic


!------------------------------------------------------------------------
!--- GR sonic point flow solution ---------------------------------------
!------------------------------------------------------------------------
subroutine get_bondi_sonicpoint(rho,v,u,r,m,gamma)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,m,gamma
 real :: T,uvel,term,u0,dens,sqrtg
 real :: uc2,vc2

 n     = 1./(gamma-1.)
 mass1 = m

 uc2 = mass1/(2.*rcrit)
 vc2 = uc2/(1.-3.*uc2)

 Tc  = vc2*n/(1.+n-vc2*n*(1.+n))
 C1  = sqrt(uc2) * Tc**n * rcrit**2
 C2  = (1. + (1.+n)*Tc)**2 * (1. - 2.*mass1/rcrit + C1**2/(rcrit**4*Tc**(2.*n)))

 ! Given an r, solve eq 76 for T numerically (Hawley, Smarr and Wilson 1976a)
 call Tsolve(T,r)

 uvel = C1/(r**2 * T**n)
 dens = adiabat*T**n
 u = T*n

 !get u0 at r
 term = 1./(1.-2.*mass1/r)
 u0  = sqrt(term*(1.+term*uvel**2))
 v   = uvel/u0

 sqrtg = 1.
 rho = sqrtg*u0*dens

end subroutine get_bondi_sonicpoint

! Newton Raphson
subroutine Tsolve(T,r)
 real, intent(in)  :: r
 real, intent(out) :: T
 real    :: Tnew, diff
 logical :: converged
 integer :: its
 integer, parameter :: itsmax = 100
 real, parameter    :: tol    = 1.e-5

 ! These guess values may need to be adjusted for values of rcrit other than rcrit=8M
 if ((iswind .and. r>=rcrit) .or. (.not.iswind .and. r<rcrit)) then
    T = 0.760326*r**(-1.307)/2.   ! This guess is calibrated for rcrit=8M, and works ok up to r ~ 10^7 M
 elseif ((iswind .and. r<rcrit) .or. (.not.iswind .and. r>=rcrit)) then
    T = 100.
 endif

 converged = .false.
 its = 0
 do while (.not.converged .and. its<itsmax)
    Tnew = T - ffunc(T,r)/df(T,r)
    diff = abs(Tnew - T)/abs(T)
    converged = diff < tol
    T = Tnew
    its = its+1
 enddo

 if (.not.converged) print*,'Bondi exact solution not converged at r = ',r

end subroutine Tsolve

real function ffunc(T,r)
 real, intent(in) :: T,r
 ffunc = (1. + (1. + n)*T)**2*(1. - (2.*mass1)/r + C1**2/(r**4*T**(2.*n))) - C2
end function ffunc

real function df(T,r)
 real, intent(in) :: T,r
 df = (2.*(1. + T + n*T)*((1. + n)*r**3*(-2.*mass1 + r) - C1**2*T**(-1. - 2.*n)*(n + (-1. + n**2)*T)))/r**4
end function df

end module bondi
