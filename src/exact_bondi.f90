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

contains

subroutine exact_bondi(iplot,time,gamma,r0,rho0,m,relativistic,xpts,ypts,ierr)
 integer, intent(in)  :: iplot
 integer, intent(out) :: ierr
 real, intent(in)     :: time,gamma,r0,rho0,m
 logical, intent(in)  :: relativistic
 real, dimension(:), intent(in) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 real, parameter :: pi = 3.1415926536
 integer :: i,npts
 real :: r,vr,ur,e0,d0,den,en,cs2,rhor,mdot

 npts = size(xpts)

 print "(4(a,g8.2))",' Plotting exact Bondi solution at t = ',time

 !
 !--determine which parameter to plot
 !
 do i=1,npts
    r = xpts(i)

    vr = 0.
    ur = 0.
    rhor = 0.

    if (.not. relativistic) then

       cs2  = m/(2.*r0)
       mdot = rho0*4.*pi*r0**2*sqrt(cs2)

       if (r>=r0) then
          vr = -sqrt(-cs2*lambertw_0(-func(r,r0)))
       else
          vr = -sqrt(-cs2*lambertw_neg1(-func(r,r0)))
       endif

       rhor = mdot/(4.*pi*abs(vr)*r**2)

       ur = cs2/(gamma-1.)

       if (r<=tiny(r)) then
          vr = 0.
          ur = 0.
          rhor = 0.
       endif

    elseif (relativistic) then
      !  e0  = 0.000297118
      !  d0  = 12.
       d0 = 1.
       e0 = 1.e-9
       en = e0/((sqrt(2.*m/r)*r**2)**gamma * (1.- 2.*m/r)**((gamma + 1.)/4.))
       den = d0/(r**2*sqrt(2.*m/r*(1.- 2.*m/r)))
       ur  = en/den

       if (r0==0.) then ! Marginally bound case (r0=infinity)
          vr = -(1. - 2.*m/r)*sqrt(2.*m/r)
       elseif(r<=r0) then
          vr = -(1. - 2.*m/r)/sqrt(1.-2.*m/r0)*sqrt(2.*m*(1./r - 1./r0))
       endif

       if (r<2.) then
          vr = 0.
          ur = 0.
          rhor = 0.
       endif

     endif

    select case(iplot)
    case(1)
       ypts(i) = vr
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


! Lambert W function for the principal branch (k=0) approaching from the negative
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

! Lambert W function for the k=-1 branch
real function lambertw_neg1(x)
 real, intent(in) :: x
 real, parameter  :: M1 = 0.3361, M2 = -0.0042, M3 = -0.0201
 real :: sigma

 sigma = -1. - Log(-x)
 lambertw_neg1 = -1. - sigma - (2.*(1. - 1./(1. + (M1*Sqrt(sigma)*(1. + exp(M3*Sqrt(sigma))*M2*sigma))/Sqrt(2.))))/M1

end function lambertw_neg1

real function func(r,rc)
 real, intent(in) :: r,rc

 func = (rc/r)**4 * exp(4.*(1.-rc/r) - 1.)

end function func

end module bondi
