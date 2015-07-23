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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
! compute exact solution for Sedov-type point-like energy injection
! the solution is correct for 3D, I have attempted to fix it in 1D & 2D but
! not working yet.
!---------------------------------------------------------------------------
module sedov
  implicit none
  public :: exact_sedov
  private :: etau,rhou,pru,dudlneta,eta0

contains

subroutine exact_sedov(iplot,time,gam,rhozero,energy,rmax,rplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot
  integer, intent(out) :: ierr
  real, intent(in) :: time, gam, rhozero, energy, rmax
  real, dimension(:), intent(inout) :: rplot
  real, dimension(size(rplot)), intent(out) :: yplot

  integer, parameter :: ndim = 3
  integer :: npts
  real, parameter :: pi = 3.1415926536
  real :: rshock, dr
  real :: rhoshock, ushock, prshock
  real :: eta_0, power
  real :: ubar,ubarzero,dubar
  real :: phi,dphi
  integer :: i,ishock

  ierr = 0
  print*,' Plotting 3D Sedov similarity solution at t = ',time
  print*,' rhozero = ',rhozero,' energy = ',energy, ' rmax = ',rmax
  if (abs(time).lt.1.e-10) then
     print*,'nothing at t=0, returning'
     ierr = 1
     return
  endif

  npts = size(rplot)
  eta_0  = eta0(gam,ndim)
  print*,' eta0 = ',eta_0, ' gamma = ',gam

  power = 1./(ndim+2)
  dr = rmax/float(npts-1)
!
!--calculate radius and velocity of shock from dimensional analysis
!
  rshock = eta_0*(energy*time**2./rhozero)**power
!!  ushock = 2.*power*eta_0*((energy*time**2./rhozero)**(power-1.))*energy*time/rhozero
  ushock = 2.*power*rshock/time
  print*,' rshock = ',rshock, ' ushock = ',ushock
!
!--on x-y plots plot a circle at radius rshock
!
  if (iplot.eq.0) then
     dphi = 2.*pi/real(npts-1)
     phi = 0.
     do i=1,npts
        phi = (i-1)*dphi
        rplot(i) = rshock*cos(phi)
        yplot(i) = rshock*sin(phi)
     enddo
     return
  endif
!
!--jump conditions to find states behind shock
!
  rhoshock = rhozero*(gam+1.)/(gam-1.)
  prshock = 2./(gam+1.)*rhoshock*ushock**2

  rplot(1) = 0.0
  select case(iplot)
  case(2)
     yplot(1) = 0.0 ! shouldn't be zero for pressure
  case default
     yplot(1) = 0.0
  end select

  ishock = INT((rshock - rplot(1))/dr)

  if (ishock.gt.0) then
     ishock = min(ishock,npts)
!
!--the solution for rho is given as a function of ubar (dimensionless velocity)
!  ubar varies from (gamma+1)/2*gamma at eta = 0 to 1 at eta = 1
!  (eta is the dimensionless radius eta = r/rshock, so eta=1 is the shock position)
!
     ubarzero = 0.5*(gam+1.)/gam
     dubar = (1.0 - ubarzero)/REAL(ishock-1)  ! need to check the denominator
!
!--solution behind shock front (similarity solution)
!  (I really want to start from ubar = ubarzero, but am having problems)
!
     do i=2,ishock

        ubar = ubarzero + (i-1)*dubar ! again need to check this
        rplot(i) = etau(ubar,gam,ndim)*rshock

        select case(iplot)
        case(1)  ! rho
           yplot(i) = rhoshock*rhou(ubar,gam)
        case(2)  ! pr
           yplot(i) = prshock*(rplot(i)/rshock)**2*pru(ubar,gam)
        case(3)  ! utherm
           yplot(i) = prshock*(rplot(i)/rshock)**2*pru(ubar,gam)/ &
                ((gam-1.)*rhoshock*rhou(ubar,gam))
        case(4)  ! 1/2 v^2
           yplot(i) = 0.5*(4./(5.*(gam+1.))*rplot(i)/time*ubar)**2
        case(5)
           yplot(i) = (4./(5.*(gam+1.))*rplot(i)/time*ubar)
        end select
        !print*,'u,r, rho = ',ubar,rplot(i),rhoplot(i)

     enddo

  endif
!
!--solution ahead of shock front
!
  ishock = max(ishock,1)
  if (ishock.lt.npts) then
     do i=ishock,npts
        rplot(i) = rshock + (i-ishock)*dr
        select case(iplot)
        case(1)  ! rho
           yplot(i) = rhozero
        case default  ! pr,utherm,v
           yplot(i) = 0.
        end select
     enddo
  endif

  return
end subroutine exact_sedov

!
!--eta (dimensionless radius) as a function of u_bar (dimensionless velocity)
!
real function etau(u,gamma,ndim)
  implicit none
  integer, intent(in) :: ndim
  real, intent(in) :: u,gamma
  real :: gam1,term1,term2,power1,power2
  gam1 = gamma-1.
  power1 = (-12.+7.*gamma-13.*gamma**2)/(5.*(-1.+gamma+6.*gamma**2))
  power2 = gam1/(1.+2.*gamma)
  term1 = ((5.+5.*gamma+2.*u-6.*gamma*u)/(7.-gamma))**power1
  term2 = ((2.*gamma*u-gamma-1.)/gam1)**power2

  etau = u**(-2./(ndim+2.))*term1*term2
end function etau

!
!--rhobar (dimensionless density) as a function of u_bar (dimensionless velocity)
!
real function rhou(u,gamma)
  implicit none
  real, intent(in) :: u,gamma
  real :: gam1,power1,power2,power3,term1,term2,term3

  gam1 = gamma-1.
  power1 = (2./(gamma-2.))
  power2 = -(12.-7.*gamma+13.*gamma**2)/(2.-3.*gamma-11.*gamma**2+6.*gamma**3)
  power3 = 3./(1.+2.*gamma)
  term1 = ((1.+ gamma - 2.*u)/gam1)**power1
  term2 = ((5.+5.*gamma+2.*u-6.*gamma*u)/(7.-gamma))**power2
  term3 = ((2.*gamma*u - gamma -1.)/gam1)**power3

  rhou = term1*term2*term3

end function rhou
!
!--prbar (dimensionless pressure) as a function of u_bar (dimensionless velocity)
!
real function pru(u,gamma)
  implicit none
  real, intent(in) :: u,gamma

  pru = (gamma+1. - 2.*u)/(2.*gamma*u - gamma - 1.)*u**2*rhou(u,gamma)

end function pru
!
!--du /dln eta - required for the integral to compute eta0
!
real function dudlneta(u,gamma)
  implicit none
  real, intent(in) :: u,gamma
  real :: term1,term2

  term1 = u*(5.+5.*gamma + 2.*u - 6.*gamma*u)*(-1.-gamma+2.*gamma*u)
  term2 = 2.*(1.+gamma)*(1.+gamma - 2.*u - 2.*gamma*u + 2.*gamma*u**2)

  dudlneta = term1/term2

end function dudlneta

!
!--eta_0 as a function of gamma
!
real function eta0(gamma,ndim)
  implicit none
  integer, parameter :: ipts = 50000
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: ndim
  real, intent(in) :: gamma
  integer :: i
  real :: u0, u, du
  real :: sum, term, weight, factor

!  if (abs(gamma-5./3.).lt.1.e-3) then
!     eta0 = 1.1517
!  else
!     print*,'warning: don''t know eta0: integral not implemented'
!     eta0 = 1.0
!  endif

  u0 = 0.5*(gamma+1.)/gamma
  du = (1. - u0)/REAL(ipts)
!
!--integrate using Simpson's 1/3 rule
!
  sum = 0.
  do i=1,ipts
     weight = 1.0
     if (mod(i,2).eq.0) then
        weight = 4./3.
     else
        weight = 2./3.
     endif
     if ((i.eq.1).or.(i.eq.ipts)) weight = 1./3.
     u = u0 + i*du
     term = (pru(u,gamma) + rhou(u,gamma)*u**2)*(etau(u,gamma,ndim)**(ndim+2))/dudlneta(u,gamma)
     sum = sum + weight*du*term
  enddo

  if (ndim.eq.3) then
     factor = 4.*pi
  elseif (ndim.eq.2) then
     factor = 2.*pi
  else
     factor = 1.
  endif

  eta0 = (factor*8.*sum/(25.*(gamma**2 - 1.)))**(-1./REAL(ndim+2))

end function eta0

end module sedov
