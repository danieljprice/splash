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
! Plot exact solution for a magnetohydrodynamic oblique C-shock
! (ie. one dimensional MHD shock problem with ambipolar diffusion)
!
! Solution from Mac-Low, Norman, Konigl & Wardle (1995), ApJ 442, 726
! http://adsabs.harvard.edu/abs/1995ApJ...442..726M
!
! ----------------------------------------------------------------------
module Cshock
 implicit none
 public :: exact_Cshock

contains

subroutine exact_Cshock(iplot,time,gamma,xmin,xmax,xpts,ypts,ierr)
 integer, intent(in) :: iplot
 integer, intent(out) :: ierr
 real, intent(in) :: time,gamma,xmin,xmax !,mach_sonic,mach_alfven,theta
 real, dimension(:), intent(inout) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 real, dimension(size(xpts)) :: D
 real, parameter :: pi = 3.1415926536
 real :: machs,macha,theta,xshock
 real :: rhoi,cs,rhon0,Bfield0,b0,shockl,vs
 integer :: npts,i
 
 npts = size(xpts)
 machs = 50. ! sonic Mach number
 macha = 5.  ! Alfvenic Mach number
 theta = pi/4.
 D(npts) = 1. + 1.e-4 ! upstream
 rhoi    = 1.e-5
 cs      = 0.1
 rhon0   = 1.
 Bfield0 = 1.
 b0      = sin(theta)
 shockl  = Bfield0/(gamma*rhoi*sqrt(rhon0))
 vs      = cs*machs
 xshock  = vs*time

 print "(4(a,f6.2))",&
  ' Plotting exact C-shock at t = ',time,' M = ',machs,' M_A = ',macha,' theta = ',theta
 print "(4(a,es10.3))",' shock length L = ',shockL,' shock is at x = ',xshock

 call integrate(xmin,xmax,xshock,xpts,macha,machs,theta,shockl,D,npts)

 !
 !--determine which parameter to plot
 !
 select case(iplot)
 case(1)
    ypts(1:npts) = D(1:npts)*rhon0  ! rho (neutrals)
    !print*,' D = ',D(1:npts)
 case(2)
    do i=1,npts   ! By = B_0*B
       ypts(i) = Bfield0*get_b(b0,macha,machs,D(i))
    enddo
 case(3) ! vx (neutrals)
    ypts(1:npts) = vs/D(1:npts)
 !case(4)
    !do i=1,npts   ! vy (neutrals)
    !   ypts(i) = Bfield0*get_b(b0,macha,machs,D(i))
    !enddo
 case(5)         ! Bx
    ypts(1:npts) = Bfield0*cos(theta)
 case default
    print*,'error: unknown solution to plot'
 end select
 ierr = 0

 return
end subroutine exact_Cshock

real function rhs(D,macha,machs,theta,shockl)
 real, intent(in) :: D,macha,machs,theta,shockl
 real :: term,sintheta,costheta2,b0,b
 
 term = (1./D**2 - 1./machs**2)*shockl
 
 sintheta = sin(theta)
 costheta2 = cos(theta)**2
 b0 = sintheta
 b = get_b(b0,macha,machs,D)
 
 rhs = b/macha*(b - D*((b - b0)/macha**2*costheta2 + sintheta))/(b**2 + costheta2)/term
 
end function rhs

real function get_b(b0,macha,machs,D)
 real, intent(in) :: b0,macha,machs,D

 get_b  = sqrt(b0**2 + 2.*macha**2*(D-1.)*(1./D - 1./machs**2))

end function get_b

subroutine integrate(xmin,xmax,xshock,xpts,macha,machs,theta,shockl,D,npts)
 integer, intent(in) :: npts
 real,    intent(in) :: xmin,xmax,xshock,macha,machs,theta,shockl
 real, dimension(npts), intent(inout) :: xpts
 real, dimension(npts), intent(inout) :: D
 real :: Dhalf,dx,xminshock
 integer :: i

 !
 ! set up grid to have good resolution around the shock
 !  and just two points extending to xmin and xmax
 !
 xminshock = min(xshock - 100.*shockl,xmin)
 dx = (xshock - xminshock)/real(npts-1)
 xpts(npts) = xmax
 xpts(1)    = xmin
 do i=2,npts-1
    xpts(i) = xmin + (i-1)*dx
 enddo

 do i=npts-1,1,-1
    if (xpts(i) > xshock) then
       D(i) = D(i+1)
    else
       !--mid-point rule
       Dhalf = D(i+1) - 0.5*dx*rhs(D(i+1),macha,machs,theta,-shockl)
       D(i)  = D(i+1) - 0.5*dx*rhs(Dhalf,macha,machs,theta,-shockl)
    endif
 enddo

end subroutine integrate

end module Cshock
