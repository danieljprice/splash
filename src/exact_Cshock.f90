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

subroutine exact_Cshock(iplot,time,gamma,machs,macha,xmin,xmax,xpts,ypts,ierr)
 integer, intent(in) :: iplot
 integer, intent(out) :: ierr
 real, intent(in) :: time,gamma,machs,macha,xmin,xmax
 real, dimension(:), intent(inout) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 real, dimension(size(xpts)) :: D
 real, parameter :: pi = 3.1415926536
 real :: theta,xshock,ambi_gamma,ambi_rhoi,vx0,vy0,Bx,By,rhon_pre
 real :: cs,rhon0,Bfield0,b0,shockl,vs,va,K1,K2,By_post,P_post,rhon_post,By0,Pr0
 real :: sintheta,costheta,vx2,vx,vy,rhon,dvx,vxin
 integer :: npts,i
 logical printout
 
 printout = .false.
 npts = size(xpts)
 theta = pi/4.
 costheta = cos(theta)
 sintheta = sin(theta)
 D(npts) = 1. + 1.e-6 ! upstream
 cs      = 0.1
 rhon0   = 1.
 Bfield0 = 1. ! this gives Bx0 = By0 = 1/sqrt(2) as in Choi et al. (2009)
 ambi_gamma = 1.
 ambi_rhoi  = 1.e-5
 b0      = sintheta
 shockl  = Bfield0/(ambi_gamma*ambi_rhoi*sqrt(rhon0))
 va      = Bfield0/sqrt(rhon0)
 xshock  = 6./8.*va*time

 if ( printout ) open(unit = 625,file="Cshock_splash.dat")


 print "(4(a,g8.2))",&
  ' Plotting exact C-shock at t = ',time,' M = ',machs,' M_A = ',macha,' theta = ',theta
 print "(4(a,es10.3))",' shock length L = ',shockL,' shock is at x = ',xshock

 call integrate(xmin,xmax,xshock,xpts,macha,machs,theta,shockl,D,npts)
 
 !
 !  compute velocity jump across shock: See Mac-Low et al. (1995). This is the difference 
 !  in the velocity across the shock front since we assume that the
 !  post-shock gas is at rest
 !
! !--post-shock, assume vx = 0
 By_post   = Bfield0*get_b(b0,macha,machs,D(1))
 rhon_post = D(1)*rhon0
 P_post    = rhon_post*cs**2
! K1 = P_post + 0.5*By_post*By_post
! print*,' K1 is ',K1

 !--pre-shock
 vx0  = -5.0
 vxin = -4.45
 dvx  = vxin-vx0
 vy0 = 0.
 Bx       = Bfield0*costheta
 By0      = Bfield0*get_b(b0,macha,machs,D(npts))
 rhon_pre = D(npts)*rhon0
 Pr0      = rhon_pre*cs**2
 K1 = Pr0 + 0.5*By0*By0 + rhon_pre*vx0**2
 K2 = rhon_pre*vx0*vy0 - Bx*By0

 vx2 = (K1 - 0.5*By_post**2 - P_post)/rhon_post
 if (vx2 > 0.) then
    vx = -sqrt(vx2)
    print "(1x,a,g10.3)",'vx post-shock = ',vx
 else
    vx = 0.
    print*,'error, post-shock vx is imaginary'
 endif
 vs = vx0 - vx
 !print*,'vs = ',vs
 !
 !--determine which parameter to plot
 !
 do i=1,npts
    rhon = D(i)*rhon0
    By   = Bfield0*get_b(b0,macha,machs,D(i))
    vx2  = (K1 - 0.5*By**2 - rhon*cs**2)/rhon
    if (vx2 > epsilon(vx2)) then
       vx = -sqrt(vx2)
       vy = (K2 + Bx*By)/(rhon*vx)
    else
       vx = 0.
       vy = 0.
    endif
    vx = vx + dvx

    select case(iplot)
    case(1)
       ypts(i) = rhon ! rho (neutrals)
    case(2)
       ypts(i) = By
    case(3) ! vx (neutrals)
       ypts(i) = vx
    case(4) ! vy (neutrals)
       ypts(i) = vy
    case(5) ! Bx
       ypts(i) = Bx
    case default
       ypts(i) = 0.
    end select
    if ( printout ) write(625,*) i,xpts(i),rhon,Bx,By,vx,vy
 enddo
 ierr = 0

 if ( printout ) close(625)
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
       D(i)  = D(i+1) - dx*rhs(Dhalf,macha,machs,theta,-shockl)
    endif
 enddo

end subroutine integrate

end module Cshock
