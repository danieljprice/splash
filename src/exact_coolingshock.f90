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
!  Copyright (C) 2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
! Plot exact solution for cooling shock problem
! (ie. one dimensional shock problem with prescribed cooling function)
!
! Taken from Creasey et al. (2011), MNRAS 415, 3706-3720
! http://adsabs.harvard.edu/abs/2011MNRAS.415.3706C
!
! ----------------------------------------------------------------------
module coolingshock
 implicit none
 public :: exact_coolingshock

contains

subroutine exact_coolingshock(iplot,time,gamma,xshock,rho_0,v_0,p_0,cooling_length,xplot,yplot,ierr)
 use shock, only:get_pstar
 integer, intent(in) :: iplot
 integer, intent(out) :: ierr
 real, intent(in) :: time,gamma,xshock
 real, intent(in) :: rho_0,v_0,p_0,cooling_length
 real, dimension(:), intent(out) :: xplot
 real, dimension(size(xplot)), intent(out) :: yplot
    !integer, parameter :: nr = 100
 integer :: i,nr
 real :: n_0,u_0,cs0,pstar,vshock,r,r1,dr,a
 real :: vc,rhoc,fac,umass,mH_code,kB_T_on_Lambda_n0
 real :: lambda_min,lambda_max,r1_min,r1_max

 ierr = 0
 print*,' v_0,rho_0,P_0 = ',v_0,rho_0,p_0,' cooling length = ',cooling_length

 u_0 = p_0/((gamma-1.)*rho_0)
 cs0 = sqrt(gamma*p_0/rho_0)   ! upstream sound speed
 kB_T_on_Lambda_n0 = cooling_length/cs0   ! Eq 17 in Creasey et al
 print*,'kB*T/(lambda*n0) = ',kB_T_on_Lambda_n0,' u0 = ',u_0,' gamma = ',gamma

 vc = 0.
 !call riemannsolver()
 !call get_pstar(gamma,p_0,p_0,v_0,-v_0,cs0,cs0,pstar,vshock)
 print*,' got vshock = ',vshock
 vshock = 0.5*v_0 + sqrt((0.5*v_0)**2 + p_0/rho_0)  ! just below Eq A20 in C11
 a = rho_0*(v_0 - vshock)**2/p_0   ! square of the Mach number
 rhoc = -(v_0 - vshock)*rho_0/vshock

 fac = (vc - vshock)/(gamma-1.)*(kB_T_on_Lambda_n0) !*rho_0/rhoc
 !fac = (vc - vshock)/(Aterm*rhoc) 
 print*,' vshock = ',vshock,' rho c = ',rhoc,' a = ',a

 nr = size(xplot)
 ! define a grid in the similarity variable r1 = rho0/rho = 1./lambda
 lambda_min = a/(a + 1.)*(gamma+1)/gamma
 lambda_max = 1.
 print*,' lambda_min,lambda_max=',lambda_min,lambda_max
 r1_min = 1./lambda_max
 r1_max = 1./lambda_min
 dr = (r1_max - r1_min)/(nr-1)   ! r goes from 0 to 1
 do i=1,nr
    r1 = r1_min + (i-0.5)*dr
    xplot(i) = xshock - fac * &
              ( (gamma-a)/(a-1.)*log(r1 - 1.) &   ! note typo in formula in their appendix
              + (1.-a*gamma)/((a-1.)*a**2)*log(1. - a*r1) &
              - (a+1.)/a*r1 - 0.5*(gamma+1.)*r1**2 )
    select case(iplot)
    case(3) ! pressure
       yplot(i) = (gamma-1.)*(rho_0/r1)*(u_0*((a+1.)*r1 - a*r1*r1)) ! P = (gamma-1) rho u
    case(2) ! u
       yplot(i) = u_0*((a+1.)*r1 - a*r1*r1) ! from u/u0 (Eq A8 in C11)
    case(1) ! v
       yplot(i) = vshock + (vc - vshock)*r1 ! from Eq A7 in C11
    case default  ! rho
       yplot(i) = rhoc/r1 ! from rho/rho0 = lambda (Eq A6 in C11)
    end select
 enddo
 print*,' xmin,xmax = ',minval(xplot),maxval(xplot)

end subroutine exact_coolingshock

end module coolingshock
