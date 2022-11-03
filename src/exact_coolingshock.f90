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

subroutine exact_coolingshock(iplot,time,gamma,xshock,rho_0,v_0,p_0,Lambda,xplot,yplot,ierr)
 use physcon, only:mH
 integer, intent(in) :: iplot
 integer, intent(out) :: ierr
 real, intent(in) :: time,gamma,xshock
 real, intent(in) :: rho_0,v_0,p_0,Lambda
 real, dimension(:), intent(out) :: xplot
 real, dimension(size(xplot)), intent(out) :: yplot
 integer, parameter :: nr = 100
 integer :: i
 real :: n_0,u_0,kB_T_on_Lambda_n0,vshock,r,r1,dr,a

 ierr = 0

 ! P = n*kB*T
 n_0 = rho_0/mH
 kB_T_on_Lambda_n0 = 0.1 !p_0/(n_0**2*Lambda)
 u_0 = p_0/((gamma-1.)*rho_0)
  print*,'kB*T/(lambda*n0) = ',kB_T_on_Lambda_n0

 vshock = 0.5*v_0 - sqrt((0.5*v_0)**2 + p_0/rho_0)  ! just below Eq A20 in C11
 a = rho_0*vshock**2/p_0

 print*,' vshock = ',vshock,' a = ',a

 ! define a grid in the similarity variable r1 = rho0/rho
 dr = 1.0/nr   ! r goes from 0 to 1
 do i=1,nr
    r1 = (i-0.5)*dr
    xplot(i) = xshock - vshock/(gamma-1.)*(kB_T_on_Lambda_n0)*&
              ( (gamma-a)/(a-1.)*log(1 - r1) &
              + (1.-a*gamma)/((a-1.)*a**2)*log(1. - a*r1) &
              - (a+1.)/a*r1 - 0.5*(gamma+1.)*r1**2 )
    select case(iplot)
    case default
       yplot(i) = u_0*((a+1.)*r1 - a*r1*r1)  ! u/u_0
       print*,r1,' x = ',xplot(i),' y = ',yplot(i)
    end select
 enddo

end subroutine exact_coolingshock

end module coolingshock
