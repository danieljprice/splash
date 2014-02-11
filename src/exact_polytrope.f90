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

 !---------------------------------------
 ! numerically integrate a polytrope
 ! with the density = 1 at the centre
 ! This uses scaled variables (no sigma in the equation)
 ! then the radius is scaled to give the correct mass
 ! Based on an old subroutine from Joe Monaghan
 !-----------------------------------------

module polytrope
 implicit none
 public :: exact_polytrope

contains

subroutine exact_polytrope(gamma,polyk,totmass,rplot,denplot,npts,ierr)
  implicit none
  integer, intent(out) :: npts,ierr
  real, intent(in)     :: gamma
  real, intent(in)     :: polyk,totmass
  real, dimension(:), intent(inout) :: rplot
  real, dimension(size(rplot)), intent(out) :: denplot

  integer :: i,j
  real, parameter :: pi = 3.1415926536
  real, dimension(size(rplot)) :: r,v,den
  real :: dr,an,rhs,rstar,totmassf
  real :: rhocentre,fac,rfac,G
  
  ierr = 0
  print*,' gamma           :',gamma
  dr = 0.001
  G = 1.
  an = 1./(gamma-1.)
  v(1) = 0.0
  v(2) = dr*(1.0 - dr*dr/6. )
  r(1) = 0.
 
  i = 2
  do while (v(i).ge.0.)
    r(i) = (i-1)*dr
    rhs = - r(i)*(v(i)/r(i))**an
    v(i+1) = 2*v(i) - v(i-1) + dr*dr*rhs
    i = i + 1
    if (i+1.gt.size(rplot)) then
       dr = dr*2.
       r(2) = dr
       v(2) = dr*(1.0 - dr*dr/6. )
       i = 2
    endif
  enddo
  npts = i-1
  rstar = r(npts)

  !--------------------------------------
  ! calculate the mass out to radius r
  ! using the density without the central
  ! density multiplier- call this totmassf
  ! the true scaled totmass = 1.
  !----------------------------------------

  den(1) = 1.0
  totmassf = 0.
  do j = 2,npts
    den(j) = (v(j)/r(j))**an
    totmassf = totmassf + 4.*pi*r(j)*r(j)*den(j)*dr
  enddo

  !---------------------------------------------------
  ! rescale the central density to give desired massq
  ! then rescale the radius to match this
  !---------------------------------------------------
  fac = (gamma*polyk)/(4.*pi*G*(gamma - 1.))
  rhocentre = ((totmass/totmassf)/fac**1.5)**(2./(3.*gamma - 4.))
  rfac = sqrt(fac*rhocentre**(gamma - 2.))
  
  print*,' Rstar = ',rstar*rfac
  print*,' central density :',rhocentre
  print*,' total mass      :',totmass

  rplot = r * rfac
  denplot = rhocentre * den

  return
end subroutine exact_polytrope

end module polytrope
