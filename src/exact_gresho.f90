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
!  Copyright (C) 2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
! compute exact solution for gresho vortex problem
! ----------------------------------------------------------------------
module gresho
  implicit none

contains

subroutine exact_gresho(iplot,xplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot
  real, intent(in),  dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr

  print*,'plotting gresho vortex '
!
! check for errors
!
  ierr = 0
  select case(iplot)  
  case(2) ! pressure
     where (xplot < 0.2)
        yplot = 5. + 12.5*xplot**2
     elsewhere (xplot < 0.4)
        yplot = 9. + 12.5*xplot**2 - 20.*xplot + 4.*log(5.*xplot)
     elsewhere
        yplot = 3. + 4.*log(2.)
     end where
  case(1) ! vphi
     where (xplot < 0.2)
        yplot = 5.*xplot
     elsewhere (xplot < 0.4)
        yplot = 2. - 5.*xplot
     elsewhere
        yplot = 0.
     end where
  end select

  return
end subroutine exact_gresho

end module gresho
