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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
!  Plots various analytic density profiles
!
!  Currently implemented:
!
!   1) Plummer sphere
!   2) Hernquist sphere
!
!  Added by D. Price 5/9/05
! ----------------------------------------------------------------------
module densityprofiles
  implicit none

contains

subroutine exact_densityprofiles(iplot,iprofile,Msphere,rsoft,xplot,yplot,ierr)
  implicit none
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: iplot,iprofile
  real, intent(in), dimension(2) :: Msphere,rsoft
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  integer :: i
!
! check for errors
!
  ierr = 0
  if (all(Msphere.le.0.)) then
     print*,'error: mass <= 0 in exact_densityprofile'
     ierr = 2
     return
  endif
  if (any(rsoft.lt.0.)) then
     print*,'error: rsoft < 0 in exact_densityprofile'
     ierr = 3
     return
  endif

  select case(iprofile)
  case(1)
!
!--Plummer sphere
!
    select case(iplot)
    case(2)
!--potential
       do i=1,size(xplot)
          yplot(i) = -Msphere(1)/sqrt(rsoft(1)**2 + xplot(i)**2) &
                     -Msphere(2)/sqrt(rsoft(2)**2 + xplot(i)**2)
       enddo
!--force
    case(3)
       do i=1,size(xplot)
          yplot(i) = Msphere(1)*xplot(i)/((rsoft(1)**2 + xplot(i)**2)**1.5) &
                   + Msphere(2)*xplot(i)/((rsoft(2)**2 + xplot(i)**2)**1.5)
       enddo
!--density
    case default
       do i=1,size(xplot)
          yplot(i) = 3.*Msphere(1)*rsoft(1)**2/(4.*pi*(rsoft(1)**2 + xplot(i)**2)**2.5) &
                   + 3.*Msphere(2)*rsoft(2)**2/(4.*pi*(rsoft(2)**2 + xplot(i)**2)**2.5)
       enddo
    end select

  case(2)
!
!--Hernquist model (use tiny to prevent divergences in cusp)
!
    select case(iplot)
    case(2)
!--potential
       do i=1,size(xplot)
          yplot(i) = -Msphere(1)/(rsoft(1) + xplot(i)) &
                     -Msphere(2)/(rsoft(2) + xplot(i))
       enddo
!--force
    case(3)
       do i=1,size(xplot)
          yplot(i) = Msphere(1)/(rsoft(1) + xplot(i))**2 &
                   + Msphere(2)/(rsoft(2) + xplot(i))**2
       enddo
!--density
    case default
       do i=1,size(xplot)
          yplot(i) = Msphere(1)*rsoft(1)/ &
                     ((2.*pi*xplot(i)*(rsoft(1) + xplot(i))**3)) &
                   + Msphere(2)*rsoft(2)/ &
                     ((2.*pi*xplot(i)*(rsoft(2) + xplot(i))**3))
       enddo
    end select

  case default
     ierr = 1
  end select

  return
end subroutine exact_densityprofiles

end module densityprofiles
