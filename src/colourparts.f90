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

module colourparts
 implicit none

contains

subroutine colour_particles(dat,datmin,datmax,icolour,npart)
 use plotlib, only: plot_qcir
 implicit none
 integer, intent(in) :: npart
 real, dimension(npart), intent(in) :: dat
 real, intent(in) :: datmin, datmax
 integer, dimension(npart), intent(out) :: icolour
 integer :: i,icolourmin,icolourmax,icolourtemp
 real :: dx

 call plot_qcir(icolourmin,icolourmax)

 dx = (datmax - datmin)/real(icolourmax - icolourmin)

 do i=1,npart
    icolourtemp = int((dat(i) - datmin)/dx) + icolourmin
    if (icolourtemp.gt.icolourmax) icolourtemp = icolourmax
    if (icolourtemp.lt.icolourmin) icolourtemp = icolourmin
    icolour(i) = icolourtemp
 enddo

 return
end subroutine colour_particles

end module colourparts
