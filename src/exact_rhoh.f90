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

!----------------------------------------------------------
! plots the relation between smoothing length and density
!
! ie. h = h_fact*(pmass/rho)^(1/ndim)
!
!----------------------------------------------------------
module rhoh
 implicit none
 public :: exact_rhoh

contains

subroutine exact_rhoh(iplot,ndim,hfact,pmassval,xplot,yplot,ierr)
 implicit none
 integer, intent(in) :: iplot,ndim
 integer, intent(out) :: ierr
 real, intent(in) :: hfact,pmassval
 real, dimension(:), intent(in) :: xplot
 real, dimension(size(xplot)), intent(out) :: yplot

 if (hfact.gt.0.01) then
    ierr = 0

    if (iplot.eq.2) then ! x axis is h
       where (xplot > tiny(xplot))
         yplot(:) = pmassval*(hfact/xplot(:))**ndim
       elsewhere
         yplot(:) = huge(yplot)
       end where
    else ! y axis is h
       where (xplot > tiny(xplot))
         yplot(:) = hfact*(pmassval/xplot(:))**(1./FLOAT(ndim))
       elsewhere
         yplot(:) = huge(yplot)
       end where
    endif
    write(*,"(a,f5.2,a,es9.2,a,i1,a)") ' plotting h = ',hfact, &
                               '*(',pmassval,'/rho)**(1/',ndim,')'
 else
    print "(a)",'error: hfact = 0: can''t plot h vs rho exact solution'
    ierr = 1
 endif

 return
end subroutine exact_rhoh

end module rhoh
