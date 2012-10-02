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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Utility module containing wrapper routines for coordinate
! transformations
!-----------------------------------------------------------------
module geomutils
 implicit none
 
 public :: change_coords
 
 private
 
contains

!-----------------------------------------------------------------
! transform all columns of data for a given particle
! to new coordinate system (does both coordinates
! and components of vectors)
!-----------------------------------------------------------------
subroutine change_coords(vals,ncols,ndim,icoords,icoordsnew,x0)
 use params,        only:doub_prec
 use geometry,      only:coord_transform,vector_transform
 use labels,        only:ix,iamvec
 implicit none
 integer,                intent(in)    :: ncols,ndim,icoords,icoordsnew
 real(kind=doub_prec), dimension(ncols), intent(inout) :: vals
 real, dimension(ndim),  intent(in)    :: x0
 real, dimension(ndim) :: xcoords,xcoordsnew,vec,vecnew
 integer :: iamvecprev,icol

 !--perform transformations on coordinates
 xcoords(1:ndim) = vals(ix(1:ndim)) - x0(1:ndim)
 call coord_transform(xcoords(1:ndim),ndim,icoords,xcoordsnew(1:ndim),ndim,icoordsnew)
 vals(ix(1:ndim)) = xcoordsnew(1:ndim)

 !--transform all vector quantities to new coord system
 iamvecprev = 0
 do icol=1,ncols - ndim + 1
    if (iamvec(icol).gt.0 .and. iamvec(icol).ne.iamvecprev) then                          
       iamvecprev = iamvec(icol)
       vec(1:ndim) = vals(iamvec(icol):iamvec(icol)+ndim-1)
       call vector_transform(xcoords,vec,ndim,icoords,vecnew,ndim,icoordsnew)
       vals(iamvec(icol):iamvec(icol)+ndim-1) = vecnew(1:ndim)
    endif
 enddo
 
end subroutine change_coords

end module geomutils
