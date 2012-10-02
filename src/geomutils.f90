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
 
 public :: change_coords, changecoords, changeveccoords
 
 private
 
contains

!-----------------------------------------------------------------
! transform all columns of data for a given particle
! to new coordinate system (does both coordinates
! and components of vectors)
! this version uses DOUBLE PRECISION for vals
!-----------------------------------------------------------------
subroutine change_coords(vals,ncols,ndim,icoords,icoordsnew,x0,v0)
 use params,        only:doub_prec
 use geometry,      only:coord_transform,vector_transform
 use labels,        only:ix,iamvec,ivx
 implicit none
 integer,                intent(in)    :: ncols,ndim,icoords,icoordsnew
 real(kind=doub_prec), dimension(ncols), intent(inout) :: vals
 real, dimension(ndim),  intent(in)    :: x0,v0
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
       if (icol.eq.ivx) then
          vec(1:ndim) = vals(iamvec(icol):iamvec(icol)+ndim-1) - v0(1:ndim)       
       else
          vec(1:ndim) = vals(iamvec(icol):iamvec(icol)+ndim-1)
       endif
       call vector_transform(xcoords,vec,ndim,icoords,vecnew,ndim,icoordsnew)
       vals(iamvec(icol):iamvec(icol)+ndim-1) = vecnew(1:ndim)
    endif
 enddo
 
end subroutine change_coords

!-------------------------------------------------------------------
! interface to coordinate-system transformations
!-------------------------------------------------------------------
subroutine changecoords(iplotx,iploty,xplot,yplot,ntot,ndim,dat)
 use geometry,      only:coord_transform,labelcoordsys
 use settings_data, only:xorigin,itrackpart,icoords,icoordsnew
 use labels,        only:is_coord,ix
 implicit none
 integer, intent(in) :: iplotx,iploty,ntot,ndim
 real, dimension(:), intent(inout) :: xplot,yplot
 real, dimension(:,:), intent(in)  :: dat
 real, dimension(ndim) :: xcoords,xcoordsnew
 integer :: j,ixcoord,iycoord
 logical :: iscoordx,iscoordy

 iscoordx = is_coord(iplotx,ndim)
 iscoordy = is_coord(iploty,ndim)
 if (iscoordx .or. iscoordy) then
    print*,'changing coords from ',trim(labelcoordsys(icoords)), &
           ' to ',trim(labelcoordsys(icoordsnew))
    if (itrackpart.gt.0) print*,' (relative to particle ',itrackpart,')'

    !--get offsets in range 1->ndim for the case where particle
    !  coords are not first in plot arrays
    ixcoord = iplotx - ix(1) + 1
    if (iscoordx .and. (ixcoord.le.0 .or. ixcoord.gt.ndim)) then
       print*,'ERROR in x coordinate offset in arrays: cannot change coordinate system'
       return
    endif
    iycoord = iploty - ix(1) + 1
    if (iscoordy .and. (iycoord.le.0 .or. iycoord.gt.ndim)) then
       print*,'ERROR in y coordinate offset in arrays: cannot change coordinate system'
       return
    endif

    do j=1,ntot
       if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
          xcoords(1:ndim) = dat(j,ix(1:ndim)) - dat(itrackpart,ix(1:ndim))
       else
          xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
       endif

       call coord_transform(xcoords(1:ndim),ndim,icoords, &
                            xcoordsnew(1:ndim),ndim,icoordsnew)
       if (iscoordx) xplot(j) = xcoordsnew(ixcoord)
       if (iscoordy) yplot(j) = xcoordsnew(iycoord)
    enddo
 endif

end subroutine changecoords

!-------------------------------------------------------------------
! interface to coordinate-system transformations for vectors
!-------------------------------------------------------------------
subroutine changeveccoords(iplot,xploti,ntot,ndim,dat)
 use geometry,      only:vector_transform,labelcoordsys
 use settings_data, only:xorigin,itrackpart,icoords,icoordsnew
 use labels,        only:ivx,iamvec,ix
 implicit none
 integer, intent(in) :: iplot,ntot,ndim
 real, dimension(:), intent(inout) :: xploti
 real, dimension(ndim) :: xcoords,vecnew,vecin
 real, dimension(:,:), intent(in)  :: dat
 integer :: j

 if (iamvec(iplot).gt.0) then
    if (iplot-iamvec(iplot)+1 .le. ndim) then
       print*,'changing vector component from ', &
        trim(labelcoordsys(icoords)),' to ',trim(labelcoordsys(icoordsnew))
       if (itrackpart.gt.0 .and. iamvec(iplot).eq.ivx) then
          print*,' (velocities relative to particle ',itrackpart,')'
       endif
       do j=1,ntot
          if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
             xcoords(1:ndim) = dat(j,ix(1:ndim)) - dat(itrackpart,ix(1:ndim))
             if (iamvec(iplot).eq.ivx) then
                vecin(1:ndim) = dat(j,iamvec(iplot):iamvec(iplot)+ndim-1) &
                                - dat(itrackpart,iamvec(iplot):iamvec(iplot)+ndim-1)
             else
                vecin(1:ndim) = dat(j,iamvec(iplot):iamvec(iplot)+ndim-1)
             endif
          else
             xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
             vecin(1:ndim) = dat(j,iamvec(iplot):iamvec(iplot)+ndim-1)
          endif
          call vector_transform(xcoords(1:ndim),vecin(1:ndim), &
               ndim,icoords,vecnew(1:ndim),ndim,icoordsnew)
          xploti(j) = vecnew(iplot-iamvec(iplot)+1)
       enddo
    else
       print*,'error: can''t convert vector components with ndimV > ndim'
    endif
 endif

 return
end subroutine changeveccoords

end module geomutils
