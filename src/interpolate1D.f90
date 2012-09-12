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

!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 1D interpolation
!
!----------------------------------------------------------------------

module interpolations1D
 implicit none
 public :: interpolate1D

contains

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even 1D grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_j w_j W(r-r_j, h_j)
!
!     where _j is the quantity at the neighbouring particle j and
!     W is the smoothing kernel, for which we use the usual cubic spline.
!     For an SPH interpolation the weight for each particle should be
!     the dimensionless quantity
!
!     w_j = m_j / (rho_j * h_j**ndim)
!
!     Other weights can be used (e.g. constants), but in this case the
!     normalisation option should also be set.
!
!     Input: particle coordinates  : x      (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat    (npart)
!
!            number of pixels in x   : npixx
!            pixel width             : pixwidth
!            option to normalise interpolation : normalise (.true. or .false.)
!
!     Output: smoothed data            : datsmooth (npixx)
!
!     Written by Daniel Price 2003-2006
!--------------------------------------------------------------------------

subroutine interpolate1D(x,hh,weight,dat,itype,npart,  &
     xmin,datsmooth,npixx,pixwidth,normalise)

  use kernels, only:cnormk1D,radkernel,wfunc
  implicit none
  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,pixwidth
  real, intent(out), dimension(npixx) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx) :: datnorm

  integer :: i,ipix,ipixmin,ipixmax
  real :: hi,hi1,radkern,q2,wab,const
  real :: term,termnorm,dx,xpix

  datsmooth = 0.
  term = 0.
  if (normalise) then
     print*,'interpolating (normalised) from particles to 1D grid: npix,xmin,max=',npixx,xmin,xmin+npixx*pixwidth
  else
     print*,'interpolating (non-normalised) from particles to 1D grid: npix,xmin,max=',npixx,xmin,xmin+npixx*pixwidth
  endif
  if (pixwidth.le.0.) then
     print*,'interpolate1D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate1D: warning: ignoring some or all particles with h < 0'
  endif
  const = cnormk1D  ! normalisation constant
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = radkernel*hi   ! radius of the smoothing kernel
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth) + 1

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (ipixmax.gt.npixx) ipixmax = npixx ! to pixels in the image
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do ipix = ipixmin,ipixmax
        xpix = xmin + (ipix-0.5)*pixwidth
        dx = xpix - x(i)
        q2 = dx*dx*hi1*hi1
        !
        !--SPH kernel - standard cubic spline
        !
        wab = wfunc(q2)
        !
        !--calculate data value at this pixel using the summation interpolant
        !
        datsmooth(ipix) = datsmooth(ipix) + term*wab
        if (normalise) datnorm(ipix) = datnorm(ipix) + termnorm*wab
     enddo

  enddo over_parts
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > 0.)
        datsmooth = datsmooth/datnorm
     end where
  endif

  return

end subroutine interpolate1D

end module interpolations1D
