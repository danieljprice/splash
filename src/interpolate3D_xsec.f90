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
!  Module containing all of the routines required for cross sections
!  through 3D data
!
!----------------------------------------------------------------------

module xsections3D
 use kernels, only:cnormk3D,radkernel,radkernel2,wfunc
 implicit none
 public :: interpolate3D_fastxsec, interpolate3D_xsec_vec

contains

!--------------------------------------------------------------------------
!
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!            cross section location: zslice
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------

subroutine interpolate3D_fastxsec(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zslice
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,q2,wab,const,xi,yi,hi21
  real :: termnorm,term,dy,dy2,dz,dz2,ypix,rescalefac
  real, dimension(npixx) :: dx2i

  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised)...',zslice
  else
     print*,'taking fast cross section (non-normalised)...',zslice
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate3D_xsec: error: pixel width <= 0'
     return
  elseif (npart.le.0) then
     print*,'interpolate3D_xsec: error: npart = 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec: WARNING: ignoring some or all particles with h < 0'
  endif
  const = cnormk3D
  !
  !--renormalise dat array by first element to speed things up
  !
  if (dat(1).gt.tiny(dat)) then
     rescalefac = dat(1)
  else
     rescalefac = 1.0
  endif
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
     hi1 = 1./hi
     hi21 = hi1*hi1
     radkern = radkernel*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2*hi21
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (dz2 .lt. radkernel2) then

        xi = x(i)
        yi = y(i)
        termnorm = const*weight(i)
        term = termnorm*dat(i)/rescalefac
        !
        !--for each particle work out which pixels it contributes to
        !
        ipixmin = int((xi - radkern - xmin)/pixwidthx)
        jpixmin = int((yi - radkern - ymin)/pixwidthy)
        ipixmax = int((xi + radkern - xmin)/pixwidthx) + 1
        jpixmax = int((yi + radkern - ymin)/pixwidthy) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xmin + (ipix-0.5)*pixwidthx - xi)**2)*hi21 + dz2
        enddo
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixmin,ipixmax
              q2 = dx2i(ipix) + dy2
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.radkernel2) then
                 wab = wfunc(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                 if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo over_parts                    ! over particles
  !
  !--normalise dat array
  !
  if (normalise) then
     !--normalise everywhere (required if not using SPH weighting)
     where (datnorm > tiny(datnorm))
        datsmooth = datsmooth/datnorm
     end where
  endif
  datsmooth = datsmooth*rescalefac

  return

end subroutine interpolate3D_fastxsec

!--------------------------------------------------------------------------
!     program to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
!
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!            cross section location: zslice
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------

subroutine interpolate3D_xsec_vec(x,y,z,hh,weight,vecx,vecy,itype,npart,&
     xmin,ymin,zslice,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zslice
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,q2,wab,const
  real :: termx,termy,termnorm,dx,dy,dz,dz2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised)...',zslice
  else
     print*,'taking fast cross section (non-normalised)...',zslice
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate3D_xsec_vec: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec_vec: WARNING: ignoring some or all particles with h < 0'
  endif
  const = cnormk3D ! normalisation constant (3D)
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
     hi1 = 1./hi
     radkern = radkernel*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (abs(dz) .lt. radkern) then
        termnorm = const*weight(i)
        termx = termnorm*vecx(i)
        termy = termnorm*vecy(i)
        !
        !--for each particle work out which pixels it contributes to
        !
        ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
        jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
        ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
        jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - y(i)
           do ipix = ipixmin,ipixmax
              xpix = xmin + (ipix-0.5)*pixwidthx
              dx = xpix - x(i)
              q2 = (dx*dx + dy*dy + dz2)*hi1*hi1
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.radkernel2) then
                 wab = wfunc(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !
                 vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab
                 vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab
                 if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo over_parts                    ! over particles
  !
  !--normalise dat array(s)
  !
  if (normalise) then
     where (datnorm > tiny(datnorm))
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif

  return

end subroutine interpolate3D_xsec_vec

end module xsections3D
