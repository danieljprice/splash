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
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 2D interpolation
!
!----------------------------------------------------------------------

module interpolations2D
 use kernels, only:radkernel2,radkernel,cnormk2D,wfunc
 implicit none
 public :: interpolate2D, interpolate2D_xsec, interpolate2D_vec
 public :: interpolate_part, interpolate_part1

contains
!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_j w_j dat_j W(r-r_j, h_j)
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
!     Input: particle coordinates  : x,y    (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat    (npart)
!
!            number of pixels in x,y : npixx,npixy
!            pixel width             : pixwidth
!            option to normalise interpolation : normalise (.true. or .false.)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price 2003-2012
!--------------------------------------------------------------------------

subroutine interpolate2D(x,y,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,periodicx,periodicy)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise,periodicx,periodicy
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ipixi,jpixi
  real :: hi,hi1,radkern,q2,wab,const
  real :: term,termnorm,dx,dy,xpix,ypix

  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print "(1x,a)",'interpolating from particles to 2D grid (normalised)...'
  else
     print "(1x,a)",'interpolating from particles to 2D grid (non-normalised)...'
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print "(1x,a)",'interpolate2D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate2D: warning: ignoring some or all particles with h < 0'
  endif
  const = cnormk2D  ! normalisation constant
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
     radkern = radkernel*hi  ! radius of the smoothing kernel
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1      ! make sure they only contribute
        if (ipixmax.gt.npixx) ipixmax = npixx  ! to pixels in the image
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        jpixi = jpix
        if (periodicy) then
           if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
           if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1        
        endif
        ypix = ymin + (jpix-0.5)*pixwidthy
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           ipixi = ipix
           if (periodicx) then
              if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
              if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
           endif
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)         
           q2 = (dx*dx + dy*dy)*hi1*hi1
           !
           !--SPH kernel
           !
           wab = wfunc(q2)
           !
           !--calculate data value at this pixel using the summation interpolant
           !
           datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
           if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab

        enddo
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

end subroutine interpolate2D

!--------------------------------------------------------------------------
!
!     ** this version does vector quantities
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field    : vecsmoothx (npixx,npixy)
!                                      : vecsmoothy (npixx,npixy)
!
!     Daniel Price, University of Exeter, March 2005
!--------------------------------------------------------------------------

subroutine interpolate2D_vec(x,y,hh,weight,vecx,vecy,itype,npart, &
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,periodicx,periodicy)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx,vecsmoothy
  logical, intent(in) :: normalise,periodicx,periodicy
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ipixi,jpixi
  real :: hi,hi1,radkern,q2,wab,const
  real :: termnorm,termx,termy,dx,dy,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  if (normalise) then
     print "(1x,a)",'interpolating vector field from particles to 2D grid (normalised)...'
  else
     print "(1x,a)",'interpolating vector field from particles to 2D grid (non-normalised)...'
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate2D_vec: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate2D_vec: warning: ignoring some or all particles with h < 0'
  endif
  const = cnormk2D  ! normalisation constant
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
     radkern = radkernel*hi  ! radius of the smoothing kernel
     termx = termnorm*vecx(i)
     termy = termnorm*vecy(i)
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1
        if (ipixmax.gt.npixx) ipixmax = npixx  
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        jpixi = jpix
        if (periodicy) then
           if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
           if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1        
        endif
        ypix = ymin + (jpix-0.5)*pixwidthy
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           ipixi = ipix
           if (periodicx) then
              if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
              if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
           endif
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)
           q2 = (dx*dx + dy*dy)*hi1*hi1
           !
           !--SPH kernel
           !
           wab = wfunc(q2)
           !
           !--calculate data value at this pixel using the summation interpolant
           !
           vecsmoothx(ipixi,jpixi) = vecsmoothx(ipixi,jpixi) + termx*wab
           vecsmoothy(ipixi,jpixi) = vecsmoothy(ipixi,jpixi) + termy*wab
           if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab

        enddo
     enddo

  enddo over_parts
  !
  !--normalise dat arrays
  !
  if (normalise) then
     where (datnorm > 0.)
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif

  return

end subroutine interpolate2D_vec

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     this version takes any 1D cross section through a 2D data set
!     the 1D line is specified by two points, (x1,y1) and (x2,y2)
!     (ie. this is for arbitrary oblique cross sections)
!
!     NB: A similar version could be used for 2D oblique cross sections
!         of 3D data. In this case we would need to find the intersection
!         between the smoothing sphere and the cross section plane. However
!         in 3D it is simpler just to rotate the particles first and then take
!         a straight cross section.
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Feb 2004
!--------------------------------------------------------------------------

subroutine interpolate2D_xsec(x,y,hh,weight,dat,itype,npart,&
     x1,y1,x2,y2,datsmooth,npixx,normalise)

  implicit none
  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: x1,y1,x2,y2
  real, intent(out), dimension(npixx) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx) :: datnorm

  integer :: i,ipix,ipixmin,ipixmax
  real :: hi,hi1,radkern,q2,wab,const
  real :: term,termnorm,dx,dy,xpix,ypix,pixwidth,xpixwidth,xlength
  real :: gradient,yintercept,aa,bb,cc,determinant,det
  real :: xstart,xend,ystart,yend,rstart,rend
  real :: tol
  logical :: xsame, ysame, debug

  debug = .false.
  !
  !--check for errors in input
  !
  tol = 1.e-3
  ysame = (abs(y2 - y1).lt.tol)
  xsame = (abs(x2 - x1).lt.tol)
  if (xsame.and.ysame) then
     print*,'error: interpolate2D_xsec: zero length cross section'
     return
  endif
  if (npixx.eq.0) then
     print*,'error: interpolate2D_xsec: npix = 0 '
     return
  endif
  print*,'oblique 1D cross section through 2D data: npix =',npixx
  !
  !--work out the equation of the line y = mx + c from the two points input
  !
  gradient = 0.
  if (.not.xsame) gradient = (y2-y1)/(x2-x1)
  yintercept = y2 - gradient*x2
  print*,'cross section line: y = ',gradient,'x + ',yintercept
  !
  !--work out length of line and divide into pixels
  !
  xlength = sqrt((x2-x1)**2 + (y2-y1)**2)
  pixwidth = xlength/real(npixx)
  xpixwidth = (x2 - x1)/real(npixx)
  if (debug) then
     print*,'length of line = ',xlength
     print*,'pixel width = ',pixwidth, ' in x direction = ',xpixwidth
  endif
  !
  !--now interpolate to the line of pixels
  !
  datsmooth = 0.
  datnorm = 0.
  const = cnormk2D   ! normalisation constant
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
     radkern = radkernel*hi    ! radius of the smoothing kernel
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !  to do this we need to work out the two points at which the line
     !  intersects the particles smoothing circle
     !  given by the equation (x-xi)^2 + (y-yi)^2 = (2h)^2.
     !  The x co-ordinates of these points are the solutions to a
     !  quadratic with co-efficients:

     aa = 1. + gradient**2
     bb = 2.*gradient*(yintercept - y(i)) - 2.*x(i)
     cc = x(i)**2 + y(i)**2 - 2.*yintercept*y(i) + yintercept**2 &
          - radkern**2
     !
     !--work out whether there are any real solutions and find them
     !
     determinant = bb**2 - 4.*aa*cc
     if (determinant < 0) then
        !!print*,' particle ',i,': does not contribute ',x(i),y(i)
     else
        det = sqrt(determinant)
        xstart = (-bb - det)/(2.*aa)
        xend =  (-bb + det)/(2.*aa)
        if (xstart.lt.x1) xstart = x1
        if (xstart.gt.x2) xstart = x2
        if (xend.gt.x2) xend = x2
        if (xend.lt.x1) xend = x1
        ystart = gradient*xstart + yintercept
        yend = gradient*xend + yintercept
        !
        !--work out position in terms of distance (no. of pixels) along the line
        !
        rstart = sqrt((xstart-x1)**2 + (ystart-y1)**2)
        rend = sqrt((xend-x1)**2 + (yend-y1)**2)

        ipixmin = int(rstart/pixwidth)
        ipixmax = int(rend/pixwidth) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (ipixmax.lt.1) ipixmax = 1
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (ipixmin.gt.npixx) ipixmax = npixx
               !
        !--loop over pixels, adding the contribution from this particle
        !
        !if (debug) print*,' particle ',i,': ',ipixmin,ipixmax,xstart,x(i),xend
        do ipix = ipixmin,ipixmax

           xpix = x1 + (ipix-0.5)*xpixwidth
           ypix = gradient*xpix + yintercept
           dy = ypix - y(i)
           dx = xpix - x(i)
           q2 = (dx*dx + dy*dy)*hi1*hi1
           !
           !--SPH kernel
           !
           wab = wfunc(q2)
           !
           !--calculate data value at this pixel using the summation interpolant
           !
           datsmooth(ipix) = datsmooth(ipix) + term*wab
           if (normalise) datnorm(ipix) = datnorm(ipix) + termnorm*wab

        enddo

     endif

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

end subroutine interpolate2D_xsec

!--------------------------------------------------------------------------
!     subroutine to render particles onto a pixel array
!     at the maximum or minimum colour
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------
subroutine interpolate_part(x,y,hh,npart,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)
  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh
  real, intent(in) :: xmin,ymin,pixwidth,datval
  real, intent(inout), dimension(npixx,npixy) :: datsmooth
  real, intent(inout), dimension(npixx,npixy), optional :: brightness
  integer :: i

  if (pixwidth.le.0.) then
     print "(1x,a)",'interpolate_part: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate_part: warning: ignoring some or all particles with h < 0'
  endif
  !
  !--loop over particles
  !
  if (present(brightness)) then
     do i=1,npart
        call interpolate_part1(x(i),y(i),hh(i),xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval)
     enddo
  else
     do i=1,npart
        call interpolate_part1(x(i),y(i),hh(i),xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)
     enddo
  endif
  return

end subroutine interpolate_part

!--------------------------------------------------------------------------
!     subroutine to render a single particle onto a pixel array
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------
subroutine interpolate_part1(xi,yi,hi,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)
  implicit none
  real, intent(in) :: xi,yi,hi,xmin,ymin,pixwidth,datval
  integer, intent(in) :: npixx,npixy
  real, intent(inout), dimension(npixx,npixy) :: datsmooth
  real, intent(inout), dimension(npixx,npixy), optional :: brightness
  integer :: ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: radkern,radkern2,rab2
  real :: dx,dy2,xpix,ypix

  !
  !--skip particles with wrong h's
  !
  if (hi.le.tiny(hi)) return
  !
  !--set kernel related quantities
  !
  radkern = max(hi,2.*pixwidth)
  radkern2 = radkern*radkern  ! radius of the smoothing kernel
  !
  !--for each particle work out which pixels it contributes to
  !
  ipixmin = int((xi - radkern - xmin)/pixwidth)
  jpixmin = int((yi - radkern - ymin)/pixwidth)
  ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
  jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

  if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
  if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
  if (ipixmax.gt.npixx) ipixmax = npixx
  if (jpixmax.gt.npixy) jpixmax = npixy
  !
  !--loop over pixels, adding the contribution from this particle
  !
  do jpix = jpixmin,jpixmax
     ypix = ymin + (jpix-0.5)*pixwidth
     dy2 = (ypix - yi)**2
     do ipix = ipixmin,ipixmax
        xpix = xmin + (ipix-0.5)*pixwidth
        dx = xpix - xi
        rab2 = dx**2 + dy2
        !
        !--set data value at this pixel to maximum
        !
        if (rab2.lt.radkern2) then
           datsmooth(ipix,jpix) = datval
           if (present(brightness)) then
              brightness(ipix,jpix) = 1.0
           endif
        endif
     enddo
  enddo

  return
end subroutine interpolate_part1

!--------------------------------------------------------------------------
!     subroutine to interpolate from arbitrary data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_j dat_j W(r-r_j, \Delta) / sum_j W(r-r_j, \Delta)
!
!     where _j is the quantity at the neighbouring particle j and
!     W is the smoothing kernel. THe interpolation is normalised.
!
!     Input: data points  : x,y    (npart)
!            third scalar to use as weight : dat (npart) (optional)
!
!            number of pixels in x,y : npixx,npixy
!            pixel width             : pixwidth
!            option to normalise interpolation : normalise (.true. or .false.)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price 2016
!--------------------------------------------------------------------------

subroutine interpolate2D_pixels(x,y,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,periodicx,periodicy,dat)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise,periodicx,periodicy
  real, intent(in), dimension(npart), optional :: dat
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ipixi,jpixi
  real :: hi,hi1,radkernx,radkerny,q2,wab,const
  real :: term,termnorm,dx,dy,xpix,ypix

  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print "(1x,a)",'interpolating from particles to 2D pixels (normalised)...'
  else
     print "(1x,a)",'interpolating from particles to 2D pixels (non-normalised)...'
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print "(1x,a)",'interpolate2D: error: pixel width <= 0'
     return
  endif
  const = cnormk2D  ! normalisation constant
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
     termnorm = const
     !if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = 2.*pixwidthx
     if (hi.le.tiny(hi)) cycle over_parts
     hi1 = 1./hi
     !
     !--set kernel related quantities
     !
     radkernx = radkernel*pixwidthx  ! radius of the smoothing kernel
     radkerny = radkernel*pixwidthy  ! radius of the smoothing kernel
     if (present(dat)) then
        term = termnorm*dat(i)
     else
        term = termnorm
     endif
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkernx - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkerny - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkernx - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkerny - ymin)/pixwidthy) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1      ! make sure they only contribute
        if (ipixmax.gt.npixx) ipixmax = npixx  ! to pixels in the image
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif
     !
     !--loop over pixels, adding the contribution from this particle
     !
     if (i.eq.1) print*,jpixmin,jpixmax,ipixmin,ipixmax
     do jpix = jpixmin,jpixmax
        jpixi = jpix
        if (periodicy) then
           if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
           if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
        endif
        ypix = ymin + (jpix-0.5)*pixwidthy
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           ipixi = ipix
           if (periodicx) then
              if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
              if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
           endif
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)         
           q2 = (dx*dx/pixwidthx**2 + dy*dy/pixwidthy**2)
           !
           !--SPH kernel
           !
           wab = wfunc(q2)
           !
           !--calculate data value at this pixel using the summation interpolant
           !
           datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
           if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab

        enddo
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

end subroutine interpolate2D_pixels

end module interpolations2D
