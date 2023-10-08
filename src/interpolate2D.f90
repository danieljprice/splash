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
!  Copyright (C) 2005-2019 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 2D interpolation
!
!----------------------------------------------------------------------

module interpolations2D
 use kernels,       only:radkernel2,radkernel,cnormk2D,wfunc,pint,select_kernel,soft_func
 use timing,        only:wall_time,print_time
 use interpolation, only:iroll
 implicit none
 public :: interpolate2D, interpolate2D_xsec, interpolate2D_vec
 public :: interpolate_part, interpolate_part1
 public :: interpolate2D_pixels, interpolate2D_fromgrid

 private

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
!     Written by Daniel Price 2003-2020
!     Exact rendering implemented by Maya Petkova and Daniel Price 2018
!--------------------------------------------------------------------------

subroutine interpolate2D(x,y,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,exact,periodicx,periodicy,iverbose)

 integer, intent(in) :: npart,npixx,npixy
 real, intent(in), dimension(npart) :: x,y,hh,weight,dat
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
 real, intent(out), dimension(npixx,npixy) :: datsmooth
 logical, intent(in) :: normalise,exact,periodicx,periodicy
 integer, intent(in) :: iverbose
 real, dimension(npixx,npixy) :: datnorm

 integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
 integer :: ipixi,jpixi
 real :: hi,hi1,radkern,q2,wab,const
 real :: term,termnorm,dx,dy,dy2,xpix,ypix,hi21
 real :: t_start,t_end,t_used,xmini,ymini,denom,datnorm_min
 real, dimension(npixx) :: dx2i

! Maya
 real :: pixint,d1,d2,r0

 call wall_time(t_start)
 if (.not.associated(wfunc)) call select_kernel(0)

 datsmooth = 0.
 datnorm = 0.
 if (iverbose > 0) then
    if (exact) then
       print "(1x,a)",'interpolating from particles to 2D grid (exact)...'
    elseif (normalise) then
       print "(1x,a)",'interpolating from particles to 2D grid (normalised)...'
    else
       print "(1x,a)",'interpolating from particles to 2D grid (non-normalised)...'
    endif
 endif
 if (abs(pixwidthx) < tiny(0.) .or. abs(pixwidthy) < tiny(0.)) then
    print "(1x,a)",'interpolate2D: error: pixel size = 0'
    return
 endif
 ! allow for negative pixel widths (just flips image)
 xmini = min(xmin,xmin + npixy*pixwidthy)
 ymini = min(ymin,ymin + npixy*pixwidthx)

 if (any(hh(1:npart) <= tiny(hh)) .and. iverbose >= 0) then
    print*,'interpolate2D: warning: ignoring some or all particles with h < 0'
 endif
 const = cnormk2D  ! normalisation constant
 !
 !--loop over particles
 !
!$omp parallel do default(none) &
!$omp shared(hh,x,y,weight,dat,datsmooth,datnorm,itype,npart) &
!$omp shared(xmini,ymini,normalise,exact,radkernel,const) &
!$omp shared(pixwidthx,pixwidthy,periodicx,periodicy,npixx,npixy) &
!$omp private(hi,radkern,wab) &
!$omp private(hi1,hi21,term,termnorm,jpixi,ipixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,denom) &
!$omp private(q2,xpix,ypix,dx,dy,dx2i,dy2,r0,d1,d2,pixint) &
!$omp private(i,ipix,jpix) &
!$omp schedule (guided, 10)
 over_parts: do i=1,npart
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0) cycle over_parts
    !
    !--skip particles with zero weights
    !
    termnorm = const*weight(i)
    if (termnorm <= 0.) cycle over_parts
    !
    !--skip particles with wrong h's
    !
    hi = hh(i)
    if (hi <= tiny(hi)) cycle over_parts
    !
    !--set kernel related quantities
    !
    hi1 = 1./hi
    hi21 = hi1*hi1
    radkern = radkernel*hi  ! radius of the smoothing kernel
    term = termnorm*dat(i)
    denom = 1./abs(pixwidthx*pixwidthy)/const*hi**2
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((x(i) - radkern - xmini)/abs(pixwidthx))
    jpixmin = int((y(i) - radkern - ymini)/abs(pixwidthy))
    ipixmax = int((x(i) + radkern - xmini)/abs(pixwidthx)) + 1
    jpixmax = int((y(i) + radkern - ymini)/abs(pixwidthy)) + 1
    ipixmin = min(ipixmin,ipixmax)
    ipixmax = max(ipixmin,ipixmax)
    jpixmin = min(jpixmin,jpixmax)
    jpixmax = max(jpixmin,jpixmax)

    if (.not.periodicx) then
       if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
       if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
    endif
    if (.not.periodicy) then
       if (jpixmin < 1)     jpixmin = 1
       if (jpixmax > npixy) jpixmax = npixy
    endif

    if (exact .and. ipixmax-ipixmin < 5) then
       !
       !--loop over pixels boundaries, adding the contribution from this particle
       !
       !--first pixel row
       !
       if (jpixmax >= jpixmin) then
          jpix = jpixmin
          jpixi = jpix
          if (periodicy) jpixi = iroll(jpix,npixy)
          ypix = ymini + (jpix-0.5)*pixwidthy
          dy = ypix - y(i)

          do ipix = ipixmin,ipixmax
             ipixi = ipix
             if (periodicx) ipixi = iroll(ipix,npixx)
             xpix = xmini + (ipix-0.5)*pixwidthx
             dx = xpix - x(i)

             !--top boundary
             r0 = 0.5*pixwidthy - dy
             d1 = 0.5*pixwidthx + dx
             d2 = 0.5*pixwidthx - dx
             pixint = pint(r0, d1, d2, hi1)

             wab = pixint*denom
             !$omp atomic
             datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
             if (normalise) then
                !$omp atomic
                datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
             endif
          enddo
       endif

       !
       !--first pixel column
       !
       if (ipixmax >= ipixmin) then
          ipix = ipixmin
          ipixi = ipix
          if (periodicx) ipixi = iroll(ipix,npixx)
          xpix = xmini + (ipix-0.5)*pixwidthx
          dx = xpix - x(i)

          do jpix = jpixmin,jpixmax
             jpixi = jpix
             if (periodicy) jpixi = iroll(jpixi,npixy)
             ypix = ymini + (jpix-0.5)*pixwidthy
             dy = ypix - y(i)

             !--left boundary
             r0 = 0.5*pixwidthx - dx
             d1 = 0.5*pixwidthy - dy
             d2 = 0.5*pixwidthy + dy
             pixint = pint(r0, d1, d2, hi1)

             wab = pixint*denom
             !$omp atomic
             datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
             if (normalise) then
                !$omp atomic
                datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
             endif
          enddo
       endif

       !
       !--other pixels
       !
       do jpix = jpixmin,jpixmax
          jpixi = jpix
          if (periodicy) jpixi = iroll(jpix,npixy)
          ypix = ymini + (jpix-0.5)*pixwidthy
          dy = ypix - y(i)

          do ipix = ipixmin,ipixmax
             ipixi = ipix
             if (periodicx) ipixi = iroll(ipix,npixx)
             xpix = xmini + (ipix-0.5)*pixwidthx
             dx = xpix - x(i)
             !
             !--Kernel integral
             !
             !--bottom boundary
             r0 = 0.5*pixwidthy + dy
             d1 = 0.5*pixwidthx - dx
             d2 = 0.5*pixwidthx + dx
             pixint = pint(r0, d1, d2, hi1)

             wab = pixint*denom
             !$omp atomic
             datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
             if (normalise) then
                !$omp atomic
                datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
             endif

             if (jpix < jpixmax) then
                jpixi = jpix+1
                if (periodicy) jpixi = iroll(jpixi,npixy)

                !$omp atomic
                datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) - term*wab
                if (normalise) then
                   !$omp atomic
                   datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) - termnorm*wab
                endif

                jpixi = jpix
                if (periodicy) jpixi = iroll(jpix,npixy)
             endif

             !--right boundary
             r0 = 0.5*pixwidthx + dx
             d1 = 0.5*pixwidthy + dy
             d2 = 0.5*pixwidthy - dy
             pixint = pint(r0, d1, d2, hi1)

             wab = pixint*denom
             !$omp atomic
             datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
             if (normalise) then
                !$omp atomic
                datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
             endif

             if (ipix < ipixmax) then
                ipixi = ipix+1
                if (periodicx) ipixi = iroll(ipixi,npixx)

                !$omp atomic
                datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) - term*wab
                if (normalise) then
                   !$omp atomic
                   datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) - termnorm*wab
                endif
             endif
          enddo
       enddo
    else
       !
       !--precalculate an array of dx2 for this particle (optimisation)
       !
       do ipix=ipixmin,ipixmax
          dx2i(ipix) = ((xmini + (ipix-0.5)*pixwidthx - x(i))**2)*hi21
       enddo
       !
       !--loop over pixels, adding the contribution from this particle
       !
       do jpix = jpixmin,jpixmax
          jpixi = jpix
          if (periodicy) jpixi = iroll(jpix,npixy)
          ypix = ymini + (jpix-0.5)*pixwidthy
          dy = ypix - y(i)
          dy2 = dy*dy*hi21
          do ipix = ipixmin,ipixmax
             ipixi = ipix
             if (periodicx) ipixi = iroll(ipix,npixx)
!             xpix = xmin + (ipix-0.5)*pixwidthx
!             dx = xpix - x(i)
             q2 = dx2i(ipix) + dy2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
             !q2 = (dx*dx + dy*dy)*hi1*hi1


             !
             !--SPH kernel
             !
             wab = wkernel(q2)
             !
             !--calculate data value at this pixel using the summation interpolant
             !
             !$omp atomic
             datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
             if (normalise) then
                !$omp atomic
                datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
             endif
          enddo
       enddo
    endif

 enddo over_parts
 !$omp end parallel do

 !if (exact) then
 !    print*, 'sum of datpix = ', sum(datsmooth)/(npixx*npixy)
 !    print*, 'max of datpix = ', maxval(datsmooth)
 !    print*, 'min of datpix = ', minval(datsmooth)
 !endif
 !
 !--normalise dat array
 !
 if (normalise) then
    !where (datnorm > 0.)
    !   datsmooth = datsmooth/datnorm
    !end where
    !
    ! compute the minimum possible value for datnorm, then multiply
    ! with a kernel-softened version of 1/x to avoid dividing by zero
    !
    datnorm_min = minval(const*weight,mask=(weight > 0.))
    datsmooth = datsmooth*soft_func(datnorm,datnorm_min)
 endif

 call wall_time(t_end)
 t_used = t_end - t_start
 if (t_used > 10.) call print_time(t_used)

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
!     Exact rendering implemented by Maya Petkova and Daniel Price 2018å
!--------------------------------------------------------------------------

subroutine interpolate2D_vec(x,y,hh,weight,vecx,vecy,itype,npart, &
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,exact,periodicx,periodicy)

 integer, intent(in) :: npart,npixx,npixy
 real, intent(in), dimension(npart) :: x,y,hh,weight,vecx,vecy
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
 real, intent(out), dimension(npixx,npixy) :: vecsmoothx,vecsmoothy
 logical, intent(in) :: normalise,exact,periodicx,periodicy
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
 if (pixwidthx <= 0. .or. pixwidthy <= 0.) then
    print*,'interpolate2D_vec: error: pixel width <= 0'
    return
 endif
 if (any(hh(1:npart) <= tiny(hh))) then
    print*,'interpolate2D_vec: warning: ignoring some or all particles with h < 0'
 endif
 const = cnormk2D  ! normalisation constant
 if (.not.associated(wfunc)) call select_kernel(0)
 !
 !--loop over particles
 !
 over_parts: do i=1,npart
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0) cycle over_parts
    !
    !--skip particles with zero weights
    !
    termnorm = const*weight(i)
    if (termnorm <= 0.) cycle over_parts
    !
    !--skip particles with wrong h's
    !
    hi = hh(i)
    if (hi <= tiny(hi)) cycle over_parts
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
       if (ipixmin < 1)     ipixmin = 1
       if (ipixmax > npixx) ipixmax = npixx
    endif
    if (.not.periodicy) then
       if (jpixmin < 1)     jpixmin = 1
       if (jpixmax > npixy) jpixmax = npixy
    endif
    !
    !--loop over pixels, adding the contribution from this particle
    !
    do jpix = jpixmin,jpixmax
       jpixi = jpix
       if (periodicy) jpixi = iroll(jpix,npixy)
       ypix = ymin + (jpix-0.5)*pixwidthy
       dy = ypix - y(i)
       do ipix = ipixmin,ipixmax
          ipixi = ipix
          if (periodicx) ipixi = iroll(ipix,npixx)
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
 ysame = (abs(y2 - y1) < tol)
 xsame = (abs(x2 - x1) < tol)
 if (xsame.and.ysame) then
    print*,'error: interpolate2D_xsec: zero length cross section'
    return
 endif
 if (npixx==0) then
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
 if (.not.associated(wfunc)) call select_kernel(0)
 !
 !--loop over particles
 !
 over_parts: do i=1,npart
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0) cycle over_parts
    !
    !--skip particles with zero weights
    !
    termnorm = const*weight(i)
    if (termnorm <= 0.) cycle over_parts
    !
    !--skip particles with wrong h's
    !
    hi = hh(i)
    if (hi <= tiny(hi)) cycle over_parts
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
       if (xstart < x1) xstart = x1
       if (xstart > x2) xstart = x2
       if (xend > x2) xend = x2
       if (xend < x1) xend = x1
       ystart = gradient*xstart + yintercept
       yend = gradient*xend + yintercept
       !
       !--work out position in terms of distance (no. of pixels) along the line
       !
       rstart = sqrt((xstart-x1)**2 + (ystart-y1)**2)
       rend = sqrt((xend-x1)**2 + (yend-y1)**2)

       ipixmin = int(rstart/pixwidth)
       ipixmax = int(rend/pixwidth) + 1

       if (ipixmin < 1) ipixmin = 1 ! make sure they only contribute
       if (ipixmax < 1) ipixmax = 1
       if (ipixmax > npixx) ipixmax = npixx
       if (ipixmin > npixx) ipixmax = npixx
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

end subroutine interpolate2D_xsec

!--------------------------------------------------------------------------
!     subroutine to render particles onto a pixel array
!     at the maximum or minimum colour
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------
subroutine interpolate_part(x,y,hh,npart,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)
 integer, intent(in) :: npart,npixx,npixy
 real, intent(in), dimension(npart) :: x,y,hh
 real, intent(in) :: xmin,ymin,pixwidth,datval
 real, intent(inout), dimension(npixx,npixy) :: datsmooth
 real, intent(inout), dimension(npixx,npixy), optional :: brightness
 integer :: i

 if (pixwidth <= 0.) then
    print "(1x,a)",'interpolate_part: error: pixel width <= 0'
    return
 endif
 if (any(hh(1:npart) <= tiny(hh))) then
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

end subroutine interpolate_part

!--------------------------------------------------------------------------
!     subroutine to render a single particle onto a pixel array
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------
subroutine interpolate_part1(xi,yi,hi,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)
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
 if (hi <= tiny(hi)) return
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

 if (ipixmin < 1) ipixmin = 1 ! make sure they only contribute
 if (jpixmin < 1) jpixmin = 1 ! to pixels in the image
 if (ipixmax > npixx) ipixmax = npixx
 if (jpixmax > npixy) jpixmax = npixy
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
       if (rab2 < radkern2) then
          datsmooth(ipix,jpix) = datval
          if (present(brightness)) then
             brightness(ipix,jpix) = 1.0
          endif
       endif
    enddo
 enddo

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
!     W is the smoothing kernel. The interpolation is normalised.
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
!     Written by Daniel Price 2017
!--------------------------------------------------------------------------

subroutine interpolate2D_pixels(x,y,itype,npart, &
     xmin,ymin,xmax,ymax,datsmooth,npixx,npixy,&
     normalise,adaptive,dat,datpix2,fac,weights)

 use timing, only:wall_time,print_time
 integer, intent(in) :: npart,npixx,npixy
 real, intent(in), dimension(npart) :: x,y
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,xmax,ymax
 real, intent(out), dimension(npixx,npixy) :: datsmooth
 logical, intent(in) :: normalise,adaptive
 real, intent(in), dimension(npart), optional :: dat,weights
 real, dimension(npixx,npixy), intent(out), optional :: datpix2
 real, intent(in), optional :: fac

 real, dimension(npixx,npixy) :: datnorm,datold
 real, dimension(npixx) :: dx2i,qq2,wabi

 integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,its,itsmax
 real :: hi,hi1,radkernx,radkerny,q2,wab,const
 real :: term,termnorm,dy,xpix,ypix,ddx,ddy
 real :: xi,yi,pixwidthx,pixwidthy,dy2
 real :: t1,t2,hfac

 if (adaptive) then
    print "(1x,a)",'interpolating from particles to 2D pixels (adaptive)...'
 else
    print "(1x,a)",'interpolating from particles to 2D pixels...'
 endif

 if (adaptive) then
    itsmax = 3
 else
    itsmax = 1
 endif
 ! default smoothing length in units of pixel scale
 hfac = 1.5
 if (present(fac)) then
    if (fac > 0.) hfac = fac
 endif
 call wall_time(t1)
 if (.not.associated(wfunc)) call select_kernel(0)

 iterations: do its=1,itsmax
    datsmooth = 0.
    datnorm = 0.

    const = cnormk2D  ! normalisation constant
    !
    !--loop over particles
    !
    ddx = npixx/(xmax - xmin)
    ddy = npixy/(ymax - ymin)

    pixwidthx = 1. !/npixx
    pixwidthy = 1. !/npixy
    if (pixwidthx <= 0. .or. pixwidthy <= 0.) then
       print "(1x,a)",'interpolate2D: error: pixel width <= 0'
       return
    endif

    !$omp parallel do default(none) &
    !$omp shared(npart,itype,x,y,xmin,ymin,ddx,ddy,its,itsmax,weights) &
    !$omp shared(datold,datsmooth,datnorm,npixx,npixy,const,radkernel,radkernel2,dat) &
    !$omp shared(pixwidthx,pixwidthy,normalise,hfac) &
    !$omp private(i,xi,yi,ipix,jpix,hi,hi1) &
    !$omp private(radkernx,radkerny,ipixmin,ipixmax,jpixmin,jpixmax) &
    !$omp private(dx2i,xpix,ypix,dy,dy2,q2,wab,term,termnorm,qq2,wabi)
    over_parts: do i=1,npart
       !
       !--skip particles with itype < 0
       !
       if (itype(i) < 0) cycle over_parts

       !
       !--scale particle positions into viewport coordinates
       !
       xi = (x(i) - xmin)*ddx
       yi = (y(i) - ymin)*ddy
       hi = 1.0*pixwidthx
       if (itsmax==1) hi=hi*hfac  ! in units of pixel spacing

       ipix = int(xi)
       jpix = int(yi)
       if (its > 1 .and. ipix >= 1 .and. ipix <= npixx.and. jpix >= 1 .and. jpix <= npixy) then
          hi = max(min(hfac/sqrt(datold(ipix,jpix)),100.),1.)
       endif
       hi1 = 1./hi
       termnorm = const*hi1*hi1
       !
       !--set kernel related quantities
       !
       radkernx = radkernel*hi  ! radius of the smoothing kernel
       radkerny = radkernel*hi  ! radius of the smoothing kernel
       if (present(weights)) then
          termnorm = termnorm*weights(i)
       endif
       if (present(dat)) then
          term = termnorm*dat(i)
       else
          term = termnorm
       endif
       if (termnorm <= 0.) cycle over_parts
       !
       !--for each particle work out which pixels it contributes to
       !
       ipixmin = int((xi - radkernx))
       jpixmin = int((yi - radkerny))
       ipixmax = int((xi + radkernx)) + 1
       jpixmax = int((yi + radkerny)) + 1

       if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
       if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
       if (jpixmin < 1)     jpixmin = 1
       if (jpixmax > npixy) jpixmax = npixy
       !
       !--precalculate an array of dx2 for this particle (optimisation)
       !
       do ipix=ipixmin,ipixmax
          xpix = (ipix-0.5)*pixwidthx
          dx2i(ipix) = ((xpix - xi)**2)*hi1*hi1
       enddo
       !
       !--loop over pixels, adding the contribution from this particle
       !
       do jpix = jpixmin,jpixmax
          ypix = (jpix-0.5)*pixwidthy
          dy = ypix - yi
          dy2 = dy*dy*hi1*hi1

          do ipix = ipixmin,ipixmax
             q2 = dx2i(ipix) + dy2
             !
             !--SPH kernel
             !
             if (q2 < radkernel2) then
                wab = wkernel(q2)
                !
                !--calculate data value at this pixel using the summation interpolant
                !
                !$omp atomic
                datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                if (normalise) then
                   !$omp atomic
                   datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
                endif
             endif

          enddo
       enddo

    enddo over_parts
    !$omp end parallel do

    if (present(dat)) then
       datold = datnorm
    else
       datold = datsmooth
    endif
    !
    !--normalise dat array
    !
    if (normalise) then
       where (datnorm > 0.)
          datsmooth = datsmooth/datnorm
       end where
    endif

 enddo iterations

 if (present(datpix2)) datpix2 = datnorm
 call wall_time(t2)
 if (t2-t1 > 1.) call print_time(t2-t1)

end subroutine interpolate2D_pixels

!--------------------------------------------------------------------------
!     subroutine to interpolate from even grid of pixels TO arbitrary points
!
!     Input: data points  : x,y    (npart)
!            smoothing length : hh (npart)
!            number of pixels in x,y : npixx,npixy
!            pixel width             : pixwidth
!            data on pixels : datpix(npixx,npixy)
!
!     Output: smoothed data          : dat (npart)
!
!     Written by Daniel Price 2021
!--------------------------------------------------------------------------

subroutine interpolate2D_fromgrid(x,y,hh,dat,gradh,sigma,mask,npart, &
     xmin,ymin,datpix,npixx,npixy,pixwidthx,pixwidthy)

 use timing, only:wall_time,print_time
 integer, intent(in) :: npart,npixx,npixy
 real,    intent(in),  dimension(npart) :: x,y,hh
 real,    intent(out), dimension(npart) :: dat,gradh,sigma
 integer, intent(in),  dimension(npart) :: mask
 real,    intent(in) :: xmin,ymin,pixwidthx,pixwidthy
 real,    intent(in), dimension(npixx,npixy) :: datpix

 real, dimension(npixx) :: dx2i,qq2

 integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
 real :: hi,hi1,radkernx,radkerny,q2,wab,const,datpart
 real :: term,dy,xpix,ypix,xi,yi,dy2,dwabdh,gradhi
 real :: t1,t2

 !print "(1x,a)",'interpolating from 2D pixels to particles...'

 call wall_time(t1)
 if (.not.associated(wfunc)) call select_kernel(0)

 const = cnormk2D  ! normalisation constant

 if (pixwidthx <= 0. .or. pixwidthy <= 0.) then
    print "(1x,a)",'interpolate2D: error: pixel width <= 0'
    return
 endif
 !
 !--loop over particles
 !
 !$omp parallel do default(none) &
 !$omp shared(npart,mask,x,y,hh,xmin,ymin) &
 !$omp shared(dat,gradh,sigma,datpix,npixx,npixy,const,radkernel,radkernel2) &
 !$omp shared(pixwidthx,pixwidthy) &
 !$omp private(i,xi,yi,ipix,jpix,hi,hi1,datpart) &
 !$omp private(radkernx,radkerny,ipixmin,ipixmax,jpixmin,jpixmax) &
 !$omp private(dx2i,xpix,ypix,dy,dy2,q2,wab,term,qq2,dwabdh,gradhi)
 over_parts: do i=1,npart
    !
    !--skip particles with itype < 0
    !
    if (mask(i) < 0) cycle over_parts

    xi = x(i)
    yi = y(i)
    hi = hh(i)  ! in units of pixel spacing

    hi1 = 1./hi
    term = const*pixwidthx*pixwidthy*hi1*hi1 ! this is dx*dy/h**2
    if (term <= 0.) cycle over_parts
    !
    !--set kernel related quantities
    !
    radkernx = radkernel*hi  ! radius of the smoothing kernel
    radkerny = radkernel*hi  ! radius of the smoothing kernel
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((xi - radkernx - xmin)/abs(pixwidthx))
    jpixmin = int((yi - radkerny - ymin)/abs(pixwidthy))
    ipixmax = int((xi + radkernx - xmin)/abs(pixwidthx)) + 1
    jpixmax = int((yi + radkerny - ymin)/abs(pixwidthy)) + 1

    if (ipixmin < 1)     ipixmin = 1      ! make sure we only receive contributions
    if (ipixmax > npixx) ipixmax = npixx  ! from pixels in the image
    if (jpixmin < 1)     jpixmin = 1
    if (jpixmax > npixy) jpixmax = npixy
    !
    !--precalculate an array of dx2 for this particle (optimisation)
    !
    do ipix=ipixmin,ipixmax
       xpix = xmin + (ipix-0.5)*pixwidthx
       dx2i(ipix) = ((xpix - xi)**2)*hi1*hi1
    enddo
    !
    !--loop over pixels, adding the contribution TO this particle
    !
    datpart = 0.
    gradhi = 0.
    do jpix = jpixmin,jpixmax
       ypix = ymin + (jpix-0.5)*pixwidthy
       dy = ypix - yi
       dy2 = dy*dy*hi1*hi1

       do ipix = ipixmin,ipixmax
          q2 = dx2i(ipix) + dy2
          !
          !--SPH kernel
          !
          if (q2 < radkernel2) then
             wab = wkernel(q2)
             !
             !--interpolate to particle using the summation interpolant
             !
             datpart = datpart + datpix(ipix,jpix)*wab
             dwabdh = -(2.*wab + sqrt(q2)*dwkernel(q2))
             gradhi = gradhi + datpix(ipix,jpix)*dwabdh
          endif
       enddo
    enddo
    dat(i) = term*datpart
    gradh(i) = term*gradhi*hi1

    datpart = 0.
    do jpix = jpixmin,jpixmax
       ypix = ymin + (jpix-0.5)*pixwidthy
       dy = ypix - yi
       dy2 = dy*dy*hi1*hi1

       do ipix = ipixmin,ipixmax
          q2 = dx2i(ipix) + dy2
          !
          !--SPH kernel
          !
          if (q2 < radkernel2) then
             wab = wkernel(q2)
             !
             !--interpolate to particle using the summation interpolant
             !
             datpart = datpart + (datpix(ipix,jpix) - dat(i))**2*wab
          endif
       enddo
    enddo
    sigma(i) = sqrt(term*datpart)

 enddo over_parts
 !$omp end parallel do

 call wall_time(t2)
 if (t2-t1 > 1.) call print_time(t2-t1)

end subroutine interpolate2D_fromgrid

!------------------------------------------------------------
! interface to kernel routine to avoid problems with openMP
!-----------------------------------------------------------
real function wkernel(q2)
 use kernels, only:wfunc
 real, intent(in) :: q2

 wkernel = wfunc(q2)

end function wkernel

real function dwkernel(q2)
 use kernels, only:dwfunc
 real, intent(in) :: q2

 dwkernel = dwfunc(q2)

end function dwkernel

end module interpolations2D
