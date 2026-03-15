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
!  Module containing all of the routines required for interpolation
!  from 3D data to a 3D grid in arbitrary coordinate systems (SLOW!)
!
!----------------------------------------------------------------------

module interpolations3Dgeom
 use kernels,       only:radkernel2,radkernel,cnormk3D
 use geometry,      only:labelcoordsys,coord_is_length,igeom_cartesian,coord_transform
 use interpolation, only:doub_prec,iroll
 implicit none
 public :: interpolate3Dgeom,interpolate3Dgeom_vec

contains
!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is interpolated according to the formula
!
!     datsmooth(pixel) = sum_b weight_b dat_b W(r-r_b, h_b)
!
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline.
!
!     For a standard SPH smoothing the weight function for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     this version is written for slices through a rectangular volume, ie.
!     assumes a uniform pixel size in x,y, whilst the number of pixels
!     in the z direction can be set to the number of cross-section slices.
!
!     Input: particle coordinates  : x,y,z (npart)
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npix(1),npix(2),npix(3))
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!     Revised for "splash to grid", Monash University 02/11/09
!     Version for arbitrary coordinate systems 02/05/18
!--------------------------------------------------------------------------

subroutine interpolate3Dgeom(igeom,x,y,z,hh,weight,dat,itype,npart,&
     xmin,datsmooth,npix,pixwidth,xorigin,normalise,periodic)
 integer, intent(in) :: igeom,npart,npix(3)
 real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin(3),pixwidth(3),xorigin(3)
 real(doub_prec), intent(out), dimension(npix(1),npix(2),npix(3)) :: datsmooth
 logical, intent(in) :: normalise,periodic(3)
 real(doub_prec), dimension(npix(1),npix(2),npix(3)) :: datnorm

 integer :: i,ipix,jpix,kpix,ierr
 integer :: iprintinterval,iprintnext
 integer :: ipixmin(3),ipixmax(3)
 integer :: ipixi,jpixi,kpixi,nsub,isub,jsub,ksub,nsubpix
 integer :: maxsubpix,meansubpix,ncount
 real :: xminpix(3),hmin,dsubpix(3)
 real :: xi(3),xci(3),xcoord(3),hi,hi1,hi21,radkern,wab,q2,const
 real :: term,termnorm,xpix(3),dx(3)
 !real :: t_start,t_end
 logical :: iprintprogress
 !$ integer :: omp_get_num_threads
 integer(kind=selected_int_kind(10)) :: iprogress  ! up to 10 digits

 datsmooth = 0.
 datnorm = 0.
 if (normalise) then
    print "(1x,a)",'interpolating from particles to 3D '//trim(labelcoordsys(igeom))//' grid (normalised) ...'
 else
    print "(1x,a)",'interpolating from particles to 3D '//trim(labelcoordsys(igeom))//' grid (non-normalised) ...'
 endif
 if (any(pixwidth <= 0.)) then
    print "(1x,a)",'interpolate3D: error: pixel width <= 0'
    return
 endif
 if (any(hh(1:npart) <= tiny(hh))) then
    print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
 endif

 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = (npart  >=  100000) .or. (npix(1)*npix(2)  > 100000)
 !$ iprintprogress = .false.
 !
 !--loop over particles
 !
 iprintinterval = 25
 if (npart >= 1e6) iprintinterval = 10
 iprintnext = iprintinterval
 !
 !--get starting CPU time
 !
 !call cpu_time(t_start)

 xminpix(:) = xmin(:) - 0.5*pixwidth(:)
 !
 !--use a minimum smoothing length on the grid to make
 !  sure that particles contribute to at least one pixel
 !
 hmin = 0.
 do i=1,3
    if (coord_is_length(i,igeom)) hmin = max(hmin,0.5*pixwidth(i))
 enddo

 const = cnormk3D  ! normalisation constant (3D)
 nsub = 0
 maxsubpix = 0
 meansubpix = 0
 ncount = 0
 !
 !--loop over particles
 !
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,radkernel,radkernel2) &
!$omp shared(xminpix,pixwidth,xorigin) &
!$omp shared(npix,const,igeom) &
!$omp shared(datnorm,normalise,periodic) &
!$omp shared(iprintprogress,iprintinterval) &
!$omp shared(hmin) &
!$omp private(hi,xi,xci,xcoord,xpix,radkern,hi1,hi21) &
!$omp private(term,termnorm) &
!$omp private(iprogress) &
!$omp firstprivate(iprintnext) &
!$omp private(ipixmin,ipixmax,ierr) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx,q2,wab,isub,jsub,ksub,nsubpix,dsubpix) &
!$omp reduction(+:nsub,meansubpix,ncount) &
!$omp reduction(max:maxsubpix)
!$omp master
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end master

!$omp do schedule (guided, 2)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
    if (iprintprogress) then
       iprogress = 100*i/npart
       if (iprogress >= iprintnext) then
          write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
          iprintnext = iprintnext + iprintinterval
       endif
    endif
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle over_parts

    hi = hh(i)
    if (hi <= 0.) cycle over_parts
    hi1 = 1./hi
    hi21 = hi1*hi1
    radkern = radkernel*hi   ! radius of the smoothing kernel
    termnorm = const*weight(i)
    term = termnorm*dat(i)
    !
    !--set kernel related quantities
    !
    xci(1) = x(i) + xorigin(1)  ! xci = position in cartesians
    xci(2) = y(i) + xorigin(2)
    xci(3) = z(i) + xorigin(3)
    call get_pixel_limits(xci,xi,radkern,ipixmin,ipixmax,npix,pixwidth,xmin,periodic,igeom,ierr)
    if (ierr /= 0) cycle over_parts

    if (hi < hmin) then
       nsub = nsub + 1
       nsubpix = max(2, int(hmin/hi) + 1)
    else
       nsubpix = 1
    endif
    ! statistics on sub-pixel interpolation
    maxsubpix = max(maxsubpix,nsubpix)
    meansubpix = meansubpix + nsubpix
    ncount = ncount + 1

    dsubpix(:) = pixwidth(:) / real(nsubpix)
    !
    !--loop over pixels, adding the contribution from this particle
    !  When nsubpix > 1 we sub-divide each pixel and average the kernel
    !  so that small kernels still contribute to the coarse grid.
    !
    do kpix = ipixmin(3),ipixmax(3)
       kpixi = kpix
       if (periodic(3)) kpixi = iroll(kpix,npix(3))

       do jpix = ipixmin(2),ipixmax(2)
          jpixi = jpix
          if (periodic(2)) jpixi = iroll(jpix,npix(2))

          do ipix = ipixmin(1),ipixmax(1)
             ipixi = ipix
             if (periodic(1)) ipixi = iroll(ipix,npix(1))

             wab = 0.
             do ksub = 1,nsubpix
                xcoord(3) = xmin(3) + (kpix-1)*pixwidth(3) + (ksub-0.5)*dsubpix(3)
                do jsub = 1,nsubpix
                   xcoord(2) = xmin(2) + (jpix-1)*pixwidth(2) + (jsub-0.5)*dsubpix(2)
                   do isub = 1,nsubpix
                      xcoord(1) = xmin(1) + (ipix-1)*pixwidth(1) + (isub-0.5)*dsubpix(1)
                      call coord_transform(xcoord,3,igeom,xpix,3,igeom_cartesian)
                      dx = xpix(:) - xci(:)
                      q2 = (dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))*hi21
                      if (q2 < radkernel2) wab = wab + wkernel(q2)
                   enddo
                enddo
             enddo
             wab = wab / real(nsubpix**3)
             if (wab > 0.) then
                !$omp atomic
                datsmooth(ipixi,jpixi,kpixi) = datsmooth(ipixi,jpixi,kpixi) + term*wab
                if (normalise) then
                   !$omp atomic
                   datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                endif
             endif
          enddo
       enddo
    enddo
 enddo over_parts
!$omp enddo
!$omp end parallel

 !
 !--normalise dat array
 !
 if (normalise) then
    where (datnorm > tiny(datnorm))
       datsmooth = datsmooth/datnorm
    end where
 endif
 print "(20x,a,i0,a,i0,a,f6.1/)",'particles on subgrid: ',nsub,&
       ' max subpixels: ',maxsubpix,' mean subpixels: ',meansubpix/real(ncount)

end subroutine interpolate3Dgeom

!------------------------------------------------------------------------------
!
! Same functionality as interpolate3Dgeom, but for vector quantities
!
!------------------------------------------------------------------------------
subroutine interpolate3Dgeom_vec(igeom,x,y,z,hh,weight,datvec,itype,npart,&
     xmin,datsmooth,npix,pixwidth,xorigin,normalise,periodic)
 integer, intent(in) :: igeom,npart,npix(3)
 real, intent(in), dimension(npart) :: x,y,z,hh,weight
 real, intent(in), dimension(npart,3)  :: datvec
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin(3),pixwidth(3),xorigin(3)
 real(doub_prec), intent(out), dimension(3,npix(1),npix(2),npix(3)) :: datsmooth
 logical, intent(in) :: normalise,periodic(3)
 real(doub_prec), dimension(npix(1),npix(2),npix(3)) :: datnorm

 integer :: i,ipix,jpix,kpix,ierr
 integer :: iprintinterval,iprintnext
 integer :: ipixmin(3),ipixmax(3)
 integer :: ipixi,jpixi,kpixi,isub,jsub,ksub,nsubpix
 real :: xminpix(3),hmin,dsubpix(3)
 real :: xi(3),xci(3),xcoord(3),hi,hi1,hi21,radkern,wab,q2,const
 real :: term(3),termnorm,xpix(3),dx(3),ddatnorm
 !real :: t_start,t_end
 logical :: iprintprogress
 !$ integer :: omp_get_num_threads
 integer(kind=selected_int_kind(10)) :: iprogress  ! up to 10 digits

 datsmooth = 0.
 datnorm = 0.
 if (normalise) then
    print "(1x,a)",'interpolating to 3D '//trim(labelcoordsys(igeom))//' grid (normalised) ...'
 else
    print "(1x,a)",'interpolating to 3D '//trim(labelcoordsys(igeom))//' grid (non-normalised) ...'
 endif
 if (any(pixwidth <= 0.)) then
    print "(1x,a)",'interpolate3D: error: pixel width <= 0'
    return
 endif
 if (any(hh(1:npart) <= tiny(hh))) then
    print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
 endif

 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = (npart  >=  100000) .or. (npix(1)*npix(2)  > 100000)
 !$ iprintprogress = .false.
 !
 !--loop over particles
 !
 iprintinterval = 25
 if (npart >= 1e6) iprintinterval = 10
 iprintnext = iprintinterval
 !
 !--get starting CPU time
 !
 !call cpu_time(t_start)

 xminpix(:) = xmin(:) - 0.5*pixwidth(:)
 !
 !--use a minimum smoothing length on the grid to make
 !  sure that particles contribute to at least one pixel
 !
 hmin = 0.
 do i=1,3
    if (coord_is_length(i,igeom)) hmin = max(hmin,0.5*pixwidth(i))
 enddo

 const = cnormk3D  ! normalisation constant (3D)
 !
 !--loop over particles
 !
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,datvec,itype,datsmooth,npart) &
!$omp shared(xmin,radkernel,radkernel2) &
!$omp shared(xminpix,pixwidth,xorigin) &
!$omp shared(npix,const,igeom) &
!$omp shared(datnorm,normalise,periodic) &
!$omp shared(hmin) &
!$omp shared(iprintprogress,iprintinterval) &
!$omp private(hi,xi,xci,xcoord,xpix,radkern,hi1,hi21) &
!$omp private(term,termnorm) &
!$omp private(ipixmin,ipixmax,ierr) &
!$omp private(iprogress) &
!$omp firstprivate(iprintnext) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx,q2,wab,isub,jsub,ksub,nsubpix,dsubpix)
!$omp master
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end master

!$omp do schedule (guided, 2)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
    if (iprintprogress) then
       iprogress = 100*i/npart
       if (iprogress >= iprintnext) then
          write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
          iprintnext = iprintnext + iprintinterval
       endif
    endif
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle over_parts

    hi = hh(i)
    if (hi <= 0.) cycle over_parts
    termnorm = const*weight(i)
    hi1 = 1./hi
    hi21 = hi1*hi1
    radkern = radkernel*hi   ! radius of the smoothing kernel
    term(:) = termnorm*datvec(i,:)
    !
    !--set kernel related quantities
    !
    xci(1) = x(i) + xorigin(1)  ! xci = position in cartesians
    xci(2) = y(i) + xorigin(2)
    xci(3) = z(i) + xorigin(3)
    call get_pixel_limits(xci,xi,radkern,ipixmin,ipixmax,npix,pixwidth,xmin,periodic,igeom,ierr)

    if (ierr /= 0) cycle over_parts

    if (hi < hmin) then
       nsubpix = max(2, int(hmin/hi) + 1)
    else
       nsubpix = 1
    endif
    dsubpix(:) = pixwidth(:) / real(nsubpix)
    !
    !--loop over pixels, adding the contribution from this particle
    !
    do kpix = ipixmin(3),ipixmax(3)
       kpixi = kpix
       if (periodic(3)) kpixi = iroll(kpix,npix(3))

       do jpix = ipixmin(2),ipixmax(2)
          jpixi = jpix
          if (periodic(2)) jpixi = iroll(jpix,npix(2))

          do ipix = ipixmin(1),ipixmax(1)
             ipixi = ipix
             if (periodic(1)) ipixi = iroll(ipix,npix(1))

             wab = 0.
             do ksub = 1,nsubpix
                xcoord(3) = xmin(3) + (kpix-1)*pixwidth(3) + (ksub-0.5)*dsubpix(3)
                do jsub = 1,nsubpix
                   xcoord(2) = xmin(2) + (jpix-1)*pixwidth(2) + (jsub-0.5)*dsubpix(2)
                   do isub = 1,nsubpix
                      xcoord(1) = xmin(1) + (ipix-1)*pixwidth(1) + (isub-0.5)*dsubpix(1)
                      call coord_transform(xcoord,3,igeom,xpix,3,igeom_cartesian)
                      dx = xpix(:) - xci(:)
                      q2 = (dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))*hi21
                      if (q2 < radkernel2) wab = wab + wkernel(q2)
                   enddo
                enddo
             enddo
             wab = wab / real(nsubpix**3)
             if (wab > 0.) then
                !$omp atomic
                datsmooth(1,ipixi,jpixi,kpixi) = datsmooth(1,ipixi,jpixi,kpixi) + term(1)*wab
                !$omp atomic
                datsmooth(2,ipixi,jpixi,kpixi) = datsmooth(2,ipixi,jpixi,kpixi) + term(2)*wab
                !$omp atomic
                datsmooth(3,ipixi,jpixi,kpixi) = datsmooth(3,ipixi,jpixi,kpixi) + term(3)*wab
                if (normalise) then
                   !$omp atomic
                   datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                endif
             endif
          enddo
       enddo
    enddo
 enddo over_parts
!$omp enddo
!$omp end parallel

 !
 !--normalise dat array
 !
 if (normalise) then
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(datsmooth,datnorm,npix) &
    !$omp private(kpix,jpix,ipix,ddatnorm)
    do kpix=1,npix(3)
       do jpix=1,npix(2)
          do ipix=1,npix(1)
             if (datnorm(ipix,jpix,kpix) > tiny(datnorm)) then
                ddatnorm = 1./datnorm(ipix,jpix,kpix)
                datsmooth(1,ipix,jpix,kpix) = datsmooth(1,ipix,jpix,kpix)*ddatnorm
                datsmooth(2,ipix,jpix,kpix) = datsmooth(2,ipix,jpix,kpix)*ddatnorm
                datsmooth(3,ipix,jpix,kpix) = datsmooth(3,ipix,jpix,kpix)*ddatnorm
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
 endif

end subroutine interpolate3Dgeom_vec

!------------------------------------------------------------
! interface to kernel routine to avoid problems with openMP
!-----------------------------------------------------------
real function wkernel(q2)
 use kernels, only:wfunc
 real, intent(in) :: q2

 wkernel = wfunc(q2)

end function wkernel

!--------------------------------------------------------------------------
!
!  utility routine for use in above routines
!  IN:
!    xci - coordinates of particle in cartesians
!    radkern - radius of kernel (e.g. 2h)
!  OUT:
!    xi - non-cartesian coordinates of particle (e.g. r, phi, z)
!    ipixmin,ipixmax,jpixmin,jpixmax - pixel limits
!
!--------------------------------------------------------------------------
subroutine get_pixel_limits(xci,xi,radkern,ipixmin,ipixmax,npix,pixwidth,xmin,periodic,igeom,ierr)
 use geometry, only:coord_is_periodic,get_coord_limits
 real, intent(in)  :: xci(3),radkern,pixwidth(3),xmin(3)
 real, intent(out) :: xi(3)
 integer, intent(out) :: ipixmin(3),ipixmax(3),ierr
 integer, intent(in)  :: igeom,npix(3)
 logical, intent(in)  :: periodic(3)
 real :: xpixmin(3),xpixmax(3)
 integer :: i

 ierr = 0
 !
 !--get limits of rendering in new coordinate system
 !
 call get_coord_limits(radkern,xci,xi,xpixmin,xpixmax,igeom)
 !print*,' R min = ',xpixmin(1),' Rmax = ',xpixmax(1)
 !
 !--now work out contributions to pixels in the transformed space.
 !  For periodic dimensions use the bounding-box range (may extend
 !  outside [1,npix]); the caller wraps with iroll().  Clamp only
 !  if the range exceeds a full period.
 !
 do i = 1, 3
    ipixmin(i) = int((xpixmin(i) - xmin(i))/pixwidth(i)) + 1
    ipixmax(i) = int((xpixmax(i) - xmin(i))/pixwidth(i)) + 1
    if (periodic(i)) then
       if (ipixmax(i) - ipixmin(i) >= npix(i)) then
          ipixmin(i) = 1
          ipixmax(i) = npix(i)
       endif
    else
       if (ipixmin(i) < 1)       ipixmin(i) = 1
       if (ipixmax(i) > npix(i)) ipixmax(i) = npix(i)
       if (ipixmin(i) > npix(i) .or. ipixmax(i) < 1) ierr = ierr + 1
    endif
 enddo

end subroutine get_pixel_limits

end module interpolations3Dgeom
