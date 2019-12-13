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
!  from 3D data to a 3D grid (SLOW!)
!
!----------------------------------------------------------------------

module interpolations3D
 use kernels,       only:radkernel2,radkernel,cnormk3D
 use interpolation, only:iroll
 implicit none
 integer, parameter :: doub_prec = kind(0.d0)
 public :: interpolate3D,interpolate3D_vec

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
!     Output: smoothed data            : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!     Revised for "splash to grid", Monash University 02/11/09
!     Maya Petkova contributed exact subgrid interpolation, April 2019
!--------------------------------------------------------------------------

subroutine interpolate3D(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidthx,pixwidthy,pixwidthz,&
     normalise,periodicx,periodicy,periodicz)

 integer, intent(in) :: npart,npixx,npixy,npixz
 real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
 real(doub_prec), intent(out), dimension(npixx,npixy,npixz) :: datsmooth
 logical, intent(in) :: normalise,periodicx,periodicy,periodicz
 real(doub_prec), dimension(npixx,npixy,npixz) :: datnorm

 integer :: i,ipix,jpix,kpix
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nwarn
 real :: xminpix,yminpix,zminpix,hmin !,dhmin3
 real, dimension(npixx) :: dx2i
 real :: xi,yi,zi,hi,hi1,hi21,radkern,wab,q2,const,dyz2,dz2
 real :: term,termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax,dfac
 !real :: t_start,t_end
 logical :: iprintprogress

! Exact rendering
 real :: pixint, wint
 logical, parameter :: exact_rendering = .true.   ! use exact rendering y/n
 integer :: usedpart, negflag

#ifdef _OPENMP
 integer :: omp_get_num_threads
#else
 integer(kind=selected_int_kind(10)) :: iprogress  ! up to 10 digits
#endif

 datsmooth = 0.
 datnorm = 0.
 if (exact_rendering) then
    print "(1x,a)",'interpolating to 3D grid (exact/Petkova+2018 on subgrid) ...'
 elseif (normalise) then
    print "(1x,a)",'interpolating to 3D grid (normalised) ...'
 else
    print "(1x,a)",'interpolating to 3D grid (non-normalised) ...'
 endif
 if (pixwidthx <= 0. .or. pixwidthy <= 0 .or. pixwidthz <= 0) then
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
 iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000) .or. exact_rendering
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

 usedpart = 0

 xminpix = xmin - 0.5*pixwidthx
 yminpix = ymin - 0.5*pixwidthy
 zminpix = zmin - 0.5*pixwidthz
 pixwidthmax = max(pixwidthx,pixwidthy,pixwidthz)
 !
 !--use a minimum smoothing length on the grid to make
 !  sure that particles contribute to at least one pixel
 !
 hmin = 0.5*pixwidthmax
 !dhmin3 = 1./(hmin*hmin*hmin)

 const = cnormk3D  ! normalisation constant (3D)
 nwarn = 0
 !
 !--loop over particles
 !
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
!$omp shared(xminpix,yminpix,zminpix,pixwidthx,pixwidthy,pixwidthz) &
!$omp shared(npixx,npixy,npixz,const) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
!$omp shared(hmin,pixwidthmax) & !,dhmin3) &
!$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
!$omp private(term,termnorm,xpixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
!$omp private(pixint,wint,negflag,dfac) &
!$omp reduction(+:nwarn,usedpart)
!$omp master
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end master

!$omp do schedule (guided, 2)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
#ifndef _OPENMP
    if (iprintprogress) then
       iprogress = 100*i/npart
       if (iprogress >= iprintnext) then
          write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
          iprintnext = iprintnext + iprintinterval
       endif
    endif
#endif
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0) cycle over_parts

    hi = hh(i)
    if (hi <= 0.) then
       cycle over_parts
    elseif (hi < hmin) then
       !
       !--use minimum h to capture subgrid particles
       !  (get better results *without* adjusting weights)
       !
       termnorm = const*weight(i) !*(hi*hi*hi)*dhmin3
       if (.not.exact_rendering) hi = hmin
    else
       termnorm = const*weight(i)
    endif

    !
    !--set kernel related quantities
    !
    xi = x(i)
    yi = y(i)
    zi = z(i)

    hi1 = 1./hi
    hi21 = hi1*hi1
    radkern = radkernel*hi   ! radius of the smoothing kernel
    !termnorm = const*weight(i)
    term = termnorm*dat(i)
    dfac = hi**3/(pixwidthx*pixwidthy*pixwidthz*const)
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((xi - radkern - xmin)/pixwidthx)
    jpixmin = int((yi - radkern - ymin)/pixwidthy)
    kpixmin = int((zi - radkern - zmin)/pixwidthz)
    ipixmax = int((xi + radkern - xmin)/pixwidthx) + 1
    jpixmax = int((yi + radkern - ymin)/pixwidthy) + 1
    kpixmax = int((zi + radkern - zmin)/pixwidthz) + 1

    if (.not.periodicx) then
       if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
       if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
    endif
    if (.not.periodicy) then
       if (jpixmin < 1)     jpixmin = 1
       if (jpixmax > npixy) jpixmax = npixy
    endif
    if (.not.periodicz) then
       if (kpixmin < 1)     kpixmin = 1
       if (kpixmax > npixz) kpixmax = npixz
    endif

    negflag = 0

    !
    !--precalculate an array of dx2 for this particle (optimisation)
    !
    nxpix = 0
    do ipix=ipixmin,ipixmax
       nxpix = nxpix + 1
       ipixi = ipix
       if (periodicx) ipixi = iroll(ipix,npixx)
       xpixi = xminpix + ipix*pixwidthx
       !--watch out for errors with periodic wrapping...
       if (nxpix <= size(dx2i)) then
          dx2i(nxpix) = ((xpixi - xi)**2)*hi21
       endif
    enddo

    !--if particle contributes to more than npixx pixels
    !  (i.e. periodic boundaries wrap more than once)
    !  truncate the contribution and give warning
    if (nxpix > npixx) then
       nwarn = nwarn + 1
       ipixmax = ipixmin + npixx - 1
    endif
    !
    !--loop over pixels, adding the contribution from this particle
    !
    do kpix = kpixmin,kpixmax
       kpixi = kpix
       if (periodicz) kpixi = iroll(kpix,npixz)

       zpix = zminpix + kpix*pixwidthz
       dz = zpix - zi
       dz2 = dz*dz*hi21

       do jpix = jpixmin,jpixmax
          jpixi = jpix
          if (periodicy) jpixi = iroll(jpix,npixy)

          ypix = yminpix + jpix*pixwidthy
          dy = ypix - yi
          dyz2 = dy*dy*hi21 + dz2

          nxpix = 0
          do ipix = ipixmin,ipixmax
             if ((kpix==kpixmin).and.(jpix==jpixmin).and.(ipix==ipixmin)) then
                usedpart = usedpart + 1
             endif

             nxpix = nxpix + 1
             ipixi = ipix
             if (periodicx) ipixi = iroll(ipix,npixx)

             q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21

             if (exact_rendering .and. ipixmax-ipixmin <= 4) then
                if (q2 < radkernel2 + 3.*pixwidthmax**2*hi21) then
                   xpixi = xminpix + ipix*pixwidthx

                   ! Contribution of the cell walls in the xy-plane
                   pixint = 0.0
                   wint = wallint(zpix-zi+0.5*pixwidthz,xi,yi,xpixi,ypix,pixwidthx,pixwidthy,hi)
                   pixint = pixint + wint

                   wint = wallint(zi-zpix+0.5*pixwidthz,xi,yi,xpixi,ypix,pixwidthx,pixwidthy,hi)
                   pixint = pixint + wint

                   ! Contribution of the cell walls in the xz-plane
                   wint = wallint(ypix-yi+0.5*pixwidthy,xi,zi,xpixi,zpix,pixwidthx,pixwidthz,hi)
                   pixint = pixint + wint

                   wint = wallint(yi-ypix+0.5*pixwidthy,xi,zi,xpixi,zpix,pixwidthx,pixwidthz,hi)
                   pixint = pixint + wint

                   ! Contribution of the cell walls in the yz-plane
                   wint = wallint(xpixi-xi+0.5*pixwidthx,zi,yi,zpix,ypix,pixwidthz,pixwidthy,hi)
                   pixint = pixint + wint

                   wint = wallint(xi-xpixi+0.5*pixwidthx,zi,yi,zpix,ypix,pixwidthz,pixwidthy,hi)
                   pixint = pixint + wint

                   wab = pixint*dfac ! /(pixwidthx*pixwidthy*pixwidthz*const)*hi**3

                   if (pixint < -0.01d0) then
                      print*, "Error: (",ipixi,jpixi,kpixi,") -> ", pixint, term*wab
                   endif

                   !
                   !--calculate data value at this pixel using the summation interpolant
                   !
                   !$omp atomic
                   datsmooth(ipixi,jpixi,kpixi) = datsmooth(ipixi,jpixi,kpixi) + term*wab
                   if (normalise) then
                      !$omp atomic
                      datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                   endif
                endif
             else
                if (q2 < radkernel2) then
                   !
                   !--SPH kernel - standard cubic spline
                   !
                   wab = wkernel(q2)
                   !
                   !--calculate data value at this pixel using the summation interpolant
                   !
                   !$omp atomic
                   datsmooth(ipixi,jpixi,kpixi) = datsmooth(ipixi,jpixi,kpixi) + term*wab
                   if (normalise) then
                      !$omp atomic
                      datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                   endif
                endif
             endif
          enddo
       enddo
    enddo
 enddo over_parts
!$omp enddo
!$omp end parallel

 if (nwarn > 0) then
    print "(a,i11,a,/,a)",' interpolate3D: WARNING: contributions truncated from ',nwarn,' particles',&
                           '                that wrap periodic boundaries more than once'
 endif
 !
 !--normalise dat array
 !
 if (normalise) then
    where (datnorm > tiny(datnorm))
       datsmooth = datsmooth/datnorm
    end where
 endif

 print*, 'Number of particles in the volume: ', usedpart

end subroutine interpolate3D

subroutine interpolate3D_vec(x,y,z,hh,weight,datvec,itype,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidthx,pixwidthy,pixwidthz,&
     normalise,periodicx,periodicy,periodicz)

 integer, intent(in) :: npart,npixx,npixy,npixz
 real, intent(in), dimension(npart)    :: x,y,z,hh,weight
 real, intent(in), dimension(npart,3)  :: datvec
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
 real(doub_prec), intent(out), dimension(3,npixx,npixy,npixz) :: datsmooth
 logical, intent(in) :: normalise,periodicx,periodicy,periodicz
 real(doub_prec), dimension(npixx,npixy,npixz) :: datnorm

 integer :: i,ipix,jpix,kpix
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nwarn
 real :: xminpix,yminpix,zminpix
 real, dimension(npixx) :: dx2i
 real :: xi,yi,zi,hi,hi1,hi21,radkern,wab,q2,const,dyz2,dz2
 real :: termnorm,dy,dz,ypix,zpix,xpixi,ddatnorm
 real, dimension(3) :: term
 !real :: t_start,t_end
 logical :: iprintprogress
#ifdef _OPENMP
 integer :: omp_get_num_threads
#else
 integer(kind=selected_int_kind(10)) :: iprogress  ! up to 10 digits
#endif

 datsmooth = 0.
 datnorm = 0.
 if (normalise) then
    print "(1x,a)",'interpolating to 3D grid (normalised) ...'
 else
    print "(1x,a)",'interpolating to 3D grid (non-normalised) ...'
 endif
 if (pixwidthx <= 0. .or. pixwidthy <= 0. .or. pixwidthz <= 0.) then
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
 iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000)
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

 xminpix = xmin - 0.5*pixwidthx
 yminpix = ymin - 0.5*pixwidthy
 zminpix = zmin - 0.5*pixwidthz

 const = cnormk3D  ! normalisation constant (3D)
 nwarn = 0

!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,datvec,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
!$omp shared(xminpix,yminpix,zminpix,pixwidthx,pixwidthy,pixwidthz) &
!$omp shared(npixx,npixy,npixz,const) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
!$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
!$omp private(term,termnorm,xpixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
!$omp reduction(+:nwarn)
!$omp master
!$  print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end master
 !
 !--loop over particles
 !
!$omp do schedule (guided, 2)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
#ifndef _OPENMP
    if (iprintprogress) then
       iprogress = 100*i/npart
       if (iprogress >= iprintnext) then
          write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
          iprintnext = iprintnext + iprintinterval
       endif
    endif
#endif
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0) cycle over_parts

    hi = hh(i)
    if (hi <= 0.) cycle over_parts

    !
    !--set kernel related quantities
    !
    xi = x(i)
    yi = y(i)
    zi = z(i)

    hi1 = 1./hi
    hi21 = hi1*hi1
    radkern = radkernel*hi   ! radius of the smoothing kernel
    termnorm = const*weight(i)
    term(:) = termnorm*datvec(i,:)
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((xi - radkern - xmin)/pixwidthx)
    jpixmin = int((yi - radkern - ymin)/pixwidthy)
    kpixmin = int((zi - radkern - zmin)/pixwidthz)
    ipixmax = int((xi + radkern - xmin)/pixwidthx) + 1
    jpixmax = int((yi + radkern - ymin)/pixwidthy) + 1
    kpixmax = int((zi + radkern - zmin)/pixwidthz) + 1

    if (.not.periodicx) then
       if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
       if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
    endif
    if (.not.periodicy) then
       if (jpixmin < 1)     jpixmin = 1
       if (jpixmax > npixy) jpixmax = npixy
    endif
    if (.not.periodicz) then
       if (kpixmin < 1)     kpixmin = 1
       if (kpixmax > npixz) kpixmax = npixz
    endif
    !
    !--precalculate an array of dx2 for this particle (optimisation)
    !
    nxpix = 0
    do ipix=ipixmin,ipixmax
       nxpix = nxpix + 1
       ipixi = ipix
       if (periodicx) ipixi = iroll(ipix,npixx)
       xpixi = xminpix + ipix*pixwidthx
       !--watch out for errors with perioic wrapping...
       if (nxpix <= size(dx2i)) then
          dx2i(nxpix) = ((xpixi - xi)**2)*hi21
       endif
    enddo

    !--if particle contributes to more than npixx pixels
    !  (i.e. periodic boundaries wrap more than once)
    !  truncate the contribution and give warning
    if (nxpix > npixx) then
       nwarn = nwarn + 1
       ipixmax = ipixmin + npixx - 1
    endif
    !
    !--loop over pixels, adding the contribution from this particle
    !
    do kpix = kpixmin,kpixmax
       kpixi = kpix
       if (periodicz) kpixi = iroll(kpix,npixz)
       zpix = zminpix + kpix*pixwidthz
       dz = zpix - zi
       dz2 = dz*dz*hi21

       do jpix = jpixmin,jpixmax
          jpixi = jpix
          if (periodicy) jpixi = iroll(jpix,npixy)
          ypix = yminpix + jpix*pixwidthy
          dy = ypix - yi
          dyz2 = dy*dy*hi21 + dz2

          nxpix = 0
          do ipix = ipixmin,ipixmax
             ipixi = ipix
             if (periodicx) ipixi = iroll(ipix,npixx)
             nxpix = nxpix + 1
             q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
             !
             !--SPH kernel - standard cubic spline
             !
             if (q2 < radkernel2) then
                wab = wkernel(q2)
                !
                !--calculate data value at this pixel using the summation interpolant
                !
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

 if (nwarn > 0) then
    print "(a,i11,a,/,a)",' interpolate3D: WARNING: contributions truncated from ',nwarn,' particles',&
                           '                that wrap periodic boundaries more than once'
 endif
 !
 !--normalise dat array
 !
 if (normalise) then
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(datsmooth,datnorm,npixz,npixy,npixx) &
    !$omp private(kpix,jpix,ipix,ddatnorm)
    do kpix=1,npixz
       do jpix=1,npixy
          do ipix=1,npixx
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

 return

end subroutine interpolate3D_vec

!------------------------------------------------------------
! interface to kernel routine to avoid problems with openMP
!-----------------------------------------------------------
real function wkernel(q2)
 use kernels, only:wfunc
 real, intent(in) :: q2

 wkernel = wfunc(q2)

end function wkernel

!------------------------------------------------------------
! functions to evaluate exact overlap of kernel with wall boundaries
! see Petkova, Laibe & Bonnell (2018), J. Comp. Phys
!------------------------------------------------------------
real function wallint(r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi)
 real, intent(in) :: r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi
 real(doub_prec) :: R_0, d1, d2, dx, dy, h

 wallint = 0.0
 dx = xc - xp
 dy = yc - yp
 h = hi

 !
 ! Contributions from each of the 4 sides of a cell wall
 !
 R_0 = 0.5*pixwidthy + dy
 d1 = 0.5*pixwidthx - dx
 d2 = 0.5*pixwidthx + dx
 wallint = wallint + pint(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthy - dy
 d1 = 0.5*pixwidthx + dx
 d2 = 0.5*pixwidthx - dx
 wallint = wallint + pint(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthx + dx
 d1 = 0.5*pixwidthy + dy
 d2 = 0.5*pixwidthy - dy
 wallint = wallint + pint(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthx - dx
 d1 = 0.5*pixwidthy - dy
 d2 = 0.5*pixwidthy + dy
 wallint = wallint + pint(r0, R_0, d1, d2, h)

end function wallint


real function pint(r0, R_0, d1, d2, hi)

 real(doub_prec), intent(in) :: R_0, d1, d2, hi
 real, intent(in) :: r0
 real(doub_prec) :: ar0, aR_0
 real(doub_prec) :: int1, int2
 integer :: fflag = 0

 if (abs(r0) < tiny(0.)) then
    pint = 0.d0
    return
 endif

 if (r0  >  0.d0) then
    pint = 1.d0
    ar0 = r0
 else
    pint = -1.d0
    ar0 = -r0
 endif

 if (R_0  >  0.d0) then
    aR_0 = R_0
 else
    pint = -pint
    aR_0 = -R_0
 endif

 int1 = full_integral_3D(d1, ar0, aR_0, hi)
 int2 = full_integral_3D(d2, ar0, aR_0, hi)

 if (int1 < 0.d0) int1 = 0.d0
 if (int2 < 0.d0) int2 = 0.d0

 if (d1*d2  >=  0) then
    pint = pint*(int1 + int2)
    if (int1 + int2 < 0.d0) print*, 'Error: int1 + int2 < 0'
 elseif (abs(d1)  <  abs(d2)) then
    pint = pint*(int2 - int1)
    if (int2 - int1 < 0.d0) print*, 'Error: int2 - int1 < 0: ', int1, int2, '(', d1, d2,')'
 else
    pint = pint*(int1 - int2)
    if (int1 - int2 < 0.d0) print*, 'Error: int1 - int2 < 0: ', int1, int2, '(', d1, d2,')'
 endif

end function pint

real(doub_prec) function full_integral_3D(d, r0, R_0, h)

 real(doub_prec), intent(in) :: d, r0, R_0, h
 real(doub_prec) :: B1, B2, B3, a, logs, u, u2, h2
 real(doub_prec), parameter :: pi = 4.*atan(1.)
 real(doub_prec) :: tanphi, phi, a2, cosp, cosp2, mu2, mu2_1, r0h, r03, r0h2, r0h3, r0h_2, r0h_3, tanp
 real(doub_prec) :: r2, R_, linedist2, phi1, phi2, cosphi, sinphi
 real(doub_prec) :: I0, I1, I_1, I_2, I_3, I_4, I_5
 real(doub_prec) :: J_1, J_2, J_3, J_4, J_5
 real(doub_prec) :: D1, D2, D3

 r0h = r0/h
 tanphi = abs(d)/R_0
 phi = atan(tanphi)

 if (abs(r0h) < tiny(0.) .or. abs(R_0/h) < tiny(0.) .or. abs(phi) < tiny(0.)) then
    full_integral_3D = 0.0
    return
 endif

 h2 = h*h
 r03 = r0*r0*r0
 r0h2 = r0h*r0h
 r0h3 = r0h2*r0h
 r0h_2 = 1./r0h2
 r0h_3 = 1./r0h3

 if (r0 >= 2.0*h) then
    B3 = 0.25*h2*h
 elseif (r0 > h) then
    B3 = 0.25*r03 *(-4./3. + (r0h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2)
    B2 = 0.25*r03 *(-4./3. + (r0h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3)
 else
    B3 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2)
    B2 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2)
    B1 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3)
 endif

 a = R_0/r0
 a2 = a*a

 linedist2 = (r0*r0 + R_0*R_0)
 cosphi = cos(phi)
 R_ = R_0/cosphi
 r2 = (r0*r0 + R_*R_)

 D2 = 0.0
 D3 = 0.0

 if (linedist2 < h2) then
    !////// phi1 business /////
    cosp = R_0/sqrt(h2-r0*r0)
    call get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)

    D2 = -1./6.*I_2 + 0.25*(r0h) *I_3 - 0.15*r0h2 *I_4 + 1./30.*r0h3 *I_5 - 1./60. *r0h_3 *I1 + (B1-B2)/r03 *I0
 endif
 if (linedist2 < 4.*h2) then
    !////// phi2 business /////
    cosp = R_0/sqrt(4.0*h2-r0*r0)
    call get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)

    D3 = 1./3.*I_2 - 0.25*(r0h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2
 endif

 !//////////////////////////////
 call get_I_terms(cosphi,a2,a,I0,I1,I_2,I_3,I_4,I_5)

 if (r2 < h2) then
    full_integral_3D = r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0)
 elseif (r2 < 4.*h2) then
    full_integral_3D=  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - &
    &   1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2)
 else
    full_integral_3D = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3)
 endif

end function full_integral_3D

subroutine get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)
 real(doub_prec), intent(in) :: cosp,a2,a
 real(doub_prec), intent(out) :: I0,I1,I_2,I_3,I_4,I_5
 real(doub_prec) :: cosp2,p,tanp,u2,u,logs,I_1,mu2_1,fac

 p = acos(cosp)
 cosp2 = cosp*cosp
 mu2_1 = 1. / (1. + cosp2/a2)
 tanp = tan(p)

 I0  = p
 I_2 = p +    a2 * tanp
 I_4 = p + 2.*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2)

 u2 = (1.-cosp2)*mu2_1
 u = sqrt(u2)
 logs = log((1.+u)/(1.-u))
 I1 = atan2(u,a)

 fac = 1./(1.-u2)
 I_1 = 0.5*a*logs + I1
 I_3 = I_1 + a*0.25*(1.+a2)*(2.*u*fac + logs)
 I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10.*u - 6.*u*u2)*fac*fac + 3.*logs)

end subroutine get_I_terms

end module interpolations3D
