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
 use kernels,       only:radkernel2,radkernel,cnormk3D,wallint
 use interpolation, only:iroll
 use timing,        only:wall_time,print_time
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
 real(doub_prec), allocatable :: datnorm(:,:,:)

 integer :: i,ipix,jpix,kpix
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nwarn,threadid
 real :: xminpix,yminpix,zminpix,hmin !,dhmin3
 real, dimension(npixx) :: dx2i
 real :: xi,yi,zi,hi,hi1,hi21,radkern,wab,q2,const,dyz2,dz2
 real :: term,termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax,dfac
 real :: t_start,t_end,t_used
 logical :: iprintprogress

! Exact rendering
 real :: pixint, wint
 logical, parameter :: exact_rendering = .true.   ! use exact rendering y/n
 integer :: usedpart, negflag

!$ integer :: omp_get_num_threads,omp_get_thread_num
 integer(kind=selected_int_kind(10)) :: iprogress,j  ! up to 10 digits

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
 !if (any(hh(1:npart) <= tiny(hh))) then
!    print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
 !endif

 call wall_time(t_start)

 datsmooth = 0.
 if (normalise) then
    allocate(datnorm(npixx,npixy,npixz))
    datnorm = 0.
 endif
 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000) !.or. exact_rendering
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
 j = 0_8
 threadid = 1
 !
 !--loop over particles
 !
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
!$omp shared(xminpix,yminpix,zminpix,pixwidthx,pixwidthy,pixwidthz) &
!$omp shared(npixx,npixy,npixz,const) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
!$omp shared(hmin,pixwidthmax) &
!$omp shared(iprintprogress,iprintinterval,j) &
!$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
!$omp private(term,termnorm,xpixi,iprogress) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
!$omp private(pixint,wint,negflag,dfac,threadid) &
!$omp firstprivate(iprintnext) &
!$omp reduction(+:nwarn,usedpart)
!$omp master
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end master

!$omp do schedule (guided, 2)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
    if (iprintprogress) then
       !$omp atomic
       j=j+1_8
       !$ threadid = omp_get_thread_num()
       iprogress = 100*j/npart
       if (iprogress >= iprintnext .and. threadid==1) then
          write(*,"(i3,'%.')",advance='no') iprogress
          iprintnext = iprintnext + iprintinterval
       endif
    endif
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle over_parts

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
    !dfac = hi**3/(pixwidthx*pixwidthy*const)
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
 if (allocated(datnorm)) deallocate(datnorm)

 call wall_time(t_end)
 t_used = t_end - t_start
 if (t_used > 10.) call print_time(t_used)

 !print*, 'Number of particles in the volume: ', usedpart

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
 !$ integer :: omp_get_num_threads
 integer(kind=selected_int_kind(10)) :: iprogress  ! up to 10 digits

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
!$omp shared(iprintprogress,iprintinterval) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
!$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
!$omp private(term,termnorm,xpixi) &
!$omp private(iprogress,iprintnext) &
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

end module interpolations3D
