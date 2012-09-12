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
!  Module containing all of the routines required for interpolation
!  from 3D data to a 3D grid (SLOW!)
!
!----------------------------------------------------------------------

module interpolations3D
 use kernels, only:radkernel2,radkernel,cnormk3D
 implicit none
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
!--------------------------------------------------------------------------

subroutine interpolate3D(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth,zpixwidth,&
     normalise,periodicx,periodicy,periodicz)
  implicit none
  integer, intent(in) :: npart,npixx,npixy,npixz
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,zmin,pixwidth,zpixwidth
  real, intent(out), dimension(npixx,npixy,npixz) :: datsmooth
  logical, intent(in) :: normalise,periodicx,periodicy,periodicz
  real, dimension(npixx,npixy,npixz) :: datnorm

  integer :: i,ipix,jpix,kpix
  integer :: iprintinterval,iprintnext
  integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
  integer :: ipixi,jpixi,kpixi,nxpix,nwarn
  real :: xminpix,yminpix,zminpix,hmin !,dhmin3
  real, dimension(npixx) :: dx2i
  real :: xi,yi,zi,hi,hi1,hi21,radkern,wab,q2,const,dyz2,dz2
  real :: term,termnorm,dy,dz,ypix,zpix,xpixi
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
     print "(1x,a)",'interpolating from particles to 3D grid (normalised) ...'  
  else
     print "(1x,a)",'interpolating from particles to 3D grid (non-normalised) ...'
  endif
  if (pixwidth.le.0.) then
     print "(1x,a)",'interpolate3D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
  endif

  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval
  !
  !--get starting CPU time
  !
  !call cpu_time(t_start)

  xminpix = xmin - 0.5*pixwidth
  yminpix = ymin - 0.5*pixwidth
  zminpix = zmin - 0.5*zpixwidth
  !
  !--use a minimum smoothing length on the grid to make
  !  sure that particles contribute to at least one pixel
  !
  hmin = 0.5*pixwidth
  !dhmin3 = 1./(hmin*hmin*hmin)

  const = cnormk3D  ! normalisation constant (3D)
  nwarn = 0
  !
  !--loop over particles
  !
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
!$omp shared(xminpix,yminpix,zminpix,pixwidth,zpixwidth) &
!$omp shared(npixx,npixy,npixz,const) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
!$omp shared(hmin) & !,dhmin3) &
!$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
!$omp private(term,termnorm,xpixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
!$omp reduction(+:nwarn)
!$omp master
#ifdef _OPENMP
  print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
#endif
!$omp end master

!$omp do schedule (guided, 2)
  over_parts: do i=1,npart
     !
     !--report on progress
     !
#ifndef _OPENMP
     if (iprintprogress) then
        iprogress = 100*i/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
#endif
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts

     hi = hh(i)
     if (hi.le.0.) then
        cycle over_parts
     elseif (hi.lt.hmin) then
     !
     !--use minimum h to capture subgrid particles
     !  (get better results *without* adjusting weights)
     !
        termnorm = const*weight(i) !*(hi*hi*hi)*dhmin3
        hi = hmin
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
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((xi - radkern - xmin)/pixwidth)
     jpixmin = int((yi - radkern - ymin)/pixwidth)
     kpixmin = int((zi - radkern - zmin)/zpixwidth)
     ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
     jpixmax = int((yi + radkern - ymin)/pixwidth) + 1
     kpixmax = int((zi + radkern - zmin)/zpixwidth) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1      ! make sure they only contribute
        if (ipixmax.gt.npixx) ipixmax = npixx  ! to pixels in the image
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif
     if (.not.periodicz) then
        if (kpixmin.lt.1)     kpixmin = 1
        if (kpixmax.gt.npixz) kpixmax = npixz
     endif
     !
     !--precalculate an array of dx2 for this particle (optimisation)
     !
     nxpix = 0
     do ipix=ipixmin,ipixmax
        nxpix = nxpix + 1
        ipixi = ipix
        if (periodicx) then
           if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
           if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
        endif
        xpixi = xminpix + ipix*pixwidth
        !--watch out for errors with perioic wrapping...
        if (nxpix.le.size(dx2i)) then
           dx2i(nxpix) = ((xpixi - xi)**2)*hi21
        endif
     enddo
     
     !--if particle contributes to more than npixx pixels
     !  (i.e. periodic boundaries wrap more than once)
     !  truncate the contribution and give warning
     if (nxpix.gt.npixx) then
        nwarn = nwarn + 1
        ipixmax = ipixmin + npixx - 1
     endif
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do kpix = kpixmin,kpixmax
        kpixi = kpix
        if (periodicz) then
           if (kpixi.lt.1)     kpixi = mod(kpixi,npixz) + npixz
           if (kpixi.gt.npixz) kpixi = mod(kpixi-1,npixz) + 1
        endif
        zpix = zminpix + kpix*zpixwidth
        dz = zpix - zi
        dz2 = dz*dz*hi21
        
        do jpix = jpixmin,jpixmax
           jpixi = jpix
           if (periodicy) then
              if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
              if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
           endif
           ypix = yminpix + jpix*pixwidth
           dy = ypix - yi
           dyz2 = dy*dy*hi21 + dz2
           
           nxpix = 0
           do ipix = ipixmin,ipixmax
              nxpix = nxpix + 1
              ipixi = ipix
              if (periodicx) then
                 if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
                 if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
              endif
              q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.radkernel2) then                  
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
           enddo
        enddo
     enddo
  enddo over_parts
!$omp end do
!$omp end parallel

  if (nwarn.gt.0) then
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
  
  return

end subroutine interpolate3D

subroutine interpolate3D_vec(x,y,z,hh,weight,datvec,itype,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth,zpixwidth,&
     normalise,periodicx,periodicy,periodicz)
  implicit none
  integer, intent(in) :: npart,npixx,npixy,npixz
  real, intent(in), dimension(npart)    :: x,y,z,hh,weight
  real, intent(in), dimension(npart,3)  :: datvec
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,zmin,pixwidth,zpixwidth
  real, intent(out), dimension(3,npixx,npixy,npixz) :: datsmooth
  logical, intent(in) :: normalise,periodicx,periodicy,periodicz
  real, dimension(npixx,npixy,npixz) :: datnorm

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
     print "(1x,a)",'interpolating from particles to 3D grid (normalised) ...'  
  else
     print "(1x,a)",'interpolating from particles to 3D grid (non-normalised) ...'
  endif
  if (pixwidth.le.0.) then
     print "(1x,a)",'interpolate3D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
  endif

  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval
  !
  !--get starting CPU time
  !
  !call cpu_time(t_start)

  xminpix = xmin - 0.5*pixwidth
  yminpix = ymin - 0.5*pixwidth
  zminpix = zmin - 0.5*zpixwidth
!  xmax = xmin + npixx*pixwidth
!  ymax = ymin + npixy*pixwidth

  const = cnormk3D  ! normalisation constant (3D)
  nwarn = 0

!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,datvec,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
!$omp shared(xminpix,yminpix,zminpix,pixwidth,zpixwidth) &
!$omp shared(npixx,npixy,npixz,const) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
!$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
!$omp private(term,termnorm,xpixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
!$omp reduction(+:nwarn)
!$omp master
#ifdef _OPENMP
  print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
#endif
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
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
#endif
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts

     hi = hh(i)
     if (hi.le.0.) cycle over_parts

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
     ipixmin = int((xi - radkern - xmin)/pixwidth)
     jpixmin = int((yi - radkern - ymin)/pixwidth)
     kpixmin = int((zi - radkern - zmin)/zpixwidth)
     ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
     jpixmax = int((yi + radkern - ymin)/pixwidth) + 1
     kpixmax = int((zi + radkern - zmin)/zpixwidth) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1      ! make sure they only contribute
        if (ipixmax.gt.npixx) ipixmax = npixx  ! to pixels in the image
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif
     if (.not.periodicz) then
        if (kpixmin.lt.1)     kpixmin = 1
        if (kpixmax.gt.npixz) kpixmax = npixz
     endif
     !
     !--precalculate an array of dx2 for this particle (optimisation)
     !
     nxpix = 0
     do ipix=ipixmin,ipixmax
        nxpix = nxpix + 1
        ipixi = ipix
        if (periodicx) then
           if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
           if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
        endif
        xpixi = xminpix + ipix*pixwidth
        !--watch out for errors with perioic wrapping...
        if (nxpix.le.size(dx2i)) then
           dx2i(nxpix) = ((xpixi - xi)**2)*hi21
        endif
     enddo
     
     !--if particle contributes to more than npixx pixels
     !  (i.e. periodic boundaries wrap more than once)
     !  truncate the contribution and give warning
     if (nxpix.gt.npixx) then
        nwarn = nwarn + 1
        ipixmax = ipixmin + npixx - 1
     endif
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do kpix = kpixmin,kpixmax
        kpixi = kpix
        if (periodicz) then
           if (kpixi.lt.1)     kpixi = mod(kpixi,npixz) + npixz
           if (kpixi.gt.npixz) kpixi = mod(kpixi-1,npixz) + 1
        endif
        zpix = zminpix + kpix*zpixwidth
        dz = zpix - zi
        dz2 = dz*dz*hi21
        
        do jpix = jpixmin,jpixmax
           jpixi = jpix
           if (periodicy) then
              if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
              if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
           endif
           ypix = yminpix + jpix*pixwidth
           dy = ypix - yi
           dyz2 = dy*dy*hi21 + dz2
           
           nxpix = 0
           do ipix = ipixmin,ipixmax
              ipixi = ipix
              if (periodicx) then
                 if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
                 if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
              endif
              nxpix = nxpix + 1
              q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.radkernel2) then                  
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
!$omp end do
!$omp end parallel
  
  if (nwarn.gt.0) then
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
              if (datnorm(ipix,jpix,kpix).gt.tiny(datnorm)) then
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
 implicit none
 real, intent(in) :: q2

 wkernel = wfunc(q2)

end function wkernel

end module interpolations3D
