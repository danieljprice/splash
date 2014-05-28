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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 3D projections
!  (where rendered quantity is integrated along the line of sight)
!
!----------------------------------------------------------------------

module projections3D
 implicit none

 integer, parameter :: maxcoltable = 1000
 real, dimension(0:maxcoltable) :: coltable
 real, private :: dq2table = 4./maxcoltable
 real, private :: ddq2table = maxcoltable/4.

 public :: setup_integratedkernel
 public :: interpolate3D_projection
 public :: interpolate3D_proj_vec,interp3D_proj_vec_synctron
 public :: wfromtable

contains

subroutine setup_integratedkernel
!-------------------------------------------------------------
!     tabulates the integral through the cubic spline kernel
!     tabulated in (r/h)**2 so that sqrt is not necessary
!-------------------------------------------------------------
 use kernels, only:wfunc,radkernel2,cnormk3D
 implicit none
 integer :: i,j
 real :: rxy2,deltaz,dz,z,q2,wkern,coldens
 integer, parameter :: npts = 100

 !print "(1x,a)",'setting up integrated kernel table...'
 dq2table = radkernel2/maxcoltable
 ddq2table = 1./dq2table

 do i=0,maxcoltable-1
!
!--tabulate for (cylindrical) r**2 between 0 and radkernel**2
!
    rxy2 = i*dq2table
!
!--integrate z between 0 and sqrt(radkernel^2 - rxy^2)
!
    deltaz = sqrt(radkernel2 - rxy2)
    dz = deltaz/real(npts-1)
    coldens = 0.
    if (deltaz.ne.deltaz) print "(a)",'WARNING: NaN in kernel table setup'
    do j=1,npts
       z = (j-1)*dz
       q2 = rxy2 + z*z
       wkern = wfunc(q2)
       if (j.eq.1 .or. j.eq.npts) then
          coldens = coldens + 0.5*wkern*dz ! trapezoidal rule
       else
          coldens = coldens + wkern*dz
       endif
    enddo
    coltable(i)=2.0*coldens*cnormk3D
 end do
 coltable(maxcoltable) = 0.
 
 return
end subroutine setup_integratedkernel
!
! This function interpolates from the table of integrated kernel values
! to give w(q)
!
real function wfromtable(q2)
 implicit none
 real, intent(in) :: q2
 real :: dxx,dwdx
 integer :: index, index1
 !
 !--find nearest index in table
 !
 index = max(int(q2*ddq2table),0) ! the max prevents seg faults on NaNs for q2
 !index = min(index,maxcoltable) ! should be unnecessary if q2 < radkernel checked
 index1 = min(index + 1,maxcoltable)
 !
 !--find increment along from this index
 !
 dxx = q2 - index*dq2table
 !
 !--find gradient
 !
 dwdx = (coltable(index1) - coltable(index))*ddq2table
 !
 !--compute value of integrated kernel
 !
 wfromtable = coltable(index) + dwdx*dxx

end function wfromtable

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b weight_b dat_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     The (dimensionless) weight for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     the interface is written in this form to avoid floating exceptions
!     on physically scaled data.
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!     3D perspective added Nov 2005
!--------------------------------------------------------------------------

subroutine interpolate3D_projection(x,y,z,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise,zobserver,dscreen, &
     useaccelerate)

  use kernels, only:radkernel,radkernel2
  use timing,  only:wall_time,print_time
  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zobserver,dscreen
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm
  logical, intent(in) :: useaccelerate
  real :: row(npixx)

  integer :: ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,npixpartx,npixparty
  integer :: iprintinterval, iprintnext,ipixi,jpixi,jpixcopy
  integer :: nsubgrid,nfull,nok
#ifdef _OPENMP
  integer :: omp_get_num_threads,i
#else
  integer(kind=selected_int_kind(10)) :: iprogress,i  ! up to 10 digits
#endif
  real :: hi,hi1,hi21,radkern,wab,q2,xi,yi,xminpix,yminpix
  real :: term,termnorm,dy,dy2,ypix,zfrac,hsmooth,horigi
  real :: xpixmin,xpixmax,xmax,ypixmin,ypixmax,ymax
  real :: hmin,fac,hminall !,dhmin3
  real, dimension(npixx) :: xpix,dx2i
  real :: t_start,t_end,t_used
  logical :: iprintprogress,use3Dperspective,accelerate
  
  datsmooth = 0.
  term = 0.
  if (normalise) then
     print "(1x,a)",'projecting (normalised) from particles to pixels...'
     datnorm = 0.
  elseif (useaccelerate) then
     print "(1x,a)",'projecting (fast) from particles to pixels...'    
  else
     print "(1x,a)",'projecting from particles to pixels...'  
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0) then
     print "(1x,a)",'interpolate3D_proj: error: pixel width <= 0'
     return
  endif
  !nout = count(hh(1:npart).le.0.)
  !if (nout.gt.0) then
  !   print*,'interpolate3D_projection: warning: ignoring ',nout,' particles with h <= 0'
  !endif
  !
  !--check column density table has actually been setup
  !
  if (abs(coltable(1)).le.1.e-5) then
     call setup_integratedkernel
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
  use3Dperspective = abs(dscreen).gt.tiny(dscreen)
!
!--get starting CPU time
!
  call wall_time(t_start)

  xminpix = xmin - 0.5*pixwidthx
  yminpix = ymin - 0.5*pixwidthy
  xmax = xmin + npixx*pixwidthx
  ymax = ymin + npixy*pixwidthy
!
!--use a minimum smoothing length on the grid to make
!  sure that particles contribute to at least one pixel
!
  hmin = 0.5*max(pixwidthx,pixwidthy)
  !dhmin3 = 1./(hmin*hmin*hmin)
!
!--store x value for each pixel (for optimisation)
!  
  do ipix=1,npixx
     xpix(ipix) = xminpix + ipix*pixwidthx
  enddo
  nsubgrid = 0
  nok = 0
  hminall = huge(hminall)

!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,xmax,ymax,xminpix,yminpix,xpix,pixwidthx,pixwidthy) &
!$omp shared(npixx,npixy,dscreen,zobserver,use3dperspective,useaccelerate) &
!$omp shared(datnorm,normalise,radkernel,radkernel2) &
!$omp firstprivate(hmin) & !,dhmin3) &
!$omp private(hi,zfrac,xi,yi,radkern,xpixmin,xpixmax,ypixmin,ypixmax) &
!$omp private(hsmooth,horigi,hi1,hi21,term,termnorm,npixpartx,npixparty,jpixi,ipixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,accelerate) &
!$omp private(dx2i,row,q2,ypix,dy,dy2,wab) &
!$omp private(i,ipix,jpix,jpixcopy,fac) &
!$omp reduction(+:nsubgrid,nok) &
!$omp reduction(min:hminall)
!$omp master
#ifdef _OPENMP
  print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
#endif
!$omp end master

!$omp do schedule (guided, 2)
  over_particles: do i=1,npart
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
     if (itype(i).lt.0) cycle over_particles
     !
     !--set h related quantities
     !
     hi = hh(i)
     horigi = hi
     if (hi.le.0.) then
        cycle over_particles
     elseif (use3Dperspective) then
        if (z(i).gt.zobserver) cycle over_particles
        zfrac = abs(dscreen/(z(i)-zobserver))
        hi = hi*zfrac
     endif
     
     radkern = radkernel*hi ! radius of the smoothing kernel

     !--cycle as soon as we know the particle does not contribute
     xi = x(i)
     xpixmin = xi - radkern
     if (xpixmin.gt.xmax) cycle over_particles
     xpixmax = xi + radkern
     if (xpixmax.lt.xmin) cycle over_particles
     
     yi = y(i)
     ypixmin = yi - radkern
     if (ypixmin.gt.ymax) cycle over_particles
     ypixmax = yi + radkern
     if (ypixmax.lt.ymin) cycle over_particles

     !--take resolution length as max of h and 1/2 pixel width
     if (hi.lt.hmin) then
        hminall = min(hi,hminall)
        nsubgrid = nsubgrid + 1
        hsmooth = hmin
        fac = 1. !(horigi*horigi*horigi)*dhmin3      ! factor by which to adjust the weight
     else
        fac = 1.
        hsmooth = hi
        nok = nok + 1
     endif
     radkern = radkernel*hsmooth
     !
     !--set kernel related quantities
     !
     hi1 = 1./hsmooth
     hi21 = hi1*hi1
     termnorm = weight(i)*fac*horigi
     term = termnorm*dat(i) ! h gives the z length scale (NB: no perspective)
     !
     !--for each particle work out which pixels it contributes to
     !
     npixpartx = int(radkern/pixwidthx) + 1
     npixparty = int(radkern/pixwidthy) + 1
     jpixi = int((yi-ymin)/pixwidthy) + 1
     ipixi = int((xi-xmin)/pixwidthx) + 1
     ipixmin = ipixi - npixpartx
     ipixmax = ipixi + npixpartx
     jpixmin = jpixi - npixparty
     jpixmax = jpixi + npixparty

!     ipixmin = int((xi - radkern - xmin)/pixwidth)
!     jpixmin = int((yi - radkern - ymin)/pixwidth)
!     ipixmax = ipixmin + npixpart !!int((xi + radkern - xmin)/pixwidth) + 1
!     jpixmax = jpixmin + npixpart !!int((yi + radkern - ymin)/pixwidth) + 1
     !
     !--loop over pixels, adding the contribution from this particle
     !  copy by quarters if all pixels within domain
     !
     accelerate = useaccelerate .and. (.not.normalise) &
                 .and. npixpartx.gt.5 .and. npixparty.gt.5 &
                 .and. ipixmin.ge.1 .and. ipixmax.le.npixx &
                 .and. jpixmin.ge.1 .and. jpixmax.le.npixy
     
     if (accelerate) then
        !--adjust xi, yi to centre of pixel
        xi = xminpix + ipixi*pixwidthx
        yi = yminpix + jpixi*pixwidthy
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
        enddo
        do jpix = jpixi,jpixmax
           ypix = yminpix + jpix*pixwidthy
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixi,ipixmax
              q2 = dx2i(ipix) + dy2
              !
              !--SPH kernel - integral through cubic spline
              !  interpolate from a pre-calculated table
              !
              if (q2.lt.radkernel2) then
                 wab = wfromtable(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !                  
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                 row(ipix) = term*wab
              else
                 row(ipix) = 0.
              endif
           enddo
           !--NB: the following actions can and should be vectorized (but I don't know how...)
           !--copy top right -> top left
           do ipix=ipixmin,ipixi-1
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + row(ipixmax-(ipix-ipixmin))
           enddo
           if (jpix.ne.jpixi) then
              jpixcopy = jpixi - (jpix-jpixi)
              !--copy top right -> bottom left 
              do ipix=ipixmin,ipixi-1
                 datsmooth(ipix,jpixcopy) = datsmooth(ipix,jpixcopy) + row(ipixmax-(ipix-ipixmin))
              enddo
              !--copy top right -> bottom right
              do ipix=ipixi,ipixmax
                 datsmooth(ipix,jpixcopy) = datsmooth(ipix,jpixcopy) + row(ipix)
              enddo
           endif
        enddo
          
     else
        ipixmin = int((xi - radkern - xmin)/pixwidthx)
        ipixmax = int((xi + radkern - xmin)/pixwidthx)
        jpixmin = int((yi - radkern - ymin)/pixwidthy)
        jpixmax = int((yi + radkern - ymin)/pixwidthy)

        if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
        if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
        enddo

        do jpix = jpixmin,jpixmax
           ypix = yminpix + jpix*pixwidthy
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixmin,ipixmax
              !xpix = xminpix + ipix*pixwidthx
              !dx = xpix - xi
              !rab2 = (xminpix + ipix*pixwidthx - xi)**2 + dy2
              q2 = dx2i(ipix) + dy2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
              !
              !--SPH kernel - integral through cubic spline
              !  interpolate from a pre-calculated table
              !
              if (q2.lt.radkernel2) then
                 wab = wfromtable(q2)
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
     endif

  enddo over_particles
!$omp end do
!$omp end parallel
!
!--normalise dat array
!
  if (normalise) then
     !--normalise everywhere (required if not using SPH weighting)
     where (datnorm > tiny(datnorm))
        datsmooth = datsmooth/datnorm
     end where
  endif
!
!--warn about subgrid interpolation
!
  if (nsubgrid.gt.1) then
     nfull = int((xmax-xmin)/(hminall)) + 1
     if (nsubgrid.gt.0.1*nok) &
     print "(a,i9,a,/,a,i6,a)",' Warning: pixel size > 2h for ',nsubgrid,' particles', &
                               '          need',nfull,' pixels for full resolution'
  endif
!
!--get/print timings
!
  call wall_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.10.) call print_time(t_used)
  
  return

end subroutine interpolate3D_projection

!--------------------------------------------------------------------------
!
!     Same as previous but for a vector quantity
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price 23/12/04
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_vec(x,y,z,hh,weight,vecx,vecy,itype,npart,&
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,normalise,zobserver,dscreen)

  use kernels, only:radkernel,radkernel2
  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zobserver,dscreen
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy
  logical, intent(in) :: normalise
  real, dimension(:,:), allocatable :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,ierr
  real :: hi,hi1,hi21,radkern,q2,wab,rab2,const,zfrac,hsmooth
  real :: termx,termy,termnorm,dx,dy,dy2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  if (normalise) then
     print "(1x,a)",'projecting vector (normalised) from particles to pixels...'
     allocate(datnorm(npixx,npixy),stat=ierr)
     if (ierr /= 0) then
        print "(a)",'interpolate3D_proj_vec: error allocating memory'
        return
     endif
     datnorm = 0.
  else
     print "(1x,a)",'projecting vector from particles to pixels...'  
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print "(a)",'interpolate3D_proj_vec: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,vecx,vecy,itype,vecsmoothx,vecsmoothy,npart) &
!$omp shared(xmin,ymin,pixwidthx,pixwidthy,zobserver,dscreen,datnorm) &
!$omp shared(npixx,npixy,normalise,radkernel,radkernel2) &
!$omp private(hi,radkern,const,zfrac,ypix,xpix) &
!$omp private(hsmooth,hi1,hi21,termx,termy,termnorm) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax) &
!$omp private(dy,dy2,dx,rab2,q2,wab) &
!$omp private(i,ipix,jpix)

!$omp do schedule(guided, 2)
  over_particles: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_particles
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     const = weight(i)*hi ! h gives the z length scale (NB: no perspective)
     if (hi.le.0.) then
        cycle over_particles
     elseif (abs(dscreen).gt.tiny(dscreen)) then
        if (z(i).gt.zobserver) cycle over_particles
        zfrac = abs(dscreen/(z(i)-zobserver))
        hi = hi*zfrac
     endif

     !--take resolution length as max of h and 1/2 pixel width
     hsmooth = max(hi,0.5*min(pixwidthx,pixwidthy))

     radkern = radkernel*hsmooth    ! radius of the smoothing kernel
     hi1 = 1./hsmooth
     hi21 = hi1*hi1
        
     termx = const*vecx(i)
     termy = const*vecy(i)
     termnorm = const
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

     ! PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     ! PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

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
        dy2 = dy*dy
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)
           rab2 = dx**2 + dy2
           q2 = rab2*hi21
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !  
              vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab
              vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab
              if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
           endif

        enddo
     enddo

  enddo over_particles
!$omp end do
!$omp end parallel

  if (normalise .and. allocated(datnorm)) then
     !--normalise everywhere
     where (datnorm > tiny(datnorm))
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif
  if (allocated(datnorm)) deallocate(datnorm)

  return

end subroutine interpolate3D_proj_vec

!--------------------------------------------------------------------------
!
!     Computes synchrotron emission (Stokes Q, U) for a given B field
!     at present assuming no faraday rotation
!
!     For references see:
!
!      Urbanik et al. (1997), A&A 326, 465
!      Sokoloff et al. (1998), MNRAS, 299, 189
!      Gomez & Cox (2004), ApJ 615, 744 (for Cosmic Ray distribution esp.)
!
!     Faraday rotation could be included easily but
!     I have not yet done so.
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field   : stokesQ (npixx,npixy)
!                                     : stokesU (npixx,npixy)
!                                     : stokesI (npixx,npixy)
!
!     DOES NOT WORK FOR ROTATED CONFIGURATIONS YET!!
!     (ie. z is assumed to be z_galaxy in the cosmic ray distribution)
!
!     Daniel Price 14/03/07
!--------------------------------------------------------------------------

subroutine interp3D_proj_vec_synctron(x,y,z,hh,weight,vecx,vecy,itype,npart,&
     xmin,ymin,stokesQ,stokesU,stokesI,npixx,npixy,pixwidth,rcrit,zcrit,alpha, &
     qpixwidth,getIonly,utherm,uthermcutoff)

  use kernels, only:radkernel,radkernel2
  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidth,rcrit,zcrit,alpha,qpixwidth
  logical, intent(in) :: getIonly
  real, intent(out), dimension(npixx,npixy) :: stokesQ,stokesU,stokesI
  real, intent(in), dimension(npart), optional :: utherm
  real, intent(in), optional :: uthermcutoff

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,hi21,radkern,q2,wab,const,hsmooth
  real :: termx,termy,term,dy,dy2,ypix,xi,yi,zi
  real :: crdens,emissivity,Bperp,angle,pintrinsic,rcyl
  real, dimension(npixx) :: dx2i
  
  if (getIonly) then
     stokesI = 0.
  else
     stokesU = 0.
     stokesQ = 0.
  endif
  termx = 0.
  termy = 0.
  term = 0.
  pintrinsic = (3. + 3.*alpha)/(5. + 3.*alpha)

  if (getIonly) then
     print "(1x,a)",'getting synchrotron intensity map from B field...'
  else
     print "(1x,a)",'getting synchrotron polarisation map from B field...'
     print*,' assuming cosmic ray electron distribution exp(-r/',rcrit,' -z/',zcrit,') (kpc)'
     print*,' synchrotron spectral index I_nu = nu^-',alpha
     print*,' intrinsic polarisation fraction = ',pintrinsic
  endif
  if (present(utherm) .and. present(uthermcutoff)) then
     print*,' using only particles with utherm > ',uthermcutoff
  endif
  
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj_vec_synchrotron: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,vecx,vecy,itype,stokesq,stokesu,stokesi,npart) &
!$omp shared(xmin,ymin,pixwidth,rcrit,zcrit,alpha,radkernel,radkernel2) &
!$omp shared(npixx,npixy,pintrinsic,qpixwidth,getionly,utherm,uthermcutoff) &
!$omp private(hi,xi,yi,zi,radkern,const) &
!$omp private(hsmooth,hi1,hi21,term,termx,termy) &
!$omp private(rcyl,crdens,bperp,emissivity,angle) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax) &
!$omp private(dy,dy2,dx2i,ypix,q2,wab) &
!$omp private(i,ipix,jpix)

!$omp do schedule(guided, 2)
  over_particles: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_particles
     !
     !--skip particles with utherm < uthermcutoff
     !
     if (present(utherm) .and. present(uthermcutoff)) then
        if (utherm(i).lt.uthermcutoff) cycle over_particles
     endif
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_particles
     const = weight(i)*hi ! h gives the z length scale (NB: no perspective)
     zi = z(i)

     !--take resolution length as max of h and 1/2 pixel width
     !  (for intensity calculation, qpixwidth is pixel width of Q,U calculation)
     hsmooth = max(hi,0.5*pixwidth,0.5*qpixwidth)
     hi1 = 1./hsmooth
     hi21 = hi1*hi1
     radkern = radkernel*hsmooth    ! radius of the smoothing kernel
     xi = x(i)
     yi = y(i)
     
     !--assumed distribution of cosmic ray electrons in galaxy
     !  (should use UNROTATED x,y if rotation added)
     rcyl = sqrt(xi**2 + yi**2)
     crdens = exp(-rcyl/rcrit - abs(zi)/zcrit)
     
     !--calculate synchrotron emissivity based on Bperp and a spectral index alpha
     Bperp = sqrt(vecx(i)**2 + vecy(i)**2)
     emissivity = crdens*Bperp**(1. + alpha)
     
     if (getIonly) then
        term = emissivity*const
        termx = 0.
        termy = 0.
     else
        term = 0.
        !--faraday rotation would change angle here
        angle = atan2(vecy(i),vecx(i))     
        termx = pintrinsic*emissivity*const*COS(angle)
        termy = pintrinsic*emissivity*const*SIN(angle)
     endif
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((xi - radkern - xmin)/pixwidth)
     jpixmin = int((yi - radkern - ymin)/pixwidth)
     ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
     jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

     ! PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     ! PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy

     !
     !--precalculate an array of dx2 for this particle (optimisation)
     !
     do ipix=ipixmin,ipixmax
        dx2i(ipix) = ((xmin + (ipix-0.5)*pixwidth - xi)**2)*hi21
     enddo
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix-0.5)*pixwidth
        dy = ypix - yi
        dy2 = dy*dy*hi21
        do ipix = ipixmin,ipixmax
           !xpix = xmin + (ipix-0.5)*pixwidth
           !dx = xpix - xi
           q2 = dx2i(ipix) + dy2
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !
              if (getIonly) then
                 stokesI(ipix,jpix) = stokesI(ipix,jpix) + term*wab           
              else
                 stokesQ(ipix,jpix) = stokesQ(ipix,jpix) + termx*wab
                 stokesU(ipix,jpix) = stokesU(ipix,jpix) + termy*wab
              endif
           endif

        enddo
     enddo

  enddo over_particles
!$omp end do
!$omp end parallel

  return

end subroutine interp3D_proj_vec_synctron


end module projections3D
