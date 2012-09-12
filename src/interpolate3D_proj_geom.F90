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
!  Module containing routines required for 3D projections
!  in different coordinate systems
!
!----------------------------------------------------------------------

module projections3Dgeom
 use projections3D, only:setup_integratedkernel,wfromtable,coltable
 use kernels,       only:radkernel,radkernel2
 implicit none

 public :: interpolate3D_proj_geom
! public :: interpolate3D_proj_geom_vec

#ifdef _OPENMP
 character(len=5), parameter :: str = 'cpu s'
#else
 character(len=1), parameter :: str = 's'
#endif

contains

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
!     ** This version is for 3D projections in alternative coordinate
!     ** systems, e.g. \Int rho d\phi
!
!     The (dimensionless) weight for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     the interface is written in this form to avoid floating exceptions
!     on physically scaled data.
!
!     Input: particle coordinates      : x,y,z (npart)
!            smoothing lengths         : hh    (npart)
!            weight for each particle  : weight (npart)
!            scalar data to smooth     : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price July 2011
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_geom(x,y,z,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise,igeom,&
     iplotx,iploty,iplotz,ix)

  use geometry, only:igeom_cartesian,coord_transform,coord_is_length
  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real,    intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real,    intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real,    intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  integer, intent(in) :: igeom,iplotx,iploty,iplotz
  integer, dimension(3), intent(in) :: ix
  real, dimension(npixx,npixy) :: datnorm

  integer :: ipix,jpix,ixcoord,iycoord,izcoord
  integer :: iprintinterval, iprintnext, itmin
#ifdef _OPENMP
  integer :: omp_get_num_threads,i
#else
  integer(kind=selected_int_kind(10)) :: iprogress,i  ! up to 10 digits
#endif
  real, dimension(3) :: xcoord, xpix
  real :: hi,hi1,hi21,radkern,wab,q2,xci(3),xminpix,yminpix
  real :: term,termnorm,dx,dx2,dy,dy2
  real :: xmax,ymax
  real :: t_start,t_end,t_used,tsec
  logical :: iprintprogress,islengthx,islengthy,islengthz
  
  datsmooth = 0.
  term = 0.
  if (normalise) then
     print "(1x,a)",'projecting (normalised, non-cartesian) from particles to pixels...'
     datnorm = 0.
  else
     print "(1x,a)",'projecting (non-cartesian) from particles to pixels...'  
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0) then
     print "(1x,a)",'interpolate3D_proj_geom: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.0.)) then
     print*,'interpolate3D_proj_geom: warning: ignoring some or all particles with h <= 0'
  endif
  !
  !--get plotted coordinates in range 1->ndim
  !
  ixcoord = iplotx - ix(1) + 1
  iycoord = iploty - ix(1) + 1
  izcoord = iplotz - ix(1) + 1
  if (ixcoord.le.0 .or. ixcoord.gt.3) then
     print*,' ERROR finding x coordinate offset, cannot render'
     return
  endif
  if (iycoord.le.0 .or. iycoord.gt.3) then
     print*,' ERROR finding y coordinate offset, cannot render'
     return
  endif
  if (izcoord.le.0 .or. izcoord.gt.3) then
     print*,' ERROR finding y coordinate offset, cannot render'
     return
  endif
  !
  !--check if coordinate is a length (i.e., not an angle)
  !
  islengthx = coord_is_length(ixcoord,igeom)
  islengthy = coord_is_length(iycoord,igeom)
  islengthz = coord_is_length(izcoord,igeom)
  !print*,' islength = ',islengthx,islengthy,islengthz
  !print*,' y axis is ',iycoord
  !print*,' x axis is ',ixcoord
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
!
!--get starting CPU time
!
  call cpu_time(t_start)

  xminpix = xmin - 0.5*pixwidthx
  yminpix = ymin - 0.5*pixwidthy
  xmax = xmin + npixx*pixwidthx
  ymax = ymin + npixy*pixwidthy
!
!--use a minimum smoothing length on the grid to make
!  sure that particles contribute to at least one pixel
!
!  hmin = 0.5*max(pixwidthx,pixwidthy)

!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,xmax,ymax,xminpix,yminpix,pixwidthx,pixwidthy) &
!$omp shared(npixx,npixy,ixcoord,iycoord,izcoord,islengthx,islengthy,islengthz,igeom) &
!$omp shared(datnorm,normalise) &
!$omp private(hi,xi,xci,radkern) &
!$omp private(hi1,hi21,term,termnorm) &
!$omp private(q2,dx,dx2,dy,dy2,wab,xcoord,xpix) &
!$omp private(i,ipix,jpix)
!$omp master
!$    print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
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
     if (hi.le.0.) cycle over_particles

     radkern = radkernel*hi ! radius of the smoothing kernel

     xci(1) = sqrt(x(i)**2 + z(i)**2)
     if (xci(1).lt.0.) print*,'error r < 0',xci(1)
     xci(2) = y(i)
     xci(3) = z(i)
     !--get particle coords in new coord system
     !call coord_transform(xci(:),3,igeom_cartesian,xi(:),3,igeom)
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     hi21 = hi1*hi1

     ! h gives the z length scale (NB: no perspective)
     if (islengthz) then
        termnorm = weight(i)*hi
     else
        termnorm = weight(i)
     endif
     term = termnorm*dat(i)

     !
     !--loop over pixels, adding the contribution from this particle
     !
     if (islengthz) then
        xcoord(izcoord) = 1. !xci(3)     
     else
        xcoord(izcoord) = 0. ! use phi=0 so get x = r cos(phi) = r
     endif
     do jpix = 1,npixy
        xcoord(iycoord) = yminpix + jpix*pixwidthy
        do ipix = 1,npixx
           xcoord(ixcoord) = xminpix + ipix*pixwidthx
           xpix(:) = xcoord(:)
           !call coord_transform(xcoord(:),3,igeom,xpix(:),3,igeom_cartesian)
           
!           !--this is in cartesians
           dy   = xpix(iycoord) - xci(2)
           dx   = xpix(ixcoord) - xci(1)

           dx2  = dx*dx
           dy2  = dy*dy
           q2   = (dx2 + dy2)*hi21
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
!--get ending CPU time
!
  call cpu_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.60.) then
     itmin = int(t_used/60.)
     tsec = t_used - (itmin*60.)
     print "(1x,a,i4,a,f5.2,1x,a)",'completed in',itmin,' min ',tsec,trim(str)
  else
     print "(1x,a,f5.2,1x,a)",'completed in ',t_used,trim(str)
  endif
  
  return

end subroutine interpolate3D_proj_geom

end module projections3Dgeom
