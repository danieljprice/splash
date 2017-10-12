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
 use geometry,      only:igeom_cartesian,coord_transform,coord_is_length, &
                         coord_transform_limits,igeom_cylindrical,get_coord_limits
 implicit none

 public :: interpolate3D_proj_geom, interpolate3D_xsec_geom
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

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real,    intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real,    intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real,    intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(inout) :: normalise
  integer, intent(in) :: igeom,iplotx,iploty,iplotz
  integer, dimension(3), intent(in) :: ix
  real, dimension(npixx,npixy) :: datnorm
  real, parameter :: pi = 3.1415926536

  integer :: ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ixcoord,iycoord,izcoord,ierr
  integer :: iprintinterval, iprintnext, itmin
#ifdef _OPENMP
  integer :: omp_get_num_threads,i
#else
  integer(kind=selected_int_kind(10)) :: iprogress,i  ! up to 10 digits
#endif
  real, dimension(3) :: xcoord, xpix
  real, dimension(3), save :: xci, xi
!$omp threadprivate(xci,xi)
  real :: hi,hi1,hi21,radkern,wab,q2,xminpix,yminpix
  real :: term,termnorm,dx,dx2,dy,dy2
  real :: xmax,ymax,hmin,horigi
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
  !--get information about the coordinates
  !
  call get_coord_info(iplotx,iploty,iplotz,ix(1),igeom,ixcoord,iycoord,izcoord,&
                      islengthx,islengthy,islengthz,ierr)
  if (ierr /= 0) return
  !
  !--if z coordinate is not a length, use normalised interpolation
  !  (e.g. azimuthally averaged density)
  !
  if (.not.islengthz) normalise = .true.
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
  hmin = 0.
  if (islengthx) hmin = 0.5*pixwidthx
  if (islengthy) hmin = max(hmin,0.5*pixwidthy)

!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart,hmin) &
!$omp shared(xmin,ymin,xmax,ymax,xminpix,yminpix,pixwidthx,pixwidthy) &
!$omp shared(npixx,npixy,ixcoord,iycoord,izcoord,islengthx,islengthy,islengthz,igeom) &
!$omp shared(datnorm,normalise,radkernel,radkernel2) &
!$omp private(hi,horigi,radkern) &
!$omp private(hi1,hi21,term,termnorm) &
!$omp private(q2,dx,dx2,dy,dy2,wab,xcoord,xpix) &
!$omp private(i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,ierr)
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
     horigi = hh(i)
     if (horigi.le.0.) cycle over_particles
     hi = max(horigi,hmin)

     radkern = radkernel*hi ! radius of the smoothing kernel
     
     !
     !--get limits of contribution from particle in cartesian space
     !
     xci(1) = x(i)
     xci(2) = y(i)
     xci(3) = z(i)
     call get_pixel_limits(xci,xi,radkern,ipixmin,ipixmax,jpixmin,jpixmax,igeom,&
                           npixx,npixy,pixwidthx,pixwidthy,xmin,ymin,ixcoord,iycoord,ierr)
     if (ierr /= 0) cycle over_particles  
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     hi21 = hi1*hi1

     ! h gives the z length scale (NB: no perspective)
     if (islengthz) then
        termnorm = weight(i)*horigi
     elseif (igeom.eq.igeom_cylindrical) then
        termnorm = weight(i)*atan(radkern/xi(ixcoord))/pi
     else
        termnorm = weight(i)
     endif
     term = termnorm*dat(i)

     !
     !--loop over pixels, adding the contribution from this particle
     !
     if (islengthz) then
        xcoord(izcoord) = xi(izcoord) ! assume all pixels at same r as particlefor theta-phi
     else
        xcoord(izcoord) = 0. ! use phi=0 so get x = r cos(phi) = r
     endif
     do jpix = jpixmin,jpixmax
        xcoord(iycoord) = yminpix + jpix*pixwidthy
        do ipix = ipixmin,ipixmax
           xcoord(ixcoord) = xminpix + ipix*pixwidthx
           
           !--now transform to get location of pixel in cartesians
           call coord_transform(xcoord,3,igeom,xpix,3,igeom_cartesian)
           
           !--find distances using cartesians and perform interpolation
           dy   = xpix(iycoord) - xci(iycoord)
           dx   = xpix(ixcoord) - xci(ixcoord)

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

!--------------------------------------------------------------------------
!
!   ** In this version 3D data is interpolated to a single 2D cross section
!
!   ** Note that the cross section is always taken in the z co-ordinate
!   ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x1,x2,x3 (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!            cross section location: zslice
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Daniel Price, Monash University 2nd May 2017
!--------------------------------------------------------------------------

subroutine interpolate3D_xsec_geom(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise,igeom, &
     iplotx,iploty,iplotz,ix)

  use kernels, only:cnormk3D,wfunc
  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zslice
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  integer, intent(in) :: igeom,iplotx,iploty,iplotz
  integer, dimension(3), intent(in) :: ix
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ixcoord,iycoord,izcoord,ierr
  real :: hi,hi1,radkern,q2,wab,const,hi21
  real :: termnorm,term,dx,dx2,dy,dy2,dz,dz2
  real :: xmax,ymax,xminpix,yminpix
  real, dimension(3) :: xcoord,xpix,xci,xi
  logical :: islengthx,islengthy,islengthz

  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised, non-cartesian)...',zslice
  else
     print*,'taking fast cross section (non-cartesian)...',zslice
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate3D_xsec: error: pixel width <= 0'
     return
  elseif (npart.le.0) then
     print*,'interpolate3D_xsec: error: npart = 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec_geom: WARNING: ignoring some or all particles with h < 0'
  endif
  const = cnormk3D

  !
  !--get information about the coordinates
  !
  call get_coord_info(iplotx,iploty,iplotz,ix(1),igeom,ixcoord,iycoord,izcoord,&
                      islengthx,islengthy,islengthz,ierr)
  if (ierr /= 0) return
  !
  !--if z coordinate is not a length, quit
  !
  if (.not.islengthz) then
     print*,'interpolate3D_xsec_geom: ERROR xsec not implemented when z is an angle'
     return
  endif

  xminpix = xmin - 0.5*pixwidthx
  yminpix = ymin - 0.5*pixwidthy
  xmax = xmin + npixx*pixwidthx
  ymax = ymin + npixy*pixwidthy
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set h related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
     !horigi = hh(i)
     !if (horigi.le.0.) cycle over_parts
     !hi = max(horigi,hmin)

     radkern = radkernel*hi ! radius of the smoothing kernel
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     hi21 = hi1*hi1
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2
     xcoord(izcoord) = 1.
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (dz2 .lt. radkernel2) then
        !
        !--get limits of contribution from particle in cartesian space
        !
        xci(1) = x(i)
        xci(2) = y(i)
        xci(3) = z(i)
        call get_pixel_limits(xci,xi,radkern,ipixmin,ipixmax,jpixmin,jpixmax,igeom,&
                              npixx,npixy,pixwidthx,pixwidthy,xmin,ymin,ixcoord,iycoord,ierr)
        if (ierr /= 0) cycle over_parts

        termnorm = const*weight(i)
        term = termnorm*dat(i) !/rescalefac
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           xcoord(iycoord) = yminpix + jpix*pixwidthy
           do ipix = ipixmin,ipixmax
              xcoord(ixcoord) = xminpix + ipix*pixwidthx

              !--now transform to get location of pixel in cartesians
              call coord_transform(xcoord,3,igeom,xpix,3,igeom_cartesian)

              !--find distances using cartesians and perform interpolation
              dy   = xpix(iycoord) - xci(iycoord)
              dx   = xpix(ixcoord) - xci(ixcoord)

              dx2  = dx*dx
              dy2  = dy*dy
              q2 = (dx*dx + dy*dy + dz2)*hi1*hi1
              !
              !--SPH kernel - integral through cubic spline
              !  interpolate from a pre-calculated table
              !
              if (q2.lt.radkernel2) then
                 wab = wfunc(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab

                 if (normalise) then
                    datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
                 endif
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
  !datsmooth = datsmooth*rescalefac

  return

end subroutine interpolate3D_xsec_geom

!--------------------------------------------------------------------------
!
!  utility routine for use in above routines
!
!--------------------------------------------------------------------------
subroutine get_coord_info(iplotx,iploty,iplotz,ix1,igeom,ixcoord,iycoord,izcoord,&
                          islengthx,islengthy,islengthz,ierr)
  integer, intent(in)  :: iplotx,iploty,iplotz,ix1,igeom
  integer, intent(out) :: ixcoord,iycoord,izcoord,ierr
  logical, intent(out) :: islengthx,islengthy,islengthz

  ierr = 0
  !
  !--get plotted coordinates in range 1->ndim
  !
  ixcoord = iplotx - ix1 + 1
  iycoord = iploty - ix1 + 1
  izcoord = iplotz - ix1 + 1
  if (ixcoord.le.0 .or. ixcoord.gt.3) then
     print*,' ERROR finding x coordinate offset, cannot render'
     ierr = 1
  endif
  if (iycoord.le.0 .or. iycoord.gt.3) then
     print*,' ERROR finding y coordinate offset, cannot render'
     ierr = 2
  endif
  if (izcoord.le.0 .or. izcoord.gt.3) then
     print*,' ERROR finding y coordinate offset, cannot render'
     ierr = 3
  endif
  !
  !--check if coordinate is a length (i.e., not an angle)
  !
  islengthx = coord_is_length(ixcoord,igeom)
  islengthy = coord_is_length(iycoord,igeom)
  islengthz = coord_is_length(izcoord,igeom)

end subroutine get_coord_info

!--------------------------------------------------------------------------
!
!  utility routine for use in above routines
!
!--------------------------------------------------------------------------
subroutine get_pixel_limits(xci,xi,radkern,ipixmin,ipixmax,jpixmin,jpixmax,igeom,&
                            npixx,npixy,pixwidthx,pixwidthy,xmin,ymin,ixcoord,iycoord,ierr)  
  real, intent(in)  :: xci(3),radkern,pixwidthx,pixwidthy,xmin,ymin
  real, intent(out) :: xi(3)
  integer, intent(out) :: ipixmin,ipixmax,jpixmin,jpixmax,ierr
  integer, intent(in)  :: igeom,npixx,npixy,ixcoord,iycoord
  real :: xpixmin(3),xpixmax(3)
  
  ierr = 0
  
  xpixmin(:) = xci(:) - radkern
  xpixmax(:) = xci(:) + radkern
  !
  !--transform these into limits of the contributions 
  !  in the new coordinate system
  !
  call get_coord_limits(radkern,xci,xi,xpixmin,xpixmax,igeom)
  !
  !--now work out contributions to pixels in the the transformed space
  !
  ipixmax = int((xpixmax(ixcoord) - xmin)/pixwidthx)
  if (ipixmax.lt.1) ierr = 1
  jpixmax = int((xpixmax(iycoord) - ymin)/pixwidthy)
  if (jpixmax.lt.1) ierr = 2

  ipixmin = int((xpixmin(ixcoord) - xmin)/pixwidthx)
  if (ipixmin.gt.npixx) ierr = 3
  jpixmin = int((xpixmin(iycoord) - ymin)/pixwidthy)
  if (jpixmin.gt.npixy) ierr = 4

  if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
  if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
  if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
  if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)

end subroutine get_pixel_limits

end module projections3Dgeom
