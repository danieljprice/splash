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

module interpolate3D_opacity
 use projections3D, only:wfromtable,coltable
 use kernels,       only:radkernel2
 use sort,          only:indexx
 use interpolation, only:weight_sink
 implicit none

contains
!--------------------------------------------------------------------------
! $Id: interpolate3D_opacity.f90,v 1.16 2007/11/20 17:05:35 dprice Exp $
!
!     subroutine to do a ray trace through the particle data
!
!     we use the radiation transport equation along a ray, that is
!     the change in intensity from one side of a particle to the other is
!     given by:
!
!     I_nu = I_nu(0) exp(-tau_i) + S_nu (1 - exp(-tau_i))
!
!     where tau_i is the integrated optical depth through the particle,
!     and S_nu is the colour calculated from a colour table for the rendered data.
!     We calculate an intensity in red, green and blue for colour plots.
!
!     tau_i = kappa \int rho dz
!
!     this is calculated using the SPH kernel for rho, so for each pixel
!     the optical depth is incremented as the sum
!
!     tau_i = kappa \sum_j m_j \int W dz
!
!     where \int W dz is the SPH kernel integrated along one spatial dimension.
!     This is interpolated from a pre-calculated table (see module projections3D for this).
!
!     kappa is the monochromatic mass extinction coefficient
!     (particle cross section per unit mass) and is a constant for all particles
!     which must be given as input (although see below for calculations of a
!     meaningful values for kappa in terms of "surface depth in units of smoothing lengths")
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            particle masses       : pmass (npmass)
!            smoothing lengths     : hh    (npart)
!            weight                : m/(h^3 rho) (not used, but skips particles with w <= 0)
!            scalar data to smooth : dat   (npart)
!
!     Particle masses can be sent in as either a single scalar (npmass = 1)
!      or as an array of length npart (npmass=npart)
!
!     Settings: zobs, dz1 : settings for 3D projection
!               rkappa    : particle cross section per unit mass
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!             brightness array         : brightness (npixx,npixy)
!
!--------------------------------------------------------------------------

subroutine interp3D_proj_opacity(x,y,z,pmass,npmass,hh,weight,dat,zorig,itype,npart, &
     xmin,ymin,datsmooth,brightness,npixx,npixy,pixwidth,zobserver,dscreenfromobserver, &
     rkappa,zcut,iverbose)

 implicit none
 real, parameter :: pi=3.1415926536
 integer, intent(in) :: npart,npixx,npixy,npmass,iverbose
 real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat,zorig
 real, intent(in), dimension(npmass) :: pmass
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,pixwidth,zobserver,dscreenfromobserver, &
                      zcut,rkappa
 real, dimension(npixx,npixy), intent(out) :: datsmooth, brightness

 integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,nused,nsink
 integer :: iprintinterval, iprintnext,itmin
 integer, dimension(npart) :: iorder
 integer(kind=selected_int_kind(12)) :: ipart
 real :: hi,hi1,hi21,radkern,q2,wab,pmassav
 real :: term,dy,dy2,ypix,zfrac,hav,zcutoff
 real :: fopacity,tau,rkappatemp,termi,xi,yi
 real :: t_start,t_end,t_used,tsec
 logical :: iprintprogress,adjustzperspective,rendersink
 real, dimension(npixx) :: xpix,dx2i
 real :: xminpix,yminpix
!#ifdef _OPENMP
!  integer :: OMP_GET_NUM_THREADS
!#else
 integer(kind=selected_int_kind(12)) :: iprogress
!#endif

 datsmooth = 0.
 term = 0.
 brightness = 0.
 if (pixwidth <= 0.) then
    if (iverbose >= -1) print "(a)",'interpolate3D_opacity: error: pixel width <= 0'
    return
 endif
 if (any(hh(1:npart) <= tiny(hh)) .and. iverbose >= -1) then
    print*,'interpolate3D_opacity: warning: ignoring some or all particles with h < 0'
 endif
 !--check that npmass is sensible
 if (npmass < 1 .or. npmass > npart) then
    if (iverbose >= -1) print*,'interpolate3D_opacity: ERROR in input number of particle masses '
    return
 endif
 !--these values for npmass are not sensible but the routine will still work
 if (npmass /= 1 .and. npmass /= npart) then
    if (iverbose >= -1) print*,'WARNING: interpolate3D_opacity: number of particle masses input =',npmass
 endif

 if (abs(dscreenfromobserver) > tiny(dscreenfromobserver)) then
    adjustzperspective = .true.
    zcutoff = zobserver
 else
    adjustzperspective = .false.
    zcutoff = huge(zobserver)
 endif

!
!--kappa is the opacity in units of length^2/mass
!  sent as an input parameter as it should be kept constant throughout the simulation
!
!  However we compute a reasonable estimate below based on the current plot so that
!  we can give the "actual" optical depth for the current frame in terms of number of
!  smoothing lengths. This is purely for diagnostic purposes only.
!
!--calculate average h
 hav = sum(hh(1:npart))/real(npart)
!--average particle mass
 pmassav = sum(pmass(1:npmass))/real(npmass)
 rkappatemp = pi*hav*hav/(pmassav*coltable(0))
 if (iverbose >= 0) print "(1x,a,g8.2,a)",'ray tracing: surface depth ~ ',rkappatemp/rkappa,' smoothing lengths'
 !print "(1x,a,f6.2,a)",'typical surface optical depth is ~',rkappatemp/rkappa,' smoothing lengths'
 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = ((npart  >=  1000000) .or. (npixx*npixy  > 500000)) .and. (iverbose >= 0)
 !
 !--loop over particles
 !
 iprintinterval = 25
 if (npart >= 1e7) iprintinterval = 10
 iprintnext = iprintinterval
!
!--get starting CPU time
!
 call cpu_time(t_start)
!
!--first sort the particles in z so that we do the opacity in the correct order
!
 call indexx(npart,z,iorder)
!
!--store x value for each pixel (for optimisation)
!
 xminpix = xmin - 0.5*pixwidth
 yminpix = ymin - 0.5*pixwidth
 do ipix=1,npixx
    xpix(ipix) = xminpix + ipix*pixwidth
 enddo

 nused = 0
 nsink = 0

!!$omp parallel default(none) &
!!$omp SHARED(hh,z,x,y,zorig,pmass,dat,itype,datsmooth,npmass,npart) &
!!$omp SHARED(xmin,ymin,xminpix,yminpix,xpix,pixwidth) &
!!$omp SHARED(npixx,npixy,dscreenfromobserver,zobserver,adjustzperspective) &
!!$omp SHARED(zcut,zcutoff,iorder,rkappa,brightness) &
!!$omp PRIVATE(hi,zfrac,xi,yi,radkern) &
!!$omp PRIVATE(hi1,hi21,term,termi) &
!!$omp PRIVATE(ipixmin,ipixmax,jpixmin,jpixmax) &
!!$omp PRIVATE(dx2i,q2,ypix,dy,dy2,wab) &
!!$omp PRIVATE(ipart,i,ipix,jpix,tau,fopacity) &
!!$omp REDUCTION(+:nused)
!!$omp MASTER
!#ifdef _OPENMP
!  print "(1x,a,i3,a)",'Using ',OMP_GET_NUM_THREADS(),' cpus'
!#endif
!!$omp END MASTER

!!$omp DO ORDERED SCHEDULE(dynamic)
 over_particles: do ipart=1,npart
    !
    !--report on progress
    !
!#ifndef _OPENMP
    if (iprintprogress) then
       iprogress = 100*(ipart/real(npart))
       if (iprogress >= iprintnext) then
          write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,ipart
          iprintnext = iprintnext + iprintinterval
       endif
    endif
!#endif
    !
    !--render in order from back to front
    !
    i = iorder(ipart)
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0) cycle over_particles

    !
    !--skip particles with weight < 0
    !  but not if weight == weight_sink (=-1)
    !
    rendersink = .false.
    if (abs(weight(i) - weight_sink) < tiny(0.)) then
       rendersink = .true.
    elseif (weight(i) <= 0.) then
       cycle over_particles
    endif

    !
    !--allow slicing [take only particles with z(unrotated) < zcut]
    !
    particle_within_zcut: if (zorig(i) < zcut .and. z(i) < zcutoff) then

       !  count particles within slice
       nused = nused + 1
       !
       !--adjust h according to 3D perspective
       !  need to be careful -- the kernel quantities
       !  change with z (e.g. radkern, r^2/h^2)
       !  but *not* the 1/h^2 in tau (because the change in 1/h^2 in tau
       !  would be cancelled by the corresponding change to h^2 in kappa)
       !
       hi = hh(i)
       if (hi <= 0.) then
          cycle over_particles
       elseif (adjustzperspective) then
          zfrac = abs(dscreenfromobserver/(z(i)-zobserver))
          hi = hi*zfrac
       endif

       !--these are the quantities used in the kernel r^2/h^2
       radkern = 2.*hi
       hi1 = 1./hi
       hi21 = hi1*hi1
       !--this is the term which multiplies tau
       if (npmass==npart) then
          term = pmass(i)/(hh(i)*hh(i))
       else
          term = pmass(1)/(hh(i)*hh(i))
       endif
       !
       !--determine colour contribution of current point
       !  (work out position in colour table)
       !
!     dati = dat(i)
       xi = x(i)
       yi = y(i)
       termi = dat(i)
       !
       !--sink particles can have weight set to -1
       !  indicating that we should include them in the rendering
       !
       if (rendersink) then
          termi = pmass(i)/(4./3.*pi*hh(i)**3)  ! define "density" of a sink
          nsink = nsink + 1
       endif
       !
       !--for each particle work out which pixels it contributes to
       !
       ipixmin = int((xi - radkern - xmin)/pixwidth)
       jpixmin = int((yi - radkern - ymin)/pixwidth)
       ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
       jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

       if (ipixmin < 1) ipixmin = 1  ! make sure they only contribute
       if (jpixmin < 1) jpixmin = 1  ! to pixels in the image
       if (ipixmax > npixx) ipixmax = npixx ! (note that this optimises
       if (jpixmax > npixy) jpixmax = npixy !  much better than using min/max)
       !
       !--precalculate an array of dx2 for this particle (optimisation)
       !
       do ipix=ipixmin,ipixmax
          dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
       enddo

       !
       !--loop over pixels, adding the contribution from this particle
       !
       do jpix = jpixmin,jpixmax
          ypix = yminpix + jpix*pixwidth
          dy = ypix - yi
          dy2 = dy*dy*hi21
          do ipix = ipixmin,ipixmax
             q2 = dx2i(ipix) + dy2
             !
             !--SPH kernel - integral through cubic spline
             !  interpolate from a pre-calculated table
             !
             if (q2 < radkernel2) then
                wab = wfromtable(q2)
                !
                !--get incremental tau for this pixel from the integrated SPH kernel
                !
                tau = rkappa*wab*term
                fopacity = 1. - exp(-tau)
                !
                !--render, obscuring previously drawn pixels by relevant amount
                !  also calculate total brightness (`transparency') of each pixel
                !
                datsmooth(ipix,jpix) = (1.-fopacity)*datsmooth(ipix,jpix) + fopacity*termi
                brightness(ipix,jpix) = brightness(ipix,jpix) + fopacity
             endif

          enddo
       enddo

    endif particle_within_zcut

 enddo over_particles
!!$omp end do
!!$omp end parallel

!
!--get ending CPU time
!
 if (nsink > 99) then
    if (iverbose >= 0) print*,'rendered ',nsink,' sink particles'
 elseif (nsink > 0) then
    if (iverbose >= 0) print "(1x,a,i2,a)",'rendered ',nsink,' sink particles'
 endif
 call cpu_time(t_end)
 t_used = t_end - t_start
 if (t_used > 60.) then
    itmin = int(t_used/60.)
    tsec = t_used - (itmin*60.)
    if (iverbose >= 0) print "(1x,a,i4,a,f5.2,1x,a)",'completed in',itmin,' min ',tsec,'s'
 elseif (t_used > 10.) then
    if (iverbose >= 0) print "(1x,a,f5.2,1x,a)",'completed in ',t_used,'s'
 endif
 if (zcut < huge(zcut) .and. iverbose >= 0) print*,'slice contains ',nused,' of ',npart,' particles'

 return

end subroutine interp3D_proj_opacity

end module interpolate3D_opacity
