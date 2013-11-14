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
     rkappa,zcut)

  implicit none
  real, parameter :: pi=3.1415926536
  integer, intent(in) :: npart,npixx,npixy,npmass
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
  print "(1x,a)",'ray tracing from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_opacity: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_opacity: warning: ignoring some or all particles with h < 0'
  endif
  !--check that npmass is sensible
  if (npmass.lt.1 .or. npmass.gt.npart) then
     print*,'interpolate3D_opacity: ERROR in input number of particle masses '
     return
  endif
  !--these values for npmass are not sensible but the routine will still work
  if (npmass.ne.1 .and. npmass.ne.npart) then
     print*,'WARNING: interpolate3D_opacity: number of particle masses input =',npmass
  endif

  if (abs(dscreenfromobserver).gt.tiny(dscreenfromobserver)) then
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
  print*,'average h = ',hav,' average mass = ',pmassav
  print "(1x,a,f6.2,a)",'typical surface optical depth is ~',rkappatemp/rkappa,' smoothing lengths'
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

!!$OMP PARALLEL default(none) &
!!$OMP SHARED(hh,z,x,y,zorig,pmass,dat,itype,datsmooth,npmass,npart) &
!!$OMP SHARED(xmin,ymin,xminpix,yminpix,xpix,pixwidth) &
!!$OMP SHARED(npixx,npixy,dscreenfromobserver,zobserver,adjustzperspective) &
!!$OMP SHARED(zcut,zcutoff,iorder,rkappa,brightness) &
!!$OMP PRIVATE(hi,zfrac,xi,yi,radkern) &
!!$OMP PRIVATE(hi1,hi21,term,termi) &
!!$OMP PRIVATE(ipixmin,ipixmax,jpixmin,jpixmax) &
!!$OMP PRIVATE(dx2i,q2,ypix,dy,dy2,wab) &
!!$OMP PRIVATE(ipart,i,ipix,jpix,tau,fopacity) &
!!$OMP REDUCTION(+:nused)
!!$OMP MASTER
!#ifdef _OPENMP
!  print "(1x,a,i3,a)",'Using ',OMP_GET_NUM_THREADS(),' cpus'
!#endif
!!$OMP END MASTER

!!$OMP DO ORDERED SCHEDULE(dynamic)
  over_particles: do ipart=1,npart
     !
     !--report on progress
     !
!#ifndef _OPENMP
     if (iprintprogress) then
        iprogress = 100*(ipart/npart)
        if (iprogress.ge.iprintnext) then
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
     particle_within_zcut: if (zorig(i).lt.zcut .and. z(i).lt.zcutoff) then

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
     if (hi.le.0.) then
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
     if (npmass.eq.npart) then
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
           if (q2.lt.radkernel2) then
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
!!$OMP END DO
!!$OMP END PARALLEL

!
!--get ending CPU time
!
  if (nsink > 99) then
     print*,'rendered ',nsink,' sink particles'
  elseif (nsink > 0) then
     print "(1x,a,i2,a)",'rendered ',nsink,' sink particles'
  endif
  call cpu_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.60.) then
     itmin = int(t_used/60.)
     tsec = t_used - (itmin*60.)
     print "(1x,a,i4,a,f5.2,1x,a)",'completed in',itmin,' min ',tsec,'s'
  else
     print "(1x,a,f5.2,1x,a)",'completed in ',t_used,'s'
  endif
  if (zcut.lt.huge(zcut)) print*,'slice contains ',nused,' of ',npart,' particles'

  return

end subroutine interp3D_proj_opacity

subroutine interp3D_proj_opacity_writeppm(datsmooth,brightness,npixx,npixy,datmin,datmax,istep)
  use colours, only:rgbtable,ncolours
  implicit none
  integer, intent(in) :: npixx,npixy
  real, intent(in), dimension(npixx,npixy) :: datsmooth,brightness
  real, intent(in) :: datmin,datmax
  integer, intent(in) :: istep
  character(len=120) :: filename
!  real, dimension(3,npixx,npixy) :: rgb
  real, dimension(3) :: rgbi,drgb
  real :: dati,ddatrange,datfraci,ftable
  integer :: ipix,jpix,ir,ib,ig,ierr,maxcolour,indexi
!
!--check for errors
!
  if (abs(datmax-datmin).gt.tiny(datmin)) then
     ddatrange = 1./abs(datmax-datmin)
  else
     print "(a)",'error: datmin=datmax : pointless writing ppm file'
     return
  endif
!
!--write PPM--
!
  write(filename,"(a,i5.5,a)") 'splash_',istep,'.ppm'
  open(unit=78,file=filename,status='replace',form='formatted',iostat=ierr)
  if (ierr /=0) then
     print*,'error opening ppm file'
     return
  endif
  print "(a)", 'writing to file '//trim(filename)
!
!--PPM header
!
  maxcolour = 255
  write(78,"(a)") 'P3'
  write(78,"(a)") '# splash.ppm created by splash (c) 2005-2007 Daniel Price'
  write(78,"(i4,1x,i4)") npixx, npixy
  write(78,"(i3)") maxcolour
!--pixel information
  do jpix = npixy,1,-1
     do ipix = 1,npixx

        dati = datsmooth(ipix,jpix)
        datfraci = (dati - datmin)*ddatrange
        datfraci = max(datfraci,0.)
        datfraci = min(datfraci,1.)
        !--define colour for current particle
        ftable = datfraci*ncolours
        indexi = int(ftable) + 1
        indexi = min(indexi,ncolours)
        if (indexi.lt.ncolours) then
        !--do linear interpolation from colour table
           drgb(:) = rgbtable(:,indexi+1) - rgbtable(:,indexi)
           rgbi(:) = rgbtable(:,indexi) + (ftable - int(ftable))*drgb(:)
        else
           rgbi(:) = rgbtable(:,indexi)
        endif
        rgbi(:) = rgbi(:)*min(brightness(ipix,jpix),1.0)
        ir = max(min(int(rgbi(1)*maxcolour),maxcolour),0)
        ig = max(min(int(rgbi(2)*maxcolour),maxcolour),0)
        ib = max(min(int(rgbi(3)*maxcolour),maxcolour),0)
        write(78,"(i3,1x,i3,1x,i3,2x)") ir,ig,ib
     enddo
  enddo
  close(unit=78)

  return
end subroutine interp3D_proj_opacity_writeppm

end module interpolate3D_opacity
