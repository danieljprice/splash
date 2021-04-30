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
!  Copyright (C) 2021- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module lightcurve
 use params, only:int1,doub_prec
 implicit none

 public :: get_lightcurve
 public :: get_temp_from_u

 private

contains

!---------------------------------------------------------
! routine to to compute luminosity, effective temperature
! effective blackbody radius, Reff, vs time from
! SPH particle data
!
! We solve the equation of radiative transfer along a
! ray for each pixel in the image, performing the ray
! trace through particles assuming grey blackbody
! emission (i.e. each particle emits sigma*T^4)
!
! We then compute the emitting area (area of optically
! thick material) and the effective temperature, putting
! these together to give a total luminosity
!
! Used to generate synthetic lightcurves
!---------------------------------------------------------
subroutine get_lightcurve(time,ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,lum,rphoto,temp)
 use labels,                only:ix,ih,irho,ipmass,itemp,ikappa
 use limits,                only:lim,get_particle_subset
 use interpolate3D_opacity, only:interp3D_proj_opacity
 use particle_data,         only:icolourme
 use interpolation,         only:get_n_interp,set_interpolation_weights
 use settings_data,         only:iRescale,iverbose,required,UseTypeInRenderings
 use settings_part,         only:iplotpartoftype
 use settings_render,       only:npix,inormalise=>inormalise_interpolations,&
                                 idensityweightedinterpolation,exact_rendering
 use settings_units,        only:units,unit_interp
 use physcon,               only:steboltz,pi,au,rsun=>solarrcgs,Lsun
 use write_pixmap,          only:write_pixmap_ascii
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: time
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: lum,rphoto,temp
 integer :: n,isinktype,npixx,npixy,ierr,j,nfreq
 real, dimension(3) :: xmin,xmax
 real, dimension(:),   allocatable :: weight,x,y,z,flux,opacity,freq
 real, dimension(:,:), allocatable :: datpix,taupix
 real :: zobs,dzobs,dx,dy,area,dlogfreq,freqmin,freqmax

 lum = 0.
 rphoto = 0.
 temp = 0.
 if (ndim /= 3) then
    print "(a)",' ERROR: lightcurve only works with 3 dimensional data'
    return
 endif
 if (.not. (ih > 0 .and. ipmass > 0 .and. irho > 0 .and. itemp > 0)) then
    print "(a)",' ERROR: could not locate h,mass,rho or temperature in data'
    return
 endif
 xmin(1:ndim) = lim(ix(1:ndim),1)
 xmax(1:ndim) = lim(ix(1:ndim),2)
 !
 !--set number of particles to use in the interpolation routines
 !  and allocate memory for weights
 !
 n = get_n_interp(ntypes,npartoftype,UseTypeInRenderings,iplotpartoftype,size(itype),.false.)
 allocate(weight(n),x(n),y(n),z(n),flux(n),opacity(n),stat=ierr)
 if (ierr /= 0) then
    print*,' ERROR allocating memory for interpolation weights, aborting...'
    return
 endif
 x(1:n) = dat(1:n,ix(1))
 y(1:n) = dat(1:n,ix(2))
 z(1:n) = dat(1:n,ix(3))
 !
 !--allocate memory for image
 !
 !do j=6,12
 npixx = npix !2**j
 if (npixx < 8) npixx = 1024
 dx = (xmax(1)-xmin(1))/npixx
 npixy = int((xmax(2)-xmin(2) - 0.5*dx)/dx) + 1
 dy = (xmax(2)-xmin(2))/npixy
 print "(a,i0,a,i0,a)",' Using ',npixx,' x ',npixy,' pixels'
 print "(2(1x,a,es10.3,'->',es10.3,a,/))",'x = [',xmin(1),xmax(1),']','y = [',xmin(2),xmax(2),']'

 if (allocated(datpix) .or. allocated(taupix)) deallocate(datpix,taupix)
 allocate(datpix(npixx,npixy),taupix(npixx,npixy))
 !
 !--set interpolation weights (w = m/(rho*h^ndim)
 !
 isinktype = 0 !get_sink_type(ntypes)
 call set_interpolation_weights(weight,dat,itype,(iplotpartoftype .and. UseTypeInRenderings),&
      n,npartoftype,masstype,ntypes,ncolumns,irho,ipmass,ih,ndim,iRescale,&
      idensityweightedinterpolation,inormalise,units,unit_interp,required,.false.,isinktype)
 !
 !--set default mask and apply range restrictions to data
 !
 icolourme(:) = 1
 call get_particle_subset(icolourme,dat,ncolumns)

 !
 ! specify opacity
 !
 if (ikappa > 0) then
    opacity = dat(1:n,ikappa)
 else
    print*,' WARNING: using fixed opacity kappa = 0.3 cm^2/g for lightcurve'
    opacity = 0.3
 endif
 !
 ! specify source function for each particle
 !
 flux = steboltz*dat(1:n,itemp)**4
 !flux = steboltz*5.778e3**4  ! check solar luminosity

 !nfreq = 100
 !freqmin = 1e4
 !freqmax = 1e22
 !freq = logspace(nfreq,freqmin,freqmax)  ! frequency grid in Hz
 !flux(1:nfreq) = B_nu(5.778e6,freq)
 !zobs = integrate_log(flux(1:nfreq),freq,freqmin,freqmax)
 !print "(3(es10.3))",steboltz*5.778e6**4,pi*zobs
 !
 ! raytrace SPH data to 2D image to get flux
 !
 zobs = huge(zobs)  ! no 3D perspective
 dzobs = 0.
 call interp3D_proj_opacity(x,y,z,&
      dat(1:n,ipmass),n,dat(1:n,ih),weight, &
      flux,z,icolourme(1:n), &
      n,xmin(1),xmin(2),datpix,taupix,npixx,npixy,&
      dx,dy,zobs,dzobs,opacity,huge(zobs),iverbose,.false.)

 lum = 4.*sum(datpix)*dx*dy

 ! luminosity is integrated flux
 print "(/,a,2(es10.3,a))",' luminosity = ',lum,' erg/s = ',lum/Lsun,' L_sun'

 area = count(taupix >= 1.)*dx*dy
 print "(a,1pg10.3,a)",' emitting area = ',area/au**2,' au^2'
 print "(a,1pg10.3,a)",' maximum Temp  = ',(maxval(datpix)/steboltz)**0.25,' K'

 ! effective temperature: total flux equals that of a blackbody at T=Teff
 temp = (lum/area/(4.*steboltz))**0.25
 print "(/,a,1pg10.3,a)",' Teff = ',temp,' K'

 rphoto = sqrt(lum/(4.*pi*steboltz*temp**4))
 print "(a,2(es10.3,a))",' radius = ',rphoto/au,' au = ',rphoto/rsun,' rsun'

 !if (j==12) call write_pixmap_ascii(datpix,npixx,npixy,xmin(1),xmin(2),dx,minval(datpix),maxval(datpix),'intensity','lum.pix',time)
 !enddo

end subroutine get_lightcurve

!---------------------------------------------------------
! routine to to compute temperature from
! internal energy assuming a mix of gas and radiation
! pressure, where Trad = Tgas. That is, we solve the
! quartic equation
!
!  a*T^4 + 3/2*rho*kb*T/mu = rho*u
!
! to determine the temperature from the supplied density
! and internal energy (rho, u).
! INPUT:
!    rho - density [g/cm^3]
!    u - internal energy [erg/g]
! OUTPUT:
!    temp - temperature [K]
!---------------------------------------------------------
real elemental function get_temp_from_u(rho,u) result(temp)
 use physcon, only:kb_on_mh,radconst
 real(doub_prec), intent(in) :: rho,u
 real(doub_prec) :: ft,dft,dt
 real(doub_prec), parameter :: tol = 1.e-8
 real(doub_prec), parameter :: mu = 0.6
 integer :: its

 ! Take minimum of gas and radiation temperatures as initial guess
 temp = min(u*mu/(1.5*kb_on_mh),(u*rho/radconst)**0.25)

 dt = huge(0.)
 its = 0
 do while (abs(dt) > tol*temp .and. its < 500)
    its = its + 1
    ft = u*rho - 1.5*kb_on_mh*temp*rho/mu - radconst*temp**4
    dft = - 1.5*kb_on_mh*rho/mu - 4.*radconst*temp**3
    dt = ft/dft ! Newton-Raphson
    if (temp - dt > 1.2*temp) then
       temp = 1.2*temp
    elseif (temp - dt < 0.8*temp) then
       temp = 0.8*temp
    else
       temp = temp - dt
    endif
 enddo

end function get_temp_from_u

!---------------------------------------------------------
! Planck function
! INPUT:
!    temp - temperature [K]
!    nu - frequency [Hz]
! OUTPUT:
!    B_nu - Planck function erg/s/cm^2/Hz/steradian
!---------------------------------------------------------
real elemental function B_nu(temp,nu)
 use physcon, only:c,hplanck,kboltz
 real(doub_prec), intent(in) :: temp,nu
 real(doub_prec) :: hnu_on_kT,hnu3_on_c2

 hnu_on_kT  = hplanck*nu/(kboltz*temp)
 hnu3_on_c2 = hplanck*nu**3/c**2

 if (hnu_on_kT < 300.) then
    B_nu = 2.*hnu3_on_c2/(exp(hnu_on_kT) - 1.d0)
 else
    B_nu = epsilon(0.)
 endif

end function B_nu

!--------------------------------------------------------
!+
!  Function to fill an array with equally log-spaced points
!+
!--------------------------------------------------------
function logspace(n,xmin,xmax) result(x)
 integer, intent(in) :: n
 real, intent(in)    :: xmin,xmax
 real, allocatable   :: x(:)
 integer :: i
 real    :: dx

 allocate(x(n))
 dx = log10(xmax/xmin)/real(n-1)
 do i=1,n
    x(i) = log10(xmin) + (i-1)*dx
 enddo

 x = 10.**x

end function logspace

!--------------------------------------------------------
!+
!  Integrate function on evenly spaced logarithmic grid
!  i.e. \int f(x) dx = \int x f(x) d(ln x)
!+
!--------------------------------------------------------
real(doub_prec) function integrate_log(f,x,xmin,xmax) result(fint)
 real(doub_prec), intent(in) :: xmin,xmax
 real(doub_prec), intent(in) :: x(:),f(:)
 real(doub_prec) :: dlogx
 integer :: n,i

 n = size(f)
 dlogx = log(xmax/xmin)/(n-1)
 fint = 0.
 do i=2,n
    fint = fint + 0.5*(f(i)*x(i) + f(i-1)*x(i-1))*dlogx
 enddo

end function integrate_log

end module lightcurve
