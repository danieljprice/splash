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
subroutine get_lightcurve(ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,&
                          lum,rphoto,temp,specfile)
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
 use physcon,               only:steboltz,pi,au,rsun=>solarrcgs,Lsun,c,&
                                 cm_to_nm,keV_to_Hz,keV_to_erg,kboltz
 use write_pixmap,          only:write_pixmap_ascii
 use filenames,             only:tagline
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: lum,rphoto,temp
 character(len=*), intent(in) :: specfile
 integer :: n,isinktype,npixx,npixy,ierr,j,i,nfreq,iu1,iu2,iu3
 real, dimension(3) :: xmin,xmax
 real, dimension(:),   allocatable :: weight,x,y,z,flux,opacity
 real, dimension(:),   allocatable :: freq,spectrum,bb_spectrum
 real, dimension(:,:), allocatable :: img,taupix,flux_nu
 real, dimension(:,:,:), allocatable :: img_nu
 real :: zobs,dzobs,dx,dy,area,freqmin,freqmax,lam_max,Tc,freq_max,bb_scale
 real :: lum_x,lum_bb,r_bb,Tx,bb_scalex

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
 !--set number of pixels and pixel scale in each direction
 !
 !do j=6,12
 npixx = npix !2**j
 if (npixx < 8) npixx = 1024
 dx = (xmax(1)-xmin(1))/npixx
 npixy = int((xmax(2)-xmin(2) - 0.5*dx)/dx) + 1
 dy = (xmax(2)-xmin(2))/npixy
 print "(a,i0,a,i0,a)",' Using ',npixx,' x ',npixy,' pixels'
 print "(2(1x,a,es10.3,'->',es10.3,a,/))",'x = [',xmin(1),xmax(1),']','y = [',xmin(2),xmax(2),']'
 !
 !--allocate memory for image
 !
 if (allocated(img) .or. allocated(taupix)) deallocate(img,taupix)
 allocate(img(npixx,npixy),taupix(npixx,npixy))
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

 ! frequency-dependent version
 nfreq = 128
 freqmin = 1e8
 freqmax = 1e22
 freq = logspace(nfreq,freqmin,freqmax)  ! frequency grid in Hz
 allocate(flux_nu(nfreq,n))
 do i=1,n
    flux_nu(:,i) = B_nu(dat(i,itemp),freq)
!    flux_nu(:,i) = B_nu(5.778e3,freq)
 enddo
 if (allocated(img_nu)) deallocate(img_nu)
 allocate(img_nu(nfreq,npixx,npixy))
 !
 ! raytrace SPH data to 2D image to get flux
 !
 zobs = huge(zobs)  ! no 3D perspective
 dzobs = 0.
 call interp3D_proj_opacity(x,y,z,&
      dat(1:n,ipmass),n,dat(1:n,ih),weight, &
      flux,z,icolourme(1:n), &
      n,xmin(1),xmin(2),img,taupix,npixx,npixy,&
      dx,dy,zobs,dzobs,opacity,huge(zobs),iverbose,.false.,datv=flux_nu,datvpix=img_nu)

 lum = 4.*sum(img)*dx*dy
 print*,'grey luminosity = ',lum,' erg/s'

 ! integrate flux over all frequencies to give Flux = \int F_\nu d\nu = pi \int B_nu dnu
 do j=1,npixy
    do i=1,npixx
       img(i,j) = pi*integrate_log(img_nu(1:nfreq,i,j),freq,freqmin,freqmax)
    enddo
 enddo

 lum = 4.*sum(img)*dx*dy

 ! luminosity is integrated flux
 print "(/,a,2(es10.3,a))",' L_bol = ',lum,' erg/s = ',lum/Lsun,' L_sun'

 area = count(taupix >= 1.)*dx*dy
 print "(a,1pg10.3,a)",' emitting area = ',area/au**2,' au^2'
 print "(/,a,1pg10.3,a)",' Tmax  = ',(maxval(img)/steboltz)**0.25,' K'

 ! effective temperature: total flux equals that of a blackbody at T=Teff
 temp = (lum/area/(4.*steboltz))**0.25
 freq_max = Wien_nu_from_T(temp)
 lam_max = c/freq_max*cm_to_nm
 print "(a,3(1pg10.3,a))",' Teff  = ',temp,' K: Blackbody peak at ',freq_max,' Hz / ',lam_max,' nm'

 ! get integrated spectrum from integrating over all pixels in the image
 allocate(spectrum(nfreq),bb_spectrum(nfreq))
 do i=1,nfreq
    spectrum(i) = sum(img_nu(i,1:npixx,1:npixy))*dx*dy
 enddo

 ! get colour temperature by fitting the blackbody peak
 call get_colour_temperature(spectrum,freq,Tc,freq_max,bb_scale)
 print "(a,3(1pg10.3,a))",' Tc    = ',Tc,' K: Blackbody peak at ',freq_max,' Hz / ',nu_to_lam(freq_max),' nm'

 ! effective photospheric radius, using Teff
 rphoto = sqrt(lum/(4.*pi*steboltz*temp**4))
 print "(a,2(es10.3,a))",' R_eff = ',rphoto/au,' au = ',rphoto/rsun,' rsun'

 ! L_bb and effective photospheric radius, using blackbody at T=Tc
 bb_spectrum = B_nu(Tc,freq(:))*bb_scale
 lum_bb = 4.*pi*integrate_log(bb_spectrum(1:nfreq),freq,freqmin,freqmax)
 print "(a,2(es10.3,a))",' L_bb  = ',lum_bb,' L_bb /L_x = ',lum_bb/lum_x
 r_bb = sqrt(lum_bb/(4.*pi*steboltz*Tc**4))
 print "(a,2(es10.3,a))",' R_bb  = ',r_bb/au,' au = ',r_bb/rsun,' rsun'

 ! fit blackbody to optical spectrum
 print*,' previous best fit T = ',Tc,' with bbscale = ',bb_scale
 call bb_fit(bb_spectrum,freq,lam_to_nu(400.),lam_to_nu(700.),Tc,bb_scale,min=0.5*Tc,max=10.*Tc)
 print*,' Best fit T = ',Tc,' with bbscale = ',bb_scale

 ! fit blackbody to x-ray spectrum
 Tx = Tc
 call bb_fit(spectrum,freq,real(0.3*keV_to_Hz),real(10.*keV_to_Hz),Tx,bb_scalex,min=0.1*Tc,max=100.*Tc)
 print*,' Tx = ',Tx,' kTx = ',kboltz*Tx/keV_to_erg,' keV'

 ! finally, compute x-ray luminosity
 bb_spectrum = B_nu(Tx,freq(:))*bb_scale
 lum_x = 4.*pi*integrate_log(bb_spectrum(1:nfreq),freq,freqmin,freqmax)
 print "(a,2(es10.3,a))",' L_x  = ',lum_x,' L_bb /L_x = ',lum_bb/lum_x
 r_bb = sqrt(lum_x/(4.*pi*steboltz*Tx**4))
 print "(a,2(es10.3,a))",' R_bbx = ',r_bb/au,' au = ',r_bb/rsun,' rsun'

 open(newunit=iu1,file=trim(specfile)//'.spec',status='replace',iostat=ierr)
 open(newunit=iu2,file=trim(specfile)//'.bbfit',status='replace',iostat=ierr)
 open(newunit=iu3,file=trim(specfile)//'.bbfitx',status='replace',iostat=ierr)
 write(iu1,"(a)") '# model spectrum, computed with '//trim(tagline)
 write(iu1,"(a)") '# wavelength [nm], F_\lambda'
 write(iu2,"(a)") '# best fit blackbody spectrum'
 write(iu2,"(a)") '# wavelength [nm], F_\lambda'
 write(iu3,"(a)") '# best fit blackbody spectrum to x-ray band'
 write(iu3,"(a)") '# wavelength [nm], F_\lambda'
 do i=1,nfreq
    write(iu1,*) nu_to_lam(freq(i)),spectrum(i)
    write(iu2,*) nu_to_lam(freq(i)),B_nu(Tc,freq(i))*bb_scale
    write(iu3,*) nu_to_lam(freq(i)),B_nu(Tx,freq(i))*bb_scalex
 enddo
 close(iu1)
 close(iu2)
 close(iu3)
 !print*,' bolometric luminosity = ',4.*pi*integrate_log(spectrum(1:nfreq),freq,freqmin,freqmax)

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
 real, intent(in) :: temp,nu
 real(doub_prec) :: hnu_on_kT,hnu3_on_c2

 hnu_on_kT  = hplanck*nu/(kboltz*temp)
 hnu3_on_c2 = hplanck*nu**3/c**2

 if (hnu_on_kT < 300.) then
    B_nu = 2.*hnu3_on_c2/(exp(hnu_on_kT) - 1.d0)
 else
    B_nu = tiny(0.)
 endif

end function B_nu

!---------------------------------------------------------
! Planck function derivative
! INPUT:
!    temp - temperature [K]
!    nu - frequency [Hz]
! OUTPUT:
!    dB_nu/dT - Planck function derivative w.r.t. T
!---------------------------------------------------------
real elemental function dBnu_dT(temp,nu)
 use physcon, only:c,hplanck,kboltz
 real, intent(in) :: temp,nu
 real(doub_prec) :: hnu_on_kT,hnu3_on_c2
 real :: Bnu,expterm

 hnu_on_kT  = hplanck*nu/(kboltz*temp)
 hnu3_on_c2 = hplanck*nu**3/c**2

 Bnu = B_nu(temp,nu)
 if (hnu_on_kT < 300.) then
    expterm = exp(hnu_on_kT)
    dBnu_dT = Bnu*expterm/(expterm - 1.)*hnu_on_kT/temp
 else
    dBnu_dT = tiny(0.)
 endif

end function dBnu_dT

!---------------------------------------------------------
! Wien's displacement law, in frequency
!---------------------------------------------------------
real function Wien_T_from_nu(nu) result(T)
 real, intent(in) :: nu

 T = nu/5.88e10

end function Wien_T_from_nu

!---------------------------------------------------------
! Wien's displacement law, frequency from temperature
!---------------------------------------------------------
real function Wien_nu_from_T(T) result(nu)
 real, intent(in) :: T

 nu = 5.88e10*T

end function Wien_nu_from_T

!---------------------------------------------------------
! Convert frequency in Hz to wavelength in nm
!---------------------------------------------------------
real function nu_to_lam(nu) result(lam)
 use physcon, only:c,cm_to_nm
 real, intent(in) :: nu

 lam = (c/nu)*cm_to_nm

end function nu_to_lam

!---------------------------------------------------------
! Convert frequency in Hz to wavelength in nm
!---------------------------------------------------------
real function lam_to_nu(lam) result(nu)
 use physcon, only:c,cm_to_nm
 real, intent(in) :: lam

 nu = (c/lam)*cm_to_nm

end function lam_to_nu

!---------------------------------------------------------
! colour temperature, from fitting peak of blackbody
!---------------------------------------------------------
subroutine get_colour_temperature(spectrum,freq,Tc,freq_max,bb_scale)
 real, intent(in) :: spectrum(:),freq(:)
 real, intent(out) :: Tc,freq_max,bb_scale
 integer :: imax(1),j

 imax = maxloc(spectrum)
 j = imax(1)
 freq_max = freq(j)
 Tc = Wien_T_from_nu(freq_max)

 ! scaling factor by which to shift Blackbody
 ! in y direction to match the peak flux
 bb_scale = spectrum(j)/B_nu(Tc,freq_max)

end subroutine get_colour_temperature

!--------------------------------------------------------
!+
!  Fit blackbody to spectrum over a frequency range
!  by least squares minimisation
!+
!--------------------------------------------------------
subroutine bb_fit(spec,freq,fmin,fmax,T,bbscale,min,max)
 real, intent(in) :: fmin,fmax
 real, intent(in) :: freq(:),spec(:)
 real, intent(inout) :: T
 real, intent(out) :: bbscale
 real, intent(in), optional :: min,max
 real :: Tmax,Tmin,func
 integer :: its

 Tmin = 1e3
 Tmax = 2.*T - Tmin
! if (present(min)) Tmin = min
! if (present(max)) Tmax = max
 print*,'initial min/max = ',Tmin,Tmax,0.5*(Tmin + Tmax)

 ! rootfind, using bisection method
 its = 0
 func = huge(func)
 do while(abs(func) > 1.e-20 .and. its < 100)
    T = 0.5*(Tmin + Tmax)
    func = dchisq_dT(spec,freq,fmin,fmax,T,bbscale)
    if (func > 0.) then
       Tmax = T
    else
       Tmin = T
    endif
    its = its + 1
    if (abs(func) < tiny(func)) func = -huge(func)
    !print*,its,func,'T = ',T,' min/max = ',Tmin,Tmax,bbscale
 enddo
 print*,its,func,'T = ',T,' min/max = ',Tmin,Tmax,bbscale

end subroutine bb_fit

!--------------------------------------------------------
!+
!  derivative of the chisquared error between the
!  input data (spec) and the blackbody spectrum at T
!  this is the function which we find the root of
!+
!--------------------------------------------------------
real function dchisq_dT(spec,freq,fmin,fmax,T,bbscale)
 real, intent(in) :: fmin,fmax
 real, intent(in) :: freq(:),spec(:)
 real, intent(in) :: T
 real, intent(out) :: bbscale
 integer :: i,iloc(1)
 real :: err,fmatch,Bnu

 ! fit bb_scale to match maximum value of the spectrum
 ! in the specified frequency range
 iloc = maxloc(spec,mask=(freq > fmin .and. freq < fmax))
 i = max(iloc(1),1)
 !print*,'matching at temperature ',i,Wien_T_from_nu(freq(i)),' spec = ', &
!        spec(i),' b_nu = ',B_nu(T,freq(i)),' scale=',bbscale

 bbscale = 1.e-5
 if (B_nu(T,freq(i)) > epsilon(0.)) bbscale = spec(i)/B_nu(T,freq(i))
 bbscale = max(bbscale,1e-5)

 dchisq_dT = 0.
 do i=1,size(freq)
    if (freq(i) > fmin .and. freq(i) < fmax) then
       Bnu = B_nu(T,freq(i))
       err = spec(i)/bbscale - Bnu
       dchisq_dT = dchisq_dT - 2.*err/dBnu_dT(T,freq(i))
    endif
 enddo

end function dchisq_dT

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
real function integrate_log(f,x,xmin,xmax) result(fint)
 real, intent(in) :: xmin,xmax
 real, intent(in) :: x(:),f(:)
 real :: dlogx
 integer :: n,i

 n = size(f)
 dlogx = log(xmax/xmin)/(n-1)
 fint = 0.
 do i=2,n
    if (x(i-1) > xmin .and. x(i) < xmax) then
       fint = fint + 0.5*(f(i)*x(i) + f(i-1)*x(i-1))*dlogx
    endif
 enddo

end function integrate_log

end module lightcurve
