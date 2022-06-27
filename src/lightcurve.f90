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
                          lum,rphoto,temp,lum_bb,r_bb,Tc,specfile)
 use labels,                only:ix,ih,irho,ipmass,itemp,ikappa,ivx,ipmomx
 use limits,                only:lim,get_particle_subset
 use lightcurve_utils,      only:get_temp_from_u
 use interpolate3D_opacity, only:interp3D_proj_opacity
 use particle_data,         only:icolourme
 use interpolation,         only:get_n_interp,set_interpolation_weights
 use settings_data,         only:iRescale,iverbose,required,UseTypeInRenderings,ndimV
 use settings_part,         only:iplotpartoftype
 use settings_render,       only:npix,inormalise=>inormalise_interpolations,&
                                 idensityweightedinterpolation,exact_rendering
 use settings_units,        only:units,unit_interp
 use physcon,               only:steboltz,pi,au,rsun=>solarrcgs,Lsun,c,cm_to_nm
 use write_pixmap,          only:write_pixmap_ascii
 use filenames,             only:tagline
 use blackbody,             only:B_nu,logspace,Wien_nu_from_T,nu_to_lam,&
                                 integrate_log,get_colour_temperature
 use settings_xsecrot,      only:anglex,angley,anglez,taupartdepth
 use rotation,              only:rotate3D
 use system_utils,          only:get_command_flag
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: lum,rphoto,temp,lum_bb,r_bb,Tc
 character(len=*), intent(in) :: specfile
 integer :: n,isinktype,npixx,npixy,ierr,j,i,k,nfreq
 integer, parameter :: iu1 = 45
 real, dimension(3) :: xmin,xmax
 real, dimension(:),   allocatable :: weight,x,y,z,flux,opacity
 real, dimension(:),   allocatable :: freq,spectrum,bb_spectrum
 real, dimension(:,:), allocatable :: img,taupix,flux_nu,v_on_c
 real, dimension(:,:,:), allocatable :: img_nu
 real :: zobs,dzobs,dx,dy,area,freqmin,freqmax,lam_max,freq_max,bb_scale,opacity_factor
 real :: ax,ay,az,xi(3),betaz,lorentz,doppler_factor,doppler_factor_max,taui,lprev
 logical :: relativistic
 character(len=20) :: tmpfile

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
 !
 !--allow for reduction in opacity by factor of v/c for photon trapping
 !
 relativistic = .not.get_command_flag('nonrel') .and. (ivx > 0 .and. ndimV >= 3 .and. ipmomx > 0)
 if (relativistic) print "(a)",' relativistic corrections ON (use --nonrel to switch off)'

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

 if (relativistic) then
    allocate(v_on_c(3,n),stat=ierr)
    do i=1,3
       v_on_c(i,:) = dat(1:n,ivx+i-1)/c ! velocity in units of speed of light
    enddo
 endif
 if (abs(anglez)>0. .or. abs(angley)>0. .or. abs(anglex)>0.) then
    print*, 'Rotating particles around (z,y,x) by',anglez,angley,anglex
    ax = anglex*pi/180.0 ! convert degrees to radians to pass into rotate
    ay = angley*pi/180.0
    az = anglez*pi/180.0
    do i=1,n
       xi = (/x(i),y(i),z(i)/)
       call rotate3D(xi,ax,ay,az,0.,0.)
       x(i) = xi(1)
       y(i) = xi(2)
       z(i) = xi(3)
       if (allocated(v_on_c)) call rotate3D(v_on_c(:,i),ax,ay,az,0.,0.)
    enddo
 endif

 !
 !--set number of pixels and pixel scale in each direction
 !
 npixx = npix
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
    if (abs(taupartdepth - 1.0) > tiny(0.) .and. taupartdepth > 0.) then
       print "(a,es10.2)",' opacity = opacity x ',taupartdepth
       opacity = opacity*taupartdepth
    endif
 else
    opacity = 0.35
    if (taupartdepth > 0.) opacity = taupartdepth
    print "(a,1pg10.2,a)",' WARNING: using fixed opacity kappa = ',maxval(opacity),' cm^2/g for lightcurve'
 endif
 !
 ! specify source function for each particle
 !
 flux = steboltz*dat(1:n,itemp)**4 ! grey version

 ! frequency-dependent version
 nfreq = 128
 freqmin = 1e8
 freqmax = 1e22
 freq = logspace(nfreq,freqmin,freqmax)  ! frequency grid in Hz
 if (allocated(flux_nu)) deallocate(flux_nu)
 allocate(flux_nu(nfreq,n))
 doppler_factor = 1.
 opacity_factor = 1.
 doppler_factor_max = 0.
 do i=1,n
    if (relativistic) then
       betaz = 1. + v_on_c(3,i)
       lorentz = 1./sqrt(1. - dot_product(v_on_c(:,i),v_on_c(:,i)))
       doppler_factor = betaz*lorentz               ! nu / nu_0, Lorentz transform
       opacity_factor = lorentz*(1. - v_on_c(3,i))  ! nu_0 / nu, inverse Lorentz transform
       doppler_factor_max = max(doppler_factor,doppler_factor_max)
    endif
    opacity(i) = opacity(i)*opacity_factor
    flux_nu(:,i) = B_nu(dat(i,itemp),freq*doppler_factor)
 enddo
 if (relativistic) print*,' max relativistic correction=',doppler_factor_max

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
 if (allocated(spectrum)) deallocate(spectrum,bb_spectrum)
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
 print "(a,2(es10.3,a))",' L_bb  = ',lum_bb
 r_bb = sqrt(lum_bb/(4.*pi*steboltz*Tc**4))
 print "(a,2(es10.3,a))",' R_bb  = ',r_bb/au,' au = ',r_bb/rsun,' rsun'

 print "(a)",' WRITING '//trim(specfile)//'.spec'
 open(unit=iu1,file=trim(specfile)//'.spec',status='replace',iostat=ierr)
 write(iu1,"(a)") '# model spectrum, computed with '//trim(tagline)
 write(iu1,"(a)") '# frequency [Hz], F_\nu'
 do i=1,nfreq
    write(iu1,*) freq(i),spectrum(i)
 enddo
 close(iu1)

 if (get_command_flag('writepix')) then
    call write_pixmap_ascii(img,npixx,npixy,xmin(1),xmin(2),dx,&
         minval(img),maxval(img),'intensity','lum.pix',0.)
 endif

end subroutine get_lightcurve

end module lightcurve
