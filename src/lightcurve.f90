module lightcurve
 use params, only:int1,doub_prec
 implicit none

 public :: get_lightcurve
 public :: get_temp_from_u

 private

contains

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
 use physcon,               only:steboltz,pi,au
 use write_pixmap,          only:write_pixmap_ascii
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: time
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: lum,rphoto,temp
 integer :: n,isinktype,npixx,npixy,ierr
 real, dimension(3) :: xmin,xmax
 real, dimension(:),   allocatable :: weight,x,y,z,flux,opacity
 real, dimension(:,:), allocatable :: datpix,taupix
 real :: zobs,dzobs,dx,dy,area

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
 npixx = 512 !npix
 if (npixx < 8) npixx = 512
 dx = (xmax(1)-xmin(1))/npixx
 npixy = int((xmax(2)-xmin(2) - 0.5*dx)/dx) + 1
 dy = (xmax(2)-xmin(2))/npixy
 print "(a,i0,a,i0,a)",' Using ',npixx,' x ',npixy,' pixels'
 print "(2(1x,a,es10.3,'->',es10.3,a,/))",'x = [',xmin(1),xmax(1),']','y = [',xmin(2),xmax(2),']'

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
 !
 ! raytrace SPH data to 2D image to get flux
 !
 zobs = huge(zobs)  ! no 3D perspective
 dzobs = 0.
 call interp3D_proj_opacity(x,y,z,&
      dat(1:n,ipmass),n,dat(1:n,ih),weight, &
      flux,z,icolourme(1:n), &
      n,xmin(1),xmin(2),datpix,taupix,npixx,npixy,&
      dx,dy,zobs,dzobs,opacity,huge(zobs),iverbose,exact_rendering)

 area = count(taupix >= 1.)*dx*dy
 print "(/,a,1pg10.3,a)",' emitting area = ',area/au**2,' au^2'
 print "(a,1pf10.3,a)",' maximum Temp  = ',(maxval(datpix)/steboltz)**0.25,' K'

 ! effective temperature: total flux equals that of a blackbody at T=Teff
 temp = (sum(datpix)/count(taupix >= 1.)/steboltz)**0.25
 print "(/,a,1pf10.3,a)",' Teff = ',temp,' K'

 ! luminosity is integrated flux
 lum = area*steboltz*temp**4
 print "(a,es10.3,a)",' luminosity = ',lum,' erg/s'

 rphoto = sqrt(lum/(4.*pi*steboltz*temp**4))
 print "(a,es10.3,a)",' radius = ',rphoto/1.49e13,' au'

 call write_pixmap_ascii(datpix,npixx,npixy,xmin(1),xmin(2),dx,minval(datpix),maxval(datpix),'temp','lum.pix',time)

end subroutine get_lightcurve

real elemental function get_temp_from_u(rho,u) result(temp)
 real(doub_prec), intent(in) :: rho,u
 real(doub_prec) :: ft,dft,dt
 real(doub_prec), parameter :: tol = 1.e-8
 real(doub_prec), parameter :: radconst = 7.5646d-15
 real(doub_prec), parameter :: kb_on_mh = real(1.38066d-16/1.67262158d-24)
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

end module lightcurve
