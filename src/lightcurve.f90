module lightcurve
 use params, only:int1,doub_prec
 implicit none

 public :: get_lightcurve

 private

contains

subroutine get_lightcurve(time,ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,lum,rphoto,temp)
 use convert_grid, only:convert_to_grid
 use labels,       only:iutherm,ix
 use limits,       only:lim
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: time
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: lum,rphoto,temp
 real, allocatable :: rhogrid(:,:,:),ugrid(:,:,:),tgrid(:,:,:)
 integer :: nlos,idim,nx,ny,nz
 real :: xmin(ndim),xmax(ndim)

 if (ndim /= 3) then
    print "(a)",' ERROR: lightcurve only works with 3 dimensional data'
    return
 endif
 !
 ! interpolate SPH data to 3D grid to get density and thermal energy
 !
 call convert_to_grid(time,dat,ntypes,npartoftype,masstype,itype,ncolumns,'none',&
                      'none',.false.,icols=(/iutherm/),rhogrid=rhogrid,dat3D=ugrid)
 !
 ! convert thermal energy to temperature
 !
 print*,' max density = ',maxval(rhogrid)
 print*,' max utherm  = ',maxval(ugrid)
 tgrid = get_temp_from_u(rhogrid,ugrid)
 print*,' max temp    = ',maxval(tgrid)
 if (allocated(ugrid)) deallocate(ugrid)
 !
 ! find the photosphere and get the temperature and luminosity
 !
 nx = size(rhogrid,1)
 ny = size(rhogrid,2)
 nz = size(rhogrid,3)
 xmin(1:ndim) = lim(ix(1:ndim),1)
 xmax(1:ndim) = lim(ix(1:ndim),2)
 do idim=3,3
    call get_photosphere(idim,1,nx,ny,nz,rhogrid,tgrid,xmin,xmax,lum,rphoto,temp)
    !call get_photosphere(idim,-1,rhogrid,ugrid,lum,rphoto,temp)
 enddo
 !
 ! clean up memory
 !
 if (allocated(rhogrid)) deallocate(rhogrid)
 if (allocated(tgrid)) deallocate(tgrid)

end subroutine get_lightcurve

subroutine get_photosphere(idim,idir,nx,ny,nz,rho,temp,xmin,xmax,lum,rphoto,tphoto)
 integer, intent(in)  :: idim,idir,nx,ny,nz
 real,    intent(in)  :: rho(nx,ny,nz),temp(nx,ny,nz),xmin(:),xmax(:)
 real,    intent(out) :: lum,rphoto,tphoto
 real, allocatable :: tlast(:,:),tau(:,:)
 real :: dx,dy,dz,tau_crit,kappa
 integer :: k,j,i,ierr,npixels_in_photosphere
 real, parameter :: steboltz = 5.67051e-5
 real, parameter :: pi = 4.*atan(1.)

 tau_crit = 1.
 kappa = 0.2

 select case(idim)
 case(3)
    dz = (xmax(3) - xmin(3))/nz
    dy = (xmax(2) - xmin(2))/ny
    dx = (xmax(1) - xmin(1))/ny
    allocate(tlast(nx,ny),tau(nx,ny),stat=ierr)
    tau = 0.
    tlast = 0.
    do k=nz,1,-1
       where(tau < tau_crit)
          tlast(:,:) = temp(:,:,k)
          tau(:,:)   = tau(:,:) + kappa*rho(:,:,k)*dz
       end where
    enddo
    npixels_in_photosphere = count(tau > 1.)
    lum = sum(steboltz*tlast**4*dx*dy,mask=(tau > 1.))
    tphoto = sum(tlast,mask=(tau > 1.))/npixels_in_photosphere
    print "(1x,g8.2,a)",100.*npixels_in_photosphere/real(nx*ny),'% of image is optically thick'
    rphoto = sqrt(lum/(4.*pi*steboltz*tphoto**4))
    print*,' tau (max,ave,min) = ',maxval(tau),sum(tau)/(nx*ny),minval(tau)
    deallocate(tlast,tau)
 end select

end subroutine get_photosphere

real elemental function get_temp_from_u(rho,u) result(temp)
 real, intent(in) :: rho,u
 real :: ft,dft,dt
 real, parameter :: tol = 1.e-8
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
