module lightcurve
 use params, only:int1,doub_prec
 implicit none

 ! public :: get_lightcurve
 public :: get_temp_from_u,ionisation_fraction

 private

contains

! subroutine get_lightcurve(time,ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,lum,rphoto,temp)
!  use convert_grid, only:convert_to_grid
!  use labels,       only:iutherm,ix
!  use limits,       only:lim
!  integer, intent(in)  :: ncolumns,ntypes,ndim
!  integer, intent(in)  :: npartoftype(:)
!  integer(kind=int1), intent(in) :: itype(:)
!  real,    intent(in)  :: time
!  real,    intent(in)  :: masstype(:)
!  real,    intent(in)  :: dat(:,:)
!  real,    intent(out) :: lum,rphoto,temp
!  real, allocatable :: rhogrid(:,:,:),ugrid(:,:,:),tgrid(:,:,:)
!  integer :: nlos,idim,nx,ny,nz
!  real :: xmin(ndim),xmax(ndim)
!
!  if (ndim /= 3) then
!     print "(a)",' ERROR: lightcurve only works with 3 dimensional data'
!     return
!  endif
!  !
!  ! interpolate SPH data to 3D grid to get density and thermal energy
!  !
!  call convert_to_grid(time,dat,ntypes,npartoftype,masstype,itype,ncolumns,'none',&
!                       'none',.false.,icols=(/iutherm/),rhogrid=rhogrid,dat3D=ugrid)
!  !
!  ! convert thermal energy to temperature
!  !
!  print*,' max density = ',maxval(rhogrid)
!  print*,' max utherm  = ',maxval(ugrid)
!  tgrid = get_temp_from_u(rhogrid,ugrid)
!  print*,' max temp    = ',maxval(tgrid)
!  if (allocated(ugrid)) deallocate(ugrid)
!  !
!  ! find the photosphere and get the temperature and luminosity
!  !
!  nx = size(rhogrid,1)
!  ny = size(rhogrid,2)
!  nz = size(rhogrid,3)
!  xmin(1:ndim) = lim(ix(1:ndim),1)
!  xmax(1:ndim) = lim(ix(1:ndim),2)
!  do idim=3,3
!     call get_photosphere(idim,1,nx,ny,nz,rhogrid,tgrid,xmin,xmax,lum,rphoto,temp)
!     !call get_photosphere(idim,-1,rhogrid,ugrid,lum,rphoto,temp)
!  enddo
!  !
!  ! clean up memory
!  !
!  if (allocated(rhogrid)) deallocate(rhogrid)
!  if (allocated(tgrid)) deallocate(tgrid)
!
! end subroutine get_lightcurve

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


! subroutine calc_energy(particlemass,poten,xyzh,vxyzu,xyzmh_ptmass,phii,epoti,ekini,einti,etoti)
!  ! Warning: Do not sum epoti or etoti as it is to obtain a total energy; this would not give the correct
!  !          total energy due to complications related to double-counting.
!  use ptmass, only:get_accel_sink_gas
!  use part,   only:nptmass
!  real, intent(in)                       :: particlemass
!  real(4), intent(in)                    :: poten
!  real, dimension(4), intent(in)         :: xyzh,vxyzu
!  real, dimension(5,nptmass), intent(in) :: xyzmh_ptmass
!  real, intent(out)                      :: phii,epoti,ekini,einti,etoti
!  real                                   :: fxi,fyi,fzi

!  phii = 0.0

!  call get_accel_sink_gas(nptmass,xyzh(1),xyzh(2),xyzh(3),xyzh(4),xyzmh_ptmass,fxi,fyi,fzi,phii)

!  epoti = 2.*poten + particlemass * phii ! For individual particles, need to multiply 2 to poten to get GmM/r
!  ekini = particlemass * 0.5 * dot_product(vxyzu(1:3),vxyzu(1:3))
!  einti = particlemass * vxyzu(4)
!  etoti = epoti + ekini + einti
! end subroutine calc_energy


!----------------------------------------------------------------
!+
!  Solves three Saha equations simultaneously to return ion
!  fractions of hydrogen and helium. Assumes inputs in cgs units
!+
!----------------------------------------------------------------
subroutine ionisation_fraction(dens,temp,X,Y,xh0,xh1,xhe0,xhe1,xhe2)
 real, intent(in) :: dens, temp, X, Y
 real, intent(out):: xh0, xh1, xhe0, xhe1, xhe2
 real             :: n, nh, nhe
 real             :: A, B, C, const
 real             :: xh1g, xhe1g, xhe2g
 real             :: f, g, h
 real, parameter  :: chih0 = 13.6, chihe0 = 24.6, chihe1 = 54.4
 real, dimension(3,3) :: M, M_inv
 real, dimension(3) :: dx
 integer          :: i
 real, parameter :: twopi=6.2831853072d0,kboltz=1.38066d-16,eV=1.60219d-12,&
                    planckh=6.6260755d-27,mass_electron_cgs=9.10938291d-28,mass_proton_cgs=1.67262158d-24

 nh = X * dens / mass_proton_cgs
 nhe = Y * dens / (4. * mass_proton_cgs)
 n = nh + nhe

 const = (sqrt(twopi * mass_electron_cgs * kboltz) / planckh)**3 / n

 A = 1. * const * temp**(1.5) * exp(-chih0 * eV / (kboltz * temp))
 B = 4. * const * temp**(1.5) * exp(-chihe0 * eV / (kboltz * temp))
 C = 1. * const * temp**(1.5) * exp(-chihe1 * eV / (kboltz * temp))

 xh1g = 0.4
 xhe1g = 0.3
 xhe2g = 0.2

 do i=1,50
    f = xh1g * (xh1g + xhe1g + 2*xhe2g) - A * ((nh/n) - xh1g)
    g = xhe1g * (xh1g + xhe1g + 2*xhe2g) - B * ((nhe/n) - xhe1g - xhe2g)
    h = xhe2g * (xh1g + xhe1g + 2*xhe2g) - C * xhe1g

    M(1,:) = (/ 2*xh1g + xhe1g + 2*xhe2g + A, xh1g, 2*xh1g /)
    M(2,:) = (/ xhe1g, xh1g + 2*xhe1g + 2*xhe2g + B, 2*xhe1g + B /)
    M(3,:) = (/ xhe2g, xhe2g - C, xh1g + xhe1g + 4*xhe2g /)

    call minv(M, M_inv)

    dx = matmul(M_inv, (/ -f, -g, -h/))

    xh1g = xh1g + dx(1)
    xhe1g = xhe1g + dx(2)
    xhe2g = xhe2g + dx(3)
 enddo

 xh1 = xh1g * n / nh
 xhe1 = xhe1g * n / nhe
 xhe2 = xhe2g * n / nhe
 xh0 = ((nh/n) - xh1g) * n / nh
 xhe0 = ((nhe/n) - xhe1g - xhe2g) * n / nhe
end subroutine ionisation_fraction



subroutine minv (M, M_inv)

 implicit none

 real, dimension(3,3), intent(in)  :: M
 real, dimension(3,3), intent(out) :: M_inv

 real :: det
 real, dimension(3,3) :: cofactor


 det =   M(1,1)*M(2,2)*M(3,3)  &
       - M(1,1)*M(2,3)*M(3,2)  &
       - M(1,2)*M(2,1)*M(3,3)  &
       + M(1,2)*M(2,3)*M(3,1)  &
       + M(1,3)*M(2,1)*M(3,2)  &
       - M(1,3)*M(2,2)*M(3,1)

 cofactor(1,1) = +(M(2,2)*M(3,3)-M(2,3)*M(3,2))
 cofactor(1,2) = -(M(2,1)*M(3,3)-M(2,3)*M(3,1))
 cofactor(1,3) = +(M(2,1)*M(3,2)-M(2,2)*M(3,1))
 cofactor(2,1) = -(M(1,2)*M(3,3)-M(1,3)*M(3,2))
 cofactor(2,2) = +(M(1,1)*M(3,3)-M(1,3)*M(3,1))
 cofactor(2,3) = -(M(1,1)*M(3,2)-M(1,2)*M(3,1))
 cofactor(3,1) = +(M(1,2)*M(2,3)-M(1,3)*M(2,2))
 cofactor(3,2) = -(M(1,1)*M(2,3)-M(1,3)*M(2,1))
 cofactor(3,3) = +(M(1,1)*M(2,2)-M(1,2)*M(2,1))

 M_inv = transpose(cofactor) / det

 return

end subroutine minv

end module lightcurve
