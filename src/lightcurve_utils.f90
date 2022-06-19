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
!  Copyright (C) 2021-2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module lightcurve_utils
 use params, only:doub_prec
 implicit none

 public :: get_temp_from_u,ionisation_fraction,ionisation_fraction_Honly,get_opacity

 private

contains
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
! routine to return simple gas opacities
! from electron scattering and free-free absorption
!
! INPUT:
!    rho - density [g/cm^3]
!    T -  temperature [K]
! OUTPUT:
!    kappa - cm^2/g
!---------------------------------------------------------

real elemental function get_opacity(rho,T,X,Y) result(kappa)
 real(doub_prec), intent(in) :: rho,T,X,Y
 real(doub_prec) :: kappa_ff,kappa_es,kappa_H,kappa_mol,kappa_abs
 !real(doub_prec) :: Z
 real(doub_prec) :: xfrac,ne,xh0,xh1,xhe0,xhe1,xhe2
 real(doub_prec), parameter :: sigma_e = 6.652e-25 ! Thomson cross section

 call ionisation_fraction_Honly(rho,T,xfrac,ne)
 !call ionisation_fraction(rho,T,X,Y,xh0,xh1,xhe0,xhe1,xhe2,ne)

 !Z = max(1. - X - Y,0.)  ! metallicity

 ! free-free emission (Kramer's law)
! kappa_ff = 0.64e23*rho*T**(-3.5)
 !kappa_ff = 4.e25*(1. + X)*(Z + 0.001)*rho*T**(-3.5)

 ! opacity due to negative Hydrogen
 !kappa_H = 1.1e-25*sqrt(Z)*sqrt(rho)*T**7.7
 !kappa_H = 1.1e-25*sqrt(Z)*sqrt(rho)*T**7.7

 ! opacity due to molecules
 !kappa_mol = 0.1*Z

 ! electron scattering
! if (T > 7000.) then
!    kappa_es = 0.2*(1. + X)
! else
!    kappa_es = 0.
 !endif
 kappa_es = sigma_e*ne/rho !*0.5*(1. + X) ! 0.2*(1 + X) if fully ionised

 ! total opacity, valid between T = 1.5 x 10^3 and 10^9 K
 !kappa = kappa_mol + 1./(1./kappa_H + 1./(kappa_es + kappa_ff))

 !kappa = 1./(1./kappa_H + 1./(kappa_es + kappa_ff))
 kappa = kappa_es

 ! absorption opacity
 !kappa_abs = kappa_ff + kappa_H !+ kappa_ff
 !kappa = sqrt((kappa_ff + kappa_es + kappa_H)*kappa_abs)
 !kappa = sqrt(kappa*(kappa_abs))

end function get_opacity

!----------------------------------------------------------------
!+
!  Solves three Saha equations simultaneously to return ion
!  fractions of hydrogen and helium. Assumes inputs in cgs units
!+
!----------------------------------------------------------------
pure subroutine ionisation_fraction(dens,temp,X,Y,xh0,xh1,xhe0,xhe1,xhe2,ne)
 use vectorutils, only:matrixinvert3D
 use physcon,     only:pi,kboltz,hplanck,mh
 real, intent(in) :: dens,temp,X,Y
 real, intent(out):: xh0,xh1,xhe0,xhe1,xhe2,ne
 real             :: n,nh,nhe,A,B,C,const,xh1g,xhe1g,xhe2g,f,g,h
 real, dimension(3,3) :: M,M_inv
 real, dimension(3) :: dx
 integer          :: i,ierr
 real, parameter  :: twopi=2.*pi,eV=1.60219d-12,mass_electron_cgs=9.10938291d-28,&
                     chih0=13.6,chihe0=24.6,chihe1=54.4

 nh = X * dens / mh
 nhe = Y * dens / (4. * mh)
 n = nh + nhe

 const = (sqrt(twopi * mass_electron_cgs * kboltz) / hplanck)**3 / n

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

    call matrixinvert3D(M,M_inv,ierr)
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

 ne = xh1*nh + xhe1*nhe + xhe2*nhe

end subroutine ionisation_fraction

!----------------------------------------------------------------
!+
!  Same as above but for H-only, can be solved analytically
!  without any matrix inversion or iterations.
!  Assumes inputs in cgs units.
!+
!----------------------------------------------------------------
pure subroutine ionisation_fraction_Honly(dens,temp,xh1,ne)
 use physcon,     only:pi,kboltz,hplanck,mh
 real(doub_prec), intent(in) :: dens,temp
 real(doub_prec), intent(out):: xh1,ne
 real(doub_prec)             :: n,nh,nhe,A,const
 real(doub_prec), parameter  :: twopi=2.*pi,eV=1.60219d-12,mass_electron_cgs=9.10938291d-28,&
                                chih0=13.6

 nh = dens / mh
 n = nh

 const = (sqrt(twopi * mass_electron_cgs * kboltz) / hplanck)**3 / n
 A = const * temp**(1.5) * exp(-chih0 * eV / (kboltz * temp))

 ! solve quadratic equation x^2/(1-x) = A, i.e. x^2 + Ax - A = 0
 xh1 = 0.5*(-A + sqrt(A**2 + 4.*A))

 ! electron density
 ne = xh1*n

end subroutine ionisation_fraction_Honly

end module lightcurve_utils
