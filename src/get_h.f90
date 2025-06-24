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
!  Copyright (C) 2005-2025 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module get_h
!
! this module contains routines to compute a smoothing length 
! for post-processing of simulations where h is not output
! in the file (e.g. dark matter, stars etc)
!
 implicit none
 public :: compute_h_and_density

 real, parameter :: tolh = 1.e-4
 integer, parameter :: its_max = 50

 private

contains

!
! public-facing routine to compute h and density
! from a set of positions. No guess for h is needed
!
subroutine compute_h_and_density(ndim,n,hfac,x,h,dens)
 integer, intent(in)  :: ndim,n
 real,    intent(in)  :: hfac
 real,    intent(in)  :: x(ndim,n)
 real,    intent(out) :: h(n),dens(n)
 real :: hguess

 print*,' >>> Computing smoothing lengths for ',n,' particles...'

 hguess = abs(maxval(x(1,1:n)) - minval(x(1,1:n)))/n**(1./3.)
 h = hguess
 print*, ' guessing h = ',hguess

 ! call guess_h(ndim,x,h)
 call densityiterate(ndim,n,hfac,x,h,dens)

end subroutine compute_h_and_density

!
! routine to perform the simultaneous
! (number) density and smoothing length calculation
!
subroutine densityiterate(ndim,n,hfac,x,h,dens)
 use kernels, only:wfunc,dwfunc,radkernel2,select_kernel_by_name,&
                   cnormk1D,cnormk2D,cnormk3D
 use timing,  only:wall_time,print_time
 integer, intent(in)    :: ndim,n
 real,    intent(in)    :: hfac
 real,    intent(in)    :: x(ndim,n)
 real,    intent(inout) :: h(n)
 real,    intent(out)   :: dens(n)
 real :: xi(ndim),dx(ndim)
 real :: rhoi,gradhi,hi,hi_old,hi1,hi21,hterm
 real :: q2,q,wi,dwi,dwdhi,t1,t2
 logical :: converged
 integer, parameter :: its_max = 50
 integer :: its,i,j,k

 call wall_time(t1)
 call select_kernel_by_name('cubic')
 k = 0
 !$omp parallel do schedule(dynamic) default(none)&
 !$omp shared(n,x,h,dens,ndim,radkernel2,wfunc,dwfunc,hfac,k) &
 !$omp shared(cnormk1D,cnormk2D,cnormk3D) &
 !$omp private(xi,dx,rhoi,gradhi) &
 !$omp private(hi,hi_old,hi1,hi21,hterm) &
 !$omp private(q2,q,wi,dwi,dwdhi,converged,its)
 do i=1,n
    !$omp atomic
    k = k + 1
    if (mod(k,10000)==0) then
       !$omp critical
       print*,k
       !$omp end critical
    endif
    xi     = x(1:ndim,i)
    hi     = h(i)
    hi_old = hi
    converged = .false.
    its = 0
    do while(.not.converged .and. its < its_max)
       its = its + 1
       rhoi   = 0.
       gradhi = 0.
       hi1    = 1./hi
       hi21   = hi1*hi1
       hterm  = hi1**ndim
       do j=1,n
          dx = xi - x(:,j)
          q2 = dot_product(dx,dx)*hi21
          if (q2 < radkernel2) then
             q = sqrt(q2)
             wi = wfunc(q2)*hterm
             dwi = dwfunc(q2)*hterm*hi1
             rhoi = rhoi + wi
             dwdhi = (-q*dwi - ndim*wi*hi1)
             gradhi = gradhi + dwdhi
             !print*,'j=',j,wi,hterm
          endif
       enddo
       select case(ndim)
       case(1)
          rhoi = rhoi*cnormk1D
          gradhi = gradhi*cnormk1D
       case(2)
          rhoi = rhoi*cnormk2D
          gradhi = gradhi*cnormk2D
       case default
          rhoi = rhoi*cnormk3D
          gradhi = gradhi*cnormk3D
       end select
       call update_h(hfac,hi_old,rhoi,hi,gradhi,ndim,converged)
    enddo

    ! store converged solution
    dens(i) = rhoi
    h(i) = hi
 enddo
 !$omp end parallel do
 call wall_time(t2)
 call print_time(t2-t1)

end subroutine densityiterate

!
! cleanup and iterate h once the sums over
! neighbours are done
!
subroutine update_h(hfac,h_old,rho,h,gradh,ndim,converged)
 real, intent(in)     :: hfac,h_old,rho
 real, intent(inout)  :: h,gradh
 integer, intent(in)  :: ndim
 logical, intent(out) :: converged
 real :: dhdrho,rhoh,omega,func,dfdh1,hnew

 rhoh   = (hfac/abs(h))**ndim
 dhdrho = -h/(ndim*rho)
 omega  = 1. - dhdrho*gradh
 gradh  = 1./omega

 ! Newton Raphson solver, see Price & Monaghan (2007)
 func = rhoh - rho
 if (omega > tiny(omega)) then
    dfdh1 = dhdrho/omega
 else
    dfdh1 = dhdrho/abs(omega + epsilon(omega))
 endif

 ! damped newton-raphson iteration
 hnew = h - func*dfdh1
 if (hnew > 1.2*h) then
    hnew = 1.2*h
 elseif (hnew < 0.8*h) then
    hnew = 0.8*h
 endif

 ! convergence condition
 converged = ((abs(hnew-h)/h_old) < tolh .and. omega > 0. .and. h > 0.)
 h = hnew

end subroutine update_h

end module get_h