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
!  Copyright (C) 2026 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Solvers for small systems of linear and nonlinear equations
!
! solve_linear : Gaussian elimination with partial pivoting
! newton_nd    : damped Newton-Raphson iteration with a
!                finite difference Jacobian
!
! This module has no dependencies, so it can be used standalone
!
!-----------------------------------------------------------------------
module nonlinsolve
 implicit none
 integer, parameter, private :: dp = kind(0.d0)
 public :: solve_linear,newton_nd,resid_func

 abstract interface
  subroutine resid_func(x,f,n)
   import :: dp
   integer, intent(in) :: n
   real(dp), intent(in) :: x(n)
   real(dp), intent(out) :: f(n)
  end subroutine resid_func
 end interface

 private

contains

!-----------------------------------------------------------------------
!
! solve the linear system a x = b by Gaussian elimination with
! partial pivoting. On output b contains the solution x and
! a is overwritten. ierr = 1 if the matrix is singular
!
!-----------------------------------------------------------------------
subroutine solve_linear(a,b,n,ierr)
 integer,  intent(in)    :: n
 real(dp), intent(inout) :: a(n,n),b(n)
 integer,  intent(out)   :: ierr
 integer  :: i,k,ipiv
 real(dp) :: fac,rowtmp(n),btmp

 ierr = 0
 do k=1,n
    !--choose the largest available pivot in column k
    ipiv = k - 1 + maxloc(abs(a(k:n,k)),dim=1)
    if (abs(a(ipiv,k)) <= tiny(1._dp)) then
       ierr = 1
       return
    endif
    if (ipiv /= k) then
       rowtmp(:)  = a(k,:)
       a(k,:)     = a(ipiv,:)
       a(ipiv,:)  = rowtmp(:)
       btmp       = b(k)
       b(k)       = b(ipiv)
       b(ipiv)    = btmp
    endif
    !--eliminate entries below the pivot
    do i=k+1,n
       fac    = a(i,k)/a(k,k)
       a(i,k:n) = a(i,k:n) - fac*a(k,k:n)
       b(i)   = b(i) - fac*b(k)
    enddo
 enddo

 !--back substitution
 do i=n,1,-1
    b(i) = (b(i) - dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
 enddo

end subroutine solve_linear

!-----------------------------------------------------------------------
!
! solve f(x) = 0 for n unknowns using Newton-Raphson iterations
! with a finite difference Jacobian. The step is halved until the
! residual decreases, which helps convergence from poor guesses
!
! input  : func   : subroutine returning residuals f(x)
!          x      : initial guess
!          tol    : convergence tolerance on the step and residual
!          maxits : maximum number of iterations
! output : x      : solution
!          its    : number of iterations taken
!          ierr   : 0 on success, 1 if singular, 2 if not converged,
!                   3 if no step reduces the residual
!
!-----------------------------------------------------------------------
subroutine newton_nd(func,x,n,tol,maxits,its,ierr)
 procedure(resid_func)   :: func
 integer,  intent(in)    :: n,maxits
 real(dp), intent(inout) :: x(n)
 real(dp), intent(in)    :: tol
 integer,  intent(out)   :: its,ierr
 integer, parameter :: maxhalve = 30
 real(dp) :: f(n),fnew(n),xnew(n),dx(n),jac(n,n),fnorm,lambda
 integer  :: ihalve
 logical  :: converged

 ierr = 0
 its  = 0
 call func(x,f,n)
 fnorm = norm2(f)
 converged = (fnorm < tol)

 do while (.not.converged .and. its < maxits)
    its = its + 1
    call get_jacobian(func,x,f,n,jac)
    dx = -f
    call solve_linear(jac,dx,n,ierr)
    if (ierr /= 0) return

    !--damped step: halve the step until the residual decreases,
    !  giving up (with x unchanged) if no step reduces it
    lambda = 1._dp
    do ihalve=1,maxhalve
       xnew = x + lambda*dx
       call func(xnew,fnew,n)
       if (norm2(fnew) < fnorm) exit
       lambda = 0.5_dp*lambda
    enddo
    if (ihalve > maxhalve) then
       ierr = 3
       return
    endif

    !--a small step only indicates convergence if it was a full Newton step
    converged = (lambda >= 1._dp .and. maxval(abs(dx)/max(1._dp,abs(x))) < tol)
    x = xnew
    f = fnew
    fnorm = norm2(f)
    if (fnorm < tol) converged = .true.
 enddo

 if (.not.converged) ierr = 2

end subroutine newton_nd

!-----------------------------------------------------------------------
!
! Jacobian of func at x by one-sided finite differences,
! given the residual f = func(x)
!
!-----------------------------------------------------------------------
subroutine get_jacobian(func,x,f,n,jac)
 procedure(resid_func)  :: func
 integer,  intent(in)   :: n
 real(dp), intent(in)   :: x(n),f(n)
 real(dp), intent(out)  :: jac(n,n)
 real(dp) :: xp(n),fp(n),h
 integer  :: j

 do j=1,n
    xp = x
    h  = sqrt(epsilon(1._dp))*max(1._dp,abs(x(j)))
    xp(j) = x(j) + h
    call func(xp,fp,n)
    jac(:,j) = (fp - f)/h
 enddo

end subroutine get_jacobian

end module nonlinsolve
