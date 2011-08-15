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
!  Copyright (C) 2005-2011 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------
module cubic
 implicit none

contains
!-------------------------------------------------------------
! this subroutine finds the real solutions to
! a cubic equation of the form
!
! a*x^3 + b*x^2 + c*x + d
!
! formulae taken from:
! Woan, The Cambridge Handbook of Physics Formulas, 2000, p51
!
! input  : a,b,c,d : coefficients of cubic polynomial
! output : x(3)    : array containing up to 3 real solutions
!                    => x(1:nreal) non-zero, rest set to zero
!          nreal   : number of real solutions
!
! Daniel Price, 22/2/07
! dprice@astro.ex.ac.uk
!-------------------------------------------------------------
subroutine cubicsolve(a,b,c,d,x,nreal,check)
 implicit none
 real, intent(in) :: a,b,c,d
 real, intent(out), dimension(3) :: x
 integer, intent(out) :: nreal
 logical, intent(in), optional :: check
 real :: p,q,det,sqrtdet
 real :: a2,b2,u,v,y1,y2,y3,term,phi
 real, parameter :: eps = 1000.*epsilon(0.)
 real, parameter :: pi = 3.14159265358979323846
 integer :: i

 x = 0.
!
!--handle all trivial cases (quadratic, linear, all zero)
!
 if (abs(a).lt.eps) then
    det = c**2 - 4.*b*d
    if (det.lt.0.) then ! no solutions to quadratic
       nreal = 0
    else
       if (abs(b).lt.eps) then
       !--no solutions if a = 0, b = 0, c = 0
          if (abs(c).lt.eps) then
             nreal = 0
          else
       !--solve linear equation if a = 0, b = 0
             nreal = 1
             x(1) = -d/c
          endif
       else
       !--solve quadratic for a = 0
          nreal = 2
          sqrtdet = sqrt(det)
          x(1) = 0.5*(-c + sqrtdet)/b
          x(2) = 0.5*(-c - sqrtdet)/b
       endif
    endif
 else
!
!--cubic solution
!
    a2 = a**2
    b2 = b**2
    p = (c/a - b2/(3.*a2))
    q = (2.*b**3/(27.*a2*a) - b*c/(3.*a2) + d/a)
    det = (p**3)/27. + 0.25*q**2
!
!--determine number of solutions
!
    if (det.lt.0.) then
    !--3 distinct real roots
       nreal = 3
       term = sqrt(abs(p)/3.)
       phi = ACOS(-0.5*q*term**(-3))

       !--these are the solutions to the reduced cubic
       !  y^3 + py + q = 0
       y1 = 2.*term*COS(phi/3.)
       y2 = -2.*term*COS((phi + pi)/3.)
       y3 = -2.*term*COS((phi - pi)/3.)
    else
    !--1 real, 2 complex roots
       term = -0.5*q + sqrt(det)
       !--must take cube root of positive quantity, then give sign later
       !  (otherwise gives NaNs)
       u = (abs(term))**(1/3.)*SIGN(1.0,term)
       term = -0.5*q - sqrt(det)
       v = (abs(term))**(1/3.)*SIGN(1.0,term)
       nreal = 1
       y1 = u + v
    !--if det=0, 3 real roots, but at least 2 equal, so max of 2 unique roots)
       if (abs(det).lt.tiny(det)) then
          nreal = 2
          y2 = -(u + v)/2.
       endif
       y3 = 0.
    endif
    !--return solutions to original cubic, not reduced cubic
    term = b/(3.*a)
    if (nreal.ge.1) x(1) = y1 - term
    if (nreal.ge.2) x(2) = y2 - term
    if (nreal.ge.3) x(3) = y3 - term

 endif

 if (present(check)) then
    if (check) then
       !--verify the cubic solution
       print*,'verifying: ',a,'x^3 + ',b,'x^2 + ',c,'x + ',d
       do i=1,nreal
          term = a*x(i)**3 + b*x(i)**2 + c*x(i) + d
          if (abs(term).lt.eps) then
             print*,'root ',i,':',x(i),'f=',term,': OK'
          else
             print*,'root ',i,':',x(i),'f=',term,': FAILED',eps
          endif
       enddo
    endif
 endif
 return

end subroutine cubicsolve

!-------------------------------------------------------------
! this subroutine returns both the real and complex
! solutions to a cubic equation of the form
!
! x^3 + b*x^2 + c*x + d
!
! input  : b,c,d : coefficients of cubic polynomial
! output : x(3)  : array of 3 COMPLEX solutions
!          nreal : number of real solutions
!
! The form of the equation above means that we
! do not need to handle trivial cases (quadratic, etc.)
! and that there will always be 3 solutions.
!
! Daniel Price, daniel.price@monash.edu 21/01/2011
!
!-------------------------------------------------------------
subroutine cubicsolve_complex(b,c,d,x,nreal,check)
 implicit none
 real,    intent(in) :: b,c,d
 complex, intent(out), dimension(3) :: x
 integer, intent(out), optional :: nreal
 logical, intent(in), optional :: check
 double precision :: p,q,q2,xi
 double precision :: b2,term,termA,det,phi
 real :: termr,termi
 double precision :: fx,dfx
 real, parameter :: eps = 1000.*epsilon(0.)
 double precision, parameter :: pi = 3.14159265358979323846d0
 integer :: i,j,nroots

 x = (0.,0.)
!
!--preliminaries
!
 b2 = b*b
 p = (c - b2/3.)
 q = (2.*b2*b - 9.*b*c + 27.*d)/27.
 q2 = q*q
 det = (p*p*p)/27. + 0.25*q2
 if (det < 0) then
 !--3 distinct real roots
    nroots = 3
    term = sqrt(abs(p)/3.)
    phi = ACOS(-0.5*q*term**(-3))

    !--these are the solutions to the reduced cubic
    !  y^3 + py + q = 0
    x(1) = real(2.d0*term*COS(phi/3.d0))
    x(2) = real(-2.d0*term*COS((phi + pi)/3.d0))
    x(3) = real(-2.d0*term*COS((phi - pi)/3.d0))
 else
 !--1 real, two complex
    nroots = 1
    if (abs(det).lt.tiny(det)) nroots = 2
    term = -0.5*q + sqrt(det)
    termA = (abs(term))**(1.d0/3.d0)*SIGN(1.0d0,term)

    x(1) = real(termA - p/(3.*termA))
    termr = real(-0.5*termA + p/(6.*termA)) ! convert from double prec.
    termi = real(0.5*sqrt(3.d0)*(termA + p/(3.d0*termA)))
    x(2) = cmplx(termr,termi)
    x(3) = cmplx(termr,-termi)
 endif

 !--return solutions to original cubic, not reduced cubic
 x(:) = x(:) - b/3.

 !--if determinant is small, take a couple of Newton-Raphson iterations
 !  to beat down the error
 if (abs(det).lt.eps) then
    do i=1,nroots
       xi = dble(x(i))
       do j=1,3
          fx = xi*(xi*(xi + b) + c) + d
          dfx = xi*(3.d0*xi + 2.d0*b) + c
          if (abs(dfx).gt.0.) xi = xi - fx/dfx
       enddo
       x(i) = real(xi)
    enddo
 endif

 if (present(nreal)) nreal = nroots

 !--the following lines can be used for debugging
 if (present(check)) then
    if (check) then
       !--verify the cubic solution
       print*,'verifying: x^3 + ',b,'x^2 + ',c,'x + ',d
       do i=1,3
          term = real(x(i)**3 + b*x(i)**2 + c*x(i) + d)
          if (abs(term).lt.eps) then
             print*,'root ',i,':',x(i),'f=',term,': OK'
          else
             print*,'root ',i,':',x(i),'f=',term,': FAILED',eps
          endif
       enddo
    endif
 endif

 return
end subroutine cubicsolve_complex

end module cubic
