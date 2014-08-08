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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------------------
! Plots Lynden-Bell & Pringle (1978) solution for viscous
! ring spreading in an accretion disk.
!
! Input: radius in disk
! Output: surface density
!
! D. Price, UofExeter 25.2.08
!----------------------------------------------------------------------
module ringspread
 implicit none
 public :: exact_ringspread, ringspreadfunc
 private

contains

subroutine exact_ringspread(iplot,time,Mdisk,Rdisk,viscnu,xplot,yplot,ierr)
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: time,Mdisk,Rdisk,viscnu
 real, intent(in), dimension(:) :: xplot
 real, intent(out), dimension(size(xplot)) :: yplot
 integer, intent(out) :: ierr
 real, parameter :: pi = 3.1415926536
 integer :: i
 real :: R2,sigma,tvisc
 double precision :: tau,x
!
! check for errors in input parameters
!
 ierr = 0
 if (Mdisk.le.0.) then
    print*,'error: mass <= 0 in exact_ringspread'
    ierr = 2
    return
 elseif (Rdisk.le.0.) then
    print*,'error: rdisk < 0 in exact_ringspread'
    ierr = 3
    return
 elseif (viscnu.le.tiny(viscnu)) then
    print*,'error: viscosity <= 0 in ringspreading solution'
    ierr = 4
    return
 endif

 R2 = Rdisk*Rdisk
 tvisc = R2/(12.*viscnu)
 tau = time/tvisc

 print "(a,1pe9.2,a,1pe9.2,a,0pf6.2,a,f6.2)", &
 ' Plotting ring spreading solution: tau = ',tau,' nu = ',viscnu,' R0 = ',Rdisk,' M = ',Mdisk

 do i=1,size(xplot)
    x = xplot(i)/Rdisk
    sigma = Mdisk/real((pi*R2))*ringspreadfunc(x,tau)
    !print*,'x = ',xplot(i),Rdisk,tau,sigma

    select case(iplot)
    case(1)
    !--density
       yplot(i) = sigma
    case default
    !--pressure
       yplot(i) = 0.
    end select
 enddo

 return
end subroutine exact_ringspread

!----------------------------------------------------------------------
! evaluates the surface density as a function of x and tau
!----------------------------------------------------------------------
double precision function ringspreadfunc(x,tau)
 implicit none
 double precision, intent(in) :: x, tau
 double precision :: xfunc,besfunc,dummy,term

 if (tau.le.epsilon(tau) .or. x.le.tiny(x)) then
    ringspreadfunc = 0.
 else
    xfunc = 2.*x/tau
    term = exp(-(1.+x*x)/tau)

    !--prevent blowups at t=0: no point evaluating
    !  the Bessel function if the exp term is zero.
    if (term.gt.tiny(term)) then
       call bessik(xfunc,0.25d0,besfunc,dummy,dummy,dummy)
    else
       besfunc = 0.
    endif
    ringspreadfunc = 1./(tau*x**0.25)*term*besfunc
    !print*,'ringspreadfunc = ',x,xfunc,-(1.+x*x)/tau,exp(-(1.+x*x)/tau),ringspreadfunc,xfunc,tau,besfunc
 endif

 return
end function ringspreadfunc


!----------------------------------------------------------------------
! remainder of this file are routines for evaluating the modified
! Bessel function in the above solution
!----------------------------------------------------------------------

subroutine bessik(x,xnu,ri,rk,rip,rkp)
 implicit none
 integer, parameter :: maxit = 10000
 double precision, intent(in) :: x, xnu
 double precision, intent(out) :: ri,rip,rk,rkp
 double precision, parameter :: eps = 1.e-10, fpmin = 1.e-30, &
   xmin = 2., pi = 3.141592653589793d0

! Returns the modified Bessel functions ri = I\nu, rk = K\nu and their
! derivatives rip = I'\nu and rkp = K'\nu, for positive x and for
! xn = \nu .ge. 0. The relative accuracy is within one or two significant
! digits of eps.
!
! All internal arithmetic in double precision
!
! This routine written by Daniel Price, UofExeter 25/2/08
! Adapted from Press et. al. (1992) Numerical Recipes in FORTRAN 77
!
 integer:: i,l,nl
 double precision :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2, &
  ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl, &
  ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2

 if (x <= 0. .or. xnu < 0.) then
    print*,' bad arguments in bessik ',x,xnu
!    return
 endif

 nl = int(xnu+0.5d0) ! nl is the number of downward recurrences of the I's
 xmu = xnu - nl      ! and upward recurrences of K's. xmu lies between
 xmu2 = xmu*xmu      ! -1/2 and 1/2
 xi = 1.d0/x
 xi2 = 2.d0*xi
 h = xnu*xi
 if (h.lt.fpmin) h = fpmin
 b = xi2*xnu
 d = 0.d0
 c = h
 do i=1,maxit
    b = b + xi2
    d = 1.d0/(b + d)
    c = b + 1.d0/c
    del = c*d
    h = del*h
    if (abs(del-1.d0).lt.eps) goto 1
 enddo
 print*,'x too large in bessik; try asymptotic expansion'
1 continue
 ril = fpmin
 ripl = h*ril
 ril1 = ril
 rip1 = ripl
 fact = xnu*xi
 do l=nl,1,-1
    ritemp = fact*ril + ripl
    fact = fact - xi
    ripl = fact*ritemp + ril
    ril = ritemp
 enddo
 f = ripl/ril
 if (x < xmin) then
    x2 = 0.5d0*x
    pimu = pi*xmu
    if (abs(pimu).lt.eps) then
       fact = 1.d0
    else
       fact = pimu/sin(pimu)
    endif
    d = -log(x2)
    e = xmu*d
    if (abs(e).lt.eps) then
       fact2 = 1.d0
    else
       fact2 = sinh(e)/e
    endif
    ! Chebyshev evaluation of Gamma1, Gamma2
    call beschb(xmu,gam1,gam2,gampl,gammi)
    ff = fact*(gam1*cosh(e) + gam2*fact2*d)  ! f0
    sum = ff
    e = exp(e)
    p = 0.5d0*e/gampl
    q = 0.5d0/(e*gammi)
    c = 1.d0
    d = x2*x2
    sum1 = p
    do i=1,maxit
       ff = (i*ff + p + q)/(i*i - xmu2)
       c = c*d/i
       p = p/(i - xmu)
       q = q/(i + xmu)
       del = c*ff
       sum = sum + del
       del1 = c*(p - i*ff)
       sum1 = sum1 + del1
       if (abs(del).lt.abs(sum)*eps) goto 2
    enddo
    print*,' bessk series failed to converge'
2   continue
    rkmu = sum
    rk1 = sum1*xi2
 else
    b = 2.d0*(1.d0 + x)
    d = 1.d0/b
    delh = d
    h = delh
    q1 = 0.d0
    q2 = 1.d0
    a1 = 0.25d0 - xmu2
    c = a1
    q = c
    a = -a1
    s = 1.d0 + q*delh
    do i=2,maxit
       a = a - 2*(i-1)
       c = -a*c/i
       qnew = (q1 - b*q2)/a
       q1 = q2
       q2 = qnew
       q = q + c*qnew
       b = b + 2.d0
       d = 1.d0/(b + a*d)
       delh = (b*d - 1.d0)*delh
       h = h + delh
       dels = q*delh
       s = s + dels
       if (abs(dels/s).lt.eps) goto 3
    enddo
    print*,'bessik: failure to converge in cf2'
3   continue
    h = a1*h
    rkmu = sqrt(pi/(2.d0*x))*exp(-x)/s
    rk1 = rkmu*(xmu + x + 0.5d0 - h)*xi
 endif
 rkmup = xmu*xi*rkmu - rk1
 rimu = xi/(f*rkmu - rkmup)
 ri = (rimu*ril1)/ril
 rip = (rimu*rip1)/ril
 do i=1,nl
    rktemp = (xmu + i)*xi2*rk1 + rkmu
    rkmu = rk1
    rk1 = rktemp
 enddo
 rk = rkmu
 rkp = xnu*xi*rkmu - rk1

 return
end subroutine bessik

subroutine beschb(x,gam1,gam2,gampl,gammi)
 implicit none
 integer, parameter :: nuse1=7,nuse2=8
 double precision, intent(in) :: x
 double precision, intent(out) :: gam1,gam2,gammi,gampl
!
! Evaluates Gamma_1 and Gamma_2 by Chebyshev expansion for |x| < 1/2.
! Also returns 1/Gamma(1 + x) and 1/Gamma(1-x).
!
! In double precision, set NUSE1 = 7, NUSE2 = 8.
! In single precision, set NUSE1 = 2, NUSE2 = 5.
!
! This routine written by Daniel Price, UofExeter 25/2/08
! Adapted from Press et. al. (1992) Numerical Recipes in FORTRAN 77
!
 double precision :: xx,c1(7),c2(8)
 save c1,c2
 data c1/-1.142022680371168d0, 6.5165112670737d-3, &
         3.087090173086d-4, -3.4706269649d-6,6.9437664d-9, &
         3.67795d-11,-1.356d-13/
 data c2/1.843740587300905d0, -7.68528408447867d-2, &
         1.2719271366546d-3, -4.9717367042d-6, -3.31261198d-8, &
         2.423096d-10, -1.702d-13, -1.49d-15/
 xx = 8.d0*x*x - 1.d0               ! multiply x by 2 to make range be -1 to 1,
 gam1 = chebev(-1.d0,1.d0,c1,NUSE1,xx)  ! and then apply transformation for
 gam2 = chebev(-1.d0,1.d0,c2,NUSE2,xx)  ! evaluating even Chebyshev series
 gampl = gam2 - x*gam1
 gammi = gam2 + x*gam1

 return
end subroutine beschb

double precision function chebev(a,b,c,m,x)
 implicit none
 integer, intent(in) :: m
 double precision, intent(in) :: a,b,x,c(m)
!
! Chebyshev evaluation: All arguments are input. c(1:m) is an array of
! Chebyshev coefficients, the first m elements of c output from chebft
! (which must have been called with the same a and b). The Chebyshev
! polynomial \sum_{k=1}^m c_k T_{k-1}(y) - c1/2 is evaluated at a
! point y = [ x - (b + a)/2 ] / [(b - a)/2], and the result is returned
! as the function value.
!
! This routine written by Daniel Price, UofExeter 25/2/08
! Adapted from Press et. al. (1992) Numerical Recipes in FORTRAN 77
!
 integer :: j
 double precision :: d, dd, sv, y, y2

 if ((x-a)*(x-b).gt.0.) then
    print*,'error: x not in range in chebev'
    chebev = 0.
    return
 endif
 d = 0.
 dd = 0.
 y = (2.*x - a - b)/(b - a)  ! change of variable
 y2 = 2.*y
 do j=m,2,-1  ! Clenshaw's recurrence
    sv = d
    d = y2*d - dd + c(j)
    dd = sv
 enddo
 chebev = y*d - dd + 0.5*c(1) ! last step is different

end function chebev

end module ringspread
