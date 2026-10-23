!-----------------------------------------------------------------------
!
! unit tests of the nonlinsolve module: linear solve by Gaussian
! elimination and Newton-Raphson for systems of nonlinear equations
!
!-----------------------------------------------------------------------
module testnonlin
 use nonlinsolve, only:solve_linear,newton_nd
 implicit none
 integer, parameter, private :: dp = kind(0.d0)

contains

!-----------------------------------------------------------------------
! check a pass/fail condition and report
!-----------------------------------------------------------------------
subroutine checkval(label,val,ref,tol,nfail)
 character(len=*), intent(in) :: label
 real(dp), intent(in) :: val,ref,tol
 integer, intent(inout) :: nfail

 if (abs(val-ref) <= tol*max(1._dp,abs(ref))) then
    print "(a,a,es14.6)",' PASS: ',label,val
 else
    print "(a,a,2es14.6)",' FAIL: ',label,val,ref
    nfail = nfail + 1
 endif

end subroutine checkval

!-----------------------------------------------------------------------
! linear system needing pivoting (zero on the leading diagonal)
!-----------------------------------------------------------------------
subroutine test_linear(nfail)
 integer, intent(inout) :: nfail
 real(dp) :: a(3,3),b(3)
 integer :: ierr

 a = reshape((/0._dp,2._dp,1._dp, 1._dp,1._dp,3._dp, 2._dp,1._dp,1._dp/),(/3,3/))
 ! solution x = (1,2,3) => b = a x
 b = matmul(a,(/1._dp,2._dp,3._dp/))
 call solve_linear(a,b,3,ierr)
 if (ierr /= 0) nfail = nfail + 1
 call checkval('linear x(1)',b(1),1._dp,1.e-12_dp,nfail)
 call checkval('linear x(2)',b(2),2._dp,1.e-12_dp,nfail)
 call checkval('linear x(3)',b(3),3._dp,1.e-12_dp,nfail)

end subroutine test_linear

!-----------------------------------------------------------------------
! intersection of the circle x^2 + y^2 = 4 with the line x = y
!-----------------------------------------------------------------------
subroutine test_newton_circle(nfail)
 integer, intent(inout) :: nfail
 real(dp) :: x(2)
 integer :: its,ierr

 x = (/1._dp,0.5_dp/)
 call newton_nd(circle,x,2,1.e-12_dp,50,its,ierr)
 if (ierr /= 0) then
    print*,' FAIL: newton circle ierr = ',ierr
    nfail = nfail + 1
 endif
 call checkval('newton circle x',x(1),sqrt(2._dp),1.e-10_dp,nfail)
 call checkval('newton circle y',x(2),sqrt(2._dp),1.e-10_dp,nfail)

end subroutine test_newton_circle

!-----------------------------------------------------------------------
! a poor initial guess should still converge with step damping
!-----------------------------------------------------------------------
subroutine test_newton_badguess(nfail)
 integer, intent(inout) :: nfail
 real(dp) :: x(2)
 integer :: its,ierr

 x = (/30._dp,-0.1_dp/)
 call newton_nd(circle,x,2,1.e-12_dp,100,its,ierr)
 if (ierr /= 0) then
    print*,' FAIL: newton bad guess ierr = ',ierr
    nfail = nfail + 1
 endif
 call checkval('newton bad guess x',abs(x(1)),sqrt(2._dp),1.e-10_dp,nfail)

end subroutine test_newton_badguess

!-----------------------------------------------------------------------
! f = x^2 + 1 has no root and from x = 0 no damped Newton step reduces
! the residual. Newton should then fail, leaving x unchanged
!-----------------------------------------------------------------------
subroutine test_newton_nodescent(nfail)
 integer, intent(inout) :: nfail
 real(dp) :: x(1)
 integer :: its,ierr

 x = 0._dp
 call newton_nd(noroot,x,1,1.e-12_dp,50,its,ierr)
 if (ierr == 0) then
    print*,' FAIL: newton with no descent returned ierr = 0'
    nfail = nfail + 1
 endif
 call checkval('newton with no descent leaves x unchanged',x(1),0._dp,1.e-14_dp,nfail)

end subroutine test_newton_nodescent

!-----------------------------------------------------------------------
! from x = 0.3 damped Newton on f = x^2 + 1 takes ever smaller steps
! towards the minimum at x = 0; this must not be reported as convergence
!-----------------------------------------------------------------------
subroutine test_newton_noroot(nfail)
 integer, intent(inout) :: nfail
 real(dp) :: x(1)
 integer :: its,ierr

 x = 0.3_dp
 call newton_nd(noroot,x,1,1.e-3_dp,200,its,ierr)
 if (ierr == 0) then
    print*,' FAIL: newton with no root returned ierr = 0 at x = ',x
    nfail = nfail + 1
 else
    print*,' PASS: newton with no root returned ierr = ',ierr
 endif

end subroutine test_newton_noroot

subroutine noroot(x,f,n)
 integer, intent(in) :: n
 real(dp), intent(in) :: x(n)
 real(dp), intent(out) :: f(n)

 f = x**2 + 1._dp

end subroutine noroot

subroutine circle(x,f,n)
 integer, intent(in) :: n
 real(dp), intent(in) :: x(n)
 real(dp), intent(out) :: f(n)

 f(1) = x(1)**2 + x(2)**2 - 4._dp
 f(2) = x(1) - x(2)

end subroutine circle

end module testnonlin

program test_nonlinsolve
 use testnonlin, only:test_linear,test_newton_circle,test_newton_badguess,test_newton_nodescent,&
                      test_newton_noroot
 implicit none
 integer :: nfail

 nfail = 0
 call test_linear(nfail)
 call test_newton_circle(nfail)
 call test_newton_badguess(nfail)
 call test_newton_nodescent(nfail)
 call test_newton_noroot(nfail)

 if (nfail == 0) then
    print "(a)",' ALL NONLINSOLVE TESTS PASSED'
 else
    print "(a,i3,a)",' FAILED ',nfail,' NONLINSOLVE TESTS'
    stop 1
 endif

end program test_nonlinsolve
