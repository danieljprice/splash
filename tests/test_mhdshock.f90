!-----------------------------------------------------------------------
!
! unit tests of the exact MHD Riemann solver in exact_mhdshock.f90
!
! Reference values for the preset problems are the plateau states
! of the tabulated solutions previously hardwired in exact_mhdshock.f90
! (Ryu & Jones 1995; Balsara 1998; Dai & Woodward 1994; Brio & Wu 1988),
! sampled midway between the tabulated wave positions
!
!-----------------------------------------------------------------------
module testmhdshock
 use mhdshock, only:exact_mhdshock,get_mhdshock_preset
 implicit none
 integer, parameter :: maxpts = 5000
 integer, parameter :: irho = 1, ipr = 2, ivx = 3, ivy = 4, ivz = 5, iby = 6, ibz = 7
 real, parameter :: pi = 4.*atan(1.)
 real, parameter :: const = 1./sqrt(4.*pi)
 real :: xpts(maxpts),ypts(maxpts),xr(maxpts),yr(maxpts),xb(maxpts),yb(maxpts)

contains

!-----------------------------------------------------------------------
! check a value against a reference, with tolerance relative to
! max(|ref|,floor) so that small reference values are not over-tested
!-----------------------------------------------------------------------
subroutine checkval(label,val,ref,tol,nfail)
 character(len=*), intent(in) :: label
 real, intent(in) :: val,ref,tol
 integer, intent(inout) :: nfail
 real, parameter :: floor = 0.1

 if (abs(val-ref) <= tol*max(abs(ref),floor)) then
    print "(a,a,2es14.6)",' PASS: ',label,val,ref
 else
    print "(a,a,2es14.6)",' FAIL: ',label,val,ref
    nfail = nfail + 1
 endif

end subroutine checkval

!-----------------------------------------------------------------------
! linear interpolation of the plotted solution at position x
!-----------------------------------------------------------------------
real function value_at(x,xpts,ypts,npts)
 real, intent(in) :: x,xpts(:),ypts(:)
 integer, intent(in) :: npts
 integer :: i

 value_at = ypts(npts)
 do i=2,npts
    if (xpts(i) >= x) then
       if (xpts(i) > xpts(i-1)) then
          value_at = ypts(i-1) + (ypts(i)-ypts(i-1))*(x-xpts(i-1))/(xpts(i)-xpts(i-1))
       else
          value_at = ypts(i)
       endif
       return
    endif
 enddo

end function value_at

!-----------------------------------------------------------------------
! plot the solution for one of the preset initial states
!-----------------------------------------------------------------------
subroutine plot_preset(iv,ishk,time,gamma,x,y,n,ierr)
 integer, intent(in)  :: iv,ishk
 real,    intent(in)  :: time,gamma
 real,    intent(inout) :: x(:)
 real,    intent(out) :: y(size(x))
 integer, intent(out) :: n,ierr
 real :: bx,primL(7),primR(7)

 call get_mhdshock_preset(ishk,bx,primL,primR,ierr)
 call exact_mhdshock(iv,time,gamma,bx,primL,primR,-0.5,0.5,0.,x,y,n,ierr)

end subroutine plot_preset

!-----------------------------------------------------------------------
! compare the solution for a preset problem to reference plateau
! states ref(1:7,j) = (rho,p,vx,vy,vz,By,Bz) at positions xsample(j)
!-----------------------------------------------------------------------
subroutine check_preset(ishk,label,time,gamma,xsample,ref,tol,nfail,ivars)
 integer, intent(in) :: ishk
 character(len=*), intent(in) :: label
 real, intent(in) :: time,gamma,xsample(:),ref(:,:),tol
 integer, intent(inout) :: nfail
 integer, intent(in) :: ivars(:)
 integer :: npts,ierr,j,k,iv
 character(len=60) :: string

 do k=1,size(ivars)
    iv = ivars(k)
    call plot_preset(iv,ishk,time,gamma,xpts,ypts,npts,ierr)
    if (ierr /= 0) then
       print*,' FAIL: '//trim(label)//' ierr = ',ierr
       nfail = nfail + 1
       return
    endif
    do j=1,size(xsample)
       write(string,"(a,' var ',i1,' at x=',f7.4)") trim(label),iv,xsample(j)
       call checkval(trim(string),value_at(xsample(j),xpts,ypts,npts),ref(iv,j),tol,nfail)
    enddo
 enddo

end subroutine check_preset

!-----------------------------------------------------------------------
! Brio & Wu (1988), gamma = 2, with the compound wave
!-----------------------------------------------------------------------
subroutine test_briowu(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,4)

 ref(:,1) = (/0.67623,0.447,0.63721,-0.23345,0.,2.1*const,0./)
 ref(:,2) = (/0.6962,0.516,0.600,-1.584,0.,-1.9*const,0./)
 ref(:,3) = (/0.2352,0.516,0.600,-1.584,0.,-1.9*const,0./)
 ref(:,4) = (/0.117,0.0876,-0.24,-0.166,0.,-3.25*const,0./)
 call check_preset(1,'Brio-Wu',0.1,2.,(/-0.055,0.03,0.10,0.24/),ref,0.03,nfail,(/1,2,3,4,6/))

end subroutine test_briowu

!-----------------------------------------------------------------------
! fast switch-on shock and slow shock (Ryu & Jones 1995)
!-----------------------------------------------------------------------
subroutine test_switchon(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,4)

 ! two densities in the old table are corrected to be consistent with
 ! the other tabulated values: rho = p^(1/gamma) behind the fast
 ! rarefaction (table had 0.5955), and u = cA behind the switch-on
 ! shock given the tabulated vx (table had 0.2337)
 ref(:,1) = (/0.59953,0.42629,0.81237,-0.59961,0.,0.28431,0./)
 ref(:,2) = (/0.55151,0.37090,0.89416,-0.5447,0.,0.31528,0./)
 ref(:,3) = (/0.41272,0.37090,0.89416,-0.5447,0.,0.31528,0./)
 ref(:,4) = (/0.22337,0.12402,0.24722,-0.91164,0.,0.43086,0./)
 call check_preset(2,'switch-on',0.15,5./3.,(/-0.06,0.06,0.19,0.30/),ref,5.e-3,nfail,(/1,2,3,4,6/))

end subroutine test_switchon

!-----------------------------------------------------------------------
! seven discontinuities (Ryu & Jones 1995 test 2a)
!-----------------------------------------------------------------------
subroutine test_sevenjump(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,4)

 ref(:,1) = (/1.4903,1.6558,0.60588,0.11235,0.55686,1.4383,0.79907/)
 ref(:,2) = (/1.6343,1.9317,0.57538,0.047602,0.24734,1.4126,0.43772/)
 ref(:,3) = (/1.4735,1.9317,0.57538,0.047601,0.24734,1.4126,0.43772/)
 ref(:,4) = (/1.3090,1.5844,0.53432,-0.094572,-0.047286,1.5078,0.75392/)
 call check_preset(3,'7-jump',0.2,5./3.,(/-0.08,0.085,0.15,0.3/),ref,5.e-3,nfail,(/1,2,3,4,5,6,7/))

end subroutine test_sevenjump

!-----------------------------------------------------------------------
! isothermal MHD (Balsara 1998); tabulated values are approximate
!-----------------------------------------------------------------------
subroutine test_isothermal(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,5)

 ref(:,1) = (/1.515,1.515,0.65,0.13,0.57,5.2*const,2.885*const/)
 ref(:,2) = (/1.515,1.515,0.65,0.24,0.31,5.7*const,1.76*const/)
 ref(:,3) = (/1.745,1.745,0.62,0.071,0.255,5.22*const,1.62*const/)
 ref(:,4) = (/1.36,1.36,0.54,-0.215,0.165,5.96*const,1.85*const/)
 ref(:,5) = (/1.36,1.36,0.54,-0.125,-0.06,5.58*const,2.79*const/)
 call check_preset(4,'isothermal',0.2,1.,(/-0.06,0.0525,0.12,0.185,0.3/),ref,0.03,nfail,(/1,3,4,5,6,7/))

end subroutine test_isothermal

!-----------------------------------------------------------------------
! rarefactions with no normal field (Ryu & Jones 1995 test 5a)
!-----------------------------------------------------------------------
subroutine test_rarefaction(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,1)

 ref(:,1) = (/0.49653,0.31134,0.,0.,0.,0.49638,0./)
 call check_preset(5,'rarefaction',0.1,5./3.,(/0./),ref,1.e-3,nfail,(/1,2,3,6/))

end subroutine test_rarefaction

!-----------------------------------------------------------------------
! Mach 25 colliding flows (Dai & Woodward 1994)
!-----------------------------------------------------------------------
subroutine test_mach25(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,1)

 ref(:,1) = (/3.982,1806.,0.,-0.07727,-0.01932,15.95*const,3.988*const/)
 call check_preset(6,'Mach 25',0.03,5./3.,(/0./),ref,5.e-3,nfail,(/1,2,3,4,5,6,7/))

end subroutine test_mach25

!-----------------------------------------------------------------------
! fast and slow shocks (Ryu & Jones 1995 test 1a)
!-----------------------------------------------------------------------
subroutine test_rj1a(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,4)

 ref(:,1) = (/2.6797,150.98,0.72113,0.23139,0.,3.8389,0./)
 ref(:,2) = (/2.6713,150.19,0.72376,0.35684,0.,4.0380,0./)
 ref(:,3) = (/3.8508,150.19,0.72376,0.35684,0.,4.0380,0./)
 ref(:,4) = (/3.7481,143.57,0.70505,-0.38804,0.,5.4272,0./)
 call check_preset(7,'RJ95 1a',0.08,5./3.,(/-0.2,0.02,0.085,0.245/),ref,5.e-3,nfail,(/1,2,3,4,6/))

end subroutine test_rj1a

!-----------------------------------------------------------------------
! with zero magnetic field the solution should reduce to the
! hydrodynamic Riemann problem (Sod shock tube)
!-----------------------------------------------------------------------
subroutine test_hydro_limit(nfail)
 use shock, only:exact_shock
 integer, intent(inout) :: nfail
 real :: xs(6),yref(6)
 real :: primL(7),primR(7)
 integer :: npts,ierr,iv,j
 character(len=40) :: string

 primL = (/1.,1.,0.,0.,0.,0.,0./)
 primR = (/0.125,0.1,0.,0.,0.,0.,0./)
 xs = (/-0.25,-0.1,0.,0.1,0.2,0.3/)
 do iv=1,3
    call exact_mhdshock(iv,0.2,1.4,0.,primL,primR,-0.5,0.5,0.,xpts,ypts,npts,ierr)
    call exact_shock(iv,0.2,1.4,0.,primL(1),primR(1),primL(2),primR(2),0.,0.,0.,xs,yref,ierr)
    do j=1,size(xs)
       write(string,"('Sod var ',i1,' at x=',f6.3)") iv,xs(j)
       call checkval(trim(string),value_at(xs(j),xpts,ypts,npts),yref(j),1.e-4,nfail)
    enddo
 enddo

end subroutine test_hydro_limit

!-----------------------------------------------------------------------
! with Bx = 0 the quantities By/rho and p/rho^gamma are constant
! across the rarefaction fans
!-----------------------------------------------------------------------
subroutine test_invariants(nfail)
 integer, intent(inout) :: nfail
 real :: x
 integer :: nr,np,nb,ierr,j
 character(len=40) :: string
 real, parameter :: gam = 5./3.

 call plot_preset(1,5,0.1,gam,xr,yr,nr,ierr)
 call plot_preset(2,5,0.1,gam,xpts,ypts,np,ierr)
 call plot_preset(6,5,0.1,gam,xb,yb,nb,ierr)
 do j=1,5
    x = -0.25 + 0.04*j
    write(string,"('By/rho in fan at x=',f6.3)") x
    call checkval(trim(string),value_at(x,xb,yb,nb)/value_at(x,xr,yr,nr),1.,1.e-4,nfail)
    write(string,"('p/rho^gam in fan at x=',f6.3)") x
    call checkval(trim(string),value_at(x,xpts,ypts,np)/value_at(x,xr,yr,nr)**gam,1.,1.e-4,nfail)
 enddo

end subroutine test_invariants

end module testmhdshock

program test_mhdshock
 use testmhdshock
 implicit none
 integer :: nfail

 nfail = 0
 call test_hydro_limit(nfail)
 call test_invariants(nfail)
 call test_rarefaction(nfail)
 call test_rj1a(nfail)
 call test_sevenjump(nfail)
 call test_mach25(nfail)
 call test_switchon(nfail)
 call test_briowu(nfail)
 call test_isothermal(nfail)

 if (nfail == 0) then
    print "(a)",' ALL MHD SHOCK TESTS PASSED'
 else
    print "(a,i4,a)",' FAILED ',nfail,' MHD SHOCK TESTS'
    stop 1
 endif

end program test_mhdshock
