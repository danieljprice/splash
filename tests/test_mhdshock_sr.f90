!-----------------------------------------------------------------------
!
! unit tests of the exact special relativistic MHD Riemann solver
! in exact_mhdshock_sr.f90
!
! Reference values are the states in each region of the Riemann fan
! tabulated by Giacomazzo & Rezzolla (2006), J. Fluid Mech. 562, 223,
! Tables 2, 3 and 5-12, given to four significant figures. As in
! those tables, p is the total (gas + magnetic) pressure
!
!-----------------------------------------------------------------------
module testmhdshocksr
 use mhdshock_sr, only:exact_mhdshock_sr,get_mhdshock_sr_preset,solve_mhdshock_sr,&
                       total_pressure_sr
 implicit none
 integer, parameter :: maxpts = 5000
 real :: xpts(maxpts),ypts(maxpts)

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
! compare the states in each region of the Riemann fan for a preset
! problem with ref(1:7,i) = (rho,p_total,vx,vy,vz,By,Bz) in region i
!-----------------------------------------------------------------------
subroutine check_regions(ishk,label,gamma,ref,tol,nfail,tolreg)
 integer, intent(in) :: ishk
 character(len=*), intent(in) :: label
 real, intent(in) :: gamma,ref(:,:),tol
 integer, intent(inout) :: nfail
 real, intent(in), optional :: tolreg(:)
 real :: bx,primL(7),primR(7),ureg(7,8),val(7),toli
 integer :: nreg,ierr,i,k
 character(len=60) :: string

 call get_mhdshock_sr_preset(ishk,bx,primL,primR,ierr)
 call solve_mhdshock_sr(gamma,bx,primL,primR,nreg,ureg,ierr)
 if (ierr /= 0 .or. nreg /= size(ref,2)) then
    print*,' FAIL: '//trim(label)//' ierr = ',ierr,' nreg = ',nreg
    nfail = nfail + 1
    return
 endif
 do i=1,nreg
    val = ureg(:,i)
    val(2) = total_pressure_sr(ureg(:,i),bx)
    toli = tol
    if (present(tolreg)) toli = tolreg(i)
    do k=1,7
       write(string,"(a,' R',i1,' var ',i1)") trim(label),i,k
       call checkval(trim(string),val(k),ref(k,i),toli,nfail)
    enddo
 enddo

end subroutine check_regions

!-----------------------------------------------------------------------
! tests with zero normal field (GR06 Tables 2 and 3)
!-----------------------------------------------------------------------
subroutine test_tangential(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,4)

 ref(:,1) = (/1.0,230.0,0.0,0.0,0.0,20.0,0.0/)
 ref(:,2) = (/0.2410,16.11,0.8497,0.0,0.0,9.141,0.0/)
 ref(:,3) = (/0.6426,16.11,0.8497,0.0,0.0,0.0,0.0/)
 ref(:,4) = (/0.1,1.0,0.0,0.0,0.0,0.0,0.0/)
 call check_regions(2,'Komissarov ST2',4./3.,ref,1.e-3,nfail)

 ref(:,1) = (/1.0,18.19,0.1,0.3,0.4,6.0,2.0/)
 ref(:,2) = (/1.581,44.59,-0.3073,0.3082,0.2927,9.582,3.194/)
 ref(:,3) = (/5.489e-4,44.59,-0.3073,0.7488,0.5556,1.023,4.092/)
 ref(:,4) = (/0.01,5138.,0.5,0.4,0.3,5.0,20.0/)
 call check_regions(10,'generic shock tube',5./3.,ref,1.e-3,nfail)

end subroutine test_tangential

!-----------------------------------------------------------------------
! Komissarov (1999) shock tube 1 and collision (GR06 Tables 5 and 6)
!-----------------------------------------------------------------------
subroutine test_komissarov(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,8)

 ref(:,1) = (/1.0,1001.,0.,0.,0.,0.,0./)
 ref(:,2) = (/0.06984,29.27,0.9115,0.,0.,0.,0./)
 ref(:,3) = ref(:,2)
 ref(:,4) = ref(:,2)
 ref(:,5) = (/0.8846,29.27,0.9115,0.,0.,0.,0./)
 ref(:,6) = ref(:,5)
 ref(:,7) = ref(:,5)
 ref(:,8) = (/0.1,1.5,0.,0.,0.,0.,0./)
 call check_regions(1,'Komissarov ST1',4./3.,ref,1.e-3,nfail)

 ref(:,1) = (/1.0,52.92,0.9806,0.,0.,10.0,0./)
 ref(:,2) = (/6.331,257.1,0.4380,0.4069,0.,19.60,0./)
 ref(:,3) = ref(:,2)
 ref(:,4) = (/27.42,281.9,0.,-0.6811,0.,0.,0./)
 ref(:,5) = ref(:,4)
 ref(:,6) = (/6.331,257.1,-0.4380,0.4069,0.,-19.60,0./)
 ref(:,7) = ref(:,6)
 ref(:,8) = (/1.0,52.92,-0.9806,0.,0.,-10.0,0./)
 call check_regions(3,'collision',4./3.,ref,1.e-3,nfail)

end subroutine test_komissarov

!-----------------------------------------------------------------------
! Balsara (2001) tests 1-5 (GR06 Tables 7-11)
!-----------------------------------------------------------------------
subroutine test_balsara(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,8)

 ref(:,1) = (/1.0,1.625,0.,0.,0.,1.0,0./)
 ref(:,2) = (/0.6257,0.6989,0.3742,-0.03561,0.,0.6594,0./)
 ref(:,3) = ref(:,2)
 ref(:,4) = (/0.7092,0.7062,0.2555,-0.6804,0.,-0.4285,0./)
 ref(:,5) = (/0.2695,0.7062,0.2555,-0.6804,0.,-0.4285,0./)
 ref(:,6) = (/0.1223,0.6976,-0.02080,-0.003460,0.,-0.9769,0./)
 ref(:,7) = ref(:,6)
 ref(:,8) = (/0.125,0.725,0.,0.,0.,-1.0,0./)
 call check_regions(4,'Balsara 1',2.,ref,1.e-3,nfail)

 ref(:,1) = (/1.0,78.50,0.,0.,0.,6.0,6.0/)
 ref(:,2) = (/0.4300,23.21,0.6344,-0.09981,-0.09981,3.045,3.045/)
 ref(:,3) = ref(:,2)
 ref(:,4) = (/0.3830,22.84,0.6770,-0.05566,-0.05566,3.205,3.205/)
 ref(:,5) = (/2.828,22.84,0.6770,-0.05566,-0.05566,3.205,3.205/)
 ref(:,6) = (/1.582,20.72,0.4688,-0.2538,-0.2538,3.971,3.971/)
 ref(:,7) = ref(:,6)
 ref(:,8) = (/1.0,13.99,0.,0.,0.,0.7,0.7/)
 call check_regions(5,'Balsara 2',5./3.,ref,1.e-3,nfail)

 ref(:,1) = (/1.0,1099.,0.,0.,0.,7.0,7.0/)
 ref(:,2) = (/0.1381,86.04,0.9246,-0.03513,-0.03513,2.238,2.238/)
 ref(:,3) = ref(:,2)
 ref(:,4) = (/0.09798,76.53,0.9529,0.04366,0.04366,4.670,4.670/)
 ref(:,5) = (/10.10,76.53,0.9529,0.04366,0.04366,4.670,4.670/)
 ref(:,6) = (/1.218,63.63,0.4670,-0.4270,-0.4270,9.408,9.408/)
 ref(:,7) = ref(:,6)
 ref(:,8) = (/1.0,50.59,0.,0.,0.,0.7,0.7/)
 call check_regions(6,'Balsara 3',5./3.,ref,1.e-3,nfail)

 ref(:,1) = (/1.0,50.20,0.999,0.,0.,7.0,7.0/)
 ref(:,2) = (/51.75,1184.,0.04408,0.03263,0.03263,16.68,16.68/)
 ref(:,3) = ref(:,2)
 ref(:,4) = (/61.48,1188.,0.,-0.2877,-0.2877,0.,0./)
 ref(:,5) = ref(:,4)
 ref(:,6) = (/51.75,1184.,-0.04408,0.03263,0.03263,-16.68,-16.68/)
 ref(:,7) = ref(:,6)
 ref(:,8) = (/1.0,50.20,-0.999,0.,0.,-7.0,-7.0/)
 call check_regions(7,'Balsara 4',5./3.,ref,1.e-3,nfail)

 !--GR06 quote a reduced accuracy of 3e-4 for this test
 ref(:,1) = (/1.08,2.885,0.4,0.3,0.2,0.3,0.3/)
 ref(:,2) = (/2.447,5.908,-0.1331,0.2111,0.1751,0.2662,0.5076/)
 ref(:,3) = (/2.447,5.908,-0.1215,0.1264,0.1158,-0.1182,0.2302/)
 ref(:,4) = (/2.050,5.616,-0.04547,-0.1463,0.2146,-1.175,0.5852/)
 ref(:,5) = (/1.884,5.616,-0.04543,-0.1462,0.2149,-1.175,0.5850/)
 ref(:,6) = (/1.642,5.488,-0.1129,-0.04606,0.1601,-1.429,0.7320/)
 ref(:,7) = (/1.642,5.488,-0.1155,-0.08536,0.1027,-1.272,0.9468/)
 ref(:,8) = (/1.0,2.918,-0.45,-0.2,0.2,-0.7,0.5/)
 call check_regions(8,'Balsara 5',5./3.,ref,5.e-3,nfail)

end subroutine test_balsara

!-----------------------------------------------------------------------
! generic Alfven test (GR06 Table 12)
!-----------------------------------------------------------------------
subroutine test_alfven(nfail)
 integer, intent(inout) :: nfail
 real :: ref(7,8)

 ref(:,1) = (/1.0,23.76,0.,0.3,0.4,6.0,2.0/)
 ref(:,2) = (/0.9219,20.83,0.06232,0.3050,0.4193,5.622,1.892/)
 ref(:,3) = (/0.9219,20.83,0.07109,0.3669,0.2429,5.691,0.8502/)
 ref(:,4) = (/1.263,20.87,0.03886,0.1147,0.2054,5.130,0.7680/)
 ref(:,5) = (/1.099,20.87,0.03886,0.1147,0.2054,5.130,0.7680/)
 ref(:,6) = (/0.9130,20.85,0.01607,-0.05009,0.1813,5.505,0.8195/)
 ref(:,7) = (/0.9130,20.85,0.01341,-6.599e-4,-2.640e-4,5.073,2.029/)
 ref(:,8) = (/0.9,20.30,0.,0.,0.,5.0,2.0/)
 !--the tabulated R6 state does not satisfy the tangential momentum jump
 !  condition across the right slow shock (residual ~2e-4, ten times that
 !  from rounding to four figures), so R6 is compared with a looser tolerance
 call check_regions(9,'generic Alfven',5./3.,ref,1.e-3,nfail, &
                    tolreg=(/1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-2,1.e-3,1.e-3/))

end subroutine test_alfven

!-----------------------------------------------------------------------
! with zero magnetic field the solution should reduce to the
! special relativistic hydrodynamic Riemann problem
!-----------------------------------------------------------------------
subroutine test_hydro_limit(nfail)
 use shock_sr, only:exact_shock_sr
 integer, intent(inout) :: nfail
 real :: xs(6),yref(6),primL(7),primR(7),x
 integer :: npts,ierr,iv,j,k,ivsr(3)
 character(len=40) :: string

 primL = (/10.,13.33,0.,0.,0.,0.,0./)
 primR = (/1.,1.e-6,0.,0.,0.,0.,0./)
 ivsr  = (/1,2,3/)
 !--exact_shock_sr places the initial discontinuity at x = 0
 xs = (/-0.3,-0.15,-0.05,0.1,0.25,0.4/)
 do k=1,3
    iv = ivsr(k)
    call exact_mhdshock_sr(iv,0.4,5./3.,0.,primL,primR,-0.5,0.5,0.,xpts,ypts,npts,ierr)
    call exact_shock_sr(iv,0.4,5./3.,primL(1),primR(1),primL(2),primR(2),0.,0.,xs,yref,ierr)
    do j=1,size(xs)
       x = xs(j)
       write(string,"('SR hydro var ',i1,' at x=',f6.3)") iv,x
       call checkval(trim(string),value_at(x,xpts,ypts,npts),yref(j),2.e-3,nfail)
    enddo
 enddo

end subroutine test_hydro_limit

!-----------------------------------------------------------------------
! linear interpolation of the plotted solution at position x
!-----------------------------------------------------------------------
real function value_at(x,xp,yp,np)
 real, intent(in) :: x,xp(:),yp(:)
 integer, intent(in) :: np
 integer :: i

 value_at = yp(np)
 do i=2,np
    if (xp(i) >= x) then
       if (xp(i) > xp(i-1)) then
          value_at = yp(i-1) + (yp(i)-yp(i-1))*(x-xp(i-1))/(xp(i)-xp(i-1))
       else
          value_at = yp(i)
       endif
       return
    endif
 enddo

end function value_at

end module testmhdshocksr

program test_mhdshock_sr
 use testmhdshocksr
 implicit none
 integer :: nfail

 nfail = 0
 call test_tangential(nfail)
 call test_komissarov(nfail)
 call test_balsara(nfail)
 call test_alfven(nfail)
 call test_hydro_limit(nfail)

 if (nfail == 0) then
    print "(a)",' ALL SRMHD SHOCK TESTS PASSED'
 else
    print "(a,i4,a)",' FAILED ',nfail,' SRMHD SHOCK TESTS'
    stop 1
 endif

end program test_mhdshock_sr
