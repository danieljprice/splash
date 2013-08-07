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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------
!
!  Module containing kernel functions used for interpolation
!
!-------------------------------------------------------------
module kernels
 implicit none
 integer, parameter, public :: nkernels = 6
 character(len=24), dimension(0:nkernels), parameter, public :: kernelname = &
    (/'default [cubic]         ', &
      'M4 cubic spline   (2h)  ', &
      'M5 quartic        (2.5h)', &
      'M6 quintic spline (3h)  ', &
      'Wendland C2       (2h)  ', &
      'Wendland C4       (2h)  ', &
      'Wendland C6       (2h)  '/)

 real, parameter :: pi = 3.14159236
 integer, public :: ikernel = 0
 real, public  :: radkernel = 2.
 real, public  :: radkernel2 = 4.
 real, public  :: cnormk1D = 2./3.
 real, public  :: cnormk2D = 10./(7.*pi)
 real, public  :: cnormk3D = 1./pi
 procedure(real), pointer :: wfunc
 
 public :: wfunc, select_kernel, select_kernel_by_name
 
 private

contains

!-----------------------------------------------
!
! Kernel selection routine, sets radkernel,
! cnormk values and pointer to kernel function
!
!--------------------------------------
subroutine select_kernel(j)
 integer, intent(in) :: j

 if (j.ge.1 .and. j.le.nkernels) then
    !--print only if NOT using the default kernel
    print "(a,/)",' Using '//trim(kernelname(j))//' kernel' 
 endif

 select case(j)
 case(6) ! Wendland 3D C6
    ikernel = 6
    radkernel = 2.
    cnormk1D = 15./16.
    cnormk2D = 39./(14.*pi)
    cnormk3D = 1365./(512.*pi)
    wfunc => w_wendlandc6
 case(5) ! Wendland 3D C4
    ikernel = 5
    radkernel = 2.
    cnormk1D = 27./32.
    cnormk2D = 9./(4.*pi)
    cnormk3D = 495./(256.*pi)
    wfunc => w_wendlandc4
 case(4) ! Wendland 3D C2
    ikernel = 4
    radkernel = 2.
    cnormk1D = 0.75
    cnormk2D = 7./(4.*pi)
    cnormk3D = 21./(16.*pi)
    wfunc => w_wendlandc2
 case(3) ! M6 quintic, 3h
    ikernel = 3
    radkernel = 3.0
    cnormk1D = 1./120.
    cnormk2D = 7./(478*pi)
    cnormk3D = 1./(120.*pi)
    wfunc => w_quintic
 case(2) ! M5 quartic, 2.5h
    ikernel = 2
    radkernel = 2.5
    cnormk1D = 1./24.
    cnormk2D = 96./(1199.*pi)
    cnormk3D = 1./(20.*pi)
    wfunc => w_quartic
 case default  !-- cubic spline kernel
    if (j.eq.1) then
       ikernel = 1 ! deliberately chose cubic spline
    else
       ikernel = 0 ! just whatever is the default
    endif
    radkernel = 2.0
    cnormk1D = 2./3.
    cnormk2D = 10./(7.*pi)
    cnormk3D = 1./pi
    wfunc => w_cubic
 end select
 radkernel2 = radkernel*radkernel

end subroutine select_kernel

!--------------------------------------
!
! Kernel selection based on string
!
!--------------------------------------
subroutine select_kernel_by_name(string)
 use asciiutils, only:lcase
 character(len=*), intent(in) :: string
 integer :: i,jkern
 
 jkern = 0
 !
 !--check if string exactly matches a kernel name
 !
 do i=1,nkernels
    if (trim(adjustl(lcase(string)))==trim(adjustl(lcase(kernelname(i))))) then
       jkern = i
    endif
 enddo
 !
 !--if no match to a kernel name, look for other possible strings
 !
 if (ikernel.eq.0) then
    select case(trim(adjustl(lcase(string))))
    case('wendlandc6','wendland c6','6th order wendland','wendland 3d c6','w6','wendland6')
       jkern = 6
    case('wendlandc4','wendland c4','4th order wendland','wendland 3d c4','w4','wendland4')
       jkern = 5
    case('wendlandc2','wendland c2','2nd order wendland','wendland 3d c2','w2','wendland2')
       jkern = 4
    case('quintic','quintic spline','m6','quintic b-spline')
       jkern = 3
    case('quartic','quartic spline','m5','quartic b-spline')
       jkern = 2
    case('cubic','cubic spline','m4','cubic b-spline')
       jkern = 1
    end select
 endif

 call select_kernel(jkern)

end subroutine select_kernel_by_name

!---------------------------------------
!
!  Functional forms of various kernels
!
!--------------------------------------
real function w_cubic(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2.lt.1.0) then
    q = sqrt(q2)
    w_cubic = 1.-1.5*q2 + 0.75*q2*q
 elseif (q2.lt.4.0) then
    q = sqrt(q2)
    w_cubic = 0.25*(2.-q)**3
 else
    w_cubic = 0.
 endif

end function w_cubic

real function w_quartic(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 q = sqrt(q2)
 if (q.lt.0.5) then
    w_quartic = (2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4
 elseif (q.lt.1.5) then
    w_quartic = (2.5-q)**4 - 5.*(1.5-q)**4
 elseif (q.lt.2.5) then
    w_quartic = (2.5-q)**4
 else
    w_quartic = 0.
 endif

end function w_quartic

real function w_quintic(q2)
 implicit none
 real, intent(in) :: q2
 real :: q,q4

 if (q2.lt.1.0) then
    q = sqrt(q2)
    q4 = q2*q2
    w_quintic = 66.-60.*q2 + 30.*q4 - 10.*q4*q
 elseif ((q2.ge.1.0).and.(q2.lt.4.0)) then
    q = sqrt(q2)
    w_quintic = (3.-q)**5 - 6.*(2.-q)**5
 elseif ((q2.ge.4.0).and.(q2.lt.9.0)) then
    q = sqrt(q2)
    w_quintic = (3.-q)**5
 else
    w_quintic = 0.0
 endif

end function w_quintic

real function w_quartic2h(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 q = sqrt(q2)
 if (q.lt.0.4) then
    w_quartic2h = (2.-q)**4 - 5.*(1.2-q)**4 + 10.*(0.4-q)**4
 elseif (q.lt.1.2) then
    w_quartic2h = (2.-q)**4 - 5.*(1.2-q)**4
 elseif (q.lt.2.) then
    w_quartic2h = (2.-q)**4
 else
    w_quartic2h = 0.
 endif

end function w_quartic2h

real function w_wendlandc2(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc2 = (1. - 0.5*q)**4*(2.*q + 1.)
 else
    w_wendlandc2 = 0.
 endif

end function w_wendlandc2

real function w_wendlandc4(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc4 = (1. - 0.5*q)**6*(35./12.*q2 + 3.*q + 1.)
 else
    w_wendlandc4 = 0.
 endif

end function w_wendlandc4

real function w_wendlandc6(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc6 = (1. - 0.5*q)**8*(4.*q2*q + 25./4.*q2 + 4.*q + 1.)
 else
    w_wendlandc6 = 0.
 endif

end function w_wendlandc6

end module kernels
