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

!-------------------------------------------------------------
!
!  Module containing kernel functions used for interpolation
!
!-------------------------------------------------------------
module kernels
 implicit none
 integer, parameter, public :: nkernels = 3
 character(len=22), dimension(0:nkernels), parameter, public :: kernelname = &
    (/'default [cubic]       ', &
      'M4 cubic spline   (2h)', &
      'M5 quartic        (2h)', &
      'M6 quintic spline (3h)'/)

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
    print "(a,/)",' Using '//trim(kernelname(j))//' kernel' 
 else
    print "(a,/)",' Using default '//trim(kernelname(1))//' kernel'
 endif

 select case(j)
 case(3) ! M6 quintic
    ikernel = 3
    radkernel = 3.0
    cnormk1D = 1./120.
    cnormk2D = 7./(478*pi)
    cnormk3D = 1./(120.*pi)
    wfunc => w_quintic
 case(2) ! M5 quintic, squashed to 2h
    ikernel = 2
    radkernel = 2.0
    cnormk1D = 3125./24576.
    cnormk2D = 46875./(153472.*pi)
    cnormk3D = 15625./(65536.*pi)
    wfunc => w_quartic2h
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
    case('quintic','quintic spline','m6')
       jkern = 3
    case('quartic','quartic2h','m5','squashed quartic','quartic spline')
       jkern = 2
    case('cubic','cubic spline','m4')
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
pure real function w_cubic(q2)
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

pure real function w_quintic(q2)
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

pure real function w_quartic2h(q2)
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

end module kernels
