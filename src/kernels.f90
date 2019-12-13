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
!  Copyright (C) 2005-2018 Daniel Price. All rights reserved.
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

 real, parameter :: pi = 4.*atan(1.)
 integer, public :: ikernel = 0
 real, public  :: radkernel = 2.
 real, public  :: radkernel2 = 4.
 real, public  :: cnormk1D = 2./3.
 real, public  :: cnormk2D = 10./(7.*pi)
 real, public  :: cnormk3D = 1./pi
 procedure(k_func), pointer, public :: wfunc

 abstract interface
  pure function k_func(q)
   real, intent(in) :: q
   real :: k_func
  end function k_func
 end interface

 public :: select_kernel, select_kernel_by_name
 public :: pint

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

 if (j >= 1 .and. j <= nkernels) then
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
    if (j==1) then
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
 if (ikernel==0) then
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
pure real function w_cubic(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2 < 1.0) then
    q = sqrt(q2)
    w_cubic = 1.-1.5*q2 + 0.75*q2*q
 elseif (q2 < 4.0) then
    q = sqrt(q2)
    w_cubic = 0.25*(2.-q)**3
 else
    w_cubic = 0.
 endif

end function w_cubic

pure real function w_quartic(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 q = sqrt(q2)
 if (q < 0.5) then
    w_quartic = (2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4
 elseif (q < 1.5) then
    w_quartic = (2.5-q)**4 - 5.*(1.5-q)**4
 elseif (q < 2.5) then
    w_quartic = (2.5-q)**4
 else
    w_quartic = 0.
 endif

end function w_quartic

pure real function w_quintic(q2)
 implicit none
 real, intent(in) :: q2
 real :: q,q4

 if (q2 < 1.0) then
    q = sqrt(q2)
    q4 = q2*q2
    w_quintic = 66.-60.*q2 + 30.*q4 - 10.*q4*q
 elseif ((q2 >= 1.0).and.(q2 < 4.0)) then
    q = sqrt(q2)
    w_quintic = (3.-q)**5 - 6.*(2.-q)**5
 elseif ((q2 >= 4.0).and.(q2 < 9.0)) then
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
 if (q < 0.4) then
    w_quartic2h = (2.-q)**4 - 5.*(1.2-q)**4 + 10.*(0.4-q)**4
 elseif (q < 1.2) then
    w_quartic2h = (2.-q)**4 - 5.*(1.2-q)**4
 elseif (q < 2.) then
    w_quartic2h = (2.-q)**4
 else
    w_quartic2h = 0.
 endif

end function w_quartic2h

pure real function w_wendlandc2(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2 < 4.) then
    q = sqrt(q2)
    w_wendlandc2 = (1. - 0.5*q)**4*(2.*q + 1.)
 else
    w_wendlandc2 = 0.
 endif

end function w_wendlandc2

pure real function w_wendlandc4(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2 < 4.) then
    q = sqrt(q2)
    w_wendlandc4 = (1. - 0.5*q)**6*(35./12.*q2 + 3.*q + 1.)
 else
    w_wendlandc4 = 0.
 endif

end function w_wendlandc4

pure real function w_wendlandc6(q2)
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2 < 4.) then
    q = sqrt(q2)
    w_wendlandc6 = (1. - 0.5*q)**8*(4.*q2*q + 25./4.*q2 + 4.*q + 1.)
 else
    w_wendlandc6 = 0.
 endif

end function w_wendlandc6

!---------------------------------------
!
!  Integrals of the kernel, used
!  when performing exact interpolation
!  See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
pure real function pint(r0, d1, d2, hi1)
 real, intent(in) :: r0, d1, d2, hi1
 real :: ar0,q0,tphi1,tphi2,phi1,phi2

 if (abs(r0) < tiny(0.)) then
    pint = 0.
    return
 elseif (r0 > 0.) then
    pint = 1.
    ar0 = r0
 else
    pint = -1.
    ar0 = -r0
 endif
 q0 = ar0*hi1
 tphi1 = abs(d1)/ar0
 tphi2 = abs(d2)/ar0
 phi1 = atan(tphi1)
 phi2 = atan(tphi2)

 if (d1*d2 >= 0.) then
    pint = pint*(full_2d_mod(phi1,tphi1,q0) + full_2d_mod(phi2,tphi2,q0))
 else
    if (abs(d1) < abs(d2)) then
       pint = pint*(full_2d_mod(phi2,tphi2,q0) - full_2d_mod(phi1,tphi1,q0))
    else
       pint = pint*(full_2d_mod(phi1,tphi1,q0) - full_2d_mod(phi2,tphi2,q0))
    endif
 endif

end function pint

!---------------------------------------
!
! Helper functions for kernel integrals
! See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
pure real function full_2d_mod(phi,tphi,q0)
 real, intent(in) :: phi,tphi,q0
 real :: q, phi1, phi2, tphi1, tphi2, cphi, cphi1, cphi2

 if (q0 <= 1.0) then
    cphi = cos(phi)
    q = q0/cphi

    if (q <= 1.0) then
       full_2d_mod = F1_2d(phi,tphi,cphi,q0)
    elseif (q <= 2.0) then
       cphi1 = q0
       phi1 = acos(q0)
       tphi1 = tan(phi1)
       full_2d_mod = F2_2d(phi,tphi,cphi,q0) - F2_2d(phi1,tphi1,cphi1,q0) &
                   + F1_2d(phi1,tphi1,cphi1,q0)
    else
       cphi1 = q0
       phi1 = acos(q0)
       cphi2 = 0.5*q0
       phi2 = acos(0.5*q0)
       tphi1 = tan(phi1)
       tphi2 = tan(phi2)
       full_2d_mod = F3_2d(phi) - F3_2d(phi2) + F2_2d(phi2,tphi2,cphi2,q0) &
                                - F2_2d(phi1,tphi1,cphi1,q0) &
                                + F1_2d(phi1,tphi1,cphi1,q0)
    endif
 elseif (q0 <= 2.0) then
    cphi = cos(phi)
    q = q0/cphi

    if (q <= 2.0) then
       full_2d_mod = F2_2d(phi,tphi,cphi,q0)
    else
       cphi2 = 0.5*q0
       phi2 = acos(0.5*q0)
       tphi2 = tan(phi2)
       full_2d_mod = F3_2d(phi) - F3_2d(phi2) + F2_2d(phi2,tphi2,cphi2,q0)
    endif
 else ! q0 > 2
    full_2d_mod = F3_2d(phi)
 endif

end function full_2d_mod

pure real function F1_2d(phi,tphi,cphi,q0)
 real, intent(in) :: phi,tphi,cphi,q0
 real :: I2, I4, I5, logs, cphi2, q02, q03

 !tphi = tan(phi)
 !cphi = cos(phi)
 cphi2 = cphi*cphi

 q02 = q0*q0
 q03 = q02*q0

 logs = log(tan(phi/2.+pi/4.))

 I2 = tphi
 I4 = 1./3. * tphi*(2. + 1./cphi2)

 I5 = 1./16. * (0.5*(11.*sin(phi) + 3.*sin(3.*phi))/cphi2/cphi2 + 6.*logs)

 F1_2d =  5./7.*q02/pi  * (I2 - 3./4.*q02 *I4 + 0.3*q03 *I5)

end function F1_2d

pure real function F2_2d(phi, tphi, cphi, q0)
 real, intent(in) :: phi, tphi, cphi, q0
 real :: I0, I2, I3, I4, I5, logs, cphi2, q02, q03

! tphi = tan(phi)
! cphi = cos(phi)
 cphi2 = cphi*cphi

 q02 = q0*q0
 q03 = q02*q0

 logs = log(tan(phi/2.+pi/4.))

 I0 = phi
 I2 = tphi
 I4 = 1./3. * tphi*(2. + 1./(cphi2))

 I3 = 1./2.  * (tphi/cphi + logs)
 I5 = 1./16. * (0.5*(11.*sin(phi) + 3.*sin(3.*phi))/cphi2/cphi2 + 6.*logs)

 F2_2d =  5./7.*q02/pi  * (2.*I2 - 2.*q0 *I3 + 3./4.*q02 *I4 - 1./10.*q03 *I5 &
                           - 1./10./q02 *I0)

end function F2_2d

pure real function F3_2d(phi)
 real, intent(in) :: phi
 real :: I0

 I0 = phi
 F3_2d = 0.5/pi  *I0

end function F3_2d

end module kernels
