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
 use kernel_integrals, only:maxpieces,maxpow,get_kernel_integral_coeffs, &
                            edge_int_2d,edge_int_3d
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
 integer, parameter :: doub_prec = kind(0.d0)

 ! coefficient tables for exact kernel integrals (set in select_kernel)
 integer :: npieces_int = 0, nmax_int = 0
 real(doub_prec) :: qb_int(maxpieces)
 real(doub_prec) :: g2_int(-1:maxpow,maxpieces+1),a3_int(-1:maxpow,maxpieces+1)

 procedure(k_func), pointer, public :: wfunc
 procedure(k_func), pointer, public :: dwfunc

 abstract interface
  pure function k_func(q)
   real, intent(in) :: q
   real :: k_func
  end function k_func
 end interface

 public :: select_kernel, select_kernel_by_name
 public :: pint,wallint,soft_func

 private

contains

!-----------------------------------------------
!
! Kernel selection routine, sets radkernel,
! cnormk values and pointer to kernel function
!
!--------------------------------------
subroutine select_kernel(j,verbose)
 integer, intent(in) :: j
 logical, intent(in), optional :: verbose
 logical :: isverbose

 isverbose = .true.
 if (present(verbose)) isverbose = verbose

 if (j >= 1 .and. j <= nkernels) then
    !--print only if NOT using the default kernel
    if (isverbose) print "(a,/)",' Using '//trim(kernelname(j))//' kernel'
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
    dwfunc => dw_cubic
 end select
 radkernel2 = radkernel*radkernel
 call get_kernel_integral_coeffs(ikernel,npieces_int,nmax_int,qb_int,g2_int,a3_int)

end subroutine select_kernel

!--------------------------------------
!
! Kernel selection based on string
!
!--------------------------------------
subroutine select_kernel_by_name(string,verbose)
 use asciiutils, only:lcase
 character(len=*), intent(in) :: string
 logical, intent(in), optional :: verbose
 integer :: i,jkern
 logical :: isverbose

 isverbose = .true.
 if (present(verbose)) isverbose = verbose

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
 if (jkern==0) then
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

 if (jkern == 0) then
    ! keep the configured kernel rather than silently falling back to cubic
    if (isverbose) print "(a)",' ERROR: unknown kernel "'//trim(string)// &
         '", keeping current setting'
    call select_kernel(ikernel,verbose=isverbose)
 else
    call select_kernel(jkern,verbose=isverbose)
 endif

end subroutine select_kernel_by_name

!---------------------------------------
!
!  Functional forms of various kernels
!
!--------------------------------------
pure real function w_cubic(q2)
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

pure real function dw_cubic(q2)
 real, intent(in) :: q2
 real :: q

 if (q2 < 1.0) then
    q = sqrt(q2)
    dw_cubic = q*(2.25*q - 3.)
 elseif (q2 < 4.0) then
    q = sqrt(q2)
    dw_cubic = -0.75*(2.-q)**2
 else
    dw_cubic = 0.
 endif

end function dw_cubic

pure real function w_quartic(q2)
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
!  generalised to any kernel in the
!  kernel_integrals module
!
!---------------------------------------
pure real function pint(r0, d1, d2, hi1)
 real, intent(in) :: r0, d1, d2, hi1
 real(doub_prec) :: q0, v1, v2

 if (abs(r0) < tiny(0.)) then
    pint = 0.
    return
 elseif (r0 > 0.) then
    pint = 1.
 else
    pint = -1.
 endif
 q0 = abs(r0)*hi1
 call edge_int_2d(abs(d1)*real(hi1,doub_prec),abs(d2)*real(hi1,doub_prec),q0, &
                  npieces_int,nmax_int,qb_int,g2_int,v1,v2)
 pint = pint*real(combine_vertices(d1,d2,v1,v2))

end function pint

!---------------------------------------
!
! combine the two vertex integrals on
! an edge depending on whether the foot
! of the perpendicular lies between them
!
!---------------------------------------
pure real(doub_prec) function combine_vertices(d1,d2,v1,v2) result(v)
 real,            intent(in) :: d1, d2
 real(doub_prec), intent(in) :: v1, v2

 if (d1*d2 >= 0.) then
    v = v1 + v2
 elseif (abs(d1) < abs(d2)) then
    v = v2 - v1
 else
    v = v1 - v2
 endif

end function combine_vertices

!------------------------------------------------------------
! 3D functions to evaluate exact overlap of kernel with wall boundaries
! see Petkova, Laibe & Bonnell (2018), J. Comp. Phys
!------------------------------------------------------------
pure real function wallint(r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi)
 real, intent(in) :: r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi
 real :: dx, dy

 dx = xc - xp
 dy = yc - yp

 !
 ! Contributions from each of the 4 sides of a cell wall
 !
 wallint = pint3D(r0, 0.5*pixwidthy + dy, 0.5*pixwidthx - dx, 0.5*pixwidthx + dx, hi) &
         + pint3D(r0, 0.5*pixwidthy - dy, 0.5*pixwidthx + dx, 0.5*pixwidthx - dx, hi) &
         + pint3D(r0, 0.5*pixwidthx + dx, 0.5*pixwidthy + dy, 0.5*pixwidthy - dy, hi) &
         + pint3D(r0, 0.5*pixwidthx - dx, 0.5*pixwidthy - dy, 0.5*pixwidthy + dy, hi)

end function wallint

!------------------------------------------------------------
! contribution from one side of a cell wall: r0 is the distance
! to the wall plane, R_0 the distance in the plane to the edge
! and d1, d2 the distances along the edge to each vertex
!------------------------------------------------------------
pure real function pint3D(r0, R_0, d1, d2, h)
 real, intent(in) :: r0, R_0, d1, d2, h
 real(doub_prec) :: c, p, hi1, int1, int2

 if (abs(r0) < tiny(0.)) then
    pint3D = 0.
    return
 endif

 pint3D = sign(1.,r0)*sign(1.,R_0)
 hi1 = 1.d0/h
 c = abs(r0)*hi1
 p = abs(R_0)*hi1
 call edge_int_3d(abs(d1)*hi1,abs(d2)*hi1,p,c,npieces_int,nmax_int,qb_int,a3_int,int1,int2)
 int1 = max(int1,0.d0)
 int2 = max(int2,0.d0)
 pint3D = pint3D*real(max(combine_vertices(d1,d2,int1,int2),0.d0))

end function pint3D

!------------------------------------------------------------
! function to return a soft maximum for 1/x with no bias
! for x >> eps using the cubic spline kernel softening
! i.e. something equivalent to 1/sqrt(x**2 + eps**2) but
! with compact support, i.e. f=1/x when x > 2*eps
!------------------------------------------------------------
pure elemental real function soft_func(x,eps) result(f)
 real, intent(in)  :: x,eps
 real :: q,q2,q4

 q = x/eps
 q2 = q*q
 if (q < 1.) then
    q4 = q2*q2
    f = (1./eps)*(q4*q/10. - 3.*q4/10. + 2.*q2/3. - 7./5.)
 elseif (q < 2.) then
    q4 = q2*q2
    f = (1./eps)*(q*(-q4*q + 9.*q4 - 30.*q2*q + 40.*q2 - 48.) + 2.)/(30.*q)
 else
    f = -1./x
 endif
 f = -f

end function soft_func

end module kernels
