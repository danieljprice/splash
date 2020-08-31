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
!  Copyright (C) 2005-2020 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module providing library version of splash read_data routines
! specifies c interfaces to corresponding Fortran subroutines
!-------------------------------------------------------------------------
module libutils
  implicit none

  public

 contains

 subroutine check_argcv()
  include 'libinclude.f90'
 end subroutine check_argcv

function ctypes_to_fstring(array)
 use, intrinsic :: iso_c_binding, only:c_char
 character(kind=c_char), dimension(:), intent(in) :: array
 character(len=size(array)) :: ctypes_to_fstring
 integer :: i

 ctypes_to_fstring = ''
 do i=1,size(array)
    if (array(i)==achar(0)) exit
    ctypes_to_fstring(i:i) = array(i)
 enddo

end function ctypes_to_fstring

end module libutils
