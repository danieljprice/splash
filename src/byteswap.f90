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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module byteswap
 implicit none
 integer, parameter :: little_endian_order8(8) = [8,7,6,5,4,3,2,1]
 integer, parameter :: little_endian_order4(4) = [4,3,2,1]

 public :: bs

 interface bs
  module procedure reverse_bytes_r4,reverse_bytes_r8,&
                   reverse_bytes_int,reverse_bytes_int8
 end interface bs

 private

contains
!--------------------------------------------------------
! function to swap byte order (real 4)
!--------------------------------------------------------
real(kind=4) elemental function reverse_bytes_r4(x) result(y)
 real(kind=4), intent(in) :: x
 character(len=1)         :: src_arr(4), dst_arr(4)

 ! put whatever the variable was into an 8 x 1 character string
 src_arr = transfer(x,src_arr)
 ! swap the bytes around into the desired order
 dst_arr = src_arr(little_endian_order4)
 ! copy the 8 x 1 character string back to represent the original type
 y = transfer(dst_arr,y)

end function reverse_bytes_r4

!--------------------------------------------------------
! function to swap byte order (real 8)
!--------------------------------------------------------
real(kind=8) elemental function reverse_bytes_r8(x) result(y)
 real(kind=8), intent(in) :: x
 character(len=1)         :: src_arr(8), dst_arr(8)

 ! put whatever the variable was into an 8 x 1 character string
 src_arr = transfer(x,src_arr)
 ! swap the bytes around into the desired order
 dst_arr = src_arr(little_endian_order8)
 ! copy the 8 x 1 character string back to represent the original type
 y = transfer(dst_arr,y)

end function reverse_bytes_r8

!--------------------------------------------------------
! function to swap byte order (int 4)
!--------------------------------------------------------
integer(kind=4) elemental function reverse_bytes_int(x) result(y)
 integer(kind=4), intent(in) :: x
 character(len=1)            :: src_arr(4), dst_arr(4)

 ! put whatever the variable was into an 8 x 1 character string
 src_arr = transfer(x,src_arr)
 ! swap the bytes around into the desired order
 dst_arr = src_arr(little_endian_order4)
 ! copy the 8 x 1 character string back to represent the original type
 y = transfer(dst_arr,y)

end function reverse_bytes_int

!--------------------------------------------------------
! function to swap byte order (int 4)
!--------------------------------------------------------
integer(kind=8) elemental function reverse_bytes_int8(x) result(y)
 integer(kind=8), intent(in) :: x
 character(len=1)            :: src_arr(8), dst_arr(8)

 ! put whatever the variable was into an 8 x 1 character string
 src_arr = transfer(x,src_arr)
 ! swap the bytes around into the desired order
 dst_arr = src_arr(little_endian_order8)
 ! copy the 8 x 1 character string back to represent the original type
 y = transfer(dst_arr,y)

end function reverse_bytes_int8

end module byteswap
