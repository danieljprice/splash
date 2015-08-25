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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------
!
!  utility routines used during data read
!
!-----------------------------------------------------------
module dataread_utils
 use params, only:doub_prec
 implicit none

 public :: check_range
 integer, private :: iverbose_level = 1 ! can be changed
 
 private

 ! generic interface check_range
 interface check_range
  module procedure check_range_int, check_range_intarr,check_range_double
 end interface check_range

contains

!---------------------------------------------
!  set verboseness for remaining routines
!---------------------------------------------
subroutine set_check_range_verboseness(ilevel)
 integer, intent(in) :: ilevel

 iverbose_level = ilevel

end subroutine set_check_range_verboseness

!-----------------------------------------------
! standardised print statement for range errors
!-----------------------------------------------
subroutine handle_range_error(tag,string,ierror)
 character(len=*), intent(in)    :: tag,string
 integer,          intent(inout) :: ierror
 
 if (iverbose_level > 0) then
    print "(1x,5a)",'ERROR: ',trim(tag),' value of ',trim(string),' out of range'
 endif
 ierror = ierror + 1

end subroutine handle_range_error

!----------------------------------------------------
! check that an integer is within a prescribed range
!----------------------------------------------------
subroutine check_range_int(ivar,tag,min,max,err)
 integer,          intent(in) :: ivar
 character(len=*), intent(in) :: tag
 integer,          intent(in),  optional :: min,max
 integer,          intent(out), optional :: err
 integer :: ierror
 character(len=12) :: string
 
 write(string,"(i12)") ivar
 string = trim(adjustl(string))

 ierror = 0
 if (present(min)) then
    if (ivar < min) call handle_range_error(tag,string,ierror)
 endif
 if (present(max)) then
    if (ivar > max) call handle_range_error(tag,string,ierror)
 endif
 if (present(err)) then
    err = ierror
 endif

end subroutine check_range_int

!------------------------------------------------------------------------
! check that all values of an integer array is within a prescribed range
!------------------------------------------------------------------------
subroutine check_range_intarr(ivar,tag,min,max,err)
 integer,          intent(in) :: ivar(:)
 character(len=*), intent(in) :: tag
 integer,          intent(in),  optional :: min,max
 integer,          intent(out), optional :: err
 integer :: ierror,i
 character(len=12) :: string

 ierror = 0
 do i=1,size(ivar)
    write(string,"(i12)") ivar(i)
    string = trim(adjustl(string))

    if (present(min)) then
       if (ivar(i) < min) call handle_range_error(tag,string,ierror)
    endif
    if (present(max)) then
       if (ivar(i) > max) call handle_range_error(tag,string,ierror)
    endif
 enddo
 if (present(err)) then
    err = ierror
 endif

end subroutine check_range_intarr

!----------------------------------------------------
! check that a real*8 is within a prescribed range
!----------------------------------------------------
subroutine check_range_double(dvar,tag,min,max,err)
 real(doub_prec),  intent(in) :: dvar
 character(len=*), intent(in) :: tag
 real(doub_prec),  intent(in),  optional :: min,max
 integer,          intent(out), optional :: err
 integer :: ierror
 character(len=12) :: string
 
 write(string,"(1pg12.3)") dvar
 string = trim(adjustl(string))
 
 ierror = 0 
 if (present(min)) then
    if (dvar < min) call handle_range_error(tag,string,ierror)
 endif
 if (present(max)) then
    if (dvar > max) call handle_range_error(tag,string,ierror)
 endif
 if (present(err)) then
    err = ierror
 endif

end subroutine check_range_double

end module dataread_utils
