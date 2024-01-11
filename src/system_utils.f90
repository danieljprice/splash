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
!  Copyright (C) 2005-2021 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!
! this module contains certain useful utilities which
! depend on system commands (in the system_commands module)
!
module system_utils
 use asciiutils, only:lcase
 implicit none
 public :: ienvironment,lenvironment,renvironment,lenvstring,ienvstring
 public :: envlist,ienvlist,lenvlist,renvlist,get_command_option,count_matching_args
 public :: get_command_flag,get_user,get_copyright,get_environment_or_flag

 private

contains

!
!--this routine returns a variable
!  from EITHER a command line flag --foo=bar
!  OR from an environment variable MY_FOO=bar
!
!  The flag takes priority over the environment variable
!
subroutine get_environment_or_flag(variable,string)
 character(len=*), intent(in)  :: variable
 character(len=*), intent(out) :: string
 integer :: ierr,iloc

 ! try as command flag, excluding everything before the last
 ! underscore, e.g. MY_GOOD_FOO becomes --foo
 iloc = index(variable,'_',back=.true. )
 call get_option(variable(iloc+1:),string,ierr)

 ! then try as environment variable, including underscores
 ! e.g. MY_GOOD_FOO=bar
 if (ierr /= 0) call get_environment_variable(variable,string)
 !print*,' GOT ',trim(variable),' = ',trim(string)

end subroutine get_environment_or_flag

 !
 !--this routine returns an integer variable
 !  from an environment variable setting
 !
 !  if the errval argument is present then this is the
 !  value assigned when an error has occurred
 !  (default is zero)
 !
integer function ienvironment(variable,errval)
 character(len=*), intent(in) :: variable
 character(len=30) :: string
 integer, intent(in), optional :: errval

 call get_environment_or_flag(variable,string)
 if (present(errval)) then
    ienvironment = ienvstring(string,errval)
 else
    ienvironment = ienvstring(string)
 endif

end function ienvironment

 !
 !--this routine returns a floating point (real) variable
 !  from an environment variable setting
 !
 !  if the errval argument is present then this is the
 !  value assigned when an error has occurred
 !  (default is zero)
 !
real function renvironment(variable,errval)
 character(len=*), intent(in) :: variable
 character(len=30) :: string
 real, intent(in), optional :: errval
 integer :: ierr

 call get_environment_or_flag(variable,string)
 if (len_trim(string) > 0) then
    read(string,*,iostat=ierr) renvironment
 else
    ierr = 1
 endif

 if (ierr /= 0) then
    if (present(errval)) then
       renvironment = errval
    else
       renvironment = 0.
    endif
 endif

end function renvironment

 !
 !--this routine returns a logical variable
 !  from an environment variable setting
 !
logical function lenvironment(variable)
 character(len=*), intent(in) :: variable
 character(len=30) :: string

 call get_environment_or_flag(variable,string)
 lenvironment = lenvstring(string)

end function lenvironment

 !
 !--utility routine to determine whether a string
 !  should be interpreted as true or false
 !
logical function lenvstring(string)
 character(len=*), intent(in) :: string

 if (string(1:1)=='y'.or.string(1:1)=='Y' &
    .or.string(1:1)=='t'.or.string(1:1)=='T' &
    .or.trim(string)=='on'.or.trim(string)=='ON' &
    .or.trim(string)=='1') then
    lenvstring = .true.
 else
    lenvstring = .false.
 endif

end function lenvstring

 !
 !--utility routine to extract integer value from string
 !
integer function ienvstring(string,errval)
 character(len=*), intent(in)  :: string
 integer, intent(in), optional :: errval
 character(len=5) :: fmtstring
 integer :: ierr

 if (len_trim(string) > 0) then
    !--use a formatted read - this is to avoid a compiler bug
    !  should in general be more robust anyway
    write(fmtstring,"(a,i2,a)",iostat=ierr) '(i',len_trim(string),')'
    read(string,fmtstring,iostat=ierr) ienvstring
 else
    ierr = 1
 endif

 if (ierr /= 0) then
    if (present(errval)) then
       ienvstring = errval
    else
       ienvstring = 0
    endif
 endif

end function ienvstring
!
 !--utility routine to extract real value from string
 !
real function renvstring(string,errval)
 character(len=*), intent(in)  :: string
 real, intent(in), optional :: errval
 character(len=5) :: fmtstring
 integer :: ierr

 if (len_trim(string) > 0) then
    read(string,*,iostat=ierr) renvstring
 else
    ierr = 1
 endif

 if (ierr /= 0) then
    if (present(errval)) then
       renvstring = errval
    else
       renvstring = 0.
    endif
 endif

end function renvstring
 !
 !--this routine returns an arbitrary number of
 !  comma separated strings
 !
subroutine envlist(variable,nlist,list)
 character(len=*), intent(in) :: variable
 integer, intent(out) :: nlist
 character(len=*), dimension(:), intent(out), optional :: list
 character(len=120) :: string
 character(len=10) :: dummy
 integer :: i1,i2,ierr
 logical :: notlistfull

 !--set list to blank strings if argument is present
 if (present(list)) then
    list = ' '
 endif

 !--get envlist from the environment
 call get_environment_or_flag(variable,string)

 !--split the string on commas
 i1 = 1
 i2 = index(string,',')-1
 if (i2==-1) i2 = len_trim(string)
 nlist = 0
 ierr = 0
 notlistfull = .true.

 !--for each comma separated string, add a list member
 do while(i2 >= i1 .and. notlistfull .and. ierr==0)
    nlist = nlist + 1
    !print*,'i1,i2,stringfrag= ',i1,i2,trim(string(i1:))
    if (present(list)) then
       read(string(i1:i2),"(a)",iostat=ierr) list(nlist)
       notlistfull = (nlist < size(list))
    else
       read(string(i1:i2),"(a)",iostat=ierr) dummy ! to get ierr at end of string
       notlistfull = .true.
    endif
    i1 = i2 + 2
    i2 = index(string(i1:),',')
    if (i2==0) then
       i2 = len_trim(string)
    else
       i2 = i2 + i1 - 2
    endif
 enddo

end subroutine envlist
!
!--return comma separated list of integers
!
function ienvlist(variable,nlist,errval)
 character(len=*), intent(in) :: variable
 integer, intent(in) :: nlist
 integer, intent(in), optional :: errval
 character(len=30), dimension(nlist) :: list
 integer :: ienvlist(nlist),i,ngot

 ngot = nlist
 call envlist(variable,ngot,list)

 do i=1,nlist
    if (present(errval)) then
       ienvlist(i) = ienvstring(list(i),errval=errval)
    else
       ienvlist(i) = ienvstring(list(i))
    endif
 enddo

end function ienvlist
!
!--return comma separated list of logicals
!
function lenvlist(variable,nlist)
 character(len=*), intent(in) :: variable
 integer, intent(in) :: nlist
 character(len=30), dimension(nlist) :: list
 logical :: lenvlist(nlist)
 integer :: i,ngot

 ngot = nlist
 call envlist(variable,ngot,list)

 do i=1,nlist
    lenvlist(i) = lenvstring(list(i))
 enddo

end function lenvlist
!
!--return comma separated list of reals
!
function renvlist(variable,nlist,errval)
 character(len=*), intent(in) :: variable
 integer, intent(in) :: nlist
 real, intent(in), optional :: errval
 character(len=30), dimension(nlist) :: list
 real :: renvlist(nlist)
 integer :: i,ngot

 ngot = nlist
 call envlist(variable,ngot,list)

 do i=1,nlist
    if (present(errval)) then
       renvlist(i) = renvstring(list(i),errval=errval)
    else
       renvlist(i) = renvstring(list(i))
    endif
 enddo

end function renvlist

!
!--find logical-valued option from command line arguments
!  as in --arg (true if present, false if not)
!
subroutine get_option(variable,value,err)
 character(len=*), intent(in) :: variable
 character(len=*), intent(out) :: value
 character(len=80) :: string
 integer, intent(out) :: err
 integer :: nargs,iarg,ieq

 err = 1
 nargs = command_argument_count()
 do iarg=1,nargs
    call get_command_argument(iarg,string)
    if (string(1:2)=='--' .and. index(lcase(string),lcase(variable)) > 0) then
       err = 0
       ieq = index(string,'=',back=.true.)
       if (ieq > 0) then
          value = string(ieq+1:)
       else
          value = 'True'  ! raw flag --foo, equivalent to --foo=True
       endif
    endif
 enddo

end subroutine get_option
!
!--find real-valued option from command line arguments
!  as in --arg=blah
!
real function get_command_option(variable,default) result(val)
 character(len=*), intent(in) :: variable
 real, intent(in), optional   :: default
 character(len=80) :: string
 integer :: ierr,nargs,ieq,iarg

 val = 0.
 if (present(default)) val = default
 nargs = command_argument_count()
 do iarg=1,nargs
    call get_command_argument(iarg,string)
    ieq = index(string,'=')
    if (string(1:1)=='-' .and. index(string,variable) > 0 .and. ieq > 0) then
       read(string(ieq+1:),*,iostat=ierr) val
    endif
 enddo

end function get_command_option

!
!--find logical-valued option from command line arguments
!  as in --arg (true if present, false if not)
!
logical function get_command_flag(variable) result(val)
 character(len=*), intent(in) :: variable
 character(len=80) :: string
 integer :: nargs,iarg

 val = .false.
 nargs = command_argument_count()
 do iarg=1,nargs
    call get_command_argument(iarg,string)
    if (string(1:1)=='-' .and. index(string,variable) > 0) val = .true.
 enddo

end function get_command_flag

!
!--count the number of arguments matching a certain substring
!  e.g. with particular filename extension
!
integer function count_matching_args(string,id) result(n)
 character(len=*), intent(in) :: string
 integer, intent(out), optional :: id(:)
 character(len=1024) :: myarg
 integer :: iarg,nargs

 n = 0
 nargs = command_argument_count()
 do iarg=1,nargs
    call get_command_argument(iarg,myarg)
    if (index(myarg,string) > 0) then
       n = n + 1
       if (present(id)) then
          if (size(id) >= n) id(n) = iarg
       endif
    endif
 enddo

end function count_matching_args

!
!--get the name of the logged in user
!
subroutine get_user(string)
 character(len=*), intent(out) :: string
 character(len=*), parameter :: tempfile = '/tmp/splash.username'
 logical :: iexist
 integer :: ierr,iu

 string = ''
 inquire(file='/proc/cpuinfo',exist=iexist) ! exists only in Linux
 if (.not.iexist) then
    call system('id -F $USER > '//trim(tempfile)) ! Mac
 else
    call system('getent passwd "$USER" | cut -d: -f5 | cut -d, -f1 > '//trim(tempfile)) ! Linux
 endif
 inquire(file=tempfile,exist=iexist)
 if (iexist) then
    open(newunit=iu,file=tempfile,action='read',status='old',iostat=ierr)
    read(iu,"(a)",iostat=ierr) string
    close(iu,status='delete',iostat=ierr)
 endif

end subroutine get_user

!
!--get copyright string involving logged in user and current year
!
function get_copyright() result(string)
 character(len=30) :: string
 character(len=8) :: year

 call get_user(string)
 call date_and_time(date=year)

 string = '(c) '//year(1:4)//' '//trim(string)

end function get_copyright

end module system_utils
