!
! this module contains certain useful utilities which
! depend on system commands (in the system_commands module)
!
module system_utils
 use system_commands, only:get_environment
 implicit none
 public :: lenvironment,renvironment

contains

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
    
    call get_environment(variable,string)
    read(string,*,iostat=ierr) renvironment
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
    
    call get_environment(variable,string)
    if (trim(string).eq.'yes'.or.trim(string).eq.'YES' &
    .or.trim(string).eq.'true'.or.trim(string).eq.'TRUE' &
    .or.trim(string).eq.'on'.or.trim(string).eq.'ON') then
       lenvironment = .true.
    else
       lenvironment = .false.
    endif
 
 end function lenvironment
 
end module system_utils
