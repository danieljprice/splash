!
! this module contains wrappers for all of the 
! system and compiler dependent routines
!
! these are called from the main program by their generic names,
! and in here the actual call to the system is performed
!
! THIS VERSION IS FOR THE NAG f95 COMPILER
!
module system_commands
 use f90_unix  ! uncomment this for NAG f95 compiler
 implicit none

contains
 
 subroutine get_number_arguments(nargs)
    integer, intent(out) :: nargs
    
    nargs = iargc()
        
 end subroutine get_number_arguments

 subroutine get_argument(iarg,argstring)
    integer, intent(in) :: iarg
    character(len=*), intent(out) :: argstring
    
    call getarg(iarg,argstring)
        
 end subroutine get_argument
 
 subroutine get_environment(variable,value)
    character(len=*), intent(in) :: variable
    character(len=*), intent(out) :: value
 
    call getenv(variable,value)
 
 end subroutine get_environment

 !--this routine returns a logical variable
 !  from an environment variable setting
 !  and should not require modification
 ! (just acts as an interface to get_environment)
 !
 logical function lenvironment(variable)
    character(len=*), intent(in) :: variable
    character(len=10) :: string
    
    call get_environment(variable,string)
    if (trim(string).eq.'yes'.or.trim(string).eq.'YES' &
    .or.trim(string).eq.'true'.or.trim(string).eq.'TRUE' &
    .or.trim(string).eq.'on'.or.trim(string).eq.'ON') then
       lenvironment = .true.
    else
       lenvironment = .false.
    endif
 
 end function lenvironment
 
end module system_commands
