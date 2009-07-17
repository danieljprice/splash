!
! this module contains wrappers for all of the 
! system and compiler dependent routines
!
! these are called from the main program by their generic names,
! and in here the actual call to the system is performed
!
! THIS ONE IS FOR MANY UNIX COMPILERS
!
module system_commands
 implicit none

contains
 
 subroutine get_number_arguments(nargs)
    integer, intent(out) :: nargs
    integer :: iargc
    
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

end module system_commands
