!
! this module contains wrappers for all of the 
! system and compiler dependent routines
!
! these are called from the main program by their generic names,
! and in here the actual call to the system is performed
!
module system_commands
 !!use f90_unix  ! uncomment this for NAG f95 compiler
 implicit none

contains
 
 subroutine get_number_arguments(nargs)
    integer, intent(out) :: nargs
    integer :: iargc
    
    nargs = iargc()
        
 end subroutine get_number_arguments

 subroutine get_command_argument(iarg,argstring)
    integer, intent(in) :: iarg
    character(len=*), intent(out) :: argstring
    
    call getarg(iarg,argstring)
        
 end subroutine get_command_argument
 
end module system_commands
