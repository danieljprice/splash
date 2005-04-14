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
 
end module system_commands
