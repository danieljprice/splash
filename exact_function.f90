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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!----------------------------------------------------------------------
! Plots arbitrary analytic function y = f(x)
! Uses the function parser module by Roland Schmehl
!----------------------------------------------------------------------
module exactfunction
  implicit none

contains

subroutine exact_function(string,xplot,yplot,ierr)
  use fparser, only:initf,parsef,evalf,endf,EvalErrType,EvalErrMsg,rn
  implicit none
  character(len=*), intent(in) :: string
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  integer :: i
  real(kind=rn), dimension(1) :: xi
  character(len=1), dimension(1), parameter :: var = (/'x'/)
  
  
  print "(a)",' Plotting function f(x) = '//trim(string)
  if (len_trim(string).le.0) then
     print "(a)",' *** ERROR: blank function string in exact_function call'
     ierr = 1
     return
  endif
  ierr = 0
  call initf(1)
  call parsef(1,string,var)
  
  if (EvalErrType.ne.0) then
     print "(a)",' *** ERROR parsing function: '//trim(EvalerrMsg())//' ***'
     ierr = EvalErrType
  else
     do i=1,size(xplot)
        xi = xplot(i)                 ! type conversion here
        yplot(i) = real(evalf(1,xi))  ! type conversion back
        if (EvalErrType > 0) ierr = EvalErrType
     enddo
     if (ierr.ne.0) then
        print "(a)",' *** WARNING: errors during function evaluation: '//trim(EvalerrMsg())
        !--set exit error to zero so we plot the results anyway
        ierr = 0
     endif
  endif
  call endf
    
  return
end subroutine exact_function

!----------------------------------------------------------------
! check syntax in the function string - this subroutine 
! mainly just an interface to checking routines in fparser
!----------------------------------------------------------------
subroutine check_function(string,ierr)
 use fparser, only:checkf
 implicit none
 character(len=*), intent(in) :: string
 integer, intent(out)         :: ierr
 
 ierr = checkf(string,(/'x'/))

end subroutine check_function

end module exactfunction
