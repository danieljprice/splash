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
  use fparser, only:initf,evalf,endf,EvalErrType,EvalErrMsg,rn
  implicit none
  character(len=*), intent(in) :: string
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  integer :: i,j,nfunc
  real(kind=rn), dimension(:), allocatable     :: val
  
  print "(a)",' Plotting function f(x) = '//trim(string)
  if (len_trim(string).le.0) then
     print "(a)",' *** ERROR: blank function string in exact_function call'
     ierr = 1
     return
  endif
  ierr = 0
  !
  !--work out how many subfunctions the string contains 
  !  and allocate memory for the sub function values appropriately
  !
  call get_nfunc(string,nfunc)
  allocate(val(nfunc),stat=ierr)
  if (ierr /= 0) then
     print "(a)",' ERROR allocating memory for ',nfunc,' sub-functions in exact_function'
     if (allocated(val)) deallocate(val)
     return
  endif

  call initf(nfunc)
  call parse_subfunctions(string,nfunc,.false.,ierr)
  
  if (EvalErrType.ne.0) then
     print "(a)",' *** ERROR parsing function: '//trim(EvalerrMsg())//' ***'
     ierr = EvalErrType
  else
     do i=1,size(xplot)
        val(1) = xplot(i)                 ! type conversion here
        
        !--evaluate sub-functions in order of dependency
        do j=2,nfunc
           val(j) = evalf(j,val(1:j-1))
        enddo
        yplot(i) = real(evalf(1,val(1:nfunc)))  ! type conversion back
        if (EvalErrType /= 0) ierr = EvalErrType
     enddo
     if (ierr.ne.0) then
        print "(a)",' *** WARNING: errors during function evaluation: '//trim(EvalerrMsg())
        !--set exit error to zero so we plot the results anyway
        ierr = 0
     endif
  endif

  call endf
  if (allocated(val)) deallocate(val)
    
  return
end subroutine exact_function

!----------------------------------------------------------------
! check syntax in the function string - this subroutine 
! mainly just an interface to checking routines in fparser
!----------------------------------------------------------------
subroutine check_function(string,ierr)
! use fparser, only:checkf
 implicit none
 character(len=*), intent(in) :: string
 integer, intent(out)         :: ierr
 integer :: nfunc
 
 call get_nfunc(string,nfunc)
 call parse_subfunctions(string,nfunc,.true.,ierr)
! ierr = checkf(string,(/'x'/))

end subroutine check_function

!----------------------------------------------------------------
! allow sub-function syntax (f(x) = y, y = 24*x)
!----------------------------------------------------------------
subroutine parse_subfunctions(string,nfunc,check,ierr)
 use fparser, only:checkf,parsef,EvalErrMsg,EvalErrType
 implicit none
 character(len=*), intent(in) :: string
 integer, intent(in)          :: nfunc
 logical, intent(in)          :: check
 integer, intent(out)         :: ierr
 
 character(len=len(string)), dimension(nfunc) :: var
 integer :: ieq,ivars,lstr,j,icommaprev

 var(1) = 'x'
 ivars = 1
 lstr = len_trim(string)
 icommaprev = lstr+1
 do j=lstr,1,-1
!
!--split the string according to commas
!
    if (string(j:j)==',') then
       !--sub functions must be of the form f(var) = val
       ieq = j + index(string(j+1:lstr),'=')
       if (ieq.eq.j) then
          print "(a)",'*** Error in sub-function syntax, missing equals sign in comma-separated function list'
          ierr = 4
          return
       endif
       !--variable is what lies to left of equals sign
       ivars = ivars + 1
       var(ivars) = string(j+1:ieq-1)
       if (len_trim(var(ivars)).le.0) then
          print "(a)",'*** Error in sub-function syntax, blank variable '
          ierr = 3
          return
       endif
       !--function is what lies to right of equals sign
       if (check) then
          if (ivars.eq.2) print "(a)",'Evaluating sub-functions in the order:'
          print "(a)",trim(var(ivars))//' = '//string(ieq+1:icommaprev-1)
          ierr = checkf(string(ieq+1:icommaprev-1),var(1:ivars-1))
          if (ierr /= 0) return
       else
          call parsef(ivars,string(ieq+1:icommaprev-1),var(1:ivars-1))
          if (EvalErrType.ne.0) then
             print "(a)",' *** ERROR parsing function: '//trim(EvalerrMsg())//' ***'
             ierr = EvalErrType
             return
          endif
       endif
       icommaprev = j
    endif
 enddo
 if (ivars.ne.nfunc) then
    print "(a)",' Internal consistency error in parse_subfunctions:'
    print*,' nsubfunctions ',ivars,' not equal to that obtained in get_nfunc, ',nfunc
 endif
!
!--finally, check/parse combined function
!
 if (check) then
    if (ivars.ge.2) print "(a)",'f('//trim(var(1))//') = '//string(1:icommaprev-1)
    ierr = checkf(string(1:icommaprev-1),var(1:ivars))
 else
    call parsef(1,string(1:icommaprev-1),var(1:ivars))
    if (EvalErrType.ne.0) then
       print "(a)",' *** ERROR parsing function: '//trim(EvalerrMsg())//' ***'
       ierr = EvalErrType
    endif
 endif
 
end subroutine parse_subfunctions

!----------------------------------------------------------------
! query the number of sub-functions (number of commas)
!----------------------------------------------------------------
subroutine get_nfunc(string,nfunc)
 implicit none
 character(len=*), intent(in) :: string
 integer, intent(out) :: nfunc
 integer :: j
 
 nfunc = 1
 do j=1,len_trim(string)
    if (string(j:j)==',') nfunc = nfunc + 1
 enddo

end subroutine get_nfunc

end module exactfunction
