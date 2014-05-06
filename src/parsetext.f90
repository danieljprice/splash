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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  module containing routines to parse text strings containing
!  variables.
!
!  For example:
!   %(t + 10)
!   %(t*2)
!   t = %t
!   time = %t.5   ! specifying formatting to 5 sig-figs
!
!  Uses the function parser module to evaluate the functions
!  once extracted from the text. Functions that cannot be
!  correctly evaluated are left intact in the text
!
!  Dependencies:
!     fparser module
!     plot_numb from the plotting library
!-----------------------------------------------------------------
module parsetext
 use fparser, only:rn
 implicit none

contains

subroutine parse_text(string,vars,vals)
 use asciiutils, only:string_replace,string_sub,lcase
 character(len=*), intent(inout) :: string
 real(kind=rn),    dimension(:), intent(in) :: vals
 character(len=*), dimension(:), intent(in) :: vars
 character(len=1) :: ch
 character(len=len(string)+128) :: newstring
 integer :: ia,iz,i0,i9,lenstr,npar,i,istart,iend,ndecimal,ierr,i1,i2
 logical :: in_variable,parse
 real :: r
 integer, parameter :: ndecimal_default = 3
 character(len=32) :: varstring
!
!--look for strings of the form:
!
!  %var
!  %(var)
!  %var.n
!  %(var + 10).n
!  %(var1*var2 + 10)
!  %(var1*var2)
!
 ia = iachar('a')
 iz = iachar('z')
 i0 = iachar('0')
 i9 = iachar('9')
 lenstr = len(newstring)
 parse = .false.
 in_variable = .false.
! print*,'parsing ',trim(string)
 newstring = string
 npar = 0
 ndecimal = ndecimal_default
 i1 = 0
 i2 = 0
 !
 ! In the following we extract two substrings:
 !
 ! string(istart:iend) is the string to replace, i.e. %(blah).5
 ! string(i1:i2) is the variable/function to evaluate, i.e. blah
 !
 i = 0
 do while (i < lenstr)
    i = i + 1
    ch = lcase(newstring(i:i))
    select case(ch)
    case('%')
       in_variable = .true.
       if (i.gt.1) then
          if (newstring(i-1:i-1).eq.'\') in_variable = .false.
       endif
       if (in_variable) then
          istart = i
          iend = 0
          i1 = i + 1
          i2 = 0
       endif
    case('(')
       if (in_variable) then
          npar = npar + 1
       endif
    case(')')
       if (in_variable) then
          npar = max(npar - 1,0)
          if (i.ge.lenstr) then
             iend = i
             if (i2 < i1) i2 = iend
             parse = .true.
          elseif (npar.eq.0) then
             if (newstring(i+1:i+1).ne.'.') then
                iend = i
                if (i2 <= i1) i2 = iend
                parse = .true.
             endif
          endif
       endif
    case('.')
       if (in_variable .and. npar <= 0 .and. i < lenstr) then
          read(newstring(i+1:i+1),"(i1)",iostat=ierr) ndecimal
          if (ierr.ne.0) ndecimal = 3
          iend = i+1
          if (i2 < i1) i2 = i - 1
          parse = .true.
       endif
    case default
       if ((.not.((iachar(ch) >= ia .and. iachar(ch) <= iz)  &
              .or.(iachar(ch) >= i0 .and. iachar(ch) <= i9)) &
            .or. i.eq.lenstr) .and. npar <= 0) then
          if (in_variable) then
             if (i.eq.lenstr) then
                iend = i
             else
                iend = i - 1
             endif
             if (i2 < i1) i2 = iend
             parse = .true.
          endif
       endif
    end select
    if (parse) then
       in_variable = .false.
       !print*,'variable = ',newstring(istart:iend), ', ndecimal = ',ndecimal
       !print*,'formula = ',newstring(i1:i2),i1,i2
       r = parse_formula(newstring(i1:i2),vars,vals,ierr)
       if (ierr.eq.0) then
          !print*,' r = ',r,' ierr = ',ierr
          call get_varstring(r,ndecimal,varstring)
          !print*,'varstring: "',varstring,'"'
          call string_sub(newstring,istart,iend,trim(varstring))
       endif
       i = i + (len_trim(varstring) - (iend - istart)) - 1
       !print*,'newstring = ',newstring  !(1:i),len_trim(varstring),iend-istart
       parse = .false.
       ndecimal = ndecimal_default
    endif
 enddo
!
! get rid of escape sequence on %
!
 call string_replace(newstring,'\%','%')
!
! replace original string (possibly truncated)
!
 string = trim(newstring)
! print*,' string: "',trim(newstring),'"'

end subroutine parse_text

!---------------------------------------------------------------------------
!
! write the real number r to a string
! with ndec decimal places
!
! uses plot_numb to do the formatting
!
!---------------------------------------------------------------------------
subroutine get_varstring(r,ndec,string)
 use plotlib, only:plot_numb
 real,    intent(in) :: r
 integer, intent(in)  :: ndec
 character(len=*), intent(out) :: string
 real :: rtmp
 integer :: mm,pp,nc

 if (abs(r) < tiny(r)) then
    string = '0'
    nc = 1
 else
    rtmp = abs(r)
    mm = nint(r/10.**(int(log10(rtmp)-ndec)))
    pp = int(log10(rtmp) - ndec)
    call plot_numb(mm,pp,1,string,nc)
 endif

end subroutine get_varstring

!---------------------------------------------------------------------------
!
! evaluate the variable or function via the function parser
!
! i.e. (t + 10) or (t*10) etc.
!
! OUTPUT: a real number
!
! unknown variables or un-parsable syntax return an error and a zero value
!
!---------------------------------------------------------------------------
real function parse_formula(string,vars,vals,ierr)
 use fparser, only:initf,evalf,endf,checkf,parsef
 character(len=*), intent(in) :: string
 character(len=*), dimension(:), intent(in) :: vars
 real(kind=rn),    dimension(:), intent(in) :: vals
 integer, intent(out) :: ierr
 
 call initf(1)
 ierr = checkf(string,vars,verbose=.false.)
 if (ierr.eq.0) then
    call parsef(1,string,vars)
    parse_formula = real(evalf(1,vals))
 else
    parse_formula = 0.
 endif
 call endf

end function parse_formula

end module parsetext
