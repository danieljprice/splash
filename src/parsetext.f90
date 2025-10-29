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
!  Copyright (C) 2005-2025 Daniel Price. All rights reserved.
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
!   %t:datetime   ! formats as YYYY-MM-DD (rounded to days)
!   %(t + 2013-04-10):datetime  ! formats as YYYY-MM-DD from base date (t in seconds)
!   %(t + 2013/04/10):datetime  ! formats as YYYY/MM/DD from base date (t in seconds)
!   %(t + 2013-04-10 12:00):dt  ! formats as YYYY-MM-DD HH:MM:SS from base date/time (t in seconds)
!   %(t + 2013/04/10 12:00):dt  ! formats as YYYY/MM/DD HH:MM:SS from base date/time (t in seconds)
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
 character(len=len(string)+128) :: temp_formula
 integer :: ia,iz,i0,i9,lenstr,npar,i,istart,iend,ndecimal,ierr,i1,i2
 logical :: in_variable,parse,is_datetime
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
!  %(t + 2013-04-10):datetime
!  %(t + 2013-04-10 12:00):dt
!  %t:datetime
!  %t:dt
!
 ia = iachar('a')
 iz = iachar('z')
 i0 = iachar('0')
 i9 = iachar('9')
 lenstr = len(newstring)
 parse = .false.
 in_variable = .false.
 is_datetime = .false.
! print*,'parsing ',trim(string)
 newstring = string
 npar = 0
 ndecimal = ndecimal_default
 i1 = 0
 i2 = 0
 !
 ! In the following we extract two substrings:
 !
 ! string(istart:iend) is the string to replace, i.e. %(blah).5 or %(blah):datetime
 ! string(i1:i2) is the variable/function to evaluate, i.e. blah
 !
 i = 0
 do while (i < lenstr)
    i = i + 1
    ch = lcase(newstring(i:i))
    select case(ch)
    case('%')
       in_variable = .true.
       if (i > 1) then
          if (newstring(i-1:i-1)=='\') in_variable = .false.
       endif
       if (in_variable) then
          istart = i
          iend = 0
          i1 = i + 1
          i2 = 0
          is_datetime = .false.
       endif
    case('(')
       if (in_variable) npar = npar + 1
    case(')')
       if (in_variable) then
          npar = max(npar - 1,0)
          if (i >= lenstr) then
             iend = i
             if (i2 < i1) i2 = iend
             parse = .true.
          elseif (npar==0) then
             if (newstring(i+1:i+1) /= '.' .and. newstring(i+1:i+1) /= ':') then
                iend = i
                if (i2 <= i1) i2 = iend
                parse = .true.
             endif
          endif
       endif
    case('.')
       if (in_variable .and. npar <= 0 .and. i < lenstr) then
          read(newstring(i+1:i+1),"(i1)",iostat=ierr) ndecimal
          if (ierr /= 0) ndecimal = 3
          iend = i+1
          if (i2 < i1) i2 = i - 1
          parse = .true.
       endif
    case(':')
       if (in_variable .and. npar <= 0 .and. i < lenstr) then
          ! check for datetime format
          if (i+8 <= lenstr .and. lcase(newstring(i+1:i+8)) == 'datetime') then
             is_datetime = .true.; iend = i+8
          elseif (i+2 <= lenstr .and. lcase(newstring(i+1:i+2)) == 'dt') then
             is_datetime = .true.; iend = i+2
          else
             ! not a datetime format, continue parsing as normal
             is_datetime = .false.
          endif
          if (is_datetime) then
             if (i2 < i1) i2 = i - 1
             parse = .true.
          endif
       endif
    case default
       if ((.not.((iachar(ch) >= ia .and. iachar(ch) <= iz)  &
              .or.(iachar(ch) >= i0 .and. iachar(ch) <= i9)) &
            .or. i==lenstr) .and. npar <= 0) then
          if (in_variable) then
             if (i==lenstr) then
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
       if (is_datetime) then
          call create_parseable_formula(newstring(i1:i2), temp_formula)
          r = parse_formula(temp_formula,vars,vals,ierr)
          if (ierr==0) then
             call get_datetime_string(r, varstring, newstring(i1:i2))
          else
             varstring = 'ERROR'
          endif
       else
          r = parse_formula(newstring(i1:i2),vars,vals,ierr)
          if (ierr==0) then
             call number_to_string(r,ndecimal,varstring)
          else
             varstring = 'ERROR'
          endif
       endif

       call string_sub(newstring,istart,iend,trim(varstring))
       i = i + (len_trim(varstring) - (iend - istart)) - 1
       parse = .false.; ndecimal = ndecimal_default; is_datetime = .false.
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
subroutine number_to_string(r,ndec,string)
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

end subroutine number_to_string

!---------------------------------------------------------------------------
!
! write the real number r (time in seconds) to a datetime string
! format can be 'datetime' or 'dt'
!
! converts seconds since epoch to datetime format
! base_date is optional and can be in format 'YYYY-MM-DD' or 'YYYY-MM-DD HH:MM' or 'YYYY/MM/DD HH:MM'
! If base_date has no time component, the output will be rounded to days
! original_format preserves the separators and format from the input
!
!---------------------------------------------------------------------------
subroutine get_datetime_string(r,string,original_format)
 real,    intent(in) :: r
 character(len=*), intent(out) :: string
 character(len=*), intent(in), optional :: original_format
 integer :: year, month, day, hour, minute, second
 character(len=1)  :: date_separator
 integer :: seconds_since_epoch, total_seconds
 logical :: is_datetime_format

 ! convert real to integer seconds
 seconds_since_epoch = nint(r)

 ! determine separators and format from original format if provided
 date_separator = '-'
 is_datetime_format = .true.
 if (present(original_format)) then
    call determine_separator(original_format, date_separator)
    is_datetime_format = index(original_format, ':') > 0 ! if has time component HH:MM
 endif

 ! convert seconds since Unix epoch to date components
 year = 1970; month = 1; day = 1
 total_seconds = seconds_since_epoch
 day = day + total_seconds / 86400
 total_seconds = mod(total_seconds, 86400)

 ! handle negative remaining seconds
 if (total_seconds < 0) then
    total_seconds = total_seconds + 86400
    day = day - 1
 endif

 ! convert remaining seconds to time components
 hour = total_seconds / 3600
 minute = (total_seconds - hour * 3600) / 60
 second = total_seconds - hour * 3600 - minute * 60

 ! normalize the date by working backwards from the epoch
 do while (day > days_in_year(year))
    day = day - days_in_year(year)
    year = year + 1
 enddo

 do while (day > days_in_month(month, year))
    day = day - days_in_month(month, year)
    month = month + 1
    if (month > 12) then
       month = 1
       year = year + 1
    endif
 enddo

 ! format the date string
 write(string,"(i4,a1,i2.2,a1,i2.2)") year, date_separator, month, date_separator, day

 ! add time if datetime format is requested
 if (is_datetime_format) then
    write(string,"(a,1x,i2.2,':',i2.2)") trim(string), hour, minute
 endif

end subroutine get_datetime_string

!---------------------------------------------------------------------------
!
! Create a parseable formula by replacing datetime strings with time
! since a reference epoch in seconds
! For example: "t + 2013-04-10 12:00" becomes "t + 1365595200"
!
!---------------------------------------------------------------------------
subroutine create_parseable_formula(original_formula, parseable_formula)
 use asciiutils, only:string_replace
 character(len=*), intent(in) :: original_formula
 character(len=*), intent(out) :: parseable_formula
 integer :: i, len_formula
 integer :: year, month, day, hour, minute, second
 integer :: seconds_since_epoch
 character(len=32) :: temp_str
 character(len=32) :: datetime_str

 parseable_formula = original_formula
 len_formula = len_trim(parseable_formula)

 i = 1
 do while (i <= len_formula)
    if (parse_datetime_pattern(original_formula,i,year,month,day,hour,minute,second,datetime_str)) then
       ! convert to seconds since Unix epoch (1970-01-01)
       seconds_since_epoch = datetime_to_seconds(year, month, day, hour, minute, second)

       ! replace the datetime string with the value in seconds since 1970
       write(temp_str, '(i0)') seconds_since_epoch
       call string_replace(parseable_formula, trim(datetime_str), trim(temp_str))
       len_formula = len_trim(parseable_formula)
    endif
    i = i + 1
 enddo

end subroutine create_parseable_formula

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
 if (ierr==0) then
    call parsef(1,string,vars)
    parse_formula = real(evalf(1,vals))
 else
    parse_formula = 0.
 endif
 call endf

end function parse_formula

!---------------------------------------------------------------------------
!
! Determine date separators from the original format string
!
!---------------------------------------------------------------------------
subroutine determine_separator(format_str, date_sep)
 character(len=*), intent(in) :: format_str
 character(len=1), intent(out) :: date_sep

 date_sep = '-'
 if (index(format_str, '/') > 0) date_sep = '/'

end subroutine determine_separator

!---------------------------------------------------------------------------
!
! Convert datetime components to seconds since Unix epoch (1970-01-01 00:00:00)
!
!---------------------------------------------------------------------------
integer function datetime_to_seconds(year,month,day,hour,minute,second)
 integer, intent(in) :: year,month,day,hour,minute,second
 integer :: days_since_epoch,i

 ! days since Unix epoch (1970-01-01)
 days_since_epoch = 0

 ! add days for each year from 1970 to year-1
 do i = 1970, year - 1
    days_since_epoch = days_since_epoch + days_in_year(i)
 enddo

 ! add days for months in current year (up to month-1)
 do i = 1, month - 1
    days_since_epoch = days_since_epoch + days_in_month(i, year)
 enddo

 ! add days in current month (day-1 because we start from day 1)
 days_since_epoch = days_since_epoch + day - 1

 ! convert to seconds and add time components
 datetime_to_seconds = days_since_epoch * 86400 + hour * 3600 + minute * 60 + second

end function datetime_to_seconds

!---------------------------------------------------------------------------
!
! Parse a datetime string and extract year, month, day, hour, minute, second
! returns .true. if a valid datetime pattern is found
!
!---------------------------------------------------------------------------
logical function parse_datetime_pattern(str,start_pos,year,month,day,hour,minute,second,datetime_str)
 character(len=*), intent(in) :: str
 integer, intent(in) :: start_pos
 integer, intent(out) :: year, month, day, hour, minute, second
 character(len=*), intent(out) :: datetime_str
 integer :: len_str,ierr
 character(len=1) :: separator,sep2

 parse_datetime_pattern = .false.
 len_str = len_trim(str)

 if (start_pos + 9 > len_str) return

 ! Check for YYYY-MM-DD or YYYY/MM/DD pattern
 separator = str(start_pos+4:start_pos+4)
 sep2      = str(start_pos+7:start_pos+7)
 if ((separator /= '-' .and. separator /= '/') .or. sep2 /= separator) return

 ! parse date components using formatted read
 read(str(start_pos:start_pos+9), '(i4,1x,i2,1x,i2)',iostat=ierr) year, month, day
 if (ierr /= 0) return

 ! Check if there's a time component
 hour = 0; minute = 0; second = 0
 datetime_str = str(start_pos:start_pos+9)

 if (start_pos + 15 <= len_str .and. str(start_pos+10:start_pos+10) == ' ') then
    if (str(start_pos+13:start_pos+13) == ':') then
       ! parse time components HH:MM using formatted read
       read(str(start_pos+11:start_pos+15), '(i2,1x,i2)', iostat=ierr) hour, minute
       if (ierr /= 0) return
       datetime_str = str(start_pos:start_pos+15)
       ! read seconds if HH:MM:SS
       if (str(start_pos+16:start_pos+16) == ':' .and. start_pos + 18 <= len_str) then
          read(str(start_pos+17:start_pos+18), *, iostat=ierr) second
          if (ierr /= 0) return
          datetime_str = str(start_pos:start_pos+18)
       endif
    endif
 endif

 parse_datetime_pattern = .true.

end function parse_datetime_pattern

!---------------------------------------------------------------------------
!
! check if a year is a leap year
!
!---------------------------------------------------------------------------
logical function is_leap_year(year)
 integer, intent(in) :: year

 is_leap_year = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. (mod(year, 400) == 0)

end function is_leap_year

!---------------------------------------------------------------------------
!
! get number of days in a month (accounting for leap years)
!
!---------------------------------------------------------------------------
integer function days_in_month(month,year)
 integer, intent(in) :: month, year
 integer, parameter :: days_in_month_array(12) = &
          (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

 days_in_month = days_in_month_array(month)
 if (month == 2 .and. is_leap_year(year)) days_in_month = 29

end function days_in_month

!---------------------------------------------------------------------------
!
! get number of days in a year (accounting for leap years)
!
!---------------------------------------------------------------------------
integer function days_in_year(year)
 integer, intent(in) :: year

 days_in_year = 365 + merge(1, 0, is_leap_year(year))

end function days_in_year

end module parsetext
