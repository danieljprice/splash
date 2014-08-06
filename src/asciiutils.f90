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

!---------------------------------------------------------------------------
! module containing various utility subroutines
! related to reading from ascii files and dealing with string variables
!
! written by Daniel Price, University of Exeter 2007 24th April '07
! revised at Monash University, Nov '08.
! daniel.price@monash.edu
!
! this is a standalone module with no dependencies
!---------------------------------------------------------------------------
module asciiutils
 implicit none
 public :: read_asciifile,get_ncolumns,get_nrows,ncolumnsline,safename,basename
 public :: cstring,fstring
 public :: string_replace, string_delete, nheaderlines, string_sub
 public :: ucase,lcase
 public :: get_line_containing

 private

!--------------------------------------------------
! Generic interface to ascii file read for either
! character arrays (ie. each line is an element)
! or an array of real numbers
!--------------------------------------------------
 interface read_asciifile
   module procedure read_asciifile_char, read_asciifile_real,&
                    read_asciifile_real_string
 end interface read_asciifile

contains

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns array of character strings (one per line)
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile_char(filename,nlinesread,charline,ierror)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 character(len=*), dimension(:), intent(out) :: charline
 integer, intent(out), optional :: ierror
 integer, parameter :: iunit = 66 ! logical unit number for read operation
 integer :: ierr,i,maxlines
 logical :: iexist

 nlinesread = 0
 if (present(ierror)) ierror = 0

 !--if file does not exist, do nothing and return
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    if (present(ierror)) ierror = -1
    return
 endif

 open(unit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
 !--error opening file (but file does exist)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
    if (present(ierror)) ierror = ierr
    return
 endif

 maxlines = size(charline)
 do i=1,maxlines
    read(iunit,"(a)",err=66,end=99) charline(i)
 enddo
 !--end of array limits
 !  check to see if there is anything more in the file. Report error if there is.
 read(iunit,"(a)",iostat=ierr)
 if (ierr.eq.0) then
    print "(a,i6)",' WARNING: array limits reached reading '//trim(filename)//', max = ',maxlines
 endif
 nlinesread = maxlines
 close(unit=iunit)
 return

 !--error encountered
66 continue
  print "(a,i6)",' ERROR reading '//trim(filename)//' at line ',i-1
  if (present(ierror)) ierror = 1
  nlinesread = i-1
  close(unit=iunit)
  return

 !--reached end of file (the expected behaviour)
99 continue
  nlinesread = i-1
  close(unit=iunit)
  return

end subroutine read_asciifile_char

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns array of real numbers (either one per line or all on same line)
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile_real(filename,nlinesread,realarr,ierror)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 real, dimension(:), intent(out) :: realarr
 integer, intent(out), optional :: ierror
 integer, parameter :: iunit = 66 ! logical unit number for read operation
 integer :: ierr,i,maxlines
 logical :: iexist

 nlinesread = 0
 if (present(ierror)) ierror = 0

 !--if file does not exist, do nothing and return
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    if (present(ierror)) ierror = -1
    return
 endif

 open(unit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
 !--error opening file (but file does exist)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
    if (present(ierror)) then
       ierror = ierr
    endif
    return
 endif

 realarr(:) = -666.
 maxlines = size(realarr)
 read(iunit,*,err=66,end=99) (realarr(i),i=1,maxlines)

 !--end of array limits
 print "(a,i6)",' WARNING: array limits reached reading '//trim(filename)//', max = ',maxlines
 nlinesread = maxlines
 close(unit=iunit)
 return

 !--error encountered
66 continue
  print "(a,i6)",' ERROR reading '//trim(filename)//' at line ',i-1
  if (present(ierror)) ierror = 1
  do i=1,maxlines
     if (abs(realarr(i)+666.).gt.tiny(0.)) nlinesread = nlinesread + 1
  enddo
  close(unit=iunit)
  return

 !--reached end of file (the expected behaviour)
99 continue
  do i=1,maxlines
     if (abs(realarr(i)+666.).gt.tiny(0.)) nlinesread = nlinesread + 1
  enddo
  close(unit=iunit)
  return

end subroutine read_asciifile_real

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns array of real numbers and corresponding string
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile_real_string(filename,nlinesread,realarr,charline,ierror)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 real, dimension(:), intent(out) :: realarr
 character(len=*), dimension(:), intent(out) :: charline
 integer, intent(out), optional :: ierror
 integer, parameter :: iunit = 66 ! logical unit number for read operation
 integer :: ierr,i,maxlines
 logical :: iexist

 nlinesread = 0
 if (present(ierror)) ierror = 0

 !--if file does not exist, do nothing and return
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    if (present(ierror)) ierror = -1
    return
 endif

 open(unit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
 !--error opening file (but file does exist)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
    if (present(ierror)) then
       ierror = ierr
    endif
    return
 endif

 if (size(realarr) /= size(charline)) then
    print "(a)",' WARNING: array size mismatch in call to read_asciifile'
 endif

 realarr(:) = -666.
 maxlines = min(size(realarr),size(charline))
 read(iunit,*,err=66,end=99) (realarr(i),charline(i),i=1,maxlines)

 !--end of array limits
 print "(a,i6)",' WARNING: array limits reached reading '//trim(filename)//', max = ',maxlines
 nlinesread = maxlines
 close(unit=iunit)
 return

 !--error encountered
66 continue
  print "(a,i6)",' ERROR reading '//trim(filename)//' at line ',i-1
  if (present(ierror)) ierror = 1
  do i=1,maxlines
     if (abs(realarr(i)+666.).gt.tiny(0.)) nlinesread = nlinesread + 1
  enddo
  close(unit=iunit)
  return

 !--reached end of file (the expected behaviour)
99 continue
  do i=1,maxlines
     if (abs(realarr(i)+666.).gt.tiny(0.)) nlinesread = nlinesread + 1
  enddo

  close(unit=iunit)
  return

end subroutine read_asciifile_real_string

!---------------------------------------------------------------------------
! utility to work out number of columns of real numbers
! in an ascii file
!
! file must already be open and at the start
! slightly ad-hoc but its the best way I could think of!
!---------------------------------------------------------------------------
subroutine get_ncolumns(lunit,ncolumns,nheaderlines)
 implicit none
 integer, intent(in) :: lunit
 integer, intent(out) :: ncolumns,nheaderlines
 integer :: ierr,ncolprev,ncolsthisline
 character(len=5000) :: line
 logical :: nansinfile,infsinfile

 nheaderlines = 0
 line = ' '
 ierr = 0
 ncolumns = 0
 ncolprev = -100
 ncolsthisline = 0
 nansinfile = .false.
 infsinfile = .false.
!
!--loop until we find two consecutive lines with the same number of columns (but non zero)
!
 do while ((len_trim(line).eq.0 .or. ncolsthisline.ne.ncolprev .or. ncolumns.le.0) .and. ierr.eq.0)
    ncolprev = ncolumns
    read(lunit,"(a)",iostat=ierr) line
    if (index(line,'NaN').gt.0) nansinfile = .true.
    if (index(line,'Inf').gt.0) infsinfile = .true.
    if (len_trim(line).eq.0) then
       ncolsthisline = -1
    else
       if (ierr.eq.0) ncolsthisline = ncolumnsline(line)
       ncolumns = ncolsthisline
    endif
    nheaderlines = nheaderlines + 1
    !print*,'DEBUG: header line ',nheaderlines,' ncols = ',ncolsthisline,'"'//trim(line)//'"'
 enddo
 !--subtract 2 from the header line count (the last two lines which were the same)
 nheaderlines = max(nheaderlines - 2,0)
 if (ierr .gt.0 .or. ncolumns.le.0) then
    ncolumns = 0
 elseif (ierr .lt. 0) then
    !print*,ncolumns,ncolprev
 endif
 if (nansinfile) print "(a)",' INDIAN BREAD WARNING!! NaNs in file!!'
 if (infsinfile) print "(a)",' WARNING!! Infs in file!!'
 rewind(lunit)

 if (ncolumns.eq.0) print "(a)",' ERROR: no columns of real numbers found'

end subroutine get_ncolumns

!---------------------------------------------------------------------------
! utility to work out number of rows in file
!---------------------------------------------------------------------------
subroutine get_nrows(lunit,nheaderlines,nlines)
 implicit none
 integer, intent(in)  :: lunit,nheaderlines
 integer, intent(out) :: nlines
 integer :: ierr,i

 rewind(lunit)
 ierr = 0
 do i=1,nheaderlines
    read(lunit,*,iostat=ierr)
 enddo
 nlines = 0
 do while (ierr==0)
    read(lunit,*,iostat=ierr)
    if (ierr==0) nlines = nlines + 1
 enddo

end subroutine get_nrows

!---------------------------------------------------------------------------
!
! function returning the number of columns of real numbers from a given line
!
!---------------------------------------------------------------------------
integer function ncolumnsline(line)
 implicit none
 character(len=*), intent(in) :: line
 real :: dummyreal(1000)
 integer :: ierr,i

 dummyreal = -666666.0

 ierr = 0
 read(line,*,iostat=ierr) (dummyreal(i),i=1,size(dummyreal))

 i = 1
 ncolumnsline = 0
 do while(abs(dummyreal(i)+666666.).gt.tiny(0.))
    ncolumnsline = ncolumnsline + 1
    i = i + 1
    if (i.gt.size(dummyreal)) then
       print "(a)",'*** ERROR: too many columns in file'
       ncolumnsline = size(dummyreal)
       return
    endif
 enddo

end function ncolumnsline

!----------------------------------------------------------------------
!
! Small utility to return the number of comment lines in an ascii
! file. These are lines that do not begin with a number.
! 
! This is slightly different to what is done in the get_ncolumns
! routine, where header lines are any lines not having the same number
! of columns. Here we do not attempt to evaluate the number of data
! columns.
!
! File must be open and at the desired starting position
!----------------------------------------------------------------------
integer function nheaderlines(lunit)
 integer, intent(in) :: lunit
 real    :: dum
 integer :: ierr

 dum = -666.
 nheaderlines = 0
 ierr = -1
 do while (abs(dum+666.).lt.tiny(0.) .or. ierr.ne.0)
    nheaderlines = nheaderlines + 1
    read(lunit,*,iostat=ierr) dum
 enddo
 nheaderlines = nheaderlines - 1

end function nheaderlines

!---------------------------------------------------------------------------
!
! function stripping '/', '\' and spaces out of filenames
!
!---------------------------------------------------------------------------
function safename(string)
 implicit none
 character(len=*), intent(in) :: string
 character(len=len(string)) :: safename
 integer :: ipos

 safename = string

 !--remove forward slashes which can be mistaken for directories: replace with '_'
 call string_replace(safename,'/','_')
 call string_replace(safename,' ','_')
 
 !--delete brackets and operators of all kinds
 call string_delete(safename,'{')
 call string_delete(safename,'}')
 call string_delete(safename,'(')
 call string_delete(safename,')')
 call string_delete(safename,'[')
 call string_delete(safename,']')
 call string_delete(safename,'<')
 call string_delete(safename,'>')
 call string_delete(safename,'*')
 call string_delete(safename,'?')
 call string_delete(safename,'^')
 call string_delete(safename,'''')
 call string_delete(safename,'"')
 call string_delete(safename,'&')
 call string_delete(safename,'#')
 call string_delete(safename,'|')

 !--remove escape sequences: remove '\' and position following
 ipos = index(trim(safename),'\')
 do while (ipos.ne.0)
    safename = safename(1:ipos-1)//safename(ipos+2:len_trim(safename))
    ipos = index(trim(safename),'\')
 enddo

end function safename

!---------------------------------------------------------------------------
!
! function stripping the directory off a filename
!
!---------------------------------------------------------------------------
function basename(string)
 implicit none
 character(len=*), intent(in) :: string
 character(len=len(string)) :: basename
 integer :: i,iposmax

 basename = string

 !--find the last forward slash
 iposmax = 0
 i = len_trim(string)
 do while(i.ge.2 .and. iposmax.eq.0)
    i = i - 1
    if (string(i:i).eq.'/') iposmax = i
 enddo
 basename = trim(string(iposmax+1:))

end function basename

!---------------------------------------------------------------------------
!
! function to safely convert a string to c format (ie. with a terminating
! ascii null character)
!
!---------------------------------------------------------------------------
function cstring(string)
 implicit none
 character(len=*), intent(in) :: string
 character(len=len(string)+1) :: cstring

 cstring = trim(string)//achar(0)

end function cstring

!---------------------------------------------------------------------------
!
! function to safely convert a string from c format (ie. with a terminating
! ascii null character) back to a normal Fortran string
!
!---------------------------------------------------------------------------
function fstring(array)
 use, intrinsic :: iso_c_binding, only:c_char
 implicit none
 character(kind=c_char), dimension(:), intent(in) :: array
 character(len=size(array)-1) :: fstring
 integer :: i

 fstring = ''
 do i=1,size(array)
    if (array(i).eq.achar(0)) exit
    fstring(i:i) = array(i)
 enddo

end function fstring

!---------------------------------------------------------------------------
!
! subroutine to replace a matching section of a string with another
! string, possibly of differing length
!
!---------------------------------------------------------------------------
subroutine string_replace(string,skey,sreplacewith)
 implicit none
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skey,sreplacewith
 integer :: ipos,lensub

 ipos = index(trim(string),skey)
 lensub = len(skey)
 do while(ipos.gt.0)
    string = string(1:ipos-1)//sreplacewith//string(ipos+lensub:len_trim(string))
    ipos = index(trim(string),skey)
 enddo

end subroutine string_replace

!---------------------------------------------------------------------------
!
! subroutine to replace a specified section of a string with a
! replacement string, possibly of differing length
!
!---------------------------------------------------------------------------
subroutine string_sub(string,i1,i2,sreplacewith)
 implicit none
 character(len=*), intent(inout) :: string
 integer, intent(in)             :: i1,i2
 character(len=*), intent(in)    :: sreplacewith
 character(len=len(string))      :: oldstring

 oldstring = string
 if (i2 < len_trim(string)) then
    string = oldstring(1:i1-1)//sreplacewith//oldstring(i2+1:len_trim(oldstring))
 else
    string = oldstring(1:i1-1)//sreplacewith
 endif

end subroutine string_sub

!---------------------------------------------------------------------------
!
! subroutine to delete all matching occurrences of key from string
!
!---------------------------------------------------------------------------
pure subroutine string_delete(string,skey)
 implicit none
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skey
 integer :: ipos,lensub

 ipos = index(trim(string),skey)
 lensub = len(skey)
 do while(ipos.gt.0)
    string = string(1:ipos-1)//string(ipos+lensub:len_trim(string))
    ipos = index(trim(string),skey)
 enddo

end subroutine string_delete

!---------------------------------------------------------------------------
!
! Converts a string to upper case
!
!---------------------------------------------------------------------------
function ucase(string)
 implicit none
 character(len=*), intent(in) :: string
 character(len=len(string))   :: ucase
 integer :: is,ia
 integer, parameter           :: aoffset = 32

 ucase = string
 do is = 1, len(ucase)
    ia = iachar(ucase(is:is))
    if (ia >= iachar('a').and.ia <= iachar('z')) &
        ucase(is:is) = achar(ia-aoffset)
 enddo

end function ucase

!---------------------------------------------------------------------------
!
! Converts a string to lower case
!
!---------------------------------------------------------------------------
function lcase(string)
 implicit none
 character(len=*), intent(in) :: string
 character(len=len(string))   :: lcase
 integer :: is,ia
 integer, parameter           :: aoffset = 32

 lcase = string
 do is = 1, len(lcase)
    ia = iachar(lcase(is:is))
    if (ia >= iachar('A').and.ia <= iachar('Z')) &
        lcase(is:is) = achar(ia+aoffset)
 enddo

end function lcase

!---------------------------------------------------------------------------
!
! Converts a string to lower case
!
!---------------------------------------------------------------------------
integer function get_line_containing(filename,string)
 character(len=*), intent(in) :: filename, string
 character(len=130) :: line
 integer :: i,ierr
 integer, parameter :: lu=95
 
 get_line_containing = 0
 open(unit=lu,file=filename,status='old',iostat=ierr)
 i = 0
 do while(ierr.eq.0)
    i = i + 1
    read(lu,"(a)",iostat=ierr) line
    if (index(line,string).ne.0) get_line_containing = i
 enddo
 close(lu)
 
end function get_line_containing

end module asciiutils
