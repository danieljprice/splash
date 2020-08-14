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
!  Copyright (C) 2005-2018 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
! module containing various utility subroutines
! related to reading from ascii files and dealing with string variables
!
! written by Daniel Price, University of Exeter 2007 24th April '07
! revised at Monash University, 2008 -
! daniel.price@monash.edu
!
! this is a standalone module with no dependencies
!---------------------------------------------------------------------------
module asciiutils
 implicit none
 public :: read_asciifile,get_ncolumns,get_nrows,ncolumnsline,safename,basename
 public :: cstring,fstring,add_escape_chars
 public :: string_replace, string_delete, nheaderlines, string_sub
 public :: ucase,lcase,strip
 public :: get_line_containing
 public :: enumerate,isdigit,split
 public :: get_column_labels
 public :: match_tag,match_taglist,append_number,make_tags_unique
 public :: count_non_blank,find_repeated_tags

 private

!--------------------------------------------------
! Generic interface to ascii file read for either
! character arrays (ie. each line is an element)
! or an array of real numbers
!--------------------------------------------------
 interface read_asciifile
  module procedure read_asciifile_char, read_asciifile_real,&
                    read_asciifile_real_string, read_asciifile_realarr
 end interface read_asciifile

 interface string_delete
  module procedure string_delete1,string_delete_array
 end interface string_delete

contains

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns array of character strings (one per line)
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile_char(filename,nlinesread,charline,ierror)
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
 if (ierr==0) then
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
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 real, dimension(:), intent(out) :: realarr
 integer, intent(out), optional :: ierror
 integer, parameter :: iunit = 66 ! logical unit number for read operation
 integer :: ierr,i,maxlines
 logical :: iexist

 i = 0
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
    if (abs(realarr(i)+666.) > tiny(0.)) nlinesread = nlinesread + 1
 enddo
 close(unit=iunit)
 return

 !--reached end of file (the expected behaviour)
99 continue
 do i=1,maxlines
    if (abs(realarr(i)+666.) > tiny(0.)) nlinesread = nlinesread + 1
 enddo
 close(unit=iunit)
 return

end subroutine read_asciifile_real

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns 2D array of real numbers (i.e. tabulated data)
!---------------------------------------------------------------------------
subroutine read_asciifile_realarr(filename,nlinesread,realarr,ierror)
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 real, dimension(:,:), intent(out) :: realarr
 integer, intent(out), optional :: ierror
 integer, parameter :: iunit = 66 ! logical unit number for read operation
 integer :: ierr,i,ncols,ncolsfile,nheader
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
 else
    ! get number of columns
    call get_ncolumns(iunit,ncolsfile,nheader)
    ! skip header lines
    do i=1,nheader
       read(iunit,*,iostat=ierr)
    enddo
    ! read 2D array from file
    ncols = min(ncolsfile,size(realarr(:,1)))
    nlinesread = 0
    do while (ierr==0)
       nlinesread = nlinesread + 1
       read(iunit,*,iostat=ierr) realarr(1:ncols,nlinesread)
    enddo
    nlinesread = max(nlinesread - 1,0)
    close(iunit)
 endif

end subroutine read_asciifile_realarr

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns array of real numbers and corresponding string
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile_real_string(filename,nlinesread,realarr,charline,ierror)
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
    if (abs(realarr(i)+666.) > tiny(0.)) nlinesread = nlinesread + 1
 enddo
 close(unit=iunit)
 return

 !--reached end of file (the expected behaviour)
99 continue
 do i=1,maxlines
    if (abs(realarr(i)+666.) > tiny(0.)) nlinesread = nlinesread + 1
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
subroutine get_ncolumns(lunit,ncolumns,nheaderlines,maxheaderlines)
 integer, intent(in) :: lunit
 integer, intent(out) :: ncolumns,nheaderlines
 integer, intent(in), optional :: maxheaderlines
 integer :: ierr,ncolprev,ncolsthisline,maxlines
 character(len=5000) :: line
 logical :: nansinfile,infsinfile

 if (present(maxheaderlines)) then
    maxlines = maxheaderlines
 else
    maxlines = 1000
 endif
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
 do while ((len_trim(line)==0 .or. ncolsthisline /= ncolprev .or. ncolumns <= 0) &
           .and. ierr==0 .and. nheaderlines <= maxlines)
    ncolprev = ncolumns
    read(lunit,"(a)",iostat=ierr) line
    if (index(line,'NaN') > 0) nansinfile = .true.
    if (index(line,'Inf') > 0) infsinfile = .true.
    if (len_trim(line)==0) then
       ncolsthisline = -1
    else
       if (ierr==0) ncolsthisline = ncolumnsline(line)
       ncolumns = ncolsthisline
    endif
    nheaderlines = nheaderlines + 1
    !print*,'DEBUG: header line ',nheaderlines,' ncols = ',ncolsthisline,'"'//trim(line)//'"'
 enddo
 !--subtract 2 from the header line count (the last two lines which were the same)
 nheaderlines = max(nheaderlines - 2,0)
 if (ierr  > 0 .or. ncolumns <= 0) then
    ncolumns = 0
 elseif (ierr  <  0) then
    !print*,ncolumns,ncolprev
 endif
 if (nansinfile) print "(a)",' INDIAN BREAD WARNING!! NaNs in file!!'
 if (infsinfile) print "(a)",' WARNING!! Infs in file!!'
 rewind(lunit)

 if (ncolumns==0) print "(a)",' ERROR: no columns of real numbers found'

end subroutine get_ncolumns

!---------------------------------------------------------------------------
! utility to work out number of rows in file
!---------------------------------------------------------------------------
subroutine get_nrows(lunit,nheaderlines,nlines)
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
 character(len=*), intent(in) :: line
 real :: dummyreal(1000)
 integer :: ierr,i

 dummyreal = -666666.0

 ierr = 0
 read(line,*,iostat=ierr) (dummyreal(i),i=1,size(dummyreal))

 i = 1
 ncolumnsline = 0
 do while(abs(dummyreal(i)+666666.) > tiny(0.) .or. dummyreal(i) /= dummyreal(i))
    ncolumnsline = ncolumnsline + 1
    i = i + 1
    if (i > size(dummyreal)) then
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
 do while (abs(dum+666.) < tiny(0.) .or. ierr /= 0)
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
 do while (ipos /= 0)
    safename = safename(1:ipos-1)//safename(ipos+2:len_trim(safename))
    ipos = index(trim(safename),'\')
 enddo

end function safename

!---------------------------------------------------------------------------
!
! function to insert escape characters so filenames appear correctly in legend
!
!---------------------------------------------------------------------------
function add_escape_chars(string)
 character(len=*), intent(in) :: string
 character(len=len(string)) :: add_escape_chars

 add_escape_chars = string
 call string_replace(add_escape_chars,'_','\_')
 call string_replace(add_escape_chars,'^','\^')

end function add_escape_chars

!---------------------------------------------------------------------------
!
! function to strip spaces out of a string
!
!---------------------------------------------------------------------------
function nospaces(string)
 character(len=*), intent(in) :: string
 character(len=len(string)) :: nospaces

 nospaces = string
 call string_delete(nospaces,' ')

end function nospaces

!---------------------------------------------------------------------------
!
! function stripping the directory off a filename
!
!---------------------------------------------------------------------------
function basename(string)
 character(len=*), intent(in) :: string
 character(len=len(string)) :: basename
 integer :: i,iposmax

 basename = string

 !--find the last forward slash
 iposmax = 0
 i = len_trim(string)
 do while(i >= 2 .and. iposmax==0)
    i = i - 1
    if (string(i:i)=='/') iposmax = i
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
 character(kind=c_char), dimension(:), intent(in) :: array
 character(len=size(array)-1) :: fstring
 integer :: i

 fstring = ''
 do i=1,size(array)
    if (array(i)==achar(0)) exit
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
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skey,sreplacewith
 character(len=len(string)) :: remstring
 integer :: ipos,ioffset,lensub

 ipos = index(trim(string),skey)
 lensub = len(skey)
 do while(ipos > 0)
    remstring = string(ipos+lensub:len_trim(string))
    ioffset = ipos - 1 + len(sreplacewith)
    string = string(1:ipos-1)//sreplacewith//remstring
    ipos = index(trim(remstring),skey)
    if (ipos > 0) ipos = ipos + ioffset
 enddo

end subroutine string_replace

!---------------------------------------------------------------------------
!
! subroutine to replace a specified section of a string with a
! replacement string, possibly of differing length
!
!---------------------------------------------------------------------------
subroutine string_sub(string,i1,i2,sreplacewith)
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
pure subroutine string_delete1(string,skey)
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skey
 integer :: ipos,lensub

 ipos = index(string,skey)
 lensub = len(skey)
 do while(ipos > 0)
    string = string(1:ipos-1)//string(ipos+lensub:len_trim(string))
    ipos = index(trim(string),skey)
 enddo

end subroutine string_delete1

!---------------------------------------------------------------------------
!
! alternate version of above for arrays of keys
!
!---------------------------------------------------------------------------
pure subroutine string_delete_array(string,skeys)
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skeys(:)
 integer :: i

 do i=1,size(skeys)
    call string_delete(string,skeys(i))
 enddo

end subroutine string_delete_array
!---------------------------------------------------------------------------
!
! function version of string_delete
!
!---------------------------------------------------------------------------
pure elemental function strip(string,skey)
 character(len=*), intent(in) :: string,skey
 character(len=len(string))   :: strip

 strip = string
 call string_delete(strip,skey)

end function strip

!---------------------------------------------------------------------------
!
! Converts a string to upper case
!
!---------------------------------------------------------------------------
pure elemental function ucase(string)
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
pure elemental function lcase(string)
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
! indicate if a character is a digit (number) or not
!
!---------------------------------------------------------------------------
pure elemental logical function isdigit(string)
 character(len=1), intent(in) :: string
 integer :: ia

 isdigit = .false.
 ia = iachar(string)
 if (ia >= iachar('0').and.ia <= iachar('9')) isdigit = .true.

end function isdigit

!---------------------------------------------------------------------------
!
! search a file for the line containing a particular string
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
 do while(ierr==0)
    i = i + 1
    read(lu,"(a)",iostat=ierr) line
    if (index(line,string) /= 0) get_line_containing = i
 enddo
 close(lu)

end function get_line_containing

!---------------------------------------------------------------------------
!
! Convert an integer into the corresponding entry in a list of strings
!
!---------------------------------------------------------------------------
function enumerate(i,stringarr,default) result(string)
 integer, intent(in) :: i
 character(len=*), intent(in), dimension(:) :: stringarr
 integer, intent(in), optional :: default
 character(len=len(stringarr)) :: string

 string = ''
 if (i >= 1 .and. i <= size(stringarr)) then
    string = trim(stringarr(i))
 elseif (present(default)) then
    if (default >= 1 .and. i <= size(stringarr)) then
       string = trim(stringarr(default))
    endif
 endif

end function enumerate

!---------------------------------------------------------------------------
!
! Split a string into substrings based on a delimiter
!
!---------------------------------------------------------------------------
pure subroutine split(string,delim,stringarr,nsplit)
 character(len=*), intent(in)  :: string
 character(len=*), intent(in)  :: delim
 character(len=*), intent(out), dimension(:) :: stringarr
 integer,          intent(out) :: nsplit
 integer :: i,j,imax,iend

 i = 1
 nsplit = 0
 imax = len(string)
 do while(nsplit < size(stringarr) .and. i <= imax)
    ! find next non-blank character
    if (string(i:i)==' ') then
       do while (string(i:i)==' ')
          i = i + 1
          if (i > imax) exit
       enddo
       if (i > imax) exit
    endif

    ! look for next occurrence of delimiter
    j = index(string(i:),delim) - 1
    ! if no delimiter found, use whole rest of string
    if (j < 0) j = imax
    ! set end of substring
    iend = min(i+j-1,imax)
    ! extract the substring
    nsplit = nsplit + 1
    if (nsplit <= size(stringarr)) then
       stringarr(nsplit) = string(i:iend)
    endif
    i = iend + len(delim) + 1
 enddo

end subroutine split

!---------------------------------------------------------------------------
!
! extract a list of labels from the header line of a file
!
!---------------------------------------------------------------------------
subroutine get_column_labels(line,nlabels,labels,method)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: nlabels
 character(len=*), dimension(:), intent(out) :: labels
 integer,          intent(out), optional :: method
 integer :: i1,i2,i,nlabelstmp,istyle
 character(len=1) :: leadingchar

 nlabels = 0
 i1 = 1
 istyle = 0
 !
 ! strip leading comment character ('#')
 !
 leadingchar = trim(adjustl(line))
 if (leadingchar=='#') then
    i1 = index(line,'#') + 1
 endif
 ! strip anything preceding an equals sign
 i1 = max(i1,index(line,'=')+1)
 i2 = i1

 if (index(nospaces(line),'][') > 0) then
    !
    ! format style 1: # [ mylabel1 ] [ mylabel2 ] [ mylabel3 ]
    !
    istyle = 1
    call split(line(i1:),']',labels,nlabels)
 elseif (index(line,',') > 1) then
    !
    ! format style 2: mylabel1,mylabel2,mylabel3
    !
    istyle = 2
    call split(line(i1:),',',labels,nlabelstmp)
    nlabels = count_sensible_labels(nlabelstmp,labels)
 else
    !
    ! format style 3: #     mylabel1     mylabel2     mylabel3
    !
    istyle = 3
    call split(line(i1:),'  ',labels,nlabelstmp)
    !
    ! this style is dangerous, so perform sanity checks
    ! on the labels to ensure they are sensible
    !
    nlabels = count_sensible_labels(nlabelstmp,labels)
    if (nlabels <= 1) then
       !
       ! format style 4: x y z vx vy vz
       ! (this style is also dangerous)
       !
       istyle = 4
       call split(line(i1:),' ',labels,nlabelstmp)
       nlabels = count_sensible_labels(nlabelstmp,labels)
    endif
 endif
 if (present(method)) method = istyle
 !
 ! clean up
 !
 do i=1,nlabels
    ! delete brackets
    if (nlabels <= size(labels)) then
       call string_delete(labels(i),',')
       if (istyle==1) then
          call string_delete(labels(i),'[')
          call string_delete(labels(i),']')
       endif
       if (istyle==1 .or. istyle==2) then
          labels(i) = trim(adjustl(labels(i)))
          ! delete leading numbers
          i1 = 1
          do while (isdigit(labels(i)(i1:i1)))
             labels(i)(i1:i1) = ' '
             i1 = i1 + 1
          enddo
       endif
       labels(i) = trim(adjustl(labels(i)))
    endif
 enddo

end subroutine get_column_labels

!---------------------------------------------------------------------------
!
! count the number of sensible labels in a list of possible labels
!
!---------------------------------------------------------------------------
integer function count_sensible_labels(n,labels) result(m)
 integer, intent(in) :: n
 character(len=*), dimension(n), intent(in) :: labels
 integer :: i

 m = 0
 do i=1,n
    if (is_sensible_label(labels(i))) m = m + 1
 enddo

end function count_sensible_labels

!---------------------------------------------------------------------------
!
! determine if a particular string makes sense as a column label or not
!
!---------------------------------------------------------------------------
logical function is_sensible_label(string)
 character(len=*), intent(in) :: string
 real    :: dum
 integer :: ierr
 real, parameter :: dum_prev =  -66666666.

 is_sensible_label = .true.

 ! should not start with a decimal point
 if (string(1:1)=='.') is_sensible_label = .false.

 ! should not contain equals sign
 !if (index(string,'=') > 0) is_sensible_label = .false.

 dum = dum_prev
 ! should not be able to read it as a real number
 read(string,*,iostat=ierr) dum
 if (ierr==0 .and. abs(dum-dum_prev) > tiny(dum)) is_sensible_label = .false.

end function is_sensible_label

!------------------------------------------
! match tag against a list of string_sub
!------------------------------------------
integer function match_tag(tags,tag)
 character(len=*), intent(in) :: tags(:)
 character(len=*), intent(in) :: tag
 integer :: i

 match_tag = 0 ! default if not found
 do i=1,size(tags)
    if (trim(tags(i))==trim(adjustl(tag))) then
       match_tag = i
       exit  ! only match first occurrence
    endif
 enddo

end function match_tag

!-----------------------------------------------
! match multiple tags against a list of strings
! e.g. find 'x','y','z' in list of labels
!-----------------------------------------------
subroutine match_taglist(taglist,tags,istartmatch,nmatch)
 character(len=*), intent(in)  :: taglist(:)
 character(len=*), intent(in)  :: tags(:)
 integer,          intent(out) :: istartmatch,nmatch
 integer :: i,j

 istartmatch = 0
 nmatch = 0
 if (size(taglist) < 1) return
 do i=1,size(tags)
    if (nmatch <= 1 .and. trim(tags(i))==trim(taglist(1))) then
       nmatch = 1
       istartmatch = i
       do j=2,size(taglist)
          if (trim(tags(i+j-1))==trim(taglist(j))) then
             nmatch = nmatch + 1
          endif
       enddo
    endif
 enddo

end subroutine match_taglist

!------------------------------
! Append a number to a string
! e.g. string,2 -> string2
!------------------------------
subroutine append_number(string,j)
 character(len=*), intent(inout) :: string
 integer,          intent(in)    :: j
 character(len=12) :: strj

 write(strj,"(i12)") j
 string = trim(string)//trim(adjustl(strj))

end subroutine append_number

!----------------------------------------------------------------------
! Append numbers to otherwise identical tags to make them unique
! e.g. massoftype1, massoftype2, massoftype3, etc.
!----------------------------------------------------------------------
subroutine make_tags_unique(ntags,tags)
 integer, intent(in) :: ntags
 character(len=*), dimension(ntags), intent(inout) :: tags
 character(len=len(tags)) :: tagprev
 integer :: i,j

 if (ntags < 1) return
 j = 0
 tagprev = tags(1)
 do i=2,ntags
    if (tags(i)==tagprev) then
       j = j + 1
       if (j==1) then
          call append_number(tags(i-1),j)
          j = j + 1
       endif
       call append_number(tags(i),j)
    else
       tagprev = tags(i)
       j = 0
    endif
 enddo

end subroutine make_tags_unique

!-----------------------------------------------------------------
!
!  utility to count number of non-blank strings in a list
!
!-----------------------------------------------------------------
integer function count_non_blank(string)
 character(len=*), dimension(:), intent(in) :: string
 integer :: i

 count_non_blank = 0
 do i=1,size(string)
    if (len_trim(string(i))<=0) exit
    count_non_blank = count_non_blank + 1
 enddo

end function count_non_blank

!---------------------------------------------------------------------
!  utility to identify repeated tags in a list, written
!  for use when labels were made by "make_tags_unique"
!
!  IN:
!    tag - string to match, e.g. "massoftype"
!    ntags - length of array to search
!    tags  - sequence of character strings to search
!
!  the routine then searches for consecutive entries matching
!  the tag, e.g. massoftype1, massoftype2 etc and returns:
!
!  OUT:
!    istartlist - start of list (e.g. location of massoftype1)
!    nlist - number of tags matching string (i.e. N in massoftypeN)
!----------------------------------------------------------------------
subroutine find_repeated_tags(tag,ntags,tags,istartlist,nlist)
 character(len=*), intent(in)  :: tag
 integer,          intent(in)  :: ntags
 character(len=*), intent(in)  :: tags(ntags)
 integer,          intent(out) :: istartlist,nlist
 integer :: i
 logical :: consecutive

 istartlist = 0
 nlist = 0
 consecutive = .false.
 do i=1,ntags
    if (trim(tag)==tags(i)(1:len_trim(tag))) then
       if (nlist==0) then
          istartlist = i
          consecutive = .true.
       endif
       if (consecutive) nlist = nlist + 1
    else
       consecutive = .false.
    endif
 enddo

end subroutine find_repeated_tags

end module asciiutils
