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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
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
 public :: read_asciifile,get_ncolumns,get_nrows,ncolumnsline,safename,basename,numfromfile
 public :: cstring,fstring,add_escape_chars
 public :: string_replace,string_delete,get_nheaderlines,string_sub
 public :: ucase,lcase,strip
 public :: get_line_containing
 public :: enumerate,isdigit,get_digits,integer_to_string,split
 public :: get_column_labels,read_column_labels
 public :: match_tag,match_taglist,append_number,make_tags_unique,get_value
 public :: match_column,match_tag_start,match_integer,match_lists
 public :: count_non_blank,find_repeated_tags,count_char
 public :: get_extensions,readline_csv,extension
 public :: sort_filenames_for_comparison
 public :: read_var_from_file
 integer, parameter :: max_line_length = 10000 ! for finding number of columns

 private

!--------------------------------------------------
! Generic interface to ascii file read for either
! character arrays (ie. each line is an element)
! or an array of real numbers
!--------------------------------------------------
 interface read_asciifile
  module procedure read_asciifile_char, read_asciifile_real,&
                   read_asciifile_real_string, read_asciifile_realarr, &
                   read_asciifile_int
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
subroutine read_asciifile_char(filename,nlinesread,charline,ierror,skip)
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 character(len=*), dimension(:), intent(out) :: charline
 integer, intent(out), optional :: ierror
 logical, intent(in), optional :: skip
 integer :: ierr,i,j,maxlines,iunit
 logical :: iexist,do_skip
 character(len=1) :: temp

 nlinesread = 0
 if (present(ierror)) ierror = 0

 !--if file does not exist, do nothing and return
 inquire(file=filename,exist=iexist)
 if (.not.iexist) then
    if (present(ierror)) ierror = -1
    return
 endif

 open(newunit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
 !--error opening file (but file does exist)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
    if (present(ierror)) ierror = ierr
    return
 endif

 ! read lines from file, skipping blank lines
 maxlines = size(charline)
 i = 0
 j = 1
 ierr = 0
 do_skip = .false.
 if (present(skip)) do_skip = skip
 over_lines: do while(j <= maxlines .and. ierr == 0)
    i = i + 1
    read(iunit,"(a)",iostat=ierr) charline(j)
    ! skip blank and comment lines
    if (ierr == 0) then
       if (do_skip) then
          temp = adjustl(charline(j))
          if (len_trim(charline(j)) > 0 .and. temp(1:1) /= '#') j = j + 1
       else
          j = j + 1
       endif
    endif
 enddo over_lines
 nlinesread = j-1

 ! emit warnings if errors or reached array limits
 if (nlinesread >= maxlines) then
    !--end of array limits
    !  check to see if there is anything more in the file. Report error if there is.
    read(iunit,"(a)",iostat=ierr)
    if (ierr==0) then
       print "(/,a,i6,/)",' WARNING: array limits reached reading '//trim(filename)//', max = ',maxlines
    endif
    nlinesread = min(maxlines,nlinesread-1)
 elseif (ierr > 0) then
    print "(a,i6)",' ERROR reading '//trim(filename)//' at line ',i
    if (present(ierror)) ierror = 1
 endif
 close(unit=iunit)

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
! returns array of real numbers (either one per line or all on same line)
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile_int(filename,nlinesread,intarr,ierror)
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 integer, dimension(:), intent(out) :: intarr
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

 intarr(:) = -666
 maxlines = size(intarr)
 read(iunit,*,err=66,end=99) (intarr(i),i=1,maxlines)

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
    if (abs(intarr(i)+666) > 0) nlinesread = nlinesread + 1
 enddo
 close(unit=iunit)
 return

 !--reached end of file (the expected behaviour)
99 continue
 do i=1,maxlines
    if (abs(intarr(i)+666) > 0) nlinesread = nlinesread + 1
 enddo
 close(unit=iunit)
 return

end subroutine read_asciifile_int

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
subroutine get_ncolumns(lunit,ncolumns,nheaderlines,csv,maxheaderlines)
 integer, intent(in) :: lunit
 integer, intent(out) :: ncolumns,nheaderlines
 integer, intent(in), optional :: maxheaderlines
 logical, intent(in), optional :: csv
 integer :: ierr,ncolprev,ncolprev2,ncolsthisline,maxlines,ncolstot
 character(len=max_line_length) :: line
 logical :: nansinfile,infsinfile,is_csv

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
 ncolprev2 = -200
 ncolsthisline = 0
 ncolstot = 0
 nansinfile = .false.
 infsinfile = .false.
 is_csv = .false.
 if (present(csv)) is_csv = csv
!
!--loop until we find two consecutive lines with the same number of columns (but non zero)
!  if ncolumns==1 then we must find 3 consecutive lines
!
 do while ((len_trim(line)==0 .or. ncolsthisline /= ncolprev .or. ncolumns < 1 .or. &
           (ncolumns==1 .and. ncolsthisline /= ncolprev2)) &
           .and. ierr==0 .and. nheaderlines <= maxlines)
    ncolprev2 = ncolprev
    ncolprev = ncolumns
    read(lunit,"(a)",iostat=ierr) line
    if (index(line,'NaN') > 0) nansinfile = .true.
    if (index(line,'Inf') > 0) infsinfile = .true.
    if (len_trim(line)==0) then
       ncolsthisline = -1
    else
       if (ierr==0) ncolsthisline = ncolumnsline(line,csv=is_csv,ntot=ncolstot)
       ncolumns = ncolsthisline
    endif
    nheaderlines = nheaderlines + 1
    !print*,'DEBUG: header line ',nheaderlines,' ncols = ',ncolsthisline,'"'//trim(line)//'"'
 enddo
 !--subtract 2 from the header line count (the last two lines which were the same)
 nheaderlines = max(nheaderlines - 2,0)
 if (ncolumns==1) nheaderlines = max(nheaderlines - 1,0)
 if (is_csv) ncolumns = ncolstot

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
integer function ncolumnsline(line,csv,ntot)
 character(len=*), intent(in)   :: line
 logical, intent(in),  optional :: csv
 integer, intent(out), optional :: ntot
 real :: dummyreal(1000)
 integer :: ierr,i
 logical :: use_commas

 use_commas= .false.
 if (present(csv)) use_commas = csv
 if (use_commas) then
    if (present(ntot)) then
       ncolumnsline = ncolumnsline_csv(line,ntot)
    else
       ncolumnsline = ncolumnsline_csv(line)
    endif
    return
 endif

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

!---------------------------------------------------------------------------
!
! function returning the number of columns of real numbers from a given line
!
!---------------------------------------------------------------------------
integer function ncolumnsline_csv(line,ntot) result(ncols)
 character(len=*), intent(in) :: line
 integer, parameter :: lenf = 15
 character(len=lenf) :: fields(len(line)/lenf)
 integer, intent(out), optional :: ntot
 integer :: i,ierr,nfields
 real :: dum

 ! split line by commas
 call split(line,',',fields,nfields)

 ! report how many columns contain real numbers
 ! or blank (non-text) entries
 ncols = 0
 do i=1,nfields
    if (len_trim(fields(i))==0) then
       ncols = ncols + 1
    else
       dum = -666666.
       read(fields(i),*,iostat=ierr) dum
       if (ierr==0) then
          ncols = ncols + 1
       endif
    endif
 enddo

 if (present(ntot)) ntot = nfields

end function ncolumnsline_csv

!---------------------------------------------------------------------------
!
! read a line from a csv file and parse for real numbers
!
!---------------------------------------------------------------------------
subroutine readline_csv(line,ncols,datcol)
 character(len=*), intent(in) :: line
 integer, intent(in)  :: ncols
 real,    intent(out) :: datcol(ncols)
 integer, parameter :: lenf = 15
 character(len=lenf) :: fields(ncols)
 !logical, intent(in)  :: mask(ncols)
 integer :: nfields,i,icol,ierr

 ! split line by commas
 call split(line,',',fields,nfields)

 ! read only columns that contain real numbers
 icol = 0
 do i=1,min(nfields,ncols)
    icol = icol + 1
    read(fields(i),*,iostat=ierr) datcol(icol)
 enddo

end subroutine readline_csv

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
integer function get_nheaderlines(lunit) result(nheaderlines)
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

end function get_nheaderlines

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
! function to get file extension
!
!---------------------------------------------------------------------------
function extension(string) result(ext)
 character(len=*), intent(in) :: string
 character(len=len(string))   :: ext
 integer :: idot

 idot = index(string,'.',back=.true.)
 if (idot > 2 .and. idot+1 <= len(string)) then
    ext = string(idot:)
 else
    ext = ''
 endif

end function extension

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
 integer :: ipos,imax,lensub,i,iposnext

 ipos = index(trim(string),skey)
 lensub = len(skey)
 imax   = len(string)
 i = 0
 do while (ipos > 0 .and. i <= imax)
    iposnext = ipos + len(sreplacewith)
    i = i + 1  !  only allow as many replacements as characters
    remstring = string(ipos+lensub:len_trim(string))
    string = string(1:ipos-1)//sreplacewith//remstring
    ipos = index(trim(string(iposnext:)),skey)
    if (ipos > 0) ipos = ipos + iposnext
    !print*,ipos,' string = ',trim(string),'iposnext = ',trim(string(iposnext:))
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

!------------------------------------------------------------------------
!     get_digits: for an integer i returns number of digits it contains
!     and a list of these *without* using write statements
!
!     i            : integer to split into digits
!     nmax           : dimensions of digits array
!     digits(nmax) : array of digits
!     ndigits      : number of digits in i
!------------------------------------------------------------------------
pure subroutine get_digits(i,digits,ndigits)
 integer, intent(in) :: i
 integer, intent(out) :: ndigits
 integer, intent(out), dimension(:) :: digits
 integer :: j,isubtract,idigit

 ndigits = 0

 isubtract = 0

 do j=size(digits),0,-1
    if (i >= 10**j) then
       ndigits = ndigits + 1
       idigit = (i - isubtract)/10**j
       digits(ndigits) = idigit
       isubtract = isubtract + digits(ndigits)*10**j
    endif
 enddo

end subroutine get_digits

!---------------------------------------------------------------------------
!
! convert a string to an integer WITHOUT using write statement
! (so this can be used in a write or print statement)
!
!---------------------------------------------------------------------------
pure elemental function integer_to_string(i) result(string)
 integer, intent(in) :: i
 character(len=12) :: string
 integer :: i0,ndigits,j
 integer :: idigit(12)

 string = ''
 i0 = iachar('0')
 call get_digits(i,idigit,ndigits)
 do j=2,ndigits
    string(j:j) = achar(i0+idigit(j))
 enddo

end function integer_to_string

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
 character(len=*), intent(out), dimension(:), optional :: stringarr
 integer,          intent(out) :: nsplit
 integer :: i,j,imax,iend,nmax

 i = 1
 nsplit = 0
 imax = len(string)
 nmax = imax
 if (present(stringarr)) nmax = size(stringarr)

 do while(nsplit < nmax .and. i <= imax)
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
    if (nsplit <= nmax .and. present(stringarr)) then
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
subroutine get_column_labels(line,nlabels,labels,method,ndesired,csv)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: nlabels
 character(len=*), dimension(:), intent(out) :: labels
 integer,          intent(out), optional :: method
 integer,          intent(in),  optional :: ndesired
 logical,          intent(in),  optional :: csv
 integer :: i1,i2,i,nlabelstmp,nlabels_prev,istyle,ntarget
 character(len=1) :: leadingchar
 character(len=4), parameter :: spaces = '    '
 logical :: is_csv

 nlabels = 0
 i1 = 1
 istyle = 0
 ntarget = -1
 is_csv = .false.
 if (present(csv)) is_csv = csv
 if (present(ndesired)) ntarget = ndesired
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

 if (index(nospaces(line),'][') > 0 .and. .not.is_csv) then
    !
    ! format style 1: # [ mylabel1 ] [ mylabel2 ] [ mylabel3 ]
    !
    istyle = 1
    i1 = max(index(line,'[')+1,i1)    ! strip leading square bracket
    ! try with different number of spaces between brackets (if labels not found)
    over_spaces1: do i=4,0,-1
       call split(line(i1:),']'//spaces(1:i)//'[',labels,nlabels)
       if (nlabels > 1) exit over_spaces1
    enddo over_spaces1
 elseif (index(line,',') > 1 .or. is_csv) then
    !
    ! format style 2: mylabel1,mylabel2,mylabel3
    !
    istyle = 2
    call split(line(i1:),',',labels,nlabelstmp)
    if (is_csv) then
       nlabels = nlabelstmp  ! allow blank/arbitrary labels in csv format
    else
       nlabels = count_sensible_labels(nlabelstmp,labels)
    endif
 else
    !
    ! format style 3: #     mylabel1     mylabel2     mylabel3
    !
    istyle = 3
    ! try splitting with 4, then 3, then 2 spaces until the number of labels decreases
    nlabels_prev = 0
    over_spaces: do i=4,2,-1
       call split(line(i1:),spaces(1:i),labels,nlabelstmp)
       ! quit if we already have the target number of labels
       if (nlabelstmp == ntarget) exit over_spaces

       ! if the number of labels is > 1 but has decreased, quit, unless nlabels
       ! still exceeds the number of labels we are hoping for (ntarget)
       if ((nlabelstmp < nlabels_prev .or. nlabelstmp >= max(nlabels_prev,2)  &
            .and. i < 4 .and. .not. (ntarget > 0 .and. nlabelstmp > ntarget))) then
          ! take the answer with the previous number of spaces
          call split(line(i1:),spaces(1:i+1),labels,nlabelstmp)
          exit over_spaces
       endif
       nlabels_prev = nlabelstmp
    enddo over_spaces
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
 ! delete loose trailing square bracket but only if not matching
 if (istyle==1) then
    if (index(labels(nlabels),']') > 0) then
       i1 = count_char(labels(nlabels),'[') ! number of open brackets
       i2 = count_char(labels(nlabels),']') ! number of closed brackets
       if (i2 > i1) then ! if brackets do not match
          ! find last trailing bracket
          i2 = index(labels(nlabels),']',back=.true.)
          ! delete it, but only if followed by spaces
          if (i2==len_trim(labels(nlabels))) then
             labels(nlabels) = labels(nlabels)(1:i2-1)
          endif
       endif
    endif
 endif

end subroutine get_column_labels

!---------------------------------------------------------------------------
!
! interface to the above routine that also searches for the line
! containing the column labels in the list of header lines
!
!---------------------------------------------------------------------------
subroutine read_column_labels(iunit,nheaderlines,ncols,nlabels,labels,csv,debug)
 integer,          intent(in)  :: iunit,nheaderlines,ncols
 integer,          intent(out) :: nlabels
 character(len=*), dimension(:), intent(out) :: labels
 logical, intent(in), optional :: csv,debug
 character(len=len(labels(1))), dimension(size(labels)) :: tmplabel
 character(len=max_line_length) :: line
 logical :: is_csv,verbose,got_labels
 integer :: i,imethod,ierr,nwanted

 is_csv = .false.
 verbose = .false.
 if (present(csv)) is_csv = csv
 if (present(debug)) verbose = debug
 got_labels = .false.
 nlabels = 0
 nwanted = min(ncols,size(labels)) ! can either retrieve all labels or completely fill the labels array
 labels = ''
 rewind(iunit)
 do i=1,nheaderlines
    read(iunit,"(a)",iostat=ierr) line
    !--try to match column labels from this header line, if not already matched (or dubious match)
    call get_column_labels(trim(line),nlabels,tmplabel,method=imethod,ndesired=nwanted,csv=csv)
    !--if we get nlabels > ncolumns, use them, but keep trying for a better match
    if ((got_labels .and. nlabels == nwanted) .or. &
        (.not.got_labels .and. nlabels >= nwanted  & ! only allow single-spaced labels if == ncols
         .and. (.not.(imethod>=4) .or. nlabels==nwanted))) then
       labels(1:nwanted) = tmplabel(1:nwanted)
       got_labels = .true.
    endif
    if (verbose) print "(5(1x,a,i0))",'DEBUG: line ',i,'nlabels = ',nlabels,&
                 'want ',ncols,'method=',imethod,'len_trim(line)=',len_trim(line) !,' LABELS= '//tmplabel(1:ncols)
 enddo

end subroutine read_column_labels

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
! match tag against a list of tags
! returns index of matching tag in the list
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

!--------------------------------------------
! as above but only match first N characters
! where N=5 by default
!--------------------------------------------
integer function match_tag_start(tags,tag,n)
 character(len=*), intent(in) :: tags(:)
 character(len=*), intent(in) :: tag
 integer, intent(in), optional :: n
 integer :: i,ilen
 character(len=len(tag)) :: str1,str2

 ilen = 5
 if (present(n)) ilen = n

 match_tag_start = 0 ! default if not found
 do i=1,size(tags)
    str1 = tags(i)(1:ilen)
    str2 = tag(1:ilen)
    if (trim(lcase(str2))==trim(lcase(str1))) then
       match_tag_start = i
       exit  ! only match first occurrence
    endif
 enddo

end function match_tag_start

!------------------------------------------
! match tag against a list of tags
! or by giving the column number explicitly
! returns index of matching tag in the list
!------------------------------------------
integer function match_column(tags,tag)
 character(len=*), intent(in) :: tags(:)
 character(len=*), intent(in) :: tag
 integer :: ierr

 ! try to match the string tag first
 match_column = match_tag(tags,tag)
 if (match_column == 0) then
    ! try to read it as an integer from the string
    read(tag,*,iostat=ierr) match_column
 endif

end function match_column

!----------------------------------------------
! match tag against a list of tags
! and extract the value from an array of reals
!----------------------------------------------
real function get_value(tag,tags,vals,default)
 character(len=*), intent(in) :: tag
 character(len=*), intent(in) :: tags(:)
 real, intent(in) :: vals(:)
 real, intent(in), optional :: default
 integer :: itag

 itag = match_tag(tags,tag)
 if (itag > 0 .and. itag <= size(vals)) then
    get_value = vals(itag)
 else
    if (present(default)) then
       get_value = default
    else
       get_value = 0.
    endif
 endif

end function get_value

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
       do j=2,min(size(taglist),size(tags)-i+1)
          if (trim(tags(i+j-1))==trim(taglist(j))) then
             nmatch = nmatch + 1
          endif
       enddo
    endif
 enddo

end subroutine match_taglist

!------------------------------------------
! find first integer that matches in a
! list of integers
!------------------------------------------
integer function match_integer(ivals,i)
 integer, intent(in) :: ivals(:)
 integer, intent(in) :: i
 integer :: k

 match_integer = 0
 do k=1,size(ivals)
    if (ivals(k)==i) then
       match_integer = k
       exit
    endif
 enddo

end function match_integer

!------------------------------------------
! find labels that match between two lists
! (match == first N characters are the same)
! output is a list of indices of labels
! from list1 that match list2
!------------------------------------------
function match_lists(list1,list2) result(imap)
 character(len=*), intent(in) :: list1(:),list2(:)
 integer :: imap(size(list1))
 integer :: j,icol

 do j=1,size(list1)
    icol = match_tag_start(list2,list1(j))
    if (icol > 0) imap(j) = icol
 enddo

end function match_lists

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

!-----------------------------------------------------------------
!
!  utility to count number of times a character appears in a string
!
!-----------------------------------------------------------------
integer function count_char(string,mychar)
 character(len=*), intent(in) :: string
 character(len=1), intent(in) :: mychar
 integer :: i

 count_char = 0
 do i=1,len(string)
    if (string(i:i)==mychar) count_char = count_char + 1
 enddo

end function count_char

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

!------------------------------------------------------------
! utility to return up to five file extensions
!------------------------------------------------------------
subroutine get_extensions(string,extensions)
 character(len=*), intent(in) :: string
 character(len=12), dimension(5), intent(out) :: extensions(5)
 character(:), allocatable :: tmp_string

 integer :: ppos_new
 integer :: ppos_old
 integer :: i

 ppos_new = scan(trim(string),".", BACK= .true.)
 ppos_old = len(string)
 tmp_string = lcase(string)

 do i=1,5
    if (ppos_new > 0) then
       extensions(i) = trim(tmp_string(ppos_new:ppos_old))
       tmp_string=tmp_string(1:ppos_new-1)
       ppos_old=ppos_new-1
       ppos_new=scan(trim(tmp_string),".",BACK=.true.)
    else
       extensions(i)=""
    endif
 enddo

end subroutine get_extensions

!---------------------------------------------------------------------------
!+
!  extract the start of the file extension, if the filename does not
!  end with digits
!+
!---------------------------------------------------------------------------
pure integer function get_idot(string)
 character(len=*), intent(in) :: string
 integer :: ilen

 ilen = len_trim(string)
 get_idot = 0
 !
 ! if file ends in at least two numbers then use the numbers at the end
 ! (two is to avoid problems with .hdf5 etc)
 !
 if (ilen >= 2) then
    if (isdigit(string(ilen:ilen)) .and. isdigit(string(ilen-1:ilen-1))) then
       get_idot = ilen + 1
    endif
 endif
 !
 ! otherwise, look for numbers before the file extension (e.g. _0000.dat)
 !
 if (get_idot==0) then
    get_idot = index(string,'.',back=.true.)
    if (get_idot==0) get_idot = len_trim(string) + 1
 endif

end function get_idot

!----------------------------------------------------------------
!+
!  this function extracts the number at the end of the filename
!+
!----------------------------------------------------------------
integer function numfromfile(filename)
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: string
 integer :: idot,istartnum,ilen,i,ierr
!
!--extract current number from filename
!
 string = basename(filename)
 idot = get_idot(string)
 istartnum = 0
 do i=idot-1,1,-1
    if (istartnum==0) then
       if (.not.isdigit(string(i:i))) istartnum = i
    endif
 enddo
 if (istartnum /= 0) istartnum = istartnum + 1
 ilen = idot - istartnum

 if (ilen > 0) then
    read(string(istartnum:istartnum+ilen-1),*,iostat=ierr) numfromfile
    if (ierr /= 0) then
       !print*,'internal error in numfromfilename'
       numfromfile = -1
    endif
 else
    numfromfile = 0
 endif

end function numfromfile

!------------------------------------------------------------
! utility to reorder a list of files. We use an insertion
! sort algorithm which preserves relative order of filenames
! so a list of files like:
!
!  disc_00000 disc_00001 disc_00002 disc1_00000 disc1_00002
!
! will be reordered as:
!
!  disc_00000 disc1_00000 disc_00001 disc1_00001 disc_00002
!
!------------------------------------------------------------
subroutine sort_filenames_for_comparison(nfiles,filenames)
 integer, intent(in) :: nfiles
 character(len=*), intent(inout) :: filenames(nfiles)
 character(len=100) :: key
 integer :: i,j,num

 do i = 2, nfiles
    key = filenames(i)
    num = numfromfile(key)
    j = i - 1

    ! move elements of the filenames array with number
    ! greater than the key to one position ahead 
    ! of their current position
    do while(j >= 1 .and. numfromfile(filenames(j)) > num)
       filenames(j+1) = filenames(j)
       j = j - 1
    enddo
    filenames(j+1) = key
 enddo

end subroutine sort_filenames_for_comparison

!------------------------------------------------------------
! utility to read a variable from an ascii file
! in the form:
!
!   var = val   ! comment
!
! val is returned as a string, if not found leaves the
! input value unmodified
!------------------------------------------------------------
subroutine read_var_from_file(var,val,filename,ierr)
 character(len=*), intent(in) :: var
 character(len=*), intent(inout) :: val
 character(len=*), intent(in) :: filename
 integer, intent(out) :: ierr
 character(len=130) :: line
 integer :: j,lu,ieq,istart
 logical :: match

 open(newunit=lu,file=filename,status='old',iostat=ierr)
 match = .false.
 do while(ierr==0 .and. .not.match)
    read(lu,"(a)",iostat=ierr) line
    if (index(line,var) /= 0) then
       istart = index(line,var)
       ieq = index(line(istart:),'=')
       val = line(istart+ieq+1:)
       match = .true.
    endif
 enddo
 close(lu)

 if (match) then
    val = trim(adjustl(val))
    ! now cull the variable at a space or comma
    do j=1,len_trim(val)
       if (any((/' ',',',';',':','=','!'/)==val(j:j))) then
          val = val(1:j)
          exit
       endif
    enddo
 !else
 !print*,trim(var)//' not found in '//trim(filename)
 endif

end subroutine read_var_from_file

end module asciiutils
