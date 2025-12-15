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
!
!  Helper utilities for reading simple JSON files in SPLASH
!
!-----------------------------------------------------------------
module json_utils
 use iso_fortran_env, only:int64
 use params,          only:doub_prec
 implicit none
 private

 integer, parameter, public :: json_kind_unknown = 0
 integer, parameter, public :: json_kind_object  = 1
 integer, parameter, public :: json_kind_array   = 2
 integer, parameter, public :: json_kind_string  = 3
 integer, parameter, public :: json_kind_number  = 4
 integer, parameter, public :: json_kind_boolean = 5
 integer, parameter, public :: json_kind_null    = 6

 integer, parameter, public :: json_success      = 0
 integer, parameter, public :: json_err_not_found = 1
 integer, parameter, public :: json_err_parse     = 2
 integer, parameter, public :: json_err_type      = 3
 integer, parameter, public :: json_err_io        = 4

 type :: json_path_token
    character(:), allocatable :: key
    logical :: has_index = .false.
    integer :: index = -1
 end type json_path_token

 public :: json_read_file
 public :: json_get_raw_value
 public :: json_get_string_or_array
 public :: json_get_real
 public :: json_get_integer
 public :: json_get_logical
 public :: json_get_value_by_path
 public :: json_array_get_element
 public :: json_array_length
 public :: json_read

 interface json_read
  module procedure json_get_string_or_array, json_get_real, json_get_integer, json_get_logical,&
                   json_get_value_by_path, json_get_array_int64, json_get_array_int
 end interface json_read

contains

!-----------------------------------------------------------------
! read a JSON file into a string
!-----------------------------------------------------------------
subroutine json_read_file(filename, content, ierr)
 character(len=*), intent(in)  :: filename
 character(:), allocatable, intent(out) :: content
 integer, intent(out) :: ierr
 integer :: iu, ios
 logical :: exists
 integer(int64) :: file_size

 inquire(file=filename, exist=exists)
 if (.not. exists) then
    ierr = json_err_io
    allocate(character(len=0) :: content)
    return
 endif

 open(newunit=iu,file=trim(filename),status='old',action='read',access='stream', &
      form='unformatted',iostat=ios)
 if (ios /= 0) then
    ierr = json_err_io
    allocate(character(len=0) :: content)
    return
 endif

 inquire(unit=iu, size=file_size, iostat=ios)
 if (ios /= 0 .or. file_size < 0_int64) then
    close(iu)
    ierr = json_err_io
    allocate(character(len=0) :: content)
    return
 endif

 if (file_size == 0_int64) then
    allocate(character(len=0) :: content)
 else
    allocate(character(len=int(file_size)) :: content)
    read(iu,iostat=ios) content
    if (ios /= 0) then
       close(iu)
       ierr = json_err_io
       deallocate(content)
       allocate(character(len=0) :: content)
       return
    endif
 endif

 close(iu)
 ierr = json_success

end subroutine json_read_file

!-----------------------------------------------------------------
! get a string value from a JSON text
!-----------------------------------------------------------------
subroutine json_get_string_or_array(json_text, key, value, ierr)
 character(len=*), intent(in)  :: json_text
 character(len=*), intent(in)  :: key
 character(:), allocatable, intent(out) :: value
 integer, intent(out) :: ierr
 integer :: kind

 call json_get_raw_value(json_text,key,value,kind,ierr)
 if (ierr /= json_success) return

 if (kind /= json_kind_string .and. kind /= json_kind_array) then
    if (allocated(value)) deallocate(value)
    ierr = json_err_type
    allocate(character(len=0) :: value)
 else
    ierr = json_success
 endif

end subroutine json_get_string_or_array

!-----------------------------------------------------------------
! get a real value from a JSON text
!-----------------------------------------------------------------
subroutine json_get_real(json_text, key, value, ierr)
 character(len=*), intent(in) :: json_text
 character(len=*), intent(in) :: key
 real(kind=doub_prec), intent(out) :: value
 integer, intent(out) :: ierr
 character(:), allocatable :: raw
 integer :: kind

 call json_get_raw_value(json_text,key,raw,kind,ierr)
 if (ierr /= json_success) return

 if (kind /= json_kind_number) then
    ierr = json_err_type
    value = 0.0d0
    if (allocated(raw)) deallocate(raw)
    return
 endif

 read(raw,*,iostat=ierr) value
 if (ierr /= 0) then
    ierr = json_err_parse
    value = 0.0d0
 else
    ierr = json_success
 endif
 if (allocated(raw)) deallocate(raw)

end subroutine json_get_real

!-----------------------------------------------------------------
! get an integer value from a JSON text
!-----------------------------------------------------------------
subroutine json_get_integer(json_text, key, value, ierr)
 character(len=*), intent(in) :: json_text
 character(len=*), intent(in) :: key
 integer, intent(out) :: value
 integer, intent(out) :: ierr
 character(:), allocatable :: raw
 integer :: kind

 call json_get_raw_value(json_text,key,raw,kind,ierr)
 if (ierr /= json_success) return

 if (kind /= json_kind_number) then
    ierr = json_err_type
    value = 0
    if (allocated(raw)) deallocate(raw)
    return
 endif

 read(raw,*,iostat=ierr) value
 if (ierr /= 0) then
    ierr = json_err_parse
    value = 0
 else
    ierr = json_success
 endif
 if (allocated(raw)) deallocate(raw)

end subroutine json_get_integer

!-----------------------------------------------------------------
! get an integer(kind=int64) array from a JSON text
!-----------------------------------------------------------------
subroutine json_get_array_int64(json_text, key, values, ierr)
 character(len=*), intent(in) :: json_text
 character(len=*), intent(in) :: key
 integer(kind=int64), allocatable, intent(out) :: values(:)
 integer, intent(out) :: ierr
 character(:), allocatable :: raw
 integer :: kind, parse_err

 ierr = json_success
 call json_get_raw_value(json_text, key, raw, kind, ierr)
 if (ierr /= json_success) then
    allocate(values(0))
    return
 endif

 if (kind /= json_kind_array) then
    ierr = json_err_type
    allocate(values(0))
    if (allocated(raw)) deallocate(raw)
    return
 endif

 call parse_int64_array(raw, values, parse_err)
 if (allocated(raw)) deallocate(raw)
 if (parse_err /= 0) then
    ierr = json_err_parse
    if (.not. allocated(values)) allocate(values(0))
 endif

end subroutine json_get_array_int64

!-----------------------------------------------------------------
! get an integer array from a JSON text
!-----------------------------------------------------------------
subroutine json_get_array_int(json_text, key, values, ierr)
 character(len=*), intent(in) :: json_text
 character(len=*), intent(in) :: key
 integer, allocatable, intent(out) :: values(:)
 integer, intent(out) :: ierr
 character(:), allocatable :: raw
 integer :: kind, parse_err

 ierr = json_success
 call json_get_raw_value(json_text, key, raw, kind, ierr)
 if (ierr /= json_success) then
    allocate(values(0))
    return
 endif

 if (kind /= json_kind_array) then
    ierr = json_err_type
    allocate(values(0))
    if (allocated(raw)) deallocate(raw)
    return
 endif

 call parse_int_array(raw, values, parse_err)
 if (allocated(raw)) deallocate(raw)
 if (parse_err /= 0) then
    ierr = json_err_parse
    if (.not. allocated(values)) allocate(values(0))
 endif

end subroutine json_get_array_int

!-----------------------------------------------------------------
! get a logical value from a JSON text
!-----------------------------------------------------------------
subroutine json_get_logical(json_text, key, value, ierr)
 character(len=*), intent(in) :: json_text
 character(len=*), intent(in) :: key
 logical, intent(out) :: value
 integer, intent(out) :: ierr
 character(:), allocatable :: raw
 integer :: kind

 call json_get_raw_value(json_text,key,raw,kind,ierr)
 if (ierr /= json_success) return

 if (kind /= json_kind_boolean) then
    ierr = json_err_type
    value = .false.
    if (allocated(raw)) deallocate(raw)
    return
 endif

 select case (trim(raw))
 case('true')
    value = .true.
    ierr = json_success
 case('false')
    value = .false.
    ierr = json_success
 case default
    value = .false.
    ierr = json_err_parse
 end select
 if (allocated(raw)) deallocate(raw)

end subroutine json_get_logical

!-----------------------------------------------------------------
! get a value from a JSON text by specified path
!-----------------------------------------------------------------
subroutine json_get_value_by_path(json_text, path, value_text, kind, ierr)
 character(len=*), intent(in)  :: json_text
 character(len=*), intent(in)  :: path
 character(:), allocatable, intent(out) :: value_text
 integer, intent(out) :: kind
 integer, intent(out) :: ierr
 type(json_path_token), allocatable :: tokens(:)
 character(:), allocatable :: current
 character(:), allocatable :: raw
 integer :: current_kind
 integer :: i, ntokens

 ierr = json_success
 current_kind = detect_value_kind(json_text)
 current = json_text

 call parse_path(path,tokens,ntokens,ierr)
 if (ierr /= json_success) then
    if (allocated(current)) deallocate(current)
    allocate(character(len=0) :: value_text)
    kind = json_kind_unknown
    return
 endif

 do i = 1, ntokens
    if (len_trim(tokens(i)%key) > 0) then
       call json_get_raw_value(current,trim(tokens(i)%key),raw,kind,ierr)
       if (ierr /= json_success) exit
    else
       raw = current
       kind = current_kind
    endif

    if (tokens(i)%has_index) then
       call json_array_get_element(raw,tokens(i)%index,current,current_kind,ierr)
       if (allocated(raw)) deallocate(raw)
       if (ierr /= json_success) exit
    else
       if (allocated(current)) deallocate(current)
       call move_alloc(raw,current)
       current_kind = kind
    endif
 enddo

 if (allocated(raw)) then
    deallocate(raw)
 end if

 if (ierr == json_success) then
    kind = current_kind
    call move_alloc(current,value_text)
 else
    if (allocated(current)) deallocate(current)
    allocate(character(len=0) :: value_text)
 end if

 if (allocated(tokens)) deallocate(tokens)

end subroutine json_get_value_by_path

!-----------------------------------------------------------------
! get an element from a JSON array by giving the index
!-----------------------------------------------------------------
subroutine json_array_get_element(array_text, index, value_text, kind, ierr)
 character(len=*), intent(in)  :: array_text
 integer, intent(in) :: index
 character(:), allocatable, intent(out) :: value_text
 integer, intent(out) :: kind
 integer, intent(out) :: ierr
 integer :: len_arr, i, depth, current_index
 integer :: start_pos, end_pos
 logical :: in_string, escaped
 character :: ch

 ierr = json_success
 kind = json_kind_unknown
 if (index < 0) then
    ierr = json_err_not_found
    allocate(character(len=0) :: value_text)
    return
 endif

 len_arr = len_trim(array_text)
 if (len_arr < 2) then
    ierr = json_err_parse
    allocate(character(len=0) :: value_text)
    return
 endif
 if (array_text(1:1) /= '[') then
    ierr = json_err_type
    allocate(character(len=0) :: value_text)
    return
 endif

 in_string = .false.
 escaped = .false.
 depth = 0
 current_index = -1
 start_pos = 0
 end_pos = 0

 do i = 2, len_arr
    if (start_pos == 0) then
       if (.not. is_whitespace_char(array_text(i:i))) start_pos = i
    endif

    ch = array_text(i:i)
    if (in_string) then
       if (escaped) then
          escaped = .false.
       else if (ch == '\') then
          escaped = .true.
       else if (ch == '"') then
          in_string = .false.
       endif
    else
       select case (ch)
       case('"')
          in_string = .true.
       case('[')
          depth = depth + 1
       case('{')
          depth = depth + 1
       case(']')
          if (depth == 0) then
             if (start_pos /= 0) then
                end_pos = i - 1
                current_index = current_index + 1
                if (current_index == index) exit
                start_pos = 0
             endif
             exit
          else
             depth = depth - 1
          endif
       case('}')
          depth = depth - 1
       case(',')
          if (depth == 0) then
             if (start_pos /= 0) then
                end_pos = i - 1
                current_index = current_index + 1
                if (current_index == index) exit
                start_pos = 0
             endif
          endif
       end select
    endif
 enddo

 if (current_index /= index) then
    if (end_pos == 0 .and. start_pos /= 0) then
       end_pos = len_arr - 1
       current_index = current_index + 1
    endif
 endif

 if (current_index /= index .or. start_pos == 0) then
    ierr = json_err_not_found
    allocate(character(len=0) :: value_text)
    kind = json_kind_unknown
    return
 endif

 call extract_trimmed_span(array_text,start_pos,end_pos,value_text)
 kind = detect_value_kind(value_text)

end subroutine json_array_get_element

!-----------------------------------------------------------------
! get the number of items in a JSON array
!-----------------------------------------------------------------
subroutine json_array_length(array_text, nitems, ierr)
 character(len=*), intent(in) :: array_text
 integer, intent(out) :: nitems
 integer, intent(out) :: ierr
 integer :: len_arr, i, depth
 logical :: in_string, escaped
 character :: ch
 logical :: expecting_value

 ierr = json_success
 nitems = 0

 len_arr = len_trim(array_text)
 if (len_arr < 2 .or. array_text(1:1) /= '[') then
    ierr = json_err_type
    return
 endif

 in_string = .false.
 escaped = .false.
 depth = 0
 expecting_value = .true.

 do i = 2, len_arr
    ch = array_text(i:i)
    if (in_string) then
       if (escaped) then
          escaped = .false.
       else if (ch == '\') then
          escaped = .true.
       else if (ch == '"') then
          in_string = .false.
       endif
       cycle
    endif

    select case (ch)
    case('"')
       in_string = .true.
       if (depth == 0 .and. expecting_value) then
          nitems = nitems + 1
          expecting_value = .false.
       endif
    case('[','{')
       if (depth == 0 .and. expecting_value) then
          nitems = nitems + 1
          expecting_value = .false.
       endif
       depth = depth + 1
    case(']', '}')
       depth = depth - 1
    case(',')
       if (depth == 0) expecting_value = .true.
    case default
       if (depth == 0 .and. expecting_value .and. .not. is_whitespace_char(ch)) then
          nitems = nitems + 1
          expecting_value = .false.
       endif
    end select
 enddo

end subroutine json_array_length

!-----------------------------------------------------------------
! get the raw value from a JSON key : value pair
!-----------------------------------------------------------------
subroutine json_get_raw_value(json_text, key, raw_value, kind, ierr)
 character(len=*), intent(in)  :: json_text
 character(len=*), intent(in)  :: key
 character(:), allocatable, intent(out) :: raw_value
 integer, intent(out) :: kind
 integer, intent(out) :: ierr
 integer :: key_pos, colon_pos, value_start, value_end
 character :: first_char
 character(:), allocatable :: string_value

 ierr = json_success
 kind = json_kind_unknown
 key_pos = find_key_position(json_text,trim(key))
 if (key_pos <= 0) then
    ierr = json_err_not_found
    allocate(character(len=0) :: raw_value)
    return
 endif

 colon_pos = find_colon_after(json_text,key_pos,len_trim(key))
 if (colon_pos <= 0) then
    ierr = json_err_parse
    allocate(character(len=0) :: raw_value)
    return
 endif

 value_start = colon_pos + 1
 call skip_whitespace(json_text,value_start)
 if (value_start > len(json_text)) then
    ierr = json_err_parse
    allocate(character(len=0) :: raw_value)
    return
 endif

 first_char = json_text(value_start:value_start)
 select case (first_char)
 case('"')
    value_end = value_start
    call parse_string(json_text,value_start,string_value,ierr)
    if (ierr /= json_success) then
       allocate(character(len=0) :: raw_value)
       return
    endif
    call move_alloc(string_value,raw_value)
    kind = json_kind_string
    return
 case('{')
    value_end = find_matching_delimiter(json_text,value_start,'{','}')
    if (value_end <= 0) then
       ierr = json_err_parse
       allocate(character(len=0) :: raw_value)
       return
    endif
    call extract_trimmed_span(json_text,value_start,value_end,raw_value)
    kind = json_kind_object
    return
 case('[')
    value_end = find_matching_delimiter(json_text,value_start,'[',']')
    if (value_end <= 0) then
        ierr = json_err_parse
        allocate(character(len=0) :: raw_value)
        return
    endif
    call extract_trimmed_span(json_text,value_start,value_end,raw_value)
    kind = json_kind_array
    return
 case('t','f')
    call parse_literal(json_text,value_start,raw_value,value_end,ierr)
    if (ierr /= json_success) return
    kind = json_kind_boolean
    return
 case('n')
    call parse_literal(json_text,value_start,raw_value,value_end,ierr)
    if (ierr /= json_success) return
    kind = json_kind_null
    return
 case default
    call parse_number(json_text,value_start,raw_value,value_end,ierr)
    if (ierr /= json_success) return
    kind = json_kind_number
    return
 end select

end subroutine json_get_raw_value

!-----------------------------------------------------------------
! parse a JSON path into a list of tokens
!-----------------------------------------------------------------
subroutine parse_path(path, tokens, ntokens, ierr)
 character(len=*), intent(in) :: path
 type(json_path_token), allocatable, intent(out) :: tokens(:)
 integer, intent(out) :: ntokens
 integer, intent(out) :: ierr
 integer :: i, n, start
 type(json_path_token) :: token
 character(len=:), allocatable :: segment
 integer :: last
 integer :: ios, idx_value

 ierr = json_success
 ntokens = 0
 if (len_trim(path) == 0) then
    ierr = json_err_parse
    allocate(tokens(0))
    return
 endif

 i = 1
 n = len_trim(path)
 allocate(tokens(0))

 do while (i <= n)
    call skip_whitespace(path,i)
    if (i > n) exit
    start = i
    do while (i <= n)
       if (path(i:i) == '.' .or. path(i:i) == '[') exit
       i = i + 1
    enddo
    last = i - 1
    if (last >= start) then
       segment = path(start:last)
       token%key = trim(segment)
       token%has_index = .false.
       token%index = -1
       call append_token(tokens,ntokens,token)
    else
       token%key = ''
       token%has_index = .false.
       token%index = -1
       call append_token(tokens,ntokens,token)
    endif

    do while (i <= n)
       if (path(i:i) /= '[') exit
       i = i + 1
       start = i
       do while (i <= n)
          if (path(i:i) == ']') exit
          i = i + 1
       enddo
       if (i > n) then
          ierr = json_err_parse
          deallocate(tokens)
          allocate(tokens(0))
          ntokens = 0
          return
       endif
       read(path(start:i-1),*,iostat=ios) idx_value
       if (ios /= 0) then
          ierr = json_err_parse
          deallocate(tokens)
          allocate(tokens(0))
          ntokens = 0
          return
       endif
       if (.not. tokens(ntokens)%has_index) then
          tokens(ntokens)%has_index = .true.
          tokens(ntokens)%index = idx_value
       else
          token%key = ''
          token%has_index = .true.
          token%index = idx_value
          call append_token(tokens,ntokens,token)
       endif
       i = i + 1
    enddo

   if (i <= n) then
      if (path(i:i) == '.') i = i + 1
   endif
 enddo

end subroutine parse_path

!-----------------------------------------------------------------
! append a token to a list of tokens
!-----------------------------------------------------------------
subroutine append_token(tokens, ntokens, token)
 type(json_path_token), allocatable, intent(inout) :: tokens(:)
 integer, intent(inout) :: ntokens
 type(json_path_token), intent(in) :: token
 type(json_path_token), allocatable :: tmp(:)

 allocate(tmp(ntokens+1))
 if (ntokens > 0) tmp(1:ntokens) = tokens
 tmp(ntokens+1) = token
 ntokens = ntokens + 1
 call move_alloc(tmp,tokens)

end subroutine append_token

!-----------------------------------------------------------------
! parse a literal value from a JSON text
!-----------------------------------------------------------------
subroutine parse_literal(text, pos, value, end_pos, ierr)
 character(len=*), intent(in) :: text
 integer, intent(inout) :: pos
 character(:), allocatable, intent(out) :: value
 integer, intent(out) :: end_pos
 integer, intent(out) :: ierr
 integer :: len_text, start

 ierr = json_success
 len_text = len(text)
 start = pos

 do while (pos <= len_text)
    if (.not. is_identifier_char(text(pos:pos))) exit
    pos = pos + 1
 enddo
 end_pos = pos - 1
 call extract_trimmed_span(text,start,end_pos,value)

 select case (trim(value))
 case('true','false','null')
    ! valid literal
 case default
    ierr = json_err_parse
 end select

end subroutine parse_literal

!-----------------------------------------------------------------
! parse a number value from a JSON text
!-----------------------------------------------------------------
subroutine parse_number(text, pos, value, end_pos, ierr)
 character(len=*), intent(in) :: text
 integer, intent(inout) :: pos
 character(:), allocatable, intent(out) :: value
 integer, intent(out) :: end_pos
 integer, intent(out) :: ierr
 integer :: len_text, start
 character :: ch

 ierr = json_success
 len_text = len(text)
 start = pos

 do while (pos <= len_text)
    ch = text(pos:pos)
    if (.not. is_number_char(ch)) exit
    pos = pos + 1
 enddo
 end_pos = pos - 1
 call extract_trimmed_span(text,start,end_pos,value)

end subroutine parse_number

!-----------------------------------------------------------------
! parse a string value from a JSON text
!-----------------------------------------------------------------
subroutine parse_string(text, pos, value, ierr)
 character(len=*), intent(in) :: text
 integer, intent(inout) :: pos
 character(:), allocatable, intent(out) :: value
 integer, intent(out) :: ierr
 integer :: len_text, i, count
 logical :: escaped
 character :: ch
 character(:), allocatable :: buffer

 ierr = json_success
 len_text = len(text)
 if (pos > len_text .or. text(pos:pos) /= '"') then
    ierr = json_err_parse
    allocate(character(len=0) :: value)
    return
 endif

 pos = pos + 1
 allocate(character(len=len_text) :: buffer)
 count = 0
 escaped = .false.

i = pos
do while (i <= len_text)
   ch = text(i:i)
   if (escaped) then
      count = count + 1
      select case (ch)
      case('"','\','/')
         buffer(count:count) = ch
      case('b')
         buffer(count:count) = char(8)
      case('f')
         buffer(count:count) = char(12)
      case('n')
         buffer(count:count) = char(10)
      case('r')
         buffer(count:count) = char(13)
      case('t')
         buffer(count:count) = char(9)
      case('u')
         ! basic handling: skip unicode escape and insert placeholder
         buffer(count:count) = '?'
         if (i+4 <= len_text) i = i + 4
      case default
         buffer(count:count) = ch
      end select
      escaped = .false.
   else
      select case (ch)
      case('\')
         escaped = .true.
      case('"')
         pos = i + 1
         exit
      case default
         count = count + 1
         buffer(count:count) = ch
      end select
   endif
   i = i + 1
enddo

 if (escaped) then
    ierr = json_err_parse
    deallocate(buffer)
    allocate(character(len=0) :: value)
    return
 endif

 allocate(character(len=count) :: value)
 if (count > 0) value = buffer(1:count)
 deallocate(buffer)

end subroutine parse_string

!-----------------------------------------------------------------
! extract a trimmed span from a JSON text
!-----------------------------------------------------------------
subroutine extract_trimmed_span(text, start_pos, end_pos, span)
 character(len=*), intent(in) :: text
 integer, intent(in) :: start_pos, end_pos
 character(:), allocatable, intent(out) :: span
 integer :: s, e

 s = start_pos
 e = end_pos
 do while (s <= e .and. is_whitespace_char(text(s:s)))
    s = s + 1
 enddo
 do while (e >= s .and. is_whitespace_char(text(e:e)))
    e = e - 1
 enddo
 if (e < s) then
    allocate(character(len=0) :: span)
 else
    allocate(character(len=e-s+1) :: span)
    span = text(s:e)
 endif

end subroutine extract_trimmed_span

!-----------------------------------------------------------------
! find the position of the matching delimiter in a JSON text
!-----------------------------------------------------------------
integer function find_matching_delimiter(text, start_pos, open_char, close_char) result(pos)
 character(len=*), intent(in) :: text
 integer, intent(in) :: start_pos
 character, intent(in) :: open_char, close_char
 integer :: depth, i, len_text
 logical :: in_string, escaped
 character :: ch

 len_text = len(text)
 depth = 0
 in_string = .false.
 escaped = .false.
 pos = 0

 do i = start_pos, len_text
    ch = text(i:i)
    if (in_string) then
       if (escaped) then
          escaped = .false.
       else if (ch == '\') then
          escaped = .true.
       else if (ch == '"') then
          in_string = .false.
       endif
    else
       if (ch == '"') then
          in_string = .true.
       else if (ch == open_char) then
          depth = depth + 1
       else if (ch == close_char) then
          depth = depth - 1
          if (depth == 0) then
             pos = i
             return
          endif
       endif
    endif
 enddo

end function find_matching_delimiter

!-----------------------------------------------------------------
! find the position of a key in a JSON text
!-----------------------------------------------------------------
integer function find_key_position(text, key) result(pos)
 character(len=*), intent(in) :: text
 character(len=*), intent(in) :: key
 integer :: len_text, len_key, i
 logical :: in_string, escaped
 character :: ch

 len_text = len(text)
 len_key = len_trim(key)
 pos = 0
 in_string = .false.
 escaped = .false.

 do i = 1, len_text - len_key - 1
    ch = text(i:i)
    if (in_string) then
       if (escaped) then
          escaped = .false.
       else if (ch == '\') then
          escaped = .true.
       else if (ch == '"') then
          in_string = .false.
       endif
    else
       if (ch == '"') then
          if (i + len_key < len_text) then
             if (text(i+1:i+len_key) == key .and. text(i+len_key+1:i+len_key+1) == '"') then
                pos = i
                return
             endif
          endif
          in_string = .true.
       endif
    endif
 enddo
end function find_key_position

!-----------------------------------------------------------------
! find the position of the colon after a key in a JSON text
!-----------------------------------------------------------------
integer function find_colon_after(text, key_pos, key_len) result(pos)
 character(len=*), intent(in) :: text
 integer, intent(in) :: key_pos, key_len
 integer :: i, len_text

 len_text = len(text)
 pos = 0
 i = key_pos + key_len + 2

 do while (i <= len_text)
    if (.not. is_whitespace_char(text(i:i))) then
       if (text(i:i) == ':') then
          pos = i
          return
       else
          return
       endif
    endif
    i = i + 1
 enddo

end function find_colon_after

!-----------------------------------------------------------------
! skip whitespace in a JSON text
!-----------------------------------------------------------------
subroutine skip_whitespace(text, pos)
 character(len=*), intent(in) :: text
 integer, intent(inout) :: pos
 integer :: len_text

 len_text = len(text)
 do while (pos <= len_text)
    if (.not. is_whitespace_char(text(pos:pos))) exit
    pos = pos + 1
 enddo

end subroutine skip_whitespace

!-----------------------------------------------------------------
! check if a character is a whitespace character
!-----------------------------------------------------------------
logical function is_whitespace_char(ch) result(res)
 character, intent(in) :: ch

 res = (ch == ' ' .or. ch == char(9) .or. ch == char(10) .or. ch == char(13))

end function is_whitespace_char

!-----------------------------------------------------------------
! check if a character is an identifier character
!-----------------------------------------------------------------
logical function is_identifier_char(ch) result(res)
 character, intent(in) :: ch

 res = ((ch >= 'a' .and. ch <= 'z') .or. (ch >= 'A' .and. ch <= 'Z') .or. (ch >= '0' .and. ch <= '9') .or. ch == '_')

end function is_identifier_char

!-----------------------------------------------------------------
! check if a character is a number character
!-----------------------------------------------------------------
logical function is_number_char(ch) result(res)
 character, intent(in) :: ch

 res = (ch == '+' .or. ch == '-' .or. ch == '.' .or. ch == 'e' .or. ch == 'E' .or. (ch >= '0' .and. ch <= '9'))

end function is_number_char

!-----------------------------------------------------------------
! detect the kind of a value from a JSON text
!-----------------------------------------------------------------
integer function detect_value_kind(value_text) result(kind)
 character(len=*), intent(in) :: value_text
 character :: first_char
 integer :: i, n

 n = len(value_text)
 first_char = ' '
 do i = 1, n
    if (.not. is_whitespace_char(value_text(i:i))) then
       first_char = value_text(i:i)
       exit
    endif
 enddo

 if (first_char == ' ') then
    kind = json_kind_unknown
    return
 endif
 select case (first_char)
 case('{')
    kind = json_kind_object
 case('[')
    kind = json_kind_array
 case('"')
    kind = json_kind_string
 case('t','f')
    if (trim(value_text) == 'true' .or. trim(value_text) == 'false') then
       kind = json_kind_boolean
    else
       kind = json_kind_unknown
    endif
 case('n')
    if (trim(value_text) == 'null') then
       kind = json_kind_null
    else
       kind = json_kind_unknown
    endif
 case default
    kind = json_kind_number
 end select

end function detect_value_kind

!-----------------------------------------------------------------
! parse an array of 64 bit integers from a JSON text
!-----------------------------------------------------------------
subroutine parse_int64_array(array_text,values,ierr)
 character(len=*), intent(in) :: array_text
 integer(kind=int64), allocatable, intent(out) :: values(:)
 integer, intent(out) :: ierr
 integer :: n,i,value_kind
 character(len=:), allocatable :: elem

 ierr = 0
 call json_array_length(array_text,n,ierr)
 if (ierr /= json_success) then
    ierr = 1
    allocate(values(0))
    return
 endif
 allocate(values(n))
 do i = 1, n
    call json_array_get_element(array_text,i-1,elem,value_kind,ierr)
    if (ierr == json_success) then
       read(elem,*,iostat=ierr) values(i)
    else
       ierr = 1
    endif
    if (allocated(elem)) deallocate(elem)
    if (ierr /= 0) then
       values = 0_int64
       return
    endif
 enddo

end subroutine parse_int64_array

!-----------------------------------------------------------------
! parse an array of integers from a JSON text
!-----------------------------------------------------------------
subroutine parse_int_array(array_text,values,ierr)
 character(len=*), intent(in) :: array_text
 integer, allocatable, intent(out) :: values(:)
 integer, intent(out) :: ierr
 integer :: n,i,value_kind
 character(len=:), allocatable :: elem

 ierr = 0
 call json_array_length(array_text,n,ierr)
 if (ierr /= json_success) then
    ierr = 1
    allocate(values(0))
    return
 endif
 allocate(values(n))
 do i = 1, n
    call json_array_get_element(array_text,i-1,elem,value_kind,ierr)
    if (ierr == json_success) then
       read(elem,*,iostat=ierr) values(i)
    else
       ierr = 1
    endif
    if (allocated(elem)) deallocate(elem)
    if (ierr /= 0) then
       values = 0
       return
    endif
 enddo

end subroutine parse_int_array

end module json_utils

