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
!  Copyright (C) 2005-2020 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!
!  Shamrock native binary reader for SPLASH
!
!-----------------------------------------------------------------
module readdata_shamrock
 use params,          only:maxplot
 use iso_fortran_env, only:int64,int8
 use json_utils,      only:json_read,json_get_value_by_path, &
                           json_array_length,json_array_get_element,json_success, &
                           json_kind_number,json_kind_string,json_kind_object, &
                           json_kind_array,parse_int64_array,parse_int_array
 implicit none
 private

 type :: shamrock_field
    character(len=64) :: name = ''
    character(len=16) :: ftype = ''
    integer           :: ncomp = 0
    integer           :: col_start = 0
    logical           :: is_vector = .false.
 end type shamrock_field

 character(len=32) :: tagarr(maxplot)

 public :: read_data_shamrock, set_labels_shamrock

contains

!-----------------------------------------------------------------
! read the data from a Shamrock file
!-----------------------------------------------------------------
subroutine read_data_shamrock(rootname,istepstart,ipos,nstepsread)
 use particle_data,  only:dat,npartoftype,masstype,maxcol,maxpart,headervals
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread,iverbose,debugmode
 use mem_allocation, only:alloc
 use labels,         only:headertags,irho,ih
 use asciiutils,     only:numfromfile,match_tag
 character(len=*), intent(in)  :: rootname
 integer,          intent(in)  :: istepstart,ipos
 integer,          intent(out) :: nstepsread

 type(shamrock_field), allocatable :: fields(:)
 logical :: iexist,reallocate
 integer :: nfields,ncolstep,npart,nsteps_to_read,idx_xyz,nheader,npatches,ierr,i,iu
 integer, allocatable :: pids(:)
 integer(kind=int64) :: base_pos
 integer(kind=int64), allocatable :: offsets(:), bytecounts(:), field_counts(:)
 character(len=:),    allocatable :: user_header,sched_header,file_header
 character(len=len(rootname)+16) :: datfile
 real :: kernel_hfact, particle_mass

 nstepsread = 0
 nsteps_to_read = 1
 npart = 0
 ncolstep = 4

 if (len_trim(rootname) == 0) then
    print*,' **** no data read **** '
    return
 endif

 datfile = trim(rootname)
 inquire(file=datfile,exist=iexist)
 if (.not. iexist) then
    datfile = trim(rootname)//'.sham'
    inquire(file=datfile,exist=iexist)
 endif
 if (.not. iexist) then
    print "(a)",' *** error: '//trim(rootname)//': file not found ***'
    return
 endif

 if (iverbose==1 .and. ipos==1) print "(1x,a)",'reading Shamrock native format'
 open(newunit=iu,file=trim(datfile),status='old',access='stream',&
      form='unformatted',action='read',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR opening '//trim(datfile)
    return
 endif

 call read_int64_and_string(iu,user_header,ierr)
 if (ierr /= 0) print "(a)",' ERROR reading user header'
 call read_int64_and_string(iu,sched_header,ierr)
 if (ierr /= 0) print "(a)",' ERROR reading scheduler header'
 call read_int64_and_string(iu,file_header,ierr)
 if (ierr /= 0) print "(a)",' ERROR reading file header'

 inquire(unit=iu,pos=base_pos,iostat=ierr)
 if (ierr /= 0) return

 call read_shamrock_header(user_header,sched_header,file_header, &
      fields,offsets,bytecounts,pids,kernel_hfact,particle_mass,ierr)
 if (ierr /= 0) return

 nfields = size(fields)
 tagarr(:) = ''
 call initialise_columns(fields,nfields,ncolstep,ierr)
 if (ierr /= 0) return

 npatches = size(pids)
 allocate(field_counts(nfields))

 idx_xyz = match_tag(fields(:)%name,'xyz')
 if (idx_xyz <= 0) then
    print*,' ERROR: xyz field missing in Shamrock layout'
    return
 endif

 npart = 0
 do i = 1,npatches
    call read_patch_counts(iu,base_pos+offsets(i),nfields,field_counts,ierr)
    if (ierr /= 0) exit
    npart = npart + int(field_counts(idx_xyz))
 enddo

 if (npart <= 0) then
    print*,' ERROR: no particles found in Shamrock dump'
    return
 endif

 ncolumns = ncolstep
 reallocate = (npart > maxpart)
 if (reallocate .or. .not.(allocated(dat))) then
    call alloc(npart,nsteps_to_read,max(ncolumns+ncalc,maxcol),mixedtypes=.false.)
 endif

 dat(1:npart,1:ncolumns,istepstart) = 0.0
 ipartialread = .false.
 ndim = 3
 ndimV = 3

 call read_shamrock_data(iu,base_pos,npatches,fields,offsets, &
      istepstart,npart,field_counts,idx_xyz,ierr)

 if (ierr /= 0) print "(a)",' ERROR reading Shamrock data'
 close(iu)

 call set_labels_shamrock

 ! extract all possible values from the user header
 nheader = 4
 headertags(1:4) = (/'npart','pmass','hfact','nfile'/)
 headervals(1,istepstart) = real(npart)
 headervals(2,istepstart) = real(particle_mass)
 headervals(3,istepstart) = real(kernel_hfact)
 headervals(4,istepstart) = real(numfromfile(datfile))
 if (debugmode) print*,user_header
 call get_header_pairs(user_header,headertags,headervals(:,istepstart),nheader)

 masstype(1,istepstart) = real(particle_mass)
 npartoftype(1,istepstart) = npart

 ! compute the density from the mass and the smoothing length
 if (ih > 0 .and. irho > 0) then
    do i=1,npart
       dat(i,irho,istepstart) = particle_mass*(kernel_hfact/abs(dat(i,ih,istepstart)))**3
    enddo
 endif
 nstepsread = 1

 if (allocated(field_counts)) deallocate(field_counts)
 if (allocated(fields)) deallocate(fields)
 if (allocated(offsets)) deallocate(offsets)
 if (allocated(bytecounts)) deallocate(bytecounts)
 if (allocated(pids)) deallocate(pids)

end subroutine read_data_shamrock

!-----------------------------------------------------------------
! set the labels for the data
!-----------------------------------------------------------------
subroutine set_labels_shamrock
 use labels,        only:label,labeltype,ix,iamvec,labelvec,irho,ih,ipmass,  &
                          set_vector_labels,label_synonym,ivx
 use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings,ncolumns
 use geometry,      only:labelcoord
 integer :: i

 ix(1) = 1
 ix(2) = 2
 if (ndim >= 3) ix(3) = 3

 do i = 1, ncolumns
    label(i) = trim(tagarr(i))
    select case (trim(label(i)))
    case('rho')
       irho = i
    case('h','hpart')
       ih = i
    case('mass')
       ipmass = i
    case('vx')
       ivx = i
    end select
 enddo

 call set_vector_labels(ncolumns,ndimV,iamvec,labelvec,label,labelcoord(:,1))

 if (ix(1) > 0) label(ix(1:ndim)) = labelcoord(1:ndim,1)
 if (irho > 0)  label(irho) = 'density'
 if (ipmass > 0) label(ipmass) = 'particle mass'
 if (ih > 0)   label(ih) = 'h'

 ntypes = 1
 labeltype(1) = 'gas'
 UseTypeInRenderings(:) = .true.

end subroutine set_labels_shamrock

!-----------------------------------------------------------------
! initialise the columns
!-----------------------------------------------------------------
subroutine initialise_columns(fields,nfields,ncolstep,ierr)
 use labels, only:iamvec
 type(shamrock_field), intent(inout) :: fields(:)
 integer, intent(in) :: nfields
 integer, intent(inout) :: ncolstep
 integer, intent(out) :: ierr
 integer :: i

 ierr = 0
 tagarr(:) = ''
 ncolstep = 4
 tagarr(1) = 'x'
 tagarr(2) = 'y'
 tagarr(3) = 'z'
 tagarr(4) = 'rho'
 iamvec(:) = 0

 do i = 1, nfields
    select case (trim(fields(i)%ftype))
    case('f64_3')
       fields(i)%ncomp = 3
       if (trim(fields(i)%name) == 'xyz') then
          fields(i)%col_start = 1
          fields(i)%is_vector = .false.
       else
          call register_vector_field(fields(i),ncolstep,ierr)
          if (ierr /= 0) return
       endif
    case('f64')
       fields(i)%ncomp = 1
       call register_scalar_field(fields(i),ncolstep,ierr)
       if (ierr /= 0) return
    case default
       print*,' WARNING: unsupported field type ',trim(fields(i)%ftype)
    end select
 enddo

end subroutine initialise_columns

!-----------------------------------------------------------------
! register a scalar field
!-----------------------------------------------------------------
subroutine register_scalar_field(field,ncolstep,ierr)
 type(shamrock_field), intent(inout) :: field
 integer, intent(inout) :: ncolstep
 integer, intent(out) :: ierr
 character(len=32) :: label

 ierr = 0
 label = field%name(1:min(len_trim(field%name),len(label)))
 field%col_start = find_or_add_label(label,ncolstep,ierr)
 field%is_vector = .false.

end subroutine register_scalar_field

!-----------------------------------------------------------------
! register a vector field
!-----------------------------------------------------------------
subroutine register_vector_field(field,ncolstep,ierr)
 use labels, only:iamvec
 type(shamrock_field), intent(inout) :: field
 integer, intent(inout) :: ncolstep
 integer, intent(out) :: ierr
 character(len=32) :: vec_labels(3),base_label
 integer :: base, name_len

 ierr = 0
 vec_labels(:) = ''
 base_label = ''
 select case (trim(field%name))
 case('vxyz')
    base_label = 'v'
 case('axyz')
    base_label = 'a'
 case('axyz_ext')
    base_label = 'aext'
 case default
    name_len = min(len_trim(field%name), len(base_label))
    if (name_len > 0) then
       base_label = field%name(1:name_len)
    else
       base_label = 'vec'
    endif
 end select

 vec_labels = (/trim(base_label),trim(base_label),trim(base_label)/)

 base = find_or_add_vector(vec_labels,ncolstep,ierr)
 if (ierr /= 0) return
 field%col_start = base
 field%is_vector = .true.
 iamvec(base:base+2) = base

end subroutine register_vector_field

!-----------------------------------------------------------------
! find the index of a label
!-----------------------------------------------------------------
integer function find_or_add_label(label,ncolstep,ierr) result(idx)
 character(len=*), intent(in) :: label
 integer, intent(inout) :: ncolstep
 integer, intent(out) :: ierr
 integer :: i

 ierr = 0
 do i = 1, ncolstep
    if (trim(tagarr(i)) == trim(label)) then
       idx = i
       return
    endif
 enddo
 ncolstep = ncolstep + 1
 if (ncolstep > maxplot) then
    ierr = 1
    idx = 0
    return
 endif
 tagarr(ncolstep) = trim(label)
 idx = ncolstep

end function find_or_add_label

!-----------------------------------------------------------------
! find the index of a vector
!-----------------------------------------------------------------
integer function find_or_add_vector(labeli,ncolstep,ierr) result(idx)
 character(len=*), intent(in) :: labeli(3)
 integer, intent(inout) :: ncolstep
 integer, intent(out) :: ierr
 integer :: i, base

 ierr = 0
 do i = 1, ncolstep-2
    if (trim(tagarr(i)) == trim(labeli(1)) .and. trim(tagarr(i+1)) == trim(labeli(2)) &
        .and. trim(tagarr(i+2)) == trim(labeli(3))) then
       idx = i
       return
    endif
 enddo

 base = ncolstep + 1
 if (base+2 > maxplot) then
    ierr = 1
    idx = 0
    return
 endif
 tagarr(base) = trim(labeli(1))
 tagarr(base+1) = trim(labeli(2))
 tagarr(base+2) = trim(labeli(3))
 ncolstep = base + 2
 idx = base

end function find_or_add_vector

!-----------------------------------------------------------------
! read the header of the Shamrock file
!-----------------------------------------------------------------
subroutine read_shamrock_header(user_header,sched_header,file_header,fields,offsets, &
                                bytecounts,pids,kernel_hfact,particle_mass,ierr)
 character(len=*),                  intent(in)  :: user_header,sched_header,file_header
 type(shamrock_field), allocatable, intent(out) :: fields(:)
 integer(kind=int64),  allocatable, intent(out) :: offsets(:), bytecounts(:)
 integer, allocatable,              intent(out) :: pids(:)
 real,                              intent(out) :: kernel_hfact, particle_mass
 integer,                           intent(out) :: ierr
 character(len=:), allocatable :: layout_text
 character(len=:), allocatable :: elem_text
 character(len=:), allocatable :: tmp
 character(len=:), allocatable :: arr_text
 integer :: nfields,kind,i,istat,arr_len

 ierr = 0
 kernel_hfact = 1.2
 particle_mass = 0.0

 call extract_real(user_header,'solver_config.gpart_mass',particle_mass)
 call extract_kernel_hfact(user_header,kernel_hfact)

 istat = json_success
 call json_read(sched_header,'patchdata_layout',layout_text,istat)
 if (istat /= json_success) then
    ierr = 1
    allocate(fields(0), offsets(0), bytecounts(0), pids(0))
    return
 endif

 call json_array_length(layout_text,nfields,istat)
 if (istat /= json_success .or. nfields <= 0) then
    ierr = 1
    if (allocated(layout_text)) deallocate(layout_text)
    allocate(fields(0), offsets(0), bytecounts(0), pids(0))
    return
 endif

 allocate(fields(nfields))
 do i = 1, nfields
    call json_array_get_element(layout_text,i-1,elem_text,kind,istat)
    if (istat /= json_success .or. kind /= json_kind_object) then
       ierr = 1
       exit
    endif
    call json_read(elem_text,'field_name',tmp,istat)
    if (istat /= json_success) then
       ierr = 1
       exit
    endif
    fields(i)%name = trim(tmp)
    if (allocated(tmp)) deallocate(tmp)
    call json_read(elem_text,'type',tmp,istat)
    if (istat /= json_success) then
       ierr = 1
       exit
    endif
    fields(i)%ftype = trim(tmp)
    if (allocated(tmp)) deallocate(tmp)
    if (allocated(elem_text)) deallocate(elem_text)
 enddo
 if (allocated(tmp)) deallocate(tmp)
 if (allocated(elem_text)) deallocate(elem_text)
 if (allocated(layout_text)) deallocate(layout_text)
 if (ierr /= 0) then
    if (allocated(fields)) deallocate(fields)
    allocate(fields(0), offsets(0), bytecounts(0), pids(0))
    return
 endif

 call json_read(file_header,'offsets',arr_text,istat)
 if (istat /= json_success) then
    ierr = 1
    allocate(offsets(0), bytecounts(0), pids(0))
    return
 endif
 call parse_int64_array(arr_text,offsets,ierr)
 if (allocated(arr_text)) deallocate(arr_text)
 if (ierr /= 0) return

 call json_read(file_header,'bytecounts',arr_text,istat)
 if (istat /= json_success) then
    ierr = 1
    allocate(bytecounts(0), pids(0))
    return
 endif
 call parse_int64_array(arr_text,bytecounts,ierr)
 if (allocated(arr_text)) deallocate(arr_text)
 if (ierr /= 0) return

 call json_read(file_header,'pids',arr_text,istat)
 if (istat /= json_success) then
    ierr = 1
    allocate(pids(0))
    return
 endif
 call parse_int_array(arr_text,pids,ierr)
 if (allocated(arr_text)) deallocate(arr_text)
 if (ierr /= 0) return

 arr_len = size(pids)
 if (size(offsets) /= arr_len .or. size(bytecounts) /= arr_len) then
    ierr = 1
 endif

end subroutine read_shamrock_header

!-----------------------------------------------------------------
! get the header pairs from the user header by parsing the JSON
! string for key:value pairs
!-----------------------------------------------------------------
subroutine get_header_pairs(user_header,headertags,headervals,nheader)
 use settings_data, only:debugmode
 character(len=*), intent(in) :: user_header
 character(len=*), intent(out) :: headertags(:)
 real, intent(out) :: headervals(:)
 integer, intent(inout) :: nheader
 character(len=len(user_header)) :: string
 character(len=len(headertags(1))) :: key_string
 character(len=32) :: val_string
 integer :: icolon,ierr,i1,i2

 string = user_header
 icolon = index(user_header,':')
 do while (icolon > 0)
    ! key is the thing between the previous two double quotes
    i2 = index(string(1:icolon-1),'"',back=.true.)
    i1 = index(string(1:i2-1),'"',back=.true.)
    key_string = string(i1+1:min(i2-1,i1+len(key_string)))

    ! val is the thing after the colon
    val_string = string(icolon+1:icolon+len(val_string))

    ! try to the value as a real, if this fails, skip it
    read(val_string,*,iostat=ierr) headervals(nheader)
    if (ierr == 0) then
       nheader = nheader + 1
       headertags(nheader) = trim(key_string)
       if (debugmode) print*,trim(key_string)//':',headervals(nheader)
    endif
    string = string(icolon+1:) ! move to the next colon
    icolon = index(string,':')
 enddo

end subroutine get_header_pairs

!-----------------------------------------------------------------
! read 64 bit integer and the string following it from file
!-----------------------------------------------------------------
subroutine read_int64_and_string(iu,str,ierr)
 integer, intent(in) :: iu
 character(len=:), allocatable, intent(out) :: str
 integer, intent(out) :: ierr
 integer(kind=int64) :: len_str

 ierr = 0
 read(iu,iostat=ierr) len_str
 allocate(character(len=int(len_str)) :: str)
 read(iu,iostat=ierr) str

end subroutine read_int64_and_string

!-----------------------------------------------------------------
! read the patch counts from the file
!-----------------------------------------------------------------
subroutine read_patch_counts(iu,pos_start,nfields,counts,ierr)
 integer, intent(in) :: iu
 integer(kind=int64), intent(in) :: pos_start
 integer, intent(in) :: nfields
 integer(kind=int64), intent(out) :: counts(:)
 integer, intent(out) :: ierr
 integer(kind=int64) :: prehead
 integer :: j

 ierr = 0
 read(iu,pos=pos_start,iostat=ierr) prehead
 if (ierr /= 0) return
 do j = 1, nfields
    read(iu,iostat=ierr) counts(j)
    if (ierr /= 0) return
 enddo

end subroutine read_patch_counts

!-----------------------------------------------------------------
! read the data from the Shamrock file into the data arrays
!-----------------------------------------------------------------
subroutine read_shamrock_data(iu,base_pos,npatches,fields,offsets, &
                              istepstart,npart,counts_buffer,idx_xyz,ierr)
 type(shamrock_field), intent(in) :: fields(:)
 integer, intent(in) :: npatches
 integer(kind=int64), intent(in) :: offsets(:),base_pos
 integer(kind=int64), intent(inout) :: counts_buffer(:)
 integer, intent(in) :: iu, istepstart, npart
 integer, intent(in) :: idx_xyz
 integer, intent(out) :: ierr
 integer :: patch,nfields,i1,nobj,ios,npart_got
 integer(kind=int64) :: pos

 ierr = 0
 nfields = size(fields)
 i1 = 0
 npart_got = 0

   do patch = 1,npatches
      pos = base_pos + offsets(patch)
      call read_patch_counts(iu,pos,nfields,counts_buffer,ios)
      if (ios /= 0) then
         ierr = 1
         return
      endif
      nobj = int(counts_buffer(idx_xyz))
      call read_patch_data_block(iu,nfields,fields,counts_buffer, &
           istepstart,i1,ierr)
      if (ierr /= 0) return
      i1 = i1 + nobj
      npart_got = npart_got + nobj
   enddo

 if (npart_got /= npart) then
    print*,' WARNING: particle count mismatch in Shamrock reader'
 endif

end subroutine read_shamrock_data

!-----------------------------------------------------------------
! read the data from a single patch into the data arrays
!-----------------------------------------------------------------
subroutine read_patch_data_block(iu,nfields,fields,counts,istepstart,i1,ierr)
 type(shamrock_field), intent(in)  :: fields(:)
 integer,              intent(in)  :: nfields,istepstart,i1,iu
   integer(kind=int64),  intent(in)  :: counts(:)
 integer,              intent(out) :: ierr
   integer :: idx

 ierr = 0
   ! stream is already positioned at the first field payload by read_patch_counts

 do idx = 1, nfields
    select case (trim(fields(idx)%ftype))
    case('f64_3')
       call read_vector_field(iu,fields(idx),counts(idx),i1,istepstart,ierr)
    case('f64')
       call read_scalar_field(iu,fields(idx),counts(idx),i1,istepstart,ierr)
    case default
       ! skip unsupported field types for now
    end select
    if (ierr /= 0) return
 enddo

end subroutine read_patch_data_block

!-----------------------------------------------------------------
! read the data from a vector field into the data arrays
!-----------------------------------------------------------------
subroutine read_vector_field(iu,field,count,i1,istepstart,ierr)
 use particle_data, only:dat
 type(shamrock_field), intent(in) :: field
 integer(kind=int64),  intent(in) :: count
 integer, intent(in) :: iu, i1, istepstart
 integer, intent(out) :: ierr
 real, allocatable :: flat(:)
 integer :: p,comp,base_col
 integer(kind=int64) :: size_bytes, aligned_bytes, pad_bytes
 integer(kind=int8), allocatable :: pad(:)

 ierr = 0
 if (count <= 0) return

 base_col = field%col_start
 allocate(flat(int(count)*field%ncomp))
 read(iu,iostat=ierr) flat
 if (ierr /= 0) then
    ierr = 1
    deallocate(flat)
    return
 endif
 do p = 1,int(count)
    do comp = 1, field%ncomp
       dat(i1+p,base_col+comp-1,istepstart) = real(flat((p-1)*field%ncomp+comp))
    enddo
 enddo
 deallocate(flat)

 size_bytes = count*field%ncomp*8_int64
 aligned_bytes = align8(size_bytes)
 pad_bytes = aligned_bytes - size_bytes
 if (pad_bytes > 0_int64) then
    allocate(pad(pad_bytes))
    read(iu,iostat=ierr) pad
    deallocate(pad)
    if (ierr /= 0) ierr = 1
 endif

end subroutine read_vector_field

!-----------------------------------------------------------------
! read the data from a scalar field into the data arrays
!-----------------------------------------------------------------
subroutine read_scalar_field(iu,field,count,i1,istepstart,ierr)
 use particle_data, only:dat
 type(shamrock_field), intent(in) :: field
 integer(kind=int64),  intent(in) :: count
 integer, intent(in)  :: iu,i1,istepstart
 integer, intent(out) :: ierr
 real, allocatable :: values(:)
 integer(kind=int64) :: size_bytes, aligned_bytes, pad_bytes
 integer(kind=int8), allocatable :: pad(:)

 ierr = 0
 if (count <= 0) return
 allocate(values(count))
 read(iu,iostat=ierr) values
 if (ierr /= 0) then
    ierr = 1
    deallocate(values)
    return
 endif
 dat(i1+1:i1+count,field%col_start,istepstart) = real(values)
 deallocate(values)

 size_bytes = count*8_int64
 aligned_bytes = align8(size_bytes)
 pad_bytes = aligned_bytes - size_bytes
 if (pad_bytes > 0_int64) then
    allocate(pad(pad_bytes))
    read(iu,iostat=ierr) pad
    deallocate(pad)
    if (ierr /= 0) ierr = 1
 endif

end subroutine read_scalar_field

!-----------------------------------------------------------------
! align the number of bytes to the next multiple of 8
!-----------------------------------------------------------------
integer(kind=int64) function align8(nbytes) result(res)
 integer(kind=int64), intent(in) :: nbytes
 integer(kind=int64) :: rem

 res = nbytes
 rem = modulo(res,8_int64)
 if (rem /= 0_int64) res = res + (8_int64 - rem)

end function align8

!-----------------------------------------------------------------
! extract a real value from a JSON text
!-----------------------------------------------------------------
subroutine extract_real(json_text,path,value)
 use json_utils, only:json_get_value_by_path
 character(len=*), intent(in) :: json_text
 character(len=*), intent(in) :: path
 real, intent(out) :: value
 character(len=:), allocatable :: raw
 integer :: value_kind, ierr, ios

 value = 0.
 call json_get_value_by_path(json_text,path,raw,value_kind,ierr)
 if (ierr == json_success .and. value_kind == json_kind_number) then
    read(raw,*,iostat=ios) value
    if (ios /= 0) value = 0.
 endif
 if (allocated(raw)) deallocate(raw)

end subroutine extract_real

!-----------------------------------------------------------------
! extract the kernel hfact from a JSON text
!-----------------------------------------------------------------
subroutine extract_kernel_hfact(json_text,hfact)
 use json_utils, only:json_get_value_by_path
 character(len=*), intent(in) :: json_text
 real, intent(out) :: hfact
 character(len=:), allocatable :: kernel
 integer :: value_kind, ierr

 hfact = 1.2
 call json_get_value_by_path(json_text,'solver_config.kernel_id',kernel,value_kind,ierr)
 if (ierr == json_success .and. value_kind == json_kind_string) then
    select case (kernel(1:min(2,len_trim(kernel))))
    case('M4','M5')
       hfact = 1.2
    case('M6')
       hfact = 1.0
    case default
       hfact = 1.2
    end select
 endif
 if (allocated(kernel)) deallocate(kernel)

end subroutine extract_kernel_hfact

end module readdata_shamrock

