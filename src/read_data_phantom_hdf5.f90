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
!
! Standalone reader for Phantom / sphNG HDF5 dumps (C HDF5 via iso_c_binding).
!
!-----------------------------------------------------------------
module phantomhdf5read
 use params, only:maxplot
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none

 character(len=lenlabel), dimension(maxplot) :: blocklabel
 integer, dimension(maxplot) :: blocksize
 integer :: phantom_hdf5_j = 1

 interface
  subroutine phantom_hdf5_header(filename,time,npart,ncol,has_header,ndim,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   real(kind=c_double), intent(out) :: time
   integer(kind=c_int), intent(out) :: npart,ncol,has_header,ndim,ierr
  end subroutine phantom_hdf5_header

  subroutine phantom_hdf5_data(filename,npart,ncol,isrequired,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   integer(kind=c_int), intent(in) :: npart,ncol
   integer(kind=c_int), intent(out) :: ierr
   integer(kind=c_int), dimension(ncol), intent(in) :: isrequired
  end subroutine phantom_hdf5_data

  integer(c_int) function phantom_hdf5_is_phantom_file(filename) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
  end function phantom_hdf5_is_phantom_file
 end interface

end module phantomhdf5read

!-----------------------------------------------------------------
!
! read_data and set_labels for Phantom HDF5 format (-f phantom_hdf5)
!
!-----------------------------------------------------------------
module readdata_phantom_hdf5
 use phantomhdf5read
 implicit none

 public :: read_data_phantom_hdf5, set_labels_phantom_hdf5, file_format_is_phantom_hdf5

 private :: column_looks_like_sph_smoothing_length

contains

!-----------------------------------------------------------------
! read Phantom HDF5 dump into dat (standalone format reader)
!-----------------------------------------------------------------
subroutine read_data_phantom_hdf5(rootname,indexstart,ipos,nstepsread)
 use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol
 use params,         only:maxplot
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,lowmemorymode
 use mem_allocation, only:alloc
 use asciiutils,     only:cstring
 use phantomhdf5read, only:phantom_hdf5_header,phantom_hdf5_data
 integer, intent(in)  :: indexstart,ipos
 integer, intent(out) :: nstepsread
 character(len=*), intent(in) :: rootname
 character(len=len(rootname)+10) :: datfile
 integer :: i,j,ncolstep,ilastrequired
 integer(c_int) :: npart_c,ncol_c,has_header_c,ndim_c,ierr_c
 real(c_double) :: time_c
 logical :: iexist
 integer, dimension(maxplot) :: isrequired

 nstepsread = 0
 if (len_trim(rootname) == 0) return

 datfile = trim(rootname)
 inquire(file=datfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(rootname)//': file not found ***'
    return
 endif

 j = indexstart
 phantom_hdf5_j = j
 blocklabel = ' '
 blocksize = 0
 ndim = 3
 ndimV = 3
 print "(1x,a)",'reading Phantom HDF5 format'
 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)

 call phantom_hdf5_header(cstring(datfile),time_c,npart_c,ncol_c,has_header_c,ndim_c,ierr_c)
 if (ierr_c /= 0) then
    print "(a)",' *** ERROR READING PHANTOM HDF5 HEADER ***'
    return
 endif

 ncolstep = ncol_c
 ndim = max(1,ndim_c)

 print "(a,i10,a,es12.5)",' npart = ',npart_c,' time = ',time_c

 if (ncolstep > maxplot) then
    print "(1x,a,i0,a,i0,a)",'ERROR: file has ',ncolstep,' columns but maxplot=',maxplot
    return
 endif

 if (npart_c > maxpart .or. (ncolstep+ncalc) > maxcol .or. .not.allocated(dat)) then
    if (lowmemorymode) then
       ilastrequired = 0
       do i=1,ncolstep+ncalc
          if (required(i)) ilastrequired = i
       enddo
       call alloc(max(int(npart_c),maxpart),j,ilastrequired,mixedtypes=.true.)
    else
       call alloc(max(int(npart_c),maxpart),max(j,1),max(ncolstep+ncalc,maxcol),mixedtypes=.true.)
    endif
 endif

 ncolumns = ncolstep
 nstepsread = 1
 npartoftype(:,j) = 0
 npartoftype(1,j) = npart_c
 time(j) = real(time_c, kind=kind(time(j)))
 gamma(j) = real(5./3., kind=kind(gamma(j)))
 masstype(1,j) = 0.

 isrequired(:) = 0
 where (required(1:ncolstep)) isrequired(1:ncolstep) = 1
 if (.not.all(required(1:ncolstep))) ipartialread = .true.

 call phantom_hdf5_data(cstring(datfile),npart_c,ncol_c,isrequired(1:ncolstep),ierr_c)
 if (ierr_c /= 0) then
    print "(a)",' *** ERROR READING PHANTOM HDF5 DATA ***'
    return
 endif

 call set_labels_phantom_hdf5

 if (nstepsread > 0) then
    print "(a,i10,a)",' >> read ',npartoftype(1,j),' particles'
 endif

end subroutine read_data_phantom_hdf5

!-----------------------------------------------------------------
! C callback: copy one column from HDF5 read buffer into dat
!-----------------------------------------------------------------
subroutine read_phantom_hdf5_data_fromc(icol,npart,temparr,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 use particle_data, only:dat
 integer(kind=c_int), intent(in) :: icol,npart
 real(kind=c_double), intent(in) :: temparr(npart)
 character(kind=c_char), intent(in) :: name(256)
 integer :: i

 if (icol < 1 .or. icol > size(dat(1,:,1))) return
 if (phantom_hdf5_j < 1 .or. phantom_hdf5_j > size(dat(1,1,:))) return
 do i=1,min(npart,size(dat(:,1,1)))
    dat(i,icol,phantom_hdf5_j) = real(temparr(i), kind=kind(dat(i,icol,phantom_hdf5_j)))
 enddo

end subroutine read_phantom_hdf5_data_fromc

!-----------------------------------------------------------------
! C callback: store HDF5 dataset name and rank for column labelling
!-----------------------------------------------------------------
subroutine set_blocklabel_phantom(icol,irank,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int,c_char
 use phantomhdf5read, only:blocklabel,blocksize
 use asciiutils, only:fstring
 integer(kind=c_int), intent(in) :: icol,irank
 character(kind=c_char), intent(in) :: name(256)

 if (icol < 1 .or. icol > size(blocklabel)) return
 if (name(1) == achar(0)) return
 blocklabel(icol) = trim(fstring(name))
 blocksize(icol) = irank

end subroutine set_blocklabel_phantom

!-----------------------------------------------------------------
! return true if column values look like SPH smoothing length (h)
!-----------------------------------------------------------------
logical function column_looks_like_sph_smoothing_length(icol)
 use particle_data, only:dat,npartoftype
 use phantomhdf5read, only:phantom_hdf5_j
 integer, intent(in) :: icol
 integer :: i,n,npos

 column_looks_like_sph_smoothing_length = .false.
 if (phantom_hdf5_j < 1 .or. phantom_hdf5_j > size(dat,3)) return
 if (icol < 1 .or. icol > size(dat,2)) return
 n = npartoftype(1,phantom_hdf5_j)
 if (n <= 0) return
 npos = 0
 do i=1,min(n,size(dat,1))
    if (dat(i,icol,phantom_hdf5_j) > 0.) npos = npos + 1
 enddo
 column_looks_like_sph_smoothing_length = (npos > n/2)

end function column_looks_like_sph_smoothing_length

!-----------------------------------------------------------------
! set column labels and index variables for Phantom HDF5 read
!-----------------------------------------------------------------
subroutine set_labels_phantom_hdf5
 use labels, only:label,labeltype,ix,ivx,irho,ih,ipmass,iamvec,labelvec,make_vector_label
 use params
 use settings_data, only:ndim,ndimV,UseTypeInRenderings,ntypes,ncolumns
 use geometry, only:labelcoord
 use phantomhdf5read, only:blocklabel,blocksize
 use asciiutils, only:lcase
 integer :: icol,n

 ix = 0
 ivx = 0
 irho = 0
 ih = 0
 ipmass = 0

 do icol=1,ncolumns
    label(icol) = trim(blocklabel(icol))
    select case(trim(lcase(blocklabel(icol))))
    case('x')
       ix(1) = icol
    case('y')
       ix(2) = icol
    case('z')
       ix(3) = icol
    case('vx','vxyz')
       ivx = icol
    case('rho','density')
       irho = icol
    case('h')
       if (column_looks_like_sph_smoothing_length(icol)) ih = icol
    case('pmass','mass')
       ipmass = icol
    end select
    n = blocksize(icol)
    if (n >= 3 .and. ivx == 0 .and. index(trim(lcase(blocklabel(icol))),'v') > 0) then
       ivx = icol
    endif
 enddo

 if (ix(1) > 0 .and. all(ix(1:ndim) > 0)) label(ix(1:ndim)) = labelcoord(1:ndim,1)
 if (ivx > 0) call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))

 ntypes = 1
 labeltype(1) = 'gas'
 UseTypeInRenderings(1) = .true.

end subroutine set_labels_phantom_hdf5

!-----------------------------------------------------------------
! return true if filename is a Phantom / sphNG HDF5 dump
!-----------------------------------------------------------------
logical function file_format_is_phantom_hdf5(filename) result(is_phantom_hdf5)
 use asciiutils,     only:cstring
 use phantomhdf5read, only:phantom_hdf5_is_phantom_file
 use, intrinsic :: iso_c_binding, only:c_int
 character(len=*), intent(in) :: filename
 integer(c_int) :: is_phantom

 is_phantom_hdf5 = .false.
 if (index(filename,'.h5') == 0 .and. index(filename,'.hdf5') == 0) return
 is_phantom = phantom_hdf5_is_phantom_file(cstring(filename))
 is_phantom_hdf5 = (is_phantom == 1)

end function file_format_is_phantom_hdf5

end module readdata_phantom_hdf5
