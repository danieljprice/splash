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
!  Copyright (C) 2023- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING VTK FILES
!
!-------------------------------------------------------------------------
module readdata_vtk
 use params, only:maxplot
 implicit none

 public :: read_data_vtk, set_labels_vtk
 character(len=32) :: tagarr(maxplot)

 private

contains

subroutine read_data_vtk(rootname,istepstart,ipos,nstepsread)
 use particle_data,    only:dat,npartoftype,masstype,maxcol,maxpart,headervals
 use settings_data,    only:ndim,ndimV,ncolumns,ncalc,ipartialread,iverbose
 use mem_allocation,   only:alloc
 use labels,           only:ih,irho,ipmass,headertags
 use system_utils,     only:get_command_option
 use asciiutils,       only:get_value,numfromfile
 integer, intent(in)                :: istepstart,ipos
 integer, intent(out)               :: nstepsread
 character(len=*), intent(in)       :: rootname
 character(len=len(rootname)+10)    :: datfile
 integer               :: i,ierr
 integer               :: ncolstep,npart,nsteps_to_read
 logical               :: iexist,reallocate,is_binary
 real, allocatable     :: pmass(:)
 real :: hfac

 nstepsread = 0
 nsteps_to_read = 1

 if (len_trim(rootname) > 0) then
    datfile = trim(rootname)
 else
    print*,' **** no data read **** '
    return
 endif
!
!--check if first data file exists
!
 if (iverbose==1 .and. ipos==1) print "(1x,a)",'reading vtk format'
 inquire(file=datfile,exist=iexist)
 if (.not.iexist) then
    !
    !--append .vtk on the endif not already present
    !
    datfile=trim(rootname)//'.vtk'
    inquire(file=datfile,exist=iexist)
 endif
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(rootname)//': file not found ***'
    return
 endif
!
!--read data from snapshots
!
 i = istepstart
 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
 !
 !--open file and read header information
 !
 call read_vtk_header(trim(datfile),ndim,ndimV,ncolstep,npart,is_binary,ierr)
 if (ierr > 0) then
    print*,'ERROR reading vtk file: ierr=',ierr
    return
 else
    print "(1x,a,i0,a,i0)",'npart = ',npart,' ncolumns = ',ncolstep
 endif
 !
 !--allocate or reallocate memory for main data array
 !
 ncolumns = ncolstep
 nstepsread = 1
 reallocate = (npart > maxpart)
 if (reallocate .or. .not.(allocated(dat))) then
    call alloc(npart,nsteps_to_read,max(ncolumns+ncalc,maxcol),mixedtypes=.false.)
 endif
!
!--now memory has been allocated, read main data arrays
!
 ipartialread = .false.
 call read_vtk_legacy_binary(trim(datfile),npart,ncolstep,ierr,dat(:,:,i))
 call set_labels_vtk

 ! set header variables from useful information from the file
 headertags(1) = 'npart'
 headervals(1,i) = real(npart)
 headertags(2) = 'nfile'
 headervals(2,i) = real(numfromfile(datfile))

!
!--set particle mass and number of particles of type 1
!
 if (npart >= 1) then
    !
    !--reconstruct particle mass from h and rho
    !
    if (ih > 0 .and. irho > 0 .and. ipmass==0) then
       allocate(pmass(npart))
       hfac = 1.0
       pmass = dat(1:npart,irho,i)*(dat(1:npart,ih,i)/hfac)**ndim
       if (any(abs(pmass - pmass(1)) > 1e-6)) print*,' WARNING: unequal particle masses not supported'
       !print*,pmass(1:10)
       masstype(1,i) = pmass(1)
       deallocate(pmass)
    else
       masstype(1,i) = 1./npart
    endif
    npartoftype(1,i) = npart
 endif

end subroutine read_data_vtk

!----------------------------------------------------
! read header information
!----------------------------------------------------
subroutine read_vtk_header(filename,ndim,ndimV,ncolstep,npart,is_binary,ierr)
 character(len=*)     :: filename
 integer, intent(out) :: ndim,ndimV,ncolstep,npart,ierr
 logical, intent(out) :: is_binary
 character(len=256)   :: version,header
 character(len=6)     :: fileformat
 character(len=30)    :: dataset,datatype
 integer :: iu
 
 ierr = 0
 npart = 0
 ndim = 3
 ndimV = 3
 ncolstep = 0
 tagarr(:) = ''
 
 open(newunit=iu,file=filename,status='old',action='read',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print*,' ERROR opening '//trim(filename)
    return
 endif

 ! read version information from first line
 read(iu,"(a)",iostat=ierr) version
 if (ierr /= 0) then
    print*,' ERROR reading vtk version from first line'
 endif

 ! read data description from second line
 read(iu,"(a)",iostat=ierr) header
 if (ierr /= 0) then
    print*,' ERROR reading vtk header from second line'
 else
    print "(a)",trim(version)//': '//trim(header)
 endif

 ! read file format from third line
 is_binary = .false.
 read(iu,"(a)",iostat=ierr) fileformat
 if (ierr /= 0) then
    print*,' ERROR reading vtk file format from third line'
 else
    if (trim(fileformat)=='BINARY' .or. trim(fileformat)=='binary') then
       is_binary = .true.
    else
       print "(2x,a)",trim(fileformat)
    endif
 endif

 ! read dataset type
 read(iu,*,iostat=ierr) dataset,datatype
 if (ierr /= 0) then
    print*,trim(dataset),trim(datatype)
    print*,' ERROR reading vtk dataset type ',ierr
 else
    if (trim(datatype) /= 'UNSTRUCTURED_GRID') then
       print "(a)",' ERROR: got '//trim(datatype)//' but splash reader only supports UNSTRUCTURED_GRID'
       ierr = 1
       return
    endif
 endif
 close(iu)

 if (is_binary) then
    call read_vtk_legacy_binary(filename,npart,ncolstep,ierr)
 else
    print "(a)",' ERROR: got fileformat='//trim(fileformat)//' but splash reader only supports BINARY'
 endif

end subroutine read_vtk_header

!----------------------------------------------------
! read vtk legacy binary format
!----------------------------------------------------
subroutine read_vtk_legacy_binary(filename,npart,ncolstep,ierr,dat)
 use settings_data, only:debugmode
 use byteswap,      only:bs
 use labels,        only:iamvec
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: npart,ncolstep,ierr
 real,             intent(out), optional :: dat(:,:)
 character(len=32)  :: dataset,datatype,tmp
 character(len=256) :: line
 character(len=1)   :: newline
 real(kind=4), allocatable :: xyz(:,:),rval(:,:)
 real(kind=8), allocatable :: r8val(:,:)
 integer :: iu,nline,n,i,nvec,nfields,npts

 ! open the file with stream access and use get_next_line to parse the ascii lines
 open(newunit=iu,file=filename,status='old',action='read',form='unformatted',iostat=ierr,access='stream')
 nline = 0
 ncolstep = 0
 iamvec(:) = 0
 do while(ierr == 0)
    call get_next_line(iu,line,n,nline,ierr)

    ! Process the line as needed, e.g., print it
    select case (line(1:6))
    case('POINTS')
       read(line(1:n),*,iostat=ierr) dataset,npart,datatype
       !print*,trim(dataset),npart,trim(datatype)
       allocate(xyz(3,npart),rval(1,npart))
       read(iu,iostat=ierr) xyz,newline
       if (present(dat)) then
          dat(1:npart,1:3) = bs(transpose(xyz)) ! byteswap from big endian -> little
          !print*,' x vals = ',dat(1:4,1)
          !print*,' y vals = ',dat(1:4,2)
          !print*,' z vals = ',dat(1:4,3)
       endif
       deallocate(xyz,rval)
       ncolstep = ncolstep + 3
    case('FIELD ')
       nvec = 0
       read(line(1:n),*,iostat=ierr) tmp,tmp,nfields
       !print*,' number of fields = ',nfields
       do i=1,nfields
          call get_next_line(iu,line,n,nline,ierr)
          if (ierr /= 0 .or. n < 1) exit
          read(line(1:n),*,iostat=ierr) tagarr(ncolstep+1),nvec,npts,datatype
          !print*,'GOT: ',tagarr(ncolstep+1),': nvec=',nvec,'npts=',npts,trim(datatype)
          select case(trim(datatype))
          case('float')
             ! read/convert floats to working precision
             allocate(rval(nvec,npts))
             read(iu,iostat=ierr) rval,newline
             if (debugmode) print*,bs(rval(1,1:min(npts,10)))
             if (npts==npart) then
                if (present(dat)) dat(1:npart,ncolstep+1:ncolstep+nvec) = bs(transpose(rval)) ! byteswap from big endian -> little
                if (nvec > 1) iamvec(ncolstep+1:ncolstep+nvec) = ncolstep
                ncolstep = ncolstep + nvec
             else
                if (present(dat)) print*,'WARNING: skipping float '//trim(tagarr(ncolstep+1))//' as npts /= npart'
             endif
             deallocate(rval)
          case('double')
             ! read/convert doubles to working precision
             allocate(r8val(nvec,npts))
             read(iu,iostat=ierr) r8val,newline
             if (debugmode) print*,bs(r8val(1,1:min(npts,10)))
             if (npts==npart) then
                if (present(dat)) dat(1:npart,ncolstep+1:ncolstep+nvec) = bs(transpose(r8val)) ! byteswap from big endian -> little
                if (nvec > 1) iamvec(ncolstep+1:ncolstep+nvec) = ncolstep
                ncolstep = ncolstep + nvec
             else
                if (present(dat)) print*,'WARNING: skipping double '//trim(tagarr(ncolstep+1))//' as npts /= npart'
             endif
             deallocate(r8val)
          case('int')
             ! just silently skip integers
             allocate(rval(nvec,npts))
             read(iu,iostat=ierr) rval,newline
             deallocate(rval)
          case default
             if (present(dat)) print*,'Warning: unknown data type '//trim(datatype)//' assuming 4-bytes and skipping'
             allocate(rval(nvec,npts))
             read(iu,iostat=ierr) rval,newline
             deallocate(rval)
          end select
       enddo
    end select
 enddo

end subroutine read_vtk_legacy_binary

!------------------------------------------------------
! read a line of ascii text from an unformatted stream
!------------------------------------------------------
subroutine get_next_line(iu,line,n,nline,ierr)
 use settings_data, only:debugmode
 integer,          intent(in)    :: iu
 character(len=*), intent(out)   :: line
 integer,          intent(inout) :: nline
 integer,          intent(out)   :: n,ierr

 n = 1
 line = ''
 do while(n < len(line))
    read(iu, iostat=ierr) line(n:n)
    if (ierr /= 0) exit
    if (line(n:n) == char(10)) exit
    n = n + 1
 end do
 n = n - 1  ! n is the line length

 if (ierr == 0) then
    nline = nline + 1
    if (debugmode) print "(a,i2,a)",' DEBUG: ',nline,' : '//trim(line(1:n))
 endif

end subroutine get_next_line

!------------------------------------------------------------
! set labels for each column of data
!------------------------------------------------------------
subroutine set_labels_vtk
 use labels,         only:label,labeltype,ix,ipmass,ih,irho,labelvec,iamvec,&
                          set_vector_labels,label_synonym
 use settings_data,  only:ndim,ndimV,ntypes,UseTypeInRenderings,ncolumns
 use geometry,       only:labelcoord
 integer :: i

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_fits ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_fits ***'
    return
 endif

 ix(1) = 1
 ix(2) = 2
 if (ndim >= 3) ix(3) = 3
 !ipmass = ndim+3
 do i=4,ncolumns
    label(i) = trim(tagarr(i))
    if (label_synonym(label(i))=='density') irho = i
    if (label_synonym(label(i))=='h') ih = i
 enddo
 call set_vector_labels(ncolumns,ndimV,iamvec,labelvec,label,labelcoord(:,1))

 ! set labels of the quantities read in
 if (ix(1) > 0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
 if (irho > 0)    label(irho)       = 'density'
 if (ipmass > 0)  label(ipmass)     = 'particle mass'
 if (ih > 0)      label(ih)         = 'h'

 ! set labels for each particle type
 ntypes = 1
 labeltype(1) = 'gas'
 UseTypeInRenderings(:) = .true.

end subroutine set_labels_vtk

end module readdata_vtk
