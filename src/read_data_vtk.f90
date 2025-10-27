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
! Supports both UNSTRUCTURED_GRID and STRUCTURED_GRID formats
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
 integer               :: ncolstep,npart,nsteps_to_read,nheader
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
 nheader = 2
 call set_labels_vtk
 call read_vtk_legacy_binary(trim(datfile),npart,ncolstep,ierr,dat(:,:,i),headertags,headervals(:,i),nheader)
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
    if (trim(datatype) /= 'UNSTRUCTURED_GRID' .and. trim(datatype) /= 'STRUCTURED_GRID') then
       print "(a)",' ERROR: got '//trim(datatype)//' but splash reader only supports UNSTRUCTURED_GRID or STRUCTURED_GRID'
       ierr = 1
       return
    endif
 endif
 close(iu)

 if (is_binary) then
    ! For both UNSTRUCTURED_GRID and STRUCTURED_GRID, we need to read the file
    ! to determine npart and ncolstep from the actual data sections
    call read_vtk_legacy_binary(filename,npart,ncolstep,ierr)
 else
    print "(a)",' ERROR: got fileformat='//trim(fileformat)//' but splash reader only supports BINARY'
 endif

end subroutine read_vtk_header

!----------------------------------------------------
! read vtk legacy binary format
!----------------------------------------------------
subroutine read_vtk_legacy_binary(filename,npart,ncolstep,ierr,dat,headertags,headervals,nheader)
 use settings_data, only:debugmode
 use byteswap,      only:bs
 use labels,        only:iamvec,ih,ipmass,irho
 use geometry,      only:igeom_cartesian,igeom_spherical,coord_transform
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: npart,ncolstep,ierr
 real,             intent(out), optional :: dat(:,:)
 character(len=*), intent(inout), optional :: headertags(:)
 real,             intent(inout), optional :: headervals(:)
 integer,          intent(inout), optional :: nheader
 character(len=32)  :: dataset,datatype,tmp
 character(len=256) :: line
 character(len=1)   :: newline
 real(kind=4), allocatable :: xyz(:,:),rval(:,:)
 real(kind=8), allocatable :: r8val(:,:)
 real(kind=4), allocatable :: x1c(:),x2c(:),x3c(:)
 integer :: iu,nline,n,i,nvec,nfields,npts,ispace
 integer :: nx,ny,nz,grid_type,ncells,icell,j,k,nvertices
 logical :: found_dimensions,found_points
 logical :: in_cell_data,in_point_data
 real    :: xin(3),xout(3)
 integer :: igeom

 ! open the file with stream access and use get_next_line to parse the ascii lines
 open(newunit=iu,file=filename,status='old',action='read',form='unformatted',iostat=ierr,access='stream')
 nline = 0
 ncolstep = 0
 iamvec(:) = 0
 found_dimensions = .false.
 found_points = .false.
 in_cell_data = .false.
 in_point_data = .false.
 grid_type = 0  ! 0=unknown, 1=unstructured, 2=structured
 ncells = 0

 do while(ierr == 0)
    call get_next_line(iu,line,n,nline,ierr)

    ! Process the line as needed, e.g., print it
    select case (line(1:6))
    case('DIMENS')
       ! Handle DIMENSIONS section for STRUCTURED_GRID
       read(line(1:n),*,iostat=ierr) dataset,nx,ny,nz
       if (ierr == 0) then
          nvertices = nx * ny * nz
          ncells = (nx-1) * (ny-1) * (nz-1)  ! number of cells in structured grid
          npart = ncells  ! use cells as particles for structured grids
          found_dimensions = .true.
          grid_type = 2  ! structured grid
          if (.not.present(headervals)) print "(a,i0,' x ',i0,' x ',i0,' ncells = ',i0)",&
                ' STRUCTURED_GRID: ',nx-1,ny-1,nz-1,ncells
       else
          print*,' ERROR reading DIMENSIONS section'
       endif
    case('POINTS')
       read(line(1:n),*,iostat=ierr) dataset,npts,datatype
       found_points = .true.

       ! For STRUCTURED_GRID, verify that npts matches our calculated npart
       if (found_dimensions) then
          if (npts /= nvertices) then
             print*,' ERROR: POINTS count (',npts,') does not match DIMENSIONS product (',nvertices,')'
             ierr = 1
             exit
          endif
       else
          ! For UNSTRUCTURED_GRID, use the npts from POINTS section
          npart = npts
          grid_type = 1  ! unstructured grid
          if (debugmode) print*,' UNSTRUCTURED_GRID: npart = ',npart
       endif

       !print*,trim(dataset),npart,trim(datatype)
       allocate(xyz(3,npts),rval(1,npts))
       read(iu,iostat=ierr) xyz,newline
       if (ierr /= 0) then
          print*,' ERROR reading POINTS data'
          deallocate(xyz,rval)
          exit
       endif
       if (present(dat) .and. npart==npts) then
          dat(1:npart,1:3) = bs(transpose(xyz)) ! byteswap from big endian -> little
       endif
       deallocate(xyz,rval)
       ncolstep = ncolstep + 3
    case('FIELD ')
       nvec = 0
       read(line(1:n),*,iostat=ierr) tmp,tmp,nfields
       if (debugmode) print*,' number of fields = ',nfields
       do i=1,nfields
          ! get the next line
          call get_next_line(iu,line,n,nline,ierr)
          if (ierr /= 0 .or. n < 1) exit
          ! extract the tag and number of vectors and points
          ispace = index(line,' ')                ! position of first space
          tagarr(ncolstep+1) = line(1:ispace-1)  ! tag is everything up to first space
          read(line(ispace:n),*,iostat=ierr) nvec,npts,datatype
          !print*,'LINE:',line(1:n)!
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
             elseif (npts==1 .and. present(headervals) .and. present(headertags) .and. present(nheader)) then
                ! put single scalar values in the header
                call add_to_header(headervals,headertags,nheader,real(bs(rval(1,1))),tagarr(ncolstep+1))
             elseif (trim(tagarr(ncolstep+1)) == 'X1C_NATIVE_COORDINATES') then
                if (.not. allocated(x1c)) allocate(x1c(npts))
                if (ierr == 0) x1c = bs(rval(1,:))
             elseif (trim(tagarr(ncolstep+1)) == 'X2C_NATIVE_COORDINATES') then
                if (.not. allocated(x2c)) allocate(x2c(npts))
                if (ierr == 0) x2c = bs(rval(1,:))
             elseif (trim(tagarr(ncolstep+1)) == 'X3C_NATIVE_COORDINATES') then
                if (.not. allocated(x3c)) allocate(x3c(npts))
                if (ierr == 0) x3c = bs(rval(1,:))
             else
                if (present(dat) .and. debugmode) print*,'skipping float '//trim(tagarr(ncolstep+1))//' as npts /= ncells'
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
             elseif (npts==1 .and. present(headervals) .and. present(headertags) .and. present(nheader)) then
                ! put single scalar values in the header
                call add_to_header(headervals,headertags,nheader,real(bs(r8val(1,1))),tagarr(ncolstep+1))
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
    case('CELL_D')
       ! CELL_DATA section - data applies to cells, not points
       read(line(1:n),*,iostat=ierr) dataset,npts
       in_cell_data = .true.
       in_point_data = .false.
       if (debugmode) print*,' Entering CELL_DATA section with ',npts,' cells'
       if (found_dimensions .and. npts /= ncells) then
          if (debugmode) print*,' WARNING: CELL_DATA count (',npts,') does not match expected cells (',ncells,')'
       endif
    case('POINT_')
       ! POINT_DATA section - data applies to points
       read(line(1:n),*,iostat=ierr) dataset,npts
       in_point_data = .true.
       in_cell_data = .false.
       if (debugmode) print*,' Entering POINT_DATA section with ',npts,' points'
       if (npts /= npart) then
          if (debugmode) print*,' WARNING: POINT_DATA count (',npts,') does not match npart (',npart,')'
       endif
    case('SCALAR')
       ! SCALARS section - read scalar data
       read(line(1:n),*,iostat=ierr) dataset,tmp,datatype
       if (ierr == 0) then
          tagarr(ncolstep+1) = tmp
          if (debugmode) print*,' Reading SCALARS: '//trim(tmp)//' type: '//trim(datatype)

          ! Check if LOOKUP_TABLE is on the same line or next line
          if (index(line,'LOOKUP_TABLE') > 0) then
             ! LOOKUP_TABLE is on the same line, skip it
             if (debugmode) print*,' LOOKUP_TABLE on same line, skipping'
          else
             ! Read the next line which should be LOOKUP_TABLE
             call get_next_line(iu,line,n,nline,ierr)
             if (ierr == 0 .and. index(line(1:n),'LOOKUP_TABLE') > 0) then
                if (debugmode) print*,' Skipping LOOKUP_TABLE line'
             endif
          endif

          ! Determine how many data points to read based on context
          if (in_cell_data) then
             npts = ncells
             if (debugmode) print*,' Reading SCALARS for cells: ',npts
          else
             npts = npart
             if (debugmode) print*,' Reading SCALARS for points: ',npts
          endif

          ! Read the binary data directly (not as text lines)
          select case(trim(datatype))
          case('float')
             allocate(rval(1,npts))
             read(iu,iostat=ierr) rval,newline
             ncolstep = ncolstep + 1
             if (ierr == 0) then
                if (present(dat)) dat(1:npts,ncolstep+1) = bs(rval(1,:))
                if (debugmode) print*,' Successfully read SCALARS data for '//trim(tmp)
             else
                if (debugmode) print*,' ERROR reading SCALARS binary data for '//trim(tmp)
             endif
             deallocate(rval)
          case('double')
             allocate(r8val(1,npts))
             read(iu,iostat=ierr) r8val,newline
             ncolstep = ncolstep + 1
             if (ierr == 0) then
                if (present(dat)) dat(1:npts,ncolstep+1) = bs(r8val(1,:))
                if (debugmode) print*,' Successfully read SCALARS data for '//trim(tmp)
             else
                if (debugmode) print*,' ERROR reading SCALARS binary data for '//trim(tmp)
             endif
             deallocate(r8val)
          case default
             if (present(dat)) print*,'Warning: unknown SCALARS data type '//trim(datatype)//' skipping'
          end select
       endif
    case('VECTOR')
       ! VECTORS section - read vector data
       read(line(1:n),*,iostat=ierr) dataset,tmp,datatype
       if (ierr == 0) then
          tagarr(ncolstep+1) = tmp
          if (debugmode) print*,' Reading VECTORS: '//trim(tmp)//' type: '//trim(datatype)

          ! Determine how many data points to read based on context
          if (in_cell_data) then
             npts = ncells
             if (debugmode) print*,' Reading VECTORS for cells: ',npts
          else
             npts = npart
             if (debugmode) print*,' Reading VECTORS for points: ',npts
          endif

          ! Read the data based on type
          select case(trim(datatype))
          case('float')
             allocate(rval(3,npts))
             read(iu,iostat=ierr) rval,newline
             if (ierr == 0) then
                if (present(dat)) dat(1:npts,ncolstep+1:ncolstep+3) = bs(transpose(rval))
                iamvec(ncolstep+1:ncolstep+3) = ncolstep
                ncolstep = ncolstep + 3
             endif
             deallocate(rval)
          case('double')
             allocate(r8val(3,npts))
             read(iu,iostat=ierr) r8val,newline
             if (ierr == 0) then
                if (present(dat)) dat(1:npts,ncolstep+1:ncolstep+3) = bs(transpose(r8val))
                iamvec(ncolstep+1:ncolstep+3) = ncolstep
                ncolstep = ncolstep + 3
             endif
             deallocate(r8val)
          case default
             if (present(dat)) print*,'Warning: unknown VECTORS data type '//trim(datatype)//' skipping'
          end select
       endif
    case('LOOKUP')
       ! LOOKUP_TABLE section - just skip this line
       if (debugmode) print*,' Skipping LOOKUP_TABLE section'
    end select
 enddo

 ! For structured grids, add fake h and mass columns
 if (grid_type == 2 .and. found_dimensions) then
    if (debugmode) print*,' Adding fake h and mass columns for structured grid'

    ! Add h column (cell size)
    ncolstep = ncolstep + 1
    tagarr(ncolstep) = 'h'
    ih = ncolstep

    if (present(dat)) dat(1:npart,ncolstep) = 1.0  ! Placeholder - should be calculated from actual coordinates
    if (debugmode) print*,' Added fake h column (cell size)'

    ! Add mass column (cell mass = density * cell volume)
    ncolstep = ncolstep + 1
    tagarr(ncolstep) = 'mass'
    ipmass = ncolstep
    print*,'Added mass column (cell mass)',ipmass

    if (present(dat)) dat(1:npart,ncolstep) = 1.0  ! Placeholder - should be calculated from density * cell_volume
    if (debugmode) print*,' Added fake mass column (cell mass)'

    ! Construct position data from cell coordinates
    if (present(dat) .and. allocated(x1c) .and. allocated(x2c) .and. allocated(x3c)) then
       !print*,'npts = ',npts,size(x1c),size(x2c),size(x3c),nx-1,ny-1,nz-1
       igeom = igeom_spherical
       do k = 1, nz-1
          xin(2) = x3c(k)
          do j = 1, ny-1
             xin(3) = x2c(j)
             do i = 1, nx-1
                xin(1) = x1c(i)
                icell = i + (j-1)*(nx-1) + (k-1)*(nx-1)*(ny-1)
                call coord_transform(xin,3,igeom,xout,3,igeom_cartesian)
                dat(icell,1) = xout(1)
                dat(icell,2) = xout(2)
                dat(icell,3) = xout(3)
                if (i < nx-1) then
                   dat(icell,ih) = x1c(i+1) - x1c(i)
                else
                   dat(icell,ih) = x1c(i) - x1c(i-1)
                endif
                dat(icell,ipmass) = dat(icell,ih)**3 * dat(icell,irho)
             enddo
          enddo
       enddo
       if (debugmode) print*,' Constructed position data from cell coordinates'
    endif
 endif

 ! Final validation
 if (.not. found_points) then
    print*,' ERROR: no POINTS section found in VTK file'
    ierr = 1
 elseif (grid_type == 2 .and. .not. found_dimensions) then
    print*,' ERROR: STRUCTURED_GRID format but no DIMENSIONS section found'
    ierr = 1
 elseif (npart <= 0) then
    print*,' ERROR: invalid number of points: ',npart
    ierr = 1
 endif

 close(iu)

end subroutine read_vtk_legacy_binary

!------------------------------------------------------
! add a value to the header
!------------------------------------------------------
subroutine add_to_header(headervals,headertags,nheader,value,tag)
 real, intent(inout) :: headervals(:)
 character(len=*), intent(inout) :: headertags(:)
 integer, intent(inout) :: nheader
 real, intent(in) :: value
 character(len=*), intent(in) :: tag

 nheader = nheader + 1
 headervals(nheader) = value
 headertags(nheader) = trim(tag)

end subroutine add_to_header

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
                          set_vector_labels,label_synonym,ivx
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
    if (iamvec(i) > 0 .and.trim(label(i))=='v') ivx = i
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
