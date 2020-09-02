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
!  Copyright (C) 2005-2011 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE DRAGON CODE
! HANDLES BOTH ASCII AND BINARY FILES
!
! THE FOLLOWING ENVIRONMENT VARIABLES AFFECT THIS FORMAT:
!
! DSPLASH_EXTRACOLS : set to number of extra columns to read (after itype)
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxpart,maxplot,maxstep) : main data array
!
! npartoftype(maxstep): number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!
! Partial data read implemented means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------

module unit_constants
 integer, parameter :: DP = selected_real_kind(p=15) ! double precision

 ! Length units in metres
 real(kind=DP),parameter :: r_pc    = 3.08568E16_DP     ! parsec
 real(kind=DP),parameter :: r_au    = 1.49597870E11_DP  ! astronomical unit
 real(kind=DP),parameter :: r_sun   = 6.96E8_DP         ! solar radius
 real(kind=DP),parameter :: r_earth = 6.371E6_DP        ! Earth radius

 ! Mass units in kilograms
 real(kind=DP),parameter :: m_sun   = 1.98892E30_DP     ! solar mass
 real(kind=DP),parameter :: m_jup   = 1.8986E27_DP      ! Jupiter mass
 real(kind=DP),parameter :: m_earth = 5.9736E24_DP      ! Earth mass

 ! Time units in seconds
 real(kind=DP),parameter :: myr = 3.1556952E13_DP       ! megayear
 real(kind=DP),parameter :: yr  = 3.1556952E7_DP        ! year
 real(kind=DP),parameter :: day = 8.64E4_DP             ! day

end module unit_constants

module readdata_dragon
 implicit none
 
 public :: read_data_dragon, set_labels_dragon
 
 private 
contains

subroutine read_data_dragon(rootname,istepstart,ipos,nstepsread)
 use particle_data, only:dat,iamtype,npartoftype,time,gamma,maxpart,maxcol,maxstep
 use params
 use settings_data, only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,ntypes
 use settings_units, only:unitzintegration, unit_interp
 use mem_allocation, only:alloc
 use labels, only:label,labeltype,labelzintegration
 use system_utils, only:ienvironment
 implicit none
 integer, intent(in) :: istepstart,ipos
 integer, intent(out) :: nstepsread
 character(len=*), intent(in) :: rootname
 character(len=len(rootname)+10) :: datfile
 integer, parameter :: iunit = 16
 integer :: i,j,icol,ierr,iambinaryfile,itype
 integer :: ncolstep,npart_max,nstep_max,ntoti,nlastcol,nextracols
 logical :: iexist,reallocate,doubleprec
 character(len=11) :: fmt
 character(len=50) :: string

 integer :: nei_want,nei_min,nmax
 integer, dimension(:), allocatable :: iparttype
 integer, dimension(20) :: idata
 real, dimension(50) :: rdata
 real(doub_prec), dimension(50) :: rdatadb
 real :: timetemp,gammatemp,runit,massunit,sinksoft,sinkrad

 nstepsread = 0
 if (len_trim(rootname) > 0) then
    datfile = trim(rootname)
 else
    print*,' **** no data read **** '
    return
 endif
!
!--check if first data file exists
!
 inquire(file=datfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(datfile)//': file not found ***'
    return
 endif
!
!--set parameters which do not vary between timesteps
!
 ndim = 3
 ndimV = 3
 ncolstep = ndim + ndimV + 4 ! pos x 3, vel x 3, temp, h, rho, mass
 nlastcol = ncolstep
 nextracols = ienvironment('DSPLASH_EXTRACOLS',0)
 if (nextracols > 0 .and. nextracols <= 99) then
    print "(a,i2,a)",' ASSUMING ',nextracols,' EXTRA COLUMNS BEYOND ITYPE'
    ncolstep = ncolstep + nextracols
 else
    nextracols = 0
 endif
 ncolumns = ncolstep
 call set_labels_dragon
!
!--read data from snapshots
!
 j = istepstart

 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
 !
 !--open data file and read data
 !
 !
 !--determine whether file is binary or ascii, open it and read the header
 !
 inquire(file=datfile,form=fmt)
 !print*,'fmt = ',fmt

!   select case(trim(adjustl(fmt)))
!   case('UNFORMATTED')
!      iambinaryfile = 1
!      open(unit=iunit,file=datfile,status='old',form='unformatted',iostat=ierr)
!      write (6,*) "Compiler identified as UNFORMATTED"
!   case('FORMATTED')
!      iambinaryfile = 0
!      open(unit=iunit,file=datfile,status='old',form='formatted',iostat=ierr)
!      write (6,*) "Compiler identified as FORMATTED"
!   case default
 !--if compiler cannot distinguish the two, try binary first, then ascii
 iambinaryfile = -1
 open(unit=iunit,file=datfile,status='old',form='unformatted',iostat=ierr)
!   end select

 if (ierr /= 0) then
    print "(a)",'*** ERROR OPENING '//trim(datfile)//' ***'
    return
 endif
 !
 !--read the file header
 !  try binary format first, and if unsuccessful try ascii
 !
 doubleprec = .true.
 if (iambinaryfile==1) then
    print "(a)",' reading binary dragon format '
    call read_dragonheader_binary(iunit,ierr)
 elseif (iambinaryfile==0) then
    print "(a)",' reading ascii dragon format '
    call read_dragonheader_ascii(iunit,ierr,iambinaryfile)
 else
    call read_dragonheader_binary(iunit,ierr)
    if (ierr==0) then
       !--if successful binary header read, file is doubleprec binary
       iambinaryfile = 1
       print "(a)",' reading binary dragon format '
       print "(a)",' Double precision file'
    else
       !--otherwise, close binary file, and assume file is single precision binary
       doubleprec = .false.
       close(unit=iunit)
       iambinaryfile = 1
       open(unit=iunit,file=datfile,status='old',form='unformatted',iostat=ierr)
       call read_dragonheader_binary(iunit,ierr)
       if (ierr==0) then
          print "(a)",' reading binary dragon format '
          print "(a)",' Single precision file'
       else
          print "(a)",' reading ascii dragon format '
          iambinaryfile = 0
          close(unit=iunit)
          open(unit=iunit,file=datfile,status='old',form='formatted',iostat=ierr)
          call read_dragonheader_ascii(iunit,ierr,iambinaryfile)
          if (ierr/=0) then
             print "(a)",' ERROR reading ascii file header: wrong endian binary? '
             close (iunit)
             ndim = 0
             ncolumns = 0
             return
          endif
       endif
    endif
 endif
 !
 !--get values of quantities from the header
 !
 ntoti = idata(1)
 nmax = idata(3)
 nei_want = idata(10)
 nei_min = idata(12)

 !--check for errors in integer header (either from corrupt file or wrong endian)
 if (ntoti <= 0 .or. ntoti > 1.e10 .or. nmax < 0 &
      .or. nei_want < 0 .or. nei_want > 1e6 .or. nei_min < 0) then
    if (iambinaryfile==1) then
       print "(a)",' ERROR reading binary file header: wrong endian? '
    else
       print "(a)",' ERROR reading ascii file header '
    endif
    ndim = 0
    ncolumns = 0
    close(unit=iunit)
    return
 endif

 if (doubleprec) then
    timetemp = rdatadb(1)
    runit = rdatadb(21)
    massunit = rdatadb(22)
    gammatemp = rdatadb(26)
    sinksoft = rdatadb(27)
    sinkrad = rdatadb(38)
 else
    timetemp = rdata(1)
    runit = rdata(21)
    massunit = rdata(22)
    gammatemp = rdata(26)
    sinksoft = rdata(27)
    sinkrad = rdata(38)
 endif

 !--assume first that the file is single precision, check values are sensible, if not try double
!   if (iambinaryfile==1) then
!      if (timetemp < 0. .or. runit < 0. .or. massunit < 0. .or. gammatemp < 0. &
!       .or. gammatemp > 6.) then
!         print "(a)",' double precision file'
!         doubleprec = .true.
!         rewind(iunit)
!         call read_dragonheader_binary(iunit,ierr)
!      else
!         print "(a)",' single precision file'
!      endif
!   endif

 print*,'time             : ',timetemp
 print*,'gamma            : ',gammatemp
 print*,'n_total          : ',ntoti

 if (ierr /= 0) then
    if (iambinaryfile==1) then
       print "(a)",' ERROR reading real part of binary file header '
    else
       print "(a)",' ERROR reading real part of ascii file header '
    endif
    ndim = 0
    ncolumns = 0
    close(unit=iunit)
    return
 endif
 !
 !--if successfully read header, increment the nstepsread counter
 !
 nstepsread = nstepsread + 1
 !
 !-- now work out dimensionless weight unit and z integration unit
 !
 call find_weights(unit_interp,unitzintegration,labelzintegration)
 !
 !--now read data
 !
 reallocate = .false.
 npart_max = maxpart
 nstep_max = max(maxstep,1)

 if (ntoti > maxpart) then
    reallocate = .true.
    if (maxpart > 0) then
       ! if we are reallocating, try not to do it again
       npart_max = int(1.1*ntoti)
    else
       ! if first time, save on memory
       npart_max = int(ntoti)
    endif
 endif
 if (j >= maxstep .and. j /= 1) then
    nstep_max = j + max(10,INT(0.1*nstep_max))
    reallocate = .true.
 endif
 !
 !--reallocate memory for main data array
 !
 if (reallocate .or. .not.(allocated(dat))) then
    call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol),mixedtypes=.true.)
 endif
 !
 !--copy header into header arrays
 !
 npartoftype(:,j) = 0
 npartoftype(1,j) = ntoti
 time(j) = timetemp
 gamma(j) = gammatemp
 !
 !--read particle data
 !
 if (ntoti > 0) then
    if (iambinaryfile==1) then
       call read_dragonbody_binary(iunit,ierr)
    else
       call read_dragonbody_ascii(iunit,ierr)
    endif
 else
    ntoti = 0
    npartoftype(1,i) = 0
    dat(:,:,i) = 0.
 endif

 if (allocated(iamtype)) then
    !--relabel particle types
    call set_types(iamtype(:,j),ntoti,npartoftype(:,j))
 endif
 if (any(npartoftype(2:,j) /= 0)) then
    do itype=1,ntypes
       if (npartoftype(itype,j) > 0) then
          string = ' '
          write(string,"(a)") 'n_'//trim(labeltype(itype))
          write(string(18:len(string)),"(a)") ':'
          print*,trim(string),' ',npartoftype(itype,j)
       endif
    enddo
 endif
!
!--set flag to indicate that only part of this file has been read
!
 if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--close data file and return
!
 if (allocated(iparttype)) deallocate(iparttype)
 close(unit=iunit)

 return

contains

!----------------------------------------------------
! binary header read
!----------------------------------------------------
subroutine read_dragonheader_binary(iunitb,ierr)
 implicit none
 integer, intent(in) :: iunitb
 integer, intent(out) :: ierr

 read(iunitb,end=55,iostat=ierr) idata
 if (doubleprec) then
    read(iunitb,end=55,iostat=ierr) rdatadb
 else
    read(iunitb,end=55,iostat=ierr) rdata
 endif

 return

55 continue
 !print "(a)",' ERROR: end of file in binary header read'
 ierr = -1
 return

end subroutine read_dragonheader_binary

!----------------------------------------------------
! ascii header read
!----------------------------------------------------
subroutine read_dragonheader_ascii(iunita,ierr,iwarn)
 implicit none
 integer, intent(in) :: iunita,iwarn
 integer, intent(out) :: ierr

 do i=1,size(idata)
    read(iunita,*,end=55,iostat=ierr) idata(i)
 enddo

 do i=1,size(rdata)
    read(iunita,*,end=55,iostat=ierr) rdata(i)
 enddo
 doubleprec = .false.
 return

55 continue
 if (iwarn >= 0) print "(a)",' ERROR: end of file in binary header read'
 ierr = -1
 return

end subroutine read_dragonheader_ascii

!----------------------------------------------------
! binary body read
!----------------------------------------------------
subroutine read_dragonbody_binary(iunitb,ierr)
 implicit none
 integer, intent(in) :: iunitb
 integer, intent(out) :: ierr
 real(doub_prec), dimension(:,:), allocatable :: dummyx
 real(doub_prec), dimension(:), allocatable :: dummy
 integer, dimension(:), allocatable :: idumtype
 integer :: icol

 if (doubleprec .and. any(required(1:ndim+ndimV))) then
    allocate(dummyx(3,ntoti),stat=ierr)
    if (ierr /= 0) then
       print *,' ERROR allocating memory'
       goto 56
    endif
 endif

 !--positions
 if (any(required(1:ndim))) then
    if (doubleprec) then
       read(iunitb,end=55,iostat=ierr) dummyx(1:ndim,1:ntoti)
       do i=1,ntoti
          dat(i,1:ndim,j) = real(dummyx(1:ndim,i))
       enddo
    else
       read(iunitb,end=55,iostat=ierr) (dat(i,1:ndim,j),i=1,ntoti)
    endif
    if (ierr /= 0) print*,' WARNING: errors reading positions '
 else
    read(iunitb,end=55,iostat=ierr)
    if (ierr /= 0) print*,' WARNING: error skipping positions '
 endif

 !--velocities
 if (any(required(ndim+1:ndim+ndimV))) then
    if (doubleprec) then
       read(iunitb,end=55,iostat=ierr) dummyx(1:ndimV,1:ntoti)
       do i=1,ntoti
          dat(i,ndim+1:ndim+ndimV,j) = real(dummyx(1:ndimV,i))
       enddo
    else
       read(iunitb,end=55,iostat=ierr) (dat(i,ndim+1:ndim+ndimV,j),i=1,ntoti)
    endif
    if (ierr /= 0) print*,' WARNING: errors reading velocities '
 else
    read(iunitb,end=55,iostat=ierr)
    if (ierr /= 0) print*,' WARNING: error skipping velocities '
 endif

 if (doubleprec .and. any(required(ndim+ndimV+1:ncolstep))) then
    allocate(dummy(ntoti),stat=ierr)
    if (ierr /= 0) then
       print*,' ERROR allocating memory'
       goto 56
    endif
 endif

 !--the rest
 do icol = ndim+ndimV+1,nlastcol
    if (required(icol)) then
       if (doubleprec) then
          read(iunitb,end=55,iostat=ierr) dummy(1:ntoti)
          dat(1:ntoti,icol,j) = real(dummy(1:ntoti))
       else
          read(iunitb,end=55,iostat=ierr) dat(1:ntoti,icol,j)
       endif
       if (ierr /= 0) print*,' WARNING: errors reading '//trim(label(icol))
    else
       read(iunitb,end=55,iostat=ierr)
       if (ierr /= 0) print*,' WARNING: error skipping '//trim(label(icol))
    endif
 enddo

 if (size(iamtype(:,j)) > 1) then
    allocate(idumtype(ntoti),stat=ierr)
    if (ierr /= 0) then
       print*,'error reading type, assuming all gas'
       iamtype(1:ntoti,j) = 1
    else
       read(iunitb,end=55,iostat=ierr) idumtype(1:ntoti)
       iamtype(1:ntoti,j) = idumtype(1:ntoti)
    endif
    deallocate(idumtype)
    if (ierr /= 0) print*,' WARNING: error reading itype'
 endif

 !--extra columns beyond itype
 do icol = nlastcol+1,nlastcol+nextracols
    if (required(icol)) then
       if (doubleprec) then
          read(iunitb,end=55,iostat=ierr) dummy(1:ntoti)
          dat(1:ntoti,icol,j) = real(dummy(1:ntoti))
       else
          read(iunitb,end=55,iostat=ierr) dat(1:ntoti,icol,j)
       endif
       if (ierr /= 0) print*,' WARNING: errors reading '//trim(label(icol))
    else
       read(iunitb,end=55,iostat=ierr)
       if (ierr /= 0) print*,' WARNING: error skipping '//trim(label(icol))
    endif
 enddo

 if (allocated(dummyx)) deallocate(dummyx)
 if (allocated(dummy)) deallocate(dummy)
 return

55 continue
 print "(a)",' ERROR: end of file in binary read'
56 continue
 ierr = -1

 if (allocated(dummyx)) deallocate(dummyx)
 if (allocated(dummy)) deallocate(dummy)
 return

end subroutine read_dragonbody_binary

!----------------------------------------------------
! ascii body read
!----------------------------------------------------
subroutine read_dragonbody_ascii(iunita,ierr)
 implicit none
 integer, intent(in) :: iunita
 integer, intent(out) :: ierr
 integer :: nerr,idumtype

 !--positions
 nerr = 0
 do i=1,ntoti
    read(iunita,*,end=55,iostat=ierr) dat(i,1:ndim,j)
    if (ierr /= 0) nerr = nerr + 1
 enddo
 if (nerr > 0) print*,' WARNING: ',nerr,' errors reading positions '

 !--velocities
 nerr = 0
 do i=1,ntoti
    read(iunita,*,end=55,iostat=ierr) dat(i,ndim+1:ndim+ndimV,j)
    if (ierr /= 0) nerr = nerr + 1
 enddo
 if (nerr > 0) print*,' WARNING: ',nerr,' errors reading velocities '

 !--the rest
 if (any(required(ndim+ndimV+1:nlastcol))) then
    do icol = ndim+ndimV+1,nlastcol
       nerr = 0
       do i=1,ntoti
          read(iunita,*,end=55,iostat=ierr) dat(i,icol,j)
          if (ierr /= 0) nerr = nerr + 1
       enddo
       if (nerr > 0) print*,' WARNING: ',nerr,' errors reading '//trim(label(icol))
    enddo
 endif

 !--particle type
 if (size(iamtype(:,j)) > 1) then
    nerr = 0
    do i=1,ntoti
       read(iunita,*,end=55,iostat=ierr) idumtype
       iamtype(i,j) = idumtype
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print*,' WARNING: ',nerr,' errors reading itype'
 endif

 !--the rest
 if (any(required(nlastcol+1:nlastcol+nextracols))) then
    do icol = nlastcol+1,nlastcol+nextracols
       nerr = 0
       do i=1,ntoti
          read(iunita,*,end=55,iostat=ierr) dat(i,icol,j)
          if (ierr /= 0) nerr = nerr + 1
       enddo
       if (nerr > 0) print*,' WARNING: ',nerr,' errors reading '//trim(label(icol))
    enddo
 endif

 return

55 continue
 print "(a)",' ERROR: end of file in ascii read'
 ierr = -1
 return

end subroutine read_dragonbody_ascii

!----------------------------------------------------
! translate types into order (for old dragon read)
!----------------------------------------------------
subroutine set_types(itypei,ntotal,noftype)
 implicit none
 integer(kind=int1), dimension(:), intent(inout) :: itypei
 integer, intent(in) :: ntotal
 integer, dimension(:), intent(out) :: noftype
 integer :: ngas,nsink,nbnd,ncloud,nsplit,nunknown,nstar

!--types
!  1 gas
! -1 sink
! -2 star
!  6 boundary (fixed)
!  9 intercloud (hydro only)
!  4 split particle (obsolete?)
!
!--we translate these into
!  1 gas
!  2 boundary
!  3 sink
!  4 intercloud
!  5 split
!  6 unknown / the rest
!  7 star
!
 ngas = 0
 nsink = 0
 nbnd = 0
 ncloud = 0
 nsplit = 0
 nunknown = 0
 nstar = 0
 do i=1,ntotal
    select case(itypei(i))
    case(1)
       ngas = ngas + 1
       itypei(i) = 1
    case(-1)
       nsink = nsink + 1
       itypei(i) = 3
    case(6)
       nbnd = nbnd + 1
       itypei(i) = 2
    case(9)
       ncloud = ncloud + 1
       itypei(i) = 4
    case(4)
       nsplit = nsplit + 1
       itypei(i) = 5
    case(-2)
       nstar = nstar + 1
       itypei(i) = 6
    case default
       nunknown = nunknown + 1
!        itypei(i) = 6
       write (6,*) "Unknown particle type ", itypei(i), "!!"
       stop
    end select
 enddo

 noftype(1) = ngas
 noftype(2) = nbnd
 noftype(3) = nsink
 noftype(4) = ncloud
 noftype(5) = nsplit
 noftype(6) = nstar
 if (sum(noftype(1:6)) /= ntotal) then
    print "(a)",' INTERNAL ERROR setting number in each type in dragon read'
 endif

 return
end subroutine set_types

end subroutine read_data_dragon

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels_dragon
 use labels, only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,ih,irho
 use params
 use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings
 use geometry, only:labelcoord
 implicit none
 integer :: i

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_dragon ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_dragon ***'
    return
 endif

 do i=1,ndim
    ix(i) = i
 enddo
 ivx = ndim+1
 label(ivx+ndimV) = 'temperature'
 ih = ivx + ndimV + 1
 irho = ih + 1        ! location of rho in data array
 ipmass = irho + 1
 !
 !--set labels of the quantities read in
 !
 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 label(irho) = 'density'
 !label(iutherm) = 'u'
 label(ipmass) = 'particle mass'
 label(ih) = 'h'
 !
 !--set labels for vector quantities
 !
 iamvec(ivx:ivx+ndimV-1) = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'
 do i=1,ndimV
    label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
 enddo

 !--set labels for each particle type
 !
 ntypes = 6
 labeltype(1) = 'gas'
 labeltype(2) = 'boundary'
 labeltype(3) = 'sink'
 labeltype(4) = 'cloud'
 labeltype(5) = 'split'
 labeltype(6) = 'star'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .true.
 UseTypeInRenderings(3) = .false.
 UseTypeInRenderings(4) = .true.
 UseTypeInRenderings(5) = .true.
 UseTypeInRenderings(6) = .false.

!-----------------------------------------------------------
 return
end subroutine set_labels_dragon

subroutine find_weights(out_unit_interp,out_unitzintegration,out_labelzintegration)
 use labels, only:ipmass,ih,irho
 use params
 use settings_data, only:ndim
 use unit_constants
 use system_commands, only:get_environment
 implicit none

 real(doub_prec), intent(out)      :: out_unit_interp
 real, intent(out)                 :: out_unitzintegration
 character(len=20), intent(out)    :: out_labelzintegration
 real(doub_prec)                   :: dm_unit, dh_unit, drho_unit, dr_unit
 logical                           :: do_dimweight, do_zintegration
 character(len=20)                 :: rho_length_label
 real(doub_prec)                   :: rho_length
 character(len=20) :: r_unit                  ! length unit
 character(len=20) :: m_unit                  ! mass unit
 character(len=20) :: rho_unit                ! density unit
 character(len=20) :: h_unit                  ! smoothing length unit

 call get_environment("DRAGON_R_UNIT",r_unit)
 call get_environment("DRAGON_M_UNIT",m_unit)
 call get_environment("DRAGON_RHO_UNIT",rho_unit)
 call get_environment("DRAGON_H_UNIT",H_unit)

 out_unit_interp = 1.0
 out_unitzintegration = 1.0
 out_labelzintegration = ""

 do_dimweight = .TRUE.
 do_zintegration = .TRUE.

! Length unit in S.I. units (m)
 if (r_unit=="") then
    print*,'No positions or no position units!'
    print*,'Set environment variable DRAGON_R_UNIT to:'
    print*,'  pc, au, r_sun, r_earth, km, m, cm or 1 (dimensionless)'
    do_zintegration = .FALSE.
    dr_unit = 1._DP
 elseif (r_unit=="pc") then
    dr_unit = r_pc
 elseif (r_unit=="au") then
    dr_unit = r_au
 elseif (r_unit=="r_sun") then
    dr_unit = r_sun
 elseif (r_unit=="r_earth") then
    dr_unit = r_earth
 elseif (r_unit=="km") then
    dr_unit = 1000.0_DP
 elseif (r_unit=="m") then
    dr_unit = 1.0_DP
 elseif (r_unit=="cm") then
    dr_unit = 0.01_DP
 elseif (r_unit=="1") then
    dr_unit = 1._DP
 else
    print*,'Unknown position unit ', r_unit, '!'
    do_zintegration = .FALSE.
    dr_unit = 1._DP
 endif

! Length unit in S.I. units (m)
 if (h_unit=="") then
    print*,'No smoothing lengths or no smoothing length units!'
    print*,'Set environment variable DRAGON_H_UNIT to: '
    print*,'  pc, au, r_sun, r_earth, km, m, cm or 1 (dimensionless)'
    do_dimweight = .FALSE.
    dh_unit = 1._DP
 elseif (h_unit=="pc") then
    dh_unit = r_pc
 elseif (h_unit=="au") then
    dh_unit = r_au
 elseif (h_unit=="r_sun") then
    dh_unit = r_sun
 elseif (h_unit=="r_earth") then
    dh_unit = r_earth
 elseif (h_unit=="km") then
    dh_unit = 1000.0_DP
 elseif (h_unit=="m") then
    dh_unit = 1.0_DP
 elseif (h_unit=="cm") then
    dh_unit = 0.01_DP
 elseif (h_unit=="1") then
    dh_unit = 1._DP
 else
    print*,'Unknown smoothing length unit ', h_unit, '!'
    do_dimweight = .FALSE.
    dh_unit = 1._DP
 endif

! Mass units in S.I. units (kg)
 if (m_unit=="") then
    print*,'No masses or no mass units!'
    print*,'Set environment variable DRAGON_M_UNIT to:'
    print*,'  m_sun, m_jup, m_earth, kg, g or 1 (dimensionless)'
    do_dimweight = .FALSE.
    dm_unit = 1._DP
 elseif (m_unit=="m_sun") then
    dm_unit = m_sun
 elseif (m_unit=="m_jup") then
    dm_unit = m_jup
 elseif (m_unit=="m_earth") then
    dm_unit = m_earth
 elseif (m_unit=="kg") then
    dm_unit = 1._DP
 elseif (m_unit=="g") then
    dm_unit = 1.0E-3_DP
 elseif (m_unit=="1") then
    dm_unit = 1._DP
 else
    print*,'Unknown mass unit ', m_unit, '!'
    do_dimweight = .FALSE.
    dm_unit = 1._DP
 endif

 ! Density units in S.I. units (i.e. kg/m^3)
 if (rho_unit=="") then
    print*,'No densities or no density units!'
    print*,'Set environment variable DRAGON_RHO_UNIT to:'
    if (ndim==3) print*,'  m_sun_pc3, kg_m3, g_cm3 or 1 (dimensionless)'
    if (ndim==2) print*,'  m_sun_pc2, kg_m2, g_cm2 or 1 (dimensionless)'
    if (ndim==1) print*,'  1 (dimensionless)'
    do_dimweight = .FALSE.
    do_zintegration = .FALSE.
    rho_length = 1._DP
 elseif (rho_unit=="m_sun_pc3") then
    drho_unit = m_sun / (r_pc**3)
    rho_length = r_pc
    rho_length_label = "pc"
 elseif (rho_unit=="m_sun_pc2") then
    drho_unit = m_sun / (r_pc**2)
    rho_length = r_pc
    rho_length_label = "pc"
 elseif (rho_unit=="kg_m3") then
    drho_unit = 1.0_DP
    rho_length = 1.0_DP
    rho_length_label = "m"
 elseif (rho_unit=="kg_m2") then
    drho_unit = 1.0_DP
    rho_length = 1.0_DP
    rho_length_label = "m"
 elseif (rho_unit=="g_cm3") then
    drho_unit = 1.0E3_DP
    rho_length = 0.01_DP
    rho_length_label = "cm"
 elseif (rho_unit=="g_cm2") then
    drho_unit = 10.0_DP
    rho_length = 0.01_DP
    rho_length_label = "cm"
 elseif (rho_unit=="1") then
    drho_unit = 1._DP
    rho_length = 1._DP
    rho_length_label = ""
 else
    print*,'Unknown density unit ', rho_unit, '!'
    do_dimweight = .FALSE.
    do_zintegration = .FALSE.
    rho_length = 1._DP
 endif

 if (do_dimweight) then
    out_unit_interp = dm_unit/(drho_unit*dh_unit**ndim)
 else
    print*,'Cannot create dimensionless weight'
    print*,'(unnormalised rendered plots may be incorrect)'
 endif

 if (do_zintegration) then
    out_unitzintegration = dr_unit / rho_length
    out_labelzintegration = rho_length_label
 else
    print*,'Cannot set unitzintegration'
    print*,'(column density plots may be incorrect)'
 endif

 return
end subroutine find_weights
end module readdata_dragon
