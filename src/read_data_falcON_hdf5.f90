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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR HDF5 OUTPUT FROM W. DEHNEN'S FALCON CODE
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
! Columns with the 'required' flag set to false are not read 
!-------------------------------------------------------------------------
!
!  The module below contains interface routines to c functions
!  that perform the actual calls to the HDF5 libs
!
!-------------------------------------------------------------------------
module falcONhdf5read
 use params, only:maxplot,doub_prec
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none
 character(len=lenlabel), dimension(maxplot) :: blocklabel
 integer, parameter :: maxtypes = 6
 integer :: i_current_step

 interface
   ! opens a falcON HDF5 snapshot file
   subroutine open_falcON_file(filename,ierr) bind(c,name="open_falcON_file")
    import
    character(c_char), intent(in)  :: filename(*)
    integer(c_int),    intent(out) :: ierr
   end subroutine open_falcON_file
   
   ! queries whether a file is open
   function falcON_file_is_open() bind(c,name="falcON_file_is_open")
    import
    integer(c_int) :: falcON_file_is_open
   end function falcON_file_is_open

   ! closes the currently open (if any) falcON HDF5 snapshot file.
   subroutine close_falcON_file() bind(c,name="close_falcON_file")
    ! no arguments
   end subroutine close_falcON_file
   
   ! queries if there is another snapshot present the currently open file
   function num_falcON_snapshots(ierr) bind(c,name="num_falcON_snapshots")
    import
    integer(c_int) :: num_falcON_snapshots
    integer(c_int), intent(out) :: ierr
   end function num_falcON_snapshots

   ! set falcON debugging level
   subroutine set_falcON_debugging_level(level) bind(c,name="set_falcON_debugging_level")
    import
    integer(c_int), intent(in) :: level
   end subroutine set_falcON_debugging_level

   ! read falcON header
   subroutine open_falcON_snapshot(ntype,npart,ncol,dimX,dimV,time,hper,ierr) &
     bind(c,name="open_falcON_snapshot")
    import
    integer(c_int), intent(out) :: ntype,ncol,dimX,dimV,ierr
    integer(c_int), intent(out) :: npart(*)
    real(c_double), intent(out) :: time,hper(3)
   end subroutine open_falcON_snapshot
   
   ! read falcON data
   subroutine read_falcON_snapshot(ierr) bind(c,name="read_falcON_snapshot")
    import
    integer(c_int), intent(out) :: ierr   
   end subroutine read_falcON_snapshot
 end interface

contains

 ! map types from falcON to splash
 integer function itypemap_falcON(itype)
  integer, intent(in) :: itype

  select case(itype)
  case(1) ! sinks
     itypemap_falcON = 2
  case(2) ! gas
     itypemap_falcON = 1
  case(3:4)
     itypemap_falcON = itype
  case(5:maxtypes)
     itypemap_falcON = itype+1
  case default
     itypemap_falcON = 5 ! unknown
  end select
 
 end function itypemap_falcON

 ! get starting position in particle array
 integer function ioffset(itype,npartoftype)
  integer, intent(in) :: itype
  integer, intent(in) :: npartoftype(:)
  integer :: i
  
  ioffset = 0
  do i=1,size(npartoftype)
     if (i < itype) ioffset = ioffset + npartoftype(i)
  enddo
 
 end function ioffset

end module falcONhdf5read

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------
subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol
  use params,         only:doub_prec,maxparttypes !,maxplot
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread, &
                           ntypes,debugmode,iverbose,buffer_steps_in_file
  use mem_allocation, only:alloc
  use labels,         only:print_types,labeltype
  use system_utils,   only:lenvironment
  use asciiutils,     only:cstring
  use dataread_utils, only:check_range
  use falcONhdf5read
  implicit none
  integer, intent(in)                :: istepstart,ipos
  integer, intent(out)               :: nstepsread
  character(len=*), intent(in)       :: rootname
  character(len=len(rootname)+10)    :: datfile
  integer               :: i,j,ierr,ierror(8),istep,nsteps_to_read
  integer               :: ncolstep,npart_max,nstep_max,ntoti,ntotall
  integer               :: npartoftypei(maxparttypes)
  logical               :: iexist,reallocate,debug,goterrors
  real(doub_prec)       :: timetemp,hperiodic(3)
  !integer, dimension(maxplot) :: isrequired

  nstepsread = 0
  goterrors  = .false.

  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** '
     return
  endif
!
!--check if first data file exists
!
  if (iverbose==1 .and. ipos==1) print "(1x,a)",'reading FalcON hdf5 format'
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     !
     !--append .h5 on the end if not already present
     !
     datfile=trim(rootname)//'.h5'
     inquire(file=datfile,exist=iexist)
     if (.not.iexist) then
        print "(a)",' *** error: '//trim(rootname)//': file not found ***'
        return
     endif
  endif
  !
  ! set parameters which do not vary between timesteps
  !
  ndim  = 3
  ndimV = 3
  debug = (debugmode .or. lenvironment('FSPLASH_DEBUG'))
  ! 
  ! read data from snapshots
  !
  i = istepstart
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  ! open file and read header information
  !
  if (debug) print*,'DEBUG: reading header...'

  call open_falcON_file(cstring(datfile),ierr)
  if (ierr /= 0) then
     print "(a)", '*** ERROR OPENING FALCON FILE ***'
     return
  endif
  
  if (falcON_file_is_open() /= 1) then
     print "(a)", '*** ERROR: falcON_file_is_open /= 1 after opening ***'
     return
  endif

  if (debug) call set_falcON_debugging_level(3);
  
  nstep_max = num_falcON_snapshots(ierr);
  if (debug) print*,'got ',nstep_max,' falcON snapshots in file'
  if (nstep_max <= 0) then
     print "(a)",'*** ERROR: no falcON snapshots found in file ***'
     return
  endif
  
  ntotall = 0
  if (buffer_steps_in_file) then
     nsteps_to_read = nstep_max
  else
     nsteps_to_read = 1
  endif

  i = istepstart
  over_snapshots: do istep=1,nstep_max
  !
  ! read falcON header
  !
  if (debug) print*,'DEBUG: opening snapshot ',i
  npartoftypei(:) = 0
  call open_falcON_snapshot(ntypes, npartoftypei, &
                            ncolstep, ndim, ndimV, timetemp, &
                            hperiodic, ierr)
  !
  ! error checking on header info
  !
  ierror(:) = 0
  call check_range(ntypes,'ntypes',min=1,err=ierror(1))
  call check_range(npartoftypei(1:ntypes),'npartoftype',min=0,err=ierror(2))
  call check_range(sum(npartoftypei(1:ntypes)),'ntot',min=1,err=ierror(3))
  call check_range(ndim,'ndim',min=1,max=3,err=ierror(4))
  call check_range(ndimV,'ndimV',min=ndim,max=3,err=ierror(5))
  call check_range(timetemp,'time',min=0.d0,err=ierror(6))
  call check_range(ierr,'error during header read',min=0,max=0,err=ierror(7))
  if (any(ierror(1:7) > 0)) then
     print*,'*** ERROR during falcON header read ***'
     return
  endif
  ncolumns = ncolstep
  ntoti = sum(npartoftypei(1:ntypes))
  ntotall = max(ntoti,ntotall)
  !
  ! print header information
  !
  if (iverbose >= 1 .and. buffer_steps_in_file .or. istep.eq.ipos) then
     !print "(2(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep
     !print "(a,i2)",' ntypes: ',ntypes,' 
     !print*,' npartoftype = ',(npartoftypei(itypemap_falcON(j)),j=1,ntypes)
     !print*,' ncolstep = ',ncolstep,' ndim = ',ndim,ndimV
     print*,' time = ',timetemp !,' hper = ',hperiodic(:)
  endif
  !
  ! now read data
  !
  reallocate = .false.
  npart_max = maxpart

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*ntotall)
     else
        ! if first time, save on memory
        npart_max = int(ntoti)
     endif
  endif
  !
  ! reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nsteps_to_read,max(ncolumns+ncalc,maxcol))
  endif

  !
  ! copy header data into allocated arrays
  !
  if (buffer_steps_in_file .or. istep.eq.ipos) then
     do j=1,ntypes
        npartoftype(itypemap_falcON(j),i) = npartoftypei(j)
     enddo
     time(i) = real(timetemp)
     masstype(:,i) = 0. ! all masses read from file
  endif
  !
  ! read particle data
  !
  got_particles: if (ntoti > 0) then
     !isrequired(:) = 0
     !where (required(1:ncolumns)) isrequired(1:ncolumns) = 1
     if (buffer_steps_in_file .or. istep.eq.ipos) then
        i_current_step = i
        call read_falcON_snapshot(ierr);
        if (ierr /= 0) then
           print "(/,1x,a,/)",' *** ERROR reading falcON snapshot ***'
           print*,'Press any key to continue (but there is likely something wrong with the file...)'
           read*
        endif
        call print_types(npartoftype(:,i),labeltype)
        i = i + 1
     endif
     nstepsread = nstepsread + 1
  else
     !
     ! cover the special case where no particles have been read
     !
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.

  endif got_particles
  
  enddo over_snapshots
  !
  ! now memory has been allocated, set arrays which are constant for all time
  !
  gamma = 5./3.
  !
  ! set flag to indicate that only part of this file has been read
  !
  ipartialread = .false.
  !if (.not.all(required(1:ncolstep))) ipartialread = .true.
  !
  ! call set labels to identify location of smoothing length
  !
  call set_labels

  !if (nstepsread.gt.0) then
  !   print "(a,i10,a)",' >> read ',sum(npartoftype(:,i)),' particles'
  !endif
  call close_falcON_file()
  return

end subroutine read_data

subroutine read_falcON_data_into_splash(icol,npartoftypei,temparr,itypec) bind(c,name="read_falcON_data_into_splash")
  use, intrinsic :: iso_c_binding, only:c_int,c_double
  use particle_data,  only:dat,iamtype,npartoftype
  use settings_data,  only:debugmode
  use labels,         only:label
  use falcONhdf5read, only:itypemap_falcON,ioffset,i_current_step
  implicit none
  integer(kind=c_int), intent(in) :: icol,npartoftypei,itypec
  real(kind=c_double), intent(in) :: temparr(*)
  integer(kind=c_int) :: icolput
  integer :: istart,iend,nmax,itype,i
  
  itype = itypec + 1 ! convert from c to Fortran indexing

  icolput = icol + 1
  if (debugmode) print "(3(a,i2),a,i8)",'DEBUG: Step ',i_current_step,' column ',icol,&
    ' type ',itypemap_falcON(itype),' -> '//trim(label(icolput))
  
  ! check column is within array limits
  if (icolput.gt.size(dat(1,:,1)) .or. icolput.eq.0) then
     print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in receive_data_fromc'
     return
  endif

  ! ensure no array overflows
  istart = ioffset(itypemap_falcON(itype),npartoftype(:,i_current_step)) + 1
  iend   = min(istart + npartoftypei - 1,size(dat(:,1,1)))
  nmax   = iend - istart + 1

  ! copy data into main splash array
  if (debugmode) print*,'DEBUG: COPYING TO ',istart,iend,' total = ',1,nmax
  
  ! this should never happen
  if (i_current_step < 1 .or. i_current_step > size(dat(1,1,:))) then
     print*,'INTERNAL ERROR in indexing during falcON read'
     return
  endif
  dat(istart:iend,icolput,i_current_step) = real(temparr(1:nmax),kind=kind(dat))

  ! set particle type
  if (size(iamtype(:,1)).gt.1) then
     print*,' SETTING TYPES ',istart,iend
     do i=istart,iend
        iamtype(i,i_current_step) = int(itypemap_falcON(itype),kind=kind(iamtype))
     enddo
  endif

  return
end subroutine read_falcON_data_into_splash

!------------------------------------------------------------
! set labels for each column of data
!------------------------------------------------------------
subroutine set_labels
  use labels,        only:label,iamvec,labelvec,ix,ivx,ipmass, &
                          ih,irho,iax,iutherm !ipr,iutherm
  use settings_data,  only:ndim,ndimV,UseTypeInRenderings
  use geometry,       only:labelcoord
  use falcONhdf5read, only:blocklabel
  use asciiutils,     only:lcase
  implicit none
  integer :: i,icol

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif

  ix = 0
  iutherm = 0
  do icol=1,size(blocklabel)
     select case(trim(lcase(blocklabel(icol))))
     case('x')
        ix(1) = icol
     case('y')
        ix(2) = icol
     case('z')
        ix(3) = icol
     case('vx')
        ivx = icol
     case('ax')
        iax = icol
     case('h')
        ih = icol
     case('mass')
        ipmass = icol
     case('srho')
        irho = icol
     end select
     label(icol) = trim(blocklabel(icol))
  enddo

  ! set labels of the quantities read in
  if (ix(1).gt.0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
  if (irho.gt.0)    label(irho)       = 'density'
  !if (iutherm.gt.0) label(iutherm)    = 'u'
  if (ipmass.gt.0)  label(ipmass)     = 'particle mass'

  ! set labels for vector quantities
  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = trim(labelvec(ivx))//'_'//labelcoord(i,1)
     enddo
  endif

  if (iax.gt.0) then
     iamvec(iax:iax+ndimV-1) = iax
     labelvec(iax:iax+ndimV-1) = 'a'
     do i=1,ndimV
        label(iax+i-1) = trim(labelvec(iax))//'_'//labelcoord(i,1)
     enddo
  endif

  ! labels for each particle type already set
  UseTypeInRenderings(:) = .false.
  UseTypeInRenderings(1) = .true.

!-----------------------------------------------------------
  return
end subroutine set_labels

subroutine set_splash_block_label(icol,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 use falcONhdf5read, only:blocklabel
 use asciiutils,     only:fstring
 implicit none
 integer(kind=c_int),    intent(in) :: icol
 character(kind=c_char), intent(in) :: name(256)

 blocklabel(icol+1) = trim(fstring(name))
 !print*,icol,' name = ',trim(blocklabel(icol))

end subroutine set_splash_block_label

subroutine set_splash_particle_label(itypec,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 use asciiutils,     only:fstring
 use labels,         only:labeltype
 use falcONhdf5read, only:itypemap_falcON
 implicit none
 integer(kind=c_int),    intent(in) :: itypec
 character(kind=c_char), intent(in) :: name(256)

 !print*,' got type = ',itypec,' setting ',itypemap_falcON(itypec+1),'= ',trim(fstring(name))
 labeltype(itypemap_falcON(itypec+1)) = trim(fstring(name))

end subroutine set_splash_particle_label
