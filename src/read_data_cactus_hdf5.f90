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
!  Copyright (C) 2005-2017 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR HDF5 OUTPUT FROM THE CACTUS CODE
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
module cactushdf5read
 use params, only:maxplot,doub_prec
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none
 real :: hsoft
 character(len=lenlabel), dimension(maxplot) :: blocklabel
 logical :: havewarned = .false.
 integer, parameter :: maxtypes = 6

 interface
   subroutine read_cactus_hdf5_header(filename,npart,ncol,ndim,ndimV,time,ierr) bind(c)
    import
    character(kind=c_char), dimension(*), intent(in) :: filename
    integer(kind=c_int), intent(out) :: npart,ncol,ndim,ndimV,ierr
    real(kind=c_double), intent(out) :: time
   end subroutine read_cactus_hdf5_header

   subroutine read_cactus_hdf5_data(filename,maxtypes,npartoftypei,&
                                    ncol,isrequired,ierr) bind(c)
    import
    implicit none
    character(kind=c_char), dimension(*), intent(in)  :: filename
    integer(kind=c_int), intent(in), value :: maxtypes
    integer(kind=c_int), dimension(6), intent(in) :: npartoftypei
    integer(kind=c_int), intent(in), value  :: ncol
    integer(kind=c_int), intent(out) :: ierr
    integer(kind=c_int), dimension(ncol), intent(in)  :: isrequired
   end subroutine read_cactus_hdf5_data
 end interface

end module cactushdf5read

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------
subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep,iamtype
  use params,         only:doub_prec,maxparttypes,maxplot
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread, &
                           ntypes,debugmode,iverbose
  use settings_page,  only:legendtext
  use mem_allocation, only:alloc
  use labels,         only:ih,irho,ipmass,labeltype
  use system_utils,   only:renvironment,lenvironment,ienvironment,envlist
  use asciiutils,     only:cstring
  use cactushdf5read, only:hsoft,blocklabel,havewarned,read_cactus_hdf5_header, &
                           read_cactus_hdf5_data,maxtypes
  use dataread_utils, only:count_types
  implicit none
  integer, intent(in)                :: istepstart,ipos
  integer, intent(out)               :: nstepsread
  character(len=*), intent(in)       :: rootname
  character(len=len(rootname)+10)    :: datfile,densfile,hfile
  character(len=20)                  :: string
  integer               :: i,j,itype,ierr
  integer               :: nunknown
  integer               :: ncolstep,npart_max,nstep_max,ntoti,ntotall,idot
  integer, parameter    :: iunit = 11
  logical               :: iexist,reallocate,debug,goterrors
  integer, dimension(maxplot) :: isrequired
  real(doub_prec) :: timetemp

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
  print "(1x,a)",'reading CACTUS HDF5 format'
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
!--set parameters which do not vary between timesteps
!
  ndim  = 3
  ndimV = 3
! 
!--read data from snapshots
!
  i = istepstart
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open file and read header information
  !
  if (debug) print*,'DEBUG: reading header...'
  call read_cactus_hdf5_header(cstring(datfile),ntoti,ncolstep,ndim,ndimV,timetemp,ierr)
  if (ierr /= 0) then
     print "(a)", '*** ERROR READING HEADER ***'
     return
  endif
  ncolumns = ncolstep

  if (iverbose >= 1) print "(2(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

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
  if (i.ge.maxstep .and. i.ne.1) then
     nstep_max = i + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif
  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)
  endif

  !
  !--copy header data into allocated arrays
  !
  npartoftype(1,i) = ntoti
  time(i) = real(timetemp)
  masstype(:,i) = 0. ! all masses read from file
  !
  !--read particle data
  !
  got_particles: if (ntoti > 0) then
     
     isrequired(:) = 0
     where (required(1:ncolumns)) isrequired(1:ncolumns) = 1
     
     call read_cactus_hdf5_data(cstring(datfile),ntypes,npartoftype(:,i),ncolumns,isrequired,ierr)

     nstepsread = 1
  endif got_particles
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--count particles by type
!
  call count_types(ntoti,iamtype(:,i),npartoftype(:,i),nunknown)
!
!--set flag to indicate that only part of this file has been read
!
  ipartialread = .false.
!
!--call set labels to identify location of smoothing length
!
  call set_labels
!
!--cover the special case where no particles have been read
!
  if (ntoti.le.0) then
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif

  if (nstepsread.gt.0) then
     print "(a,i10,a)",' >> read ',sum(npartoftype(:,istepstart+nstepsread-1)),' particles'
  endif
  return

end subroutine read_data

subroutine read_cactus_hdf5_data_fromc(icol,ntot,np,temparr) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int,c_double
  use particle_data,  only:dat
  use settings_data,  only:debugmode
  use labels,         only:label
  implicit none
  integer(kind=c_int), intent(in) :: icol,ntot,np
  real(kind=c_double), intent(in) :: temparr(np)
  integer(kind=c_int) :: i,icolput
  integer :: nmax,nerr,idi,i1,i2
  logical :: useids

  icolput = icol
  i1 = ntot-np+1
  i2 = ntot
  
  if (debugmode) print "(a,i2,a,i8,a,i8)",'DEBUG: reading column ',icol,' -> '//trim(label(icolput))//' parts ',i1,' to ',i2
  
  ! check column is within array limits
  if (icolput.gt.size(dat(1,:,1)) .or. icolput.eq.0) then
     print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in receive_data_fromc'
     return
  endif

  ! ensure no array overflows
  nmax = min(i2,size(dat(:,1,1)))

  ! copy data into main splash array
  dat(i1:i2,icolput,1) = real(temparr(1:np))

  return
end subroutine read_cactus_hdf5_data_fromc

subroutine read_cactus_itype_fromc(ntot,np,itype) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int
  use particle_data,  only:iamtype
  implicit none
  integer(kind=c_int), intent(in) :: ntot,np
  integer(kind=c_int), intent(in) :: itype(np)
  integer :: i,i1,i2

  i1 = ntot-np+1
  i2 = ntot
  ! set particle type
  if (size(iamtype(:,1)).gt.1) then
     iamtype(i1:i2,1) = itype(1:np)
  endif

  return
end subroutine read_cactus_itype_fromc

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass, &
                          ih,irho,ipr,iutherm,iBfirst,idivB,iax
  use params
  use settings_data,  only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings,iformat
  use geometry,       only:labelcoord
  use system_utils,   only:envlist,ienvironment
  use cactushdf5read, only:hsoft,blocklabel
  use asciiutils,     only:lcase
  implicit none
  integer :: i,j,icol,irank

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif
  blocklabel(1:5) = (/'x','y','z','h','m'/)

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
     case('mass','m')
        ipmass = icol
     case('dens','density')
        irho = icol
     end select
     label(icol) = trim(blocklabel(icol))
  enddo

  ! set labels of the quantities read in
  if (ix(1).gt.0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
  !if (irho.gt.0)    label(irho)       = 'density'
  !if (iutherm.gt.0) label(iutherm)    = 'u'
  !if (ipmass.gt.0)  label(ipmass)     = 'particle mass'
  !if (ih.gt.0)      label(ih)         = 'h'

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

  ! set labels for each particle type
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'ghost'
  UseTypeInRenderings(:) = .true.

!-----------------------------------------------------------
  return
end subroutine set_labels

subroutine set_blocklabel(icol,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 use cactushdf5read, only:blocklabel
 use asciiutils,    only:fstring
 implicit none
 integer(kind=c_int),    intent(in) :: icol
 character(kind=c_char), intent(in) :: name(12)

 blocklabel(icol) = trim(fstring(name))
 !print*,icol,' name = ',trim(blocklabel(icol))

end subroutine set_blocklabel
