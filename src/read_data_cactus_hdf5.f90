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
 character(len=lenlabel), dimension(maxplot) :: blocklabel
 character(len=130) :: datfileprev = ' '
 logical :: file_is_open = .false.
 integer :: ntoti_prev,ncol_prev,nstep_prev

 interface
   subroutine open_cactus_hdf5_file(filename,istep,npart,ncol,nstep_max,ndim,ndimV,time,ierr) bind(c)
    import
    character(kind=c_char), dimension(*), intent(in) :: filename
    integer(kind=c_int), intent(in), value :: istep
    integer(kind=c_int), intent(out) :: npart,ncol,nstep_max,ndim,ndimV,ierr
    real(kind=c_double), intent(out) :: time
   end subroutine open_cactus_hdf5_file

   subroutine read_cactus_hdf5_data(filename,istep,npart,time,dx,ierr) bind(c)
    import
    character(kind=c_char), dimension(*), intent(in)  :: filename
    integer(kind=c_int), intent(in), value :: istep
    integer(kind=c_int), intent(out) :: npart,ierr
    real(kind=c_double), intent(out) :: time,dx
   end subroutine read_cactus_hdf5_data

   subroutine close_cactus_hdf5_file(ierr) bind(c)
    import
    integer(kind=c_int), intent(out) :: ierr
   end subroutine close_cactus_hdf5_file

 end interface

end module cactushdf5read

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------
subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep,iamtype
  use params,         only:doub_prec
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread,iverbose,buffer_steps_in_file
  use settings_page,  only:legendtext
  use mem_allocation, only:alloc
  use labels,         only:ih,irho,ipmass
  use system_utils,   only:renvironment,lenvironment,ienvironment,envlist
  use asciiutils,     only:cstring
  use cactushdf5read, only:open_cactus_hdf5_file,read_cactus_hdf5_data,close_cactus_hdf5_file,&
                           datfileprev,file_is_open,ntoti_prev,ncol_prev,nstep_prev
  use dataread_utils, only:count_types
  implicit none
  integer, intent(in)                :: istepstart,ipos
  integer, intent(out)               :: nstepsread
  character(len=*), intent(in)       :: rootname
  character(len=len(rootname)+10)    :: datfile
  integer               :: i,istep,ierr
  integer               :: nunknown
  integer               :: ncolstep,npart_max,nstep_max,nsteps_to_read,ntoti
  logical               :: iexist,reallocate,goterrors
  real(doub_prec) :: timetemp,dx,vol

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
  if (iverbose==1 .and. ipos==1) print "(1x,a)",'reading CACTUS HDF5 format'
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     !
     !--append .h5 on the end if not already present
     !
     datfile=trim(rootname)//'.h5'
     inquire(file=datfile,exist=iexist)
  endif
!
!--close previous file if filenames do not match
!
  if (trim(datfile)/=trim(datfileprev) .and. file_is_open) then
     call close_cactus_hdf5_file(ierr)
     file_is_open = .false.
  endif
  
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(rootname)//': file not found ***'
     return
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
  if (.not.file_is_open) write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open file and read header information
  !
  if (file_is_open) then
     ntoti     = ntoti_prev
     ncolstep  = ncol_prev
     nstep_max = nstep_prev
  else
     call open_cactus_hdf5_file(cstring(datfile),ipos,ntoti,ncolstep,nstep_max,ndim,ndimV,timetemp,ierr)
     if (ierr /= 0) then
        print "(a)", '*** ERROR READING HEADER ***'
        call close_cactus_hdf5_file(ierr)
        return
     endif
     file_is_open = .true.
     datfileprev = datfile
     ntoti_prev = ntoti
     ncol_prev  = ncolstep
     nstep_prev = nstep_max
  endif
  ncolumns = ncolstep
  if (iverbose >= 1) print "(3(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep,' nsteps: ',nstep_max

  istep = 1 
  over_snapshots: do istep=1,nstep_max
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  if (buffer_steps_in_file) then
     nsteps_to_read = nstep_max
  else
     nsteps_to_read = max(maxstep,1)
  endif
  if (nsteps_to_read > maxstep) reallocate = .true.

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*ntoti)
     else
        ! if first time, save on memory
        npart_max = int(ntoti)
     endif
  endif
  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nsteps_to_read,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)
  endif

  !
  !--read particle data
  !
  got_particles: if (ntoti > 0) then

     if (buffer_steps_in_file .or. ipos.eq.istep) then
        call read_cactus_hdf5_data(cstring(datfile),istep,ntoti,timetemp,dx,ierr)
        call set_labels
        ! set smoothing length and particle mass
        if (ih > 0) dat(:,ih,i) = real(dx)
        vol = dx**ndim
        if (ipmass > 0 .and. irho > 0) dat(:,ipmass,i) = dat(:,irho,i)*real(vol)
        call count_types(ntoti,iamtype(:,i),npartoftype(:,i),nunknown)
        masstype(:,i) = 0. ! all masses read from file
        time(i) = real(timetemp)
        i = i + 1
     endif

     nstepsread = nstepsread + 1

  endif got_particles
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read
!
  ipartialread = .false.
!
!--call set labels to identify location of smoothing length
!
  call set_labels
  
  enddo over_snapshots

  if (nstepsread.gt.0) then
     print "(a,i10,a)",' >> read ',sum(npartoftype(:,istepstart)),' cells'
  endif

  !if (ipos==nstep_max) then
  !   call close_cactus_hdf5_file(ierr)
  !   file_is_open = .false.
  !endif

end subroutine read_data

subroutine read_cactus_hdf5_data_fromc(icol,ntot,np,temparr) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int,c_double
  use particle_data,  only:dat
  use settings_data,  only:debugmode
  use labels,         only:label
  implicit none
  integer(kind=c_int), intent(in) :: icol,ntot,np
  real(kind=c_double), intent(in) :: temparr(np)
  integer(kind=c_int) :: icolput
  integer :: nmax,i1,i2

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
  use params,         only:int1
  implicit none
  integer(kind=c_int), intent(in) :: ntot,np
  integer(kind=c_int), intent(in) :: itype(np)
  integer :: i1,i2

  i1 = ntot-np+1
  i2 = ntot
  ! set particle type
  if (size(iamtype(:,1)).gt.1) then
     iamtype(i1:i2,1) = int(itype(1:np),kind=int1)
  endif

  return
end subroutine read_cactus_itype_fromc

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,iutherm,ih,irho
  use params
  use settings_data,  only:ndim,ndimV,ntypes,UseTypeInRenderings
  use geometry,       only:labelcoord
  use system_utils,   only:envlist,ienvironment
  use cactushdf5read, only:blocklabel
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
  blocklabel(1:5) = (/'x','y','z','h','m'/)

  ix = 0
  iutherm = 0
  irho = 0
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
     case('h')
        ih = icol
     case('mass','m')
        ipmass = icol
     case('dens','density','rho')
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

  ! set labels for each particle type
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'ghost'
  UseTypeInRenderings(:) = .true.
  UseTypeInRenderings(2) = .false.
  

!-----------------------------------------------------------
  return
end subroutine set_labels

subroutine set_blocklabel(icol,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 use cactushdf5read, only:blocklabel
 use asciiutils,    only:fstring
 implicit none
 integer(kind=c_int),    intent(in) :: icol
 character(kind=c_char), intent(in) :: name(24)
 character(len=24) :: temp
 integer :: ivar
 
 temp = fstring(name)
 ivar = index(temp,'::')
 if (ivar > 0) temp = temp(ivar+2:)
 blocklabel(icol) = trim(temp)
 !print*,icol,' name = ',trim(blocklabel(icol))

end subroutine set_blocklabel

subroutine sort_cactus_data(n,iter,iorder) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int
 use sort, only:indexxi
 implicit none
 integer(kind=c_int), intent(in)  :: n
 integer(kind=c_int), intent(in)  :: iter(n)
 integer(kind=c_int), intent(out) :: iorder(n)

 call indexxi(n,iter,iorder)

end subroutine sort_cactus_data
