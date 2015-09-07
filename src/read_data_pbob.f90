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
! THIS VERSION IS FOR .PBOB FILES BY DAVID BROWN
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
!-------------------------------------------------------------------------
!
!  The module below contains interface routines to c functions
!  that perform the actual read of the .pbob file information
!
!-------------------------------------------------------------------------
module pbobread
 use params, only:maxplot,doub_prec
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none
 character(len=lenlabel), dimension(maxplot) :: blocklabel
 integer, parameter :: maxtypes = 6

 interface
   subroutine read_pbob_header(filename,npart,ncol,ndim,ndimV,time,ierr) bind(c)
    import
    character(kind=c_char), dimension(*), intent(in) :: filename
    integer(kind=c_int), intent(out) :: npart,ncol,ndim,ndimV,ierr
    real(kind=c_double), intent(out) :: time
   end subroutine read_pbob_header

   subroutine read_pbob_data(filename,np,time_slice,ierr) bind(c)
    import
    implicit none
    character(kind=c_char), dimension(*), intent(in)  :: filename
    integer(kind=c_int), intent(in), value :: np
    integer(kind=c_int), intent(in), value :: time_slice
    integer(kind=c_int), intent(out) :: ierr
   end subroutine read_pbob_data
 end interface

end module pbobread

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------
subroutine read_data(rootname,istepstart,nstepsread)
  use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep,iamtype
  use params,         only:doub_prec,maxparttypes,maxplot
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread, &
                           ntypes,debugmode,iverbose
  use mem_allocation, only:alloc
  use labels,         only:labeltype,print_types
  use asciiutils,     only:cstring
  use pbobread,       only:blocklabel,read_pbob_header,read_pbob_data,maxtypes
  implicit none
  integer, intent(in)                :: istepstart
  integer, intent(out)               :: nstepsread
  character(len=*), intent(in)       :: rootname
  character(len=len(rootname)+10)    :: datfile,densfile,hfile
  character(len=20)                  :: string
  integer               :: i,j,itype,ierr
  integer               :: index1,index2,nhfac
  integer               :: ncolstep,npart_max,nstep_max,ntoti,ntotall,idot
  integer, parameter    :: iunit = 11
  logical               :: iexist,reallocate,usez,debug,goterrors
  real(doub_prec)       :: timetemp,ztemp
  real :: hfact,hfactmean,pmassi
  real, parameter :: pi = 3.1415926536

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
  print "(1x,a)",'reading PBOB format'
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     !
     !--append .silo on the end if not already present
     !
     datfile=trim(rootname)//'.pbob'
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
  call read_pbob_header(cstring(datfile),ntoti,ncolstep,ndim,ndimV,timetemp,ierr)
  if (ierr /= 0) then
     print "(a)", '*** ERROR READING HEADER ***'
     return
  endif
  ncolumns = ncolstep
  !if (iverbose >= 1) print "(a,1x,i10,a,es10.3)",'  ndim: ',ndim,'     time: ',timetemp
  !if (iverbose >= 1) print "(2(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep
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
     
     !do while(ierr == 0)
        if (i.ge.maxstep .and. i.ne.1) then
           nstep_max = i + max(10,INT(0.1*nstep_max))
           reallocate = .true.
        endif
        !
        !--reallocate memory for main data array
        !
        if (reallocate) call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)

        call read_pbob_data(cstring(datfile),ntoti,i,ierr)
        if (ierr == 0) then
           call set_labels ! sets ntypes and labeltype
           if (size(iamtype(:,i)).gt.1) then
              npartoftype(:,i) = 0
              do j=1,ntoti
                 itype = iamtype(j,i)
                 if (itype > 0 .and. itype <= ntypes) then
                    npartoftype(itype,i) = npartoftype(itype,i) + 1
                 else ! catch-all "unknown" type
                    npartoftype(ntypes+1,i) = npartoftype(ntypes+1,i) + 1
                 endif
              enddo
           endif
           !i = i + 1
           nstepsread = nstepsread + 1
           call print_types(npartoftype(:,i),labeltype)
        endif
        
     !enddo
  endif got_particles
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
  ipartialread = .false.
!
!--cover the special case where no particles have been read
!
  if (ntoti.le.0) then
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif
 
  return
end subroutine read_data

subroutine read_pbob_data_fromc(icol,np,temparr,itype,tag) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
  use particle_data,  only:dat,iamtype
  use settings_data,  only:debugmode
  use labels,         only:label
  use asciiutils,     only:fstring
  use pbobread,       only:blocklabel
  integer(kind=c_int), intent(in),value :: icol,np
  integer(kind=c_int), intent(in) :: itype(np)
  real(kind=c_double), intent(in) :: temparr(np)
  character(kind=c_char), intent(in) :: tag(256)
  integer(kind=c_int) :: i,icolput
  integer :: nmax,nerr,idi

  icolput = icol
  if (debugmode) print "(a,i2,a,i2,a,i8)",'DEBUG: reading column ',icol,' type ',itype,' -> '//trim(label(icolput))

  ! check column is within array limits
  if (icolput.gt.size(dat(1,:,1)) .or. icolput.eq.0) then
     print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in read_data_fromc'
     return
  endif
  
  blocklabel(icol) = trim(fstring(tag))

  ! ensure no array overflows
  nmax = min(np,size(dat(:,1,1)))

  ! copy data into main splash array
  dat(1:nmax,icolput,1) = real(temparr(1:nmax))

  ! set particle type
  if (size(iamtype(:,1)).gt.1) then
     do i=1,nmax
        iamtype(i,1) = int(itype(i),kind=kind(iamtype))
     enddo
  endif

  return
end subroutine read_pbob_data_fromc

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
  use pbobread,       only:blocklabel
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
     case('p')
        ipr = icol
     case('mass','m')
        ipmass = icol
     case('density','rho')
        irho = icol
     end select
     label(icol) = trim(blocklabel(icol))
  enddo

  ! set labels of the quantities read in
  if (ix(1).gt.0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)

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
  ntypes = 3
  labeltype(1) = 'interior'
  labeltype(2) = 'ghost'
  labeltype(3) = 'boundary'
  UseTypeInRenderings(:) = .false.
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .true.
  UseTypeInRenderings(3) = .false.

!-----------------------------------------------------------
  return
end subroutine set_labels
