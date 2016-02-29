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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM ANDREA URBAN'S CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,maxpart,maxcol,maxstep,time,gamma
  use params
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation, only:alloc
  !use system_utils,   only:lenvironment,ienvironment
  use asciiutils,     only:get_ncolumns
  implicit none
  integer, intent(in)          :: indexstart,ipos
  integer, intent(out)         :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,lab,icol,ilen
  integer :: nprint,npart_max,nstep_max
  integer :: ncol,nerr,nheaderlines
  logical :: iexist,timeset,gammaset
  real    :: dummyreal
  character(len=len(rootname)+5) :: dumpfile
  integer, parameter :: iunit = 15

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  dumpfile = trim(rootname)
  !
  !--check if first data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
     return
  endif
  !
  !--fix number of spatial dimensions
  !
  ndim = 3
  ndimV = 3

  !--number of columns to read from file
  !  this is determined automatically
  !ncol = 13

  j = indexstart
  nstepsread = 0
  print "(a)",' reading Andrea Urban ascii file format'

  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the file and read the number of particles
  !
  open(unit=iunit,file=dumpfile,status='old',form='formatted',iostat=ierr)
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  else
     call get_ncolumns(iunit,ncol,nheaderlines)
     ncol = max(ncol - 1,0)
     if (ncol.le.0) then
        print "(a,/)",' *** no data read from file ***'
        return
     endif
     !
     !--allocate memory initially
     !
     nprint = 10001
     nstep_max = max(nstep_max,indexstart,1)
     if (.not.allocated(dat) .or. (nprint.gt.npart_max) .or. (ncol+ncalc).gt.maxcol) then
        npart_max = max(npart_max,nprint)
        call alloc(npart_max,nstep_max,ncol+ncalc)
     endif
  endif

  npart_max = max(npart_max,nprint)
  ncolumns = ncol
!
!--allocate/reallocate memory if j > maxstep
!
  if (j.gt.maxstep) then
     call alloc(maxpart,j+1,maxcol)
  endif

!
!--read header lines, try to use it to set time
!
  timeset = .false.
  gammaset = .false.

  if (nheaderlines.gt.0) then
     print*,'skipping ',nheaderlines,' header lines'
     do i=1,nheaderlines
        read(iunit,*,iostat=ierr) dummyreal
        if (timeset .and. .not.gammaset .and. ierr.eq.0 &
           .and. dummyreal.gt.0.999999 .and. dummyreal.lt.2.000001) then
           print*,'setting gamma = ',dummyreal,' from header line ',i
           gamma(j) = dummyreal
           gammaset = .true.
        endif
        if (ierr.eq.0 .and. .not. timeset) then
           time(j) = dummyreal
           timeset = .true.
           print*,'setting time = ',dummyreal,' from header line ',i
        endif
     enddo
  endif
!
!--now read the timestep data in the dumpfile
!
  dat(:,:,j) = 0.
!  time(j) = -1.0 ! time not read

!
!--now read the timestep data in the dumpfile
!
  i = 0
  ierr = 0
  nerr = 0
  overparts: do while (ierr >= 0)
     i = i + 1
     if (i.gt.npart_max) then ! reallocate memory if necessary
        npart_max = 10*npart_max
        call alloc(npart_max,nstep_max,ncol+ncalc)
     endif
     read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,10),lab,(dat(i,icol,j),icol=11,ncol)
     if (ierr > 0) then
        nerr = nerr + 1
        if (nerr .le. 10) print "(a,i8,a)",' ERROR reading data from line ',i+nheaderlines,', skipping'
        i = i - 1 ! ignore lines with errors
     endif
  enddo overparts
  close(iunit)

  nprint = i - 1
  nstepsread = nstepsread + 1

  if (nerr > 10) then
     print "(a,i8,a)",' *** WARNING: errors whilst reading file on ',nerr,' lines: skipped these ***'
  endif
  if (ierr < 0) then
     print*,'end of file: npart = ',nprint
  endif


  npartoftype(:,j) = 0
  npartoftype(1,j) = nprint

!
!--now open the sink particle file and read it
!
  !--find the last underscore in the file name
  ilen = index(rootname,'_',back=.true.)
  if (ilen.le.0) ilen = len_trim(rootname) + 1
  dumpfile = rootname(1:ilen-1)//'_S'
  inquire(file=trim(dumpfile),exist=iexist)
  if (iexist) then
     open(unit=iunit+1,file=trim(dumpfile),form='formatted',status='old',iostat=ierr)
     if (ierr.ne.0) then
        print "(a)",' ERROR: could not open sink particle file '//trim(dumpfile)
     else
        i = npartoftype(1,j)
        ierr = 0
        nerr = 0
        oversinks: do while (ierr >= 0)
           i = i + 1
           if (i.gt.npart_max) then ! reallocate memory if necessary
              npart_max = npart_max + 1000
              call alloc(npart_max,nstep_max,ncol+ncalc)
           endif
           read(iunit+1,*,iostat=ierr) (dat(i,icol,j),icol = 1,10)
           if (ierr > 0) then
              nerr = nerr + 1
              if (nerr .le. 10) print "(a,i8,a)",' ERROR reading sink data from line ',i+nheaderlines,', skipping'
              i = i - 1 ! ignore lines with errors
           endif
        enddo oversinks
     endif
     npartoftype(2,j) = i - 1 - npartoftype(1,j)
     print "(a,i8,a)",' read ',npartoftype(2,j),' sink particles from '//trim(dumpfile)
     close(iunit+1)
  else
     print "(a)",' sink particle file ('//trim(dumpfile)//') not present'
  endif
!
!--look for a _t file for the time (interim measure)
!
  dumpfile = rootname(1:ilen-1)//'_t'
  inquire(file=trim(dumpfile),exist=iexist)
  if (iexist) then
     open(unit=iunit+2,file=trim(dumpfile),form='formatted',status='old',iostat=ierr)
     if (ierr.ne.0) then
        print "(a)",' ERROR: could not open time file '//trim(dumpfile)
     else
        read(iunit+2,*,iostat=ierr) time(j)
        if (ierr.ne.0) then
           print "(a)",' ERROR reading time from file '//trim(dumpfile)
        else
           print*,' got time = ',time(j),' from file '//trim(dumpfile)
        endif
     endif
     close(iunit+2)
  endif

return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels,        only:label,labelvec,labeltype,iamvec,&
                           ix,ivx,ih,irho,iutherm,ipmass
  use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings
  use geometry,      only:labelcoord
  !use settings_units, only:units,unitslabel
  implicit none
  integer :: i

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif

  do i=1,ndim
     ix(i) = i
  enddo
  ivx     = 4
  ipmass  = 7
  ih      = 8
  irho    = 9
  iutherm = 10

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(ipmass) = 'particle mass'
  label(ih)     = 'h'
  label(irho)   = 'density'
  if (iutherm.gt.0) label(iutherm) = 'u'
  label(11)     = 't\ddust\u'
  label(12)     = 'N\dcol\u'
  label(13)     = 'N\dloc\u'

  if (ivx.ne.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'sink'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.

!-----------------------------------------------------------

  return
end subroutine set_labels
