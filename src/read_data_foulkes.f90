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
! THIS VERSION IS FOR STEVE FOULKES' ASCII DATA FORMAT
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
! ntot(maxstep)       : total number of particles in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data, only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,xorigin
  use mem_allocation, only:alloc
  use system_utils, only:lenvironment
  use labels, only:ih
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,nerr,iunit,ncolstep
  integer :: nprint,npart_max,nstep_max,icol
  integer :: nmodel,nstar,ncolread
  logical :: iexist
  real :: tread,hmax,dtmin,tdg,hfac
  real, dimension(3) :: xptmass,yptmass,vxptmass,vyptmass
  character(len=len(rootname)+4) :: dumpfile

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  iunit = 15  ! logical unit number for input

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
  !--fix number of spatial dimensions (0 means no particle coords)
  !
  ndim = 3
  ndimV = 3
  nstar = 2

  j = indexstart
  nstepsread = 0
  print "(a)",' Steve Foulkes/Carol Haswell/James Murray ascii data format'
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

  !
  !--open the file and read the number of particles
  !
  open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='formatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  endif
!
!--read header line, set time
!
  read(iunit,*,iostat=ierr) nmodel,nprint,hmax,tread,dtmin,tdg
  if (ierr /= 0) print "(a)",' WARNING: error(s) reading first header line'
  print "(a,1pe10.3,a,i10,a,i10)",' time = ',tread,' npart =',nprint,' model number = ',nmodel
  print "(3(a,1pe10.4))",' hmax = ',hmax,' dtmin = ',dtmin,' tdg = ',tdg
!
!--determine how many columns we are going to read
!
  if (lenvironment('FSPLASH_READALL')) then
     ncolstep = 26
     print "(a)",' reading all columns'
  else
     ncolstep = 15
     print "(a,i2)",' reading only to column 15 (setenv FSPLASH_READALL ''yes'' to read all)'
  endif
!
!--(re)allocate memory
!
  nstep_max = max(maxstep,nstep_max,indexstart,1)
  if (.not.allocated(dat) .or. (nprint.gt.maxpart) .or. (ncolstep+ncalc).gt.maxcol) then
     npart_max = max(npart_max,nprint+nstar,maxpart)
     !--allow extra room if reallocating
     if (allocated(dat)) npart_max = max(npart_max,INT(1.1*(nprint+nstar)),maxpart)
     call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol))
  endif
!
!--set the necessary parameters
!
  ncolumns = ncolstep
  nstepsread = nstepsread + 1
  npartoftype(:,j) = 0
  npartoftype(1,j) = nprint
  npartoftype(2,j) = nstar
  time(j) = tread
  gamma(j) = 5./3.
!
!--now read the timestep data in the dumpfile
!
  nerr = 0
!
!--only read required columns
!
  ncolread = 0
  do i=1,ncolstep
     if (required(i)) ncolread = i
  enddo
  if (ncolread.ne.ncolstep) then
     ipartialread = .true.
     print*,' reading only up to column ',ncolread
  endif

  do i=1,nprint
     read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,ncolread)
     if (ierr.ne.0) nerr = nerr + 1
  enddo
  if (nerr > 0) print *,' ERRORS reading particle data on ',nerr,' lines'
!
!--read point mass information
!
  read(iunit,*,iostat=ierr) xptmass(1),yptmass(1),xptmass(2),yptmass(2),hfac
!--point mass velocities not read, though would put them here if they were
  vxptmass(1) = 0.
  vyptmass(1) = 0.
  vxptmass(2) = 0.
  vyptmass(2) = 0.
  if (ierr /= 0) print *,' ERROR reading primary and secondary positions'

  close(iunit)

!--set labels to get ih for setting smoothing length of stars
  call set_labels

!--copy star particle properties into main data array
  do i=nprint+1,nprint+nstar
     dat(i,1,j) = xptmass(i-nprint)
     dat(i,2,j) = yptmass(i-nprint)
     dat(i,3,j) = 0.
     dat(i,4,j) = vxptmass(i-nprint)
     dat(i,5,j) = vyptmass(i-nprint)
     dat(i,6,j) = 0.
     dat(i,7:ncolstep,j) = 0.
     if (ih.gt.0) dat(i,ih,j) = epsilon(0.) ! small but non-zero smoothing length
  enddo
  print "(' primary     (x,y) = (',1pe10.2,',',1pe10.2,')')",xptmass(1),yptmass(1)
  print "(' secondary   (x,y) = (',1pe10.2,',',1pe10.2,')')",xptmass(2),yptmass(2)
  print "(a)",' setting origin to primary position... '
  xorigin(1) = xptmass(1)
  xorigin(2) = yptmass(1)
  xorigin(3) = 0.


return
end subroutine read_data

!!-------------------------------------------------------------------
!! set labels for each column of data
!!
!! read these from a file called 'columns' in the current directory
!! then take sensible guesses as to which quantities are which
!! from the column labels
!!
!!-------------------------------------------------------------------

subroutine set_labels
  use labels, only:label,labeltype,ix,irho,ipmass,ih,ivx,iamvec,labelvec
  use params
  use settings_data, only:ntypes,ndim,ndimV,UseTypeInRenderings
  use geometry, only:labelcoord
  implicit none
  integer :: i,ifx

  do i=1,ndim
     ix(i) = i
  enddo
  ivx = ndim+1
  ipmass = ivx+ndimV
  ifx = ipmass+1
  irho = ifx+ndimV
  label(irho+1) = 'du/dt'
  label(irho+2) = 'C\d\s\u'
  label(irho+3) = 'alpha'
  ih = irho+4
  label(irho+5) = 'kpc'
  label(irho+6) = 'schb'
  label(irho+7) = 'dtp'
  label(irho+8) = 'sigma'
  label(irho+9) = 'pdr2'
  label(irho+10) = 'av_sep'
  label(irho+11) = 'radius'
  label(irho+12) = 'viscosity flag'
  label(irho+13) = 'neighbour number'
  label(irho+14) = 'iwas'
  label(irho+15) = 'll'
  if (irho.gt.0) label(irho) = 'density'
  if (ih.gt.0) label(ih) = 'smoothing length'
  if (ipmass.gt.0) label(ipmass) = 'particle mass'

  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
       label(ivx+i-1) = 'v\d'//labelcoord(i,1)
     enddo
  endif
  if (ifx.gt.0) then
     iamvec(ifx:ifx+ndimV-1) = ifx
     labelvec(ifx:ifx+ndimV-1) = 'f'
     do i=1,ndimV
       label(ifx+i-1) = 'f\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'star'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.

!-----------------------------------------------------------

  return
end subroutine set_labels
