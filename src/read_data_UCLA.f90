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
! THIS VERSION IS FOR SKY KING / MARK MORRIS' (UCLA) ASCII DATA FORMAT
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
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation, only:alloc
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,nerr,iunit,ncolstep,ncolenv
  integer :: nprint,npart_max,nstep_max,icol
  integer :: igatherb,ntot,ninit,ninit1,nstar
  logical :: iexist
  real :: tread,pmass,hbav1,dttot,xmacc,xlxacc,xlyacc,xlzacc
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
  ncolstep = 10  ! create one column for particle mass
  nstar = 3

  j = indexstart
  nstepsread = 0
  print "(a)",' reading Sky King/Mark Morris (UCLA) ascii data format '

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
!--read header lines, try to use it to set time
!
  read(iunit,*,iostat=ierr) nprint,tread,pmass,xptmass(1),yptmass(1),xptmass(2),yptmass(2),&
                vxptmass(1),vyptmass(1),vxptmass(2),vyptmass(2), &
                xptmass(3),yptmass(3),vxptmass(3),vyptmass(3)
  if (ierr /= 0) print "(a)",' WARNING: error(s) reading first header line'
  print "(a,i10,a,1pe10.3)",' npart = ',nprint, ' time = ',tread

  read(iunit,*,iostat=ierr) igatherb,ntot,ninit,ninit1,hbav1,dttot,xmacc,xlxacc,xlyacc,xlzacc
  if (ierr /= 0) print "(a)",' WARNING: error(s) reading second header line'
  print "(a)",' header info: '
  print "(2(a11,i10))",' igatherb: ',igatherb,' ntot: ',ntot
  print "(2(a11,i10))",' ninit: ',ninit,' ninit1: ',ninit1
  print "(3(a11,1pe10.4))",' hbav1: ',hbav1,' dttot: ',dttot,' xmacc: ',xmacc
  print "(3(a11,1pe10.4))",' accelx: ',xlxacc,' accely: ',xlyacc,' accelz: ',xlzacc
!
!--(re)allocate memory
!
  nstep_max = max(nstep_max,indexstart,1)
  if (.not.allocated(dat) .or. (nprint.gt.maxpart) .or. (ncolstep+ncalc).gt.maxcol) then
     npart_max = max(npart_max,INT(1.1*(nprint)),maxpart)
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
!
!--now read the timestep data in the dumpfile
!
  nerr = 0
  do i=1,nprint
     read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,9)
     if (ierr.ne.0) nerr = nerr + 1
  enddo
  if (nerr > 0) print *,' ERRORS reading particle data on ',nerr,' lines'
  close(iunit)

!--set particle mass from column
  dat(1:nprint,10,j) = pmass
!--copy star particle properties into main data array
  do i=nprint+1,nprint+nstar
     dat(i,1,j) = xptmass(i-nprint)
     dat(i,2,j) = yptmass(i-nprint)
     dat(i,3,j) = 0.
     dat(i,4,j) = vxptmass(i-nprint)
     dat(i,5,j) = vyptmass(i-nprint)
     dat(i,6,j) = 0.
  enddo


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
  use labels, only:label,labeltype,ix,irho,ipmass,ih,ipr,ivx,iamvec,labelvec
  use params
  use settings_data, only:ncolumns,ntypes,ndim,ndimV,UseTypeInRenderings
  use geometry, only:labelcoord
  implicit none
  integer :: i,ierr,ndimVtemp

  do i=1,ndim
     ix(i) = i
  enddo
  ivx = 4
  irho = 7
  ipr = 8
  ih = 9
  ipmass = 10
  label(irho) = 'density'
  label(ipr) = 'pressure'
  label(ih) = 'smoothing length'
  label(ipmass) = 'particle mass'

  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
       label(ivx+i-1) = 'v\d'//labelcoord(i,1)
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
