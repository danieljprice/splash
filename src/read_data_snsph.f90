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
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM THE SNSPH CODE
! USING THE SELF-DESCRIBING FORMAT
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
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname

  integer :: j,ierr,ntoti,getcol
  integer :: npart_max,nstep_max,ncolstep
  real    :: timei, gammai
  logical :: iexist

  character(len=len(rootname)+10) :: dumpfile

  nstepsread = 0
  npart_max = maxpart

  dumpfile = trim(rootname)
  !
  !--check if first data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: ',trim(dumpfile),' file not found ***'
     return
  endif
  !
  !--fix number of spatial dimensions
  !
  ndim = 3
  ndimV = 3
  !
  !--allocate memory initially
  !
  nstep_max = max(indexstart,1)

  j = indexstart
  nstepsread = 0

  print "(1x,a)",'reading SNSPH format'
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--read number of columns and number of particles
  !
  ncolstep = getcol()
  call getdata(dumpfile,len_trim(dumpfile),npart_max,timei,gammai)

  !
  !--get number of particles from header and allocate memory
  !
  ntoti = npart_max
  ncolumns = ncolstep

  if (.not.allocated(dat) .or. ntoti.gt.npart_max) then
     npart_max = max(npart_max,INT(1.1*ntoti))
     call alloc(npart_max,nstep_max,ncolstep+ncalc)
  endif

  npart_max = max(npart_max,ntoti)
!
!--allocate/reallocate memory if j > maxstep
!
  if (j.gt.maxstep) then
     call alloc(maxpart,j+2*nstepsread,maxcol)
  endif
!
!--now read the timestep data in the dumpfile
!
  time(j) = timei

  print "(a,i5,a,f10.3,a,i8)",'| step ',j,': t = ',time(j),' ntotal = ',ntoti

  nstepsread = nstepsread + 1
!
!--set particle numbers
!
  npartoftype(:,j) = 0
  npartoftype(1,j) = ntoti
!
!--now read data
!
  call readsdf(dumpfile,len_trim(dumpfile),dat(:,:,j),maxpart,maxcol,ierr)

  if (nstepsread .gt. 0 .and. j.gt.0) then
     print*,'>> end of dump file: ntotal = ',sum(npartoftype(:,j))
  endif

  return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use settings_data
  use geometry, only:labelcoord
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
  ivx = 4
  ih = 10        !  smoothing length
  iutherm = 8  !  thermal energy
  ipmass = 7   !  particle mass
  irho = 9    ! location of rho in data array

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = 'density'
  label(iutherm) = 'temperature'
  label(ih) = 'h'
  label(ipmass) = 'particle mass'
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  !
  !--set labels for each particle type
  !
  ntypes = 1
  labeltype(1) = 'gas'

!-----------------------------------------------------------

  return
end subroutine set_labels
