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
! THIS VERSION IS FOR JOE'S 2D SPH CODE
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
! ntot(maxstep)       : total number of particles in each timestep
! iam(maxpart,maxstep): integer identification of particle type
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
  use settings_data, only:ndim,ndimV,ncolumns
  use mem_allocation
  implicit none
  integer, intent(IN) :: indexstart,ipos
  integer, intent(OUT) :: nstepsread
  character(LEN=*), intent(IN) :: rootname
  integer :: i,j,ifile,ierr
  integer :: istep,nprint,npart_max,nstep_max,icol
  logical :: iexist
  character(LEN=LEN(rootname)+4) :: dumpfile
  real :: timei,dti,hi

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  ifile = 1

  dumpfile = trim(rootname)
!  if (index(dumpfile,'.plt').eq.0) dumpfile = trim(rootname)//'.plt'
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
  ndim = 2
  ndimV = 2
  ncolumns = 7  ! number of columns in file
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,2)

  j = indexstart
  nstepsread = 0

  print "(1x,a)",'reading Joe Monaghan ascii format'
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the file and read the number of particles
  !
     open(unit=15,iostat=ierr,file=dumpfile,status='old',form='formatted')
     if (ierr /= 0) then
        print*,'*** ERROR OPENING ',trim(dumpfile),' ***'
     else
        !
        !--read the number of particles in the first step,
        !  allocate memory and rewind
        !
        read(15,*,end=55,iostat=ierr) istep,nprint,timei,dti
        print*,'first time = ',timei,nprint
        if (.not.allocated(dat) .or. (nprint.gt.npart_max)) then
           npart_max = max(npart_max,INT(1.1*(nprint)))
           call alloc(npart_max,nstep_max,ncolumns)
        endif
        rewind(15)
     endif
     if (ierr /= 0) then
        print*,'*** ERROR READING TIMESTEP HEADER ***'
     else

        oversteps: do
    !
    !--loop over the timesteps in this file
    !
           npart_max = max(npart_max,nprint)
    !
    !--allocate/reallocate memory if j > maxstep
    !
           if (j.gt.maxstep) then
              call alloc(maxpart,2*j,maxcol)
           endif
    !
    !--now read the timestep data in the dumpfile
    !
           read(15,*,end=55,iostat=ierr) istep,nprint,hi,time(j),dti
           do i=1,nprint
              read(15,*,end=55,iostat=ierr) (dat(i,icol,j),icol = 1,ncolumns)
           enddo
           masstype(1,j) = 1. 

           if (ierr /= 0) then
              print "(a)",'|*** ERROR READING TIMESTEP ***'
              return
           else
              nstepsread = nstepsread + 1
           endif

           npartoftype(:,j) = 0
           npartoftype(1,j) = nprint

           print*,j,' time = ',time(j)
           gamma(j) = 1.666666666667
           j = j + 1

        enddo oversteps
     endif

55 continue
  !
  !--reached end of file
  !
  close(15)

  print*,'nstepsread = ',nstepsread
  print*,'>> end of dump file: nsteps =',j-1,'ntot = ',npartoftype(1,j-1),'nptmass=',npartoftype(2,j-1)

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
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  ivx = ndim+1
  ipr = ndim+ndimV+2
  ih = ndim+ndimV+3        !  smoothing length
  label(ih) = 'h'
  irho = ndim+ndimV+1     ! location of rho in data array
  label(irho) = 'density'
  if (ipr.gt.0) label(ipr) = 'pressure'
  iutherm = 0  !  thermal energy
!  label(iutherm) = 'u'
  ipmass = 0  !  particle mass
!  label(ipmass) = 'particle mass'

  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  !
  !--set labels for each particle type
  !
  ntypes = 1 !!maxparttypes
  labeltype(1) = 'gas'
  UseTypeInRenderings(1) = .true.

!-----------------------------------------------------------

  return
end subroutine set_labels
