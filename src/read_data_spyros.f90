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
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM THE VINE CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! *** CONVERTS TO SINGLE PRECISION ***
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
  use particle_data
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname

  integer :: j,ierr,ntoti,irec
  integer :: npart_max,nstep_max,ncolstep
  logical :: iexist

  character(len=len(rootname)+10) :: dumpfile

  !--we are assuming dump is double precision
  real(doub_prec) :: variable1,variable2
  real(doub_prec), dimension(:), allocatable :: dattemp

  nstepsread = 0
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
  !
  !--allocate memory initially
  !
  nstep_max = max(indexstart,1)

  j = indexstart
  nstepsread = 0

  print "(1x,a)",'reading Spyros Kitsonias'' format'
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted',&
       access='direct',recl=4000000)
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
  else
     !
     !--read timestep header (integers only)
     !
     read(15,iostat=ierr) variable1,variable2
     if (ierr < 0) then
        print "(a)",'*** END OF FILE IN TIMESTEP HEADER ***'
        return
     elseif (ierr /= 0) then
        print "(a)",'*** ERROR READING TIMESTEP HEADER ***'
        return
     endif
     !
     !--get number of particles from header and allocate memory
     !
     ntoti = int(variable1)
     ncolstep = int(variable2)
     ncolumns = ncolstep

     if (.not.allocated(dat) .or. ntoti.gt.npart_max) then
        npart_max = max(npart_max,INT(1.1*ntoti))
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
     !
     !--rewind file
     !
     rewind(15)
  endif

  npart_max = max(npart_max,ntoti)
!
!--allocate/reallocate memory if j > maxstep
!
  if (j.gt.maxstep) then
     call alloc(maxpart,j+2*nstepsread,maxcol)
  endif
!
!--allocate a temporary array for double precision variables
!
  if (allocated(dattemp)) deallocate(dattemp)
  allocate(dattemp(npart_max),stat=ierr)
  dattemp = 0.
  if (ierr /= 0) print*,'not enough memory in read_data'
!
!--now read the timestep data in the dumpfile
!
  print "(a,i5,a,f8.3,a,i8)",'| step ',j,': t = ',time(j),' ntotal = ',ntoti

  nstepsread = nstepsread + 1
!
!--set particle numbers
!
  npartoftype(:,j) = 0
  npartoftype(1,j) = ntoti
!
!--read all records
!
  overcolumns: do irec=1,ncolstep
     read(15,rec=irec,iostat=ierr) variable1,variable2,dattemp(1:ntoti)
     if (ierr < 0) then
        print "(a)",'*** END OF FILE IN READ DATA (CHECK PRECISION) ***'
        exit overcolumns
     elseif (ierr /= 0) then
        print "(a)",'*** ERROR READING DATA ***'
        exit overcolumns
     endif
!
!--time and gamma
!
     if (irec.eq.2) time(j) = real(variable1)
     if (irec.eq.5) gamma(j) = real(variable2)
!
!--convert to single precision
!
     dat(1:ntoti,irec,j) = real(dattemp(1:ntoti))

  enddo overcolumns
!
!--clean up
!
  if (allocated(dattemp)) deallocate(dattemp)
!
!--close file
!
 close(15)
 if (nstepsread .gt. 0) then
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
  iutherm = 9  !  thermal energy
  ipmass = 8   !  particle mass
  irho = 7    ! location of rho in data array

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
  UseTypeInRenderings(1) = .true.

!-----------------------------------------------------------

  return
end subroutine set_labels
