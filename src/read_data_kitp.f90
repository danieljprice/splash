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
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT IN KITP FORMAT
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
  use particle_data, only:dat,time,npartoftype,gamma,maxpart,maxcol
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation, only:alloc
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr
  integer :: np
  integer :: ncol,nread,nstep_max
  logical :: iexist
  character(len=len(rootname)) :: dumpfile
  real :: wp,timei

  nstepsread = 0
  nstep_max = 0
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
  ncol = 7
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,1)
  j = indexstart
  nstepsread = 0

  print "(1x,a)",'reading KITP SPH format'
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,file=dumpfile,status='old',form='unformatted',iostat=ierr)
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  else
     !
     !--read the number of particles in the header and allocate memory
     !
     read(15,iostat=ierr) np,wp
     timei = 0.
     print "(a,f10.2,a,i10,a,f10.4)",' time: ',timei,' npart: ',np,' wp: ',wp
     !--barf if stupid values read
     if (np.le.0 .or. np.gt.1e10) then
        print "(a)",' *** ERRORS IN TIMESTEP HEADER: WRONG ENDIAN? ***'
        close(15)
        return
     elseif (ierr /= 0) then
        print "(a)",'*** WARNING: ERRORS READING HEADER ***'
        close(15)
        return
     endif
     ncolumns = ncol

     if (.not.allocated(dat) .or. np.gt.maxpart) then
        call alloc(np,nstep_max,ncol+ncalc)
     endif
     !
     !--now read the timestep data in the dumpfile
     !
     dat(:,:,j) = 0.
     time(j) = 0.

     nread = 0
     do i=1,ncol
        read(15,end=44,iostat=ierr) dat(1:np,i,j)
        if (ierr /= 0) print*,' error reading column ',i
        nread = nread + 1
     enddo

44   continue

     if (nread.lt.ncol) then
        print "(a)",' WARNING: END OF FILE: read to column ',nread
     endif

     nstepsread = nstepsread + 1
     npartoftype(1,j) = np
     gamma(j) = 1.666666666667
     j = j + 1

  endif

  close(15)

  if (allocated(npartoftype)) then
     print*,'>> end of dump file: nsteps =',j-1,'ntot = ',sum(npartoftype(:,j-1))
  endif
return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,labelvec,labeltype,iamvec,&
              ix,ivx,ih,irho,iutherm,ipmass
  use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings
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
  irho = 7
  label(ix(1:ndim)) = labelcoord(1:ndim,1)

  if (ivx.ne.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = labelvec(ivx)//'\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 1
  labeltype(1) = 'gas'
  UseTypeInRenderings(1) = .true.

!-----------------------------------------------------------

  return
end subroutine set_labels
