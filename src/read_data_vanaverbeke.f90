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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT
! FROM SIGFRIED VANAVERBEKE'S CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! VSPLASH_SINGLEPREC if 'YES' or 'TRUE' then assumes data is single precision
! VSPLASH_NCOL to change the number of columns read from the file,
! e.g. setenv VSPLASH_NCOL=13
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
  use particle_data, only:dat,time,npartoftype,gamma,maxpart
  use params
  use settings_data, only:ndim,ndimV,ncolumns
  use mem_allocation, only:alloc
  use system_utils, only:lenvironment,ienvironment
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,k,ierr
  integer :: nprint,ntotal,npart_max,nstep_max
  integer :: ncol,nread,nerr,ncoltemp,nsink,nsinkcol,nacc
  logical :: iexist,doubleprec
  integer, parameter :: iunit = 15
  character(len=len(rootname)) :: dumpfile
  real :: timei,gammai
  real(doub_prec), dimension(maxplot) :: datdb
  real(doub_prec) :: timedb,gammadb

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
  ncol = 10
  nsinkcol = 7
  doubleprec = .true.

  !--can override these settings with environment variables
  if (lenvironment('VSPLASH_SINGLEPREC')) doubleprec = .false.
  ncoltemp = ienvironment('VSPLASH_NCOL')
  if (ncoltemp.gt.0) ncol = ncoltemp
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,1)
  j = indexstart
  nstepsread = 0
  print "(a)",' reading Vanaverbeke format...'

  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=iunit,file=dumpfile,status='old',form='unformatted',iostat=ierr)
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  else
     timei = 0.
     read(iunit,iostat=ierr) nprint,nacc,nsink
     print "(3(a,i10))",'npart:',nprint,' naccreted:',nacc,' nsinks:',nsink
     if (doubleprec) then
        read(iunit,iostat=ierr) timedb,gammadb
        timei = real(timedb)
        gammai = real(gammadb)
     else
        read(iunit,iostat=ierr) timei,gammai
     endif
     print "(2(a,1pe12.3))",'time:',timei,' gamma:',gammai

     !--barf if stupid values read
     if (nprint.lt.0 .or. nacc.lt.0 .or. nsink.lt.0 .or. nsink.gt.1e7) then
        print "(a)",' *** ERROR IN TIMESTEP HEADER: WRONG ENDIAN? (or old header format)?'
        close(iunit)
        return
     elseif (ierr /= 0) then
        print "(a)",'*** ERROR READING TIMESTEP HEADER: WRONG ENDIAN? ***'
        close(iunit)
        return
     endif
     if (timei.lt.0. .or. gammai.lt.1.0 .or. gammai.gt.2.0) then
        print*,'*** ERROR IN HEADER: strange time and/or gamma read: wrong precision?'
     endif
     ncolumns = ncol

     ntotal = nprint + nsink
     if (.not.allocated(dat) .or. ntotal.gt.npart_max) then
        npart_max = max(npart_max,ntotal)
        call alloc(npart_max,nstep_max,ncolumns)
     endif
     !
     !--now read the timestep data in the dumpfile
     !
     dat(:,:,j) = 0.
     time(j) = timei
     gamma(j) = gammai

     if (doubleprec) then
        nerr = 0
        nread = 0
        do i=1,nprint
           nread = nread + 1
           read(iunit,end=44,iostat=ierr) (datdb(k),k=1,ncol)
           if (ierr /= 0) then
              nerr = nerr + 1
           else
              dat(i,1:ncol,j) = real(datdb(1:ncol))
           endif
        enddo
        do i=nprint+1,nprint+nsink
           nread = nread + 1
           read(iunit,end=44,iostat=ierr) (datdb(k),k=1,nsinkcol)
           if (ierr /= 0) then
              nerr = nerr + 1
           else
              dat(i,1:nsinkcol,j) = real(datdb(1:nsinkcol))
           endif
        enddo
     else
        nerr = 0
        nread = 0
        do i=1,nprint
           nread = nread + 1
           read(iunit,end=44,iostat=ierr) dat(i,1:ncol,j)
           if (ierr /= 0) nerr = nerr + 1
        enddo
        do i=nprint+1,nprint+nsink
           nread = nread + 1
           read(iunit,end=44,iostat=ierr) dat(i,1:nsinkcol,j)
           if (ierr /= 0) nerr = nerr + 1
        enddo
     endif

     goto 45
44   continue
     print "(a,i10)",' WARNING: END-OF-FILE AT LINE ',nread
45   continue
     if (nerr.gt.0) print *,'*** WARNING: ERRORS DURING READ ON ',nerr,' LINES'

     nstepsread = nstepsread + 1
     npartoftype(1,j) = nprint - nacc
     npartoftype(2,j) = nacc
     npartoftype(3,j) = nsink
     j = j + 1

  endif

  close(iunit)

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
              ix,ivx,ih,irho,iutherm,ipmass,ispsound
  use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings
  use geometry, only:labelcoord
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
  ivx = 4
  ipmass = 7
  ih = ipmass + 1
  iutherm = ih + 1
  irho = iutherm + 1
  ispsound = 0

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(ih) = 'h'
  if (irho.gt.0) label(irho) = '\gr'
  if (ipmass.gt.0) label(ipmass) = 'particle mass'
  if (iutherm.gt.0) label(iutherm) = 'u'
  if (ispsound.gt.0) label(ispsound) = 'c_s'

  if (ivx.ne.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = trim(labelvec(ivx))//'\d'//trim(labelcoord(i,1))
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 3
  labeltype(1) = 'gas'
  labeltype(2) = 'accreted/dead'
  labeltype(3) = 'sink'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .true.
  UseTypeInRenderings(3) = .false.

!-----------------------------------------------------------

  return
end subroutine set_labels
