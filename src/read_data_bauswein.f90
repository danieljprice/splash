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
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM ANDREAS BAUSWEIN'S CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! BSPLASH_R8 if 'YES' or 'TRUE' then assumes data is double precision
! BSPLASH_NCOL to change the number of columns read from the file,
! e.g. setenv BSPLASH_NCOL=22
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
  use settings_data, only:ndim,ndimV,ncolumns
  use mem_allocation, only:alloc
  use system_utils, only:lenvironment,ienvironment
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,k,ierr
  integer :: nprint,n1,npart_max,nstep_max
  integer :: ncol,nread,nerr,ncoltemp
  logical :: iexist,doubleprec
  character(len=len(rootname)) :: dumpfile
  real :: timei,dti
  real(doub_prec), dimension(maxcol) :: datdb
  real(doub_prec) :: timedb,dtdb

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
  ncol = 21
  doubleprec = .false.

  !--can override these settings with environment variables
  if (lenvironment('BSPLASH_R8')) doubleprec = .true.
  ncoltemp = ienvironment('BSPLASH_NCOL')
  if (ncoltemp.gt.0) ncol = ncoltemp
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,1)
  j = indexstart
  nstepsread = 0

  print "(1x,a)",'reading Andreas Bauswein format'
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
     if (doubleprec) then
        read(15,end=55,iostat=ierr) nprint,n1,timedb,dtdb
        timei = real(timedb)
        dti = real(dtdb)
     else
        read(15,end=55,iostat=ierr) nprint,n1,timei,dti
     endif
     print "(a,f10.2,a,i10,a,i10,a,f10.4)",' time: ',timei,' npart: ',nprint,' n1: ',n1,' dt = ',dti
     !--barf if stupid values read
     if (nprint.le.0 .or. nprint.gt.1e10) then
        print "(a)",' *** ERRORS IN TIMESTEP HEADER: WRONG ENDIAN? ***'
        close(15)
        return
     elseif (ierr /= 0) then
        print "(a)",'*** WARNING: ERRORS READING HEADER ***'
     endif
     if (timei.lt.0. .or. dti.lt.0.) print "(a)",'*** ERROR: t < 0: use setenv BSPLASH_R8=TRUE for double precision'
     ncolumns = ncol

     if (.not.allocated(dat) .or. nprint.gt.npart_max) then
        npart_max = max(npart_max,nprint)
        call alloc(npart_max,nstep_max,ncolumns)
     endif
     !
     !--now read the timestep data in the dumpfile
     !
     dat(:,:,j) = 0.
     time(j) = timei

     if (doubleprec) then
        nread = 0
        nerr = 0
        do i=1,nprint
           nread = nread + 1
           read(15,end=44,iostat=ierr) (datdb(k),k=1,ncol)
           if (ierr /= 0) nerr = nerr + 1
           dat(i,4,j) = real(datdb(1))
           dat(i,1:3,j) = real(datdb(2:4))
           dat(i,5:ncol,j) = real(datdb(5:ncol))
        enddo
     else
        nread = 0
        nerr = 0
        do i=1,nprint
           nread = nread + 1
           read(15,end=44,iostat=ierr) dat(i,4,j),dat(i,1:3,j),(dat(i,k,j),k=5,ncol)
           if (ierr /= 0) nerr = nerr + 1
        enddo
     endif

44   continue
     if (nerr.gt.0) print *,'*** WARNING: ERRORS DURING READ ON ',nerr,' LINES'

     if (nread.lt.nprint) then
        print "(a)",' WARNING: END OF FILE: read to particle ',nread
        nprint = nread
     endif

     nstepsread = nstepsread + 1
     npartoftype(1,j) = nprint
     gamma(j) = 1.666666666667
     j = j + 1

  endif

55 continue
  !
  !--reached end of file during header read
  !
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
  !use settings_units, only:units,unitslabel
  implicit none
  integer :: i,ipmom

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
  ih = 4
  ivx = 5
  ipmom = 8
  iutherm = 11
  ipmass = 14
  irho = 17

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(ih) = 'h'
  if (iutherm.gt.0) label(iutherm) = 'u'
  label(12) = 'psi\dpot\u'
  label(13) = 'alpha\dpot\u'
  label(ipmass) = 'particle mass'
  label(15) = '\gr'
  label(16) = 'P\deff'
  label(irho) = '\gr*'
  label(18) = 'tau'
  label(19) = '\ga\dvisc'
  label(20) = 'Ye'
  label(21) = 'temperature'

  if (ivx.ne.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = labelvec(ivx)//'\d'//labelcoord(i,1)
     enddo
  endif

  if (ipmom.ne.0) then
     iamvec(ipmom:ipmom+ndimV-1) = ipmom
     labelvec(ipmom:ipmom+ndimV-1) = 'momentum'
     do i=1,ndimV
        label(ipmom+i-1) = labelvec(ipmom)//'\d'//labelcoord(i,1)
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
