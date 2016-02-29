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
! THIS VERSION IS FOR FLASH HDF5 TRACER PARTICLES
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
module flash_hdf5read
  use, intrinsic :: iso_c_binding, only:c_int,c_float,c_char
 !interface to the c versions
 interface
   subroutine read_flash_hdf5_header(filename,time,npart,ncol,ierr) bind(c)
    import
    character(kind=c_char,len=1), intent(in) :: filename
    real(kind=c_float), intent(out) :: time
    integer(kind=c_int), intent(out) :: npart,ncol,ierr
   end subroutine read_flash_hdf5_header

   subroutine read_flash_hdf5_data(filename,npart,ncol,isrequired,ierr) bind(c)
    import
    character(kind=c_char,len=1), intent(in) :: filename
    integer(kind=c_int), intent(in) :: npart,ncol
    integer(kind=c_int), dimension(ncol), intent(in) :: isrequired
    integer(kind=c_int), intent(out) :: ierr
   end subroutine read_flash_hdf5_data
 end interface

contains

!
! function which maps from the order in which columns
! are read from the HDF5 file to the order in which they
! are stored in SPLASH. Differs because there are a couple
! of useless arrays that we do not read/store (ie. first column
! is on/off tag, 5th column is particle ID which we use to order
! the particles)
!
 integer function icolshuffle(icol)
   implicit none
   integer, intent(in) :: icol

   select case(icol)
   case(1)
      icolshuffle = 4
   case(2,3,4)
      icolshuffle = icol - 1
   case(5)
      icolshuffle = 0
   case default
      icolshuffle = icol
   end select

 end function icolshuffle

end module flash_hdf5read

subroutine read_data(dumpfile,indexstart,ipos,nstepsread)
  use particle_data, only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,lowmemorymode
  use mem_allocation, only:alloc
  use flash_hdf5read
  use asciiutils, only:cstring
  use labels, only:ih,irho
  use system_utils, only:renvironment
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: dumpfile

  integer :: i,j,ncolstep,ilastrequired
  integer :: nprint,npart_max,nstep_max,ierr
  integer, dimension(0:maxplot) :: isrequired
  logical :: iexist
  real :: tread,hfact,totmass

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
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

  j = indexstart
  nstepsread = 0
  print "(a)",' reading FLASH tracer particles (HDF5) data format '

  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

  call read_flash_hdf5_header(cstring(dumpfile),tread,nprint,ncolstep,ierr)
  ncolstep = ncolstep - 1   ! subtract particle ID column
  print "(a,i10,a,es10.3,a,i2)",' npart = ',nprint,' time = ',tread

  call set_labels
  if (ih.gt.0 .and. required(ih)) required(irho) = .true.
!
!--(re)allocate memory
!
  nstep_max = max(nstep_max,indexstart,1)
  if (.not.allocated(dat) .or. (nprint.gt.maxpart) .or. (ncolstep+ncalc).gt.maxcol) then
     npart_max = max(npart_max,nprint,maxpart)
     if (lowmemorymode) then
        ilastrequired = 0
        do i=1,ncolstep+ncalc
           if (required(i)) ilastrequired = i
        enddo
        call alloc(npart_max,j,ilastrequired)
     else
        call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol))
     endif
  endif
!
!--set the necessary parameters
!
  ncolumns = ncolstep
  nstepsread = nstepsread + 1
  npartoftype(:,j) = 0
  npartoftype(1,j) = nprint

  totmass = renvironment('FSPLASH_TOTMASS',-1.0)
  if (totmass.gt.0.) then
     print "(a,1pe10.3)",' setting total mass for all particles using FSPLASH_TOTMASS=',totmass
  else
     print "(a)",' FSPLASH_TOTMASS not set, assuming total mass of all particles is 1.0'
     totmass = 1.0
  endif
  masstype(1,j) = totmass/real(nprint)
  time(j) = tread
  gamma(j) = 5./3.
!
!--map "required" array to integers
!  also remap to the order read from the c data read
!
  isrequired(:) = 0
  do i=1,ncolstep
     if (icolshuffle(i).ne.0 .and. required(icolshuffle(i))) then
        !print*,'required '//trim(label(icolshuffle(i)))//' so must read ',i
        isrequired(i) = 1
     endif
  enddo

  if (.not.all(required(1:ncolstep))) then
     ipartialread = .true.
  else
     ipartialread = .false.
  endif
!
!--now read the timestep data in the dumpfile
!  (to avoid Fortran calling C with the array, we don't actually
!   pass the dat array here - instead we get c to
!   "call back" to fill the dat array, below)
!
  call read_flash_hdf5_data(cstring(dumpfile),nprint,ncolstep+1,isrequired(1:ncolstep+1),ierr)

  if (required(ih)) then
     hfact = 1.2
     hfact = renvironment('FSPLASH_HFACT',1.2)
     print "(a,i2,a,f5.2,a)",' creating smoothing length in column ',ih,' using h =',hfact,'(m/rho)^(1/3)'
     dat(1:nprint,ih,j) = hfact*(masstype(1,j)/dat(1:nprint,irho,j))**(1./3.)
  endif

return
end subroutine read_data

subroutine receive_data_fromc(icol,npart,temparr,id) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int,c_double
  use particle_data,  only:dat
  use flash_hdf5read, only:icolshuffle
  use labels,         only:label
  implicit none
  integer(kind=c_int), intent(in) :: icol,npart
  real(kind=c_double), dimension(npart), intent(in) :: temparr
  integer(kind=c_int), dimension(npart), intent(in) :: id
  integer(kind=c_int) :: i,icolput

  icolput = icolshuffle(icol)
  if (icolput.gt.size(dat(1,:,1)) .or. icolput.eq.0) then
     print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in receive_data_fromc'
     return
  endif
  print "(a,i2,a)",' reading column ',icol,' -> '//trim(label(icolput))

  do i=1,npart
     if (id(i).lt.1 .or. id(i).gt.npart) then
        print*,' ERROR in particle id = ',id(i)
     else
        dat(id(i),icolput,1) = real(temparr(i))
     endif
  enddo

  return
end subroutine receive_data_fromc

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
  use settings_data, only:ntypes,ndim,ndimV,UseTypeInRenderings,ncolumns
  use geometry, only:labelcoord
  implicit none
  integer :: i

  do i=1,ndim
     ix(i) = i
     label(i) = labelcoord(i,1)
  enddo
  ih = 5
  ivx = 6
  ipmass = 0
  label(4) = 'density (from grid)'
  label(ih) = 'smoothing length'
  irho = 4
  if (ncolumns.ge.9) then
     irho = 9
     label(9) = 'density (on particles)'
  else
     irho = 4
  endif

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
  ntypes = 1
  labeltype(1) = 'tracer'
  UseTypeInRenderings(1) = .true.

!-----------------------------------------------------------

  return
end subroutine set_labels
