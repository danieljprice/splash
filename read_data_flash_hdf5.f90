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

 !interface to the c versions
 interface
   subroutine read_flash_hdf5_header(filename,time,npart,ncol,ierr) bind(c)
    character(len=*), intent(in) :: filename
    real, intent(out) :: time
    integer, intent(out) :: npart,ncol,ierr
   end subroutine read_flash_hdf5_header

   subroutine read_flash_hdf5_data(filename,npart,ncol,ierr) bind(c)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: npart,ncol
    integer, intent(out) :: ierr
   end subroutine read_flash_hdf5_data
 end interface

end module flash_hdf5read

subroutine read_data(dumpfile,indexstart,nstepsread)
  use particle_data, only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation, only:alloc
  use flash_hdf5read
  use asciiutils, only:cstring
  use labels, only:ih,irho
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: dumpfile
  
  integer :: j,ncolstep
  integer :: nprint,npart_max,nstep_max,ierr
  logical :: iexist
  real :: tread,hfact

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
  print "(a,i10,a,1pe10.3)",' npart = ',nprint,' time = ',tread
!
!--(re)allocate memory
!
  nstep_max = max(nstep_max,indexstart,1)
  if (.not.allocated(dat) .or. (nprint.gt.maxpart) .or. (ncolstep+ncalc).gt.maxcol) then
     npart_max = max(npart_max,nprint,maxpart)
     call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol))
  endif
!
!--set the necessary parameters
!
  ncolumns = ncolstep
  nstepsread = nstepsread + 1
  npartoftype(:,j) = 0
  npartoftype(1,j) = nprint
  masstype(1,j) = 1.0/real(nprint)
  time(j) = tread
  gamma(j) = 5./3.
!
!--now read the timestep data in the dumpfile
!  (to avoid Fortran calling C with the array, we don't actually
!   pass the dat array here - instead we get c to
!   "call back" to fill the dat array, below)
!
  call read_flash_hdf5_data(cstring(dumpfile),nprint,ncolstep,ierr)
  
  call set_labels
  hfact = 1.2
  print "(a,i2)",' creating smoothing length in column ',ih
  dat(1:nprint,ih,j) = hfact*(masstype(1,j)/dat(1:nprint,irho,j))**(1./3.)
  
     
return
end subroutine read_data

subroutine receive_data_fromc(icol,npart,temparr,id) bind(c)
  use particle_data, only:dat
  implicit none
  integer, intent(in) :: icol,npart
  double precision, dimension(npart), intent(in) :: temparr
  integer, dimension(npart), intent(in) :: id
  integer :: i,icolput
  
  select case(icol)
  case(1)
     icolput = 4
  case(2,3,4)
     icolput = icol - 1
  case default
     icolput = icol
  end select
  print*,'reading column = ',icol+1,'->',icolput

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
  use settings_data, only:ntypes,ndim,ndimV,UseTypeInRenderings
  use geometry, only:labelcoord
  implicit none
  integer :: i

  do i=1,ndim
     ix(i) = i
  enddo
  irho = 4
  ih = 9
  ivx = 6
  ipmass = 0
  label(irho) = 'density'
  label(ih) = 'smoothing length'
  
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
