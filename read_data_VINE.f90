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

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname

  integer, parameter :: iheadlength = 1000
  integer :: i,j,ierr,nparti,ntoti
  integer :: npart_max,nstep_max,ncolstep
  logical :: iexist
    
  character(len=len(rootname)+10) :: dumpfile
  integer, dimension(iheadlength) :: iheader
  integer, dimension(:), allocatable :: ipindx
  
  !--we are assuming dump is double precision
  real(doub_prec), dimension(iheadlength) :: dheader
  real(doub_prec), dimension(:,:), allocatable :: dattemp, dattempvec

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
  ncolstep = 12  ! number of columns in file
  ncolumns = ncolstep
  !
  !--allocate memory initially
  !
  nstep_max = max(indexstart,1)

  j = indexstart
  nstepsread = 0
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
  else
     !
     !--read timestep header (integers only)
     !
     read(15,iostat=ierr) (iheader(i),i=1,iheadlength)
     !
     !--get number of particles from header and allocate memory
     !
     ntoti = iheader(2)
     nparti = iheader(3)
     
     if (.not.allocated(dat) .or. ntoti.gt.npart_max) then
        npart_max = max(npart_max,INT(1.1*ntoti))
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
     !
     !--rewind file
     !
     rewind(15)
  endif
  if (ierr /= 0) then
     print "(a)",'*** ERROR READING TIMESTEP HEADER ***'
  else

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
       allocate(dattemp(npart_max,ncolstep),stat=ierr)
       dattemp = 0.
       if (ierr /= 0) print*,'not enough memory in read_data'
!
!--allocate a temporary array for vectors
!
       if (allocated(dattempvec)) deallocate(dattempvec)
       allocate(dattempvec(7,npart_max),stat=ierr)
       dattempvec = 0.
       if (ierr /= 0) print*,'not enough memory in read_data'
!
!--allocate a temporary array for particle index
!
       if (allocated(ipindx)) deallocate(ipindx)
       allocate(ipindx(npart_max),stat=ierr)
       ipindx = 0
       if (ierr /= 0) print*,'not enough memory in read_data'
!
!--now read the timestep data in the dumpfile
!
       write(*,"(a,i5,a)",advance="no") '| step ',j,': '

       read(15,iostat=ierr), &
            (iheader(i),i=1,iheadlength), &
            (dheader(i),i=1,iheadlength), &
            (dattempvec(1:4,i),i=1,ntoti), &
            (dattempvec(5:7,i),i=1,ntoti), &
            (dattemp(i,8), i=1,ntoti), &
            (dattemp(i,9), i=1,nparti), &
            (dattemp(i,10), i=1,nparti), &
            (dattemp(i,11), i=1,nparti), &
            (dattemp(i,12), i=1,ntoti), &
            (ipindx(i), i=1,ntoti)

       if (ierr < 0) then
          print "(a)",'*** END OF FILE IN READ DATA (CHECK PRECISION) ***'
          close(15)
          return
       elseif (ierr /= 0) then
          print "(a)",'*** ERROR READING DATA ***'
          close(15)
          return
       else
          nstepsread = nstepsread + 1
       endif
!
!--spit out time
!
       time(j) = real(dheader(1))
       gamma(j) = real(dheader(4))
       print "(a,f8.3,a,i8)",'t = ',time(j),' ntotal = ',ntoti
!
!--convert vectors to columns and double to single precision
!
       do i=1,7
          dat(ipindx(1:ntoti),i,j) = real(dattempvec(i,1:ntoti))
       enddo
!
!--now convert scalars
!
       dat(ipindx(1:ntoti),8:ncolstep,j) = real(dattemp(1:ntoti,8:ncolstep))

!
!--set particle numbers
!
       npartoftype(1,j) = nparti
       npartoftype(2,j) = ntoti - nparti
!
!--clean up
!
       if (allocated(dattemp)) deallocate(dattemp)
       if (allocated(dattempvec)) deallocate(dattempvec)
       if (allocated(ipindx)) deallocate(ipindx)

  endif

 !
 !--reached end of file
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
  ivx = 5
  ih = 8        !  smoothing length
  iutherm = 9  !  thermal energy
  ipmass = 4   !  particle mass      
  irho = 10     ! location of rho in data array
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = '\gr'
  label(iutherm) = 'u'
  label(ih) = 'h'
  label(ipmass) = 'particle mass'
  label(11) = 'alpha'
  label(12) = 'poten'
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
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'Nbody'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
