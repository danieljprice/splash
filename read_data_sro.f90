!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM MATTHEW BATE'S CODE
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
! ntot(maxstep)       : total number of particles in each timestep
! iam(maxpart,maxstep): integer identification of particle type
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
  use settings_data  
  use mem_allocation
  implicit none
  integer, intent(IN) :: indexstart
  integer, intent(OUT) :: nstepsread
  character(LEN=*), intent(IN) :: rootname
  integer, parameter :: maxptmass = 10
  integer :: i,j,ifile,ierr
  integer :: nprint,nptmass,npart_max,nstep_max
  logical :: iexist,magfield   
  character(LEN=LEN(rootname)) :: dumpfile
  real :: timei

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  ifile = 1
  magfield = .true.

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
  if (magfield) then
     ncolumns = 11
  else
     ncolumns = 7  ! number of columns in file  
  endif
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,2)

  j = indexstart
  nstepsread = 0
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
     open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
     if (ierr /= 0) then
        print*,'*** ERROR OPENING ',trim(dumpfile),' ***'
     else
        !
        !--read the number of particles in the first step,
        !  allocate memory and rewind
        !
        read(15,end=55,iostat=ierr) timei,nprint,nptmass
        print*,'first time = ',timei,nprint,nptmass
        if (.not.allocated(dat) .or. (nprint+nptmass).gt.npart_max) then
           npart_max = max(npart_max,INT(1.1*(nprint+nptmass)))
           call alloc(npart_max,nstep_max,ncolumns)
        endif
        rewind(15)
     endif
     if (ierr /= 0) then
        print*,'*** ERROR READING TIMESTEP HEADER ***'
     else
!
!--loop over the timesteps in this file
!     
        npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
        if (j.gt.maxstep) then
           call alloc(maxpart,j+1,maxcol)
        endif
!
!--now read the timestep data in the dumpfile
!
        if (magfield) then
           read(15,end=55,iostat=ierr) time(j),nprint,nptmass, &
             (dat(i,1,j),i=1,nprint),(dat(i,2,j),i=1,nprint),  &
             (dat(i,3,j),i=1,nprint),(dat(i,4,j),i=1,nprint),  &
             (dat(i,5,j),i=1,nprint),(dat(i,6,j),i=1, nprint), &
             (dat(i,8,j),i=1,nprint),(dat(i,9,j),i=1,nprint),  &
             (dat(i,10,j),i=1,nprint),(dat(i,11,j),i=1,nprint),&             
             (dat(i,7,j), i=nprint+1, nprint+nptmass), &
             (dat(i,1,j), i=nprint+1, nprint+nptmass), &
             (dat(i,2,j), i=nprint+1, nprint+nptmass), &
             (dat(i,3,j), i=nprint+1, nprint+nptmass)        
        else
           read(15,end=55,iostat=ierr) time(j),nprint,nptmass, &
             (dat(i,1,j), i=1, nprint), (dat(i,2,j), i=1,nprint), &
             (dat(i,3,j), i=1, nprint), (dat(i,4,j), i=1,nprint), &
             (dat(i,5,j), i=1, nprint), (dat(i,6,j), i=1, nprint), &
             (dat(i,7,j), i=nprint+1, nprint+nptmass), &
             (dat(i,1,j), i=nprint+1, nprint+nptmass), &
             (dat(i,2,j), i=nprint+1, nprint+nptmass), &
             (dat(i,3,j), i=nprint+1, nprint+nptmass)
        endif
             
        dat(1:nprint,7,j) = 1.4/real(2971627)
             
        if (ierr /= 0) then
           print "(a)",'|*** ERROR READING TIMESTEP ***'
           return
        else
           nstepsread = nstepsread + 1
        endif

        ntot(j) = nprint+nptmass
        npartoftype(1,j) = nprint
        npartoftype(2,j) = nptmass

        print*,j,' time = ',time(j)
        gamma(j) = 1.666666666667
        j = j + 1
     
     endif

55 continue
  !
  !--reached end of file
  !
  close(15)

  print*,'>> end of dump file: nsteps =',j-1,'ntot = ',ntot(j-1),'nptmass=',npartoftype(2,j-1)
   
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
  ivx = 0
  ih = 4        !  smoothing length
  irho = 5     ! location of rho in data array
  iutherm = 6  !  thermal energy
  ipmass = 7  !  particle mass
  if (ncolumns.gt.7) then
     iBfirst = 8
     idivB = 11 
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = labelvec(iBfirst)//'\d'//labelcoord(i,1)
     enddo
     label(11) = 'div B'
  endif
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
!  do i=1,ndimV
!     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
!  enddo
  label(irho) = '\gr'      
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'
  !
  !--set labels for each particle type
  !
  ntypes = 2 !!maxparttypes
  labeltype(1) = 'gas'
  labeltype(2) = 'point mass'
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
