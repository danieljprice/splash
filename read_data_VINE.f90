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
  use particle_data, only:npartoftype,dat,time,gamma,maxcol,maxpart,maxstep
  use params, only:doub_prec
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use labels, only:ivx, iBfirst
  use mem_allocation, only:alloc
  use system_commands, only:lenvironment
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname

  integer :: iheadlength
  integer :: i,j,ierr,nparti,ntoti,i1,icol
  integer :: npart_max,nstep_max,ncolstep
  logical :: iexist,mhdread
    
  character(len=len(rootname)+10) :: dumpfile
  integer, parameter :: maxheadlength = 1000
  integer, dimension(maxheadlength) :: iheader
  integer, dimension(:), allocatable :: ipindx,itstepbin
  
  !--we are assuming dump is double precision
  real(doub_prec), dimension(maxheadlength) :: dheader
  real(doub_prec), dimension(:,:), allocatable :: dattemp, dattempvec

  nstepsread = 0
  npart_max = maxpart
  !--this is the default header length
  iheadlength = maxheadlength

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
  mhdread = .false.
  if (lenvironment('VINE_MHD')) then
     mhdread = .true.
  endif
  
  nstep_max = max(indexstart,1)

  j = indexstart
  nstepsread = 0
  nparti = 0
  ncolstep = 0
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
  else
     if (mhdread) print "(a)",' assuming MHD file from VINE_MHD setting'
     !
     !--read timestep header (integers only)
     !
     read(15,iostat=ierr) (iheader(i),i=1,iheadlength)
     !
     !--get number of particles from header and allocate memory
     !
     iheadlength = iheader(1)
     if (iheadlength.gt.maxheadlength) print "(a)",' ERROR: header length too big!'
     ntoti = iheader(2)
     nparti = iheader(3)
     ndim = iheader(30)
     ndimV = ndim
     if (mhdread) then
        ncolstep = 2*ndim + 6 + ndim
     else
        ncolstep = 2*ndim + 6
     endif
     ncolumns = ncolstep
     if (.not.allocated(dat) .or. ntoti.gt.npart_max) then
        if (.not.allocated(dat)) then
           npart_max = ntoti
        else
           npart_max = max(npart_max,INT(1.1*ntoti))
        endif
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
       if (mhdread) then
          allocate(dattempvec(3*ndim+1,npart_max),stat=ierr)
       else
          allocate(dattempvec(2*ndim+1,npart_max),stat=ierr)       
       endif
       dattempvec = 0.
       if (ierr /= 0) print*,'not enough memory in read_data'
!
!--allocate a temporary array for particle index
!
       if (allocated(ipindx)) deallocate(ipindx)
       allocate(ipindx(npart_max),stat=ierr)
       !ipindx = 0
       if (ierr /= 0) print*,'not enough memory in read_data'
!
!--allocate a temporary array for itstepbin (MHD only)
!
       if (mhdread) then
          if (allocated(itstepbin)) deallocate(itstepbin)
          allocate(itstepbin(npart_max),stat=ierr)
          !itstepbin = 0
          if (ierr /= 0) print*,'not enough memory in read_data'
       endif
!
!--now read the timestep data in the dumpfile
!
       write(*,"(a,i5,a)",advance="no") '| step ',j,': '

       ivx = ndim + 2 ! location of vx in 'columns'
!      starting point for non position and velocity columns
       icol = ndim + 1 + ndimV + 1

       if (mhdread) then
          iBfirst = icol+5
          read(15,iostat=ierr), &
               (iheader(i),i=1,iheadlength), &
               (dheader(i),i=1,iheadlength), &
               (dattempvec(1:ndim+1,i),i=1,ntoti), &
               (dattempvec(ivx:ivx+ndimV-1,i),i=1,ntoti), &
               (dattemp(i,icol), i=1,ntoti), &
               (dattemp(i,icol+1), i=1,nparti), &
               (dattemp(i,icol+2), i=1,nparti), &
               (dattemp(i,icol+3), i=1,nparti), &
               (dattemp(i,icol+4), i=1,ntoti), &
               (ipindx(i), i=1,ntoti), &
               (itstepbin(i),i=1,ntoti), &
               (dattempvec(ivx+ndimV:ivx+2*ndimV-1,i),i=1,nparti)       
       else
          read(15,iostat=ierr), &
               (iheader(i),i=1,iheadlength), &
               (dheader(i),i=1,iheadlength), &
               (dattempvec(1:ndim+1,i),i=1,ntoti), &
               (dattempvec(ivx:ivx+ndimV-1,i),i=1,ntoti), &
               (dattemp(i,icol), i=1,ntoti), &
               (dattemp(i,icol+1), i=1,nparti), &
               (dattemp(i,icol+2), i=1,nparti), &
               (dattemp(i,icol+3), i=1,nparti), &
               (dattemp(i,icol+4), i=1,ntoti), &
               (ipindx(i), i=1,ntoti)
       endif

       if (ierr < 0) then
          print "(a)",'*** END OF FILE IN READ DATA ***'
       elseif (ierr /= 0) then
          if (mhdread) then
             print "(a)",'*** ERROR READING DATA: MAYBE NOT AN MHD FILE?? ***'          
          else
             print "(a)",'*** ERROR READING DATA ***'
          endif
       endif
       nstepsread = nstepsread + 1
!
!--spit out time
!
       time(j) = real(dheader(1))
       gamma(j) = real(dheader(4))
       print "(a,f8.3,2(a,i8))",'t = ',time(j),' n(SPH) = ',ntoti,' n(Nbody) = ',ntoti-nparti
!
!--convert posm and velocity vectors to columns and double to single precision
!
       do i=1,2*ndim+1
          dat(ipindx(1:ntoti),i,j) = real(dattempvec(i,1:ntoti))
       enddo
!
!--convert B vectors to columns and double to single precision
!
       if (mhdread) then
          i1 = iBfirst - 1
          do i=ivx+ndimV,ivx+2*ndimV-1
             i1 = i1 + 1
             dat(ipindx(1:nparti),i1,j) = real(dattempvec(i,1:nparti))
          enddo
       endif
!
!--now convert scalars
!
       dat(ipindx(1:ntoti),icol:ncolstep,j) = real(dattemp(1:ntoti,icol:ncolstep))

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
       if (allocated(itstepbin)) deallocate(itstepbin)

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
  use labels, only:label,ih,ipmass,ivx,iutherm,irho,ix,iBfirst, &
                   labelvec,iamvec,labeltype
  use params
  use settings_data, only:ndim,ndimV,UseTypeInRenderings,ntypes
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
  ipmass = ndim+1   !  particle mass      
  ivx = ndim+2
  
  ih = ndim + 1 + ndimV + 1       !  smoothing length
  iutherm = ih+1  !  thermal energy
  irho = iutherm+1     ! location of rho in data array
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = 'density'
  label(iutherm) = 'u'
  label(ih) = 'h'
  label(ipmass) = 'particle mass'
  label(irho+1) = 'alpha'
  label(irho+2) = 'poten'
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
     enddo
  endif
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
