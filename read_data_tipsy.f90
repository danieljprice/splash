!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING TIPSY FILES
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
  use particle_data, only:dat,time,npartoftype,gamma,maxpart
  use params
  use settings_data, only:ndim,ndimV,ncolumns
  use mem_allocation, only:alloc
  use labels, only:label,ih
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,ic,icol
  integer :: nprint,ngas,ndark,nptmass,npart_max,nstep_max
  integer :: ncol,nerr,nread
  logical :: iexist
  character(len=len(rootname)) :: dumpfile
  real :: timei

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

  nstep_max = max(nstep_max,indexstart,1)
  j = indexstart
  nstepsread = 0
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,file=dumpfile,status='old',form='formatted',iostat=ierr)
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  else
     !
     !--read the number of particles in the header and allocate memory
     !
     read(15,*,end=55,iostat=ierr) nprint,ngas,nptmass
     read(15,*,end=55,iostat=ierr) ndim
     read(15,*,end=55,iostat=ierr) timei
     ndark = nprint - ngas - nptmass
     print "(a,f10.2,1(a,i1))",' time: ',timei,' ndim: ',ndim
     print "(4(a,i10))",' ntot: ',nprint,' ngas: ',ngas,' ndark: ',ndark,' nstar: ',nptmass
     !--barf if stupid values read
     if (nprint.le.0 .or. nprint.gt.1e10 .or. ndim.le.0 .or. ndim.gt.3 .or. ndark.lt.0) then
        print "(a)",' *** ERROR READING TIMESTEP HEADER ***'
        close(15)
        return
     elseif (ierr /= 0) then
        print "(a)",'*** WARNING: ERRORS READING HEADER ***'  
     endif
     ndimV = ndim
     ncol = 2*ndim + 4
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

     nread = 0
     call set_labels
     
     !--pmass,x,y,z,vx,vy,vz
     do ic=1,2*ndim+1
        nerr = 0
        nread = nread + 1
        if (ic.eq.1) then ! pmass
           icol = ndim + 1
        elseif (ic.ge.2 .and. ic.le.ndim+1) then ! x, y, z
           icol = ic - 1
        else ! everything after
           icol = ic
        endif
        !print "(1x,a)",trim(label(icol))
        nerr = 0
        do i=1,nprint
           read(15,*,end=44,iostat=ierr) dat(i,icol,j)
           if (ierr /= 0) nerr = nerr + 1
        enddo
        if (nerr.gt.0) print "(/,a)",'*** WARNING: ERRORS READING '//trim(label(icol))//' ON ',nerr,' LINES'
     enddo
     !--h dark matter
     if (ndark.gt.0) then
        nerr = 0
        do i=ngas+1,ngas+ndark-1
           read(15,*,end=44,iostat=ierr) dat(i,ih,j)
           if (ierr /= 0) nerr = nerr + 1
        enddo
        if (nerr.gt.0) print *,'*** WARNING: ERRORS READING DARK MATTER H ON ',nerr,' LINES'
     endif
     !--h star particles
     if (nptmass.gt.0) then
        nerr = 0
        do i=ngas+ndark+1,ngas+ndark+nptmass
           read(15,*,end=44,iostat=ierr) dat(i,ih,j)
           if (ierr /= 0) nerr = nerr + 1
        enddo
        if (nerr.gt.0) print *,'*** WARNING: ERRORS READING PTMASS H ON ',nerr,' LINES'
     endif
     !--density, temperature, sph smoothing length
     do icol=2*ndim+2,ncol
        nread = nread + 1
        !print "(1x,a)",trim(label(icol))
        do i=1,ngas
           read(15,*,end=44,iostat=ierr) dat(i,icol,j)
           if (ierr /= 0) nerr = nerr + 1
        enddo
        if (nerr.gt.0) print "(/,a)",'*** WARNING: ERRORS READING '//trim(label(icol))//' ON ',nerr,' LINES'
     enddo
     
     !--multiply sph smoothing length x 2
     dat(1:ngas,ih,j) = dat(1:ngas,ih,j)*2.

44   continue
     
     if (nread.lt.ncol) then
        print "(a,i2)",' WARNING: END OF FILE: READ TO COLUMN ',nread
        ncolumns = nread
     endif

     nstepsread = nstepsread + 1
     npartoftype(1,j) = ngas
     npartoftype(2,j) = ndark
     npartoftype(3,j) = nptmass
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
  ipmass = ndim + 1
  ivx = ndim + 2
  irho = ivx + ndim
  iutherm = irho + 1
  ih = iutherm + 1
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(ih) = 'h'
  if (iutherm.gt.0) label(iutherm) = 'u'
  label(ipmass) = 'particle mass'
  label(irho) = 'density'

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
  ntypes = 3
  labeltype(1) = 'gas'
  labeltype(2) = 'dark matter'
  labeltype(3) = 'star'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.
  UseTypeInRenderings(3) = .false.
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
