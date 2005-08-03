!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR GENERAL ASCII DATA FORMATS
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

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data, only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,iunit,ncolstep
  integer :: nprint,npart_max,nstep_max,icol
  logical :: iexist
  character(len=len(rootname)+4) :: dumpfile

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  iunit = 15  ! logical unit number for input

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
  !--fix number of spatial dimensions (0 means no particle coords)
  !
  ndim = 0
  ndimV = 0

  j = indexstart
  nstepsread = 0
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the file and read the number of particles
  !
  open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='formatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
  else
     call get_ncolumns(iunit,ncolstep)
     if (ncolstep.le.0) return
     !
     !--allocate memory initially
     !
     nprint = 101
     nstep_max = max(nstep_max,indexstart,1)
     if (.not.allocated(dat) .or. (nprint.gt.npart_max)) then
        npart_max = max(npart_max,INT(1.1*(nprint)))
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
  endif

  npart_max = max(npart_max,nprint)
  ncolumns = ncolstep
!
!--allocate/reallocate memory if j > maxstep
!
  if (j.gt.maxstep) then
     call alloc(maxpart,j+1,maxcol)
  endif

!
!--now read the timestep data in the dumpfile
!
  i = 0
  ierr = 0
  overparts: do while (ierr == 0)
     i = i + 1
     if (i.gt.npart_max) then ! reallocate memory if necessary
        npart_max = 10*npart_max
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
     read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,ncolstep)
  enddo overparts

  nprint = i - 1
  nstepsread = nstepsread + 1
  if (ierr < 0) then
     print*,' end of file: npts = ',nprint
  elseif (ierr > 0) then
     print*,' *** error reading file, npts = ',nprint,' ***'
  endif


  npartoftype(:,j) = 0
  npartoftype(1,j) = nprint

  time(j) = 0. !!real(j)
  !!print*,' setting "time" = dump number = ',time(j)
  gamma(j) = 1.666666666667

  close(iunit)
     
return

contains
!
! utility to work out number of columns of real numbers
! in an ascii output file
!
! file must already be open and at the start
! slightly ad-hoc but its the best way I could think of!
!
subroutine get_ncolumns(lunit,ncolumns)
 implicit none
 integer, intent(in) :: lunit
 integer, intent(out) :: ncolumns
 integer :: ierr,i,nblanklines
 character(len=2000) :: line
 real :: dummyreal(100)

 nblanklines = 0
 line = ' '
 ierr = 0
 do while (len_trim(line).eq.0 .and. ierr.eq.0)
    read(lunit,"(a)",iostat=ierr) line
    nblanklines = nblanklines + 1
 enddo
 if (ierr .ne.0 ) then
    ncolumns = 0
    return
 else
    if (nblanklines.gt.1) print*,'skipped ',nblanklines-1,' blank lines'
    rewind(lunit)
 endif
 dummyreal = -666.0
 
 ierr = 0
 read(line,*,iostat=ierr) (dummyreal(i),i=1,size(dummyreal))
 if (ierr .gt. 0) then
    print "(a)",' WARNING: not all columns contain real numbers '
 endif

 i = 1
 ncolumns = 0
 do while(abs(dummyreal(i)+666.).gt.1.e-10)
    ncolumns = ncolumns + 1
    i = i + 1
    if (i.gt.size(dummyreal)) then
       print "(a)",'*** ERROR: too many columns in file'
       return
    endif
 enddo
 if (ncolumns.eq.0) then
    print "(a)",' ERROR: no columns of real numbers found'
 else
    print "(a,i3)",' number of data columns = ',ncolumns
 endif
end subroutine get_ncolumns
                   
end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!
!! * basically in this case we do nothing except guess that
!!   the first 3 columns are coordinates
!!
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,labeltype,ix
  use params
  use settings_data, only:ncolumns,ntypes
  use geometry, only:labelcoord
  implicit none
  integer :: i,ierr

  ix(:) = 0
!  do i=1,ndim
!     ix(i) = i
!  enddo
!  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  
!
!--read column labels from the columns file if it exists
!  
  open(unit=51,file='columns',status='old',iostat=ierr)
  if (ierr /=0) then
     print "(a)",' columns file not found: using default labels'
  else
     overcols: do i=1,ncolumns
        read(51,"(a)",iostat=ierr) label(i)
        if (ierr < 0) then
           print "(a,i3)",' ERROR: end of file in columns file: read to column ',i-1
           exit overcols
        elseif (ierr > 0) then
           print "(a)",' *** error reading from columns file ***'
           exit overcols
        endif
     enddo overcols
     close(unit=51)
  endif
  print "(3(/,a),/)",'WARNING: Rendering capabilities cannot be enabled', &
                 '         until positions of rho, h, pmass etc are', &
                 '         known (see read_data_ascii.f90 for details)'
  
!!  ivx = ndim+1
!!  ih = ndim+ndimV+2        !  smoothing length
!!  label(ih) = 'h'
!!  irho = ndim+ndimV+1     ! location of rho in data array
!!  label(irho) = 'density'      
!!  iutherm = 0  !  thermal energy
!!  label(iutherm) = 'u'
!!  ipmass = 0  !  particle mass
!!  label(ipmass) = 'particle mass'
  
!!  iamvec(ivx:ivx+ndimV-1) = ivx
!!  labelvec(ivx:ivx+ndimV-1) = 'v'
!!  do i=1,ndimV
!!     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
!!  enddo
  !
  !--set labels for each particle type
  !
  ntypes = 1 !!maxparttypes
  labeltype(1) = 'gas'
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
