!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nfilesteps  : number of steps read from this file
! hfact       : constant relating smoothing length to particle spacing
! ivegotdata  : flag which indicates successful data read
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npart(maxstep)      : number of particles in each timestep
! ntot(maxstep)       : total number of particles in each timestep
! nghost(maxstep)     : number of ghost particles in each timestep
! iam(maxpart,maxstep): integer identification of particle type
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,nfilesteps)
  use particle_data
  use params
  use labels
  use settings
  implicit none
  integer, intent(OUT) :: nfilesteps
  character(LEN=20) :: datfile
  character(LEN=2) :: fileno
  character(LEN=*), intent(IN) :: rootname
  integer :: i,j,k,ifile
  integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
  integer :: npartin,ntotin
  logical :: iexist,reallocate
  real :: timein,gammain

  if (rootname(1:1).ne.' ') then
     !
     !--if rootname does not contain .dat, make it end in .dat
     !
     if (index(rootname,'.dat').eq.0) then
        datfile = trim(rootname)//'.dat'
     else
        datfile = trim(rootname)  
     endif
     ifile = 1
     print*,'rootname = ',rootname
  else
     print*,' **** no data read **** ' 
     return
  endif

  print *,' opening ',datfile
  k=1      

50 continue
  ivegotdata = .false.
  !
  !--set everything to zero initially
  !
  dat = 0.
  time = 0.
  gamma = 0.
  hfact = 0.
  iam = 0
  npart = 0
  ntot = 0
  nghost = 0
  ndim = 0
  ndimV = 0
  ncolumns = 0
  ncalc = 0
  nfilesteps = 0
  !
  !--open data file and read data
  !
  open(unit=11,ERR=81,file=datfile,status='old',form='formatted')

  nfilesteps = 100000
!
!--read first header line
!
  read(11,*,ERR=78,end=80) timein,npartin,ntotin,gammain, &
       hfact,ndim,ndimV,ncol_max	 
!  print*,'reading time = ',timein,npartin,ntotin,gammain, &
!       ndim_max,ndimV_max,ncol_max     
!
!--allocate memory for data arrays (initially for 11 timesteps)
!
  if (ntotin.lt.5500) then
     nstep_max = 111
  elseif (ntotin.lt.111111) then
     nstep_max = 11
  else 
     nstep_max = 5
  endif
  npart_max = ntotin
  if (.not.allocated(dat) .or. ntotin.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncol_max.gt.maxcol) then
     call alloc(ntotin,nstep_max,ncol_max)
  endif
!
!--rewind file
!
  rewind(11)

  do i=1,nfilesteps
     reallocate = .false.
     !
     !--read header line for this timestep
     !
     read(11,*,ERR=78,end=67) time(i),npart(i),ntot(i),gamma(i), &
          hfact,ndim,ndimV,ncolumns	 
     print*,'reading time = ',time(i),npart(i),ntot(i),gamma(i), &
          ndim,ndimV,ncolumns     
     if (ncolumns.ne.ncol_max) then
        print*,'*** Warning number of columns not equal for timesteps'
        print*,'ncolumns = ',ncolumns,ncol_max
     endif
     if (ncolumns.gt.ncol_max) then
        reallocate = .true.
        ncol_max = ncolumns
     endif

     if (ndim.gt.ndim_max) ndim_max = ndim
     if (ndim_max.gt.3) stop 'error: ndim in file> 3'
     if (ndimV.gt.ndimV_max) ndimV_max = ndimV  
     if (ndimV_max.gt.3) stop 'error: ndimV in file> 3' 
     nghost(i) = ntot(i) - npart(i)
     if (ntot(i).gt.maxpart) then
        !print*, 'ntot greater than array limits!!'    
        reallocate = .true.
	npart_max = ntot(i)
     endif
     if (i.eq.nstep_max) then
        nstep_max = i + 10
        reallocate = .true.
     endif
     !
     !--reallocate memory for main data array
     !
     if (reallocate) then
        call alloc(npart_max,nstep_max,ncol_max)
     endif

  
     if (ntot(i).gt.0) then

        read (11,*, end=66,ERR=77) (dat(1:ncolumns,j,i),j=1,ntot(i))

     else
        ntot(i) = 1
        npart(i) = 1
        nghost(i) = 0
        dat(:,:,i) = 0.
     endif
     iam(:,i) = 0
  enddo

  print*,' REACHED ARRAY LIMITS IN READFILE'

  nfilesteps = i-1		! this is if reached array limits
  ntot(i-1) = j-1
  nghost(i-1) = ntot(i-1) - npart(i-1)
  goto 68

66 continue
  nfilesteps = i		! timestep there but data incomplete
  ntot(i) = j-1
  nghost(i) = ntot(i) - npart(i)
  goto 68

67 continue
  nfilesteps = i-1		! no timestep there at all

68 continue
  !
  !--close data file and return
  !      	      
  close(unit=11)

  ivegotdata = .true.
  ncolumns = ncol_max
  ndim = ndim_max
  ndimV = ndimV_max
  print*,'ncolumns = ',ncolumns

  print*,'>> READ all steps =',nfilesteps,'last step ntot = ',ntot(nfilesteps)

  !!------------------------------------------------------------
  !! set labels for each column of data

  do i=1,ndim
     ix(i) = i
  enddo
  ivx = ndim + 1
  ivlast = ndim + ndimV
  irho = ndim + ndimV + 1		! location of rho in data array
  ipr = ndim + ndimV + 2		!  pressure 
  iutherm = ndim + ndimV + 3	!  thermal energy
  ih = ndim + ndimV + 4		!  smoothing length
  ipmass = ndim + ndimV + 5	!  particle mass      

  label(ix(1:ndim)) = labelcoord(1:ndim)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i)
  enddo
  label(irho) = '\gr'
  label(ipr) = 'P      '
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'      

  label(ndim + ndimV+6) = '\ga'
  if (ncolumns.gt.ndim+ndimV+6) then
     iBfirst = ndim + ndimV+6+1	! location of Bx
     iBlast = ndim + ndimV+6+ndimV	! location of Bz      
     do i=1,ndimV
        label(ndim + ndimV+6+i) = 'B\d'//labelcoord(i) !' (x10\u-3\d)'	!//'/rho'
     enddo
     idivB = ndim+ndimV+ndimV+7	 
     label(idivB) = 'div B'
     do i=1,ndimV
        label(ndim + ndimV+ndimV+7 + i) = 'J'//labelcoord(i)
     enddo
  else	 
     iBfirst = 0
     iBlast = 0
  endif
  !--these are here for backwards compatibility -- could be removed
  if (ncolumns.gt.ndim+3*ndimV+7) then
     label(ndim + 3*ndimV+8) = 'v_parallel'
     label(ndim + 3*ndimV+9) = 'v_perp'
     label(ndim + 3*ndimV+10) = 'B_parallel'
     label(ndim + 3*ndimV+11) = 'B_perp'
  endif


!-----------------------------------------------------------

  return    
!
!--errors
!
77 continue
  print*,' *** Error encountered while reading file ***'
  print*,' -> Check that magnetic field is toggled correctly'
  return

78 continue
  print*,' *** Error encountered while reading timestep ***'
  print*,' -> number of columns in data file not equal to'
  print*,'    that set as a parameter - edit and recompile'
  return

79 continue
  print*,' *** Error reading data file header: check format ***'
  return

80 continue
  print*,' *** data file empty, no steps read ***'
  return

81 continue
  print*,' *** Error: can''t open data file ***'
  return

end subroutine read_data
