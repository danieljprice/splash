!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! THIS VERSION FOR DAN'S SPMHD CODE (BINARY DUMPS)
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

subroutine read_data(rootname,istart,nfilesteps)
  use exact
  use particle_data
  use params
  use labels
  use settings
  implicit none
  integer, intent(IN) :: istart
  integer, intent(OUT) :: nfilesteps
  character(LEN=*), intent(IN) :: rootname
  character(LEN=LEN(rootname)+4) :: datfile
  character(LEN=2) :: fileno
  integer :: i,j,k,ifile,icol,ipos
  integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
  integer :: npartin,ntotin
  logical :: iexist,reallocate
  real*8 :: timein,gammain,hfactin
  real*8, dimension(:), allocatable :: dattemp
  real*8, dimension(:,:), allocatable :: dattempvec

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
     !print*,'rootname = ',rootname
  else
     print*,' **** no data read **** ' 
     return
  endif

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  k=1      

50 continue
  ivegotdata = .false.
  !
  !--open data file and read data
  !
  open(unit=11,ERR=81,file=datfile,status='old',form='unformatted')

  nfilesteps = 100000
!
!--read first header line
!
  read(11,ERR=78,end=80) timein,npartin,ntotin,gammain, &
       hfactin,ndim_max,ndimV_max,ncol_max	 
!  print*,'reading time = ',timein,npartin,ntotin,gammain, &
!       ndim_max,ndimV_max,ncol_max     
!
!--allocate memory for data arrays (initially for 11 timesteps)
!
  if (ntotin.lt.5500) then
     nstep_max = max(111,maxstep)
  elseif (ntotin.lt.111111) then
     nstep_max = max(11,maxstep)
  else 
     nstep_max = max(5,maxstep)
  endif
  npart_max = max(ntotin,maxpart)
  if (.not.allocated(dat) .or. ntotin.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncol_max.gt.maxcol) then
     call alloc(npart_max,nstep_max,ncol_max)
  endif
!
!--rewind file
!
  rewind(11)

  do i=istart,nfilesteps
     !!print*,' reading step ',i
     reallocate = .false.
     npart_max = maxpart
     nstep_max = maxstep
     ncol_max = maxcol
     !
     !--read header line for this timestep
     !
     read(11,ERR=67,end=67) timein,npart(i),ntot(i),gammain, &
          hfactin,ndim,ndimV,ncolumns	 
     time(i) = real(timein)
     gamma(i) = real(gammain)
     hfact = real(hfactin)
     npartoftype(1,i) = npart(i)
     npartoftype(2,i) = ntot(i)-npart(i)
     print*,'reading time = ',time(i),npart(i),ntot(i),gamma(i), &
          hfact,ndim,ndimV,ncolumns
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
	npart_max = int(1.1*ntot(i))
     endif
     if (i.eq.maxstep) then
        nstep_max = i + max(10,INT(0.1*nstep_max))
        reallocate = .true.
     endif
     !
     !--reallocate memory for main data array
     !
     if (reallocate) then
        call alloc(npart_max,nstep_max,ncol_max)
     endif

  
     if (ntot(i).gt.0) then
        !
	!--read position vector
        !
	icol = 1
	if (allocated(dattempvec)) deallocate(dattempvec)
	allocate(dattempvec(3,ntot(i)))
	read (11,end=66,ERR=67) dattempvec(1:ndim,1:ntot(i))
	do ipos = 1,ndim
	   dat(1:ntot(i),ipos,i) = real(dattempvec(ipos,1:ntot(i)))
	enddo
	icol = icol + ndim
	!
        !--read velocity vector
        !
	read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))	
	do ipos = icol,icol+ndimV-1
	   dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
	enddo
	icol = icol + ndimV
	!
        !--read scalar variables
        !
	if (allocated(dattemp)) deallocate(dattemp)
        allocate(dattemp(ntot(i)))
	do j = 1,5
           read (11,end=66,ERR=67) dattemp(1:ntot(i))
	   dat(1:ntot(i),icol,i) = real(dattemp(1:ntot(i)))
	   icol = icol + 1
	enddo
	!
        !--read alpha
        !
	read (11,end=66,ERR=67) dattempvec(1:3,1:ntot(i))
	do ipos = icol,icol+2
	   dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
	enddo
	icol = icol + 3
	!
        !--Bfield
        !
	read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))
	do ipos = icol, icol+ndimV-1
	   dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
	enddo
	icol = icol + ndimV
	!
	!--curl B
	!
        read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))
	do ipos = icol, icol+ndimV-1
	   dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
        enddo
	icol = icol + ndimV
	!
        !--div B
	!
        read (11,end=66,ERR=67) dattemp(1:ntot(i))
	dat(1:ntot(i),icol,i) = real(dattemp(1:ntot(i)))
        icol = icol + 1
	!
        !--psi
	!
        read (11,end=66,ERR=67) dattemp(1:ntot(i))
	dat(1:ntot(i),icol,i) = real(dattemp(1:ntot(i)))
        icol = icol + 1	
	!
	!--force
	!
        read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))
	do ipos = icol, icol+ndimV-1
	   dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
	enddo
	icol = icol+ndimV-1
	!!print*,'columns read = ',icol,' should be = ',ncolumns

	deallocate(dattemp,dattempvec)
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

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use settings
  implicit none
  integer :: i

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

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  label(irho) = '\gr'
  label(ipr) = 'P      '
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'      

  label(ndim + ndimV+6) = '\ga'
  label(ndim + ndimV+7) = '\ga\du'
  if (ncolumns.gt.ndim+ndimV+7) then
     label(ndim + ndimV+8) = '\ga\dB'
     iBfirst = ndim + ndimV+8+1	! location of Bx
     iBlast = ndim + ndimV+8+ndimV	! location of Bz      
     do i=1,ndimV
        label(ndim + ndimV+8+i) = 'B\d'//labelcoord(i,1) !' (x10\u-3\d)'	!//'/rho'
     enddo
     do i=1,ndimV
        label(ndim + ndimV+ndimV+8 + i) = 'J'//labelcoord(i,1)
     enddo
     iJfirst = ndim+2*ndimV+8+1
     idivB = ndim+3*ndimV+9 
     label(idivB) = 'div B'
  else	 
     iBfirst = 0
     iBlast = 0
  endif
  if (ncolumns.gt.ndim+3*ndimV+9) then
     label(ndim+3*ndimV+10) = 'psi'
  endif
   if (ncolumns.gt.ndim+3*ndimV+10) then
     label(ndim+3*ndimV+11) = 'f_visc_x'
     label(ndim+3*ndimV+12) = 'f_visc_y'
     label(ndim+3*ndimV+13) = 'f_x'
     label(ndim+3*ndimV+14) = 'f_y'
  endif 
  
  !--these are here for backwards compatibility -- could be removed
!  if (ncolumns.gt.ndim+3*ndimV+7) then
!     label(ndim + 3*ndimV+8) = 'v_parallel'
!     label(ndim + 3*ndimV+9) = 'v_perp'
!     label(ndim + 3*ndimV+10) = 'B_parallel'
!     label(ndim + 3*ndimV+11) = 'B_perp'
!  endif

 !
 !--set labels for each type of particles
 !
 ntypes = 2
 labeltype(1) = 'gas'
 labeltype(2) = 'ghost'

!-----------------------------------------------------------

  return 
end subroutine set_labels
