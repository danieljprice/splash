!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE (VERSION 2.0)
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
! dat(maxpart,maxplot,maxstep) : main data array
!
! npart(maxstep)      : number of particles in each timestep
! ntot(maxstep)       : total number of particles in each timestep
! nghost(maxstep)     : number of ghost particles in each timestep
! iam(maxpart,maxstep): integer identification of particle type
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,istart,nfilesteps)
  use particle_data
  use params
  use labels
  use settings
  implicit none
  integer, intent(IN) :: istart
  integer, intent(OUT) :: nfilesteps
  character(LEN=*), intent(IN) :: rootname
  character(LEN=LEN(rootname)+10) :: datfile
  character(LEN=2) :: fileno
  integer*4, dimension(6) :: Npartoftypetemp
  integer, dimension(:), allocatable :: iamtemp
  integer :: i,j,itype,icol,ifile,idashpos,ioftype
  integer :: index1,index2,indexstart,indexend,Nmassesdumped
  integer :: ncol_max,npart_max,nstep_max,int_from_string
  logical :: iexist,reallocate
  real*8 :: timetemp
  real*8, dimension(6) :: Massoftype
  real*4, dimension(:), allocatable :: dattemp1
  real*4, dimension(:,:), allocatable :: dattemp

  if (rootname(1:1).ne.' ') then
     idashpos = index(rootname,'_') ! position of '_' in string
     if (idashpos.eq.0) then
        !
        !--for rootnames without the '_000', read all files starting at #1
        !
	ifile = 1
        !--work out the first filename
        write(datfile,"(a,'_','00',i1)") trim(rootname),ifile
     else
        !
	!--otherwise just read this dump
	!
        ifile = int_from_string(rootname(idashpos+1:idashpos+3))
	datfile = trim(rootname)
     endif
  else
     print*,' **** no data read **** ' 
     return
  endif
!
!--check if first data file exists
!
  ivegotdata = .false.  
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     print*,' *** error: ',trim(datfile),' file not found ***'    
     return
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim = 3
  ndimV = 3
  ncol_max = 12 ! 3 x pos, 3 x vel, utherm, rho, Ne, h, pmass
!
!--read data from snapshots
!  
  i = istart
  !
  !--allocate memory for data arrays (initially for 11 timesteps)
  !
  if (i.eq.1) then  ! on first step, allocate for several timesteps
     npart_max = 1
     nstep_max = 1
     ncol_max = 12
     if (.not.allocated(dat)) then
        call alloc(npart_max,nstep_max,ncol_max)
     endif
  endif

  do while (iexist)
     write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
     !
     !--open data file and read data
     !
     open(11,ERR=81,file=datfile,status='old',form='unformatted')
     !
     !--read header for this timestep
     !
     read(11,ERR=70,end=80) Npartoftype(:,i),Massoftype,timetemp 
     print*,'Npartoftype =',Npartoftype(:,i)
     ntot(i) = int(sum(Npartoftype))
     npart(i) = int(Npartoftype(1,i))
     print*,'Npartoftype =',Npartoftype(:,i)
     print*,'Massoftype = ',Massoftype
     time(i) = real(timetemp)
     print*,'t = ',time(i),' npart, ntot = ',npart(i),ntot(i)

     reallocate = .false.
     npart_max = maxpart
     nstep_max = maxstep
     ncol_max = maxcol

     if (ntot(i).gt.maxpart) then
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
        if (allocated(dattemp)) deallocate(dattemp)
        allocate(dattemp(3,ntot(i)))
        !
        !--read positions of all particles
        !
	print*,'positions ',ntot(i)
        read (11, end=66) dattemp(1:3,1:ntot(i))	
	do icol=1,3
	   dat(1:ntot(i),icol,i) = dattemp(icol,1:ntot(i))
	enddo
	!
        !--same for velocities
        !
	print*,'velocities ',ntot(i)
        read (11, end=66,ERR=72) dattemp(1:3,1:ntot(i))		
	do icol=4,6
	   dat(1:ntot(i),icol,i) = dattemp(icol-3,1:ntot(i))
	enddo
	!
        !--read particle ID
	!
	print*,'particle ID ',ntot(i)
	if (allocated(iamtemp)) deallocate(iamtemp)
	allocate(iamtemp(npart_max))
        read (11, end=66,ERR=73) iamtemp(1:ntot(i))
	iam(1:ntot(i),i) = int(iamtemp(1:ntot(i)))
	deallocate(iamtemp)
        !
        !--read particle masses
        !
        !--work out total number of masses dumped 
	Nmassesdumped = 0
	do itype = 1,6
	   if (abs(Massoftype(itype).lt.1.e-8)) then
	      Nmassesdumped = Nmassesdumped + Npartoftype(itype,i)
	   endif
	enddo
	print*,'particle masses ',Nmassesdumped

	!--read this number of entries
        if (allocated(dattemp1)) deallocate(dattemp1)
        allocate(dattemp1(Nmassesdumped))
	read(11,end=66,err=74) dattemp1(1:Nmassesdumped)
	!--now copy to the appropriate sections of the .dat array
	indexstart = 1
	index1 = 1
	
	do itype=1,6
	   if (Npartoftype(itype,i).ne.0) then
	      index2 = index1 + Npartoftype(itype,i) -1
	      if (abs(Massoftype(itype)).lt.1.e-8) then ! masses dumped
	         indexend = indexstart + Npartoftype(itype,i) - 1
	         print*,'read ',Npartoftype(itype,i),' masses for type ', &
		        itype,index1,'->',index2,indexstart,'->',indexend
	         dat(index1:index2,7,i) = dattemp1(indexstart:indexend)
	      else  ! masses not dumped
	         print*,'setting masses for type ',itype,' = ', &
		        real(Massoftype(itype)),index1,'->',index2
	         dat(index1:index2,7,i) = real(Massoftype(itype))
	      endif
	      index1 = index2 + 1
	      indexstart = indexend + 1
	   endif
	enddo
	deallocate(dattemp1)
	!
        !--read other quantities for rest of particles
        !
	print*,'gas properties ',npart(i)
        do icol=8,12
	   !!print*,icol
	   read (11, end=66,ERR=78) dat(1:npart(i),icol,i)
	enddo
	
     else
        ntot(i) = 1
        npart(i) = 1
        dat(:,:,i) = 0.
     endif
     iam(:,i) = 0
     !
     !--if just the rootname has been input, 
     !  set next filename and see if it exists
     !
     ifile = ifile + 1
     i = i + 1
     if (idashpos.eq.0) then
        if (ifile.lt.10) then
           write(datfile,"(a,'_','00',i1)") trim(rootname),ifile
        elseif (ifile.lt.100) then
           write(datfile,"(a,'_','0',i2)") trim(rootname),ifile     
        elseif (ifile.lt.1000) then
           write(datfile,"(a,'_',i3)") trim(rootname),ifile
        else
           print*,'error: ifile > 1000 in filename'
	   return
        endif
        inquire(file=datfile,exist=iexist)
     else
        iexist = .false. ! exit loop
     endif
  enddo

  nfilesteps = i-1		! this is if reached array limits
  !!ntot(i-1) = j-1
!
!--now memory has been allocated, set arrays which are constant for all time
!
  nghost = 0
  gamma = 5./3.
  goto 68

66 continue
  print*,'*** end of file reached in ',trim(datfile),' ***'
  nfilesteps = i		! timestep there but data incomplete
  ntot(i) = j-1
  nghost(i) = 0.
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
  print*,'ncolumns = ',ncolumns

  print*,'>> Finished reading: steps =',nfilesteps,'last step ntot = ',ntot(nfilesteps)

  !!------------------------------------------------------------
  !! set locations of particular items of data used in the plotting program

  do i=1,ndim
     ix(i) = i
  enddo
  ivx = 4
  ivlast = 6
  ipmass = 7
  irho = 9	! location of rho in data array
  ipr = 0
  iutherm = 8	        !  thermal energy
  ih = 12		!  smoothing length
  !
  !--set labels of the quantities read in
  !
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = '\gr'
  label(iutherm) = 'u'
  label(10) = 'Ne'
  label(11) = 'N\dH'
  label(ih) = 'h'
  label(ipmass) = 'particle mass'
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  
  !--set labels for each particle type
  !
  ntypes = maxparttypes
  labeltype(1) = 'gas'
  labeltype(2) = 'dark matter'
  labeltype(5) = 'star'

!-----------------------------------------------------------

  return    
!
!--errors
!
70 continue
  print*,' *** Error encountered while reading timestep header ***'
  print*,' Npartoftype = ',Npartoftype(:,i)
  print*,' Massoftype = ',Massoftype
  return

71 continue
  print*,' *** Error encountered while reading positions ***'
  return

72 continue
  print*,' *** Error encountered while reading velocities ***'
  return

73 continue
  print*,' *** Error encountered while reading particle ID ***'
  return

74 continue
  print*,' *** Error encountered while reading particle masses ***'
  return

75 continue
  print*,' *** Error encountered while reading star particle masses ***'
  return

78 continue
  print*,' *** Error encountered while reading gas particle properties ***'
  return

80 continue
  print*,' *** data file empty, no steps read ***'
  return

81 continue
  print*,' *** Error: can''t open data file ***'
  return

end subroutine read_data
