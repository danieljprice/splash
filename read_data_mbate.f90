!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING FORMATTED OUTPUT FROM MATTHEW BATE'S CODE
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

subroutine read_data(rootname,istart,nfilesteps)
  use particle_data
  use params
  use labels
  use settings
  implicit none
  integer, intent(IN) :: istart
  integer, intent(OUT) :: nfilesteps
  character(LEN=2) :: fileno
  character(LEN=*), intent(IN) :: rootname
  integer :: i,j,k,ifile
  integer :: ncol_max,nstep_max
  integer :: ntotin
  logical :: iexist
  real :: gammain,timeff
  
  character(LEN=LEN(rootname)+10), dimension(1000) :: filename,sinkfile
  character(LEN=LEN(rootname)+10) :: imagefile
  integer :: nsinkcolumns, int_from_string
  logical :: magfield

  print*,'entering read_data'
!
!--assume MHD if filename starts with m
!
  magfield = .false.
  if (rootname(1:1).EQ.'m') magfield = .true.
!
!--fix number of spatial dimensions
!
  ndim = 3
  ndimV = 3
  if (magfield) then
     ncolumns = 18	! number of columns in file
     nsinkcolumns = 10
  else 
     ncolumns = 11
     nsinkcolumns = 10
  endif
  ivegotdata = .false.
  iexist = .false.
  
  print*,'rootname'
  ifile = int_from_string(rootname(6:7))
  print*,' image file starting at number ',ifile
  
  write(fileno,"(i1,i1)") ifile/10,mod(ifile,10)
  !!  fileno = achar(48+ifile/10)//achar(48+mod(ifile,10))
  imagefile = 'IP'//rootname(1:5)//fileno
  print *,' opening ',imagefile
  k=1
  !
  !--allocate memory initially
  !
  ncol_max = 1
  ntotin = 1
  nstep_max = 100
  call alloc(ntotin,nstep_max,ncol_max)
  
50 continue
  i = k
  open(unit=15,file=imagefile,status='old',form='formatted',ERR=81)
     read(15,*,end=55) gammain
  !      print*,'Number of particles = ',ntot
     print*,'gamma = ',gammain
     do i=k,100000
	read(15,*, end=55,ERR=79) timeff,time(i),ntot(i),nghost(i),filename(i)
	!!         filename(i) = '../'//filename(i)
	print*,filename(i)
	gamma(i) = gammain
	!
	!--reallocate if exceeded max timesteps
	!
	if (i.eq.maxstep) then
	   nstep_max = maxstep + 10
	   call alloc(ntotin,nstep_max,ncol_max)
	endif
     enddo
55   continue
     print*,'end of image file, nsteps=',i-1
     ifile = ifile + 1
     write(fileno,"(i1,i1)") ifile/10,mod(ifile,10)
     !!fileno = achar(48+ifile/10)//achar(48+mod(ifile,10))
     imagefile = 'IP'//rootname(1:5)//fileno
     inquire (file=imagefile, exist=iexist)
     print*,imagefile,' exist = ',iexist
     k = i
     if (iexist) then
	print*,' opening ',imagefile
	goto 50
     endif
      
     nfilesteps = k-1
56   continue
  close(15)

  ntotin = maxval(ntot)
  nstep_max = nfilesteps
  ncol_max = ncolumns
  !
  !--allocate memory for main data array before reading data
  !
  if (.not.allocated(dat) .or. ntotin.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncol_max.gt.maxcol) then
     call alloc(ntotin,nstep_max,ncol_max)
  endif
  !
  !--now read data
  !
  do i=1,nfilesteps
     if (ntot(i).gt.maxpart) print*,'ntot > array limits!!'      
     !         READ (11,*,END=66) time(i)
     print*,'t = ',time(i), '  file = ',trim(filename(i))
     !
     !--read data from QG file (gas particles)
     !
     open(unit=11,file=trim(filename(i)),status='old',form='formatted')
        read (11,*, end=66, ERR=77) (dat(1:ncolumns,j,i),iam(j,i),j=1,maxpart)
66      continue
  
        ntot(i) = j-1
	if (ntot(i)-nghost(i).gt.0) then
	   npart(i) = ntot(i) - nghost(i)	! assumes always more ghosts
	else					! than particles
	   npart(i) = ntot(i)
	endif
	print*,'Number of particles,ghosts = ', &
	     npart(i),nghost(i),':',ntot(i)-npart(i),' ghosts output'
     close(unit=11)
     !
     !--read data from QS file (sink particles)
     !
     sinkfile(i) = filename(i)(1:1)//'S'//filename(i)(3:)
     inquire (file=sinkfile(i), exist=iexist)

     if (iexist) then
	print*,'reading sink file ',sinkfile(i)
	open(unit=12,file=sinkfile(i),status='old',form='formatted')
	read(12,*,end=68,ERR=67) (dat(1:nsinkcolumns,j,i),j=ntot(i)+1,maxpart)
67      continue
        print*,' Error reading sinkfile - no sinks read'	    
68      continue	    
        print*,' sinks = ',j-1 - ntot(i)
        do k = ntot(i),j-1
           iam(k,i) = 1
        enddo
        ntot(i) = j-1
        print*,'ntotal = ',ntot(i) 
     else
	print*,'sink file not found'                
     endif

  enddo      

  print*,'>> READ all steps =',i-1,'ntot = ',ntot(i-1),'nghost=',ntot(i-1)-npart(i-1)
      
  ivegotdata = .true.
  
  !!------------------------------------------------------------
  !! set labels for each column of data
  
  do i=1,ndim
     ix(i) = i
  enddo
  ivx = ndim + 1
  ivlast = ndim + ndimV
  irho = ndim + ndimV + 1		! location of rho in data array
  iutherm = ndim + ndimV + 2	!  thermal energy
  ih = ndim + ndimV + 3		!  smoothing length
  ipmass = ndim + ndimV + 4	!  particle mass      
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  label(irho) = '\gr'      
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'      
  
  label(ndim + ndimV+5) = '\ga'
  if (magfield) then
     iBfirst = ndim + ndimV+5+1	! location of Bx
     iBlast = ndim + ndimV+5+ndimV	! location of Bz      
     do i=1,ndimV
	label(ndim + ndimV+5+i) = 'B\d'//labelcoord(i,1) !' (x10\u-3\d)'	!//'/rho'
     enddo
     idivB = ndim + ndimV+ndimV+6
     label(ndim + ndimV+ndimV+6) = 'div B'
     do i=1,ndimV
	label(ndim + ndimV+ndimV+6 + i) = 'J'//labelcoord(i,1)
     enddo
  else	 
     iBfirst = 0
     iBlast = 0
  endif
  !      label(ndim + ndimV+2*ndimV+8) = 'v_parallel'
  !      label(ndim + ndimV+2*ndimV+9) = 'v_perp'
  !      label(ndim + ndimV+2*ndimV+10) = 'B_parallel'
  !      label(ndim + ndimV+2*ndimV+11) = 'B_perp'
  !
  !--specify which of the possible quantities you would like to calculate
  !  (0 = not calculated)
  ncalc = 9	! specify number to calculate
  ipr = ncolumns + 1
  ientrop = ncolumns + 2      
  irad = ncolumns + 3
  irad2 = 0
  ipmag = ncolumns + 4
  ibeta = ncolumns + 5
  itotpr = ncolumns + 6      
  ike = ncolumns + 7
  idivBerr = ncolumns + 8
  itimestep = ncolumns + 9
  if (ipr.NE.0) label(ipr) = 'P      '
  
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
print*,' *** Error reading IP file: check format ***'
return

80 continue
print*,' *** data file empty, no steps read ***'
return

81 continue
print*,' *** Error: can''t open IP file ***'
return
                    
end subroutine read_data
