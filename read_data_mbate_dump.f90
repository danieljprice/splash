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
    
  character(LEN=20), dimension(1000) :: filename,sinkfile
  character(LEN=20) :: dumpfile
  integer :: nsinkcolumns, int_from_string
  integer :: nprint, nghosti, n1, n2, rhozero, RK2
  logical :: magfield
  real*8, dimension(:,:), allocatable :: dattemp
  real*8 :: udisti,umassi,utimei, umagfdi, timei, gammai
  real*8 :: escap,tkin,tgrav,tterm,tmag

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
     ncolumns = 19	! number of columns in file
     nsinkcolumns = 10
  else 
     ncolumns = 11
     nsinkcolumns = 10
  endif
  ivegotdata = .false.
  iexist = .false.
  
  print*,'rootname'
  ifile = int_from_string(rootname(6:7))
  print*,' dumps starting at number ',ifile
  
  write(fileno,"(i1,i1)") ifile/10,mod(ifile,10)
  dumpfile = rootname(1:5)//fileno
  print *,' opening ',dumpfile
  k=1
  !
  !--allocate memory initially
  !
  ncol_max = 19
  ntotin = 30000
  nstep_max = 100
  nprint = 1
  call alloc(ntotin,nstep_max,ncol_max)
 !
 !--open the (unformatted) binary file and read the number of particles
 !
50 continue
  j = k
  open(unit=15,file=dumpfile,status='old',form='unformatted',ERR=81)
!
!--loop over the timesteps in this file
!     
     do j=k,100000
     
        ntotin = max(ntotin,nprint)
!
!--read the initial information to determine number of particles
!
!        read(15,end=55) udisti, umassi, utimei, umagfdi,  &
!            nprint, nghosti, n1, n2, gt, gammai, rhozero, RK2
!
!--allocate/reallocate memory now that we know the number of particles
!
	print*,ntotin
	if (ntotin.gt.maxpart) then
	   call alloc(ntotin,nstep_max,ncol_max)
	endif
!
!--allocate a temporary array for double precision variables
!
        if (allocated(dattemp)) deallocate(dattemp)
	allocate(dattemp(ncol_max,ntotin))
!
!--now read the timestep data in the dumpfile
!

        read(15,end=55) udisti, umassi, utimei, umagfdi,  &
          nprint, nghosti, n1, n2, timei, gammai, rhozero, RK2, &
	  (dattemp(7,i), i=1, nprint), (dattemp(8,i), i=1,nprint), &
          escap, tkin, tgrav, tterm, tmag, &
          (dattemp(1,i), i=1, nprint), (dattemp(2,i), i=1, nprint), &
          (dattemp(3,i), i=1, nprint), (dattemp(4,i), i=1, nprint), &
          (dattemp(5,i), i=1, nprint), (dattemp(6,i), i=1, nprint), &
          (dattemp(9,i), i=1, nprint), (dattemp(10,i), i=1, nprint), &
          (dattemp(11,i), i=1, nprint), (dattemp(12,i), i=1, nprint), &  
          (dattemp(13,i), i=1, nprint), (dattemp(14,i), i=1, nprint), &
          (dattemp(15,i), i=1, nprint), (dattemp(16,i), i=1, nprint), &
          (dattemp(17,i), i=1, nprint), (dattemp(18,i), i=1, nprint), &
          (dattemp(19,i), i=1, nprint)
!
!--convert to single precision
!     
          print*,'converting to single precision ',ncol_max
	  print*,'x,y(1) = ',dattemp(1,1),dattemp(2,1)
          dat(1:ncol_max,1:nprint,j) = real(dattemp(1:ncol_max,1:nprint))
          deallocate(dattemp)

	  ntot(j) = nprint
	  nghost(j) = nghosti
	  npart(j) = ntot(j) - nghost(j)
	  gamma(j) = real(gammai)
	  time(j) = real(timei)
          if (ntotin.eq.30000) ntotin = nprint
     enddo

55   continue
     print*,'end of dump file, nsteps=',j-1
     ifile = ifile + 1
     write(fileno,"(i1,i1)") ifile/10,mod(ifile,10)
     dumpfile = rootname(1:5)//fileno
     inquire (file=dumpfile, exist=iexist)
     print*,dumpfile,' exist = ',iexist
     k = j
     if (iexist) then
	print*,' opening ',dumpfile
	goto 50
     endif
      
     nfilesteps = k-1
56   continue
  close(15)

  print*,'>> READ all steps =',j-1,'ntot = ',ntot(j-1),'nghost=',ntot(j-1)-npart(j-1)
      
  ivegotdata = .true.
  
  !!------------------------------------------------------------
  !! set labels for each column of data
  
  do i=1,ndim
     ix(i) = i
  enddo
  ivx = 4
  ivlast = 6
  irho = 18		! location of rho in data array
  iutherm = 16	!  thermal energy
  ih = 7		!  smoothing length
  ipmass = 17	!  particle mass      
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  label(irho) = '\gr'      
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'     
  label(8) = 'alpha'
  label(19) = 'psi' 
  
  label(ndim + ndimV+5) = '\ga'
  if (magfield) then
     iBfirst = 9	! location of Bx
     iBlast = 11	! location of Bz      
     do i=1,ndimV
	label(iBfirst + i-1) = 'B\d'//labelcoord(i,1) !' (x10\u-3\d)'	!//'/rho'
     enddo
     idivB = 12
     label(idivB) = 'div B'
     do i=1,ndimV
	label(13 + i-1) = 'J'//labelcoord(i,1)
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
print*,' *** Error: can''t open dump file ***'
return
                    
end subroutine read_data
