!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!-------------------------------------------------------------------------
     
      SUBROUTINE read_data(rootname,ifile,gamma,time,
     &  dat1,iam,ncolumns,ncalc,nfilesteps,ntot,npart,
     &  nghost,ndim,ndimV,hfact,ivegotdata)
      USE params
      USE labels
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ncolumns,ncalc,nfilesteps
      INTEGER, INTENT(OUT) :: ndim,ndimV
      INTEGER, INTENT(OUT), DIMENSION(max) :: iam
      INTEGER, INTENT(OUT), DIMENSION(maxstep) :: npart,ntot,nghost
      REAL, INTENT(OUT), DIMENSION(maxplot,max,maxstep) :: dat1
      REAL, INTENT(OUT), DIMENSION(maxstep) :: time,gamma
      REAL, INTENT(OUT) :: hfact
      REAL, DIMENSION(maxstep) :: timeff
      REAL :: gamma_temp
      CHARACTER*20, DIMENSION(maxstep) :: filename,sinkfile
      CHARACTER*20 :: imagefile,datfile
      CHARACTER fileno*2
      CHARACTER(LEN=*), INTENT(IN) :: rootname
      INTEGER i,j,k,ifile,ifile10,ifile1
      INTEGER :: nsinkcolumns
      LOGICAL :: iexist,ivegotdata,magfield      

!
!--assume MHD if filename starts with m
!
      magfield = .false.
      IF (rootname(1:1).EQ.'m') magfield = .true.
!
!--fix number of spatial dimensions
!      
      ndim = 3
      ndimV = 3
      IF (magfield) THEN
         ncolumns = 18	! number of columns in file
	 nsinkcolumns = 10
      ELSE 
         ncolumns = 11
	 nsinkcolumns = 10
      ENDIF	 
      ivegotdata = .false.
      
      PRINT*,'rootname'
      ifile10 = IACHAR(rootname(6:6))-48
      ifile1 = IACHAR(rootname(7:7))-48
      PRINT*,ifile10,ifile1
      ifile = 10*ifile10 + ifile1
      
      PRINT*,' image file starting at number ',ifile
      fileno = ACHAR(48+ifile/10)//ACHAR(48+mod(ifile,10))
      imagefile = 'IP'//rootname(1:5)//fileno
      PRINT *,' opening ',imagefile
      k=1
      
50    CONTINUE
      i = k
      OPEN(unit=15,file=imagefile,status='old',form='formatted')
      READ(15,*,END=55) gamma_temp
c      PRINT*,'Number of particles = ',ntot
      PRINT*,'gamma = ',gamma_temp
      DO i=k,maxstep
         READ(15,*, END=55,ERR=79) time(i),timeff(i),
     &                      ntot(i),nghost(i),filename(i)
!         filename(i) = '../'//filename(i)
	 PRINT*,filename(i)
	 gamma(i) = gamma_temp
      ENDDO
55    CONTINUE
      PRINT*,'end of image file, nsteps=',i-1
      ifile = ifile + 1
      fileno = ACHAR(48+ifile/10)//ACHAR(48+mod(ifile,10))
      imagefile = 'IP'//rootname(1:5)//fileno
      INQUIRE (file=imagefile, exist=iexist)
      k = i
      IF (iexist) THEN
         PRINT*,' opening ',imagefile
         GOTO 50
      ENDIF
      
      nfilesteps = k-1      
56    CONTINUE
      CLOSE(15)
            
      DO i=1,nfilesteps
         IF (ntot(i).gt.max) PRINT*,'ntot greater than array limits!!'      
!         READ (11,*,END=66) time(i)
	 PRINT*,'t = ',time(i), '  file = ',filename(i)
!
!--read data from QG file (gas particles)
!
         OPEN(unit=11,file=filename(i),status='old',form='formatted')
         READ (11,*, END=66, ERR=77)
     &	      (dat1(1:ncolumns,j,i),iam(j),j=1,max)
66       CONTINUE

	 ntot(i) = j-1
	 IF (ntot(i)-nghost(i).gt.0) THEN
	    npart(i) = ntot(i) - nghost(i)	! assumes always more ghosts
	 ELSE					! than particles
	    npart(i) = ntot(i)
	 ENDIF
	    PRINT*,'Number of particles,ghosts = ',
     &        npart(i),nghost(i),':',ntot(i)-npart(i),' ghosts output'
	 CLOSE(unit=11)
!
!--read data from QS file (sink particles)
!
         sinkfile(i) = filename(i)(1:1)//'S'//filename(i)(3:)
	 INQUIRE (file=sinkfile(i), exist=iexist)
	 IF (iexist) THEN
            PRINT*,'reading sink file ',sinkfile(i)
	    OPEN(unit=12,file=sinkfile(i),status='old',form='formatted')
	    READ(12,*,END=68,ERR=67)
     &           (dat1(1:nsinkcolumns,j,i),j=ntot(i)+1,max)
67          CONTINUE
            PRINT*,' Error reading sinkfile - no sinks read'	    
68          CONTINUE	    
            PRINT*,' sinks = ',j-1 - ntot(i)
            DO k = ntot(i),j-1
               iam(k) = 1
            ENDDO
            ntot(i) = j-1
	    PRINT*,'ntotal = ',ntot(i) 
	 ELSE
	    PRINT*,'sink file not found'                
	 ENDIF

      ENDDO      

      PRINT*,'>> READ all steps =',i-1,'ntot = ',ntot(i-1),
     &          'nghost=',ntot(i-1)-npart(i-1)
      
      ivegotdata = .true.

!!------------------------------------------------------------
!! set labels for each column of data

      DO i=1,ndim
         ix(i) = i
      ENDDO
      ivx = ndim + 1
      ivlast = ndim + ndimV
      irho = ndim + ndimV + 1		! location of rho in data array
      iutherm = ndim + ndimV + 2	!  thermal energy
      ih = ndim + ndimV + 3		!  smoothing length
      ipmass = ndim + ndimV + 4	!  particle mass      

      label(ix(1:ndim)) = labelcoord(1:ndim)
      DO i=1,ndimV
         label(ivx+i-1) = 'v\d'//labelcoord(i)
      ENDDO
      label(irho) = '\gr'      
      label(iutherm) = 'u'
      label(ih) = 'h       '
      label(ipmass) = 'particle mass'      
      
      label(ndim + ndimV+5) = '\ga'
      IF (magfield) THEN
	 iBfirst = ndim + ndimV+5+1	! location of Bx
	 iBlast = ndim + ndimV+5+ndimV	! location of Bz      
         DO i=1,ndimV
	    label(ndim + ndimV+5+i) = 'B\d'//labelcoord(i) !' (x10\u-3\d)'	!//'/rho'
	 ENDDO
	 idivB = ndim + ndimV+ndimV+6
	 label(ndim + ndimV+ndimV+6) = 'div B'
	 DO i=1,ndimV
	    label(ndim + ndimV+ndimV+6 + i) = 'J'//labelcoord(i)
	 ENDDO
      ELSE	 
         iBfirst = 0
	 iBlast = 0
      ENDIF
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
      IF (ipr.NE.0) label(ipr) = 'P      '

!-----------------------------------------------------------

      RETURN    
!
!--errors
!
77    CONTINUE
      PRINT*,' *** Error encountered while reading file ***'
      PRINT*,' -> Check that magnetic field is toggled correctly'
      RETURN      
      
78    CONTINUE
      PRINT*,' *** Error encountered while reading timestep ***'
      PRINT*,' -> number of columns in data file not equal to'
      PRINT*,'    that set as a parameter - edit and recompile'
      RETURN

79    CONTINUE
      PRINT*,' *** Error reading data file header: check format ***'
      RETURN

80    CONTINUE
      PRINT*,' *** data file empty, no steps read ***'
      RETURN

81    CONTINUE
      PRINT*,' *** Error: can''t open data file ***'
      RETURN
                    
      END SUBROUTINE read_data
