!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!-------------------------------------------------------------------------
     
      SUBROUTINE read_data(rootname,ifile,gamma,time,
     &  dat1,iam,ncolumns,ncalc,nfilesteps,ntot,npart,nghost,
     &  ndim,ndimV,hfact,ivegotdata)
      USE params
      USE labels
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ncolumns,ncalc,nfilesteps
      INTEGER, INTENT(OUT) :: ndim,ndimV
      INTEGER, INTENT(OUT), DIMENSION(max) :: iam
      INTEGER, INTENT(OUT), DIMENSION(maxstep) :: npart,ntot,nghost
      REAL, INTENT(OUT), DIMENSION(maxplot,max,maxstep) :: dat1
      REAL, INTENT(OUT), DIMENSION(maxstep) :: time, gamma
      REAL, INTENT(OUT) :: hfact
      CHARACTER*20 imagefile,filename(maxstep),datfile
      CHARACTER fileno*2
      CHARACTER(LEN=*), INTENT(IN) :: rootname
      INTEGER i,j,k,ifile
      INTEGER ncol_max,ndim_max,ndimV_max
      LOGICAL :: iexist,ivegotdata,magfield
      
!
!--assume MHD if filename starts with m
!
      magfield = .false.
      IF (rootname(1:1).NE.' ') THEN
!
!--if rootname does not contain .dat, make it end in .dat
!
	 IF (INDEX(rootname,'.dat').EQ.0) THEN
           datfile = rootname(1:LEN_TRIM(rootname))//'.dat'
	 ELSE
	   datfile = rootname(1:LEN_TRIM(rootname))  
         ENDIF
         IF (rootname(1:1).EQ.'m') magfield = .true.	 
         ifile = 1
         PRINT*,'rootname = ',rootname, ' mhd = ',magfield
      ELSE
         PRINT*,' **** no data read **** ' 
	 RETURN
      ENDIF

      PRINT *,' opening ',datfile
      k=1      
      
50    CONTINUE
!
!--open data file and read data
!
      ivegotdata = .false.

      OPEN(unit=11,ERR=81,file=datfile,status='old',form='formatted')
!      READ(11,*,ERR=79,END=80) npart(1),gamma,hfact
      
      nfilesteps = maxstep-1
      
      DO i=1,nfilesteps
!
!--read header line for this timestep
!
         READ(11,*,ERR=78,END=67) time(i),npart(i),ntot(i),gamma(i),
     &                            hfact,ndim,ndimV,ncolumns	 
         PRINT*,'reading time = ',time(i),npart(i),ntot(i),gamma(i),
     &	                          ndim,ndimV,ncolumns     
         IF (ncolumns.GT.ncol_max) ncol_max = ncolumns
	 IF (ncolumns.NE.ncol_max) THEN
	    PRINT*,'*** Warning number of columns not equal for timesteps'
	 ENDIF
	 IF (ndim.GT.ndim_max) ndim_max = ndim
	 IF (ndimV.GT.ndimV_max) ndimV_max = ndimV   
!
!--allocate memory for main data array here
!     
c	 PRINT*,'data columns = ',nplot
	 IF (ntot(i).GT.max) STOP 'ntot greater than array limits!!'      
         IF (ntot(i).GT.0) THEN
	    READ (11,*, END=66,ERR=77) (dat1(1:ncolumns,j,i),j=1,ntot(i))
         ELSE
	    ntot(i) = 1
	    npart(i) = 1
	    nghost(i) = 0
	    dat1(:,:,i) = 0.
	 ENDIF
      ENDDO 
      
      PRINT*,' REACHED ARRAY LIMITS IN READFILE'

66    CONTINUE
      nfilesteps = i		! timestep there but data incomplete
      ntot(i) = j-1
      nghost(i) = ntot(i) - npart(i)
      GOTO 68
        
67    CONTINUE
      nfilesteps = i-1		! no timestep there at all
     
68    CONTINUE
!
!--close data file and return
!      	      
      CLOSE(unit=11)
      
      ivegotdata = .true.
      ncolumns = ncol_max
      ndim = ndim_max
      ndimV = ndimV_max
      PRINT*,'ncolumns = ',ncolumns

      PRINT*,'>> READ all steps =',nfilesteps,'last step ntot = ',
     &    ntot(nfilesteps)

!!------------------------------------------------------------
!! set labels for each column of data

      DO i=1,ndim
         ix(i) = i
      ENDDO
      ivx = ndim + 1
      ivlast = ndim + ndimV
      irho = ndim + ndimV + 1		! location of rho in data array
      ipr = ndim + ndimV + 2		!  pressure 
      iutherm = ndim + ndimV + 3	!  thermal energy
      ih = ndim + ndimV + 4		!  smoothing length
      ipmass = ndim + ndimV + 5	!  particle mass      

      label(ix(1:ndim)) = labelcoord(1:ndim)
      DO i=1,ndimV
         label(ivx+i-1) = 'v\d'//labelcoord(i)
      ENDDO
      label(irho) = '\gr'
      label(ipr) = 'P      '
      label(iutherm) = 'u'
      label(ih) = 'h       '
      label(ipmass) = 'particle mass'      
      
      label(ndim + ndimV+6) = '\ga'
      IF (magfield) THEN
	 iBfirst = ndim + ndimV+6+1	! location of Bx
	 iBlast = ndim + ndimV+6+ndimV	! location of Bz      
         DO i=1,ndimV
	    label(ndim + ndimV+6+i) = 'B\d'//labelcoord(i) !' (x10\u-3\d)'	!//'/rho'
	 ENDDO
	 idivB = ndim+ndimV+ndimV+7	 
	 label(idivB) = 'div B'
	 DO i=1,ndimV
	    label(ndim + ndimV+ndimV+7 + i) = 'J'//labelcoord(i)
	 ENDDO
      ELSE	 
         iBfirst = 0
	 iBlast = 0
      ENDIF
      label(ndim + 3*ndimV+8) = 'v_parallel'
      label(ndim + 3*ndimV+9) = 'v_perp'
      label(ndim + 3*ndimV+10) = 'B_parallel'
      label(ndim + 3*ndimV+11) = 'B_perp'
!
!--specify which of the possible quantities you would like to calculate
!  (0 = not calculated)
      ncalc = 8	! specify number to calculate
      ientrop = ncolumns + 1      
      irad = ncolumns + 2
      irad2 = ncolumns + 3
      ipmag = ncolumns + 4
      ibeta = ncolumns + 5
      itotpr = ncolumns + 6      
      ike = ncolumns + 7
      idivBerr = ncolumns + 8
      itimestep = 0

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
