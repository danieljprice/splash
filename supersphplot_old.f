      PROGRAM supersphplot
!--------------------------------------------------------------------------
!     plots me a graph or two from sph output
!
!     file format is specified in the subroutine read_data   
!     also plots density renderings and vector maps
!     for 3D data does cross-sections as render plots
!
!     uses PGPLOT routines to plot graphs
!
!     Written by: Daniel Price, Institute of Astronomy, Cambridge UK
!          email: dprice@ast.cam.ac.uk
!
!     This version for both NDSPMHD and Matthew Bate's code 2003
!     Changes log:
!      18/07/03 - transformations (log, 1/x etc)
!      15/07/03 - interpolate2D,3D - much simpler than smooth_pixels
!      25/06/03 - clever Makefile - makes for dansph or M. Bate SPH
!               - subroutines in different files
!      20/06/03 - rendering can handle zero density, prints if min=max
!      19/06/03 - multiplot with rendering, no x array
!      16/06/03 - labels in module, specified in read_data
!      09/06/03 - ndim, ndimV changeable, reads ncolumns in data file
!      26/05/03 - calculated quantities no longer in read_data
!      16/05/03 - colour schemes + rendering using SPH summation
!                 output format has changed, reads pmass array also
!      13/05/03 - read polytrope in menu option
!      01/05/03 - can manually enter plot limits 
!      29/04/03 - bug in initial plot limits fixed
!--------------------------------------------------------------------------
      
      USE params
      USE labels      
      IMPLICIT NONE      
      REAL, PARAMETER :: pi = 3.1415926536
      INTEGER :: ndim,ndimV
      INTEGER nstep,i,j,k,numplot,ncolumns,len,ipolyc
      INTEGER ipickx,ipicky,iplotx,iploty
      INTEGER ilist,ihoc,ibin,listsize
      INTEGER ncolours,nstart,n_end,nfilesteps,nfreq,ifile
      INTEGER, DIMENSION(maxstep) :: npart,ntot,ntotplot
      INTEGER, DIMENSION(max) :: iam
      INTEGER ntotmin,nacross,ndown
      INTEGER ncalc,nextra,ndataplots,ihalf,iadjust
      INTEGER imark, imarkg, npart1
      INTEGER ixsec,nxsec, nbins,nc,icircpart
      INTEGER nyplot,nyplots,nyplotmulti
      INTEGER linestylein, iexact, icolours,ncontours      
      INTEGER irender,ivecplot,irenderplot,npix,npixvec
      INTEGER ivecplot_old,npix_old,npixvec_old,ncontours_old
      INTEGER, DIMENSION(maxplot) :: multiploty,multiplotx
      INTEGER, DIMENSION(maxplot) :: irendermulti,ivecplotmulti      
      INTEGER, DIMENSION(maxplot) :: npixmulti,npixvecmulti 
      INTEGER, DIMENSION(maxplot) :: ncontoursmulti,itrans
      CHARACTER(LEN=8) string	! used in PGPLOT calls
      REAL, DIMENSION(maxplot,max,maxstep) :: dat
      REAL, DIMENSION(2,max) :: vecplot
      REAL, DIMENSION(maxstep) :: time
      REAL, DIMENSION(maxplot,2) :: lim
      REAL, DIMENSION(max) :: xplot,yplot,renderplot
      REAL, DIMENSION(:,:), ALLOCATABLE :: datpix
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: datpix3D
      REAL :: xmin,xmax,ymin,ymax,zmin,zmax,xminrender,xmaxrender
      REAL :: Bmin,Bmax,vmin,vmax,Bmag,rendermin,rendermax
      REAL gamma,hfact,temp
      REAL scale,rhomax,rhomin,binsize
      REAL scalemax
      REAL zoom,diff,mid
      REAL xsecmin,xsecmax,dxsec
      REAL charheight
      REAL angledeg,anglexy,runit(ndimmax)	! to plot r at some angle
!--exact solution parameters
!--toy star
      INTEGER :: norder
      REAL :: Htstar,Atstar,Ctstar,totmass,sigma,sigma0 ! toy star parameters
!--sound wave      
      REAL delta,lambda
!--polytrope
      REAL mtot,maxrho,akfac,den(ipolycmax),rad(ipolycmax)
      
      CHARACTER(LEN=20) :: rootname,filename
      CHARACTER(LEN=60) :: title,titlex,datfile
      CHARACTER(LEN=1) :: ans,dummy
      CHARACTER(LEN=20) :: labelx,labely,labelz,labelrender
      CHARACTER(LEN=25) :: transform_label
!--plot options
      LOGICAL :: animate,iadapt,magfield,done,iexist,ihavereadfilename
      LOGICAL :: plotcirc,plotcircall,x_sec,flythru,imulti,ipagechange
      LOGICAL :: iplotline,iplotlinein,iplotav,ilabelpart
      LOGICAL :: iplotpart,iplotpartvec,iplotcont,iplotghost
      LOGICAL :: ishowopts, ivegotdata,isamexaxis
      LOGICAL, DIMENSION(maxplot) :: iplotcontmulti, iplotpartvecmulti
      LOGICAL :: iplotcont_old,iplotpartvec_old
!
!-----------------------------------------------------------------------------
!
      PRINT*,' Welcome to Dan''s ND supersphplot 2002...version 3.4 '
!
!--set default options
!
      rootname = 'blank'
      numplot=maxplot 	! reset if read from file
      ncalc = 3		! number of columns to calculate(e.g. radius)
      nextra = 0	! extra plots aside from particle data
      iextra = 0
      ncolumns=maxplot-ncalc
      ndim = ndimmax
      ndimV = ndim	! default velocity same dim as coords
      len=7
      ncolours=10
      nstart = 1
      n_end = maxstep        
      nfilesteps = maxstep
      nfreq = 1
      animate = .true.
      magfield = .true.
      iadapt = .true.
      done = .false.
      plotcirc = .false.
      plotcircall = .false.
      icircpart = 1
      x_sec = .false.
      flythru = .false.
      ipagechange = .true.
      xmin=0.0
      xmax=-1.
      scalemax = 1.0
      zoom = 1.0
      imark = 17
      imarkg = 4
      nacross = 1
      ndown = 1
      charheight = 1.0
      iplotline = .false.
      iplotlinein = .false.
      linestylein = 4
      iexact = 0
      iplotav = .false.
      nbins = 24
      ilabelpart = .false.
      iplotpart = .true.
      iplotpartvec = .true.
      npix = 20
      npixvec = 20
      ivecplot = 0
      irender = 0
      irenderplot = 0
      iplotcont = .true.
      ncontours = 30
      iplotghost = .true.
      icolours = 0
      lambda = 1.0	! sound wave exact solution
      delta = 0.005	! sound wave exact solution
      ishowopts = .false.
      ivegotdata = .false.	! flag when data has been read in
      labelcoord(1) = 'x'
      labelcoord(2) = 'y'
      labelcoord(3) = 'z'
      isamexaxis = .true.

      nyplots = 1	! number of plots each
      nyplotmulti = nyplots
      multiploty(1) = ndim+1 
      multiplotx(1) = 1
      irendermulti(1) = 0
      ivecplotmulti(1) = 0
      itrans(:) = 0

! ------------------------------------------
! read default options if file exists
!
      INQUIRE (exist=iexist, file='defaults')
      IF (iexist) THEN
         OPEN(unit=1,file='defaults',status='old',form='formatted')
	    READ(1,*,END=7,ERR=8) animate,magfield,iadapt,x_sec,
     &                             flythru,plotcirc
	    READ(1,*,END=7,ERR=8) imark, imarkg, nacross, ndown,nyplotmulti     
	    READ(1,*,END=7,ERR=8) iplotline,iplotlinein,linestylein
	    READ(1,*,END=7,ERR=8) iexact, iplotav, nbins
	    READ(1,*,END=7,ERR=8) irender,ivecplot,iplotpartvec,npix,npixvec
	    READ(1,*,END=7,ERR=8) iplotcont,ncontours,icolours,iplotghost
	    READ(1,*,END=7,ERR=8) delta,lambda
	    READ(1,*,END=7,ERR=8) Htstar,Atstar,Ctstar,sigma0,norder

            READ(1,*,END=7,ERR=8) itrans(:)
	    READ(1,*,END=7,ERR=8) multiplotx(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) multiploty(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) irendermulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) iplotcontmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) ncontoursmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) ivecplotmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) npixmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) npixvecmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) iplotpartvecmulti(1:nyplotmulti)	 
	 CLOSE(unit=1)
	 PRINT*,'read default options from file '
      ENDIF
      GOTO 9
7     CONTINUE
      PRINT*,'**** Warning: end of file in defaults ****'
      CLOSE(unit=1)
      GOTO 9
8     CONTINUE
      PRINT*,'Error reading defaults from file'
      CLOSE(unit=1)     
9     CONTINUE

c ------------------------------------------
c prompt for title

       title = ' '
!       print*,' Enter title for graphs '
!       read(*,808) title
808    FORMAT (a40)       
       print*,'title = ',title

c ------------------------------------------
c get rootname from command line/file

      CALL getarg(1, rootname)
      IF (rootname(1:1).NE.' ') ihavereadfilename = .true.
      
! ----------------------------------------------------------------
! menu - loop back to here once finished plotting/setting options
!
200   CONTINUE
!
!--numplot is the total number of data columns (read + calculated)
!   not including the particle co-ordinates
!  nextra are extra graphs to plot (e.g. convergence plots)
!
      IF (iexact.EQ.4) THEN	! toy star plot A-C plane
         nextra = 1
	 iextra = ncolumns + ncalc + 1
	 label(iextra) = 'A-C plane'
	 IF (magfield) THEN
	    sigma = sigma0
	 ELSE
	    sigma = 0.
	 ENDIF
      ELSE		! otherwise no extra plots
         nextra = 0	 
      ENDIF	 

      numplot = ncolumns + ncalc + nextra
      ndataplots = ncolumns + ncalc
      
      IF ((nacross.GT.1).OR.(ndown.GT.1)) THEN
         title = '          '
      ENDIF
!
!--these are the quantities calculated by the program
!      
      IF (ientrop.NE.0) label(ientrop) = 'entropy'
      IF (irad.NE.0) label(irad) = 'radius '
      IF (irad2.NE.0) label(irad2) = 'r_parallel'
      IF (ike.NE.0) label(ike) = 'specific KE'
      IF (ipr.NE.0) label(ipr) = 'P_gas '
      IF (ipmag.NE.0) label(ipmag) = 'P_mag'
      IF (itotpr.NE.0) label(itotpr) = 'P_gas + P_mag'
      IF (ibeta.NE.0) label(ibeta) = 'Plasma \gb'
      IF (idivBerr.NE.0) label(idivBerr) = 'h |div B| / |B|'
!
!--if first time through, read in timesteps from file(s)
!   
      IF (ihavereadfilename) THEN
	 ipicky = numplot+2
         GOTO 25
      ENDIF	 
      
!----------------------------------------------------------------------
!--print menu
!     
      CALL print_menu(ipicky,ipickx)
       
      PRINT*,' You may choose from a delectable sample of plots '
      PRINT 12
      ihalf = numplot/2		! print in two columns
      iadjust = MOD(numplot,2)
      PRINT 11, (i,transform_label(label(i),itrans(i)),
     &             ihalf + i + iadjust,
     &             transform_label(label(ihalf + i +
     &              iadjust),itrans(ihalf+i+iadjust)),i=1,ihalf)
      IF (iadjust.NE.0) THEN
         PRINT 13, ihalf + iadjust,transform_label(label(ihalf +
     &	 iadjust),itrans(ihalf+iadjust))
      ENDIF
      PRINT 12
      PRINT 18,numplot+1,
     &        'multiplot ',
     &         multiploty(1:nyplotmulti)
!
!--supersphplot options
!     
      IF (ishowopts) THEN
     
      PRINT 12
      PRINT 13,numplot+2,'read dat'
      PRINT 14,numplot+3,
     &        'Change number of timesteps read  ',(n_end-nstart+1)/nfreq 
      PRINT 15,numplot+4,
     &        'toggle animate                   ',animate
      PRINT 15,numplot+5,
     &        'toggle adaptive/fixed limits     ',iadapt
      PRINT 15,numplot+6,
     &        'toggle mag field                 ',magfield
      PRINT 15,numplot+7,
     &        'toggle cross section/projection  ',x_sec
      PRINT 15,numplot+8,
     &        'toggle circles of interaction    ',plotcirc
      PRINT 16,numplot+9,
     &        'Change graph markers             ',imark,imarkg
      PRINT 16,numplot+10,
     &        'Change plots per page            ',nacross,ndown
      PRINT 14,numplot+11,
     &        'Set multiplot                    ',nyplotmulti
      PRINT 15,numplot+12,
     &        'toggle page change               ',ipagechange
      PRINT 19,numplot+13,
     &        'toggle plot line (all, init)     ',iplotline,iplotlinein 
      PRINT 14,numplot+14,
     &        'toggle exact solution            ',iexact
      PRINT 20,numplot+15,
     &        'toggle plot average line         ',iplotav,nbins
      PRINT 15,numplot+16,
     &        'toggle label particles           ',ilabelpart
      PRINT 15,numplot+17,
     &        'toggle plot ghosts               ',iplotghost
      PRINT 16,numplot+18,
     &        'set rendering/vector plots       ',irender,ivecplot
      PRINT 13,numplot+19,
     &        'apply transformations (log10,1/x)'
      IF (iadapt) THEN
         PRINT 17,numplot+20,
     &        'Zoom out/in                      ',scalemax      
      ELSE
         PRINT 17,numplot+20,
     &        'Zoom out/in/set manual limits    ',zoom
      ENDIF
      
      ENDIF	! show/hide opts
      
      PRINT 12
      IF (ishowopts) THEN
         PRINT 13,numplot+21,'hide options'      
      ELSE
         PRINT 13,numplot+21,'supersphplot options'
      ENDIF
      PRINT 13,numplot+22,'Save defaults'
      PRINT 13,numplot+23,'Exit supersphplot'
      PRINT 12
9901  CONTINUE
      PRINT*,'Please enter your selection now (y axis): '
11    FORMAT(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12    FORMAT(1x,45('-'))
13    FORMAT(1x,i2,')',1x,a)
14    FORMAT(1x,i2,')',1x,a,'( ',i2, ' )')
15    FORMAT(1x,i2,')',1x,a,'( ',L1,' )')
16    FORMAT(1x,i2,')',1x,a,'( ',i2,',',i2,' )')
17    FORMAT(1x,i2,')',1x,a,'( ',f4.2,' )')
18    FORMAT(1x,i2,')',1x,a,'( ',10(i2,1x),' )')
19    FORMAT(1x,i2,')',1x,a,'( ',L1,',',1x,L1,' )')
20    FORMAT(1x,i2,')',1x,a,'( ',L1,',',i2,' )')
!
!--read y axis selection and if needed x axis selection
!
      READ (*,*, ERR=9901) ipicky

25    CONTINUE
      IF ((ipicky.le.(numplot-nextra)).and.(ipicky.gt.0)) THEN
         PRINT*,' (x axis):'
	 READ (*,*,ERR=25) ipickx
	 IF ((ipickx.gt.numplot) .or. ipickx.lt.1) GOTO 200
      ENDIF   
      
!-----------------------------------------------------------------------
! menu actions
!      
      imulti = .false.
      isamexaxis = .true.	! used to determine whether to plot labels
      
      IF ((ipicky.gt.numplot+22).or.(ipicky.lt.1)) THEN
         STOP
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+1) THEN
          imulti=.true.
	  IF (ANY(multiplotx(1:nyplotmulti).NE.multiplotx(1))) THEN
	     isamexaxis = .false.
	  ENDIF
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+2) THEN
          IF (.not.ihavereadfilename) THEN
	     PRINT*,'Enter filename to read:' 	!(5 characters):'
	     READ(*,30) rootname
30           FORMAT(a)	   	  
	  ENDIF
	  ihavereadfilename = .false.
	  CALL read_data(rootname,ifile,gamma,time(:),
     &         dat(:,:,:),iam(:),
     &         ncolumns,ncalc,nfilesteps,ntot(:),npart(:),
     &         ndim,ndimV,hfact,ivegotdata)
!
!--calculate various quantities
!     
          PRINT*,'calculating radius etc'
	  numplot = ncolumns + ncalc
	  DO i=nstart,n_end
	     DO j=1,npart(i)
!!--pressure if not in data array
	        IF ((ipr.GT.ncolumns)
     &		   .AND.(irho.NE.0).AND.(iutherm.NE.0))
     &		   dat(ipr,j,i) = dat(irho,j,i)*dat(iutherm,j,i)*(gamma-1.)
!!--entropy
                IF (ientrop.NE.0) THEN
  	           IF (dat(irho,j,i).GT.1.e-10) THEN
                      dat(ientrop,j,i) = dat(ipr,j,i)
     &		                        /dat(irho,j,i)**gamma
                   ELSE  
		      dat(ientrop,j,i) = 0.
		   ENDIF		   
		ENDIF   
!!--radius	 
                IF (irad.NE.0)    
     &	           dat(irad,j,i) = SQRT(DOT_PRODUCT(dat(ix(1:ndim),j,i),
     &                                             dat(ix(1:ndim),j,i)))
!!--distance along the line with angle 30deg w.r.t. x axis
                IF (irad2.NE.0) THEN
                   angledeg = 30.
		   anglexy = angledeg*pi/180.
		   runit(1) = COS(anglexy)
		   IF (ndim.GE.2) runit(2) = SIN(anglexy)
	           dat(irad2,j,i) = DOT_PRODUCT(dat(ix(1:ndim),j,i),
     &		                                runit(1:ndim))
		ENDIF
!!--specific KE
                IF ((ike.NE.0).AND.(ivx.NE.0)) 
     &             dat(ike,j,i) = 0.5*DOT_PRODUCT(dat(ivx:ivlast,j,i),
     &		                                  dat(ivx:ivlast,j,i))
!!--magnetic pressure
                IF ((ipmag.NE.0).AND.(iBfirst.NE.0))
     &             dat(ipmag,j,i) = 
     &                 0.5*DOT_PRODUCT(dat(iBfirst:iBlast,j,i),
     &		                       dat(iBfirst:iBlast,j,i))
!!--plasma beta
                IF ((ibeta.NE.0).AND.(ipmag.NE.0)) THEN
		   IF (abs(dat(ipmag,j,i)).GT.1e-10) THEN
                      dat(ibeta,j,i) = dat(ipr,j,i)/dat(ipmag,j,i)
		   ELSE  
		      dat(ibeta,j,i) = 0.
		   ENDIF
		ENDIF   
!!--total pressure (gas + magnetic)     
                IF ((itotpr.NE.0).AND.(ipr.NE.0).AND.(ipmag.NE.0)) THEN
		   dat(itotpr,j,i) = dat(ipr,j,i) + dat(ipmag,j,i)
		ENDIF
!!--div B error	(h*divB / abs(B))	
		IF ((idivB.NE.0).AND.(ih.NE.0).AND.(iBfirst.NE.0)) THEN
		   Bmag = SQRT(DOT_PRODUCT(dat(iBfirst:iBlast,j,i),
     &		                           dat(iBfirst:iBlast,j,i)))
                   IF (Bmag.GT.1e-10) THEN
		     dat(idivBerr,j,i) = 
     &		         abs(dat(idivB,j,i))*dat(ih,j,i)/Bmag
		   ELSE
		     PRINT*,'Bmag < 1e-10 !!'
		     dat(idivBerr,j,i) = 0.
		   ENDIF
		ENDIF
	     ENDDO
	  ENDDO
	  PRINT*,'Setting plot limits...'
	  n_end = nfilesteps
	  ntotplot(:) = npart(:)
	  IF (iplotghost) ntotplot(:) = ntot(:)
	  zoom = 1.0
!!--find limits of particle properties	  
	  lim(:,1) = 1.e6
	  lim(:,2) = -1.e6
	  DO j=1,numplot
 	     DO i=nstart,n_end
	        lim(j,1) = MIN(lim(j,1),MINVAL(dat(j,1:ntotplot(i),i)))
	        temp = MAXVAL(dat(j,1:ntotplot(i),i))*scalemax
		IF (lim(j,2).LT.temp) lim(j,2) = temp
	     ENDDO
	     IF (lim(j,2).EQ.lim(j,1)) THEN
	        PRINT*,label(j),' min = max = ',lim(j,1)
	     ENDIF
	  ENDDO
!!--limits of magnetic field (in all dimensions) for vector plot
	  IF (iBfirst.NE.0) THEN
	     IF (iadapt) THEN
	        Bmin = 0.0
		Bmax = maxval(dat(iBfirst:iBlast,1:ntotplot(1),1))*scalemax
	     ELSE
	        Bmax = 0.0
		Bmin = 0.0
	        DO j=iBfirst,iBlast
	          !Bmin = min(Bmin,lim(j,1))
		  IF (lim(j,2).gt.Bmax) THEN
		     Bmax=lim(j,2)
		     PRINT*,' Bmax = ',label(j),lim(j,2)
	          ENDIF
		ENDDO
	     ENDIF
	     PRINT*,'Bmin,Bmax = ',Bmin,Bmax
	  ENDIF
	  PRINT*,'plot limits set'
!
!--read toy star file for toy star solution
!
          IF (iexact.EQ.4) THEN
	     filename = TRIM(rootname)//'.tstar'
	     OPEN(UNIT=20,ERR=8801,FILE=filename,STATUS='old')
	      READ(20,*,ERR=8801) Htstar,Ctstar,Atstar
	      READ(20,*,ERR=8801) sigma0
	      READ(20,*,ERR=8801) norder
	     CLOSE(UNIT=20)
	     PRINT*,' >> read ',filename
	     PRINT*,' H,C,A,sigma,n = ',Htstar,Ctstar,Atstar,sigma0,norder
	     GOTO 200
8801         CONTINUE
             PRINT*,'Cannot open/ error reading ',filename	     
	  ENDIF
	  GOTO 200
!------------------------------------------------------------------------	  
      ELSEIF (ipicky.eq.numplot+3) THEN
          PRINT*,' Start, end at timestep numbers:'
	  READ*,nstart,n_end
	  IF (n_end.gt.nfilesteps) THEN
	     PRINT*,'n_end greater than nfilesteps, reset to ',nfilesteps
	     n_end = nfilesteps
	  ENDIF
	  PRINT*,' Frequency of steps to READ:'
	  READ*,nfreq
          PRINT *,' Steps = ',(n_end-nstart+1)/nfreq
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+4) THEN
          animate=.not.animate
          PRINT *,' Animation = ',animate
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+5) THEN
          iadapt=.not.iadapt
          PRINT *,' Adaptive plot limits = ',iadapt
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+6) THEN
          magfield=.not.magfield
          PRINT *,' Mag field = ',magfield
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+7) THEN
          x_sec =.not.x_sec
	  PRINT *,' Cross section = ',x_sec
	  flythru = .false.
	  IF (x_sec) THEN
	     PRINT*,'Do you want a fly-through (y/n)?'
	     READ*,ans
	     IF (ans.eq.'y'.or.ans.eq.'Y') flythru=.true.
	  ENDIF
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+8) THEN
	  plotcirc=.not.plotcirc
          PRINT*,' Plot circles of interaction = ',plotcirc
	  IF (plotcirc) THEN
	     PRINT*,'Plot all circles (y/n)?'
	     READ*,ans
	     plotcircall = .true.
	     IF (ans.ne.'y'.AND.ans.ne.'Y') THEN
	        plotcircall = .false.
	        PRINT*,'Enter particle number to plot circle around'
	        READ*,icircpart
	     ENDIF
	  ENDIF
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+9) THEN
          PRINT*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
	  PRINT*,' Enter PGPLOT marker # (particles):'
	  READ*,imark
	  PRINT*,' Enter PGPLOT marker # (ghosts):'
	  READ*,imarkg
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+10) THEN
	  PRINT*,'Enter number of plots across:'
	  READ*,nacross
	  PRINT*,'Enter number of plots down'
	  READ*,ndown
	  GOTO 200	    	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+11) THEN
	  PRINT*,'Enter number of plots per timestep:'
	  READ*,nyplotmulti
          IF (nyplotmulti.GT.numplot) THEN 
	     nyplotmulti = numplot
	     PRINT*,'nyplots > numplot, reset to ',numplot
	  ENDIF
	  nacross = nyplotmulti/2
	  IF (nacross.EQ.0) nacross = 1
	  ndown = nyplotmulti/nacross
	  PRINT*,'setting nacross,ndown = ',nacross,ndown  
          PRINT*,'Same x axis for all (y/n)?'
	  READ*,ans
	  IF (ans.EQ.'y'.OR.ans.EQ.'Y') THEN
 	     PRINT*,'Enter x axis for all plots '
	     READ*,multiplotx(1)
	     multiplotx(2:nyplotmulti) = multiplotx(1)	     
	  ENDIF
	  DO i=1,nyplotmulti
	     PRINT*,'Enter y axis plot number ',i	     
	     READ*,multiploty(i)
	     IF ((ans.NE.'y'.AND.ans.NE.'Y')
     &           .AND.multiploty(i).le.ndataplots) THEN
	        PRINT*,'Enter x axis plot number ',i
  	        READ*,multiplotx(i)
	     ENDIF
	     IF ((multiplotx(i).LE.ndim)
     &	         .AND.(multiploty(i).LE.ndim)) THEN
                CALL get_render_options(irendermulti(i),
     &		     npixmulti(i),icolours,iplotcontmulti(i),
     &	             ncontoursmulti(i),ivecplotmulti(i),npixvecmulti(i),
     &               iplotpartvecmulti(i),ndim,numplot)
             ENDIF	     
	  ENDDO
	  GOTO 200	    	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+12) THEN
	  ipagechange=.not.ipagechange
          PRINT*,' Page changing = ',ipagechange
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+13) THEN
	  PRINT*,' Plot initial only(i), all(a), both(b) or not (n)?'
	  READ*,ans
	  iplotline = .false.
	  iplotlinein = .false.
	  IF (ans.EQ.'i'.OR.ans.EQ.'b') iplotlinein = .true.
	  IF (ans.EQ.'a'.OR.ans.EQ.'b') iplotline = .true.
	  IF (iplotlinein) THEN
	     PRINT*,' Enter PGPLOT line style (0-5)'
	     READ*,linestylein
	  ENDIF
          PRINT*,' Plot line = ',iplotline,iplotlinein
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+14) THEN
	  IF (iexact.NE.0) THEN
	     iexact = 0
	  ELSE   
40           WRITE (*,41)
41           FORMAT(' Enter exact solution to plot ',/,
     &              ' 0) none ',/,
     &              ' 1) polytrope ',/,
     &              ' 2) soundwave ',/,     
     &              ' 3) sedov blast wave',/,     
     &              ' 4) toy star ')
             READ*,iexact
	     IF ((iexact.LT.0).OR.(iexact.GT.4)) GOTO 40
	     PRINT*,' plotting exact solution number ',iexact
!
!--enter parameters for soundwave solution
!
	     IF (iexact.EQ.4) THEN
	        PRINT*,' Toy star: enter parameters H, A, C'
		READ*,Htstar,Atstar,Ctstar
		sigma = 0.
		IF (magfield) THEN
		   PRINT*,' Enter const for By proportional to rho '
		   READ*,sigma0
		   sigma = sigma0
		ENDIF
		PRINT*,'Do you want oscillations?'
		READ*,ans
		norder = -1
		IF (ans.EQ.'y'.OR.ans.EQ.'Y') THEN
		   PRINT*,'Enter order:'
		   READ*,norder
		ENDIF
	     ELSEIF (iexact.EQ.2) THEN
	        PRINT*,' Enter wavelength of sound wave lambda '
		READ*,lambda
		PRINT*,' Enter amplitude '
		READ*,delta
	     ELSEIF (iexact.EQ.1) THEN
!
!--read exact solution for a polytrope
!
                ipolyc = ipolycmax
		INQUIRE (exist = iexist, file='polycalc.dat')
		IF (.not.iexist) THEN
                   PRINT*,' file polycalc.dat does not exist'
                ENDIF
                OPEN(unit=14,file='polycalc.dat',	
     &		     status='old',form='formatted')
                READ(14,10) dummy
10              FORMAT(a)      
                READ(14,*) maxrho,mtot
                READ(14,*) akfac
                READ(14,*, END=100) (den(i),rad(i),i=1,ipolyc) 
                CLOSE(14)
                GOTO 101
100             CONTINUE
                PRINT*,'End of polycalc.dat, i=',i-1
                ipolyc = i-1
                CLOSE(14)
101             CONTINUE	     
	     ENDIF
	  ENDIF
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+15) THEN
	  iplotav=.not.iplotav
          IF (iplotav) THEN
51	     PRINT*,' Enter number of bins for averaging '
	     READ*,nbins
	     IF ((nbins.LE.0).OR.(nbins.GT.1000)) GOTO 51
          ENDIF
	  PRINT*,' Plot average, nbins = ',iplotav,nbins
	  GOTO 200 		  	    	  	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+16) THEN
c	  label particles with particle numbers
	  ilabelpart=.not.ilabelpart
          PRINT*,' label particles = ',ilabelpart
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+17) THEN
c	  plot ghost particles?
	  iplotghost=.not.iplotghost
          PRINT*,' plot ghost particles = ',iplotghost
	  IF (iplotghost) ntotplot(:) = ntot(:)
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+18) THEN      
          CALL get_render_options(irender,npix,icolours,iplotcont,
     &	        ncontours,ivecplot,npixvec,iplotpartvec,ndim,numplot)
          GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+19) THEN
299       CONTINUE
	  PRINT*,'Enter plot number to apply transformation (0 to finish)'
	  READ*,ipicky
	  IF (ipicky.LE.numplot .AND. ipicky.GT.0) THEN
	     WRITE(*,*) (i,') ',transform_label('x',i), i=1,5)
             PRINT*,'Enter transformation to apply:'
             READ*,itrans(ipicky)
	     GOTO 299
	  ELSE
	     GOTO 200
	  ENDIF
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+20) THEN
	  IF (.not.iadapt) THEN
	     PRINT*,'Do you want to manually enter limits?'
	     READ*,ans
	     IF (ans.EQ.'y' .OR. ans.EQ.'Y') THEN
301	        PRINT*,'Enter plot number to set limits (0 to finish)'
		READ*,ipicky
		IF ((ipicky.LE.numplot).AND.(ipicky.GT.0)) THEN
		   PRINT*,' Current limits   (min,max):',
     &		          lim(ipicky,1),lim(ipicky,2)
		   PRINT*,' Enter new limits (min,max):'
		   READ*,lim(ipicky,1),lim(ipicky,2)
		   PRINT*,' >> limits set    (min,max):',
     &		         lim(ipicky,1),lim(ipicky,2)
		ELSE
		   PRINT*,'Finished, returning'
		   GOTO 200
		ENDIF
		GOTO 301
	     ELSE
	        PRINT*,'Enter zoom factor for fixed limits'
	        READ*,zoom
	        DO i=1,numplot
	           diff = lim(i,2)- lim(i,1)
	           mid = 0.5*(lim(i,1) + lim(i,2))
		   lim(i,1) = mid - 0.5*zoom*diff
		   lim(i,2) = mid + 0.5*zoom*diff
	        ENDDO
	     ENDIF	
	  ELSE
	     PRINT*,'Enter scale factor (adaptive limits)'
	     READ*,scalemax
	  ENDIF
	  GOTO 200
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+21) THEN
!	  show/hide plot options
	  ishowopts = .not.ishowopts
	  GOTO 200 	  
!------------------------------------------------------------------------
      ELSEIF (ipicky.eq.numplot+22) THEN
          OPEN(unit=1,file='defaults',status='replace',
     &         form='formatted')
	      WRITE(1,*) animate,magfield,iadapt,x_sec,flythru,plotcirc
	      WRITE(1,*) imark,imarkg,nacross,ndown,nyplotmulti
	      WRITE(1,*) iplotline,iplotlinein,linestylein
	      WRITE(1,*) iexact,iplotav,nbins
	      WRITE(1,*) irender,ivecplot,iplotpartvec,npix,npixvec
	      WRITE(1,*) iplotcont,ncontours,icolours,iplotghost
	      WRITE(1,*) delta,lambda
	      WRITE(1,*) Htstar,Atstar,Ctstar,sigma0,norder
	      
	      WRITE(1,*) itrans(:)
	      WRITE(1,*) multiplotx(1:nyplotmulti)
	      WRITE(1,*) multiploty(1:nyplotmulti)
	      WRITE(1,*) irendermulti(1:nyplotmulti)
	      WRITE(1,*) iplotcontmulti(1:nyplotmulti)
	      WRITE(1,*) ncontoursmulti(1:nyplotmulti)
	      WRITE(1,*) ivecplotmulti(1:nyplotmulti)
	      WRITE(1,*) npixmulti(1:nyplotmulti)
	      WRITE(1,*) npixvecmulti(1:nyplotmulti)
	      WRITE(1,*) iplotpartvecmulti(1:nyplotmulti)    
	  CLOSE(unit=1)
	  PRINT*,'default options saved to file'
	  GOTO 200
      ENDIF

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! if data has not been read from a file, prompt for file
      IF (.not.ivegotdata) THEN
         PRINT*,' No data '
	 ihavereadfilename = .false.
	 ipicky = numplot+2
         GOTO 25
      ENDIF
        
!------------------------------------------------------------------------
! initialise PGPLOT
!------------------------------------------------------------------------

!!--set current plot to first in multiplot array if doing multiplot
      nxsec = 1
      IF (imulti) THEN 
         ipickx = multiplotx(1)
	 ipicky = multiploty(1)
	 nyplots = nyplotmulti
      ELSE
         nyplots = 1 
      ENDIF	 
      
!------------------------------------------------------------------------
!!--co-ordinate plot initialisation

      IF (ipicky.le.ndim .and. ipickx.le.ndim) THEN

!!--if cross-section read limits for slice      
        IF (x_sec.and.flythru) THEN
	   DO j=1,3
	      IF ((j.ne.ipickx).and.(j.ne.ipicky)) ixsec = j
	   ENDDO	       
           PRINT*,'Enter number of cross-section slices:'
	   READ*,nxsec
	   dxsec = (lim(ixsec,2)-lim(ixsec,1))/float(nxsec)
	   xsecmin = lim(ixsec,1)-dxsec
	   xsecmax = lim(ixsec,1)
	   xsecpos = xsecmin
        ENDIF      

!!--set title of plot

        IF ((.not.imulti).AND.(nacross*ndown.EQ.1)) THEN
           IF (x_sec) THEN
	      titlex = 'cross-section'
           ELSE
	      titlex = 'projection'   
           ENDIF
	   titlex = TRIM(label(ipickx))//TRIM(label(ipicky))
     &         //' '//titlex	          

	   IF (irender.GT.ndim) THEN
	   titlex = TRIM(titlex)//' - '//TRIM(label(irender))
     &     //' rendering'
           ENDIF
	   IF (ivecplot.EQ.1) titlex = ' velocity map: '//titlex
	   IF (ivecplot.EQ.2) titlex = ' magnetic field map: '//titlex
	ELSE
	   titlex = ' '
	ENDIF

!!--initialise PGPLOT     
	CALL PGBEGIN(0,'?',nacross,ndown)
!
!--set colour table
!
	IF (((irender.GT.ndim).OR.ANY(irendermulti(1:nyplots).GT.ndim))
     &	    .AND.(icolours.GT.0)) THEN
           CALL setcolours(icolours)
        ENDIF

!!------------------------------------------------------------------------      
!!--not co-ordinate plot
!!
      ELSE

        CALL PGBEGIN(0,'?',nacross,ndown)

      ENDIF
!!------------------------------------------------------------------------
! general initialisations

      IF (animate) CALL PGASK(.false.)
      IF (ndim.EQ.1) x_sec = .false.

      iploty = ipicky	
      iplotx = ipickx
!
!--save settings for rendering/vector plots
!
      iplotcont_old = iplotcont
      ncontours_old = ncontours
      ivecplot_old = ivecplot
      npix_old = npix
      npixvec_old = npixvec
      iplotpartvec_old = iplotpartvec

!
!--if plotting ghost particles, set ntotplot = ntot, else ntot=npart
!
      ntotplot(:) = npart(:)
      IF (iplotghost) ntotplot = ntot(:)
!
!--set fill style for circle plots
!      
      IF (plotcirc) CALL PGSFS(2)
!
!--increase character size depending on the number of graphs on the page
!
!      PRINT*,' Enter character height '
!      READ*,charheight
      charheight = 1.0
      IF ((ndown*nacross).GT.1) charheight = 2.0
!      charheight = 0.5*(nacross+ndown)

!!------------------------------------------------------------------------      
!! plot graphs
!!------------------------------------------------------------------------            
      over_timesteps: DO i=nstart,n_end,nfreq

	  npart1 = npart(i) + 1   	 	          
            
!-------------------------------------
! to plot multiple plots per timestep
!-------------------------------------
          over_plots: DO nyplot=1,nyplots
!--make sure character height is set correctly             
	     CALL PGSCH(charheight)	  
!--for consecutive plots	     
	     iploty = ipicky + nyplot - 1
!--where plots are specified
             IF (imulti) THEN	        
                iploty = multiploty(nyplot)
		iplotx = multiplotx(nyplot)		
	     ENDIF
!--------------------------------------------------------------
!  copy from main dat array into xplot, yplot and renderplot
!  apply transformations (log, 1/x, etc) if appropriate
!--------------------------------------------------------------
             CALL transform(dat(iplotx,:,i),xplot,itrans(iplotx),max)
	     CALL transform(dat(iploty,:,i),yplot,itrans(iploty),max)
	     labelx = transform_label(label(iplotx),itrans(iplotx))
	     labely = transform_label(label(iploty),itrans(iploty))
	     CALL transform(lim(iplotx,1),xmin,itrans(iplotx),1)
	     CALL transform(lim(iplotx,2),xmax,itrans(iplotx),1)
	     CALL transform(lim(iploty,1),ymin,itrans(iploty),1)
	     CALL transform(lim(iploty,2),ymax,itrans(iploty),1)

!--------------------------------------------------------------

!--write username, date on plot
!         IF (nacross.le.2.and.ndown.le.2) CALL PGIDEN

!--adjust plot limits if adaptive plot limits set
	 IF ((ipagechange.AND.iadapt).AND.(iplotx.LE.ndataplots)
     &      .AND.(iploty.LE.ndataplots)) THEN
	    xmin = minval(xplot(1:ntotplot(i)))
	    xmax = maxval(xplot(1:ntotplot(i)))*scalemax
	    ymin = minval(yplot(1:ntotplot(i)))
	    ymax = maxval(yplot(1:ntotplot(i)))*scalemax
	 ENDIF

!-------------------------------------------
!--plots with co-ordinates as x and y axis
!-------------------------------------------

	 IF ((iploty.le.ndim).and.(iplotx.le.ndim)) THEN
	 
	    IF (imulti) THEN
               irenderplot = irendermulti(nyplot)      
	       iplotcont = iplotcontmulti(nyplot)
	       ncontours = ncontoursmulti(nyplot)
	       ivecplot = ivecplotmulti(nyplot)
	       npix = npixmulti(nyplot)
	       npixvec = npixvecmulti(nyplot)
	       iplotpartvec = iplotpartvecmulti(nyplot)
	    ELSE
	       irenderplot = irender
	       iplotcont = iplotcont_old
	       ncontours = ncontours_old
	       ivecplot = ivecplot_old
	       npix = npix_old
	       npixvec = npixvec_old
	       iplotpartvec = iplotpartvec_old	       
	    ENDIF
!
!--rendering setup and interpolation
!
             IF (irenderplot.GT.ndim) THEN
!!--do transformations on rendering array
	        CALL transform(dat(irenderplot,:,i),renderplot(:),
     &		               itrans(irenderplot),max)
		labelrender = transform_label(label(irenderplot),
     &		                         itrans(irenderplot))

!!--if adaptive limits, find limits of rendering array		
		IF (iadapt) THEN
	           rendermin = minval(renderplot(1:ntotplot(i)))*scalemax
  	           rendermax = maxval(renderplot(1:ntotplot(i)))*scalemax		
!!--or apply transformations to fixed limits
		ELSE
		   CALL transform(lim(irenderplot,1),rendermin,
     &		               itrans(irenderplot),1)
		   CALL transform(lim(irenderplot,2),rendermax,
     &		               itrans(irenderplot),1)
                ENDIF
!
!--interpolate from particles to fixed grid using SPH summation
!		
!--do not apply transformations to the co-ordinates here		
		xmin = lim(ix(1),1)
		xmax = lim(ix(1),2)
		ymin = lim(ix(2),1)
		ymax = lim(ix(2),2)
		zmin = lim(ix(3),1)
		zmax = lim(ix(3),2)
		xminrender = MIN(xmin,ymin,zmin)
		xmaxrender = xmax
		IF (xmaxrender.GT.ymax) xmaxrender = ymax
		IF (xmaxrender.GT.zmax) xmaxrender = zmax
!!--determine number of pixels in rendered image
		pixwidth = (xmax-xmin)/REAL(npix)
		npixy = INT((ymax-ymin)/pixwidth) + 1
		npixz = INT((zmax-zmin)/pixwidth) + 1
!!--allocate memory for rendering array
		ALLOCATE ( datpix(npix,npixy) )

		SELECT CASE(ndim)
		 CASE(2)
		     CALL interpolate2D(
     &		     dat(ix(1),1:ntot(i),i),dat(ix(2),1:ntot(i),i)
     &		     dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),
     &               dat(ih,1:ntot(i),i),
     &		     dat(irenderplot,1:ntot(i),i),
     &               ntot(i),xminrender,ymin,
     &	             datpix,npix,npixy,pixwidth)
		 CASE(3)  
!!--allocate memory for 3D rendering array		 
		     ALLOCATE ( datpix3D(npix,npixy,npixz) )
		     CALL interpolate3D(
     &		     dat(ix(1),1:ntot(i),i),dat(ix(2),1:ntot(i),i)
     &		     dat(ix(3),1:ntot(i),i),dat(ipmass,1:ntot(i),i),
     &               dat(irho,1:ntot(i),i),dat(ih,1:ntot(i),i),
     &		     dat(irenderplot,1:ntot(i),i),
     &               ntot(i),xminrender,ymin,zmin,
     &	             datpix3D,npix,npixy,npixz,pixwidth)
		END SELECT
	     ENDIF
!
!--if vector plot determine whether or not to plot the particles as well
!
            iplotpart = .true.
            IF (ivecplot.GT.0) iplotpart = iplotpartvec     
	    IF (irenderplot.GT.0) iplotpart = .false
!
!--loop over cross-section slices
!
          over_cross_sections: DO k=1,nxsec
        
!!--multislice cross section (flythru)
   	     IF (x_sec.AND.flythru) THEN
	        IF (iplotpart) THEN
		   xsecmin = xsecmin+dxsec
	           xsecmax = xsecmax+dxsec
		ELSE
		   xsecpos = xsecpos + dxsec
		ENDIF
	     ELSEIF (x_sec.AND.i.EQ.1) THEN

!!--if single slice cross-section read limits for slice	    
	       DO j=1,ndim
	          IF ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j
	       ENDDO	 	       
	       PRINT*,' Cross section = ',ixsec          	       
!!--if particle cross-section, plot particles within a certain co-ordinate range	       
	       IF (iplotpart) THEN
	          PRINT 33,label(ixsec)
33            FORMAT(' Enter ',a1,' min, max for cross section slice:')
	          READ*,xsecmin,xsecmax
	          IF (xsecmax.gt.lim(ixsec,2)
     &               .or.xsecmin.lt.lim(ixsec,1)) THEN
	             PRINT*,'Warning: Cross section outside data limits' 
		  ENDIF
!!--if rendering single slice, determine position of slice
               ELSE
	          PRINT 133,label(ixsec)
133	          FORMAT(' Enter ',a1,' co-ordinate of slice ')
		  READ*,xsecpos
	          IF (xsecmin.gt.lim(ixsec,2)) xsecpos = lim(ixsec,2)
		  IF (xsecmin.lt.lim(ixsec,1)) xsecpos = lim(ixsec,1)		  	          
	       ENDIF	 
	     ENDIF   
!
!--take cross section
!
            IF (ndim.EQ.3) THEN
	       SELECT CASE(ixsec)	! which direction to slice in
                CASE(1)
                 ipixxsec = INT((xsecpos - xmin)/npixx)
                 datpix = datpix3D(ipixxsec,:,:)
                CASE(2)
                 ipixxsec = INT((xsecpos - ymin)/npixy)
                 datpix = datpix3D(:,ipixxsec,:)
                CASE(3)
                 ipixxsec = INT((xsecpos - zmin)/npixz)
                 PRINT*,' cross section = ',ixsec,ipixxsec
                 datpix = datpix3D(:,:,ipixxsec)
               END SELECT	       		  	       
	    ENDIF   
	    	    
!----------------------------------------------------------------------------
!--set up PGPLOT page
!
	    IF ((ipagechange).OR.((.not.ipagechange).AND.(i.EQ.nstart))) THEN
!	       CALL PGENV(lim(iplotx,1),lim(iplotx,2),
!     &                 lim(iploty,1),lim(iploty,2),1,1)	! 0 for no axes
	       CALL PGPAGE
	       IF (nacross*ndown.GT.1) THEN
!	          IF (imulti) THEN
		     CALL PGSVP(0.2,0.8,0.2,0.98)
!		  ELSE
!		     CALL PGSVP(0.0,1.0,0.0,1.0)
!		  ENDIF   
	       ELSE
	          CALL PGSVP(0.1,0.9,0.1,0.9)
	       ENDIF	  
	       CALL PGWNAD(xmin,xmax,ymin,ymax)	!  pgwnad does equal aspect ratios
	       CALL PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)	       
            ELSEIF (nyplot.EQ.1) THEN
	       CALL PGPANL(1,1)
	    ELSE
	       CALL PGPAGE
	    ENDIF
!----------------------------------------------------------------------------	    
!--print plot limits to screen
	    PRINT 34, time(i)
	    PRINT*,TRIM(labely),'min,max = ',ymin,ymax
	    PRINT*,TRIM(labelx),'min,max = ',xmin,xmax
34          FORMAT (5('-'),' t = ',f8.4,1x,15('-'))
    
!--set plot limits	    
	    CALL PGWNAD(xmin,xmax,ymin,ymax)	! pgwnad does equal aspect ratios
!--label plot
	    IF (((nyplots-nyplot).LT.nacross).OR.(.not.isamexaxis)) THEN
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
	      CALL PGLABEL(labelx,' ',titlex)	    
	    ELSE
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
!	      CALL PGLABEL(' ',labely,titlex)
	    ENDIF

	    IF (x_sec.AND.iplotpart) PRINT 35,label(ixsec),xsecmin,
     &                          label(ixsec),xsecmax
35          FORMAT('Cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

            xminrender = MIN(xmin,ymin)
	    xmaxrender = xmax
	    IF (xmaxrender.LT.ymax) xmaxrender = ymax

cc--density/scalar field rendering	    
	    IF (irenderplot.GT.ndim) THEN

	       PRINT*,TRIM(labelrender),' min, max = ',rendermin,rendermax	       
	       SELECT CASE(ndim)
	        CASE(2)
	           CALL smooth_render(
     &	               xplot(1:ntot(i)),yplot(1:ntot(i)),
     &	               xminrender,xmaxrender,ymin,ymax,
     &                 renderplot(1:ntot(i)),
     &                 rendermin,rendermax,labelrender,
     &                 npart(i),ntot(i),npix,
     &                 icolours,iplotcont,ncontours,
     &                 dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),
     &                 dat(ih,1:ntot(i),i))
		CASE(3)
!--do not apply transformations to the co-ordinates here		
		   xmin = lim(ix(1),1)
		   xmax = lim(ix(1),2)
		   ymin = lim(ix(2),1)
		   ymax = lim(ix(2),2)
		   zmin = lim(ix(3),1)
		   zmax = lim(ix(3),2)
		   xminrender = MIN(xmin,ymin,zmin)
		   xmaxrender = xmax
		   IF (xmaxrender.GT.ymax) xmaxrender = ymax
		   IF (xmaxrender.GT.zmax) xmaxrender = zmax
		   CALL render3D(
     &		       dat(ix(1),1:ntot(i),i),dat(ix(2),1:ntot(i),i),
     &                 dat(ix(3),1:ntot(i),i),
     &	               xminrender,xmaxrender,ymin,ymax,zmin,zmax,
     &                 renderplot(1:ntot(i)),
     &                 rendermin,rendermax,labelrender,
     &                 npart(i),ntot(i),npix,
     &                 icolours,iplotcont,ncontours,
     &                 dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),
     &                 dat(ih,1:ntot(i),i),ixsec,xsecmin)
	       END SELECT
	       
	    ELSE
cc--particle plots

!
!--if particle cross section, plot particles only in a defined coordinate range
!
            IF (x_sec.AND.iplotpart) THEN
	       DO j=1,npart(i)
	          IF ((dat(ixsec,j,i).lt.xsecmax)
     &		  .and.(dat(ixsec,j,i).gt.xsecmin)) THEN
     		      CALL PGPT(1,xplot(j),yplot(j),imark)
		   IF (plotcirc) THEN
		      CALL PGCIRC(xplot(j),yplot(j),2.*dat(ih,j,i))
		   ENDIF
		  ENDIF
	       ENDDO
	       DO j=npart1,ntotplot(i)
	          IF ((dat(ixsec,j,i).lt.xsecmax)
     &		  .and.(dat(ixsec,j,i).gt.xsecmin)) THEN
     		      CALL PGPT(1,xplot(j),yplot(j),imarkg)
		  ENDIF	       
	       ENDDO
!
!--or simply plot all particles
!
	    ELSE	     	     
	     
cc--plot particle positions
	       IF (iplotpart) CALL PGPT(npart(i),xplot(1:npart(i)),
     &                                  yplot(1:npart(i)),imark)
cc--plot ghost particles with different marker
               IF ((ntotplot(i).GT.npart(i)).AND.(iplotpart)) THEN
	          CALL PGPT(ntot(i)-npart(i),xplot(npart1:ntot(i)),
     &                          yplot(npart1:ntot(i)),imarkg)
               ENDIF
cc--plot circles of interaction
	       IF (plotcirc) THEN		  
	         IF (plotcircall) THEN
		  PRINT*,'plotting circles of interaction',npart(i) 
		  DO j=1,npart(i)
		     CALL PGCIRC(xplot(j),yplot(j),2.*dat(ih,j,i))
		  ENDDO
		 ELSE 
		  PRINT*,'plotting circles of interaction',icircpart
		  CALL PGCIRC(xplot(icircpart),yplot(icircpart),2*dat(ih,icircpart,i))
                 ENDIF
	       ENDIF
	       IF (ilabelpart) THEN
cc--plot particle labels
                  PRINT*,'plotting particle labels ',ntotplot(i)
		  DO j=1,ntotplot(i)
		     CALL PGNUMB(j,0,1,string,nc)
		     CALL PGSCH(0.5*charheight)
		     CALL PGTEXT(xplot(j),yplot(j),string(1:nc))
		     CALL PGSCH(charheight)
		  ENDDO
	       ENDIF	! ilabelpart
             
	     ENDIF	! if x_sec else    
	    ENDIF	! if irender
	    
	    	    
cc--velocity vector map
	    IF (ivecplot.EQ.1 .AND. ivx.NE.0) THEN
	       DO j=1,ntotplot(i)
	          vecplot(1,j) = dat(iplotx+ivx-1,j,i)
		  vecplot(2,j) = dat(iploty+ivx-1,j,i)
	       ENDDO
	       vmax = lim(iplotx+ivx-1,2)
	       IF (lim(iploty+ivx-1,2).gt.vmax) vmax = lim(iploty+ivx-1,2)
	       vmin = min(lim(iplotx+ivx-1,1),lim(iploty+ivx-1,1))
	       CALL coarse_render(
     &              xplot(1:ntotplot(i)),
     &              yplot(1:ntotplot(i)),
     &              xminrender,xmaxrender,vecplot(:,1:ntotplot(i)),
     &              vmin,vmax,ntotplot(i),npixvec,2,icolours,iplotcont)

!	       scale = 0.08*(limx(iploty,2)-limx(iploty,1))
!	       CALL PGSCH(0.35)
!	       DO j=1,ntotplot(i)
!	       CALL PGARRO(x(iploty,j,i),x(iplotx,j,i),
!     &	                   x(iploty,j,i)+v(iploty,j,i)*scale,
!     &			   x(iplotx,j,i)+v(iplotx,j,i)*scale)
!	       ENDDO
!	       CALL PGSCH(1.0)    
	    ELSEIF ((ivecplot.EQ.2).AND.(iBfirst.NE.0)) THEN
cc--plot vector map of magnetic field
	       PRINT*,'plotting magnetic field: ',
     &          label(iBfirst+iplotx-1),label(iBfirst+iploty-1)
	       DO j=1,ntotplot(i)
	          vecplot(1,j) = dat(iBfirst+iplotx-1,j,i)
		  vecplot(2,j) = dat(iBfirst+iploty-1,j,i)
	       ENDDO

	       CALL coarse_render(
     &              xplot(1:ntotplot(i)),
     &              yplot(1:ntotplot(i)),
     &              xminrender,xmaxrender,vecplot(1:2,1:ntotplot(i)),
     &              Bmin,Bmax,ntotplot(i),npixvec,2,icolours,iplotcont)	    

	    ENDIF
	    
	    IF (nyplot.EQ.1) CALL legend(time(i))	    
	    
	    ENDDO over_cross_sections
	    
	    IF (ALLOCATED(datpix)) DEALLOCATE(datpix)
	    IF (ALLOCATED(datpix3D)) DEALLOCATE(datpix3D) 
	    
c--------------------------------
cc--not both coordinates
c--------------------------------
	 ELSEIF ((iploty.gt.ndim .or. iplotx.gt.ndim)
     &     .AND.(iploty.LE.ndataplots .AND. iplotx.LE.ndataplots)) THEN
	    
	    IF ((ipagechange).OR.
     &	       ((.not.ipagechange).AND.(i.EQ.nstart))) THEN
!	       CALL PGENV(limx(iplotx,1),limx(iplotx,2),
!     &                 lim(iploty,1),lim(iploty,2),0,0)	! 0 for no axes
	       CALL PGPAGE
	       IF ((nacross*ndown).GT.1) THEN
	          IF (imulti) THEN
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ELSE
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ENDIF
	       ELSE
	          CALL PGSVP(0.1,0.9,0.1,0.9)	       
	       ENDIF
	       CALL PGSWIN(xmin,xmax,ymin,ymax)
	       CALL PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)	       	       
            ELSEIF (nyplot.EQ.1) THEN
	       CALL PGPANL(1,1)
	    ELSE
	       CALL PGPAGE
	    ENDIF

!--print plot limits to screen
	    PRINT 34, time(i)
	    PRINT*,TRIM(labely),'min,max = ',ymin,ymax
	    PRINT*,TRIM(labelx),'min,max = ',xmin,xmax
!
!--set plot limits
!	    
	    CALL PGSWIN(xmin,xmax,ymin,ymax)

	    IF (((nyplots-nyplot).LT.nacross).OR.(.not.isamexaxis)) THEN
!--print x and y labels
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
	      CALL PGLAB(labelx,' ',title)	    
	    ELSE
!--print y labels only	    
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
!	      CALL PGLAB(' ',labely,title)
	    ENDIF
	    
	    IF ((i.EQ.nstart).AND.iplotlinein) THEN	! plot initial conditions as dotted line
	       CALL PGSLS(linestylein)
	    ELSE
!--plot time on plot
	       IF (nyplot.EQ.1) CALL legend(time(i))
!--plot particles
	       CALL PGSLS(1)
	       CALL PGSCH(1.0)	! reset character height before plotting particles
	       CALL PGPT(npart(i),xplot(1:npart(i)),yplot(1:npart(i)),imark)
	    ENDIF      
!--plot line joining the particles
	    IF (iplotline.OR.(iplotlinein.AND.(i.EQ.nstart))) THEN
	       CALL PGLINE(npart(i),xplot(1:npart(i)),yplot(1:npart(i)))     
            ENDIF
!--plot ghost particles with different marker
	    IF (ntotplot(i).GT.npart(i))
     &	       CALL PGPT(ntot(i)-npart(i),xplot(npart1:ntot(i)),
     &                                    yplot(npart1:ntot(i)),imarkg)
            
	    IF (plotcirc) THEN
	       PRINT*,'plotting circles of interaction',npart(i)
	       DO j=1,npart(i)
		  CALL PGCIRC(xplot(j),yplot(j),2.*dat(ih,j,i))
	       ENDDO
	    ENDIF
            CALL PGSLS(1)	! reset 
            CALL PGSCH(charheight)
!
!--plot average line
!
	    IF (iplotav) CALL plot_average(xplot(1:npart(i)),
     &                        yplot(1:npart(i)),npart(i),nbins)
!
!--plot particle labels
!
	    IF (ilabelpart) THEN
               PRINT*,'plotting particle labels ',ntotplot(i)
	       DO j=1,ntotplot(i)
		  CALL PGNUMB(j,0,1,string,nc)
		  CALL PGSCH(0.5*charheight)
		  CALL PGTEXT(xplot(j),yplot(j),string(1:nc))
		  CALL PGSCH(charheight)
	       ENDDO
	    ENDIF	! ilabelpart

         ELSEIF (iploty.LE.numplot) THEN
!---------------------------------------------------
! additional plots (not plots of particle data)
!---------------------------------------------------
            ! e.g. call routine to do convergence plot here
	    IF ((ipagechange).OR.((.not.ipagechange).AND.(i.EQ.nstart))) THEN
	       CALL PGPAGE
	       IF ((nacross*ndown).GT.1) THEN
	          IF (imulti) THEN
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ELSE
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ENDIF
	       ELSE
	          CALL PGSVP(0.1,0.9,0.1,0.9)	       
	       ENDIF
!	       CALL PGSWIN(lim(iplotx,1),lim(iplotx,2),
!     &	                   lim(iploty,1),lim(iploty,2))
!	       CALL PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)	       	       
            ELSEIF (nyplot.EQ.1) THEN
	       CALL PGPANL(1,1)
	    ELSE
	       CALL PGPAGE
	    ENDIF

	    IF (iexact.EQ.4) THEN
	       CALL exact_toystar_ACplane(Atstar,Ctstar,sigma,gamma)
	    ENDIF 

	    IF (nyplot.EQ.1) CALL legend(time(i))
	      
	
	 ELSE	 ! if plot not in correct range
	    CALL PGPAGE
	
	 ENDIF      ! ploty = whatever

!-----------------------------
!!--plot exact solution
!-----------------------------

         SELECT CASE(iexact)
         CASE(1)		! polytrope
	    IF ((iploty.eq.irho).and.(iplotx.eq.irad))
     &	    CALL PGLINE(ipolyc,rad(1:ipolyc),den(1:ipolyc))
	 
	 CASE(2)	! soundwave
	    IF ((iploty.eq.irho).and.(iplotx.eq.1).and.(i.ne.1)) THEN
	       CALL exact_swave(time(i),delta,lambda,gamma,
     &              xplot(1:npart(i)),yplot(1:npart(i)),
     &              dat(iutherm,1:npart(i),i),npart(i))
            ENDIF
	    
	 CASE(3)	! Sedov blast wave
	    IF ((iploty.eq.irho).and.(iplotx.eq.irad))
     &	    CALL exact_sedov(time(i),gamma,xplot(1:npart(i)),
     &                      yplot(1:npart(i)),npart(i))
	 
	 CASE(4)	! Toy star
!	    totmass = SUM(dat(ipmass,1:npart(i),i))
!	    Htstar = (0.75*totmass)**(2./3.)*Ctstar**(1./3.)
!	    Htstar = 1.0    
!	    PRINT*,' totmass,H,A,C in = ',totmass,Htstar,Atstar,Ctstar
	    IF (iplotx.EQ.1) THEN	! if x axis is x coordinate
	       IF (iploty.EQ.irho) THEN
	          CALL exact_toystar(time(i),gamma,
     &		        Htstar,Atstar,Ctstar,sigma,norder,1)
	       ELSEIF (iploty.EQ.ipr) THEN
	          CALL exact_toystar(time(i),gamma,
     &			Htstar,Atstar,Ctstar,sigma,norder,2)	       
	       ELSEIF (iploty.EQ.iutherm) THEN
	          CALL exact_toystar(time(i),gamma,
     &			Htstar,Atstar,Ctstar,sigma,norder,3)	       
	       ELSEIF (iploty.EQ.ivx) THEN
	          CALL exact_toystar(time(i),gamma,
     &			Htstar,Atstar,Ctstar,sigma,norder,4)	       
	       ELSEIF (iploty.EQ.iBfirst+1) THEN
	          CALL exact_toystar(time(i),gamma,
     &			Htstar,Atstar,Ctstar,sigma,norder,5)
	       ENDIF
	    ELSEIF (iplotx.EQ.irho) THEN
	       IF (iploty.EQ.iBfirst+1) THEN
	          CALL exact_toystar(time(i),gamma,
     &			Htstar,Atstar,Ctstar,sigma,norder,6)	       
	       ENDIF   
	    ENDIF
	    
	    IF (iploty.EQ.iextra) THEN	! plot point on A-C plane
	       CALL exact_toystar(time(i),gamma,
     &	                Htstar,Atstar,Ctstar,sigma,norder,7)
	    ENDIF
	 END SELECT
!
!--plot h = 1/rho
!
	 IF ((iploty.EQ.ih).AND.(iplotx.EQ.irho)) THEN
	    CALL plot_rhoh(hfact,ndim)
	 ENDIF
	
         ENDDO over_plots	! over plots per timestep (nyplot)
      
      ENDDO over_timesteps
      
      IF (animate) THEN
         PRINT*,'press RETURN to finish'
	 READ*
      ENDIF
c ------------------------------------------------------------------------      

300   CONTINUE
      CALL PGEND
      GOTO 200

999   CONTINUE                 
      END

!-----------------------------------------------------------------------------
! options for rendering / vector plots
!-----------------------------------------------------------------------------

      SUBROUTINE get_render_options(irender,npix,icolours,iplotcont,
     &	        ncontours,ivecplot,npixvec,iplotpartvec,ndim,numplot)     
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim, numplot
      INTEGER, INTENT(OUT) :: irender,npix,icolours,ivecplot,npixvec
      INTEGER, INTENT(OUT) :: ncontours
      LOGICAL, INTENT(OUT) :: iplotcont,iplotpartvec
      CHARACTER(LEN=1) :: ans      
!
!--rendering options
!
      PRINT*,' Enter field (from menu) for rendering (0=none): '
      READ*,irender
      IF ((irender.GT.ndim)
     &	 .AND.(irender.LE.numplot)) THEN 
      
         PRINT*,' Enter number of pixels'
         READ*,npix
      
         icolours = -1
         DO WHILE (icolours.LT.0)
	    PRINT*,' Enter colour scheme for rendering '
            PRINT*,' (-ve for demo, 0 for contours only)'
	    READ*,icolours
	    IF (icolours.LT.0) CALL colour_demo
         ENDDO
      
         IF (icolours.NE.0) THEN
	    PRINT*,' Plot contours (y/n)?'
	    READ*,ans
	    iplotcont = .false.
	    IF (ans.EQ.'y'.OR.ans.EQ.'Y') iplotcont = .true.
         ELSE
	    iplotcont = .true.
         ENDIF
 
         IF (iplotcont) THEN
	    PRINT*,' Enter number of contours between min,max'
   	    READ*,ncontours	     
         ENDIF
      ELSE
	 irender = 0   
      ENDIF
!
!--vector plot options
!	  
111   PRINT*,'vector plot: '//
     &       'velocity (1), magnetic field (2) or none(0)?'
      READ*,ivecplot
      IF ((ivecplot.GT.2).OR.(ivecplot.LT.0)) THEN
	 GOTO 111
      ELSEIF (ivecplot.GT.0) THEN	
	 PRINT*,' Enter number of pixels'
	 READ*,npixvec
         ans = ' '
         IF (irender.EQ.0) THEN
	    PRINT*,' Plot particles (y/n) '
	    READ*,ans
	 ENDIF
	 iplotpartvec = .false.
	 IF (ans.EQ.'y' .OR. ans.EQ.'Y') iplotpartvec = .true.	
      ENDIF
    
      RETURN
      END SUBROUTINE get_render_options
