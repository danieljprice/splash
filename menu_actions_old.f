!!
!!    implements menu options
!!
      SUBROUTINE menu_actions(ipicky)      
      USE exact_params
      USE filenames
      USE labels
      USE settings
      USE multiplot
      USE particle_data
      USE prompting
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ipicky
      INTEGER :: i,j,iaction,ipick
      REAL, PARAMETER :: pi = 3.1415926536
      REAL :: diff, mid, temp, Bmag
      REAL :: angledeg,anglexy,runit(ndimmax)	! to plot r at some angle
      CHARACTER(LEN=30) :: filename
      CHARACTER(LEN=25) :: transform_label
      CHARACTER(LEN=1) :: ans,dummy
      
      iaction = ipicky - numplot      
      
      SELECT CASE(iaction)
            
!------------------------------------------------------------------------
       CASE(2)
          IF (.not.ihavereadfilename) THEN
	     CALL prompt('Enter filename to read',rootname)
!	     PRINT*,'Enter filename to read:' 	!(5 characters):'
!	     READ(*,30) rootname
!30           FORMAT(a)	   	  
	  ENDIF
	  ihavereadfilename = .false.
	  nfilesteps = maxstep
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
	     RETURN
8801         CONTINUE
             PRINT*,'Cannot open/ error reading ',filename	     
	  ENDIF
	  RETURN
!------------------------------------------------------------------------	  
       CASE(3)
          PRINT*,' Start, end at timestep numbers:'
	  READ*,nstart,n_end
	  IF (n_end.gt.nfilesteps) THEN
	     PRINT*,'n_end greater than nfilesteps, reset to ',nfilesteps
	     n_end = nfilesteps
	  ENDIF
	  PRINT*,' Frequency of steps to READ:'
	  READ*,nfreq
          PRINT *,' Steps = ',(n_end-nstart+1)/nfreq
	  RETURN
!------------------------------------------------------------------------
       CASE(4)
          animate=.not.animate
          PRINT *,' Animation = ',animate
	  RETURN
!------------------------------------------------------------------------
       CASE(5)
          iadapt=.not.iadapt
          PRINT *,' Adaptive plot limits = ',iadapt
	  RETURN
!------------------------------------------------------------------------
       CASE(6)
          magfield=.not.magfield
          PRINT *,' Mag field = ',magfield
	  RETURN
!------------------------------------------------------------------------
       CASE(7)
          x_sec =.not.x_sec
	  PRINT *,' Cross section = ',x_sec
	  flythru = .false.
	  IF (x_sec) THEN
	     PRINT*,'Do you want a fly-through (y/n)?'
	     READ*,ans
	     IF (ans.eq.'y'.or.ans.eq.'Y') flythru=.true.
	  ENDIF
	  RETURN
!------------------------------------------------------------------------
       CASE(8)
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
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(9)
          PRINT*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
	  PRINT*,' Enter PGPLOT marker # (particles):'
	  READ*,imark
	  PRINT*,' Enter PGPLOT marker # (ghosts):'
	  READ*,imarkg
	  RETURN
!------------------------------------------------------------------------
       CASE(10)
	  PRINT*,'Enter number of plots across:'
	  READ*,nacross
	  PRINT*,'Enter number of plots down'
	  READ*,ndown
	  RETURN	    	  
!------------------------------------------------------------------------
       CASE(11)
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
     &               iplotpartvecmulti(i),x_secmulti(i),
     &               xsecposmulti(i),ndim,numplot)
             ENDIF	     
	  ENDDO
	  RETURN	    	  
!------------------------------------------------------------------------
       CASE(12)
	  ipagechange=.not.ipagechange
          PRINT*,' Page changing = ',ipagechange
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(13)
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
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(14)
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
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(15)
	  iplotav=.not.iplotav
          IF (iplotav) THEN
51	     PRINT*,' Enter number of bins for averaging '
	     READ*,nbins
	     IF ((nbins.LE.0).OR.(nbins.GT.1000)) GOTO 51
          ENDIF
	  PRINT*,' Plot average, nbins = ',iplotav,nbins
	  RETURN 		  	    	  	  
!------------------------------------------------------------------------
       CASE(16)
c	  label particles with particle numbers
	  ilabelpart=.not.ilabelpart
          PRINT*,' label particles = ',ilabelpart
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(17)
c	  plot ghost particles?
	  iplotghost=.not.iplotghost
          PRINT*,' plot ghost particles = ',iplotghost
	  IF (iplotghost) ntotplot(:) = ntot(:)
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(18)
          CALL get_render_options(irender,npix,icolours,iplotcont,
     &	        ncontours,ivecplot,npixvec,iplotpartvec,
     &          x_sec,xsecpos,ndim,numplot)
          RETURN
!------------------------------------------------------------------------
       CASE(19)
299       CONTINUE
	  PRINT*,'Enter plot number to apply transformation '
          PRINT*,'(0 to finish, -ve to set all)'
	  READ*,ipick
	  IF (ipick.LE.numplot .AND. ipick.NE.0) THEN
	     WRITE(*,*) (i,') ',TRIM(transform_label('x',i)), i=1,5)
             PRINT*,'Enter transformation to apply:'
	     IF (ipick.LT.0) THEN
	        READ*,ipick
		itrans(:) = ipick
		RETURN
	     ELSE
	        READ*,itrans(ipick)
	        GOTO 299
	     ENDIF	
	  ELSE
	     RETURN
	  ENDIF
!------------------------------------------------------------------------
       CASE(20)
	  IF (.not.iadapt) THEN
	     PRINT*,'Do you want to manually enter limits?'
	     READ*,ans
	     IF (ans.EQ.'y' .OR. ans.EQ.'Y') THEN
301	        PRINT*,'Enter plot number to set limits (0 to finish)'
		READ*,ipick
		IF ((ipick.LE.numplot).AND.(ipick.GT.0)) THEN
		   PRINT*,' Current limits   (min,max):',
     &		          lim(ipick,1),lim(ipick,2)
		   PRINT*,' Enter new limits (min,max):'
		   READ*,lim(ipick,1),lim(ipick,2)
		   PRINT*,' >> limits set    (min,max):',
     &		         lim(ipick,1),lim(ipick,2)
		ELSE
		   PRINT*,'Finished, returning'
		   RETURN
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
	  RETURN
!------------------------------------------------------------------------
       CASE(21)
!	  show/hide plot options
	  ishowopts = .not.ishowopts
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(22)
          CALL write_defaults
	  RETURN
       CASE DEFAULT
          STOP ' Error: menu action not defined'  
      
      END SELECT
      
      RETURN      
      END SUBROUTINE menu_actions
