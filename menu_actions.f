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
      REAL :: papersizey
      CHARACTER(LEN=30) :: filename
      CHARACTER(LEN=25) :: transform_label,crap
      CHARACTER(LEN=1) :: ans,dummy
      LOGICAL :: ians
      
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
	  CALL read_data(rootname,ifile,gamma(:),time(:),
     &         dat(:,:,:),iam(:),
     &         ncolumns,ncalc,nfilesteps,ntot(:),npart(:),
     &         nghost(:),ndim,ndimV,hfact,ivegotdata)
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
     &		   dat(ipr,j,i) = dat(irho,j,i)*dat(iutherm,j,i)*(gamma(i)-1.)
!!--entropy
                IF (ientrop.NE.0) THEN
  	           IF (dat(irho,j,i).GT.1.e-10) THEN
                      dat(ientrop,j,i) = dat(ipr,j,i)
     &		                        /dat(irho,j,i)**gamma(i)
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
		IF ((idivBerr.NE.0).AND.
     &		(idivB.NE.0).AND.(ih.NE.0).AND.(iBfirst.NE.0)) THEN
		   Bmag = SQRT(DOT_PRODUCT(dat(iBfirst:iBlast,j,i),
     &		                           dat(iBfirst:iBlast,j,i)))
                   IF (Bmag.GT.0.) THEN
		     dat(idivBerr,j,i) = 
     &		         abs(dat(idivB,j,i))*dat(ih,j,i)/Bmag
		   ELSE
!		     PRINT*,'Bmag = 0 !!'
		     dat(idivBerr,j,i) = 0.
		   ENDIF
		ENDIF
!!--h*SQRT(rho)/abs(B)	(timestep)
		IF ((itimestep.NE.0).AND.
     &		(ih.NE.0).AND.(irho.NE.0).AND.(iBfirst.NE.0)) THEN
		   Bmag = SQRT(DOT_PRODUCT(dat(iBfirst:iBlast,j,i),
     &		                           dat(iBfirst:iBlast,j,i)))
                   IF (Bmag.NE.0.) THEN
		     dat(itimestep,j,i) = 
     &		         dat(ih,j,i)*SQRT(dat(irho,j,i))/Bmag
		   ELSE
		     PRINT*,'Bmag = 0 !!'
		     dat(itimestep,j,i) = 0.
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
          CALL prompt('Start at timestep ',nstart,1,nfilesteps)
	  CALL prompt('End at timestep   ',n_end,nstart)
	  IF (n_end.gt.nfilesteps) THEN
	     PRINT*,'n_end greater than nfilesteps, reset to ',nfilesteps
	     n_end = nfilesteps
	  ENDIF
	  CALL prompt(' Frequency of steps to read',nfreq,1,nfilesteps)
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
          xsec_nomulti =.not.xsec_nomulti
	  PRINT *,' Cross section = ',xsec_nomulti
	  flythru = .false.
	  IF (xsec_nomulti) THEN
	     CALL prompt('Do you want a fly-through',flythru)
!	     READ*,ans
!	     IF (ans.eq.'y'.or.ans.eq.'Y') flythru=.true.
	  ENDIF
	  RETURN
!------------------------------------------------------------------------
       CASE(8)
	  plotcirc=.not.plotcirc
          PRINT*,' Plot circles of interaction = ',plotcirc
	  IF (plotcirc) THEN	     
	     CALL prompt('Plot all circles?',plotcircall)
	     IF (.NOT.plotcircall) THEN
	        CALL prompt('Enter particle number to plot circle around',
     &		            icircpart,1,MAXVAL(ntot))
	     ENDIF
	  ENDIF
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(9)
          PRINT*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
	  CALL prompt(' Enter PGPLOT marker # (particles):',imark)
	  CALL prompt(' Enter PGPLOT marker # (ghosts)   :',imarkg)
	  CALL prompt(' Enter PGPLOT marker # (sinks)    :',imarksink)
	  RETURN
!------------------------------------------------------------------------
       CASE(10)
          CALL prompt('Enter number of plots across:',nacross,1,numplot)
	  CALL prompt('Enter number of plots down  :',ndown,1,numplot)
	  RETURN	    	  
!------------------------------------------------------------------------
       CASE(11)
	  CALL prompt('Enter number of plots per timestep:',
     & 	              nyplotmulti,1,numplot)
!	  READ*,nyplotmulti
!          IF (nyplotmulti.GT.numplot) THEN 
!	     nyplotmulti = numplot
!	     PRINT*,'nyplots > numplot, reset to ',numplot
!	  ENDIF
	  nacross = nyplotmulti/2
	  IF (nacross.EQ.0) nacross = 1
	  ndown = nyplotmulti/nacross
	  PRINT*,'setting nacross,ndown = ',nacross,ndown  
          CALL prompt('Same x axis for all?',ians)
	  IF (ians) THEN
 	     CALL prompt('Enter x axis for all plots',multiplotx(1),
     &                   1,numplot)	     
	     multiplotx(2:nyplotmulti) = multiplotx(1)	     
	  ENDIF
	  DO i=1,nyplotmulti
	     PRINT*,'Plot number ',i,':'
	     CALL prompt(' y axis ',multiploty(i),1,numplot)
	     IF (.NOT.ians.AND.multiploty(i).le.ndataplots) THEN
	        CALL prompt(' x axis ',multiplotx(i),1,numplot)
	     ENDIF
	     IF ((multiplotx(i).LE.ndim)
     &	         .AND.(multiploty(i).LE.ndim)) THEN
                CALL get_render_options(irendermulti(i),
     &		     npixmulti(i),icolours,iplotcontmulti(i),
     &	             ncontoursmulti(i),ivecplotmulti(i),npixvecmulti(i),
     &               iplotpartvecmulti(i),x_secmulti(i),
     &               xsecposmulti(i),backgnd_vec_multi(i),ndim,numplot)
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
	  PRINT*,' 0) PGPLOT default'
	  PRINT*,' 1) small square movie '
	  PRINT*,' 2) large/multiple movie'
	  PRINT*,' 3) single small graph'
	  PRINT*,' 4) duo small graph '
	  PRINT*,' 5) Custom size ' 
          CALL prompt(' Enter option for paper size ',ipapersize,0,5)
	  SELECT CASE(ipapersize)
	     CASE(1) 
	        papersizex = 0.25*11.7
		aspectratio = 1.0
	     CASE(2)
	        papersizex = 0.5*11.7
		aspectratio = 1.0
	     CASE(3) 
	        papersizex = 0.5*11.7 
		aspectratio = 1./SQRT(2.)
	     CASE(4)
	        papersizex = 11.7
		aspectratio = 0.5/SQRT(2.)	
	     CASE(5)
	        CALL prompt(' x size (inches) ',papersizex,0.0,12.0)
		CALL prompt(' y size (inches) or aspect ratio (-ve)',
     &		            papersizey,-12.0,12.0)
                IF (papersizey.LT.0.0) THEN
		   aspectratio = abs(papersizey)
		ELSE
                   aspectratio = papersizey/papersizex
                ENDIF
	     CASE DEFAULT
	        papersizex = 0.0	! no call to PGPAP if they are zero
		aspectratio = 0.0	
	  END SELECT
	  RETURN 	  
	  
!------------------------------------------------------------------------
       CASE(14)
	  PRINT*,' Plot initial only(i), all(a), both(b) or not (n)?'
	  READ*,ans
	  iplotline = .false.
	  iplotlinein = .false.
	  IF (ans.EQ.'i'.OR.ans.EQ.'b') iplotlinein = .true.
	  IF (ans.EQ.'a'.OR.ans.EQ.'b') iplotline = .true.
	  IF (iplotlinein) THEN
	     CALL prompt('Enter PGPLOT line style',linestylein,0,5)
	  ENDIF
          PRINT*,' Plot line = ',iplotline,iplotlinein
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(15)
	  IF (iexact.NE.0) THEN
	     iexact = 0
	  ELSE   
40           WRITE (*,41)
41           FORMAT(' 0) none ',/,
     &              ' 1) polytrope ',/,
     &              ' 2) soundwave ',/,     
     &              ' 3) sedov blast wave',/,     
     &              ' 4) toy star ',/,
     &              ' 5) MHD shock tubes ')
             CALL prompt('Enter exact solution to plot',iexact,0,5)
	     PRINT*,' plotting exact solution number ',iexact
!
!--enter parameters for soundwave solution
!
	     IF (iexact.EQ.4) THEN
	        PRINT*,' Toy star: '
		CALL prompt('Enter parameter A (v = Ax) ',Atstar)
                CALL prompt('Enter parameter H (rho = H - Cx^2)',Htstar)
		CALL prompt('Enter parameter C (rho = H - Cx^2)',
     &                       Ctstar,0.0)		
		sigma = 0.
		IF (magfield) THEN
		   CALL prompt('Enter parameter sigma (By = sigma rho)',
     &                         sigma0)
		   sigma = sigma0
		ENDIF
		ians = .false.
		CALL prompt('Do you want oscillations?',ians)
		norder = -1
		IF (ians) THEN
		   CALL prompt('Enter order',norder,0)
		ENDIF
	     ELSEIF (iexact.EQ.2) THEN
	        CALL prompt('Enter wavelength of sound wave lambda ',
     &                      lambda,0.0)		
                CALL prompt('Enter amplitude ',delta,0.0)
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
       CASE(16)
	  iplotav=.not.iplotav
          IF (iplotav) THEN
	     CALL prompt('Enter number of bins for averaging ',nbins,1,1000)
          ENDIF
	  PRINT*,' Plot average, nbins = ',iplotav,nbins
	  RETURN 		  	    	  	  
!------------------------------------------------------------------------
       CASE(17)
c	  label particles with particle numbers
	  ilabelpart=.not.ilabelpart
          PRINT*,' label particles = ',ilabelpart
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(18)
c	  plot ghost particles?
	  CALL prompt('Plot ghost particles? ',iplotghost)
	  CALL prompt('Plot sink particles? ',iplotsink)
          PRINT*,' plot ghost particles = ',iplotghost
	  PRINT*,' plot sink particles = ',iplotsink
	  IF (iplotghost) ntotplot(:) = npart(:) + nghost(:)
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(19)
          CALL get_render_options(irender,npix_nomulti,icolours,
     &	       iplotcont_nomulti,ncontours_nomulti,
     &	       ivecplot_nomulti,npixvec_nomulti,iplotpartvec_nomulti,
     &         xsec_nomulti,xsecpos_nomulti,backgnd_vec_nomulti,
     &         ndim,numplot)
          RETURN
!------------------------------------------------------------------------
       CASE(20)
299       CONTINUE
          ipick = 0
	  PRINT*,'Enter plot number to apply transformation '
!          PRINT*,'(0 = finish, -ve to set all)'
	  CALL prompt('(0 = finish, -1 = set all) ',ipick)
	  IF (ipick.LE.numplot .AND. ipick.NE.0) THEN
	     WRITE(*,300) (i,TRIM(transform_label('x',i)), i=1,5)
300          FORMAT(1x,i1,') ',a)		
             PRINT*,'Enter transformation to apply'//
     &	      ' (or a combination e.g. 321)'
	     IF (ipick.LT.0) THEN
	        ipick = 0
	        CALL prompt(' ',ipick,0)
		itrans(:) = ipick
		RETURN
	     ELSE	        
	        CALL prompt(' ',itrans(ipick),0)
	        GOTO 299
	     ENDIF	
	  ELSE
	     RETURN
	  ENDIF
!------------------------------------------------------------------------
       CASE(21)
	  IF (.not.iadapt) THEN
	     ians = .false.
	     CALL prompt('Do you want to manually enter limits?',ians)
	     IF (ians) THEN
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
	        CALL prompt('Enter zoom factor for fixed limits',zoom,0.0)
	        DO i=1,numplot
	           diff = lim(i,2)- lim(i,1)
	           mid = 0.5*(lim(i,1) + lim(i,2))
		   lim(i,1) = mid - 0.5*zoom*diff
		   lim(i,2) = mid + 0.5*zoom*diff
	        ENDDO
	     ENDIF	
	  ELSE
	     CALL prompt('Enter scale factor (adaptive limits)',scalemax,0.0)
	  ENDIF
	  RETURN
!------------------------------------------------------------------------
       CASE(22)
!	  show/hide plot options
	  ishowopts = .not.ishowopts
	  RETURN 	  
!------------------------------------------------------------------------
       CASE(23)
          CALL write_defaults
	  RETURN
       CASE DEFAULT
          STOP ' Error: menu action not defined'  
      
      END SELECT
      
      RETURN      
      END SUBROUTINE menu_actions
