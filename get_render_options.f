!-----------------------------------------------------------------------------
! options for rendering / vector plots
!-----------------------------------------------------------------------------

      SUBROUTINE get_render_options(irender,npix,icolours,iplotcont,
     &	        ncontours,ivecplot,npixvec,iplotpartvec,
     &          x_sec,xsecpos,backgnd_vec,ndim,numplot)     
      USE prompting
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim, numplot
      INTEGER, INTENT(OUT) :: irender,npix,icolours,ivecplot,npixvec
      INTEGER, INTENT(OUT) :: ncontours
      LOGICAL, INTENT(OUT) :: iplotcont,iplotpartvec,x_sec,backgnd_vec
      REAL, INTENT(OUT) :: xsecpos
      CHARACTER(LEN=1) :: ans      
!
!--rendering options
!
      CALL prompt(' Enter field (from menu) for rendering (0=none)',
     &             irender,0,numplot)
      IF ((irender.GT.ndim)
     &	 .AND.(irender.LE.numplot)) THEN 
      
         CALL prompt('Enter number of pixels along x axis',npix,1,1000)
      
100      CONTINUE
         PRINT*,'(-ve = demo, 0 = contours only)'
	 CALL prompt('Enter colour scheme for rendering ',icolours,max=5)
	 IF (icolours.LT.0) THEN
	    CALL colour_demo
	    GOTO 100
	 ENDIF   
         
	 IF (ndim.EQ.3) THEN
	    PRINT*,' Take cross section (c) '//
     &	           'or projection (p) through 3D data?'
	    READ*,ans
	    x_sec = .false.
	    IF (ans.EQ.'c') x_sec = .true.
	    PRINT*,' Cross section = ',x_sec
	    xsecpos = 0.0
	    IF (x_sec) THEN
	       CALL prompt('Enter co-ordinate location of cross section slice',
     &                     xsecpos)	       
	    ENDIF   
	 ENDIF
	 
         IF (icolours.NE.0) THEN
	    CALL prompt(' Plot contours?',iplotcont)
         ELSE
	    iplotcont = .true.
         ENDIF
 
         IF (iplotcont) THEN
	    CALL prompt(' Enter number of contours between min,max',
     &                  ncontours,1,100)
         ENDIF
      ELSE
	 irender = 0   
      ENDIF
!
!--vector plot options
!	  
      CALL prompt('vector plot: '//
     &  'velocity (1), magnetic field (2) or none(0)?',ivecplot,0,2)
      IF (ivecplot.GT.0) THEN	
	 CALL prompt('Enter number of pixels',npixvec,1,1000)
         iplotpartvec = .false.
         IF (irender.EQ.0) THEN
	    CALL prompt(' Plot particles?',iplotpartvec)
	 ENDIF
	 backgnd_vec = .false.
	 CALL prompt('Use background color for vector plot?',backgnd_vec)
      ENDIF
    
      RETURN
      END SUBROUTINE get_render_options
