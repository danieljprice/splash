! ------------------------------------------------------------------------
!      defines colour schemes for rendering
!      ** add your own here **
! ------------------------------------------------------------------------     
      SUBROUTINE colour_set(icolours)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: icolours
      INTEGER i,icolourmin,icolmin,icolmax,ncolmax
      INTEGER ncolours
      REAL red,green,blue,hue,light,sat
      REAL dhue,dlight,dsat,dred,dgreen,dblue
      LOGICAL rgb
      
      red = 0.0
      green = 0.0
      blue = 0.0
      rgb = .false.
!
!--set number of colours
!      
      ncolours = 256
!
!--set first colour index (warning: colours 1-16 have presets, so
!  overwriting these means that line graphs that use colour will come
!  out funny). Best to leave 0 and 1 alone as these are black and white.
!
      icolourmin = 2
!
!--inquire as to colour range available on current device
!  adjust ncolours if necessary
!      
      CALL PGQCOL(icolmin,icolmax)
      IF (icolourmin.lt.icolmin) icolourmin = icolmin
      ncolmax = icolmax - icolourmin
      IF (ncolours.gt.ncolmax) THEN  
         ncolours = ncolmax
         PRINT*,'Warning: Device allows only ',ncolours+1,' colours'
      ENDIF
!
!--starting values
!      
      IF (icolours.EQ.2) THEN	! red
         rgb = .false.
         hue=100
         light=0.0
         sat=1.0      
	 dlight = 1.0/FLOAT(ncolours)
	 dsat = 0.0
	 dhue = 80.0/FLOAT(ncolours)	 
      ELSEIF (icolours.EQ.3) THEN	! blue
         rgb = .false.
         hue=300
         light=0.
         sat=1.0
	 dlight = 1.0/FLOAT(ncolours)
	 dsat = 0.0
	 dhue = 40.0/FLOAT(ncolours)
      ELSEIF (icolours.EQ.4) THEN	! rainbow
         rgb = .false.         
         hue=100
         light=0.5
         sat=1.0
	 dlight = 0.0	!/FLOAT(ncolours)
	 dsat = 0.0
	 dhue = 320.0/FLOAT(ncolours)
      ELSEIF (icolours.EQ.5) THEN	! green
         rgb = .false.         
         hue=200
         light=0.0
         sat=1.0
	 dlight = 1.0/FLOAT(ncolours)
	 dsat = 0.0
	 dhue = 80.0/FLOAT(ncolours)
      ELSEIF (icolours.EQ.6) THEN
         rgb = .false.
         hue=120
         light=0.0
         sat= 0.    
	 dlight = 1.0/FLOAT(ncolours)
	 dsat = 0.	!1.0/FLOAT(ncolours)
	 dhue = 80.0/FLOAT(ncolours)	 
      ELSEIF (icolours.EQ.10) THEN	! rgb attempts
         rgb = .true.
         red = 0.3333
         green = 0.0
         blue = 0.0      
	 dred = 0.0	!0.333/FLOAT(ncolours)
	 dblue = 0.3333/FLOAT(ncolours)
	 dgreen = 0.0	!0.333/FLOAT(ncolours)	 
      ENDIF
!
!--increment values in steps
!
      DO i=icolourmin,icolourmin+ncolours	! set colour indexes from 17-> 17+ncolours
!         red = red + 0.05
!         green = green + 0.015
	 
	 IF (rgb) THEN
	    red = red + dred
	    blue = blue + dblue
	    green = green + dgreen
	    CALL PGSCR(i,red,green,blue)	 
	 ELSE   
	    hue = hue + dhue
	    sat = sat + dsat
	    light = light + dlight
	    CALL PGSHLS(i,hue,light,sat)
         ENDIF
      ENDDO
!
!--set this as the range of colour indices to use
!
      CALL PGSCIR(icolourmin,icolourmin+ncolours)  
      
      END
