! ------------------------------------------------------------------------
!      defines colour schemes for rendering
!      ** add your own here **
! ------------------------------------------------------------------------     
      SUBROUTINE setcolours(icolours)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: icolours
      INTEGER i,nlevels
      REAL red,green,blue,hue,light,sat
      REAL dhue,dlight,dsat,dred,dgreen,dblue
      LOGICAL rgb
      
      red = 0.0
      green = 0.0
      blue = 0.0
      rgb = .false.
      nlevels = 50	! number of colour levels
!
!--starting values
!      
      IF (icolours.EQ.2) THEN	! red
         rgb = .false.
         hue=100
         light=0.0
         sat=1.0      
	 dlight = 1.0/FLOAT(nlevels)
	 dsat = 0.0
	 dhue = 80.0/FLOAT(nlevels)	 
      ELSEIF (icolours.EQ.3) THEN	! blue
         rgb = .false.
         hue=300
         light=0.
         sat=1.0
	 dlight = 1.0/FLOAT(nlevels)
	 dsat = 0.0
	 dhue = 40.0/FLOAT(nlevels)
      ELSEIF (icolours.EQ.4) THEN	! rainbow
         rgb = .false.         
         hue=100
         light=0.5
         sat=1.0
	 dlight = 0.0	!/FLOAT(nlevels)
	 dsat = 0.0
	 dhue = 320.0/FLOAT(nlevels)
      ELSEIF (icolours.EQ.5) THEN	! green
         rgb = .false.         
         hue=200
         light=0.0
         sat=1.0
	 dlight = 1.0/FLOAT(nlevels)
	 dsat = 0.0
	 dhue = 80.0/FLOAT(nlevels)
      ELSEIF (icolours.EQ.6) THEN
         rgb = .false.
         hue=120
         light=0.0
         sat= 0.    
	 dlight = 1.0/FLOAT(nlevels)
	 dsat = 0.	!1.0/FLOAT(nlevels)
	 dhue = 80.0/FLOAT(nlevels)	 
      ELSEIF (icolours.EQ.10) THEN	! rgb attempts
         rgb = .true.
         red = 0.3333
         green = 0.0
         blue = 0.0      
	 dred = 0.0	!0.333/FLOAT(nlevels)
	 dblue = 0.3333/FLOAT(nlevels)
	 dgreen = 0.0	!0.333/FLOAT(nlevels)	 
      ENDIF
!
!--increment values in steps
!
      DO i=20,20+nlevels	! set colour indexes from 20-> 20+nlevels
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
      CALL PGSCIR(20,20+nlevels)  
      
      END
