!--------------------------------------------------------------------------
!     program to interpolate from particle data to even grid of pixels
!     this version takes any 1D cross section through a 2D data set
!     the 1D line is specified by two points, (x1,y1) and (x2,y2)
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data 	   : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Feb 2004
!--------------------------------------------------------------------------

      SUBROUTINE interpolate2D_xsec(x,y,pmass,rho,hh,dat,npart,
     &    x1,y1,x2,y2,datsmooth,npixx)
      
      IMPLICIT NONE
      REAL, PARAMETER :: pi = 3.1415926536      
      INTEGER, INTENT(IN) :: npart,npixx
      REAL, INTENT(IN), DIMENSION(npart) :: x,y,pmass,rho,hh,dat
      REAL, INTENT(IN) :: x1,y1,x2,y2
      REAL, INTENT(OUT), DIMENSION(npixx) :: datsmooth
      
      INTEGER :: i,ipix,ipixmin,ipixmax
      REAL :: hi,hi1,h2,radkern,qq,wab,rab,const
      REAL :: term,dx,dy,xpix,ypix,pixwidth,xpixwidth,xlength
      REAL :: gradient,yintercept,aa,bb,cc,determinant,det
      REAL :: xstart,xend,ystart,yend,rstart,rend
      REAL :: tol
      LOGICAL :: xsame, ysame
!
!--check for errors in input
!
      tol = 1.e-3
      ysame = (abs(y2 - y1).lt.tol)
      xsame = (abs(x2 - x1).lt.tol)
      IF (xsame.and.ysame) THEN
         PRINT*,'error: interpolate: zero length cross section'
	 RETURN
      ENDIF
      IF (npixx.EQ.0) THEN
         PRINT*,'error: interpolate: npix = 0 '
	 RETURN
      ENDIF
      PRINT*,'oblique 1D cross section through 2D data: npix =',npixx
!
!--work out the equation of the line y = mx + c from the two points input
!
      gradient = 0.
      IF (.not.xsame) gradient = y2-y1/(x2-x1)
      yintercept = y2 - gradient*x2
      print*,'line equation: y = ',gradient,'x + ',yintercept
!
!--work out length of line and divide into pixels
!
      xlength = SQRT((x2-x1)**2 + (y2-y1)**2)
      pixwidth = xlength/REAL(npixx)
      xpixwidth = (x2 - x1)/REAL(npixx)
!      print*,' line length, pixwidth, xincrement = ',
!     &           xlength,pixwidth,xpixwidth

      datsmooth = 0.
      term = 0.
!
!--loop over particles
!      
      DO i=1,npart
!
!--set kernel related quantities
!
         hi = hh(i)
	 hi1 = 1./hi
	 h2 = hi*hi
	 radkern = 2.*hi	! radius of the smoothing kernel
         const = 10./(7.*pi*h2)	! normalisation constant
	 IF (rho(i).NE.0.) term = pmass(i)*dat(i)/rho(i) 
!
!--for each particle work out which pixels it contributes to
!  to do this we need to work out the two points at which the line
!  intersects the particle's smoothing sphere (circle).
!  given by the equation (x-xi)^2 + (y-yi)^2 = (2h)^2             
!  the x co-ordinates of the points are the solutions to a
!  quadratic with co-efficients:

	 aa = 1. + gradient**2
	 bb = 2.*gradient*(yintercept - y(i)) - 2.*x(i)
	 cc = x(i)**2 + y(i)**2 - 2.*yintercept*y(i) * yintercept**2
     &        - radkern**2
!
!--work out whether there are any real solutions and find them
!
	 determinant = bb**2 - 4.*aa*cc
	 if (determinant < 0) then
!	    print*,' particle ',i,': does not contribute ',x(i),y(i)	 
	 else
	    det = SQRT(determinant)
	    xstart = (-bb - det)/(2.*aa)
	    xend =  (-bb + det)/(2.*aa)
	    if (xstart.lt.x1) xstart = x1
	    if (xstart.gt.x2) xstart = x2
	    if (xend.gt.x2) xend = x2
	    if (xend.lt.x1) xend = x1
	    ystart = gradient*xstart + yintercept
	    yend = gradient*xend + yintercept
!	    print*,' particle ',i,': xstart, end = ',xstart,xend
!	    print*,'               : ystart, end = ',ystart,yend
	    
!
!--work out position in terms of distance (no. of pixels) along the line
!
            rstart = SQRT((xstart-x1)**2 + (ystart-y1)**2)
            rend = SQRT((xend-x1)**2 + (yend-y1)**2)
	 
            ipixmin = INT(rstart/pixwidth)
            ipixmax = INT(rend/pixwidth)
	    
            if (ipixmin.LT.1) ipixmin = 1 ! make sure they only contribute
            if (ipixmax.LT.1) ipixmax = 1
	    if (ipixmax.GT.npixx) ipixmax = npixx
	    if (ipixmin.GT.npixx) ipixmax = npixx
	    
!            print*,' ipixmin,ipixmax = ',ipixmin,ipixmax
!
!--loop over pixels, adding the contribution from this particle
!
            DO ipix = ipixmin,ipixmax
               
               xpix = x1 + (ipix)*xpixwidth - 0.5*xpixwidth
               ypix = gradient*xpix + yintercept
               dy = ypix - y(i)
               dx = xpix - x(i)
               rab = SQRT(dx**2 + dy**2)
               qq = rab*hi1
!
!--SPH kernel - standard cubic spline in 2D
!		     
               IF (qq.LT.1.0) THEN
                  wab = const*(1.-1.5*qq**2 + 0.75*qq**3)
               ELSEIF (qq.LT.2.0) THEN
                  wab = const*0.25*(2.-qq)**3
               ELSE
                  wab = 0.
               ENDIF
!
!--calculate data value at this pixel using the summation interpolant
!		  
               datsmooth(ipix) = datsmooth(ipix) + term*wab		          
               
            ENDDO
            
         ENDIF
         
      ENDDO
            
      RETURN
      
      END SUBROUTINE interpolate2D_xsec
