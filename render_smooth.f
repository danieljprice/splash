!------------------------------------------------------------------------
!  This subroutine calls the interpolation routine to smooth data
!  from the SPH particles onto a uniform grid and then calls
!  the appropriate PGPLOT routines to plot the pixel map
!
!  scalar field dat is interpolated using the SPH summation interpolant
!
!------------------------------------------------------------------------
      
      SUBROUTINE smooth_render(x,y,xminpart,xmaxpart,yminpart,ymaxpart,
     &                         dat,datmin,datmax,label,
     &                         npart,ntot,npixx,icolours,iplotcont,
     &                         nc,pmass,rho,hh)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npart,ntot,npixx,icolours,nc
      INTEGER i,j,k,npixy
      REAL, INTENT(IN), DIMENSION(ntot) :: x,y,pmass,rho,hh,dat
      REAL, INTENT(IN) :: xminpart,xmaxpart,yminpart,ymaxpart
      LOGICAL, INTENT(IN) :: iplotcont
      CHARACTER(LEN=20), INTENT(IN) :: label
            
      REAL, DIMENSION(npixx,npixx) :: datpix      
      REAL :: xmin,xmax,ymin,ymax
      REAL dx,scale,zero
      REAL trans(6),levels(nc),dcont,datmin,datmax
      
      zero= 0.0E0
!     set up grid for rendering
      
      xmax = xmaxpart + 0.000001
      xmin = xminpart - 0.000001
      ymax = ymaxpart + 0.000001
      ymin = yminpart - 0.000001
      dx = (xmax-xmin)/REAL(npixx)	! this is pixel width

      trans(1) = xmin-0.5*dx		! this is for the PGIMAG call
      trans(2) = dx			! see help for PGIMAG/PGGRAY/PGCONT
      trans(3) = 0.0
      trans(4) = xmin-0.5*dx
      trans(5) = 0.0
      trans(6) = dx

      npixy = INT((ymax-ymin)/dx) + 1	! pixels in y direction
      IF (npixy.GT.npixx) npixy = npixx
      PRINT*,'npixx,npixy = ',npixx,npixy
!
!--if doing scalar field rendering, interpolate to finer grid
!      
      CALL interpolate2D(x,y,pmass,rho,hh,dat,ntot,
     &                   xmin,ymin,datpix,npixx,npixx,dx)	 
	       
      PRINT*,'rendering...'
!
!--set contour levels
!      
      dcont = (datmax-datmin)/REAL(nc+1)   ! even contour levels
      DO i=1,nc
         levels(i) = datmin + REAL(i)*dcont
      ENDDO
!
!--nb: plots use my modification of PGWEDG which plots vertical numbers on axes
!	 
      IF (icolours.EQ.1) THEN	! greyscale
         CALL DANPGWEDG('RG',0.5,4.5,datmin,datmax,label)
         CALL PGGRAY(datpix,npixx,npixx,1,npixx,1,npixy,
     &                                   datmin,datmax,trans)
     
      ELSEIF (icolours.GT.1) THEN	! colour
         CALL DANPGWEDG('RI',0.5,4.5,datmin,datmax,label)
!         CALL PGWEDG('RI',2.0,4.0,datmin,datmax,' ')
!         CALL PGPIXL(datpix,npixx,npixx,1,npixx,1,npixx,xmin,xmax,ymin,ymax)
         CALL PGIMAG(datpix,npixx,npixx,1,npixx,1,npixy,
     &                                datmin,datmax,trans)
!         CALL PGHI2D(datpix,npixx,npixx,1,npixx,1,npixy,
!     &	             x,1,0.1,.true.,y)
	 
      ENDIF
!
!--plot contours
!
      IF (iplotcont)
     &	 CALL PGCONT(datpix,npixx,npixx,1,npixx,1,npixy,levels,nc,trans)
     
      RETURN
           
      END
