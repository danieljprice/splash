!-------------------------------------------------------------------
!  This subroutine takes an array of particle positions in x,y
!  and a scalar or 2D vector array 'dat' defined on the particles
!  and interpolates to a grid which can then be rendered.
!
!  scalar field is interpolated to a finer grid using SPH summation
!
!------------------------------------------------------------------
      
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
            
      REAL, DIMENSION(npixx,npixx) :: xpix,ypix,datpix      
      REAL :: xmin,xmax,ymin,ymax
      REAL dx,scale,zero
      REAL trans(6),levels(nc),dcont,datmin,datmax
      
      INTEGER, PARAMETER :: ncol=10
      REAL :: eye(3), lattice(3,3), light(3), rgb(3,ncol)
      REAL :: dred,dgreen,dblue
      
      
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
      PRINT*,' npixx,npixy = ',npixx,npixy
      ! need to send in grid as an array
      DO j=1,npixy
	 DO i=1,npixx
	    xpix(i,j) = xmin + (i-1)*dx - 0.5*dx
	    ypix(i,j) = ymin + (j-1)*dx - 0.5*dx
	 ENDDO
      ENDDO   
!
!--if doing scalar field rendering, interpolate to finer grid
!      
      IF (ANY (rho.EQ.0.)) THEN
         PRINT*,'Error: density zero somewhere, can''t render'
      ELSE
         CALL smooth_pixels(x,y,pmass,rho,hh,dat,npart,ntot,
     &                   xpix,ypix,datpix,npixx,npixy,
     &                   xmin,xmax,ymin,ymax)
      ENDIF
      
      PRINT*,'rendering...'
      dcont = (datmax-datmin)/REAL(nc+1)   ! even contour levels
      DO i=1,nc
         levels(i) = datmin + REAL(i)*dcont
      ENDDO
      
      DATA eye /1.0, 1.0, 3.0/
      DATA light /-1.0, -1.0, -1.0/
      lattice(1,1:3) = 0.0
      lattice(2,1:3) = -1.0
      lattice(3,1:3) = -0.5
      
      rgb(1,1) = 0.0
      rgb(2,1) = 0.0
      rgb(3,1) = 0.0
      dred = 0.1
      dgreen = 0.01
      dblue = 0.1
            
      DO i=2,ncol
         rgb(1,i) = rgb(1,i-1) + dred
         rgb(2,i) = rgb(2,i-1) + dgreen	 
	 rgb(3,i) = rgb(3,i-1) + dblue
      ENDDO	 
      
      CALL COLSRF(RGB,ncol,1.0,1,10,2,0.5,0.5,1.0)
      CALL SB2SRF(eye,lattice,datpix,npixx,npixx,
     &             datmin,datmax,datmax,1.0,10,1,light,.false.)
!!      CALL SBSURF(eye,lattice,dens
      RETURN
      
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
!         CALL PGPIXL(datpix,npix,npix,1,npix,1,npix,xmin,xmax,ymin,ymax)
         CALL PGIMAG(datpix,npixx,npixx,1,npixx,1,npixy,
     &                                datmin,datmax,trans)

!         CALL PGHI2D(datpix,npixx,npixx,1,npixx,1,npixy,
!     &	             x,1,0.1,.true.,y)
	 
      ENDIF

      IF (iplotcont)
     &	 CALL PGCONT(datpix,npixx,npixx,1,npixx,1,npixy,levels,nc,trans)
     
      RETURN
           
      END
