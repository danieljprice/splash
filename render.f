!------------------------------------------------------------------------
!  This subroutine takes a 2D grid of data and renders it using PGPLOT
!  Rendering is either greyscale (icolours = 1) or colour (icolours>1)
!  Also plots nc contours between datmin and datmax.
!------------------------------------------------------------------------
      
      SUBROUTINE render(datpix,datmin,datmax,label,npixx,npixy,
     &                  xmin,ymin,dx,icolours,iplotcont,nc,log)
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npixx,npixy,nc,icolours
      REAL, INTENT(IN) :: xmin,ymin,datmin,datmax,dx
      REAL, DIMENSION(npixx,npixy), INTENT(IN) :: datpix
      LOGICAL, INTENT(IN) :: iplotcont,log
      CHARACTER(LEN=*), INTENT(IN) :: label            
      
      INTEGER i,j,k
      REAL :: trans(6),levels(nc),dcont
      CHARACTER(LEN=1) :: clog
      
!     set up grid for rendering      

      trans(1) = xmin - 0.5*dx		! this is for the PGIMAG call
      trans(2) = dx			! see help for PGIMAG/PGGRAY/PGCONT
      trans(3) = 0.0
      trans(4) = ymin - 0.5*dx
      trans(5) = 0.0
      trans(6) = dx

!     set character to send to PGWEDG call if log (danpgwedg only)      
      clog = ' '
      IF (log) clog = 'L'
     
      PRINT*,'rendering...',npixx,'x',npixy,',array size=',SIZE(datpix)
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
         CALL DANPGWEDG('RGV'//clog,0.5,4.5,datmin,datmax,label)
         CALL PGGRAY(datpix,npixx,npixy,1,npixx,1,npixy,
     &                                   datmin,datmax,trans)
     
      ELSEIF (icolours.GT.1) THEN	! colour
         CALL DANPGWEDG('RIV'//clog,0.5,4.5,datmin,datmax,label)
!         CALL PGWEDG('RI',2.0,4.0,datmin,datmax,' ')
!         CALL PGPIXL(datpix,npixx,npixx,1,npixx,1,npixx,xmin,xmax,ymin,ymax)
         CALL PGIMAG(datpix,npixx,npixy,1,npixx,1,npixy,
     &                                datmin,datmax,trans)
!         CALL PGHI2D(datpix,npixx,npixx,1,npixx,1,npixy,
!     &	             x,1,0.1,.true.,y)
	 
      ENDIF
!
!--plot contours
!
      IF (iplotcont) THEN
         PRINT*,'plotting ',nc,' contours...'
     	 CALL PGCONT(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans)
      ENDIF
      
      RETURN
           
      END
