!------------------------------------------------
! demonstration plot of all the colour schemes
!------------------------------------------------
      SUBROUTINE colour_demo
      IMPLICIT NONE
      INTEGER :: i,j,nschemes,nc
!
!--npixx should be >= ncolours in setcolours.f
!      
      INTEGER, PARAMETER :: npixx = 256
      INTEGER, PARAMETER :: npixy = npixx/10
      REAL, DIMENSION(npixx,npixy) :: sample
      REAL :: xmin,xmax,ymin,ymax,dx,dy,trans(6)
      CHARACTER(len=10) :: STRING
      
      nschemes = 10

      CALL PGBEGIN(0,'/xw',1,1)
      CALL PGPAPER(5.0,0.25/SQRT(2.))
      
      xmin = 0.0
      xmax = 1.0
      ymin = 0.0
      ymax = 0.1
      dx = (xmax-xmin)/FLOAT(npixx)
      dy = (ymax-ymin)/FLOAT(npixy)
      trans(1) = xmin
      trans(2) = dx
      trans(3) = 0.0
      trans(4) = xmin
      trans(5) = 0.0
      trans(6) = dx
      
      DO j=1,npixy
         DO i=1,npixx
            sample(i,j) = (i-1)*dx
         ENDDO
      ENDDO
      
      CALL PGENV(xmin,xmax,ymin,ymax,1,-1)
      CALL PGSCH(3.0)
      CALL PGGRAY(sample,npixx,npixy,1,npixx,1,npixy,
     &            MINVAL(sample),MAXVAL(sample),trans)
      CALL PGNUMB(1,0,0,STRING,NC)
      CALL PGMTXT('T',1.0,0.5,0.5,STRING(1:nc))     
      
      DO i=2,nschemes
         CALL PGSCH(1.0)      
         CALL PGENV(xmin,xmax,ymin,ymax,1,-1)   	 
         CALL PGSCH(3.0)
         CALL PGNUMB(i,0,0,STRING,NC)
         CALL PGMTXT('T',1.0,0.5,0.5,STRING(1:nc))     
         CALL colour_set(i)
         CALL PGIMAG(sample,npixx,npixy,1,npixx,1,npixy,
     &               MINVAL(sample),MAXVAL(sample),trans)
      ENDDO
      
      CALL PGSCH(1.0)
      CALL PGEND 
      
      END
