!---------------------------
! plot rho \propto 1/h
!
      SUBROUTINE plot_rhoh(hfact,ndim)
      IMPLICIT NONE
      INTEGER, PARAMETER :: npts = 100
      INTEGER, INTENT(IN) :: ndim
      INTEGER i
      REAL, DIMENSION(0:npts) :: xplot,yplot
      REAL, INTENT(IN) :: hfact
      REAL hfact,xmax,dx
      
c      hfact = 1.5*massp
      xmax = 3.0
      dx = (xmax-xplot(0))/float(npts)
      xplot(0) = 0.01
      
      DO i=1,npts
         xplot(i) = xplot(0)+dx*i
         yplot(i) = hfact/(xplot(i))**(1./FLOAT(ndim))
c	 print*,i,' x,y = ',xplot(i),yplot(i)
      ENDDO	      
      PRINT*,' plotting h = ',hfact,'/rho**ndim'
      
      CALL PGLINE(npts+1,xplot,yplot)    
      
      RETURN
      END
