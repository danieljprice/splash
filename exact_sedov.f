!---------------------------------------------------------------------------
! compute exact solution for Sedov-type point-like energy injection
!---------------------------------------------------------------------------

      SUBROUTINE exact_sedov(time,gam,r,rho,npart)
      IMPLICIT NONE
      INTEGER, PARAMETER :: npts=100
      REAL, PARAMETER :: pi = 3.1415926536
      INTEGER, INTENT(in) :: npart
      REAL, INTENT(in) :: time, gam
      REAL, DIMENSION(npart), INTENT(in) :: r,rho
      REAL, DIMENSION(0:npts) :: rplot, rhoplot
      REAL rmax, rshock, rhomin, rhomax, rhozero, dr
      REAL beta, energy
      INTEGER i,j,ishock
      
      PRINT*,' Plotting similarity solution at t = ',time
      PRINT*,' Valid for gamma = 5./3., we have gamma = ',gam      
      
      beta = 1.1517
      energy = 1.0
      rhozero = 1.0
      rmax = MAXVAL(r)
            
      dr = rmax/float(npts-1)
      rshock = beta*(energy*time**2./rhozero)**(0.2)
      
      rplot(0) = 0.0
      ishock = npts+1
      
      DO i=1,npts
         rplot(i) = i*dr
	 IF (i.LT.ishock) rhoplot(i) = rhozero*rplot(i)
	 IF (rplot(i).GT.rshock.AND.rplot(i-1).LE.rshock) THEN
	    ishock = i
	    rhoplot(i) = rhozero*(gam+1.)/(gam-1.)
	 ENDIF
	 IF (i.GT.ishock) rhoplot(i) = rhozero
      ENDDO
      
      CALL PGLINE(npts+1,rplot,rhoplot)
      
      RETURN
      END SUBROUTINE exact_sedov
