! ----------------------------------------------------------------------
! compute exact solution for a sound wave perturbation
! ----------------------------------------------------------------------

      SUBROUTINE exact_swave(time,delta,lambda,gam,x,rho,utherm,npart)
      IMPLICIT NONE
      INTEGER, PARAMETER :: npts=100
      REAL, PARAMETER :: pi = 3.1415926536
      INTEGER, INTENT(in) :: npart
      REAL, INTENT(in) :: time, delta, lambda, gam
      REAL, DIMENSION(npart), INTENT(in) :: x,rho,utherm
      REAL, DIMENSION(npts) :: xplot, rhoplot
      REAL xmin, xmax, rhomin, rhomax, rhoav, dx
      REAL period, omega, vsoundav, uthermav
      INTEGER i,j
      
      PRINT*,' Plotting sound wave at t = ',time
      xmax = MAXVAL(x)
      xmin = MINVAL(x)
      PRINT*,' xmin, xmax, gamma = ',xmin, xmax, gam      
      
      rhomax = maxval(rho)
      rhomin = minval(rho)
      rhoav = 0.0
      uthermav = 0.0
            
      DO i=1,npart
         rhoav = rhoav + rho(i)
	 uthermav = uthermav + utherm(i)
      ENDDO
            
      rhoav = rhoav/float(npart)
      PRINT*,' Average density = ',rhoav,rhomin+0.5*(rhomax-rhomin)
      uthermav = uthermav/float(npart)
      PRINT*,' Average utherm = ',uthermav
      IF (gam.eq.1.0) THEN
         vsoundav = SQRT(2./3.*uthermav)
      ELSE
         vsoundav = SQRT(uthermav*gam*(gam-1.))
      ENDIF
      
      PRINT*,' Average sound speed = ',vsoundav
      period = 1.0	!(xmax-xmin)/vsoundav
      PRINT*,' period = ',period
      omega = 2.*pi/period
      
      dx = (xmax-xmin)/float(npts-1)
      DO i=1,npts
         xplot(i) = xmin+(i-1)*dx
	 rhoplot(i) = rhoav*(1. + delta*
     &              SIN(2.*pi/lambda*(xplot(i)-xmin) - omega*time))
      ENDDO
      
      CALL PGLINE(npts,xplot,rhoplot)
      
!      CALL plot_average(x,rho,npart)
            
      RETURN
      END SUBROUTINE exact_swave
