! ----------------------------------------------------------------------
! compute exact solution for a hydrodynamic shock
! (ie. one dimensional Riemann problem)
!
! input parameters are initial left and right states of 
! density, pressure and velocity
!
! Computes shock profile at time t
! ----------------------------------------------------------------------

      SUBROUTINE exact_shock(iplot,time,gamma,
     &           rho_L,rho_R,pr_L,pr_R,vx_L,vx_R,xmin,xmax)
      IMPLICIT NONE
      INTEGER, PARAMETER :: npts=100
      INTEGER, INTENT(IN) :: iplot
      REAL, INTENT(in) :: time,gamma,xmin,xmax
      REAL, INTENT(in) :: rho_L,rho_R,pr_L,pr_R,vx_L,vx_R
      REAL, DIMENSION(npts) :: xplot, yplot
      REAL :: dx,vsound_L, vsound_R
      INTEGER i,j
      
      PRINT*,' Plotting exact Riemann solution at t = ',time
!
!--set up grid for exact solution
!          
      dx = (xmax-xmin)/REAL(npts)
      DO i=1,npts
         xplot(i) = xmin + (i-1)*dx
      ENDDO	 
!
!--calculate various quantities
!      
      vsound_L = SQRT(gamma*pr_L/rho_L)
      vsound_R = SQRT(gamma*pr_R/rho_R)
!
!--compute desired quantity at time t
!
      
!
!--plot exact line using PGPLOT
!
      CALL PGLINE(npts,xplot,yplot)
            
      RETURN
      END SUBROUTINE exact_shock
