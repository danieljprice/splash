!------------------------------------------------------------
! plot exact solution for toystar
!
! linear solution uses Gegenbauer/Legendre polynomials
!
! non-linear solution solves ODEs, assumes linear velocity
!
! For details see Monaghan and Price (2003) MNRAS
!------------------------------------------------------------

      SUBROUTINE exact_toystar(time,gamma,H0,A0,C0,sigma,
     &                         norder,iplot)
      IMPLICIT NONE
      INTEGER, PARAMETER :: npts = 100
      INTEGER, INTENT(IN) :: iplot,norder
      INTEGER i,j,its,nsteps
      REAL, DIMENSION(0:npts) :: xplot,yplot
      REAL, INTENT(IN) :: time,gamma,sigma
      REAL, INTENT(IN) :: H0, C0, A0	! parameters for toy star
      REAL const	! parameter for toy star
      REAL Hprev, Cprev, Aprev, A2, Htemp, Ctemp, Atemp, H, C, A
      REAL radstar,dx,dt,tnow
      REAL func, fderiv,rhoplot,deltarho
      REAL fprevC,fprevA,fprevH,ftempC,ftempA,ftempH
      REAL gamp1,gamm1,gam1,fact,constK,omega
      REAL Gn,Pm,fnorm
      LOGICAL linear

      IF (norder.GE.0) linear = .true.
!
      gamp1 = gamma + 1.
      gamm1 = gamma - 1.
      gam1 = 1./gamm1
      constK = 0.25

      IF (linear) THEN
!---------------------------------------------------------------------------
!--linear solution uses Gegenbauer & Legendre Polynomials
!  (this is for the toy star oscillations)

      omega = SQRT(0.5*(norder+1.)*(norder+2.))
      fnorm = 2.*(norder+1)*(norder+2)/REAL(2.*norder + 3.)

      PRINT*,' Plotting toy star oscills: time, norder, omega = ',
     &       time,norder,omega,H0,C0,A0

      IF (C0.LE.0.) THEN 
	 PRINT*,'*** C = 0 = illegal'
	 RETURN
      ELSE	 
	 radstar = SQRT(H0/C0)
      ENDIF	 
      xplot(0) = -radstar
      dx = (radstar-xplot(0))/float(npts)

      DO i=0,npts
         xplot(i) = xplot(0)+dx*i
c	 print*,i,' x,y = ',xplot(i),yplot(i)
         rhoplot = (H0 - C0*xplot(i)**2)
         IF (rhoplot.LE.0.) rhoplot = 0.
	 deltarho = Pm(xplot(i),norder+1)*SIN(omega*time)
         rhoplot = rhoplot**gam1 + 2.*omega*A0*deltarho/fnorm
         IF (iplot.EQ.1) THEN	! plot solution for density
            yplot(i) = rhoplot
         ELSEIF (iplot.EQ.2) THEN	! plot solution for pressure
            yplot(i) = constK*rhoplot**gamma
         ELSEIF (iplot.EQ.3) THEN  ! plot solution for utherm
            yplot(i) = constK*(rhoplot**gamm1)/gamm1
         ELSEIF (iplot.EQ.4) THEN  ! plot solution for vx
            yplot(i) = A0*Gn(xplot(i),norder)*COS(omega*time)
	 ELSEIF (iplot.EQ.5) THEN ! plot solution for By
	    yplot(i) = sigma*rhoplot
         ENDIF
      
      ENDDO      
      
      IF (iplot.EQ.6) THEN	! plot By \propto rho
         dx = (H0**gam1)/float(npts) ! ie (rhomax - 0)/npts
	 xplot(0) = 0.
	 yplot(0) = sigma*xplot(0)
	 DO i=1,npts
	    xplot(i) = xplot(0) + dx*i
	    yplot(i) = sigma*(xplot(i))
	 ENDDO
      ENDIF

      IF (iplot.EQ.7) THEN	! plot current point on A-C plane
         CALL PGPT(1,C0,A0*COS(omega*time),4)
      ELSE			! plot normal exact solution line
         CALL PGLINE(npts+1,xplot,yplot)
      ENDIF
       
!---------------------------------------------------------------------------
!  non-linear solution
!
      ELSE
      
!  solve for H, C and A given initial conditions on v, rho and the time.
!
      nsteps = 1000*(INT(time) + 1)
      its = 10
      
      Hprev = H0
      Cprev = C0
      Aprev = A0

!      PRINT*,' nsteps,H,C,A in = ',nsteps,H0,C0,A0
      dt = time/nsteps

      fact = 2.*(constK + 0.5*sigma**2)*gamma*gam1

      const = (A0**2 + 1. + 2.*fact*C0*gam1)*C0**(-2./gamp1)
      
      PRINT*,' Plotting toy star: time, H0, C0, A0, k = ',
     &       time,Hprev,Cprev,Aprev,const
     
      tnow = 0.
      DO i = 1,nsteps
         tnow = tnow + dt
! integrate using improved Euler         
	 fprevC = -Cprev*Aprev*gamp1
	 fprevA = fact*Cprev -1.-Aprev**2
	 fprevH = -Aprev*Hprev*gamm1
! predictor	 
	 Ctemp = Cprev + dt*(-Cprev*Aprev*gamp1)
	 Atemp = Aprev + dt*(fact*Cprev-1.-Aprev**2)
	 Htemp = Hprev + dt*(-Aprev*Hprev*gamm1)
	 
	 ftempC = -Ctemp*Atemp*gamp1
	 ftempA = fact*Ctemp -1. -Aprev**2
	 ftempH = -Atemp*Htemp*gamm1
! corrector	 
	 C = Cprev + 0.5*dt*(fprevC + ftempC)
	 A = Aprev + 0.5*dt*(fprevA + ftempA)
	 H = Hprev + 0.5*dt*(fprevH + ftempH)
	 
         Cprev = C
	 Aprev = A
	 Hprev = H
	 
!	 print*,' time = ',tnow
!	 IF ((abs(C-C0).LT.5.e-3).AND.
!     &	     (abs(A-A0).LT.5.e-3).AND.(tnow.GT.5e-3)) THEN
!	     PRINT*,'*** period, t = ',tnow,' err = ',abs(C-C0)+abs(A-A0)
!	 ENDIF
	 
      ENDDO
      
      const = (A**2 + 1. + 2.*fact*C*gam1)*C**(-2./gamp1)

      print*,' C, A, H, k = ',C,A,H,const

      IF (C.LE.0.) THEN 
         radstar = 0.5
	 STOP '*** C = 0 = illegal'
      ELSE	 
	 radstar = SQRT(H/C)
      ENDIF	 
      xplot(0) = -radstar
      dx = (radstar-xplot(0))/float(npts)
            
      DO i=0,npts
         xplot(i) = xplot(0)+dx*i
c	 print*,i,' x,y = ',xplot(i),yplot(i)
         rhoplot = (H - C*xplot(i)**2)
         IF (rhoplot.LE.0.) rhoplot = 0.
         rhoplot = rhoplot**gam1
         IF (iplot.EQ.1) THEN	! plot solution for density
            yplot(i) = rhoplot
         ELSEIF (iplot.EQ.2) THEN	! plot solution for pressure
            yplot(i) = constK*rhoplot**gamma
         ELSEIF (iplot.EQ.3) THEN  ! plot solution for utherm
            yplot(i) = constK*(rhoplot**gamm1)/gamm1
         ELSEIF (iplot.EQ.4) THEN  ! plot solution for vx
	    yplot(i) = A*xplot(i)
	 ELSEIF (iplot.EQ.5) THEN ! plot solution for By
	    yplot(i) = sigma*rhoplot
         ENDIF
      
      ENDDO      
      
      IF (iplot.EQ.6) THEN	! plot By \propto rho
         dx = (H**gam1)/float(npts) ! ie (rhomax - 0)/npts
	 xplot(0) = 0.
	 yplot(0) = sigma*xplot(0)
	 DO i=1,npts
	    xplot(i) = xplot(0) + dx*i
	    yplot(i) = sigma*(xplot(i))
	 ENDDO
      ENDIF
      
      IF (iplot.EQ.7) THEN	! plot current point on A-C plane
         CALL PGPT(1,C,A,4)
      ELSE			! plot normal exact solution line
         CALL PGLINE(npts+1,xplot,yplot)
      ENDIF
!
!------------------------------------------------------------------------
!      
      ENDIF	 
      
      RETURN
      END SUBROUTINE exact_toystar
 !
!--function to evaluate the Gegenbauer polynomial of index n given x
!
      FUNCTION Gn(x,n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: x
      INTEGER :: i
      REAL :: Gn,Gnminus1,Gnminus2
      REAL :: fnorm
 
      fnorm = 2.*(n+1)*(n+2)/REAL(2.*n + 3.)
!      PRINT*,' fnorm = ',fnorm
!
!--specify first two Gegenbauer polynomials
!
      Gnminus2 = 1.
      Gnminus1 = 3.*x
!
!--use recurrence relation to calculate the rest
!
      SELECT CASE (n)
         CASE (0) 
            Gn = Gnminus2
         CASE (1)
            Gn = Gnminus1
         CASE (2:)
	    DO i=2,n
	       Gn = ((2*i+1)*x*Gnminus1 - (i+1)*Gnminus2)/REAL(i)
	       Gnminus2 = Gnminus1
	       Gnminus1 = Gn
	    ENDDO
      END SELECT
            
      Gn = Gn/fnorm    
	    
      END FUNCTION Gn     

!
!--function to calculate a Legendre Polynomial of order m
!
      FUNCTION Pm(x,m)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL, INTENT(IN) :: x
      INTEGER :: i      
      REAL :: Pmminus1,Pmminus2,Pm
!
!--specify first two Legendre polynomials
!
      Pmminus2 = 1.
      Pmminus1 = x
 
      SELECT CASE(m)
         CASE (0)
            Pm = 1.
         CASE (1)
            Pm = x
         CASE (2:)	! use recurrence relation to calculate the rest
            DO i=2,m
               Pm = ((2.*(i-1.)+1.)*x*Pmminus1 
     &	                     - (i-1.)*Pmminus2)/REAL(i)
    	       Pmminus2 = Pmminus1
   	       Pmminus1 = Pm
            ENDDO
      END SELECT
  
      END FUNCTION Pm      
