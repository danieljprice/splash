!------------------------------------------------------------
! plot exact solution for toystar in one dimension
!
! linear solution uses Gegenbauer/Legendre polynomials
!
! non-linear solution solves ODEs, assumes linear velocity
!
! For details see Monaghan and Price (2004) MNRAS
!------------------------------------------------------------

subroutine exact_toystar(time,gamma,H0,A0,C0,sigma,norder,iplot)
  implicit none
  integer, parameter :: npts = 100
  integer, intent(in) :: iplot,norder
  integer i,j,its,nsteps
  real, dimension(0:npts) :: xplot,yplot
  real, intent(in) :: time,gamma,sigma
  real, intent(in) :: H0, C0, A0    ! parameters for toy star
  real const ! parameter for toy star
  real Hprev, Cprev, Aprev, A2, Htemp, Ctemp, Atemp, H, C, A
  real radstar,dx,dt,tnow
  real func, fderiv,rhoplot,deltarho
  real fprevC,fprevA,fprevH,ftempC,ftempA,ftempH
  real gamp1,gamm1,gam1,fact,constK,omega
  real Gn,Pm,fnorm
  logical linear

  if (norder.ge.0) linear = .true.
  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam1 = 1./gamm1
  constK = 0.25

  if (linear) then
!---------------------------------------------------------------------------
!  linear solution uses Gegenbauer & Legendre Polynomials
!  (this is for the toy star oscillations)

     omega = sqrt(0.5*(norder+1.)*(norder+2.))
     fnorm = 2.*(norder+1)*(norder+2)/real(2.*norder + 3.)

     print*,' Plotting toy star oscills: time, norder, omega = ', &
          time,norder,omega,H0,C0,A0

     if (C0.le.0.) then 
        print*,'*** C = 0 = illegal in input'
        return
     else         
        radstar = sqrt(H0/C0)
     endif
     xplot(0) = -radstar
     dx = (radstar-xplot(0))/float(npts)

     do i=0,npts
        xplot(i) = xplot(0)+dx*i
        !         print*,i,' x,y = ',xplot(i),yplot(i)
        rhoplot = (H0 - C0*xplot(i)**2)
        if (rhoplot.le.0.) rhoplot = 0.
        deltarho = Pm(xplot(i),norder+1)*sin(omega*time)
        rhoplot = rhoplot**gam1 + 2.*omega*A0*deltarho/fnorm
        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for vx
           yplot(i) = A0*Gn(xplot(i),norder)*cos(omega*time)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     if (iplot.eq.6) then        ! plot By \propto rho
        dx = (H0**gam1)/float(npts) ! ie (rhomax - 0)/npts
        xplot(0) = 0.
        yplot(0) = sigma*xplot(0)
        do i=1,npts
           xplot(i) = xplot(0) + dx*i
           yplot(i) = sigma*(xplot(i))
        enddo
     endif

     if (iplot.eq.7) then        ! plot current point on A-C plane
        call PGPT(1,C0,A0*cos(omega*time),4)
     else                        ! plot normal exact solution line
        call PGLINE(npts+1,xplot,yplot)
     endif

!---------------------------------------------------------------------------
!  non-linear solution for the fundamental (n=1) mode
!
  else
     !
     !  solve for H, C and A given initial conditions on v, rho and the time.
     !
     nsteps = 1000*(int(time) + 1)
     its = 10

     Hprev = H0
     Cprev = C0
     Aprev = A0

     !      PRINT*,' nsteps,H,C,A in = ',nsteps,H0,C0,A0
     dt = time/nsteps

     fact = 2.*(constK + 0.5*sigma**2)*gamma*gam1

     const = (A0**2 + 1. + 2.*fact*C0*gam1)*C0**(-2./gamp1)

     print*,' Plotting toy star: time, H0, C0, A0, k = ', &
          time,Hprev,Cprev,Aprev,const

     tnow = 0.
     do i = 1,nsteps
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

        !         print*,' time = ',tnow
        !         IF ((abs(C-C0).LT.5.e-3).AND. &
        !                  (abs(A-A0).LT.5.e-3).AND.(tnow.GT.5e-3)) THEN
        !             PRINT*,'*** period, t = ',tnow,' err = ',abs(C-C0)+abs(A-A0)
        !         ENDIF

     enddo

     const = (A**2 + 1. + 2.*fact*C*gam1)*C**(-2./gamp1)

     print*,' C, A, H, k = ',C,A,H,const

     if (C.le.0.) then 
        radstar = 0.5
        stop '*** C = 0 = illegal'
     else         
        radstar = sqrt(H/C)
     endif
     xplot(0) = -radstar
     dx = (radstar-xplot(0))/float(npts)

     do i=0,npts
        xplot(i) = xplot(0)+dx*i
        !         print*,i,' x,y = ',xplot(i),yplot(i)
        rhoplot = (H - C*xplot(i)**2)
        if (rhoplot.le.0.) rhoplot = 0.
        rhoplot = rhoplot**gam1
        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for vx
           yplot(i) = A*xplot(i)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     if (iplot.eq.6) then        ! plot By \propto rho
        dx = (H**gam1)/float(npts) ! ie (rhomax - 0)/npts
        xplot(0) = 0.
        yplot(0) = sigma*xplot(0)
        do i=1,npts
           xplot(i) = xplot(0) + dx*i
           yplot(i) = sigma*(xplot(i))
        enddo
     endif

     if (iplot.eq.7) then        ! plot current point on A-C plane
        call PGPT(1,C,A,4)
     else                        ! plot normal exact solution line
        call PGLINE(npts+1,xplot,yplot)
     endif
!
!------------------------------------------------------------------------
!      
  endif

  return
end subroutine exact_toystar
!
!--function to evaluate the Gegenbauer polynomial of index n given x
!
function Gn(x,n)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x
  integer :: i
  real :: Gn,Gnminus1,Gnminus2
  real :: fnorm

  fnorm = 2.*(n+1)*(n+2)/real(2.*n + 3.)
  !      PRINT*,' fnorm = ',fnorm
  !
  !--specify first two Gegenbauer polynomials
  !
  Gnminus2 = 1.
  Gnminus1 = 3.*x
  !
  !--use recurrence relation to calculate the rest
  !
  select case (n)
  case (0) 
     Gn = Gnminus2
  case (1)
     Gn = Gnminus1
  case (2:)
     do i=2,n
        Gn = ((2*i+1)*x*Gnminus1 - (i+1)*Gnminus2)/real(i)
        Gnminus2 = Gnminus1
        Gnminus1 = Gn
     enddo
  end select

  Gn = Gn/fnorm    

end function Gn

!
!--function to calculate a Legendre Polynomial of order m
!
function Pm(x,m)
  implicit none
  integer, intent(in) :: m
  real, intent(in) :: x
  integer :: i      
  real :: Pmminus1,Pmminus2,Pm
  !
  !--specify first two Legendre polynomials
  !
  Pmminus2 = 1.
  Pmminus1 = x

  select case(m)
  case (0)
     Pm = 1.
  case (1)
     Pm = x
  case (2:)        ! use recurrence relation to calculate the rest
     do i=2,m
        Pm = ((2.*(i-1.)+1.)*x*Pmminus1 - (i-1.)*Pmminus2)/real(i)
        Pmminus2 = Pmminus1
        Pmminus1 = Pm
     enddo
  end select

end function Pm
