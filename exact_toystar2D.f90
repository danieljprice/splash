!------------------------------------------------------------
! plot exact solution for toystar in two dimensions
!
! non-linear solution solves ODEs, assumes linear velocity
!
! the solutions are all plots against radius
!
! For details see Monaghan and Price (2004), in prep.
!------------------------------------------------------------

subroutine exact_toystar2D(time,gamma,H0,A0,C0,sigma,norder,iplot)
  implicit none
  integer, parameter :: npts = 100
  integer, intent(in) :: iplot,norder
  integer :: i,j,k,its
  integer :: jmode,smode
  real, parameter :: pi = 3.1415926536
  real, dimension(0:npts) :: xplot,yplot
  real, intent(in) :: time,gamma,sigma
  real, intent(in) :: H0, C0, A0	! parameters for toy star
  real :: Aprev, A,H,C, term,const
  real :: radstar,dx,dt,tnow
  real :: func, fderiv,rhoplot,deltarho
  real :: gamp1,gamm1,gam1,fact,constK,omega
  real :: drhor
  logical linear

  if (norder.ge.0) linear = .true.
  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam1 = 1./gamm1
  constK = 0.25   ! this is K from P = K*rho**gamma

  if (linear) then
!---------------------------------------------------------------------------
!  linear solution

     jmode = norder   ! radial mode
     smode = 0        ! non-axisymmetric modes (theta)
     ! omega is the frequency of oscillation
     omega = sqrt((jmode + smode)*(jmode+smode + 2./gamm1) - smode**2)
     print*,' Plotting toy star linear solution '
     print*,' Amplitude = ',A0,' period = ',2*pi/omega,' H,C = ',H0,C0
     read*
     if (C0.le.0.) then 
        radstar = 0.5
        print*,'*** C = 0 = illegal'
        return
     else	 
        radstar = sqrt(H0/C0)
     endif
     xplot(0) = -radstar
     dx = (radstar-xplot(0))/float(npts)

     do i=0,npts
        xplot(i) = xplot(0)+dx*i
        !	 print*,i,' x,y = ',xplot(i),yplot(i)
        rhoplot = (H0 - C0*xplot(i)**2)
        if (rhoplot.le.0.) rhoplot = 0.
        deltarho = drhor(jmode,smode,xplot(i),gamma)  ! functional form of rho(r)
        print*,'deltarho = ',rhoplot,deltarho,xplot(i)
        rhoplot = (rhoplot + deltarho*A0*SIN(omega*time))**gam1

        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for vx,vy
           yplot(i) = A0*xplot(i)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     call PGLINE(npts+1,xplot,yplot)

!---------------------------------------------------------------------------
!  non-linear solution
!
  else

     !  solve for H, C and A given initial conditions on v, rho and the time.
     !
     H = H0
     C = C0
     Aprev = A0

     omega = 1.0
     const = 4.*omega**2 + 4.*A0**2 
     term = 1.-4.*omega**2/const
     if (term.le.0.) then
        print*,'error: const or omega wrong, sqrt < 0'
        return
     else
        term = sqrt(term)
        !
        !--this is the solution to the 2nd order ODE for alpha
        !
        A = COS(2.*omega*time)*term/(1. + SIN(2.*omega*time)*term)
     endif

     print*,' Plotting toy star: time, A = ',time,A

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
        !	 print*,i,' x,y = ',xplot(i),yplot(i)
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
        case(4)                 ! plot solution for vx,vy
           yplot(i) = A*xplot(i)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     call PGLINE(npts+1,xplot,yplot)
!
!------------------------------------------------------------------------
!      
  endif

  return
end subroutine exact_toystar2D

!
!--function that evaluates the polynomial for rho(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
real function drhor(j,m,rad,gamma)
  implicit none
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma
  real :: ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     drhor = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!  
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  drhor = 0.
  akprev = 1.0  ! this is a_0 which is the amplitude
  print*,'mode = ',j,' nu**2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     ak = akprev*(kprev**2 + 2.*kprev/gamm1 - freqsq)/REAL(k**2)
     print*,'coeff ',k,' = ',ak,k**2,2.*k/gamm1
     drhor = drhor + ak*rad**k
     akprev = ak
  enddo

end function drhor

!
!--function that evaluates the polynomial for v(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
!real function dvr(j,m,rad)
!  implicit none
!  integer :: j,m  ! j is the radial mode, m is the theta mode
!  real :: rad
!  
!end function dvr
