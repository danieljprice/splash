!------------------------------------------------------------
! plot exact solution for toystar in two dimensions
!
! non-linear solution solves ODEs, assumes linear velocity
!
! the solutions are all plots against radius
!
! For details see Monaghan and Price (2005), in prep.
!------------------------------------------------------------

subroutine exact_toystar2D(time,gamma,polyk,totmass, &
                           H0,A0,C0,sigma,norder,iplot)
  implicit none
  integer, intent(in) :: iplot,norder
  real, intent(in) :: time,gamma,polyk,totmass,sigma
  real, intent(in) :: H0, C0, A0        ! parameters for toy star
  real :: B0
  integer, parameter :: npts = 100
  real, parameter :: pi = 3.141592653589
  integer :: i
  integer :: jmode,smode
  real, dimension(0:npts) :: xplot,yplot
  real :: Aprev, A,H,C, term,const,omeg
  real :: radstar,dx
  real :: rhoplot,deltarho
  real :: gamp1,gamm1,gam1,constK,omega,omega2
  real :: drhor
  logical linear

  linear = (norder.ge.0)
  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'Error: no toy star solution for isothermal yet'
     return
  endif
  gam1 = 1./gamm1
  if (polyK.le.0.) then
     print*,'Error: polytropic K <= 0 on input: using 0.25 by default'
     constK = 0.25
  else
     constK = polyK  !!0.25   ! this is K from P = K*rho**gamma
  endif
  
  if (linear) then
!---------------------------------------------------------------------------
!  linear solution

     print*,' Plotting 2D toy star: linear solution '
     jmode = norder   ! radial mode
     smode = 0        ! non-axisymmetric modes (theta)
     ! omega is the frequency of oscillation
     omega2 = (jmode + smode)*(jmode+smode + 2./gamm1) - smode**2
     if (omega2.le.0.) then
        print*,'Error: sqrt < 0 in linear toy star  ',omega2
        return
     else
        omega = sqrt(omega2)
     endif
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
        !         print*,i,' x,y = ',xplot(i),yplot(i)
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

     print*,'Plotting 2D toy star: non-linear'
     !  solve for H, C and A given initial conditions on v, rho and the time.
     !
     H = H0
     C = C0
     Aprev = A0
     B0 = 0.
!
!--this is the static solution, determined from the total mass, polyk, gamma and omega
!
     omeg = 1.0  ! this is omega from the main code (ie. from potential)
     
     radstar = sqrt(gamma*totmass/(pi*gamm1))
     H = omeg**2*gamm1*radstar**2/(2.*polyk*gamma)
     C = 0.5*gamm1*omeg**2/(gamma*polyk)
     print*,'r_star = ',radstar,' rho = (',H,'-',C,'^2)**',gamm1
!
!--work out period of oscillation
!
     omega = 4.*(B0**2 + C*polyk*gamma**2/gamm1)
     if (omega.le.1.e-5) then
        print*,'ERROR: sqrt < 0 in omega'
        return
     else
        omega = sqrt(omega)
     endif
     print*,'period = ',2.*pi/omega
!
!--solve for alpha(t)
!    
     const = 4.*omega**2 + 4.*A0**2 
     term = 1.-4.*omega**2/const
     if (term.le.0.) then
        if (abs(A0).gt.1.e-3) print*,'warning: const or omega wrong, sqrt < 0 : assuming static solution'
        A = 0.
     else
        term = sqrt(term)
        !
        !--this is the solution to the 2nd order ODE for alpha
        !
        A = omega*COS(2.*omega*time)*term/(1. + SIN(2.*omega*time)*term)
     endif

     print*,' Plotting toy star: time, A = ',time,A

     !if (C.le.0.) then 
     !   radstar = 0.5
     !   stop '*** C = 0 = illegal'
     !elseif (A.le.1.e-5) then
     !else
     !   radstar = sqrt(H/C)
     !endif
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
        case(4)                 ! plot solution for vx,vy
           yplot(i) = A*xplot(i)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     if (iplot.gt.0 .and. iplot.le.5) then
        call pgline(npts+1,xplot,yplot)
     elseif (iplot.eq.0) then
        call pgsfs(2)
        call pgcirc(0.0,0.0,radstar)
     endif
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
