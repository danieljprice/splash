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
  integer i,j,its,nsteps
  real, dimension(0:npts) :: xplot,yplot
  real, intent(in) :: time,gamma,sigma
  real, intent(in) :: H0, C0, A0	! parameters for toy star
  real Aprev, A,H,C, term,const
  real radstar,dx,dt,tnow
  real func, fderiv,rhoplot,deltarho
  real gamp1,gamm1,gam1,fact,constK,omega
  logical linear

  if (norder.ge.0) linear = .true.
  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam1 = 1./gamm1
  constK = 0.25   ! this is K from P = K*rho**gamma

  if (linear) then
!---------------------------------------------------------------------------
!  linear solution not done yet
     print*,' linear solution not implemented '

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
