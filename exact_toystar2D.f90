!------------------------------------------------------------
! plot exact solution for toystar in two dimensions
!
! non-linear solution solves ODEs, assumes linear velocity
!
! the solutions are all plots against radius
!
! For details see Monaghan and Price (2005), in prep.
!
!
! iplot = 0 gives x vs y
!
! iplot = 1->5 gives rho, pr, u, vx, vy vs r
!------------------------------------------------------------

module toystar2D
 implicit none

contains

subroutine exact_toystar2D(time,gamma,polyk,totmass, &
                           H0,A0,C0,Brhofac,jorder,morder,iplot)
  use toystar2D_utils
  implicit none
  integer, intent(in) :: iplot,jorder,morder
  real, intent(in) :: time,gamma,polyk,totmass,Brhofac
  real, intent(in) :: H0, C0, A0        ! parameters for toy star
  real :: B0
  integer, parameter :: npts = 100
  real, parameter :: pi = 3.141592653589
  integer :: i
  integer :: jmode,smode
  real, dimension(0:npts) :: xplot,yplot
  real :: Aprev, A,H,C, term,const,omega,omegasq
  real :: radstar,dx,nu2
  real :: rhoplot,deltarho,vplot,deltav
  real :: gamp1,gamm1,gam1,constK,sigma
  logical linear

  linear = (jorder.ge.0 .or. morder.ge.0)
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

  omega = 1.0  ! this is omega from the main code (ie. from potential)
  omegasq = omega**2
  
  if (linear) then
!---------------------------------------------------------------------------
!  linear solution

     print*,' Plotting 2D toy star: linear solution '
     jmode = jorder   ! radial mode
     smode = morder        ! non-axisymmetric modes (theta)
     
     ! sigma is the frequency of oscillation
     nu2 = (jmode + smode)*(jmode+smode + 2./gamm1) - smode**2
     if (nu2.le.0.) then
        print*,'Error: nu^2 < 0 in linear toy star  ',nu2
        print*,' radial mode = ',jmode,' theta mode = ',smode
        return
     else
        sigma = sqrt(0.5*omegasq*gamm1*nu2)
     endif
     print*,' Amplitude = ',A0,' period = ',2*pi/sigma,' H,C = ',H0,C0

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
        deltarho = etar(jmode,smode,xplot(i)/radstar,gamma)  ! functional form of rho(r)
        print*,'deltarho = ',rhoplot,deltarho,xplot(i)
        rhoplot = (rhoplot + deltarho*A0*SIN(sigma*time))**gam1
        
        vplot = deltav*COS(sigma*time)

        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for vx,vy
           yplot(i) = vplot
        case(5)                 ! plot solution for By
           yplot(i) = Brhofac*rhoplot
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
     
     radstar = sqrt(gamma*totmass/(pi*gamm1))
     H = omegasq*gamm1*radstar**2/(2.*polyk*gamma)
     C = 0.5*gamm1*omegasq/(gamma*polyk)
     print*,'r_star = ',radstar,' rho = (',H,'-',C,'^2)**',gamm1
!
!--work out period of oscillation
!
     sigma = 4.*(B0**2 + C*polyk*gamma**2/gamm1)
     if (sigma.le.1.e-5) then
        print*,'ERROR: sqrt < 0 in sigma'
        return
     else
        sigma = sqrt(sigma)
     endif
     print*,'period = ',2.*pi/sigma
!
!--solve for alpha(t)
!    
     const = 4.*sigma**2 + 4.*A0**2 
     term = 1.-4.*sigma**2/const
     if (term.le.0.) then
        if (abs(A0).gt.1.e-3) print*,'warning: const or omega wrong, sqrt < 0 : assuming static solution'
        A = 0.
     else
        term = sqrt(term)
        !
        !--this is the solution to the 2nd order ODE for alpha
        !
        A = omega*COS(2.*sigma*time)*term/(1. + SIN(2.*sigma*time)*term)
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

end module toystar2D
