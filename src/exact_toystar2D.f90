!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!------------------------------------------------------------
! Exact solutions for toystar in two dimensions
!
! For details see Monaghan and Price (2005), in prep.
!
!------------------------------------------------------------
module toystar2D
 implicit none
 public :: exact_toystar2D
 public :: etar, detadr  ! public because they are used in setup

contains

!------------------------------------------------------------
! calculate exact solution for toystar in two dimensions
!
! non-linear solution solves ODEs, assumes linear velocity
!
! the solutions are all plots against radius
!
! iplot = 0 gives x vs y
!
! iplot = 1->5 gives rho, pr, u, vx, vy vs r
!------------------------------------------------------------

subroutine exact_toystar2D(iplot,time,gamma,polyk,totmass, &
                           ampl,denscentre,C0,jorder,morder, &
                           V11,V22,V12,V21,xplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot,jorder,morder
  integer, intent(out) :: ierr
  real, intent(in) :: time,gamma,polyk,totmass
  real, intent(in) :: C0, ampl, denscentre     ! parameters for toy star
  real, intent(in) :: V11,V22,V12,V21
  real, dimension(:), intent(inout) :: xplot
  real, dimension(:), intent(out) :: yplot

  real, parameter :: pi = 3.141592653589
  integer :: i,npts,nsteps
  integer :: jmode,smode
  real :: H,C,D,B,omega,omegasq
  real :: radstar,dx,nu2 !!,scalefac
  real :: rhoplot,deltarho,vplot,deltav
  real :: gamm1,gam1,constK,sigma,period
  real, dimension(8) :: params, paramsp,dparams
  real :: fac,dt,t,phi,dphi,cosphi,sinphi,denom
  real :: massbefore,massafter
  logical linear

  ierr = 1
  npts = size(xplot)

  linear = (jorder.ge.0 .and. morder.ge.0)
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'Error: no toy star solution for isothermal yet'
     ierr = 1
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
  B = 0.
  D = 0.

  if (linear) then
!---------------------------------------------------------------------------
!  linear solution

     print*,' Plotting 2D toy star: linear solution r mode = ',jorder,' phi mode = ',morder
     jmode = jorder   ! radial mode
     smode = morder        ! non-axisymmetric modes (theta)

     ! sigma is the frequency of oscillation
     nu2 = (jmode + smode)*(jmode+smode + 2./gamm1) - smode**2
     if (nu2.le.0.) then
        print*,'Error: nu^2 < 0 in linear toy star  ',nu2
        print*,' radial mode = ',jmode,' theta mode = ',smode
        ierr = 2
        return
     else
        sigma = sqrt(0.5*omegasq*gamm1*nu2)
     endif
     print*,' Amplitude = ',ampl,' period = ',2*pi/sigma,' H,C = ',denscentre,C0

     !!scalefac = polyk*gamma/(sigma*gamm1)
     radstar = sqrt((2.*polyk*gamma*denscentre**gamm1)/gamm1)

     xplot(1) = 0.
     dx = (radstar-xplot(1))/float(npts-1)

     do i=1,npts
        xplot(i) = xplot(1)+dx*(i-1)
        !         print*,i,' x,y = ',xplot(i),yplot(i)
        rhoplot = (denscentre - C0*xplot(i)**2)
        if (rhoplot.le.0.) rhoplot = 0.
        deltarho = etar(jmode,smode,xplot(i)/radstar,gamma)  ! functional form of rho(r)
        !!print*,'deltarho = ',rhoplot,deltarho,xplot(i)
        rhoplot = (rhoplot + deltarho*ampl*SIN(sigma*time))**gam1

        deltav = ampl*detadr(jmode,smode,xplot(i)/radstar,gamma)
        vplot = deltav*COS(sigma*time)

        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for v_r
           yplot(i) = vplot
        case(5)
           yplot(i) = vplot**2
        end select

     enddo

!---------------------------------------------------------------------------
!  non-linear solution
!
  else

     if (jorder.lt.0 .and. smode.lt.0) then
        smode = 2
        jmode = 0
     else
        smode = 0
        jmode = 2
     endif
     print*,'Plotting 2D toy star: non-linear'
     !  solve for H, C and A given initial conditions on v, rho and the time.
     !
     H = denscentre
!
!--this is the static solution, determined from the total mass, polyk, gamma and omega
!
     radstar = sqrt(gamma*totmass/(pi*gamm1))
     H = omegasq*gamm1*radstar**2/(2.*polyk*gamma)
     C = 0.5*gamm1*omegasq/(gamma*polyk)
     print*,' r_star = ',radstar,' rho = (',H,'-',C,'^2)**',gamm1
     D = C
     B = 0.
!
!--now solve the ODEs for V11,V22,V12,V21,H,C,D & B
!  (we use a simple modified euler method)
!
     params(1)= V11
     params(2)= V22
     params(3)= V12
     params(4)= V21
     params(5)= H
     params(6)= C
     params(7)= D
     params(8)= B
     fac= 2.*gamma*polyk*gam1

     massbefore = pi*gamm1/gamma*H**(gamma*gam1)/(sqrt(C*D - B**2))
!
!--get frequency to determine timestep
!
     nu2 = (jmode + smode)*(jmode+smode + 2./gamm1) - smode**2
     if (nu2.le.0.) then
        print*,'Error: nu^2 < 0 in exact toy star  ',nu2
        print*,' radial mode = ',jmode,' theta mode = ',smode
        ierr = 2
        return
     else
        sigma = sqrt(0.5*omegasq*gamm1*nu2)
     endif
!
!--solve ODE's
!
     period = 2.*pi/sigma
     dt = 0.001*period
     nsteps = int(time/dt)
     dt = time/real(nsteps)
     t = 0.
     do i=1,nsteps
        t = t + dt
        call param_derivs(params,dparams,fac,gamm1,omegasq)
        paramsp(:)= params(:) + 0.5*dt*dparams(:)
        call param_derivs(paramsp,dparams,fac,gamm1,omegasq)
        params(:)= params(:) + dt*dparams(:)
     enddo
!
!--have now got solution at current time
!
     H= params(5)
     C= params(6)
     D= params(7)
     B= params(8)

     print*,' solved ODEs to time = ',t,' in ',nsteps,' nsteps '
     massafter = pi*gamm1/gamma*H**(gamma*gam1)/(sqrt(C*D - B**2))
     print*,' conserved mass before = ',massbefore,' after =',massafter

     if (C.le.0.) then
        radstar = 0.5
        stop '*** C = 0 = illegal'
     !!elseif (A.le.1.e-5) then
     else
        radstar = sqrt(H/C)
     endif
     xplot(1) = 0.
     dx = (radstar-xplot(1))/float(npts-1)

     do i=1,npts
        xplot(i) = xplot(1)+dx*(i-1)
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
        case(4)                 ! plot solution for v_r
           yplot(i) = ampl*xplot(i)
        case(5)
           yplot(i) = (ampl*xplot(i))**2
        end select
     enddo
!
!------------------------------------------------------------------------
!
  endif

  if (iplot.gt.0 .and. iplot.le.5) then
     ierr = 0
  elseif (iplot.eq.0) then
     print*,' plotting non-axisymmetric boundary'
!
!--for x-y plots we plot the rho=0 curve (ie. toy star boundary)
!
     dphi = 2.*pi/real(npts-1)
     phi = 0.
     do i=1,npts
        phi = (i-1)*dphi
        cosphi = cos(phi)
        sinphi = sin(phi)
        denom = C*cosphi**2 + 2.*B*cosphi*sinphi + D*sinphi**2
        radstar = SQRT(H/denom)
        xplot(i) = radstar*cosphi
        yplot(i) = radstar*sinphi
     enddo
     ierr = 0
!     call pgsfs(2)
!     call pgcirc(0.0,0.0,radstar)
!     ierr = 3
  endif

  return
end subroutine exact_toystar2D

!
!--function that evaluates the polynomial for rho(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
!  rad = r/r_star
!  j = radial (axisymmetric) mode
!  m = theta mode
!
!  solution is for delta(rho**(gamma-1))
!  ie. rho**(gamma-1) = rho_0**(gamma-1) + etar
!
!  and takes the form
!
!  etar = rad**m sum_k a_k rad**k
!
real function etar(j,m,rad,gamma)
  implicit none
  integer, intent(in) :: j,m ! j is the radial mode, m is the theta mode
  integer :: k,kprev
  real, intent(in) :: rad,gamma
  real :: denom,ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     etar = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  akprev = 1.0  ! this is a_0 which is the amplitude
  etar = akprev
  !!print*,'mode = ',j,m,' nu^2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     denom = real((kprev + 2 + m)**2 - m**2)
     ak = akprev*(kprev**2 + 2.*kprev*m + 2.*(kprev+m)/gamm1 - freqsq)/denom
     !!print*,'coeff ',k,' = ',ak,k**2,2.*k/gamm1
     etar = etar + ak*rad**k
     akprev = ak
  enddo

  etar = etar * rad**m

end function etar

!
!--function that evaluates the polynomial for v(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
real function detadr(j,m,rad,gamma)
  implicit none
  integer, intent(in) :: j, m  ! j is the radial mode, m is the theta mode
  integer :: k,kprev
  real, intent(in) :: rad,gamma
  real :: denom,term1,term2
  real :: ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     detadr = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  detadr = 0.
  akprev = 1.0  ! this is a_0 which is the amplitude
  term1 = akprev
  term2 = 0.
!  print*,'mode = ',j,m,' nu^2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     denom = real((kprev + 2 + m)**2 - m**2)
     ak = akprev*(kprev**2 + 2.*kprev*m + 2.*(kprev+m)/gamm1 - freqsq)/denom
     !!print*,'coeff ',k,' = ',ak,k*ak,rad,(k-1)
     term1 = term1 + ak*rad**k
     term2 = term2 + k*ak*rad**(k-1)
     akprev = ak
  enddo

  if (m.eq.0) then
     detadr = term2
  else
     detadr = m*rad**(m-1)*term1 + rad**m*term2
  endif

end function detadr

subroutine param_derivs(func,dfunc,fac,gamm1,omegasq)
  implicit none
  real, intent(in), dimension(8) :: func
  real, intent(out), dimension(8) :: dfunc
  real, intent(in) :: fac, gamm1,omegasq
  real :: term, gamma

  term = func(1) + func(2)
  gamma = gamm1 + 1.
  dfunc(1)= fac*func(6) - func(1)*func(1) - func(3)*func(4) - omegasq
  dfunc(2)= fac*func(7) - func(2)*func(2) - func(3)*func(4) - omegasq
  dfunc(3)= fac*func(8) - func(3)*term
  dfunc(4)= fac*func(8) - func(4)*term
  dfunc(5)= -gamm1*term*func(5)
  dfunc(6)= -2.*func(6)*func(1) - gamm1*func(6)*term - 2.*func(8)*func(4)
  dfunc(7)= -2.*func(7)*func(2) - gamm1*func(7)*term - 2.*func(8)*func(3)
  dfunc(8)= -func(6)*func(3) - func(7)*func(4) - gamma*func(8)*term

  return
end subroutine param_derivs

!----------------------------------------------------------------------
!
! this subroutine plots the alpha-beta relation in the 2D Toy star solution
!
!----------------------------------------------------------------------
subroutine exact_toystar_ACplane2D(astart,bstart,sigmain,gamma)
  use plotlib, only:plot_swin,plot_box,plot_label,plot_line
  implicit none
  real, intent(in) :: astart,bstart,sigmain,gamma
  integer, parameter :: npts = 2000
  integer :: i
  real :: gamm1,gam1,sigma
  real :: polyk,Omega2,constk
  real :: xstart,xend,ymin,ymax,xi,term,extra
  real, dimension(npts) :: xplot, yplot

  print*,' plotting alpha-beta plane...'

  gamm1 = gamma - 1.
  gam1 = 1./gamm1
  polyk = 0.25
  Omega2 = 1.0
  sigma = 1.0
  print*,' alpha = ',astart,' beta = ',bstart, 'sigma = ',sigma,sigmain
  !
  !--find integration constant from starting values of alpha and beta
  !
  constk = (astart**2 + bstart**2 + Omega2 &
            + 2.*polyk*gamma*(sigma*bstart**gamma)*gam1**2)/bstart

  print*,' integration constant = ',constk

  !
  !--find limits of plot (ie. where alpha = 0)
  !
  xstart = 0.25
  xend = 2.0

  !      print*,'plotting k = ',k,' cstart = ',cstart,' astart = ',astart
  !      print*,'min c = ',xstart,' max c = ',xend

  xstart = xstart + 0.000001
  xend = xend - 0.000001

  extra = 0.1*(xend-xstart)
  !!xcentre = 0.5*(xstart + xend)
  ymax = 2.0
  ymin = -2.0

  call plot_swin(xstart-extra,xend+extra,ymin,ymax)
  call plot_box('bcnst',0.0,0,'1bvcnst',0.0,0)
  call plot_label ('beta','alpha',' ')

  do i=1,npts
     xi = xstart + (i-1)*npts
     xplot(i) = xi
     term = -(xi**2 + Omega2 + 2.*polyk*gamma*(sigma*xi**gamma)*gam1**2 + constk*xi)
     if (term.le.0) then
        yplot(i) = 0.
     else
        yplot(i) = sqrt(term)
     endif
  enddo
  call plot_line(npts,xplot,yplot)

  return

end subroutine exact_toystar_ACplane2D
end module toystar2D
