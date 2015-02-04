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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
! compute exact solution for the one dimensional Riemann problem
! (hydrodynamic shock)
!
! input parameters are initial left and right states of
! density, pressure and velocity
!
! Computes shock profile at time t
!
! Calls a separate subroutine to calculate the post-shock pressure
! and velocity (this is the difficult bit).
!
! Handles all cases of left and right-going shocks and rarefactions
!
! Daniel Price, Institute of Astronomy, Cambridge, 2004
!               University of Exeter 2004-2008
!               Monash University 2008-
!
! dprice@astro.ex.ac.uk
!-----------------------------------------------------------------------
module shock
 implicit none
 public :: exact_shock
 private :: get_pstar, get_pstar_isothermal, f_and_df

contains

subroutine exact_shock(iplot,time,gammain,rho_L,rho_R,p_L,p_R,v_L,v_R,rdust_to_gas,xplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot
  integer, intent(out) :: ierr
  real, intent(in) :: time,gammain
  real, intent(in) :: rho_L,rho_R,p_L,p_R,v_L,v_R,rdust_to_gas
  real, dimension(:), intent(in) :: xplot
  real, dimension(size(xplot)), intent(out) :: yplot

  integer :: i
  real, dimension(size(xplot)) :: dens, pr, vel
  real :: cs_L,cs_R, gamfac
  real :: ppost, vpost, vleft, vright, gamma
  real :: xzero,xleft,xleftleft,xcontact,xright,xrightright
  logical :: useisothermal,leftisshock,rightisshock

  gamma = gammain
  print*,'Plotting exact Riemann solution at t = ',time,' gamma = ',gamma
!
! check for errors in input
!
  ierr = 0
  if (rho_L.le.0. .or. rho_R.le.0.) then
     print*,'error: rho <= 0 on input : ',rho_L,rho_R
     ierr = 1
     return
  elseif (p_L .le.0. .or. p_R .le.0.) then
     print*,'error: pr <= 0 on input ',p_L, p_R
     ierr = 2
     return
  endif
  if (gamma < 1.) then
     print*,'Error: gamma = ',gamma,' setting to 5/3'
     gamma = 5./3.
  endif
!
!  xzero is the position of the shock at t=0
!
  xzero = 0.
!
!  define sound speeds to left and right of shock tube
!
  cs_L = sqrt(gamma*p_L/rho_L)
  cs_R = sqrt(gamma*p_R/rho_R)
  if (rdust_to_gas .gt.epsilon(rdust_to_gas)) then
     cs_L = cs_L*sqrt(1./(1.+rdust_to_gas))
     cs_R = cs_R*sqrt(1./(1.+rdust_to_gas))   
  endif
  gamfac = (gamma-1.)/(gamma + 1.)

!------------------------------------------------------------
! find post-shock pressure via the Riemann solver
! (this version also returns the post shock speed vpost
!  although this can be calculated from ppost)
!------------------------------------------------------------

  if (gamma.gt.1.0001) then
     call get_pstar(gamma,p_L,p_R,v_L,v_R,cs_L,cs_R,ppost,vpost)
     useisothermal = .false.
  else
     print*,'using isothermal solver...',p_L/rho_L, p_R/rho_R
     useisothermal = .true.
     call get_pstar_isothermal(cs_L*cs_L,v_L,v_R,rho_L,rho_R,ppost,vpost)
  endif

!------------------------------------------------------------
!  using this, calculate various speeds needed in order to
!  reconstruct the shock profile
!------------------------------------------------------------

  !
  !  check whether solutions are shocks or rarefactions
  !
  if (ppost .gt. p_L) then
     leftisshock = .true.
  else
     leftisshock = .false.
  endif
  if (ppost .gt. p_R) then
     rightisshock = .true.
  else
     rightisshock = .false.
  endif
  if (leftisshock) then
     write(*,"(a)",advance='no') ' left-going wave is a shock; '
  else
     write(*,"(a)",advance='no') ' left-going wave is a rarefaction; '
  endif
  if (rightisshock) then
     write(*,"(a)") 'right-going wave is a shock'
  else
     write(*,"(a)") 'right-going wave is a rarefaction'
  endif

  if (rightisshock) then ! right hand wave is a shock
  !
  !  speed of the shock front
  !
     vright = v_R + cs_R**2*(ppost/p_R - 1.)/(gamma*(vpost-v_R))
  else ! right hand wave is a rarefaction
  !
  !  speed at which the right end of rarefaction fan moves
  !
     vright = cs_R + 0.5*(gamma+1.)*vpost - 0.5*(gamma-1.)*v_R
  endif
  !
  !  repeat for left-going wave
  !
  if (leftisshock) then
     vleft = -(v_L + cs_L**2*(ppost/p_L - 1.)/(gamma*(vpost-v_L)))
  else
     vleft = cs_L - 0.5*(gamma+1.)*vpost + 0.5*(gamma-1.)*v_L
  endif

!-------------------------------------------------------------
! now work out the locations of various features in the shock
!-------------------------------------------------------------

!
! position of left-going shock or back end of left-going rarefaction
!
  xleft = xzero - vleft*time
  if (leftisshock) then
     xleftleft = xleft
  else
!
! front end of left-going expansion fan (propagates at sound speed into "left" fluid)
!
     xleftleft = xzero - (cs_L - v_L)*time
  endif
!
! position of the interface between the "left" fluid from the "right" fluid
! (the contact discontinuity)
!
  xcontact = xzero + vpost*time

!
! position of right-going shock or back end of right-going rarefaction
!
  xright = xzero + vright*time
  if (rightisshock) then
     xrightright = xright
  else
!
! right end of right-going expansion fan (propagates at sound speed into "right" fluid)
!
     xrightright = xzero + (cs_R + v_R)*time
  endif

!--------------------------------------------------------------
! reconstruct the shock profile for all x
!--------------------------------------------------------------

!--here is a cheap, dirty f90 version for crap compilers
  do i=1,size(xplot)
     if (xplot(i) <= xleftleft) then
!       undisturbed medium to the left
        pr(i) = p_L
        dens(i) = rho_L
        vel(i) = v_L
     elseif (xplot(i) < xleft) then
        if (leftisshock) then
           pr(i) = ppost
           dens(i) = rho_L*(gamfac+ppost/p_L)/(1+gamfac*ppost/p_L)
!           dens(i) = rho_L*(ppost/p_L)**(1./gamma)
           vel(i) = vpost
        else
!       inside expansion fan
           if (useisothermal) then ! this is a bit of a guess
              dens(i) = rho_L*exp((xleftleft-xplot(i))/(cs_L*time) + v_L/cs_L)
           else
              dens(i) = rho_L*(gamfac*(xzero-xplot(i))/(cs_L*time) + gamfac*v_L/cs_L + (1.-gamfac))**(2./(gamma-1.))
           endif
           pr(i) = p_L*(dens(i)/rho_L)**gamma
           vel(i) = (1.-gamfac)*(cs_L -(xzero-xplot(i))/time) + gamfac*v_L
        endif
     elseif (xplot(i) < xcontact) then
!       between left expansion fan/shock and contact discontinuity
!       post-shock, ahead of contact discontinuity but before right going wave
        pr(i) = ppost
        if (leftisshock) then
           dens(i) = rho_L*(gamfac+ppost/p_L)/(1+gamfac*ppost/p_L)
        else
           dens(i) = rho_L*(ppost/p_L)**(1./gamma)
        endif
        vel(i) = vpost
     elseif (xplot(i) < xright) then
!       post-shock, ahead of contact discontinuity but before right going wave
        pr(i) = ppost
        if (rightisshock) then
           dens(i) = rho_R*(gamfac+ppost/p_R)/(1+gamfac*ppost/p_R)
        else
           dens(i) = rho_R*(ppost/p_R)**(1./gamma)
        endif
        vel(i) = vpost
     elseif (xplot(i) < xrightright) then
        if (rightisshock) then
!       irrelevant as in this case xrightright = xright
        else
!        inside expansion fan to right
           if (useisothermal) then ! this is a bit of a guess
              dens(i) = rho_R*exp(-(xrightright-xplot(i))/(cs_R*time) - v_R/cs_R)
           else
              dens(i) = rho_R*(gamfac*(xplot(i)-xzero)/(cs_R*time) - gamfac*v_R/cs_R + (1.-gamfac))**(2./(gamma-1.))
           endif
           pr(i) = p_R*(dens(i)/rho_R)**gamma
           vel(i) = (1.-gamfac)*(-cs_R - (xzero-xplot(i))/time) + gamfac*v_R
        endif
     else
!       undisturbed medium to the right
        pr(i) = p_R
        dens(i) = rho_R
        vel(i) = v_R
     endif
  enddo

!--this is the beautiful, f95 version (which won't compile on pgf90)
!  where(xplot <= xleft)  ! <= otherwise problems at t=0
!!    undisturbed medium to the left
!     pr = p_L
!     dens = rho_L
!     vel = v_L
!  elsewhere(xplot < xfan)
!!    inside expansion fan
!     dens = rho_L*(gamfac*(xzero-xplot)/(cs_L*time) + (1.-gamfac))**(2./(gamma-1.))
!     pr = p_L*(dens/rho_L)**gamma
!     vel = (1.-gamfac)*(cs_L -(xzero-xplot)/time)
!  elsewhere(xplot < xcontact)
!!    between expansion fan and contact discontinuity
!     pr = ppost
!     dens = rho_L*(ppost/p_L)**(1./gamma)
!     vel = vpost
!  elsewhere(xplot < xshock)
!!    post-shock, ahead of contact discontinuity
!     pr = ppost
!     dens = rho_R*(gamfac+ppost/p_R)/(1+gamfac*ppost/p_R)
!     vel = vpost
!  elsewhere
!!    undisturbed medium to the right
!     pr = p_R
!     dens = rho_R
!     vel = v_R
!  end where

!------------------------------------
!  determine which solution to plot
!------------------------------------
  select case(iplot)
  case(1)
     yplot = dens
  case(2)
     yplot = pr
  case(3)
     yplot = vel
  case(4)
     if (gamma.gt.1.0001) then
        yplot = pr/((gamma-1.)*dens)
     else
        yplot = pr/dens
     endif
  case(5) ! deltav, where vd = 0
     yplot = -vel
  case(6) ! eps, where rhod = const
     yplot = rho_R/(dens + rho_R)
  end select

  return
end subroutine exact_shock

!-------------------------------------------------------------------
! Implementation of the exact Riemann solver given in Toro (1992)
!
! Solves for the post-shock pressure (pr) and velocity (vstar)
! given the initial left and right states
!
! Does not matter if high P / high rho is on left or right
!
! Daniel Price, Institute of Astronomy, Cambridge, UK, 2004
! dprice@ast.cam.ac.uk
!-------------------------------------------------------------------
subroutine get_pstar(gamma,p_L,p_R,v_L,v_R,c_L,c_R,pr,vstar)
  implicit none
  real, parameter :: tol = 1.5e-6
  real, intent(in) :: gamma,p_L,p_R,v_L,v_R,c_L,c_R
  real, intent(out) :: pr,vstar
  integer, parameter :: maxits = 30
  integer :: its
  real :: prnew, f_L, f_R, dfdp_L, dfdp_R, f, df, dp
  real :: power, denom
!
!--get an initial starting estimate of intermediate pressure
!  this one is from Toro(1992) - gives basically the right answer
!  for pressure jumps below about 4
!
  power = (gamma-1.)/(2.*gamma)
  denom = c_L/p_L**power + c_R/p_R**power
  prnew = ((c_L + c_R + (v_L - v_R)*0.5*(gamma-1.))/denom)**(1./power)
  pr = p_L
  its = 0

  !!print*,'initial guess = ',prnew

  do while (abs(prnew-pr).gt.tol .and. its.lt.maxits)

     its = its + 1
     pr = prnew
!
!--evaluate the function and its derivatives
!
     call f_and_df(pr,p_L,c_L,gamma,f_L,dfdp_L)
     call f_and_df(pr,p_R,c_R,gamma,f_R,dfdp_R)
!
!--then get new estimate of pr
!
     f = f_L + f_R + (v_R - v_L)
     df = dfdp_L + dfdp_R
!
!--Newton-Raphson iterations
!
     dp =  -f/df
     prnew = pr + dp

  enddo

  if (its.eq.maxits) print*,'WARNING: its not converged in riemann solver'
  pr = prnew
  vstar = v_L - f_L

  print*,'its =',its,' p* =',prnew,'v* =',vstar,v_R + f_R

end subroutine get_pstar

!
!--pressure function
!  H is pstar/p_L or pstar/p_R
!
subroutine f_and_df(prstar,pr,cs,gam,fp,dfdp)
  implicit none
  real, intent(in) :: prstar, pr, gam, cs
  real, intent(out) :: fp,dfdp
  real :: H,term, power, gamm1, denom

  H = prstar/pr
  gamm1 = gam - 1.

  if (H.gt.1.) then  ! shock
     denom = gam*((gam+1.)*H + gamm1)
     term = sqrt(2./denom)
     fp = (H - 1.)*cs*term

     dfdp = cs*term/pr + (H - 1.)*cs/term*(-1./denom**2)*gam*(gam+1.)/pr
  else               ! rarefaction
     power = gamm1/(2.*gam)
     fp = (H**power - 1.)*(2.*cs/gamm1)

     dfdp = 2.*cs/gamm1*power*H**(power-1.)/pr
  endif

end subroutine f_and_df

!-------------------------------------------------------------
! Non-iterative isothermal Riemann solver
! from Balsara (1994), ApJ 420, 197-212
!
! See also Cha & Whitworth (2003), MNRAS 340, 73-90
!-------------------------------------------------------------
subroutine get_pstar_isothermal(cs2,v_L,v_R,rho_L,rho_R,pstar,vstar)
  implicit none
  real, intent(in) :: cs2,v_L,v_R,rho_L,rho_R
  real, intent(out) :: pstar,vstar
  real :: sqrtrho_L, sqrtrho_R, X, vdiff, determinant, vstar2

  sqrtrho_L = sqrt(rho_L)
  sqrtrho_R = sqrt(rho_R)

  X = sqrtrho_L*sqrtrho_R/(sqrtrho_L + sqrtrho_R)
  vdiff = v_L - v_R
  determinant = (X*vdiff)**2 + 4.*cs2*X*(sqrtrho_L + sqrtrho_R)

  pstar = 0.25*(X*vdiff + sqrt(determinant))**2
  vstar = v_L - (pstar - cs2*rho_L)/(sqrt(pstar*rho_L))
  vstar2 = v_R + (pstar - cs2*rho_R)/(sqrt(pstar*rho_R))
  print*,' pstar = ',pstar,' vstar = ',vstar,vstar2

end subroutine get_pstar_isothermal

end module shock
