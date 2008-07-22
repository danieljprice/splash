! ----------------------------------------------------------------------
! compute exact solution for the one dimensional Riemann problem
! in Special Relativity
!
! input parameters are initial left and right states of 
! density, pressure and velocity
!
! Computes shock profile at time t
!
! Calls the exact solution routine provided by Marti & Mueller
!
! Daniel Price, Institute of Astronomy, Cambridge, 2004
!               University of Exeter 2004-2006
!
! dprice@astro.ex.ac.uk
!-----------------------------------------------------------------------
module shock
 implicit none
 public :: exact_shock
 
contains

subroutine exact_shock(iplot,time,gamma,rho_L,rho_R,p_L,p_R,v_L,v_R,xplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot
  integer, intent(out) :: ierr
  real, intent(in) :: time,gamma
  real, intent(in) :: rho_L,rho_R,p_L,p_R,v_L,v_R
  real, dimension(:), intent(inout) :: xplot
  real, dimension(size(xplot)), intent(out) :: yplot
  real*8, dimension(size(xplot)) :: rad,dens,pr,vel,uu
  real*8 :: rhol,rhor,pl,prr,vl,vr,gam,t
  
  print*,'Plotting Special Relativistic Riemann solution at t = ',time,' gamma = ',gamma
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

  rhol = rho_L
  rhor = rho_R
  pl = p_L
  prr = p_R
  vl = v_L
  vr = v_R
  gam = gamma
  t = time
  call riemann(size(xplot),rad,dens,pr,vel,uu,rhol,rhor,pl,prr,vl,vr,gam,t,0.0d0)
  xplot = real(rad)

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
  real, parameter :: tol = 1.5e-4
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

  print*,'its = ',its,' pr = ',prnew,'v = ',vstar,v_R + f_R

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
