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
subroutine riemannsolver(gamma,p_L,p_R,v_L,v_R,c_L,c_R,pr,vstar)
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

end subroutine riemannsolver

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
