!
! modules containing routines for exact toystar solutions in 2D
! contains two modules :
!
! toystar2D_utils (containing functions utilised in the exact solution)
! exact_toystar2D (containing routine for plotting the exact solution)
!
module toystar2D_utils
 implicit none

contains
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
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma,denom
  real :: ak,akprev,gamm1,freqsq
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

  etar = 0.
  akprev = 1.0  ! this is a_0 which is the amplitude
  !!print*,'mode = ',j,' nu^2 = ',freqsq,' a_0 = ',akprev
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
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma,denom,term1,term2
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
  term1 = 0.
  term2 = 0.
  akprev = 1.0  ! this is a_0 which is the amplitude
!  print*,'mode = ',j,' nu^2 = ',freqsq,' a_0 = ',akprev
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

end module toystar2D_utils
