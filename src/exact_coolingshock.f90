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
!  Copyright (C) 2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
! Plot exact solution for the cooling shock problem: two identical
! streams (rho_0, p_0, +/- v_0) collide at x0 and form a cold dense
! layer bounded by two radiative shocks, with a linear cooling function
! in the post-shock gas.
!
! Taken from Creasey et al. (2011), MNRAS 415, 3706-3720
! http://adsabs.harvard.edu/abs/2011MNRAS.415.3706C
!
! The post-shock cooling region is the stationary similarity solution
! of Appendix A1, in terms of lambda = rho/rho_c. The shock speed and
! cold state follow from the isothermal jump conditions (A17-A20).
! The shock position follows from mass conservation: the cold layer
! plus the cooling region must contain all the gas swept up, giving
!
!   |x_s - x0| = |u_s| t + D/(rho_c - rho_0),
!
! where D = int (rho_c - rho) dx is the mass deficit of the cooling
! region, computed here in closed form. This is exact once the cooling
! region is stationary, i.e. for t >> cooling_length/c_0.
! ----------------------------------------------------------------------
module coolingshock
 implicit none
 public :: exact_coolingshock

 private

contains

!-----------------------------------------------------------------------
!+
!  Exact solution for the cooling shock problem
!
!  iplot  : 0=density, 1=velocity, 2=thermal energy, 3=pressure
!  xshock : position where the two streams collide
!  v_0    : speed of the left stream (> 0), right stream has -v_0
!  cooling_length : kB T0 c0/(Lambda n0), eq. 17 of Creasey et al.
!+
!-----------------------------------------------------------------------
subroutine exact_coolingshock(iplot,time,gamma,xshock,rho_0,v_0,p_0,cooling_length,xplot,yplot,ierr)
 integer, intent(in)  :: iplot
 integer, intent(out) :: ierr
 real,    intent(in)  :: time,gamma,xshock
 real,    intent(in)  :: rho_0,v_0,p_0,cooling_length
 real, dimension(:), intent(out) :: xplot
 real, dimension(size(xplot)), intent(out) :: yplot
 real, parameter :: xfar = 1.e4  ! distance used for the far upstream point
 integer :: i,nr,nh
 real :: c0,u_0,us,rhoc,a,mach,rhos,lambda_s,prefac,deficit,dshock
 real :: us_ad,rhoc_ad,pc_ad

 ierr = 0
 nr = size(xplot)
 if (gamma <= 1.) then
    print "(a)",' ERROR: gamma must be > 1 for cooling shock solution'
    ierr = 1
    return
 endif
 if (rho_0 <= 0. .or. p_0 <= 0. .or. v_0 <= 0. .or. cooling_length <= 0.) then
    print "(a)",' ERROR: density, pressure, velocity and cooling length must be > 0 for cooling shock'
    ierr = 2
    return
 endif
 if (nr < 10) then
    print "(a)",' ERROR: too few points for cooling shock solution'
    ierr = 3
    return
 endif

 c0  = sqrt(gamma*p_0/rho_0)        ! upstream sound speed
 u_0 = p_0/((gamma-1.)*rho_0)       ! upstream (= cold layer) thermal energy

 ! velocity of the left-moving shock (A20) and density of the cold central gas (A17)
 us   = 0.5*v_0 - sqrt(0.25*v_0**2 + p_0/rho_0)
 rhoc = (us - v_0)*rho_0/us
 a    = rho_0*us**2/p_0             ! eq. A9 with the cold gas (v=0) as reference state

 ! density immediately behind the shock from the Rankine-Hugoniot relations (A21-A23)
 mach = abs(v_0 - us)/c0
 if (mach <= 1.) then
    print "(a)",' ERROR: injection velocity too low to form a shock in cooling shock solution'
    ierr = 4
    return
 endif
 rhos = rho_0/((gamma-1.)/(gamma+1.) + 2./((gamma+1.)*mach**2))
 lambda_s = rhos/rhoc
 ! the logarithms in F need a < lambda < 1 throughout the cooling region
 if (lambda_s <= a .or. lambda_s >= 1.) then
    print "(a)",' ERROR: no valid cooling region for these parameters in cooling shock solution'
    ierr = 5
    return
 endif

 ! distance behind the shock is xi = prefac*(F(lambda) - F(lambda_s)), from eq. A11
 prefac = -us*cooling_length*rho_0/((gamma-1.)*c0*rhoc)

 ! mass deficit of the cooling region, int (rho_c - rho) dxi from lambda_s to 1,
 ! obtained by integrating by parts and using the antiderivative G of F
 deficit = rhoc*prefac*(-(1. - lambda_s)*creasey_F(lambda_s,a,gamma) &
                        + creasey_G1(a,gamma) - creasey_G(lambda_s,a,gamma))

 ! distance of each shock from the collision point once the cooling region is stationary
 dshock = abs(us)*time + deficit/(rhoc - rho_0)

 ! shock speed, density and pressure of the hot central gas without cooling (A12-A16),
 ! which is the exact solution at early times before the cooling region forms
 us_ad   = -0.25*(gamma - 3.)*v_0 - 0.25*sqrt((gamma + 1.)**2*v_0**2 + 16.*p_0/rho_0)
 rhoc_ad = (us_ad - v_0)*rho_0/us_ad
 pc_ad   = rhoc_ad*(p_0/rho_0 + 0.5*(gamma - 1.)*v_0**2)

 ! fill the left half (x <= xshock) and mirror it about xshock
 nh = nr/2
 xplot(1) = xshock - xfar
 yplot(1) = select_quantity(iplot,rho_0,v_0,u_0,gamma)
 if (abs(us_ad)*time < dshock) then
    ! adiabatic shock is behind the stationary cooling shock: cooling not yet established
    call fill_plateau(iplot,gamma,xshock-abs(us_ad)*time,xshock,rho_0,v_0,u_0, &
                      rhoc_ad,pc_ad/((gamma-1.)*rhoc_ad),xplot(2:nh),yplot(2:nh))
 else
    call fill_cooling(iplot,gamma,xshock-dshock,xshock,rho_0,v_0,u_0,us,rhoc,a,lambda_s,prefac, &
                      xplot(2:nh),yplot(2:nh))
 endif

 ! right half is the mirror image, with velocity reversed
 do i=1,nh
    xplot(nr+1-i) = 2.*xshock - xplot(i)
    yplot(nr+1-i) = yplot(i)
    if (iplot==1) yplot(nr+1-i) = -yplot(i)
 enddo
 if (mod(nr,2)==1) then
    xplot(nh+1) = xshock
    yplot(nh+1) = yplot(nh)
 endif

end subroutine exact_coolingshock

!-----------------------------------------------------------------------
!+
!  Left half of the profile with a uniform state behind the shock at xs,
!  extending to the collision point x0 (with xs = x0 this is the initial jump)
!+
!-----------------------------------------------------------------------
subroutine fill_plateau(iplot,gamma,xs,x0,rho_0,v_0,u_0,rhoc,uc,xplot,yplot)
 integer, intent(in)  :: iplot
 real,    intent(in)  :: gamma,xs,x0,rho_0,v_0,u_0,rhoc,uc
 real,    intent(out) :: xplot(:),yplot(:)
 integer :: i,n

 n = size(xplot)
 xplot(1) = xs
 yplot(1) = select_quantity(iplot,rho_0,v_0,u_0,gamma)
 do i=2,n
    xplot(i) = xs + (x0 - xs)*real(i-2)/real(n-2)
    if (xs < x0) then
       yplot(i) = select_quantity(iplot,rhoc,0.,uc,gamma)
    else
       yplot(i) = yplot(1)
    endif
 enddo

end subroutine fill_plateau

!-----------------------------------------------------------------------
!+
!  Left half of the profile with the stationary cooling region behind
!  the shock at xs, followed by cold gas at rest up to the collision
!  point x0. Points are spaced evenly in log(1-lambda), since
!  x ~ -log(1-lambda) near lambda = 1
!+
!-----------------------------------------------------------------------
subroutine fill_cooling(iplot,gamma,xs,x0,rho_0,v_0,u_0,us,rhoc,a,lambda_s,prefac,xplot,yplot)
 integer, intent(in)  :: iplot
 real,    intent(in)  :: gamma,xs,x0,rho_0,v_0,u_0,us,rhoc,a,lambda_s,prefac
 real,    intent(out) :: xplot(:),yplot(:)
 real, parameter :: lambda_max = 1. - 1.e-4  ! end of the cooling region (rho = 0.9999 rho_c)
 integer :: i,n
 real :: lambda,rho,v,u

 n = size(xplot)
 xplot(1) = xs
 yplot(1) = select_quantity(iplot,rho_0,v_0,u_0,gamma)
 do i=2,n-1
    lambda = 1. - (1. - lambda_s)*((1. - lambda_max)/(1. - lambda_s))**(real(i-2)/real(n-3))
    rho = rhoc*lambda                             ! eq. A6
    v   = us*(1. - 1./lambda)                     ! eq. A7
    u   = u_0*((a + 1.)/lambda - a/lambda**2)     ! eq. A8
    xplot(i) = xs + prefac*(creasey_F(lambda,a,gamma) - creasey_F(lambda_s,a,gamma))
    yplot(i) = select_quantity(iplot,rho,v,u,gamma)
    ! at early times the cooling region is truncated at the collision point
    if (xplot(i) >= x0) then
       xplot(i:n) = x0
       yplot(i:n) = yplot(i)
       return
    endif
 enddo
 xplot(n) = x0
 yplot(n) = select_quantity(iplot,rhoc,0.,u_0,gamma)

end subroutine fill_cooling

!-----------------------------------------------------------------------
!+
!  Return the requested quantity given rho, v and u
!+
!-----------------------------------------------------------------------
real function select_quantity(iplot,rho,v,u,gamma)
 integer, intent(in) :: iplot
 real,    intent(in) :: rho,v,u,gamma

 select case(iplot)
 case(3)
    select_quantity = (gamma-1.)*rho*u
 case(2)
    select_quantity = u
 case(1)
    select_quantity = v
 case default
    select_quantity = rho
 end select

end function select_quantity

!-----------------------------------------------------------------------
!+
!  Bracket in x(lambda) from integrating eq. A11:
!  F = A ln(1/lambda-1) + B ln(1-a/lambda) - C/lambda - E/lambda^2
!+
!-----------------------------------------------------------------------
real function creasey_F(lambda,a,gamma)
 real, intent(in) :: lambda,a,gamma
 real :: ca,cb,cc,ce

 call creasey_coeffs(a,gamma,ca,cb,cc,ce)
 creasey_F = ca*log(1./lambda - 1.) + cb*log(1. - a/lambda) - cc/lambda - ce/lambda**2

end function creasey_F

!-----------------------------------------------------------------------
!+
!  Antiderivative of F with respect to lambda (a < lambda < 1)
!+
!-----------------------------------------------------------------------
real function creasey_G(lambda,a,gamma)
 real, intent(in) :: lambda,a,gamma
 real :: ca,cb,cc,ce

 call creasey_coeffs(a,gamma,ca,cb,cc,ce)
 creasey_G = ca*(-(1. - lambda)*log(1. - lambda) - lambda*log(lambda)) &
           + cb*((lambda - a)*log(lambda - a) - lambda*log(lambda)) &
           - cc*log(lambda) + ce/lambda

end function creasey_G

!-----------------------------------------------------------------------
!+
!  Limit of the antiderivative of F as lambda -> 1
!+
!-----------------------------------------------------------------------
real function creasey_G1(a,gamma)
 real, intent(in) :: a,gamma
 real :: ca,cb,cc,ce

 call creasey_coeffs(a,gamma,ca,cb,cc,ce)
 creasey_G1 = cb*(1. - a)*log(1. - a) + ce

end function creasey_G1

!-----------------------------------------------------------------------
!+
!  Coefficients A, B, C and E of the terms in F
!+
!-----------------------------------------------------------------------
subroutine creasey_coeffs(a,gamma,ca,cb,cc,ce)
 real, intent(in)  :: a,gamma
 real, intent(out) :: ca,cb,cc,ce

 ca = (gamma - a)/(a - 1.)
 cb = (1. - a*gamma)/((a - 1.)*a**2)
 cc = (a + 1.)/a
 ce = 0.5*(gamma + 1.)

end subroutine creasey_coeffs

end module coolingshock
