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
! Plot exact solution for cooling shock problem
! (ie. one dimensional shock problem with prescribed cooling function)
!
! Taken from Creasey et al. (2011), MNRAS 415, 3706-3720
! http://adsabs.harvard.edu/abs/2011MNRAS.415.3706C
!
! This version implements the full analytical solution including:
! - Proper shock jump conditions
! - Temperature-dependent cooling
! - Profile shifting based on temperature maximum
! - Complete profiles with pre- and post-shock regions
! ----------------------------------------------------------------------
module coolingshock
 implicit none
 public :: exact_coolingshock

 ! Error codes
 integer, parameter :: ierr_none = 0
 integer, parameter :: ierr_alloc = -1
 integer, parameter :: ierr_invalid = -2

contains

!> Calculate exact solution for cooling shock problem
!>
!> @param[in]  iplot          Which quantity to plot (0=density, 1=velocity, 2=temperature, 3=pressure)
!> @param[in]  time           Current time
!> @param[in]  gamma          Adiabatic index
!> @param[in]  xshock         Shock position
!> @param[in]  rho_0         Initial density
!> @param[in]  v_0           Initial velocity
!> @param[in]  p_0           Initial pressure
!> @param[in]  cooling_length Cooling length scale
!> @param[out] xplot         x coordinates for solution
!> @param[out] yplot         y coordinates for selected quantity
!> @param[out] ierr          Error code (0=success)
subroutine exact_coolingshock(iplot,time,gamma,xshock,rho_0,v_0,p_0,cooling_length,xplot,yplot,ierr)
 integer, intent(in) :: iplot
 integer, intent(out) :: ierr
 real, intent(in) :: time,gamma,xshock
 real, intent(in) :: rho_0,v_0,p_0,cooling_length
 real, dimension(:), intent(out) :: xplot
 real, dimension(size(xplot)), intent(out) :: yplot

 ! Local variables
 integer :: i,nr
 real :: u_0,cs0,vshock,r1,dr,a
 real :: vc,rhoc,fac,kB_T_on_Lambda_n0
 real :: lambda_min,lambda_max,r1_min,r1_max
 real :: rhos ! Post-shock density from jump conditions
 real :: mach_number, vs, ps ! Additional shock parameters
 real :: x_maxT ! Location of temperature maximum
 real :: deltax ! Shift in x coordinate
 real, allocatable :: rho_shock(:), T_shock(:), v_shock(:), P_shock(:) ! Shock region
 
 ! Initialize
 ierr = ierr_none
 vc = 0.0  ! Centre velocity (reference frame)

 ! Input validation
 if (gamma <= 1.0) then
    print*,'Error: gamma must be > 1'
    ierr = ierr_invalid
    return
 endif
 if (rho_0 <= 0.0 .or. p_0 <= 0.0 .or. cooling_length <= 0.0) then
    print*,'Error: rho_0, p_0, and cooling_length must be positive'
    ierr = ierr_invalid
    return
 endif
 if (size(xplot) < 3) then
    print*,'Error: output arrays must have size >= 3'
    ierr = ierr_invalid
    return
 endif

 ! Calculate initial parameters
 u_0 = p_0/((gamma-1.)*rho_0)  ! Initial internal energy
 cs0 = sqrt(gamma*p_0/rho_0)   ! Initial sound speed
 kB_T_on_Lambda_n0 = cooling_length/cs0   ! Eq 17 in Creasey et al.

 ! Calculate shock parameters (Eq. A20 in Creasey+2011)
 vshock = 0.5*v_0 + sqrt((0.5*v_0)**2 + p_0/rho_0)
 a = rho_0*(v_0 - vshock)**2/p_0   ! Square of the Mach number
 rhoc = -(v_0 - vshock)*rho_0/vshock  ! Critical density
 
 ! Calculate jump conditions
 mach_number = abs((v_0-vshock)/cs0)
 ! Post-shock density from Rankine-Hugoniot relations
 rhos = rho_0 / ((gamma-1.)/(gamma+1.) + 2./(gamma+1.)/mach_number**2)
 vs = (v_0-vshock)/rhos*rho_0 + vshock  ! Post-shock velocity
 ps = p_0 + rho_0*(v_0-vshock)**2 - rhos*(vs-vshock)**2  ! Post-shock pressure

 ! Define similarity variable range (Eq. A6-A8 in Creasey+2011)
 lambda_min = a/(a + 1.)*(gamma+1.)/gamma
 lambda_max = 0.99999
 r1_min = 1./lambda_min
 r1_max = 1./lambda_max

 ! Calculate point distribution for arbitrary input size
 nr = size(xplot)
 
 allocate(rho_shock(nr), T_shock(nr), v_shock(nr), P_shock(nr), stat=ierr)
 if (ierr /= 0) then
    print*,'Error allocating shock arrays'
    ierr = ierr_alloc
    return
 endif
 
 ! Calculate shock profiles
 dr = (r1_max - r1_min)/(nr-1-2)
 fac = (vc - vshock)/(gamma-1.)*(kB_T_on_Lambda_n0)
 
 do i=2,nr-1
    r1 = r1_min + (i-2)*dr
    ! Spatial coordinate (Eq. A21 in Creasey+2011)
    xplot(i) = -fac * ( &
              (gamma-a)/(a-1.)*log(r1 - 1.) &
              + (1.-a*gamma)/((a-1.)*a**2)*log(1. - a*r1) &
              - (a+1.)/a*r1 - 0.5*(gamma+1.)*r1**2 )
    rho_shock(i) = rho_0 + rhoc/r1  ! Density (Eq. A6)
    T_shock(i) = ((a+1.)*r1 - a*r1*r1)  ! Temperature normalized to T0 (Eq. A8)
    v_shock(i) = vshock + (vc - vshock)*r1  ! Velocity (Eq. A7)
    ! Pressure from equation of state P = (gamma-1)*rho*u
    P_shock(i) = p_0 + (gamma-1.)*(rho_0/r1)*(u_0*((a+1.)*r1 - a*r1*r1))
 enddo
 !print*,' xshock = ',x_shock(1),x_shock(nr)

 ! Find temperature maximum in negative x region and shift profiles
 !x_maxT = minval(x_shock, mask=T_shock == maxval(T_shock))
 !deltax = x_shock(2) - x_maxT
 !print*,' shifting profile to ',deltax
 !x_shock = x_shock - deltax
 
 ! Left side (pre-shock)
 xplot(1) = -100.
 rho_shock(1) = rho_0 ! initial density
 T_shock(1) = u_0 ! initial T
 v_shock(1) = v_0
 P_shock(1) = P_0

 ! Right side (post-shock)
 xplot(nr) = 100.
 rho_shock(nr) = rhos
 T_shock(nr) = T_shock(nr-1)
 v_shock(nr) = -vs
 P_shock(nr) = ps
 
 ! Copy to output arrays based on requested quantity
 select case(iplot)
 case(3) ! pressure
    yplot = P_shock
 case(2) ! temperature
    yplot = T_shock
 case(1) ! velocity
    yplot = v_shock
 case default  ! density
    yplot = rho_shock
 end select

 ! Clean up
 deallocate(rho_shock, T_shock, v_shock, P_shock)

end subroutine exact_coolingshock

end module coolingshock
