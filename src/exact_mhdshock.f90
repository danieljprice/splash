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
!  Copyright (C) 2005-2026 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Exact solution of the Riemann problem in ideal magnetohydrodynamics
! (one dimensional MHD shock tubes), for adiabatic or isothermal gas
!
! The solution is found by iterating on the strengths of the fast,
! rotational (Alfven) and slow waves either side of the contact until
! pressure, velocity and transverse field match at the contact, following
!
!  Torrilhon (2002), Research report 2002-06, SAM, ETH Zurich
!  Takahashi & Yamada (2013), J. Plasma Phys. 79, 335
!  Takahashi & Yamada (2014), J. Plasma Phys. 80, 255
!
! The physically relevant solution is not always unique. For problems
! where the transverse fields and velocities lie in a plane (e.g. the
! Brio & Wu 1988 problem) we return the solution of the coplanar
! ("reduced") MHD system, where the slow wave can be an intermediate
! shock or a compound wave (Falle & Komissarov 2001; Takahashi &
! Yamada 2013 sec. 5.1), as found in numerical simulations. Otherwise
! the regular solution is returned. Zero transverse field on one side
! is handled with switch-on shocks and rarefactions.
!
! The preset problems are from Brio & Wu (1988), Ryu & Jones (1995),
! Balsara (1998) and Dai & Woodward (1994)
!
! Daniel Price, Monash University, 2026
! daniel.price@monash.edu
!
!-----------------------------------------------------------------------
module mhdshock
 use cubic,       only:cubicsolve
 use nonlinsolve, only:newton_nd
 implicit none
 public :: exact_mhdshock,get_mhdshock_preset
 integer, parameter, public :: nmhdshocksolns = 7
 character(len=30), dimension(nmhdshocksolns), parameter, public :: mhdprob = &
    (/'Brio & Wu (1988), gamma=2     ', &
      'fast switch-on shock (RJ95)   ', &
      '7 jump shock (RJ95 2a)        ', &
      'isothermal (Balsara 1998)     ', &
      'rarefaction, Bx=0 (RJ95 5a)   ', &
      'Mach 25 (Dai & Woodward 1994) ', &
      'fast/slow shocks (RJ95 1a)    '/)

 private

 integer, parameter :: dp = kind(0.d0)
 integer, parameter :: ifast = 1, islow = 2
 integer, parameter :: imode_general = 1, imode_coplanar = 2, imode_tangential = 3
 integer, parameter :: nstep_raref = 4000  ! integration steps per rarefaction
 integer, parameter :: nrec = 100          ! number of nodes recorded per rarefaction
 integer, parameter :: maxnodes = 4096
 real(dp), parameter :: pi = 4._dp*atan(1._dp)
 real(dp), parameter :: tol_newton = 1.e-11_dp
 real(dp), parameter :: big = huge(1._dp)

 !--settings for the current Riemann problem, shared with the residual function
 real(dp) :: gam,bn,csiso2
 logical  :: isothermal
 integer  :: imode
 real(dp) :: ul0(7),ur0(7),edir(2),smaxl,smaxr,pscale,vscale,bscale

 !--solution nodes (xi = x/t, primitive state)
 logical  :: record = .false.
 integer  :: nnodes
 real(dp) :: xinode(maxnodes),unode(7,maxnodes)

contains

!-----------------------------------------------------------------------
!
! initial left and right states (rho,p,vx,vy,vz,By,Bz) and Bx
! for the preset MHD shock tube problems. Problem 4 requires an
! isothermal equation of state (gamma = 1)
!
!-----------------------------------------------------------------------
subroutine get_mhdshock_preset(ishk,bx,primL,primR,ierr)
 integer, intent(in)  :: ishk
 real,    intent(out) :: bx,primL(7),primR(7)
 integer, intent(out) :: ierr
 real, parameter :: c = 1./sqrt(4.*real(pi))

 ierr = 0
 select case(ishk)
 case(1) ! Brio & Wu (1988)
    bx = 0.75
    primL = (/1.,1.,0.,0.,0.,1.,0./)
    primR = (/0.125,0.1,0.,0.,0.,-1.,0./)
 case(2) ! fast switch-on shock and slow shock (Ryu & Jones 1995)
    bx = 1.
    primL = (/1.,1.,0.,0.,0.,1.,0./)
    primR = (/0.2,0.1,0.,0.,0.,0.,0./)
 case(3) ! seven discontinuities (Ryu & Jones 1995 test 2a)
    bx = 2.*c
    primL = (/1.08,0.95,1.2,0.01,0.5,3.6*c,2.*c/)
    primR = (/1.,1.,0.,0.,0.,4.*c,2.*c/)
 case(4) ! isothermal version of the above (Balsara 1998)
    bx = 2.*c
    primL = (/1.08,1.08,1.2,0.01,0.5,3.6*c,2.*c/)
    primR = (/1.,1.,0.,0.,0.,4.*c,2.*c/)
 case(5) ! rarefactions with no normal field (Ryu & Jones 1995 test 5a)
    bx = 0.
    primL = (/1.,1.,-1.,0.,0.,1.,0./)
    primR = (/1.,1.,1.,0.,0.,1.,0./)
 case(6) ! Mach 25 colliding flows (Dai & Woodward 1994)
    bx = 4.*c
    primL = (/1.,1.,36.87,-0.1546,-0.03864,4.*c,1.*c/)
    primR = (/1.,1.,-36.87,0.,0.,4.*c,1.*c/)
 case(7) ! fast and slow shocks (Ryu & Jones 1995 test 1a)
    bx = 5.*c
    primL = (/1.,20.,10.,0.,0.,5.*c,0./)
    primR = (/1.,1.,-10.,0.,0.,5.*c,0./)
 case default
    print*,'error: unknown MHD shock tube ',ishk
    bx = 0.
    primL = 0.
    primR = 0.
    ierr = 1
 end select

end subroutine get_mhdshock_preset

!-----------------------------------------------------------------------
!
! plot the exact solution of the MHD Riemann problem for arbitrary
! left and right states (rho,p,vx,vy,vz,By,Bz) with normal field Bx.
! gamma = 1 selects an isothermal equation of state, with the
! sound speed given by p/rho in the left state
!
! input  : iplot : variable to plot (1=rho,2=p,3=vx,4=vy,5=vz,6=By,7=Bz,
!                  8=utherm,9=Bx)
! output : xpts,ypts,npts : the solution as a piecewise linear curve
!
!-----------------------------------------------------------------------
subroutine exact_mhdshock(iplot,time,gamma,bx,primL,primR,xmin,xmax,xshock,xpts,ypts,npts,ierr)
 integer, intent(in)  :: iplot
 integer, intent(out) :: npts,ierr
 real,    intent(in)  :: time,gamma,bx,primL(7),primR(7),xmin,xmax,xshock
 real, dimension(:), intent(inout) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 real(dp), allocatable :: uplot(:,:),xplot(:)
 integer  :: i,n

 print "(a,es10.3,a)",' Plotting exact MHD shock tube at t = ',time,' from Takahashi & Yamada (2014)'
 ierr = 0
 npts = 0
 xpts = 0.
 ypts = 0.
 if (primL(1) <= 0. .or. primR(1) <= 0. .or. primL(2) <= 0. .or. primR(2) <= 0.) then
    print*,'error: rho or p <= 0 in initial states'
    ierr = 1
    return
 endif

 allocate(uplot(7,maxnodes+2),xplot(maxnodes+2))
 if (time <= 0.) then
    !--just plot the initial conditions
    n = 4
    xplot(1:4) = (/xmin,xshock,xshock,xmax/)
    uplot(:,1) = primL
    uplot(:,2) = primL
    uplot(:,3) = primR
    uplot(:,4) = primR
 else
    call solve_mhd_riemann(real(gamma,dp),real(bx,dp),real(primL,dp),real(primR,dp),ierr)
    if (ierr /= 0) then
       print*,'WARNING: MHD Riemann solver failed to converge, ierr = ',ierr
       return
    endif
    !--convert xi = x/t to x, adding the undisturbed states at either end
    n = 0
    if (xmin < xshock + xinode(1)*time) call add_point(n,xplot,uplot,real(xmin,dp),unode(:,1))
    do i=1,nnodes
       call add_point(n,xplot,uplot,xshock + xinode(i)*time,unode(:,i))
    enddo
    if (xmax > xplot(n)) call add_point(n,xplot,uplot,real(xmax,dp),unode(:,nnodes))
 endif

 npts = min(n,size(xpts))
 xpts(1:npts) = real(xplot(1:npts))
 select case(iplot)
 case(1:7)
    ypts(1:npts) = real(uplot(iplot,1:npts))
 case(8)
    if (abs(gamma-1.) > 1.e-4) then
       ypts(1:npts) = real(uplot(2,1:npts)/((gamma-1.)*uplot(1,1:npts)))
    else
       ypts(1:npts) = real(uplot(2,1:npts)/uplot(1,1:npts))
    endif
 case(9)
    ypts(1:npts) = bx
 case default
    print*,'error: unknown solution to plot'
    ierr = 2
 end select

end subroutine exact_mhdshock

!-----------------------------------------------------------------------
! append a point (x,u) to the plotted solution
!-----------------------------------------------------------------------
subroutine add_point(n,xplot,uplot,x,u)
 integer,  intent(inout) :: n
 real(dp), intent(inout) :: xplot(:),uplot(:,:)
 real(dp), intent(in)    :: x,u(7)

 n = n + 1
 xplot(n)   = x
 uplot(:,n) = u

end subroutine add_point

!-----------------------------------------------------------------------
!
! solve the MHD Riemann problem, returning the solution in the
! module variables xinode, unode, nnodes
!
!-----------------------------------------------------------------------
subroutine solve_mhd_riemann(gamma,bx,primL,primR,ierr)
 real(dp), intent(in)  :: gamma,bx,primL(7),primR(7)
 integer,  intent(out) :: ierr
 real(dp) :: x(5),urtarget(7)
 integer  :: n,its,istep,nsteps,itry
 integer, parameter :: nsteps_try(2) = (/10,50/)

 gam = gamma
 bn  = bx
 isothermal = (abs(gamma - 1._dp) < 1.e-4_dp)
 csiso2 = primL(2)/primL(1)
 ul0 = primL
 ur0 = primR
 if (isothermal) then
    if (abs(primR(2)/primR(1) - csiso2) > 1.e-6_dp*csiso2) &
       print*,'WARNING: isothermal sound speed differs between left and right states'
    ul0(2) = csiso2*ul0(1)
    ur0(2) = csiso2*ur0(1)
 endif
 call setup_problem(n)

 x = 0._dp
 call newton_nd(residual,x(1:n),n,tol_newton,200,its,ierr)

 !--if this fails, march the right state from the left state to its final value
 urtarget = ur0
 do itry=1,size(nsteps_try)
    if (ierr == 0) exit
    nsteps = nsteps_try(itry)
    x = 0._dp
    do istep=1,nsteps
       ur0 = ul0 + (urtarget - ul0)*real(istep,dp)/real(nsteps,dp)
       call setup_problem(n)
       call newton_nd(residual,x(1:n),n,tol_newton,200,its,ierr)
       if (ierr /= 0) exit
    enddo
 enddo
 if (ierr /= 0) return

 !--evaluate the solution once more, recording the nodes
 record = .true.
 nnodes = 0
 call build_solution(x(1:n),n)
 record = .false.
 if (any(xinode(2:nnodes) < xinode(1:nnodes-1) - 1.e-8_dp*maxval(abs(xinode(1:nnodes))))) &
    print*,'WARNING: waves in MHD Riemann solution are not correctly ordered'

end subroutine solve_mhd_riemann

!-----------------------------------------------------------------------
!
! set the solution mode, scales and default field direction
! for the current left and right states
!
!-----------------------------------------------------------------------
subroutine setup_problem(n)
 integer, intent(out) :: n
 real(dp) :: perp(2),vdiff(2),tol

 pscale = max(ptot(ul0),ptot(ur0))
 bscale = sqrt(pscale)
 vscale = sqrt(pscale/min(ul0(1),ur0(1))) + maxval(abs(ul0(3:5))) + maxval(abs(ur0(3:5)))

 !--default direction for the transverse field
 vdiff = ur0(4:5) - ul0(4:5)
 if (btrans(ul0) > 0._dp) then
    edir = ul0(6:7)/btrans(ul0)
 elseif (btrans(ur0) > 0._dp) then
    edir = ur0(6:7)/btrans(ur0)
 elseif (norm2(vdiff) > 0._dp) then
    edir = vdiff/norm2(vdiff)
 else
    edir = (/1._dp,0._dp/)
 endif

 !--choose between the tangential, coplanar and general solutions
 perp = (/-edir(2),edir(1)/)
 tol  = 1.e-10_dp
 if (abs(bn) <= tiny(bn)) then
    imode = imode_tangential
    n = 2
 elseif (abs(dot_product(ul0(6:7),perp)) <= tol*bscale .and. &
         abs(dot_product(ur0(6:7),perp)) <= tol*bscale .and. &
         abs(dot_product(vdiff,perp)) <= tol*vscale) then
    imode = imode_coplanar
    n = 4
 else
    imode = imode_general
    n = 5
 endif

 smaxl = get_smax(ul0,-1._dp,edir)
 smaxr = get_smax(ur0, 1._dp,edir)

end subroutine setup_problem

!-----------------------------------------------------------------------
!
! residual of the matching conditions at the contact discontinuity,
! as a function of the wave parameters x
!
!-----------------------------------------------------------------------
subroutine residual(x,f,n)
 integer,  intent(in)  :: n
 real(dp), intent(in)  :: x(n)
 real(dp), intent(out) :: f(n)
 real(dp) :: ul(7),ur(7)

 call get_contact_states(x,n,ul,ur,.false.,.false.)
 select case(imode)
 case(imode_tangential)
    f(1) = (ptot(ul) - ptot(ur))/pscale
    f(2) = (ul(3) - ur(3))/vscale
 case(imode_coplanar)
    f(1) = (ul(2) - ur(2))/pscale
    f(2) = (ul(3) - ur(3))/vscale
    f(3) = dot_product(ul(4:5) - ur(4:5),edir)/vscale
    f(4) = dot_product(ul(6:7) - ur(6:7),edir)/bscale
 case default
    f(1) = (ul(2) - ur(2))/pscale
    f(2) = (ul(3) - ur(3))/vscale
    f(3) = (ul(4) - ur(4))/vscale
    f(4) = (ul(5) - ur(5))/vscale
    f(5) = (btrans(ul) - btrans(ur))/bscale
 end select

end subroutine residual

!-----------------------------------------------------------------------
!
! construct the solution nodes, ordered by xi = x/t, by applying
! the waves on the left side and then those on the right side
!
!-----------------------------------------------------------------------
subroutine build_solution(x,n)
 integer,  intent(in) :: n
 real(dp), intent(in) :: x(n)
 real(dp) :: ul(7),ur(7),vcontact
 integer  :: nleft

 nnodes = 0
 !--left-going waves are recorded in order of increasing xi
 call get_contact_states(x,n,ul,ur,.true.,.false.)
 vcontact = 0.5_dp*(ul(3) + ur(3))
 call add_node(vcontact,ul)
 nleft = nnodes
 !--right-going waves are recorded from the right state inwards,
 !  so reverse them to get increasing xi
 call get_contact_states(x,n,ul,ur,.false.,.true.)
 call add_node(vcontact,ur)
 xinode(nleft+1:nnodes)  = xinode(nnodes:nleft+1:-1)
 unode(:,nleft+1:nnodes) = unode(:,nnodes:nleft+1:-1)

end subroutine build_solution

!-----------------------------------------------------------------------
!
! states either side of the contact discontinuity, given the
! wave parameters x. The parameters are the strengths of the left
! fast and slow waves, the right fast and slow waves and, in the
! general case, the angle of the transverse field at the contact.
! recl, recr control whether the left/right waves are recorded as nodes
!
!-----------------------------------------------------------------------
subroutine get_contact_states(x,n,ul,ur,recl,recr)
 integer,  intent(in)  :: n
 real(dp), intent(in)  :: x(n)
 real(dp), intent(out) :: ul(7),ur(7)
 logical,  intent(in)  :: recl,recr
 real(dp) :: estar(2),theta

 ul = ul0
 ur = ur0
 select case(imode)
 case(imode_tangential)
    record = recl
    call fast_wave(ul,x(1),-1._dp,edir,smaxl)
    record = recr
    call fast_wave(ur,x(2), 1._dp,edir,smaxr)
 case(imode_coplanar)
    record = recl
    call fast_wave(ul,x(1),-1._dp,edir,smaxl)
    call slow_wave(ul,x(2),-1._dp,edir,.true.)
    record = recr
    call fast_wave(ur,x(3), 1._dp,edir,smaxr)
    call slow_wave(ur,x(4), 1._dp,edir,.true.)
 case default
    theta = atan2(edir(2),edir(1)) + x(5)
    estar = (/cos(theta),sin(theta)/)
    record = recl
    call fast_wave(ul,x(1),-1._dp,estar,smaxl)
    call rotational_wave(ul,estar,-1._dp)
    call slow_wave(ul,x(2),-1._dp,estar,.false.)
    record = recr
    call fast_wave(ur,x(3), 1._dp,estar,smaxr)
    call rotational_wave(ur,estar, 1._dp)
    call slow_wave(ur,x(4), 1._dp,estar,.false.)
 end select
 record = .false.

end subroutine get_contact_states

!-----------------------------------------------------------------------
! append a node to the solution
!-----------------------------------------------------------------------
subroutine add_node(xi,u)
 real(dp), intent(in) :: xi,u(7)

 if (nnodes >= maxnodes) return
 nnodes = nnodes + 1
 xinode(nnodes)  = xi
 unode(:,nnodes) = u

end subroutine add_node

!-----------------------------------------------------------------------
!
! fast wave with strength parameter psi. psi >= 0 gives a shock with
! upstream Mach number M0 = cf/a + psi; psi < 0 gives a rarefaction of
! length -psi, which ends in a switch-off rarefaction at smax (Takahashi
! & Yamada 2014 eq. 4.1). sig = -1 for left-going, +1 for right-going
!
!-----------------------------------------------------------------------
subroutine fast_wave(u,psi,sig,edef,smax)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: psi,sig,edef(2),smax
 real(dp) :: a2,ca2,cf2,cs2

 if (psi >= 0._dp) then
    call get_speeds(u,a2,ca2,cf2,cs2)
    call fast_shock(u,sqrt(cf2/a2)+psi,sig,edef)
 else
    call rarefaction(u,min(-psi,smax),sig,ifast,edef)
 endif

end subroutine fast_wave

!-----------------------------------------------------------------------
!
! slow wave with strength parameter psi. psi <= 0 gives a rarefaction
! of length -psi. In the general case psi > 0 gives a regular slow
! shock with downstream field b1 = b0 (1 - tanh(psi)) (Takahashi &
! Yamada 2014 eq. 3.4). In the coplanar case the downstream field is
! b1 = b0 - psi sqrt(p0), continuing past switch-off into intermediate
! shocks, and beyond the slow-sonic point into a compound wave
!
!-----------------------------------------------------------------------
subroutine slow_wave(u,psi,sig,edef,coplanar)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: psi,sig,edef(2)
 logical,  intent(in)    :: coplanar
 real(dp) :: b0,b1,bstar,sqrtp0,a2,ca2,cf2,cs2
 logical  :: ok

 b0 = btrans(u)
 if (b0 <= small_field(u)) then
    !--no transverse field: switch-on slow rarefaction or ordinary slow wave
    call get_speeds(u,a2,ca2,cf2,cs2)
    if (a2 > ca2) then
       call rarefaction(u,abs(psi),sig,islow,-sign(1._dp,psi)*edef)
    elseif (psi > 0._dp) then
       call euler_shock(u,1._dp + (sqrt(ca2/a2) - 1._dp)*tanh(psi),sig)
    else
       call rarefaction(u,-psi,sig,islow,edef)
    endif
    return
 endif
 if (psi <= 0._dp) then
    call rarefaction(u,-psi,sig,islow,edef)
    return
 endif

 sqrtp0 = sqrt(u(2))
 if (.not.coplanar) then
    call slow_shock(u,b0*(1._dp - tanh(psi)),sig,ok)
    return
 endif
 b1 = b0 - psi*sqrtp0
 if (b1 >= 0._dp) then
    call slow_shock(u,b1,sig,ok)
 else
    bstar = get_compound_point(u)
    if (b1 >= bstar) then
       call slow_shock(u,b1,sig,ok)
    else
       call slow_shock(u,bstar,sig,ok)
       call rarefaction(u,(bstar - b1)/sqrtp0,sig,islow,edef)
    endif
 endif

end subroutine slow_wave

!-----------------------------------------------------------------------
!
! rotational discontinuity rotating the transverse field to lie along
! the direction estar, keeping its magnitude (Alfven wave)
!
!-----------------------------------------------------------------------
subroutine rotational_wave(u,estar,sig)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: estar(2),sig
 real(dp) :: bt0(2),speed

 if (btrans(u) <= small_field(u)) return
 bt0 = u(6:7)
 speed = u(3) + sig*abs(bn)/sqrt(u(1))
 if (record) call add_node(speed,u)
 u(6:7) = btrans(u)*estar
 u(4:5) = u(4:5) - sig*sign(1._dp,bn)*(u(6:7) - bt0)/sqrt(u(1))
 if (record) call add_node(speed,u)

end subroutine rotational_wave

!-----------------------------------------------------------------------
!
! fast shock with upstream Mach number mach (relative to the sound speed).
! The compression ratio X = rho1/rho0 is the smallest root X > 1 of a
! cubic derived from the jump conditions (Takahashi & Yamada 2014 App. B)
! For zero upstream transverse field the shock is a switch-on shock
! where possible, or otherwise an ordinary (Euler) shock
!
!-----------------------------------------------------------------------
subroutine fast_shock(u,mach,sig,edef)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: mach,sig,edef(2)
 real(dp) :: rho0,p0,b0,q,w,x,b1,c(4),roots(3),a2,ca2,cf2,cs2
 integer  :: nroots

 rho0 = u(1)
 p0   = u(2)
 b0   = btrans(u)
 q    = bn**2
 call get_speeds(u,a2,ca2,cf2,cs2)
 w    = rho0*mach**2*a2   ! upstream momentum flux rho0 u0^2

 if (b0 <= small_field(u)) then
    !--switch-on shock, with the downstream flow speed equal to the Alfven speed
    x  = w/max(q,tiny(q))
    b1 = switch_on_field2(rho0,p0,w,x)
    if (ca2 > a2 .and. x > 1._dp .and. b1 > 0._dp) then
       call apply_shock(u,x,w,0._dp,sqrt(b1),sig,edef)
    else
       call euler_shock(u,mach,sig)
    endif
    return
 endif

 if (isothermal) then
    c(1) = 2._dp*csiso2*q**2*rho0
    c(2) = -w*(2._dp*b0**2*q - b0**2*w + 4._dp*csiso2*q*rho0 + 2._dp*q**2)
    c(3) = w**2*(b0**2 + 2._dp*csiso2*rho0 + 4._dp*q)
    c(4) = -2._dp*w**3
 else
    c(1) = 2._dp*q*(b0**2*gam*w - b0**2*w + 2._dp*gam*p0*q + gam*q*w - q*w)
    c(2) = -2._dp*w*(b0**2*gam*q + b0**2*gam*w + b0**2*q - 2._dp*b0**2*w + 4._dp*gam*p0*q &
                     + gam*q**2 + 2._dp*gam*q*w + q**2 - 2._dp*q*w)
    c(3) = 2._dp*w**2*(b0**2*gam + 2._dp*gam*p0 + 2._dp*gam*q + gam*w + 2._dp*q - w)
    c(4) = -2._dp*w**3*(gam + 1._dp)
 endif
 call get_cubic_roots(c,roots,nroots)
 x = minval(roots(1:nroots),mask=roots(1:nroots) > 1._dp)
 if (nroots == 0 .or. x > big/2._dp) return
 b1 = b0*(w - q)*x/(w - q*x)
 call apply_shock(u,x,w,b0,b1,sig,u(6:7)/b0)

end subroutine fast_shock

!-----------------------------------------------------------------------
!
! square of the downstream transverse field of a switch-on shock with
! compression ratio x, from the jump conditions with zero upstream field
!
!-----------------------------------------------------------------------
real(dp) function switch_on_field2(rho0,p0,w,x)
 real(dp), intent(in) :: rho0,p0,w,x
 real(dp) :: ptot1,v0,v1

 ptot1 = p0 + w*(1._dp - 1._dp/x)  ! downstream gas + magnetic pressure
 if (isothermal) then
    switch_on_field2 = 2._dp*(ptot1 - csiso2*rho0*x)
 else
    v0 = 1._dp/rho0
    v1 = v0/x
    switch_on_field2 = -2._dp*(gam - 1._dp)/v1*((p0*v0 - ptot1*v1)/(gam - 1._dp) &
                       + 0.5_dp*(p0 + ptot1)*(v0 - v1))
 endif

end function switch_on_field2

!-----------------------------------------------------------------------
!
! ordinary (hydrodynamic) shock with upstream Mach number mach,
! leaving the transverse field unchanged (zero)
!
!-----------------------------------------------------------------------
subroutine euler_shock(u,mach,sig)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: mach,sig
 real(dp) :: w,x,b0

 b0 = btrans(u)
 if (isothermal) then
    w = u(1)*mach**2*csiso2
    x = mach**2
 else
    w = gam*u(2)*mach**2
    x = (gam + 1._dp)*w/(2._dp*gam*u(2) + (gam - 1._dp)*w)
 endif
 if (x <= 1._dp) return
 call apply_shock(u,x,w,b0,b0,sig,edir)

end subroutine euler_shock

!-----------------------------------------------------------------------
!
! slow (or intermediate) shock giving downstream transverse field b1,
! measured along the upstream field direction. Writing X = 1 + (b1-b0) z,
! the jump conditions give a quadratic for z with no cancellation for
! weak shocks. The root with the smaller compression is taken, giving
! the "plus branch" of the slow Hugoniot locus (Torrilhon 2002)
!
!-----------------------------------------------------------------------
subroutine slow_shock(u,b1,sig,ok)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: b1,sig
 logical,  intent(out)   :: ok
 real(dp) :: x,w

 call slow_shock_jump(u,b1,x,w,ok)
 if (ok) call apply_shock(u,x,w,btrans(u),b1,sig,u(6:7)/btrans(u))

end subroutine slow_shock

!-----------------------------------------------------------------------
! compression ratio x and momentum flux w for a slow shock to field b1
!-----------------------------------------------------------------------
subroutine slow_shock_jump(u,b1,x,w,ok)
 real(dp), intent(in)  :: u(7),b1
 real(dp), intent(out) :: x,w
 logical,  intent(out) :: ok
 real(dp) :: rho0,p0,b0,q,d,a,b,c,disc,qq,z(2)
 integer  :: i

 rho0 = u(1)
 p0   = u(2)
 b0   = btrans(u)
 q    = bn**2
 d    = b1 - b0
 ok   = .false.
 x    = 1._dp
 w    = 0._dp
 if (abs(d) <= epsilon(d)*b0) return

 if (isothermal) then
    a = -2._dp*b0*csiso2*rho0
    b = -(2._dp*b0**2 + b0*d - 2._dp*csiso2*rho0 + 2._dp*q)
    c = 2._dp*b0 + d
 else
    a = -2._dp*((1._dp - gam)*b0**2*d + 2._dp*gam*b0*p0 + (1._dp - gam)*q*d)
    b = 2._dp*(-2._dp*b0**2 - gam*b0*d + 2._dp*gam*p0 - 2._dp*q)
    c = 2._dp*(2._dp*b0 + d)
 endif
 disc = b**2 - 4._dp*a*c
 if (disc < 0._dp .or. abs(a) <= tiny(a)) return
 qq = -0.5_dp*(b + sign(sqrt(disc),b))
 z(1) = qq/a
 z(2) = c/qq

 x = big
 do i=1,2
    if (1._dp + d*z(i) > 1._dp .and. 1._dp - b0*z(i) > 0._dp) then
       if (1._dp + d*z(i) < x) then
          x = 1._dp + d*z(i)
          w = q*x/(1._dp - b0*z(i))
          ok = .true.
       endif
    endif
 enddo

end subroutine slow_shock_jump

!-----------------------------------------------------------------------
!
! apply a shock with compression ratio x and upstream momentum flux
! w = rho0 u0^2 to the state u, recording nodes if required. b0 and b1
! are the upstream and downstream transverse fields along beta
!
!-----------------------------------------------------------------------
subroutine apply_shock(u,x,w,b0,b1,sig,beta)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: x,w,b0,b1,sig,beta(2)
 real(dp) :: speed

 if (record) call add_node(u(3) + sig*sqrt(w/u(1)),u)
 call downstream(u,x,w,b0,b1,sig,beta,speed)
 if (record) call add_node(speed,u)

end subroutine apply_shock

!-----------------------------------------------------------------------
!
! downstream state of a shock from the jump conditions in the shock
! frame, and the shock speed in the lab frame
!
!-----------------------------------------------------------------------
subroutine downstream(u,x,w,b0,b1,sig,beta,speed)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: x,w,b0,b1,sig,beta(2)
 real(dp), intent(out)   :: speed
 real(dp) :: rho0,un0

 rho0  = u(1)
 un0   = sqrt(w/rho0)
 speed = u(3) + sig*un0
 u(1) = rho0*x
 if (isothermal) then
    u(2) = csiso2*u(1)
 else
    u(2) = u(2) + w*(1._dp - 1._dp/x) + 0.5_dp*(b0**2 - b1**2)
 endif
 u(3)   = u(3) + sig*un0*(1._dp - 1._dp/x)
 u(4:5) = u(4:5) - sig*bn*(b1 - b0)*beta/(rho0*un0)
 u(6:7) = b1*beta

end subroutine downstream

!-----------------------------------------------------------------------
!
! downstream field b* on the slow Hugoniot locus where the downstream
! flow speed equals the slow speed (2 -> 3,4 intermediate shock),
! found by bisection. Beyond this point the slow wave is a compound
! wave (Takahashi & Yamada 2013 sec. 5.1)
!
!-----------------------------------------------------------------------
real(dp) function get_compound_point(u)
 real(dp), intent(in) :: u(7)
 real(dp) :: blo,bhi,bmid,b0
 integer  :: i

 b0  = btrans(u)
 bhi = 0._dp
 blo = -b0
 do i=1,100
    if (beyond_sonic(u,blo)) exit
    blo = 2._dp*blo
 enddo
 do i=1,100
    bmid = 0.5_dp*(blo + bhi)
    if (beyond_sonic(u,bmid)) then
       blo = bmid
    else
       bhi = bmid
    endif
 enddo
 get_compound_point = bhi

end function get_compound_point

!-----------------------------------------------------------------------
! true if the slow shock to field b1 is not possible, or has a
! downstream flow speed exceeding the downstream slow speed
!-----------------------------------------------------------------------
logical function beyond_sonic(u,b1)
 real(dp), intent(in) :: u(7),b1
 real(dp) :: u1(7),x,w,a2,ca2,cf2,cs2,speed
 logical  :: ok

 call slow_shock_jump(u,b1,x,w,ok)
 beyond_sonic = .true.
 if (.not.ok) return
 u1 = u
 call downstream(u1,x,w,btrans(u),b1,1._dp,u(6:7)/btrans(u),speed)
 call get_speeds(u1,a2,ca2,cf2,cs2)
 beyond_sonic = (w/(u(1)*x**2) > cs2)

end function beyond_sonic

!-----------------------------------------------------------------------
!
! rarefaction of length s in the fast or slow family, integrating
! du/ds = -r/rho along the right eigenvector r with the midpoint method.
! A fast rarefaction stops if the transverse field switches off
!
!-----------------------------------------------------------------------
subroutine rarefaction(u,s,sig,ifam,edef)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: s,sig,edef(2)
 integer,  intent(in)    :: ifam
 real(dp) :: ds,bdir0(2)
 logical  :: check_switchoff
 integer  :: i,ierr

 if (s <= 0._dp) return
 ds = s/nstep_raref
 check_switchoff = (ifam == ifast .and. btrans(u) > small_field(u))
 bdir0 = field_dir(u,edef)
 if (record) call add_node(wave_speed(u,sig,ifam),u)
 do i=1,nstep_raref
    call midpoint_step(u,ds,sig,ifam,edef,ierr)
    if (ierr /= 0) exit  ! vacuum reached
    if (check_switchoff .and. dot_product(u(6:7),bdir0) <= 0._dp) then
       u(6:7) = 0._dp
       exit
    endif
    if (record .and. mod(i,nstep_raref/nrec) == 0) call add_node(wave_speed(u,sig,ifam),u)
 enddo
 if (record) call add_node(wave_speed(u,sig,ifam),u)

end subroutine rarefaction

!-----------------------------------------------------------------------
!
! maximum length of a fast rarefaction before the transverse field
! switches off (Takahashi & Yamada 2014 sec. 5.2)
!
!-----------------------------------------------------------------------
real(dp) function get_smax(u0,sig,edef)
 real(dp), intent(in) :: u0(7),sig,edef(2)
 real(dp), parameter :: ds = 2.e-3_dp, smax_limit = 40._dp
 real(dp) :: u(7),uold(7),bdir0(2),s,dold,dnew
 integer  :: ierr
 real(dp) :: a2,ca2,cf2,cs2

 get_smax = big
 if (abs(bn) <= tiny(bn)) return
 if (btrans(u0) <= small_field(u0)) then
    call get_speeds(u0,a2,ca2,cf2,cs2)
    if (ca2 > a2) get_smax = 0._dp
    return
 endif
 u = u0
 bdir0 = u0(6:7)/btrans(u0)
 s = 0._dp
 do while (s < smax_limit)
    uold = u
    call midpoint_step(u,ds,sig,ifast,edef,ierr)
    if (ierr /= 0) return
    s = s + ds
    dnew = dot_product(u(6:7),bdir0)
    if (dnew <= 0._dp) then
       dold = dot_product(uold(6:7),bdir0)
       get_smax = s - ds*dnew/(dnew - dold)
       return
    endif
 enddo

end function get_smax

!-----------------------------------------------------------------------
!
! one midpoint step of length ds along a rarefaction. ierr = 1 and u
! is unchanged if the step would reach vacuum (rho or p <= 0)
!
!-----------------------------------------------------------------------
subroutine midpoint_step(u,ds,sig,ifam,edef,ierr)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: ds,sig,edef(2)
 integer,  intent(in)    :: ifam
 integer,  intent(out)   :: ierr
 real(dp) :: umid(7),unew(7),r(7)

 ierr = 1
 call raref_rhs(u,sig,ifam,edef,r)
 umid = u + 0.5_dp*ds*r
 if (umid(1) <= 0._dp .or. umid(2) <= 0._dp) return
 call raref_rhs(umid,sig,ifam,edef,r)
 unew = u + ds*r
 if (unew(1) <= 0._dp .or. unew(2) <= 0._dp) return
 if (isothermal) unew(2) = csiso2*unew(1)
 u = unew
 ierr = 0

end subroutine midpoint_step

!-----------------------------------------------------------------------
! right hand side of the rarefaction equations, du/ds = -r/rho
!-----------------------------------------------------------------------
subroutine raref_rhs(u,sig,ifam,edef,rhs)
 real(dp), intent(in)  :: u(7),sig,edef(2)
 integer,  intent(in)  :: ifam
 real(dp), intent(out) :: rhs(7)

 call eigenvector(u,sig,ifam,edef,rhs)
 rhs = -rhs/u(1)

end subroutine raref_rhs

!-----------------------------------------------------------------------
!
! right eigenvector for fast or slow waves in primitive variables
! (rho,p,vx,vy,vz,By,Bz), normalised as in Roe & Balsara (1996)
! so that it remains well defined in the degenerate limits
!
!-----------------------------------------------------------------------
subroutine eigenvector(u,sig,ifam,edef,r)
 real(dp), intent(in)  :: u(7),sig,edef(2)
 integer,  intent(in)  :: ifam
 real(dp), intent(out) :: r(7)
 real(dp) :: rho,a2,ca2,cf2,cs2,alphaf,alphas,beta(2),sgnbn,a,cf,cs

 rho = u(1)
 call get_speeds(u,a2,ca2,cf2,cs2)
 if (cf2 - cs2 <= 1.e-12_dp*cf2) then
    alphaf = sqrt(0.5_dp)
    alphas = sqrt(0.5_dp)
 else
    alphaf = sqrt(max(a2 - cs2,0._dp)/(cf2 - cs2))
    alphas = sqrt(max(cf2 - a2,0._dp)/(cf2 - cs2))
 endif
 beta  = field_dir(u,edef)
 sgnbn = sign(1._dp,bn)
 a  = sqrt(a2)
 cf = sqrt(cf2)
 cs = sqrt(cs2)
 if (ifam == ifast) then
    r(1)   = alphaf*rho
    r(2)   = alphaf*rho*a2
    r(3)   = sig*alphaf*cf
    r(4:5) = -sig*alphas*cs*beta*sgnbn
    r(6:7) = alphas*sqrt(rho)*a*beta
 else
    r(1)   = alphas*rho
    r(2)   = alphas*rho*a2
    r(3)   = sig*alphas*cs
    r(4:5) = sig*alphaf*cf*beta*sgnbn
    r(6:7) = -alphaf*sqrt(rho)*a*beta
 endif

end subroutine eigenvector

!-----------------------------------------------------------------------
! characteristic speed of the fast or slow wave
!-----------------------------------------------------------------------
real(dp) function wave_speed(u,sig,ifam)
 real(dp), intent(in) :: u(7),sig
 integer,  intent(in) :: ifam
 real(dp) :: a2,ca2,cf2,cs2

 call get_speeds(u,a2,ca2,cf2,cs2)
 if (ifam == ifast) then
    wave_speed = u(3) + sig*sqrt(cf2)
 else
    wave_speed = u(3) + sig*sqrt(cs2)
 endif

end function wave_speed

!-----------------------------------------------------------------------
! squared sound, Alfven, fast and slow speeds
!-----------------------------------------------------------------------
subroutine get_speeds(u,a2,ca2,cf2,cs2)
 real(dp), intent(in)  :: u(7)
 real(dp), intent(out) :: a2,ca2,cf2,cs2
 real(dp) :: tot

 if (isothermal) then
    a2 = csiso2
 else
    a2 = gam*u(2)/u(1)
 endif
 ca2 = bn**2/u(1)
 tot = a2 + (bn**2 + u(6)**2 + u(7)**2)/u(1)
 cf2 = 0.5_dp*(tot + sqrt(max(tot**2 - 4._dp*a2*ca2,0._dp)))
 cs2 = 0._dp
 if (cf2 > 0._dp) cs2 = a2*ca2/cf2

end subroutine get_speeds

!-----------------------------------------------------------------------
!
! real roots of the cubic c(1) x^3 + c(2) x^2 + c(3) x + c(4) = 0.
! cubicsolve (in default precision) gives one root, which is refined
! in double precision; the remaining roots come from the deflated quadratic
!
!-----------------------------------------------------------------------
subroutine get_cubic_roots(c,roots,nroots)
 real(dp), intent(in)  :: c(4)
 real(dp), intent(out) :: roots(3)
 integer,  intent(out) :: nroots
 real(dp) :: cn(4),x,f,df,dx,qa,qb,qc
 real     :: xs(3)
 integer  :: nreal,its

 roots  = big
 nroots = 0
 cn = c/maxval(abs(c))
 if (abs(cn(1)) < 1.e-12_dp) then
    call quadratic_roots(cn(2),cn(3),cn(4),roots,nroots)
    return
 endif
 call cubicsolve(real(cn(1)),real(cn(2)),real(cn(3)),real(cn(4)),xs,nreal)
 x = real(xs(1),dp)
 do its=1,50
    f  = ((cn(1)*x + cn(2))*x + cn(3))*x + cn(4)
    df = (3._dp*cn(1)*x + 2._dp*cn(2))*x + cn(3)
    if (abs(df) <= tiny(df)) exit
    dx = -f/df
    x  = x + dx
    if (abs(dx) <= 4._dp*epsilon(x)*abs(x)) exit
 enddo
 nroots = 1
 roots(1) = x
 !--deflate to a quadratic for the other two roots
 qa = cn(1)
 qb = cn(2) + qa*x
 qc = cn(3) + qb*x
 call quadratic_roots(qa,qb,qc,roots(2:3),nreal)
 nroots = nroots + nreal

end subroutine get_cubic_roots

!-----------------------------------------------------------------------
! real roots of a x^2 + b x + c = 0, avoiding cancellation
!-----------------------------------------------------------------------
subroutine quadratic_roots(a,b,c,roots,nroots)
 real(dp), intent(in)  :: a,b,c
 real(dp), intent(out) :: roots(:)
 integer,  intent(out) :: nroots
 real(dp) :: disc,qq

 roots = big
 nroots = 0
 if (abs(a) <= tiny(a)) then
    if (abs(b) > tiny(b)) then
       nroots = 1
       roots(1) = -c/b
    endif
    return
 endif
 disc = b**2 - 4._dp*a*c
 if (disc < 0._dp) return
 qq = -0.5_dp*(b + sign(sqrt(disc),b))
 nroots = 2
 roots(1) = qq/a
 if (abs(qq) > tiny(qq)) then
    roots(2) = c/qq
 else
    roots(2) = roots(1)
 endif

end subroutine quadratic_roots

!-----------------------------------------------------------------------
! magnitude of the transverse field and total pressure
!-----------------------------------------------------------------------
real(dp) function btrans(u)
 real(dp), intent(in) :: u(7)

 btrans = sqrt(u(6)**2 + u(7)**2)

end function btrans

real(dp) function ptot(u)
 real(dp), intent(in) :: u(7)

 ptot = u(2) + 0.5_dp*(u(6)**2 + u(7)**2)

end function ptot

!-----------------------------------------------------------------------
! threshold below which the transverse field is treated as zero
!-----------------------------------------------------------------------
real(dp) function small_field(u)
 real(dp), intent(in) :: u(7)

 small_field = 1.e-12_dp*sqrt(u(2))

end function small_field

!-----------------------------------------------------------------------
! unit vector along the transverse field, or edef if it is zero
!-----------------------------------------------------------------------
function field_dir(u,edef) result(beta)
 real(dp), intent(in) :: u(7),edef(2)
 real(dp) :: beta(2)

 if (btrans(u) > small_field(u)) then
    beta = u(6:7)/btrans(u)
 else
    beta = edef
 endif

end function field_dir

end module mhdshock
