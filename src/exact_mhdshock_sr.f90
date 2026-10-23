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
!  Copyright (C) 2026 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Exact solution of the Riemann problem in special relativistic
! magnetohydrodynamics (one dimensional SRMHD shock tubes), following
!
!  Giacomazzo & Rezzolla (2006), J. Fluid Mech. 562, 223
!
! As in that paper, the solution is regular: fast, Alfven and slow
! waves either side of a contact discontinuity (or two fast waves and
! a tangential discontinuity when Bx = 0), with no compound waves.
!
! Each wave is parameterised by a single strength and the strengths
! are iterated until the total pressure, velocity and transverse field
! match at the contact. Shocks and Alfven discontinuities are found
! by solving the jump conditions [[F]] = Vs [[U]] numerically,
! with one extra condition setting the strength: the jump in total
! pressure (fast shocks), the downstream transverse field along the
! upstream field direction (slow shocks) or the rotation angle of
! the transverse field (Alfven waves). Rarefactions are integrated
! along the right eigenvectors of the flux Jacobian, with the
! magnetosonic speeds from the quartic dispersion relation
! (e.g. Mignone & Bodo 2006, MNRAS 368, 1040)
!
! The preset problems are from Komissarov (1999), Balsara (2001)
! and Giacomazzo & Rezzolla (2006)
!
! Daniel Price, Monash University, 2026
! daniel.price@monash.edu
!
!-----------------------------------------------------------------------
module mhdshock_sr
 use nonlinsolve, only:newton_nd,solve_linear
 implicit none
 public :: exact_mhdshock_sr,get_mhdshock_sr_preset,solve_mhdshock_sr,total_pressure_sr
 integer, parameter, public :: nmhdshocksrsolns = 10
 character(len=40), dimension(nmhdshocksrsolns), parameter, public :: mhdsrprob = &
    (/'Komissarov (1999) shock tube 1, G=4/3  ', &
      'Komissarov (1999) shock tube 2, G=4/3  ', &
      'Komissarov (1999) collision, G=4/3     ', &
      'Balsara (2001) test 1, G=2             ', &
      'Balsara (2001) test 2, G=5/3           ', &
      'Balsara (2001) test 3, G=5/3           ', &
      'Balsara (2001) test 4, G=5/3           ', &
      'Balsara (2001) test 5, G=5/3           ', &
      'generic Alfven test (GR06), G=5/3      ', &
      'generic shock tube, Bx=0 (GR06), G=5/3 '/)

 private

 integer, parameter :: dp = kind(0.d0)
 integer, parameter :: ifast = 1, ialfven = 2, islow = 3
 integer, parameter :: icon_ptot = 1, icon_bpar = 2, icon_angle = 3
 integer, parameter :: nstep_raref = 100   ! RK4 integration steps per rarefaction
 integer, parameter :: nrec = 100          ! nodes recorded per rarefaction
 integer, parameter :: maxnodes = 4096
 real(dp), parameter :: tol_newton = 1.e-10_dp
 real(dp), parameter :: big = huge(1._dp)
 real(dp), parameter :: bad_residual = 1.e30_dp

 !--settings for the current Riemann problem, shared with the residual function
 real(dp) :: gam,bn,ul0(7),ur0(7),pscale
 logical  :: tangential,wave_error

 !--last solution of each discontinuity (fast-, Alfven-, slow-, slow+, Alfven+,
 !  fast+), used as the initial guess when the wave strengths change slightly
 logical  :: cache_ok(6)
 real(dp) :: cache_ua(7,6),cache_ub(7,6),cache_vs(6),cache_target(6)

 !--inputs for the last solution, which is reused if the inputs are unchanged
 logical  :: have_solution = .false.
 real(dp) :: last_inputs(16)

 !--solution nodes (xi = x/t, primitive state) and states in each region
 logical  :: record = .false.
 integer  :: nnodes,nnodes_left,nregion
 real(dp) :: xinode(maxnodes),unode(7,maxnodes),uregion(7,8)

contains

!-----------------------------------------------------------------------
!
! plot the exact solution of the SRMHD Riemann problem for left and
! right states (rho,p,vx,vy,vz,By,Bz), where rho is the rest-frame
! density and p the gas pressure, with normal field Bx
!
! input  : iplot : variable to plot (1=rho,2=p,3=vx,4=vy,5=vz,6=By,7=Bz,
!                  8=utherm,9=Bx,10=rho*=rho*W)
! output : xpts,ypts,npts : the solution as a piecewise linear curve
!
!-----------------------------------------------------------------------
subroutine exact_mhdshock_sr(iplot,time,gamma,bx,primL,primR,xmin,xmax,xshock,xpts,ypts,npts,ierr)
 integer, intent(in)  :: iplot
 integer, intent(out) :: npts,ierr
 real,    intent(in)  :: time,gamma,bx,primL(7),primR(7),xmin,xmax,xshock
 real, dimension(:), intent(inout) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 real(dp), allocatable :: uplot(:,:),xplot(:)
 integer  :: i,n

 print "(a,es10.3,a)",' Plotting exact SRMHD shock tube at t = ',time,' from Giacomazzo & Rezzolla (2006)'
 npts = 0
 xpts = 0.
 ypts = 0.
 call check_states(primL,primR,ierr)
 if (ierr /= 0) return

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
    call solve_riemann_sr(real(gamma,dp),real(bx,dp),real(primL,dp),real(primR,dp),ierr)
    if (ierr /= 0) then
       print*,'WARNING: SRMHD Riemann solver failed to converge, ierr = ',ierr
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
    ypts(1:npts) = real(uplot(2,1:npts)/((gamma - 1.)*uplot(1,1:npts)))
 case(9)
    ypts(1:npts) = bx
 case(10)
    ypts(1:npts) = real(uplot(1,1:npts)/sqrt(1._dp - sum(uplot(3:5,1:npts)**2,dim=1)))
 case default
    print*,'error: unknown solution to plot'
    ierr = 2
 end select

end subroutine exact_mhdshock_sr

!-----------------------------------------------------------------------
!
! solve the SRMHD Riemann problem and return the states in each
! region of the Riemann fan, R1-R8 when Bx /= 0 (R1-R4 when Bx = 0),
! in the order of Giacomazzo & Rezzolla (2006) Figs. 1 and 2
!
!-----------------------------------------------------------------------
subroutine solve_mhdshock_sr(gamma,bx,primL,primR,nreg,ureg,ierr)
 real,    intent(in)  :: gamma,bx,primL(7),primR(7)
 integer, intent(out) :: nreg,ierr
 real,    intent(out) :: ureg(7,8)

 nreg = 0
 ureg = 0.
 call check_states(primL,primR,ierr)
 if (ierr /= 0) return
 call solve_riemann_sr(real(gamma,dp),real(bx,dp),real(primL,dp),real(primR,dp),ierr)
 if (ierr /= 0) return
 nreg = nregion
 ureg(:,1:nreg) = real(uregion(:,1:nreg))

end subroutine solve_mhdshock_sr

!-----------------------------------------------------------------------
! total (gas + magnetic) pressure for primitive state prim
!-----------------------------------------------------------------------
real function total_pressure_sr(prim,bx)
 real, intent(in) :: prim(7),bx
 real(dp) :: bn_saved

 bn_saved = bn
 bn = bx
 total_pressure_sr = real(ptot(real(prim,dp)))
 bn = bn_saved

end function total_pressure_sr

!-----------------------------------------------------------------------
!
! initial left and right states (rho,p,vx,vy,vz,By,Bz) and Bx for
! the preset problems. The adiabatic index for each is given in mhdsrprob
!
!-----------------------------------------------------------------------
subroutine get_mhdshock_sr_preset(ishk,bx,primL,primR,ierr)
 integer, intent(in)  :: ishk
 real,    intent(out) :: bx,primL(7),primR(7)
 integer, intent(out) :: ierr

 ierr = 0
 select case(ishk)
 case(1) ! Komissarov (1999) shock tube 1
    bx = 1.
    primL = (/1.,1000.,0.,0.,0.,0.,0./)
    primR = (/0.1,1.,0.,0.,0.,0.,0./)
 case(2) ! Komissarov (1999) shock tube 2
    bx = 0.
    primL = (/1.,30.,0.,0.,0.,20.,0./)
    primR = (/0.1,1.,0.,0.,0.,0.,0./)
 case(3) ! Komissarov (1999) collision
    bx = 10.
    primL = (/1.,1., 5./sqrt(26.),0.,0., 10.,0./)
    primR = (/1.,1.,-5./sqrt(26.),0.,0.,-10.,0./)
 case(4) ! Balsara (2001) test 1, relativistic version of Brio & Wu
    bx = 0.5
    primL = (/1.,1.,0.,0.,0.,1.,0./)
    primR = (/0.125,0.1,0.,0.,0.,-1.,0./)
 case(5) ! Balsara (2001) test 2
    bx = 5.
    primL = (/1.,30.,0.,0.,0.,6.,6./)
    primR = (/1.,1.,0.,0.,0.,0.7,0.7/)
 case(6) ! Balsara (2001) test 3
    bx = 10.
    primL = (/1.,1000.,0.,0.,0.,7.,7./)
    primR = (/1.,0.1,0.,0.,0.,0.7,0.7/)
 case(7) ! Balsara (2001) test 4
    bx = 10.
    primL = (/1.,0.1, 0.999,0.,0., 7., 7./)
    primR = (/1.,0.1,-0.999,0.,0.,-7.,-7./)
 case(8) ! Balsara (2001) test 5
    bx = 2.
    primL = (/1.08,0.95,0.4,0.3,0.2,0.3,0.3/)
    primR = (/1.,1.,-0.45,-0.2,0.2,-0.7,0.5/)
 case(9) ! generic Alfven test (Giacomazzo & Rezzolla 2006)
    bx = 1.
    primL = (/1.,5.,0.,0.3,0.4,6.,2./)
    primR = (/0.9,5.3,0.,0.,0.,5.,2./)
 case(10) ! generic shock tube with Bx = 0 (Giacomazzo & Rezzolla 2006)
    bx = 0.
    primL = (/1.,0.01,0.1,0.3,0.4,6.,2./)
    primR = (/0.01,5000.,0.5,0.4,0.3,5.,20./)
 case default
    print*,'error: unknown SRMHD shock tube ',ishk
    bx = 0.
    primL = 0.
    primR = 0.
    ierr = 1
 end select

end subroutine get_mhdshock_sr_preset

!-----------------------------------------------------------------------
! check that the initial states are physical
!-----------------------------------------------------------------------
subroutine check_states(primL,primR,ierr)
 real,    intent(in)  :: primL(7),primR(7)
 integer, intent(out) :: ierr

 ierr = 0
 if (primL(1) <= 0. .or. primR(1) <= 0. .or. primL(2) <= 0. .or. primR(2) <= 0.) then
    print*,'error: rho or p <= 0 in initial states'
    ierr = 1
 elseif (sum(primL(3:5)**2) >= 1. .or. sum(primR(3:5)**2) >= 1.) then
    print*,'error: |v| >= 1 in initial states'
    ierr = 1
 endif

end subroutine check_states

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
! solve the Riemann problem, returning the solution in the module
! variables xinode, unode, nnodes and the region states uregion
!
!-----------------------------------------------------------------------
subroutine solve_riemann_sr(gamma,bx,primL,primR,ierr)
 real(dp), intent(in)  :: gamma,bx,primL(7),primR(7)
 integer,  intent(out) :: ierr
 real(dp) :: x(6)
 integer  :: n,its,istep,nsteps,itry
 integer, parameter :: nsteps_try(2) = (/10,40/)

 ierr = 0
 if (have_solution) then
    if (all(abs(last_inputs - (/gamma,bx,primL,primR/)) <= 0._dp)) return
 endif
 have_solution = .false.
 gam = gamma
 bn  = bx
 ul0 = primL
 ur0 = primR
 call setup_problem(n)

 x = 0._dp
 call newton_nd(residual,x(1:n),n,tol_newton,100,its,ierr)

 !--if this fails, march the right state from the left state to its final value
 do itry=1,size(nsteps_try)
    if (ierr == 0) exit
    nsteps = nsteps_try(itry)
    x = 0._dp
    do istep=1,nsteps
       ur0 = interpolate_state(primL,primR,real(istep,dp)/real(nsteps,dp))
       call setup_problem(n)
       call newton_nd(residual,x(1:n),n,tol_newton,100,its,ierr)
       if (ierr /= 0) exit
    enddo
 enddo
 if (ierr /= 0) return

 call build_solution(x(1:n),n)
 have_solution = .true.
 last_inputs = (/gamma,bx,primL,primR/)

end subroutine solve_riemann_sr

!-----------------------------------------------------------------------
!
! choose between the full solution and the reduced solution with only
! fast waves, used when Bx = 0 or when there are no transverse fields
! or velocity jumps (the hydrodynamic limit)
!
!-----------------------------------------------------------------------
subroutine setup_problem(n)
 integer, intent(out) :: n
 logical :: notrans

 cache_ok = .false.
 pscale = max(ptot(ul0),ptot(ur0))
 notrans = all(abs(ul0(6:7)) <= 0._dp) .and. all(abs(ur0(6:7)) <= 0._dp) .and. &
           all(abs(ul0(4:5) - ur0(4:5)) <= 0._dp)
 tangential = (abs(bn) <= 0._dp .or. notrans)
 if (tangential) then
    n = 2
 else
    n = 6
 endif

end subroutine setup_problem

!-----------------------------------------------------------------------
!
! intermediate state between uL and uR for the continuation method,
! interpolating the transverse field in magnitude and angle so that
! it does not pass through zero
!
!-----------------------------------------------------------------------
function interpolate_state(uL,uR,f) result(u)
 real(dp), intent(in) :: uL(7),uR(7),f
 real(dp) :: u(7),bl,br,thl,thr,dth,th

 u(1:5) = uL(1:5) + f*(uR(1:5) - uL(1:5))
 bl  = norm2(uL(6:7))
 br  = norm2(uR(6:7))
 thl = atan2(uL(7),uL(6))
 thr = atan2(uR(7),uR(6))
 if (bl <= 0._dp) thl = thr
 if (br <= 0._dp) thr = thl
 dth = modulo(thr - thl + acos(-1._dp),2._dp*acos(-1._dp)) - acos(-1._dp)
 if (abs(abs(dth) - acos(-1._dp)) < 1.e-8_dp) dth = acos(-1._dp)
 th  = thl + f*dth
 u(6:7) = (bl + f*(br - bl))*(/cos(th),sin(th)/)

end function interpolate_state

!-----------------------------------------------------------------------
!
! residual of the matching conditions at the contact discontinuity
! given the wave parameters x = (fast-, Alfven-, slow-, slow+, Alfven+,
! fast+), or x = (fast-, fast+) for the reduced problem
!
!-----------------------------------------------------------------------
subroutine residual(x,f,n)
 integer,  intent(in)  :: n
 real(dp), intent(in)  :: x(n)
 real(dp), intent(out) :: f(n)
 real(dp) :: ul(7),ur(7)

 call get_contact_states(x,n,ul,ur)
 if (wave_error) then
    f = bad_residual
    return
 endif
 f(1) = (ptot(ul) - ptot(ur))/pscale
 f(2) = ul(3) - ur(3)
 if (n > 2) then
    f(3:4) = ul(4:5) - ur(4:5)
    f(5:6) = (ul(6:7) - ur(6:7))/sqrt(pscale)
 endif

end subroutine residual

!-----------------------------------------------------------------------
!
! states either side of the contact discontinuity given the wave
! parameters x, recording the region states (and nodes if required)
!
!-----------------------------------------------------------------------
subroutine get_contact_states(x,n,ul,ur)
 integer,  intent(in)  :: n
 real(dp), intent(in)  :: x(n)
 real(dp), intent(out) :: ul(7),ur(7)
 real(dp) :: xfull(6)

 if (n == 2) then
    xfull = (/x(1),0._dp,0._dp,0._dp,0._dp,x(2)/)
 else
    xfull = x
 endif
 wave_error = .false.
 ul = ul0
 ur = ur0
 uregion(:,1) = ul0
 uregion(:,8) = ur0
 call fast_wave(ul,xfull(1),-1._dp)
 uregion(:,2) = ul
 call alfven_wave(ul,xfull(2),-1._dp)
 uregion(:,3) = ul
 call slow_wave(ul,xfull(3),-1._dp)
 uregion(:,4) = ul
 nnodes_left = nnodes

 !--right-going waves are recorded separately so they can be reversed
 call fast_wave(ur,xfull(6),1._dp)
 uregion(:,7) = ur
 call alfven_wave(ur,xfull(5),1._dp)
 uregion(:,6) = ur
 call slow_wave(ur,xfull(4),1._dp)
 uregion(:,5) = ur

end subroutine get_contact_states

!-----------------------------------------------------------------------
!
! evaluate the solution recording the nodes, ordered by xi = x/t,
! and the states in each region of the Riemann fan
!
!-----------------------------------------------------------------------
subroutine build_solution(x,n)
 integer,  intent(in) :: n
 real(dp), intent(in) :: x(n)
 real(dp) :: ul(7),ur(7),vcontact
 real(dp), allocatable :: xinode_r(:),unode_r(:,:)
 integer  :: nleft,nright

 !--evaluate the waves on each side, recording the nodes
 record = .true.
 nnodes = 0
 call get_contact_states(x,n,ul,ur)
 record = .false.
 !--split the nodes: left-going waves are first, right-going waves are
 !  recorded from the right state inwards
 nleft = nnodes_left
 nright = nnodes - nleft
 allocate(xinode_r(nright),unode_r(7,nright))
 xinode_r(1:nright)  = xinode(nnodes:nleft+1:-1)
 unode_r(:,1:nright) = unode(:,nnodes:nleft+1:-1)
 nnodes = nleft
 vcontact = 0.5_dp*(ul(3) + ur(3))
 call add_node(vcontact,ul)
 call add_node(vcontact,ur)
 xinode(nnodes+1:nnodes+nright)  = xinode_r(1:nright)
 unode(:,nnodes+1:nnodes+nright) = unode_r(:,1:nright)
 nnodes = nnodes + nright

 if (tangential .and. abs(bn) <= 0._dp) then
    nregion = 4
    uregion(:,2) = uregion(:,4)
    uregion(:,3) = uregion(:,5)
    uregion(:,4) = uregion(:,8)
 else
    nregion = 8
 endif

end subroutine build_solution

!-----------------------------------------------------------------------
!
! fast wave with strength psi: psi > 0 gives a shock with total
! pressure jump ln(ptot_b/ptot_a) = psi, psi < 0 a rarefaction of
! length -psi. sig = -1 for left-going, +1 for right-going waves
!
!-----------------------------------------------------------------------
subroutine fast_wave(u,psi,sig)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: psi,sig

 if (wave_error) return
 if (psi > 0._dp) then
    call jump(u,ifast,sig,icon_ptot,psi)
 elseif (psi < 0._dp) then
    call rarefaction(u,-psi,ifast,sig)
 endif

end subroutine fast_wave

!-----------------------------------------------------------------------
! Alfven discontinuity rotating the transverse field by angle phi
!-----------------------------------------------------------------------
subroutine alfven_wave(u,phi,sig)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: phi,sig

 if (wave_error) return
 if (abs(phi) > 0._dp) call jump(u,ialfven,sig,icon_angle,phi)

end subroutine alfven_wave

!-----------------------------------------------------------------------
!
! slow wave with strength psi: psi > 0 gives a shock with downstream
! transverse field along the upstream field direction b0 - psi sqrt(ptot),
! which may pass through zero, psi < 0 a rarefaction of length -psi
!
!-----------------------------------------------------------------------
subroutine slow_wave(u,psi,sig)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: psi,sig

 if (wave_error) return
 if (psi > 0._dp) then
    call jump(u,islow,sig,icon_bpar,norm2(u(6:7)) - psi*sqrt(ptot(u)))
 elseif (psi < 0._dp) then
    call rarefaction(u,-psi,islow,sig)
 endif

end subroutine slow_wave

!-----------------------------------------------------------------------
!
! shock or Alfven discontinuity in family ifam, with strength set by
! the constraint icon = target. The strength is increased by
! continuation, doubling the step after each success and halving it
! after a failure, starting either from the last solution for this
! wave (if the upstream state is similar) or from the linearised
! solution for a weak discontinuity
!
!-----------------------------------------------------------------------
subroutine jump(u,ifam,sig,icon,target)
 real(dp), intent(inout) :: u(7)
 integer,  intent(in)    :: ifam,icon
 real(dp), intent(in)    :: sig,target
 real(dp), parameter :: eps = 1.e-7_dp
 real(dp) :: ua(7),ub(7),r(7),rnone(7),lam,vs,ref(3),t0,tk,c0,c1
 integer  :: ierr,islot

 ua = u
 ref(1) = ptot(ua)
 ref(2) = atan2(ua(7),ua(6))
 ref(3) = sqrt(ref(1))
 t0 = 0._dp
 if (icon == icon_bpar) t0 = norm2(ua(6:7))
 islot = ifam
 if (sig > 0._dp) islot = 7 - ifam

 !--first try starting from the last solution for this wave, if similar
 if (cache_ok(islot)) then
    if (maxval(abs(ua - cache_ua(:,islot))/state_scale(ua)) < 0.05_dp .and. &
        abs(target - cache_target(islot)) < 0.05_dp*max(abs(target - t0),1.e-3_dp)) then
       ub = cache_ub(:,islot)
       vs = cache_vs(islot)
       call jump_newton(ua,icon,ref,target,vs,ub,ierr)
       if (ierr == 0) then
          call finish_jump(u,ua,ub,vs,target,islot)
          return
       endif
    endif
 endif

 !--otherwise start from the linearised solution for a weak discontinuity
 rnone = 0._dp
 call eigenvector(ua,ifam,sig,rnone,.false.,r,lam,ierr)
 if (ierr /= 0) then
    wave_error = .true.
    return
 endif
 tk = t0 + (target - t0)*1.e-3_dp
 c0 = constraint(icon,ua,ref,tk)
 c1 = constraint(icon,ua - eps*r,ref,tk)
 if (abs(c1 - c0) <= tiny(c0)) then
    wave_error = .true.
    return
 endif
 ub = ua - (c0*eps/(c1 - c0))*r
 vs = lam
 call continuation(ua,icon,ref,t0,target,1.e-3_dp,vs,ub,ierr)
 if (ierr /= 0) then
    wave_error = .true.
    return
 endif
 call finish_jump(u,ua,ub,vs,target,islot)

end subroutine jump

!-----------------------------------------------------------------------
!
! solve the jump conditions for constraint values from tstart to
! target, starting with a step df (as a fraction of the total), doubling
! the step after each success and halving it after a failure. For fast
! shocks the initial guess for each step is extrapolated from the
! previous two solutions
!
!-----------------------------------------------------------------------
subroutine continuation(ua,icon,ref,tstart,target,dfstart,vs,ub,ierr)
 real(dp), intent(in)    :: ua(7),ref(3),tstart,target,dfstart
 integer,  intent(in)    :: icon
 real(dp), intent(inout) :: vs,ub(7)
 integer,  intent(out)   :: ierr
 real(dp), parameter :: df_min = 1.e-8_dp
 real(dp) :: ubtry(7),vstry,frac,df,tk,ubold(7),vsold,fracold,fac
 logical  :: have_old

 frac = 0._dp
 df = dfstart
 ierr = 0
 have_old = .false.
 do while (frac < 1._dp)
    df = min(df,1._dp - frac)
    tk = tstart + (target - tstart)*(frac + df)
    ubtry = ub
    vstry = vs
    if (have_old .and. icon == icon_ptot) then
       !--secant predictor from the previous two solutions, used only for
       !  fast shocks where there is a single branch of solutions
       fac = df/(frac - fracold)
       ubtry = ub + fac*(ub - ubold)
       vstry = vs + fac*(vs - vsold)
       if (.not.valid_jump((/vstry,ubtry/))) then
          ubtry = ub
          vstry = vs
       endif
    endif
    call jump_newton(ua,icon,ref,tk,vstry,ubtry,ierr)
    if (ierr == 0) then
       ubold = ub
       vsold = vs
       fracold = frac
       have_old = .true.
       frac = frac + df
       ub = ubtry
       vs = vstry
       df = 2._dp*df
    else
       df = 0.5_dp*df
       if (df < df_min) return
    endif
 enddo
 ierr = 0

end subroutine continuation

!-----------------------------------------------------------------------
! store the solution of a discontinuity, recording nodes if required
!-----------------------------------------------------------------------
subroutine finish_jump(u,ua,ub,vs,target,islot)
 real(dp), intent(out) :: u(7)
 real(dp), intent(in)  :: ua(7),ub(7),vs,target
 integer,  intent(in)  :: islot

 if (record) then
    call add_node(vs,ua)
    call add_node(vs,ub)
 endif
 cache_ok(islot) = .true.
 cache_ua(:,islot) = ua
 cache_ub(:,islot) = ub
 cache_vs(islot) = vs
 cache_target(islot) = target
 u = ub

end subroutine finish_jump

!-----------------------------------------------------------------------
!
! the extra condition fixing the strength of a discontinuity with
! downstream state ub. ref = (ptot, field angle, sqrt(ptot)) upstream
!
!-----------------------------------------------------------------------
real(dp) function constraint(icon,ub,ref,target)
 integer,  intent(in) :: icon
 real(dp), intent(in) :: ub(7),ref(3),target

 select case(icon)
 case(icon_ptot)
    constraint = log(ptot(ub)/ref(1)) - target
 case(icon_bpar)
    constraint = (ub(6)*cos(ref(2)) + ub(7)*sin(ref(2)) - target)/ref(3)
 case default
    constraint = (-ub(6)*sin(ref(2)+target) + ub(7)*cos(ref(2)+target))/ref(3)
 end select

end function constraint

!-----------------------------------------------------------------------
!
! solve the jump conditions [[F]] = vs [[U]] across a discontinuity
! together with the constraint for the downstream state ub and the
! speed vs, by damped Newton-Raphson iteration from the initial guess
!
!-----------------------------------------------------------------------
subroutine jump_newton(ua,icon,ref,target,vs,ub,ierr)
 real(dp), intent(in)    :: ua(7),ref(3),target
 integer,  intent(in)    :: icon
 real(dp), intent(inout) :: vs,ub(7)
 integer,  intent(out)   :: ierr
 integer, parameter :: maxits = 20, maxhalve = 40
 real(dp) :: z(8),zn(8),zp(8),f(8),fn(8),fp(8),dz(8),jac(8,8)
 real(dp) :: uca(7),fa(7),sc(7),fnorm,lambda,h
 integer  :: its,j,ihalve

 call get_cons_flux(ua,uca,fa)
 sc = max(abs(uca),abs(fa)) + 1.e-10_dp*(sum(abs(uca)) + sum(abs(fa)))
 z = (/vs,ub/)
 call jump_residual(z,uca,fa,sc,icon,ref,target,f)
 ierr = 0
 do its=1,maxits
    fnorm = norm2(f)
    if (fnorm < 1.e-13_dp) exit
    do j=1,8
       h  = 1.e-8_dp*max(1._dp,abs(z(j)))
       zp = z
       zp(j) = z(j) + h
       call jump_residual(zp,uca,fa,sc,icon,ref,target,fp)
       jac(:,j) = (fp - f)/h
    enddo
    dz = -f
    call solve_linear(jac,dz,8,ierr)
    if (ierr /= 0) exit
    lambda = 1._dp
    do ihalve=1,maxhalve
       zn = z + lambda*dz
       if (valid_jump(zn)) then
          call jump_residual(zn,uca,fa,sc,icon,ref,target,fn)
          if (norm2(fn) < fnorm) exit
       endif
       lambda = 0.5_dp*lambda
    enddo
    if (ihalve > maxhalve) exit
    z = zn
    f = fn
    if (maxval(abs(lambda*dz)) < 1.e-14_dp*max(1._dp,maxval(abs(z)))) exit
 enddo
 ierr = 0
 if (norm2(f) > 1.e-9_dp) ierr = 2
 vs = z(1)
 ub = z(2:8)

end subroutine jump_newton

!-----------------------------------------------------------------------
! residual of the jump conditions and constraint, z = (vs,ub)
!-----------------------------------------------------------------------
subroutine jump_residual(z,uca,fa,sc,icon,ref,target,f)
 real(dp), intent(in)  :: z(8),uca(7),fa(7),sc(7),ref(3),target
 integer,  intent(in)  :: icon
 real(dp), intent(out) :: f(8)
 real(dp) :: ucb(7),fb(7)

 if (.not.valid_jump(z)) then
    f = bad_residual
    return
 endif
 call get_cons_flux(z(2:8),ucb,fb)
 f(1:7) = ((fb - fa) - z(1)*(ucb - uca))/sc
 f(8)   = constraint(icon,z(2:8),ref,target)

end subroutine jump_residual

!-----------------------------------------------------------------------
! check that the shock speed and downstream state are physical
!-----------------------------------------------------------------------
logical function valid_jump(z)
 real(dp), intent(in) :: z(8)

 valid_jump = (abs(z(1)) < 1._dp .and. valid_state(z(2:8)))

end function valid_jump

!-----------------------------------------------------------------------
!
! rarefaction of length s in the fast or slow family, integrating
! du/ds = r along the right eigenvector r (normalised, with density
! decreasing) using the fourth order Runge-Kutta method
!
!-----------------------------------------------------------------------
subroutine rarefaction(u,s,ifam,sig)
 real(dp), intent(inout) :: u(7)
 real(dp), intent(in)    :: s,sig
 integer,  intent(in)    :: ifam
 real(dp) :: ds,r(7),k1(7),k2(7),k3(7),k4(7),un(7),lam
 logical  :: have_r
 integer  :: i,ierr

 ds = s/nstep_raref
 have_r = .false.
 r = 0._dp
 do i=1,nstep_raref
    call eigenvector(u,ifam,sig,r,have_r,k1,lam,ierr)
    if (record) call add_node(lam,u)
    if (ierr == 0) call raref_stage(u,0.5_dp*ds,k1,ifam,sig,k2,ierr)
    if (ierr == 0) call raref_stage(u,0.5_dp*ds,k2,ifam,sig,k3,ierr)
    if (ierr == 0) call raref_stage(u,ds,k3,ifam,sig,k4,ierr)
    if (ierr > 0) wave_error = .true.  ! no eigenvector found
    if (ierr /= 0) exit                ! or vacuum reached (ierr < 0)
    un = u + ds*(k1 + 2._dp*k2 + 2._dp*k3 + k4)/6._dp
    if (.not.valid_state(un)) exit
    u = un
    r = k4
    have_r = .true.
 enddo
 if (record) call add_node(wave_speed(u,ifam,sig),u)

end subroutine rarefaction

!-----------------------------------------------------------------------
! eigenvector at u + h*k for a Runge-Kutta stage, oriented along k.
! ierr = -1 if the stage reaches vacuum, > 0 if the eigenvector fails
!-----------------------------------------------------------------------
subroutine raref_stage(u,h,k,ifam,sig,knew,ierr)
 real(dp), intent(in)  :: u(7),h,k(7),sig
 integer,  intent(in)  :: ifam
 real(dp), intent(out) :: knew(7)
 integer,  intent(out) :: ierr
 real(dp) :: ustage(7),lam

 ierr = -1
 knew = k
 ustage = u + h*k
 if (.not.valid_state(ustage)) return
 call eigenvector(ustage,ifam,sig,k,.true.,knew,lam,ierr)

end subroutine raref_stage

!-----------------------------------------------------------------------
!
! right eigenvector r of the flux Jacobian for the wave family ifam
! moving in direction sig, found by inverse iteration on
! dF/du - lambda dU/du. r is normalised and oriented either along rprev
! (if have_prev) or so that the density decreases along r
!
!-----------------------------------------------------------------------
subroutine eigenvector(u,ifam,sig,rprev,have_prev,r,lam,ierr)
 real(dp), intent(in)  :: u(7),sig,rprev(7)
 integer,  intent(in)  :: ifam
 logical,  intent(in)  :: have_prev
 real(dp), intent(out) :: r(7),lam
 integer,  intent(out) :: ierr
 real(dp) :: dudw(7,7),dfdw(7,7),m(7,7),mc(7,7),y(7),sc(7)
 integer  :: it,i

 lam = wave_speed(u,ifam,sig)
 call get_jacobians(u,dudw,dfdw)
 m  = dfdw - lam*dudw
 sc = state_scale(u)
 if (have_prev) then
    r = rprev
 else
    r = (/1._dp,1.1_dp,0.9_dp,1.2_dp,0.8_dp,1.3_dp,0.7_dp/)*sc
 endif
 do it=1,merge(1,2,have_prev)
    mc = m
    y  = r
    call solve_linear(mc,y,7,ierr)
    if (ierr /= 0) then
       !--exactly singular: shift the matrix slightly
       mc = m
       do i=1,7
          mc(i,i) = mc(i,i) + 1.e-13_dp*maxval(abs(m))
       enddo
       y = r
       call solve_linear(mc,y,7,ierr)
       if (ierr /= 0) return
    endif
    r = y/norm2(y/sc)
 enddo
 if (have_prev) then
    if (dot_product(r/sc,rprev/sc) < 0._dp) r = -r
 elseif (r(1) > 0._dp) then
    r = -r
 endif

end subroutine eigenvector

!-----------------------------------------------------------------------
! characteristic speed of wave family ifam moving in direction sig
!-----------------------------------------------------------------------
real(dp) function wave_speed(u,ifam,sig)
 real(dp), intent(in) :: u(7),sig
 integer,  intent(in) :: ifam
 real(dp) :: lam(6)
 integer  :: i

 call get_speeds(u,lam)
 select case(ifam)
 case(ifast)
    i = 1
 case(ialfven)
    i = 2
 case default
    i = 3
 end select
 if (sig > 0._dp) i = 7 - i
 wave_speed = lam(i)

end function wave_speed

!-----------------------------------------------------------------------
!
! characteristic speeds (fast-, Alfven-, slow-, slow+, Alfven+, fast+).
! The magnetosonic speeds are the roots of the quartic
!  rho h (1-cs^2) a^4 - (1 - lam^2)[(b^2 + rho h cs^2) a^2 - cs^2 B^2] = 0
! with a = W (lam - vx) and B = b^x - lam b^0. The fast speeds are found
! by Newton iteration from lam = -1 and 1 (which converges monotonically
! to the outermost roots since all roots are real) and the slow speeds
! from the quadratic left after dividing out the fast roots
!
!-----------------------------------------------------------------------
subroutine get_speeds(u,lam)
 real(dp), intent(in)  :: u(7)
 real(dp), intent(out) :: lam(6)
 real(dp) :: w,b0,bvec(3),b2,e,cs2,et,q1(3),q4(5),t2(3),c(5),quad(3),tmp(4)
 real(dp) :: sw,disc,qq,r1,r2
 integer  :: i

 call get_aux(u,w,b0,bvec,b2)
 e   = u(1) + gam/(gam - 1._dp)*u(2)
 cs2 = gam*u(2)/e
 q1  = (/u(3)**2,-2._dp*u(3),1._dp/)
 q4  = poly_mul(q1,q1)
 t2  = (b2 + e*cs2)*w**2*q1 - cs2*(/bvec(1)**2,-2._dp*bvec(1)*b0,b0**2/)
 c   = e*(1._dp - cs2)*w**4*q4
 c(1:5) = c(1:5) - poly_mul((/1._dp,0._dp,-1._dp/),t2)

 lam(6) = poly_root_newton(c, 1._dp)
 lam(1) = poly_root_newton(c,-1._dp)
 !--divide out the fast roots to leave a quadratic for the slow speeds
 call poly_deflate(c,lam(6),tmp)
 call poly_deflate(tmp,lam(1),quad)
 disc = max(quad(2)**2 - 4._dp*quad(3)*quad(1),0._dp)
 qq = -0.5_dp*(quad(2) + sign(sqrt(disc),quad(2)))
 r1 = qq/quad(3)
 r2 = r1
 if (abs(qq) > tiny(qq)) r2 = quad(1)/qq
 lam(3) = poly_root_newton(c,min(r1,r2))
 lam(4) = poly_root_newton(c,max(r1,r2))

 !--Alfven speeds
 et = e + b2
 sw = sqrt(et)*w
 r1 = (bvec(1) + sw*u(3))/(b0 + sw)
 r2 = (bvec(1) - sw*u(3))/(b0 - sw)
 lam(2) = min(r1,r2)
 lam(5) = max(r1,r2)
 do i=1,6
    lam(i) = max(min(lam(i),1._dp),-1._dp)
 enddo

end subroutine get_speeds

!-----------------------------------------------------------------------
! product of polynomials with coefficients in ascending powers
!-----------------------------------------------------------------------
function poly_mul(a,b) result(c)
 real(dp), intent(in) :: a(:),b(:)
 real(dp) :: c(size(a)+size(b)-1)
 integer  :: i,j

 c = 0._dp
 do i=1,size(a)
    do j=1,size(b)
       c(i+j-1) = c(i+j-1) + a(i)*b(j)
    enddo
 enddo

end function poly_mul

!-----------------------------------------------------------------------
! divide polynomial c (ascending powers) by (x - x0), dropping the remainder
!-----------------------------------------------------------------------
subroutine poly_deflate(c,x0,q)
 real(dp), intent(in)  :: c(:),x0
 real(dp), intent(out) :: q(size(c)-1)
 integer  :: n,k

 n = size(c)
 q(n-1) = c(n)
 do k=n-2,1,-1
    q(k) = c(k+1) + x0*q(k+1)
 enddo

end subroutine poly_deflate

!-----------------------------------------------------------------------
! root of the polynomial c (ascending powers) by Newton iteration from x0
!-----------------------------------------------------------------------
real(dp) function poly_root_newton(c,x0) result(x)
 real(dp), intent(in) :: c(:),x0
 real(dp) :: f,df,dx
 integer  :: its,k

 x = x0
 do its=1,100
    f  = c(size(c))
    df = 0._dp
    do k=size(c)-1,1,-1
       df = df*x + f
       f  = f*x + c(k)
    enddo
    if (abs(df) <= tiny(df)) exit
    dx = -f/df
    x  = x + dx
    if (abs(dx) <= 4._dp*epsilon(x)*max(abs(x),1.e-3_dp)) exit
 enddo

end function poly_root_newton

!-----------------------------------------------------------------------
! Jacobians of the conserved variables and fluxes with respect to the
! primitive variables, by centred finite differences
!-----------------------------------------------------------------------
subroutine get_jacobians(u,dudw,dfdw)
 real(dp), intent(in)  :: u(7)
 real(dp), intent(out) :: dudw(7,7),dfdw(7,7)
 real(dp) :: up(7),um(7),ucp(7),ucm(7),fp(7),fm(7),sc(7),h
 integer  :: j

 sc = state_scale(u)
 do j=1,7
    h = 1.e-6_dp*max(abs(u(j)),sc(j))
    up = u
    um = u
    up(j) = u(j) + h
    um(j) = u(j) - h
    call get_cons_flux(up,ucp,fp)
    call get_cons_flux(um,ucm,fm)
    dudw(:,j) = (ucp - ucm)/(2._dp*h)
    dfdw(:,j) = (fp - fm)/(2._dp*h)
 enddo

end subroutine get_jacobians

!-----------------------------------------------------------------------
! typical magnitude of each primitive variable, used for scaling
!-----------------------------------------------------------------------
function state_scale(u) result(sc)
 real(dp), intent(in) :: u(7)
 real(dp) :: sc(7)

 sc(1)   = u(1)
 sc(2)   = ptot(u)
 sc(3:5) = 1._dp
 sc(6:7) = sqrt(ptot(u))

end function state_scale

!-----------------------------------------------------------------------
!
! conserved variables U = (D, E, Sx, Sy, Sz, By, Bz) and their fluxes
! in the x direction, for primitive variables (rho,p,vx,vy,vz,By,Bz)
! (e.g. Giacomazzo & Rezzolla 2006 eq. 2.8)
!
!-----------------------------------------------------------------------
subroutine get_cons_flux(u,ucons,flux)
 real(dp), intent(in)  :: u(7)
 real(dp), intent(out) :: ucons(7),flux(7)
 real(dp) :: w,b0,bvec(3),b2,wtot,pt,s(3)

 call get_aux(u,w,b0,bvec,b2)
 wtot = u(1) + gam/(gam - 1._dp)*u(2) + b2
 pt   = u(2) + 0.5_dp*b2
 s    = wtot*w**2*u(3:5) - b0*bvec
 ucons(1)   = u(1)*w
 ucons(2)   = wtot*w**2 - pt - b0**2
 ucons(3:5) = s
 ucons(6:7) = u(6:7)
 flux(1)    = u(1)*w*u(3)
 flux(2)    = s(1)
 flux(3:5)  = s*u(3) - bvec*bn/w
 flux(3)    = flux(3) + pt
 flux(6:7)  = u(6:7)*u(3) - bn*u(4:5)

end subroutine get_cons_flux

!-----------------------------------------------------------------------
! Lorentz factor, magnetic four-vector b (time and space parts) and b^2
!-----------------------------------------------------------------------
subroutine get_aux(u,w,b0,bvec,b2)
 real(dp), intent(in)  :: u(7)
 real(dp), intent(out) :: w,b0,bvec(3),b2
 real(dp) :: bfield(3),vdotb

 bfield = (/bn,u(6),u(7)/)
 w      = 1._dp/sqrt(1._dp - sum(u(3:5)**2))
 vdotb  = dot_product(u(3:5),bfield)
 b0     = w*vdotb
 bvec   = bfield/w + b0*u(3:5)
 b2     = sum(bfield**2)/w**2 + vdotb**2

end subroutine get_aux

!-----------------------------------------------------------------------
! total (gas + magnetic) pressure
!-----------------------------------------------------------------------
real(dp) function ptot(u)
 real(dp), intent(in) :: u(7)
 real(dp) :: w,b0,bvec(3),b2

 call get_aux(u,w,b0,bvec,b2)
 ptot = u(2) + 0.5_dp*b2

end function ptot

!-----------------------------------------------------------------------
! check that a state has positive density and pressure and |v| < 1
!-----------------------------------------------------------------------
logical function valid_state(u)
 real(dp), intent(in) :: u(7)

 valid_state = (u(1) > 0._dp .and. u(2) > 0._dp .and. sum(u(3:5)**2) < 1._dp)

end function valid_state

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

end module mhdshock_sr
