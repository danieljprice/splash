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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

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
!               University of Exeter 2004-2008
!               Monash University 2008-
!
! daniel.price@monash.edu
!-----------------------------------------------------------------------
module shock_sr
 implicit none
 public :: exact_shock_sr

contains

subroutine exact_shock_sr(iplot,time,gamma,rho_L,rho_R,p_L,p_R,v_L,v_R,xplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot
  integer, intent(out) :: ierr
  real, intent(in) :: time,gamma
  real, intent(in) :: rho_L,rho_R,p_L,p_R,v_L,v_R
  real, dimension(:), intent(inout) :: xplot
  real, dimension(size(xplot)), intent(out) :: yplot
  double precision, dimension(size(xplot)) :: rad,dens,pr,vel,uu
  double precision :: rhol,rhor,pl,prr,vl,vr,gam,t

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
  rad(:) = xplot(:)
  call riemann(size(xplot),rad,dens,pr,vel,uu,rhol,rhor,pl,prr,vl,vr,gam,t,0.0d0)

!------------------------------------
!  determine which solution to plot
!------------------------------------
  select case(iplot)
  case(1)
     yplot = real(dens)
  case(2)
     yplot = real(pr)
  case(3)
     yplot = real(vel)
  case(4)
     if (gamma.gt.1.0001) then
        yplot = real(pr/((gamma-1.)*dens))
     else
        yplot = real(pr/dens)
     endif
  case(5) ! rho*
     yplot = real(dens/dsqrt(1.d0-vel**2))
  end select

  return
end subroutine exact_shock_sr

! ------------
!n    name: r i e m a n n
! ------------

!p    purpose:
!p    this program computes the solution of a 1d
!p    relativistic riemann problem (for constant-gamma ideal gases) with
!p    initial data ul if x<0.5 and ur if x>0.5
!p    in the whole spatial domain [0, 1]
!

!c    comments:
!c    see marti and mueller, jfm, 1994
!c
!c    written by:     jose-maria marti
!c                    departamento de astronomia y astrofisica
!c                    universidad de valencia
!c                    46100 burjassot (valencia), spain
!c                    jose-maria.marti@uv.es
!c    and
!c                    ewald mueller
!c                    max-planck-institut fuer astrophysik
!c                    karl-schwarzschild-str. 1
!c                    85741 garching, germany
!c                    emueller@mpa-garching.mpg.de
!c
!c    Modifications by D. Price (daniel.price@sci.monash.edu.au):
!c    07/2007: used as subroutine instead of program
!c    09/2008: reformatted in free-form and included in exact_shock_sr module
!c

  subroutine riemann(mn,rad,rhoa,pa,vela,ua,rholin,rhorin,plin,prin, &
                     vlin,vrin,gamin,tin,x0)
  implicit none

! -------
! common blocks
! -------

  double precision rhol, pl, ul, hl, csl, vell, wl, &
                   rhor, pr, ur, hr, csr, velr, wr
  common /states/  rhol, pl, ul, hl, csl, vell, wl, &
                   rhor, pr, ur, hr, csr, velr, wr

  double precision rhols, uls, hls, csls, vells, vshockl
  common /ls/      rhols, uls, hls, csls, vells, vshockl

  double precision rhors, urs, hrs, csrs, velrs, vshockr
  common /rs/      rhors, urs, hrs, csrs, velrs, vshockr

  double precision gamma
  common /adind/   gamma

! ---------
! internal variables
! ---------

  integer, intent(in) :: mn
  integer         n, i, iloop
!  parameter      (mn = 400)

  double precision tol, pmin, pmax, dvel1, dvel2, check !, xmin
  double precision, intent(in) :: rholin,rhorin,plin,prin,vlin,vrin,gamin,tin

  double precision ps, vels

  double precision, intent(out) :: rhoa(mn), pa(mn), vela(mn), ua(mn)

  double precision, intent(in) :: x0
  double precision xi

  double precision, intent(in) :: rad(mn)
  double precision x1, x2, x3, x4, x5, t

! -------
! initial states
! -------

!  write(*,*) ' adiabatic index of the gas: '
! read (*,*)   gamma
  gamma = gamin

!  write(*,*) ' time for the solution: '
!  read (*,*)   t
  t = tin
!  xmin = x0 - 0.5d0

! -----
! left state
! -----

!  write(*,*) ' -- left state -- '
!  write(*,*) '      pressure     : '
!  read (*,*) pl
!  write(*,*) '      density      : '
!  read (*,*) rhol
!  write(*,*) '      flow velocity: '
!  read (*,*) vell
  pl = plin
  rhol = rholin
  vell = vlin

  pr = prin
  rhor = rhorin
  velr = vrin

! ------
! right state
! ------

!  write(*,*) ' -- right state -- '
!  write(*,*) '      pressure     : '
!  read (*,*) pr
!  write(*,*) '      density      : '
!  read (*,*) rhor
!  write(*,*) '      flow velocity: '
!  read (*,*) velr

! ------------------------------
! specific internal energy, specific enthalpy, sound speed and
! flow lorentz factors in the initial states
! ------------------------------

  ul  = pl/(gamma-1.d0)/rhol
  ur  = pr/(gamma-1.d0)/rhor

  hl  = 1.d0+ul+pl/rhol
  hr  = 1.d0+ur+pr/rhor

  csl = dsqrt(gamma*pl/rhol/hl)
  csr = dsqrt(gamma*pr/rhor/hr)

  wl  = 1.d0/dsqrt(1.d0-vell**2)
  wr  = 1.d0/dsqrt(1.d0-velr**2)

! --------
! number of points
! --------

  n   = mn

! -------------
! tolerance for the solution
! -------------

  tol = 0.d0

!

  iloop = 0

  pmin  = (pl + pr)/2.d0
  pmax  = pmin

5 iloop = iloop + 1

  pmin  = 0.5d0*max(pmin,0.d0)
  pmax  = 2.d0*pmax

  call getdvel(pmin, dvel1)

  call getdvel(pmax, dvel2)

  check = dvel1*dvel2
  if (check.gt.0.d0) goto 5

! ---------------------------
! pressure and flow velocity in the intermediate states
! ---------------------------

  call getp(pmin, pmax, tol, ps)

  vels = 0.5d0*(vells + velrs)

! ---------------
! solution on the numerical mesh
! ---------------

! -----------
! positions of the waves
! -----------

  if (pl.ge.ps) then

    x1 = x0 + (vell - csl )/(1.d0 - vell*csl )*t
    x2 = x0 + (vels - csls)/(1.d0 - vels*csls)*t

  else

    x1 = x0 + vshockl*t
    x2 = x1

  end if

  x3 = x0 + vels*t

  if (pr.ge.ps) then

    x4 = x0 + (vels + csrs)/(1.d0 + vels*csrs)*t
    x5 = x0 + (velr + csr )/(1.d0 + velr*csr )*t

  else

    x4 = x0 + vshockr*t
    x5 = x4

  end if

! ----------
! solution on the mesh
! ----------

  !do i=1,n

  !  rad(i) = xmin + dfloat(i-1)/dfloat(n-1)

  !enddo

  do i=1,n

    if (rad(i).le.x1) then

      pa(i)   = pl
      rhoa(i) = rhol
      vela(i) = vell
      ua(i)   = ul

    else if (rad(i).le.x2) then

      xi = (rad(i) - x0)/t

      call raref(xi, rhol,  csl,  vell,  'l', &
                     rhoa(i), pa(i), ua(i),       vela(i))

    else if (rad(i).le.x3) then

      pa(i)   = ps
      rhoa(i) = rhols
      vela(i) = vels
      ua(i)   = uls

    else if (rad(i).le.x4) then

      pa(i)   = ps
      rhoa(i) = rhors
      vela(i) = vels
      ua(i)   = urs

    else if (rad(i).le.x5) then

      xi = (rad(i) - x0)/t

      call raref(xi, rhor, csr,  velr,  'r', &
                     rhoa(i), pa(i), ua(i),       vela(i))

    else

      pa(i)   = pr
      rhoa(i) = rhor
      vela(i) = velr
      ua(i)   = ur

    end if

   enddo

!  open (3,file='solution.dat',form='formatted',status='new')

!  write(3,150) n, t
! 150  format(i5,1x,f10.5)

!  do 60 i=1,n
!    write(3,200) rad(i),pa(i),rhoa(i),vela(i),ua(i)
! 60   continue

! 200  format(5(e15.8,1x))

!  close(3)

  return
  end subroutine riemann

! ----------
!n    name: g e t d v e l
! ----------

!p    purpose:
!p    compute the difference in flow speed between left and right intermediate
!p    states for given left and right states and pressure
!

!c    comments
!c    none

  subroutine getdvel( p, dvel )

  implicit none

! -----
! arguments
! -----

  doubleprecision, intent(in)  :: p
  doubleprecision, intent(out) :: dvel

! -------
! common blocks
! -------

  double precision rhols,uls,hls,csls,vells,vshockl
  common /ls/      rhols,uls,hls,csls,vells,vshockl

  double precision rhors,urs,hrs,csrs,velrs,vshockr
  common /rs/      rhors,urs,hrs,csrs,velrs,vshockr

  double precision rhol, pl, ul, hl, csl, vell, wl, &
                   rhor, pr, ur, hr, csr, velr, wr
  common /states/  rhol, pl, ul, hl, csl, vell, wl, &
                   rhor, pr, ur, hr, csr, velr, wr

  double precision gamma
  common /adind/   gamma

! -----
! left wave
! -----

  call getvel(p, rhol, pl, hl,  csl,  vell,  wl, 'l', &
                 rhols,    uls, hls, csls, vells, vshockl )

! -----
! right wave
! -----

  call getvel(p, rhor, pr, hr,  csr,  velr,  wr, 'r', &
                 rhors,    urs, hrs, csrs, velrs, vshockr )

  dvel = vells - velrs

  return
  end subroutine getdvel

! -------
!n    name: g e t p
! -------

!p    purpose:
!p    find the pressure in the intermediate state of a riemann problem in
!p    relativistic hydrodynamics
!

!c    comments:
!c    this routine uses a combination of interval bisection and inverse
!c    quadratic interpolation to find the root in a specified interval.
!c    it is assumed that dvel(pmin) and dvel(pmax) have opposite signs without
!c    a check.
!c    adapted from "computer methods for mathematical computation",
!c    by g. e. forsythe, m. a. malcolm, and c. b. moler,
!c    prentice-hall, englewood cliffs n.j.
!
  subroutine getp( pmin, pmax, tol, ps )

  implicit none

! -----
! arguments
! -----

  doubleprecision, intent(in) :: pmin, pmax, tol
  doubleprecision, intent(out) :: ps

! -------
! common blocks
! -------

  doubleprecision gamma
  common /adind/  gamma

  doubleprecision rhol, pl, ul, hl, csl, vell, wl, &
                  rhor, pr, ur, hr, csr, velr, wr
  common /states/ rhol, pl, ul, hl, csl, vell, wl, &
                  rhor, pr, ur, hr, csr, velr, wr

! ---------
! internal variables
! ---------

  doubleprecision a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, q, r, s

! -------------
! compute machine precision
! -------------

  eps  = 1.d0
10    eps  = eps/2.d0
  tol1 = 1.d0 + eps
  if( tol1 .gt. 1.d0 ) go to 10

! -------
! initialization
! -------

  a  = pmin
  b  = pmax
  call getdvel(a,fa)
  call getdvel(b,fb)

! -----
! begin step
! -----

20    c  = a
  fc = fa
  d  = b - a
  e  = d
30    if( dabs(fc) .ge. dabs(fb) )go to 40
  a  = b
  b  = c
  c  = a
  fa = fb
  fb = fc
  fc = fa

! --------
! convergence test
! --------

40    tol1 = 2.d0*eps*dabs(b) + 0.5d0*tol
  xm   = 0.5d0*(c - b)
  if( dabs(xm) .le. tol1 ) go to 90
  if( fb .eq. 0.d0 ) go to 90

! ------------
! is bisection necessary?
! ------------

  if( dabs(e) .lt. tol1 ) go to 70
  if( dabs(fa) .le. dabs(fb) ) go to 70

! ------------------
! is quadratic interpolation possible?
! ------------------

  if( a .ne. c ) go to 50

! ----------
! linear interpolation
! ----------

  s = fb/fa
  p = 2.d0*xm*s
  q = 1.d0 - s
  go to 60

! ----------------
! inverse quadratic interpolation
! ----------------

50 q = fa/fc
  r = fb/fc
  s = fb/fa
  p = s*(2.d0*xm*q*(q - r) - (b - a)*(r - 1.d0))
  q = (q - 1.d0)*(r - 1.d0)*(s - 1.d0)

! ------
! adjust signs
! ------

60 if( p .gt. 0.d0 ) q = -q
  p = dabs(p)

! --------------
! is interpolation acceptable?
! --------------

  if( (2.d0*p) .ge. (3.d0*xm*q-dabs(tol1*q)) ) go to 70
  if( p .ge. dabs(0.5d0*e*q) ) go to 70
  e = d
  d = p/q
  go to 80

! -----
! bisection
! -----

70 d = xm
  e = d

! -------
! complete step
! -------

80 a  = b
  fa = fb
  if( dabs(d) .gt. tol1 ) b = b+d
  if( dabs(d) .le. tol1 ) b = b+dsign(tol1,xm)
  call getdvel(b,fb)
  if( (fb*(fc/dabs(fc))) .gt. 0.d0) go to 20
  go to 30

! --
! done
! --

90 ps = b

  return
  end subroutine getp

! ---------
!n    name: g e t v e l
! ---------

!p    purpose:
!p    compute the flow velocity behind a rarefaction or shock in terms of the
!p    post-wave pressure for a given state ahead the wave in a relativistic
!p    flow
!

!c    comments:
!c    this routine closely follows the expressions in Marti and Mueller,
!c    J. fluid mech., (1994)

  subroutine getvel( p, rhoa, pa, ha, csa, vela, wa, s,  &
                     rho, u,  h,  cs,  vel,  vshock )

  implicit none

! -----
! arguments
! -----

  double precision, intent(in) :: p, rhoa, pa, ha, csa, vela, wa
  character(len=1), intent(in) :: s
  double precision, intent(out) :: rho, u, h, cs, vel, vshock

! -------
! common blocks
! -------

  double precision gamma
  common /adind/   gamma

! ---------
! internal variables
! ---------

  double precision a, b, c, sign
  double precision j, wshock
  double precision k, sqgl1

! ---------------
! left or right propagating wave
! ---------------

  sign = 0.d0
  if (s.eq.'l') sign = -1.d0

  if (s.eq.'r') sign =  1.d0

!

  if (p.gt.pa) then

!   ---
!   shock
!   ---

    a  = 1.d0+(gamma-1.d0)*(pa-p)/gamma/p
    b  = 1.d0-a
    c  = ha*(pa-p)/rhoa-ha**2

!   ----------------
!   check for unphysical enthalpies
!   ----------------

    if (c.gt.(b**2/4.d0/a)) then
       print*,'getvel: unphysical specific enthalpy in intermediate state'
       return
    endif

!   -----------------------------
!   specific enthalpy in the post-wave state
!   (from the equation of state and the taub adiabat,
!   eq.(74), mm94)
!   -----------------------------

    h = (-b+dsqrt(b**2-4.d0*a*c))/2.d0/a

!   ---------------
!   density in the post-wave state
!   (from eq.(73), mm94)
!   ---------------

    rho = gamma*p/(gamma-1.d0)/(h-1.d0)

!   ------------------------
!   specific internal energy in the post-wave state
!   (from the equation of state)
!   ------------------------

    u = p/(gamma-1.d0)/rho

!   --------------------------
!   mass flux across the wave
!   (from the rankine-hugoniot relations, eq.(71), mm94)
!   --------------------------

    j = sign*dsqrt((p-pa)/(ha/rhoa-h/rho))

!   ----------
!   shock velocity
!   (from eq.(86), mm94
!   ----------

    a      = j**2+(rhoa*wa)**2
    b      = -vela*rhoa**2*wa**2
    vshock = (-b+sign*j**2*dsqrt(1.d0+rhoa**2/j**2))/a
    wshock = 1.d0/dsqrt(1.d0-vshock**2)

!   -------------------
!   flow velocity in the post-shock state
!   (from eq.(67), mm94)
!   -------------------

    a = wshock*(p-pa)/j+ha*wa*vela
    b = ha*wa+(p-pa)*(wshock*vela/j+1.d0/rhoa/wa)

    vel = a/b

!   ---------------------
!   local sound speed in the post-shock state
!   (from the equation of state)
!   ---------------------

    cs = dsqrt(gamma*p/rho/h)

  else

!   ------
!   rarefaction
!   ------

!   ---------------------------
!   politropic constant of the gas across the rarefaction
!   ---------------------------

    k = pa/rhoa**gamma

!   ---------------
!   density behind the rarefaction
!   ---------------

    rho = (p/k)**(1.d0/gamma)

!   ------------------------
!   specific internal energy behind the rarefaction
!   (from the equation of state)
!   ------------------------

    u = p/(gamma-1.d0)/rho

!   --------------------
!   local sound speed behind the rarefaction
!   (from the equation of state)
!   --------------------

    cs = dsqrt(gamma*p/(rho+gamma*p/(gamma-1.d0)))

!   ------------------
!   flow velocity behind the rarefaction
!   ------------------

    sqgl1 = dsqrt(gamma-1.d0)
    a   = (1.d0+vela)/(1.d0-vela)*((sqgl1+csa)/(sqgl1-csa)* &
          (sqgl1-cs )/(sqgl1+cs ))**(-sign*2.d0/sqgl1)

    vel = (a-1.d0)/(a+1.d0)

  end if

  end subroutine getvel

! --------
!n    name: r a r e f
! --------

!p    purpose:
!p    compute the flow state in a rarefaction for given pre-wave state
!

!c    comments:
!c    this routine closely follows the expressions in marti and mueller,
!c    j. fluid mech., (1994)

  subroutine raref( xi, rhoa, csa, vela, s, rho, p, u, vel )

  implicit none

! -----
! arguments
! -----

  double precision, intent(in) :: xi

  double precision, intent(in) :: rhoa, csa, vela

  character, intent(in) ::  s

  double precision, intent(out) :: rho, p, u, vel

! -------
! common blocks
! -------

  double precision gamma
  common /adind/   gamma

! ---------
! internal variables
! ---------

  double precision b, c, d, k, l, v, ocs2, fcs2, dfdcs2, cs2, sign

! ---------------
! left or right propagating wave
! ---------------

  sign = 0.d0
  if (s.eq.'l') sign =  1.d0

  if (s.eq.'r') sign = -1.d0

  b    = dsqrt(gamma - 1.d0)
  c    = (b + csa)/(b - csa)
  d    = -sign*b/2.d0
  k    = (1.d0 + xi)/(1.d0 - xi)
  l    = c*k**d
  v    = ((1.d0 - vela)/(1.d0 + vela))**d

  ocs2 = csa

25 fcs2   = l*v*(1.d0 + sign*ocs2)**d*(ocs2 - b) + (1.d0 - sign*ocs2)**d*(ocs2 + b)

  dfdcs2 = l*v*(1.d0 + sign*ocs2)**d* &
   (1.d0 + sign*d*(ocs2 - b)/(1.d0 + sign*ocs2)) + &
   (1.d0 - sign*ocs2)**d* &
   (1.d0 - sign*d*(ocs2 + b)/(1.d0 - sign*ocs2))

  cs2 = ocs2 - fcs2/dfdcs2

  if (abs(cs2 - ocs2)/ocs2.gt.5.e-7)then
    ocs2 = cs2
    goto 25
  end if

  vel = (xi + sign*cs2)/(1.d0 + sign*xi*cs2)

  rho = rhoa*((cs2**2*(gamma - 1.d0 - csa**2))/ &
   (csa**2*(gamma - 1.d0 - cs2**2)))**(1.d0/(gamma - 1.d0))

  p   = cs2**2*(gamma - 1.d0)*rho/(gamma - 1.d0 - cs2**2)/gamma

  u   = p/(gamma - 1.d0)/rho

  return
  end subroutine raref

end module shock_sr
