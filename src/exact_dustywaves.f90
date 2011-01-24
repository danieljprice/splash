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
!  Copyright (C) 2005-2011 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

! ----------------------------------------------------------------------
! compute exact solution for a linear wave
! plots a sine function with a given amplitude, period and wavelength
! ----------------------------------------------------------------------
module dustywaves
  implicit none

contains

subroutine exact_dustywave(iplot,time,ampl,cs,Kdragin,lambda,x0,ymean,xplot,yplot,ierr)
  use cubic,   only:cubicsolve_complex
  use plotlib, only:plot_line,plot_sls
  implicit none
  integer :: i
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: iplot
  real, intent(in)    :: time, ampl, cs, Kdragin, lambda, x0, ymean 
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  real, dimension(size(xplot)) :: yplot1
  real :: rhodeq,rhogeq,rhodsol,rhogsol,vdeq,vgeq,vgsol,vdsol
  real :: aa,bb,cc,w1r,w2r,w3r,w1i,w2i,w3i
  real :: k,xk,arg1,arg2,arg3,vgas,vdust,rhogas,rhodust
  real :: vd1r,vd1i,vd2r,vd2i,vd3r,vd3i
  real :: vg1r,vg1i,vg2r,vg2i,vg3r,vg3i
  real :: rhod1r,rhod1i,rhod2r,rhod2i,rhod3r,rhod3i
  real :: rhog1r,rhog1i,rhog2r,rhog2i,rhog3r,rhog3i
  real :: tgas1,tdust1,Kdrag
  complex :: xc(3)

  Kdrag = Kdragin

  print*,'plotting two-fluid gas/dust linear wave solution ... '
  print*,' lambda = ',lambda,' ampl = ',ampl,' cs = ',cs,' Kdrag = ',Kdrag
!
! check for errors
!
  ierr = 0
  if (ampl.lt.0.) then
     print*,'error: amplitude < 0 on input'
     ierr = 1
     return
  endif
  if (lambda <= 0.) then
     print*,'error: lambda <= 0 on input'
     ierr = 2
     return
  endif
  if (cs <= 0) then
     print*,'error: sound speed <= 0 on input'
     ierr = 3
     return
  endif
  if (ymean <= 0) then
     print*,'error: mean density <= 0 on input'
     ierr = 4
     return
  endif
  if (Kdrag < 0) then
     print*,'error: drag coefficient < 0 on input'
     ierr = 5
     return
  elseif (abs(Kdrag).lt.1.e-8) then
     print*,' WARNING: Kdrag = 0 on input; using tiny to avoid divergence '
     Kdrag = 0.
  endif

  rhodeq  = ymean ! initial dust density
  rhogeq  = ymean ! initial gas density
  print*,' rho(dust),0 = ',ymean,' rho(gas),0 = ',ymean
  rhodsol = ampl  ! amplitude of dust density perturbation
  rhogsol = ampl  ! amplitude of gas density perturbation
  vdeq    = 0.
  vgeq    = 0.
  vgsol   = ampl  ! amplitude of gas velocity perturbation
  vdsol   = ampl  ! amplitude of dust velocity perturbation
  !cs      = 1.0
  !Kdrag   = 1.0
  k       = 2.*pi/lambda ! wavenumber
  
  vd1r = 0.
  vd1i = 0.
  vd2r = 0.
  vd2i = 0.
  vd3r = 0.
  vd3i = 0.
  rhod1r = 0.
  rhod1i = 0.
  rhod2r = 0.
  rhod2i = 0.
  rhod3r = 0.
  rhod3i = 0.
  !
  !--solve cubic to get the 3 solutions for omega
  !  (these each have both real and imaginary components,
  !   labelled w1r, w1i etc.)
  !
  tdust1 = Kdrag/rhodeq
  tgas1  = Kdrag/rhogeq
  aa = -(tdust1 + tgas1)
  bb = k**2*cs**2
  cc = bb*tdust1
  
  call cubicsolve_complex(aa,bb,cc,xc)
  !--get solutions for (w = -iy instead of y)
  xc = xc*cmplx(0.,1.)
  print*,' roots are ',xc
  
  w1r = real(xc(1))
  w2r = real(xc(2))
  w3r = real(xc(3))
  
  w1i = aimag(xc(1))
  w2i = aimag(xc(2))
  w3i = aimag(xc(3))

!-------------------------------
! G A S  V E L O C I T I E S
!-------------------------------
  vg3r = - ( - w3r**2*w1i*rhogsol*w2i**2 + w1i*w3r**3*w2r*rhogsol + w1i*w2i*w3i*w3r**2*rhogsol +&
        rhogeq*w3i*w3r**2*w2r*k*vgsol + w2i*w3r**3*w1r*rhogsol - w3r*w3i**2*k*Kdrag*vgsol - w3r*w2r*w1r*k*Kdrag*vgsol +&
        w3r*w1i*w3i**2*w2r*rhogsol - w3r*w2i*rhogeq*w3i**2*k*vgsol + w3r*w2i*rhogeq*w1r**2*k*vgsol +&
        w3r*w2i*w3i**2*w1r*rhogsol + w3r*w2r*w1r*k*Kdrag*vdsol + rhogeq*w3i*w3r**2*k*w1r*vgsol -&
        w2i*rhogeq*w3r**3*k*vgsol - w3r**3*k*Kdrag*vgsol - w3i*w3r**2*k**2*cs**2*rhogsol + w3r**3*k*Kdrag*vdsol -&
        w3i*w3r**2*w2r*w1r*rhogsol + w3r*w1i*w2i**2*rhogeq*k*vgsol - w3i*w2i**2*rhogeq*k*w1r*vgsol +&
        w3r*w3i**2*k*Kdrag*vdsol - w3r*w1i*w2r*k**2*cs**2*rhogsol + w3r*w1i*w2i*k*Kdrag*vgsol +&
        w3i*w2i*w1r*k*Kdrag*vdsol - w3i**2*w2r*k*Kdrag*vdsol + w3i**3*rhogeq*w2r*k*vgsol -&
        w3r*w1i*rhogeq*w3i**2*k*vgsol - w3i*w2i*Kdrag*w1r*k*vgsol - w3r*w1i*w2i*k*Kdrag*vdsol +&
        w3r*w1i*w2r**2*rhogeq*k*vgsol - w3r*w2i*k**2*cs**2*w1r*rhogsol + w3r*w1i**2*w2i*rhogeq*k*vgsol -&
        w1i*rhogeq*w3r**3*k*vgsol - w3i*rhogeq*w2r*w1r**2*k*vgsol + w3i**2*w2r*k*Kdrag*vgsol - w3i**2*k*Kdrag*vdsol*w1r&
        + w3i**2*k*Kdrag*vgsol*w1r + w3i**2*w2i*k**2*cs**2*rhogsol - w3i*w1i**2*rhogeq*w2r*k*vgsol -&
        w3r**2*w2i*w1r**2*rhogsol + w3i**2*w1i*k**2*cs**2*rhogsol + w3i*w1i*w2r*k*Kdrag*vdsol +&
        w3r**2*w2i*k**2*cs**2*rhogsol - w3i*w1i*w2r*k*Kdrag*vgsol - w3i*w1i*w2i*k**2*cs**2*rhogsol -&
        w3i*rhogeq*w2r**2*k*w1r*vgsol + w3i**3*rhogeq*k*w1r*vgsol + w3i*w2r*w1r*k**2*cs**2*rhogsol +&
        w3r**2*w2r*k*Kdrag*vgsol - w3r**2*w2r*k*Kdrag*vdsol - w3r**2*k*Kdrag*vdsol*w1r + w3r**2*k*Kdrag*vgsol*w1r -&
        w3r**2*w1i*w2r**2*rhogsol + w3r**2*w1i*k**2*cs**2*rhogsol - w3r**2*w1i**2*rhogsol*w2i -&
        w3i**2*w1i*w2r**2*rhogsol - w3i**2*w1i*rhogsol*w2i**2 - w3i**3*k**2*cs**2*rhogsol + w3i*w2i**2*w1r**2*rhogsol +&
        w3i*w2r**2*w1r**2*rhogsol - w3i**2*w2i*w1r**2*rhogsol + w3i*w1i**2*rhogsol*w2i**2 + w3i**3*w1i*w2i*rhogsol -&
        w3i**2*w1i**2*rhogsol*w2i + w3i*w1i**2*w2r**2*rhogsol - w3i**3*w2r*w1r*rhogsol)/rhogeq/k/( - 2*w2r*w3r + w2r**2&
        + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  vg3i =(w3i*w2r*w1r*k*Kdrag*vgsol - w2r*w3r**2*k**2*cs**2*rhogsol - w3i**2*w2r*k**2*cs**2*rhogsol +&
        w2r*w1r*w3r*k**2*cs**2*rhogsol + w3r**2*k*w3i*Kdrag*vdsol + w3i**2*w3r*w2r*w1r*rhogsol + w3i**3*k*Kdrag*vdsol +&
        w3r**3*w2r*w1r*rhogsol - w1r*w3r**2*w2r**2*rhogsol + w3r**3*k**2*cs**2*rhogsol + w3r*w2r**2*w1r**2*rhogsol -&
        w3i**2*w1r*w2r**2*rhogsol - w3i**3*k*Kdrag*vgsol - w2i**2*w3i**2*w1r*rhogsol + w2i**2*w3r*w1r**2*rhogsol -&
        w2i**2*w1r*w3r**2*rhogsol + w1i*w3i**3*w2r*rhogsol + w1i**2*w2i**2*w3r*rhogsol - w1i**2*w3i**2*w2r*rhogsol +&
        w1i**2*w2r**2*w3r*rhogsol - w1i**2*w3r**2*w2r*rhogsol - w3r**2*k*w3i*Kdrag*vgsol +&
        w3i**2*w3r*k**2*cs**2*rhogsol - w3i*w2r*w1r*k*Kdrag*vdsol - w3i**2*w1r*k**2*cs**2*rhogsol -&
        w1r*w3r**2*k**2*cs**2*rhogsol - w1r**2*w3r**2*w2r*rhogsol - w3i**2*w1r**2*w2r*rhogsol +&
        w1i*w3i**2*k*Kdrag*vgsol - w1i**2*w3i*w2i*rhogeq*k*vgsol + w1i**2*rhogeq*w3i**2*k*vgsol -&
        w1i**2*rhogeq*w2r*w3r*k*vgsol + w1i**2*rhogeq*w3r**2*k*vgsol - w1i*w2i*w3r**3*rhogsol -&
        w2i*rhogeq*w3i**3*k*vgsol - w2i*rhogeq*w3i*w1r**2*k*vgsol - w2i*rhogeq*w3r**2*k*w3i*vgsol +&
        w2i*w3i*w3r**2*w1r*rhogsol - w2i*w3r**2*k*Kdrag*vdsol + w2i*w3i**3*w1r*rhogsol + w2i*w1r*w3r*k*Kdrag*vdsol +&
        w2i*w3i**2*k*Kdrag*vgsol + w2i*k**2*cs**2*w3i*w1r*rhogsol - w2i*w3r*Kdrag*w1r*k*vgsol -&
        w2i*w3i**2*k*Kdrag*vdsol + w2i**2*rhogeq*w3i**2*k*vgsol - w2i**2*rhogeq*w3r*k*w1r*vgsol +&
        w2i**2*rhogeq*w3r**2*k*vgsol + w2i*w3r**2*k*Kdrag*vgsol + rhogeq*w2r**2*w3i**2*k*vgsol -&
        rhogeq*w3r*w2r**2*k*w1r*vgsol - rhogeq*w3i**2*w3r*w2r*k*vgsol - rhogeq*w3r**3*k*w1r*vgsol -&
        rhogeq*w3r*w2r*w1r**2*k*vgsol - rhogeq*w3r**3*w2r*k*vgsol + rhogeq*w3i**2*w1r**2*k*vgsol +&
        2*rhogeq*w3i**2*k*vgsol*w1r*w2r + 2*rhogeq*w3r**2*k*vgsol*w1r*w2r - rhogeq*w3i**2*w3r*k*w1r*vgsol +&
        rhogeq*w2r**2*w3r**2*k*vgsol + rhogeq*w1r**2*w3r**2*k*vgsol - w1i*w2r**2*w3i*rhogeq*k*vgsol -&
        w1i*w2i**2*w3i*rhogeq*k*vgsol + w1i*w3r**2*k*Kdrag*vgsol + w1i*w3r**2*w2r*w3i*rhogsol -&
        w1i*w3r**2*k*Kdrag*vdsol - w1i*rhogeq*w3i**3*k*vgsol - w1i*w2i*w3i**2*w3r*rhogsol +&
        2*w1i*w2i*rhogeq*w3i**2*k*vgsol + 2*w1i*w2i*rhogeq*w3r**2*k*vgsol + w1i*w2i*k*w3i*Kdrag*vdsol -&
        w1i*w2i*w3r*k**2*cs**2*rhogsol + w1i*w3i*w2r*k**2*cs**2*rhogsol + w1i*w2r*w3r*k*Kdrag*vdsol -&
        w1i*w2r*w3r*k*Kdrag*vgsol - w1i*rhogeq*w3r**2*k*w3i*vgsol - w1i*w2i*k*w3i*Kdrag*vgsol -&
        w1i*w3i**2*k*Kdrag*vdsol)/rhogeq/k/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i&
        + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  vg2r =(w3r**2*w1i*rhogsol*w2i**2 + w2r**3*k*Kdrag*vgsol + w3r*w2r*w1r*k*Kdrag*vgsol + w3r*w2i*rhogeq*w1r**2*k*vgsol -&
        w3r*w2r*w1r*k*Kdrag*vdsol + w3r*w1i*w2r*k**2*cs**2*rhogsol + w3r*w1i*w2i*k*Kdrag*vgsol -&
        w3i*w2i*w1r*k*Kdrag*vdsol + w3i*w2i*Kdrag*w1r*k*vgsol - w3r*w1i*w2i*k*Kdrag*vdsol -&
        w3r*w2i*k**2*cs**2*w1r*rhogsol + w3r*w1i**2*w2i*rhogeq*k*vgsol - w3i*rhogeq*w2r*w1r**2*k*vgsol -&
        w3i*w1i**2*rhogeq*w2r*k*vgsol - w3r**2*w2i*w1r**2*rhogsol + w3i*w1i*w2r*k*Kdrag*vdsol -&
        w3i*w1i*w2r*k*Kdrag*vgsol + w3i*w1i*w2i*k**2*cs**2*rhogsol + w3i*w2r*w1r*k**2*cs**2*rhogsol +&
        w3r**2*w1i*w2r**2*rhogsol - w3r**2*w1i**2*rhogsol*w2i + w3i**2*w1i*w2r**2*rhogsol + w3i**2*w1i*rhogsol*w2i**2 +&
        w3i*w2i**2*w1r**2*rhogsol + w3i*w2r**2*w1r**2*rhogsol - w3i**2*w2i*w1r**2*rhogsol + w3i*w1i**2*rhogsol*w2i**2 -&
        w3i**2*w1i**2*rhogsol*w2i + w3i*w1i**2*w2r**2*rhogsol - w2i**3*w1i*w3i*rhogsol + w2i**3*k**2*cs**2*rhogsol +&
        w3i**2*w2i*rhogeq*k*w1r*vgsol + w2i*rhogeq*w3r**2*k*w1r*vgsol - w2i**2*Kdrag*w1r*k*vgsol -&
        w2i**2*w3i*k**2*cs**2*rhogsol - w2i**3*rhogeq*k*w1r*vgsol - w2i**2*w3r*k*Kdrag*vgsol + w2i**2*w1r*k*Kdrag*vdsol&
        + w2i**2*w3r*k*Kdrag*vdsol + w2i**3*w3r*w1r*rhogsol - w2i**2*w1i*k**2*cs**2*rhogsol - w2i**3*w3r*rhogeq*k*vgsol&
        - w2r*w2i**2*k*Kdrag*vdsol + w2r*w2i**2*k*Kdrag*vgsol - w1i*w2r**2*w3i*rhogsol*w2i +&
        w2r*w2i**2*rhogeq*w3i*k*vgsol - w2r*w1i*rhogeq*w3i**2*k*vgsol - w2r*w1i*rhogeq*w3r**2*k*vgsol +&
        w2r*w1i*w2i**2*rhogeq*k*vgsol - w2r*w3i*w2i**2*w1r*rhogsol - w2r*w1i*w2i**2*w3r*rhogsol +&
        w1i*w2r**3*rhogeq*k*vgsol - w3i*w1r*rhogsol*w2r**3 + w2r**2*w1r*w2i*w3r*rhogsol - w2i*w1r*rhogeq*w2r**2*k*vgsol&
        + w3i*w2r**3*rhogeq*k*vgsol - w1i*w2r**3*w3r*rhogsol - w2r**2*w1i*k**2*cs**2*rhogsol + w2r**2*w1r*k*Kdrag*vdsol&
        + w2r**2*w3r*k*Kdrag*vdsol - w2r**2*w3i*k**2*cs**2*rhogsol - w2r**2*w3r*k*Kdrag*vgsol -&
        w2r**2*Kdrag*w1r*k*vgsol - w2r**3*k*Kdrag*vdsol + w2i*w2r**2*k**2*cs**2*rhogsol -&
        w2i*rhogeq*w2r**2*w3r*k*vgsol)/rhogeq/k/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/(w2i**2&
        + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  vg2i = - ( - w2i**3*k*Kdrag*vdsol - w2r**3*k**2*cs**2*rhogsol + w2r**2*w2i*k*Kdrag*vgsol + w2i**2*w3i*k*Kdrag*vdsol -&
        w2i**2*w2r*k**2*cs**2*rhogsol + w2i**3*rhogeq*w3i*k*vgsol - w2i**2*rhogeq*w1r**2*k*vgsol +&
        w2r**2*w3r*k**2*cs**2*rhogsol + w2i**2*rhogeq*w2r*w3r*k*vgsol + w2i**2*k**2*cs**2*w1r*rhogsol -&
        w2r**2*w2i*k*Kdrag*vdsol + w2r**3*rhogeq*w3r*k*vgsol + w3i*w2r**2*w2i*rhogeq*k*vgsol + k*w3i*Kdrag*vdsol*w2r**2&
        - w3i*w2i*w1r*rhogsol*w2r**2 - k*w3i*Kdrag*vgsol*w2r**2 - w2r*w1r*w2i**2*w3r*rhogsol - w2i**2*w3i*k*Kdrag*vgsol&
        + w2i**2*w1r*rhogeq*w2r*k*vgsol + w2i**2*w3r*k**2*cs**2*rhogsol + w2i**3*k*Kdrag*vgsol +&
        w2r**3*rhogeq*k*w1r*vgsol + w2r**2*w1r*k**2*cs**2*rhogsol - w2r**2*rhogeq*w1r**2*k*vgsol -&
        w1i**2*w2i**2*rhogeq*k*vgsol - w1i**2*w2r**2*rhogeq*k*vgsol - w1i*w2i*w2r**2*w3r*rhogsol +&
        w1i*w2i**2*k*Kdrag*vdsol + w1i*w2r**2*w2i*rhogeq*k*vgsol - w1i*w2i**2*k*Kdrag*vgsol + w1i*w2i**3*rhogeq*k*vgsol&
        + w1i*w2r*w3i*rhogsol*w2i**2 + w1i*w2r**2*k*Kdrag*vdsol - w1i*w2r**2*k*Kdrag*vgsol + w1i*w2r**3*w3i*rhogsol -&
        w1i*w2i**3*w3r*rhogsol + w3i*w2r*w1r*k*Kdrag*vgsol - w2r*w1r*w3r*k**2*cs**2*rhogsol + w1r*w3r**2*w2r**2*rhogsol&
        + w3r*w2r**2*w1r**2*rhogsol + w3i**2*w1r*w2r**2*rhogsol + w2i**2*w3i**2*w1r*rhogsol + w2i**2*w3r*w1r**2*rhogsol&
        + w2i**2*w1r*w3r**2*rhogsol + w1i**2*w2i**2*w3r*rhogsol - w1i**2*w3i**2*w2r*rhogsol + w1i**2*w2r**2*w3r*rhogsol&
        - w1i**2*w3r**2*w2r*rhogsol - w3i*w2r*w1r*k*Kdrag*vdsol - w1r**2*w3r**2*w2r*rhogsol - w3i**2*w1r**2*w2r*rhogsol&
        + w1i**2*w3i*w2i*rhogeq*k*vgsol + w1i**2*rhogeq*w2r*w3r*k*vgsol + w2i*rhogeq*w3i*w1r**2*k*vgsol +&
        w2i*w1r*w3r*k*Kdrag*vdsol - w2i*k**2*cs**2*w3i*w1r*rhogsol - w2i*w3r*Kdrag*w1r*k*vgsol -&
        w2i**2*rhogeq*w3i**2*k*vgsol - 2*w2i**2*rhogeq*w3r*k*w1r*vgsol - w2i**2*rhogeq*w3r**2*k*vgsol -&
        rhogeq*w2r**2*w3i**2*k*vgsol - 2*rhogeq*w3r*w2r**2*k*w1r*vgsol + rhogeq*w3r*w2r*w1r**2*k*vgsol +&
        rhogeq*w3i**2*k*vgsol*w1r*w2r + rhogeq*w3r**2*k*vgsol*w1r*w2r - rhogeq*w2r**2*w3r**2*k*vgsol -&
        2*w1i*w2r**2*w3i*rhogeq*k*vgsol - 2*w1i*w2i**2*w3i*rhogeq*k*vgsol + w1i*w2i*rhogeq*w3i**2*k*vgsol +&
        w1i*w2i*rhogeq*w3r**2*k*vgsol - w1i*w2i*k*w3i*Kdrag*vdsol - w1i*w2i*w3r*k**2*cs**2*rhogsol +&
        w1i*w3i*w2r*k**2*cs**2*rhogsol - w1i*w2r*w3r*k*Kdrag*vdsol + w1i*w2r*w3r*k*Kdrag*vgsol +&
        w1i*w2i*k*w3i*Kdrag*vgsol - w2r**3*w1r*w3r*rhogsol - w3i*w2i**3*w1r*rhogsol)/rhogeq/k/( - 2*w2r*w3r + w2r**2 +&
        w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  vg1r = - ( - rhogeq*k*w1r**3*vgsol*w3i - w1r**2*k*Kdrag*vdsol*w3r + w2r*k*Kdrag*vgsol*w1i**2 +&
        w3r*k*Kdrag*vgsol*w1i**2 + w1r*k*Kdrag*vdsol*w1i**2 - w3r*k*Kdrag*vdsol*w1i**2 - Kdrag*w1r*k*vgsol*w1i**2 -&
        Kdrag*w1r**3*k*vgsol - rhogeq*k*w1r*vgsol*w1i**2*w3i + w2r*w3i*w1r*rhogsol*w1i**2 - w2r*w1i*w3r*rhogsol*w1r**2&
        + w2r*k*Kdrag*vgsol*w1r**2 + w2r*w1i*rhogeq*k*vgsol*w1r**2 + w2r*w1i**3*rhogeq*k*vgsol + w2r*w3i*w1r**3*rhogsol&
        - w1i*k**2*cs**2*rhogsol*w1r**2 - w1i**3*k**2*cs**2*rhogsol - w2r*w1i**3*w3r*rhogsol - w2r*k*Kdrag*vdsol*w1r**2&
        + w3r**2*w1i*rhogsol*w2i**2 - w2r*k*Kdrag*vdsol*w1i**2 + w3r*w1i*rhogeq*k*vgsol*w1r**2 +&
        w3r*w1i**3*rhogeq*k*vgsol + w1r**3*k*Kdrag*vdsol - w3r*w2r*w1r*k*Kdrag*vgsol + w3r*w2r*w1r*k*Kdrag*vdsol -&
        w3r*w1i*w2i**2*rhogeq*k*vgsol + w3i*w2i**2*rhogeq*k*w1r*vgsol + w3r*w1i*w2r*k**2*cs**2*rhogsol -&
        w3r*w1i*w2i*k*Kdrag*vgsol - w3i*w2i*w1r*k*Kdrag*vdsol + w3i*w2i*Kdrag*w1r*k*vgsol + w3r*w1i*w2i*k*Kdrag*vdsol -&
        w3r*w1i*w2r**2*rhogeq*k*vgsol - w3r*w2i*k**2*cs**2*w1r*rhogsol - w3r**2*w2i*w1r**2*rhogsol +&
        w3i*w1i*w2r*k*Kdrag*vdsol - w3i*w1i*w2r*k*Kdrag*vgsol - w3i*w1i*w2i*k**2*cs**2*rhogsol +&
        w3i*rhogeq*w2r**2*k*w1r*vgsol - w3i*w2r*w1r*k**2*cs**2*rhogsol + w3r**2*w1i*w2r**2*rhogsol -&
        w3r**2*w1i**2*rhogsol*w2i + w3i**2*w1i*w2r**2*rhogsol + w3i**2*w1i*rhogsol*w2i**2 - w3i*w2i**2*w1r**2*rhogsol -&
        w3i*w2r**2*w1r**2*rhogsol - w3i**2*w2i*w1r**2*rhogsol - w3i*w1i**2*rhogsol*w2i**2 - w3i**2*w1i**2*rhogsol*w2i -&
        w3i*w1i**2*w2r**2*rhogsol - w2i*rhogeq*k*w1r*vgsol*w1i**2 - w2i*rhogeq*k*w1r**3*vgsol +&
        w3i**2*w2i*rhogeq*k*w1r*vgsol + w2i*rhogeq*w3r**2*k*w1r*vgsol - w2r*w1i*rhogeq*w3i**2*k*vgsol -&
        w2r*w1i*rhogeq*w3r**2*k*vgsol + w2i*w3r*w1r**3*rhogsol + w2i*k**2*cs**2*rhogsol*w1r**2 +&
        w2i*w1i*w3i*rhogsol*w1r**2 + w2i*k**2*cs**2*rhogsol*w1i**2 + w2i*w3r*w1r*rhogsol*w1i**2 +&
        w2i*w1i**3*w3i*rhogsol + k**2*cs**2*rhogsol*w3i*w1r**2 + Kdrag*w1r**2*k*vgsol*w3r +&
        k**2*cs**2*rhogsol*w1i**2*w3i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w2i**2 + w1r**2&
        - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r)/rhogeq/k 

  vg1i =(w2i**2*rhogeq*w1r**2*k*vgsol + w2r**2*rhogeq*w1r**2*k*vgsol + w2r*w1r**3*w3r*rhogsol +&
        k**2*cs**2*w1r**3*rhogsol - w2i*w3i*w1r**3*rhogsol + w2i*w1i**3*w3r*rhogsol + w1i**2*w2i**2*rhogeq*k*vgsol +&
        w1i**2*w2r**2*rhogeq*k*vgsol - w2i*w3i*w1r*rhogsol*w1i**2 - w2i*k*Kdrag*vdsol*w1i**2 - w2i*k*Kdrag*vdsol*w1r**2&
        + w2i*w1i*w3r*rhogsol*w1r**2 - w2i*w1i*rhogeq*k*vgsol*w1r**2 - w2i*w1i**3*rhogeq*k*vgsol +&
        w2i*k*Kdrag*vgsol*w1r**2 + w2i*k*Kdrag*vgsol*w1i**2 - w3i*w2r*w1r*k*Kdrag*vgsol +&
        w2r*w1r*w3r*k**2*cs**2*rhogsol + w1r*w3r**2*w2r**2*rhogsol - w3r*w2r**2*w1r**2*rhogsol +&
        w3i**2*w1r*w2r**2*rhogsol + w2i**2*w3i**2*w1r*rhogsol - w2i**2*w3r*w1r**2*rhogsol + w2i**2*w1r*w3r**2*rhogsol -&
        w1i**2*w2i**2*w3r*rhogsol - w1i**2*w3i**2*w2r*rhogsol - w1i**2*w2r**2*w3r*rhogsol - w1i**2*w3r**2*w2r*rhogsol +&
        w3i*w2r*w1r*k*Kdrag*vdsol - w1r**2*w3r**2*w2r*rhogsol - w3i**2*w1r**2*w2r*rhogsol +&
        2*w1i**2*w3i*w2i*rhogeq*k*vgsol + w1i**2*rhogeq*w3i**2*k*vgsol + 2*w1i**2*rhogeq*w2r*w3r*k*vgsol +&
        w1i**2*rhogeq*w3r**2*k*vgsol + 2*w2i*rhogeq*w3i*w1r**2*k*vgsol + w2i*w1r*w3r*k*Kdrag*vdsol -&
        w2i*k**2*cs**2*w3i*w1r*rhogsol - w2i*w3r*Kdrag*w1r*k*vgsol - w2i**2*rhogeq*w3r*k*w1r*vgsol -&
        rhogeq*w3r*w2r**2*k*w1r*vgsol + 2*rhogeq*w3r*w2r*w1r**2*k*vgsol + rhogeq*w3i**2*w1r**2*k*vgsol -&
        rhogeq*w3i**2*k*vgsol*w1r*w2r - rhogeq*w3r**2*k*vgsol*w1r*w2r + rhogeq*w1r**2*w3r**2*k*vgsol -&
        w1i*w2r**2*w3i*rhogeq*k*vgsol - w1i*w2i**2*w3i*rhogeq*k*vgsol - w1i*w2i*rhogeq*w3i**2*k*vgsol -&
        w1i*w2i*rhogeq*w3r**2*k*vgsol + w1i*w2i*k*w3i*Kdrag*vdsol + w1i*w2i*w3r*k**2*cs**2*rhogsol +&
        w1i*w3i*w2r*k**2*cs**2*rhogsol - w1i*w2r*w3r*k*Kdrag*vdsol + w1i*w2r*w3r*k*Kdrag*vgsol -&
        w1i*w2i*k*w3i*Kdrag*vgsol - k*Kdrag*vdsol*w3i*w1r**2 - rhogeq*w1r**3*k*vgsol*w3r -&
        k**2*cs**2*w1r**2*rhogsol*w3r + k**2*cs**2*w1r*rhogsol*w1i**2 - k*Kdrag*vdsol*w1i**2*w3i -&
        w1i**2*rhogeq*k*vgsol*w1r*w3r - w1i*rhogeq*k*vgsol*w3i*w1r**2 + k*Kdrag*vgsol*w3i*w1r**2 -&
        w3r*k**2*cs**2*rhogsol*w1i**2 - w2r*k**2*cs**2*rhogsol*w1r**2 - w2r*k**2*cs**2*rhogsol*w1i**2 -&
        w1i*k*Kdrag*vgsol*w1r**2 + w1i*w2r*w3i*rhogsol*w1r**2 + w1i*k*Kdrag*vdsol*w1r**2 + w1i**3*k*Kdrag*vdsol +&
        w1i**3*w2r*w3i*rhogsol - w1i**3*k*Kdrag*vgsol - w1r**3*rhogeq*w2r*k*vgsol + w3i*k*Kdrag*vgsol*w1i**2 +&
        w2r*w1r*w3r*rhogsol*w1i**2 - w1r*rhogeq*w2r*k*vgsol*w1i**2 - w1i**3*rhogeq*k*vgsol*w3i)/( - 2*w1i*w3i + w1i**2&
        + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r)/rhogeq/k 

!-------------------------------
! D U S T  V E L O C I T I E S
!-------------------------------
  if (Kdrag.gt.0.) then
  vd3r =( - rhogsol*Kdrag*w2r**2*w1r**2*w3i + w3r*w2i*Kdrag**2*w1i*k*vdsol - w3r*w2i*Kdrag**2*w1i*k*vgsol +&
        w3r*Kdrag**2*w2r*w1r*k*vgsol - rhogsol*cs**4*k**4*rhogeq*w2r*w1r + w3r*w2i*Kdrag*k**2*cs**2*w1r*rhogsol -&
        w3r*w2i*Kdrag*w1i**2*rhogeq*k*vgsol - w3r*Kdrag*w1i*w2r**2*rhogeq*k*vgsol + w3r*rhogeq**2*w3i**2*k*vgsol*w2r**2&
        + rhogeq**2*w2i**2*w3r**3*k*vgsol - w3r*rhogeq*w2r*w1i**2*w3i**2*rhogsol + vgsol*Kdrag*k*rhogeq*w3i*w1r*w2i**2&
        + Kdrag*w3i*w1r*rhogsol*w2r*w3r**2 + Kdrag*w3i**2*w1i*rhogsol*w2i**2 - w3i*Kdrag*w1i*w2i*w3r**2*rhogsol +&
        Kdrag*w3r**2*k**2*cs**2*rhogsol*w3i - rhogeq*w1r*w3r**3*k**2*cs**2*rhogsol -&
        w3r*rhogeq*w2r*w3i**2*rhogsol*w1r**2 - w3r*rhogeq*w1i*k*w3i**2*Kdrag*vdsol +&
        w3r*rhogeq**2*w3i**2*w1r**2*k*vgsol + w3r*rhogeq**2*w1i**2*w3i**2*k*vgsol - w3r*w2i*Kdrag*rhogeq*w1r**2*k*vgsol&
        - rhogeq*w2r*rhogsol*w1r**2*w3r**3 - rhogeq*w1r*rhogsol*w2r**2*w3r**3 - rhogeq*w2r*w1i**2*rhogsol*w3r**3 +&
        2*w3r*rhogeq**2*w3i**2*k*vgsol*w1r*w2r - w3r*rhogeq*w2r*k**2*cs**2*rhogsol*w3i**2 +&
        2*w3r*rhogeq**2*w3i**2*w1i*w2i*k*vgsol - w3r*rhogeq*w3i**2*w2i**2*w1r*rhogsol -&
        w3r*rhogeq*w3i**2*w1r*rhogsol*w2r**2 + w3r*rhogeq**2*w3i**2*w2i**2*k*vgsol -&
        w3r*rhogeq*w3i**2*w2i*k*Kdrag*vdsol - rhogeq*w2i*w3r**3*k*Kdrag*vdsol + w3r*Kdrag*w1i*w2r*k**2*cs**2*rhogsol -&
        w3r*rhogeq*w2i**2*w1i*k*Kdrag*vgsol - vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r +&
        rhogsol*cs**2*k**2*rhogeq*w2i*w3i*w1r**2 - rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**2 -&
        rhogeq*w2i**2*w1r*w3r**3*rhogsol - rhogsol*Kdrag*w3i*w2i**2*w1r**2 + w1i*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2 +&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r - rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**2 +&
        rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**2 + w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i -&
        w1i**2*rhogsol*Kdrag*w3i*w2i**2 - w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2 - rhogsol*Kdrag*w1i**2*w3i*w2r**2 -&
        rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2 + vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r -&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r + vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r +&
        rhogeq**2*w3i**2*w2i**2*k*w1r*vgsol + vgsol*Kdrag*k*rhogeq*w3i*w1r*w2r**2 + w1i*rhogsol*cs**4*k**4*rhogeq*w2i -&
        Kdrag**2*w3r**2*w2r*k*vgsol + rhogeq*w3i**4*k**2*cs**2*rhogsol + Kdrag**2*w3r**2*k*vdsol*w1r -&
        w3i*rhogeq*k**4*cs**4*w1i*rhogsol + rhogeq*w3i**3*w2i*w1r**2*rhogsol - rhogeq*w3i**4*w1i*w2i*rhogsol +&
        rhogeq*w3i**3*w1i**2*rhogsol*w2i + rhogeq*w3i**3*w1i*w2r**2*rhogsol - Kdrag*w3r**3*w2i*w1r*rhogsol +&
        Kdrag*w3r**2*w1i**2*rhogsol*w2i + Kdrag**2*w3r**2*w2r*k*vdsol - Kdrag*w3r**2*w2i*k**2*cs**2*rhogsol +&
        Kdrag**2*w3r**3*k*vgsol - w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol + w3i*rhogeq*k**3*cs**2*Kdrag*vdsol*w1r -&
        w3i**2*rhogeq**2*k**3*cs**2*w2r*vgsol - w3i**2*rhogeq**2*k**3*cs**2*w1r*vgsol +&
        w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol + rhogeq*w3i**2*w1i*w2r*k*Kdrag*vgsol + rhogeq*w3i**3*w1i*rhogsol*w2i**2&
        - rhogeq*w3i**2*w1i**2*w2r**2*rhogsol + rhogeq*w3i**3*k*Kdrag*vdsol*w1r - rhogeq*w3i**2*w1i**2*rhogsol*w2i**2 -&
        rhogeq*w3i**2*w2i*w1r*k*Kdrag*vdsol + rhogeq**2*w3i**2*w1i**2*w2r*k*vgsol +&
        rhogeq*w3i*w3r**2*w1i*rhogsol*w2i**2 - w3i*rhogeq*k**3*cs**2*Kdrag*vgsol*w1r -&
        rhogeq*w3i*w3r**2*w1i*k**2*cs**2*rhogsol + rhogeq*w3i*w3r**2*k*Kdrag*vdsol*w1r +&
        rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol + Kdrag*w3r**2*w1i*w2r**2*rhogsol + Kdrag*w3r**2*w2i*w1r**2*rhogsol -&
        rhogeq**2*w3i**4*w2r*k*vgsol - rhogeq**2*w3i**4*k*w1r*vgsol - rhogeq*w3i**3*w1i*k**2*cs**2*rhogsol -&
        rhogeq*w3i**2*w1i*w2r*k*Kdrag*vdsol + w3i**2*rhogeq*k**4*cs**4*rhogsol - Kdrag*w3r**3*w1i*w2r*rhogsol -&
        Kdrag**2*w3r**2*k*vgsol*w1r + Kdrag*w3r**2*w1i*rhogsol*w2i**2 - Kdrag*w3r**2*w1i*k**2*cs**2*rhogsol +&
        rhogeq**2*w3i**2*w2r**2*k*w1r*vgsol + rhogeq**2*w3i**2*w2r*w1r**2*k*vgsol + rhogeq*w3i**3*w2r*k*Kdrag*vdsol +&
        rhogeq*w3i**2*w2i*Kdrag*w1r*k*vgsol - rhogeq*w3i**3*w2i*k**2*cs**2*rhogsol -&
        2*Kdrag*w3i*w3r**2*rhogeq*w2r*k*vgsol - 2*Kdrag*w3i*w3r**2*rhogeq*k*w1r*vgsol +&
        Kdrag*w3i*w1i*w2i*k**2*cs**2*rhogsol + Kdrag**2*w3i**2*w3r*k*vgsol - Kdrag**2*w3i**2*w3r*k*vdsol +&
        Kdrag**2*w3i**2*k*vdsol*w1r - Kdrag*w3i**2*w2i*k**2*cs**2*rhogsol - Kdrag**2*w3i*w2i*w1r*k*vdsol -&
        Kdrag*w3i**2*w2i*w3r*w1r*rhogsol + Kdrag**2*w3i**2*w2r*k*vdsol + Kdrag*w3i**2*w2i*w1r**2*rhogsol +&
        Kdrag*w3i**2*w1i**2*rhogsol*w2i + Kdrag**2*w3i*w2i*w1r*k*vgsol + Kdrag*w3i**2*w1i*w2r**2*rhogsol +&
        Kdrag*w3i*rhogeq*w2r*w1r**2*k*vgsol - Kdrag*w3i**2*w1i*w3r*w2r*rhogsol + Kdrag*w3i*w1i**2*rhogeq*w2r*k*vgsol +&
        rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol - w3i*rhogeq*k**4*cs**4*w2i*rhogsol +&
        rhogeq*w3r*w2r**2*w1r*k**2*cs**2*rhogsol - Kdrag*w3i**3*w1i*w2i*rhogsol - Kdrag*w3i*w2r*w1r*k**2*cs**2*rhogsol&
        - 2*Kdrag*w3i**3*rhogeq*k*w1r*vgsol + Kdrag*w3i**3*k**2*cs**2*rhogsol - 2*rhogeq*w3r**2*w3i**2*w1i*w2i*rhogsol&
        - rhogeq*w3i**2*w2i**2*w1r**2*rhogsol + Kdrag*w3i**3*w2r*w1r*rhogsol - rhogeq*w3i**2*w2r**2*w1r**2*rhogsol -&
        2*rhogeq**2*w3r**2*w3i**2*w2r*k*vgsol - rhogeq*w3r**2*w1r*w2i*k*Kdrag*vgsol + rhogeq*w3i**4*w2r*w1r*rhogsol +&
        Kdrag**2*w3i*w1i*w2r*k*vgsol - 2*rhogeq**2*w3r**2*w3i**2*k*w1r*vgsol + rhogeq*w3r**2*w1r*w2i*k*Kdrag*vdsol +&
        2*rhogeq*w3r**2*w3i**2*k**2*cs**2*rhogsol + rhogeq*w3r*w2r*w1r**2*k**2*cs**2*rhogsol -&
        2*Kdrag*w3i**3*rhogeq*w2r*k*vgsol + 2*rhogeq**2*w3r*cs**2*k**3*w1i*w3i*vgsol -&
        rhogeq*w3r*cs**2*k**3*w1i*Kdrag*vgsol - Kdrag**2*w3i**2*k*vgsol*w1r - rhogeq**2*w3r*cs**2*k**3*w2i**2*vgsol -&
        rhogeq**2*w3r*cs**2*k**3*w1r**2*vgsol + rhogeq**2*w3r**2*cs**2*k**3*w1r*vgsol - Kdrag**2*w3r**3*k*vdsol +&
        rhogeq*w3i*w3r**2*w1i*w2r**2*rhogsol + rhogeq*w3r*cs**2*k**3*w1i*Kdrag*vdsol +&
        rhogeq*w3i*w3r**2*w2r*k*Kdrag*vdsol - rhogeq*w3i*w3r**2*w2i*k**2*cs**2*rhogsol -&
        2*rhogeq**2*w3r*cs**2*k**3*vgsol*w1r*w2r + rhogeq*w3i*w3r**2*w1i**2*rhogsol*w2i -&
        2*rhogeq**2*w3r*cs**2*k**3*w1i*w2i*vgsol + rhogeq*w3r*cs**2*k**3*w2i*Kdrag*vdsol +&
        rhogeq**2*w3r**2*cs**2*k**3*w2r*vgsol + rhogeq*w3i*w3r**2*w2i*w1r**2*rhogsol -&
        Kdrag*w3i**2*w1i*k**2*cs**2*rhogsol - Kdrag**2*w3i*w1i*w2r*k*vdsol + rhogeq*w3r*cs**4*k**4*w1r*rhogsol -&
        2*rhogeq**2*w3r*w3i*w1i*w2i**2*k*vgsol + rhogeq*w3r**2*w2r**2*w1r**2*rhogsol +&
        rhogeq*w3r*cs**4*k**4*w2r*rhogsol - 2*rhogeq*w3r*cs**2*k**3*w3i*Kdrag*vdsol +&
        2*rhogeq*w3r*cs**2*k**3*w3i*Kdrag*vgsol - rhogeq*w3r*cs**2*k**3*w2i*Kdrag*vgsol -&
        rhogeq**2*w3r**2*w1r*w2i**2*k*vgsol - rhogeq**2*w3r*cs**2*k**3*w1i**2*vgsol +&
        2*rhogeq**2*w3r*cs**2*k**3*w2i*w3i*vgsol - rhogeq**2*w3r*cs**2*k**3*w2r**2*vgsol -&
        2*rhogeq**2*w3r*w3i*w2i*w1r**2*k*vgsol + 2*rhogeq*w3r**3*w1i*k*Kdrag*vgsol + 2*rhogeq*w3r**3*w2i*k*Kdrag*vgsol&
        + rhogeq*w3r*w2i**2*k**2*cs**2*w1r*rhogsol + rhogeq*w3r**2*w2r**2*w1i**2*rhogsol -&
        rhogeq*w3r**2*w2r*w1i*k*Kdrag*vgsol - rhogeq*w3r**2*cs**4*k**4*rhogsol - 2*rhogeq**2*w3r*w1i*w3i*w2r**2*k*vgsol&
        - rhogeq**2*w3r**2*w2r**2*k*vgsol*w1r - 2*rhogeq*w3r*w3i*w2r*w1r*k*Kdrag*vdsol +&
        rhogeq*w3r**2*w1r**2*w2i**2*rhogsol + 2*rhogeq*w3r*w3i*w2r*w1r*k*Kdrag*vgsol -&
        rhogeq**2*w3r**2*w2r*w1i**2*k*vgsol - 2*rhogeq**2*w3r*w3i*w2i*w1i**2*k*vgsol +&
        2*rhogeq*w3r**2*w3i**2*w2r*w1r*rhogsol + rhogeq*w3r**2*w2i**2*w1i**2*rhogsol - Kdrag**2*w3i**2*w2r*k*vgsol +&
        rhogeq*w3r*w1i**2*w2r*k**2*cs**2*rhogsol - rhogeq**2*w3r**2*w2r*w1r**2*k*vgsol +&
        rhogeq*w3r**2*w2r*w1i*k*Kdrag*vdsol + 2*rhogeq*w3r*w3i**2*w2i*k*Kdrag*vgsol - rhogeq*w3r**4*w1i*w2i*rhogsol -&
        rhogeq**2*w3r**4*k*w1r*vgsol - rhogeq**2*w3r**4*w2r*k*vgsol + vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r +&
        rhogeq*w3r**4*k**2*cs**2*rhogsol + 2*rhogeq*w3r*w1i*w3i**2*k*Kdrag*vgsol -&
        2*rhogeq*w3r*w3i*w2i*w1i*k*Kdrag*vgsol + 2*rhogeq*w3r*w3i*w2i*w1i*k*Kdrag*vdsol + rhogeq*w3r**4*w2r*w1r*rhogsol&
        + 2*rhogeq**2*w1i*w2i*w3r**3*k*vgsol - w3r*Kdrag**2*w2r*w1r*k*vdsol - w3r*rhogeq*w3i**2*w1r*k**2*cs**2*rhogsol&
        + rhogeq**2*k*vgsol*w3r**3*w2r**2 - rhogeq*w3r**3*k**2*cs**2*rhogsol*w2r + 2*rhogeq**2*w3r**3*k*w1r*vgsol*w2r -&
        rhogeq*w1i*k*Kdrag*vdsol*w3r**3 + rhogeq**2*w1i**2*k*vgsol*w3r**3 +&
        rhogeq**2*w1r**2*k*vgsol*w3r**3)/k/Kdrag/rhogeq/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 -&
        2*w3i*w2i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  vd3i =(2*rhogeq*k*Kdrag*vdsol*w3i**2*w3r**2 - rhogeq**2*w3i**2*w2i*w1r**2*k*vgsol + rhogeq*w3i**4*w2i*w1r*rhogsol -&
        2*rhogeq**2*w3i**2*w2i*w3r**2*k*vgsol - rhogeq**2*w3i**4*w2i*k*vgsol - rhogeq*w3r*w2i*w1i**2*w3i**2*rhogsol -&
        rhogeq*w3r**3*w2i*w1r**2*rhogsol - rhogeq*w3i**3*w2i**2*w1r*rhogsol +&
        4*rhogeq*w3i*w2r*w1r*w3r*k**2*cs**2*rhogsol + 2*rhogeq*w3i*w1i**2*w2i**2*w3r*rhogsol +&
        2*rhogeq*w3i*w2i**2*w3r*w1r**2*rhogsol - rhogeq*w2r*k**2*cs**2*rhogsol*w3i**3 +&
        w3r**2*rhogeq**2*w1i**2*w2i*k*vgsol - 2*cs**2*k**3*rhogeq**2*w3i*w3r*w1r*vgsol +&
        2*cs**4*k**4*rhogeq*w3i*w3r*rhogsol + Kdrag**2*w2i*w3r**2*k*vgsol - w3i**2*Kdrag*w1i*w2i*w3r*rhogsol +&
        w3i**2*Kdrag**2*w2i*k*vgsol - w3i**2*Kdrag**2*w2i*k*vdsol - Kdrag*w1i*w2i*w3r**3*rhogsol +&
        w3i*Kdrag*w2i*w3r**2*w1r*rhogsol - Kdrag**2*w2i*w3r**2*k*vdsol - 4*rhogeq*w3i*w1i*w2i*w3r*k**2*cs**2*rhogsol +&
        rhogeq*w3i**2*w1i*w2i*k*Kdrag*vdsol - rhogeq*w3i*w2i**2*w1r*w3r**2*rhogsol +&
        2*rhogeq**2*w3i*w1i*w2i*w3r**2*k*vgsol + 2*rhogeq**2*w3i**3*w1i*w2i*k*vgsol -&
        2*rhogeq**2*w3i**2*w1i*w3r**2*k*vgsol - rhogeq**2*w3i**4*w1i*k*vgsol + w3i**3*Kdrag*w2i*w1r*rhogsol -&
        rhogeq*w3i*w2i*w3r**2*k*Kdrag*vdsol + 2*rhogeq*w3i*w2i*w1r*w3r*k*Kdrag*vdsol +&
        2*rhogeq*w3i**2*w2i*w3r**2*w1r*rhogsol + rhogeq**2*w3i*w2i**2*w3r**2*k*vgsol -&
        2*rhogeq**2*w3i*w2i**2*w3r*k*w1r*vgsol + rhogeq**2*w3i**3*w2i**2*k*vgsol - rhogeq*w3i**3*w2i*k*Kdrag*vdsol -&
        2*rhogeq*w3i*w2i*w3r*Kdrag*w1r*k*vgsol - rhogeq*w3i**3*w1r*k**2*cs**2*rhogsol +&
        2*rhogeq**2*w3i**3*k*vgsol*w1r*w2r - rhogeq**2*w3i**2*w1i**2*w2i*k*vgsol -&
        rhogeq*w3i*w1r*w3r**2*k**2*cs**2*rhogsol - Kdrag*rhogeq*w3r**4*k*vgsol - rhogeq*w3r*w2i**2*w1i*w3i**2*rhogsol -&
        rhogeq*w3r**3*w2i*w1i**2*rhogsol - rhogeq*w3r**3*w2i**2*w1i*rhogsol + rhogeq*w3r**3*w2i*k**2*cs**2*rhogsol -&
        w3r**4*rhogeq**2*w1i*k*vgsol + w3r**4*rhogeq*w2i*w1r*rhogsol + w3r**2*rhogeq**2*w2i*w1r**2*k*vgsol -&
        w3r**2*rhogeq*w1i*w2i*k*Kdrag*vdsol - w3r**4*rhogeq**2*w2i*k*vgsol + 2*cs**2*k**3*rhogeq**2*w3i*vgsol*w2r*w1r -&
        Kdrag**2*w3i**3*k*vgsol + Kdrag**2*w3i**3*k*vdsol - rhogeq*w3r*w2i*w3i**2*rhogsol*w1r**2 +&
        rhogeq*w3r*w2i*w3i**2*k**2*cs**2*rhogsol + Kdrag*w3i**2*w1r*rhogsol*w2r*w3r + Kdrag*w1r*w3r**3*rhogsol*w2r +&
        cs**2*k**3*rhogeq*Kdrag*vdsol*w3i**2 - cs**2*k**3*rhogeq*Kdrag*vdsol*w3r**2 - rhogeq*w3i**3*w1r*rhogsol*w2r**2&
        - rhogeq*w1i*w3r**3*rhogsol*w2r**2 + rhogeq*w1i*w3r**4*rhogsol*w2r - rhogeq*w2i**2*k**2*cs**2*w3i*w1r*rhogsol -&
        rhogeq*w2i**2*w3r*Kdrag*w1r*k*vgsol + Kdrag*w2r*w1r*w3r*k**2*cs**2*rhogsol - Kdrag*w1r*w3r**2*w2r**2*rhogsol +&
        Kdrag*w3r*w2r**2*w1r**2*rhogsol - Kdrag*w3i**2*w1r*w2r**2*rhogsol - w2i**2*Kdrag*w3i**2*w1r*rhogsol +&
        w2i**2*Kdrag*w3r*w1r**2*rhogsol - w2i**2*Kdrag*w1r*w3r**2*rhogsol + w2i**2*Kdrag*w1i**2*w3r*rhogsol -&
        Kdrag*w1i**2*w3i**2*w2r*rhogsol + Kdrag*w1i**2*w2r**2*w3r*rhogsol - Kdrag*w1i**2*w3r**2*w2r*rhogsol -&
        Kdrag**2*w3i*w2r*w1r*k*vdsol - Kdrag*w1r**2*w3r**2*w2r*rhogsol + 2*cs**2*k**3*rhogeq**2*w1i*w2i*w3i*vgsol -&
        cs**2*k**3*rhogeq**2*w1i*w3i**2*vgsol + cs**2*k**3*rhogeq**2*w1i*w3r**2*vgsol -&
        cs**2*k**3*rhogeq*w1i*w3i*Kdrag*vdsol - cs**4*k**4*rhogeq*w1i*w3r*rhogsol +&
        cs**2*k**3*rhogeq*w1i*w3i*Kdrag*vgsol + cs**2*k**3*rhogeq**2*w3i*w2r**2*vgsol -&
        cs**2*k**2*rhogeq*w3i*w1r*rhogsol*w2r**2 + cs**2*k**3*rhogeq*w2i*w3i*Kdrag*vgsol -&
        cs**4*k**4*rhogeq*w2i*w3r*rhogsol - cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol +&
        cs**2*k**2*rhogeq*w1i*w2r**2*w3r*rhogsol + cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol -&
        cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol - cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol -&
        cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol + 2*cs**2*k**2*rhogeq*w2i*w3i**2*w1r*rhogsol +&
        cs**2*k**2*rhogeq*w2i*w3r*w1r**2*rhogsol - 2*cs**2*k**2*rhogeq*w2i*w1r*w3r**2*rhogsol +&
        cs**2*k**2*rhogeq*w1i**2*w2i*w3r*rhogsol + cs**2*k**3*rhogeq**2*w1i**2*w3i*vgsol +&
        cs**2*k**3*rhogeq**2*w3i*w1r**2*vgsol + cs**2*k**3*rhogeq*w1r*w3r*Kdrag*vdsol -&
        cs**4*k**4*rhogeq*w3i*w1r*rhogsol - cs**2*k**3*rhogeq*w3r*Kdrag*w1r*vgsol -&
        cs**2*k**3*rhogeq**2*w2i*w3i**2*vgsol + cs**2*k**3*rhogeq**2*w2i*w3r**2*vgsol -&
        rhogeq**2*w2i**2*w1i*w3i**2*k*vgsol + rhogeq**2*w2i**2*w1i*w3r**2*k*vgsol +&
        rhogeq*w2i**2*w1i*w3r*k**2*cs**2*rhogsol - rhogeq*w2i**2*w1i*k*w3i*Kdrag*vgsol -&
        cs**2*k**3*rhogeq*w2i*w3i*Kdrag*vdsol + cs**2*k**3*rhogeq**2*w2i**2*w3i*vgsol -&
        cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol + cs**4*k**4*rhogeq*w2i*w1r*rhogsol + Kdrag**2*w3i*w2r*w1r*k*vgsol -&
        Kdrag*w3i**2*w1r**2*w2r*rhogsol + 3*w2i*Kdrag*w1i*rhogeq*w3r**2*k*vgsol + w2i*Kdrag**2*w1i*k*w3i*vdsol -&
        w2i*Kdrag*w1i*w3r*k**2*cs**2*rhogsol + Kdrag*w1i*w3i*w2r*k**2*cs**2*rhogsol + Kdrag**2*w1i*w2r*w3r*k*vdsol -&
        Kdrag**2*w1i*w2r*w3r*k*vgsol - w2i*Kdrag**2*w1i*k*w3i*vgsol - w2r**2*rhogeq**2*w1i*w3i**2*k*vgsol +&
        w2r**2*rhogeq**2*w1i*w3r**2*k*vgsol - rhogeq*w2r*k**2*cs**2*rhogsol*w3i*w1r**2 -&
        2*rhogeq*w2r*w1i*w3r**2*k**2*cs**2*rhogsol + 2*rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w3i**2 -&
        rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i - rhogeq*w2r*cs**2*k**3*w3r*Kdrag*vgsol +&
        rhogeq*w2r*cs**4*k**4*w1i*rhogsol + rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol - rhogeq*w2r*cs**4*k**4*w3i*rhogsol -&
        rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol + rhogeq*w2r*cs**2*k**3*w3r*Kdrag*vdsol + Kdrag**2*w3i*k*vdsol*w3r**2 -&
        w2i*Kdrag*w1i**2*w3i*rhogeq*k*vgsol - Kdrag*w1i**2*rhogeq*w2r*w3r*k*vgsol - w2i*Kdrag*rhogeq*w3i*w1r**2*k*vgsol&
        + w2i*Kdrag**2*w1r*w3r*k*vdsol + w2i*Kdrag*k**2*cs**2*w3i*w1r*rhogsol - w2i*Kdrag**2*w3r*w1r*k*vgsol +&
        w2i**2*Kdrag*rhogeq*w3i**2*k*vgsol + w2i**2*Kdrag*rhogeq*w3r**2*k*vgsol + Kdrag*rhogeq*w2r**2*w3i**2*k*vgsol -&
        Kdrag*rhogeq*w3r*w2r**2*k*w1r*vgsol - Kdrag*rhogeq*w3r*w2r*w1r**2*k*vgsol +&
        3*Kdrag*rhogeq*w3i**2*k*vgsol*w1r*w2r + Kdrag*rhogeq*w3r**2*k*vgsol*w1r*w2r +&
        Kdrag*rhogeq*w2r**2*w3r**2*k*vgsol - Kdrag*w1i*w2r**2*w3i*rhogeq*k*vgsol + w2i*Kdrag*w1i*rhogeq*w3i**2*k*vgsol&
        + rhogeq**2*w3i**3*k*vgsol*w2r**2 + cs**2*k**3*rhogeq*Kdrag*vgsol*w3r**2 + Kdrag**2*w1i*k*vgsol*w3r**2 +&
        Kdrag**2*w1i*k*vgsol*w3i**2 - Kdrag**2*w1i*k*vdsol*w3i**2 - Kdrag**2*w1i*k*vdsol*w3r**2 -&
        cs**2*k**3*rhogeq*Kdrag*vgsol*w3i**2 + Kdrag*w1i**2*rhogeq*k*vgsol*w3i**2 + Kdrag*w1i**2*rhogeq*k*vgsol*w3r**2&
        + Kdrag*w3r**3*k**2*cs**2*rhogsol + Kdrag*w1i*w2r*w3i**3*rhogsol - Kdrag**2*w3i*k*vgsol*w3r**2 +&
        Kdrag*w1i*w2r*w3i*rhogsol*w3r**2 + Kdrag*w3r*k**2*cs**2*rhogsol*w3i**2 - Kdrag*w3r**2*k**2*cs**2*rhogsol*w2r -&
        Kdrag*k**2*cs**2*w1r*rhogsol*w3r**2 + rhogeq**2*w3i**3*w1r**2*k*vgsol + rhogeq**2*w1i**2*w3i**3*k*vgsol +&
        Kdrag*rhogeq*w1r**2*k*vgsol*w3r**2 - Kdrag*k**2*cs**2*w1r*rhogsol*w3i**2 - Kdrag*w2r*k**2*cs**2*rhogsol*w3i**2&
        - 2*rhogeq**2*w1i**2*w3i*k*vgsol*w2r*w3r - 2*rhogeq**2*w3i*w1r**2*k*vgsol*w2r*w3r -&
        rhogeq*w1r*w3r*k*Kdrag*vdsol*w3i**2 + rhogeq*w1r*w3r**2*k*Kdrag*vdsol*w2r + Kdrag*rhogeq*w1r**2*k*vgsol*w3i**2&
        - rhogeq*w2r*w1r*k*Kdrag*vdsol*w3i**2 + rhogeq*w2r*w1i*rhogsol*w3i**4 - rhogeq*w2r*w1i**2*w3i**3*rhogsol -&
        rhogeq*w1i*w2r**2*w3r*rhogsol*w3i**2 - rhogeq*w2r*w3r*k*Kdrag*vdsol*w3i**2 -&
        rhogeq*w1i*k*w3i*Kdrag*vdsol*w3r**2 + rhogeq*w1i*w3r*k**2*cs**2*rhogsol*w3i**2 -&
        rhogeq*w1r*w3r**3*k*Kdrag*vdsol - rhogeq*w2r*w3i**3*rhogsol*w1r**2 + rhogeq*w1i*w3r**3*k**2*cs**2*rhogsol +&
        rhogeq**2*w1i**2*w3i*k*vgsol*w3r**2 - Kdrag*rhogeq*w3i**4*k*vgsol + rhogeq**2*w3i*w1r**2*k*vgsol*w3r**2 +&
        rhogeq*k*Kdrag*vdsol*w3i**4 + 2*rhogeq*w1i*k*w3i*Kdrag*vdsol*w2r*w3r - rhogeq*w1i*k*w3i**3*Kdrag*vdsol +&
        2*rhogeq*w2r*w1i*rhogsol*w3r**2*w3i**2 + rhogeq*k*Kdrag*vdsol*w3r**4 - rhogeq*w3i*w1r*rhogsol*w2r**2*w3r**2 +&
        2*rhogeq*w2r**2*w1i**2*w3i*rhogsol*w3r - w3i*rhogeq*w3r**2*k**2*cs**2*rhogsol*w2r -&
        2*w3i*rhogeq**2*w3r*k*w1r*vgsol*w2r**2 + 2*w3i*rhogeq**2*w3r**2*k*w1r*vgsol*w2r -&
        rhogeq*w2r*w1i**2*w3i*rhogsol*w3r**2 - 2*rhogeq*w1i*k*w3i*Kdrag*vgsol*w2r*w3r +&
        rhogeq**2*w3i*k*vgsol*w3r**2*w2r**2 - rhogeq*w2r*w3i*rhogsol*w1r**2*w3r**2 +&
        2*rhogeq*w2r**2*w3i*rhogsol*w1r**2*w3r - rhogeq*w2r*w3r**3*k*Kdrag*vdsol -&
        2*cs**2*k**3*rhogeq**2*w3i*vgsol*w2r*w3r - 2*Kdrag*rhogeq*w3i**2*k*vgsol*w3r**2)/k/Kdrag/rhogeq/( - 2*w2r*w3r +&
        w2r**2 + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  vd2r =(w2i**3*w1i**2*rhogeq*rhogsol*w3i + w2i**3*Kdrag*k**2*cs**2*rhogsol + w2i**3*rhogsol*Kdrag*w3r*w1r +&
        w2i**3*vdsol*Kdrag*k*rhogeq*w3r + Kdrag**2*w2r**3*k*vgsol - Kdrag**2*k*vdsol*w2r**3 -&
        w2i**3*rhogsol*cs**2*k**2*rhogeq*w3i - w2i**4*vgsol*k*rhogeq**2*w3r - w2r*cs**2*k**2*rhogeq*w2i**2*w3r*rhogsol&
        - cs**2*k**2*rhogeq*w3r*rhogsol*w2r**3 + 2*w2r*w2i**2*Kdrag*w3i*rhogeq*k*vgsol -&
        2*w2i**3*vgsol*Kdrag*k*rhogeq*w3r - rhogeq*w2r**2*k**2*cs**2*rhogsol*w3i*w2i -&
        2*Kdrag*rhogeq*w2r**2*w3r*k*vgsol*w2i - rhogeq**2*w2r**4*w3r*k*vgsol + 2*Kdrag*w3i*rhogeq*k*vgsol*w2r**3 -&
        vgsol*k*Kdrag**2*w2r**2*w1r - rhogeq*w3i**2*w1r*w2r**3*rhogsol - rhogeq*w3r*w2r**3*w1r**2*rhogsol -&
        rhogeq**2*w2r**4*k*vgsol*w1r + rhogeq**2*w1r**2*k*vgsol*w2r**3 - rhogeq*k**2*cs**2*w1r*rhogsol*w2r**3 +&
        rhogeq**2*w3r**2*k*vgsol*w2r**3 - rhogeq*w3i*k*Kdrag*vdsol*w2r**3 - rhogeq*w1r*w3r**2*w2r**3*rhogsol +&
        2*rhogeq**2*w3r*k*w1r*vgsol*w2r**3 + rhogeq*w2r**4*k**2*cs**2*rhogsol + rhogeq*w2r**4*w3r*rhogsol*w1r +&
        Kdrag**2*w2r**2*w3r*k*vdsol - Kdrag**2*w2r**2*w3r*k*vgsol - Kdrag*w3i*w2r**2*k**2*cs**2*rhogsol +&
        w2r*cs**2*k**2*rhogeq*w1r*w3r**2*rhogsol + w2r*cs**2*k**2*rhogeq*w3i**2*w1r*rhogsol +&
        rhogeq**2*w3i**2*k*vgsol*w2r**3 - Kdrag*w3i*w1r*rhogsol*w2r**3 - rhogsol*cs**2*k**2*rhogeq*w3i**2*w1r**2 -&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i**2 - 2*w1i*vgsol*cs**2*k**3*rhogeq**2*w3i*w2r -&
        w1i*vdsol*Kdrag*k*rhogeq*w2r**3 + 2*w1i*vgsol*k*rhogeq**2*w3i*w2r**3 + 2*w1i*vgsol*Kdrag*k*rhogeq*w2r**3 -&
        w1i*Kdrag*w3r*rhogsol*w2r**3 - w1i*rhogsol*cs**2*k**2*Kdrag*w2r**2 - w1i*rhogeq*w2r**4*w3i*rhogsol -&
        w2i**3*w1i*rhogsol*Kdrag*w3i + w2i**4*rhogeq*rhogsol*w3r*w1r + w2i**4*rhogeq*k**2*cs**2*rhogsol -&
        w2i**4*w1i*rhogeq*rhogsol*w3i + w2i**3*w1i*rhogeq*rhogsol*w3r**2 - w2i**3*w1i*rhogsol*cs**2*k**2*rhogeq +&
        w2i**3*vdsol*Kdrag*k*rhogeq*w1r + w2i**3*rhogeq*rhogsol*w3i*w1r**2 + w2i*w1i*rhogsol*cs**2*k**2*rhogeq*w3i**2 +&
        2*w2i*cs**2*k**3*rhogeq**2*w3i*vgsol*w2r - 2*w2i**3*vgsol*Kdrag*k*rhogeq*w1r - w2i*w1i*Kdrag*w2r**2*w3i*rhogsol&
        + 2*w2i*w1i*w2r*cs**2*k**3*rhogeq**2*vgsol - w2i*w1i*rhogsol*cs**2*k**2*rhogeq*w2r**2 +&
        w2i*w1i*rhogsol*cs**2*k**2*rhogeq*w3r**2 + w2i*w1i*rhogeq*w3i**2*w2r**2*rhogsol -&
        2*w2i*w1i**2*w2r*rhogeq**2*w3i*k*vgsol - 2*w2i*w1i*vgsol*Kdrag*k*rhogeq*w3i*w2r +&
        w2i*w1i*rhogeq*w3r**2*w2r**2*rhogsol - 2*w2i*w1i*w2r*rhogeq**2*w3i**2*k*vgsol +&
        2*w2i*w1i*w2r*rhogeq*k*w3i*Kdrag*vdsol - 2*w2i*w1i*w2r*rhogeq**2*w3r**2*k*vgsol -&
        2*w2i*w2r*cs**2*k**3*rhogeq*Kdrag*vdsol - 2*w2i*w2r*rhogeq*w1r*w3r*k*Kdrag*vdsol +&
        w2i*w1i**2*rhogeq*w2r**2*w3i*rhogsol + w2i*rhogeq*w2r**2*w3r*k*Kdrag*vdsol +&
        2*w2i*w2r*cs**2*k**3*rhogeq*Kdrag*vgsol - 2*w2i*w2r*rhogeq**2*w3i*w1r**2*k*vgsol +&
        w2i*Kdrag*w2r**2*k**2*cs**2*rhogsol + w2i*Kdrag*w2r**2*w1r*w3r*rhogsol + w2i*rhogeq*w2r**2*w3i*rhogsol*w1r**2 +&
        w2i*rhogeq*w2r**2*w1r*k*Kdrag*vdsol + 2*w2i*vgsol*Kdrag*k*rhogeq*w3r*w2r*w1r -&
        2*w2i*Kdrag*rhogeq*k*vgsol*w2r**2*w1r - w2i**4*vgsol*k*rhogeq**2*w1r - w2r**2*vgsol*k*rhogeq**2*w1r**2*w3r -&
        w2r**2*cs**4*k**4*rhogeq*rhogsol - w2r**2*w1i**2*vgsol*k*rhogeq**2*w3r + w2r**2*w1i*vdsol*Kdrag*k*rhogeq*w3r +&
        w2r**2*vdsol*Kdrag*k*rhogeq*w3i*w1r - rhogsol*cs**2*k**2*rhogeq*w1r**2*w3r**2 -&
        2*w2i**2*w2r**2*vgsol*k*rhogeq**2*w3r - w2i**2*vdsol*Kdrag*k*rhogeq*w3i*w1r + w2i**2*cs**4*k**4*rhogeq*rhogsol&
        + 2*w2i**2*w1i*vgsol*Kdrag*k*rhogeq*w2r - 2*w2i**2*w1i*rhogeq*w2r**2*w3i*rhogsol -&
        w2i**2*w1i*w2r*Kdrag*w3r*rhogsol + 2*w2i**2*w1i*vgsol*k*rhogeq**2*w3i*w2r - w2i**2*w1i*rhogsol*cs**2*k**2*Kdrag&
        + w2i**2*w1i**2*w2r*rhogeq**2*k*vgsol - w2i**2*w1i*vdsol*Kdrag*k*rhogeq*w3r -&
        w2i**2*w1i*vdsol*Kdrag*k*rhogeq*w2r - w2i**2*w2r*rhogeq*k**2*cs**2*w1r*rhogsol -&
        w2i**2*w2r*Kdrag*w3i*w1r*rhogsol + w2i**2*w2r*rhogeq**2*w1r**2*k*vgsol + 2*w2i**2*vgsol*k*rhogeq**2*w3r*w2r*w1r&
        - w2i**2*w1i**2*w2r*rhogeq*w3r*rhogsol + w2i**2*w1i**2*vgsol*k*rhogeq**2*w3r -&
        w2i**2*w2r*rhogeq*w3i*k*Kdrag*vdsol + w2i**2*w2r*rhogeq**2*w3i**2*k*vgsol -&
        w2i**2*w2r*rhogeq*w1r*w3r**2*rhogsol - w2i**2*w2r*rhogeq*w3i**2*w1r*rhogsol -&
        w2i**2*w2r*rhogeq*w3r*w1r**2*rhogsol + w2i**2*w2r*Kdrag**2*k*vgsol + w2i**2*w2r*rhogeq**2*w3r**2*k*vgsol -&
        2*w2i**2*rhogeq**2*k*vgsol*w2r**2*w1r + 2*w2i**2*rhogeq*w2r**2*k**2*cs**2*rhogsol +&
        2*w2i**2*rhogeq*w3r*rhogsol*w2r**2*w1r + w2i**2*vdsol*k*Kdrag**2*w1r - w2i**2*w2r*Kdrag**2*k*vdsol -&
        w2i**2*vgsol*k*Kdrag**2*w3r - w2i**2*vgsol*k*Kdrag**2*w1r + w2i**2*vdsol*k*Kdrag**2*w3r +&
        w2i**2*vgsol*k*rhogeq**2*w1r**2*w3r - w2i**2*rhogsol*cs**2*k**2*Kdrag*w3i + w2i**3*w1i*rhogeq*rhogsol*w3i**2 +&
        rhogsol*Kdrag*w2r**2*w1r**2*w3i - w3r*w2i*Kdrag**2*w1i*k*vdsol + w3r*w2i*Kdrag**2*w1i*k*vgsol +&
        w3r*Kdrag**2*w2r*w1r*k*vgsol + rhogsol*cs**4*k**4*rhogeq*w2r*w1r - w3r*w2i*Kdrag*k**2*cs**2*w1r*rhogsol +&
        w3r*w2i*Kdrag*w1i**2*rhogeq*k*vgsol - w3r*Kdrag*w1i*w2r**2*rhogeq*k*vgsol + vgsol*Kdrag*k*rhogeq*w3i*w1r*w2i**2&
        + Kdrag*w3i**2*w1i*rhogsol*w2i**2 + w3r*w2i*Kdrag*rhogeq*w1r**2*k*vgsol + vdsol*k*Kdrag**2*w2r**2*w1r +&
        w1i**2*rhogeq**2*k*vgsol*w2r**3 - w1i**2*rhogeq*w3r*rhogsol*w2r**3 - w1i**2*rhogsol*cs**2*k**2*rhogeq*w3r**2 +&
        w3r*Kdrag*w1i*w2r*k**2*cs**2*rhogsol + w3r*rhogeq*w2i**2*w1i*k*Kdrag*vgsol +&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r + rhogsol*cs**2*k**2*rhogeq*w2i*w3i*w1r**2 +&
        rhogsol*Kdrag*w3i*w2i**2*w1r**2 - vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r +&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i + w1i**2*rhogsol*Kdrag*w3i*w2i**2 + rhogsol*Kdrag*w1i**2*w3i*w2r**2 -&
        vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r + vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r -&
        vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r + rhogeq**2*w3i**2*w2i**2*k*w1r*vgsol -&
        vgsol*Kdrag*k*rhogeq*w3i*w1r*w2r**2 - w1i*rhogsol*cs**4*k**4*rhogeq*w2i + w3i*rhogeq*k**4*cs**4*w1i*rhogsol -&
        Kdrag*w3r**2*w1i**2*rhogsol*w2i - w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol - w3i*rhogeq*k**3*cs**2*Kdrag*vdsol*w1r&
        - w3i**2*rhogeq**2*k**3*cs**2*w2r*vgsol + w3i**2*rhogeq**2*k**3*cs**2*w1r*vgsol +&
        w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol - rhogeq*w3i**2*w1i*w2r*k*Kdrag*vgsol +&
        rhogeq*w3i**2*w1i**2*w2r**2*rhogsol - rhogeq*w3i**2*w1i**2*rhogsol*w2i**2 +&
        w3i*rhogeq*k**3*cs**2*Kdrag*vgsol*w1r + rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol + Kdrag*w3r**2*w1i*w2r**2*rhogsol&
        - Kdrag*w3r**2*w2i*w1r**2*rhogsol + Kdrag*w3r**2*w1i*rhogsol*w2i**2 - rhogeq**2*w3i**2*w2r**2*k*w1r*vgsol +&
        rhogeq*w3i**2*w2i*Kdrag*w1r*k*vgsol + Kdrag*w3i*w1i*w2i*k**2*cs**2*rhogsol - Kdrag**2*w3i*w2i*w1r*k*vdsol -&
        Kdrag*w3i**2*w2i*w1r**2*rhogsol - Kdrag*w3i**2*w1i**2*rhogsol*w2i + Kdrag**2*w3i*w2i*w1r*k*vgsol +&
        Kdrag*w3i**2*w1i*w2r**2*rhogsol - Kdrag*w3i*rhogeq*w2r*w1r**2*k*vgsol - Kdrag*w3i*w1i**2*rhogeq*w2r*k*vgsol -&
        rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol - w3i*rhogeq*k**4*cs**4*w2i*rhogsol +&
        Kdrag*w3i*w2r*w1r*k**2*cs**2*rhogsol - rhogeq*w3i**2*w2i**2*w1r**2*rhogsol +&
        rhogeq*w3i**2*w2r**2*w1r**2*rhogsol + rhogeq*w3r**2*w1r*w2i*k*Kdrag*vgsol - Kdrag**2*w3i*w1i*w2r*k*vgsol +&
        rhogeq*w3r*w2r*w1r**2*k**2*cs**2*rhogsol + rhogeq*w3r*cs**2*k**3*w1i*Kdrag*vgsol -&
        rhogeq**2*w3r*cs**2*k**3*w2i**2*vgsol + rhogeq**2*w3r*cs**2*k**3*w1r**2*vgsol +&
        rhogeq**2*w3r**2*cs**2*k**3*w1r*vgsol - rhogeq*w3r*cs**2*k**3*w1i*Kdrag*vdsol -&
        2*rhogeq**2*w3r*cs**2*k**3*vgsol*w1r*w2r + rhogeq*w3r*cs**2*k**3*w2i*Kdrag*vdsol -&
        rhogeq**2*w3r**2*cs**2*k**3*w2r*vgsol + Kdrag**2*w3i*w1i*w2r*k*vdsol - rhogeq*w3r*cs**4*k**4*w1r*rhogsol +&
        rhogeq*w3r**2*w2r**2*w1r**2*rhogsol + rhogeq*w3r*cs**4*k**4*w2r*rhogsol - rhogeq*w3r*cs**2*k**3*w2i*Kdrag*vgsol&
        + rhogeq**2*w3r**2*w1r*w2i**2*k*vgsol + rhogeq**2*w3r*cs**2*k**3*w1i**2*vgsol +&
        rhogeq**2*w3r*cs**2*k**3*w2r**2*vgsol + rhogeq*w3r**2*w2r**2*w1i**2*rhogsol -&
        rhogeq*w3r**2*w2r*w1i*k*Kdrag*vgsol - rhogeq**2*w3r**2*w2r**2*k*vgsol*w1r - rhogeq*w3r**2*w1r**2*w2i**2*rhogsol&
        - rhogeq*w3r**2*w2i**2*w1i**2*rhogsol + rhogeq*w3r*w1i**2*w2r*k**2*cs**2*rhogsol -&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r - w3r*Kdrag**2*w2r*w1r*k*vdsol)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 +&
        w2r**2 - 2*w2r*w1r)/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/rhogeq/Kdrag/k 

  vd2i = - 1/k*( - w2i**3*Kdrag**2*k*vdsol + rhogeq**2*w2i**2*w1i**2*w3i*k*vgsol +&
        2*rhogeq**2*w2i*w1i**2*w2r*w3r*k*vgsol + rhogeq**2*w2i**2*w3i*w1r**2*k*vgsol +&
        rhogeq*w2i**2*w1r*w3r*k*Kdrag*vdsol - 2*rhogeq*w2i**2*k**2*cs**2*w3i*w1r*rhogsol -&
        3*rhogeq*w2i**2*w3r*Kdrag*w1r*k*vgsol - rhogeq**2*w2i**3*w3i**2*k*vgsol - 2*rhogeq**2*w2i**3*w3r*k*w1r*vgsol -&
        rhogeq**2*w2i**3*w3r**2*k*vgsol - rhogeq**2*w2i*w2r**2*w3i**2*k*vgsol - 2*rhogeq**2*w2i*w3r*w2r**2*k*w1r*vgsol&
        - Kdrag*w2r*w1r*w3r*k**2*cs**2*rhogsol + Kdrag*w1r*w3r**2*w2r**2*rhogsol + Kdrag*w3r*w2r**2*w1r**2*rhogsol +&
        Kdrag*w3i**2*w1r*w2r**2*rhogsol + w2i**2*Kdrag*w3i**2*w1r*rhogsol + w2i**2*Kdrag*w3r*w1r**2*rhogsol +&
        w2i**2*Kdrag*w1r*w3r**2*rhogsol + w2i**2*Kdrag*w1i**2*w3r*rhogsol - Kdrag*w1i**2*w3i**2*w2r*rhogsol +&
        Kdrag*w1i**2*w2r**2*w3r*rhogsol - Kdrag*w1i**2*w3r**2*w2r*rhogsol - Kdrag**2*w3i*w2r*w1r*k*vdsol -&
        Kdrag*w1r**2*w3r**2*w2r*rhogsol - w2i**3*Kdrag*w3i*w1r*rhogsol - rhogeq**2*w2i**3*w1i**2*k*vgsol -&
        rhogeq**2*w2i*w1i**2*w2r**2*k*vgsol - 2*rhogeq*w2i**2*w1i*w2r**2*w3r*rhogsol + rhogeq*w2i**3*w1i*k*Kdrag*vdsol&
        + 2*rhogeq**2*w2i**2*w1i*w2r**2*k*vgsol + rhogeq**2*w2i**4*w1i*k*vgsol + rhogeq*w2i*w1i*w2r**2*k*Kdrag*vdsol -&
        rhogeq*w2i**4*w1i*w3r*rhogsol + 2*rhogeq*w2i*w3i*w2r*w1r*k*Kdrag*vgsol -&
        4*rhogeq*w2i*w2r*w1r*w3r*k**2*cs**2*rhogsol + rhogeq*w2i*w1r*w3r**2*w2r**2*rhogsol +&
        rhogeq*w2i*w3r*w2r**2*w1r**2*rhogsol + rhogeq*w2i*w3i**2*w1r*w2r**2*rhogsol + rhogeq*w2i**3*w3i**2*w1r*rhogsol&
        + rhogeq*w2i**3*w3r*w1r**2*rhogsol + rhogeq*w2i**3*w1r*w3r**2*rhogsol + rhogeq*w2i**3*w1i**2*w3r*rhogsol -&
        2*rhogeq*w2i*w1i**2*w3i**2*w2r*rhogsol + rhogeq*w2i*w1i**2*w2r**2*w3r*rhogsol -&
        2*rhogeq*w2i*w1i**2*w3r**2*w2r*rhogsol - 2*rhogeq*w2i*w3i*w2r*w1r*k*Kdrag*vdsol -&
        2*rhogeq*w2i*w1r**2*w3r**2*w2r*rhogsol - rhogeq*w2i**4*k*Kdrag*vdsol + 2*rhogeq*w2i**2*w2r**2*k*Kdrag*vgsol +&
        rhogeq*w2i**3*w3i*k*Kdrag*vdsol + rhogeq**2*w2i**4*w3i*k*vgsol - rhogeq**2*w2i**3*w1r**2*k*vgsol +&
        rhogeq*w2i*w2r**2*w3r*k**2*cs**2*rhogsol + rhogeq*w2i**3*k**2*cs**2*w1r*rhogsol +&
        rhogeq*w2i*k*w3i*Kdrag*vdsol*w2r**2 - 2*rhogeq*w2i**2*w2r**2*k*Kdrag*vdsol -&
        2*rhogeq*w2i*w3i**2*w1r**2*w2r*rhogsol - 2*cs**2*k**3*rhogeq**2*w1i*w2i*w3i*vgsol +&
        cs**2*k**3*rhogeq**2*w1i*w3i**2*vgsol + cs**2*k**3*rhogeq**2*w1i*w3r**2*vgsol -&
        cs**2*k**3*rhogeq*w1i*w3i*Kdrag*vdsol - cs**4*k**4*rhogeq*w1i*w3r*rhogsol +&
        cs**2*k**3*rhogeq*w1i*w3i*Kdrag*vgsol - Kdrag*w2r**3*k**2*cs**2*rhogsol + w2i*Kdrag**2*w2r**2*k*vgsol +&
        w2i**2*Kdrag**2*w3i*k*vdsol - w2i**2*Kdrag*w2r*k**2*cs**2*rhogsol - w2i**2*Kdrag*rhogeq*w1r**2*k*vgsol +&
        Kdrag*w2r**2*w3r*k**2*cs**2*rhogsol + w2i**2*Kdrag*k**2*cs**2*w1r*rhogsol - w2i*Kdrag**2*w2r**2*k*vdsol +&
        Kdrag**2*k*w3i*vdsol*w2r**2 - w2i*Kdrag*w3i*w1r*rhogsol*w2r**2 - Kdrag**2*k*w3i*vgsol*w2r**2 -&
        w2i**2*Kdrag*w2r*w1r*w3r*rhogsol - w2i**2*Kdrag**2*w3i*k*vgsol + w2i**2*Kdrag*w3r*k**2*cs**2*rhogsol +&
        Kdrag*w2r**2*w1r*k**2*cs**2*rhogsol - Kdrag*w2r**2*rhogeq*w1r**2*k*vgsol - w2i**2*Kdrag*w1i**2*rhogeq*k*vgsol -&
        Kdrag*w1i**2*w2r**2*rhogeq*k*vgsol - w2i*Kdrag*w1i*w2r**2*w3r*rhogsol + w2i**2*Kdrag**2*w1i*k*vdsol -&
        w2i**2*Kdrag**2*w1i*k*vgsol + w2i**2*Kdrag*w1i*w2r*w3i*rhogsol + Kdrag**2*w1i*w2r**2*k*vdsol -&
        Kdrag**2*w1i*w2r**2*k*vgsol + Kdrag*w1i*w2r**3*w3i*rhogsol - w2i**3*Kdrag*w1i*w3r*rhogsol +&
        cs**2*k**3*rhogeq*w2r**2*Kdrag*vdsol - cs**2*k**3*rhogeq**2*w3i*w2r**2*vgsol +&
        2*cs**2*k**2*rhogeq*w3i*w1r*rhogsol*w2r**2 - cs**2*k**3*rhogeq*w2i*w3i*Kdrag*vgsol +&
        2*cs**2*k**3*rhogeq**2*w2i*w1r*w2r*vgsol + cs**4*k**4*rhogeq*w2i*w3r*rhogsol +&
        cs**2*k**3*rhogeq*w2i**2*Kdrag*vgsol - cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol +&
        2*cs**2*k**2*rhogeq*w1i*w2r**2*w3r*rhogsol + cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol -&
        cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol - cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol +&
        cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol + cs**2*k**2*rhogeq*w2i*w3i**2*w1r*rhogsol +&
        cs**2*k**2*rhogeq*w2i*w3r*w1r**2*rhogsol + cs**2*k**2*rhogeq*w2i*w1r*w3r**2*rhogsol +&
        cs**2*k**2*rhogeq*w1i**2*w2i*w3r*rhogsol + cs**2*k**3*rhogeq**2*w1i**2*w3i*vgsol +&
        cs**2*k**3*rhogeq**2*w3i*w1r**2*vgsol + cs**2*k**3*rhogeq*w1r*w3r*Kdrag*vdsol -&
        cs**4*k**4*rhogeq*w3i*w1r*rhogsol - cs**2*k**3*rhogeq*w3r*Kdrag*w1r*vgsol -&
        cs**2*k**3*rhogeq**2*w2i*w3i**2*vgsol - 2*cs**2*k**3*rhogeq**2*w2i*w3r*w1r*vgsol -&
        cs**2*k**3*rhogeq**2*w2i*w3r**2*vgsol + 2*rhogeq**2*w2i**2*w3i*w2r**2*k*vgsol -&
        2*rhogeq*w2i**2*w3i*w1r*rhogsol*w2r**2 + rhogeq*w2i**3*w3r*k**2*cs**2*rhogsol + rhogeq*w2i**4*k*Kdrag*vgsol +&
        rhogeq*w2i*w2r**2*w1r*k**2*cs**2*rhogsol + 2*rhogeq**2*w2i*w3r*w2r*w1r**2*k*vgsol +&
        2*rhogeq**2*w2i*w3i**2*k*vgsol*w1r*w2r + 2*rhogeq**2*w2i*w3r**2*k*vgsol*w1r*w2r -&
        rhogeq**2*w2i*w2r**2*w3r**2*k*vgsol - 2*rhogeq**2*w2i*w1i*w2r**2*w3i*k*vgsol -&
        2*rhogeq**2*w2i**3*w1i*w3i*k*vgsol + rhogeq**2*w2i**2*w1i*w3i**2*k*vgsol + rhogeq**2*w2i**2*w1i*w3r**2*k*vgsol&
        - rhogeq*w2i**2*w1i*k*w3i*Kdrag*vdsol - 2*rhogeq*w2i**2*w1i*w3r*k**2*cs**2*rhogsol +&
        4*rhogeq*w2i*w1i*w3i*w2r*k**2*cs**2*rhogsol - 2*rhogeq*w2i*w1i*w2r*w3r*k*Kdrag*vdsol +&
        2*rhogeq*w2i*w1i*w2r*w3r*k*Kdrag*vgsol - rhogeq*w2i**2*w1i*k*w3i*Kdrag*vgsol - rhogeq*w2i**4*w3i*w1r*rhogsol -&
        cs**2*k**3*rhogeq*w2i**2*Kdrag*vdsol - cs**2*k**3*rhogeq*w2r**2*Kdrag*vgsol +&
        cs**2*k**3*rhogeq*w2i*w3i*Kdrag*vdsol - 2*cs**4*k**4*rhogeq*w2i*w2r*rhogsol +&
        cs**2*k**3*rhogeq**2*w2i**2*w3i*vgsol - cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol +&
        2*cs**2*k**3*rhogeq**2*w2i*w2r*w3r*vgsol + w2i**3*Kdrag**2*k*vgsol + cs**4*k**4*rhogeq*w2i*w1r*rhogsol -&
        rhogeq**2*w2i*w2r**2*w1r**2*k*vgsol + Kdrag**2*w3i*w2r*w1r*k*vgsol - Kdrag*w3i**2*w1r**2*w2r*rhogsol +&
        w2i*Kdrag*w1i*rhogeq*w3r**2*k*vgsol - w2i*Kdrag**2*w1i*k*w3i*vdsol - w2i*Kdrag*w1i*w3r*k**2*cs**2*rhogsol +&
        Kdrag*w1i*w3i*w2r*k**2*cs**2*rhogsol - Kdrag**2*w1i*w2r*w3r*k*vdsol + Kdrag**2*w1i*w2r*w3r*k*vgsol +&
        w2i*Kdrag**2*w1i*k*w3i*vgsol - Kdrag*w2r**3*w1r*w3r*rhogsol + w2r**4*rhogeq*k*Kdrag*vgsol -&
        w2r**4*rhogeq*k*Kdrag*vdsol + w2r**4*rhogeq**2*w3i*k*vgsol - w2r**4*rhogeq*w3i*w1r*rhogsol -&
        w2r**4*rhogeq*w1i*w3r*rhogsol + w2r**4*rhogeq**2*w1i*k*vgsol - w2r**2*rhogeq**2*w1i**2*w3i*k*vgsol -&
        w2r**2*rhogeq**2*w3i*w1r**2*k*vgsol - w2r**2*rhogeq*w1r*w3r*k*Kdrag*vdsol - w2r**2*rhogeq**2*w1i*w3i**2*k*vgsol&
        - w2r**2*rhogeq**2*w1i*w3r**2*k*vgsol + w2r**2*rhogeq*w1i*k*w3i*Kdrag*vdsol -&
        rhogeq*w2r*k**2*cs**2*rhogsol*w3i*w1r**2 - rhogeq*w2r*w1i*w3r**2*k**2*cs**2*rhogsol -&
        rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w3i**2 - rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i +&
        rhogeq*w2r*w3i*rhogsol*w2i**2*w1r**2 + rhogeq*w2r**3*w3i*w1r**2*rhogsol -&
        rhogeq*w2r*w2i**2*w3i*k**2*cs**2*rhogsol - rhogeq*w2r**3*w3i*k**2*cs**2*rhogsol +&
        rhogeq*w2r*w2i**2*w1r*k*Kdrag*vdsol + rhogeq*w2r**3*w1r*k*Kdrag*vdsol + rhogeq*w2r*w2i**2*w3r*k*Kdrag*vdsol +&
        rhogeq*w2r**3*w3r*k*Kdrag*vdsol - rhogeq*w2r*w2i**2*w1i*k**2*cs**2*rhogsol +&
        rhogeq*w2r*w2i**2*w1i*rhogsol*w3i**2 + rhogeq*w2r*w2i**2*w1i*rhogsol*w3r**2 +&
        rhogeq*w2r*w2i**2*w1i**2*w3i*rhogsol + rhogeq*w2r*cs**2*k**3*w3r*Kdrag*vgsol +&
        rhogeq*w2r*cs**4*k**4*w1i*rhogsol + rhogeq*w2r**3*w1i*rhogsol*w3i**2 + rhogeq*w2r**3*w1i*rhogsol*w3r**2 +&
        rhogeq*w2r**3*w1i**2*w3i*rhogsol - rhogeq*w2r**3*w1i*k**2*cs**2*rhogsol + rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol&
        + rhogeq*w2r*cs**4*k**4*w3i*rhogsol - rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol -&
        rhogeq*w2r*cs**2*k**3*w3r*Kdrag*vdsol + w2i*Kdrag*w1i**2*w3i*rhogeq*k*vgsol +&
        Kdrag*w1i**2*rhogeq*w2r*w3r*k*vgsol + w2i*Kdrag*rhogeq*w3i*w1r**2*k*vgsol + w2i*Kdrag**2*w1r*w3r*k*vdsol -&
        w2i*Kdrag*k**2*cs**2*w3i*w1r*rhogsol - w2i*Kdrag**2*w3r*w1r*k*vgsol - w2i**2*Kdrag*rhogeq*w3i**2*k*vgsol -&
        w2i**2*Kdrag*rhogeq*w3r**2*k*vgsol - Kdrag*rhogeq*w2r**2*w3i**2*k*vgsol - Kdrag*rhogeq*w3r*w2r**2*k*w1r*vgsol +&
        Kdrag*rhogeq*w3r*w2r*w1r**2*k*vgsol + Kdrag*rhogeq*w3i**2*k*vgsol*w1r*w2r + Kdrag*rhogeq*w3r**2*k*vgsol*w1r*w2r&
        - Kdrag*rhogeq*w2r**2*w3r**2*k*vgsol - 3*Kdrag*w1i*w2r**2*w3i*rhogeq*k*vgsol +&
        w2i*Kdrag*w1i*rhogeq*w3i**2*k*vgsol)/Kdrag/rhogeq/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 -&
        2*w3i*w2i)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  vd1r = - (rhogsol*cs**2*k**2*rhogeq*w3r**2*w2i**2 + rhogsol*cs**2*k**2*rhogeq*w3r**2*w2r**2 +&
        rhogsol*cs**2*k**2*rhogeq*w3i**2*w2r**2 + rhogsol*cs**2*k**2*rhogeq*w3i**2*w2i**2 +&
        Kdrag*w3i*w1i**2*k**2*cs**2*rhogsol - 2*vgsol*k*rhogeq**2*w1r*w1i**2*w2r*w3r -&
        2*vgsol*k*rhogeq**2*w1r**3*w2r*w3r - w2r*cs**2*k**2*rhogeq*w1r*w3r**2*rhogsol -&
        w2r*cs**2*k**2*rhogeq*w3i**2*w1r*rhogsol - w2i*w1i*rhogsol*cs**2*k**2*rhogeq*w3i**2 -&
        w2i*w1i*rhogsol*cs**2*k**2*rhogeq*w3r**2 + w2r**2*vgsol*k*rhogeq**2*w1r**2*w3r -&
        w2r**2*w1i**2*vgsol*k*rhogeq**2*w3r - w2i**2*w1i**2*vgsol*k*rhogeq**2*w3r + w2i**2*vgsol*k*rhogeq**2*w1r**2*w3r&
        - rhogsol*Kdrag*w2r**2*w1r**2*w3i + w3r*w2i*Kdrag**2*w1i*k*vdsol - w3r*w2i*Kdrag**2*w1i*k*vgsol -&
        w3r*Kdrag**2*w2r*w1r*k*vgsol - rhogsol*cs**4*k**4*rhogeq*w2r*w1r - w3r*w2i*Kdrag*k**2*cs**2*w1r*rhogsol -&
        w3r*w2i*Kdrag*w1i**2*rhogeq*k*vgsol - w3r*Kdrag*w1i*w2r**2*rhogeq*k*vgsol + vgsol*Kdrag*k*rhogeq*w3i*w1r*w2i**2&
        + Kdrag*w3i**2*w1i*rhogsol*w2i**2 + w3r*w2i*Kdrag*rhogeq*w1r**2*k*vgsol + w3r*Kdrag*w1i*w2r*k**2*cs**2*rhogsol&
        - w3r*rhogeq*w2i**2*w1i*k*Kdrag*vgsol - vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r - rhogsol*Kdrag*w3i*w2i**2*w1r**2&
        - w1i*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2 + vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r -&
        rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**2 - w1i**2*rhogsol*Kdrag*w3i*w2i**2 - rhogsol*Kdrag*w1i**2*w3i*w2r**2 +&
        vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r - vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r -&
        vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r + vgsol*Kdrag*k*rhogeq*w3i*w1r*w2r**2 + w1i*rhogsol*cs**4*k**4*rhogeq*w2i&
        + w3i*rhogeq*k**4*cs**4*w1i*rhogsol - Kdrag*w3r**2*w1i**2*rhogsol*w2i - w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol -&
        w3i*rhogeq*k**3*cs**2*Kdrag*vdsol*w1r - w3i**2*rhogeq**2*k**3*cs**2*w2r*vgsol +&
        w3i**2*rhogeq**2*k**3*cs**2*w1r*vgsol + w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol -&
        rhogeq*w3i**2*w1i*w2r*k*Kdrag*vgsol + rhogeq*w3i**2*w1i**2*w2r**2*rhogsol + rhogeq*w3i**2*w1i**2*rhogsol*w2i**2&
        - rhogeq**2*w3i**2*w1i**2*w2r*k*vgsol + w3i*rhogeq*k**3*cs**2*Kdrag*vgsol*w1r +&
        rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol + Kdrag*w3r**2*w1i*w2r**2*rhogsol - Kdrag*w3r**2*w2i*w1r**2*rhogsol +&
        Kdrag*w3r**2*w1i*rhogsol*w2i**2 + rhogeq**2*w3i**2*w2r*w1r**2*k*vgsol + rhogeq*w3i**2*w2i*Kdrag*w1r*k*vgsol -&
        Kdrag*w3i*w1i*w2i*k**2*cs**2*rhogsol - Kdrag**2*w3i*w2i*w1r*k*vdsol - Kdrag*w3i**2*w2i*w1r**2*rhogsol -&
        Kdrag*w3i**2*w1i**2*rhogsol*w2i + Kdrag**2*w3i*w2i*w1r*k*vgsol + Kdrag*w3i**2*w1i*w2r**2*rhogsol +&
        Kdrag*w3i*rhogeq*w2r*w1r**2*k*vgsol - Kdrag*w3i*w1i**2*rhogeq*w2r*k*vgsol +&
        rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol - w3i*rhogeq*k**4*cs**4*w2i*rhogsol -&
        rhogeq*w3r*w2r**2*w1r*k**2*cs**2*rhogsol - Kdrag*w3i*w2r*w1r*k**2*cs**2*rhogsol -&
        rhogeq*w3i**2*w2i**2*w1r**2*rhogsol - rhogeq*w3i**2*w2r**2*w1r**2*rhogsol + rhogeq*w3r**2*w1r*w2i*k*Kdrag*vgsol&
        - Kdrag**2*w3i*w1i*w2r*k*vgsol + rhogeq*w3r*cs**2*k**3*w1i*Kdrag*vgsol - rhogeq**2*w3r*cs**2*k**3*w2i**2*vgsol&
        - rhogeq**2*w3r*cs**2*k**3*w1r**2*vgsol + rhogeq**2*w3r**2*cs**2*k**3*w1r*vgsol -&
        rhogeq*w3r*cs**2*k**3*w1i*Kdrag*vdsol + 2*rhogeq**2*w3r*cs**2*k**3*vgsol*w1r*w2r +&
        rhogeq*w3r*cs**2*k**3*w2i*Kdrag*vdsol - rhogeq**2*w3r**2*cs**2*k**3*w2r*vgsol + Kdrag**2*w3i*w1i*w2r*k*vdsol -&
        rhogeq*w3r*cs**4*k**4*w1r*rhogsol - rhogeq*w3r**2*w2r**2*w1r**2*rhogsol + rhogeq*w3r*cs**4*k**4*w2r*rhogsol -&
        rhogeq*w3r*cs**2*k**3*w2i*Kdrag*vgsol + rhogeq**2*w3r*cs**2*k**3*w1i**2*vgsol -&
        rhogeq**2*w3r*cs**2*k**3*w2r**2*vgsol - rhogeq*w3r*w2i**2*k**2*cs**2*w1r*rhogsol +&
        rhogeq*w3r**2*w2r**2*w1i**2*rhogsol - rhogeq*w3r**2*w2r*w1i*k*Kdrag*vgsol - rhogeq*w3r**2*w1r**2*w2i**2*rhogsol&
        - rhogeq**2*w3r**2*w2r*w1i**2*k*vgsol + rhogeq*w3r**2*w2i**2*w1i**2*rhogsol +&
        rhogeq**2*w3r**2*w2r*w1r**2*k*vgsol + vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r + w3r*Kdrag**2*w2r*w1r*k*vdsol -&
        2*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r**2 + w1i**3*rhogsol*cs**2*k**2*rhogeq*w3i -&
        rhogeq**2*w3i**2*k*w1r**3*vgsol - w1i**2*rhogsol*cs**4*k**4*rhogeq + 2*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i +&
        rhogsol*cs**2*k**2*rhogeq*w3i*w1r**2*w1i - rhogsol*cs**2*k**2*rhogeq*w1r**4 +&
        2*w3r*Kdrag*rhogeq*w1r**2*k*vgsol*w1i + 2*w3r*Kdrag*w1i**3*rhogeq*k*vgsol - 2*vgsol*Kdrag*k*rhogeq*w3i*w1r**3 +&
        w3r*Kdrag**2*w1i**2*k*vgsol - w3r*Kdrag**2*w1i**2*k*vdsol + vdsol*k*Kdrag**2*w1r**3 +&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w2r*w1r - vgsol*k*Kdrag**2*w1r**3 - rhogeq**2*w3i**2*k*w1r*vgsol*w1i**2 +&
        rhogsol*cs**2*k**2*rhogeq*w1r**3*w2r + w1i**2*rhogsol*Kdrag*w3i*w2r*w1r + rhogsol*Kdrag*w3i*w1r**3*w2r -&
        w1i**4*rhogsol*cs**2*k**2*rhogeq - 2*w3r*rhogeq*w1i*k*Kdrag*vgsol*w2r*w1r -&
        2*vgsol*Kdrag*k*rhogeq*w3i*w1r*w1i**2 + w2i*rhogeq*rhogsol*w3i*w1r**4 - w2i*rhogeq*w3i**2*w1i**3*rhogsol +&
        rhogeq*w3i**2*w1i**2*rhogsol*w2r*w1r - 2*w2i*vgsol*Kdrag*k*rhogeq*w1r*w1i**2 +&
        w2i*Kdrag*k**2*cs**2*rhogsol*w1i**2 - 2*w2i*vgsol*Kdrag*k*rhogeq*w1r**3 +&
        2*w2i*w1i**2*rhogeq*rhogsol*w3i*w1r**2 + w2i*vdsol*Kdrag*k*rhogeq*w1r*w1i**2 + w2i*rhogsol*Kdrag*w3r*w1r*w1i**2&
        + w2i*vdsol*Kdrag*k*rhogeq*w1r**3 + w2i*w1i**4*rhogeq*rhogsol*w3i + w2i*vdsol*Kdrag*k*rhogeq*w3r*w1i**2 -&
        rhogeq**2*w3r**2*w1r**3*k*vgsol + 2*w2i*vgsol*Kdrag*k*rhogeq*w3i*w1r*w1i +&
        w2i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w1i + w2i*rhogsol*Kdrag*w3i*w1r**2*w1i +&
        2*w2i*rhogeq**2*w3i**2*k*w1r*vgsol*w1i - 2*w2i*w3i*vgsol*k*rhogeq**2*w1r*w1i**2 +&
        w2i*w1i**3*rhogsol*cs**2*k**2*rhogeq + 2*w2i*w3i*rhogeq**2*k**3*cs**2*w1r*vgsol -&
        2*w2i*w3i*vgsol*k*rhogeq**2*w1r**3 - 2*w2i*vdsol*Kdrag*k*rhogeq*w1r*w1i*w3i - w2i*rhogeq*w3r**2*w1i**3*rhogsol&
        + w2i*Kdrag*k**2*cs**2*rhogsol*w1r**2 - w2i*vdsol*Kdrag*k*rhogeq*w1r**2*w3r - w2i**2*vgsol*k*rhogeq**2*w1r**3 +&
        w2i**2*rhogeq*rhogsol*w3r*w1r**3 - w2i**2*w1i**3*rhogeq*rhogsol*w3i - w2i**2*vgsol*k*rhogeq**2*w1r*w1i**2 -&
        w2i**2*w1i*rhogeq*rhogsol*w3i*w1r**2 + 2*w2i**2*vgsol*k*rhogeq**2*w1r*w1i*w3i +&
        w2i**2*rhogeq*rhogsol*w3r*w1r*w1i**2 - w2i*rhogeq*w3r**2*w1r**2*rhogsol*w1i -&
        2*w2i*rhogeq**2*k**3*cs**2*w1r*vgsol*w1i + 2*w2i*rhogeq**2*w3r**2*w1r*k*vgsol*w1i -&
        w2i*rhogeq*w3i**2*w1r**2*rhogsol*w1i + w2i*w1i**3*rhogsol*Kdrag*w3i + w2i*rhogsol*Kdrag*w3r*w1r**3 +&
        rhogeq*w3i**2*w1r**3*rhogsol*w2r - 2*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i -&
        rhogeq**2*w3r**2*w1r*k*vgsol*w1i**2 + rhogeq*w3r*k**2*cs**2*w1r*rhogsol*w1i**2 +&
        rhogeq*w3r**2*w1r**3*rhogsol*w2r + rhogeq*w3r**2*w1i**2*rhogsol*w2r*w1r + rhogeq*w3r*k**2*cs**2*w1r**3*rhogsol&
        + 2*w1i*vgsol*Kdrag*k*rhogeq*w2r*w1r**2 + w1i**4*w2r*rhogeq**2*k*vgsol - w1i**3*w2r*Kdrag*w3r*rhogsol +&
        w1i**4*vgsol*k*rhogeq**2*w3r - w2r*rhogeq*w3i*k*Kdrag*vdsol*w1r**2 + 2*w1i**3*vgsol*Kdrag*k*rhogeq*w2r -&
        2*w1i**2*w2r*rhogeq*w3r*rhogsol*w1r**2 - w1i**4*w2r*rhogeq*w3r*rhogsol - w1i*vdsol*Kdrag*k*rhogeq*w2r*w1r**2 -&
        2*rhogeq**2*k**3*cs**2*w1r*vgsol*w1i*w3i - w1i**3*vdsol*Kdrag*k*rhogeq*w2r -&
        w1i*vdsol*Kdrag*k*rhogeq*w3r*w1r**2 - w1i**3*vdsol*Kdrag*k*rhogeq*w3r + cs**4*k**4*rhogeq*rhogsol*w1r**2 +&
        2*w1i*vdsol*Kdrag*k*rhogeq*w2r*w1r*w3r - w1i**3*rhogsol*cs**2*k**2*Kdrag +&
        2*w1i**2*vgsol*k*rhogeq**2*w3r*w1r**2 - w1i*rhogsol*cs**2*k**2*Kdrag*w1r**2 +&
        2*w1i**2*w2r*rhogeq**2*k*vgsol*w1r**2 - w1i*w2r*Kdrag*w3r*rhogsol*w1r**2 + w2r*rhogeq**2*w1r**4*k*vgsol +&
        vgsol*k*rhogeq**2*w1r**4*w3r + vdsol*Kdrag*k*rhogeq*w1r**3*w3i - w1i*rhogeq*rhogsol*w3i*w2r**2*w1r**2 +&
        vdsol*Kdrag*k*rhogeq*w1r*w1i**2*w3i - vdsol*k*Kdrag**2*w1r**2*w3r + rhogeq*rhogsol*w3r*w1r**3*w2r**2 +&
        vgsol*k*Kdrag**2*w1r**2*w3r - w1i**3*rhogeq*rhogsol*w3i*w2r**2 + 2*vgsol*k*rhogeq**2*w1r*w1i*w3i*w2r**2 -&
        vgsol*k*rhogeq**2*w1r**3*w2r**2 + Kdrag*k**2*cs**2*rhogsol*w3i*w1r**2 - vgsol*k*rhogeq**2*w1r*w1i**2*w2r**2 -&
        vgsol*k*Kdrag**2*w1r*w1i**2 + w2r*Kdrag**2*k*vgsol*w1r**2 + w1i**2*vdsol*Kdrag*k*rhogeq*w2r*w3i +&
        rhogeq*rhogsol*w3r*w1r*w1i**2*w2r**2 - w2r*rhogeq*w3r*w1r**4*rhogsol - w2r*Kdrag**2*k*vdsol*w1r**2 +&
        vdsol*k*Kdrag**2*w1r*w1i**2 - w2r*Kdrag**2*k*vdsol*w1i**2 + w2r*Kdrag**2*k*vgsol*w1i**2)/( - 2*w1i*w3i + w1i**2&
        + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/k/Kdrag/rhogeq/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 -&
        2*w2r*w1r) 

  vd1i =(rhogeq**2*w3i**2*w2i*w1r**2*k*vgsol - w3r**2*rhogeq**2*w1i**2*w2i*k*vgsol -&
        rhogeq**2*w3i**2*w1i**2*w2i*k*vgsol + w3r**2*rhogeq**2*w2i*w1r**2*k*vgsol - rhogeq**2*w1r**4*k*vgsol*w3i +&
        rhogeq*k**2*cs**2*w1r**3*rhogsol*w3i + Kdrag*w2r*w1r*w3r*rhogsol*w1i**2 - Kdrag*k**2*cs**2*w1r**2*rhogsol*w3r +&
        2*cs**4*k**4*rhogeq*w1r*rhogsol*w1i + Kdrag**2*w3i*k*vgsol*w1r**2 - Kdrag*w3r*k**2*cs**2*rhogsol*w1i**2 -&
        Kdrag*w1i**4*rhogeq*k*vgsol + Kdrag*k**2*cs**2*w1r*rhogsol*w1i**2 + vdsol*k*Kdrag*rhogeq*w2r*w1r**2*w3r +&
        2*vdsol*k*Kdrag*rhogeq*w1i**2*w1r**2 - rhogeq*w2r*w1i**3*w3i**2*rhogsol - 2*rhogeq**2*w1i**2*k*vgsol*w3i*w1r**2&
        + vdsol*k*Kdrag*rhogeq*w1i**4 + 2*rhogeq*w2r*w1r*k*Kdrag*vdsol*w1i*w3i - rhogeq*w1i*k*Kdrag*vdsol*w3i*w1r**2 +&
        rhogeq*w2r*w1i**4*w3i*rhogsol + 2*rhogeq*w1i*w2r**2*w3r**2*rhogsol*w1r - rhogeq*w1i*w3r*rhogsol*w2r**2*w1r**2 -&
        rhogeq*w1i**3*w3r**2*rhogsol*w2r - 2*rhogeq**2*w1i*w2r**2*k*vgsol*w1r*w3r + Kdrag**2*w1i**3*k*vdsol -&
        2*rhogeq**2*w1i*w3i**2*k*vgsol*w2r*w1r - 2*cs**2*k**3*rhogeq**2*w1i*vgsol*w1r*w3r -&
        2*cs**2*k**3*rhogeq**2*w1i*vgsol*w2r*w1r - Kdrag*w2r*k**2*cs**2*rhogsol*w1i**2 -&
        2*Kdrag*rhogeq*w1r**2*k*vgsol*w1i**2 - Kdrag*w2r*k**2*cs**2*rhogsol*w1r**2 + rhogeq**2*w1i**3*w3i**2*k*vgsol -&
        Kdrag**2*k*vdsol*w3i*w1r**2 + 4*rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w1r*w3r +&
        2*rhogeq*w3i**2*w1r*rhogsol*w2r**2*w1i + rhogeq**2*w1i*k*vgsol*w3i**2*w1r**2 +&
        rhogeq**2*w1i*k*vgsol*w2r**2*w1r**2 - rhogeq*w3i*w1r**3*rhogsol*w2r**2 - rhogeq**2*w1i**4*k*vgsol*w3i -&
        rhogeq*w1i**3*k*Kdrag*vdsol*w3i - Kdrag**2*k*vdsol*w1i**2*w3i + cs**2*k**3*rhogeq*Kdrag*vdsol*w1i**2 -&
        2*rhogeq*w1i*k*w3i*Kdrag*vgsol*w2r*w1r - rhogeq*w1i*w3r*k**2*cs**2*rhogsol*w1r**2 -&
        rhogeq*w1i**3*w3r*k**2*cs**2*rhogsol - cs**2*k**3*rhogeq*Kdrag*vdsol*w1r**2 + vdsol*k*Kdrag*rhogeq*w1r**4 -&
        rhogeq*w3i*w1r*rhogsol*w2r**2*w1i**2 - rhogeq*w2r*w3i**2*rhogsol*w1r**2*w1i -&
        rhogeq*w1r*w3r*k*Kdrag*vdsol*w1i**2 - rhogeq*w1r**3*w3r*k*Kdrag*vdsol +&
        rhogeq*k**2*cs**2*w3i*w1r*rhogsol*w1i**2 - cs**2*k**3*rhogeq*Kdrag*vgsol*w1i**2 +&
        rhogeq**2*w1i*k*vgsol*w1r**2*w3r**2 + Kdrag*k**2*cs**2*w1r**3*rhogsol + Kdrag*w2r*w1r**3*w3r*rhogsol +&
        Kdrag**2*w1i*k*vdsol*w1r**2 + cs**2*k**3*rhogeq*Kdrag*vgsol*w1r**2 - rhogeq*w2r*w1r**3*k*Kdrag*vdsol -&
        rhogeq*w2r*w1i**3*k**2*cs**2*rhogsol - rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w1r**2 -&
        rhogeq*w1i**3*w2r**2*w3r*rhogsol + 2*rhogeq*w2r*w3i*rhogsol*w1r**2*w1i**2 + rhogeq**2*w1i**3*k*vgsol*w3r**2 -&
        rhogeq*w1i*w3r**2*rhogsol*w2r*w1r**2 - rhogeq*k*Kdrag*vdsol*w1i**2*w2r*w3r -&
        2*rhogeq**2*w1i*k*vgsol*w2r*w1r*w3r**2 + Kdrag**2*w1i**2*k*vgsol*w3i - rhogeq*w2r*w1r*k*Kdrag*vdsol*w1i**2 -&
        Kdrag*rhogeq*w1r**4*k*vgsol + rhogeq*w2r*w3i*rhogsol*w1r**4 + rhogeq**2*w1i**3*w2r**2*k*vgsol -&
        Kdrag**2*w1i*k*vgsol*w1r**2 + Kdrag*w1i*w2r*w3i*rhogsol*w1r**2 + Kdrag*w1i**3*w2r*w3i*rhogsol +&
        2*rhogeq**2*w1i**3*k*vgsol*w2r*w3r + 2*rhogeq**2*w1i*k*vgsol*w1r**2*w2r*w3r +&
        2*cs**2*k**3*rhogeq**2*w1i*vgsol*w2r*w3r - rhogeq**2*w2i**2*w1i**2*w3i*k*vgsol +&
        rhogeq**2*w2i**2*w3i*w1r**2*k*vgsol + rhogeq*w2i**2*k**2*cs**2*w3i*w1r*rhogsol -&
        rhogeq*w2i**2*w3r*Kdrag*w1r*k*vgsol + Kdrag*w2r*w1r*w3r*k**2*cs**2*rhogsol + Kdrag*w1r*w3r**2*w2r**2*rhogsol -&
        Kdrag*w3r*w2r**2*w1r**2*rhogsol + Kdrag*w3i**2*w1r*w2r**2*rhogsol + w2i**2*Kdrag*w3i**2*w1r*rhogsol -&
        w2i**2*Kdrag*w3r*w1r**2*rhogsol + w2i**2*Kdrag*w1r*w3r**2*rhogsol - w2i**2*Kdrag*w1i**2*w3r*rhogsol -&
        Kdrag*w1i**2*w3i**2*w2r*rhogsol - Kdrag*w1i**2*w2r**2*w3r*rhogsol - Kdrag*w1i**2*w3r**2*w2r*rhogsol +&
        Kdrag**2*w3i*w2r*w1r*k*vdsol - Kdrag*w1r**2*w3r**2*w2r*rhogsol + 2*cs**2*k**3*rhogeq**2*w1i*w2i*w3i*vgsol +&
        cs**2*k**3*rhogeq**2*w1i*w3i**2*vgsol + cs**2*k**3*rhogeq**2*w1i*w3r**2*vgsol -&
        cs**2*k**3*rhogeq*w1i*w3i*Kdrag*vdsol - cs**4*k**4*rhogeq*w1i*w3r*rhogsol +&
        cs**2*k**3*rhogeq*w1i*w3i*Kdrag*vgsol + w2i**2*Kdrag*rhogeq*w1r**2*k*vgsol + Kdrag*w2r**2*rhogeq*w1r**2*k*vgsol&
        + w2i**2*Kdrag*w1i**2*rhogeq*k*vgsol + Kdrag*w1i**2*w2r**2*rhogeq*k*vgsol -&
        cs**2*k**3*rhogeq**2*w3i*w2r**2*vgsol + cs**2*k**2*rhogeq*w3i*w1r*rhogsol*w2r**2 -&
        cs**2*k**3*rhogeq*w2i*w3i*Kdrag*vgsol + cs**4*k**4*rhogeq*w2i*w3r*rhogsol -&
        cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol - cs**2*k**2*rhogeq*w1i*w2r**2*w3r*rhogsol -&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol + cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol +&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol + cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol +&
        cs**2*k**2*rhogeq*w2i*w3i**2*w1r*rhogsol - 2*cs**2*k**2*rhogeq*w2i*w3r*w1r**2*rhogsol +&
        cs**2*k**2*rhogeq*w2i*w1r*w3r**2*rhogsol + 2*cs**2*k**2*rhogeq*w1i**2*w2i*w3r*rhogsol -&
        cs**2*k**3*rhogeq**2*w1i**2*w3i*vgsol + cs**2*k**3*rhogeq**2*w3i*w1r**2*vgsol +&
        cs**2*k**3*rhogeq*w1r*w3r*Kdrag*vdsol - cs**4*k**4*rhogeq*w3i*w1r*rhogsol -&
        cs**2*k**3*rhogeq*w3r*Kdrag*w1r*vgsol - cs**2*k**3*rhogeq**2*w2i*w3i**2*vgsol -&
        cs**2*k**3*rhogeq**2*w2i*w3r**2*vgsol - rhogeq*w2i**2*w1i*w3r*k**2*cs**2*rhogsol -&
        rhogeq*w2i**2*w1i*k*w3i*Kdrag*vgsol + cs**2*k**3*rhogeq*w2i*w3i*Kdrag*vdsol -&
        cs**2*k**3*rhogeq**2*w2i**2*w3i*vgsol + cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol -&
        cs**4*k**4*rhogeq*w2i*w1r*rhogsol - Kdrag**2*w3i*w2r*w1r*k*vgsol - Kdrag*w3i**2*w1r**2*w2r*rhogsol -&
        w2i*Kdrag*w1i*rhogeq*w3r**2*k*vgsol + w2i*Kdrag**2*w1i*k*w3i*vdsol + w2i*Kdrag*w1i*w3r*k**2*cs**2*rhogsol +&
        Kdrag*w1i*w3i*w2r*k**2*cs**2*rhogsol - Kdrag**2*w1i*w2r*w3r*k*vdsol + Kdrag**2*w1i*w2r*w3r*k*vgsol -&
        w2i*Kdrag**2*w1i*k*w3i*vgsol - w2r**2*rhogeq**2*w1i**2*w3i*k*vgsol + w2r**2*rhogeq**2*w3i*w1r**2*k*vgsol -&
        2*rhogeq*w2r*k**2*cs**2*rhogsol*w3i*w1r**2 - rhogeq*w2r*w1i*w3r**2*k**2*cs**2*rhogsol -&
        rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w3i**2 + 2*rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i +&
        rhogeq*w2r*cs**2*k**3*w3r*Kdrag*vgsol - rhogeq*w2r*cs**4*k**4*w1i*rhogsol -&
        rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol + rhogeq*w2r*cs**4*k**4*w3i*rhogsol +&
        rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol - rhogeq*w2r*cs**2*k**3*w3r*Kdrag*vdsol - w2i*rhogeq**2*w1i**4*k*vgsol +&
        w2i*Kdrag*w1i**2*w3i*rhogeq*k*vgsol + 3*Kdrag*w1i**2*rhogeq*w2r*w3r*k*vgsol +&
        3*w2i*Kdrag*rhogeq*w3i*w1r**2*k*vgsol + w2i*Kdrag**2*w1r*w3r*k*vdsol - w2i*Kdrag*k**2*cs**2*w3i*w1r*rhogsol -&
        w2i*Kdrag**2*w3r*w1r*k*vgsol - Kdrag*rhogeq*w3r*w2r**2*k*w1r*vgsol + Kdrag*rhogeq*w3r*w2r*w1r**2*k*vgsol -&
        Kdrag*rhogeq*w3i**2*k*vgsol*w1r*w2r - Kdrag*rhogeq*w3r**2*k*vgsol*w1r*w2r - Kdrag*w1i*w2r**2*w3i*rhogeq*k*vgsol&
        - w2i*Kdrag*w1i*rhogeq*w3i**2*k*vgsol + Kdrag*w1i**2*rhogeq*k*vgsol*w3i**2 + Kdrag*w1i**2*rhogeq*k*vgsol*w3r**2&
        + Kdrag*rhogeq*w1r**2*k*vgsol*w3r**2 - w2i*vdsol*k*Kdrag*rhogeq*w3i*w1r**2 + Kdrag*rhogeq*w1r**2*k*vgsol*w3i**2&
        + 2*w2i*w3i*rhogeq**2*w1i**3*k*vgsol + 2*w2i*w3i*rhogeq**2*w1i*k*vgsol*w1r**2 + w2i*Kdrag**2*k*vgsol*w1r**2 +&
        w2i*rhogeq*k**2*cs**2*w1r*rhogsol*w1i**2 + w2i*rhogeq*k**2*cs**2*w1r**3*rhogsol +&
        2*w2i*rhogeq*w3r*w1r**2*rhogsol*w1i**2 - w2i*rhogeq*w1i**3*k*Kdrag*vdsol - w2i*rhogeq*w1i*k*Kdrag*vdsol*w1r**2&
        - 2*w2i*rhogeq**2*w1i**2*k*vgsol*w1r**2 - w2i*rhogeq*w3i**2*w1r*rhogsol*w1i**2 -&
        2*w2i*rhogeq*w3r*Kdrag*w1r*k*vgsol*w1i + w2i*Kdrag**2*k*vgsol*w1i**2 - w2i*Kdrag**2*k*vdsol*w1r**2 +&
        w2i*rhogeq*w1i**2*k*Kdrag*vdsol*w3i + 2*w2i*rhogeq*w1i*k*Kdrag*vdsol*w1r*w3r - w2i*Kdrag*w3i*w1r*rhogsol*w1i**2&
        - w2i*rhogeq*w3i**2*w1r**3*rhogsol - w2i*rhogeq*w1i**2*w3r**2*rhogsol*w1r - w2i*rhogeq*w3r**2*w1r**3*rhogsol +&
        w2i*rhogeq*w3r*w1r**4*rhogsol - 2*w2i**2*rhogeq**2*w1i*k*vgsol*w1r*w3r - w2i**2*rhogeq*w1i*w3r*rhogsol*w1r**2 -&
        w2i**2*rhogeq*w3i*w1r**3*rhogsol + 2*w2i**2*rhogeq*w1i*w3r**2*rhogsol*w1r - w2i**2*rhogeq*w1i**3*w3r*rhogsol +&
        w2i**2*rhogeq**2*w1i**3*k*vgsol + w2i**2*rhogeq**2*w1i*k*vgsol*w1r**2 - w2i**2*rhogeq*w3i*w1r*rhogsol*w1i**2 +&
        w2i*Kdrag*w1i*w3r*rhogsol*w1r**2 - 4*w2i*rhogeq*k**2*cs**2*w1r*rhogsol*w1i*w3i - w2i*Kdrag*w3i*w1r**3*rhogsol +&
        2*w2i**2*rhogeq*w3i**2*w1r*rhogsol*w1i + w2i*rhogeq*w1i**4*w3r*rhogsol + w2i*Kdrag*w1i**3*w3r*rhogsol -&
        w2i*rhogeq**2*w1r**4*k*vgsol - w2i*Kdrag**2*k*vdsol*w1i**2 - Kdrag**2*w1i**3*k*vgsol)/( - 2*w1i*w3i + w1i**2 +&
        w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 -&
        2*w2r*w1r)/rhogeq/Kdrag/k 
  endif

!-------------------------------
! G A S  D E N S I T I E S
!-------------------------------
  rhog3r =(w2i*w3i**2*w1r*rhogsol + w1i*w3i**2*w2r*rhogsol - w3i**2*k*Kdrag*vgsol + w3i**2*k*Kdrag*vdsol +&
        2*w1r*w2r*w3i*rhogeq*k*vgsol + w2r**2*w3i*rhogeq*k*vgsol + w2i**2*w3i*rhogeq*k*vgsol - w1i**2*w2r*w3i*rhogsol -&
        w1r**2*w2r*w3i*rhogsol - w2r**2*w3i*w1r*rhogsol - w2i**2*w3i*w1r*rhogsol + 2*k*w3i*w2i*rhogeq*w1i*vgsol +&
        w1i**2*rhogeq*w3i*k*vgsol + 2*w3i*w3r*w2r*w1r*rhogsol + w1i*w3i*k*Kdrag*vgsol - w3i*w1r*k**2*cs**2*rhogsol +&
        2*w3i*w3r*k**2*cs**2*rhogsol + w3r**2*k*Kdrag*vgsol - w3r**2*k*Kdrag*vdsol - w3i*w2r*k**2*cs**2*rhogsol +&
        w2i*w3i*k*Kdrag*vgsol - w2i*rhogeq*w3i**2*k*vgsol - w1i*rhogeq*w3i**2*k*vgsol - w1i*w3r**2*w2r*rhogsol -&
        w2i*w3r**2*w1r*rhogsol - w2i*w3i*k*Kdrag*vdsol - w2i*rhogeq*w1r**2*k*vgsol - 2*w1i*w2i*w3i*w3r*rhogsol -&
        w1i**2*w2i*rhogeq*k*vgsol - w2r*w1r*k*Kdrag*vdsol + w2r*w1r*k*Kdrag*vgsol + w3r*w2i*w1r**2*rhogsol +&
        w3r*w1i*w2r**2*rhogsol + w3r*w1i*rhogsol*w2i**2 + w3r*w1i**2*rhogsol*w2i - w1i*w3i*k*Kdrag*vdsol +&
        rhogeq*w3i*w1r**2*k*vgsol + w2i*rhogeq*w3r**2*k*vgsol - 2*rhogeq*w3i*w3r*k*w1r*vgsol -&
        2*rhogeq*w3i*w3r*w2r*k*vgsol + w2i*k**2*cs**2*w1r*rhogsol - w1i*w2i**2*rhogeq*k*vgsol -&
        w1i*w2r**2*rhogeq*k*vgsol + w1i*w2i*k*Kdrag*vdsol + w3r*w2r*k*Kdrag*vdsol - w1i*w2i*k*Kdrag*vgsol +&
        w1i*rhogeq*w3r**2*k*vgsol + w1i*w2r*k**2*cs**2*rhogsol - w3r*w1i*k**2*cs**2*rhogsol - w3r*w2r*k*Kdrag*vgsol +&
        w3r*k*Kdrag*vdsol*w1r - w3r*k*Kdrag*vgsol*w1r - w3r*w2i*k**2*cs**2*rhogsol)/( - 2*w2r*w3r + w2r**2 + w3r**2 +&
        w3i**2 + w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  rhog3i =(w2i**2*w1r**2*rhogsol + w2r**2*w1r**2*rhogsol - w3i**2*w2r*w1r*rhogsol - w3r*w1r*w2r**2*rhogsol -&
        w3r*w1r**2*w2r*rhogsol + w3r**2*k**2*cs**2*rhogsol + w3r**2*w2r*w1r*rhogsol - w3i**2*k**2*cs**2*rhogsol -&
        w2i**2*w3r*w1r*rhogsol - w2i*w3i*w1r**2*rhogsol + w1i*w2i*rhogsol*w3i**2 - w1i*w2i*rhogsol*w3r**2 -&
        w1i**2*w3i*rhogsol*w2i - w1i*w3i*w2r**2*rhogsol - w1i**2*w3r*w2r*rhogsol - w2i**2*rhogeq*k*w1r*vgsol -&
        2*w2i*w3r*rhogeq*k*w3i*vgsol + 2*w2i*w3r*w3i*w1r*rhogsol + w2i*w3r*k*Kdrag*vgsol - w2i*Kdrag*w1r*k*vgsol +&
        w2i*w3i*k**2*cs**2*rhogsol + w2i*w1r*k*Kdrag*vdsol - w2i*w3r*k*Kdrag*vdsol + 2*w3r*rhogeq*k*vgsol*w1r*w2r -&
        w3r*w2r*k**2*cs**2*rhogsol + w3r*rhogeq*w2r**2*k*vgsol + 2*w3r*k*w3i*Kdrag*vdsol + w3r*rhogeq*w1r**2*k*vgsol -&
        2*w3r*k*w3i*Kdrag*vgsol - w3r*w1r*k**2*cs**2*rhogsol - w3i*w2r*k*Kdrag*vdsol + w2r*w1r*k**2*cs**2*rhogsol -&
        rhogeq*w2r**2*k*w1r*vgsol + rhogeq*w3i**2*w2r*k*vgsol + rhogeq*w3i**2*k*w1r*vgsol - rhogeq*w2r*w1r**2*k*vgsol +&
        w3i*w2r*k*Kdrag*vgsol - k*Kdrag*vdsol*w3i*w1r + k*Kdrag*vgsol*w3i*w1r - w3r**2*rhogeq*w2r*k*vgsol -&
        w3r**2*rhogeq*k*w1r*vgsol - w1i*w2i*k**2*cs**2*rhogsol + 2*w1i*w3r*w2r*w3i*rhogsol -&
        2*w1i*w3r*rhogeq*k*w3i*vgsol + w1i*w3r*k*Kdrag*vgsol + w1i**2*rhogsol*w2i**2 - w1i*w3i*rhogsol*w2i**2 +&
        w1i**2*w2r**2*rhogsol + 2*w1i*w2i*w3r*rhogeq*k*vgsol + w2i**2*w3r*rhogeq*k*vgsol - w1i*w3r*k*Kdrag*vdsol +&
        w1i*w3i*k**2*cs**2*rhogsol + w1i*w2r*k*Kdrag*vdsol - w1i*w2r*k*Kdrag*vgsol + w1i**2*w3r*rhogeq*k*vgsol -&
        w1i**2*rhogeq*w2r*k*vgsol)/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i +&
        w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  rhog2r = - (w2i*w3i**2*w1r*rhogsol - w1i*w3i**2*w2r*rhogsol - w2r**2*w3i*rhogeq*k*vgsol + w2i**2*w3i*rhogeq*k*vgsol -&
        w1i**2*w2r*w3i*rhogsol - w1r**2*w2r*w3i*rhogsol + w2r**2*w3i*w1r*rhogsol - w2i**2*w3i*w1r*rhogsol -&
        2*k*w3i*w2i*rhogeq*w1i*vgsol - w2i**2*k*Kdrag*vdsol + w1i**2*rhogeq*w3i*k*vgsol + w1i*w3i*k*Kdrag*vgsol -&
        w3i*w1r*k**2*cs**2*rhogsol + w3i*w2r*k**2*cs**2*rhogsol - w2i*w3i*k*Kdrag*vgsol - w2i*rhogeq*w3i**2*k*vgsol +&
        w1i*rhogeq*w3i**2*k*vgsol - w1i*w3r**2*w2r*rhogsol + w2i*w3r**2*w1r*rhogsol - w2r**2*k*Kdrag*vgsol +&
        w2i*w3i*k*Kdrag*vdsol + w2i**2*k*Kdrag*vgsol - w2i*rhogeq*w1r**2*k*vgsol + 2*w2i*w1r*rhogeq*w2r*k*vgsol -&
        2*w2r*w1r*w2i*w3r*rhogsol - w1i**2*w2i*rhogeq*k*vgsol - w2r*w1r*k*Kdrag*vdsol + w2r*w1r*k*Kdrag*vgsol +&
        w2r**2*k*Kdrag*vdsol + 2*w2i*rhogeq*w2r*w3r*k*vgsol + w3r*w2i*w1r**2*rhogsol - 2*w2i*w2r*k**2*cs**2*rhogsol +&
        w3r*w1i*w2r**2*rhogsol - w3r*w1i*rhogsol*w2i**2 + w3r*w1i**2*rhogsol*w2i - w1i*w3i*k*Kdrag*vdsol +&
        rhogeq*w3i*w1r**2*k*vgsol - w2i*rhogeq*w3r**2*k*vgsol + w2i*k**2*cs**2*w1r*rhogsol + 2*w1i*w2r*w3i*rhogsol*w2i&
        + w1i*w2i**2*rhogeq*k*vgsol - w1i*w2r**2*rhogeq*k*vgsol + w1i*w2i*k*Kdrag*vdsol - w3r*w2r*k*Kdrag*vdsol -&
        w1i*w2i*k*Kdrag*vgsol + w1i*rhogeq*w3r**2*k*vgsol - 2*w2i*rhogeq*w3r*k*w1r*vgsol + w1i*w2r*k**2*cs**2*rhogsol -&
        w3r*w1i*k**2*cs**2*rhogsol + w3r*w2r*k*Kdrag*vgsol + w3r*k*Kdrag*vdsol*w1r - w3r*k*Kdrag*vgsol*w1r +&
        w3r*w2i*k**2*cs**2*rhogsol)/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 - 2*w3i*w2i)/(w2i**2 + w1r**2 -&
        2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  rhog2i = - ( - w3r**2*rhogsol*w1r**2 - w1i**2*w3i**2*rhogsol - w3i**2*rhogsol*w1r**2 + w2i**2*k**2*cs**2*rhogsol +&
        w3i**2*w2r*w1r*rhogsol - w3r*w1r*w2r**2*rhogsol + w3r*w1r**2*w2r*rhogsol + w3r**2*w2r*w1r*rhogsol +&
        w2i**2*w3r*w1r*rhogsol + w2i*w3i*w1r**2*rhogsol - k**2*cs**2*rhogsol*w2r**2 - 2*w2i*w2r*k*Kdrag*vdsol +&
        w1i*w2i*rhogsol*w3i**2 + 2*w2i*w2r*k*Kdrag*vgsol - 2*w2i*w2r*w3i*w1r*rhogsol + w1i*w2i*rhogsol*w3r**2 +&
        w1i**2*w3i*rhogsol*w2i + w1i*w3i*w2r**2*rhogsol - w1i**2*w3r**2*rhogsol - 2*w1i*w2i*w2r*w3r*rhogsol -&
        2*w1i*w2r*rhogeq*w3i*k*vgsol + w1i**2*w3r*w2r*rhogsol - w2i**2*rhogeq*k*w1r*vgsol - w2i*w3r*k*Kdrag*vgsol -&
        w2i*Kdrag*w1r*k*vgsol - w2i*w3i*k**2*cs**2*rhogsol + w2i*w1r*k*Kdrag*vdsol + w2i*w3r*k*Kdrag*vdsol -&
        2*w3r*rhogeq*k*vgsol*w1r*w2r + w3r*w2r*k**2*cs**2*rhogsol + w3r*rhogeq*w2r**2*k*vgsol +&
        w3r*rhogeq*w1r**2*k*vgsol - w3r*w1r*k**2*cs**2*rhogsol + w3i*w2r*k*Kdrag*vdsol + w2r*w1r*k**2*cs**2*rhogsol +&
        rhogeq*w2r**2*k*w1r*vgsol - rhogeq*w3i**2*w2r*k*vgsol + rhogeq*w3i**2*k*w1r*vgsol - rhogeq*w2r*w1r**2*k*vgsol -&
        w3i*w2r*k*Kdrag*vgsol - k*Kdrag*vdsol*w3i*w1r + k*Kdrag*vgsol*w3i*w1r - w3r**2*rhogeq*w2r*k*vgsol +&
        w3r**2*rhogeq*k*w1r*vgsol + 2*w1i*w2i*w2r*rhogeq*k*vgsol + 2*w2i*w2r*rhogeq*w3i*k*vgsol -&
        w1i*w2i*k**2*cs**2*rhogsol + w1i*w3r*k*Kdrag*vgsol - w1i*w3i*rhogsol*w2i**2 - w2i**2*w3r*rhogeq*k*vgsol -&
        w1i*w3r*k*Kdrag*vdsol + w1i*w3i*k**2*cs**2*rhogsol + w1i*w2r*k*Kdrag*vdsol - w1i*w2r*k*Kdrag*vgsol +&
        w1i**2*w3r*rhogeq*k*vgsol - w1i**2*rhogeq*w2r*k*vgsol)/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 -&
        2*w3i*w2i)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  rhog1r =(w2i*w3i**2*w1r*rhogsol - w1i*w3i**2*w2r*rhogsol + k*Kdrag*vdsol*w1i**2 - k*Kdrag*vdsol*w1r**2 -&
        w2r**2*w3i*rhogeq*k*vgsol - w2i**2*w3i*rhogeq*k*vgsol + w1i**2*w2r*w3i*rhogsol - w1r**2*w2r*w3i*rhogsol +&
        w2r**2*w3i*w1r*rhogsol + w2i**2*w3i*w1r*rhogsol - k*Kdrag*vgsol*w1i**2 + k*Kdrag*vgsol*w1r**2 +&
        2*k*w3i*w2i*rhogeq*w1i*vgsol - w1i**2*rhogeq*w3i*k*vgsol + w1i*w3i*k*Kdrag*vgsol - w3i*w1r*k**2*cs**2*rhogsol +&
        w3i*w2r*k**2*cs**2*rhogsol - w2i*w3i*k*Kdrag*vgsol - w2i*rhogeq*w3i**2*k*vgsol + w1i*rhogeq*w3i**2*k*vgsol -&
        w1i*w3r**2*w2r*rhogsol + w2i*w3r**2*w1r*rhogsol + w2i*w3i*k*Kdrag*vdsol + w2i*rhogeq*w1r**2*k*vgsol +&
        2*k**2*cs**2*w1r*rhogsol*w1i - 2*w1i*rhogeq*k*vgsol*w2r*w1r - w1i**2*w2i*rhogeq*k*vgsol + w2r*w1r*k*Kdrag*vdsol&
        - w2r*w1r*k*Kdrag*vgsol + 2*w1i*rhogeq*k*vgsol*w2r*w3r - w3r*w2i*w1r**2*rhogsol - w3r*w1i*w2r**2*rhogsol -&
        2*w1i*rhogeq*k*vgsol*w1r*w3r - w3r*w1i*rhogsol*w2i**2 + w3r*w1i**2*rhogsol*w2i - w1i*w3i*k*Kdrag*vdsol +&
        rhogeq*w3i*w1r**2*k*vgsol + 2*w3r*w1i*rhogsol*w2r*w1r - w2i*rhogeq*w3r**2*k*vgsol - w2i*k**2*cs**2*w1r*rhogsol&
        + w1i*w2i**2*rhogeq*k*vgsol + w1i*w2r**2*rhogeq*k*vgsol - w1i*w2i*k*Kdrag*vdsol - w3r*w2r*k*Kdrag*vdsol +&
        w1i*w2i*k*Kdrag*vgsol + w1i*rhogeq*w3r**2*k*vgsol - w1i*w2r*k**2*cs**2*rhogsol - w3r*w1i*k**2*cs**2*rhogsol +&
        w3r*w2r*k*Kdrag*vgsol + w3r*k*Kdrag*vdsol*w1r - w3r*k*Kdrag*vgsol*w1r + w3r*w2i*k**2*cs**2*rhogsol -&
        2*w2i*w3i*w1r*rhogsol*w1i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w2i**2 + w1r**2 -&
        2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  rhog1i =(w3i**2*rhogsol*w2i**2 + 2*w2i*w3i*rhogeq*k*w1r*vgsol - k**2*cs**2*rhogsol*w1i**2 -&
        2*rhogeq*k*w1r*vgsol*w1i*w3i + 2*w1i*w3i*rhogsol*w2r*w1r - 2*Kdrag*w1r*k*vgsol*w1i + 2*w1r*k*Kdrag*vdsol*w1i -&
        w3i**2*w2r*w1r*rhogsol - w3r*w1r*w2r**2*rhogsol + w3r*w1r**2*w2r*rhogsol - w3r**2*w2r*w1r*rhogsol -&
        w2i**2*w3r*w1r*rhogsol - w2i*w3i*w1r**2*rhogsol + 2*w2i*w3r*w1r*rhogsol*w1i - 2*w2i*rhogeq*k*w1r*vgsol*w1i -&
        w1i*w2i*rhogsol*w3i**2 - w1i*w2i*rhogsol*w3r**2 + w1i**2*w3i*rhogsol*w2i - w1i*w3i*w2r**2*rhogsol -&
        w1i**2*w3r*w2r*rhogsol + w2i**2*rhogeq*k*w1r*vgsol - w2i*w3r*k*Kdrag*vgsol + w2i*Kdrag*w1r*k*vgsol -&
        w2i*w3i*k**2*cs**2*rhogsol - w2i*w1r*k*Kdrag*vdsol + w2i*w3r*k*Kdrag*vdsol + 2*w3r*rhogeq*k*vgsol*w1r*w2r +&
        w3r*w2r*k**2*cs**2*rhogsol - w3r*rhogeq*w2r**2*k*vgsol - w3r*rhogeq*w1r**2*k*vgsol - w3r*w1r*k**2*cs**2*rhogsol&
        + w3i*w2r*k*Kdrag*vdsol - w2r*w1r*k**2*cs**2*rhogsol + rhogeq*w2r**2*k*w1r*vgsol - rhogeq*w3i**2*w2r*k*vgsol +&
        rhogeq*w3i**2*k*w1r*vgsol - rhogeq*w2r*w1r**2*k*vgsol - w3i*w2r*k*Kdrag*vgsol - k*Kdrag*vdsol*w3i*w1r +&
        k*Kdrag*vgsol*w3i*w1r - w3r**2*rhogeq*w2r*k*vgsol + w3r**2*rhogeq*k*w1r*vgsol + w1i*w2i*k**2*cs**2*rhogsol +&
        w1i*w3r*k*Kdrag*vgsol - w1i*w3i*rhogsol*w2i**2 - w2i**2*w3r*rhogeq*k*vgsol - w1i*w3r*k*Kdrag*vdsol +&
        w1i*w3i*k**2*cs**2*rhogsol - w1i*w2r*k*Kdrag*vdsol + w1i*w2r*k*Kdrag*vgsol + w1i**2*w3r*rhogeq*k*vgsol +&
        w1i**2*rhogeq*w2r*k*vgsol + rhogsol*w3i**2*w2r**2 + rhogsol*w2r**2*w3r**2 + w2i**2*w3r**2*rhogsol +&
        k**2*cs**2*rhogsol*w1r**2)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w2i**2 + w1r**2 -&
        2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

!-------------------------------
! D U S T  D E N S I T I E S
!-------------------------------
  if (Kdrag.gt.0.) then
  rhod3r =rhodeq*(rhogeq*w3i**5*w2i*w1r*rhogsol + Kdrag**2*w3i**4*k*vdsol - rhogeq*w2r*k**2*cs**2*rhogsol*w3i**4 -&
        2*rhogeq*w2r*k**2*cs**2*rhogsol*w3i**2*w3r**2 + w3r*rhogeq*w3i**2*w1i*w2r*k*Kdrag*vdsol -&
        w3r*rhogeq**2*w3i**2*w2r*w1r**2*k*vgsol - w3r*rhogeq*w3i**2*w2i*Kdrag*w1r*k*vgsol +&
        w3r*rhogeq*w3i**2*w2i*w1r*k*Kdrag*vdsol - w3r*rhogeq*w3i**2*w1i*w2r*k*Kdrag*vgsol +&
        w3r*rhogeq*w3i**2*w1i**2*w2r**2*rhogsol - w3r*rhogeq**2*w3i**4*k*w1r*vgsol -&
        w3r*rhogeq**2*w3i**2*w1i**2*w2r*k*vgsol - 3*w3r*w3i**2*rhogeq**2*k**3*cs**2*w1r*vgsol +&
        w3r*rhogeq*w3i**2*w1i**2*rhogsol*w2i**2 + 2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol -&
        2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol + 2*w3r*w3i*rhogeq*k**3*cs**2*Kdrag*vdsol*w1r -&
        2*w3r*w3i*rhogeq*k**3*cs**2*Kdrag*vgsol*w1r + w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol -&
        3*w3r*w3i**2*rhogeq**2*k**3*cs**2*w2r*vgsol + w3r*rhogeq*w3i**4*k**2*cs**2*rhogsol +&
        2*rhogeq*w3r**2*w1i*w3i**3*w2r*rhogsol - 2*rhogeq**2*w3r**2*w1i*w3i**3*k*vgsol +&
        rhogeq*w3r**2*w1i**2*w2r*k**2*cs**2*rhogsol - rhogeq**2*w3r**3*w2r*w1r**2*k*vgsol +&
        rhogeq*w3r**3*w2r*w1i*k*Kdrag*vdsol - rhogeq**2*w3r**3*w2r*w1i**2*k*vgsol -&
        rhogeq**2*w3r**2*w3i*w2i*w1i**2*k*vgsol + 2*rhogeq*w3r**3*w3i**2*w2r*w1r*rhogsol - rhogeq**2*w3r**5*k*w1r*vgsol&
        - rhogeq**2*w3r**5*w2r*k*vgsol + 2*rhogeq*w3r**2*w1i*w3i**2*k*Kdrag*vgsol +&
        2*rhogeq*w3r**2*w3i**2*w2i*k*Kdrag*vgsol - rhogeq*w3r**5*w1i*w2i*rhogsol +&
        rhogeq*w3r**2*w3i*w2i*w1i*k*Kdrag*vgsol + rhogeq*w3r**2*w3i*w2i*w1i*k*Kdrag*vdsol +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**2 - 2*rhogeq*w3i**2*w2i**2*w1r*rhogsol*w3r**2 -&
        Kdrag**2*w3i**4*k*vgsol + 3*rhogeq**2*w3r**2*cs**2*k**3*w1i*w3i*vgsol + rhogeq*w3r**3*w2r**2*w1r**2*rhogsol -&
        rhogeq*w3r**2*cs**2*k**3*w1i*Kdrag*vgsol + Kdrag*w3i*w1i**2*w3r**2*rhogeq*k*vgsol +&
        rhogeq*w3r**2*w2r**2*w1r*k**2*cs**2*rhogsol - Kdrag*w3r**3*w1i*k**2*cs**2*rhogsol +&
        Kdrag*w3r**3*w1i*w2r**2*rhogsol - Kdrag**2*w3r**3*w2r*k*vgsol + Kdrag**2*w3r**3*k*vdsol*w1r -&
        2*w3i*rhogeq*k**2*cs**2*w1i*w3r**2*w2r*rhogsol - rhogeq*w3r**3*w2r*w1i*k*Kdrag*vgsol +&
        2*rhogeq*w3r**2*w3i**3*w2i*w1r*rhogsol - rhogeq*w3r**3*cs**4*k**4*rhogsol +&
        rhogeq*w3r**2*cs**4*k**4*w2r*rhogsol + rhogeq*w3r**2*cs**4*k**4*w1r*rhogsol -&
        rhogeq**2*w3r**2*w3i*w1i*w2i**2*k*vgsol - 3*rhogeq*w3r**2*cs**2*k**3*w3i*Kdrag*vdsol +&
        3*rhogeq*w3r**2*cs**2*k**3*w3i*Kdrag*vgsol - rhogeq*w3r**2*cs**2*k**3*w2i*Kdrag*vgsol -&
        rhogeq**2*w3r**3*w1r*w2i**2*k*vgsol + rhogeq*w3r**2*w2i**2*k**2*cs**2*w1r*rhogsol +&
        rhogeq*w3r**3*w2r**2*w1i**2*rhogsol - 2*w3i*rhogeq*k**2*cs**2*w2i*w3r**2*w1r*rhogsol -&
        Kdrag*w3r**4*w1i*w2r*rhogsol - Kdrag**2*w3r**3*k*vgsol*w1r + Kdrag*w3r**3*w1i*rhogsol*w2i**2 -&
        2*Kdrag*w3i*w3r**3*rhogeq*w2r*k*vgsol + Kdrag**2*w3i*w1i*w3r**2*k*vgsol - 2*Kdrag*w3i*w3r**3*rhogeq*k*w1r*vgsol&
        + Kdrag*w3i*w2i**2*w3r**2*rhogeq*k*vgsol - rhogeq**2*w3r**4*w3i*w1i*k*vgsol -&
        Kdrag*w3i*w3r**2*w1r**2*w2r*rhogsol - Kdrag*w3i*w3r**2*w1r*w2r**2*rhogsol - Kdrag*w3i*w1i**2*w3r**2*w2r*rhogsol&
        - Kdrag*w3i*w2i**2*w3r**2*w1r*rhogsol + Kdrag**2*w3i*w2i*w3r**2*k*vgsol - Kdrag**2*w3i*w2i*w3r**2*k*vdsol -&
        Kdrag*w3i*w3r**2*w1r*k**2*cs**2*rhogsol + rhogeq*w3r**4*w3i*w2i*w1r*rhogsol + rhogeq*w3r**4*w1i*w3i*w2r*rhogsol&
        - 2*rhogeq**2*w3r**3*w3i**2*w2r*k*vgsol - Kdrag**2*w3r**4*k*vdsol + Kdrag*w3i*w3r**2*rhogeq*w2r**2*k*vgsol +&
        Kdrag*w3i*w3r**2*rhogeq*w1r**2*k*vgsol - rhogeq*w3r**3*w1r*w2i*k*Kdrag*vgsol -&
        2*rhogeq**2*w3r**3*w3i**2*k*w1r*vgsol - 2*rhogeq*w3r**3*w3i**2*w1i*w2i*rhogsol +&
        rhogeq*w3r**3*w1r*w2i*k*Kdrag*vdsol - rhogeq**2*w3r**4*w3i*w2i*k*vgsol - rhogeq*w3r**4*w3i*k*Kdrag*vgsol +&
        2*rhogeq*w3r**3*w3i**2*k**2*cs**2*rhogsol + rhogeq*w3r**2*w2r*w1r**2*k**2*cs**2*rhogsol +&
        rhogeq**2*w3r**3*cs**2*k**3*w2r*vgsol - Kdrag**2*w3i*w1i*w3r**2*k*vdsol +&
        3*rhogeq**2*w3r**2*cs**2*k**3*w2i*w3i*vgsol - 2*rhogeq**2*w3r**2*cs**2*k**3*vgsol*w1r*w2r -&
        rhogeq**2*w3r**2*cs**2*k**3*w2r**2*vgsol + rhogeq*w3r**2*cs**2*k**3*w1i*Kdrag*vdsol -&
        2*rhogeq**2*w3r**2*cs**2*k**3*w1i*w2i*vgsol - rhogeq**2*w3r**2*cs**2*k**3*w2i**2*vgsol -&
        rhogeq**2*w3r**2*cs**2*k**3*w1r**2*vgsol + rhogeq**2*w3r**3*cs**2*k**3*w1r*vgsol +&
        2*w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i - w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2 +&
        w3r*w1i*rhogsol*cs**4*k**4*rhogeq*w2i + 2*w3r*w1i*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2 +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r - rhogeq**2*w3i**3*w2i*w1r**2*k*vgsol -&
        w3r*rhogeq*w3i**4*w1i*w2i*rhogsol - w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2 +&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r - w2i**2*Kdrag*w3i**3*w1r*rhogsol - rhogeq*w2i*w3r**4*k*Kdrag*vdsol +&
        2*cs**2*k**3*rhogeq**2*w3i**2*vgsol*w2r*w1r + 2*Kdrag*w3i*w1r*rhogsol*w2r*w3r**3 +&
        cs**2*k**3*rhogeq*Kdrag*vdsol*w3i**3 - 2*w3i*Kdrag*w1i*w2i*w3r**3*rhogsol + w3i**3*Kdrag**2*w2i*k*vgsol -&
        w3i**3*Kdrag**2*w2i*k*vdsol - 4*rhogeq*w3i**2*w1i*w2i*w3r*k**2*cs**2*rhogsol +&
        rhogeq*w3i**3*w1i*w2i*k*Kdrag*vdsol - rhogeq*w2i**2*w1r*w3r**4*rhogsol +&
        4*rhogeq**2*w3i**2*w1i*w2i*w3r**2*k*vgsol + 2*rhogeq**2*w1i*w2i*w3r**4*k*vgsol -&
        Kdrag*w1i**2*w3i**3*w2r*rhogsol - cs**4*k**4*rhogeq*w3i**2*w1r*rhogsol + Kdrag**2*w3i**2*w2r*w1r*k*vgsol +&
        Kdrag**2*w2r*w1r*k*vgsol*w3r**2 + w2i*Kdrag**2*w1i*k*w3i**2*vdsol + w2i*Kdrag**2*w1i*k*vdsol*w3r**2 +&
        Kdrag*w1i*w3i**2*w2r*k**2*cs**2*rhogsol - rhogeq*w3i**4*w1r*rhogsol*w2r**2 -&
        2*rhogeq*w3i**2*w1r*rhogsol*w2r**2*w3r**2 - rhogeq*w2i**2*k**2*cs**2*w3i**2*w1r*rhogsol -&
        Kdrag*w3i**3*w1r*w2r**2*rhogsol - rhogeq*w3i**4*w1r*k**2*cs**2*rhogsol -&
        2*rhogeq*w3i**2*w1r*k**2*cs**2*rhogsol*w3r**2 + 2*rhogeq**2*w3i**4*k*vgsol*w1r*w2r +&
        4*rhogeq**2*w3i**2*k*vgsol*w1r*w2r*w3r**2 - rhogeq**2*w3i**3*w1i**2*w2i*k*vgsol -&
        rhogeq*w1r*w3r**4*k**2*cs**2*rhogsol + 2*rhogeq**2*w3i**2*w2i**2*w3r**2*k*vgsol +&
        rhogeq**2*w2i**2*w3r**4*k*vgsol + rhogeq**2*w3i**4*w2i**2*k*vgsol - rhogeq*w3i**4*w2i*k*Kdrag*vdsol -&
        rhogeq**2*w3i**5*w1i*k*vgsol + w3i**4*Kdrag*w2i*w1r*rhogsol - w3r*Kdrag*w3i**2*w1i*k**2*cs**2*rhogsol -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r + 2*w3r*rhogsol*cs**2*k**2*rhogeq*w2i*w3i*w1r**2 -&
        rhogeq**2*w3i**5*w2i*k*vgsol + Kdrag*w3r**3*w2i*w1r**2*rhogsol - Kdrag*w3r**4*w2i*w1r*rhogsol +&
        Kdrag*w3r**3*w1i**2*rhogsol*w2i + Kdrag**2*w3r**3*w2r*k*vdsol - Kdrag*w3r**3*w2i*k**2*cs**2*rhogsol -&
        w3r*rhogeq**2*w3i**2*w2i**2*k*w1r*vgsol + w3r*rhogeq*w3i**2*w2i**2*w1r**2*rhogsol +&
        2*w3r*Kdrag*w3i**3*w2r*w1r*rhogsol + w3r*rhogeq*w3i**2*w2r**2*w1r**2*rhogsol -&
        2*w3r*Kdrag*w3i**3*w1i*w2i*rhogsol - 2*w3r*Kdrag*w3i**3*rhogeq*w2r*k*vgsol -&
        2*w3r*Kdrag*w3i**3*rhogeq*k*w1r*vgsol + 2*w3r*Kdrag*w3i**3*k**2*cs**2*rhogsol -&
        w3r*Kdrag*w3i**2*w2i*k**2*cs**2*rhogsol + w3r*rhogeq*w3i**4*w2r*w1r*rhogsol + w3r*Kdrag**2*w3i**2*w2r*k*vdsol +&
        w3r*Kdrag*w3i**2*w2i*w1r**2*rhogsol + w3r*Kdrag*w3i**2*w1i**2*rhogsol*w2i + w3r*Kdrag*w3i**2*w1i*w2r**2*rhogsol&
        + w3r*Kdrag**2*w3i**2*k*vdsol*w1r + w3r*rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol -&
        2*w3r*w3i*rhogeq*k**4*cs**4*w2i*rhogsol + 3*w3r*w3i**2*rhogeq*k**4*cs**4*rhogsol -&
        w3r*rhogeq**2*w3i**2*w2r**2*k*w1r*vgsol - w3r*rhogeq**2*w3i**4*w2r*k*vgsol -&
        w3r*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**2 + w3r*Kdrag*w3i**2*w1i*rhogsol*w2i**2 -&
        w3r*rhogsol*cs**4*k**4*rhogeq*w2r*w1r + w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r -&
        w3r*Kdrag**2*w3i**2*w2r*k*vgsol - cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol*w3i -&
        cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol*w3i + 2*cs**2*k**2*rhogeq*w2i*w3i**3*w1r*rhogsol +&
        rhogeq**2*w3i**4*k*vgsol*w2r**2 - cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol*w3i +&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol*w3i - cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol*w3i +&
        2*rhogeq**2*w3r**4*k*w1r*vgsol*w2r - rhogeq*w2r*w1i**2*rhogsol*w3r**4 + rhogeq**2*k*vgsol*w3r**4*w2r**2 -&
        rhogeq*w2r*rhogsol*w1r**2*w3r**4 - rhogeq*w3i**4*w2i**2*w1r*rhogsol - rhogeq*w1i*k*Kdrag*vdsol*w3r**4 -&
        rhogeq*w2r*w3i**4*rhogsol*w1r**2 - 2*rhogeq*w2r*w3i**2*rhogsol*w1r**2*w3r**2 - w2i*Kdrag**2*w1i*k*w3i**2*vgsol&
        - w2r**2*rhogeq**2*w1i*w3i**3*k*vgsol - rhogeq*w2r*k**2*cs**2*rhogsol*w3i**2*w1r**2 -&
        rhogeq*w2i**2*w1i*k*w3i**2*Kdrag*vgsol - rhogeq*w2i**2*w1i*k*Kdrag*vgsol*w3r**2 -&
        cs**2*k**3*rhogeq*w2i*w3i**2*Kdrag*vdsol + cs**2*k**3*rhogeq**2*w2i**2*w3i**2*vgsol -&
        cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol*w3i + cs**4*k**4*rhogeq*w2i*w1r*rhogsol*w3i -&
        Kdrag*w3i**3*w1r**2*w2r*rhogsol - cs**2*k**3*rhogeq**2*w2i*w3i**3*vgsol - rhogeq**2*w2i**2*w1i*w3i**3*k*vgsol +&
        cs**2*k**3*rhogeq**2*w1i**2*w3i**2*vgsol - rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol*w3i -&
        w2i*Kdrag*w1i**2*rhogeq*k*vgsol*w3r**2 - Kdrag*w1i*w2r**2*w3i**2*rhogeq*k*vgsol -&
        Kdrag*w1i*w2r**2*rhogeq*k*vgsol*w3r**2 + w2i*Kdrag*w1i*rhogeq*w3i**3*k*vgsol +&
        2*rhogeq**2*w3i**2*k*vgsol*w2r**2*w3r**2 + Kdrag**2*w1i*k*vgsol*w3i**3 - Kdrag**2*w1i*k*vdsol*w3i**3 +&
        2*rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w3i**3 - rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i**2 +&
        rhogeq*w2r*cs**4*k**4*w1i*rhogsol*w3i + rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol*w3i -&
        rhogeq*w2r*cs**4*k**4*w3i**2*rhogsol - w2i*Kdrag**2*w1i*k*vgsol*w3r**2 +&
        2*rhogeq**2*w3i**2*w1r**2*k*vgsol*w3r**2 - cs**2*k**3*rhogeq*Kdrag*vgsol*w3i**3 +&
        Kdrag*w1i**2*rhogeq*k*vgsol*w3i**3 + 2*Kdrag*w3r**3*k**2*cs**2*rhogsol*w3i + Kdrag*w1i*w2r*w3i**4*rhogsol +&
        rhogeq**2*w1i**2*w3i**4*k*vgsol + 2*rhogeq**2*w1i**2*w3i**2*k*vgsol*w3r**2 -&
        Kdrag*k**2*cs**2*w1r*rhogsol*w3i**3 - Kdrag*w2r*k**2*cs**2*rhogsol*w3i**3 + rhogeq*w2r*w1i*rhogsol*w3i**5 +&
        Kdrag*rhogeq*w2r**2*w3i**3*k*vgsol + 3*Kdrag*rhogeq*w3i**3*k*vgsol*w1r*w2r -&
        w2i*Kdrag*rhogeq*w3i**2*w1r**2*k*vgsol - w2i*Kdrag*rhogeq*w1r**2*k*vgsol*w3r**2 +&
        w2i*Kdrag*k**2*cs**2*w3i**2*w1r*rhogsol + w2i*Kdrag*k**2*cs**2*w1r*rhogsol*w3r**2 -&
        w2i*Kdrag*w1i**2*w3i**2*rhogeq*k*vgsol - rhogeq*w3r**4*k**2*cs**2*rhogsol*w2r +&
        w2i**2*Kdrag*rhogeq*w3i**3*k*vgsol + Kdrag*w1i*w2r*k**2*cs**2*rhogsol*w3r**2 +&
        cs**2*k**3*rhogeq**2*w3i**2*w1r**2*vgsol - 2*rhogeq*w3r**2*w3i**3*k*Kdrag*vgsol + rhogeq*w3r**5*w2r*w1r*rhogsol&
        + rhogeq*w3r**3*w2i**2*w1i**2*rhogsol - rhogeq*w3r**2*w3i*w2r*w1r*k*Kdrag*vdsol +&
        rhogeq*w3r**3*w1r**2*w2i**2*rhogsol + 3*rhogeq*w3r**2*w3i*w2r*w1r*k*Kdrag*vgsol +&
        rhogeq*w3r**4*w3i*k*Kdrag*vdsol - rhogeq**2*w3r**3*w2r**2*k*vgsol*w1r - rhogeq**2*w3r**2*w3i*w2i*w1r**2*k*vgsol&
        + 2*rhogeq*w3r**4*w1i*k*Kdrag*vgsol - 2*rhogeq**2*w3r**2*w3i**3*w2i*k*vgsol +&
        2*rhogeq*w3r**2*w3i**3*k*Kdrag*vdsol + rhogeq*w3r**5*k**2*cs**2*rhogsol -&
        rhogeq**2*w3r**2*w1i*w3i*w2r**2*k*vgsol + 2*rhogeq*w3r**4*w2i*k*Kdrag*vgsol -&
        Kdrag*w3i*w3r**2*w2r*k**2*cs**2*rhogsol + rhogeq*w3r**2*cs**2*k**3*w2i*Kdrag*vdsol -&
        rhogeq**2*w3r**2*cs**2*k**3*w1i**2*vgsol + 4*rhogeq*w3i**2*w2r*w1r*w3r*k**2*cs**2*rhogsol -&
        2*w3r*w3i*rhogeq*k**4*cs**4*w1i*rhogsol - w3r*Kdrag**2*w3i**2*k*vgsol*w1r + Kdrag**2*w3r**4*k*vgsol -&
        2*rhogeq*w3i**2*w2i*w3r**2*k*Kdrag*vdsol + 2*rhogeq**2*w3i**4*w1i*w2i*k*vgsol -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r + w3r*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r -&
        w3r*rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**2 - cs**2*k**3*rhogeq*w1i*w3i**2*Kdrag*vdsol +&
        cs**2*k**3*rhogeq*w1i*w3i**2*Kdrag*vgsol + cs**2*k**3*rhogeq**2*w3i**2*w2r**2*vgsol -&
        cs**2*k**2*rhogeq*w3i**2*w1r*rhogsol*w2r**2 + cs**2*k**3*rhogeq*w2i*w3i**2*Kdrag*vgsol -&
        Kdrag**2*w3i**2*w2r*w1r*k*vdsol - Kdrag**2*w2r*w1r*k*vdsol*w3r**2 + 2*cs**2*k**3*rhogeq**2*w1i*w2i*w3i**2*vgsol&
        - cs**2*k**3*rhogeq**2*w1i*w3i**3*vgsol - rhogeq*w2r*w1r*k*Kdrag*vdsol*w3i**3 -&
        rhogeq*w2r*w1i**2*w3i**4*rhogsol - 2*rhogeq*w2r*w1i**2*w3i**2*rhogsol*w3r**2 +&
        Kdrag*rhogeq*w1r**2*k*vgsol*w3i**3 - 2*rhogeq*w1i*k*w3i**2*Kdrag*vdsol*w3r**2 + rhogeq**2*w1i**2*k*vgsol*w3r**4&
        - Kdrag*rhogeq*w3i**5*k*vgsol + rhogeq**2*w1r**2*k*vgsol*w3r**4 + rhogeq*k*Kdrag*vdsol*w3i**5 -&
        rhogeq*w1i*k*w3i**4*Kdrag*vdsol - rhogeq*w1r*rhogsol*w2r**2*w3r**4 +&
        rhogeq**2*w3i**4*w1r**2*k*vgsol)/rhogeq/Kdrag/(w3i**2 + w3r**2)/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 +&
        w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  rhod3i =rhodeq*(vdsol*Kdrag*cs**2*k**3*rhogeq*w3i*w2i*w1r - rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2*w1r**2 +&
        rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i*w1r**2 + rhogsol*cs**2*k**2*rhogeq*w2i*w1r**2*w3r**2 -&
        vgsol*Kdrag*k*rhogeq*w1r*w3r**2*w2i**2 - vgsol*Kdrag*k*rhogeq*w3i**2*w1r*w2i**2 -&
        Kdrag*w3i**3*w1i*rhogsol*w2i**2 + rhogsol*cs**4*k**4*rhogeq*w2i*w3r*w1r + rhogsol*cs**4*k**4*rhogeq*w3i*w2r*w1r&
        - rhogeq**2*w3r*w3i**4*w2i*k*vgsol + rhogeq*w3r*w3i**4*k*Kdrag*vdsol - rhogeq*w3r**4*k**2*cs**2*rhogsol*w3i -&
        rhogeq*w3r*w3i**4*k*Kdrag*vgsol + rhogeq**2*w3r**4*w3i*k*w1r*vgsol + rhogeq**2*w3r**4*w3i*w2r*k*vgsol -&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w3i*w2i*w1r + rhogsol*Kdrag*w3i**2*w2i**2*w1r**2 +&
        rhogsol*Kdrag*w2i**2*w1r**2*w3r**2 + rhogsol*Kdrag*w2r**2*w1r**2*w3i**2 - 2*rhogeq*w3r*w1i*w3i**3*k*Kdrag*vgsol&
        - 2*rhogeq*w3r*w3i**3*w2i*k*Kdrag*vgsol + rhogeq*w3r**4*w3i*w1i*w2i*rhogsol +&
        3*rhogeq*w3r*w3i**2*w2i*w1i*k*Kdrag*vgsol - rhogeq*w3r*w3i**2*w2i*w1i*k*Kdrag*vdsol +&
        rhogeq**2*w3r**3*cs**2*k**3*w1i*vgsol - rhogeq*w3r**4*w3i*w2r*w1r*rhogsol + rhogeq*w3r*w1i*w3i**4*w2r*rhogsol -&
        rhogeq**2*w3r*w1i*w3i**4*k*vgsol - 2*rhogeq*w3r*w1i**2*w3i*w2r*k**2*cs**2*rhogsol -&
        rhogeq**2*w3r**2*w2r*w3i*w1r**2*k*vgsol + rhogeq*w3r**2*w2r*w1i*w3i*k*Kdrag*vdsol -&
        rhogeq**2*w3r**2*w2r*w1i**2*w3i*k*vgsol + rhogeq**2*w3r**3*w2i**2*w1i*k*vgsol +&
        3*rhogeq*w3r**3*w2i*w1i*k*Kdrag*vgsol - rhogeq*w3r**3*w2i*w1i*k*Kdrag*vdsol +&
        rhogeq**2*w3r*w3i**2*w2i*w1i**2*k*vgsol - 2*rhogeq*w3r**2*w3i**3*w2r*w1r*rhogsol +&
        rhogeq*w3r**2*w3i*w2i**2*w1i**2*rhogsol + rhogeq*w3r*w3i**2*w2r*w1r*k*Kdrag*vdsol +&
        rhogeq*w3r**2*w1r**2*w2i**2*w3i*rhogsol + rhogeq*w3r*w3i**2*w2r*w1r*k*Kdrag*vgsol +&
        2*rhogeq*w3r**3*w3i**2*k*Kdrag*vdsol + rhogeq**2*w3r**3*w1i*w2r**2*k*vgsol +&
        rhogeq**2*w3r*w1i*w3i**2*w2r**2*k*vgsol - rhogeq**2*w3r**2*w2r**2*w3i*k*vgsol*w1r + Kdrag**2*w3i**3*w2r*k*vgsol&
        + rhogeq**2*w3r*w3i**2*w2i*w1r**2*k*vgsol - 2*rhogeq*w3r**3*w1i*w3i*k*Kdrag*vgsol -&
        2*rhogeq*w3r**3*w3i*w2i*k*Kdrag*vgsol - 2*rhogeq*w3r*w3i*w2i**2*k**2*cs**2*w1r*rhogsol +&
        rhogeq*w3r**2*w2r**2*w1i**2*w3i*rhogsol - rhogeq*w3r**2*w2r*w1i*w3i*k*Kdrag*vgsol +&
        rhogeq*w3r*w3i**4*w2i*w1r*rhogsol - vgsol*cs**2*k**3*rhogeq**2*w2i*w1r**2*w3r +&
        3*rhogeq*w3r**2*cs**4*k**4*w3i*rhogsol - 2*rhogeq*w3r*cs**4*k**4*w3i*w2r*rhogsol -&
        2*rhogeq*w3r*cs**4*k**4*w3i*w1r*rhogsol - rhogeq*w3r**3*cs**2*k**3*Kdrag*vdsol +&
        rhogeq**2*w3r*w3i**2*w1i*w2i**2*k*vgsol + 3*rhogeq*w3r*cs**2*k**3*w3i**2*Kdrag*vdsol +&
        rhogeq**2*w3r**3*w2i*w1i**2*k*vgsol - 3*rhogeq*w3r*cs**2*k**3*w3i**2*Kdrag*vgsol +&
        2*rhogeq*w3r*cs**2*k**3*w2i*w3i*Kdrag*vgsol - rhogeq**2*w3r**2*w1r*w2i**2*w3i*k*vgsol +&
        rhogeq*w3r**2*w2r**2*w3i*w1r**2*rhogsol + 2*rhogeq**2*w3r*cs**2*k**3*w1i**2*w3i*vgsol -&
        3*rhogeq**2*w3r*cs**2*k**3*w2i*w3i**2*vgsol + 4*rhogeq**2*w3r*cs**2*k**3*w3i*vgsol*w1r*w2r +&
        2*rhogeq**2*w3r*cs**2*k**3*w2r**2*w3i*vgsol - 2*rhogeq*w3i**2*w3r**2*w2r*k*Kdrag*vdsol +&
        4*rhogeq*w3i*w3r**2*w2r*w1r*k**2*cs**2*rhogsol + 2*rhogeq*w3i**2*w3r**2*w2i*k**2*cs**2*rhogsol -&
        2*rhogeq*w3i**2*w3r**2*w1i**2*rhogsol*w2i - 2*rhogeq*w3i**2*w3r**2*w1i*w2r**2*rhogsol +&
        rhogeq**2*w3r**3*cs**2*k**3*w2i*vgsol - 2*rhogeq*w3r*cs**2*k**3*w1i*w3i*Kdrag*vdsol +&
        4*rhogeq**2*w3r*cs**2*k**3*w1i*w2i*w3i*vgsol + 2*rhogeq**2*w3r*cs**2*k**3*w2i**2*w3i*vgsol +&
        2*rhogeq**2*w3r*cs**2*k**3*w3i*w1r**2*vgsol - 3*rhogeq**2*w3r**2*cs**2*k**3*w3i*w1r*vgsol -&
        2*rhogeq*w3r*cs**2*k**3*w2i*w3i*Kdrag*vdsol - 3*rhogeq**2*w3r*cs**2*k**3*w1i*w3i**2*vgsol -&
        rhogeq**2*w3r*cs**2*k**3*w1i*w2r**2*vgsol - rhogeq**2*w3r*cs**2*k**3*w1i*w2i**2*vgsol +&
        2*w3i*Kdrag**2*w3r**3*k*vdsol + 2*rhogeq*w3r*cs**2*k**3*w1i*w3i*Kdrag*vgsol +&
        rhogeq*w3r**3*cs**2*k**3*Kdrag*vgsol - 3*rhogeq**2*w3r**2*cs**2*k**3*w3i*w2r*vgsol -&
        2*rhogeq*w3i**2*w3r**2*w2i*w1r**2*rhogsol - Kdrag**2*w3i**2*w1i*w3r*k*vdsol + Kdrag**2*w3i**3*k*vgsol*w1r +&
        Kdrag*w3i**3*w1i*k**2*cs**2*rhogsol + Kdrag**2*w3i**2*w1i*w2r*k*vdsol - Kdrag**2*w3i**2*w1i*w2r*k*vgsol +&
        Kdrag*w3i**2*w1i**2*w3r*rhogeq*k*vgsol + 2*rhogeq**2*w3r**2*w3i**3*w2r*k*vgsol -&
        rhogeq*w3r**2*w1r*w2i*w3i*k*Kdrag*vgsol - rhogeq*w3i**5*w2r*w1r*rhogsol + 2*rhogeq**2*w3r**2*w3i**3*k*w1r*vgsol&
        + rhogeq**2*w3r**3*w2i*w1r**2*k*vgsol + rhogeq*w3r**5*w1i*w2r*rhogsol + 2*rhogeq*w3r**2*w3i**3*w1i*w2i*rhogsol&
        - rhogeq**2*w3r**5*w1i*k*vgsol + rhogeq*w3i**3*w2i**2*w1r**2*rhogsol - Kdrag*w3i**4*w2r*w1r*rhogsol +&
        rhogeq*w3r**2*w1r*w2i*w3i*k*Kdrag*vdsol + rhogeq*w3r**3*w2r*w1r*k*Kdrag*vdsol + rhogeq*w3r**5*w2i*w1r*rhogsol -&
        2*rhogeq**2*w3r**3*w3i**2*w2i*k*vgsol - 2*rhogeq*w3r**3*w3i**2*k*Kdrag*vgsol -&
        2*rhogeq*w3r**2*w3i**3*k**2*cs**2*rhogsol - rhogeq*w3r**5*k*Kdrag*vgsol -&
        2*rhogeq*w3r*w2r*w1r**2*w3i*k**2*cs**2*rhogsol - 2*rhogeq*w3r*w2r**2*w1r*w3i*k**2*cs**2*rhogsol -&
        rhogeq**2*w3r**5*w2i*k*vgsol - 2*rhogeq**2*w3r**3*w3i**2*w1i*k*vgsol + rhogeq*w3i**3*w2r**2*w1r**2*rhogsol -&
        Kdrag*w3i**2*w3r*w1r**2*w2r*rhogsol - Kdrag*w3i**2*w3r*w1r*w2r**2*rhogsol + Kdrag*w3i**4*w1i*w2i*rhogsol -&
        Kdrag*w3i**2*w1i**2*w3r*w2r*rhogsol + Kdrag*w3i**2*w2r*w1r*k**2*cs**2*rhogsol +&
        2*Kdrag*w3i**4*rhogeq*w2r*k*vgsol + 2*Kdrag*w3i**4*rhogeq*k*w1r*vgsol - Kdrag*w3i**4*k**2*cs**2*rhogsol -&
        Kdrag*w3i**2*w2i**2*w3r*w1r*rhogsol + Kdrag**2*w3i**2*w2i*w3r*k*vgsol - Kdrag**2*w3i**2*w2i*w1r*k*vgsol +&
        Kdrag*w3i**3*w2i*k**2*cs**2*rhogsol + Kdrag**2*w3i**2*w2i*w1r*k*vdsol - Kdrag**2*w3i**2*w2i*w3r*k*vdsol -&
        Kdrag*w3i**2*w3r*w1r*k**2*cs**2*rhogsol + 2*Kdrag*w3i**3*w2i*w3r*w1r*rhogsol - Kdrag**2*w3i**3*w2r*k*vdsol -&
        Kdrag*w3i**3*w2i*w1r**2*rhogsol - Kdrag*w3i**3*w1i**2*rhogsol*w2i - Kdrag*w3i**3*w1i*w2r**2*rhogsol +&
        2*rhogeq*w3r**3*w3i**2*w2i*w1r*rhogsol + rhogeq*w3r**3*w2r*w1r*k*Kdrag*vgsol + rhogeq*w3r**5*k*Kdrag*vdsol +&
        2*rhogeq*w3r**3*w1i*w3i**2*w2r*rhogsol - Kdrag*w3i**2*w3r*w2r*k**2*cs**2*rhogsol +&
        Kdrag*w3i**2*w3r*rhogeq*w2r**2*k*vgsol + Kdrag*w3i**2*w3r*rhogeq*w1r**2*k*vgsol - 2*Kdrag**2*w3i**3*w3r*k*vgsol&
        + 2*Kdrag**2*w3i**3*w3r*k*vdsol - Kdrag*w3i**2*rhogeq*w2r*w1r**2*k*vgsol +&
        2*Kdrag*w3i**2*w3r**2*rhogeq*w2r*k*vgsol + Kdrag**2*w3i**2*w1i*w3r*k*vgsol +&
        2*Kdrag*w3i**2*w3r**2*rhogeq*k*w1r*vgsol - Kdrag*w3i**2*w1i*w2i*k**2*cs**2*rhogsol -&
        Kdrag**2*w3i**3*k*vdsol*w1r + 2*Kdrag*w3i**3*w1i*w3r*w2r*rhogsol - Kdrag*w3i**2*w1i**2*rhogeq*w2r*k*vgsol +&
        Kdrag*w3i**2*w2i**2*w3r*rhogeq*k*vgsol - w3i*rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol +&
        w3i**2*rhogeq*k**4*cs**4*w2i*rhogsol + w3i*Kdrag*w3r**2*w1i*k**2*cs**2*rhogsol +&
        Kdrag*w3r**3*w1i**2*rhogeq*k*vgsol - Kdrag*w3r**2*w1i**2*rhogeq*w2r*k*vgsol +&
        2*w3i**2*rhogeq*k**2*cs**2*w2i*w3r*w1r*rhogsol + Kdrag**2*w3r**2*w1i*w2r*k*vdsol -&
        w3i**3*rhogeq*k**4*cs**4*rhogsol - Kdrag**2*w3r**2*w1i*w2r*k*vgsol + 2*w3i*Kdrag*w3r**3*w1i*w2r*rhogsol +&
        Kdrag**2*w3r**3*w1i*k*vgsol + w3i*Kdrag**2*w3r**2*k*vgsol*w1r - Kdrag*w3r**2*rhogeq*w2r*w1r**2*k*vgsol -&
        w3i*Kdrag*w3r**2*w1i*rhogsol*w2i**2 - Kdrag**2*w3r**3*w1i*k*vdsol + Kdrag*w3r**3*rhogeq*w2r**2*k*vgsol +&
        Kdrag*w3r**3*rhogeq*w1r**2*k*vgsol - Kdrag*w3r**3*w1r*k**2*cs**2*rhogsol - rhogeq**2*w3i**3*w2r**2*k*w1r*vgsol&
        + rhogeq**2*w3i**5*w2r*k*vgsol + rhogeq**2*w3i**5*k*w1r*vgsol + rhogeq*w3i**4*w1i*k**2*cs**2*rhogsol +&
        rhogeq*w3i**3*w1i*w2r*k*Kdrag*vdsol - rhogeq**2*w3i**3*w2r*w1r**2*k*vgsol - rhogeq*w3i**4*w2r*k*Kdrag*vdsol -&
        rhogeq*w3i**3*w2i*Kdrag*w1r*k*vgsol + rhogeq*w3i**4*w2i*k**2*cs**2*rhogsol +&
        rhogeq*w3i**3*w2i*w1r*k*Kdrag*vdsol - rhogeq*w3i**3*w1i*w2r*k*Kdrag*vgsol - rhogeq*w3i**4*w1i*rhogsol*w2i**2 +&
        rhogeq*w3i**3*w1i**2*w2r**2*rhogsol - rhogeq*w3i**4*k*Kdrag*vdsol*w1r +&
        w3r**2*rhogeq*cs**2*k**3*w2r*Kdrag*vdsol - 2*w3r**3*rhogeq*cs**2*k**2*w2i*w1r*rhogsol -&
        2*w3r**3*rhogeq*cs**2*k**2*w1i*w2r*rhogsol - w3r**2*rhogeq*cs**2*k**3*w2r*Kdrag*vgsol -&
        rhogeq**2*w3i**3*w1i**2*w2r*k*vgsol - w3r**2*rhogeq*cs**2*k**3*Kdrag*vgsol*w1r +&
        w3r**2*rhogeq*cs**2*k**3*Kdrag*vdsol*w1r - w3r**2*rhogeq*cs**4*k**4*w2i*rhogsol +&
        w3i**3*rhogeq**2*k**3*cs**2*w1r*vgsol + rhogeq*w3i**3*w1i**2*rhogsol*w2i**2 -&
        w3i**2*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol - 2*rhogeq*w3i**2*w3r**2*w1i*rhogsol*w2i**2 +&
        w3i**2*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol - w3i**2*rhogeq*k**3*cs**2*Kdrag*vdsol*w1r +&
        w3i**2*rhogeq*k**3*cs**2*Kdrag*vgsol*w1r - Kdrag*w3r**3*w1r*w2r**2*rhogsol - Kdrag*w3r**3*w1r**2*w2r*rhogsol -&
        Kdrag**2*w3r**2*w2i*w1r*k*vgsol + 2*rhogeq*w3i**2*w3r**2*w1i*k**2*cs**2*rhogsol -&
        2*rhogeq*w3i**2*w3r**2*k*Kdrag*vdsol*w1r + Kdrag**2*w3r**3*w2i*k*vgsol -&
        4*rhogeq*w3i*w3r**2*w1i*w2i*k**2*cs**2*rhogsol - w3i*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol +&
        w3i**3*rhogeq**2*k**3*cs**2*w2r*vgsol - w3r**2*rhogeq*cs**4*k**4*w1i*rhogsol - Kdrag*w3r**3*w1i**2*w2r*rhogsol&
        + Kdrag*w3r**4*w2r*w1r*rhogsol - Kdrag*w3r**3*w2i**2*w1r*rhogsol - w3i*Kdrag*w3r**2*w1i*w2r**2*rhogsol -&
        w3i*Kdrag*w3r**2*w2i*w1r**2*rhogsol + Kdrag*w3r**4*k**2*cs**2*rhogsol - Kdrag*w3r**4*w1i*w2i*rhogsol +&
        2*w3i*Kdrag*w3r**3*w2i*w1r*rhogsol + Kdrag**2*w3r**2*w2i*w1r*k*vdsol - Kdrag**2*w3r**3*w2i*k*vdsol -&
        Kdrag*w3r**3*w2r*k**2*cs**2*rhogsol - w3i*Kdrag*w3r**2*w1i**2*rhogsol*w2i - w3i*Kdrag**2*w3r**2*w2r*k*vdsol +&
        Kdrag*w3r**2*w2r*w1r*k**2*cs**2*rhogsol + w3i*Kdrag*w3r**2*w2i*k**2*cs**2*rhogsol +&
        Kdrag*w3r**3*w2i**2*rhogeq*k*vgsol - 2*w3i*Kdrag**2*w3r**3*k*vgsol + w3i*Kdrag**2*w3r**2*w2r*k*vgsol -&
        rhogeq*w3i**5*k**2*cs**2*rhogsol - w3i*Kdrag**2*w3r**2*k*vdsol*w1r - Kdrag*w3r**2*w1i*w2i*k**2*cs**2*rhogsol -&
        w3r**4*rhogeq*w2i*w1r**2*rhogsol - w3r**4*rhogeq*w1i**2*rhogsol*w2i - w3r**4*rhogeq*w1i*rhogsol*w2i**2 -&
        w3r**4*rhogeq*k*Kdrag*vdsol*w1r - w3r**4*rhogeq*w2r*k*Kdrag*vdsol + w3r**4*rhogeq*w2i*k**2*cs**2*rhogsol +&
        w3r**4*rhogeq*w1i*k**2*cs**2*rhogsol - w3r**4*rhogeq*w1i*w2r**2*rhogsol +&
        2*w3i**2*rhogeq*k**2*cs**2*w1i*w3r*w2r*rhogsol + w3i**2*rhogeq*k**4*cs**4*w1i*rhogsol -&
        rhogeq*w3i**4*w2i*w1r**2*rhogsol + rhogeq*w3i**5*w1i*w2i*rhogsol - rhogeq*w3i**4*w1i**2*rhogsol*w2i -&
        rhogeq*w3i**4*w1i*w2r**2*rhogsol - rhogeq**2*w3i**3*w2i**2*k*w1r*vgsol +&
        rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2*w3i - vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3i*w2r +&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r - vdsol*Kdrag*cs**2*k**3*rhogeq*w2r*w1r*w3r -&
        vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r*w3i - vgsol*Kdrag*k*rhogeq*w3r**2*w2r**2*w1r -&
        vgsol*Kdrag*k*rhogeq*w3i**2*w1r*w2r**2 + rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**2*w3i +&
        rhogsol*cs**2*k**2*rhogeq*w1i*w2r**2*w3r**2 - rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2r**2 +&
        rhogsol*Kdrag*w1r**2*w3r**2*w2r**2 + w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r +&
        w1i**2*rhogsol*Kdrag*w2i**2*w3r**2 - w1i**2*vgsol*cs**2*k**3*rhogeq**2*w2i*w3r -&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2 + w1i**2*rhogsol*Kdrag*w3i**2*w2i**2 +&
        w1i**2*rhogsol*Kdrag*w2r**2*w3r**2 + w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2 +&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3r**2 - w1i*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r -&
        w1i*rhogsol*cs**4*k**4*rhogeq*w3i*w2i + w1i*rhogsol*cs**4*k**4*rhogeq*w2r*w3r -&
        w1i*rhogsol*cs**2*k**2*rhogeq*w3i**2*w2i**2 + w1i*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3r**2 +&
        rhogsol*Kdrag*w1i**2*w3i**2*w2r**2 + vgsol*Kdrag*cs**2*k**3*rhogeq*w2r*w1r*w3r -&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r)/rhogeq/Kdrag/(w3i**2 + w3r**2)/( - 2*w2r*w3r + w2r**2 + w3r**2 +&
        w3i**2 + w2i**2 - 2*w3i*w2i)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r) 

  rhod2r = - ( - rhogeq*rhogsol*w2r**5*w1r*w3r - w1i**2*rhogsol*Kdrag*w2r**3*w3i - w1i*rhogsol*Kdrag*w2r**3*w3r**2 +&
        w1i*rhogeq*rhogsol*w3i*w2r**5 + w1i**2*vgsol*cs**2*k**3*rhogeq**2*w2r**2 - vdsol*Kdrag*k*rhogeq*w3i*w1r*w2r**3&
        + cs**2*k**3*rhogeq*w3i*Kdrag*vgsol*w2r**2 + rhogeq*rhogsol*w2r**4*w1r**2*w3r - rhogsol*Kdrag*w3i*w1r**2*w2r**3&
        + rhogsol*Kdrag*w2r**4*w3i*w1r - vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1i + rhogsol*cs**4*k**4*rhogeq*w2r**3 +&
        vgsol*Kdrag*k*rhogeq*w3i*w1r*w2r**3 - 2*vgsol*Kdrag*k*rhogeq*w3i*w2r**4 + vgsol*Kdrag*k*rhogeq*w3r*w2r**3*w1i -&
        rhogsol*cs**2*k**2*rhogeq*w2r**5 - cs**2*k**3*rhogeq*w3i*Kdrag*vdsol*w2r**2 + vdsol*Kdrag*k*rhogeq*w2r**4*w1i +&
        vgsol*k*rhogeq**2*w2r**3*w1r*w3r**2 - 2*vgsol*k*rhogeq**2*w3i*w2r**4*w1i + w1i*rhogsol*Kdrag*w2r**4*w3r +&
        w1i*rhogsol*cs**2*k**2*Kdrag*w2r**3 - w1i**2*rhogeq*rhogsol*w2r**3*w3i**2 - w1i*rhogsol*Kdrag*w2r**3*w3i**2 +&
        vdsol*k*Kdrag**2*w3r*w2r**2*w1r + vgsol*Kdrag*k*rhogeq*w3i*w2r**2*w1i**2 -&
        rhogsol*cs**2*k**2*Kdrag*w3i*w2r**2*w1r - cs**4*k**4*rhogeq*w3r*rhogsol*w2r**2 +&
        w1i**2*rhogeq*rhogsol*w2r**4*w3r - vgsol*cs**2*k**3*rhogeq**2*w2r**3*w1r - vdsol*k*Kdrag**2*w3r*w2r**3 +&
        w2i**4*Kdrag**2*k*vgsol - rhogeq**2*w3r**2*k*vgsol*w2r**4 - rhogeq**2*w2r**4*w3i**2*k*vgsol +&
        rhogeq**2*w3r*w2r**3*k*w1r**2*vgsol + rhogeq**2*w1i**2*w2r**3*w3r*k*vgsol + rhogeq*w3i**2*w1r*w2r**4*rhogsol +&
        rhogeq*w1r*w3r**2*w2r**4*rhogsol + rhogeq*w3i*w2r**2*w1r**2*k*Kdrag*vgsol + rhogeq*w3i*k*Kdrag*vdsol*w2r**4 -&
        2*rhogeq*w2r**4*k*Kdrag*vgsol*w1i + rhogeq*w2r**4*w3r*k**2*cs**2*rhogsol +&
        2*w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol + w2i**2*Kdrag*rhogeq*w3i*w1r**2*k*vgsol +&
        2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w2r**2*w3r - w2i**2*Kdrag*k**2*cs**2*w3i*w1r*rhogsol +&
        w2i**2*Kdrag**2*w1r*w3r*k*vdsol - w2i**3*Kdrag*rhogeq*w3i**2*k*vgsol +&
        rhogeq*w3r**2*w1i**2*w2r*k**2*cs**2*rhogsol - cs**2*k**3*rhogeq*w2i**3*Kdrag*vdsol -&
        3*w2r*rhogsol*cs**4*k**4*rhogeq*w2i**2 + cs**2*k**3*rhogeq*w2i**3*Kdrag*vgsol -&
        w2i**3*Kdrag*rhogeq*w1r**2*k*vgsol + w2i**3*Kdrag*k**2*cs**2*w1r*rhogsol + w2i**3*Kdrag**2*w3i*k*vdsol -&
        w2i**3*Kdrag**2*w3i*k*vgsol - rhogeq*w3r**2*w2r**2*w1r*k**2*cs**2*rhogsol +&
        rhogeq*w3r**2*w2i**2*k**2*cs**2*w1r*rhogsol - cs**2*k**3*rhogeq*w3i*Kdrag*vgsol*w2i**2 +&
        cs**4*k**4*rhogeq*w3r*rhogsol*w2i**2 + cs**2*k**3*rhogeq**2*w2i**3*w3i*vgsol +&
        cs**2*k**3*rhogeq*w3i*Kdrag*vdsol*w2i**2 + 2*Kdrag*w2r**3*rhogeq*k*vgsol*w2i*w1r +&
        2*rhogeq**2*w2i**2*w2r**3*w3r*k*vgsol + w2r*vgsol*k*rhogeq**2*w2i**4*w1r -&
        cs**2*k**3*rhogeq**2*w2i**2*vgsol*w1r**2 + 2*rhogeq*w2i**2*w3r*w2r**2*w1r**2*rhogsol +&
        rhogeq*w2i**4*w1r*w3r**2*rhogsol - 2*rhogeq**2*w2i**2*w1r**2*k*vgsol*w2r**2 +&
        2*rhogeq**2*w2i**3*w3i*k*vgsol*w2r**2 + rhogeq*w2i**4*w3i**2*w1r*rhogsol - w2i**4*Kdrag*w3i*w1r*rhogsol -&
        4*rhogeq*w2i**2*w2r*w1r*w3r*k**2*cs**2*rhogsol - rhogeq*w2i**5*k*Kdrag*vdsol + rhogeq*w2i**4*w3r*w1r**2*rhogsol&
        - rhogeq**2*w2i**4*w3r**2*k*vgsol - rhogeq**2*w2i**4*w3i**2*k*vgsol - 3*rhogeq*w2i**3*w3r*Kdrag*w1r*k*vgsol +&
        w2i**3*Kdrag*w3r*w1r**2*rhogsol + w2i**3*Kdrag*w1r*w3r**2*rhogsol + w2i**3*Kdrag*w3i**2*w1r*rhogsol +&
        rhogeq*w3r**2*w2r*w1r**2*k**2*cs**2*rhogsol - rhogeq**2*w3r**2*cs**2*k**3*vgsol*w1r*w2r +&
        rhogeq**2*w3r**2*cs**2*k**3*w2r**2*vgsol + rhogeq**2*w3r**2*cs**2*k**3*w1i*w2i*vgsol -&
        rhogeq**2*w3r**2*cs**2*k**3*w2i**2*vgsol - 2*rhogeq**2*w2i**4*w3r*k*w1r*vgsol +&
        w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2 - w3r*w1i*rhogsol*cs**4*k**4*rhogeq*w2i -&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r - w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2 -&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r - cs**2*k**3*rhogeq**2*w3i**2*vgsol*w2r*w1r +&
        rhogeq*w2i**2*k**2*cs**2*w3i**2*w1r*rhogsol + w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r -&
        2*w3r*rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol + w3r*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**2 +&
        w3r*rhogsol*cs**4*k**4*rhogeq*w2r*w1r - w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r -&
        w2i**2*Kdrag**2*w3r*w1r*k*vgsol - w2i**3*Kdrag*rhogeq*w3r**2*k*vgsol + cs**2*k**2*rhogeq*w2i**4*w3r*rhogsol +&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol*w3i - 2*cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol*w3i +&
        cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol*w3i - cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol*w3i +&
        2*cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol*w3i + rhogeq*w2r*k**2*cs**2*rhogsol*w3i**2*w1r**2 +&
        rhogeq*w2i**2*w1i*k*w3i**2*Kdrag*vgsol + rhogeq*w2i**2*w1i*k*Kdrag*vgsol*w3r**2 -&
        cs**2*k**3*rhogeq**2*w2i**2*w3i**2*vgsol + cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol*w3i -&
        cs**4*k**4*rhogeq*w2i*w1r*rhogsol*w3i + rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol*w3i +&
        Kdrag*w1i*w2r**2*w3i**2*rhogeq*k*vgsol + Kdrag*w1i*w2r**2*rhogeq*k*vgsol*w3r**2 +&
        rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i**2 - rhogeq*w2r*cs**4*k**4*w1i*rhogsol*w3i -&
        rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol*w3i + w2i**5*Kdrag*rhogeq*k*vgsol + 2*w2i*vgsol*Kdrag*k*rhogeq*w3r*w2r**3&
        + 2*w2r*vgsol*Kdrag*k*rhogeq*w2i**3*w3r - w1i**2*w2i**3*Kdrag*rhogeq*k*vgsol -&
        w1i**2*cs**2*k**3*rhogeq**2*w2i**2*vgsol - w1i**2*w2r*rhogeq*rhogsol*w2i**2*w3r**2 -&
        w1i**2*w2r*rhogsol*Kdrag*w3i*w2i**2 + 2*w1i**2*rhogeq*w2i**2*w3r*rhogsol*w2r**2 +&
        w1i**2*rhogeq*w2i**4*w3r*rhogsol - 2*w1i**2*w2r*rhogsol*cs**2*k**2*rhogeq*w3i*w2i -&
        w1i**2*rhogeq**2*w2i**4*k*vgsol + w1i**2*rhogeq**2*w2i**3*w3i*k*vgsol -&
        2*w1i**2*rhogeq**2*w2i**2*k*vgsol*w2r**2 - w1i**2*w2r*rhogeq*rhogsol*w2i**2*w3i**2 -&
        w1i**2*vgsol*Kdrag*k*rhogeq*w2r**2*w2i + w1i**2*vgsol*k*rhogeq**2*w2i*w3i*w2r**2 +&
        w1i**2*w2i**3*Kdrag*w3r*rhogsol + w1i**2*w2i**2*Kdrag*w3i*rhogeq*k*vgsol +&
        w1i**2*w2r*vgsol*k*rhogeq**2*w2i**2*w3r + w1i**2*rhogsol*Kdrag*w2i*w2r**2*w3r + w1i*rhogeq**2*w2i**5*k*vgsol +&
        w1i*w2i**2*Kdrag**2*w3i*k*vgsol - w1i*w2i**2*Kdrag*w3r*k**2*cs**2*rhogsol -&
        w1i*cs**2*k**3*rhogeq*w2i**2*Kdrag*vgsol - 4*w1i*rhogeq**2*w2i**2*w3i*w2r**2*k*vgsol -&
        2*w1i*rhogeq*w2i**3*w3r*k**2*cs**2*rhogsol + w1i*rhogeq**2*w2i*w2r**2*w3r**2*k*vgsol +&
        w1i*cs**2*k**3*rhogeq*w2i**2*Kdrag*vdsol + 2*w1i*cs**4*k**4*rhogeq*w2i*w2r*rhogsol -&
        w1i*w2i**3*Kdrag**2*k*vgsol + 4*w1i*rhogeq*w2r*w2i**2*w3i*k**2*cs**2*rhogsol -&
        w1i*rhogeq*w2r*w2i**2*w3r*k*Kdrag*vdsol - 2*vgsol*k*rhogeq**2*w2r**4*w3r*w1r - vgsol*k*Kdrag**2*w3r*w2r**2*w1r&
        - 2*vgsol*k*rhogeq**2*w2i**2*w2r**2*w3r**2 - 3*vgsol*cs**2*k**3*rhogeq**2*w2i*w3i*w2r**2 +&
        vdsol*Kdrag*k*rhogeq*w2r**2*w1r*w3r*w2i - vgsol*k*Kdrag**2*w2i*w3i*w2r**2 - rhogeq*rhogsol*w3i*w2r**4*w1r*w2i -&
        2*vdsol*Kdrag*k*rhogeq*w2i**3*w2r**2 + vgsol*k*Kdrag**2*w3r*w2r**3 - w1i**2*vgsol*k*rhogeq**2*w2r**4 -&
        w2r*rhogsol*Kdrag*w3i*w2i**2*w1r**2 - w2r*rhogeq*rhogsol*w2i**2*w3r**2*w1r**2 -&
        w2r*rhogeq*rhogsol*w2i**2*w3i**2*w1r**2 - w2r*rhogsol*cs**2*k**2*rhogeq*w2i**4 -&
        w2r*rhogeq*rhogsol*w2i**4*w3r*w1r + w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r -&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r - w3r*rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**2 +&
        cs**2*k**3*rhogeq**2*w3i**2*w2r**2*vgsol - cs**2*k**2*rhogeq*w3i**2*w1r*rhogsol*w2r**2 +&
        cs**2*k**3*rhogeq**2*w1i*w2i*w3i**2*vgsol + rhogsol*cs**2*k**2*rhogeq*w2r**4*w1r +&
        cs**2*k**3*rhogeq**2*w1r**2*w2r**2*vgsol - cs**4*k**4*rhogeq*w2r**2*rhogsol*w1r -&
        vdsol*k*Kdrag**2*w3i*w2r**2*w1i - rhogeq*rhogsol*w3r**2*w2r**3*w1r**2 - vdsol*k*Kdrag**2*w2r**3*w1r -&
        rhogeq*rhogsol*w3i**2*w1r**2*w2r**3 + vgsol*k*Kdrag**2*w2r**3*w1r - vdsol*Kdrag*k*rhogeq*w3r*w2r**3*w1i +&
        vgsol*k*Kdrag**2*w3i*w2r**2*w1i + vgsol*k*rhogeq**2*w2r**3*w3i**2*w1r + vgsol*k*rhogeq**2*w3r*w2r**5 +&
        vgsol*k*rhogeq**2*w1r*w2r**5 - rhogsol*cs**2*k**2*Kdrag*w2r**2*w3r*w1i + rhogsol*cs**2*k**2*Kdrag*w2r**3*w3i +&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1i - w1i**2*rhogeq*rhogsol*w3r**2*w2r**3 -&
        vgsol*cs**2*k**3*rhogeq**2*w3r*w2r**3 - vgsol*k*rhogeq**2*w2r**4*w1r**2 -&
        4*vgsol*k*rhogeq**2*w2i**2*w3r*w2r**2*w1r - w2r*vdsol*Kdrag*k*rhogeq*w2i**2*w3i*w1r +&
        w2r*vgsol*k*rhogeq**2*w2i**2*w1r*w3r**2 - 2*w2r*rhogsol*Kdrag*w2i**3*w3r*w1r +&
        3*w2r*vgsol*cs**2*k**3*rhogeq**2*w2i**2*w1r - 2*rhogeq*w2i**3*k**2*cs**2*w3i*w1r*rhogsol +&
        rhogeq**2*w2i**5*w3i*k*vgsol + rhogeq*w2i**3*w1r*w3r*k*Kdrag*vdsol + rhogeq**2*w2i**3*w3i*w1r**2*k*vgsol +&
        rhogeq**2*w2i**4*w2r*w3r*k*vgsol + 2*rhogeq*rhogsol*w2i**2*w3i**2*w2r**2*w1r -&
        2*rhogsol*Kdrag*w2i*w2r**3*w3r*w1r + 2*rhogeq*rhogsol*w2i**2*w3r**2*w2r**2*w1r +&
        vdsol*k*Kdrag**2*w2i*w3i*w2r**2 + rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2*w1r -&
        3*vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w2i + rhogsol*Kdrag*w2i*w3i**2*w2r**2*w1r +&
        rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2*w3r + rhogsol*Kdrag*w2r**2*w3r*w2i*w1r**2 +&
        rhogeq**2*w2i*w3i*k*vgsol*w2r**4 + w2i**3*Kdrag*w3r*k**2*cs**2*rhogsol +&
        3*cs**2*k**3*rhogeq**2*w2i**2*w2r*w3r*vgsol - rhogeq*w2i**5*w3i*w1r*rhogsol +&
        cs**4*k**4*rhogeq*w2i**2*w1r*rhogsol + 2*vgsol*k*rhogeq**2*w2i**2*w2r**3*w1r +&
        2*vgsol*Kdrag*k*rhogeq*w2i**3*w2r**2 - vgsol*Kdrag*k*rhogeq*w3i**2*w2r**2*w2i -&
        3*vgsol*Kdrag*k*rhogeq*w2r**2*w1r*w3r*w2i + vgsol*Kdrag*k*rhogeq*w2i*w2r**4 -&
        vgsol*Kdrag*k*rhogeq*w2i*w2r**2*w1r**2 - 2*vgsol*Kdrag*k*rhogeq*w2i**2*w3i*w2r**2 +&
        rhogsol*Kdrag*w2r**2*w3r**2*w1r*w2i + 2*w1i*vgsol*k*rhogeq**2*w2i**3*w2r**2 -&
        2*w1i*rhogeq*rhogsol*w2i**3*w2r**2*w3r - 3*w1i*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r**2 +&
        2*w1i*rhogsol*Kdrag*w2r**3*w2i*w3i + 2*w1i*rhogeq*rhogsol*w2i**2*w2r**3*w3i -&
        w1i*rhogeq*w2i**3*k*w3i*Kdrag*vgsol - 2*w1i*w2r*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2 +&
        w1i*vgsol*k*rhogeq**2*w2i*w2r**4 + w1i*w2r*rhogeq*rhogsol*w3i*w2i**4 + w1i*rhogeq**2*w2i**3*w3i**2*k*vgsol +&
        w1i*rhogeq**2*w2i**3*w3r**2*k*vgsol + w1i*rhogeq**2*w2i*w2r**2*w3i**2*k*vgsol + w1i*rhogeq*w2i**4*k*Kdrag*vdsol&
        - 2*w1i*rhogeq*w2i**2*w2r**2*k*Kdrag*vgsol - w1i*rhogeq*w2i**3*w3i*k*Kdrag*vdsol -&
        2*w1i*rhogeq**2*w2i**4*w3i*k*vgsol + 2*w1i*rhogeq*w2i*w2r**2*w3r*k**2*cs**2*rhogsol -&
        w1i*rhogeq*w2i*k*w3i*Kdrag*vdsol*w2r**2 + 2*w1i*rhogeq*w2i**2*w2r**2*k*Kdrag*vdsol -&
        w1i*w2i*Kdrag**2*w2r**2*k*vgsol - w1i*w2i**2*Kdrag**2*w3i*k*vdsol + w1i*w2i**2*Kdrag*w2r*k**2*cs**2*rhogsol +&
        w1i*w2i*Kdrag**2*w2r**2*k*vdsol - w1i*rhogeq*w2i**5*w3r*rhogsol - w1i*w2i**4*Kdrag*w3r*rhogsol -&
        w1i*w2r*rhogsol*Kdrag*w2i**2*w3r**2 + 2*w1i*w2r*rhogsol*Kdrag*w3i*w2i**3 +&
        w1i*cs**2*k**3*rhogeq**2*w2i**3*vgsol - w1i*w2r*rhogsol*Kdrag*w3i**2*w2i**2 +&
        w1i*w2r*vgsol*Kdrag*k*rhogeq*w2i**2*w3r - 2*w1i*w2r*rhogsol*cs**2*k**2*rhogeq*w2i*w3r**2 -&
        w1i*rhogeq*rhogsol*w2r**4*w3r*w2i - w1i*vgsol*Kdrag*k*rhogeq*w2i*w3i*w2r**2 + w1i*w2i**3*Kdrag**2*k*vdsol +&
        2*rhogeq*w2i**2*k**2*cs**2*w1r*rhogsol*w2r**2 + rhogeq*w2i**4*w3i*k*Kdrag*vdsol -&
        rhogeq**2*w2i**4*w1r**2*k*vgsol + rhogeq*w2i**4*k**2*cs**2*w1r*rhogsol -&
        2*vgsol*k*rhogeq**2*w2i**2*w3i**2*w2r**2 - 2*rhogsol*cs**2*k**2*Kdrag*w2i*w2r**3 -&
        2*rhogeq*rhogsol*w2r**3*w1r*w2i**2*w3r + 2*vdsol*Kdrag*k*rhogeq*w2i**2*w3i*w2r**2 -&
        2*w2r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r + w2r*vgsol*k*rhogeq**2*w2i**2*w1r*w3i**2 -&
        2*w2r*rhogsol*cs**2*k**2*rhogeq*w2i*w3i*w1r**2 + w2r*vgsol*k*Kdrag**2*w2i**2*w1r -&
        2*w2r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r + 2*w2r*vgsol*Kdrag*k*rhogeq*w2i**3*w1r +&
        w2r*vgsol*k*rhogeq**2*w1r**2*w2i**2*w3r + 2*w2r*rhogsol*cs**4*k**4*rhogeq*w3i*w2i +&
        w2r*vgsol*Kdrag*k*rhogeq*w2i**2*w3i*w1r + 2*w2r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r +&
        w2r*vgsol*k*Kdrag**2*w2i**2*w3r - w2r*vdsol*k*Kdrag**2*w2i**2*w1r - 2*w2r*rhogsol*cs**2*k**2*Kdrag*w2i**3 -&
        w2r*vdsol*k*Kdrag**2*w2i**2*w3r + 2*w2r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r +&
        w2r*rhogsol*cs**2*k**2*Kdrag*w3i*w2i**2 - 2*rhogeq*rhogsol*w2i**3*w3i*w2r**2*w1r -&
        vdsol*Kdrag*k*rhogeq*w2i*w2r**4 + vgsol*k*rhogeq**2*w1r**2*w2i*w3i*w2r**2 -&
        vgsol*Kdrag*k*rhogeq*w2r**2*w3r**2*w2i + 2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i*w2r**2*w1r -&
        2*rhogsol*cs**2*k**2*rhogeq*w2r**3*w2i**2 + 3*vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w2i -&
        w2i**4*Kdrag**2*k*vdsol + vdsol*k*Kdrag**2*w2r**4 - vgsol*k*Kdrag**2*w2r**4)*rhodeq/(w2i**2 + w1r**2 -&
        2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r)/Kdrag/(w2i**2 + w2r**2)/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 +&
        w2i**2 - 2*w3i*w2i)/rhogeq 

  rhod2i = - rhodeq*(rhogeq*rhogsol*w2r**4*w2i*w1r*w3r + 2*rhogeq*rhogsol*w2r**2*w1r*w2i**3*w3r -&
        rhogeq*rhogsol*w2i**4*w3i*w2r*w1r - rhogeq*rhogsol*w3i**2*w2i*w1r**2*w2r**2 -&
        rhogeq*rhogsol*w2i*w3r**2*w2r**2*w1r**2 + vdsol*k*Kdrag**2*w2i*w2r**2*w1r + vdsol*k*Kdrag**2*w2i*w3r*w2r**2 +&
        3*vgsol*cs**2*k**3*rhogeq**2*w2i*w3r*w2r**2 + 3*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r**2*w1r -&
        vdsol*k*Kdrag**2*w2r**2*w3i*w1r + 3*vgsol*cs**2*k**3*rhogeq**2*w2i**2*w3i*w2r - vgsol*k*Kdrag**2*w2i**3*w1r +&
        2*vdsol*Kdrag*k*rhogeq*w2i**2*w3r*w2r**2 + vdsol*Kdrag*k*rhogeq*w2i**4*w1r + vdsol*Kdrag*k*rhogeq*w1r*w2r**4 -&
        vdsol*Kdrag*k*rhogeq*w2r**5 - vdsol*Kdrag*k*rhogeq*w2r*w1r*w3r*w2i**2 - vdsol*Kdrag*k*rhogeq*w2i**3*w3i*w1r -&
        vdsol*Kdrag*k*rhogeq*w2i*w3i*w1r*w2r**2 + vdsol*Kdrag*k*rhogeq*w2i**4*w3r - vgsol*k*Kdrag**2*w2i*w2r**2*w1r +&
        vgsol*k*Kdrag**2*w2r**2*w3i*w1r - vgsol*k*Kdrag**2*w2i*w3r*w2r**2 - vgsol*k*Kdrag**2*w2i**2*w3i*w2r +&
        vgsol*k*Kdrag**2*w2i**2*w3i*w1r - 2*rhogeq*rhogsol*w3i*w2r**3*w1r*w2i**2 +&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w3r - vdsol*Kdrag*cs**2*k**3*rhogeq*w3i*w2i*w1r +&
        2*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i*w2r + vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**3 +&
        vgsol*k*rhogeq**2*w1r**2*w2i**3*w3r + vgsol*k*rhogeq**2*w1r**2*w2r**2*w2i*w3r +&
        vgsol*k*rhogeq**2*w2i*w2r**2*w3i**2*w1r + vgsol*k*rhogeq**2*w2i*w2r**2*w1r*w3r**2 +&
        vgsol*k*rhogeq**2*w2i**3*w1r*w3r**2 - vgsol*k*rhogeq**2*w1r**2*w2i**2*w3i*w2r +&
        vgsol*k*rhogeq**2*w2i**4*w3i*w2r + vgsol*k*rhogeq**2*w2i**3*w1r*w3i**2 - vgsol*k*rhogeq**2*w2i*w1r*w2r**4 -&
        vgsol*k*rhogeq**2*w2i*w3r*w2r**4 - 2*vgsol*k*rhogeq**2*w2i**3*w2r**2*w3r -&
        2*vgsol*k*rhogeq**2*w2i**3*w2r**2*w1r + 2*vgsol*k*rhogeq**2*w2i**2*w3i*w2r**3 +&
        2*vdsol*Kdrag*k*rhogeq*w2i**2*w2r**2*w1r - vdsol*Kdrag*k*rhogeq*w2i**4*w2r -&
        2*vdsol*Kdrag*k*rhogeq*w2i**2*w2r**3 + vdsol*Kdrag*k*rhogeq*w3r*w2r**4 -&
        rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2*w1r**2 + rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i*w1r**2 -&
        rhogsol*cs**2*k**2*rhogeq*w2i*w1r**2*w3r**2 + 2*rhogsol*cs**2*k**2*rhogeq*w2r**3*w3i*w1r -&
        2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w2r**2*w3i + rhogsol*cs**2*k**2*rhogeq*w2i*w2r**4 -&
        rhogsol*cs**2*k**2*rhogeq*w3i*w2i**4 - rhogsol*cs**2*k**2*rhogeq*w2r**4*w3i +&
        2*rhogsol*cs**2*k**2*rhogeq*w1r*w2r*w2i*w3r**2 + 2*rhogsol*cs**2*k**2*rhogeq*w2i*w2r*w1r**2*w3r +&
        2*rhogsol*cs**2*k**2*rhogeq*w2r**2*w2i**3 - 2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i*w2r*w1r +&
        2*rhogsol*cs**2*k**2*rhogeq*w3i**2*w1r*w2r*w2i - 4*rhogsol*cs**2*k**2*rhogeq*w2r**2*w2i*w3r*w1r -&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1r - 3*vdsol*Kdrag*cs**2*k**3*rhogeq*w2r*w2i**2 +&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w1r - vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w3r -&
        vgsol*Kdrag*k*rhogeq*w2r*w3r**2*w2i**2 - 2*vgsol*Kdrag*k*rhogeq*w2i**4*w1r + vgsol*Kdrag*k*rhogeq*w2i**4*w2r +&
        2*vgsol*Kdrag*k*rhogeq*w2i**2*w2r**3 - 2*vgsol*Kdrag*k*rhogeq*w2i**2*w3r*w2r**2 +&
        2*vgsol*Kdrag*k*rhogeq*w2i*w3i*w2r**3 + 2*vgsol*Kdrag*k*rhogeq*w2i**3*w3i*w2r + vgsol*Kdrag*k*rhogeq*w2r**5 +&
        vgsol*Kdrag*k*rhogeq*w2i**3*w3i*w1r - vgsol*Kdrag*k*rhogeq*w2r*w1r*w3r*w2i**2 -&
        vgsol*Kdrag*k*rhogeq*w3i**2*w2r*w2i**2 + vgsol*Kdrag*k*rhogeq*w1r*w3r**2*w2i**2 -&
        2*vgsol*Kdrag*k*rhogeq*w2i**4*w3r + vgsol*Kdrag*k*rhogeq*w3i**2*w1r*w2i**2 +&
        vgsol*Kdrag*k*rhogeq*w2i*w3i*w1r*w2r**2 - 2*vgsol*Kdrag*k*rhogeq*w2i**2*w2r**2*w1r -&
        vgsol*Kdrag*k*rhogeq*w2i**2*w2r*w1r**2 + vgsol*Kdrag*k*rhogeq*w2i**2*w3r*w1r**2 -&
        rhogsol*cs**4*k**4*rhogeq*w3i*w2i**2 + rhogsol*cs**4*k**4*rhogeq*w2r**2*w3i -&
        3*rhogsol*cs**4*k**4*rhogeq*w2i*w2r**2 + 2*rhogsol*cs**4*k**4*rhogeq*w2i*w2r*w3r -&
        rhogsol*cs**4*k**4*rhogeq*w2i*w3r*w1r + 2*rhogsol*cs**4*k**4*rhogeq*w2i*w2r*w1r -&
        rhogsol*cs**4*k**4*rhogeq*w3i*w2r*w1r + rhogsol*cs**2*k**2*Kdrag*w2i**4 - rhogsol*cs**2*k**2*Kdrag*w2r**4 +&
        rhogsol*cs**2*k**2*Kdrag*w2i**2*w2r*w1r - rhogsol*cs**2*k**2*Kdrag*w2r**2*w2i*w3i -&
        rhogsol*cs**2*k**2*Kdrag*w1r*w2i**2*w3r + rhogsol*cs**2*k**2*Kdrag*w2i**2*w2r*w3r +&
        rhogsol*cs**2*k**2*Kdrag*w3r*w2r**3 + rhogsol*cs**2*k**2*Kdrag*w2r**3*w1r - rhogsol*cs**2*k**2*Kdrag*w3i*w2i**3&
        + vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1r + rhogsol*Kdrag*w2r**3*w3i**2*w1r +&
        rhogsol*Kdrag*w2i**2*w3i**2*w2r*w1r + 3*vgsol*Kdrag*cs**2*k**3*rhogeq*w2r*w2i**2 -&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w1r + vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w3r -&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w3r - 2*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i*w2r +&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w3i*w2i*w1r - vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**3 +&
        rhogsol*Kdrag*w2r*w3r*w2i**2*w1r**2 + rhogsol*Kdrag*w2r*w3r**2*w1r*w2i**2 + rhogsol*Kdrag*w2r**3*w3r*w1r**2 +&
        rhogsol*Kdrag*w2r**3*w3r**2*w1r + rhogsol*Kdrag*w3i*w2i*w1r**2*w2r**2 + rhogsol*Kdrag*w3i*w2i**3*w1r**2 -&
        rhogsol*Kdrag*w3i**2*w2i**2*w1r**2 + 2*vgsol*k*Kdrag**2*w2i**3*w2r - 2*rhogsol*Kdrag*w2i*w2r**3*w3i*w1r -&
        2*rhogsol*Kdrag*w2i**3*w3i*w2r*w1r + vdsol*k*Kdrag**2*w2i**3*w1r - 2*vdsol*k*Kdrag**2*w2i**3*w2r +&
        vdsol*k*Kdrag**2*w2i**3*w3r - vgsol*k*Kdrag**2*w2i**3*w3r + rhogsol*cs**4*k**4*rhogeq*w2i**3 +&
        vdsol*k*Kdrag**2*w3i*w2r**3 - rhogsol*Kdrag*w2r**4*w1r*w3r - vgsol*cs**2*k**3*rhogeq**2*w2i**3*w3r -&
        vgsol*cs**2*k**3*rhogeq**2*w2i**3*w1r - rhogsol*Kdrag*w2i**2*w1r**2*w3r**2 + rhogsol*Kdrag*w2i**4*w3r*w1r -&
        rhogsol*Kdrag*w2r**2*w1r**2*w3i**2 - 2*vdsol*k*Kdrag**2*w2i*w2r**3 - vgsol*cs**2*k**3*rhogeq**2*w3i*w2r**3 +&
        rhogeq*rhogsol*w1r**2*w2i**4*w3i - 2*vgsol*cs**2*k**3*rhogeq**2*w2r*w3r**2*w2i +&
        vgsol*cs**2*k**3*rhogeq**2*w2i*w1r*w3r**2 + vgsol*cs**2*k**3*rhogeq**2*w2i*w1r**2*w3r -&
        2*vgsol*cs**2*k**3*rhogeq**2*w3i**2*w2i*w2r + vgsol*cs**2*k**3*rhogeq**2*w2i*w1r*w3i**2 -&
        2*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r*w1r**2 - 4*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r*w3r*w1r -&
        rhogeq*rhogsol*w3i**2*w2i**3*w1r**2 - rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2*w3i +&
        vgsol*cs**2*k**3*rhogeq**2*w1i*w2r*w3r**2 + vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3i*w2r +&
        vgsol*cs**2*k**3*rhogeq**2*w1i*w2r*w3i**2 - vgsol*k*rhogeq**2*w1i*w3i**2*w2r**3 -&
        vgsol*k*rhogeq**2*w1i*w2r**3*w3r**2 - vgsol*k*rhogeq**2*w3i*w2r**3*w1r**2 - vgsol*k*rhogeq**2*w1i**2*w3i*w2r**3&
        + vdsol*Kdrag*k*rhogeq*w3i*w1i*w2r**3 - vdsol*Kdrag*k*rhogeq*w2r**3*w1r*w3r -&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r + vdsol*Kdrag*cs**2*k**3*rhogeq*w2r*w1r*w3r +&
        vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r*w3i + vgsol*Kdrag*k*rhogeq*w3r**2*w2r**2*w1r +&
        vgsol*Kdrag*k*rhogeq*w3i**2*w1r*w2r**2 - vgsol*Kdrag*k*rhogeq*w1i**2*w2r**3 -&
        vgsol*Kdrag*k*rhogeq*w2r**3*w1r**2 + vgsol*Kdrag*k*rhogeq*w2r**2*w1r**2*w3r -&
        vgsol*Kdrag*k*rhogeq*w2r**3*w1r*w3r + vgsol*Kdrag*k*rhogeq*w1i**2*w3r*w2r**2 -&
        rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**2*w3i - rhogsol*cs**2*k**2*rhogeq*w1i*w2r**2*w3r**2 -&
        rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2r**2 + rhogsol*cs**2*k**2*Kdrag*w1i*w3i*w2r**2 -&
        rhogsol*cs**2*k**2*Kdrag*w3r*w1r*w2r**2 - rhogsol*Kdrag*w1r**2*w3r**2*w2r**2 + vdsol*k*Kdrag**2*w2i**2*w3i*w2r&
        - vdsol*k*Kdrag**2*w2i**2*w3i*w1r + 2*rhogeq*rhogsol*w2i**2*w2r**2*w1r**2*w3i + w1i*vdsol*k*Kdrag**2*w2r**3 -&
        w1i*vdsol*Kdrag*k*rhogeq*w2i*w3r*w2r**2 + w1i*vdsol*Kdrag*k*rhogeq*w2i**2*w3i*w2r -&
        w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r + 2*w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w2r -&
        w1i*vgsol*k*rhogeq**2*w2i**2*w2r*w3r**2 - w1i*rhogeq*rhogsol*w3i*w2i**5 + w1i*rhogeq*rhogsol*w2i**4*w3i**2 -&
        w1i*vgsol*k*Kdrag**2*w2r**3 - rhogeq*rhogsol*w3i*w1r*w2r**5 - w1i**2*rhogeq*rhogsol*w2i*w3r**2*w2r**2 -&
        w1i**2*rhogeq*rhogsol*w2i*w2r**2*w3i**2 + 2*w1i**2*rhogeq*rhogsol*w2i**2*w2r**2*w3i +&
        w1i**2*rhogeq*rhogsol*w2r**4*w3i - w1i**2*rhogeq*rhogsol*w2i**3*w3i**2 -&
        2*w1i**2*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r + w1i**2*rhogsol*Kdrag*w3r*w2r**3 -&
        w1i**2*rhogsol*Kdrag*w2i**2*w3r**2 - w1i**2*rhogeq*rhogsol*w2i**3*w3r**2 + w1i**2*rhogsol*Kdrag*w2r**2*w2i*w3i&
        + w1i**2*vgsol*cs**2*k**3*rhogeq**2*w2i*w3r + rhogeq*rhogsol*w2i**5*w3r*w1r +&
        w1i**2*vgsol*Kdrag*k*rhogeq*w2i**2*w3r - w1i**2*vgsol*Kdrag*k*rhogeq*w2r*w2i**2 -&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2 - w1i**2*rhogsol*Kdrag*w3i**2*w2i**2 -&
        w1i**2*rhogsol*Kdrag*w2r**2*w3r**2 + w1i**2*rhogsol*Kdrag*w2i**2*w2r*w3r +&
        2*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w2r*w3r + w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2 -&
        w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3r**2 + w1i**2*vgsol*k*rhogeq**2*w2i**3*w3r +&
        w1i**2*vgsol*k*rhogeq**2*w2r**2*w2i*w3r - w1i**2*vgsol*k*rhogeq**2*w2i**2*w3i*w2r +&
        w1i**2*rhogsol*Kdrag*w3i*w2i**3 + w1i*rhogeq*rhogsol*w2r**4*w3i**2 + w1i**2*rhogeq*rhogsol*w3i*w2i**4 -&
        2*w1i*rhogeq*rhogsol*w2r**3*w3r*w2i**2 - w1i*vgsol*k*Kdrag**2*w2r*w2i**2 + w1i*vgsol*k*Kdrag**2*w2i**2*w3r -&
        w1i*rhogeq*rhogsol*w2i**4*w2r*w3r + w1i*vgsol*k*Kdrag**2*w2r**2*w3r + w1i*vgsol*k*rhogeq**2*w2i**4*w2r -&
        w1i*rhogsol*cs**2*k**2*rhogeq*w2i**4 - w1i*rhogsol*cs**2*k**2*rhogeq*w2r**4 +&
        2*w1i*rhogeq*rhogsol*w2r**2*w2i**2*w3r**2 - 2*w1i*rhogeq*rhogsol*w3i*w2i**3*w2r**2 -&
        w1i*rhogeq*rhogsol*w3i*w2r**4*w2i - 4*w1i*w3i*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r -&
        w1i*vdsol*k*Kdrag**2*w2i**2*w3r + 3*w1i*vgsol*cs**2*k**3*rhogeq**2*w2r*w2i**2 -&
        2*w1i*rhogsol*Kdrag*w2i**3*w2r*w3r + w1i*rhogsol*Kdrag*w2r**2*w2i*w3r**2 - 2*w1i*rhogsol*Kdrag*w2i*w2r**3*w3r +&
        w1i*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r - w1i*rhogsol*cs**4*k**4*rhogeq*w2i**2 +&
        w1i*rhogeq*rhogsol*w2i**4*w3r**2 + w1i*rhogsol*cs**2*k**2*Kdrag*w3i*w2i**2 -&
        w1i*rhogsol*cs**2*k**2*Kdrag*w2i**3 + w1i*rhogeq*rhogsol*w2r**4*w3r**2 + w1i*rhogsol*cs**4*k**4*rhogeq*w3i*w2i&
        - w1i*rhogsol*cs**4*k**4*rhogeq*w2r*w3r + w1i*rhogsol*Kdrag*w3i*w2r**4 + 2*w1i*vgsol*Kdrag*k*rhogeq*w2i**3*w2r&
        + 2*w1i*vgsol*Kdrag*k*rhogeq*w2i*w2r**3 + w1i*vgsol*Kdrag*k*rhogeq*w2i*w3r*w2r**2 -&
        3*w1i*vgsol*Kdrag*k*rhogeq*w2i**2*w3i*w2r - w1i*rhogsol*Kdrag*w3i*w2i**4 + w1i*vgsol*Kdrag*k*rhogeq*w2i**3*w3r&
        - w1i*rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2 + w1i*rhogsol*cs**4*k**4*rhogeq*w2r**2 -&
        2*w1i*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w2r + w1i*rhogsol*Kdrag*w2r**2*w2i*w3i**2 -&
        w1i*vgsol*cs**2*k**3*rhogeq**2*w2r**3 - w1i*vdsol*k*Kdrag**2*w2r**2*w3r + w1i*vdsol*k*Kdrag**2*w2r*w2i**2 +&
        w1i*vgsol*k*rhogeq**2*w2r**5 + w1i*rhogsol*Kdrag*w2i**3*w3r**2 + 4*w1i*rhogsol*cs**2*k**2*rhogeq*w2r**2*w2i*w3i&
        + w1i*rhogsol*Kdrag*w2i**3*w3i**2 - w1i*vgsol*k*rhogeq**2*w2i**2*w3i**2*w2r -&
        2*w1i*rhogsol*cs**2*k**2*rhogeq*w2i**2*w2r*w3r + 2*w1i*rhogsol*cs**2*k**2*rhogeq*w3r*w2r**3 -&
        2*w1i*rhogsol*cs**2*k**2*rhogeq*w2i**2*w2r**2 + w1i*rhogsol*cs**2*k**2*rhogeq*w3i**2*w2i**2 +&
        w1i*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3r**2 - w1i*rhogeq*rhogsol*w2r**5*w3r +&
        2*w1i*rhogeq*rhogsol*w3i**2*w2r**2*w2i**2 - w1i*vdsol*Kdrag*k*rhogeq*w2i**3*w3r +&
        2*w1i*vgsol*k*rhogeq**2*w2i**2*w2r**3 - rhogeq*rhogsol*w1r**2*w3r**2*w2i**3 + rhogeq*rhogsol*w1r**2*w2r**4*w3i&
        + 2*vgsol*k*Kdrag**2*w2i*w2r**3 - vgsol*k*Kdrag**2*w3i*w2r**3 + vgsol*k*rhogeq**2*w2r**5*w3i -&
        vgsol*k*rhogeq**2*w2i**5*w3r - vgsol*k*rhogeq**2*w2i**5*w1r + rhogsol*cs**2*k**2*rhogeq*w2i**5 -&
        3*vgsol*Kdrag*k*rhogeq*w3i*w1i*w2r**3 - vgsol*Kdrag*k*rhogeq*w3i**2*w2r**3 - vgsol*Kdrag*k*rhogeq*w2r**3*w3r**2&
        - rhogsol*Kdrag*w1i**2*w3i**2*w2r**2 - vgsol*Kdrag*cs**2*k**3*rhogeq*w2r*w1r*w3r +&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r)/rhogeq/( - 2*w2r*w3r + w2r**2 + w3r**2 + w3i**2 + w2i**2 -&
        2*w3i*w2i)/(w2i**2 + w2r**2)/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r)/Kdrag 

  rhod1r = - rhodeq*( - rhogeq*w3r**3*w2r**3*w1r**3*rhogsol - rhogeq*w1i**3*w3i**3*w2r**3*rhogsol +&
        rhogsol*cs**4*k**4*rhogeq*w2r**3*w3r**2 + vgsol*Kdrag*k*rhogeq*w3i**3*w2r**4 +&
        vgsol*Kdrag*k*rhogeq*w3i*w2r**4*w3r**2 - 3*vgsol*Kdrag*k*rhogeq*w3r*w2r**3*w1i*w3i**2 -&
        3*vgsol*Kdrag*k*rhogeq*w3r**3*w2r**3*w1i - rhogsol*Kdrag*w2r**4*w3i**3*w1r -&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1i*w3i**2 - vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1i*w3r**2 -&
        rhogsol*cs**4*k**4*rhogeq*w2r**3*w3i**2 - vdsol*Kdrag*k*rhogeq*w3i*w1r*w2r**3*w3r**2 -&
        cs**2*k**3*rhogeq*w3i*Kdrag*vgsol*w2r**2*w3r**2 + w1i*rhogsol*Kdrag*w2r**3*w3r**4 +&
        2*w1i**2*vgsol*cs**2*k**3*rhogeq**2*w2r**2*w3i**2 + 2*w1i*rhogsol*Kdrag*w2r**3*w3r**2*w3i**2 -&
        vdsol*Kdrag*k*rhogeq*w3i**3*w1r*w2r**3 - cs**2*k**3*rhogeq*w3i**3*Kdrag*vgsol*w2r**2 +&
        rhogeq*rhogsol*w2r**4*w1r**2*w3r*w3i**2 + rhogeq*rhogsol*w2r**4*w1r**2*w3r**3 +&
        rhogsol*Kdrag*w3i**3*w1r**2*w2r**3 + rhogsol*Kdrag*w3i*w1r**2*w2r**3*w3r**2 -&
        rhogsol*Kdrag*w2r**4*w3i*w1r*w3r**2 - w1i**2*rhogsol*Kdrag*w2r**3*w3i**3 -&
        w1i**2*rhogsol*Kdrag*w2r**3*w3i*w3r**2 - w2i**4*Kdrag*w3i**3*w1r*rhogsol - w2i**4*Kdrag*w3i*w1r*rhogsol*w3r**2&
        - rhogeq*w2i**4*w3i**4*w1r*rhogsol - rhogeq*w2i**4*w1r*w3r**4*rhogsol + rhogeq*w2i**4*w3r*w1r**2*rhogsol*w3i**2&
        + rhogeq**2*w2i**4*w3r**4*k*vgsol + rhogeq**2*w2i**4*w3i**4*k*vgsol -&
        cs**2*k**3*rhogeq*w3i**3*Kdrag*vdsol*w2i**2 - cs**2*k**3*rhogeq*w3i*Kdrag*vdsol*w2i**2*w3r**2 +&
        cs**2*k**3*rhogeq*w3i**3*Kdrag*vgsol*w2i**2 + cs**2*k**3*rhogeq*w3i*Kdrag*vgsol*w2i**2*w3r**2 -&
        cs**4*k**4*rhogeq*w3r*rhogsol*w2i**2*w3i**2 - cs**4*k**4*rhogeq*w3r**3*rhogsol*w2i**2 +&
        cs**2*k**3*rhogeq**2*w2i**3*w3i**3*vgsol - rhogeq*w3r**4*w2i**2*k**2*cs**2*w1r*rhogsol +&
        cs**2*k**3*rhogeq**2*w2i**3*w3i*vgsol*w3r**2 - 2*cs**2*k**3*rhogeq**2*w2i**2*vgsol*w1r**2*w3r**2 +&
        2*rhogeq*w2i**2*w3r*w2r**2*w1r**2*rhogsol*w3i**2 + 2*rhogeq*w2i**2*w3r**3*w2r**2*w1r**2*rhogsol -&
        2*rhogeq*w2i**4*w1r*w3r**2*rhogsol*w3i**2 + w2i**3*Kdrag*k**2*cs**2*w1r*rhogsol*w3r**2 +&
        w2i**3*Kdrag**2*w3i*k*vgsol*w3r**2 + 2*rhogeq*w3r**2*w2r**2*w1r*k**2*cs**2*rhogsol*w3i**2 +&
        rhogeq*w3r**4*w2r**2*w1r*k**2*cs**2*rhogsol - 2*rhogeq*w3r**2*w2i**2*k**2*cs**2*w1r*rhogsol*w3i**2 -&
        2*w2i**3*Kdrag*rhogeq*w1r**2*k*vgsol*w3i**2 - 2*w2i**3*Kdrag*rhogeq*w1r**2*k*vgsol*w3r**2 +&
        w2i**3*Kdrag*k**2*cs**2*w1r*rhogsol*w3i**2 - w2i**3*Kdrag**2*w3i*k*vdsol*w3r**2 -&
        cs**2*k**3*rhogeq*w2i**3*Kdrag*vdsol*w3i**2 + cs**2*k**3*rhogeq*w2i**3*Kdrag*vdsol*w3r**2 -&
        w2r*rhogsol*cs**4*k**4*rhogeq*w2i**2*w3i**2 - w2i**3*Kdrag**2*w3i**3*k*vdsol -&
        2*rhogeq*w3r**2*w1i**2*w2r*k**2*cs**2*rhogsol*w3i**2 + w2r*rhogsol*cs**4*k**4*rhogeq*w2i**2*w3r**2 +&
        cs**2*k**3*rhogeq*w2i**3*Kdrag*vgsol*w3i**2 - cs**2*k**3*rhogeq*w2i**3*Kdrag*vgsol*w3r**2 +&
        w2i**3*Kdrag**2*w3i**3*k*vgsol + w2i**2*Kdrag*k**2*cs**2*w3i**3*w1r*rhogsol +&
        2*w3r**3*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol - 2*w2i**2*Kdrag*rhogeq*w3i**3*w1r**2*k*vgsol -&
        2*w2i**2*Kdrag*rhogeq*w3i*w1r**2*k*vgsol*w3r**2 + 2*w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol*w3i**2 +&
        w2i**2*Kdrag*k**2*cs**2*w3i*w1r*rhogsol*w3r**2 - w2i**2*Kdrag**2*w1r*w3r*k*vdsol*w3i**2 -&
        w2i**2*Kdrag**2*w1r*w3r**3*k*vdsol + w2i**3*Kdrag*rhogeq*w3i**4*k*vgsol +&
        2*w2i**3*Kdrag*rhogeq*w3i**2*k*vgsol*w3r**2 - rhogeq*w3r**4*w1i**2*w2r*k**2*cs**2*rhogsol -&
        rhogeq*w2r**4*k*Kdrag*vgsol*w1i*w3i**2 - rhogeq*w2r**4*k*Kdrag*vgsol*w1i*w3r**2 -&
        2*rhogeq*w3i**2*w1r*w2r**4*rhogsol*w3r**2 - 2*rhogeq**2*w1i**2*w2r**3*w3r*k*vgsol*w3i**2 -&
        2*rhogeq**2*w1i**2*w2r**3*w3r**3*k*vgsol + rhogeq**2*w2r**4*w3i**4*k*vgsol +&
        2*rhogeq**2*w3r*w2r**3*k*w1r**2*vgsol*w3i**2 + rhogeq**2*w3r**4*k*vgsol*w2r**4 +&
        2*rhogeq**2*w3r**3*w2r**3*k*w1r**2*vgsol - rhogeq*w3i**4*w1r*w2r**4*rhogsol - rhogeq*w1r*w3r**4*w2r**4*rhogsol&
        - 2*rhogeq*w3i**3*w2r**2*w1r**2*k*Kdrag*vgsol - 2*rhogeq*w3i*w2r**2*w1r**2*k*Kdrag*vgsol*w3r**2 +&
        2*vgsol*cs**2*k**3*rhogeq**2*w2r**3*w1r*w3r**2 + vdsol*k*Kdrag**2*w3r*w2r**3*w3i**2 +&
        rhogsol*cs**2*k**2*Kdrag*w3i*w2r**2*w1r*w3r**2 + cs**4*k**4*rhogeq*w3r*rhogsol*w2r**2*w3i**2 +&
        cs**4*k**4*rhogeq*w3r**3*rhogsol*w2r**2 + w1i**2*rhogeq*rhogsol*w2r**4*w3r*w3i**2 +&
        rhogsol*cs**2*k**2*Kdrag*w3i**3*w2r**2*w1r + vdsol*k*Kdrag**2*w3r**3*w2r**3 +&
        2*rhogeq**2*w3r**2*k*vgsol*w2r**4*w3i**2 - vdsol*k*Kdrag**2*w3r*w2r**2*w1r*w3i**2 -&
        vdsol*k*Kdrag**2*w3r**3*w2r**2*w1r + w1i**2*rhogeq*rhogsol*w2r**3*w3i**4 +&
        2*w1i**2*rhogeq*rhogsol*w2r**3*w3i**2*w3r**2 + w1i*rhogsol*Kdrag*w2r**3*w3i**4 +&
        w1i**2*rhogeq*rhogsol*w2r**4*w3r**3 + w1i*rhogsol*cs**2*k**2*Kdrag*w2r**3*w3r**2 -&
        vgsol*k*rhogeq**2*w3i*w2r**4*w1i*w3r**2 + w1i*rhogsol*Kdrag*w2r**4*w3r*w3i**2 + w1i*rhogsol*Kdrag*w2r**4*w3r**3&
        + w1i*rhogsol*cs**2*k**2*Kdrag*w2r**3*w3i**2 - vgsol*k*rhogeq**2*w2r**3*w1r*w3r**4 -&
        vgsol*k*rhogeq**2*w3i**3*w2r**4*w1i - 2*vgsol*k*rhogeq**2*w2r**3*w1r*w3r**2*w3i**2 +&
        cs**2*k**3*rhogeq*w3i**3*Kdrag*vdsol*w2r**2 + cs**2*k**3*rhogeq*w3i*Kdrag*vdsol*w2r**2*w3r**2 +&
        vgsol*Kdrag*k*rhogeq*w3i**3*w1r*w2r**3 + vgsol*Kdrag*k*rhogeq*w3i*w1r*w2r**3*w3r**2 +&
        w1i*rhogeq*w2r*w2i**2*w3r*k*Kdrag*vdsol*w3i**2 - vgsol*k*rhogeq**2*w2r**4*w3r*w1r*w3i**2 -&
        vgsol*k*rhogeq**2*w2r**4*w3r**3*w1r + 4*w1i*cs**4*k**4*rhogeq*w2i*w2r*rhogsol*w3i**2 -&
        w1i*w2i**3*Kdrag**2*k*vgsol*w3i**2 - w1i*w2i**3*Kdrag**2*k*vgsol*w3r**2 +&
        w1i*rhogeq*w2r*w2i**2*w3r**3*k*Kdrag*vdsol - 2*w1i*rhogeq**2*w2i*w2r**2*w3r**2*k*vgsol*w3i**2 -&
        w1i*rhogeq**2*w2i*w2r**2*w3r**4*k*vgsol + 3*w1i*cs**2*k**3*rhogeq*w2i**2*Kdrag*vdsol*w3i**2 -&
        w1i*cs**2*k**3*rhogeq*w2i**2*Kdrag*vdsol*w3r**2 + vgsol*k*Kdrag**2*w3r*w2r**2*w1r*w3i**2 +&
        2*vgsol*k*rhogeq**2*w2i**2*w2r**2*w3r**4 + vgsol*cs**2*k**3*rhogeq**2*w2i*w3i**3*w2r**2 +&
        vgsol*cs**2*k**3*rhogeq**2*w2i*w3i*w2r**2*w3r**2 - vdsol*Kdrag*k*rhogeq*w2r**2*w1r*w3r**3*w2i +&
        w1i*w2i**2*Kdrag*w3r*k**2*cs**2*rhogsol*w3i**2 + w1i*w2i**2*Kdrag*w3r**3*k**2*cs**2*rhogsol -&
        2*w1i*rhogeq**2*w2i**2*w3i**3*w2r**2*k*vgsol - w1i*w2i**2*Kdrag**2*w3i*k*vgsol*w3r**2 -&
        3*w1i*cs**2*k**3*rhogeq*w2i**2*Kdrag*vgsol*w3i**2 + w1i*cs**2*k**3*rhogeq*w2i**2*Kdrag*vgsol*w3r**2 -&
        2*w1i*rhogeq**2*w2i**2*w3i*w2r**2*k*vgsol*w3r**2 - 2*w1i**2*w2r*vgsol*k*rhogeq**2*w2i**2*w3r*w3i**2 -&
        2*w1i**2*w2r*vgsol*k*rhogeq**2*w2i**2*w3r**3 - vgsol*k*Kdrag**2*w3r**3*w2r**3 -&
        w1i**2*rhogsol*Kdrag*w2i*w2r**2*w3r*w3i**2 - w1i**2*rhogsol*Kdrag*w2i*w2r**2*w3r**3 -&
        w1i*w2i**2*Kdrag**2*w3i**3*k*vgsol + 2*w1i**2*vgsol*k*rhogeq**2*w2i*w3i*w2r**2*w3r**2 -&
        w1i**2*w2i**3*Kdrag*w3r*rhogsol*w3i**2 + w1i**2*w2r*rhogeq*rhogsol*w2i**2*w3i**4 +&
        2*w1i**2*vgsol*k*rhogeq**2*w2i*w3i**3*w2r**2 + 2*w1i**2*rhogeq**2*w2i**3*w3i*k*vgsol*w3r**2 -&
        4*w1i**2*w2r*rhogsol*cs**2*k**2*rhogeq*w3i**3*w2i - 4*w1i**2*w2r*rhogsol*cs**2*k**2*rhogeq*w3i*w2i*w3r**2 +&
        w1i**2*rhogeq*w2i**4*w3r*rhogsol*w3i**2 + 2*w1i**2*rhogeq**2*w2i**3*w3i**3*k*vgsol -&
        w1i**2*w2r*rhogsol*Kdrag*w3i*w2i**2*w3r**2 + 2*w1i**2*rhogeq*w2i**2*w3r*rhogsol*w2r**2*w3i**2 +&
        2*w1i**2*rhogeq*w2i**2*w3r**3*rhogsol*w2r**2 + 2*w1i**2*cs**2*k**3*rhogeq**2*w2i**2*vgsol*w3r**2 +&
        2*w1i**2*w2r*rhogeq*rhogsol*w2i**2*w3r**2*w3i**2 + w1i**2*w2r*rhogeq*rhogsol*w2i**2*w3r**4 -&
        w1i**2*w2r*rhogsol*Kdrag*w3i**3*w2i**2 + rhogeq*w2r*cs**4*k**4*w1i*rhogsol*w3i**3 +&
        rhogeq*w2r*cs**4*k**4*w1i*rhogsol*w3i*w3r**2 + rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol*w3i**3 +&
        rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol*w3i*w3r**2 - Kdrag*w1i*w2r**2*rhogeq*k*vgsol*w3r**4 -&
        rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i**4 + cs**4*k**4*rhogeq*w2i*w1r*rhogsol*w3i*w3r**2 -&
        rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol*w3i**3 + cs**4*k**4*rhogeq*w2i*w1r*rhogsol*w3i**3 -&
        rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol*w3i*w3r**2 - Kdrag*w1i*w2r**2*w3i**4*rhogeq*k*vgsol -&
        2*Kdrag*w1i*w2r**2*w3i**2*rhogeq*k*vgsol*w3r**2 - rhogeq*w2i**2*w1i*k*w3i**4*Kdrag*vgsol -&
        2*rhogeq*w2i**2*w1i*k*w3i**2*Kdrag*vgsol*w3r**2 + cs**2*k**3*rhogeq**2*w2i**2*w3i**4*vgsol -&
        rhogeq*w2r*k**2*cs**2*rhogsol*w3i**4*w1r**2 - rhogeq*w2i**2*w1i*k*Kdrag*vgsol*w3r**4 -&
        cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol*w3i**3 - cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol*w3i*w3r**2 -&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol*w3i*w3r**2 - 2*cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol*w3i**3 +&
        w2i**3*Kdrag*rhogeq*w3r**4*k*vgsol + w1i**2*rhogeq*w2i**4*w3r**3*rhogsol +&
        w2i**2*Kdrag**2*w3r*w1r*k*vgsol*w3i**2 - cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol*w3i**3 -&
        2*cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol*w3i*w3r**2 + cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol*w3i**3 +&
        cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol*w3i*w3r**2 + cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol*w3i**3 +&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol*w3i*w3r**2 - w3r**3*rhogsol*cs**4*k**4*rhogeq*w2r*w1r +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w3i**2 + w3r**3*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**2*w3i**2 - w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w3i**2 -&
        w3r**3*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r + 2*w3r**3*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**2 -&
        w3r*rhogsol*cs**4*k**4*rhogeq*w2r*w1r*w3i**2 - w1i**2*w2i**3*Kdrag*w3r**3*rhogsol +&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r*w3i**2 + w3r**3*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r -&
        2*w3r**3*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2 + cs**2*k**3*rhogeq**2*w3i**4*vgsol*w2r*w1r -&
        rhogeq*w2i**2*k**2*cs**2*w3i**4*w1r*rhogsol + w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r*w3i**2 -&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2*w3i**2 + w2i**2*Kdrag**2*w3r**3*w1r*k*vgsol -&
        rhogeq**2*w2i**4*w3r**3*k*w1r*vgsol + w3r*w1i*rhogsol*cs**4*k**4*rhogeq*w2i*w3i**2 +&
        w3r**3*w1i*rhogsol*cs**4*k**4*rhogeq*w2i + w3r**3*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r -&
        rhogeq**2*w3r**4*cs**2*k**3*w2r**2*vgsol - 2*rhogeq**2*w3r**2*cs**2*k**3*w1i*w2i*vgsol*w3i**2 -&
        2*rhogeq**2*w3r**2*cs**2*k**3*w2r**2*vgsol*w3i**2 - rhogeq**2*w3r**4*cs**2*k**3*w1i*w2i*vgsol +&
        2*rhogeq**2*w3r**2*cs**2*k**3*w2i**2*vgsol*w3i**2 + rhogeq**2*w3r**4*cs**2*k**3*w2i**2*vgsol -&
        rhogeq**2*w2i**4*w3r*k*w1r*vgsol*w3i**2 - 2*rhogeq*w3r**2*w2r*w1r**2*k**2*cs**2*rhogsol*w3i**2 -&
        rhogeq*w3r**4*w2r*w1r**2*k**2*cs**2*rhogsol + 2*rhogeq**2*w3r**2*cs**2*k**3*vgsol*w1r*w2r*w3i**2 +&
        rhogeq**2*w3r**4*cs**2*k**3*vgsol*w1r*w2r - 2*w2i**3*Kdrag*w1r*w3r**2*rhogsol*w3i**2 -&
        w2i**3*Kdrag*w1r*w3r**4*rhogsol - w2i**3*Kdrag*w3i**4*w1r*rhogsol + rhogeq*w2i**3*w3r*Kdrag*w1r*k*vgsol*w3i**2&
        + rhogeq*w2i**3*w3r**3*Kdrag*w1r*k*vgsol + w2i**3*Kdrag*w3r*w1r**2*rhogsol*w3i**2 +&
        w2i**3*Kdrag*w3r**3*w1r**2*rhogsol + rhogeq*w2i**4*w3r**3*w1r**2*rhogsol +&
        2*rhogeq**2*w2i**4*w3r**2*k*vgsol*w3i**2 + 2*vgsol*Kdrag*k*rhogeq*w2i**2*w3i**3*w2r**2 +&
        2*vgsol*Kdrag*k*rhogeq*w2i**2*w3i*w2r**2*w3r**2 - 2*vgsol*Kdrag*k*rhogeq*w2i*w2r**2*w1r**2*w3i**2 -&
        2*vgsol*Kdrag*k*rhogeq*w2i*w2r**2*w1r**2*w3r**2 - w1r**2*rhogsol*rhogeq*w1i*w2r**3*w3i*w3r**2 +&
        vgsol*Kdrag*k*rhogeq*w2r**2*w1r*w3r*w2i*w3i**2 + vgsol*Kdrag*k*rhogeq*w2r**2*w1r*w3r**3*w2i +&
        vgsol*Kdrag*k*rhogeq*w3i**4*w2r**2*w2i - rhogsol*Kdrag*w2r**2*w3r**4*w1r*w2i +&
        cs**4*k**4*rhogeq*w2i**2*w1r*rhogsol*w3i**2 + cs**4*k**4*rhogeq*w2i**2*w1r*rhogsol*w3r**2 +&
        2*vgsol*Kdrag*k*rhogeq*w3i**2*w2r**2*w2i*w3r**2 - cs**2*k**3*rhogeq**2*w2i**2*w2r*w3r*vgsol*w3i**2 -&
        cs**2*k**3*rhogeq**2*w2i**2*w2r*w3r**3*vgsol + rhogsol*Kdrag*w2r**2*w3r**3*w2i*w1r**2 +&
        rhogsol*Kdrag*w2r**2*w3r*w2i*w1r**2*w3i**2 - w2i**3*Kdrag*w3r*k**2*cs**2*rhogsol*w3i**2 -&
        w2i**3*Kdrag*w3r**3*k**2*cs**2*rhogsol - vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w2i*w3r**2 -&
        rhogsol*Kdrag*w2i*w3i**4*w2r**2*w1r - 2*rhogsol*Kdrag*w2i*w3i**2*w2r**2*w1r*w3r**2 -&
        rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2*w3r*w3i**2 - rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2*w3r**3 -&
        2*rhogeq*rhogsol*w2i**2*w3r**4*w2r**2*w1r - 2*rhogeq*rhogsol*w2i**2*w3i**4*w2r**2*w1r -&
        vdsol*k*Kdrag**2*w2i*w3i**3*w2r**2 - vdsol*k*Kdrag**2*w2i*w3i*w2r**2*w3r**2 +&
        rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2*w1r*w3i**2 + rhogsol*cs**2*k**2*Kdrag*w2i*w2r**2*w1r*w3r**2 +&
        vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w2i*w3i**2 - w2i*w1r**2*w1i*rhogeq**2*k*vgsol*w3i**2*w2r**2 -&
        rhogeq*w2i**3*w1r*w3r*k*Kdrag*vdsol*w3i**2 - rhogeq*w2i**3*w1r*w3r**3*k*Kdrag*vdsol -&
        2*rhogeq**2*w2i**3*w3i**3*w1r**2*k*vgsol - 2*rhogeq**2*w2i**3*w3i*w1r**2*k*vgsol*w3r**2 -&
        4*rhogeq*rhogsol*w2i**2*w3i**2*w2r**2*w1r*w3r**2 - w2r*vdsol*Kdrag*k*rhogeq*w2i**2*w3i**3*w1r -&
        2*w2r*vgsol*k*rhogeq**2*w2i**2*w1r*w3r**2*w3i**2 - 2*vgsol*k*rhogeq**2*w2i**2*w3r**3*w2r**2*w1r -&
        w2r*vdsol*Kdrag*k*rhogeq*w2i**2*w3i*w1r*w3r**2 - w2r*vgsol*k*rhogeq**2*w2i**2*w1r*w3r**4 +&
        2*w2r*vgsol*cs**2*k**3*rhogeq**2*w2i**2*w1r*w3r**2 - vgsol*cs**2*k**3*rhogeq**2*w3r**3*w2r**3 -&
        2*vgsol*k*rhogeq**2*w2i**2*w3r*w2r**2*w1r*w3i**2 + vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1i*w3r**2 +&
        w1i**2*rhogeq*rhogsol*w3r**4*w2r**3 - vgsol*cs**2*k**3*rhogeq**2*w3r*w2r**3*w3i**2 +&
        rhogsol*cs**2*k**2*Kdrag*w2r**2*w3r**3*w1i - rhogsol*cs**2*k**2*Kdrag*w2r**3*w3i**3 -&
        rhogsol*cs**2*k**2*Kdrag*w2r**3*w3i*w3r**2 + vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w1i*w3i**2 -&
        vgsol*k*rhogeq**2*w2r**3*w3i**4*w1r - vgsol*k*Kdrag**2*w3i*w2r**2*w1i*w3r**2 +&
        rhogsol*cs**2*k**2*Kdrag*w2r**2*w3r*w1i*w3i**2 + vgsol*k*Kdrag**2*w2r**3*w1r*w3i**2 +&
        vgsol*k*Kdrag**2*w2r**3*w1r*w3r**2 + vdsol*Kdrag*k*rhogeq*w3r*w2r**3*w1i*w3i**2 +&
        vdsol*Kdrag*k*rhogeq*w3r**3*w2r**3*w1i - vgsol*k*Kdrag**2*w3i**3*w2r**2*w1i -&
        3*cs**4*k**4*rhogeq*w2r**2*rhogsol*w1r*w3r**2 + vdsol*k*Kdrag**2*w3i*w2r**2*w1i*w3r**2 +&
        2*rhogeq*rhogsol*w3r**2*w2r**3*w1r**2*w3i**2 - vdsol*k*Kdrag**2*w2r**3*w1r*w3i**2 -&
        vdsol*k*Kdrag**2*w2r**3*w1r*w3r**2 + rhogsol*cs**2*k**2*rhogeq*w2r**4*w1r*w3r**2 -&
        2*cs**2*k**3*rhogeq**2*w1r**2*w2r**2*vgsol*w3i**2 + cs**4*k**4*rhogeq*w2r**2*rhogsol*w1r*w3i**2 +&
        vdsol*k*Kdrag**2*w3i**3*w2r**2*w1i - cs**2*k**3*rhogeq**2*w1i*w2i*w3i**4*vgsol -&
        rhogsol*cs**2*k**2*rhogeq*w2r**4*w1r*w3i**2 - w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r*w3i**2 -&
        w3r**3*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r - w3r*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r*w3i**2 +&
        rhogeq*rhogsol*w3r**4*w2r**3*w1r**2 - cs**2*k**3*rhogeq**2*w3i**4*w2r**2*vgsol +&
        cs**2*k**2*rhogeq*w3i**4*w1r*rhogsol*w2r**2 + w2r*rhogsol*Kdrag*w3i*w2i**2*w1r**2*w3r**2 +&
        2*w2r*rhogeq*rhogsol*w2i**2*w3r**2*w1r**2*w3i**2 + w2r*rhogeq*rhogsol*w2i**2*w3r**4*w1r**2 +&
        w2r*rhogeq*rhogsol*w2i**2*w3i**4*w1r**2 + rhogeq*rhogsol*w3i**4*w1r**2*w2r**3 +&
        w2r*rhogsol*Kdrag*w3i**3*w2i**2*w1r**2 - w3r**3*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r -&
        vgsol*k*Kdrag**2*w3r*w2r**3*w3i**2 - vdsol*Kdrag*k*rhogeq*w2r**2*w1r*w3r*w2i*w3i**2 +&
        vgsol*k*Kdrag**2*w2i*w3i**3*w2r**2 + vgsol*k*Kdrag**2*w2i*w3i*w2r**2*w3r**2 +&
        vgsol*k*Kdrag**2*w3r**3*w2r**2*w1r + 4*vgsol*k*rhogeq**2*w2i**2*w2r**2*w3r**2*w3i**2 -&
        4*w2r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w3r**2 + w2r*vgsol*k*Kdrag**2*w2i**2*w1r*w3i**2 +&
        w2r*vgsol*k*Kdrag**2*w2i**2*w1r*w3r**2 + 2*w2r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r*w3i**2 +&
        2*w2r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r**3 - w2r*vgsol*k*rhogeq**2*w2i**2*w1r*w3i**4 -&
        w2r*vgsol*k*Kdrag**2*w2i**2*w3r**3 - rhogeq*w2i**4*k**2*cs**2*w1r*rhogsol*w3i**2 +&
        rhogeq*w2i**4*k**2*cs**2*w1r*rhogsol*w3r**2 + 2*vgsol*k*rhogeq**2*w2i**2*w3i**4*w2r**2 +&
        2*rhogeq*w2i**2*k**2*cs**2*w1r*rhogsol*w2r**2*w3r**2 - w1i*vgsol*Kdrag*k*rhogeq*w2i*w3i**3*w2r**2 -&
        w1i*vgsol*Kdrag*k*rhogeq*w2i*w3i*w2r**2*w3r**2 + 2*w1i*w2r*rhogsol*cs**2*k**2*rhogeq*w2i*w3r**4 +&
        w1i*w2i**3*Kdrag**2*k*vdsol*w3i**2 + w1i*w2i**3*Kdrag**2*k*vdsol*w3r**2 -&
        2*rhogeq*w2i**2*k**2*cs**2*w1r*rhogsol*w2r**2*w3i**2 - 2*w1i*cs**2*k**3*rhogeq**2*w2i**3*vgsol*w3i**2 +&
        w1i*w2i**4*Kdrag*w3r*rhogsol*w3i**2 + 2*w1i*w2r*rhogsol*Kdrag*w2i**2*w3r**2*w3i**2 +&
        w1i*w2i*Kdrag**2*w2r**2*k*vdsol*w3r**2 + w1i*w2r*rhogsol*Kdrag*w2i**2*w3r**4 +&
        w1i*w2r*rhogsol*Kdrag*w3i**4*w2i**2 - 3*w1i*w2r*vgsol*Kdrag*k*rhogeq*w2i**2*w3r*w3i**2 -&
        3*w1i*w2r*vgsol*Kdrag*k*rhogeq*w2i**2*w3r**3 + w1i*w2i**2*Kdrag**2*w3i*k*vdsol*w3r**2 +&
        w1i*w2i**2*Kdrag*w2r*k**2*cs**2*rhogsol*w3i**2 + w1i*w2i**2*Kdrag*w2r*k**2*cs**2*rhogsol*w3r**2 +&
        w1i*w2i*Kdrag**2*w2r**2*k*vdsol*w3i**2 - w1i*rhogeq*w2i*k*w3i*Kdrag*vdsol*w2r**2*w3r**2 -&
        w1i*w2i*Kdrag**2*w2r**2*k*vgsol*w3i**2 - w1i*rhogeq*w2i*k*w3i**3*Kdrag*vdsol*w2r**2 -&
        w1i*w2i*Kdrag**2*w2r**2*k*vgsol*w3r**2 - w1i*rhogeq*w2i**3*w3i*k*Kdrag*vdsol*w3r**2 -&
        w1i*rhogeq**2*w2i**4*w3i**3*k*vgsol - w1i*rhogeq**2*w2i**4*w3i*k*vgsol*w3r**2 -&
        w1i*rhogeq*w2i**3*w3i**3*k*Kdrag*vdsol + w1i*w2i**2*Kdrag**2*w3i**3*k*vdsol -&
        2*w1i*rhogeq*w2i**2*w2r**2*k*Kdrag*vgsol*w3i**2 - 2*w1i*rhogeq*w2i**2*w2r**2*k*Kdrag*vgsol*w3r**2 -&
        w1i*rhogeq**2*w2i**3*w3r**4*k*vgsol - w1i*rhogeq**2*w2i*w2r**2*w3i**4*k*vgsol + w1i*w2i**4*Kdrag*w3r**3*rhogsol&
        + 2*w1i*w2r*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**4 + 4*w1i*w2r*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2*w3r**2 -&
        w1i*rhogeq*w2i**3*k*w3i*Kdrag*vgsol*w3r**2 - w1i*rhogeq**2*w2i**3*w3i**4*k*vgsol -&
        2*w1i*rhogeq**2*w2i**3*w3i**2*k*vgsol*w3r**2 - 2*w1i*vgsol*cs**2*k**3*rhogeq**2*w2i*w2r**2*w3i**2 -&
        w1i*rhogeq*w2i**3*k*w3i**3*Kdrag*vgsol - w1r**2*vdsol*Kdrag*k*rhogeq*w1i*w3r**2*w2r**2 -&
        4*w2i*cs**2*k**3*rhogeq**2*w3i*w1r**2*vgsol*w2r*w3r + 4*w2i*cs**2*k**3*rhogeq**2*w1i**2*w3i*vgsol*w2r*w3r -&
        w1r*rhogeq**2*k*vgsol*w1i**2*w2r**3*w3r**2 + vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w2i*w3r**2 -&
        vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w2i*w3i**2 - 2*vgsol*k*rhogeq**2*w1r**2*w2i*w3i*w2r**2*w3r**2 +&
        vgsol*Kdrag*k*rhogeq*w2r**2*w3r**4*w2i - 2*vgsol*k*rhogeq**2*w1r**2*w2i*w3i**3*w2r**2 -&
        w2r*rhogsol*cs**2*k**2*Kdrag*w3i*w2i**2*w3r**2 + w2r*vdsol*k*Kdrag**2*w2i**2*w3r*w3i**2 +&
        w2r*vdsol*k*Kdrag**2*w2i**2*w3r**3 + 4*w2r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w3r**2 -&
        w2r*rhogsol*cs**2*k**2*Kdrag*w3i**3*w2i**2 - w2r*vgsol*k*Kdrag**2*w2i**2*w3r*w3i**2 -&
        w2r*vdsol*k*Kdrag**2*w2i**2*w1r*w3i**2 - w2r*vdsol*k*Kdrag**2*w2i**2*w1r*w3r**2 +&
        w2r*vgsol*Kdrag*k*rhogeq*w2i**2*w3i**3*w1r + w2r*vgsol*Kdrag*k*rhogeq*w2i**2*w3i*w1r*w3r**2 -&
        2*w2r*rhogsol*cs**4*k**4*rhogeq*w3i**3*w2i - 2*w2r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r*w3i**2 -&
        2*w2r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r**3 + 2*w2r*vgsol*k*rhogeq**2*w1r**2*w2i**2*w3r*w3i**2 +&
        2*w2r*vgsol*k*rhogeq**2*w1r**2*w2i**2*w3r**3 - 2*w2r*rhogsol*cs**4*k**4*rhogeq*w3i*w2i*w3r**2 +&
        w2i*vdsol*Kdrag*k*rhogeq*w1i**2*w3i**2*w2r**2 + 2*w2i*w3i*w3r*rhogeq**2*k**3*cs**2*w1r**3*vgsol +&
        2*w2i*w3i*w3r*rhogeq**2*k**3*cs**2*w1r*vgsol*w1i**2 - w2i*w1i*w1r**2*rhogeq**2*k**3*cs**2*vgsol*w3r**2 +&
        w2i*w1i*w1r**2*rhogeq**2*k**3*cs**2*vgsol*w3i**2 + w2i*rhogsol*w1r**3*rhogeq*w2r**2*w3i**3 -&
        2*cs**2*k**3*rhogeq**2*w1i*w2r**3*vgsol*w3i*w3r + w2i*rhogsol*w1r**3*rhogeq*w2r**2*w3r**2*w3i -&
        2*cs**2*k**3*rhogeq**2*w1i*vgsol*w3i*w1r**2*w2r*w3r + 2*w3r*rhogeq*w3i**2*w2i**2*Kdrag*w1r*k*vgsol*w1i +&
        rhogeq*w3r**2*w3i*w1i**2*k*Kdrag*vdsol*w2r**2 + 4*w3i*cs**2*k**3*rhogeq**2*w1i*vgsol*w1r*w3r*w2r**2 +&
        2*w3r*w3i**2*rhogeq**2*k**3*cs**2*w1r*vgsol*w2i*w1i - w3r*rhogeq*w3i**2*w1i**3*w2r**2*rhogsol*w2i -&
        2*cs**2*k**3*rhogeq**2*w1i**3*vgsol*w3i*w2r*w3r + w2i*w1r**2*rhogeq*k*Kdrag*vdsol*w3r**2*w2r**2 +&
        w2i*w1r**2*rhogeq*k*Kdrag*vdsol*w3i**2*w2r**2 + w2i*w1r*rhogsol*rhogeq*w1i**2*w2r**2*w3i**3 -&
        w2i*rhogeq**2*k*vgsol*w1i**3*w3r**2*w2r**2 - w2i*rhogeq**2*k*vgsol*w1i**3*w3i**2*w2r**2 -&
        4*w2i*rhogeq*w3r**3*k**2*cs**2*w1r*rhogsol*w1i*w2r + 2*w2i*cs**2*k**3*rhogeq**2*w3i*w3r*w1r*vgsol*w2r**2 +&
        w2i*w1r*rhogsol*rhogeq*w1i**2*w2r**2*w3r**2*w3i + 4*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol*w2i*w1i -&
        4*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol*w2i*w1i - w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r**3*vgsol -&
        w2i*w1r**2*w1i*rhogeq**2*k*vgsol*w3r**2*w2r**2 + w2i*vdsol*Kdrag*k*rhogeq*w1i**2*w3r**2*w2r**2 +&
        2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol*w1r**2 - 2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol*w2i**2 -&
        2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol*w1r**2 + 2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol*w2i**2 +&
        2*w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol*w2i**2 - 2*w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol*w2i*w1i -&
        2*rhogeq*w3r**2*w1i**2*w2r*k**2*cs**2*rhogsol*w1r**2 + 2*rhogeq*w3r**2*w1i**3*w2r*k**2*cs**2*rhogsol*w2i -&
        2*w3r*w3i**2*rhogeq**2*k**3*cs**2*w2r*vgsol*w2i*w1i - 4*w3r*rhogsol*cs**2*k**2*rhogeq*w1i**2*w3i*w2r**2*w2i +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**2*w1r**2 + rhogeq*w3r**2*w3i*w2i**2*w1i**2*k*Kdrag*vdsol +&
        4*w3r*rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**2*w2i**2 - rhogeq**2*w3r**2*w3i*w2i**2*w1i**3*k*vgsol -&
        rhogsol*rhogeq*w1i**3*w2r**3*w3i*w3r**2 + w1r**4*rhogeq**2*k*vgsol*w3r**2*w2r**2 +&
        2*rhogeq*w3r**2*w2r**2*w1r**3*k**2*cs**2*rhogsol - 2*rhogeq*w3r**2*w2r**2*w1r*k**2*cs**2*rhogsol*w2i*w1i +&
        2*w3i*rhogeq*k**2*cs**2*w1i*w3r**2*w2r*rhogsol*w1r**2 - rhogeq*w3r**3*w2r**2*w1r**2*rhogsol*w2i*w1i +&
        2*rhogeq*w3r**2*cs**4*k**4*w2r*rhogsol*w1r**2 + 2*Kdrag*w3r**3*w1i*w2r**2*rhogsol*w2i**2 -&
        rhogeq**2*w3r**2*w3i*w1i*w2i**2*k*vgsol*w1r**2 - 2*rhogeq*w3r**2*cs**4*k**4*w1r*rhogsol*w2i*w1i +&
        2*w3i*rhogeq*k**2*cs**2*w2i**2*w3r**2*w1r*rhogsol*w1i - 2*rhogeq*w3r**2*w2i**3*k**2*cs**2*w1r*rhogsol*w1i +&
        4*rhogeq**2*w3r**3*w1r*w2i**3*k*vgsol*w1i - 2*rhogeq*w3r**2*cs**2*k**3*w2i*Kdrag*vgsol*w1r**2 -&
        w2i**2*w1r**2*rhogsol*rhogeq*w1i*w2r*w3r**2*w3i - w2i**2*w1r**2*rhogsol*rhogeq*w1i*w2r*w3i**3 +&
        Kdrag*w3i*w2i**4*w3r**2*rhogeq*k*vgsol - w1r**3*rhogeq**2*k*vgsol*w2r**3*w3r**2 -&
        4*w2i**2*w3i*cs**2*k**3*rhogeq**2*w1i*vgsol*w1r*w3r - w2i**2*w1r*rhogeq**2*k*vgsol*w1i**2*w2r*w3r**2 -&
        w2i**2*w1r**2*vdsol*Kdrag*k*rhogeq*w1i*w3r**2 + w2i**2*w1r**4*rhogeq**2*k*vgsol*w3r**2 -&
        2*w2i**2*cs**2*k**3*rhogeq**2*w1i*vgsol*w3i*w2r*w3r + 2*Kdrag*w3i*w2i**3*w3r**2*w1r*rhogsol*w1i -&
        2*Kdrag*w3i*w3r**2*w1r*w2r**2*rhogsol*w2i**2 + 2*Kdrag*w3i*w3r**2*w1r*w2r**2*rhogsol*w2i*w1i -&
        2*w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**3*w2r*w2i - 2*w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2*w2i**2 -&
        w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**4 + 2*w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2*w2i*w1i +&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r*w2i**2 + 2*w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r*w1r**2 +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r*w1r**2 - 2*w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i**2*w2r*w2i +&
        2*w3r*w1i*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**4 + w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r*w2i**2 -&
        rhogeq**2*w3i**3*w2i**2*w1r**2*k*vgsol*w1i + w3r*w1i*rhogsol*cs**4*k**4*rhogeq*w2i**3 +&
        w3r*w1i*rhogsol*cs**4*k**4*rhogeq*w2i*w1r**2 - w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**4 -&
        2*w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**2 + 2*w3r*w1i**3*rhogsol*cs**2*k**2*rhogeq*w2i**3 +&
        2*rhogeq**2*w3r**3*cs**2*k**3*w1r*vgsol*w2i*w1i - 4*w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**3*w3i +&
        2*w3r*w1i**3*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i + w2i**2*w1r**2*Kdrag**2*k*vdsol*w3r**2 -&
        w2i**2*w1r**3*rhogeq**2*k*vgsol*w2r*w3r**2 - w2i**2*w1r**3*rhogeq**2*k*vgsol*w2r*w3i**2 -&
        rhogeq**2*w3r**2*cs**2*k**3*vgsol*w1r**3*w2r + w2i**2*w1r**2*Kdrag**2*k*vdsol*w3i**2 -&
        w2i**2*rhogsol*rhogeq*w1i**3*w2r*w3r**2*w3i - w2i**2*rhogsol*rhogeq*w1i**3*w2r*w3i**3 +&
        4*rhogeq**2*w3r**2*cs**2*k**3*vgsol*w1r*w2r*w2i*w1i - 2*rhogeq**2*w3r**2*cs**2*k**3*w2r**2*vgsol*w2i**2 -&
        2*rhogeq**2*w3r**3*cs**2*k**3*w2r*vgsol*w2i*w1i - 2*rhogeq*w3r**2*w2r*w1r**2*k**2*cs**2*rhogsol*w2i**2 -&
        rhogeq*w3r**2*w2r*w1r**4*k**2*cs**2*rhogsol + 2*rhogeq*w3r**2*w2r*w1r**2*k**2*cs**2*rhogsol*w2i*w1i +&
        2*rhogeq*w3r**3*w1r*w2i**2*k*Kdrag*vgsol*w1i - w3r*rhogeq*w3i**2*w1i**2*rhogsol*w2i**2*w2r*w1r +&
        2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vdsol*w1i**2 + 2*w3r*w3i*rhogeq*k**3*cs**2*w2r**3*Kdrag*vdsol -&
        2*w3r*w1i**2*rhogsol*cs**4*k**4*rhogeq*w2i**2 + 2*w3r*w1i*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2*w1r**2 -&
        w3r*rhogeq**2*w3i**2*w1i**2*w2r**2*k*vgsol*w1r + 2*w3r*rhogeq*w3i**2*w1i*w2r**2*k*Kdrag*vgsol*w1r -&
        w3r*rhogeq*w3i**2*w1i**2*w2r**3*rhogsol*w1r - w3r*rhogeq**2*w3i**2*w2r**2*w1r**3*k*vgsol -&
        2*rhogeq*w2i**2*k**2*cs**2*w3i**2*w1r**3*rhogsol + 2*rhogeq*w2i**3*k**2*cs**2*w3i**2*w1r*rhogsol*w1i +&
        2*Kdrag*w3i**3*w1r*w2r**2*rhogsol*w2i*w1i - 2*cs**4*k**4*rhogeq*w3i**2*w1r*rhogsol*w2i*w1i -&
        rhogeq*w3r**3*w2i**2*w1i**2*rhogsol*w2r*w1r + rhogeq*w3i**3*w1i**2*w2i**2*k*Kdrag*vdsol +&
        cs**2*k**3*rhogeq**2*w3i**2*vgsol*w2r*w1r**3 - 4*cs**2*k**3*rhogeq**2*w3i**2*vgsol*w2r*w1r*w2i*w1i -&
        rhogeq**2*w3r**3*w2r**2*w1i**2*k*vgsol*w1r - rhogeq*w3r**2*w1i**4*w2r*k**2*cs**2*rhogsol +&
        2*rhogeq*w3r**2*w1i**2*w2r**2*k**2*cs**2*rhogsol*w1r - w3r*rhogeq**2*k**3*cs**2*w2r**2*w1r*vgsol*w1i**2 -&
        4*w3r*w3i*rhogeq*k**3*cs**2*w2r**2*Kdrag*vdsol*w1r - 2*w3r*w3i*rhogeq*k**3*cs**2*w2r**3*Kdrag*vgsol +&
        4*w3r*w3i*rhogeq*k**3*cs**2*w2r**2*Kdrag*vgsol*w1r - w3r*rhogeq**2*k**3*cs**2*w2r**3*w1r**2*vgsol +&
        w2i**2*w1r**4*rhogeq**2*k*vgsol*w3i**2 + w2i**2*vgsol*k*rhogeq**2*w1i**4*w3r**2 +&
        w2i**2*vgsol*k*rhogeq**2*w1i**4*w3i**2 - 4*w3r*rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**3*w1r -&
        2*rhogeq**2*w3r**2*cs**2*k**3*w1i*w3i*vgsol*w2r*w1r - 2*rhogeq*w3r**2*cs**2*k**3*w1i*Kdrag*vgsol*w2r*w1r +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w1i**3*w3i*w2r**2 + 2*w3r*rhogsol*cs**2*k**2*rhogeq*w1i*w3i*w2r**4 -&
        2*w3r*w3i*rhogeq*k**3*cs**2*w2r*Kdrag*vgsol*w1i**2 + w3r*rhogeq**2*k**3*cs**2*w2r**4*w1r*vgsol -&
        w2i**2*vdsol*Kdrag*k*rhogeq*w1i**3*w3r**2 - w2i**2*w1r*rhogeq**2*k*vgsol*w1i**2*w2r*w3i**2 -&
        w2i**2*vdsol*Kdrag*k*rhogeq*w1i**3*w3i**2 + 2*w2i**2*w1r**2*rhogeq**2*k*vgsol*w1i**2*w3r**2 +&
        2*w2i**2*w1r**2*rhogeq**2*k*vgsol*w1i**2*w3i**2 - w2i**2*w1r**2*Kdrag**2*k*vgsol*w3r**2 -&
        w2i**2*w1r**2*Kdrag**2*k*vgsol*w3i**2 + w2i**2*w1r**2*rhogeq*k*Kdrag*vdsol*w3r**2*w3i +&
        w2i**2*w1r**2*rhogeq*k*Kdrag*vdsol*w3i**3 - w2i**2*w1r**2*vdsol*Kdrag*k*rhogeq*w1i*w3i**2 -&
        rhogeq**2*w3r**3*w2r**2*w1r**3*k*vgsol + 4*rhogeq**2*w3r**2*w3i*w1i*w2i**2*k*vgsol*w2r*w1r -&
        rhogeq**2*w3r**3*w1r*w2i**2*k*vgsol*w1i**2 - w2i**3*rhogeq**2*w1i*k*vgsol*w1r**2*w3r**2 -&
        w2i**3*rhogeq**2*w1i*k*vgsol*w3i**2*w1r**2 - w2i**3*rhogeq**2*w1i**3*w3i**2*k*vgsol -&
        rhogeq*w3r**3*w2r**3*w1i**2*rhogsol*w1r - w2i**3*rhogeq**2*w1i**3*k*vgsol*w3r**2 +&
        2*rhogeq*w3r**3*w2r**2*w1i*k*Kdrag*vgsol*w1r - 2*w3i*rhogeq*k**2*cs**2*w1i*w3r**2*w2r**2*rhogsol*w1r +&
        2*w3i*rhogeq*k**2*cs**2*w1i**3*w3r**2*w2r*rhogsol - 2*Kdrag*w3r**3*w1i*w2r**3*rhogsol*w1r -&
        2*rhogeq*w3r**2*w2r**3*w1r**2*k**2*cs**2*rhogsol + 2*w2i**3*cs**2*k**3*rhogeq**2*w3i*w3r*w1r*vgsol -&
        w3i**2*w1r*rhogeq**2*k*vgsol*w1i**2*w2r**3 - w3i**2*w1r**3*rhogeq**2*k*vgsol*w2r**3 -&
        2*Kdrag*w3r**3*w1i*rhogsol*w2i**2*w2r*w1r + w2i**3*rhogsol*w1r**3*rhogeq*w3r**2*w3i +&
        Kdrag**2*w1i**2*k*vgsol*w3r**2*w2r**2 + w2i**3*w1r*rhogsol*rhogeq*w1i**2*w3r**2*w3i +&
        w2i**3*w1r**2*rhogeq*k*Kdrag*vdsol*w3r**2 + w2i**3*w1r*rhogsol*rhogeq*w1i**2*w3i**3 +&
        w2i**3*w1r**2*rhogeq*k*Kdrag*vdsol*w3i**2 + w2i**3*vdsol*Kdrag*k*rhogeq*w1i**2*w3r**2 +&
        w2i**3*vdsol*Kdrag*k*rhogeq*w1i**2*w3i**2 + w2i**3*rhogsol*w1r**3*rhogeq*w3i**3 -&
        4*w3r*w1i*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2*w2r*w1r - w3r*w1i**4*rhogsol*cs**2*k**2*rhogeq*w2i**2 -&
        2*w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w2r**2 + 2*w3r*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w2r*w1r +&
        w3r*w1i**3*rhogsol*cs**4*k**4*rhogeq*w2i + w3r*w1i*rhogsol*cs**4*k**4*rhogeq*w2i*w2r**2 -&
        rhogeq**2*w3r**2*cs**2*k**3*w1i**3*w2i*vgsol + 2*rhogeq*w3r**2*cs**2*k**3*w1i*Kdrag*vdsol*w2r*w1r -&
        rhogeq**2*w3r**2*cs**2*k**3*vgsol*w1r*w2r*w1i**2 - rhogeq**2*w3r**2*cs**2*k**3*w2r**4*vgsol +&
        2*rhogeq**2*w3r**2*cs**2*k**3*w2i*w3i*vgsol*w2r*w1r - 2*rhogeq*w2i**2*k**2*cs**2*w3i**2*w1r*rhogsol*w1i**2 -&
        2*Kdrag*w1i*w3i**2*w2r**2*k**2*cs**2*rhogsol*w1r + 2*rhogeq*w2i**2*k**2*cs**2*w3i**2*w1r**2*rhogsol*w2r -&
        Kdrag**2*w3i**2*w2r**2*w1r**2*k*vgsol - Kdrag**2*w2r**2*w1r**2*k*vgsol*w3r**2 -&
        4*rhogeq*w3i**2*w1i*w2i*w3r*k**2*cs**2*rhogsol*w2r*w1r - Kdrag**2*w1i**2*k*vdsol*w3r**2*w2r**2 +&
        cs**2*k**3*rhogeq**2*w3i**2*vgsol*w2r*w1r*w1i**2 + w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2r**3 -&
        w3r*rhogsol*cs**2*k**2*rhogeq*w2r**4*w1r**2 + 2*w3r*rhogsol*cs**2*k**2*rhogeq*w2r**3*w1r**3 +&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1i**4*w2r - 2*w3r*rhogsol*cs**2*k**2*rhogeq*w2r**2*w1r**2*w1i**2 +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r**3 - 2*w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r**2*w1r +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i**3*w2r - w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r**3 -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**3*w1r + 2*w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w1r*w1i -&
        w3r*rhogeq*w3i**2*w2r**2*w1r**2*rhogsol*w2i*w1i - rhogeq*w3r**3*w2i**3*w1i**3*rhogsol -&
        w3r*rhogeq*w3i**2*w2i**3*w1r**2*rhogsol*w1i + 4*w3r*rhogeq**2*w3i**2*w2i**3*k*w1r*vgsol*w1i -&
        w3r*rhogeq**2*w3i**2*w2i**2*k*w1r**3*vgsol - 2*cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol*w3i*w1r**2 +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**3*w1r + w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r**3 -&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol*w3i*w1r**2 + 2*cs**2*k**3*rhogeq*w1i**2*w2i**2*Kdrag*vgsol*w3i -&
        cs**2*k**3*rhogeq*w1i*w2i**3*Kdrag*vgsol*w3i - cs**2*k**3*rhogeq**2*w1i*w2i**4*vgsol*w3i +&
        cs**2*k**3*rhogeq**2*w1i**2*w2i**3*vgsol*w3i - 2*w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w1r*w1i -&
        w3r*rhogsol*cs**4*k**4*rhogeq*w2r*w1r*w2i**2 - w3r*rhogsol*cs**4*k**4*rhogeq*w2r*w1r**3 -&
        w3r*rhogsol*cs**2*k**2*rhogeq*w2i**4*w1r**2 - w3r*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**4 +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w2i**3*w1r**2*w1i + 4*w3r*rhogeq**2*w3i**2*w2r**2*k*w1r*vgsol*w2i*w1i -&
        2*w3r*w3i*rhogeq*k**4*cs**4*w2i**3*rhogsol - 2*w3r*w3i*rhogeq*k**4*cs**4*w2i*rhogsol*w1r**2 +&
        w3r*rhogeq**2*k**3*cs**2*w2i**4*w1r*vgsol + w3r*rhogeq**2*k**3*cs**2*w2i**2*w1r**3*vgsol +&
        2*w3r*Kdrag*w3i**2*w1i*w2r**2*rhogsol*w2i**2 - 2*w3r*rhogeq**2*k**3*cs**2*w2i**3*w1r*vgsol*w1i +&
        4*w3r*w3i*rhogeq*k**4*cs**4*w2i**2*rhogsol*w1i - 2*cs**4*k**4*rhogeq*w2i**2*w1r*rhogsol*w3i*w1i +&
        cs**4*k**4*rhogeq*w2i**3*w1r*rhogsol*w3i + cs**4*k**4*rhogeq*w2i*w1r**3*rhogsol*w3i +&
        cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol*w3i*w1r**2 + 2*cs**2*k**2*rhogeq*w2i**2*w3i**3*w1r*rhogsol*w1i -&
        rhogeq*w2i**4*w1i*k*w3i**2*Kdrag*vgsol + rhogeq*w2i**2*w1i*k*w3i**2*Kdrag*vgsol*w1r**2 +&
        cs**2*k**3*rhogeq**2*w2i**4*w3i**2*vgsol + 2*rhogeq*w2r*k**2*cs**2*rhogsol*w3i**2*w1r**2*w2i*w1i -&
        rhogeq*w2i**4*w1i*k*Kdrag*vgsol*w3r**2 - cs**2*k**3*rhogeq**2*w2i**3*w1r**2*vgsol*w3i -&
        cs**2*k**3*rhogeq**2*w2i*w1r**4*vgsol*w3i + w2i**2*Kdrag**2*w1i**2*k*w3i**2*vgsol -&
        w2r**2*rhogeq**2*w1i*w3i**3*k*vgsol*w1r**2 - rhogeq*w2r*k**2*cs**2*rhogsol*w3i**2*w1r**4 -&
        2*cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol*w3i*w2i**2 - cs**2*k**3*rhogeq**2*w1i*w2r**2*vgsol*w3i*w1r**2 +&
        cs**2*k**3*rhogeq**2*w1i**2*w2r**2*vgsol*w3i*w2i + cs**2*k**3*rhogeq*w1i*w2i**3*Kdrag*vdsol*w3i -&
        2*cs**2*k**3*rhogeq*w1i**2*w2i**2*Kdrag*vdsol*w3i + cs**2*k**3*rhogeq**2*w1i**3*w2i**2*vgsol*w3i +&
        cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol*w3i*w1r**2 - 2*rhogeq*w2r*w1i**2*k**2*cs**2*rhogsol*w3i**2*w1r**2 +&
        Kdrag*w1i*w2r**2*w3i**2*rhogeq*k*vgsol*w1r**2 + Kdrag*w1i*w2r**2*rhogeq*k*vgsol*w3r**2*w1r**2 +&
        w2i**2*Kdrag*w1i**3*rhogeq*k*vgsol*w3r**2 - w3r*rhogeq*w3i**2*w1i**3*rhogsol*w2i**3 +&
        rhogeq*w2i**2*w1i*k*Kdrag*vgsol*w3r**2*w1r**2 + cs**2*k**3*rhogeq**2*w1i**3*w3i**2*vgsol*w2i -&
        rhogeq*w2r*cs**2*k**3*w1r**3*Kdrag*vdsol*w3i + rhogeq*w2r*cs**4*k**4*w1i*rhogsol*w3i*w1r**2 -&
        2*rhogeq*w2r*cs**4*k**4*w1i**2*rhogsol*w3i*w2i + 2*rhogeq*w2r*w1i**3*k**2*cs**2*rhogsol*w3i**2*w2i +&
        rhogeq*w2r*cs**4*k**4*w1i*rhogsol*w3i*w2i**2 + 2*rhogeq*w2r*w1i*k**2*cs**2*rhogsol*w3i**3*w1r**2 -&
        rhogeq**2*w3r**3*w1r**3*w2i**2*k*vgsol + w2i**2*Kdrag**2*w1i**2*k*vgsol*w3r**2 +&
        rhogeq*w2r*cs**2*k**3*Kdrag*w1r**3*vgsol*w3i + rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol*w3i*w2i**2 -&
        rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol*w3i*w2i**2 - 2*w2i**2*Kdrag*k**2*cs**2*w1r*rhogsol*w3r**2*w1i -&
        2*w2i**2*Kdrag*k**2*cs**2*w3i**2*w1r*rhogsol*w1i - rhogeq*w3r**3*w2r**2*w1i**3*rhogsol*w2i +&
        rhogeq*w3i**3*k*Kdrag*vdsol*w1i**2*w2r**2 - rhogeq*w3r**3*w1r**2*w2i**3*rhogsol*w1i +&
        w2i**4*Kdrag*rhogeq*w3i**3*k*vgsol - vdsol*k*Kdrag**2*w1i**2*w2r**2*w3i**2 +&
        2*cs**2*k**2*rhogeq*w3i**2*w1r*rhogsol*w2r**2*w2i*w1i + 2*cs**2*k**3*rhogeq*w1i**2*w3i**2*Kdrag*vgsol*w2i -&
        2*cs**2*k**3*rhogeq*w1i**2*w3i**2*Kdrag*vdsol*w2i + 2*cs**2*k**3*rhogeq**2*w3i**2*w2r**2*vgsol*w2i**2 +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w1i**3*w2r**2*w2i + 2*w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i**2*w2r*w2i +&
        w3r*vgsol*cs**2*k**3*rhogeq**2*w1r**4*w2r - w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r*w2i**2 -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r*w1r**2 - w3r*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r*w2i**2 -&
        2*w3r*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2r*w2i*w1i + 2*rhogeq*w3i**2*w2r*w1r**3*w3r*k**2*cs**2*rhogsol -&
        2*w3r*w3i*rhogeq*k**4*cs**4*w1i**2*rhogsol*w2i - rhogeq**2*w3r**2*w1i*w3i*w2r**2*k*vgsol*w1r**2 +&
        2*rhogeq*w3r**2*cs**2*k**3*w2i*Kdrag*vdsol*w1r**2 + w2i**2*Kdrag*w1i**3*w3i**2*rhogeq*k*vgsol +&
        4*rhogeq**2*w3r**3*w2r**2*k*vgsol*w1r*w2i*w1i - w3r*rhogeq*w3i**2*w2r**3*w1r**3*rhogsol -&
        w3r*rhogeq*w3i**2*w2i**2*w1r**3*rhogsol*w2r - w3r*rhogeq**2*w3i**2*w2i**2*k*w1r*vgsol*w1i**2 -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w2r**2 + 2*w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r**2*w2r -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w1i**2 - rhogeq**2*w3r**2*cs**2*k**3*w2i**4*vgsol -&
        2*w3r*Kdrag*w3i**2*w1i*rhogsol*w2i**2*w2r*w1r - w3r*rhogsol*cs**4*k**4*rhogeq*w2r*w1r*w1i**2 -&
        2*w3r*w3i*rhogeq*k**4*cs**4*w2i*rhogsol*w2r**2 - 2*w3r*Kdrag*w3i**2*w1i*w2r**3*rhogsol*w1r +&
        w3r*rhogeq**2*k**3*cs**2*w2i**2*w1r*vgsol*w1i**2 + 4*w3r*w3i*rhogeq*k**4*cs**4*w2i*rhogsol*w2r*w1r -&
        w3r*rhogsol*cs**4*k**4*rhogeq*w2r**3*w1r + 2*w3r*rhogsol*cs**4*k**4*rhogeq*w2r**2*w1r**2 +&
        w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w1i**2 + 2*w3r*rhogsol*cs**2*k**2*rhogeq*w2i**2*w1r**3*w2r +&
        rhogeq**2*w3i**2*w1r**4*k*vgsol*w2r**2 + 2*rhogeq**2*w1i**2*w3i**2*k*vgsol*w1r**2*w2r**2 +&
        vgsol*k*Kdrag**2*w1i**2*w2r**2*w3i**2 + w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r*w2r**2 -&
        2*w3r*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w1r**2*w2r - cs**2*k**3*rhogeq**2*w1i**3*w2r**2*vgsol*w3i +&
        2*cs**2*k**3*rhogeq**2*w1i**2*w2i*vgsol*w3i*w2r*w1r + cs**2*k**3*rhogeq*w1i**3*w2i*Kdrag*vdsol*w3i -&
        cs**2*k**3*rhogeq**2*w1i**4*w2i*vgsol*w3i + cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vdsol*w3i*w2r**2 -&
        cs**2*k**3*rhogeq**2*w1i*w2r**4*vgsol*w3i - cs**2*k**3*rhogeq*w1i*w2i*Kdrag*vgsol*w3i*w2r**2 +&
        2*cs**2*k**3*rhogeq**2*w1i*w2i**2*vgsol*w3i*w2r*w1r + 2*cs**2*k**3*rhogeq**2*w1i*w2r**3*vgsol*w3i*w1r +&
        2*w2i**3*Kdrag*w3i**3*w1r*rhogsol*w1i + 2*rhogeq*w2r**2*cs**2*k**3*w1r**2*Kdrag*vdsol*w3i +&
        Kdrag*w1i**3*w2r**2*w3i**2*rhogeq*k*vgsol + 2*cs**2*k**3*rhogeq**2*w2i*w3i**3*vgsol*w2r*w1r +&
        2*cs**2*k**3*rhogeq**2*w2i*w1r**3*vgsol*w3i*w2r + cs**4*k**4*rhogeq*w2i*w1r*rhogsol*w3i*w1i**2 +&
        cs**4*k**4*rhogeq*w2i*w1r*rhogsol*w3i*w2r**2 - cs**2*k**3*rhogeq**2*w2i*w1r**2*vgsol*w3i*w2r**2 -&
        2*cs**4*k**4*rhogeq*w2i*w1r**2*rhogsol*w3i*w2r + 4*rhogeq**2*w2i**2*w1i*w3i**3*k*vgsol*w2r*w1r -&
        rhogeq*w2r*cs**2*k**3*w1r*Kdrag*vdsol*w3i*w1i**2 - rhogeq*w2r**3*cs**2*k**3*w1r*Kdrag*vdsol*w3i -&
        cs**2*k**3*rhogeq*w1i**3*w2i*Kdrag*vgsol*w3i + 2*rhogeq*w2i**2*w1i*k*Kdrag*vgsol*w3r**2*w2r*w1r +&
        2*rhogeq*w2r**3*k**2*cs**2*rhogsol*w3i**2*w1r**2 + 2*rhogeq*w2i**2*w1i*k*w3i**2*Kdrag*vgsol*w2r*w1r -&
        w2r**2*rhogeq**2*w1i**3*w3i**3*k*vgsol + 4*w2r**3*rhogeq**2*w1i*w3i**3*k*vgsol*w1r -&
        2*rhogeq*w2r**2*w1i*k**2*cs**2*rhogsol*w3i**3*w1r + 2*rhogeq*w2r*w1i**3*k**2*cs**2*rhogsol*w3i**3 -&
        rhogeq*w2r*w1i**4*k**2*cs**2*rhogsol*w3i**2 + 2*Kdrag*w1i*w2r**3*w3i**2*rhogeq*k*vgsol*w1r +&
        Kdrag*w1i**3*w2r**2*rhogeq*k*vgsol*w3r**2 + 2*Kdrag*w1i*w2r**3*rhogeq*k*vgsol*w3r**2*w1r -&
        2*rhogeq*w2r*cs**4*k**4*w3i**2*rhogsol*w1i**2 + rhogeq*w2r**3*cs**2*k**3*Kdrag*w1r*vgsol*w3i -&
        2*rhogeq*w2r**2*cs**2*k**3*Kdrag*w1r**2*vgsol*w3i - 2*rhogeq*w2r**2*cs**4*k**4*w1i*rhogsol*w3i*w1r +&
        rhogeq*w2r*cs**2*k**3*Kdrag*w1r*vgsol*w3i*w1i**2 + rhogeq*w2r*cs**4*k**4*w1i**3*rhogsol*w3i +&
        rhogeq*w2r**3*cs**4*k**4*w1i*rhogsol*w3i - w2i**2*Kdrag**2*w1i**2*k*vdsol*w3r**2 -&
        w2i**2*Kdrag**2*w1i**2*k*w3i**2*vdsol + 2*w1r**3*rhogeq*k**2*cs**2*rhogsol*w3r**3*w2r -&
        2*Kdrag*w3i**3*w1r*w2r**2*rhogsol*w2i**2 + 2*w1r**2*rhogeq**2*k*vgsol*w1i**2*w3r**2*w2r**2 +&
        rhogeq*w3r**2*w3i*w2r**2*w1r**2*k*Kdrag*vdsol - rhogeq*w3r**3*w1r**3*w2i**2*rhogsol*w2r -&
        rhogeq**2*w3r**2*w1i**3*w3i*w2r**2*k*vgsol + 4*rhogeq**2*w3r**2*w1i*w3i*w2r**3*k*vgsol*w1r +&
        2*rhogeq*w3i**2*w2r*w1r*w3r*k**2*cs**2*rhogsol*w1i**2 - w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i**3*w2r -&
        w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r**3 + 2*w3r*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2r**2*w1r -&
        w3r*rhogsol*cs**2*k**2*rhogeq*w1i**4*w2r**2 - w3r*rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**4 +&
        2*w3r*rhogsol*cs**2*k**2*rhogeq*w1i**2*w2r**3*w1r + 2*cs**2*k**3*rhogeq*w1i*w3i**2*Kdrag*vdsol*w2r*w1r -&
        2*cs**2*k**3*rhogeq*w1i*w3i**2*Kdrag*vgsol*w2r*w1r + cs**2*k**3*rhogeq**2*w3i**2*w2r**4*vgsol +&
        Kdrag**2*w3i**2*w2r**2*w1r**2*k*vdsol + Kdrag**2*w2r**2*w1r**2*k*vdsol*w3r**2 -&
        2*cs**2*k**3*rhogeq**2*w1i*w3i**3*vgsol*w2r*w1r + rhogeq*w2r**2*w1r**2*k*Kdrag*vdsol*w3i**3 +&
        rhogeq**2*w1i**4*w3i**2*k*vgsol*w2r**2 - rhogeq*w1i*k*w3i**2*Kdrag*vdsol*w1r**2*w2r**2 -&
        rhogeq*w1i**3*k*w3i**2*Kdrag*vdsol*w2r**2 - rhogeq*w1i*w3i**3*w2r**3*rhogsol*w1r**2 -&
        2*Kdrag*w1i*w2r**2*k**2*cs**2*rhogsol*w3r**2*w1r - rhogeq**2*w3i**3*w1i**3*w2i**2*k*vgsol +&
        2*w1r*rhogsol*cs**2*k**2*w1i**2*rhogeq*w3r**3*w2r - vdsol*Kdrag*k*rhogeq*w1i**3*w3r**2*w2r**2 +&
        vgsol*k*rhogeq**2*w1i**4*w3r**2*w2r**2)/( - 2*w1i*w3i + w1i**2 + w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w3i**2&
        + w3r**2)/rhogeq/(w2i**2 + w2r**2)/Kdrag/(w2i**2 + w1r**2 - 2*w2i*w1i + w1i**2 + w2r**2 - 2*w2r*w1r) 

  rhod1i = - (rhodeq*vgsol*cs**2*k**3*rhogeq**2*w3i**3*w2r**3 - rhodsol*rhogeq*Kdrag*w2i**4*w3i**4 +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w2i**3*w3r**3 - rhodsol*rhogeq*Kdrag*w2i**4*w3r**4 -&
        rhodeq*rhogeq*rhogsol*w1r**2*w1i*w3i**3*w2i**3 - rhodsol*rhogeq*Kdrag*w1r**2*w3r**4*w2i**2 -&
        4*rhodsol*rhogeq*Kdrag*w1r**2*w2r**3*w3r**3 - rhodsol*rhogeq*Kdrag*w1r**2*w3i**2*w2r**4 -&
        rhodsol*rhogeq*Kdrag*w1r**2*w3r**2*w2r**4 - rhodsol*rhogeq*Kdrag*w1r**2*w3i**4*w2r**2 -&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w3i**2*w2r**2*w3r**2 - 4*rhodsol*rhogeq*Kdrag*w1r**2*w2r**3*w3r*w3i**2 -&
        rhodsol*rhogeq*Kdrag*w1r**2*w3r**4*w2r**2 + 2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2i**2*w3r**3 -&
        rhodsol*rhogeq*Kdrag*w1r**2*w3r**2*w2i**4 - 4*rhodsol*rhogeq*Kdrag*w1r**2*w2r*w3r*w2i**2*w3i**2 -&
        4*rhodsol*rhogeq*Kdrag*w1r**2*w2r*w3r**3*w2i**2 - 2*rhodsol*rhogeq*Kdrag*w1r**2*w3i**2*w2i**2*w3r**2 -&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w3i**2*w2i**2*w2r**2 - rhodsol*rhogeq*Kdrag*w1r**2*w3i**4*w2i**2 -&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w3r**2*w2i**2*w2r**2 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**3*w3i*w2i**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**3*w2i*w2r*w3r - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**3*w2i*w3i**2 +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**3*w2i*w3r**2 + 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**3*w3r*w3i*w2r +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**3*w3i*w2r**2 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w2i**2*w3r*w3i -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w2i**2*w3i*w2r -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w2i*w2r*w3r**2 -&
        rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2i**2*w2r*w3i**2 - rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w3i**2*w2r**3 -&
        rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2r**3*w3r**2 - rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2r**2*w3r**3 -&
        rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2r**2*w3r*w3i**2 + rhodeq*vgsol*k*rhogeq**2*w1r**3*w2i**3*w3i**2 +&
        rhodeq*vgsol*k*rhogeq**2*w1r**3*w2i**3*w3r**2 + rhodeq*vgsol*k*rhogeq**2*w1r**3*w2i**2*w3i**3 -&
        rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2i**2*w3r*w3i**2 + rhodeq*vgsol*k*rhogeq**2*w1r**3*w2i**2*w3i*w3r**2 +&
        4*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w3i**2*w2i**2 -&
        8*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w2i*w3r*w3i*w2r +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w3i*w2i*w3r**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w3i**3*w2i +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w2i*w3i*w2r**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w3i**2*w3r*w2r +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w3r**3*w2r + 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w2r**3*w3r&
        + 4*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w2r**2*w3r**2 - rhodeq*vdsol*k*Kdrag**2*w2i**3*w3r**3 +&
        2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2i**2*w2r*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2i**2*w3r*w3i**2 +&
        2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2i**2*w2r*w3i**2 + 2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w3i**2*w2r**3 +&
        2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2r**3*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2r**2*w3r**3 +&
        2*rhodsol*rhogeq*Kdrag*w1i**2*w1r*w2r**2*w3r*w3i**2 - rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2i**2*w2r*w3r**2 -&
        rhodeq*vgsol*k*rhogeq**2*w1i*w1r**2*w2i**2*w3r**3 - rhodsol*rhogeq*Kdrag*w1r**4*w3i**2*w2i**2 -&
        rhodeq*w1i*vgsol*k*Kdrag**2*w2i**2*w3r**3 + rhodeq*rhogsol*cs**2*k**2*Kdrag*w3i**3*w2i**3 +&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w3i*w2i**3*w3r**2 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w2r*w3r**3*w2i**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w2r*w3r*w2i**2*w3i**2 + rhodeq*rhogsol*cs**2*k**2*Kdrag*w3i**3*w2i*w2r**2 +&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w3i*w2i*w2r**2*w3r**2 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w2r**3*w3r**3 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w2r**3*w3r*w3i**2 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w3r*w3i**2 +&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w3i**2*w2r**3 + rhodeq*vgsol*Kdrag*k*rhogeq*w2i**4*w3r**3 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w2i**4*w3r*w3i**2 + 2*rhodeq*vgsol*Kdrag*k*rhogeq*w2r**2*w2i**2*w3i**2*w3r +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w2i**2*w3r**4*w2r + rhodeq*vgsol*Kdrag*k*rhogeq*w2i**2*w3i**4*w2r +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w2r**2*w2i**2*w3r**3 + 2*rhodeq*vgsol*Kdrag*k*rhogeq*w2i**2*w3r**2*w3i**2*w2r +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w2r**4*w3r**3 + rhodeq*vgsol*Kdrag*k*rhogeq*w2r**3*w3r**4 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w2r**4*w3i**2*w3r + rhodeq*vgsol*Kdrag*k*rhogeq*w2r**3*w3i**4 +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w2r**3*w3i**2*w3r**2 + rhodeq*vgsol*k*Kdrag**2*w2i**3*w3r**3 +&
        rhodeq*vgsol*k*Kdrag**2*w2i**3*w3r*w3i**2 + rhodeq*vgsol*k*Kdrag**2*w3i*w2r*w3r**2*w2i**2 +&
        rhodeq*vgsol*k*Kdrag**2*w2r*w3i**3*w2i**2 + rhodeq*vgsol*k*Kdrag**2*w2i*w2r**2*w3i**2*w3r +&
        rhodeq*vgsol*k*Kdrag**2*w2i*w2r**2*w3r**3 + rhodeq*vgsol*k*Kdrag**2*w3i**3*w2r**3 +&
        rhodeq*vgsol*k*Kdrag**2*w3i*w2r**3*w3r**2 + 2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**3*w3r*w3i +&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w3r*w3i**2 + rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w2r*w3i**2&
        + rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w3r**3 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w2r*w3r**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r*w3i*w2r**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i*w2r*w3r**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i**3*w2r - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**3*w3r**2 -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i**2*w3r*w3i*w2r - rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i*w2r**2*w3r**2 +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i*w3i**2*w2r**2 - 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i*w3r**3*w2r -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i*w3i**2*w3r*w2r - rhodeq*rhogsol*cs**4*k**4*rhogeq*w3i**3*w2r**2 -&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w3i*w2r**2*w3r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i**2*w3i*w3r**2 -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w3r*w3i*w2r**3 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w3r**3 -&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2r**4 - rhodeq*w1i**4*vgsol*cs**2*k**3*rhogeq**2*w3r*w2i -&
        rhodeq*w1i*rhogsol*Kdrag*w2i**3*w3r**4 + rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**4*w2i*w3i**2 +&
        rhodeq*w1i**2*rhogsol*Kdrag*w2r*w3r*w2i**2*w3i**2 + 2*rhodsol*rhogeq*Kdrag*w1i**3*w2i**3*w3i**2 +&
        2*rhodsol*rhogeq*Kdrag*w1i**3*w2i**3*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1i**3*w2i**2*w3i**3 +&
        2*rhodsol*rhogeq*Kdrag*w1i**3*w2i**2*w3i*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1i**3*w2i*w3i**2*w2r**2 +&
        2*rhodsol*rhogeq*Kdrag*w1i**3*w2i*w2r**2*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1i**3*w3i**3*w2r**2 +&
        2*rhodsol*rhogeq*Kdrag*w1i**3*w3i*w2r**2*w3r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w3r*w2i**3 +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w2i**2*w3i*w2r + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w2i*w3r**3 +&
        4*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w2i*w2r*w3r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w2i*w3r*w3i**2 +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w2i*w3r*w2r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w3i*w2r*w3r**2 +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w3i**3*w2r + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w3i*w2r**3 +&
        4*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w3r*w3i*w2r**2 - rhodeq*vdsol*k*Kdrag**2*w2i**3*w3r*w3i**2 -&
        rhodeq*vdsol*k*Kdrag**2*w3i*w2r*w3r**2*w2i**2 - rhodeq*vdsol*k*Kdrag**2*w2r*w3i**3*w2i**2 -&
        rhodeq*vdsol*k*Kdrag**2*w2i*w2r**2*w3i**2*w3r - rhodeq*vdsol*k*Kdrag**2*w2i*w2r**2*w3r**3 -&
        rhodeq*vdsol*k*Kdrag**2*w3i**3*w2r**3 - rhodeq*vdsol*k*Kdrag**2*w3i*w2r**3*w3r**2 -&
        rhodeq*rhogsol*Kdrag*w1r*w2i**4*w3r**3 - rhodeq*rhogsol*Kdrag*w1r*w2i**4*w3r*w3i**2 -&
        2*rhodeq*rhogsol*Kdrag*w1r*w2r**2*w2i**2*w3i**2*w3r - rhodeq*rhogsol*Kdrag*w1r*w2i**2*w3r**4*w2r -&
        rhodeq*rhogsol*Kdrag*w1r*w2i**2*w3i**4*w2r - 2*rhodeq*rhogsol*Kdrag*w1r*w2r**2*w2i**2*w3r**3 -&
        2*rhodeq*rhogsol*Kdrag*w1r*w2i**2*w3r**2*w3i**2*w2r - rhodeq*rhogsol*Kdrag*w1r*w2r**4*w3r**3 -&
        rhodeq*rhogsol*Kdrag*w1r*w2r**3*w3r**4 - rhodeq*rhogsol*Kdrag*w1r*w2r**4*w3i**2*w3r -&
        rhodeq*rhogsol*Kdrag*w1r*w2r**3*w3i**4 - 2*rhodeq*rhogsol*Kdrag*w1r*w2r**3*w3i**2*w3r**2 -&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i**3*w3r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i**3*w3i**2 +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w2i**2*w3i**3 + rhodeq*vgsol*k*rhogeq**2*w1r**3*w2i*w3i**2*w2r**2 +&
        rhodeq*vgsol*k*rhogeq**2*w1r**3*w2i*w2r**2*w3r**2 + rhodeq*vgsol*k*rhogeq**2*w1r**3*w3i**3*w2r**2 +&
        rhodeq*vgsol*k*rhogeq**2*w1r**3*w3i*w2r**2*w3r**2 + 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w3i*w2i**3 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r*w2i**2*w2r*w3r - rhodsol*rhogeq*Kdrag*w1r**4*w2r**2*w3r**2 -&
        rhodeq*w1i*vgsol*k*Kdrag**2*w2i**2*w2r*w3r**2 - rhodeq*w1i*vgsol*k*Kdrag**2*w2i**2*w3r*w3i**2 -&
        rhodeq*w1i*vgsol*k*Kdrag**2*w2i**2*w2r*w3i**2 - rhodeq*w1i*vgsol*k*Kdrag**2*w3i**2*w2r**3 -&
        rhodeq*w1i*vgsol*k*Kdrag**2*w2r**3*w3r**2 - rhodeq*w1i*vgsol*k*Kdrag**2*w2r**2*w3r**3 -&
        rhodeq*w1i*vgsol*k*Kdrag**2*w2r**2*w3r*w3i**2 + rhodeq*rhogsol*Kdrag*w1r**4*w3i**2*w2i**2 +&
        rhodeq*rhogsol*Kdrag*w1r**4*w2i**2*w3r**2 + rhodeq*rhogsol*Kdrag*w1r**4*w3i**2*w2r**2 -&
        3*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i*w2i*w2r**2*w3r**2 - 3*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**3*w2i*w2r**2 -&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w2r**3*w3r**3 - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**2*w2r**4 -&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3r**2*w2r**4 - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**4*w2r**2 -&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**2*w2r**2*w3r**2 - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w2r**3*w3r*w3i**2 -&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3r**4*w2r**2 - rhodsol*rhogeq*Kdrag*w1r**4*w2i**2*w3r**2 -&
        rhodsol*rhogeq*Kdrag*w1r**4*w3i**2*w2r**2 - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3r**2*w2i**4 -&
        3*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**3*w2i**3 - 3*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i*w2i**3*w3r**2 -&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**2*w2i**2*w3r**2 - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w2r*w3r**3*w2i**2 -&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3r**4*w2i**2 - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3r**2*w2i**2*w2r**2 -&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w2r*w3r*w2i**2*w3i**2 - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**2*w2i**2*w2r**2&
        - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**2*w2i**4 - rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w3i**4*w2i**2 +&
        rhodeq*rhogsol*Kdrag*w1r**4*w2r**2*w3r**2 - 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w1i**2*w3r*w2i -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w1i**2*w3i*w2r + 2*rhodsol*rhogeq*Kdrag*w1i*w3i**3*w2i**4 +&
        2*rhodsol*rhogeq*Kdrag*w1i*w2i**4*w3i*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1i*w2i**3*w3r**4 +&
        2*rhodsol*rhogeq*Kdrag*w1i*w3i**4*w2i**3 + 4*rhodsol*rhogeq*Kdrag*w1i*w3i**2*w2i**3*w3r**2 +&
        4*rhodsol*rhogeq*Kdrag*w1i*w3i**3*w2i**2*w2r**2 + 4*rhodsol*rhogeq*Kdrag*w1i*w2i**2*w2r**2*w3r**2*w3i +&
        2*rhodsol*rhogeq*Kdrag*w1i*w2r**2*w3i**4*w2i + 2*rhodsol*rhogeq*Kdrag*w1i*w2i*w2r**2*w3r**4 +&
        4*rhodsol*rhogeq*Kdrag*w1i*w2r**2*w3i**2*w3r**2*w2i + 2*rhodsol*rhogeq*Kdrag*w1i*w2r**4*w3r**2*w3i +&
        2*rhodsol*rhogeq*Kdrag*w1i*w2r**4*w3i**3 + 2*rhodeq*w1i**2*rhogsol*Kdrag*w3i**2*w2i**2*w2r**2 +&
        3*rhodeq*w1i**2*rhogsol*Kdrag*w3i*w2i*w2r**2*w3r**2 + 3*rhodeq*w1i**2*rhogsol*Kdrag*w3i**3*w2i*w2r**2 +&
        rhodeq*w1i**2*rhogsol*Kdrag*w2r**3*w3r**3 + rhodeq*w1i**2*rhogsol*Kdrag*w3i**2*w2r**4 +&
        rhodeq*w1i**2*rhogsol*Kdrag*w3r**2*w2r**4 + rhodeq*w1i**2*rhogsol*Kdrag*w3i**4*w2r**2 +&
        2*rhodeq*w1i**2*rhogsol*Kdrag*w3i**2*w2r**2*w3r**2 + rhodeq*w1i**2*rhogsol*Kdrag*w2r**3*w3r*w3i**2 +&
        rhodeq*w1i**2*rhogsol*Kdrag*w3r**4*w2r**2 + rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**4*w2i*w3r**2 +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**4*w3i*w2r**2 + rhodeq*w1i**2*rhogsol*Kdrag*w3r**2*w2i**4 +&
        3*rhodeq*w1i**2*rhogsol*Kdrag*w3i**3*w2i**3 + 3*rhodeq*w1i**2*rhogsol*Kdrag*w3i*w2i**3*w3r**2 +&
        2*rhodeq*w1i**2*rhogsol*Kdrag*w3i**2*w2i**2*w3r**2 + rhodeq*w1i**2*rhogsol*Kdrag*w2r*w3r**3*w2i**2 +&
        rhodeq*w1i**2*rhogsol*Kdrag*w3r**4*w2i**2 + 2*rhodeq*w1i**2*rhogsol*Kdrag*w3r**2*w2i**2*w2r**2 +&
        rhodeq*w1i**2*rhogsol*Kdrag*w3i**4*w2i**2 - rhodeq*w1i*rhogsol*Kdrag*w3i**4*w2i**3 -&
        2*rhodeq*w1i*rhogsol*Kdrag*w3i**2*w2i**3*w3r**2 - 2*rhodeq*w1i*rhogsol*Kdrag*w3i**3*w2i**2*w2r**2 -&
        2*rhodeq*w1i*rhogsol*Kdrag*w2i**2*w2r**2*w3r**2*w3i - rhodeq*w1i*rhogsol*Kdrag*w2r**2*w3i**4*w2i -&
        rhodeq*w1i*rhogsol*Kdrag*w2i*w2r**2*w3r**4 - 2*rhodeq*w1i*rhogsol*Kdrag*w2r**2*w3i**2*w3r**2*w2i -&
        rhodeq*w1i*rhogsol*Kdrag*w2r**4*w3r**2*w3i - rhodeq*w1i*rhogsol*Kdrag*w2r**4*w3i**3 +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**4*w3i*w2i**2 + rhodeq*w1i**2*rhogsol*Kdrag*w3i**2*w2i**4 -&
        rhodeq*w1i**4*vgsol*cs**2*k**3*rhogeq**2*w3i*w2r - 2*rhodeq*w1i**3*rhogsol*Kdrag*w2i**3*w3i**2 -&
        2*rhodeq*w1i**3*rhogsol*Kdrag*w2i**3*w3r**2 - 2*rhodeq*w1i**3*rhogsol*Kdrag*w2i**2*w3i**3 -&
        2*rhodeq*w1i**3*rhogsol*Kdrag*w2i**2*w3i*w3r**2 - 2*rhodeq*w1i**3*rhogsol*Kdrag*w2i*w3i**2*w2r**2 -&
        2*rhodeq*w1i**3*rhogsol*Kdrag*w2i*w2r**2*w3r**2 - 2*rhodeq*w1i**3*rhogsol*Kdrag*w3i**3*w2r**2 -&
        2*rhodeq*w1i**3*rhogsol*Kdrag*w3i*w2r**2*w3r**2 + rhodeq*w1i**4*rhogsol*Kdrag*w3i**2*w2i**2 +&
        rhodeq*w1i**4*rhogsol*Kdrag*w2i**2*w3r**2 + rhodeq*w1i**4*rhogsol*Kdrag*w3i**2*w2r**2 +&
        rhodeq*w1i**4*rhogsol*Kdrag*w2r**2*w3r**2 - rhodeq*w1i*rhogsol*Kdrag*w3i**3*w2i**4 -&
        rhodeq*w1i*rhogsol*Kdrag*w2i**4*w3i*w3r**2 + rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3r**2*w2r**4 +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**4*w2r**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2r**2*w3r**2 +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3r**4*w2r**2 + rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3r**2*w2i**4 -&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2i**4 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2i**2*w3r**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**2*w2i**2*w2r**2 -&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3i**4*w2i**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3r**2*w2i**2*w2r**2 -&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w3r**4*w2i**2 + 2*rhodeq*rhogsol*Kdrag*w1r**2*w1i**2*w3i**2*w2i**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r**2*w1i**2*w2i**2*w3r**2 + 2*rhodeq*rhogsol*Kdrag*w1r**2*w1i**2*w3i**2*w2r**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r**2*w1i**2*w2r**2*w3r**2 + rhodeq*vdsol*k*Kdrag**2*w1r*w3i*w2r**2*w3r**2 -&
        rhodeq*w1i**3*vgsol*k*rhogeq**2*w2i**2*w2r*w3r**2 - rhodeq*w1i**3*vgsol*k*rhogeq**2*w2i**2*w3r*w3i**2 -&
        rhodeq*w1i**3*vgsol*k*rhogeq**2*w2i**2*w2r*w3i**2 - rhodeq*w1i**3*vgsol*k*rhogeq**2*w2i**2*w3r**3 -&
        rhodeq*w1i**3*vgsol*k*rhogeq**2*w3i**2*w2r**3 - rhodeq*w1i**3*vgsol*k*rhogeq**2*w2r**3*w3r**2 -&
        rhodeq*w1i**3*vgsol*k*rhogeq**2*w2r**2*w3r**3 - rhodeq*w1i**3*vgsol*k*rhogeq**2*w2r**2*w3r*w3i**2 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**4*w3r*w2i - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**4*w3i*w2r -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i*w2i**3 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i**2*w3r**2 +&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i**2*w2r*w3r -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**2*w2i**2 -&
        4*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i*w3r*w3i*w2r -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i*w3i*w2r**2 -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i*w2i*w3r**2 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**3*w2i +&
        rhodeq*vgsol*k*rhogeq**2*w1i*w2r**3*w3i**4 + 2*rhodeq*vgsol*k*rhogeq**2*w1i*w2r**3*w3i**2*w3r**2 +&
        rhodeq*vdsol*k*Kdrag**2*w1r*w2i**3*w3i**2 + rhodeq*vdsol*k*Kdrag**2*w1r*w2i**3*w3r**2 +&
        rhodeq*vdsol*k*Kdrag**2*w1r*w2i**2*w3i**3 + rhodeq*vdsol*k*Kdrag**2*w1r*w2i**2*w3i*w3r**2 +&
        rhodeq*vdsol*k*Kdrag**2*w1r*w2i*w3i**2*w2r**2 + rhodeq*vdsol*k*Kdrag**2*w1r*w3i**3*w2r**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i*w3i**3*w2r**3 + 2*rhodeq*rhogsol*Kdrag*w1r*w1i*w3i*w2r**3*w3r**2 +&
        rhodeq*vgsol*k*rhogeq**2*w1i*w2i**4*w3r**3 + rhodeq*vgsol*k*rhogeq**2*w1i*w2i**4*w3r*w3i**2 +&
        2*rhodeq*vgsol*k*rhogeq**2*w1i*w2r**2*w2i**2*w3i**2*w3r + rhodeq*vgsol*k*rhogeq**2*w1i*w2i**2*w3r**4*w2r +&
        rhodeq*vgsol*k*rhogeq**2*w1i*w2i**2*w3i**4*w2r + 2*rhodeq*vgsol*k*rhogeq**2*w1i*w2r**2*w2i**2*w3r**3 +&
        2*rhodeq*vgsol*k*rhogeq**2*w1i*w2i**2*w3r**2*w3i**2*w2r + rhodeq*vgsol*k*rhogeq**2*w1i*w2r**4*w3r**3 +&
        rhodeq*vgsol*k*rhogeq**2*w1i*w2r**3*w3r**4 - rhodeq*rhogeq*rhogsol*w1r**3*w2i**3*w3r*w3i**2 -&
        rhodeq*rhogeq*rhogsol*w1r**3*w3i*w2r*w3r**2*w2i**2 - rhodeq*rhogeq*rhogsol*w1r**3*w2r*w3i**3*w2i**2 -&
        rhodeq*rhogeq*rhogsol*w1r**3*w2i*w2r**2*w3i**2*w3r - rhodeq*rhogeq*rhogsol*w1r**3*w2i*w2r**2*w3r**3 -&
        rhodeq*rhogeq*rhogsol*w1r**3*w3i**3*w2r**3 - rhodeq*rhogeq*rhogsol*w1r**3*w3i*w2r**3*w3r**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i*w2i**3*w3r**3 + 2*rhodeq*rhogsol*Kdrag*w1r*w1i*w2i**3*w3r*w3i**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i*w3i*w2r*w3r**2*w2i**2 + 2*rhodeq*rhogsol*Kdrag*w1r*w1i*w2r*w3i**3*w2i**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i*w2i*w2r**2*w3i**2*w3r - rhodeq*rhogeq*rhogsol*w1r*w1i**2*w2r*w3i**3*w2i**2 -&
        rhodeq*rhogeq*rhogsol*w1r*w1i**2*w2i*w2r**2*w3i**2*w3r - rhodeq*rhogeq*rhogsol*w1r*w1i**2*w2i*w2r**2*w3r**3 -&
        rhodeq*rhogeq*rhogsol*w1r*w1i**2*w3i**3*w2r**3 - rhodeq*rhogeq*rhogsol*w1r*w1i**2*w3i*w2r**3*w3r**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w2i**3*w3i**2 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w2i**3*w3r**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w2i**2*w3i**3 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w2i**2*w3i*w3r**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w2i*w3i**2*w2r**2 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w2i*w2r**2*w3r**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w3i**3*w2r**2 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i*w3i*w2r**2*w3r**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r*w3r*w3i*w2i**4 -&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r*w2i**2*w3r*w3i*w2r**2 -&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r*w2i*w2r*w3r**2*w3i**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r*w2i*w3r**4*w2r - 2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r*w3i**4*w2r*w2i -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r*w3i*w2r**4*w3r - rhodeq*rhogeq*rhogsol*w1r*w1i**2*w2i**3*w3r**3 -&
        rhodeq*rhogeq*rhogsol*w1r*w1i**2*w3i*w2r*w3r**2*w2i**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w1r**2*w3r*w2i +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w1r**2*w3i*w2r - 2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w2i**3*w3i**2 -&
        2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w2i**3*w3r**2 - 2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w2i**2*w3i**3 -&
        2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w2i**2*w3i*w3r**2 - 2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w2i*w3i**2*w2r**2 -&
        2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w2i*w2r**2*w3r**2 - 2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w3i**3*w2r**2 -&
        2*rhodeq*w1i*rhogsol*Kdrag*w1r**2*w3i*w2r**2*w3r**2 + 2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w2i*w3r*w3i&
        + 2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w2i*w3i*w2r -&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w2r*w3r**2 -&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w3r*w2r**2 -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w1r**2*w3r*w2i -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w1r**2*w3i*w2r + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2i**2*w2r*w3r**2&
        + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2i**2*w3r*w3i**2 + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2i**2*w2r*w3i**2&
        + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2i**2*w3r**3 + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w3i**2*w2r**3 +&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2r**3*w3r**2 + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2r**2*w3r**3 +&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r*w2r**2*w3r*w3i**2 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w3r*w3i**2&
        - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w3r*w2r**2 -&
        4*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w2r*w3i**2 -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r*w3r**2 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i**3*w2r -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r**3 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w1i**2*w3i*w2i**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w1i**2*w2i*w3i**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w1i**2*w2i*w3r**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w1i**2*w3i*w2r**2 +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3i**3*w2r + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3i*w2r*w3r**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3r*w3i*w2r**2 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3r*w2i**3&
        - 4*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i**2*w3r*w3i -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i**2*w3i*w2r + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3r*w2i**3 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i**2*w3r*w3i -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i**2*w3i*w2r +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i*w3r*w2r**2 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i*w3r*w3i**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i*w2r*w3r**2 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i*w3r**3&
        - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w3i*w2r**3 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w2i*w2r*w3i**2 - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w3r**3&
        - rhodeq*rhogeq*rhogsol*w1r*w1i**2*w2i**3*w3r*w3i**2 - rhodeq*rhogeq*rhogsol*w1r**3*w2i**3*w3r**3 +&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i*w2i*w2r**2*w3r**3 + rhodeq*vgsol*k*rhogeq**2*w1i*w2r**4*w3i**2*w3r +&
        rhodeq*vdsol*k*Kdrag**2*w1r*w2i*w2r**2*w3r**2 + rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**2*w3r*w2r +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i**2*w2i*w3i*w2r + rhodeq*w1i*rhogsol*cs**4*k**4*rhogeq*w1r**2*w2r*w3r&
        - 2*rhodeq*w1i**2*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i*w2r + rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w2i**3*w3r*w3i**2&
        + rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i**2*w2r**2*w3r**2 +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w2i*w2r**2*w3r**2 -&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w3r*w3i**2 + rhodeq*rhogeq*rhogsol*w1r**2*w2i*w2r**2*w3r**4 -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r**3*w3i*w2i + rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r**3*w2r*w3r +&
        rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**4 - 2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**3*w3r**2 +&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**3*w3i**2 +&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i**3 +&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i*w3r**2 +&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i**2*w3i*w2r**2 +&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2*w3r**2 +&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2*w2r**2 +&
        rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**4 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w2r**2*w3r**2 +&
        rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w2i*w3r**4 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i*w2r**2*w3r**2 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i**3*w2r**2 + rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w3i*w2r**4 -&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3i*w2i**3 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w2i**2*w2r*w3r +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w2i**2*w3r**2 - 3*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3i**2*w2i**2 +&
        4*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w2i*w3r*w3i*w2r + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i**2*w2i**2*w3r**2 +&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i**2*w3i**2*w2r**2 + rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i**2*w2r**2*w3r**2 +&
        rhodeq*w1i**4*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**2 + rhodeq*w1i**4*rhogsol*cs**2*k**2*rhogeq*w2i*w3i**2 +&
        rhodeq*w1i**4*rhogsol*cs**2*k**2*rhogeq*w2i*w3r**2 + rhodeq*w1i**4*rhogsol*cs**2*k**2*rhogeq*w3i*w2r**2 +&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1i**2*w3i**2*w2i**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i*w3i*w2i**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i*w2i*w3i**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i*w2i*w3r**2 +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i*w3i*w2r**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1r**3*w3i**2*w2i**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r**3*w2i**2*w3r**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1r**3*w3i**2*w2r**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r**3*w2r**2*w3r**2 + rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w2i**3*w3i**2 +&
        rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w2i**3*w3r**2 + rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w2i**2*w3i**3 +&
        rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w2i**2*w3i*w3r**2 + rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w2i*w3i**2*w2r**2 +&
        rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w2i*w2r**2*w3r**2 + rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w3i**3*w2r**2 +&
        rhodeq*w1i**2*vgsol*k*rhogeq**2*w1r*w3i*w2r**2*w3r**2 - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w2r*w3i**2&
        - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w3i**2*w2r**3 - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2r**3*w3r**2 -&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2r**2*w3r**3 - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2r**2*w3r*w3i**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r**3*w3i*w2i - rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r**3*w2r*w3r -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r**2*w2i*w3r**2 - 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r**2*w2i*w2r*w3r -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r**2*w3i*w2r**2 - 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r**2*w3r*w3i*w2r -&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w2r*w3r**2 - 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w3r**3 -&
        2*rhodeq*w1i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3i*w2i**3 -&
        2*rhodeq*w1i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3i**2*w2i**2 -&
        2*rhodeq*w1i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i*w3i*w2r**2 -&
        2*rhodeq*w1i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3i*w2i*w3r**2 -&
        2*rhodeq*w1i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3i**3*w2i +&
        2*rhodeq*w1i*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2r**2*w3r**2 + rhodeq*vgsol*Kdrag*k*rhogeq*w1r**3*w3i**2*w2i**2&
        + rhodeq*vgsol*Kdrag*k*rhogeq*w1r**3*w2i**2*w3r**2 + rhodeq*vgsol*Kdrag*k*rhogeq*w1r**3*w3i**2*w2r**2 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r**3*w2r**2*w3r**2 + 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w3i**3*w2r**2 +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w3i*w2r**2*w3r**2 -&
        2*rhodeq*w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i*w2i**2 -&
        2*rhodeq*w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i*w3i**2 -&
        2*rhodeq*w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i*w3r**2 -&
        2*rhodeq*w1i*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i*w2r**2 - rhodeq*vgsol*k*rhogeq**2*w1r*w3i**3*w2i**4 -&
        rhodeq*vgsol*k*rhogeq**2*w1r*w2i**4*w3i*w3r**2 - rhodeq*vgsol*k*rhogeq**2*w1r*w2i**3*w3r**4 -&
        rhodeq*vgsol*k*rhogeq**2*w1r*w3i**4*w2i**3 - 2*rhodeq*vgsol*k*rhogeq**2*w1r*w3i**2*w2i**3*w3r**2 -&
        2*rhodeq*vgsol*k*rhogeq**2*w1r*w3i**3*w2i**2*w2r**2 - 2*rhodeq*vgsol*k*rhogeq**2*w1r*w2i**2*w2r**2*w3r**2*w3i -&
        rhodeq*vgsol*k*rhogeq**2*w1r*w2r**2*w3i**4*w2i - rhodeq*vgsol*k*rhogeq**2*w1r*w2i*w2r**2*w3r**4 -&
        2*rhodeq*vgsol*k*rhogeq**2*w1r*w2r**2*w3i**2*w3r**2*w2i - rhodeq*vgsol*k*rhogeq**2*w1r*w2r**4*w3r**2*w3i -&
        rhodeq*vgsol*k*rhogeq**2*w1r*w2r**4*w3i**3 + 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w2i**3*w3i**2 +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w2i**3*w3r**2 + 2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w2i**2*w3i**3 +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w2i**2*w3i*w3r**2 +&
        2*rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i*w2i*w3i**2*w2r**2 + rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w2r*w3i**3*w2i**2 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w2i*w2r**2*w3i**2*w3r + rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w2i*w2r**2*w3r**3 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w3i**3*w2r**3 + rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w3i*w2r**3*w3r**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i**3*w3r*w2i + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i**3*w3i*w2r +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i**2*w3i**2*w2i**2 + rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i**2*w2i**2*w3r**2 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1r*w1i**2*w3i**2*w2r**2 - 2*rhodeq*w1i**2*vdsol*Kdrag*cs**2*k**3*rhogeq*w2r*w3i**2&
        + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i**2*w3i*w2i -&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i**2*w2r*w3r + rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w2i**3*w3r**3 +&
        rhodeq*vgsol*Kdrag*k*rhogeq*w1i*w3i*w2r*w3r**2*w2i**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w1i**2*w3i**2*w2i**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w1i**2*w2i**2*w3r**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w1i**2*w3i**2*w2r**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w1i**2*w2r**2*w3r**2 - 2*rhodeq*w1i**2*vdsol*Kdrag*cs**2*k**3*rhogeq*w3r*w2i**2&
        - 2*rhodeq*w1i**2*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r*w3i +&
        2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i**2*w2r*w3i**2 + 2*rhodeq*vgsol*k*Kdrag**2*w1r*w1i*w3i**2*w2i**2 +&
        2*rhodeq*vgsol*k*Kdrag**2*w1r*w1i*w2i**2*w3r**2 + 2*rhodeq*vgsol*k*Kdrag**2*w1r*w1i*w3i**2*w2r**2 +&
        2*rhodeq*vgsol*k*Kdrag**2*w1r*w1i*w2r**2*w3r**2 - 4*rhodsol*rhogeq*Kdrag*w1r*w1i*w2i**3*w3r**3 -&
        4*rhodsol*rhogeq*Kdrag*w1r*w1i*w2i**3*w3r*w3i**2 - 4*rhodsol*rhogeq*Kdrag*w1r*w1i*w3i*w2r*w3r**2*w2i**2 -&
        4*rhodsol*rhogeq*Kdrag*w1r*w1i*w2r*w3i**3*w2i**2 - 4*rhodsol*rhogeq*Kdrag*w1r*w1i*w2i*w2r**2*w3i**2*w3r -&
        4*rhodsol*rhogeq*Kdrag*w1r*w1i*w2i*w2r**2*w3r**3 - 4*rhodsol*rhogeq*Kdrag*w1r*w1i*w3i**3*w2r**3 -&
        4*rhodsol*rhogeq*Kdrag*w1r*w1i*w3i*w2r**3*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w2i**3*w3i**2 +&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w2i**3*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w2i**2*w3i**3 +&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w2i**2*w3i*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w2i*w3i**2*w2r**2 +&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w2i*w2r**2*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w3i**3*w2r**2 +&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w1i*w3i*w2r**2*w3r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w1i**2*w3r*w2i +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w1i**2*w3i*w2r - rhodeq*w1i*rhogsol*cs**4*k**4*rhogeq*w1r**2*w3i*w2i +&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3r**3*w2r + rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2r**3*w3r +&
        3*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2r**2*w3r**2 -&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**2*w2r**2 + rhodeq*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w3r*w2i**2 +&
        2*rhodeq*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w2i*w3i*w2r + 2*rhodeq*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w2i*w3r*w3i&
        + rhodeq*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w2r*w3i**2 - rhodeq*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w2r*w3r**2 -&
        rhodeq*rhogeq**2*cs**2*k**3*w1i**3*vgsol*w3r*w2r**2 + 2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i**2*w3r*w2i**2&
        + 2*rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1i**2*w2i*w3r*w3i -&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3i*w2i*w3r**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2r**2*w3r*w3i**2 +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3i*w2i**4 - 2*rhodsol*rhogeq*Kdrag*w1i**2*w3i**2*w2i**2*w3r**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r**2*w3i**2*w2i**2 - rhodeq*vgsol*k*Kdrag**2*w1r*w2i**2*w3i*w3r**2 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2i**2*w2r*w3i**2 + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**2*w2i**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w2r*w3r*w2i**2*w3i**2 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w3r*w2r**4 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i*w3i**4 + rhodeq*rhogeq*rhogsol*w1r**2*w2r**4*w3r**2*w3i +&
        rhodeq*rhogeq*rhogsol*w1r**2*w2r**4*w3i**3 - 2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2i**2*w2r*w3r**2 -&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2i**2*w3r*w3i**2 - 2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2i**2*w2r*w3i**2 -&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2i**2*w3r**3 - 2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w3i**2*w2r**3 -&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2r**3*w3r**2 - 2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2r**2*w3r**3 -&
        2*rhodeq*rhogsol*Kdrag*w1r*w1i**2*w2r**2*w3r*w3i**2 - rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3i**3*w2i -&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w2i*w3i*w2r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w2r**3*w3r +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w2r**2*w3r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3r**3*w2r +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3i**2*w3r*w2r + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i*w3i**2*w2r**2 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w2r*w3r**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w3r*w3i**2 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w2r*w3i**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2i**2*w3r**3 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w3i**2*w2r**3 + rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2r**3*w3r**2 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i**2*w2r**2*w3r**3 + 4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i**2*w3r*w3i*w2r&
        + 2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i**2*w3i*w2r**2 +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i*w3i**4 + 4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i*w3r**3*w2r +&
        rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i*w3r**4 +&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i*w3i**2*w3r*w2r +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w2i*w3i**2*w3r**2 +&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3r*w3i*w2r**3 + rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**2*w3i*w2r**4 -&
        rhodeq*w1i**3*rhogeq*rhogsol*w3i*w2i**3*w3r**2 + rhodeq*w1i**3*rhogeq*rhogsol*w2r*w3r**3*w2i**2 +&
        rhodeq*w1i**3*rhogeq*rhogsol*w2r*w3r*w2i**2*w3i**2 - rhodeq*w1i**3*rhogeq*rhogsol*w3i**3*w2i*w2r**2 -&
        rhodeq*w1i**3*rhogeq*rhogsol*w3i*w2i*w2r**2*w3r**2 + rhodeq*w1i**3*rhogeq*rhogsol*w2r**3*w3r**3 +&
        rhodeq*w1i**3*rhogeq*rhogsol*w2r**3*w3r*w3i**2 - rhodsol*rhogeq*Kdrag*w1i**2*w3i**2*w2i**4 -&
        rhodsol*rhogeq*Kdrag*w1i**2*w3r**2*w2i**4 - 4*rhodsol*rhogeq*Kdrag*w1i**2*w3i**3*w2i**3 -&
        4*rhodsol*rhogeq*Kdrag*w1i**2*w3i*w2i**3*w3r**2 - rhodsol*rhogeq*Kdrag*w1i**2*w3i**4*w2i**2 -&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**3*w3r*w3i - rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w3r*w3i**2&
        - rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w2r*w3i**2 - rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w3r**3 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i**2*w2r*w3r**2 -&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3r*w3i*w2r**2 -&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i*w2r*w3r**2 -&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2i*w3i**3*w2r + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**3*w3r**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w3r**3 + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w2r**2*w3r*w3i**2 -&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w3i**2*w2r**3 - 2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w2i**2*w3r*w3i -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w2i**2*w3i*w2r -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w2i*w2r*w3r**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w2i*w3r*w3i**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w2i*w2r*w3i**2 - 2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w2i*w3r**3&
        - 2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w3i*w2r**3 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1r**3*w3r*w3i*w2r**2 + rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r**3*w3r*w2i +&
        rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r**3*w3i*w2r - rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r**2*w2i**2*w3r**2 -&
        rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r**2*w3i**2*w2r**2 - rhodeq*rhogsol*cs**2*k**2*Kdrag*w1r**2*w2r**2*w3r**2 -&
        rhodeq*vgsol*k*Kdrag**2*w1r*w2i**3*w3i**2 - rhodeq*vgsol*k*Kdrag**2*w1r*w2i**3*w3r**2 -&
        rhodeq*vgsol*k*Kdrag**2*w1r*w2i**2*w3i**3 - 2*rhodsol*rhogeq*Kdrag*w1i**2*w3i**2*w2i**2*w2r**2 -&
        2*rhodsol*rhogeq*Kdrag*w1i**2*w3r**2*w2i**2*w2r**2 - rhodsol*rhogeq*Kdrag*w1i**2*w3r**4*w2i**2 -&
        4*rhodsol*rhogeq*Kdrag*w1i**2*w3i**3*w2i*w2r**2 - 4*rhodsol*rhogeq*Kdrag*w1i**2*w3i*w2i*w2r**2*w3r**2 -&
        rhodsol*rhogeq*Kdrag*w1i**2*w3i**2*w2r**4 - rhodsol*rhogeq*Kdrag*w1i**2*w3r**2*w2r**4 -&
        rhodsol*rhogeq*Kdrag*w1i**2*w3i**4*w2r**2 - 2*rhodsol*rhogeq*Kdrag*w1i**2*w3i**2*w2r**2*w3r**2 -&
        rhodsol*rhogeq*Kdrag*w1i**2*w3r**4*w2r**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2i**2*w3r**3 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w3i**2*w2r**3 + rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2r**3*w3r**2 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2r**2*w3r**3 + rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2r**2*w3r*w3i**2 +&
        rhodeq*w1i*vdsol*k*Kdrag**2*w2i**2*w2r*w3r**2 + rhodeq*w1i*vdsol*k*Kdrag**2*w2i**2*w3r*w3i**2 +&
        rhodeq*w1i*vdsol*k*Kdrag**2*w2i**2*w2r*w3i**2 + rhodeq*w1i*vdsol*k*Kdrag**2*w2i**2*w3r**3 +&
        rhodeq*w1i*vdsol*k*Kdrag**2*w3i**2*w2r**3 + rhodeq*w1i*vdsol*k*Kdrag**2*w2r**3*w3r**2 +&
        rhodeq*w1i*vdsol*k*Kdrag**2*w2r**2*w3r**3 + rhodeq*w1i*vdsol*k*Kdrag**2*w2r**2*w3r*w3i**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3r*w2i**3 + 4*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i**2*w3r*w3i&
        + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i**2*w3i*w2r + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w3r**3&
        + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w3r*w3i**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w3r*w2r**2 +&
        4*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w2i*w2r*w3i**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r*w3r**2 + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i**3*w2r +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1i*w3i*w2r**3 + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i*w2i**3 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i**2*w3r**2 -&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i**2*w2r*w3r - rhodeq*vgsol*k*Kdrag**2*w1r*w2i*w3i**2*w2r**2 -&
        rhodeq*vgsol*k*Kdrag**2*w1r*w2i*w2r**2*w3r**2 - rhodeq*vgsol*k*Kdrag**2*w1r*w3i**3*w2r**2 -&
        rhodeq*vgsol*k*Kdrag**2*w1r*w3i*w2r**2*w3r**2 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r**2*w3r*w2i**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r**2*w2i*w3i*w2r +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r**2*w2i*w3r*w3i +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r**2*w2r*w3i**2 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r**2*w2r*w3r**2 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w1r**2*w3r*w2r**2 -&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w2i*w3r*w3i -&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w2i*w3i*w2r +&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w2r*w3r**2 +&
        2*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r**2*w3r*w2r**2 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w1r*w3i*w2i**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w1r*w2i*w2r*w3r -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w1r*w2i*w3i**2 +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w1r*w2i*w3r**2 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w1r*w3r*w3i*w2r +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i**2*w1r*w3i*w2r**2 - 2*rhodeq*w1i**3*rhogsol*cs**2*k**2*rhogeq*w3i*w2i**3&
        - 2*rhodeq*w1i**3*rhogsol*cs**2*k**2*rhogeq*w3i**2*w2i**2 -&
        2*rhodeq*w1i**3*rhogsol*cs**2*k**2*rhogeq*w2i*w3i*w2r**2 -&
        2*rhodeq*w1i**3*rhogsol*cs**2*k**2*rhogeq*w3i*w2i*w3r**2 - 2*rhodeq*w1i**3*rhogsol*cs**2*k**2*rhogeq*w3i**3*w2i&
        + 2*rhodeq*w1i**3*rhogsol*cs**2*k**2*rhogeq*w2r**2*w3r**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w2i**3*w3r**3 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w2i**3*w3r*w3i**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w3i*w2r*w3r**2*w2i**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w2r*w3i**3*w2i**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w2i*w2r**2*w3i**2*w3r -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w2i*w2r**2*w3r**3 - rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w3i**3*w2r**3 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1i*w3i*w2r**3*w3r**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2i**2*w2r*w3r**2 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r**2*w2i**2*w3r*w3i**2 - rhodeq*w1i*rhogeq*rhogsol*w3r**4*w2i**4 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w3i**3*w2i*w2r**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w3i*w2i*w2r**2*w3r**2 -&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w2r**3*w3r**3 - rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w2r**3*w3r*w3i**2 -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w3r*w2i**4 - 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2i**3*w3r*w3i -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2i**2*w3r*w2r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2i**2*w2r*w3r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2i*w3r*w3i*w2r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2i*w3i*w2r*w3r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2i*w3i**3*w2r -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2r**2*w3r*w3i**2 +&
        4*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i*w3r*w3i*w2r +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2i*w3i*w2r**2 +&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i*w2i*w3r**2 + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**3*w2i -&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**2*w3r*w2r - rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3r**3*w2r -&
        rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2r**3*w3r - 3*rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w2r**2*w3r**2&
        + rhodeq*vdsol*Kdrag*cs**2*k**3*rhogeq*w1r*w3i**2*w2r**2 - rhodeq*w1i**3*rhogsol*cs**4*k**4*rhogeq*w3i*w2i +&
        rhodeq*w1i**3*rhogsol*cs**4*k**4*rhogeq*w2r*w3r - rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i**2*w3i*w2i +&
        rhodeq*vgsol*Kdrag*cs**2*k**3*rhogeq*w1r*w1i**2*w2r*w3r - rhodeq*w1i*rhogeq*rhogsol*w3i**4*w2i**4 -&
        2*rhodeq*w1i*rhogeq*rhogsol*w3r**2*w2i**4*w3i**2 - 2*rhodeq*w1i*rhogeq*rhogsol*w3i**4*w2i**2*w2r**2 -&
        2*rhodeq*w1i*rhogeq*rhogsol*w3r**4*w2i**2*w2r**2 - 4*rhodeq*w1i*rhogeq*rhogsol*w3i**2*w2i**2*w2r**2*w3r**2 -&
        2*rhodeq*w1i*rhogeq*rhogsol*w3i**2*w2r**4*w3r**2 - rhodeq*w1i*rhogeq*rhogsol*w3r**4*w2r**4 -&
        rhodeq*w1i*rhogeq*rhogsol*w3i**4*w2r**4 - 2*rhodeq*vdsol*k*Kdrag**2*w1r*w1i*w3i**2*w2i**2 -&
        2*rhodeq*vdsol*k*Kdrag**2*w1r*w1i*w2i**2*w3r**2 - 2*rhodeq*vdsol*k*Kdrag**2*w1r*w1i*w3i**2*w2r**2 -&
        2*rhodeq*vdsol*k*Kdrag**2*w1r*w1i*w2r**2*w3r**2 + rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w3i**3*w2i**3 +&
        rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w3i*w2i**3*w3r**2 - rhodeq*vdsol*Kdrag*k*rhogeq*w1r*w2r*w3r**3*w2i**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i*w3i**2*w3r*w2r -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i*w3r**3*w2r -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i*w3i**2*w2r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i*w3i**2*w3r**2 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w3i*w2r**4 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w3r*w3i*w2r**3 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w3r**4*w2r -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2r*w3i**4 - 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2r**3*w3r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2r**2*w3r**3 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1i*w2r*w3r**2*w3i**2 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w3i*w2i**4 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i**3*w3i**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i**2*w3i*w3r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i**2*w3i**3 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i**2*w3i*w2r**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i**2*w3r*w3i*w2r - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r*w2i*w3r**4 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3i**3*w2i**3 + rhodeq*rhogsol*Kdrag*w1r**2*w3i*w2i**3*w3r**2 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3r**2*w2i**4 + 3*rhodeq*rhogsol*Kdrag*w1r**2*w2r*w3r**3*w2i**2 -&
        rhodeq*w1i**3*rhogeq*rhogsol*w3i**3*w2i**3 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w2i*w3r**3*w2r**2 +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w2i*w3r*w3i**2*w2r**2 + rhodeq*rhogeq*rhogsol*w1r**2*w3i**3*w2i**4 +&
        rhodeq*rhogeq*rhogsol*w1r**2*w2i**4*w3i*w3r**2 + rhodeq*rhogeq*rhogsol*w1r**2*w2i**3*w3r**4 +&
        rhodeq*rhogeq*rhogsol*w1r**2*w3i**4*w2i**3 + 2*rhodeq*rhogeq*rhogsol*w1r**2*w3i**2*w2i**3*w3r**2 +&
        2*rhodeq*rhogeq*rhogsol*w1r**2*w3i**3*w2i**2*w2r**2 + 2*rhodeq*rhogeq*rhogsol*w1r**2*w2i**2*w2r**2*w3r**2*w3i +&
        rhodeq*rhogeq*rhogsol*w1r**2*w2r**2*w3i**4*w2i + 2*rhodeq*rhogeq*rhogsol*w1r**2*w2r**2*w3i**2*w3r**2*w2i +&
        4*rhodsol*rhogeq*Kdrag*w1r*w2r**2*w2i**2*w3r**3 + 2*rhodsol*rhogeq*Kdrag*w1r*w2r**3*w3i**4 +&
        4*rhodsol*rhogeq*Kdrag*w1r*w2r**3*w3i**2*w3r**2 - 2*rhodsol*rhogeq*Kdrag*w1r**2*w1i**2*w3i**2*w2i**2 -&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w1i**2*w2i**2*w3r**2 - 2*rhodsol*rhogeq*Kdrag*w1r**2*w1i**2*w3i**2*w2r**2 -&
        2*rhodsol*rhogeq*Kdrag*w1r**2*w1i**2*w2r**2*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1r*w2i**4*w3r**3 -&
        2*rhodeq*rhogsol*Kdrag*w1r**3*w2r**2*w3r*w3i**2 + 2*rhodsol*rhogeq*Kdrag*w1r*w2i**4*w3r*w3i**2 +&
        4*rhodsol*rhogeq*Kdrag*w1r*w2r**2*w2i**2*w3i**2*w3r + 2*rhodsol*rhogeq*Kdrag*w1r*w2i**2*w3r**4*w2r +&
        2*rhodsol*rhogeq*Kdrag*w1r*w2i**2*w3i**4*w2r - 2*rhodeq*rhogsol*Kdrag*w1r**3*w2i**2*w2r*w3i**2 -&
        2*rhodeq*rhogsol*Kdrag*w1r**3*w3i**2*w2r**3 - 2*rhodeq*rhogsol*Kdrag*w1r**3*w2r**3*w3r**2 -&
        2*rhodeq*rhogsol*Kdrag*w1r**3*w2r**2*w3r**3 + rhodeq*rhogeq*rhogsol*w1r**2*w1i*w2r**3*w3r*w3i**2 -&
        rhodeq*w1i**3*vgsol*Kdrag*cs**2*k**3*rhogeq*w3r*w2i - rhodeq*w1i**3*vgsol*Kdrag*cs**2*k**3*rhogeq*w3i*w2r -&
        2*rhodeq*rhogsol*Kdrag*w1r**3*w2i**2*w2r*w3r**2 - 2*rhodeq*rhogsol*Kdrag*w1r**3*w2i**2*w3r*w3i**2 -&
        2*rhodeq*rhogsol*Kdrag*w1r**3*w2i**2*w3r**3 - rhodeq*rhogeq*rhogsol*w1r**2*w1i*w3i*w2i**3*w3r**2 +&
        rhodeq*rhogeq*rhogsol*w1r**2*w1i*w2r*w3r**3*w2i**2 + rhodeq*rhogeq*rhogsol*w1r**2*w1i*w2r*w3r*w2i**2*w3i**2 -&
        rhodeq*rhogeq*rhogsol*w1r**2*w1i*w3i**3*w2i*w2r**2 - rhodeq*rhogeq*rhogsol*w1r**2*w1i*w3i*w2i*w2r**2*w3r**2 +&
        rhodeq*rhogeq*rhogsol*w1r**2*w1i*w2r**3*w3r**3 - 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i**2*w2i*w2r*w3r -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i**2*w3r*w3i*w2r -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w3r*w3i*w2r**2 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w3i*w2r**3 + 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i**2*w3i*w2i**2&
        + 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1i**2*w2i*w3i**2 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w2i*w3r*w3i**2 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w2i*w2r*w3i**2 -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w1i*w2r*w3i**2 -&
        2*rhodeq*w1i**2*rhogsol*cs**2*k**2*rhogeq*w1r*w2i*w3r**3 -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w1i*w3r*w2i**2 - 2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w1i*w3r*w2r**2 -&
        2*rhodeq*rhogsol*cs**4*k**4*rhogeq*w1r*w1i*w2r*w3r**2 + rhodeq*w1i**2*rhogeq*rhogsol*w2i*w2r**2*w3r**4 +&
        2*rhodeq*w1i**2*rhogeq*rhogsol*w2r**2*w3i**2*w3r**2*w2i + rhodeq*w1i**2*rhogeq*rhogsol*w2r**4*w3r**2*w3i +&
        2*rhodeq*w1i**2*rhogeq*rhogsol*w3i**3*w2i**2*w2r**2 + 2*rhodeq*w1i**2*rhogeq*rhogsol*w2i**2*w2r**2*w3r**2*w3i +&
        2*rhodeq*w1i**2*rhogeq*rhogsol*w3i**2*w2i**3*w3r**2 + rhodeq*w1i**2*rhogeq*rhogsol*w2r**2*w3i**4*w2i +&
        rhodeq*w1i**2*rhogeq*rhogsol*w2r**4*w3i**3 + 2*rhodeq*rhogsol*Kdrag*w1r**2*w3r**2*w2i**2*w2r**2 +&
        2*rhodeq*rhogsol*Kdrag*w1r**2*w3i**2*w2i**2*w2r**2 - 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w3r*w3i*w2r**2&
        + rhodeq*w1i**2*rhogeq*rhogsol*w2i**4*w3i*w3r**2 + rhodeq*w1i**2*rhogeq*rhogsol*w2i**3*w3r**4 +&
        rhodeq*w1i**2*rhogeq*rhogsol*w3i**4*w2i**3 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i*w3r**3 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i*w2r*w3i**2 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w3i*w2r**3&
        - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w3i**3*w2r - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w3i*w2r*w3r**2&
        + rhodeq*w1i**2*rhogeq*rhogsol*w3i**3*w2i**4 - rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w3r*w2i**3 +&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i**2*w3r*w3i +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i**2*w3i*w2r -&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i*w3r*w2r**2 +&
        rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i*w3r*w3i**2 -&
        2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w1r**2*w2i*w2r*w3r**2 + 4*rhodsol*rhogeq*Kdrag*w1r*w2i**2*w3r**2*w3i**2*w2r&
        + 2*rhodsol*rhogeq*Kdrag*w1r*w2r**4*w3r**3 + 2*rhodsol*rhogeq*Kdrag*w1r*w2r**3*w3r**4 +&
        2*rhodsol*rhogeq*Kdrag*w1r*w2r**4*w3i**2*w3r + 2*rhodsol*rhogeq*Kdrag*w1r**3*w2i**2*w3r**3 +&
        2*rhodsol*rhogeq*Kdrag*w1r**3*w2r**2*w3r**3 + 2*rhodsol*rhogeq*Kdrag*w1r**3*w2r**2*w3r*w3i**2 +&
        2*rhodsol*rhogeq*Kdrag*w1r**3*w2i**2*w2r*w3r**2 + 2*rhodsol*rhogeq*Kdrag*w1r**3*w2i**2*w3r*w3i**2 +&
        2*rhodsol*rhogeq*Kdrag*w1r**3*w2i**2*w2r*w3i**2 + 3*rhodeq*rhogsol*Kdrag*w1r**2*w2r**3*w3r**3 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3i**2*w2r**4 + rhodeq*rhogsol*Kdrag*w1r**2*w3r**2*w2r**4 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3i**4*w2r**2 + 3*rhodeq*rhogsol*Kdrag*w1r**2*w2r**3*w3r*w3i**2 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3r**4*w2r**2 + 2*rhodeq*rhogsol*Kdrag*w1r**2*w3i**2*w2r**2*w3r**2 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3i**4*w2i**2 + 2*rhodeq*rhogsol*Kdrag*w1r**2*w3i**2*w2i**2*w3r**2 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3i**3*w2i*w2r**2 + rhodeq*rhogsol*Kdrag*w1r**2*w3i*w2i*w2r**2*w3r**2 +&
        3*rhodeq*rhogsol*Kdrag*w1r**2*w2r*w3r*w2i**2*w3i**2 + rhodeq*rhogsol*Kdrag*w1r**2*w3r**4*w2i**2 +&
        rhodeq*rhogsol*Kdrag*w1r**2*w3i**2*w2i**4 + 2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w3i**2*w2r**3 -&
        rhodsol*rhogeq*Kdrag*w1i**4*w3i**2*w2r**2 + 2*rhodsol*rhogeq*Kdrag*w1r**3*w3i**2*w2r**3 +&
        2*rhodsol*rhogeq*Kdrag*w1r**3*w2r**3*w3r**2 + 4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i*w3i*w2r*w3r**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2r**2*w3r**3 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2r**2*w3r*w3i**2 - rhodsol*rhogeq*Kdrag*w1i**4*w3i**2*w2i**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i**2*w3r*w3i**2 +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i**2*w2r*w3i**2 +&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i**3*w3r*w3i +&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i**2*w3r**3 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i**2*w2r*w3r**2 +&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i*w3r*w3i*w2r**2 -&
        2*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2r**3*w3r**2 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w3i**3*w2r*w2i**2&
        + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w3i*w2r*w3r**2*w2i**2 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w3i*w2r**3*w3r**2&
        + 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w2i*w3r**4*w2r + 2*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w3i**4*w2r*w2i +&
        4*rhodeq*vgsol*cs**2*k**3*rhogeq**2*w3i**2*w2r*w3r**2*w2i - rhodsol*rhogeq*Kdrag*w3i**4*w2r**4 -&
        2*rhodsol*rhogeq*Kdrag*w3i**4*w2r**2*w2i**2 - 2*rhodsol*rhogeq*Kdrag*w2i**4*w3i**2*w3r**2 -&
        rhodsol*rhogeq*Kdrag*w2r**4*w3r**4 - 2*rhodsol*rhogeq*Kdrag*w2r**4*w3r**2*w3i**2 +&
        4*rhodeq*rhogeq**2*w3r*cs**2*k**3*w3i*vgsol*w2r**2*w2i**2 + 2*rhodeq*rhogeq**2*w3r*cs**2*k**3*w3i*vgsol*w2r**4&
        - 2*rhodsol*rhogeq*Kdrag*w3r**4*w2i**2*w2r**2 + rhodeq*vgsol*cs**2*k**3*rhogeq**2*w2i**3*w3r*w3i**2 +&
        2*rhodeq*rhogeq**2*w3r*cs**2*k**3*w3i*vgsol*w2i**4 - 4*rhodsol*rhogeq*Kdrag*w3i**2*w2r**2*w3r**2*w2i**2 +&
        4*rhodeq*rhogsol*cs**2*k**2*rhogeq*w1i*w1r*w2i*w3i**3*w2r - rhodsol*rhogeq*Kdrag*w1i**4*w2i**2*w3r**2 -&
        rhodsol*rhogeq*Kdrag*w1i**4*w2r**2*w3r**2 - rhodsol*rhogeq*Kdrag*w1r**2*w3i**2*w2i**4)/( - 2*w1i*w3i + w1i**2 +&
        w3r**2 + w1r**2 + w3i**2 - 2*w1r*w3r)/(w3i**2 + w3r**2)/Kdrag/(w2i**2 + w2r**2)/(w2i**2 + w1r**2 - 2*w2i*w1i +&
        w1i**2 + w2r**2 - 2*w2r*w1r)/rhogeq 
  endif
  
  print*,'w1 = ',w1r,w1i
  print*,'w2 = ',w2r,w2i
  print*,'w3 = ',w3r,w3i
  print*,' vgas1 = ',vg1r,vg1i
  print*,' vgas2 = ',vg2r,vg2i
  print*,' vgas3 = ',vg3r,vg3i
  print*,' vdust1 = ',vd1r,vd1i
  print*,' vdust2 = ',vd2r,vd2i
  print*,' vdust3 = ',vd3r,vd3i

!-------------------------------
! F I N A L  S O L U T I O N
!-------------------------------
  do i=1,size(xplot)
     xk =  2.*pi/lambda*(xplot(i)-x0)
     arg1 = xk - w1r*time
     arg2 = xk - w2r*time
     arg3 = xk - w3r*time
     vgas = vgeq &
          + vg1r*exp(-w1i*time)*cos(arg1) - vg1i*exp(-w1i*time)*sin(arg1) &
          + vg2r*exp(-w2i*time)*cos(arg2) - vg2i*exp(-w2i*time)*sin(arg2) &
          + vg3r*exp(-w3i*time)*cos(arg3) - vg3i*exp(-w3i*time)*sin(arg3)

     vdust = vdeq &
            + vd1r*exp(-w1i*time)*cos(arg1) - vd1i*exp(-w1i*time)*sin(arg1) &
            + vd2r*exp(-w2i*time)*cos(arg2) - vd2i*exp(-w2i*time)*sin(arg2) &
            + vd3r*exp(-w3i*time)*cos(arg3) - vd3i*exp(-w3i*time)*sin(arg3)

     rhogas = rhogeq &
            + rhog1r*exp(-w1i*time)*cos(arg1) - rhog1i*exp(-w1i*time)*sin(arg1) &
            + rhog2r*exp(-w2i*time)*cos(arg2) - rhog2i*exp(-w2i*time)*sin(arg2) &
            + rhog3r*exp(-w3i*time)*cos(arg3) - rhog3i*exp(-w3i*time)*sin(arg3)

     rhodust = rhodeq &
             + rhod1r*exp(-w1i*time)*cos(arg1) - rhod1i*exp(-w1i*time)*sin(arg1) &
             + rhod2r*exp(-w2i*time)*cos(arg2) - rhod2i*exp(-w2i*time)*sin(arg2) &
             + rhod3r*exp(-w3i*time)*cos(arg3) - rhod3i*exp(-w3i*time)*sin(arg3)

     select case(iplot)
     case(2)
        yplot1(i) = vdust
        yplot(i) = vgas     
     case default
        yplot1(i) = rhodust
        yplot(i) = rhogas
     end select
  enddo
  
  !
  !--plot dust solution, but only if Kdrag > 0
  !
  if (abs(Kdrag).gt.tiny(Kdrag)) then
     call plot_sls(2)
     call plot_line(size(xplot),xplot,yplot1)
     call plot_sls(1)
  endif
    
  return
end subroutine exact_dustywave

end module dustywaves
