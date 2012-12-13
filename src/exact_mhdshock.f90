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
! Plot exact solution for a magnetohydrodynamic shock
! (ie. one dimensional MHD Riemann problem)
!
! For want of a better solution these are just taken from the tables in
! Ryu & Jones (1995), ApJ 442, 228 or from ruler and pencil
! on the results in Balsara (1998), or from running the Athena code via
! http://rainman.astro.illinois.edu/ddr/oned/cgi-bin/athena.pl
!
! ----------------------------------------------------------------------
module mhdshock
 implicit none
 public :: exact_mhdshock
 integer, parameter, public :: nmhdshocksolns = 7
 character(len=23), dimension(nmhdshocksolns), parameter, public :: mhdprob = &
    (/'Brio/Wu (gamma=2)      ', &
      'fast/slow shock  (RJ95)', &
      '7 jump shock     (RJ95)', &
      'isothermal shock (B98) ', &
      'rarefaction wave (RJ95)', &
      'Mach 25 shock    (DW94)', &
      'Toth shock   (RJ95/T00)'/)


contains

subroutine exact_mhdshock(iplot,ishk,time,gamma,xmin,xmax,xshock,xpts,ypts,npts,ierr)
  implicit none
  integer, intent(in) :: iplot,ishk
  integer, intent(out) :: npts,ierr
  real, intent(in) :: time,gamma,xmin,xmax,xshock
  real, dimension(:), intent(inout) :: xpts
  real, dimension(size(xpts)), intent(out) :: ypts

  real, dimension(16) :: rho,pr,vx,vy,vz,By,Bz
  real :: const, Bxzero,tfac

  print*,' Plotting exact mhd shock #',ishk,' at t = ',time
  !
  !--set up grid for exact solution
  !
  const = 1./SQRT(4.*3.1415926536)
  ierr = 0

  select case(ishk)   ! which solution to plot
  case(1)
     !
     !--Brio & Wu problem with gamma = 2
     !
     tfac = time/0.1
     vz = 0.
     Bz = 0.
     Bxzero = 0.75
     npts = 14
     xpts(1) = xmin
     xpts(2) = -0.18*tfac
     xpts(3) = -0.08*tfac
     xpts(4:6) = -0.03*tfac
     xpts(7) = -0.005*tfac
     xpts(8:9) = 0.06*tfac
     xpts(10:11) = 0.147*tfac
     xpts(12) = 0.33*tfac
     xpts(13) = 0.36*tfac
     xpts(14) = xmax

     rho(1:2) = 1.0
     rho(3:4) = 0.67623
     rho(5) = 0.827
     rho(6) = 0.775
     rho(7:8) = 0.6962
     rho(9:10) = 0.2352
     rho(11:12) = 0.117
     rho(13:14) = 0.125

     pr(1:2) = 1.0
     pr(3:4) = 0.447
     pr(5) = 0.727219
     pr(6) = 0.6
     pr(7:10) = 0.5160
     pr(11:12) = 0.0876
     pr(13:14) = 0.1

     vx(1:2) = 0.0
     vx(3:4) = 0.63721
     vx(5) = 0.48
     vx(6) = 0.52
     vx(7:10) = 0.600
     vx(11:12) = -0.24
     vx(13:14) = 0.0

     vy(1:2) = 0.0
     vy(3:4) = -0.23345
     vy(5) = -1.3
     vy(6) = -1.4
     vy(7:10) = -1.584
     vy(11:12) = -0.166
     vy(13:14) = 0.

     By(1:2) = 1.0
     By(3:4) = 2.1*const
     By(5) = -1.2*const
     By(6) = -1.3*const
     By(7:10) = -1.9*const
     By(11:12) = -3.25*const
     By(13:14) = -1.0

  case(2)
     !
     !--fast/slow shock from RJ95
     !
     tfac = time/0.15
     vz = 0.
     Bz = 0.
     Bxzero = 1.0
     npts = 12
     xpts(1) = xmin
     xpts(2) = -0.27*tfac
     xpts(3) = -0.09*tfac
     xpts(4) = -0.03*tfac
     xpts(5) = -0.01*tfac
     xpts(6:7) = 0.135*tfac
     xpts(8:9) = 0.25*tfac
     xpts(10:11) = 0.35*tfac
     xpts(12) = xmax

     rho(1:2) = 1.0
     rho(3:4) = 0.5955
     rho(5:6) = 0.55151
     rho(7:8) = 0.41272
     rho(9:10) = 0.2337
     rho(11:12) = 0.2

     pr(1:2) = 1.0
     pr(3:4) = 0.42629
     pr(5:8) = 0.37090
     pr(9:10) = 0.12402
     pr(11:12) = 0.1

     vx(1:2) = 0.0
     vx(3:4) = 0.81237
     vx(5:8) = 0.89416
     vx(9:10) = 0.24722
     vx(11:12) = 0.0

     vy(1:2) = 0.0
     vy(3:4) = -0.59961
     vy(5:8) = -0.5447
     vy(9:10) = -0.91164
     vy(11:12) = 0.

     By(1:2) = 1.0
     By(3:4) = 0.28431
     By(5:8) = 0.31528
     By(9:10) = 0.43086
     By(11:12) = 0.0

  case(3)
     !
     !--problem with 7 discontinuities from RJ95
     !
     tfac = time/0.2
     Bxzero = 2.*const
     npts = 16
     xpts(1) = xmin
     xpts(2:3) = -0.19*tfac
     xpts(4:5) = 0.03*tfac
     xpts(6:7) = 0.051*tfac
     xpts(8:9) = 0.12*tfac    ! contact discontinuity
     xpts(10:11) = 0.18*tfac
     xpts(12:13) = 0.205*tfac
     xpts(14:15) = 0.45*tfac
     xpts(16) = xmax

     rho(1:2) = 1.08
     rho(3:4) = 1.4903
     rho(5:8) = 1.6343
     rho(9:10) = 1.4735
     rho(11:14) = 1.3090
     rho(15:16) = 1.0

     pr(1:2) = 0.95
     pr(3:4) = 1.6558
     pr(5:10) = 1.9317
     pr(11:14) = 1.5844
     pr(15:16) = 1.0

     vx(1:2) = 1.2
     vx(3:5) = 0.60588
     vx(6:11) = 0.57538
     vx(12:14) = 0.53432
     vx(15:16) = 0.0

     vy(1:2) = 0.01
     vy(3:4) = 0.11235
     vy(5:6) = 0.22157
     vy(7:8) = 0.047602
     vy(9:10) = 0.047601
     vy(11:12) = -0.18411
     vy(13:14) = -0.094572
     vy(15:16) = 0.0

     vz(1:2) = 0.5
     vz(3:4) = 0.55686
     vz(5:6) = 0.30125
     vz(7:10) = 0.24734
     vz(11:12) = 0.17554
     vz(13:14) = -0.047286
     vz(15:16) = 0.0

     By(1:2) = 1.0155
     By(3:4) = 1.4383
     By(5:6) = 1.5716
     By(7:10) = 1.4126
     By(11:12) = 1.6103
     By(13:14) = 1.5078
     By(15:16) = 1.1284

     Bz(1:2) = 0.56419
     Bz(3:4) = 0.79907
     Bz(5:6) = 0.48702
     Bz(7:10) = 0.43772
     Bz(11:12) = 0.49899
     Bz(13:14) = 0.75392
     Bz(15:16) = 0.56419

  case(4)
     !
     !--isothermal MHD problem from Balsara (1998)
     !
     tfac = time/0.2
     Bxzero = 2.*const
     npts = 14
     xpts(1) = xmin
     xpts(2:3) = -0.15*tfac
     xpts(4:5) = 0.035*tfac
     xpts(6:7) = 0.07*tfac
     xpts(8:9) = 0.17*tfac
     xpts(10:11) = 0.2*tfac
     xpts(12:13) = 0.41*tfac
     xpts(14) = xmax

     rho(1:2) = 1.08
     rho(3:6) = 1.515
     rho(7:8) = 1.745
     rho(9:12) = 1.36
     rho(13:14) = 1.0

     vx(1:2) = 1.2
     vx(3:6) = 0.65
     vx(7:8) = 0.62
     vx(9:12) = 0.54
     vx(13:14) = 0.0

     vy(1:2) = 0.01
     vy(3:4) = 0.13
     vy(5:6) = 0.24
     vy(7:8) = 0.071
     vy(9:10) = -0.215
     vy(11:12) = -0.125
     vy(13:14) = 0.0

     vz(1:2) = 0.5
     vz(3:4) = 0.57
     vz(5:6) = 0.31
     vz(7:8) = 0.255
     vz(9:10) = 0.165
     vz(11:12) = -0.06
     vz(13:14) = 0.0

     By(1:2) = 3.6*const
     By(3:4) = 5.2*const
     By(5:6) = 5.7*const
     By(7:8) = 5.22*const
     By(9:10) = 5.96*const
     By(11:12) = 5.58*const
     By(13:14) = 4.0*const

     Bz(1:2) = 2.0*const
     Bz(3:4) = 2.885*const
     Bz(5:6) = 1.76*const
     Bz(7:8) = 1.62*const
     Bz(9:10) = 1.85*const
     Bz(11:12) = 2.79*const
     Bz(13:14) = 2.0*const

     pr = rho

  case(5)
     !
     !--rarefaction from RJ95
     !
     tfac = time/0.1
     npts = 6
     vy = 0.
     vz = 0.
     Bz = 0.
     Bxzero = 0.
     xpts(1) = xmin
     xpts(2) = -0.27*tfac
     xpts(3) = -0.12*tfac
     xpts(4) = 0.12*tfac
     xpts(5) = 0.27*tfac
     xpts(6) = xmax

     rho(1:2) = 1.0
     rho(3:4) = 0.49653
     rho(5:6) = 1.0

     pr(1:2) = 1.0
     pr(3:4) = 0.31134
     pr(5:6) = 1.0

     vx(1:2) = -1.0
     vx(3:4) = 0.    ! this is approximate (to 10-7)
     vx(5:6) = 1.0

     By(1:2) = 1.0
     By(3:4) = 0.49638
     By(5:6) = 1.0

  case(6)
     !
     !--mach 25 shocks from Dai and Woodward (1994)
     !
     tfac = time/0.03
     Bxzero = 4.*const
     npts = 6
     rho(1:2) = 1.0
     rho(3:4) = 3.982
     rho(5:6) = 0.1

     pr(1:2) = 1.0
     pr(3:4) = 1806.0
     pr(5:6) = 1.0

     !       machno = 0.5*25.5
     !       vs = SQRT(gamma*pr(1)/rho(1))
     !
     xpts(1) = xmin
     xpts(2:3) = -0.35*tfac
     xpts(4:5) = 0.35*tfac
     xpts(6) = xmax

     vx(1:2) = 36.87
     vx(3:4) = 0.0
     vx(5:6) = -36.87

     vy(1:2) = -0.1546
     vy(3:4) = -0.07727
     vy(5:6) = 0.0

     vz(1:2) = -0.03864
     vz(3:4) = -0.01932
     vz(5:6) = 0.0

     By(1:2) = 4.0*const
     By(3:4) = 15.95*const
     By(5:6) = 4.0*const

     Bz(1:2) = 1.0*const
     Bz(3:4) = 3.988*const
     Bz(5:6) = 1.0*const

  case(7)
     !
     !--Problem 1A in Ryu and Jones (1995)
     !
     Bxzero = 5.*const
     npts = 12
     tfac = time/0.08
     xpts(1) = xmin
     xpts(2:3) = -0.386*tfac
     xpts(4:5) = -0.01*tfac
     xpts(6:7) = 0.0505*tfac
     xpts(8:9) = 0.12*tfac
     xpts(10:11) = 0.37*tfac
     xpts(12) = xmax

     rho(1:2) = 1.0
     rho(3:4) = 2.6797
     rho(5:6) = 2.6713
     rho(7:8) = 3.8508
     rho(9:10) = 3.7481
     rho(11:12) = 1.0

     pr(1:2) = 20.0
     pr(3:4) = 150.98
     pr(5:8) = 150.19
     pr(9:10) = 143.57
     pr(11:12) = 1.0

     vx(1:2) = 10.0
     vx(3:4) = 0.72113
     vx(5:8) = 0.72376
     vx(9:10) = 0.70505
     vx(11:12) = -10.0

     vy(1:2) = 0.0
     vy(3:4) = 0.23139
     vy(5:8) = 0.35684
     vy(9:10) = -0.38804
     vy(11:12) = 0.0

     vz(1:12) = 0.0

     By(1:2) = 1.4105
     By(3:4) = 3.8389
     By(5:8) = 4.0380
     By(9:10) = 5.4272
     By(11:12) = 1.4105

     Bz(1:12) = 0.0

  case default
     ierr = 1
     npts = 0
     ypts = 0.
     xpts = 0.
     return

  end select

  !
  !--plot just the initial conditions at t=0
  !
  if (abs(time).le.0.) then
     rho(1:2) = rho(1)
     pr(1:2)  = pr(1)
     vx(1:2)  = vx(1)
     vy(1:2)  = vy(1)
     vz(1:2)  = vz(1)
     By(1:2)  = By(1)
     Bz(1:2)  = Bz(1)
     xpts(3)  = 0.
     xpts(4)  = xpts(npts)
     rho(3:4) = rho(npts)
     pr(3:4)  = pr(npts)
     vx(3:4)  = vx(npts)
     vy(3:4)  = vy(npts)
     vz(3:4)  = vz(npts)
     By(3:4)  = By(npts)
     Bz(3:4)  = Bz(npts)
     npts   = 4
  endif

  !
  !--translate positions if shock is not at x=0
  !
  xpts(:) = xpts(:) + xshock

  !
  !--determine which parameter to plot
  !
  select case(iplot)
  case(1)
     ypts(1:npts) = rho(1:npts)
  case(2)
     ypts(1:npts) = pr(1:npts)
  case(3)
     ypts(1:npts) = vx(1:npts)
  case(4)
     ypts(1:npts) = vy(1:npts)
  case(5)
     ypts(1:npts) = vz(1:npts)
  case(6)
     ypts(1:npts) = By(1:npts)
  case(7)
     ypts(1:npts) = Bz(1:npts)
  case(8)
     print*,'gamma = ',gamma
     if (abs(gamma-1.).gt.1.e-5) then
        where (abs(rho(1:npts)) > 0.)
           ypts(1:npts) = pr(1:npts) / ((gamma-1.)*rho(1:npts))
        end where
     else
        print*,' ***isothermal: utherm solution not valid'
        ypts(1:npts) = 0.
     endif
  case(9)
     ypts(1:npts) = Bxzero
  case default
     print*,'error: unknown solution to plot'
  end select

  return
end subroutine exact_mhdshock

end module mhdshock
