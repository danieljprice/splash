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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!------------------------------------------------------------
! plot exact solution for toystar in one dimension
!
! linear solution uses Gegenbauer/Legendre polynomials
!
! non-linear solution solves ODEs, assumes linear velocity
!
! For details see Monaghan and Price (2004) MNRAS
!------------------------------------------------------------
module toystar1D
 implicit none
 public :: exact_toystar1D, exact_toystar_ACplane
 private :: Pm, Gn

contains

subroutine exact_toystar1D(iplot,time,gamma,H0,A0,C0, &
                           sigma,norder,xplot,yplot,npts,ierr)
  use plotlib, only:plot_pt1
  implicit none
  integer, intent(in) :: iplot,norder,npts
  integer, intent(out) :: ierr
  real, intent(in) :: time,gamma,sigma
  real, intent(in) :: H0, C0, A0    ! parameters for toy star
  real, dimension(:), intent(inout) :: xplot
  real, dimension(size(xplot)), intent(out) :: yplot

  integer :: i,nsteps
  real :: const ! parameter for toy star
  real :: Hprev, Cprev, Aprev, Htemp, Ctemp, Atemp, H, C, A
  real :: radstar,dx,dt,tnow
  real :: rhoplot,deltarho
  real :: fprevC,fprevA,fprevH,ftempC,ftempA,ftempH
  real :: gamp1,gamm1,gam1,fact,constK,omega
  real :: fnorm
  logical :: linear

  linear = .false.
  if (norder.ge.0) linear = .true.
  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam1 = 1./gamm1
  constK = 0.25

  ierr = 0

  if (linear) then
!---------------------------------------------------------------------------
!  linear solution uses Gegenbauer & Legendre Polynomials
!  (this is for the toy star oscillations)

     omega = sqrt(0.5*(norder+1.)*(norder+2.))
     fnorm = 2.*(norder+1)*(norder+2)/real(2.*norder + 3.)

     print*,' Plotting toy star oscills: time, norder, omega = ', &
          time,norder,omega,H0,C0,A0

     if (C0.le.0.) then
        print*,'*** C = 0 = illegal in input'
        ierr = 1
        return
     else
        radstar = sqrt(H0/C0)
     endif
     xplot(1) = -radstar
     dx = (radstar-xplot(1))/float(npts-1)

     do i=2,npts
        xplot(i) = xplot(1)+dx*(i-1)
        !         print*,i,' x,y = ',xplot(i),yplot(i)
        rhoplot = (H0 - C0*xplot(i)**2)
        if (rhoplot.le.0.) rhoplot = 0.
        deltarho = Pm(xplot(i),norder+1)*sin(omega*time)
        rhoplot = rhoplot**gam1 + 2.*omega*A0*deltarho/fnorm
        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for vx
           yplot(i) = A0*Gn(xplot(i),norder)*cos(omega*time)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     if (iplot.eq.6) then        ! plot By \propto rho
        dx = (H0**gam1)/float(npts-1) ! ie (rhomax - 0)/npts
        xplot(1) = 0.
        yplot(1) = sigma*xplot(1)
        do i=2,npts
           xplot(i) = xplot(1) + dx*(i-1)
           yplot(i) = sigma*(xplot(i))
        enddo
     endif

     if (iplot.eq.7) then        ! plot current point on A-C plane
        call plot_pt1(C0,A0*cos(omega*time),4)
        ierr = 2 ! do not plot again
     else                        ! plot normal exact solution line
        ierr = 0
     endif

!---------------------------------------------------------------------------
!  non-linear solution for the fundamental (n=1) mode
!
  else
     !
     !  solve for H, C and A given initial conditions on v, rho and the time.
     !
     nsteps = 1000*(int(time) + 1)

     Hprev = H0
     Cprev = C0
     Aprev = A0

     !      PRINT*,' nsteps,H,C,A in = ',nsteps,H0,C0,A0
     dt = time/nsteps

     fact = 2.*(constK + 0.5*sigma**2)*gamma*gam1

     const = (A0**2 + 1. + 2.*fact*C0*gam1)*C0**(-2./gamp1)

     print*,' Plotting toy star: time, H0, C0, A0, k = ', &
          time,Hprev,Cprev,Aprev,const

     tnow = 0.
     do i = 1,nsteps
        tnow = tnow + dt
        ! integrate using improved Euler
        fprevC = -Cprev*Aprev*gamp1
        fprevA = fact*Cprev -1.-Aprev**2
        fprevH = -Aprev*Hprev*gamm1
        ! predictor
        Ctemp = Cprev + dt*(-Cprev*Aprev*gamp1)
        Atemp = Aprev + dt*(fact*Cprev-1.-Aprev**2)
        Htemp = Hprev + dt*(-Aprev*Hprev*gamm1)

        ftempC = -Ctemp*Atemp*gamp1
        ftempA = fact*Ctemp -1. -Aprev**2
        ftempH = -Atemp*Htemp*gamm1
        ! corrector
        C = Cprev + 0.5*dt*(fprevC + ftempC)
        A = Aprev + 0.5*dt*(fprevA + ftempA)
        H = Hprev + 0.5*dt*(fprevH + ftempH)

        Cprev = C
        Aprev = A
        Hprev = H

        !         print*,' time = ',tnow
        !         IF ((abs(C-C0).LT.5.e-3).AND. &
        !                  (abs(A-A0).LT.5.e-3).AND.(tnow.GT.5e-3)) THEN
        !             PRINT*,'*** period, t = ',tnow,' err = ',abs(C-C0)+abs(A-A0)
        !         ENDIF

     enddo

     const = (A**2 + 1. + 2.*fact*C*gam1)*C**(-2./gamp1)

     print*,' C, A, H, k = ',C,A,H,const

     if (C.le.0.) then
        radstar = 0.5
        print*,'*** C = 0 = illegal'
        ierr = 1
        return
     else
        radstar = sqrt(H/C)
     endif
     xplot(1) = -radstar
     dx = (radstar-xplot(1))/float(npts-1)

     do i=2,npts
        xplot(i) = xplot(1)+dx*(i-1)
        !         print*,i,' x,y = ',xplot(i),yplot(i)
        rhoplot = (H - C*xplot(i)**2)
        if (rhoplot.le.0.) rhoplot = 0.
        rhoplot = rhoplot**gam1
        select case(iplot)
        case(1)                 ! plot solution for density
           yplot(i) = rhoplot
        case(2)                 ! plot solution for pressure
           yplot(i) = constK*rhoplot**gamma
        case(3)                 ! plot solution for utherm
           yplot(i) = constK*(rhoplot**gamm1)/gamm1
        case(4)                 ! plot solution for vx
           yplot(i) = A*xplot(i)
        case(5)                 ! plot solution for By
           yplot(i) = sigma*rhoplot
        end select

     enddo

     if (iplot.eq.6) then        ! plot By \propto rho
        dx = (H**gam1)/float(npts-1) ! ie (rhomax - 0)/npts
        xplot(1) = 0.
        yplot(1) = sigma*xplot(1)
        do i=2,npts
           xplot(i) = xplot(1) + dx*(i-1)
           yplot(i) = sigma*(xplot(i))
        enddo
     endif

     if (iplot.eq.7) then        ! plot current point on A-C plane
        call plot_pt1(C,A,4)
        ierr = 2 ! do not plot again
     else                        ! plot normal exact solution line
        ierr = 0
     endif
!
!------------------------------------------------------------------------
!
  endif

  return
end subroutine exact_toystar1D
!
!--function to evaluate the Gegenbauer polynomial of index n given x
!
real function Gn(x,n)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x
  integer :: i
  real :: Gnminus1,Gnminus2
  real :: fnorm

  fnorm = 2.*(n+1)*(n+2)/real(2.*n + 3.)
  !      PRINT*,' fnorm = ',fnorm
  !
  !--specify first two Gegenbauer polynomials
  !
  Gnminus2 = 1.
  Gnminus1 = 3.*x
  Gn = 0. ! avoid compiler warning
  !
  !--use recurrence relation to calculate the rest
  !
  select case (n)
  case (0)
     Gn = Gnminus2
  case (1)
     Gn = Gnminus1
  case (2:)
     do i=2,n
        Gn = ((2*i+1)*x*Gnminus1 - (i+1)*Gnminus2)/real(i)
        Gnminus2 = Gnminus1
        Gnminus1 = Gn
     enddo
  end select

  Gn = Gn/fnorm

end function Gn

!
!--function to calculate a Legendre Polynomial of order m
!
real function Pm(x,m)
  implicit none
  integer, intent(in) :: m
  real, intent(in) :: x
  integer :: i
  real :: Pmminus1,Pmminus2
  !
  !--specify first two Legendre polynomials
  !
  Pmminus2 = 1.
  Pmminus1 = x
  Pm = 0. ! avoid compiler warning

  select case(m)
  case (0)
     Pm = 1.
  case (1)
     Pm = x
  case (2:)        ! use recurrence relation to calculate the rest
     do i=2,m
        Pm = ((2.*(i-1.)+1.)*x*Pmminus1 - (i-1.)*Pmminus2)/real(i)
        Pmminus2 = Pmminus1
        Pmminus1 = Pm
     enddo
  end select

end function Pm


!----------------------------------------------------------------------
!
! this subroutine plots the A-C relation in the 1D Toy star solution
!
!----------------------------------------------------------------------
subroutine exact_toystar_ACplane(astart,cstart,sigma,gamma)
  use plotlib, only:plot_swin,plot_funx,plot_label,plot_box
  implicit none
  real, intent(in) :: astart,cstart,sigma,gamma
  real :: constk,gam1,gamm1,gamp1,fact
  real :: xstart,xend,xcentre,c,cnew,k
  real :: func,func2,funct,fderiv,ymin,ymax,extra
  external func,func2
  common /kconst/ k,fact,gam1,gamp1

  print*,' plotting a-c plane...'

  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam1 = 1./gamm1
  constk = 0.25
  !      print*,' k, kdash = ',constk,constk + 0.5*sigma**2
  fact = 2.*(constk + 0.5*sigma**2)*gamma*gam1

  k = (astart**2 + 1. + 2.*fact*cstart*gam1)*cstart**(-2./gamp1)

  !      print*,' k,fact = ',k,fact
  !
  !--find limits of plot (ie. where a = 0)
  !
  c = 1.e6
  cnew = 0.25

  do while (abs(c-cnew).gt.1.e-5)

     c = cnew

     funct = k*c**(2./gamp1) - 2.*fact*c*gam1 - 1.
     fderiv = 2.*k/gamp1*c**(-gamm1/gamp1) - 2.*fact*gam1

     cnew = c - funct/fderiv

     if (cnew.lt.0.) print*,'eek c < 0'

  enddo

  xstart = cnew

  c = 1.e6
  cnew = 6.37935

  do while (abs(c-cnew).gt.1.e-5)

     c = cnew

     funct = k*c**(2./gamp1) - 2.*fact*c*gam1 - 1.
     fderiv = 2.*k/gamp1*c**(-gamm1/gamp1) - 2.*fact*gam1

     cnew = c - funct/fderiv

  enddo

  xend = cnew

  !      print*,'plotting k = ',k,' cstart = ',cstart,' astart = ',astart
  !      print*,'min c = ',xstart,' max c = ',xend

  xstart = xstart + 0.000001
  xend = xend - 0.000001

  extra = 0.1*(xend-xstart)
  xcentre = 0.5*(xstart + xend)
  ymax = 1.5*func(xcentre)
  ymin = 1.5*func2(xcentre)

  call plot_swin(xstart-extra,xend+extra,ymin,ymax)
  call plot_box('bcnst',0.0,0,'1bvcnst',0.0,0)
  call plot_funx(func,10000,xstart,xend,1)
  call plot_funx(func2,10000,xstart,xend,1)

  call plot_label ('c','a',' ')
  return

end subroutine exact_toystar_ACplane

end module toystar1D

!------------------------------------
!
! these functions must be external
!
!------------------------------------
real function func(x)
  real, intent(in) :: x
  real :: k,term,fact,gam1,gamp1
  common /kconst/ k,fact,gam1,gamp1

  !      print*,'k = ',k

  term = -1 -2.*fact*x*gam1 + k*x**(2./gamp1)
  if (term.le.0.) then
  !         print*,' warning: func < 0 ',term
     func = 0.
  else
     func = sqrt(term)
  endif

end function func

real function func2(x)
  implicit none
  real, intent(in) :: x
  real :: k,term,fact,gam1,gamp1
  common /kconst/ k,fact,gam1,gamp1

  !      print*,' k = ',k

  term = -1 -2.*fact*x*gam1 + k*x**(2./gamp1)
  if (term.le.0.) then
     func2 = 0.
  else
     func2 = -sqrt(term)
  endif

end function func2
