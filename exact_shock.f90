! ----------------------------------------------------------------------
! compute exact solution for the one dimensional Riemann problem
! (hydrodynamic shock)
!
! input parameters are initial left and right states of 
! density, pressure and velocity
!
! Computes shock profile at time t
!
! Calls a separate subroutine to calculate the post-shock pressure
! and velocity (this is the difficult bit).
!
! Daniel Price, Institute of Astronomy, Cambridge, 2004
! dprice@ast.cam.ac.uk
!-----------------------------------------------------------------------

subroutine exact_shock(iplot,time,gamma,rho_L,rho_R,p_L,p_R,v_L,v_R,xmin,xmax)
  implicit none
  integer, parameter :: npts=2000
  integer, intent(in) :: iplot
  real, intent(in) :: time,gamma,xmin,xmax
  real, intent(in) :: rho_L,rho_R,p_L,p_R,v_L,v_R
  
  integer :: i,j
  real, dimension(npts) :: xplot, yplot, dens, pr, vel
  real :: dx,cs_L,cs_R, gamfac
  real :: ppost, vpost, vfan, vshock
  real :: xzero,xleft,xfan,xcontact,xshock

  print*,'Plotting exact Riemann solution at t = ',time,' gamma = ',gamma
!
! check for errors in input
!
  if (rho_L.le.0. .or. rho_R.le.0.) then
     print*,'error: rho <= 0 on input : ',rho_L,rho_R
     return
  elseif (p_L .le.0. .or. p_R .le.0.) then
     print*,'error: pr <= 0 on input ',p_L, p_R
     return
  elseif (gamma.lt.1.0001) then
     print*,'error: isothermal solver not implemented'
     return
  endif
!
! set up grid for exact solution
!          
  dx = (xmax-xmin)/real(npts)
  do i=1,npts
     xplot(i) = xmin + (i-1)*dx
  enddo
!
!  xzero is the position of the shock at t=0
!
  xzero = 0.
!
!  define sound speeds to left and right of shock tube
!      
  cs_L = sqrt(gamma*p_L/rho_L)
  cs_R = sqrt(gamma*p_R/rho_R)
  gamfac = (gamma-1.)/(gamma + 1.)

!------------------------------------------------------------
! find post-shock pressure via the Riemann solver
! (this version also returns the post shock speed vpost
!  although this can be calculated from ppost)
!------------------------------------------------------------

  call riemannsolver(gamma,p_L,p_R,v_L,v_R,cs_L,cs_R,ppost,vpost)

!------------------------------------------------------------
!  using this, calculate various speeds needed in order to 
!  reconstruct the shock profile
!------------------------------------------------------------

  !  post shock velocity (speed at which fluid behind shock is 
  !  moving to the right)
  !
!!  vpost = 
  !
  !  speed of the shock front
  !
  vshock = v_R + cs_R**2*(ppost/p_R - 1.)/(gamma*(vpost-v_R))
  !
  !  speed at which the right end of rarefaction fan moves
  !
  vfan = cs_L - 0.5*(gamma+1.)*vpost + 0.5*(gamma-1.)*v_L

!-------------------------------------------------------------
! now work out the locations of various features in the shock
!-------------------------------------------------------------
!
! left end of expansion fan (propagates at the sound speed into the "left" fluid)
!
  xleft = xzero - (cs_L - v_L)*time 
!
! right end of expansion fan
!
  xfan = xzero - vfan*time
!
! position of the interface between the "left" fluid from the "right" fluid
! (the contact discontinuity)
!
  xcontact = xzero + vpost*time
!
! shock position (distance the shock has travelled into the "right" fluid)
!
  xshock = xzero + vshock*time

!--------------------------------------------------------------
! reconstruct the shock profile for all x
!--------------------------------------------------------------
  
  where(xplot <= xleft)  ! <= otherwise problems at t=0
!    undisturbed medium to the left
     pr = p_L
     dens = rho_L
     vel = v_L
  elsewhere(xplot < xfan)
!    inside expansion fan
     dens = rho_L*(gamfac*(xzero-xplot)/(cs_L*time) + (1.-gamfac))**(2./(gamma-1.))
     pr = p_L*(dens/rho_L)**gamma
     vel = (1.-gamfac)*(cs_L -(xzero-xplot)/time)
  elsewhere(xplot < xcontact)
!    between expansion fan and contact discontinuity
     pr = ppost
     dens = rho_L*(ppost/p_L)**(1./gamma)
     vel = vpost
  elsewhere(xplot < xshock)
!    post-shock, ahead of contact discontinuity
     pr = ppost
     dens = rho_R*(gamfac+ppost/p_R)/(1+gamfac*ppost/p_R)
     vel = vpost
  elsewhere
!    undisturbed medium to the right
     pr = p_R
     dens = rho_R
     vel = v_R
  end where

!------------------------------------
!  determine which solution to plot
!------------------------------------
  select case(iplot)
  case(1)
     yplot = dens
  case(2)
     yplot = pr
  case(3)
     yplot = vel
  case(4)
     yplot = pr/((gamma-1.)*dens)
  end select
!----------------------------------------------------------
!  plot this as a line on the current graph using PGPLOT
!----------------------------------------------------------
  call pgline(npts,xplot,yplot)

  return
end subroutine exact_shock
