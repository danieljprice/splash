!
! subroutine to plot the interaction region for a given particle
! in general coordinate systems (e.g. cylindrical coordinates)
!
! input:  igeom : coordinate system (0=cartesian, 1=cylindrical, 2=spherical)
!         x,y   : particle location in cartesian space
!         h     : size of smoothing sphere 
!                 (assumed isotropic in coordinate space)
!
! PGPLOT page must already be set up - this just draws the "circle"
!
subroutine plot_kernel_gr(igeom,x,y,h)
  implicit none
  integer, parameter :: npts = 100 ! big enough to give a smooth circle
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: igeom
  integer :: i
  real, intent(in) :: x,y,h
  real, dimension(2) :: xcoord, xtemp, xin
  real, dimension(2,npts) :: xpts
  real :: angle, dangle  

  xin(1) = x
  xin(2) = y
!
!--translate the particle's co-ordinate to the appropriate
!  co-ordinate system
!  (NB could do this the other way round ie. input in the co-ordinate space
!   and translate back to cartesian space)
!
  call coord_transform(xin,2,1,xcoord,2,igeom)
!
!--now step around a circle in co-ordinate space of radius h and store the 
!  location of the points in cartesian space in the 2D array xpts
!
  dangle = 2.*pi/REAL(npts-1)
  do i=1,npts
     angle = (i-1)*dangle
     xtemp(1) = xcoord(1) + h*COS(angle) 
     xtemp(2) = xcoord(2) + h*SIN(angle)
!
!--translate back to cartesian space
!
     call coord_transform(xtemp,2,igeom,xpts(:,i),2,1)
  enddo
!
!--now plot the circle using pgline
!
  call pgline(npts,xpts(1,:),xpts(2,:))

  return
end subroutine plot_kernel_gr
