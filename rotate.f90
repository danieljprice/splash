!
! This module contains all the routines for rotating the particles
! and plotting the rotated axes
!
module rotation
 implicit none
!
!--2D rotation (about z axis)
!
contains

subroutine rotate2D(xcoords,anglez)
  implicit none
  real, intent(inout) :: xcoords(2)
  real, intent(in) :: anglez
  real :: x, y, r, phi
  
  x = xcoords(1)
  y = xcoords(2)
!
!--rotate about z
!  
  r = sqrt(x**2 + y**2)
  phi = ATAN2(y,x)
  phi = phi - anglez
  xcoords(1) = r*COS(phi)
  xcoords(2) = r*SIN(phi)
  
  return 
end subroutine rotate2D

!
!--3D rotation (about x, y and z axes)
!
subroutine rotate3D(xcoords,anglex,angley,anglez)
  implicit none
  real, intent(inout) :: xcoords(3)
  real, intent(in) :: anglex, angley, anglez
  real :: x, y, z, r, phi
  
  x = xcoords(1)
  y = xcoords(2)
  z = xcoords(3)
!
!--rotate about z
!  
  r = sqrt(x**2 + y**2)
  phi = ATAN2(y,x)
  phi = phi - anglez
  x = r*COS(phi)
  y = r*SIN(phi)
!
!--rotate about y
!
  r = sqrt(x**2 + z**2)
  phi = ATAN2(x,z)
  phi = phi - angley  
  z = r*COS(phi)
  x = r*SIN(phi)
!
!--rotate about x
!
  r = sqrt(y**2 + z**2)
  phi = ATAN2(y,z)
  phi = phi - anglex  
  z = r*COS(phi)
  y = r*SIN(phi)
  
  xcoords(1) = x
  xcoords(2) = y
  xcoords(3) = z
  
  return
end subroutine rotate3D

!
!--plots rotated plot axes
!
subroutine rotate_axes2D(ioption,xmin,xmax,anglez)
  implicit none
  integer, intent(in) :: ioption
  real, intent(in), dimension(2) :: xmin,xmax
  real, intent(in) :: anglez
  integer :: i
  real, dimension(2,4) :: xpt
  print*,'plotting rotated axes...'

!--front face (pts 1->4)
  xpt(:,1) = xmin(:)  ! xmin, ymin

  xpt(1,2) = xmin(1)  ! xmin
  xpt(2,2) = xmax(2)  ! ymax

  xpt(1,3) = xmax(1)  ! xmax
  xpt(2,3) = xmax(2)  ! ymax

  xpt(1,4) = xmax(1)  ! xmax
  xpt(2,4) = xmin(2)  ! ymin
!
!--now rotate each of these coordinates
!
  do i=1,4
     call rotate2D(xpt(:,i),anglez)
  enddo
!
!--now plot boxes appropriately using points
!
  call pgsfs(2)
  call pgpoly(4,xpt(1,1:4),xpt(2,1:4))
  
  return
end subroutine rotate_axes2D

subroutine rotate_axes3D(ioption,xmin,xmax,anglez,angley,anglex)
  implicit none
  integer, intent(in) :: ioption
  real, intent(in), dimension(3) :: xmin,xmax
  real, intent(in) :: anglez, angley, anglex
  integer :: i
  real, dimension(3,8) :: xpt
  print*,'plotting rotated axes...'

!--front face (pts 1->4)
  xpt(:,1) = xmin(:)  ! xmin, ymin

  xpt(1,2) = xmin(1)  ! xmin
  xpt(2,2) = xmax(2)  ! ymax

  xpt(1,3) = xmax(1)  ! xmax
  xpt(2,3) = xmax(2)  ! ymax

  xpt(1,4) = xmax(1)  ! xmax
  xpt(2,4) = xmin(2)  ! ymin
!
!--now rotate each of these coordinates
!
  do i=1,4
     call rotate2D(xpt(:,i),anglez)
  enddo
!
!--now plot boxes appropriately using points
!
  call pgsfs(2)
  call pgpoly(4,xpt(1,1:4),xpt(2,1:4))
  
  return
end subroutine rotate_axes3D

end module rotation
