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

subroutine rotate2D(xcoords,anglez,angley)
  implicit none
  real, intent(inout) :: xcoords(2)
  real, intent(in) :: anglez,angley
  real :: x, y, r, phi
  
  x = xcoords(1)
  y = xcoords(2)
!
!--rotate about z
! 
  r = sqrt(x**2 + y**2)
  phi = ATAN2(y,x)
  phi = phi - anglez
  x = r*COS(phi)
  y = r*SIN(phi)

  xcoords(1) = x
  xcoords(2) = y
  
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
  phi = ATAN2(z,y)
  phi = phi - anglex  
  y = r*COS(phi)
  z = r*SIN(phi)
  
  xcoords(1) = x !! *tan(angley)
  xcoords(2) = y !!*tan(anglex)
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
  print*,'plotting rotated (2D) axes...'

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
     call rotate2D(xpt(:,i),anglez,0.0)
  enddo
!
!--now plot box appropriately using points
!
  call pgsfs(2)
  call pgpoly(4,xpt(1,1:4),xpt(2,1:4))
  
  return
end subroutine rotate_axes2D

subroutine rotate_axes3D(ioption,iplotx,iploty,xmin,xmax,anglex,angley,anglez)
  implicit none
  integer, intent(in) :: ioption,iplotx,iploty
  real, intent(in), dimension(3) :: xmin,xmax
  real, intent(in) :: anglex, angley, anglez
  integer :: i,idim,iline
  integer, parameter :: nlines = 10
  real, dimension(3,8) :: xpt
  real, dimension(2) :: xline,yline
  real :: dx

  !
  ! plot various options for the 3D axes
  !
  select case(ioption)
  case(1)
  !--rotated axes
     print*,'plotting rotated 3D axes...'
     xpt = 0.
     !--origin
     xpt(1:3,1) = 0.
     
     do idim=1,3
        !--plot to max of each axis
        xpt(:,2) = 0.
        xpt(idim,2) = xmax(idim)
        do i=1,2
           call rotate3D(xpt(:,i),anglex,angley,anglez)
        enddo
        !--plot each axis as an arrow
        call pgarro(xpt(iplotx,1),xpt(iploty,1),xpt(iplotx,2),xpt(iploty,2))
!!        call pgline(2,xpt(iplotx,1:2),xpt(iploty,1:2))
     enddo
  case(2)
  !--rotated box
     print*,'plotting rotated 3D box...',iplotx,iploty
     
     !--front face (pts 1->4)
     xpt(:,1) = xmin(:)  ! xmin, ymin

     xpt(1,2) = xmin(1)  ! xmin
     xpt(2,2) = xmax(2)  ! ymax

     xpt(1,3) = xmax(1)  ! xmax
     xpt(2,3) = xmax(2)  ! ymax

     xpt(1,4) = xmax(1)  ! xmax
     xpt(2,4) = xmin(2)  ! ymin
     
     xpt(3,1:4) = xmin(3) ! zmin
   
     !--back face (pts 5->8)
     do i=1,4
        xpt(1:2,i+4) = xpt(1:2,i)
     enddo
     xpt(3,5:8) = xmax(3)  
     !
     !--now rotate each of these coordinates
     !
     do i=1,8
        call rotate3D(xpt(:,i),anglex,angley,anglez)
     enddo
     !
     !--now draw lines appropriately through points
     !
     call pgsfs(2)    
     !--front face
     call pgpoly(4,xpt(iplotx,1:4),xpt(iploty,1:4))
     !--back face
     call pgpoly(4,xpt(iplotx,5:8),xpt(iploty,5:8))
     !--connecting lines ( 1->5, 2->6, 3->7, 4->8 )
     do i=1,4
        xline(1) = xpt(iplotx,i)
        yline(1) = xpt(iploty,i)
        xline(2) = xpt(iplotx,i+4)
        yline(2) = xpt(iploty,i+4)
        print*,xline,yline
        call pgline(2,xline,yline)
     enddo

  case(3)
  !--gridded x-y plane
     print*,'plotting rotated x-y plane...'

     !--lines of constant x

     dx = (xmax(1) - xmin(1))/real(nlines-1)
     do iline=1,nlines
        !--all pts at z = 0
        xpt(3,:) = 0.

        !--start from xmin, plot line from ymin to ymax
        xpt(1,1:2) = xmin(1) + (iline-1)*dx
        xpt(2,1) = xmin(2)
        xpt(2,2) = xmax(2)
        do i=1,2
           call rotate3D(xpt(:,i),anglex,angley,anglez)
        enddo
        call pgline(2,xpt(iplotx,1:2),xpt(iploty,1:2))
     enddo
     
     !--lines of constant y
     dx = (xmax(2) - xmin(2))/real(nlines-1)
     do iline=1,nlines
        !--all pts at z = 0
        xpt(3,:) = 0.

        !--start from ymin, plot line from xmin to xmax
        xpt(2,1:2) = xmin(2) + (iline-1)*dx
        xpt(1,1) = xmin(1)
        xpt(1,2) = xmax(1)
        do i=1,2
           call rotate3D(xpt(:,i),anglex,angley,anglez)
        enddo
        call pgline(2,xpt(iplotx,1:2),xpt(iploty,1:2))
     enddo     
  case default
  !--do nothing
  end select
  
  return
end subroutine rotate_axes3D

end module rotation
