!
!--rotates the particles around the z (and y) axes in 2 (and 3) dimensions
!
!  arguments:
!       angles(ndim-1) : azimuthal and (in 3D) tilt angle in radians
!       xin(ndim)      : input co-ordinates
!       xout(ndim)     : output (rotated) co-ordinates
!       xorigin(ndim)  : input co-ordinates of the origin
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

subroutine rotate(angles,xin,xout,xorigin,ndim)
  implicit none
  integer, intent(in) :: ndim
  real, intent(in), dimension(ndim-1) :: angles
  real, intent(in), dimension(ndim) :: xin,xorigin
  real, intent(out), dimension(ndim) :: xout
  real, dimension(ndim) :: xintemp,xtemp,xouttemp

  if (ndim.lt.2 .or. ndim.gt.3) then
     print*,'error: rotate: ndim < 2 or > 3'
     return
  endif
  !
  !--adjust x-y positions according to location of rotation axis
  !
  xintemp(:) = xin(:) - xorigin(:)
  !
  !--convert to cylindrical polar co-ordinates
  !
  call coord_transform(xintemp,ndim,1,xtemp,ndim,2)
  !
  !--increase rotation and tilt angles appropriately
  !
  xtemp(2) = xtemp(2) - angles(1)
  !
  !--now convert back to cartesians
  !
  call coord_transform(xtemp,ndim,2,xout,ndim,1)
  !
  !--adjust positions back from origin
  !
  xout(:) = xout(:) + xorigin(:)

  return
end subroutine rotate
