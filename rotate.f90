!
!--rotates the particles around the z (and y) axes in 2 (and 3) dimensions
!
!  arguments:
!       angles(ndim-1) : azimuthal and (in 3D) tilt angle
!       xin(ndim)      : input co-ordinates
!       xout(ndim)     : output (rotated) co-ordinates
!
subroutine rotate(angles,xin,xout,ndim)
  implicit none
  integer, intent(in) :: ndim
  real, intent(in), dimension(ndim-1) :: angles
  real, intent(in), dimension(ndim) :: xin
  real, intent(out), dimension(ndim) :: xout
  real, dimension(ndim) :: xtemp

  if (ndim.lt.2 .or. ndim.gt.3) then
     print*,'error: rotate: ndim < 2 or > 3'
     return
  endif
  !
  !--work out radius, azimuthal and tilt angle
  !  (ie. convert to spherical co-ordinates)
  !
  call coord_transform(xin,ndim,1,xtemp,ndim,3)

  !--increase angles appropriately
  xtemp(2:ndim) = xtemp(2:ndim) + angles

  !--now convert back to cartesians
  call coord_transform(xtemp,ndim,3,xout,ndim,1)

  return
end subroutine rotate
