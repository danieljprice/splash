!
!--rotates the particles around the z (and y) axes in 2 (and 3) dimensions
!
!  arguments:
!       angles(ndim-1) : azimuthal and (in 3D) tilt angle in radians
!       xin(ndim)      : input co-ordinates
!       xout(ndim)     : output (rotated) co-ordinates
!       xorigin(ndim)  : input co-ordinates of the origin
!
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
  !--convert to spherical polar co-ordinates
  !
  call coord_transform(xintemp,ndim,1,xtemp,ndim,3)
  !
  !--increase rotation and tilt angles appropriately
  !
  xtemp(2:ndim) = xtemp(2:ndim) - angles(1:2)
  !
  !--now convert back to cartesians
  !
  call coord_transform(xtemp,ndim,3,xout,ndim,1)
  !
  !--adjust positions back from origin
  !
  xout(:) = xout(:) + xorigin(:)

  return
end subroutine rotate
