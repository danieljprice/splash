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
  real, dimension(ndim) :: xintemp,xtemp

  if (ndim.lt.2 .or. ndim.gt.3) then
     print*,'error: rotate: ndim < 2 or > 3'
     return
  endif
  !
  !--rotation about z axis
  !
  !
  !  adjust x-y positions according to location of rotation axis
  !
  xintemp(1:2) = xin(1:2) - xorigin(1:2)
  if (ndim.eq.3) xintemp(ndim) = xin(ndim)
  !
  !--convert to cylindrical co-ordinates
  !
  call coord_transform(xintemp,ndim,1,xtemp,ndim,2)

  !--increase rotation angle appropriately
  xtemp(2) = xtemp(2) + angles(1)

  !--now convert back to cartesians
  call coord_transform(xtemp,ndim,2,xout,ndim,1)
!  print*,'----- xin = ',xin
!  print*,'xin-xorigin = ',xintemp,' xtemp = ',xtemp
!  print*,'xouttemp = ',xout

  xout(1:2) = xout(1:2) + xorigin(1:2)
!  print*,'xout = ',xout

  !
  !--rotation about x axis (tilt)
  !

  return
end subroutine rotate
