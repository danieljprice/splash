!
!--plots three dimensional axes, as box or as axes
!
!  arguments:
!       angles(ndim-1) : azimuthal and (in 3D) tilt angle in radians
!       xin(ndim)      : input co-ordinates
!       xout(ndim)     : output (rotated) co-ordinates
!       xorigin(ndim)  : input co-ordinates of the origin
!
subroutine rotate_axes(angles,xmin,xmax,xorigin,ndim)
  implicit none
  integer, intent(in) :: ndim
  real, intent(in), dimension(ndim-1) :: angles
  real, intent(in), dimension(ndim) :: xmin,xmax,xorigin
  integer :: i
  real, dimension(ndim,2**ndim) :: xpt
  real, dimension(2) :: xline, yline

  if (ndim.lt.2 .or. ndim.gt.3) then
     print*,'error: rotate_axes: ndim < 2 or > 3'
     return
  endif

  print*,'plotting rotated axes...'
!
!--set coordinates of all 4(8) box corners
!

!--front face (pts 1->4)
  xpt(:,1) = xmin(:)  ! xmin, ymin

  xpt(1,2) = xmin(1)  ! xmin
  xpt(2,2) = xmax(2)  ! ymax

  xpt(1,3) = xmax(1)  ! xmax
  xpt(2,3) = xmax(2)  ! ymax

  xpt(1,4) = xmax(1)  ! xmax
  xpt(2,4) = xmin(2)  ! ymin

  if (ndim.eq.3) xpt(ndim,1:4) = xmin(ndim)
  
!--back face (pts 5->8)
  do i=1,4
     xpt(1:2,i+4) = xpt(1:2,i)
  enddo
  if (ndim.eq.3) xpt(ndim,5:8) = xmax(ndim)
!
!--now rotate each of these coordinates
!
  do i=1,2**ndim
     !!print*,' corner ',i,' x,y,z = ',xpt(:,i)
     call rotate(angles,xpt(:,i),xpt(:,i),xorigin,ndim)
  enddo
!
!--now plot boxes appropriately using points
!
  call pgsfs(2)

  !--front face (points 1->4)
  call pgpoly(4,xpt(1,1:4),xpt(2,1:4))

  if (ndim.eq.3) then
     !--back face (points 5->8)
     call pgpoly(4,xpt(1,5:8),xpt(2,5:8))
     
     !--connecting lines ( 1->5, 2->6, 3->7, 4->8 )
     do i=1,4
        xline(1) = xpt(1,i)
        yline(1) = xpt(2,i)
        xline(2) = xpt(1,i+4)
        yline(2) = xpt(2,i+4)
        call pgline(2,xline,yline)
     enddo
  endif
  
  return
end subroutine rotate_axes
