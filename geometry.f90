!-----------------------------------------------------------------
! Module containing subroutines to transform between 
! different co-ordinate systems, for co-ordinates and vectors
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! itype must be one of the following:
!  itype = 1    : cartesian (default)
!  itype = 2    : cylindrical
!  itype = 3    : spherical
!
! Currently handles:
!
!  cartesian -> cylindrical, spherical polar
!  cylindrical -> cartesian
!  spherical polar -> cartesian
!
!-----------------------------------------------------------------
module geometry
 implicit none
 integer, parameter :: maxcoordsys = 3
 real, parameter, private :: pi = 3.1415926536
 character(len=20), dimension(maxcoordsys), parameter :: labelcoordsys = &
    (/'cartesian         ', &
      'cylindrical polars', &
      'spherical polars  '/)
 character(len=*), dimension(3,maxcoordsys), parameter :: labelcoord = &
    reshape((/'x    ','y    ','z    ', &
              'r    ','phi  ','z    ', &
              'r    ','theta','phi  '/),shape=(/3,maxcoordsys/))
 
contains
!-----------------------------------------------------------------
! Subroutine to transform between different co-ordinate systems
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! xin(ndimin)   : input co-ordinates, in ndimin dimensions
! itypein       : input co-ordinate type
!
! xout(ndimout) : output co-ordinates, in ndimout dimensions
! itypeout      : output co-ordinate type
!

!
!-----------------------------------------------------------------
subroutine coord_transform(xin,ndimin,itypein,xout,ndimout,itypeout)
  implicit none
  integer, intent(in) :: ndimin,ndimout,itypein,itypeout
  real, intent(in), dimension(ndimin) :: xin
  real, intent(out), dimension(ndimout) :: xout
!
!--check for errors in input
!
  if (itypeout.eq.itypein) then
     xout(1:ndimout) = xin(1:ndimout)
     return
  elseif (ndimin.lt.1.or.ndimin.gt.3) then
     print*,'Error: coord transform: invalid number of dimensions on input'
     return
  elseif (ndimout.lt.1.or.ndimout.gt.3) then
     print*,'Error: coord transform: invalid number of dimensions on output'
     return
  elseif (ndimout.gt.ndimin) then
     print*,'Error: coord transform: ndimout must be <= ndimin'
     return
  elseif (abs(xin(1)).lt.1e-8 .and. ndimout.ge.2 .and. &
       (itypein.eq.2 .or. itypein.eq.3)) then
     print*,'Warning: coord transform: r=0 on input: cannot return angle'
     xout(1:ndimout) = xin(1:ndimout)
     return
  endif
!
!--now do transformation
!
  select case(itypein)
!
!--input is cylindrical polars
!
  case(2)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
        if (itypeout.ne.1) print*,'warning: using default cartesian output'
        if (ndimout.eq.1) then
           xout(1) = xin(1)
        else  ! r,phi,z -> x,y,z
           xout(1) = xin(1)*COS(xin(2))
           xout(2) = xin(1)*SIN(xin(2))
           if (ndimout.gt.2) xout(3) = xin(3)
        endif
     end select
!
!--input is spherical polars
!
  case(3)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
        if (itypeout.ne.1) print*,'warning: using default cartesian output'
        select case(ndimout)
           case(1) ! r -> x
              xout(1) = xin(1)
           case(2) ! r,phi -> x,y
              xout(1) = xin(1)*COS(xin(2))
              xout(2) = xin(1)*SIN(xin(2))
           case(3) ! r,phi,theta -> x,y,z
              xout(1) = xin(1)*COS(xin(2))*SIN(xin(3))
              xout(2) = xin(1)*SIN(xin(2))*SIN(xin(3))
              xout(3) = xin(1)*COS(xin(3))
        end select
     end select
!
!--input is cartesian co-ordinates
!
  case default
     select case(itypeout)
     case(2)
        !
        !--output is cylindrical
        !
        if (ndimin.eq.1) then
           xout(1) = abs(xin(1))   ! cylindrical r
        else
           xout(1) = SQRT(DOT_PRODUCT(xin(1:2),xin(1:2)))
           if (ndimout.ge.2) xout(2) = ATAN2(xin(2),xin(1)) ! phi
           if (ndimout.eq.3) xout(3) = xin(3) ! z
        endif
     case(3)
        !
        ! output is spherical
        !
        xout(1) = SQRT(DOT_PRODUCT(xin,xin))! r  
        if (ndimout.ge.2) xout(2) = ATAN2(xin(2),xin(1)) ! phi
        if (ndimout.ge.3) then
           !
           ! theta = ACOS(z/r)
           !
           xout(3) = ACOS(xin(3)/xout(1))
           !--sort out which quadrant for theta
           !if (xin(3).lt.0. .and. xin(1).lt.0.) then
           !  xout(3) = xout(3) + pi
           !lseif (xin(3).lt.0.) then
           !!   xout(3) = pi - xout(3)
           !elseif (xin(1).lt.0.) then
           !   xout(3) = -xout(3)
           !endif    
        endif
     case default
        !
        ! just copy
        !
        xout(1:ndimout) = xin(1:ndimout)
     end select
  end select

  return
end subroutine coord_transform

!-----------------------------------------------------------------
! Subroutine to transform vector components 
!  between different co-ordinate systems
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! Arguments:
!   xin(ndimin)   : input co-ordinates, in ndimin dimensions
!   vecin(ndimin) : components of vector in input co-ordinate basis
!   itypein      : input co-ordinate type
!
!   vecout(ndimout) : components of vector in output co-ordinate basis
!   itypeout      : output co-ordinate type
!
! coords must be one of the following:
!   'cartesian' (default)
!   'cylindrical'
!   'spherical'
!
! Currently handles:
!
!  cartesian -> cylindrical, spherical polar
!  cylindrical -> cartesian
!  spherical polar -> cartesian
!
!-----------------------------------------------------------------
subroutine vector_transform(xin,vecin,ndimin,itypein,vecout,ndimout,itypeout)
  implicit none
  integer, intent(in) :: ndimin,ndimout,itypein,itypeout
  real, intent(in), dimension(ndimin) :: xin,vecin
  real, intent(out), dimension(ndimout) :: vecout
  integer :: i
  real, dimension(3,3) :: dxdx
  real :: sinphi, cosphi
  real :: rr,rr1
!
!--check for errors in input
!
  if (ndimout.gt.ndimin) then
     print*,'Error: vec transform: ndimout must be <= ndimin'
     return
  elseif (itypein.eq.itypeout) then
     vecout(1:ndimout) = vecin(1:ndimout)
     return
  elseif (ndimin.lt.1.or.ndimin.gt.3) then
     print*,'Error: vec transform: invalid number of dimensions on input'
     return
  elseif (ndimout.lt.1.or.ndimout.gt.3) then
     print*,'Error: vec transform: invalid number of dimensions on output'
     return
  elseif (abs(xin(1)).lt.1e-8 .and. &
       (itypein.eq.2 .or. itypein.eq.3)) then
     print*,'Error: vec transform: r=0 on input'
     return
  endif
!
!--set Jacobian matrix to zero
!
  dxdx = 0.
!
!--calculate non-zero components of Jacobian matrix for the transformation
!
  select case(itypein)
!
!--input is spherical polars
!
  case(3)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
        dxdx(1,1) = SIN(xin(2))*COS(xin(3))         ! dx/dr
        dxdx(1,2) = xin(1)*COS(xin(2))*COS(xin(3))  ! dx/dtheta
        dxdx(1,3) = -xin(1)*SIN(xin(2))*SIN(xin(3)) ! dx/dphi
        dxdx(2,1) = SIN(xin(2))*SIN(xin(3))         ! dy/dr
        dxdx(2,2) = xin(1)*COS(xin(2))*SIN(xin(3))  ! dy/dtheta
        dxdx(2,3) = xin(1)*SIN(xin(2))*COS(xin(3))  ! dy/dphi
        dxdx(3,1) = COS(xin(2))                     ! dz/dr
        dxdx(3,2) = -xin(1)*SIN(xin(2))             ! dz/dtheta
     end select
!
!--input is cylindrical polars
!
  case(2)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
	sinphi = SIN(xin(2))
	cosphi = COS(xin(2))
        dxdx(1,1) = cosphi            ! dx/dr
        dxdx(1,2) = -xin(1)*sinphi    ! dx/dphi
        dxdx(2,1) = sinphi            ! dy/dr
        dxdx(2,2) = xin(1)*cosphi     ! dy/dphi
        dxdx(3,3) = 1.                ! dz/dz
     end select
!
!--input is cartesian co-ordinates (default)
!
  case default
     select case(itypeout)
     case(3)
        !
        ! output is spherical
        !
        rr = sqrt(dot_product(xin,xin))
	if (rr.ne.0.) then
	   rr1 = 1./rr
        else
	   rr1 = 0.
	endif
           dxdx(1,1) = xin(1)*rr1  ! dr/dx
           if (ndimin.ge.2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
           if (ndimin.eq.3) dxdx(1,3) = 1.          ! dr/dz 
           if (ndimout.ge.2) then
              ! FILL THESE IN
              dxdx(2,1) = 0. ! dtheta/dx
              dxdx(2,2) = 0. ! dtheta/dy
              dxdx(2,3) = 0. ! dtheta/dz
              dxdx(3,1) = 0. ! dphi/dx
              dxdx(3,2) = 0. ! dphi/dy
              dxdx(3,3) = 0. ! dphi/dz
           endif
     case(2)
        !
        !--output is cylindrical
        !
        rr = sqrt(dot_product(xin(1:max(ndimin,2)),xin(1:max(ndimin,2))))
	if (rr.ne.0.) then
	   rr1 = 1./rr
        else
	   rr1 = 0.
	endif        
        dxdx(1,1) = xin(1)*rr1  ! dr/dx
        if (ndimin.ge.2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
        if (ndimout.ge.2) then
           dxdx(2,1) = -xin(2)*rr1**2 ! dphi/dx
           dxdx(2,2) = xin(1)*rr1**2  ! dphi/dy
           if (ndimout.eq.3) dxdx(3,3) = 1.  ! dz/dz 
        endif
     case default
        print*,'coord transform: invalid co-ordinate type on output'
        vecout(1:ndimout) = vecin(1:ndimout)
        return
     end select
  end select
!
!--now perform transformation using Jacobian matrix
!
  do i=1,ndimout
     vecout(i) = dot_product(dxdx(i,1:ndimin),vecin(1:ndimin))
  enddo

  return
end subroutine vector_transform

!------------------------------------------------------------------
! this subroutine attempts to switch plot limits / boundaries 
! between various co-ordinate systems.
!------------------------------------------------------------------
subroutine coord_transform_limits(xmin,xmax,itypein,itypeout,ndim)
 implicit none
 integer, intent(in) :: itypein,itypeout,ndim
 real, dimension(ndim), intent(inout) :: xmin,xmax
 real, dimension(ndim) :: xmaxtemp,xmintemp
!
!--check for errors in input
!
 if (ndim.lt.1 .or. ndim.gt.3) then
    print*,'Error: limits coord transform: ndim invalid on input'
    return
 endif 
!
!--by default do nothing
!
 xmintemp(1:ndim) = xmin(1:ndim)
 xmaxtemp(1:ndim) = xmax(1:ndim)

 select case(itypein)
!
!--input is spherical
!
 case(3)
    select case(itypeout)
    case default
    !
    !--cartesian output
    !
    xmintemp(1:ndim) = -xmax(1)
    xmaxtemp(1:ndim) = xmax(1)

    end select
!
!--input is cylindrical
!
 case(2)
    select case(itypeout)
    case default
    !
    !--cartesian output
    !
    xmintemp(1:max(ndim,2)) = -xmax(1)
    xmaxtemp(1:max(ndim,2)) = xmax(1)
        
    end select
!
!--input is cartesian
!
 case default
    select case(itypeout)
    !
    !--output is spherical
    !
    case(3)
       !--rmin, rmax
       xmintemp(1) = 0.
       xmaxtemp(1) = max(maxval(abs(xmin(1:ndim))), &
                         maxval(abs(xmax(1:ndim))))
       if (ndim.ge.2) then
          xmintemp(2) = -pi
          xmaxtemp(2) = pi
          if (ndim.ge.3) then
             xmintemp(3) = -0.5*pi
             xmaxtemp(3) = 0.5*pi
          endif
       endif
    !
    !--output is cylindrical
    !
    case(2)
       !--rmin, rmax
       xmintemp(1) = 0.
       xmaxtemp(1) = max(maxval(abs(xmin(1:max(2,ndim)))), &
                         maxval(abs(xmax(1:max(2,ndim)))))
       if (ndim.ge.2) then
          xmintemp(2) = -pi
          xmaxtemp(2) = pi
       endif
    end select
 end select

 xmin(:) = min(xmintemp(:),xmaxtemp(:))
 xmax(:) = max(xmintemp(:),xmaxtemp(:))
 
 return
end subroutine coord_transform_limits

end module geometry
