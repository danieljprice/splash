!-----------------------------------------------------------------
! Standalone module containing subroutines to transform between 
! different co-ordinate systems, for co-ordinates and vectors
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! itype must be one of the following:
!  itype = 1    : cartesian (default)
!  itype = 2    : cylindrical
!  itype = 3    : spherical
!  itype = 4    : toroidal
!
! Currently handles:
!
!  cartesian -> cylindrical, spherical polar
!  cylindrical -> cartesian
!  spherical polar -> cartesian
!  toroidal r,theta,phi <-> cartesian
!
! written by Daniel Price 2004-2007
! as part of the SPLASH SPH visualisation package
!-----------------------------------------------------------------
module geometry
 implicit none
 integer, parameter, public :: maxcoordsys = 4

 character(len=*), dimension(maxcoordsys), parameter, public :: labelcoordsys = &
    (/'cartesian   x,y,z      ', &
      'cylindrical r,phi,z    ', &
      'spherical   r,phi,theta', &
      'toroidal    r,theta,phi'/)
 character(len=*), dimension(3,maxcoordsys), parameter, public :: labelcoord = &
    reshape((/'x    ','y    ','z    ', &
              'r    ','phi  ','z    ', &
              'r    ','phi  ','theta', &
              'r_tor','theta','phi  '/),shape=(/3,maxcoordsys/))

 public :: coord_transform, vector_transform, coord_transform_limits

 real, parameter, private :: pi = 3.1415926536
 real, parameter, private :: Rtorus = 1.0
 
 private
 
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
  real :: rcyl
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
!--input is torus co-ordinates
!
  case(4)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
        if (itypeout.ne.1) print*,'warning: using default cartesian output'
        if (ndimin.ne.3) then
           xout(1:ndimout) = xin(1:ndimout)
        else
           rcyl = xin(1)*COS(xin(2)) + Rtorus
           xout(1) = rcyl*COS(xin(3))
           if (ndimout.ge.2) xout(2) = rcyl*SIN(xin(3))
           if (ndimout.ge.3) xout(3) = xin(1)*SIN(xin(2))
        endif
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
        !--output is spherical
        !
        xout(1) = SQRT(DOT_PRODUCT(xin,xin))! r  
        if (ndimout.ge.2) xout(2) = ATAN2(xin(2),xin(1)) ! phi
        if (ndimout.ge.3) then
           ! theta = ACOS(z/r)
           xout(3) = ACOS(xin(3)/xout(1))
        endif
     case(4)
        !
        !--output is torus r,theta,phi co-ordinates
        !
        if (ndimin.ne.3) then
           ! not applicable if ndim < 3
           xout(1:ndimout) = xin(1:ndimout)
        else
           rcyl = SQRT(xin(1)**2 + xin(2)**2)
           xout(1) = SQRT(xin(3)**2 + (rcyl - Rtorus)**2)
           if (ndimout.ge.2) xout(2) = ASIN(xin(3)/xout(1))
           if (ndimout.ge.3) xout(3) = ATAN2(xin(2),xin(1))
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
  real :: rr,rr1,rcyl,rcyl2,rcyl1
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
     print*,'Warning: vec transform: r=0 on input, setting vec = 0 at origin'
     vecout = 0.
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
!--input is toroidal
!
  case(4)
     select case(itypeout)
     case default  
        dxdx(1,1) = COS(xin(2))*COS(xin(3))         ! dx/dr
        dxdx(1,2) = -SIN(xin(2))*COS(xin(3))        ! 1/r dx/dtheta
        dxdx(1,3) = SIN(xin(3))                     ! 1/rcyl dx/dphi
        dxdx(2,1) = COS(xin(2))*SIN(xin(3))         ! dy/dr
        dxdx(2,2) = -SIN(xin(2))*SIN(xin(3))        ! 1/r dy/dtheta
        dxdx(2,3) = COS(xin(3))                     ! 1/rcyl dy/dphi
        dxdx(3,1) = SIN(xin(3))                     ! dz/dr
        dxdx(3,2) = COS(xin(3))                     ! 1/r dz/dtheta
!        dxdx(3,3) = 0.                             ! dz/dphi
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
        dxdx(1,1) = COS(xin(2))*SIN(xin(3))         ! dx/dr
        dxdx(1,2) = -SIN(xin(2))                    ! 1/rcyl dx/dphi
        dxdx(1,3) = COS(xin(2))*COS(xin(3))         ! 1/r dx/dtheta
        dxdx(2,1) = SIN(xin(2))*SIN(xin(3))         ! dy/dr
        dxdx(2,2) = COS(xin(2))                     ! 1/rcyl dy/dphi
        dxdx(2,3) = SIN(xin(2))*COS(xin(3))         ! 1/r dy/dtheta
        dxdx(3,1) = COS(xin(3))                     ! dz/dr
        dxdx(3,3) = -SIN(xin(3))                    ! 1/r dz/dtheta
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
        dxdx(1,2) = -sinphi           ! 1/r*dx/dphi
        dxdx(2,1) = sinphi            ! dy/dr
        dxdx(2,2) = cosphi            ! 1/r*dy/dphi
        dxdx(3,3) = 1.                ! dz/dz
     end select
!
!--input is cartesian co-ordinates (default)
!
  case default
     select case(itypeout)
     case(4)
        !
        ! output is toroidal
        !
        rcyl = sqrt(xin(1)**2 + xin(2)**2)
        if (rcyl.gt.tiny(rcyl)) then
           rcyl1 = 1./rcyl
        else
           rcyl1 = 0.
        endif
        rr = sqrt((rcyl - Rtorus)**2 + xin(3)**2)
        if (rr.gt.tiny(rr)) then
           rr1 = 1./rr
        else
           rr1 = 0.
        endif
        dxdx(1,1) = (rcyl - Rtorus)*xin(1)*rr1*rcyl1 ! dr/dx
        dxdx(1,2) = (rcyl - Rtorus)*xin(2)*rr1*rcyl1 ! dr/dy
        dxdx(1,3) = xin(3)*rr1                       ! dr/dz
        dxdx(2,1) = -xin(3)*xin(1)*rr1*rcyl1         ! dtheta/dx
        dxdx(2,2) = -xin(3)*xin(2)*rr1*rcyl1         ! dtheta/dy
        dxdx(2,3) = (rcyl - Rtorus)*rr1              ! dtheta/dz
        dxdx(3,1) = -xin(2)*rcyl1                    ! dphi/dx
        dxdx(3,2) = xin(1)*rcyl1                     ! dphi/dy
!        dxdx(3,3) = 0.                              ! dphi/dz
     case(3)
        !
        ! output is spherical
        !
        rr = sqrt(dot_product(xin,xin))
        if (rr.gt.tiny(rr)) then
           rr1 = 1./rr
        else
           rr1 = 0.
        endif
        dxdx(1,1) = xin(1)*rr1  ! dr/dx
        if (ndimin.ge.2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
        if (ndimin.eq.3) dxdx(1,3) = xin(3)*rr1  ! dr/dz 
        if (ndimin.ge.2) then
           rcyl2 = dot_product(xin(1:2),xin(1:2))
           rcyl = sqrt(rcyl2)
           if (rcyl.gt.tiny(rcyl)) then
              rcyl1 = 1./rcyl
           else
              rcyl1 = 0.
           endif
           dxdx(2,1) = -xin(2)*rcyl1 ! rcyl dphi/dx
           dxdx(2,2) = xin(1)*rcyl1  ! rcyl dphi/dy
           dxdx(2,3) = 0.
           if (ndimin.ge.3) then
              dxdx(3,1) = xin(1)*xin(3)*rr1*rcyl1 ! r dtheta/dx
              dxdx(3,2) = xin(2)*xin(3)*rr1*rcyl1 ! r dtheta/dy
              dxdx(3,3) = -rcyl2*rr1*rcyl1 ! r dtheta/dz
           endif
        endif
     case(2)
        !
        !--output is cylindrical
        !
        rr = sqrt(dot_product(xin(1:min(ndimin,2)),xin(1:min(ndimin,2))))
        if (rr.ne.0.) then
           rr1 = 1./rr
        else
           rr1 = 0.
        endif        
        dxdx(1,1) = xin(1)*rr1  ! dr/dx
        if (ndimin.ge.2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
        if (ndimout.ge.2) then
           dxdx(2,1) = -xin(2)*rr1 ! r*dphi/dx
           dxdx(2,2) = xin(1)*rr1  ! r*dphi/dy
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
 print*,'modifying plot limits for new coordinate system'
!
!--by default do nothing
!
 xmintemp(1:ndim) = xmin(1:ndim)
 xmaxtemp(1:ndim) = xmax(1:ndim)

 select case(itypein)
!
!--input is toroidal
!
 case(4)
    select case(itypeout)
    case default
    !
    !--cartesian output
    !
    xmintemp(1:min(ndim,2)) = -Rtorus - xmax(1)
    xmaxtemp(1:min(ndim,2)) = Rtorus + xmax(1)
    if (ndim.eq.3) then
       xmintemp(3) = -xmax(1)
       xmaxtemp(3) = xmax(1)
    endif
    
    end select
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
    case(4)
    !
    !--output is toroidal
    !
    xmintemp(1) = 0.
    xmaxtemp(1) = max(maxval(abs(xmax(1:min(ndim,2))))-Rtorus, &
                      maxval(abs(xmin(1:min(ndim,2))))-Rtorus)
    if (ndim.ge.2) then
       xmintemp(2) = -0.5*pi
       xmaxtemp(2) = 0.5*pi
       if (ndim.ge.3) then
          xmintemp(3) = -pi
          xmaxtemp(3) = pi
       endif
    endif
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
             xmintemp(3) = 0.
             xmaxtemp(3) = pi
          endif
       endif
    !
    !--output is cylindrical
    !
    case(2)
       !--rmin, rmax
       xmintemp(1) = 0.
       if (ndim.ge.2) then
          xmaxtemp(1) = max(maxval(abs(xmin(1:max(2,ndim)))), &
                         maxval(abs(xmax(1:max(2,ndim)))))
          xmintemp(2) = -pi
          xmaxtemp(2) = pi
       else
          xmaxtemp(1) = max(abs(xmin(1)),abs(xmax(1)))
       endif
    end select
 end select

 xmin(:) = min(xmintemp(:),xmaxtemp(:))
 xmax(:) = max(xmintemp(:),xmaxtemp(:))
 
 return
end subroutine coord_transform_limits

end module geometry
