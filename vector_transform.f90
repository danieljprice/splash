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
  real, intent(out), dimension(ndimout) :: xout,vecout
  real, dimension(3,3) :: dxdx
  real :: sinphi, cosphi, sintheta, costheta
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
        if (itypeout(1:3).ne.'car') then
           print*,'output option not implemented: using cartesian output'
        endif
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
        if (itypeout(1:3).ne.'car') then
           print*,'output option not implemented: using cartesian output'
        endif
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
     case(2)
        !
        !--output is cylindrical
        !
        rr = sqrt(dot_product(xin(1:max(ndimin,2)),xin(1:max(ndimin,2)))
	if (rr.ne.0.) then
	   rr1 = 1./rr
        else
	   rr1 = 0.
	endif        
        dxdx(1,1) = xin(1)*rr1  ! dr/dx
        if (ndimin.ge.2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
        if (ndimout.ge.2) then
           term = xin(1)*(1. + (xin(2)/xin(1))**2)
           dxdx(2,1) = -xin(2)/(xin(1)*term) ! dphi/dx
           dxdx(2,2) = 1./term               ! dphi/dy
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
