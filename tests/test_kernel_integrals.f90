!-----------------------------------------------------------------------
!+
!  Original cubic-spline-only implementation of the exact kernel
!  integrals (Petkova, Laibe & Bonnell 2018), kept as a reference
!  for regression testing of the generic implementation
!+
!-----------------------------------------------------------------------
module legacy_cubic
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer, parameter :: doub_prec = kind(0.d0)
 public :: pint_legacy,wallint_legacy
 private

contains

!---------------------------------------
!
!  Integrals of the kernel, used
!  when performing exact interpolation
!  See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
pure real function pint_legacy(r0, d1, d2, hi1)
 real, intent(in) :: r0, d1, d2, hi1
 real :: ar0,q0,tphi1,tphi2,phi1,phi2

 if (abs(r0) < tiny(0.)) then
    pint_legacy = 0.
    return
 elseif (r0 > 0.) then
    pint_legacy = 1.
    ar0 = r0
 else
    pint_legacy = -1.
    ar0 = -r0
 endif
 q0 = ar0*hi1
 tphi1 = abs(d1)/ar0
 tphi2 = abs(d2)/ar0
 phi1 = atan(tphi1)
 phi2 = atan(tphi2)

 if (d1*d2 >= 0.) then
    pint_legacy = pint_legacy*(full_2d_mod(phi1,tphi1,q0) + full_2d_mod(phi2,tphi2,q0))
 else
    if (abs(d1) < abs(d2)) then
       pint_legacy = pint_legacy*(full_2d_mod(phi2,tphi2,q0) - full_2d_mod(phi1,tphi1,q0))
    else
       pint_legacy = pint_legacy*(full_2d_mod(phi1,tphi1,q0) - full_2d_mod(phi2,tphi2,q0))
    endif
 endif

end function pint_legacy

!---------------------------------------
!
! Helper functions for kernel integrals
! See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
pure real function full_2d_mod(phi,tphi,q0)
 real, intent(in) :: phi,tphi,q0
 real :: q, phi1, phi2, tphi1, tphi2, cphi, cphi1, cphi2

 if (q0 <= 1.0) then
    cphi = cos(phi)
    q = q0/cphi

    if (q <= 1.0) then
       full_2d_mod = F1_2d(phi,tphi,cphi,q0)
    elseif (q <= 2.0) then
       cphi1 = q0
       phi1 = acos(q0)
       tphi1 = tan(phi1)  !  tan (acos x) = sqrt(1-x**2)/x
       full_2d_mod = F2_2d(phi,tphi,cphi,q0) - F2_2d(phi1,tphi1,cphi1,q0) &
                   + F1_2d(phi1,tphi1,cphi1,q0)
    else
       cphi1 = q0
       phi1 = acos(q0)
       cphi2 = 0.5*q0
       phi2 = acos(0.5*q0)
       tphi1 = tan(phi1)
       tphi2 = tan(phi2)
       full_2d_mod = F3_2d(phi) - F3_2d(phi2) + F2_2d(phi2,tphi2,cphi2,q0) &
                                - F2_2d(phi1,tphi1,cphi1,q0) &
                                + F1_2d(phi1,tphi1,cphi1,q0)
    endif
 elseif (q0 <= 2.0) then
    cphi = cos(phi)
    q = q0/cphi

    if (q <= 2.0) then
       full_2d_mod = F2_2d(phi,tphi,cphi,q0)
    else
       cphi2 = 0.5*q0
       phi2 = acos(0.5*q0)
       tphi2 = tan(phi2)
       full_2d_mod = F3_2d(phi) - F3_2d(phi2) + F2_2d(phi2,tphi2,cphi2,q0)
    endif
 else ! q0 > 2
    full_2d_mod = F3_2d(phi)
 endif

end function full_2d_mod

pure real function F1_2d(phi,tphi,cphi,q0)
 real, intent(in) :: phi,tphi,cphi,q0
 real :: I2, I4, I5, logs, cphi2, q02, q03

 !tphi = tan(phi)
 !cphi = cos(phi)
 cphi2 = cphi*cphi

 q02 = q0*q0
 q03 = q02*q0

 logs = log(tan(phi/2.+pi/4.))

 I2 = tphi
 I4 = 1./3. * tphi*(2. + 1./cphi2)

 I5 = 1./16. * (0.5*(11.*sin(phi) + 3.*sin(3.*phi))/cphi2/cphi2 + 6.*logs)

 F1_2d =  5./7.*q02/pi  * (I2 - 3./4.*q02 *I4 + 0.3*q03 *I5)

end function F1_2d

pure real function F2_2d(phi, tphi, cphi, q0)
 real, intent(in) :: phi, tphi, cphi, q0
 real :: I0, I2, I3, I4, I5, logs, cphi2, q02, q03

! tphi = tan(phi)
! cphi = cos(phi)
 cphi2 = cphi*cphi

 q02 = q0*q0
 q03 = q02*q0

 logs = log(tan(phi/2.+pi/4.))

 I0 = phi
 I2 = tphi
 I4 = 1./3. * tphi*(2. + 1./(cphi2))

 I3 = 1./2.  * (tphi/cphi + logs)
 I5 = 1./16. * (0.5*(11.*sin(phi) + 3.*sin(3.*phi))/cphi2/cphi2 + 6.*logs)

 F2_2d =  5./7.*q02/pi  * (2.*I2 - 2.*q0 *I3 + 3./4.*q02 *I4 - 1./10.*q03 *I5 &
                           - 1./10./q02 *I0)

end function F2_2d

pure real function F3_2d(phi)
 real, intent(in) :: phi
 real :: I0

 I0 = phi
 F3_2d = 0.5/pi  *I0

end function F3_2d

!------------------------------------------------------------
! 3D functions to evaluate exact overlap of kernel with wall boundaries
! see Petkova, Laibe & Bonnell (2018), J. Comp. Phys
!------------------------------------------------------------
real function wallint_legacy(r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi)
 real, intent(in) :: r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi
 real(doub_prec) :: R_0, d1, d2, dx, dy, h

 wallint_legacy = 0.0
 dx = xc - xp
 dy = yc - yp
 h = hi

 !
 ! Contributions from each of the 4 sides of a cell wall
 !
 R_0 = 0.5*pixwidthy + dy
 d1 = 0.5*pixwidthx - dx
 d2 = 0.5*pixwidthx + dx
 wallint_legacy = wallint_legacy + pint3D_legacy(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthy - dy
 d1 = 0.5*pixwidthx + dx
 d2 = 0.5*pixwidthx - dx
 wallint_legacy = wallint_legacy + pint3D_legacy(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthx + dx
 d1 = 0.5*pixwidthy + dy
 d2 = 0.5*pixwidthy - dy
 wallint_legacy = wallint_legacy + pint3D_legacy(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthx - dx
 d1 = 0.5*pixwidthy - dy
 d2 = 0.5*pixwidthy + dy
 wallint_legacy = wallint_legacy + pint3D_legacy(r0, R_0, d1, d2, h)

end function wallint_legacy


real function pint3D_legacy(r0, R_0, d1, d2, hi)

 real(doub_prec), intent(in) :: R_0, d1, d2, hi
 real, intent(in) :: r0
 real(doub_prec) :: ar0, aR_0
 real(doub_prec) :: int1, int2

 if (abs(r0) < tiny(0.)) then
    pint3D_legacy = 0.d0
    return
 endif

 if (r0  >  0.d0) then
    pint3D_legacy = 1.d0
    ar0 = r0
 else
    pint3D_legacy = -1.d0
    ar0 = -r0
 endif

 if (R_0  >  0.d0) then
    aR_0 = R_0
 else
    pint3D_legacy = -pint3D_legacy
    aR_0 = -R_0
 endif

 int1 = full_integral_3D(d1, ar0, aR_0, hi)
 int2 = full_integral_3D(d2, ar0, aR_0, hi)

 if (int1 < 0.d0) int1 = 0.d0
 if (int2 < 0.d0) int2 = 0.d0

 if (d1*d2  >=  0) then
    pint3D_legacy = pint3D_legacy*real(int1 + int2)
    if (int1 + int2 < 0.d0) pint3D_legacy = 0. !print*, 'Error: int1 + int2 < 0'
 elseif (abs(d1)  <  abs(d2)) then
    pint3D_legacy = pint3D_legacy*real(int2 - int1)
    if (int2 - int1 < 0.d0) pint3D_legacy = 0.! print*, 'Error: int2 - int1 < 0: ', int1, int2, '(', d1, d2,')'
 else
    pint3D_legacy = pint3D_legacy*real(int1 - int2)
    if (int1 - int2 < 0.d0) pint3D_legacy = 0. !print*, 'Error: int1 - int2 < 0: ', int1, int2, '(', d1, d2,')'
 endif

end function pint3D_legacy

real(doub_prec) function full_integral_3D(d, r0, R_0, h)

 real(doub_prec), intent(in) :: d, r0, R_0, h
 real(doub_prec) :: B1, B2, B3, a, h2
 real(doub_prec), parameter :: pi = 4.*atan(1.)
 real(doub_prec) :: tanphi, phi, a2, cosp, r0h, r03, r0h2, r0h3, r0h_2, r0h_3
 real(doub_prec) :: r2, R_, linedist2, cosphi
 real(doub_prec) :: I0, I1, I_2, I_3, I_4, I_5
 real(doub_prec) :: D2, D3

 r0h = r0/h
 tanphi = abs(d)/R_0
 phi = atan(tanphi)
 B1 = 0.d0
 B2 = 0.d0
 B3 = 0.d0
 a = 0.d0
 h2 = 0.d0

 if (abs(r0h) < tiny(0.) .or. abs(R_0/h) < tiny(0.) .or. abs(phi) < tiny(0.)) then
    full_integral_3D = 0.0
    return
 endif

 h2 = h*h
 r03 = r0*r0*r0
 r0h2 = r0h*r0h
 r0h3 = r0h2*r0h
 r0h_2 = 1./r0h2
 r0h_3 = 1./r0h3

 if (r0 >= 2.0*h) then
    B3 = 0.25*h2*h
 elseif (r0 > h) then
    B3 = 0.25*r03 *(-4./3. + (r0h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2)
    B2 = 0.25*r03 *(-4./3. + (r0h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3)
 else
    B3 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2)
    B2 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2)
    B1 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3)
 endif

 a = R_0/r0
 a2 = a*a

 linedist2 = (r0*r0 + R_0*R_0)
 cosphi = cos(phi)
 R_ = R_0/cosphi
 r2 = (r0*r0 + R_*R_)

 D2 = 0.0
 D3 = 0.0

 if (linedist2 < h2) then
    !////// phi1 business /////
    cosp = R_0/sqrt(h2-r0*r0)
    call get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)

    D2 = -1./6.*I_2 + 0.25*(r0h) *I_3 - 0.15*r0h2 *I_4 + 1./30.*r0h3 *I_5 - 1./60. *r0h_3 *I1 + (B1-B2)/r03 *I0
 endif
 if (linedist2 < 4.*h2) then
    !////// phi2 business /////
    cosp = R_0/sqrt(4.0*h2-r0*r0)
    call get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)

    D3 = 1./3.*I_2 - 0.25*(r0h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2
 endif

 !//////////////////////////////
 call get_I_terms(cosphi,a2,a,I0,I1,I_2,I_3,I_4,I_5,phi=phi,tanphi=tanphi)

 if (r2 < h2) then
    full_integral_3D = r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0)
 elseif (r2 < 4.*h2) then
    full_integral_3D=  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - &
    &   1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2)
 else
    full_integral_3D = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3)
 endif

end function full_integral_3D

subroutine get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5,phi,tanphi)
 real(doub_prec), intent(in) :: cosp,a2,a
 real(doub_prec), intent(out) :: I0,I1,I_2,I_3,I_4,I_5
 real(doub_prec), intent(in), optional :: phi,tanphi
 real(doub_prec) :: cosp2,p,tanp,u2,u,logs,I_1,mu2_1,fac

 cosp2 = cosp*cosp
 if (present(phi)) then
    p = phi
    tanp = tanphi
 else
    p = acos(cosp)
    tanp = sqrt(1.-cosp2)/cosp ! tan(p)
 endif

 mu2_1 = 1. / (1. + cosp2/a2)
 I0  = p
 I_2 = p +    a2 * tanp
 I_4 = p + 2.*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2)

 u2 = (1.-cosp2)*mu2_1
 u = sqrt(u2)
 logs = log((1.+u)/(1.-u))
 I1 = atan2(u,a)

 fac = 1./(1.-u2)
 I_1 = 0.5*a*logs + I1
 I_3 = I_1 + a*0.25*(1.+a2)*(2.*u*fac + logs)
 I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10.*u - 6.*u*u2)*fac*fac + 3.*logs)

end subroutine get_I_terms

end module legacy_cubic

!-----------------------------------------------------------------------
!+
!  Helper routines for testing the exact kernel integrals: cell
!  integrals built from the edge/wall integrals, reference values
!  by numerical quadrature, and pass/fail bookkeeping
!+
!-----------------------------------------------------------------------
module kernel_integral_tests
 use, intrinsic :: ieee_arithmetic, only:ieee_is_nan
 use kernels,       only:pint,wallint,wfunc,radkernel,cnormk2D,cnormk3D,ikernel,get_kernel_pieces
 use kernel_integrals, only:maxpieces,maxdeg,maxpow,set_kernel_integral_coeffs
 use legacy_cubic,  only:pint_legacy,wallint_legacy
 use projections3D, only:interpolate3D_projection,setup_integratedkernel
 use projections3Dgeom,     only:interpolate3D_proj_geom
 use interpolate3D_opacity, only:interp3D_proj_opacity
 use geometry,              only:igeom_cartesian
 implicit none
 integer, parameter :: dp = kind(0.d0)
 integer :: ntests = 0, nfail = 0
 ! 4-point Gauss-Legendre nodes and weights on [-1,1]
 real, parameter :: xgl(4) = (/-0.8611363115940526,-0.3399810435848563, &
                                0.3399810435848563, 0.8611363115940526/)
 real, parameter :: wgl(4) = (/ 0.3478548451374538, 0.6521451548625461, &
                                0.6521451548625461, 0.3478548451374538/)

contains

!-----------------------------------------------------------------------
!+
!  check a value against a reference, recording pass/fail
!+
!-----------------------------------------------------------------------
subroutine checkval(label,val,ref,tol)
 character(len=*), intent(in) :: label
 real, intent(in) :: val,ref,tol
 real :: err,tolp

 ! do not ask for more than the working precision allows
 tolp = max(tol,1000.*epsilon(1.)*max(abs(ref),1.))
 ntests = ntests + 1
 err = abs(val - ref)
 if (err > tolp .or. ieee_is_nan(val)) then
    nfail = nfail + 1
    print "(3x,a,a,es23.15,a,es23.15,a,es10.3,a)",label,' = ',val,' should be ',ref,' (err=',err,') FAILED'
 else
    print "(3x,a,a,es10.3,a)",label,' err =',err,' OK'
 endif

end subroutine checkval

!-----------------------------------------------------------------------
!+
!  integral of the 2D kernel of a particle at (x0,y0) over a pixel
!  centred on (xc,yc), summed from the four pixel edges
!+
!-----------------------------------------------------------------------
real function cell2d(x0,y0,xc,yc,wx,wy,h,legacy)
 real,    intent(in) :: x0,y0,xc,yc,wx,wy,h
 logical, intent(in) :: legacy
 real :: dx,dy

 dx = xc - x0
 dy = yc - y0
 cell2d = pint2(0.5*wy - dy,0.5*wx + dx,0.5*wx - dx,1./h,legacy) &
        + pint2(0.5*wy + dy,0.5*wx - dx,0.5*wx + dx,1./h,legacy) &
        + pint2(0.5*wx - dx,0.5*wy - dy,0.5*wy + dy,1./h,legacy) &
        + pint2(0.5*wx + dx,0.5*wy + dy,0.5*wy - dy,1./h,legacy)

end function cell2d

!-----------------------------------------------------------------------
!+
!  call either the generic or legacy 2D edge integral
!+
!-----------------------------------------------------------------------
real function pint2(r0,d1,d2,hi1,legacy)
 real,    intent(in) :: r0,d1,d2,hi1
 logical, intent(in) :: legacy

 if (legacy) then
    pint2 = pint_legacy(r0,d1,d2,hi1)
 else
    pint2 = pint(r0,d1,d2,hi1)
 endif

end function pint2

!-----------------------------------------------------------------------
!+
!  integral of the 3D kernel of a particle at (x0,y0,z0) over a cell
!  centred on (xc,yc,zc), summed from the six cell walls
!+
!-----------------------------------------------------------------------
real function cell3d(x0,y0,z0,xc,yc,zc,wx,wy,wz,h,legacy)
 real,    intent(in) :: x0,y0,z0,xc,yc,zc,wx,wy,wz,h
 logical, intent(in) :: legacy

 cell3d = wall2(0.5*wz + (zc-z0),x0,y0,xc,yc,wx,wy,h,legacy) &
        + wall2(0.5*wz - (zc-z0),x0,y0,xc,yc,wx,wy,h,legacy) &
        + wall2(0.5*wy + (yc-y0),x0,z0,xc,zc,wx,wz,h,legacy) &
        + wall2(0.5*wy - (yc-y0),x0,z0,xc,zc,wx,wz,h,legacy) &
        + wall2(0.5*wx + (xc-x0),z0,y0,zc,yc,wz,wy,h,legacy) &
        + wall2(0.5*wx - (xc-x0),z0,y0,zc,yc,wz,wy,h,legacy)

end function cell3d

!-----------------------------------------------------------------------
!+
!  call either the generic or legacy 3D wall integral
!+
!-----------------------------------------------------------------------
real function wall2(r0,xp,yp,xc,yc,wx,wy,h,legacy)
 real,    intent(in) :: r0,xp,yp,xc,yc,wx,wy,h
 logical, intent(in) :: legacy

 if (legacy) then
    wall2 = wallint_legacy(r0,xp,yp,xc,yc,wx,wy,h)
 else
    wall2 = wallint(r0,xp,yp,xc,yc,wx,wy,h)
 endif

end function wall2

!-----------------------------------------------------------------------
!+
!  composite Gauss-Legendre nodes and weights on the part of [a,b]
!  overlapping [lo,hi], using npanel panels
!+
!-----------------------------------------------------------------------
subroutine get_nodes(a,b,lo,hi,npanel,x,w,n)
 real,    intent(in)  :: a,b,lo,hi
 integer, intent(in)  :: npanel
 real,    intent(out) :: x(:),w(:)
 integer, intent(out) :: n
 real :: x1,x2,dx,xm
 integer :: i,k

 n = 0
 x1 = max(a,lo)
 x2 = min(b,hi)
 if (x2 <= x1) return
 dx = (x2 - x1)/npanel
 do i=1,npanel
    xm = x1 + (i-0.5)*dx
    do k=1,4
       n = n + 1
       x(n) = xm + 0.5*dx*xgl(k)
       w(n) = 0.5*dx*wgl(k)
    enddo
 enddo

end subroutine get_nodes

!-----------------------------------------------------------------------
!+
!  reference 2D cell integral by numerical quadrature
!+
!-----------------------------------------------------------------------
real function quad2d(xc,yc,wx,wy,npanel)
 real,    intent(in) :: xc,yc,wx,wy
 integer, intent(in) :: npanel
 real    :: x(4*npanel),wxq(4*npanel),y(4*npanel),wyq(4*npanel)
 real(dp) :: sum
 integer :: i,j,nx,ny

 call get_nodes(xc-0.5*wx,xc+0.5*wx,-radkernel,radkernel,npanel,x,wxq,nx)
 call get_nodes(yc-0.5*wy,yc+0.5*wy,-radkernel,radkernel,npanel,y,wyq,ny)
 sum = 0._dp
 do j=1,ny
    do i=1,nx
       sum = sum + wxq(i)*wyq(j)*wfunc(x(i)**2 + y(j)**2)
    enddo
 enddo
 quad2d = real(cnormk2D*sum)

end function quad2d

!-----------------------------------------------------------------------
!+
!  reference 3D cell integral by numerical quadrature
!+
!-----------------------------------------------------------------------
real function quad3d(xc,yc,zc,wc,npanel)
 real,    intent(in) :: xc,yc,zc,wc
 integer, intent(in) :: npanel
 real    :: x(4*npanel),wxq(4*npanel),y(4*npanel),wyq(4*npanel),z(4*npanel),wzq(4*npanel)
 real(dp) :: sum
 integer :: i,j,k,nx,ny,nz

 call get_nodes(xc-0.5*wc,xc+0.5*wc,-radkernel,radkernel,npanel,x,wxq,nx)
 call get_nodes(yc-0.5*wc,yc+0.5*wc,-radkernel,radkernel,npanel,y,wyq,ny)
 call get_nodes(zc-0.5*wc,zc+0.5*wc,-radkernel,radkernel,npanel,z,wzq,nz)
 sum = 0._dp
 do k=1,nz
    do j=1,ny
       do i=1,nx
          sum = sum + wxq(i)*wyq(j)*wzq(k)*wfunc(x(i)**2 + y(j)**2 + z(k)**2)
       enddo
    enddo
 enddo
 quad3d = real(cnormk3D*sum)

end function quad3d

!-----------------------------------------------------------------------
!+
!  a cell containing the whole kernel should integrate to unity
!+
!-----------------------------------------------------------------------
subroutine test_normalisation
 real :: w

 w = 2.*radkernel + 1.
 call checkval('2D normalisation, centred ',cell2d(0.,0.,0.,0.,w,w,1.,.false.),1.,1.e-12)
 call checkval('2D normalisation, offset  ',cell2d(0.3,-0.2,0.,0.,w,w,1.,.false.),1.,1.e-12)
 call checkval('3D normalisation, centred ',cell3d(0.,0.,0.,0.,0.,0.,w,w,w,1.,.false.),1.,1.e-12)
 call checkval('3D normalisation, offset  ',cell3d(0.3,-0.2,0.1,0.,0.,0.,w,w,w,1.,.false.),1.,1.e-12)

end subroutine test_normalisation

!-----------------------------------------------------------------------
!+
!  integrals over a grid of cells covering the kernel sum to unity
!+
!-----------------------------------------------------------------------
subroutine test_partition(w)
 real, intent(in) :: w
 real, parameter :: x0 = 0.137, y0 = 0.291, z0 = -0.073
 real(dp) :: sum2,sum3
 real :: xc,yc,zc
 integer :: i,j,k,n
 character(len=40) :: label

 n = int(radkernel/w) + 2
 sum2 = 0._dp
 sum3 = 0._dp
 do j=-n,n
    yc = j*w
    do i=-n,n
       xc = i*w
       sum2 = sum2 + cell2d(x0,y0,xc,yc,w,w,1.,.false.)
       do k=-n,n
          zc = k*w
          sum3 = sum3 + cell3d(x0,y0,z0,xc,yc,zc,w,w,w,1.,.false.)
       enddo
    enddo
 enddo
 write(label,"(a,f4.2)") '2D partition of unity, w/h=',w
 call checkval(trim(label),real(sum2),1.,1.e-11)
 write(label,"(a,f4.2)") '3D partition of unity, w/h=',w
 call checkval(trim(label),real(sum3),1.,1.e-11)

end subroutine test_partition

!-----------------------------------------------------------------------
!+
!  single cells partially overlapping the kernel, compared to
!  numerical quadrature (particle at the origin, h=1)
!+
!-----------------------------------------------------------------------
subroutine test_quadrature
 integer, parameter :: ncase = 6
 ! cell centre (x,y,z) and width: particle inside, outside, at a corner,
 ! on an edge/face, and cell faces lying exactly at q=0.5 and q=1
 real, parameter :: cells(4,ncase) = reshape((/0.4,0.3,0.2,1.0, &
                                               1.2,0.5,0.3,0.8, &
                                               0.5,0.5,0.5,1.0, &
                                               0.5,0.1,0.2,1.0, &
                                               0.75,0.2,-0.1,0.5, &
                                               1.25,0.6,0.4,0.5/),(/4,ncase/))
 real :: xc,yc,zc,w
 integer :: i
 character(len=40) :: label

 do i=1,ncase
    xc = cells(1,i)
    yc = cells(2,i)
    zc = cells(3,i)
    w  = cells(4,i)
    write(label,"(a,i1)") '2D vs quadrature, case ',i
    call checkval(trim(label),cell2d(0.,0.,xc,yc,w,w,1.,.false.),quad2d(xc,yc,w,w,400),1.e-9)
    write(label,"(a,i1)") '3D vs quadrature, case ',i
    call checkval(trim(label),cell3d(0.,0.,0.,xc,yc,zc,w,w,w,1.,.false.),quad3d(xc,yc,zc,w,60),1.e-8)
 enddo

end subroutine test_quadrature

!-----------------------------------------------------------------------
!+
!  degenerate geometries should give finite, sensible answers
!+
!-----------------------------------------------------------------------
subroutine test_edge_cases
 real :: v

 call checkval('2D edge through particle  ',pint(0.,0.5,0.5,1.),0.,0.)
 call checkval('2D zero length edge       ',pint(0.5,0.,0.,1.),0.,0.)
 call checkval('3D wall through particle  ',wallint(0.,0.,0.,0.,0.,1.,1.,1.),0.,0.)
 call checkval('2D h << pixel             ',cell2d(0.,0.,0.,0.,1.,1.,1.e-8,.false.),1.,1.e-12)
 call checkval('3D h << pixel             ',cell3d(0.,0.,0.,0.,0.,0.,1.,1.,1.,1.e-8,.false.),1.,1.e-12)
 call checkval('2D h >> pixel (scaled)    ',1.e16*cell2d(0.,0.,0.,0.,1.,1.,1.e8,.false.),&
                                            cnormk2D*wfunc(0.),1.e-10)
 v = cell3d(0.,0.,0.,0.5,0.5,0.5,1.,1.,1.,1.,.false.)
 call checkval('3D particle at corner     ',v,quad3d(0.5,0.5,0.5,1.,60),1.e-8)

end subroutine test_edge_cases

!-----------------------------------------------------------------------
!+
!  generic and legacy cubic implementations should agree
!+
!-----------------------------------------------------------------------
subroutine test_legacy
 integer, parameter :: n = 2000
 real :: r(7),err2,err3,x0,y0,z0,w,h
 integer :: i

 err2 = 0.
 err3 = 0.
 do i=1,n
    call random_number(r)
    x0 = 2.*r(1) - 1.
    y0 = 2.*r(2) - 1.
    z0 = 2.*r(3) - 1.
    w  = 0.05 + 2.*r(4)
    h  = 0.1 + 1.5*r(5)
    err2 = max(err2,abs(cell2d(x0,y0,0.,0.,w,w,h,.false.) - cell2d(x0,y0,0.,0.,w,w,h,.true.)))
    err3 = max(err3,abs(cell3d(x0,y0,z0,0.,0.,0.,w,w,w,h,.false.) - cell3d(x0,y0,z0,0.,0.,0.,w,w,w,h,.true.)))
 enddo
 call checkval('2D generic vs legacy cubic',err2,0.,1.e-10)
 call checkval('3D generic vs legacy cubic',err3,0.,1.e-10)

end subroutine test_legacy

!-----------------------------------------------------------------------
!+
!  time the cell integrals for random cells overlapping the kernel
!+
!-----------------------------------------------------------------------
subroutine time_cells(legacy,t2,t3)
 logical, intent(in)  :: legacy
 real,    intent(out) :: t2,t3
 integer, parameter :: n = 200000
 real, allocatable :: xp(:,:)
 real :: t1,t0,sum2,sum3
 integer :: i

 allocate(xp(3,n))
 call random_number(xp)
 xp = 2.*radkernel*xp - radkernel
 sum2 = 0.
 sum3 = 0.
 call cpu_time(t0)
 do i=1,n
    sum2 = sum2 + cell2d(xp(1,i),xp(2,i),0.,0.,0.5,0.5,1.,legacy)
 enddo
 call cpu_time(t1)
 t2 = (t1 - t0)/n*1.e9
 call cpu_time(t0)
 do i=1,n/10
    sum3 = sum3 + cell3d(xp(1,i),xp(2,i),xp(3,i),0.,0.,0.,0.5,0.5,0.5,1.,legacy)
 enddo
 call cpu_time(t1)
 t3 = (t1 - t0)/(n/10)*1.e9
 if (sum2 < 0. .or. sum3 < 0.) print*,' negative sums in timing test'

end subroutine time_cells

!-----------------------------------------------------------------------
!+
!  projected mass of a single unit-mass particle at (x0,y0) with
!  smoothing length h, on a grid of unit pixels
!+
!-----------------------------------------------------------------------
real function projmass(x0,y0,h,exact)
 real,    intent(in) :: x0,y0,h
 logical, intent(in) :: exact
 real, allocatable :: pix(:,:)
 real :: x(1),y(1),z(1),hh(1),w(1),dat(1),xmin
 integer :: itype(1),n

 n = 2*int(radkernel*h + 4)
 allocate(pix(n,n))
 x = x0; y = y0; z = 0.; hh = h; dat = 1.; itype = 1
 w = 1./h**3   ! so that weight*dat*h^3 = m = 1
 xmin = -0.5*n
 call interpolate3D_projection(x,y,z,hh,w,dat,itype,1,xmin,xmin,pix,n,n,1.,1., &
      .false.,0.,0.,.false.,exact,-10)
 projmass = real(sum(real(pix,kind=dp)))

end function projmass

!-----------------------------------------------------------------------
!+
!  projections should conserve mass: exactly (to round-off) for
!  exact rendering when h < 2 pixels, and on average over particle
!  positions for the standard method when h is below the pixel scale
!+
!-----------------------------------------------------------------------
subroutine test_projection_mass
 integer, parameter :: nh = 6, noff = 8
 real, parameter :: hvals(nh) = (/0.1,0.3,0.45,0.8,1.5,1.9/)
 real :: x0,y0,m,errmax,mmean
 integer :: i,j,k
 character(len=40) :: label

 call setup_integratedkernel()
 do i=1,nh
    errmax = 0.
    mmean = 0.
    ! average over a regular grid of particle offsets within a pixel
    do k=1,noff
       y0 = (k - 0.5)/noff - 0.5
       do j=1,noff
          x0 = (j - 0.5)/noff - 0.5
          m = projmass(x0,y0,hvals(i),.true.)
          errmax = max(errmax,abs(m - 1.))
          mmean = mmean + projmass(x0,y0,hvals(i),.false.)/noff**2
       enddo
    enddo
    write(label,"(a,f4.2)") 'exact projected mass, h/pw=',hvals(i)
    call checkval(trim(label),errmax,0.,1.e-10)
    write(label,"(a,f4.2)") 'mean projected mass,  h/pw=',hvals(i)
    call checkval(trim(label),mmean,1.,0.02)
 enddo

end subroutine test_projection_mass

!-----------------------------------------------------------------------
!+
!  place noff*noff well-separated particles on a grid of unit pixels,
!  each at a different sub-pixel offset (stratified over the pixel)
!+
!-----------------------------------------------------------------------
subroutine place_particles(noff,h,x,y,n)
 integer, intent(in)  :: noff
 real,    intent(in)  :: h
 real,    intent(out) :: x(noff*noff),y(noff*noff)
 integer, intent(out) :: n
 integer :: a,b,k,ns

 ns = 2*int(radkernel*h + 3)   ! even spacing in pixels, wider than the kernel
 n = noff*ns
 k = 0
 do b=1,noff
    do a=1,noff
       k = k + 1
       x(k) = (a - 0.5)*ns + (a - 0.5)/noff - 0.5
       y(k) = (b - 0.5)*ns + (b - 0.5)/noff - 0.5
    enddo
 enddo

end subroutine place_particles

!-----------------------------------------------------------------------
!+
!  geometry projection (cartesian) and opacity rendering should also
!  conserve mass and optical depth respectively
!+
!-----------------------------------------------------------------------
subroutine test_geom_opacity_mass
 integer, parameter :: nh = 6, noff = 8, np = noff*noff
 real, parameter :: hvals(nh) = (/0.1,0.3,0.45,0.8,1.5,1.9/)
 real :: x(np),y(np),z(np),hh(np),w(np),dat(np),rkappa(np),pmass(1)
 real, allocatable :: pix(:,:),tau(:,:)
 integer :: itype(np),i,n
 logical :: normalise
 character(len=40) :: label

 z = 0.; dat = 1.; itype = 1; rkappa = 1.; pmass = 1.
 do i=1,nh
    call place_particles(noff,hvals(i),x,y,n)
    hh = hvals(i)
    w = 1./hh**3   ! so that weight*dat*h^3 = m = 1
    allocate(pix(n,n),tau(n,n))
    normalise = .false.
    call interpolate3D_proj_geom(x,y,z,hh,w,dat,itype,np,0.,0.,pix,n,n,1.,1., &
         normalise,igeom_cartesian,1,2,3,(/1,2,3/),(/0.,0.,0./))
    write(label,"(a,f4.2)") 'geom projected mass,  h/pw=',hvals(i)
    call checkval(trim(label),real(sum(real(pix,kind=dp)))/np,1.,0.02)
    call interp3D_proj_opacity(x,y,z,pmass,1,hh,w,dat,z,itype,np,0.,0.,pix,tau,n,n,1.,1., &
         0.,0.,rkappa,huge(0.),-10,.true.)
    write(label,"(a,f4.2)") 'exact optical depth,  h/pw=',hvals(i)
    call checkval(trim(label),real(sum(real(tau,kind=dp)))/np,1.,1.e-10)
    call interp3D_proj_opacity(x,y,z,pmass,1,hh,w,dat,z,itype,np,0.,0.,pix,tau,n,n,1.,1., &
         0.,0.,rkappa,huge(0.),-10,.false.)
    write(label,"(a,f4.2)") 'mean optical depth,   h/pw=',hvals(i)
    call checkval(trim(label),real(sum(real(tau,kind=dp)))/np,1.,0.02)
    write(label,"(a,f4.2)") 'all particles visible,h/pw=',hvals(i)
    call checkval(trim(label),merge(1.,0.,min_block_sum(tau,n/noff,noff) > 0.),1.,0.)
    deallocate(pix,tau)
 enddo

end subroutine test_geom_opacity_mass

!-----------------------------------------------------------------------
!+
!  minimum over the noff*noff blocks of size ns of the summed pixel
!  values (each block contains the contribution of one particle)
!+
!-----------------------------------------------------------------------
real function min_block_sum(pix,ns,noff)
 integer, intent(in) :: ns,noff
 real,    intent(in) :: pix(ns*noff,ns*noff)
 integer :: a,b

 min_block_sum = huge(0.)
 do b=1,noff
    do a=1,noff
       min_block_sum = min(min_block_sum,sum(pix((a-1)*ns+1:a*ns,(b-1)*ns+1:b*ns)))
    enddo
 enddo

end function min_block_sum

!-----------------------------------------------------------------------
!+
!  the polynomial pieces used for the exact integrals must match the
!  kernel function used elsewhere, including its normalisation
!+
!-----------------------------------------------------------------------
subroutine test_kernel_pieces
 integer, parameter :: nq = 1000
 real(dp) :: qb(maxpieces),w(0:maxdeg,maxpieces),g2(-1:maxpow,maxpieces+1)
 real(dp) :: a3(-1:maxpow,maxpieces+1),c2,c3,q,wpoly,errmax
 integer  :: npieces,nmax,i,j,k

 call get_kernel_pieces(ikernel,npieces,qb,w)
 call set_kernel_integral_coeffs(npieces,qb,w,nmax,g2,a3,c2,c3)
 errmax = 0._dp
 do i=0,nq
    q = (i*radkernel)/nq
    ! evaluate the polynomial on the piece containing q
    j = 1
    do while (j < npieces .and. q >= qb(j))
       j = j + 1
    enddo
    wpoly = 0._dp
    do k=maxdeg,0,-1
       wpoly = wpoly*q + w(k,j)
    enddo
    if (q >= qb(npieces)) wpoly = 0._dp
    errmax = max(errmax,abs(wpoly - wfunc(real(q*q))))
 enddo
 call checkval('kernel pieces match wfunc  ',real(errmax),0.,1.e-12)
 ! tolerance allows for cancellation when integrating the expanded
 ! polynomials, e.g. (1 - q/2)^8 for the Wendland C6 kernel
 call checkval('2D normalisation constant  ',real(c2),cnormk2D,1.e-11*cnormk2D)
 call checkval('3D normalisation constant  ',real(c3),cnormk3D,1.e-11*cnormk3D)

end subroutine test_kernel_pieces

end module kernel_integral_tests

!-----------------------------------------------------------------------
!+
!  Unit tests of the exact kernel-cell overlap integrals
!  for all kernels, with timing against the legacy cubic code
!+
!-----------------------------------------------------------------------
program test_kernel_integrals
 use kernels, only:select_kernel,kernelname,nkernels
 use kernel_integral_tests
 implicit none
 integer :: ik
 real :: t2,t3,t2l,t3l

 do ik=1,nkernels
    call select_kernel(ik,verbose=.false.)
    print "(/,a)",' Testing exact integrals for '//trim(kernelname(ik))
    call test_kernel_pieces
    call test_normalisation
    call test_partition(0.37)
    call test_partition(1.3)
    call test_quadrature
    call test_edge_cases
    call test_projection_mass
    call test_geom_opacity_mass
 enddo

 print "(/,a)",' Testing against legacy cubic spline implementation'
 call select_kernel(1,verbose=.false.)
 call test_legacy

 print "(/,a)",' Timing (ns per cell integral, 2D pixel / 3D cell):'
 call time_cells(.true.,t2l,t3l)
 print "(3x,a24,2f10.1)",'legacy cubic',t2l,t3l
 do ik=1,nkernels
    call select_kernel(ik,verbose=.false.)
    call time_cells(.false.,t2,t3)
    print "(3x,a24,2f10.1)",kernelname(ik),t2,t3
 enddo

 print "(/,1x,i4,a,i4,a)",ntests-nfail,' of ',ntests,' tests passed'
 if (nfail > 0) stop 1

end program test_kernel_integrals
