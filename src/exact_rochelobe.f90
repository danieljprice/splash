!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------------
! Plots Roche Lobes (equipotential surface) for a binary system
! Guts of it adapted from an original C routine by Simon Portugies-Zwart
!-----------------------------------------------------------------------
module rochelobe
 implicit none
 real, parameter :: roche_accuracy = 1.e-6
 integer, parameter :: maxits = 100

contains

!-----------------------------------------------------------------------
! Plot Roche Lobe
!
! INPUT:
!   x1, y1 : position of primary
!   x2, y2 : position of secondary
!   m1, m2 : masses of components
!
! OUTPUT:
!   xplot, yplot: contains (half of) the Roche lobe solution
!   ierr : error condition to splash indicating plotting was done here
!-----------------------------------------------------------------------
subroutine exact_rochelobe(x1,y1,x2,y2,m1,m2,xplot,yplot,ierr)
 use plotlib, only:plot_line
 real, intent(in) :: x1,y1,x2,y2,m1,m2
 real, dimension(:), intent(inout) :: xplot,yplot
 integer, intent(out) :: ierr
 real    :: roche_radius1,roche_radius2,q,L1
 real    :: rlimit,llimit,xmax,xmin,ymax,xtmp
 real    :: angle,cosangle,sinangle,dx,dy,sep
 integer :: i,npts

 npts = (size(xplot)-1)/2
 if (npts < 1) return
 
 sep = sqrt((x2 - x1)**2 + (y2 - y1)**2)
 print "(4(a,es10.3))",' plotting Roche potential, m1 = ',m1,' m2 = ',m2,' sep = ',sep

 roche_radius1 = roche_radius(m1, m2)
 roche_radius2 = roche_radius(m2, m1)
 q  = m2/m1
 L1 = first_Lagrangian_point(1./q)

 ! We assume primary on the left
 if (m1 < m2) then
    q = 1/q
    L1 = first_Lagrangian_point(1./q)
 endif

 rlimit = right_limit(q, L1)
 llimit = left_limit(q, L1)

 call compute_lobes(q, L1, npts, xplot, yplot)

 xmin = xplot(1)
 xmax = xplot(2*npts)
 if (m1 < m2) then
    xtmp = xmin
    xmin = 1.-xmax
    xmax = 1.-xtmp
 endif

 ymax = -1000.
 do i=1,2*npts
    ymax = max(ymax,yplot(i))
 enddo

 ! some tedious fiddling for q>1 case
 if (m1 < m2) then
    xplot(:) = 1. - xplot(:)
 endif
 
 ! scale to actual separation
 xplot = xplot*sep
 yplot = yplot*sep

 ! work out angle needed to rotate into corotating frame
 dx   = x2 - x1
 dy   = y2 - y1
 angle = -atan2(dy,dx)
 cosangle = cos(angle)
 sinangle = sin(angle)
 
 ! lobes are computed assuming primary is at the origin, so shift to xprim,yprim
 ! unrotated, this is just plot_line(xplot,yplot) and plot_line(xplot,-yplot)
 call plot_line(2*npts+1,xplot*cosangle + yplot*sinangle + x1,-xplot*sinangle + yplot*cosangle + y1)
 call plot_line(2*npts+1,xplot*cosangle - yplot*sinangle + x1,-xplot*sinangle - yplot*cosangle + y1)

 !--return non-zero ierr value as we do the plotting here
 ierr = 1

end subroutine exact_rochelobe
!
! calculates outer limit of roche-lobe
!
subroutine rlimit(q, L, x, f, df, dummy)
 real, intent(in) :: q,L,x,dummy
 real, intent(out) :: f,df
 real :: qi,q11,cnst,r1,r2,r3

 qi = 1./q
 q11 = 1./(1.+qi)
 cnst = qi/L+1./(1.-L)+0.5*(1.+qi)*(L-q11)**2

 r1 = abs(x)
 r2 = abs(1-x)
 r3 = abs(x-q11)

 f = qi/r1+1./r2+0.5*(1.+qi)*r3**2 - cnst
 df = -qi*x/r1**3 +(1.-x)/r2**3 + (1.+qi)*(x-q11)

end subroutine rlimit

real function rtsafe(func,q,L,x1,x2,xll,xacc)
 real, intent(in) :: q,L,x1,x2,xll,xacc
 external :: func
 integer :: j
 real :: df,dx,dxold,f,fh,fl
 real :: temp,xh,xl,rts

 call func(q,L,x1,fl,df,xll)
 call func(q,L,x2,fh,df,xll)

 if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
   !print*,'Error occured in rtsafe, exiting...',q,L,x1,x2,fl,fh
   rtsafe = 0.
   return
 endif

 if (abs(fl) < tiny(fl)) then
    rtsafe = x1
    return
 endif

 if (abs(fh) < tiny(fh)) then
    rtsafe = x2
    return
 endif

 call func(q, L, x1, f, df, xll)
 call func(q, L, x2, f, df, xll)

 if (fl < 0.0) then
    xl=x1
    xh=x2
 else
    xh=x1
    xl=x2
 endif

 rts = 0.5*(x1+x2)
 dxold = abs(x2-x1)
 dx = dxold
 call func(q, L, rts, f, df, xll)
 do j=1,maxits
    if ((((rts-xh)*df - f)*((rts-xl)*df - f) >= 0.0) &
       .or. (abs(2.0*f) > abs(dxold*df))) then
       dxold = dx
       dx = 0.5*(xh-xl)
       rts = xl+dx
       if (abs(xl-rts) < tiny(rts)) then
          rtsafe = rts
          return
       endif
    else
       dxold = dx
       dx = f/df
       temp = rts
       rts = rts - dx
       if (abs(temp-rts) < tiny(rts)) then
          rtsafe = rts
          return
       endif
    endif
    if (abs(dx) < xacc) then
       rtsafe = rts
       return
    endif

    call func(q, L, rts, f, df, xll)
    if (f < 0.0) then
       xl = rts
    else
       xh = rts
    endif
 enddo
 rtsafe = 0.
 return
end function rtsafe

real function left_limit(q, L)
 real, intent(in) :: q,L

 left_limit = rtsafe(rlimit,q,L,-0.5*L,-L,0., roche_accuracy);

end function left_limit

real function right_limit(q, L)
 real, intent(in) :: q,L

 right_limit = rtsafe(rlimit,q,L,1.5-0.5*L,2.0-L,0., roche_accuracy);

end function right_limit

!
! return roche radius as fraction of the semi-major axis.
! So to obtain the true roche_radius call:
! real Rl = semi_major_axis * roche_radius(m1, m2);
! Eggleton PP., ApJ, 1983, 268, 368.
!
real function roche_radius(mthis, mother)
 real, intent(in) :: mthis,mother
 real :: mr,q1_3,q2_3

 mr = mthis/mother
 q1_3 = mr**(1./3.)
 q2_3 = q1_3**2

 roche_radius = 0.49*q2_3/(0.6*q2_3 + log(1 + q1_3))

end function roche_radius

real function first_Lagrangian_point(qinv)
 real, intent(in) :: qinv
 real :: fL, dfL, dL, L, q11

 q11 = 1./(1.+qinv)
 L = 0.5 + 0.2222222*log10(qinv)

 dL = 1.e7
 do while (abs(dL)>1.e-6)
   fL = qinv/L**2- 1./(1.-L)**2 - (1.+qinv)*L + 1.
   dfL=-2*qinv/L**3 - 2./(1.-L)**3 - (1.+qinv)
   dL = -fL/(dfL*L)
   L = L*(1.+dL)
 enddo

 first_Lagrangian_point = L

end function first_Lagrangian_point

subroutine rline(q, L, y, f, df, xl)
 real, intent(in)  :: q, L, y, xl
 real, intent(out) :: f, df
 real :: xsq,onexsq,qi,q11,cnst,cnst2,r1,r2

 xsq=xl*xl
 onexsq=(1.-xl)**2

 qi=q
 q11=1./(1.+qi)
 cnst =qi/L+1./(1.-L) + 0.5*(1.+qi)*(L-q11)**2
 cnst2=0.5*(1.+qi)*(xl-q11)**2 - cnst

 r1=sqrt(xsq+y)
 r2=sqrt(onexsq+y)
 f=qi/r1+1./r2+cnst2
 df =-0.5*qi/r1**3 - 0.5/r2**3

end subroutine rline

subroutine compute_lobes(q, L, npts, xplot, yplot)
 real, intent(in) :: q, L
 integer, intent(in) :: npts
 real, intent(out), dimension(2*npts+1) :: xplot,yplot
 real :: qi,q11,cnst,lrl,rrl,y1,y2,ysq,dxl,dxr
 integer :: i

 qi = 1/q
 q11 = 1./(1.+qi)
 cnst = qi/L+1./(1.-L) + 0.5*(1.+qi)*(L-q11)**2

 lrl = left_limit(q, L)
 xplot(1) = lrl
 yplot(1) = 0.

 xplot(npts+1) = L
 yplot(npts+1) = 0.

 rrl = right_limit(q, L)
 xplot(2*npts) = rrl
 yplot(2*npts) = 0.

 y1 = 0.
 y2 = L*L

 !--left lobe
 dxl = (xplot(npts+1)-xplot(1))/real(npts)
 do i=1,npts
    xplot(i+1) = xplot(1) + i*dxl
    ysq = rtsafe(rline,qi,L,y1,y2,xplot(i+1),roche_accuracy)
    yplot(i+1) = sqrt(ysq)
 enddo

 !--right lobe
 dxr = (xplot(2*npts)-xplot(npts+1))/real(npts)
 do i=1,npts
    xplot(npts+i+1) = xplot(npts+1) + i*dxr
    ysq = rtsafe(rline,qi,L,y1,y2,xplot(npts+i+1),roche_accuracy)
    yplot(npts+i+1) = sqrt(ysq)
 enddo

end subroutine compute_lobes

end module rochelobe
