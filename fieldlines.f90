!
! module for field line / flux tube plotting in 2 and 3 dimensions
!
module fieldlines
 implicit none
 public :: streamlines
 private :: trace2D,interpolate_pt

 private
 
contains

!---------------------------------------------------------------------------
!
! This subroutine integrates a 2D vector field to give the stream function
! Plotting the contours of this function gives the field/stream lines
!
! The solution is given by
!
!  A(x,y) = \int v_x dy - \int v_y dx
!
! which we compute, knowing v_x and v_y on a fixed grid of square pixels,
! by a simple trapezoidal integration for each component.
!
! For SPH, this means we first interpolate the vector field to
! get v_x and v_y on the two-dimensional grid. This then also works for
! cross sections / projections of 3D vector fields by first interpolating
! to the 2D grid.
!
!
! Inputs: vecpixx(npixx,npixy) : x component of vector field on fixed grid
!         vecpixy(npixx,npixy) : y component of vector field on fixed grid
!         xmin, ymin           : xmin and ymin of grid
!         pixwidth             : grid cell size (pixels are square)
!         npixx,npixy          : number of grid cells (pixels) in x,y
!
! Output: datpix(npixx,npixy) : stream function on fixed grid
!
! written by Daniel Price dprice@astro.ex.ac.uk
!
! Last modified: 5th Oct. 2006
!
!---------------------------------------------------------------------------
subroutine streamlines(vecpixx,vecpixy,datpix,npixx,npixy,xmin,ymin,pixwidth)
 implicit none
 real, intent(in), dimension(npixx,npixy) :: vecpixx,vecpixy
 real, intent(in) :: xmin,ymin,pixwidth
 real, intent(out), dimension(npixx,npixy) :: datpix
 integer, intent(in) :: npixx,npixy
 real :: weighty,weightx,termx,termy,fyi
 real, dimension(npixx) :: fx
 integer :: i,j
 !
 !--check for errors in input
 !
 if (pixwidth.le.0.) then
    print "(1x,a)",'streamlines: error: pixel width <= 0'
    datpix = 0.
    return
 endif
 !
 !--initialise quantities
 !
 do i=1,npixx
    fx(i) = 0.
 enddo
 !
 !--perform the integration
 !
 do j=1,npixy
    fyi = 0.
    do i=1,npixx
       termx = vecpixx(i,j)*pixwidth
       termy = -vecpixy(i,j)*pixwidth
       if (i.eq.1) then ! trapezoidal rule
          fyi = fyi + 0.5*termy
          weighty = -0.5
       else
          weighty = 0.5
          fyi = fyi + termy
       endif
       if (j.eq.1 .or. j.eq.npixy) then ! trapezoidal rule
          fx(i) = fx(i) + 0.5*termx
          weightx = -0.5
       else
          fx(i) = fx(i) + termx
          weightx = 0.5
       endif
       datpix(i,j) = fx(i) + fyi - weighty*termy - weightx*termx
    enddo
 enddo
 
 return
end subroutine streamlines

!---------------------------------------------------------------------------------
!
! THE REST OF THIS MODULE IS EITHER OLD, *VERY* EXPERIMENTAL AND/OR UNDOCUMENTED
!
!---------------------------------------------------------------------------------

!
! we want to trace the curve through a 2D vector field
!
subroutine fieldlines2D(npart,x,y,vecx,vecy,h,pmass,rho,xmin,xmax,ymin,ymax)
 implicit none
 integer, intent(in) :: npart
 real, intent(in), dimension(npart) :: x,y,vecx,vecy,h,pmass,rho
 real, intent(in) :: xmin,xmax,ymin,ymax
 integer :: i,nlines
 real :: dx,dy,xstart,ystart,ymaxline

 nlines = 10
 dx = (xmax-xmin)/nlines
 dy = (ymax-ymin)/nlines
 ystart = ymin
 ymaxline = ymin

 do while (ymaxline.lt.ymax)
    do i=1,nlines
       xstart = xmin + (i-1)*dx + 0.5*dx
       print*,' tracing field line ',i,' x, y = ',xstart,ystart
       call trace2D(xstart,ystart,xmin,xmax,ymin,ymax,ymaxline, &
                 x,y,vecx,vecy,h,pmass,rho,npart)
    enddo
    ystart = ymaxline + 0.5*dy
 enddo

end subroutine fieldlines2D

subroutine trace2D(xstart,ystart,xmin,xmax,ymin,ymax,ymaxline, &
                   x,y,vecx,vecy,h,pmass,rho,npart)
 implicit none
 integer, intent(in) :: npart
 real, intent(in) :: xstart,ystart,xmin,xmax,ymin,ymax
 real, intent(inout) :: ymaxline
 real, dimension(npart), intent(in) :: x,y,vecx,vecy,h,pmass,rho
 integer :: ipt,npix
 real, dimension(2) :: xline, yline
 real :: runit,runit1,dx,dy,vx,vy,pixwidth,sign

 xline(1) = xstart
 yline(1) = ystart
 npix = 100
 pixwidth = (xmax - xmin)/npix
 ipt = 0
 sign = 1.0
  
 do while ((xline(1).ge.xmin .and. xline(1).le.xmax) .and. &
           (yline(1).ge.ymin .and. yline(1).le.ymax) .and. ipt.lt.10*npix)
    
    ipt = ipt + 1
    !
    !--get dx and dy from interpolation of vector field from particles
    !
    call interpolate_pt(xline(1),yline(1),vx,vy, &
                        x,y,vecx,vecy,h,pmass,rho,npart)
    !
    !--get unit vector in direction of vector
    !
    runit = sqrt(vx**2 + vy**2)
    if (runit.gt.0) then
       runit1 = 1./runit
       dx = vx*runit1*pixwidth
       dy = vy*runit1*pixwidth
    else
       dx = 0.
       dy = 0.
    endif
    
    if (ipt.eq.1 .and. dy.lt.0.) sign = -1.0
    
    xline(2) = xline(1) + sign*dx
    yline(2) = yline(1) + sign*dy

    !!print*,'x, y = ',xline(2),yline(2)
    ymaxline = max(ymaxline,yline(2))
    !
    !--plot line segment
    !
    call pgline(2,xline,yline)
    
    xline(1) = xline(2)
    yline(1) = yline(2)
    
 enddo
 
 if (ipt.ge.10*npix) print*,'WARNING: infinite field line'

end subroutine trace2D

!
!--interpolate from particles to single point
!  would be nice to know neighbours
!
subroutine interpolate_pt(xpt,ypt,vxpt,vypt,x,y,vecx,vecy,h,pmass,rho,npart)
 implicit none
 integer, intent(in) :: npart
 real, dimension(npart), intent(in) :: x,y,vecx,vecy,h,pmass,rho
 real, intent(in) :: xpt, ypt
 real, intent(out) :: vxpt, vypt
 real, parameter :: pi = 3.141592653589
 real :: rho1i,term,const,dx,dy,rr,qq,wab
 integer :: i

 vxpt = 0.
 vypt = 0.

 do i=1,npart
    dx = xpt - x(i)
    dy = ypt - y(i)
    rr = sqrt(dx**2 + dy**2)
    qq = rr/h(i)
    !
    !--if particles are within range, calculate contribution to this pt
    !
    if (qq.lt.2.0) then
       const = 10./(7.*pi*h(i)**2)
       if (rho(i).ne.0.) then
          rho1i = 1./rho(i)
       else
          rho1i = 0.
       endif
       term = const*pmass(i)*rho1i
       if (qq.lt.1.0) then
          wab = (1.-1.5*qq**2 + 0.75*qq**3)
       else
          wab = 0.25*(2.-qq)**3
       endif
       
       vxpt = vxpt + term*vecx(i)*wab
       vypt = vypt + term*vecy(i)*wab
    endif
 enddo

end subroutine interpolate_pt

end module fieldlines
