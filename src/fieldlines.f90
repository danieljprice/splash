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
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!
! module for field line / flux tube plotting in 2 and 3 dimensions
!
module fieldlines
 implicit none
 public :: streamlines,vecplot3D_proj
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
!  A(x,y) = \int v_x(x,y) dy - \int v_y(x,y) dx
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
! Oct 2007: uses Simpson's rule instead of trapezoidal
!
!---------------------------------------------------------------------------
subroutine streamlines(vecpixx,vecpixy,datpix,npixx,npixy,pixwidth)
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in), dimension(npixx,npixy) :: vecpixx,vecpixy
 real, intent(in) :: pixwidth
 real, intent(out), dimension(npixx,npixy) :: datpix
 real, dimension(npixx,npixy) :: datpix2
 double precision :: fyj,fyjhalf,term,termi,termj,fxprevi,fyjprev
 double precision, dimension(npixx) :: fx,fxhalf
 integer :: i,j
 !
 !--check for errors in input
 !
 if (pixwidth.le.0.) then
    print "(1x,a)",'streamlines: error: pixel width <= 0'
    datpix = 0.
    return
 endif

 fyj = 0.
 fyjhalf = 0.
 !
 !--perform the integration forwards
 !
 do j=1,npixy
    do i=1,npixx
       term = 0.
       if (i.eq.1) then
          fyj = 0.
          fyjhalf = 0.
       else
          fyjprev = fyj
          !--trapezoidal rule in x
          termj = 0.5*pixwidth*(vecpixy(i-1,j)+vecpixy(i,j))
          fyj = fyj - termj
          if (mod(i-1,2).eq.0) then ! 3, 5, 7, 9 ...
             !
             !--for odd points, use trapezoidal solution at half grid points
             !  to get Simpson's rule
             !
             fyjhalf = fyjhalf - pixwidth*(vecpixy(i-2,j)+vecpixy(i,j))
             term = term + 4./3.*fyj - 1./3.*fyjhalf
          else
             !
             !--for even points, use Simpson's rule up to last odd point
             !  then finish with a trapezoidal integration over last two points
             !
             term = term + 4./3.*fyjprev - 1./3.*fyjhalf - termj
          endif
       endif
       !
       !--same as above but for integration in y
       !
       if (j.eq.1) then
          fx(i) = 0.
          fxhalf(i) = 0.
          fxprevi = 0.
       else
          fxprevi = fx(i)
          termi = 0.5*pixwidth*(vecpixx(i,j-1)+vecpixx(i,j))
          fx(i) = fx(i) + termi
          if (mod(j-1,2).eq.0) then
             fxhalf(i) = fxhalf(i) + pixwidth*(vecpixx(i,j-2)+vecpixx(i,j))
             term = term + 4./3.*fx(i) - 1./3.*fxhalf(i)
          else
             term = term + 4./3.*fxprevi - 1./3.*fxhalf(i) + termi
          endif
       endif

       datpix(i,j) = real(term)
    enddo
 enddo
 !
 !--perform the integration backwards
 !
 datpix2 = 0.
 do j=npixy,1,-1
    do i=npixx,1,-1
       term = 0.
       if (i.eq.npixx) then
          fyj = 0.
          fyjhalf = 0.
       else
          fyjprev = fyj
          !--trapezoidal rule in x
          termj = 0.5*pixwidth*(vecpixy(i+1,j)+vecpixy(i,j))
          fyj = fyj + termj
          if (mod(npixx-i,2).eq.0) then ! 3, 5, 7, 9 ...
             !
             !--for odd points, use trapezoidal solution at half grid points
             !  to get Simpson's rule
             !
             fyjhalf = fyjhalf + pixwidth*(vecpixy(i+2,j)+vecpixy(i,j))
             term = term + 4./3.*fyj - 1./3.*fyjhalf
          else
             !
             !--for even points, use Simpson's rule up to last odd point
             !  then finish with a trapezoidal integration over last two points
             !
             term = term + 4./3.*fyjprev - 1./3.*fyjhalf + termj
          endif
       endif
       !
       !--same as above but for integration in y
       !
       if (j.eq.npixy) then
          fx(i) = 0.
          fxhalf(i) = 0.
          fxprevi = 0.
       else
          fxprevi = fx(i)
          termi = 0.5*pixwidth*(vecpixx(i,j+1)+vecpixx(i,j))
          fx(i) = fx(i) - termi
          if (mod(npixy-j,2).eq.0) then
             fxhalf(i) = fxhalf(i) - pixwidth*(vecpixx(i,j+2)+vecpixx(i,j))
             term = term + 4./3.*fx(i) - 1./3.*fxhalf(i)
          else
             term = term + 4./3.*fxprevi - 1./3.*fxhalf(i) - termi
          endif
       endif

       datpix2(i,j) = real(term)
    enddo
 enddo

 !
 !--average the two
 !
 datpix = 0.5*(datpix + datpix2)

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
 use plotlib, only:plot_line
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
    call plot_line(2,xline,yline)

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
 use kernels, only:radkernel2,wfunc,cnormk2D
 implicit none
 integer, intent(in) :: npart
 real, dimension(npart), intent(in) :: x,y,vecx,vecy,h,pmass,rho
 real, intent(in) :: xpt, ypt
 real, intent(out) :: vxpt, vypt
 real :: rho1i,term,const,dx,dy,hi1,q2,wab
 integer :: i

 vxpt = 0.
 vypt = 0.
 const = cnormk2D

 do i=1,npart
    dx  = xpt - x(i)
    dy  = ypt - y(i)
    hi1 = 1./h(i)
    q2  = (dx*dx + dy*dy)*hi1*hi1
    !
    !--if particles are within range, calculate contribution to this pt
    !
    if (q2.lt.radkernel2) then
       if (rho(i) > 0.) then
          rho1i = 1./rho(i)
       else
          rho1i = 0.
       endif
       term = const*pmass(i)*rho1i
       wab = wfunc(q2)

       vxpt = vxpt + term*vecx(i)*wab
       vypt = vypt + term*vecy(i)*wab
    endif
 enddo

end subroutine interpolate_pt

!--------------------------------------------------------------------------
!   Visualisation of a 3D vector field in projection
!   by means of "iron filings" drawn on particles
!
!   We draw a line on each particle in the direction of the vector field
!   This line is illuminated by reflections from a lighting source, and
!   is drawn with an opacity proportional to the field strength, such that
!   strong field regions are highlighted.
!
!   For details of the lighting algorithm, see e.g.
!   Stalling, Zoeckler and Hege, 1997, IEEE Trans. Viz. Comp. Graphics, 3, 118-128
!   
!   Added by D. Price, Dec 2011
!--------------------------------------------------------------------------
subroutine vecplot3D_proj(x,y,z,vx,vy,vz,vecmax,weight,itype,n,dx,zobs,dscreen)
 use plotlib, only:plot_line,plot_bbuf,plot_ebuf,plot_slw,plot_sci,plot_set_opacity
 use plotlib, only:plot_qcr,plot_scr,plot_qlw,plot_arro,plot_sah
 use sort,    only:indexx
 implicit none
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: x,y,z,vx,vy,vz,weight
 integer, dimension(n), intent(in) :: itype
 real, intent(inout) :: vecmax
 real, intent(in) :: dx,zobs,dscreen
 integer, dimension(n) :: iorder
 integer :: i,ipart,np
 real, dimension(2) :: xpts,ypts
 real :: vxi,vyi,vzi,dvmag,zfrac,vmax,vmag,frac,ri,gi,bi,term,lw
 real :: toti,fambient,diffuse,specular,fdiff,fspec,ldotn,vdotr,ldott,vdott
 integer :: pdiff,nspec,lwold
 real, dimension(3) :: vunit,lighting,viewangle
 logical :: white_bg,use3Dperspective

 !
 !--get the max adaptively if it is not already set
 !
 if (vecmax.le.0. .or. vecmax.gt.0.5*huge(vecmax)) then
    vmax = 0.
    do i=1,n
       if (itype(i).ge.0 .and. weight(i).gt.0.) then
          vmax = max(vx(i)**2 + vy(i)**2 + vz(i)**2,vmax)
       endif
    enddo
    vmax = sqrt(vmax)
    vecmax = vmax
 else
    vmax = vecmax
 endif
 
 use3Dperspective = abs(dscreen).gt.tiny(dscreen)
 
 !
 !--work out whether or not we have a white or black
 !  background colour
 !
 !call plot_sah(1,20.,1.0)
 call plot_qcr(0,ri,gi,bi)
 white_bg = (ri + gi + bi > 1.5)
 !
 !--specify the parameters in the lighting algorithm
 !  these should differ depending on whether we are drawing
 !  on a white or black background
 ! 
 if (white_bg) then
    fambient = 0.
    fdiff = 0.1
    fspec = 0.5
 else
    fambient = 0.3
    fdiff = 0.7
    fspec = 0.8
 endif
 pdiff = 4
 nspec = 12
 !
 !--specify the viewing and lighting angles
 ! 
 viewangle = (/0.,0.,1./)
 !lighting = (/0.3,0.3,1./)
 lighting = (/0.,0.,1./)
 
 !--make sure these are normalised
 !lighting = lighting/sqrt(dot_product(lighting,lighting))

 print*,'plotting 3D field structure: min,max = ',1.e-3*vmax,vmax
!
!--first sort the particles in z so that we do the opacity in the correct order
!
 call indexx(n,z,iorder)
 call plot_bbuf
 call plot_qcr(1,ri,gi,bi)

 np = 0
 zfrac = 1.
 call plot_qlw(lwold)
 lw = 2.*lwold
 over_particles: do ipart=1,n
    i = iorder(ipart)
    if (itype(i).ge.0 .and. weight(i).gt.0.) then
       if (use3Dperspective) then
          if (z(i).gt.zobs) cycle over_particles
          zfrac = abs(dscreen/(z(i)-zobs))
       endif
       !if (mod(ipart,10)/=0) cycle over_particles
!       lw = min(zfrac,2.5)
       
       vxi = vx(i)
       vyi = vy(i)
       vzi = vz(i)
       !
       !--we draw lines on each particle with an 
       !  opacity proportional to the field strength
       !
       vmag = sqrt(vxi**2 + vyi**2 + vzi**2)
       dvmag = 1./vmag
       vunit = abs((/vxi,vyi,vzi/)*dvmag)
       frac = min(vmag/vmax,1.0)

       if (frac.ge.1.e-3) then
          !--specify the length of line to draw
          term = 1.5*dx*dvmag*zfrac
          !term = term*(vmag/vmax)**0.2
          xpts(1) = x(i) - vxi*term
          xpts(2) = x(i) + vxi*term
          ypts(1) = y(i) - vyi*term
          ypts(2) = y(i) + vyi*term
          
          !--draw "halo" in background colour with
          !  twice the thickness, same opacity
          call plot_slw(2.*lw)
          call plot_sci(0)
          call plot_set_opacity(frac)
          call plot_line(2,xpts,ypts)

          !--Phong lighting
          ldott = dot_product(lighting,vunit)
          ldotn = sqrt(1. - ldott**2)
          diffuse = fdiff*(ldotn)**pdiff

          vdott = dot_product(viewangle,vunit)
          vdotr = ldotn*sqrt(1. - vdott**2) - ldott*vdott
          specular = fspec*(vdotr)**nspec
          toti = (fambient + diffuse + specular)

          !--draw line with intensity proportional
          !  to the amount of lighting
          !call plot_scr(1,toti,toti,toti,max(frac,0.15))
          call plot_scr(1,toti,toti,toti,max(frac,0.05))
          call plot_sci(1)
          call plot_slw(lw)
          call plot_line(2,xpts,ypts)
          !call plot_arro(xpts(1),ypts(1),xpts(2),ypts(2))
          np = np + 1
       endif
    endif
 enddo over_particles
 !--reset opacity for both foreground and background colour indices
 call plot_sci(0)
 call plot_set_opacity(1.0)
 call plot_sci(1)
 call plot_scr(1,ri,gi,bi)
 call plot_set_opacity(1.0)
 call plot_ebuf
 call plot_slw(lwold)
 print*,' plotted ',np,' of ',n,' particles'

end subroutine vecplot3D_proj

end module fieldlines
