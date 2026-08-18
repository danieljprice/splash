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
! module for 3D iron-filings vector field visualisation
!
module fieldlines
 implicit none
 public :: vecplot3D_proj,interpolate_pt

contains

!
!--interpolate from particles to single point
!  would be nice to know neighbours
!
subroutine interpolate_pt(xpt,ypt,vxpt,vypt,x,y,vecx,vecy,h,pmass,rho,npart)
 use kernels, only:radkernel2,wfunc,cnormk2D
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
    if (q2 < radkernel2) then
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
subroutine vecplot3D_proj(x,y,z,vx,vy,vz,vecmax,weight,itype,n,dx,zobs,dscreen, &
                          use_zslice,zslicemin,zslicemax)
 use plotlib, only:plot_line,plot_bbuf,plot_ebuf,plot_slw,plot_sci,plot_set_opacity
 use plotlib, only:plot_qcr,plot_scr,plot_qlw,plot_arro,plot_sah
 use sort,    only:indexx
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: x,y,z,vx,vy,vz,weight
 integer, dimension(n), intent(in) :: itype
 real, intent(inout) :: vecmax
 real, intent(in) :: dx,zobs,dscreen,zslicemin,zslicemax
 logical, intent(in) :: use_zslice
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
 if (vecmax <= 0. .or. vecmax > 0.5*huge(vecmax)) then
    vmax = 0.
    do i=1,n
       if (itype(i) >= 0 .and. weight(i) > 0.) then
          vmax = max(vx(i)**2 + vy(i)**2 + vz(i)**2,vmax)
       endif
    enddo
    vmax = sqrt(vmax)
    vecmax = vmax
 else
    vmax = vecmax
 endif

 use3Dperspective = abs(dscreen) > tiny(dscreen)

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
    if (itype(i) >= 0 .and. weight(i) > 0.) then
       if (use_zslice) then
          if (z(i) < zslicemin .or. z(i) > zslicemax) cycle over_particles
       endif
       if (use3Dperspective) then
          if (z(i) > zobs) cycle over_particles
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

       if (frac >= 1.e-3) then
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
