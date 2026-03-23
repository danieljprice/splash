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
!  Copyright (C) 2005-2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module interactive_utils
 implicit none
 public :: inslice,inrectangle,incircle,inpoly
 public :: adapt_limits_interactive,get_vptxy,get_posxy,get_button_anchor
 public :: get_pix2viewport_factor, get_xw_margin_bounds
 ! from GIZA_XW_MARGIN in giza-driver-xw.c
 integer, parameter, private :: giza_margin_xw = 20

contains

!--------------------------------------------------------------------
! utilities to determine whether a point is in or out of a selection
!--------------------------------------------------------------------
logical function inslice(x,xmin,xmax)
 real, intent(in) :: x,xmin,xmax

 inslice = (x >= xmin .and. x <= xmax)

end function inslice

logical function inrectangle(x,y,xmin,xmax,ymin,ymax)
 real, intent(in) :: x,y,xmin,xmax,ymin,ymax

 inrectangle = (x >= xmin .and. x <= xmax .and. y >= ymin .and. y <= ymax)

end function inrectangle

logical function incircle(x,y,r2)
 real, intent(in) :: x,y,r2

 incircle = ((x*x + y*y) <= r2)

end function incircle

!
! Point in polygon
! See: http://en.wikipedia.org/wiki/Even-odd_rule
!
logical function inpoly(x,y,xpts,ypts,npts)
 real, intent(in) :: x,y
 real, dimension(:), intent(in) :: xpts,ypts
 integer, intent(in) :: npts
 integer :: i,j

 inpoly = .false.
 j = npts
 do i=1,npts
    if (((ypts(i) > y) .neqv. (ypts(j) > y)) .and. &
        (x < (xpts(j) - xpts(i))*(y-ypts(i))/(ypts(j) - ypts(i)) + xpts(i))) then
       inpoly = .not. inpoly
    endif
    j = i
 enddo

end function inpoly

!------------------------------------------------------------
! utility which adapts plot limits based only on the
! particles being plotted
!------------------------------------------------------------
subroutine adapt_limits_interactive(labeli,np,xarr,xmin,xmax,icolourpart,iamtype,iusetype)
 use params, only:int1
 use limits, only:assert_sensible_limits
 character(len=*), intent(in)    :: labeli
 integer, intent(in)             :: np
 real, dimension(np), intent(in) :: xarr
 real, intent(out)               :: xmin,xmax
 integer(kind=int1), dimension(:) , intent(in) :: iamtype
 integer,            dimension(np), intent(in) :: icolourpart
 logical,            dimension(:),  intent(in) :: iusetype
 integer :: itype,i
 logical :: mixedtypes

 xmin =  huge(xmin)
 xmax = -huge(xmax)
 mixedtypes = size(iamtype) >= np

 if (mixedtypes) then
    do i=1,np
       itype = int(iamtype(i))
       if (itype > 0 .and. itype <= np) then
          if (iusetype(itype) .and. icolourpart(i) > 0) then
             xmin = min(xmin,xarr(i))
             xmax = max(xmax,xarr(i))
          endif
       endif
    enddo
 else
    xmin = minval(xarr,mask=(icolourpart >= 0))
    xmax = maxval(xarr,mask=(icolourpart >= 0))
 endif
 call assert_sensible_limits(xmin,xmax)

 !print "(1x,a)",' resetting '//trim(labeli)//' limits'

end subroutine adapt_limits_interactive

!------------------------------------------------------------
! utility which translates between world co-ordinates (x,y)
! and viewport co-ordinates (relative to the whole viewport)
!------------------------------------------------------------
subroutine get_vptxy(x,y,vptx,vpty)
 use plotlib, only:plot_qvp,plot_qwin
 real, intent(in) :: x,y
 real, intent(out) :: vptx,vpty
 real :: xmini,xmaxi,ymini,ymaxi
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi

 call plot_qvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
 call plot_qwin(xmini,xmaxi,ymini,ymaxi)
 vptx = vptxmini + (x-xmini)/(xmaxi-xmini)*(vptxmaxi-vptxmini)
 vpty = vptymini + (y-ymini)/(ymaxi-ymini)*(vptymaxi-vptymini)

end subroutine get_vptxy

!------------------------------------------------------------
! utility to return x,y coordinates given viewport coords
! (only works for single-panelled plots)
!------------------------------------------------------------
subroutine get_posxy(vptx,vpty,x,y,xmini,xmaxi,ymini,ymaxi)
 use plotlib, only:plot_qvp
 real, intent(in) :: vptx,vpty
 real, intent(out) :: x,y
 real, intent(in) :: xmini,xmaxi,ymini,ymaxi
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi

 call plot_qvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
 x = xmini + (vptx-vptxmini)/(vptxmaxi-vptxmini)*(xmaxi-xmini)
 y = ymini + (vpty-vptymini)/(vptymaxi-vptymini)*(ymaxi-ymini)

end subroutine get_posxy

!---------------------------------------------
! get the anchor position for the button set
!---------------------------------------------
subroutine get_button_anchor(xleft_vp0,ybottom_vp0)
 real, intent(out) :: xleft_vp0,ybottom_vp0
 real :: pix_to_viewport_x,pix_to_viewport_y

 call get_pix2viewport_factor(pix_to_viewport_x,pix_to_viewport_y)
 xleft_vp0 = -18.*pix_to_viewport_x   ! start 18 pixels left of viewport edge
 ybottom_vp0 = -18.*pix_to_viewport_y ! start 18 pixels below bottom margin

end subroutine get_button_anchor

!---------------------------------------------
! get the y bounds of the X-window margin
!---------------------------------------------
subroutine get_xw_margin_bounds(y1,y2)
 real, intent(out) :: y1,y2
 real :: pix_to_viewport_x,pix_to_viewport_y

 call get_pix2viewport_factor(pix_to_viewport_x,pix_to_viewport_y)
 y1 = -20.*pix_to_viewport_y
 y2 = -1.*pix_to_viewport_y

end subroutine get_xw_margin_bounds

!---------------------------------------------
! get conversion factor from viewport coordinates to pixels
!---------------------------------------------
subroutine get_pix2viewport_factor(xfac,yfac)
 use plotlib, only:plot_qcs
 real, intent(out) :: xfac,yfac
 real :: xch,ych,xch_pix,ych_pix

 call plot_qcs(0,xch,ych)
 call plot_qcs(3,xch_pix,ych_pix)
 xfac = xch / xch_pix * giza_margin_xw / 20.
 yfac = ych / ych_pix * giza_margin_xw / 20.

end subroutine get_pix2viewport_factor

end module interactive_utils
