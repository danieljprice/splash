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

module pagesetup
 implicit none
 public :: redraw_axes, setpage2
 real, parameter, public :: xlabeloffset = 2.5, ylabeloffset = 4.5

 private

contains
!
!--this subroutine determines the setup of the plotting page
!  sorts out labelling of axes, positioning of windows etc
!  can be used as a replacement for PGENV and PGLABEL
!
!  divides up a single page into subpanels
!
!
!  option to tile graphs appropriately on a page:
!  divides up a single panel into subpanels, with a margin at the edge
!  should replace the call to pgenv and pglabel
!
!  for tiled plots the page setup looks like this:
!
!    |   |   |   |
!  --+---+---+---+--
!    | 1 | 2 | 3 |
!  --+---+---+---+--
!    | 4 | 5 | 6 |
!  --+---+---+---+--
!    |   |   |   |
!
!  (ie. with margins in x and y)
!  note that we divide up a single panel, so pgbeg should be called with nx=1,ny=1
!
!  arguments:
!   iplot  : current plot number
!   nx     : number of panels in x direction
!   ny     : number of panels in y direction
!   xmin,  : xmax, ymin, ymax : plot limits (if tiled should be same for all plots)
!   labelx : x axis label (if tiled should be same for all plots)
!   labely : y axis label (if tiled should be same for all plots)
!   title  : current plot title (can differ between plots)
!   just   : just=1 gives equal aspect ratios (same as in pgenv)
!   axis   : axes options (same as in pgenv with a few extra)
!   vmarginleft,right,bottom,top : initial margin sizes (% of page (if tiled) or panel (if not))
!            (default should be zero for these)
!   tile   : assumes all plots can be tiled
!
!  This version by Daniel Price, July 2006
!
subroutine setpage2(iplotin,nx,ny,xmin,xmax,ymin,ymax,labelx,labely,title,just,axis, &
                    vmarginleftin,vmarginrightin,vmarginbottomin,vmargintopin, &
                    colourbarwidth,titleoffset,isamexaxis,tile,adjustlimits,lastrow,lastplot,&
                    yscale,labelyalt,itransy)
  use plotlib,only:plot_svp,plot_swin,plot_box,plot_qvsz,plot_annotate, &
                   plot_page,plot_qcs,plot_wnad,plot_set_exactpixelboundaries, &
                   plot_qvp
  use asciiutils, only:string_delete
  implicit none
  integer, intent(in) :: iplotin,nx,ny,just,axis,itransy
  real, intent(inout) :: xmin, xmax, ymin, ymax
  real, intent(in)    :: colourbarwidth, titleoffset
  real, intent(in)    :: vmarginleftin,vmarginrightin,vmargintopin,vmarginbottomin
  real, intent(in)    :: yscale
  character(len=*), intent(in) :: labelx,labely,title,labelyalt
  logical, intent(in) :: isamexaxis,tile,adjustlimits,lastrow,lastplot
  integer iplot,ix,iy
  real vptsizeeffx,vptsizeeffy,panelsizex,panelsizey
  real vmargintop,vmarginbottom,vmarginleft,vmarginright
  real vptxmin,vptxmax,vptymin,vptymax
  real aspectratio,devaspectratio,x1,x2,y1,y2
  real xch,ych,dx,dy,xcen,ycen
  character(len=10)  :: xopts, yopts
  logical, parameter :: useexactpixelboundaries = .true.
  logical :: plot_alt_y_axis

  if (axis.eq.4) then
     plot_alt_y_axis = .true.
  else
     plot_alt_y_axis = .false.
  endif
!
! new page if iplot > number of plots on page
!
  if (iplotin.gt.nx*ny) then
     if (mod(iplotin,nx*ny).eq.1)  call plot_page
     iplot = iplotin - (nx*ny)*((iplotin-1)/(nx*ny))
  elseif (iplotin.le.0) then
     return
  else
     iplot = iplotin
  endif
!
! check for errors in input
!
  if (nx.le.0 .or. ny.le.0) return
!
! for tiled plots, adjust effective viewport size if just=1 and graphs are not square
!
  if (tile .and. just.eq.1) then
     if (abs(ymax-ymin) < tiny(ymin)) then
        print*,'setpage: error tiling plots: ymax=ymin'
        return
     endif
!
! query the current aspect ratio of the device and set aspect ratio appropriately
!
     call plot_qvsz(3,x1,x2,y1,y2)
     devaspectratio = (x2-x1)/(y2-y1)
     aspectratio = ((xmax-xmin)*nx)/((ymax-ymin)*ny)/devaspectratio
  else
     aspectratio = 1.0
  endif
!
! set positions of x and y labels in units of character height from edge
!
!  xlabeloffset = 3.0
!  ylabeloffset = 4.5
!
! query the character height as fraction of viewport
!
  call plot_qcs(0,xch,ych)
!
! set margin size in units of viewport dimensions
! allow enough room for the plot labels if they are drawn
! nb: pgplot sets the character height as some fraction of the smallest
!     dimension
!
! for tiled plots, these margins apply to the whole page
! otherwise, these are applied to each panel individually
!
  vmargintop = vmargintopin
  vmarginright = vmarginrightin
  if (axis.ge.0) then
     !--if we are drawing an axis
     !  leave a minimum of half a character
     !  spacing (or 0.7 for antialiasing)
     !  so that axis numbers are not chopped
     !  in half at the edges of the viewport
     vmargintop = max(vmargintop,0.7*ych)
     vmarginright = max(vmarginright,0.7*xch)
     !--leave space for labels
     if (axis.ne.3) then
        vmarginleft = vmarginleftin + (ylabeloffset+1.5)*xch
        vmarginbottom = vmarginbottomin + (xlabeloffset+1.0)*ych
     else
        vmarginleft = vmarginleftin + 2.0*xch
        vmarginbottom = vmarginbottomin + 1.5*ych
     endif
     if (plot_alt_y_axis) vmarginright = vmarginleft

     if (.not.tile) then
        if (ny.gt.1 .and. .not.isamexaxis) then
           vmarginbottom = vmarginbottom + 0.5*ych
        elseif (ny.gt.1) then
           vmarginbottom = vmarginbottom + 0.25*ych
        endif
     endif
  else
     vmarginleft = vmarginleftin
     vmarginbottom = vmarginbottomin
  endif
!
!--set size of each panel
!
  ix = iplot - ((iplot-1)/nx)*nx
  iy = (iplot-1)/nx + 1

  if (tile) then
     !--also leave room for title if necessary
     if (titleoffset.ge.0.) then
        vmargintop = vmargintop + (titleoffset+1.)*ych
     endif

     !
     ! effective viewport size = size - margins (only used for tiled
     !
     vptsizeeffx = 1.0 - vmarginright - vmarginleft
     vptsizeeffy = 1.0 - vmargintop - vmarginbottom
     !     reduce x or y size if just=1 to get right aspect ratio
     if (aspectratio.le.1.0 .and. just.eq.1) then
        if (aspectratio*vptsizeeffy.lt.vptsizeeffx) then
           vptsizeeffx = aspectratio*vptsizeeffy
        !  but this could still be bigger than the margins allow...
        else
           vptsizeeffy = vptsizeeffx/aspectratio
        endif
     elseif (aspectratio.gt.1.0 .and. just.eq.1) then
        if (vptsizeeffx/aspectratio.lt.vptsizeeffy) then
           vptsizeeffy = vptsizeeffx/aspectratio
        !  but this could still be bigger than the margins allow...
        else
           vptsizeeffx = vptsizeeffy*aspectratio
        endif
     endif

     panelsizex = vptsizeeffx/nx
     panelsizey = vptsizeeffy/ny
!         print*,ix,iy,nx,ny
!         print*,panelsizex,panelsizey,vptsizeeffx,vptsizeeffy

!     print*,'margins = ',vmarginleft,vmarginright
     vptxmin = vmarginleft + (ix-1)*panelsizex
     vptxmax = vptxmin + panelsizex
     vptymax = (1.0 - vmargintop) - (iy-1)*panelsizey
     vptymin = vptymax - panelsizey
  else
     !--use full page for non-tiled plots, then set margins inside each panel
     panelsizex = 1.0/nx
     panelsizey = 1.0/ny
     vptxmin = (ix-1)*panelsizex + vmarginleft
     vptxmax = ix*panelsizex - vmarginright
     vptymax = 1.0 - (iy-1)*panelsizey - vmargintop
     vptymin = 1.0 - iy*panelsizey + vmarginbottom


     !--also leave room for title if necessary
     if (titleoffset.ge.0.) then
        vptymax = vptymax - (titleoffset+1.)*ych
     endif

     !--also leave room for colour bar if necessary
     if (colourbarwidth.GT.0.) then
        vptxmax = vptxmax - (colourbarwidth + 1.6)*xch
     endif

  endif
  !    print*,vptxmin,vptxmax,vptymin,vptymax
!
! set viewport
!
 !print*,'setting ',vptxmin,vptxmax,vptymin,vptymax
 call plot_svp(vptxmin,vptxmax,vptymin,vptymax)
!
! set axes
!
 if (just.eq.1) then
    if (nx*ny.eq.1 .and. adjustlimits) then
       !--query viewport aspect ratio
       call plot_qvp(3,x1,x2,y1,y2)
       devaspectratio = (x2-x1)/(y2-y1)
       
       !--adjust limits to match viewport aspect ratio
       dx = xmax - xmin
       dy = ymax - ymin
       if (devaspectratio*dy/dx.ge.1.) then
          xcen = 0.5*(xmin + xmax)
          xmin = xcen - 0.5*devaspectratio*dy
          xmax = xcen + 0.5*devaspectratio*dy
          print*,' auto-adjusting xmin = ',xmin,' xmax = ',xmax
       else
          ycen = 0.5*(ymin + ymax)
          ymin = ycen - 0.5*dx/devaspectratio
          ymax = ycen + 0.5*dx/devaspectratio
          print*,' auto-adjusting ymin = ',ymin,' ymax = ',ymax
       endif
    endif
    call plot_wnad(xmin,xmax,ymin,ymax)
 else
    call plot_swin(xmin,xmax,ymin,ymax)
 endif
!
! adjust viewport to lie exactly on pixel boundaries
!
 if (useexactpixelboundaries) call plot_set_exactpixelboundaries()
!
! option to return before actually doing anything
!
  if (trim(title).eq.'NOPGBOX') return
!
! set options for call to pgbox (draws axes) and label axes where appropriate
! (options are exactly as in pgenv apart from axis=-3,-4 which i have added)
!
  yopts = '*'
  select case(axis)
  case(-4)
     xopts = 'BCT'
  case(-3)
     xopts = 'BCST'
  case(-2)
     xopts = ' '
  case(-1)
    xopts = 'BC'
  case(0,4)
    xopts = 'BCST'
  case(1)
    xopts = 'ABCST'
  case(2)
    xopts = 'ABCGST'
  case(3)
    xopts = 'BCST'
  case(10)
    xopts = 'BCSTL'
    yopts = 'BCST'
  case(20)
    xopts = 'BCST'
    yopts = 'BCSTL'
  case(30)
    xopts = 'BCSTL'
    yopts = 'BCSTL'
  case default
    print*,'setpage: illegal axis argument.'
    xopts = 'BCNST'
  end select
  if (yopts.eq.'*') yopts = xopts
  
  if (plot_alt_y_axis) call string_delete(yopts,'C')
!
! label plot
!
  if (tile) then
     !
     ! decide whether to number and label the y axis
     !
     if (ix.eq.nx .and. axis.ge.0) then
        !
        !--apply label to right hand side axis if used
        !
        if (plot_alt_y_axis) then
           call plot_second_y_axis(yopts,just,axis,itransy,yscale,ylabeloffset,labelyalt)
        endif
     endif
     if (ix.eq.1 .and. axis.ge.0) then
        !
        !--label "normal" y axis
        !
        if (axis.eq.3) then
           yopts = '1N'//trim(yopts)
        else
           yopts = '1VN'//trim(yopts)
           call plot_annotate('L',ylabeloffset,0.5,0.5,labely)
        endif
     elseif (axis.ge.0) then
        !yopts = trim(yopts)//'N'
     endif
     !
     ! decide whether to number and label the x axis
     !
     if ((iy.eq.ny .or. lastplot .or. lastrow) .and. axis.ge.0) then
        xopts = 'N'//trim(xopts)
        if (axis.ne.3) call plot_annotate('B',xlabeloffset,0.5,0.5,labelx)
     endif
     !
     ! plot the title if inside the plot boundaries
     !
     if (titleoffset.lt.0.) call plot_annotate('t',-titleoffset,0.96,1.0,title)

  elseif (axis.ge.0) then
     !
     !--label x axis only if on last row
     !  or if x axis quantities are different
     !
     if (((ny*nx-iplot).lt.nx).or.(.not.isamexaxis).or.lastplot) then
       if (axis.ne.3) call plot_annotate('B',xlabeloffset,0.5,0.5,labelx)
     endif
     !--always plot numbers
     xopts = 'N'//trim(xopts)
     !
     !--apply label to right hand side axis if used
     !
     if (plot_alt_y_axis) then
        call plot_second_y_axis(yopts,just,axis,itransy,yscale,ylabeloffset,labelyalt)
     endif
     !
     !--always label y axis
     !
     if (axis.eq.3) then
        yopts = '1N'//trim(yopts)
     else
        yopts = '1VN'//trim(yopts)
        call plot_annotate('L',ylabeloffset,0.5,0.5,labely)
     endif
     !
     !--always plot title
     !
     call plot_annotate('T',-titleoffset,0.5,0.5,title)

  endif

  call plot_box(xopts,0.0,0,yopts,0.0,0)

  return
end subroutine

!
!--this subroutine is a cut down version of the above, which ONLY redraws the axes
!  (so that axes can be redrawn on *top* of what has been plotted).
!
!  inputs:
!         axis   : axes options (same as in PGENV, with axis=-4,-3,+3 added)
!

subroutine redraw_axes(iaxis,just,yscale,itransy)
  use plotlib, only:plot_box
  implicit none
  integer, intent(in) :: iaxis,just,itransy
  character(len=10) :: xopts, yopts
  real, intent(in)  :: yscale
!
!--set plot axes (options are exactly as in PGENV, with axis=-4,-3,+3 added)
!
  yopts = '*'
  select case(iaxis)
  case(-4)
     xopts = 'BCT'
  case(-3)
     xopts = 'BCST'
  case(-2)
     xopts = ' '
  case(-1)
     xopts = 'BC'
  case(0)
     xopts = 'BCST'
  case(1)
     xopts = 'ABCST'
  case(2)
     xopts = 'ABCGST'
  case(3)
     xopts = 'BCST'
  case(4)
     xopts = 'BCST'
     yopts = 'BST'
  case(10)
     xopts = 'BCSTL'
     yopts = 'BCST'
  case(20)
     xopts = 'BCST'
     yopts = 'BCSTL'
  case(30)
     xopts = 'BCSTL'
     yopts = 'BCSTL'
  case default
     print*,'redraw_axes: illegal AXIS argument.'
     xopts = 'BCST'
  end select
  if (yopts.eq.'*') yopts = xopts

  if (iaxis.eq.4) call plot_second_y_axis(yopts,just,iaxis,itransy,yscale)
  call plot_box(xopts,0.0,0,yopts,0.0,0)

  return
end subroutine redraw_axes

subroutine plot_second_y_axis(yopts,just,iaxis,itransy,yscale,ylabeloffset,labely)
 use plotlib,    only:plot_box,plot_annotate,plot_qwin,plot_swin,plot_wnad
 use asciiutils, only:string_delete
 use transforms, only:transform,transform_inverse,transform_label
 implicit none
 character(len=*), intent(in) :: yopts
 real,             intent(in) :: yscale
 real,             intent(in), optional :: ylabeloffset
 character(len=*), intent(in), optional :: labely
 integer,          intent(in) :: just,iaxis,itransy
 character(len=10)  :: yoptsi
 character(len=120) :: labelyalt
 real :: xmin,xmax,ymin,ymax,yminalt,ymaxalt
 
 yoptsi = yopts
 call string_delete(yoptsi,'B')
 call string_delete(yoptsi,'N')

 !--save plot window settings
 call plot_qwin(xmin,xmax,ymin,ymax)

 !--scaling of y axis: multiplication in un-transformed space
 yminalt = ymin
 ymaxalt = ymax
 if (itransy.gt.0) call transform_inverse(yminalt,ymaxalt,itransy)
 yminalt = yminalt*yscale
 ymaxalt = ymaxalt*yscale
 if (itransy.gt.0) call transform(yminalt,ymaxalt,itransy)
 
 !--set plot window to new scaled y axis
 if (just.eq.1) then
    call plot_wnad(xmin,xmax,yminalt,ymaxalt)
 else
    call plot_swin(xmin,xmax,yminalt,ymaxalt)
 endif
 !--draw axes and label on right hand side of box
 if (iaxis.eq.3) then
    call plot_box(' ',0.0,0,'1MC'//trim(yoptsi),0.0,0) 
 else
    call plot_box(' ',0.0,0,'1VMC'//trim(yoptsi),0.0,0)
 endif
 if (present(labely) .and. present(ylabeloffset)) then
    labelyalt = labely
    if (itransy.gt.0) labelyalt = transform_label(labely,itransy)
    call plot_annotate('R',ylabeloffset,0.5,0.5,labelyalt)
 endif

 !--reset plot window
 if (just.eq.1) then
    call plot_wnad(xmin,xmax,ymin,ymax)
 else
    call plot_swin(xmin,xmax,ymin,ymax)
 endif

end subroutine plot_second_y_axis

end module pagesetup
