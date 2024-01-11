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

!------------------------------------------------------------------------
!  Module containing routines related to plotting the colour bar
!  in various styles
!------------------------------------------------------------------------
module colourbar
 implicit none
 integer, parameter, public :: maxcolourbarstyles = 12
 character(len=28), dimension(0:maxcolourbarstyles), parameter, public :: &
   labelcolourbarstyles = (/'no colour bar               ', &
                            'vertical (right hand side)  ', &
                            'horizontal (underneath plot)', &
                            'plot-hugging vertical       ', &
                            'plot-hugging horizontal     ', &
                            'one-sided vertical          ', &
                            'one-sided horizontal        ', &
                            'floating/inset vertical     ', &
                            'floating/inset horizontal   ', &
                            'vertical (left hand side)   ', &
                            'horizontal (on top)         ', &
                            'custom vertical             ', &
                            'custom horizontal           '/)

 integer, parameter, public :: maxfloatingstyles = 5
 character(len=*), dimension(maxfloatingstyles), parameter, public :: &
   labelfloatingstyles = (/' 1) Top left    ', &
                           ' 2) Top right   ', &
                           ' 3) Bottom left ', &
                           ' 4) Bottom right', &
                           ' 5) Custom      '/)
 !
 !--these are settings that have default values but can
 !  be changed if required
 !
 real, public :: ColourBarDisp = 5.0
 real, public :: ColourBarWidth = 2. ! width in character heights
 logical, public :: iplotcolourbarlabel = .true.
 logical, public :: ilogcolourbaraxis = .false.

 public :: plotcolourbar,incolourbar,incolourbarlabel,barisvertical
 public :: get_colourbarmargins,isfloating,adjustcolourbar,iscustombar
 public :: set_floating_bar_style
 real, private, save :: xlabeloffsetsave = 0.
 real, parameter, private :: dispall = 0.25
 real, public :: ColourBarPosx = 0.01   ! default x pos of short/fat bars
 real, public :: ColourBarPosy = 0.05    ! default y pos of short/fat bars
 real, public :: ColourBarLen = 0.25    ! default length of short/fat bars
 character(len=10), public :: ColourBarFmtStr = 'BCMSTV    '
 private

contains

!--------------------------------------------------------
! this subroutine plots the colour bar in various styles
!--------------------------------------------------------
subroutine plotcolourbar(istyle,icolours,datmin,datmax,label,log, &
           xlabeloffset,vptxminfull,vptxmaxfull,vptyminfull,vptymaxfull)
 use plotlib,   only:plot_set_exactpixelboundaries,plotlib_is_pgplot
 use plotlib,   only:plot_bbuf,plot_ebuf,plot_qwin,plot_qvp,plot_qcs,&
                     plot_svp,plot_swin,plot_imag,plot_box,plot_annotate,plot_gray,&
                     plotlib_extend_pad
 integer, intent(in) :: istyle,icolours
 real, intent(in) :: datmin,datmax,xlabeloffset
 character(len=*), intent(in) :: label
 logical, intent(in) :: log
 real, intent(in), optional :: vptxminfull,vptxmaxfull,vptyminfull,vptymaxfull
 integer, parameter :: maxpixwedg = 400
 real, dimension(6), parameter :: trans = (/-0.5,1.0,0.0,-0.5,0.0,1.0/)
 real, dimension(1,maxpixwedg) :: sampley
 real, dimension(maxpixwedg,1) :: samplex
 character(len=1) :: logc
 integer :: i,npixwedg
 real :: disp,width,xch,ych,dx
 real :: xmin,xmax,ymin,ymax,vptxmin,vptxmax,vptymin,vptymax
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi
 real :: vptxminp,vptxmaxp,vptyminp,vptymaxp
 real :: xmaxpix,xminpix,yminpix,ymaxpix
!
!--return on style 0
!
 if (istyle <= 0) return
!
!--set colour bar displacement and width in character heights
!
 width = ColourBarWidth
 xlabeloffsetsave = xlabeloffset
!
!--start buffering
!
 call plot_bbuf
!
!--query and save current viewport, page and character height settings
!
 call plot_qwin(xmin,xmax,ymin,ymax)
 call plot_qvp(0,vptxmin,vptxmax,vptymin,vptymax)
 call plot_qcs(0,xch,ych)

 !--if colour bar stretches across multiple plots,
 !  override settings for vptymin and vptymax with input values
 if (present(vptxminfull) .and. present(vptxmaxfull) .and. &
     present(vptyminfull) .and. present(vptymaxfull)) then
    vptxmini = vptxminfull
    vptxmaxi = vptxmaxfull
    vptymini = vptyminfull
    vptymaxi = vptymaxfull
 else
    vptxmini = vptxmin
    vptxmaxi = vptxmax
    vptymini = vptymin
    vptymaxi = vptymax
 endif
 disp = dispall

 select case(istyle)
    !------------------------
    ! horizontal colour bar
    !------------------------
 case(2,4,6,8,10,12)

    if (istyle==4) disp = 0. ! plot-hugging
    !
    !--set viewport for the wedge
    !
    if (isfloating(istyle)) then
       vptxminp = vptxmini  ! to
       vptxmaxp = vptxmaxi  ! avoid
       vptyminp = vptymini  ! compiler
       vptymaxp = vptymaxi  ! warnings
       call barlimits(vptxmini,vptxmaxi,vptxminp,vptxmaxp,ColourBarPosx,ColourBarLen)
       call barlimits(vptymini,vptymaxi,vptyminp,vptymaxp,ColourBarPosy,ColourBarLen)
       vptymaxi = vptymini + width*ych
    elseif (istyle==10) then
       ! on top of plot
       vptymini = vptymaxi + (disp+0.1)*ych
       vptymaxi = vptymini + width*ych
    else
       vptymaxi = vptymini - (disp + xlabeloffset)*ych
       vptymini = vptymaxi - width*ych
    endif
    call plot_svp(vptxmini,vptxmaxi,vptymini,vptymaxi)
    call plot_set_exactpixelboundaries()
    !--check number of pixels in colour bar
    call plot_qvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
    npixwedg = max(min(int(xmaxpix-xminpix),maxpixwedg),2)
    !
    !--fill array with all values from datmin to datmax
    !
    dx = (datmax-datmin)/real(npixwedg-1)
    do i=1,npixwedg
       samplex(i,:) = datmin + (i-1)*dx
    enddo
    !
    !--draw colour bar, by cleverly setting window size
    !
    call plot_swin(1.0,real(npixwedg),0.0,1.0)

    if (abs(icolours) > 0) then        ! colour
       !--check if the colour bar will be more than 1024 pixels
       if ((xmaxpix-xminpix) <= 1024 .or. .not.plotlib_is_pgplot) then
          !
          !--the standard way is to use the default line below
          !
          if (icolours==1) then
             call plot_gray(samplex,npixwedg,1,1,npixwedg,1,1,datmin,datmax,trans,iextend=plotlib_extend_pad)
          else
             call plot_imag(samplex,npixwedg,1,1,npixwedg,1,1,datmin,datmax,trans,iextend=plotlib_extend_pad)
          endif
       else
          !
          !--if > 1024 pixels, we instead use the following:
          !  this is a workaround for a PGPLOT bug with large colour bars
          !  (> 1024 device pixels long) - plot colour bar in two halves.
          !  this works up to 2048 pixels, really should divide by n.
          !
          call plot_svp(vptxmini,vptxmaxi-0.5*(vptxmaxi-vptxmini),vptymini,vptymaxi)
          call plot_swin(1.0,real(npixwedg/2),0.0,1.0)
          call plot_set_exactpixelboundaries()
          call plot_imag(samplex,npixwedg,1,1,npixwedg/2,1,1,datmin,datmax,trans)

          call plot_svp(vptxmaxi-0.5*(vptxmaxi-vptxmini)-0.001,vptxmaxi,vptymini,vptymaxi)
          call plot_swin(real(npixwedg/2 + 1),real(npixwedg),0.0,1.0)
          call plot_set_exactpixelboundaries()
          call plot_imag(samplex,npixwedg,1,npixwedg/2+1,npixwedg,1,1,datmin,datmax,trans)
          call plot_svp(vptxmini,vptxmaxi,vptymini,vptymaxi)
          call plot_set_exactpixelboundaries()
       endif
    endif
    call plot_swin(datmin,datmax,0.0,1.0)
    !
    !--draw labelled frame around the wedge
    !
    if (istyle==12) then
       call plot_box(ColourBarFmtStr,0.0,0,'BC',0.0,0)
    elseif (istyle==10) then
       call plot_box('BCSMT',0.0,0,'BC',0.0,0)
    elseif (istyle==4 .or. istyle==6 .or. istyle==8) then
       call plot_box('BNST',0.0,0,'BC',0.0,0)
       if (istyle==6 .or. istyle==8) call plot_box('C',0.0,0,' ',0.0,0)
    else
       call plot_box('BCNST',0.0,0,'BC',0.0,0)
    endif
    !
    !--write the units label: the position is relative to the bottom of
    !  the wedge because of the way we have defined the viewport.
    !  For the horizontal colour bar this never needs to change
    !  (0.25 space + 1 character height for numeric labels + 0.25 space
    !   + 1 character height for actual label = 2.5 character heights)
    !
    if (len_trim(label) > 0 .and. iplotcolourbarlabel) then
       if (istyle==10) then
          call plot_annotate('T',2.5,0.5,0.5,trim(label))
       else
          call plot_annotate('B',2.5,0.5,0.5,trim(label))
       endif
    endif

    !-------------------------------
    ! vertical colour bar (default)
    !-------------------------------
 case default
    if (istyle==3) disp = 0. ! plot-hugging
    !
    !--set viewport for the wedge
    !
    if (isfloating(istyle)) then
       vptxminp = vptxmini  ! to
       vptxmaxp = vptxmaxi  ! avoid
       vptyminp = vptymini  ! compiler
       vptymaxp = vptymaxi  ! warnings
       call barlimits(vptxmini,vptxmaxi,vptxminp,vptxmaxp,ColourBarPosx,ColourBarLen)
       call barlimits(vptymini,vptymaxi,vptyminp,vptymaxp,ColourBarPosy,ColourBarLen)
       vptxmaxi = vptxmini + width*xch
    elseif (istyle==9) then
       vptxmaxi = vptxmini - disp*xch
       vptxmini = vptxmaxi - width*xch
    else
       vptxmini = vptxmaxi + disp*xch
       vptxmaxi = vptxmini + width*xch
    endif
    call plot_svp(vptxmini,vptxmaxi,vptymini,vptymaxi)
    call plot_set_exactpixelboundaries()

    !--check number of pixels in colour bar
    call plot_qvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
    npixwedg = max(min(int(ymaxpix-yminpix),maxpixwedg),2)
    dx = (datmax-datmin)/real(npixwedg-1)
    do i=1,npixwedg
       sampley(:,i) = datmin + (i-1)*dx
    enddo
    !
    !--use log ticks
    !
    logc = ''
    if (ilogcolourbaraxis) logc = 'L'
    !
    !--draw colour bar, by cleverly setting window size
    !
    call plot_swin(0.0,1.0,0.0,real(npixwedg))
    if (icolours==1) then        ! greyscale
       call plot_gray(sampley,1,npixwedg,1,1,1,npixwedg,datmin,datmax,trans,iextend=plotlib_extend_pad)
    elseif (abs(icolours) > 0) then        ! colour
       call plot_imag(sampley,1,npixwedg,1,1,1,npixwedg,datmin,datmax,trans,iextend=plotlib_extend_pad)
    endif
    call plot_swin(0.0,1.0,datmin,datmax)
    !
    !--draw labelled frame around the wedge
    !
    if (istyle==11) then
       call plot_box('BC',0.0,0,ColourBarFmtStr,0.0,0)
    elseif (istyle==9) then
       call plot_box('BC',0.0,0,'BCNSTV',0.0,0)
    elseif (istyle==3 .or. istyle==5 .or. istyle==7) then
       call plot_box('BC',0.0,0,'CMSTV',0.0,0)
       if (istyle==5 .or. istyle==7) call plot_box(' ',0.0,0,'B',0.0,0)
    else
       call plot_box('BC',0.0,0,'BCMSTV'//logc,0.0,0)
    endif
    !
    !--write the units label: the position is relative to the edge of
    !  the wedge because of the way we have defined the viewport.
    !  For the vertical colour bar ColourBarDisp is a set by default to
    !  the maximum size for the numeric label (written horizontally) -
    !  this is about 4 character heights for something like "-5 x 10^10"
    !  We allow the user to adjust this parameter to bring the label
    !  closer where the numeric labels are smaller (e.g. "-5").
    !
    if (len_trim(label) > 0 .and. iplotcolourbarlabel) then
       if (istyle==9) then
          call plot_annotate('L',ColourBarDisp+0.75,1.0,1.0,trim(label))
       else
          call plot_annotate('R',ColourBarDisp+0.75,1.0,1.0,trim(label))
       endif
    endif
 end select
!
!--reset window and viewport
!
 call plot_svp(vptxmin,vptxmax,vptymin,vptymax)
 call plot_swin(xmin,xmax,ymin,ymax)
 call plot_ebuf

 return
end subroutine plotcolourbar

!-------------------------------------------------------
! query function to see if colour bar is plotted
! vertically or horizontally for a given style
!-------------------------------------------------------
logical function barisvertical(istyle)
 integer, intent(in) :: istyle

 barisvertical = .true.
 if (istyle <= 0) return

 select case(istyle)
 case(2,4,6,8,10,12)
    barisvertical = .false.
 case default
    barisvertical = .true.
 end select

end function barisvertical

!-------------------------------------------------------
! query function to see if a given position on
! the plot should lie within the colour bar or not
!-------------------------------------------------------
logical function incolourbar(istyle,iunits,xpt,ypt,xmin,xmax,ymin,ymax)
 use plotlib, only:plot_qcs
 integer, intent(in) :: istyle,iunits
 real, intent(in) :: xpt,ypt,xmin,xmax,ymin,ymax
 real :: xminbar,xmaxbar,yminbar,ymaxbar,xch,ych,barwidth

 incolourbar = .false.
 if (istyle <= 0) return

 select case(istyle)
 case(8,12)
    call barlimits(xminbar,xmaxbar,xmin,xmax,ColourBarPosx,ColourBarLen)
    call barlimits(yminbar,ymaxbar,ymin,ymax,ColourBarPosy,ColourBarLen)
    call plot_qcs(iunits,xch,ych)
    ymaxbar = yminbar + 2.*ColourBarWidth*ych
    if (iplotcolourbarlabel) then
       yminbar = yminbar - 3.0*ych
    else
       yminbar = yminbar - 2.0*ych
    endif
    if ((xpt >= xminbar .and. xpt <= xmaxbar) .and. &
        (ypt >= yminbar .and. ypt <= ymaxbar)) then
       incolourbar = .true.
    endif
 case(7,11)
    call barlimits(xminbar,xmaxbar,xmin,xmax,ColourBarPosx,ColourBarLen)
    call barlimits(yminbar,ymaxbar,ymin,ymax,ColourBarPosy,ColourBarLen)
    call plot_qcs(iunits,xch,ych)
    if (iplotcolourbarlabel) then
       barwidth = (2.*ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch
    else
       barwidth = (2.*ColourBarWidth+0.75 + 5.0)*xch
    endif
    xmaxbar = xminbar + barwidth
    if ((xpt >= xminbar .and. xpt <= xmaxbar) .and. &
        (ypt >= yminbar .and. ypt <= ymaxbar)) then
       incolourbar = .true.
    endif
 case(9)  ! colour bar on left
    if (xpt < xmin) incolourbar = .true.
 case(10) ! colour bar on top
    if (ypt > ymax) incolourbar = .true.
 case(2,4,6)
    if (ypt < ymin) incolourbar = .true.
 case default
    if (xpt > xmax) incolourbar = .true.
 end select

 return
end function incolourbar

!-------------------------------------------------------
! query function to see if a given position on
! the plot should lie within the colour bar label or not
!-------------------------------------------------------
logical function incolourbarlabel(istyle,iunits,xpt,ypt,xmin,xmax,ymin,ymax)
 use plotlib, only:plot_qcs
 integer, intent(in) :: istyle,iunits
 real, intent(in) :: xpt,ypt,xmin,xmax,ymin,ymax
 real :: xch,ych,disp,xminbar,xmaxbar,yminbar,ymaxbar

 incolourbarlabel = .false.
 if (iplotcolourbarlabel) then
    call plot_qcs(iunits,xch,ych)
    disp = dispall
    if (istyle==3 .or. istyle==4) disp = 0.
    select case(istyle)
    case(8,12)
       call barlimits(xminbar,xmaxbar,xmin,xmax,ColourBarPosx,ColourBarLen)
       if (ypt < (ymin-(disp + xlabeloffsetsave + ColourBarWidth+2.0)*ych) .and. &
           ypt > (ymin-(disp + xlabeloffsetsave + ColourBarWidth+3.0)*ych) .and. &
           xpt > xminbar .and. xpt < xmaxbar) incolourbarlabel = .true.
    case(7,11)
       call barlimits(yminbar,ymaxbar,ymin,ymax,ColourBarPosy,ColourBarLen)
       if (xpt > (xmax+(disp + ColourBarWidth-0.25 + max(ColourBarDisp-0.25,0.0))*xch) .and. &
           xpt < (xmax+(disp + ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch) .and. &
           ypt > yminbar .and. ypt < ymaxbar) incolourbarlabel = .true.
    case(2,4,6)
       if (ypt < (ymin-(disp + xlabeloffsetsave + ColourBarWidth+2.0)*ych) .and. &
           ypt > (ymin-(disp + xlabeloffsetsave + ColourBarWidth+3.0)*ych)) incolourbarlabel = .true.
    case(9)
       if (xpt < (xmin-(disp + ColourBarWidth-0.25 + max(ColourBarDisp-0.25,0.0))*xch) .and. &
           xpt > (xmin-(disp + ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch)) incolourbarlabel = .true.
    case(10)
       if (ypt > (ymax+(disp + xlabeloffsetsave + ColourBarWidth+2.0)*ych) .and. &
           ypt < (ymax+(disp + xlabeloffsetsave + ColourBarWidth+3.0)*ych)) incolourbarlabel = .true.
    case default
       if (xpt > (xmax+(disp + ColourBarWidth-0.25 + max(ColourBarDisp-0.25,0.0))*xch) .and. &
           xpt < (xmax+(disp + ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch)) incolourbarlabel = .true.
    end select
 endif

 return
end function incolourbarlabel

!------------------------------------------
! utility function to avoid repeated code
!------------------------------------------
subroutine barlimits(barmin,barmax,posmin,posmax,pos,barlen)
 real, intent(out) :: barmin,barmax
 real, intent(in)  :: posmin,posmax,pos,barlen
 real :: dpos

 dpos = (posmax - posmin) ! in case posmin and barmin are same variable
 barmin = posmin +    pos*dpos
 barmax = barmin + barlen*dpos

end subroutine barlimits

!-------------------------------------------------------
! query function to get margins which should
! be allowed on the page in order to later plot
! the colour bar
!-------------------------------------------------------
subroutine get_colourbarmargins(istyle,xminmargin,xmaxmargin,yminmargin,ymaxmargin,barwidth)
 use plotlib, only:plot_qcs,plot_qvp
 integer, intent(in) :: istyle
 real, intent(inout) :: xminmargin,xmaxmargin,yminmargin,ymaxmargin
 real, intent(out) :: barwidth
 real :: xch,ych,vptxmin,vptxmax,vptymin,vptymax

 barwidth = 0.
 if (istyle <= 0) return
 call plot_qcs(0,xch,ych)
 call plot_qvp(0,vptxmin,vptxmax,vptymin,vptymax)

 if (barisvertical(istyle)) then
    if (iplotcolourbarlabel) then
       barwidth = (ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch
    else
       barwidth = (ColourBarWidth+0.75 + 5.0)*xch
    endif
    if (isfloating(istyle)) then
       barwidth = max((ColourBarPosx-1.) + barwidth,0.)
    endif
    if (istyle==9) then
       xminmargin = xminmargin + barwidth
    else
       xmaxmargin = xmaxmargin + barwidth
    endif
 else
    if (iplotcolourbarlabel) then
       barwidth = (ColourBarWidth+3.0)*ych  ! ie. width + 2.5 + 0.5 margin
    else
       barwidth = (ColourBarWidth+2.0)*ych  ! ie. width + 1.5 + 0.5 margin
    endif
    if (isfloating(istyle)) then
       barwidth = max(-(ColourBarPosy - (barwidth - ColourBarWidth*ych)),0.)
    endif
    if (istyle==10) then
       ymaxmargin = ymaxmargin + barwidth + 0.35*ych
    else
       yminmargin = yminmargin + barwidth
    endif
 endif

 return
end subroutine get_colourbarmargins

!-------------------------------------------------------
! query function for floating colour bar styles
!-------------------------------------------------------
logical function isfloating(istyle)
 integer, intent(in) :: istyle

 select case(istyle)
 case(7,8,11,12)
    isfloating = .true.
 case default
    isfloating = .false.
 end select

end function isfloating

!-------------------------------------------------------
! query function for custom colour bar styles
!-------------------------------------------------------
logical function iscustombar(istyle)
 integer, intent(in) :: istyle

 if (istyle==12 .or. istyle==11) then
    iscustombar = .true.
 else
    iscustombar = .false.
 endif

end function iscustombar

!---------------------------------------------------------------------
! utility function used when interactively changing colour bar limits
!---------------------------------------------------------------------
subroutine adjustcolourbar(istyle,xpt1,ypt1,xpt2,ypt2,&
                           xmin,xmax,ymin,ymax,barmin,barmax)
 integer, intent(in)    :: istyle
 real,    intent(in)    :: xpt1,ypt1,xpt2,ypt2,xmin,xmax,ymin,ymax
 real,    intent(inout) :: barmin,barmax
 real :: dbar,xminbar,xmaxbar,yminbar,ymaxbar

 if (istyle==8 .or. istyle==12) then
    !--floating horizontal bar
    xminbar = xmin    + ColourBarPosx*(xmax - xmin)
    xmaxbar = xminbar + ColourBarLen*(xmax - xmin)
 else
    xminbar = xmin
    xmaxbar = xmax
 endif
 if (istyle==7 .or. istyle==11) then
    !--floating vertical bar
    yminbar = ymin    + ColourBarPosy*(ymax - ymin)
    ymaxbar = yminbar + ColourBarLen*(ymax - ymin)
 else
    yminbar = ymin
    ymaxbar = ymax
 endif

 if (barisvertical(istyle)) then
    if ((ymaxbar-yminbar) > 0.) then
       dbar = (barmax-barmin)/(ymaxbar-yminbar)
    else
       dbar = 0.
    endif
    barmax = barmin + (max(ypt1,ypt2)-yminbar)*dbar
    barmin = barmin + (min(ypt1,ypt2)-yminbar)*dbar
 else
    if ((xmaxbar-xminbar) > 0.) then
       dbar = (barmax-barmin)/(xmaxbar-xminbar)
    else
       dbar = 0.
    endif
    barmax = barmin + (max(xpt1,xpt2)-xminbar)*dbar
    barmin = barmin + (min(xpt1,xpt2)-xminbar)*dbar
 endif

end subroutine adjustcolourbar

subroutine set_floating_bar_style(iColourBarStyle,iColourBarPos)
 integer, intent(in) :: iColourBarStyle,iColourBarPos

 select case(iColourBarPos)
 case(1)
    if (barisvertical(iColourBarStyle)) then
       ColourBarPosx = 0.01
       ColourBarPosy = 0.74
    else
       ColourBarPosx = 0.01
       ColourBarPosy = 0.95
    endif
 case(2)
    if (barisvertical(iColourBarStyle)) then
       ColourBarPosx = 0.82 ! minus width in ch
       ColourBarPosy = 0.74
    else
       ColourBarPosx = 0.73
       ColourBarPosy = 0.95
    endif
 case(3)
    if (barisvertical(iColourBarStyle)) then
       ColourBarPosx = 0.01
       ColourBarPosy = 0.01
    else
       ColourBarPosx = 0.015
       ColourBarPosy = 0.075
    endif
 case(4)
    if (barisvertical(iColourBarStyle)) then
       ColourBarPosx = 0.82 ! minus width in ch
       ColourBarPosy = 0.01
    else
       ColourBarPosx = 0.73
       ColourBarPosy = 0.075
    endif
 end select

end subroutine set_floating_bar_style

end module colourbar
