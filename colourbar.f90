!------------------------------------------------------------------------
!  Module containing routines related to plotting the colour bar
!  in various styles
!
!  DJP 16/10/2007
!
!  Distributed as part of SPLASH: a visualisation tool for SPH data
!------------------------------------------------------------------------
module colourbar
 implicit none
 integer, parameter, public :: maxcolourbarstyles = 6
 character(len=28), dimension(0:maxcolourbarstyles), parameter, public :: &
   labelcolourbarstyles = (/'no colour bar               ', &
                            'vertical (right hand side)  ', &
                            'horizontal (underneath plot)', &
                            'plot-hugging vertical       ', &
                            'plot-hugging horizontal     ', &
                            'one-sided vertical          ', &
                            'one-sided horizontal        '/)
 !
 !--these are settings that have default values but can
 !  be changed if required
 !
 real, public :: ColourBarDisp = 5.0
 real, parameter, public :: ColourBarWidth = 2. ! width in character heights
 logical, public :: iplotcolourbarlabel = .true.

 public :: plotcolourbar,incolourbar,incolourbarlabel,barisvertical,get_colourbarmargins
 real, private :: xlabeloffsetsave = 0.
 real, parameter, private :: dispall = 0.25
 private
 
contains

!-------------------------------------------------------
! this subroutine plots the colour bar 
! (similar to PGWEDG: differences are that
! text is written vertically, txtsep is a
! changeable parameter and the character height
! is not changed)
!-------------------------------------------------------
subroutine plotcolourbar(istyle,icolours,datmin,datmax,label,log, &
           xlabeloffset,vptxminfull,vptxmaxfull,vptyminfull,vptymaxfull)
 use pagesetup, only:set_exactpixelboundaries
 implicit none
 integer, intent(in) :: istyle,icolours
 real, intent(in) :: datmin,datmax,xlabeloffset
 character(len=*), intent(in) :: label
 logical, intent(in) :: log
 real, intent(in), optional :: vptxminfull,vptxmaxfull,vptyminfull,vptymaxfull
 integer, parameter :: npixwedg = 400
 real, dimension(6), parameter :: trans = (/-0.5,1.0,0.0,-0.5,0.0,1.0/)
 real, dimension(1,npixwedg) :: sampley
 real, dimension(npixwedg,1) :: samplex
 integer :: i
 real :: disp,width,xch,ych,dx
 real :: xmin,xmax,ymin,ymax,vptxmin,vptxmax,vptymin,vptymax
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi
 real :: xmaxpix,xminpix,yminpix,ymaxpix
!
!--return on style 0
!
 if (istyle.le.0) return
!
!--set colour bar displacement and width in character heights
!
 width = ColourBarWidth
 xlabeloffsetsave = xlabeloffset
!
!--start buffering
!
 call pgbbuf
!
!--query and save current viewport, page and character height settings
!
 call pgqwin(xmin,xmax,ymin,ymax)
 call pgqvp(0,vptxmin,vptxmax,vptymin,vptymax)
 call pgqcs(0,xch,ych)
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
!
!--fill array with all values from datmin to datmax
!
 dx = (datmax-datmin)/real(npixwedg-1)
 do i=1,npixwedg
    sampley(1,i) = datmin + (i-1)*dx
    samplex(i,1) = datmin + (i-1)*dx
 enddo
 disp = dispall

 select case(istyle)
 !------------------------
 ! horizontal colour bar
 !------------------------
 case(2,4,6)

   if (istyle.eq.4) disp = 0. ! plot-hugging
   !
   !--set viewport for the wedge
   !
   vptymaxi = vptymini - (disp + xlabeloffset)*ych
   vptymini = vptymaxi - width*ych
   call pgsvp(vptxmini,vptxmaxi,vptymini,vptymaxi)
   call set_exactpixelboundaries()
   !
   !--draw colour bar, by cleverly setting window size
   !
   call pgswin(1.0,real(npixwedg),0.0,1.0)

   !--check number of pixels in colour bar
   call pgqvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
   
   if (abs(icolours).gt.0) then        ! colour
   !--check if the colour bar will be more than 1024 pixels
       if ((xmaxpix-xminpix).le.1024) then
    !   
    !--the standard way is to use the default line below
    !
          call pgimag(samplex,npixwedg,1,1,npixwedg,1,1,datmin,datmax,trans)   
       else
    !
    !--if > 1024 pixels, we instead use the following: 
    !  this is a workaround for a PGPLOT bug with large colour bars
    !  (> 1024 device pixels long) - plot colour bar in two halves.
    !  this works up to 2048 pixels, really should divide by n.
    !
          call pgsvp(vptxmini,vptxmaxi-0.5*(vptxmaxi-vptxmini),vptymini,vptymaxi)
          call pgswin(1.0,real(npixwedg/2),0.0,1.0)
          call set_exactpixelboundaries()
          call pgimag(samplex,npixwedg,1,1,npixwedg/2,1,1,datmin,datmax,trans)

          call pgsvp(vptxmaxi-0.5*(vptxmaxi-vptxmini)-0.001,vptxmaxi,vptymini,vptymaxi)
          call pgswin(real(npixwedg/2 + 1),real(npixwedg),0.0,1.0)
          call set_exactpixelboundaries()
          call pgimag(samplex,npixwedg,1,npixwedg/2+1,npixwedg,1,1,datmin,datmax,trans)
          call pgsvp(vptxmini,vptxmaxi,vptymini,vptymaxi)
          call set_exactpixelboundaries()
       endif
    endif
    call pgswin(datmin,datmax,0.0,1.0)
   !
   !--draw labelled frame around the wedge
   !
    if (istyle.eq.4 .or. istyle.eq.6) then
       call pgbox('BNST',0.0,0,'BC',0.0,0)
       if (istyle.eq.6) call pgbox('C',0.0,0,' ',0.0,0)    
    else
       call pgbox('BCNST',0.0,0,'BC',0.0,0)
    endif
   !
   !--write the units label: the position is relative to the bottom of 
   !  the wedge because of the way we have defined the viewport. 
   !  For the horizontal colour bar this never needs to change 
   !  (0.25 space + 1 character height for numeric labels + 0.25 space
   !   + 1 character height for actual label = 2.5 character heights)
   !
    if (label.ne.' ' .and. iplotcolourbarlabel) then
       call pgmtxt('B',2.5,0.5,0.5,trim(label))
    endif

 !-------------------------------
 ! vertical colour bar (default)
 !-------------------------------
 case default
    
    if (istyle.eq.3) disp = 0. ! plot-hugging
   !
   !--set viewport for the wedge
   !
    vptxmini = vptxmaxi + disp*xch
    vptxmaxi = vptxmini + width*xch
    call pgsvp(vptxmini,vptxmaxi,vptymini,vptymaxi)
    call set_exactpixelboundaries()
   !
   !--draw colour bar, by cleverly setting window size
   !
    call pgswin(0.0,1.0,1.0,real(npixwedg))
   ! if (abs(icolours).eq.1) then        ! greyscale
   !    call pggray(sample,1,npixwedg,1,1,1,npixwedg,datmin,datmax,trans)
    if (abs(icolours).gt.0) then        ! colour
       call pgimag(sampley,1,npixwedg,1,1,1,npixwedg,datmin,datmax,trans)
    endif
    call pgswin(0.0,1.0,datmin,datmax)
   !
   !--draw labelled frame around the wedge
   !
    if (istyle.eq.3 .or. istyle.eq.5) then
       call pgbox('BC',0.0,0,'CMSTV',0.0,0)
       if (istyle.eq.5) call pgbox(' ',0.0,0,'B',0.0,0)
    else
       call pgbox('BC',0.0,0,'BCMSTV',0.0,0)
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
    if (label.ne.' ' .and. iplotcolourbarlabel) then
       call pgmtxt('R',ColourBarDisp+0.75,1.0,1.0,trim(label))
    endif
 end select
!
!--reset window and viewport
!
 call pgsvp(vptxmin,vptxmax,vptymin,vptymax)
 call pgswin(xmin,xmax,ymin,ymax)
 call pgebuf

 return
end subroutine plotcolourbar

!-------------------------------------------------------
! query function to see if colour bar is plotted
! vertically or horizontally for a given style
!-------------------------------------------------------
logical function barisvertical(istyle)
 implicit none
 integer, intent(in) :: istyle

 barisvertical = .true. 
 if (istyle.le.0) return
 
 select case(istyle)
 case(2,4,6)
    barisvertical = .false.
 case default
    barisvertical = .true.
 end select

end function barisvertical

!-------------------------------------------------------
! query function to see if a given position on
! the plot should lie within the colour bar or not
!-------------------------------------------------------
logical function incolourbar(istyle,xpt,ypt,xmin,xmax,ymin,ymax)
 implicit none
 integer, intent(in) :: istyle
 real, intent(in) :: xpt,ypt,xmin,xmax,ymin,ymax
 
 incolourbar = .false.
 if (istyle.le.0) return
 
 select case(istyle)
 case(2,4,6)
    if (ypt.lt.ymin) incolourbar = .true.
 case default
    if (xpt.gt.xmax) incolourbar = .true.
 end select
 
 return
end function incolourbar

!-------------------------------------------------------
! query function to see if a given position on
! the plot should lie within the colour bar or not
!-------------------------------------------------------
logical function incolourbarlabel(istyle,xpt,ypt,xmin,xmax,ymin,ymax)
 implicit none
 integer, intent(in) :: istyle
 real, intent(in) :: xpt,ypt,xmin,xmax,ymin,ymax
 real :: xch,ych,disp
 
 incolourbarlabel = .false.
 if (iplotcolourbarlabel) then
    call pgqcs(4,xch,ych)
    print*,'checking colourbar label ',xpt,ypt,ymin-2.5*ych,ych
    disp = dispall
    if (istyle.eq.3 .or. istyle.eq.4) disp = 0.
    select case(istyle)
    case(2,4,6)
       if (ypt.lt.(ymin-(disp + xlabeloffsetsave + ColourBarWidth+2.0)*ych) .and. &
           ypt.gt.(ymin-(disp + xlabeloffsetsave + ColourBarWidth+3.0)*ych)) incolourbarlabel = .true.
    case default
       if (xpt.gt.(xmax+(disp + ColourBarWidth-0.25 + max(ColourBarDisp-0.25,0.0))*xch) .and. &
           xpt.lt.(xmax+(disp + ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch)) incolourbarlabel = .true.
    end select
 endif
 
 return
end function incolourbarlabel

!-------------------------------------------------------
! query function to get margins which should
! be allowed on the page in order to later plot
! the colour bar
!-------------------------------------------------------
subroutine get_colourbarmargins(istyle,xmaxmargin,yminmargin,barwidth)
 implicit none
 integer, intent(in) :: istyle
 real, intent(inout) :: xmaxmargin,yminmargin
 real, intent(out) :: barwidth
 real :: xch,ych

 barwidth = 0.
 if (istyle.le.0) return
 call pgqcs(0,xch,ych)

 select case(istyle)
 case(2,4,6)
    if (iplotcolourbarlabel) then
       barwidth = (ColourBarWidth+3.0)*ych  ! ie. width + 2.5 + 0.5 margin
    else
       barwidth = (ColourBarWidth+2.0)*ych  ! ie. width + 1.5 + 0.5 margin
    endif
    yminmargin = yminmargin + barwidth
 case default
    if (iplotcolourbarlabel) then
       barwidth = (ColourBarWidth+0.75 + max(ColourBarDisp+0.75,0.0))*xch
    else
       barwidth = (ColourBarWidth+0.75 + 5.0)*xch
    endif
    xmaxmargin = xmaxmargin + barwidth
 end select
 
 return
end subroutine get_colourbarmargins

end module colourbar
