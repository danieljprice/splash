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
module interactive_buttons
 implicit none
 integer, parameter :: maxbuttons = 10
 integer, parameter :: maxtypes = 6
 integer, parameter, public :: &
    ibutton_rectangle = 1, &
    ibutton_circle    = 2, &
    ibutton_irregular = 3, &
    ibutton_text      = 4, &
    ibutton_forward   = 5, &
    ibutton_backward  = 6, &
    ibutton_plus      = 7, &
    ibutton_minus     = 8, &
    ibutton_adapt     = 9, &
    ibutton_erase     = 10 ! erase must always be the last one

 type button
    real :: x(2),y(2)
    integer :: type
    character(len=64) :: help
 end type button

 integer, parameter, public  :: max_button_instant = 5

 ! these are just temporarily set to some default values
 integer, parameter, private :: button_default = 6
 integer, public :: button_pressed = button_default
 logical, private :: buttons_drawn = .false.
 integer, private :: nbuttons = 0

 type(button) :: buttons(maxbuttons)

 public :: draw_buttons,inbutton,erase_buttons,button_type
 public :: press_button,press_button_if_inside,print_button_help

 interface press_button
   module procedure press_button,press_button_default
 end interface press_button

 private

contains

!---------------------------------------------
! draw the full set of buttons
! in principle this routine could be separate
!---------------------------------------------
subroutine draw_buttons(onclick)
 logical, intent(in), optional :: onclick
 logical :: my_onclick

 my_onclick = .true.
 if (present(onclick)) my_onclick = onclick

 if (buttons_drawn .and. .not.my_onclick) return

 ! instant-action buttons, change max_button_instant if you add more
 call draw_button(1,ibutton_backward,'previous timestep (keystroke ''b'')')    ! backward arrow
 call draw_button(2,ibutton_forward,'next timestep (space bar)')         ! forward arrow
 call draw_button(3,ibutton_plus,'zoom in (keystroke ''+'')')            ! plus
 call draw_button(4,ibutton_minus,'zoom out (keystroke ''-'')')          ! minus
 call draw_button(5,ibutton_adapt,'adapt plot limits (keystroke ''a'')') ! a

 ! deferred-action buttons
 call draw_button(6,ibutton_rectangle,'rectangle selection (left click)') ! rectangle selection
 call draw_button(7,ibutton_irregular,'irregular selection (shift + left click)') ! rectangle selection
 call draw_button(8,ibutton_circle,'circular selection (right click)')     ! circle selection
 call draw_button(9,ibutton_text,'add text annotation (ctrl-t)')      ! text annotation
 buttons_drawn = .true.

end subroutine draw_buttons

!---------------------------------------------
! erase the screen area behind the button set
!---------------------------------------------
subroutine erase_buttons()
 use plotlib, only:plot_sci,plot_sfs

 if (.not.buttons_drawn) return ! do not erase if nothing drawn
 call reset_buttons()
 call plot_sci(0)
 call draw_button(nbuttons,ibutton_erase,'')
 call plot_sci(1)
 buttons_drawn = .false.
 nbuttons = 0

end subroutine erase_buttons

!---------------------------------------------
! draw a single button with position and type
!---------------------------------------------
 subroutine draw_button(i,itype,msg)
 use plotlib, only:plot_qwin,plot_qcs,plot_qfs,plot_sfs,plot_rect,plot_sls,&
                   plot_set_opacity,plot_set_clipping,plot_get_clipping,&
                   plot_qls,plot_qch,plot_sch,plot_pt1,plot_line1,plot_poly
 use interactive_utils, only:get_posxy,get_button_anchor,get_pix2viewport_factor
 integer, intent(in) :: i,itype
 character(len=*), intent(in) :: msg
 real :: xmin,xmax,ymin,ymax
 real :: xleft_vp,ybottom_vp,oldch,x0,y0,width
 real :: xleft,xright,ybottom,ytop,xch,ych,xsize,ysize,xbuf,ybuf
 real :: pix_to_viewport_x,pix_to_viewport_y
 real :: xleft_vp0,ybottom_vp0,fac,xadapt_size
 integer :: oldfill,oldclip,oldls
 integer, parameter :: npts = 14
 real :: xpts(npts),ypts(npts)

 if (i > maxbuttons .or. i < 0) return
 nbuttons = i
 ! query and save previous settings
 call plot_qfs(oldfill)
 call plot_qls(oldls)
 call plot_qch(oldch)
 call plot_get_clipping(oldclip)

 if (itype==ibutton_erase) then
    call plot_sfs(1) ! solid
    call plot_set_opacity(1.0)
 elseif (button_pressed==i) then
    call plot_sfs(1) ! solid
    call plot_set_opacity(0.25)
 else
    call plot_sfs(2) ! hollow
    call plot_set_opacity(0.25)
 endif
 call plot_set_clipping(0)
 call plot_qwin(xmin,xmax,ymin,ymax)
 call get_pix2viewport_factor(pix_to_viewport_x,pix_to_viewport_y)

 ! we start 18 pixels below the bottom margin
 ysize = 14.0 * pix_to_viewport_y
 ybuf = 1.0 * pix_to_viewport_y   ! pixel buffer

 ! x size is less constrained, but use same buffer size
 xsize = 20. * pix_to_viewport_x  ! golden ratio times y size
 xbuf = 1.0 * pix_to_viewport_x

 ! set character height to 12 pixels
 call plot_sch(1.0)
 call plot_qcs(0,xch,ych)
 fac = 1.0
 if (ysize > 0.85714*ych) fac = 0.85714*ych/ysize  ! 12 point font
 call plot_sch(fac)

 ! find the anchor position for the button set
 call get_button_anchor(xleft_vp0,ybottom_vp0)

 if (itype==ibutton_erase) then
    xleft_vp = xleft_vp0 - xbuf
    ybottom_vp = ybottom_vp0 - ybuf
    ysize = ysize + 2.*ybuf ! erase an area slightly larger than original drawing
    xsize = xsize + 2.*xbuf
    width = i*xsize + xbuf
 else
    xleft_vp = xleft_vp0 + (i-1)*(xsize+xbuf)
    ybottom_vp = ybottom_vp0
    width = xsize
 endif
 ! translate from viewport to world coords
 call get_posxy(xleft_vp,ybottom_vp,xleft,ybottom,xmin,xmax,ymin,ymax)
 ! translate from viewport to world coords
 call get_posxy(xleft_vp+width,ybottom_vp+ysize,xright,ytop,xmin,xmax,ymin,ymax)

 call plot_qcs(4,xch,ych)

 buttons(i)%x(1) = xleft
 buttons(i)%x(2) = xright
 buttons(i)%y(1) = ybottom
 buttons(i)%y(2) = ytop
 buttons(i)%type = itype
 buttons(i)%help = trim(msg)

 if (itype==ibutton_erase) then
    call plot_rect(buttons(i)%x(1),buttons(i)%x(2),buttons(i)%y(1),buttons(i)%y(2),0.2*ych)
 else
    if (button_pressed==i) then
       call plot_sfs(2) ! hollow
       call plot_rect(buttons(i)%x(1),buttons(i)%x(2),buttons(i)%y(1),buttons(i)%y(2),0.2*ych)
       call plot_sfs(1) ! solid
    endif
    call plot_set_opacity(0.25)
    call plot_rect(buttons(i)%x(1),buttons(i)%x(2),buttons(i)%y(1),buttons(i)%y(2),0.2*ych)
 endif 
 xsize = xright-xleft
 ysize = ytop-ybottom
 x0 = 0.5*(xleft + xright)
 y0 = 0.5*(ybottom + ytop)

 ! draw the pattern inside the button
 call plot_sls(2) ! dashed
 call plot_sfs(2) ! hollow
 select case(itype)
 case(ibutton_rectangle)
    xbuf = 0.33*xch
    ybuf = 0.33*ych
    call plot_sch(0.33*fac)
    call plot_rect(buttons(i)%x(1)+xbuf,buttons(i)%x(2)-xbuf,buttons(i)%y(1)+ybuf,buttons(i)%y(2)-ybuf)
 case(ibutton_forward)
    call plot_pt1(x0,y0,29)
 case(ibutton_backward)
    call plot_pt1(x0,y0,28)
 case(ibutton_circle)
    call plot_sls(1)
    call plot_pt1(x0,y0,22)
 case(ibutton_text)
    call plot_button_text_centered(x0,y0,'a')
 case(ibutton_plus)
    call plot_button_text_centered(x0,y0,'+')
 case(ibutton_minus)
    call plot_button_text_centered(x0,y0,'-')
 case(ibutton_adapt)
    xbuf = 0.1*xsize
    ybuf = 0.1*ysize
    xadapt_size = 3.9
    call plot_sls(1)
    call plot_line1(x0+xbuf,y0+ybuf,x0+xadapt_size*xbuf,y0+xadapt_size*ybuf)
    call plot_line1(x0-xbuf,y0-ybuf,x0-xadapt_size*xbuf,y0-xadapt_size*ybuf)
    call plot_line1(x0-xbuf,y0+ybuf,x0-xadapt_size*xbuf,y0+xadapt_size*ybuf)
    call plot_line1(x0+xbuf,y0-ybuf,x0+xadapt_size*xbuf,y0-xadapt_size*ybuf)
 case(ibutton_irregular)
    call plot_sch(0.33*fac)
    xpts(:) = (/0.1358,0.1830,0.3393,0.5224,0.7149,0.8179,0.8279,&
                0.8479,0.8370,0.8127,0.7187,0.4499,0.2631,0.1487/)
    ypts(:) = (/0.5490,0.7558,0.7616,0.7888,0.7727,0.7462,0.6163,&
                0.4769,0.4346,0.3383,0.2358,0.1834,0.1787,0.3446/)
    xpts(:) = buttons(i)%x(1) + xpts(:)*(buttons(i)%x(2) - buttons(i)%x(1))
    ypts(:) = buttons(i)%y(1) + ypts(:)*(buttons(i)%y(2) - buttons(i)%y(1))
    call plot_poly(npts,xpts,ypts)
 end select

 ! restore settings
 call plot_sch(oldch)
 call plot_set_opacity(1.0)
 call plot_sfs(oldfill)
 call plot_sls(oldls)
 call plot_set_clipping(oldclip)

end subroutine draw_button

!---------------------------------------------
! plot labelled button but with the text
! in the button centred both horizontally
! and vertically at the position xc,yc
!---------------------------------------------
subroutine plot_button_text_centered(xc,yc,text)
 use plotlib, only:plot_qtxt,plot_text
 real, intent(in) :: xc,yc
 character(len=*), intent(in) :: text
 real :: xbox(4),ybox(4),xbp(4),ybp(4),x0,y0,y_plus_mid

 call plot_qtxt(xc,yc,0.0,0.0,text,xbox,ybox)
 x0 = 0.5*(minval(xbox)+maxval(xbox))
 y0 = 0.5*(minval(ybox)+maxval(ybox))
!
!--hyphen: the full bounding box is much taller than the bar, so centring
!  that box puts the bar too high; the bar aligns with the math axis, so
!  use the vertical centre of '+' at the same reference as the y anchor
!  (same as plot_text / plot_qtxt with fjust=0).
!
 if (trim(text)=='-') then
    call plot_qtxt(xc,yc,0.0,0.0,'+',xbp,ybp)
    y_plus_mid = 0.5*(minval(ybp)+maxval(ybp))
    call plot_text(2.0*xc - x0,2.0*yc - y_plus_mid,text)
 else
    call plot_text(2.0*xc - x0,2.0*yc - y0,text)
 endif

end subroutine plot_button_text_centered

!---------------------------------------------
! reset all parameters from buttons
!---------------------------------------------
subroutine reset_buttons()
 integer :: i

 do i=1,size(buttons)
    buttons(i)%x = 0.
    buttons(i)%y = 0.
 enddo
 !button_pressed = 0

end subroutine reset_buttons

!-------------------------------------------------------
! check if a point is inside a button
!-------------------------------------------------------
integer function inbutton(xpt,ypt)
 use interactive_utils, only:inrectangle
 real, intent(in) :: xpt,ypt
 integer :: i

 inbutton = 0
 do i=1,nbuttons
    if (inrectangle(xpt,ypt,buttons(i)%x(1),buttons(i)%x(2),buttons(i)%y(1),buttons(i)%y(2))) then
       inbutton = i
       return
    endif
 enddo

end function inbutton

!---------------------------------------------
! query what type of buttonw was pressed
!---------------------------------------------
integer function button_type(i)
 integer, intent(in) :: i

 button_type = 0
 if (i > 0 .and. i <= nbuttons) button_type = buttons(i)%type

end function button_type

!---------------------------------------------
! designate a button as pressed
!---------------------------------------------
subroutine press_button(i)
 integer, intent(in) :: i

 if (i > 0 .and. i <= nbuttons) button_pressed = i
 call erase_buttons()
 call draw_buttons(onclick=.true.)

end subroutine press_button

!---------------------------------------------
! press the default button
!---------------------------------------------
subroutine press_button_default()

 ! do this in case buttons are not drawn
 button_pressed = button_default

 ! should have same effect, but nbuttons=0 if buttons not drawn
 call press_button(button_default)

end subroutine press_button_default

!---------------------------------------------
! click a button if inside
!---------------------------------------------
integer function press_button_if_inside(xpt,ypt) result(i)
 real, intent(in) :: xpt,ypt

 i = inbutton(xpt,ypt)
 if (i > 0) call press_button(i)

end function press_button_if_inside

!---------------------------------------------
! print help message for a button
!---------------------------------------------
subroutine print_button_help(xpt,ypt,ibutton)
 use interactive_help, only:print_message
 real, intent(in) :: xpt,ypt
 integer, intent(out) :: ibutton

 ibutton = inbutton(xpt,ypt)
 if (ibutton > 0) call print_message(buttons(ibutton)%help)

end subroutine print_button_help

end module interactive_buttons
