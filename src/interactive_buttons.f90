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

 real, parameter, public :: xleft_vp0 = -0.01, ytop_vp0=1.02
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
 !button_pressed = 0
 buttons_drawn = .false.
 nbuttons = 0

end subroutine erase_buttons

!---------------------------------------------
! draw a single button with position and type
!---------------------------------------------
subroutine draw_button(i,itype,msg)
 use plotlib, only:plot_qwin,plot_qcs,plot_qfs,plot_sfs,plot_rect,plot_sls,&
                   plot_set_opacity,plot_set_clipping,plot_get_clipping,&
                   plot_qls,plot_qch,plot_sch,plot_pt1,plot_text,plot_line1,plot_poly
 use interactive_utils, only:get_posxy
 integer, intent(in) :: i,itype
 character(len=*), intent(in) :: msg
 real :: xmin,xmax,ymin,ymax
 real :: xleft_vp,ytop_vp,oldch,x0,y0,width
 real :: xleft,xright,ybottom,ytop,xch,ych,xsize,ysize,xbuf,ybuf
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
 call plot_qcs(0,xch,ych)
 xsize = 2.0*xch
 ysize = 1.5*ych
 xbuf = 0.075*xch
 ybuf = 0.075*ych
 if (itype==ibutton_erase) then
    xleft_vp = xleft_vp0 - xbuf
    ytop_vp = ytop_vp0 + ybuf
    ysize = ysize + 2.*ybuf ! erase an area slightly larger than original drawing
    xsize = xsize + xbuf
    width = i*xsize + xbuf
 else
    xleft_vp = xleft_vp0 + (i-1)*(xsize+xbuf)
    ytop_vp = ytop_vp0
    width = xsize
 endif
 ! translate from viewport to world coords
 call get_posxy(xleft_vp,ytop_vp,xleft,ytop,xmin,xmax,ymin,ymax)
 ! translate from viewport to world coords
 call get_posxy(xleft_vp+width,ytop_vp-ysize,xright,ybottom,xmin,xmax,ymin,ymax)

 call plot_qcs(4,xch,ych)
 buttons(i)%x(1) = xleft
 buttons(i)%x(2) = xright
 buttons(i)%y(1) = ybottom
 buttons(i)%y(2) = ytop
 buttons(i)%type = itype
 buttons(i)%help = trim(msg)
 !print*,'button ',i,buttons(i)%x,buttons(i)%y
 if (itype==ibutton_erase) then
    call plot_rect(buttons(i)%x(1),buttons(i)%x(2),buttons(i)%y(1),buttons(i)%y(2))
 else
    call plot_rect(buttons(i)%x(1),buttons(i)%x(2),buttons(i)%y(1),buttons(i)%y(2),0.2*ych)
 endif
 x0 = 0.5*(xleft + xright)
 y0 = 0.5*(ybottom + ytop)

 ! draw the pattern inside the button
 call plot_sls(2) ! dashed
 call plot_sfs(2) ! hollow
 select case(itype)
 case(ibutton_rectangle)
    xbuf = 0.33*xch
    ybuf = 0.33*ych
    call plot_sch(0.33)
    call plot_rect(buttons(i)%x(1)+xbuf,buttons(i)%x(2)-xbuf,buttons(i)%y(1)+ybuf,buttons(i)%y(2)-ybuf)
 case(ibutton_forward)
    call plot_pt1(x0,y0,29)
 case(ibutton_backward)
    call plot_pt1(x0,y0,28)
 case(ibutton_circle)
    call plot_pt1(x0,y0,23)
 case(ibutton_text)
    xbuf = 0.6*xch
    ybuf = 0.33*ych
    call plot_text(buttons(i)%x(1)+xbuf,buttons(i)%y(1)+ybuf,'T')
 case(ibutton_plus)
    xbuf = 0.6*xch
    ybuf = 0.45*ych
    call plot_text(buttons(i)%x(1)+xbuf,buttons(i)%y(1)+ybuf,'+')
 case(ibutton_minus)
    xbuf = 0.8*xch
    ybuf = 0.45*ych
    call plot_text(buttons(i)%x(1)+xbuf,buttons(i)%y(1)+ybuf,'-')
 case(ibutton_adapt)
    xbuf = 0.1*xch
    ybuf = 0.1*ych
!    call plot_sah(1,20.,1.0)
!    call plot_sah(2,45.,0.7)
   ! call plot_sah(2,45.0,0.7)
    xsize = 5.5
    call plot_sls(1)
    call plot_line1(x0+xbuf,y0+ybuf,x0+xsize*xbuf,y0+xsize*ybuf)
    call plot_line1(x0-xbuf,y0-ybuf,x0-xsize*xbuf,y0-xsize*ybuf)
    call plot_line1(x0-xbuf,y0+ybuf,x0-xsize*xbuf,y0+xsize*ybuf)
    call plot_line1(x0+xbuf,y0-ybuf,x0+xsize*xbuf,y0-xsize*ybuf)
 case(ibutton_irregular)
    call plot_sch(0.33)
    xpts(:) = (/0.1258,0.1830,0.3393,0.5224,0.7149,0.8179,0.8279,&
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
 if (ibutton > 0) call print_message('                           ',&
                                     buttons(ibutton)%help)

end subroutine print_button_help

end module interactive_buttons
