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
module interactive_help
 implicit none
 character(len=128) :: help_msg = '', last_msg = ''

 public :: print_message,print_last_message
 public :: clear_messages

 interface print_message
  module procedure print_message,print_var,print_message_xy,&
                   print_message_single,print_last_message
 end interface

 private

contains
!--------------------------------
! print message into help area
!--------------------------------
subroutine print_message(line1,line2,onclick)
 use plotlib,           only:plot_qch,plot_qwin,plot_sch,plot_qcs,plot_set_opacity,&
                             plotlib_supports_alpha,plot_text,plot_stbg,plot_ptxt,&
                             plot_get_clipping,plot_set_clipping
 use legends,           only:plot_box_around_text_xy
 use interactive_utils, only:get_posxy,get_vptxy,get_button_anchor,get_xw_margin_bounds
 character(len=*), intent(in) :: line1,line2
 logical, intent(in), optional :: onclick
 real :: xmin,xmax,ymin,ymax,x0,y0,x1,y1
 real :: xch,ych,oldch,fjust
 real, parameter :: xpos = 1.0 ! position in viewport coords
 real, parameter :: xpos_c = 0.5
 logical :: my_onclick
 logical, parameter :: centre_text = .true.
 integer :: oldclip
 real :: xpos_tmp,ypos,ym(2),ym_vp(2)

 my_onclick = .true.
 if (present(onclick)) my_onclick = onclick

 call plot_qch(oldch)
 call plot_qwin(xmin,xmax,ymin,ymax)
 call plot_get_clipping(oldclip)

 !call plot_stbg(0)
 call plot_sch(1.0)
 call get_button_anchor(xpos_tmp,ypos)

 ! get y bounds of X-window margin in world coordinates
 call get_xw_margin_bounds(ym_vp(1),ym_vp(2))
 call get_posxy(xpos_c,ym_vp(1),x1,ym(1),xmin,xmax,ymin,ymax)
 call get_posxy(xpos_c,ym_vp(2),x1,ym(2),xmin,xmax,ymin,ymax)

 ! shrink character height if it is too large to fit in the margin
 call plot_qcs(0,xch,ych)
 if (ych > 0.85*abs(ym_vp(2) - ym_vp(1))) then
    call plot_sch(0.85*abs(ym_vp(2) - ym_vp(1))/ych)
    call plot_qcs(0,xch,ych)
 endif

 if (plotlib_supports_alpha) call plot_set_opacity(0.25)
 call plot_set_clipping(0)

 ! line 2
 if (trim(line2) /= trim(help_msg) .or. my_onclick) then
    if (centre_text) then
       call get_posxy(xpos_c,ypos,x1,y1,xmin,xmax,ymin,ymax)
       fjust = 0.5
    else
       call get_posxy(xpos,ypos,x1,y1,xmin,xmax,ymin,ymax)
       fjust = 0.0
    endif

    ! erase the previous help message at the same anchor first
    if (len_trim(help_msg) > 0) then
       call plot_box_around_text_xy(x1,y1,fjust,1.0,trim(help_msg),ybounds=ym)
    endif

    ! block out pixels behind the new text by drawing a box
    call plot_box_around_text_xy(x1,y1,fjust,1.0,trim(line2),ybounds=ym)

    ! plot the help text
    if (len_trim(line2) > 0) call plot_ptxt(x1,y1,0.,fjust,trim(line2))

    ! save previous message
    help_msg = line2
 endif

 ! line 1
 call get_posxy(xpos,ypos,x0,y0,xmin,xmax,ymin,ymax)
 if (len_trim(last_msg) > 0) then
    call plot_box_around_text_xy(x0,y0,1.0,1.0,trim(last_msg),ybounds=ym)
 endif
 call plot_box_around_text_xy(x0,y0,1.0,1.0,trim(line1),ybounds=ym)
 if (len_trim(line1) > 0) call plot_ptxt(x0,y0,0.,1.0,trim(line1))
 last_msg = line1

 ! restore settings
 if (plotlib_supports_alpha) call plot_set_opacity(1.0)
 call plot_stbg(-1)
 call plot_sch(oldch)
 call plot_set_clipping(oldclip)

end subroutine print_message

!----------------------------------------------------
! print the value of a variable along with a message
!----------------------------------------------------
subroutine print_var(name,val,msg)
 character(len=*), intent(in) :: name,msg
 real, intent(in) :: val
 integer :: ierr
 character(len=64) :: string

 write(string,"(1x,a,1pg10.3)",iostat=ierr) name,val

 call print_message('',trim(string)//': '//msg)

end subroutine print_var

!----------------------------------------------------
! print x,y position along with a message
!----------------------------------------------------
subroutine print_message_xy(xpt,ypt,msg)
 use parsetext, only:number_to_string
 real, intent(in) :: xpt,ypt
 character(len=*), intent(in) :: msg
 character(len=20) :: strx,stry,strxy
 character(len=64) :: string

 ! plot x,y position as cursor moves
 call number_to_string(xpt,3,strx)
 call number_to_string(ypt,3,stry)
 strxy = '('//trim(strx)//', '//trim(stry)//')'

 string = trim(msg)
 call print_message(strxy,string)

end subroutine print_message_xy

!----------------------------------------------------
! print message with no xy box
!----------------------------------------------------
subroutine print_message_single(msg)
 character(len=*), intent(in) :: msg

 call print_message('',msg)

end subroutine print_message_single

!----------------------------------------------------
! print last saved help message
!----------------------------------------------------
subroutine print_last_message()

 call print_message('',help_msg)

end subroutine print_last_message

!----------------------------------------------------
! clear any saved help messages
!----------------------------------------------------
subroutine clear_messages()

 !help_msg = ''
 last_msg = ''

end subroutine clear_messages

end module interactive_help
