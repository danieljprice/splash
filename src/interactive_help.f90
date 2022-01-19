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
 character(len=128) :: help_msg = ''

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
                             plotlib_supports_alpha,plot_text,plot_stbg
 use legends,           only:plot_box_around_text_xy
 use interactive_utils, only:get_posxy
 character(len=*), intent(in) :: line1,line2
 logical, intent(in), optional :: onclick
 real :: xmin,xmax,ymin,ymax,x0,y0,x1,y1
 real :: xch,ych,oldch
 real, parameter :: xpos = -0.01, ypos = -0.02 ! position in viewport coords
 logical :: my_onclick

 my_onclick = .true.
 if (present(onclick)) my_onclick = onclick

 call plot_qch(oldch)
 call plot_qwin(xmin,xmax,ymin,ymax)

 call plot_stbg(0)
 call plot_sch(1.0)
 call plot_qcs(0,xch,ych)
 if (plotlib_supports_alpha) call plot_set_opacity(0.25)

 ! line 2
 if (trim(line2) /= trim(help_msg) .or. my_onclick) then
    call get_posxy(xpos,ypos,x1,y1,xmin,xmax,ymin,ymax)

    ! block out pixels behind text by drawing a box
    call plot_box_around_text_xy(x1,y1,0.,1.0,line2)

    ! plot the help text
    if (len_trim(line2) > 0) call plot_text(x1,y1,line2)

    ! save previous message
    help_msg = line2
 endif

 ! line 1
 call get_posxy(xpos,ypos+1.5*ych,x0,y0,xmin,xmax,ymin,ymax)
 call plot_box_around_text_xy(x0,y0,0.,1.0,line1)
 if (len_trim(line1) > 0) call plot_text(x0,y0,line1)

 ! restore settings
 if (plotlib_supports_alpha) call plot_set_opacity(1.0)
 call plot_stbg(-1)
 call plot_sch(oldch)

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

 help_msg = ''

end subroutine clear_messages

end module interactive_help
