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
!  Copyright (C) 2005-2010 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!----------------------------------------------------------
!  module handling page colour schemes, for generic
!  plotting library
!----------------------------------------------------------
module pagecolours
 implicit none
 integer, parameter :: maxpagecolours = 2

contains

!----------------------------------------------------------
!  query function for the colour scheme
!----------------------------------------------------------
function pagecolourscheme(ischeme,short)
 implicit none
 integer, intent(in)            :: ischeme
 character(len=21)              :: pagecolourscheme
 logical, intent(in), optional  :: short
 logical :: use_short

 use_short = .false.
 if (present(short)) use_short = short

 select case(ischeme)
 case(2)
    pagecolourscheme = 'white-on-black'
 case(1)
    pagecolourscheme = 'black-on-white'
 case default
    if (use_short) then
       pagecolourscheme = 'default'
    else
       pagecolourscheme = 'plot library default'
    endif
 end select

end function pagecolourscheme

!----------------------------------------------------------
!  set the colour index 1 and 0 of the plotting library
!  corresponding to the foreground and background colours
!  (must be called after the plot library is initialised)
!----------------------------------------------------------
subroutine set_pagecolours(ischeme)
 use plotlib, only:plot_scr
 implicit none
 integer, intent(in) :: ischeme

 select case(ischeme)
 case(2) !--white-on-black
    call plot_scr(0,0.,0.,0.)
    call plot_scr(1,1.,1.,1.)
 case(1) !--black-on-white
    call plot_scr(0,1.,1.,1.)
    call plot_scr(1,0.,0.,0.)
 end select

end subroutine set_pagecolours

!----------------------------------------------------------
!  query function for the name of the foreground colour
!----------------------------------------------------------
function colour_fore(ischeme)
 implicit none
 integer, intent(in) :: ischeme
 character(len=5)    :: colour_fore

 select case(ischeme)
 case(2)
    colour_fore = 'white'
 case(1)
    colour_fore = 'black'
 case default
    colour_fore = ' '
 end select

end function colour_fore

!----------------------------------------------------------
!  query function for the name of the background colour
!----------------------------------------------------------
function colour_back(ischeme)
 implicit none
 integer, intent(in) :: ischeme
 character(len=5)    :: colour_back

 select case(ischeme)
 case(2)
    colour_back = 'black'
 case(1)
    colour_back = 'white'
 case default
    colour_back = ' '
 end select

end function colour_back

end module pagecolours
