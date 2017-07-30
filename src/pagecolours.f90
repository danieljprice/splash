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

!----------------------------------------------------------
!  write line colours to file
!----------------------------------------------------------
subroutine write_coloursfile(filename,nc,linecolours)
 character(len=*), intent(in) :: filename
 integer, intent(in) :: nc
 real,    intent(in) :: linecolours(3,nc)
 integer :: i,ierr
 integer, parameter :: lu = 37
 
 open(unit=lu,file=trim(filename),status='replace',iostat=ierr)
 do i=1,nc
    write(lu,"(3(f8.3,1x))") linecolours(:,i)
 enddo
 close(lu)

end subroutine write_coloursfile

!----------------------------------------------------------
!  read line colours from file
!----------------------------------------------------------
subroutine read_coloursfile(filename,maxc,linecolours,nc,ierr)
 character(len=*), intent(in) :: filename
 integer, intent(in)  :: maxc
 real,    intent(out) :: linecolours(3,maxc)
 integer, intent(out) :: nc,ierr
 integer :: i
 integer, parameter :: lu = 38
 
 ierr = 0
 nc = 0
 open(unit=lu,file=trim(filename),status='old',iostat=ierr)
 if (ierr /= 0) return
 do i=1,maxc
    read(lu,*,iostat=ierr) linecolours(:,i)
    if (ierr==0) nc = nc + 1
 enddo
 close(lu)
 
end subroutine read_coloursfile

!----------------------------------------------------------
!  read line colours from file
!----------------------------------------------------------
subroutine set_linecolours(linepalette,filename,maxc,linecolours)
 use plotlib, only:plot_scr,plot_set_palette,plotlib_maxpalette
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: linepalette,maxc
 real,             intent(out) :: linecolours(3,maxc)
 integer :: i,nc,ierr

 if (linepalette > 0 .and. linepalette <= plotlib_maxpalette) then
    call plot_set_palette(linepalette)
 elseif (linepalette < 0) then 
    call read_coloursfile(filename,maxc,linecolours,nc,ierr)
    do i=1,nc
       call plot_scr(i+1,linecolours(1,i),linecolours(2,i),linecolours(3,i))
    enddo
 endif
 
end subroutine set_linecolours

end module pagecolours
