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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
! module containing application programming interfaces for basic
! plotting functions. The idea is to add more to this module to
! eventually use it to be able to change backends more easily.
!---------------------------------------------------------------------------
module plotutils
 use plotlib, only:plot_line,plot_bins
 implicit none
 public :: plotline,plotbins,formatreal

 private

contains

!
!  line plotting, with blanking
!
subroutine plotline(npts,xline,yline,blank)
 integer, intent(in) :: npts
 real, intent(in), dimension(:) :: xline,yline
 real, intent(in), optional :: blank
 integer :: i,nseg,istart

 if (present(blank)) then
    nseg = 0
    istart = 1
    !--plot line in segments, leaving blank segments where y=blank
    do i=1,npts
       if (abs(yline(i)-blank) < tiny(yline) .or. i==npts) then
          if (nseg > 0) call plot_line(nseg,xline(istart:istart+nseg),yline(istart:istart+nseg))
          istart = i+1
          nseg = 0
       else
          nseg = min(nseg + 1,npts-1)
       endif
    enddo
 else
    call plot_line(npts,xline,yline)
 endif

 return
end subroutine plotline

!
!  binned histogram plotting, with blanking
!
subroutine plotbins(nbins,xbins,ybins,blank)
 integer, intent(in) :: nbins
 real, intent(in), dimension(:) :: xbins,ybins
 real, intent(in), optional :: blank
 integer :: i,nseg,istart

 if (present(blank)) then
    nseg = 0
    istart = 1
    !--plot line in segments, leaving blank segments where y=blank
    do i=1,nbins
       if (abs(ybins(i)-blank) < tiny(ybins) .or. i==nbins) then
          if (nseg > 0) call plot_bins(nseg,xbins(istart:istart+nseg),ybins(istart:istart+nseg),.true.)
          istart = i+1
          nseg = 0
       else
          nseg = min(nseg + 1,nbins-1)
       endif
    enddo
 else
    call plot_bins(nbins,xbins,ybins,.true.)
 endif

 return
end subroutine plotbins

!
!  formatting of real variables into strings (like PGNUMB)
!
subroutine formatreal(val,string,ierror)
 real, intent(in) :: val
 character(len=*), intent(out) :: string
 integer, intent(out), optional :: ierror
 integer :: ierr,i,idot
 logical :: nonzero

 if (abs(val) >= 1.d99) then
    write(string,"(1pe10.3)",iostat=ierr) val
 elseif (abs(val) < 1.e-3 .or. abs(val) >= 1.e4) then
    write(string,"(1pe9.2)",iostat=ierr) val
 elseif (abs(val) < 0.1) then
    write(string,"(f8.3)",iostat=ierr) val
 elseif (abs(val) >= 100.) then
    write(string,"(f8.0)",iostat=ierr) val
 else
    write(string,"(f8.2)",iostat=ierr) val
 endif
 string = adjustl(trim(string))

 if (present(ierror)) ierror = ierr

 !
 !--strip trailing zeros after the decimal place
 !  (and the decimal place if it is the last character)
 !
 idot = index(string,'.')
 if (idot > 0) then
    nonzero = .false.
    do i = len_trim(string),idot,-1
       if (.not.nonzero .and. string(i:i)=='0') then
          string(i:i) = ' '
       elseif (.not.nonzero .and. string(i:i)=='.') then
          string(i:i) = ' '
          nonzero = .true.
       else
          nonzero = .true.
       endif
    enddo
 endif
 string = trim(string)

 return
end subroutine formatreal

end module plotutils
