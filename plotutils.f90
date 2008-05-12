!---------------------------------------------------------------------------
! module containing application programming interfaces for basic
! plotting functions. The idea is to add more to this module to
! eventually use it to be able to change backends more easily.
!---------------------------------------------------------------------------
module plotutils
 implicit none
 public :: plotline,plotbins
 
 private

contains

!
!  line plotting, with blanking
!
subroutine plotline(npts,xline,yline,blank)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(:) :: xline,yline
 real, intent(in), optional :: blank
 integer :: i,nseg,istart
 
 if (present(blank)) then
    nseg = 0
    istart = 1
    !--plot line in segments, leaving blank segments where y=blank
    do i=1,npts
       if (abs(yline(i)-blank).lt.tiny(yline) .or. i.eq.npts) then
          if (nseg.gt.0) call pgline(nseg,xline(istart:istart+nseg),yline(istart:istart+nseg))
          istart = i+1
          nseg = 0
       else
          nseg = min(nseg + 1,npts-1)
       endif
    enddo
 else
    call pgline(npts,xline,yline) 
 endif

 return
end subroutine plotline

!
!  binned histogram plotting, with blanking
!
subroutine plotbins(nbins,xbins,ybins,blank)
 implicit none
 integer, intent(in) :: nbins
 real, intent(in), dimension(:) :: xbins,ybins
 real, intent(in), optional :: blank
 integer :: i,nseg,istart
 
 if (present(blank)) then
    nseg = 0
    istart = 1
    !--plot line in segments, leaving blank segments where y=blank
    do i=1,nbins
       if (abs(ybins(i)-blank).lt.tiny(ybins) .or. i.eq.nbins) then
          if (nseg.gt.0) call pgbin(nseg,xbins(istart:istart+nseg),ybins(istart:istart+nseg),.true.)
          istart = i+1
          nseg = 0
       else
          nseg = min(nseg + 1,nbins-1)
       endif
    enddo
 else
    call pgbin(nbins,xbins,ybins,.true.)
 endif
 
 return
end subroutine plotbins

end module plotutils
