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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module sort
 implicit none
 public :: indexx
 private

contains

subroutine indexx(n, arr, indx)
!************************************************************
!                                                           *
!  This is INDEXX using the quicksort algorithm.            *
!                                                           *
!************************************************************
 implicit none
 integer, parameter :: m=7, nstack=500
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: arr
 integer, dimension(n), intent(out) :: indx

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer, dimension(nstack) :: istack
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l.lt.m) then
     do j = l + 1, ir
        indxt = indx(j)
        a = arr(indxt)
        do i = j - 1, 1, -1
           if (arr(indx(i)).le.a) goto 2
           indx(i + 1) = indx(i)
        end do
        i = 0
2       indx(i + 1) = indxt
     end do
     if (jstack.eq.0) return
     ir = istack(jstack)
     l = istack(jstack - 1)
     jstack = jstack - 2
  else
     k = (l + ir)/2
     itemp = indx(k)
     indx(k) = indx(l + 1)
     indx(l + 1) = itemp
     if (arr(indx(l + 1)).gt.arr(indx(ir))) then
        itemp = indx(l + 1)
        indx(l + 1) = indx(ir)
        indx(ir) = itemp
     endif
     if (arr(indx(l)).gt.arr(indx(ir))) then
        itemp = indx(l)
        indx(l) = indx(ir)
        indx(ir) = itemp
     endif
     if (arr(indx(l + 1)).gt.arr(indx(l))) then
        itemp = indx(l + 1)
        indx(l + 1) = indx(l)
        indx(l) = itemp
     endif
     i = l + 1
     j = ir
     indxt = indx(l)
     a = arr(indxt)
3    continue
     i = i + 1
     if (arr(indx(i)).lt.a) goto 3
4    continue
     j = j - 1
     if (arr(indx(j)).gt.a) goto 4
     if (j.lt.i) goto 5
     itemp = indx(i)
     indx(i) = indx(j)
     indx(j) = itemp
     goto 3

5    indx(l) = indx(j)
     indx(j) = indxt
     jstack = jstack + 2
     if (jstack.gt.nstack) then
        print*,'fatal error!!! stacksize exceeded in sort'
        print*,'(need to set parameter nstack higher in subroutine indexx '
        print*,' this is in the file sort.f90)'
        stop
     endif
     if (ir - i + 1.ge.j - l) then
        istack(jstack) = ir
        istack(jstack - 1) = i
        ir = j - 1
     else
        istack(jstack) = j - 1
        istack(jstack - 1) = l
        l = i
     endif
  endif

goto 1
end subroutine indexx

end module sort
