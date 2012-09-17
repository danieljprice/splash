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
!  Copyright (C) 2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module implementing generic menu where objects can be added
! to/deleted from a list
!-------------------------------------------------------------------------
module promptlist
 implicit none
 interface
  subroutine print_objlist(nobj)
   integer, intent(in) :: nobj
  end subroutine
 end interface
 
 interface
  subroutine add_obj(istart,iend,nobj)
   integer, intent(in) :: istart,iend
   integer, intent(inout) :: nobj
  end subroutine add_obj
 end interface
 
 interface
  subroutine delete_obj(iobj,nobj)
   integer, intent(in) :: iobj,nobj
  end subroutine delete_obj
 end interface

contains

subroutine prompt_list(nobj,maxobj,objname,print_objlist,add_obj,delete_obj)
 use prompting, only:prompt
 implicit none
 integer, intent(inout) :: nobj
 integer, intent(in)    :: maxobj
 character(len=*), intent(in) :: objname
 procedure(), pointer, intent(in) :: print_objlist,add_obj,delete_obj
 character(len=1) :: charp
 logical          :: done,first
 integer          :: istart,iend,ipick
 
 ipick = nobj + 1
 done  = .false.
 first = .true.
 charp = 'a'
 objmenu: do while(.not.done)
    call print_objlist(nobj)
    iend = maxobj
    if (nobj.gt.0 .or. .not.first) then
       charp='q'
       print*
       call prompt(' a)dd '//trim(objname)//', e)dit, d)elete, c)lear all or q)uit/finish?',&
                   charp,list=(/'a','e','d','c','q','s','S','Q'/),noblank=.true.)
       select case(charp)
       case('a')
          istart = nobj
          iend = nobj + 1
       case('e')
          if (nobj.gt.0) then
             ipick = 0
             call prompt(' pick a '//objname//' to edit ',ipick,0,nobj)
             if (ipick.gt.0) then
                istart = ipick - 1
                iend   = istart + 1
             else
                istart = 0
                iend = 1
                first = .false.
                cycle objmenu             
             endif
          else
             istart = 0
             iend   = 1
          endif
          first = .false.
       case('d')
          if (nobj.gt.0) then
             ipick = 0
             call prompt(' pick a '//objname//' to delete ',ipick,0,nobj)
             call delete_obj(ipick,nobj)
          else
             print*,'nothing to delete!'
          endif
          first = .false.
          cycle objmenu
       case('c')
          nobj = 0
          first = .false.
          cycle objmenu
       case('q','Q','s','S')
          done = .true.
       case default
          istart = 0
          iend = maxobj
       end select
    else
       istart = 0
       iend   = 1
    endif
    if (.not.done) call add_obj(istart,iend,nobj)
    first = .false.
 enddo objmenu

 return
end subroutine prompt_list

end module promptlist
