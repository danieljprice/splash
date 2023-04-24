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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! module which handles prompt for mapping of exact solution
! columns to the data columns in splash
!-----------------------------------------------------------------
module map_columns
 use asciiutils, only:match_integer
 implicit none
 integer :: nlab
 character(len=120), allocatable :: label_local(:)
 character(len=120), allocatable :: exact_label_local(:)
 integer, allocatable :: imap(:)

 procedure(print_map), pointer, private :: printmap => null()
 procedure(add_map),    pointer, private :: addmap => null()
 procedure(delete_map), pointer, private :: delmap => null()

contains

!-----------------------------------------------------------------
! interactive prompting for mapping list
!-----------------------------------------------------------------
subroutine map_columns_interactive(imapexact,label,exact_label,nlab_exact)
 use promptlist, only:prompt_list
 integer, intent(inout) :: imapexact(:)
 character(len=*), intent(in) :: label(:),exact_label(:)
 integer, intent(in) :: nlab_exact
 integer :: nmap

 ! allocate memory and copy input arrays
 imap = imapexact
 label_local = label
 exact_label_local = exact_label
 nlab = nlab_exact
 !maxlabels = size(label)
 if (nlab <= 2) return

 nmap = nlab

 ! prompt for edited map
 printmap => print_map
 addmap => add_map
 delmap => delete_map
 call prompt_list(nmap,nlab,'mapping',printmap,addmap,delmap)

 ! copy back to actual array
 imapexact = imap

 ! clean up
 if (allocated(label_local)) deallocate(label_local)
 if (allocated(exact_label_local)) deallocate(exact_label_local)

end subroutine map_columns_interactive

!-------------------------------------------
! print the current list of column mappings
!-------------------------------------------
subroutine print_map(nmap)
 integer, intent(in) :: nmap
 integer :: j,k,icol

 print "(/,a,i2)", ' Current mapping:'
 if (nmap==0) imap = 0 ! reset if nmap = 0, happens after clear operation

 do j=1,nlab
    if (imap(j) > 0) then
       print "(i2,': ',a12,a,i2,a)",j,exact_label_local(j),' -> ',imap(j),') '//trim(label_local(imap(j)))
    else
       print "(i2,': ',a)",j,trim(exact_label_local(j))
    endif
 enddo

end subroutine print_map

!------------------------------------------
! add new mapping
!------------------------------------------
subroutine add_map(istart,iend,nmap)
 use prompting,  only:prompt
 use asciiutils, only:match_tag_start
 integer, intent(in)    :: istart,iend
 integer, intent(inout) :: nmap
 integer :: i,icol

 if (istart == nmap) then
    !do i=1,nlab
   !    print "(i2,': ',a)",i,trim(exact_label_local(i))
    !enddo
    ! default prompt is the first column not already mapped
    icol = 1
    do while (imap(icol) > 0)
       icol = icol + 1
    enddo
    call prompt('which exact column to map?',icol,0,nlab)
 else
    icol = istart+1
 endif
 if (icol <= 0) return

 if (imap(icol)==0) then
    imap(icol) = match_tag_start(label_local,exact_label_local(icol))
 endif
 call prompt('which splash column to map '//trim(exact_label_local(icol))//' to?',imap(icol),0)

 nmap = nlab

end subroutine add_map

!------------------------------------------
! delete column mapping
!------------------------------------------
subroutine delete_map(icol,nmap)
 integer, intent(in)    :: icol
 integer, intent(inout) :: nmap

 if (icol > 0 .and. icol <= size(imap)) imap(icol) = 0
 print*,' deleting ',icol

 nmap = nlab

end subroutine delete_map

end module map_columns
