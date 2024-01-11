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

 public :: map_labels,map_columns_in_file,map_columns_interactive
 public :: print_mapping

 private

contains

!---------------------------------------------------
! map labels in one list onto corresponding labels
! in another list
!---------------------------------------------------
function map_labels(list1,list2) result(imap)
 use asciiutils, only:match_lists
 use labels,     only:label_synonym
 character(len=*), intent(in) :: list1(:),list2(:)
 integer :: imap(size(list1)),i

 imap = match_lists(label_synonym(list1),label_synonym(list2))

 ! make mapping unique
 do i=2,size(imap)
    if (any(imap(1:i-1)==imap(i))) imap(i) = 0
 enddo

end function map_labels

!---------------------------------------------------
! map labels in the ascii file read as the exact
! solution onto corresponding labels in the data
! (i.e. find which columns correspond to which)
!---------------------------------------------------
subroutine map_columns_in_file(filename,ncols,nrows,imap,label,exact_labels,nlabels,ierr)
 use asciiutils, only:get_ncolumns,get_nrows,read_column_labels
 use labels,     only:set_default_labels
 character(len=*), intent(in) :: filename
 integer, intent(out) :: ncols,nrows,nlabels
 integer, intent(out) :: imap(:), ierr
 character(len=*), intent(in)  :: label(:)
 character(len=*), intent(out) :: exact_labels(:)
 integer :: nheaderlines,ncolumns,iu

 ncolumns = size(label)
 nlabels = 0; ncols = 0; nrows = 0
 call set_default_labels(exact_labels(:))

 open(newunit=iu,file=filename,status='old',iostat=ierr)
 if (ierr == 0) then
    call get_ncolumns(iu,ncols,nheaderlines)
    call get_nrows(iu,nheaderlines,nrows)
    call read_column_labels(iu,nheaderlines,ncols,nlabels,exact_labels) !,debug=.true.
    close(iu)
    !if (nlabels > 0) print "(/,1x,a,2(i0,a),/)",trim(filename)//': ',nrows,' lines, ',nlabels,' columns'
 endif
 if (nlabels > 0) then
    nlabels = min(nlabels,size(imap),size(exact_labels),ncols)
    imap(1:nlabels) = map_labels(exact_labels(1:nlabels),label(1:ncolumns))
 else
    nlabels = min(ncols,size(exact_labels))
    imap = 0
 endif

end subroutine map_columns_in_file

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
 nlab = min(nlab_exact,size(imap))
 !maxlabels = size(label)
 if (nlab < 2) return

 nmap = min(nlab,size(imap))

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
subroutine print_mapping(ncols1,imap1,list1,list2)
 integer, intent(in) :: ncols1
 integer, intent(in) :: imap1(ncols1)
 character(len=*), intent(in) :: list1(:),list2(:)
 integer :: j

 print "(/,a)", ' Current mapping:'

 do j=1,ncols1
    if (imap1(j) > 0 .and. imap1(j) <= size(list2)) then
       print "(i3,': ',a12,a,i2,a)",j,list1(j),' -> ',imap1(j),') '//trim(list2(imap1(j)))
    else
       print "(i3,': ',a12,a)",j,list1(j),' -> None'
    endif
 enddo

end subroutine print_mapping

!-------------------------------------------
! print the current list of column mappings
!-------------------------------------------
subroutine print_map(nmap)
 integer, intent(in) :: nmap

 if (nmap==0) imap = 0 ! reset if nmap = 0, happens after clear operation
 call print_mapping(nlab,imap,exact_label_local,label_local)

end subroutine print_map

!------------------------------------------
! add new mapping
!------------------------------------------
subroutine add_map(istart,iend,nmap)
 use prompting,  only:prompt
 use asciiutils, only:match_tag_start
 use labels,     only:label_synonym
 integer, intent(in)    :: istart,iend
 integer, intent(inout) :: nmap
 integer :: icol

 if (istart == nmap) then
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
    imap(icol) = match_tag_start(label_synonym(label_local),&
                                 label_synonym(exact_label_local(icol)))
 endif
 call prompt('which splash column to map '//trim(exact_label_local(icol))// &
             ' to?',imap(icol),0,size(label_local))

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
