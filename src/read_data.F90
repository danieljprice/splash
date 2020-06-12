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
!  Copyright (C) 2005-2020 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------
!
!  Moldule used to consolidate all of the read_data options
!
!-------------------------------------------------------------

module readdata
 implicit none
 public :: select_data_format, guess_format, print_available_formats
 
 abstract interface
 subroutine set_labels_subroutine
 end subroutine
 end interface
 
 abstract interface
 subroutine read_data_subroutine(rootname,indexstart,iposn,nstepsread)
  integer,            intent(in)  :: indexstart,iposn
  integer,            intent(out) :: nstepsread
  character(len=*),   intent(in)  :: rootname
 
 end subroutine
 end interface
 
 procedure(read_data_subroutine),  pointer, public  :: read_data
 procedure(set_labels_subroutine), pointer, public  :: set_labels
 
 private
 
contains 

subroutine select_data_format(string,ierr)
 use readdata_sphNG,    only:read_data_sphNG,   set_labels_sphNG
 use readdata_ascii,    only:read_data_ascii,   set_labels_ascii
 use readdata_ndspmhd,  only:read_data_ndspmhd, set_labels_ndspmhd
 use readdata_gadget,   only:read_data_gadet,   set_labels_gadet
 use readdata_VINE,     only:read_data_VINE,    set_labels_VINE
 use readdata_sro,      only:read_data_sro,     set_labels_sro
 use readdata_dragon,   only:read_data_dragon,  set_labels_dragon
 use readdata_seren,    only:read_data_seren,   set_labels_seren
 use readdata_tipsy,    only:read_data_tipsy,   set_labels_tipsy
 use readdata_mhutch,   only:read_data_sro,     set_labels_mhutch
 
 
 use asciiutils,        only:lcase
 
 character(len=*),  intent(in)  :: string
 integer,           intent(out) :: ierr
 
 logical  :: selected_format
 
 selected_format = .false.
 
 ! Search input string for matching supported formats
 

 select case(trim(adjustl(lcase(string))))
 case('sphng', 'phantom')
   read_data=>read_data_sphNG
   set_labels=>set_labels_sphNG
   selected_format = .true.
 case('ascii')
   read_data=>read_data_ascii
   set_labels=>set_labels_ascii
   selected_format = .true.
 end select
 
 
 if (.not. selected_format) then
   print "(a)",' *** WARNING: file format '//trim(string)//' not found ***'
   ierr=1
 endif
 
end subroutine select_data_format

subroutine print_available_formats

 print "(a,/)",'Supported file formats:'
 print "(a)",' -ascii            : ascii file format (default)'
 print "(a)",' -phantom -sphng   : Phantom and sphNG '
 print "(a)"
 
end subroutine print_available_formats

subroutine guess_format(nfiles,filenames)
 integer, intent(in)          :: nfiles
 character(len=*), intent(in) :: filenames(:)

 print "(a)",' *** guess_format not yet implemented ***'

end subroutine guess_format

end module readdata
