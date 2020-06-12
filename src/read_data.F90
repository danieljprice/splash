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
 public :: select_data_format, guess_format
 public :: print_available_formats_short, print_available_formats
 
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
 ! This list is the reason why we really need a standard file format for SPH
 use readdata_sphNG,        only:read_data_sphNG,        set_labels_sphNG
 use readdata_ascii,        only:read_data_ascii,        set_labels_ascii
 use readdata_ndspmhd,      only:read_data_ndspmhd,      set_labels_ndspmhd
 use readdata_gadget,       only:read_data_gadget,       set_labels_gadget
 use readdata_VINE,         only:read_data_VINE,         set_labels_VINE
 use readdata_sro,          only:read_data_sro,          set_labels_sro
 use readdata_dragon,       only:read_data_dragon,       set_labels_dragon
 use readdata_seren,        only:read_data_seren,        set_labels_seren
 use readdata_tipsy,        only:read_data_tipsy,        set_labels_tipsy
 use readdata_mhutch,       only:read_data_mhutch,       set_labels_mhutch
 use readdata_UCLA,         only:read_data_UCLA,         set_labels_UCLA
 use readdata_aly,          only:read_data_aly,          set_labels_aly
 use readdata_bauswein,     only:read_data_bauswein,     set_labels_bauswein
 use readdata_dansph_old,   only:read_data_dansph_old,   set_labels_dansph_old
 use readdata_egaburov,     only:read_data_egaburov,     set_labels_egaburov
! !use readdata_fits,         only:read_data_fits,         set_labels_fits
 use readdata_foulkes,      only:read_data_foulkes,      set_labels_foulkes
 use readdata_gadget_jsb,   only:read_data_gadget_jsb,   set_labels_gadget_jsb
 use readdata_jjm,          only:read_data_jjm,          set_labels_jjm
 use readdata_jjmmulti,     only:read_data_jjmmulti,     set_labels_jjmmulti
 use readdata_kitp,         only:read_data_kitp,         set_labels_kitp
 use readdata_mbate,        only:read_data_mbate,        set_labels_mbate
 use readdata_mbate_hydro,  only:read_data_mbate_hydro,  set_labels_mbate_hydro
 use readdata_mbate_mhd,    only:read_data_mbate_mhd,    set_labels_mbate_mhd
 use readdata_oilonwater,   only:read_data_oilonwater,   set_labels_oilonwater
! use readdata_pbob,         only:read_data_pbob,         set_labels_pbob
 use readdata_rsph,         only:read_data_rsph,         set_labels_rsph
 use readdata_scw,         only:read_data_scw,         set_labels_scw
 use readdata_vanaverbeke,  only:read_data_vanaverbeke,  set_labels_vanaverbeke
! !use readdata_sphysics,     only:read_data_sphysics,     set_labels_sphysics
 use readdata_spyros,       only:read_data_spyros,       set_labels_spyros
 use readdata_urban,        only:read_data_urban,        set_labels_urban
 use readdata_jules,        only:read_data_jules,        set_labels_jules
! use readdata_snsph,        only:read_data_snsph,        set_labels_snsph

 use asciiutils,        only:lcase
 
 character(len=*),  intent(in)  :: string
 integer,           intent(out) :: ierr
 
 logical  :: selected_format
 
 selected_format = .false.
 
 !----------------------------------------------------
 ! Search input string for matching supported formats
 !----------------------------------------------------

 select case(trim(adjustl(lcase(string))))
 
 case('sphng', 'phantom', 'phantomsph')
   read_data=>read_data_sphNG
   set_labels=>set_labels_sphNG
   selected_format = .true.
   
 case('ascii')
   read_data=>read_data_ascii
   set_labels=>set_labels_ascii
   selected_format = .true.
 
 case('ndspmhd')
   read_data=>read_data_ndspmhd
   set_labels=>set_labels_ndspmhd
   selected_format = .true.
 
 case('gadget')
   read_data=>read_data_gadget
   set_labels=>set_labels_gadget
   selected_format = .true.
 
 case('vine')
   read_data=>read_data_VINE
   set_labels=>set_labels_VINE
   selected_format = .true.
  
 case('sro', 'srosph')
   read_data=>read_data_sro
   set_labels=>set_labels_sro
   selected_format = .true.
 
 case('dragon')
   read_data=>read_data_dragon
   set_labels=>set_labels_dragon
   selected_format = .true.
 
 case('seren')
   read_data=>read_data_seren
   set_labels=>set_labels_seren
   selected_format = .true.
 
 case('tipsy', 'gasoline')
   read_data=>read_data_tipsy
   set_labels=>set_labels_tipsy
   selected_format = .true.
  
 case('mhutch')
   read_data=>read_data_mhutch
   set_labels=>set_labels_mhutch
   selected_format = .true.
   
 case('ucla', 'ascii_ucla', 'ucla_ascii')
   read_data=>read_data_ucla
   set_labels=>set_labels_ucla
   selected_format = .true.
   
 case('aly')
   read_data=>read_data_aly
   set_labels=>set_labels_aly
   selected_format = .true.
   
 case('bauswein')
   read_data=>read_data_bauswein
   set_labels=>set_labels_bauswein
   selected_format = .true.
   
! case('fits')
!   read_data=>read_data_fits
!   set_labels=>set_labels_fits
!   selected_format = .true.
! 
 case('dansph_old')
   read_data=>read_data_dansph_old
   set_labels=>set_labels_dansph_old
   selected_format = .true.
   
 case('egaburov')
   read_data=>read_data_egaburov
   set_labels=>set_labels_egaburov
   selected_format = .true.
 
 case('foulkes')
   read_data=>read_data_foulkes
   set_labels=>set_labels_foulkes
   selected_format = .true.
 
 case('gadget_jsb', 'gadget-jsb')
   read_data=>read_data_gadget_jsb
   set_labels=>set_labels_gadget_jsb
   selected_format = .true.
 
 case('jjm')
   read_data=>read_data_jjm
   set_labels=>set_labels_jjm
   selected_format = .true.
   
 case('jjmmulti', 'jjm_multi', 'jjm_multiphase')
   read_data=>read_data_jjmmulti
   set_labels=>set_labels_jjmmulti
   selected_format = .true.
   
 case('jules')
   read_data=>read_data_jules
   set_labels=>set_labels_jules
   selected_format = .true.
!

 case('kitp')
   read_data=>read_data_kitp
   set_labels=>set_labels_kitp
   selected_format = .true.
   
 case('mbate')
   read_data=>read_data_mbate
   set_labels=>set_labels_mbate
   selected_format = .true.
   
 case('mbate_hydro', 'mbate-hydro')
   read_data=>read_data_mbate_hydro
   set_labels=>set_labels_mbate_hydro
   selected_format = .true.
   
 case('mbate_mhd', 'mbate-mhd')
   read_data=>read_data_mbate_mhd
   set_labels=>set_labels_mbate_mhd
   selected_format = .true.
   
 case('oilonwater')
   read_data=>read_data_oilonwater
   set_labels=>set_labels_oilonwater
   selected_format = .true.
   
! case('pbob')
!   read_data=>read_data_pbob
!   set_labels=>set_labels_pbob
!   selected_format = .true.
   
 case('rsph')
   read_data=>read_data_rsph
   set_labels=>set_labels_rsph
   selected_format = .true.

 case('scw')
   read_data=>read_data_scw
   set_labels=>set_labels_scw
   selected_format = .true.
      
! case('snsph')
!   read_data=>read_data_snsph
!   set_labels=>set_labels_snsph
!   selected_format = .true.
!   
! case('sphysics')
!   read_data=>read_data_sphysics
!   set_labels=>set_labels_sphysics
!   selected_format = .true.
!  
 case('spyros')
   read_data=>read_data_spyros
   set_labels=>set_labels_spyros
   selected_format = .true.
    
 case('urban')
   read_data=>read_data_urban
   set_labels=>set_labels_urban
   selected_format = .true.
    
 case('vanaverbeke')
   read_data=>read_data_vanaverbeke
   set_labels=>set_labels_vanaverbeke
   selected_format = .true.
 end select
 
 
 if (.not. selected_format) then
   print "(a)",' *** WARNING: file format '//trim(string)//' not found ***'
   ierr=1
 endif
 
end subroutine select_data_format

subroutine print_available_formats_short

 print "(a)"
 print "(a)",'To select data formats, use the short cuts below, or'
 print "(a)",'use the -f or --format command line options'
 print "(a,/)",'Supported file formats:'
 print "(a)",' -ascii            : ascii file format (default)'
 print "(a)",' -phantom -sphng   : Phantom and sphNG '
 print "(a)",' -ndspmhd          : ndsphmd '
 print "(a)",' -gadet            : Gadget  '
 print "(a)",' -seren            : Seren '
 print "(a)",' ..plus many others. Type --help for a full list '
 print "(a)"
 
end subroutine print_available_formats_short


subroutine print_available_formats
 
 print "(a)"
 print "(a)",'To select data formats, use the short cuts below, or'
 print "(a)",'use the -f or --format command line options'
 print "(a,/)",'Supported file formats:'
 print "(a)",' -ascii            : ascii file format (default)'
 print "(a)",' -phantom -sphng   : Phantom and sphNG '
 print "(a)",' -ndspmhd          : ndsphmd '
 print "(a)",' -gadet            : Gadget  '
 print "(a)",' -seren            : Seren '
 print "(a)"
 
end subroutine print_available_formats

subroutine guess_format(nfiles,filenames)
 integer, intent(in)          :: nfiles
 character(len=*), intent(in) :: filenames(:)

 print "(a)",' *** guess_format not yet implemented ***'

end subroutine guess_format

end module readdata
