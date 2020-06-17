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
 public :: print_available_formats
 
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

!----------------------------------------------------------------------
! subroutine to select string specified data format 
!----------------------------------------------------------------------

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
 use readdata_egaburov,     only:read_data_egaburov,     set_labels_egaburov
 use readdata_foulkes,      only:read_data_foulkes,      set_labels_foulkes
 use readdata_gadget_jsb,   only:read_data_gadget_jsb,   set_labels_gadget_jsb
 use readdata_jjm,          only:read_data_jjm,          set_labels_jjm
 use readdata_jjmmulti,     only:read_data_jjmmulti,     set_labels_jjmmulti
 use readdata_mbate,        only:read_data_mbate,        set_labels_mbate
 use readdata_oilonwater,   only:read_data_oilonwater,   set_labels_oilonwater
! use readdata_pbob,         only:read_data_pbob,         set_labels_pbob
 use readdata_rsph,         only:read_data_rsph,         set_labels_rsph
 use readdata_vanaverbeke,  only:read_data_vanaverbeke,  set_labels_vanaverbeke
! !use readdata_sphysics,     only:read_data_sphysics,     set_labels_sphysics
 use readdata_spyros,       only:read_data_spyros,       set_labels_spyros
 use readdata_urban,        only:read_data_urban,        set_labels_urban
! use readdata_snsph,        only:read_data_snsph,        set_labels_snsph

 ! Make hdf5 fortran/c modules available if compiled with hdf5
#ifdef HDF5
 use readdata_amuse_hdf5,   only:read_data_amuse_hdf5,   set_labels_amuse_hdf5
 use readdata_cactus_hdf5,  only:read_data_cactus_hdf5,  set_labels_cactus_hdf5
! use readdata_falcON_hdf5,  only:read_data_falcON_hdf5,  set_labels_falcON_hdf5
 use readdata_flash_hdf5,   only:read_data_flash_hdf5,   set_labels_flash_hdf5
 use readdata_gadget_hdf5,  only:read_data_gadget_hdf5,  set_labels_gadget_hdf5
#endif

 ! Same for FITS files
#ifdef FITS
 use readdata_fits,         only:read_data_fits,         set_labels_fits
#endif

 use asciiutils,        only:lcase
 
 character(len=*),  intent(in)  :: string
 integer,           intent(out) :: ierr
 
 
 ! Check if SPLASH has been compiled with hdf5, if hdf5 format is requested
 if ((index(string, 'hdf5') > 0) .or. (index(string, '.h5') > 0)) then
#ifndef HDF5
   print "(a)", '*** ERROR: hdf5 file format requested, but SPLASH'
   print "(a)", '           has not been compiled with HDF5 ***   '
   stop
#else
   continue
#endif
 end if

 ! Check if SPLASH has been compiled with the fits library, if fits format is requested
 if (string == 'fits') then
#ifndef FITS
   print "(a)", '*** ERROR: .fits file given, but SPLASH has not been compiled'
   print "(a)", '           with the fits library. ***'
   stop
#else
   continue
#endif
 end if

 !----------------------------------------------------
 ! Search input string for matching supported formats
 !----------------------------------------------------

 select case(trim(adjustl(lcase(string))))
 
 case('sphng', 'phantom', 'phantomsph')
   read_data=>read_data_sphNG
   set_labels=>set_labels_sphNG

 case('ascii')
   read_data=>read_data_ascii
   set_labels=>set_labels_ascii

 case('ndspmhd')
   read_data=>read_data_ndspmhd
   set_labels=>set_labels_ndspmhd

 case('gadget')
   read_data=>read_data_gadget
   set_labels=>set_labels_gadget

 case('vine')
   read_data=>read_data_VINE
   set_labels=>set_labels_VINE

 case('sro', 'srosph')
   read_data=>read_data_sro
   set_labels=>set_labels_sro

 case('dragon')
   read_data=>read_data_dragon
   set_labels=>set_labels_dragon

 case('seren')
   read_data=>read_data_seren
   set_labels=>set_labels_seren

 case('tipsy', 'gasoline')
   read_data=>read_data_tipsy
   set_labels=>set_labels_tipsy

 case('mhutch')
   read_data=>read_data_mhutch
   set_labels=>set_labels_mhutch

 case('ucla', 'ascii_ucla', 'ucla_ascii')
   read_data=>read_data_ucla
   set_labels=>set_labels_ucla

 case('aly')
   read_data=>read_data_aly
   set_labels=>set_labels_aly

 case('bauswein')
   read_data=>read_data_bauswein
   set_labels=>set_labels_bauswein

 case('egaburov')
   read_data=>read_data_egaburov
   set_labels=>set_labels_egaburov

 case('foulkes')
   read_data=>read_data_foulkes
   set_labels=>set_labels_foulkes

 case('gadget_jsb', 'gadget-jsb')
   read_data=>read_data_gadget_jsb
   set_labels=>set_labels_gadget_jsb

 case('jjm')
   read_data=>read_data_jjm
   set_labels=>set_labels_jjm

 case('jjmmulti', 'jjm_multi', 'jjm_multiphase')
   read_data=>read_data_jjmmulti
   set_labels=>set_labels_jjmmulti
   
 case('mbate')
   read_data=>read_data_mbate
   set_labels=>set_labels_mbate

 case('oilonwater')
   read_data=>read_data_oilonwater
   set_labels=>set_labels_oilonwater

! case('pbob')
!   read_data=>read_data_pbob
!   set_labels=>set_labels_pbob

 case('rsph')
   read_data=>read_data_rsph
   set_labels=>set_labels_rsph

! case('snsph')
!   read_data=>read_data_snsph
!   set_labels=>set_labels_snsph

! case('sphysics')
!   read_data=>read_data_sphysics
!   set_labels=>set_labels_sphysics

 case('spyros')
   read_data=>read_data_spyros
   set_labels=>set_labels_spyros

 case('urban')
   read_data=>read_data_urban
   set_labels=>set_labels_urban

 case('vanaverbeke')
   read_data=>read_data_vanaverbeke
   set_labels=>set_labels_vanaverbeke

 ! Make the hdf5 data formats available if SPLASH has been compiled with HDF5
#ifdef HDF5
 case('phantom_hdf5', 'sphng_hdf5', 'phantomsph_hdf5')
   print "(a)",  'Phantom HDF5 files are currently not supported :( '
   stop

 case('gadget_hdf5')
   read_data=>read_data_gadget_hdf5
   set_labels=>set_labels_gadget_hdf5
   
 case('amuse_hdf5')
   read_data=>read_data_amuse_hdf5
   set_labels=>set_labels_amuse_hdf5
   
! case('falcon_hdf5')
!   read_data=>read_data_falcON_hdf5
!   set_labels=>set_labels_falcON_hdf5
   
 case('flash_hdf5')
   read_data=>read_data_flash_hdf5
   set_labels=>set_labels_flash_hdf5
 
 case('cactus', 'cactus_hdf5')
   read_data=>read_data_cactus_hdf5
   set_labels=>set_labels_cactus_hdf5
#endif

 ! Make fits routines available if SPLASH has been compiled with FITS
#ifdef FITS
 case('fits')
   read_data=>read_data_fits
   set_labels=>set_labels_fits
#endif

 case default
   print "(a)",' *** WARNING: file format '//trim(string)//' not found ***'
   ierr=1
 end select
 

 
 

 
end subroutine select_data_format

!----------------------------------------------------------------------
!  subroutine for printing help messages associated with reading data
!----------------------------------------------------------------------

subroutine print_available_formats(string)
 character(len=*), intent(in), optional :: string

 print "(/,a)",'To select data formats, use the short cuts below, or'
 print "(a)"  ,'use the -f or --format command line options.'
 print "(a)"  ,'Multiple data formats are not support in a single instance.'
 print "(a,/)",'Supported data formats:'
 print "(a)"  ,' -ascii            : ascii file format (default)'
 print "(a)"  ,' -phantom -sphng   : Phantom and sphNG codes'
 print "(a)"  ,' -ndspmhd          : ndsphmd code'
 print "(a)"  ,' -gadget           : Gadget  code'
 print "(a)"  ,' -seren            : Seren code'
 
 if (string=='short') then
   print "(a,/)",' ..plus many others. Type --help for a full list '
 else
   print "(a)"  ,' -flash            : FLASH code'
   print "(a)"  ,' -tispy -gasoline  : Gasoline code'
   print "(a)"  ,' -ucla             : UCLA ascii format'
 end if

#ifdef HDF5
 print "(a)"  ,'This build of SPLASH supports the HDF5 file format.'
 print "(a)"  ,'HDF5 files will be automatically recognised if they'
 print "(a)"  ,'end with .h5, however you must specify a supported'
 print "(a)"  ,'data format from above.'
#else
 print "(a)"  ,'This build of SPLASH does not support HDF5. '
#endif

end subroutine print_available_formats

!-----------------------------------------------------------------------------------
! subroutine for guessing the file format if not specified, or full info not given
!-----------------------------------------------------------------------------------

subroutine guess_format(nfiles,filenames,ierr,string)
 integer, intent(in)                     :: nfiles
 character(len=*), intent(in)            :: filenames(:)
 character(len=*), intent(in), optional  :: string
 integer, intent(out)                    :: ierr
 
 character(len=5), dimension(3), parameter :: extensions = &
           (/'.fits','.h5  ','.pb  '/)
 
 logical    :: selected_format
 integer    :: i
 
 
 selected_format = .false.

 do i = 1, size(extensions)
   
   if (any((index(filenames, extensions(i)) > 0))) then
     call select_data_format('fits',ierr)
   end if
   
 end do
 
end subroutine guess_format

end module readdata
