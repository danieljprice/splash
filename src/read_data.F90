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

!-------------------------------------------------------------
!
!  Moldule used to consolidate all of the read_data options
!
!-------------------------------------------------------------

module readdata
 ! This list is the reason why we really need a standard file format for SPH
 use readdata_sphNG,        only:read_data_sphNG,        set_labels_sphNG,   file_format_is_sphNG
 use readdata_ascii,        only:read_data_ascii,        set_labels_ascii
 use readdata_ndspmhd,      only:read_data_ndspmhd,      set_labels_ndspmhd, file_format_is_ndspmhd
 use readdata_gadget,       only:read_data_gadget,       set_labels_gadget,  file_format_is_gadget
 use readdata_VINE,         only:read_data_VINE,         set_labels_VINE
 use readdata_sro,          only:read_data_sro,          set_labels_sro
 use readdata_dragon,       only:read_data_dragon,       set_labels_dragon
 use readdata_seren,        only:read_data_seren,        set_labels_seren
 use readdata_tipsy,        only:read_data_tipsy,        set_labels_tipsy, file_format_is_tipsy
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
 use readdata_rsph,         only:read_data_rsph,         set_labels_rsph
 use readdata_vanaverbeke,  only:read_data_vanaverbeke,  set_labels_vanaverbeke
 use readdata_spyros,       only:read_data_spyros,       set_labels_spyros
 use readdata_urban,        only:read_data_urban,        set_labels_urban
 use readdata_starsmasher,  only:read_data_starsmasher,  set_labels_starsmasher
 use readdata_vtk,          only:read_data_vtk,          set_labels_vtk

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

 ! If the PBOB_DIR is given, then also include this
#ifdef PBOB_DIR
 use readdata_pbob,         only:read_data_pbob,         set_labels_pbob
#endif

#ifdef H5PART_DIR
 use readdata_h5part,       only:read_data_h5part,       set_labels_h5part
#endif

! use readdata_snsph,        only:read_data_snsph,        set_labels_snsph

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

subroutine select_data_format(string_in,ierr)
 use asciiutils,        only:lcase

 character(len=*),  intent(in)  :: string_in
 integer,           intent(out) :: ierr
 character(len=len(string_in)) :: string

 ! allow both -ascii and --ascii in flags
 if (string_in(1:1)=='-') then
    string = string_in(2:)
 else
    string = string_in
 endif

 !----------------------------
 !  Checks for dependencies
 !----------------------------

 ! Check if SPLASH has been compiled with hdf5, if hdf5 format is requested
 if ((index(string, 'hdf5') > 0) .or. (index(string, '.h5') > 0)) then
#ifndef HDF5
   print "(/,a)", ' *** ERROR: hdf5 file format requested, but SPLASH not compiled with HDF5 ***'
   print "(a)", '           Try make HDF5=yes   '
   stop
#else
   continue
#endif
 end if

 ! Check if H5PART is being requested
 if (string == 'h5part') then
#ifndef H5PART_DIR
   print "(/,a)", ' *** ERROR: H5PART file given, but SPLASH not compiled with the H5PART library. ***'
   print "(a)", '            Try make H5PART_DIR=/path/to/h5part/ Note: You must also compile with HDF5=yes'
   stop
#else
   continue
#endif
 end if

 ! Check if SPLASH has been compiled with the fits library, if fits format is requested
 if (string == 'fits') then
#ifndef FITS
   print "(/,a)", ' *** ERROR: .fits file given, but SPLASH not compiled with the fits library. ***'
   print "(a)", '            Try make FITS=yes '
#endif
 end if

 ! Check if PBOB is being requested
 if (string == 'pbob') then
#ifndef PBOB_DIR
   print "(/,a)", ' *** ERROR: .pbob file given, but SPLASH has not been compiled with the PBOB library. ***'
   print "(a)", '            Try make PBOB_DIR=/path/to/pbob/ '
#else
   continue
#endif
 end if

 ierr=0

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

 case('sro', 'srosph', 'magma')
   read_data=>read_data_sro
   set_labels=>set_labels_sro

 case('dragon')
   read_data=>read_data_dragon
   set_labels=>set_labels_dragon

 case('seren','gandalf')
   read_data=>read_data_seren
   set_labels=>set_labels_seren

 case('tipsy', 'gasoline')
   read_data=>read_data_tipsy
   set_labels=>set_labels_tipsy

 case('mhutch')
   read_data=>read_data_mhutch
   set_labels=>set_labels_mhutch

 case('ucla', 'ascii_ucla', 'ucla_ascii', 'sky')
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

 case('foulkes', 'steve')
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

 case('rsph')
   read_data=>read_data_rsph
   set_labels=>set_labels_rsph

! case('snsph', 'snsplash')
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

 case('gradsph', 'vanaverbeke', 'sigfried')
   read_data=>read_data_vanaverbeke
   set_labels=>set_labels_vanaverbeke

 case('starsmasher', 'star_smasher', 'jamiesph')
   read_data=> read_data_starsmasher
   set_labels=>set_labels_starsmasher

 case('vtk')
   read_data=> read_data_vtk
   set_labels=>set_labels_vtk

 ! Make the hdf5 data formats available if SPLASH has been compiled with HDF5
#ifdef HDF5
 case('phantom_hdf5', 'sphng_hdf5', 'phantomsph_hdf5')
   print "(a)",  'Phantom HDF5 files are currently not supported :( '
   stop

 case('gadget_hdf5','swift')
   read_data=>read_data_gadget_hdf5
   set_labels=>set_labels_gadget_hdf5

 case('amuse_hdf5')
   read_data=>read_data_amuse_hdf5
   set_labels=>set_labels_amuse_hdf5

! case('falcon_hdf5', 'falconhdf5', 'falcon')
!   read_data=>read_data_falcON_hdf5
!   set_labels=>set_labels_falcON_hdf5

 case('flash_hdf5', 'flashhdf5', 'flash')
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

 ! Make PBOB available if PBOB_DIR defined
#ifdef PBOB_DIR
 case('pbob')
   read_data=>read_data_pbob
   set_labels=>set_labels_pbob
#endif

 ! Make H5PART available of H5PART_DIR defined
#ifdef H5PART_DIR
 case('h5part')
   read_data=>read_data_h5part
   set_labels=>set_labels_h5part
#endif

 case default
   !call guess_format()
   ierr=1
 end select


end subroutine select_data_format

!----------------------------------------------------------------------
!  subroutine for printing help messages associated with reading data
!----------------------------------------------------------------------

subroutine print_available_formats(string)
 character(len=*), intent(in), optional :: string

 !print "(/,a)",' To select data formats, use the shortcuts below, or use the -f or --format command line options'
 !print "(a)"  ,' Multiple data formats are not support in a single instance.'
 if (string == 'short') then
    print "(/,a,/)",'Example data formats (type --formats for full list):'
 else
    print "(/,a,/)",'Supported data formats (auto = automatically identified from file contents):'
 endif
 print "(a)"  ,' -ascii,-csv          : ascii text/csv format (default)'
 print "(a)"  ,' -phantom -sphng      : Phantom and sphNG codes (auto)'
 print "(a)"  ,' -vtk                 : vtk legacy binary format (auto)'
 print "(a)"  ,' -ndspmhd             : ndspmhd code (auto)'
 print "(a)"  ,' -gandalf,-seren      : Gandalf/Seren code'
#ifdef HDF5
 print "(a)"  ,' -gadget -gadget_hdf5 : Gadget code (auto)'
 print "(a)"  ,' -falcon -falcon_hdf5 : FalcON code'
 print "(a)"  ,' -flash  -flash_hdf5  : FLASH code'
 print "(a)"  ,' -cactus -cactus_hdf5 : Cactus code'
 print "(a,/)",' -amuse  -amuse_hdf5  : AMUSE Framework'
#else
 print "(a)"  ,' -gadget              : Gadget code (auto)'
#endif
#ifdef FITS
 print "(a)"  ,' -fits                : FITS format (auto)'
#endif
#ifdef PBOB_DIR
 print "(a)"  ,' -pbob                : PBOB format'
#endif
#ifdef H5PART_DIR
 print "(a)"  ,' -h5part              : H5PART format'
#endif

 if (string /= 'short') then
    print "(a)",' -tipsy -gasoline     : Gasoline code (auto)'
    print "(a)",' -vine                : VINE SPH code'
    print "(a)",' -rsph                : Regularised SPH'
    print "(a)",' -starsmasher         : Star Smasher code'
    print "(a)",' -dragon              : DRAGON code'
    print "(a)",' -sro -magma          : Stephan Rosswog SPH code'
    print "(a)",' -gradsph             : GRADSPH code'
    print "(a)"," -mhutch              : Mark Hutchison's code"
    print "(a)"," -mbate               : Matthew Bate's code"
    print "(a)",' -oilonwater          : Oil-on-Water binary accretion SPH'
    print "(a)",' -ucla                : UCLA ascii format'
    print "(a)",' -urban               : Andrea Urban ascii format'
    print "(a)",' -spyros              : Spyros Kitsionas format'
    print "(a)",' -jjm  -jjmmulti      : Joe Monaghan format'
    print "(a)",' -bauswein            : Andreas Bauswein format'
    print "(a)",' -egaburov            : Evghenii Gaburov format'
    print "(a)",' -aly                 : Aly Reheam format'
    print "(a)",' -foulkes             : Foulkes ascii format'
    print "(a)",' -vanaverbeke         : Sigfried Vanaverbeke code'
    print "(a)",' -gadget_jsb          : GADGET Jamie Bolton variant'
#ifndef HDF5
    print "(a)",' -gadget_hdf5         : Gadget HDF5     [not compiled]'
    print "(a)",' -falcon -falcon_hdf5 : FalcON code     [not compiled]'
    print "(a)",' -flash  -flash_hdf5  : FLASH code      [not compiled]'
    print "(a)",' -cactus -cactus_hdf5 : Cactus code     [not compiled]'
    print "(a)",' -amuse  -amuse_hdf5  : AMUSE Framework [not compiled]'
#endif
#ifndef FITS
    print "(a)",' -fits                : FITS format     [not compiled]'
#endif
 end if

#ifdef HDF5
 print "(/,a)" ,'This build supports HDF5 formats. HDF5 files will be automatically'
 print "(a)",' recognised if they end with .h5, however you must specify a supported data format.'
 print "(a)",'  add a suffix "_hdf5" to above format if your data files do not end with .h5.'
#else
 print "(/,a)",'This build does not support HDF5. Compile with HDF5=yes to change this.'
#endif
#ifndef FITS
 print "(a)",'This build does not support FITS. Compile with FITS=yes to change this.'
#endif

end subroutine print_available_formats

!-----------------------------------------------------------------------------------
! subroutine for guessing the file format if not specified, or full info not given
!-----------------------------------------------------------------------------------
subroutine guess_format(nfiles,filenames,ierr,informat)
 use asciiutils, only:get_extensions
 integer, intent(in)                     :: nfiles
 character(len=*), intent(in)            :: filenames(:)
 character(len=*), intent(in), optional  :: informat ! This is given if --format <string> is supplied
 integer, intent(out)                    :: ierr
 character(len=12), dimension(5) :: extensions(5)

 call get_extensions(filenames(1), extensions)

 ierr = 0
 !
 ! try to guess the file format from the extension
 !
 if (any(index(extensions, '.h5') > 0) .or. any(index(extensions, '.hdf5') > 0)) then
    if (present(informat)) then
       call select_data_format(informat//"_hdf5",ierr)
    elseif (any((index(extensions, '.pb') > 0))) then
       call select_data_format("phantom_hdf5", ierr)
    else
       call select_data_format('gadget_hdf5',ierr)
    endif
 elseif (any((index(extensions, '.fits') > 0))) then
    call select_data_format('fits',ierr)
 elseif (any((index(extensions, '.vtk') > 0))) then
    call select_data_format('vtk',ierr)
 elseif (any((index(extensions, '.pb') > 0))) then
    call select_data_format('phantom', ierr)
 elseif (any((index(extensions, '.pbob') > 0))) then
    call select_data_format('pbob', ierr)
 elseif (any((index(extensions, '.csv') > 0))) then
    call select_data_format('ascii', ierr)
 else
    !
    ! if cannot guess from extension, then
    ! try to guess from the filename/header
    !
    call guess_format_from_file_header(filenames(1),ierr)
    !
    ! it is ok to not get a format, just assume ascii
    ! in this case we just return ierr /= 0
    !
 endif

end subroutine guess_format

!------------------------------------------------------------
! subroutine for guessing the file format from the filename
! and (if filename matches) the first few lines of the file
!------------------------------------------------------------
subroutine guess_format_from_file_header(filename,ierr)
 character(len=*), intent(in) :: filename
 integer, intent(out) :: ierr

 ierr = 1
 if (file_format_is_sphNG(filename)) then
    call select_data_format('sphNG',ierr)
 elseif (file_format_is_gadget(filename)) then
    call select_data_format('gadget',ierr)
 elseif (index(filename,'.hdf5') > 0) then
    call select_data_format('gadget_hdf5',ierr)
 elseif (file_format_is_ndspmhd(filename)) then
    call select_data_format('ndspmhd',ierr)
 elseif (file_format_is_tipsy(filename)) then
    call select_data_format('tipsy',ierr)
 endif

end subroutine guess_format_from_file_header

end module readdata
