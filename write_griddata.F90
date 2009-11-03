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
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Module implementing "splash to grid" operation, writing
! 3D gridded data in various output formats
!-----------------------------------------------------------------
module write_griddata
 implicit none

 public :: isgridformat,print_gridformats
 public :: open_gridfile,write_grid
 private

contains

!-----------------------------------------------------------------
! utility to check if a format selection is valid
!-----------------------------------------------------------------
logical function isgridformat(string)
 use asciiutils, only:lcase
 implicit none
 character(len=*), intent(in) :: string

 isgridformat = .false.
 select case(trim(lcase(string)))
 case('grid')
     isgridformat = .true.
 case('gridascii')
     isgridformat = .true.
 case('gridbinary')
     isgridformat = .true.
 end select

end function isgridformat

!-----------------------------------------------------------------
! print usage if format selection not valid
!-----------------------------------------------------------------
subroutine print_gridformats
 implicit none

 print "(/,a)",' grid conversion mode ("splash to X dumpfiles"): '
 print "(a)",'    splash to grid         : interpolate basic SPH data (density, plus velocity if present in data)'
 print "(a)",'                             to 3D grid, write grid data to file (using default output=ascii)'
 print "(a)",'           to gridascii    : as above, grid data written in ascii format'
 print "(a)",'           to gridbinary   : as above, grid data in simple unformatted binary format'
 print "(a)",'        allto grid         : as above, interpolating *all* columns to the grid (and output file)'
 print "(a)",'        allto gridascii    : as above, with ascii output'
 print "(a)",'        allto gridbinary   : as above, with binary output'
 
 return
end subroutine print_gridformats

!------------------------------------------------------
! open grid file for output, write header
!------------------------------------------------------
subroutine open_gridfile(iunit,filenamein,outformat,npixels,ncolumns,time,ierr)
 use asciiutils, only:lcase
 implicit none
 integer, intent(in)               :: iunit
 character(len=*), intent(in)      :: filenamein,outformat
 character(len=len(filenamein)+10) :: filename
 integer, dimension(3), intent(in) :: npixels
 integer, intent(in)               :: ncolumns
 real, intent(in)                  :: time
 integer, intent(out)              :: ierr
!  
!--Only have to do something here for formats
!  that have all columns in the same file
!
 select case(trim(lcase(outformat)))
 case('gridascii','grid')
 !
 !--ascii output uses individual files
 !
    print "(/,a,i2)",'-----> WRITING TO ASCII OUTPUT FILES'

 case('gridbinary','gridbin')
 !
 !--simple unformatted binary format
 !
    filename = trim(filenamein)//'.grid'
    print "(/,a,i2)",'----> WRITING TO '//trim(filename)//' on unit ',iunit
    print "(a)",     '      (using unformatted binary format)'
    open(unit=iunit,file=trim(filename),form='unformatted',status='replace',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(filename)//' for output!'
       return
    endif

    write(iunit,iostat=ierr) npixels(1:3),ncolumns,time
    if (ierr /= 0) then
       print "(a)",' ERROR writing header to file!'
       return
    endif
 
 case('hdf5')

 case default
    ! return error if bad format
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in open_gridfile'
    ierr = 1
    return
 end select
 
end subroutine open_gridfile

!------------------------------------------------------
! write a particular column to the grid output file
!------------------------------------------------------
subroutine write_grid(iunit,filenamein,outformat,dat,npixels,label,time,ierr)
 use asciiutils, only:ucase,lcase,safename
 use filenames,  only:tagline
 implicit none
 integer, intent(in)                :: iunit
 character(len=*), intent(in)       :: filenamein,outformat
 real, dimension(:,:,:), intent(in) :: dat
 integer, dimension(3), intent(in)  :: npixels
 character(len=*), intent(in)       :: label
 real, intent(in)                   :: time
 integer, intent(out)               :: ierr
 character(len=len(filenamein)+20)  :: filename
 integer :: i,j,k
 
 select case(trim(lcase(outformat)))
 case('gridascii','grid')
    filename = trim(filenamein)//'_'//trim(safename(label))//'_grid.dat'
    print "(a)",'-----> WRITING '//trim(ucase(label))//' to '//trim(filename)
    
    !
    !--open ascii file
    !
    open(unit=iunit,file=trim(filename),form='formatted',status='replace',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'        
       return
    endif
    write(iunit,"(a)",iostat=ierr) '# '//trim(tagline)
    write(iunit,"(a)",iostat=ierr) &
      '# '//trim(filename)//' produced using "splash to '//trim(outformat)// &
      '" on file '//trim(filenamein)
    write(iunit,"(a)",iostat=ierr) '#'
    write(iunit,"(a)",iostat=ierr) '# time:'
    write(iunit,"(a,es15.7)",iostat=ierr) '# ',time
    write(iunit,"(a)",iostat=ierr) '#'
    write(iunit,"(a)",iostat=ierr) '# file contains:'
    write(iunit,"(a)",iostat=ierr) '# '//trim(label)//' interpolated to 3D grid '
    write(iunit,"(a)",iostat=ierr) '#'
    write(iunit,"(a)",iostat=ierr) '# written in the form: '
    write(iunit,"(a)",iostat=ierr) '#   do k=1,nz'
    write(iunit,"(a)",iostat=ierr) '#      do j=1,ny'
    write(iunit,"(a)",iostat=ierr) '#         write(*,*) (dat(i,j,k),i=1,nx)'
    write(iunit,"(a)",iostat=ierr) '#      enddo'
    write(iunit,"(a)",iostat=ierr) '#   enddo'
    write(iunit,"(a)",iostat=ierr) '#'
    write(iunit,"(a)",iostat=ierr) '# grid dimensions:'
    write(iunit,"(a)",iostat=ierr) '# nx    ny    nz'
    write(iunit,*) npixels(1:3)
    do k=1,npixels(3)
       do j=1,npixels(2)
          write(iunit,"(2048(es14.6,1x))") (dat(i,j,k),i=1,npixels(1))
       enddo
    enddo
    close(unit=iunit)

 case('gridbinary','gridbin')
    print "(a)",'-----> WRITING '//trim(ucase(label))
    write(iunit,iostat=ierr) (((dat(i,j,k),i=1,npixels(1)),j=1,npixels(2)),k=1,npixels(3))
 case('hdf5')

 case default
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in open_gridfile'
    return
 end select
 
end subroutine write_grid

end module write_griddata
