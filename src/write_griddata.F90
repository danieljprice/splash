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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Module implementing "splash to grid" operation, writing
! 3D gridded data in various output formats
!-----------------------------------------------------------------
module readwrite_griddata
 implicit none

 public :: isgridformat,print_gridformats
 public :: open_gridfile_w,open_gridfile_r
 public :: write_grid,read_gridcolumn,write_gridlimits
 
 !
 !--generic interface for reading grid column data
 !  into 1D and 3D arrays
 !
 interface read_gridcolumn
  module procedure read_gridcolumn3D,read_gridcolumn1D
 end interface read_gridcolumn

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
 case('gridbinary','gridbin')
     isgridformat = .true.
 case('gridascii2')
     isgridformat = .true.
 end select

end function isgridformat

!-----------------------------------------------------------------
! print grid limits
!-----------------------------------------------------------------
subroutine write_gridlimits(ndim,xmin,xmax,labelx)
 integer, intent(in) :: ndim
 real,    intent(in) :: xmin(ndim),xmax(ndim)
 character(len=*), intent(in) :: labelx(ndim)
 integer :: i

 print "(a)",' grid dimensions:'
 do i=1,ndim
    if (maxval(abs(xmax)).lt.1.e7) then
       print "(1x,a,': ',f14.6,' -> ',f14.6)",trim(labelx(i)),xmin(i),xmax(i)
    else
       print "(1x,a,': ',es14.6,' -> ',es14.6)",trim(labelx(i)),xmin(i),xmax(i)
    endif
 enddo

end subroutine write_gridlimits

!-----------------------------------------------------------------
! print usage if format selection not valid
!-----------------------------------------------------------------
subroutine print_gridformats
 implicit none

 print "(/,a)",' Grid conversion mode ("splash to X dumpfiles"): '
 print "(a)",'    splash to grid         : interpolate basic SPH data (density, plus velocity if present in data)'
 print "(a)",'                             to 2D or 3D grid, write grid data to file (using default output=ascii)'
 print "(a)",'           to gridascii    : as above, grid data written in ascii format'
 print "(a)",'           to gridascii2   : grid data written in ascii format, all in one file'
 print "(a)",'           to gridbinary   : as above, grid data in simple unformatted binary format:'
 print "(a)",'                                write(unit) nx,ny,nz,ncolumns,time                 [ 4 bytes each ]'
 print "(a)",'                                write(unit) (((rho(i,j,k),i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]'
 print "(a)",'                                write(unit) (((vx(i,j,k), i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]'
 print "(a)",'                                write(unit) (((vy(i,j,k), i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]'
 print "(a)",'                                write(unit) (((...(i,j,k),i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]'
 print "(a)",'        allto grid         : as above, interpolating *all* columns to the grid (and output file)'
 print "(a)",'        allto gridascii    : as above, with ascii output'
 print "(a)",'        allto gridbinary   : as above, with binary output'

 return
end subroutine print_gridformats

!------------------------------------------------------
! open grid file for (write) output, write header
!------------------------------------------------------
subroutine open_gridfile_w(iunit,filenamein,outformat,ndim,ncolumns,npixels,time,ierr)
 use asciiutils, only:lcase
 implicit none
 integer, intent(in)               :: iunit
 character(len=*), intent(in)      :: filenamein,outformat
 character(len=len(filenamein)+10) :: filename
 integer, intent(in)                  :: ndim,ncolumns
 integer, dimension(ndim), intent(in) :: npixels
 real, intent(in)                     :: time
 integer, intent(out)                 :: ierr
!  
!--Only have to do something here for formats
!  that have all columns in the same file
!
 ierr = 0
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

    write(iunit,iostat=ierr) npixels(1:ndim),ncolumns,time
    if (ierr /= 0) then
       print "(a)",' ERROR writing header to file!'
       return
    endif

 case('gridascii2')

    print "(/,a,i2)",'-----> WRITING TO ASCII OUTPUT FILES (WITH X, Y, Z, COL)'
     
 case('hdf5')

 case default
    ! return error if bad format
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in open_gridfile'
    ierr = 1
    return
 end select
 
end subroutine open_gridfile_w


!------------------------------------------------------
! open grid file for reading, read header
!------------------------------------------------------
subroutine open_gridfile_r(iunit,filename,informat,ndim,ncolumns,npixels,time,ierr)
 use asciiutils, only:lcase
 implicit none
 integer, intent(in)                :: iunit,ndim
 character(len=*), intent(in)       :: filename
 character(len=*), intent(inout)    :: informat
 integer, intent(out)               :: ncolumns
 integer, dimension(ndim), intent(out) :: npixels
 real, intent(out)                     :: time
 integer, intent(out)                  :: ierr
!
!--read only implemented for binary grid format at present
!
 ierr = 0
 select case(trim(lcase(informat)))
 case('gridbinary','gridbin')
 !
 !--simple unformatted binary format
 !
    print "(/,a,i2)",'----> READING '//trim(filename)//' on unit ',iunit
    print "(a)",     '      (using unformatted binary format)'
    open(unit=iunit,file=trim(filename),form='unformatted',status='old',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(filename)//' for reading!'
       return
    endif

    read(iunit,iostat=ierr) npixels(1:ndim),ncolumns,time
    if (ierr /= 0) then
       print "(a)",' ERROR reading header to file!'
       return
    endif
 case default
    ! return error if bad format
    print "(a)",' ERROR: cannot read grid format '''//trim(informat)//''' in open_gridfile_r'
    ierr = 1
    return
 end select
 
end subroutine open_gridfile_r

!------------------------------------------------------
! write a particular column to the grid output file
!------------------------------------------------------
subroutine write_grid(iunit,filenamein,outformat,ndim,ncolgrid,npixels,label,time,&
                      pixwidth,xmin,ierr,dat,dat3D,dat2D,label3D)
 use asciiutils, only:ucase,lcase,safename
 use filenames,  only:tagline
 implicit none
 integer, intent(in)                :: iunit
 character(len=*), intent(in)       :: filenamein,outformat
 integer, intent(in)                  :: ndim,ncolgrid
 integer, dimension(ndim), intent(in) :: npixels
 character(len=*), intent(in)         :: label
 real, intent(in)                     :: time,pixwidth
 real, dimension(3), intent(in)       :: xmin
 integer, intent(out)                 :: ierr
 character(len=len(filenamein)+20)    :: filename
 real, dimension(:,:,:),   intent(in), optional :: dat
 real, dimension(:,:,:,:), intent(in), optional :: dat3D
 real, dimension(:,:),     intent(in), optional :: dat2D
 character(len=*), intent(in), optional :: label3D(ncolgrid)
 integer :: i,j,k,n
 real    :: xi,yi,zi
 
 ierr = 0
 if (ndim.eq.3 .and. .not.(present(dat3D) .or. present(dat))) then
    print "(a)",' ERROR in call to write_grid: ndim=3 but 3D grid not passed'
    ierr = 1
 elseif (ndim.eq.2 .and. .not.present(dat2D)) then
    print "(a)",' ERROR in call to write_grid: ndim=2 but 2D grid not passed'
    ierr = 1
 elseif (.not.(ndim.eq.2 .or. ndim.eq.3)) then
    print "(a,i2,a)",' ERROR in call to write_grid: cannot write grid for ',ndim,' dimensions'
    ierr = 2
 endif
 if (ierr /= 0) return
 
 select case(trim(lcase(outformat)))
 case('gridascii','grid')
    if (ncolgrid > 1) then
       filename = trim(filenamein)//'_grid.dat'
    else
       filename = trim(filenamein)//'_'//trim(safename(label))//'_grid.dat'
    endif
    print "(/,a)",'-----> WRITING to '//trim(filename)
    
    !
    !--open ascii file
    !
    open(unit=iunit,file=trim(filename),form='formatted',status='replace',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'        
       return
    endif
    write(iunit,"(a)",err=100) '# '//trim(tagline)
    write(iunit,"(a)",err=100) &
      '# '//trim(filename)//' produced using "splash to '//trim(outformat)// &
      '" on file '//trim(filenamein)
    write(iunit,"(a)",err=100) '#'
    write(iunit,"(a)",err=100) '# time:'
    write(iunit,"(a,es15.7)",iostat=ierr) '# ',time
    write(iunit,"(a)",err=100) '#'
    write(iunit,"(a)",err=100) '# file contains:'
    write(iunit,"(a,i1,a)",err=100) '# '//trim(label)//' interpolated to ',ndim,'D grid '
    write(iunit,"(a)",err=100) '#'
    write(iunit,"(a)",err=100) '# written in the form: '
    if (ndim.eq.3) then
       write(iunit,"(a)",err=100) '#   do k=1,nz'
       write(iunit,"(a)",err=100) '#      do j=1,ny'
       write(iunit,"(a)",err=100) '#         write(*,*) (dat(i,j,k),i=1,nx)'
       write(iunit,"(a)",err=100) '#      enddo'
       write(iunit,"(a)",err=100) '#   enddo'
    else
       write(iunit,"(a)",err=100) '#   do j=1,ny'
       write(iunit,"(a)",err=100) '#      write(*,*) (dat(i,j),i=1,nx)'
       write(iunit,"(a)",err=100) '#   enddo'
    endif
    write(iunit,"(a)",err=100) '#'
    write(iunit,"(a)",err=100) '# grid dimensions:'
    if (present(dat)) then
       write(iunit,"(a)",err=100) '# nx    ny    nz'
       write(iunit,*,err=100) npixels(1:ndim)
       do k=1,npixels(3)
          do j=1,npixels(2)
             write(iunit,"(2048(es14.6,1x))",err=100) (dat(i,j,k),i=1,npixels(1))
          enddo
       enddo
    elseif (present(dat3D)) then
       write(iunit,"(a)",err=100) '# nx    ny    nz'
       write(iunit,*,err=100) npixels(1:ndim)
       do k=1,npixels(3)
          do j=1,npixels(2)
             write(iunit,"(2048(es14.6,1x))",err=100) (dat3D(1,i,j,k),i=1,npixels(1))
          enddo
       enddo
    elseif (present(dat2D)) then
       write(iunit,"(a)",err=100) '# nx    ny'
       write(iunit,*,err=100) npixels(1:ndim)
       do j=1,npixels(2)
          write(iunit,"(2048(es14.6,1x))",err=100) (dat2D(i,j),i=1,npixels(1))
       enddo
    endif
    close(unit=iunit)
    return

 case('gridbinary','gridbin')
    print "(a)",'-----> WRITING '//trim(ucase(label))
    if (present(dat)) then
       write(iunit,iostat=ierr) (((dat(i,j,k),i=1,npixels(1)),j=1,npixels(2)),k=1,npixels(3))
    elseif (present(dat3D)) then
       write(iunit,iostat=ierr) ((((dat3D(n,i,j,k),i=1,npixels(1)),j=1,npixels(2)),k=1,npixels(3)),n=1,ncolgrid)
    elseif (present(dat2D)) then
       write(iunit,iostat=ierr) ((dat2D(i,j),i=1,npixels(1)),j=1,npixels(2))    
    endif
 case('gridascii2')
    if (ncolgrid > 1) then
       filename = trim(filenamein)//'_grid.dat'
       print "(a)",'-----> WRITING to '//trim(filename)
    else
       filename = trim(filenamein)//'_'//trim(safename(label))//'_grid.dat'
       print "(a)",'-----> WRITING '//trim(ucase(label))//' to '//trim(filename)
    endif
    
    !
    !--open ascii file
    !
    open(unit=iunit,file=trim(filename),form='formatted',status='replace',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'        
       return
    endif
    write(iunit,"(a)",err=100) '# '//trim(tagline)
    write(iunit,"(a)",err=100) &
      '# '//trim(filename)//' produced using "splash to '//trim(outformat)// &
      '" on file '//trim(filenamein)
    write(iunit,"(a)",err=100) '# grid dimensions:'
    if (present(dat3D) .and. present(dat)) then
       write(iunit,"(a)",err=100) '# nx    ny    nz'
       write(iunit,"(a,3(i5,1x))",err=100) '# ',npixels(1:3)
       write(iunit,"('#',64('[',a13,']'))",err=100) 'x','y','z',trim(label),(trim(label3D(n)),n=1,ncolgrid)
       do k=1,npixels(3)
          write(*,"('.')",ADVANCE='NO')
          zi = xmin(3) + (k-0.5)*pixwidth
          do j=1,npixels(2)
             yi = xmin(2) + (j-0.5)*pixwidth
             do i=1,npixels(1)
                xi = xmin(1) + (i-0.5)*pixwidth
                write(iunit,"(64(es14.6,1x))") xi,yi,zi,dat(i,j,k),(dat3D(n,i,j,k),n=1,ncolgrid)
             enddo
          enddo
       enddo
    elseif (present(dat3D)) then
       write(iunit,"(a)",err=100) '# nx    ny    nz'
       write(iunit,"(a,3(i5,1x))",err=100) '# ',npixels(1:3)
       write(iunit,"('#',64('[',a13,']'))",err=100) 'x','y','z',(trim(label3D(n)),n=1,ncolgrid)
       do k=1,npixels(3)
          write(*,"('.')",ADVANCE='NO')
          zi = xmin(3) + (k-0.5)*pixwidth
          do j=1,npixels(2)
             yi = xmin(2) + (j-0.5)*pixwidth
             do i=1,npixels(1)
                xi = xmin(1) + (i-0.5)*pixwidth
                write(iunit,"(64(es14.6,1x))") xi,yi,zi,(dat3D(n,i,j,k),n=1,ncolgrid)
             enddo
          enddo
       enddo
    elseif (present(dat2D)) then
       write(iunit,"(a)",err=100) '# nx    ny '
       write(iunit,"(a,2(i5,1x))",err=100) '# ',npixels(1:2)
       write(iunit,"('#',3('[',a13,']'))",err=100) 'x','y',trim(label)
       do j=1,npixels(2)
          write(*,"('.')",ADVANCE='NO')
          yi = xmin(2) + (j-0.5)*pixwidth
          do i=1,npixels(1)
             xi = xmin(1) + (i-0.5)*pixwidth
             write(iunit,"(3(es14.6,1x))") xi,yi,dat2D(i,j)
          enddo
       enddo
    endif
    write(*,*)

 case('hdf5')

 case default
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in write_grid'
    return
 end select

 return
!
!--error handling during write
!
100 continue
    print "(a)",' ERROR writing grid file'
    close(unit=iunit)
    return

end subroutine write_grid

!------------------------------------------------------------------
! read a particular column from the grid output file into 3D array
!------------------------------------------------------------------
subroutine read_gridcolumn3D(iunit,dat,npixels,ierr)
 implicit none
 integer, intent(in)                 :: iunit
 real, dimension(:,:,:), intent(out) :: dat
 integer, dimension(3), intent(in)   :: npixels
 integer, intent(out)                :: ierr
 integer :: i,j,k
 
 print "(a,i4,'x',i4,'x',i4,a)",'-----> READING ',npixels(:),' data points'
 read(iunit,iostat=ierr) (((dat(i,j,k),i=1,npixels(1)),j=1,npixels(2)),k=1,npixels(3))

end subroutine read_gridcolumn3D

!------------------------------------------------------------------
! read a particular column from the grid output file into 1D array
!------------------------------------------------------------------
subroutine read_gridcolumn1D(iunit,dat,ngrid,ierr)
 implicit none
 integer, intent(in)                 :: iunit
 real, dimension(:), intent(out)     :: dat
 integer, intent(in)                 :: ngrid
 integer, intent(out)                :: ierr
 integer :: i
 
 print "(a,i10,a)",'-----> READING ',ngrid,' data points'
 read(iunit,iostat=ierr) (dat(i),i=1,ngrid)

end subroutine read_gridcolumn1D

end module readwrite_griddata
