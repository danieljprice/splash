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

!-----------------------------------------------------------------
!     module containing output routines for writing pixel map
!     to output file in various formats
!
!     (c) D. Price 21/09/07
!     Added read routines June 2009
!-----------------------------------------------------------------
module write_pixmap
 use filenames, only:fileprefix,tagline
 implicit none
 logical, public :: iwritepixmap = .false.
 logical, public :: ireadpixmap = .false.
 character(len=5), public :: pixmapformat = ' '
 character(len=5), public :: readpixformat = ' '
 public :: isoutputformat,writepixmap,write_pixmap_ppm
 public :: isinputformat,readpixmap

 private

contains

!-----------------------------------------------------------------
! utility to check if an output format selection is valid
!-----------------------------------------------------------------
logical function isoutputformat(string)
 implicit none
 character(len=*), intent(in) :: string

 isoutputformat = .false.
 select case(trim(string))
 case('ascii','ppm')
     isoutputformat = .true.
 end select

 if (.not.isoutputformat) then
    print "(a)",' possible formats for -o option: '
    print "(a)",' -o ppm   : dump pixel map to ppm file'
    print "(a)",' -o ascii : dump pixel map to ascii file'
    print "(a)",' use -p to change the prefix for the filenames'
 endif

 return
end function isoutputformat

!-----------------------------------------------------------------
! utility to check if an input format selection is valid
!-----------------------------------------------------------------
logical function isinputformat(string)
 implicit none
 character(len=*), intent(in) :: string

 isinputformat = .false.
 select case(trim(string))
 case('ascii','ftn','ftn512','chf')
     isinputformat = .true.
 end select

 if (.not.isinputformat) then
    print "(a)",' possible formats for -readpix option: '
    print "(a)",' -readpix ascii  : read pixel maps from ascii file'
    print "(a)",' -readpix ftn    : read pixel maps from unformatted fortran file "read(1) dat(:,:)"'
    print "(a)",' -readpix ftn512 : read pixel maps from unformatted fortran file "read(1) dat(1:512,1:512)"'
 endif

 return
end function isinputformat

!-----------------------------------------------------------------
!  wrapper routine for all output formats
!-----------------------------------------------------------------
subroutine writepixmap(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,labu,istep,xsec,dumpfile)
 implicit none
 integer, intent(in) :: npixx,npixy
 real,    intent(in), dimension(npixx,npixy) :: datpix
 real,    intent(in) :: xmin,ymin,dx,datmin,datmax
 logical, intent(in) :: xsec
 character(len=*), intent(in) :: label,labu,dumpfile
 integer, intent(in) :: istep

 select case(trim(pixmapformat))
 case('ascii')
    call write_pixmap_ascii(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,labu,istep,xsec,dumpfile)
 case('ppm')
    call write_pixmap_ppm(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,istep)
! case('fits')
!    call write_pixmap_fits(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label)
 case default
    print "(a)",' ERROR: invalid output format for pixel map '
 end select

end subroutine writepixmap

!-----------------------------------------------------------------
!   output pixmap as an ascii file
!-----------------------------------------------------------------
subroutine write_pixmap_ascii(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,labu,istep,xsec,dumpfile)
 use labels, only:shortlabel
 implicit none
 integer, intent(in) :: npixx,npixy,istep
 real,    intent(in), dimension(npixx,npixy) :: datpix
 real,    intent(in) :: xmin,ymin,dx,datmin,datmax
 logical, intent(in) :: xsec
 character(len=*), intent(in) :: label,labu,dumpfile
 character(len=10) :: stringx,stringy
 character(len=30) :: fmtstring
 character(len=len(dumpfile)+10) :: filename
 integer :: ierr,j
 integer, parameter :: iunit = 166
!
!--write ascii file
!
 !write(filename,"(a,i5.5,a)") trim(fileprefix)//'_',istep,'.dat'
 call get_pixmap_filename(filename,dumpfile,shortlabel(label,labu),'.pix',xsec)
 open(unit=iunit,file=filename,status='replace',form='formatted',iostat=ierr)
 if (ierr /=0) then
    print*,'error opening '//trim(filename)
    return
 endif
 write(*,"(a)",ADVANCE='NO') '> writing pixel map to file '//trim(filename)//' ...'

 write(stringx,"(i10)") npixx
 write(stringy,"(i10)") npixy
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by '//trim(tagline)
 write(iunit,"(a)",err=66) '# Contains 2D pixel array '//trim(adjustl(stringx))//' x '//trim(adjustl(stringy))//' written as '
 write(iunit,"(a)",err=66) '#   do j=1,'//trim(adjustl(stringy))
 write(iunit,"(a)",err=66) '#      write(*,*) dat(1:'//trim(adjustl(stringx))//',j)'
 write(iunit,"(a)",err=66) '#   enddo'
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# '//trim(label)//': min = ',datmin,' max = ',datmax
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# x axis: min = ',xmin,' max = ',xmin+(npixx-1)*dx
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# y axis: min = ',ymin,' max = ',ymin+(npixy-1)*dx
 write(iunit,"(a)",err=66) '# '//trim(adjustl(stringx))//' '//trim(adjustl(stringy))

 write(fmtstring,"(a,i6,a)",iostat=ierr) '(',npixx,'(1pe14.6))'
 if (ierr /= 0) then
    do j=1,npixy
       write(iunit,*,err=66) datpix(1:npixx,j)
    enddo
 else
    do j=1,npixy
       write(iunit,fmtstring,err=66) datpix(1:npixx,j)
    enddo
 endif
 close(iunit)
 print "(a)",'OK'
 return

66 continue
 print "(a)",' ERROR during write '
 close(iunit)
 return

end subroutine write_pixmap_ascii

!-----------------------------------------------------------------
!   output pixmap as a raw .ppm file
!-----------------------------------------------------------------
subroutine write_pixmap_ppm(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,istep,brightness)
 use colours, only:rgbtable,ncolours
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in), dimension(npixx,npixy) :: datpix
 real, intent(in), dimension(npixx,npixy), optional :: brightness
 real, intent(in) :: xmin,ymin,dx,datmin,datmax
 character(len=*), intent(in) :: label
 integer, intent(in) :: istep
 character(len=120) :: filename
 real, dimension(3) :: rgbi,drgb
 real :: dati,ddatrange,datfraci,ftable
 integer :: ipix,jpix,ir,ib,ig,ierr,maxcolour,indexi
 integer, parameter :: iunit = 167
!
!--check for errors
!
 if (abs(datmax-datmin).gt.tiny(datmin)) then
    ddatrange = 1./abs(datmax-datmin)
 else
    print "(a)",'error: datmin=datmax : pointless writing ppm file'
    return
 endif
!
!--write PPM--
!
 write(filename,"(a,i5.5,a)") trim(fileprefix)//'_',istep,'.ppm'
 open(unit=iunit,file=filename,status='replace',form='formatted',iostat=ierr)
 if (ierr /=0) then
    print*,'error opening ppm file'
    return
 endif
 write(*,"(a)",ADVANCE='NO') '> writing pixel map to file '//trim(filename)//' ...'
!
!--PPM header
!
 maxcolour = 255
 write(iunit,"(a)",err=66) 'P3'
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by '//trim(tagline)
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# '//trim(label)//': min = ',datmin,' max = ',datmax
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# x axis: min = ',xmin,' max = ',xmin+(npixx-1)*dx
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# y axis: min = ',ymin,' max = ',ymin+(npixy-1)*dx
 write(iunit,"(i4,1x,i4)",err=66) npixx, npixy
 write(iunit,"(i3)",err=66) maxcolour
!--pixel information
 do jpix = npixy,1,-1
    do ipix = 1,npixx
       dati = datpix(ipix,jpix)
       datfraci = (dati - datmin)*ddatrange
       datfraci = max(datfraci,0.)
       datfraci = min(datfraci,1.)
       !--define colour for current particle
       ftable = datfraci*ncolours
       indexi = int(ftable) + 1
       indexi = min(indexi,ncolours)
       if (indexi.lt.ncolours) then
       !--do linear interpolation from colour table
          drgb(:) = rgbtable(:,indexi+1) - rgbtable(:,indexi)
          rgbi(:) = rgbtable(:,indexi) + (ftable - int(ftable))*drgb(:)
       else
          rgbi(:) = rgbtable(:,indexi)
       endif
       if (present(brightness)) then
          rgbi(:) = rgbi(:)*min(brightness(ipix,jpix),1.0)
       endif
       ir = max(min(int(rgbi(1)*maxcolour),maxcolour),0)
       ig = max(min(int(rgbi(2)*maxcolour),maxcolour),0)
       ib = max(min(int(rgbi(3)*maxcolour),maxcolour),0)
       write(iunit,"(i3,1x,i3,1x,i3,2x)",err=66) ir,ig,ib
    enddo
 enddo
 close(unit=iunit)
 print "(a)",'OK'
 return

66 continue
 print "(a)",' ERROR during write '
 close(iunit)
 return
end subroutine write_pixmap_ppm


!-----------------------------------------------------------------
!  read in pixels from file
!-----------------------------------------------------------------
subroutine readpixmap(datpix,npixx,npixy,dumpfile,label,istep,xsec,ierr)
 use asciiutils, only:nheaderlines
 implicit none
 real, intent(out), dimension(:,:), allocatable :: datpix
 integer, intent(out)            :: npixx,npixy,ierr
 character(len=*), intent(in)    :: dumpfile
 character(len=*), intent(in)    :: label
 integer,          intent(in)    :: istep
 logical,          intent(in)    :: xsec
 integer            :: i,j,nheader,nerr
 integer, parameter :: iunit = 168
 character(len=128) :: filename
 character(len=2)   :: char
 logical :: iexist

 ierr = 0
 select case(trim(adjustl(readpixformat)))
 case('ascii') ! splash pixmap output files
    call check_for_pixmap_files(filename,dumpfile,label,'.pix',istep,xsec,iexist)
    if (.not.iexist) then
       print "('*',a,'*',/,72('*'))",' Create a file with one of these names (or a soft link) and try again '
       ierr = 1
       return
    endif

    open(unit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
    if (ierr /=0) then
       print*,'error opening '//trim(filename)
       return
    else
       npixx = 0
       npixy = 0
       nheader = nheaderlines(iunit)
       rewind(iunit)
       do i=1,nheader-1
          read(iunit,*,iostat=ierr)
       enddo
       read(iunit,*,iostat=ierr) char,npixx,npixy
       if (ierr /= 0 .or. npixx.le.0 .or. npixy.le.0) then
          print*,'ERROR reading size of pixel map, got nx = ',npixx,' ny = ',npixy,&
                 ', skipped ',nheader,' header lines'
       else
          print "(a,i5,a,i5,a)",' reading',npixx,' x',npixy,' pixel map from '//trim(filename)
       endif
       allocate(datpix(npixx,npixy),stat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR allocating memory for pixel map'
          close(iunit)
          return
       endif
       nerr = 0
       do j=1,npixy
          read(iunit,*,iostat=ierr) datpix(1:npixx,j)
          if (ierr /= 0) nerr = nerr + 1
       enddo
       if (nerr /= 0) print "(a,i3,a,i3)",' WARNING: ',nerr,' errors reading pixel map from '//trim(filename)//' on unit ',iunit

       close(iunit)
    endif

 case('ftn','ftn512','chf') ! Christoph Federrath files
    npixx = 512
    npixy = 512
    !
    !--cycle through possible filenames
    !
    call check_for_pixmap_files(filename,dumpfile,label,'.pix',istep,xsec,iexist)
    if (.not.iexist) then
       print "('*',a,'*',/,72('*'))",' Create a file with one of these names (or a soft link) and try again '
       ierr = 1
       return
    endif

    open(unit=iunit,file=filename,status='old',form='unformatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR: cannot open '//trim(filename)
       ierr = 2
       return
    else
       print "(a)",' reading pixel map from '//trim(filename)
       allocate(datpix(npixx,npixy),stat=ierr)
       read(iunit,iostat=ierr) datpix
       if (ierr /= 0) print "(a,i3)",' WARNING: ERRORS reading pixel map from '//trim(filename)//' on unit ',iunit
       close(iunit)
    endif
 case default
    if (len_trim(readpixformat).le.0) then
       print "(a)",' ERROR: pixel format not set prior to read_pixmap call'
    else
       print "(a)",' ERROR: unknown pixmap format '//trim(adjustl(readpixformat))
    endif
    ierr = 1
 end select

end subroutine readpixmap

!-----------------------------------------------------------------
!  look for pixel map files matching a variety of naming schemes
!-----------------------------------------------------------------
subroutine get_pixmap_filename(filename,dumpfile,label,ext,xsec)
 use asciiutils, only:basename,safename
 character(len=*), intent(out) :: filename
 character(len=*), intent(in)  :: dumpfile,label,ext
 logical,          intent(in)  :: xsec

 if (xsec) then
    filename = trim(basename(dumpfile))//'_'//trim(safename(label))//'_slice'//trim(ext)
 else
    filename = trim(basename(dumpfile))//'_'//trim(safename(label))//'_proj'//trim(ext)
 endif

end subroutine get_pixmap_filename

!-----------------------------------------------------------------
!  look for pixel map files matching a variety of naming schemes
!-----------------------------------------------------------------
subroutine check_for_pixmap_files(filename,dumpfile,label,ext,istep,xsec,iexist)
 use asciiutils, only:safename,basename
 implicit none
 character(len=*), intent(out) :: filename
 character(len=*), intent(in)  :: dumpfile,label,ext
 integer,          intent(in)  :: istep
 logical,          intent(in)  :: xsec
 logical,          intent(out) :: iexist
 character(len=len(dumpfile))  :: dumpfilei
 integer :: i,maxnames
 logical :: printinfo

 dumpfilei = dumpfile
 iexist = .false.
 i = 0
 maxnames = 5
 printinfo = .false.

 do while (.not.iexist .and. i.lt.maxnames)
    i = i + 1
    select case(i)
    case(1)
       if (xsec) then
          filename = trim(dumpfilei)//'_'//trim(safename(label))//'_slice'//trim(ext)
       else
          filename = trim(dumpfilei)//'_'//trim(safename(label))//'_proj'//trim(ext)
       endif
    case(2)
       filename = trim(dumpfilei)//'_'//trim(safename(label))//trim(ext)
    case(3)
       if (xsec) then
          filename = trim(dumpfilei)//'_slice'//trim(ext)
       else
          filename = trim(dumpfilei)//'_proj'//trim(ext)
       endif
    case(4)
       filename = trim(dumpfilei)//trim(ext)
    case(5)
       write(filename,"(a,i5.5,a)") trim(fileprefix)//'_',istep,trim(ext)
    end select
    !
    !--query to see if file exists
    !
    if (printinfo) then
       print "('*',a,'*')",'  no file '//filename(1:60)//''
    else
       inquire(file=filename,exist=iexist)
    endif
    if (.not.iexist) then
       !
       !--try the same files again but in the current directory
       !  instead of the directory in which the dump files are located
       !
       if (i.eq.maxnames) then
          if (len_trim(dumpfilei).ne.len_trim(basename(dumpfile))) then
             i = 0
             dumpfilei = basename(dumpfile)
          elseif (.not.printinfo) then
             print "(72('*'),/,'*',a,12x,'*')", &
                   ' ERROR: could not find any '//trim(ext)//' files with matching names:'
             dumpfilei = dumpfile
             printinfo = .true.
             i = 0
          endif
       endif
    endif
 enddo
 
end subroutine check_for_pixmap_files

end module write_pixmap
