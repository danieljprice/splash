!-----------------------------------------------------------------
!     module containing output routines for writing pixel map
!     to output file in various formats
!
!     (c) D. Price 21/09/07
!     Added read routines June 2009
!-----------------------------------------------------------------
module write_pixmap
 use filenames, only:fileprefix
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
subroutine writepixmap(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,istep)
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in), dimension(npixx,npixy) :: datpix
 real, intent(in) :: xmin,ymin,dx,datmin,datmax
 character(len=*), intent(in) :: label
 integer, intent(in) :: istep
 
 select case(trim(pixmapformat))
 case('ascii')
    call write_pixmap_ascii(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,istep)
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
subroutine write_pixmap_ascii(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label,istep)
 implicit none
 integer, intent(in) :: npixx,npixy,istep
 real, intent(in), dimension(npixx,npixy) :: datpix
 real, intent(in) :: xmin,ymin,dx,datmin,datmax
 character(len=*), intent(in) :: label
 character(len=10) :: stringx,stringy
 character(len=30) :: fmtstring
 character(len=len_trim(fileprefix)+10) :: filename
 integer :: ierr,j
 integer, parameter :: iunit = 166
!
!--write ascii file
!
 write(filename,"(a,i5.5,a)") trim(fileprefix)//'_',istep,'.dat'
 open(unit=iunit,file=filename,status='replace',form='formatted',iostat=ierr)
 if (ierr /=0) then
    print*,'error opening '//trim(filename)
    return
 endif
 write(*,"(a)",ADVANCE='NO') '> writing pixel map to file '//trim(filename)//' ...'

 write(stringx,"(i10)") npixx
 write(stringy,"(i10)") npixy
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by splash (c) 2005-2009 Daniel Price'
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
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by splash (c) 2005-2009 Daniel Price'
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
subroutine readpixmap(datpix,npixx,npixy,dumpfile,label,icol,time,istep,xsec,ierr)
 use asciiutils, only:safename,basename
 implicit none
 real, intent(out), dimension(:,:), allocatable :: datpix
 integer, intent(out)            :: npixx,npixy,ierr
 character(len=*), intent(in)    :: dumpfile
 character(len=*), intent(inout) :: label
 integer,          intent(in)    :: icol
 real,             intent(in)    :: time
 integer,          intent(in)    :: istep
 logical,          intent(in)    :: xsec
 integer            :: i,maxnames
 integer, parameter :: iunit = 168
 character(len=128) :: filename
 character(len=len(dumpfile)) :: dumpfilei
 logical :: iexist,printinfo
 
 ierr = 0
 select case(trim(adjustl(readpixformat)))
 case('ascii') ! splash pixmap output files
    write(filename,"(a,i5.5,a)") trim(fileprefix)//'_',istep,'.dat'

    open(unit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
    if (ierr /=0) then
       print*,'error opening '//trim(filename)
       return
    else
       close(iunit)
    endif

 case('ftn','ftn512','chf') ! Christoph Federrath files
    npixx = 512
    npixy = 512
    !
    !--cycle through possible filenames
    !
    dumpfilei = dumpfile
    iexist = .false.
    i = 0
    filenametries=' '
    maxnames = 4
    printinfo = .false.
    do while (.not.iexist .and. i.lt.maxnames)
       i = i + 1
       select case(i)
       case(1)
          if (xsec) then
             filename = trim(dumpfilei)//'_'//trim(safename(label))//'_slice.pix'
          else
             filename = trim(dumpfilei)//'_'//trim(safename(label))//'_proj.pix'
          endif
       case(2)
          filename = trim(dumpfilei)//'_'//trim(safename(label))//'.pix'
       case(3)
          if (xsec) then
             filename = trim(dumpfilei)//'_slice.pix'
          else
             filename = trim(dumpfilei)//'_proj.pix'
          endif
       case(4)
          filename = trim(dumpfilei)//'.pix'
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
                print "(72('*'),/,'*',a,12x,'*')",' ERROR: could not find any .pix files with matching names:'
                dumpfilei = dumpfile
                printinfo = .true.
                i = 0
             endif
          endif
       endif
    enddo
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


end module write_pixmap
