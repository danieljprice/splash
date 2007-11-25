!-----------------------------------------------------------------
!     module containing output routines for writing pixel map
!     to output file in various formats
!
!     (c) D. Price 21/09/07
!-----------------------------------------------------------------
module write_pixmap
 use filenames, only:fileprefix
 implicit none
 logical, public :: iwritepixmap
 character(len=5), public :: pixmapformat
 public :: isoutputformat,writepixmap,write_pixmap_ppm
 
 private

contains

!-----------------------------------------------------------------
! utility to check if a format selection is valid
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
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by splash (c) 2005-2007 Daniel Price'
 write(iunit,"(a)",err=66) '# Contains 2D pixel array '//trim(adjustl(stringx))//' x '//trim(adjustl(stringy))//' written as '
 write(iunit,"(a)",err=66) '#   do j=1,'//trim(adjustl(stringy))
 write(iunit,"(a)",err=66) '#      write(*,*) dat(1:'//trim(adjustl(stringx))//',j)'
 write(iunit,"(a)",err=66) '#   enddo'
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# '//trim(label)//': min = ',datmin,' max = ',datmax
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# x axis: min = ',xmin,' max = ',xmin+(npixx-1)*dx
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# y axis: min = ',ymin,' max = ',ymin+(npixy-1)*dx
 write(iunit,"(a)",err=66) '# '//trim(adjustl(stringx))//' '//trim(adjustl(stringy))
 do j=1,npixy
    write(iunit,*,err=66) datpix(1:npixx,j)
 enddo
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
 write(iunit,"(a)",err=66) '# '//trim(adjustl(filename))//' created by splash (c) 2005-2007 Daniel Price'
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

end module write_pixmap
