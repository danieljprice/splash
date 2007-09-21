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
 public :: isoutputformat,writepixmap
 
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
 case('ascii')
     isoutputformat = .true.
 end select
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
! case('ppm')
!    call write_pixmap_ppm(datpix,npixx,npixy,xmin,ymin,dx,datmin,datmax,label)
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
 write(*,"(a)",ADVANCE='NO'), '> writing pixel map to file '//trim(filename)//' ...'

 write(stringx,"(i10)") npixx
 write(stringy,"(i10)") npixy
 write(iunit,"(a)",err=66) '# '//trim(filename)//' created by splash (c) 2005-2007 Daniel Price'
 write(iunit,"(a)",err=66) '# Contains 2D pixel array '//trim(adjustl(stringx))//' x '//trim(adjustl(stringy))//' written as '
 write(iunit,"(a)",err=66) '#   do j=1,'//trim(adjustl(stringy))
 write(iunit,"(a)",err=66) '#      write(*,*) dat(1:'//trim(adjustl(stringx))//',j)'
 write(iunit,"(a)",err=66) '#   enddo'
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# '//trim(label)//': min = ',datmin,' max = ',datmax
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# x axis: min = ',xmin,' max = ',xmin+(npixx-1)*dx
 write(iunit,"(a,1pe14.6,a,1pe14.6)",err=66) '# y axis: min = ',ymin,' max = ',ymin+(npixy-1)*dx
 write(iunit,"(a)",err=66),'# '//trim(adjustl(stringx))//' '//trim(adjustl(stringy))
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

end module write_pixmap
