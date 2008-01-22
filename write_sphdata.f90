!-----------------------------------------------------------------
!     module containing output routines for writing SPH data
!     as read by SPLASH to output file in various formats
!
!     (c) D. Price 22/01/08
!-----------------------------------------------------------------
module write_sphdata
 public :: issphformat,write_sphdump
 private

contains

!-----------------------------------------------------------------
! utility to check if a format selection is valid
!-----------------------------------------------------------------
logical function issphformat(string)
 implicit none
 character(len=*), intent(in) :: string

 issphformat = .false.
 select case(trim(string))
 case('ascii')
     issphformat = .true.
 end select
 
 if (.not.issphformat) then
    print "(a)",' possible formats for convert mode ("splash to X"): '
    print "(a)",' splash to ascii : convert SPH data to ascii file'
 endif
 
 return
end function issphformat

subroutine write_sphdump(dat,npart,ncolumns,filename,outformat)
 implicit none
 integer, intent(in) :: npart,ncolumns
 real, intent(in), dimension(npart,ncolumns) :: dat
 character(len=*), intent(in) :: filename,outformat
 integer, parameter :: iunit = 83
 integer :: ierr,i
 character(len=40) :: fmtstring
 
 write(fmtstring,"(i10,a)") ncolumns,'(1pe15.7,1x)'
 fmtstring = '('//trim(adjustl(fmtstring))//')'

 select case(trim(outformat))
 case ('ascii')
    print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//'.ascii WITH ',ncolumns,' COLUMNS'
    open(unit=iunit,file=trim(filename)//'.ascii',status='replace',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR OPENING FILE FOR WRITING'
          return
       endif

       do i=1,npart
          write(iunit,fmtstring,err=100) dat(i,1:ncolumns)
       enddo
    close(iunit)

    return
100 continue
    close(iunit)
    print*,'ERROR WRITING ASCII FILE'
    return
 case default
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in write_sphdump'
    return
 end select
 
end subroutine write_sphdump

end module write_sphdata
