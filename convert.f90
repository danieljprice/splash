!-----------------------------------------------------------------
!     interface routine to convert all of the dump files
!     on the SPLASH command line to a specified output format
!
!     (c) D. Price 22/01/08
!-----------------------------------------------------------------
module convert
 implicit none

contains

subroutine convert_all(outformat,igotfilenames)
 use particle_data, only:time,dat,npartoftype,iamtype
 use settings_data, only:ncolumns,ncalc,required,ntypes
 use filenames, only:rootname,nstepsinfile,nfiles
 use write_sphdata, only:write_sphdump
 use getdata, only:get_data
 use prompting, only:ucase
 implicit none
 character(len=*), intent(in) :: outformat
 logical, intent(in) :: igotfilenames
 integer :: ifile,idump,ntotal
 character(len=len(rootname)+4) :: filename
 
 required = .true.
 print "(/,5('-'),a,/)",'> CONVERTING DUMPFILES TO '//trim(ucase(outformat))//' FORMAT '
 do ifile=1,nfiles
    !--read data from dump file + calculate extra columns
    if (ifile.eq.1) then
       call get_data(ifile,igotfilenames)
    else
       call get_data(ifile,.true.)
    endif
    !--dump each step in file to an output file
    do idump = 1,nstepsinfile(ifile)
       if (idump.gt.1) then
          write(filename,"(a,'_',i3.3)") rootname(ifile),idump
       endif
       ntotal = sum(npartoftype(1:ntypes,idump))
       call write_sphdump(time(idump),dat(1:ntotal,1:ncolumns+ncalc,idump),ntotal,ntypes, &
                          npartoftype(1:ntypes,idump),iamtype(:,idump), &
                          ncolumns+ncalc,rootname(ifile),outformat)
    enddo
 enddo
 print "(/,5('-'),a,/)",'> FINISHED CONVERTING DUMP FILES '

 return
end subroutine convert_all

end module convert
