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
 use particle_data, only:time,dat,npartoftype,masstype,iamtype
 use settings_data, only:ncolumns,ncalc,required,ntypes,ndimV
 use filenames, only:rootname,nstepsinfile,nfiles
 use write_sphdata, only:write_sphdump
 use analysis, only:isanalysis,open_analysis,write_analysis,close_analysis
 use getdata, only:get_data
 use asciiutils, only:ucase
 implicit none
 character(len=*), intent(in) :: outformat
 logical, intent(inout) :: igotfilenames
 logical :: doanalysis
 integer :: ifile,idump,ntotal
 character(len=len(rootname)+4) :: filename
 character(len=10) :: string
 
 required = .true.
 doanalysis = isanalysis(outformat,noprint=.true.)
 
 if (.not.doanalysis) then
 !
 !--for format conversion each dump file is independent
 !
    print "(/,5('-'),a,/)",'> CONVERTING DUMPFILES TO '//trim(ucase(outformat))//' FORMAT '
 endif
 
 !
 !--if nfiles = 0 (ie. no files read from command line), then call get_data here
 !  to also get nfiles correctly prior to the loop
 !
 if (nfiles.eq.0) then
    call get_data(1,igotfilenames)
    igotfilenames = .true.
 endif
 
 do ifile=1,nfiles
    !--read data from dump file + calculate extra columns
    if (ifile.eq.1) then
       call get_data(ifile,igotfilenames)
 !
 !--for analysis we need to initialise the output file
 !  and close it at the end - do this here so we know
 !  the first filename and ndimV, labels etc.
 !
       if (doanalysis) then
          call open_analysis(rootname(1),outformat,required,ncolumns+ncalc,ndimV)
       endif
    else
       call get_data(ifile,.true.)
    endif
    !--dump each step in file to an output file
    do idump = 1,nstepsinfile(ifile)
       if (idump.gt.1) then
          write(filename,"(a,'_',i3.3)") rootname(ifile),idump
       endif
       ntotal = sum(npartoftype(1:ntypes,idump))
       if (doanalysis) then
          call write_analysis(time(idump),dat(1:ntotal,1:ncolumns+ncalc,idump),ntotal,ntypes, &
                          npartoftype(1:ntypes,idump),masstype(1:ntypes,idump),iamtype(:,idump), &
                          ncolumns+ncalc,ndimV,outformat)
       else
          call write_sphdump(time(idump),dat(1:ntotal,1:ncolumns+ncalc,idump),ntotal,ntypes, &
                          npartoftype(1:ntypes,idump),iamtype(:,idump), &
                          ncolumns+ncalc,rootname(ifile),outformat)
       endif
    enddo
 enddo

 !--for analysis we need to start and end differently
 if (doanalysis) then
    call close_analysis
    write(string,"(i10)") nfiles
    print "(/,5('-'),a,i5,a,/)",'> FINISHED CALCULATING '//trim(ucase(outformat))//' (USED '//trim(adjustl(string))//' DUMP FILES)'
 else
    print "(/,5('-'),a,/)",'> FINISHED CONVERTING DUMP FILES '
 endif
 return
end subroutine convert_all

end module convert
