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
!     interface routine to convert all of the dump files
!     on the SPLASH command line to a specified output format
!
!     (c) D. Price 22/01/08
!-----------------------------------------------------------------
module convert
 implicit none

contains

subroutine convert_all(outformat,igotfilenames,useall)
 use particle_data, only:time,gamma,dat,npartoftype,masstype,iamtype
 use settings_data, only:ncolumns,ncalc,required,ntypes,ndim,ndimV,lowmemorymode
 use filenames,     only:rootname,nstepsinfile,nfiles,limitsfile
 use write_sphdata, only:write_sphdump
 use readwrite_griddata, only:isgridformat
 use analysis,      only:isanalysis,open_analysis,write_analysis,close_analysis
 use convert_grid,  only:convert_to_grid
 use getdata,       only:get_data
 use asciiutils,    only:ucase
 use limits,        only:read_limits
 implicit none
 character(len=*), intent(in) :: outformat
 logical, intent(inout)       :: igotfilenames
 logical, intent(in)          :: useall
 logical :: doanalysis,converttogrid
 integer :: ifile,idump,ntotal,ierr
 character(len=len(rootname)+4) :: filename
 character(len=10) :: string

 required      = .true.  ! read whole dump file by default
 doanalysis    = isanalysis(outformat,noprint=.true.)
 converttogrid = isgridformat(outformat)
 lowmemorymode = .false. ! must not be true for first file

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
       call get_data(ifile,igotfilenames,firsttime=.true.)
       !
       ! read plot limits from file (overrides get_data limits settings)
       !
       call read_limits(trim(limitsfile),ierr)
 !
 !--for analysis we need to initialise the output file
 !  and close it at the end - do this here so we know
 !  the first filename and ndimV, labels etc.
 !
       if (doanalysis) then
          call open_analysis(outformat,required,ncolumns+ncalc,ndim,ndimV)
       endif
    else
       call get_data(ifile,.true.)
    endif
    !--dump each step in file to an output file
    do idump = 1,nstepsinfile(ifile)
       if (idump.gt.1) then
          write(filename,"(a,'_',i3.3)") rootname(ifile),idump
       else
          filename = trim(rootname(ifile))
       endif
       ntotal = sum(npartoftype(1:ntypes,idump))
       if (doanalysis) then
          call write_analysis(time(idump),dat(1:ntotal,:,idump),ntotal,ntypes, &
                          npartoftype(1:ntypes,idump),masstype(1:ntypes,idump),iamtype(:,idump), &
                          ncolumns+ncalc,ndim,ndimV,outformat)
       elseif (converttogrid) then
          call convert_to_grid(time(idump),dat(:,:,idump),ntypes,&
                               npartoftype(1:ntypes,idump),masstype(1:ntypes,idump),iamtype(:,idump), &
                               ncolumns+ncalc,filename,outformat,useall)
       else
          call write_sphdump(time(idump),gamma(idump),dat(1:ntotal,1:ncolumns+ncalc,idump),ntotal,ntypes, &
                          npartoftype(1:ntypes,idump),masstype(1:ntypes,idump),iamtype(:,idump), &
                          ncolumns+ncalc,filename,outformat)
       endif
    enddo
 enddo

 !--for analysis we need to start and end differently
 if (doanalysis) then
    call close_analysis(outformat)
    write(string,"(i10)") nfiles
    print "(/,5('-'),a,i5,a,/)",'> FINISHED CALCULATING '//trim(ucase(outformat))//' (USED '//trim(adjustl(string))//' DUMP FILES)'
 else
    print "(/,5('-'),a,/)",'> FINISHED CONVERTING DUMP FILES '
 endif
 return
end subroutine convert_all

end module convert
