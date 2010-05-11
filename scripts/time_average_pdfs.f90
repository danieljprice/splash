!--------------------------------------------------------------
!  Utility to get the time averaged PDF from the PDF
!  files produced by SPLASH
!--------------------------------------------------------------
program time_average_pdfs
 implicit none
 integer, parameter :: iunit = 10, iout = 7, iprint = 6
 integer, parameter :: maxbins = 2000, maxfiles= 200
 integer :: nargs,iarg,nbins,nbinsprev,nfiles,i,ierr
 character(len=120) :: filename,line
 real, dimension(maxbins) :: xval,xvalprev,ave,var
 real, dimension(maxbins,maxfiles) :: pdfval

 nargs = command_argument_count()
 if (nargs.le.0) then
    print "(a)",'A utility for computing the time averaged Probability Distribution Function'
    print "(a)",'from a bunch of PDF files as output by SPLASH'
    print "(/,a)",'Usage: time_average_pdfs pdffiles'
    stop
 elseif (nargs.gt.maxfiles) then
    print "(a)",' ERROR: number of files exceeds array limits, edit and recompile'
    stop
 endif
 
 nbinsprev = 0
 xvalprev = 0.
 ave = 0.
 pdfval = 0.
 xval = 0.
 nfiles = 0
 
 over_files: do iarg=1,nargs
    call get_command_argument(iarg,filename)
    
    open(unit=iunit,file=filename,iostat=ierr,status='old')
    if (ierr /= 0) then
       print "(a)",' ERROR: cannot open '//trim(filename)//' for read'
    else
       nbins = 0
       do while (ierr.eq.0 .and. nbins.lt.maxbins)
          read(iunit,"(a)",iostat=ierr) line
          if (ierr.eq.0 .and. adjustl(line(1:1)).ne.'#') then
             nbins = nbins + 1
             read(line,*,iostat=ierr) xval(nbins),pdfval(nbins,nfiles+1)
          endif
       enddo
       if (ierr.eq.0) print "(a,i4,a,i4,a)",' ERROR! number of bins ',nbins,' exceeds maximum (',maxbins,')'
       if (all(pdfval(1:nbins,nfiles+1).lt.tiny(0.))) then
          write(iprint,*) 'skipping '//trim(filename)//': PDF = 0'
          cycle over_files
       else
          write(iprint,*) trim(filename)//' nbins = ',nbins
       endif       
       !--error checks
       if (nbins.le.0) then
          print "(a)",' ERROR: no data read from file, skipping'
       endif
       
       if (nbinsprev.gt.0) then
          if (nbins.ne.nbinsprev) then
             print "(a)",' ERROR: number of bins has changed between files, skipping file'
             cycle over_files
          elseif(.not.all(abs(xval(1:nbins)-xvalprev(1:nbins)).lt.1.e-6)) then
             do i=1,nbins
                if (abs(xval(i)-xvalprev(i)).gt.1.e-6) print*,i,xval(i),xvalprev(i)
             enddo
             print "(a)",' ERROR: location of bins has changed between files, skipping file'
             cycle over_files
          endif
       else
          nbinsprev = nbins
          xvalprev = xval
       endif
       nfiles = nfiles + 1
       ave(1:nbins) = ave(1:nbins) + pdfval(1:nbins,nfiles)
    endif
 enddo over_files
 
 !--compute average
 if (nfiles.gt.0) then
    write(iprint,*) 'last file: '//trim(filename)//' nfiles = ',nfiles
    ave(1:nbins) = ave(1:nbins)/real(nfiles)
 else
    print*,' ERROR: nfiles = ',nfiles
    stop
 endif
 
 !--compute standard deviation
 var = 0.
 do iarg=1,nfiles
    var(1:nbins) = var(1:nbins) + (pdfval(1:nbins,iarg) -  ave(1:nbins))**2
 enddo
 var(1:nbins) = var(1:nbins)/nfiles
 
 open(unit=iout,file='averaged_pdf.dat',status='replace',form='formatted')
 write(iout,"(a)") '# [01  xval ] [02  PDF(average)] [03 st. dev] [04 variance]'
 do i=1,nbins
    if (ave(i).gt.0.) then ! only spit out non-zero bins
       write(iout,"(4(1pe10.4,2x))") xval(i),ave(i),sqrt(var(i)),var(i)
    endif
 enddo
 close(iout)

end program time_average_pdfs
