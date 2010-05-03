program grid2pdf
 use system_commands,    only:get_number_arguments,get_argument
 use readwrite_griddata, only:open_gridfile_r,read_gridcolumn
 use pdfs,               only:pdf_calc,pdf_write
 implicit none
 logical, parameter     :: usefixedbins = .true.
 integer, parameter     :: nbins = 270
 integer, parameter     :: iunit = 13
 integer, dimension(3)  :: nxgrid
 integer                :: ierr,ncolumns,itransx,ifile,nargs
 integer                :: ntot
 real                   :: time,rhologmin,rhologmax,pdfmin,pdfmax
 real, dimension(nbins) :: xbin,pdf
 character(len=120)     :: filename,tagline
 character(len=20)      :: informat
 logical                :: volweighted
 real, dimension(:), allocatable :: rhogrid
 
 tagline = 'grid2pdf: a SPLASH utility (c) 2010 Daniel Price'
 
 call get_number_arguments(nargs)
 if (nargs.le.0) then
    print "(a)",trim(tagline)
    print "(/,a,/)",'Usage: grid2pdf gridfile(s)'
    stop
 endif
 
 do ifile=1,nargs
    !
    !--get the filename off the command line
    !
    call get_argument(ifile,filename)
    !
    !--open the grid data file and read the header information
    !
    informat = 'gridbinary'
    call open_gridfile_r(iunit,filename,informat,nxgrid(:),ncolumns,time,ierr)
    !
    !--allocate memory for the grid data
    !
    ntot = product(nxgrid)
    if (.not.allocated(rhogrid) .or. ntot.ne.size(rhogrid)) then
       if (allocated(rhogrid)) deallocate(rhogrid)
       allocate(rhogrid(ntot))
    endif
    !
    !--read one particular column (in this case, the density)
    !
    call read_gridcolumn(iunit,rhogrid,ntot,ierr)
    !
    !--close the file
    !
    close(iunit)
    !
    !--set the parameters for calculating the PDF
    !  [in this case, we want the PDF of ln (rho)]
    !
    rhologmin = -15.
    rhologmax = 12.
    where (rhogrid.gt.0.)
       rhogrid = log(rhogrid)
    end where
    itransx = 6
    !
    !--calculate the PDF
    !
    call pdf_calc(ntot,rhogrid,rhologmin,rhologmax,nbins,xbin,&
                  pdf,pdfmin,pdfmax,itransx,usefixedbins,volweighted,ierr)
    !
    !--write the PDF to file
    !
    if (ierr.eq.0) then
       call pdf_write(nbins,xbin,pdf,'rhogrid',itransx,volweighted,trim(filename),trim(tagline))
    endif
 enddo
 !
 !--clean up/deallocate memory
 !
 if (allocated(rhogrid)) deallocate(rhogrid)

end program grid2pdf
