!
!  wrapper for the main data read
!  ensures that same procedure occurs on initial read as from menu option
!
!  drives reading of all files listed on command line
!
!  Arguments:
!   ireadfile : if < 0, reads from all files
!               if > 0, reads only from the filename rootname(ireadfile)
!               if = 0, no data read, just call labelling and exact_params
!
subroutine get_data(ireadfile)
  use exact
  use filenames
  use settings
  use particle_data
  use prompting
  implicit none
  integer, intent(in) :: ireadfile
  integer :: i,istart,ifinish,ierr

  if (.not.ihavereadfilename) then
     call prompt('Enter filename to read',rootname(1))
     nfiles = 1
  endif
  ihavereadfilename = .false.
  ifinish = maxstep
  !
  !--set everything to zero initially
  !
  ndim = 0
  ndimV = 0
  ncolumns = 0
  ncalc = 0
  n_end = 0
  istart = 1

  if (ireadfile.le.0) then
     !
     !--read all steps from the data file
     !
     nstepsinfile(1:nfiles) = 0
     do i=1,nfiles
        ifinish = maxstep
        call read_data(rootname(i),istart,ifinish)
        istart = ifinish + 1 ! number of next step in data array
        nstepsinfile(i) = ifinish - istart + 1
	print*,'nsteps in file = ',nstepsinfile(i)
     enddo
     nstart = 1
     n_end = ifinish
     nstepstotal = n_end
     !
     !--read plot limits from file, otherwise set plot limits
     !
     call read_limits(ierr)
     if (ierr.gt.0) call set_limits

  elseif (ireadfile.gt.0) then
     !
     !--read from a single file only
     !
     nstepsinfile(ireadfile) = 0
     print*,'reading single step'
     ifinish = maxstep
     call read_data(rootname(ireadfile),istart,ifinish)
     nstepsinfile(ireadfile) = ifinish - istart + 1
     print*,'nsteps in file = ',nstepsinfile(ireadfile)
     !
     !--assume there are the same number of steps in the other files
     !  which have not been read
     !
     do i=1,nfiles
        if (nstepsinfile(i).eq.0) then
	   nstepsinfile(i) = nstepsinfile(ireadfile)
	endif
     enddo
     nstart = 1
     n_end = sum(nstepsinfile(1:nfiles))
     nstepstotal = n_end
     print*,'nend = ',n_end
  endif
  !
  !--set labels for each column of data
  !
  print*,'setting plot labels...'
  call set_labels
  
  numplot = ncolumns

  if (ivegotdata) then
     !
     !--calculate various additional quantities
     !     
     call calc_quantities
  endif
  !
  !--read exact solution parameters from files if present
  !
  if (iexact.ne.0) call read_exactparams(iexact,ierr)
  
  return
end subroutine get_data
