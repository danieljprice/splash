!
!  wrapper for the main data read
!  ensures that same procedure occurs on initial read as from menu option
!
subroutine get_data   
  use exact_params
  use filenames
  use settings
  use particle_data
  use prompting
  implicit none
  integer :: i,istep,ierr
  character(LEN=30) :: filename

  if (.not.ihavereadfilename) then
     call prompt('Enter filename to read',rootname(1))
     nfiles = 1
  endif
  ihavereadfilename = .false.
  nfilesteps = maxstep
  !
  !--set everything to zero initially
  !
  hfact = 0.
  ndim = 0
  ndimV = 0
  ncolumns = 0
  ncalc = 0
  nfilesteps = 0
  n_end = 0
  istep = 1
  !
  !--read the data from the file
  !
  do i=1,nfiles
     nfilesteps = maxstep
     call read_data(rootname(i),istep,nfilesteps)
     istep = nfilesteps + 1 ! current location of istep in data array
  enddo

  nstart = 1
  n_end = nfilesteps

  if (ivegotdata) then
     !
     !--calculate various additional quantities
     !     
     call calc_quantities	  
     !
     !--read plot limits from file, otherwise set plot limits
     !
     call read_limits(ierr)
     if (ierr.gt.0) call set_limits
  endif
  !
  !--read exact solution parameters from files if present
  !
  if (iexact.ne.0) call read_exactparams(iexact,ierr)
  
  return
end subroutine get_data
