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
  !
  !--calculate various additional quantities
  !     
  call calc_quantities	  
  !
  !--read plot limits from file, otherwise set plot limits
  !
  call read_limits(ierr)
  if (ierr.gt.0) call set_limits
  !
  !--read toy star file for toy star solution
  !
  if (iexact.eq.4) then
     filename = trim(rootname(1))//'.tstar'
     open(UNIT=20,ERR=8801,FILE=filename,STATUS='old')
     read(20,*,ERR=8801) Htstar,Ctstar,Atstar
     read(20,*,ERR=8801) sigma0
     read(20,*,ERR=8801) norder
     close(UNIT=20)
     print*,' >> read ',filename
     print*,' H,C,A,sigma,n = ',Htstar,Ctstar,Atstar,sigma0,norder
     return
8801 continue
     print*,'Cannot open/ error reading ',filename	     
  endif
  
  return
end subroutine get_data
