!-----------------------------------------------------------------------
! read exact solution parameters from files
! (in ndspmhd these files are used in the input to the code)
!
! called after main data read and if exact solution chosen from menu
!-----------------------------------------------------------------------
subroutine read_exactparams(iexact,ierr)
  use exact_params
  use filenames
  implicit none
  integer, intent(in) :: iexact
  integer, intent(out) :: ierr
  integer :: int_from_string
  character(LEN=30) :: filename

  select case(iexact)
  case(1)
  !
  !--shock tube parameters from .shk file
  !
     filename = trim(rootname(1))//'.shk'
     open(UNIT=19,ERR=7701,FILE=filename,STATUS='old')
     read(19,*,ERR=7777) rho_L, rho_R
     read(19,*,ERR=7777) pr_L, pr_R
     read(19,*,ERR=7777) v_L, v_R
     close(UNIT=19)
     print*,'>> read ',filename
     print*,' rhoL, rho_R = ',rho_L,rho_R
     print*,' pr_L, pr_R  = ',pr_L, pr_R
     print*,' v_L,  v_R   = ',v_L, v_R
     return
7701 print*,'no file ',filename
     ierr = 1
     return
7777 print*,'error reading ',filename
     close(UNIT=19)
     ierr = 2
     return

  case(4)
  !
  !--read toy star file for toy star solution
  !
     filename = trim(rootname(1))//'.tstar'
     open(unit=20,ERR=8801,FILE=filename,STATUS='old')
     read(20,*,ERR=8888) Htstar,Ctstar,Atstar
     read(20,*,ERR=8888) sigma0
     read(20,*,ERR=8888) norder
     close(UNIT=20)
     print*,' >> read ',filename
     print*,' H,C,A,sigma,n = ',Htstar,Ctstar,Atstar,sigma0,norder
     return
8801 continue
     print*,'no file ',filename
     ierr = 1
     return
8888 print*,'error reading ',filename
     close(UNIT=20)
     ierr = 2
     return

  case(6)
  !
  !--attempt to guess which MHD shock tube has been done from filename
  !
     ishk = int_from_string(rootname(1)(5:5))
     return

  end select
  
  return
end subroutine read_exactparams
