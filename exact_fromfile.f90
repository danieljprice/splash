!------------------------------------------------
! reads an exact solution from a file
!
! file should contain two columns containing
! x axis data and y axis data
!
! this is plotted as a line on the chosen graph
!------------------------------------------------
subroutine exact_fromfile(filename,ierr)
  use exact_params
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ierr
  integer :: i

  open(unit=33,file=filename,err=20,status='old',form='formatted')
  read(33,*,END=10,ERR=30) (xexact(i), yexact(i), i=1,maxexactpts)
  print*,'WARNING: reached array limits in ',trim(filename),': partial solution read'
  ierr = -1
  return
10 continue
  iexactpts = i-1
  print*,'finished reading ',trim(filename),' : ',iexactpts,' read'
  return
20 print*,'error opening ',filename
  ierr = 1
  return
30 print*,'error reading ',trim(filename),': partial solution read'
  iexactpts = i - 1
  ierr = -2
  return

end subroutine exact_fromfile
