!------------------------------------------------
! reads an exact solution from a file
!
! file should contain two columns containing
! x axis data and y axis data
!
! this is plotted as a line on the chosen graph
!------------------------------------------------
module exactfromfile
 implicit none

contains

subroutine exact_fromfile(filename,xexact,yexact,iexactpts,ierr)
  implicit none
  character(len=*), intent(in) :: filename
  real, intent(out), dimension(:) :: xexact, yexact
  integer, intent(out) :: iexactpts, ierr
  integer :: i

  ierr = 0
  open(unit=33,file=filename,iostat=ierr,status='old',form='formatted')
  if (ierr /= 0) then
     ierr = 1
     print*,'error opening ',filename
     return
  endif
  read(33,*,end=10,err=20) (xexact(i), yexact(i), i=1,size(xexact))
  print*,'WARNING: reached array limits in ',trim(filename),': partial solution read'
  ierr = -1
  return
10 continue
  iexactpts = i-1
  print*,'finished reading ',trim(filename),' : ',iexactpts,' read'
  return
20 print*,'error reading ',trim(filename),': partial solution read'
  iexactpts = i - 1
  ierr = -2
  return

end subroutine exact_fromfile

end module exactfromfile
