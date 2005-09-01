module titles
 implicit none
 integer, parameter :: maxtitles = 50
 integer, parameter :: maxsteptitles = 1000
 character(len=60), dimension(maxtitles), public :: pagetitles
 character(len=60), dimension(maxsteptitles), public :: steptitles
 
contains
!
!--reads a list of titles (one per line), to be used to label each plot on page
!
subroutine read_titles(ntitles)
  implicit none
  integer, intent(out) :: ntitles
  integer :: i
  character(len=50) :: titlefile

  titlefile = 'titlelist'
  ntitles = 0

  open(unit=56,file=titlefile,status='old',form='formatted',ERR=997)
  print*,'reading plot titles from file ',trim(titlefile)
  do i=1,maxtitles
     read(56,"(a)",err=998,end=66) pagetitles(i)
  enddo
  print*,'WARNING: array limits reached read ',maxtitles,' titles'
  ntitles = maxtitles
  close(unit=56)
  return
66 continue
  ntitles = i-1
  close(unit=56)

  return

997 continue  ! title file does not exist, so do nothing and return
  return
998 continue
  print*,'*** error reading title file : at line ',i-1
  ntitles = i-1
  close(unit=56)
  return

end subroutine read_titles
!
!--reads a list of titles (one per line), to be used to label each timestep
!
subroutine read_steptitles(nsteptitles)
  implicit none
  integer, intent(out) :: nsteptitles
  integer :: i
  character(len=50) :: titlefile

  titlefile = 'legend'
  nsteptitles = 0

  open(unit=57,file=titlefile,status='old',form='formatted',ERR=997)
  print*,'reading legend text from file ''',trim(titlefile),''''
  do i=1,maxtitles
     read(57,"(a)",err=998,end=66) steptitles(i)
  enddo
  print*,'WARNING: array limits reached read ',maxtitles,' titles'
  nsteptitles = maxtitles
  close(unit=56)
  return
66 continue
  nsteptitles = i-1
  close(unit=57)

  return

997 continue  ! title file does not exist, so do nothing and return
  return
998 continue
  print*,'*** error reading legend file : at line ',i-1
  nsteptitles = i-1
  close(unit=57)
  return

end subroutine read_steptitles

end module titles
