module titles
 implicit none
 integer, parameter, private :: maxtitles = 50
 integer, parameter, private :: maxsteplegend = 100
 character(len=60), dimension(maxtitles), public :: pagetitles
 character(len=60), dimension(maxsteplegend), public :: steplegend
 public :: read_titles, read_steplegend
 
 private
 
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
!--reads a list of labels (one per line) to be used in the timestep legend
! (ie. for multiple timesteps on same page)
!
subroutine read_steplegend(nsteptitles)
  implicit none
  integer, intent(out) :: nsteptitles
  integer :: i
  character(len=50) :: legendfile

  legendfile = 'legend'
  nsteptitles = 0

  open(unit=57,file=legendfile,status='old',form='formatted',ERR=997)
  print*,'reading legend text from file ''',trim(legendfile),''''
  do i=1,maxtitles
     read(57,"(a)",err=998,end=66) steplegend(i)
  enddo
  print*,'WARNING: array limits reached read ',maxtitles,' titles'
  nsteptitles = maxtitles
  close(unit=56)
  return
66 continue
  nsteptitles = i-1
  close(unit=57)

  return

997 continue  ! legend file does not exist, so do nothing and return
  return
998 continue
  print*,'*** error reading legend file : at line ',i-1
  nsteptitles = i-1
  close(unit=57)
  return

end subroutine read_steplegend

end module titles
