!
!--reads a list of titles (one per line), to be used to label each timestep
!
subroutine read_titles(titlelist,ntitles,maxtitles)
  use filenames
  implicit none
  integer, intent(in) :: maxtitles
  integer, intent(out) :: ntitles
  integer :: i
  character(len=*), dimension(maxtitles), intent(out) :: titlelist
  character(len=50) :: titlefile

  titlefile = 'titlelist'  !!trim(rootname(1))//'.titles'

  open(unit=56,file=titlefile,status='old',form='formatted',ERR=997)
  print*,'reading plot titles from file ',trim(titlefile)
  do i=1,maxtitles
     read(56,"(a)",err=998,end=66) titlelist(i)
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
