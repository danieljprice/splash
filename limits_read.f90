!
!--read plot limits for all columns from a file
!
subroutine read_limits(ierr)
  use filenames
  use labels
  use limits
  use settings_data
  implicit none
  integer :: i,ierr
  character(len=26) :: limitsfile

  ierr = 0
  limitsfile = trim(rootname(1))//'.limits'

  open(unit=54,file=limitsfile,status='old',form='formatted',ERR=997)
  print*,'reading plot limits from file ',trim(limitsfile)
  do i=1,numplot
     read(54,*,err=998,end=999) lim(i,1),lim(i,2)
     if (lim(i,1).eq.lim(i,2)) then
        print*,label(i),' min = max = ',lim(i,1)
     endif
  enddo
  close(unit=54)

  return

997 continue
  print*,trim(limitsfile),' not found'
  ierr = 1
  return
998 continue
  print*,'*** error reading limits from file'
  ierr = 2
  close(unit=54)
  return
999 continue
  print*,'end of file in ',trim(limitsfile),': limits read to column ',i
  ierr = -1
  close(unit=54)
  return

end subroutine read_limits
