!
!--save plot limits for all columns to a file
!
subroutine save_limits
  use filenames
  use limits
  use settings_data
  implicit none
  integer :: i
  character(len=26) :: limitsfile

  limitsfile = trim(rootname(1))//'.limits'
  print*,'saving plot limits to file ',trim(limitsfile)

  open(unit=55,file=limitsfile,status='replace',form='formatted',ERR=998)
  do i=1,numplot
     write(55,*,err=999) lim(i,1),lim(i,2)
  enddo
  close(unit=55)

  return

998 continue
  print*,'*** error opening limits file: limits not saved'
  return
999 continue
  print*,'*** error saving limits'
  close(unit=55)
  return

end subroutine save_limits
