!
! subroutines to do with setting of plot limits from data
!
module limits
 use params
 implicit none
 real, dimension(maxplot,2) :: lim

contains
!
!--set plot limits for all columns
!
!  NB: does not differentiate between particle types at the moment
!
subroutine set_limits(ifromstep,itostep,ifromcol,itocol)
  use labels
  use particle_data
  implicit none
  integer, intent(in) :: ifromstep,itostep,ifromcol,itocol
  integer :: i,j,k

  print*,'setting plot limits...'
  !!--find limits of particle properties	  
  lim(:,1) = 1.e6
  lim(:,2) = -1.e6
  do i=ifromstep,itostep
     do j=ifromcol,itocol
        do k=1,ntot(i)
           lim(j,1) = min(lim(j,1),dat(k,j,i))
           lim(j,2) = max(lim(j,2),dat(k,j,i))
        enddo
     enddo
  enddo
  !
  !--warn if limits are the same
  ! 
  do j=ifromcol,itocol
     if (lim(j,2).eq.lim(j,1)) then
        print*,j,label(j),' min = max = ',lim(j,1)
     endif  
  enddo
  print*,'plot limits set'

end subroutine set_limits
!
!--save plot limits for all columns to a file
!
subroutine save_limits
  use filenames
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
!
!--read plot limits for all columns from a file
!
subroutine read_limits(ierr)
  use filenames
  use labels
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

end module limits
