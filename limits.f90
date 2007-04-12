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
  use labels, only:label
  use geometry, only:coord_transform_limits
  use particle_data, only:npartoftype,dat,maxcol
  use settings_data, only:ndim,icoords,icoordsnew
  implicit none
  integer, intent(in) :: ifromstep,itostep,ifromcol,itocol
  integer :: i,j,k,ntoti

  print 100,ifromstep,itostep,ifromcol,itocol
100 format(/' setting plot limits: steps ',i5,'->',i5,' cols ',i2,'->',i3)
  if (ifromcol.gt.maxcol .or. itocol.gt.maxcol) then
     print "(a)",' *** error: set_limits: column > array size ***'
     return
  endif

  !!--find limits of particle properties	  
  lim(ifromcol:itocol,1) = huge(lim)
  lim(ifromcol:itocol,2) = -huge(lim)
  do i=ifromstep,itostep
     ntoti = sum(npartoftype(:,i))
     do j=ifromcol,itocol
        do k=1,ntoti
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
        print "(a,a20,a,1pe9.2)",'  warning: ',label(j),' min = max = ',lim(j,1)
     endif  
  enddo
  print "(a/)",' plot limits set'
  
  !
  !--transform coord limits into new coordinate system if coord transform is applied
  !
  if (icoordsnew.ne.icoords .and. ndim.gt.0 .and. ifromcol.le.ndim) then
     call coord_transform_limits(lim(1:ndim,1),lim(1:ndim,2), &
                                 icoords,icoordsnew,ndim)
  endif
  
  return
end subroutine set_limits
!
!--save plot limits for all columns to a file
!
subroutine write_limits(limitsfile)
  use settings_data, only:numplot
  use prompting, only:prompt
  implicit none
  character(len=*), intent(in) :: limitsfile
  integer :: i

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

end subroutine write_limits
!
!--read plot limits for all columns from a file
!
subroutine read_limits(limitsfile,ierr)
  use labels, only:label
  use settings_data, only:numplot,ncolumns,ncalc
  use prompting, only:prompt
  implicit none
  character(len=*), intent(in) :: limitsfile
  integer, intent(out) :: ierr
  integer :: i

  ierr = 0

  open(unit=54,file=limitsfile,status='old',form='formatted',err=997)
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
  !--only give error if we really do not have enough columns
  !  (on first call nextra is not set)
  if (i.le.ncolumns+ncalc) then
     print*,'end of file in ',trim(limitsfile),': limits read to column ',i
     ierr = -1
  endif
  close(unit=54)
  return

end subroutine read_limits

end module limits
