!
!--set plot limits for all columns
!
!  NB: does not differentiate between particle types at the moment
!
subroutine set_limits(ifromstep,itostep,ifromcol,itocol)
  use labels
  use limits
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
