!
!--set plot limits for all columns
!
!  NB: does not differentiate between particle types at the moment
!
subroutine set_limits
  use labels
  use limits
  use particle_data
  use settings
  implicit none
  integer :: i,j,k

  print*,'setting plot limits...'
  zoom = 1.0
  !!--find limits of particle properties	  
  lim(:,1) = 1.e6
  lim(:,2) = -1.e6
  do i=nstart,n_end
     do j=1,numplot
        do k=1,ntot(i)
           lim(j,1) = min(lim(j,1),dat(k,j,i))
           lim(j,2) = max(lim(j,2),dat(k,j,i)*scalemax)
        enddo
     enddo
  enddo
  !
  !--warn if limits are the same
  ! 
  do j=1,numplot
     if (lim(j,2).eq.lim(j,1)) then
        print*,j,label(j),' min = max = ',lim(j,1)
     endif  
  enddo
  print*,'plot limits set'

end subroutine set_limits
