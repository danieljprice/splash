!
!--set plot limits for all columns
!
subroutine set_limits
  use labels
  use limits
  use particle_data
  use settings
  implicit none
  integer :: i,j,k

  print*,'setting plot limits...'
  ntotplot(nstart:n_end) = npart(nstart:n_end)
  if (iplotghost) ntotplot(nstart:n_end) = ntot(nstart:n_end)
  zoom = 1.0
  !!--find limits of particle properties	  
  lim(:,1) = 1.e6
  lim(:,2) = -1.e6
  do j=1,numplot
     do i=nstart,n_end
        do k=1,ntotplot(i)
           lim(j,1) = min(lim(j,1),dat(j,k,i))
           lim(j,2) = max(lim(j,2),dat(j,k,i)*scalemax)
        enddo
     enddo
     if (lim(j,2).eq.lim(j,1)) then
        print*,label(j),' min = max = ',lim(j,1)
     endif
  enddo
  !!--limits of magnetic field (in all dimensions) for vector plot
  if (iBfirst.ne.0) then
     if (iadapt) then
        Bmin = 0.0
        Bmax = maxval(dat(iBfirst:iBlast,1:ntotplot(1),1))*scalemax
     else
        Bmax = 0.0
        Bmin = 0.0
        do j=iBfirst,iBlast
           !Bmin = min(Bmin,lim(j,1))
           if (lim(j,2).gt.Bmax) then
              Bmax=lim(j,2)
              print*,' Bmax = ',label(j),lim(j,2)
           endif
        enddo
     endif
     print*,'Bmin,Bmax = ',Bmin,Bmax
  endif
  print*,'plot limits set'
  !

end subroutine set_limits
