!-------------------------------------------------------------------
!
!  Produces a vector plot on particle data.
!
!  At the moment, just does a simple averaging by binning particles
!  and taking the average of vx,vy in the cell to give a vector for
!  that cell.
!
!------------------------------------------------------------------

subroutine vectorplot(x,y,xmin,ymin,dx,vecx,vecy,vecmax,ntot,npixx,npixy)
  implicit none
  integer, intent(in) :: ntot,npixx,npixy
  integer i,j,k,ix,iy
  integer, dimension(npixx,npixy) :: ihoc,numcell
  integer, dimension(ntot) :: ll
  real, intent(in), dimension(ntot) :: x,y,vecx,vecy
  real, intent(in) :: vecmax
  real, dimension(npixx,npixy) :: vecpixx,vecpixy
  real, intent(in) :: xmin,ymin,dx
  real :: scale,zero,xmingrid,vmax
  real :: trans(6)

  zero= 0.0E0
  !
  !--set up grid for rendering
  !
  xmingrid = xmin - 0.000001

  trans(1) = xmin-0.5*dx		! this is for the PGVECT call
  trans(2) = dx			! see help for PGIMAG/PGGRAY/PGCONT
  trans(3) = 0.0
  trans(4) = ymin-0.5*dx
  trans(5) = 0.0
  trans(6) = dx
  !
  !--interpolation is to a coarser grid, so just average      
  !
  !  bin particles into cells using a link list 
  !     
  ihoc(:,:) = -1   ! head of chain
  numcell(:,:) = 0
  do i=1,ntot
     ix = int((x(i)-xmin)/dx)+1
     iy = int((y(i)-ymin)/dx)+1
     if ((ix.ge.1).and.(ix.le.npixx).and.(iy.ge.1).and.(iy.le.npixy)) then
        ll(i)=ihoc(ix,iy)   ! set link list of this particle to old head of list
        ihoc(ix,iy) = i	    ! set head of chain to this particle
     endif
  enddo
  !
  !--add up total vx,vy in each cell
  !
  vecpixx(:,:) = 0.
  vecpixy(:,:) = 0.
  do j=1,npixy
     do i=1,npixx
        k = ihoc(i,j)
        do while (k.ne.-1)
           vecpixx(i,j) = vecpixx(i,j) + vecx(k)
           vecpixy(i,j) = vecpixy(i,j) + vecy(k)
           numcell(i,j) = numcell(i,j) + 1
           k = ll(k)
        enddo
     enddo
  enddo
  !
  !--divide by number of particles in that cell to get average vx,vy
  !
  do j=1,npixy
     do i=1,npixx
        if (numcell(i,j).ne.0) then
           vecpixx(i,j) = vecpixx(i,j)/float(numcell(i,j))
           vecpixy(i,j) = vecpixy(i,j)/float(numcell(i,j))         
        endif
     enddo
  enddo

  call pgsah(2,45.0,0.7)   ! arrow style
  call pgsch(0.3)	  ! size of arrow head
  if (vecmax.le.0.0) then  ! adaptive limits
     scale = 0.0
     vmax = max(maxval(vecpixx(:,:)),maxval(vecpixy(:,:)))
     if (vmax.gt.0.) scale = 0.1/vmax
  else
     scale=0.1/vecmax
  endif
  print*,'vector map: pixels = ',npixx,npixy,' scale = ',scale
  call pgvect(vecpixx(:,:),vecpixy(:,:),npixx,npixy, &
       1,npixx,1,npixy,scale,0,trans,-1000.0)
  call pgsch(1.0)

  return

end subroutine vectorplot
