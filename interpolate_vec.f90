!--------------------------------------------------------------------------
!    Interpolates vector quantity from particles to even grid of pixels
!
!    This version just does a simple averaging by binning particles
!    and taking the average of vx,vy in the cell to give a vector for
!    that cell. This is because the interpolation of a vector quantity is
!    usually to a *coarser* grid than the particles.
!
!    Input: particle coordinates  : x,y   (npart)
!           vector data to smooth : vecx  (npart)
!                                   vecy  (npart)
!           grid setup : xmin, ymin, dx
!
!     Output: smoothed vector field   : vecpixx (npixx,npixy)
!                                     : vecpixy (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 20/8/04
!--------------------------------------------------------------------------

subroutine interpolate_vec(x,y,vecx,vecy, &
     xmin,ymin,dx,vecpixx,vecpixy,npart,npixx,npixy)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,vecx,vecy
  real, intent(in) :: xmin,ymin,dx
  real, intent(out), dimension(npixx,npixy) :: vecpixx, vecpixy
  real, parameter :: pi = 3.1415926536
  integer :: i,j,k,ix,iy
  integer, dimension(npixx,npixy) :: ihoc,numcell
  integer, dimension(npart) :: ll

  print*,'averaging vector field onto pixels...'
  if (dx.le.0.) then
     print*,'interpolate_vec: error: pixel width <= 0'
     return
  endif
  !
  !--interpolation is to a coarser grid, so just average      
  !
  !  bin particles into cells using a link list 
  !     
  ihoc(:,:) = -1   ! head of chain
  numcell(:,:) = 0
  do i=1,npart
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
  return
  
end subroutine interpolate_vec
