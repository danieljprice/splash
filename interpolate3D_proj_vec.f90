!--------------------------------------------------------------------------
!     program to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a vector field integrated along the line of sight.
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price 23/12/04
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_vec(x,y,pmass,rho,hh,vecx,vecy,npart,&
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth)

  use column
  implicit none
  real, parameter :: pi = 3.1415926536      
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,pmass,rho,hh,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: index,index1
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: rho1i,termx,termy,dx,dy,xpix,ypix
  real :: dxx,dwdx,dmaxcoltable

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  dmaxcoltable = 1./real(maxcoltable)
  print*,'projecting vector from particles to pixels...'
  if (pixwidth.le.0.) then
     print*,'interpolate3D_proj_vec: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
  over_particles: do i=1,npart
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D_proj_vec: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     radkern = 2.*hi    ! radius of the smoothing kernel
     const = hi1*hi1
     if (rho(i).ne.0.) then
        rho1i = 1./rho(i)
     else
        rho1i = 0.
     endif
        
     termx = const*pmass(i)*vecx(i)*rho1i
     termy = const*pmass(i)*vecy(i)*rho1i
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)
     jpixmax = int((y(i) + radkern - ymin)/pixwidth)

     ! PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     ! PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
           dx = xpix - x(i)
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !                     
           !!--find nearest index in table
           index = int(0.5*qq*maxcoltable) + 1
           if (index.lt.maxcoltable) then
              index1 = index + 1
              !--find increment along from this index
              dxx = qq - index*maxcoltable
              !--find gradient                  
              dwdx = (coltable(index1) - coltable(index))*dmaxcoltable
              !--compute value of integrated kernel
              wab = coltable(index) + dwdx*dxx
              !
              !--calculate data value at this pixel using the summation interpolant
              !  
              vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab          
              vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab          

           endif

        enddo
     enddo

  enddo over_particles

  return

end subroutine interpolate3D_proj_vec
