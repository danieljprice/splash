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
!     this version is written for slices through a rectangular volume, ie.
!     assumes a uniform pixel size in x,y, whilst the number of pixels
!     in the z direction can be set to the number of cross-section slices.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!--------------------------------------------------------------------------

subroutine interpolate3D(x,y,z,pmass,rho,hh,dat,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth,zpixwidth)

  implicit none
  real, parameter :: pi = 3.1415926536      
  integer, intent(IN) :: npart,npixx,npixy,npixz
  real, intent(IN), dimension(npart) :: x,y,z,pmass,rho,hh,dat
  real, intent(IN) :: xmin,ymin,zmin,pixwidth,zpixwidth
  real, intent(OUT), dimension(npixx,npixy,npixz) :: datsmooth

  integer :: i,j,ipix,jpix,kpix
  integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
  real :: hi,hi1,h3,radkern,qq,wab,rab,const
  real :: term,dx,dy,dz,xpix,ypix,zpix

  datsmooth = 0.
  term = 0.
  print*,'interpolating from particles to 3D grid...'
   if (pixwidth.le.0.) then
     print*,'interpolate3D: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
  do i=1,npart
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     h3 = hi*hi*hi
     radkern = 2.*hi   ! radius of the smoothing kernel
     const = 1./(pi*h3)  ! normalisation constant (3D)
     term = 0.
     if (rho(i).ne.0.) term = pmass(i)*dat(i)/rho(i) 
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     kpixmin = int((z(i) - radkern - zmin)/zpixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)
     jpixmax = int((y(i) + radkern - ymin)/pixwidth)
     kpixmax = int((z(i) + radkern - zmin)/zpixwidth)

     !         PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     !         PRINT*,'z slices = ',kpixmin,zmin + kpixmin*zpixwidth, !- 0.5*zpixwidth,
     !     &                 kpixmax,zmin + kpixmax*zpixwidth                !- 0.5*zpixwidth
     !        PRINT*,'should cover z = ',z(i)-radkern,' to ',z(i)+radkern         
     !         PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
     !         PRINT*,'xmin,ymin = ',xmin,ymin,zmin

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (kpixmin.lt.1) kpixmin = 1
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     if (kpixmax.gt.npixz) kpixmax = npixz
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do kpix = kpixmin,kpixmax
        zpix = zmin + (kpix)*zpixwidth  !- 0.5*zpixwidth
        dz = zpix - z(i)
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix)*pixwidth  !- 0.5*pixwidth
           dy = ypix - y(i)  
           do ipix = ipixmin,ipixmax
              xpix = xmin + (ipix)*pixwidth   !- 0.5*pixwidth
              dx = xpix - x(i)
              rab = sqrt(dx**2 + dy**2 + dz**2)
              qq = rab*hi1
              !
              !--SPH kernel - standard cubic spline
              !                     
              if (qq.lt.1.0) then
                 wab = const*(1.-1.5*qq**2 + 0.75*qq**3)
              elseif (qq.lt.2.0) then
                 wab = const*0.25*(2.-qq)**3
              else
                 wab = 0.
              endif
              !
              !--calculate data value at this pixel using the summation interpolant
              !  
              datsmooth(ipix,jpix,kpix) = datsmooth(ipix,jpix,kpix) + term*wab          

           enddo
        enddo
     enddo
  enddo

  return

end subroutine interpolate3D
