!----------------------------------------------------------------------
!
!  Module containing all of the routines required for cross sections
!  through 3D data
!
!----------------------------------------------------------------------

module xsections3D
 implicit none
 real, parameter, private :: pi = 3.1415926536      
 public :: interpolate3D, interpolate3D_fastxsec, interpolate3D_xsec_vec
 
contains

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
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
  integer, intent(in) :: npart,npixx,npixy,npixz
  real, intent(in), dimension(npart) :: x,y,z,pmass,rho,hh,dat
  real, intent(in) :: xmin,ymin,zmin,pixwidth,zpixwidth
  real, intent(out), dimension(npixx,npixy,npixz) :: datsmooth

  integer :: i,ipix,jpix,kpix
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

!--------------------------------------------------------------------------
!
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!            cross section location: zslice
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------

subroutine interpolate3D_fastxsec(x,y,z,pmass,rho,hh,dat,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,pmass,rho,hh,dat
  real, intent(in) :: xmin,ymin,pixwidth,zslice
  real, intent(out), dimension(npixx,npixy) :: datsmooth

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,h3,radkern,qq,wab,rab,const
  real :: term,dx,dy,dz,dz2,xpix,ypix

  datsmooth = 0.
  term = 0.
  print*,'taking fast cross section...',zslice
  if (pixwidth.le.0.) then
     print*,'interpolate3D_xsec: error: pixel width <= 0'
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
        print*,'interpolate3D_xsec: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     h3 = hi*hi*hi
     radkern = 2.*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (abs(dz) .lt. radkern) then

        const = 1./(pi*h3)  ! normalisation constant (3D)
        term = 0.
        if (rho(i).ne.0.) term = const*pmass(i)*dat(i)/rho(i) 
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
              rab = sqrt(dx**2 + dy**2 + dz2)
              qq = rab*hi1
              !
              !--SPH kernel - standard cubic spline
              !     
              if (qq.lt.2.0) then
                 if (qq.lt.1.0) then
                    wab = (1.-1.5*qq**2 + 0.75*qq**3)
                 else
                    wab = 0.25*(2.-qq)**3
                 endif
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !  
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab          

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo                     ! over particles

  return

end subroutine interpolate3D_fastxsec

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
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!            cross section location: zslice
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------

subroutine interpolate3D_xsec_vec(x,y,z,pmass,rho,hh,vecx,vecy,npart,&
     xmin,ymin,zslice,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,pmass,rho,hh,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth,zslice
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,h3,radkern,qq,wab,rab,const
  real :: rho1i,termx,termy,dx,dy,dz,dz2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  print*,'taking fast cross section...',zslice
  if (pixwidth.le.0.) then
     print*,'interpolate3D_xsec_vec: error: pixel width <= 0'
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
        print*,'interpolate3D_xsec_vec: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     h3 = hi*hi*hi
     radkern = 2.*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (abs(dz) .lt. radkern) then

        const = 1./(pi*h3)  ! normalisation constant (3D)
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
              rab = sqrt(dx**2 + dy**2 + dz2)
              qq = rab*hi1
              !
              !--SPH kernel - standard cubic spline
              !     
              if (qq.lt.2.0) then
                 if (qq.lt.1.0) then
                    wab = (1.-1.5*qq**2 + 0.75*qq**3)
                 else
                    wab = 0.25*(2.-qq)**3
                 endif
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !  
                 vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab          
                 vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab          

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo                     ! over particles

  return

end subroutine interpolate3D_xsec_vec

end module xsections3D
