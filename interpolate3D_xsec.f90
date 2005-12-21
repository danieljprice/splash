!----------------------------------------------------------------------
!
!  Module containing all of the routines required for cross sections
!  through 3D data
!
!----------------------------------------------------------------------

module xsections3D
 implicit none
 real, parameter, private :: dpi = 1./3.1415926536      
 public :: interpolate3D, interpolate3D_fastxsec, interpolate3D_xsec_vec
 
contains

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is interpolated according to the formula
!
!     datsmooth(pixel) = sum_b weight_b dat_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline.
!
!     For a standard SPH smoothing the weight function for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     this version is written for slices through a rectangular volume, ie.
!     assumes a uniform pixel size in x,y, whilst the number of pixels
!     in the z direction can be set to the number of cross-section slices.
!
!     Input: particle coordinates  : x,y,z (npart)
!            smoothing lengths     : hh    (npart)
!            weights for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!--------------------------------------------------------------------------

subroutine interpolate3D(x,y,z,hh,weight,dat,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth,zpixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy,npixz
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  real, intent(in) :: xmin,ymin,zmin,pixwidth,zpixwidth
  real, intent(out), dimension(npixx,npixy,npixz) :: datsmooth
  real, dimension(npixx,npixy,npixz) :: datnorm

  integer :: i,ipix,jpix,kpix
  integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,termnorm,dx,dy,dz,xpix,ypix,zpix

  datsmooth = 0.
  datnorm = 0.
  print*,'interpolating from particles to 3D grid...'
   if (pixwidth.le.0.) then
     print*,'interpolate3D: error: pixel width <= 0'
     return
  endif
  const = dpi  ! normalisation constant (3D)
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
     radkern = 2.*hi   ! radius of the smoothing kernel
     termnorm = const*weight(i)
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     kpixmin = int((z(i) - radkern - zmin)/zpixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidth) + 1
     kpixmax = int((z(i) + radkern - zmin)/zpixwidth) + 1

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
        zpix = zmin + (kpix-0.5)*zpixwidth
        dz = zpix - z(i)
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidth
           dy = ypix - y(i)  
           do ipix = ipixmin,ipixmax
              xpix = xmin + (ipix-0.5)*pixwidth
              dx = xpix - x(i)
              rab = sqrt(dx**2 + dy**2 + dz**2)
              qq = rab*hi1
              !
              !--SPH kernel - standard cubic spline
              !                     
              if (qq.lt.1.0) then
                 wab = 1.-1.5*qq**2 + 0.75*qq**3
              elseif (qq.lt.2.0) then
                 wab = 0.25*(2.-qq)**3
              else
                 wab = 0.
              endif
              !
              !--calculate data value at this pixel using the summation interpolant
              !  
              datsmooth(ipix,jpix,kpix) = datsmooth(ipix,jpix,kpix) + term*wab          
              datnorm(ipix,jpix,kpix) = datnorm(ipix,jpix,kpix) + termnorm*wab          

           enddo
        enddo
     enddo
  enddo
  !
  !--normalise dat array
  !
  where (datnorm > 0.)
     datsmooth = datsmooth/datnorm
  end where
  
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

subroutine interpolate3D_fastxsec(x,y,z,hh,weight,dat,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  real, intent(in) :: xmin,ymin,pixwidth,zslice
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: termnorm,term,dx,dy,dz,dz2,xpix,ypix,rescalefac

  datsmooth = 0.
  datnorm = 0.
  print*,'taking fast cross section...',zslice
  if (pixwidth.le.0.) then
     print*,'interpolate3D_xsec: error: pixel width <= 0'
     return
  elseif (npart.le.0) then
     print*,'interpolate3D_xsec: error: npart = 0'
     return
  endif
  const = dpi
  !
  !--renormalise dat array by first element to speed things up
  !
  if (dat(1).gt.tiny(dat)) then
     rescalefac = dat(1)
  else
     rescalefac = 1.0
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

        termnorm = const*weight(i)
        term = termnorm*dat(i)/rescalefac
        !
        !--for each particle work out which pixels it contributes to
        !               
        ipixmin = int((x(i) - radkern - xmin)/pixwidth)
        jpixmin = int((y(i) - radkern - ymin)/pixwidth)
        ipixmax = int((x(i) + radkern - xmin)/pixwidth) + 1
        jpixmax = int((y(i) + radkern - ymin)/pixwidth) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidth
           dy = ypix - y(i)
           do ipix = ipixmin,ipixmax
              xpix = xmin + (ipix-0.5)*pixwidth
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
                 datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab      

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo                     ! over particles
  !
  !--normalise dat array
  !
  where (datnorm > 0.)
     datsmooth = datsmooth/datnorm*rescalefac
  end where
  
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

subroutine interpolate3D_xsec_vec(x,y,z,hh,weight,vecx,vecy,npart,&
     xmin,ymin,zslice,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth,zslice
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: termx,termy,termnorm,dx,dy,dz,dz2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  print*,'taking fast cross section...',zslice
  if (pixwidth.le.0.) then
     print*,'interpolate3D_xsec_vec: error: pixel width <= 0'
     return
  endif
  const = dpi ! normalisation constant (3D)
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
        termnorm = const*weight(i)
        termx = termnorm*vecx(i)
        termy = termnorm*vecy(i)
        !
        !--for each particle work out which pixels it contributes to
        !               
        ipixmin = int((x(i) - radkern - xmin)/pixwidth)
        jpixmin = int((y(i) - radkern - ymin)/pixwidth)
        ipixmax = int((x(i) + radkern - xmin)/pixwidth) + 1
        jpixmax = int((y(i) + radkern - ymin)/pixwidth) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidth
           dy = ypix - y(i)
           do ipix = ipixmin,ipixmax
              xpix = xmin + (ipix-0.5)*pixwidth
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
                 datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab      

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo                     ! over particles
  !
  !--normalise dat array(s)
  !
  where (datnorm > 0.)
     vecsmoothx = vecsmoothx/datnorm
     vecsmoothy = vecsmoothy/datnorm
  end where

  return

end subroutine interpolate3D_xsec_vec

end module xsections3D
