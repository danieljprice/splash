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
!            weight for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!--------------------------------------------------------------------------

subroutine interpolate3D(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth,zpixwidth,normalise)
  implicit none
  integer, intent(in) :: npart,npixx,npixy,npixz
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,zmin,pixwidth,zpixwidth
  real, intent(out), dimension(npixx,npixy,npixz) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy,npixz) :: datnorm

  integer :: i,ipix,jpix,kpix
  integer :: iprintinterval,iprintnext,iprogress
  integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
  real :: xminpix,yminpix,zminpix
  real, dimension(npixx) :: xpix,dx2i
  real :: xi,yi,zi,hi,hi1,hi21,radkern,qq,wab,q2,const,dyz2,dz2
  real :: term,termnorm,dx,dy,dz,ypix,zpix
  real :: t_start,t_end
  logical :: iprintprogress
  
  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print "(1x,a)",'interpolating from particles to 3D grid (normalised) ...'  
  else
     print "(1x,a)",'interpolating from particles to 3D grid (non-normalised) ...'
  endif
  if (pixwidth.le.0.) then
     print "(1x,a)",'interpolate3D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
  endif

  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval
  !
  !--get starting CPU time
  !
  call cpu_time(t_start)

  xminpix = xmin - 0.5*pixwidth
  yminpix = ymin - 0.5*pixwidth
  zminpix = zmin - 0.5*zpixwidth
!  xmax = xmin + npixx*pixwidth
!  ymax = ymin + npixy*pixwidth
!
!--store x value for each pixel (for optimisation)
!  
  do ipix=1,npixx
     xpix(ipix) = xminpix + ipix*pixwidth
  enddo
  const = dpi  ! normalisation constant (3D)
  !
  !--loop over particles
  !      
  over_parts: do i=1,npart
     !
     !--report on progress
     !
     if (mod(i,10000).eq.0) then
        call cpu_time(t_end)
        print*,i,t_end-t_start
     endif
     if (iprintprogress) then
        iprogress = 100*i/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts

     hi = hh(i)
     if (hi.le.0.) cycle over_parts

     !
     !--set kernel related quantities
     !
     xi = x(i)
     yi = y(i)
     zi = z(i)

     hi1 = 1./hi
     hi21 = hi1*hi1
     radkern = 2.*hi   ! radius of the smoothing kernel
     termnorm = const*weight(i)
     term = termnorm*dat(i)


     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((xi - radkern - xmin)/pixwidth)
     jpixmin = int((yi - radkern - ymin)/pixwidth)
     kpixmin = int((zi - radkern - zmin)/zpixwidth)
     ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
     jpixmax = int((yi + radkern - ymin)/pixwidth) + 1
     kpixmax = int((zi + radkern - zmin)/zpixwidth) + 1

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (kpixmin.lt.1) kpixmin = 1
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     if (kpixmax.gt.npixz) kpixmax = npixz
     !
     !--precalculate an array of dx2 for this particle (optimisation)
     !
     do ipix=ipixmin,ipixmax
        dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
     enddo
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do kpix = kpixmin,kpixmax
        zpix = zminpix + kpix*zpixwidth
        dz = zpix - zi
        dz2 = dz*dz*hi21
        
        do jpix = jpixmin,jpixmax
           ypix = yminpix + jpix*pixwidth
           dy = ypix - yi
           dyz2 = dy*dy*hi21 + dz2
           
           do ipix = ipixmin,ipixmax
              q2 = dx2i(ipix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.4.0) then                  
                 if (qq.lt.1.0) then
                    qq = sqrt(q2)
                    wab = 1.-1.5*q2 + 0.75*q2*qq
!
!-- the following lines use a fast inverse sqrt function
!
!                    if (q2.gt.epsilon(q2)) then
!                       qq = q2*finvsqrt(q2)
!                       wab = 1.-1.5*q2 + 0.75*q2*qq
!                    else
!                       wab = 1.
!                    endif
                 else
                    qq = sqrt(q2)
!                    qq = q2*finvsqrt(q2)
                    wab = 0.25*(2.-qq)**3
                 endif
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !  
                 datsmooth(ipix,jpix,kpix) = datsmooth(ipix,jpix,kpix) + term*wab          
                 if (normalise) datnorm(ipix,jpix,kpix) = datnorm(ipix,jpix,kpix) + termnorm*wab          
              endif
           enddo
        enddo
     enddo
  enddo over_parts
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > tiny(datnorm))
        datsmooth = datsmooth/datnorm
     end where
  endif
  
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

subroutine interpolate3D_fastxsec(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidth,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidth,zslice
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,qq2,wab,const,xi,yi,hi21
  real :: termnorm,term,dy,dy2,dz,dz2,ypix,rescalefac
  real, dimension(npixx) :: dx2i

  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised)...',zslice
  else
     print*,'taking fast cross section (non-normalised)...',zslice
  endif
  if (pixwidth.le.0.) then
     print*,'interpolate3D_xsec: error: pixel width <= 0'
     return
  elseif (npart.le.0) then
     print*,'interpolate3D_xsec: error: npart = 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec: WARNING: ignoring some or all particles with h < 0'
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
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
     hi1 = 1./hi
     hi21 = hi1*hi1
     radkern = 2.*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2*hi21
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (dz2 .lt. 4.0) then

        xi = x(i)
        yi = y(i)
        termnorm = const*weight(i)
        term = termnorm*dat(i)/rescalefac
        !
        !--for each particle work out which pixels it contributes to
        !               
        ipixmin = int((xi - radkern - xmin)/pixwidth)
        jpixmin = int((yi - radkern - ymin)/pixwidth)
        ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
        jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xmin + (ipix-0.5)*pixwidth - xi)**2)*hi21 + dz2
        enddo
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidth
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixmin,ipixmax
              qq2 = dx2i(ipix) + dy2
              !
              !--SPH kernel - standard cubic spline
              !     
              if (qq2.lt.4.0) then
                 qq = sqrt(qq2)
                 if (qq.lt.1.0) then
                    wab = (1.-1.5*qq2 + 0.75*qq*qq2)
                 else
                    wab = 0.25*(2.-qq)**3
                 endif
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !  
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                 if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo over_parts                    ! over particles
  !
  !--normalise dat array
  !
  if (normalise) then
     !--normalise everywhere (required if not using SPH weighting)
     where (datnorm > tiny(datnorm))
        datsmooth = datsmooth/datnorm
     end where
  endif
  datsmooth = datsmooth*rescalefac
  
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

subroutine interpolate3D_xsec_vec(x,y,z,hh,weight,vecx,vecy,itype,npart,&
     xmin,ymin,zslice,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidth,zslice
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: termx,termy,termnorm,dx,dy,dz,dz2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised)...',zslice
  else
     print*,'taking fast cross section (non-normalised)...',zslice
  endif
  if (pixwidth.le.0.) then
     print*,'interpolate3D_xsec_vec: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec_vec: WARNING: ignoring some or all particles with h < 0'
  endif
  const = dpi ! normalisation constant (3D)
  !
  !--loop over particles
  !      
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
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
                 if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab      

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo over_parts                    ! over particles
  !
  !--normalise dat array(s)
  !
  if (normalise) then
     where (datnorm > tiny(datnorm))
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif

  return

end subroutine interpolate3D_xsec_vec

end module xsections3D
