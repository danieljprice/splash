!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 2D interpolation
!
!----------------------------------------------------------------------

module interpolations2D
 implicit none
 real, parameter, private :: pi = 3.1415926536      
 public :: interpolate2D, interpolate2D_xsec, interpolate2D_vec
 
contains
!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_j w_j dat_j W(r-r_j, h_j)
! 
!     where _j is the quantity at the neighbouring particle j and
!     W is the smoothing kernel, for which we use the usual cubic spline.
!     For an SPH interpolation the weight for each particle should be
!     the dimensionless quantity
!
!     w_j = m_j / (rho_j * h_j**ndim)
!
!     Other weights can be used (e.g. constants), but in this case the
!     normalisation option should also be set.
!
!     Input: particle coordinates  : x,y    (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat    (npart)
!
!            number of pixels in x,y : npixx,npixy
!            pixel width             : pixwidth
!            option to normalise interpolation : normalise (.true. or .false.)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price 2003-2006
!--------------------------------------------------------------------------

subroutine interpolate2D(x,y,hh,weight,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  real, intent(in) :: xmin,ymin,pixwidth
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,termnorm,dx,dy,xpix,ypix

  datsmooth = 0.
  datnorm = 0.
  print*,'interpolating from particles to 2D grid...'
  if (pixwidth.le.0.) then
     print*,'interpolate2D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate2D: warning: ignoring some or all particles with h < 0'
  endif
  const = 10./(7.*pi)  ! normalisation constant
  !
  !--loop over particles
  !      
  over_parts: do i=1,npart
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = 2.*hi  ! radius of the smoothing kernel
     term = termnorm*dat(i)
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
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - standard cubic spline
           !                     
           if (qq.lt.1.0) then
              wab = (1.-1.5*qq**2 + 0.75*qq**3)
           elseif (qq.lt.2.0) then
              wab = 0.25*(2.-qq)**3
           else
              wab = 0.
           endif
           !
           !--calculate data value at this pixel using the summation interpolant
           !  
           datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
           if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab

        enddo
     enddo

  enddo over_parts
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > 0.)
        datsmooth = datsmooth/datnorm
     end where
  endif

  return

end subroutine interpolate2D

!--------------------------------------------------------------------------
!     
!     ** this version does vector quantities
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field    : vecsmoothx (npixx,npixy)
!                                      : vecsmoothy (npixx,npixy)
!
!     Daniel Price, University of Exeter, March 2005
!--------------------------------------------------------------------------

subroutine interpolate2D_vec(x,y,hh,weight,vecx,vecy,npart, &
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx,vecsmoothy
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: termnorm,termx,termy,dx,dy,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  print*,'interpolating vector field from particles to 2D grid...'
  if (pixwidth.le.0.) then
     print*,'interpolate2D_vec: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate2D_vec: warning: ignoring some or all particles with h < 0'
  endif
  const = 10./(7.*pi)  ! normalisation constant
  !
  !--loop over particles
  !      
  over_parts: do i=1,npart
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = 2.*hi  ! radius of the smoothing kernel
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
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - standard cubic spline
           !                     
           if (qq.lt.1.0) then
              wab = (1.-1.5*qq**2 + 0.75*qq**3)
           elseif (qq.lt.2.0) then
              wab = 0.25*(2.-qq)**3
           else
              wab = 0.
           endif
           !
           !--calculate data value at this pixel using the summation interpolant
           !  
           vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab
           vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab
           if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab

        enddo
     enddo

  enddo over_parts
  !
  !--normalise dat arrays
  !
  if (normalise) then
     where (datnorm > 0.)
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif

  return

end subroutine interpolate2D_vec

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!     
!     this version takes any 1D cross section through a 2D data set
!     the 1D line is specified by two points, (x1,y1) and (x2,y2)
!     (ie. this is for arbitrary oblique cross sections)
!
!     NB: A similar version could be used for 2D oblique cross sections
!         of 3D data. In this case we would need to find the intersection
!         between the smoothing sphere and the cross section plane. However
!         in 3D it is simpler just to rotate the particles first and then take
!         a straight cross section.
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Feb 2004
!--------------------------------------------------------------------------

subroutine interpolate2D_xsec(x,y,hh,weight,dat,npart,&
     x1,y1,x2,y2,datsmooth,npixx,normalise)

  implicit none
  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  real, intent(in) :: x1,y1,x2,y2
  real, intent(out), dimension(npixx) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx) :: datnorm

  integer :: i,ipix,ipixmin,ipixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,termnorm,dx,dy,xpix,ypix,pixwidth,xpixwidth,xlength
  real :: gradient,yintercept,aa,bb,cc,determinant,det
  real :: xstart,xend,ystart,yend,rstart,rend
  real :: tol
  logical :: xsame, ysame, debug
  
  debug = .false.
  !
  !--check for errors in input
  !
  tol = 1.e-3
  ysame = (abs(y2 - y1).lt.tol)
  xsame = (abs(x2 - x1).lt.tol)
  if (xsame.and.ysame) then
     print*,'error: interpolate2D_xsec: zero length cross section'
     return
  endif
  if (npixx.eq.0) then
     print*,'error: interpolate2D_xsec: npix = 0 '
     return
  endif
  print*,'oblique 1D cross section through 2D data: npix =',npixx
  !
  !--work out the equation of the line y = mx + c from the two points input
  !
  gradient = 0.
  if (.not.xsame) gradient = (y2-y1)/(x2-x1)
  yintercept = y2 - gradient*x2
  print*,'cross section line: y = ',gradient,'x + ',yintercept
  !
  !--work out length of line and divide into pixels
  !
  xlength = sqrt((x2-x1)**2 + (y2-y1)**2)
  pixwidth = xlength/real(npixx)
  xpixwidth = (x2 - x1)/real(npixx)
  if (debug) then
     print*,'length of line = ',xlength
     print*,'pixel width = ',pixwidth, ' in x direction = ',xpixwidth
  endif
  !
  !--now interpolate to the line of pixels
  !
  datsmooth = 0.
  datnorm = 0.
  const = 10./(7.*pi)   ! normalisation constant
  !
  !--loop over particles
  !      
  over_parts: do i=1,npart
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = 2.*hi    ! radius of the smoothing kernel
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !  to do this we need to work out the two points at which the line
     !  intersects the particles smoothing circle
     !  given by the equation (x-xi)^2 + (y-yi)^2 = (2h)^2.             
     !  The x co-ordinates of these points are the solutions to a
     !  quadratic with co-efficients:

     aa = 1. + gradient**2
     bb = 2.*gradient*(yintercept - y(i)) - 2.*x(i)
     cc = x(i)**2 + y(i)**2 - 2.*yintercept*y(i) + yintercept**2 &
          - radkern**2
     !
     !--work out whether there are any real solutions and find them
     !
     determinant = bb**2 - 4.*aa*cc
     if (determinant < 0) then
        !!print*,' particle ',i,': does not contribute ',x(i),y(i) 
     else
        det = sqrt(determinant)
        xstart = (-bb - det)/(2.*aa)
        xend =  (-bb + det)/(2.*aa)
        if (xstart.lt.x1) xstart = x1
        if (xstart.gt.x2) xstart = x2
        if (xend.gt.x2) xend = x2
        if (xend.lt.x1) xend = x1
        ystart = gradient*xstart + yintercept
        yend = gradient*xend + yintercept
        !
        !--work out position in terms of distance (no. of pixels) along the line
        !
        rstart = sqrt((xstart-x1)**2 + (ystart-y1)**2)
        rend = sqrt((xend-x1)**2 + (yend-y1)**2)

        ipixmin = int(rstart/pixwidth)
        ipixmax = int(rend/pixwidth) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (ipixmax.lt.1) ipixmax = 1
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (ipixmin.gt.npixx) ipixmax = npixx
               !
        !--loop over pixels, adding the contribution from this particle
        !
        !if (debug) print*,' particle ',i,': ',ipixmin,ipixmax,xstart,x(i),xend 
        do ipix = ipixmin,ipixmax

           xpix = x1 + (ipix-0.5)*xpixwidth
           ypix = gradient*xpix + yintercept
           dy = ypix - y(i)
           dx = xpix - x(i)
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - standard cubic spline in 2D
           !     
           if (qq.lt.1.0) then
              wab = (1.-1.5*qq**2 + 0.75*qq**3)
           elseif (qq.lt.2.0) then
              wab = 0.25*(2.-qq)**3
           else
              wab = 0.
           endif
           !
           !--calculate data value at this pixel using the summation interpolant
           !  
           datsmooth(ipix) = datsmooth(ipix) + term*wab
           if (normalise) datnorm(ipix) = datnorm(ipix) + termnorm*wab      

        enddo

     endif

  enddo over_parts
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > 0.)
        datsmooth = datsmooth/datnorm
     end where
  endif

  return

end subroutine interpolate2D_xsec

end module interpolations2D
