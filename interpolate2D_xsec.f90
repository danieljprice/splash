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
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data 	   : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Feb 2004
!--------------------------------------------------------------------------

subroutine interpolate2D_xsec(x,y,pmass,rho,hh,dat,npart,&
     x1,y1,x2,y2,datsmooth,npixx)

  implicit none
  real, parameter :: pi = 3.1415926536      
  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,y,pmass,rho,hh,dat
  real, intent(in) :: x1,y1,x2,y2
  real, intent(out), dimension(npixx) :: datsmooth

  integer :: i,ipix,ipixmin,ipixmax
  real :: hi,hi1,h2,radkern,qq,wab,rab,const
  real :: term,dx,dy,xpix,ypix,pixwidth,xpixwidth,xlength
  real :: gradient,yintercept,aa,bb,cc,determinant,det
  real :: xstart,xend,ystart,yend,rstart,rend
  real :: tol
  logical :: xsame, ysame
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

  !
  !--now interpolate to the line of pixels
  !
  datsmooth = 0.
  term = 0.
  !
  !--loop over particles
  !      
  do i=1,npart
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate2D_xsec: error: h <= 0 ',i,hi
	return
     endif
     hi1 = 1./hi
     h2 = hi*hi
     radkern = 2.*hi    ! radius of the smoothing kernel
     const = 10./(7.*pi*h2)   ! normalisation constant
     if (rho(i).ne.0.) term = pmass(i)*dat(i)/rho(i) 
     !
     !--for each particle work out which pixels it contributes to
     !  to do this we need to work out the two points at which the line
     !  intersects the particles smoothing circle
     !  given by the equation (x-xi)^2 + (y-yi)^2 = (2h)^2.             
     !  The x co-ordinates of these points are the solutions to a
     !  quadratic with co-efficients:

     aa = 1. + gradient**2
     bb = 2.*gradient*(yintercept - y(i)) - 2.*x(i)
     cc = x(i)**2 + y(i)**2 - 2.*yintercept*y(i) * yintercept**2 &
          - radkern**2
     !
     !--work out whether there are any real solutions and find them
     !
     determinant = bb**2 - 4.*aa*cc
     if (determinant < 0) then
        !    print*,' particle ',i,': does not contribute ',x(i),y(i) 
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
        ipixmax = int(rend/pixwidth)

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (ipixmax.lt.1) ipixmax = 1
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (ipixmin.gt.npixx) ipixmax = npixx
       	!
        !--loop over pixels, adding the contribution from this particle
        !
        do ipix = ipixmin,ipixmax

           xpix = x1 + (ipix)*xpixwidth - 0.5*xpixwidth
           ypix = gradient*xpix + yintercept
           dy = ypix - y(i)
           dx = xpix - x(i)
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - standard cubic spline in 2D
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
           datsmooth(ipix) = datsmooth(ipix) + term*wab          

        enddo

     endif

  enddo

  return

end subroutine interpolate2D_xsec
