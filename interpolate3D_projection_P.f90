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
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!     OpenMP parallelized by DJP Nov 2004
!--------------------------------------------------------------------------

subroutine interpolate3D_projection(x,y,pmass,rho,hh,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth)

  use column
  implicit none
  real, parameter :: pi = 3.1415926536      
  integer, intent(IN) :: npart,npixx,npixy
  real, intent(IN), dimension(npart) :: x,y,pmass,rho,hh,dat
  real, intent(IN) :: xmin,ymin,pixwidth
  real, intent(OUT), dimension(npixx,npixy) :: datsmooth

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: index, index1
  integer :: iprintnext, iprogress
  real :: hi,hi1,h2,radkern,qq,wab,rab,const
  real :: term,dx,dy,xpix,ypix
  real :: dxx,dwdx
  real :: dmaxcoltable
  logical :: iprintprogress

  datsmooth = 0.
  term = 0.
  dmaxcoltable = 1./real(maxcoltable)
  print*,'projecting from particles to pixels...'
  if (pixwidth.le.0.) then
     print*,'interpolate3D_proj: error: pixel width <= 0'
     return
  endif
  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .and. (npixx*npixy .gt.11000)
  !
  !--loop over particles
  !
  iprintnext = 25
!$OMP PARALLEL default(none)
!$OMP& shared(hh,x,pmass,dat,rho,datsmooth,npart)
!$OMP& shared(xmin,xmax,ymin,ymax)
!$OMP& shared(npixx,npixy,pixwidth)
!$OMP& shared(coltable,maxcoltable)
!$OMP& private(hi,hi1,h2,radkern,const,term)
!$OMP& private(ipixmin,ipixmax,jpixmin,jpixmax)
!$OMP& private(dx,dy,rab,qq,index,index1,dxx,wab,dwdx)
!$OMP& private(i,ipix,jpix)
!$OMP DO SCHEDULE (runtime)  
  do i=1,npart
     !
     !--report on progress
     !
     !if (iprintprogress) then
     !   iprogress = 100*i/npart
     !   if (iprogress.ge.iprintnext) then
     !      write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
     !      iprintnext = iprintnext + 25
     !   endif
     !endif
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     !if (hi.le.0.) then
     !  print*,'interpolate3D_proj: error: h <= 0 ',i,hi
     !   return
     !endif
     hi1 = 1./hi
     h2 = hi*hi
     radkern = 2.*hi  !radius of the smoothing kernel
     !         const = 10./(7.*pi*h2)  ! normalisation constant
     const = hi1*hi1
     if (rho(i).ne.0.) term = const*pmass(i)*dat(i)/rho(i) 
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)
     jpixmax = int((y(i) + radkern - ymin)/pixwidth)

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
     if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        do ipix = ipixmin,ipixmax

           ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
           xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
           dy = ypix - y(i)
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
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab          
           endif

        enddo
     enddo

  enddo
!$OMP END DO
!$OMP END PARALLEL

  return

end subroutine interpolate3D_projection
