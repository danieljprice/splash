!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 3D projections
!  (where rendered quantity is integrated along the line of sight)
!
!----------------------------------------------------------------------

module projections3D
 implicit none

 integer, parameter, private :: maxcoltable = 1000
 real, parameter, private :: pi = 3.1415926536
 real, parameter, private :: radkernel = 2.0
 real, dimension(maxcoltable), private :: coltable
 real, parameter :: dmaxcoltable = 1./real(maxcoltable)

 public :: setup_integratedkernel
 public :: interpolate3D_projection
 public :: interpolate3D_proj_vec
 public :: interpolate3D_proj_opacity
 private :: indexx

contains

subroutine setup_integratedkernel
!-------------------------------------------------------------
!     tabulates the integral through the cubic spline kernel
!     subroutine originally by Matthew Bate
!-------------------------------------------------------------
 implicit none
 integer :: i,j
 real :: r, dist, step, ypos, v, v2, val
 real :: coldens, v2m

 print "(1x,a)",'setting up integrated kernel table...'

 do i=1,maxcoltable
    r=(i-1)/500.
    dist=sqrt(4.0-r*r)
    step=dist/4000.0
    ypos=0.
         
    coldens=0.0
    do j=1,4000
       v=sqrt(r*r+ypos*ypos)
       if (v.lt.1.0) then
          v2=v*v
          val=1.0-1.5*v2+0.75*v2*v
          coldens=coldens+val*step
       else
          v2m=2.0-v
          val=0.25*v2m*v2m*v2m
          coldens=coldens+val*step        
       endif
       ypos=ypos+step
    end do
    coltable(i)=2.0*coldens/pi
 end do

!      check=0.0
!      do i=1,1000
!         r=(i-1)/500.
!         write(*,*) r, table(i)
!         check=check+2*pi*r*table(i)*(1.0/500.)
!      end do
!      write(*,*)check

 return
end subroutine setup_integratedkernel

!
! This function interpolates from the table of integrated kernel values
! to give w(q)
!
real function wfromtable(qq)
 implicit none
 real :: qq,dxx,dwdx
 integer :: index, index1
 !
 !--find nearest index in table
 !
 index = int(0.5*qq*maxcoltable) + 1
 index1 = min(index + 1,maxcoltable)
 !
 !--find increment along from this index
 !
 dxx = 0.5*qq*maxcoltable - index*dmaxcoltable
 !
 !--find gradient
 !
 dwdx = (coltable(index1) - coltable(index))*dmaxcoltable
 !
 !--compute value of integrated kernel
 !
 wfromtable = coltable(index) + dwdx*dxx

end function wfromtable

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
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!     3D perspective added Nov 2005
!--------------------------------------------------------------------------

subroutine interpolate3D_projection(x,y,z,pmass,rho,hh,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth,zobs,dz1)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,pmass,rho,hh,dat
  real, intent(in) :: xmin,ymin,pixwidth,zobs,dz1
  real, intent(out), dimension(npixx,npixy) :: datsmooth

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: iprintinterval, iprintnext, iprogress, itmin
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,rho1i,dx,dy,xpix,ypix,zfrac
  real :: t_start,t_end,t_used,tsec
  logical :: iprintprogress

  datsmooth = 0.
  term = 0.
  print "(1x,a)",'projecting from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj: error: pixel width <= 0'
     return
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
  
  over_particles: do i=1,npart
     !
     !--report on progress
     !
     if (iprintprogress) then
        iprogress = 100*i/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D_proj: error: h <= 0 ',i,hi
        return
     elseif (abs(dz1).gt.tiny(dz1)) then
        zfrac = abs(dz1/(z(i)-zobs))
        hi = hi*zfrac
     endif
     hi1 = 1./hi
     radkern = 2.*hi  !radius of the smoothing kernel
     !         const = 10./(7.*pi*h2)  ! normalisation constant
     const = hi1*hi1
     if (rho(i).gt.0.) then
        rho1i = 1./rho(i)
     else
        rho1i = 0.
     endif
     term = const*pmass(i)*dat(i)*rho1i
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
           if (qq.lt.radkernel) then
              wab = wfromtable(qq)
              !
              !--calculate data value at this pixel using the summation interpolant
              !                  
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab          
           endif

        enddo
     enddo

  enddo over_particles
!
!--get ending CPU time
!
  call cpu_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.60) then
     itmin = int(t_used/60)
     tsec = t_used - (itmin*60)
     print*,'completed in ',itmin,' min ',tsec,'s'
  else
     print*,'completed in ',t_used,'s'
  endif
  
  return

end subroutine interpolate3D_projection

!--------------------------------------------------------------------------
!
!     Same as previous but for a vector quantity
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

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,pmass,rho,hh,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: rho1i,termx,termy,dx,dy,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  print "(1x,a)",'projecting vector from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj_vec: error: pixel width <= 0'
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
           if (qq.lt.radkernel) then
              wab = wfromtable(qq)
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
!     ** In this version the opacity is set according to the column density
!     ** and the colour corresponds to the rendered quantity. The effect is
!     ** to show the surface values of the quantity and to "see through"
!     ** low density regions
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price Nov 2005
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_opacity(x,y,z,pmass,rho,hh,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth,zobs,dz1,rhomin,rhomax,datmin,datmax,itrans,istep)

  use transforms
  implicit none
  integer, intent(in) :: npart,npixx,npixy,itrans,istep
  real, intent(in), dimension(npart) :: x,y,z,pmass,rho,hh,dat
  real, intent(in) :: xmin,ymin,pixwidth,zobs,dz1,rhomin,rhomax,datmin,datmax
  real, dimension(npixx,npixy), intent(out) :: datsmooth
  real, dimension(3,npixx,npixy) :: rgb

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: iprintinterval, iprintnext, iprogress, itmin
  integer, dimension(npart) :: iorder
  integer :: ipart,ifirstcolour,ilastcolour,ir,ib,ig,ierr,maxcolour,indexi
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,rho1i,dx,dy,xpix,ypix,zfrac
  real :: dwnorm
  real :: drhorange,fopacity,rhoi
  real, dimension(1) :: dati
  real :: t_start,t_end,t_used,tsec
  real :: rgbtable(3,256),ddatrange,datfraci
  logical :: iprintprogress
  character(len=120) :: filename

  datsmooth = 0.
  term = 0.
  dwnorm = 1./coltable(1)
  print "(1x,a)",'projecting (with variable opacity) from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj: error: pixel width <= 0'
     return
  endif
  if (abs(rhomax-rhomin).gt.tiny(rhomin)) then
     drhorange = 1./abs(rhomax-rhomin)
  else
     print "(a)",'error: rhomin=rhomax in opacity rendering'
     return
  endif
  if (abs(datmax-datmin).gt.tiny(datmin)) then
     ddatrange = 1./abs(datmax-datmin)
  else
     print "(a)",'error: datmin=datmax in opacity rendering'
     return
  endif
  rgb = 0.
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
!
!--get the rgb colours from the colour table
!  
  call tabulate_colours(rgbtable,ifirstcolour,ilastcolour)
!
!--first sort the particles in z so that we do the opacity in the correct order
!
  call indexx(npart,z,iorder)
  
  over_particles: do ipart=1,npart
     !
     !--report on progress
     !
     if (iprintprogress) then
        iprogress = 100*ipart/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,ipart
           iprintnext = iprintnext + iprintinterval
        endif
     endif
     !
     !--render in order from back to front
     !
     i = iorder(ipart)
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D_proj_opacity: error: h <= 0 ',i,hi
        return
     elseif (abs(dz1).gt.tiny(dz1)) then
        zfrac = abs(dz1/(z(i)-zobs))
        hi = hi*zfrac
     endif
     hi1 = 1./hi
     radkern = 2.*hi  !radius of the smoothing kernel
     const = hi1*hi1
     rhoi = rho(i)
!     if (rhoi.gt.0.) then
!        rho1i = 1./rhoi
!     else
!        rho1i = 0.
!     endif
     !term = const*pmass(i)*dat(i)*rho1i
     term = dat(i)
     dati = term
     call transform(dati,itrans)
     datfraci = (dati(1) - datmin)*ddatrange
     datfraci = max(datfraci,0.)
     datfraci = min(datfraci,1.)
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
           if (qq.lt.radkernel) then
              wab = wfromtable(qq)
              !--use normalised value here
              wab = wab*dwnorm
              !
              !--opacity is the fractional density
              !  NB this is not quite right if log is applied to density
              !
              fopacity = (rhoi - rhomin)*drhorange*wab
              fopacity = min(fopacity,1.0)
              fopacity = max(fopacity,0.0)
              !
              !--determine colour contribution of current point
              !  (work out position in colour table)
              !  NB: wab does not make sense here if we have taken the log of dat
              !
              indexi = int(datfraci*(ilastcolour-ifirstcolour)) + ifirstcolour
              !
              !--render, obscuring previously drawn pixels by relevant amount
              !
              rgb(1,ipix,jpix) = (1.-fopacity)*rgb(1,ipix,jpix) + fopacity*rgbtable(1,indexi)
              rgb(2,ipix,jpix) = (1.-fopacity)*rgb(2,ipix,jpix) + fopacity*rgbtable(2,indexi)
              rgb(3,ipix,jpix) = (1.-fopacity)*rgb(3,ipix,jpix) + fopacity*rgbtable(3,indexi)
              !
              !--this is the rendering of the colour value -- ie. position in the colour table
              !  previously drawn colours (data values) are obscured by relevant amount
              ! 
              datsmooth(ipix,jpix) = (1.-fopacity)*datsmooth(ipix,jpix) + fopacity*(term*wab)
           endif

        enddo
     enddo

  enddo over_particles
!
!--write PPM--
!  
  write(filename,"(a,i5.5,a)") 'supersphplot_',istep,'.ppm' 
  open(unit=78,file=filename,status='replace',form='formatted',iostat=ierr)
  if (ierr /=0) then
     print*,'error opening ppm file'
     return
  endif
  print "(1x,a,i5.5,a)", 'writing to file supersphplot_',istep,'.ppm' 
!
!--PPM header
!
  maxcolour = 256
  write(78,"(a)") 'P3'
  write(78,"(a)") '# supersphplot.ppm created by supersphplot (c) 2005 Daniel Price'
  write(78,"(i4,1x,i4)") npixx, npixy
  write(78,"(i3)") maxcolour
!--pixel information
  do jpix = npixy,1,-1
     do ipix = 1,npixx
        ir = max(min(int(rgb(1,ipix,jpix)*maxcolour),maxcolour),0)
        ig = max(min(int(rgb(2,ipix,jpix)*maxcolour),maxcolour),0)
        ib = max(min(int(rgb(3,ipix,jpix)*maxcolour),maxcolour),0)
        write(78,"(i3,1x,i3,1x,i3,2x)") ir,ig,ib
     enddo
  enddo
  close(unit=78)
! 
!--get ending CPU time
!
  call cpu_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.60) then
     itmin = int(t_used/60)
     tsec = t_used - (itmin*60)
     print*,'completed in ',itmin,' min ',tsec,'s'
  else
     print*,'completed in ',t_used,'s'
  endif
  
  return

end subroutine interpolate3D_proj_opacity

!
!--this subroutine queries PGPLOT for the current colour table
!  returns colour table and start/end colour indices
!
subroutine tabulate_colours(rgbtable,istart,iend)
 implicit none
 real, intent(out), dimension(:,:) :: rgbtable
 integer, intent(out) :: istart,iend
 integer :: i
 
 call pgqcir(istart,iend)
 if (iend.gt.size(rgbtable(1,:))) then
    print*,'error: too many colours for array size in tabulate_colours'
 endif
 
 do i=istart,iend
    call pgqcr(i,rgbtable(1,i),rgbtable(2,i),rgbtable(3,i))
 enddo
 
 return
end subroutine tabulate_colours


subroutine indexx(n, arr, indx)
!************************************************************
!                                                           *
!  This is INDEXX using the quicksort algorithm.            *
!                                                           *
!************************************************************
 implicit none
 integer, parameter :: m=7, nstack=500
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: arr
 integer, dimension(n), intent(out) :: indx

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer, dimension(nstack) :: istack
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l.lt.m) then
   do j = l + 1, ir
      indxt = indx(j)
      a = arr(indxt)
      do i = j - 1, 1, -1
         if (arr(indx(i)).le.a) goto 2
         indx(i + 1) = indx(i)
      end do
      i = 0
2     indx(i + 1) = indxt
   end do
   if (jstack.eq.0) return
   ir = istack(jstack)
   l = istack(jstack - 1)
   jstack = jstack - 2
  else
   k = (l + ir)/2
   itemp = indx(k)
   indx(k) = indx(l + 1)
   indx(l + 1) = itemp
   if (arr(indx(l + 1)).gt.arr(indx(ir))) then
      itemp = indx(l + 1)
      indx(l + 1) = indx(ir)
      indx(ir) = itemp
   endif
   if (arr(indx(l)).gt.arr(indx(ir))) then
      itemp = indx(l)
      indx(l) = indx(ir)
      indx(ir) = itemp
   endif
   if (arr(indx(l + 1)).gt.arr(indx(l))) then
      itemp = indx(l + 1)
      indx(l + 1) = indx(l)
      indx(l) = itemp
   endif
   i = l + 1
   j = ir
   indxt = indx(l)
   a = arr(indxt)

3  continue
   i = i + 1
   if (arr(indx(i)).lt.a) goto 3
4  continue
   j = j - 1
   if (arr(indx(j)).gt.a) goto 4
   if (j.lt.i) goto 5
   itemp = indx(i)
   indx(i) = indx(j)
   indx(j) = itemp
   goto 3

5  indx(l) = indx(j)
   indx(j) = indxt
   jstack = jstack + 2
   if (jstack.gt.nstack) then
      print*,'fatal error!!! stacksize exceeded in sort'
      print*,'(need to set parameter nstack higher in subroutine indexx '
      print*,' this is in the file interpolate3D_projection.f90)'
      stop
   endif
   if (ir - i + 1.ge.j - l) then
      istack(jstack) = ir
      istack(jstack - 1) = i
      ir = j - 1
   else
      istack(jstack) = j - 1
      istack(jstack - 1) = l
      l = i
   endif
endif

goto 1
end subroutine indexx

end module projections3D
