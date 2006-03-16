module opacityrendering3D
 use projections3D, only:interpolate3D_projection,wfromtable,radkernel2,coltable
 implicit none
 private :: indexx

contains
!--------------------------------------------------------------------------
!     subroutine to do a ray trace through the particle data
!
!     we use the radiation transport equation along a ray, that is
!     the change in intensity from one side of a particle to the other is
!     given by:
!
!     I_nu = I_nu(0) exp(-tau_i) + S_nu (1 - exp(-tau_i))
!
!     where tau_i is the integrated optical depth through the particle,
!     and S_nu is the colour calculated from a colour table for the rendered data.
!     We calculate an intensity in red, green and blue for colour plots.
!
!     tau_i = kappa \int rho dz
!
!     this is calculated using the SPH kernel for rho, so for each pixel
!     the optical depth is incremented as the sum
!
!     tau_i = kappa \sum_j m_j \int W dz
!
!     where \int W dz is the SPH kernel integrated along one spatial dimension.
!     This is interpolated from a pre-calculated table (see module projections3D for this).
!
!     kappa is the monochromatic mass extinction coefficient 
!     (particle cross section per unit mass) and is a constant for all particles
!     which must be given as input (although see below for calculations of a
!     meaningful values for kappa in terms of "surface depth in units of smoothing lengths")
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            particle masses       : pmass (npart)
!            smoothing lengths     : hh    (npart)
!            scalar data to smooth : dat   (npart)
!
!     Settings: zobs, dz1 : settings for 3D projection
!               rkappa    : particle cross section per unit mass
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     NB: AT PRESENT WE WRITE A PPM FILE DIRECTLY WITH THE RGB COLOURS
!         AND OUTPUT JUST THE MONOCHROMATIC VERSION TO SUPERSPHPLOT.
!
!     (c) 2005 Daniel Price. Last modified Dec 2005.
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_opacity(x,y,z,pmass,hh,dat,zorig,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth,zobserver,dscreenfromobserver, &
     rkappa,zcut,datmin,datmax,itrans,istep)

  use transforms
  use colours, only:rgbtable,ncolours
  implicit none
  real, parameter :: pi=3.1415926536
  integer, intent(in) :: npart,npixx,npixy,itrans,istep
  real, intent(in), dimension(npart) :: x,y,z,pmass,hh,dat,zorig
  real, intent(in) :: xmin,ymin,pixwidth,zobserver,dscreenfromobserver, &
                      zcut,datmin,datmax,rkappa
  real, dimension(npixx,npixy), intent(out) :: datsmooth
  real, dimension(npixx,npixy) :: brightness
!  real, dimension(3,npixx,npixy) :: rgb
  real, dimension(3) :: rgbi,drgb

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,nused
  integer :: iprintinterval, iprintnext, iprogress, itmin
  integer, dimension(npart) :: iorder
  integer :: ipart,ir,ib,ig,ierr,maxcolour,indexi
  real :: hi,hi1,hi21,radkern,q2,wab,rab2,pmassav
  real :: term,dx,dy,dy2,xpix,ypix,zfrac,hav
  real :: fopacity,tau,rkappatemp,termi,xi,yi
  real, dimension(1) :: dati
  real :: t_start,t_end,t_used,tsec
  real :: ddatrange,datfraci,ftable
  logical :: iprintprogress
  character(len=120) :: filename

  datsmooth = 0.
  term = 0.
  brightness = 0.
  print "(1x,a)",'ray tracing from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_opacity: error: pixel width <= 0'
     return
  endif
  if (abs(datmax-datmin).gt.tiny(datmin)) then
     ddatrange = 1./abs(datmax-datmin)
  else
     print "(a)",'error: datmin=datmax in opacity rendering'
     return
  endif

!
!--kappa is the opacity in units of length^2/mass
!  sent as an input parameter as it should be kept constant throughout the simulation
!
!  However we compute a reasonable estimate below based on the current plot so that
!  we can give the "actual" optical depth for the current frame in terms of number of
!  smoothing lengths. This is purely for diagnostic purposes only.
!
!--calculate average h
  hav = sum(hh(1:npart))/real(npart)
!--average particle mass
  pmassav = sum(pmass(1:npart))/real(npart)
  rkappatemp = pi*hav*hav/(pmassav*coltable(0))
  print*,'average h = ',hav,' average mass = ',pmassav
  print "(1x,a,f6.2,a)",'typical surface optical depth is ~',rkappatemp/rkappa,' smoothing lengths'  
!  rgb = 0.
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
!--first sort the particles in z so that we do the opacity in the correct order
!
  call indexx(npart,z,iorder)
  
  nused = 0
  
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
     !--allow slicing [take only particles with z(unrotated) < zcut]
     !
     particle_within_zcut: if (zorig(i).lt.zcut) then
    
     !  count particles within slice
     nused = nused + 1
     !
     !--adjust h according to 3D perspective
     !  need to be careful -- the kernel quantities
     !  change with z (e.g. radkern, r^2/h^2)
     !  but *not* the 1/h^2 in tau (because the change in 1/h^2 in tau
     !  would be cancelled by the corresponding change to h^2 in kappa)
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D_proj_opacity: error: h <= 0 ',i,hi
        return
     elseif (abs(dscreenfromobserver).gt.tiny(dscreenfromobserver)) then
        zfrac = abs(dscreenfromobserver/(z(i)-zobserver))
        hi = hi*zfrac
     endif
     
     !--these are the quantities used in the kernel r^2/h^2
     radkern = 2.*hi
     hi1 = 1./hi
     hi21 = hi1*hi1
     !--this is the term which multiplies tau
     term = pmass(i)/(hh(i)*hh(i))
     !
     !--determine colour contribution of current point
     !  (work out position in colour table)
     !
!     dati = dat(i)
     xi = x(i)
     yi = y(i)
     termi = dat(i)
!     call transform(dati,itrans)
!     datfraci = (dati(1) - datmin)*ddatrange
!     datfraci = max(datfraci,0.)
!     datfraci = min(datfraci,1.)
!     !--define colour for current particle
!     ftable = datfraci*ncolours
!     indexi = int(ftable) + 1
!     indexi = min(indexi,ncolours)
!     if (indexi.lt.ncolours) then
!     !--do linear interpolation from colour table
!        drgb(:) = rgbtable(:,indexi+1) - rgbtable(:,indexi)
!        rgbi(:) = rgbtable(:,indexi) + (ftable - int(ftable))*drgb(:)
!     else
!        rgbi(:) = rgbtable(:,indexi)
!     endif
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((xi - radkern - xmin)/pixwidth)
     jpixmin = int((yi - radkern - ymin)/pixwidth)
     ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
     jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
     if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix-0.5)*pixwidth
        dy = ypix - yi
        dy2 = dy*dy
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix-0.5)*pixwidth
           dx = xpix - xi
           rab2 = dx**2 + dy2
           q2 = rab2*hi21
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--get incremental tau for this pixel from the integrated SPH kernel
              !
              tau = rkappa*wab*term
              fopacity = 1. - exp(-tau)
              !
              !--render, obscuring previously drawn pixels by relevant amount
              !  also calculate total brightness (`transparency') of each pixel
              !
              datsmooth(ipix,jpix) = (1.-fopacity)*datsmooth(ipix,jpix) + fopacity*termi
              brightness(ipix,jpix) = brightness(ipix,jpix) + fopacity
           endif

        enddo
     enddo
     
     endif particle_within_zcut

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
  maxcolour = 255
  write(78,"(a)") 'P3'
  write(78,"(a)") '# supersphplot.ppm created by supersphplot (c) 2005 Daniel Price'
  write(78,"(i4,1x,i4)") npixx, npixy
  write(78,"(i3)") maxcolour
!--pixel information
  do jpix = npixy,1,-1
     do ipix = 1,npixx

        dati(1) = datsmooth(ipix,jpix)
        call transform(dati,itrans)
        datfraci = (dati(1) - datmin)*ddatrange
        datfraci = max(datfraci,0.)
        datfraci = min(datfraci,1.)
        !--define colour for current particle
        ftable = datfraci*ncolours
        indexi = int(ftable) + 1
        indexi = min(indexi,ncolours)
        if (indexi.lt.ncolours) then
        !--do linear interpolation from colour table
           drgb(:) = rgbtable(:,indexi+1) - rgbtable(:,indexi)
           rgbi(:) = rgbtable(:,indexi) + (ftable - int(ftable))*drgb(:)
        else
           rgbi(:) = rgbtable(:,indexi)
        endif
        rgbi(:) = rgbi(:)*min(brightness(ipix,jpix),1.0)
        ir = max(min(int(rgbi(1)*maxcolour),maxcolour),0)
        ig = max(min(int(rgbi(2)*maxcolour),maxcolour),0)
        ib = max(min(int(rgbi(3)*maxcolour),maxcolour),0)

!        ir = max(min(int(rgb(1,ipix,jpix)*maxcolour),maxcolour),0)
!        ig = max(min(int(rgb(2,ipix,jpix)*maxcolour),maxcolour),0)
!        ib = max(min(int(rgb(3,ipix,jpix)*maxcolour),maxcolour),0)
!!        if (rgb(1,ipix,jpix).gt.0.999) print*,rgb(1,ipix,jpix),ir
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
  if (zcut.lt.huge(zcut)) print*,'slice contains ',nused,' of ',npart,' particles'
  
  return

end subroutine interpolate3D_proj_opacity


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
      print*,' this is in the file interpolate3D_opacity.f90)'
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

end module opacityrendering3D
