!-----------------------------------------------------------------
!  module handling plotting of azimuthally-averaged quantities
!  for disc simulations
!-----------------------------------------------------------------
module disc
 implicit none
 integer, parameter, private :: nbins = 350
 integer, dimension(nbins), private :: ninbin
 real, dimension(nbins), private :: radius,sigma,spsound

 public :: disccalc,discplot

 private

contains

subroutine disccalc(iplot,npart,rpart,npmass,pmass,rminin,rmaxin,ymin,ymax,itransx,itransy,utherm)
 use transforms, only:transform_limits_inverse,transform_inverse,transform
 implicit none
 integer, intent(in) :: iplot,npart,npmass,itransx,itransy
 real, dimension(npart), intent(in) :: rpart
 real, dimension(npmass), intent(in) :: pmass
 real, dimension(npart), intent(in), optional :: utherm
 real, intent(in) :: rminin,rmaxin
 real, intent(out) :: ymin,ymax
 integer :: i,ibin
 real, parameter :: pi = 3.1415926536
 real :: pmassi,rbin,deltar,area,rmin,rmax
 real :: sigmai,toomreq,epicyclic,Omegai,spsoundi
 real, dimension(1) :: rad

 ninbin(:) = 0
 sigma(:) = 0.
 spsound(:) = 0.

 pmassi = 0
 if (npmass.le.0) then
    print*,' INTERNAL ERROR in discplot: dimension of mass array <= 0'
    return
 endif
!
!--print info
!
 select case(iplot)
 case(1)
    print "(a,i10,a,i4)",' calculating disc surface density profile: npart = ',npart, ' nbins = ',nbins
 case(2)
    if (present(utherm)) then
       print "(a)",' calculating Toomre Q parameter (assuming Mstar=1 and a Keplerian rotation profile)'
    else
       print "(a)",' ERROR: cannot calculate Toomre Q parameter: thermal energy not present in dump file'
       return
    endif
 case default
    print "(a)",' ERROR: unknown plot in discplot. '
 end select
!
!--if transformations (e.g. log) are applied to r, then limits
!  will already be set in transformed space - need to obtain
!  limits in non-transformed space.
!
 rmin = rminin
 rmax = rmaxin
 if (itransx.gt.0) call transform_limits_inverse(rmin,rmax,itransx)
!
!--set array of radius values for plotting
!
 deltar = (rmax - rmin)/(nbins - 1)
 do ibin=1,nbins
    radius(ibin) = rmin + (ibin-0.5)*deltar
 enddo
!
!--calculate surface density in each radial bin
!
!$omp parallel default(none) &
!$omp shared(npart,rpart,sigma,npmass,pmass,itransx,rmin,deltar) &
!$omp shared(ninbin,spsound,utherm) &
!$omp private(i,rad,pmassi,ibin,rbin,area)
!$omp do
 do i=1,npart
    if (itransx.eq.0) then
       rad(1) = rpart(i)
    else
       rad(1) = rpart(i)
       call transform_inverse(rad,itransx)
    endif
    if (npmass.ge.npart) then
       pmassi = pmass(i)
    else
       pmassi = pmass(1)
    endif
    ibin = int((rad(1) - rmin)/deltar) + 1
    if (ibin.gt.0 .and. ibin.le.nbins) then
       rbin = rmin + (ibin-0.5)*deltar

       area = pi*((rbin + 0.5*deltar)**2 - (rbin - 0.5*deltar)**2)
!$omp atomic
       sigma(ibin) = sigma(ibin) + pmassi/area
       if (present(utherm)) then
!$omp atomic
          spsound(ibin) = spsound(ibin) + 2./3.*utherm(i)
!$omp atomic
          ninbin(ibin) = ninbin(ibin) + 1
       endif
    endif
 enddo
!$omp end do
!$omp end parallel

!
!--calculate Toomre Q parameter in each bin using surface density
! 
 if (iplot.eq.2) then
    epicyclic = 0.
    do ibin=1,nbins
       sigmai = sigma(ibin)
!
!--for Toomre Q need the epicyclic frequency
!  in a Keplerian disc kappa = Omega
!
       Omegai = sqrt(1./radius(ibin)**3)
       epicyclic = Omegai
!
!--spsound is RMS sound speed for all particles in the annulus
!
       if (ninbin(ibin).gt.0) then
          spsoundi = sqrt(spsound(ibin)/real(ninbin(ibin)))
       else
          spsoundi = 0.
       endif
!
!--now calculate Toomre Q
!
       toomreq = spsoundi*epicyclic/(pi*sigmai)
       sigma(ibin) = toomreq
    enddo
 endif
 
 if (itransx.gt.0) call transform(radius,itransx)
 if (itransy.gt.0) call transform(sigma,itransy)
!
!--return min and max of y axis so adaptive plot limits can be set
!
 ymin = minval(sigma(1:nbins))
 ymax = maxval(sigma(1:nbins))

 return
end subroutine disccalc

subroutine discplot
 implicit none
 
 call pgline(nbins,radius,sigma)

end subroutine discplot

end module disc
