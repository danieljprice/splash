!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  module handling plotting of azimuthally-averaged quantities
!  for disc simulations
!-----------------------------------------------------------------
module disc
 implicit none
 integer, parameter, private :: maxbins = 1001
 integer, dimension(maxbins), private :: ninbin
 real, dimension(maxbins), private :: radius,yplot,sigma,spsound,omegasq
 integer, private :: nbins

 public :: disccalc,discplot

 private

contains

subroutine disccalc(iplot,npart,rpart,pmass_in,rminin,rmaxin,ymin,ymax,&
                    itransx,itransy,icolourpart,iamtype,usetype,noftype,gamma,mstar,&
                    dat,units,units_dz,iRescale,ipmass,irho,ispsound,iutherm,ix,ivx)
 use transforms, only:transform_limits_inverse,transform_inverse,transform
 use params,     only:int1,maxparttypes,doub_prec
 use physcon,    only:pi
 use part_utils, only:igettype
 integer,                          intent(in)  :: iplot,npart,itransx,itransy
 real, dimension(npart),           intent(in)  :: rpart
 real,                             intent(in)  :: pmass_in,rminin,rmaxin,gamma,mstar
 real,                             intent(out) :: ymin,ymax
 integer, dimension(npart),        intent(in)  :: icolourpart
 integer(kind=int1), dimension(:), intent(in)  :: iamtype
 logical, dimension(maxparttypes), intent(in)  :: usetype
 integer, dimension(maxparttypes), intent(in)  :: noftype
 real(doub_prec),                  intent(in)  :: units(:),units_dz
 logical,                          intent(in)  :: iRescale
 real, dimension(:,:),             intent(in)  :: dat
 integer,                          intent(in)  :: ipmass,irho,ispsound,iutherm,ix(:),ivx

 integer :: i,ibin
 real :: pmassi,rbin,deltar,area,rmin,rmax
 real :: x,y,vx,vy,r,vphi
 real(doub_prec) :: sigmai,toomreq,epicyclic,Omegai,spsoundi,dOmegasq_dr
 real(doub_prec) :: unit_mass,unit_r,unit_u,unit_dz,unit_dens,unit_cs2
 real, dimension(1) :: rad
 logical :: mixedtypes
 integer :: itype,np

 ninbin(:) = 0
 sigma(:) = 0.
 spsound(:) = 0.
 omegasq(:) = 0.
 pmassi = 0.

 ! work out the unit of mass, r needed for computing Toomre Q
 unit_mass = 1.d0
 unit_r    = 1.d0
 unit_u    = 1.d0
 unit_dz   = 1.d0
 unit_dens = 1.d0
 unit_cs2 = 1.d0
 unit_dz  = 1.d0
 if (iRescale) then
    if (ix(1) > 0) unit_r = units(ix(1))
    if (iutherm > 0) unit_u = units(iutherm)
    if (irho > 0) unit_dens = units(irho)
    if (ipmass > 0) then
       unit_mass = units(ipmass)
    elseif (irho > 0) then
       unit_mass = units(irho)*units_dz**3
    endif
    unit_dz = units_dz
    if (ispsound > 0) unit_cs2 = units(ispsound)
 endif
!
!--print info
!
 select case(iplot)
 case(1)
    print "(a)",' calculating disc surface density profile '
 case(2)
    if (iutherm > 0 .or. ispsound > 0) then
       print "(a,es10.3,a)",' calculating Toomre Q parameter (assuming Mstar=',mstar,' and Keplerian rotation)'
       if (ispsound == 0) then
          if (gamma < 1.00001) then
             print "(a)",' isothermal equation of state: using cs^2 = 2/3*utherm'
          else
             print "(a,f6.3,a,f6.3,a)",' ideal gas equation of state: using cs^2 = ',gamma*(gamma-1),'*u (gamma = ',gamma,')'
          endif
       endif
    else
       print "(a)",' ERROR: cannot calculate Toomre Q parameter: thermal energy/sound speed not present in dump file'
       return
    endif
 case(3)
    print "(a)",' calculating kappa^2/Omega^2 '
 case default
    print "(a)",' ERROR: unknown plot in discplot. '
    return
 end select
!
!--if transformations (e.g. log) are applied to r, then limits
!  will already be set in transformed space - need to obtain
!  limits in non-transformed space.
!
 rmin = rminin
 rmax = rmaxin
 if (itransx > 0) call transform_limits_inverse(rmin,rmax,itransx)
!
!--try to get appropriate value for nbins
!
 nbins = min(4*int(npart**(1./3.)) + 1,maxbins)
!
!--set array of radius values for plotting
!
 deltar = (rmax - rmin)/(nbins - 1)
 do ibin=1,nbins
    radius(ibin) = rmin + (ibin-0.5)*deltar
 enddo
 mixedtypes = size(iamtype) >= npart
!
!--calculate surface density in each radial bin
!
 np = 0
!$omp parallel do default(none) &
!$omp shared(npart,rpart,sigma,omegasq,pmass_in,itransx,icolourpart,rmin,deltar,nbins,mstar) &
!$omp shared(ninbin,spsound,gamma,iamtype,mixedtypes,usetype,noftype) &
!$omp shared(dat,unit_cs2,ispsound,iutherm,ivx,ipmass,ix) &
!$omp private(i,rad,pmassi,ibin,rbin,area,itype,x,y,vx,vy,r,vphi) &
!$omp reduction(+:np)
 over_parts: do i=1,npart
    !--skip particles with itype < 0
    if (icolourpart(i) < 0) cycle over_parts
    if (mixedtypes) then
       itype = int(iamtype(i))
    else
       itype = igettype(i,noftype)
    endif
    if (.not.usetype(itype)) cycle over_parts
    np = np + 1

    if (itransx==0) then
       rad(1) = rpart(i)
    else
       rad(1) = rpart(i)
       call transform_inverse(rad,itransx)
    endif
    if (ipmass > 0) then
       pmassi = dat(i,ipmass)
    else
       pmassi = pmass_in
    endif
    ibin = int((rad(1) - rmin)/deltar) + 1
    if (ibin > 0 .and. ibin <= nbins) then
       rbin = rmin + (ibin-0.5)*deltar

       area = pi*((rbin + 0.5*deltar)**2 - (rbin - 0.5*deltar)**2)
       !$omp atomic
       sigma(ibin) = sigma(ibin) + pmassi/area
       if (ispsound > 0) then
          !$omp atomic
          spsound(ibin) = spsound(ibin) + real((dat(i,ispsound))**2/unit_cs2)
       elseif (iutherm > 0) then
          if (gamma < 1.00001) then
             !$omp atomic
             spsound(ibin) = spsound(ibin) + real(2./3.*(dat(i,iutherm)/unit_cs2))
          else
             !$omp atomic
             spsound(ibin) = spsound(ibin) + real(gamma*(gamma-1.)*(dat(i,iutherm)/unit_cs2))
          endif
       endif
      if (ivx > 0 .and. ix(1) > 0 .and. ix(2) > 0) then
         !--angular velocity Omega = vphi / r with
         !  vphi from full x,y,vx,vy in the disc plane
         x  = dat(i,ix(1))
         y  = dat(i,ix(2))
         vx = dat(i,ivx)
         vy = dat(i,ivx+1)
         r  = sqrt(x*x + y*y)
         if (r > 0.) then
            vphi = (x*vy - y*vx)/r
            !$omp atomic
            omegasq(ibin) = omegasq(ibin) + (vphi/r)**2
         endif
      endif

       !$omp atomic
       ninbin(ibin) = ninbin(ibin) + 1
    endif
 enddo over_parts
!$omp end parallel do

 print "(1x,a,i10,a,i10,a,i4,a)",'used ',np,' of ',npart,' particles in ',nbins,' bins'
 select case(iplot)
 case(3)
!
!--average Keplerian rotation frequency in each bin
!
   do ibin=1,nbins
       omegasq(ibin) = omegasq(ibin)/real(ninbin(ibin))
   enddo
   do ibin=1,nbins
       ! compute epicyclic frequency r*dOmega^2/dr + 4 Omega^2
       if (ibin < nbins) then
          dOmegasq_dr = (omegasq(ibin+1) - omegasq(ibin))/(radius(ibin+1) - radius(ibin))
       else
          dOmegasq_dr = (omegasq(ibin) - omegasq(ibin-1))/(radius(ibin) - radius(ibin-1))
       endif
       epicyclic = radius(ibin)*dOmegasq_dr + 4*omegasq(ibin)
       yplot(ibin) = epicyclic/omegasq(ibin)
    enddo
 case(2)
!
!--calculate Toomre Q parameter in each bin using surface density
!
   epicyclic = 0.
    do ibin=1,nbins
       sigmai = sigma(ibin)*(unit_r**2/unit_mass)  ! convert back to code units
!
!--for Toomre Q need the epicyclic frequency
!  in a Keplerian disc kappa = Omega
!
       Omegai = sqrt(mstar/(radius(ibin)/unit_r)**3)
       epicyclic = Omegai
!
!--spsound is RMS sound speed for all particles in the annulus
!
       if (ninbin(ibin) > 0) then
          spsoundi = sqrt(spsound(ibin)/real(ninbin(ibin))) ! unit conversion already done
       else
          spsoundi = 0.
       endif
!
!--now calculate Toomre Q
!
       if (sigmai > 0.) then
          toomreq = spsoundi*epicyclic/(pi*sigmai)
       else
          toomreq = 0.
       endif
       yplot(ibin) = real(toomreq,kind=kind(yplot))
    enddo
 case default
!
!--return surface density in units of [g/cm^2], not [g/cm^3 au] or [Msun/au^2]
!
    yplot = sigma*(unit_r**2/unit_mass)  ! convert back to code units
    yplot = yplot*unit_dens*unit_dz      ! convert to units of unit_dens * unit_dz
 end select
!
!--give very small number instead of zero to avoid problems on log plots
!
 where (yplot(1:nbins) < tiny(0.)) yplot(1:nbins) = minval(yplot,mask=yplot>0.)

 if (itransx > 0) call transform(radius,itransx)
 if (itransy > 0) call transform(yplot,itransy)
!
!--return min and max of y axis so adaptive plot limits can be set
!
 ymin = minval(yplot(1:nbins),mask=(yplot(1:nbins) /= 0.))
 ymax = maxval(yplot(1:nbins),mask=(yplot(1:nbins) /= 0.))

end subroutine disccalc

!---------------------------------------------------
!
! subroutine to actually perform the disc plotting
!
!---------------------------------------------------
subroutine discplot
 use plotlib, only:plot_line

 call plot_line(nbins,radius,yplot)

end subroutine discplot

end module disc
