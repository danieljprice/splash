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
 real, dimension(maxbins), private :: radius,sigma,spsound
 integer, private :: nbins

 public :: disccalc,discplot

 private

contains

subroutine disccalc(iplot,npart,rpart,npmass,pmass,unit_mass,unit_r,rminin,rmaxin,ymin,ymax,&
                    itransx,itransy,icolourpart,iamtype,usetype,noftype,gamma,unit_u,u,u_is_spsound)
 use transforms, only:transform_limits_inverse,transform_inverse,transform
 use params,     only:int1,maxparttypes,doub_prec
 use part_utils, only:igettype
 implicit none
 integer,                          intent(in)  :: iplot,npart,npmass,itransx,itransy
 real, dimension(npart),           intent(in)  :: rpart
 real, dimension(npmass),          intent(in)  :: pmass
 real(doub_prec),                  intent(in)  :: unit_mass,unit_r
 real,                             intent(in)  :: rminin,rmaxin,gamma
 real,                             intent(out) :: ymin,ymax
 integer, dimension(npart),        intent(in)  :: icolourpart
 integer(kind=int1), dimension(:), intent(in)  :: iamtype
 logical, dimension(maxparttypes), intent(in)  :: usetype
 integer, dimension(maxparttypes), intent(in)  :: noftype
 real(doub_prec),                  intent(in), optional :: unit_u
 real, dimension(npart),           intent(in), optional :: u
 logical,                          intent(in), optional :: u_is_spsound

 integer :: i,ibin
 real, parameter :: pi = 3.1415926536
 real :: pmassi,rbin,deltar,area,rmin,rmax
 real(doub_prec) :: sigmai,toomreq,epicyclic,Omegai,spsoundi,unit_cs2
 real, dimension(1) :: rad
 logical :: mixedtypes,gotspsound
 integer :: itype,np

 ninbin(:) = 0
 sigma(:) = 0.
 spsound(:) = 0.

 pmassi = 0
 if (npmass.le.0) then
    print*,' INTERNAL ERROR in discplot: dimension of mass array <= 0'
    return
 endif
 gotspsound = .false.
 if (present(u_is_spsound)) gotspsound = u_is_spsound
 if (present(unit_u)) then
    unit_cs2 = unit_u
 else
    unit_cs2 = 1.d0
 endif
!
!--print info
!
 select case(iplot)
 case(1)
    print "(a,i4,a)",' calculating disc surface density profile using',nbins,' bins'
 case(2)
    if (present(u)) then
       print "(a)",' calculating Toomre Q parameter (assuming Mstar=1 and a Keplerian rotation profile)'
       if (.not.gotspsound) then
          if (gamma.lt.1.00001) then
             print "(a)",' isothermal equation of state: using cs^2 = 2/3*utherm'
          else
             print "(a,f6.3,a,f6.3,a)",' ideal gas equation of state: using cs^2 = ',gamma*(gamma-1),'*u (gamma = ',gamma,')'
          endif
       endif
    else
       print "(a)",' ERROR: cannot calculate Toomre Q parameter: thermal energy/sound speed not present in dump file'
       return
    endif
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
 if (itransx.gt.0) call transform_limits_inverse(rmin,rmax,itransx)
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
 mixedtypes = size(iamtype).ge.npart
!
!--calculate surface density in each radial bin
!
 np = 0
!$omp parallel do default(none) &
!$omp shared(npart,rpart,sigma,npmass,pmass,itransx,icolourpart,rmin,deltar,nbins) &
!$omp shared(ninbin,spsound,gamma,u,iamtype,mixedtypes,usetype,noftype,gotspsound,unit_cs2) &
!$omp private(i,rad,pmassi,ibin,rbin,area,itype) &
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
       if (present(u)) then
          if (gotspsound) then
             !$omp atomic
             spsound(ibin) = spsound(ibin) + real((u(i))**2/unit_cs2)
          else
             if (gamma.lt.1.00001) then
                !$omp atomic
                spsound(ibin) = spsound(ibin) + real(2./3.*(u(i)/unit_cs2))
             else
                !$omp atomic
                spsound(ibin) = spsound(ibin) + real(gamma*(gamma-1.)*(u(i)/unit_cs2))
             endif
          endif
          !$omp atomic
          ninbin(ibin) = ninbin(ibin) + 1
       endif
    endif
 enddo over_parts
!$omp end parallel do

 print "(1x,a,i10,a,i10,a,i4,a)",'used ',np,' of ',npart,' particles in ',nbins,' bins'

!
!--calculate Toomre Q parameter in each bin using surface density
!
 if (iplot.eq.2) then
    epicyclic = 0.
    do ibin=1,nbins
       sigmai = sigma(ibin)*(unit_r**2/unit_mass)  ! convert back to code units
!
!--for Toomre Q need the epicyclic frequency
!  in a Keplerian disc kappa = Omega
!
       Omegai = sqrt(1./(radius(ibin)/unit_r)**3)
       epicyclic = Omegai
!
!--spsound is RMS sound speed for all particles in the annulus
!
       if (ninbin(ibin).gt.0) then
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
       sigma(ibin) = real(toomreq,kind=kind(sigma))
    enddo
 endif

 if (itransx.gt.0) call transform(radius,itransx)
 if (itransy.gt.0) call transform(sigma,itransy)
!
!--return min and max of y axis so adaptive plot limits can be set
!
 ymin = minval(sigma(1:nbins),mask=(sigma(1:nbins).ne.0.))
 ymax = maxval(sigma(1:nbins),mask=(sigma(1:nbins).ne.0.))

 return
end subroutine disccalc

!---------------------------------------------------
!
! subroutine to actually perform the disc plotting
!
!---------------------------------------------------
subroutine discplot
 use plotlib, only:plot_line
 implicit none

 call plot_line(nbins,radius,sigma)

end subroutine discplot

end module disc
