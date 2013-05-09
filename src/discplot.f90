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
 integer, parameter, private :: nbins = 350
 integer, dimension(nbins), private :: ninbin
 real, dimension(nbins), private :: radius,sigma,spsound

 public :: disccalc,discplot

 private

contains

subroutine disccalc(iplot,npart,rpart,npmass,pmass,rminin,rmaxin,ymin,ymax,&
                    itransx,itransy,icolourpart,iamtype,usetype,noftype,gamma,u,u_is_spsound)
 use transforms, only:transform_limits_inverse,transform_inverse,transform
 use params,     only:int1,maxparttypes
 use part_utils, only:igettype
 implicit none
 integer,                          intent(in)  :: iplot,npart,npmass,itransx,itransy
 real, dimension(npart),           intent(in)  :: rpart
 real, dimension(npmass),          intent(in)  :: pmass
 real,                             intent(in)  :: rminin,rmaxin,gamma
 real,                             intent(out) :: ymin,ymax
 integer, dimension(npart),        intent(in)  :: icolourpart
 integer(kind=int1), dimension(:), intent(in)  :: iamtype
 logical, dimension(maxparttypes), intent(in)  :: usetype
 integer, dimension(maxparttypes), intent(in)  :: noftype
 real, dimension(npart),           intent(in), optional :: u
 logical,                          intent(in), optional :: u_is_spsound

 integer :: i,ibin
 real, parameter :: pi = 3.1415926536
 real :: pmassi,rbin,deltar,area,rmin,rmax
 real :: sigmai,toomreq,epicyclic,Omegai,spsoundi
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
!$omp shared(npart,rpart,sigma,npmass,pmass,itransx,icolourpart,rmin,deltar) &
!$omp shared(ninbin,spsound,gamma,u,iamtype,mixedtypes,usetype,noftype,gotspsound) &
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
             spsound(ibin) = spsound(ibin) + u(i)
          else
             if (gamma.lt.1.00001) then
                !$omp atomic
                spsound(ibin) = spsound(ibin) + 2./3.*u(i)
             else
                !$omp atomic
                spsound(ibin) = spsound(ibin) + gamma*(gamma-1.)*u(i)
             endif
          endif
          !$omp atomic
          ninbin(ibin) = ninbin(ibin) + 1
       endif
    endif
 enddo over_parts
!$omp end parallel do

 print "(1x,a,i10,a,i10,a)",'used ',np,' of ',npart,' particles'

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
