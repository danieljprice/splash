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
!  Copyright (C) 2005-2020 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module providing library version of splash exact routines
! specifies c interfaces to corresponding Fortran subroutines
!-------------------------------------------------------------------------
module libexact

 use shock,            only:exact_shock
 use shock_sr,         only:exact_shock_sr
 use sedov,            only:exact_sedov
 use polytrope,        only:exact_polytrope
 use toystar1D,        only:exact_toystar1D
 use toystar2D,        only:exact_toystar2D
 use gresho,           only:exact_gresho
 use mhdshock,         only:exact_mhdshock
 use rhoh,             only:exact_rhoh
 use densityprofiles,  only:exact_densityprofiles
 use torus,            only:exact_torus
 use ringspread,       only:exact_ringspread
 use dustywaves,       only:exact_dustywave
 use rochelobe,        only:exact_rochelobe
 use cshock,           only:exact_cshock
 use planetdisc,       only:exact_planetdisc
 use bondi,            only:exact_bondi

 use iso_c_binding,    only: c_float, c_int, c_bool

 implicit none

 public

contains
subroutine check_argcv_f() bind(c)
 include 'libinclude.f90'
end subroutine check_argcv_f

subroutine shock_c(&
    iplot,npart,time,gamma,xshock,rho_L,rho_R,p_L,p_R,v_L,v_R,&
    rdust_to_gas,xplot,yplot,ierr) bind(c, name='shock')
 integer(c_int), intent(in)  :: iplot, npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: time,gamma,xshock
 real(c_float),  intent(in)  :: rho_L,rho_R,p_L,p_R,v_L,v_R,rdust_to_gas
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call shock(iplot,time,gamma,xshock,rho_L,rho_R,p_L,p_R,v_L,v_R,&
            rdust_to_gas,xplot,yplot,ierr)

end subroutine shock_c

subroutine shock_sr_c(iplot,npart,time,gamma,rho_L,rho_R,p_L,p_R,v_L,v_R,&
                      xplot,yplot,ierr) bind(c, name='shock_sr')
 integer(c_int), intent(in)    :: iplot, npart
 integer(c_int), intent(out)   :: ierr
 real(c_float),  intent(in)    :: time,gamma
 real(c_float),  intent(in)    :: rho_L,rho_R,p_L,p_R,v_L,v_R
 real(c_float),  intent(inout) :: xplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call shock_sr(iplot,time,gamma,rho_L,rho_R,p_L,p_R,v_L,v_R,xplot,yplot,ierr)

end subroutine shock_sr_c

subroutine sedov_c(iplot,npart,time,gamma,rhozero,energy,rmax,&
                   rplot,yplot,ierr) bind(c, name='sedov')
 integer(c_int), intent(in)    :: iplot, npart
 integer(c_int), intent(out)   :: ierr
 real(c_float),  intent(in)    :: time,gamma
 real(c_float),  intent(in)    :: rhozero,energy,rmax
 real(c_float),  intent(inout) :: rplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call sedov(iplot,time,gamma,rhozero,energy,rmax,rplot,yplot,ierr)

end subroutine sedov_c

subroutine polytrope_c(npart,gamma,polyk,totmass,rplot,yplot,&
                       npartout,ierr) bind(c, name='polytrope')
 integer(c_int), intent(in)    :: npart
 integer(c_int), intent(out)   :: ierr,npartout
 real(c_float),  intent(in)    :: gamma,polyk,totmass
 real(c_float),  intent(inout) :: rplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call polytrope(gamma,polyk,totmass,rplot,yplot,npartout,ierr)

end subroutine polytrope_c

subroutine toystar1D_c(iplot,npart,time,gamma,H0,A0,C0,sigma,norder,&
                      xplot,yplot,ierr) bind(c, name='toystar1D')
 integer(c_int), intent(in)    :: iplot,npart,norder
 integer(c_int), intent(out)   :: ierr
 real(c_float),  intent(in)    :: time,gamma,sigma,H0,A0,C0
 real(c_float),  intent(inout) :: xplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call toystar1D(iplot,time,gamma,H0,A0,C0,sigma,&
                norder,xplot,yplot,npart,ierr)

end subroutine toystar1D_c

subroutine toystar2D_c(iplot,npart,time,gamma,polyk,totmass,ampl,&
                       denscentre,C0,jorder,morder,V11,V22,V12,V21,&
                       xplot,yplot,ierr) bind(c, name='toystar2D')
 integer(c_int), intent(in)   :: iplot,npart,jorder,morder
 integer(c_int), intent(out)  :: ierr
 real(c_float),  intent(in)   :: time,gamma,polyk,totmass,&
                                    C0, ampl, denscentre,&
                                    V11,V22,V12,V21
 real(c_float), intent(inout) :: xplot(npart)
 real(c_float), intent(out)   :: yplot(npart)

 call toystar2D(iplot,time,gamma,polyk,totmass, &
                ampl,denscentre,C0,jorder,morder, &
                V11,V22,V12,V21,xplot,yplot,ierr)

end subroutine toystar2D_c

subroutine gresho_c(iplot,npart,xplot,yplot,ierr) bind(c, name='gresho')
 integer(c_int), intent(in)  :: iplot,npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call gresho(iplot,xplot,yplot,ierr)
end subroutine gresho_c

subroutine mhdshock_c(iplot,npart,ishk,time,gamma,xmin,xmax,xshock,&
                      xplot,yplot,npts,ierr) bind(c, name='mhdshock')
 integer(c_int), intent(in)    :: iplot,npart,ishk
 integer(c_int), intent(out)   :: npts,ierr
 real(c_float),  intent(in)    :: time,gamma,xmin,xmax,xshock
 real(c_float),  intent(inout) :: xplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call mhdshock(iplot,ishk,time,gamma,xmin,xmax,xshock,xplot,yplot,npts,ierr)

end subroutine mhdshock_c

subroutine rhoh_c(iplot,npart,ndim,hfact,pmassval,&
                 xplot,yplot,ierr) bind(c, name='rhoh')
 integer(c_int), intent(in)  :: iplot,ndim,npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: hfact,pmassval
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call rhoh(iplot,ndim,hfact,pmassval,xplot,yplot,ierr)

end subroutine rhoh_c

subroutine densityprofiles_c(iplot,npart,iprofile,Mspherex,Mspherey,&
                             rsoftx,rsofty,xplot,yplot,ierr)&
                             bind(c, name='denstyprofiles')
 integer(c_int), intent(in)  :: iplot,iprofile,npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: Mspherex,Mspherey,rsoftx,rsofty
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call densityprofiles(iplot,iprofile,[Mspherex,Mspherey],&
      [rsoftx,rsofty],xplot,yplot,ierr)

end subroutine densityprofiles_c

subroutine torus_c(iplot,npart,itorus,Mstar,Rtorus,AA,distortion,&
                   gamma,xplot,yplot,ierr) bind(c, name='torus')
 integer(c_int), intent(in)  :: iplot,itorus,npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: Mstar,Rtorus,AA,gamma,distortion
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call torus(iplot,itorus,Mstar,Rtorus,AA,distortion,gamma,xplot,yplot,ierr)

end subroutine torus_c

subroutine ringspread_c(iplot,npart,time,Mdisk,Rdisk,viscnu,&
                      xplot,yplot,ierr) bind(c, name='ringspread')
 integer(c_int), intent(in)  :: iplot,npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: time,Mdisk,Rdisk,viscnu
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call ringspread(iplot,time,Mdisk,Rdisk,viscnu,xplot,yplot,ierr)

end subroutine ringspread_c

subroutine dustywave_c(iplot,npart,time,ampl,cs,Kdragin,lambda,x0,&
                       rhog0,rhod0,xplot,yplot,ierr) bind(c, name='dustywave')
 integer(c_int), intent(in)  :: iplot,npart
 integer(c_int), intent(out) :: ierr
 real(c_float),  intent(in)  :: time, ampl, cs, Kdragin, lambda, x0, rhog0, rhod0
 real(c_float),  intent(in)  :: xplot(npart)
 real(c_float),  intent(out) :: yplot(npart)

 call dustywave(iplot,time,ampl,cs,Kdragin,lambda,x0,&
                rhog0,rhod0,xplot,yplot,ierr)

end subroutine dustywave_c

subroutine rochelobe_c(npart,x1,y1,x2,y2,m1,m2,&
                       xplot,yplot,ierr) bind(c, name='rochelobe')
 integer(c_int), intent(in)    :: npart
 integer(c_int), intent(out)   :: ierr
 real(c_float),  intent(in)    :: x1,y1,x2,y2,m1,m2
 real(c_float),  intent(inout) :: xplot(npart),yplot(npart)

 call rochelobe(x1,y1,x2,y2,m1,m2,xplot,yplot,ierr)
 ierr = 0

end subroutine rochelobe_c

subroutine cshock_c(iplot,npart,time,gamma,machs,macha,xmin,xmax,&
                    xplot,yplot,ierr) bind(c, name='cshock')
 integer(c_int), intent(in)    :: iplot,npart
 integer(c_int), intent(out)   :: ierr
 real(c_float),  intent(in)    :: time,gamma,machs,macha,xmin,xmax
 real(c_float),  intent(inout) :: xplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call cshock(iplot,time,gamma,machs,macha,xmin,xmax,xplot,yplot,ierr)

end subroutine cshock_c

subroutine planetdisc_c(iplot,npart,ispiral,time,HonR,rplanet,q,narms,&
                        params,rplot,yplot,ierr) bind(c, name='planetdisc')
 integer(c_int), intent(in)    :: iplot,ispiral,narms,npart
 integer(c_int), intent(out)   :: ierr
 real(c_float),  intent(in)    :: time, HonR, rplanet, q, params(7,10)
 real(c_float),  intent(inout) :: rplot(npart)
 real(c_float),  intent(out)   :: yplot(npart)

 call planetdisc_f(iplot,ispiral,time,HonR,rplanet,q,narms,&
                   params,rplot,yplot,ierr)
 ierr = 0
end subroutine planetdisc_c

subroutine bondi_c(iplot,npart,time,gamma,const1,const2,m,relativistic,&
                   geodesic_flow,is_wind,xplot,yplot,ierr) bind(c, name='bondi')
 integer(c_int),  intent(in)  :: iplot,npart
 integer(c_int),  intent(out) :: ierr
 real(c_float),   intent(in)  :: time,gamma,const1,const2,m
 logical(c_bool), intent(in)  :: relativistic, geodesic_flow,is_wind
 real(c_float),   intent(in)  :: xplot(npart)
 real(c_float),   intent(out) :: yplot(npart)

 logical :: relativistic_f,geodesic_flow_f,is_wind_f

 relativistic_f  = relativistic
 geodesic_flow_f = geodesic_flow
 is_wind_f       = is_wind

 call bondi(iplot,time,gamma,const1,const2,m,relativistic_f,&
            geodesic_flow_f,is_wind_f,xplot,yplot,ierr)

end subroutine bondi_c
end module libexact
