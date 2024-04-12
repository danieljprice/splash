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
!  Copyright (C) 2022- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module extinction
 implicit none

 public :: get_extinction, get_extinction_los

 private

contains
!---------------------------------------------------------
! routine to compute line-of-sight extinction to a
! collection of points (e.g. positions of sink particles)
!---------------------------------------------------------
subroutine get_extinction(ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,&
                          npts,xpts,ypts,zpts,coldens,icol)
 use labels,           only:ix,ih,irho,ipmass
 use limits,           only:get_particle_subset
 use projections3D,    only:interpolate3D_proj_points
 use particle_data,    only:icolourme
 use interpolation,    only:get_n_interp,set_interpolation_weights
 use settings_data,    only:iRescale,iverbose,required,UseTypeInRenderings
 use settings_part,    only:iplotpartoftype
 use settings_render,  only:npix,inormalise=>inormalise_interpolations,&
                            idensityweightedinterpolation,exact_rendering
 use settings_xsecrot, only:anglex,angley,anglez
 use settings_units,   only:units,unit_interp,unitzintegration
 use rotation,         only:rotate_particles
 use params,           only:int1
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 integer, intent(in)  :: npts,icol
 real,    intent(inout) :: xpts(npts),ypts(npts),zpts(npts)
 real,    intent(out) :: coldens(:)
 integer :: n,ierr,mycol
 real, dimension(:), allocatable :: weight,x,y,z

 !
 !--sanity checks
 !
 if (ndim /= 3) then
    print "(a)",' ERROR: extinction only works with 3 dimensional data'
    return
 endif
 if (.not. (ih > 0 .and. ipmass > 0 .and. irho > 0)) then
    print "(a)",' ERROR: could not locate h,mass or rho in data'
    return
 endif
 if (icol < 0 .or. icol > ncolumns) then
    print "(a)",' ERROR: invalid column number for extinction'
    return
 endif
 mycol = icol
 if (mycol==0) mycol = irho
 !
 !--set number of particles to use in the interpolation routines
 !  and allocate memory for weights
 !
 n = get_n_interp(ntypes,npartoftype,UseTypeInRenderings,iplotpartoftype,size(itype),.false.)
 allocate(weight(n),x(n),y(n),z(n),stat=ierr)
 if (ierr /= 0) then
    print*,' ERROR allocating memory for interpolation weights, aborting...'
    return
 endif
 x(1:n) = dat(1:n,ix(1))
 y(1:n) = dat(1:n,ix(2))
 z(1:n) = dat(1:n,ix(3))
 !
 !--rotate positions if necessary
 !
 call rotate_particles(n,x,y,z,anglex,angley,anglez)
 !
 !--set interpolation weights (w = m/(rho*h^ndim)
 !
 call set_interpolation_weights(weight,dat,itype,(iplotpartoftype .and. UseTypeInRenderings),&
      n,npartoftype,masstype,ntypes,ncolumns,irho,ipmass,ih,ndim,iRescale,&
      idensityweightedinterpolation,inormalise,units,unit_interp,required,.false.,isinktype=0)
 !
 !--set default mask and apply range restrictions to data
 !
 icolourme(:) = 1
 call get_particle_subset(icolourme,dat,ncolumns)

 !
 !--rotate positions if necessary
 !
 call rotate_particles(npts,xpts,ypts,zpts,anglex,angley,anglez)
 !
 ! interpolate column density (or other column-integrated quantity) 
 ! to desired point locations
 !
 call interpolate3D_proj_points(x,y,z,dat(1:n,ih),weight,dat(1:n,mycol),icolourme,n,&
                                xpts,ypts,zpts,coldens,npts,inormalise,iverbose)
 !
 ! adjust units of z-integrated quantity
 !
 if (iRescale .and. units(ih) > 0. .and. .not.inormalise) then
    coldens = coldens*(unitzintegration/units(ih))
 endif

 deallocate(x,y,z,weight)

end subroutine get_extinction

!--------------------------------------------------------------
! routine to compute line-of-sight extinction from a
! given point (currently the origin) along a set of directions
!--------------------------------------------------------------
subroutine get_extinction_los(ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,&
                        ndirs,anglex_vals,angley_vals,anglez_vals,coldens,icol)
 use settings_xsecrot, only:anglex,angley,anglez
 use params,           only:int1
 integer, intent(in)  :: ncolumns,ntypes,ndim
 integer, intent(in)  :: npartoftype(:)
 integer(kind=int1), intent(in) :: itype(:)
 real,    intent(in)  :: masstype(:)
 real,    intent(in)  :: dat(:,:)
 integer, intent(in)  :: ndirs,icol
 real,    intent(in)  :: anglex_vals(ndirs),angley_vals(ndirs),anglez_vals(ndirs)
 real,    intent(out) :: coldens(ndirs)
 integer, parameter :: npts = 1
 integer :: i
 real :: xpts(npts),ypts(npts),zpts(npts),coldensi(1)

 ! set the position from which to compute extinction as the origin
 xpts = 0.; ypts = 0.; zpts = 0.

 do i=1,ndirs
    anglex = anglex_vals(i)
    angley = angley_vals(i)
    anglez = anglez_vals(i)
    call get_extinction(ncolumns,dat,npartoftype,masstype,itype,ndim,ntypes,&
                        npts,xpts,ypts,zpts,coldensi,icol)
    coldens(i) = coldensi(1)
 enddo

 end subroutine get_extinction_los

end module extinction
