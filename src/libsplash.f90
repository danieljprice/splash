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
!  Copyright (C) 2005-2019 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! module providing library version of splash interpolation routines
! specifies c interfaces to corresponding Fortran subroutines
!-------------------------------------------------------------------------
module libsplash

 use projections3D,         only: interpolate3d_projection_f => interpolate3d_projection, &
                                   interpolate3d_proj_vec_f => interpolate3d_proj_vec

 use xsections3D,           only: interpolate3d_fastxsec_f => interpolate3d_fastxsec, &
                                   interpolate3d_xsec_vec_f => interpolate3d_xsec_vec

 use interpolate3D_opacity, only: interp3d_proj_opacity_f => interp3d_proj_opacity

 use projections3Dgeom,     only: interpolate3D_proj_geom_f => interpolate3D_proj_geom, &
                                   interpolate3D_xsec_geom_f => interpolate3D_xsec_geom

 use iso_c_binding,         only: c_float, c_int, c_bool

 implicit none

contains
!-------------------------------------------------------------------------
! 3D projection (column density)
!-------------------------------------------------------------------------
subroutine interpolate3d_projection_c(                                         &
  x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx, npixy, &
  pixwidthx, pixwidthy, normalise, zobserver, dscreen, useaccelerate, iverbose &
  ) bind(c)

 real(c_float),   intent(in)  :: x(npart),      &
                                  y(npart),      &
                                  z(npart),      &
                                  hh(npart),     &
                                  weight(npart), &
                                  dat(npart),    &
                                  xmin,          &
                                  ymin,          &
                                  pixwidthx,     &
                                  pixwidthy,     &
                                  zobserver,     &
                                  dscreen
 integer(c_int),  intent(in)  :: npart,         &
                                  npixx,         &
                                  npixy,         &
                                  itype(npart),  &
                                  iverbose
 logical(c_bool), intent(in)  :: normalise,     &
                                  useaccelerate
 real(c_float),   intent(out) :: datsmooth(npixx,npixy)

 logical :: normalise_f, &
             useaccelerate_f

 normalise_f     = normalise
 useaccelerate_f = useaccelerate

 call interpolate3d_projection_f(                                        &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx, &
    npixy, pixwidthx, pixwidthy, normalise_f, zobserver, dscreen,         &
    useaccelerate_f, iverbose)

end subroutine interpolate3d_projection_c

!-------------------------------------------------------------------------
! 3D projection of vectors
!-------------------------------------------------------------------------
subroutine interpolate3d_proj_vec_c(                                     &
  x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, vecsmoothx, &
  vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise, zobserver,  &
  dscreen, iverbose                                                      &
  ) bind(c)

 integer(c_int),  intent(in)  :: npart,                   &
                                  npixx,                   &
                                  npixy,                   &
                                  itype(npart),            &
                                  iverbose
 real(c_float),   intent(in)  :: x(npart),                &
                                  y(npart),                &
                                  z(npart),                &
                                  hh(npart),               &
                                  weight(npart),           &
                                  vecx(npart),             &
                                  vecy(npart),             &
                                  xmin,                    &
                                  ymin,                    &
                                  pixwidthx,               &
                                  pixwidthy,               &
                                  zobserver,               &
                                  dscreen
 logical(c_bool), intent(in)  :: normalise
 real(c_float),   intent(out) :: vecsmoothx(npixx,npixy), &
                                  vecsmoothy(npixx,npixy)

 logical :: normalise_f

 normalise_f = normalise

 call interpolate3d_proj_vec_f(                                            &
    x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, vecsmoothx,  &
    vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise_f, zobserver, &
    dscreen, iverbose)

end subroutine interpolate3d_proj_vec_c

!-------------------------------------------------------------------------
! cross sections of 3D data
!-------------------------------------------------------------------------
subroutine interpolate3d_fastxsec_c(                                     &
  x, y, z, hh, weight, dat, itype, npart, xmin, ymin, zslice, datsmooth, &
  npixx, npixy, pixwidthx, pixwidthy, normalise, iverbose                &
  ) bind(c)

 integer(c_int),   intent(in) :: npart,         &
                                  npixx,         &
                                  npixy,         &
                                  itype(npart),  &
                                  iverbose
 real(c_float),   intent(in)  :: x(npart),      &
                                  y(npart),      &
                                  z(npart),      &
                                  hh(npart),     &
                                  weight(npart), &
                                  dat(npart),    &
                                  xmin,          &
                                  ymin,          &
                                  pixwidthx,     &
                                  pixwidthy,     &
                                  zslice
 logical(c_bool), intent(in)  :: normalise
 real(c_float),   intent(out) :: datsmooth(npixx,npixy)

 logical :: normalise_f

 normalise_f = normalise

 call interpolate3d_fastxsec_f(                                           &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, zslice, datsmooth, &
    npixx, npixy, pixwidthx, pixwidthy, normalise_f, iverbose)

end subroutine interpolate3d_fastxsec_c

!-------------------------------------------------------------------------
! cross sections of 3D vector data
!-------------------------------------------------------------------------
subroutine interpolate3d_xsec_vec_c(                                     &
  x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, zslice,     &
  vecsmoothx, vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise, &
  iverbose) bind(c)

 integer(c_int),  intent(in)  :: npart,                   &
                                  npixx,                   &
                                  npixy,                   &
                                  itype(npart),            &
                                  iverbose
 real(c_float),   intent(in)  :: x(npart),                &
                                  y(npart),                &
                                  z(npart),                &
                                  hh(npart),               &
                                  weight(npart),           &
                                  vecx(npart),             &
                                  vecy(npart),             &
                                  xmin,                    &
                                  ymin,                    &
                                  pixwidthx,               &
                                  pixwidthy,               &
                                  zslice
 logical(c_bool), intent(in)  :: normalise
 real(c_float),   intent(out) :: vecsmoothx(npixx,npixy), &
                                  vecsmoothy(npixx,npixy)

 logical :: normalise_f

 normalise_f = normalise

 call interpolate3d_xsec_vec_f(                                             &
    x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, zslice,       &
    vecsmoothx, vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise_f, &
    iverbose)

end subroutine interpolate3d_xsec_vec_c

!-------------------------------------------------------------------------
! opacity rendering of 3D data
!-------------------------------------------------------------------------
subroutine interp3d_proj_opacity_c(                                         &
  x, y, z, pmass, npmass, hh, weight, dat, zorig, itype, npart, xmin, ymin, &
  datsmooth, brightness, npixx, npixy, pixwidth, zobserver,                 &
  dscreenfromobserver, rkappa, zcut, iverbose                               &
  ) bind(c)

 integer(c_int), intent(in)  :: npart,                  &
                                 npixx,                  &
                                 npixy,                  &
                                 npmass,                 &
                                 itype(npart),           &
                                 iverbose
 real(c_float),  intent(in)  :: x(npart),               &
                                 y(npart),               &
                                 z(npart),               &
                                 hh(npart),              &
                                 weight(npart),          &
                                 dat(npart),             &
                                 zorig(npart),           &
                                 pmass(npmass),          &
                                 xmin,                   &
                                 ymin,                   &
                                 pixwidth,               &
                                 zobserver,              &
                                 dscreenfromobserver,    &
                                 zcut,                   &
                                 rkappa
 real(c_float),  intent(out) :: datsmooth(npixx,npixy), &
                                 brightness(npixx,npixy)

 call interp3d_proj_opacity_f(                                               &
    x, y, z, pmass, npmass, hh, weight, dat, zorig, itype, npart, xmin, ymin, &
    datsmooth, brightness, npixx, npixy, pixwidth, zobserver,                 &
    dscreenfromobserver, rkappa, zcut, iverbose)

end subroutine interp3d_proj_opacity_c

subroutine interpolate3D_proj_geom_c(                                          &
  x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx, npixy, &
  pixwidthx, pixwidthy, normalise, igeom, iplotx, iploty, iplotz, ix, xorigin  &
  ) bind(c)

 integer(c_int),  intent(in)    :: npart,         &
                                    npixx,         &
                                    npixy,         &
                                    itype(npart),  &
                                    igeom,         &
                                    iplotx,        &
                                    iploty,        &
                                    iplotz,        &
                                    ix(3)
 real(c_float),   intent(in)    :: x(npart),      &
                                    y(npart),      &
                                    z(npart),      &
                                    hh(npart),     &
                                    weight(npart), &
                                    dat(npart),    &
                                    xmin,          &
                                    ymin,          &
                                    pixwidthx,     &
                                    pixwidthy,     &
                                    xorigin(3)
 real(c_float),   intent(out)   :: datsmooth(npixx,npixy)
 logical(c_bool), intent(inout) :: normalise

 logical :: normalise_f

 normalise_f = normalise

 call interpolate3D_proj_geom_f(                                            &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx,    &
    npixy, pixwidthx, pixwidthy, normalise_f, igeom, iplotx, iploty, iplotz, &
    ix, xorigin)

end subroutine interpolate3D_proj_geom_c

subroutine interpolate3D_xsec_geom_c(                                          &
  x, y, z, hh, weight, dat, itype, npart, xmin, ymin, zslice, datsmooth,       &
  npixx, npixy, pixwidthx, pixwidthy, normalise, igeom, iplotx, iploty,        &
  iplotz, ix, xorigin) bind(c)

 integer(c_int),  intent(in)    :: npart,         &
                                    npixx,         &
                                    npixy,         &
                                    itype(npart),  &
                                    igeom,         &
                                    iplotx,        &
                                    iploty,        &
                                    iplotz,        &
                                    ix(3)
 real(c_float),   intent(in)    :: x(npart),      &
                                    y(npart),      &
                                    z(npart),      &
                                    hh(npart),     &
                                    weight(npart), &
                                    dat(npart),    &
                                    xmin,          &
                                    ymin,          &
                                    zslice,        &
                                    pixwidthx,     &
                                    pixwidthy,     &
                                    xorigin(3)
 real(c_float),   intent(out)   :: datsmooth(npixx,npixy)
 logical(c_bool), intent(inout) :: normalise

 logical :: normalise_f

 normalise_f = normalise

 call interpolate3D_xsec_geom_f(                                           &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, zslice, datsmooth,  &
    npixx, npixy, pixwidthx, pixwidthy, normalise_f, igeom, iplotx, iploty, &
    iplotz, ix, xorigin)

end subroutine interpolate3D_xsec_geom_c

end module libsplash
