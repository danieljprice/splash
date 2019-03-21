module libsplash

  use projections3D,         only: interpolate3d_projection, &
                                   interpolate3d_proj_vec

  use xsections3D,           only: interpolate3d_fastxsec,   &
                                   interpolate3d_xsec_vec

  use interpolate3D_opacity, only: interp3d_proj_opacity

  use iso_c_binding,         only: c_float, c_int, c_bool

  implicit none

contains

subroutine c_interpolate3d_projection(                                         &
  x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx, npixy, &
  pixwidthx, pixwidthy, normalise, zobserver, dscreen, useaccelerate           &
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
                                  itype(npart)
  logical(c_bool), intent(in)  :: normalise,     &
                                  useaccelerate
  real(c_float),   intent(out) :: datsmooth(npixx,npixy)

  logical :: normalise_f, &
             useaccelerate_f

  normalise_f     = normalise
  useaccelerate_f = useaccelerate

  call interpolate3d_projection(                                          &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx, &
    npixy, pixwidthx, pixwidthy, normalise_f, zobserver, dscreen,         &
    useaccelerate_f)

end subroutine c_interpolate3d_projection

subroutine c_interpolate3d_proj_vec(                                     &
  x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, vecsmoothx, &
  vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise, zobserver,  &
  dscreen                                                                &
  ) bind(c)

  integer(c_int),  intent(in)  :: npart,                   &
                                  npixx,                   &
                                  npixy,                   &
                                  itype(npart)
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

  call interpolate3d_proj_vec(                                              &
    x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, vecsmoothx,  &
    vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise_f, zobserver, &
    dscreen)

end subroutine c_interpolate3d_proj_vec

subroutine c_interpolate3d_fastxsec(                                     &
  x, y, z, hh, weight, dat, itype, npart, xmin, ymin, zslice, datsmooth, &
  npixx, npixy, pixwidthx, pixwidthy, normalise                          &
  ) bind(c)

  integer(c_int),   intent(in) :: npart,         &
                                  npixx,         &
                                  npixy,         &
                                  itype(npart)
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

  call interpolate3d_fastxsec(                                             &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, zslice, datsmooth, &
    npixx, npixy, pixwidthx, pixwidthy, normalise_f)

end subroutine c_interpolate3d_fastxsec

subroutine c_interpolate3d_xsec_vec(                                     &
  x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, zslice,     &
  vecsmoothx, vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise  &
  ) bind(c)

  integer(c_int),  intent(in)  :: npart,                   &
                                  npixx,                   &
                                  npixy,                   &
                                  itype(npart)
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

  call interpolate3d_xsec_vec(                                               &
    x, y, z, hh, weight, vecx, vecy, itype, npart, xmin, ymin, zslice,       &
    vecsmoothx, vecsmoothy, npixx, npixy, pixwidthx, pixwidthy, normalise_f)

end subroutine c_interpolate3d_xsec_vec

subroutine c_interp3d_proj_opacity(                                         &
  x, y, z, pmass, npmass, hh, weight, dat, zorig, itype, npart, xmin, ymin, &
  datsmooth, brightness, npixx, npixy, pixwidth, zobserver,                 &
  dscreenfromobserver, rkappa, zcut                                         &
  ) bind(c)

  integer(c_int), intent(in)  :: npart,                  &
                                 npixx,                  &
                                 npixy,                  &
                                 npmass,                 &
                                 itype(npart)
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

  call interp3d_proj_opacity(x, y, z, pmass, npmass, hh, weight, dat, zorig, &
    itype, npart, xmin, ymin, datsmooth, brightness, npixx, npixy, pixwidth, &
    zobserver, dscreenfromobserver, rkappa, zcut)

end subroutine c_interp3d_proj_opacity

end module libsplash
