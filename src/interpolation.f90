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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  Module containing utility routines for SPH kernel interpolation
!-----------------------------------------------------------------
module interpolation
 implicit none
 public :: set_interpolation_weights
 real, parameter, public :: weight_sink = -1.
 private

contains
!-------------------------------------------------------------------
! Set interpolation weights for the particles. The weights are
! calculated using:
!
!           w = m/(rho*h**ndim),
!
! where we need to handle a few special scenarios:
!
! 1) Firstly, the weight should be calculated in a consistent
!    set of units. Safest way is to use the data as originally
!    read from the dump file, before any unit scaling was applied.
!
! 2) Particle weights are set to zero for particle types not
!    used in the rendering.
!
! 3) If particle mass not read, it is still possible to perform
!    interpolations, but not using the SPH weights. These
!    interpolations *must* therefore be normalised.
!-------------------------------------------------------------------
subroutine set_interpolation_weights(weighti,dati,iamtypei,usetype, &
           ninterp,npartoftype,masstype,ntypes,ndataplots,irho,ipmass,ih,ndim, &
           iRescale,idensityweighted,inormalise,units,unit_interp,required, &
           rendersinks)
  use params, only:doub_prec,int1,maxplot
  use labels, only:get_sink_type
  implicit none
  real, dimension(:), intent(out)              :: weighti
  real, dimension(:,:), intent(in)             :: dati
  integer(kind=int1), dimension(:), intent(in) :: iamtypei
  logical, dimension(:), intent(in)            :: usetype
  logical, dimension(0:maxplot), intent(in)    :: required
  integer, intent(in)                          :: ih,irho,ipmass,ndim
  integer, intent(in)                          :: ninterp,ntypes,ndataplots
  integer, dimension(:), intent(in)            :: npartoftype
  real,    dimension(:), intent(in)            :: masstype
  logical, intent(in)                          :: iRescale,idensityweighted,rendersinks
  logical, intent(inout)                       :: inormalise
  real, dimension(0:maxplot), intent(in)       :: units
  real(doub_prec), intent(in)                  :: unit_interp
  integer         :: i2,i1,itype,ipart,isinktype
  real(doub_prec) :: dunitspmass,dunitsrho,dunitsh

  !
  !-- unit_interp is a multiplication factor that can
  !   be used to scale the "weight" in case the
  !   units read from the read_data routine
  !   are inconsistent (this is the case for the SEREN read)
  !
  dunitspmass = 1.d0
  dunitsrho   = 1.d0
  dunitsh     = 1.d0
  if (iRescale) then
     if (ipmass.gt.0) dunitspmass = 1.d0/units(ipmass)
     if (ih.gt.0)     dunitsh     = 1.d0/units(ih)
     if (irho.gt.0)   dunitsrho   = 1.d0/units(irho)
  endif
  dunitspmass = dunitspmass * unit_interp
  
  isinktype = get_sink_type(ntypes)

  if (ipmass.gt.0 .and. ipmass.le.ndataplots .and. &
      irho.gt.0 .and. irho.le.ndataplots .and. &
      ih .gt. 0 .and. ih.le.ndataplots .and. &
      required(ipmass) .and. required(irho) .and. required(ih)) then

     if (size(iamtypei) > 1) then
        !
        !--particles with mixed types
        !
        !$omp parallel do default(none) &
        !$omp shared(ninterp,iamtypei,weighti,dati,rendersinks,isinktype) &
        !$omp shared(usetype,idensityweighted,dunitsrho)    &
        !$omp shared(ipmass,ih,irho,dunitspmass,dunitsh,ndim) &
        !$omp private(ipart,itype)
        do ipart=1,ninterp
           itype = iamtypei(ipart)
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(ipart) = weight_sink
              else
                 weighti(ipart) = 0.
              endif
           elseif (idensityweighted) then
              if (dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (dati(ipart,ipmass)*dunitspmass)/ &
                                 ((dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           else
              if (dati(ipart,irho) > tiny(dati) .and. dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (dati(ipart,ipmass)*dunitspmass)/ &
                                 ((dati(ipart,irho)*dunitsrho)*(dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           endif
        enddo
        !$omp end parallel do
     else
        !
        !--particles ordered by type
        !
        i2 = 0
        over_types: do itype=1,ntypes
           i1 = i2 + 1
           i2 = i2 + npartoftype(itype)
           i2 = min(i2,ninterp)
           if (i1 > i2) cycle over_types
           !--set weights to zero for particle types not used in the rendering
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(i1:i2) = weight_sink
              else
                 weighti(i1:i2) = 0.
              endif
           elseif (idensityweighted) then
           !--for density weighted interpolation use m/h**ndim
              where(dati(i1:i2,ih) > tiny(dati))
                 weighti(i1:i2) = (dati(i1:i2,ipmass)*dunitspmass)/ &
                                  ((dati(i1:i2,ih)*dunitsh)**ndim)
              elsewhere
                 weighti(i1:i2) = 0.
              endwhere
           else
           !--usual interpolation use m/(rho h**ndim)
              where(dati(i1:i2,irho) > tiny(dati) .and. dati(i1:i2,ih) > tiny(dati))
                 weighti(i1:i2) = (dati(i1:i2,ipmass)*dunitspmass)/ &
                                  ((dati(i1:i2,irho)*dunitsrho)*(dati(i1:i2,ih)*dunitsh)**ndim)
              elsewhere
                 weighti(i1:i2) = 0.
              endwhere
           endif
        enddo over_types
     endif

     if (idensityweighted) then
        print "(a)",' USING DENSITY WEIGHTED INTERPOLATION '
        inormalise = .true.
     endif

  elseif (any(masstype(1:ntypes).gt.0.) .and. &
          irho.gt.0 .and. irho.le.ndataplots .and. &
          ih .gt. 0 .and. ih.le.ndataplots .and. &
          required(irho) .and. required(ih)) then

     if (size(iamtypei) > 1) then
        !
        !--particles with mixed types
        !
        !$omp parallel do default(none) &
        !$omp shared(ninterp,iamtypei,weighti,dati,rendersinks,isinktype) &
        !$omp shared(usetype,idensityweighted,dunitsrho,masstype) &
        !$omp shared(ih,irho,dunitspmass,dunitsh,ndim)    &
        !$omp private(ipart,itype)
        do ipart=1,ninterp
           itype = iamtypei(ipart)
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(ipart) = weight_sink              
              else
                 weighti(ipart) = 0.
              endif
           elseif (idensityweighted) then
              if (dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (masstype(itype)*dunitspmass)/ &
                                 ((dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           else
              if (dati(ipart,irho) > tiny(dati) .and. dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (masstype(itype)*dunitspmass)/ &
                                 ((dati(ipart,irho)*dunitsrho)*(dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           endif
        enddo
        !$omp end parallel do
     else
        !
        !--particles ordered by type
        !
        i2 = 0
        over_types2: do itype=1,ntypes
           i1 = i2 + 1
           i2 = i2 + npartoftype(itype)
           i2 = min(i2,ninterp)
           if (i1 > i2) cycle over_types2
           !--set weights to zero for particle types not used in the rendering
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(i1:i2) = weight_sink
              else
                 weighti(i1:i2) = 0.
              endif
           else
              where(dati(i1:i2,irho) > tiny(dati) .and. dati(i1:i2,ih) > tiny(dati))
                 weighti(i1:i2) = masstype(itype)/ &
                                ((dati(i1:i2,irho)*dunitsrho)*(dati(i1:i2,ih)*dunitsh)**ndim)
              elsewhere
                 weighti(i1:i2) = 0.
              endwhere
           endif
        enddo over_types2
     endif

     if (idensityweighted) then
        print "(a)",' USING DENSITY WEIGHTED INTERPOLATION '
        inormalise = .true.
     endif
  else
     if (required(ih) .and. required(irho) .and. ih.gt.0 .and. irho.gt.0) then
        print "(a)",' WARNING: particle mass not set: using normalised interpolations'
     endif
  !--if particle mass has not been set, then must use normalised interpolations
     weighti(1:ninterp) = 1.0
     inormalise = .true.
  endif

end subroutine set_interpolation_weights

end module interpolation
