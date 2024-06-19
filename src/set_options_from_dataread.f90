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
module set_options_from_dataread
 implicit none

contains

subroutine set_options_dataread()
 use exact,          only:read_exactparams
 use settings_part,  only:iplotpartoftype,PlotonRenderings
 use particle_data,  only:npartoftype
 use settings_data,  only:ndim,ncolumns,ntypes,iexact,iverbose,UseTypeInRenderings,idefaults_file_read
 use filenames,      only:rootname,nsteps,ifileopen
 use labels,         only:labeltype,irho,ih,get_sink_type
 use asciiutils,     only:ucase
 use settings_render,  only:icolour_particles
 use settings_xsecrot, only:xsec_nomulti
 integer :: itype,nplot,ierr
 !
 !--check for errors in data read / print warnings
 !
 if (ndim > 0 .and. ncolumns > 0 .and. nsteps > 0 .and. .not.idefaults_file_read) then
    if (sum(npartoftype(:,1)) > 0 .and. npartoftype(1,1)==0 .and. .not.any(iplotpartoftype(2:))) then
       if (iverbose >= 1) print "(/,a)",' WARNING! DATA APPEARS TO CONTAIN NO '//trim(ucase(labeltype(1)))//' PARTICLES'
       itype = 0
       nplot = 0
       do while (nplot==0 .and. itype < ntypes)
          itype = itype + 1
          if (npartoftype(itype,1) > 0) then
             iplotpartoftype(itype) = .true.
             if (UseTypeInRenderings(itype)) nplot = nplot + npartoftype(itype,1)
             if (iverbose >= 1) print "(a)",' (plotting of '//trim(labeltype(itype))//' particles turned ON)'
          endif
       enddo
       print*
    endif
    ! Set particle colouring if it can't render
    if (irho==0 .or. ih==0) then
       icolour_particles = .true.
    endif
 endif
 !
 !--turn sink particles ON by default
 !
 if (ndim > 0 .and. .not.idefaults_file_read) then  ! only on first data read
    itype = get_sink_type(ntypes)
    if (itype > 0) then
       if (npartoftype(itype,1) > 0 .and. .not.UseTypeInRenderings(itype)) then
          iplotpartoftype(itype) = .true.
          PlotOnRenderings(itype) = .true.
       endif
    endif
 endif

 !--for fits cubes, turn on cross section by default if 3D
 if (ndim >= 3 .and. .not.idefaults_file_read .and. ifileopen > 0) then  ! only on first data read
    if (index(rootname(ifileopen),'.fits') > 0) then
       print "(a)",' turning cross section ON for fits cube'
       xsec_nomulti = .true.
    endif
 endif
 !
 !--read exact solution parameters from files if present
 !
 if (iexact /= 0) then
    if (ifileopen > 0) then
       call read_exactparams(iexact,rootname(ifileopen),ierr)
    endif
 endif

end subroutine set_options_dataread

end module set_options_from_dataread
