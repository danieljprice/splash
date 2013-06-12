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

!-------------------------------------------------------------------------
! Module containing settings and options relating to vector plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_vecplot
 implicit none
 integer :: npixvec,minpartforarrow,iVecLegendOnPanel
 logical :: UseBackgndColorVecplot, iplotpartvec
 logical :: iVecplotLegend,iplotstreamlines,iplotarrowheads
 logical :: iplotsynchrotron,ihidearrowswherenoparts,iallarrowssamelength
 real :: hposlegendvec,vposlegendvec
 real :: rcrit,zcrit,synchrotronspecindex,uthermcutoff

 namelist /vectoropts/ npixvec, UseBackgndColorVecplot,iplotpartvec,&
          iVecplotLegend,hposlegendvec,vposlegendvec,iplotstreamlines, &
          iplotarrowheads,iplotsynchrotron,rcrit,zcrit,synchrotronspecindex, &
          uthermcutoff,ihidearrowswherenoparts,minpartforarrow,iallarrowssamelength,&
          iVecLegendOnPanel

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_vecplot
  implicit none

  npixvec = 40        ! pixels in x direction on vector plots
  UseBackgndColorVecplot = .false. ! plot vector plot using black/white
  iplotpartvec = .true.   ! whether to plot particles on vector plot
  iVecplotLegend = .true.
  iVecLegendOnPanel = 0   ! all panels
  hposlegendvec = 0.02
  vposlegendvec = -1.5
  iplotstreamlines = .false. ! plot stream lines instead of arrows
  iplotarrowheads = .true.
  iplotsynchrotron = .false.
  zcrit = 2.5 ! kpc
  rcrit = 13. ! kpc
  synchrotronspecindex = 0.8
  uthermcutoff = -1. ! flags this as uninitialised
  ihidearrowswherenoparts = .false.
  minpartforarrow = 1
  iallarrowssamelength = .false.

  return
end subroutine defaults_set_vecplot

!----------------------------------------------------------------------
! sets options relating to vector plots
!----------------------------------------------------------------------
subroutine submenu_vecplot(ichoose)
 use prompting,     only:prompt,print_logical
 use settings_data, only:ndim,numplot
 use labels,        only:iutherm
 use limits,        only:lim
 use legends,       only:prompt_panelselect
 implicit none
 integer, intent(in) :: ichoose
 integer :: ians

 ians = ichoose
 print "(a)",'--------------- vector plot options -------------------'

 if (ians.le.0 .or. ians.gt.7) then
    print 10,npixvec,print_logical(UseBackgndColorVecplot), &
             print_logical(iVecplotLegend),print_logical(iplotstreamlines), &
             print_logical(iplotarrowheads), &
             print_logical(ihidearrowswherenoparts), &
             print_logical(iallarrowssamelength)
10  format( &
             ' 0) exit ',/, &
             ' 1) change number of pixels                   (',i4,' )',/, &
             ' 2) use background colour for arrows          ( ',a,' )',/, &
             ' 3) vector plot legend settings               ( ',a,' )',/, &
             ' 4) plot stream/field lines instead of arrows ( ',a,' )',/, &
             ' 5) turn arrow heads on/off                   ( ',a,' )',/, &
             ' 6) hide arrows where there are no particles  ( ',a,' )',/, &
             ' 7) all arrows same length - ie. direction only ( ',a,' )')
    call prompt('enter option',ians,0,7)
 endif
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
 case(1)
    call prompt('enter number of pixels',npixvec,1,1000)
!------------------------------------------------------------------------
 case(2)
    UseBackgndColorVecplot = .not.UseBackgndColorVecplot
    print*,'use background colour on vector plots = ', &
           print_logical(UseBackgndColorVecplot)
!------------------------------------------------------------------------
 case(3)
    call prompt('plot vector legend?',iVecplotLegend)
    if (iVecplotLegend) then
       print*,'note: H key in interactive mode can also be used to set positions'
       call prompt('Enter horizontal position as fraction of viewport', &
                   hposlegendvec,0.0,1.0)
       call prompt('Enter vertical position in character heights from top', &
                    vposlegendvec)
       call prompt_panelselect('vector legend',iVecLegendOnPanel)
    endif
!------------------------------------------------------------------------
 case(4)
    iplotstreamlines = .not.iplotstreamlines
    print "(2(a,/))",' Note: the number of stream lines plotted is determined by', &
                     ' the "change number of contours" option in the r)ender menu'
    call prompt('use stream lines instead of arrows? ',iplotstreamlines)
!------------------------------------------------------------------------
 case(5)
    iplotarrowheads = .not.iplotarrowheads
    call prompt('plot arrow heads? ',iplotarrowheads)
    if (ndim.eq.3 .and. .not.iplotarrowheads) then
       call prompt(' plot synchrotron map? ',iplotsynchrotron)
       if (iplotsynchrotron) then
          if (iutherm.lt.0 .or. iutherm.gt.numplot) then
             print "(a)",' Warning: cannot use thermal energy cutoff in synchrotron plots'
             print "(a)",' (could not locate thermal energy in data columns)'
          endif
          call prompt(' enter rcrit for cosmic ray electron distribution exp(-r/rcrit -z/zcrit)',rcrit,0.)
          call prompt(' enter zcrit for cosmic ray electron distribution exp(-r/rcrit -z/zcrit)',zcrit,0.)
          call prompt(' enter synchrotron spectral index I_nu = nu^-alpha ',synchrotronspecindex,0.)
          if (iutherm.gt.0 .and. iutherm.le.numplot) then
             !--set sensible default value for uthermcutoff
             if (uthermcutoff.lt.-tiny(uthermcutoff)) then
                uthermcutoff = 0.5*(lim(iutherm,1) + lim(iutherm,2))
             endif
             call prompt(' enter threshold thermal energy in current units (u < utherm not used) ',uthermcutoff,0.)
          endif
       endif
    endif
!------------------------------------------------------------------------
 case(6)
    call prompt('hide vector arrows where there are no particles ? ',ihidearrowswherenoparts)
    if (ihidearrowswherenoparts) then
       call prompt(' enter minimum number of particles in pixel cell for arrow to be plotted ',minpartforarrow,1)
    endif
!------------------------------------------------------------------------
 case(7)
    iallarrowssamelength = .not.iallarrowssamelength
    call prompt('make all arrows same length (ie. only show direction, not magnitude) ?',iallarrowssamelength)
 end select

 return
end subroutine submenu_vecplot

end module settings_vecplot
