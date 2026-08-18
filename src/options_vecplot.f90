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
 integer, parameter :: ivecstyle_streamlines = 0
 integer, parameter :: ivecstyle_arrows = 1
 integer, parameter :: ivecstyle_ironfilings = 2
 integer :: npixvec,minpartforarrow,iVecLegendOnPanel,ivecstyle
 logical :: UseBackgndColorVecplot, iplotpartvec
 logical :: iVecplotLegend,iplotarrowheads
 logical :: iplotsynchrotron,ihidearrowswherenoparts,iallarrowssamelength
 real :: hposlegendvec,vposlegendvec,streamdensity
 real :: rcrit,zcrit,synchrotronspecindex,uthermcutoff

 namelist /vectoropts/ npixvec, UseBackgndColorVecplot,iplotpartvec,&
          iVecplotLegend,hposlegendvec,vposlegendvec, &
          iplotarrowheads,iplotsynchrotron,rcrit,zcrit,synchrotronspecindex, &
          uthermcutoff,ihidearrowswherenoparts,minpartforarrow,iallarrowssamelength,&
          iVecLegendOnPanel,ivecstyle,streamdensity

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_vecplot

 npixvec = 40        ! pixels in x direction on vector plots (arrow style)
 UseBackgndColorVecplot = .false. ! plot vector plot using black/white
 iplotpartvec = .true.   ! whether to plot particles on vector plot
 iVecplotLegend = .true.
 iVecLegendOnPanel = 0   ! all panels
 hposlegendvec = 0.02
 vposlegendvec = 2.0      ! same default as time legend; inside the viewport
 ivecstyle = ivecstyle_streamlines
 streamdensity = 1.0
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
 integer, intent(in) :: ichoose
 integer :: ians

 ians = ichoose
 print "(a)",'--------------- vector plot options -------------------'

 if (ians <= 0 .or. ians > 5) then
    print 10,npixvec,print_logical(UseBackgndColorVecplot), &
             print_logical(iVecplotLegend),ivecstyle, &
             print_logical(ihidearrowswherenoparts)
10  format( &
             ' 0) exit ',/, &
             ' 1) change number of pixels                   (',i4,' )',/, &
             ' 2) use background colour for arrows          ( ',a,' )',/, &
             ' 3) vector plot legend settings               ( ',a,' )',/, &
             ' 4) vector plot style                         (',i4,' )',/, &
             ' 5) hide arrows where there are no particles  ( ',a,' )')
    call prompt('enter option',ians,0,5)
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
    if (ndim==3) then
       call prompt('enter vector plot style (0=streamlines, 1=arrows, 2=iron filings)', &
                   ivecstyle,0,2)
    else
       call prompt('enter vector plot style (0=streamlines, 1=arrows)',ivecstyle,0,1)
    endif
    if (ivecstyle == ivecstyle_streamlines) then
       call prompt('enter streamline density',streamdensity,0.05,10.)
       call prompt('plot arrow heads on streamlines? ',iplotarrowheads)
    elseif (ivecstyle == ivecstyle_arrows) then
       call prompt('plot arrow heads? ',iplotarrowheads)
       if (ndim==3 .and. .not.iplotarrowheads) then
          call prompt(' plot synchrotron map? ',iplotsynchrotron)
          if (iplotsynchrotron) then
             if (iutherm < 0 .or. iutherm > numplot) then
                print "(a)",' Warning: cannot use thermal energy cutoff in synchrotron plots'
                print "(a)",' (could not locate thermal energy in data columns)'
             endif
             call prompt(' enter rcrit for cosmic ray electron distribution exp(-r/rcrit -z/zcrit)',rcrit,0.)
             call prompt(' enter zcrit for cosmic ray electron distribution exp(-r/rcrit -z/zcrit)',zcrit,0.)
             call prompt(' enter synchrotron spectral index I_nu = nu^-alpha ',synchrotronspecindex,0.)
             if (iutherm > 0 .and. iutherm <= numplot) then
                !--set sensible default value for uthermcutoff
                if (uthermcutoff < -tiny(uthermcutoff)) then
                   uthermcutoff = 0.5*(lim(iutherm,1) + lim(iutherm,2))
                endif
                call prompt(' enter threshold thermal energy in current units (u < utherm not used) ',uthermcutoff,0.)
             endif
          endif
       endif
       call prompt('make all arrows same length (ie. only show direction, not magnitude) ?',iallarrowssamelength)
    endif
!------------------------------------------------------------------------
 case(5)
    call prompt('hide vector arrows where there are no particles ? ',ihidearrowswherenoparts)
    if (ihidearrowswherenoparts) then
       call prompt(' enter minimum number of particles in pixel cell for arrow to be plotted ',minpartforarrow,1)
    endif
 end select

 return
end subroutine submenu_vecplot

end module settings_vecplot
