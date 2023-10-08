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
!  Copyright (C) 2005-2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module containing settings and options related to the data read
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module options_data
 implicit none
 public :: submenu_data,defaults_set_data
 private

contains

!---------------------------------------------
! set default values for these options
! (most should be set upon call to read_data)
!---------------------------------------------
subroutine defaults_set_data
 use settings_data
 use params, only:maxplot
 integer :: i

 numplot=maxplot   ! reset if read from file
 ncalc = 0         ! number of columns to calculate(e.g. radius)
 nextra = 0        ! extra plots aside from particle data
 ncolumns=maxplot-ncalc        ! number of columns in data file
 ndataplots = ncolumns
 ndim = 0          ! number of coordinate dimensions
 ndimV = ndim      ! default velocity same dim as coords
 istartatstep = 1        ! timestep to start from
 iendatstep = 1000      ! timestep to finish on
 nfreq = 1         ! frequency of timesteps to read
 icoords = 1       ! coordinate system of simulation
 iformat = 0       ! file format
 buffer_data = .false.
 buffer_steps_in_file = .false.
 iUseStepList = .false.
 do i=1,size(isteplist)
    isteplist(i) = i
 enddo
 iCalcQuantities = .false.
 DataIsBuffered = .false.
 iRescale = .false.
 enforce_code_units = .false.
 idefaults_file_read = .false.
 ivegotdata = .false.
 ntypes = 1
 xorigin = 0.
 track_string = '0' ! particle tracking limits (none)
 ipartialread = .false.  ! strictly unnecessary as set in get_data
 iverbose = 1
 UseFakeDustParticles = .false.
 iautorender = 0

end subroutine defaults_set_data

!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine submenu_data(ichoose)
 use filenames,      only:nsteps,nstepsinfile,ifileopen
 use prompting,      only:prompt,print_logical
 use getdata,        only:get_data
 use settings_data,  only:buffer_data,iCalcQuantities,iRescale, &
                          DataIsBuffered,numplot,ncalc,ncolumns,UseFakeDustParticles
 use calcquantities, only:calc_quantities,setup_calculated_quantities
 use limits,         only:set_limits,rescale_limits
 use labels,         only:label,unitslabel,labelzintegration,lenlabel,strip_units,idustfrac
 use settings_units, only:units,units_old,set_units,write_unitsfile,unitzintegration
 use fparser,        only:rn,mu0
 integer, intent(in) :: ichoose
 integer             :: ians,i,ncalcwas,maxopt,len_max
 character(len=1)    :: charp
 character(len=lenlabel) :: labeli,labelui
 logical             :: UnitsHaveChanged,iRescaleprev,oldval

 ians = ichoose

 print "(a)",'----------------- data read options -------------------'

 if (ians <= 0 .or. ians > 3) then
    print 10, print_logical(iCalcQuantities), print_logical(iRescale)
10  format( &
           ' 0) exit ',/,               &
           ' 1) calculate extra quantities             (  ',a,' )',/, &
           ' 2) physical units on/off/edit             (  ',a,' )')
    if (idustfrac > 0) then
        print &
         "(' 3) use fake dust particles                (  ',a,' )')",print_logical(UseFakeDustParticles)
       maxopt = 3
    else
       maxopt = 2
    endif
    call prompt('enter option',ians,0,maxopt)
 endif
!
!--options
!
 select case(ians)
!  case(2)
!     iUseStepList = .false.
!     call prompt('Start at timestep ',istartatstep,1,nsteps)
!     call prompt('End at timestep   ',iendatstep,istartatstep,nsteps)
!     call prompt(' Frequency of steps to read',nfreq,1,nsteps)
!     print *,' Steps = ',(iendatstep-istartatstep+1)/nfreq
! !------------------------------------------------------------------------
!  case(3)
!     iUseStepList = .true.
!     istartatstep = 1
!     nfreq = 1
!     iendatstep = min(iendatstep,size(isteplist),nsteps)
!     call prompt('Enter number of steps to plot ', &
!          iendatstep,1,size(isteplist))
!     do i=1,iendatstep
!        if (isteplist(i) <= 0 .or. isteplist(i) > nsteps) isteplist(i) = i
!        write(fmtstring,"(a,i2)") 'Enter step ',i
!        call prompt(fmtstring,isteplist(i),1,nsteps)
!     enddo
 case(1)
    ncalcwas = ncalc
    call setup_calculated_quantities(ncalc)

    iCalcQuantities = (ncalc > 0)
    if (iCalcQuantities) then
       if (DataIsBuffered) then
          call calc_quantities(1,nsteps)
          call set_limits(1,nsteps,ncolumns+ncalcwas+1,ncolumns+ncalc)
       else
          if (ifileopen > 0) then
             call calc_quantities(1,nstepsinfile(ifileopen))
             call set_limits(1,nstepsinfile(ifileopen),ncolumns+ncalcwas+1,ncolumns+ncalc)
          endif
       endif
    else
       print "(/a)",' Calculation of extra quantities is '//print_logical(iCalcQuantities)
    endif
!------------------------------------------------------------------------
 case(2)
    print "(/,a,/)",' Current settings for physical units:'
    len_max = min(maxval(len_trim(label(:))),20)
    labeli = 'dz '//trim(labelzintegration)
    print "(2x,a,es10.3,a)",labeli(1:len_max+3)//' = ',unitzintegration,' x dz'
    labeli = 'time'//trim(unitslabel(0))
    print "('  0) ',a,es10.3,a)",labeli(1:len_max)//' = ',units(0),' x time'
    do i=1,ncolumns
       labeli = strip_units(label(i),unitslabel(i))
       labelui = trim(labeli)//trim(unitslabel(i))
       print "(i3,') ',a,es10.3,a)",i,labelui(1:len_max)//' = ',units(i),' x '//trim(labeli)
    enddo
    print "(a)"

    UnitsHaveChanged = .false.
    iRescaleprev = iRescale
    units_old = 1.d0
    if (iRescale) units_old = units
    charp='y'
    if (iRescale) charp = 'n'
    call prompt('Use physical units? y)es, n)o, e)dit',charp,&
                list=(/'y','Y','n','N','e','E'/),noblank=.true.)
    select case(charp(1:1))
    case('y','Y')
       iRescale = .true.
    case('e','E')
       call set_units(ncolumns,numplot,UnitsHaveChanged)
       if (.not.iRescale .and. UnitsHaveChanged) iRescale = .true.
    case default
       iRescale = .false.
    end select

    if ((iRescale .and..not. iRescaleprev) .or. (iRescaleprev .and..not.iRescale) &
        .or. UnitsHaveChanged) then
       if (iRescale) then
          mu0 = 12.566370614359_rn
       else
          mu0 = 1.0_rn
       endif
       if (buffer_data) then
          call get_data(-1,.true.)
       else
          call get_data(1,.true.,firsttime=.false.)
       endif
       if (iRescale) then
          call rescale_limits(fac=units(1:)/units_old(1:))
       else ! switching from physical back to code units
          units_old = 1.d0
          call rescale_limits(fac=units_old(1:)/units(1:))
       endif
    endif
 case(3)
    oldval = UseFakeDustParticles
    call prompt( 'Use fake dust particles?',UseFakeDustParticles)
    ! re-read data if option has changed
    if (UseFakeDustParticles .neqv. oldval) then
       if (buffer_data) then
          call get_data(-1,.true.)
       else
          call get_data(1,.true.,firsttime=.true.)
       endif
    endif
 end select

end subroutine submenu_data

end module options_data
