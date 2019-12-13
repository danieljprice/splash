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
 implicit none
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
 iUseStepList = .false.
 do i=1,size(isteplist)
    isteplist(i) = i
 enddo
 iCalcQuantities = .false.
 DataIsBuffered = .false.
 iRescale = .false.
 ivegotdata = .false.
 ntypes = 1
 xorigin = 0.
 itracktype = 0 ! particle tracking limits (none)
 itrackoffset = 0
 ipartialread = .false.  ! strictly unnecessary as set in get_data
 iverbose = 1

 return
end subroutine defaults_set_data

!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine submenu_data(ichoose)
 use filenames,      only:nsteps,nstepsinfile,ifileopen,unitsfile
 use prompting,      only:prompt,print_logical
 use getdata,        only:get_data,get_labels
 use settings_data,  only:istartatstep,iendatstep,nfreq,iUseStepList, &
                          isteplist,buffer_data,iCalcQuantities,iRescale, &
                          DataIsBuffered,numplot,ncalc,ncolumns
 use calcquantities, only:calc_quantities,setup_calculated_quantities
 use limits,         only:set_limits
 use labels,         only:label,unitslabel,labelzintegration,lenlabel,shortstring
 use settings_units, only:units,set_units,write_unitsfile,unitzintegration
 use fparser,        only:rn,mu0
 implicit none
 integer, intent(in) :: ichoose
 integer             :: ians, i, ncalcwas
 character(len=30)   :: fmtstring
 character(len=1)    :: charp
 character(len=lenlabel) :: labeli
 logical             :: ireadnow,UnitsHaveChanged,iRescaleprev,iwriteunitsfile

 ians = ichoose

 print "(a)",'----------------- data read options -------------------'

 if (ians <= 0 .or. ians > 8) then
    if (iUseStepList) then
       print 10, iendatstep,print_logical(iUseStepList),print_logical(buffer_data), &
                 print_logical(iCalcQuantities),print_logical(iRescale)
    else
       print 10, (iendatstep-istartatstep+1)/nfreq,print_logical(iUseStepList), &
                 print_logical(buffer_data),print_logical(iCalcQuantities), &
                 print_logical(iRescale)
    endif
10  format( &
           ' 0) exit ',/,               &
           ' 1) read new data /re-read data',/,      &
           ' 2) change number of timesteps used        ( ',i5, ' )',/, &
           ' 3) plot selected steps only               (  ',a,' )',/, &
           ' 4) buffer snapshots into memory           (  ',a,' )',/, &
           ' 5) calculate extra quantities             (  ',a,' )',/, &
           ' 6) use physical units                     (  ',a,' )')
    call prompt('enter option',ians,0,6)
 endif
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
 case(1)
    if (buffer_data) then
       call get_data(-1,.false.)
    else
       call get_data(1,.false.,firsttime=.true.)
    endif
!------------------------------------------------------------------------
 case(2)
    iUseStepList = .false.
    call prompt('Start at timestep ',istartatstep,1,nsteps)
    call prompt('End at timestep   ',iendatstep,istartatstep,nsteps)
    call prompt(' Frequency of steps to read',nfreq,1,nsteps)
    print *,' Steps = ',(iendatstep-istartatstep+1)/nfreq
!------------------------------------------------------------------------
 case(3)
    iUseStepList = .true.
    istartatstep = 1
    nfreq = 1
    iendatstep = min(iendatstep,size(isteplist),nsteps)
    call prompt('Enter number of steps to plot ', &
         iendatstep,1,size(isteplist))
    do i=1,iendatstep
       if (isteplist(i) <= 0 .or. isteplist(i) > nsteps) isteplist(i) = i
       write(fmtstring,"(a,i2)") 'Enter step ',i
       call prompt(fmtstring,isteplist(i),1,nsteps)
    enddo
!------------------------------------------------------------------------
 case(4)
    buffer_data = .not.buffer_data
    print "(/a)",' Buffering of data = '//print_logical(buffer_data)
    if (buffer_data) then
       call prompt('Do you want to read all files into memory now?',ireadnow)
       if (ireadnow) then
          call get_data(-1,.true.)
       endif
    endif
!------------------------------------------------------------------------
 case(5)
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
 case(6)
    print "(/,a)",' Current settings for physical units:'
    call get_labels ! reset labels for printing
    print "(2x,a,a3,a,a3,es9.2)",'dz '//trim(labelzintegration),' = ','dz',' x ',unitzintegration
    print "('  0) ',a,a3,a,a3,es9.2)",'time'//trim(unitslabel(0)),' = ','time',' x ',units(0)
    do i=1,ncolumns
       labeli = shortstring(label(i),unitslabel(i))
       print "(i3,') ',a,a3,a,a3,es10.3)",i,trim(labeli)//trim(unitslabel(i)),' = ',trim(labeli),' x ',units(i)
    enddo
    print "(a)"

    UnitsHaveChanged = .false.
    iRescaleprev = iRescale
    charp='y'
    if (iRescale) charp = 'n'
    call prompt('Use physical units? y)es, n)o, e)dit',charp,&
                list=(/'y','Y','n','N','e','E'/),noblank=.true.)
    select case(charp(1:1))
    case('y','Y')
       iRescale = .true.
    case('e','E')
       call set_units(ncolumns,numplot,UnitsHaveChanged)
       iwriteunitsfile = UnitsHaveChanged
       call prompt(' save units to file? ',iwriteunitsfile)
       if (iwriteunitsfile) call write_unitsfile(trim(unitsfile),numplot)
       if (.not.iRescale .and. UnitsHaveChanged) iRescale = .true.
       if (.not.UnitsHaveChanged) call get_labels
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
          call get_data(1,.true.,firsttime=.true.)
       endif
    endif
!------------------------------------------------------------------------
 end select

 return
end subroutine submenu_data

end module options_data
