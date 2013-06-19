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
 use labels,         only:label,unitslabel,labelzintegration
 use settings_units, only:units,set_units,write_unitsfile,unitzintegration
 implicit none
 integer, intent(in) :: ichoose
 integer             :: ians, i
 character(len=30)   :: fmtstring
 logical             :: ireadnow,UnitsHaveChanged,iRescaleprev,iwriteunitsfile

 ians = ichoose

 print "(a)",'----------------- data read options -------------------'

 if (ians.le.0 .or. ians.gt.8) then
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
           ' 4) buffering of data on/off               (  ',a,' )',/, &
           ' 5) turn calculate extra quantities on/off (  ',a,' )',/, &
           ' 6) edit list of calculated quantities               ',/, &
           ' 7) use physical units                     (  ',a,' )',/,&
           ' 8) change physical unit settings ')
    call prompt('enter option',ians,0,8)
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
       if (isteplist(i).le.0 .or. isteplist(i).gt.nsteps) isteplist(i) = i
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
 case(5,6)
    if (ians.eq.5) iCalcQuantities = .not.iCalcQuantities

    if (iCalcQuantities .or. ians.eq.6) then
       call setup_calculated_quantities(ncalc)

       if (ians.eq.6 .and. .not.iCalcQuantities) then
          if (ncalc.gt.0) iCalcQuantities = .true.
       endif
       if (iCalcQuantities) then
          if (DataIsBuffered) then
             call calc_quantities(1,nsteps)
             call set_limits(1,nsteps,ncolumns+1,ncolumns+ncalc)
          else
             if (ifileopen.gt.0) then
                call calc_quantities(1,nstepsinfile(ifileopen))
                call set_limits(1,nstepsinfile(ifileopen),ncolumns+1,ncolumns+ncalc)
             endif
          endif
       endif
    else
       print "(/a)",' Calculation of extra quantities is '//print_logical(iCalcQuantities)
    endif
!------------------------------------------------------------------------
 case(7)
    print "(a)",'current settings for conversion to physical units are:'
    call get_labels ! reset labels for printing
    do i=1,ncolumns
       print "(a,a3,a,a3,es10.3)",trim(label(i))//trim(unitslabel(i)),' = ',trim(label(i)),' x ',units(i)
    enddo
    print "(a,a3,a,a3,es9.2)",'time'//trim(unitslabel(0)),' = ','time',' x ',units(0)
    print "(a,a3,a,a3,es9.2)",'dz '//trim(labelzintegration),' = ','dz',' x ',unitzintegration

    iRescaleprev = iRescale
    iRescale = .not.iRescale
    call prompt('Use physical units?',iRescale)

    if ((iRescale .and..not. iRescaleprev) .or. (iRescaleprev .and..not.iRescale)) then
       if (buffer_data) then
          call get_data(-1,.true.)
       else
          call get_data(1,.true.,firsttime=.true.)
       endif
    endif
!------------------------------------------------------------------------
 case(8)
    UnitsHaveChanged = .false.
    call set_units(ncolumns,numplot,UnitsHaveChanged)

    iwriteunitsfile = .true.
    call prompt(' save units to file? ',iwriteunitsfile)
    if (iwriteunitsfile) call write_unitsfile(trim(unitsfile),numplot)

    if (.not.iRescale .and. UnitsHaveChanged) call prompt('Apply physical units to data?',iRescale)

    !
    !--re-read/rescale data if units have changed
    !
    if (UnitsHaveChanged) then
       if (buffer_data) then
          call get_data(-1,.true.)
       else
          call get_data(1,.true.,firsttime=.true.)
       endif
    else
       call get_labels
    endif
!------------------------------------------------------------------------
 end select

 return
end subroutine submenu_data

end module options_data
