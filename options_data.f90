!-------------------------------------------------------------------------
! Module containing settings and options related to the data read
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module options_data
 implicit none
! integer :: numplot,ncalc,ncolumns,nextra
! integer :: ndataplots
! integer :: ndim, ndimv 
! integer :: icoords, iformat, ntypes
! integer :: istartatstep,iendatstep,nfreq
! integer, dimension(10) :: isteplist
! logical :: ivegotdata, buffer_data, iUseStepList!
!
! !namelist /dataopts/ buffer_data
 !
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
  icoords = 1       ! co-ordinate system of simulation
  iformat = 0       ! file format
  buffer_data = .false.
  iUseStepList = .false.
  do i=1,size(isteplist)
     isteplist(i) = i
  enddo
  iCalcQuantities = .false.
  DataIsBuffered = .false.
  units(:) = 1.0
  unitslabel(:) = ' '
  iRescale = .false.
  ivegotdata = .false.
  ntypes = 1
  
  return
end subroutine defaults_set_data

!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine submenu_data
 use filenames, only:nsteps,nstepsinfile,ifileopen
 use prompting
 use getdata, only:get_data
 use settings_data, only:istartatstep,iendatstep,nfreq,iUseStepList, &
     isteplist,buffer_data,iCalcQuantities,iRescale,units,unitslabel, &
     DataIsBuffered,numplot,ncalc,ncolumns
 use calcquantities, only:calc_quantities
 use limits, only:set_limits
 use labels, only:label
 implicit none
 integer :: ians, i, icol
 character(len=30) :: fmtstring
 logical :: ireadnow,UnitsHaveChanged,iRescaleprev
 real :: unitsprev,dunits
 
 ians = 0
 
 if (iUseStepList) then
    print 10, iendatstep,print_logical(iUseStepList),print_logical(buffer_data), &
              print_logical(iCalcQuantities),print_logical(iRescale)
 else
    print 10, (iendatstep-istartatstep+1)/nfreq,print_logical(iUseStepList), &
              print_logical(buffer_data),print_logical(iCalcQuantities), &
              print_logical(iRescale)
 endif
10  format('----------------- data read options -------------------',/,&
           ' 0) exit ',/,               &
           ' 1) read new data /re-read data',/,      &
           ' 2) change number of timesteps used        ( ',i5, ' )',/, &
           ' 3) plot selected steps only               (  ',a,' )',/, &
           ' 4) buffering of data on/off               (  ',a,' )',/, &
           ' 5) turn calculate extra quantities on/off (  ',a,' )',/, &
           ' 6) use physical units                     (  ',a,' )',/,&
           ' 7) change physical unit settings ')
 call prompt('enter option',ians,0,7)
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
    print "(/a,L1)",' Buffering of data = ',buffer_data
    if (buffer_data) then
       call prompt('Do you want to read all data into memory now?',ireadnow)
       if (ireadnow) then
          call get_data(-1,.true.)
       endif
    endif
!------------------------------------------------------------------------
 case(5)
    iCalcQuantities = .not.iCalcQuantities
    if (iCalcQuantities) then
       if (DataIsBuffered) then
          call calc_quantities(1,nsteps)
          call set_limits(1,nsteps,numplot-ncalc+1,numplot)
       else
          call calc_quantities(1,nstepsinfile(ifileopen))
          call set_limits(1,nstepsinfile(ifileopen),numplot-ncalc+1,numplot)
       endif
    else
       print "(/a,L1)",' Calculation of extra quantities = ',iCalcQuantities    
    endif
!------------------------------------------------------------------------
 case(6) 
    print "(a)",'current settings for conversion to physical units are:'
    call set_labels ! reset labels for printing
    do i=1,ncolumns
       print "(a,a3,a,a3,1pe10.3)",trim(label(i))//trim(unitslabel(i)),' = ',trim(label(i)),' x ',units(i)
    enddo
    print "(a,a3,a,a3,1pe8.2)",'time'//trim(unitslabel(0)),' = ','time',' x ',units(0)
    
    iRescaleprev = iRescale
    iRescale = .not.iRescale
    call prompt('Use physical units?',iRescale)
    
    if ((iRescale .and..not. iRescaleprev) .or. (iRescaleprev .and..not.iRescale)) then
       if (buffer_data) then
          call get_data(-1,.true.)
       else
          call get_data(1,.true.,firsttime=.true.)
       endif    
    elseif (iRescale) then
       do i=1,ncolumns
          label(i) = trim(label(i))//trim(unitslabel(i))
       enddo
    endif
!------------------------------------------------------------------------
 case(7)
    UnitsHaveChanged = .false.
    icol = 1
    print "(a)",' *** WARNING: if units are set in the data read, changes here have no effect ***'
    do while(icol.ge.0)
       icol = -1
       call prompt('enter column to change units (0=time,-1=quit,-2=reset all)',icol,-2,numplot)
       if (icol.ge.0) then
          unitsprev = units(icol)
          if (icol.gt.ncolumns) then
             print "(a)",' WARNING: calculated quantities are automatically calculated in physical units '
             print "(a)",' this means that units set here will be re-scalings of these physical values'
          endif
          if (icol.eq.0) then
             call prompt('enter time units (new=old*units)',units(icol),0.)
          else
             call prompt('enter '//trim(label(icol))//' units (new=old*units)',units(icol),0.)
          endif
          if (units(icol).gt.tiny(units)) then
             if (abs(units(icol) - unitsprev).gt.tiny(units)) UnitsHaveChanged = .true.
             if (len_trim(unitslabel(icol)).eq.0 .or. UnitsHaveChanged) then
             !--suggest a label amendment if none already set or if units have changed
                dunits = 1./units(icol)
                if (dunits.gt.100 .or. dunits.lt.1.e-1) then
                   write(unitslabel(icol),"(1pe8.1)") dunits
                else
                   write(unitslabel(icol),"(f5.1)") dunits                  
                endif
                unitslabel(icol) = ' [ x '//trim(adjustl(unitslabel(icol)))//' ]'
             endif
             !--label amendment can be overwritten
             call prompt('enter label amendment ',unitslabel(icol))
          else
             UnitsHaveChanged = .true.
             units(icol) = 1.0
             unitslabel(icol) = ' '
          endif
       elseif (icol.eq.-2) then
          UnitsHaveChanged = .true.
          print "(/a)",' resetting all units to unity...'
          units = 1.0
          unitslabel = ' '
       endif
       print*
    enddo
    
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
       call set_labels
       do i=1,numplot
          label(i) = trim(label(i))//trim(unitslabel(i))
       enddo
    endif
!------------------------------------------------------------------------
 end select

 return
end subroutine submenu_data

end module options_data
