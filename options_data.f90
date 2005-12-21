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
! integer :: nstart,n_end,nfreq
! integer, dimension(10) :: isteplist
! logical :: ivegotdata, buffer_data, iUseStepList!
!
! !namelist /dataopts/ buffer_data
 !
! public :: submenu_data,defaults_set_data
! private

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
  ndim = 0          ! number of coordinate dimensions
  ndimV = ndim      ! default velocity same dim as coords
  nstart = 1        ! timestep to start from
  n_end = 1000      ! timestep to finish on
  nfreq = 1         ! frequency of timesteps to read
  icoords = 1       ! co-ordinate system of simulation
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
  
  return
end subroutine defaults_set_data

!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine submenu_data
 use filenames, only:nstepstotal,nstepsinfile,ifileopen
 use prompting
 use getdata, only:get_data
 use settings_data
 use calcquantities, only:calc_quantities
 use limits, only:set_limits
 use labels, only:label
 implicit none
 integer :: ians, i, icol
 character(len=30) :: fmtstring
 logical :: ireadnow,UnitsHaveChanged,iChange
 real :: unitsprev,dunits
 
 ians = 0
 
 if (iUseStepList) then
    print 10, n_end,iUseStepList,buffer_data,iCalcQuantities,iRescale
 else
    print 10, (n_end-nstart+1)/nfreq,iUseStepList,buffer_data,iCalcQuantities,iRescale
 endif
 
10  format(' 0) exit ',/,               &
           ' 1) read new data /re-read data',/,      &
           ' 2) change number of timesteps used        ( ',i5, ' )',/, &
           ' 3) plot selected steps only               (  ',L1,' )',/, &
           ' 4) buffering of data on/off               (  ',L1, ' )',/, &
           ' 5) turn calculate extra quantities on/off (  ',L1, ' )',/, &
           ' 6) rescale data (change units)            (  ',L1, ' )')
 call prompt('enter option',ians,0,6)
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
    call prompt('Start at timestep ',nstart,1,nstepstotal)
    call prompt('End at timestep   ',n_end,nstart,nstepstotal)
    call prompt(' Frequency of steps to read',nfreq,1,nstepstotal)
    print *,' Steps = ',(n_end-nstart+1)/nfreq
!------------------------------------------------------------------------
 case(3)
    iUseStepList = .true.
    nstart = 1
    nfreq = 1
    n_end = min(n_end,size(isteplist))
    call prompt('Enter number of steps to plot ', &
         n_end,1,size(isteplist))
    do i=1,n_end
       if (isteplist(i).le.0 .or. isteplist(i).gt.nstepstotal) isteplist(i) = i
       write(fmtstring,"(a,i2)") 'Enter step ',i
       call prompt(fmtstring,isteplist(i),1,nstepstotal)
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
          call calc_quantities(nstart,n_end)
          call set_limits(nstart,n_end,numplot-ncalc+1,numplot)
       else
          call calc_quantities(1,nstepsinfile(ifileopen))
          call set_limits(1,nstepsinfile(ifileopen),numplot-ncalc+1,numplot)
       endif
    else
       print "(/a,L1)",' Calculation of extra quantities = ',iCalcQuantities    
    endif
!------------------------------------------------------------------------
 case(6)
    UnitsHaveChanged = .false.
    call prompt('Rescale data using unit arrays?',iRescale)
    if (iRescale) then
       iChange = .false.
       call prompt('Do you want to set the units? '// &
                   '(warning: may be overwritten by data read)',iChange)
       icol = 0
       if (iChange) icol = 1
       do while(icol.gt.0)
          icol = 0
          call prompt('enter column to rescale (0=quit)(-1=reset all)',icol,-1,numplot)
          if (icol.gt.0) then
             unitsprev = units(icol)          
             call prompt('enter units for this column (new=old*units)',units(icol),0.)
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
          elseif (icol.lt.0) then
             UnitsHaveChanged = .true.
             print "(/a)",' resetting all units to unity...'
             units = 1.0
             unitslabel = ' '
          endif
          print*
       enddo
    endif
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
