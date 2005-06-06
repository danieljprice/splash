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
  iCalcQuantities = .true.
  DataIsBuffered = .false.
  
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
 implicit none
 integer :: ians, i
 character(len=30) :: fmtstring
 logical :: ireadnow
 
 ians = 0
 if (iUseStepList) then
    print 10, n_end,iUseStepList,buffer_data,iCalcQuantities
 else
    print 10, (n_end-nstart+1)/nfreq,iUseStepList,buffer_data,iCalcQuantities
 endif
 
10  format(' 0) exit ',/,               &
           ' 1) read new data ',/,      &
           ' 2) change number of timesteps used        ( ',i5, ' )',/, &
           ' 3) plot selected steps only               (  ',L1,' )',/, &
           ' 4) buffering of data on/off               (  ',L1, ' )',/, &
           ' 5) turn calculate extra quantities on/off (  ',L1, ' )')
 call prompt('enter option',ians,0,5)
!
!--options
!
 select case(ians)
 case(1)
    if (buffer_data) then
       call get_data(-1,.false.)
    else
       call get_data(1,.false.,firsttime=.true.)
    endif
 case(2)
    iUseStepList = .false.
    call prompt('Start at timestep ',nstart,1,nstepstotal)
    call prompt('End at timestep   ',n_end,nstart,nstepstotal)
    call prompt(' Frequency of steps to read',nfreq,1,nstepstotal)
    print *,' Steps = ',(n_end-nstart+1)/nfreq
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
 case(4)
    buffer_data = .not.buffer_data
    print "(/a,L1)",' Buffering of data = ',buffer_data
    if (buffer_data) then
       call prompt('Do you want to read all data into memory now?',ireadnow)
       if (ireadnow) then
          call get_data(-1,.false.)
       endif
    endif
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
 end select

 return
end subroutine submenu_data

end module options_data
