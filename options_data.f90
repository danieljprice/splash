!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine options_data
 use filenames, only:nstepstotal
 use prompting
 use settings_data
 use getdata, only:get_data
 implicit none
 integer :: ians, i
 character(len=30) :: fmtstring
 
 ians = 0
 if (iUseStepList) then
    print 10, n_end,iUseStepList,buffer_data
 else
    print 10, (n_end-nstart+1)/nfreq,iUseStepList,buffer_data
 endif
 
10  format(' 0) exit ',/,               &
           ' 1) read new data ',/,      &
           ' 2) change number of timesteps used ( ',i5, ' )',/, &
           ' 3) plot selected steps only        (  ',L1,' )',/, &
           ' 4) toggle buffering of data        (  ',L1, ' )')
 call prompt('enter option',ians,0,4)
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
    print*,'buffering of data = ',buffer_data
 end select

 return
end subroutine options_data
