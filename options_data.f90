!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine options_data
 use filenames
 use prompting
 use settings_data
 implicit none
 integer :: ians
 
 ians = 0
 print 10, (n_end-nstart+1)/nfreq,buffer_data
10  format(' 0) exit ',/,               &
           ' 1) read new data ',/,      &
           ' 2) change number of timesteps read ( ',i2, ' )',/, &
           ' 3) toggle buffering of data        (  ',L1, ' )')
 call prompt('enter option',ians,0,3)
!
!--options
!
 select case(ians)
 case(1)
    if (buffer_data) then
       call get_data(-1)
    else
       call get_data(1)
    endif
 case(2)
    call prompt('Start at timestep ',nstart,1,nstepstotal)
    call prompt('End at timestep   ',n_end,nstart,nstepstotal)
    call prompt(' Frequency of steps to read',nfreq,1,nstepstotal)
    print *,' Steps = ',(n_end-nstart+1)/nfreq
 case(3)
    buffer_data = .not.buffer_data
    print*,'buffering of data = ',buffer_data
 end select

 return
end subroutine options_data
