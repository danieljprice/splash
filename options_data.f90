!----------------------------------------------------------------------
! sets options relating to current data
! (read new data or change timesteps plotted)
!----------------------------------------------------------------------
subroutine options_data
 use filenames
 use prompting
 use settings
 implicit none
 integer :: ians
  
 print 10, (n_end-nstart+1)/nfreq 
10  format(' 0) exit ',/, 		&
           ' 1) read new data ',/,     	&
           ' 2) change number of timesteps read ( ',i2, ' )')
 call prompt('enter option',ians,0,2)
!
!--options
!
 select case(ians)
 case(1)
    call get_data    
 case(2)
    call prompt('Start at timestep ',nstart,1,nfilesteps)
    call prompt('End at timestep   ',n_end,nstart,nfilesteps)
    call prompt(' Frequency of steps to read',nfreq,1,nfilesteps)
    print *,' Steps = ',(n_end-nstart+1)/nfreq
 end select

 return
end subroutine options_data
