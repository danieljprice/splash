!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write
 use exact_params
 use settings
 use multiplot
 implicit none
       
 open(unit=1,file='defaults',status='replace',form='formatted')
    write(1,NML=plotopts)
    write(1,NML=exactparams)
    write(1,NML=multi)
 close(unit=1)
 print*,'default options saved to file'
    
 return              
end subroutine defaults_write
