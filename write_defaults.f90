!
!     writes default options to file (should match read_defaults)
!
subroutine write_defaults
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
end subroutine write_defaults
