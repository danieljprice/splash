!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write
 use exact
 use settings
 use multiplot
 implicit none
 integer :: ierr
       
 open(unit=1,file='defaults',status='replace',form='formatted')
    write(1,NML=plotopts)
    write(1,NML=pageopts)
    write(1,NML=renderopts)
    write(1,NML=vectoropts)
    call defaults_write_exact(1,ierr)
    write(1,NML=multi)
 close(unit=1)
 print*,'default options saved to file'
    
 return              
end subroutine defaults_write
