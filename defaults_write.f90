!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write
 use exact
 use filenames
 use settings
 use multiplot
 implicit none
 integer :: ierr, i
       
 open(unit=1,file='defaults',status='replace',form='formatted')
    write(1,NML=plotopts)
    write(1,NML=pageopts)
    write(1,NML=renderopts)
    write(1,NML=vectoropts)
    write(1,NML=exactparams)
    write(1,NML=multi)
    do i=1,nfiles
       write(1,"(a)") trim(rootname(i))
    enddo
 close(unit=1)
 print*,'default options saved to file'
    
 return              
end subroutine defaults_write
