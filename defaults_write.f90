!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write
 use exact
 use filenames
 use settings_data
 use settings_part
 use settings_page
 use settings_render
 use settings_vecplot
 use settings_xsecrot
 use settings_powerspec
 use multiplot
 implicit none
 integer :: i
       
 open(unit=1,file='defaults',status='replace',form='formatted')
    write(1,NML=dataopts)
    write(1,NML=plotopts)
    write(1,NML=pageopts)
    write(1,NML=renderopts)
    write(1,NML=vectoropts)
    write(1,NML=xsecrotopts)
    write(1,NML=powerspecopts)
    write(1,NML=exactparams)
    write(1,NML=multi)
    do i=1,nfiles
       write(1,"(a)") trim(rootname(i))
    enddo
 close(unit=1)
 print*,'default options saved to file'
    
 return              
end subroutine defaults_write
