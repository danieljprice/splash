!
! reads default options from file
!
subroutine defaults_read
 use multiplot
 use settings
 use exact_params
 implicit none
 
 inquire (exist=iexist, file='defaults')
 if (iexist) then
    open(unit=1,file='defaults',status='old',form='formatted')
    read(1,NML=plotopts,end=7,err=8)
    read(1,NML=exactparams,end=7,err=8)
    read(1,NML=multi,end=7,err=8)
 close(unit=1)
 print*,'read default options from file '

 endif
 goto 9
7 continue
 print*,'**** warning: end of file in defaults ****'
 close(unit=1)
 goto 9
8 continue
 print*,'error reading defaults from file'
 close(unit=1)
9 continue

 return
end subroutine defaults_read
