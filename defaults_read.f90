!-----------------------------------------------
! reads default options from file
! uses namelist input to group the options
! these are specified in the modules
!-----------------------------------------------
subroutine defaults_read
 use multiplot
 use settings
 use exact_params
 implicit none
 logical :: iexist
 
 inquire (exist=iexist, file='defaults')
 if (iexist) then
    open(unit=1,file='defaults',status='old',form='formatted')
    read(1,NML=plotopts,end=77,err=8)
    goto 9
8   print*,'error reading plot options from defaults'
9   continue
    read(1,NML=exactparams,end=77,err=10)
    goto 11
10  print*,'error reading exact solution parameters from defaults'    
11  continue    
    read(1,NML=multi,end=77,err=10)
    goto 13
12  print*,'error reading multiplot options from defaults'
13  continue
 close(unit=1)
 print*,'read default options from file '
 return
 endif
77 continue
 print*,'**** warning: end of file in defaults ****'
 close(unit=1)

 return
end subroutine defaults_read
