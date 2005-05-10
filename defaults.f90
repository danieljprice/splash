!-------------------------------------------------
!
! Module containing subroutines relating to 
! setting/saving default options
!
!-------------------------------------------------
module defaults
 implicit none
 
contains

!!
!! set initial default options
!! these are used if no defaults file is found
!!
subroutine defaults_set
  use exact, only:defaults_set_exact
  use filenames
  use labels
  use limits
  use multiplot
  use settings_limits
  use options_data, only:defaults_set_data
  use settings_data, only:ndim
  use settings_part, only:defaults_set_part
  use settings_page, only:defaults_set_page
  use settings_render, only:defaults_set_render
  use settings_vecplot, only:defaults_set_vecplot
  use settings_xsecrot, only:defaults_set_xsecrotate
  use settings_powerspec, only:defaults_set_powerspec
  use particle_data, only:maxpart,maxstep,maxcol
  implicit none
  integer :: i
!
!--set defaults for submenu options
!
  call defaults_set_data
  call defaults_set_limits
  call defaults_set_page
  call defaults_set_part
  call defaults_set_render
  call defaults_set_xsecrotate
  call defaults_set_vecplot
  call defaults_set_exact
  call defaults_set_powerspec 
!
!--limits (could set them to anything but min & max must be different
!          to enable them to be reset interactively if not set elsewhere)
!
  lim(:,1) = 0.
  lim(:,2) = 1.
  itrans(:) = 0
!
!--multiplot
!
  nyplotmulti = 4           ! number of plots in multiplot
  multiploty(:) = 0
  do i=1,4
     multiploty(i) = ndim+i  ! first plot : y axis
  enddo
  multiplotx(:) = 1          ! first plot : x axis
  irendermulti(:) = 0        ! rendering
  ivecplotmulti(:) = 0       ! vector plot
  x_secmulti(:) = .false.    ! take cross section?
  xsecposmulti(:) = 0.0      ! position of cross section
  iplotcontmulti(:) = .false.
  !
  !--array positions of specific quantities
  !
  ix = 0
  !do i=1,ndim
  !   ix(i) = i       ! coords
  !enddo
  ivx = 0      ! vx
  irho = 0     ! density
  ipr = 0      ! pressure
  iutherm = 0  ! thermal energy
  ih = 0       ! smoothing length
  ipmass = 0   ! particle mass
  ipr = 0      ! pressure
  irad = 0     ! radius
  ipowerspec = 0 ! power spectrum
  !
  !--filenames
  !
  rootname = ' '
  !
  !--data array sizes
  !
  maxpart = 0
  maxcol = 0
  maxstep = 0
  !
  !--labels
  !
  do i=1,maxplot
     write(label(i),"(a,i3)") 'column ',i
  enddo
  labeltype(1) = 'gas'
  labeltype(2) = 'type 2'
  labeltype(3) = 'type 3'
  labeltype(4) = 'type 4'
  labeltype(5) = 'type 5'
  labeltype(6) = 'type 6'
  iamvec(:) = 0
  labelvec = ' '
  
  return    
end subroutine defaults_set
!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write
 use exact, only:exactopts,exactparams
 use filenames, only:rootname,nfiles
 use settings_data, only:dataopts
 use settings_part, only:plotopts
 use settings_page, only:pageopts
 use settings_render, only:renderopts
 use settings_vecplot, only:vectoropts
 use settings_xsecrot, only:xsecrotopts
 use settings_powerspec, only:powerspecopts
 use multiplot, only:multi
 implicit none
 integer :: i
       
 open(unit=1,file='defaults',status='replace',form='formatted', &
      delim='apostrophe') ! without delim namelists may not be readable
    write(1,NML=dataopts)
    write(1,NML=plotopts)
    write(1,NML=pageopts)
    write(1,NML=renderopts)
    write(1,NML=vectoropts)
    write(1,NML=xsecrotopts)
    write(1,NML=powerspecopts)
    write(1,NML=exactopts)
    write(1,NML=exactparams)
    write(1,NML=multi)
    do i=1,nfiles
       write(1,"(a)") trim(rootname(i))
    enddo
 close(unit=1)
 print*,'default options saved to file'
    
 return              
end subroutine defaults_write
!-----------------------------------------------
! reads default options from file
! uses namelist input to group the options
! these are specified in the modules
!-----------------------------------------------
subroutine defaults_read
 use filenames, only:rootname,maxfile
 use multiplot
 use settings_data, only:dataopts
 use settings_part, only:plotopts
 use settings_page, only:pageopts
 use settings_render, only:renderopts
 use settings_vecplot, only:vectoropts
 use settings_xsecrot, only:xsecrotopts
 use settings_powerspec, only:powerspecopts
 use exact, only:exactopts,exactparams
 implicit none
 logical :: iexist
 integer :: ierr,i
 
 inquire (exist=iexist, file='defaults')
 if (iexist) then
    open(unit=1,file='defaults',status='old',form='formatted')
    
    ierr = 0
    read(1,NML=dataopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading data options from defaults'    
    
    ierr = 0
    read(1,NML=plotopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading plot options from defaults'

    ierr = 0
    read(1,NML=pageopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading page options from defaults'

    ierr = 0
    read(1,NML=renderopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading render options from defaults'

    ierr = 0
    read(1,NML=vectoropts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading vector plot options from defaults'

    ierr = 0
    read(1,NML=xsecrotopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading xsec/rotation options from defaults'

    ierr = 0
    read(1,NML=powerspecopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading power spectrum options from defaults'

    ierr = 0
    read(1,NML=exactopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading exact solution options from defaults'    

    ierr = 0
    read(1,NML=exactparams,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading exact solution parameters from defaults'    
  
    ierr = 0
    read(1,NML=multi,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading multiplot options from defaults'

    do i=1,maxfile
       read(1,*,end=66,iostat=ierr) rootname(i)
    enddo
66  continue

    close(unit=1)
    print*,'read default options from file '
    return
 else
    print*,'defaults file not found: using program settings'
    return
 endif
 
77 continue
 print*,'**** warning: end of file in defaults ****'
 close(unit=1)

 return
end subroutine defaults_read

end module defaults
