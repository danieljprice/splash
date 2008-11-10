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
!! initialise some variables
!! this should only be called at code start
!!
subroutine defaults_set_initial
  use filenames, only:rootname
  use labels, only:label,labeltype,iamvec,labelvec,reset_columnids
  use limits, only:lim,range
  use particle_data, only:maxpart,maxstep,maxcol
  use settings_data, only:UseTypeInRenderings
  use settings_page, only:device
  implicit none
  integer :: i
!
!--limits (could set them to anything but min & max must be different
!          to enable them to be reset interactively if not set elsewhere)
!
  lim(:,1) = 0.
  lim(:,2) = 1.
  range(:,:) = 0.
  call reset_columnids()
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
  !  column labels
  do i=1,size(label)
     write(label(i),"(a,i3)") 'column ',i
  enddo
  !  particle types
  labeltype(1) = 'gas'
  do i=2,size(labeltype)
     write(labeltype(i),"(a,i1)") 'type ',i
  enddo
  UseTypeInRenderings(:) = .false.
  UseTypeInRenderings(1) = .true.  
  
  !  vector labels
  iamvec(:) = 0
  labelvec = ' '
  
  !  device from command line
  device = ' '

  return
end subroutine defaults_set_initial

!!
!! set initial default options
!! these are used if no defaults file is found
!!
subroutine defaults_set(use_evdefaults)
  use exact, only:defaults_set_exact
  use multiplot
  use settings_limits, only:defaults_set_limits
  use options_data, only:defaults_set_data
  use settings_part, only:defaults_set_part,defaults_set_part_ev
  use settings_page, only:defaults_set_page,defaults_set_page_ev
  use settings_render, only:defaults_set_render
  use settings_vecplot, only:defaults_set_vecplot
  use settings_xsecrot, only:defaults_set_xsecrotate
  use settings_powerspec, only:defaults_set_powerspec
  use settings_units, only:defaults_set_units
  use titles, only:pagetitles,steplegend
  implicit none
  logical, intent(in) :: use_evdefaults
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
  call defaults_set_units
!
!--if using evsplash, override some default options
!
  if (use_evdefaults) then
     print "(a)",'setting evsplash defaults'
     call defaults_set_page_ev
     call defaults_set_part_ev
  endif

  itrans(:) = 0
!
!--multiplot
!
  nyplotmulti = 4           ! number of plots in multiplot
  multiploty(:) = 0
  do i=1,4
     multiploty(i) = 1+i       ! first plot : y axis
  enddo
  multiplotx(:) = 1          ! first plot : x axis
  irendermulti(:) = 0        ! rendering
  ivecplotmulti(:) = 0       ! vector plot
  x_secmulti(:) = .false.    ! take cross section?
  xsecposmulti(:) = 0.0      ! position of cross section
  iplotcontmulti(:) = .false.
  !
  !--titles
  !
  pagetitles = ' '
  steplegend = ' '
  
  return    
end subroutine defaults_set
!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write(filename)
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
 use shapes, only:shapeopts
 implicit none
 character(len=*), intent(in) :: filename
 integer :: i,ierr
       
 open(unit=1,file=trim(adjustl(filename)),status='replace',form='formatted', &
      delim='apostrophe',iostat=ierr) ! without delim namelists may not be readable
    if (ierr /= 0) then 
       print*,'ERROR: cannot write file '//trim(filename)
       close(unit=1)
       return
    endif
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
    write(1,NML=shapeopts)
    do i=1,nfiles
       write(1,"(a)") trim(rootname(i))
    enddo
 close(unit=1)
 print "(a)",'default options saved to file '//trim(filename)
    
 return              
end subroutine defaults_write
!-----------------------------------------------
! reads default options from file
! uses namelist input to group the options
! these are specified in the modules
!-----------------------------------------------
subroutine defaults_read(filename)
 use filenames, only:rootname,maxfile
 use multiplot, only:multi
 use settings_data, only:dataopts
 use settings_part, only:plotopts
 use settings_page, only:pageopts
 use settings_render, only:renderopts
 use settings_vecplot, only:vectoropts
 use settings_xsecrot, only:xsecrotopts
 use settings_powerspec, only:powerspecopts
 use exact, only:exactopts,exactparams
 use shapes, only:shapeopts
 implicit none
 character(len=*), intent(in) :: filename
 logical :: iexist
 integer :: ierr,i
 
 inquire (exist=iexist, file=filename)
 if (iexist) then
    open(unit=1,file=filename,status='old',form='formatted',delim='apostrophe',err=88)
    
    ierr = 0
    read(1,NML=dataopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading data options from '//trim(filename)    
    
    ierr = 0
    read(1,NML=plotopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading plot options from '//trim(filename)

    ierr = 0
    read(1,NML=pageopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading page options from '//trim(filename)

    ierr = 0
    read(1,NML=renderopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading render options from '//trim(filename)

    ierr = 0
    read(1,NML=vectoropts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading vector plot options from '//trim(filename)

    ierr = 0
    read(1,NML=xsecrotopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading xsec/rotation options from '//trim(filename)

    ierr = 0
    read(1,NML=powerspecopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading power spectrum options from '//trim(filename)

    ierr = 0
    read(1,NML=exactopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading exact solution options from '//trim(filename)    

    ierr = 0
    read(1,NML=exactparams,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading exact solution parameters from '//trim(filename)    
  
    ierr = 0
    read(1,NML=multi,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading multiplot options from '//trim(filename)

    ierr = 0
    read(1,NML=shapeopts,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading shape options from '//trim(filename)

    if (len_trim(rootname(1)).eq.0) then
       do i=1,maxfile
          read(1,*,end=66,iostat=ierr) rootname(i)
       enddo
    endif
66  continue

    close(unit=1)
    print*,'read default options from '//trim(filename)
    return
 else
    print*,trim(filename)//': file not found: using program settings'
    return
 endif

88 continue
 print "(a)",' *** error opening defaults file '//trim(filename)//': using program settings'

 return
end subroutine defaults_read

end module defaults
