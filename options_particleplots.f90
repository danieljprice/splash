!-------------------------------------------------------------------------
! Module containing settings and options relating to particle plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_part
 use params
 use settings_data, only:icoordsnew
 implicit none
 integer, dimension(maxparttypes) :: imarktype,idefaultcolourtype
 integer, dimension(100) :: icircpart
 integer :: ncircpart
 integer :: linestyle, linecolour,linestylethisstep,linecolourthisstep, iexact
 logical, dimension(maxparttypes) :: iplotpartoftype,PlotOnRenderings
 logical :: iplotline,ilabelpart,ifastparticleplot

 namelist /plotopts/ iplotline,linestyle,linecolour, &
   imarktype,iplotpartoftype,PlotOnRenderings, &
   iexact,icoordsnew,ifastparticleplot,idefaultcolourtype

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_part
  use settings_data, only:icoords
  implicit none

  ncircpart = 0
  iplotline = .false.     ! plot line joining the particles
  linestyle = 1           ! PGPLOT line style for above
  linecolour = 1
  linestylethisstep = 1
  linecolourthisstep = 1
  iexact = 0              ! exact solution to plot
  ilabelpart = .false.    ! plot particle numbers
  icoordsnew = icoords
  icircpart(:) = 0
  
  iplotpartoftype(1) = .true. ! whether or not to plot particles of certain types
  iplotpartoftype(2:maxparttypes) = .false.
  PlotOnRenderings = .false.
  imarktype = 1              ! PGPLOT marker for all particles
  imarktype(2) = 4           ! PGPLOT marker for ghost/dark matter particles
  imarktype(3) = 17          ! PGPLOT marker for sink particles
  idefaultcolourtype = -1     ! default colour for each particle type
  ifastparticleplot = .true. ! allow crowded-field elimination on particle plots

  return
end subroutine defaults_set_part

!----------------------------------------------------------------------
! submenu with options relating to particle plots
!----------------------------------------------------------------------
subroutine submenu_particleplots(ichoose)
  use exact, only:options_exact,submenu_exact
  use labels, only:labeltype
  use limits, only:lim
  use settings_data, only:icoords,ntypes,ndim,UseTypeInRenderings
  use particle_data, only:npartoftype
  use prompting
  use geometry, only:maxcoordsys,labelcoordsys,coord_transform_limits
  implicit none
  integer, intent(in) :: ichoose
  integer :: i,iaction,n,itype,icoordsprev
  character(len=2) :: charntypes
  character(len=20) :: substring1,substring2
  character(len=1000) :: fmtstring

  iaction = ichoose
  
  !--we require some tricks with the format string to print only the actual number of
  !  particle types rather than the whole array
  !
  if (ntypes.gt.100) print*,'WARNING: Internal error: ntypes too large for formatting in particle plot menu'
  if (ntypes.le.0) then
     substring1 = "no types specified"
     substring2 = "not applicable" 
  elseif (ntypes.eq.1) then
     substring1 = "a"
     substring2 = "i2"
  else
     write(charntypes,"(i2)") ntypes-1
     substring1 = charntypes//"(a,',',1x),a"
     substring2 = charntypes//"(i2,',',1x),i2"
  endif
  
  fmtstring="("// &
         "' 0) exit ',/,"// &
         "' 1) turn on/off particles by type       ( ',"//trim(substring1)//",' )',/,"// &
         "' 2) change graph markers for each type  ( ',"//trim(substring2)//",' )',/,"//  &
         "' 3) set colour for each particle type   ( ',"//trim(substring2)//",' )',/,"//  &
         "' 4) plot line joining particles         ( ',a,' ) ',/,"// &
         "' 5) plot smoothing circles              ( ',i3,' ) ',/,"// &
         "' 6) use fast particle plotting          ( ',a,' ) ',/,"// &
         "' 7) change coordinate systems           ( ',i2,' ) ',/,"// &
         "' 8) plot exact solution                 ( ',i2,' ) ',/,"// &
         "' 9) exact solution plot options ')"

  print "(a)",'------------- particle plot options -------------------'
  if (iaction.le.0 .or. iaction.gt.9) then
     print fmtstring,(trim(print_logical(iplotpartoftype(i))),i=1,ntypes), &
              imarktype(1:ntypes),idefaultcolourtype(1:ntypes),print_logical(iplotline), &
              ncircpart,print_logical(ifastparticleplot),icoordsnew,iexact

     call prompt('enter option',iaction,0,9)
  endif
!
  select case(iaction)
!------------------------------------------------------------------------
  case(1)
     !          plot particles by type?
     do itype=1,ntypes
        call prompt('Plot '//trim(labeltype(itype))//' particles?',iplotpartoftype(itype))
        if (iplotpartoftype(itype) .and. itype.gt.1) then
           if (.not.UseTypeInRenderings(itype)) then
              call prompt('>> Plot '//trim(labeltype(itype))//' particles on top of rendered plots?',PlotOnRenderings(itype))
           else
              PlotonRenderings(itype) = .false.
           endif
        elseif (.not.iplotpartoftype(itype)) then
           PlotonRenderings(itype) = .false.
        endif
     enddo
     return           
!------------------------------------------------------------------------
  case(2)
     print*,'(0 Square) (1 .) (2 +) (3 *) (4 o) (5 x) (17 bold circle) (-8 bigger bold circle)'
     do itype=1,ntypes
        call prompt(' Enter PGPLOT marker for '//trim(labeltype(itype)) &
             //' particles:',imarktype(itype),-8,31)
     enddo
     return   
!------------------------------------------------------------------------
  case(3)
     print "(2(a,/),/,4(a,/))", &
           ' Warning: setting a colour for a particle type overrides', &
           '          (at each new timestep) colours set interactively ', &
           ' -1 = retain interactively set colours between timesteps', &
           '  0 = background ',&
           '  1 = foreground ',&
           '  2->10 = various colours (see PGPLOT default colour indices)'
     do itype=1,ntypes
        call prompt(' Enter default colour for '//trim(labeltype(itype)) &
             //' particles:',idefaultcolourtype(itype),-1,14)
     enddo
     return   
!------------------------------------------------------------------------
  case(4)
     call prompt('plot line joining particles?',iplotline)
     if (iplotline) then     
        call prompt('Enter PGPLOT line style to use ',linestyle,0,5)
        call prompt('Enter PGPLOT colour for line ',linecolour,0,15)
     endif
     return 
!!-----------------------------------------------------------------------
!!  case(5)
!     !          label particles with particle numbers
!     ilabelpart=.not.ilabelpart
!     print*,' label particles = ',ilabelpart
!     return           
!------------------------------------------------------------------------
  case(5)
     print*,'Note that circles of interaction can also be set interactively'
     call prompt('Enter number of circles to draw',ncircpart,0,size(icircpart))
     if (ncircpart.gt.0) then
        do n=1,ncircpart
           if (icircpart(n).eq.0) then
              if (n.gt.1) then
                 icircpart(n) = icircpart(n-1)+1
              else
                 icircpart(n) = 1
              endif
           endif
           call prompt('Enter particle number to plot circle around', &
                    icircpart(n),1,maxval(npartoftype(1,:)))
        enddo
     endif
     return           
!------------------------------------------------------------------------
  case(6)
     print "(1x,a,/,a,/)",'Fast particle plotting excludes particles in crowded regions', &
                     ' Turn this option off to always plot every particle'  
     call prompt('Allow fast particle plotting?',ifastparticleplot)
     return 
!------------------------------------------------------------------------
  case(7)
     print 20,icoords
     do i=1,maxcoordsys
        print 30,i,labelcoordsys(i)
     enddo
20   format(' 0) reset (=',i2,')')
30   format(1x,i1,')',1x,a)
     icoordsprev = icoordsnew
     call prompt(' Enter coordinate system to plot in:', &
                 icoordsnew,0,maxcoordsys)
     if (icoordsnew.eq.0) icoordsnew = icoords
     if (icoordsnew.ne.icoordsprev) then
        call coord_transform_limits(lim(1:ndim,1),lim(1:ndim,2), &
                                    icoordsprev,icoordsnew,ndim)
     endif
     return
!------------------------------------------------------------------------
  case(8)
     call submenu_exact(iexact)
     return
!------------------------------------------------------------------------
  case(9)
     call options_exact
     return     
!------------------------------------------------------------------------
  case default
     return

  end select

  return      
end subroutine submenu_particleplots

end module settings_part
