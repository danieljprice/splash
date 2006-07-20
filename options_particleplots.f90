!-------------------------------------------------------------------------
! Module containing settings and options relating to particle plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_part
 use params
 implicit none
 integer, dimension(maxparttypes) :: imarktype
 integer, dimension(100) :: icircpart
 integer :: ncircpart, icoordsnew
 integer :: linestyle, linecolour,linestylethisstep,linecolourthisstep, iexact
 logical, dimension(maxparttypes) :: iplotpartoftype,PlotOnRenderings
 logical :: iplotline,ilabelpart

 namelist /plotopts/ iplotline,linestyle,linecolour, &
   imarktype,iplotpartoftype,PlotOnRenderings, &
   iexact,icoordsnew

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

  return
end subroutine defaults_set_part

!----------------------------------------------------------------------
! submenu with options relating to particle plots
!----------------------------------------------------------------------
subroutine submenu_particleplots
  use exact, only:options_exact,submenu_exact
  use labels, only:labeltype
  use limits, only:lim
  use settings_data, only:icoords,ntypes,ndim,UseTypeInRenderings
  use particle_data, only:npartoftype
  use prompting
  use geometry, only:maxcoordsys,labelcoordsys,coord_transform_limits
  implicit none
  integer :: i,iaction,n,itype,icoordsprev

  iaction = 0
  print 10, iplotline,ilabelpart,ncircpart, &
        iplotpartoftype,imarktype,icoordsnew,iexact
10  format(' 0) exit ',/,                 &
         ' 1) plot line joining particles        ( ',L1,' ) ',/, &
         ' 2) label particles                    ( ',L1,' ) ',/,           &
         ' 3) plot circles of interaction        ( ',i3,' ) ',/,           &
         ' 4) turn on/off particles by type      ( ',6(L1,',',1x),' )',/,  &
         ' 5) change graph markers for each type ( ',6(i2,',',1x),' )',/,  &
         ' 6) change coordinate systems          ( ',i2,' ) ',/,           &
         ' 7) plot exact solution                ( ',i2,' ) ',/, &
         ' 8) exact solution options')
    call prompt('enter option',iaction,0,8)
!
  select case(iaction)

!------------------------------------------------------------------------
  case(1)
     call prompt('plot line joining particles?',iplotline)
     if (iplotline) then     
        call prompt('Enter PGPLOT line style to use ',linestyle,0,5)
        call prompt('Enter PGPLOT colour for line ',linecolour,0,15)
     endif
     return 
!-----------------------------------------------------------------------
  case(2)
     !          label particles with particle numbers
     ilabelpart=.not.ilabelpart
     print*,' label particles = ',ilabelpart
     return           
!------------------------------------------------------------------------
  case(3)
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
  case(4)
     !          plot particles by type?
     do itype=1,ntypes
        call prompt('Plot '//trim(labeltype(itype))//' particles?',iplotpartoftype(itype))
        if (iplotpartoftype(itype) .and. itype.gt.1) then
           if (.not.UseTypeInRenderings(itype)) then
              call prompt('Plot on top of rendered plots?',PlotOnRenderings(itype))
           else
              PlotonRenderings(itype) = .false.
           endif
        elseif (.not.iplotpartoftype(itype)) then
           PlotonRenderings(itype) = .false.
        endif
     enddo
     return           
!------------------------------------------------------------------------
  case(5)
     print*,'(0 Square) (1 .) (2 +) (3 *) (4 o) (5 x) (17 bold circle) (-8 bigger bold circle)'
     do itype=1,ntypes
        call prompt(' Enter PGPLOT marker for '//trim(labeltype(itype)) &
             //' particles:',imarktype(itype),-8,31)
     enddo
     return   
!------------------------------------------------------------------------
  case(6)
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
        print*,'modifying plot limits for new coordinate system'
        call coord_transform_limits(lim(1:ndim,1),lim(1:ndim,2), &
                                    icoordsprev,icoordsnew,ndim)
     endif
     return
!------------------------------------------------------------------------
  case(7)
     call submenu_exact(iexact)
     return
!------------------------------------------------------------------------
  case(8)
     call options_exact
     return     
!------------------------------------------------------------------------
  case default
     return

  end select

  return      
end subroutine submenu_particleplots

end module settings_part
