!-------------------------------------------------------------------------
! Module containing settings and options relating to particle plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_part
 use params
 implicit none
 integer, dimension(maxparttypes) :: imarktype
 integer :: ncircpart, icoordsnew
 integer, dimension(:), allocatable :: icircpart
 integer :: nc
 integer :: linestylein, iexact
 logical, dimension(maxparttypes) :: iplotpartoftype
 logical :: iplotline,iplotlinein,ilabelpart

 namelist /plotopts/ iplotline,iplotlinein,linestylein,  &
   imarktype,iplotpartoftype,iexact,icoordsnew

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_part
  implicit none

  ncircpart = 0
  iplotline = .false.     ! plot line joining the particles
  iplotlinein = .false.   ! " " but on first step only
  linestylein = 4         ! PGPLOT line style for above
  iexact = 0              ! exact solution to plot
  ilabelpart = .false.    ! plot particle numbers
  
  iplotpartoftype(1) = .true. ! whether or not to plot particles of certain types
  iplotpartoftype(2:maxparttypes) = .false.
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
  use labels
  use settings_data, only:icoords,ntypes
  use particle_data, only:npartoftype
  use prompting
  use geometry, only:maxcoordsys,labelcoordsys
  implicit none
  integer :: i,iaction,n,itype
  character(LEN=1) :: ans

  iaction = 0      
  print 10, iplotline,iplotlinein,ilabelpart,ncircpart, &
        iplotpartoftype,imarktype,icoordsnew,iexact
10  format(' 0) exit ',/,                 &
         ' 1) toggle plot line                   ( ',L1,',',1x,L1,' ) ',/, &
         ' 2) toggle label particles             ( ',L1,' ) ',/,           &
         ' 3) plot circles of interaction        ( ',i3,' ) ',/,           &
         ' 4) toggle plot particles by type      ( ',6(L1,',',1x),' )',/,  &
         ' 5) change graph markers for each type ( ',6(i2,',',1x),' )',/,  &
         ' 6) change coordinate systems          ( ',i2,' ) ',/,           &
         ' 7) plot exact solution                ( ',i2,' ) ',/, &
         ' 8) exact solution options')
    call prompt('enter option',iaction,0,8)
!
  select case(iaction)

!------------------------------------------------------------------------
  case(1)
     print*,' Plot initial only(i), all(a), both(b) or not (n)?'
     read*,ans
     iplotline = .false.
     iplotlinein = .false.
     if (ans.eq.'i'.or.ans.eq.'b') iplotlinein = .true.
     if (ans.eq.'a'.or.ans.eq.'b') iplotline = .true.
     if (iplotlinein) then
        call prompt('Enter PGPLOT line style',linestylein,0,5)
     endif
     print*,' Plot line = ',iplotline,iplotlinein
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
     enddo
     return           
!------------------------------------------------------------------------
  case(5)
     print*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
     do itype=1,ntypes
        call prompt(' Enter PGPLOT marker for '//trim(labeltype(itype)) &
             //' particles:',imarktype(itype))
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
     call prompt(' Enter coordinate system to plot in:', &
                 icoordsnew,0,maxcoordsys)
     if (icoordsnew.eq.0) icoordsnew = icoords
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
