!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2019 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module containing settings and options relating to particle plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_part
 use params
 use settings_data, only:icoordsnew,iexact
 implicit none
 integer, dimension(maxparttypes) :: imarktype,idefaultcolourtype,itypeorder
 integer, dimension(100)          :: icircpart
 integer, dimension(maxplot)      :: ilocerrbars
 logical, dimension(maxparttypes) :: iplotpartoftype,PlotOnRenderings,UseTypeInContours
 integer :: ncircpart,ismooth_particle_plots
 integer :: linestyle, linecolour,linestylethisstep,linecolourthisstep,ErrorBarType
 logical :: iplotline,ilabelpart,ifastparticleplot,iploterrbars
 real    :: hfacmarkers,rref,betaflare,mstari

 namelist /plotopts/ iplotline,linestyle,linecolour, &
   imarktype,iplotpartoftype,PlotOnRenderings, &
   iexact,icoordsnew,ifastparticleplot,idefaultcolourtype,&
   itypeorder,UseTypeInContours,iploterrbars,ilocerrbars,hfacmarkers,&
   ErrorBarType,ismooth_particle_plots,rref,betaflare,mstari

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_part
 use settings_data, only:icoords
 implicit none
 integer :: i

 ncircpart = 0
 iplotline = .false.     ! plot line joining the particles
 linestyle = 1           ! line style for above
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
 imarktype = 1              ! marker type for all particles
 imarktype(2) = 4           ! marker type for ghost/dark matter particles
 imarktype(3) = 17          ! marker type for sink particles
 imarktype(5) = 3           ! marker type for star particles (gadget)
 idefaultcolourtype = -1     ! default colour for each particle type
 ifastparticleplot = .true. ! allow crowded-field elimination on particle plots
 do i=1,maxparttypes
    itypeorder(i) = i
 enddo
 UseTypeInContours(:) = iplotpartoftype(:)
 iploterrbars = .false.    ! plot error bars for a particular column
 ilocerrbars(:) = 0     ! location of data for error bars in dat array
 hfacmarkers = 1.0
 ErrorBarType = 0
 ismooth_particle_plots = 0
 rref = 1.
 betaflare = 1.25
 mstari = 1.  ! used in Toomre Q calculation

 return
end subroutine defaults_set_part

!---------------------------------------------
! changed default values for these options
!---------------------------------------------
subroutine defaults_set_part_ev
 implicit none

 iplotline = .true.     ! plot line joining the particles
 iplotpartoftype(1:maxparttypes) = .false. ! whether or not to plot particles of certain types
 UseTypeInContours(:) = iplotpartoftype(:)

 return
end subroutine defaults_set_part_ev

!---------------------------------------------
! changed default values for these options
!---------------------------------------------
subroutine initialise_coord_transforms
 use geometry, only:set_flaring_index

 call set_flaring_index(rref,betaflare)

end subroutine initialise_coord_transforms

!----------------------------------------------------------------------
! submenu with options relating to particle plots
!----------------------------------------------------------------------
subroutine submenu_particleplots(ichoose)
 use exact,           only:options_exact,submenu_exact
 use labels,          only:labeltype,ih,label,idustfracsum,ideltavsum
 use limits,          only:lim
 use settings_data,   only:icoords,ntypes,ndim,ndimV,UseTypeInRenderings, &
                            ndataplots,ndusttypes,idustfrac_plot,ideltav_plot,ncalc
 use settings_render, only:iplotcont_nomulti
 use particle_data,   only:npartoftype,iamtype
 use prompting,       only:prompt,print_logical
 use geometry,        only:maxcoordsys,labelcoordsys,coord_transform_limits,&
                            igeom_flaredcyl,igeom_logflared,set_flaring_index
 use multiplot,       only:itrans
 use plotlib,         only:plotlib_maxlinestyle,plotlib_maxlinecolour
 use calcquantities,  only:calc_quantities
 use settings_data,   only:DataIsBuffered,numplot
 use filenames,       only:nsteps,nstepsinfile,ifileopen
 use geomutils,       only:set_coordlabels
 use calcquantities,  only:setup_calculated_quantities
 use asciiutils,      only:enumerate
 implicit none
 integer, intent(in) :: ichoose
 integer             :: i,iaction,n,itype,icoordsprev,ierr,icol
 integer             :: idustfrac_prev
 character(len=2)    :: charntypes
 character(len=20)   :: substring1,substring2,substring3
 character(len=1000) :: fmtstring
 character(len=120)  :: contline
 character(len=3)    :: idustfracsum_string
 logical :: ians

 iaction = ichoose

 !--we require some tricks with the format string to print only the actual number of
 !  particle types rather than the whole array
 !
 if (ntypes > 100) print*,'WARNING: Internal error: ntypes too large for formatting in particle plot menu'
 if (ntypes <= 0) then
    substring1 = "no types specified"
    substring2 = "not applicable"
 elseif (ntypes==1) then
    substring1 = "a"
    substring2 = "i2"
 else
    write(charntypes,"(i2)") ntypes-1
    substring1 = charntypes//"(a,',',1x),a"
    substring2 = charntypes//"(i2,',',1x),i2"
 endif
 substring3 = enumerate(ismooth_particle_plots+1,(/'OFF  ','FIXED','ADAPT'/),default=1)
 if (iplotcont_nomulti) then
    contline = "'            use in contour plots:       ( ',"//trim(substring1)//",' )',/,"
 else
    contline = ' '
 endif

 fmtstring="("// &
         "' 0) exit ',/,"// &
         "' 1) turn on/off particles by type       ( ',"//trim(substring1)//",' )',/,"//trim(contline)// &
         "' 2) change graph markers for each type  ( ',"//trim(substring2)//",' )',/,"//  &
         "' 3) set colour for each particle type   ( ',"//trim(substring2)//",' )',/,"//  &
         "' 4) smooth particle plots               ( ',a,' )',/,"//  &
         "' 5) plot line joining particles         ( ',a,' ) ',/,"// &
         "' 6) plot error bars/smoothing circles   ( ',a,' ) ',/,"// &
         "' 7) change coordinate systems           ( ',i2,' ) ',/,"// &
         "' 8) plot exact solution                 ( ',i2,' ) ',/,"// &
         "' 9) exact solution plot options ')"

 print "(a)",'------------- particle plot options -------------------'
 if (iaction <= 0 .or. iaction > 9) then
    if (iplotcont_nomulti) then
       print fmtstring,(trim(print_logical(iplotpartoftype(i))),i=1,ntypes), &
                 (trim(print_logical(UseTypeInContours(i),mask=UseTypeInRenderings(i))),i=1,ntypes), &
                 imarktype(1:ntypes),idefaultcolourtype(1:ntypes),trim(substring3), &
                 print_logical(iplotline),print_logical(ncircpart > 0 .or.iploterrbars),icoordsnew,iexact
    else
       print fmtstring,(trim(print_logical(iplotpartoftype(i))),i=1,ntypes), &
                 imarktype(1:ntypes),idefaultcolourtype(1:ntypes),trim(substring3), &
                 print_logical(iplotline),print_logical(ncircpart > 0 .or.iploterrbars),icoordsnew,iexact
    endif
    call prompt('enter option',iaction,0,9)
 endif
!
 select case(iaction)
!------------------------------------------------------------------------
 case(1)
    !          plot particles by type?
    do itype=1,ntypes
       if (UseTypeinRenderings(itype) .and. ndim > 1) then
          call prompt('Plot '//trim(labeltype(itype))//' particles / use in renderings?',iplotpartoftype(itype))
          if (iplotcont_nomulti) then
             call prompt('Use '//trim(labeltype(itype))//' particles in contour plots?',UseTypeInContours(itype))
          endif
       else
          call prompt('Plot '//trim(labeltype(itype))//' particles?',iplotpartoftype(itype))
          UseTypeInContours(itype) = .false.
       endif
       if (iplotpartoftype(itype) .and. itype > 1) then
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
    print "(/,' Marker options (for all from -8->31, see plot library userguide):',11(/,i2,') ',a))", &
           0,'square',1,'.',2,'+',3,'*',4,'o',5,'x',17,'bold circle',-8,'large bold circle', &
           32,'solid circle, size proportional to h', &
           33,'open circle,  size proportional to h', &
           34,'outlined solid circle, size prop. to h'

    !print*,'(0 Square) (1 .) (2 +) (3 *) (4 o) (5 x) (17 bold circle) (-8 bigger bold circle)'
    do itype=1,ntypes
       call prompt(' Enter marker to use for '//trim(labeltype(itype)) &
             //' particles:',imarktype(itype),-8,35)
    enddo
    if (any(imarktype(1:ntypes) >= 32)) then
       print*
       call prompt(' Enter proportionality factor for scalable markers (radius = fac*h)',hfacmarkers)
    endif

    return
!------------------------------------------------------------------------
 case(3)
    print "(2(a,/),/,4(a,/))", &
           ' Warning: setting a colour for a particle type overrides', &
           '          (at each new timestep) colours set interactively ', &
           ' -1 = retain interactively set colours between timesteps', &
           '  0 = background ',&
           '  1 = foreground ',&
           '  2->10 = various colours (see default colour indices for plot library)'
    do itype=1,ntypes
       call prompt(' Enter default colour for '//trim(labeltype(itype)) &
             //' particles:',idefaultcolourtype(itype),-1,14)
    enddo
    return
!------------------------------------------------------------------------
 case(4)
    !print "(3(/,a))",' 0) no smoothing, raw particle plots',&
    !                 ' 1) render particle plots with fixed h', &
    !                 ' 2) render particle plots with adaptive h (slower)'
    call prompt('smooth particle plots? (0=none 1=fixed 2=adaptive)',ismooth_particle_plots,0,2)
    if (ismooth_particle_plots==0) then
       if (size(iamtype(:,1)) > 1) then
          print "(3(/,a),/)", &
             ' WARNING: changing type plotting order currently has no effect ', &
             '          when particle types are mixed in the dump file', &
             '          (for sphNG read disable this using -lowmem on the command line)'
       endif
       if (ntypes > 1) then
          ians = .false.
          call prompt('set plot order of particle types manually?',ians)
          if (ians) then
             print "(9(i1,'=',a,', '))",(i,trim(labeltype(i)),i=1,ntypes)
             call prompt('enter first particle type to plot',itypeorder(1),1,ntypes)
             do i=2,ntypes
                ierr = 1
                do while (ierr /= 0)
                   itype = itypeorder(i)
                   call prompt('enter next particle type to plot',itype,1,ntypes)
                   if (any(itypeorder(1:i-1)==itype)) then
                      print "(a)",' error: cannot be same as previous type'
                      ierr = 1
                   else
                      itypeorder(i) = itype
                      ierr = 0
                   endif
                enddo
             enddo
          endif
       endif

       print "(/,a,/,a,/)",' Fast particle plotting excludes particles in crowded regions', &
                            ' Turn this option off to always plot every particle'
       call prompt('Allow fast particle plotting?',ifastparticleplot)
    endif
    return
!------------------------------------------------------------------------
 case(5)
    call prompt('plot line joining particles?',iplotline)
    if (iplotline) then
       call prompt('Enter line style to use ',linestyle,1,plotlib_maxlinestyle)
       call prompt('Enter colour for line ',linecolour,0,plotlib_maxlinecolour)
    endif
    return
!------------------------------------------------------------------------
 case(6)
    if (ndim <= 1 .or. ih <= 0) then
       icol = 0
       do icol=1,ndataplots
          if (ilocerrbars(icol) > 0) print "(a,i2,a,i2,a)", &
              'column ',ilocerrbars(icol),' contains errors for column ',icol,':'//label(icol)
       enddo
       if (any(ilocerrbars(1:ndataplots) > 0)) then
          call prompt('turn on plotting of error bars? ',iploterrbars)
          if (.not.iploterrbars) return
       else
          iploterrbars = .false.
       endif
       icol = 1
       do while(icol /= 0)
          icol = 0
          call prompt('Enter column to set location of error bars for (0=none)',icol,0,ndataplots)
          if (icol > 0) then
             call prompt('Enter location of error data for this column in the data',ilocerrbars(icol),0)
             if (ilocerrbars(icol) <= 0 .or. ilocerrbars(icol) > ndataplots) then
                print "(a,i2)",' WARNING: currently no data in column ',ilocerrbars(icol)
             else
                iploterrbars = .true.
             endif
             if (all(ilocerrbars(1:ndataplots) <= 0)) iploterrbars = .false.
          else
             if (all(ilocerrbars(1:ndataplots) <= 0)) iploterrbars = .false.
          endif
       enddo
       print "(2(/,a))",' 0) Default style |--|',&
                        ' 1) Semi-transparent shaded region'
       call prompt('Select error bar style',ErrorBarType,0,1)
    else
       print*,'Circles of interaction can also be set interactively'
       call prompt('Enter number of circles to draw',ncircpart,0,size(icircpart))
       if (ncircpart > 0) then
          do n=1,ncircpart
             if (icircpart(n)==0) then
                if (n > 1) then
                   icircpart(n) = icircpart(n-1)+1
                else
                   icircpart(n) = 1
                endif
             endif
             call prompt('Enter particle number to plot circle around', &
                       icircpart(n),1,maxval(npartoftype(1,:)))
          enddo
       endif
    endif
    return
!------------------------------------------------------------------------
 case(7)
    print "(' 0) reset (=',i2,')')",icoords
    do i=1,maxcoordsys
       print "(1x,i1,')',1x,a)",i,labelcoordsys(i)
    enddo
    icoordsprev = icoordsnew
    call prompt(' Enter coordinate system to plot in:', &
                 icoordsnew,0,maxcoordsys)
    if (icoordsnew==0) icoordsnew = icoords
    select case(icoordsnew)
       !case(igeom_rotated)
       !  call prompt('enter rotation angle a (degrees)',rot_angle_a)
       ! call prompt('enter rotation angle b (degrees)',rot_angle_b)
       !call set_rotation_angles(rot_angle_a*pi/180.,rot_angle_b*pi/180.)
    case(igeom_flaredcyl,igeom_logflared)
       call prompt('enter reference radius (Rref) in z''=z(R/Rref)**(-beta)',rref)
       print "(2(/,a),/)",' Use beta = 1.5 - q to give correct flaring index in a disc', &
                         ' where q is sound speed index i.e. cs = cs_0(R/R_0)^-q '
       call prompt('enter flaring index beta',betaflare)
       call set_flaring_index(rref,betaflare)
    end select

    if (icoordsnew /= icoordsprev) then
       itrans(1:ndim) = 0
       call coord_transform_limits(lim(1:ndim,1),lim(1:ndim,2), &
                                    icoordsprev,icoordsnew,ndim)
       call set_coordlabels(numplot)
       if (DataIsBuffered) then
          call calc_quantities(1,nsteps)
       else
          call calc_quantities(1,nstepsinfile(ifileopen))
       endif
    endif
    return
!------------------------------------------------------------------------
 case(8)
    call submenu_exact(iexact)
    return
!------------------------------------------------------------------------
 case(9)
    call options_exact(iexact)
    return
!------------------------------------------------------------------------
 case default
    return

 end select

 return
end subroutine submenu_particleplots

end module settings_part
