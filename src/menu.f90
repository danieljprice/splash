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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!--------------------
!  SPLASH MAIN MENU
!--------------------
module mainmenu
 implicit none
 public :: menu,allowrendering,set_extracols

 private

contains

subroutine menu
  use filenames,        only:defaultsfile,limitsfile,fileprefix,set_filenames
  use labels,           only:label,labelvec,iamvec,isurfdens,itoomre,ipdf,icolpixmap,is_coord,ix
  use limits,           only:write_limits,lim2,lim,reset_lim2,lim2set
  use options_data,     only:submenu_data
  use settings_data,    only:ndim,numplot,ndataplots,nextra,ncalc,ivegotdata, &
                             buffer_data,ncolumns,icoords,icoordsnew
  use settings_limits,  only:submenu_limits,iadapt
  use settings_part,    only:submenu_particleplots
  use settings_page,    only:submenu_page,submenu_legend,interactive,nacross,ndown
  use settings_render,  only:submenu_render,iplotcont_nomulti,icolours,double_rendering
  use settings_vecplot, only:submenu_vecplot,iplotpartvec
  use settings_xsecrot, only:submenu_xsecrotate,xsec_nomulti
  use multiplot
  use prompting,        only:prompt,print_logical
  use transforms,       only:transform_label
  use defaults,         only:defaults_write
  use getdata,          only:get_data
  use geomutils,        only:set_coordlabels
  use timestepping
  implicit none
  integer            :: i,icol,ihalf,iadjust,indexi,ierr
  integer            :: ipicky,ipickx,irender,ivecplot,icontourplot
  integer            :: iamvecprev, ivecplottemp,ichoose
  character(len=2)   :: ioption
  character(len=100) :: vecprompt,string
  character(len=20)  :: rprompt
  character(len=*), parameter :: sep="(55('-'))"
  logical            :: iAllowRendering

  irender = 0
  icontourplot = 0
  ivecplot = 0
  if (ndim > 1 .and. ix(1) > 0) then
     ipickx = ix(1)
  else
     ipickx = 1
  endif
  ipicky = 1

  menuloop: do
!---------------------------------------------------------------------------
!  preliminaries
!---------------------------------------------------------------------------
!
!--make sure the number of columns is set appropriately
!  (nextra can change depending on what options are set)
!
  !
  !--numplot is the total number of data columns (read + calculated)
  !   not including the particle co-ordinates
  !  nextra are extra graphs to plot (e.g. convergence plots, power spectrum)
  !
  !  note that numplot and ndataplots should *only* be set here
  !  this means that even if ncolumns changes during data reads while plotting
  !  we don't start plotting new quantities
  !
  call set_extracols(ncolumns,ncalc,nextra,numplot,ndataplots)
!
!--set the coordinate and vector labels
!  if working in a different coordinate system
!
  call set_coordlabels(numplot)

!--set contents of the vector plotting prompt
  vecprompt(1:6) = '0=none'
  indexi = 7
  iamvecprev = 0
  do icol=1,numplot
     if (iamvec(icol).ne.0 .and. iamvec(icol).ne.iamvecprev) then
        iamvecprev = iamvec(icol)
        if (iamvec(icol).ge.10) then
           write(vecprompt(indexi:),"(',',1x,i2,'=',a)") &
                 iamvec(icol),trim(labelvec(icol))
        else
           write(vecprompt(indexi:),"(',',1x,i1,'=',a)") &
                 iamvec(icol),trim(labelvec(icol))
        endif
        indexi = len_trim(vecprompt) + 1
     endif
  enddo

  ichoose = 0

!---------------------------------------------------------------------------
!  print menu
!---------------------------------------------------------------------------

  if (numplot.gt.0) then
!
!--data columns
!
     call print_header()
     print sep
     ihalf = numplot/2                ! print in two columns
     iadjust = mod(numplot,2)
     print "(1x,i2,')',1x,a20,1x,i2,')',1x,a20)", &
          (i,transform_label(label(i),itrans(i)), &
          ihalf + i + iadjust, transform_label(label(ihalf + i + iadjust), &
          itrans(ihalf+i+iadjust)),i=1,ihalf)
     if (iadjust.ne.0) then
        print "(1x,i2,')',1x,a20)", &
              ihalf + iadjust,transform_label(label(ihalf + iadjust), &
              itrans(ihalf+iadjust))
     endif
!
!--multiplot
!
     print sep
     print "(1x,i2,')',1x,a,'[ ',i2,' ]',5x,a2,') ',a)", &
           numplot+1,'multiplot ',nyplotmulti,'m','set multiplot '
  else
!
!--if no data
!
     print "(/a)",' No data: You may choose from the options below '
  endif
!
!--options
!
  print sep
  if (ndim.le.1) then
     print "(a)",' d(ata) p(age) o(pts) l(imits) le(g)end s,S(ave) q(uit)'
  else
     print "(a)",' d(ata) p(age) o(pts) l(imits) le(g)end h(elp)'
     print "(a)",' r(ender) v(ector) x(sec/rotate) s,S(ave) q(uit)'
  endif
  print sep

!
!--prompt user for selection
!
  write(*,"(a)",ADVANCE='NO') 'Please enter your selection now (y axis or option):'
  read(*,"(a)",iostat=ierr) string
  if (ierr < 0) stop 'reached end of input' ! end of input (e.g. in script)
  if (ierr > 0) stop !'error reading input'
  ioption = string(1:2)

!------------------------------------------------------------
!  if input is an integer and within range, plot data
!------------------------------------------------------------
  read(ioption,*,iostat=ierr) ipicky
  if (ierr /= 0) ipicky = -1
  
  !--try to read more integers from the string
  !  if present, use these to set up an "instant multiplot"
  if (ipicky > 0 .and. ipicky < numplot+1 .and. len_trim(string) > 2) then
     call set_instant_multiplot(string,ipicky,ipickx,numplot,nyplotmulti,&
                                multiplotx,multiploty,nacross,ndown)
  endif

  if (ipicky > 0 .and. ipicky <= numplot+1) then

     if (.not.ivegotdata) then
        !
        !--do not allow plotting if no data - instead try to read data
        !
        print*,' no data '
        if (buffer_data) then
           call get_data(-1,.false.)
        else
           call get_data(1,.false.,firsttime=.true.)
        endif
     else
        !
        !--if needed prompt for x axis selection
        !
        if (ipicky <= (numplot-nextra)) then
           if (ipickx==0) then
              if (ndim > 1 .and. ix(1) > 0) then
                 ipickx = ix(1)
              else
                 ipickx = 1 ! do not allow zero as default
              endif
           endif
           if (ipickx==ipicky) then ! do not allow x same as y by default
              if (ipickx > 1) then
                 ipickx = ipicky-1
              else
                 ipickx = ipicky+1
              endif
           endif
           call prompt(' (x axis) ',ipickx)
           !--go back to y prompt if out of range
           if (ipickx.gt.numplot .or. ipickx.le.0) cycle menuloop
           !
           !--work out whether rendering is allowed
           !
           iAllowRendering = allowrendering(ipickx,ipicky,xsec_nomulti)
           !
           !--prompt for render and vector plots
           ! -> only allow if in "natural" coord system, otherwise h's would be wrong)
           ! (a future feature might be to interpolate in icoord then translate the pixels
           !  to icoordsnew, or alternatively plot non-cartesian pixel shapes)
           ! -> also do not allow if transformations are applied
           !
           if (is_coord(ipicky,ndim) .and. is_coord(ipickx,ndim)) then
              if (iAllowRendering) then
                 call prompt('(render) (0=none)',irender,0,(numplot-nextra))
                 if (irender > 0 .and. iplotcont_nomulti .and. icolours /= 0) then
                    if (double_rendering) then
                       rprompt = '2nd render'
                    else
                       rprompt = 'contours'
                    endif
                    call prompt('('//trim(rprompt)//') (0=none)',icontourplot,0,(numplot-nextra))
                    if (icontourplot==irender) then
                       if (iadapt) then
                          print "(a)",' limits for '//trim(rprompt)//' are adaptive'
                       else
                          if (.not.lim2set(icontourplot)) lim2(icontourplot,:) = lim(icontourplot,:)
                          call prompt(' enter min for '//trim(rprompt)//':',lim2(icontourplot,1))
                          call prompt(' enter max for '//trim(rprompt)//':',lim2(icontourplot,2))
                          if (all(abs(lim2(icontourplot,:)-lim(icontourplot,:)) < tiny(lim))) then
                             call reset_lim2(icontourplot)
                          endif
                       endif
                    endif
                 endif
              else
                 irender = 0
              endif
              if (any(iamvec(1:numplot).ne.0) .and. (icoordsnew.eq.icoords)) then
                 ivecplottemp = -1
                 ierr = 1
                 do while(ierr /= 0 .and. ivecplottemp /= 0)
                    ivecplottemp = ivecplot
                    ierr = 0
                    call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
                    if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
                       print "(a)",'Error, value not in list'
                       ierr = 1
                    endif
                 enddo
                 ivecplot = ivecplottemp
              else
                 ivecplot = 0
              endif
              if (ivecplot.gt.0 .and. irender.eq.0) then
                 call prompt('plot particles?',iplotpartvec)
              endif
           else
              irender = 0
              ivecplot = 0
           endif
        elseif (ipicky > 0 .and. ipicky==itoomre .or. ipicky==isurfdens) then
            if (ipicky==isurfdens) print "(a)",' setting x axis to r for surface density plot'
            if (ipicky==itoomre) print "(a)",' setting x axis to r for Toomre Q plot'
            ipickx = 1
            irender = 0
            ivecplot = 0
        elseif (ipicky > 0 .and. ipicky==ipdf) then
            call prompt(' enter x axis for PDF calculation ',ipickx,1,ndataplots)
            irender = 0
            ivecplot = 0
        elseif (ipicky > 0 .and. ipicky==icolpixmap) then
            call prompt(' enter corresponding SPH column for particle data ',irender,0,ndataplots)
            ipickx = 0
            ivecplot = 0
        elseif (ipicky==numplot+1) then
        !
        !--for multiplots, check that options are valid. If not, re-prompt for multiplot
        !  settings
        !
            ipickx   = 0
            irender  = 0
            ivecplot = 0
            if (any(multiploty(1:nyplotmulti) <= 0) .or. &
                any(multiploty(1:nyplotmulti) > numplot) .or. &
                any(multiplotx(1:nyplotmulti) <= 0) .or. &
                any(multiplotx(1:nyplotmulti) > numplot)) then
               print "(/,a,/)",'ERROR: multiplot settings out of range, please re-enter these'
               call options_multiplot
            endif
        endif
        !
        !--call main plotting routine
        !
        call timestep_loop(ipicky,ipickx,irender,icontourplot,ivecplot)
     endif
!------------------------------------------------------------------------
!  if input is an integer > numplot+1, quit
!------------------------------------------------------------------------
  elseif (ipicky > numplot+1) then
     return
  else
!------------------------------------------------------------------------
!  if input is a string, use menu options
!------------------------------------------------------------------------
!--  Menu shortcuts; so you can type e.g. o2 and get the o)ptions menu, item 2
     read(ioption(2:2),*,iostat=ierr) ichoose
     if (ierr /= 0) ichoose = 0

     select case(ioption(1:1))
!------------------------------------------------------------------------
!+ Sets up plotting of (m)ultiple quantities per timestep
     case('m','M')
        call options_multiplot
!------------------------------------------------------------------------
!+ This submenu sets options relating to the (d)ata read
     case('d','D')
        call submenu_data(ichoose)
!------------------------------------------------------------------------
!+ This option turns (i)nteractive mode on/off
     case('i','I')
        interactive = .not.interactive
        print "(a)",' Interactive mode is '//print_logical(interactive)
!------------------------------------------------------------------------
!+ This submenu sets (p)age setup options
     case('p','P')
        call submenu_page(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets particle plot (o)ptions
     case('o','O')
        call submenu_particleplots(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets le(g)end and title options
     case('g','G')
        call submenu_legend(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets (r)endering options
     case('r','R')
        if (ndim.le.1) print "(a)",'WARNING: these options have no effect in < 2D'
        call submenu_render(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets (v)ector plotting options
     case('v','V')
        if (ndim.le.1) print "(a)",'WARNING: these options have no effect in < 2D'
        call submenu_vecplot(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets cross section and rotation options
     case('x','X')
        if (ndim.le.1) print "(a)",'WARNING: these options have no effect in < 2D'
        call submenu_xsecrotate(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets options relating to the plot limits
     case('l','L')
        call submenu_limits(ichoose)
!------------------------------------------------------------------------
!+ The (s)ave option saves the default options to a
!+ file called `splash.defaults'' in the current directory which
!+ is read automatically upon the next invocation of splash.
!+ This file uses namelist formatting and may be edited
!+ manually prior to startup if so desired. This is quite
!+ useful for setting multiplots with many plots per page
!+ The (S)ave option writes both the defaults file and
!+ also saves the current plot limits to a file called
!+ 'splash.limits' which is also read automatically
!+ at startup.
     case('s')
        if (ioption(2:2).eq.'a') then
           call prompt('enter prefix for defaults file: ',fileprefix,noblank=.true.)
           if (index(fileprefix,'.defaults').eq.0) then
              defaultsfile = trim(fileprefix)//'.defaults'
           else
              defaultsfile = trim(fileprefix)
           endif
        endif
        call defaults_write(defaultsfile)
     case('S')
        if (ioption(2:2).eq.'a' .or. ioption(2:2).eq.'A') then
           call prompt('enter prefix for filenames: ',fileprefix,noblank=.true.)
           call set_filenames(trim(fileprefix))
        endif
        call defaults_write(defaultsfile)
        call write_limits(limitsfile)
!------------------------------------------------------------------------
!+ Slightly obsolete: prints whatever help may be helpful
     case('h','H')
        print "(10(/a))",' Hint: menu items can be shortcut by typing, e.g. o2 for ',&
                 ' the o)ptions menu, item 2.',' ', &
                 ' for detailed help, consult the user guide',' ',&
                 '  (splash/docs/splash.pdf ',&
                 '   or http://users.monash.edu.au/~dprice/splash/userguide/)', &
                 ' ',' and/or the online FAQ. If you''re really stuck, email me! '
        read*
!------------------------------------------------------------------------
!+ (q)uit, unsurprisingly, quits. Typing a number greater
!+ than the number of data columns also exits the program
!+ (e.g. I often simply type 99 to exit).
     case('q','Q')
        return
!------------------------------------------------------------------------
     case DEFAULT
        print "(a)",'unknown option '//trim(ioption)
     end select

  endif

  enddo menuloop

  return

 contains

!----------------------------------------------------
! multiplot setup
!----------------------------------------------------
  subroutine options_multiplot
   use settings_page,   only:nacross, ndown
   use settings_render, only:iplotcont_nomulti
   use limits,          only:lim,lim2,lim2set,reset_lim2
   use labels,          only:is_coord,labeltype
   use params,          only:maxparttypes
   use settings_data,   only:ntypes
   implicit none
   integer :: ifac,ierr,itype,nvalues
   logical :: isamex, isamey, icoordplot, anycoordplot, imultisamepanel
   integer, dimension(maxparttypes)   :: itypelist

   call prompt('Enter number of plots per timestep:',nyplotmulti,1,numplot)

   isamey = all(multiploty(1:nyplotmulti).eq.multiploty(1))
   if (ndim.ge.2) call prompt('Same y axis for all?',isamey)
   if (isamey) then
      call prompt('Enter y axis for all plots',multiploty(1),1,numplot)
      multiploty(2:nyplotmulti) = multiploty(1)
   endif

   isamex = all(multiplotx(1:nyplotmulti).eq.multiplotx(1))
   call prompt('Same x axis for all?',isamex)
   if (isamex) then
      call prompt('Enter x axis for all plots',multiplotx(1),1,numplot)
      multiplotx(2:nyplotmulti) = multiplotx(1)
   endif

   anycoordplot = .false.
   do i=1,nyplotmulti
      print*,'-------------- Plot number ',i,' --------------'
      if (.not.isamey) then
         call prompt(' y axis ',multiploty(i),1,numplot)
      endif
      if (multiploty(i).le.ndataplots .and. .not.isamex) then
         call prompt(' x axis ',multiplotx(i),1,ndataplots)
      else
         if (multiploty(i).eq.isurfdens) then
            print "(a)",' setting x axis to r for surface density plot'
            multiplotx(i) = 1
         elseif (multiploty(i).eq.itoomre) then
            print "(a)",' setting x axis to r for Toomre Q plot'
            multiplotx(i) = 1
         elseif (multiploty(i).eq.ipdf) then
            call prompt(' enter x axis for PDF calculation ',multiplotx(i),1,ndataplots)
         elseif (multiploty(i).eq.icolpixmap) then
            call prompt(' enter corresponding SPH column for particle data ',irendermulti(i),0,ndataplots)
            multiplotx(i) = 1
         elseif(.not.isamex) then
            multiplotx(i) = multiploty(i)
         endif
      endif
      !
      !--work out whether rendering is allowed
      !
      iAllowRendering = allowrendering(multiplotx(i),multiploty(i))

      icoordplot = (is_coord(multiplotx(i),ndim) .and. is_coord(multiploty(i),ndim))
      if (icoordplot) anycoordplot = icoordplot

      if (icoordplot) then
         if (iAllowRendering) then
            call prompt('(render) (0=none)',irendermulti(i),0,numplot-nextra)
            if (irendermulti(i).gt.0 .and. iplotcont_nomulti .and. icolours.ne.0) then
               if (double_rendering) then
                  rprompt = '2nd render'
               else
                  rprompt = 'contours'
               endif
               call prompt('('//trim(rprompt)//') (0=none)',icontourmulti(i),0,numplot-nextra)
               if (icontourmulti(i).eq.irendermulti(i)) then
                  if (iadapt) then
                     print "(a)",' limits for '//trim(rprompt)//' are adaptive '
                  else
                     if (.not.lim2set(icontourmulti(i))) lim2(icontourmulti(i),:) = lim(icontourmulti(i),:)
                     call prompt(' enter min for '//trim(rprompt)//':',lim2(icontourmulti(i),1))
                     call prompt(' enter max for '//trim(rprompt)//':',lim2(icontourmulti(i),2))
                     if (all(abs(lim2(icontourmulti(i),:)-lim(icontourmulti(i),:)) < tiny(lim))) then
                        call reset_lim2(icontourmulti(i))
                     endif
                  endif
               endif
            else
               icontourmulti(i) = 0
            endif
            !iplotcontmulti(i) = iplotcont_nomulti
         endif
         if (any(iamvec(1:numplot).gt.0)) then
            ivecplottemp = -1
            ierr = 1
            do while(ierr.ne.0 .and. ivecplottemp.ne.0)
               ivecplottemp = ivecplotmulti(i)
               ierr = 0
               call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
               if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
                  print "(a)",'Error, value not in list'
                  ierr = 1
               endif
            enddo
            ivecplotmulti(i) = ivecplottemp
         else
            ivecplotmulti(i) = 0
         endif
         if (ivecplotmulti(i).gt.0 .and. irendermulti(i).eq.0) then
            call prompt('plot particles?',iplotpartvec)
         endif
      else
         !
         !--set irender, icontour and ivecplot to zero if no rendering allowed
         !
         if (multiploty(i).ne.icolpixmap) irendermulti(i) = 0
         icontourmulti(i) = 0
         ivecplotmulti(i) = 0
      endif

      if (icoordplot .and. ndim.ge.2) then
         call prompt(' is this a cross section (no=projection)? ',x_secmulti(i))
         if (x_secmulti(i)) then
            call prompt('enter co-ordinate location of cross section slice',xsecposmulti(i))
         endif
      endif
      !
      !--prompt for selection of different particle types
      !  if more than one SPH particle type is present
      !
      itypelist = 0
      if (ntypes.ge.2) then
         call prompt('use all active particle types?',iusealltypesmulti(i))

         if (iusealltypesmulti(i)) then
            nvalues = 0
            itypelist(:) = 0
         else
         !
         !--prepare list of types based on current iplotpartoftypemulti
         !
            nvalues = 0
            do itype=1,ntypes
               if (iplotpartoftypemulti(itype,i)) then
                  nvalues = nvalues + 1
                  itypelist(nvalues) = itype
               endif
            enddo
            if (nvalues.eq.0) then
               print*,'warning: internal error in type list'
               itypelist(:) = 0
               nvalues = 1
            endif
         !
         !--prompt for list of types to use
         !
            do itype=1,ntypes
               print "(i2,':',1x,a)",itype,'use '//trim(labeltype(itype))//' particles'
            enddo

            call prompt('Enter type or list of types to use',itypelist,nvalues,1,ntypes)
         !
         !--set which particle types to plot
         !
            iplotpartoftypemulti(:,i) = .false.
            iplotpartoftypemulti(itypelist(1:nvalues),i) = .true.

         endif
      else
         !
         !--if ntypes < 2 always use the (only) particle type
         !
         iusealltypesmulti(i) = .true.
      endif
   enddo

   if (isamex .and. .not.anycoordplot) then
      imultisamepanel = .false.
      !call prompt('plot all plots in same panel? (default is different panels)',imultisamepanel)
   else
      imultisamepanel = .false.
   endif

   if (nyplotmulti.eq.1 .or. imultisamepanel) then
      nacross = 1
      ndown = 1
      print*,'setting nacross,ndown = ',nacross,ndown
   elseif (mod(nacross*ndown,nyplotmulti).ne.0) then
      !--guess nacross,ndown based on largest factor
      ifac = nyplotmulti/2
      do while (mod(nyplotmulti,ifac).ne.0 .and. ifac.gt.1)
         ifac = ifac - 1
      end do
      if (ifac.le.1) then
         nacross = nyplotmulti/2
      else
         nacross = ifac
      endif
      if (nacross.le.0) nacross = 1
      ndown = nyplotmulti/nacross
      print*,'setting nacross,ndown = ',nacross,ndown
   else
      print*,'nacross = ',nacross,' ndown = ',ndown
   endif

   return
   end subroutine options_multiplot
end subroutine menu

!----------------------------------------------
! utility function which determines whether
! or not rendering is allowed or not
!----------------------------------------------
logical function allowrendering(iplotx,iploty,xsec)
 use labels,          only:ih,irho !,ipmass
 use multiplot,       only:itrans
 use settings_data,   only:ndataplots,icoords,icoordsnew
 use settings_render, only:icolour_particles
 implicit none
 integer, intent(in) :: iplotx,iploty
 logical, intent(in), optional :: xsec
 integer :: itransx,itransy
 logical :: is_xsec
 
 if (present(xsec)) then
    is_xsec = xsec
 else
    is_xsec = .true.
 endif
 itransx = 0
 itransy = 0
 if (iplotx.gt.0) itransx = itrans(iplotx)
 if (iploty.gt.0) itransy = itrans(iploty)
!
!--work out whether rendering is allowed based on presence of rho, h & m in data read
!  also must be in base coordinate system and no transformations applied
!
 if ((ih.gt.0 .and. ih.le.ndataplots) &
    .and.(irho.gt.0 .and. irho.le.ndataplots) &
    .and.(icoords.eq.icoordsnew .or. .not.is_xsec) &
    .and.(itransx.eq.0 .and. itransy.eq.0)) then

    allowrendering = .true.
 else
    allowrendering = .false.
    if (icolour_particles) allowrendering = .true.
 endif

end function allowrendering

!----------------------------------------------
! utility function which sets up the "extra"
! plot columns and returns the total number
! of allowed columns for plotting
!----------------------------------------------
subroutine set_extracols(ncolumns,ncalc,nextra,numplot,ndataplots)
 use params,        only:maxplot
 use labels,        only:ipowerspec,iacplane,isurfdens,itoomre,iutherm,ipdf,label,icolpixmap
 use settings_data, only:ndim,icoordsnew,ivegotdata,debugmode
 use settings_part, only:iexact
 use system_utils,  only:lenvironment
 use write_pixmap,  only:ireadpixmap
 implicit none
 integer, intent(in)    :: ncolumns
 integer, intent(inout) :: ncalc
 integer, intent(out)   :: nextra,numplot,ndataplots

 !
 !-add extra columns (but not if nothing read from file)
 !
 if (ncolumns.gt.0) then
    nextra = 0
    ipowerspec = 0
    iacplane = 0
    isurfdens = 0
    itoomre = 0
    if (ndim.eq.3 .and. icoordsnew.eq.2 .or. icoordsnew.eq.3) then
       nextra = nextra + 1
       isurfdens = ncolumns + ncalc + nextra
       label(isurfdens) = 'Surface density'
       if (iutherm.gt.0 .and. iutherm.le.ncolumns) then
          nextra = nextra + 1
          itoomre = ncolumns + ncalc + nextra
          label(itoomre) = 'Toomre Q parameter'
       endif
    endif
    if (ndim.eq.3 .and. lenvironment('SPLASH_TURB')) then  !--Probability Density Function
       nextra = nextra + 1
       ipdf = ncolumns + ncalc + nextra
       label(ipdf) = 'PDF'
    endif

    if (ndim.le.1 .and. lenvironment('SPLASH_TURB')) then !! .or. ndim.eq.3) then ! if 1D or no coord data (then prompts for which x)
       nextra = nextra + 1      ! one extra plot = power spectrum
       ipowerspec = ncolumns + ncalc + nextra
       label(ipowerspec) = '1D power spectrum'
    else
       ipowerspec = 0
    endif
    if (iexact.eq.6) then       ! toy star plot a-c plane
       nextra = nextra + 1
       iacplane = ncolumns + ncalc + nextra
       label(iacplane) = 'a-c plane'
    else
       iacplane = 0
    endif
    !nextra = nextra + 1
    !label(ncolumns+ncalc+nextra) = 'gwaves'
    if (ndim.ge.2) then
       if (ireadpixmap) then
          nextra = nextra + 1
          icolpixmap = ncolumns + ncalc + nextra
          label(icolpixmap) = '2D pixel map'
       endif
    endif
 endif
!
!--now that we know nextra, set the total number of allowed plots (numplot).
!
 if (ivegotdata) then
    numplot = ncolumns + ncalc + nextra
    if (numplot.gt.maxplot) then
       print "(a,i3,a)",' ERROR: total number of columns = ',numplot,' is greater '
       print "(a,i3,a)",'        than the current allowed maximum (',maxplot,').'
       print "(a)",'        This is set by the parameter "maxplot" in the params module'
       print "(a)",'        in the file globaldata.f90 -- edit this and recompile splash'
       print "(a)",'        (or email me if this happens to increase the default limit)'
       stop
    endif
    ndataplots = ncolumns + ncalc
 else
    numplot = 0
    ndataplots = 0
    ncalc = 0
 endif

 if (debugmode) print*,'DEBUG: numplot = ',numplot, ' ncalc = ',ncalc,' ndataplots = ',ndataplots

 return
end subroutine set_extracols

!----------------------------------------
! instant multiplot setup from main menu
!----------------------------------------
subroutine set_instant_multiplot(string,ipicky,ipickx,numplot,nmulti,multiplotx,multiploty,nx,ny)
 use params, only:maxplot
 use prompting, only:prompt
 character(len=*), intent(in) :: string
 integer, intent(in) :: numplot
 integer, intent(inout) :: ipicky,ipickx
 integer, intent(inout) :: nmulti,nx,ny
 integer, intent(inout) :: multiplotx(:),multiploty(:)
 integer :: ipickarr(maxplot),ierr,i
 
 ipickarr = 0
 read(string,*,iostat=ierr) ipickarr
 i = 1
 do while (i < size(ipickarr) .and. ipickarr(i) /= 0 .and. ipickarr(i) <= numplot)
    i = i + 1
 enddo
 if (i > 2) then
    nmulti = i-1
    !--make sure nmulti matches the number of panels on the page
    if (nmulti /= nx*ny) then
       nx = 1
       ny = nmulti
    endif
    multiploty(1:nmulti) = ipickarr(1:nmulti)
    ipicky = numplot + 1
    if (ipickx==0) ipickx = 1 ! do not allow zero as default
    call prompt(' (x axis) ',ipickx)
    multiplotx(1:nmulti) = ipickx
 endif

end subroutine set_instant_multiplot

!--------------------
! print menu header
!--------------------
subroutine print_header
 integer :: v(8),i
 integer, parameter :: m(48) = (/32,68,111,110,39,116,32,102,111,114,103,101,116,32,116,111,&
                      32,115,101,110,100,32,68,97,110,105,101,108,32,97,32,98,&
                      105,114,116,104,100,97,121,32,109,101,115,115,97,103,101,33/)
 integer, parameter :: c(49) = (/32,120,111,120,111,120,111,120,111,32,77,101,114,114,121,32,&
                      67,104,114,105,115,116,109,97,115,32,102,114,111,109,32,115,&
                      112,108,97,115,104,33,32,120,111,120,111,120,111,120,111,120,111/)
 integer, parameter :: d(49) = (/32,79,111,79,111,79,111,79,32,83,80,76,65,83,72,32,119, &
                      105,115,104,101,115,32,121,111,117,32,97,32,118,101,114,121,32,104,&
                      97,112,112,121,32,33,32,79,111,79,111,79,111,79/)
 call date_and_time(values=v)
 if (v(2)==m(1)/4 .and. v(3)==v(2)-2) then
    print "(/,48(a))",(achar(m(i)),i=1,48)
 elseif (v(2)==(m(1)-20) .and. v(3)>20) then
    print "(/,49(a))",(achar(c(i)),i=1,49)
 elseif (v(2)==nint(0.6) .and. v(3)==d(2)/79) then
    print "(/,40(a),i4,9(a))",(achar(d(i)),i=1,40),v(1),(achar(d(i)),i=41,49)
 else
    print "(/a)",' You may choose from a delectable sample of plots'
 endif

end subroutine print_header

end module mainmenu
