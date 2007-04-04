!--------------------
!    MAIN MENU
!--------------------
module mainmenu
 implicit none
 public :: menu
 
 private

contains

subroutine menu
  use filenames, only:defaultsfile,limitsfile,animfile
  use labels, only:label,labelvec,iamvec,iacplane,ipowerspec,ih,irho,ipmass
  use limits, only:write_limits
  use options_data, only:submenu_data
  use settings_data, only:ndim,numplot,ndataplots,nextra,ncalc,ivegotdata, &
                     icoords,buffer_data,ncolumns,iRescale
  use settings_limits, only:submenu_limits
  use settings_part, only:submenu_particleplots,iexact,icoordsnew
  use settings_page, only:submenu_page,interactive
  use settings_render, only:submenu_render
  use settings_vecplot, only:submenu_vecplot,iplotpartvec
  use settings_xsecrot, only:submenu_xsecrotate,write_animfile
  use settings_units, only:unitslabel
  use multiplot
  use prompting, only:prompt
  use transforms, only:transform_label
  use defaults, only:defaults_write
  use geometry, only:labelcoord
  use getdata, only:get_data
  use timestepping
  implicit none
  integer :: i,icol,ihalf,iadjust,index,ierr
  integer :: ipicky,ipickx,irender,ivecplot
  integer :: iamvecprev, ivecplottemp,ichoose
  character(len=2) :: ioption
  character(len=100) :: vecprompt
  logical :: iAllowRendering

  irender = 0
  ivecplot = 0
  ipickx = 1
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
  nextra = 0
  ipowerspec = 0
  iacplane = 0
  if (ndim.le.1) then ! if 1D or no coord data (then prompts for which x)
     nextra = 1      ! one extra plot = power spectrum
     ipowerspec = ncolumns + ncalc + 1
     label(ipowerspec) = '1D power spectrum'
  else
     ipowerspec = 0
  endif
  if (iexact.eq.4) then       ! toy star plot a-c plane
     nextra = nextra + 1
     iacplane = ncolumns + ncalc + nextra
     label(iacplane) = 'a-c plane'
  else
     iacplane = 0
  endif
  !nextra = nextra + 1
  !label(ncolumns+ncalc+nextra) = 'gwaves'

  if (ivegotdata) then
     numplot = ncolumns + ncalc + nextra
     if (numplot.gt.maxplot) then
        print*,numplot,ncolumns,ncalc,nextra
        stop 'ERROR: numplot > multiplot array limits: reset this in module params'
     endif
     ndataplots = ncolumns + ncalc
  else
     numplot = 0
     ndataplots = 0
     ncalc = 0
  endif
     
!--set coordinate and vector labels (depends on coordinate system)
  if (icoords.ne.0 .or. icoordsnew.ne.0) then
     do i=1,ndim
        label(i) = labelcoord(i,icoordsnew)
        if (iRescale .and. icoords.eq.icoordsnew) then
           label(i) = trim(label(i))//trim(unitslabel(i))
        endif
     enddo
     do i=1,numplot
        if (iamvec(i).ne.0) then
           label(i) = trim(labelvec(iamvec(i)))//'\d'//labelcoord(i-iamvec(i)+1,icoordsnew)
           if (iRescale) then
              label(i) = trim(label(i))//trim(unitslabel(i))
           endif
        endif
     enddo
  endif
!--set contents of the vector plotting prompt
  vecprompt(1:6) = '0=none'
  index = 7
  iamvecprev = 0
  do icol=1,numplot
     if (iamvec(icol).ne.0 .and. iamvec(icol).ne.iamvecprev) then
        iamvecprev = iamvec(icol)
        if (iamvec(icol).ge.10) then
           write(vecprompt(index:),"(',',1x,i2,'=',a)") &
                 iamvec(icol),trim(labelvec(icol))
        else
           write(vecprompt(index:),"(',',1x,i1,'=',a)") &        
                 iamvec(icol),trim(labelvec(icol))
        endif
        index = len_trim(vecprompt) + 1
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
     print "(/a)",' You may choose from a delectable sample of plots '
     print 12
     ihalf = numplot/2                ! print in two columns
     iadjust = mod(numplot,2)
     print 11, (i,transform_label(label(i),itrans(i)), &
          ihalf + i + iadjust, transform_label(label(ihalf + i + iadjust), &
          itrans(ihalf+i+iadjust)),i=1,ihalf)
     if (iadjust.ne.0) then
        print 13, ihalf + iadjust,transform_label(label(ihalf + iadjust), &
             itrans(ihalf+iadjust))
     endif 
!
!--multiplot
!  
     print 12
     print 18,numplot+1,'multiplot ',nyplotmulti,'m','set multiplot '
  else
!
!--if no data
!
     print "(/a)",' No data: You may choose from the options below '
  endif
  
11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a20)
12 format(55('-'))
13 format(1x,i2,')',1x,a20)
18 format(1x,i2,')',1x,a,'[ ',i2,' ]',5x,a2,') ',a)

!
!--options 
! 
  print 12
  print "(a)",' d(ata) i(nteractive) p(age) o(pts) l(imits) h(elp)'
  print "(a)",' r(ender) v(ector) x(sec/rotate) s,S(ave) q(uit)' 
  print 12

!
!--prompt user for selection
!
  write(*,"(a)",ADVANCE='NO') 'Please enter your selection now (y axis or option):'
  read(*,*,iostat=ierr) ioption
  if (ierr < 0) stop 'reached end of input' ! end of input (e.g. in script)
  if (ierr > 0) stop 'error reading input' 

!------------------------------------------------------------
!  if input is an integer and within range, plot data
!------------------------------------------------------------
  read(ioption,*,iostat=ierr) ipicky
  if (ierr /= 0) ipicky = -1

  if (ipicky.gt.0 .and. ipicky.le.numplot+1) then
     
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
        if (ipicky.le.(numplot-nextra)) then
           if (ipickx.eq.0) ipickx = 1 ! do not allow zero as default
           call prompt(' (x axis) ',ipickx)
           !--go back to y prompt if out of range
           if (ipickx.gt.numplot .or. ipickx.le.0) cycle menuloop
           !
           !--work out whether rendering is allowed based on presence of rho, h & m in data read
           !  also must be in base coordinate system and no transformations applied
           !
           iAllowRendering = (ih.gt.0 .and. ih.le.ndataplots) &
                        .and.(irho.gt.0 .and. irho.le.ndataplots) &
                        .and.(ipmass.gt.0 .and. ipmass.le.ndataplots)  &
                        .and.(icoords.eq.icoordsnew) &
                        .and.(itrans(ipickx).eq.0 .and. itrans(ipicky).eq.0)
           !
           !--prompt for render and vector plots 
           ! -> only allow if in "natural" coord system, otherwise h's would be wrong)
           ! (a future feature might be to interpolate in icoord then translate the pixels
           !  to icoordsnew, or alternatively plot non-cartesian pixel shapes)
           ! -> also do not allow if transformations are applied
           !
           if (ipicky.le.ndim .and. ipickx.le.ndim .and. iAllowRendering) then
              call prompt('(render) (0=none)',irender,0,numplot)
              if (any(iamvec(1:numplot).ne.0)) then
                 ivecplottemp = -1
                 ierr = 1
                 do while(ierr.ne.0 .and. ivecplottemp.ne.0)
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
        endif
        !
        !--call main plotting routine
        !
        call timestep_loop(ipicky,ipickx,irender,ivecplot)
     endif
!------------------------------------------------------------------------
!  if input is an integer > numplot+1, quit
!------------------------------------------------------------------------
  elseif (ipicky.gt.numplot+1) then
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
        print*,' Interactive mode = ',interactive
!------------------------------------------------------------------------
!+ This submenu sets (p)age setup options
     case('p','P')
        call submenu_page(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets particle plot (o)ptions
     case('o','O')
        call submenu_particleplots(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets (r)endering options
     case('r','R')
        call submenu_render(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets (v)ector plotting options
     case('v','V')
        call submenu_vecplot(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets cross section and rotation options
     case('x','X')
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
        call defaults_write(defaultsfile)
     case('S')
        call defaults_write(defaultsfile)
        call write_limits(limitsfile)
        call write_animfile(animfile)
!------------------------------------------------------------------------
!+ Slightly obsolete: prints whatever help may be helpful
     case('h','H')
        print "(10(/a))",' Hint: menu items can be shortcut by typing, e.g. o2 for ',&
                 ' the o)ptions menu, item 2.',&
                 '   ', &
                 ' for detailed help, consult the user guide',&
                 '  (splash/docs/splash.pdf ',&
                 '   or http://www.astro.ex.ac.uk/people/dprice/splash/userguide/)', &
                 ' and/or the online FAQ. If you''re really stuck, email me! '
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
   use settings_page, only: nacross, ndown
   use settings_render, only: iplotcont_nomulti
   implicit none
   integer :: ifac,ierr
   logical :: isamex, isamey, icoordplot
   
   call prompt('Enter number of plots per timestep:',nyplotmulti,1,numplot)
   if (nyplotmulti.eq.1) then
      nacross = 1
      ndown = 1
   else
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
   endif
   print*,'setting nacross,ndown = ',nacross,ndown 
   isamex = all(multiplotx(1:nyplotmulti).eq.multiplotx(1))
   call prompt('Same x axis for all?',isamex)
   if (isamex) then
      call prompt('Enter x axis for all plots',multiplotx(1),1,numplot)
      multiplotx(2:nyplotmulti) = multiplotx(1)
   endif
   isamey = all(multiploty(1:nyplotmulti).eq.multiploty(1))
   if (ndim.ge.2) call prompt('Same y axis for all?',isamey)
   if (isamey) then
      call prompt('Enter y axis for all plots',multiploty(1),1,numplot)
      multiploty(2:nyplotmulti) = multiploty(1)
   endif

   do i=1,nyplotmulti
      print*,'-------------- Plot number ',i,' --------------'
      if (.not.isamey .or. multiploty(i).gt.ndataplots .or. multiploty(i).le.0) then
         call prompt(' y axis ',multiploty(i),1,numplot)
      endif
      if (.not.isamex) then
         if (multiploty(i).le.ndataplots) then
            call prompt(' x axis ',multiplotx(i),1,ndataplots)
         else
            multiplotx(i) = multiploty(i)
         endif
      endif
      !
      !--work out whether rendering is allowed based on presence of rho, h & m in data read
      !  also must be in base coordinate system and no transformations applied
      !
      iAllowRendering = (ih.gt.0 .and. ih.le.ndataplots) &
                   .and.(irho.gt.0 .and. irho.le.ndataplots) &
                   .and.(ipmass.gt.0 .and. ipmass.le.ndataplots)  &
                   .and.(icoords.eq.icoordsnew) &
                   .and.(itrans(multiplotx(i)).eq.0 .and. itrans(multiploty(i)).eq.0)
      
      icoordplot = (multiplotx(i).le.ndim .and. multiploty(i).le.ndim)
      
      if (icoordplot .and.iAllowRendering) then
         call prompt('(render) (0=none)',irendermulti(i),0,numplot)
         iplotcontmulti(i) = iplotcont_nomulti

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
         irendermulti(i) = 0
         ivecplotmulti(i) = 0
      endif

      if (icoordplot .and. ndim.ge.2) then
         call prompt(' is this a cross section (no=projection)? ',x_secmulti(i))
         if (x_secmulti(i)) then
            call prompt('enter co-ordinate location of cross section slice',xsecposmulti(i))
         endif
      endif
      
   enddo
   
   return
   end subroutine options_multiplot
end subroutine menu

end module mainmenu
