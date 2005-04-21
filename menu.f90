!--------------------
!    MAIN MENU
!--------------------
module mainmenu
 implicit none
 public :: menu
 
 private

contains

subroutine menu
  use filenames
  use labels, only: label,labelvec,iamvec,iacplane,ipowerspec,ih,irho,ipmass
  use options_data, only:submenu_data
  use settings_data, only:ndim,numplot,ndataplots,nextra,ncalc,ivegotdata, &
                     icoords,buffer_data,ncolumns
  use settings_limits, only:submenu_limits
  use settings_part, only:submenu_particleplots,iexact,icoordsnew
  use settings_page, only:submenu_page,interactive
  use settings_render, only:submenu_render
  use settings_vecplot, only:submenu_vecplot,iplotpartvec
  use settings_xsecrot, only:submenu_xsecrotate
  use multiplot
  use prompting
  use transforms
  use defaults, only:defaults_write
  use geometry, only:labelcoord
  use getdata, only:get_data
  use timestepping
  implicit none
  integer :: i,icol,ihalf,iadjust,index,ierr
  integer :: ipicky,ipickx,irender,ivecplot
  integer :: iamvecprev, ivecplottemp
  character(len=2) :: ioption
  character(len=50) :: vecprompt
  logical :: iAllowRendering

  irender = 0
  ivecplot = 0
  ipowerspec = 0
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
  nextra = 0
  if (ndim.eq.1) then
     nextra = 1      ! one extra plot = power spectrum
     ipowerspec = ncolumns + ncalc + 1
     label(ipowerspec) = '1D power spectrum'
  endif
  if (iexact.eq.4) then       ! toy star plot a-c plane
     nextra = nextra + 1
     iacplane = ncolumns + ncalc + nextra
     label(iacplane) = 'a-c plane'
  endif

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
     label(1:ndim) = labelcoord(1:ndim,icoordsnew)
     do i=1,numplot
        if (iamvec(i).ne.0) then
           label(i) = trim(labelvec(iamvec(i)))//'\d'//labelcoord(i-iamvec(i)+1,icoordsnew)
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
  
11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(55('-'))
13 format(1x,i2,')',1x,a)
18 format(1x,i2,')',1x,a,'[ ',i2,' ]',5x,a2,') ',a)

!
!--options 
! 
  print 12
  print "(a)",' d(ata) i(nteractive) p(age) o(pts) l(imits) h(elp)'
  print "(a)",' r(ender) v(ector) x(sec/rotate) s(ave) q(uit)' 
  print 12

!
!--prompt user for selection
!
9901 continue
  write(*,"(a)",ADVANCE='NO') 'Please enter your selection now (y axis or option):'
  read(*,*,ERR=9901) ioption 

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
           ! (apologies for the gratuitous use of goto - I learnt BASIC as a kid)
           if (ipickx.gt.numplot .or. ipickx.le.0) goto 9901
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
              ivecplottemp = -1
              do while(.not.any(iamvec(1:numplot).eq.ivecplottemp).and.ivecplottemp.ne.0)
                 ivecplottemp = ivecplot
                 call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
                 if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
                    print "(a)",'Error, value not in list' 
                 endif
              enddo
              ivecplot = ivecplottemp
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
     select case(adjustl(ioption))
!------------------------------------------------------------------------
     case('m','M')
        call options_multiplot
     case('?m','?M')
        print "(/a)",' Sets up plotting of (m)ultiple quantities per timestep'
!------------------------------------------------------------------------
     case('d','D')
        call submenu_data
     case('?d','?D')
        print "(/a)",' This submenu sets options relating to the (d)ata read'
!------------------------------------------------------------------------
     case('i','I')
        interactive = .not.interactive
        print*,' Interactive mode = ',interactive
     case('?i','?I')
        print "(/a)",' This option turns (i)nteractive mode on/off'
!------------------------------------------------------------------------
     case('p','P')
        call submenu_page
     case('?p','?P')
        print "(/a)",' This submenu sets (p)age setup options'
!------------------------------------------------------------------------
     case('o','O')
        call submenu_particleplots
     case('?o','?O')
        print "(/a)",' This submenu sets particle plot (o)ptions'
!------------------------------------------------------------------------
     case('r','R')
        call submenu_render
     case('?r','?R')
        print "(/a)",' This submenu sets (r)endering options'
!------------------------------------------------------------------------
     case('v','V')
        call submenu_vecplot
     case('?v','?V')
        print "(/a)",' This submenu sets (v)ector plotting options'
!------------------------------------------------------------------------
     case('x','X')
        call submenu_xsecrotate
     case('?x','?X')
        print "(/a)",' This submenu sets cross section and rotation options'
!------------------------------------------------------------------------
     case('l','L')
        call submenu_limits
     case('?l','?L')
        print "(/a)",' This submenu sets options relating to the plot limits'
        call submenu_limits(help=.true.)
!------------------------------------------------------------------------
     case('s','S')
        call defaults_write
     case('?s','?S')
        print "(6(/a))",' The (s)ave option saves the default options to a ', &
                    ' file called `defaults'' in the current directory which', &
                    ' is read automatically upon the next invocation of supersphplot.',&
                    ' This file uses namelist formatting and may be edited ', &
                    ' manually prior to startup if so desired. This is quite',&
                    ' useful for setting multiplots with many plots per page'
!------------------------------------------------------------------------
     case('h','H')
        print "(4(/a))",' For help on any menu item type a question mark ',&
                 ' preceding the appropriate letter.',&
                 ' ie. ?d gives help on the data submenu',&
                 '     ?x gives help on cross section/rotation', &
                 '     ...you get the idea'
     case('?h','?H')
        print "(2(/a))",' You mean to tell me you need help on help???',&
                      ' You need help.'
!------------------------------------------------------------------------
     case('q','Q')
        return
     case('?q','?Q')
        print "(3(/a))",' (q)uit, unsurprisingly, quits. Typing a number greater',&
                    ' than the number of data columns also exits the program',&
                    ' (e.g. I often simply type 99 to exit).'
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
   !use multiplot
   use settings_page, only: nacross, ndown
   use settings_render, only: iplotcont_nomulti
   use settings_xsecrot
   !use prompting
   implicit none
   logical :: iansx, iansy, ichange
   
   call prompt('Enter number of plots per timestep:',nyplotmulti,1,numplot)
   nacross = nyplotmulti/2
   if (nacross.eq.0) nacross = 1
   ndown = nyplotmulti/nacross
   print*,'setting nacross,ndown = ',nacross,ndown 
   iansx = .true.
   call prompt('Same x axis for all?',iansx)
   if (iansx) then
      call prompt('Enter x axis for all plots',multiplotx(1),1,numplot)
      multiplotx(2:nyplotmulti) = multiplotx(1)        
   endif
   iansy = .false.
   if (ndim.ge.2) call prompt('Same y axis for all?',iansy)
   if (iansy) then
      call prompt('Enter y axis for all plots',multiploty(1),1,numplot)
      multiploty(2:nyplotmulti) = multiploty(1)
   endif

   do i=1,nyplotmulti
      print*,'Plot number ',i,':'
      if (.not.iansy .or. multiploty(i).gt.ndataplots .or. multiploty(i).le.0) then
         call prompt(' y axis ',multiploty(i),1,numplot)
      endif
      if (.not.iansx.and.multiploty(i).le.ndataplots) then
         call prompt(' x axis ',multiplotx(i),1,numplot)
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
      
      if ((multiplotx(i).le.ndim).and.(multiploty(i).le.ndim) .and.iAllowRendering) then
         call prompt('(render) (0=none)',irendermulti(i),0,numplot)
         if (irendermulti(i).ne.0) then
            ichange = .false.
            call prompt(' change rendering options for this plot? ',ichange)
            if (ichange) then
               call prompt('plot contours? ',iplotcontmulti(i))
               if (ndim.ge.2) then
                  call prompt(' cross section (no=projection)? ',x_secmulti(i))
                  if (x_secmulti(i)) then
                     call prompt('enter co-ordinate location of cross section slice',xsecposmulti(i))
                  endif
               endif
            elseif (i.eq.1) then
               print*,'copying options from rendering settings'
               iplotcontmulti(i) = iplotcont_nomulti
               x_secmulti(i) = xsec_nomulti
               xsecposmulti(i) = xsecpos_nomulti
            else  
               print*,'using same rendering options as plot 1'       
               iplotcontmulti(i) = iplotcontmulti(1)
               x_secmulti(i) = x_secmulti(1)
               xsecposmulti(i) = xsecposmulti(1)
            endif
         endif
         ivecplottemp = -1
         do while(.not.any(iamvec(1:numplot).eq.ivecplottemp).and.ivecplottemp.ne.0)
            ivecplottemp = ivecplot
            call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
            if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
               print "(a)",'Error, value not in list' 
            endif
         enddo
         ivecplotmulti(i) = ivecplottemp
         if (ivecplotmulti(i).gt.0 .and. irendermulti(i).eq.0) then
            call prompt('plot particles?',iplotpartvec)
         endif
      else
         irendermulti(i) = 0
         ivecplotmulti(i) = 0
      endif
   enddo
   
   return
   end subroutine options_multiplot
end subroutine menu

end module mainmenu
