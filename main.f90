!
! This subroutine drives the main plotting loop
!
subroutine main(ipicky,ipickx,irender)
  use params
  use exact_params
  use filenames
  use labels
  use multiplot
  use particle_data
  use prompting
  use settings   
  implicit none
  integer, intent(in) :: ipicky, ipickx, irender

  integer :: i,j,k,n
  integer :: iplotx,iploty,ivecx,ivecy
  integer :: nyplot,nyplots      
  integer :: npart1
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: ivecplot,npix,npixvec,npixyvec,ncontours
  integer :: irenderprev, istepprev
  integer :: isizex,isizey      ! for sending datpix to transform
  integer :: nsink,nsinkstart,nsinkend,nghoststart,nghostend
  integer :: ishk,int_from_string
  integer :: igrid, ngrid

  character(len=8) :: string     ! used in pgplot calls
  real, dimension(2,maxpart) :: vecplot
  real, dimension(maxpart) :: xplot,yplot
  real, dimension(:), allocatable :: datpix1D, xgrid
  real, dimension(:,:), allocatable :: datpix
  real, dimension(:,:,:), allocatable :: datpix3D
  real :: xmin,xmax,ymin,ymax,zmin,zmax
  real :: vecmax,rendermin,rendermax
  real :: xsecmin,xsecmax,dxsec,xsecpos
  real :: pixwidth
  real :: charheight
  real :: dxgrid
  real :: xpt1,ypt1,xpt2,ypt2

  logical :: iplotcont,iplotpartvec,x_sec
  logical :: log,use_backgnd_color_vecplot

  character(len=60) :: title,titlex
  character(len=1) :: logx,logy
  character(len=20) :: labelx,labely,labelrender
  character(len=25) :: transform_label

  !------------------------------------------------------------------------
  ! initialisations
  !------------------------------------------------------------------------


  x_sec = xsec_nomulti
  ivecplot = ivecplot_nomulti
  iplotcont = iplotcont_nomulti
  iplotpartvec = iplotpartvec_nomulti
  npix = npix_nomulti
  npixvec = npixvec_nomulti
  ncontours = ncontours_nomulti
  iplotcont = iplotcont_nomulti
  use_backgnd_color_vecplot = backgnd_vec_nomulti

  if (ndim.ne.3) x_sec = .false.

  !!--set current plot to first in multiplot array if doing multiplot
  nxsec = 1
  if (ipicky.eq.numplot+1) then   ! multiplot
     imulti=.true.
     if (any(multiplotx(1:nyplotmulti).ne.multiplotx(1))) then
        isamexaxis = .false.
     endif
     iplotx = multiplotx(1)
     iploty = multiploty(1)
     nyplots = nyplotmulti
  else
     nyplots = 1 
     iploty = ipicky
     iplotx = ipickx
  endif
  !!--if doing multiplot can only take a single cross section slice      
  if (imulti) then
     flythru = .false.
     nxsec = 1
  endif

  !------------------------------------------------------------------------
  ! co-ordinate plot initialisation

  if (ipicky.le.ndim .and. ipickx.le.ndim) then

     !!--work out coordinate that is not being plotted	 
     ixsec = 0
     do j=1,ndim
        if ((j.ne.ipickx).and.(j.ne.ipicky)) ixsec = j
     enddo
     if (ixsec.eq.0) x_sec = .false.   ! ie can only have x_sec in 3D

     !!--if series of cross sections (flythru), set position of first one      
     if (x_sec.and.flythru) then
        print 32,label(ixsec)
32      format('enter number of ',a1,' cross-section slices')
        read*,nxsec
        !!--dxsec is the distance between slices	    
        dxsec = (lim(ixsec,2)-lim(ixsec,1))/float(nxsec)
        xsecpos = lim(ixsec,1) - 0.5*dxsec
        xsecpos_nomulti = xsecpos

        !!--if single cross-section, read position of cross-section slice

     elseif (x_sec.and.iplotpart.and.irender.le.ndim) then
        print 33,label(ixsec)
33      format(' enter ',a1,' position for cross section slice:')
        read*,xsecpos
        if (xsecpos.gt.lim(ixsec,2).or.xsecpos.lt.lim(ixsec,1)) then
           print*,'warning: cross section outside data limits' 
        endif
        xsecpos_nomulti = xsecpos
     endif

     !!--set title of plot

     if ((.not.imulti).and.(nacross*ndown.eq.1)) then
        !        if (x_sec) then
        !           titlex = 'cross-section'
        !        else
        !           titlex = 'projection'   
        !        endif	
        !        titlex = trim(label(ipickx))//trim(label(ipicky))//' '//titlex	          
        !        
        !        if (irender.gt.ndim) then
        !           titlex = trim(titlex)//' - '//trim(label(irender))//' rendering'
        !        endif
        titlex = ' '
        if (ivecplot.eq.1) titlex = trim(' velocity map: '//titlex)
        if (ivecplot.eq.2) titlex = trim(' magnetic field map: '//titlex)
     else
        titlex = ' '
     endif

     !!--initialise pgplot     
     call pgbegin(0,'?',nacross,ndown)
     !
     !--set colour table
     !
     if (((irender.gt.ndim).or.any(irendermulti(1:nyplots).gt.ndim)) &
          .and.(icolours.gt.0)) then
        call colour_set(icolours)
     endif

     !!------------------------------------------------------------------------      
     ! non- co-ordinate plot initialisations
     !
  else
     !!--prompt for options if plotting power spectrum      
     if (ipicky.eq.ipowerspec) call options_powerspec
     !!--no title if more than one plot on the page
     if ((nacross.gt.1).or.(ndown.gt.1)) title = '          '

     call pgbegin(0,'?',nacross,ndown)    !  initialise PGPLOT

  endif
  !!------------------------------------------------------------------------
  ! general initialisations

  !!--set paper size
  if (nacross.eq.2 .and. ndown.eq.1) then
     call pgpaper(11.7,0.5/sqrt(2.))
  elseif (ipapersize.gt.0 .and. papersizex.gt.0.0 .and. aspectratio.gt.0.0 ) then
     call pgpaper(papersizex,aspectratio)
  endif
  !!--turn on/off page prompting
  if (animate) call pgask(.false.)

  !
  !--if plotting ghost particles, set ntotplot = ntot, else ntot=npart
  !
  ntotplot(:) = npart(:)
  if (iplotghost) ntotplot = npart(:) + nghost(:)
  if (iplotsink) ntotplot = ntot(:)
  !
  !--set fill style for circle plots
  !      
  if (plotcirc) call pgsfs(2)
  !
  !--increase character size depending on the number of graphs on the page
  !
  !      print*,' enter character height '
  !      read*,charheight
  charheight = 1.0
  if ((ndown*nacross).gt.1) charheight = 2.0
  !      charheight = 0.5*(nacross+ndown)

  !------------------------------------------------------------------------      
  ! loop over timesteps 
  !------------------------------------------------------------------------            
  over_timesteps: do i=nstart,n_end,nfreq

     npart1 = npart(i) + 1
     irenderprev = 0
     istepprev = 0  
     !-------------------------------------
     ! loop over plots per timestep
     !-------------------------------------
     over_plots: do nyplot=1,nyplots
        !--make sure character height is set correctly
        call pgsch(charheight)
        !--for consecutive plots (ie. if not multi but nyplots > 1 plots consecutive numbers)	     
        iploty = ipicky + nyplot - 1
        !--set current x, y plot from multiplot array
        if (imulti) then
           iploty = multiploty(nyplot)
           iplotx = multiplotx(nyplot)
        endif
        !--------------------------------------------------------------
        !  copy from main dat array into xplot, yplot 
        !  apply transformations (log, 1/x, etc) if appropriate
        !--------------------------------------------------------------
        if (iploty.le.numplot-nextra .and. iplotx.le.numplot-nextra) then
           call transform(dat(iplotx,:,i),xplot,itrans(iplotx),maxpart)
           call transform(dat(iploty,:,i),yplot,itrans(iploty),maxpart)
           !--set axis labels, applying transformation if appropriate
           labelx = transform_label(label(iplotx),itrans(iplotx))
           labely = transform_label(label(iploty),itrans(iploty))
           !--set x,y plot limits, applying transformation if appropriate
           call transform_limits(lim(iplotx,1),lim(iplotx,2),  &
                xmin,xmax,itrans(iplotx))
           call transform_limits(lim(iploty,1),lim(iploty,2),  &
                ymin,ymax,itrans(iplotx))
           !--work out whether to use log axes - this is for the call to pgbox
           logx = ' '
           logy = ' '
           !if (itrans(iplotx).eq.1) logx = 'l'
           !if (itrans(iploty).eq.1) logy = 'l'

           !--write username, date on plot
           !         if (nacross.le.2.and.ndown.le.2) call pgiden

           !--adjust plot limits if adaptive plot limits set
           if ((ipagechange.and.iadapt).and.(iplotx.le.ndataplots) &
                .and.(iploty.le.ndataplots)) then
              xmin = minval(xplot(1:ntotplot(i)))
              xmax = maxval(xplot(1:ntotplot(i)))*scalemax
              ymin = minval(yplot(1:ntotplot(i)))
              ymax = maxval(yplot(1:ntotplot(i)))*scalemax
           endif

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! plots with co-ordinates as x and y axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ((iploty.le.ndim).and.(iplotx.le.ndim)) then

           !!--set rendering options equal to settings in multiplot	 
           if (imulti) then
              irenderplot = irendermulti(nyplot)      
              iplotcont = iplotcontmulti(nyplot)
              ncontours = ncontoursmulti(nyplot)
              ivecplot = ivecplotmulti(nyplot)
              npix = npixmulti(nyplot)
              npixvec = npixvecmulti(nyplot)
              iplotpartvec = iplotpartvecmulti(nyplot)
              x_sec = x_secmulti(nyplot)
              xsecpos = xsecposmulti(nyplot)
              use_backgnd_color_vecplot = backgnd_vec_multi(nyplot)
           else
              irenderplot = irender
              iplotcont = iplotcont_nomulti
              ncontours = ncontours_nomulti
              ivecplot = ivecplot_nomulti
              npix = npix_nomulti
              npixvec = npixvec_nomulti
              iplotpartvec = iplotpartvec_nomulti
              use_backgnd_color_vecplot = backgnd_vec_nomulti
              x_sec = xsec_nomulti
              xsecpos = xsecpos_nomulti
           endif
           npixx = npix

           !!--work out coordinate that is not being plotted	 
           do j=1,ndim
              if ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j
           enddo
           !!           if (ixsec.eq.0) x_sec = .false.   ! ie can only have x_sec in 3D	   

           !------------------------------------------------------------------
           !  rendering setup and interpolation (this is the rendering done
           !  *before* the cross sections are taken, e.g. to 3D grid)
           !------------------------------------------------------------------
           if ((irenderplot.gt.ndim).and. &
                ((ndim.eq.3).or.(ndim.eq.2.and..not.x_sec))) then
              !
              !--interpolate from particles to fixed grid using sph summation
              !		
              !--do not apply any transformations to the co-ordinates
              xmin = lim(iplotx,1)
              xmax = lim(iplotx,2)
              ymin = lim(iploty,1)
              ymax = lim(iploty,2)
              !!--determine number of pixels in rendered image (npix = pixels in x direction)
              pixwidth = (xmax-xmin)/real(npix)
              npixx = int((xmax-xmin)/pixwidth) + 1
              npixy = int((ymax-ymin)/pixwidth) + 1
              isizex = npixx
              isizey = npixy
              print*,'npixx, npixy = ',npixx,npixy
              !!--only need z pixels if working with interpolation to 3D grid
              if ((ndim.ge.3).and.(x_sec.and.nxsec.gt.2)) then
                 zmin = lim(ixsec,1)
                 zmax = lim(ixsec,2)
                 !		   npixz = int((zmax-zmin)/pixwidth) + 1		 
                 !!--number of z pixels is equal to number of cross sections
                 npixz = nxsec
                 print*,'npixz = ',npixz
              endif

              !!--if rendering array is the same as the previous plot, reuse the array		
              if (irenderplot.eq.irenderprev .and. i.eq.istepprev) then
                 print*,'same rendering, using previous array...'
              else
                 if (allocated(datpix)) deallocate(datpix)
                 if (allocated(datpix3D)) deallocate(datpix3D)
                 select case(ndim)
                 case(2)
                    !!--interpolate to 2D grid
                    !!  allocate memory for rendering array
                    if (.not. x_sec) then
                       isizex = npixx
                       isizey = npixy
                       allocate ( datpix(npixx,npixy) )
                       call interpolate2D( &
                            dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i), &
                            dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i), &
                            dat(ih,1:ntot(i),i),dat(irenderplot,1:ntot(i),i), &
                            ntot(i),xmin,ymin,datpix,npixx,npixy,pixwidth)
                    endif
                 case(3)
                    !!--interpolation to 3D grid - then take multiple cross sections/projections
                    !!  do this if taking more than 2 cross sections, otherwise use fast xsec
                    if (x_sec.and.nxsec.gt.2) then
                       !!--allocate memory for 3D rendering array
                       allocate ( datpix3D(npixx,npixy,npixz) )
                       !!--interpolate from particles to 3D grid
                       call interpolate3D( &
                            dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i), &
                            dat(ixsec,1:ntot(i),i),dat(ipmass,1:ntot(i),i),  &
                            dat(irho,1:ntot(i),i),dat(ih,1:ntot(i),i), &
                            dat(irenderplot,1:ntot(i),i), &
                            ntot(i),xmin,ymin,zmin,datpix3D,npixx,npixy,npixz,pixwidth,dxsec)
                    endif
                 end select
              endif

              irenderprev = irenderplot
              istepprev = i
           endif
           !
           !--if vector plot determine whether or not to plot the particles as well
           !
           iplotpart = .true.
           if (ivecplot.gt.0) iplotpart = iplotpartvec
           if (irenderplot.gt.0) iplotpart = .false.

           !
           !%%%%%%%%%%%%%%% loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
           !
           over_cross_sections: do k=1,nxsec

              !------------------------------------------------------------
              ! for multislice cross section (flythru)
              ! increment the position of the current cross section slice
              !-------------------------------------------------------------

              if (x_sec.and.flythru) then
                 xsecpos = xsecpos + dxsec
                 !!--for cross sections of particle plots, need range of co-ordinates in which
                 !!  particles may lie
                 if (iplotpart) then
                    xsecmin = xsecpos-0.5*dxsec
                    xsecmax = xsecpos+0.5*dxsec
                 endif
              endif

              !------------take projections/cross sections through 3D data-----------------!
              if (irenderplot.gt.ndim .and. ndim.eq.3) then

                 !!--allocate memory for 2D rendered array
                 if (allocated(datpix)) deallocate(datpix)
                 allocate ( datpix(npixx,npixy) )

                 !------------------------------------------------------------------------
                 ! if we have rendered to a 3D grid, take cross sections from this array
                 !------------------------------------------------------------------------
                 if (x_sec .and. nxsec.gt.2) then
                    ipixxsec = int((xsecpos-zmin)/dxsec) + 1
                    if (ipixxsec.gt.npixz) ipixxsec = npixz
                    print*,TRIM(label(ixsec)),' = ',xsecpos, &
                         ' cross section, pixel ',ipixxsec
                    datpix = datpix3D(:,:,ipixxsec)    ! slices are in 3rd dimension

                 else
                    !-------------------------------------------------------------------
                    !  or do a fast projection/cross section of 3D data to 2D array
                    !-------------------------------------------------------------------
                    if (x_sec) then
                       !!--do fast cross-section
                       print*,trim(label(ix(ixsec))),' = ',xsecpos,  &
                            ' : fast cross section', xmin,ymin
                       call interpolate3D_fastxsec( &
                            dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i), &
                            dat(ixsec,1:ntot(i),i), &
                            dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),    &
                            dat(ih,1:ntot(i),i),dat(irenderplot,1:ntot(i),i), &
                            ntot(i),xmin,ymin,xsecpos,datpix,npixx,npixy,pixwidth)
                    else
                       !!--do fast projection
                       call interpolate3D_projection( &
                            dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i), &
                            dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),   &
                            dat(ih,1:ntot(i),i), dat(irenderplot,1:ntot(i),i), &
                            ntot(i),xmin,ymin,datpix,npixx,npixy,pixwidth)
                    endif

                 endif ! whether 3D grid or fast renderings

                 !-------------take cross sections through 2D data------------------!
              elseif (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) then
                 !-------------------------------------------------------------------
                 !  or do a fast cross section through 2D data to 1D array
                 !-------------------------------------------------------------------
                 !!--interpolate from 2D data to 1D line
                 !!  line is specified by giving two points, (x1,y1) and (x2,y2)
                 call prompt('enter xmin of cross section line',xpt1,xmin,xmax)
                 call prompt('enter xmax of cross section line',xpt2,xmin,xmax)
                 call prompt('enter ymin of cross section line',ypt1,ymin,ymax)
                 call prompt('enter ymax of cross section line',ypt2,ymin,ymax)
                 isizex = npixx
                 !!--set up 1D grid
                 if (allocated(xgrid)) deallocate(xgrid)
                 allocate ( xgrid(npixx) )
                 xmin = 0.   ! distance (r) along cross section
                 xmax = SQRT((ypt2-ypt1)**2 + (xpt2-xpt1)**2)
                 dxgrid = (xmax-xmin)/REAL(npixx)
                 do igrid = 1,npixx
                    xgrid(igrid) = xmin + igrid*dxgrid - 0.5*dxgrid
                 enddo
                 !!--interpolate to 1D cross section
                 if (allocated(datpix1D)) deallocate(datpix1D)
                 allocate ( datpix1D(npixx) )
                 call interpolate2D_xsec( &
                      dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i), &
                      dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),    &
                      dat(ih,1:ntot(i),i),dat(irenderplot,1:ntot(i),i), &
                      ntot(i),xpt1,ypt1,xpt2,ypt2,datpix1D,npixx)
                 !
                 !--find limits of datpix1D for plotting
                 !  do transformations on rendered array where appropriate
                 !  set these as ymin,ymax and set labels of plot
                 !
                 call transform(datpix1D,datpix1D,itrans(irenderplot),npixx)
                 labely = transform_label(label(irenderplot),itrans(irenderplot))
                 labelx = 'cross section'
                 !!--if adaptive limits, find limits of datpix
                 if (iadapt) then
                    ymin = minval(datpix1D)
                    ymax = maxval(datpix1D)
                    !!--or apply transformations to fixed limits
                 else
                    call transform_limits(lim(irenderplot,1),lim(irenderplot,2), &
                         ymin,ymax,itrans(irenderplot))
                 endif

              endif ! 2 or 3D and rendering
              !-----end of preliminary muff for 2D/3D cross sections/renderings ------------------

              !-----------------------
              ! set up pgplot page
              !-----------------------
              if ((ipagechange).or.((.not.ipagechange).and.(i.eq.nstart))) then
                 !	       call pgenv(lim(iplotx,1),lim(iplotx,2), &
                 !                     lim(iploty,1),lim(iploty,2),1,1)	! 0 for no axes
                 call pgpage
                 if (nacross*ndown.gt.1) then
                    !	         if (imulti) then
                    if (axes) then
                       call pgsvp(0.2,0.8,0.2,0.98)
                    else    ! if no axes use full viewport
                       call pgsvp(0.02,0.98,0.02,0.98)
                    endif
                    !		  else
                    !		     call pgsvp(0.0,1.0,0.0,1.0)
                    !		  endif
                 else
                    if (axes) then
                       call pgsvp(0.1,0.9,0.1,0.9)
                    else      ! if no axes use full viewport
                       call pgsvp(0.02,0.98,0.02,0.98)
                    endif
                 endif
                 if (ndim.eq.2 .and. x_sec) then
                    call pgswin(xmin,xmax,ymin,ymax)
                 else
                    call pgwnad(xmin,xmax,ymin,ymax)   !  pgwnad does equal aspect ratios
                 endif
                 !!--plot axes (log if appropriate)
                 if (axes) then
                    call pgbox('bcnst'//logx,0.0,0,'1bvcnst'//logy,0.0,0)
                 elseif (ivecplot.ne.0) then
                    call pgbox('bc',0.0,0,'bc',0.0,0)  ! draw box only for vector plots
                 endif
              elseif (nyplot.eq.1) then
                 call pgpanl(1,1)
              else
                 call pgpage
              endif

              !---------------------------------
              ! set plot limits and label plot
              !---------------------------------

              !--print plot limits to screen
              print 34, time(i),i
              print*,trim(labely),'min,max = ',ymin,ymax
              print*,trim(labelx),'min,max = ',xmin,xmax
34            format (5('-'),' t = ',f8.4,', dump #',i3,1x,10('-'))
              if (x_sec.and.iplotpart) print 35,label(ixsec),xsecmin,label(ixsec),xsecmax
35            format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

              if (ndim.eq.2 .and. x_sec) then
                 call pgswin(xmin,xmax,ymin,ymax)
              else
                 call pgwnad(xmin,xmax,ymin,ymax)   !  pgwnad does equal aspect ratios
              endif

              !--label plot
              if (axes) then
                 if (((nyplots-nyplot).lt.nacross).or.(.not.isamexaxis)) then
                    call pgmtxt('l',3.0,0.5,1.0,labely)
                    call pglabel(labelx,' ',trim(titlex))
                 else
                    call pgmtxt('l',3.0,0.5,1.0,labely)
                    !	     call pglabel(' ',labely,trim(titlex))
                 endif
              else
                 !--if multiple plots showing contours only, label with a,b,c etc
                 if ((nacross*ndown.gt.1).and.(irenderplot.gt.ndim) &
                      .and.(icolours.eq.0).and.imulti) then
                    select case(nyplot)
                    case(1)
                       call pgmtext('T',-1.5,0.05,0.0,'a)') 
                    case(2)
                       call pgmtext('T',-1.5,0.05,0.0,'b)')
                    case(3)
                       call pgmtext('T',-1.5,0.05,0.0,'c)')
                    case(4)
                       call pgmtext('T',-1.5,0.05,0.0,'d)')
                    case(5)
                       call pgmtext('T',-1.5,0.05,0.0,'e)')
                    case(6)
                       call pgmtext('T',-1.5,0.05,0.0,'f)')
                    case(7)
                       call pgmtext('T',-1.5,0.05,0.0,'g)')
                    case(8)
                       call pgmtext('T',-1.5,0.05,0.0,'h)')
                    end select
                 endif
              endif

              !------------------------------
              ! now actually plot the data
              !------------------------------
              !---------------------------------------------------------------
              ! density/scalar field rendering
              ! having got our 2D array of gridded data (datpix), we plot it
              !---------------------------------------------------------------
              if (irenderplot.gt.ndim .and.    &
                   ((ndim.eq.3).or.(ndim.eq.2.and. .not.x_sec)) ) then
                 !!--do transformations on rendered array  
                 call transform2(datpix,itrans(irenderplot),isizex,isizey)
                 labelrender = label(irenderplot)
                 !!--set label for column density (projection) plots (2268 or 2412 for integral sign)
                 if (ndim.eq.3 .and..not. x_sec) then         
                    labelrender = '\(2268) '//trim(labelrender)//' d'//trim(label(ix(ixsec)))
                 endif
                 labelrender = transform_label(labelrender,itrans(irenderplot))
                 !!--set log axes for call to render
                 log = .false.
                 if (itrans(irenderplot).eq.1) log = .true.

                 !!--if adaptive limits, find limits of rendered array
                 if (iadapt) then
                    rendermin = minval(datpix)
                    rendermax = maxval(datpix)
                    !!--or apply transformations to fixed limits
                 else
                    call transform_limits(lim(irenderplot,1),lim(irenderplot,2), &
                         rendermin,rendermax,itrans(irenderplot))
                 endif
                 !!--print plot limits to screen
                 print*,trim(labelrender),' min, max = ',rendermin,rendermax       
                 !!--call subroutine to actually render the image       
                 call render(datpix,rendermin,rendermax,trim(labelrender),  &
                      npixx,npixy,xmin,ymin,pixwidth,    &
                      icolours,iplotcont,ncontours,log)

              elseif (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) then
                 !---------------------------------------------------------------
                 ! plot 1D cross section through 2D data    
                 !---------------------------------------------------------------

                 !!--plot 1D cross section (contents of datpix)       
                 call pgline(npixx,xgrid,datpix1D) 
              else
                 !-----------------------
                 ! particle plots
                 !-----------------------
                 !
                 !--if particle cross section, plot particles only in a defined coordinate range
                 !
                 if (x_sec.and.iplotpart) then
                    do j=1,npart(i)
                       if ((dat(ixsec,j,i).lt.xsecmax) &
                            .and.(dat(ixsec,j,i).gt.xsecmin)) then
                          call pgpt(1,xplot(j),yplot(j),imark)
                          !!--plot circles of interaction
                          if (plotcirc) then
                             call pgcirc(xplot(j),yplot(j),2.*dat(ih,j,i))
                          endif
                          if (ilabelpart) then
                             !!--plot particle labels
                             call pgnumb(j,0,1,string,nc)
                             call pgsch(0.5*charheight)
                             call pgtext(xplot(j),yplot(j),string(1:nc))
                             call pgsch(charheight)
                          endif! ilabelpart
                       endif
                    enddo
                    !!--plot ghosts using different marker
                    do j=npart1,ntotplot(i)
                       if ((dat(ixsec,j,i).lt.xsecmax) &
                            .and.(dat(ixsec,j,i).gt.xsecmin)) then
                          call pgpt(1,xplot(j),yplot(j),imarkg)
                          if (ilabelpart) then
                             !!--plot ghost labels
                             call pgnumb(j,0,1,string,nc)
                             call pgsch(0.5*charheight)
                             call pgtext(xplot(j),yplot(j),string(1:nc))
                             call pgsch(charheight)
                          endif! ilabelpart
                       endif
                    enddo

                 else          
                    !
                    !--or simply plot all particles
                    !
                    !!--plot particle positions
                    if (iplotpart) call pgpt(npart(i),xplot(1:npart(i)),yplot(1:npart(i)),imark)
                    !!--plot ghost particles with different marker
                    if (iplotpart .and. iplotghost .and. nghost(i).gt.0) then
                       nghoststart = npart(i) + 1
                       nghostend = npart(i) + nghost(i)
                       call pgpt(nghost(i),xplot(nghoststart:nghostend), &
                            yplot(nghoststart:nghostend),imarkg)         
                    endif
                    !!--plot circles of interaction (circles of radius 2h around each particle)
                    if (plotcirc) then  
                       if (plotcircall) then
                          print*,'plotting circles of interaction',npart(i) 
                          do j=1,npart(i)
                             call pgcirc(xplot(j),yplot(j),2.*dat(ih,j,i))
                          enddo
                       else 
                          icoords = 2
                          print*,'plotting circles of interaction',ncircpart
                          do n = 1,ncircpart   
                             if (icoords.gt.1) then
                                print*,'coordinate system = ',icoords     
                                call plot_kernel_gr(icoords,xplot(icircpart(n)),  &
                                     yplot(icircpart(n)),2*dat(ih,icircpart(n),i))
                             else
                                call pgcirc(xplot(icircpart(n)),  &
                                     yplot(icircpart(n)),2*dat(ih,icircpart(n),i))
                             endif
                          enddo
                       endif
                    endif
                    if (ilabelpart) then
                       !!--plot particle labels
                       print*,'plotting particle labels ',ntotplot(i)
                       do j=1,ntotplot(i)
                          call pgnumb(j,0,1,string,nc)
                          call pgsch(0.5*charheight)
                          call pgtext(xplot(j),yplot(j),string(1:nc))
                          call pgsch(charheight)
                       enddo
                    endif! ilabelpart

                 endif! if x_sec else    
              endif! if irender
              !-----------------------------------------------------------------------------
              ! sink particles (want these to appear on both particle plots and renderings)
              !-----------------------------------------------------------------------------

              !--plot sink particles with different marker again

              nsink = ntot(i) - nghost(i) - npart(i)
              if (iplotsink .and. nsink.gt.0) then
                 nsinkstart = npart(i) + nghost(i) + 1
                 nsinkend = ntot(i)
                 print*,' plotting ',nsink,' sink particles...'
                 call pgsci(2)
                 call pgsch(2.*charheight)
                 call pgpt(nsink,xplot(nsinkstart:nsinkend), &
                      yplot(nsinkstart:nsinkend),imarksink) 
                 call pgsch(charheight)
                 call pgsci(1)
              endif

              !----------------------------
              ! vector maps
              !----------------------------    

              !!--velocity vector map
              if  (ivecplot.ne.0) then
                 if (ivecplot.eq.1 .and. ivx.ne.0) then
                    ivecx = ivx + iplotx - 1 ! 
                    ivecy = ivx + iploty - 1
                    print*,'plotting velocity field'        
                 elseif (ivecplot.eq.2 .and. iBfirst.ne.0) then
                    ivecx = iBfirst + iplotx - 1 ! 
                    ivecy = iBfirst + iploty - 1
                    print*,'plotting magnetic field'
                    vecmax = Bmax
                 endif
                 !!--check for errors
                 if ((ivecx.le.ndim).or.(ivecx.gt.ndataplots) &
                      .or.(ivecy.le.ndim).or.(ivecy.gt.ndataplots)) then
                    print*,'error finding location of vector plot in array'
                 else
                    !!--determine number of pixels in rendered image (npix = pixels in x direction)
                    pixwidth = (xmax-xmin)/real(npixvec)
                    npixyvec = int((ymax-ymin)/pixwidth) + 1
                    !--copy appropriate velocity data to a 2D array
                    !do j=1,ntotplot(i)
                    !   vecplot(1,j) = dat(ivecx,j,i)
                    !   vecplot(2,j) = dat(ivecy,j,i)
                    !enddo
                    if (iadapt) then
                       vecmax = -1.0  ! plot limits then set in vectorplot
                    else                    
                       vecmax = max(lim(ivecx,2),lim(ivecy,2))
                    endif
                    
                    !!--plot arrows in either background or foreground colour
                    if (use_backgnd_color_vecplot) call pgsci(0)
                    !!--call routine to do vector plot off particles      
                    call vectorplot(xplot(1:ntotplot(i)),yplot(1:ntotplot(i)),  &
                         xmin,ymin,pixwidth, &
                         dat(ivecx,1:ntotplot(i),i),dat(ivecy,1:ntotplot(i),i), &
                         vecmax,ntotplot(i),npixvec,npixyvec)
                    if (use_backgnd_color_vecplot) call pgsci(1)
                    !!--old stuff here is to plot arrows on the particles themselves
                    !scale = 0.08*(lim(iploty,2)-lim(iploty,1))
                    !call pgsch(0.35)! character height (size of arrow head)
                    !do j=1,ntotplot(i)
                    !   call pgarro(yplot(i),xplot(i), &
                    !   yplot(i)+vecplot(2,i)*scale,   &
                    !   xplot(i)+vecplot(1,i)*scale)
                    !enddo
                    !       call pgsch(1.0)    ! reset character height
                 endif
              endif
              !
              !--print legend if this is the first plot on the page
              !    
              if (nyplot.eq.1) call legend(time(i))    

              !
              !--%%%%%%%%%%%%% end loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
              !
           enddo over_cross_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! not both coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !-----------------------
           ! set up pgplot page
           !-----------------------
        elseif ((iploty.gt.ndim .or. iplotx.gt.ndim)  &
             .and.(iploty.le.ndataplots .and. iplotx.le.ndataplots)) then

           if ((ipagechange).or.((.not.ipagechange).and.(i.eq.nstart))) then
              ! call pgenv(limx(iplotx,1),limx(iplotx,2), &
              !            lim(iploty,1),lim(iploty,2),0,0)! 0 for no axes
              call pgpage
              if ((nacross*ndown).gt.1) then
                 if (axes) then
                    call pgsvp(0.2,0.99,0.2,0.98)
                 else! if no axes use full viewport
                    call pgsvp(0.02,0.98,0.02,0.98)
                 endif
              else
                 if (axes) then
                    call pgsvp(0.1,0.9,0.1,0.9)       
                 else
                    call pgsvp(0.02,0.98,0.02,0.98)
                 endif
              endif
              call pgswin(xmin,xmax,ymin,ymax)
              if (axes) call pgbox('bcnst'//logx,0.0,0,'1bvcnst'//logy,0.0,0)              
           elseif (nyplot.eq.1) then
              call pgpanl(1,1)
           else
              call pgpage
           endif

           !---------------------------------
           ! set plot limits and label plot
           !---------------------------------


           !--print plot limits to screen
           print 34, time(i),i
           print*,trim(labely),'min,max = ',ymin,ymax
           print*,trim(labelx),'min,max = ',xmin,xmax
           !
           !--set plot limits
           !    
           call pgswin(xmin,xmax,ymin,ymax)

           if (axes) then
              if (((nyplots-nyplot).lt.nacross).or.(.not.isamexaxis)) then
                 !--print x and y labels
                 call pgmtxt('l',3.0,0.5,1.0,labely)
                 call pglab(labelx,' ',' ' )   !!trim(title)) 
              else
                 !--print y labels only    
                 call pgmtxt('l',3.0,0.5,1.0,labely)
                 !      call pglab(' ',labely,trim(title))
              endif
           endif

           !--------------------------------
           ! now plot particles
           !--------------------------------

           if ((i.eq.nstart).and.iplotlinein) then! plot initial conditions as dotted line
              call pgsls(linestylein)
           else
              !--plot time on plot
              if (nyplot.eq.1) call legend(time(i))
              !--plot particles
              call pgsls(1)
              call pgsch(1.0)! reset character height before plotting particles
              call pgpt(npart(i),xplot(1:npart(i)),yplot(1:npart(i)),imark)
           endif
           !--plot line joining the particles
           if (iplotline.or.(iplotlinein.and.(i.eq.nstart))) then
              call pgline(npart(i),xplot(1:npart(i)),yplot(1:npart(i)))     
           endif
           !--plot ghost particles with different marker
           if (iplotghost.and.nghost(i).gt.0) then
              nghoststart = npart(i) + 1
              nghostend = npart(i) + nghost(i)
              call pgpt(nghost(i),xplot(nghoststart:nghostend), &
                   yplot(nghoststart:nghostend),imarkg)
           endif
           !--plot sink particles with different marker again
           nsink = ntot(i) - npart(i) - nghost(i)
           nsinkstart = npart(i) + nghost(i) + 1
           nsinkend = ntot(i)
           if (iplotsink .and. nsink.gt.0) then
              print*,'plotting ',nsink,' sinks...'
              call pgpt(nsink,xplot(nsinkstart:nsinkend), &
                   yplot(nsinkstart:nsinkend),imarksink) 
           endif
           !--plot circles of interaction (error bar of length 2h on co-ordinate axis)
           if (plotcirc) then
              !!--on all particles      
              if (plotcircall) then
                 if (iplotx.le.ndim) then
                    print*,'plotting error bars x axis',npart(i) 
                    call pgerrb(5,npart(i),xplot(1:npart(i)), &
                         yplot(1:npart(i)),2.*dat(ih,1:npart(i),i),1.0)
                 elseif (iploty.le.ndim) then
                    print*,'plotting error bars y axis',npart(i) 
                    call pgerrb(6,npart(i),xplot(1:npart(i)), &
                         yplot(1:npart(i)),2.*dat(ih,1:npart(i),i),1.0)
                 endif
              else 
                 !!--only on specified particles
                 do n=1,ncircpart
                    if (iplotx.le.ndim) then
                       print*,'plotting error bar x axis',icircpart(n)
                       call pgerrb(5,1,xplot(icircpart(n)),yplot(icircpart(n)), &
                            2.*dat(ih,icircpart(n),i),1.0)
                    elseif (iploty.le.ndim) then
                       print*,'plotting error bar y axis',icircpart(n)
                       call pgerrb(6,1,xplot(icircpart(n)),yplot(icircpart(n)), &
                            2.*dat(ih,icircpart(n),i),1.0)      
                    endif
                 enddo
              endif
           endif

           call pgsls(1)! reset 
           call pgsch(charheight)
           !
           !--plot average line
           !
           if (iplotav) call plot_average(xplot(1:npart(i)), &
                yplot(1:npart(i)),npart(i),nbins)
           !
           !--plot particle labels
           !
           if (ilabelpart) then
              print*,'plotting particle labels ',ntotplot(i)
              do j=1,ntotplot(i)
                 call pgnumb(j,0,1,string,nc)
                 call pgsch(0.5*charheight)
                 call pgtext(xplot(j),yplot(j),string(1:nc))
                 call pgsch(charheight)
              enddo
           endif! ilabelpart

        elseif (iploty.le.numplot) then! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! additional plots (not plots of particle data - e.g where some additional 
! information is read from a file and plotted on the same page as the 
! particle plots, or where some additional plot is calculated
! from the particle data, such as errors etc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !
           !--setup plotting page as normal
           !
           if ((ipagechange).or.((.not.ipagechange).and.(i.eq.nstart))) then
              call pgpage
              if ((nacross*ndown).gt.1) then
                 if (imulti) then
                    call pgsvp(0.2,0.99,0.2,0.98)
                 else
                    call pgsvp(0.2,0.99,0.2,0.98)
                 endif
              else
                 call pgsvp(0.1,0.9,0.1,0.9)       
              endif
              ! call pgswin(lim(iplotx,1),lim(iplotx,2), &
              !             lim(iploty,1),lim(iploty,2))
              ! call pgbox('bcnst',0.0,0,'1bvcnst',0.0,0)              
           elseif (nyplot.eq.1) then
              call pgpanl(1,1)
           else
              call pgpage
           endif
           !
           !--then call subroutine to plot the additional plot
           !
           ! e.g. call routine to do convergence plot here
           if (iexact.eq.4) then
              call exact_toystar_acplane(atstar,ctstar,sigma,gamma(i))
           endif
           !
           !--power spectrum plots (uses x and data as yet unspecified)
           !
           if (iploty.eq.ipowerspec) then 

              if (.not.idisordered) then! interpolate first
                 !!--allocate memory for 1D grid (size = 2*npart)
                 ngrid = 2*npart(i)
                 if (allocated(datpix1D)) deallocate(datpix1D)
                 if (allocated(xgrid)) deallocate(xgrid)
                 allocate (datpix1D(ngrid))
                 allocate (xgrid(ngrid))
                 !!--set up 1D grid
                 xmin = lim(ix(1),1)
                 xmax = lim(ix(1),2)
                 dxgrid = (xmax-xmin)/ngrid
                 do igrid = 1,ngrid
                    xgrid(igrid) = xmin + igrid*dxgrid - 0.5*dxgrid
                 enddo
                 !!--interpolate to 1D grid  
                 call interpolate1D(dat(ix(1),1:npart(i),i), & 
                      dat(ipmass,1:npart(i),i),dat(irho,1:npart(i),i), &
                      dat(ih,1:npart(i),i),dat(ipowerspecy,1:npart(i),i), & 
                      npart(i),xmin,datpix1D,ngrid,dxgrid)
                 !!--plot interpolated 1D data to check it
                 !print*,'plotting interpolated data...'
                 !call pgswin(xmin,xmax,minval(datpix1D),maxval(datpix1D),0,1)
                 !call pgbox('BCNST',0.0,0,'1BVCNST',0.0,0)      
                 !call pglabel('x',label(ipowerspecy),'1D interpolation')
                 !call pgline(ngrid,xgrid,datpix1D)
                 !read*
                 !call pgpage! change page

                 !!--call power spectrum calculation on the even grid
                 call plot_powerspectrum(ngrid,xgrid,datpix1D,idisordered)              
              else
                 !!--or else call power spectrum calculation on the particles themselves    
                 call plot_powerspectrum(npart(i), dat(ix(1),1:npart(i),i), &
                      dat(ipowerspecy,1:npart(i),i),idisordered)
              endif
           endif
           !
           !--if this is the first plot on the page, print legend
           !
           if (nyplot.eq.1) call legend(time(i))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! if plot not in correct range
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
           call pgpage! just skip to next plot

        endif   ! ploty = whatever

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! plot exact solution on top of the plot already on the page
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        select case(iexact)
        case(1)! polytrope
           if ((iploty.eq.irho).and.(iplotx.eq.irad)) &
                call pgline(ipolyc,rad(1:ipolyc),den(1:ipolyc))

        case(2)! soundwave
           if ((iploty.eq.irho).and.(iplotx.eq.1).and.(i.ne.1)) then
              call exact_swave(time(i),delta,lambda,gamma(i), &
                   xplot(1:npart(i)),yplot(1:npart(i)), &
                   dat(iutherm,1:npart(i),i),npart(i))
           endif

        case(3)! sedov blast wave
           if ((iploty.eq.irho).and.(iplotx.eq.irad)) &
                call exact_sedov(time(i),gamma(i),xplot(1:npart(i)), &
                yplot(1:npart(i)),npart(i))

        case(4)! toy star
           !    totmass = sum(dat(ipmass,1:npart(i),i))
           !    htstar = (0.75*totmass)**(2./3.)*ctstar**(1./3.)
           !    htstar = 1.0    
           !    print*,' totmass,h,a,c in = ',totmass,htstar,atstar,ctstar
           if (iBfirst.ne.0) then
              sigma = sigma0
           else
              sigma = 0.
           endif
           if (ndim.eq.1) then
              !
              !--1D toy star solutions
              !
              if (iplotx.eq.1 .or. iplotx.eq.irad) then! if x axis is x or r
                 if (iploty.eq.irho) then
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,1)
                 elseif (iploty.eq.ipr) then
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,2)       
                 elseif (iploty.eq.iutherm) then
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,3)       
                 elseif (iploty.eq.ivx) then
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,4)       
                 elseif (iploty.eq.ibfirst+1) then
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,5)
                 endif
              elseif (iplotx.eq.irho) then
                 if (iploty.eq.ibfirst+1) then
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,6)       
                 endif
              endif

              if (iploty.eq.iacplane) then! plot point on a-c plane
                 call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,7)
              endif
           else
              !
              !--2D and 3D toy star solutions
              !
              if ((iplotx.eq.1 .and. iploty.eq.ivx) &
                   .or. (iplotx.eq.2 .and. iploty.eq.ivx+1)) then
                 call exact_toystar2D(time(i),gamma(i), &
                      htstar,atstar,ctstar,sigma,norder,4)
              endif
              if (iplotx.eq.irad) then
                 if (iploty.eq.irho) then
                    call exact_toystar2D(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,1)
                 elseif (iploty.eq.ipr) then
                    call exact_toystar2D(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,2)
                 elseif (iploty.eq.iutherm) then
                    call exact_toystar2D(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,3)
                 elseif (iploty.eq.ivx .or. iploty.eq.ivx+1) then
                    call exact_toystar2D(time(i),gamma(i), &
                         htstar,atstar,ctstar,sigma,norder,4)
                 elseif (iploty.eq.ike) then
                    call exact_toystar2D(time(i),gamma(i), &
                         htstar,atstar,ctstar,sigma,norder,4)
                 endif
              endif
           endif

        case(5) ! mhd shock tubes
           if (iplotx.eq.1) then
              !       print*,'rootname = ',rootname,rootname(5:5)
              !--if not already set, try to determine solution to plot from filename
              if (ishk.eq.0) ishk = int_from_string(rootname(5:5))
              !--otherwise prompt for shock type       
              if (ishk.eq.0) then ! prompt
                 call prompt('enter shock solution to plot',ishk,0,6)
              endif
              if (iploty.eq.irho) then
                 call exact_mhdshock(1,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.ipr) then
                 call exact_mhdshock(2,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.ivx) then
                 call exact_mhdshock(3,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.ivx+1) then
                 call exact_mhdshock(4,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.ivlast) then
                 call exact_mhdshock(5,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.ibfirst+1) then
                 call exact_mhdshock(6,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.iblast) then
                 call exact_mhdshock(7,ishk,time(i),gamma(i),xmin,xmax)
              elseif (iploty.eq.iutherm) then
                 call exact_mhdshock(8,ishk,time(i),gamma(i),xmin,xmax)
              endif
           endif
        end select
        !
        !--plot h = (1/rho)^(1/ndim)
        !
        if ((iploty.eq.ih).and.(iplotx.eq.irho)) then
           call exact_rhoh(hfact,ndim)
        endif

     enddo over_plots ! over plots per timestep (nyplot)

  enddo over_timesteps

  if (animate) then
     print*,'press return to finish'
     read*
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

300 continue
  call pgend

  return
end subroutine main
