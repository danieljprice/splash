!
! This subroutine drives the main plotting loop
!
subroutine main(ipicky,ipickx,irender,ivecplot)
  use params
  use exact_params
  use filenames
  use labels
  use limits
  use multiplot
  use particle_data
  use prompting
  use settings 
  use transforms
  implicit none
  integer, intent(in) :: ipicky, ipickx, irender, ivecplot

  integer, parameter :: maxtitles = 50
  integer :: i,j,k,n,ierr
  integer :: iplotx,iploty,irenderplot,ivectorplot,ivecx,ivecy
  integer :: nyplot,nyplots      
  integer :: ninterp,npart1,npartdim
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: npixyvec
  integer :: irenderprev, istepprev, iadvance
  integer :: nsink,nsinkstart,nsinkend,nghoststart,nghostend
  integer :: int_from_string
  integer :: ngrid
  integer :: just, ntitles
  integer :: iplots,iplotsonpage

  character(len=8) :: string     ! used in pgplot calls
  real, parameter :: pi = 3.1415926536
  real, dimension(maxpart) :: xplot,yplot
  real, dimension(:), allocatable :: datpix1D, xgrid
  real, dimension(:,:), allocatable :: datpix,vecpixx,vecpixy
  real, dimension(:,:,:), allocatable :: datpix3D
  real, dimension(ndim) :: xcoords
  real :: xmin,xmax,ymin,ymax,zmin,zmax,ymean
  real :: vecmax,rendermin,rendermax
  real :: xsecmin,xsecmax,dxsec,xsecpos
  real :: pixwidth
  real :: charheight, charheightmm
  real :: dxgrid
  real, dimension(2) :: angles

  logical :: iplotcont,x_sec,isamexaxis,isameyaxis
  logical :: log, inewpage, tile_plots, debug

  character(len=60) :: title,titlex
  character(len=20) :: labelx,labely,labelrender
  character(len=60), dimension(maxtitles) :: titlelist

  debug = .false.

  !------------------------------------------------------------------------
  ! initialisations
  !------------------------------------------------------------------------

  x_sec = xsec_nomulti
  iplotcont = iplotcont_nomulti
  title = ' '
  titlex = ' '
  isamexaxis = .true.  ! same x axis on all plots? (only relevant for >1 plots per page)
  isameyaxis = .true.  ! same y axis on all plots?
  tile_plots = .false.
  iplots = 0 ! counter for how many plots have been plotted in total
  iplotsonpage = 0  ! counter for how many plots on page

  if (ndim.eq.1) x_sec = .false. ! can't have xsec in 1D
  nxsec = 1

  imulti = .false.
  if (ipicky.eq.numplot+1) imulti = .true.

  if (imulti) then   ! multiplot
     !
     !--for a multiplot, set current plot to first in multiplot array
     !
     imulti=.true.
     iplotx = multiplotx(1)
     iploty = multiploty(1)
     nyplots = nyplotmulti
     !
     !--if doing multiplot can only take a single cross section slice      
     !
     flythru = .false.
     nxsec = 1
     !
     !--work out whether to tile plots and make labelling decisions
     !
     if (any(multiplotx(1:nyplotmulti).ne.multiplotx(1))) then
        isamexaxis = .false.
     endif
     if (any(multiploty(1:nyplotmulti).ne.multiploty(1))) then
        isameyaxis = .false.
     endif
  else
     !
     !--or else set number of plots = 1 and use ipicky and ipickx
     !
     nyplots = 1 
     iploty = ipicky
     iplotx = ipickx
  endif
  !
  !--work out whether or not to tile plots on the page
  !
  tile_plots = tile .and. isamexaxis.and.isameyaxis.and. .not. iadapt

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
	call prompt(' enter '//trim(label(ixsec))//' position for cross-section slice:', &
	             xsecpos_nomulti,lim(ixsec,1),lim(ixsec,2))
        !!--default thickness is half of the average particle spacing
	npartdim = int(maxval(npart(nstart:n_end))**(1./real(ndim)))
	print*,'average # of particles in each dimension = ',npartdim
	dxsec = (lim(ixsec,2)-lim(ixsec,1))/float(npartdim)
	call prompt(' enter thickness of cross section slice:', &
	             dxsec,0.0,lim(ixsec,2)-lim(ixsec,1))  
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
     call read_titles(titlelist,ntitles,maxtitles)

     !!--initialise pgplot
     if (tile_plots) then
        call pgbegin(0,'?',1,1)
     else
        call pgbegin(0,'?',nacross,ndown)
     endif
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
     if ((nacross.gt.1).or.(ndown.gt.1)) title = ' '

     if (tile_plots) then
        call pgbegin(0,'?',1,1)
     else
        call pgbegin(0,'?',nacross,ndown)
     endif

  endif
  !!------------------------------------------------------------------------
  ! general initialisations

  !!--set paper size if necessary
  if (ipapersize.gt.0 .and. papersizex.gt.0.0 .and. aspectratio.gt.0.0 ) then
     call pgpaper(papersizex,aspectratio)
  endif
  !!--turn off page prompting
  if (animate .or. interactive) call pgask(.false.)
  !!if (animate .and. .not. interactive) call pgbbuf !! start buffering output

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
  !--set character height in mm
  !
  charheightmm = 4.0
  !!if ((ndown*nacross).gt.1 .and..not. tile_plots) charheight = 2.0
  !      charheight = 0.5*(nacross+ndown)

  !------------------------------------------------------------------------      
  ! loop over timesteps (flexible to allow going forwards/backwards in
  !                      interactive mode)
  !------------------------------------------------------------------------            
  i = nstart
  iadvance = nfreq   ! amount to increment timestep by (changed in interactive)

  over_timesteps: do while (i.le.n_end)
     npart1 = npart(i) + 1
     irenderprev = 0
     istepprev = 0  
     !-------------------------------------
     ! loop over plots per timestep
     !-------------------------------------
     over_plots: do nyplot=1,nyplots
        !--make sure character height is set correctly
        call danpgsch(charheightmm,2) ! set in mm
	call pgqch(charheight) ! in PGPLOT scaled units
        !--for consecutive plots (ie. if not multi but nyplots > 1 plots consecutive numbers)	     
        iploty = ipicky + nyplot - 1
        !--set current x, y plot from multiplot array
        if (imulti) then
           iploty = multiploty(nyplot)
           iplotx = multiplotx(nyplot)
        endif
        !--------------------------------------------------------------
        !  copy from main dat array into xplot, yplot 
	!  also set labels and plot limits
        !--------------------------------------------------------------
        if (iploty.le.ndataplots .and. iplotx.le.ndataplots) then
           xplot = dat(:,iplotx,i)
           yplot = dat(:,iploty,i)
           labelx = label(iplotx)
           labely = label(iploty)
           if (iadvance.ne.0) then
	      xmin = lim(iplotx,1)
              xmax = lim(iplotx,2)
              ymin = lim(iploty,1)
              ymax = lim(iploty,2)
	   endif
           
           !--change coordinate system if relevant
           if (icoordsnew.ne.icoords) then
              !--do this if one is a coord but not if rendering
              if ((iplotx.le.ndim .or. iploty.le.ndim) &
                   .and..not.(iplotx.le.ndim.and.iploty.le.ndim &
                   .and.irenderplot.gt.ndim)) then
                 print*,'changing to new coordinate system',icoords,icoordsnew
                 do j=1,ntotplot(i)
                    call coord_transform(dat(j,ix(:),i),ndim,icoords, &
                                         xcoords(:),ndim,icoordsnew)
                    if (iplotx.le.ndim) xplot(j) = xcoords(iplotx)
                    if (iploty.le.ndim) yplot(j) = xcoords(iploty)
                 enddo
              endif
           endif
           !--apply transformations (log, 1/x etc) if appropriate
           !  also change labels and limits appropriately
           if (itrans(iplotx).ne.0) then
              call transform(xplot,xplot,itrans(iplotx),maxpart)
              labelx = transform_label(labelx,itrans(iplotx))
              call transform_limits(xmin,xmax,xmin,xmax,itrans(iplotx))
           endif
           if (itrans(iploty).ne.0) then
              call transform(yplot,yplot,itrans(iploty),maxpart)
              labely = transform_label(labely,itrans(iploty))
              call transform_limits(ymin,ymax,ymin,ymax,itrans(iploty))
           endif
           
           !--write username, date on plot
           !         if (nacross.le.2.and.ndown.le.2) call pgiden

           !--adjust plot limits if adaptive plot limits set
           if (ipagechange .and. iadapt .and. iadvance.ne.0) then
              xmin = minval(xplot(1:ntotplot(i)))
              xmax = maxval(xplot(1:ntotplot(i)))*scalemax
              ymin = minval(yplot(1:ntotplot(i)))
              ymax = maxval(yplot(1:ntotplot(i)))*scalemax
	   endif

           !!-reset co-ordinate plot limits if particle tracking           
	   if (itrackpart.gt.0 .and. iadvance.ne.0) then
              if (iplotx.le.ndim) then
                 xmin = xplot(itrackpart) - xminoffset_track(iplotx)
                 xmax = xplot(itrackpart) + xmaxoffset_track(iplotx)
              endif
              if (iploty.le.ndim) then
                 ymin = yplot(itrackpart) - xminoffset_track(iploty)
                 ymax = yplot(itrackpart) + xmaxoffset_track(iploty)	   
              endif
           endif

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! plots with co-ordinates as x and y axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ((iploty.le.ndim).and.(iplotx.le.ndim)) then

           !!--set rendering options equal to settings in multiplot	 
           if (imulti) then
              irenderplot = irendermulti(nyplot)      
              ivectorplot = ivecplotmulti(nyplot)
              iplotcont = iplotcontmulti(nyplot)
              x_sec = x_secmulti(nyplot)
              xsecpos = xsecposmulti(nyplot)
           else
              irenderplot = irender
              ivectorplot = ivecplot
              iplotcont = iplotcont_nomulti
              x_sec = xsec_nomulti
              xsecpos = xsecpos_nomulti
           endif
           npixx = npix

           !!--work out coordinate that is not being plotted	 
           do j=1,ndim
              if ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j
           enddo
	   !
           !--set number of particles to use in the interpolation routines
           !  (ie. including only gas particles and ghosts)
	   ninterp = npart(i) + nghost(i)	

           !
	   !--rotate the particles about the z (and y) axes
	   !  only applies to particle plots at the moment
	   !
	   if (ndim.ge.2 .and. irotate) then
	       !
	       !--convert angles to radians
	       !
	       angles(1) = anglerot*pi/180.
	       angles(2) = angletilt*pi/180.
               print*,'rotating particles',irotate,angles
               do j=1,ntotplot(i)
                  call rotate(angles(1:ndim-1),dat(j,ix(1:ndim),i), &
		              xcoords(:),xorigin(1:ndim),ndim)
                  xplot(j) = xcoords(iplotx)
                  yplot(j) = xcoords(iploty)
               enddo
	   endif 
           
	   !------------------------------------------------------------------
           !  rendering setup and interpolation (this is the rendering done
           !  *before* the cross sections are taken, e.g. to 3D grid)
           !------------------------------------------------------------------
           if ((irenderplot.gt.ndim).and. &
                ((ndim.eq.3).or.(ndim.eq.2.and..not.x_sec))) then
              !
              !--interpolate from particles to fixed grid using sph summation
              !		
              !!--use the un-transformed co-ordinates
	      if (itrackpart.le.0 .and. iadvance.ne.0) then ! do not reset limits if particle tracking
	         xmin = lim(iplotx,1)
                 xmax = lim(iplotx,2)
                 ymin = lim(iploty,1)
                 ymax = lim(iploty,2)
	      endif
	      
              !!--determine number of pixels in rendered image (npix = pixels in x direction)
              pixwidth = (xmax-xmin)/real(npix)
              npixx = int((xmax-xmin)/pixwidth) + 1
              npixy = int((ymax-ymin)/pixwidth) + 1
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
                       allocate ( datpix(npixx,npixy) )
                       call interpolate2D( &
                            dat(1:ninterp,iplotx,i),dat(1:ninterp,iploty,i), &
                            dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i), &
                            dat(1:ninterp,ih,i),dat(1:ninterp,irenderplot,i), &
                            ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth)
                    endif
                 case(3)
                    !!--interpolation to 3D grid - then take multiple cross sections/projections
                    !!  do this if taking more than 2 cross sections, otherwise use fast xsec
                    if (x_sec.and.nxsec.gt.2) then
                       !!--allocate memory for 3D rendering array
                       allocate ( datpix3D(npixx,npixy,npixz) )
                       !!--interpolate from particles to 3D grid
                       call interpolate3D( &
                            dat(1:ninterp,iplotx,i),dat(1:ninterp,iploty,i), &
                            dat(1:ninterp,ixsec,i),dat(1:ninterp,ipmass,i),  &
                            dat(1:ninterp,irho,i),dat(1:ninterp,ih,i), &
                            dat(1:ninterp,irenderplot,i), &
                            ninterp,xmin,ymin,zmin,datpix3D,npixx,npixy,npixz,pixwidth,dxsec)
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
           if (ivectorplot.gt.0) iplotpart = iplotpartvec
           if (irenderplot.gt.0) iplotpart = .false.

           !
           !%%%%%%%%%%%%%%% loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
           !
           over_cross_sections: do k=1,nxsec

              if (x_sec) then
                 !!--for multislice cross section (flythru)
                 !!  increment the position of the current cross section slice          
		 if (flythru) xsecpos = xsecpos + dxsec
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
                            dat(1:ninterp,iplotx,i),dat(1:ninterp,iploty,i), &
                            dat(1:ninterp,ixsec,i), &
                            dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i),    &
                            dat(1:ninterp,ih,i),dat(1:ninterp,irenderplot,i), &
                            ninterp,xmin,ymin,xsecpos,datpix,npixx,npixy,pixwidth)
                    else
                       !!--do fast projection
                       call interpolate3D_projection( &
                            dat(1:ninterp,iplotx,i),dat(1:ninterp,iploty,i), &
                            dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i),   &
                            dat(1:ninterp,ih,i), dat(1:ninterp,irenderplot,i), &
                            ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth)
                    endif

                 endif ! whether 3D grid or fast renderings

                 !-------------take cross sections through 2D data------------------!
              elseif (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) then
                 !-------------------------------------------------------------------
                 !  or do a fast cross section through 2D data to 1D array
                 !-------------------------------------------------------------------
                 !!--interpolate from 2D data to 1D line
                 !!  line is specified by giving two points, (x1,y1) and (x2,y2)
                 !--set up 1D grid and allocate memory for datpix1D
                 if (iadvance.ne.0) then
		    xmin = 0.   ! distance (r) along cross section
                    xmax = SQRT((xseclineY2-xseclineY1)**2 + (xseclineX2-xseclineX1)**2)
                 endif
		 dxgrid = (xmax-xmin)/REAL(npixx)
                 call set_grid1D(xmin,dxgrid,npixx)

                 call interpolate2D_xsec( &
                      dat(1:ninterp,iplotx,i),dat(1:ninterp,iploty,i), &
                      dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i),    &
                      dat(1:ninterp,ih,i),dat(1:ninterp,irenderplot,i), &
                      ninterp,xseclineX1,xseclineY1,xseclineX2,xseclineY2, &
		      datpix1D,npixx)
                 !
                 !--find limits of datpix1D for plotting
                 !  do transformations on rendered array where appropriate
                 !  set these as ymin,ymax and set labels of plot
                 !
                 call transform(datpix1D,datpix1D,itrans(irenderplot),npixx)
                 labely = transform_label(label(irenderplot),itrans(irenderplot))
                 labelx = 'cross section'
                 !!--if adaptive limits, find limits of datpix
		 if (iadvance.ne.0) then               
		    if (iadapt) then
                       ymin = minval(datpix1D)
                       ymax = maxval(datpix1D)
                       !!--or apply transformations to fixed limits
                    else
                       call transform_limits(lim(irenderplot,1),lim(irenderplot,2), &
                            ymin,ymax,itrans(irenderplot))
                    endif
		 endif

              endif ! 2 or 3D and rendering
              !-----end of preliminary muff for 2D/3D cross sections/renderings ------------------

	      !---------------------------------
	      ! output some muff to the screen
	      !---------------------------------
	      
              print 34, time(i),i
              print*,trim(labely),'min,max = ',ymin,ymax
              print*,trim(labelx),'min,max = ',xmin,xmax
34            format (5('-'),' t = ',f8.4,', dump #',i4,1x,10('-'))
              if (x_sec.and.iplotpart) print 35,label(ixsec),xsecmin,label(ixsec),xsecmax
35            format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)
	      
              !--------------------------------------------------------------
              ! set up pgplot page (this is my version of PGENV and PGLABEL)
              !--------------------------------------------------------------

	      iplots = iplots + 1
	      iplotsonpage = iplotsonpage + 1
	      if (iplotsonpage.gt.nacross*ndown) iplotsonpage = 1
	      
	      just = 1  ! x and y axis have same scale
	      ! unless 1D xsec through 2D data
	      if (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) just = 0 
	      
	      if (tile_plots) then
	         if (iplotsonpage.eq.1 .and. ipagechange) call pgpage
                 call danpgtile(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                                labelx,labely,titlex,just,iaxis)
	      else
	         !--change the page if pagechange set
	         !  or, if turned off, between plots on first page only
	         inewpage = ipagechange .or. (iplots.le.nacross*ndown)
	         call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
	           labelx,labely,titlex,just,iaxis, &
		   isamexaxis,isameyaxis,inewpage)
	      endif
	      
	      if (iaxis.lt.0) then
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
              if (irenderplot.gt.ndim .and.    &
                   ((ndim.eq.3).or.(ndim.eq.2.and. .not.x_sec)) ) then
                 !---------------------------------------------------------------
                 ! scalar quantity which has been rendered to a 2D pixel array (datpix)
                 !---------------------------------------------------------------
                 !!--do transformations on rendered array  
                 call transform2(datpix,itrans(irenderplot),npixx,npixy)
                 labelrender = label(irenderplot)
                 !!--set label for column density (projection) plots (2268 or 2412 for integral sign)
                 if (ndim.eq.3 .and..not. x_sec) then         
                    labelrender = '\(2268) '//trim(labelrender)//' d'//trim(label(ix(ixsec)))
                 endif
                 !!--apply transformations to the label for the rendered quantity 
		 !!  but don't do this for log as we use a logarithmic axis instead
                 log = .false.
                 if (itrans(irenderplot).eq.1) then
		    log = .true.
		 else 
                    labelrender = transform_label(labelrender,itrans(irenderplot))
		 endif
		 !!--limits for rendered quantity
                 if (iadapt) then
                    !!--if adaptive limits, find limits of rendered array
                    rendermin = minval(datpix)
                    rendermax = maxval(datpix)
                 else
                    !!--or apply transformations to fixed limits
                    call transform_limits(lim(irenderplot,1),lim(irenderplot,2), &
                         rendermin,rendermax,itrans(irenderplot))
                 endif
                 !!--print plot limits to screen
                 print*,trim(labelrender),' min, max = ',rendermin,rendermax       
                 !!--call subroutine to actually render the image       
                 call render(datpix,rendermin,rendermax,trim(labelrender),  &
                      npixx,npixy,xmin,ymin,pixwidth,    &
                      icolours,iplotcont,iPlotColourBar,ncontours,log)

              elseif (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) then
                 !---------------------------------------------------------------
                 ! plot 1D cross section through 2D data (contents of datpix) 
                 !---------------------------------------------------------------  
                 call pgline(npixx,xgrid,datpix1D) 
              else
                 !-----------------------
                 ! particle plots
                 !-----------------------
                 !
                 !--if particle cross section, plot particles only in a defined coordinate range
                 !
                 if (ndim.gt.2 .and. x_sec.and.iplotpart) then
                    do j=1,npart(i)
                       if ((dat(j,ixsec,i).lt.xsecmax) &
                            .and.(dat(j,ixsec,i).gt.xsecmin)) then
                          call pgpt(1,xplot(j),yplot(j),imark)
                          !!--plot circles of interaction
                          if (plotcirc) then
                             call pgcirc(xplot(j),yplot(j),2.*dat(j,ih,i))
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
                       if ((dat(j,ixsec,i).lt.xsecmax) &
                            .and.(dat(j,ixsec,i).gt.xsecmin)) then
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
                    if (iplotpart) then
		       if (debug) print*,'plotting particles'
                       call pgpt(npart(i),xplot(1:npart(i)),yplot(1:npart(i)),imark)
                    endif
                    !!--plot ghost particles with different marker
                    if (iplotpart .and. iplotghost .and. nghost(i).gt.0) then
		       if (debug) print*,'plotting ',nghost(i),'ghosts'
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
                             call pgcirc(xplot(j),yplot(j),2.*dat(j,ih,i))
                          enddo
                       else
                          print*,'plotting circles of interaction',ncircpart
                          do n = 1,ncircpart   
                             if (icoords.gt.1) then   
                                call plot_kernel_gr(icoords,xplot(icircpart(n)),  &
                                     yplot(icircpart(n)),2*dat(icircpart(n),ih,i))
                             else
                                call pgcirc(xplot(icircpart(n)),  &
                                     yplot(icircpart(n)),2*dat(icircpart(n),ih,i))
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
                    !
		    !--enter interactive mode
		    !
                    !if (iplotpart.and.interactive) then
		    !   call interactive_part(ntotplot(i),xplot(1:ntotplot(i)), &
		    !                         yplot(1:ntotplot(i)),iadvance) 
                    !endif

                 endif! if x_sec else    
              endif! if irender
              !-----------------------------------------------------------------------------
              ! sink particles (want these to appear on both particle plots and renderings)
              !-----------------------------------------------------------------------------

              !--plot sink particles with different marker again

              nsink = ntot(i) - nghost(i) - npart(i)
              if (debug) print*,'nsink = ',nsink
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

              !--------------------------------------------------------------
              ! vector maps (can be on top of particle plots and renderings)
              !--------------------------------------------------------------

              if (ivectorplot.ne.0 .and. ndim.ge.2) then
	         !!--choose quantity to be plotted
                 if (ivectorplot.eq.1 .and. ivx.ne.0) then
                    ivecx = ivx + iplotx - 1 ! 
                    ivecy = ivx + iploty - 1
                    print*,'plotting velocity field'        
                 elseif (ivectorplot.eq.2 .and. iBfirst.ne.0) then
                    ivecx = iBfirst + iplotx - 1 ! 
                    ivecy = iBfirst + iploty - 1
                    print*,'plotting magnetic field'
                 endif
                 !!--check for errors
                 if ((ivecx.le.ndim).or.(ivecx.gt.ndataplots) &
                      .or.(ivecy.le.ndim).or.(ivecy.gt.ndataplots)) then
                    print*,'error finding location of vector plot in array'
                 else
                    !!--determine number of pixels in rendered image (npix = pixels in x direction)
                    pixwidth = (xmax-xmin)/real(npixvec)
                    npixyvec = int((ymax-ymin)/pixwidth) + 1

                    if (iadapt) then
                       vecmax = -1.0  ! plot limits then set in vectorplot
                    else                    
                       vecmax = max(lim(ivecx,2),lim(ivecy,2))
                    endif
                    
                    !!--plot arrows in either background or foreground colour
                    if (UseBackgndColorVecplot) call pgsci(0)
                    
		    if (allocated(vecpixx)) deallocate(vecpixx)
		    if (allocated(vecpixy)) deallocate(vecpixy)
		    allocate(vecpixx(npixvec,npixyvec),stat=ierr)
		    allocate(vecpixy(npixvec,npixyvec),stat=ierr)
		    if (ierr.ne.0) then
		       print*,'error allocating memory for vector cross section'
		    else
		       if (x_sec .and. ndim.eq.3) then ! take vector plot in cross section
		          !
			  !--interpolate vector from particles to cross section
			  !
			  call interpolate3D_xsec_vec(dat(1:ninterp,iplotx,i),dat(1:ninterp,iploty,i), &
                            dat(1:ninterp,ixsec,i),dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i),  &
                            dat(1:ninterp,ih,i),dat(1:ninterp,ivecx,i),dat(1:ninterp,ivecy,i), &
			    ninterp,xmin,ymin,xsecpos, &
			    vecpixx,vecpixy,npixvec,npixyvec,pixwidth)
		       else
		          !
			  !--or interpolate (via averaging) to coarser grid
			  !
                          call interpolate_vec(xplot(1:ntotplot(i)),yplot(1:ntotplot(i)), &
                            dat(1:ntotplot(i),ivecx,i),dat(1:ntotplot(i),ivecy,i), &
			    xmin,ymin,pixwidth,vecpixx,vecpixy, &
                            ntotplot(i),npixvec,npixyvec)		       
		       endif
		       !
		       !--plot rendered vector map
		       !
		       log = .false.
		       call render_vec(vecpixx,vecpixy,vecmax, &
			   npixvec,npixyvec,xmin,ymin,pixwidth,log)
		       deallocate(vecpixx,vecpixy)
                    endif
                    if (UseBackgndColorVecplot) call pgsci(1)
		 endif
              endif
              !
              !--print legend if this is the first plot on the page
              !    
              if (nyplot.eq.1 .and. i.le.nacross) then
	         call legend(time(i),hposlegend,vposlegend)
	      endif
	      if (i.le.ntitles) then
	         if (titlelist(i)(1:1).ne.' ') then
		    call pgmtxt('T',vpostitle,hpostitle,fjusttitle,trim(titlelist(i)))
	         endif
	      endif

              !
	      !--enter interactive mode
	      !	     
             if (interactive) then
	        iadvance = nfreq
		call interactive_part(ntotplot(i),iplotx,iploty,irenderplot, &
		     xplot(1:ntotplot(i)),yplot(1:ntotplot(i)), &
		     xmin,xmax,ymin,ymax,iadvance) 
             endif


              !
              !--%%%%%%%%%%%%% end loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
              !
           enddo over_cross_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! not both coordinates - these are just particle plots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        elseif ((iploty.gt.ndim .or. iplotx.gt.ndim)  &
             .and.(iploty.le.ndataplots .and. iplotx.le.ndataplots)) then

	   !---------------------------------
	   ! output some muff to the screen
	   !---------------------------------
	   
           print 34, time(i),i
           print*,trim(labely),'min,max = ',ymin,ymax
           print*,trim(labelx),'min,max = ',xmin,xmax
	      
           !--------------------------------------------------------------
           ! set up pgplot page (this is my version of PGENV and PGLABEL)
           !--------------------------------------------------------------

	   iplots = iplots + 1
	   iplotsonpage = iplotsonpage + 1
	   if (iplotsonpage.gt.nacross*ndown) iplotsonpage = 1

           just = 0  ! x and y axis have different scales
	   
	   if (tile_plots) then
	      if (iplotsonpage.eq.1 .and. ipagechange) call pgpage
              call danpgtile(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                             labelx,labely,title,just,iaxis)
	   else
 	      !--change the page if pagechange set
	      !  or, if turned off, between plots on first page only
	      inewpage = ipagechange .or. (iplots.le.nacross*ndown)
	      call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
	        labelx,labely,title,just,iaxis, &
		isamexaxis,isameyaxis,inewpage)
	   endif

           !--------------------------------
           ! now plot particles
           !--------------------------------

           if ((i.eq.nstart).and.iplotlinein) then! plot initial conditions as dotted line
              call pgsls(linestylein)
           else
              !--plot time on plot
              if (nyplot.eq.1) call legend(time(i),hposlegend,vposlegend)
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
                         yplot(1:npart(i)),2.*dat(1:npart(i),ih,i),1.0)
                 elseif (iploty.le.ndim) then
                    print*,'plotting error bars y axis',npart(i) 
                    call pgerrb(6,npart(i),xplot(1:npart(i)), &
                         yplot(1:npart(i)),2.*dat(1:npart(i),ih,i),1.0)
                 endif
              else 
                 !!--only on specified particles
                 do n=1,ncircpart
                    if (iplotx.le.ndim) then
                       print*,'plotting error bar x axis',icircpart(n)
                       call pgerrb(5,1,xplot(icircpart(n)), &
		            yplot(icircpart(n)), &
			    (2.*dat(icircpart(n):icircpart(n),ih,i)),1.0)
                    elseif (iploty.le.ndim) then
                       print*,'plotting error bar y axis',icircpart(n)
                       call pgerrb(6,1,xplot(icircpart(n)),yplot(icircpart(n)), &
                            (2.*dat(icircpart(n):icircpart(n),ih,i)),1.0)      
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
           !
           !--enter interactive mode
           !
	   iadvance = nfreq
           if (interactive) then
	      call interactive_part(npart(i),iplotx,iploty,0,xplot(1:npart(i)), &
	                            yplot(1:npart(i)),xmin,xmax,ymin,ymax,iadvance)
           endif
	   
        elseif (iploty.le.numplot) then! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! additional plots (not plots of particle data - e.g where some additional 
! information is read from a file and plotted on the same page as the 
! particle plots, or where some additional plot is calculated
! from the particle data, such as errors etc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   
           !--------------------------------------------------------------
           ! set up pgplot page, but do not plot axes or
	   ! labels
           !--------------------------------------------------------------

	   iplots = iplots + 1
	   iplotsonpage = iplotsonpage + 1
	   if (iplotsonpage.gt.nacross*ndown) iplotsonpage = 1

           just = 0  ! this is irrelevant since axes are not plotted
	   xmin = 0.0 ! these are also irrelevant
	   xmax = 1.0
	   ymin = 0.0
	   ymax = 1.0
	   
	   if (tile_plots) then
	      if (iplotsonpage.eq.1 .and. ipagechange) call pgpage
              call danpgtile(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                             ' ',' ',' ',just,-2)
	   else
 	      !--change the page if pagechange set
	      !  or, if turned off, between plots on first page only
	      inewpage = ipagechange .or. (iplots.le.nacross*ndown)
	      call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
	        ' ',' ',' ',just,-2,isamexaxis,isameyaxis,inewpage)
	   endif	   

           !--------------------------------------------------------------
           !  then call subroutine to plot the additional plot
           ! e.g. call routine to do convergence plot here
           !--------------------------------------------------------------
           !
	   !--A vs C for exact toystar solution
	   !
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
                 !!--set up 1D grid
                 xmin = lim(ix(1),1)
                 xmax = lim(ix(1),2)
                 dxgrid = (xmax-xmin)/ngrid
                 call set_grid1D(xmin,dxgrid,ngrid)

                 !!--interpolate to 1D grid  
                 call interpolate1D(dat(1:npart(i),ix(1),i), & 
                      dat(1:npart(i),ipmass,i),dat(1:npart(i),irho,i), &
                      dat(1:npart(i),ih,i),dat(1:npart(i),ipowerspecy,i), & 
                      npart(i),xmin,datpix1D,ngrid,dxgrid)
                 !!--plot interpolated 1D data to check it
                 print*,'plotting interpolated data...'
                 call pgswin(xmin,xmax,minval(datpix1D),maxval(datpix1D),0,1)
                 call pgbox('BCNST',0.0,0,'1BVCNST',0.0,0)      
                 call pglabel('x',label(ipowerspecy),'1D interpolation')
                 call pgline(ngrid,xgrid,datpix1D)
                 read*
                 call pgpage! change page

                 !!--call power spectrum calculation on the even grid
                 call plot_powerspectrum(ngrid,nfreqspec,wavelengthmax, &
		      xgrid,datpix1D,idisordered,itrans(iploty))              
              else
                 !!--or else call power spectrum calculation on the particles themselves    
                 call plot_powerspectrum(npart(i),nfreqspec,wavelengthmax, &
		      dat(1:npart(i),ix(1),i), &
                      dat(1:npart(i),ipowerspecy,i),idisordered,itrans(iploty))
              endif
           endif
           !
           !--if this is the first plot on the page, print legend
           !
           if (iplotsonpage.eq.1) call legend(time(i),hposlegend,vposlegend)

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
        case(1)! shock tube
           if (iplotx.eq.ix(1)) then
              if (iploty.eq.irho) then
                 call exact_shock(1,time(i),gamma(i),rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
              elseif (iploty.eq.ipr) then
                 call exact_shock(2,time(i),gamma(i),rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
              elseif (iploty.eq.ivx) then
                 call exact_shock(3,time(i),gamma(i),rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
              elseif (iploty.eq.iutherm) then
                 call exact_shock(4,time(i),gamma(i),rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
              endif
           endif

        case(2)! sedov blast wave
           if (iplotx.eq.irad) then
              if (iploty.eq.irho) then
                 call exact_sedov(time(i),gamma(i),rhosedov,esedov,lim(irad,2),1)
              elseif (iploty.eq.ipr) then
                 call exact_sedov(time(i),gamma(i),rhosedov,esedov,lim(irad,2),2)                 
              elseif (iploty.eq.iutherm) then
                 call exact_sedov(time(i),gamma(i),rhosedov,esedov,lim(irad,2),3)                
              elseif (iploty.eq.ike) then
                 call exact_sedov(time(i),gamma(i),rhosedov,esedov,lim(irad,2),4)                 
              endif
           endif

        case(3)! polytrope
           if (iploty.eq.irho .and. iplotx.eq.irad) call exact_polytrope(gamma(i))

        case(4)! toy star
           if (iBfirst.ne.0) then
              sigma = sigma0
           else
              sigma = 0.
           endif
           if (ndim.eq.1) then
              !
              !--1D toy star solutions
              !
              if (iplotx.eq.ix(1) .or. iplotx.eq.irad) then! if x axis is x or r
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
              if ((iplotx.eq.ix(1) .and. iploty.eq.ivx) &
                   .or. (iplotx.eq.ix(2) .and. iploty.eq.ivx+1)) then
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

        case(5)! linear wave
           if ((iploty.eq.iwaveploty).and.(iplotx.eq.iwaveplotx)) then 
	      ymean = SUM(yplot(1:npart(i)))/REAL(npart(i)) 
              call exact_wave(time(i),ampl,period,lambda,xmin,xmax,ymean)
           endif

        case(6) ! mhd shock tubes
           if (iplotx.eq.ix(1)) then
              !--prompt for shock type if not set  
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
	
	case(7) ! exact solution read from file
	   if (iplotx.eq.iexactplotx .and. iploty.eq.iexactploty) then   
	      call pgline(iexactpts,xexact,yexact)
	   endif
        end select
        !
        !--plot h = (1/rho)^(1/ndim)
        !
        if ((iploty.eq.ih).and.(iplotx.eq.irho)) then
           call exact_rhoh(hfact,ndim,dat(1:npart(i),ipmass,i),npart(i),xmin,xmax)
        endif

     enddo over_plots ! over plots per timestep (nyplot)
!
!--increment timestep
!
     i = i + iadvance
     if (i.lt.1) then
        print*,'reached first step: can''t go back'
	i = 1
     elseif (i.lt.nstart) then
        print*,'warning: i < nstart'
     endif

  enddo over_timesteps

  if (.not.interactive) then
     !!if (animate) call pgebuf
     print*,'press return to finish'
     read*
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

300 continue
  call pgend

  return

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------
! sets up a one dimensional grid of pixels
! and allocates memory for datpix1D
!--------------------------------------------
  subroutine set_grid1D(xmin1D,dxgrid1D,ngridpts)
    implicit none
    integer, intent(in) :: ngridpts
    real, intent(in) :: xmin1D, dxgrid1D
    integer :: igrid

    if (allocated(datpix1D)) deallocate(datpix1D)
    if (allocated(xgrid)) deallocate(xgrid)
    allocate (datpix1D(ngridpts))
    allocate (xgrid(ngridpts))

    do igrid = 1,ngridpts
       xgrid(igrid) = xmin1D + igrid*dxgrid1D - 0.5*dxgrid1D
    enddo

  end subroutine set_grid1D

end subroutine main
