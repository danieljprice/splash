!
! This subroutine drives the main plotting loop
!
subroutine mainloop(ipicky,ipickx,irender,ivecplot)
  use params
  use colours, only:colour_set
  use labels
  use limits, only:lim
  use multiplot
  use particle_data, only:npartoftype,time
  use prompting
  use settings_data, only:ndim,nstart,n_end,nfreq,numplot,buffer_data,imulti
  use settings_page
  use settings_render, only:icolours,iplotcont_nomulti
  use settings_vecplot, only:iplotpartvec
  use settings_xsecrot, only:ixsec,xsec_nomulti,xsecpos_nomulti,flythru,nxsec
  implicit none
  integer, intent(in) :: ipicky, ipickx, irender, ivecplot

  integer, parameter :: maxtitles = 50
  integer :: i,j,ierr,ifile
  integer :: iplotx,iploty,irenderplot,ivectorplot,ivecx,ivecy
  integer :: nyplots,npartdim      
  integer :: irenderprev, istepprev, iadvance
  integer :: ngrid
  integer :: just, ntitles
  integer :: iplots,iplotsonpage
  integer :: index1,index2,itype

  real, dimension(:), allocatable :: datpix1D, xgrid
  real :: xmin,xmax,ymin,ymax,zmin,zmax,ymean
  real :: vecmax,rendermin,rendermax,dummymin,dummymax
  real :: dxsec,xsecpos
  real :: pixwidth
  real :: charheight, charheightmm
  real :: dxgrid,xmingrid,xmaxgrid
  real :: angletempx, angletempy, angletempz

  logical :: iplotpart,iplotcont,x_sec,isamexaxis,isameyaxis
  logical :: log, inewpage, tile_plots, isave, lastplot
  logical :: initialise_xsec

  character(len=60), dimension(maxtitles) :: titlelist

  !------------------------------------------------------------------------
  ! initialisations
  !------------------------------------------------------------------------

  titlelist = ' '
  isamexaxis = .true.  ! same x axis on all plots? (only relevant for >1 plots per page)
  isameyaxis = .true.  ! same y axis on all plots?
  tile_plots = .false.
  iplots = 0 ! counter for how many plots have been plotted in total
  iplotsonpage = 0  ! counter for how many plots on page
  irenderplot = 0
  ivectorplot = 0
  x_sec = xsec_nomulti
  iplotcont = iplotcont_nomulti
  lastplot = .false.
  iplotpart = .true.
  if (ivecplot.ne.0) iplotpart = iplotpartvec
  xmin = 0.
  xmax = 0.
  ymin = 0.
  ymax = 0.

  if (ndim.eq.1) x_sec = .false. ! can't have xsec in 1D
  nxsec = 1

  if (ipicky.eq.numplot+1) then   ! multiplot
     imulti = .true.
     nyplots = nyplotmulti
     !
     !--if doing multiplot can only take a single cross section slice      
     !
     flythru = .false.
     nxsec = 1
     !
     !--work out whether to tile plots and make labelling decisions
     !
     if (any(multiplotx(1:nyplotmulti).ne.multiplotx(1))) isamexaxis = .false.
     if (any(multiploty(1:nyplotmulti).ne.multiploty(1))) isameyaxis = .false.
  else
     !
     !--or else set number of plots = 1 and use ipicky and ipickx
     !
     imulti = .false.
     nyplots = 1 
     iploty = ipicky
     iplotx = ipickx
  endif
  !
  !--work out whether or not to tile plots on the page
  !
  tile_plots = tile .and. isamexaxis.and.isameyaxis.and. .not. iadapt

  !------------------------------------------------------------------------
  ! initialise options to be set before plotting

  initialise_xsec = iploty.le.ndim .and. iplotx.le.ndim
  if (imulti) then
     do i=1,nyplotmulti
        if (multiplotx(i).le.ndim .and. multiploty(i).le.ndim) then
           initialise_xsec = .true.
        endif
     enddo
  endif
  
  if (initialise_xsec) then
     !!--work out coordinate that is not being plotted 
     ixsec = 0
     if (x_sec) then
        do j=1,ndim
           if ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j
        enddo
     endif
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
        npartdim = int(maxval(npartoftype(1,nstart:n_end))**(1./real(ndim)))
        print*,'average # of particles in each dimension = ',npartdim
        dxsec = (lim(ixsec,2)-lim(ixsec,1))/float(npartdim)
        call prompt(' enter thickness of cross section slice:', &
                     dxsec,0.0,lim(ixsec,2)-lim(ixsec,1))  
     endif
  endif

  !!--prompt for options if plotting power spectrum      
  if (iploty.eq.ipowerspec &
     .or. (imulti.and.any(multiploty(1:nyplotmulti).eq.ipowerspec))) then
     call options_powerspec
  endif

  !!--set plot titles
  ntitles = 0
  call read_titles(titlelist,ntitles,maxtitles)

  !!------------------------------------------------------------------------
  ! initialise PGPLOT

  !!--start PGPLOT
  if (tile_plots) then
     call pgbegin(0,'?',1,1)
  else
     call pgbegin(0,'?',nacross,ndown)
  endif

  !!--set paper size if necessary
  if (ipapersize.gt.0 .and. papersizex.gt.0.0 .and. aspectratio.gt.0.0 ) then
     call pgpaper(papersizex,aspectratio)
  endif
  !!--turn off page prompting
  call pgask(.false.)
  !!if (.not. interactive) call pgbbuf !! start buffering output
  
  !!--set colour table
  if (((irender.gt.ndim).or.any(irendermulti(1:nyplots).gt.ndim)) &
       .and.(icolours.gt.0)) then
     call colour_set(icolours)
  endif
  
  !!--set character height in mm
  charheightmm = 4.0
  !!if ((ndown*nacross).gt.1 .and..not. tile_plots) charheight = 2.0
  !      charheight = 0.5*(nacross+ndown)


  !------------------------------------------------------------------------      
  ! loop over timesteps (flexible to allow going forwards/backwards in
  !                      interactive mode)
  !------------------------------------------------------------------------            
  i = nstart
  iadvance = nfreq   ! amount to increment timestep by (changed in interactive)
  ifile = 1

  over_timesteps: do while (i.le.n_end)

     if (.not.buffer_data) then    
        !
        !--make sure we have data for this timestep
        !
        call get_nextstep(i,ifile)
        if (i.eq.-666) exit over_timesteps
     endif
     !
     !--check timestepping
     !
     if (i.lt.1) then
        print*,'reached first step: can''t go back'
        i = 1
     endif
     if (i.lt.nstart) then
        print*,'warning: i < nstart'
     endif        

     print 33, time(i),i
33   format (5('-'),' t = ',f9.4,', dump #',i5,1x,18('-'))
     irenderprev = 0
     istepprev = 0  

     call plotstep
!
!--increment timestep
!
     if (iadvance.eq.-666) exit over_timesteps
     i = i + iadvance

  enddo over_timesteps

  if (.not.interactive) then
     !!call pgebuf
     print*,'press return to finish'
     read*
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pgend

  return

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plotstep
  use exact, only:exact_solution, &
             atstar,ctstar,htstar,sigma,iwaveplotx,iwaveploty
  use particle_data
  use rotation
  use settings_data
  use settings_limits
  use settings_part
  use settings_page
  use settings_render
  use settings_vecplot
  use settings_xsecrot
  use settings_powerspec

  use transforms
  use interactive_routines
  use geometry
  use legends, only:legend
  use fieldlines
  use particleplots
  use powerspectrums

  implicit none
  integer :: j,k
  integer :: nyplot
  integer :: ninterp
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: npixyvec

  real, parameter :: pi = 3.1415926536
  real, parameter :: tol = 1.e-10 ! used to compare real numbers
  real, dimension(:,:), allocatable :: datpix,vecpixx,vecpixy
  real, dimension(:,:,:), allocatable :: datpix3D
  real, dimension(ndim) :: xcoords
  real, dimension(max(maxpart,2000)) :: xplot,yplot,zplot
  real :: angleradx, anglerady, angleradz
  real :: xsecmin,xsecmax

  character(len=len(label(1))+20) :: labelx,labely,labelz,labelrender,labelvecplot

34   format (25(' -'))

  !-------------------------------------
  ! loop over plots per timestep
  !-------------------------------------
  over_plots: do nyplot=1,nyplots

     if (nyplot.gt.1) print 34 
     !--make sure character height is set correctly
     call danpgsch(charheightmm,2) ! set in mm
     call pgqch(charheight) ! in PGPLOT scaled units

     !--set current x, y, render and vector plot from multiplot array
     if (imulti) then
        iploty = multiploty(nyplot)
        iplotx = multiplotx(nyplot)
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
     !--------------------------------------------------------------
     !  copy from main dat array into xplot, yplot 
     !  also set labels and plot limits
     !--------------------------------------------------------------
     if (iploty.le.ndataplots .and. iplotx.le.ndataplots) then
        xplot = dat(:,iplotx,i)
        yplot = dat(:,iploty,i)
        zplot = 0. !--reset later if x-sec
        xsecmin = 0. !-- " " 
        xsecmax = 0.
        labelx = label(iplotx)
        labely = label(iploty)
        if (iadvance.ne.0) then
           xmin = lim(iplotx,1)
           xmax = lim(iplotx,2)
           ymin = lim(iploty,1)
           ymax = lim(iploty,2)
           angletempx = anglex
           angletempy = angley
           angletempz = anglez
        endif

        !--change coordinate system if relevant
        if (icoordsnew.ne.icoords) then
           !--do this if one is a coord but not if rendering
           if ((iplotx.le.ndim .or. iploty.le.ndim) &
                .and..not.(iplotx.le.ndim.and.iploty.le.ndim &
                .and.irenderplot.gt.ndim)) then
              print*,'changing to new coordinate system',icoords,icoordsnew
              do j=1,ntot(i)
                 call coord_transform(dat(j,ix(1:ndim),i),ndim,icoords, &
                                      xcoords(1:ndim),ndim,icoordsnew)
                 if (iplotx.le.ndim) xplot(j) = xcoords(iplotx)
                 if (iploty.le.ndim) yplot(j) = xcoords(iploty)
              enddo
           endif
        endif
        !--apply transformations (log, 1/x etc) if appropriate
        !  also change labels and limits appropriately
        if (.not.(iplotx.le.ndim .and. iploty.le.ndim)) then
           if (itrans(iplotx).ne.0) then
              call transform(xplot,itrans(iplotx))
              labelx = transform_label(labelx,itrans(iplotx))
              call transform_limits(xmin,xmax,itrans(iplotx))
           endif
           if (itrans(iploty).ne.0) then
              call transform(yplot,itrans(iploty))
              labely = transform_label(labely,itrans(iploty))
              call transform_limits(ymin,ymax,itrans(iploty))
           endif
        endif

        !--write username, date on plot
        !         if (nacross.le.2.and.ndown.le.2) call pgiden

        !--adjust plot limits if adaptive plot limits set
        if (ipagechange .and. iadapt .and. iadvance.ne.0) then
           xmin = 1.e12
           xmax = -1.e12
           ymin = 1.e12
           ymax = -1.e12
           !--find maximum over all particle types being plotted
           index1 = 1
           do itype=1,maxparttypes
              index2 = index1 + npartoftype(itype,i) - 1
              if (iplotpartoftype(itype).and.npartoftype(itype,i).gt.0) then
                 xmin = min(xmin,minval(xplot(index1:index2)))
                 xmax = max(xmax,maxval(xplot(index1:index2))*scalemax)
                 ymin = min(ymin,minval(yplot(index1:index2)))
                 ymax = max(ymax,maxval(yplot(index1:index2))*scalemax)
              endif
              index1 = index2 + 1
           enddo
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

        npixx = npix

        !!--work out coordinate that is not being plotted         
        if (x_sec) then
           do j=1,ndim
              if ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j
           enddo
        endif
        if (ixsec.ne.0) then
           zplot(:) = dat(:,ixsec,i)
           labelz = label(ixsec)
        endif
        !
        !--set number of particles to use in the interpolation routines
        !  (ie. including only gas particles and ghosts)
        !--if plotting ghost particles, set ntotplot = ntot, else ntot=npart
        !
        ninterp = npartoftype(1,i)
        if (labeltype(2)(1:5).eq.'ghost' .and. iplotpartoftype(2)) then
           ninterp = ninterp + npartoftype(2,i)
        endif
        !
        !--rotate the particles about the z (and y) axes
        !  only applies to particle plots at the moment
        !
        if (ndim.ge.2 .and. irotate) then
            !
            !--convert angles to radians
            !
            angleradz = angletempz*pi/180.
            print*,'rotating particles about z by ',angletempz
            if (ndim.eq.3) then
               anglerady = angletempy*pi/180.
               angleradx = angletempx*pi/180.
               print*,'rotating particles about y by ',angletempy
               print*,'rotating particles about x by ',angletempx
            endif
            do j=1,ntot(i)
               xcoords(1:ndim) = dat(j,ix(1:ndim),i) - xorigin(1:ndim)
               if (ndim.eq.2) then
                  call rotate2D(xcoords(:),angleradz)
               elseif (ndim.eq.3) then
                  call rotate3D(xcoords(:),angleradx,anglerady,angleradz)
               endif
               xplot(j) = xcoords(iplotx) + xorigin(iplotx)
               yplot(j) = xcoords(iploty) + xorigin(iploty)
               if (ixsec.gt.0) then
                  zplot(j) = xcoords(ixsec) + xorigin(ixsec)
               endif
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
              !                   npixz = int((zmax-zmin)/pixwidth) + 1                 
              !!--number of z pixels is equal to number of cross sections
              npixz = nxsec
              print*,'npixz = ',npixz
           endif

           !!--if rendering array is the same as the previous plot, reuse the array                
           if (irenderplot.eq.irenderprev .and. i.eq.istepprev) then
              if (.not.x_sec .or. (x_sec.and.ndim.eq.3.and.nxsec.gt.2)) then
                 print*,'same rendering, using previous array...'
              endif
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
                         xplot(1:ninterp),yplot(1:ninterp), &
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
                         xplot(1:ninterp),yplot(1:ninterp), &
                         zplot(1:ninterp),dat(1:ninterp,ipmass,i),  &
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

           if (k.gt.1) print 34 

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
                         xplot(1:ninterp),yplot(1:ninterp), &
                         zplot(1:ninterp), &
                         dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i),    &
                         dat(1:ninterp,ih,i),dat(1:ninterp,irenderplot,i), &
                         ninterp,xmin,ymin,xsecpos,datpix,npixx,npixy,pixwidth)
                 else
                    !!--do fast projection
                    call interpolate3D_projection( &
                         xplot(1:ninterp),yplot(1:ninterp), &
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
              call transform(datpix1D,itrans(irenderplot))
              labely = transform_label(label(irenderplot),itrans(irenderplot))
              labelx = 'cross section'
              !!--if adaptive limits, find limits of datpix
              if (iadvance.ne.0) then               
                 if (iadapt) then
                    ymin = minval(datpix1D)
                    ymax = maxval(datpix1D)
                 else
                    !!--or apply transformations to fixed limits
                    ymin = lim(irenderplot,1)
                    ymax = lim(irenderplot,2)
                    call transform_limits(ymin,ymax,itrans(irenderplot))
                 endif
              endif

           endif ! 2 or 3D and rendering
           !-----end of preliminary muff for 2D/3D cross sections/renderings ------------------

           !---------------------------------
           ! output some muff to the screen
           !---------------------------------

           print*,trim(labely),'min,max = ',ymin,ymax
           print*,trim(labelx),'min,max = ',xmin,xmax
           if (x_sec.and.iplotpart) print 35,label(ixsec),xsecmin,label(ixsec),xsecmax
35            format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

           !--------------------------------------------------------------
           ! set up pgplot page (this is my version of PGENV and PGLABEL)
           !--------------------------------------------------------------

           iplots = iplots + 1
           iplotsonpage = iplotsonpage + 1
           if (iplotsonpage.gt.nacross*ndown) iplotsonpage = 1

           just = 1  ! x and y axis have same scale
           ! unless 1D xsec through 2D data or non-cartesian
           if ((irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) &
               .or.(icoordsnew.ne.1)) then
              just = 0 
           endif

           if (tile_plots) then
              if (iplotsonpage.eq.1 .and. ipagechange) call pgpage
              call danpgtile(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                             trim(labelx),trim(labely),' ',just,iaxis)
           else
              !--change the page if pagechange set
              !  or, if turned off, between plots on first page only
              inewpage = ipagechange .or. (iplots.le.nacross*ndown)
              call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                trim(labelx),trim(labely),' ',just,iaxis, &
                isamexaxis,isameyaxis,inewpage)
           endif

           if (irotate .and. irotateaxes.gt.0) then
              if (ndim.eq.3) then
                 call rotate_axes3D(irotateaxes,iplotx,iploty, &
                      lim(ix(:),1),lim(ix(:),2),angleradx,anglerady,angleradz)
              elseif (ndim.eq.2) then
                 call rotate_axes2D(irotateaxes,lim(ix(1:ndim),1), &
                                   lim(ix(1:ndim),2),angleradz)
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
              if (iadapt .and. iadvance.ne.0) then
                 !!--if adaptive limits, find limits of rendered array
                 rendermin = minval(datpix)
                 rendermax = maxval(datpix)
              elseif (iadvance.ne.0) then                   
                 !!--or apply transformations to fixed limits
                 rendermin = lim(irenderplot,1)
                 rendermax = lim(irenderplot,2)
                 call transform_limits(rendermin,rendermax,itrans(irenderplot))
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
              if (iplotpart) then
                 call particleplot(xplot(1:ntot(i)),yplot(1:ntot(i)), &
                   zplot(1:ntot(i)),dat(1:ntot(i),ih,i),ntot(i),iplotx,iploty, &
                   icolourme(1:ntot(i)),npartoftype(:,i), &
                   x_sec,xsecmin,xsecmax,labelz)
              endif
           endif

           !--------------------------------------------------------------
           ! vector maps (can be on top of particle plots and renderings)
           !--------------------------------------------------------------
           if (ivectorplot.ne.0 .and. ndim.ge.2) then
             if (iamvec(ivectorplot).ne.0) then
              !!--choose quantity to be plotted
              ivecx = iamvec(ivectorplot) + iplotx - 1
              ivecy = iamvec(ivectorplot) + iploty - 1

              labelvecplot = trim(labelvec(ivectorplot))
              print*,'plotting vector field ',trim(labelvecplot)
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
                 !!--allocate memory for interpolated vector components
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
                       call interpolate3D_xsec_vec(xplot(1:ninterp), &
                         yplot(1:ninterp),zplot(1:ninterp), &
                         dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i),  &
                         dat(1:ninterp,ih,i),dat(1:ninterp,ivecx,i),dat(1:ninterp,ivecy,i), &
                         ninterp,xmin,ymin,xsecpos, &
                         vecpixx,vecpixy,npixvec,npixyvec,pixwidth)
                    else
                       !
                       !--or interpolate (via averaging) to coarser grid
                       !
                       !call fieldlines2D(ninterp,xplot(1:ninterp),yplot(1:ninterp), &
                       !     dat(1:ninterp,ivecx,i),dat(1:ninterp,ivecy,i), &
                       !     dat(1:ninterp,ih,i),dat(1:ninterp,ipmass,i), &
                       !     dat(1:ninterp,irho,i),xmin,xmax,ymin,ymax)

                       call interpolate_vec(xplot(1:ninterp),yplot(1:ninterp), &
                         dat(1:ninterp,ivecx,i),dat(1:ninterp,ivecy,i), &
                         xmin,ymin,pixwidth,vecpixx,vecpixy, &
                         ninterp,npixvec,npixyvec)                       
                    endif
                    !
                    !--plot it
                    !
                    call render_vec(vecpixx,vecpixy,vecmax, &
                         npixvec,npixyvec,xmin,ymin,pixwidth,labelvecplot)
                    deallocate(vecpixx,vecpixy)
                 endif
                 if (UseBackgndColorVecplot) call pgsci(1)
                endif
              endif
           endif
           !
           !--print legend if this is the first plot on the page
           !    
           if (nyplot.eq.1) then
              call legend(time(i),hposlegend,vposlegend)
           endif
           !
           !--print title if appropriate
           !
           if (i.le.ntitles) then
              if (titlelist(i)(1:1).ne.' ') then
                 call pgmtxt('T',vpostitle,hpostitle,fjusttitle,trim(titlelist(i)))
              endif
           endif
           !
           !--plot exact solution if relevant (before going interactive)
           !
           if (iexact.ne.0) call exact_solution(iplotx,iploty,iexact,ndim,ndimV, &
                            time(i),xmin,xmax,0.0,gamma(i),dat(1:ntot(i),ipmass,i), &
                            ntot(i))
           !
           !--enter interactive mode
           !
           lastplot = (i.eq.n_end .and. nyplot.eq.nyplots .and. k.eq.nxsec)

           if (interactive .and. (nacross*ndown.eq.1)) then
              iadvance = nfreq
              call interactive_part(ninterp,iplotx,iploty,ixsec,irenderplot, &
                   xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                   dat(1:ninterp,ih,i),icolourme(1:ninterp), &
                   xmin,xmax,ymin,ymax,xsecmin,xsecmax,rendermin,rendermax, &
                   angletempx,angletempy,angletempz,ndim,iadvance,isave)
              !--turn rotation on if necessary
              if (abs(angletempx-anglex).gt.tol) irotate = .true.
              if (abs(angletempy-angley).gt.tol) irotate = .true.
              if (abs(angletempz-anglez).gt.tol) irotate = .true.
              if (iadvance.eq.-666) return
           elseif (iplotsonpage.eq.nacross*ndown .or. lastplot) then
              !
              !--timestep control only if multiple plots on page
              !
              iadvance = nfreq
              call interactive_step(iadvance,xmin,xmax,ymin,ymax)
              if (iadvance.eq.-666) return
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

        print*,trim(labely),' min,max = ',ymin,ymax
        print*,trim(labelx),' min,max = ',xmin,xmax

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
                          trim(labelx),trim(labely),' ',just,iaxis)
        else
            !--change the page if pagechange set
           !  or, if turned off, between plots on first page only
           inewpage = ipagechange .or. (iplots.le.nacross*ndown)
           call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
             trim(labelx),trim(labely),' ',just,iaxis, &
             isamexaxis,isameyaxis,inewpage)
        endif

        !--------------------------------
        ! now plot particles
        !--------------------------------

        !--plot time on plot
        if (nyplot.eq.1) call legend(time(i),hposlegend,vposlegend)

        call particleplot(xplot(1:ntot(i)),yplot(1:ntot(i)), &
             zplot(1:ntot(i)),dat(1:ntot(i),ih,i),ntot(i),iplotx,iploty, &
             icolourme(1:ntot(i)),npartoftype(:,i),.false.,0.0,0.0,' ')

        if ((i.eq.nstart).and.iplotlinein) then! plot initial conditions as dotted line
           call pgsls(linestylein)
        endif

        !--plot line joining the particles
        if (iplotline.or.(iplotlinein.and.(i.eq.nstart))) then
           call pgline(npartoftype(1,i),xplot(1:npartoftype(1,i)), &
                       yplot(1:npartoftype(1,i)))     
        endif
        call pgsls(1)! reset 
        call pgsch(charheight)
        !
        !--plot exact solution if relevant
        !
        if (iexact.eq.5 .and.(iploty.eq.iwaveploty).and.(iplotx.eq.iwaveplotx)) then 
           ymean = SUM(yplot(1:npartoftype(1,i)))/REAL(npartoftype(1,i)) 
        else
           ymean = 0.
        endif

        if (iexact.ne.0 .or. (iploty.eq.irho .and. iplotx.eq.ih) .or. &
           (iplotx.eq.irho .and. iploty.eq.ih)) then
           call exact_solution(iplotx,iploty,iexact,ndim,ndimV,  &
                               time(i),xmin,xmax,ymean,gamma(i), &
                               dat(1:ntot(i),ipmass,i),ntot(i))
        endif
        !
        !--enter interactive mode
        !
        lastplot = (i.eq.n_end .and. nyplot.eq.nyplots)

        if (interactive .and. (nacross*ndown.eq.1)) then
           iadvance = nfreq
           call interactive_part(ntot(i),iplotx,iploty,0,irenderplot, &
                xplot(1:ntot(i)),yplot(1:ntot(i)),zplot(1:ntot(i)), &
                dat(1:ntot(i),ih,i),icolourme(1:ntot(i)), &
                xmin,xmax,ymin,ymax,dummymin,dummymax,dummymin,dummymax, &
                angletempx,angletempy,angletempz,ndim,iadvance,isave)
           if (iadvance.eq.-666) return
        elseif (iplotsonpage.eq.nacross*ndown .or. lastplot) then
           !
           !--timestep control only if multiple plots on page
           !
           iadvance = nfreq
           call interactive_step(iadvance,xmin,xmax,ymin,ymax)
           if (iadvance.eq.-666) return
        endif

     elseif (iploty.le.numplot) then! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! additional plots (not plots of particle data - e.g where some additional 
! information is read from a file and plotted on the same page as the 
! particle plots, or where some additional plot is calculated
! from the particle data, such as errors etc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
           labelx = 'frequency'
           labely = 'power'
           if (iadvance.ne.0) then
              xmin = 1./wavelengthmax  ! freq min
              xmax = nfreqspec*xmin
           endif

           if (.not.idisordered .and. iadvance.ne.0) then! interpolate first
              !!--allocate memory for 1D grid (size = 2*npart)
              ngrid = 2*npartoftype(1,i)
              !!--set up 1D grid
              xmingrid = lim(ix(1),1)
              xmaxgrid = lim(ix(1),2)
              dxgrid = (xmaxgrid-xmingrid)/ngrid
              call set_grid1D(xmingrid,dxgrid,ngrid)

              ninterp = ntot(i)
              !!--interpolate to 1D grid  
              call interpolate1D(dat(1:ninterp,ix(1),i), & 
                   dat(1:ninterp,ipmass,i),dat(1:ninterp,irho,i), &
                   dat(1:ninterp,ih,i),dat(1:ninterp,ipowerspecy,i), & 
                   ninterp,xmingrid,datpix1D,ngrid,dxgrid)
              !!--plot interpolated 1D data to check it
              !!print*,minval(datpix1D),maxval(datpix1D)

              !call pgswin(xmin,xmax,minval(datpix1D),maxval(datpix1D),0,1)
              !call pgbox('BCNST',0.0,0,'1BVCNST',0.0,0)      
              !call pglabel('x',label(ipowerspecy),'1D interpolation')
              !call pgline(ngrid,xgrid,datpix1D)
              !read*
              !call pgpage! change page

              if (iadvance.ne.0) then
                 xmin = 1./wavelengthmax  ! freq min
                 xmax = nfreqspec*xmin
              endif
              !!--call power spectrum calculation on the even grid
              call powerspectrum_fourier(ngrid,xgrid,datpix1D,nfreqspec, &
                   xplot(1:nfreqspec),xmin,xmax,yplot(1:nfreqspec))
              if (allocated(datpix1D)) deallocate(datpix1D)              
           elseif (iadvance.ne.0) then
              !!--or else call power spectrum calculation on the particles themselves    
              call powerspectrum_lomb(ntot(i),dat(1:ntot(i),ix(1),i), &
                   dat(1:ntot(i),ipowerspecy,i),nfreqspec, &
                   xplot(1:nfreqspec),xmin,xmax,yplot(1:nfreqspec))
           endif

           if (iadvance.ne.0) then
              ymin = minval(yplot(1:nfreqspec))
              ymax = maxval(yplot(1:nfreqspec))
           endif

           !!--uncomment next few lines to plot wavelengths instead
           labelx = 'wavelength'
           if (iadvance.ne.0) then
              zplot(1:nfreqspec) = 1./xplot(1:nfreqspec)
              xplot(1:nfreqspec) = zplot(1:nfreqspec)
              xmin = minval(xplot(1:nfreqspec))
              xmax = maxval(xplot(1:nfreqspec))
           endif

           if (itrans(iploty).ne.0) then
              call transform(xplot,itrans(iploty))
              labelx = transform_label(labelx,itrans(iploty))

              call transform(yplot,itrans(iploty))
              labely = transform_label(labely,itrans(iploty))
              if (iadvance.ne.0) then
                 call transform_limits(xmin,xmax,itrans(iploty))
                 call transform_limits(ymin,ymax,itrans(iploty))
              endif
           endif

           !--------------------------------------------------------------
           ! output some muff to the screen
           !--------------------------------------------------------------

           print*,trim(labelx),'min,max = ',xmin,xmax
           print*,trim(labely),'min,max = ',ymin,ymax

           !--------------------------------------------------------------
           ! set up pgplot page
           !--------------------------------------------------------------

           iplots = iplots + 1
           iplotsonpage = iplotsonpage + 1
           if (iplotsonpage.gt.nacross*ndown) iplotsonpage = 1
           just = 0

           if (tile_plots) then
              if (iplotsonpage.eq.1 .and. ipagechange) call pgpage
              call danpgtile(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                             trim(labelx),trim(labely),'Power spectrum',just,iaxis)
           else
               !--change the page if pagechange set
              !  or, if turned off, between plots on first page only
              inewpage = ipagechange .or. (iplots.le.nacross*ndown)
              call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                trim(labelx),trim(labely),'Power spectrum', &
                just,iaxis,isamexaxis,isameyaxis,inewpage)
           endif 

           call pgline(nfreqspec,xplot,yplot)

        endif
        !
        !--if this is the first plot on the page, print legend
        !
        if (iplotsonpage.eq.1) call legend(time(i),hposlegend,vposlegend)

        lastplot = (i.eq.n_end)

        if (interactive .and. (iplotsonpage.eq.nacross*ndown .or. lastplot)) then
           iadvance = nfreq
           call interactive_step(iadvance,xmin,xmax,ymin,ymax)
           if (iadvance.eq.-666) return
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if plot not in correct range
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else
        call pgpage! just skip to next plot

     endif   ! ploty = whatever


  enddo over_plots ! over plots per timestep (nyplot)
  
end subroutine plotstep

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
  
!-------------------------------------------------------------------
! works out whether or not we need to read another dump into memory
!-------------------------------------------------------------------
  subroutine get_nextstep(i,ifile)
   use filenames
   implicit none
   integer, intent(inout) :: i,ifile
   integer :: iskipfiles
   
   !
   !--if data is not stored in memory, read next step from file
   !  skip files if necessary. At the moment assumes number of steps in
   !  each file are the same
   !
   if (i.gt.nstepsinfile(ifile)) then
      if (nstepsinfile(i).ge.1) then
         iskipfiles = (i-nstepsinfile(ifile))/nstepsinfile(ifile)
      else
         print*,'*** error in timestepping: file contains zero timesteps'
         iskipfiles = 1
      endif
      if (iskipfiles.gt.1) then
         print*,'skipping ',iskipfiles,' files '
      elseif (iskipfiles.le.0) then
         print*,'error with iskipfiles = ',iskipfiles
         iskipfiles = 1
      endif
      ifile = ifile+iskipfiles
      if (ifile.le.nfiles) then
         call get_data(ifile,.true.)
         print*,'starting at step ',MOD(i-1,nstepsinfile(ifile))+1
         i = MOD(i-1,nstepsinfile(ifile)) + 1
      else
         i=-666
      endif
   elseif (i.lt.1) then
      ifile = ifile-1
      if (ifile.ge.1) then
         iskipfiles = (i-1)/nstepsinfile(ifile)
         if (abs(iskipfiles).gt.0) print*,'skipping back ',abs(iskipfiles),' files'
         ifile = ifile + iskipfiles + 1
         if (ifile.lt.1) ifile = 1
         call get_data(ifile,.true.)
      else
         ifile = 1
      endif
      i = 1
   endif
  
   return
  end subroutine get_nextstep     
        
end subroutine mainloop
