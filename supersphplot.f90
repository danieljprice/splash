program supersphplot
!---------------------------------------------------------------------------------
!     plotting utility for sph data in 1, 2 and 3 dimensions.
!
!     uses pgplot routines to plot graphs and utilises the rendering
!     tools to plot density renderings and vector plots of 2D and 3D data.
!
!     subroutines as follows (in alphabetical order):
!
!     calc_quantities    : calculates additional quantities from particle data
!     colour_demo        : demonstration of colour schemes for rendering
!     danpgwedg          : my very minor modification of pgwedg
!     exact_mhdshock     : some "exact" solutions for mhd shocks 
!     exact_polytrope    : exact solution for a polytrope
!     exact_rhoh	 : plots exact relation between density and smoothing length
!     exact_sedov        : exact solution for sedov blast wave
!     exact_swave        : exact solution for a linear sound wave
!     exact_toystar      : exact solution for the toy star problem
!     get_render_options : prompts user for options for render plots
!     interpolate1D	 : interpolation of 1D sph data to 1D grid using sph kernel
!     interpolate2D	 : interpolation of 2D sph data to 2D grid using sph kernel     
!     interpolate3D	 : interpolation of 3D sph data to 3D grid using sph kernel
!     interpolate3D_fastxsec   : fast cross section through 3D data using sph kernel
!     interpolate3D_projection : fast projection of 3D data to 2D grid using integrated sph kernel
!     legend		       : plots legend on plot (time)
!     menu_actions	 : plot options in menu format
!     modules		 : contains all shared (global) variables
!     plot_average	 : bins particles along x-axis and plots average line
!     plot_powerspectrum : calls powerspectrum and plots it
!     powerspectrum_fourier : calculates power spectrum of 1D data on ordered pts
!     powerspectrum_lomb : calculates power spectrum of 1D data on disordered pts
!     print_menu	 : prints menu
!     read_data_dansph   : reads data from my format of data files
!     read_data_mrbsph   : reads data from matthew bate's format of data files
!     read_defaults	 : reads default plot options from file
!     render_coarse	 : interpolates to grid by averaging from particles (used for vector plots)
!     render_smooth	 : calls interpolate and plots render maps
!     set_defaults	 : sets default plot options if not read from file
!     setcolours	 : sets up pgplot colour table for rendering
!     supersphplot	 : main program, does plots
!     transform	 	 : applies various transformations to data (log10, 1/x, etc)
!     write_defaults	 : writes default plot options to file
!
!     file format is specified in the subroutine read_data   
!
!     written by: daniel price, institute of astronomy, cambridge uk
!          email: dprice@ast.cam.ac.uk
!
!     this version for both ndspmhd and matthew bate's code 2003
!     changes log:
!      17/12/03 - 1D interpolation and crap power spectrum
!		- some options moved to separate subroutines
!      16/12/03 - labels on particle cross sections
!      15/12/03 - namelist input/output, freeform source in modules
! 		- bug fix in read_data (nghosts) and interpolation routines
!      09/12/03 - power spectrum plotting in 1D
!      24/11/03 - calc_quantities in separate subroutine, rhoh moved
!      28/10/03 - bug fix when no data 
!      14/10/03 - new colour schemes
!      09/10/03 - gamma different betw. timesteps, paper size, menuitems
!      23/09/03 - plot error bars, fast cross-section, multi-transform
!      18/09/03 - exact mhd shocks from papers
!      16/09/03 - sink particle plotting, fixed bug in render_coarse
!      12/09/03 - fast column density plots, several bug fixes
!      11/09/03 - bug in vector plots (xminrender,xmaxrender)
!      11/08/03 - prompting, bug fix in resetting options after multiplot
!      30/07/03 - 3D cross sections/projections
!      29/07/03 - split into more subroutines (menu etc)
!      18/07/03 - transformations (log, 1/x etc)
!      15/07/03 - interpolate2D,3D - much simpler than smooth_pixels
!      25/06/03 - clever makefile - makes for dansph or m. bate sph
!               - subroutines in different files
!      20/06/03 - rendering can handle zero density, prints if min=max
!      19/06/03 - multiplot with rendering, no x array
!      16/06/03 - labels in module, specified in read_data
!      09/06/03 - ndim, ndimv changeable, reads ncolumns in data file
!      26/05/03 - calculated quantities no longer in read_data
!      16/05/03 - colour schemes + rendering using sph summation
!                 output format has changed, reads pmass array also
!      13/05/03 - read polytrope in menu option
!      01/05/03 - can manually enter plot limits 
!      29/04/03 - bug in initial plot limits fixed
!
!
!      plots can be of two types: co-ordinate plots or not
!      1) co-ordinate plots have co-ordinates as x and y axis
!         these plots can be rendered with any scalar or vector array.
!         
!         the rendering routines interpolate from the particles to either
!         a 2D or 3D grid. in 3D you can either render to a 3D grid and take
!         cross sections, or render to a 2D grid using a table of the integrated
!         sph kernel. this 2D rendering results in a map of the quantity
!         integrated through the third co-ordinate. 
!         rendering to a 3D grid can be quite slow - it is only efficient
!         if many cross sections are taken all at once from the same data.
!
!      2) other plots have a variety of options, with lines joining the particles
!         and various exact solutions. plot limits can be fixed or adaptive.
!
!      multiplot enables you to set up multiple plots per page, mixing from any type.
!
!      improvements wishlist: allocatable arrays
!----------------------------------------------------------------------------------
      
  use params
  use exact_params
  use filenames
  use labels
  use multiplot
  use particle_data
  use prompting
  use settings      
  implicit none      
  integer :: nstep,i,j,k,ipix,len
  integer :: ipickx,ipicky,iplotx,iploty
  integer :: ilist,ihoc,ibin,listsize
  integer :: ntotmin
  integer :: nyplot,nyplots      
  integer :: npart1
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: ivecplot,npix,npixvec,ncontours
  integer :: irenderprev, istepprev
  integer :: isizex,isizey	! for sending datpix to transform
  integer :: nsink,nsinkstart,nsinkend,nghoststart,nghostend
  integer :: ishk,int_from_string
  integer :: igrid, ngrid
  
  character(len=8) :: string	! used in pgplot calls
  real, dimension(2,max) :: vecplot
  real, dimension(max) :: xplot,yplot,renderplot
  real, dimension(:), allocatable :: datpix1D, xgrid
  real, dimension(:,:), allocatable :: datpix
  real, dimension(:,:,:), allocatable :: datpix3D
  real :: xmin,xmax,ymin,ymax,zmin,zmax,xminrender,xmaxrender
  real :: vmin,vmax,rendermin,rendermax
  real :: scale,rhomax,rhomin,binsize
  real :: xsecmin,xsecmax,dxsec,xsecpos
  real :: pixwidth
  real :: charheight
  real :: dxgrid
  
  logical :: iplotcont,iplotpartvec,x_sec
  logical :: log,use_backgnd_color_vecplot
  
  character(len=20) :: filename
  character(len=60) :: title,titlex,datfile
  character(len=1) :: ans,dummy,logx,logy
  character(len=20) :: labelx,labely,labelz,labelrender
  character(len=25) :: transform_label
  !
  !-----------------------------------------------------------------------------
  !
  call print_header
  print*,'( version 4.2 )'
  !
  !--set default options
  !
  call set_defaults
  !
  !--initialise variables
  !      
  rootname = 'blank'
  xmin = 0.0
  xmax = -1.
  charheight = 1.0
  len = 7
  nyplots = 1
  labelcoord(1) = 'x'
  labelcoord(2) = 'y'
  labelcoord(3) = 'z'
  
  isamexaxis = .true.
  irenderplot = 0
  ishowopts = .false.
  ivegotdata = .false.
  
  ! ---------------------------------------------
  ! read default options from file if it exists
  !
  call read_defaults
  
  ! ------------------------------------------
  ! prompt for title
  
  title = ' '
  !       print*,' enter title for graphs '
  !       read(*,808) title
808 format (a40)       
  print*,'title = ',title

  ! ---------------------------------------------------
  ! get rootname from command line/file and read file
  
  call getarg(1, rootname)
  if (rootname(1:1).ne.' ') then
     ihavereadfilename = .true.
     call menu_actions(numplot+2)
  endif
  !------------------------------------------------------------
  ! setup kernel table for fast column density plots in 3D
  call setup_integratedkernel
  
  ! ----------------------------------------------------------------
  ! menu - loop back to here once finished plotting/setting options
  !
  ipicky = 1
  menuloop: do while (ipicky.gt.0 .and. ipicky.le.numplot+menuitems)
     
100  continue
     
     !
     !--numplot is the total number of data columns (read + calculated)
     !   not including the particle co-ordinates
     !  nextra are extra graphs to plot (e.g. convergence plots, power spectrum)
     !
     nextra = 0
     if (ndim.eq.1) then
        nextra = 1	! one extra plot = power spectrum
        ipowerspec = ncolumns + ncalc + 1
        label(ipowerspec) = '1D power spectrum'
     endif
     if (iexact.eq.4) then	! toy star plot a-c plane
        nextra = nextra + 1
        iacplane = ncolumns + ncalc + nextra
        label(iacplane) = 'a-c plane'
        if (magfield) then
           sigma = sigma0
        else
           sigma = 0.
        endif
     endif
     
     if (ivegotdata) then
        numplot = ncolumns + ncalc + nextra
        if (numplot.gt.maxplot) then
           print*,numplot,ncolumns,ncalc,nextra
           stop ' numplot > array limits, see modules.f'
        endif
        ndataplots = ncolumns + ncalc
     else
        numplot = 0
        ndataplots = 0
     endif
     !
     !--these are the quantities calculated by the program
     !      
     if (ientrop.ne.0) label(ientrop) = 'entropy'
     if (irad.ne.0) label(irad) = 'radius '
     if (irad2.ne.0) label(irad2) = 'r_parallel'
     if (ike.ne.0) label(ike) = 'specific KE'
     if (ipr.ne.0) label(ipr) = 'P'	!'p_gas '
     if (ipmag.ne.0) label(ipmag) = 'P_mag'
     if (itotpr.ne.0) label(itotpr) = 'P_gas + P_mag'
     if (ibeta.ne.0) label(ibeta) = 'plasma \gb'
     if (idivberr.ne.0) label(idivberr) = 'h |div B| / |B|'
     if (itimestep.ne.0) label(itimestep) = 'h sqrt(\gr) / |B|'
     
     !----------------------------------------------------------------------
     !  print menu
     !     
     call print_menu(ipicky,ipickx)
     
     !-----------------------------------------------------------------------
     ! set plot options from menu
     !      
     imulti = .false.
     isamexaxis = .true.	! used to determine whether to plot labels
     
     if ((ipicky.gt.numplot+1).and.(ipicky.le.numplot+menuitems)) then
        !-------------------------------------------------------
        !     adjust plot settings via menu options
        !-------------------------------------------------------
        
        call menu_actions(ipicky)
     else 	! do plot 
        
        !------------------------------------------------------
        !     or else plot data
        !-------------------------------------------------------
        
        !!--if data has not been read from a file, prompt for file
        if (.not.ivegotdata) then
           print*,' no data '
           ihavereadfilename = .false.
           ipicky = numplot+2
           call menu_actions(ipicky)
           goto 100	! return to menu
        endif
        
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
        if (ipicky.eq.numplot+1) then	! multiplot
           imulti=.true.
           if (any(multiplotx(1:nyplotmulti).ne.multiplotx(1))) then
              isamexaxis = .false.
           endif
           ipickx = multiplotx(1)
           ipicky = multiploty(1)
           nyplots = nyplotmulti
        else
           nyplots = 1 
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
           do j=1,ndim
              if ((j.ne.ipickx).and.(j.ne.ipicky)) ixsec = j	 
           enddo
           if (ixsec.eq.0) x_sec = .false.   ! ie can only have x_sec in 3D
           
           !!--if series of cross sections (flythru), set position of first one      
           if (x_sec.and.flythru) then
              print 32,label(ixsec)
32            format('enter number of ',a1,' cross-section slices')
              read*,nxsec
              !!--dxsec is the distance between slices	    
              dxsec = (lim(ixsec,2)-lim(ixsec,1))/float(nxsec)
              xsecpos = lim(ixsec,1) - 0.5*dxsec
              xsecpos_nomulti = xsecpos
              
              !!--if single cross-section, read position of cross-section slice
              
           elseif (x_sec.and.iplotpart.and.irender.le.ndim) then
              print 33,label(ixsec)
33            format(' enter ',a1,' position for cross section slice:')
              read*,xsecpos
              if (xsecpos.gt.lim(ixsec,2).or.xsecpos.lt.lim(ixsec,1)) then
                 print*,'warning: cross section outside data limits' 
              endif
              xsecpos_nomulti = xsecpos
           endif
           
           !!--set title of plot
           
           if ((.not.imulti).and.(nacross*ndown.eq.1)) then
              if (x_sec) then
                 titlex = 'cross-section'
              else
                 titlex = 'projection'   
              endif
              titlex = trim(label(ipickx))//trim(label(ipicky))//' '//titlex	          
              
              if (irender.gt.ndim) then
                 titlex = trim(titlex)//' - '//trim(label(irender))//' rendering'
              endif
              if (ivecplot.eq.1) titlex = ' velocity map: '//titlex
              if (ivecplot.eq.2) titlex = ' magnetic field map: '//titlex
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
              call setcolours(icolours)
           endif
           
           !!------------------------------------------------------------------------      
           ! non- co-ordinate plot initialisations
           !
        else
           !!--prompt for options if plotting power spectrum      
           if (ipicky.eq.ipowerspec) call options_powerspec
           !!--no title if more than one plot on the page
           if ((nacross.gt.1).or.(ndown.gt.1)) title = '          '
           
           call pgbegin(0,'?',nacross,ndown)	!  initialise PGPLOT
           
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
        
        iploty = ipicky	
        iplotx = ipickx
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
                 call transform(dat(iplotx,:,i),xplot,itrans(iplotx),max)
                 call transform(dat(iploty,:,i),yplot,itrans(iploty),max)
                 !--set axis labels, applying transformation if appropriate
                 labelx = transform_label(label(iplotx),itrans(iplotx))
                 labely = transform_label(label(iploty),itrans(iploty))
                 !--set x,y plot limits, applying transformation if appropriate
                 call transform(lim(iplotx,1),xmin,itrans(iplotx),1)
                 call transform(lim(iplotx,2),xmax,itrans(iplotx),1)
                 call transform(lim(iploty,1),ymin,itrans(iploty),1)
                 call transform(lim(iploty,2),ymax,itrans(iploty),1)
                 !--work out whether to use log axes - this is for the call to pgbox
                 logx = ' '
                 logy = ' '
                 if (itrans(iplotx).eq.1) logx = 'l'
                 if (itrans(iploty).eq.1) logy = 'l'
                 
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
                 
                 !!--work out coordinate that is not being plotted	 
                 do j=1,ndim
                    if ((j.ne.ipickx).and.(j.ne.ipicky)) ixsec = j	 
                 enddo
                 if (ixsec.eq.0) x_sec = .false.   ! ie can only have x_sec in 3D	   
                 !!--set limits for rendering area
                 xminrender = MINVAL(lim(ix(1:ndim),1))
                 xmaxrender = MAXVAL(lim(ix(1:ndim),2))	    
                 
                 !------------------------------------------------------------------
                 !  rendering setup and interpolation (this is the rendering done
                 !  *before* the cross sections are taken, e.g. to 3D grid)
                 !------------------------------------------------------------------
                 if ((irenderplot.gt.ndim).and.(ndim.ge.2)) then
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
                          !!--allocate memory for rendering array
                          isizex = npixx
                          isizey = npixy
                          allocate ( datpix(npixx,npixy) )
                          call interpolate2D( &
                               dat(ix(1),1:ntot(i),i),dat(ix(2),1:ntot(i),i), &
                               dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i), &
                               dat(ih,1:ntot(i),i),dat(irenderplot,1:ntot(i),i), &
                               ntot(i),xmin,ymin,datpix,npixx,npixy,pixwidth)
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
                    
                    !---------------preliminary muff for 3D renderings-----------------------!
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
                          datpix = datpix3D(:,:,ipixxsec)	! slices are in 3rd dimension
                          
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
                       
                    endif ! 3D and rendering
                    !--------------end of preliminary muff for 3D renderings ------------------
	    	    
                    !-----------------------
                    ! set up pgplot page
                    !-----------------------
                    if ((ipagechange).or.((.not.ipagechange).and.(i.eq.nstart))) then
                       !	       call pgenv(lim(iplotx,1),lim(iplotx,2), &
                       !                     lim(iploty,1),lim(iploty,2),1,1)	! 0 for no axes
                       call pgpage
                       if (nacross*ndown.gt.1) then
                          !	          if (imulti) then
                          call pgsvp(0.2,0.8,0.2,0.98)
                          !		  else
                          !		     call pgsvp(0.0,1.0,0.0,1.0)
                          !		  endif   
                       else
                          call pgsvp(0.1,0.9,0.1,0.9)
                       endif
                       call pgwnad(xmin,xmax,ymin,ymax)	!  pgwnad does equal aspect ratios
                       !	       call pgswin(xmin,xmax,ymin,ymax)	!  not equal	
                       !!--plot axes (log if appropriate)
                       call pgbox('bcnst'//logx,0.0,0,'1bvcnst'//logy,0.0,0)	       
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
34                  format (5('-'),' t = ',f8.4,', dump #',i3,1x,10('-'))
                    if (x_sec.and.iplotpart) print 35,label(ixsec),xsecmin,label(ixsec),xsecmax
35                  format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)
                    
                    !--set plot limits	    
                    call pgwnad(xmin,xmax,ymin,ymax)	! pgwnad does equal aspect ratios
                    !--label plot
                    if (((nyplots-nyplot).lt.nacross).or.(.not.isamexaxis)) then
                       call pgmtxt('l',3.0,0.5,1.0,labely)
                       call pglabel(labelx,' ',titlex)	    
                    else
                       call pgmtxt('l',3.0,0.5,1.0,labely)
                       !	      call pglabel(' ',labely,titlex)
                    endif
                    
                    !------------------------------
                    ! now actually plot the data
                    !------------------------------
                    !---------------------------------------------------------------
                    ! density/scalar field rendering
                    ! having got our 2D array of gridded data (datpix), we plot it	    
                    !---------------------------------------------------------------
                    if (irenderplot.gt.ndim) then	    
                       !!--do transformations on rendered array	       
                       call transform2(datpix,datpix,itrans(irenderplot),isizex,isizey)
                       labelrender = label(irenderplot)
                       !!--set label for column density (projection) plots (2268 or 2412 for integral sign)
                       if (ndim.eq.3 .and..not. x_sec) then	       	  
                          labelrender = '\(2268) '//trim(labelrender)//' d'//trim(label(ix(ixsec)))
                       endif
                       labelrender = transform_label(labelrender,itrans(irenderplot))
                       !!--set log axes for call to render
                       if (itrans(irenderplot).eq.1) log = .true.
                       
                       !!--if adaptive limits, find limits of rendered array		
                       if (iadapt) then
                          rendermin = minval(datpix)
                          rendermax = maxval(datpix)	
                          !!--or apply transformations to fixed limits
                       else
                          call transform(lim(irenderplot,1),rendermin,itrans(irenderplot),1)
                          call transform(lim(irenderplot,2),rendermax,itrans(irenderplot),1)
                       endif
                       !!--print plot limits to screen
                       print*,trim(labelrender),' min, max = ',rendermin,rendermax	       
                       !!--call subroutine to actually render the image	       
                       call render(datpix,rendermin,rendermax,trim(labelrender),  &
                            npixx,npixy,xmin,ymin,pixwidth,    &
                            icolours,iplotcont,ncontours,log)
                       
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
                                endif	! ilabelpart
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
                                endif	! ilabelpart
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
                                print*,'plotting circle of interaction',icircpart
                                call pgcirc(xplot(icircpart),yplot(icircpart),2*dat(ih,icircpart,i))
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
                          endif	! ilabelpart
             
                       endif	! if x_sec else    
                    endif	! if irender
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
                    if (ivecplot.eq.1 .and. ivx.ne.0) then
                       print*,'plotting velocity field'
                       !!--copy appropriate velocity data to a 2D array
                       do j=1,ntotplot(i)
                          vecplot(1,j) = dat(iplotx+ivx-1,j,i)
                          vecplot(2,j) = dat(iploty+ivx-1,j,i)
                       enddo
                       vmax = lim(iplotx+ivx-1,2)
                       if (lim(iploty+ivx-1,2).gt.vmax) vmax = lim(iploty+ivx-1,2)
                       vmin = min(lim(iplotx+ivx-1,1),lim(iploty+ivx-1,1))
                       !!--plot arrows in either background or foreground colour
                       if (use_backgnd_color_vecplot) call pgsci(0)
                       !!--render to a grid by binning particles and taking average vx, vy in each cell      
                       call coarse_render(xplot(1:ntotplot(i)),yplot(1:ntotplot(i)),  &
                            xminrender,xmaxrender,vecplot(:,1:ntotplot(i)), &
                            vmin,vmax,ntotplot(i),npixvec,2,icolours,iplotcont)
                       if (use_backgnd_color_vecplot) call pgsci(1)
                       !!--old stuff here is to plot arrows on the particles themselves
                       !scale = 0.08*(lim(iploty,2)-lim(iploty,1))
                       !call pgsch(0.35)	! character height (size of arrow head)
                       !do j=1,ntotplot(i)
                       !   call pgarro(yplot(i),xplot(i), &
                       !   yplot(i)+vecplot(2,i)*scale,   &
                       !   xplot(i)+vecplot(1,i)*scale)
                       !enddo
                       !	       call pgsch(1.0)    ! reset character height
                    elseif ((ivecplot.eq.2).and.(ibfirst.ne.0)) then
                       !!--plot vector map of magnetic field
                       print*,'plotting magnetic field: ', &
                            label(ibfirst+iplotx-1),label(ibfirst+iploty-1)
                       do j=1,ntotplot(i)
                          vecplot(1,j) = dat(ibfirst+iplotx-1,j,i)
                          vecplot(2,j) = dat(ibfirst+iploty-1,j,i)
                       enddo
                       if (use_backgnd_color_vecplot) call pgsci(0)
                       call coarse_render(xplot(1:ntotplot(i)),yplot(1:ntotplot(i)), &
                            xminrender,xmaxrender,vecplot(1:2,1:ntotplot(i)), &
                            bmin,bmax,ntotplot(i),npixvec,2,icolours,iplotcont)	    
                       if (use_backgnd_color_vecplot) call pgsci(1)
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
                    !            lim(iploty,1),lim(iploty,2),0,0)	! 0 for no axes
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
                    call pgswin(xmin,xmax,ymin,ymax)
                    call pgbox('bcnst'//logx,0.0,0,'1bvcnst'//logy,0.0,0)	       	       
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
                 
                 if (((nyplots-nyplot).lt.nacross).or.(.not.isamexaxis)) then
                    !--print x and y labels
                    call pgmtxt('l',3.0,0.5,1.0,labely)
                    call pglab(labelx,' ',title)	    
                 else
                    !--print y labels only	    
                    call pgmtxt('l',3.0,0.5,1.0,labely)
                    !	      call pglab(' ',labely,title)
                 endif
                 
                 !--------------------------------
                 ! now plot particles
                 !--------------------------------
                 
                 if ((i.eq.nstart).and.iplotlinein) then	! plot initial conditions as dotted line
                    call pgsls(linestylein)
                 else
                    !--plot time on plot
                    if (nyplot.eq.1) call legend(time(i))
                    !--plot particles
                    call pgsls(1)
                    call pgsch(1.0)	! reset character height before plotting particles
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
                       !!--only on a specified particle
                       if (iplotx.le.ndim) then
                          print*,'plotting error bar x axis',npart(i) 
                          call pgerrb(5,1,xplot(icircpart),yplot(icircpart), &
                               2.*dat(ih,icircpart,i),1.0)
                       elseif (iploty.le.ndim) then
                          print*,'plotting error bar y axis',icircpart
                          call pgerrb(6,1,xplot(icircpart),yplot(icircpart), &
                               2.*dat(ih,icircpart,i),1.0)		      
                       endif
                    endif
                 endif

                 call pgsls(1)	! reset 
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
                 endif	! ilabelpart
                 
              elseif (iploty.le.numplot) then	! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! additional plots (not plots of particle data - e.g where some additional 
                 ! information is read from a file and plotted on the same page as the 
                 ! particle plots, or where some additional plot is calculated
                 ! from the particle data)
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
                    
                    if (.not.idisordered) then	! interpolate first
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
                       !call pgpage	! change page
                       
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
                 call pgpage	! just skip to next plot
                 
              endif   ! ploty = whatever
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! plot exact solution on top of the plot already on the page
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              select case(iexact)
              case(1)		! polytrope
                 if ((iploty.eq.irho).and.(iplotx.eq.irad)) &
                      call pgline(ipolyc,rad(1:ipolyc),den(1:ipolyc))
                 
              case(2)	! soundwave
                 if ((iploty.eq.irho).and.(iplotx.eq.1).and.(i.ne.1)) then
                    call exact_swave(time(i),delta,lambda,gamma(i), &
                         xplot(1:npart(i)),yplot(1:npart(i)), &
                         dat(iutherm,1:npart(i),i),npart(i))
                 endif
	    
              case(3)	! sedov blast wave
                 if ((iploty.eq.irho).and.(iplotx.eq.irad)) &
                      call exact_sedov(time(i),gamma(i),xplot(1:npart(i)), &
                      yplot(1:npart(i)),npart(i))
                 
              case(4)	! toy star
                 !	    totmass = sum(dat(ipmass,1:npart(i),i))
                 !	    htstar = (0.75*totmass)**(2./3.)*ctstar**(1./3.)
                 !	    htstar = 1.0    
                 !	    print*,' totmass,h,a,c in = ',totmass,htstar,atstar,ctstar
                 if (iplotx.eq.1) then	! if x axis is x coordinate
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
                 
                 if (iploty.eq.iacplane) then	! plot point on a-c plane
                    call exact_toystar(time(i),gamma(i),htstar,atstar,ctstar,sigma,norder,7)
                 endif
                 
              case(5) 	! mhd shock tubes
                 if (iplotx.eq.1) then
                    !	       print*,'rootname = ',rootname,rootname(5:5)
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
              
           enddo over_plots	! over plots per timestep (nyplot)
      
        enddo over_timesteps
      
        if (animate) then
           print*,'press return to finish'
           read*
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

300     continue
        call pgend
        
     endif 	! if plot or not
     
  enddo menuloop

999 continue                 
end program supersphplot
