      PROGRAM supersphplot
!---------------------------------------------------------------------------------
!     Plotting utility for SPH data in 1, 2 and 3 dimensions.
!
!     Uses PGPLOT routines to plot graphs and utilises the rendering
!     tools to plot density renderings and vector plots of 2D and 3D data.
!
!     Subroutines as follows (in alphabetical order):
!
!     calc_quantities    : calculates additional quantities from particle data
!     colour_demo        : demonstration of colour schemes for rendering
!     danpgwedg          : my very minor modification of pgwedg
!     exact_mhdshock     : some "exact" solutions for MHD shocks 
!     exact_polytrope    : exact solution for a polytrope
!     exact_rhoh	 : plots exact relation between density and smoothing length
!     exact_sedov        : exact solution for sedov blast wave
!     exact_swave        : exact solution for a linear sound wave
!     exact_toystar      : exact solution for the toy star problem
!     get_render_options : prompts user for options for render plots
!     interpolate2D	 : interpolation of 2D SPH data to 2D grid using SPH kernel     
!     interpolate3D	 : interpolation of 3D SPH data to 3D grid using SPH kernel
!     interpolate3D_fastxsec   : fast cross section through 3D data using SPH kernel
!     interpolate3D_projection : fast projection of 3D data to 2D grid using integrated SPH kernel
!     legend		       : plots legend on plot (time)
!     lomb_powerspectrum1D     : calculates power spectrum of data on particles
!     menu_actions	 : plot options in menu format
!     modules		 : contains all shared (global) variables
!     plot_average	 : bins particles along x-axis and plots average line
!     plot_powerspectrum : calls powerspectrum and plots it
!     print_menu	 : prints menu
!     read_data_dansph   : reads data from my format of data files
!     read_data_mrbsph   : reads data from Matthew Bate's format of data files
!     read_defaults	 : reads default plot options from file
!     render_coarse	 : interpolates to grid by averaging from particles (used for vector plots)
!     render_smooth	 : calls interpolate and plots render maps
!     set_defaults	 : sets default plot options if not read from file
!     setcolours	 : sets up PGPLOT colour table for rendering
!     supersphplot	 : main program, does plots
!     transform	 	 : applies various transformations to data (log10, 1/x, etc)
!     write_defaults	 : writes default plot options to file
!
!     File format is specified in the subroutine read_data   
!
!     Written by: Daniel Price, Institute of Astronomy, Cambridge UK
!          email: dprice@ast.cam.ac.uk
!
!     This version for both NDSPMHD and Matthew Bate's code 2003
!     Changes log:
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
!      25/06/03 - clever Makefile - makes for dansph or M. Bate SPH
!               - subroutines in different files
!      20/06/03 - rendering can handle zero density, prints if min=max
!      19/06/03 - multiplot with rendering, no x array
!      16/06/03 - labels in module, specified in read_data
!      09/06/03 - ndim, ndimV changeable, reads ncolumns in data file
!      26/05/03 - calculated quantities no longer in read_data
!      16/05/03 - colour schemes + rendering using SPH summation
!                 output format has changed, reads pmass array also
!      13/05/03 - read polytrope in menu option
!      01/05/03 - can manually enter plot limits 
!      29/04/03 - bug in initial plot limits fixed
!
!
!      Plots can be of two types: co-ordinate plots or not
!      1) Co-ordinate plots have co-ordinates as x and y axis
!         These plots can be rendered with any scalar or vector array.
!         
!         The rendering routines interpolate from the particles to either
!         a 2D or 3D grid. In 3D you can either render to a 3D grid and take
!         cross sections, or render to a 2D grid using a table of the integrated
!         SPH kernel. This 2D rendering results in a map of the quantity
!         integrated through the third co-ordinate. 
!         Rendering to a 3D grid can be quite slow - it is only efficient
!         if many cross sections are taken all at once from the same data.
!
!      2) Other plots have a variety of options, with lines joining the particles
!         and various exact solutions. Plot limits can be fixed or adaptive.
!
!      Multiplot enables you to set up multiple plots per page, mixing from any type.
!
!      Improvements wishlist: allocatable arrays
!----------------------------------------------------------------------------------
      
      USE params
      USE exact_params
      USE filenames
      USE labels
      USE multiplot
      USE particle_data
      USE prompting	! for shock type only
      USE settings      
      IMPLICIT NONE      
      INTEGER :: nstep,i,j,k,ipix,len
      INTEGER :: ipickx,ipicky,iplotx,iploty
      INTEGER :: ilist,ihoc,ibin,listsize
      INTEGER :: ntotmin
      INTEGER :: nyplot,nyplots      
      INTEGER :: npart1
      INTEGER :: npixy,npixz,ipixxsec
      INTEGER :: ivecplot,npix,npixvec,ncontours
      INTEGER :: irenderprev, istepprev
      INTEGER :: isizex,isizey	! for sending datpix to transform
      INTEGER :: nsink,nsinkstart,nsinkend,nghoststart,nghostend
      INTEGER :: ishk,int_from_string,ichoosey

      CHARACTER(LEN=8) :: string	! used in PGPLOT calls
      REAL, DIMENSION(2,max) :: vecplot
      REAL, DIMENSION(max) :: xplot,yplot,renderplot
      REAL, DIMENSION(:,:), ALLOCATABLE :: datpix
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: datpix3D
      REAL :: xmin,xmax,ymin,ymax,zmin,zmax,xminrender,xmaxrender
      REAL :: vmin,vmax,rendermin,rendermax
      REAL :: scale,rhomax,rhomin,binsize
      REAL :: xsecmin,xsecmax,dxsec,xsecpos
      REAL :: pixwidth
      REAL :: charheight

      LOGICAL :: iplotcont,iplotpartvec,x_sec
      LOGICAL :: log,use_backgnd_color_vecplot      
      
      CHARACTER(LEN=20) :: filename
      CHARACTER(LEN=60) :: title,titlex,datfile
      CHARACTER(LEN=1) :: ans,dummy,logx,logy
      CHARACTER(LEN=20) :: labelx,labely,labelz,labelrender
      CHARACTER(LEN=25) :: transform_label
!
!-----------------------------------------------------------------------------
!
      CALL print_header
      PRINT*,'( version 4.1 )'
!
!--set default options
!
      CALL set_defaults
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
      CALL read_defaults

c ------------------------------------------
c prompt for title

       title = ' '
!       print*,' Enter title for graphs '
!       read(*,808) title
808    FORMAT (a40)       
       print*,'title = ',title

c ---------------------------------------------------
c get rootname from command line/file and read file

      CALL getarg(1, rootname)
      IF (rootname(1:1).NE.' ') THEN
         ihavereadfilename = .true.
         CALL menu_actions(numplot+2)
      ENDIF	 
!------------------------------------------------------------
! setup kernel table for fast column density plots in 3D
      CALL setup_integratedkernel
      
! ----------------------------------------------------------------
! menu - loop back to here once finished plotting/setting options
!
      ipicky = 1
      menuloop: DO WHILE (ipicky.GT.0 .AND. ipicky.LE.numplot+menuitems)

100   CONTINUE

!
!--numplot is the total number of data columns (read + calculated)
!   not including the particle co-ordinates
!  nextra are extra graphs to plot (e.g. convergence plots, power spectrum)
!
      nextra = 0
      IF (ndim.EQ.1) THEN
         nextra = 1	! one extra plot = power spectrum
         ipowerspec = ncolumns + ncalc + 1
         label(ipowerspec) = '1D Power spectrum'
      ENDIF
      IF (iexact.EQ.4) THEN	! toy star plot A-C plane
         nextra = nextra + 1
	 iACplane = ncolumns + ncalc + nextra
	 label(iACplane) = 'A-C plane'
	 IF (magfield) THEN
	    sigma = sigma0
	 ELSE
	    sigma = 0.
	 ENDIF
      ENDIF	 

      IF (ivegotdata) THEN
         numplot = ncolumns + ncalc + nextra
         IF (numplot.GT.maxplot) THEN
            PRINT*,numplot,ncolumns,ncalc,nextra
            STOP ' numplot > array limits, see modules.f'
         ENDIF	 
         ndataplots = ncolumns + ncalc
      ELSE
         numplot = 0
	 ndataplots = 0
      ENDIF
!
!--these are the quantities calculated by the program
!      
      IF (ientrop.NE.0) label(ientrop) = 'entropy'
      IF (irad.NE.0) label(irad) = 'radius '
      IF (irad2.NE.0) label(irad2) = 'r_parallel'
      IF (ike.NE.0) label(ike) = 'specific KE'
      IF (ipr.NE.0) label(ipr) = 'P'	!'P_gas '
      IF (ipmag.NE.0) label(ipmag) = 'P_mag'
      IF (itotpr.NE.0) label(itotpr) = 'P_gas + P_mag'
      IF (ibeta.NE.0) label(ibeta) = 'Plasma \gb'
      IF (idivBerr.NE.0) label(idivBerr) = 'h |div B| / |B|'
      IF (itimestep.NE.0) label(itimestep) = 'h SQRT(\gr) / |B|'
      
!----------------------------------------------------------------------
!  print menu
!     
      CALL print_menu(ipicky,ipickx)
       
!-----------------------------------------------------------------------
! set plot options from menu
!      
      imulti = .false.
      isamexaxis = .true.	! used to determine whether to plot labels
      
      IF ((ipicky.GT.numplot+1).AND.
     &    (ipicky.LE.numplot+menuitems)) THEN
!-------------------------------------------------------
!     adjust plot settings via menu options
!-------------------------------------------------------

         CALL menu_actions(ipicky)
      ELSE 	! do plot 

!------------------------------------------------------
!     or else plot data
!-------------------------------------------------------

!!--if data has not been read from a file, prompt for file
       IF (.not.ivegotdata) THEN
          PRINT*,' No data '
	  ihavereadfilename = .false.
	  ipicky = numplot+2
	  CALL menu_actions(ipicky)
	  GOTO 100	! return to menu
       ENDIF

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
      
      IF (ndim.NE.3) x_sec = .false.
        
!!--set current plot to first in multiplot array if doing multiplot
      nxsec = 1
      IF (ipicky.EQ.numplot+1) THEN	! multiplot
         imulti=.true.
	 IF (ANY(multiplotx(1:nyplotmulti).NE.multiplotx(1))) THEN
	    isamexaxis = .false.
	 ENDIF
         ipickx = multiplotx(1)
	 ipicky = multiploty(1)
	 nyplots = nyplotmulti
      ELSE
         nyplots = 1 
      ENDIF	 
!!--if doing multiplot can only take a single cross section slice      
      IF (imulti) THEN
         flythru = .false.
	 nxsec = 1
      ENDIF	 
      
!------------------------------------------------------------------------
! co-ordinate plot initialisation

      IF (ipicky.le.ndim .and. ipickx.le.ndim) THEN

!!--work out coordinate that is not being plotted	 
	 DO j=1,ndim
            IF ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j	 
         ENDDO	 	       
         IF (ixsec.EQ.0) x_sec = .false.   ! ie can only have x_sec in 3D

!!--if series of cross sections (flythru), set position of first one      
         IF (x_sec.and.flythru) THEN
            PRINT*,'Enter number of cross-section slices:'
	    READ*,nxsec
	    dxsec = (lim(ixsec,2)-lim(ixsec,1))/float(nxsec)
	    xsecmin = lim(ixsec,1)-dxsec
	    xsecmax = lim(ixsec,1)
	    xsecpos = xsecmin
	    xsecpos_nomulti = xsecpos

!!--if particle cross-section, read coordinate range to plot particles in

         ELSEIF (x_sec.AND.iplotpart.AND.irender.LE.ndim) THEN
            PRINT 33,label(ixsec)
33          FORMAT(' Enter ',a1,' min, max for cross section slice:')
            READ*,xsecmin,xsecmax
	    IF (xsecmax.gt.lim(ixsec,2)
     &	    	.or.xsecmin.lt.lim(ixsec,1)) THEN
	        PRINT*,'Warning: Cross section outside data limits' 
	    ENDIF
         ENDIF
	 
!!--set title of plot

        IF ((.not.imulti).AND.(nacross*ndown.EQ.1)) THEN
           IF (x_sec) THEN
	      titlex = 'cross-section'
           ELSE
	      titlex = 'projection'   
           ENDIF
	   titlex = TRIM(label(ipickx))//TRIM(label(ipicky))
     &         //' '//titlex	          

	   IF (irender.GT.ndim) THEN
	   titlex = TRIM(titlex)//' - '//TRIM(label(irender))
     &     //' rendering'
           ENDIF
	   IF (ivecplot.EQ.1) titlex = ' velocity map: '//titlex
	   IF (ivecplot.EQ.2) titlex = ' magnetic field map: '//titlex
	ELSE
	   titlex = ' '
	ENDIF

!!--initialise PGPLOT     
	CALL PGBEGIN(0,'?',nacross,ndown)	
!
!--set colour table
!
	IF (((irender.GT.ndim).OR.ANY(irendermulti(1:nyplots).GT.ndim))
     &	    .AND.(icolours.GT.0)) THEN
           CALL setcolours(icolours)
        ENDIF

!!------------------------------------------------------------------------      
! non- co-ordinate plot initialisations
!
      ELSE

        IF ((nacross.GT.1).OR.(ndown.GT.1)) THEN
           title = '          '
        ENDIF

        CALL PGBEGIN(0,'?',nacross,ndown)

      ENDIF
!!------------------------------------------------------------------------
! general initialisations

!!--set paper size
      IF (nacross.EQ.2 .AND. ndown.EQ.1) THEN
         CALL PGPAPER(11.7,0.5/SQRT(2.))
      ELSEIF (ipapersize.GT.0 .AND. papersizex.GT.0.0
     &        .AND. aspectratio.GT.0.0 ) THEN
         CALL PGPAPER(papersizex,aspectratio)
      ENDIF	 


      IF (animate) CALL PGASK(.false.)

      iploty = ipicky	
      iplotx = ipickx
!
!--if plotting ghost particles, set ntotplot = ntot, else ntot=npart
!
      ntotplot(:) = npart(:)
      IF (iplotghost) ntotplot = npart(:) + nghost(:)
      IF (iplotsink) ntotplot = ntot(:)
!
!--set fill style for circle plots
!      
      IF (plotcirc) CALL PGSFS(2)
!
!--increase character size depending on the number of graphs on the page
!
!      PRINT*,' Enter character height '
!      READ*,charheight
      charheight = 1.0
      IF ((ndown*nacross).GT.1) charheight = 2.0
!      charheight = 0.5*(nacross+ndown)

!------------------------------------------------------------------------      
! loop over timesteps 
!------------------------------------------------------------------------            
      over_timesteps: DO i=nstart,n_end,nfreq

	  npart1 = npart(i) + 1   	 	          
          irenderprev = 0
	  istepprev = 0  
!-------------------------------------
! loop over plots per timestep
!-------------------------------------
          over_plots: DO nyplot=1,nyplots
!--make sure character height is set correctly             
	     CALL PGSCH(charheight)	  
!--for consecutive plots (ie. if not multi but nyplots > 1 plots consecutive numbers)	     
	     iploty = ipicky + nyplot - 1
!--set current x, y plot from multiplot array
             IF (imulti) THEN	        
                iploty = multiploty(nyplot)
		iplotx = multiplotx(nyplot)		
	     ENDIF
!--------------------------------------------------------------
!  copy from main dat array into xplot, yplot 
!  apply transformations (log, 1/x, etc) if appropriate
!--------------------------------------------------------------
             IF (iploty.LE.numplot-nextra 
     &     .AND. iplotx.LE.numplot-nextra) THEN
	        CALL transform(dat(iplotx,:,i),xplot,itrans(iplotx),max)
	        CALL transform(dat(iploty,:,i),yplot,itrans(iploty),max)
	        labelx = transform_label(label(iplotx),itrans(iplotx))
	        labely = transform_label(label(iploty),itrans(iploty))
	        CALL transform(lim(iplotx,1),xmin,itrans(iplotx),1)
	        CALL transform(lim(iplotx,2),xmax,itrans(iplotx),1)
	        CALL transform(lim(iploty,1),ymin,itrans(iploty),1)
	        CALL transform(lim(iploty,2),ymax,itrans(iploty),1)
!--work out whether to use log axes - this is for the call to PGBOX
                logx = ' '
                logy = ' '
                IF (itrans(iplotx).EQ.1) logx = 'L'
	        IF (itrans(iploty).EQ.1) logy = 'L'

!--write username, date on plot
!         IF (nacross.le.2.and.ndown.le.2) CALL PGIDEN

!--adjust plot limits if adaptive plot limits set
	       IF ((ipagechange.AND.iadapt).AND.(iplotx.LE.ndataplots)
     &            .AND.(iploty.LE.ndataplots)) THEN
	          xmin = minval(xplot(1:ntotplot(i)))
	          xmax = maxval(xplot(1:ntotplot(i)))*scalemax
	          ymin = minval(yplot(1:ntotplot(i)))
	          ymax = maxval(yplot(1:ntotplot(i)))*scalemax
	       ENDIF

            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! plots with co-ordinates as x and y axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 IF ((iploty.le.ndim).and.(iplotx.le.ndim)) THEN

!!--set rendering options equal to settings in multiplot	 
	    IF (imulti) THEN
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
	    ELSE
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
	    ENDIF

!!--work out coordinate that is not being plotted	 
	    DO j=1,ndim
               IF ((j.ne.iplotx).and.(j.ne.iploty)) ixsec = j	 
            ENDDO	 	       	   
            IF (ixsec.EQ.0) x_sec = .false.   ! ie can only have x_sec in 3D	   
!!--set limits for rendering area
            xminrender = lim(ix(1),1)
	    xmaxrender = lim(ix(1),2)
	    DO j=1,ndim
	       IF (xminrender.GT.lim(ix(j),1)) xminrender = lim(ix(j),1)
	       IF (xmaxrender.LT.lim(ix(j),2)) xmaxrender = lim(ix(j),2)
	    ENDDO
	    
	    
!------------------------------------------------------------------
!  rendering setup and interpolation (this is the rendering done
!  *before* the cross sections are taken, e.g. to 3D grid)
!------------------------------------------------------------------
             IF ((irenderplot.GT.ndim).AND.(ndim.GE.2)) THEN
!
!--interpolate from particles to fixed grid using SPH summation
!		
!--do not apply any transformations to the co-ordinates
		xmin = lim(ix(1),1)
		xmax = lim(ix(1),2)
		xminrender = xmin
		xmaxrender = xmax
!!--determine number of pixels in rendered image (npix = pixels in x direction)
		pixwidth = (xmax-xmin)/REAL(npix)
                print*,'npix = ',npix
		IF (ndim.GE.2) THEN
		   ymin = lim(ix(2),1)
		   ymax = lim(ix(2),2)
		   xminrender = MIN(xmin,ymin)
		   IF (ymax.GT.xmaxrender) xmaxrender = ymax
		   npixy = INT((ymax-ymin)/pixwidth) + 1
		   print*,'npixy = ',npixy
!!--only need z pixels if working with interpolation to 3D grid
		   IF ((ndim.GE.3).AND.(x_sec.AND.nxsec.GT.2)) THEN
		      zmin = lim(ix(3),1)
		      zmax = lim(ix(3),2)
		      xminrender = MIN(xminrender,zmin)
		      IF (zmax.GT.xmaxrender) xmaxrender = zmax
		      npixz = INT((zmax-zmin)/pixwidth) + 1		 
		      print*,'npixz = ',npixz
		   ENDIF
		   
                ENDIF		
		
!!--if rendering array is the same as the previous plot, reuse the array		
		IF (irenderplot.EQ.irenderprev 
     &		   .AND. i.EQ.istepprev) THEN
		   PRINT*,'same rendering, using previous array...'		
		ELSE   
	         IF (ALLOCATED(datpix)) DEALLOCATE(datpix)		
 	         IF (ALLOCATED(datpix3D)) DEALLOCATE(datpix3D) 
		 SELECT CASE(ndim)
		  CASE(2)
!!--allocate memory for rendering array
                      isizex = npix
		      isizey = npixy
		      ALLOCATE ( datpix(npix,npixy) )
		      CALL interpolate2D(
     &		      dat(ix(1),1:ntot(i),i),dat(ix(2),1:ntot(i),i),
     &		      dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),
     &                dat(ih,1:ntot(i),i),
     &		      dat(irenderplot,1:ntot(i),i),
     &                ntot(i),xminrender,ymin,
     &	              datpix,npix,npixy,pixwidth)
		  CASE(3)
!!--interpolation to 3D grid - then take multiple cross sections/projections		  
!!  do this if taking more than 2 cross sections, otherwise use fast xsec
		      IF (x_sec.AND.nxsec.GT.2) THEN  
!!--allocate memory for 3D rendering array		 
		         ALLOCATE ( datpix3D(npix,npixy,npixz) )
!!--interpolate from particles to 3D grid
		         CALL interpolate3D(
     &		         dat(ix(1),1:ntot(i),i),dat(ix(2),1:ntot(i),i),
     &		         dat(ix(3),1:ntot(i),i),dat(ipmass,1:ntot(i),i),
     &                   dat(irho,1:ntot(i),i),dat(ih,1:ntot(i),i),
     &		         dat(irenderplot,1:ntot(i),i),
     &                   ntot(i),xminrender,ymin,zmin,
     &	                 datpix3D,npix,npixy,npixz,pixwidth)
		      ENDIF
		  END SELECT                 
		ENDIF
		
		irenderprev = irenderplot
		istepprev = i
	    ENDIF
!
!--if vector plot determine whether or not to plot the particles as well
!
            iplotpart = .true.
            IF (ivecplot.GT.0) iplotpart = iplotpartvec     
	    IF (irenderplot.GT.0) iplotpart = .false.

!
!--%%%%%%%%%%%%% loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
!
          over_cross_sections: DO k=1,nxsec

!------------------------------------------------------------
! for multislice cross section (flythru)
! increment the position of the current cross section slice	   
!-------------------------------------------------------------

   	    IF (x_sec.AND.flythru) THEN
!!--for cross sections of particle plots, need range of co-ordinates in which
!!  particles may lie
	       IF (iplotpart) THEN
		  xsecmin = xsecmin+dxsec
	          xsecmax = xsecmax+dxsec
	       ELSE
!!--for cross sections through rendered data, need only co-ordinate of slice	       
		  xsecpos = xsecpos + dxsec
	       ENDIF		   
            ENDIF

!------------------------------------------------------------------------
! preliminary muff for 3D renderings
!------------------------------------------------------------------------
	    IF (irenderplot.GT.ndim .AND. ndim.EQ.3) THEN	    

!!--the line below is the condition for doing a full interpolation to a 3D grid
!!  - need to make sure it is the same everywhere
!!  - could be a bit more clever about when we do this
!!    (at the moment all the projections are done via the fast projection,
!!     although the ability to do projections from the 3D grid is there)
!!
             IF (x_sec .AND. nxsec.GT.2) THEN
!------------------------------------------------------------------------
! if we have rendered to a 3D grid, take cross sections 
! or projections through this grid
!------------------------------------------------------------------------
	       IF (ALLOCATED(datpix)) DEALLOCATE(datpix)
	       SELECT CASE(ixsec)	! which direction to slice in
                 CASE(1)
!!--allocate memory for rendering array
                     isizex = npixy
		     isizey = npixz
		     ALLOCATE ( datpix(npixy,npixz) )
		     IF (x_sec) THEN	! take cross section
                        ipixxsec = INT((xsecpos - xmin)/pixwidth) + 1
			IF (ipixxsec.GT.npix) ipixxsec = npix
                        PRINT*,'x = ',xsecpos,
     &			       ' cross section, pixel ',ipixxsec
                        datpix = datpix3D(ipixxsec,:,:)
		     ELSE  		! take column density 
		        PRINT*,'taking projection through x data...'		     
		        datpix = 0.
		        DO ipix=1,npix
		           datpix(:,:) = datpix(:,:) 
     &			   + datpix3D(ipix,:,:)*pixwidth
		        ENDDO
		     ENDIF
                 CASE(2)
!!--allocate memory for rendering array
                     isizex = npix
		     isizey = npixz
		     ALLOCATE ( datpix(npix,npixz) )
		     IF (x_sec) THEN
                        ipixxsec = INT((xsecpos - ymin)/pixwidth) + 1
			IF (ipixxsec.GT.npixy) ipixxsec = npixy
                        PRINT*,'y = ',xsecpos,
     &			       ' cross section, pixel ',ipixxsec
                        datpix = datpix3D(:,ipixxsec,:)
		     ELSE
		        PRINT*,'taking projection through y data...'		     
		        datpix = 0.
			DO ipix=1,npixy
			   datpix(:,:) = datpix(:,:) 
     &			   + datpix3D(:,ipix,:)*pixwidth
			ENDDO
		     ENDIF		     	
                 CASE(3)
!!--allocate memory for rendering array
                     isizex = npix
		     isizey = npixy
		     ALLOCATE ( datpix(npix,npixy) )	
                     IF (x_sec) THEN
		        ipixxsec = INT((xsecpos - zmin)/pixwidth) + 1
			IF (ipixxsec.GT.npixz) ipixxsec = npixz
                        PRINT*,'z = ',xsecpos,
     &			       ' cross section, pixel ',ipixxsec
                        datpix = datpix3D(:,:,ipixxsec)
		     ELSE
		        PRINT*,'taking projection through z data...'
		        datpix = 0.
			DO ipix=1,npixz			
			   datpix(:,:) = datpix(:,:) 
     &			   + datpix3D(:,:,ipix)*pixwidth
			ENDDO
		     ENDIF  
                 END SELECT	       		  	       	     

	      ELSE
!-------------------------------------------------------------------
!  or do a fast projection/cross section of 3D data to 2D array
!-------------------------------------------------------------------

!!--determine limits of rendering plot
		 xminrender = lim(ix(iplotx),1)
		 xmaxrender = lim(ix(iplotx),2)
		 ymin = lim(ix(iploty),1)
		 ymax = lim(ix(iploty),2)		 
!!--determine number of pixels in rendered image (npix = pixels in x direction)
		 pixwidth = (xmaxrender-xminrender)/REAL(npix)
		 npixy = INT((ymax-ymin)/pixwidth) + 1
                 print*,'npix,npixy = ',npix,npixy
                 isizex = npix
		 isizey = npixy
!!--allocate memory for the 2D array		 
		 IF (ALLOCATED(datpix)) DEALLOCATE(datpix)
                 ALLOCATE ( datpix(npix,npixy) )
!!--do fast cross-section
		 IF (x_sec) THEN
                    PRINT*,TRIM(label(ix(ixsec))),' = ',xsecpos,
     &			       ' : fast cross section'

		    CALL interpolate3D_fastxsec(
     &		         dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i),
     &                   dat(ixsec,1:ntot(i),i),
     &		         dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),
     &                   dat(ih,1:ntot(i),i),
     &		         dat(irenderplot,1:ntot(i),i),
     &                   ntot(i),xminrender,ymin,xsecpos,
     &	                 datpix,npix,npixy,pixwidth)			 
		 ELSE
!!--do fast projection		 
		    CALL interpolate3D_projection(
     &		         dat(iplotx,1:ntot(i),i),dat(iploty,1:ntot(i),i),
     &		         dat(ipmass,1:ntot(i),i),dat(irho,1:ntot(i),i),
     &                   dat(ih,1:ntot(i),i),
     &		         dat(irenderplot,1:ntot(i),i),
     &                   ntot(i),xminrender,ymin,
     &	                 datpix,npix,npixy,pixwidth)
                 ENDIF	      
	      
	      ENDIF ! whether 3D grid or fast renderings
	      	 
	    ENDIF 
!--------------end of preliminary muff for 3D renderings ------------------
	    	    
!-----------------------
! set up PGPLOT page
!-----------------------
	    IF ((ipagechange).OR.((.not.ipagechange).AND.(i.EQ.nstart))) THEN
!	       CALL PGENV(lim(iplotx,1),lim(iplotx,2),
!     &                 lim(iploty,1),lim(iploty,2),1,1)	! 0 for no axes
	       CALL PGPAGE
	       IF (nacross*ndown.GT.1) THEN
!	          IF (imulti) THEN
		     CALL PGSVP(0.2,0.8,0.2,0.98)
!		  ELSE
!		     CALL PGSVP(0.0,1.0,0.0,1.0)
!		  ENDIF   
	       ELSE
	          CALL PGSVP(0.1,0.9,0.1,0.9)
	       ENDIF	  
	       CALL PGWNAD(xmin,xmax,ymin,ymax)	!  pgwnad does equal aspect ratios
!!--plot axes (log if appropriate)
	       CALL PGBOX('BCNST'//logx,0.0,0,'1BVCNST'//logy,0.0,0)	       
            ELSEIF (nyplot.EQ.1) THEN
	       CALL PGPANL(1,1)
	    ELSE
	       CALL PGPAGE
	    ENDIF

!---------------------------------
! set plot limits and label plot
!---------------------------------

!--print plot limits to screen
	    PRINT 34, time(i),i
	    PRINT*,TRIM(labely),'min,max = ',ymin,ymax
	    PRINT*,TRIM(labelx),'min,max = ',xmin,xmax
34          FORMAT (5('-'),' t = ',f8.4,', dump #',i3,1x,10('-'))
    
!--set plot limits	    
	    CALL PGWNAD(xmin,xmax,ymin,ymax)	! pgwnad does equal aspect ratios
!--label plot
	    IF (((nyplots-nyplot).LT.nacross).OR.(.not.isamexaxis)) THEN
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
	      CALL PGLABEL(labelx,' ',titlex)	    
	    ELSE
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
!	      CALL PGLABEL(' ',labely,titlex)
	    ENDIF

	    IF (x_sec.AND.iplotpart) PRINT 35,label(ixsec),xsecmin,
     &                          label(ixsec),xsecmax
35          FORMAT('Cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

!------------------------------
! Now actually plot the data
!------------------------------
!---------------------------------------------------------------
! density/scalar field rendering
! having got our 2D array of gridded data (datpix), we plot it	    
!---------------------------------------------------------------
	    IF (irenderplot.GT.ndim) THEN	    
!!--do transformations on rendered array
c               PRINT*,' SIZE = ',SIZE(datpix(:,1)),SIZE(datpix(1,:))	       
	       CALL transform2(datpix,datpix,itrans(irenderplot),
     &		               isizex,isizey)
	       labelrender = label(irenderplot)
!!--set label for column density plots (2268 or 2412 for integral sign)
	       IF (ndim.EQ.3 .AND..NOT. x_sec) THEN	       	  
	          labelrender = '\(2268) '//TRIM(labelrender)//
     &		                ' d'//TRIM(label(ix(ixsec)))
	       ENDIF
	       labelrender = transform_label(labelrender,
     &		                         itrans(irenderplot))
!!--set log axes for call to render
	       IF (itrans(irenderplot).EQ.1) log = .true.

!!--if adaptive limits, find limits of rendered array		
	       IF (iadapt) THEN
	          rendermin = minval(datpix)
  	          rendermax = maxval(datpix)	
!!--or apply transformations to fixed limits
	       ELSE
		  CALL transform(lim(irenderplot,1),rendermin,
     &		              itrans(irenderplot),1)
		  CALL transform(lim(irenderplot,2),rendermax,
     &		              itrans(irenderplot),1)
               ENDIF
!!--print plot limits to screen
	       PRINT*,TRIM(labelrender),' min, max = ',rendermin,rendermax	       
!!--call subroutine to actually render the image	       
	       CALL render(datpix,rendermin,rendermax,TRIM(labelrender),
     &	               npix,npixy,xmin,ymin,pixwidth,
     &                 icolours,iplotcont,ncontours,log)
	       
	    ELSE
!-----------------------
! particle plots
!-----------------------
!
!--if particle cross section, plot particles only in a defined coordinate range
!
            IF (x_sec.AND.iplotpart) THEN
	       DO j=1,npart(i)
	          IF ((dat(ixsec,j,i).lt.xsecmax)
     &		  .and.(dat(ixsec,j,i).gt.xsecmin)) THEN
     		      CALL PGPT(1,xplot(j),yplot(j),imark)
		   IF (plotcirc) THEN
		      CALL PGCIRC(xplot(j),yplot(j),2.*dat(ih,j,i))
		   ENDIF
		  ENDIF
	       ENDDO
	       DO j=npart1,ntotplot(i)
	          IF ((dat(ixsec,j,i).lt.xsecmax)
     &		  .and.(dat(ixsec,j,i).gt.xsecmin)) THEN
     		      CALL PGPT(1,xplot(j),yplot(j),imarkg)
		  ENDIF	       
	       ENDDO

	    ELSE	     	     
!
!--or simply plot all particles
!
cc--plot particle positions
	       IF (iplotpart) CALL PGPT(npart(i),xplot(1:npart(i)),
     &                                  yplot(1:npart(i)),imark)
cc--plot ghost particles with different marker
               IF (iplotpart .AND. iplotghost .AND. nghost(i).GT.0) THEN
	          nghoststart = npart(i) + 1
	          nghostend = npart(i) + nghost(i)
     	          CALL PGPT(nghost(i),xplot(nghoststart:nghostend),
     &                              yplot(nghoststart:nghostend),imarkg)  	       
               ENDIF
cc--plot circles of interaction (circles of radius 2h around each particle)
	       IF (plotcirc) THEN		  
	         IF (plotcircall) THEN
		  PRINT*,'plotting circles of interaction',npart(i) 
		  DO j=1,npart(i)
		     CALL PGCIRC(xplot(j),yplot(j),2.*dat(ih,j,i))
		  ENDDO
		 ELSE 
		  PRINT*,'plotting circle of interaction',icircpart
		  CALL PGCIRC(xplot(icircpart),yplot(icircpart),2*dat(ih,icircpart,i))
                 ENDIF
	       ENDIF
	       IF (ilabelpart) THEN
cc--plot particle labels
                  PRINT*,'plotting particle labels ',ntotplot(i)
		  DO j=1,ntotplot(i)
		     CALL PGNUMB(j,0,1,string,nc)
		     CALL PGSCH(0.5*charheight)
		     CALL PGTEXT(xplot(j),yplot(j),string(1:nc))
		     CALL PGSCH(charheight)
		  ENDDO
	       ENDIF	! ilabelpart
             
	     ENDIF	! if x_sec else    
	    ENDIF	! if irender
!-----------------------------------------------------------------------------
! sink particles (want these to appear on both particle plots and renderings)
!-----------------------------------------------------------------------------
        
!--plot sink particles with different marker again

	    nsink = ntot(i) - nghost(i) - npart(i)
            IF (iplotsink .AND. nsink.GT.0) THEN
	       nsinkstart = npart(i) + nghost(i) + 1
	       nsinkend = ntot(i)
	       PRINT*,' plotting ',nsink,' sink particles...'
               CALL PGSCI(2)
	       CALL PGSCH(2.*charheight)
	       CALL PGPT(nsink,xplot(nsinkstart:nsinkend),
     &                         yplot(nsinkstart:nsinkend),imarksink) 
               CALL PGSCH(charheight)
               CALL PGSCI(1)
	    ENDIF

!----------------------------
! vector maps
!----------------------------	    
	    	    
cc--velocity vector map
	    IF (ivecplot.EQ.1 .AND. ivx.NE.0) THEN
	       PRINT*,'plotting velocity field'
cc--copy appropriate velocity data to a 2D array
	       DO j=1,ntotplot(i)
	          vecplot(1,j) = dat(iplotx+ivx-1,j,i)
		  vecplot(2,j) = dat(iploty+ivx-1,j,i)
	       ENDDO
	       vmax = lim(iplotx+ivx-1,2)
	       IF (lim(iploty+ivx-1,2).gt.vmax) vmax = lim(iploty+ivx-1,2)
	       vmin = min(lim(iplotx+ivx-1,1),lim(iploty+ivx-1,1))
cc--plot arrows in either background or foreground colour
               IF (use_backgnd_color_vecplot) CALL PGSCI(0)
cc--render to a grid by binning the particles and taking average vx, vy in each cell	       
	       CALL coarse_render(
     &              xplot(1:ntotplot(i)),
     &              yplot(1:ntotplot(i)),
     &              xminrender,xmaxrender,vecplot(:,1:ntotplot(i)),
     &              vmin,vmax,ntotplot(i),npixvec,2,icolours,iplotcont)
               IF (use_backgnd_color_vecplot) CALL PGSCI(1)
cc--old stuff here is to plot arrows on the particles themselves
!	       scale = 0.08*(lim(iploty,2)-lim(iploty,1))
!	       CALL PGSCH(0.35)	! character height (size of arrow head)
!	       DO j=1,ntotplot(i)
!	       CALL PGARRO(yplot(i),xplot(i),
!     &	                   yplot(i)+vecplot(2,i)*scale,
!     &			   xplot(i)+vecplot(1,i)*scale)
!	       ENDDO
!	       CALL PGSCH(1.0)    ! reset character height
	    ELSEIF ((ivecplot.EQ.2).AND.(iBfirst.NE.0)) THEN
cc--plot vector map of magnetic field
	       PRINT*,'plotting magnetic field: ',
     &          label(iBfirst+iplotx-1),label(iBfirst+iploty-1)
	       DO j=1,ntotplot(i)
	          vecplot(1,j) = dat(iBfirst+iplotx-1,j,i)
		  vecplot(2,j) = dat(iBfirst+iploty-1,j,i)
	       ENDDO
               IF (use_backgnd_color_vecplot) CALL PGSCI(0)
	       CALL coarse_render(
     &              xplot(1:ntotplot(i)),
     &              yplot(1:ntotplot(i)),
     &              xminrender,xmaxrender,vecplot(1:2,1:ntotplot(i)),
     &              Bmin,Bmax,ntotplot(i),npixvec,2,icolours,iplotcont)	    
               IF (use_backgnd_color_vecplot) CALL PGSCI(1)
	    ENDIF
!
!--print legend if this is the first plot on the page
!	    
	    IF (nyplot.EQ.1) CALL legend(time(i))	    
	    
!
!--%%%%%%%%%%%%% end loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
!
	    ENDDO over_cross_sections
	    	    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! not both coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------
! set up PGPLOT page
!-----------------------
	 ELSEIF ((iploty.gt.ndim .or. iplotx.gt.ndim)
     &     .AND.(iploty.LE.ndataplots .AND. iplotx.LE.ndataplots)) THEN
	    
	    IF ((ipagechange).OR.
     &	       ((.not.ipagechange).AND.(i.EQ.nstart))) THEN
!	       CALL PGENV(limx(iplotx,1),limx(iplotx,2),
!     &                 lim(iploty,1),lim(iploty,2),0,0)	! 0 for no axes
	       CALL PGPAGE
	       IF ((nacross*ndown).GT.1) THEN
	          IF (imulti) THEN
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ELSE
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ENDIF
	       ELSE
	          CALL PGSVP(0.1,0.9,0.1,0.9)	       
	       ENDIF
	       CALL PGSWIN(xmin,xmax,ymin,ymax)
	       CALL PGBOX('BCNST'//logx,0.0,0,'1BVCNST'//logy,0.0,0)	       	       
            ELSEIF (nyplot.EQ.1) THEN
	       CALL PGPANL(1,1)
	    ELSE
	       CALL PGPAGE
	    ENDIF

!---------------------------------
! set plot limits and label plot
!---------------------------------


!--print plot limits to screen
	    PRINT 34, time(i),i
	    PRINT*,TRIM(labely),'min,max = ',ymin,ymax
	    PRINT*,TRIM(labelx),'min,max = ',xmin,xmax
!
!--set plot limits
!	    
	    CALL PGSWIN(xmin,xmax,ymin,ymax)

	    IF (((nyplots-nyplot).LT.nacross).OR.(.not.isamexaxis)) THEN
!--print x and y labels
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
	      CALL PGLAB(labelx,' ',title)	    
	    ELSE
!--print y labels only	    
	      CALL PGMTXT('L',3.0,0.5,1.0,labely)
!	      CALL PGLAB(' ',labely,title)
	    ENDIF

!--------------------------------
! now plot particles
!--------------------------------
	    
	    IF ((i.EQ.nstart).AND.iplotlinein) THEN	! plot initial conditions as dotted line
	       CALL PGSLS(linestylein)
	    ELSE
!--plot time on plot
	       IF (nyplot.EQ.1) CALL legend(time(i))
!--plot particles
	       CALL PGSLS(1)
	       CALL PGSCH(1.0)	! reset character height before plotting particles
	       CALL PGPT(npart(i),xplot(1:npart(i)),yplot(1:npart(i)),imark)
	    ENDIF      
!--plot line joining the particles
	    IF (iplotline.OR.(iplotlinein.AND.(i.EQ.nstart))) THEN
	       CALL PGLINE(npart(i),xplot(1:npart(i)),yplot(1:npart(i)))     
            ENDIF
!--plot ghost particles with different marker
	    IF (iplotghost.AND.nghost(i).GT.0) THEN
	       nghoststart = npart(i) + 1
	       nghostend = npart(i) + nghost(i)
     	       CALL PGPT(nghost(i),xplot(nghoststart:nghostend),
     &                             yplot(nghoststart:nghostend),imarkg)
            ENDIF
!--plot sink particles with different marker again
	    nsink = ntot(i) - npart(i) - nghost(i)
	    nsinkstart = npart(i) + nghost(i) + 1
	    nsinkend = ntot(i)
            IF (iplotsink .AND. nsink.GT.0) THEN
	       PRINT*,'plotting ',nsink,' sinks...'
               CALL PGPT(nsink,xplot(nsinkstart:nsinkend),
     &                         yplot(nsinkstart:nsinkend),imarksink) 
	    ENDIF     
!--plot circles of interaction (error bar of length 2h on co-ordinate axis)
	    IF (plotcirc) THEN	
!!--on all particles	    	  
	       IF (plotcircall) THEN
	          IF (iplotx.LE.ndim) THEN
 		     PRINT*,'plotting error bars x axis',npart(i) 
		     CALL PGERRB(5,npart(i),xplot(1:npart(i)),
     &		      yplot(1:npart(i)),2.*dat(ih,1:npart(i),i),1.0)
		  ELSEIF (iploty.LE.ndim) THEN
 		     PRINT*,'plotting error bars y axis',npart(i) 
		     CALL PGERRB(6,npart(i),xplot(1:npart(i)),
     &		     yplot(1:npart(i)),2.*dat(ih,1:npart(i),i),1.0)
		  ENDIF
	       ELSE 
!!--only on a specified particle
	          IF (iplotx.LE.ndim) THEN
 		     PRINT*,'plotting error bar x axis',npart(i) 
		     CALL PGERRB(5,1,xplot(icircpart),yplot(icircpart),
     &		                   2.*dat(ih,icircpart,i),1.0)
                  ELSEIF (iploty.LE.ndim) THEN
		     PRINT*,'plotting error bar y axis',icircpart
		     CALL PGERRB(6,1,xplot(icircpart),yplot(icircpart),
     &		                   2.*dat(ih,icircpart,i),1.0)		      
		  ENDIF
               ENDIF
	    ENDIF

            CALL PGSLS(1)	! reset 
            CALL PGSCH(charheight)
!
!--plot average line
!
	    IF (iplotav) CALL plot_average(xplot(1:npart(i)),
     &                        yplot(1:npart(i)),npart(i),nbins)
!
!--plot particle labels
!
	    IF (ilabelpart) THEN
               PRINT*,'plotting particle labels ',ntotplot(i)
	       DO j=1,ntotplot(i)
		  CALL PGNUMB(j,0,1,string,nc)
		  CALL PGSCH(0.5*charheight)
		  CALL PGTEXT(xplot(j),yplot(j),string(1:nc))
		  CALL PGSCH(charheight)
	       ENDDO
	    ENDIF	! ilabelpart

         ELSEIF (iploty.LE.numplot) THEN	! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! additional plots (not plots of particle data - e.g where some additional 
! information is read from a file and plotted on the same page as the 
! particle plots)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!--setup plotting page as normal
!
	    IF ((ipagechange).OR.((.not.ipagechange).AND.(i.EQ.nstart))) THEN
	       CALL PGPAGE
	       IF ((nacross*ndown).GT.1) THEN
	          IF (imulti) THEN
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ELSE
		     CALL PGSVP(0.2,0.99,0.2,0.98)
		  ENDIF
	       ELSE
	          CALL PGSVP(0.1,0.9,0.1,0.9)	       
	       ENDIF
!	       CALL PGSWIN(lim(iplotx,1),lim(iplotx,2),
!     &	                   lim(iploty,1),lim(iploty,2))
!	       CALL PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)	       	       
            ELSEIF (nyplot.EQ.1) THEN
	       CALL PGPANL(1,1)
	    ELSE
	       CALL PGPAGE
	    ENDIF
!
!--then call subroutine to plot the additional plot
!
            ! e.g. call routine to do convergence plot here
	    IF (iexact.EQ.4) THEN
	       CALL exact_toystar_ACplane(Atstar,Ctstar,sigma,gamma(i))
	    ENDIF
!
!--power spectrum plots (uses x and data as yet unspecified)
!
	    IF (iploty.EQ.ipowerspec) THEN 
	    !
            ! prompt for data to take power spectrum of 
	    ! 
	       CALL prompt('Enter data to take power spectrum of',
     &	                   ichoosey,ndim+1,numplot-nextra)
               CALL plot_powerspectrum(npart(i),dat(ix(1),1:npart(i),i),
     &                                 dat(ichoosey,1:npart(i),i))
            ENDIF
!
!--if this is the first plot on the page, print legend
!
	    IF (nyplot.EQ.1) CALL legend(time(i))
	      	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if plot not in correct range
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 ELSE
	    CALL PGPAGE	! just skip to next plot
	
	 ENDIF   ! ploty = whatever

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! plot exact solution on top of the plot already on the page
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SELECT CASE(iexact)
         CASE(1)		! polytrope
	    IF ((iploty.eq.irho).and.(iplotx.eq.irad))
     &	    CALL PGLINE(ipolyc,rad(1:ipolyc),den(1:ipolyc))
	 
	 CASE(2)	! soundwave
	    IF ((iploty.eq.irho).and.(iplotx.eq.1).and.(i.ne.1)) THEN
	       CALL exact_swave(time(i),delta,lambda,gamma(i),
     &              xplot(1:npart(i)),yplot(1:npart(i)),
     &              dat(iutherm,1:npart(i),i),npart(i))
            ENDIF
	    
	 CASE(3)	! Sedov blast wave
	    IF ((iploty.eq.irho).and.(iplotx.eq.irad))
     &	    CALL exact_sedov(time(i),gamma(i),xplot(1:npart(i)),
     &                      yplot(1:npart(i)),npart(i))
	 
	 CASE(4)	! Toy star
!	    totmass = SUM(dat(ipmass,1:npart(i),i))
!	    Htstar = (0.75*totmass)**(2./3.)*Ctstar**(1./3.)
!	    Htstar = 1.0    
!	    PRINT*,' totmass,H,A,C in = ',totmass,Htstar,Atstar,Ctstar
	    IF (iplotx.EQ.1) THEN	! if x axis is x coordinate
	       IF (iploty.EQ.irho) THEN
	          CALL exact_toystar(time(i),gamma(i),
     &		        Htstar,Atstar,Ctstar,sigma,norder,1)
	       ELSEIF (iploty.EQ.ipr) THEN
	          CALL exact_toystar(time(i),gamma(i),
     &			Htstar,Atstar,Ctstar,sigma,norder,2)	       
	       ELSEIF (iploty.EQ.iutherm) THEN
	          CALL exact_toystar(time(i),gamma(i),
     &			Htstar,Atstar,Ctstar,sigma,norder,3)	       
	       ELSEIF (iploty.EQ.ivx) THEN
	          CALL exact_toystar(time(i),gamma(i),
     &			Htstar,Atstar,Ctstar,sigma,norder,4)	       
	       ELSEIF (iploty.EQ.iBfirst+1) THEN
	          CALL exact_toystar(time(i),gamma(i),
     &			Htstar,Atstar,Ctstar,sigma,norder,5)
	       ENDIF
	    ELSEIF (iplotx.EQ.irho) THEN
	       IF (iploty.EQ.iBfirst+1) THEN
	          CALL exact_toystar(time(i),gamma(i),
     &			Htstar,Atstar,Ctstar,sigma,norder,6)	       
	       ENDIF   
	    ENDIF
	    
	    IF (iploty.EQ.iACplane) THEN	! plot point on A-C plane
	       CALL exact_toystar(time(i),gamma(i),
     &	                Htstar,Atstar,Ctstar,sigma,norder,7)
	    ENDIF
	    
	 CASE(5) 	! MHD shock tubes
	    IF (iplotx.EQ.1) THEN
!	       PRINT*,'rootname = ',rootname,rootname(5:5)
!--if not already set, try to determine solution to plot from filename
	       IF (ishk.EQ.0) ishk = int_from_string(rootname(5:5))
!--otherwise prompt for shock type	       
	       IF (ishk.EQ.0) THEN ! prompt
	          CALL prompt('Enter shock solution to plot',ishk,0,6)
	       ENDIF
	       IF (iploty.EQ.irho) THEN
	          CALL exact_mhdshock(1,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.ipr) THEN
	          CALL exact_mhdshock(2,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.ivx) THEN
	          CALL exact_mhdshock(3,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.ivx+1) THEN
	          CALL exact_mhdshock(4,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.ivlast) THEN
	          CALL exact_mhdshock(5,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.iBfirst+1) THEN
	          CALL exact_mhdshock(6,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.iBlast) THEN
	          CALL exact_mhdshock(7,ishk,time(i),gamma(i),xmin,xmax)
	       ELSEIF (iploty.EQ.iutherm) THEN
	          CALL exact_mhdshock(8,ishk,time(i),gamma(i),xmin,xmax)
	       ENDIF  
	    ENDIF   
	 END SELECT
!
!--plot h = (1/rho)^(1/ndim)
!
	 IF ((iploty.EQ.ih).AND.(iplotx.EQ.irho)) THEN
	    CALL exact_rhoh(hfact,ndim)
	 ENDIF
	
         ENDDO over_plots	! over plots per timestep (nyplot)
      
      ENDDO over_timesteps
      
      IF (animate) THEN
         PRINT*,'press RETURN to finish'
	 READ*
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

300   CONTINUE
      CALL PGEND
      
      ENDIF 	! if plot or not
      
      ENDDO menuloop

999   CONTINUE                 
      END
