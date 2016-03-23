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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!------------------------------------------------------------------------
!
! This is the core routine for the whole code.
! Drives the plotting pipeline, ie. calls all the routines which
! do the work.
!
! I am gradually trying to make this routine more modular...
!
!------------------------------------------------------------------------
module timestep_plotting
  use params, only:maxplot
  implicit none

  integer, private :: ninterp
  integer, private :: iplotx,iploty,iplotz,irender,irenderplot,icontourplot
  integer, private :: ivectorplot,ivecx,ivecy
  integer, private :: nyplots,npartdim,nyplotfirstonpage,ifirststeponpage
  integer, private :: ngrid,nframefirstonpage
  integer, private :: just,ntitles,nsteplegendlines
  integer, private :: iplots,ipanel
  integer, private :: iframesave
  integer, private :: npixx,npixy,npixz

  real, dimension(:),     allocatable, private :: datpix1D, xgrid
  real, dimension(:,:),   allocatable, private :: datpix,datpixcont,brightness
  real, dimension(:,:,:), allocatable, private :: datpix3D,datpixcont3D
  real, private :: xmin,xmax,ymin,ymax,zmin
  real, private :: rendermin,rendermax,vecmax,contmin,contmax
  real, private :: dz,zslicepos,zobservertemp,dzscreentemp,taupartdepthtemp,rkappafac
  real, private :: dxgrid,xmingrid,xmaxgrid
  real, private :: angletempx, angletempy, angletempz
  !--buffer for interactive mode on multiplots
  integer, dimension(maxplot) :: iplotxtemp,iplotytemp,irendertemp,icontourtemp,ivecplottemp
  real,    dimension(maxplot) :: xminmulti,xmaxmulti,xminadapt,xmaxadapt
  real,    dimension(maxplot) :: vptxmin,vptxmax,vptymin,vptymax,barwmulti
  real, private :: xminadapti,xmaxadapti,yminadapti,ymaxadapti,renderminadapt,rendermaxadapt
  real, private :: contminadapt,contmaxadapt
  real, parameter, private :: pi = 3.1415926536

  logical, private :: iplotpart,iplotcont,x_sec,isamexaxis,isameyaxis,iamrendering,idoingvecplot
  logical, private :: inewpage, tile_plots, lastplot
  logical, private :: imulti,irerender,iAllowspaceforcolourbar,ihavesetweights
  logical, private :: interactivereplot,ihavesetcolours,vectordevice,gotcontours
  logical, private :: OneColourBarPerRow

  public :: initialise_plotting, plotstep
  private

contains

!
! initialise plotting options
! called once for all steps
!
subroutine initialise_plotting(ipicky,ipickx,irender_nomulti,icontour_nomulti,ivecplot)
  use params
  use colours,            only:colour_set
  use labels,             only:label,ipowerspec,ih,ipmass,irho,iamvec,isurfdens,&
                               is_coord,itoomre,iutherm,ipdf,ix,icolpixmap
  use limits,             only:lim,rangeset
  use multiplot,          only:multiplotx,multiploty,irendermulti,icontourmulti, &
                               nyplotmulti,x_secmulti,ivecplotmulti
  use prompting,          only:prompt
  use titles,             only:read_titles,read_steplegend
  use settings_data,      only:ndim,ndimV,numplot,ncolumns,ncalc,ndataplots,required, &
                               icoords,icoordsnew,debugmode,ntypes,usetypeinrenderings
  use settings_page,      only:nacross,ndown,ipapersize,tile,papersizex,aspectratio,&
                               iPageColours,iadapt,iadaptcoords,linewidth,device,nomenu,&
                               interactive,ipapersizeunits,usecolumnorder
  use pagecolours,        only:set_pagecolours
  use settings_part,      only:linecolourthisstep,linecolour,linestylethisstep,linestyle,iexact,iplotpartoftype
  use settings_render,    only:icolours,iplotcont_nomulti,iColourBarStyle,icolour_particles
  use settings_xsecrot,   only:xsec_nomulti,xsecpos_nomulti,flythru,nxsec,irotate, &
                               xseclineX1,xseclineX2,xseclineY1,xseclineY2,xsecwidth, &
                               use3Dperspective,use3Dopacityrendering,zobserver,dzscreenfromobserver,taupartdepth
  use settings_powerspec, only:options_powerspec,options_pdf
  use particle_data,      only:npartoftype,masstype
  use projections3D,      only:coltable
  use plotlib,            only:plot_init,plot_qcur,plot_slw,plot_env,plot_curs,plot_band, &
                               plot_close,plot_qinf
  use system_utils,       only:renvironment
  use calcquantities,     only:get_calc_data_dependencies
  implicit none
  real, parameter     :: pi=3.1415926536
  integer, intent(in) :: ipicky,ipickx,irender_nomulti,icontour_nomulti,ivecplot
  integer             :: i,j,ifirst,iplotzprev,ilen,ierr,irow,ntries
  logical             :: iadapting,icoordplot,iallrendered,ians
  real                :: hav,pmassav,dzsuggest
  integer, dimension(:), allocatable :: ifirstinrow
  character(len=1)    :: char
  character(len=20)   :: devstring

  !------------------------------------------------------------------------
  ! initialisations
  ! should initialise all saved variables here
  !------------------------------------------------------------------------

  isamexaxis = .true.  ! same x axis on all plots? (only relevant for >1 plots per page)
  isameyaxis = .true.  ! same y axis on all plots?
  tile_plots = .false.
  iplots = 0 ! counter for how many plots have been plotted in total
  ipanel = 0  ! counter for which panel we are in on plotting page
  irender     = 0
  irenderplot = 0
  ivectorplot = 0
  x_sec = xsec_nomulti
  iplotcont = iplotcont_nomulti
  lastplot = .false.
  iplotpart = .true.
  irerender = .false.
  interactivereplot = .false.
  nyplotfirstonpage = 1 ! should be unnecessary, but to be on the safe side
  ifirststeponpage = 1  ! again, should be unnecessary
  nframefirstonpage = 1

  iplotxtemp(:)   = 1  ! this is just to be safe, so any spurious array access
  iplotytemp(:)   = 2  ! does not give an out-of-bounds error
  irendertemp(:)  = 0
  ivecplottemp(:) = 0

  xmin = 0.
  xmax = 0.
  ymin = 0.
  ymax = 0.
  rendermin = 0.
  rendermax = 0.
  contmin   = 0.
  contmax   = 0.
  vecmax    = 0.
  xminadapt = huge(xminadapt)
  xmaxadapt = -huge(xmaxadapt)
  renderminadapt = huge(renderminadapt)
  rendermaxadapt = -huge(rendermaxadapt)
  contminadapt   = huge(contminadapt)
  contmaxadapt   = -huge(contmaxadapt)

  !xminpagemargin = renvironment('SPLASH_MARGIN_XMIN')
  !xmaxpagemargin = renvironment('SPLASH_MARGIN_XMAX')
  !yminpagemargin = renvironment('SPLASH_MARGIN_YMIN')
  !ymaxpagemargin = renvironment('SPLASH_MARGIN_YMAX')

  if (ndim.eq.1) x_sec = .false. ! can't have xsec in 1D
  nxsec = 1

  iamrendering = .false.
  idoingvecplot = .false.
  if (ipicky.eq.numplot+1) then   ! multiplot
     imulti = .true.
     nyplots = nyplotmulti
     iplotx = 0  ! set these to zero by default for multiplots
     iploty = 0  ! (they should not be used in that case)
     !
     !--if doing multiplot can only take a single cross section slice
     !
     flythru = .false.
     nxsec = 1
     !
     !--work out whether to tile plots and make labelling decisions
     !
     if (any(multiplotx(1:nyplotmulti).ne.multiplotx(1))) isamexaxis = .false.
     if (any(multiploty(1:nyplotmulti).ne.multiploty(1))) then
        isameyaxis = .false.
        if (any(multiploty(1:nyplotmulti).eq.icolpixmap)) then
           isameyaxis = .true.
           do i=1,nyplotmulti
              if (.not.(multiploty(i).eq.icolpixmap .or. multiploty(i).eq.multiploty(1))) isameyaxis = .false.
           enddo
        endif
     endif
     if (any(irendermulti(1:nyplotmulti).gt.0)) iamrendering = .true.
     if (any(x_secmulti(1:nyplotmulti))) x_sec = .true.
     if (any(ivecplotmulti(1:nyplotmulti).gt.0)) idoingvecplot = .true.
  else
     !
     !--or else set number of plots = 1 and use ipicky and ipickx
     !
     imulti = .false.
     nyplots = 1
     iploty = ipicky
     iplotx = ipickx
     if (irender_nomulti.gt.0) iamrendering = .true.
     if (ivecplot.gt.0) idoingvecplot = .true.
  endif

  !------------------------------------------------------------------------
  ! initialise options to be set before plotting

  !
  !--determine z coordinate for 3D plots
  !
  icoordplot = is_coord(iploty,ndim) .and. is_coord(iplotx,ndim)
  iallrendered = iamrendering
  iplotz = 0
  if (imulti) then
     do i=1,nyplotmulti
        if (is_coord(multiplotx(i),ndim) .and. is_coord(multiploty(i),ndim)) then
           icoordplot = .true.
           !--this check is to see if any co-ordinate plots involve just particles
           !  (if so need to initialise the cross section slice width)
           if (irendermulti(i).le.0) iallrendered = .false.
           iplotzprev = iplotz
           !!--work out coordinate that is not being plotted on cross-section/ 3D plots
           iplotz = 0
           if (ndim.ge.3 .and. (x_sec .or. use3Dperspective .or. irotate)) then
              do j=1,ndim
                 if ((multiplotx(i).ne.multiploty(i)).and. &
                     (ix(j).ne.multiplotx(i)).and.(ix(j).ne.multiploty(i))) iplotz = ix(j)
              enddo
              !--use only first iplotz in the case of multiple slices
              !  (only effect is on default values for slice thickness etc below)
              if (iplotzprev.gt.0) iplotz = iplotzprev
           endif
        endif
     enddo
  elseif (icoordplot) then
     !!--work out coordinate that is not being plotted
     if (ndim.ge.3) then
        do j=1,ndim
           if ((iplotx.ne.iploty).and. &
               (ix(j).ne.iplotx).and.(ix(j).ne.iploty)) iplotz = ix(j)
        enddo
     endif
  endif
  if (debugmode) print*,'DEBUG: iplotz = ',iplotz

  !
  !--work out whether or not to tile plots on the page
  !  if plots are coord plots, make tiling decisions based on iadaptcoords
  !  otherwise use iadapt
  !
  if (icoordplot) then
     iadapting = iadaptcoords
  else
     iadapting = iadapt
  endif
  tile_plots = tile .and. (isamexaxis.and.isameyaxis .or. isameyaxis.and.ndown.eq.1  &
                      .or. isamexaxis.and.nacross.eq.1) .and. (nacross*ndown.gt.1)
  !--do not tile if limits are adaptive
  if (tile_plots .and. (iadapting .or. (iamrendering .and. iadapt .and. iColourbarStyle.gt.0))) then
     print "(a)",'WARNING: cannot tile plots because limits are set to adaptive'
     tile_plots = .false.
  endif

  !--( a further constraint on plot tiling is required in the case of
  !    multiple renderings which would involve different colour bars )
  OneColourBarPerRow = .false.
  if (iamrendering .and. icolours.ne.0 .and. iColourbarStyle.gt.0) then
     !--this option means that a margin is set aside for a colour bar on tiled plots
     iAllowspaceforcolourbar = .true.
     !--do not allow tiled plots if multiple (different) colour bars are plotted
     if (tile_plots) then
        ifirst = 0
        do i=1,nyplots
           if (irendermulti(i).gt.0 .and. ifirst.eq.0) ifirst = i
           if (ifirst.gt.0) then
              if (irendermulti(i).gt.0 .and. irendermulti(i).ne.irendermulti(ifirst)) then
                 tile_plots = .false.
              endif
           endif
        enddo
        !--this means colour bars are not the same, but we can still tile if
        !  all the colour bars in each row are the same
        if (.not.tile_plots .and. mod(nacross*ndown,nyplots).eq.0) then
           OneColourBarPerRow = .true.
           allocate(ifirstinrow(ndown))
           ifirstinrow(:) = 0
           do i=1,nyplots
              if (usecolumnorder) then
                 irow = (i-1)/nacross + 1
              else
                 irow = i - ((i-1)/ndown)*ndown
              endif
              if (ifirstinrow(irow).eq.0) ifirstinrow(irow) = i
              if (irendermulti(i).ne.irendermulti(ifirstinrow(irow))) then
                 OneColourBarPerRow = .false.
              endif
           enddo
           deallocate(ifirstinrow)
           if (OneColourBarPerRow) tile_plots = .true.
        endif
        if (.not.tile_plots) print "(a)",'WARNING: cannot tile plots because of multiple colour bars'
     endif
  else
     iAllowspaceforcolourbar = .false.
  endif

  if (icoordplot) then
     if (x_sec .and. iplotz.gt.0) then
!
!--if series of cross sections (flythru), set position of first one
!
        if (flythru) then
           print 32,label(iplotz)
32         format('enter number of ',a1,' cross-section slices')
           read*,nxsec
           !!--dz is the distance between slices
           dz = (lim(iplotz,2)-lim(iplotz,1))/float(nxsec)
           zslicepos = lim(iplotz,1) - 0.5*dz
           xsecpos_nomulti = zslicepos
        else
!
!--if single cross-section, read position of cross-section slice
!
           if (.not.imulti) then
              !--make sure position falls within the limits
              if (xsecpos_nomulti.lt.lim(iplotz,1) &
              .or.xsecpos_nomulti.gt.lim(iplotz,2)) then
                  xsecpos_nomulti = (lim(iplotz,2)+lim(iplotz,1))/2.
              endif
              call prompt(' enter '//trim(label(iplotz))// &
                       ' position for cross-section slice:', &
                       xsecpos_nomulti,lim(iplotz,1),lim(iplotz,2))
           endif
!
!--set thickness if plotting particles
!  (default thickness is half of the average particle spacing)
!
           if (.not.iallrendered .or. icolour_particles) then
              if (xsecwidth.gt.0. .and. xsecwidth.lt.(lim(iplotz,2)-lim(iplotz,1))) then
                 !--already set
                 dzsuggest = xsecwidth
              else
                 !--xsecwidth not set; suggest a good value
                 npartdim = int(maxval(npartoftype(:,1))**(1./real(ndim)))
                 print*,'average # of particles in each dimension = ',npartdim
                 if (npartdim.gt.0) then
                    dzsuggest = (lim(iplotz,2)-lim(iplotz,1))/float(npartdim)
                 else
                    dzsuggest = 0.01*(lim(iplotz,2)-lim(iplotz,1))
                 endif
              endif
              dz = dzsuggest

              if (imulti) then
                 call prompt(' enter thickness for cross section slice(s):', &
                           dz,0.0,lim(iplotz,2)-lim(iplotz,1))
              else
                 call prompt(' enter thickness of cross section slice:', &
                           dz,0.0,lim(iplotz,2)-lim(iplotz,1))
              endif
              !--if dz has been set from the prompt, save the setting,
              !  otherwise suggest (possibly different) value again next time
              if (abs(dz-dzsuggest).gt.tiny(dz)) xsecwidth = dz

           elseif (ndim.eq.3) then
!
!--for rendered cross sections in 3D, set thickness to 10%
!  this is the distance slices are moved up and down in interactive mode
!
              dz = 0.1*(lim(iplotz,2)-lim(iplotz,1))
           endif
        endif ! flythru or single
!
!--set up for 1D cross sections through 2D data
!
    elseif (ndim.eq.2 .and. x_sec) then
       ians = .false.
       call prompt('set cross section position interactively?',ians)

       if (ians) then
       !
       !--set cross section position interactively
       !
          call plot_init('/xw',ierr)
          call plot_env(lim(1,1),lim(1,2),lim(2,1),lim(2,2),1,0)
          ierr = plot_curs(xseclineX1,xseclineY1,char)
          print*,'please select cross section line'
          ierr = plot_band(1,1,xseclineX1,xseclineY1,xseclineX2,xseclineY2,char)
          print*,'cross section line: xmin = ',xseclineX1,' xmax = ',xseclineX2
          print*,'                    ymin = ',xseclineY1,' ymax = ',xseclineY2
          call plot_close
       else
       !
       !--set position manually
       !
          if (abs(xseclineX2-xseclineX1).lt.1.e-5 .and. &
              abs(xseclineY2-xseclineY1).lt.1.e-5) then
          !--if not already set (ie. if all = 0.0)
          !  then set default line to diagonal across the domain
             xseclineX1 = lim(1,1)
             xseclineX2 = lim(1,2)
             xseclineY1 = lim(2,1)
             xseclineY2 = lim(2,2)
          endif
          print*,'enter position of cross section through 2D data:'
          call prompt('enter xmin of cross section line',xseclineX1)
          call prompt('enter xmax of cross section line',xseclineX2)
          call prompt('enter ymin of cross section line',xseclineY1)
          call prompt('enter ymax of cross section line',xseclineY2)
       endif
     endif

     if (iplotz.gt.0 .and. use3Dperspective) then
!
!--initialise 3D perspective
!
       !--set default values if none set
       if (abs(zobserver).lt.tiny(zobserver)) zobserver = 10.*lim(iplotz,2)
       if (abs(dzscreenfromobserver).lt.tiny(dzscreenfromobserver)) dzscreenfromobserver = zobserver
       call prompt('enter z coordinate of observer ',zobserver)
       dzscreenfromobserver = zobserver
!       call prompt('enter distance for unit magnification ',dzscreenfromobserver,0.)
!
!--initialise opacity for 3D opacity rendering
!
       if (use3Dopacityrendering .and. (iamrendering .or. idoingvecplot)) then
          hav = lim(ih,1) !! 0.5*(lim(ih,2) + lim(ih,1))
          if (hav.le.epsilon(hav)) hav = 0.5*lim(ih,2) ! take 0.5*max if min is zero
          if (ipmass.gt.0) then
             pmassav = lim(ipmass,1)
             if (pmassav.le.epsilon(hav)) pmassav = 0.5*lim(ipmass,2) ! take 0.5*max if min is zero
          else  ! handle case where mass is not a data column
             pmassav = maxval(masstype)
             do i=1,ntypes
                if (iplotpartoftype(i) .and. usetypeinrenderings(i) &
                    .and. any(masstype(i,:).gt.0.)) pmassav = min(pmassav,maxval(masstype(i,:)))
             enddo
          endif
          print*,'using current h and pmass limits to calculate kappa (cross section/unit mass)'
          print*,'min h = ',hav,' min particle mass = ',pmassav
          print*,'[ kappa = pi*h_min**2/(particle_mass*n_smoothing_lengths) ]'
          call prompt('enter approximate surface depth (number of smoothing lengths):',taupartdepth,0.)
          rkappafac = pi*hav*hav/(pmassav*coltable(0))
          print*,'kappa (particle cross section per unit mass) = ',rkappafac/taupartdepth
       endif
    endif

  endif

  !!--prompt for options if plotting power spectrum
  if (iploty.eq.ipowerspec .and. .not. imulti &
     .or. (imulti.and.any(multiploty(1:nyplotmulti).eq.ipowerspec))) then
     call options_powerspec
  endif

  !!--prompt for options if plotting PDFs
  if (iploty.eq.ipdf .and. .not. imulti &
     .or. (imulti.and.any(multiploty(1:nyplotmulti).eq.ipdf))) then
     call options_pdf
  endif

  !!--for fast data read, set which columns are required from the file
  !   (note that required(0)= whatever is a valid statement, just has no effect)
  required = .false.  ! by default, no columns required
  if (debugmode) print*,'DEBUG: imulti = ',imulti,' iamrendering = ',iamrendering
!  if (fastdataread) then
     if (imulti) then
        required(multiplotx(1:nyplotmulti)) = .true.
        required(multiploty(1:nyplotmulti)) = .true.
        required(irendermulti(1:nyplotmulti)) = .true.
        required(icontourmulti(1:nyplotmulti)) = .true.
     else
        if (iploty.ne.icolpixmap) required(iplotx) = .true.
        required(iploty) = .true.
     endif
     required(iplotz) = .true.
     if ((iamrendering .or. idoingvecplot) .and. &
        (iploty.ne.icolpixmap .or. imulti .or. iploty.eq.0)) then
        required(ipmass) = .true.
        required(irho) = .true.
        required(ih) = .true.
        required(irender_nomulti) = .true.
        required(icontour_nomulti) = .true.
     endif
  !!--need to read columns used for range restrictions
     do i=1,ndataplots
        if (rangeset(i)) required(i) = .true.
     enddo
  !!--need mass for some exact solutions
     if (iexact.eq.7 .or. iploty.eq.isurfdens) required(ipmass) = .true.
     if (iploty.eq.itoomre .and. iploty.gt.0) required(iutherm) = .true.
  !!--only require actual dependencies of calculated quantities
     if (any(required(ncolumns+1:ncolumns+ncalc))) call get_calc_data_dependencies(required)
     !if (any(required(ncolumns+1:ncolumns+ncalc))) required = .true.
  !!--vectors
     if (imulti) then
        do i=1,nyplotmulti
           if (ivecplotmulti(i).gt.0) then
              required(ivecplotmulti(i):ivecplotmulti(i)+ndimV-1) = .true.
           endif
        enddo
     elseif (ivecplot.gt.0) then
        required(ivecplot:ivecplot+ndimV-1) = .true.
     endif
   !!--if geometry is not default must read all coords
   !!  and if we are plotting a vector component, all components
     if (icoordsnew.ne.icoords) then
        required(ix(1:ndim)) = .true.
        if (iplotx.gt.0 .and. iplotx.le.numplot) then
           if (iamvec(iplotx).gt.0) required(iamvec(iplotx):iamvec(iplotx)+ndimV-1) = .true.
        endif
        if (iploty.gt.0 .and. iploty.le.numplot) then
           if (iamvec(iploty).gt.0) required(iamvec(iploty):iamvec(iploty)+ndimV-1) = .true.
        endif
     endif
!  endif
  if (debugmode) print*,'DEBUG: required(1:ncolumns) = ',required(1:ncolumns+ncalc)

  !!--read step titles (don't need to store ntitles for this)
  nsteplegendlines = 0
  call read_steplegend(nsteplegendlines)
  !!--read plot titles
  ntitles = 0
  call read_titles(ntitles)

  !!------------------------------------------------------------------------
  !! initialise the plotting library

  if (len_trim(device).le.0) then
     devstring = '?'           ! prompt for device
  else
     devstring = trim(device)  ! device specified on command line
  endif

  ierr = 1
  ntries = 0
  do while(ierr.ne.0)
     ntries = ntries + 1
     if (ipapersize.gt.0 .and. papersizex.gt.0.0 .and. aspectratio.gt.0.0 ) then
        call plot_init(trim(devstring),ierr,papersizex,aspectratio,ipapersizeunits)
     else
        call plot_init(trim(devstring),ierr)  ! use default paper size
     endif
     !--abort if device specified on command line returns an error
     if (len_trim(device).gt.0 .and. ierr.ne.0) then  ! -ve indicates an error
        print "(a)",' ERROR: unknown device "'//trim(device)//'"'
        stop
     endif
     if (ierr.ne.0) print "(a)",' ERROR opening plotting device'
     if (ntries.gt.10) stop
  enddo

  !--query whether or not device is interactive
  if (plot_qcur()) then
     !--turn menu and interactive mode on if
     !  interactive device invoked from the command line
     if (nomenu) then
        interactive = .true.
        nomenu      = .false.
     endif
     !--smoothing length is required if interactive device and coordinate plot
     if (icoordplot) required(ih) = .true.
  endif

  !!--set background/foreground colours
  call set_pagecolours(iPageColours)

  !!--set colour table
  !   do this regardless of whether rendering or not
  if (iamrendering .and. icolours.ne.0) then
     call colour_set(icolours)
     ihavesetcolours = .true.
  else
     ihavesetcolours = .false.
  endif

  !!--determine whether or not the device is vector or not
  !   this affects the choice of line with (if auto line width is used -- see below)
  !   and also the automatic resolution determination (should not apply to vector devices)
  !
  call plot_qinf('TYPE',devstring,ilen)
  select case(devstring(1:ilen))
  case('PS','CPS','VPS','VCPS','NULL','LATEX')
     vectordevice = .true.
  case default
     vectordevice = .false.
  end select

  !!--set line width (0=auto based on whether device is vector or not)
  if (linewidth.le.0) then
     if (vectordevice) then
        print "(a)",' setting line width = 2 for '//devstring(1:ilen)//' device'
        call plot_slw(2)
     else
        call plot_slw(1)
     endif
  else
     call plot_slw(linewidth)
  endif

  linecolourthisstep = linecolour
  linestylethisstep = linestyle

end subroutine initialise_plotting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plotstep(ipos,istep,istepsonpage,irender_nomulti,icontour_nomulti,ivecplot, &
                    iamtype,npartoftype,masstype,dat,timei,gammai,ipagechange,iadvance)
  use params,             only:int1,maxparttypes,doub_prec
  use colours,            only:colour_set
  use filenames,          only:nsteps,rootname,ifileopen,tagline
  use exact,              only:exact_solution,atstar,ctstar,sigma,iPlotExactOnlyOnPanel
  use toystar1D,          only:exact_toystar_ACplane
  use toystar2D,          only:exact_toystar_ACplane2D
  use labels,             only:label,shortlabel,labelvec,iamvec,lenlabel,lenunitslabel,ih,irho,ipmass,ix,iacplane, &
                               ipowerspec,isurfdens,itoomre,ispsound,iutherm,ipdf,icolpixmap,is_coord,labeltype,&
                               labelzintegration,unitslabel,integrate_label,get_sink_type
  use limits,             only:lim,get_particle_subset,lim2,lim2set
  use multiplot,          only:multiplotx,multiploty,irendermulti,ivecplotmulti,itrans, &
                               icontourmulti,x_secmulti,xsecposmulti,iusealltypesmulti,iplotpartoftypemulti
  use particle_data,      only:maxpart,maxcol,icolourme
  use settings_data,      only:numplot,ndataplots,icoords,icoordsnew,ndim,ndimV,nfreq,iRescale, &
                               iendatstep,ntypes,UseTypeInRenderings,itracktype,itrackoffset,&
                               required,ipartialread,xorigin,lowmemorymode,debugmode,iverbose
  use settings_limits,    only:iadapt
  use settings_part,      only:iexact,iplotpartoftype,imarktype,PlotOnRenderings,UseTypeInContours, &
                               iplotline,linecolourthisstep,linestylethisstep,ifastparticleplot, &
                               iploterrbars,ilocerrbars
  use settings_page,      only:nacross,ndown,interactive,iaxis,usesquarexy,yscalealt,labelyalt, &
                               charheight,iPlotTitles,vpostitle,hpostitle,fjusttitle,nstepsperpage
  use settings_render,    only:npix,ncontours,icolours,iColourBarStyle,icolour_particles,&
                               inormalise_interpolations,ifastrender,ilabelcont,double_rendering,&
                               projlabelformat,iapplyprojformat
  use settings_vecplot,   only:npixvec,iplotpartvec
  use settings_xsecrot,   only:nxsec,irotateaxes,xsec_nomulti,irotate,flythru,use3Dperspective, &
                               use3Dopacityrendering,writeppm,anglex,angley,anglez,zobserver,&
                               dzscreenfromobserver,taupartdepth,xsecpos_nomulti, &
                               xseclineX1,xseclineX2,xseclineY1,xseclineY2, &
                               nseq,nframes,getsequencepos,insidesequence,rendersinks
  use settings_powerspec, only:nfreqspec,freqmin,freqmax,ipowerspecx,ipowerspecy,&
                               idisordered,npdfbins
  use settings_units,     only:units,unitzintegration
!
!--subroutines called from this routine
!
  use colourparts
  use transforms,            only:transform,transform_limits,transform_label,transform_inverse,islogged
  use interactive_routines
  use part_utils,            only:get_tracked_particle,locate_first_two_of_type,get_binary
  use particleplots,         only:particleplot,plot_errorbarsx,plot_errorbarsy
  use powerspectrums,        only:powerspectrum,powerspec3D_sph
  use interpolations1D,      only:interpolate1D
  use interpolations2D,      only:interpolate2D, interpolate2D_xsec !, interpolate2D_pixels
  use interpolations3D,      only:interpolate3D
  use projections3D,         only:interpolate3D_projection
  use projections3Dgeom,     only:interpolate3D_proj_geom
  use interpolate3D_opacity, only:interp3D_proj_opacity !,interp3D_proj_opacity_writeppm
  use xsections3D,           only:interpolate3D_fastxsec,interpolate3D_xsec_vec
  use render,                only:render_pix
  use pagesetup,             only:redraw_axes
  use disc,                  only:disccalc,discplot
  use exactfromfile,         only:exact_fromfile
  use exact,                 only:iexact_rochelobe,use_sink_data,mprim,msec,xprim,xsec
  use write_pixmap,          only:iwritepixmap,writepixmap,write_pixmap_ppm,readpixmap
  use pdfs,                  only:pdf_calc,pdf_write
  use plotutils,             only:plotline
  use geometry,              only:coord_is_length
  use geomutils,             only:changecoords,changeveccoords
  use legends,               only:ipanelselect
  use plotlib,               only:plot_sci,plot_page,plot_sch,plot_qci,plot_qls,plot_sls, &
                                  plot_line,plot_pt1,plotlib_is_pgplot,plotlib_supports_alpha
  implicit none
  integer, intent(inout) :: ipos, istepsonpage
  integer, intent(in)    :: istep,irender_nomulti,icontour_nomulti,ivecplot
  integer(kind=int1), dimension(:), intent(in) :: iamtype
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real,    dimension(maxparttypes), intent(in) :: masstype
  real,    dimension(:,:),          intent(in) :: dat
  real,       intent(in) :: timei,gammai
  logical,    intent(in) :: ipagechange
  integer, intent(inout) :: iadvance

  logical, dimension(maxparttypes) :: iusetype

  integer :: ntoti,iz,iseqpos,itrackpart
  integer :: i,j,k,icolumn,irow
  integer :: nyplot,nframesloop
  integer :: irenderpart
  integer :: npixyvec,nfreqpts,ipixxsec
  integer :: icolourprev,linestyleprev
  integer :: ierr,ipt,nplots,nyplotstart,iaxisy,iaxistemp,icol
  integer :: ivectemp,iamvecx,iamvecy,itransx,itransy,itemp
  integer :: iframe,isize,isinktype,isink1,isink2

  real, parameter :: tol = 1.e-10 ! used to compare real numbers
  real, parameter :: error_in_log = -666. ! magic number used to flag error with log(0.)
  real(doub_prec) :: unit_mass, unit_r, unit_u
  real, dimension(:), allocatable    :: xplot,yplot,zplot
  real, dimension(:), allocatable    :: hh,weight
  real, dimension(:), allocatable    :: renderplot
  real, dimension(:,:), allocatable  :: vecplot
  real :: rkappa
  real :: zslicemin,zslicemax,dummy,pmassmin,pmassmax,pmassav(1)
  real :: pixwidth,pixwidthy,pixwidthvec,pixwidthvecy,dxfreq

  character(len=lenlabel+lenunitslabel) :: labelx,labely,labelz,labelrender,labelvecplot,labelcont
  character(len=lenunitslabel) :: labeltimeunits
  character(len=12) :: string

  logical :: iPlotColourBar, rendering, inormalise, logged, loggedcont
  logical :: dumxsec, isetrenderlimits, iscoordplot
  logical :: ichangesize, initx, inity, isameweights, volweightedpdf, got_h
  logical, parameter :: isperiodicx = .false. ! feature not implemented
  logical, parameter :: isperiodicy = .false. ! feature not implemented
  logical, parameter :: isperiodicz = .false. ! feature not implemented
  logical, dimension(maxparttypes) :: PlotonRender_tmp

34   format (25(' -'))

  !--set labels to blank (just in case)
  labelx = ' '
  labely = ' '
  labelz = ' '
  labelrender = ' '
  labelvecplot = ' '

  !
  !--allocate temporary memory required for plotting
  !
  isize = max(maxpart,2000)
  !--do not allocate the temporary arrays if the dat array has not been allocated
  if (lowmemorymode .and. maxcol.eq.0) then
     isize = 2000
  endif
  if (debugmode) print*,'DEBUG: in plotstep, allocating local memory...'
  ierr = 0
  allocate(xplot(isize),stat=ierr)
  if (ierr /= 0) stop 'out of memory in plotstep allocating temporary x array'
  allocate(yplot(isize),stat=ierr)
  if (ierr /= 0) stop 'out of memory in plotstep allocating temporary y array'
  allocate(zplot(isize),stat=ierr)
  if (ierr /= 0) stop 'out of memory in plotstep allocating temporary z array'
  if (allocated(xplot)) xplot = 0.
  if (allocated(yplot)) yplot = 0.
  if (allocated(zplot)) zplot = 0.

  allocate(hh(maxpart),weight(maxpart),stat=ierr)
  if (ierr /= 0) stop 'out of memory in plotstep allocating temporary h,weight arrays'
  hh = 0.
  if (debugmode) print*,'DEBUG: in plotstep, allocated local memory successfully'

  dummy = 0.
  labeltimeunits = ' '
  dumxsec = .false.
  isetrenderlimits = .false.
  k = nxsec ! matters for lastplot in page_setup for non-coord plots
  if (iReScale) labeltimeunits = unitslabel(0)
  iaxistemp = iaxis

  !--set the arrays needed for rendering if they are present
  if (ih.gt.0 .and. ih.le.ndataplots .and. (required(ih) .or. .not.ipartialread)) then
     hh(:) = dat(:,ih)
     got_h = .true.
  else
     got_h = .false.
  endif

  if (ipmass.gt.0 .and. ipmass.le.ndataplots) then
     if (required(ipmass)) then
        pmassmin = minval(dat(:,ipmass))
        pmassmax = maxval(dat(:,ipmass))
     else
        pmassmin = 0.
        pmassmax = 0.
     endif
  else
     pmassmin = minval(masstype,mask=(masstype.gt.0.))
     pmassmax = maxval(masstype)
     pmassav = masstype(1)
     if (pmassav(1).le.0.) then
        do i=ntypes,2,-1
           if (UseTypeInRenderings(i) .and. iplotpartoftype(i) .and. masstype(i).gt.0.) then
              pmassav = masstype(i)
           endif
        enddo
     endif
  endif
  !
  !--set number of particles to use in the interpolation routines
  !  (by default, only the gas particles)
  !
  ntoti = sum(npartoftype)
  ninterp = npartoftype(1)
  if (any(UseTypeInRenderings(2:ntypes).and.iplotpartoftype(2:ntypes)) &
      .or. size(iamtype).gt.1 .or. (use3Dopacityrendering .and. rendersinks)) ninterp = ntoti

  !--work out the identity of a particle being tracked
  if (debugmode) print*,'DEBUG: itracktype = ',itracktype,' itrackoffset = ',itrackoffset
  itrackpart = get_tracked_particle(itracktype,itrackoffset,npartoftype,iamtype)
  if (itrackpart.eq.0) then
     write(string,"(i12)") itrackoffset
     string = adjustl(string)
     if (itracktype.gt.0 .and. itracktype.le.ntypes) then
        print "(/,a,/)",' WARNING: tracked '//trim(labeltype(itracktype))//' particle #'//trim(string)//' not found in data'
     elseif (itrackoffset.gt.0) then
        print "(/,a,/)",' WARNING: tracked particle #'//trim(string)//' not found in data'
     endif
  else
     write(string,"(i12)") itrackpart
     print "(/,a,/)",' Tracking particle #'//trim(adjustl(string))
  endif

  !--non-SPH particle types cannot be used in contours
  where (.not.UseTypeInRenderings(:))
     UseTypeInContours(:) = .false.
  end where

  !--check for consistency that if particles are not plotted,
  !  they are also not plotted on renderings
  do i=1,ntypes
     if (.not.iplotpartoftype(i)) PlotOnRenderings(i) = .false.
  enddo
  !
  !--check whether or not the particle types used for contouring are
  !  the same as the particle types used for rendering
  !
  isameweights = .true.
  do i=1,ntypes
     if (UseTypeInRenderings(i) .and. &
        .not.(iplotpartoftype(i).eqv.UseTypeInContours(i))) isameweights = .false.
  enddo
  !
  !--set weight factor for interpolation routines
  !
  ihavesetweights = .false.
  if (iamrendering .or. idoingvecplot) then
     if (debugmode) print*,'DEBUG: setting interpolation weights...'
     call set_weights(weight,dat,iamtype,(iplotpartoftype .and. UseTypeInRenderings))
  else
     if (debugmode) print*,'DEBUG: interpolation weights not set because no rendering...'
  endif

  !--set the colour table if it has not been set and particles have been coloured previously
  if (any(icolourme(1:ntoti).gt.16) .and. .not.ihavesetcolours) call colour_set(icolours)

  !
  !--exclude subset of particles if parameter range restrictions are set
  !
  call get_particle_subset(icolourme,dat,ndataplots)
  !
  !--add a loop over frames for animation sequences
  !  but only generate extra frames if we are inside a sequence
  !
  iseqpos = (ipos-1)/(nacross*ndown) + 1
  !print*,' iseqpos = ',iseqpos,ipos

  iframe = 0
  if (nseq.gt.0 .and. insidesequence(iseqpos)) then
     if (nacross*ndown.eq.1) then
        nframesloop = nframes
     else
        nframesloop = max(iframesave+1,1)
        iframe = iframesave
        !print*,'iframe=',iframesave,'iseqpos = ',iseqpos,insidesequence(iseqpos)
     endif
  else
     nframesloop = 1
  endif
  if (debugmode) print*,'DEBUG: starting frame loop...'

  !--loop over frames: flexible to allow forwards/backwards in interactive mode
  over_frames: do while (iframe.lt.nframesloop)

  if (interactivereplot .and. ipos.eq.ifirststeponpage .and. iframe.eq.0) then
     iframe = min(nframefirstonpage,nframesloop)
  else
     iframe = iframe + 1
  endif

  !print*,'iframe = ',iframe, ipagechange,nstepsperpage
  !--sanity check on frame number, should never happen...
  if (iframe.eq.0) then
     print*,' Internal error in iframe, setting to 1 '
     iframe = 1
  endif
  !-------------------------------------
  ! loop over plots per timestep
  ! (jump to first on the page if replotting in interactivemode)
  !-------------------------------------
  if (interactivereplot .and. ipos.eq.ifirststeponpage .and. iframe.eq.nframefirstonpage) then
     nyplotstart = nyplotfirstonpage
     ipanel = 0
  else
     nyplotstart = 1
  endif

  over_plots: do nyplot=nyplotstart,nyplots

     if (nyplot.gt.1 .or. iframe.gt.1) print 34
     !--make sure character height is set correctly
     call plot_sch(charheight) ! in PGPLOT scaled units

     iPlotColourBar = .false.   ! should be false by default until set to true
     iaxistemp = iaxis

     !--set current x, y, render and vector plot from multiplot array
     if (imulti) then
        iploty = multiploty(nyplot)
        iplotx = multiplotx(nyplot)
        irender = irendermulti(nyplot)
        ivectorplot = ivecplotmulti(nyplot)
        icontourplot = icontourmulti(nyplot)
        iplotcont = .false. !iplotcontmulti(nyplot)
        x_sec = x_secmulti(nyplot)
        zslicepos = xsecposmulti(nyplot)
        if (iusealltypesmulti(nyplot)) then
           iusetype(:) = iplotpartoftype(:)
        else
           iusetype(:) = iplotpartoftypemulti(:,nyplot)
        endif
        if (irender.gt.0 .or. icontourplot.gt.0 .or. ivectorplot.gt.0) then
           if (debugmode) print*,'DEBUG: resetting interpolation weights for multiplot...'
           call set_weights(weight,dat,iamtype,(iusetype .and. UseTypeInRenderings))
        endif
     else
        if (.not.interactivereplot) irender = irender_nomulti
        ivectorplot = ivecplot
        icontourplot = icontour_nomulti
        iplotcont = .false. !iplotcont_nomulti
        if (.not.interactivereplot) x_sec = xsec_nomulti
        if (.not.interactivereplot .and. x_sec) zslicepos = xsecpos_nomulti
        iusetype(:) = iplotpartoftype(:)
     endif

     !--if the contour plot is the same as the rendered plot,
     !  do not interpolate twice. Instead simply plot the contours
     !  of the rendered quantity when plotting the render plot.
     if (irender.gt.0 .and. irender.le.numplot) then
        if (icontourplot.eq.irender .and. isameweights) then
           icontourplot = 0
           iplotcont = .true.
           !print "(a)",' contouring same as rendering'
        elseif (icolours.eq.0) then
           iplotcont = .true.
        endif
     else
        !--contours not allowed if not rendering
        !  (because this can be achieved by rendering with colour scheme 0)
        icontourplot = 0
     endif
     !--flag to indicate that we have actually got the contoured quantity,
     !  set to true once interpolation has been done.
     if (.not.interactivereplot .or. irerender) gotcontours = .false.

     if (icolour_particles) then
        irenderpart = irender
        irenderplot = 0
     else
        irenderpart = 0
        irenderplot = irender
     endif

     if (ivectorplot.gt.0) iplotpart = iplotpartvec

     !--if replotting in interactive mode, use the temporarily stored plot limits
     !  (check iplot values are sensible though, otherwise will seg fault here)
     if (interactivereplot .and. (nacross*ndown.gt.1 .or. iploty.gt.ndataplots)) then
        if (iplotx.gt.0 .and. iplotx.le.numplot) then
           xmin = xminmulti(iplotx)
           xmax = xmaxmulti(iplotx)
        endif
        if (iploty.gt.0 .and. iploty.le.numplot) then
           ymin = xminmulti(iploty)
           ymax = xmaxmulti(iploty)
        endif
        if (irender.gt.0 .and. irender.le.numplot) then
           rendermin = xminmulti(irender)
           rendermax = xmaxmulti(irender)
           if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
              contmin = xminmulti(icontourplot)
              contmax = xmaxmulti(icontourplot)
           endif
        endif
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! initialisation for plots of particle data
     ! copy from main dat array into xplot, yplot
     ! also set labels and plot limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     initx = (iplotx.gt.0 .and. iplotx.le.ndataplots)
     inity = (iploty.gt.0 .and. iploty.le.ndataplots)
     if (iplotx.gt.0 .and. iplotx.le.numplot) labelx = label(iplotx)
     if (iploty.gt.0 .and. iploty.le.numplot) labely = label(iploty)
     iscoordplot = (is_coord(iplotx,ndim) .and. is_coord(iploty,ndim))

     initdataplots: if (initx .or. inity) then
        if (debugmode) print*,'DEBUG: initialising data plots...',initx,inity,iplotx,iploty,ntoti,size(xplot)
        if (initx) then
           !--check for errors
           if (iplotx.gt.size(dat(1,:)) .or. iplotx.lt.1) then
              print*,'ERROR: Internal error with out-of-bounds x column = ',iplotx
              exit over_plots
           endif
           xplot(1:ntoti) = dat(1:ntoti,iplotx)
           iamvecx = iamvec(iplotx)
        else
           iamvecx = 0
        endif
        if (inity) then
           !--check for errors
           if (iploty.gt.size(dat(1,:)) .or. iploty.lt.1) then
              print*,'ERROR: Internal error with out-of-bounds y column = ',iploty
              exit over_plots
           endif
           yplot(1:ntoti) = dat(1:ntoti,iploty)
           iamvecy = iamvec(iploty)
        else
           iamvecy = 0
        endif
        if (debugmode) print*,'DEBUG: setting itrans...'
        itransx = 0
        itransy = 0
        if (iplotx.gt.0 .and. iplotx.le.numplot) itransx = itrans(iplotx)
        if (iploty.gt.0 .and. iploty.le.numplot) itransy = itrans(iploty)

        zplot = 0.   !--set later if x-sec
        zslicemin = -huge(zslicemax) !-- " "
        zslicemax = huge(zslicemax)
        if (.not.interactivereplot) then
!           if (iplotx.gt.0 .and. iplotx.le.numplot .and. ipos.eq.ifirststeponpage) then
           if (iplotx.gt.0 .and. iplotx.le.numplot) then
              xmin = lim(iplotx,1)
              xmax = lim(iplotx,2)
           endif
!           if (iploty.gt.0 .and. iploty.le.numplot .and. ipos.eq.ifirststeponpage) then
           if (iploty.gt.0 .and. iploty.le.numplot) then
              ymin = lim(iploty,1)
              ymax = lim(iploty,2)
           endif
           angletempx = anglex
           angletempy = angley
           angletempz = anglez
           if (ndim.eq.3 .and. use3Dperspective) then
              dzscreentemp = dzscreenfromobserver
              zobservertemp = zobserver
              taupartdepthtemp = taupartdepth
           else
              dzscreentemp = 0.
              zobservertemp = 0.
              taupartdepthtemp = 0.
           endif
        else
           if (ndim.eq.3 .and. use3Dperspective) dzscreentemp = zobservertemp
        endif
        !
        !--flag for whether or not we have raw particle plot or not
        !  (not allowed to use transformations on coords otherwise)
        !
        rendering = (iscoordplot .and.(irenderplot.gt.0 .or. ivectorplot.gt.0) &
                     .and.(.not.icolour_particles))
        !
        !--change coordinate system if relevant
        !
        if (icoordsnew.ne.icoords) then
           !--do this if one is a coord but not if rendering
           call changecoords(iplotx,iploty,xplot,yplot,ntoti,ndim,itrackpart,dat)
           if (iamvecx.gt.0) call changeveccoords(iplotx,xplot,ntoti,ndim,itrackpart,dat)
           if (iamvecy.gt.0) call changeveccoords(iploty,yplot,ntoti,ndim,itrackpart,dat)
        endif

        !--apply transformations (log, 1/x etc) if appropriate
        !  also change labels and limits appropriately
        if (.not.(rendering)) then
           if (itransx.ne.0) call applytrans(xplot,xmin,xmax,labelx,itransx,'x',iplotx,iaxis,interactivereplot)
           if (itransy.ne.0) call applytrans(yplot,ymin,ymax,labely,itransy,'y',iploty,iaxis,interactivereplot)
        endif

        !--write username, date on plot
        !         if (nacross.le.2.and.ndown.le.2) call pgiden
        !
        !--adjust plot limits if adaptive plot limits set
        !  (find minimum/maximum only on particle types actually plotted)
        !
        if (itrackpart.le.0 .and. .not.(iscoordplot .and. irotate)) then
           if (initx) call adapt_limits(iplotx,xplot,xmin,xmax,xminadapti,xmaxadapti,'x',&
                                        iamtype,ntoti,npartoftype,iusetype,ipagechange)
           if (inity) call adapt_limits(iploty,yplot,ymin,ymax,yminadapti,ymaxadapti,'y',&
                                        iamtype,ntoti,npartoftype,iusetype,ipagechange)
        endif

        !!-reset co-ordinate plot limits if particle tracking
        if (itrackpart.gt.0 .and. .not.interactivereplot) then
           if (initx) call settrackinglimits(itrackpart,iplotx,xplot,xmin,xmax)
           if (inity) call settrackinglimits(itrackpart,iploty,yplot,ymin,ymax)
        endif

        !--override settings based on positions in sequence
        if (nseq.gt.0) then
           call getsequencepos(iseqpos,iframe,iplotx,iploty,irender, &
                angletempx,angletempy,angletempz,zobservertemp,dzscreentemp,taupartdepthtemp,&
                zslicepos,xmin,xmax,ymin,ymax,rendermin,rendermax,isetrenderlimits)
        endif
        !--for 3D perspective, do not plot particles behind the observer
        if (ndim.eq.3.and.use3Dperspective) then
           dzscreenfromobserver = zobserver
           zslicemax = zobservertemp
           if (use3Dopacityrendering) rkappa = rkappafac/taupartdepthtemp
        endif

     endif initdataplots
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! plots with co-ordinates as x and y axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (iscoordplot) then

        if (debugmode) print*,'DEBUG: starting coord plot...'
        if (.not.interactivereplot .or. irerender) then
           npixx = npix
           if (npixx.le.0) npixx = 500 ! default for other uses of npixx if auto pixels are used
        endif

        !!--page setup preliminaries
        if (usesquarexy) then
           just = 1  ! x and y axis have same scale
           ! unless 1D xsec through 2D data or non-cartesian
           if ((irender.gt.ndim .and. ndim.eq.2 .and. x_sec) &
               .or.(icoordsnew.gt.1 .and. .not.(coord_is_length(iplotx,icoordsnew) &
                                          .and. coord_is_length(iploty,icoordsnew)))) then
              just = 0
           endif
        else
           just = 0
        endif
        !--work out if colour bar is going to be plotted
        !  (leave space in page setup if so)
        iPlotColourBar = .false.
        if (irender.gt.ndim .and..not.(ndim.eq.2.and.x_sec)) iPlotColourBar = (iColourBarStyle.gt.0)


        !!--work out coordinate that is not being plotted
        iz = 0
        if (ndim.ge.3) then
           do j=1,ndim
              if ((iplotx.ne.iploty).and. &
                  (ix(j).ne.iplotx).and.(ix(j).ne.iploty)) iz = ix(j)
           enddo
        endif
        iplotz = iz ! this is used as cross sectioned quantity
        if (iplotz.gt.0 .and. iplotz.le.ndataplots) then
           zplot(1:ntoti) = dat(1:ntoti,iplotz)
           labelz = label(iplotz)
        endif
        if (debugmode) print*,'DEBUG: iplotz = ',iplotz

        if (.not.interactivereplot) then
           irerender = .false.
        endif
        !
        !--rotate the particles about the z (and y and x) axes
        !  only applies to particle plots at the moment
        !
        if (ndim.ge.2 .and. (irotate .or. (ndim.eq.3 .and.use3Dperspective)) &
            .and. icoordsnew.eq.1) then
           if ((irotate .and. (angletempx.gt.tiny(0.) .or. angletempy.gt.tiny(0.))) &
               .or.(ndim.eq.3 .and.use3Dperspective .and. dzscreentemp.gt.tiny(0.))) then
              if (iaxis.ge.0) iaxistemp = -3
           endif
           ivectemp = 0

           !--for vector plots with rotation, need to allocate temporary
           !  arrays to hold the rotated vector components
           if (ivectorplot.gt.0) then
              ichangesize = .false.
              if (allocated(vecplot)) then
                 if (size(vecplot(1,:)).lt.ninterp) ichangesize = .true.
              endif
              if (.not.allocated(vecplot) .or. ichangesize) then
                 if (allocated(vecplot)) deallocate(vecplot)
                 allocate(vecplot(ndim,ninterp),stat=ierr)
                 if (ierr /= 0) then
                    print "(a)",' ERROR allocating memory for vector plot + rotation '
                    stop
                 endif
              endif
              ivectemp = ivectorplot
           endif
           call rotationandperspective(angletempx,angletempy,angletempz,dzscreentemp,zobservertemp, &
                xplot,yplot,zplot,ntoti,iplotx,iploty,iplotz,dat,ivectemp,vecplot,itrackpart)
           !--adapt plot limits after rotations have been done
           if (.not.interactivereplot) then
              call adapt_limits(iplotx,xplot,xmin,xmax,xminadapti,xmaxadapti,'x',&
                                iamtype,ntoti,npartoftype,iusetype,ipagechange)
              call adapt_limits(iploty,yplot,ymin,ymax,yminadapti,ymaxadapti,'y',&
                                iamtype,ntoti,npartoftype,iusetype,ipagechange)
           endif
           !!-reset co-ordinate plot limits if particle tracking
           if (itrackpart.gt.0 .and. .not.interactivereplot) then
              call settrackinglimits(itrackpart,iplotx,xplot,xmin,xmax)
              call settrackinglimits(itrackpart,iploty,yplot,ymin,ymax)
           endif
        endif

        !------------------------------------------------------------------
        !  rendering setup and interpolation (this is the rendering done
        !  *before* the cross sections are taken, e.g. to 3D grid)
        !------------------------------------------------------------------
        if ((irenderplot.gt.ndim).and. &
             ((ndim.eq.3).or.(ndim.eq.2.and..not.x_sec))) then

           !!--determine number of pixels in rendered image (npix = pixels in x direction)
           if (npix.gt.0) then
              npixx = npix
              call page_setup(dummy=.true.) ! do this here in case limits are auto-adjusted
              pixwidth  = (xmax-xmin)/real(npix)
              if (just.eq.1) then
                 pixwidthy = pixwidth
              else
                 pixwidthy = pixwidth*(ymax-ymin)/(xmax - xmin)
              endif
           else
           !!--automatically reset the pixel number to match the device
              call page_setup(dummy=.true.) !--npixx and npixy are determined here
              pixwidth = (xmax-xmin)/real(npixx)
              if (just.eq.1) then
                 pixwidthy = pixwidth
              else
                 pixwidthy = (ymax-ymin)/real(npixy)
              endif
           endif
           npixx = max(int((1. - epsilon(0.))*(xmax-xmin)/pixwidth) + 1,1)
           npixy = max(int((1. - epsilon(0.))*(ymax-ymin)/pixwidthy) + 1,1)

           !!--only need z pixels if working with interpolation to 3D grid
           !  (then number of z pixels is equal to number of cross sections)
           if ((ndim.ge.3).and.(x_sec.and.nxsec.gt.2)) then
              zmin = lim(iplotz,1)
              npixz = nxsec
           endif

           if (.not.interactivereplot .or. irerender) then
              if (allocated(datpix)) then
                 if (npixx.ne.size(datpix(:,1)) .or. npixy.ne.size(datpix(1,:))) then
                    deallocate(datpix)
                    allocate (datpix(npixx,npixy))
                    if (ndim.eq.3 .and. use3Dperspective .and. use3Dopacityrendering) then
                       if (allocated(brightness)) deallocate(brightness)
                       allocate(brightness(npixx,npixy))
                    endif
                    if (icontourplot.gt.ndim) then
                       if (allocated(datpixcont)) deallocate(datpixcont)
                       allocate(datpixcont(npixx,npixy))
                    endif
                 endif
              else
                 allocate (datpix(npixx,npixy))
                 if (ndim.eq.3 .and. use3Dperspective .and. use3Dopacityrendering) then
                    if (allocated(brightness)) deallocate(brightness)
                    allocate(brightness(npixx,npixy))
                 endif
                 if (icontourplot.gt.ndim) then
                    if (allocated(datpixcont)) deallocate(datpixcont)
                    allocate(datpixcont(npixx,npixy))
                 endif
              endif

              select case(ndim)
              case(2)
                 !!--interpolate to 2D grid
                 !!  allocate memory for rendering array
                 if (.not. x_sec) then
                    call interpolate2D(xplot(1:ninterp),yplot(1:ninterp), &
                         hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                         icolourme(1:ninterp),ninterp,xmin,ymin,datpix,npixx,npixy, &
                         pixwidth,pixwidthy,inormalise,isperiodicx,isperiodicy)
                    !--also get contour plot data
                    if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                       if (.not.isameweights) & ! set contouring weights as necessary
                          call set_weights(weight,dat,iamtype,UseTypeInContours)

                       call interpolate2D(xplot(1:ninterp),yplot(1:ninterp), &
                            hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,icontourplot), &
                            icolourme(1:ninterp),ninterp,xmin,ymin,datpixcont,npixx,npixy, &
                            pixwidth,pixwidthy,inormalise,isperiodicx,isperiodicy)
                       gotcontours = .true.

                       if (.not.isameweights) & ! reset weights
                          call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                    endif
                 endif
              case(3)
                 !!--interpolation to 3D grid - then take multiple cross sections/projections
                 !!  do this if taking more than 2 cross sections, otherwise use fast xsec
                 if (x_sec.and.nxsec.gt.2) then
                    !!--allocate memory for 3D rendering array
                    if (allocated(datpix3D)) deallocate(datpix3D)
                    allocate ( datpix3D(npixx,npixy,npixz) )
                    !!--interpolate from particles to 3D grid
                    call interpolate3D(xplot(1:ninterp),yplot(1:ninterp), &
                         zplot(1:ninterp),hh(1:ninterp),weight(1:ninterp), &
                         dat(1:ninterp,irenderplot),icolourme(1:ninterp), &
                         ninterp,xmin,ymin,zmin,datpix3D,npixx,npixy,npixz,pixwidth,dz, &
                         inormalise,isperiodicx,isperiodicy,isperiodicz)

                    if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                       !!--allocate memory for 3D contouring array
                       if (allocated(datpixcont3D)) deallocate(datpixcont3D)
                       allocate ( datpixcont3D(npixx,npixy,npixz) )

                       if (.not.isameweights) & ! set contouring weights as necessary
                          call set_weights(weight,dat,iamtype,UseTypeInContours)

                       !!--interpolate from particles to 3D grid
                       call interpolate3D(xplot(1:ninterp),yplot(1:ninterp), &
                            zplot(1:ninterp),hh(1:ninterp),weight(1:ninterp), &
                            dat(1:ninterp,icontourplot),icolourme(1:ninterp), &
                            ninterp,xmin,ymin,zmin,datpixcont3D,npixx,npixy,npixz,pixwidth,dz, &
                            inormalise,isperiodicx,isperiodicy,isperiodicz)
                       gotcontours = .true.

                       if (.not.isameweights) & ! reset weights
                          call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                    endif
                 endif
              end select
           endif

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
              if (flythru) zslicepos = zslicepos + dz
              !!--for cross sections of particle plots, need range of co-ordinates in which
              !!  particles may lie
              zslicemin = zslicepos-0.5*dz
              zslicemax = zslicepos+0.5*dz
           endif

           !------------take projections/cross sections through 3D data-----------------!
           if (irenderplot.gt.0 .and. ndim.eq.3) then

              !!--allocate memory for 2D rendered array
              
              if (.not.interactivereplot) then
                 if (allocated(datpix)) then
                    if (npixx.ne.size(datpix(:,1)) .or. npixy.ne.size(datpix(1,:))) then
                       deallocate(datpix)
                       if (debugmode) print*,'reallocating datpix...'
                       allocate ( datpix(npixx,npixy) )
                    endif
                 else
                    if (debugmode) print*,'allocating datpix...'
                    allocate ( datpix(npixx,npixy) )
                 endif
                 if (icontourplot.gt.ndim) then
                    if (allocated(datpixcont)) then
                       if (npixx.ne.size(datpixcont(:,1)) .or. npixy.ne.size(datpixcont(1,:))) then
                          deallocate(datpixcont)
                          if (debugmode) print*,'reallocating datpixcont...'
                          allocate ( datpixcont(npixx,npixy) )
                       endif
                    else
                       if (debugmode) print*,'allocating datpixcont...'
                       allocate ( datpixcont(npixx,npixy) )                    
                    endif
                 endif
              endif

              !------------------------------------------------------------------------
              ! if we have rendered to a 3D grid, take cross sections from this array
              !------------------------------------------------------------------------
              if (x_sec .and. nxsec.gt.2) then
                 ipixxsec = int(0.99999*(zslicepos-zmin)/dz) + 1
                 if (ipixxsec.gt.npixz) ipixxsec = npixz
                 print*,TRIM(label(iplotz)),' = ',zslicepos, &
                      ' cross section, pixel ',ipixxsec
                 datpix = datpix3D(:,:,ipixxsec)    ! slices are in 3rd dimension
                 if (gotcontours) then
                    datpixcont = datpixcont3D(:,:,ipixxsec)
                 endif
              else
                 !-------------------------------------------------------------------
                 !  or do a fast projection/cross section of 3D data to 2D array
                 !-------------------------------------------------------------------

                 !--only rerender if absolutely necessary
                 if (.not.interactivereplot .or. irerender) then
                    if (x_sec) then
                       if (use3Dperspective .and. use3Dopacityrendering) then
                          !!--do surface-rendered cross-section with opacity
                          if (iverbose > 0) print*,trim(label(ix(iplotz))),' = ',zslicepos,  &
                                ' : opacity-rendered cross section', xmin,ymin
                          if (ipmass.gt.0) then
                             if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),weight(1:ninterp),&
                                  dat(1:ninterp,icontourplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,zslicepos)
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),weight(1:ninterp), &
                               dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,zslicepos)
                          else
                             if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  pmassav,1,hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,icontourplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,zslicepos)
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               pmassav,1,hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,zslicepos)
                          endif
                       elseif (use3Dperspective) then
                          print*,'ERROR: X_SEC WITH 3D PERSPECTIVE NOT IMPLEMENTED'
                          datpix = 0.
                       else
                          !!--do fast cross-section
                          if (iverbose > 0) print*,trim(label(ix(iplotz))),' = ',zslicepos,  &
                                ' : fast cross section', xmin,ymin
                          call interpolate3D_fastxsec( &
                               xplot(1:ninterp),yplot(1:ninterp), &
                               zplot(1:ninterp),hh(1:ninterp), &
                               weight(1:ninterp),dat(1:ninterp,irenderplot),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,zslicepos,datpix,npixx,npixy,pixwidth, &
                               pixwidthy,inormalise)
                          !!--same but for contour plot
                          if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                             if (.not.isameweights) & ! set contouring weights as necessary
                                call set_weights(weight,dat,iamtype,UseTypeInContours)

                             call interpolate3D_fastxsec( &
                                  xplot(1:ninterp),yplot(1:ninterp), &
                                  zplot(1:ninterp),hh(1:ninterp), &
                                  weight(1:ninterp),dat(1:ninterp,icontourplot),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,zslicepos,datpixcont,npixx,npixy,pixwidth, &
                                  pixwidthy,inormalise)
                             gotcontours = .true.

                             if (.not.isameweights) & ! reset weights
                                call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                          endif
                       endif
                    else
                       if (use3Dperspective .and. use3Dopacityrendering) then
                          !!--do fast projection with opacity
                          if (ipmass.gt.0) then
                             !--contour plot first
                             if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),weight(1:ninterp),&
                                  dat(1:ninterp,icontourplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,huge(zslicepos))
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),&
                               weight(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,huge(zslicepos))
                          else
                             !--do contour plot first so brightness corresponds to render plot
                             if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  pmassav,1,hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,huge(zslicepos))
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               pmassav,1,hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,huge(zslicepos))
                          endif
                       else
                          !!--do fast projection of z integrated data (e.g. column density)
                          if (icoordsnew.ne.icoords) then
                             call interpolate3D_proj_geom( &
                                  dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)), &
                                  hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                                  icolourme(1:ninterp),ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth, &
                                  pixwidthy,inormalise,icoordsnew,iplotx,iploty,iplotz,ix)
                          else
                             call interpolate3D_projection( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                                  icolourme(1:ninterp),ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth, &
                                  pixwidthy,inormalise,zobservertemp,dzscreentemp,ifastrender)
                          endif
                          !!--same but for contour plot
                          if (icontourplot.gt.0 .and. icontourplot.le.numplot) then
                             if (.not.isameweights) & ! set contouring weights as necessary
                                call set_weights(weight,dat,iamtype,UseTypeInContours)

                             call interpolate3D_projection( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,icontourplot), &
                                  icolourme(1:ninterp),ninterp,xmin,ymin,datpixcont,npixx,npixy,pixwidth, &
                                  pixwidthy,inormalise,zobservertemp,dzscreentemp,ifastrender)
                             gotcontours = .true.

                             if (.not.isameweights) & ! reset weights
                                call set_weights(weight,dat,iamtype,UseTypeInRenderings)
                          endif
                          !!--adjust the units of the z-integrated quantity
                          if (iRescale .and. units(ih).gt.0. .and. .not.inormalise) then
                             datpix = datpix*(unitzintegration/units(ih))
                             if (gotcontours) then
                                datpixcont = datpixcont*(unitzintegration/units(ih))
                             endif
                          endif
                       endif
                    endif
                 endif

              endif ! whether 3D grid or fast renderings

              !-------------take cross sections through 2D data------------------!
           elseif (irenderplot.gt.0 .and. ndim.eq.2 .and. x_sec) then
              !-------------------------------------------------------------------
              !  or do a fast cross section through 2D data to 1D array
              !-------------------------------------------------------------------
              !!--interpolate from 2D data to 1D line
              !!  line is specified by giving two points, (x1,y1) and (x2,y2)
              !--set up 1D grid and allocate memory for datpix1D
              if (.not.interactivereplot) then
                 xmin = 0.   ! distance (r) along cross section
                 xmax = SQRT((xseclineY2-xseclineY1)**2 + (xseclineX2-xseclineX1)**2)
              endif
              dxgrid = (xmax-xmin)/REAL(npixx)
              call set_grid1D(xmin,dxgrid,npixx)

              call interpolate2D_xsec( &
                   dat(1:ninterp,iplotx),dat(1:ninterp,iploty),&
                   hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                   icolourme(1:ninterp),ninterp,xseclineX1,xseclineY1,xseclineX2,xseclineY2, &
                   datpix1D,npixx,inormalise)
              !
              !--find limits of datpix1D for plotting
              !  do transformations on rendered array where appropriate
              !  set these as ymin,ymax and set labels of plot
              !
              call transform(datpix1D,itrans(irenderplot))
              labely = transform_label(label(irenderplot),itrans(irenderplot))
              if (abs(xseclineY2-xseclineY1).gt.epsilon(0.)) then
                 labelx = 'cross section' ! only if cross-section is oblique (otherwise keep x axis label)
              endif
              !!--if adaptive limits, find limits of datpix
              if (.not.interactivereplot) then
                 ymin = minval(datpix1D)
                 ymax = maxval(datpix1D)
                 xminadapt(irenderplot) = min(ymin,xminadapt(irenderplot))
                 xmaxadapt(irenderplot) = max(ymax,xmaxadapt(irenderplot))
                 if (iadapt) then
                    print*,' adapting y limits'
                 else
                    !!--or use fixed limits and apply transformations
                    ymin = lim(irenderplot,1)
                    ymax = lim(irenderplot,2)
                    call transform_limits(ymin,ymax,itrans(irenderplot))
                 endif
              endif

           endif ! 2 or 3D and rendering

           !------------------------------------------------------------------
           !   apply transformations to, and find limits for the 2D
           !   pixel array datpix resulting from the interpolation operations
           !   do this *before* the page setup so that rendermin,max
           !   can be stored in page_setup for interactive plots
           !------------------------------------------------------------------
           if (irenderplot.gt.0 .and. irenderplot.le.numplot) then
              if (ndim.eq.3 .or. (ndim.eq.2 .and..not.x_sec)) then
                 !!--determine whether rendered quantity is logged or not
                 logged = islogged(itrans(irenderplot))
                 if (gotcontours) then
                    loggedcont = islogged(itrans(icontourplot))
                 else
                    loggedcont = .false.
                 endif
                 !!--do transformations on rendered array (but only the first time!)
                 if (.not.interactivereplot .or. irerender) then
                    if (logged) then
                       !!--if log, then set zero values to some large negative number
                       !   but exclude this value from adaptive limits determination
                       call transform(datpix,itrans(irenderplot),errval=error_in_log)
                    else
                       call transform(datpix,itrans(irenderplot))
                    endif
                    if (gotcontours) then
                       if (loggedcont) then
                          call transform(datpixcont,itrans(icontourplot),errval=error_in_log)
                       else
                          call transform(datpixcont,itrans(icontourplot))
                       endif
                    endif
                 endif

                 !!--set label for rendered quantity
                 labelrender = label(irenderplot)
                 if (gotcontours) labelcont = label(icontourplot)

                 !!--set label for column density (projection) plots
                 if (ndim.eq.3 .and..not. x_sec .and..not.(use3Dperspective.and.use3Dopacityrendering)) then
                    labelrender = integrate_label(labelrender,irender,iz,inormalise,iRescale,&
                                                  labelzintegration,projlabelformat,iapplyprojformat)
                    if (gotcontours) labelcont = integrate_label(labelcont,icontourplot,iz,inormalise,&
                                                 iRescale,labelzintegration,projlabelformat,iapplyprojformat)
                 endif
                 !!--apply transformations to the label(s) for the rendered and contoured quantit(y,ies)
                 labelrender = transform_label(labelrender,itrans(irenderplot))
                 if (gotcontours) labelcont = transform_label(labelcont,itrans(icontourplot))

                 !!--limits for rendered quantity
                 if (.not.interactivereplot .or. irerender) then
                    !!--find (adaptive) limits of rendered array
                    if (logged) then
                       renderminadapt = minval(datpix,mask=abs(datpix-error_in_log).gt.tiny(datpix)) ! see above
                    else
                       renderminadapt = minval(datpix)
                    endif
                    rendermaxadapt = maxval(datpix)
                    !--fix case where no limits are set due to NaNs etc.
                    if (renderminadapt.gt.rendermaxadapt) then
                       print "(a)",' WARNING: NaNs in rendered quantity'
                       renderminadapt = 0.
                       rendermaxadapt = 0.
                    endif

                    if (gotcontours) then
                       if (loggedcont) then
                          contminadapt = minval(datpixcont,mask=abs(datpixcont+666.).gt.tiny(datpixcont))
                       else
                          contminadapt = minval(datpixcont)
                       endif
                       contmaxadapt = maxval(datpixcont)
                       !--fix case where no limits are set due to NaNs etc.
                       if (contminadapt.gt.contmaxadapt) then
                          print "(a)",' WARNING: NaNs in contoured quantity'
                          contminadapt = 0.
                          contmaxadapt = 1.
                       endif
                    endif

                    if (.not.interactivereplot .and. .not.isetrenderlimits) then
                       if (iadapt) then
                          print*,'adapting render limits'
                          rendermin = renderminadapt
                          rendermax = rendermaxadapt
                       else
                          !!--or apply transformations to fixed limits
                          rendermin = lim(irenderplot,1)
                          rendermax = lim(irenderplot,2)
                          call transform_limits(rendermin,rendermax,itrans(irenderplot))
                       endif
                       if (gotcontours) then
                          if (iadapt) then
                             print*,'adapting contour limits'
                             contmin = contminadapt
                             contmax = contmaxadapt
                          elseif (icontourplot.eq.irenderplot .and. lim2set(icontourplot)) then
                             contmin = lim2(icontourplot,1)
                             contmax = lim2(icontourplot,2)
                             call transform_limits(contmin,contmax,itrans(icontourplot))
                          else
                             contmin = lim(icontourplot,1)
                             contmax = lim(icontourplot,2)
                             call transform_limits(contmin,contmax,itrans(icontourplot))
                          endif
                       endif
                    endif
                 endif
                 if (iplotcont .and. .not.gotcontours) then
                 !
                 !  this is the case where contoured quantity=rendered quantity
                 !  don't need to recalculate the pixel array but limits can be independent
                 !  => do this even during interactive replotting as rendermin,max can be changed
                 !     but contour limits should copy changes unless separate contour limits are set
                 !
                    contmin = rendermin
                    contmax = rendermax
                    if (lim2set(irenderplot) .and. .not.iadapt) then
                       contmin = lim2(irenderplot,1)
                       contmax = lim2(irenderplot,2)
                       call transform_limits(contmin,contmax,itrans(irenderplot))
                    endif
                 endif

                 !!  do not let max=0 on log plots as this is suspiciously wrong
                 if (logged) then
                    if (iadapt .and. abs(rendermax).lt.tiny(datpix)) then
                       !!print*,'max=0 on log plot, fixing'
                       rendermax = maxval(datpix)
                    endif
                 endif
                 if (gotcontours .and. loggedcont) then
                    if (iadapt .and. abs(contmax).lt.tiny(datpixcont)) then
                       contmax = maxval(datpixcont)
                    endif
                 endif
              endif

           !-------------------------------------------------------------------------
           !   similar but where particle colouring is used instead of interpolation
           !-------------------------------------------------------------------------
           elseif (irenderpart.gt.0 .and. iplotpart) then
              !--allocate memory for particle colouring
              if (.not.allocated(renderplot)) then
                 allocate(renderplot(ntoti),stat=ierr)
                 if (ierr /= 0) stop 'error allocating temporary array for particle colouring'
              endif

              !--apply transformations to render array and set label
              renderplot(1:ntoti) = dat(1:ntoti,irenderpart)
              call transform(renderplot(1:ntoti),itrans(irenderpart))
              labelrender = label(irenderpart)
              labelrender = transform_label(labelrender,itrans(irenderpart))

              call adapt_limits(irenderpart,renderplot(1:ntoti),rendermin,rendermax, &
                                renderminadapt,rendermaxadapt,trim(labelrender),&
                                iamtype,ntoti,npartoftype,iusetype,ipagechange)

              !!--limits for rendered quantity
              if (.not.interactivereplot .and. .not.isetrenderlimits) then
                 !!--find (adaptive) limits of rendered array
                 !   (note: something may be not quite right here with adapt during anim sequences)

                 if (iadapt) then
                    rendermin = renderminadapt
                    rendermax = rendermaxadapt
                 else
                    !!--use fixed limits and apply transformations
                    rendermin = lim(irenderpart,1)
                    rendermax = lim(irenderpart,2)
                    call transform_limits(rendermin,rendermax,itrans(irenderpart))
                 endif
              endif
              !
              !--actually colour the particles
              !
              call colour_particles(renderplot(1:ntoti), &
                   rendermin,rendermax,icolourme(1:ntoti),ntoti)

              !--deallocate memory
              if (allocated(renderplot)) deallocate(renderplot)

           endif

           !-------------------------------------------------------------------------
           !   similarly, get limits to be used in the vector plot so we can store
           !   them during page_setup for interactive plots
           !-------------------------------------------------------------------------
           if (ivectorplot.ne.0 .and. ndim.ge.2) then
             if (iamvec(ivectorplot).ne.0) then
                !!--choose quantity to be plotted
                ivecx = iamvec(ivectorplot) + iplotx - 1
                ivecy = iamvec(ivectorplot) + iploty - 1

                if (.not.interactivereplot) then ! not if vecmax changed interactively
                   if (iadapt) then
                      vecmax = -1.0  ! plot limits then set in vectorplot
                   else
                      vecmax = max(lim(ivecx,2),lim(ivecy,2))
                   endif
                endif
             endif
           endif

           !-----end of preliminary muff for 2D/3D cross sections/renderings ------------------

           !---------------------------------
           ! setup page
           !---------------------------------

           call page_setup

           !--add to log
           if (x_sec.and.iplotpart .and. iplotz.gt.0 .and. iverbose > 1) then
              print "(' cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)",&
                    label(iplotz),zslicemin,label(iplotz),zslicemax
           endif

           !------------------------------
           ! now actually plot the data
           !------------------------------
           if (irenderplot.gt.0) then
              if ((ndim.eq.3).or.(ndim.eq.2.and. .not.x_sec)) then

                 !!--call subroutine to actually render the image
                 if (gotcontours .and. double_rendering) then
                    !--if double rendering, plot first image in greyscale
                    call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                      npixx,npixy,xmin,ymin,pixwidth,pixwidthy, &
                      1,iplotcont,0,ncontours,.false.,ilabelcont,contmin,contmax)
                 else
                    if (use3Dperspective .and. use3Dopacityrendering .and. ndim.eq.3 .and. writeppm) then
                       call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                         npixx,npixy,xmin,ymin,pixwidth,pixwidthy,    &
                         icolours,iplotcont,0,ncontours,.false.,ilabelcont,contmin,contmax,alpha=brightness)                    
                    else
                       call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                         npixx,npixy,xmin,ymin,pixwidth,pixwidthy,    &
                         icolours,iplotcont,0,ncontours,.false.,ilabelcont,contmin,contmax)
                    endif
                 endif

                 !!--contour/2nd render plot of different quantity on top of 1st rendering
                 if (gotcontours) then
                    if (double_rendering) then
                       call render_pix(datpixcont,contmin,contmax,trim(labelcont), &
                            npixx,npixy,xmin,ymin,pixwidth,pixwidthy,icolours,.false.,&
                            0,ncontours,.false.,ilabelcont,transparent=.true.)
                    else
                       call render_pix(datpixcont,contmin,contmax,trim(labelcont), &
                            npixx,npixy,xmin,ymin,pixwidth,pixwidthy,0,.true.,0,ncontours,&
                            .false.,ilabelcont)
                    endif
                 endif

                 PlotOnRender_tmp(:) = PlotOnRenderings(:)
                 isinktype = get_sink_type(ntypes)
                 if (use3Dperspective .and. use3Dopacityrendering .and. rendersinks .and. isinktype > 0) then
                    PlotOnRender_tmp(isinktype) = .false.
                 endif
                 
                 !!--write ppm if interpolate3D_opacity
                 if ((.not. plotlib_supports_alpha) .and. &
                     (use3Dperspective .and. use3Dopacityrendering .and. ndim.eq.3 .and. writeppm)) then
                    !!--plot non-gas particle types (e.g. sink particles) on top (and to ppm)
                    call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                      zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                      icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRender_tmp(:), &
                      (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                      xmin,xmax,ymin,ymax,ifastparticleplot,datpix,npixx,npixy,rendermax,brightness)

                    call write_pixmap_ppm(datpix,npixx,npixy,xmin,ymin,pixwidth,rendermin,rendermax, &
                                       trim(labelrender),((istep-1)*nframesloop + iframe),brightness)
                 !!--dump pixmap to file if option set
                 elseif (iwritepixmap) then

                    !!--plot non-gas particle types (e.g. sink particles) on top (and to ppm)
                    call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                      zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                      icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRender_tmp(:), &
                      (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                      xmin,xmax,ymin,ymax,ifastparticleplot,datpix,npixx,npixy,rendermax)

                    call writepixmap(datpix,npixx,npixy,xmin,ymin,pixwidth,rendermin,rendermax,labelrender,&
                                     unitslabel(irenderplot),((istep-1)*nframesloop + iframe),x_sec,rootname(ifileopen))
                 !!--no ppm write
                 else
                    !!--plot non-gas particle types (e.g. sink particles) on top
                    call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                      zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                      icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRender_tmp(:), &
                      (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                      xmin,xmax,ymin,ymax,ifastparticleplot)
                 endif

              elseif (ndim.eq.2 .and. x_sec) then
                 !---------------------------------------------------------------
                 ! plot 1D cross section through 2D data (contents of datpix)
                 !---------------------------------------------------------------
                 call plot_line(npixx,xgrid,datpix1D)
              endif
           else
              !-----------------------
              ! particle plots
              !-----------------------
              if (iplotpart) then
                 if (debugmode .and. size(icolourme).ge.10) &
                    print*,'DEBUG: starting particle plot with ',ntoti,' particles ',&
                           zplot(1:10),icolourme(1:10),npartoftype(:),iusetype(:)
                 !!--plot all particle types
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),iamtype,npartoftype(:),iusetype(:), &
                   (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                   xmin,xmax,ymin,ymax,ifastparticleplot)
              else
                 !!--plot non-gas particle types on top of vector plots (e.g. sinks)
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRenderings(:), &
                   (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                   xmin,xmax,ymin,ymax,ifastparticleplot)

              endif
           endif

           !--------------------------------------------------------------
           ! vector maps (can be on top of particle plots and renderings)
           !--------------------------------------------------------------
           if (ivecx.gt.0 .and. ivecy.gt.0 .and. ivectorplot.gt.0) then
              if (iRescale) then
                 labelvecplot = trim(labelvec(ivectorplot))//trim(unitslabel(ivectorplot))
              else
                 labelvecplot = trim(labelvec(ivectorplot))
              endif
              !!--set label for projection plots (2268 or 2412 for integral sign)
              if (ndim.eq.3 .and..not. x_sec) then
                 if (iRescale) then
                    labelvecplot = '\(2268) '//trim(labelvecplot)//' d'// &
                      trim(label(iz)(1:index(label(iz),unitslabel(iz))-1))//' [code units]'
                 else
                    labelvecplot = '\(2268) '//trim(labelvecplot)//' d'//trim(label(iz))
                 endif
              endif
              pixwidthvec  = (xmax-xmin)/real(npixvec)
              if (just.eq.1) then
                 pixwidthvecy = pixwidthvec !(xmax-xmin)/real(npixvec)
              else
                 pixwidthvecy = pixwidthvec
              endif
              npixyvec = int(0.999*(ymax-ymin)/pixwidthvecy) + 1
              pixwidth = (xmax-xmin)/real(npixx) ! used in synchrotron plots

              if (.not.ihavesetweights) then
                 call set_weights(weight,dat,iamtype,(iusetype.and.UseTypeInRenderings))
              endif

              call vector_plot(ivecx,ivecy,npixvec,npixyvec,pixwidthvec,pixwidthvecy,vecmax,labelvecplot,got_h)

              !--vecmax is returned with the adaptive value if sent in -ve
              !  store this for use in interactive_multi
              if (xmaxmulti(ivectorplot).lt.0.) xmaxmulti(ivectorplot) = vecmax
           endif

           !---------------------------------
           ! plot rotated axes
           !---------------------------------
           if (irotate .and. irotateaxes.gt.0 .and. icoordsnew.eq.1) then
              call rotatedaxes(irotateaxes,iplotx,iploty,angletempx,angletempy,angletempz, &
                               dzscreentemp,zobservertemp)
           endif
           !
           !--redraw axes over what has been plotted
           !
           if (irenderplot.gt.0 .or. plotlib_is_pgplot) then
              call redraw_axes(iaxistemp,just,yscalealt,itransy)
           endif
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title
           !
           !--plot exact solution if relevant (before going interactive)
           !
           if (iexact.ne.0 .and.nyplot.le.nacross*ndown .and. &
               ipanelselect(iPlotExactOnlyOnPanel,ipanel,irow,icolumn)) then
              iaxisy = iaxis
              if (tile_plots .and. icolumn.ne.1) iaxisy = -1
              if (iexact.eq.iexact_rochelobe .and. use_sink_data .and. ipmass > 0 .and. ndim >= 2) then
                 isinktype = get_sink_type(ntypes)
                 call locate_first_two_of_type(isink1,isink2,isinktype,iamtype,npartoftype,ntoti)
                 mprim = dat(isink1,ipmass)
                 msec  = dat(isink2,ipmass)
                 xprim(1) = xplot(isink1)
                 xprim(2) = yplot(isink1)
                 xsec(1) = xplot(isink2)
                 xsec(2) = yplot(isink2)
              endif
              call exact_solution(iexact,iplotx,iploty, &
                   itrans(iplotx),itrans(iploty),icoordsnew, &
                   ndim,ndimV,timei,xmin,xmax,gammai, &
                   xplot(1:ntoti),yplot(1:ntoti),icolourme(1:ntoti),iamtype,npartoftype,iusetype, &
                   pmassmin,pmassmax,ntoti,imarktype(1), &
                   units(iplotx),units(iploty),irescale,iaxisy)
           endif

           !--the following line sets the number of steps on page to nstepsonpage
           !  in the case where we reach the last timestep before nstepsonpage is reached
           !  (makes interactive replotting behave better)
           if (lastplot) istepsonpage = nstepsperpage

           !
           !--enter interactive mode
           !

           if (interactive) then
              if (nacross*ndown.eq.1 .and. (nstepsperpage.eq.1 .or. nsteps.eq.1)) then
                 iadvance = nfreq
                 call interactive_part(ntoti,iplotx,iploty,iplotz,irender,icontourplot,ivecx,ivecy, &
                      xplot(1:ntoti),yplot(1:ntoti),zplot(1:ntoti), &
                      hh(1:ntoti),icolourme(1:ntoti),iamtype,iusetype,npartoftype, &
                      xmin,xmax,ymin,ymax,rendermin,rendermax,renderminadapt,rendermaxadapt,contmin,contmax,&
                      contminadapt,contmaxadapt,vecmax, &
                      angletempx,angletempy,angletempz,ndim,xorigin(1:ndim),x_sec,zslicepos,dz, &
                      zobservertemp,dzscreentemp,use3Dopacityrendering,taupartdepthtemp,&
                      (double_rendering .and. gotcontours),irerender,itrackpart,icolours,&
                      iColourBarStyle,labelrender,iadvance,ipos,iendatstep,iframe,nframesloop,interactivereplot)
                 !--turn rotation on if necessary
                 if (abs(angletempx-anglex).gt.tol) irotate = .true.
                 if (abs(angletempy-angley).gt.tol) irotate = .true.
                 if (abs(angletempz-anglez).gt.tol) irotate = .true.
                 if (iadvance.eq.-666 .or. interactivereplot) exit over_frames
              elseif ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot) then
                 !
                 !--slightly different interactive mode if multiple plots on page
                 !
                 iadvance = nfreq
!                 call interactive_step(iadvance,ipos,iendatstep,xmin,xmax,ymin,ymax)
                 nplots = ipanel
                 irerender = .true.
                 call interactive_multi(iadvance,ipos,ifirststeponpage,iendatstep,iframe,nframefirstonpage, &
                      nframesloop,ipanel,iplotxtemp(1:nplots),iplotytemp(1:nplots),irendertemp(1:nplots),&
                      icontourtemp(1:nplots),ivecplottemp(1:nplots),double_rendering,xminmulti(:),xmaxmulti(:),&
                      vptxmin(1:nplots),vptxmax(1:nplots),vptymin(1:nplots),vptymax(1:nplots),barwmulti(1:nplots), &
                      xminadapt(:),xmaxadapt(:),nacross,ndim,xorigin(1:ndim),icolours,iColourBarStyle,interactivereplot)
                 if (iadvance.eq.-666 .or. interactivereplot) exit over_frames
              endif
           endif

           !
           !--%%%%%%%%%%%%% end loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
           !
        enddo over_cross_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! not both coordinates - these are just particle plots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     elseif ((.not.iscoordplot).and.(iploty.le.ndataplots .and. iplotx.le.ndataplots)) then

        if (debugmode) print*,'DEBUG: starting particle plot...'
        !
        !--sort out particle colouring
        !  (at present this is NOT used -can't render if not co-ord plot)
        !
        if (irenderpart.gt.0 .and. irenderpart.le.numplot) then
           !--allocate memory for particle colouring
           if (.not.allocated(renderplot)) then
              allocate(renderplot(ntoti),stat=ierr)
              if (ierr /= 0) stop 'error allocating temporary array for particle colouring'
           endif

           !--apply transformations to render array and set label
           renderplot(1:ntoti) = dat(1:ntoti,irenderpart)
           call transform(renderplot(1:ntoti),itrans(irenderpart))
           labelrender = label(irenderpart)
           labelrender = transform_label(labelrender,itrans(irenderpart))

           !--limits for rendered quantity
           if (.not.interactivereplot) then
              !--find (adaptive) limits of rendered array
              call adapt_limits(irenderpart,renderplot(1:ntoti),rendermin,rendermax, &
                                renderminadapt,rendermaxadapt,trim(labelrender),&
                                iamtype,ntoti,npartoftype,iusetype,ipagechange)
              if (.not.iadapt) then
                 !!--use fixed limits and apply transformations
                 rendermin = lim(irenderpart,1)
                 rendermax = lim(irenderpart,2)
                 call transform_limits(rendermin,rendermax,itrans(irenderpart))
              endif
           endif

           !--actually colour the particles
           call colour_particles(renderplot(1:ntoti), &
                rendermin,rendermax,icolourme(1:ntoti),ntoti)

           !--deallocate memory
           if (allocated(renderplot)) deallocate(renderplot)
        endif

        !--------------------------------
        ! setup page
        !--------------------------------
        just = 0
        call page_setup

        !--------------------------------
        ! now plot particles
        !--------------------------------
        !npixx = 512
        !npixy = 512
        !if (.not.allocated(datpix)) allocate(datpix(npixx,npixy))
        !pixwidth = (xmax - xmin)/npixx
        !pixwidthy = (ymax - ymin)/npixy
        !call interpolate2D_pixels(xplot,yplot,icolourme,ntoti,xmin,ymin,datpix,npixx,npixy,&
        !     pixwidth,pixwidthy,.true.,.false.,.false.)
        !print*,maxval(datpix),minval(datpix)
        !call render_pix(datpix,0.,maxval(datpix),'blah', &
        !     npixx,npixy,xmin,ymin,pixwidth,pixwidthy,2,.false.,,ncontours,&
        !     .false.,.false.)
        !deallocate(datpix)

        call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
             zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
             icolourme(1:ntoti),iamtype,npartoftype(:),iusetype,.false., &
             zslicemin,zslicemax,' ',xmin,xmax,ymin,ymax,ifastparticleplot)

        !--------------------------------
        ! plot error bars
        !--------------------------------
        if (iploterrbars) then
           call plot_qci(icolourprev)        ! query line style and colour
           call plot_sci(linecolourthisstep) ! set colour to current line
           !--y error bars
           if (ilocerrbars(iploty).gt.0 .and. ilocerrbars(iploty).le.ndataplots) then
              call plot_errorbarsy(ntoti,xplot,yplot,dat(:,ilocerrbars(iploty)),itransy)
           endif
           !--x error bars
           if (ilocerrbars(iplotx).gt.0 .and. ilocerrbars(iplotx).le.ndataplots) then
              call plot_errorbarsx(ntoti,xplot,yplot,dat(:,ilocerrbars(iplotx)),itransx)
           endif
           call plot_sci(icolourprev)        ! restore line colour
        endif

        !
        !--redraw axes over what has been plotted
        !
        if (plotlib_is_pgplot) call redraw_axes(iaxis,just,yscalealt,itransy)
        !
        !--annotate with time / marker legend and title
        !
        call legends_and_title
        !
        !--plot exact solution (after redrawn axis for residual plots)
        !
        if (iexact.ne.0 .and.nyplot.le.nacross*ndown .and. &
            ipanelselect(iPlotExactOnlyOnPanel,ipanel,irow,icolumn)) then
           iaxisy = iaxis
           if (tile_plots .and. icolumn.ne.1) iaxisy = -1
           call exact_solution(iexact,iplotx,iploty,itrans(iplotx),itrans(iploty), &
                icoordsnew,ndim,ndimV,timei,xmin,xmax,gammai, &
                xplot(1:ntoti),yplot(1:ntoti),icolourme(1:ntoti),iamtype,npartoftype,iusetype, &
                pmassmin,pmassmax,ntoti,imarktype(1), &
                units(iplotx),units(iploty),irescale,iaxisy)
        endif
        !
        !--enter interactive mode
        !--the following line sets the number of steps on page to nstepsonpage
        !  in the case where we reach the last timestep before nstepsonpage is reached
        !  (makes interactive replotting behave better)
        if (lastplot) istepsonpage = nstepsperpage

        if (interactive) then
           if (nacross*ndown.eq.1 .and. (nstepsperpage.eq.1 .or. nsteps.eq.1)) then
              iadvance = nfreq
              call interactive_part(ntoti,iplotx,iploty,0,irenderpart,0,0,0, &
                   xplot(1:ntoti),yplot(1:ntoti),zplot(1:ntoti), &
                   hh(1:ntoti),icolourme(1:ntoti),iamtype,iusetype,npartoftype, &
                   xmin,xmax,ymin,ymax,rendermin,rendermax,renderminadapt,rendermaxadapt,&
                   contmin,contmax,contminadapt,contmaxadapt,vecmax, &
                   angletempx,angletempy,angletempz,ndim,xorigin(1:ndim), &
                   dumxsec,dummy,dummy,dummy,dummy,.false.,dummy,.false.,irerender, &
                   itrackpart,icolours,iColourBarStyle,labelrender,iadvance,ipos,iendatstep,iframe,nframesloop,interactivereplot)
              if (iadvance.eq.-666 .or. interactivereplot) exit over_frames ! this should be unnecessary
           elseif ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot) then
              !
              !--timestep control only if multiple plots on page
              !
              iadvance = nfreq
!              call interactive_step(iadvance,ipos,iendatstep,xmin,xmax,ymin,ymax)
              nplots = ipanel
              irerender = .true.
              call interactive_multi(iadvance,ipos,ifirststeponpage,iendatstep,iframe,nframefirstonpage, &
                   nframesloop,ipanel,iplotxtemp(1:nplots),iplotytemp(1:nplots),irendertemp(1:nplots),&
                   icontourtemp(1:nplots),ivecplottemp(1:nplots),.false.,xminmulti(:),xmaxmulti(:),&
                   vptxmin(1:nplots),vptxmax(1:nplots),vptymin(1:nplots),vptymax(1:nplots),barwmulti(1:nplots), &
                   xminadapt(:),xmaxadapt(:),nacross,ndim,xorigin(1:ndim),icolours,iColourBarStyle,interactivereplot)
              if (iadvance.eq.-666 .or. interactivereplot) exit over_frames
           endif
        endif

     elseif (iploty.le.numplot) then! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! additional plots (not plots of particle data - e.g where some additional
! information is read from a file and plotted on the same page as the
! particle plots, or where some additional plot is calculated
! from the particle data, such as errors etc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (debugmode) print*,'DEBUG: starting extra plot...'
        !--------------------------------------------------------------
        !  plot surface density, Toomre Q parameter
        !  or Probability Distribution Function
        !  => these all involve a new "y column"
        !     but use a particle property as the x axis
        !--------------------------------------------------------------
        if (iploty.eq.isurfdens .or. iploty.eq.itoomre .or. iploty.eq.ipdf) then
           just = 0
           if (iploty.eq.itoomre) then
              itemp = 2
              label(iploty) = 'Q\dToomre\u'
              labely = trim(label(iploty))
           elseif (iploty.eq.isurfdens) then
              itemp = 1
              label(iploty) = '\gS'
              labely = trim(label(iploty))
           elseif (iploty.eq.ipdf) then
              label(iploty) = 'PDF ('//trim(label(iplotx))//')'
              labely = trim(label(iploty))
           endif
           yplot(:) = 0.

           if (itrans(iploty).ne.0) then
              labely = trim(transform_label(label(iploty),itrans(iploty)))
           endif

           if ((.not.interactivereplot) .or. irerender) then
              !--call routines which actually calculate disc properties from the particles
              if (iploty.eq.isurfdens .or. iploty.eq.itoomre) then
                 if (ispsound.gt.0 .and. ispsound.le.ndataplots) then
                    icol = ispsound   ! sound speed is present in data
                 elseif (iutherm.gt.0 .and. iutherm.le.ndataplots) then
                    icol = iutherm    ! use thermal energy if spsound not present
                 else
                    icol = 0
                 endif
                 ! work out the unit of mass, r needed for computing Toomre Q
                 unit_mass = 1.d0
                 unit_r    = 1.d0
                 unit_u    = 1.d0
                 if (iRescale) then
                    if (ix(1) > 0) unit_r = units(ix(1))
                    if (icol > 0)  unit_u = units(icol)
                 endif
                 if (ipmass.gt.0 .and. ipmass.le.ndataplots) then
                    if (iRescale) unit_mass = units(ipmass)
                    if (icol.gt.0) then
                       call disccalc(itemp,ntoti,xplot(1:ntoti),ntoti,dat(1:ntoti,ipmass),unit_mass,unit_r, &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty), &
                            icolourme(1:ntoti),iamtype,iusetype,npartoftype,gammai,unit_u,dat(1:ntoti,icol),icol.eq.ispsound)
                    else
                       call disccalc(itemp,ntoti,xplot(1:ntoti),ntoti,dat(1:ntoti,ipmass),unit_mass,unit_r, &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty), &
                            icolourme(1:ntoti),iamtype,iusetype,npartoftype,gammai)
                    endif
                 else
                    if (iRescale .and. irho > 0) unit_mass = units(irho)*unitzintegration**3
                    if (icol.gt.0) then
                       call disccalc(itemp,ntoti,xplot(1:ntoti),1,masstype(1),unit_mass,unit_r, &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty), &
                            icolourme(1:ntoti),iamtype,iusetype,npartoftype,gammai,unit_u,dat(1:ntoti,icol),icol.eq.ispsound)
                    else
                       call disccalc(itemp,ntoti,xplot(1:ntoti),1,masstype(1),unit_mass,unit_r, &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty),&
                            icolourme(1:ntoti),iamtype,iusetype,npartoftype,gammai)
                    endif
                 endif
              elseif (iploty.eq.ipdf) then
                 if (npdfbins.gt.0) then
                    ngrid = npdfbins
                 else  ! automatic number of bins determination
                    ngrid = int(0.75*ntoti**(1./3.))+1
                 endif
                 call set_grid1D(xmin,1.,ngrid)
              !--call routine which calculates pdf on the particles
                 !--compute PDF on raw (un-transformed) data
                 xplot(1:ntoti) = dat(1:ntoti,iplotx)
                 if (irho.gt.0 .and. irho.le.ndataplots .and. ipmass.gt.0 .and. ipmass.le.ndataplots) then
                    call pdf_calc(ntoti,xplot(1:ntoti),xmin,xmax,ngrid,xgrid,datpix1D, &
                         yminadapti,ymaxadapti,(npdfbins.gt.0),volweightedpdf, &
                         ierr,icolourme(1:ntoti),dat(1:ntoti,irho),dat(1:ntoti,ipmass))
                 else
                    call pdf_calc(ntoti,xplot(1:ntoti),xmin,xmax,ngrid,xgrid,datpix1D, &
                         yminadapti,ymaxadapti,(npdfbins.gt.0),volweightedpdf, &
                         ierr,icolourme(1:ntoti))
                 endif
                 !
                 !--write PDF to file
                 !
                 if (ierr.eq.0) then
                    call pdf_write(ngrid,xgrid,datpix1D,label(iplotx), &
                                   volweightedpdf,rootname(ifileopen),tagline)
                 endif
                 !
                 !--apply transformations to PDF data
                 !
                 if (itrans(iplotx).gt.0) then
                    !--reapply the x transform
                    call transform(xplot,itrans(iplotx))
                 endif
                 if (itrans(iploty).gt.0) then
                    call transform(datpix1D,itrans(iploty))
                    call transform_limits(yminadapti,ymaxadapti,itrans(iploty))
                 endif
              endif
           endif
           if (iadapt .and. .not.interactivereplot) then
              print "(1x,a)",'adapting '//trim(labely)//' limits'
              ymin = yminadapti
              ymax = ymaxadapti
           endif

           call page_setup

           call plot_qci(icolourprev)    ! query line style and colour
           call plot_qls(linestyleprev)
           ! set appropriate colour and style if multiple steps per page
           if (nstepsperpage.gt.1) then
              call plot_sci(linecolourthisstep)
              call plot_sls(linestylethisstep)
           endif
           if (iploty.eq.itoomre .or. iploty.eq.isurfdens) then
              call discplot()
           elseif (iploty.eq.ipdf) then
              !
              !--plot PDF as line segment, with blanking at zero
              !
              call plotline(size(xgrid),xgrid,datpix1D,blank=0.)
           endif

           !--restore line size and colour
           call plot_sci(icolourprev)
           call plot_sls(linestyleprev)

           if (plotlib_is_pgplot) call redraw_axes(iaxis,just,yscalealt,itransy)
           call legends_and_title

           !
           !--plot exact solution (after redrawn axis for residual plots)
           !
           if (iexact.ne.0 .and.nyplot.le.nacross*ndown .and. &
               ipanelselect(iPlotExactOnlyOnPanel,ipanel,irow,icolumn)) then
              iaxisy = iaxis
              if (tile_plots .and. icolumn.ne.1) iaxisy = -1
              call exact_solution(iexact,iplotx,iploty,itrans(iplotx),itrans(iploty), &
                   icoordsnew,ndim,ndimV,timei,xmin,xmax,gammai, &
                   xplot(1:ntoti),yplot(1:ntoti),icolourme(1:ntoti),iamtype,npartoftype,iusetype, &
                   pmassmin,pmassmax,ntoti,imarktype(1), &
                   units(iplotx),units(iploty),irescale,iaxisy)
           endif

           if (lastplot) istepsonpage = nstepsperpage
           if (interactive .and. ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot)) then
              iadvance = nfreq
              nplots = ipanel
              irerender = .true.
              call interactive_multi(iadvance,ipos,ifirststeponpage,iendatstep,iframe,nframefirstonpage, &
                   nframesloop,ipanel,iplotxtemp(1:nplots),iplotytemp(1:nplots),irendertemp(1:nplots),&
                   icontourtemp(1:nplots),ivecplottemp(1:nplots),.false.,xminmulti(:),xmaxmulti(:),&
                   vptxmin(1:nplots),vptxmax(1:nplots),vptymin(1:nplots),vptymax(1:nplots),barwmulti(1:nplots), &
                   xminadapt(:),xmaxadapt(:),nacross,ndim,xorigin(1:ndim),icolours,iColourBarStyle,interactivereplot)
              if (iadvance.eq.-666 .or. interactivereplot) exit over_frames
           endif
           cycle over_plots
        !--------------------------------------------------------------
        !  plot Toy star A-C plane solution
        !--------------------------------------------------------------
        elseif (iexact.eq.4 .and. iploty.eq.iacplane) then
        !
        !--A vs C for exact toystar solution
        !
           if (ndim.eq.1) then
              call exact_toystar_acplane(atstar,ctstar,sigma,gammai)
           elseif (ndim.eq.2) then
              call exact_toystar_acplane2D(atstar,ctstar,sigma,gammai)
           endif
           !--increment page counter as setpage is not called
           iplots = iplots + 1
           ipanel = ipanel + 1
           if (ipanel.gt.nacross*ndown) ipanel = 1

        !--------------------------------------------------------------
        !  power spectrum plots (uses x and data as yet unspecified)
        !--------------------------------------------------------------
        elseif (iploty.eq.ipowerspec) then

           labelx = 'frequency'
           labely = 'power'
        !
        !--3D: use FFT routines
        !
           if (ndim.eq.3) then
              if (.not.ihavesetweights) then
                 call set_weights(weight,dat,iamtype,iusetype)
              endif

              call powerspec3D_sph(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)), &
                   hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,ipowerspecy),icolourme(1:ninterp), &
                   ninterp,nfreqspec,lim(ipowerspecx,1),lim(ipowerspecx,2),xplot(1:nfreqspec), &
                   yplot(1:nfreqspec),inormalise)
              xmin = max(minval(xplot(1:nfreqspec)),1.0)
              xmax = maxval(xplot(1:nfreqspec))
              nfreqpts = nfreqspec
           else
        !
        !--1D: use slow FT routines or Lomb periodogram
        !
              if (.not.interactivereplot) then
                 xmin = freqmin  ! freq min
                 xmax = freqmax  ! freq max
              endif
              if (.not.interactivereplot .and. itrans(iploty).gt.0) then
                 call transform_limits(xmin,xmax,itrans(iploty))
              endif
              !
              !--setup frequency grid (evenly spaced in transformed grid)
              !
              nfreqpts = nfreqspec
              if (nfreqpts.ge.size(xplot)) then
                 nfreqpts = size(xplot)
                 print*,' WARNING: nfreqpts > array size, restricting to ',nfreqpts
              else
                 print "(a,i6)",' number of frequency points = ',nfreqpts
              endif
              dxfreq = (xmax - xmin)/real(nfreqpts)
              do i=1,nfreqpts
                 xplot(i) = xmin + (i-1)*dxfreq
              enddo
              !
              !--transform back to frequency space
              !
              if (itrans(iploty).gt.0) &
                 call transform_inverse(xplot(1:nfreqpts),itrans(iploty))

              if (.not.idisordered) then! interpolate first
                 !!--allocate memory for 1D grid (size = 2*npart)
                 ngrid = 2*npartoftype(1)
                 !!--set up 1D grid
                 xmingrid = lim(ipowerspecx,1)
                 xmaxgrid = lim(ipowerspecx,2)
                 dxgrid = (xmaxgrid-xmingrid)/ngrid
                 call set_grid1D(xmingrid,dxgrid,ngrid)

                 ninterp = ntoti
                 !!--interpolate to 1D grid
                 if (.not.ihavesetweights) then
                    call set_weights(weight,dat,iamtype,iusetype)
                 endif

                 call interpolate1D(dat(1:ninterp,ipowerspecx),hh(1:ninterp), &
                      weight(1:ninterp),dat(1:ninterp,ipowerspecy),icolourme(1:ninterp), &
                      ninterp,xmingrid,datpix1D,ngrid,dxgrid,inormalise)
                 !!--plot interpolated 1D data to check it
                 !!print*,minval(datpix1D),maxval(datpix1D)
                 !call pgswin(xmin,xmax,minval(datpix1D),maxval(datpix1D),0,1)
                 !call pgbox('BCNST',0.0,0,'1BVCNST',0.0,0)
                 !call pglabel('x',label(ipowerspecy),'1D interpolation')
                 !call pgline(ngrid,xgrid,datpix1D)
                 !read*
                 !call pgpage! change page

                 !!--call power spectrum calculation on the even grid
                 call powerspectrum(ngrid,xgrid,datpix1D,nfreqpts,xplot(1:nfreqpts), &
                                    yplot(1:nfreqpts),idisordered)
                 if (allocated(datpix1D)) deallocate(datpix1D)
                 if (allocated(xgrid)) deallocate(xgrid)
              else
                 !!--or else call power spectrum calculation on the particles themselves
                 call powerspectrum(ntoti,dat(1:ntoti,ipowerspecx), &
                      dat(1:ntoti,ipowerspecy),nfreqpts, &
                      xplot(1:nfreqpts),yplot(1:nfreqpts),idisordered)
              endif

           endif

           if (.not.interactivereplot) then
              ymin = minval(yplot(1:nfreqspec))
              ymax = maxval(yplot(1:nfreqspec))
           endif

           !!--uncomment next few lines to plot wavelengths instead
           !labelx = 'wavelength'
           !zplot(1:nfreqspec) = 1./xplot(1:nfreqspec)
           !xplot(1:nfreqspec) = zplot(1:nfreqspec)
           !if (.not.interactivereplot) then
           !   xmin = minval(xplot(1:nfreqspec))
           !   xmax = maxval(xplot(1:nfreqspec))
           !endif

           if (itrans(iploty).ne.0) then
              call transform(xplot(1:nfreqpts),itrans(iploty))
              labelx = transform_label(labelx,itrans(iploty))

              call transform(yplot(1:nfreqpts),itrans(iploty))
              labely     = transform_label(labely,itrans(iploty))
              if (.not.interactivereplot) then
                 call transform_limits(xmin,xmax,itrans(iploty))
                 call transform_limits(ymin,ymax,itrans(iploty))
              endif
           endif

           just = 0
           call page_setup

           call plot_qci(icolourprev)    ! query line style and colour
           call plot_qls(linestyleprev)
           if (nstepsperpage.gt.1) then
              call plot_sci(linecolourthisstep) ! set appropriate colour and style if multiple steps per page
              call plot_sls(linestylethisstep)
           endif

           call plot_line(nfreqpts,xplot(1:nfreqpts),yplot(1:nfreqpts))
           print*,' maximum power at '//trim(labelx)//' = ',xplot(maxloc(yplot(1:nfreqpts)))

           call plot_sci(icolourprev)
           call plot_sls(linestyleprev)

           !
           !--redraw axes over what has been plotted
           !
           if (plotlib_is_pgplot) call redraw_axes(iaxis,just,yscalealt,itransy)
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title

        elseif (iploty.eq.icolpixmap) then
        !--------------------------------------------------------------
        !  plot the contents of a pixel map read from a file
        !--------------------------------------------------------------
           !
           !--irender should already be set, associating the pixmap
           !  with a column from the SPH data. Then we can use the
           !  limit settings from the SPH data. Otherwise just
           !  treat it like a separate column.
           !
           if (irender.eq.0) irender = icolpixmap

           !--datpix is allocated inside the readpixmap routine
           if (allocated(datpix)) deallocate(datpix)
           
           if (irender.eq.icolpixmap) then
              labelrender = '|B_\phi|/|B_p|'
           else
              labelrender = label(irender)
           endif
           call readpixmap(datpix,npixx,npixy,rootname(ifileopen),&
                shortlabel(labelrender,unitslabel(irender)),istep,x_sec,ierr)

           if (.not.interactivereplot) then
              if (ndim.ge.1) then
                 xmin = lim(ix(1),1)
                 xmax = lim(ix(1),2)
              else
                 xmin = 0.
                 xmax = 1.
              endif
              if (ndim.ge.2) then
                 ymin = lim(ix(3),1)
                 ymax = lim(ix(3),2)
              else
                 ymin = 0.
                 ymax = 1.
              endif
           endif
           if (ndim.ge.1) iplotx = ix(1)
           if (ndim.ge.2) then
              iploty = ix(3)
              labely = label(ix(3))
           endif
           pixwidth = (xmax-xmin)/real(npixx)

           if (itrans(irender).ne.0 .and. allocated(datpix)) then
              call transform(datpix,itrans(irender),errval=error_in_log)
           endif
           labelrender = transform_label(labelrender,itrans(irender))

           !--find (adaptive) limits of rendered array
           if (allocated(datpix)) then
              renderminadapt = minval(datpix,mask=abs(datpix-error_in_log).gt.tiny(datpix))
              rendermaxadapt = maxval(datpix)
           endif
           !--limits for rendered quantity
           if (.not.interactivereplot) then
              if (.not.iadapt) then
                 !!--use fixed limits and apply transformations
                 rendermin = lim(irender,1)
                 rendermax = lim(irender,2)
                 call transform_limits(rendermin,rendermax,itrans(irender))
              endif
           endif

           just = 1
           iPlotColourBar = .true.
           call page_setup

           if (ierr.eq.0 .and. allocated(datpix)) then
              !!--call subroutine to actually render the image
              call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                npixx,npixy,xmin,ymin,pixwidth,pixwidth,    &
                icolours,iplotcont,0,0,.false.,.false.)
           endif
           !
           !--redraw axes over what has been plotted
           !
           if ((allocated(datpix) .and. ierr.eq.0) .or. plotlib_is_pgplot) then
              call redraw_axes(iaxis,just,yscalealt,itransy)
           endif
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title
           irender = 0
           iploty = icolpixmap
           iplotx = 0

        else
        !--------------------------------------------------------------
        !  plot the contents of an extra two-column ascii file
        !--------------------------------------------------------------

           call exact_fromfile('gwaves1.dat',xplot,yplot,1,2,nfreqpts,ierr)
           just = 0
           labelx = 't [ms]'
           labely = 'h'
           if (.not.interactivereplot) then
              xmin = minval(xplot(1:nfreqpts))
              xmax = maxval(xplot(1:nfreqpts))
              ymin = minval(yplot(1:nfreqpts))
              ymax = maxval(yplot(1:nfreqpts))
              !--adjust y axes
              ymin = (ymin + ymax)/2. - 0.55*(ymax-ymin)
              ymax = (ymin + ymax)/2. + 0.55*(ymax-ymin)
           endif
           call page_setup
           !--plot extra point corresponding to current time
           ipt = 0
           do i=1,nfreqpts-1
              if (xplot(i).le.timei .and. xplot(i+1).gt.timei) ipt = i
           enddo
           if (ipt.ne.0) then
              call plot_pt1(xplot(ipt),yplot(ipt),4)
              call plot_line(ipt,xplot(1:ipt),yplot(1:ipt))
           endif

           if (plotlib_is_pgplot) call redraw_axes(iaxis,just,yscalealt,itransy)
           call legends_and_title

        endif

        !--the following line sets the number of steps on page to nstepsonpage
        !  in the case where we reach the last timestep before nstepsonpage is reached
        !  (makes interactive replotting behave better)
        if (lastplot) istepsonpage = nstepsperpage

        if (interactive .and.((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) &
           .or. lastplot)) then
           iadvance = nfreq
           call interactive_step(iadvance,ipos,iendatstep,xmin,xmax,ymin,ymax,interactivereplot)
           irerender = .true.
           if (iadvance.eq.-666) exit over_frames
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if plot not in correct range
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else
        print*,' error in plotting : iplotx = ',iplotx,' iploty =',iploty, 'numplot =',numplot
        call plot_page! just skip to next plot

     endif ! ploty = whatever


  enddo over_plots ! over plots per timestep (nyplot)
  enddo over_frames ! over nframes for animation sequences

  !
  !--frame changing for multiple steps on the page (ie. page split into panels)
  !
  iframesave = 0
  if (iadvance.ne.-666 .and..not.interactivereplot) then
     if (nacross*ndown.gt.1 .and. insidesequence(iseqpos)) then
        !--for the last panel on the page, reset the sequence for next time
        if (ipanel.eq.nacross*ndown) then
           if (iframe.lt.nframes) then
              ipos = ipos - nacross*ndown
           endif
           iframesave = iframe
        else
           iframesave = iframe - 1
        endif
     endif
  endif

  !--free all dynamically allocated memory
  if (.not.interactivereplot) then
     if (allocated(datpix1D)) deallocate(datpix1D)
     if (allocated(datpix)) deallocate(datpix)
     if (allocated(brightness)) deallocate(brightness)
     if (allocated(datpix3D)) deallocate(datpix3D)
     if (allocated(xgrid)) deallocate(xgrid)
     if (allocated(vecplot)) deallocate(vecplot)
     if (allocated(datpixcont)) deallocate(datpixcont)
     if (allocated(datpixcont3D)) deallocate(datpixcont3D)
  endif

  !--free temporary arrays
  if (allocated(xplot)) deallocate(xplot)
  if (allocated(yplot)) deallocate(yplot)
  if (allocated(zplot)) deallocate(zplot)
  if (allocated(hh)) deallocate(hh)
  if (allocated(weight)) deallocate(weight)

  return

contains

!----------------------------------------------
! interfaces to the page setup routines
! this is called just before a plot is
! actually plotted
!----------------------------------------------
  subroutine page_setup(dummy)
    use colourbar,     only:get_colourbarmargins
    use pagesetup,     only:setpage2
    use settings_page, only:nstepsperpage,iUseBackgroundColourForAxes, &
                       vposlegend,iPlotLegend,usecolumnorder,&
                       xminpagemargin,xmaxpagemargin,yminpagemargin,ymaxpagemargin
    use settings_limits, only:adjustlimitstodevice
    use plotlib,       only:plot_qvp,plot_sci,plot_page,plotlib_is_pgplot,plot_set_opacity
    implicit none
    integer :: iplotsave,ipanelsave,ipanelpos,npanels_remaining
    real    :: barwidth, TitleOffset,xminmargin,xmaxmargin,yminmargin,ymaxmargin
    real    :: xminpix,xmaxpix,yminpix,ymaxpix,dxpix
    logical :: ipanelchange,dum,iprint_axes,lastrow
    logical, intent(in), optional :: dummy
    character(len=7) :: string
    !--------------------------------------------
    ! whether or not this is a dummy call or not
    !--------------------------------------------
    if (present(dummy)) then
       dum = dummy
       if (debugmode) print*,'DEBUG: entering page setup (dummy)'
    else
       dum = .false.
       if (debugmode) print*,'DEBUG: entering page setup'
    endif

    !---------------------
    ! increment counters
    !---------------------
    iplotsave = iplots
    ipanelsave = ipanel

    iplots = iplots + 1
    ipanelchange = .true.
    if (nstepsperpage.eq.0 .and. iplots.gt.1) ipanelchange = .false. ! this is an option to never change panels
    if (iplots.gt.1 .and. nyplots.eq.1 .and. nacross*ndown.gt.1.and..not.ipagechange) ipanelchange = .false.
    if (ipanelchange) ipanel = ipanel + 1
    if (ipanel.gt.nacross*ndown) ipanel = 1
    ipanel = max(ipanel,1) ! catch panel=0 if panel is not changing
    !--set counter for where we are in row, col
    if (.not.usecolumnorder) then
       irow      = ipanel - ((ipanel-1)/ndown)*ndown
       icolumn   = (ipanel-1)/ndown + 1
       ipanelpos = (irow-1)*nacross + icolumn
    else
       icolumn   = ipanel - ((ipanel-1)/nacross)*nacross
       irow      = (ipanel-1)/nacross + 1
       ipanelpos = ipanel
    endif
    !--if we are in interactive mode, use the currently buffered plot limits
    if (interactivereplot .and. (nacross*ndown.gt.1 .or. (nstepsperpage.gt.1 .and. nsteps.gt.1))) then
       xmin = xminmulti(iplotx)
       xmax = xmaxmulti(iplotx)
       ymin = xminmulti(iploty)
       ymax = xmaxmulti(iploty)
       if (ivectorplot.gt.0 .and. ivectorplot.le.numplot) then
          vecmax = xmaxmulti(ivectorplot)
       endif
    endif

    !nsteps_remaining = (nsteps - istep)/nstepsperpage
    !nplots_remaining = (nyplots - nyplot)/nstepsperpage
    npanels_remaining = (nsteps - istep)/nstepsperpage !nplots_remaining*nsteps_remaining

!    npanels_remaining = (nsteps - istep)*nyplots/nstepsperpage
    lastrow  = (npanels_remaining < nacross)

    lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) &
                .and. nyplot.eq.nyplots .and. k.eq.nxsec)

    !--------------------------------------------------------------
    ! output some muff to the screen
    !--------------------------------------------------------------
    if (((interactive .and. ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot)) &
        .or. (iadapt .and. (istepsonpage.eq.nstepsperpage .or. lastplot))) .and. .not.dum) then
       print*,trim(labelx),' min, max = ',xmin,xmax
       print*,trim(labely),' min, max = ',ymin,ymax
       if (irender.gt.0 .and. .not.(ndim.eq.2 .and. x_sec)) then
          print*,trim(labelrender),' min, max = ',rendermin,rendermax
          if (gotcontours) then
             print*,trim(labelcont),' min, max = ',contmin,contmax
          endif
       endif
    endif
    !--------------------------------------------------------------
    ! set up pgplot page
    !--------------------------------------------------------------
    !--use foreground colour
    if (.not.dum) call plot_sci(1)

    !--page margins: zero if no box is drawn
   ! xminmargin = 0.0
   ! xmaxmargin = 0.0
   ! yminmargin = 0.0
   ! ymaxmargin = 0.0
    xminmargin = xminpagemargin
    xmaxmargin = xmaxpagemargin
    yminmargin = yminpagemargin
    ymaxmargin = ymaxpagemargin

    !--leave space for colour bar if necessary (at end of row only on tiled plots)
    if ((tile_plots .and. iAllowspaceforcolourbar).or.(.not.tile_plots.and.iPlotColourBar)) then
       call get_colourbarmargins(iColourBarStyle,xmaxmargin,yminmargin,barwidth)
    else
       barwidth = 0.
    endif
    !--work out whether or not to leave space above plots for titles/legends
    TitleOffset = -tiny(Titleoffset)
    if (iPlotTitles .and. nstepsperpage.eq.1 .and. vpostitle.gt.0.) TitleOffset = vpostitle
    if (iPlotLegend .and. nstepsperpage.eq.1 .and. vposlegend.lt.0.) TitleOffset = max(Titleoffset,-vposlegend)

    inewpage = ipanel.eq.1 .and. ipanelchange .and. ipagechange
    if (inewpage .and. .not.dum) then
       call plot_page
       !--store ipos and nyplot positions for first on page
       !  as starting point for interactive replotting
       nyplotfirstonpage = nyplot
       ifirststeponpage = ipos
       nframefirstonpage = iframe
    endif
    !
    !--do not allow limits to be the same
    !
    if (abs(xmax-xmin).lt.tiny(xmax)) then
       if (.not.dum) print "(a)",' WARNING: '//trim(labelx)//'min='//trim(labelx)//'max '
       xmax = xmax + 1.0
       if (xmin.gt.0.) then
          xmin = max(xmin - 1.0,xmin,0.)
       else
          xmin = xmin - 1.0
       endif
    endif
    if (abs(ymax-ymin).lt.tiny(ymax)) then
       if (.not.dum) print "(a)",' WARNING: '//trim(labely)//'min='//trim(labely)//'max '
       ymax = ymax + 1.0
       if (ymin.gt.0.) then
          ymin = max(ymin - 1.0,ymin,0.)
       else
          ymin = ymin - 1.0
       endif
    endif
    if (irender.gt.0 .and. abs(rendermax-rendermin).lt.tiny(rendermax) .or.rendermax.ne.rendermax) then
       if (.not.dum) print "(a)",' WARNING: '//trim(labelrender)//'min='//trim(labelrender)//'max '
       rendermax = rendermax + 1.0
       if (rendermin.gt.0.) then
          rendermin = max(rendermin - 1.0,rendermin,0.)
       else
          rendermin = rendermin - 1.0
       endif
    endif
    if (debugmode) print*,'DEBUG: calling setpage...',nstepsperpage
    if (nstepsperpage.gt.0 .or. inewpage) then
       if (dum) then !--fake the page setup, then return
          if (.not.(interactivereplot .and. .not.irerender)) then
             call setpage2(ipanelpos,nacross,ndown,xmin,xmax,ymin,ymax, &
                     trim(labelx),trim(labely),'NOPGBOX',just,iaxistemp, &
                     xminmargin,xmaxmargin,yminmargin,ymaxmargin, &
                     0.0,TitleOffset,isamexaxis,tile_plots,adjustlimitstodevice, &
                     lastrow,lastplot,yscalealt,labelyalt,itransy)
             call plot_qvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
             if (debugmode) print*,'DEBUG: viewport xpix=',xminpix,'->',xmaxpix,' ypix=',yminpix,'->',ymaxpix

             npixx = max(nint(xmaxpix-xminpix),1)
             npixy = max(nint(ymaxpix-yminpix),1)
             if (debugmode) print*,'DEBUG: dx = ',xmax-xmin,' dy = ',ymax-ymin
             if (debugmode) print*,'DEBUG: dxpix = ',xmaxpix-xminpix,' dypix = ',ymaxpix-yminpix
             if (vectordevice .and. npixx.gt.1024) then
                npixx = 1024/nacross
                dxpix = (xmax-xmin)/npixx
                npixy = int(0.999*(ymax-ymin)/real(dxpix)) + 1
                print "(a,i4,a,i4,a)",' auto-selecting resolution of ',npixx,' x ',npixy,' for vector device'
                print "(a)",' => set the number of pixels manually if you want more (or less) than this.'
             else
                if (npix.eq.0 .and. debugmode) &
                   print "(a,i4,a,i4)",' auto-selecting device resolution = ',npixx,' x ',npixy
                !
                !--warn about PGPLOT limitations
                !
                if (plotlib_is_pgplot) then
                   if ((xmaxpix-xminpix).gt.1024. .or. (ymaxpix-yminpix).gt.1024) then
                      print "(/,75('*'))"
                      print "(a)",'!! WARNING: PGPLOT will truncate image to 1024 pixels on pixel devices.'
                      print "(a)",'!! To fix this, change line 18 of file grimg2.f in the PGPLOT source code:'
                      print "(a)",'!!          REAL     BUFFER(1026)'
                      print "(a)",'!! changing 1026 to something much bigger, then recompile PGPLOT.'
                      print "(75('*'),/)"
                   endif
                endif
             endif

          endif
          !--restore saved attributes
          iplots = iplotsave
          ipanel = ipanelsave
          if (debugmode) print*,'DEBUG: finished dummy page setup'
          return
       else
       
          !--if we are not changing page, do not reprint the axes
          iprint_axes = ipagechange .or. inewpage .or. &
                        ((iplots.le.nacross*ndown) .and. (nyplot.le.nacross*ndown .and. istepsonpage.eq.1))
         
          if (iprint_axes) then
             if (debugmode) print*,'DEBUG: printing axes ',ipagechange,inewpage,iplots,nyplot,istepsonpage
             string = ' '
          else
             if (debugmode) print*,'DEBUG: not printing axes ',ipagechange,inewpage,iplots,nyplot,istepsonpage
             string = 'NOPGBOX'
          endif
          call setpage2(ipanelpos,nacross,ndown,xmin,xmax,ymin,ymax, &
                  trim(labelx),trim(labely),string,just,iaxistemp, &
                  xminmargin,xmaxmargin,yminmargin,ymaxmargin, &
                  0.0,TitleOffset,isamexaxis,tile_plots,adjustlimitstodevice, &
                  lastrow,lastplot,yscalealt,labelyalt,itransy)
       endif
    endif

    if (debugmode) print*,'DEBUG: setpage ok, querying and saving viewport...'
    !--query and save viewport co-ordinates set up for this panel
    call plot_qvp(0,vptxmin(ipanel),vptxmax(ipanel),vptymin(ipanel),vptymax(ipanel))

    !--------------------------------------------------------------
    ! store current page setup for interactive mode on multiplots
    !--------------------------------------------------------------
    if (tile_plots) then
       barwmulti(ipanel) = 0.
    else
       barwmulti(ipanel) = barwidth
    endif
    iplotxtemp(ipanel) = iplotx
    iplotytemp(ipanel) = iploty
    irendertemp(ipanel) = irender
    icontourtemp(ipanel) = icontourplot
    ivecplottemp(ipanel) = ivectorplot
    xminmulti(iplotx) = xmin
    xmaxmulti(iplotx) = xmax
    xminmulti(iploty) = ymin
    xmaxmulti(iploty) = ymax
    if (irender.gt.0 .and. irender.le.size(xmaxmulti)) then
       xminmulti(irender) = rendermin
       xmaxmulti(irender) = rendermax
       if (icontourplot.gt.0 .and. gotcontours .and. icontourplot.le.size(xmaxmulti)) then
          xminmulti(icontourplot) = contmin
          xmaxmulti(icontourplot) = contmax
       endif
    endif
    if (ivectorplot.gt.0 .and. ivectorplot.le.numplot) then
       xmaxmulti(ivectorplot) = vecmax
    endif

    !
    ! store adaptive plot limits for a) in interactive mode
    ! on multiple plots per page
    !
    !--adaptive plot limits are allowed to change even during
    !  interactive replotting
    !if (.not.interactivereplot) then
       if (inewpage) then
          xminadapt = huge(xminadapt)
          xmaxadapt = -huge(xmaxadapt)
       endif
       xminadapt(iplotx) = min(xminadapt(iplotx),xminadapti)
       xmaxadapt(iplotx) = max(xmaxadapt(iplotx),xmaxadapti)
       xminadapt(iploty) = min(xminadapt(iploty),yminadapti)
       xmaxadapt(iploty) = max(xmaxadapt(iploty),ymaxadapti)
       if (irender.gt.0) then
          xminadapt(irender) = min(xminadapt(irender),renderminadapt)
          xmaxadapt(irender) = max(xmaxadapt(irender),rendermaxadapt)
          if (icontourplot.gt.0) then
             xminadapt(icontourplot) = min(xminadapt(icontourplot),contminadapt)
             xmaxadapt(icontourplot) = max(xmaxadapt(icontourplot),contmaxadapt)
          endif
       endif
    !endif

    !--change to background colour index for overlaid text and axes
    if (iUseBackGroundColourForAxes) then
       call plot_sci(0)
       call plot_set_opacity(1.) ! ensure background colour is opaque
    endif
    if (debugmode) print*,'DEBUG: finished page setup'

    return
  end subroutine page_setup

!------------------------------------------------------
! draws legend(s), titles etc
! (must be called after rendering otherwise rendering
!  will overwrite plot area)
!------------------------------------------------------
  subroutine legends_and_title
    use colourbar,     only:plotcolourbar
    use legends,       only:legend,legend_markers,legend_scale,ipanelselect
    use titles,        only:pagetitles,steplegend,lensteplegend
    use filenames,     only:nstepsinfile,nfiles,rootname
    use settings_page, only:iPlotLegend,iPlotStepLegend, &
        hposlegend,vposlegend,fjustlegend,legendtext,iPlotLegendOnlyOnPanel, &
        iPlotScale,iscalepanel,dxscale,hposscale,vposscale,scaletext,&
        alphalegend,iUseBackGroundColourForAxes
    use shapes,        only:nshapes,plot_shapes
    use pagesetup,     only:xlabeloffset
    use plotlib,       only:plot_qci,plot_sci,plot_annotate,plot_set_opacity
    use labels,        only:is_coord
    implicit none
    integer :: icoloursave
    character(len=lensteplegend) :: steplegendtext
    real :: xlabeloffsettemp
    integer :: ititle
    logical :: usebox

    !--save colour index
    call plot_qci(icoloursave)

    !--use foreground colour by default for legends
    call plot_sci(1)

    !--------------------------------------------------------------
    ! plot colour bar for rendered plots (use currently set colour)
    ! do this here so it always appears OVERLAID on the renderings
    !--------------------------------------------------------------
    if (irender.gt.0) then
       !--only plot colour bar at the end of first row on tiled plots
       if (tile_plots .and..not.(ipanel.eq.nacross*ndown .or. lastplot .or. &
           (OneColourBarPerRow.and.icolumn.eq.nacross))) iPlotColourBar = .false.

       if (iPlotColourBar) then
          xlabeloffsettemp = xlabeloffset + 1.0
          if (iaxistemp.lt.0) xlabeloffsettemp = 0.

          !--for tiled plots only on last plot in first row,
          !  and use full viewport size in the y direction
          if (tile_plots .and. .not.OneColourBarPerRow) then
             if (double_rendering .and. gotcontours) then
                call plotcolourbar(iColourBarStyle,icolours,contmin,contmax, &
                     trim(labelcont),.false.,xlabeloffsettemp, &
                     minval(vptxmin(1:ipanel)),maxval(vptxmax(1:ipanel)), &
                     minval(vptymin(1:ipanel)),maxval(vptymax(1:ipanel)))
             else
                call plotcolourbar(iColourBarStyle,icolours,rendermin,rendermax, &
                     trim(labelrender),.false.,xlabeloffsettemp, &
                     minval(vptxmin(1:ipanel)),maxval(vptxmax(1:ipanel)), &
                     minval(vptymin(1:ipanel)),maxval(vptymax(1:ipanel)))
             endif
          elseif (.not.tile_plots .or. (OneColourBarPerRow .and. icolumn.eq.nacross)) then
             !!--plot colour bar
             if (double_rendering .and. gotcontours) then
                !--for double rendering, plot the colour bar in the 2nd quantity
                call plotcolourbar(iColourBarStyle,icolours,contmin,contmax, &
                               trim(labelcont),.false.,xlabeloffsettemp)
             else
                call plotcolourbar(iColourBarStyle,icolours,rendermin,rendermax, &
                               trim(labelrender),.false.,xlabeloffsettemp)
             endif
          endif
       endif
    endif

    !--plot time on plot
    if (iPlotLegend .and. nyplot.eq.1 &
        .and. ipanelselect(iPlotLegendOnlyOnPanel,ipanel,irow,icolumn) &
        .and. timei.gt.-0.5*huge(timei)) then  ! but not if time has not been read from dump

       !--change to background colour index for legend text if overlaid
       if (iUseBackGroundColourForAxes .and. vposlegend.gt.0.) then
          call plot_sci(0)
          call plot_set_opacity(alphalegend)
       endif
       usebox = (ivectorplot.gt.0)
       if (istepsonpage.eq.1) then
          call legend(legendtext,timei,labeltimeunits,hposlegend,vposlegend,fjustlegend,usebox)
       endif
    endif

    !--line/marker style/colour legend for multiple timesteps on same page
    if (iPlotStepLegend .and. istepsonpage.gt.0 &
        .and.((nyplot.eq.1 .and. iPlotLegendOnlyOnPanel.eq.0) &
        .or. ipanelselect(iPlotLegendOnlyOnPanel,ipanel,irow,icolumn))) then

       !--change to background colour index for overlaid text and axes
       if (iUseBackGroundColourForAxes .and. vposlegend.gt.0.) call plot_sci(0)
       !
       !--use filenames in legend if none set
       !
       if (nstepsperpage.ge.1 .and. nsteplegendlines.ge.nstepsperpage*nacross*ndown) then
          steplegendtext = steplegend(istepsonpage + (ipanel-1)*nstepsperpage)
       elseif (istepsonpage.le.nsteplegendlines) then
          steplegendtext = steplegend(istepsonpage)
       elseif (all(nstepsinfile(1:nfiles).le.1)) then
          steplegendtext = trim(rootname(istep))
       else
          write(steplegendtext,"(a,i4)") 'step ',istep
       endif
       if (debugmode) print "(a,i2,a)",&
          ' DEBUG: plotting step legend (step ',istepsonpage,': "'//trim(steplegendtext)//'")'

       if (iploty.gt.ndataplots) then
          call legend_markers(istepsonpage,linecolourthisstep,imarktype(1),linestylethisstep, &
            .false.,.true.,trim(steplegendtext),hposlegend,vposlegend,1.0)
       else
          call legend_markers(istepsonpage,linecolourthisstep,imarktype(1),linestylethisstep, &
            iusetype(1),iplotline,trim(steplegendtext),hposlegend,vposlegend,1.0)
       endif
    endif

    !--use foreground colour by default for title
    call plot_sci(1)

    !--print title if appropriate
    if (iPlotTitles .and. istepsonpage.eq.1 .and. ipanel.le.ntitles) then
       if (ntitles.gt.nacross*ndown) then
          ititle = (ipos - 1)/nstepsperpage + 1
          if (ititle.gt.ntitles) ititle = ipanel
       else
          ititle = ipanel
       endif
       
       if (len_trim(pagetitles(ititle)).gt.0) then

          !--change to background colour index if title is overlaid
          if (iUseBackGroundColourForAxes .and. vpostitle.lt.0.) then
             call plot_sci(0)
             call plot_set_opacity(alphalegend)
          endif

          call plot_annotate('T',vpostitle,hpostitle,fjusttitle,trim(pagetitles(ititle)))
       endif
    endif

    !--use foreground colour by default for scale
    call plot_sci(1)

    !--scale on co-ordinate plots
    if (iPlotScale .and. (iscalepanel.eq.0 .or. ipanel.eq.iscalepanel) &
                   .and. is_coord(iplotx,ndim) .and. is_coord(iploty,ndim)) then

       !--change to background colour index if title is overlaid
       if (iUseBackGroundColourForAxes .and. vposscale.gt.0.) then
          call plot_sci(0)
          call plot_set_opacity(alphalegend)
       endif
       call legend_scale(dxscale,hposscale,vposscale,scaletext)
    endif

    !--plot shapes
    if (nshapes.gt.0 .and. istepsonpage.eq.1) &
       call plot_shapes(ipanel,irow,icolumn,itrans(iplotx),itrans(iploty),timei)

    !--restore colour index
    call plot_sci(icoloursave)
    call plot_set_opacity(1.0)

    return
  end subroutine legends_and_title

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
       xgrid(igrid) = xmin1D + (igrid-0.5)*dxgrid1D
    enddo

  end subroutine set_grid1D

!-------------------------------------------------------------------
! interface for setting limits when using particle tracking limits
!-------------------------------------------------------------------
  subroutine settrackinglimits(itrackpart,iplot,xploti,xmini,xmaxi)
    use labels,          only:is_coord
    use settings_limits, only:xminoffset_track,xmaxoffset_track
    implicit none
    integer, intent(in) :: itrackpart,iplot
    real, dimension(:), intent(in) :: xploti
    real, intent(inout) :: xmini,xmaxi

    !--particle tracking limits only apply to co-ordinate axes
    if (is_coord(iplot,ndim)) then
       xmini = xploti(itrackpart) - xminoffset_track(iplot)
       xmaxi = xploti(itrackpart) + xmaxoffset_track(iplot)
    endif

    return
   end subroutine settrackinglimits

!-------------------------------------------------------------------
! interface for setting interpolation weights
! (to make calls above neater)
!-------------------------------------------------------------------
  subroutine set_weights(weighti,dati,iamtypei,usetype)
    use settings_render,  only:idensityweightedinterpolation
    use interpolation,    only:set_interpolation_weights
    use settings_units,   only:unit_interp
    use settings_xsecrot, only:rendersinks,use3Dopacityrendering
    implicit none
    real, dimension(:), intent(out) :: weighti
    real, dimension(:,:), intent(in) :: dati
    integer(kind=int1), dimension(:), intent(in) :: iamtypei
    logical, dimension(:), intent(in) :: usetype

    ihavesetweights = .true.
    inormalise = inormalise_interpolations
    call set_interpolation_weights(weighti,dati,iamtypei,usetype,&
         ninterp,npartoftype,masstype,ntypes,ndataplots,irho,ipmass,ih,ndim,&
         iRescale,idensityweightedinterpolation,inormalise,units,unit_interp,required,&
         (use3Dopacityrendering .and. rendersinks))

    return
 end subroutine set_weights

!-------------------------------------------------------------------
! interface to vector plotting routines
! so that pixel arrays are allocated appropriately
!-------------------------------------------------------------------
  subroutine vector_plot(ivecx,ivecy,numpixx,numpixy,pixwidthvec,pixwidthvecy,vmax,label,got_h)
   use settings_vecplot, only:UseBackgndColorVecplot,iplotstreamlines,iplotarrowheads, &
       iplotsynchrotron,rcrit,zcrit,synchrotronspecindex,uthermcutoff, &
       ihidearrowswherenoparts,minpartforarrow,iVecplotLegend,iVecLegendOnPanel
   use interpolations2D, only:interpolate2D_vec
   use projections3D,    only:interpolate3D_proj_vec,interp3D_proj_vec_synctron
   use interpolate_vec,  only:mask_vectors,interpolate_vec_average
   use render,           only:render_vec
   use fieldlines,       only:streamlines,vecplot3D_proj
   use labels,           only:iutherm,is_coord
   use plotlib,          only:plot_qci,plot_qlw,plot_sci,plot_slw
   use system_utils,     only:lenvironment
   use legends,          only:ipanelselect
   implicit none
   integer,          intent(in) :: ivecx,ivecy,numpixx,numpixy
   real,             intent(in) :: pixwidthvec,pixwidthvecy
   real,          intent(inout) :: vmax
   character(len=*), intent(in) :: label
   logical,          intent(in) :: got_h
   real, dimension(numpixx,numpixy) :: vecpixx, vecpixy
   real, dimension(max(npixx,numpixx),max(npixy,numpixy)) :: datpixvec
   integer :: i,j,icoloursav,linewidthprev,ivecz
   real    :: vmag
   real    :: blankval,datmax
   logical :: usevecplot,use3Dstreamlines,plotlegend

   !--query colour index and line width
   call plot_qci(icoloursav)
   call plot_qlw(linewidthprev)

   !print*,'plotting vector field ',trim(label)
   if ((is_coord(ivecx,ndim) .or. ivecx.lt.0 .or.(ivecx.gt.ndataplots)) .or. &
       (is_coord(ivecy,ndim) .or. ivecy.lt.0 .or.(ivecy.gt.ndataplots))) then
      print*,'error finding location of vector plot in array'
   else
      use3Dstreamlines = (ndim.eq.3) .and. .not.x_sec !lenvironment('SPLASH_3DSTREAMLINES')

      !--plot arrows in either background or foreground colour
      if (UseBackgndColorVecplot) then
         call plot_sci(0)
      else
         call plot_sci(1)
      endif
      usevecplot = .false.
      if (irotate) then
         if (allocated(vecplot)) usevecplot = .true.
         if (debugmode) print*,'DEBUG: using vecplot' ! this is to indicate (to me) that extra memory is in use
      endif
      !
      !--interpolate using appropriate routine for number of dimensions
      !
      select case(ndim)
      case(3)
         if (x_sec) then ! take vector plot in cross section
            if (got_h) then
               if (usevecplot) then ! using rotation
                  call interpolate3D_xsec_vec(xplot(1:ninterp), &
                    yplot(1:ninterp),zplot(1:ninterp), &
                    hh(1:ninterp),weight(1:ninterp), &
                    vecplot(1,1:ninterp),vecplot(2,1:ninterp), &
                    icolourme(1:ninterp),ninterp,xmin,ymin,zslicepos, &
                    vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,pixwidthvecy,inormalise)
               else
                  call interpolate3D_xsec_vec(xplot(1:ninterp), &
                    yplot(1:ninterp),zplot(1:ninterp), &
                    hh(1:ninterp),weight(1:ninterp), &
                    dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                    icolourme(1:ninterp),ninterp,xmin,ymin,zslicepos, &
                    vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,pixwidthvecy,inormalise)
               endif
            else
               ! don't have smoothing length, use averaging
               if (usevecplot) then
                  call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
                       vecplot(1,1:ninterp),vecplot(2,1:ninterp),icolourme(1:ninterp), &
                       xmin,ymin,pixwidthvec,pixwidthvecy,vecpixx,vecpixy, &
                       ninterp,numpixx,numpixy,zplot(1:ninterp),zslicemin,zslicemax)
               else
                  call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
                       dat(1:ninterp,ivecx),dat(1:ninterp,ivecy),icolourme(1:ninterp), &
                       xmin,ymin,pixwidthvec,pixwidthvecy,vecpixx,vecpixy, &
                       ninterp,numpixx,numpixy,zplot(1:ninterp),zslicemin,zslicemax)
               endif            
            endif
         else
            if (iplotsynchrotron .and. .not.iplotstreamlines .and. .not.iplotarrowheads) then
               !--get synchrotron polarisation vectors
               if (iutherm.gt.0 .and. iutherm.le.numplot .and. uthermcutoff.gt.0.) then
                  if (usevecplot) then
                     call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),vecplot(1,1:ninterp),vecplot(2,1:ninterp), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,pixwidthvec, &
                       rcrit,zcrit,synchrotronspecindex,pixwidthvec,.false., &
                       dat(1:ninterp,iutherm),uthermcutoff)
                  else
                     call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,pixwidthvec, &
                       rcrit,zcrit,synchrotronspecindex,pixwidthvec,.false., &
                       dat(1:ninterp,iutherm),uthermcutoff)
                  endif
               else
                  if (usevecplot) then
                     call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),vecplot(1,1:ninterp),vecplot(2,1:ninterp), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,pixwidthvec, &
                       rcrit,zcrit,synchrotronspecindex,pixwidthvec,.false.)
                  elseif (.not.iplotstreamlines) then
                     call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,pixwidthvec, &
                       rcrit,zcrit,synchrotronspecindex,pixwidthvec,.false.)
                  endif
               endif
            elseif (.not.(iplotstreamlines .and. use3Dstreamlines)) then
               if (got_h) then
                  if (usevecplot) then
                     if (.not.allocated(vecplot)) stop 'internal error: vecplot not allocated'
                     call interpolate3D_proj_vec(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),vecplot(1,1:ninterp),vecplot(2,1:ninterp), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,pixwidthvecy,&
                       .false.,zobservertemp,dzscreentemp)
                  else
                     call interpolate3D_proj_vec(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,pixwidthvecy, &
                       .false.,zobservertemp,dzscreentemp)
                  endif
               else
               ! don't have smoothing length, use averaging
                  if (usevecplot) then
                     call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
                          vecplot(1,1:ninterp),vecplot(2,1:ninterp),icolourme(1:ninterp), &
                          xmin,ymin,pixwidthvec,pixwidthvecy,vecpixx,vecpixy, &
                          ninterp,numpixx,numpixy)
                  else
                     call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
                          dat(1:ninterp,ivecx),dat(1:ninterp,ivecy),icolourme(1:ninterp), &
                          xmin,ymin,pixwidthvec,pixwidthvecy,vecpixx,vecpixy, &
                          ninterp,numpixx,numpixy)
                  endif
               endif
            endif
            !--adjust the units of the z-integrated quantity
            !if (iRescale .and. units(ih).gt.0.) then
            !   vecpixx = vecpixx*(unitzintegration/units(ih))
            !   vecpixy = vecpixy*(unitzintegration/units(ih))
            !endif
         endif
      case(2)
         !
         !--two dimensions
         !
         if (got_h) then
            if (usevecplot) then
               call interpolate2D_vec(xplot(1:ninterp),yplot(1:ninterp), &
                 hh(1:ninterp),weight(1:ninterp),vecplot(1,1:ninterp), &
                 vecplot(2,1:ninterp),icolourme(1:ninterp),ninterp,xmin,ymin, &
                 vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,pixwidthvecy,inormalise,&
                 isperiodicx,isperiodicy)
            else
               call interpolate2D_vec(xplot(1:ninterp),yplot(1:ninterp), &
                 hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,ivecx), &
                 dat(1:ninterp,ivecy),icolourme(1:ninterp),ninterp,xmin,ymin, &
                 vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,pixwidthvecy,inormalise,&
                 isperiodicx,isperiodicy)
            endif
         else
            ! don't have smoothing length, use averaging
            if (usevecplot) then
               call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
                    vecplot(1,1:ninterp),vecplot(2,1:ninterp),icolourme(1:ninterp), &
                    xmin,ymin,pixwidthvec,pixwidthvecy,vecpixx,vecpixy, &
                    ninterp,numpixx,numpixy)
            else
               call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
                    dat(1:ninterp,ivecx),dat(1:ninterp,ivecy),icolourme(1:ninterp), &
                    xmin,ymin,pixwidthvec,pixwidthvecy,vecpixx,vecpixy, &
                    ninterp,numpixx,numpixy)
            endif         
         endif

      case default
         print "(a,i1,a)",'ERROR: Cannot do vector plotting in ',ndim,' dimensions'
         return
      end select
      !
      !--plot it, either as streamlines or arrows
      !
      if (iplotstreamlines) then
         if (ndim.eq.3) then
            !--normalise the 3D vector field
            do j=1,numpixy
               do i=1,numpixx
                  vmag = sqrt(vecpixx(i,j)**2 + vecpixy(i,j)**2)
                  if (vmag.gt.tiny(vmag)) then
                     vecpixx(i,j) = vecpixx(i,j)/vmag
                     vecpixy(i,j) = vecpixy(i,j)/vmag
                  endif
              enddo
            enddo
         endif

         if (ndim.eq.3 .and. use3Dstreamlines .and. .not.x_sec) then
            if (usevecplot) then
               if (.not.allocated(vecplot)) stop 'vecplot not allocated'
               call vecplot3D_proj(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp), &
                       vecplot(1,1:ninterp),vecplot(2,1:ninterp),vecplot(3,1:ninterp),vmax, &
                       weight(1:ninterp),icolourme(1:ninterp),ninterp,pixwidthvec,zobservertemp,dzscreentemp)            
            else
               ivecz = ivecx + (iplotz - ix(1))
               call vecplot3D_proj(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp), &
                       dat(1:ninterp,ivecx),dat(1:ninterp,ivecy),dat(1:ninterp,ivecz),vmax, &
                       weight(1:ninterp),icolourme(1:ninterp),ninterp,pixwidthvec,zobservertemp,dzscreentemp)

            endif
         else
            call streamlines(vecpixx,vecpixy,datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,pixwidthvec)

            if (ihidearrowswherenoparts) then
               datmax = maxval(datpixvec(1:numpixx,1:numpixy))
               blankval = 2.*datmax
               call mask_vectors(xplot(1:ninterp),yplot(1:ninterp),icolourme(1:ninterp),ninterp, &
                                 xmin,xmax,ymin,ymax,datpixvec(1:numpixx,1:numpixy), &
                                 datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,minpartforarrow,blankval)

               !--use blanking for values of zero
               call render_pix(datpixvec(1:numpixx,1:numpixy), &
                            minval(datpixvec(1:numpixx,1:numpixy)), &
                            datmax, &
                            'crap',numpixx,numpixy,xmin,ymin,pixwidthvec,pixwidthvecy,    &
                            0,.true.,0,ncontours,.false.,ilabelcont,blank=blankval)
            else
               call render_pix(datpixvec(1:numpixx,1:numpixy), &
                            minval(datpixvec(1:numpixx,1:numpixy)), &
                            maxval(datpixvec(1:numpixx,1:numpixy)), &
                            'crap',numpixx,numpixy,xmin,ymin,pixwidthvec,pixwidthvecy,    &
                            0,.true.,0,ncontours,.false.,ilabelcont)
            endif
         endif

      else
         if (ihidearrowswherenoparts) then
            call mask_vectors(xplot(1:ninterp),yplot(1:ninterp),icolourme(1:ninterp),ninterp, &
                              xmin,xmax,ymin,ymax,vecpixx,vecpixy,numpixx,numpixy,minpartforarrow,0.)
         endif

         plotlegend = iVecplotLegend .and. ipanelselect(iVecLegendOnPanel,ipanel,irow,icolumn)
         call render_vec(vecpixx,vecpixy,vmax, &
              numpixx,numpixy,xmin,ymin,pixwidthvec,pixwidthvecy,trim(label),' ',plotlegend)

         if (iplotsynchrotron .and. .not. iplotarrowheads) then
            !--get synchrotron polarisation intensity using more pixels
            if (iutherm.gt.0 .and. iutherm.le.numplot .and. uthermcutoff.gt.0.) then
               call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                 yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                 weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                 icolourme(1:ninterp),ninterp,xmin,ymin, &
                 datpixvec(1:npixx,1:npixy),datpixvec(1:npixx,1:npixy), & ! these are just dummy arguments
                 datpixvec(1:npixx,1:npixy),npixx,npixy,pixwidth, &
                 rcrit,zcrit,synchrotronspecindex,pixwidthvec,.true., &
                 dat(1:ninterp,iutherm),uthermcutoff)
            else
               call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                 yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                 weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                 icolourme(1:ninterp),ninterp,xmin,ymin, &
                 datpixvec(1:npixx,1:npixy),datpixvec(1:npixx,1:npixy), & ! these are just dummy arguments
                 datpixvec(1:npixx,1:npixy),npixx,npixy,pixwidth, &
                 rcrit,zcrit,synchrotronspecindex,pixwidthvec,.true.)
            endif

            !--adjust the units of the z-integrated quantity
            !if (iRescale .and. units(ih).gt.0. .and..not.inormalise) then
            !   datpix = datpix*(unitzintegration/units(ih))
            !endif

            !--plot contours of synchrotron intensity
            call render_pix(datpixvec(1:npixx,1:npixy),minval(datpixvec(1:npixx,1:npixy)), &
              maxval(datpixvec(1:npixx,1:npixy)),'crap', &
              npixx,npixy,xmin,ymin,pixwidth,pixwidthy,0,.true.,0,ncontours,.false.,ilabelcont)
         endif

      endif

   endif

   !--restore colour index and line width
   call plot_sci(icoloursav)
   call plot_slw(linewidthprev)

  end subroutine vector_plot

end subroutine plotstep

!----------------------------------------------------
! adapt the (particle plot) limits to include all
! particles which are to be plotted on the page
!----------------------------------------------------
subroutine adapt_limits(iplot,xploti,xmini,xmaxi,xminadaptive,xmaxadaptive,labeli,&
                        iamtype,ntoti,npartoftype,iusetype,ipagechange)
  use params,          only:int1,maxparttypes
  use labels,          only:is_coord
  use limits,          only:assert_sensible_limits
  use settings_limits, only:scalemax,iadapt,iadaptcoords
  use settings_data,   only:debugmode,ndim
  use settings_part,   only:iplotline
  implicit none
  integer,            intent(in)    :: iplot
  real, dimension(:), intent(in)    :: xploti
  real,               intent(inout) :: xmini,xmaxi,xminadaptive,xmaxadaptive
  character(len=*),   intent(in)    :: labeli
  integer(kind=int1), dimension(:), intent(in) :: iamtype
  integer, intent(in)               :: ntoti
  integer, dimension(:), intent(in) :: npartoftype
  logical, dimension(:), intent(in) :: iusetype
  logical,               intent(in) :: ipagechange
  integer :: index1,index2,itype,i
  logical :: mixedtypes

  !--calculate adaptive limits for this quantity
  xminadaptive = huge(xminadaptive)
  xmaxadaptive = -huge(xmaxadaptive)

  mixedtypes = size(iamtype).gt.1
  if (mixedtypes) then
     do i=1,ntoti
        itype = iamtype(i)
        if (iusetype(itype) .or. (iplotline.and.itype.eq.1)) then
           xminadaptive = min(xminadaptive,xploti(i))
           xmaxadaptive = max(xmaxadaptive,xploti(i))*scalemax
        endif
     enddo
  else
     index1 = 1
     do itype=1,maxparttypes
        index2 = index1 + npartoftype(itype) - 1
        if (iusetype(itype).and.npartoftype(itype).gt.0 &
           .or. (iplotline.and.itype.eq.1)) then
           xminadaptive = min(xminadaptive,minval(xploti(index1:index2)))
           xmaxadaptive = max(xmaxadaptive,maxval(xploti(index1:index2))*scalemax)
        endif
        index1 = index2 + 1
     enddo
  endif
  !--avoid infs and NaNs
  call assert_sensible_limits(xminadaptive,xmaxadaptive)

  if (debugmode) print*,'DEBUG: ',iplot,': '//trim(labeli)// &
                 'min,max adaptive = ',xminadaptive,xmaxadaptive

  !--set these as limits if adaptive limits are on
  if (.not.interactivereplot) then
     if (((is_coord(iplot,ndim) .and. iadaptcoords) &
     .or.(.not.is_coord(iplot,ndim) .and. iadapt)) .and. ipagechange) then
        print "(1x,a)",'adapting '//trim(labeli)//' limits'
        xmini = xminadaptive
        xmaxi = xmaxadaptive
     endif
  endif

end subroutine adapt_limits

!-------------------------------------------------------------------
! interface to log, inverse transformations:
! also adjusts label (depending on
! whether log axes are also set or not).
!  (independent)
!-------------------------------------------------------------------
subroutine applytrans(xploti,xmini,xmaxi,labelxi,itransxi,chaxis,iplotxi,iaxis,intreplot)
  use transforms,    only:transform,transform_label,transform_limits
  use settings_data, only:numplot
  implicit none
  integer, intent(in)               :: itransxi,iplotxi,iaxis
  real, dimension(:), intent(inout) :: xploti
  real, intent(inout)               :: xmini,xmaxi
  character(len=*), intent(inout)   :: labelxi
  character(len=1), intent(in)      :: chaxis
  logical, intent(in)               :: intreplot
  integer           :: itranstemp,lstr
  character(len=20) :: string

  if (itransxi.ne.0) then
     if (iplotxi.gt.0 .and. iplotxi.le.numplot) call transform(xploti(:),itransxi)
     if ((chaxis.eq.'x' .and. (iaxis.eq.10 .or. iaxis.eq.30)).or. &
         (chaxis.eq.'y' .and. (iaxis.eq.20 .or. iaxis.eq.30))) then ! logarithmic axes
        write(string,*) itransxi
        string = adjustl(string)
        itranstemp = 0
        lstr = len_trim(string)
        if (string(lstr:lstr).eq.'1') then
           if (lstr.gt.1) read(string(1:lstr-1),*) itranstemp
           labelxi = transform_label(labelxi,itranstemp)
        else
           labelxi = transform_label(labelxi,itransxi)
        endif
     else
        labelxi = transform_label(labelxi,itransxi)
     endif
     if (.not.intreplot) call transform_limits(xmini,xmaxi,itransxi)
  endif
end subroutine applytrans

!-------------------------------------------------------------------
! interface for adding rotation and perspective
! (completely independent)
!-------------------------------------------------------------------
subroutine rotationandperspective(anglexi,angleyi,anglezi,dzscreen,zobs,xploti,yploti,zploti, &
                                  ntot,iplotx,iploty,iplotz,dat,ivecstart,vecploti,itrackpart)
  use labels,           only:ix
  use settings_data,    only:ndim,xorigin,debugmode
  use settings_xsecrot, only:use3Dperspective
  use rotation,         only:rotate2D,rotate3D
  implicit none
  real,                 intent(in)  :: anglexi,angleyi,anglezi,dzscreen,zobs
  real, dimension(:), intent(inout) :: xploti,yploti,zploti
  real, dimension(:,:), intent(in)  :: dat
  real, dimension(:,:), intent(out) :: vecploti
  integer,              intent(in)  :: ntot,iplotx,iploty,iplotz,ivecstart,itrackpart
  integer               :: j,iposx,iposy,iposz,i
  real                  :: angleradx,anglerady,angleradz
  real, dimension(ndim) :: xcoords,veci
  !
  !--convert angles to radians
  !
  angleradz = anglezi*pi/180.
  anglerady = angleyi*pi/180.
  angleradx = anglexi*pi/180.
  if (ndim.eq.3) then
     print "(1x,a,2(f6.2,1x),f6.2,a)",'rotation: (z, y, x) = (',anglezi,angleyi,anglexi,')'
  else
     print "(1x,a,f6.2)",'rotating particles about z by ',anglezi
  endif
  if (ndim.eq.3 .and. use3Dperspective) then
     print*,' observer height = ',zobs,', screen at ',zobs-dzscreen
  elseif (ndim.eq.3) then
     if (abs(zobs).gt.tiny(zobs) .or. abs(dzscreen).gt.tiny(dzscreen)) then
        print "(a)",' INTERNAL ERROR: no 3D perspective but observer set'
     endif
  endif
  if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
     print*,'rotating about tracked particle ',itrackpart,' x,y,z = ',dat(itrackpart,ix(1:ndim))
  elseif (any(abs(xorigin).ge.tiny(xorigin))) then
     print*,'rotating about x,y,z = ',xorigin(1:ndim)
  endif

  if (debugmode .and. ivecstart.gt.0) print "(1x,a)",'(also rotating vector components)'
  !
  !--set location of x,y and z
  !  such that:
  !  ix(iposx) = iplotx
  !  ix(iposy) = iploty
  !  ix(iposz) = iplotz
  !
  iposx = 1  ! this is "just in case"
  iposy = 2
  iposz = 3
  do i=1,ndim
     if (ix(i).eq.iplotx) then
        iposx = i
     elseif (ix(i).eq.iploty) then
        iposy = i
     elseif (ix(i).eq.iplotz) then
        iposz = i
     else
        print "(a)",' WARNING: internal error in ix setting for rotation: ix = ',ix(:)
     endif
  enddo
  if (debugmode) print*,'DEBUG: in rotation, iplotz = ',iplotz,' iposz = ',iposz, xorigin(:)

!$omp parallel default(none) &
!$omp shared(dat,xorigin,ndim,angleradx,anglerady,angleradz,zobs,dzscreen) &
!$omp shared(xploti,yploti,zploti,iposx,iposy,iposz,iplotz,ntot,ix,itrackpart) &
!$omp shared(vecploti,ivecstart) &
!$omp private(j,xcoords,veci)
!$omp do
  do j=1,ntot
     if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
        xcoords(1:ndim) = dat(j,ix(1:ndim)) - dat(itrackpart,ix(1:ndim))
     else
        xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
     endif
     if (ndim.eq.2) then
        call rotate2D(xcoords(:),angleradz)
     elseif (ndim.eq.3) then
        call rotate3D(xcoords(1:ndim),angleradx,anglerady,angleradz,zobs,dzscreen)
     endif
     if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
        xploti(j) = xcoords(iposx) + dat(itrackpart,ix(iposx))
        yploti(j) = xcoords(iposy) + dat(itrackpart,ix(iposy))
        if (iplotz.gt.0) then
           zploti(j) = xcoords(iposz) + dat(itrackpart,ix(iposz))
        endif
     else
        xploti(j) = xcoords(iposx) + xorigin(iposx)
        yploti(j) = xcoords(iposy) + xorigin(iposy)
        if (iplotz.gt.0) then
           zploti(j) = xcoords(iposz) + xorigin(iposz)
        endif
     endif
!
!--rotate vector components
!
     if (ivecstart.gt.0) then
        veci(1:ndim) = dat(j,ivecstart:ivecstart+ndim-1)
        if (ndim.eq.2) then
           call rotate2D(veci(:),angleradz)
        elseif (ndim.eq.3) then
           call rotate3D(veci(1:ndim),angleradx,anglerady,angleradz,zobs,dzscreen)
        endif
        vecploti(1,j) = veci(iposx)
        vecploti(2,j) = veci(iposy)
        if (ndim.ge.3) vecploti(3,j) = veci(iposz)
     endif
  enddo
!$omp end do
!$omp end parallel

  return
end subroutine rotationandperspective

!-------------------------------------------------------------------
! interface for plotting rotated axes
!-------------------------------------------------------------------
subroutine rotatedaxes(irotateaxes,iplotx,iploty,anglexi,angleyi,anglezi,dzscreen,zobs)
  use labels,   only:ix
  use rotation, only:rotate_axes3D,rotate_axes2D
  use settings_data, only:ndim,xorigin
  use settings_xsecrot, only:xminrotaxes,xmaxrotaxes,use3Dperspective
  implicit none
  integer, intent(in) :: irotateaxes,iplotx,iploty
  real, intent(in) :: anglexi,angleyi,anglezi
  real, intent(inout) :: dzscreen,zobs
  real :: angleradx,anglerady,angleradz
  !
  !--convert angles to radians
  !
  angleradz = anglezi*pi/180.
  anglerady = angleyi*pi/180.
  angleradx = anglexi*pi/180.

  if (ndim.eq.3) then
     if (.not.use3Dperspective .and. dzscreen.gt.tiny(zobs)) then
        print "(a)",' INTERNAL ERROR: no 3D perspective but observer set'
        zobs = 0.
        dzscreen = 0.
     endif
     call rotate_axes3D(irotateaxes,iplotx-ix(1)+1,iploty-ix(1)+1, &
          xminrotaxes(1:ndim),xmaxrotaxes(1:ndim),xorigin(1:ndim), &
          angleradx,anglerady,angleradz,zobs,dzscreen)
  elseif (ndim.eq.2) then
     call rotate_axes2D(irotateaxes,xminrotaxes(1:ndim), &
                       xmaxrotaxes(1:ndim),xorigin(1:ndim),angleradz)
  endif

  return
end subroutine rotatedaxes

end module timestep_plotting
