!------------------------------------------------------------------------
!
! This file is part of the SPLASH visualisation package for SPH.
! (c) 2002-2008 Daniel Price
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
  integer, private :: iplotx,iploty,iplotz,irenderplot,icontourplot
  integer, private :: ivectorplot,ivecx,ivecy
  integer, private :: nyplots,npartdim,nyplotfirstonpage,ifirststeponpage    
  integer, private :: ngrid,nframefirstonpage
  integer, private :: just, ntitles,nsteplegendlines
  integer, private :: iplots,ipanel
  integer, private :: iframesave

  real, dimension(:), allocatable, private :: datpix1D, xgrid
  real, dimension(:,:), allocatable, private :: datpix,datpixcont,brightness
  real, dimension(:,:,:), allocatable, private :: datpix3D,datpixcont3D
  real, private :: xmin,xmax,ymin,ymax,zmin
  real, private :: rendermin,rendermax,vecmax,contmin,contmax
  real, private :: dz,zslicepos,zobservertemp,dzscreentemp,taupartdepthtemp,rkappafac
  real, private :: dxgrid,xmingrid,xmaxgrid
  real, private :: angletempx, angletempy, angletempz
  !--buffer for interactive mode on multiplots
  integer, dimension(maxplot) :: iplotxtemp,iplotytemp,irendertemp,ivecplottemp
  real, dimension(maxplot) :: xminmulti,xmaxmulti,xminadapt,xmaxadapt
  real, dimension(maxplot) :: vptxmin,vptxmax,vptymin,vptymax,barwmulti
  real, private :: xminadapti,xmaxadapti,yminadapti,ymaxadapti,renderminadapt,rendermaxadapt
  real, private :: contminadapt,contmaxadapt
  real, parameter, private :: pi = 3.1415926536

  logical, private :: iplotpart,iplotcont,x_sec,isamexaxis,isameyaxis
  logical, private :: inewpage, tile_plots, lastplot
  logical, private :: imulti,irerender,iAllowspaceforcolourbar
  logical, private :: interactivereplot,ihavesetcolours,vectordevice
  
  public :: initialise_plotting, plotstep
  private

contains

!
! initialise plotting options
! called once for all steps
!
subroutine initialise_plotting(ipicky,ipickx,irender_nomulti,icontour_nomulti,ivecplot)
  use params
  use colours, only:colour_set
  use labels, only:label,ipowerspec,ih,ipmass,irho,iamvec,isurfdens,itoomre,iutherm,ipdf
  use limits, only:lim,rangeset
  use multiplot, only:multiplotx,multiploty,irendermulti,icontourmulti, &
                 nyplotmulti,x_secmulti,ivecplotmulti
  use prompting, only:prompt
  use titles, only:read_titles,read_steplegend
  use settings_data, only:ndim,ndimV,numplot,ncolumns,ncalc,ndataplots,required,icoords,icoordsnew
  use settings_page, only:nacross,ndown,ipapersize,tile,papersizex,aspectratio,&
      colour_fore,colour_back,iadapt,iadaptcoords,linewidth,device,nomenu,interactive
  use settings_part, only:linecolourthisstep,linecolour,linestylethisstep,linestyle,iexact
  use settings_render, only:icolours,iplotcont_nomulti,iColourBarStyle,icolour_particles
  use settings_xsecrot, only:xsec_nomulti,xsecpos_nomulti,flythru,nxsec, &
                        xseclineX1,xseclineX2,xseclineY1,xseclineY2,xsecwidth, &
                        use3Dperspective,use3Dopacityrendering,zobserver,dzscreenfromobserver,taupartdepth
  use settings_powerspec, only:options_powerspec
  use particle_data, only:npartoftype
  use projections3D, only:coltable
  use pdfs, only:options_pdf
  implicit none
  real, parameter :: pi=3.1415926536
  integer, intent(in) :: ipicky,ipickx,irender_nomulti,icontour_nomulti,ivecplot
  integer :: i,j,ierr,ifirst,iplotzprev,ilen,pgopen
  logical :: iadapting,iamrendering,icoordplot,iallrendered,ians
  real :: hav,pmassav,dzsuggest
  character(len=1) :: char
  character(len=20) :: string,devstring
  
  !------------------------------------------------------------------------
  ! initialisations
  !------------------------------------------------------------------------

  isamexaxis = .true.  ! same x axis on all plots? (only relevant for >1 plots per page)
  isameyaxis = .true.  ! same y axis on all plots?
  tile_plots = .false.
  iplots = 0 ! counter for how many plots have been plotted in total
  ipanel = 0  ! counter for which panel we are in on plotting page
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
  
  xmin = 0.
  xmax = 0.
  ymin = 0.
  ymax = 0.
  xminadapt = huge(xminadapt)
  xmaxadapt = -huge(xmaxadapt)

  if (ndim.eq.1) x_sec = .false. ! can't have xsec in 1D
  nxsec = 1

  iamrendering = .false.
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
     if (any(irendermulti(1:nyplotmulti).gt.ndim)) iamrendering = .true.
     if (any(x_secmulti(1:nyplotmulti))) x_sec = .true.
  else
     !
     !--or else set number of plots = 1 and use ipicky and ipickx
     !
     imulti = .false.
     nyplots = 1 
     iploty = ipicky
     iplotx = ipickx
     if (irender_nomulti.gt.ndim) iamrendering = .true.
  endif

  !------------------------------------------------------------------------
  ! initialise options to be set before plotting

  icoordplot = iploty.le.ndim .and. iplotx.le.ndim
  iallrendered = iamrendering
  iplotz = 0
  if (imulti) then
     do i=1,nyplotmulti
        if (multiplotx(i).le.ndim .and. multiploty(i).le.ndim) then
           icoordplot = .true.
           !--this check is to see if any co-ordinate plots involve just particles
           !  (if so need to initialise the cross section slice width)
           if (irendermulti(i).le.ndim) iallrendered = .false.
           iplotzprev = iplotz
           !!--work out coordinate that is not being plotted on cross-section/ 3D plots
           iplotz = 0
           if (ndim.ge.3 .and. (x_sec .or. use3Dperspective)) then
              do j=1,ndim
                 if ((multiplotx(i).ne.multiploty(i)).and. &
                     (j.ne.multiplotx(i)).and.(j.ne.multiploty(i))) iplotz = j
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
               (j.ne.iplotx).and.(j.ne.iploty)) iplotz = j
        enddo
     endif
  endif

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
  if (iamrendering .and. icolours.ne.0 .and. iColourbarStyle.gt.0) then
     !--this option means that a margin is set aside for a colour bar on tiled plots
     iAllowspaceforcolourbar = .true.
     !--do not allow tiled plots if multiple (different) colour bars are plotted
     if (tile_plots) then
        ifirst = 0
        do i=1,nyplots
           if (irendermulti(i).gt.ndim .and. ifirst.eq.0) ifirst = i
           if (ifirst.gt.0) then
              if (irendermulti(i).gt.ndim .and. irendermulti(i).ne.irendermulti(ifirst)) then
                 if (tile_plots) print "(a)",'WARNING: cannot tile plots because of multiple colour bars'
                 tile_plots = .false.
              endif
           endif
        enddo
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
          call pgbegin(0,'/xw',1,1)
          call pgenv(lim(1,1),lim(1,2),lim(2,1),lim(2,2),1,0)
          call pgcurs(xseclineX1,xseclineY1,char)
          print*,'please select cross section line'
          call pgband(1,1,xseclineX1,xseclineY1,xseclineX2,xseclineY2,char)
          print*,'cross section line: xmin = ',xseclineX1,' xmax = ',xseclineX2
          print*,'                    ymin = ',xseclineY1,' ymax = ',xseclineY2
          call pgend       
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
       if (abs(dzscreenfromobserver).lt.tiny(dzscreenfromobserver)) dzscreenfromobserver = lim(iplotz,2)
       call prompt('enter z coordinate of observer ',zobserver)
       call prompt('enter distance for unit magnification ',dzscreenfromobserver,0.)
!
!--initialise opacity for 3D opacity rendering
!       
       if (use3Dopacityrendering .and. iamrendering) then
          hav = lim(ih,1) !! 0.5*(lim(ih,2) + lim(ih,1))
          if (hav.le.epsilon(hav)) hav = 0.5*lim(ih,2) ! take 0.5*max if min is zero
          pmassav = lim(ipmass,1)
          if (pmassav.le.epsilon(hav)) pmassav = 0.5*lim(ipmass,2) ! take 0.5*max if min is zero
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

  !!--prompt for options if plotting power spectrum      
  if (iploty.eq.ipdf .and. .not. imulti &
     .or. (imulti.and.any(multiploty(1:nyplotmulti).eq.ipdf))) then
     call options_pdf
  endif

  !!--for fast data read, set which columns are required from the file
  !   (note that required(0)= whatever is a valid statement, just has no effect)
  required = .false.  ! by default, no columns required
!  if (fastdataread) then
     if (imulti) then
        required(multiplotx(1:nyplotmulti)) = .true.
        required(multiploty(1:nyplotmulti)) = .true.
        required(irendermulti(1:nyplotmulti)) = .true.
        required(icontourmulti(1:nyplotmulti)) = .true.
     else
        required(iplotx) = .true.
        required(iploty) = .true.
     endif
     required(iplotz) = .true.
     if (iamrendering) then
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
     if (iploty.eq.itoomre) required(iutherm) = .true.
  !!--must read everything if we are plotting a calculated quantity
     if (any(required(ncolumns+1:ncolumns+ncalc))) required = .true.
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
        required(1:ndim) = .true.
        if (iplotx.gt.0 .and. iplotx.le.numplot) then
           if (iamvec(iplotx).gt.0) required(iamvec(iplotx):iamvec(iplotx)+ndimV-1) = .true.
        endif
        if (iploty.gt.0 .and. iploty.le.numplot) then
           if (iamvec(iploty).gt.0) required(iamvec(iploty):iamvec(iploty)+ndimV-1) = .true.
        endif
     endif
!  endif

  !!--read step titles (don't need to store ntitles for this)
  nsteplegendlines = 0
  call read_steplegend(nsteplegendlines)
  !!--read plot titles
  ntitles = 0
  call read_titles(ntitles)

  !!------------------------------------------------------------------------
  ! initialise PGPLOT

  !!--start PGPLOT (prompt for device)
  if (len_trim(device).eq.0) then
     call pgbegin(0,'?',1,1)
  else
     ierr = pgopen(trim(device))
     if (ierr.le.0) then  ! zero or negative indicates an error
        print "(a)",' ERROR: unknown PGPLOT device "'//trim(device)//'"'
        stop
     endif
  endif
  
  !--query whether or not device is interactive
  call pgqinf('CURSOR',string,ilen)
  
  if (string(1:ilen).eq.'YES') then
     !--turn menu and interactive mode on if
     !  interactive device invoked from the command line
     if (nomenu) then
        interactive = .true.
        nomenu=.false.
     endif
     !--smoothing length is required if interactive device and coordinate plot
     if (icoordplot) required(ih) = .true.
  endif

  !!--set paper size if necessary
  if (ipapersize.gt.0 .and. papersizex.gt.0.0 .and. aspectratio.gt.0.0 ) then
     call pgpaper(papersizex,aspectratio)
  endif
  !!--turn off page prompting
  call pgask(.false.)
  !!if (.not. interactive) call pgbbuf !! start buffering output
  
  !!--set background/foreground colours
  ierr = 0
  if (len_trim(colour_back).gt.0) call pgscrn(0,colour_back,ierr)
  if (ierr /= 0) print 10,'background',trim(colour_back)
  if (len_trim(colour_fore).gt.0) call pgscrn(1,colour_fore,ierr)
  if (ierr /= 0) print 10,'foreground',trim(colour_fore)
10 format(' error: ',a,' colour "',a,'" not found in table')
  
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
  call pgqinf('TYPE',devstring,ilen)
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
        call pgslw(2)
     else
        print "(a)",' setting line width = 1 for '//devstring(1:ilen)//' device'
        call pgslw(1)
     endif  
  else
     call pgslw(linewidth)
  endif

  linecolourthisstep = linecolour
  linestylethisstep = linestyle

end subroutine initialise_plotting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plotstep(ipos,istep,istepsonpage,irender_nomulti,icontour_nomulti,ivecplot, &
                    iamtype,npartoftype,masstype,dat,timei,gammai,ipagechange,iadvance)
  use params
  use colours, only:colour_set
  use filenames, only:nsteps
  use exact, only:exact_solution,atstar,ctstar,sigma
  use toystar1D, only:exact_toystar_ACplane
  use toystar2D, only:exact_toystar_ACplane2D
  use labels, only:label,labelvec,iamvec,lenlabel, &
              ih,irho,ipmass,ix,iacplane,ipowerspec,isurfdens,itoomre,iutherm,ipdf
  use limits, only:lim,get_particle_subset,lim2,lim2set
  use multiplot,only:multiplotx,multiploty,irendermulti,ivecplotmulti,itrans, &
                icontourmulti,x_secmulti,xsecposmulti
  use particle_data, only:maxpart,icolourme
  use settings_data, only:numplot,ndataplots,icoords,icoordsnew,ndim,ndimV,nfreq,iRescale, &
                     iendatstep,ntypes,UseTypeInRenderings,itrackpart,required,ipartialread,xorigin
  use settings_limits, only:iadapt,iadaptcoords,scalemax
  use settings_part, only:iexact,iplotpartoftype,imarktype,PlotOnRenderings,UseTypeInContours, &
                     iplotline,linecolourthisstep,linestylethisstep,ifastparticleplot
  use settings_page, only:nacross,ndown,iadapt,interactive,iaxis,usesquarexy, &
                     charheight,iPlotTitles,vpostitle,hpostitle,fjusttitle,nstepsperpage
  use settings_render, only:npix,ncontours,icolours, &
      iColourBarStyle,icolour_particles,inormalise_interpolations,ifastrender,ilabelcont
  use settings_vecplot, only:npixvec, iplotpartvec
  use settings_xsecrot, only:nxsec,irotateaxes,xsec_nomulti,irotate,flythru,use3Dperspective, &
      use3Dopacityrendering,writeppm,anglex,angley,anglez,zobserver,dzscreenfromobserver,taupartdepth, &
      xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2,nseq,nframes,getsequencepos,insidesequence
  use settings_powerspec, only:nfreqspec,wavelengthmin,wavelengthmax,ipowerspecx,ipowerspecy,idisordered
  use settings_units, only:units,unitslabel,unitzintegration
!
!--subroutines called from this routine
!
  use colourparts
  use transforms, only:transform,transform_limits,transform_label,transform_inverse,islogged
  use interactive_routines
  use particleplots, only:particleplot
  use powerspectrums, only:powerspectrum,powerspec3D_sph
  use interpolations1D, only:interpolate1D
  use interpolations2D, only:interpolate2D, interpolate2D_xsec
  use interpolations3D, only:interpolate3D
  use projections3D, only:interpolate3D_projection
  use interpolate3D_opacity, only:interp3D_proj_opacity !,interp3D_proj_opacity_writeppm
  use xsections3D, only:interpolate3D_fastxsec,interpolate3D_xsec_vec
  use render, only:render_pix
  use pagesetup, only:redraw_axes
  use disc, only:disccalc,discplot
  use exactfromfile, only:exact_fromfile
  use write_pixmap, only:iwritepixmap,writepixmap,write_pixmap_ppm
  use pdfs, only:pdfcalc,pdfplot,npdfbins

  implicit none
  integer, intent(inout) :: ipos, istepsonpage
  integer, intent(in) :: istep,irender_nomulti,icontour_nomulti,ivecplot
  integer(kind=int1), dimension(:), intent(in) :: iamtype
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(maxparttypes), intent(in) :: masstype
  real, dimension(:,:), intent(in) :: dat
  real, intent(in) :: timei,gammai
  logical, intent(in) :: ipagechange
  integer, intent(inout) :: iadvance
  
  integer :: ntoti,iz,iseqpos
  integer :: i,j,k,icolumn,irow
  integer :: nyplot,nframesloop
  integer :: irender,irenderpart
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: npixyvec,nfreqpts
  integer :: icolourprev,linestyleprev
  integer :: ierr,ipt,nplots,nyplotstart,iaxisy,iaxistemp
  integer :: ivectemp,iamvecx,iamvecy,itransx,itransy,itemp
  integer :: iframe

  real, parameter :: tol = 1.e-10 ! used to compare real numbers
  real, dimension(max(maxpart,2000)) :: xplot,yplot,zplot
  real, dimension(maxpart) :: hh,weight
  real, dimension(:), allocatable :: renderplot
  real, dimension(:,:), allocatable :: vecplot
  real :: rkappa
  real :: zslicemin,zslicemax,dummy,pmassmin,pmassmax
  real :: pixwidth,pixwidthvec,dxfreq

  character(len=lenlabel+20) :: labelx,labely,labelz,labelrender,labelvecplot,labelcont
  character(len=20) :: labeltimeunits
  
  logical :: iPlotColourBar, rendering, inormalise, logged, loggedcont
  logical :: dumxsec, isetrenderlimits, gotcontours
  logical :: ichangesize, initx, inity, isameweights
  
34   format (25(' -'))

  !--set labels to blank (just in case)
  labelx = ' '
  labely = ' '
  labelz = ' '
  labelrender = ' '
  labelvecplot = ' '
  xplot = 0.
  yplot = 0.
  zplot = 0.
  dummy = 0.
  hh = 0.
  labeltimeunits = ' '
  dumxsec = .false.
  isetrenderlimits = .false.
  k = nxsec ! matters for lastplot in page_setup for non-coord plots
  if (iReScale) labeltimeunits = unitslabel(0)
  iaxistemp = iaxis
  
  !--set the arrays needed for rendering if they are present
  if (ih.gt.0 .and. ih.le.ndataplots .and. (required(ih) .or. .not.ipartialread)) hh(:) = dat(:,ih)
  if (ipmass.gt.0 .and. ipmass.le.ndataplots) then
     if (required(ipmass)) then
        pmassmin = minval(dat(:,ipmass))
        pmassmax = maxval(dat(:,ipmass))
     else
        pmassmin = 0.
        pmassmax = 0.
     endif
  else
     pmassmin = masstype(1)
     pmassmax = masstype(1)
  endif
  !
  !--set number of particles to use in the interpolation routines
  !  (by default, only the gas particles)
  !
  ntoti = sum(npartoftype)
  ninterp = npartoftype(1)
  if (any(UseTypeInRenderings(2:ntypes).and.iplotpartoftype(2:ntypes)) &
      .or. size(iamtype).gt.1) ninterp = ntoti
  
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
  
  !--set the colour table if it has not been set and particles have been coloured previously
  if (any(icolourme(1:ntoti).gt.16) .and. .not.ihavesetcolours) call colour_set(icolours)
  !
  !--set weight factor for interpolation routines
  !
  call set_interpolation_weights(weight,dat,iamtype,(iplotpartoftype .and. UseTypeInRenderings))
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
     call pgsch(charheight) ! in PGPLOT scaled units

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
     else
        irender = irender_nomulti
        ivectorplot = ivecplot
        icontourplot = icontour_nomulti
        iplotcont = .false. !iplotcont_nomulti
        if (.not.interactivereplot) x_sec = xsec_nomulti
        if (.not.interactivereplot .and. x_sec) zslicepos = xsecpos_nomulti
     endif
     
     !--if the contour plot is the same as the rendered plot,
     !  do not interpolate twice. Instead simply plot the contours
     !  of the rendered quantity when plotting the render plot.
     if (irender.gt.ndim .and. irender.le.numplot) then
        if (icontourplot.eq.irender .and. isameweights .and..not.lim2set(icontourplot)) then
           icontourplot = 0
           iplotcont = .true.
           !print "(a)",' contouring same as rendering'
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

     initdataplots: if (initx .or. inity) then
        if (initx) then
           xplot(1:ntoti) = dat(1:ntoti,iplotx)
           iamvecx = iamvec(iplotx)
        else
           iamvecx = 0
        endif
        if (inity) then
           yplot(1:ntoti) = dat(1:ntoti,iploty)
           iamvecy = iamvec(iploty)
        else
           iamvecy = 0
        endif
        itransx = 0
        itransy = 0
        if (iplotx.gt.0 .and. iplotx.le.numplot) itransx = itrans(iplotx)
        if (iploty.gt.0 .and. iploty.le.numplot) itransy = itrans(iploty)

        zplot = 0.   !--set later if x-sec
        zslicemin = -huge(zslicemax) !-- " " 
        zslicemax = huge(zslicemax)
        if (.not.interactivereplot) then
           if (iplotx.gt.0 .and. iplotx.le.numplot) then
              xmin = lim(iplotx,1)
              xmax = lim(iplotx,2)
           endif
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
        endif
        !
        !--flag for whether or not we have raw particle plot or not
        !  (not allowed to use transformations on coords otherwise)
        !
        rendering = (iplotx.le.ndim .and. iploty.le.ndim .and. &
                     (irenderplot.gt.ndim .or. ivectorplot.gt.0) .and. &
                     (.not.icolour_particles))
        !
        !--change coordinate system if relevant
        !        
        if (icoordsnew.ne.icoords) then
           !--do this if one is a coord but not if rendering
           if (.not.rendering) call changecoords(iplotx,iploty,xplot,yplot,ntoti)
           if (iamvecx.gt.0) call changeveccoords(iplotx,xplot,ntoti)
           if (iamvecy.gt.0) call changeveccoords(iploty,yplot,ntoti)
        endif

        !--apply transformations (log, 1/x etc) if appropriate
        !  also change labels and limits appropriately
        if (.not.(rendering)) then
           if (itransx.ne.0) call applytrans(xplot,xmin,xmax,labelx,itransx,'x',iplotx)
           if (itransy.ne.0) call applytrans(yplot,ymin,ymax,labely,itransy,'y',iploty)
        endif
 
        !--write username, date on plot
        !         if (nacross.le.2.and.ndown.le.2) call pgiden
        !
        !--adjust plot limits if adaptive plot limits set
        !  (find minimum/maximum only on particle types actually plotted)
        !
        if (itrackpart.le.0 .and. .not.irotate) then
           if (initx) call adapt_limits(iplotx,xplot,xmin,xmax,xminadapti,xmaxadapti,'x')
           if (inity) call adapt_limits(iploty,yplot,ymin,ymax,yminadapti,ymaxadapti,'y')
        endif

        !!-reset co-ordinate plot limits if particle tracking           
        if (itrackpart.gt.0 .and. .not.interactivereplot) then
           if (initx) call settrackinglimits(itrackpart,iplotx,xplot,xmin,xmax)
           if (inity) call settrackinglimits(itrackpart,iploty,yplot,ymin,ymax)
        endif

        !--override settings based on positions in sequence
        if (nseq.gt.0) then
           call getsequencepos(iseqpos,iframe,iplotx,iploty,irender, &
                angletempx,angletempy,angletempz,zobservertemp,taupartdepthtemp,&
                zslicepos,xmin,xmax,ymin,ymax,rendermin,rendermax,isetrenderlimits)
        endif
        !--for 3D perspective, do not plot particles behind the observer
        if (ndim.eq.3.and.use3Dperspective) then
           zslicemax = zobservertemp
           if (use3Dopacityrendering) rkappa = rkappafac/taupartdepthtemp
        endif
        
     endif initdataplots

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! plots with co-ordinates as x and y axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if ((iploty.le.ndim).and.(iplotx.le.ndim)) then

        if (.not.interactivereplot .or. irerender) then
           npixx = npix
           if (npixx.le.0) npixx = 500 ! default for other uses of npixx if auto pixels are used
        endif
        
        !!--page setup preliminaries
        if (usesquarexy) then
           just = 1  ! x and y axis have same scale
           ! unless 1D xsec through 2D data or non-cartesian
           if ((irender.gt.ndim .and. ndim.eq.2 .and. x_sec) &
               .or.(icoordsnew.gt.1)) then
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
                  (j.ne.iplotx).and.(j.ne.iploty)) iz = j
           enddo
        endif
        iplotz = iz ! this is used as cross sectioned quantity
        if (iplotz.gt.0 .and. iplotz.le.ndataplots) then
           zplot(1:ntoti) = dat(1:ntoti,iplotz)
           labelz = label(iplotz)
        endif

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
                 allocate(vecplot(2,ninterp),stat=ierr)
                 if (ierr /= 0) then
                    print "(a)",' ERROR allocating memory for vector plot + rotation '
                    stop
                 endif
              endif
              ivectemp = ivectorplot
           endif
           call rotationandperspective(angletempx,angletempy,angletempz,dzscreentemp,zobservertemp, &
                                       xplot,yplot,zplot,ntoti,iplotx,iploty,iplotz,dat,ivectemp,vecplot)
           !--adapt plot limits after rotations have been done
           if (.not.interactivereplot) then
              call adapt_limits(iplotx,xplot,xmin,xmax,xminadapti,xmaxadapti,'x')
              call adapt_limits(iploty,yplot,ymin,ymax,yminadapti,ymaxadapti,'y')
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
              pixwidth = (xmax-xmin)/real(npix)
              npixx = max(int(0.999*(xmax-xmin)/pixwidth) + 1,1)
              npixy = max(int(0.999*(ymax-ymin)/pixwidth) + 1,1)
           else
           !!--automatically reset the pixel number to match the device
              call page_setup(dummy=.true.)
              pixwidth = (xmax-xmin)/real(npixx)
              npixy = max(int(0.999*(ymax-ymin)/pixwidth) + 1,1)
              npixx = max(int(0.999*(xmax-xmin)/pixwidth) + 1,1)
           endif

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
                         pixwidth,inormalise)
                    !--also get contour plot data
                    if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                       if (.not.isameweights) & ! set contouring weights as necessary
                          call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)
 
                       call interpolate2D(xplot(1:ninterp),yplot(1:ninterp), &
                            hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,icontourplot), &
                            icolourme(1:ninterp),ninterp,xmin,ymin,datpixcont,npixx,npixy, &
                            pixwidth,inormalise)
                       gotcontours = .true.

                       if (.not.isameweights) & ! reset weights
                          call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
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
                         inormalise)
                    
                    if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                       !!--allocate memory for 3D contouring array
                       if (allocated(datpixcont3D)) deallocate(datpixcont3D)
                       allocate ( datpixcont3D(npixx,npixy,npixz) )

                       if (.not.isameweights) & ! set contouring weights as necessary
                          call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                       !!--interpolate from particles to 3D grid
                       call interpolate3D(xplot(1:ninterp),yplot(1:ninterp), &
                            zplot(1:ninterp),hh(1:ninterp),weight(1:ninterp), &
                            dat(1:ninterp,icontourplot),icolourme(1:ninterp), &
                            ninterp,xmin,ymin,zmin,datpixcont3D,npixx,npixy,npixz,pixwidth,dz, &
                            inormalise)
                       gotcontours = .true.

                       if (.not.isameweights) & ! reset weights
                          call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
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
              if (iplotpart) then
                 zslicemin = zslicepos-0.5*dz
                 zslicemax = zslicepos+0.5*dz
              endif
           endif

           !------------take projections/cross sections through 3D data-----------------!
           if (irenderplot.gt.ndim .and. ndim.eq.3) then

              !!--allocate memory for 2D rendered array
              if (.not.interactivereplot) then
                 if (allocated(datpix)) then
                    if (npixx.ne.size(datpix(:,1)) .or. npixy.ne.size(datpix(1,:))) then
                       deallocate(datpix)
                       if (allocated(datpixcont)) deallocate(datpixcont)
                       print*,'reallocating...'
                       allocate ( datpix(npixx,npixy) )
                       if (icontourplot.gt.ndim) allocate ( datpixcont(npixx,npixy) )
                    endif
                 else
                    print*,'allocating...'
                    allocate ( datpix(npixx,npixy) )
                    if (icontourplot.gt.ndim) allocate ( datpixcont(npixx,npixy) )
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
                          print*,trim(label(ix(iplotz))),' = ',zslicepos,  &
                               ' : opacity-rendered cross section', xmin,ymin                    
                          if (ipmass.gt.0) then
                             if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),dat(1:ninterp,icontourplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,zslicepos)
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,zslicepos)
                          else
                             if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  masstype(1),1,hh(1:ninterp),dat(1:ninterp,icontourplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,zslicepos)
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               masstype(1),1,hh(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,zslicepos)                          
                          endif
                       elseif (use3Dperspective) then
                          print*,'ERROR: X_SEC WITH 3D PERSPECTIVE NOT IMPLEMENTED'
                          datpix = 0.
                       else
                          !!--do fast cross-section
                          print*,trim(label(ix(iplotz))),' = ',zslicepos,  &
                               ' : fast cross section', xmin,ymin
                          call interpolate3D_fastxsec( &
                               xplot(1:ninterp),yplot(1:ninterp), &
                               zplot(1:ninterp),hh(1:ninterp), &
                               weight(1:ninterp),dat(1:ninterp,irenderplot),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,zslicepos,datpix,npixx,npixy,pixwidth, &
                               inormalise)
                          !!--same but for contour plot
                          if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                             if (.not.isameweights) & ! set contouring weights as necessary
                                call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                             call interpolate3D_fastxsec( &
                                  xplot(1:ninterp),yplot(1:ninterp), &
                                  zplot(1:ninterp),hh(1:ninterp), &
                                  weight(1:ninterp),dat(1:ninterp,icontourplot),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,zslicepos,datpixcont,npixx,npixy,pixwidth, &
                                  inormalise)
                             gotcontours = .true.

                             if (.not.isameweights) & ! reset weights
                                call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
                          endif
                       endif
                    else                 
                       if (use3Dperspective .and. use3Dopacityrendering) then
                          !!--do fast projection with opacity
                          if (ipmass.gt.0) then
                             !--contour plot first
                             if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),dat(1:ninterp,icontourplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,huge(zslicepos))
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               dat(1:ninterp,ipmass),ninterp,hh(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,huge(zslicepos))
                          else
                             !--do contour plot first so brightness corresponds to render plot
                             if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                                if (.not.isameweights) & ! set contouring weights as necessary
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                                call interp3D_proj_opacity( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  masstype(1),1,hh(1:ninterp),dat(1:ninterp,irenderplot), &
                                  dat(1:ninterp,iz),icolourme(1:ninterp), &
                                  ninterp,xmin,ymin,datpixcont,brightness,npixx,npixy,pixwidth,zobservertemp, &
                                  dzscreentemp,rkappa,huge(zslicepos))
                                gotcontours = .true.

                                if (.not.isameweights) & ! reset weights
                                   call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
                             endif
                             call interp3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               masstype(1),1,hh(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz),icolourme(1:ninterp), &
                               ninterp,xmin,ymin,datpix,brightness,npixx,npixy,pixwidth,zobservertemp, &
                               dzscreentemp,rkappa,huge(zslicepos))
                          endif
                       else
                          !!--do fast projection of z integrated data (e.g. column density)
                          call interpolate3D_projection( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                               icolourme(1:ninterp),ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth, &
                               inormalise,zobservertemp,dzscreentemp,ifastrender)
                          !!--same but for contour plot
                          if (icontourplot.gt.ndim .and. icontourplot.le.numplot) then
                             if (.not.isameweights) & ! set contouring weights as necessary
                                call set_interpolation_weights(weight,dat,iamtype,UseTypeInContours)

                             call interpolate3D_projection( &
                                  xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                                  hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,icontourplot), &
                                  icolourme(1:ninterp),ninterp,xmin,ymin,datpixcont,npixx,npixy,pixwidth, &
                                  inormalise,zobservertemp,dzscreentemp,ifastrender)
                             gotcontours = .true.

                             if (.not.isameweights) & ! reset weights
                                call set_interpolation_weights(weight,dat,iamtype,UseTypeInRenderings)
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
           elseif (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) then
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
                       call transform(datpix,itrans(irenderplot),errval=-666.)
                    else
                       call transform(datpix,itrans(irenderplot))
                    endif
                    if (gotcontours) then
                       if (loggedcont) then
                          call transform(datpixcont,itrans(icontourplot),errval=-666.)
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
                    labelrender = integrate_label(labelrender,irender,ix(iz),inormalise)
                    if (gotcontours) labelcont = integrate_label(labelcont,icontourplot,ix(iz),inormalise)
                 endif
                 !!--apply transformations to the label(s) for the rendered and contoured quantit(y,ies)
                 labelrender = transform_label(labelrender,itrans(irenderplot))
                 if (gotcontours) labelcont = transform_label(labelcont,itrans(icontourplot))

                 !!--limits for rendered quantity
                 if (.not.interactivereplot .or. irerender) then
                    !!--find (adaptive) limits of rendered array
                    if (logged) then
!                          rendermin = minval(datpix,mask=datpix.ne.-666.) ! see above
                       renderminadapt = minval(datpix,mask=abs(datpix+666.).gt.tiny(datpix)) ! see above
                    else
                       renderminadapt = minval(datpix)
                    endif
                    rendermaxadapt = maxval(datpix)
                    if (gotcontours) then
                       if (loggedcont) then
                          contminadapt = minval(datpixcont,mask=abs(datpixcont+666.).gt.tiny(datpixcont))
                       else
                          contminadapt = minval(datpixcont)
                       endif
                       contmaxadapt = maxval(datpixcont)
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
                                renderminadapt,rendermaxadapt,trim(labelrender))
              
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

           lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) &
                       .and. nyplot.eq.nyplots .and. k.eq.nxsec)

           !--add to log
           if (x_sec.and.iplotpart.and.iplotz.gt.0) print 35,label(iplotz),zslicemin,label(iplotz),zslicemax
35            format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

           !------------------------------
           ! now actually plot the data
           !------------------------------
           if (irenderplot.gt.ndim) then
              if ((ndim.eq.3).or.(ndim.eq.2.and. .not.x_sec)) then

                 !!--call subroutine to actually render the image
                 call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                   npixx,npixy,xmin,ymin,pixwidth,    &
                   icolours,iplotcont,0,ncontours,.false.,ilabelcont)

                 !!--contour plot of different quantity on top of rendering
                 if (gotcontours) then
                    call render_pix(datpixcont,contmin,contmax,trim(labelcont), &
                         npixx,npixy,xmin,ymin,pixwidth,0,.true.,0,ncontours,.false.,ilabelcont)
                 endif
                              
                 !!--write ppm if interpolate3D_opacity
                 if (use3Dperspective .and. use3Dopacityrendering .and. ndim.eq.3 .and. writeppm) then
                    !!--plot non-gas particle types (e.g. sink particles) on top (and to ppm)
                    call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                      zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                      icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRenderings(:), &
                      (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                      xmin,xmax,ymin,ymax,ifastparticleplot,datpix,npixx,npixy,rendermax,brightness)

                    call write_pixmap_ppm(datpix,npixx,npixy,xmin,ymin,pixwidth,rendermin,rendermax, &
                                          trim(labelrender),((istep-1)*nframesloop + iframe),brightness)                 
                 !!--dump pixmap to file if option set
                 elseif (iwritepixmap) then

                    !!--plot non-gas particle types (e.g. sink particles) on top (and to ppm)
                    call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                      zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                      icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRenderings(:), &
                      (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                      xmin,xmax,ymin,ymax,ifastparticleplot,datpix,npixx,npixy,rendermax)

                    call writepixmap(datpix,npixx,npixy,xmin,ymin,pixwidth,rendermin,rendermax,trim(labelrender),&
                                     ((istep-1)*nframesloop + iframe))
                 !!--no ppm write
                 else
                    !!--plot non-gas particle types (e.g. sink particles) on top
                    call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                      zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                      icolourme(1:ntoti),iamtype,npartoftype(:),PlotOnRenderings(:), &
                      (x_sec.or.use3Dperspective),zslicemin,zslicemax,labelz, &
                      xmin,xmax,ymin,ymax,ifastparticleplot)
                 endif

              elseif (ndim.eq.2 .and. x_sec) then
                 !---------------------------------------------------------------
                 ! plot 1D cross section through 2D data (contents of datpix) 
                 !---------------------------------------------------------------  
                 call pgline(npixx,xgrid,datpix1D)
              endif
           else
              !-----------------------
              ! particle plots
              !-----------------------
              if (iplotpart) then
                 !!--plot all particle types
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),iamtype,npartoftype(:),iplotpartoftype(:), &
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
                      trim(label(ix(iz))(1:index(label(ix(iz)),unitslabel(ix(iz)))-1))//' [code units]'
                 else
                    labelvecplot = '\(2268) '//trim(labelvecplot)//' d'//trim(label(ix(iz)))
                 endif
              endif
              pixwidthvec = (xmax-xmin)/real(npixvec)
              npixyvec = int(0.999*(ymax-ymin)/pixwidthvec) + 1
              pixwidth = (xmax-xmin)/real(npixx) ! used in synchrotron plots

              call vector_plot(ivecx,ivecy,npixvec,npixyvec,pixwidthvec,vecmax,labelvecplot)

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
           call redraw_axes(iaxistemp)
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title
           !
           !--plot exact solution if relevant (before going interactive)
           !
           if (iexact.ne.0) then
              iaxisy = iaxis
              if (tile_plots .and. icolumn.ne.1) iaxisy = -1
              call exact_solution(iexact,iplotx,iploty, &
                   itrans(iplotx),itrans(iploty),icoordsnew, &
                   ndim,ndimV,timei,xmin,xmax,gammai, &
                   xplot(1:npartoftype(1)),yplot(1:npartoftype(1)), &
                   pmassmin,pmassmax,npartoftype(1),imarktype(1), &
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
              if (nacross*ndown.eq.1 .and. nstepsperpage.eq.1) then
                 iadvance = nfreq
                 call interactive_part(ninterp,iplotx,iploty,iplotz,irender,icontourplot,ivecx,ivecy, &
                      xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                      hh(1:ninterp),icolourme(1:ninterp), &
                      xmin,xmax,ymin,ymax,rendermin,rendermax,renderminadapt,rendermaxadapt,contmin,contmax,vecmax, &
                      angletempx,angletempy,angletempz,ndim,xorigin(1:ndim),x_sec,zslicepos,dz, &
                      zobservertemp,dzscreentemp,use3Dopacityrendering,taupartdepthtemp,irerender, &
                      itrackpart,icolours,iColourBarStyle,labelrender,iadvance,ipos,iendatstep,iframe,nframesloop,interactivereplot)
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
                      ivecplottemp(1:nplots),xminmulti(:),xmaxmulti(:),vptxmin(1:nplots),vptxmax(1:nplots), &
                      vptymin(1:nplots),vptymax(1:nplots),barwmulti(1:nplots), &
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

     elseif ((iploty.gt.ndim .or. iplotx.gt.ndim)  &
          .and.(iploty.le.ndataplots .and. iplotx.le.ndataplots)) then
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
                                renderminadapt,rendermaxadapt,trim(labelrender))
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
        call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
             zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
             icolourme(1:ntoti),iamtype,npartoftype(:),iplotpartoftype,.false., &
             zslicemin,zslicemax,' ',xmin,xmax,ymin,ymax,ifastparticleplot)
        !
        !--redraw axes over what has been plotted
        !
        call redraw_axes(iaxis)
        !
        !--annotate with time / marker legend and title
        !
        call legends_and_title
        !
        !--plot exact solution (after redrawn axis for residual plots)
        !
        if (iexact.ne.0) then
           iaxisy = iaxis
           if (tile_plots .and. icolumn.ne.1) iaxisy = -1
           call exact_solution(iexact,iplotx,iploty,itrans(iplotx),itrans(iploty), &
                icoordsnew,ndim,ndimV,timei,xmin,xmax,gammai, &
                xplot(1:npartoftype(1)),yplot(1:npartoftype(1)), &
                pmassmin,pmassmax,npartoftype(1),imarktype(1), &
                units(iplotx),units(iploty),irescale,iaxisy)
        endif
        !
        !--enter interactive mode
        !
        lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) .and. nyplot.eq.nyplots)
        !--the following line sets the number of steps on page to nstepsonpage
        !  in the case where we reach the last timestep before nstepsonpage is reached
        !  (makes interactive replotting behave better)
        if (lastplot) istepsonpage = nstepsperpage

        if (interactive) then
           if (nacross*ndown.eq.1 .and. nstepsperpage.eq.1) then
              iadvance = nfreq
              call interactive_part(ntoti,iplotx,iploty,0,irenderpart,0,0,0, &
                   xplot(1:ntoti),yplot(1:ntoti),zplot(1:ntoti), &
                   hh(1:ntoti),icolourme(1:ntoti), &
                   xmin,xmax,ymin,ymax,rendermin,rendermax,renderminadapt,rendermaxadapt,contmin,contmax,vecmax, &
                   angletempx,angletempy,angletempz,ndim,xorigin(1:ndim), &
                   dumxsec,dummy,dummy,dummy,dummy,.false.,dummy,irerender, &
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
                   ivecplottemp(1:nplots),xminmulti(:),xmaxmulti(:),vptxmin(1:nplots),vptxmax(1:nplots), &
                   vptymin(1:nplots),vptymax(1:nplots),barwmulti(1:nplots), &
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
              label(iploty) = 'PDF '//trim(label(iplotx))
              labely = trim(label(iploty))
           endif

           if (itrans(iploty).ne.0) labely = trim(transform_label(label(iploty),itrans(iploty)))

           if ((.not.interactivereplot) .or. irerender) then
              !--call routines which actually calculate disc properties from the particles
              if (iploty.eq.isurfdens .or. iploty.eq.itoomre) then
                 if (ipmass.gt.0 .and. ipmass.le.ndataplots) then
                    if (iutherm.gt.0 .and. iutherm.le.ndataplots) then
                       call disccalc(itemp,ntoti,xplot(1:ntoti),ntoti,dat(1:ntoti,ipmass), &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty), &
                            gammai,dat(1:ntoti,iutherm))
                    else
                       call disccalc(itemp,ntoti,xplot(1:ntoti),ntoti,dat(1:ntoti,ipmass), &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty),gammai)
                    endif
                 else
                    if (iutherm.gt.0 .and. iutherm.le.ndataplots) then
                       call disccalc(itemp,ntoti,xplot(1:ntoti),1,masstype(1), &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty), &
                            gammai,dat(1:ntoti,iutherm))
                    else
                       call disccalc(itemp,ntoti,xplot(1:ntoti),1,masstype(1), &
                            xmin,xmax,yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty),gammai)
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
                 if (irho.gt.0 .and. irho.le.ndataplots) then
                    call pdfcalc(ntoti,xplot(1:ntoti),xmin,xmax,ngrid,xgrid,datpix1D, &
                         yminadapti,ymaxadapti,itrans(iplotx),itrans(iploty),label(iplotx), &
                         (npdfbins.gt.0),icolourme(1:ntoti),dat(1:ntoti,irho))          
                 else
                    call pdfcalc(ntoti,xplot(1:ntoti),xmin,xmax,ngrid,xgrid,datpix1D,yminadapti,ymaxadapti, &
                         itrans(iplotx),itrans(iploty),label(iplotx),(npdfbins.gt.0),icolourme(1:ntoti))
                 endif
              endif
           endif
           if (iadapt .and. .not.interactivereplot) then
              print "(1x,a)",'adapting '//trim(labely)//' limits'
              ymin = yminadapti
              ymax = ymaxadapti
           endif

           call page_setup

           call pgqci(icolourprev)    ! query line style and colour
           call pgqls(linestyleprev)
           ! set appropriate colour and style if multiple steps per page
           if (nstepsperpage.gt.1) then
              call pgsci(linecolourthisstep)
              call pgsls(linestylethisstep)
           endif
           if (iploty.eq.itoomre .or. iploty.eq.isurfdens) then
              call discplot()
           elseif (iploty.eq.ipdf) then
              call pdfplot(size(xgrid),xgrid,datpix1D)
           endif

           !--restore line size and colour
           call pgsci(icolourprev)
           call pgsls(linestyleprev)

           call redraw_axes(iaxis)
           call legends_and_title

           lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) .and. nyplot.eq.nyplots)
           if (lastplot) istepsonpage = nstepsperpage
           if (interactive .and. ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot)) then
              iadvance = nfreq
              nplots = ipanel
              irerender = .true.
              call interactive_multi(iadvance,ipos,ifirststeponpage,iendatstep,iframe,nframefirstonpage, &
                   nframesloop,ipanel,iplotxtemp(1:nplots),iplotytemp(1:nplots),irendertemp(1:nplots),&
                   ivecplottemp(1:nplots),xminmulti(:),xmaxmulti(:),vptxmin(1:nplots),vptxmax(1:nplots), &
                   vptymin(1:nplots),vptymax(1:nplots),barwmulti(1:nplots), &
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
                 xmin = 1./wavelengthmax  ! freq min
                 xmax = 1./wavelengthmin  ! freq max
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
              labely = transform_label(labely,itrans(iploty))
              if (.not.interactivereplot) then
                 call transform_limits(xmin,xmax,itrans(iploty))
                 call transform_limits(ymin,ymax,itrans(iploty))
              endif
           endif
           
           just = 0
           call page_setup

           call pgqci(icolourprev)    ! query line style and colour
           call pgqls(linestyleprev)
           if (nstepsperpage.gt.1) then
              call pgsci(linecolourthisstep) ! set appropriate colour and style if multiple steps per page
              call pgsls(linestylethisstep)
           endif
           
           call pgline(nfreqpts,xplot(1:nfreqpts),yplot(1:nfreqpts))
           print*,' maximum power at '//trim(labelx)//' = ',xplot(maxloc(yplot(1:nfreqpts)))

           call pgsci(icolourprev)
           call pgsls(linestyleprev)
           
           !
           !--redraw axes over what has been plotted
           !
           call redraw_axes(iaxis)           
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title

        else
        !--------------------------------------------------------------
        !  plot the contents of an extra two-column ascii file
        !--------------------------------------------------------------           

           call exact_fromfile('gwaves1.dat',xplot,yplot,nfreqpts,ierr)
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
              call pgpt1(xplot(ipt),yplot(ipt),4)
              call pgline(ipt,xplot(1:ipt),yplot(1:ipt))
           endif
           
           call redraw_axes(iaxis)
           call legends_and_title
        
        endif

        lastplot = (ipos.eq.iendatstep .or. istep.eq.nsteps)
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
        call pgpage! just skip to next plot

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
  
  return
    
contains

!----------------------------------------------
! interfaces to the page setup routines
! this is called just before a plot is
! actually plotted
!----------------------------------------------
  subroutine page_setup(dummy)
    use colourbar, only:get_colourbarmargins
    use pagesetup, only:setpage2
    use settings_page, only:nstepsperpage,iUseBackgroundColourForAxes,vposlegend,iPlotLegend
    implicit none
    integer :: iplotsave,ipanelsave
    real :: barwidth, TitleOffset,xminmargin,xmaxmargin,yminmargin,ymaxmargin
    real :: xminpix,xmaxpix,yminpix,ymaxpix,dxpix
    logical :: ipanelchange,dum
    logical, intent(in), optional :: dummy

    !--------------------------------------------
    ! whether or not this is a dummy call or not
    !--------------------------------------------
    if (present(dummy)) then
       dum = dummy
    else
       dum = .false.
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
    icolumn = ipanel - ((ipanel-1)/nacross)*nacross
    irow = (ipanel-1)/nacross + 1

    !--if we are in interactive mode, use the currently buffered plot limits
    if (interactivereplot .and. (nacross*ndown.gt.1 .or. nstepsperpage.gt.1)) then
       xmin = xminmulti(iplotx)
       xmax = xmaxmulti(iplotx)
       ymin = xminmulti(iploty)
       ymax = xmaxmulti(iploty)
       if (ivectorplot.gt.0 .and. ivectorplot.le.numplot) then
          vecmax = xmaxmulti(ivectorplot)
       endif
    endif

    !--------------------------------------------------------------
    ! output some muff to the screen
    !--------------------------------------------------------------

    if (interactive .and. .not.dum) then
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
    if (.not.dum) call pgsci(1)

    !--page margins: zero if no box is drawn
    if (iaxistemp.eq.-2) then
       xminmargin = 0.0
       xmaxmargin = 0.0
       yminmargin = 0.0
       ymaxmargin = 0.0
    else ! leave a small buffer to allow for line width changes
       xminmargin = 0.001
       xmaxmargin = 0.001
       yminmargin = 0.001
       ymaxmargin = 0.001
    endif
    
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
       call pgpage
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
    if (nstepsperpage.ne.0 .or. inewpage) then
       if (dum) then !--fake the page setup, then return
          if (.not.(interactivereplot .and. .not.irerender)) then
             call setpage2(ipanel,nacross,ndown,xmin,xmax,ymin,ymax, &
                     trim(labelx),trim(labely),'NOPGBOX',just,iaxistemp, &
                     xminmargin,xmaxmargin,yminmargin,ymaxmargin, &
                     0.0,TitleOffset,isamexaxis,tile_plots)
             call pgqvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
             npixx = nint(xmaxpix-xminpix)
             npixy = nint(ymaxpix-yminpix)
             if (vectordevice .and. npixx.gt.1024) then
                npixx = 1024/nacross
                dxpix = (xmax-xmin)/npixx
                npixy = int(0.999*(ymax-ymin)/real(dxpix)) + 1
                print "(a,i4,a,i4,a)",' auto-selecting resolution of ',npixx,' x ',npixy,' for vector device'
                print "(a)",' => set the number of pixels manually if you want more (or less) than this.'
             else
                print "(a,i4,a,i4)",' auto-selecting device resolution = ',npixx,' x ',npixy
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
          !--restore saved attributes
          iplots = iplotsave
          ipanel = ipanelsave
          return
       else
          call setpage2(ipanel,nacross,ndown,xmin,xmax,ymin,ymax, &
                  trim(labelx),trim(labely),' ',just,iaxistemp, &
                  xminmargin,xmaxmargin,yminmargin,ymaxmargin, &
                  0.0,TitleOffset,isamexaxis,tile_plots)
       endif
    endif
    
    !--query and save viewport co-ordinates set up for this panel
    call pgqvp(0,vptxmin(ipanel),vptxmax(ipanel),vptymin(ipanel),vptymax(ipanel))

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
    ivecplottemp(ipanel) = ivectorplot
    xminmulti(iplotx) = xmin
    xmaxmulti(iplotx) = xmax
    xminmulti(iploty) = ymin
    xmaxmulti(iploty) = ymax
    if (irender.gt.0 .and. irender.le.numplot) then
       xminmulti(irender) = rendermin
       xmaxmulti(irender) = rendermax
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
       if (irender.gt.0 .and. irender.le.numplot) then
          xminadapt(irender) = min(xminadapt(irender),renderminadapt)
          xmaxadapt(irender) = max(xmaxadapt(irender),rendermaxadapt)
       endif
    !endif
    
    !--change to background colour index for overlaid text and axes
    if (iUseBackGroundColourForAxes) call pgsci(0)
    
    return
  end subroutine page_setup
  
!------------------------------------------------------
! draws legend(s), titles etc
! (must be called after rendering otherwise rendering
!  will overwrite plot area)
!------------------------------------------------------
  subroutine legends_and_title
    use colourbar, only:plotcolourbar
    use legends, only:legend,legend_markers,legend_scale
    use titles, only:pagetitles,steplegend,lensteplegend
    use filenames, only:nstepsinfile,nfiles,rootname
    use settings_page, only:iPlotLegend,iPlotStepLegend, &
        hposlegend,vposlegend,fjustlegend,legendtext,iPlotLegendOnlyOnPanel, &
        iPlotScale,iscalepanel,dxscale,hposscale,vposscale,scaletext,iUseBackGroundColourForAxes
    use shapes, only:nshapes,plot_shapes
    use pagesetup, only:xlabeloffset
    implicit none
    integer :: icoloursave
    character(len=lensteplegend) :: steplegendtext
    real :: xlabeloffsettemp
    
    !--save colour index
    call pgqci(icoloursave)

    !--use foreground colour by default for legends
    call pgsci(1)

    !--------------------------------------------------------------
    ! plot colour bar for rendered plots (use currently set colour)
    ! do this here so it always appears OVERLAID on the renderings
    !--------------------------------------------------------------
    if (irender.gt.0 .and. irender.le.numplot) then
       lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) &
                          .and. nyplot.eq.nyplots .and. k.eq.nxsec)
       !--only plot colour bar at the end of first row on tiled plots
       if (tile_plots .and..not.(ipanel.eq.nacross*ndown .or. lastplot)) iPlotColourBar = .false.

       if (iPlotColourBar) then
          xlabeloffsettemp = xlabeloffset + 1.0
          if (iaxistemp.lt.0) xlabeloffsettemp = 0.

          !--for tiled plots only on last plot in first row,
          !  and use full viewport size in the y direction
          if (tile_plots) then
             call plotcolourbar(iColourBarStyle,icolours,rendermin,rendermax, &
                  trim(labelrender),.false.,xlabeloffsettemp, &
                  minval(vptxmin(1:ipanel)),maxval(vptxmax(1:ipanel)), &
                  minval(vptymin(1:ipanel)),maxval(vptymax(1:ipanel)))
          else
             !!--plot colour bar, but only if last in row
             call plotcolourbar(iColourBarStyle,icolours,rendermin,rendermax, &
                            trim(labelrender),.false.,xlabeloffsettemp)
          endif
       endif
    endif
    
    !--plot time on plot
    if (iPlotLegend .and. nyplot.eq.1 &
        .and..not.(iPlotLegendOnlyOnPanel.gt.0 .and. ipanel.ne.iPlotLegendOnlyOnPanel) &
        .and..not.(iPlotLegendOnlyOnPanel.eq.-1 .and. irow.gt.1) &
        .and..not.(iPlotLegendOnlyOnPanel.eq.-2 .and. icolumn.gt.1) &
        .and. timei.gt.-0.5*huge(timei)) then  ! but not if time has not been read from dump

       !--change to background colour index for legend text if overlaid
       if (iUseBackGroundColourForAxes .and. vposlegend.gt.0.) call pgsci(0)

       call legend(legendtext,timei,labeltimeunits,hposlegend,vposlegend,fjustlegend)
    endif

    !--line/marker style/colour legend for multiple timesteps on same page
    if (iPlotStepLegend .and. istepsonpage.gt.0 &
        .and.((nyplot.eq.1 .and. iPlotLegendOnlyOnPanel.eq.0) &
        .or.(iPlotLegendOnlyOnPanel.gt.0 .and. ipanel.eq.iPlotLegendOnlyOnPanel) &
        .or.(iPlotLegendOnlyOnPanel.eq.-1 .and. irow.eq.1) &
        .or.(iPlotLegendOnlyOnPanel.eq.-2 .and. icolumn.eq.1))) then

       !--change to background colour index for overlaid text and axes
       if (iUseBackGroundColourForAxes .and. vposlegend.gt.0.) call pgsci(0)
       !
       !--use filenames in legend if none set
       !
       if (istepsonpage.le.nsteplegendlines) then
          steplegendtext = steplegend(istepsonpage)
       elseif (all(nstepsinfile(1:nfiles).le.1)) then
          steplegendtext = trim(rootname(istep))
       else
          write(steplegendtext,"(a,i4)") 'step ',istep
       endif
       if (iploty.gt.ndataplots) then
          call legend_markers(istepsonpage,linecolourthisstep,imarktype(1),linestylethisstep, &
            .false.,.true.,trim(steplegendtext),hposlegend,vposlegend)           
       else
          call legend_markers(istepsonpage,linecolourthisstep,imarktype(1),linestylethisstep, &
            iplotpartoftype(1),iplotline,trim(steplegendtext),hposlegend,vposlegend)
       endif
    endif
    
    !--use foreground colour by default for title
    call pgsci(1)

    !--print title if appropriate
    if (iPlotTitles .and. istepsonpage.eq.1 .and. ipanel.le.ntitles) then
       if (len_trim(pagetitles(ipanel)).gt.0) then
          
          !--change to background colour index if title is overlaid
          if (iUseBackGroundColourForAxes .and. vpostitle.lt.0.) call pgsci(0)

          call pgmtxt('T',vpostitle,hpostitle,fjusttitle,trim(pagetitles(ipanel)))
       endif
    endif

    !--use foreground colour by default for scale
    call pgsci(1)
    
    !--scale on co-ordinate plots
    if (iPlotScale .and. (iscalepanel.eq.0 .or. ipanel.eq.iscalepanel) &
                   .and. any(ix(1:ndim).eq.iplotx) .and. any(ix(1:ndim).eq.iploty)) then

       !--change to background colour index if title is overlaid
       if (iUseBackGroundColourForAxes .and. vposscale.gt.0.) call pgsci(0)

       call legend_scale(dxscale,hposscale,vposscale,scaletext)
    endif
    
    !--plot shapes
    if (nshapes.gt.0) call plot_shapes(ipanel,irow,icolumn,itrans(iplotx),itrans(iploty))
    
    !--restore colour index
    call pgsci(icoloursave)
    
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

!----------------------------------------------------
! adapt the (particle plot) limits to include all
! particles which are to be plotted on the page
!---------------------------------------------------
  
  subroutine adapt_limits(iplot,xploti,xmini,xmaxi,xminadaptive,xmaxadaptive,labeli)
    implicit none
    integer, intent(in) :: iplot
    real, dimension(:), intent(in) :: xploti
    real, intent(inout) :: xmini,xmaxi,xminadaptive,xmaxadaptive
    character(len=*), intent(in) :: labeli
    integer :: index1,index2,itype    
    !--calculate adaptive limits for this quantity
    xminadaptive = huge(xminadaptive)
    xmaxadaptive = -huge(xmaxadaptive)
    index1 = 1
    do itype=1,maxparttypes
       index2 = index1 + npartoftype(itype) - 1
       if (iplotpartoftype(itype).and.npartoftype(itype).gt.0 &
          .or. (iplotline.and.itype.eq.1)) then
          xminadaptive = min(xminadaptive,minval(xploti(index1:index2)))
          xmaxadaptive = max(xmaxadaptive,maxval(xploti(index1:index2))*scalemax)
       endif
       index1 = index2 + 1
    enddo
    
    !--set these as limits if adaptive limits are on
    if (.not.interactivereplot) then
       if ((iplot.le.ndim .and. iadaptcoords) &
       .or.(iplot.gt.ndim .and. iadapt) .and. ipagechange) then
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
!-------------------------------------------------------------------
  subroutine applytrans(xploti,xmini,xmaxi,labelxi,itransxi,chaxis,iplotxi)
    implicit none
    integer, intent(in) :: itransxi,iplotxi
    real, dimension(:), intent(inout) :: xploti
    real, intent(inout) :: xmini,xmaxi
    character(len=*), intent(inout) :: labelxi
    character(len=1), intent(in) :: chaxis
    integer :: itranstemp
    character(len=20) :: string
       
    if (itransxi.ne.0) then
        if (iplotxi.gt.0 .and. iplotxi.le.numplot) call transform(xploti(:),itransxi)
        if ((chaxis.eq.'x' .and. iaxis.eq.10 .or. iaxis.eq.30).or. &
            (chaxis.eq.'y' .and. iaxis.eq.20 .or. iaxis.eq.30)) then ! logarithmic axes
           write(string,*) itransx
           string = adjustl(string)
           itranstemp = 0
           if (string(len_trim(string):len_trim(string)).eq.'1') then  
              if (len_trim(string).gt.1) read(string(1:len_trim(string)-1),*) itranstemp
              labelxi = transform_label(labelxi,itranstemp)
           else
              labelxi = transform_label(labelxi,itransxi)
           endif
        else
           labelxi = transform_label(labelxi,itransxi)
        endif
        if (.not.interactivereplot) call transform_limits(xmini,xmaxi,itransxi)
     endif
  end subroutine applytrans

!-------------------------------------------------------------------
! interface to coordinate-system transformations
!-------------------------------------------------------------------
  subroutine changecoords(iplotx,iploty,xplot,yplot,ntot)
   use geometry, only:coord_transform,labelcoordsys
   use settings_data, only:xorigin
   implicit none
   integer, intent(in) :: iplotx,iploty,ntot
   real, dimension(:), intent(inout) :: xplot,yplot
   real, dimension(ndim) :: xcoords,xcoordsnew
   integer :: j

   if (iplotx.le.ndim .or. iploty.le.ndim) then
      print*,'changing coords from ',trim(labelcoordsys(icoords)), &
             ' to ',trim(labelcoordsys(icoordsnew))
      if (itrackpart.gt.0) print*,' (relative to particle ',itrackpart,')'

      do j=1,ntot
         if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
            xcoords(1:ndim) = dat(j,ix(1:ndim)) - dat(itrackpart,ix(1:ndim))
         else
            xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
         endif

         call coord_transform(xcoords(1:ndim),ndim,icoords, &
                              xcoordsnew(1:ndim),ndim,icoordsnew)
         if (iplotx.le.ndim) xplot(j) = xcoordsnew(iplotx)
         if (iploty.le.ndim) yplot(j) = xcoordsnew(iploty)
      enddo
   endif
   
  end subroutine changecoords
  
!-------------------------------------------------------------------
! interface to coordinate-system transformations for vectors
!-------------------------------------------------------------------
  subroutine changeveccoords(iplot,xploti,ntot)
   use geometry, only:vector_transform,labelcoordsys
   use settings_data, only:xorigin
   use labels, only:ivx
   implicit none
   integer, intent(in) :: iplot,ntot
   real, dimension(:), intent(inout) :: xploti
   integer :: j
   real, dimension(ndim) :: xcoords,vecnew,vecin
   
   if (iamvec(iplot).gt.0) then
      if (iplot-iamvec(iplot)+1 .le. ndim) then
         print*,'changing vector component from ', &
          trim(labelcoordsys(icoords)),' to ',trim(labelcoordsys(icoordsnew))
         if (itrackpart.gt.0 .and. iamvec(iplot).eq.ivx) print*,' (velocities relative to particle ',itrackpart,')'
         do j=1,ntot
            if (itrackpart.gt.0 .and. itrackpart.le.ntot) then
               xcoords(1:ndim) = dat(j,ix(1:ndim)) - dat(itrackpart,ix(1:ndim))
               if (iamvec(iplot).eq.ivx) then
                  vecin(1:ndim) = dat(j,iamvec(iplot):iamvec(iplot)+ndim-1) - dat(itrackpart,iamvec(iplot):iamvec(iplot)+ndim-1)
               else
                  vecin(1:ndim) = dat(j,iamvec(iplot):iamvec(iplot)+ndim-1)
               endif
            else
               xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
               vecin(1:ndim) = dat(j,iamvec(iplot):iamvec(iplot)+ndim-1)
            endif
            call vector_transform(xcoords(1:ndim),vecin(1:ndim), &
                 ndim,icoords,vecnew(1:ndim),ndim,icoordsnew)
            xploti(j) = vecnew(iplot-iamvec(iplot)+1)
         enddo
      else
         print*,'error: can''t convert vector components with ndimV > ndim'
      endif
   endif
   
   return
  end subroutine changeveccoords
  
!-------------------------------------------------------------------
! interface for setting limits when using particle tracking limits
!-------------------------------------------------------------------
  subroutine settrackinglimits(itrackpart,iplot,xploti,xmini,xmaxi)
    use settings_limits, only:xminoffset_track,xmaxoffset_track
    implicit none
    integer, intent(in) :: itrackpart,iplot
    real, dimension(:), intent(in) :: xploti
    real, intent(inout) :: xmini,xmaxi
    
    !--particle tracking limits only apply to co-ordinate axes
    if (iplot.le.ndim) then
       xmini = xploti(itrackpart) - xminoffset_track(iplot)
       xmaxi = xploti(itrackpart) + xmaxoffset_track(iplot)
    endif
    
    return
   end subroutine settrackinglimits

!-------------------------------------------------------------------
! interface for setting interpolation weights
!-------------------------------------------------------------------
  subroutine set_interpolation_weights(weighti,dati,iamtypei,usetype)
    use settings_render, only:idensityweightedinterpolation
    implicit none
    real, dimension(:), intent(out) :: weighti
    real, dimension(:,:), intent(in) :: dati
    integer(kind=int1), dimension(:), intent(in) :: iamtypei
    logical, dimension(:), intent(in) :: usetype
    integer :: i2,i1,itype,ipart
    real(doub_prec) :: dunitspmass,dunitsrho,dunitsh

    dunitspmass = 1.d0
    dunitsrho = 1.d0
    dunitsh = 1.d0
    if (iRescale) then
       if (ipmass.gt.0) dunitspmass = 1.d0/units(ipmass)
       if (ih.gt.0) dunitsh = 1.d0/units(ih)
       if (irho.gt.0) dunitsrho = 1.d0/units(irho)
    endif

    if (ipmass.gt.0 .and. ipmass.le.ndataplots .and. &
        irho.gt.0 .and. irho.le.ndataplots .and. &
        ih .gt. 0 .and. ih.le.ndataplots .and. &
        required(ipmass) .and. required(irho) .and. required(ih)) then
       
       if (size(iamtype).gt.1) then
          !
          !--particles with mixed types
          !
          do ipart=1,ninterp
             itype = iamtypei(ipart)
             if (.not.usetype(itype)) then
                weighti(ipart) = 0.
             elseif (idensityweightedinterpolation) then
                if (dati(ipart,ih) > tiny(dati)) then
                   weighti(ipart) = (dati(ipart,ipmass)*dunitspmass)/ &
                                   ((dati(ipart,ih)*dunitsh)**ndim)
                else
                   weighti(ipart) = 0.
                endif
             else
                if (dati(ipart,irho) > tiny(dati) .and. dati(ipart,ih) > tiny(dati)) then
                   weighti(ipart) = (dati(ipart,ipmass)*dunitspmass)/ &
                                   ((dati(ipart,irho)*dunitsrho)*(dati(ipart,ih)*dunitsh)**ndim)
                else
                   weighti(ipart) = 0.
                endif
             endif
          enddo
       else
          !
          !--particles ordered by type
          !
          i2 = 0
          do itype=1,ntypes
             i1 = i2 + 1
             i2 = i2 + npartoftype(itype)
             !--set weights to zero for particle types not used in the rendering
             if (.not.usetype(itype)) then
                weighti(i1:i2) = 0.
             elseif (idensityweightedinterpolation) then
             !--for density weighted interpolation use m/h**ndim
                where(dati(i1:i2,ih) > tiny(dati))
                   weighti(i1:i2) = (dati(i1:i2,ipmass)*dunitspmass)/ &
                                    ((dati(i1:i2,ih)*dunitsh)**ndim)
                elsewhere
                   weighti(i1:i2) = 0.
                endwhere
             else
             !--usual interpolation use m/(rho h**ndim)
                where(dati(i1:i2,irho) > tiny(dati) .and. dati(i1:i2,ih) > tiny(dati))
                   weighti(i1:i2) = (dati(i1:i2,ipmass)*dunitspmass)/ &
                                    ((dati(i1:i2,irho)*dunitsrho)*(dati(i1:i2,ih)*dunitsh)**ndim)
                elsewhere
                   weighti(i1:i2) = 0.
                endwhere
             endif
          enddo
       endif
       
       if (idensityweightedinterpolation) then
          print "(a)",' USING DENSITY WEIGHTED INTERPOLATION '
          inormalise = .true.
       else
          inormalise = inormalise_interpolations
       endif

    elseif (masstype(1).gt.0. .and. &
            irho.gt.0 .and. irho.le.ndataplots .and. &
            ih .gt. 0 .and. ih.le.ndataplots .and. &
            required(irho) .and. required(ih)) then
       
       if (size(iamtype).gt.1) then
          !
          !--particles with mixed types
          !
          do ipart=1,ninterp
             itype = iamtypei(ipart)
             if (.not.usetype(itype)) then
                weighti(ipart) = 0.
             elseif (idensityweightedinterpolation) then
                if (dati(ipart,ih) > tiny(dati)) then
                   weighti(ipart) = (masstype(itype)*dunitspmass)/ &
                                   ((dati(ipart,ih)*dunitsh)**ndim)
                else
                   weighti(ipart) = 0.
                endif
             else
                if (dati(ipart,irho) > tiny(dati) .and. dati(ipart,ih) > tiny(dati)) then
                   weighti(ipart) = (masstype(itype)*dunitspmass)/ &
                                   ((dati(ipart,irho)*dunitsrho)*(dati(ipart,ih)*dunitsh)**ndim)
                else
                   weighti(ipart) = 0.
                endif
             endif
          enddo
       else
          !
          !--particles ordered by type
          !
          i2 = 0
          do itype=1,ntypes
             i1 = i2 + 1
             i2 = i2 + npartoftype(itype)
             !--set weights to zero for particle types not used in the rendering
             if (.not.usetype(itype)) then
                weighti(i1:i2) = 0.
             else
                where(dati(i1:i2,irho) > tiny(dati) .and. dati(i1:i2,ih) > tiny(dati))
                   weighti(i1:i2) = masstype(itype)/ &
                                  ((dati(i1:i2,irho)*dunitsrho)*(dati(i1:i2,ih)*dunitsh)**ndim)
                elsewhere
                   weighti(i1:i2) = 0.
                endwhere
             endif
          enddo
       endif
       
       if (idensityweightedinterpolation) then
          print "(a)",' USING DENSITY WEIGHTED INTERPOLATION '
          inormalise = .true.
       else
          inormalise = inormalise_interpolations
       endif
    else
       if (required(ih) .and. required(irho)) then 
          print "(a)",' WARNING: particle mass not set: using normalised interpolations'
       endif
    !--if particle mass has not been set, then must use normalised interpolations
       weight(1:ninterp) = 1.0
       inormalise = .true.
    endif

  end subroutine set_interpolation_weights

!-------------------------------------------------------------------
! interface to vector plotting routines
! so that pixel arrays are allocated appropriately
!-------------------------------------------------------------------
  subroutine vector_plot(ivecx,ivecy,numpixx,numpixy,pixwidthvec,vmax,label)
   use settings_vecplot, only:UseBackgndColorVecplot,iplotstreamlines,iplotarrowheads, &
       iplotsynchrotron,rcrit,zcrit,synchrotronspecindex,uthermcutoff, &
       ihidearrowswherenoparts,minpartforarrow
   use interpolations2D, only:interpolate2D_vec
   use projections3D, only:interpolate3D_proj_vec,interp3D_proj_vec_synctron
   use interpolate_vec, only:mask_vectors
   use render, only:render_vec
   use fieldlines, only:streamlines
   use labels, only:iutherm
   implicit none
   integer, intent(in) :: ivecx,ivecy,numpixx,numpixy
   real, intent(in) :: pixwidthvec
   real, intent(inout) :: vmax
   character(len=*), intent(in) :: label
   real, dimension(numpixx,numpixy) :: vecpixx, vecpixy
   real, dimension(max(npixx,numpixx),max(npixy,numpixy)) :: datpixvec
   integer :: i,j,icoloursav,linewidthprev
   real :: vmag
   real :: blankval,datmax
   logical :: usevecplot
   
   !--query colour index and line width
   call pgqci(icoloursav)
   call pgqlw(linewidthprev)

   !print*,'plotting vector field ',trim(label)
   if ((ivecx.le.ndim).or.(ivecx.gt.ndataplots) &
        .or.(ivecy.le.ndim).or.(ivecy.gt.ndataplots)) then
      print*,'error finding location of vector plot in array'
   else
      !--plot arrows in either background or foreground colour
      if (UseBackgndColorVecplot) then
         call pgsci(0)
      else
         call pgsci(1)
      endif
      usevecplot = .false.
      if (irotate) then
         if (allocated(vecplot)) usevecplot = .true.
         print*,'using vecplot' ! this is to indicate (to me) that extra memory is in use
      endif
      !
      !--interpolate using appropriate routine for number of dimensions
      !
      select case(ndim)
      case(3)
         if (x_sec) then ! take vector plot in cross section
            if (usevecplot) then ! using rotation
               call interpolate3D_xsec_vec(xplot(1:ninterp), &
                 yplot(1:ninterp),zplot(1:ninterp), &
                 hh(1:ninterp),weight(1:ninterp), &
                 vecplot(1,1:ninterp),vecplot(2,1:ninterp), &
                 icolourme(1:ninterp),ninterp,xmin,ymin,zslicepos, &
                 vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,inormalise)            
            else
               call interpolate3D_xsec_vec(xplot(1:ninterp), &
                 yplot(1:ninterp),zplot(1:ninterp), &
                 hh(1:ninterp),weight(1:ninterp), &
                 dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                 icolourme(1:ninterp),ninterp,xmin,ymin,zslicepos, &
                 vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,inormalise)
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
                  else
                     call interp3D_proj_vec_synctron(xplot(1:ninterp), &
                       yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                       weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                       icolourme(1:ninterp),ninterp,xmin,ymin, &
                       vecpixx,vecpixy,datpixvec(1:numpixx,1:numpixy),numpixx,numpixy,pixwidthvec, &
                       rcrit,zcrit,synchrotronspecindex,pixwidthvec,.false.)
                  endif
               endif
            else
            !   call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
            !     dat(1:ninterp,ivecx),dat(1:ninterp,ivecy),icolourme(1:ninterp), &
            !     xmin,ymin,pixwidth,vecpixx,vecpixy, &
            !     ninterp,numpixx,numpixy)

               if (usevecplot) then
                  call interpolate3D_proj_vec(xplot(1:ninterp), &
                    yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                    weight(1:ninterp),vecplot(1,1:ninterp),vecplot(2,1:ninterp), &
                    icolourme(1:ninterp),ninterp,xmin,ymin, &
                    vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,.false.,zobservertemp,dzscreentemp)
               else
                  call interpolate3D_proj_vec(xplot(1:ninterp), &
                    yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
                    weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
                    icolourme(1:ninterp),ninterp,xmin,ymin, &
                    vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,.false.,zobservertemp,dzscreentemp)               
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
         !--or interpolate (via averaging) to coarser grid
         !
         !call fieldlines2D(ninterp,xplot(1:ninterp),yplot(1:ninterp), &
         !     dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
         !     hh(1:ninterp),pmass(1:ninterp), &
         !     rho(1:ninterp),xmin,xmax,ymin,ymax)
         !call interpolate_vec_average(xplot(1:ninterp),yplot(1:ninterp), &
         !  dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
         !  xmin,ymin,pixwidthvec,vecpixx,vecpixy, &
         !  ninterp,numpixx,numpixy)
         
         if (usevecplot) then
            call interpolate2D_vec(xplot(1:ninterp),yplot(1:ninterp), &
              hh(1:ninterp),weight(1:ninterp),vecplot(2,1:ninterp), &
              vecplot(2,1:ninterp),icolourme(1:ninterp),ninterp,xmin,ymin, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,inormalise)         
         else
            call interpolate2D_vec(xplot(1:ninterp),yplot(1:ninterp), &
              hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,ivecx), &
              dat(1:ninterp,ivecy),icolourme(1:ninterp),ninterp,xmin,ymin, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidthvec,inormalise)
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
                         'crap',numpixx,numpixy,xmin,ymin,pixwidthvec,    &
                         0,.true.,0,ncontours,.false.,ilabelcont,blank=blankval)
         else
            call render_pix(datpixvec(1:numpixx,1:numpixy), &
                         minval(datpixvec(1:numpixx,1:numpixy)), &
                         maxval(datpixvec(1:numpixx,1:numpixy)), &
                         'crap',numpixx,numpixy,xmin,ymin,pixwidthvec,    &
                         0,.true.,0,ncontours,.false.,ilabelcont)
         endif

      else
         if (ihidearrowswherenoparts) then
            call mask_vectors(xplot(1:ninterp),yplot(1:ninterp),icolourme(1:ninterp),ninterp, &
                              xmin,xmax,ymin,ymax,vecpixx,vecpixy,numpixx,numpixy,minpartforarrow,0.)
         endif
      
         call render_vec(vecpixx,vecpixy,vmax, &
              numpixx,numpixy,xmin,ymin,pixwidthvec,trim(label),' ')
         
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
              npixx,npixy,xmin,ymin,pixwidth,0,.true.,0,ncontours,.false.,ilabelcont)
         endif

      endif

   endif
   
   !--restore colour index and line width
   call pgsci(icoloursav)
   call pgslw(linewidthprev)
  
  end subroutine vector_plot

end subroutine plotstep

!-------------------------------------------------------------------
! interface for adding rotation and perspective
! (completely independent)
!-------------------------------------------------------------------
subroutine rotationandperspective(anglexi,angleyi,anglezi,dzscreen,zobs,xploti,yploti,zploti, &
                                  ntot,iplotx,iploty,iplotz,dat,ivecstart,vecploti)
  use labels, only:ix
  use settings_data, only:ndim,xorigin,itrackpart
  use settings_xsecrot, only:use3Dperspective
  use rotation, only:rotate2D,rotate3D
  implicit none
  real, intent(in) :: anglexi,angleyi,anglezi,dzscreen,zobs
  real, dimension(:), intent(inout) :: xploti,yploti,zploti
  real, dimension(:,:), intent(in) :: dat
  real, dimension(:,:), intent(out) :: vecploti
  integer, intent(in) :: ntot,iplotx,iploty,iplotz,ivecstart
  integer :: j
  real :: angleradx,anglerady,angleradz
  real, dimension(ndim) :: xcoords,veci
  !
  !--convert angles to radians
  !
  angleradz = anglezi*pi/180.
  anglerady = angleyi*pi/180.
  angleradx = anglexi*pi/180.
  print "(1x,a,f6.2)",'rotating particles about z by ',anglezi
  if (ndim.eq.3) then
     print "(1x,a,f6.2)",'rotating particles about y by ',angleyi
     print "(1x,a,f6.2)",'rotating particles about x by ',anglexi
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

  if (ivecstart.gt.0) print "(1x,a)",'(also rotating vector components)'

!$omp parallel default(none) &
!$omp shared(dat,xorigin,ndim,angleradx,anglerady,angleradz,zobs,dzscreen) &
!$omp shared(xploti,yploti,zploti,iplotx,iploty,iplotz,ntot,ix,itrackpart) &
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
        xploti(j) = xcoords(iplotx) + dat(itrackpart,ix(iplotx))
        yploti(j) = xcoords(iploty) + dat(itrackpart,ix(iploty))
        if (iplotz.gt.0) then
           zploti(j) = xcoords(iplotz) + dat(itrackpart,ix(iplotz))
        endif     
     else
        xploti(j) = xcoords(iplotx) + xorigin(iplotx)
        yploti(j) = xcoords(iploty) + xorigin(iploty)
        if (iplotz.gt.0) then
           zploti(j) = xcoords(iplotz) + xorigin(iplotz)
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
        vecploti(1,j) = veci(iplotx)
        vecploti(2,j) = veci(iploty)
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
     call rotate_axes3D(irotateaxes,iplotx,iploty, &
          xminrotaxes(1:ndim),xmaxrotaxes(1:ndim),xorigin(1:ndim), &
          angleradx,anglerady,angleradz,zobs,dzscreen)
  elseif (ndim.eq.2) then
     call rotate_axes2D(irotateaxes,xminrotaxes(1:ndim), &
                       xmaxrotaxes(1:ndim),xorigin(1:ndim),angleradz)
  endif

  return
end subroutine rotatedaxes

!---------------------------------------------------------------
! interface for adjusting the label for column-integrated plots
!---------------------------------------------------------------
function integrate_label(labelin,iplot,izcol,normalise)
  use settings_data, only:iRescale
  use settings_units, only:labelzintegration,unitslabel
  use settings_render, only:projlabelformat,iapplyprojformat
  use labels, only:irho,label
  use asciiutils, only:string_replace
  implicit none
  character(len=*), intent(in) :: labelin
  integer, intent(in) :: iplot,izcol
  logical, intent(in) :: normalise
  character(len=len(label)+20) :: integrate_label

  if (len_trim(projlabelformat).ne.0 .and. (iapplyprojformat.eq.0 .or. iapplyprojformat.eq.iplot)) then
     integrate_label = projlabelformat
     call string_replace(integrate_label,'%l',trim(labelin))
     if (iRescale) then
        call string_replace(integrate_label,'%z',trim(label(izcol)(1:index(label(izcol),unitslabel(izcol))-1)))
        call string_replace(integrate_label,'%uz',trim(unitslabel(izcol)))
     else
        call string_replace(integrate_label,'%z',trim(label(izcol)))
     endif
  else
     if (normalise) then
        integrate_label = '< '//trim(labelin)//' >'
     else
        if (iRescale) then
           integrate_label = '\(2268) '//trim(labelin)//' d'// &
              trim(label(izcol)(1:index(label(izcol),unitslabel(izcol))-1))//trim(labelzintegration)
        else
           integrate_label = '\(2268) '//trim(labelin)//' d'//trim(label(izcol))
        endif
        if (iplot.eq.irho .and. (index(labelin,'density').ne.0 .or. index(labelin,'rho').ne.0)) then
           integrate_label = 'column density'
           !--try to get units label right for column density
           !  would be nice to have a more robust way of knowing what the units mean
           if (iRescale .and. index(labelzintegration,'cm').gt.0  &
                        .and. trim(adjustl(unitslabel(irho))).eq.'[g/cm\u3\d]') then
              integrate_label = trim(integrate_label)//' [g/cm\u2\d]'
           endif
        endif
     endif
  endif
end function integrate_label

end module timestep_plotting
