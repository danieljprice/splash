module timestep_plotting
  implicit none

  integer, parameter, private :: maxtitles = 50
  integer, private :: ninterp
  integer, private :: iplotx,iploty,iplotz,irenderplot,ivectorplot,ivecx,ivecy
  integer, private :: nyplots,npartdim      
  integer, private :: ngrid
  integer, private :: just, ntitles
  integer, private :: iplots,iplotsonpage

  real, dimension(:), allocatable, private :: datpix1D, xgrid
  real, private :: xmin,xmax,ymin,ymax,zmin,ymean
  real, private :: rendermin,rendermax
  real, private :: dxsec,xsecpos
  real, private :: charheight, charheightmm
  real, private :: dxgrid,xmingrid,xmaxgrid
  real, private :: angletempx, angletempy, angletempz

  logical, private :: iplotpart,iplotcont,x_sec,isamexaxis,isameyaxis
  logical, private :: log, inewpage, tile_plots, isave, lastplot
  logical, private :: initialise_xsec
  logical, private :: imulti,iChangeRenderLimits

  character(len=60), dimension(maxtitles), private :: titlelist

contains

!
! initialise plotting options
! called once for all steps
!
subroutine initialise_plotting(ipicky,ipickx,irender)
  use params
  use colours, only:colour_set
  use labels, only:label,ix,ipowerspec
  use limits, only:lim
  use multiplot
  use prompting
  use settings_data, only:ndim,numplot
  use settings_page
  use settings_render, only:icolours,iplotcont_nomulti
  use settings_vecplot, only:iplotpartvec
  use settings_xsecrot, only:xsec_nomulti,xsecpos_nomulti,flythru,nxsec
  use settings_powerspec, only:options_powerspec
  use particle_data, only:npartoftype
  implicit none
  integer, intent(in) :: ipicky,ipickx,irender
  integer :: i,j,ierr
  
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
  iChangeRenderLimits = .false.
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
  
  iplotz = 0
  if (initialise_xsec) then
     !!--work out coordinate that is not being plotted 
     iplotz = 0
     if (ndim.ge.3 .and. x_sec) then
        do j=1,ndim
           if ((iplotx.ne.iploty).and. &
               (j.ne.iplotx).and.(j.ne.iploty)) iplotz = j
        enddo
     endif
     
     if (x_sec .and. iplotz.gt.0) then
!
!--if series of cross sections (flythru), set position of first one
!
        if (flythru) then
           print 32,label(iplotz)
32         format('enter number of ',a1,' cross-section slices')
           read*,nxsec
           !!--dxsec is the distance between slices            
           dxsec = (lim(iplotz,2)-lim(iplotz,1))/float(nxsec)
           xsecpos = lim(iplotz,1) - 0.5*dxsec
           xsecpos_nomulti = xsecpos
        else
!
!--if single cross-section, read position of cross-section slice
!
           call prompt(' enter '//trim(label(iplotz))// &
                       ' position for cross-section slice:', &
                       xsecpos_nomulti,lim(iplotz,1),lim(iplotz,2))
!
!--set thickness if plotting particles
!  (default thickness is half of the average particle spacing)
!
           if (irender.le.0 .or. irender.gt.numplot) then
              npartdim = int(maxval(npartoftype(:,1))**(1./real(ndim)))
              print*,'average # of particles in each dimension = ',npartdim
              if (npartdim.gt.0) then
                 dxsec = (lim(iplotz,2)-lim(iplotz,1))/float(npartdim)
              else
                 dxsec = 0.
              endif
              call prompt(' enter thickness of cross section slice:', &
                           dxsec,0.0,lim(iplotz,2)-lim(iplotz,1))
           elseif (ndim.eq.3) then
!
!--for rendered cross sections in 3D, set thickness to 10%
!  this is the distance slices are moved up and down in interactive mode
!           
              dxsec = 0.1*(lim(iplotz,2)-lim(iplotz,1))
           endif
        endif ! flythru or single
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
  
  !!--set background/foreground colours
  ierr = 0
  if (len_trim(colour_back).gt.0) call pgscrn(0,colour_back,ierr)
  if (ierr /= 0) print 10,'background',trim(colour_back)
  if (len_trim(colour_fore).gt.0) call pgscrn(1,colour_fore,ierr)
  if (ierr /= 0) print 10,'foreground',trim(colour_fore)
10 format(' error: ',a,' colour "',a,'" not found in table')
  
  !!--set colour table
  if (((irender.gt.ndim).or.any(irendermulti(1:nyplots).gt.ndim)) &
       .and.(icolours.gt.0)) then
     call colour_set(icolours)
  endif
  
  !!--set character height in mm
  charheightmm = 4.0
  !!if ((ndown*nacross).gt.1 .and..not. tile_plots) charheight = 2.0
  !      charheight = 0.5*(nacross+ndown)
  
  !!--set line width to something visible
  call pgslw(3)

end subroutine initialise_plotting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plotstep(istep,irender,ivecplot, &
                    npartoftype,dat,timei,gammai,ipagechange,iadvance)
  use params
  use exact, only:exact_solution, &
             atstar,ctstar,sigma,iwaveplotx,iwaveploty
  use toystar1D, only:exact_toystar_ACplane
  use labels
  use limits
  use multiplot
  use particle_data, only:maxpart,icolourme
  use rotation
  use settings_data, only:numplot,ndataplots,icoords,ndim,ndimv,nstart,n_end,nfreq
  use settings_limits
  use settings_part, only:icoordsnew,iexact,iplotlinein,linestylein, &
                     iplotline,iplotpartoftype
  use settings_page, only:nacross,ndown,iadapt,interactive,iaxis, &
                     hpostitle,vpostitle,fjusttitle
  use settings_render, only:npix,ncontours,icolours,iplotcont_nomulti, &
                       iPlotColourBar,icolour_particles
  use settings_vecplot, only:npixvec, iplotpartvec
  use settings_xsecrot
  use settings_powerspec
!
!--subroutines called from this routine
!
  use colourparts
  use transforms
  use interactive_routines
  use geometry
  use legends, only:legend
  use particleplots
  use powerspectrums
  use interpolations1D, only:interpolate1D
  use interpolations2D, only:interpolate2D, interpolate2D_xsec
  use projections3D, only:interpolate3D_projection
  use xsections3D, only:interpolate3D, interpolate3D_fastxsec, &
                        interpolate3D_xsec_vec
  use render, only:render_pix,colourbar

  implicit none
  integer, intent(in) :: istep, irender, ivecplot
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(:,:), intent(in) :: dat
  real, intent(in) :: timei,gammai
  logical, intent(in) :: ipagechange
  integer, intent(inout) :: iadvance
  
  integer :: ntoti,iz
  integer :: j,k
  integer :: nyplot
  integer :: irenderpart,irendered
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: npixyvec
  integer :: index1,index2,itype

  real, parameter :: pi = 3.1415926536
  real, parameter :: tol = 1.e-10 ! used to compare real numbers
  real, dimension(:,:), allocatable :: datpix
  real, dimension(:,:,:), allocatable :: datpix3D
  real, dimension(ndim) :: xcoords,vecnew,xmintemp,xmaxtemp
  real, dimension(max(maxpart,2000)) :: xplot,yplot,zplot,renderplot
  real :: angleradx, anglerady, angleradz
  real :: rendermintemp,rendermaxtemp
  real :: xsecmin,xsecmax,dummy
  real :: pixwidth

  character(len=len(label(1))+20) :: labelx,labely,labelz,labelrender,labelvecplot
  character(len=120) :: title
  character(len=5) :: string
  
34   format (25(' -'))

  !--set labels to blank (just in case)
  labelx = ' '
  labely = ' '
  labelz = ' '
  labelrender = ' '
  labelvecplot = ' '
  xplot = 0.
  yplot = 0.
  dummy = 0.
  !
  !--set number of particles to use in the interpolation routines
  !  (ie. including only gas particles and ghosts)
  !--if plotting ghost particles, set ntotplot = ntot, else ntot=npart
  !
  ntoti = sum(npartoftype)
  ninterp = npartoftype(1)
  if (labeltype(2)(1:5).eq.'ghost' .and. iplotpartoftype(2)) then
     ninterp = ninterp + npartoftype(2)
  endif
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
        if (icolour_particles) then
           irenderpart = irendermulti(nyplot)
           irenderplot = 0
        else
           irenderpart = 0
           irenderplot = irendermulti(nyplot)
        endif
        ivectorplot = ivecplotmulti(nyplot)
        iplotcont = iplotcontmulti(nyplot)
        x_sec = x_secmulti(nyplot)
        xsecpos = xsecposmulti(nyplot)
     else
        if (icolour_particles) then
           irenderpart = irender
           irenderplot = 0
        else
           irenderpart = 0
           irenderplot = irender
        endif
        ivectorplot = ivecplot
        iplotcont = iplotcont_nomulti
        x_sec = xsec_nomulti
        if (iadvance.ne.0) xsecpos = xsecpos_nomulti        
     endif
     if (ivectorplot.gt.0) iplotpart = iplotpartvec

     !--------------------------------------------------------------
     !  copy from main dat array into xplot, yplot 
     !  also set labels and plot limits
     !--------------------------------------------------------------
     if (iploty.le.ndataplots .and. iplotx.le.ndataplots) then
        xplot(1:ntoti) = dat(1:ntoti,iplotx)
        yplot(1:ntoti) = dat(1:ntoti,iploty)
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
        !
        !--change coordinate system if relevant
        !
        if (icoordsnew.ne.icoords) then
           !--do this if one is a coord but not if rendering
           if ((iplotx.le.ndim .or. iploty.le.ndim) &
                .and..not.(iplotx.le.ndim.and.iploty.le.ndim &
                .and.irenderplot.gt.ndim)) then
              print*,'changing coords from ',trim(labelcoordsys(icoords)), &
                     ' to ',trim(labelcoordsys(icoordsnew))
              do j=1,ntoti
                 call coord_transform(dat(j,ix(1:ndim)),ndim,icoords, &
                                      xcoords(1:ndim),ndim,icoordsnew)
                 if (iplotx.le.ndim) xplot(j) = xcoords(iplotx)
                 if (iploty.le.ndim) yplot(j) = xcoords(iploty)
              enddo
              if (iadvance.ne.0) then
                 xmintemp = lim(ix(1:ndim),1)
                 xmaxtemp = lim(ix(1:ndim),2)
                 call coord_transform_limits(xmintemp,xmaxtemp, &
                                             icoords,icoordsnew,ndim)
                 if (iplotx.le.ndim) then
                    xmin = xmintemp(iplotx)
                    xmax = xmaxtemp(iplotx)
                 endif
                 if (iploty.le.ndim) then
                    ymin = xmintemp(iploty)
                    ymax = xmaxtemp(iploty)
                 endif
              endif
           endif
           if (iamvec(iplotx).gt.0) then
              if (iplotx-iamvec(iplotx)+1 .le. ndim) then
                 print*,'changing vector component from ', &
                  trim(labelcoordsys(icoords)),' to ',trim(labelcoordsys(icoordsnew))
                 do j=1,ntoti
                    call vector_transform(dat(j,ix(1:ndim)), &
                         dat(j,iamvec(iplotx):iamvec(iplotx)+ndim-1), &
                         ndim,icoords,vecnew(1:ndim),ndim,icoordsnew)
                    xplot(j) = vecnew(iplotx-iamvec(iplotx)+1)
                 enddo
              else
                 print*,'error: can''t convert vector components with ndimV > ndim'
              endif
           endif
           if (iamvec(iploty).gt.0) then
              if (iploty-iamvec(iploty)+1 .le.ndim) then
                 print*,'changing vector component from ', &
                  trim(labelcoordsys(icoords)),' to ',trim(labelcoordsys(icoordsnew))
                 do j=1,ntoti
                    call vector_transform(dat(j,ix(1:ndim)), &
                         dat(j,iamvec(iploty):iamvec(iploty)+ndim-1), &
                         ndim,icoords,vecnew(1:ndim),ndim,icoordsnew)
                    yplot(j) = vecnew(iploty-iamvec(iploty)+1)
                 enddo
              else
                 print*,'error: can''t convert vector components with ndimV > ndim'
              endif
           endif
        endif
        !--apply transformations (log, 1/x etc) if appropriate
        !  also change labels and limits appropriately
        if (.not.(iplotx.le.ndim .and. iploty.le.ndim)) then
           if (itrans(iplotx).ne.0) then
              call transform(xplot(1:ntoti),itrans(iplotx))
              labelx = transform_label(labelx,itrans(iplotx))
              if (iadvance.ne.0) call transform_limits(xmin,xmax,itrans(iplotx))
           endif
           if (itrans(iploty).ne.0) then
              call transform(yplot(1:ntoti),itrans(iploty))
              labely = transform_label(labely,itrans(iploty))
              if (iadvance.ne.0) call transform_limits(ymin,ymax,itrans(iploty))
           endif
        endif

        !--write username, date on plot
        !         if (nacross.le.2.and.ndown.le.2) call pgiden

        !--adjust plot limits if adaptive plot limits set
        if (ipagechange .and. iadapt .and. iadvance.ne.0) then
           xmin = huge(xmin)
           xmax = -huge(xmax)
           ymin = huge(ymin)
           ymax = -huge(ymax)
           !--find maximum over all particle types being plotted
           index1 = 1
           print*,'adapting plot limits'
           do itype=1,maxparttypes
              index2 = index1 + npartoftype(itype) - 1
              if (iplotpartoftype(itype).and.npartoftype(itype).gt.0) then
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
        iz = 0
        if (ndim.ge.3) then
           do j=1,ndim
              if ((iplotx.ne.iploty).and. &
                  (j.ne.iplotx).and.(j.ne.iploty)) iz = j
           enddo
        endif
        
        iplotz = 0
        if (x_sec) iplotz = iz ! this is used as cross sectioned quantity
        if (iplotz.gt.0 .and. iplotz.le.ndataplots) then
           zplot(1:ntoti) = dat(1:ntoti,iplotz)
           labelz = label(iplotz)
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
            anglerady = angletempy*pi/180.
            angleradx = angletempx*pi/180.
            print "(a,f6.2)",'rotating particles about z by ',angletempz
            if (ndim.eq.3) then
               print "(a,f6.2)",'rotating particles about y by ',angletempy
               print "(a,f6.2)",'rotating particles about x by ',angletempx
            endif
            do j=1,ntoti
               xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
!               if (ndim.eq.2) then
!                  call rotate2D(xcoords(:),angleradz,anglerady)
!               elseif (ndim.eq.3) then
                  call rotate3D(xcoords(:),angleradx,anglerady,angleradz)
!               endif
               xplot(j) = xcoords(iplotx) + xorigin(iplotx)
               yplot(j) = xcoords(iploty) + xorigin(iploty)
               if (iplotz.gt.0) then
                  zplot(j) = xcoords(iplotz) + xorigin(iplotz)
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
           !--interpolate from particles to fixed grid using SPH summation
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

           !!--only need z pixels if working with interpolation to 3D grid
           !  (then number of z pixels is equal to number of cross sections)
           if ((ndim.ge.3).and.(x_sec.and.nxsec.gt.2)) then
              zmin = lim(iplotz,1)
              npixz = nxsec
           endif

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
                      dat(1:ninterp,ipmass),dat(1:ninterp,irho), &
                      dat(1:ninterp,ih),dat(1:ninterp,irenderplot), &
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
                      zplot(1:ninterp),dat(1:ninterp,ipmass),  &
                      dat(1:ninterp,irho),dat(1:ninterp,ih), &
                      dat(1:ninterp,irenderplot), &
                      ninterp,xmin,ymin,zmin,datpix3D,npixx,npixy,npixz,pixwidth,dxsec)
              endif
           end select

        endif
        !
        !--if vector plot determine whether or not to plot the particles as well
        !
        iplotpart = .true.
        if (ivectorplot.gt.0) iplotpart = iplotpartvec
        if (irenderplot.gt.0) iplotpart = .false.
        !
        !--this is a flag to say whether or not rendered limits have been changed
        !  interactively. False by default, but must retain value whilst
        !  iadvance = 0
        !
        if (iadvance.ne.0) iChangeRenderLimits = .false.

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
                 print*,TRIM(label(iplotz)),' = ',xsecpos, &
                      ' cross section, pixel ',ipixxsec
                 datpix = datpix3D(:,:,ipixxsec)    ! slices are in 3rd dimension

              else
                 !-------------------------------------------------------------------
                 !  or do a fast projection/cross section of 3D data to 2D array
                 !-------------------------------------------------------------------
                 if (x_sec) then
                    !!--do fast cross-section
                    print*,trim(label(ix(iplotz))),' = ',xsecpos,  &
                         ' : fast cross section', xmin,ymin
                    call interpolate3D_fastxsec( &
                         xplot(1:ninterp),yplot(1:ninterp), &
                         zplot(1:ninterp), &
                         dat(1:ninterp,ipmass),dat(1:ninterp,irho),    &
                         dat(1:ninterp,ih),dat(1:ninterp,irenderplot), &
                         ninterp,xmin,ymin,xsecpos,datpix,npixx,npixy,pixwidth)
                 else
                    !!--do fast projection
                    call interpolate3D_projection( &
                         xplot(1:ninterp),yplot(1:ninterp), &
                         dat(1:ninterp,ipmass),dat(1:ninterp,irho),   &
                         dat(1:ninterp,ih), dat(1:ninterp,irenderplot), &
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
                   dat(1:ninterp,iplotx),dat(1:ninterp,iploty), &
                   dat(1:ninterp,ipmass),dat(1:ninterp,irho),    &
                   dat(1:ninterp,ih),dat(1:ninterp,irenderplot), &
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
           ! setup page
           !---------------------------------

           just = 1  ! x and y axis have same scale
           ! unless 1D xsec through 2D data or non-cartesian
           if ((irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) &
               .or.(icoordsnew.ne.1)) then
              just = 0 
           endif
           title = ' '

           call page_setup

           !--add to log
           if (x_sec.and.iplotpart.and.iplotz.gt.0) print 35,label(iplotz),xsecmin,label(iplotz),xsecmax
35            format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

           !------------------------------
           ! now actually plot the data
           !------------------------------
           if (irenderplot.gt.ndim) then
              if ((ndim.eq.3).or.(ndim.eq.2.and. .not.x_sec)) then
                 !---------------------------------------------------------------
                 ! scalar quantity which has been rendered to a 2D pixel array (datpix)
                 !---------------------------------------------------------------
                 !!--do transformations on rendered array  
                 call transform2(datpix,itrans(irenderplot))
                 labelrender = label(irenderplot)
                 !!--set label for column density (projection) plots (2268 or 2412 for integral sign)
                 if (ndim.eq.3 .and..not. x_sec) then
                    labelrender = '\(2268) '//trim(labelrender)//' d'//trim(label(ix(iz)))
                 endif
                 !!--apply transformations to the label for the rendered quantity 
                 !!  but don't do this for log as we use a logarithmic axis instead
                 log = .false.
                 !if (itrans(irenderplot).eq.1) then
                 !   log = .true.
                 !else 
                    labelrender = transform_label(labelrender,itrans(irenderplot))
                 !endif
                 !!--limits for rendered quantity
                 if (iadvance.ne.0 .or. .not.iChangeRenderLimits) then
                    if (iadapt) then
                       !!--if adaptive limits, find limits of rendered array
                       rendermin = minval(datpix)
                       rendermax = maxval(datpix)
                       print*,'adapting render limits'
                    else
                       !!--or apply transformations to fixed limits
                       rendermin = lim(irenderplot,1)
                       rendermax = lim(irenderplot,2)
                       call transform_limits(rendermin,rendermax,itrans(irenderplot))
                    endif
                 endif
                 
                 !!--if log, then set zero values to minimum
                 !!  (must be done after minimum is known)
                 !!
                 write(string,*) itrans(irenderplot)
                 if (index(string,'1').ne.0) then
                    !!print*,'setting zero values to ',rendermin,' on log array'
                    where (abs(datpix).lt.tiny(datpix))
                       datpix = rendermin
                    end where
                    !!  also do not let max=0 as this is suspiciously wrong
                    if (iadapt .and. abs(rendermax).lt.tiny(datpix)) then
                       !!print*,'max=0 on log plot, fixing'
                       rendermax = maxval(datpix)
                    endif
                 endif

                 !!--print plot limits to screen
                 print*,trim(labelrender),' min, max = ',rendermin,rendermax
                      
                 !!--call subroutine to actually render the image       
                 call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                      npixx,npixy,xmin,ymin,pixwidth,    &
                      icolours,iplotcont,ncontours,log)

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
                 !
                 !--sort out particle colouring
                 !
                 if (irenderpart.gt.0 .and. irenderpart.le.numplot) then
                    renderplot(1:ntoti) = dat(1:ntoti,irenderpart)
                    call transform(renderplot(1:ntoti),itrans(irenderpart))
                    labelrender = label(irenderpart)
                    labelrender = transform_label(labelrender,itrans(irenderpart))
                    
                    !!--limits for rendered quantity
                    if (iadvance.ne.0 .or. .not.iChangeRenderLimits) then
                       if (iadapt) then
                          !!--if adaptive limits, find limits of rendered array
                          rendermin = minval(renderplot(1:ntoti))
                          rendermax = maxval(renderplot(1:ntoti))
                       else
                          !!--or apply transformations to fixed limits
                          rendermin = lim(irenderpart,1)
                          rendermax = lim(irenderpart,2)
                          call transform_limits(rendermin,rendermax,itrans(irenderpart))
                       endif
                    endif
                    !!--print plot limits to screen
                    print*,trim(labelrender),' min, max = ',rendermin,rendermax       

                    call colour_particles(renderplot(1:ntoti), &
                         rendermin,rendermax, &
                         icolourme(1:ntoti),ntoti)
                    !!--plot colour bar
                    if (iPlotColourBar) call colourbar(icolours,rendermin,rendermax, &
                                                       trim(labelrender),.false.)
                 endif
                 !
                 !--do particle plot
                 !
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),dat(1:ntoti,ih),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),npartoftype(:), &
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
                pixwidth = (xmax-xmin)/real(npixvec)
                npixyvec = int((ymax-ymin)/pixwidth) + 1

                call vector_plot(ivecx,ivecy,npixvec,npixyvec,pixwidth,labelvecplot)
             endif
           endif
           !---------------------------------
           ! plot rotated axes
           !---------------------------------
           if (irotate .and. irotateaxes.gt.0) then
              if (ndim.eq.3) then
                 call rotate_axes3D(irotateaxes,iplotx,iploty, &
                      xminrotaxes(1:ndim),xmaxrotaxes(1:ndim),xorigin(1:ndim), &
                      angleradx,anglerady,angleradz)
              elseif (ndim.eq.2) then
                 call rotate_axes2D(irotateaxes,xminrotaxes(1:ndim), &
                                   xmaxrotaxes(1:ndim),xorigin(1:ndim),angleradz)
              endif
           endif
           !
           !--print legend if this is the first plot on the page
           !    
           if (nyplot.eq.1) then
              call legend(timei)
           endif
           !
           !--print title if appropriate
           !
           if (istep.le.ntitles) then
              if (titlelist(istep)(1:1).ne.' ') then
                 call pgmtxt('T',vpostitle,hpostitle,fjusttitle,trim(titlelist(istep)))
              endif
           endif
           !
           !--plot exact solution if relevant (before going interactive)
           !
           if (iexact.ne.0) then
               call exact_solution(iexact,iplotx,iploty, &
                    itrans(iplotx),itrans(iploty),icoordsnew, &
                    ndim,ndimV,timei,xmin,xmax,0.0,gammai, &
                    dat(1:npartoftype(1),ipmass),npartoftype(1))
           endif
           !
           !--enter interactive mode
           !
           lastplot = (istep.eq.n_end .and. nyplot.eq.nyplots .and. k.eq.nxsec)

           if (interactive) then
              if (nacross*ndown.eq.1) then
                 iadvance = nfreq
                 rendermintemp = rendermin
                 rendermaxtemp = rendermax
                 if (icolour_particles) then
                    irendered = irenderpart
                 else
                    irendered = irenderplot
                 endif
                 call interactive_part(ninterp,iplotx,iploty,iplotz,irendered, &
                      xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                      dat(1:ninterp,ih),icolourme(1:ninterp), &
                      xmin,xmax,ymin,ymax,xsecpos,dxsec,rendermin,rendermax, &
                      angletempx,angletempy,angletempz,ndim,iadvance,isave)
                 !--turn rotation on if necessary
                 if (abs(angletempx-anglex).gt.tol) irotate = .true.
                 if (abs(angletempy-angley).gt.tol) irotate = .true.
                 if (abs(angletempz-anglez).gt.tol) irotate = .true.
                 if (abs(rendermintemp-rendermin).gt.tol) iChangeRenderLimits = .true.
                 if (abs(rendermaxtemp-rendermax).gt.tol) iChangeRenderLimits = .true.
                 if (iadvance.eq.-666) return
              elseif (iplotsonpage.eq.nacross*ndown .or. lastplot) then
                 !
                 !--timestep control only if multiple plots on page
                 !
                 iadvance = nfreq
                 call interactive_step(iadvance,xmin,xmax,ymin,ymax)
                 if (iadvance.eq.-666) return
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

        !--------------------------------
        ! setup page
        !--------------------------------
        just = 0
        title = ' '
        call page_setup

        !--------------------------------
        ! now plot particles
        !--------------------------------

        !--plot time on plot
        if (nyplot.eq.1) call legend(timei)
        !
        !--sort out particle colouring
        !
        if (irenderpart.gt.0 .and. irenderpart.le.numplot) then
           renderplot(1:ntoti) = dat(1:ntoti,irenderpart)
           call transform(renderplot(1:ntoti),itrans(irenderpart))

           !!--limits for rendered quantity
           if (iadapt .and. iadvance.ne.0) then
              !!--if adaptive limits, find limits of rendered array
              rendermin = minval(renderplot(1:ntoti))
              rendermax = maxval(renderplot(1:ntoti))
           elseif (iadvance.ne.0) then                   
              !!--or apply transformations to fixed limits
              rendermin = lim(irenderpart,1)
              rendermax = lim(irenderpart,2)
              call transform_limits(rendermin,rendermax,itrans(irenderpart))
           endif
           !!--print plot limits to screen
           print*,trim(labelrender),' min, max = ',rendermin,rendermax       

           call colour_particles(renderplot(1:ntoti), &
                rendermin,rendermax, &
                icolourme(1:ntoti),ntoti)
        endif

        call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
             zplot(1:ntoti),dat(1:ntoti,ih),ntoti,iplotx,iploty, &
             icolourme(1:ntoti),npartoftype(:),.false.,0.0,0.0,' ')

        if ((istep.eq.nstart).and.iplotlinein) then! plot initial conditions as dotted line
           call pgsls(linestylein)
        endif

        !--plot line joining the particles
        if (iplotline.or.(iplotlinein.and.(istep.eq.nstart))) then
           call pgline(npartoftype(1),xplot(1:npartoftype(1)), &
                       yplot(1:npartoftype(1)))     
        endif
        call pgsls(1)! reset 
        call pgsch(charheight)
        !
        !--plot exact solution if relevant
        !
        if (iexact.eq.5 .and.(iploty.eq.iwaveploty).and.(iplotx.eq.iwaveplotx)) then 
           ymean = SUM(yplot(1:npartoftype(1)))/REAL(npartoftype(1)) 
        else
           ymean = 0.
        endif

        if (iexact.ne.0 .or. (iploty.eq.irho .and. iplotx.eq.ih) .or. &
           (iplotx.eq.irho .and. iploty.eq.ih)) then
           call exact_solution(iexact,iplotx,iploty,itrans(iplotx),itrans(iploty), &
                icoordsnew,ndim,ndimV,timei,xmin,xmax,ymean,gammai, &
                dat(1:npartoftype(1),ipmass),npartoftype(1))
        endif
        !
        !--enter interactive mode
        !
        lastplot = (istep.eq.n_end .and. nyplot.eq.nyplots)

        if (interactive) then
           if (nacross*ndown.eq.1) then
              iadvance = nfreq
              call interactive_part(ntoti,iplotx,iploty,0,irenderpart, &
                   xplot(1:ntoti),yplot(1:ntoti),zplot(1:ntoti), &
                   dat(1:ntoti,ih),icolourme(1:ntoti), &
                   xmin,xmax,ymin,ymax,dummy,dummy,rendermin,rendermax, &
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
           call exact_toystar_acplane(atstar,ctstar,sigma,gammai)
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
              ngrid = 2*npartoftype(1)
              !!--set up 1D grid
              xmingrid = lim(ix(1),1)
              xmaxgrid = lim(ix(1),2)
              dxgrid = (xmaxgrid-xmingrid)/ngrid
              call set_grid1D(xmingrid,dxgrid,ngrid)

              ninterp = ntoti
              !!--interpolate to 1D grid  
              call interpolate1D(dat(1:ninterp,ix(1)), & 
                   dat(1:ninterp,ipmass),dat(1:ninterp,irho), &
                   dat(1:ninterp,ih),dat(1:ninterp,ipowerspecy), & 
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
              call powerspectrum_lomb(ntoti,dat(1:ntoti,ix(1)), &
                   dat(1:ntoti,ipowerspecy),nfreqspec, &
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
              call transform(xplot(1:nfreqspec),itrans(iploty))
              labelx = transform_label(labelx,itrans(iploty))

              call transform(yplot(1:nfreqspec),itrans(iploty))
              labely = transform_label(labely,itrans(iploty))
              if (iadvance.ne.0) then
                 call transform_limits(xmin,xmax,itrans(iploty))
                 call transform_limits(ymin,ymax,itrans(iploty))
              endif
           endif
           
           just = 0
           title = 'Power Spectrum'
           call page_setup

           call pgline(nfreqspec,xplot(1:nfreqspec),yplot(1:nfreqspec))

        endif
        !
        !--if this is the first plot on the page, print legend
        !
        if (iplotsonpage.eq.1) call legend(timei)

        lastplot = (istep.eq.n_end)

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
  
contains

!----------------------------------------------
! interfaces to the page setup routines
! this is called just before a plot is
! actually plotted
!----------------------------------------------
  subroutine page_setup
    implicit none
    
    !--------------------------------------------------------------
    ! output some muff to the screen
    !--------------------------------------------------------------

    if (interactive) then
       print*,trim(labelx),'min,max = ',xmin,xmax
       print*,trim(labely),'min,max = ',ymin,ymax
    endif

    !---------------------
    ! increment counters
    !---------------------

    iplots = iplots + 1
    iplotsonpage = iplotsonpage + 1
    if (iplotsonpage.gt.nacross*ndown) iplotsonpage = 1

    !--------------------------------------------------------------
    ! set up pgplot page
    !--------------------------------------------------------------

    if (tile_plots) then
       if (iplotsonpage.eq.1 .and. ipagechange) call pgpage
       call danpgtile(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
                      trim(labelx),trim(labely),trim(title),just,iaxis)
    else
        !--change the page if pagechange set
       !  or, if turned off, between plots on first page only
       inewpage = ipagechange .or. (iplots.le.nacross*ndown)
       call setpage(iplotsonpage,nacross,ndown,xmin,xmax,ymin,ymax, &
         trim(labelx),trim(labely),trim(title), &
         just,iaxis,isamexaxis,isameyaxis,inewpage)
    endif 

    return
  end subroutine page_setup

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
! interface to vector plotting routines
! so that pixel arrays are allocated appropriately
!-------------------------------------------------------------------
  subroutine vector_plot(ivecx,ivecy,numpixx,numpixy,pixwidth,label)
   use fieldlines
   use settings_vecplot
   use interpolations2D, only:interpolate2D_vec
   use projections3D, only:interpolate3D_proj_vec
   use render, only:render_vec
   implicit none
   integer, intent(in) :: ivecx,ivecy,numpixx,numpixy
   real, intent(in) :: pixwidth
   character(len=*), intent(in) :: label
   real, dimension(numpixx,numpixy) :: vecpixx, vecpixy
   real :: vecmax

   print*,'plotting vector field ',trim(label)
   if ((ivecx.le.ndim).or.(ivecx.gt.ndataplots) &
        .or.(ivecy.le.ndim).or.(ivecy.gt.ndataplots)) then
      print*,'error finding location of vector plot in array'
   else
      if (iadapt) then
         vecmax = -1.0  ! plot limits then set in vectorplot
      else                    
         vecmax = max(lim(ivecx,2),lim(ivecy,2))
      endif

      !!--plot arrows in either background or foreground colour
      if (UseBackgndColorVecplot) call pgsci(0)

      !
      !--interpolate using appropriate routine for number of dimensions
      !
      select case(ndim)
      case(3)
         if (x_sec) then ! take vector plot in cross section
            call interpolate3D_xsec_vec(xplot(1:ninterp), &
              yplot(1:ninterp),zplot(1:ninterp), &
              dat(1:ninterp,ipmass),dat(1:ninterp,irho),  &
              dat(1:ninterp,ih),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
              ninterp,xmin,ymin,xsecpos, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidth)
         else
            call interpolate3D_proj_vec(xplot(1:ninterp), &
              yplot(1:ninterp),dat(1:ninterp,ipmass), &
              dat(1:ninterp,irho),dat(1:ninterp,ih), &
              dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
              ninterp,xmin,ymin, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidth)
         endif
      case(2)
         !
         !--or interpolate (via averaging) to coarser grid
         !
         !call fieldlines2D(ninterp,xplot(1:ninterp),yplot(1:ninterp), &
         !     dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
         !     dat(1:ninterp,ih),dat(1:ninterp,ipmass), &
         !     dat(1:ninterp,irho),xmin,xmax,ymin,ymax)
         !call interpolate_vec(xplot(1:ninterp),yplot(1:ninterp), &
         !  dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
         !  xmin,ymin,pixwidth,vecpixx,vecpixy, &
         !  ninterp,numpixx,numpixy)
         
         call interpolate2D_vec(xplot(1:ninterp),yplot(1:ninterp), &
              dat(1:ninterp,ipmass),dat(1:ninterp,irho), &
              dat(1:ninterp,ih),dat(1:ninterp,ivecx), &
              dat(1:ninterp,ivecy),ninterp,xmin,ymin, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidth)
      
      case default
         print "(a,i1,a)",'ERROR: Cannot do vector plotting in ',ndim,' dimensions'
         return
      end select
      !
      !--plot it
      !
      call render_vec(vecpixx,vecpixy,vecmax, &
           numpixx,numpixy,xmin,ymin,pixwidth,trim(label))

      if (UseBackgndColorVecplot) call pgsci(1)

   endif
  
  end subroutine vector_plot

end subroutine plotstep
  

end module timestep_plotting
