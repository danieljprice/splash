module interactive_routines
 use colourbar, only:barisvertical,incolourbar,incolourbarlabel
 implicit none
 public :: interactive_part,interactive_step,interactive_multi
 private :: mvlegend,mvtitle,save_limits,save_rotation
 private :: get_vptxy
 real, private :: xcursor = 0.5
 real, private :: ycursor = 0.5
 
 private
 
contains
!
!--interactive tools on particle plots
!  allows user to change settings interactively
!
!  Arguments:
!
!  INPUT:
!   npart   : number of particles plotted
!   iplotx  : quantity plotted as x axis
!   iploty  : quantity plotted as y axis 
!   iplotz  : quantity to use in selecting particles
!   irender : quantity rendered
!   xcoords(npart) : x coordinates of particles
!   ycoords(npart) : y coordinates of particles
!   zcoords(npart) : z coordinates (or third quantity) of particles
!   hi(npart)      : smoothing lengths of particles
!   zmin, zmax     : range of z within which to plot particles
!   istep          : current step position
!   ilaststep      : position of last timestep
!
! CHANGEABLE:
!   icolourpart(npart) : flag indicating colour of particles
!   xmin, xmax, ymin, ymax : current plot limits
!   rendermin, rendermax   : current rendering limits
!   vecmax : maximum vector limits
!
! OUTPUT:
!   iadvance : integer telling the loop how to advance the timestep
!   irerender : if set, redo rendering. Anything which requires rendering
!               to be recalculated must set this
!
subroutine interactive_part(npart,iplotx,iploty,iplotz,irender,icontour,ivecx,ivecy, &
  xcoords,ycoords,zcoords,hi,icolourpart,xmin,xmax,ymin,ymax, &
  rendermin,rendermax,renderminadapt,rendermaxadapt,contmin,contmax,vecmax, &
  anglex,angley,anglez,ndim,xorigin,x_sec,zslicepos,dzslice, &
  zobserver,dscreen,use3Dopacity,taupartdepth,irerender,itrackpart,icolourscheme, &
  iColourBarStyle,labelrender,iadvance,istep,ilaststep,iframe,nframes,interactivereplot)
  use settings_xsecrot, only:setsequenceend
  use shapes, only:inshape,edit_shape,edit_textbox,delete_shape
  use multiplot, only:itrans
  use settings_render, only:projlabelformat,iapplyprojformat
  implicit none
  integer, intent(in) :: npart,irender,icontour,ndim,iplotz,ivecx,ivecy,istep,ilaststep,iframe,nframes
  integer, intent(inout) :: iColourBarStyle
  integer, intent(inout) :: iplotx,iploty,itrackpart,icolourscheme
  integer, intent(out) :: iadvance
  integer, dimension(npart), intent(inout) :: icolourpart
  real, dimension(npart), intent(in) :: xcoords,ycoords,zcoords,hi
  real, intent(inout) :: xmin,xmax,ymin,ymax,rendermin,rendermax,vecmax,contmin,contmax,taupartdepth
  real, intent(inout) :: anglex,angley,anglez,zslicepos,dzslice,zobserver,dscreen
  real, intent(in) :: renderminadapt,rendermaxadapt
  real, intent(in), dimension(ndim) :: xorigin
  character(len=*), intent(inout) :: labelrender
  logical, intent(inout) :: x_sec
  logical, intent(out) :: irerender,interactivereplot
  logical, intent(in) :: use3Dopacity
  real, parameter :: pi=3.141592653589
  integer :: i,iclosest,nc,ierr,ixsec,ishape
  integer :: nmarked,ncircpart,itrackparttemp
  integer, dimension(1000) :: icircpart
  real :: xpt,ypt
  real :: xpt2,ypt2,xcen,ycen,xminwin,xmaxwin,yminwin,ymaxwin
  real :: xptmin,xptmax,yptmin,yptmax,zptmin,zptmax
  real :: rmin,rr,gradient,yint,dx,dy,dr,anglerad
  real :: xlength,ylength,renderlength,renderpt,drender,zoomfac
  real :: dxlength,dylength,xmaxin,ymaxin
  real, dimension(4) :: xline,yline
  character(len=1) :: char,char2
  character(len=20) :: string
  logical :: iexit, rotation, verticalbar, iamincolourbar

  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  xline = 0.
  yline = 0.
  !
  !--convert saved cursor position (saved in viewport coords)
  !  back to coordinates
  !
  call pgqwin(xminwin,xmaxwin,yminwin,ymaxwin)
  call get_posxy(xcursor,ycursor,xpt,ypt,xminwin,xmaxwin,yminwin,ymaxwin)
  
!  xpt = 0.
!  ypt = 0.
  xpt2 = 0.
  ypt2 = 0.
  xmaxin = xmax
  ymaxin = ymax
  zoomfac = 1.0
  nc = 0
  ncircpart = 0
  itrackparttemp = itrackpart
  iexit = .false.
  rotation = .false.
  irerender = .false.
  interactivereplot = .false.
  if (iplotx.le.ndim .and. iploty.le.ndim .and. ndim.ge.2) rotation = .true.
  verticalbar = barisvertical(iColourBarStyle)
  
  if (iplotz.gt.0 .and. x_sec) then
     zptmin = zslicepos - 0.5*dzslice
     zptmax = zslicepos + 0.5*dzslice
  else
  !--if not using z range, make it encompass all the particles
     zptmin = -huge(zptmin)
     zptmax = huge(zptmax)
  endif
 
  interactiveloop: do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
     !
     !--find closest particle
     !  
     rmin = 1.e6
     iclosest = 0
     xlength = xmax - xmin
     ylength = ymax - ymin
     if (xlength.gt.tiny(xlength)) then
        dxlength = 1./xlength
     else
        dxlength = 0.
     endif
     if (ylength.gt.tiny(ylength)) then
        dylength = 1./ylength
     else
        dylength = 0.
     endif
          
     do i=1,npart
        !rr = (xcoords(i)-xpt)**2 + (ycoords(i)-ypt)**2
        !--use distance normalised on screen
        rr = ((xcoords(i)-xpt)*dxlength)**2 + ((ycoords(i)-ypt)*dylength)**2
        if (rr.lt.rmin) then
           iclosest = i
           rmin = rr
        endif
     enddo

     !--query the position of the colour bar
     iamincolourbar = incolourbar(iColourBarStyle,xpt,ypt,xmin,xmax,ymin,ymax)

     select case(char)
     !
     !--particle plot stuff
     !
     case('p')
        if (iclosest.gt.0 .and. iclosest.le.npart) then
           print*,' closest particle = ',iclosest,'x = ',xcoords(iclosest),' y =',ycoords(iclosest)
           call plot_number(iclosest,xcoords(iclosest),ycoords(iclosest))
        else
           print*,'error: could not determine closest particle'
        endif
     case('c')
        if (iclosest.gt.0 .and. iclosest.le.npart) then
           print*,'plotting circle of interaction on particle ',iclosest, &
                  ' h = ',hi(iclosest)
           !--save settings for these
           ncircpart = ncircpart + 1
           if (ncircpart.gt.size(icircpart)) then
              print*,'WARNING: ncircles > array limits, cannot save'
              ncircpart = size(icircpart)
           else
              icircpart(ncircpart) = iclosest
           endif
           call pgsfs(2)
           call pgcirc(xcoords(iclosest),ycoords(iclosest),2.*hi(iclosest))
        else
           print*,'error: could not determine closest particle'
        endif
     case('t')
     !--track closest particle (must save to activate)
        if (itrackpart.ne.0 .and. itrackparttemp.eq.itrackpart) then
           itrackpart = 0
           itrackparttemp = 0
           print*,' particle tracking limits OFF'
        else
           if (iclosest.gt.0 .and. iclosest.le.npart) then
              itrackparttemp = iclosest
              call plot_number(iclosest,xcoords(iclosest),ycoords(iclosest))
              print*,' limits set to track particle ',itrackparttemp,' x, y = ', &
                     xcoords(iclosest),ycoords(iclosest)
              print*,' save settings to activate '
           else
              print*,'error: could not determine closest particle'
           endif
        endif
     case('g')   ! draw a line between two points
        xline(2) = xpt
        yline(2) = ypt
        !--mark first point
        call pgpt1(xpt,ypt,4)
        !--select second point
        print*,' select another point (using left click or g) to plot line '
        call pgband(1,1,xline(2),yline(2),xline(3),yline(3),char2)
        !--draw line if left click or g
        select case(char2)
        case('A','g')
           print*
           !--mark second point
           call pgpt1(xline(3),yline(3),4)
           xlength = xline(3)-xline(2)
           if (abs(xlength).lt.tiny(xlength)) then
              xline(1) = xline(2)
              xline(4) = xline(2)
              yline(1) = ymin
              yline(4) = ymax
              print*,' error: gradient = infinite'
           elseif (xline(2).lt.xline(3)) then 
              xline(1) = xmin
              xline(4) = xmax
           else
              xline(1) = xmax
              xline(4) = xmin
           endif           
           ylength = yline(3)-yline(2)
           dr = sqrt(xlength**2 + ylength**2)
           print*,' (x1,y1) = (',xline(2),',',yline(2),')'
           print*,' (x2,y2) = (',xline(3),',',yline(3),')'
           print*,' dr = ',dr,' dx = ',xlength,' dy = ',ylength
           if (abs(xlength).gt.tiny(xlength)) then
              gradient = ylength/xlength
              yint = yline(3) - gradient*xline(3)
              print*,' gradient = ',gradient,' y intercept = ',yint
              yline(1) = gradient*xline(1) + yint
              yline(4) = gradient*xline(4) + yint
           endif
           !--plot line joining the two points
           call pgline(4,xline,yline)
        case default
           print*,' action cancelled'       
        end select
     !
     !--help
     !
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' left click (or A) : select region'
        print*,'      left click again to zoom on selection'
        print*,'      or press 1-9 to mark selected particles with colour 1-9'
        print*,' left click on colour bar : change rendering limits'
        print*,' +   : zoom in by 10%'
        print*,' -(_): zoom out by 10(20)%'
        print*,' a   : (a)djust/reset plot limits to fit '
        print*,'      NB for these options cursor inside the plot changes both x and y,'
        print*,'       whereas cursor over a specific axis zooms on that axis only'
        if (irender.ne.0) then
           print*,'     (applies to colour bar if mouse is over colour bar)'
        endif
        print*,' l: (l)og / unlog axis (with cursor over the axis to change)'  
        if (irender.ne.0) then
           print*,'  (applies to colour bar if mouse is over colour bar)'
        endif
        print*,' o: re-centre plot on (o)rigin'
        print*,' backspace: delete annotation - over axes       : removes axes'
        print*,'                              - once no axes    : removes legends'
        print*,'                              - once no legend  : removes titles'
        print*,'                              - over colour bar : removes colour bar'
        print*,' r: (r)eplot current plot'
        print*,' R: (R)eset/remove all range restrictions'
        print*,' p: label closest (p)article'
        print*,' t: t)rack closest particle/turn tracking off (coord plots only)'
        print*,' c: plot (c)ircle of interaction for closest particle'
        print*,' C: centre plot on (C)ursor position (keeping aspect ratio)'
        print*,' g: plot a line and find its g)radient'
        print*,' G: move le(G)end to current position'
        print*,' T: move (T)itle to current position'
        print*,' H: move vector plot legend H(ere)'
        print*,' m: change colour m)ap to next (rendered plots only)'
        print*,' M: change colour M)ap to previous (rendered plots only)'
        print*,' i: i)nvert colour map (rendered plots only)'
        print*,' v: decrease arrow size on vector plots (Z for x10)'
        print*,' V: increase arrow size on vector plots (Z for x10)'
        print*,' w: adapt vector arrow size to maximum value'
        if (ndim.ge.3) then
           print*,' k: decrease opacity on opacity-rendered plots (Z for x10)'
           print*,' K: increase opacity on opacity-rendered plots (Z for x10)'
        endif
        print*,' e, E: use current frame/settings as end point to animation sequence'
        if (ndim.gt.1) then
           print*,' , .: rotate about z axis by +(-) 15 degrees (coord plots only)'
           print*,' < >: rotate about z axis by +(-) 30 degrees (coord plots only)'
           if (ndim.ge.3) then
              print*,'[ ]: rotate about x axis by +/- 15 degrees (coord plots only)'
              print*,'{ }: rotate about x axis by +/- 30 degrees (coord plots only)'
              print*,'/ \: rotate about y axis by +/- 15 degrees (coord plots only)'
              print*,'? |: rotate about y axis by +/- 30 degrees (coord plots only)'
              print*,' x: take cross section (coord plots only)'
              if (iplotz.gt.0) then
                 print*,' u: move cross section/perspective position up (towards observer)'
                 print*,' U: move cross section/perspective position up more'
                 print*,' d: move cross section/perspective position down (away from observer)'
                 print*,' D: move cross section/perspective position down more'
              endif
           endif
        endif
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' h: (h)elp'
        print*,' s: (s)ave current settings for all steps'
        print*,' q,Q,esc: (q)uit plotting'
        print*
        print*,' z,Z(oom) : timstepping, zoom and limits-changing options '
        print*,'          are multiplied by a factor of 10'        
        print*,'-------------------------------------------------------'
     case('s','S')
        itrackpart = itrackparttemp
        if (itrackpart.eq.0) then
           call save_limits(iplotx,xmin,xmax)
           call save_limits(iploty,ymin,ymax)
        else
           print*,'tracking particle ',itrackpart,'x,y = ',xcoords(itrackpart),ycoords(itrackpart)
           if (iplotx.le.ndim) call save_limits_track(iplotx,xmin,xmax,xcoords(itrackpart))
           if (iploty.le.ndim) call save_limits_track(iploty,ymin,ymax,ycoords(itrackpart))
           call save_itrackpart_recalcradius()
        endif
        if (irender.gt.0) call save_limits(irender,rendermin,rendermax)
        if (icontour.gt.0) then
           if (icontour.eq.irender &
               .and. abs(rendermin-contmin).le.tiny(0.) &
               .and. abs(rendermax-contmax).le.tiny(0.)) then
              call reset_limits2(icontour)
           elseif (icontour.eq.irender) then
              call save_limits(icontour,contmin,contmax,setlim2=.true.)
           else
              call save_limits(icontour,contmin,contmax)
           endif
        endif
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           call save_limits(ivecx,-vecmax,vecmax)
           call save_limits(ivecy,-vecmax,vecmax)
        endif
        if (ndim.eq.3 .and. iplotz.gt.0 .and. irender.gt.0 .and. use3Dopacity) then
           call save_opacity(taupartdepth)
        endif
        if (ncircpart.gt.0) call save_circles(ncircpart,icircpart)
        if (rotation) call save_rotation(ndim,anglex,angley,anglez)
        if (iplotz.gt.0) then
           if (x_sec) then
              call save_xsecpos(zslicepos,x_sec)
           else
              call save_perspective(zobserver,dscreen)
           endif
        endif
        print*,'> plot settings saved (...but NOT to disk, use ''S'' from the main menu for that) <'
     !
     !--actions on left click
     !
     case('A') ! left click
        print*,'select area: '
        print*,'left click : zoom'
        !
        !--change colour bar limits
        !
        ishape = inshape(xpt,ypt,itrans(iplotx),itrans(iploty))
        if (ishape.gt.0) then
           call edit_shape(ishape,xpt,ypt,itrans(iplotx),itrans(iploty))
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        elseif (iamincolourbar .and. irender.gt.0) then
           if (incolourbarlabel(iColourBarStyle,xpt,ypt,xmin,xmax,ymin,ymax)) then
              if (verticalbar) then
                 call edit_textbox(xpt,ypt,90.,labelrender)
                 projlabelformat = trim(labelrender)
                 iapplyprojformat = irender
              else
                 call edit_textbox(xpt,ypt,0.,labelrender)
                 projlabelformat = trim(labelrender)
                 iapplyprojformat = irender
              endif
              iadvance = 0
              interactivereplot = .true.
              iexit = .true.
           else
              print*,'click to set rendering limits'
              if (verticalbar) then
                 call pgband(3,1,xpt,ypt,xpt2,ypt2,char2)
              else
                 call pgband(4,1,xpt,ypt,xpt2,ypt2,char2)           
              endif
              if (char2 == 'A') then
                 if (verticalbar) then
                    drender = (rendermax-rendermin)/(ymax-ymin)
                    rendermax = rendermin + (max(ypt,ypt2)-ymin)*drender
                    rendermin = rendermin + (min(ypt,ypt2)-ymin)*drender
                 else
                    drender = (rendermax-rendermin)/(xmax-xmin)
                    rendermax = rendermin + (max(xpt,xpt2)-xmin)*drender
                    rendermin = rendermin + (min(xpt,xpt2)-xmin)*drender              
                 endif
                 print*,'setting render min = ',rendermin
                 print*,'setting render max = ',rendermax
                 iadvance = 0
                 interactivereplot = .true.
                 iexit = .true.
              endif
           endif
        !
        !--zoom or mark particles
        !
        else
           if (irender.le.0) then
              print*,'1-9 = mark selected particles with colours 1-9'
           endif
           print*,'0 = hide selected particles'
           print*,'p = plot selected particles only'
           print*,'c = plot circles of interaction on selected parts'
           print*,'x = use particles within x parameter range only'
           print*,'y = use particles within y parameter range only'
           print*,'r = use particles within x and y parameter range only'
           print*,'R = remove all range restrictions'
           call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
           xptmin = min(xpt,xpt2)
           xptmax = max(xpt,xpt2)
           yptmin = min(ypt,ypt2)
           yptmax = max(ypt,ypt2)
           print*,'xrange = ',xptmin,'->',xptmax
           print*,'yrange = ',yptmin,'->',yptmax
           if (iplotz.ne.0 .and. x_sec) then
              print*,'(zrange = ',zptmin,'->',zptmax,')'
           endif
           select case (char2)
           case('A')   ! zoom if another left click
              call pgsfs(2)
              call pgrect(xpt,xpt2,ypt,ypt2)
              xmin = xptmin
              xmax = xptmax
              ymin = yptmin
              ymax = yptmax
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           case('0','1','2','3','4','5','6','7','8','9') ! mark particles
              if (irender.le.0 .or. char2.eq.'0' .or. char2.eq.'1') then
                 nmarked = 0
                 do i=1,npart
                    if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                    .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax) &
                    .and.(zcoords(i).ge.zptmin .and. zcoords(i).le.zptmax)) then
                        read(char2,*,iostat=ierr) icolourpart(i)
                        if (ierr /=0) then
                           print*,'*** error marking particle' 
                           icolourpart(i) = 1
                        !--translate 0 to icolourpart = -1 for non-plotted particles
                        elseif (icolourpart(i).eq.0) then
                           icolourpart(i) = -1
                        endif
                        nmarked = nmarked + 1
                    endif
                 enddo
                 print*,'marked ',nmarked,' particles in selected region'
              endif
              iadvance = 0
              if (nmarked.gt.0) irerender = .true.
              interactivereplot = .true.
              iexit = .true.
           case('p') ! plot selected particles only
              nmarked = 0
              do i=1,npart
                 if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                 .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax) &
                 .and.(zcoords(i).ge.zptmin .and. zcoords(i).le.zptmax)) then
                    nmarked = nmarked + 1
                    if (icolourpart(i).le.0) icolourpart(i) = 1
                 else
                    icolourpart(i) = -1
                 endif
              enddo
              print*,'plotting selected ',nmarked,' particles only'
              iadvance = 0
              irerender = .true.
              interactivereplot = .true.
              iexit = .true.           
           case('c') ! set circles of interaction in marked region
              ncircpart = 0
              do i=1,npart
                 if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                 .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax) &
                 .and.(zcoords(i).ge.zptmin .and. zcoords(i).le.zptmax)) then
                     if (ncircpart.lt.size(icircpart)) then
                        ncircpart = ncircpart + 1                        
                        icircpart(ncircpart) = i
                        call pgsfs(2)
                        call pgcirc(xcoords(i),ycoords(i),2.*hi(i))
                     endif
                 endif
              enddo
              print*,'set ',ncircpart,' circles of interaction in selected region'
              if (ncircpart.eq.size(icircpart)) print*,' (first ',size(icircpart),' only)'
           case('x')
              call restrict_range(iplotx,xptmin,xptmax)
              iadvance = 0
              irerender = .true.
              interactivereplot = .true.
              iexit = .true.
           case('y')
              call restrict_range(iploty,yptmin,yptmax)
              iadvance = 0
              irerender = .true.
              interactivereplot = .true.
              iexit = .true.
           case('r')
              call restrict_range(iplotx,xptmin,xptmax)
              call restrict_range(iploty,yptmin,yptmax)
              iadvance = 0
              irerender = .true.
              interactivereplot = .true.
              iexit = .true.
           case('R')
              call reset_ranges
              iadvance = 0
              irerender = .true.
              interactivereplot = .true.
              iexit = .true.
           case default
              print*,' action cancelled'
           end select 
        endif   
     !
     !--zooming
     !
     case('-','_','+','o','C') ! zoom in/out
        xlength = xmax - xmin
        ylength = ymax - ymin
        xcen = 0.5*(xmax + xmin)
        ycen = 0.5*(ymax + ymin)
        renderlength = rendermax - rendermin
        select case(char)
        case('-')
           xlength = 1.1*zoomfac*xlength
           ylength = 1.1*zoomfac*ylength
           renderlength = 1.1*zoomfac*renderlength
        case('_')
           xlength = 1.2*zoomfac*xlength
           ylength = 1.2*zoomfac*ylength
           renderlength = 1.2*zoomfac*renderlength
        case('+')
           xlength = 0.9/zoomfac*xlength
           ylength = 0.9/zoomfac*ylength
           renderlength = 0.9/zoomfac*renderlength
        case('o')
           if (itrackpart.gt.0) then
               print*,'centreing limits on tracked particle ',itrackpart,'x,y = ',xcoords(itrackpart),ycoords(itrackpart)
               xcen = xcoords(itrackpart)
               ycen = ycoords(itrackpart)
           else
              if (iplotx.le.ndim) then
                 xcen = xorigin(iplotx)
              else
                 xcen = 0.
              endif
              if (iploty.le.ndim) then
                 ycen = xorigin(iploty)
              else
                 ycen = 0.
              endif
              print*,' centreing plot on origin x,y = ',xcen,ycen
           endif
        case('C')
           xcen = xpt
           ycen = ypt
        end select
        if (iamincolourbar .and. irender.gt.0) then
           !--rendering zoom does not allow pan - renderpt is always centre of axis
           renderpt = 0.5*(rendermin + rendermax)
           rendermin = renderpt - 0.5*renderlength
           rendermax = renderpt + 0.5*renderlength
           print*,'zooming on colour bar: min, max = ',rendermin,rendermax
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        else
           if (xpt.ge.xmin .and. xpt.le.xmax .and. ypt.le.ymaxin) then
              xmin = xcen - 0.5*xlength
              xmax = xcen + 0.5*xlength
              print*,'zooming on x axis: min, max = ',xmin,xmax
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           endif
           if (ypt.ge.ymin .and. ypt.le.ymax .and. xpt.le.xmaxin) then
              ymin = ycen - 0.5*ylength
              ymax = ycen + 0.5*ylength
              print*,'zooming on y axis: min, max = ',ymin,ymax
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           endif
        endif
        if (char.eq.'o') then
           xpt = xcen
           ypt = ycen
        endif
     case('a') ! reset plot limits
        if (iamincolourbar .and. irender.gt.0) then
           rendermin = renderminadapt
           rendermax = rendermaxadapt
           iadvance = 0              ! that it should change the render limits
           interactivereplot = .true.
           iexit = .true.
        else
           xmaxin = xmax
           if (xpt.ge.xmin .and. xpt.le.xmax .and. ypt.le.ymaxin) then
              print*,'resetting x limits'
              xmin = minval(xcoords,mask=(icolourpart.ge.0))
              xmax = maxval(xcoords,mask=(icolourpart.ge.0))
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           endif
           if (ypt.ge.ymin .and. ypt.le.ymax .and. xpt.le.xmaxin) then
              print*,'resetting y limits'
              ymin = minval(ycoords,mask=(icolourpart.ge.0))
              ymax = maxval(ycoords,mask=(icolourpart.ge.0))
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           endif
        endif
     !
     !--zoom in/out on vector plots (arrow size)
     !
     case('v')
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'decreasing vector arrow size'
           vecmax = 1.2*zoomfac*vecmax
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        endif
     case('V')
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'increasing vector arrow size'
           vecmax = 0.8/zoomfac*vecmax
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        endif
     case('w','W')
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'adapting vector arrow size'
           vecmax = -1.0
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        endif
     !
     !--change opacity on 3D opacity rendered plots
     !
     case('k')
        if (ndim.eq.3 .and. iplotz.gt.0 .and. irender.ne.0 .and. use3Dopacity) then
           print*,'decreasing opacity by factor of ',1.5*zoomfac
           !--opacity goes like 1/taupartdepth
           taupartdepth = 1.5*zoomfac*taupartdepth
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('K')
        if (ndim.eq.3 .and. iplotz.gt.0 .and. irender.ne.0 .and. use3Dopacity) then
           print*,'increasing opacity by factor of ',1.5*zoomfac
           !--opacity goes like 1/taupartdepth
           taupartdepth = taupartdepth/(1.5*zoomfac)
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     !
     !--set/unset log axes
     !
     case('l')
        !
        !--change colour bar, y and x itrans between log / not logged
        !
        if (iamincolourbar .and. irender.gt.0) then
           if (icontour.eq.irender) then
              call change_itrans2(irender,rendermin,rendermax,contmin,contmax)
           else
              call change_itrans(irender,rendermin,rendermax)
           endif
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        elseif (xpt.lt.xmin) then
           if (iploty.le.ndim .and. irender.gt.0) then
              print "(a)",'error: cannot log coordinate axes with rendering'
           else
              call change_itrans(iploty,ymin,ymax)
              iadvance = 0
              interactivereplot = .true.
              iexit = .true.
           endif
        elseif (ypt.lt.ymin) then
           if (iplotx.le.ndim .and. irender.gt.0) then
              print "(a)",'error: cannot log coordinate axes with rendering'
           else
              call change_itrans(iplotx,xmin,xmax)
              iadvance = 0
              interactivereplot = .true.
              iexit = .true.
           endif
        endif
     !
     !--reset all range restrictions
     !
     case('R')
        call reset_ranges
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     !
     !--save as end point of animation sequence
     !
     case('e','E')
        call setsequenceend(istep,iplotx,iploty,irender,rotation, &
                          anglex,angley,anglez,zobserver,use3Dopacity,taupartdepth, &
                          x_sec,zslicepos,xmin,xmax,ymin,ymax,rendermin,rendermax)
     !
     !--rotation
     !
     case(',')
        if (rotation) then
           print*,'changing z rotation angle by -15 degrees...'
           anglez = anglez - 15.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('<')
        if (rotation) then
           print*,'changing z rotation angle by -30 degrees...'
           anglez = anglez - 30.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('.')
        if (rotation) then
           print*,'changing z rotation angle by 15 degrees...'
           anglez = anglez + 15.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('>')
        if (rotation) then
           print*,'changing z rotation angle by 30 degrees...'
           anglez = anglez + 30.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('/')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by -15 degrees...'
           angley = angley - 15.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('?')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by -30 degrees...'
           angley = angley - 30.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('\')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by 15 degrees...'
           angley = angley + 15.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('|')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by 30 degrees...'
           angley = angley + 30.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('[')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by -15 degrees...'
           anglex = anglex - 15.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('{')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by -30 degrees...'
           anglex = anglex - 30.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case(']')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by 15 degrees...'
           anglex = anglex + 15.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     case('}')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by 30 degrees...'
           anglex = anglex + 30.
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     !
     !--set cross section position
     !
     case('x')
        if (rotation .and. ndim.ge.3) then
           xline(1) = xpt
           yline(1) = ypt
           !--work out which is the third dimension
           do i=1,3
              if (i.ne.iplotx .and. i.ne.iploty) ixsec = i
           enddo
           print*,' select cross section position (using left click or x)'
           call pgband(1,1,xline(1),yline(1),xline(2),yline(2),char2)
           !--work out cross section if left click or x again
           select case(char2)
           case('A','x')
              !--plot the cross section line
              call pgline(2,xline(1:2),yline(1:2))
              !--work out angle with the x axis
              !  and offset of line from origin
              dx = xline(2) - xline(1)
              dy = yline(2) - yline(1)
              anglerad = ATAN2(dy,dx)
              select case(ixsec)
              case(1)
                 anglex = 180.*anglerad/pi + anglex
                 print*,'setting angle x = ',anglex
              case(2)
                 angley = 180.*anglerad/pi + angley
                 print*,'setting angle y = ',angley
              case(3)
                 anglez = 180.*anglerad/pi + anglez
                 print*,'setting angle z = ',anglez
              end select
              iploty = ixsec
              !--work out offset of cross section line
              ! y intercept
              yint = yline(2) - (dy/dx)*xline(2)
              zslicepos = yint/COS(anglerad)
              !--if we are in column density mode, change back to cross-section mode
              x_sec = .true.
              print*,'iploty = ',ixsec, ' xsecpos = ',zslicepos
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           case default
              print*,' action cancelled'
           end select
        endif
     !
     !--cross sections
     !
     case('u') ! move cross section up by dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position up by ',dzslice
              zslicepos = zslicepos + dzslice
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           else
              print*,'shifting perspective position up ',dscreen
              zobserver = zobserver + dscreen
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.           
           endif
        endif
     case('U') ! move cross section up by 2*dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position up by ',2.*dzslice
              zslicepos = zslicepos + 2.*dzslice
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.
           else
              print*,'shifting perspective position up by ',2.*dscreen
              zobserver = zobserver + 2.*dscreen
              iadvance = 0
              interactivereplot = .true.
              irerender = .true.
              iexit = .true.           
           endif
        endif
     case('d') ! move cross section down by dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position down by ',dzslice
              zslicepos = zslicepos - dzslice
           else
              print*,'shifting perspective position down by ',dscreen
              zobserver = zobserver - dscreen
           endif
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.           
        endif     
     case('D') ! move cross section down by 2*dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position down by ',2.*dzslice
              zslicepos = zslicepos - 2.*dzslice
           else
              print*,'shifting perspective position down by ',2.*dscreen
              zobserver = zobserver - 2.*dscreen
           endif
           iadvance = 0
           interactivereplot = .true.
           irerender = .true.
           iexit = .true.
        endif
     !
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        call mvlegend(xpt,ypt,xmin,xmax,ymax)
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case('T') ! move title here
        print*,'setting title position to current location...'
        call mvtitle(xpt,ypt,xmin,xmax,ymax)
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case('H') ! move vector legend here
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'setting vector plot legend to current location...'
           call mvlegendvec(xpt,ypt,xmin,xmax,ymax)
        endif
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case('m') ! change colour map (next scheme)
        call change_colourmap(icolourscheme,1)
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case('M') ! change colour map (previous scheme)
        call change_colourmap(icolourscheme,-1)
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case('i') ! invert colour map
        icolourscheme = -icolourscheme
        call change_colourmap(icolourscheme,0)
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case(achar(8)) ! delete plot annotation / colour bar (backspace)
        ishape = inshape(xpt,ypt,itrans(iplotx),itrans(iploty))
        if (ishape.gt.0) then
           call delete_shape(ishape)
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        elseif (iamincolourbar .and. irender.gt.0) then
           iColourBarStyle = 0
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        elseif (xpt.lt.xmin .or. xpt.gt.xmax .or. ypt.lt.ymin .or. ypt.gt.ymax) then
           call deleteaxes()
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        endif
     !
     !--timestepping
     !
     case('q','Q',achar(27),achar(3))
        iadvance = -666
        print*,'quitting...'
        iexit = .true.
     case('X','b','B') ! right click -> go back
        iadvance = -abs(iadvance)
        iexit = .true.
     case('r') ! replot
        iadvance = 0
        interactivereplot = .true.
        irerender = .true.
        iexit = .true.
     case(' ','n','N') ! space
        iadvance = abs(iadvance)
        iexit = .true.
     case('0','1','2','3','4','5','6','7','8','9')
        read(char,*,iostat=ierr) iadvance
        if (ierr /=0) then
           print*,'*** internal error setting timestep jump' 
           iadvance = 1
        endif
        iadvance = int(zoomfac*iadvance)
        print*,' setting timestep jump = ',iadvance
     case(')')
        iadvance = int(zoomfac*10)
        print*,' setting timestep jump = ',iadvance
     !
     !--multiply everything by a factor of 10     
     !
     case('z','Z')
        zoomfac = 10.*zoomfac
        if (zoomfac.gt.1000000.) then
           zoomfac = 1.0
        endif
        print*,' LIMITS/TIMESTEPPING CHANGES NOW x ',zoomfac
     !
     !--unknown
     !
     case default
        print*,' x, y = ',xpt,ypt,'; unknown option "',trim(char), '" ',iachar(char)
     end select
     
     !
     !--save cursor position relative to the viewport
     !
     call pgqwin(xminwin,xmaxwin,yminwin,ymaxwin)
     call get_vptxy(xpt,ypt,xcursor,ycursor,xminwin,xmaxwin,yminwin,ymaxwin)

     if (rotation) then
        if (anglez.ge.360.) anglez = anglez - 360.
        if (anglez.lt.0.) anglez = anglez + 360.
        if (ndim.gt.2) then
           if (angley.ge.360.) angley = angley - 360.
           if (angley.lt.0.) angley = angley + 360.
           if (anglex.ge.360.) anglex = anglex - 360.
           if (anglex.lt.0.) anglex = anglex + 360.        
        endif
     endif
     !
     !--do not let timestep go outside of bounds
     !  if we are at the first/last step, just print message and do nothing
     !  if iadvance trips over the bounds, jump to last/first step
     !
     if (iadvance.ne.-666 .and. iexit) then
        if (istep + iadvance .gt. ilaststep .and. iframe.eq.nframes) then
           print "(1x,a)",'reached last timestep'
           if (ilaststep-istep .gt.0) then
              iadvance= ilaststep - istep
           else
              iexit = .false.
           endif
        elseif (istep + iadvance .lt. 1 .and. iframe.eq.1) then
           print "(1x,a)",'reached first timestep: can''t go back'
           if (1-istep .lt.0) then
              iadvance= 1 - istep
           else
              iexit = .false.
           endif
        endif
     endif

  enddo interactiveloop
  return
end subroutine interactive_part

!
! cut down version of interactive mode -> controls timestepping only
! used in powerspectrum / extra plots
! THIS IS NOW LARGELY OBSOLETE (superseded by interactive_multi)
!  AND WILL BE REMOVED IN FUTURE VERSIONS
!
subroutine interactive_step(iadvance,istep,ilaststep,xmin,xmax,ymin,ymax,interactivereplot)
 implicit none
 integer, intent(inout) :: iadvance
 integer, intent(in) :: istep,ilaststep
 real, intent(inout) :: xmin,xmax,ymin,ymax
 logical, intent(out) :: interactivereplot
 integer :: nc,ierr
 real :: xpt,ypt,xpt2,ypt2
 real :: xlength, ylength, zoomfac
 character(len=1) :: char,char2
 character(len=5) :: string
 logical :: iexit
 
  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  xpt = 0.
  ypt = 0.
  zoomfac = 1.0
  iexit = .false.
  interactivereplot = .false.
  
  do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
     print*,'x, y = ',xpt,ypt,' function = ',char
     
     select case(char)
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' select area and zoom : left click (or A)'
        print*,' zoom in by 10%       : +'
        print*,' zoom out by 10(20)%      : - (_)'
        print*,' (r)eplot current plot        : r'
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' G : move legend to current position'
        print*,' T : move title to current position'
        print*,' (h)elp                       : h'
        print*,' (q)uit plotting              : q, Q, esc'             
        print*,'-------------------------------------------------------'

     case('A') ! left click
        !
        !--draw rectangle from the point and reset the limits
        !
        print*,'select area: '
        print*,'left click : zoom'
        call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
        print*,xpt,ypt,xpt2,ypt2,char2
        select case (char2)
        case('A')   ! zoom if another left click
           call pgrect(xpt,xpt2,ypt,ypt2)
           xmin = min(xpt,xpt2)
           xmax = max(xpt,xpt2)
           ymin = min(ypt,ypt2)
           ymax = max(ypt,ypt2)
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        case default
           print*,' action cancelled'
        end select    
     !
     !--zooming
     !
     case('-','_','+','o') ! zoom out by 10 or 20%
        xlength = xmax - xmin
        ylength = ymax - ymin
        select case(char)
        case('-')
           xlength = 1.1*zoomfac*xlength
           ylength = 1.1*zoomfac*ylength
        case('_')
           xlength = 1.2*zoomfac*xlength
           ylength = 1.2*zoomfac*ylength
        case('+')
           xlength = 0.9/zoomfac*xlength
           ylength = 0.9/zoomfac*ylength
        case('o') !--reset cursor to origin
           xpt = 0.
           ypt = 0.
        end select
        if (xpt.ge.xmin .and. xpt.le.xmax .and. ypt.le.ymax) then
           print*,'zooming on x axis'
           xmin = xpt - 0.5*xlength
           xmax = xpt + 0.5*xlength
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        endif
        if (ypt.ge.ymin .and. ypt.le.ymax .and. xpt.le.xmax) then
           print*,'zooming on y axis'
           ymin = ypt - 0.5*ylength
           ymax = ypt + 0.5*ylength
           iadvance = 0
           interactivereplot = .true.
           iexit = .true.
        endif
     !
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        call mvlegend(xpt,ypt,xmin,xmax,ymax)
     case('T') ! move title here
        print*,'setting title position to current location...'
        call mvtitle(xpt,ypt,xmin,xmax,ymax)
     !
     !--timestepping
     !
     case('q','Q',achar(27),achar(3))
        iadvance = -666
        print*,'quitting...'
        iexit = .true.
     case('X','b','B') ! right click -> go back
        iadvance = -abs(iadvance)
        iexit = .true.
     case('r','R') ! replot
        iadvance = 0
        interactivereplot = .true.
        iexit = .true.
     case(' ','n','N') ! space
        iadvance = abs(iadvance)
        iexit = .true.
     case('0','1','2','3','4','5','6','7','8','9')
        read(char,*,iostat=ierr) iadvance
        if (ierr /=0) then
           print*,'*** internal error setting timestep jump' 
           iadvance = 1
        endif
        iadvance = int(zoomfac*iadvance)
        print*,' setting timestep jump = ',iadvance
     case(')')
        iadvance = int(zoomfac*10)
        print*,' setting timestep jump = ',iadvance
     !
     !--multiply everything by a factor of 10     
     !
     case('z','Z')
        zoomfac = 10.*zoomfac
        if (zoomfac.gt.1000000.) then
           zoomfac = 1.0
        endif
        print*,' LIMITS/TIMESTEPPING CHANGES NOW x ',zoomfac
     !
     !--unknown
     !
     case default
        print*,' x, y = ',xpt,ypt,'; unknown option "',trim(char), '" ',iachar(char)
     end select
     !
     !--do not let timestep go outside of bounds
     !  if we are at the first/last step, just print message and do nothing
     !  if iadvance trips over the bounds, jump to last/first step
     !
     if (iadvance.ne.-666 .and. iexit) then
        if (istep + iadvance .gt. ilaststep) then
           print "(1x,a)",'reached last timestep'
           if (ilaststep-istep .gt.0) then
              iadvance= ilaststep - istep
           else
              iexit = .false.
           endif
        elseif (istep + iadvance .lt. 1) then
           print "(1x,a)",'reached first timestep: can''t go back'
           if (1-istep .lt.0) then
              iadvance= 1 - istep
           else
              iexit = .false.
           endif
        endif
     endif
  enddo
  return
end subroutine interactive_step

!
! interactive mode for multiple plots per page - requires determination of which plot/panel
!  a mouse-click refers to from stored settings for the viewport and limits for each plot.
! (this could be made into the only subroutine required)
!
subroutine interactive_multi(iadvance,istep,ifirststeponpage,ilaststep,iframe,ifirstframeonpage,nframes, &
                             lastpanel,iplotxarr,iplotyarr,irenderarr,ivecarr,xmin,xmax,vptxmin,vptxmax,vptymin,vptymax, &
                             barwmulti,xminadapt,xmaxadapt,nacross,ndim,xorigin,icolourscheme, &
                             iColourBarStyle,interactivereplot)
 implicit none
 integer, intent(inout) :: iadvance
 integer, intent(inout) :: istep,iframe,lastpanel,iColourBarStyle
 integer, intent(in) :: ifirststeponpage,ilaststep,nacross,ndim,ifirstframeonpage,nframes
 integer, intent(inout) :: icolourscheme
 integer, intent(in), dimension(:) :: iplotxarr,iplotyarr,irenderarr,ivecarr
 real, dimension(:), intent(in) :: vptxmin,vptxmax,vptymin,vptymax,barwmulti
 real, dimension(:), intent(inout) :: xmin,xmax,xminadapt,xmaxadapt
 real, intent(in), dimension(ndim) :: xorigin
 logical, intent(out) :: interactivereplot
 integer :: nc,ierr,ipanel,ipanel2,istepin,istepnew,i,istepjump,istepsonpage
 real :: xpt,ypt,xpt2,ypt2,xpti,ypti,renderpt,xptmin,xptmax,yptmin,yptmax
 real :: xlength,ylength,renderlength,drender,zoomfac
 real :: vptxi,vptyi,vptx2i,vpty2i,vptxceni,vptyceni
 real :: xmini,xmaxi,ymini,ymaxi,xcen,ycen,gradient,dr,yint,xmaxin
 real, dimension(4) :: xline,yline
 character(len=1) :: char,char2
 character(len=5) :: string
 logical :: iexit,iamincolourbar,verticalbar
 
  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  zoomfac = 1.0
  xpt2 = 0.
  ypt2 = 0.
  !
  !--convert saved cursor position (saved in viewport coords)
  !  back to world coordinates:
  !
  !--query window limits in world coords
  call pgqwin(xmini,xmaxi,ymini,ymaxi)

  !--determine which plot the cursor falls on
  ipanel = getpanel(xcursor,ycursor)

  !--set the window to correspond to the panel we last left the cursor in
  call set_panel(ipanel)

  !--set the position in x,y relative to this panel
  call getxy(xcursor,ycursor,xpt,ypt,ipanel)

  iexit = .false.
  interactivereplot = .false.
  istepin = istep
  istepnew = ifirststeponpage - iadvance
  istepsonpage = abs(istep - ifirststeponpage)/iadvance + 1
  istepjump = 1
!  print*,'istep = ',istepnew
!  print*,'steps on page = ',istepsonpage
  verticalbar = barisvertical(iColourBarStyle)
  
  interactive_loop: do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
!     print*,'x, y = ',xpt,ypt,' function = ',char
     !
     !--determine which plot the cursor falls on
     !
     call get_vptxy(xpt,ypt,vptxi,vptyi,xmini,xmaxi,ymini,ymaxi)
     ipanel = getpanel(vptxi,vptyi)
     !print*,'xpt,ypt = ',xpt,ypt,vptxi,vptyi,ipanel

     !--translate vpt co-ords to x,y in current panel
     call getxy(vptxi,vptyi,xpti,ypti,ipanel)

     !--query the position of the colour bar
     if (ipanel.gt.0) then
        iamincolourbar = incolourbar(iColourBarStyle,xpti,ypti,xmin(iplotxarr(ipanel)), &
                         xmax(iplotxarr(ipanel)),xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
     else
        iamincolourbar = .false.
     endif

     select case(char)
     case('h')
        print*,'----- interactive mode commands (multiple plots per page) -----'
        print*,' left click (or A) : select region'
        print*,'  left click again to zoom on selection'
        print*,'  (select colour bar to change rendering limits)'
        print*,' +   : zoom in by 10%'
        print*,' -(_): zoom out by 10(20)%'
        print*,' a   : (a)djust/reset plot limits to fit '
        print*,'      NB for these options cursor inside the plot changes both x and y,'
        print*,'       whereas cursor over a specific axis zooms on that axis only'
        print*,'       (applies to colour bar if mouse is over colour bar)'
        print*,' l: (l)og / unlog axis (with cursor over the axis to change)'  
        print*,'    (applies to colour bar if mouse is over colour bar)'
        print*,' o: re-centre plot on (o)rigin'
        print*,' C: centre plot on (C)ursor position (keeping aspect ratio)'        
        print*,' backspace: delete annotation - over axes       : removes axes'
        print*,'                              - once no axes    : removes legends'
        print*,'                              - once no legend  : removes titles'
        print*,'                              - over colour bar : removes colour bar'
        print*,' r: (r)eplot current plot'
        print*,' R: (R)eset/remove all range restrictions'
        print*,' g: plot a line and find its g)radient'
        print*,' G: move le(G)end to current position'
        print*,' T: move (T)itle to current position'
        print*,' H: move vector plot legend H(ere)'
        print*,' m: change colour m)ap to next'
        print*,' M: change colour M)ap to previous'
        print*,' i: i)nvert colour map'
        print*,' v: decrease arrow size on vector plots'
        print*,' V: increase arrow size on vector plots'
        print*,' w: adapt vector arrow size to maximum value'
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' h: (h)elp'
        print*,' s: (s)ave current settings for all steps'
        print*,' q,Q: (q)uit plotting'
        print*
        print*,' z,Z(oom) : timstepping, zoom and limits-changing options '
        print*,'          are multiplied by a factor of 10'        
        print*,'---------------------------------------------------------------'
     case('g')   ! draw a line between two points
        xline(2) = xpti
        yline(2) = ypti
        call set_panel(ipanel)
        !--mark first point
        call pgpt1(xpti,ypti,4)
        !--select second point
        print*,' select another point (using left click or g) to plot line '
        call pgband(1,1,xline(2),yline(2),xline(3),yline(3),char2)
        !--draw line if left click or g
        select case(char2)
        case('A','g')
           print*
           !--mark second point
           call pgpt1(xline(3),yline(3),4)
           xlength = xline(3)-xline(2)
           if (abs(xlength).lt.tiny(xlength)) then
              xline(1) = xline(2)
              xline(4) = xline(2)
              yline(1) = xmin(iplotyarr(ipanel))
              yline(4) = xmax(iplotyarr(ipanel))
              print*,' error: gradient = infinite'
           elseif (xline(2).lt.xline(3)) then 
              xline(1) = xmin(iplotxarr(ipanel))
              xline(4) = xmax(iplotxarr(ipanel))
           else
              xline(1) = xmax(iplotxarr(ipanel))
              xline(4) = xmin(iplotxarr(ipanel))
           endif           
           ylength = yline(3)-yline(2)
           dr = sqrt(xlength**2 + ylength**2)
           print*,' (x1,y1) = (',xline(2),',',yline(2),')'
           print*,' (x2,y2) = (',xline(3),',',yline(3),')'
           print*,' dr = ',dr,' dx = ',xlength,' dy = ',ylength
           if (abs(xlength).gt.tiny(xlength)) then
              gradient = ylength/xlength
              yint = yline(3) - gradient*xline(3)
              print*,' gradient = ',gradient,' y intercept = ',yint
              yline(1) = gradient*xline(1) + yint
              yline(4) = gradient*xline(4) + yint
           endif
           !--plot line joining the two points
           call pgline(4,xline,yline)
           call reset_panel
        case default
           print*,' action cancelled'       
        end select
     case('s','S')
        do i=1,size(vptxmin)
           call save_limits(iplotxarr(i),xmin(iplotxarr(i)),xmax(iplotxarr(i)))
           call save_limits(iplotyarr(i),xmin(iplotyarr(i)),xmax(iplotyarr(i)))
           if (irenderarr(i).gt.0) call save_limits(irenderarr(i),xmin(irenderarr(i)),xmax(irenderarr(i)))
        enddo
        print*,'> plot settings saved <'
     case('A') ! left click
        !
        !--draw rectangle from the point and reset the limits
        !
        print*,'select area: '
        print*,'left click : zoom'
        print*,'x = use particles within x parameter range only'
        print*,'y = use particles within y parameter range only'
        print*,'r = use particles within x and y parameter range only'
        print*,'R = remove all range restrictions'

        !
        !--change colour bar limits
        !
        if (ipanel.gt.0 .and. iamincolourbar .and. irenderarr(ipanel).gt.0) then
           print*,'click to set rendering limits'
           if (verticalbar) then
              call pgband(3,1,xpt,ypt,xpt2,ypt2,char2)
           else
              call pgband(4,1,xpt,ypt,xpt2,ypt2,char2)           
           endif
           if (char2 == 'A') then
              call get_vptxy(xpt2,ypt2,vptx2i,vpty2i,xmini,xmaxi,ymini,ymaxi)
              if (barwmulti(ipanel).gt.tiny(barwmulti)) then
                 if (verticalbar) then
                    drender = (xmax(irenderarr(ipanel))-xmin(irenderarr(ipanel)))/ &
                              (vptymax(ipanel) -vptymin(ipanel))
                    xmax(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (max(vptyi,vpty2i)-vptymin(ipanel))*drender
                    xmin(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (min(vptyi,vpty2i)-vptymin(ipanel))*drender               
                 else
                    drender = (xmax(irenderarr(ipanel))-xmin(irenderarr(ipanel)))/ &
                              (vptxmax(ipanel) -vptxmin(ipanel))
                    xmax(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (max(vptxi,vptx2i)-vptxmin(ipanel))*drender
                    xmin(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (min(vptxi,vptx2i)-vptxmin(ipanel))*drender
                 endif
              else
              !--for global colour bars (ie. on tiled plots) use viewport co-ordinates to set render limits
                 if (verticalbar) then
                    drender = (xmax(irenderarr(ipanel))-xmin(irenderarr(ipanel)))/ &
                              (maxval(vptymax) - minval(vptymin))
                    xmax(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (max(vptyi,vpty2i)-minval(vptymin))*drender
                    xmin(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (min(vptyi,vpty2i)-minval(vptymin))*drender
                 else
                    drender = (xmax(irenderarr(ipanel))-xmin(irenderarr(ipanel)))/ &
                              (maxval(vptxmax) - minval(vptxmin))
                    xmax(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (max(vptxi,vptx2i)-minval(vptxmin))*drender
                    xmin(irenderarr(ipanel)) = xmin(irenderarr(ipanel)) + (min(vptxi,vptx2i)-minval(vptxmin))*drender
                 endif             
              endif
              print*,'setting render min, max = ',xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
        else
           call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
           !call pgrect(xpt,xpt2,ypt,ypt2)
           call get_vptxy(xpt2,ypt2,vptx2i,vpty2i,xmini,xmaxi,ymini,ymaxi)
           !--use centre point of first click and current click to
           !  better determine panel
           vptxceni = 0.5*(vptxi + vptx2i)
           vptyceni = 0.5*(vptyi + vpty2i)
           ipanel2 = getpanel(vptxceni,vptyceni)
           if (ipanel2.gt.0 .and. ipanel2.ne.ipanel) then
              ipanel = ipanel2
              print*,'panel = ',ipanel
           endif
           if (ipanel.le.0) cycle interactive_loop
           call getxy(vptx2i,vpty2i,xpt2,ypt2,ipanel)
           !--reset first point according to current panel
           call getxy(vptxi,vptyi,xpti,ypti,ipanel)
           xptmin = min(xpti,xpt2)
           xptmax = max(xpti,xpt2)
           yptmin = min(ypti,ypt2)
           yptmax = max(ypti,ypt2)
              
           select case (char2)
           case('A')   ! zoom if another left click
              xmin(iplotxarr(ipanel)) = xptmin
              xmax(iplotxarr(ipanel)) = xptmax
              xmin(iplotyarr(ipanel)) = yptmin
              xmax(iplotyarr(ipanel)) = yptmax
              print*,'setting limits: xmin = ',xmin(iplotxarr(ipanel)),' xmax = ',xmax(iplotxarr(ipanel))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           case('x')
              call restrict_range(iplotxarr(ipanel),xptmin,xptmax)
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           case('y')
              call restrict_range(iplotyarr(ipanel),yptmin,yptmax)
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           case('r')
              call restrict_range(iplotxarr(ipanel),xptmin,xptmax)
              call restrict_range(iplotyarr(ipanel),yptmin,yptmax)
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           case('R')
              call reset_ranges
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           case default
              print*,' action cancelled'
           end select
        endif
     !
     !--zooming
     !
     case('-','_','+','o','C') ! zoom in/out
        if (ipanel.le.0) cycle interactive_loop
        xlength = xmax(iplotxarr(ipanel)) - xmin(iplotxarr(ipanel))
        ylength = xmax(iplotyarr(ipanel)) - xmin(iplotyarr(ipanel))
        xcen = 0.5*(xmax(iplotxarr(ipanel)) + xmin(iplotxarr(ipanel)))
        ycen = 0.5*(xmax(iplotyarr(ipanel)) + xmin(iplotyarr(ipanel)))
        if (irenderarr(ipanel).gt.0) then
           renderlength = xmax(irenderarr(ipanel)) - xmin(irenderarr(ipanel))
        else
           renderlength = 0.
        endif
        select case(char)
        case('-')
           xlength = 1.1*zoomfac*xlength
           ylength = 1.1*zoomfac*ylength
           renderlength = 1.1*zoomfac*renderlength
        case('_')
           xlength = 1.2*zoomfac*xlength
           ylength = 1.2*zoomfac*ylength
           renderlength = 1.2*zoomfac*renderlength
        case('+')
           xlength = 0.9/zoomfac*xlength
           ylength = 0.9/zoomfac*ylength
           renderlength = 0.9/zoomfac*renderlength
        case('o')
           if (iplotxarr(ipanel).le.ndim) then
              xcen = xorigin(iplotxarr(ipanel))
           else
              xcen = 0.
           endif
           if (iplotyarr(ipanel).le.ndim) then
              ycen = xorigin(iplotyarr(ipanel))
           else
              ycen = 0.
           endif
           print*,' centreing plot on origin x,y = ',xcen,ycen
        case('C')
           xcen = xpti
           ycen = ypti
        end select
        xmaxin = xmax(iplotxarr(ipanel))
        if (iamincolourbar .and. irenderarr(ipanel).gt.0) then
           !--rendering zoom does not allow pan - renderpt is always centre of axis
           renderpt = 0.5*(xmin(irenderarr(ipanel)) + xmax(irenderarr(ipanel)))
           xmin(irenderarr(ipanel)) = renderpt - 0.5*renderlength
           xmax(irenderarr(ipanel)) = renderpt + 0.5*renderlength
           print*,'zooming on colour bar: min, max = ',xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        else
           if (xpti.ge.xmin(iplotxarr(ipanel)) .and. xpti.le.xmax(iplotxarr(ipanel)) .and. ypti.le.xmax(iplotyarr(ipanel))) then
              xmin(iplotxarr(ipanel)) = xcen - 0.5*xlength
              xmax(iplotxarr(ipanel)) = xcen + 0.5*xlength
              print*,'zooming on x axis: min, max = ',xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
           if (ypti.ge.xmin(iplotyarr(ipanel)) .and. ypti.le.xmax(iplotyarr(ipanel)) .and. xpti.le.xmaxin) then
              xmin(iplotyarr(ipanel)) = ycen - 0.5*ylength
              xmax(iplotyarr(ipanel)) = ycen + 0.5*ylength
              print*,'zooming on y axis: min, max = ',xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
        endif

     case('a') ! adapt plot limits
        if (iamincolourbar .and. irenderarr(ipanel).gt.0) then
           print*,'adapting render limits ',xminadapt(irenderarr(ipanel)),xmaxadapt(irenderarr(ipanel))
           xmin(irenderarr(ipanel)) = xminadapt(irenderarr(ipanel))
           xmax(irenderarr(ipanel)) = xmaxadapt(irenderarr(ipanel))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        else
           if (xpti.ge.xmin(iplotxarr(ipanel)) .and. xpti.le.xmax(iplotxarr(ipanel)) &
              .and. ypti.le.xmax(iplotyarr(ipanel))) then
              print*,'adapting x limits ',xminadapt(iplotxarr(ipanel)),xmaxadapt(iplotxarr(ipanel))
              xmin(iplotxarr(ipanel)) = xminadapt(iplotxarr(ipanel))
              xmax(iplotxarr(ipanel)) = xmaxadapt(iplotxarr(ipanel))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
           if (ypti.ge.xmin(iplotyarr(ipanel)) .and. ypti.le.xmax(iplotyarr(ipanel)) &
              .and. xpti.le.xmax(iplotxarr(ipanel))) then
              print*,'adapting y limits ',xminadapt(iplotyarr(ipanel)),xmaxadapt(iplotyarr(ipanel))
              xmin(iplotyarr(ipanel)) = xminadapt(iplotyarr(ipanel))
              xmax(iplotyarr(ipanel)) = xmaxadapt(iplotyarr(ipanel))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
        endif
     !
     !--zoom in/out on vector plots (arrow size)
     !
     case('v')
        if (ivecarr(ipanel).gt.0) then
           print*,'decreasing vector arrow size'
           xmax(ivecarr(ipanel)) = 1.2*zoomfac*xmax(ivecarr(ipanel))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        endif
     case('V')
        if (ivecarr(ipanel).gt.0) then
           print*,'increasing vector arrow size'
           xmax(ivecarr(ipanel)) = 0.8/zoomfac*xmax(ivecarr(ipanel))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        endif
     case('w','W')
        if (ivecarr(ipanel).gt.0) then
           print*,'adapting vector arrow size'
           xmax(ivecarr(ipanel)) = -1.0
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        endif
     !
     !--set/unset log axes
     !
     case('l')
        !
        !--change colour bar, y and x itrans between log / not logged
        !
        if (iamincolourbar .and. irenderarr(ipanel).gt.0) then
           call change_itrans2(irenderarr(ipanel),xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel)),&
                               xminadapt(irenderarr(ipanel)),xmaxadapt(irenderarr(ipanel)))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        elseif (xpti.lt.xmin(iplotxarr(ipanel))) then
           if (iplotyarr(ipanel).le.ndim .and. irenderarr(ipanel).gt.0) then
              print "(a)",'error: cannot log coordinate axes with rendering'
           else
              call change_itrans2(iplotyarr(ipanel),xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)),&
                                  xminadapt(iplotyarr(ipanel)),xmaxadapt(iplotyarr(ipanel)))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
        elseif (ypti.lt.xmin(iplotyarr(ipanel))) then
           if (iplotxarr(ipanel).le.ndim .and. irenderarr(ipanel).gt.0) then
              print "(a)",'error: cannot log coordinate axes with rendering'
           else
              call change_itrans2(iplotxarr(ipanel),xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),&
                                  xminadapt(iplotxarr(ipanel)),xmaxadapt(iplotxarr(ipanel)))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
        endif
     !
     !--reset all range restrictions
     !
     case('R')
        call reset_ranges
        interactivereplot = .true.
        istep = istepnew
        iexit = .true.
     !
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        if (ipanel.gt.0) then
           call mvlegend(xpti,ypti,xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmax(iplotyarr(ipanel)))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        endif
     case('T') ! move title here
        if (ipanel.gt.0) then
           print*,'setting title position to current location...'
           call mvtitle(xpti,ypti,xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmax(iplotyarr(ipanel)))
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        endif
     case('H') ! move vector legend here
        if (ipanel.gt.0) then
           if (ivecarr(ipanel).gt.0) then
              print*,'setting vector plot legend to current location...'
              call mvlegendvec(xpti,ypti,xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmax(iplotyarr(ipanel)))
              istep = istepnew
              interactivereplot = .true.
              iexit = .true.
           endif
        endif
     case('m') ! change colour map (next scheme)
        call change_colourmap(icolourscheme,1)
        istep = istepnew
        interactivereplot = .true.
        iexit = .true.
     case('M') ! change colour map (previous scheme)
        call change_colourmap(icolourscheme,-1)
        istep = istepnew
        interactivereplot = .true.
        iexit = .true.
     case('i') ! invert colour map
        icolourscheme = -icolourscheme
        call change_colourmap(icolourscheme,0)
        istep = istepnew
        interactivereplot = .true.
        iexit = .true.
     case(achar(8)) ! delete plot annotation / colour bar (backspace)
        if (iamincolourbar .and. irenderarr(ipanel).gt.0) then
           iColourBarStyle = 0
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        elseif (xpti.lt.xmin(iplotxarr(ipanel)) .or. xpti.gt.xmax(iplotxarr(ipanel)) &
           .or. ypti.lt.xmin(iplotyarr(ipanel)) .or. ypti.gt.xmax(iplotyarr(ipanel))) then
           call deleteaxes()
           istep = istepnew
           interactivereplot = .true.
           iexit = .true.
        endif

     !
     !--timestepping
     !
     case('q','Q',achar(27),achar(3))
        iadvance = -666
        print*,'quitting...'
        iexit = .true.
     case('X','b','B') ! right click -> go back
!        iadvance = -abs(iadvance)
        istep = istepin - (istepjump)*istepsonpage - iadvance*istepsonpage
        lastpanel = 0
        iexit = .true.
     case('r') ! replot
        interactivereplot = .true.
        istep = istepnew
        iexit = .true.
     case(' ','n','N') ! space
        !iadvance = abs(iadvance)
        istep = istepin + (istepjump-1)*istepsonpage
        lastpanel = 0
        iexit = .true.
     case('0','1','2','3','4','5','6','7','8','9')
        read(char,*,iostat=ierr) istepjump
        if (ierr /=0) then
           print*,'*** internal error setting timestep jump' 
           istepjump = 1
        endif
        istepjump = int(zoomfac*istepjump)
        print*,' setting timestep jump = ',istepjump
     case(')')
        istepjump = int(zoomfac*10)
        print*,' setting timestep jump = ',istepjump
     !
     !--multiply everything by a factor of 10     
     !
     case('z','Z')
        zoomfac = 10.*zoomfac
        if (zoomfac.gt.1000000.) then
           zoomfac = 1.0
        endif
        print*,' LIMITS/TIMESTEPPING CHANGES NOW x ',zoomfac
     !
     !--unknown
     !
     case default
        print*,' x, y = ',xpti,ypti,'; unknown option "',trim(char),'" ',iachar(char)
     end select

     !
     !--save cursor position relative to the viewport
     !
     call get_vptxy(xpt,ypt,xcursor,ycursor,xmini,xmaxi,ymini,ymaxi)
     call reset_panel

     !
     !--do not let timestep go outside of bounds
     !  if we are at the first/last step, just print message and do nothing
     !  if iadvance trips over the bounds, jump to last/first step
     !
     if (iadvance.ne.-666 .and. iexit) then
        if (istep + iadvance .gt. ilaststep .and. iframe.eq.nframes) then
           print "(1x,a)",'reached last timestep'
           if (istepin.ne.ilaststep) then
              istep = ilaststep - istepsonpage*iadvance
           else
              istep = istepin
              iexit = .false.
           endif
        elseif (istep + iadvance .lt. 1 .and. ifirstframeonpage.eq.1) then
           print "(1x,a)",'reached first timestep: can''t go back'
           if (ifirststeponpage.ne.1) then
              istep = 1 - iadvance
           else
              istep = istepin
              iexit = .false.
           endif
        endif
     endif
  enddo interactive_loop
  return
  
  contains
  
   !---
   ! utility which translates between world co-ordinates (x,y)
   ! and viewport co-ordinates (relative to the whole viewport)
   !---
   !subroutine get_vptxy(x,y,vptx,vpty)
   ! implicit none
   ! real, intent(in) :: x,y
   ! real, intent(out) :: vptx,vpty
   ! real :: vptxmini,vptxmaxi,vptymini,vptymaxi
   !
   ! call pgqvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
   ! vptx = vptxmini + (x-xmini)/(xmaxi-xmini)*(vptxmaxi-vptxmini)
   ! vpty = vptymini + (y-ymini)/(ymaxi-ymini)*(vptymaxi-vptymini)
   !
   !end subroutine get_vptxy

   !--------
   ! utility to return which panel we are in given a point on the viewport
   ! and the viewport limits for each panel.
   !--------
   integer function getpanel(vptx,vpty)
    implicit none
    real, intent(in) :: vptx,vpty
    real :: vptxmini,vptxmaxi,vptymini,vptymaxi
    integer :: i,icol
    
    getpanel = 0
    icol = 0
    do i=1,size(vptxmin)
       icol = icol + 1
       if (icol.gt.nacross) icol = 1
       if (icol.gt.1) then
          ! if column>1 assign panel by being to the right of previous panel
          vptxmini = vptxmax(i-1)+barwmulti(i-1)
       else
          vptxmini = -0.1 ! allow for some error
       endif
       !--if last column extend xmax to right of page
       if (icol.eq.nacross) then
          vptxmaxi = 1.1
       elseif (verticalbar) then ! otherwise use max of current panel + space containing the colour bar
          vptxmaxi = vptxmax(i) + barwmulti(i)
       else
          vptxmaxi = vptxmax(i)
       endif
       
       !--if first row extend ymax to top of page
       if (i.le.nacross) then
          vptymaxi = 1.1
       else
          vptymaxi = vptymax(i)
       endif
       !--if last row then allow ymin to extend to bottom of page
       if (i.gt.(size(vptxmin)-nacross)) then
          vptymini = -0.1
       ! if not last row assign panel by being above row below
       elseif (i+nacross.le.size(vptxmin)) then
          vptymini = vptymax(i+nacross)
       elseif (.not.verticalbar) then
          vptymini = vptymin(i) - barwmulti(i)
       else
          vptymini = vptymin(i)
       endif
       if (vptx.gt.vptxmini .and. vptx.lt.vptxmaxi .and. &
           vpty.gt.vptymini .and. vpty.lt.vptymaxi) then
          if (getpanel.ne.0) print*,'Warning: multiple matching panels found ',getpanel,i
          getpanel = i
       endif
    enddo
    if (getpanel.eq.0) then
       print*,'Error determining panel: assuming last'
       getpanel = size(vptxmin)
    endif

    end function getpanel
    
   !--------
   ! utility to return x,y coordinates in a given panel given viewport coords
   !--------
    
    subroutine getxy(vptx,vpty,x,y,ipanel)
     implicit none
     real, intent(in) :: vptx,vpty
     real, intent(out) :: x,y
     integer, intent(in) :: ipanel
     
     if (ipanel.gt.0) then
        x = xmin(iplotxarr(ipanel)) + (vptx-vptxmin(ipanel))/(vptxmax(ipanel)-vptxmin(ipanel)) &
                          *(xmax(iplotxarr(ipanel))-xmin(iplotxarr(ipanel)))
        y = xmin(iplotyarr(ipanel)) + (vpty-vptymin(ipanel))/(vptymax(ipanel)-vptymin(ipanel)) &
                          *(xmax(iplotyarr(ipanel))-xmin(iplotyarr(ipanel)))
     endif
     
     return   
    end subroutine getxy

   !--------
   ! utility to reset the drawing surface so we can draw in a panel
   !--------
    
    subroutine set_panel(ipanel)
     implicit none
     integer, intent(in) :: ipanel
     
     if (ipanel.gt.0) then
        call pgswin(xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
        !--really should save viewport setting here, but doesn't matter 
        !  so long as interactive mode is the last thing called
        call pgsvp(vptxmin(ipanel),vptxmax(ipanel),vptymin(ipanel),vptymax(ipanel))
     endif
     
     return   
    end subroutine set_panel
    
    subroutine reset_panel
     implicit none
     
     call pgswin(xmini,xmaxi,ymini,ymaxi)
     
    end subroutine reset_panel

end subroutine interactive_multi

!------------------------------------------------------------
! utility which translates between world co-ordinates (x,y)
! and viewport co-ordinates (relative to the whole viewport)
!------------------------------------------------------------
subroutine get_vptxy(x,y,vptx,vpty,xmini,xmaxi,ymini,ymaxi)
 implicit none
 real, intent(in) :: x,y
 real, intent(out) :: vptx,vpty
 real, intent(in) :: xmini,xmaxi,ymini,ymaxi
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi

 call pgqvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
 vptx = vptxmini + (x-xmini)/(xmaxi-xmini)*(vptxmaxi-vptxmini)
 vpty = vptymini + (y-ymini)/(ymaxi-ymini)*(vptymaxi-vptymini)

end subroutine get_vptxy

!------------------------------------------------------------
! utility to return x,y coordinates given viewport coords
! (only works for single-panelled plots)
!------------------------------------------------------------
subroutine get_posxy(vptx,vpty,x,y,xmini,xmaxi,ymini,ymaxi)
 implicit none
 real, intent(in) :: vptx,vpty
 real, intent(out) :: x,y
 real, intent(in) :: xmini,xmaxi,ymini,ymaxi
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi

 call pgqvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
 x = xmini + (vptx-vptxmini)/(vptxmaxi-vptxmini)*(xmaxi-xmini)
 y = ymini + (vpty-vptymini)/(vptymaxi-vptymini)*(ymaxi-ymini)

 return   
end subroutine get_posxy

!-----------------------------------------------------------
! These subroutines interface to the actual plot settings
!-----------------------------------------------------------

!
!--plot a label showing the particle ID on the plot
!
subroutine plot_number(i,xi,yi)
 implicit none
 integer, intent(in) :: i
 real, intent(in) :: xi,yi
 integer :: nc
 real :: charheight
 character(len=20) :: string

 !--convert number to text string
 call pgnumb(i,0,1,string,nc)
 !--query and store character height
 call pgqch(charheight)
 !--change character height
 call pgsch(2.0)
 !--plot text string
 call pgtext(xi,yi,string(1:nc))
 !--reset character height
 call pgsch(charheight)

 return
end subroutine plot_number


subroutine deleteaxes()
 use settings_page, only:iaxis,iPlotLegend,& !iPlotStepLegend, &
                    iPlotTitles,iPlotScale 
 use settings_vecplot, only:iVecplotLegend
 implicit none

 if (iaxis.eq.-2) then
    !
    !-- would be better to do this properly by
    !   determining whether or not the cursor is over
    !   the legend, shape, title or whatever annotation the user
    !   wishes to be deleted. Instead at the moment we just
    !   delete the legends once the axes are already gone, and
    !   then in a somewhat arbitrary order.
    !
    iVecplotLegend = .false.
    if (iPlotLegend) then
       iPlotLegend = .false.
    elseif (iPlotTitles) then
       iPlotTitles = .false.
    elseif (iPlotScale) then
       iPlotScale = .false.
    endif
 elseif (iaxis.le.2 .and. iaxis.gt.-2) then
    iaxis = iaxis - 1
 elseif (iaxis.gt.2) then
    iaxis = -1
 elseif (iaxis.lt.-2) then
    iaxis = -2
 endif

end subroutine deleteaxes
!
!--move the legend to the current position
!
subroutine mvlegend(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hposlegend,vposlegend,fjustlegend,iPlotLegend
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 iPlotLegend = .true.
 hposlegend = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vposlegend = (ymax - yi)/ych
! !--automatically change justification
! if (hposlegend < 0.25) then
    fjustlegend = 0.0
! elseif (hposlegend > 0.75) then
!    fjustlegend = 1.0
! else
!    fjustlegend = 0.5
! endif
 print*,'hpos = ',hposlegend,' vpos = ',vposlegend,' just = ',fjustlegend
 
 return
end subroutine mvlegend
!
!--move the vector legend to the current position
!
subroutine mvlegendvec(xi,yi,xmin,xmax,ymax)
 use settings_vecplot, only:hposlegendvec,vposlegendvec,iVecplotLegend
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 iVecplotLegend = .true.
 hposlegendvec = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vposlegendvec = (ymax - yi)/ych
 print*,'hpos = ',hposlegendvec,' vpos = ',vposlegendvec
 
 return
end subroutine mvlegendvec
!
!--move the title to the current position
!
subroutine mvtitle(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hpostitle,vpostitle,fjusttitle,iPlotTitles
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 iPlotTitles = .true.
 hpostitle = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vpostitle = (yi - ymax)/ych

 !--automatically change justification
 if (hpostitle < 0.25) then
    fjusttitle = 0.0
 elseif (hpostitle > 0.75) then
    fjusttitle = 1.0
 else
    fjusttitle = 0.5
 endif
 print*,'hpos = ',hpostitle,' vpos = ',vpostitle,' just = ',fjusttitle
 
 return
end subroutine mvtitle

!
!--saves current plot limits
!
subroutine save_limits(iplot,xmin,xmax,setlim2)
 use limits, only:lim,lim2
 use multiplot, only:itrans
 use settings_data, only:ndim
 use settings_limits, only:iadapt,iadaptcoords
 use transforms, only:transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax
 logical, intent(in), optional :: setlim2
 logical :: uselim2
 real :: xmintemp,xmaxtemp
 
 uselim2 = .false.
 if (present(setlim2)) uselim2 = setlim2
 
 if (itrans(iplot).ne.0) then
    xmintemp = xmin
    xmaxtemp = xmax
    call transform_limits_inverse(xmintemp,xmaxtemp,itrans(iplot))
    if (uselim2) then
       lim2(iplot,1) = xmintemp
       lim2(iplot,2) = xmaxtemp
    else
       lim(iplot,1) = xmintemp
       lim(iplot,2) = xmaxtemp
    endif
 else
    if (uselim2) then
       lim2(iplot,1) = xmin
       lim2(iplot,2) = xmax
    else
       lim(iplot,1) = xmin
       lim(iplot,2) = xmax
    endif
 endif
 !
 !--change appropriate plot limits to fixed (not adaptive)
 !
 if (iplot.le.ndim) then
    iadaptcoords = .false.
 else
    iadapt = .false.
 endif
 
 return
end subroutine save_limits

!
!--implements parameter range restriction
!
subroutine restrict_range(iplot,xmin,xmax)
 use limits, only:range
 use multiplot, only:itrans
 use transforms, only:transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax
 real :: xmintemp,xmaxtemp
 
 if (itrans(iplot).ne.0) then
    xmintemp = xmin
    xmaxtemp = xmax
    call transform_limits_inverse(xmintemp,xmaxtemp,itrans(iplot))
    range(iplot,1) = xmintemp
    range(iplot,2) = xmaxtemp
 else
    range(iplot,1) = xmin
    range(iplot,2) = xmax
 endif
 
 return
end subroutine restrict_range

!
!--interface to routine which removes all parameter range restrictions
!
subroutine reset_ranges()
 use limits, only:reset_all_ranges
 implicit none
 
 call reset_all_ranges()
 
 return
end subroutine reset_ranges

!
!--interface to routine which resets second set of limits
!
subroutine reset_limits2(icol)
 use limits, only:reset_lim2
 implicit none
 integer, intent(in) :: icol
 
 call reset_lim2(icol)
 
 return
end subroutine reset_limits2

!
!--saves current plot limits for particle tracking
!
subroutine save_limits_track(iplot,xmin,xmax,xi)
 use multiplot, only:itrans
 use settings_data, only:ndim
 use settings_limits, only:xminoffset_track,xmaxoffset_track
 use transforms, only:transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax,xi
 real :: xmintemp,xmaxtemp
 
 if (iplot.gt.ndim) then
    print*,'ERROR in save_limits_track: iplot>ndim'
    return
 elseif (itrans(iplot).ne.0) then
    xmintemp = xmin
    xmaxtemp = xmax
    call transform_limits_inverse(xmintemp,xmaxtemp,itrans(iplot))
    xminoffset_track(iplot) = abs(xi - xmintemp)
    xminoffset_track(iplot) = abs(xmaxtemp - xi)
 else
    xminoffset_track(iplot) = abs(xi - xmin)
    xmaxoffset_track(iplot) = abs(xmax - xi)
 endif
 
 return
end subroutine save_limits_track
!
!--recalculates radius
!
subroutine save_itrackpart_recalcradius()
 use filenames, only:nsteps,nstepsinfile,ifileopen
 use settings_data, only:numplot,DataIsBuffered,iCalcQuantities
 use calcquantities, only:calc_quantities
 use labels, only:irad
 implicit none

 if (iCalcQuantities .and. irad.gt.0 .and. irad.le.numplot) then
    print "(a)",' Recalculating radius relative to tracked particle'
    if (DataIsBuffered) then
       call calc_quantities(1,nsteps)
    else
       call calc_quantities(1,nstepsinfile(ifileopen))
    endif
 endif
  
 return
end subroutine save_itrackpart_recalcradius
!
!--toggles log/unlog
!  note this only changes a pure log transform: will not change combinations
!
subroutine change_itrans(iplot,xmin,xmax)
 use multiplot, only:itrans
 use settings_data, only:numplot
 use transforms, only:transform_limits,transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(inout) :: xmin, xmax
 
 if (iplot.le.numplot) then
    if (itrans(iplot).eq.1) then
       itrans(iplot) = 0
       !!--untransform the plot limits
       call transform_limits_inverse(xmin,xmax,1)
    else
       itrans(iplot) = 1
       !!--transform the plot limits
       call transform_limits(xmin,xmax,itrans(iplot))
    endif
 endif
 
end subroutine change_itrans

subroutine change_itrans2(iplot,xmin,xmax,xmina,xmaxa)
 use multiplot, only:itrans
 use settings_data, only:numplot
 use transforms, only:transform_limits,transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(inout) :: xmin, xmax, xmina, xmaxa
 
 if (iplot.le.numplot) then
    if (itrans(iplot).eq.1) then
       itrans(iplot) = 0
       !!--untransform the plot limits
       call transform_limits_inverse(xmin,xmax,1)
       call transform_limits_inverse(xmina,xmaxa,1)
    else
       itrans(iplot) = 1
       !!--transform the plot limits
       call transform_limits(xmin,xmax,itrans(iplot))
       call transform_limits(xmina,xmaxa,itrans(iplot))
    endif
 endif
 
end subroutine change_itrans2

!
!--saves rotation options
!
subroutine save_rotation(ndim,anglexi,angleyi,anglezi)
 use settings_xsecrot, only:anglex,angley,anglez
 implicit none
 integer, intent(in) :: ndim
 real, intent(in) :: anglexi,angleyi,anglezi
 
 anglez = anglezi
 if (ndim.ge.3) then
    anglex = anglexi
    angley = angleyi
 endif
 
 return
end subroutine save_rotation

!
!--saves cross section position
!
subroutine save_xsecpos(xsecpos,xsec)
 use settings_xsecrot, only:xsecpos_nomulti,xsec_nomulti
 implicit none
 real, intent(in) :: xsecpos
 logical, intent(in) :: xsec
 
 xsecpos_nomulti = xsecpos
 xsec_nomulti = xsec
 
 return
end subroutine save_xsecpos

!
!--saves 3D perspective
!
subroutine save_perspective(zpos,dz)
 use settings_xsecrot, only:zobserver,dzscreenfromobserver
 implicit none
 real, intent(in) :: zpos,dz
 
 zobserver = zpos
 dzscreenfromobserver = dz
 
 return
end subroutine save_perspective
!
!--saves 3D opacity
!
subroutine save_opacity(taupartdepthi)
 use settings_xsecrot, only:taupartdepth
 implicit none
 real, intent(in) :: taupartdepthi
 
 taupartdepth = taupartdepthi
 
 return
end subroutine save_opacity
!
!--saves circles of interaction
!
subroutine save_circles(ncircpartset,icircpartset)
 use settings_part, only:ncircpart,icircpart
 implicit none
 integer, intent(in) :: ncircpartset
 integer, intent(in), dimension(:) :: icircpartset
 integer :: imax
 
 imax = min(size(icircpartset),size(icircpart),ncircpartset)
 ncircpart = imax
 icircpart(1:imax) = icircpartset(1:imax)
 print*,'saving ',imax,' circles of interaction only'
 
end subroutine save_circles
!
!--change colour map
!
subroutine change_colourmap(imap,istep)
 use colours, only:colour_set,ncolourschemes
 implicit none
 integer, intent(inout) :: imap
 integer, intent(in) :: istep
 
 imap = imap + istep
 if (abs(imap).gt.ncolourschemes) imap = 1
 if (abs(imap).lt.1) imap = ncolourschemes
 call colour_set(imap)
 
end subroutine change_colourmap

end module interactive_routines
