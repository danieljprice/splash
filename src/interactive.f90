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
!  Copyright (C) 2005-2019 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

module interactive_routines
 use colourbar, only:barisvertical,incolourbar,incolourbarlabel,adjustcolourbar
 implicit none
 public :: interactive_part,interactive_step,interactive_multi
 public :: set_movie_mode,save_limits
 private :: mvlegend,mvtitle,save_rotation
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
  xcoords,ycoords,zcoords,hi,icolourpart,iamtype,usetype,npartoftype,xmin,xmax,ymin,ymax, &
  rendermin,rendermax,renderminadapt,rendermaxadapt,contmin,contmax,&
  contminadapt,contmaxadapt,vecmax, &
  anglex,angley,anglez,ndim,xorigin,x_sec,zslicepos,dzslice, &
  zobserver,dscreen,use3Dopacity,rkappa,double_rendering,irerender,itrackpart,icolourscheme, &
  iColourBarStyle,labelrender,iadvance,istep,ilaststep,iframe,nframes,interactivereplot)
 use settings_xsecrot, only:setsequenceend
 use shapes,           only:inshape,edit_shape,edit_textbox,delete_shape,nshapes,add_shape_interactive
 use multiplot,        only:itrans
 use labels,           only:is_coord,ix,get_sink_type
 use limits,           only:assert_sensible_limits
 use settings_render,  only:projlabelformat,iapplyprojformat
 use settings_data,    only:ndataplots,ntypes,icoords,icoordsnew
 use plotlib,          only:plot_qwin,plot_curs,plot_sfs,plot_circ,plot_line,plot_pt1, &
                             plot_rect,plot_band,plot_sfs,plot_qcur,plot_left_click,plot_right_click,&
                             plot_scroll_left,plot_scroll_right,plotlib_is_pgplot,&
                             plot_shift_click,plot_lcur,plot_poly
 use params,           only:int1,maxparttypes
 use part_utils,       only:igettype
 use particleplots,    only:plot_kernel_gr
 use legends,          only:in_legend
 integer, intent(in) :: npart,icontour,ndim,iplotz,ivecx,ivecy,istep,ilaststep,iframe,nframes
 integer, intent(inout) :: irender,iColourBarStyle
 integer, intent(inout) :: iplotx,iploty,itrackpart,icolourscheme
 integer, intent(inout) :: iadvance
 integer, dimension(npart), intent(inout) :: icolourpart
 real, dimension(npart), intent(in) :: xcoords,ycoords,zcoords,hi
 integer(kind=int1), dimension(:), intent(in) :: iamtype
 logical, dimension(maxparttypes), intent(in) :: usetype
 integer, dimension(maxparttypes), intent(in) :: npartoftype
 real, intent(inout) :: xmin,xmax,ymin,ymax,rendermin,rendermax,vecmax,contmin,contmax,rkappa
 real, intent(inout) :: anglex,angley,anglez,zslicepos,dzslice,zobserver,dscreen
 real, intent(in) :: renderminadapt,rendermaxadapt,contminadapt,contmaxadapt
 real, intent(in), dimension(ndim) :: xorigin
 character(len=*), intent(inout) :: labelrender
 logical, intent(inout) :: x_sec
 logical, intent(out) :: irerender,interactivereplot
 logical, intent(in) :: use3Dopacity, double_rendering
 real,    parameter :: pi=4.*atan(1.)
 integer, parameter :: maxpts = 64
 integer :: i,iclosest,iclosestsink,ierr,ixsec,ishape,ilegend,itype,npts
 integer :: nmarked,ncircpart,itrackparttemp,iadvancenew,itypesink
 integer, dimension(1000) :: icircpart
 real :: xpt,ypt
 real :: xpt2,ypt2,xcen,ycen,xminwin,xmaxwin,yminwin,ymaxwin
 real :: xptmin,xptmax,yptmin,yptmax,zptmin,zptmax,rptmax2
 real :: rmin,rminsink,rr,gradient,yint,dx,dy,dr,anglerad
 real :: xlength,ylength,renderlength,renderpt,zoomfac,scalefac
 real :: dxlength,dylength,xmaxin,ymaxin,contlength
 real, dimension(4)      :: xline,yline
 real, dimension(maxpts) :: xpts,ypts
 character(len=1) :: char,char2
 logical :: iexit, rotation, verticalbar, iamincolourbar, mixedtypes, use3Dperspective
 logical :: iadvanceset, leftclick, iselectpoly, iselectcircle
 logical, save :: print_help = .true.
 logical, save :: in_movie_mode = .false.

 if (plot_qcur()) then
    if (.not.print_help) print*,'entering interactive mode...press h in plot window for help'
    !print*, plot_left_click
 else
    !print*,'cannot enter interactive mode: device has no cursor'
    return
 endif

 mixedtypes = size(iamtype) >= npart
 use3Dperspective = abs(dscreen) > tiny(0.)
 iadvanceset = .false.

 char = 'A'
 xline = 0.
 yline = 0.
 !
 !--convert saved cursor position (saved in viewport coords)
 !  back to coordinates
 !
 call plot_qwin(xminwin,xmaxwin,yminwin,ymaxwin)
 call get_posxy(xcursor,ycursor,xpt,ypt,xminwin,xmaxwin,yminwin,ymaxwin)

!  xpt = 0.
!  ypt = 0.
 xpt2 = 0.
 ypt2 = 0.
 xmaxin = xmax
 ymaxin = ymax
 zoomfac = 1.0
 scalefac = 1.1
 ncircpart = 0
 itrackparttemp = itrackpart
 iexit = .false.
 rotation = .false.
 irerender = .false.
 interactivereplot = .false.
 if (is_coord(iplotx,ndim) .and. is_coord(iploty,ndim) .and. ndim >= 2) rotation = .true.
 verticalbar = barisvertical(iColourBarStyle)

 if (iplotz > 0 .and. x_sec) then
    zptmin = zslicepos - 0.5*dzslice
    zptmax = zslicepos + 0.5*dzslice
 else
    !--if not using z range, make it encompass all the particles
    zptmin = -huge(zptmin)
    zptmax = huge(zptmax)
 endif

 interactiveloop: do while (.not.iexit)
    if (print_help) then
       char = 'h'
       ierr = 0
       print_help = .false.
    else
       ierr = plot_curs(xpt,ypt,char)
    endif
    !
    !--exit if the device is not interactive
    !
    if (ierr==1) return
    !print*,'x,y = ',xpt,ypt,' function = ',char,iachar(char)

    !
    !--find closest particle
    !
    rmin = huge(rmin)
    rminsink = huge(rmin)
    iclosest = 0
    iclosestsink = 0
    itypesink = get_sink_type(ntypes)
    xlength = xmax - xmin
    ylength = ymax - ymin
    if (xlength > tiny(xlength)) then
       dxlength = 1./xlength
    else
       dxlength = 0.
    endif
    if (ylength > tiny(ylength)) then
       dylength = 1./ylength
    else
       dylength = 0.
    endif
    itype = 1
    over_npart: do i=1,npart
       if (ntypes > 1) then
          if (mixedtypes) then
             itype = int(iamtype(i))
          else
             itype = igettype(i,npartoftype)
          endif
          if (.not.usetype(itype)) cycle over_npart
       endif
       !--use distance normalised on screen
       rr = ((xcoords(i)-xpt)*dxlength)**2 + ((ycoords(i)-ypt)*dylength)**2
       if (rr < rmin) then
          iclosest = i
          rmin = rr
       endif
       if (itype==itypesink .and. rr < rminsink) then
          iclosestsink = i
          rminsink = rr
       endif
    enddo over_npart

    !--query the position of the colour bar
    iamincolourbar = incolourbar(iColourBarStyle,4,xpt,ypt,xmin,xmax,ymin,ymax)

    select case(char)
       !
       !--particle plot stuff
       !
    case('p')
       if (iclosest > 0 .and. iclosest <= npart) then
          print*,' closest particle = ',iclosest,'x = ',xcoords(iclosest),' y =',ycoords(iclosest)
          call plot_number(iclosest,xcoords(iclosest),ycoords(iclosest))
       else
          print*,'error: could not determine closest particle'
       endif
    case('c')
       if (iclosest > 0 .and. iclosest <= npart) then
          print*,'plotting circle of interaction on particle ',iclosest, &
                  ' h = ',hi(iclosest)
          !--save settings for these
          ncircpart = ncircpart + 1
          if (ncircpart > size(icircpart)) then
             print*,'WARNING: ncircles > array limits, cannot save'
             ncircpart = size(icircpart)
          else
             icircpart(ncircpart) = iclosest
          endif
          call plot_sfs(2)
          if (icoordsnew /= icoords) then
             call plot_kernel_gr(icoordsnew,icoords,iplotx,iploty,&
                   xcoords(iclosest),ycoords(iclosest),zcoords(iclosest),2.*hi(iclosest))
          else
             call plot_circ(xcoords(iclosest),ycoords(iclosest),2.*hi(iclosest))
          endif
       else
          print*,'error: could not determine closest particle'
       endif
    case('t')
       !--track closest particle (must save to activate)
       if (itrackpart /= 0 .and. itrackparttemp==itrackpart) then
          itrackpart = 0
          itrackparttemp = 0
          print*,' particle tracking limits OFF'
       else
          if (iclosest > 0 .and. iclosest <= npart) then
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
       call plot_pt1(xpt,ypt,4)
       !--select second point
       print*,' select another point (using left click or g) to plot line '
       ierr = plot_band(1,1,xline(2),yline(2),xline(3),yline(3),char2)
       !--draw line if left click or g
       select case(char2)
       case(plot_left_click,'g')
          print*
          !--mark second point
          call plot_pt1(xline(3),yline(3),4)
          xlength = xline(3)-xline(2)
          if (abs(xlength) < tiny(xlength)) then
             xline(1) = xline(2)
             xline(4) = xline(2)
             yline(1) = ymin
             yline(4) = ymax
             print*,' error: gradient = infinite'
          elseif (xline(2) < xline(3)) then
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
          if (abs(xlength) > tiny(xlength)) then
             gradient = ylength/xlength
             yint = yline(3) - gradient*xline(3)
             print*,' gradient = ',gradient,' y intercept = ',yint
             yline(1) = gradient*xline(1) + yint
             yline(4) = gradient*xline(4) + yint
          endif
          !--plot line joining the two points
          call plot_line(4,xline,yline)
       case default
          print*,' action cancelled'
       end select
       !
       !--help
       !
    case('h')
       print "(/,a)",' -------------- interactive mode commands -----------------------------'
       print*,' SPACE BAR                : skip to next timestep/file'
       print*,' 0,1,2,3..9 and click     : go forward/back n timesteps (back=r.click)'
       print*,' left click (or A)        : zoom/select'
       print*,' right click (or X or b)  : previous timestep'
       print*,' shift+left click         : IRREGULAR particle selection'
       print*,' left click on colour bar : change rendering limits'
       print*,' +/-      : zoom IN/OUT (_ for out by 20%) '
       print*,' a        : (a)dapt plot limits (inside box, over axes or colour bar)'
       print*,' l        : (l)og / unlog axis  (over x/y axis or colour bar)'
       print*,' o/C/n    : re-centre plot on (o)rigin/(C)ursor/(n)earest sink position'
       print*,' backspace: delete annotation  (over axes, legend, title or colour bar)'
       print*,' r        : (r)eplot current plot'
       print*,' R        : (R)eset/remove all range restrictions'
       print*,' p/c      : label closest (p)article/plot (c)ircle of interaction'
       print*,' t        : t)rack closest particle/turn tracking off (coord plots only)'
       print*,' g        : plot a line and find its g)radient'
       print*,' ctrl-t, ^: add text or arrow(^) annotation at current position'
       print*,' ctrl-m   : toggle Hollywood mode'
       print*,' G/T/H    : move le(G)end, (T)itle or (H) vector legend to current position'
       print*,' m/M/i    : change colour map (m=next,M=previous,i=invert) (rendered plots only)'
       if (irender > 0) print*,' f/F      : f)lip to next/previous column in rendering (rendered plots only)'
       print*,' v/V/w    : decrease/increase/adapt arrow size on vector plots (Z for x10)'
       if (ndim >= 3) then
          print*,' k/K      : decrease/increase opacity on opacity-rendered plots (Z for x10)'
       endif
       print*,' e/E      : use current frame/settings as end point to animation sequence'
       if (ndim > 1) then
          print*,' , . < >  : rotate about z axis by +(-) 15,30 degrees (coord plots only)'
          if (ndim >= 3) then
             print*,' [ ] { }  : rotate about x axis by +/- 15,30 degrees (coord plots only)'
             print*,' / \ ? |  : rotate about y axis by +/- 15,30 degrees (coord plots only)'
             print*,' x        : take cross section (coord plots only)'
             if (iplotz > 0) then
                print*,' u/U/d/D  : move cross section/perspective pos. up/down (towards/away from observer)'
             endif
          endif
       endif
       print*,' s        : (s)ave current settings for all steps'
       print*,' q/Q/esc  : (q)uit plotting'
       print*,' z/Z(oom) : timestepping, zoom and limits-changing multiplied by 10'
       print*,'----------------------------------------------------------------------'
    case('s','S')
       itrackpart = itrackparttemp
       if (itrackpart<=0 .or. itrackpart > size(xcoords)) then
          call save_limits(iplotx,xmin,xmax)
          call save_limits(iploty,ymin,ymax)
          call save_itrackpart_recalcradius(itrackpart) ! set saved value to zero
       else
          print*,'tracking particle ',itrackpart,'x,y = ',xcoords(itrackpart),ycoords(itrackpart)
          if (is_coord(iplotx,ndim)) then
             call save_limits_track(iplotx,xmin,xmax,xcoords(itrackpart))
          else
             call save_limits(iplotx,xmin,xmax)
          endif
          if (is_coord(iploty,ndim)) then
             call save_limits_track(iploty,ymin,ymax,ycoords(itrackpart))
          else
             call save_limits(iploty,ymin,ymax)
          endif
          call save_itrackpart_recalcradius(itrackpart)
       endif
       if (irender > 0) call save_limits(irender,rendermin,rendermax)
       if (icontour > 0) then
          if (icontour==irender &
               .and. abs(rendermin-contmin) <= tiny(0.) &
               .and. abs(rendermax-contmax) <= tiny(0.)) then
             call reset_limits2(icontour)
          elseif (icontour==irender .and. .not.double_rendering) then
             call save_limits(icontour,contmin,contmax,setlim2=.true.)
          else
             call save_limits(icontour,contmin,contmax)
          endif
       endif
       if (ivecx > 0 .and. ivecy > 0) then
          call save_limits(ivecx,-vecmax,vecmax)
          call save_limits(ivecy,-vecmax,vecmax)
       endif
       if (ndim==3 .and. iplotz > 0 .and. irender > 0 .and. use3Dopacity) then
          call save_opacity(rkappa)
       endif
       if (ncircpart > 0) call save_circles(ncircpart,icircpart)
       if (rotation) call save_rotation(ndim,anglex,angley,anglez)
       if (iplotz > 0) then
          if (x_sec) then
             call save_xsecpos(zslicepos,x_sec)
          else
             call save_perspective(zobserver,dscreen)
          endif
       endif
       call save_windowsize()
       print*,'> interactively set limits saved <'
       !
       !--actions on left click
       !
    case(plot_left_click,plot_right_click,plot_shift_click) ! left click
       !
       !--change colour bar limits
       !
       leftclick = (char == plot_left_click)
       ishape = inshape(xpt,ypt,itrans(iplotx),itrans(iploty),xmin,xmax,ymin,ymax)
       if (ishape > 0) then
          call edit_shape(ishape,xpt,ypt,itrans(iplotx),itrans(iploty),first=.false.)
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       elseif (iamincolourbar .and. irender > 0 .and. leftclick) then
          if (incolourbarlabel(iColourBarStyle,4,xpt,ypt,xmin,xmax,ymin,ymax)) then
             if (verticalbar) then
                call edit_textbox(xpt,ypt,90.,0.,labelrender)
                projlabelformat = trim(labelrender)
                iapplyprojformat = irender
             else
                call edit_textbox(xpt,ypt,0.,0.,labelrender)
                projlabelformat = trim(labelrender)
                iapplyprojformat = irender
             endif
             iadvance = 0
             interactivereplot = .true.
             iexit = .true.
          else
             print*,'click to set rendering limits'
             if (verticalbar) then
                ierr = plot_band(3,1,xpt,ypt,xpt2,ypt2,char2)
             else
                ierr = plot_band(4,1,xpt,ypt,xpt2,ypt2,char2)
             endif
             if (char2 == plot_left_click) then
                if (double_rendering) then
                   call adjustcolourbar(iColourBarStyle,xpt,ypt,xpt2,ypt2,&
                                         xmin,xmax,ymin,ymax,contmin,contmax)
                   !print*,'setting doublerender min = ',contmin
                   !print*,'setting doublerender max = ',contmax
                else
                   call adjustcolourbar(iColourBarStyle,xpt,ypt,xpt2,ypt2,&
                                         xmin,xmax,ymin,ymax,rendermin,rendermax)
                   !print*,'setting render min = ',rendermin
                   !print*,'setting render max = ',rendermax
                endif
                iadvance = 0
                interactivereplot = .true.
                iexit = .true.
             endif
          endif
          !
          !--zoom or mark particles
          !
       else
          print*,'select area: '
          !--Note: circle selection is not implemented in PGPLOT
          iselectpoly = (char==plot_shift_click .or.((.not.leftclick).and.plotlib_is_pgplot))
          iselectcircle = (char==plot_right_click .and. .not.plotlib_is_pgplot)
          if (iselectpoly) then
             if (plotlib_is_pgplot) then
                print*,'left click/a)dd points; middle click/d)elete points; X/x)finish'
             else
                print*,'left click/a)dd points; middle click/d)elete points; q)uit/abort'
                if (irender <= 0) print*,'1-9 = close polygon and mark particles with colours 1-9'
                print*,'0 = close polygon and hide selected particles'
                print*,'p = close polygon and plot selected particles only'
                print*,'c = close polygon and plot circles of interaction on selected parts'
             endif
          else
             print*,'left click : zoom'
             if (irender <= 0) print*,'1-9 = mark selected particles with colours 1-9'
             print*,'0 = hide selected particles'
             print*,'p = plot selected particles only'
             print*,'c = plot circles of interaction on selected parts'
          endif
          if (leftclick .or. iselectcircle) then
             print*,'x = use particles within x parameter range only'
             print*,'y = use particles within y parameter range only'
             print*,'r = use particles within x and y parameter range only'
             print*,'R = remove all range restrictions'
          endif

          npts = 1
          xpts(1) = xpt   ! to avoid problems with uninitialised variables
          ypts(1) = ypt
          if (iselectpoly) then
             call plot_lcur(maxpts,npts,xpts,ypts,char2)
             if (plotlib_is_pgplot) then
                if (irender <= 0) print*,'1-9 = close polygon and mark particles with colours 1-9'
                print*,'0 = close polygon and hide selected particles'
                print*,'p = close polygon and plot selected particles only'
                print*,'c = close polygon and plot circles of interaction on selected parts'
                ierr = plot_band(0,1,xpt,ypt,xpt2,ypt2,char2)
             endif
             xptmin = minval(xpts(1:npts))
             xptmax = maxval(xpts(1:npts))
             yptmin = minval(ypts(1:npts))
             yptmax = maxval(ypts(1:npts))
          elseif (iselectcircle) then
             ierr = plot_band(8,1,xpt,ypt,xpt2,ypt2,char2)
             rptmax2 = (xpt2-xpt)**2 + (ypt2-ypt)**2
             rr = sqrt(rptmax2)
             xptmin = xpt - rr
             xptmax = xpt + rr
             yptmin = ypt - rr
             yptmax = ypt + rr
          else ! left click: rectangle selection
             ierr = plot_band(2,1,xpt,ypt,xpt2,ypt2,char2)
             xptmin = min(xpt,xpt2)
             xptmax = max(xpt,xpt2)
             yptmin = min(ypt,ypt2)
             yptmax = max(ypt,ypt2)
          endif
          if (.not.(iselectcircle.or.iselectpoly) .or. char2==plot_left_click) then ! rectangle selection
             !print*,'xrange = ',xptmin,'->',xptmax
             !print*,'yrange = ',yptmin,'->',yptmax
             if (iplotz /= 0 .and. x_sec) then
                print*,'(zrange = ',zptmin,'->',zptmax,')'
             endif
          endif

          select case (char2)
          case(plot_left_click) ! zoom if another left click
             call plot_sfs(2)
             !--draw the selected shape in case zooming is slow
             if (iselectpoly) then
                call plot_poly(npts,xpts,ypts)
             elseif (iselectcircle) then
                call plot_circ(xpt,ypt,sqrt(rptmax2))
             else
                call plot_rect(xpt,xpt2,ypt,ypt2)
             endif
             !--zoom is identical for rectangle, poly and circle selection
             xmin = xptmin
             xmax = xptmax
             ymin = yptmin
             ymax = yptmax
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          case('0','1','2','3','4','5','6','7','8','9') ! mark particles
             nmarked = 0
             if (irender <= 0 .or. char2=='0' .or. char2=='1') then
                do i=1,npart
                   if (inslice(zcoords(i),zptmin,zptmax) .and. &
                        (leftclick .and. inrectangle(xcoords(i),ycoords(i),xptmin,xptmax,yptmin,yptmax) &
                    .or.(iselectcircle .and. incircle(xcoords(i)-xpt,ycoords(i)-ypt,rptmax2)) &
                    .or.(iselectpoly .and. inpoly(xcoords(i),ycoords(i),xpts,ypts,npts)))) then
                      read(char2,*,iostat=ierr) icolourpart(i)
                      if (ierr /=0) then
                         print*,'*** error marking particle'
                         icolourpart(i) = 1
                         !--translate 0 to icolourpart = -1 for non-plotted particles
                      elseif (icolourpart(i)==0) then
                         icolourpart(i) = -1
                      endif
                      nmarked = nmarked + 1
                   endif
                enddo
                print*,'marked ',nmarked,' particles in selected region'
             endif
             iadvance = 0
             if (nmarked > 0) irerender = .true.
             interactivereplot = .true.
             iexit = .true.
          case('p') ! plot selected particles only
             nmarked = 0
             do i=1,npart
                if (inslice(zcoords(i),zptmin,zptmax) .and. &
                     (leftclick .and. inrectangle(xcoords(i),ycoords(i),xptmin,xptmax,yptmin,yptmax) &
                 .or.(iselectcircle .and. incircle(xcoords(i)-xpt,ycoords(i)-ypt,rptmax2)) &
                 .or.(iselectpoly .and. inpoly(xcoords(i),ycoords(i),xpts,ypts,npts)))) then
                   nmarked = nmarked + 1
                   if (icolourpart(i) <= 0) icolourpart(i) = 1
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
                if (inslice(zcoords(i),zptmin,zptmax) .and. &
                     (leftclick .and. inrectangle(xcoords(i),ycoords(i),xptmin,xptmax,yptmin,yptmax) &
                 .or.(iselectcircle .and. incircle(xcoords(i)-xpt,ycoords(i)-ypt,rptmax2)) &
                 .or.(iselectpoly .and. inpoly(xcoords(i),ycoords(i),xpts,ypts,npts)))) then
                   if (ncircpart < size(icircpart)) then
                      ncircpart = ncircpart + 1
                      icircpart(ncircpart) = i
                      call plot_sfs(2)
                      call plot_circ(xcoords(i),ycoords(i),2.*hi(i))
                   endif
                endif
             enddo
             print*,'set ',ncircpart,' circles of interaction in selected region'
             if (ncircpart==size(icircpart)) print*,' (first ',size(icircpart),' only)'
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
    case('-','_','+','o','C','n','N') ! zoom in/out
       xlength = xmax - xmin
       ylength = ymax - ymin
       xcen = 0.5*(xmax + xmin)
       ycen = 0.5*(ymax + ymin)
       renderlength = rendermax - rendermin
       contlength = contmax - contmin
       select case(char)
       case('-')
          xlength = zoomfac*scalefac*xlength
          ylength = zoomfac*scalefac*ylength
          renderlength = zoomfac*scalefac*renderlength
          contlength = zoomfac*scalefac*contlength
       case('_')
          xlength = 1.2*zoomfac*scalefac*xlength
          ylength = 1.2*zoomfac*scalefac*ylength
          renderlength = 1.2*zoomfac*scalefac*renderlength
          contlength = 1.2*zoomfac*scalefac*contlength
       case('+')
          xlength = xlength/(zoomfac*scalefac)
          ylength = ylength/(zoomfac*scalefac)
          renderlength = renderlength/(zoomfac*scalefac)
          contlength = contlength/(zoomfac*scalefac)
       case('o')
          if (itrackpart > 0) then
             print*,'centreing limits on tracked particle ',itrackpart,'x,y = ',xcoords(itrackpart),ycoords(itrackpart)
             xcen = xcoords(itrackpart)
             ycen = ycoords(itrackpart)
          else
             if (is_coord(iplotx,ndim)) then
                xcen = xorigin(iplotx)
             else
                xcen = 0.
             endif
             if (is_coord(iploty,ndim)) then
                ycen = xorigin(iploty)
             else
                ycen = 0.
             endif
             print*,' centreing plot on origin x,y = ',xcen,ycen
          endif
       case('C')
          xcen = xpt
          ycen = ypt
       case('n','N')
          if (iclosestsink > 0) then
             xcen = xcoords(iclosestsink)
             ycen = ycoords(iclosestsink)
             print*,'centreing limits on sink particle ',iclosestsink,'x,y = ',xcen,ycen
          else
             print*,'error: could not find closest sink particle, using origin instead'
             xcen = 0.
             ycen = 0.
          endif
       end select
       if (iamincolourbar .and. irender > 0) then
          !--rendering zoom does not allow pan - renderpt is always centre of axis
          if (double_rendering) then
             renderpt = 0.5*(contmin + contmax)
             contmin = renderpt - 0.5*contlength
             contmax = renderpt + 0.5*contlength
             call assert_sensible_limits(contmin,contmax)
             !print*,'zooming on colour bar: min, max = ',contmin,contmax
          else
             renderpt = 0.5*(rendermin + rendermax)
             rendermin = renderpt - 0.5*renderlength
             rendermax = renderpt + 0.5*renderlength
             call assert_sensible_limits(rendermin,rendermax)
             !print*,'zooming on colour bar: min, max = ',rendermin,rendermax
          endif
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       else
          if (xpt >= xmin .and. xpt <= xmax .and. ypt <= ymaxin) then
             xmin = xcen - 0.5*xlength
             xmax = xcen + 0.5*xlength
             call assert_sensible_limits(xmin,xmax)
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          endif
          if (ypt >= ymin .and. ypt <= ymax .and. xpt <= xmaxin) then
             ymin = ycen - 0.5*ylength
             ymax = ycen + 0.5*ylength
             call assert_sensible_limits(ymin,ymax)
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          endif
       endif
       if (char=='o') then
          xpt = xcen
          ypt = ycen
       endif
    case('a') ! reset plot limits
       if (iamincolourbar .and. irender > 0) then
          if (double_rendering) then
             contmin = contminadapt
             contmax = contmaxadapt
             call assert_sensible_limits(contmin,contmax)
          else
             rendermax = rendermaxadapt
             if (itrans(irender)==1 .and. renderminadapt < rendermaxadapt-5.) then ! if logged
                if (abs(rendermin-(rendermaxadapt-4.)) > epsilon(0.)) then
                   rendermin = rendermaxadapt - 4. ! if logged, do not give 20 orders of mag
                   print "(a)",' *** MIN SET 4 DEX FROM MAX, PRESS ''a'' AGAIN TO GIVE FULL RANGE ***'
                else
                   rendermin = renderminadapt
                endif
             else
                rendermin = renderminadapt
             endif
             call assert_sensible_limits(rendermin,rendermax)
          endif
          iadvance = 0              ! that it should change the render limits
          interactivereplot = .true.
          iexit = .true.
       else
          xmaxin = xmax
          ymaxin = ymax
          if (xpt >= xmin .and. xpt <= xmax .and. ypt <= ymaxin) then
             call adapt_limits_interactive('x',npart,xcoords,xmin,xmax,icolourpart,iamtype,usetype)
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          endif
          if (ypt >= ymin .and. ypt <= ymax .and. xpt <= xmaxin) then
             call adapt_limits_interactive('y',npart,ycoords,ymin,ymax,icolourpart,iamtype,usetype)
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
       if (ivecx > 0 .and. ivecy > 0) then
          !print*,'decreasing vector arrow size'
          vecmax = zoomfac*scalefac*vecmax
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
    case('V')
       if (ivecx > 0 .and. ivecy > 0) then
          !print*,'increasing vector arrow size'
          vecmax = vecmax/(zoomfac*scalefac)
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
    case('w','W')
       if (ivecx > 0 .and. ivecy > 0) then
          !print*,'adapting vector arrow size'
          vecmax = -1.0
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
       !
       !--change opacity on 3D opacity rendered plots
       !
    case('k')
       if (ndim==3 .and. iplotz > 0 .and. irender /= 0 .and. use3Dopacity) then
          print*,'decreasing opacity by factor of ',1.5*zoomfac*scalefac
          rkappa = rkappa/(zoomfac*scalefac)
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('K')
       if (ndim==3 .and. iplotz > 0 .and. irender /= 0 .and. use3Dopacity) then
          print*,'increasing opacity by factor of ',1.5*zoomfac*scalefac
          rkappa = rkappa*(zoomfac*scalefac)
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
       if (iamincolourbar .and. irender > 0) then
          if (double_rendering) then
             !call change_itrans(irender,rendermin,rendermax)
             call change_itrans(icontour,contmin,contmax)
          else
             if (icontour==irender) then
                call change_itrans2(irender,rendermin,rendermax,contmin,contmax)
             else
                call change_itrans(irender,rendermin,rendermax)
             endif
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       elseif (xpt < xmin) then
          if (is_coord(iploty,ndim) .and. irender > 0) then
             print "(a)",'error: cannot log coordinate axes with rendering'
          else
             call change_itrans(iploty,ymin,ymax)
             iadvance = 0
             interactivereplot = .true.
             iexit = .true.
          endif
       elseif (ypt < ymin) then
          if (is_coord(iplotx,ndim) .and. irender > 0) then
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
                          anglex,angley,anglez,zobserver,use3Dopacity,rkappa, &
                          x_sec,zslicepos,xmin,xmax,ymin,ymax,rendermin,rendermax)
       !
       !--rotation
       !
    case(',')
       if (rotation) then
          !print*,'changing z rotation angle by -15 degrees...'
          if (int(scalefac) > 1) then
             anglez = anglez - int(scalefac)
          else
             anglez = anglez - 15.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('<')
       if (rotation) then
          !print*,'changing z rotation angle by -30 degrees...'
          if (int(scalefac) > 1) then
             anglez = anglez - int(scalefac)
          else
             anglez = anglez - 30.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('.')
       if (rotation) then
          !print*,'changing z rotation angle by 15 degrees...'
          if (int(scalefac) > 1) then
             anglez = anglez + int(scalefac)
          else
             anglez = anglez + 15.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('>')
       if (rotation) then
          !print*,'changing z rotation angle by 30 degrees...'
          if (int(scalefac) > 1) then
             anglez = anglez + int(scalefac)
          else
             anglez = anglez + 30.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('/')
       if (rotation .and. ndim >= 2) then
          !print*,'changing y rotation angle by -15 degrees...'
          if (int(scalefac) > 1) then
             angley = angley - int(scalefac)
          else
             angley = angley - 15.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('?')
       if (rotation .and. ndim >= 2) then
          !print*,'changing y rotation angle by -30 degrees...'
          if (int(scalefac) > 1) then
             angley = angley - int(scalefac)
          else
             angley = angley - 30.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('\')
       if (rotation .and. ndim >= 2) then
          !print*,'changing y rotation angle by 15 degrees...'
          if (int(scalefac) > 1) then
             angley = angley + int(scalefac)
          else
             angley = angley + 15.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('|')
       if (rotation .and. ndim >= 2) then
          !print*,'changing y rotation angle by 30 degrees...'
          !print*,'changing y rotation angle by 15 degrees...'
          if (int(scalefac) > 1) then
             angley = angley + int(scalefac)
          else
             angley = angley + 30.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('[')
       if (rotation .and. ndim >= 3) then
          !print*,'changing x rotation angle by -15 degrees...'
          if (int(scalefac) > 1) then
             anglex = anglex - int(scalefac)
          else
             anglex = anglex - 15.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('{')
       if (rotation .and. ndim >= 3) then
          !print*,'changing x rotation angle by -30 degrees...'
          if (int(scalefac) > 1) then
             anglex = anglex - int(scalefac)
          else
             anglex = anglex - 30.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case(']')
       if (rotation .and. ndim >= 3) then
          !print*,'changing x rotation angle by 15 degrees...'
          if (int(scalefac) > 1) then
             anglex = anglex + int(scalefac)
          else
             anglex = anglex + 15.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('}')
       if (rotation .and. ndim >= 3) then
          !print*,'changing x rotation angle by 30 degrees...'
          if (int(scalefac) > 1) then
             anglex = anglex + int(scalefac)
          else
             anglex = anglex + 30.
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
       !
       !--set cross section position
       !
    case('x')
       if (rotation .and. ndim >= 3) then
          xline(1) = xpt
          yline(1) = ypt
          !--work out which is the third dimension
          do i=1,3
             if (i /= iplotx .and. i /= iploty) ixsec = i
          enddo
          print*,' select cross section position (using left click or x)'
          ierr = plot_band(1,1,xline(1),yline(1),xline(2),yline(2),char2)
          !--work out cross section if left click or x again
          select case(char2)
          case(plot_left_click,'x')
             !--plot the cross section line
             call plot_line(2,xline(1:2),yline(1:2))
             !--work out angle with the x axis
             !  and offset of line from origin
             dx = xline(2) - xline(1)
             dy = yline(2) - yline(1)
             anglerad = atan2(dy,dx)
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
             zslicepos = yint/cos(anglerad)
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
       if (iplotz > 0 .and. ndim==3) then
          if (x_sec) then
             if (int(scalefac) > 1) dzslice = scalefac
             print*,'shifting cross section position up by ',dzslice
             zslicepos = zslicepos + dzslice
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          elseif (use3Dperspective) then
             if (abs(zobserver) < tiny(0.)) then
                print*,'resetting z position'
                zobserver = 1.
             else
                print*,'shifting perspective position up by factor of ',scalefac
                zobserver = scalefac*zoomfac*zobserver
             endif
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          endif
       endif
    case('U') ! move cross section up by 2*dxsec
       if (iplotz > 0 .and. ndim==3) then
          if (x_sec) then
             print*,'shifting cross section position up by ',2.*dzslice
             zslicepos = zslicepos + 2.*dzslice
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          elseif (use3Dperspective) then
             if (abs(zobserver) < tiny(0.)) then
                print*,'resetting z position'
                zobserver = 1.
             else
                print*,'shifting perspective position up by factor of 2'
                zobserver = 2.*scalefac*zoomfac*zobserver
             endif
             iadvance = 0
             interactivereplot = .true.
             irerender = .true.
             iexit = .true.
          endif
       endif
    case('d') ! move cross section down by dxsec
       if (iplotz > 0 .and. ndim==3) then
          if (x_sec) then
             if (int(scalefac) > 1) dzslice = scalefac
             print*,'shifting cross section position down by ',dzslice
             zslicepos = zslicepos - dzslice
          elseif (use3Dperspective) then
             if (abs(zobserver) < tiny(0.)) then
                print*,'resetting z position'
                zobserver = 1.
             else
                print*,'shifting perspective position down'
                zobserver = zobserver/(zoomfac*scalefac)
             endif
          endif
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       endif
    case('D') ! move cross section down by 2*dxsec
       if (iplotz > 0 .and. ndim==3) then
          if (x_sec) then
             print*,'shifting cross section position down by ',2.*dzslice
             zslicepos = zslicepos - 2.*dzslice
          elseif (use3Dperspective) then
             if (abs(zobserver) < tiny(0.)) then
                print*,'resetting z position'
                zobserver = 1.
             else
                print*,'shifting perspective position down by factor of 2'
                zobserver = zobserver/(2.*zoomfac*scalefac)
             endif
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
       if (ivecx > 0 .and. ivecy > 0) then
          print*,'setting vector plot legend to current location...'
          call mvlegendvec(xpt,ypt,xmin,xmax,ymax)
       endif
       iadvance = 0
       interactivereplot = .true.
       iexit = .true.
    case('f') ! change rendered quantity (next)
       if (irender /= 0 .and. ndim > 0) then
          irender = irender + 1
          if (irender > ndataplots) irender = 1
          if (is_coord(irender,ndim)) irender = ix(ndim) + 1
          !if (irender==ndataplots+1) irender = 0
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       else  ! flip y axis
          iploty = iploty + 1
          if (iploty > ndataplots) iploty = 1
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
    case('F') ! change rendered quantity (previous)
       if (irender /= 0 .and. ndim > 0) then
          irender = irender - 1
          if (is_coord(irender,ndim)) irender = ix(1) - 1
          if (irender < 1) irender = ndataplots
          !if (irender==ndim) irender = 0
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       else  ! flip y axis
          iploty = iploty - 1
          if (iploty < 1) iploty = ndataplots
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
    case(achar(13))
       if (in_movie_mode) then
          call unset_movie_mode()
          in_movie_mode = .false.
       else
          call set_movie_mode(.true.)
          in_movie_mode = .true.
       endif
       iadvance = 0
       interactivereplot = .true.
       irerender = .true.
       iexit = .true.
    case('m') ! change colour map (next scheme)
       call change_colourmap(icolourscheme,1)
       iadvance = 0
       interactivereplot = .true.
       !irerender = .true.
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
    case('^') ! add arrow
       call add_shape_interactive(xpt,ypt,itrans(iplotx),itrans(iploty),0,ierr,shape_type=3)
       if (ierr==0) then
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
    case(achar(20)) ! add text shape
       call add_shape_interactive(xpt,ypt,itrans(iplotx),itrans(iploty),0,ierr,shape_type=6)
       if (ierr==0) then
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
    case(achar(8)) ! delete plot annotation / colour bar (backspace)
       ishape = inshape(xpt,ypt,itrans(iplotx),itrans(iploty),xmin,xmax,ymin,ymax)
       if (ishape > 0) then
          call delete_shape(ishape,nshapes)
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       elseif (iamincolourbar .and. irender > 0) then
          iColourBarStyle = 0
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       elseif (xpt < xmin .or. xpt > xmax .or. ypt < ymin .or. ypt > ymax) then
          call deleteaxes()
          iadvance = 0
          interactivereplot = .true.
          irerender = .true.
          iexit = .true.
       else
          ilegend = in_legend(xpt,ypt)
          if (legend_is_plotted(ilegend)) then
             call delete_legend(ilegend)
             iadvance = 0
             interactivereplot = .true.
             iexit = .true.
          else
             print*,' nothing to delete at x,y =',xpt,',',ypt
          endif
       endif
       !
       !--timestepping
       !
    case('q','Q',achar(27),achar(3))
       iadvance = -666
       !print*,'quitting...'
       iexit = .true.
    case('b','B',plot_scroll_left) ! right click -> go back
       iadvance = -abs(iadvance)
       iexit = .true.
    case('r') ! replot
       iadvance = 0
       interactivereplot = .true.
       irerender = .true.
       iexit = .true.
    case(' ',plot_scroll_right) ! space
       iadvance = abs(iadvance)
       iexit = .true.
    case('0','1','2','3','4','5','6','7','8','9')
       read(char,*,iostat=ierr) iadvancenew
       if (ierr /=0) then
          print*,'*** internal error setting timestep jump'
          iadvancenew = 1
       endif
       if ((iadvance > 1 .or. iadvanceset) .and. iadvance <= 9999) then
          iadvance = 10*iadvance + iadvancenew
          if (iadvance > 9999) iadvance = 1
       elseif (iadvancenew==0) then
          iadvance = 10
       else
          iadvance = iadvancenew
       endif
       iadvance = int(zoomfac*iadvance)
       print*,' setting timestep jump / zoom factor = ',iadvance
       iadvanceset = .true.
       scalefac = iadvance
    case(')')
       iadvance = int(zoomfac*10)
       print*,' setting timestep jump = ',iadvance
       !
       !--multiply everything by a factor of 10
       !
    case('z','Z')
       zoomfac = 10.*zoomfac
       if (zoomfac > 1000000.) then
          zoomfac = 1.0
       endif
       print*,' LIMITS/TIMESTEPPING CHANGES NOW x ',zoomfac
       !
       !--unknown
       !
    case default
       if (iachar(char) >= iachar('a')) then
          print*,' x, y = ',xpt,ypt,'; unknown option "',trim(char), '" ',iachar(char)
       endif
       print*,' x, y = ',xpt,ypt,'; GOT ',iachar(char)
    end select

    !
    !--save cursor position relative to the viewport
    !
    call plot_qwin(xminwin,xmaxwin,yminwin,ymaxwin)
    call get_vptxy(xpt,ypt,xcursor,ycursor)

    if (rotation) then
       if (anglez >= 360.) anglez = anglez - 360.
       if (anglez < 0.) anglez = anglez + 360.
       if (ndim > 2) then
          if (angley >= 360.) angley = angley - 360.
          if (angley < 0.) angley = angley + 360.
          if (anglex >= 360.) anglex = anglex - 360.
          if (anglex < 0.) anglex = anglex + 360.
       endif
    endif
    !
    !--do not let timestep go outside of bounds
    !  if we are at the first/last step, just print message and do nothing
    !  if iadvance trips over the bounds, jump to last/first step
    !
    if (iadvance /= -666 .and. iexit) then
       if (istep + iadvance  >  ilaststep .and. iframe==nframes) then
          print "(1x,a)",'reached last timestep'
          if (ilaststep-istep  > 0) then
             iadvance= ilaststep - istep
          else
             iexit = .false.
          endif
       elseif (istep + iadvance  <  1 .and. iframe==1) then
          print "(1x,a)",'reached first timestep: can''t go back'
          if (1-istep  < 0) then
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
 use plotlib, only:plot_qcur,plot_curs,plot_band,plot_rect,plot_left_click
 integer, intent(inout) :: iadvance
 integer, intent(in) :: istep,ilaststep
 real, intent(inout) :: xmin,xmax,ymin,ymax
 logical, intent(out) :: interactivereplot
 integer :: ierr
 real :: xpt,ypt,xpt2,ypt2
 real :: xlength, ylength, zoomfac
 character(len=1) :: char,char2
 logical :: iexit

 if (plot_qcur()) then
    print*,'entering interactive mode...press h in plot window for help'
 else
    !print*,'cannot enter interactive mode: device has no cursor'
    return
 endif
 char = 'A'
 xpt = 0.
 ypt = 0.
 zoomfac = 1.0
 iexit = .false.
 interactivereplot = .false.

 do while (.not.iexit)
    ierr = plot_curs(xpt,ypt,char)
    !
    !--exit if the device is not interactive
    !
    if (ierr==1) return

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

    case(plot_left_click) ! left click
       !
       !--draw rectangle from the point and reset the limits
       !
       print*,'select area: '
       print*,'left click : zoom'
       ierr = plot_band(2,1,xpt,ypt,xpt2,ypt2,char2)
       print*,xpt,ypt,xpt2,ypt2,char2
       select case (char2)
       case(plot_left_click)   ! zoom if another left click
          call plot_rect(xpt,xpt2,ypt,ypt2)
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
       if (xpt >= xmin .and. xpt <= xmax .and. ypt <= ymax) then
          print*,'zooming on x axis'
          xmin = xpt - 0.5*xlength
          xmax = xpt + 0.5*xlength
          iadvance = 0
          interactivereplot = .true.
          iexit = .true.
       endif
       if (ypt >= ymin .and. ypt <= ymax .and. xpt <= xmax) then
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
       if (zoomfac > 1000000.) then
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
    if (iadvance /= -666 .and. iexit) then
       if (istep + iadvance  >  ilaststep) then
          print "(1x,a)",'reached last timestep'
          if (ilaststep-istep  > 0) then
             iadvance= ilaststep - istep
          else
             iexit = .false.
          endif
       elseif (istep + iadvance  <  1) then
          print "(1x,a)",'reached first timestep: can''t go back'
          if (1-istep  < 0) then
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
                             lastpanel,iplotxarr,iplotyarr,irenderarr,icontourarr,ivecarr,&
                             use_double_rendering,xmin,xmax,vptxmin,vptxmax,vptymin,vptymax, &
                             barwmulti,xminadapt,xmaxadapt,nacross,ndim,xorigin,icolourscheme, &
                             iColourBarStyle,interactivereplot)
 use labels,    only:is_coord,iamvec
 use limits,    only:assert_sensible_limits
 use multiplot, only:itrans
 use shapes,    only:add_shape_interactive,inshape,edit_shape,delete_shape,nshapes
 use plotlib,   only:plot_qcur,plot_band,plot_qwin,plot_pt1,plot_curs,plot_line,plot_left_click
 use legends,   only:in_legend
 integer, intent(inout) :: iadvance
 integer, intent(inout) :: istep,iframe,lastpanel,iColourBarStyle
 integer, intent(in) :: ifirststeponpage,ilaststep,nacross,ndim,ifirstframeonpage,nframes
 integer, intent(inout) :: icolourscheme
 integer, intent(in), dimension(:) :: iplotxarr,iplotyarr,irenderarr,icontourarr,ivecarr
 real, dimension(:), intent(in) :: vptxmin,vptxmax,vptymin,vptymax,barwmulti
 real, dimension(:), intent(inout) :: xmin,xmax,xminadapt,xmaxadapt
 real, intent(in), dimension(ndim) :: xorigin
 logical, intent(in)  :: use_double_rendering
 logical, intent(out) :: interactivereplot
 integer :: ierr,ipanel,ipanel2,istepin,istepnew,i,istepjump,istepsonpage,ishape,ilegend
 integer :: istepjumpnew,ivecx,ivecy
 real :: xpt,ypt,xpt2,ypt2,xpti,ypti,renderpt,xptmin,xptmax,yptmin,yptmax
 real :: xlength,ylength,renderlength,contlength,zoomfac,scalefac
 real :: vptxi,vptyi,vptx2i,vpty2i,vptxceni,vptyceni
 real :: xmini,xmaxi,ymini,ymaxi,xcen,ycen,gradient,dr,yint,xmaxin
 real, dimension(4) :: xline,yline
 character(len=1) :: char,char2
 logical :: iexit,iamincolourbar,verticalbar,double_render,istepjumpset
 logical, save :: print_help = .true.

 if (plot_qcur()) then
    if (.not.print_help) print*,'entering interactive mode...press h in plot window for help'
 else
    !print*,'cannot enter interactive mode: device has no cursor'
    return
 endif
 char = 'A'
 zoomfac = 1.0
 scalefac = 1.1
 xpt2 = 0.
 ypt2 = 0.
 verticalbar = barisvertical(iColourBarStyle)

 !
 !--convert saved cursor position (saved in viewport coords)
 !  back to world coordinates:
 !
 !--query window limits in world coords
 call plot_qwin(xmini,xmaxi,ymini,ymaxi)

 !--determine which plot the cursor falls on
 !print*,' saved xcursor,ycursor = ',xcursor,ycursor,vptxmin,vptxmax,vptymin,vptymax
 ipanel = getpanel(xcursor,ycursor)
 !print*,' saved panel = ',ipanel

 !--set the window to correspond to the panel we last left the cursor in
 call set_panel(ipanel)

 !--set the position in x,y relative to this panel
 call getxy(xcursor,ycursor,xpt,ypt,ipanel)
 !print*,' saved x,y = ',xpt,ypt
 call get_vptxy(xpt,ypt,vptxi,vptyi)
 !print*,'saved vptx,y = ',vptxi,vptyi,ipanel

 iexit = .false.
 interactivereplot = .false.
 istepin = istep
 istepnew = ifirststeponpage - iadvance
 istepsonpage = abs(istep - ifirststeponpage)/iadvance + 1
 istepjump = 1
 istepjumpset = .false.
!  print*,'istep = ',istepnew
!  print*,'steps on page = ',istepsonpage

 interactive_loop: do while (.not.iexit)
    if (print_help) then
       print_help = .false.
       char = 'h'
       ierr = 0
    else
       ierr = plot_curs(xpt,ypt,char)
    endif
    !
    !--exit if the device is not interactive
    !
    if (ierr==1) return

    !print*,'x, y = ',xpt,ypt,' function = ',char
    !
    !--determine which plot the cursor falls on
    !
    call get_vptxy(xpt,ypt,vptxi,vptyi)
    ipanel = getpanel(vptxi,vptyi)
    !print*,'xpt,ypt = ',xpt,ypt,vptxi,vptyi,ipanel

    !--translate vpt co-ords to x,y in current panel
    call getxy(vptxi,vptyi,xpti,ypti,ipanel)

    !print*,'vptx,y = ',xpti,ypti,vptxi,vptyi,ipanel
    !--query the position of the colour bar
    if (ipanel > 0) then
       if (barwmulti(ipanel) > tiny(barwmulti)) then
          iamincolourbar = incolourbar(iColourBarStyle,4,xpti,ypti,xmin(iplotxarr(ipanel)), &
                            xmax(iplotxarr(ipanel)),xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
       else
          !--for colour bars on tiled plots, use viewport coords
          iamincolourbar = incolourbar(iColourBarStyle,0,vptxi,vptyi,&
                            minval(vptxmin),maxval(vptxmax),minval(vptymin),maxval(vptymax))
       endif
    else
       iamincolourbar = .false.
    endif

    !--work out if this plot is double rendered or not
    double_render = (icontourarr(ipanel) > 0 .and. use_double_rendering)

    select case(char)
    case('h')
       print "(/,a)",' ------- interactive mode commands (multiple plots per page) --------'
       print*,' SPACE BAR (or n)         : skip to next timestep/file'
       print*,' 0,1,2,3..9 and click     : go forward/back n timesteps (back=r.click)'
       print*,' left click (or A)        : zoom/select'
       print*,' right click (or X or b)  : previous timestep'
       print*,' left click on colour bar : change rendering limits'
       print*,' +/-      : zoom IN/OUT (_ for out by 20%) '
       print*,' a        : (a)dapt plot limits (inside box, over axes/colour bar)'
       print*,' l        : (l)og / unlog axis  (over x/y axis or colour bar)'
       print*,' o/C      : re-centre plot on (o)rigin/(C)ursor position'
       print*,' backspace: delete annotation  (over axes/legend/title/colour bar)'
       print*,' r        : (r)eplot current plot'
       print*,' R        : (R)eset/remove all range restrictions'
       print*,' g        : plot a line and find its g)radient'
       print*,' ctrl-t, ^: add text or arrow(^) at current position'
       print*,' G/T/H    : move le(G)end, (T)itle or (H) vector legend to current position'
       print*,' m/M/i    : change colour map (m=next,M=previous,i=invert) (rendered plots only)'
       print*,' v/V/w    : decrease/increase/adapt arrow size on vector plots (Z for x10)'
       print*,' s        : (s)ave current settings for all steps'
       print*,' q/Q/esc  : (q)uit plotting'
       print*,' z/Z(oom) : timestepping, zoom and limits-changing multiplied by 10'
       print*,'--------------------------------------------------------------------'
    case('g')   ! draw a line between two points
       xline(2) = xpti
       yline(2) = ypti
       call set_panel(ipanel)
       !--mark first point
       call plot_pt1(xpti,ypti,4)
       !--select second point
       print*,' select another point (using left click or g) to plot line '
       ierr = plot_band(1,1,xline(2),yline(2),xline(3),yline(3),char2)
       !--draw line if left click or g
       select case(char2)
       case(plot_left_click,'g')
          print*
          !--mark second point
          call plot_pt1(xline(3),yline(3),4)
          xlength = xline(3)-xline(2)
          if (abs(xlength) < tiny(xlength)) then
             xline(1) = xline(2)
             xline(4) = xline(2)
             yline(1) = xmin(iplotyarr(ipanel))
             yline(4) = xmax(iplotyarr(ipanel))
             print*,' error: gradient = infinite'
          elseif (xline(2) < xline(3)) then
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
          if (abs(xlength) > tiny(xlength)) then
             gradient = ylength/xlength
             yint = yline(3) - gradient*xline(3)
             print*,' gradient = ',gradient,' y intercept = ',yint
             yline(1) = gradient*xline(1) + yint
             yline(4) = gradient*xline(4) + yint
          endif
          !--plot line joining the two points
          call plot_line(4,xline,yline)
          call reset_panel
       case default
          print*,' action cancelled'
       end select
    case('s','S')
       do i=1,size(vptxmin)
          call save_limits(iplotxarr(i),xmin(iplotxarr(i)),xmax(iplotxarr(i)))
          call save_limits(iplotyarr(i),xmin(iplotyarr(i)),xmax(iplotyarr(i)))
          if (irenderarr(i) > 0) call save_limits(irenderarr(i),xmin(irenderarr(i)),xmax(irenderarr(i)))
          if (icontourarr(i) > 0) call save_limits(icontourarr(i),xmin(icontourarr(i)),xmax(icontourarr(i)))
          if (ivecarr(i) > 0) then
             ivecx = iamvec(ivecarr(i)) + iplotxarr(i) - 1
             ivecy = iamvec(ivecarr(i)) + iplotyarr(i) - 1
             if (ivecx > 0) call save_limits(ivecx,-xmax(ivecx),xmax(ivecx))
             if (ivecy > 0) call save_limits(ivecy,-xmax(ivecy),xmax(ivecy))
          endif
       enddo
       call save_windowsize()
       print*,'> interactively set limits saved <'
    case(plot_left_click) ! left click
       !ishape = inshape(xpt,ypt,itrans(iplotx),itrans(iploty))
       !if (ishape > 0) then
      !    call edit_shape(ishape,xpt,ypt,itrans(iplotx),itrans(iploty),first=.false.)
      !    iadvance = 0
      !    interactivereplot = .true.
      !!    iexit = .true.
       !endif
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
       if (ipanel > 0 .and. iamincolourbar .and. irenderarr(ipanel) > 0) then
          print*,'click to set rendering limits'
          if (verticalbar) then
             ierr = plot_band(3,1,xpt,ypt,xpt2,ypt2,char2)
          else
             ierr = plot_band(4,1,xpt,ypt,xpt2,ypt2,char2)
          endif
          if (char2 == plot_left_click) then
             call get_vptxy(xpt2,ypt2,vptx2i,vpty2i)
             !--use centre point of first click and current click to
             !  better determine panel
             vptxceni = 0.5*(vptxi + vptx2i)
             vptyceni = 0.5*(vptyi + vpty2i)
             ipanel2 = getpanel(vptxceni,vptyceni)
             if (ipanel2 > 0 .and. ipanel2 /= ipanel) then
                print*,'panel = ',ipanel2,' was ',ipanel
                ipanel = ipanel2
             endif
             call getxy(vptx2i,vpty2i,xpt2,ypt2,ipanel)
             !--reset first point according to current panel
             call getxy(vptxi,vptyi,xpti,ypti,ipanel)
             double_render = (icontourarr(ipanel) > 0 .and. use_double_rendering)

             if (barwmulti(ipanel) > tiny(barwmulti)) then
                if (double_render) then
                   call adjustcolourbar(iColourBarStyle,vptxi,vptyi,vptx2i,vpty2i,&
                         vptxmin(ipanel),vptxmax(ipanel),vptymin(ipanel),vptymax(ipanel),&
                         xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel)))
                else
                   call adjustcolourbar(iColourBarStyle,vptxi,vptyi,vptx2i,vpty2i,&
                         vptxmin(ipanel),vptxmax(ipanel),vptymin(ipanel),vptymax(ipanel),&
                         xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel)))
                endif
             else
                !--for global colour bars (ie. on tiled plots) use viewport co-ordinates to set render limits
                if (double_render) then
                   call adjustcolourbar(iColourBarStyle,vptxi,vptyi,vptx2i,vpty2i,&
                         minval(vptxmin),maxval(vptxmax),minval(vptymin),maxval(vptymax),&
                         xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel)))
                else
                   call adjustcolourbar(iColourBarStyle,vptxi,vptyi,vptx2i,vpty2i,&
                         minval(vptxmin),maxval(vptxmax),minval(vptymin),maxval(vptymax),&
                         xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel)))
                endif
             endif
             if (double_render) then
                print*,'setting double-render min, max = ',xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel))
             else
                print*,'setting render min, max = ',xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel))
             endif
             istep = istepnew
             interactivereplot = .true.
             iexit = .true.
          endif
       else
          ierr = plot_band(2,1,xpt,ypt,xpt2,ypt2,char2)
          !call pgrect(xpt,xpt2,ypt,ypt2)
          call get_vptxy(xpt2,ypt2,vptx2i,vpty2i)
          !--use centre point of first click and current click to
          !  better determine panel
          vptxceni = 0.5*(vptxi + vptx2i)
          vptyceni = 0.5*(vptyi + vpty2i)
          ipanel2 = getpanel(vptxceni,vptyceni)
          if (ipanel2 > 0 .and. ipanel2 /= ipanel) then
             ipanel = ipanel2
             print*,'panel = ',ipanel
          endif
          if (ipanel <= 0) cycle interactive_loop
          call getxy(vptx2i,vpty2i,xpt2,ypt2,ipanel)
          !--reset first point according to current panel
          call getxy(vptxi,vptyi,xpti,ypti,ipanel)
          xptmin = min(xpti,xpt2)
          xptmax = max(xpti,xpt2)
          yptmin = min(ypti,ypt2)
          yptmax = max(ypti,ypt2)

          select case (char2)
          case(plot_left_click)   ! zoom if another left click
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
       if (ipanel <= 0) cycle interactive_loop
       xlength = xmax(iplotxarr(ipanel)) - xmin(iplotxarr(ipanel))
       ylength = xmax(iplotyarr(ipanel)) - xmin(iplotyarr(ipanel))
       xcen = 0.5*(xmax(iplotxarr(ipanel)) + xmin(iplotxarr(ipanel)))
       ycen = 0.5*(xmax(iplotyarr(ipanel)) + xmin(iplotyarr(ipanel)))
       if (irenderarr(ipanel) > 0) then
          renderlength = xmax(irenderarr(ipanel)) - xmin(irenderarr(ipanel))
       else
          renderlength = 0.
       endif
       if (icontourarr(ipanel) > 0) then
          contlength = xmax(icontourarr(ipanel)) - xmin(icontourarr(ipanel))
       else
          contlength = 0.
       endif
       select case(char)
       case('-')
          xlength = zoomfac*scalefac*xlength
          ylength = zoomfac*scalefac*ylength
          renderlength = zoomfac*scalefac*renderlength
          contlength = zoomfac*scalefac*contlength
       case('_')
          xlength = 1.2*zoomfac*scalefac*xlength
          ylength = 1.2*zoomfac*scalefac*ylength
          renderlength = 1.2*zoomfac*scalefac*renderlength
          contlength = 1.2*zoomfac*scalefac*contlength
       case('+')
          xlength = xlength/(zoomfac*scalefac)
          ylength = ylength/(zoomfac*scalefac)
          renderlength = renderlength/(zoomfac*scalefac)
          contlength = contlength/(zoomfac*scalefac)
       case('o')
          if (is_coord(iplotxarr(ipanel),ndim)) then
             xcen = xorigin(iplotxarr(ipanel))
          else
             xcen = 0.
          endif
          if (is_coord(iplotyarr(ipanel),ndim)) then
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
       if (iamincolourbar .and. irenderarr(ipanel) > 0) then
          if (double_render) then
             !--rendering zoom does not allow pan - renderpt is always centre of axis
             renderpt = 0.5*(xmin(icontourarr(ipanel)) + xmax(icontourarr(ipanel)))
             xmin(icontourarr(ipanel)) = renderpt - 0.5*contlength
             xmax(icontourarr(ipanel)) = renderpt + 0.5*contlength
             call assert_sensible_limits(xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel)))
             !print*,'zooming on colour bar: min, max = ',xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel))
          else
             !--rendering zoom does not allow pan - renderpt is always centre of axis
             renderpt = 0.5*(xmin(irenderarr(ipanel)) + xmax(irenderarr(ipanel)))
             xmin(irenderarr(ipanel)) = renderpt - 0.5*renderlength
             xmax(irenderarr(ipanel)) = renderpt + 0.5*renderlength
             call assert_sensible_limits(xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel)))
             !print*,'zooming on colour bar: min, max = ',xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel))
          endif
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       else
          if (xpti >= xmin(iplotxarr(ipanel)) .and. xpti <= xmax(iplotxarr(ipanel)) .and. ypti <= xmax(iplotyarr(ipanel))) then
             xmin(iplotxarr(ipanel)) = xcen - 0.5*xlength
             xmax(iplotxarr(ipanel)) = xcen + 0.5*xlength
             call assert_sensible_limits(xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)))
             !print*,'zooming on x axis: min, max = ',xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel))
             istep = istepnew
             interactivereplot = .true.
             iexit = .true.
          endif
          if (ypti >= xmin(iplotyarr(ipanel)) .and. ypti <= xmax(iplotyarr(ipanel)) .and. xpti <= xmaxin) then
             xmin(iplotyarr(ipanel)) = ycen - 0.5*ylength
             xmax(iplotyarr(ipanel)) = ycen + 0.5*ylength
             call assert_sensible_limits(xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
             !print*,'zooming on y axis: min, max = ',xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel))
             istep = istepnew
             interactivereplot = .true.
             iexit = .true.
          endif
       endif

    case('a') ! adapt plot limits
       if (iamincolourbar .and. irenderarr(ipanel) > 0) then
          if (double_render) then
             !print*,'adapting double-render limits ',xminadapt(icontourarr(ipanel)),xmaxadapt(icontourarr(ipanel))
             xmin(icontourarr(ipanel)) = xminadapt(icontourarr(ipanel))
             xmax(icontourarr(ipanel)) = xmaxadapt(icontourarr(ipanel))
             call assert_sensible_limits(xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel)))
          else
             !print*,'adapting render limits ',xminadapt(irenderarr(ipanel)),xmaxadapt(irenderarr(ipanel))
             xmax(irenderarr(ipanel)) = xmaxadapt(irenderarr(ipanel))
             if (itrans(irenderarr(ipanel))==1 .and. &
                 xminadapt(irenderarr(ipanel)) < xmaxadapt(irenderarr(ipanel))-5.) then ! if logged
                if (abs(xmin(irenderarr(ipanel))-(xmaxadapt(irenderarr(ipanel)) - 4.)) > epsilon(0.)) then
                   xmin(irenderarr(ipanel)) = xmaxadapt(irenderarr(ipanel)) - 4.
                   print "(a)",' *** MIN SET 4 DEX FROM MAX, PRESS ''a'' AGAIN TO GIVE FULL RANGE ***'
                else
                   xmin(irenderarr(ipanel)) = xminadapt(irenderarr(ipanel))
                endif
             else
                xmin(irenderarr(ipanel)) = xminadapt(irenderarr(ipanel))
             endif
             call assert_sensible_limits(xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel)))
          endif
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       else
          !--save xmax before we go changing it so can check the y axis
          xmaxin = xmax(iplotxarr(ipanel))
          if (xpti >= xmin(iplotxarr(ipanel)) .and. xpti <= xmax(iplotxarr(ipanel)) &
              .and. ypti <= xmax(iplotyarr(ipanel))) then
             print*,'adapting x limits ',xminadapt(iplotxarr(ipanel)),xmaxadapt(iplotxarr(ipanel))
             xmin(iplotxarr(ipanel)) = xminadapt(iplotxarr(ipanel))
             xmax(iplotxarr(ipanel)) = xmaxadapt(iplotxarr(ipanel))
             call assert_sensible_limits(xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)))
             istep = istepnew
             interactivereplot = .true.
             iexit = .true.
          endif
          if (ypti >= xmin(iplotyarr(ipanel)) .and. ypti <= xmax(iplotyarr(ipanel)) &
              .and. xpti <= xmaxin) then
             print*,'adapting y limits ',xminadapt(iplotyarr(ipanel)),xmaxadapt(iplotyarr(ipanel))
             xmin(iplotyarr(ipanel)) = xminadapt(iplotyarr(ipanel))
             xmax(iplotyarr(ipanel)) = xmaxadapt(iplotyarr(ipanel))
             call assert_sensible_limits(xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
             istep = istepnew
             interactivereplot = .true.
             iexit = .true.
          endif
       endif
       !
       !--zoom in/out on vector plots (arrow size)
       !
    case('v')
       if (ivecarr(ipanel) > 0) then
          !print*,'decreasing vector arrow size'
          xmax(ivecarr(ipanel)) = zoomfac*scalefac*xmax(ivecarr(ipanel))
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       endif
    case('V')
       if (ivecarr(ipanel) > 0) then
          !print*,'increasing vector arrow size'
          xmax(ivecarr(ipanel)) = xmax(ivecarr(ipanel))/(zoomfac*scalefac)
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       endif
    case('w','W')
       if (ivecarr(ipanel) > 0) then
          !print*,'adapting vector arrow size'
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
       if (iamincolourbar .and. irenderarr(ipanel) > 0) then
          if (double_render) then
             call change_itrans2(icontourarr(ipanel),xmin(icontourarr(ipanel)),xmax(icontourarr(ipanel)),&
                               xminadapt(icontourarr(ipanel)),xmaxadapt(icontourarr(ipanel)))
          else
             call change_itrans2(irenderarr(ipanel),xmin(irenderarr(ipanel)),xmax(irenderarr(ipanel)),&
                               xminadapt(irenderarr(ipanel)),xmaxadapt(irenderarr(ipanel)))
          endif
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       elseif (xpti < xmin(iplotxarr(ipanel))) then
          if (is_coord(iplotyarr(ipanel),ndim) .and. irenderarr(ipanel) > 0) then
             print "(a)",'error: cannot log coordinate axes with rendering'
          else
             call change_itrans2(iplotyarr(ipanel),xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)),&
                                  xminadapt(iplotyarr(ipanel)),xmaxadapt(iplotyarr(ipanel)))
             istep = istepnew
             interactivereplot = .true.
             iexit = .true.
          endif
       elseif (ypti < xmin(iplotyarr(ipanel))) then
          if (is_coord(iplotxarr(ipanel),ndim) .and. irenderarr(ipanel) > 0) then
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
       if (ipanel > 0) then
          call mvlegend(xpti,ypti,xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmax(iplotyarr(ipanel)),ipanel)
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       endif
    case('T') ! move title here
       if (ipanel > 0) then
          print*,'setting title position to current location...'
          call mvtitle(xpti,ypti,xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmax(iplotyarr(ipanel)))
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       endif
    case('H') ! move vector legend here
       if (ipanel > 0) then
          if (ivecarr(ipanel) > 0) then
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
    case('^') ! add arrow shape
       call set_panel(ipanel)
       print*,' adding arrow in panel ',ipanel
       call add_shape_interactive(xpti,ypti,itrans(iplotxarr(ipanel)),itrans(iplotyarr(ipanel)),ipanel,ierr,shape_type=3)
       if (ierr==0) then
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       endif
    case(achar(20)) ! add text shape
       call set_panel(ipanel)
       print*,' adding text in panel ',ipanel
       call add_shape_interactive(xpti,ypti,itrans(iplotxarr(ipanel)),itrans(iplotyarr(ipanel)),ipanel,ierr)
       if (ierr==0) then
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       endif
    case(achar(8)) ! delete plot annotation / colour bar (backspace)
       ishape = inshape(xpti,ypti,itrans(iplotxarr(ipanel)),itrans(iplotxarr(ipanel)),&
                        xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),&
                        xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
       if (ishape > 0) then
          call delete_shape(ishape,nshapes)
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       elseif (iamincolourbar .and. irenderarr(ipanel) > 0) then
          iColourBarStyle = 0
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       elseif (xpti < xmin(iplotxarr(ipanel)) .or. xpti > xmax(iplotxarr(ipanel)) &
           .or. ypti < xmin(iplotyarr(ipanel)) .or. ypti > xmax(iplotyarr(ipanel))) then
          call deleteaxes()
          istep = istepnew
          interactivereplot = .true.
          iexit = .true.
       else
          ilegend = in_legend(xpt,ypt)
          if (legend_is_plotted(ilegend)) then
             call delete_legend(ilegend)
             iadvance = 0
             interactivereplot = .true.
             iexit = .true.
          else
             print*,' nothing to delete at x,y =',xpti,',',ypti
          endif
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
       read(char,*,iostat=ierr) istepjumpnew
       if (ierr /=0) then
          print*,'*** internal error setting timestep jump'
          istepjumpnew = 1
       endif
       if ((istepjump > 1 .or. istepjumpset) .and. istepjump <= 9999) then
          istepjump = 10*istepjump + istepjumpnew
          if (istepjump > 9999) istepjump = 1
       elseif (istepjumpnew==0) then
          istepjump = 10
       else
          istepjump = istepjumpnew
       endif
       istepjump = int(zoomfac*istepjump)
       print*,' setting timestep jump / zoom factor = ',istepjump
       istepjumpset = .true.
       scalefac = istepjump
    case(')')
       istepjump = int(zoomfac*10)
       print*,' setting timestep jump = ',istepjump
       !
       !--multiply everything by a factor of 10
       !
    case('z','Z')
       zoomfac = 10.*zoomfac
       if (zoomfac > 1000000.) then
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
    call get_vptxy(xpt,ypt,xcursor,ycursor)
    call reset_panel

    !
    !--do not let timestep go outside of bounds
    !  if we are at the first/last step, just print message and do nothing
    !  if iadvance trips over the bounds, jump to last/first step
    !
    if (iadvance /= -666 .and. iexit) then
       if (istep + iadvance  >  ilaststep .and. iframe==nframes) then
          print "(1x,a)",'reached last timestep'
          if (istepin /= ilaststep) then
             istep = ilaststep - istepsonpage*iadvance
          else
             istep = istepin
             iexit = .false.
          endif
       elseif (istep + iadvance  <  1 .and. ifirstframeonpage==1) then
          print "(1x,a)",'reached first timestep: can''t go back'
          if (ifirststeponpage /= 1) then
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

 !--------
 ! utility to return which panel we are in given a point on the viewport
 ! and the viewport limits for each panel.
 !--------
integer function getpanel(vptx,vpty)
 real, intent(in) :: vptx,vpty
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi
 integer :: i,icol

 getpanel = 0
 !
 ! first try the basic procedure of checking if the
 ! cursor falls within the viewport of a particular panel
 !
 do i=1,size(vptxmin)
    if (vptx > vptxmin(i) .and. vptx < vptxmax(i) .and. &
        vpty > vptymin(i) .and. vpty < vptymax(i)) getpanel = i
 enddo
 !
 ! if this fails, use more generous limits around the margins
 !
 icol = 0
 if (getpanel==0) then
    do i=1,size(vptxmin)
       icol = icol + 1
       if (icol > nacross) icol = 1
       if (icol > 1 .and. i > 1) then
          ! if column>1 assign panel by being to the right of previous panel
          vptxmini = vptxmax(i-1)+barwmulti(i-1)
       else
          vptxmini = -0.1 ! allow for some error
       endif
       !--if last column extend xmax to right of page
       if (icol==nacross) then
          vptxmaxi = 1.1
       elseif (verticalbar) then ! otherwise use max of current panel + space containing the colour bar
          vptxmaxi = vptxmax(i) + barwmulti(i)
       else
          vptxmaxi = vptxmax(i)
       endif

       !--if first row extend ymax to top of page
       if (i <= nacross) then
          vptymaxi = 1.1
       else
          vptymaxi = vptymax(i)
       endif
       !--if last row then allow ymin to extend to bottom of page
       if (i > (size(vptxmin)-nacross)) then
          vptymini = -0.1
          ! if not last row assign panel by being above row below
       elseif (i+nacross <= size(vptxmin)) then
          vptymini = vptymax(i+nacross)
       elseif (.not.verticalbar) then
          vptymini = vptymin(i) - barwmulti(i)
       else
          vptymini = vptymin(i)
       endif
       if (vptx > vptxmini .and. vptx < vptxmaxi .and. &
              vpty > vptymini .and. vpty < vptymaxi) then
          !print*,'matching panel ',i,vptx,vpty,vptxmini,vptxmaxi,vptymini,vptymaxi
          !if (getpanel /= 0) print*,'Warning: multiple matching panels found ',getpanel,i,vptx,vpty
          getpanel = i
       endif
    enddo
 endif
 !
 ! if we still fail, just take the last panel
 !
 if (getpanel <= 0 .or. getpanel > size(vptxmin)) then
    !print*,' vptx,y = ',vptx,vpty,vptxmin(:),vptxmax(:)
    !print*,'Error determining panel: assuming last '
    getpanel = size(vptxmin)
 endif

end function getpanel

 !--------
 ! utility to return x,y coordinates in a given panel given viewport coords
 !--------

subroutine getxy(vptx,vpty,x,y,ipanel)
 real, intent(in) :: vptx,vpty
 real, intent(out) :: x,y
 integer, intent(in) :: ipanel

 if (ipanel > 0) then
    x = xmin(iplotxarr(ipanel)) + (vptx-vptxmin(ipanel))/(vptxmax(ipanel)-vptxmin(ipanel)) &
                          *(xmax(iplotxarr(ipanel))-xmin(iplotxarr(ipanel)))
    y = xmin(iplotyarr(ipanel)) + (vpty-vptymin(ipanel))/(vptymax(ipanel)-vptymin(ipanel)) &
                          *(xmax(iplotyarr(ipanel))-xmin(iplotyarr(ipanel)))
 else
    x = 0.
    y = 0.
 endif

 return
end subroutine getxy

 !---
 ! utility which translates between world co-ordinates (x,y)
 ! and viewport co-ordinates (relative to the whole viewport)
 !---
 !subroutine get_vptxy(x,y,vptx,vpty,ipanel)
 ! implicit none
 ! real, intent(in) :: x,y
 ! real, intent(out) :: vptx,vpty
 ! integer, intent(in) :: ipanel
 !
 ! if (ipanel > 0) then
 !    vptx = vptxmin(ipanel) + (x-xmin(iplotxarr(ipanel)))/&
 !    (xmax(iplotxarr(ipanel))-xmin(iplotxarr(ipanel)))*(vptxmax(ipanel)-vptxmini(ipanel))
 !    vpty = vptymin(ipanel) + (y-xmin(iplotyarr(ipanel)))/&(xmax(iplotyarr(ipanel))-xmin(iplotyarr(ipanel)))*(vptymax(ipanel)-vptymini(ipanel))
 ! else
 !   vptx = 0.5
 !    vpty = 0.5
 ! endif
 !
 !end subroutine get_vptxy


 !--------
 ! utility to reset the drawing surface so we can draw in a panel
 !--------

subroutine set_panel(ipanel)
 use plotlib, only:plot_svp,plot_swin
 integer, intent(in) :: ipanel

 if (ipanel > 0) then
    call plot_swin(xmin(iplotxarr(ipanel)),xmax(iplotxarr(ipanel)),xmin(iplotyarr(ipanel)),xmax(iplotyarr(ipanel)))
    !--really should save viewport setting here, but doesn't matter
    !  so long as interactive mode is the last thing called
    call plot_svp(vptxmin(ipanel),vptxmax(ipanel),vptymin(ipanel),vptymax(ipanel))
 endif

 return
end subroutine set_panel

subroutine reset_panel
 use plotlib, only:plot_swin

 call plot_swin(xmini,xmaxi,ymini,ymaxi)

end subroutine reset_panel

end subroutine interactive_multi

!--------------------------------------------------------------------
! utilities to determine whether a point is in or out of a selection
!--------------------------------------------------------------------
logical function inslice(x,xmin,xmax)
 real, intent(in) :: x,xmin,xmax

 inslice = (x >= xmin .and. x <= xmax)

end function inslice

logical function inrectangle(x,y,xmin,xmax,ymin,ymax)
 real, intent(in) :: x,y,xmin,xmax,ymin,ymax

 inrectangle = (x >= xmin .and. x <= xmax .and. y >= ymin .and. y <= ymax)

end function inrectangle

logical function incircle(x,y,r2)
 real, intent(in) :: x,y,r2

 incircle = ((x*x + y*y) <= r2)

end function incircle

!
! Point in polygon
! See: http://en.wikipedia.org/wiki/Even-odd_rule
!
logical function inpoly(x,y,xpts,ypts,npts)
 real, intent(in) :: x,y
 real, dimension(:), intent(in) :: xpts,ypts
 integer, intent(in) :: npts
 integer :: i,j

 inpoly = .false.
 j = npts
 do i=1,npts
    if (((ypts(i) > y) .neqv. (ypts(j) > y)) .and. &
        (x < (xpts(j) - xpts(i))*(y-ypts(i))/(ypts(j) - ypts(i)) + xpts(i))) then
       inpoly = .not. inpoly
    endif
    j = i
 enddo

end function inpoly

!------------------------------------------------------------
! utility which adapts plot limits based only on the
! particles being plotted
!------------------------------------------------------------
subroutine adapt_limits_interactive(labeli,np,xarr,xmin,xmax,icolourpart,iamtype,iusetype)
 use params, only:int1
 use limits, only:assert_sensible_limits
 character(len=*), intent(in)    :: labeli
 integer, intent(in)             :: np
 real, dimension(np), intent(in) :: xarr
 real, intent(out)               :: xmin,xmax
 integer(kind=int1), dimension(:) , intent(in) :: iamtype
 integer,            dimension(np), intent(in) :: icolourpart
 logical,            dimension(:),  intent(in) :: iusetype
 integer :: itype,i
 logical :: mixedtypes

 xmin =  huge(xmin)
 xmax = -huge(xmax)
 mixedtypes = size(iamtype) >= np

 if (mixedtypes) then
    do i=1,np
       itype = int(iamtype(i))
       if (itype > 0 .and. itype <= np) then
          if (iusetype(itype) .and. icolourpart(i) > 0) then
             xmin = min(xmin,xarr(i))
             xmax = max(xmax,xarr(i))
          endif
       endif
    enddo
 else
    xmin = minval(xarr,mask=(icolourpart >= 0))
    xmax = maxval(xarr,mask=(icolourpart >= 0))
 endif
 call assert_sensible_limits(xmin,xmax)

 !print "(1x,a)",' resetting '//trim(labeli)//' limits'

end subroutine adapt_limits_interactive

!------------------------------------------------------------
! utility which translates between world co-ordinates (x,y)
! and viewport co-ordinates (relative to the whole viewport)
!------------------------------------------------------------
subroutine get_vptxy(x,y,vptx,vpty)
 use plotlib, only:plot_qvp,plot_qwin
 real, intent(in) :: x,y
 real, intent(out) :: vptx,vpty
 real :: xmini,xmaxi,ymini,ymaxi
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi

 call plot_qvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
 call plot_qwin(xmini,xmaxi,ymini,ymaxi)
 vptx = vptxmini + (x-xmini)/(xmaxi-xmini)*(vptxmaxi-vptxmini)
 vpty = vptymini + (y-ymini)/(ymaxi-ymini)*(vptymaxi-vptymini)

end subroutine get_vptxy

!------------------------------------------------------------
! utility to return x,y coordinates given viewport coords
! (only works for single-panelled plots)
!------------------------------------------------------------
subroutine get_posxy(vptx,vpty,x,y,xmini,xmaxi,ymini,ymaxi)
 use plotlib, only:plot_qvp
 real, intent(in) :: vptx,vpty
 real, intent(out) :: x,y
 real, intent(in) :: xmini,xmaxi,ymini,ymaxi
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi

 call plot_qvp(0,vptxmini,vptxmaxi,vptymini,vptymaxi)
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
 use plotlib, only:plot_numb,plot_qch,plot_sch,plot_text
 integer, intent(in) :: i
 real, intent(in) :: xi,yi
 integer :: nc
 real :: charheight
 character(len=20) :: string

 !--convert number to text string
 call plot_numb(i,0,1,string,nc)
 !--query and store character height
 call plot_qch(charheight)
 !--change character height
 call plot_sch(2.0)
 !--plot text string
 call plot_text(xi,yi,string(1:nc))
 !--reset character height
 call plot_sch(charheight)

 return
end subroutine plot_number

subroutine deleteaxes()
 use settings_page, only:iaxis,iPlotLegend,& !iPlotStepLegend, &
                    iPlotTitles,iPlotScale
 use settings_vecplot, only:iVecplotLegend

 if (iaxis==-2) then
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
 elseif (iaxis <= 2 .and. iaxis > -2) then
    iaxis = iaxis - 1
 elseif (iaxis > 2) then
    iaxis = -1
 elseif (iaxis < -2) then
    iaxis = -2
 endif

end subroutine deleteaxes

!
!--delete specific legends
!
subroutine delete_legend(id)
 use legends,          only:ilegend,ilegend_vec,ilegend_markers,ilegend_scale
 use settings_page,    only:iPlotLegend,iPlotStepLegend,iPlotScale
 use settings_vecplot, only:iVecplotLegend
 integer, intent(in) :: id

 select case(id)
 case(ilegend_scale)
    if (iPlotScale) print*,'> deleting scale bar'
    iPlotScale = .false.
 case(ilegend_markers)
    if (iPlotStepLegend) print*,'> deleting marker legend'
    iPlotStepLegend = .false.
 case(ilegend_vec)
    if (iVecplotLegend) print*,'> deleting vector legend'
    iVecplotLegend = .false.
 case(ilegend)
    if (iPlotLegend) print*,'> deleting legend'
    iPlotLegend = .false.
 end select

end subroutine delete_legend

!
!--check if legend is actually being plotted
!
logical function legend_is_plotted(id)
 use legends,          only:ilegend,ilegend_vec,ilegend_markers,ilegend_scale
 use settings_page,    only:iPlotLegend,iPlotStepLegend,iPlotScale
 use settings_vecplot, only:iVecplotLegend
 integer, intent(in) :: id

 legend_is_plotted = .false.
 select case(id)
 case(ilegend_scale)
    legend_is_plotted = iPlotScale
 case(ilegend_markers)
    legend_is_plotted = iPlotStepLegend
 case(ilegend_vec)
    legend_is_plotted = iVecplotLegend
 case(ilegend)
    legend_is_plotted = iPlotLegend
 end select

end function legend_is_plotted
!
!--move the legend to the current position
!
subroutine mvlegend(xi,yi,xmin,xmax,ymax,ipanel)
 use settings_page, only:hposlegend,vposlegend,fjustlegend,iPlotLegend,iPlotLegendOnlyOnPanel
 use plotlib, only:plot_qcs
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 integer, intent(in), optional :: ipanel
 real :: xch,ych

 iPlotLegend = .true.
 hposlegend = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call plot_qcs(4,xch,ych)
 vposlegend = (ymax - yi)/ych
! !--automatically change justification
! if (hposlegend < 0.25) then
 fjustlegend = 0.0
! elseif (hposlegend > 0.75) then
!    fjustlegend = 1.0
! else
!    fjustlegend = 0.5
! endif
 if (present(ipanel)) then
    if (ipanel > 0 .and. iPlotLegendOnlyOnPanel > 0) iPlotLegendOnlyOnPanel = ipanel
 endif
 print*,'hpos = ',hposlegend,' vpos = ',vposlegend,' just = ',fjustlegend

end subroutine mvlegend
!
!--move the vector legend to the current position
!
subroutine mvlegendvec(xi,yi,xmin,xmax,ymax)
 use settings_vecplot, only:hposlegendvec,vposlegendvec,iVecplotLegend
 use plotlib, only:plot_qcs
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych

 iVecplotLegend = .true.
 hposlegendvec = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call plot_qcs(4,xch,ych)
 vposlegendvec = (ymax - yi)/ych
 print*,'hpos = ',hposlegendvec,' vpos = ',vposlegendvec

end subroutine mvlegendvec
!
!--move the title to the current position
!
subroutine mvtitle(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hpostitle,vpostitle,fjusttitle,iPlotTitles
 use plotlib, only:plot_qcs
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych

 iPlotTitles = .true.
 hpostitle = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call plot_qcs(4,xch,ych)
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
 use limits,          only:lim,lim2
 use labels,          only:is_coord
 use multiplot,       only:itrans
 use settings_data,   only:ndim
 use settings_limits, only:iadapt,iadaptcoords
 use transforms,      only:transform_limits_inverse
 integer, intent(in) :: iplot
 real,    intent(in) :: xmin,xmax
 logical, intent(in), optional :: setlim2
 logical :: uselim2
 real    :: xmintemp,xmaxtemp

 uselim2 = .false.
 if (present(setlim2)) uselim2 = setlim2

 if (itrans(iplot) /= 0) then
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
 if (is_coord(iplot,ndim)) then
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
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax
 real :: xmintemp,xmaxtemp

 if (itrans(iplot) /= 0) then
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

 call reset_all_ranges()

 return
end subroutine reset_ranges

!
!--interface to routine which resets second set of limits
!
subroutine reset_limits2(icol)
 use limits, only:reset_lim2
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
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax,xi
 real :: xmintemp,xmaxtemp

 if (iplot > ndim) then
    print*,'ERROR in save_limits_track: iplot>ndim'
    return
 elseif (itrans(iplot) /= 0) then
    xmintemp = xmin
    xmaxtemp = xmax
    call transform_limits_inverse(xmintemp,xmaxtemp,itrans(iplot))
    xminoffset_track(iplot) = xi - xmintemp
    xmaxoffset_track(iplot) = xmaxtemp - xi
 else
    xminoffset_track(iplot) = abs(xi - xmin)
    xmaxoffset_track(iplot) = abs(xmax - xi)
 endif

end subroutine save_limits_track
!
!--recalculates radius
!
subroutine save_itrackpart_recalcradius(itrackpart)
 use filenames,      only:nsteps,nstepsinfile,ifileopen
 use settings_data,  only:ncalc,DataIsBuffered,iCalcQuantities,track_string
 use calcquantities, only:calc_quantities,calc_quantities_use_x0
 use part_utils,     only:is_trackstring
 integer, intent(in) :: itrackpart

 if (is_trackstring(track_string)) then
    return  ! do not overwrite strings like "maxdens"
 else
    write(track_string,"(i12)") itrackpart
 endif
 if (iCalcQuantities .and. itrackpart > 0) then
    if (ncalc > 0 .and. calc_quantities_use_x0()) then
       print "(a)",' Recalculating radius relative to tracked particle'
       if (DataIsBuffered) then
          call calc_quantities(1,nsteps)
       else
          call calc_quantities(1,nstepsinfile(ifileopen))
       endif
    endif
 endif

end subroutine save_itrackpart_recalcradius
!
!--toggles log/unlog
!  note this only changes a pure log transform: will not change combinations
!
subroutine change_itrans(iplot,xmin,xmax)
 use multiplot, only:itrans
 use settings_data, only:numplot
 use transforms, only:transform_limits,transform_limits_inverse
 integer, intent(in) :: iplot
 real, intent(inout) :: xmin, xmax

 if (iplot <= numplot) then
    if (itrans(iplot)==1) then
       itrans(iplot) = 0
       !!--untransform the plot limits
       call transform_limits_inverse(xmin,xmax,1)
    else
       itrans(iplot) = 1
       !!--transform the plot limits
       call transform_limits(xmin,xmax,itrans(iplot))
       xmin = max(xmax-4.,xmin) ! no more than 4 dex by default
    endif
 endif

end subroutine change_itrans

subroutine change_itrans2(iplot,xmin,xmax,xmina,xmaxa)
 use multiplot, only:itrans
 use settings_data, only:numplot
 use transforms, only:transform_limits,transform_limits_inverse
 integer, intent(in) :: iplot
 real, intent(inout) :: xmin, xmax, xmina, xmaxa

 if (iplot <= numplot) then
    if (itrans(iplot)==1) then
       itrans(iplot) = 0
       !!--untransform the plot limits
       call transform_limits_inverse(xmin,xmax,1)
       call transform_limits_inverse(xmina,xmaxa,1)
    else
       itrans(iplot) = 1
       !!--transform the plot limits
       call transform_limits(xmin,xmax,itrans(iplot))
       call transform_limits(xmina,xmaxa,itrans(iplot))
       xmin = max(xmax-4.,xmin) ! no more than 4 dex by default
    endif
 endif

end subroutine change_itrans2

!
!--saves rotation options
!
subroutine save_rotation(ndim,anglexi,angleyi,anglezi)
 use settings_xsecrot, only:anglex,angley,anglez
 integer, intent(in) :: ndim
 real, intent(in) :: anglexi,angleyi,anglezi

 anglez = anglezi
 if (ndim >= 3) then
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
 real, intent(in) :: zpos,dz

 zobserver = zpos
 dzscreenfromobserver = dz

 return
end subroutine save_perspective
!
!--saves 3D opacity
!
subroutine save_opacity(rkappai)
 use settings_xsecrot, only:taupartdepth
 real, intent(in) :: rkappai

 taupartdepth = rkappai

 return
end subroutine save_opacity
!
!--save the current paper size
!
subroutine save_windowsize()
 use settings_page, only:ipapersize,ipapersizeunits,papersizex,aspectratio
 use plotlib,       only:plot_qvsz
 real :: x1,x2,y1,y2,papersizey

 call plot_qvsz(0,x1,x2,y1,y2)

 if (abs(x2-x1 - 800.) > 0. .and. abs(y2-y1 - 600.) > 0.) then
    !print*,' saving paper size = ',x2-x1,' x ',y2-y1
    ipapersize = 24  ! custom
    ipapersizeunits = 0
    papersizex = x2-x1
    papersizey = y2-y1
    aspectratio = papersizey/papersizex
 endif

end subroutine save_windowsize

!
!--saves circles of interaction
!
subroutine save_circles(ncircpartset,icircpartset)
 use settings_part, only:ncircpart,icircpart
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
 use colours, only:colour_set,ncolourschemes,icustom
 integer, intent(inout) :: imap
 integer, intent(in) :: istep

 imap = imap + istep
 if (abs(imap) > ncolourschemes .and. abs(imap) /= icustom) imap = 1
 if (abs(imap) < 1) imap = ncolourschemes
 call colour_set(imap)

end subroutine change_colourmap

!
!--set movie mode
!
subroutine set_movie_mode(live)
 use settings_page,   only:iaxis,papersizex,aspectratio,ipapersize,ipapersizeunits,iPageColours
 use settings_limits, only:adjustlimitstodevice
 use settings_render, only:iColourBarStyle
 use pagecolours,     only:set_pagecolours
 use plotlib,         only:plotlib_is_pgplot,plot_pap
 use colourbar,       only:set_floating_bar_style
 use system_utils,    only:get_copyright
 use shapes,          only:add_text
 logical, intent(in) :: live

 iaxis = -2
 iPageColours = 2
 if (.not.plotlib_is_pgplot) then
    ipapersize      = 9
    ipapersizeunits = 0
    papersizex      = 1280.
    aspectratio     = 0.5625
    if (live) call plot_pap(papersizex,aspectratio,ipapersizeunits)
    iColourBarStyle = 8
    call set_floating_bar_style(iColourBarStyle,4)
    if (live) call set_pagecolours(iPageColours)
    adjustlimitstodevice = .true.
 endif
 call add_text(0.025,0.05,get_copyright())

end subroutine set_movie_mode

!
!--unset movie mode
!
subroutine unset_movie_mode()
 use settings_page,   only:iaxis,papersizex,aspectratio,ipapersize,iPageColours
 use settings_limits, only:adjustlimitstodevice
 use settings_render, only:iColourBarStyle
 use pagecolours,     only:set_pagecolours
 use plotlib,         only:plotlib_is_pgplot,plot_pap
 use system_utils,    only:get_copyright
 use shapes,          only:delete_text

 iaxis = 0
 iPageColours = 1
 if (.not.plotlib_is_pgplot) then
    ipapersize   = 0
    papersizex   = 800.
    aspectratio  = 600./800.
    iColourBarStyle   = 1
    call plot_pap(papersizex,aspectratio,0)
    call set_pagecolours(iPageColours)
    adjustlimitstodevice = .true.
 endif
 call delete_text(get_copyright())

end subroutine unset_movie_mode

end module interactive_routines
