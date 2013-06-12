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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!------------------------------------------------------------------------
!  Module containing "interface" routines between the calculated
!  pixel arrays and the plot library routines which do the actual rendering
!------------------------------------------------------------------------
module render
 use colourbar, only:plotcolourbar
 implicit none
 public :: render_pix, render_vec
 private

contains

!------------------------------------------------------------------------
!  this subroutine takes a 2D grid of data and renders it via the
!  plotting library. Rendering is either:
!  - contours (icolouropt = 0)
!  - greyscale (icolouropt = 1)
!  - colour (icolouropt>1)
!  contouring plots nc contours between datmin and datmax.
!------------------------------------------------------------------------
subroutine render_pix(datpix,datmin,datmax,label,npixx,npixy, &
                  xmin,ymin,dx,dy,icolouropt,iplotcont,iColourBarStyle,ncontours,log, &
                  ilabelcont,contmin,contmax,blank,transparent,alpha)
 use plotutils,       only:formatreal
 use plotlib,         only:plot_imag,plot_conb,plot_cons,plot_qch,plot_sch,&
                           plot_qch,plot_sch,plot_conl,plot_gray,plot_imag_transparent,&
                           plot_imag_alpha
 use contours_module, only:read_contours,contours_list,contourtitles
 implicit none
 integer, intent(in) :: npixx,npixy,ncontours,icolouropt
 real, intent(in) :: xmin,ymin,datmin,datmax,dx,dy
 real, dimension(npixx,npixy), intent(in) :: datpix
 logical, intent(in) :: iplotcont,log,ilabelcont
 integer, intent(in) :: iColourBarStyle
 character(len=*), intent(in) :: label
 real, intent(in), optional :: contmin,contmax,blank
 logical, intent(in), optional :: transparent
 real, dimension(npixx,npixy), intent(in), optional :: alpha

 integer :: i,ierr,nc
 real :: trans(6),levels(ncontours),dcont,charheight,cmin,cmax
 character(len=12) :: string
 logical :: iuse_transparent,ifixed_contours
!
!--set up grid for rendering
!
 trans(1) = xmin - 0.5*dx      ! this is for the pgimag call
 trans(2) = dx                 ! see help for pgimag/pggray/pgcont
 trans(3) = 0.0
 trans(4) = ymin - 0.5*dy
 trans(5) = 0.0
 trans(6) = dy

 iuse_transparent = .false.
 if (present(transparent)) iuse_transparent = transparent

 print*,'rendering...',npixx,'x',npixy,'=',size(datpix),' pixels'

 if (abs(icolouropt).eq.1) then        ! greyscale

    if (iColourBarStyle.gt.0) call plotcolourbar(iColourBarstyle,icolouropt,datmin,datmax,trim(label),log,0.)

    if (icolouropt.eq.1) then
       call plot_gray(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
    else !--allow inverse greyscale
       call plot_imag(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
    endif
 elseif (abs(icolouropt).gt.0) then        ! colour
    !
    !--plot colour bar
    !
    if (iColourBarStyle.gt.0) call plotcolourbar(iColourBarstyle,icolouropt,datmin,datmax,trim(label),log,0.)
    !
    !--plot pixel map
    !
    if (iuse_transparent) then
       call plot_imag_transparent(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
    else
       if (present(alpha)) then
          call plot_imag_alpha(datpix,alpha,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
       else
          call plot_imag(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
       endif
   endif
 endif
!
!--contours
!
 if (iplotcont) then
    nc = ncontours
    if (present(contmin)) then
       cmin = contmin
    else
       cmin = datmin
    endif
    if (present(contmax)) then
       cmax = contmax
    else
       cmax = datmax
    endif
!
!--set contour levels: first attempt to read these
!  from a file. If file does not exist or errors during read
!  then we construct the default levels as usual.
!
    call read_contours(nc,ierr)
    if (ierr.eq.0 .and. nc.gt.0) then
       ifixed_contours = .true.
    else
       nc = ncontours
       ifixed_contours = .false.
    endif

    if (ifixed_contours) then
       do i=1,min(nc,ncontours)
          print*,"contour @ ", contours_list(i), ": ", trim(contourtitles(i))
          levels(i) = contours_list(i)
       enddo
       dcont = 0.
    elseif (nc.le.0) then
       print*,'ERROR: cannot plot contours with ',nc,' levels'
       return
    elseif (nc.eq.1) then
       levels(1) = cmin
       dcont = 0.
    else
       dcont = (cmax-cmin)/real(nc-1)   ! even contour levels
       do i=1,nc
          levels(i) = cmin + real(i-1)*dcont
       enddo
    endif
!
!--plot contours (use pgcont if pgcons causes trouble)
!  with blanking if blank is input
!
    if (present(blank)) then
       if (.not.ifixed_contours) then
          print 10,nc,' contours (with blanking)',levels(1),levels(nc),dcont
          print 20,levels(1:nc)
       endif
       !print*,' blanking = ',blank,'min,max = ',datmin,datmax
       call plot_conb(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans,blank)
    else
       if (.not.ifixed_contours) then
          print 10,nc,' contours',levels(1),levels(nc),dcont
          print 20,levels(1:nc)
       endif
       call plot_cons(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans)
    endif
10  format(1x,'plotting ',i4,a,' between ',es10.2,' and ',es10.2,', every ',es10.2,':')
20  format(10(6(1x,es9.2),/))
!
!--labelling of contour levels
!
    if (ilabelcont) then
       call plot_qch(charheight)       ! query character height
       call plot_sch(0.75*charheight)   ! shrink character height

       do i=1,nc
          if (ifixed_contours) then
             string=adjustl(contourtitles(i))
          else
             call formatreal(levels(i),string)
          endif
          call plot_conl(datpix,npixx,npixy,1,npixx,1,npixy,levels(i),trans,trim(string),npixx/2,30)
       enddo
       call plot_sch(charheight) ! restore character height
    endif
!
!--this line prints the label inside the contour plot
!  (now obsolete-- this functionality can be achieved using plot titles)
!    call pgmtxt('T',-2.0,0.05,0.0,trim(label))

 endif

 return

end subroutine render_pix

!--------------------------------------------------------------------------
!  this subroutine takes a 2D grid of vector data (ie. x and y components)
!  and plots an arrow map of it
!--------------------------------------------------------------------------

subroutine render_vec(vecpixx,vecpixy,vecmax,npixx,npixy, &
                      xmin,ymin,dx,dy,label,unitslabel,plotlegend)
 use legends,          only:legend_vec
 use settings_vecplot, only:hposlegendvec,vposlegendvec,&
                            iplotarrowheads,iallarrowssamelength
 use plotlib,          only:plot_sah,plot_qch,plot_sch,plot_vect
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in) :: xmin,ymin,dx,dy
 real, intent(inout) :: vecmax
 real, dimension(npixx,npixy), intent(in) :: vecpixx,vecpixy
 real, dimension(npixx,npixy) :: dvmag
 character(len=*), intent(in) :: label,unitslabel
 logical,          intent(in) :: plotlegend
 real :: trans(6),scale
 real :: charheight

!set up grid for rendering

 trans(1) = xmin - 0.5*dx                ! this is for the pgimag call
 trans(2) = dx                        ! see help for pgimag/pggray/pgcont
 trans(3) = 0.0
 trans(4) = ymin - 0.5*dy
 trans(5) = 0.0
 trans(6) = dy

 print*,'vector plot..',npixx,'x',npixy,'=',size(vecpixx),' pixels'
 !!print*,'max(x component) = ',maxval(vecpixx),'max(y component) = ',maxval(vecpixy)

 if (iplotarrowheads) then
    call plot_sah(2,45.0,0.7)   ! arrow style
 else
    call plot_sah(2,0.0,1.0)
 endif
 call plot_qch(charheight)
 call plot_sch(0.3)          ! size of arrow head

 if (iallarrowssamelength) then
    !!if (vecmax.le.0.0) vecmax = 1.0 ! adaptive limits
    scale=0.9*dx !!/vecmax
    print*,trim(label),' showing direction only: max = ',vecmax

    where (abs(vecpixx).gt.tiny(vecpixx) .and. abs(vecpixy).gt.tiny(vecpixy))
       dvmag(:,:) = 1./sqrt(vecpixx**2 + vecpixy**2)
    elsewhere
       dvmag(:,:) = 0.
    end where

    call plot_vect(vecpixx(:,:)*dvmag(:,:),vecpixy(:,:)*dvmag(:,:),npixx,npixy, &
         1,npixx,1,npixy,scale,0,trans,0.0)
 else
    if (vecmax.le.0.0) then  ! adaptive limits
       scale = 0.0
       vecmax = max(maxval(vecpixx(:,:)),maxval(vecpixy(:,:)))
       if (vecmax.gt.0.) scale = dx/vecmax
    else
       scale=dx/vecmax
    endif
    print*,trim(label),' max = ',vecmax

    call plot_vect(vecpixx(:,:),vecpixy(:,:),npixx,npixy, &
         1,npixx,1,npixy,scale,0,trans,0.0)

    if (plotlegend) then
       call legend_vec(label,unitslabel,vecmax,dx,hposlegendvec,vposlegendvec,charheight)
    endif
 endif

 call plot_sch(charheight)

 return

end subroutine render_vec

end module render
