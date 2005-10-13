!------------------------------------------------------------------------
!  Module containing "interface" routines between the calculated
!  pixel arrays and the PGPLOT routines which do the actual rendering
!------------------------------------------------------------------------
module render
 implicit none
 public :: render_pix, render_vec, colourbar
 private

contains

!------------------------------------------------------------------------
!  this subroutine takes a 2D grid of data and renders it using pgplot
!  rendering is either greyscale (icolours = 1) or colour (icolours>1)
!  also plots nc contours between datmin and datmax.
!------------------------------------------------------------------------
 
subroutine render_pix(datpix,datmin,datmax,label,npixx,npixy, &
                  xmin,ymin,dx,icolours,iplotcont,iPlotColourBar,nc,log)
 implicit none
 integer, intent(in) :: npixx,npixy,nc,icolours
 real, intent(in) :: xmin,ymin,datmin,datmax,dx
 real, dimension(npixx,npixy), intent(in) :: datpix
 logical, intent(in) :: iplotcont,iPlotColourBar,log
 character(len=*), intent(in) :: label  
 
 integer :: i
 real :: trans(6),levels(nc),dcont
! 
!--set up grid for rendering 
!
 trans(1) = xmin - 0.5*dx                ! this is for the pgimag call
 trans(2) = dx                        ! see help for pgimag/pggray/pgcont
 trans(3) = 0.0
 trans(4) = ymin - 0.5*dx
 trans(5) = 0.0
 trans(6) = dx

 print*,'rendering...',npixx,'x',npixy,',array size=',size(datpix),minval(datpix)

 if (icolours.eq.1) then        ! greyscale
    if (iPlotColourBar) call colourbar(icolours,datmin,datmax,trim(label),log)
    call pggray(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)

 elseif (icolours.gt.1) then        ! colour
    if (iPlotColourBar) call colourbar(icolours,datmin,datmax,trim(label),log)
!    call pgwedg('ri',2.0,4.0,datmin,datmax,' ')
!    call pgpixl(datpix,npixx,npixx,1,npixx,1,npixx,xmin,xmax,ymin,ymax)
    call pgimag(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
!    call pghi2d(datpix,npixx,npixx,1,npixx,1,npixx,1,0.1,.true.,y) 
 endif
!
!--contours
!
 if (iplotcont) then
    print*,'plotting ',nc,' contours...'
!
!--set contour levels
! 
    dcont = (datmax-datmin)/real(nc+1)   ! even contour levels
    do i=1,nc
       levels(i) = datmin + real(i)*dcont
    enddo
!
!--plot contours (use pgcont if pgcons causes trouble)
!
    call pgcons(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans)
    call pgmtxt('T',-2.0,0.05,0.0,trim(label))

 endif
 
 return
 
end subroutine render_pix

!-------------------------------------------------------
! this subroutine interfaces to my version of PGWEDG
! which plots the colour bar (differences are that
! text is written vertically, txtsep is a
! changeable parameter and the character height
! is not changed)
!-------------------------------------------------------
subroutine colourbar(icolours,datmin,datmax,label,log)
 use settings_render, only:ColourBarDisp, ColourBarWidth
 implicit none
 integer, intent(in) :: icolours
 real, intent(in) :: datmin,datmax
 character(len=*), intent(in) :: label
 logical, intent(in) :: log
 character(len=1) :: clog
 real :: disp, width
!
!--set colour bar displacement and width in character heights
!
 disp = 0.5
 width = ColourBarWidth
!
!--set character to send to pgwedg call if log (danpgwedg only) 
!
 clog = ' '
 if (log) clog = 'l'
!
!--Note that plots use my modification of pgwedg which plots vertical numbers on axes
!          
 if (icolours.eq.1) then        ! greyscale
    call danpgwedg('rgv'//clog,disp,width,datmin,datmax,trim(label),ColourBarDisp)
 elseif (icolours.gt.1) then        ! colour
    call danpgwedg('riv'//clog,disp,width,datmin,datmax,trim(label),ColourBarDisp)
 endif

 return
end subroutine colourbar

!--------------------------------------------------------------------------
!  this subroutine takes a 2D grid of vector data (ie. x and y components)
!  and plots an arrow map of it
!--------------------------------------------------------------------------
 
subroutine render_vec(vecpixx,vecpixy,vecmax,npixx,npixy,        &
                  xmin,ymin,dx,label) 
 use legends, only:legend_vec
 use settings_vecplot, only:iVecplotLegend,hposlegendvec,vposlegendvec
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in) :: xmin,ymin,dx
 real, intent(inout) :: vecmax
 real, dimension(npixx,npixy), intent(in) :: vecpixx,vecpixy
 character(len=*), intent(in) :: label
 real :: trans(6),scale
 real :: charheight
 
!set up grid for rendering 

 trans(1) = xmin !- 0.5*dx                ! this is for the pgimag call
 trans(2) = dx                        ! see help for pgimag/pggray/pgcont
 trans(3) = 0.0
 trans(4) = ymin !- 0.5*dx
 trans(5) = 0.0
 trans(6) = dx

 print*,trim(label),' vector plot..',npixx,'x',npixy,',array size=',size(vecpixx)
 !!print*,'max(x component) = ',maxval(vecpixx),'max(y component) = ',maxval(vecpixy)

 call pgsah(2,45.0,0.7)   ! arrow style
 call pgqch(charheight)
 call pgsch(0.3)          ! size of arrow head
 if (vecmax.le.0.0) then  ! adaptive limits
    scale = 0.0
    vecmax = max(maxval(vecpixx(:,:)),maxval(vecpixy(:,:)))
    if (vecmax.gt.0.) scale = dx/vecmax
 else
    scale=dx/vecmax
 endif
 
 call pgvect(vecpixx(:,:),vecpixy(:,:),npixx,npixy, &
      1,npixx,1,npixy,scale,0,trans,-1000.0)

 if (iVecplotLegend) then
    call legend_vec(label,vecmax,dx,hposlegendvec,vposlegendvec,charheight)
 endif
 call pgsch(charheight)
 
 return
 
end subroutine render_vec

end module render
