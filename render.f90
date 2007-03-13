!------------------------------------------------------------------------
!  Module containing "interface" routines between the calculated
!  pixel arrays and the PGPLOT routines which do the actual rendering
!------------------------------------------------------------------------
module render
 implicit none
 public :: render_pix, render_vec, render_opacity, colourbar
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

 print*,'rendering...',npixx,'x',npixy,'=',size(datpix),' pixels'

 if (abs(icolours).eq.1) then        ! greyscale
    if (iPlotColourBar) call colourbar(icolours,datmin,datmax,trim(label),log)
    call pggray(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)

 elseif (abs(icolours).gt.1) then        ! colour
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
!--this line prints the label inside the contour plot
!  (now obsolete-- this functionality can be achieved using plot titles)
!    call pgmtxt('T',-2.0,0.05,0.0,trim(label))

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
subroutine colourbar(icolours,datmin,datmax,label,log, &
                     vptxmaxfull,vptyminfull,vptymaxfull)
 use settings_render, only:ColourBarDisp, ColourBarWidth
 implicit none
 integer, intent(in) :: icolours
 real, intent(in) :: datmin,datmax
 character(len=*), intent(in) :: label
 logical, intent(in) :: log
 real, intent(in), optional :: vptxmaxfull,vptyminfull,vptymaxfull
 integer, parameter :: npixwedg = 400
 real, dimension(6), parameter :: trans = (/0.0,1.0,0.0,0.0,0.0,1.0/)
 real, dimension(npixwedg) :: sample
 integer :: i
 character(len=1) :: clog
 real :: disp,width,xch,ych,dx
 real :: xmin,xmax,ymin,ymax,vptxmin,vptxmax,vptymin,vptymax
 real :: vptxmini,vptxmaxi,vptymini,vptymaxi
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

 call pgbbuf
 call pgqwin(xmin,xmax,ymin,ymax)
 call pgqvp(0,vptxmin,vptxmax,vptymin,vptymax)
 call pgqcs(0,xch,ych)
 !--if colour bar stretches across multiple plots,
 !  override settings for vptymin and vptymax with input values
 if (present(vptxmaxfull) .and. present(vptyminfull) .and. present(vptymaxfull)) then
    vptxmaxi = vptxmaxfull
    vptymini = vptyminfull
    vptymaxi = vptymaxfull
 else
    vptxmaxi = vptxmax
    vptymini = vptymin
    vptymaxi = vptymax
 endif
!
!--translate width and displacement to viewport co-ordinates
!
 width = width*xch
 disp = disp*xch
!
!--set viewport for the wedge
! 
 vptxmini = vptxmaxi + disp
 vptxmaxi = vptxmini + width*0.4
 call pgsvp(vptxmini,vptxmaxi,vptymini,vptymaxi)
!
!--fill array with all values from datmin to datmax
!
 dx = (datmax-datmin)/real(npixwedg-1)
 do i=1,npixwedg
    sample(i) = datmin + (i-1)*dx
 enddo
!
!--draw colour bar, by cleverly setting window size
!
 call pgswin(0.9,1.1,1.0,real(npixwedg))
 if (abs(icolours).eq.1) then        ! greyscale
    call pggray(sample,1,npixwedg,1,1,1,npixwedg,datmin,datmax,trans)
 elseif (abs(icolours).gt.1) then        ! colour
    call pgimag(sample,1,npixwedg,1,1,1,npixwedg,datmin,datmax,trans)
 endif
 call pgswin(0.0,1.0,datmin,datmax)
!
!--draw labelled frame around the wedge
!
 call pgbox('BC',0.0,0,'BCMSTVRV',0.0,0)
!
!--write the units label
!
 if (label.ne.' ') then
    call pgmtxt('R',ColourBarDisp+1.0,1.0,1.0,trim(label))
 endif
!
!--reset window and viewport
!
 call pgsvp(vptxmin,vptxmax,vptymin,vptymax)
 call pgswin(xmin,xmax,ymin,ymax)
 call pgebuf

 return
end subroutine colourbar

!--------------------------------------------------------------------------
!  this subroutine takes a 2D grid of vector data (ie. x and y components)
!  and plots an arrow map of it
!--------------------------------------------------------------------------
 
subroutine render_vec(vecpixx,vecpixy,vecmax,npixx,npixy,        &
                  xmin,ymin,dx,label,unitslabel) 
 use legends, only:legend_vec
 use settings_vecplot, only:iVecplotLegend,hposlegendvec,vposlegendvec,iplotarrowheads
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in) :: xmin,ymin,dx
 real, intent(inout) :: vecmax
 real, dimension(npixx,npixy), intent(in) :: vecpixx,vecpixy
 character(len=*), intent(in) :: label,unitslabel
 real :: trans(6),scale
 real :: charheight
 
!set up grid for rendering 

 trans(1) = xmin - 0.5*dx                ! this is for the pgimag call
 trans(2) = dx                        ! see help for pgimag/pggray/pgcont
 trans(3) = 0.0
 trans(4) = ymin - 0.5*dx
 trans(5) = 0.0
 trans(6) = dx

 print*,'vector plot..',npixx,'x',npixy,'=',size(vecpixx),' pixels'
 !!print*,'max(x component) = ',maxval(vecpixx),'max(y component) = ',maxval(vecpixy)

 if (iplotarrowheads) then
    call pgsah(2,45.0,0.7)   ! arrow style
 else
    call pgsah(2,0.0,1.0)
 endif
 call pgqch(charheight)
 call pgsch(0.3)          ! size of arrow head
 if (vecmax.le.0.0) then  ! adaptive limits
    scale = 0.0
    vecmax = max(maxval(vecpixx(:,:)),maxval(vecpixy(:,:)))
    if (vecmax.gt.0.) scale = dx/vecmax
 else
    scale=dx/vecmax
 endif
 print*,trim(label),' max = ',vecmax
 
 call pgvect(vecpixx(:,:),vecpixy(:,:),npixx,npixy, &
      1,npixx,1,npixy,scale,0,trans,0.0)

 if (iVecplotLegend) then
    call legend_vec(label,unitslabel,vecmax,dx,hposlegendvec,vposlegendvec,charheight)
 endif
 call pgsch(charheight)
 
 return
 
end subroutine render_vec

!
!--attempt to render an array of red, green and blue colours
!  using PGPLOT
!
subroutine render_opacity(rgbcolours,npixx,npixy,xmin,xmax,ymin,ymax, &
                          iPlotColourBar,icolours,datmin,datmax,label)
 implicit none
 integer, intent(in) :: npixx,npixy
 real, dimension(3,npixx,npixy), intent(in) :: rgbcolours
 real, intent(in) :: xmin,xmax,ymin,ymax,datmin,datmax
 logical, intent(in) :: iPlotColourBar
 integer, intent(in) :: icolours
 character(len=*), intent(in) :: label
 
 integer, dimension(npixx,npixy) :: icolourarray
 integer :: ncolours,ir,ig,ib,nshades,nshades2,index,ipix,jpix
 integer :: indexmax,indexmin
 real :: denom,red,green,blue

 if (iPlotColourBar) call colourbar(icolours,datmin,datmax,trim(label),.false.) 
!
!--set the colour table corresponding to all possible combinations
!  of red, green and blue
!
 call pgqcol(indexmin,indexmax)
 ncolours = indexmax - indexmin + 1
 nshades = int(ncolours**(1./3.))
 print*,'ncolours = ',ncolours,'nshades = ',nshades,nshades**3
 denom = 1./real(nshades-1)
 nshades2 = nshades*nshades
 !ncolours = nshades*nshades*nshades

 do ir=1,nshades
    red = (ir-1)*denom
    do ig=1,nshades
       green = (ig-1)*denom
       do ib=1,nshades
          index = (ir-1)*nshades2 + (ig-1)*nshades + (ib-1) + indexmin
          blue = (ib-1)*denom
          call pgscr(index,red,green,blue)
       enddo
    enddo
 enddo

 do jpix=1,npixy
    do ipix=1,npixx
       ir = int(rgbcolours(1,ipix,jpix)*nshades)
       ig = int(rgbcolours(2,ipix,jpix)*nshades)
       ib = int(rgbcolours(3,ipix,jpix)*nshades)
       icolourarray(ipix,jpix) = (ir-1)*nshades2 + (ig-1)*nshades + (ib-1) + indexmin
    enddo
 enddo
 
 call pgpixl(icolourarray,npixx,npixy,1,npixx,1,npixy,xmin,xmax,ymin,ymax)

end subroutine render_opacity

end module render
