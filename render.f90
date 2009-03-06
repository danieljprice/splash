!------------------------------------------------------------------------
!  Module containing "interface" routines between the calculated
!  pixel arrays and the PGPLOT routines which do the actual rendering
!------------------------------------------------------------------------
module render
 use colourbar, only:plotcolourbar
 implicit none
 public :: render_pix, render_vec, render_opacity
 private

contains

!------------------------------------------------------------------------
!  this subroutine takes a 2D grid of data and renders it using pgplot
!  rendering is either greyscale (icolours = 1) or colour (icolours>1)
!  also plots nc contours between datmin and datmax.
!------------------------------------------------------------------------
 
subroutine render_pix(datpix,datmin,datmax,label,npixx,npixy, &
                  xmin,ymin,dx,icolours,iplotcont,iColourBarStyle,nc,log,ilabelcont,blank)
 use plotutils, only:formatreal
 implicit none
 integer, intent(in) :: npixx,npixy,nc,icolours
 real, intent(in) :: xmin,ymin,datmin,datmax,dx
 real, dimension(npixx,npixy), intent(in) :: datpix
 logical, intent(in) :: iplotcont,log,ilabelcont
 integer, intent(in) :: iColourBarStyle
 character(len=*), intent(in) :: label
 real, intent(in), optional :: blank
 
 integer :: i,ierr
 real :: trans(6),levels(nc),dcont,charheight
 character(len=12) :: string
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

! if (abs(icolours).eq.1) then        ! greyscale
!    if (iPlotColourBar) call colourbar(icolours,datmin,datmax,trim(label),log)
!    call pggray(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)

 if (abs(icolours).gt.0) then        ! colour
    if (iColourBarStyle.gt.0) call plotcolourbar(iColourBarstyle,icolours,datmin,datmax,trim(label),log,0.)
!    call pgwedg('ri',2.0,4.0,datmin,datmax,' ')
    call pgimag(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
!    call pghi2d(datpix,npixx,npixx,1,npixx,1,npixx,1,0.1,.true.,y) 
 endif
!
!--contours
!
 if (iplotcont) then
!
!--set contour levels
! 
    if (nc.le.0) then
       print*,'ERROR: cannot plot contours with ',nc,' levels'
       return
    elseif (nc.eq.1) then
       levels(1) = datmin
       dcont = 0.
    else
       dcont = (datmax-datmin)/real(nc-1)   ! even contour levels
       do i=1,nc
          levels(i) = datmin + real(i-1)*dcont
       enddo
    endif
!
!--plot contours (use pgcont if pgcons causes trouble)
!  with blanking if blank is input
!
    if (present(blank)) then
       print 10,nc,' contours (with blanking)',levels(1),levels(nc),dcont
       print 20,levels(1:nc)
       !print*,' blanking = ',blank,'min,max = ',datmin,datmax
       call pgconb(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans,blank)
    else
       print 10,nc,' contours',levels(1),levels(nc),dcont
       print 20,levels(1:nc)
       call pgcons(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans)
    endif
10  format(1x,'plotting ',i4,a,' between ',1pe8.2,' and ',1pe8.2,', every ',1pe8.2,':')
20  format(10(6(1x,1pe9.2),/))
!
!--labelling of contour levels
!
    if (ilabelcont) then
       call pgqch(charheight)       ! query character height
       call pgsch(0.75*charheight)   ! shrink character height

       do i=1,nc
          call formatreal(levels(i),string)
          call pgconl(datpix,npixx,npixy,1,npixx,1,npixy,levels(i),trans,trim(string),npixx/2,30)
       enddo
       call pgsch(charheight) ! restore character height
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
 
subroutine render_vec(vecpixx,vecpixy,vecmax,npixx,npixy,        &
                  xmin,ymin,dx,label,unitslabel) 
 use legends, only:legend_vec
 use settings_vecplot, only:iVecplotLegend,hposlegendvec,vposlegendvec,iplotarrowheads,&
                            iallarrowssamelength
 implicit none
 integer, intent(in) :: npixx,npixy
 real, intent(in) :: xmin,ymin,dx
 real, intent(inout) :: vecmax
 real, dimension(npixx,npixy), intent(in) :: vecpixx,vecpixy
 real, dimension(npixx,npixy) :: dvmag
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
 
 if (iallarrowssamelength) then
    !!if (vecmax.le.0.0) vecmax = 1.0 ! adaptive limits
    scale=0.9*dx !!/vecmax
    print*,trim(label),' showing direction only: max = ',vecmax

    where (abs(vecpixx).gt.tiny(vecpixx) .and. abs(vecpixy).gt.tiny(vecpixy)) 
       dvmag(:,:) = 1./sqrt(vecpixx**2 + vecpixy**2)
    elsewhere
       dvmag(:,:) = 0.
    end where
    
    call pgvect(vecpixx(:,:)*dvmag(:,:),vecpixy(:,:)*dvmag(:,:),npixx,npixy, &
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

    call pgvect(vecpixx(:,:),vecpixy(:,:),npixx,npixy, &
         1,npixx,1,npixy,scale,0,trans,0.0)

    if (iVecplotLegend) then
       call legend_vec(label,unitslabel,vecmax,dx,hposlegendvec,vposlegendvec,charheight)
    endif
 endif
 
 call pgsch(charheight)
 
 return
 
end subroutine render_vec

!
!--attempt to render an array of red, green and blue colours
!  using PGPLOT
!
subroutine render_opacity(rgbcolours,npixx,npixy,xmin,xmax,ymin,ymax, &
                          iColourBarStyle,icolours,datmin,datmax,label)
 implicit none
 integer, intent(in) :: npixx,npixy
 real, dimension(3,npixx,npixy), intent(in) :: rgbcolours
 real, intent(in) :: xmin,xmax,ymin,ymax,datmin,datmax
 integer, intent(in) :: iColourBarStyle
 integer, intent(in) :: icolours
 character(len=*), intent(in) :: label
 
 integer, dimension(npixx,npixy) :: icolourarray
 integer :: ncolours,ir,ig,ib,nshades,nshades2,index,ipix,jpix
 integer :: indexmax,indexmin
 real :: denom,red,green,blue

 if (iColourBarStyle.gt.0) call plotcolourbar(iColourBarStyle,icolours,datmin,datmax,trim(label),.false.,0.) 
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
