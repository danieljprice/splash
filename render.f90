!------------------------------------------------------------------------
!  this subroutine takes a 2D grid of data and renders it using pgplot
!  rendering is either greyscale (icolours = 1) or colour (icolours>1)
!  also plots nc contours between datmin and datmax.
!------------------------------------------------------------------------
 
subroutine render(datpix,datmin,datmax,label,npixx,npixy,	&
	          xmin,ymin,dx,icolours,iplotcont,nc,log)  
 implicit none
 integer, intent(in) :: npixx,npixy,nc,icolours
 real, intent(in) :: xmin,ymin,datmin,datmax,dx
 real, dimension(npixx,npixy), intent(in) :: datpix
 logical, intent(in) :: iplotcont,log
 character(len=*), intent(in) :: label  
 
 integer i,j,k
 real :: trans(6),levels(nc),dcont
 character(len=1) :: clog
 
!set up grid for rendering 

 trans(1) = xmin - 0.5*dx		! this is for the pgimag call
 trans(2) = dx			! see help for pgimag/pggray/pgcont
 trans(3) = 0.0
 trans(4) = ymin - 0.5*dx
 trans(5) = 0.0
 trans(6) = dx

!set character to send to pgwedg call if log (danpgwedg only) 
 clog = ' '
 if (log) clog = 'l'

 print*,'rendering...',npixx,'x',npixy,',array size=',size(datpix),minval(datpix)
!
!--set contour levels
! 
 dcont = (datmax-datmin)/real(nc+1)   ! even contour levels
 do i=1,nc
    levels(i) = datmin + real(i)*dcont
 enddo
!
!--nb: plots use my modification of pgwedg which plots vertical numbers on axes
!	 
 if (icolours.eq.1) then	! greyscale
    call danpgwedg('rgv'//clog,0.5,4.5,datmin,datmax,label)
    call pggray(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)

 elseif (icolours.gt.1) then	! colour
    call danpgwedg('riv'//clog,0.5,4.5,datmin,datmax,label)
!    call pgwedg('ri',2.0,4.0,datmin,datmax,' ')
!    call pgpixl(datpix,npixx,npixx,1,npixx,1,npixx,xmin,xmax,ymin,ymax)
    call pgimag(datpix,npixx,npixy,1,npixx,1,npixy,datmin,datmax,trans)
!    call pghi2d(datpix,npixx,npixx,1,npixx,1,npixx,1,0.1,.true.,y)
	 
 endif
!
!--plot contours
!
 if (iplotcont) then
    print*,'plotting ',nc,' contours...'
!--use pgcont if pgcons causes trouble
    call pgcons(datpix,npixx,npixy,1,npixx,1,npixy,levels,nc,trans)
 endif
 
 return
 
end subroutine render
