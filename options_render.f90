!-----------------------------------------------------------------------------
! options for rendering / vector plots
!-----------------------------------------------------------------------------

subroutine options_render(npix,icolours,iplotcont, 		&
 			  ncontours,ivecplot,npixvec,iplotpartvec,	&
			  x_sec,xsecpos,backgnd_vec,ndim,numplot)
 use prompting
 implicit none
 integer, intent(in) :: ndim, numplot
 integer, intent(inout) :: npix,icolours,ivecplot,npixvec
 integer, intent(inout) :: ncontours
 logical, intent(inout) :: iplotcont,iplotpartvec,x_sec,backgnd_vec
 real, intent(inout) :: xsecpos
 character(len=1) :: ans 
!
!--rendering options
!
 call prompt('enter number of pixels along x axis',npix,1,1000)
 
100 continue
 print*,'(-ve = demo, 0 = contours only)'
 call prompt('enter colour scheme for rendering ',icolours,max=5)
 if (icolours.lt.0) then
    call colour_demo
    goto 100
 endif   
    
 if (ndim.eq.3) then
    print*,' take cross section (c) or projection (p) through 3d data?'
    read*,ans
    x_sec = .false.
    if (ans.eq.'c') x_sec = .true.
    print*,' cross section = ',x_sec
    xsecpos = 0.0
    if (x_sec) then
       !call prompt('Do you want a fly-through',flythru)
       !if (.not.flythru) then
       call prompt('enter co-ordinate location of cross section slice',xsecpos)
       !endif
    endif   
 endif
	 
 if (icolours.ne.0) then
    call prompt(' plot contours?',iplotcont)
 else
    iplotcont = .true.
 endif
 
 if (iplotcont) then
    call prompt(' enter number of contours between min,max',ncontours,1,100)
 endif
!
!--vector plot options
!	  
 call prompt('vector plot: velocity (1), magnetic field (2) or none(0)?',ivecplot,0,2)
 if (ivecplot.gt.0) then	
    call prompt('enter number of pixels',npixvec,1,1000)
    iplotpartvec = .false.
    call prompt(' plot particles if no renderings?',iplotpartvec)
    backgnd_vec = .false.
    call prompt('use background color for vector plot?',backgnd_vec)
 endif
    
 return
end subroutine options_render
