!--------------------
!    MAIN MENU
!--------------------
subroutine menu(quit)
  use filenames
  use labels
  use settings
  use multiplot
  use prompting
  use transforms
  implicit none
  logical, intent(out) :: quit
  integer :: i,icol,ihalf,iadjust,iaction,istep,ierr,index
  integer :: ipicky,ipickx,irender,ivecplot,int_from_string
  integer :: iamvecprev, ivecplottemp
  character(LEN=2) :: ioption
  character(LEN=30) :: filename,vecprompt
  logical :: iansx, iansy, ichange

  quit = .false.
!--set coordinate and vector labels (depends on coordinate system)
  if (icoords.ne.0) then
     label(1:ndim) = labelcoord(1:ndim,icoordsnew)
     do i=1,numplot
        if (iamvec(i).ne.0) then
	   label(i) = trim(labelvec(iamvec(i)))//'\d'//labelcoord(i-iamvec(i)+1,icoordsnew)
	endif
     enddo
  endif
!--set contents of the vector plotting prompt
  vecprompt(1:6) = '0=none'
  index = 7
  iamvecprev = 0
  do icol=1,numplot
     if (iamvec(icol).ne.0 .and. iamvec(icol).ne.iamvecprev) then
	iamvecprev = iamvec(icol)
	if (iamvec(icol).ge.10) then
	   write(vecprompt(index:),"(',',1x,i2,'=',a)") &
	         iamvec(icol),trim(labelvec(icol))
	else
	   write(vecprompt(index:),"(',',1x,i1,'=',a)") &	
	         iamvec(icol),trim(labelvec(icol))
	endif
	index = len_trim(vecprompt) + 1
     endif
  enddo 

!---------------------------------------------------------------------------
!  print menu
!---------------------------------------------------------------------------
!
!  data columns
!
  print*,' You may choose from a delectable sample of plots '
  print 12
  ihalf = numplot/2		! print in two columns
  iadjust = mod(numplot,2)
  print 11, (i,transform_label(label(i),itrans(i)), &
       ihalf + i + iadjust, transform_label(label(ihalf + i + iadjust), &
       itrans(ihalf+i+iadjust)),i=1,ihalf)
  if (iadjust.ne.0) then
     print 13, ihalf + iadjust,transform_label(label(ihalf + iadjust), &
          itrans(ihalf+iadjust))
  endif
!
!  multiplot
!  
  print 12
  print 18,numplot+1,'multiplot ',nyplotmulti,'m','set multiplot '
18 format(1x,i2,')',1x,a,'[ ',i2,' ]',5x,a2,') ',a)
!
!  options 
! 
  print 12
  if (ishowopts) then
     print 14,'d','data options'
!!     print 15,'t','change number of timesteps read  ',(n_end-nstart+1)/nfreq 
     print 16,'i','toggle interactive mode          ',interactive
     print 14,'p','page options'
     print 14,'o','particle plot options'
     print 14,'l','change plot limits'
     print 14,'r','rendering options'
     print 14,'v','vector plot options'
     print 14,'x','cross sectioning / rotation options'
     print 14,'h','hide options'
     print 14,'s','save defaults'
     print 14,'q','exit supersphplot'
  else
     print*,' d(ata) i(nteractive) p(age) o(pts) l(imits) h(elp)'
     print*,' r(ender) v(ector) x(sec/rotate) s(ave) q(uit)'
  endif
 
  print 12

11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(1x,55('-'))
13 format(1x,i2,')',1x,a)
14 format(1x,a2,')',1x,a)
15 format(1x,a2,')',1x,a,'( ',i2, ' )')
16 format(1x,a2,')',1x,a,'( ',L1,' )')
17 format(1x,a2,')',1x,a20,1x,a2,')',1x,a)
!
!--prompt user for selection
!
9901 continue
  write(*,"(a)",ADVANCE='NO') 'Please enter your selection now (y axis or option):'
  read(*,*,ERR=9901) ioption 

!------------------------------------------------------------
!  if input is an integer and within range, plot data
!------------------------------------------------------------
  ipicky = int_from_string(ioption)
  if (ipicky.gt.0 .and. ipicky.le.numplot+1) then
     
     if (.not.ivegotdata) then
        print*,' no data '
        ihavereadfilename = .false.
        call get_data
     else	
        !
        !--if needed prompt for x axis selection
        !
        if (ipicky.le.(numplot-nextra)) then
	   if (ipickx.eq.0) ipickx = 1 ! do not allow zero as default
           call prompt(' (x axis) ',ipickx)
           if (ipickx.gt.numplot .or. ipickx.le.0) then
              goto 9901
           elseif (ipicky.le.ndim .and. ipickx.le.ndim) then
              call prompt('(render) (0=none)',irender,0,numplot)
	      ivecplottemp = -1
	      do while(.not.any(iamvec(1:numplot).eq.ivecplottemp).and.ivecplottemp.ne.0)
	         ivecplottemp = ivecplot
	         call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
		 if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
		    print "(a)",'Error, value not in list' 
		 endif
	      enddo
	      ivecplot = ivecplottemp
           else
              irender = 0
	      ivecplot = 0
           endif
	endif
	!
	!--call main plotting routine
        !
        call main(ipicky,ipickx,irender,ivecplot)
     endif
!------------------------------------------------------------------------
!  if input is an integer > numplot+1, quit
!------------------------------------------------------------------------
  elseif (ipicky.gt.numplot+1) then
     quit = .true.
  else
!------------------------------------------------------------------------
!  if input is a string, use menu options
!------------------------------------------------------------------------
     select case(adjustl(ioption))
!------------------------------------------------------------------------
     case('m','M')
        call prompt('Enter number of plots per timestep:',nyplotmulti,1,numplot)
        nacross = nyplotmulti/2
        if (nacross.eq.0) nacross = 1
        ndown = nyplotmulti/nacross
        print*,'setting nacross,ndown = ',nacross,ndown 
        iansx = .true.
        call prompt('Same x axis for all?',iansx)
        if (iansx) then
           call prompt('Enter x axis for all plots',multiplotx(1),1,numplot)
           multiplotx(2:nyplotmulti) = multiplotx(1)	     
        endif
        iansy = .false.
        if (ndim.ge.2) call prompt('Same y axis for all?',iansy)
        if (iansy) then
           call prompt('Enter y axis for all plots',multiploty(1),1,numplot)
           multiploty(2:nyplotmulti) = multiploty(1)
        endif
        
        do i=1,nyplotmulti
           print*,'Plot number ',i,':'
           if (.not.iansy .or. multiploty(i).gt.ndataplots .or. multiploty(i).le.0) then
              call prompt(' y axis ',multiploty(i),1,numplot)
           endif
           if (.not.iansx.and.multiploty(i).le.ndataplots) then
              call prompt(' x axis ',multiplotx(i),1,numplot)
           endif
           if ((multiplotx(i).le.ndim).and.(multiploty(i).le.ndim)) then
              call prompt('(render) (0=none)',irendermulti(i),0,numplot)
	      if (irendermulti(i).ne.0) then
                 ichange = .false.
		 call prompt(' change rendering options for this plot? ',ichange)
                 if (ichange) then
	            call prompt('plot contours? ',iplotcontmulti(i))
	            if (ndim.ge.2) then
		       call prompt(' cross section (no=projection)? ',x_secmulti(i))
		       if (x_secmulti(i)) then
		          call prompt('enter co-ordinate location of cross section slice',xsecposmulti(i))
		       endif
	            endif
                 elseif (i.eq.1) then
                    print*,'copying options from rendering settings'
                    iplotcontmulti(i) = iplotcont_nomulti
                    x_secmulti(i) = xsec_nomulti
                    xsecposmulti(i) = xsecpos_nomulti
                 else  
                    print*,'using same rendering options as plot 1'       
                    iplotcontmulti(i) = iplotcontmulti(1)
                    x_secmulti(i) = x_secmulti(1)
                    xsecposmulti(i) = xsecposmulti(1)
		 endif
 	      endif
              call prompt('(render) (0=none)',irender,0,numplot)
	      ivecplottemp = -1
	      do while(.not.any(iamvec(1:numplot).eq.ivecplottemp).and.ivecplottemp.ne.0)
	         ivecplottemp = ivecplot
	         call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
		 if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
		    print "(a)",'Error, value not in list' 
		 endif
	      enddo
	      ivecplotmulti(i) = ivecplottemp
           endif
        enddo
        return	    	  	  
!------------------------------------------------------------------------
     case('d','D')
        call options_data
!------------------------------------------------------------------------
     case('i','I')
        interactive = .not.interactive
        print*,' Interactive mode = ',interactive
        return
!------------------------------------------------------------------------
     case('p','P')
        call options_page	 
!------------------------------------------------------------------------
     case('o','O')
        call options_particleplots
        return
!------------------------------------------------------------------------
     case('r','R')
        call options_render
        return
!------------------------------------------------------------------------
     case('v','V')
        call options_vecplot
        return
!------------------------------------------------------------------------
     case('x','X')
        call options_xsecrotate
        return	
!------------------------------------------------------------------------
     case('l','L')
        call options_limits
        return	  
!------------------------------------------------------------------------
     case('s','S')
        call defaults_write
        return
!------------------------------------------------------------------------
     case('h','H')
        ishowopts = .not.ishowopts
        return
!------------------------------------------------------------------------
     case('q','Q')
        quit = .true.
        return
!------------------------------------------------------------------------
     case DEFAULT
        goto 9901 
     end select
     
  endif
  
  return      
end subroutine menu
