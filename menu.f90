!--------------------
!    MAIN MENU
!--------------------
subroutine menu(quit)
  use filenames
  use labels
  use settings
  use multiplot
  use prompting
  implicit none
  logical, intent(out) :: quit
  integer :: i,ihalf,iadjust,iaction,istep,ierr
  integer :: ipicky,ipickx,irender,int_from_string
  character(LEN=2) :: ioption
  character(LEN=30) :: filename
  character(LEN=25) :: transform_label      
  logical :: iansx, iansy, ichange

  quit = .false.

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
     print 14,'d','read new data'
     print 15,'t','change number of timesteps read  ',(n_end-nstart+1)/nfreq 
     print 16,'i','toggle interactive mode          ',interactive
     print 14,'p','page options'
     print 14,'o','particle plot options'
     print 15,'r','rendering/vector plot options    ',ivecplot_nomulti
     print 14,'l','change plot limits'
     print 14,'h','hide options'
     print 14,'s','save defaults'
     print 14,'q','exit supersphplot'
  else
     print*,' d(ata) t(imesteps) i(nteractive) p(age) o(pts)'
     print*,' r(endering) l(imits) s(ave) h(elp) q(uit)'
  endif
 
  print 12

11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(1x,50('-'))
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
           else
              irender = 0
           endif
	endif
	!
	!--call main plotting routine
        !
        call main(ipicky,ipickx,irender)
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
     select case(ioption(1:1))
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
              call prompt(' enter field (from menu) for rendering (0=none)', &
                   irendermulti(i),0,numplot)
              if (irendermulti(i).ne.0) then
                 ichange = .false.
                 call prompt(' change rendering options for this plot? ',ichange)
                 if (ichange) then
                    call options_render( &
                         npixmulti(i),icolours,iplotcontmulti(i), &
                         ncontoursmulti(i),ivecplotmulti(i),npixvecmulti(i), &
                         iplotpartvecmulti(i),x_secmulti(i), &
                         xsecposmulti(i),backgnd_vec_multi(i),ndim,numplot)
                 elseif (i.eq.1) then
                    print*,'copying options from rendering settings'
                    npixmulti(i) = npix_nomulti	      
                    iplotcontmulti(i) = iplotcont_nomulti
                    ncontoursmulti(i) = ncontours_nomulti
                    ivecplotmulti(i) = ivecplot_nomulti
                    iplotpartvecmulti(i) = iplotpartvec_nomulti
                    x_secmulti(i) = xsec_nomulti
                    xsecposmulti(i) = xsecpos_nomulti
                    backgnd_vec_multi(i) = backgnd_vec_nomulti              
                 else  
                    print*,'using same rendering options as plot 1'       
                    npixmulti(i) = npixmulti(1)
                    iplotcontmulti(i) = iplotcontmulti(1)
                    ncontoursmulti(i) = ncontoursmulti(1)
                    ivecplotmulti(i) = ivecplotmulti(1)
                    iplotpartvecmulti(i) = iplotpartvecmulti(1)
                    x_secmulti(i) = x_secmulti(1)
                    xsecposmulti(i) = xsecposmulti(1)
                    backgnd_vec_multi(i) = backgnd_vec_multi(1)              
                 endif
              endif
           endif
        enddo
        return	    	  	  
!------------------------------------------------------------------------
     case('d','D')
        call get_data
!------------------------------------------------------------------------	  
     case('t','T')
        call prompt('Start at timestep ',nstart,1,nfilesteps)
        call prompt('End at timestep   ',n_end,nstart,nfilesteps)
        call prompt(' Frequency of steps to read',nfreq,1,nfilesteps)
        print *,' Steps = ',(n_end-nstart+1)/nfreq
        return
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
        call options_render(npix_nomulti,icolours, &
             iplotcont_nomulti,ncontours_nomulti,     &
             ivecplot_nomulti,npixvec_nomulti,iplotpartvec_nomulti, &
             xsec_nomulti,xsecpos_nomulti,backgnd_vec_nomulti,      &
             ndim,numplot)
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
