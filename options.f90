!!
!!    implements menu options
!!
subroutine options(ipicky)      
  use exact_params
  use filenames
  use labels
  use settings
  use multiplot
  use particle_data
  use prompting
  implicit none
  integer, intent(IN) :: ipicky
  integer :: i,j,k,n,iaction,ipick
  real :: diff, mid, temp
  real :: papersizey
  character(LEN=30) :: filename
  character(LEN=25) :: transform_label,crap
  character(LEN=1) :: ans
  logical :: ians, ichange

  iaction = ipicky - numplot      

  select case(iaction)

!------------------------------------------------------------------------
  case(2)
     if (.not.ihavereadfilename) then
        call prompt('Enter filename to read',rootname)
     endif
     ihavereadfilename = .false.
     nfilesteps = maxstep
     !
     !--read the data from the file
     !
     call read_data(rootname,ifile,gamma(:),time(:), &
          dat(:,:,:),iam(:), &
          ncolumns,ncalc,nfilesteps,ntot(:),npart(:), &
          nghost(:),ndim,ndimV,hfact,ivegotdata)
     !
     !--calculate various additional quantities
     !     
     if (ncalc.gt.0) call calc_quantities	  
     !
     !--set plot limits
     !
     print*,'Setting plot limits...'
     n_end = nfilesteps
     ntotplot(:) = npart(:)
     if (iplotghost) ntotplot(:) = ntot(:)
     zoom = 1.0
     !!--find limits of particle properties	  
     lim(:,1) = 1.e6
     lim(:,2) = -1.e6
     do j=1,numplot
        do i=nstart,n_end
           do k=1,ntotplot(i)
              lim(j,1) = min(lim(j,1),dat(j,k,i))
              lim(j,2) = max(lim(j,2),dat(j,k,i)*scalemax)
           enddo
        enddo
        if (lim(j,2).eq.lim(j,1)) then
           print*,label(j),' min = max = ',lim(j,1)
        endif
     enddo
     !!--limits of magnetic field (in all dimensions) for vector plot
     if (iBfirst.ne.0) then
        if (iadapt) then
           Bmin = 0.0
           Bmax = maxval(dat(iBfirst:iBlast,1:ntotplot(1),1))*scalemax
        else
           Bmax = 0.0
           Bmin = 0.0
           do j=iBfirst,iBlast
              !Bmin = min(Bmin,lim(j,1))
              if (lim(j,2).gt.Bmax) then
                 Bmax=lim(j,2)
                 print*,' Bmax = ',label(j),lim(j,2)
              endif
           enddo
        endif
        print*,'Bmin,Bmax = ',Bmin,Bmax
     endif
     print*,'plot limits set'
     !
     !--read toy star file for toy star solution
     !
     if (iexact.eq.4) then
        filename = trim(rootname)//'.tstar'
        open(UNIT=20,ERR=8801,FILE=filename,STATUS='old')
        read(20,*,ERR=8801) Htstar,Ctstar,Atstar
        read(20,*,ERR=8801) sigma0
        read(20,*,ERR=8801) norder
        close(UNIT=20)
        print*,' >> read ',filename
        print*,' H,C,A,sigma,n = ',Htstar,Ctstar,Atstar,sigma0,norder
        return
8801    continue
        print*,'Cannot open/ error reading ',filename	     
     endif
     return
!------------------------------------------------------------------------	  
  case(3)
     call prompt('Start at timestep ',nstart,1,nfilesteps)
     call prompt('End at timestep   ',n_end,nstart)
     if (n_end.gt.nfilesteps) then
        print*,'n_end greater than nfilesteps, reset to ',nfilesteps
        n_end = nfilesteps
     endif
     call prompt(' Frequency of steps to read',nfreq,1,nfilesteps)
     print *,' Steps = ',(n_end-nstart+1)/nfreq
     return
!------------------------------------------------------------------------
  case(4)
     axes=.not.axes
     print *,' Axes = ',axes
     return	  
!------------------------------------------------------------------------
  case(5)
     animate=.not.animate
     print *,' Animation = ',animate
     return
!------------------------------------------------------------------------
  case(6)
     iadapt=.not.iadapt
     print *,' Adaptive plot limits = ',iadapt
     return
!------------------------------------------------------------------------
  case(7)
     xsec_nomulti =.not.xsec_nomulti
     print *,' Cross section = ',xsec_nomulti
     flythru = .false.
     if (xsec_nomulti) then
        call prompt('Do you want a fly-through',flythru)
        !	     READ*,ans
        !	     IF (ans.eq.'y'.or.ans.eq.'Y') flythru=.true.
     endif
     return
!------------------------------------------------------------------------
  case(8)
     plotcirc=.not.plotcirc
     print*,' Plot circles of interaction = ',plotcirc
     if (plotcirc) then	     
        call prompt('Plot all circles?',plotcircall)	     
        if (.not.plotcircall) then
           call prompt('Enter number of circles to draw',ncircpart, &
                1,size(icircpart))
           do n=1,ncircpart
              call prompt('Enter particle number to plot circle around', &
                   icircpart(n),1,maxval(ntot))
           enddo
        endif
     endif
     return 	  
!------------------------------------------------------------------------
  case(9)
     print*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
     call prompt(' Enter PGPLOT marker # (particles):',imark)
     call prompt(' Enter PGPLOT marker # (ghosts)   :',imarkg)
     call prompt(' Enter PGPLOT marker # (sinks)    :',imarksink)
     return
!------------------------------------------------------------------------
  case(10)
     call prompt('Enter number of plots across:',nacross,1,numplot)
     call prompt('Enter number of plots down  :',ndown,1,numplot)
     return	    	  
!------------------------------------------------------------------------
  case(11)
     call prompt('Enter number of plots per timestep:',nyplotmulti,1,numplot)
     !	  READ*,nyplotmulti
     !          IF (nyplotmulti.GT.numplot) THEN 
     !	     nyplotmulti = numplot
     !	     PRINT*,'nyplots > numplot, reset to ',numplot
     !	  ENDIF
     nacross = nyplotmulti/2
     if (nacross.eq.0) nacross = 1
     ndown = nyplotmulti/nacross
     print*,'setting nacross,ndown = ',nacross,ndown 
     ians = .true.
     call prompt('Same x axis for all?',ians)
     if (ians) then
        call prompt('Enter x axis for all plots',multiplotx(1),1,numplot)
        multiplotx(2:nyplotmulti) = multiplotx(1)	     
     endif
     do i=1,nyplotmulti
        print*,'Plot number ',i,':'
        call prompt(' y axis ',multiploty(i),1,numplot)
        if (.not.ians.and.multiploty(i).le.ndataplots) then
           call prompt(' x axis ',multiplotx(i),1,numplot)
        endif
        if ((multiplotx(i).le.ndim).and.(multiploty(i).le.ndim)) then
           call prompt(' enter field (from menu) for rendering (0=none)', &
                irendermulti(i),0,numplot)
           ichange = .false.
           call prompt(' change options for this plot? ',ichange)
           if (ichange) then
              call options_render(irendermulti(i),  &
                   npixmulti(i),icolours,iplotcontmulti(i), &
                   ncontoursmulti(i),ivecplotmulti(i),npixvecmulti(i), &
                   iplotpartvecmulti(i),x_secmulti(i), &
                   xsecposmulti(i),backgnd_vec_multi(i),ndim,numplot)
           elseif (i.gt.1) then          
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
     enddo
     return	    	  
!------------------------------------------------------------------------
  case(12)
     ipagechange=.not.ipagechange
     print*,' Page changing = ',ipagechange
     return 	  
!------------------------------------------------------------------------
  case(13)
     print*,' 0) PGPLOT default'
     print*,' 1) small square movie '
     print*,' 2) large/multiple movie'
     print*,' 3) single small graph'
     print*,' 4) duo small graph '
     print*,' 5) Custom size ' 
     call prompt(' Enter option for paper size ',ipapersize,0,5)
     select case(ipapersize)
     case(1) 
        papersizex = 0.25*11.7
        aspectratio = 1.0
     case(2)
        papersizex = 0.5*11.7
        aspectratio = 1.0
     case(3) 
        papersizex = 0.5*11.7 
        aspectratio = 1./sqrt(2.)
     case(4)
        papersizex = 11.7
        aspectratio = 0.5/sqrt(2.)	
     case(5)
        call prompt(' x size (inches) ',papersizex,0.0,12.0)
        call prompt(' y size (inches) or aspect ratio (-ve)', &
             papersizey,-12.0,12.0)
        if (papersizey.lt.0.0) then
           aspectratio = abs(papersizey)
        else
           aspectratio = papersizey/papersizex
        endif
     case DEFAULT
        papersizex = 0.0	! no call to PGPAP if they are zero
        aspectratio = 0.0	
     end select
     return 	  

!------------------------------------------------------------------------
  case(14)
     print*,' Plot initial only(i), all(a), both(b) or not (n)?'
     read*,ans
     iplotline = .false.
     iplotlinein = .false.
     if (ans.eq.'i'.or.ans.eq.'b') iplotlinein = .true.
     if (ans.eq.'a'.or.ans.eq.'b') iplotline = .true.
     if (iplotlinein) then
        call prompt('Enter PGPLOT line style',linestylein,0,5)
     endif
     print*,' Plot line = ',iplotline,iplotlinein
     return 	  
!------------------------------------------------------------------------
  case(15)
     call options_exact(iexact)
!------------------------------------------------------------------------
  case(16)
     iplotav=.not.iplotav
     if (iplotav) then
        call prompt('Enter no. of bins for averaging ',nbins,1,1000)
     endif
     print*,' Plot average, nbins = ',iplotav,nbins
     return 		  	    	  	  
!-----------------------------------------------------------------------
  case(17)
     !	  label particles with particle numbers
     ilabelpart=.not.ilabelpart
     print*,' label particles = ',ilabelpart
     return 	  
!------------------------------------------------------------------------
  case(18)
     !	  plot ghost particles?
     call prompt('Plot ghost particles? ',iplotghost)
     call prompt('Plot sink particles? ',iplotsink)
     print*,' plot ghost particles = ',iplotghost
     print*,' plot sink particles = ',iplotsink
     if (iplotghost) ntotplot(:) = npart(:) + nghost(:)
     return 	  
!------------------------------------------------------------------------
  case(19)
     call options_render(npix_nomulti,icolours, &
          iplotcont_nomulti,ncontours_nomulti,     &
          ivecplot_nomulti,npixvec_nomulti,iplotpartvec_nomulti, &
          xsec_nomulti,xsecpos_nomulti,backgnd_vec_nomulti,      &
          ndim,numplot)
     return
!------------------------------------------------------------------------
  case(20)
     ipick = 1
     do while (ipick.gt.0 .and. ipick.le.numplot)
        ipick = 0
        print*,'Enter plot number to apply transformation '
        call prompt('(0 = finish, -1 = set all) ',ipick)
        if (ipick.le.numplot .and. ipick.ne.0) then
           write(*,300) (i,trim(transform_label('x',i)), i=1,5)
300        format(1x,i1,') ',a)		
           print*,'Enter transformation to apply (or a combination e.g. 321)'
           if (ipick.lt.0) then
              ipick = 0
              call prompt(' ',ipick,0)
              itrans(:) = ipick
              ipick = -99
           else
              call prompt(' ',itrans(ipick),0)
           endif
        endif
     enddo
     return
!------------------------------------------------------------------------
  case(21)
     if (.not.iadapt) then
        ians = .false.
        call prompt('Do you want to manually enter limits?',ians)
        if (ians) then
301        print*,'Enter plot number to set limits (0 to finish)'
           read*,ipick
           if ((ipick.le.numplot).and.(ipick.gt.0)) then
              print*,' Current limits   (min,max):', &
                   lim(ipick,1),lim(ipick,2)
              print*,' Enter new limits (min,max):'
              read*,lim(ipick,1),lim(ipick,2)
              print*,' >> limits set    (min,max):', &
                   lim(ipick,1),lim(ipick,2)
           else
              print*,'Finished, returning'
              return
           endif
           goto 301
        else
           call prompt('Enter zoom factor for fixed limits',zoom,0.0)
           do i=1,numplot
              diff = lim(i,2)- lim(i,1)
              mid = 0.5*(lim(i,1) + lim(i,2))
              lim(i,1) = mid - 0.5*zoom*diff
              lim(i,2) = mid + 0.5*zoom*diff
           enddo
        endif
     else
        call prompt('Enter scale factor (adaptive limits)',scalemax,0.0)
     endif
     return
!------------------------------------------------------------------------
  case(22)
     !	  show/hide plot options
     ishowopts = .not.ishowopts
     return 	  
!------------------------------------------------------------------------
  case(23)
     call defaults_write
     return
  case DEFAULT
     stop ' Error: menu action not defined'  

  end select

  return      
end subroutine options
