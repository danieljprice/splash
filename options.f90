!!
!!    implements menu options
!!
subroutine options(ipicky)      
  use exact_params
  use filenames
  use settings
  use multiplot
  use particle_data
  use prompting
  implicit none
  integer, intent(IN) :: ipicky
  integer :: i,j,k,n,iaction,ipick
  real :: diff, mid, temp
  character(LEN=30) :: filename
  character(LEN=1) :: ans
  logical :: ians, iansx, iansy, ichange

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
     call read_data(rootname,nfilesteps)

     nstart = 1
     n_end = nfilesteps
     !
     !--calculate various additional quantities
     !     
     call calc_quantities	  
     !
     !--set plot limits
     !
     call limits
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
     interactive = .not.interactive
     print*,' Interactive mode = ',interactive
     return
!------------------------------------------------------------------------
  case(5)
     animate = .not.animate
     print*,' Animate = ',animate
     return
!------------------------------------------------------------------------
  case(6)
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
     enddo
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
     call options_exact(iexact)
     return
!------------------------------------------------------------------------
  case(9)
     call options_page	 
!------------------------------------------------------------------------
  case(10)
     call options_particleplots
     return
!------------------------------------------------------------------------
  case(11)
     call options_render(npix_nomulti,icolours, &
          iplotcont_nomulti,ncontours_nomulti,     &
          ivecplot_nomulti,npixvec_nomulti,iplotpartvec_nomulti, &
          xsec_nomulti,xsecpos_nomulti,backgnd_vec_nomulti,      &
          ndim,numplot)
     return
!------------------------------------------------------------------------
  case(12)
     call options_limits
     return
!------------------------------------------------------------------------
  case(13)
     !	  show/hide plot options
     ishowopts = .not.ishowopts
     return 	  
!------------------------------------------------------------------------
  case(14)
     call defaults_write
     return
  case DEFAULT
     stop ' Error: menu action not defined'  

  end select

  return      
end subroutine options
