!!
!! sub-menu with utilities relating to particle plots
!!
subroutine options_particleplots   
  use settings
  use particle_data
  use prompting
  implicit none
  integer :: iaction,n
  character(LEN=1) :: ans

  iaction = 0      
    print 10, iplotline,iplotlinein,iplotav,ilabelpart,plotcirc, &
         iplotghost,iplotsink,imark,imarkg
10  format(' 0) exit ',/, 		&
         ' 1) toggle plot line              ( ',L1,',',1x,L1,' ) ',/, &
         ' 2) toggle plot average line      ( ',L1,' ) ',/,           &
         ' 3) toggle label particles        ( ',L1,' ) ',/,           &
         ' 4) toggle circles of interaction ( ',L1,' ) ',/,           &
         ' 5) toggle plot ghosts/sinks      ( ',L1,',',1x,L1,' )',/,  &
         ' 6) change graph markers          ( ',i2,',',1x,i2,' )')
    call prompt('enter option',iaction,0,6)
!
  select case(iaction)

!------------------------------------------------------------------------
  case(1)
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
  case(2)
     iplotav=.not.iplotav
     if (iplotav) then
        call prompt('Enter no. of bins for averaging ',nbins,1,1000)
     endif
     print*,' Plot average, nbins = ',iplotav,nbins
     return 		  	    	  	  
!-----------------------------------------------------------------------
  case(3)
     !	  label particles with particle numbers
     ilabelpart=.not.ilabelpart
     print*,' label particles = ',ilabelpart
     return 	  
!------------------------------------------------------------------------
  case(4)
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
  case(5)
     !	  plot ghost particles?
     call prompt('Plot ghost particles? ',iplotghost)
     call prompt('Plot sink particles? ',iplotsink)
     print*,' plot ghost particles = ',iplotghost
     print*,' plot sink particles = ',iplotsink
     if (iplotghost) ntotplot(:) = npart(:) + nghost(:)
     return 	  
!------------------------------------------------------------------------
  case(6)
     print*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
     call prompt(' Enter PGPLOT marker # (particles):',imark)
     call prompt(' Enter PGPLOT marker # (ghosts)   :',imarkg)
     call prompt(' Enter PGPLOT marker # (sinks)    :',imarksink)
     return
!------------------------------------------------------------------------
  case default
     return

  end select

  return      
end subroutine options_particleplots
