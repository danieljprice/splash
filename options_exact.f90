!----------------------------------------------------------------------
! sets options and parameters for exact solution calculation/plotting
!----------------------------------------------------------------------
subroutine options_exact(iexact)
 use exact_params
 use prompting
 implicit none
 integer, intent(inout) :: iexact
 integer :: i
 logical :: ians, iexist
 character(len=1) :: ans,dummy
 
 if (iexact.ne.0) then
    iexact = 0
 else   
    print 10
10  format(' 0) none ',/, 		&
           ' 1) polytrope ',/,		&
           ' 2) soundwave ',/,     	&
           ' 3) sedov blast wave',/,     &
           ' 4) toy star ',/,		&
           ' 5) mhd shock tubes ')
    call prompt('enter exact solution to plot',iexact,0,5)
    print*,' plotting exact solution number ',iexact
!
!--enter parameters for various exact solutions
!
    select case(iexact)
       case(1)
!
!--read exact solution for a polytrope
!
          ipolyc = ipolycmax
	  inquire (exist = iexist, file='polycalc.dat')
	  if (.not.iexist) then
	     print*,'ERROR: file polycalc.dat does not exist'
             return
	  endif
	  open(unit=14,file='polycalc.dat',status='old',form='formatted')
             read(14,20) dummy
20           format(a)      
             read(14,*) maxrho,mtot
             read(14,*) akfac
             read(14,*, end=100) (den(i),rad(i),i=1,ipolyc) 
          close(14)
          goto 101
100       continue
          print*,'end of polycalc.dat, i=',i-1
          ipolyc = i-1
          close(14)
101       continue         
       case(2)
          call prompt('enter wavelength of sound wave lambda ',lambda,0.0)		
          call prompt('enter amplitude ',delta,0.0)  
       case(3)
          call prompt('enter density of ambient medium ',rhosedov,0.0)
          call prompt('enter blast wave energy E ',esedov,0.0)
       case(4)
          print*,' toy star: '
	  call prompt('enter parameter a (v = ax) ',atstar)
          call prompt('enter parameter h (\rho = h - cx^2)',htstar)
	  call prompt('enter parameter c (\rho = h - cx^2)',ctstar,0.0)		
          sigma = 0.
	  call prompt('enter parameter sigma (By = sigma \rho)',sigma0)
	  sigma = sigma0
	  ians = .false.
	  call prompt('do you want oscillations?',ians)
	  norder = -1
	  if (ians) call prompt('enter order',norder,0)
    end select
 endif
 return
end subroutine options_exact
