!----------------------------------------------------------------------
! sets options and parameters for exact solution calculation/plotting
!----------------------------------------------------------------------
subroutine options_exact(iexact)
 use exact_params
 use prompting
 implicit none
 integer, intent(inout) :: iexact
 integer :: i, ierr
 logical :: ians, iexist
 character(len=1) :: ans,dummy
 character(len=30) :: filename
 
 print 10
10  format(' 0) none ',/, 		&
           ' 1) shock tube ',/,     	&
           ' 2) sedov blast wave ',/,   &
           ' 3) polytrope ',/,		&
           ' 4) toy star ',/,		&
	   ' 5) linear wave ',/,        &
           ' 6) mhd shock tubes (tabulated) ',/,    &
	   ' 7) read from file ')
 call prompt('enter exact solution to plot',iexact,0,7)
 print*,' plotting exact solution number ',iexact
!
!--enter parameters for various exact solutions
!
 select case(iexact)
 case(1)
!
!--read shock parameters from the .shk file
!
    call read_exactparams(iexact,ierr)
    if (ierr.ne.0) then
       call prompt('enter density to left of shock   ',rho_L,0.0)
       call prompt('enter density to right of shock  ',rho_L,0.0)   
       call prompt('enter pressure to left of shock  ',pr_L,0.0)
       call prompt('enter pressure to right of shock ',pr_R,0.0)
       call prompt('enter velocity to left of shock  ',v_L)
       call prompt('enter velocity to right of shock ',v_R)
    endif
 case(2)
    call prompt('enter density of ambient medium ',rhosedov,0.0)
    call prompt('enter blast wave energy E ',esedov,0.0)
 case(3)
    call prompt('enter polytropic k ',polyk) 
 case(4)
    print*,' toy star: '
    call read_exactparams(iexact,ierr)
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
 case(5)
    call prompt('enter y-plot to place sine wave on',iwaveploty,1)
    call prompt('enter x-plot to place sine wave on',iwaveplotx,1)
    call prompt('enter wavelength lambda ',lambda,0.0)		
    call prompt('enter amplitude ',ampl,0.0)
    call prompt('enter period ',period)
 case(6)
    print*,' MHD shock tube tables: '
    call prompt('enter solution to plot ',ishk,0,7)
 case(7)
    call prompt('enter filename: ',filename)
    call exact_fromfile(filename,ierr)
    if (ierr.gt.0) then
       iexact = 0
    else
       call prompt('enter x axis of exact solution: ',iexactplotx,1)
       call prompt('enter x axis of exact solution: ',iexactploty,1)
    endif   
 end select

 return
end subroutine options_exact
