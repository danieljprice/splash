!
! This module handles all of the settings relating to the exact solution
! plotting and calls the appropriate routines to change these settings and
! plot the actual solutions.
!
! The only things to do with exact solutions which are not entirely handled
! by this module are the h vs rho relation (because it takes all of the
! particle masses as input) and the toy star AC plane solution (because
! it is called under different circumstances to the other solutions).
! Hence nearly all the parameters are private to the module, but some are
! declared public for these reasons.
!
module exact
  implicit none
  !
  !--declare all of the parameters required for the various exact solutions
  !
  !--toy star
  integer, public :: iACplane ! label position of toy star AC plane plot
  integer, private :: norder ! for toy star
  real, public :: htstar,atstar,ctstar,sigma
  real, private :: totmass, sigma0
  !--sound wave
  integer, public :: iwaveploty,iwaveplotx ! linear wave
  real, private :: ampl,lambda,period
  !--sedov blast wave
  real, private :: rhosedov,esedov
  !--polytrope
  real, private :: polyk
  !--mhd shock solutions
  integer, private :: ishk
  !--from file
  integer, private, parameter :: maxexactpts = 1001
  integer, private :: iexactpts, iexactplotx, iexactploty
  real, private, dimension(maxexactpts) :: xexact,yexact
  !--shock tube
  real, private :: rho_L, rho_R, pr_L, pr_R, v_L, v_R
  !--rho vs h
  real, public :: hfact
  !
  !--sort these into a namelist for input/output
  !
  namelist /exactparams/ ampl,lambda,period,iwaveploty,iwaveplotx, &
       htstar,atstar,ctstar,sigma0,norder,rhosedov,esedov, &
       rho_L, rho_R, pr_L, pr_R, v_L, v_R, hfact, &
       iexactplotx,iexactploty
  private :: exactparams

contains
  !----------------------------------------------------------------------
  ! sets default values of the exact solution parameters
  !----------------------------------------------------------------------
  subroutine defaults_set_exact
    implicit none

    lambda = 1.0    ! sound wave exact solution : wavelength
    ampl = 0.005    ! sound wave exact solution : amplitude
    period = 1.0
    iwaveploty = 7
    iwaveplotx = 1
    htstar = 1.     ! toy star crap
    atstar = 1.
    ctstar = 1.
    norder = 0
    sigma0 = 0.
    rhosedov = 1.0  ! sedov blast wave
    esedov = 1.0    ! blast wave energy
    polyk = 1.0     ! polytropic k
    rho_L = 1.0     ! shock tube (default is sod problem)
    rho_R = 0.125
    pr_L = 1.0
    pr_R = 0.1
    v_L = 0.0
    v_R = 0.0
    iexactpts = maxexactpts
    iexactplotx = 0
    iexactploty = 0
    ishk = 0
    
    return
  end subroutine defaults_set_exact

  !----------------------------------------------------------------------
  ! reads values of the exact solution parameters from the defaults file
  !----------------------------------------------------------------------
  subroutine defaults_read_exact(iunit,ierr)
    implicit none
    integer, intent(in) :: iunit
    integer, intent(out) :: ierr
    
    read(iunit,NML=exactparams,iostat=ierr)

    return
  end subroutine defaults_read_exact
  
  !----------------------------------------------------------------------
  ! writes values of the exact solution parameters to the defaults file
  !----------------------------------------------------------------------
  subroutine defaults_write_exact(iunit,ierr)
    implicit none
    integer, intent(in) :: iunit
    integer, intent(out) :: ierr
    
    write(iunit,NML=exactparams,iostat=ierr)

    return
  end subroutine defaults_write_exact

  !----------------------------------------------------------------------
  ! sets options and parameters for exact solution calculation/plotting
  !----------------------------------------------------------------------
  subroutine submenu_exact(iexact)
    use prompting
    implicit none
    integer, intent(inout) :: iexact
    integer :: i, ierr
    logical :: ians, iexist
    character(len=1) :: ans,dummy
    character(len=30) :: filename

    print 10
10  format(' 0) none ',/,               &
         ' 1) shock tube ',/,           &
         ' 2) sedov blast wave ',/,     &
         ' 3) polytrope ',/,            &
         ' 4) toy star ',/,             &
         ' 5) linear wave ',/,          &
         ' 6) mhd shock tubes (tabulated) ',/,  &
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
          call prompt('enter density to right of shock  ',rho_R,0.0)   
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
       call prompt('enter parameter h (rho = h - cx^2)',htstar)
       call prompt('enter parameter c (rho = h - cx^2)',ctstar,0.0)
       sigma = 0.
       call prompt('enter parameter sigma (By = sigma*rho)',sigma0)
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
       call exact_fromfile(filename,xexact,yexact,maxexactpts,iexactpts,ierr)
       if (ierr.gt.0) then
          iexact = 0
       else
          call prompt('enter x axis of exact solution: ',iexactplotx,1)
          call prompt('enter x axis of exact solution: ',iexactploty,1)
       endif
    end select

    return
  end subroutine submenu_exact

  !-----------------------------------------------------------------------
  ! read exact solution parameters from files
  ! (in ndspmhd these files are used in the input to the code)
  !
  ! called after main data read and if exact solution chosen from menu
  !-----------------------------------------------------------------------
  subroutine read_exactparams(iexact,ierr)
    use filenames
    implicit none
    integer, intent(in) :: iexact
    integer, intent(out) :: ierr
    integer :: int_from_string
    character(LEN=30) :: filename

    select case(iexact)
    case(1)
       !
       !--shock tube parameters from .shk file
       !
       filename = trim(rootname(1))//'.shk'
       open(UNIT=19,ERR=7701,FILE=filename,STATUS='old')
       read(19,*,ERR=7777) rho_L, rho_R
       read(19,*,ERR=7777) pr_L, pr_R
       read(19,*,ERR=7777) v_L, v_R
       close(UNIT=19)
       print*,'>> read ',filename
       print*,' rhoL, rho_R = ',rho_L,rho_R
       print*,' pr_L, pr_R  = ',pr_L, pr_R
       print*,' v_L,  v_R   = ',v_L, v_R
       return
7701   print*,'no file ',filename
       ierr = 1
       return
7777   print*,'error reading ',filename
       close(UNIT=19)
       ierr = 2
       return

    case(4)
       !
       !--read toy star file for toy star solution
       !
       filename = trim(rootname(1))//'.tstar'
       open(unit=20,ERR=8801,FILE=filename,STATUS='old')
       read(20,*,ERR=8888) Htstar,Ctstar,Atstar
       read(20,*,ERR=8888) sigma0
       read(20,*,ERR=8888) norder
       close(UNIT=20)
       print*,' >> read ',filename
       print*,' H,C,A,sigma,n = ',Htstar,Ctstar,Atstar,sigma0,norder
       return
8801   continue
       print*,'no file ',filename
       ierr = 1
       return
8888   print*,'error reading ',filename
       close(UNIT=20)
       ierr = 2
       return

    case(6)
       !
       !--attempt to guess which MHD shock tube has been done from filename
       !
       ishk = int_from_string(rootname(1)(5:5))
       return

    end select

    return
  end subroutine read_exactparams

  !-----------------------------------------------------------------------
  ! this subroutine drives the exact solution plotting using the
  ! parameters which have been set
  !-----------------------------------------------------------------------

  subroutine exact_solution(iplotx,iploty,iexact,ndim,time,xmin,xmax,ymean,gamma)
    use labels
    use limits
    use prompting
    implicit none
    integer, intent(in) :: iplotx,iploty,iexact,ndim
    real, intent(in) :: time,xmin,xmax,ymean,gamma

    select case(iexact)
    case(1)! shock tube
       if (iplotx.eq.ix(1)) then
          if (iploty.eq.irho) then
             call exact_shock(1,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
          elseif (iploty.eq.ipr) then
             call exact_shock(2,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
          elseif (iploty.eq.ivx) then
             call exact_shock(3,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
          elseif (iploty.eq.iutherm) then
             call exact_shock(4,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xmin,xmax)
          endif
       endif

    case(2)! sedov blast wave
       if (iplotx.eq.irad) then
          if (iploty.eq.irho) then
             call exact_sedov(time,gamma,rhosedov,esedov,lim(irad,2),1)
          elseif (iploty.eq.ipr) then
             call exact_sedov(time,gamma,rhosedov,esedov,lim(irad,2),2)                 
          elseif (iploty.eq.iutherm) then
             call exact_sedov(time,gamma,rhosedov,esedov,lim(irad,2),3)                
          elseif (iploty.eq.ike) then
             call exact_sedov(time,gamma,rhosedov,esedov,lim(irad,2),4)                 
          endif
       endif

    case(3)! polytrope
       if (iploty.eq.irho .and. iplotx.eq.irad) call exact_polytrope(gamma)

    case(4)! toy star
       if (iBfirst.ne.0) then
          sigma = sigma0
       else
          sigma = 0.
       endif
       if (ndim.eq.1) then
          !
          !--1D toy star solutions
          !
          if (iplotx.eq.ix(1) .or. iplotx.eq.irad) then! if x axis is x or r
             if (iploty.eq.irho) then
                call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,1)
             elseif (iploty.eq.ipr) then
                call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,2)       
             elseif (iploty.eq.iutherm) then
                call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,3)       
             elseif (iploty.eq.ivx) then
                call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,4)       
             elseif (iploty.eq.ibfirst+1) then
                call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,5)
             endif
          elseif (iplotx.eq.irho) then
             if (iploty.eq.ibfirst+1) then
                call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,6)       
             endif
          endif

          if (iploty.eq.iacplane) then! plot point on a-c plane
             call exact_toystar(time,gamma,htstar,atstar,ctstar,sigma,norder,7)
          endif
       else
          !
          !--2D and 3D toy star solutions
          !
          if ((iplotx.eq.ix(1) .and. iploty.eq.ivx) &
               .or. (iplotx.eq.ix(2) .and. iploty.eq.ivx+1)) then
             call exact_toystar2D(time,gamma, &
                  htstar,atstar,ctstar,sigma,norder,4)
          endif
          if (iplotx.eq.irad) then
             if (iploty.eq.irho) then
                call exact_toystar2D(time,gamma,htstar,atstar,ctstar,sigma,norder,1)
             elseif (iploty.eq.ipr) then
                call exact_toystar2D(time,gamma,htstar,atstar,ctstar,sigma,norder,2)
             elseif (iploty.eq.iutherm) then
                call exact_toystar2D(time,gamma,htstar,atstar,ctstar,sigma,norder,3)
             elseif (iploty.eq.ivx .or. iploty.eq.ivx+1) then
                call exact_toystar2D(time,gamma, &
                     htstar,atstar,ctstar,sigma,norder,4)
             elseif (iploty.eq.ike) then
                call exact_toystar2D(time,gamma, &
                     htstar,atstar,ctstar,sigma,norder,4)
             endif
          endif
       endif

    case(5)! linear wave
       if ((iploty.eq.iwaveploty).and.(iplotx.eq.iwaveplotx)) then
          call exact_wave(time,ampl,period,lambda,xmin,xmax,ymean)
       endif

    case(6) ! mhd shock tubes
       if (iplotx.eq.ix(1)) then
          !--prompt for shock type if not set  
          if (ishk.eq.0) then ! prompt
             call prompt('enter shock solution to plot',ishk,0,6)
          endif
          if (iploty.eq.irho) then
             call exact_mhdshock(1,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.ipr) then
             call exact_mhdshock(2,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.ivx) then
             call exact_mhdshock(3,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.ivx+1) then
             call exact_mhdshock(4,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.ivlast) then
             call exact_mhdshock(5,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.ibfirst+1) then
             call exact_mhdshock(6,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.iblast) then
             call exact_mhdshock(7,ishk,time,gamma,xmin,xmax)
          elseif (iploty.eq.iutherm) then
             call exact_mhdshock(8,ishk,time,gamma,xmin,xmax)
          endif
       endif

    case(7) ! exact solution read from file
       if (iplotx.eq.iexactplotx .and. iploty.eq.iexactploty) then   
          call pgline(iexactpts,xexact,yexact)
       endif
    end select

    return

  end subroutine exact_solution

end module exact
