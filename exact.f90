!
! This module handles all of the settings relating to the exact solution
! plotting and calls the appropriate routines to change these settings and
! plot the actual solutions.
!
! The only thing to do with exact solutions that is not entirely handled
! by this module is the toy star AC plane solution (because
! it is called under different circumstances to the other solutions).
!
module exact
  implicit none
  !
  !--options used to plot the exact solution line
  !
  integer :: maxexactpts, iExactLineColour, iExactLineStyle
  logical :: iApplyTransExactFile,iCalculateExactErrors,iPlotResiduals
  real :: fracinsetResiduals,residualmax
  !
  !--declare all of the parameters required for the various exact solutions
  !
  !--toy star
  integer :: iACplane ! label position of toy star AC plane plot
  integer :: norder,morder ! for toy star
  real, public :: atstar,ctstar,sigma
  real :: htstar,alphatstar,betatstar,ctstar1,ctstar2
  real :: sigma0,totmass
  !--sound wave
  integer :: iwaveploty,iwaveplotx ! linear wave
  real :: ampl,lambda,period,xzero
  !--sedov blast wave
  real :: rhosedov,esedov
  !--polytrope
  real :: polyk
  !--mhd shock solutions
  integer :: ishk
  !--density profiles
  integer :: iprofile,icolpoten,icolfgrav
  real, dimension(2) :: Msphere,rsoft
  !--from file
  integer :: iexactplotx, iexactploty
  !--shock tube
  real :: rho_L, rho_R, pr_L, pr_R, v_L, v_R
  !--rho vs h
  real :: hfact
  character(len=120) :: filename_exact
  !--equilibrium torus
  real :: Mstar,Rtorus,distortion
  
  !
  !--sort these into a namelist for input/output
  !
  namelist /exactopts/ iexactplotx,iexactploty,filename_exact,maxexactpts, &
       iExactLineColour,iExactLineStyle,iApplyTransExactFile,iCalculateExactErrors, &
       iPlotResiduals,fracinsetResiduals,residualmax

  namelist /exactparams/ ampl,lambda,period,iwaveploty,iwaveplotx,xzero, &
       htstar,atstar,ctstar,alphatstar,betatstar,ctstar1,ctstar2, &
       polyk,sigma0,norder,morder,rhosedov,esedov, &
       rho_L, rho_R, pr_L, pr_R, v_L, v_R,ishk,hfact, &
       iprofile,Msphere,rsoft,icolpoten,icolfgrav,Mstar,Rtorus,distortion
       
  public :: defaults_set_exact,submenu_exact,options_exact,read_exactparams
  public :: exact_solution
  public :: exactopts,exactparams

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
    xzero = 0.
    htstar = 1.     ! toy star crap
    atstar = 1.
    ctstar = 1.
    alphatstar = 0.
    betatstar = 0.
    ctstar1 = 0.
    ctstar2 = 0.
    totmass = 1.
    norder = -1
    morder = 0
    sigma0 = 0.
    rhosedov = 1.0  ! sedov blast wave
    esedov = 1.0    ! blast wave energy
    polyk = 1.0     ! polytropic k
!   shock tube (default is sod problem)
    rho_L = 1.0
    rho_R = 0.125
    pr_L = 1.0
    pr_R = 0.1
    v_L = 0.0
    v_R = 0.0
    iexactplotx = 0
    iexactploty = 0
    ishk = 1
    hfact = 1.2
    filename_exact = ' '
!   density profile parameters
    iprofile = 1
    rsoft(1) = 1.0
    rsoft(2) = 0.1
    Msphere(1) = 1.0
    Msphere(2) = 0.0
    icolpoten = 0
    icolfgrav = 0
!   equilibrium torus
    Mstar = 1.0
    Rtorus = 1.0
    distortion = 1.1

    maxexactpts = 1001      ! points in exact solution plot
    iExactLineColour = 1    ! foreground
    iExactLineStyle = 1     ! solid
    iApplyTransExactFile = .true. ! false if exact from file is already logged
    iCalculateExactErrors = .true.
    iPlotResiduals = .false.
    fracinsetResiduals = 0.15
    residualmax = 0.0
    
    return
  end subroutine defaults_set_exact

  !----------------------------------------------------------------------
  ! sets which exact solution to calculate + parameters for this
  !----------------------------------------------------------------------
  subroutine submenu_exact(iexact)
    use settings_data, only:ndim
    use prompting
    use filenames, only:rootname
    implicit none
    integer, intent(inout) :: iexact
    integer :: ierr
    logical :: ians,iexist

    print 10
10  format(' 0) none ',/,               &
         ' 1) shock tube ',/,           &
         ' 2) sedov blast wave ',/,     &
         ' 3) polytrope ',/,            &
         ' 4) toy star ',/,             &
         ' 5) linear wave ',/,          &
         ' 6) mhd shock tubes (tabulated) ',/,  &
         ' 7) h vs rho ',/, &
         ' 8) radial density profiles ',/, &
         ' 9) papaloizou & pringle torus ',/, &
         '10) read from file ')
    call prompt('enter exact solution to plot',iexact,0,10)
    print*,' plotting exact solution number ',iexact
    !
    !--enter parameters for various exact solutions
    !
    select case(iexact)
    case(1)
       !
       !--read shock parameters from the .shk file
       !
       call read_exactparams(iexact,trim(rootname(1)),ierr)
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
       print "(a)",' toy star: '
       call read_exactparams(iexact,trim(rootname(1)),ierr)
       call prompt('enter polytropic k ',polyk)
       call prompt('enter total mass   ',totmass)
       call prompt('enter central density rho_0 (rho = rho_0 - cr^2)',htstar)
       call prompt('enter parameter c (rho = rho_0 - cr^2)',ctstar,0.0)
       sigma = 0.
       call prompt('enter parameter sigma (By = sigma*rho)',sigma0)
       sigma = sigma0
       ians = .false.
       call prompt('linear oscillations?',ians)
       if (ians) then
          call prompt('enter order of radial mode',norder,0)
          if (ndim.ge.2) call prompt('enter order of angular mode',morder,0)
          call prompt('enter velocity amplitude a (v = a*r) ',atstar)
       else
          print "(a)",'using exact non-linear solution:'
          ians = .true.
          if (norder.lt.0 .and. morder.lt.0) ians = .false.
          if (ndim.ge.2) call prompt('axisymmetric?',ians)
          if (ians .or. ndim.eq.1) then
             norder = -1
             morder = 0
             call prompt('enter v_r amplitude ',alphatstar)
             if (ndim.ge.2) call prompt('enter v_phi amplitude ',betatstar)
          else
             norder = -1
             morder = -1
             call prompt('enter vxx amplitude ',alphatstar)
             call prompt('enter vyy amplitude ',betatstar)
             call prompt('enter vxy amplitude ',ctstar1)
             call prompt('enter vyx amplitude ',ctstar2)
          endif
       endif
    case(5)
       call prompt('enter y-plot to place sine wave on',iwaveploty,1)
       call prompt('enter x-plot to place sine wave on',iwaveplotx,1)
       call prompt('enter starting x position',xzero)
       call prompt('enter wavelength lambda ',lambda,0.0)
       call prompt('enter amplitude ',ampl,0.0)
       call prompt('enter period ',period)
    case(6)
       print "(a)",' MHD shock tube tables: '
       if (ishk.le.0) ishk = 1
       call prompt('enter solution to plot ',ishk,1,7)
    case(7)
       call prompt('enter hfact [h = hfact*(m/rho)**1/ndim]',hfact,0.)
    case(8)
       print 20
20     format(' 1) Plummer sphere  [ rho = 3M r_s**2 /(4 pi (r**2 + r_s**2)**5/2) ]',/, &
              ' 2) Hernquist model [ rho =     M r_s /(2 pi r (r_s + r)**3        ]')
       call prompt('enter density profile to plot',iprofile,1,2)
       call prompt('enter total mass of sphere M',Msphere(1),0.)
       call prompt('enter scale length length r_s,',rsoft(1),0.)
       ians = .false.
       if (icolpoten.gt.0) ians = .true.
       call prompt('Are the gravitational potential and/or force dumped?',ians)
       if (ians) then
          call prompt('enter column containing grav. potential',icolpoten,0)
          call prompt('enter column containing grav. force',icolfgrav,0)
       endif
       call prompt('enter mass of 2nd component',Msphere(2),0.)
       call prompt('enter scale length r_s for 2nd component,',rsoft(2),0.)
    case(9)
       call prompt('enter mass of central object',Mstar,0.)
       call prompt('enter radius of torus centre',Rtorus,0.)
       call prompt('enter distortion parameter ',distortion,1.,2.)
       if (abs(polyk-1.0).lt.tiny(polyk)) polyk = 0.0764
       call prompt('enter K in P= K*rho^gamma',polyk,0.)
    case(10)
       iexist = .false.
       do while(.not.iexist)
          call prompt('enter filename ',filename_exact)
          inquire(file=filename_exact,exist=iexist)
          if (iexist) then
             print "(a)",'file seems OK'
          else
             ians = .true.
             call prompt('file does not exist: try again? ',ians)
             if (.not.ians) return
          endif
       enddo
       call prompt('enter x axis of exact solution ',iexactplotx,1)
       call prompt('enter y axis of exact solution ',iexactploty,1)
       print "(a)",'apply column transformations to exact solution?'
       call prompt(' (no if file contains e.g. log y vs log x)',iApplyTransExactFile)
    end select

    return
  end subroutine submenu_exact
  
  !---------------------------------------------------
  ! sets options relating to exact solution plotting
  !---------------------------------------------------
  subroutine options_exact
    use prompting
    implicit none
    
    call prompt('enter number of exact solution points ',maxexactpts,10,1000000)
    call prompt('enter PGPLOT line colour ',iExactLineColour,1,16)
    call prompt('enter PGPLOT line style  ',iExactLineStyle,1,5)
    call prompt('calculate error norms? ',iCalculateExactErrors)
    if (iCalculateExactErrors) then
       call prompt('plot residuals (as inset in main plot)?',iPlotResiduals)
       if (iPlotResiduals) then
          call prompt('enter fraction of plot to use for inset', &
                      fracinsetResiduals,0.1,0.9)
          call prompt('enter max residual (0 for adaptive)',residualmax,0.)
       endif
    endif
  
    return
  end subroutine options_exact

  !-----------------------------------------------------------------------
  ! read exact solution parameters from files
  ! (in ndspmhd these files are used in the input to the code)
  !
  ! called after main data read and if exact solution chosen from menu
  !-----------------------------------------------------------------------
  subroutine read_exactparams(iexact,rootname,ierr)
    use settings_data, only:ndim
    use prompting
    implicit none
    integer, intent(in) :: iexact
    character(len=*), intent(in) :: rootname
    integer, intent(out) :: ierr
    
    integer :: idash
    character(len=len_trim(rootname)+8) :: filename

    idash = index(rootname,'_')
    if (idash.eq.0) idash = len_trim(rootname)+1

    select case(iexact)
    case(1)
       !
       !--shock tube parameters from .shk file
       !
       filename = trim(rootname(1:idash-1))//'.shk'
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
       select case(ndim)
       case(1)
          filename = trim(rootname(1:idash-1))//'.tstar'
          open(unit=20,ERR=8801,FILE=filename,STATUS='old')
          read(20,*,ERR=8888) Htstar,Ctstar,Atstar
          read(20,*,ERR=8888) sigma0
          read(20,*,ERR=8888) norder
          close(UNIT=20)
          print*,' >> read ',filename
          print*,' H,C,A,sigma,n = ',Htstar,Ctstar,Atstar,sigma0,norder
          return
8801      continue
          print*,'no file ',filename
          ierr = 1
          return
8888      print*,'error reading ',filename
          close(UNIT=20)
          ierr = 2
          return
       case(2)
          filename = trim(rootname(1:idash-1))//'.tstar2D'
          open(unit=20,ERR=9901,FILE=filename,STATUS='old')
          read(20,*,ERR=9902) Htstar,Ctstar,Atstar
          read(20,*,ERR=9902) alphatstar,betatstar,ctstar1,ctstar2
          read(20,*,ERR=9902) norder,morder
          close(UNIT=20)
          print*,' >> read ',filename
          print*,' j,m = ',norder,morder
          print*,' rho_0 = ',Htstar,' - ',Ctstar,' r^2' 
          if (norder.ge.0 .and. morder.ge.0) then
             print*,' v = ',Atstar,' r'
          else
             print*,' vx = ',alphatstar,'x +',ctstar1,'y'
             print*,' vy = ',ctstar2,'x +',betatstar,'y'
          endif
          return
9901      continue
          print*,'no file ',filename
          ierr = 1
          return
9902      print*,'error reading ',filename
          close(UNIT=20)
          ierr = 2
          return
       end select
    case(6)
       !
       !--attempt to guess which MHD shock tube has been done from filename
       !
          !read(rootname(5:5),*,iostat=ios) ishk
          !if (ios.ne.0) ishk = 1
       !
       !--prompt for shock type if not set  
       !
       if (ishk.le.0) then ! prompt
          ishk = 1
          call prompt('enter shock solution to plot',ishk,1,7)
       endif
       return

    end select

    return
  end subroutine read_exactparams

  !-----------------------------------------------------------------------
  ! this subroutine drives the exact solution plotting using the
  ! parameters which have been set
  !
  ! acts as an interface between the main plotting loop and the
  ! exact solution calculation subroutines
  !
  ! The exact solution is returned from the calculation via the arrays
  ! xexact and yexact. This means that the appropriate transformations
  ! can be applied (e.g. if the graph is logarithmic) and also ensures
  ! that the line style and colour settings are applied properly.
  !
  ! Note that we attempt to space the solution evenly in the transformed
  ! space (ie. in the current plot window), but this can be overwritten
  ! in the subroutines (for example if an uneven sampling is desired or
  ! the plotting is via some similarity variable as in the Sedov solution).
  ! In these cases the resulting arrays are then transformed, possibly leading
  ! to poor sampling in some regions (e.g. an evenly spaced array will become
  ! highly uneven in logarithmic space).
  !
  ! Note that any subroutine could in principle do its own plotting, 
  ! provided that it returns ierr > 0 which means that the generic line 
  ! is not plotted. Obviously transformations could not be applied in
  ! this case.
  !
  !-----------------------------------------------------------------------

  subroutine exact_solution(iexact,iplotx,iploty,itransx,itransy,igeom, &
                            ndim,ndimV,time,xmin,xmax,gamma,xplot,yplot, &
                            pmass,npart,imarker,unitsx,unitsy,irescale,iaxisy)
    use labels, only:ix,irad,iBfirst,ivx,irho,ike,iutherm,ih,ipr
    use prompting
    use exactfromfile, only:exact_fromfile
    use mhdshock, only:exact_mhdshock
    use polytrope, only:exact_polytrope
    use rhoh, only:exact_rhoh
    use sedov, only:exact_sedov
    use shock, only:exact_shock
    use torus, only:exact_torus
    use toystar1D, only:exact_toystar1D, exact_toystar_ACplane
    use toystar2D, only:exact_toystar2D
    use wave, only:exact_wave
    use densityprofiles, only:exact_densityprofiles
    use transforms
    implicit none
    integer, intent(in) :: iexact,iplotx,iploty,itransx,itransy,igeom
    integer, intent(in) :: ndim,ndimV,npart,imarker,iaxisy
    real, intent(in) :: time,xmin,xmax,gamma,unitsx,unitsy
    real, intent(in), dimension(npart) :: xplot,yplot,pmass
    logical, intent(in) :: irescale
    real, dimension(npart) :: residuals,ypart
    
    real, parameter :: zero = 1.e-10
    integer :: i,ierr,iexactpts,iCurrentColour,iCurrentLineStyle
    real, dimension(maxexactpts) :: xexact,yexact,xtemp
    real :: pmassmin,pmassmax,dx,ymean,errL1,errL2,errLinf

    !
    !--change line style and colour settings, but save old ones
    !
    call pgqci(iCurrentColour)
    call pgqls(iCurrentLineStyle)
    call pgsci(iExactLineColour)
    call pgsls(iExactLineStyle)

    !
    !--set x axis (can be overwritten)
    !  Need to space x in transformed space (e.g. in log space)
    !  but send the values of x in *real* space to the calculation routines
    !  then need to plot x in transformed space
    !
    !  Best solution is to set x grid initially, and inverse transform to get x values.
    !  These values can then be overwritten, if required in the exact subroutines
    !  We then re-transform the x array to plot it, which means that if spacing is
    !  overwritten the resulting array can still be transformed into log space
    !  but spacing will not be even
    !
    
    !--note that xmin and xmax will already have been transformed prior to input
    !  as these were the limits used for plotting the particles
    !
    dx = (xmax - xmin)/real(maxexactpts)
    do i=1,maxexactpts
       xexact(i) = xmin + (i-1)*dx
    enddo
    if (itransx.gt.0) call transform_inverse(xexact,itransx)
    
    iexactpts = maxexactpts
    !
    !--exact solution plots must return a zero or negative value of ierr to be plotted
    !  (-ve ierr indicates a partial solution)
    !
    ierr = 666

    select case(iexact)
    case(1)! shock tube
       if (iplotx.eq.ix(1) .and. igeom.le.1) then
          if (iploty.eq.irho) then
             call exact_shock(1,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
          elseif (iploty.eq.ipr) then
             call exact_shock(2,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
          elseif (iploty.eq.ivx) then
             call exact_shock(3,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
          elseif (iploty.eq.iutherm) then
             call exact_shock(4,time,gamma,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
          endif
       endif

    case(2)! sedov blast wave
       ! this subroutine does change xexact
       if (iplotx.eq.irad .or. (igeom.eq.3 .and. iplotx.eq.ix(1))) then
          if (iploty.eq.irho) then
             call exact_sedov(1,time,gamma,rhosedov,esedov,xmax,xexact,yexact,ierr)
          elseif (iploty.eq.ipr) then
             call exact_sedov(2,time,gamma,rhosedov,esedov,xmax,xexact,yexact,ierr)
          elseif (iploty.eq.iutherm) then
             call exact_sedov(3,time,gamma,rhosedov,esedov,xmax,xexact,yexact,ierr)
          elseif (iploty.eq.ike) then
             call exact_sedov(4,time,gamma,rhosedov,esedov,xmax,xexact,yexact,ierr)
          endif
       endif

    case(3)! polytrope
       if (iploty.eq.irho .and. (iplotx.eq.irad .or.(igeom.eq.3 .and. iplotx.eq.ix(1)))) then
          call exact_polytrope(gamma,polyk,xexact,yexact,iexactpts,ierr)
       endif

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
                call exact_toystar1D(1,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
             elseif (iploty.eq.ipr) then
                call exact_toystar1D(2,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)       
             elseif (iploty.eq.iutherm) then
                call exact_toystar1D(3,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)       
             elseif (iploty.eq.ivx) then
                call exact_toystar1D(4,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)       
             elseif (iploty.eq.iBfirst+1) then
                call exact_toystar1D(5,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
             endif
          elseif (iplotx.eq.irho) then
             if (iploty.eq.iBfirst+1) then
                call exact_toystar1D(6,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)       
             endif
          endif

          if (iploty.eq.iacplane) then! plot point on a-c plane
             call exact_toystar1D(7,time,gamma,htstar,atstar,ctstar,sigma,norder, &
                                xexact,yexact,iexactpts,ierr)
          endif
       else
          !
          !--2D toy star solutions
          !  these routines change xexact
          !
          if (igeom.eq.1 .and.((iplotx.eq.ix(1) .and. iploty.eq.ivx) &
               .or. (iplotx.eq.ix(2) .and. iploty.eq.ivx+1))) then
             call exact_toystar2D(4,time,gamma,polyk,totmass, &
                  atstar,htstar,ctstar,norder,morder, &
                  alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          endif
          if (iplotx.eq.irad .or. (igeom.eq.2 .and. iplotx.eq.ix(1))) then
             if (iploty.eq.irho) then
                call exact_toystar2D(1,time,gamma,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
             elseif (iploty.eq.ipr) then
                call exact_toystar2D(2,time,gamma,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
             elseif (iploty.eq.iutherm) then
                call exact_toystar2D(3,time,gamma,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
             elseif (igeom.eq.2 .and. iploty.eq.ivx) then
                call exact_toystar2D(4,time,gamma,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
             elseif (iploty.eq.ike) then
                call exact_toystar2D(5,time,gamma,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
             endif
          elseif (iplotx.le.ndim .and. iploty.le.ndim .and. igeom.eq.1) then
             call exact_toystar2D(0,time,gamma,polyk,totmass, &
                  atstar,htstar,ctstar,norder,morder, &
                  alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          endif
       endif

    case(5)! linear wave
       if ((iploty.eq.iwaveploty).and.(iplotx.eq.iwaveplotx)) then
          ymean = SUM(yplot(1:npart))/REAL(npart)
          call exact_wave(time,ampl,period,lambda,xzero,ymean,xexact,yexact,ierr)
       endif

    case(6) ! mhd shock tubes
       ! this subroutine modifies xexact
       if (iplotx.eq.ix(1) .and. igeom.le.1) then
          if (iploty.eq.irho) then
             call exact_mhdshock(1,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.ipr) then
             call exact_mhdshock(2,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.ivx) then
             call exact_mhdshock(3,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.ivx+1 .and. ndimV.gt.1) then
             call exact_mhdshock(4,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.ivx+ndimV-1 .and. ndimV.gt.2) then
             call exact_mhdshock(5,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.iBfirst+1 .and. ndimV.gt.1) then
             call exact_mhdshock(6,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.iBfirst+ndimV-1 .and. ndimV.gt.2) then
             call exact_mhdshock(7,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.iutherm) then
             call exact_mhdshock(8,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          elseif (iploty.eq.iBfirst) then
             call exact_mhdshock(9,ishk,time,gamma,xmin,xmax, &
                                 xexact,yexact,iexactpts,ierr)
          endif
       endif
    case(7) 
       !--h = (1/rho)^(1/ndim)
       if ((iploty.eq.ih).and.(iplotx.eq.irho)) then
          !--if variable particle masses, plot one for each pmass value
          pmassmin = minval(pmass)
          pmassmax = maxval(pmass)
          call exact_rhoh(ndim,hfact,pmassmin,xexact,yexact,ierr)

          if (abs(pmassmin-pmassmax).gt.zero .and. pmassmin.gt.zero) then
             !--plot first line
             if (ierr.le.0) then
                xtemp = xexact ! must not transform xexact as this is done again below
                if (itransx.gt.0) call transform(xtemp,itransx)
                if (itransy.gt.0) call transform(yexact,itransy)
                call pgline(iexactpts,xtemp(1:iexactpts),yexact(1:iexactpts))
             endif
             !--leave this one to be plotted below  
             call exact_rhoh(ndim,hfact,pmassmax,xexact,yexact,ierr)
          endif
       endif
    case(8) ! density profiles
       if (iplotx.eq.irad .or.(igeom.eq.3 .and. iplotx.eq.ix(1))) then
          if (iploty.eq.irho) then
             call exact_densityprofiles(1,iprofile,Msphere,rsoft,xexact,yexact,ierr)
          elseif (iploty.eq.icolpoten) then
             call exact_densityprofiles(2,iprofile,Msphere,rsoft,xexact,yexact,ierr)
          elseif (iploty.eq.icolfgrav) then
             call exact_densityprofiles(3,iprofile,Msphere,rsoft,xexact,yexact,ierr)          
          endif
       endif
    case(9) ! torus
       if (iplotx.eq.irad .or.(igeom.eq.3 .and. iplotx.eq.ix(1))) then
          if (iploty.eq.irho) then
             call exact_torus(1,Mstar,Rtorus,polyk,distortion,gamma,xexact,yexact,ierr)
          elseif (iploty.eq.ipr) then
             call exact_torus(2,Mstar,Rtorus,polyk,distortion,gamma,xexact,yexact,ierr)      
          elseif (iploty.eq.iutherm) then
             call exact_torus(3,Mstar,Rtorus,polyk,distortion,gamma,xexact,yexact,ierr)      
          endif
       !--pr vs z at r=Rtorus
       elseif (igeom.eq.2 .and. iplotx.eq.ix(3) .and.iploty.eq.ipr) then
          call exact_torus(4,Mstar,Rtorus,polyk,distortion,gamma,xexact,yexact,ierr)      
       endif
    case(10) ! exact solution read from file
       if (iplotx.eq.iexactplotx .and. iploty.eq.iexactploty) then   
          call exact_fromfile(filename_exact,xexact,yexact,iexactpts,ierr)
          !--plot this untransformed (as may already be in log space)
          if (ierr.le.0 .and. .not.iApplyTransExactFile) then
             call pgline(iexactpts,xexact(1:iexactpts),yexact(1:iexactpts))
             ierr = 1
          endif
          !--change into physical units if appropriate
          if (iRescale .and. iApplyTransExactFile) then
             xexact(1:iexactpts) = xexact(1:iexactpts)*unitsx
             yexact(1:iexactpts) = yexact(1:iexactpts)*unitsy
          endif
       endif
    end select
    
    !----------------------------------------------------------
    !  plot this as a line on the current graph using PGPLOT
    !----------------------------------------------------------
    if (ierr.le.0) then
       if (itransx.gt.0) call transform(xexact(1:iexactpts),itransx)
       if (itransy.gt.0) call transform(yexact(1:iexactpts),itransy)
       call pgline(iexactpts,xexact(1:iexactpts),yexact(1:iexactpts))
       !
       !--calculate errors
       !
       if (iCalculateExactErrors) then
          !--untransform y axis again for error calculation
          if (itransy.gt.0) call transform_inverse(yexact(1:iexactpts),itransy)
          !--untransform particle y axis also
          ypart(1:npart) = yplot(1:npart)
          if (itransy.gt.0) call transform_inverse(ypart(1:npart),itransy)          
          !--calculate errors
          call calculate_errors(xexact(1:iexactpts),yexact(1:iexactpts), &
                                xplot(1:npart),ypart(1:npart),residuals(1:npart), &
                                errL1,errL2,errLinf)
          print "(3(a,1pe10.3,1x))",' L1 error = ',errL1,' L2 error = ',errL2, &
                                   ' L(infinity) error = ',errLinf
          if (iPlotResiduals) call plot_residuals(xplot,residuals,imarker,iaxisy)
       endif
    endif
    !
    !--reset line and colour settings
    !   
    call pgsci(iCurrentColour)
    call pgsls(iCurrentLineStyle)

    return

  end subroutine exact_solution

  subroutine calculate_errors(xexact,yexact,xpts,ypts,residual,errL1,errL2,errLinf)
   implicit none
   real, dimension(:), intent(in) :: xexact,yexact,xpts,ypts
   real, dimension(size(xpts)), intent(out) :: residual
   real, intent(out) :: errL1,errL2,errLinf
   integer :: i,j,npart,iused
   real :: xi,dy,dx,yexacti,err1,ymax

   errL1 = 0.
   errL2 = 0.
   errLinf = 0.
   residual = 0.
   npart = size(xpts)
   iused = 0
   ymax = -huge(ymax)
   
   do i=1,npart
      xi = xpts(i)
      yexacti = 0.
      !
      !--find nearest point in exact solution table
      !
      do j=1,size(xexact)-1
         if (xexact(j).lt.xi .and. xexact(j+1).gt.xi) then
            if (abs(residual(i)).gt.tiny(residual)) print*,'already used ',i
            !--linear interpolation from tabulated exact solution
            dy = yexact(j+1) - yexact(j)
            dx = xexact(j+1) - xexact(j)
            if (dx.gt.0.) then
               yexacti = yexact(j) + dy/dx*(xi - xexact(j))
               residual(i) = ypts(i) - yexacti
            elseif (dy.gt.0.) then
               yexacti = yexact(j)
               residual(i) = ypts(i) - yexacti
            else
               print "(a)",'error in residual calculation'
               residual(i) = 0.
            endif
            iused = iused + 1
            ymax = max(ymax,abs(yexacti))
         endif
      enddo
      err1 = abs(residual(i))
      errL1 = errL1 + err1
      errL2 = errL2 + err1**2
      errLinf = max(errLinf,err1)
      if (yexacti.gt.tiny(yexacti)) residual(i) = residual(i)/abs(yexacti)
   enddo
   !
   !--normalise errors (use maximum y value)
   !
   if (ymax.gt.tiny(ymax)) then
      errL1 = errL1/(npart*ymax)
      errL2 = sqrt(errL2/(npart*ymax**2))
      errLinf = errLinf/ymax
   else
      print "(a)",'error normalising errors'
      errL1 = 0.
      errL2 = 0.
      errLinf = 0.
   endif
   
   if (iused.ne.npart) print*,'errors calculated using ',iused,' of ',npart, 'particles'
   
   return
  end subroutine calculate_errors
  
  subroutine plot_residuals(xpts,residuals,imarker,iaxisy)
   implicit none
   real, dimension(:), intent(in) :: xpts,residuals
   integer, intent(in) :: imarker,iaxisy
   real :: vptxminold,vptxmaxold,vptyminold,vptymaxold
   real :: vptxmin,vptxmax,vptymin,vptymax
   real :: xminold,xmaxold,yminold,ymaxold,ymin,ymax
   real :: xch,ych
   integer :: ioldcolour,ioldfill

   !--query old viewport and window size
   call pgqvp(0,vptxminold,vptxmaxold,vptyminold,vptymaxold)
   call pgqwin(xminold,xmaxold,yminold,ymaxold)

   !--use specified bottom % of viewport
   vptxmin = vptxminold
   vptxmax = vptxmaxold
   vptymin = vptyminold
   vptymax = vptyminold + FracinsetResiduals*(vptymaxold - vptyminold)
   call pgsvp(vptxmin,vptxmax,vptymin,vptymax)
 
   !--set window
   if (residualmax.lt.tiny(residualmax)) then
      ymax = maxval(abs(residuals))
      print*,'max residual = ',ymax
   else
      ymax = residualmax
   endif
   ymin = -ymax

   !--erase space for residual plot
   call pgqci(ioldcolour)
   call pgqfs(ioldfill)
   call pgqcs(0,xch,ych)
   call pgsci(0)
   call pgsfs(1)
   if (iaxisy.lt.0) then
      call pgsvp(vptxmin,vptxmax,vptymin,vptymax)   
   else
      call pgsvp(vptxmin - 3.*xch,vptxmax,vptymin,vptymax)
   endif
   call pgswin(xminold,xmaxold,ymin,ymax)
   call pgrect(xminold,xmaxold,ymin,ymax)
   call pgsci(ioldcolour)
   call pgsfs(ioldfill)
   !--set window and draw axes
   call pgsvp(vptxmin,vptxmax,vptymin,vptymax)
   call pgswin(xminold,xmaxold,ymin,ymax)
   if (iaxisy.lt.0) then
      call pgbox('ABCST',0.0,0,'BCST',0.0,0)   
   else
      call pgbox('ABCST',0.0,0,'BVNCST',0.0,0)
   endif
   
   !--plot residuals
   call pgpt(size(xpts),xpts,residuals,imarker)
   
   !--restore old viewport and window
   call pgsvp(vptxminold,vptxmaxold,vptyminold,vptymaxold)
   call pgswin(xminold,xmaxold,yminold,ymaxold)
   
  end subroutine plot_residuals

end module exact
