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
  logical :: iApplyTransExactFile
  !
  !--declare all of the parameters required for the various exact solutions
  !
  !--toy star
  integer :: iACplane ! label position of toy star AC plane plot
  integer :: norder,morder ! for toy star
  real :: htstar,atstar,ctstar,sigma,alphatstar,betatstar,ctstar1,ctstar2
  real :: sigma0
  !--sound wave
  integer :: iwaveploty,iwaveplotx ! linear wave
  real :: ampl,lambda,period,xzero
  !--sedov blast wave
  real :: rhosedov,esedov
  !--polytrope
  real :: polyk
  !--mhd shock solutions
  integer :: ishk
  !--from file
  integer :: iexactplotx, iexactploty
  !--shock tube
  real :: rho_L, rho_R, pr_L, pr_R, v_L, v_R
  !--rho vs h
  real :: hfact
  character(len=120) :: filename_exact
  !
  !--sort these into a namelist for input/output
  !
  namelist /exactopts/ iexactplotx,iexactploty,filename_exact,maxexactpts, &
       iExactLineColour,iExactLineStyle,iApplyTransExactFile

  namelist /exactparams/ ampl,lambda,period,iwaveploty,iwaveplotx,xzero, &
       htstar,atstar,ctstar,alphatstar,betatstar,ctstar1,ctstar2, &
       polyk,sigma0,norder,morder,rhosedov,esedov, &
       rho_L, rho_R, pr_L, pr_R, v_L, v_R, hfact

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
    norder = -1
    morder = 0
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
    iexactplotx = 0
    iexactploty = 0
    ishk = 0
    hfact = 1.2
    filename_exact = ' '

    maxexactpts = 1001      ! points in exact solution plot
    iExactLineColour = 1    ! foreground
    iExactLineStyle = 1     ! solid
    iApplyTransExactFile = .true. ! false if exact from file is already logged
    
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
    logical :: ians

    print 10
10  format(' 0) none ',/,               &
         ' 1) shock tube ',/,           &
         ' 2) sedov blast wave ',/,     &
         ' 3) polytrope ',/,            &
         ' 4) toy star ',/,             &
         ' 5) linear wave ',/,          &
         ' 6) mhd shock tubes (tabulated) ',/,  &
         ' 7) h vs rho ',/, &
         ' 8) read from file ')
    call prompt('enter exact solution to plot',iexact,0,8)
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
       print*,' MHD shock tube tables: '
       call prompt('enter solution to plot ',ishk,0,7)
    case(7)
       call prompt('enter hfact [h = hfact*(m/rho)**1/ndim]',hfact,0.)
    case(8)
       call prompt('enter filename ',filename_exact)
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
    implicit none
    integer, intent(in) :: iexact
    character(len=*), intent(in) :: rootname
    integer, intent(out) :: ierr
    
    integer :: ios,idash
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
       read(rootname(5:5),*,iostat=ios) ishk
       if (ios.ne.0) ishk = 1
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
                            ndim,ndimV,time,xmin,xmax,ymean,gamma,pmass,npart)
    use labels, only:ix,irad,iBfirst,ivx,irho,ike,iutherm,ih,ipr
    use prompting
    use exactfromfile, only:exact_fromfile
    use mhdshock, only:exact_mhdshock
    use polytrope, only:exact_polytrope
    use rhoh, only:exact_rhoh
    use sedov, only:exact_sedov
    use shock, only:exact_shock
    use toystar1D, only:exact_toystar1D, exact_toystar_ACplane
    use toystar2D, only:exact_toystar2D
    use wave, only:exact_wave
    use transforms
    implicit none
    integer, intent(in) :: iexact,iplotx,iploty,itransx,itransy,igeom
    integer, intent(in) :: ndim,ndimV,npart
    real, intent(in) :: time,xmin,xmax,ymean,gamma
    real, intent(in), dimension(npart) :: pmass
    
    real, parameter :: zero = 1.e-10
    integer :: i,ierr,iexactpts,iCurrentColour,iCurrentLineStyle
    real, dimension(maxexactpts) :: xexact,yexact,xtemp
    real :: totmass,pmassmin,pmassmax,dx

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
       if (iplotx.eq.irad .or. (igeom.eq.2 .and. iplotx.eq.ix(1))) then
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
       if (iploty.eq.irho .and. iplotx.eq.irad) then
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
          totmass = SUM(pmass(1:npart))
          print*,'summing masses of ',npart,' particles, mass = ',totmass
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
          call exact_wave(time,ampl,period,lambda,xzero,ymean,xexact,yexact,ierr)
       endif

    case(6) ! mhd shock tubes
       ! this subroutine modifies xexact
       if (iplotx.eq.ix(1) .and. igeom.le.1) then
          !--prompt for shock type if not set  
          if (ishk.eq.0) then ! prompt
             call prompt('enter shock solution to plot',ishk,0,6)
          endif
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
    case(8) ! exact solution read from file
       if (iplotx.eq.iexactplotx .and. iploty.eq.iexactploty) then   
          call exact_fromfile(filename_exact,xexact,yexact,iexactpts,ierr)
          !--plot this untransformed (as may already be in log space)
          if (ierr.le.0 .and. .not.iApplyTransExactFile) then
             call pgline(iexactpts,xexact(1:iexactpts),yexact(1:iexactpts))
             ierr = 1
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
    endif
    !
    !--reset line and colour settings
    !   
    call pgsci(iCurrentColour)
    call pgsls(iCurrentLineStyle)

    return

  end subroutine exact_solution

end module exact
