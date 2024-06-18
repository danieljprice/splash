!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

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
 use params, only:maxplot
 implicit none
 !
 !--maximum number of solutions in any one plot
 !
 integer, parameter :: maxexact=32
 integer, parameter :: maxlabels = maxplot
 character(len=30)  :: exact_labels(maxlabels)
 !
 !--options used to plot the exact solution line
 !
 integer :: maxexactpts, iExactLineColour(maxexact), iExactLineStyle(maxexact),iPlotExactOnlyOnPanel
 integer :: iNormaliseErrors
 logical :: iApplyTransExactFile,iCalculateExactErrors,iPlotResiduals
 logical :: iApplyUnitsExactFile,iPlotExactUnder
 real :: fracinsetResiduals,residualmax,ExactAlpha
 character(len=60) :: ExactLegendText(maxexact)
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
 real :: xshock
 !--density profiles
 integer :: iprofile,icolpoten,icolfgrav
 real, dimension(2) :: Msphere,rsoft
 !--from file
 integer :: iexactplotx(maxexact), iexactploty(maxexact)
 logical :: iauto_map_columns
 !--shock tube
 real :: rho_L, rho_R, pr_L, pr_R, v_L, v_R
 !--rho vs h
 real :: hfact
 !--read from file
 integer :: imapexact(maxlabels),nfiles
 character(len=120) :: filename_exact(maxexact)
 !--equilibrium torus
 real :: Mstar,Rtorus,distortion
 !--ring spreading
 real :: Mring,Rring,viscnu
 !--dusty waves
 real :: cs,Kdrag,rhozero,rdust_to_gas
 !--arbitrary function
 integer :: nfunc
 character(len=120), dimension(maxexact) :: funcstring
 !--Roche potential
 real :: mprim,msec
 real :: xprim(3),xsec(3)
 logical :: use_sink_data
 integer, parameter :: iexact_rochelobe = 15
 !--C-shock
 real :: machs,macha
 !--planet-disc interaction
 real :: HonR,rplanet,q_index,phase
 real :: spiral_params(7,maxexact)
 integer :: ispiral,narms
 logical :: iread_wakeparams,iclockwise
 !--bondi flow
 logical :: relativistic, geodesic_flow,is_wind
 real    :: const1,const2
 !--gamma, for manual setting
 real :: gamma_exact
 logical :: use_gamma_exact
 !
 !--sort these into a namelist for input/output
 !
 namelist /exactopts/ iexactplotx,iexactploty,filename_exact,maxexactpts, &
       iExactLineColour,iExactLineStyle,iApplyTransExactFile,iCalculateExactErrors, &
       iPlotResiduals,fracinsetResiduals,residualmax,iPlotExactOnlyOnPanel,&
       iApplyUnitsExactFile,iNormaliseErrors,iPlotExactUnder,ExactAlpha

 namelist /exactparams/ ampl,lambda,period,iwaveploty,iwaveplotx,xzero, &
       htstar,atstar,ctstar,alphatstar,betatstar,ctstar1,ctstar2, &
       polyk,sigma0,norder,morder,rhosedov,esedov, &
       rho_L, rho_R, pr_L, pr_R, v_L, v_R,ishk,hfact, &
       iprofile,Msphere,rsoft,icolpoten,icolfgrav,Mstar,Rtorus,distortion, &
       Mring,Rring,viscnu,nfunc,funcstring,cs,Kdrag,rhozero,rdust_to_gas, &
       mprim,msec,imapexact,iauto_map_columns,xshock,totmass,machs,macha,&
       use_sink_data,xprim,xsec,nfiles,gamma_exact,use_gamma_exact,&
       HonR,rplanet,q_index,relativistic,geodesic_flow,is_wind,&
       const1,const2,ispiral,narms,spiral_params,phase,iread_wakeparams,iclockwise

 public :: defaults_set_exact,submenu_exact,options_exact,read_exactparams
 public :: exact_solution,get_nexact
 public :: exactopts,exactparams

contains
 !----------------------------------------------------------------------
 ! sets default values of the exact solution parameters
 !----------------------------------------------------------------------
subroutine defaults_set_exact
 use plotlib, only:plotlib_maxlinestyle
 integer :: i

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
 ishk = 1
 xshock = 0.
 hfact = 1.2
!   read from file
 nfiles = 1
 filename_exact = ' '
 imapexact = 0
 iauto_map_columns = .true.
 iexactplotx = 0
 iexactploty = 0
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
!   ring spreading
 Mring = 1.0
 Rring = 1.0
 viscnu = 1.e-3
!   dusty waves
 Kdrag = 1.0
 cs    = 1.0
 rhozero = 1.0
 rdust_to_gas = 0.0
!   Roche lobes
 use_sink_data = .true.
 mprim = 1.
 msec  = 1.
 xprim = 0.
 xsec  = 0.
 xsec(1) = 1.
!   C-shock
 machs = 50. ! sonic Mach number
 macha = 5.  ! Alfvenic Mach number
!   planet-disc interaction
 HonR = 0.05
 rplanet = 1.
 q_index = 0.25
 ispiral = 1
 narms = 1
 phase = 0.
 spiral_params = 0.
 spiral_params(2,:) = 360.
 iread_wakeparams = .true.
 iclockwise = .false.
!   Bondi
 relativistic  = .true.
 geodesic_flow = .false.
 is_wind       = .true.
 const1 = 8.
 const2 = 1.

!   gamma, if not read from file
 gamma_exact = 5./3.
 use_gamma_exact = .false.

!   arbitrary function
 nfunc = 1
 funcstring = ' '

 maxexactpts = 1001      ! points in exact solution plot
 iExactLineColour = 1    ! foreground
 do i=1,size(iExactLineStyle)
    iExactLineStyle(i) = mod(i,plotlib_maxlinestyle)
 enddo
 iApplyTransExactFile = .true. ! false if exact from file is already logged
 iApplyUnitsExactFile = .false.
 iCalculateExactErrors = .true.
 iNormaliseErrors = 1
 iPlotResiduals = .false.
 fracinsetResiduals = 0.15
 residualmax = 0.0
 iPlotExactOnlyOnPanel = 0
 iPlotExactUnder = .false.
 ExactLegendText = ''
 ExactAlpha = 1.0

end subroutine defaults_set_exact

 !----------------------------------------------------------------------
 ! sets which exact solution to calculate + parameters for this
 !----------------------------------------------------------------------
subroutine submenu_exact(iexact)
 use settings_data, only:ndim,ncolumns,ncalc
 use prompting,     only:prompt,print_logical
 use filenames,     only:rootname,ifileopen,fileprefix
 use exactfunction, only:check_function
 use mhdshock,      only:nmhdshocksolns,mhdprob
 use planetdisc,    only:maxspirals,labelspiral
 use asciiutils,    only:get_ncolumns,get_nrows,string_replace,match_tag_start,&
                         add_escape_chars,read_asciifile,read_column_labels
 use labels,        only:label
 use map_columns,   only:map_columns_in_file,map_columns_interactive,print_mapping
 integer, intent(inout) :: iexact
 integer :: ierr,itry,i,ncols,nadjust,nrows,nlab_exact
 integer, allocatable :: imapauto(:)
 logical :: ians,iexist,ltmp,prompt_for_gamma,apply_to_all
 character(len=len(filename_exact)) :: filename_tmp
 character(len=4) :: str
 character(len=len(funcstring)) :: func_prev

 print 10
10 format(' 0) none ',/,               &
           ' 1) ANY function f(x,t)',/, &
           ' 2) read from file ',/,       &
           ' 3) shock tube ',/,           &
           ' 4) sedov blast wave ',/,     &
           ' 5) polytrope ',/,            &
           ' 6) toy star ',/,             &
           ' 7) gresho vortex ',/,          &
           ' 8) mhd shock tubes (tabulated) ',/,  &
           ' 9) h vs rho ',/, &
           '10) Plummer/Hernquist spheres ',/, &
           '11) torus ',/, &
           '12) ring spreading ',/, &
           '13) special relativistic shock tube', /, &
           '14) dusty waves', /, &
           '15) Roche lobes/potential ',/, &
           '16) C-shock ',/, &
           '17) Spiral arms (planet-disc interaction)',/, &
           '18) Bondi flow')
 call prompt('enter exact solution to plot',iexact,0,18)
 print "(a,i2)",'plotting exact solution number ',iexact
 !
 !--enter parameters for various exact solutions
 !
 prompt_for_gamma = .false.
 apply_to_all = .true.
 select case(iexact)
 case(1)
    call prompt('enter number of functions to plot ',nfunc,1,maxexact)
    print "(/,a,6(/,11x,a))",' Examples: sin(2*pi*x - 0.1*t)','sqrt(0.5*x)','x^2', &
             'exp(-2*x**2 + 0.1*t)','log10(x/2)','exp(y),y=sin(pi*x)','cos(z/y),z=acos(y),y=x^2'
    ierr = 0
    overfunc: do i=1,nfunc
       ierr = 1
       itry = 0
       func_prev = funcstring(i)
       do while(ierr /= 0 .and. itry < 10)
          if (nfunc > 1) print "(/,a,i2,/,11('-'),/)",'Function ',i
          call prompt('enter function f(x,t) to plot ',funcstring(i),noblank=.true.)
          call check_function(funcstring(i),ierr)
          if (ierr /= 0 .and. len(funcstring(i))==len_trim(funcstring(i))) then
             print "(a,i3,a)",&
                     ' (errors are probably because string is too long, max length = ',&
                     len(funcstring(i)),')'
          endif
          itry = itry + 1
       enddo
       if (itry >= 10) then
          print "(a)",' *** too many tries, aborting ***'
          ierr = i-1
          exit overfunc
       endif
       print*
       call prompt('enter y axis of exact solution (0=all plots)',iexactploty(i),0)
       if (iexactploty(i) > 0) then
          call prompt('enter x axis of exact solution ',iexactplotx(i),1)
       endif
       if (len_trim(ExactLegendText(i))==0 .or. trim(funcstring(i)) /= trim(func_prev)) then
          ExactLegendText(i) = trim(funcstring(i))
       endif
       call prompt('enter text to display in legend (blank=do not show)',ExactLegendText(i))
    enddo overfunc
    if (ierr /= 0) nfunc = ierr
 case(2)
    !
    ! try to read filenames from .exactfiles if it exists
    !
    call read_exactparams(iexact,trim(rootname(1)),ierr)
    !
    ! then prompt user
    !
    call prompt('enter number of files to read per plot ',nfiles,1,maxexact)
    if (ierr /= 0 .and. nfiles > 1) then
       print "(2(/,a))",' Hint: create a file called '//trim(fileprefix)//'.exactfiles', &
                                     ' with one filename per line to automatically read the filenames'
    else
       print "(a,i2,a)",' got ',nfiles,' filenames from '//trim(fileprefix)//'.exactfiles'
    endif

    nadjust = -1
    over_files: do i=1,nfiles
       iexist = .false.
       do while(.not.iexist)
          if (i==1) print "(/,a)",'Use %f to represent current dump file, e.g. %f.exact looks for dump_000.exact'
          write(str,"(i4)") i
          call prompt('enter filename #'//trim(adjustl(str)),filename_exact(i))
          !--substitute %f for filename
          filename_tmp = filename_exact(i)
          call string_replace(filename_tmp,'%f',trim(rootname(ifileopen)))
          !--check the first file for errors
          inquire(file=filename_tmp,exist=iexist)
          if (iexist) then
             allocate(imapauto(size(imapexact)))
             imapauto = 0
             call map_columns_in_file(filename_tmp,ncols,nrows,imapauto,&
                                      label(1:ncolumns+ncalc),exact_labels,nlab_exact,ierr)
             if (ierr==0) then
                ! if more than two columns mapped successfully, prompt to keep this automatic
                if (all(imapexact==0)) imapexact = imapauto
                ! manual assignment of column mappings
                call map_columns_interactive(imapexact,label(1:ncolumns+ncalc),&
                                             exact_labels,nlab_exact)
                ! disable automated mapping if edits were made
                if (any(imapexact(1:nlab_exact) /= imapauto(1:nlab_exact))) iauto_map_columns = .false.
                print "(a)",' automated mapping is '//print_logical(iauto_map_columns)

                if (nrows > maxexactpts) maxexactpts = nrows
                if (count(imapexact > 0) < 2) print "(/,a,/)",' WARNING: not enough columns mapped to plot anything'
                ! re-prompt if there are less than two useable columns in the file
                if (ncols < 2) then
                   iexist = .false.
                   call prompt('Error: file contains < 2 readable columns: try again?',ians)
                   if (.not.ians) then
                      nadjust = i-1
                      exit over_files
                   endif
                endif
             else
                iexist = .false.
                call prompt('Error opening '//trim(filename_tmp)//': try again?',ians)
                if (.not.ians) then
                   nadjust = i-1
                   exit over_files
                endif
             endif
             if (allocated(imapauto)) deallocate(imapauto)
          else
             ians = .true.
             call prompt('file does not exist: try again? ',ians)
             if (.not.ians) then
                nadjust = i-1
                exit over_files
             endif
          endif
       enddo
       !if (i==1) then
       !    call prompt('Apply above settings to all files?',apply_to_all)
       !    if (apply_to_all) then
       !       iexactploty(:) = iexactploty(i)
       !       iexactplotx(:) = iexactplotx(i)
       !    endif
       !endif

       if (len_trim(ExactLegendText(i))==0) ExactLegendText(i) = add_escape_chars(filename_exact(i))
       !call prompt('enter text to display in legend (blank=do not show)',ExactLegendText(i))
    enddo over_files
    if (nadjust >= 0) then
       nfiles = nadjust
       if (nfiles == 1) then
          print "(/,a)",'Using 1 file only'
       elseif (nfiles > 0) then
          write(str,"(i4)") nfiles
          print "(/,a)",'Using '//trim(adjustl(str))//' files'
       endif
    endif
    if (nfiles > 0) then ! only ask if filename was read OK
       ltmp = .not.iApplyTransExactFile
       if (.not.any(exact_labels(:)(1:3) == 'log')) then
          call prompt(' are exact solutions already logged?',ltmp)
       endif
       iApplyTransExactFile = .not.ltmp
       ltmp = .not.iApplyUnitsExactFile
       call prompt(' are exact solutions in physical units?',ltmp)
       iApplyUnitsExactFile = .not.ltmp
    endif
 case(3,13)
    !
    !--read shock parameters from the .shk file
    !
    call read_exactparams(iexact,trim(rootname(1)),ierr)
    if (ierr /= 0) then
       call prompt('enter initial x position of shock',xshock)
       call prompt('enter density to left of shock   ',rho_L,0.0)
       call prompt('enter density to right of shock  ',rho_R,0.0)
       call prompt('enter pressure to left of shock  ',pr_L,0.0)
       call prompt('enter pressure to right of shock ',pr_R,0.0)
       if (iexact==13) then
          call prompt('enter velocity to left of shock  ',v_L,max=1.0)
          call prompt('enter velocity to right of shock ',v_R,max=1.0)
       else
          call prompt('enter velocity to left of shock  ',v_L)
          call prompt('enter velocity to right of shock ',v_R)
          call prompt('enter dust-to-gas ratio ',rdust_to_gas,0.)
       endif
    endif
    prompt_for_gamma = .true.
 case(4)
    call prompt('enter density of ambient medium ',rhosedov,0.0)
    call prompt('enter blast wave energy E ',esedov,0.0)
    prompt_for_gamma = .true.
 case(5)
    call prompt('enter polytropic k ',polyk)
    call prompt('enter total mass ',totmass)
    prompt_for_gamma = .true.
 case(6)
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
       if (ndim >= 2) call prompt('enter order of angular mode',morder,0)
       call prompt('enter velocity amplitude a (v = a*r) ',atstar)
    else
       print "(a)",'using exact non-linear solution:'
       ians = .true.
       if (norder < 0 .and. morder < 0) ians = .false.
       if (ndim >= 2) call prompt('axisymmetric?',ians)
       if (ians .or. ndim==1) then
          norder = -1
          morder = 0
          call prompt('enter v_r amplitude ',alphatstar)
          if (ndim >= 2) call prompt('enter v_phi amplitude ',betatstar)
       else
          norder = -1
          morder = -1
          call prompt('enter vxx amplitude ',alphatstar)
          call prompt('enter vyy amplitude ',betatstar)
          call prompt('enter vxy amplitude ',ctstar1)
          call prompt('enter vyx amplitude ',ctstar2)
       endif
    endif
    prompt_for_gamma = .true.
    !case(7)
    !   call prompt('enter y-plot to place sine wave on',iwaveploty,1)
    !   call prompt('enter x-plot to place sine wave on',iwaveplotx,1)
    !   call prompt('enter starting x position',xzero)
    !   call prompt('enter wavelength lambda ',lambda,0.0)
    !   call prompt('enter amplitude ',ampl,0.0)
    !   call prompt('enter period ',period)
 case(8)
    print "(a)",' MHD shock tube tables: '
    if (ishk <= 0) ishk = 1
    do i=1,nmhdshocksolns
       print "(i2,') ',a)",i,trim(mhdprob(i))
    enddo
    call prompt('enter solution to plot ',ishk,1,7)
    call prompt('enter initial x position of shock ',xshock)
 case(9)
    call prompt('enter hfact [h = hfact*(m/rho)**1/ndim]',hfact,0.)
 case(10)
    print 20
20  format(' 1) Plummer sphere  [ rho = 3M r_s**2 /(4 pi (r**2 + r_s**2)**5/2) ]',/, &
              ' 2) Hernquist model [ rho =     M r_s /(2 pi r (r_s + r)**3        ]')
    call prompt('enter density profile to plot',iprofile,1,2)
    call prompt('enter total mass of sphere M',Msphere(1),0.)
    call prompt('enter scale length length r_s,',rsoft(1),0.)
    ians = .false.
    if (icolpoten > 0) ians = .true.
    call prompt('Are the gravitational potential and/or force dumped?',ians)
    if (ians) then
       call prompt('enter column containing grav. potential',icolpoten,0)
       call prompt('enter column containing grav. force',icolfgrav,0)
    endif
    call prompt('enter mass of 2nd component',Msphere(2),0.)
    call prompt('enter scale length r_s for 2nd component,',rsoft(2),0.)
 case(11)
    call prompt('enter mass of central object',Mstar,0.)
    call prompt('enter radius of torus centre',Rtorus,0.)
    call prompt('enter distortion parameter ',distortion,1.,2.)
    if (abs(polyk-1.0) < tiny(polyk)) polyk = 0.0764
    call prompt('enter K in P= K*rho^gamma',polyk,0.)
    prompt_for_gamma = .true.
 case(12)
    call prompt('enter mass of ring',Mring,0.)
    call prompt('enter radius of ring centre R0',Rring,0.)
    call prompt('enter viscosity parameter nu',viscnu,0.)
 case(14)
    call prompt('enter starting x position',xzero)
    call prompt('enter wavelength lambda ',lambda,0.)
    call prompt('enter amplitude of perturbation',ampl,0.)
    call prompt('enter sound speed in gas ',cs,0.)
    call prompt('enter initial gas density ',rhozero,0.)
    call prompt('enter dust-to-gas ratio ',rdust_to_gas,0.)
    call prompt('enter drag coefficient K ',Kdrag,0.)
 case(15)
    call prompt('use data from sink particles?',use_sink_data)
    if (.not.use_sink_data) then
       call prompt('enter mass of primary star ',mprim,0.)
       call prompt('enter mass of secondary star ',msec,0.,mprim)
       call prompt('enter x position of primary',xprim(1))
       call prompt('enter y position of primary',xprim(2))
       call prompt('enter x position of secondary',xsec(1))
       call prompt('enter y position of secondary',xsec(2))
    endif
 case(16)
    call prompt('enter sonic Mach number',machs,0.)
    call prompt('enter Alfvenic Mach number ',macha,0.)
 case(17)
    do i=1,maxspirals
       print "(1x,i1,')',1x,a)",i,trim(labelspiral(i))
    enddo
    call prompt('Which spiral arm solution to plot?',ispiral,0,maxspirals)
    select case(ispiral)
    case(2)
       call prompt('enter number of spiral arms to plot',narms,1,maxexact)
       print "(3(/,a),/)",' Spiral data can be read from '//trim(fileprefix)//'.spirals file', &
                             ' with one line per arm and columns corresponding to:', &
                             ' PAmin, PAmax, a0, a1, a2, a3, a4'
       do i=1,narms
          print "(a,i1,a)",' -- arm ',i,' --'
          call prompt('enter starting position angle (deg)',spiral_params(1,i))
          call prompt('enter finishing position angle (deg)',spiral_params(2,i))
          call prompt('enter a0 in r = \sum a_i phi^i',spiral_params(3,i))
          call prompt('enter a1 in r = \sum a_i phi^i',spiral_params(4,i))
          call prompt('enter a2 in r = \sum a_i phi^i',spiral_params(5,i))
          call prompt('enter a3 in r = \sum a_i phi^i',spiral_params(6,i))
          call prompt('enter a4 in r = \sum a_i phi^i',spiral_params(7,i))
       enddo
    case default
       if (ifileopen > 0) call get_planetdisc_parameters_from_data(rplanet,HonR,q_index,phase,ierr)
       if (ierr == 0) call prompt(' use wake parameters above?',iread_wakeparams)
       if (ierr /= 0) iread_wakeparams = .false.
       if (.not.iread_wakeparams) then
          call prompt('enter disc aspect ratio at planet location (H/R)',HonR,0.,1.)
          call prompt('enter planet orbital radius ',rplanet,0.)
          call prompt('enter power-law index of sound speed cs ~ R^-q',q_index)
          call prompt('enter orbital phase in degrees ',phase)
          call prompt('clockwise rotation?',iclockwise)
       endif
    end select
 case(18)
    prompt_for_gamma = .true.
    call prompt('is it a wind (instead of accretion)?',is_wind)
    call prompt('enter mass of central object',Mstar,0.)
    call prompt('do you want a relativistic solution?',relativistic)
    if (.not.relativistic) then
       prompt_for_gamma = .false.
       call prompt('enter Parker/Bondi critical radius (rcrit)',const1)
       call prompt('enter density at critical radius (rhocrit)',const2)
    elseif (relativistic) then
       call prompt('is it the geodesic flow?',geodesic_flow)
       if (geodesic_flow) then
          const1 = 1.
          const2 = 1.e-9
          call prompt('enter constant den0',const1,min=0.)
          call prompt('enter constant en0 ',const2,min=0.)
       else
          call prompt('enter the critical radius in units of the central mass M',const1,min=2.)
          const1 = const1*Mstar
          const2 = 1.
          call prompt('enter adiabat (entropy normalisation)',const2,min=0.)
       endif
    endif

 end select

 if (prompt_for_gamma) then
    call prompt('set adiabatic gamma manually? (no=read from dumps)',use_gamma_exact)
    if (use_gamma_exact) then
       call prompt('enter gamma',gamma_exact)
    endif
 endif

 return
end subroutine submenu_exact

!---------------------------------------------------
! return number of exact solutions being plotted
!---------------------------------------------------
integer function get_nexact(iexact) result(nexact)
 integer, intent(in) :: iexact

 select case(iexact)
 case(2)
    nexact = nfiles
 case(1)
    nexact = nfunc
 case default
    nexact = 1
 end select

end function get_nexact

!---------------------------------------------------
! sets options relating to exact solution plotting
!---------------------------------------------------
subroutine options_exact(iexact)
 use prompting, only:prompt
 use plotlib,   only:plotlib_maxlinestyle,plotlib_maxlinecolour
 integer, intent(in) :: iexact
 character(len=12) :: string
 integer :: nexact,i

 call prompt('enter number of exact solution points ',maxexactpts,10,1000000)
 nexact = get_nexact(iexact)
 string = ''
 do i=1,nexact
    if (nexact > 1) write(string,"(a,i2)") 'for line ',i
    call prompt('enter line colour '//trim(string),iExactLineColour(i),0,plotlib_maxlinecolour)
 enddo
 do i=1,nexact
    if (nexact > 1) write(string,"(a,i2)") 'for line ',i
    call prompt('enter line style '//trim(string),iExactLineStyle(i),1,plotlib_maxlinestyle)
 enddo
 call prompt('Plot exact solution under data (default is over)? ',iPlotExactUnder)
 call prompt('Enter opacity for exact solution lines',ExactAlpha,0.0,1.0)
 do i=1,nexact
    call prompt('enter text to display in legend (blank=do not show)',ExactLegendText(i))
 enddo

 call prompt('calculate error norms?',iCalculateExactErrors)
 if (iCalculateExactErrors) then
    print "(/,' 0 : not normalised  L1 = 1/N \sum |y - y_exact|',/,"// &
                "' 1 : normalised      L1 = 1/N 1/max(y_exact) \sum |y - y_exact|')"
    call prompt('enter choice of error norm',iNormaliseErrors)
 endif
 if (iCalculateExactErrors) then
    call prompt('plot residuals (as inset in main plot)?',iPlotResiduals)
    if (iPlotResiduals) then
       call prompt('enter fraction of plot to use for inset', &
                      fracinsetResiduals,0.1,0.9)
       call prompt('enter max residual (0 for adaptive)',residualmax,0.)
    endif
 endif
 print "(/,'  0 : plot exact solution (where available) on every panel ',/,"// &
           "' -1 : plot exact solution on first row only ',/,"// &
           "' -2 : plot exact solution on first column only ',/,"// &
           "'  n : plot exact solution on nth panel only ')"

 call prompt('Enter selection ',iPlotExactOnlyOnPanel,-2)

end subroutine options_exact

 !-----------------------------------------------------------------------
 ! read exact solution parameters from files
 ! (in ndspmhd these files are used in the input to the code)
 !
 ! called after main data read and if exact solution chosen from menu
 !-----------------------------------------------------------------------
subroutine read_exactparams(iexact,rootname,ierr)
 use settings_data,  only:ndim,iverbose
 use prompting,      only:prompt
 use exactfunction,  only:check_function
 use filenames,      only:fileprefix,ifileopen
 use asciiutils,     only:read_asciifile,get_line_containing,read_var_from_file
 integer,          intent(in)  :: iexact
 character(len=*), intent(in)  :: rootname
 integer,          intent(out) :: ierr

 integer :: idash,nf,i,j,idrag,idum,linenum,k,ieq,ierrs(6),narmsread,nfiles_got
 character(len=len_trim(rootname)+8) :: filename
 character(len=120) :: line
 character(len=30)  :: var
 logical            :: iexist

 idash = index(rootname,'_',back=.true.)
 if (idash==0) idash = len_trim(rootname)+1

 select case(iexact)
 case(1)
    !
    !--read functions from file
    !
    filename=trim(rootname)//'.func'
    call read_asciifile(trim(filename),nf,funcstring,ierr)
    if (ierr==-1) then
       if (iverbose > 0) write(*,"(a)",advance='no') ' no file '//trim(filename)//'; '
       filename = trim(fileprefix)//'.func'
       call read_asciifile(trim(filename),nf,funcstring,ierr)
       if (ierr==-1) then
          if (iverbose > 0) print "(a)",' no file '//trim(filename)
          return
       endif
    endif

    if (nf > 0) then
       i = 0
       do while(i < nf)
          i = i + 1
          call check_function(funcstring(i),ierr,verbose=.false.)
          if (ierr /= 0) then
             print "(a)",' error parsing function '//trim(funcstring(i))//', skipping...'
             do j=i+1,nf
                funcstring(j-1) = funcstring(j)
             enddo
             funcstring(nf) = ' '
             nf = nf - 1
             i = i - 1
          endif
       enddo
       nfunc = nf
       if (iverbose > 0) print "(a,i2,a)",' read ',nfunc,' functions from '//trim(filename)
    else
       print "(a)",' *** NO FUNCTIONS READ: none will be plotted ***'
       ierr = 2
    endif
 case(2)
    filename = trim(fileprefix)//'.exactfiles'
    call read_asciifile(trim(filename),nfiles_got,filename_exact,ierr)
    if (ierr==-1 .and. nfiles <= 0) then ! look in .setup file only if no filenames set
       if (iverbose > 1) print "(a)",' no file '//trim(filename)
       !
       ! if no .exactfiles, see if the phantom .setup file exists
       ! and have a look for "outputfilename"
       !
       filename = trim(rootname(1:idash-1))//'.setup'
       call read_var_from_file('outputfilename',filename_exact(1),filename,ierr)
       if (ierr == 0) nfiles = 1
       return
    elseif (nfiles_got > 0) then
       nfiles = nfiles_got
    endif
 case(3,13)
    !
    !--shock tube parameters from .shk file
    !
    filename = trim(rootname(1:idash-1))//'.shk'
    inquire(file=filename,exist=iexist)
    if (iexist) then
       open(unit=19,file=filename,status='old',iostat=ierr)
       if (ierr==0) then
          read(19,*,iostat=ierrs(1)) rho_L, rho_R
          read(19,*,iostat=ierrs(2)) pr_L, pr_R
          read(19,*,iostat=ierrs(3)) v_L, v_R
          if (any(ierrs(1:3)/=0)) then
             print*,'error reading ',filename
             ierr = 1
          endif
       endif
       close(unit=19)
    else
       print*,'no file ',filename
       !
       ! look for .setup file for Phantom
       !
       filename= trim(rootname(1:idash-1))//'.setup'
       inquire(file=filename,exist=iexist)
       open(unit=19,file=filename,status='old',iostat=ierr)
       ierrs(:) = 0
       do while(ierr==0)
          read(19,"(a)",iostat=ierr) line
          ieq = index(line,'=')
          if (ierr==0 .and. ieq > 1) then
             var = trim(adjustl(line(1:ieq-1)))
             select case(trim(var))
             case('densleft')
                read(line(ieq+1:),*,iostat=ierrs(1)) rho_L
             case('densright')
                read(line(ieq+1:),*,iostat=ierrs(1)) rho_R
             case('prleft')
                read(line(ieq+1:),*,iostat=ierrs(1)) pr_L
             case('prright')
                read(line(ieq+1:),*,iostat=ierrs(1)) pr_R
             case('vxleft')
                read(line(ieq+1:),*,iostat=ierrs(1)) v_L
             case('vxright')
                read(line(ieq+1:),*,iostat=ierrs(1)) v_R
             case('dtg')
                read(line(ieq+1:),*,iostat=ierrs(1)) rdust_to_gas
             end select
          endif
       enddo
       ierr = 0
       if (.not.iexist .or. any(ierrs(1:6)/=0)) ierr = 1
    endif
    if (iexist) then
       print*,'>> read ',filename
       print*,' rhoL, rho_R = ',rho_L,rho_R
       print*,' pr_L, pr_R  = ',pr_L, pr_R
       print*,' v_L,  v_R   = ',v_L, v_R
    endif


 case(6)
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
8801   continue
       print*,'no file ',filename
       ierr = 1
       return
8888   print*,'error reading ',filename
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
       if (norder >= 0 .and. morder >= 0) then
          print*,' v = ',Atstar,' r'
       else
          print*,' vx = ',alphatstar,'x +',ctstar1,'y'
          print*,' vy = ',ctstar2,'x +',betatstar,'y'
       endif
       return
9901   continue
       print*,'no file ',filename
       ierr = 1
       return
9902   print*,'error reading ',filename
       close(UNIT=20)
       ierr = 2
       return
    end select
 case(8)
    !
    !--attempt to guess which MHD shock tube has been done from filename
    !
    !read(rootname(5:5),*,iostat=ios) ishk
    !if (ios /= 0) ishk = 1
    !
    !--prompt for shock type if not set
    !
    if (ishk <= 0) then ! prompt
       ishk = 1
       call prompt('enter shock solution to plot',ishk,1,7)
    endif
    return
 case(14)
    !
    !--dustywave parameters from ndspmhd or Phantom input file
    !
    filename = trim(rootname(1:idash-1))//'.in'
    var = 'Kdrag' !ndspmhd
    linenum = get_line_containing(filename,trim(var))
    if (linenum==0) then
       var = 'K_code' !phantom
       linenum = get_line_containing(filename,trim(var))
    endif
    open(unit=19,file=filename,status='old',iostat=ierr)
    if (ierr==0) then
       do i=1,linenum-1
          read(19,*,iostat=ierr)
       enddo
       read(19,"(a)",iostat=ierr) line
       if (ierr==0) then
          k = index(line,trim(var))
          if (k > 0) read(line(k+len_trim(var)+2:),*,iostat=ierr) Kdrag
          if (ierr==0) then
             print*,'>> read '//trim(var)//' = ',Kdrag,' from '//trim(filename)
          else
             read(line,*,iostat=ierr) idrag, idum, idum, Kdrag
             if (ierr==0) then
                print*,'>> read '//trim(var)//' = ',Kdrag,' from old-style '//trim(filename)
             else
                print*,'>> error reading Kdrag from '//trim(filename)
             endif
          endif
       else
          print*,'>> error reading Kdrag from '//trim(filename)
       endif
    endif
    close(unit=19)
    return
 case(17)
    !
    !--spiral arm parameters from .spirals file
    !
    select case(ispiral)
    case(2)
       filename=trim(rootname)//'.spirals'
       narmsread = 0
       call read_asciifile(trim(filename),narmsread,spiral_params,ierr)
       if (ierr==-1) then
          if (iverbose > 0) write(*,"(a)",advance='no') ' no file '//trim(filename)//'; '
          filename = trim(fileprefix)//'.spirals'
          call read_asciifile(trim(filename),narmsread,spiral_params,ierr)
          if (ierr==-1) then
             if (iverbose > 0) print "(a)",' no file '//trim(filename)
          else
             if (iverbose > 0) print*,trim(filename)//' read ',narmsread,' arms, err = ',ierr
             if (narmsread >= 0) narms = narmsread
          endif
       else
          if (iverbose > 0) print*,trim(filename)//' read ',narmsread,' arms, err = ',ierr
          narms = max(narmsread,0)
       endif
    case default
       if (ifileopen > 0 .and. iread_wakeparams) call get_planetdisc_parameters_from_data(rplanet,HonR,q_index,phase,ierr)
    end select
 end select

end subroutine read_exactparams

!-----------------------------------------------------------------------
! subroutine to automatically extract parameters needed for
! planet wake exact solution from the data
!-----------------------------------------------------------------------
subroutine get_planetdisc_parameters_from_data(rp,HR,q,angle,ierr)
 use physcon,       only:pi
 use filenames,     only:ifileopen
 use settings_data, only:ntypes,ndim,ndimV,ncolumns
 use particle_data, only:headervals,dat,iamtype,npartoftype
 use part_utils,    only:locate_nth_particle_of_type,get_binary
 use labels,        only:get_sink_type,ix,ivx,ipmass,headertags
 use asciiutils,    only:get_value
 use system_utils,  only:ienvlist
 real, intent(inout) :: rp,HR,q,angle
 integer, intent(out) :: ierr
 integer :: isinklist(2),isink1,isink2,itype,j,ntot
 logical :: got_sinks
 real :: x0(3),v0(3),domega,omega,polyk,cs

 ! do nothing if there is no data
 j = 1  ! use first data in memory
 if (ifileopen <= 0) return
 if (size(npartoftype(1,:)) < j) return

 ! either use --wake to get which planet, or --corotate
 isinklist = ienvlist('SPLASH_WAKE',2)
 if (all(isinklist==0)) isinklist = ienvlist('SPLASH_COROTATE',2)
 if (all(isinklist==0)) isinklist = (/1,2/)

 itype = get_sink_type(ntypes)
 got_sinks = .false.
 if (itype > 0) got_sinks = all(isinklist > 0) .and. (npartoftype(itype,j) > 0)
 if (got_sinks .and. all(ix > 0)) then
    call locate_nth_particle_of_type(isinklist(1),isink1,itype,iamtype(:,j),npartoftype(:,j),ntot)
    call locate_nth_particle_of_type(isinklist(2),isink2,itype,iamtype(:,j),npartoftype(:,j),ntot)
    call get_binary(isink1,isink2,dat(:,:,j),x0,v0,angle,domega,ndim,ndimV,ncolumns,ix,ivx,ipmass,0,ierr)
    if (ierr /= 0) return
    ! phase angle for planet
    angle = -angle*180./pi ! convert to degrees
    ! radial location of planet
    rp = norm2(dat(isink2,ix(:),j) - x0)
    ! power law index of sound speed profile, if found in file header
    q = get_value('qfacdisc',headertags,headervals(:,j))
    if (q <= 0.) q = q_index  ! revert to previous value if not read
    ! disc scale height, convert to aspect ratio
    polyk = 2./3.*get_value('RK2',headertags,headervals(:,j))  ! polytropic K
    cs = sqrt(polyk*(rp**2)**(-q))  ! sound speed in code units
    omega = sqrt(dat(isink1,ipmass,j)/rp**3) ! Keplerian speed
    HR = (cs/omega)/rp                       ! pressure scale height
    if (HR < 0.001 .or. HR > 0.5) HR = HonR  ! do not grab stupid numbers
    ! report results
    print "(4(a,g10.3))",' planet wake parameters: rp is ',rp,' phase is ',angle,' q is ',q,' HonR is ',HR
 endif

end subroutine get_planetdisc_parameters_from_data

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
                            ndim,ndimV,time,xmin,xmax,gamma,xplot,yplot,&
                            itag,iamtype,noftype,iplot_type, &
                            pmassmin,pmassmax,npart,imarker,unitsx,unitsy,irescale,iaxisy)
 use params,          only:int1,maxparttypes
 use settings_data,   only:ncolumns,ncalc
 use labels,          only:ix,irad,iBfirst,ivx,irho,ike,iutherm,ih,ipr,iJfirst,&
                              irhorestframe,is_coord,ideltav,idustfrac,label
 use filenames,       only:ifileopen,rootname
 use asciiutils,      only:string_replace,match_integer
 use prompting,       only:prompt
 use exactfromfile,   only:exact_fromfile
 use mhdshock,        only:exact_mhdshock
 use polytrope,       only:exact_polytrope
 use rhoh,            only:exact_rhoh
 use sedov,           only:exact_sedov
 use shock,           only:exact_shock
 use shock_sr,        only:exact_shock_sr
 use torus,           only:exact_torus
 use toystar1D,       only:exact_toystar1D  !, exact_toystar_ACplane
 use toystar2D,       only:exact_toystar2D
 use wave,            only:exact_wave
 use densityprofiles, only:exact_densityprofiles
 use exactfunction,   only:exact_function
 use ringspread,      only:exact_ringspread
 use dustywaves,      only:exact_dustywave
 use rochelobe,       only:exact_rochelobe
 use gresho,          only:exact_gresho
 use Cshock,          only:exact_Cshock
 use planetdisc,      only:exact_planetdisc
 use bondi,           only:exact_bondi
 use transforms,      only:transform,transform_inverse
 use map_columns,     only:map_columns_in_file
 use plotlib,         only:plot_qci,plot_qls,plot_sci,plot_sls,plot_line,plotlib_maxlinestyle,plot_set_opacity
 integer, intent(in) :: iexact,iplotx,iploty,itransx,itransy,igeom
 integer, intent(in) :: ndim,ndimV,npart,imarker,iaxisy
 real,    intent(in) :: time,xmin,xmax,gamma,unitsx,unitsy
 real,    intent(in) :: pmassmin,pmassmax
 real,    intent(in) :: xplot(npart),yplot(npart)
 integer, intent(in) :: itag(npart)
 integer(int1), intent(in) :: iamtype(:)
 integer,       intent(in) :: noftype(maxparttypes)
 logical,       intent(in) :: iplot_type(maxparttypes)
 logical,       intent(in) :: irescale
 integer :: mytransx,mytransy,imapx,imapy

 real, parameter :: zero = 1.e-10
 integer :: i,ierr,iexactpts,iCurrentColour,iCurrentLineStyle,LineStyle
 integer :: ncols,nrows,nlab_exact
 real, allocatable :: xexact(:),yexact(:),xtemp(:)
 real :: dx,timei,gammai
 character(len=len(filename_exact)) :: filename_tmp

 !
 !--change line style and colour settings, but save old ones
 !
 call plot_qci(iCurrentColour)
 call plot_qls(iCurrentLineStyle)
 call plot_sci(iExactLineColour(1))
 call plot_sls(iExactLineStyle(1))
 call plot_set_opacity(ExactAlpha)
 !
 !--get exact solution parameters from code input files, if present
 !
 if (ifileopen > 0) call read_exactparams(iexact,rootname(ifileopen),ierr)
 !
 !--allocate memory
 !
 allocate(xexact(maxexactpts),yexact(maxexactpts),xtemp(maxexactpts),stat=ierr)
 if (ierr /= 0) then
    print "(a)",'*** ERROR allocating memory for exact solution plotting, skipping ***'
    if (allocated(xexact)) deallocate(xexact)
    if (allocated(yexact)) deallocate(yexact)
    if (allocated(xtemp)) deallocate(xtemp)
    return
 endif

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
 xtemp = xexact
 if (itransx > 0) call transform_inverse(xexact,itransx)

 iexactpts = maxexactpts
 !
 !--exact solution plots must return a zero or negative value of ierr to be plotted
 !  (-ve ierr indicates a partial solution)
 !
 ierr = 666

 !
 !--use time=0 if time has not been read from dump file (indicated by t < 0)
 !
 if (time > 0) then
    timei = time
 else
    timei = 0.
 endif
 if (use_gamma_exact .and. gamma_exact >= 1.) then
    gammai = gamma_exact
 elseif (gamma >= 1.) then
    gammai = gamma
 else
    gammai = 5./3.
 endif

 select case(iexact)
 case(1) ! arbitrary function parsing
    do i=1,nfunc
       if ((iplotx==iexactplotx(i) .and. iploty==iexactploty(i)) .or. iexactploty(i)==0) then
          !--allow user-specified colours and line styles
          LineStyle = mod(iExactLineStyle(i),plotlib_maxlinestyle)

          call exact_function(funcstring(i),xexact,yexact,timei,ierr)
          !--plot each solution separately and calculate errors
          call plot_exact_solution(itransx,itransy,iexactpts,npart,xexact,yexact,xplot,yplot, &
                                   itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,&
                                   ls=LineStyle,lc=iExactLineColour(i))
          ierr = 1 ! indicate that we have already plotted the solution
       endif
    enddo
 case(2) ! exact solution read from file
    do i=1,nfiles
       !--substitute %f for filename
       filename_tmp = filename_exact(i)
       call string_replace(filename_tmp,'%f',trim(rootname(ifileopen)))
       !--auto-magically map columns onto splash columns (if not mapped manually)
       if (iauto_map_columns) then
          call map_columns_in_file(filename_tmp,ncols,nrows,imapexact,&
               label(1:ncolumns+ncalc),exact_labels,nlab_exact,ierr)
       endif
       ! see if there are columns in the exact file that match those of the current plot
       imapx = match_integer(imapexact,iplotx)
       imapy = match_integer(imapexact,iploty)
       ! if so, go ahead and read the data and plot the exact solution
       if (imapx > 0 .and. imapy > 0) then
          !--read exact solution from file
          print "(a)",'> reading '//trim(exact_labels(imapx))//' and '//trim(exact_labels(imapy))//' from '//trim(filename_tmp)
          call exact_fromfile(filename_tmp,xexact,yexact,imapx,imapy,iexactpts,ierr)
          !--plot this untransformed (as may already be in log space)
          if (ierr <= 0) then
             if (iApplyTransExactFile) then
                !--change into physical units if appropriate
                if (iRescale .and. iApplyUnitsExactFile) then
                   xexact(1:iexactpts) = xexact(1:iexactpts)*unitsx
                   yexact(1:iexactpts) = yexact(1:iexactpts)*unitsy
                endif
             endif
             !--change line style between files
             LineStyle = mod(iExactLineStyle(i),plotlib_maxlinestyle)
             !--do not apply log or other transformations if option set for this
             !  or if the label starts with 'log'...
             mytransy = get_transform(itransy,iApplyTransExactFile,exact_labels(imapy)(1:3))
             mytransx = get_transform(itransx,iApplyTransExactFile,exact_labels(imapx)(1:3))

             !--plot each solution separately and calculate errors
             call plot_exact_solution(mytransx,mytransy,iexactpts,npart,xexact,yexact,xplot,yplot,&
                                      itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,&
                                      ls=LineStyle,lc=iExactLineColour(i))
             ierr = 1 ! indicate that we have already plotted the solution
          endif
       else
          ! if nothing matched, make sure we do not plot anything
          ierr = 2
       endif
    enddo
 case(3)! shock tube
    if (iplotx==ix(1) .and. igeom <= 1) then
       do i=1,1  ! make 2 to plot both adiabatic and isothermal solutions
          if (i==2) then
             call plot_sls(4)
             gammai = 1.0
          endif
          if (iploty==irho) then
             call exact_shock(1,timei,gammai,xshock,rho_L,rho_R,pr_L,pr_R,v_L,v_R, &
                                 rdust_to_gas,xexact,yexact,ierr)
          elseif (iploty==ipr) then
             call exact_shock(2,timei,gammai,xshock,rho_L,rho_R,pr_L,pr_R,v_L,v_R, &
                                 rdust_to_gas,xexact,yexact,ierr)
          elseif (iploty==ivx) then
             call exact_shock(3,timei,gammai,xshock,rho_L,rho_R,pr_L,pr_R,v_L,v_R, &
                                 rdust_to_gas,xexact,yexact,ierr)
          elseif (iploty==iutherm) then
             call exact_shock(4,timei,gammai,xshock,rho_L,rho_R,pr_L,pr_R,v_L,v_R, &
                                 rdust_to_gas,xexact,yexact,ierr)
          elseif (iploty==ideltav) then
             call exact_shock(5,timei,gammai,xshock,rho_L,rho_R,pr_L,pr_R,v_L,v_R, &
                                 rdust_to_gas,xexact,yexact,ierr)
          elseif (iploty==idustfrac) then
             call exact_shock(6,timei,gammai,xshock,rho_L,rho_R,pr_L,pr_R,v_L,v_R, &
                                 rdust_to_gas,xexact,yexact,ierr)
          endif
          if (i==1 .and. ierr==0) call plot_line(iexactpts,xexact,yexact)
       enddo
    endif

 case(4)! sedov blast wave
    ! this subroutine does change xexact
    if (iplotx==irad .or. (igeom==3 .and. iplotx==ix(1))) then
       if (iploty==irho) then
          call exact_sedov(1,timei,gammai,rhosedov,esedov,xmax,xexact,yexact,ierr)
       elseif (iploty==ipr) then
          call exact_sedov(2,timei,gammai,rhosedov,esedov,xmax,xexact,yexact,ierr)
       elseif (iploty==iutherm) then
          call exact_sedov(3,timei,gammai,rhosedov,esedov,xmax,xexact,yexact,ierr)
       elseif (iploty==ike) then
          call exact_sedov(4,timei,gammai,rhosedov,esedov,xmax,xexact,yexact,ierr)
       elseif (iploty==ivx .and. igeom==3) then
          call exact_sedov(5,timei,gammai,rhosedov,esedov,xmax,xexact,yexact,ierr)
       endif
    elseif (igeom <= 1 .and. is_coord(iplotx,ndim) .and. is_coord(iploty,ndim)) then
       call exact_sedov(0,timei,gammai,rhosedov,esedov,xmax,xexact,yexact,ierr)
    endif

 case(5)! polytrope
    if (iploty==irho .and. (iplotx==irad .or.(igeom==3 .and. iplotx==ix(1)))) then
       call exact_polytrope(gammai,polyk,totmass,xexact,yexact,iexactpts,ierr)
    endif

 case(6)! toy star
    if (iBfirst /= 0) then
       sigma = sigma0
    else
       sigma = 0.
    endif
    if (ndim==1) then
       !
       !--1D toy star solutions
       !
       if (iplotx==ix(1) .or. iplotx==irad) then! if x axis is x or r
          if (iploty==irho) then
             call exact_toystar1D(1,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
          elseif (iploty==ipr) then
             call exact_toystar1D(2,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
          elseif (iploty==iutherm) then
             call exact_toystar1D(3,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
          elseif (iploty==ivx) then
             call exact_toystar1D(4,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
          elseif (iploty==iBfirst+1) then
             call exact_toystar1D(5,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
          endif
       elseif (iplotx==irho) then
          if (iploty==iBfirst+1) then
             call exact_toystar1D(6,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                   xexact,yexact,iexactpts,ierr)
          endif
       endif

       if (iploty==iacplane) then! plot point on a-c plane
          call exact_toystar1D(7,timei,gammai,htstar,atstar,ctstar,sigma,norder, &
                                xexact,yexact,iexactpts,ierr)
       endif
    else
       !
       !--2D toy star solutions
       !  these routines change xexact
       !
       if (igeom==1 .and.((iplotx==ix(1) .and. iploty==ivx) &
               .or. (iplotx==ix(2) .and. iploty==ivx+1))) then
          call exact_toystar2D(4,timei,gammai,polyk,totmass, &
                  atstar,htstar,ctstar,norder,morder, &
                  alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
       endif
       if (iplotx==irad .or. (igeom==2 .and. iplotx==ix(1))) then
          if (iploty==irho) then
             call exact_toystar2D(1,timei,gammai,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          elseif (iploty==ipr) then
             call exact_toystar2D(2,timei,gammai,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          elseif (iploty==iutherm) then
             call exact_toystar2D(3,timei,gammai,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          elseif (igeom==2 .and. iploty==ivx) then
             call exact_toystar2D(4,timei,gammai,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          elseif (iploty==ike) then
             call exact_toystar2D(5,timei,gammai,polyk,totmass, &
                     atstar,htstar,ctstar,norder,morder, &
                     alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
          endif
       elseif (is_coord(iplotx,ndim) .and. is_coord(iploty,ndim) .and. igeom==1) then
          call exact_toystar2D(0,timei,gammai,polyk,totmass, &
                  atstar,htstar,ctstar,norder,morder, &
                  alphatstar,betatstar,ctstar1,ctstar2,xexact,yexact,ierr)
       endif
    endif

 case(7)! Gresho vortex
    !if ((iploty==iwaveploty).and.(iplotx==iwaveplotx)) then
    !   ymean = SUM(yplot(1:npart))/REAL(npart)
    !   call exact_wave(timei,ampl,period,lambda,xzero,ymean,xexact,yexact,ierr)
    !endif
    if (igeom==2 .and. ndim >= 2) then
       if (iploty==ivx+1) then
          call exact_gresho(1,xexact,yexact,ierr)
       elseif (iploty==ipr) then
          call exact_gresho(2,xexact,yexact,ierr)
       endif
    endif
 case(8) ! mhd shock tubes
    ! this subroutine modifies xexact
    if (iplotx==ix(1) .and. igeom <= 1) then
       if (iploty==irho) then
          call exact_mhdshock(1,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==ipr) then
          call exact_mhdshock(2,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==ivx) then
          call exact_mhdshock(3,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==ivx+1 .and. ndimV > 1) then
          call exact_mhdshock(4,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==ivx+ndimV-1 .and. ndimV > 2) then
          call exact_mhdshock(5,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==iBfirst+1 .and. ndimV > 1) then
          call exact_mhdshock(6,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==iBfirst+ndimV-1 .and. ndimV > 2) then
          call exact_mhdshock(7,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==iutherm) then
          call exact_mhdshock(8,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       elseif (iploty==iBfirst) then
          call exact_mhdshock(9,ishk,timei,gammai,xmin,xmax,xshock, &
                                 xexact,yexact,iexactpts,ierr)
       endif
    endif
 case(9)
    !--h = (1/rho)^(1/ndim)
    if (((iploty==ih).and.(iplotx==irho)) .or. &
           ((iplotx==ih).and.(iploty==irho))) then
       if (iplotx==ih) then
          call exact_rhoh(2,ndim,hfact,pmassmin,xexact,yexact,ierr)
       else
          call exact_rhoh(1,ndim,hfact,pmassmin,xexact,yexact,ierr)
       endif

       !--if variable particle masses, plot one for each pmass value
       if (abs(pmassmin-pmassmax) > zero .and. pmassmin > zero) then
          !--plot first line
          if (ierr <= 0) then
             xtemp = xexact ! must not transform xexact as this is done again below
             if (itransx > 0) call transform(xtemp,itransx)
             if (itransy > 0) call transform(yexact,itransy)
             call plot_line(iexactpts,xtemp(1:iexactpts),yexact(1:iexactpts))
          endif
          !--leave this one to be plotted below
          if (iplotx==ih) then
             call exact_rhoh(2,ndim,hfact,pmassmax,xexact,yexact,ierr)
          else
             call exact_rhoh(1,ndim,hfact,pmassmax,xexact,yexact,ierr)
          endif
       endif
    endif
 case(10) ! density profiles
    if (iplotx==irad .or.(igeom==3 .and. iplotx==ix(1))) then
       if (iploty==irho) then
          call exact_densityprofiles(1,iprofile,Msphere,rsoft,xexact,yexact,ierr)
       elseif (iploty==icolpoten) then
          call exact_densityprofiles(2,iprofile,Msphere,rsoft,xexact,yexact,ierr)
       elseif (iploty==icolfgrav) then
          call exact_densityprofiles(3,iprofile,Msphere,rsoft,xexact,yexact,ierr)
       endif
    endif
 case(11) ! torus
    if (iplotx==irad .or.(igeom==3 .and. iplotx==ix(1))) then
       if (iploty==irho) then
          call exact_torus(1,1,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       elseif (iploty==ipr) then
          call exact_torus(2,1,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       elseif (iploty==iutherm) then
          call exact_torus(3,1,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       endif
       !--pr vs z at r=Rtorus
    elseif (igeom==2 .and. iplotx==ix(3) .and.iploty==ipr) then
       call exact_torus(4,1,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
    endif
    !--solutions for tokamak torus
    if (igeom==4 .and. iplotx==ix(1)) then
       if (iploty==irho) then
          call exact_torus(1,2,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       elseif (iploty==ipr) then
          call exact_torus(2,2,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       elseif (iploty==iutherm) then
          call exact_torus(3,2,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       elseif (iploty==iBfirst+1 .and. iBfirst > 0) then
          call exact_torus(4,2,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       elseif (iploty==iJfirst+2 .and. iJfirst > 0) then
          call exact_torus(5,2,Mstar,Rtorus,polyk,distortion,gammai,xexact,yexact,ierr)
       endif
    endif
 case(12)
    if (iplotx==irad .or.((igeom==3 .or. igeom==2) .and. iplotx==ix(1))) then
       if (iploty==irho) then
          call exact_ringspread(1,timei,Mring,Rring,viscnu,xexact,yexact,ierr)
       endif
    endif
 case(13) ! special relativistic shock tube
    if (iplotx==ix(1) .and. igeom <= 1) then
       if (iploty==irhorestframe) then
          call exact_shock_sr(1,timei,gammai,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
       elseif (iploty==ipr) then
          call exact_shock_sr(2,timei,gammai,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
       elseif (iploty==ivx) then
          call exact_shock_sr(3,timei,gammai,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
       elseif (iploty==iutherm) then
          call exact_shock_sr(4,timei,gammai,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
       elseif (iploty==irho) then
          call exact_shock_sr(5,timei,gammai,rho_L,rho_R,pr_L,pr_R,v_L,v_R,xexact,yexact,ierr)
       endif
    endif
 case(14) ! dusty wave exact solution
    if (iplotx==ix(1) .and. igeom <= 1) then
       if (iploty==ivx) then
          !--plot gas solution and calculate errors
          call exact_dustywave(1,timei,ampl,cs,Kdrag,lambda,xzero,rhozero,rhozero*rdust_to_gas,xexact,yexact,ierr)
          call plot_exact_solution(itransx,itransy,iexactpts,npart,xexact,yexact,xplot,yplot, &
                                      itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,ls=1,matchtype=1)
          !--plot dust solution
          if (Kdrag > 0.) then
             call exact_dustywave(2,timei,ampl,cs,Kdrag,lambda,xzero,rhozero,rhozero*rdust_to_gas,xexact,yexact,ierr)
             call plot_exact_solution(itransx,itransy,iexactpts,npart,xexact,yexact,xplot,yplot, &
                                         itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,ls=2,matchtype=2)
          endif
          ierr = 1
       elseif (iploty==irho) then
          !--plot gas solution and calculate errors
          call exact_dustywave(3,timei,ampl,cs,Kdrag,lambda,xzero,rhozero,rhozero*rdust_to_gas,xexact,yexact,ierr)
          call plot_exact_solution(itransx,itransy,iexactpts,npart,xexact,yexact,xplot,yplot,&
                                      itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,ls=1,matchtype=1)
          !--plot dust solution
          if (Kdrag > 0.) then
             call exact_dustywave(4,timei,ampl,cs,Kdrag,lambda,xzero,rhozero,rhozero*rdust_to_gas,xexact,yexact,ierr)
             call plot_exact_solution(itransx,itransy,iexactpts,npart,xexact,yexact,xplot,yplot,&
                                         itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,ls=2,matchtype=2)
          endif
          ierr = 1
       endif
    endif
 case(15) ! Roche potential
    if (igeom==1 .and. ndim >= 2 .and. iplotx==ix(1) .and. iploty==ix(2)) then
       call exact_rochelobe(xprim(1),xprim(2),xsec(1),xsec(2),mprim,msec,xexact,yexact,ierr)
    endif
 case(16) ! C-shock
    if (ndim >= 1 .and. iplotx==ix(1) .and. igeom <= 1) then
       if (iploty==irho) then
          call exact_Cshock(1,timei,gammai,machs,macha,xmin,xmax,xexact,yexact,ierr)
       elseif (iploty==iBfirst+1 .and. iBfirst > 0) then
          call exact_Cshock(2,timei,gammai,machs,macha,xmin,xmax,xexact,yexact,ierr)
       elseif (iploty==ivx) then
          call exact_Cshock(3,timei,gammai,machs,macha,xmin,xmax,xexact,yexact,ierr)
       elseif (iploty==ivx+1 .and. ndimV > 1) then
          call exact_Cshock(4,timei,gammai,machs,macha,xmin,xmax,xexact,yexact,ierr)
       elseif (iploty==iBfirst .and. iBfirst > 0) then
          call exact_Cshock(5,timei,gammai,machs,macha,xmin,xmax,xexact,yexact,ierr)
       endif
    endif
 case(17) ! planet-disc interaction
    if (ndim >= 2 .and. iplotx==ix(1) .and. iploty==ix(2)) then
       call exact_planetdisc(igeom,ispiral,iclockwise,timei,HonR,rplanet,q_index,phase,narms,&
                             spiral_params,xexact,yexact,ierr)
    endif
 case(18)
    if (iplotx==irad .or. (igeom==3 .and. iplotx==ix(1))) then
       if (iploty==ivx .and. igeom==3) then
          call exact_bondi(1,timei,gammai,const1,const2,Mstar,relativistic,geodesic_flow,is_wind,xexact,yexact,ierr)
       elseif (iploty==iutherm .and. igeom==3) then
          call exact_bondi(2,timei,gammai,const1,const2,Mstar,relativistic,geodesic_flow,is_wind,xexact,yexact,ierr)
       elseif (iploty==irho .and. igeom==3) then
          call exact_bondi(3,timei,gammai,const1,const2,Mstar,relativistic,geodesic_flow,is_wind,xexact,yexact,ierr)
       endif
    endif
 end select

 !----------------------------------------------------------
 !  plot this as a line on the current graph
 !----------------------------------------------------------
 if (ierr <= 0) call plot_exact_solution(itransx,itransy,iexactpts,npart,xexact,yexact,xplot,yplot, &
                                            itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy)
 !
 !--reset line and colour settings
 !
 call plot_sci(iCurrentColour)
 call plot_sls(iCurrentLineStyle)
 call plot_set_opacity(1.0)
 !
 !--deallocate memory
 !
 if (allocated(xexact)) deallocate(xexact)
 if (allocated(yexact)) deallocate(yexact)
 if (allocated(xtemp)) deallocate(xtemp)

end subroutine exact_solution

integer function get_transform(itrans,apply_trans,label_start) result(mytrans)
 use transforms, only:islogged
 use asciiutils, only:lcase
 integer, intent(in) :: itrans
 logical, intent(in) :: apply_trans
 character(len=3), intent(in) :: label_start

 ! by default, the data transformation is also applied to the exact solution
 mytrans = itrans
 ! disable this if the flag to switch it off is set
 if (.not.apply_trans) mytrans = 0

 if (lcase(label_start)=='log') then
    if (islogged(itrans)) then
       ! disable exact solution transform if solution is already logged
       mytrans = 0
    else
       ! if exact solution is logged but data isn't, apply inverse transform
       mytrans = -1
    endif
 endif

end function get_transform

!------------------------------------------------------------------
! Wrapper routine to plot the exact solution line on current graph
! and calculate errors with respect to the data
!------------------------------------------------------------------
subroutine plot_exact_solution(itransx,itransy,iexactpts,np,xexact,yexact,xplot,yplot,&
                               itag,iamtype,noftype,iplot_type,xmin,xmax,imarker,iaxisy,&
                               ls,lc,matchtype,err)
 use transforms, only:transform,transform_inverse
 use plotlib,    only:plot_line,plot_sls,plot_sci,plot_qci,plot_set_opacity
 use params,     only:int1,maxparttypes
 integer, intent(in)    :: itransx,itransy,iexactpts,np,imarker,iaxisy
 real,    intent(inout) :: xexact(:),yexact(:)
 real,    intent(in)    :: xplot(:),yplot(:),xmin,xmax
 integer, intent(in)    :: itag(:)
 integer(int1), intent(in) :: iamtype(:)
 integer,       intent(in) :: noftype(maxparttypes)
 logical,       intent(in) :: iplot_type(maxparttypes)
 integer,       intent(in), optional :: ls, matchtype, lc
 logical,       intent(in), optional :: err
 real :: residuals(np),ypart(np)
 real :: errL1,errL2,errLinf
 integer :: iused,ierr,iCurrentColour,imatchtype
 logical :: plot_err
 character(len=12) :: str1,str2

 call plot_qci(iCurrentColour)
 if (present(lc)) call plot_sci(lc)
 call plot_set_opacity(ExactAlpha)

 if (itransx /= 0) call transform(xexact(1:iexactpts),itransx)
 if (itransy /= 0) call transform(yexact(1:iexactpts),itransy)

 if (present(ls)) call plot_sls(ls)

 call plot_line(iexactpts,xexact(1:iexactpts),yexact(1:iexactpts))
 !
 !--calculate errors
 !
 if (present(err)) then
    plot_err = err
 else
    plot_err = .true.
 endif
 if (present(matchtype)) then
    imatchtype = matchtype
 else
    imatchtype = 0
 endif

 if (iCalculateExactErrors .and. plot_err) then
    !--untransform y axis again for error calculation
    if (itransy /= 0) call transform_inverse(yexact(1:iexactpts),itransy)
    !--untransform particle y axis also
    ypart(1:np) = yplot(1:np)
    if (itransy /= 0) call transform_inverse(ypart(1:np),itransy)
    !--calculate errors
    call calculate_errors(xexact(1:iexactpts),yexact(1:iexactpts),xplot(1:np),ypart,&
                            itag,iamtype,noftype,iplot_type,xmin,xmax,residuals, &
                            errL1,errL2,errLinf,iused,imatchtype)
    if (iused /= np) then
       write(str1,"(i12)",iostat=ierr) iused
       write(str2,"(i12)",iostat=ierr) np
       print "(3(a,es12.5,1x),'(used ',a,'/',a,' parts)')",' L1 err = ',errL1,'L2 err = ',errL2, &
                               'L(inf) err = ',errLinf,trim(adjustl(str1)),trim(adjustl(str2))
    else
       print "(3(a,es12.5,1x))",' L1 err = ',errL1,'L2 err = ',errL2, &
                               'L(inf) err = ',errLinf
    endif
    if (iPlotResiduals) then
       call plot_sci(1)
       call plot_residuals(xplot,residuals,imarker,iaxisy)
    endif
 endif
 call plot_sci(iCurrentColour)
 call plot_set_opacity(1.0)

end subroutine plot_exact_solution

 !--------------------------------
 ! Calculate various error norms
 !--------------------------------
subroutine calculate_errors(xexact,yexact,xpts,ypts,itag,iamtype,noftype,iplot_type,&
                              xmin,xmax,residual,errL1,errL2,errLinf,iused,imatchtype)
 use part_utils, only:igettype
 use params,     only:int1,maxparttypes
 real,    intent(in)  :: xexact(:),yexact(:),xpts(:),ypts(:),xmin,xmax
 integer, intent(in)  :: itag(:)
 integer(int1), intent(in) :: iamtype(:)
 integer,       intent(in) :: noftype(maxparttypes)
 logical,       intent(in) :: iplot_type(maxparttypes)
 real,    intent(out) :: residual(size(xpts))
 real,    intent(out) :: errL1,errL2,errLinf
 integer, intent(out) :: iused
 integer, intent(in)  :: imatchtype
 integer :: i,j,npart,nerr,nused,itype
 real    :: xi,dy,dx,yexacti,err1,ymax
 logical :: mixedtypes

 errL1 = 0.
 errL2 = 0.
 errLinf = 0.
 residual = 0.
 npart = size(xpts)
 iused = 0
 ymax = 0.
 nerr = 0
 nused = 0
 mixedtypes = size(iamtype) > 1

 do i=1,npart
    xi = xpts(i)
    yexacti = 0.
    if (xi >= xmin .and. xi <= xmax .and. itag(i) > 0) then
       if (mixedtypes) then
          itype = min(max(int(iamtype(i)),1),maxparttypes)
       else
          itype = igettype(i,noftype)
       endif
       if (iplot_type(itype) .and. ((imatchtype==0).or.(itype==imatchtype))) then
          !
          !--find nearest point in exact solution table
          !
          do j=1,size(xexact)-1
             if ((xexact(j) <= xi .and. xexact(j+1) > xi) .or. &
                 (xexact(j+1) <= xi .and. xexact(j) > xi)) then
                if (abs(residual(i)) > tiny(residual)) nerr = nerr + 1
                !--linear interpolation from tabulated exact solution
                dy = yexact(j+1) - yexact(j)
                dx = xexact(j+1) - xexact(j)
                if (dx > 0.) then
                   yexacti = yexact(j) + dy/dx*(xi - xexact(j))
                   residual(i) = ypts(i) - yexacti
                elseif (dy >= 0.) then
                   yexacti = yexact(j)
                   residual(i) = ypts(i) - yexacti
                else
                   nerr = nerr + 1
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
          if (yexacti > tiny(yexacti)) residual(i) = residual(i)/abs(yexacti)
          nused = nused + 1
       endif
    endif
 enddo
 !
 !--normalise errors (use maximum y value)
 !
 if (nused > 0) then
    if (iNormaliseErrors==1 .and. ymax > tiny(ymax)) then
       errL1 = errL1/(nused*ymax)
       errL2 = sqrt(errL2/(nused*ymax**2))
       errLinf = errLinf/ymax
    else
       errL1 = errL1/nused
       errL2 = sqrt(errL2/(nused))
    endif
 else
    !print "(a)",' error normalising errors'
    errL1 = 0.
    errL2 = 0.
    errLinf = 0.
 endif
 if (nerr > 0) print*,'WARNING: ',nerr,' errors in residual calculation'

end subroutine calculate_errors

 !------------------------------------
 ! Plot residual errors as inset plot
 !------------------------------------
subroutine plot_residuals(xpts,residuals,imarker,iaxisy)
 use plotlib, only:plot_qvp,plot_qwin,plot_svp,plot_qci,plot_qfs, &
                     plot_qcs,plot_sci,plot_sfs,plot_svp,plot_box, &
                     plot_pt,plot_swin,plot_rect
 real,    intent(in) :: xpts(:),residuals(:)
 integer, intent(in) :: imarker,iaxisy
 real :: vptxminold,vptxmaxold,vptyminold,vptymaxold
 real :: vptxmin,vptxmax,vptymin,vptymax
 real :: xminold,xmaxold,yminold,ymaxold,ymin,ymax
 real :: xch,ych
 integer :: ioldcolour,ioldfill

 !--query old viewport and window size
 call plot_qvp(0,vptxminold,vptxmaxold,vptyminold,vptymaxold)
 call plot_qwin(xminold,xmaxold,yminold,ymaxold)

 !--use specified bottom % of viewport
 vptxmin = vptxminold
 vptxmax = vptxmaxold
 vptymin = vptyminold
 vptymax = vptyminold + FracinsetResiduals*(vptymaxold - vptyminold)
 call plot_svp(vptxmin,vptxmax,vptymin,vptymax)

 !--set window
 if (residualmax < tiny(residualmax)) then
    ymax = maxval(abs(residuals))
    print*,'max residual = ',ymax
 else
    ymax = residualmax
 endif
 ymin = -ymax

 !--erase space for residual plot
 call plot_qci(ioldcolour)
 call plot_qfs(ioldfill)
 call plot_qcs(0,xch,ych)
 call plot_sci(0)
 call plot_sfs(1)
 if (iaxisy < 0) then
    call plot_svp(vptxmin,vptxmax,vptymin,vptymax)
 else
    call plot_svp(vptxmin - 3.*xch,vptxmax,vptymin,vptymax)
 endif
 call plot_swin(xminold,xmaxold,ymin,ymax)
 call plot_rect(xminold,xmaxold,ymin,ymax)
 !--restore fill style
 call plot_sfs(ioldfill)
 call plot_sci(1)

 !--set window and draw axes
 call plot_svp(vptxmin,vptxmax,vptymin,vptymax)
 call plot_swin(xminold,xmaxold,ymin,ymax)
 if (iaxisy < 0) then
    call plot_box('ABCST',0.0,0,'BCST',0.0,0)
 else
    call plot_box('ABCST',0.0,0,'BVNCST',0.0,0)
 endif

 !--plot residuals
 call plot_sci(ioldcolour)
 call plot_pt(size(xpts),xpts,residuals,imarker)

 !--restore old viewport, window and colour index
 call plot_svp(vptxminold,vptxmaxold,vptyminold,vptymaxold)
 call plot_swin(xminold,xmaxold,yminold,ymaxold)

end subroutine plot_residuals

end module exact
