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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------
!     module implementing the ability to use SPLASH to produce
!     evolution files from a sequence of SPH dump files
!     (ie. in order to produce plots of certain quantities vs time)
!     Command is "splash calc X" where X is analysis type.
!
!     (c) D. Price 06/11/08
!-------------------------------------------------------------------
module analysis
 public :: isanalysis,open_analysis,close_analysis,write_analysis
 private
 integer, private, parameter :: iunit = 89
 integer, private, parameter :: maxlevels = 20
!
! default settings for the density thresholds for massaboverho output
!
 integer, private :: nlevels,nfilesread
 real, dimension(maxlevels), private :: rholevels
 character(len=64), private :: fileout
 real, dimension(:,:), allocatable :: datmean,datvar

contains

!-----------------------------------------------------------------
! utility to check if the choice of analysis type is valid
! and if not, to print the available options
!-----------------------------------------------------------------
logical function isanalysis(string,noprint)
 implicit none
 character(len=*), intent(in) :: string
 logical, intent(in), optional :: noprint
 logical :: doprint

 isanalysis = .false.
 select case(trim(string))
 case('energies','energy')
     isanalysis = .true.
 case('massaboverho')
     isanalysis = .true.
 case('max','maxvals')
     isanalysis = .true.
 case('min','minvals')
     isanalysis = .true.
 case('mean','meanvals')
     isanalysis = .true.
 case('rms','rmsvals')
     isanalysis = .true.
 case('vrms','vrmsvals','vwrms','rmsvw')
     isanalysis = .true.
 case('rhovar','rhomach')
     isanalysis = .true.
 case('timeaverage','timeav')
     isanalysis = .true.
 end select
 
 if (present(noprint)) then
    doprint = .not.noprint
 else
    doprint = .true.
 endif
 
 if (.not.isanalysis .and. doprint) then
    print "(a)",' Analysis mode ("splash calc X dumpfiles") on a sequence of dump files: '
    print "(a)",'  splash calc energies     : calculate KE,PE,total energy vs time'
    print "(a)",'                             output to file called ''energy.out'''
    print "(a)",'         calc massaboverho : mass above a series of density thresholds vs time'
    print "(a)",'                             output to file called ''massaboverho.out'''
    print "(a)",'         calc rhomach      : density variance and RMS velocity dispersion vs. time'
    print "(a)",'                             output to file called ''rhomach.out'''
    print "(a)",'         calc max          : maximum of each column vs. time'
    print "(a)",'                             output to file called ''maxvals.out'''
    print "(a)",'         calc min          : minimum of each column vs. time'
    print "(a)",'                             output to file called ''minvals.out'''
    print "(a)",'         calc mean         : mean of each column vs. time'
    print "(a)",'                             output to file called ''meanvals.out'''
    print "(a)",'         calc rms          : (mass weighted) root mean square of each column vs. time'
    print "(a)",'                             output to file called ''rmsvals.out'''
    print "(a)",'         calc vrms         : volume weighted root mean square of each column vs. time'
    print "(a)",'                             output to file called ''rmsvals-vw.out'''
    print "(/,a)",'  the above options all produce a small ascii file with one row per input file.'
    print "(a)",'  the following option produces a file equivalent in size to one input file (in ascii format):'
    print "(/,a)",'         calc timeaverage  : time average of *all* entries for every particle'
    print "(a)",'                             output to file called ''time_average.out'''
 endif
 
 return
end function isanalysis

!----------------------------------------------------------------
!  open output file/ initialise quantities needed for analysis 
!  over all dump files
!----------------------------------------------------------------
subroutine open_analysis(dumpfile,analysistype,required,ncolumns,ndimV)
 use labels,     only:ivx,iBfirst,iutherm,irho,ipmass,label
 use asciiutils, only:read_asciifile
 use filenames,  only:rootname,nfiles,tagline
 implicit none
 integer, intent(in) :: ncolumns,ndimV
 character(len=*), intent(in) :: dumpfile,analysistype
 logical, dimension(0:ncolumns), intent(out) :: required
 character(len=1170) :: headerline   ! len=64 x 18 characters
 character(len=64) :: levelsfile
 character(len=40) :: fmtstring
 logical :: iexist,standardheader
 integer :: ierr,i
!
!--the 'required' array is used by the data reads (where implemented)
!  to determine whether or not we actually need to read a given column
!  from the file -- if not it can be skipped, leading to a faster
!  data read. Here we want to specify which columns are required
!  for the analysis in question.
!
 print "(/,5('-'),a,/)",'> CALCULATING '//trim(analysistype)//' vs time for all dump files'
 required(:)=.false.
 headerline = ' '
 standardheader = .false.

 select case(trim(analysistype))
 case('energy','energies')
    !
    !--for energies need to read particle mass, velocity, utherm and if present,
    !  magnetic field and density. The obvious limitation here is that we
    !  cannot calculate the potential energy unless it is dumped and labelled
    !  (which is not currently implemented).
    !
    required(ivx:ivx+ndimV-1) = .true.
    required(iBfirst:iBfirst+ndimV-1) = .true.
    required(iutherm) = .true. 
    required(ipmass) = .true.
    if (iBfirst.gt.0) required(irho) = .true.
    !
    !--set filename and header line
    !
    fileout = 'energy.out'
    write(headerline,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',2,'ekin',3,'etherm',4,'emag',5,'epot',6,'etot',7,'totmom'

 case('massaboverho')
    !
    !--only need to read mass and density from dump files
    !
    required(ipmass) = .true.
    required(irho) = .true.
    !
    !--need a user-configurable way of setting the density thresholds:
    !  implemented by setting them in a file which is read here
    !  (a new one is created if it doesn't exist)
    !
    levelsfile = 'massaboverho.levels'
    inquire(file=levelsfile,exist=iexist)
    if (iexist) then
       call read_asciifile(trim(levelsfile),nlevels,rholevels)
       print "(a)",' read '//trim(levelsfile)//':'
       do i=1,nlevels
          print "(a,i2,a,es8.2)",' level ',i,': rho > ',rholevels(i)
       enddo
       print*
    else
       print "(a)",' SPLASH ANALYSIS: levels file '//trim(levelsfile)//' not found'
       print "(a)",' creating one with default levels for mass > rho'
       print "(a)",' edit this file to set the density levels'
       open(unit=iunit+1,file=levelsfile,status='new',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          stop 'ERROR creating levels file'
       else
          nlevels = 10
          rholevels(1:nlevels) = (/1.e-20,1.e-19,1.e-18,1.e-17,1.e-16, &
                                   1.e-15,1.e-14,1.e-13,1.e-12,1.e-11/)
          write(iunit+1,*) rholevels(1:nlevels)
          close(iunit+1)
       endif
       stop
    endif
    !
    !--set filename and header line
    !
    fileout = 'massaboverho.out'
    write(headerline,"('#',1x,'[',i2.2,1x,a12,']',1x,20('[',i2.2,1x,a4,es8.1,a1,']',1x))") &
          1,'time',(i+1,'M(r>',rholevels(i),')',i=1,nlevels)

 case('max','maxvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'maxvals.out'
    standardheader = .true.

 case('min','minvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'minvals.out'
    standardheader = .true.

 case('mean','meanvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'meanvals.out'
    standardheader = .true.

 case('rms','rmsvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'rmsvals.out'
    standardheader = .true.

 case('vrms','vrmsvals','vwrms','rmsvw')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'rmsvals-vw.out'
    standardheader = .true.

 case('rhovar','rhomach')
    !
    !--read density, velocity info
    !
    required(ipmass) = .true.
    required(irho) = .true.
    required(ivx:ivx+ndimV-1) = .true.
    !
    !--set filename and header line
    !
    fileout = 'rhomach.out'
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) 17
    write(headerline,fmtstring) 1,'time',2,'rhomean(vw)',3,'rhomean(mw)',4,'varrho(vw)',5,'varrho(mw)',&
          6,'stddevrho(vw)',7,'stddevrho(mw)',8,'rms v (vw)',9,'rms v (mw)',10,'b (vw)',11,'b (mw)',&
          12,'s mean(vw)',13,'s mean(mw)',14,'s var(vw)',15,'s var(mw)',16,'s stddev(vw)',17,'s stddev(mw)'

 case('timeaverage','timeav')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'time_average.out'
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']''))')",iostat=ierr) 2*ncolumns
    write(headerline,fmtstring) (i,label(i)(1:12),i=1,ncolumns),&
                                (ncolumns+i,'err'//label(i)(1:9),i=1,ncolumns)
 end select

 if (standardheader) then
!
!--standard header is time in column 1, with an entry for each column following
!  (this is to avoid repeated code above)
!
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) ncolumns+1
    write(headerline,fmtstring) 1,'time',(i+1,label(i)(1:12),i=1,ncolumns)
 endif
 
!
!--do not replace the file if it already exists
!
 inquire(file=trim(fileout),exist=iexist)
 if (iexist) then
    print "(2(a,/))",' ERROR: analysis file '//trim(fileout)//' already exists', &
                    '        delete, move or rename this file and try again'
    stop
 endif
!
!--open the file for output
!
 open(unit=iunit,file=trim(fileout),status='new',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening file '//trim(fileout)//' for output'
    stop
 endif
 print "(a)",' WRITING '//trim(analysistype)//' vs time '//' TO FILE '//trim(fileout)
!
!--write header if the headerline is set
!  (no header is written if headerline is blank)
!
 if (len_trim(headerline).gt.0) then
    write(iunit,"(a)") '# '//trim(tagline)
    write(iunit,"(a)") '# '//trim(fileout)//' produced using "splash calc '//trim(analysistype)// &
                       '" on dump files '//trim(rootname(1))//'->'//trim(rootname(nfiles))
    write(iunit,"(a)") '# use asplash -e '//trim(fileout)//' to plot the contents of this file '
    write(iunit,"(a)") '#'
    write(iunit,"(a)") trim(headerline)
 endif

 nfilesread = 0

 return
end subroutine open_analysis

!----------------------------------------------------------------
!  this is the routine which actually calculates the quantities
!  required from each dump file and spits out a line to the
!  analysis file. Called once for each dump file.
!----------------------------------------------------------------
subroutine write_analysis(time,dat,ntot,ntypes,npartoftype,massoftype,iamtype,ncolumns,ndimV,analysistype)
 use labels,       only:ivx,iBfirst,iutherm,irho,ipmass,labeltype,label
 use params,       only:int1,doub_prec,maxplot
 use asciiutils,   only:ucase
 use system_utils, only:renvironment
 implicit none
 integer, intent(in)               :: ntot,ntypes,ncolumns,ndimV
 integer, intent(in), dimension(:) :: npartoftype
 real, intent(in), dimension(:)    :: massoftype
 integer(kind=int1), intent(in), dimension(:) :: iamtype
 real, intent(in)                  :: time
 real, intent(in), dimension(:,:)  :: dat
 character(len=*), intent(in)      :: analysistype
 real(kind=doub_prec), dimension(maxlevels) :: massaboverho
 integer              :: itype,i,j,ierr,ntot1,ncol1
 real(kind=doub_prec) :: ekin,emag,etherm,epot,etot,totmom,pmassi
 real(kind=doub_prec) :: rmsval,totvol,voli,rhoi,rmsvmw,v2i
 real(kind=doub_prec) :: rhomeanmw,rhomeanvw,rhovarmw,rhovarvw,bval,bvalmw
 real(kind=doub_prec) :: smeanmw,smeanvw,svarmw,svarvw,si
 real(kind=doub_prec), dimension(3) :: xmom
 real                 :: delta,dn
 character(len=20)    :: fmtstring
!
! array with one value for each column
!
 real, dimension(maxplot) :: coltemp

 nfilesread = nfilesread + 1
 if (time.ge.0.) then
    print "(/,5('-'),a,', TIME=',es8.2,' FILE #',i4,/)",&
          '> CALCULATING '//trim(ucase(analysistype)),time,nfilesread
 else
    print "(/,5('-'),a,', FILE #',i4,' (TIME NOT READ)'/)",&
          '> CALCULATING '//trim(ucase(analysistype)),nfilesread
 endif
 
 select case(trim(analysistype))
 case('energy','energies')
    ekin = 0.
    emag = 0.
    epot = 0.
    etherm = 0.
    etot = 0.
    do i=1,ntot
       itype = igettype(i)
       pmassi = particlemass(i,itype)
       
       !--kinetic energy
       if (ivx.gt.0 .and. ivx+ndimV-1.le.ncolumns) then
          ekin = ekin + pmassi*dot_product(dat(i,ivx:ivx+ndimV-1),dat(i,ivx:ivx+ndimV-1))
          xmom(1:ndimV) = xmom(1:ndimV) + dat(i,ivx:ivx+ndimV-1)
       endif
       
       !--thermal energy
       if (iutherm.gt.0 .and. iutherm.le.ncolumns) then
          etherm = etherm + pmassi*dat(i,iutherm)
       endif
       
       !--magnetic energy
       if (iBfirst.gt.0 .and. iBfirst+ndimV-1.le.ncolumns) then
          emag = emag + pmassi*dot_product(dat(i,iBfirst:iBfirst+ndimV-1),dat(i,iBfirst:iBfirst+ndimV-1))/dat(i,irho)
       endif
    enddo
    ekin = 0.5*ekin
    emag = 0.5*emag
    etot = ekin + etherm + epot + emag
    totmom = sqrt(dot_product(xmom(1:ndimV),xmom(1:ndimV)))

    print "(6(/,1x,a6,' = ',es8.2))",'etot',etot,'ekin',ekin,'etherm',etherm,'epot',epot,'emag',emag,'totmom',totmom
    !
    !--write line to output file
    !
    write(iunit,"(64(es18.10,1x))") time,ekin,etherm,emag,epot,etot,totmom

 case('massaboverho')
    massaboverho(:) = 0.
    if (irho.gt.0 .and. irho.le.ncolumns) then
       !
       !--warn if particle masses not found
       !
       if (ipmass.le.0 .or. ipmass.gt.ncolumns .and. all(massoftype.eq.0)) then
          print "(a)",' WARNING in massaboverho analysis!'// &
                      ' masses not read or are zero from dump file'
       endif
       !
       !--calculate mass above each density threshold
       !
       do i=1,ntot
          itype = igettype(i)
          pmassi = particlemass(i,itype)

          if (itype.eq.1) then
          !
          !--gas particles contribute if they are above rho
          !
             where(dat(i,irho).ge.rholevels(1:nlevels))
                massaboverho(1:nlevels) = massaboverho(1:nlevels) + pmassi
             end where
          elseif (labeltype(itype).eq.'sink') then
          !
          !--sink particles always contribute (ie. they are assumed to
          !  be above every density threshold)
          !
             massaboverho(1:nlevels) = massaboverho(1:nlevels) + pmassi
          endif
       enddo
       
       !
       !--write output to screen/terminal
       !
       do i=1,nlevels
          print "(1x,'M(rho > ',es8.2,') = ',es8.2)",rholevels(i),massaboverho(i)
       enddo
       !
       !--write line to output file
       !
       write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) nlevels+1
       write(iunit,fmtstring) time,massaboverho(1:nlevels)
       
    else
       print "(a)",' ERROR in massaboverho analysis!'// &
                   ' either mass or density not found in dump file'
       return
    endif
 case('max','maxvals')
    !
    !--calculate maximum for each column
    !
    do i=1,ncolumns
       coltemp(i) = maxval(dat(1:ntot,i))
    enddo
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'max = ',es9.2)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) time,coltemp(1:ncolumns)

 case('min','minvals')
    !
    !--calculate minimum for each column
    !
    do i=1,ncolumns
       coltemp(i) = minval(dat(1:ntot,i))
    enddo
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'min = ',es9.2)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) time,coltemp(1:ncolumns)

 case('mean','meanvals')
    !
    !--calculate mean for each column
    !
    do i=1,ncolumns
       coltemp(i) = sum(dat(1:ntot,i))/real(ntot)
    enddo
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'mean = ',es9.2)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) time,coltemp(1:ncolumns)

 case('rms','rmsvals')
    !
    !--calculate RMS for each column
    !
    do i=1,ncolumns
       coltemp(i) = sqrt(sum(dat(1:ntot,i)**2)/real(ntot))
    enddo
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'rms (mass weighted) = ',es9.2)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) time,coltemp(1:ncolumns)

 case('vrms','vrmsvals','vwrms','rmsvw')
    if (irho.le.0 .or. irho.gt.ncolumns) then
       print "(a)",' ERROR in volume weighted rms calculation!'// &
                   ' density not present / not labelled in dump file, skipping...'
       return
    endif
    !
    !--warn if particle masses not found
    !
    if (ipmass.le.0 .or. ipmass.gt.ncolumns .and. all(massoftype.eq.0)) then
       print "(a)",' WARNING in volume weighted rms calculation!'// &
                   ' masses not read or are zero from dump file'
    endif
    !
    !--calculate volume-weighted RMS for each column
    !
    do j=1,ncolumns
       rmsval = 0.
       totvol = 0.
       do i=1,ntot
          itype  = igettype(i)
          pmassi = particlemass(i,itype)
          rhoi   = dat(i,irho)
          if (rhoi.gt.0.) then
             voli = pmassi/rhoi
          else
             voli = 0.
          endif
          
          rmsval = rmsval + voli*dat(i,j)**2
          totvol = totvol + voli
       enddo
       coltemp(j) = real(sqrt(rmsval/totvol))
    enddo
    print "(1x,a,es9.2)",'volume = ',totvol
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'rms (volume weighted) = ',es9.2)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) time,coltemp(1:ncolumns)

 case('rhovar','rhomach')
    if (irho.le.0 .or. irho.gt.ncolumns) then
       print "(a)",' ERROR in density variance--rms velocity field calculation!'// &
                   ' density not present / not labelled in dump file, skipping...'
       return
    endif
    !
    !--warn if particle masses not found
    !
    if (ipmass.le.0 .or. ipmass.gt.ncolumns .and. all(massoftype.eq.0)) then
       print "(a)",' WARNING in volume weighted rms calculation!'// &
                   ' masses not read or are zero from dump file'
    endif

    if (ivx.le.0 .or. ivx.gt.ncolumns) then
       print "(a)",' WARNING in volume weighted rms calculation!'// &
                   ' velocities not present / not labelled in dump file'
    endif
    !
    !--calculate mean density and rms velocity values on first pass
    !
    rmsval = 0.
    rmsvmw = 0.
    rhomeanvw = 0.
    rhomeanmw = 0.
    totvol = 0.
    smeanvw = 0.
    smeanmw = 0.
    do i=1,ntot
       itype  = igettype(i)
       pmassi = particlemass(i,itype)
       rhoi   = dat(i,irho)
       if (rhoi.gt.0.) then
          voli = pmassi/rhoi
       else
          voli = 0.
       endif
       rhomeanmw = rhomeanmw + rhoi
       rhomeanvw = rhomeanvw + pmassi
       si        = log(rhoi)
       smeanmw   = smeanmw + si
       smeanvw   = smeanvw + voli*si
       totvol = totvol + voli

       !
       !--mean squared velocity
       !
       if (ivx.gt.0 .and. ivx.le.ncolumns) then
          v2i    = dot_product(dat(i,ivx:ivx+ndimV-1),dat(i,ivx:ivx+ndimV-1))
          rmsval = rmsval + voli*v2i
          rmsvmw = rmsvmw + v2i
       endif
    enddo
    !
    !--use the computed volume for velocity, otherwise won't be normalised correctly
    !
    rmsval = sqrt(rmsval/totvol)
    rmsvmw = sqrt(rmsvmw/dble(ntot))   

    !
    !--option to override volume from sum with environment variable
    !
    voli = renvironment('SPLASH_CALC_VOLUME',errval=-1.0)
    if (voli.gt.0.) then
       print "(1x,a,es9.2)",&
                  'volume from sum(m/rho)           = ',totvol
       totvol = voli
       print "(1x,a,es9.2,/)",&
                  '**overridden with SPLASH_CALC_VOLUME = ',totvol
    else
       print "(1x,a,es9.2,/,1x,a,/)",&
                  'volume from sum(m/rho)           = ',totvol,&
                  '(override this using SPLASH_CALC_VOLUME environment variable)'
    endif

    rhomeanmw = rhomeanmw/real(ntot)
    rhomeanvw = rhomeanvw/totvol
    smeanmw   = smeanmw/real(ntot)
    smeanvw   = smeanvw/totvol
    !
    !--calculate variance on second pass
    !
    rhovarvw = 0.
    rhovarmw = 0.
    svarvw   = 0.
    svarmw   = 0.
    totvol   = 0.
    do i=1,ntot
       itype  = igettype(i)
       pmassi = particlemass(i,itype)
       rhoi   = dat(i,irho)
       if (rhoi.gt.0.) then
          voli = pmassi/rhoi
          si = log(rhoi)
       else
          voli = 0.
          si   = 0.
       endif
       totvol   = totvol + voli
       rhovarmw = rhovarmw + (rhoi - rhomeanmw)**2
       rhovarvw = rhovarvw + voli*(rhoi - rhomeanvw)**2
       svarmw   = svarmw   + (si - smeanmw)**2
       svarvw   = svarvw   + voli*(si - smeanvw)**2
    enddo
    rhovarmw = rhovarmw/real(ntot)
    rhovarvw = rhovarvw/totvol
    svarmw   = svarmw/real(ntot)
    svarvw   = svarvw/totvol

    !
    !--write output to screen/terminal
    !
    print "(1x,'mean density (vol. weighted)     = ',es11.4,' +/- ',es11.4)",rhomeanvw,sqrt(rhovarvw)
    print "(1x,'mean density (mass weighted)     = ',es11.4,' +/- ',es11.4)",rhomeanmw,sqrt(rhovarmw)
    print "(1x,'density variance (vol. weighted) = ',es11.4)",rhovarvw
    print "(1x,'density variance (mass weighted) = ',es11.4)",rhovarmw
    print "(1x,'mean ln density (vol. weighted)     = ',es11.4,' +/- ',es11.4)",smeanvw,sqrt(svarvw)
    print "(1x,'           -0.5*var(ln density)     = ',es11.4)",-0.5*svarvw
    print "(1x,'mean ln density (mass weighted)     = ',es11.4,' +/- ',es11.4)",smeanmw,sqrt(svarmw)
    print "(1x,'ln density variance (vol. weighted) = ',es11.4)",svarvw
    print "(1x,'ln density variance (mass weighted) = ',es11.4)",svarmw
    print "(1x,'rms velocity     (vol. weighted) = ',es11.4)",rmsval
    print "(1x,'rms velocity     (mass weighted) = ',es11.4)",rmsvmw
    if (rmsval.gt.0.) then
       bval = sqrt(svarvw)/rmsval
    else
       bval = 0.
    endif
    if (rmsvmw.gt.0.) then
       bvalmw = sqrt(svarmw)/rmsvmw
    else
       bvalmw = 0.
    endif    
    print "(1x,'sqrt(sigma^2/v^2)(vol. weighted) = ',f11.3)",bval
    print "(1x,'sqrt(sigma^2/v^2)(mass weighted) = ',f11.3)",bvalmw
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) 17
    write(iunit,fmtstring) time,rhomeanvw,rhomeanmw,rhovarvw,rhovarmw,sqrt(rhovarvw),sqrt(rhovarmw),&
                           rmsval,rmsvmw,bval,bvalmw,smeanvw,smeanmw,svarvw,svarmw,sqrt(svarvw),sqrt(svarmw)

 case('timeaverage','timeav')
    if (.not.allocated(datmean)) then
       allocate(datmean(ntot,ncolumns),stat=ierr)
       if (ierr /= 0) stop 'error allocating memory for mean sum in calc'
       datmean = 0.
    endif
    if (.not.allocated(datvar)) then
       allocate(datvar(ntot,ncolumns),stat=ierr)
       if (ierr /= 0) stop 'error allocating memory for variance sum in calc'
       datvar = 0.
    endif
    ntot1 = size(datmean(:,1))
    if (ntot.gt.ntot1) then
       print*,' WARNING: nrows = ',ntot,' > nrows from previous dumpfile =',ntot1
       print*,'          ignoring all rows/particles greater than ',ntot1
    elseif (ntot.lt.ntot1) then
       print*,' WARNING: nrows = ',ntot,' < nrows from previous dumpfile =',ntot1
       print*,'          assuming zeros for rows/particles greater than ',ntot
    endif
    ncol1 = size(datmean(1,:))
    if (ncolumns.gt.ncol1) then
       print*,' WARNING: ncolumns = ',ncolumns,' > ncolumns from previous dumpfile =',ncol1
       print*,'          ignoring all rows/particles greater than ',ncol1
    elseif (ncolumns.lt.ncol1) then
       print*,' WARNING: ncolumns = ',ntot,' < ncolumns from previous dumpfile =',ncol1
       print*,'          assuming zeros for columns greater than ',ncolumns
    endif
    ntot1 = min(ntot1,ntot)
    ncol1 = min(ncol1,ncolumns)
    
    dn = 1./real(nfilesread)
    !
    !--compute the mean and variance using Knuth/Welford's compensated
    !  summation algorithm to minimise round-off error
    !  (see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
    !
    do j=1,ncol1
       do i=1,ntot1
          delta = dat(i,j) - datmean(i,j)
          datmean(i,j) = datmean(i,j) + delta*dn
          datvar(i,j) = datvar(i,j) + delta*(dat(i,j) - datmean(i,j))
       enddo
    enddo
    return
    !sum1(:,:) = sum1(:,:) + dat(1:ntot1,1:ncol1)
    !sum2(:,:) = sum2(:,:) + dat(1:ntot1,1:ncol1)**2

 case default
    print "(a)",' ERROR: unknown analysis type in write_analysis routine'
    return 
 end select

 print "(/,1x,'>>> ',a,' <<<')",'written to '//trim(fileout)

 return

contains
!
!--small internal utility to work out the particle type
!  (which depends whether or not mixed types are stored)
!
integer function igettype(i)
 implicit none
 integer :: np
 integer, intent(in) :: i

 if (size(iamtype).gt.1) then
    igettype = int(iamtype(i))
 else
    np = 0
    igettype = 0
    do while (i.gt.np .and. igettype.le.ntypes)
       igettype = igettype + 1
       np = np + npartoftype(igettype)
    enddo
 endif       

end function igettype
!
!--small internal utility to get particle mass
!  (depends on whether or not mass is stored for each particle
!   or only for each type)
!
real function particlemass(i,iparttype)
 implicit none
 integer, intent(in) :: i,iparttype

 if (ipmass.gt.0 .and. ipmass.le.ncolumns) then
    particlemass = dat(i,ipmass)
 else
    particlemass = massoftype(iparttype)
 endif

end function particlemass

end subroutine write_analysis

!---------------------
! close output file
!---------------------
subroutine close_analysis(analysistype)
 implicit none
 character(len=*), intent(in) :: analysistype
 integer :: i

 select case(trim(analysistype))
 case('timeaverage','timeav')
    print "(a)",'----> WRITING time_average.out ...'
    if (allocated(datmean) .and. allocated(datvar) .and. nfilesread.gt.0) then
       !--get standard deviation from variance (also normalise with 1/n)
       datvar(:,:) = sqrt(datvar(:,:))/sqrt(real(nfilesread))
       do i=1,size(datmean(:,1))
          write(iunit,"(1x,99(es15.7,2x))") datmean(i,:),datvar(i,:)
       enddo
    endif
 end select

 if (allocated(datmean)) deallocate(datmean)
 if (allocated(datvar)) deallocate(datvar)

 close(unit=iunit)

 return
end subroutine close_analysis

end module analysis
