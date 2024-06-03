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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------
!  module implementing the ability to use SPLASH to produce
!  evolution files from a sequence of SPH dump files
!  (ie. in order to produce plots of certain quantities vs time)
!  Command is "splash calc X" where X is analysis type.
!-------------------------------------------------------------------
module analysis
 implicit none
 public :: isanalysis,open_analysis,close_analysis,write_analysis
 private
 integer, private, parameter :: iunit = 89
 integer, private, parameter :: maxlevels = 20
 integer, private, parameter :: maxtrack = 32
!
! default settings for the density thresholds for massaboverho output
!
 integer, private :: nlevels,nfilesread,nfileout
 real, dimension(maxlevels), private :: rholevels
 character(len=64), private :: fileout(maxtrack)
 real, dimension(:,:), allocatable :: datmean,datvar
 integer, private :: itracks(maxtrack)

 !
 ! definitions needed for tdiffuse output
 !
 integer, parameter :: ndirs = 6
 character(len=2) :: dir_label(ndirs) = (/'-x','+x','-y','+y','-z','+z'/)

 ! angles required for direction:   -x  +x  -y  +y  -z  +z
 real, parameter :: anglex_vals(ndirs) = (/0.,    0., -90., 90., 180., 0./)
 real, parameter :: angley_vals(ndirs) = (/-90., 90.,   0.,  0.,   0., 0./)
 real, parameter :: anglez_vals(ndirs) = (/0.,    0.,   0.,  0.,   0., 0./)

contains

!-----------------------------------------------------------------
! utility to check if the choice of analysis type is valid
! and if not, to print the available options
!-----------------------------------------------------------------
logical function isanalysis(string,noprint)
 character(len=*), intent(in) :: string
 logical, intent(in), optional :: noprint
 logical :: doprint,verbose

 isanalysis = .false.
 verbose = .true.
 select case(trim(string))
 case('energies','energy')
    isanalysis = .true.
 case('massaboverho')
    isanalysis = .true.
 case('max','maxvals')
    isanalysis = .true.
 case('min','minvals')
    isanalysis = .true.
 case('diff','diffvals')
    isanalysis = .true.
 case('delta','deltavals')
    isanalysis = .true.
 case('amp','ampvals')
    isanalysis = .true.
 case('mean','meanvals')
    isanalysis = .true.
 case('rms','rmsvals')
    isanalysis = .true.
 case('vrms','vrmsvals','vwrms','rmsvw')
    isanalysis = .true.
 case('rhovar','rhomach')
    isanalysis = .true.
 case('kh')
    isanalysis = .true.
 case('timeaverage','timeav')
    isanalysis = .true.
 case('ratio','minus','plus')
    isanalysis = .true.
 case('tracks','track')
    isanalysis = .true.
 case('lightcurve')
    isanalysis = .true.
 case('extinction')
    isanalysis = .true.
 case('tdiffuse')
    isanalysis = .true.
 case('none')
    verbose = .false.
 end select

 if (present(noprint)) then
    doprint = .not.noprint
 else
    doprint = .true.
 endif

 if (.not.isanalysis .and. doprint .and. verbose) then
    print "(a)",' Analysis mode ("splash calc X dumpfiles") on a sequence of dump files: '
    print "(a)",'  splash calc energies     : calculate KE,PE,total energy vs time'
    print "(a)",'                             output to file called ''energy.out'''
    print "(a)",'         calc massaboverho : mass above a series of density thresholds vs time'
    print "(a)",'                             output to file called ''massaboverho.out'''
    print "(a)",'         calc extinction   : column density to all sink particles vs time'
    print "(a)",'                             output to file called ''extinction.out'''
!    print "(a)",'         calc rhomach      : density variance and RMS velocity dispersion vs. time'
!    print "(a)",'                             output to file called ''rhomach.out'''
    print "(a)",'         calc max          : maximum of each column vs. time'
    print "(a)",'                             output to file called ''maxvals.out'''
    print "(a)",'         calc min          : minimum of each column vs. time'
    print "(a)",'                             output to file called ''minvals.out'''
    print "(a)",'         calc diff         : (max - min) of each column vs. time'
    print "(a)",'                             output to file called ''diffvals.out'''
    print "(a)",'         calc amp          : 0.5*(max - min) of each column vs. time'
    print "(a)",'                             output to file called ''ampvals.out'''
    print "(a)",'         calc delta        : 0.5*(max - min)/mean of each column vs. time'
    print "(a)",'                             output to file called ''deltavals.out'''
    print "(a)",'         calc mean         : mean of each column vs. time'
    print "(a)",'                             output to file called ''meanvals.out'''
    print "(a)",'         calc rms          : (mass weighted) root mean square of each column vs. time'
    print "(a)",'                             output to file called ''rmsvals.out'''
    print "(a)",'         calc tracks       : track particle data vs time for selected particles,'
    print "(a)",'           --track=1,2,3    output to tracks-1.out,tracks-2.out,tracks-3.out'
!    print "(a)",'         calc vrms         : volume weighted root mean square of each column vs. time'
!    print "(a)",'                             output to file called ''rmsvals-vw.out'''
    print "(/,a)",'  the above options all produce a small ascii file with one row per input file.'
    print "(a)",'  the following option produces a file equivalent in size to one input file (in ascii format):'
    print "(/,a)",'         calc timeaverage  : time average of *all* entries for every particle'
    print "(a)",'                             output to file called ''time_average.out'''
    print "(/,a)",'         calc ratio        : ratio of *all* entries in each file compared to first'
    print "(a)",'                             output to file called ''ratio.out'''
    print "(/,a)",'         calc plus        : add two snapshots together'
    print "(a)",'                             output to file called ''plus.out'''
 elseif (.not.isanalysis .and. doprint) then
    print "(a)",'Analysis mode:'
    print "(a)",'   splash calc <mode>  : type "splash calc" for details'
 endif

 return
end function isanalysis

!----------------------------------------------------------------
!  open output file/ initialise quantities needed for analysis
!  over all dump files
!----------------------------------------------------------------
subroutine open_analysis(analysistype,required,ncolumns,ndim,ndimV,nsinks)
 use labels,        only:ix,ivx,ih,iBfirst,iutherm,irho,ipmass,itemp,ikappa,label,&
                         lenunitslabel,unitslabel,labelzintegration,get_unitlabel_coldens
 use asciiutils,    only:read_asciifile,basename,integer_to_string
 use filenames,     only:rootname,nfiles,tagline,fileprefix,ifileopen
 use params,        only:maxplot
 use system_utils,  only:ienvlist
 use settings_data, only:iRescale
 integer, intent(in) :: ncolumns,ndim,ndimV,nsinks
 character(len=*), intent(in) :: analysistype
 logical, dimension(0:ncolumns), intent(out) :: required
 character(len=maxplot*18) :: headerline   ! len=maxplot x 18 characters
 character(len=64) :: levelsfile
 character(len=maxplot*12) :: fmtstring
 character(len=lenunitslabel) :: labelt,labelc
 logical :: iexist,standardheader
 integer :: ierr,i,lunit
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
 nfileout = 1
 labelt = ''
 if (iRescale) then
    labelt = unitslabel(0)
 endif

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
    if (iBfirst > 0) required(irho) = .true.
    required(ix(1:ndim)) = .true.
    !
    !--set filename and header line
    !
    fileout(1) = 'energy.out'
    write(headerline,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time',2,'ekin',3,'etherm',4,'emag',5,'epot',6,'etot',7,'totmom',8,'totang'

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
          print "(a,i2,a,es9.2)",' level ',i,': rho > ',rholevels(i)
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
          rholevels(1:nlevels) = (/1e-20,1e-19,1e-18,1e-17,1e-16, &
                                   1e-15,1e-14,1e-13,1e-12,1e-11/)
          write(iunit+1,*) rholevels(1:nlevels)
          close(iunit+1)
       endif
       stop
    endif
    !
    !--set filename and header line
    !
    fileout(1) = 'massaboverho.out'
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
    fileout(1) = 'maxvals.out'
    standardheader = .true.

 case('min','minvals')

    required(:) = .true.
    fileout(1) = 'minvals.out'
    standardheader = .true.

 case('diff','diffvals')

    required(:) = .true.
    fileout(1) = 'diffvals.out'
    standardheader = .true.

 case('amp','ampvals')

    required(:) = .true.
    fileout(1) = 'ampvals.out'
    standardheader = .true.

 case('delta','deltavals','deltas')

    required(:) = .true.
    fileout(1) = 'deltavals.out'
    standardheader = .true.

 case('mean','meanvals')

    required(:) = .true.
    fileout(1) = 'meanvals.out'
    standardheader = .true.

 case('rms','rmsvals')

    required(:) = .true.
    fileout(1) = 'rmsvals.out'
    standardheader = .true.

 case('vrms','vrmsvals','vwrms','rmsvw')

    required(:) = .true.
    fileout(1) = 'rmsvals-vw.out'
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
    fileout(1) = 'rhomach.out'
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) 17
    write(headerline,fmtstring) 1,'time',2,'rhomean(vw)',3,'rhomean(mw)',4,'varrho(vw)',5,'varrho(mw)',&
          6,'stddevrho(vw)',7,'stddevrho(mw)',8,'rms v (vw)',9,'rms v (mw)',10,'b (vw)',11,'b (mw)',&
          12,'s mean(vw)',13,'s mean(mw)',14,'s var(vw)',15,'s var(mw)',16,'s stddev(vw)',17,'s stddev(mw)'

 case('kh')
    !
    !--read all columns from dump file
    !
    required(irho) = .true.
    required(ivx:ivx+ndimV-1) = .true.
    !
    !--set filename and header line
    !
    fileout(1) = 'kh.out'
    standardheader = .true.
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) 2
    write(headerline,fmtstring) 1,'time',2,'max(ekiny)'

 case('timeaverage','timeav')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout(1) = 'time_average.out'
    if (ncolumns > 0) then
       write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']''))')",iostat=ierr) 2*ncolumns
       write(headerline,fmtstring,iostat=ierr) (i,label(i)(1:12),i=1,ncolumns),&
                                               (ncolumns+i,'err'//label(i)(1:9),i=1,ncolumns)
    endif
 case('ratio','minus','add')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout(1) = trim(analysistype)//'.out'
    if (ncolumns > 0 .and. ncolumns /= maxplot) then
       write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']''))')",iostat=ierr) 2*ncolumns
       write(headerline,fmtstring,iostat=ierr) (i,label(i)(1:12),i=1,ncolumns),&
                                               (ncolumns+i,'err'//label(i)(1:9),i=1,ncolumns)
    endif
 case('tracks','track')
    !
    !--read all of dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout(1) = 'tracks.out' ! in case ntracks = 0
    !
    !--look for a file "splash.tracks" for particle ids to track
    !
    call read_asciifile(trim(fileprefix)//'.tracks',nfileout,itracks,ierr)
    if (ierr /= 0) then
       !
       !--otherwise use --tracks=1,2,3
       !
       itracks = ienvlist('SPLASH_TRACK',size(itracks))
       nfileout = count(itracks > 0)
    endif
    if (nfileout > 0) then
       print "(a,i0,a)",' TRACKING ',nfileout,' PARTICLES'
    else
       print "(a)",' Use --track=1,2,3 or splash.tracks file to specify particle IDs to track'
    endif
    do i=1,nfileout
       write(fileout(i),"(a,i0,a)") 'tracks-',itracks(i),'.out'
    enddo
    nfileout = 1
    standardheader = .true.

 case('lightcurve')
    !
    !--for lightcurve need to h, mass and utherm
    !
    required(ix(1:ndim)) = .true.
    required(ih) = .true.
    required(iutherm) = .true.
    required(ipmass) = .true.
    required(itemp) = .true.
    required(ikappa) = .true.
    !
    !--set filename and header line
    !
    if (nfiles==1) then
       fileout = 'lightcurve_'//trim(basename(rootname(ifileopen)))//'.out'
    else
       fileout = 'lightcurve.out'
    endif
    write(headerline,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
           1,'time'//trim(labelt),2,'Luminosity',3,'R_{eff}',4,'T_{eff}',&
           5,'L_{bol}',6,'R_{bb}',7,'T_c',8,'bad pix %'

 case('extinction')
    !
    !--for extinction of stars need h, mass and density
    !
    required(ix(1:ndim)) = .true.
    required(ih) = .true.
    required(ipmass) = .true.
    required(irho) = .true.
    !
    !--set filename and header line
    !
    if (nfiles==1) then
       fileout = 'extinction_'//trim(basename(rootname(ifileopen)))//'.out'
    else
       fileout = 'extinction.out'
    endif
    labelc = get_unitlabel_coldens(iRescale,labelzintegration,unitslabel(irho))

    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a15,'']'',2x))')",iostat=ierr) nsinks+1
    write(headerline,fmtstring) &
           1,'time'//trim(labelt),(i+1,'sink'//trim(adjustl(integer_to_string(i)))//trim(labelc),i=1,nsinks)
 case('tdiffuse')
    !
    !--similar to extinction, but could be any column
    !
    !required(ix(1:ndim)) = .true.
    !required(ih) = .true.
    !required(ipmass) = .true.
    !required(irho) = .true.
    required(:) = .true.
    !
    !--set filename and header line
    !
    if (nfiles==1) then
       fileout = 'tdiffuse_'//trim(basename(rootname(ifileopen)))//'.out'
    else
       fileout = 'tdiffuse.out'
    endif

    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a15,'']'',2x))')",iostat=ierr) ndirs+1
    write(headerline,fmtstring) &
           1,'time'//trim(labelt),(i+1,dir_label(i),i=1,ndirs)
 end select

 if (standardheader) then
!
!--standard header is time in column 1, with an entry for each column following
!  (this is to avoid repeated code above)
!
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) ncolumns+1
    write(headerline,fmtstring) 1,'time'//trim(labelt),(i+1,label(i)(1:12),i=1,ncolumns)
 endif

!
!--do not replace the file if it already exists
!
 do i=1,nfileout
    inquire(file=trim(fileout(i)),exist=iexist)
    if (iexist) then
       print "(2(a,/))",' ERROR: analysis file '//trim(fileout(i))//' already exists', &
                        '        delete, move or rename this file and try again'
       stop
    endif
    !
    !--open the file for output
    !
    lunit = iunit+i-1
    if (nfileout > 1) then
       open(unit=lunit,file=trim(fileout(i)),status='replace',form='formatted',iostat=ierr)
    else
       open(unit=lunit,file=trim(fileout(i)),status='new',form='formatted',iostat=ierr)
    endif
    if (ierr /= 0) then
       print "(a)",' ERROR opening file '//trim(fileout(i))//' for output'
       stop
    endif
    print "(a)",' WRITING '//trim(analysistype)//' vs time to file '//trim(fileout(i))
    !
    !--write header if the headerline is set
    !  (no header is written if headerline is blank)
    !
    if (len_trim(headerline) > 0) then
       write(lunit,"(a)") '# '//trim(tagline)
       write(lunit,"(a)") '# '//trim(fileout(i))//' produced using "splash calc '//trim(analysistype)// &
                          '" on dump files '//trim(rootname(1))//'->'//trim(rootname(nfiles))
       write(lunit,"(a)") '# use splash -ev '//trim(fileout(i))//' to plot the contents of this file '
       write(lunit,"(a)") '#'
       write(lunit,"(a)") trim(headerline)
    endif
 enddo
 !if (nfileout == 0) then
 !print "(a)",' ERROR: no output from analysis, missing options?'
 !stop
 !endif
 nfilesread = 0

end subroutine open_analysis

!----------------------------------------------------------------
!  this is the routine which actually calculates the quantities
!  required from each dump file and spits out a line to the
!  analysis file. Called once for each dump file.
!----------------------------------------------------------------
subroutine write_analysis(time,dat,ntot,ntypes,npartoftype,massoftype,&
                          iamtype,ncolumns,ndim,ndimV,analysistype)
 use labels,        only:ix,ivx,iBfirst,iutherm,irho,ipmass,labeltype,label,get_sink_type,&
                         lenunitslabel,get_unitlabel_coldens,unitslabel,labelzintegration
 use params,        only:int1,doub_prec,maxplot
 use asciiutils,    only:ucase,basename
 use system_utils,  only:renvironment,ienvironment
 use settings_part, only:iplotpartoftype
 use particle_data, only:time_was_read
 use part_utils,    only:get_positions_of_type
 use settings_data, only:xorigin,icoords,icoordsnew,track_string,iRescale
 use geomutils,     only:change_coords
 use part_utils,    only:get_tracked_particle
 use lightcurve,    only:get_lightcurve
 use extinction,    only:get_extinction,get_extinction_los
 use filenames,     only:rootname,ifileopen
 use vectorutils,   only:cross_product3D
 integer, intent(in)               :: ntot,ntypes,ncolumns,ndim,ndimV
 integer, intent(in), dimension(:) :: npartoftype
 real, intent(in), dimension(:)    :: massoftype
 integer(kind=int1), intent(in), dimension(:) :: iamtype
 real, intent(in)                  :: time
 real, intent(in), dimension(:,:)  :: dat
 character(len=*), intent(in)      :: analysistype
 real(kind=doub_prec), dimension(maxlevels) :: massaboverho
 integer              :: itype,i,j,ierr,ntot1,ncol1,nused,itrack,ifile,npts,isinktype,icol
 real(kind=doub_prec) :: ekin,emag,etherm,epot,etot,totmom,pmassi,totang
 real(kind=doub_prec) :: totvol,voli,rhoi,rmsvmw,v2i
 real(kind=doub_prec) :: rhomeanmw,rhomeanvw,rhovarmw,rhovarvw,bval,bvalmw
 real(kind=doub_prec) :: smeanmw,smeanvw,svarmw,svarvw,si,ekiny,ekinymax
 real(kind=doub_prec) :: lmin(maxplot),lmax(maxplot),lmean(maxplot),rmsvali
 real(kind=doub_prec), dimension(3) :: xmom,angmom,angmomi,ri,vi
 real                 :: delta,dn,valmin,valmax,valmean,timei
 real                 :: lum,rphoto,tphoto,l_bb,r_bb,t_bb,badfrac
 character(len=20)    :: fmtstring
 logical              :: change_coordsys
 real                 :: x0(3),v0(3)
 real, allocatable    :: xpts(:),ypts(:),zpts(:)
 character(len=lenunitslabel) :: labelt,labelc
!
! array with one value for each column
!
 real(kind=doub_prec) :: coltemp(maxplot), vals(maxplot), rmsval(maxplot)
 real :: coltemps(maxplot)

 labelt = ''
 labelc = ''
 if (iRescale) labelt = unitslabel(0)

 nfilesread = nfilesread + 1
 if (time_was_read(time)) then
    timei = time
    print "(/,5('-'),a,', TIME=',es9.2,a,' FILE #',i5,/)",&
          '> CALCULATING '//trim(ucase(analysistype)),time,trim(labelt),nfilesread
 else
    timei = 0.
    print "(/,5('-'),a,', FILE #',i5,' (TIME NOT READ)'/)",&
          '> CALCULATING '//trim(ucase(analysistype)),nfilesread
 endif

 change_coordsys = (icoordsnew /= icoords .and. ndim > 0 .and. all(ix(1:ndim) > 0))
 x0 = xorigin(:)  ! note that it is not currently possible to do splash to ascii
 v0 = 0.          ! with coords set relative to a tracked particle, so just use xorigin
 ! instead, one can use the --origin flag to make positions and velocities relative to a particle

 if (itracks(1) > 0) then
    itrack = itracks(1)  ! override particle id saved to splash.defaults file if --tracks specified
 else
    itrack = get_tracked_particle(track_string,npartoftype,iamtype,dat,irho)
 endif
 if (itrack==0) itrack = 1

 select case(trim(analysistype))
 case('energy','energies')
    ekin = 0.
    emag = 0.
    epot = 0.
    etherm = 0.
    etot = 0.
    xmom = 0.
    angmom = 0.
    nused = 0
    do i=1,ntot
       itype = igettype(i)
       if (iplotpartoftype(itype)) then
          pmassi = particlemass(i,itype)

          !--kinetic energy
          if (ivx > 0 .and. ivx+ndimV-1 <= ncolumns) then
             vi(:) = 0.
             vi(1:ndimV) = dat(i,ivx:ivx+ndimV-1)

             ekin = ekin + pmassi*dot_product(vi,vi)

             !--linear momentum
             xmom = xmom + pmassi*vi

             !--angular momentum
             if (ndim >= 1 .and. all(ix(1:ndim) > 0)) then
                ri(:) = 0.
                ri(1:ndim) = dat(i,ix(1):ix(ndim))
                call cross_product3D(ri,vi,angmomi)
                angmom(:) = angmom(:) + pmassi*angmomi(:)
             endif
          endif

          !--thermal energy
          if (iutherm > 0 .and. iutherm <= ncolumns) then
             etherm = etherm + pmassi*dat(i,iutherm)
          endif

          !--magnetic energy
          if (iBfirst > 0 .and. iBfirst+ndimV-1 <= ncolumns) then
             emag = emag + pmassi*dot_product(dat(i,iBfirst:iBfirst+ndimV-1),&
                                              dat(i,iBfirst:iBfirst+ndimV-1))/dat(i,irho)
          endif
          nused = nused + 1
       endif
    enddo
    ekin = 0.5*ekin
    emag = 0.5*emag
    etot = ekin + etherm + epot + emag
    totmom = sqrt(dot_product(xmom(1:ndimV),xmom(1:ndimV)))
    totang = sqrt(dot_product(angmom,angmom))

    print "(7(/,1x,a6,' = ',es9.2))",'etot',etot,'ekin',ekin,'etherm',etherm,'epot',epot,'emag',emag,'totmom',totmom,'totang',totang
    !
    !--write line to output file
    !
    write(iunit,"(64(es18.10,1x))") timei,ekin,etherm,emag,epot,etot,totmom,totang
    if (nused /= ntot) print*,'energies calculated using ',nused,' of ',ntot,' particles'

 case('massaboverho')
    massaboverho(:) = 0.
    if (irho > 0 .and. irho <= ncolumns) then
       !
       !--warn if particle masses not found
       !
       if (ipmass <= 0 .or. ipmass > ncolumns .and. all(massoftype < tiny(massoftype))) then
          print "(a)",' WARNING in massaboverho analysis!'// &
                      ' masses not read or are zero from dump file'
       endif
       !
       !--calculate mass above each density threshold
       !
       do i=1,ntot
          itype = igettype(i)
          pmassi = particlemass(i,itype)

          if (itype==1) then
             !
             !--gas particles contribute if they are above rho
             !
             where(dat(i,irho) >= rholevels(1:nlevels))
                massaboverho(1:nlevels) = massaboverho(1:nlevels) + pmassi
             end where
          elseif (labeltype(itype)=='sink') then
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
          print "(1x,'M(rho > ',es9.2,') = ',es9.2)",rholevels(i),massaboverho(i)
       enddo
       !
       !--write line to output file
       !
       write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) nlevels+1
       write(iunit,fmtstring) timei,massaboverho(1:nlevels)

    else
       print "(a)",' ERROR in massaboverho analysis!'// &
                   ' either mass or density not found in dump file'
       return
    endif
 case('max','maxvals')
    !
    !--calculate maximum for each column
    !
    coltemp(:) = -huge(0.d0) !maxval(dat(1:ntot,i))
    nused = 0
    do j=1,ntot
       itype = igettype(j)
       if (iplotpartoftype(itype)) then
          vals(1:ncolumns) = real(dat(j,1:ncolumns),kind=doub_prec)
          if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
          nused = nused + 1
          do i=1,ncolumns
             coltemp(i) = max(coltemp(i),vals(i))
          enddo
       endif
    enddo
    where (coltemp(:) < -0.5*huge(0.)) coltemp(:) = 0.
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'max = ',es18.10)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) timei,coltemp(1:ncolumns)
    if (nused /= ntot) print*,'max calculated using ',nused,' of ',ntot,' particles'

 case('min','minvals')
    !
    !--calculate minimum for each column
    !
    coltemp(:) = huge(0.d0) !minval(dat(1:ntot,i))
    nused = 0
    do j=1,ntot
       itype = igettype(j)
       if (iplotpartoftype(itype)) then
          vals(1:ncolumns) = real(dat(j,1:ncolumns),kind=doub_prec)
          if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
          nused = nused + 1
          do i=1,ncolumns
             coltemp(i) = min(coltemp(i),vals(i))
          enddo
       endif
    enddo
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'min = ',es18.10)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) timei,coltemp(1:ncolumns)
    if (nused /= ntot) print*,'min calculated using ',nused,' of ',ntot,' particles'

 case('diff','diffvals','amp','ampvals','delta','deltavals','deltas','tracks')
    !
    !--calculate difference between max and min for each column
    !
    lmean(:) = 0.
    lmin(:) = huge(0.d0)
    lmax(:) = -huge(0.d0)
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1

    do ifile=1,nfileout
       if (nfileout > 1) itrack = itracks(ifile)
       nused = 0
       coltemp = 0.
       do j=1,ntot
          itype = igettype(j)
          if (iplotpartoftype(itype) .or. j==itrack) then
             vals(1:ncolumns) = real(dat(j,1:ncolumns),kind=doub_prec)
             if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
             nused = nused + 1
             do i=1,ncolumns
                lmin(i) = min(lmin(i), vals(i))
                lmax(i) = max(lmax(i), vals(i))
                lmean(i) = lmean(i) + vals(i)
             enddo
             if (j==itrack) coltemp = vals
          endif
       enddo

       if (trim(analysistype)=='tracks' .or. trim(analysistype)=='track') then
          if (ifile==1) then
             do i=1,ncolumns
                print "(1x,' particle ',i8,': ',a20,' = ',es18.10)",itrack,label(i),coltemp(i)
             enddo
          endif
          !
          !--write line to output file
          !
          write(iunit+ifile-1,fmtstring) timei,coltemp(1:ncolumns)
          if (nfileout > 1) print "(1x,'>>> ',a,' <<<')",'written to '//trim(fileout(ifile))
       endif
    enddo
    if (nused > 0) lmean(:) = lmean(:)/real(nused)

    select case(trim(analysistype))
    case('amp','ampvals')
       do i=1,ncolumns
          coltemp(i) = 0.5*(lmax(i) - lmin(i))
          print "(1x,a20,'0.5*(max - min) = ',es18.10)",label(i),coltemp(i)
       enddo
    case('delta','deltavals','deltas')
       do i=1,ncolumns
          valmean = real(lmean(i))
          if (valmean > 0.) then
             coltemp(i) = 0.5*(lmax(i) - lmin(i))/valmean
          else
             coltemp(i) = 0.5*(lmax(i) - lmin(i))
          endif
          print "(1x,a20,'0.5*(max - min)/mean = ',es18.10)",label(i),coltemp(i)
       enddo
    case('diff','diffvals') ! diff, diffvals
       do i=1,ncolumns
          coltemp(i) = lmax(i) - lmin(i)
          print "(1x,a20,'(max - min) = ',es18.10)",label(i),coltemp(i)
       enddo
    end select
    !
    !--write line to output file
    !
    if (trim(analysistype) /= 'tracks' .and. trim(analysistype) /= 'track') then
       write(iunit,fmtstring) timei,coltemp(1:ncolumns)
    endif

    if (nused /= ntot) then
       select case(trim(analysistype))
       case('diff', 'diffvals')
          print*,'diff calculated using ',nused,' of ',ntot,' particles'
       case('amp','ampvals')
          print*,'amp calculated using ',nused,' of ',ntot,' particles'
       case('delta','deltavals','deltas')
          print*,'deltas calculated using ',nused,' of ',ntot,' particles'
       end select
    endif

 case('mean','meanvals')
    !
    !--calculate mean for each column
    !
    coltemp(:) = 0.
    nused = 0
    do j=1,ntot
       itype = igettype(j)
       if (iplotpartoftype(itype)) then
          vals(1:ncolumns) = real(dat(j,1:ncolumns),kind=doub_prec)
          if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
          nused = nused + 1
          do i=1,ncolumns
             coltemp(i) = coltemp(i) + vals(i)
          enddo
       endif
    enddo
    if (nused > 0) then
       coltemp(:) = coltemp(:)/real(nused)
    else
       coltemp(:) = 0.
    endif
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'mean = ',es18.10)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) timei,coltemp(1:ncolumns)
    if (nused /= ntot) print*,'mean calculated using ',nused,' of ',ntot,' particles'

 case('rms','rmsvals')
    !
    !--calculate RMS for each column
    !
    coltemp(:) = 0.
    nused = 0
    do j=1,ntot
       itype = igettype(j)
       if (iplotpartoftype(itype)) then
          vals(1:ncolumns) = real(dat(j,1:ncolumns),kind=doub_prec)
          if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
          nused = nused + 1
          do i=1,ncolumns
             coltemp(i) = coltemp(i) + vals(i)**2
          enddo
       endif
    enddo
    if (nused > 0) then
       coltemp(:) = sqrt(coltemp(:)/real(nused))
    else
       coltemp(:) = 0.
    endif
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'rms (mass weighted) = ',es18.10)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) timei,coltemp(1:ncolumns)
    if (nused /= ntot) print*,'rms calculated using ',nused,' of ',ntot,' particles'

 case('vrms','vrmsvals','vwrms','rmsvw')
    if (irho <= 0 .or. irho > ncolumns) then
       print "(a)",' ERROR in volume weighted rms calculation!'// &
                   ' density not present / not labelled in dump file, skipping...'
       return
    endif
    !
    !--warn if particle masses not found
    !
    if (ipmass <= 0 .or. ipmass > ncolumns .and. all(massoftype < tiny(massoftype))) then
       print "(a)",' WARNING in volume weighted rms calculation!'// &
                   ' masses not read or are zero from dump file'
    endif
    !
    !--calculate volume-weighted RMS for each column
    !
    rmsval(:) = 0.
    totvol = 0.
    do j=1,ntot
       itype  = igettype(j)
       if (iplotpartoftype(itype)) then
          vals(1:ncolumns) = real(dat(j,1:ncolumns),kind=doub_prec)
          if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)

          pmassi = particlemass(j,itype)
          rhoi   = dat(j,irho)
          if (rhoi > 0.) then
             voli = pmassi/rhoi
          else
             voli = 0.
          endif

          do i=1,ncolumns
             rmsval(i) = rmsval(i) + voli*vals(i)**2
          enddo
          totvol = totvol + voli
       endif
    enddo
    coltemp(:) = real(sqrt(rmsval(:)/totvol))
    print "(1x,a,es9.2)",'volume = ',totvol
    !
    !--write output to screen/terminal
    !
    do i=1,ncolumns
       print "(1x,a20,'rms (volume weighted) = ',es18.10)",label(i),coltemp(i)
    enddo
    !
    !--write line to output file
    !
    write(fmtstring,"('(',i3,'(es18.10,1x))')",iostat=ierr) ncolumns+1
    write(iunit,fmtstring) timei,coltemp(1:ncolumns)

 case('rhovar','rhomach')
    if (irho <= 0 .or. irho > ncolumns) then
       print "(a)",' ERROR in density variance--rms velocity field calculation!'// &
                   ' density not present / not labelled in dump file, skipping...'
       return
    endif
    !
    !--warn if particle masses not found
    !
    if (ipmass <= 0 .or. ipmass > ncolumns .and. all(massoftype < tiny(massoftype))) then
       print "(a)",' WARNING in volume weighted rms calculation!'// &
                   ' masses not read or are zero from dump file'
    endif

    if (ivx <= 0 .or. ivx > ncolumns) then
       print "(a)",' WARNING in volume weighted rms calculation!'// &
                   ' velocities not present / not labelled in dump file'
    endif
    !
    !--calculate mean density and rms velocity values on first pass
    !
    rmsvali = 0.
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
       if (rhoi > 0.) then
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
       if (ivx > 0 .and. ivx <= ncolumns) then
          v2i    = dot_product(dat(i,ivx:ivx+ndimV-1),dat(i,ivx:ivx+ndimV-1))
          rmsvali = rmsvali + voli*v2i
          rmsvmw = rmsvmw + v2i
       endif
    enddo
    !
    !--use the computed volume for velocity, otherwise won't be normalised correctly
    !
    rmsvali = sqrt(rmsvali/totvol)
    rmsvmw = sqrt(rmsvmw/dble(ntot))

    !
    !--option to override volume from sum with environment variable
    !
    voli = renvironment('SPLASH_CALC_VOLUME',errval=-1.0)
    if (voli > 0.) then
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
       if (rhoi > 0.) then
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
    print "(1x,'rms velocity     (vol. weighted) = ',es11.4)",rmsvali
    print "(1x,'rms velocity     (mass weighted) = ',es11.4)",rmsvmw
    if (rmsvali > 0.) then
       bval = sqrt(svarvw)/rmsvali
    else
       bval = 0.
    endif
    if (rmsvmw > 0.) then
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
    write(iunit,fmtstring) timei,rhomeanvw,rhomeanmw,rhovarvw,rhovarmw,sqrt(rhovarvw),sqrt(rhovarmw),&
                           rmsvali,rmsvmw,bval,bvalmw,smeanvw,smeanmw,svarvw,svarmw,sqrt(svarvw),sqrt(svarmw)
 case('kh')

    if (irho <= 0 .or. irho > ncolumns) then
       print "(a)",' ERROR in kh calculation!'// &
                   ' density not present / not labelled in dump file, skipping...'
       return
    endif
    if (ivx <= 0 .or. ivx > ncolumns) then
       print "(a)",' WARNING in kh calculation!'// &
                   ' velocities not present / not labelled in dump file'
    endif
    !
    !--calculate volume-weighted RMS for each column
    !
    ekinymax = 0.
    do i=1,ntot
       ekiny    = 0.5*dat(i,irho)*dat(i,ivx+1)**2
       ekinymax = max(ekiny,ekinymax)
    enddo
    print "(1x,a,es9.2)",'ekiny(max) = ',ekinymax
    write(iunit,"(2(es18.10,1x))") timei,ekinymax

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
    if (ntot > ntot1) then
       print*,' WARNING: nrows = ',ntot,' > nrows from previous dumpfile =',ntot1
       print*,'          ignoring all rows/particles greater than ',ntot1
    elseif (ntot < ntot1) then
       print*,' WARNING: nrows = ',ntot,' < nrows from previous dumpfile =',ntot1
       print*,'          assuming zeros for rows/particles greater than ',ntot
    endif
    ncol1 = size(datmean(1,:))
    if (ncolumns > ncol1) then
       print*,' WARNING: ncolumns = ',ncolumns,' > ncolumns from previous dumpfile =',ncol1
       print*,'          ignoring all rows/particles greater than ',ncol1
    elseif (ncolumns < ncol1) then
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

 case('ratio','minus','add')
    if (.not.allocated(datmean)) then
       allocate(datmean(size(dat(:,1)),size(dat(1,:))),stat=ierr)
       if (ierr /= 0) stop 'error allocating temporary memory in calc'
       datmean = 0.
    endif
    if (.not.allocated(datvar)) then
       allocate(datvar(size(dat(:,1)),size(dat(1,:))),stat=ierr)
       if (ierr /= 0) stop 'error allocating memory in calc'
       datvar = 0.
    endif
    ntot1 = size(datmean(:,1))
    if (ntot > ntot1) then
       print*,' WARNING: nrows = ',ntot,' > nrows from previous dumpfile =',ntot1
       print*,'          ignoring all rows/particles greater than ',ntot1
    elseif (ntot < ntot1) then
       print*,' WARNING: nrows = ',ntot,' < nrows from previous dumpfile =',ntot1
       print*,'          assuming zeros for rows/particles greater than ',ntot
    endif
    ncol1 = size(datmean(1,:))
    if (size(dat(1,:)) > ncol1) then
       print*,' WARNING: ncolumns = ',ncolumns,' > ncolumns from previous dumpfile =',ncol1
       print*,'          ignoring all rows/particles greater than ',ncol1
    elseif (ncolumns < ncol1) then
       print*,' WARNING: ncolumns = ',ntot,' < ncolumns from previous dumpfile =',ncol1
       print*,'          assuming zeros for columns greater than ',ncolumns
    endif
    ntot1 = min(ntot1,ntot)
    ncol1 = min(ncol1,size(dat(1,:)))
    if (ntot1 <= 0 .or. ncol1 <= 0) then
       print "(a,i2,a,i2,a)",' ERROR: nrows = ',ntot1,' ncolumns = ',ncol1,' aborting...'
       return
    endif

    if (nfilesread <= 1) then
       !--store first dump
       datmean(1:ntot1,1:ncol1) = dat(1:ntot1,1:ncol1)
    else
       select case(trim(analysistype))
       case('add')
          datvar(1:ntot1,1:ncol1) = dat(1:ntot1,1:ncol1) + datmean(1:ntot1,1:ncol1)    ! sum of current data and first step
       case('minus')
          datvar(1:ntot1,1:ncol1) = datmean(1:ntot1,1:ncol1) - dat(1:ntot1,1:ncol1)   ! first step minus current data
       case default ! ratio
          where (abs(datmean(1:ntot1,1:ncol1)) > epsilon(0.))
             datvar(1:ntot1,1:ncol1) = dat(1:ntot1,1:ncol1)/datmean(1:ntot1,1:ncol1)    ! ratio of current data to first step
          elsewhere
             datvar(1:ntot1,1:ncol1) = dat(1:ntot1,1:ncol1)/(datmean(1:ntot1,1:ncol1) + epsilon(0.))
          end where
       end select
       valmin = datvar(1,1)
       valmax = datvar(1,1)
       valmean = 0.
       do j=1,ncol1
          do i=1,ntot1
             valmin = min(datvar(i,j),valmin)
             valmax = max(datvar(i,j),valmin)
             valmean = valmean + datvar(i,j)
          enddo
       enddo
       valmean = valmean/real(ntot1*ncol1)
       print "(/,a,es10.3)",' max '//trim(analysistype)//' = ',valmax
       print "(a,es10.3)",' min '//trim(analysistype)//' = ',valmin
       print "(a,es10.3,/)",' mean '//trim(analysistype)//' = ',valmean
       print "(a)",'----> WRITING '//trim(analysistype)//'.out ...'
       if (allocated(datmean) .and. allocated(datvar)) then
          write(iunit,"('# ',i4,1x,i4)") ntot1,ncol1
          write(fmtstring,"(a,i6,a)",iostat=ierr) '(',ncol1,'(es14.6,1x,))'
          do i=1,ntot1
             write(iunit,fmtstring) datvar(i,:)
          enddo
       endif
    endif
    return
 case('lightcurve')
    call get_lightcurve(ncolumns,dat,npartoftype,massoftype,iamtype,ndim,ntypes,&
         lum,rphoto,tphoto,l_bb,r_bb,t_bb,badfrac,basename(rootname(ifileopen)),ierr)
    print "(4(/,1x,a20,' = ',es9.2))",'Luminosity',lum,'photospheric radius ',rphoto,'photospheric temperature',tphoto
    !
    !--write line to output file
    !
    if (ierr == 0) write(iunit,"(8(es18.10,1x))") timei,lum,rphoto,tphoto,l_bb,r_bb,t_bb,badfrac

 case('extinction')
    !
    !--get list of sink particle positions
    !
    isinktype = get_sink_type(ntypes)
    call get_positions_of_type(dat,npartoftype,iamtype,isinktype,ix,npts,xpts,ypts,zpts,ierr)
    if (ierr /= 0 .or. npts <= 0) then
       print*,' ERROR obtaining sink particle positions, aborting...'
       return
    endif

    call get_extinction(ncolumns,dat,npartoftype,massoftype,iamtype,ndim,ntypes,&
                        npts,xpts,ypts,zpts,coltemps,irho)

    labelc = get_unitlabel_coldens(iRescale,labelzintegration,unitslabel(irho))
    print "(100(/,1x,a20,i0,' = ',es9.2,a))",('Sigma to sink ',i,coltemps(i),trim(labelc),i=1,npts)
    !
    !--write line to output file
    !
    write(iunit,"(100(es18.10,5x))") timei,coltemps(1:npts)

 case('tdiffuse')

    icol = ienvironment('SPLASH_ANALYSIS_COL')
    if (icol <= 0 .or. icol > ncolumns) then
       icol = irho
    endif

    call get_extinction_los(ncolumns,dat,npartoftype,massoftype,iamtype,ndim,ntypes,&
                            ndirs,anglex_vals,angley_vals,anglez_vals,coltemps,icol)

    labelc = get_unitlabel_coldens(iRescale,labelzintegration,unitslabel(icol))
    print "(100(/,1x,a20,' = ',es12.4,a))",('Sigma in '//trim(dir_label(i)),coltemps(i),trim(labelc),i=1,ndirs)
    !
    !--write line to output file
    !
    write(iunit,"(100(es18.10,5x))") timei,coltemps(1:ndirs)

 case default
    print "(a)",' ERROR: unknown analysis type in write_analysis routine'
    return
 end select

 if (nfileout==1) print "(/,1x,'>>> ',a,' <<<')",'written to '//trim(fileout(1))

 return

contains
!
!--small internal utility to work out the particle type
!  (which depends whether or not mixed types are stored)
!
integer function igettype(i)
 integer :: np
 integer, intent(in) :: i

 if (size(iamtype) > 1) then
    igettype = int(iamtype(i))
 else
    np = 0
    igettype = 0
    do while (i > np .and. igettype <= ntypes)
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
 integer, intent(in) :: i,iparttype

 if (ipmass > 0 .and. ipmass <= ncolumns) then
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
 character(len=*), intent(in) :: analysistype
 integer :: i

 select case(trim(analysistype))
 case('timeaverage','timeav')
    print "(a)",'----> WRITING time_average.out ...'
    if (allocated(datmean) .and. allocated(datvar) .and. nfilesread > 0) then
       !--get standard deviation from variance (also normalise with 1/n)
       datvar(:,:) = sqrt(datvar(:,:))/sqrt(real(nfilesread))
       do i=1,size(datmean(:,1))
          write(iunit,"(1x,99(es15.7,2x))") datmean(i,:),datvar(i,:)
       enddo
    endif
 end select

 if (allocated(datmean)) deallocate(datmean)
 if (allocated(datvar)) deallocate(datvar)

 do i=1,nfileout
    close(iunit+i-1)
 enddo

 return
end subroutine close_analysis

end module analysis
