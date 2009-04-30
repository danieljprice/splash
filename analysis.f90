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
 integer, private :: nlevels
 real, dimension(maxlevels), private :: rholevels
 character(len=64), private :: fileout

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
 end select
 
 if (present(noprint)) then
    doprint = .not.noprint
 else
    doprint = .true.
 endif
 
 if (.not.isanalysis .and. doprint) then
    print "(a)",' analysis mode ("splash calc X dumpfiles") on a sequence of dump files: '
    print "(a)",' splash calc energies     : calculate KE,PE,total energy vs time'
    print "(a)",'                            output to file called ''energy.out'''
    print "(a)",'        calc massaboverho : mass above a series of density thresholds vs time'
    print "(a)",'                            output to file called ''massaboverho.out'''
    print "(a)",'        calc max          : maximum of each column vs. time'
    print "(a)",'                            output to file called ''maxvals.out'''
    print "(a)",'        calc min          : minimum of each column vs. time'
    print "(a)",'                            output to file called ''minvals.out'''
    print "(a)",'        calc mean         : mean of each column vs. time'
    print "(a)",'                            output to file called ''meanvals.out'''
 endif
 
 return
end function isanalysis

!----------------------------------------------------------------
!  open output file/ initialise quantities needed for analysis 
!  over all dump files
!----------------------------------------------------------------
subroutine open_analysis(dumpfile,analysistype,required,ncolumns,ndimV)
 use labels, only:ivx,iBfirst,iutherm,irho,ipmass,label
 use asciiutils, only:read_asciifile
 use filenames, only:rootname,nfiles
 implicit none
 character(len=*), intent(in) :: dumpfile,analysistype
 logical, dimension(0:ncolumns), intent(out) :: required
 character(len=1170) :: headerline   ! len=64 x 18 characters
 character(len=64) :: levelsfile
 character(len=40) :: fmtstring
 integer, intent(in) :: ncolumns,ndimV
 logical :: iexist
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
          1,'time',(i+1,'M(d>',rholevels(i),')',i=1,nlevels)

 case('max','maxvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'maxvals.out'
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) ncolumns+1
    write(headerline,fmtstring) 1,'time',(i+1,label(i)(1:12),i=1,ncolumns)

 case('min','minvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'minvals.out'
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) ncolumns+1
    write(headerline,fmtstring) 1,'time',(i+1,label(i)(1:12),i=1,ncolumns)

 case('mean','meanvals')
    !
    !--read all columns from dump file
    !
    required(:) = .true.
    !
    !--set filename and header line
    !
    fileout = 'meanvals.out'
    write(fmtstring,"('(''#'',1x,',i3,'(''['',i2.2,1x,a12,'']'',2x))')",iostat=ierr) ncolumns+1
    write(headerline,fmtstring) 1,'time',(i+1,label(i)(1:12),i=1,ncolumns)

 end select
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
    write(iunit,"(a)") '# SPLASH: A visualisation tool for SPH data (c)2004-2008 Daniel Price'
    write(iunit,"(a)") '# '//trim(fileout)//' produced using "splash calc '//trim(analysistype)// &
                       '" on dump files '//trim(rootname(1))//'->'//trim(rootname(nfiles))
    write(iunit,"(a)") '# use asplash -e '//trim(fileout)//' to plot the contents of this file '
    write(iunit,"(a)") '#'
    write(iunit,"(a)") trim(headerline)
 endif
 
 return
end subroutine open_analysis

!----------------------------------------------------------------
!  this is the routine which actually calculates the quantities
!  required from each dump file and spits out a line to the
!  analysis file. Called once for each dump file.
!----------------------------------------------------------------
subroutine write_analysis(time,dat,ntot,ntypes,npartoftype,massoftype,iamtype,ncolumns,ndimV,analysistype)
 use labels, only:ivx,iBfirst,iutherm,irho,ipmass,labeltype,label
 use params, only:int1,doub_prec,maxplot
 use asciiutils, only:ucase
 implicit none
 integer, intent(in) :: ntot,ntypes,ncolumns,ndimV
 integer, intent(in), dimension(:) :: npartoftype
 real, intent(in), dimension(:) :: massoftype
 integer(kind=int1), intent(in), dimension(:) :: iamtype
 real, intent(in) :: time
 real, intent(in), dimension(:,:) :: dat
 character(len=*), intent(in) :: analysistype
 real(kind=doub_prec), dimension(maxlevels) :: massaboverho
 integer :: itype,i,ierr
 real(kind=doub_prec) :: ekin,emag,etherm,epot,etot,totmom,pmassi
 real(kind=doub_prec), dimension(3) :: xmom
 character(len=20) :: fmtstring
!
! array with one value for each column
!
 real, dimension(maxplot) :: coltemp

 print "(/,5('-'),a,/)",'> CALCULATING '//trim(ucase(analysistype))
 print "(1x,'time = ',es8.2)",time
 
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
subroutine close_analysis
 implicit none

 close(unit=iunit)

 return
end subroutine close_analysis

end module analysis
