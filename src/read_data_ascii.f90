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
!  Copyright (C) 2005-2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR GENERAL ASCII DATA FORMATS
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING COMMAND LINE FLAGS:
!
! --columnsfile gives the location of the default 'columns' file
! (overridden by the presence of a `columns' file in the working directory)
!
! --ncolumns can be used to override the automatic ncolumns choice
!
! e.g. --ncolumns=10
!
! --nheaderlines can be used to override the automatic number of header line determination
!
! e.g. --nheaderlines=1
!
! --time=1.0 can be used to set the time (fixed for all files)
! --gamma=1.6 can be used to set gamma (fixed for all files)
! --timeheader=3 can be used to set the header line where the time is listed
! --gammaheader=5 can be used to set the header line where gamma is listed
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
! ntot(maxstep)       : total number of particles in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------
module asciiread
 use labels, only:lenlabel,label
 integer :: icoltype
 character(len=lenlabel), dimension(size(label)) :: label_orig

end module asciiread


module readdata_ascii
 implicit none

 public :: read_data_ascii, set_labels_ascii

 private
contains

subroutine read_data_ascii(rootname,indexstart,ipos,nstepsread)
 use particle_data,  only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep,iamtype
 use params
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,iverbose,ntypes
 use mem_allocation, only:alloc
 use asciiutils,     only:get_ncolumns,get_column_labels,isdigit,readline_csv
 use system_utils,   only:ienvironment,renvironment
 use asciiread,      only:icoltype,label_orig
 use labels,         only:lenlabel,labeltype,print_types,label
 use, intrinsic :: ieee_arithmetic
 integer, intent(in)          :: indexstart,ipos
 integer, intent(out)         :: nstepsread
 character(len=*), intent(in) :: rootname
 integer :: i,j,ierr,iunit,ncolstep,ncolenv,nerr,iheader_time,iheader_gamma
 integer :: nprint,npart_max,nstep_max,nheaderlines,nheaderenv,itype,nlabels
 integer :: noftype(maxparttypes),iverbose_was,imethod
 logical :: iexist,timeset,gammaset,got_labels,csv
 real    :: dummyreal
 real, allocatable :: dattemp(:)
 character(len=len(rootname)+4) :: dumpfile
 character(len=4096)  :: line
 character(len=lenlabel), dimension(size(label)) :: tmplabel
 character(len=10)  :: str,strc
 integer, parameter :: notset = -66

 nstepsread = 0
 nstep_max = 0
 npart_max = maxpart
 iunit = 15  ! logical unit number for input

 dumpfile = trim(rootname)

 if (iverbose > 1) print "(1x,a)",'reading ascii format'
 print "(26('>'),1x,a,1x,26('<'))",trim(dumpfile)
 !
 !--check if first data file exists
 !
 inquire(file=dumpfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
    return
 endif
 !
 !--fix number of spatial dimensions (0 means no particle coords)
 !
 ndim = 0
 ndimV = 0

 j = indexstart
 nstepsread = 0
 icoltype = 0       ! no particle type defined by default
 got_labels = .false.
 label_orig = ''
 csv = index(dumpfile,'.csv') > 0  ! if filename contains .csv
 !
 !--open the file and read the number of particles
 !
 open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='formatted')
 if (ierr /= 0) then
    print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
    return
 else
    call get_ncolumns(iunit,ncolstep,nheaderlines,csv=csv)
    !--override header lines setting
    nheaderenv = ienvironment('ASPLASH_NHEADERLINES',-1)
    if (nheaderenv >= 0) then
       if (iverbose > 0) print*,' setting nheader lines = ',nheaderenv,' from --nheaderlines flag'
       nheaderlines = nheaderenv
    endif
    !--override columns setting with environment variable
    ncolenv = ienvironment('ASPLASH_NCOLUMNS',-1)
    if (ncolenv > 0) then
       if (iverbose > 0) print "(a,i3,a)",' setting ncolumns = ',ncolenv,' from --ncolumns flag'
       ncolstep = ncolenv
    endif
    if (ncolstep <= 1) then
       !print "(a)",'*** ERROR: could not determine number of columns in file ***'
       print "(/,a)",' Are you trying to read a non-ascii file? If so, use:'
       print "(' ',/,a,/,' ')",'  splash -f <format> '//trim(dumpfile)
       print "(a)",' Use splash --formats for list of supported data formats '
       return
    endif

    if (ncolstep > 1) then
       !--search through header for column labels
       do i=1,nheaderlines
          read(iunit,"(a)",iostat=ierr) line
          !--try to match column labels from this header line, if not already matched (or dubious match)
          call get_column_labels(trim(line),nlabels,tmplabel,&
               method=imethod,ndesired=ncolstep,csv=csv)
          !--if we get nlabels > ncolumns, use them, but keep trying for a better match
          if ((got_labels .and. nlabels == ncolstep) .or. &
              (.not.got_labels .and. nlabels >= ncolstep  & ! only allow single-spaced labels if == ncols
               .and. (.not.(imethod>=4).or.nlabels==ncolstep))) then
             label_orig(1:ncolstep) = tmplabel(1:ncolstep)
             got_labels = .true.
             !print*,'DEBUG: line ',i,' nlabels = ',nlabels,' LABELS= '//tmplabel(1:ncolstep)
          endif
       enddo
    endif
    rewind(iunit)

    iverbose_was = iverbose
    iverbose = 0
    ncolumns = ncolstep
    call set_labels_ascii()  ! to see if types are defined
    iverbose = iverbose_was

    !
    !--allocate memory initially
    !
    nprint = 101
    nstep_max = max(nstep_max,indexstart,1)
    if (.not.allocated(dat) .or. (nprint > npart_max) .or. (ncolstep+ncalc) > maxcol) then
       npart_max = max(npart_max,INT(1.1*(nprint)))
       call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=(icoltype > 0))
    endif
 endif

 npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
 if (j > maxstep) then
    call alloc(maxpart,j+1,maxcol,mixedtypes=(icoltype > 0))
 endif
!
!--can set either set the time and gamma explicitly
!  using environment variables (fixed for all files)
!  or can specify on which header line the time appears
!
 timeset = .false.
 gammaset = .false.
 dummyreal = renvironment('ASPLASH_TIME',errval=-1.)
 if (dummyreal > 0.) then
    time(j) = dummyreal
    timeset = .true.
 endif
 dummyreal = renvironment('ASPLASH_GAMMA',errval=-1.)
 if (dummyreal > 0.) then
    gamma(j) = dummyreal
    gammaset = .true.
 endif
 iheader_time  = ienvironment('ASPLASH_TIMEHEADER',errval=notset)
 iheader_gamma = ienvironment('ASPLASH_GAMMAHEADER',errval=notset)
!
!--read header lines, try to use it to set time
!
 !if (nheaderlines > 0 .and. iverbose > 0) print*,'skipping ',nheaderlines,' header lines'
 do i=1,nheaderlines
    !--read header lines as character strings
    !  so that blank lines are counted in nheaderlines
    read(iunit,"(a)",iostat=ierr) line
    if (line(1:1)=='#') then
       read(line(2:),*,iostat=ierr) dummyreal
    else
       read(line,*,iostat=ierr) dummyreal
    endif
    if (i==iheader_time .and. .not.timeset) then
       if (ierr==0) then
          time(j) = dummyreal
          timeset = .true.
          print*,'setting time = ',dummyreal,' from header line ',i
          print*,'(determined from --timeheader setting)'
       else
          print "(a,i2,a)",' ** ERROR reading time from header line ',i, &
                            ' (using --timeheader)'
       endif
    elseif (i==iheader_gamma .and. .not.gammaset) then
       if (ierr==0) then
          gamma(j) = dummyreal
          gammaset = .true.
          print*,'setting gamma = ',dummyreal,' from header line ',i
          print*,'(determined from --gammaheader setting)'
       else
          print "(a,i2,a)",' ** ERROR reading gamma from header line ',i, &
                            ' (using --gammaheader)'
       endif
    elseif (timeset .and. .not.gammaset .and. ierr==0 .and. iheader_gamma==notset &
        .and. dummyreal > 1.0 .and. dummyreal < 2.000001) then
       print*,'setting gamma = ',dummyreal,' from header line ',i
       gamma(j) = dummyreal
       gammaset = .true.
    elseif (ierr==0 .and. .not. timeset .and. iheader_time==notset .and. index(line,'.') /= 0) then
       time(j) = dummyreal
       timeset = .true.
       print*,'setting time = ',dummyreal,' from header line ',i
    endif
 enddo
!
! allocate temporary array to read each line
!
 allocate(dattemp(ncolstep))
!
!--now read the timestep data in the dumpfile
!
 i = 0
 ierr = 0
 nerr = 0
 noftype(:) = 0
 ntypes = 1
 overparts: do while (ierr >= 0)
    i = i + 1
    if (i > npart_max) then ! reallocate memory if necessary
       npart_max = 10*npart_max
       call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=(icoltype > 0))
    endif
    dattemp(1:ncolstep) = ieee_value(1., ieee_quiet_nan)  ! NaN if not read
    if (csv) then
       read(iunit,"(a)",iostat=ierr) line
       call readline_csv(line,ncolstep,dattemp)
    else
       read(iunit,*,iostat=ierr) dattemp(1:ncolstep)
    endif
    !print*,ncolstep,nheaderlines,'line ',i,' got',dattemp(1:10)
    !read*
    dat(i,1:ncolstep,j) = dattemp(1:ncolstep)
    if (icoltype > 0 .and. icoltype <= ncolstep .and. ierr==0 .and. (size(iamtype(:,j)) > 1)) then
       !--set particle type from type column
       itype = nint(dat(i,icoltype,j))
       if (itype > 0 .and. itype < maxparttypes) then
          iamtype(i,j) = int(itype,kind=1)
       else
          iamtype(i,j) = 1
       endif
       itype = iamtype(i,j)
       noftype(itype) = noftype(itype) + 1
       ntypes = max(itype,ntypes)
    endif
    if (ierr > 0) then
       nerr = nerr + 1
       if (nerr  <=  10) print "(a,i8,a)",' ERROR reading data from line ',i+nheaderlines,', skipping'
       i = i - 1 ! ignore lines with errors
    endif
 enddo overparts
 if (allocated(dattemp)) deallocate(dattemp)

 nprint = i - 1
 nstepsread = nstepsread + 1

 if (nerr > 10) then
    print "(a,i8,a)",' *** WARNING: errors whilst reading file on ',nerr,' lines: skipped these ***'
 endif
 if (ierr < 0 .and. (icoltype <=0 .or. icoltype > ncolstep)) then
    write(str,"(i10)") nprint
    str = adjustl(str)
    write(strc,"(i10)") ncolstep
    strc = adjustl(strc)
    if (nheaderlines > 0 .and. iverbose > 0) then
       if (nheaderlines > 10) then
          print "(a,i3,a)",' npts = '//trim(str)//', ncols = '//trim(strc)//', skipped ',nheaderlines,' header lines'
       else
          print "(a,i1,a)",' npts = '//trim(str)//', ncols = '//trim(strc)//', skipped ',nheaderlines,' header lines'
       endif
    else
       print "(a)",' npts = '//trim(str)//', ncols = '//trim(strc)
    endif
 endif

 npartoftype(:,j) = 0
 if (icoltype > 0 .and. icoltype <= ncolstep) then
    npartoftype(1:ntypes,j) = noftype(1:ntypes)
    call print_types(npartoftype(:,j),labeltype)
 else
    npartoftype(1,j) = nprint
 endif

 close(iunit)

end subroutine read_data_ascii

!-------------------------------------------------------------------
! set labels for each column of data
!
! read these from a file called 'columns' in the current directory
! then take sensible guesses as to which quantities are which
! from the column labels
!
!-------------------------------------------------------------------
subroutine set_labels_ascii
 use asciiutils,      only:lcase,match_taglist,find_repeated_tags,add_escape_chars
 use labels,          only:label,labeltype,ix,irho,ipmass,ih,iutherm, &
                            ipr,ivx,iBfirst,iamvec,labelvec,lenlabel, &
                            make_vector_label
 use settings_data,   only:ncolumns,ndim,ndimV,UseTypeInRenderings,iverbose
 use geometry,        only:labelcoord
 use system_utils,    only:get_environment_or_flag
 use filenames,       only:fileprefix
 use asciiread,       only:icoltype,label_orig
 integer                 :: i,ierr,ndimVtemp
 character(len=120)      :: columnfile
 character(len=lenlabel) :: labeli
 logical                 :: iexist,got_time
!
!--read column labels from the columns file if it exists
!
!  first look for a columns file in the current directory
!  either called splash.columns or just 'columns'
!
 columnfile=trim(fileprefix)//'.columns'
 inquire(file=trim(columnfile),exist=iexist)
 if (.not.iexist) then
    columnfile='columns'
    inquire(file=trim(columnfile),exist=iexist)
 endif
!
!  if it does not exist see if the environment variable is set
!  and the corresponding file exists
!
 if (.not.iexist) then
    call get_environment_or_flag('ASPLASH_COLUMNSFILE',columnfile)
    if (len_trim(columnfile) > 0) then
       inquire(file=trim(columnfile),exist=iexist)
       if (iexist) then
          if (iverbose > 0) print "(a)",' using --columnsfile='//trim(columnfile)
       else
          print "(a)",' ERROR: --columnsfile='//trim(columnfile)//' DOES NOT EXIST'
          columnfile = 'columns'
       endif
    else
       columnfile = 'columns'
    endif
 endif

 open(unit=51,file=trim(columnfile),status='old',iostat=ierr)
 if (ierr /=0) then
!     print*,'HERE ',label(1)
    if (iverbose > 0 .and. len_trim(label_orig(1))==0) then
       print "(3(/,a))",' WARNING: column labels not found in file header:',&
                         ' To change the labels, create a file called ''columns'' ',&
                         '  in the current directory with one label per line'
    endif
 else
    overcols: do i=1,ncolumns
       read(51,"(a)",iostat=ierr) label_orig(i)
       if (ierr < 0) then
          if (iverbose > 0) print "(a,i3)",' end of file in columns file: read to column ',i-1
          exit overcols
       elseif (ierr > 0) then
          if (iverbose > 0) print "(a)",' *** error reading from columns file ***'
          exit overcols
       endif
    enddo overcols
    close(unit=51)
 endif
!
!--re-copy the labels from their originals
!  this is to avoid trying to match tags on labels
!  which have already been modified by splash (e.g. vectors)
!
 do i=1,ncolumns
    if (len_trim(label_orig(i)) > 0) label(i) = trim(label_orig(i))
 enddo
!
!--first, look for 'x','y','z' as consecutive labels
!  to determine the number of dimensions
!
 call match_taglist((/'x','y','z'/),label(1:ncolumns),ix(1),ndim)
 do i=2,ndim
    ix(i) = ix(1)+i-1
 enddo
 call match_taglist((/'vx','vy','vz'/),lcase(label(1:ncolumns)),ivx,ndimV)
 call match_taglist((/'bx','by','bz'/),lcase(label(1:ncolumns)),iBfirst,ndimVtemp)
 if (ndimV==0 .and. ivx==0) call match_taglist((/'ux','uy','uz'/),lcase(label(1:ncolumns)),ivx,ndimV)
!
!--make labels safe for plotting
!
 do i=1,ncolumns
    if (len_trim(label_orig(i)) > 0) label(i) = trim(add_escape_chars(label_orig(i)))
 enddo

 got_time = .false.
 do i=1,ncolumns
!
!--compare all strings in lower case, trimmed and with no preceding spaces
!
    labeli = trim(adjustl(lcase(label(i))))
    if (trim(labeli)=='t' .or. trim(labeli)=='time') got_time = .true.
!
!--guess positions of various quantities from the column labels
!
    if (.not.got_time) then
       if (ndim <= 0 .and. (trim(labeli)=='x' .or. trim(labeli)=='r' .or. labeli(1:3)=='rad')) then
          ndim = 1
          ix(1) = i
       endif
       if (ndim==1 .and. i==ix(1)+1 .and. (labeli(1:1)=='y' .or. labeli(1:1)=='z')) then
          ndim = 2
          ix(2) = i
       endif
       if (ndim==2 .and. i==ix(2)+1 .and. labeli(1:1)=='z') then
          ndim = 3
          ix(3) = i
       endif
    endif
    if (labeli(1:3)=='den' .or. index(labeli,'rho') /= 0 .or. labeli(1:3)=='\gr' .or. &
         (index(labeli,'density') /= 0 .and. irho==0)) then
       irho = i
    elseif (labeli(1:5)=='pmass' .or. labeli(1:13)=='particle mass' &
             .or. index(labeli,'mass') /= 0) then
       ipmass = i
    elseif (ipmass==0 .and. trim(labeli)=='m') then
       ipmass = i
       !--use first column labelled h as smoothing length
    elseif (ih==0 .and. (labeli(1:1)=='h' &
             .or. labeli(1:6)=='smooth')) then
       ih = i
    elseif (trim(labeli)=='u'.or.labeli(1:6)=='utherm' &
         .or.(index(labeli,'internal energy') /= 0 .and. iutherm==0)) then
       iutherm = i
    elseif (labeli(1:2)=='pr' .or. trim(labeli)=='p' .or. &
            (index(labeli,'pressure') /= 0 .and. ipr==0)) then
       ipr = i
    elseif (icoltype==0 .and. index(labeli,'type') /= 0) then
       icoltype = i
    elseif (ivx==0 .and. ndim==1 .and. trim(labeli)=='v') then
       ivx = i
    endif
 enddo
!
!--try to find vectors by looking for multiple entries starting with 'v'
!
 if (ndim > 0 .and. ivx==0) call find_repeated_tags('v_',ncolumns,label,ivx,ndimV)
 if (ndimVtemp > ndimV .and. iverbose > 0) &
     print "(a)",' WARNING: possible confusion with vector dimensions'
!
!--ignore label identifications if no spatial coordinates found
!
 if (ndim < 1) then
    ndimV = 0
    irho = 0
    ipmass = 0
    ih = 0
    iutherm = 0
    ipr = 0
    ivx = 0
 endif
 if (iverbose > 0) then
    if (ndim > 0) print "(a,i1,a,i2,a,i2)",' Assuming ',ndim,' dimensions, coords in cols ',ix(1),' to ',ix(ndim)
    !if (ndimV > 0) print "(a,i1)",' Assuming vectors have dimension = ',ndimV
    if (ndim > 0 .and. (irho > 0 .or. ipmass>0 .or. ih > 0)) write(*,"(a)",advance='no') ' Assuming'
    if (irho > 0) write(*,"(a,i2)",advance='no') ' density in column ',irho
    if (ipmass > 0) write(*,"(a,i2)",advance='no') ', mass in ',ipmass
    if (ih > 0) write(*,"(a,i2)",advance='no') ', h in ',ih
    if (iutherm > 0) write(*,"(/,a,i2)") ' Assuming thermal energy in ',iutherm
    if (iutherm==0 .and. ipr > 0) write(*,"(/,a)",advance='no') ' Assuming'
    if (ipr > 0) write(*,"(a,i2)",advance='no') ' pressure in column ',ipr
    if (ivx > 0) then
       if (ipr==0 .and. iutherm==0) write(*,*)
       if (ndimV > 1) then
          print "(a,i2,a,i2)",' Assuming velocity in cols ',ivx,' to ',ivx+ndimV-1
       else
          print "(a,i2)",' Assuming velocity in column ',ivx
       endif
    endif
    if (icoltype > 0) then
       if (ipr==0 .and. iutherm==0 .and. ivx==0) write(*,*)
       write(*,"(a,i2)") ' Assuming particle type in column ',icoltype
    endif
    if ((irho > 0 .or. ih > 0 .or. ipmass > 0) &
        .and. ipr==0 .and. iutherm==0 .and. ivx==0 .and. icoltype==0) write(*,*)

    if (ndim > 0 .and. (irho==0 .or. ipmass==0 .or. ih==0)) then
       print "(2(/,a))",' NOTE: Rendering disabled until density, h and mass columns known', &
                    ' (i.e. label relevant columns in file header or columns file)'
    endif
 endif

 call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('B',iBfirst,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 !
 !--set labels for each particle type
 !
 labeltype(1) = ''
 UseTypeInRenderings(1) = .true.

end subroutine set_labels_ascii
end module readdata_ascii
