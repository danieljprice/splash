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

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR GENERAL ASCII DATA FORMATS
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! ASPLASH_COLUMNSFILE gives the location of the default 'columns' file
! (overridden by the presence of a `columns' file in the working directory)
!
! ASPLASH_NCOLUMNS can be used to override the automatic ncolumns choice
!
! e.g. setenv ASPLASH_NCOLUMNS=10
!
! ASPLASH_NHEADERLINES can be used to override the automatic number of header line determination
!
! e.g. setenv ASPLASH_NHEADERLINES=1
!
! ASPLASH_TIMEVAL can be used to set the time (fixed for all files)
! ASPLASH_GAMMAVAL can be used to set gamma (fixed for all files)
! ASPLASH_HEADERLINE_TIME can be used to set the header line where the time is listed
! ASPLASH_HEADERLINE_GAMMA can be used to set the header line where gamma is listed
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
 integer :: icoltype

end module asciiread

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep,iamtype
  use params
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,iverbose,ntypes
  use mem_allocation, only:alloc
  use asciiutils,     only:get_ncolumns
  use system_utils,   only:ienvironment,renvironment
  use asciiread,      only:icoltype
  use labels,         only:labeltype,print_types
  implicit none
  integer, intent(in)          :: indexstart,ipos
  integer, intent(out)         :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,iunit,ncolstep,ncolenv,nerr,iheader_time,iheader_gamma
  integer :: nprint,npart_max,nstep_max,icol,nheaderlines,nheaderenv,itype
  integer :: noftype(maxparttypes),iverbose_was
  logical :: iexist,timeset,gammaset
  real    :: dummyreal
  character(len=len(rootname)+4) :: dumpfile
  character(len=40)  :: line
  integer, parameter :: notset = -66

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  iunit = 15  ! logical unit number for input

  dumpfile = trim(rootname)

  if (iverbose.gt.0) print "(1x,a)",'reading ascii format'
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

  !
  !--open the file and read the number of particles
  !
  open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='formatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  else
     call get_ncolumns(iunit,ncolstep,nheaderlines)
     !--override header lines setting
     nheaderenv = ienvironment('ASPLASH_NHEADERLINES',-1)
     if (nheaderenv.ge.0) then
        if (iverbose.gt.0) print*,' setting nheader lines = ',nheaderenv,' from ASPLASH_NHEADERLINES environment variable'
        nheaderlines = nheaderenv
     endif
     !--override columns setting with environment variable
     ncolenv = ienvironment('ASPLASH_NCOLUMNS',-1)
     if (ncolenv.gt.0) then
        if (iverbose.gt.0) print "(a,i3,a)",' setting ncolumns = ',ncolenv,' from ASPLASH_NCOLUMNS environment variable'
        ncolstep = ncolenv
     endif
     if (ncolstep.le.0) then
        print "(a)",'*** ERROR: zero/undetermined number of columns in file ***'
        return
     endif
     iverbose_was = iverbose
     iverbose = 0
     ncolumns = ncolstep
     call set_labels()  ! to see if types are defined
     iverbose = iverbose_was

     !
     !--allocate memory initially
     !
     nprint = 101
     nstep_max = max(nstep_max,indexstart,1)
     if (.not.allocated(dat) .or. (nprint.gt.npart_max) .or. (ncolstep+ncalc).gt.maxcol) then
        npart_max = max(npart_max,INT(1.1*(nprint)))
        call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=(icoltype > 0))
     endif
  endif

  npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
  if (j.gt.maxstep) then
     call alloc(maxpart,j+1,maxcol,mixedtypes=(icoltype > 0))
  endif

!
!--can set either set the time and gamma explicitly
!  using environment variables (fixed for all files)
!  or can specify on which header line the time appears
!
  timeset = .false.
  gammaset = .false.
  dummyreal = renvironment('ASPLASH_TIMEVAL',errval=-1.)
  if (dummyreal.gt.0.) then
     time(j) = dummyreal
     timeset = .true.
  endif
  dummyreal = renvironment('ASPLASH_GAMMAVAL',errval=-1.)
  if (dummyreal.gt.0.) then
     gamma(j) = dummyreal
     gammaset = .true.
  endif
  iheader_time  = ienvironment('ASPLASH_HEADERLINE_TIME',errval=notset)
  iheader_gamma = ienvironment('ASPLASH_HEADERLINE_GAMMA',errval=notset)
!
!--read header lines, try to use it to set time
!
  if (nheaderlines.gt.0 .and. iverbose.gt.0) print*,'skipping ',nheaderlines,' header lines'

  do i=1,nheaderlines
     !--read header lines as character strings
     !  so that blank lines are counted in nheaderlines
     read(iunit,"(a)",iostat=ierr) line
     read(line,*,iostat=ierr) dummyreal

     if (i.eq.iheader_time .and. .not.timeset) then
        if (ierr.eq.0) then
           time(j) = dummyreal
           timeset = .true.
           print*,'setting time = ',dummyreal,' from header line ',i
           print*,'(determined from ASPLASH_HEADERLINE_TIME setting)'
        else
           print "(a,i2,a)",' ** ERROR reading time from header line ',i, &
                            ' (using ASPLASH_HEADERLINE_TIME)'
        endif
     elseif (i.eq.iheader_gamma .and. .not.gammaset) then
        if (ierr.eq.0) then
           gamma(j) = dummyreal
           gammaset = .true.
           print*,'setting gamma = ',dummyreal,' from header line ',i
           print*,'(determined from ASPLASH_HEADERLINE_GAMMA setting)'
        else
           print "(a,i2,a)",' ** ERROR reading gamma from header line ',i, &
                            ' (using ASPLASH_HEADERLINE_GAMMA)'
        endif
     elseif (timeset .and. .not.gammaset .and. ierr.eq.0 .and. iheader_gamma.eq.notset &
        .and. dummyreal.gt.0.999999 .and. dummyreal.lt.2.000001) then
        print*,'setting gamma = ',dummyreal,' from header line ',i
        gamma(j) = dummyreal
        gammaset = .true.
     elseif (ierr.eq.0 .and. .not. timeset .and. iheader_time.eq.notset) then
        time(j) = dummyreal
        timeset = .true.
        print*,'setting time = ',dummyreal,' from header line ',i
     endif
  enddo
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
     if (i.gt.npart_max) then ! reallocate memory if necessary
        npart_max = 10*npart_max
        call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=(icoltype > 0))        
     endif
     read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,ncolstep)
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
        if (nerr .le. 10) print "(a,i8,a)",' ERROR reading data from line ',i+nheaderlines,', skipping'
        i = i - 1 ! ignore lines with errors
     endif
  enddo overparts

  nprint = i - 1
  nstepsread = nstepsread + 1

  if (nerr > 10) then
     print "(a,i8,a)",' *** WARNING: errors whilst reading file on ',nerr,' lines: skipped these ***'
  endif
  if (ierr < 0) then
     print "(2(a,i10))",' read npts = ',nprint,' ncolumns = ',ncolstep
  endif

  npartoftype(:,j) = 0
  if (icoltype > 0 .and. icoltype <= ncolstep) then
     npartoftype(1:ntypes,j) = noftype(1:ntypes)
     call print_types(npartoftype(:,j),labeltype)
  else
     npartoftype(1,j) = nprint
  endif

  close(iunit)

return
end subroutine read_data

!!-------------------------------------------------------------------
!! set labels for each column of data
!!
!! read these from a file called 'columns' in the current directory
!! then take sensible guesses as to which quantities are which
!! from the column labels
!!
!!-------------------------------------------------------------------

subroutine set_labels
  use asciiutils,      only:lcase
  use labels,          only:label,labeltype,ix,irho,ipmass,ih,iutherm, &
                            ipr,ivx,iBfirst,iamvec,labelvec,lenlabel
  !use params,          only:maxparttypes
  use settings_data,   only:ncolumns,ndim,ndimV,UseTypeInRenderings,iverbose
  use geometry,        only:labelcoord
  use system_commands, only:get_environment
  use filenames,       only:fileprefix
  use asciiread,       only:icoltype
  implicit none
  integer                 :: i,ierr,ndimVtemp
  character(len=120)      :: columnfile
  character(len=lenlabel) :: labeli
  logical                 :: iexist
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
     call get_environment('ASPLASH_COLUMNSFILE',columnfile)
     if (len_trim(columnfile).gt.0) then
        inquire(file=trim(columnfile),exist=iexist)
        if (iexist) then
           if (iverbose.gt.0) print "(a)",' using ASPLASH_COLUMNSFILE='//trim(columnfile)
        else
           print "(a)",' ERROR: ASPLASH_COLUMNSFILE='//trim(columnfile)//' DOES NOT EXIST'
           columnfile = 'columns'
        endif
     else
        columnfile = 'columns'
     endif
  endif

  open(unit=51,file=trim(columnfile),status='old',iostat=ierr)
  if (ierr /=0) then
     if (iverbose > 0) then
        print "(3(/,a))",' WARNING: columns file not found: using default labels',&
                         ' To change the labels, create a file called ''columns'' ',&
                         '  in the current directory with one label per line'
     endif
  else
     overcols: do i=1,ncolumns
        read(51,"(a)",iostat=ierr) label(i)
!
!--compare all strings in lower case, trimmed and with no preceding spaces
!
        labeli = trim(adjustl(lcase(label(i))))
!
!--guess positions of various quantities from the column labels
!
        if (ndim.le.0 .and. (labeli(1:1).eq.'x' .or. labeli(1:1).eq.'r')) then
           ndim = 1
           ix(1) = i
        endif
        if (ndim.eq.1 .and. i.eq.ix(1)+1 .and. (labeli(1:1).eq.'y' .or. labeli(1:1).eq.'z')) then
           ndim = 2
           ix(2) = i
        endif
        if (ndim.eq.2 .and. i.eq.ix(2)+1 .and. labeli(1:1).eq.'z') then
           ndim = 3
           ix(3) = i
        endif
        if (labeli(1:3).eq.'den' .or. index(labeli,'rho').ne.0 .or. labeli(1:3).eq.'\gr') then
           irho = i
        elseif (labeli(1:5).eq.'pmass' .or. labeli(1:13).eq.'particle mass' &
                .or. index(labeli,'mass').ne.0) then
           ipmass = i
        elseif (ipmass.eq.0 .and. trim(labeli).eq.'m') then
           ipmass = i
        !--use first column labelled h as smoothing length
        elseif (ih.eq.0 .and. (labeli(1:1).eq.'h' &
                .or. labeli(1:6).eq.'smooth')) then
           ih = i
        elseif (trim(labeli).eq.'u'.or.labeli(1:6).eq.'utherm' &
            .or.trim(labeli).eq.'internal energy') then
           iutherm = i
        elseif (labeli(1:2).eq.'pr' .or. trim(labeli).eq.'p') then
           ipr = i
        elseif (ivx.eq.0 .and. labeli(1:1).eq.'v') then
           ivx = i
           ndimV = 1
        elseif (icoltype==0 .and. index(labeli,'type').ne.0) then
           icoltype = i
        endif
        !--set ndimV as number of columns with v as label
        if (ivx.gt.0 .and. i.gt.ivx .and. i.le.ivx+2) then
           if (labeli(1:1).eq.'v') ndimV = i - ivx + 1
        endif
        if (iBfirst.eq.0 .and. (labeli(1:2).eq.'bx')) then
           iBfirst = i
        endif
        !--set ndimV as number of columns with v as label
        if (iBfirst.gt.0 .and. i.gt.iBfirst .and. i.le.iBfirst+2) then
           if (labeli(1:1).eq.'b') then
              ndimVtemp = i - iBfirst + 1
              if (ndimV.gt.0 .and. ndimVtemp.gt.ndimV) then
                 if (iverbose > 0) print "(a)",' WARNING: possible confusion with vector dimensions'
                 ndimV = ndimVtemp
              endif
           endif
        endif
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

  if (ndim.lt.1) ndimV = 0
  if (iverbose > 0) then
     if (ndim.gt.0) print "(a,i1)",' Assuming number of dimensions = ',ndim
     if (ndim.gt.0) print "(a,i2,a,i2)",' Assuming positions in columns ',ix(1),' to ',ix(ndim)

     if (ndimV.gt.0) print "(a,i1)",' Assuming vectors have dimension = ',ndimV
     if (irho.gt.0) print "(a,i2)",' Assuming density in column ',irho
     if (ipmass.gt.0) print "(a,i2)",' Assuming particle mass in column ',ipmass
     if (ih.gt.0) print "(a,i2)",' Assuming smoothing length in column ',ih
     if (iutherm.gt.0) print "(a,i2)",' Assuming thermal energy in column ',iutherm
     if (ipr.gt.0) print "(a,i2)",' Assuming pressure in column ',ipr
     if (ivx.gt.0) then
        if (ndimV.gt.1) then
           print "(a,i2,a,i2)",' Assuming velocity in columns ',ivx,' to ',ivx+ndimV-1
        else
           print "(a,i2)",' Assuming velocity in column ',ivx
        endif
     endif
     if (icoltype.gt.0) print "(a,i2)",' Assuming particle type in column ',icoltype

     if (ndim.eq.0 .or. irho.eq.0 .or. ipmass.eq.0 .or. ih.eq.0) then
        print "(4(/,a))",' NOTE: Rendering capabilities cannot be enabled', &
                    '  until positions of density, smoothing length and particle', &
                    '  mass are known (for the ascii read the simplest way is to ', &
                    '  label the relevant columns appropriately in the columns file)'
     endif
  endif

  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
       label(ivx+i-1) = 'v\d'//labelcoord(i,1)
     enddo
  endif
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = ivx
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
       label(iBfirst+i-1) = 'B\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  !ntypes = 1 !!maxparttypes
  labeltype(1) = 'gas'
  UseTypeInRenderings(1) = .true.


!-----------------------------------------------------------

  return
end subroutine set_labels
