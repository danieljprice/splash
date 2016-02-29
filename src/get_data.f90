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
!  wrapper for the main data read
!  ensures that same procedure occurs on initial read as from menu option
!
!  drives reading of all files listed on command line
!
!  Arguments:
!   ireadfile : if < 0, reads from all files
!               if > 0, reads only from the filename rootname(ireadfile)
!               if = 0, no data read, just call labelling and exact_params
!
module getdata
 implicit none
 public           :: get_data, get_labels
 integer, private :: ncolumnsfirst

 private

contains

subroutine get_data(ireadfile,gotfilenames,firsttime,iposinfile)
  use asciiutils,     only:ucase
  use exact,          only:read_exactparams
  use filenames,      only:rootname,nstepsinfile,nfiles,nsteps,maxfile,ifileopen
  use limits,         only:set_limits
  use settings_data,  only:ncolumns,iendatstep,ncalc,ivegotdata,    &
                      DataisBuffered,iCalcQuantities,ndim,iverbose,ntypes, &
                      iRescale,required,ipartialread,lowmemorymode,debugmode
  use settings_data,  only:iexact,buffer_steps_in_file
  use particle_data,  only:dat,time,npartoftype,maxcol
  use prompting,      only:prompt
  use labels,         only:labeltype
  use calcquantities, only:calc_quantities
  use settings_units, only:units
  use timing,         only:wall_time,print_time
  use settings_part,  only:iplotpartoftype
  use geomutils,      only:set_coordlabels
  use adjustdata,     only:adjust_data_codeunits
  implicit none
  integer, intent(in) :: ireadfile
  logical, intent(in) :: gotfilenames
  logical, intent(in), optional :: firsttime
  integer, intent(in), optional :: iposinfile
  logical :: setlimits,isfirsttime
  logical, parameter  :: dotiming = .true.
  integer :: i,istart,ierr,itype,nplot,ipos,nsteps_read
  real    :: t1,t2

  if (.not.gotfilenames) then
     if (nfiles.le.0 .or. nfiles.gt.maxfile) nfiles = 1
     call prompt('Enter number of files to read ',nfiles,1,maxfile)
     do i=1,nfiles
        call prompt('Enter filename to read',rootname(i),noblank=.true.)
     enddo
  endif
  !
  !--set everything to zero initially
  !
  ncolumns = 0
  ncalc = 0
  nsteps = 0
  istart = 1
  ivegotdata = .false.
  ifileopen = ireadfile
  DataIsBuffered = .false.
  ipartialread = .false.
  isfirsttime = .false.
  if (present(firsttime)) isfirsttime = firsttime
  if (isfirsttime) then
     iverbose = 1
  else
     iverbose = 0
  endif

  ! ipos is the offset for *which* timestep to read from files containing multiple steps
  ipos = 1
  if (present(iposinfile)) then
     if (iposinfile > 0) ipos = iposinfile
  endif
  !
  !--nstepsinfile is initialised to negative
  !  this is set progressively as files are read
  !  for non-buffered data file 1 is read and the rest are assumed to be the same
  !  then these files are corrected as they are read. By initialising nstepsinfile
  !  to negative, this means that if we get dud files (with nstepsinfile=0) we
  !  know that this is really the file contents (not just an initialised value of nstepsinfile)
  !  and can skip the file on the second encounter (see timestepping.f90)
  !
  if (isfirsttime) then
     nstepsinfile(:) = -1
     ncolumnsfirst = 0
     required = .true.
     if (lowmemorymode) required = .false.
     call endian_info()
  endif

  if (ireadfile.le.0) then
     !
     !--read all steps from the data file
     !
     nstepsinfile(1:nfiles) = 0
     required = .true.
     print "(/a)",' reading ALL dumpfiles into memory'
     !call endian_info()

     do i=1,nfiles
        call read_data(rootname(i),istart,ipos,nstepsinfile(i))

        istart = istart + nstepsinfile(i) ! number of next step in data array
        if (nstepsinfile(i).gt.0 .and. ncolumnsfirst.eq.0 .and. ncolumns.gt.0) then
           ncolumnsfirst = ncolumns
        elseif (nstepsinfile(i).gt.0 .and. ncolumns.ne.ncolumnsfirst) then
           print "(a,i2,a,i2,a)",' WARNING: file contains ',ncolumns, &
           ' columns, which differs from ',ncolumnsfirst,' read previously'
           ncolumns = max(ncolumns,ncolumnsfirst)
        endif
     enddo
     nsteps = istart - 1
     if (nsteps.gt.0) then
        ivegotdata = .true.
        DataIsBuffered = .true.
     else
        ncolumns = 0
     endif
     print "(a,i6,a,i3)",' >> Finished data read, nsteps = ',nsteps,' ncolumns = ',ncolumns

     !
     !--set labels (and units) for each column of data
     !
     !print "(/a)",' setting plot labels...'
     if (ivegotdata .and. ncolumns.gt.0) then
        call get_labels
        call adjust_data_codeunits
        !
        !--do some basic sanity checks
        !
        call check_data_read()
     endif

     if (iRescale .and. any(abs(units(0:ncolumns)-1.0).gt.tiny(units))) then
        !write(*,"(/a)") ' rescaling data...'
        do i=1,ncolumns
           if (abs(units(i)-1.0).gt.tiny(units) .and. abs(units(i)).gt.tiny(units)) then
              dat(:,i,1:nsteps) = dat(:,i,1:nsteps)*units(i)
           endif
        enddo
        time(1:nsteps) = time(1:nsteps)*units(0)
     endif
     !
     !--reset coordinate and vector labels (depending on coordinate system)
     !  Need to do this BEFORE calculating quantities
     !
     if (ivegotdata) call set_coordlabels(ncolumns)
     !
     !--calculate various additional quantities
     !
     if (nsteps.ge.1 .and. iCalcQuantities) then
        call calc_quantities(1,nsteps)
     endif
     !
     !--set plot limits
     !
     if (ierr.gt.0 .and. ivegotdata .and. nstepsinfile(1).ge.1) then
        call set_limits(1,nsteps,1,ncolumns+ncalc)
     endif

  elseif (ireadfile.gt.0) then
     !
     !--read from a single file only
     !
     nstepsinfile(ireadfile) = 0
     !if (isfirsttime) print "(/a)",' reading single dumpfile'

     if (dotiming) call wall_time(t1)
     call read_data(rootname(ireadfile),istart,ipos,nstepsinfile(ireadfile))
     !
     !--do some basic sanity checks
     !
     if (debugmode) print "(a,i3)",' DEBUG: ncolumns from data read = ',ncolumns
     if (debugmode) print "(a,i3)",' DEBUG: nsteps in file   = ',nstepsinfile(ireadfile)
     if (buffer_steps_in_file) then
        nsteps_read = nstepsinfile(ireadfile)
     else
        nsteps_read = 1
     endif

!--try different endian if failed the first time
     !if (nstepsinfile(ireadfile).eq.0) then
     !   print "(a)",' trying different endian'
     !   call read_data_otherendian(rootname(ireadfile),istart,nstepsinfile(ireadfile))
     !endif
     if (dotiming) then
        call wall_time(t2)
        if (t2-t1.gt.1.) then
           if (ipartialread) then
              call print_time(t2-t1,'time for (partial) data read = ')
              print*
           else
              call print_time(t2-t1,'time for data read = ')
              print*
           endif
        endif
!        do i=1,ncolumns+ncalc
!           print*,' required(',i,') = ',required(i)
!        enddo
     endif
     !!print*,'nsteps in file = ',nstepsinfile(ireadfile)
     if (ANY(nstepsinfile(1:ireadfile).gt.0)) ivegotdata = .true.
     if (.not.ivegotdata) ncolumns = 0
     !
     !--set ncolumns on first step only
     !
     if (ivegotdata .and. ncolumnsfirst.eq.0 .and. ncolumns.gt.0) then
        ncolumnsfirst = ncolumns
     endif
     !--override ncolumns from file and warn if different to first file
     if (ncolumnsfirst.gt.0 .and. nstepsinfile(ireadfile).gt.0) then
        if (ncolumns.ne.ncolumnsfirst) then
           print "(1x,a,i2,a,i2,a)",'WARNING: file contains ',ncolumns, &
           ' columns, which differs from ',ncolumnsfirst,' read previously'
           if (ncolumns.lt.ncolumnsfirst) then
              print "(10x,a,i2,/)",'setting data = 0 for columns > ',ncolumns
              dat(:,ncolumns+1:min(ncolumnsfirst,maxcol),1:nstepsinfile(ireadfile)) = 0.
           elseif (ncolumns.gt.ncolumnsfirst) then
              print "(10x,a,i2,a)",'extra data beyond column ',ncolumnsfirst,' will be ignored'
              print "(10x,a,/)",'(read this file first to use this data)'
           endif
           ncolumns = ncolumnsfirst
        endif
     endif

     !
     !--assume there are the same number of steps in the other files
     !  which have not been read
     !
     do i=1,nfiles
        if (nstepsinfile(i).eq.-1) then
           nstepsinfile(i) = nstepsinfile(ireadfile)
        endif
     enddo
     nsteps = sum(nstepsinfile(1:nfiles))
     !
     !--set labels (and units) for each column of data
     !  allow this to be overridden by the presence of a splash.columns file
     !
     !!print "(/a)",' setting plot labels...'
     if (ivegotdata .and. ncolumns.gt.0) then
        call get_labels
        call adjust_data_codeunits
        call check_data_read()
     endif

     if (iRescale .and. any(abs(units(0:ncolumns)-1.0).gt.tiny(units))) then
        if (debugmode) write(*,"(a)") ' rescaling data...'
        do i=1,min(ncolumns,maxcol)
           if (abs(units(i)-1.0).gt.tiny(units) .and. abs(units(i)).gt.tiny(units)) then
              dat(:,i,1:nsteps_read) = dat(:,i,1:nsteps_read)*units(i)
           endif
        enddo
        do i=1,nsteps_read
           if (time(i).gt.-0.5*huge(0.)) time(i) = time(i)*units(0)        
        enddo
     endif
     !
     !--reset coordinate and vector labels (depending on coordinate system)
     !  Need to do this BEFORE calculating quantities
     !
     if (ivegotdata) call set_coordlabels(ncolumns)

     !
     !--calculate various additional quantities
     !
     if (nsteps_read.gt.0 .and. iCalcQuantities) then
        if (ipartialread .and. .not.any(required(ncolumns+1:))) then
           !--for partial data reads do a "pretend" call to calc quantities
           !  just to get ncalc and column labels right
           call calc_quantities(1,nsteps_read,dontcalculate=.true.)
        else
           call calc_quantities(1,nsteps_read)
        endif
     endif
     !
     !--only set limits if reading the first file for the first time
     !
     setlimits = (ireadfile.eq.1 .and. ivegotdata .and. nstepsinfile(1).ge.1)
     if (.not.present(firsttime)) then
        setlimits = .false.
     elseif (.not.firsttime) then
        setlimits = .false.
     endif

     if (setlimits) then
        call set_limits(1,nsteps_read,1,ncolumns+ncalc)
        !--also set iendatstep the first time around
        iendatstep = nsteps
     endif
  endif
  !
  !--check for errors in data read / print warnings
  !
  if (ndim.ne.0 .and. ncolumns.gt.0 .and. nsteps.gt.0 .and. iverbose.eq.1) then
     if (sum(npartoftype(:,1)).gt.0 .and. npartoftype(1,1).eq.0) then
        print "(a)",' WARNING! DATA APPEARS TO CONTAIN NO '//trim(ucase(labeltype(1)))//' PARTICLES'
        itype = 0
        nplot = 0
        do while (nplot.eq.0 .and. itype < ntypes)
           itype = itype + 1
           if (npartoftype(itype,1).gt.0) then
              iplotpartoftype(itype) = .true.
              nplot = nplot + npartoftype(itype,1)
              print "(a)",' (plotting of '//trim(labeltype(itype))//' particles turned ON)'
           endif
        enddo
        print*
     endif
  endif

  !
  !--read exact solution parameters from files if present
  !
  if (iexact.ne.0) then
     if (ireadfile.lt.0) then
        call read_exactparams(iexact,rootname(1),ierr)
     else
        call read_exactparams(iexact,rootname(ireadfile),ierr)
     endif
  endif

  return
end subroutine get_data

!----------------------------------------------------------------------
!
! The following is a wrapper routine for the call to set_labels which
! overrides the label setting from the splash.columns file if present.
! Also adds the units label if the data has been rescaled.
!
!----------------------------------------------------------------------
subroutine get_labels
 use asciiutils,     only:read_asciifile
 use filenames,      only:fileprefix,unitsfile
 use labels,         only:label,unitslabel
 use settings_data,  only:ncolumns,iRescale,iverbose
 use settings_units, only:read_unitsfile
 use particle_data,  only:maxcol
 use params,         only:maxplot
 implicit none
 logical :: iexist
 integer :: nlabelsread,ierr,i

 call set_labels
 !
 !--check that label settings are sensible, fix where possible
 !
 call check_labels

 !
 !--look for a .columns file to override the default column labelling
 !
 inquire(file=trim(fileprefix)//'.columns',exist=iexist)
 nlabelsread = 0
 if (iexist) then
    call read_asciifile(trim(fileprefix)//'.columns',nlabelsread,label(1:min(ncolumns,maxcol,maxplot)))
    if (nlabelsread.lt.ncolumns) &
       print "(a,i3)",' end of file in '//trim(fileprefix)//'.columns file: labels read to column ',nlabelsread
 endif
 !
 !--read units file and change units if necessary
 !
 call read_unitsfile(trim(unitsfile),ncolumns,ierr,iverbose)
 !
 !--add units labels to labels
 !
 if (iRescale) then
    do i=1,min(ncolumns,maxcol,maxplot)
       if (index(label(i),trim(unitslabel(i))).eq.0) label(i) = trim(label(i))//trim(unitslabel(i))
    enddo
 endif

end subroutine get_labels

!----------------------------------------------------------------
!
!  utility to check that label settings are sensible
!
!----------------------------------------------------------------
subroutine check_labels
 use settings_data,   only:ndim,ndimV,ncolumns,iverbose
 use labels,          only:ix,irho,ih,ipmass
 use particle_data,   only:masstype
 use settings_render, only:icolour_particles
 implicit none
 integer :: i,ndimset

 if (ndim.ne.0 .and. ncolumns.gt.0) then
    if (ndim.lt.0 .or. ndim.gt.3) then
       print "(a)",' ERROR with ndim setting in data read, using ndim=3'
       ndim = 3
    endif
    if (ndimV.lt.0 .or. ndimV.gt.3) then
       print "(a)",' ERROR with ndimV setting in data read, using ndimV=3'
       ndimV = 3
    endif
    if (ndim.ge.2  .and. any(ix(2:ndim).eq.ix(1))) then
       print "(a)",' WARNING: error in ix setting in set_labels: fixing '
       ix(1) = max(ix(1),1)
       do i=2,ndim
          ix(i) = i
       enddo
    endif
    if (ndim.ge.1) then
       do i=1,ndim
          if (ix(i).le.0) then
             ix(i) = i
             print "(a)",' WARNING: ndim > 0 but zero ix setting in set_labels: fixing '
          endif
       enddo
    endif
    ndimset = 0
    do i=1,3
       if (ix(i).ne.0) ndimset = ndimset + 1
    enddo
    if (ndimset.ne.ndim) then
       print "(2(a,i2))",' ERROR: labels for ',ndimset,&
                         ' coordinates set but got ndim = ',ndim
    endif
    if (irho.gt.ncolumns .or. irho.lt.0) then
       print "(a)",' ERROR with irho setting in data read'
       irho = 0
    endif
    if (ih.gt.ncolumns .or. ih.lt.0) then
       print "(a)",' ERROR with ih setting in data read '
       ih = 0
    endif
    if (ipmass.gt.ncolumns .or. ipmass.lt.0) then
       print "(a)",' ERROR with ipmass setting in data read'
       ipmass = 0
    endif
    if (iverbose.ge.1) then
       if (irho.eq.0 .or. ih.eq.0) then
          print "(4(/,a))",' WARNING: Rendering capabilities cannot be enabled', &
                   '  until positions of density, smoothing length and particle', &
                   '  masses are known (specified using the integer variables ', &
                   '  irho,ih and ipmass in the read_data routine)'
          icolour_particles = .true.
       elseif (irho.gt.0 .and. ih.gt.0 .and. ipmass.eq.0 .and. all(masstype(:,:).lt.tiny(0.))) then
          print "(2(/,a))",' WARNING: Particle masses not read as array but mass not set:', &
                           '          RENDERING WILL NOT WORK! '
       endif
    endif
 endif

end subroutine check_labels

!----------------------------------------------------------------
!
!  utility to check things about the data read
!
!----------------------------------------------------------------
subroutine check_data_read
 use params,        only:maxplot,maxparttypes
 use settings_data, only:ncolumns,ndim,ndimV,ntypes,ivegotdata
 use particle_data, only:npartoftype,iamtype,dat
 use labels,        only:labeltype
 implicit none
 integer :: i,j,ntoti,nunknown,itype
 integer, dimension(maxparttypes) :: noftype

 if (ncolumns.lt.0) then
    print "(a)",' ERROR: ncolumns < 0 in data read'
    ncolumns = 0
 elseif (ncolumns.gt.maxplot) then
    print "(/,71('*'),/,'*',a,i3,a,'*',/,71('*'))",&
          ' ERROR: ncolumns > ',maxplot,' in data read: cannot list all columns in menu '
    ncolumns = maxplot
 endif

 if (ndim.gt.3) then; print "(a)",' ERROR: ndim  > 3 in data read, setting ndim = 3'; ndim = 3; endif
 if (ndim.lt.0) then; print "(a)",' ERROR: ndim  < 0 in data read, setting ndim = 0'; ndim = 0; endif
 if (ndimV.gt.3) then; print "(a)",' ERROR: ndimV > 3 in data read, setting ndimV = 3'; ndimV = 3; endif
 if (ndimV.lt.0) then; print "(a)",' ERROR: ndimV < 0 in data read, setting ndimV = 0'; ndimV = 0; endif
 if (ntypes.lt.0) then; print "(a)",' ERROR: ntypes < 0 in data read'; ntypes = 0; endif

 if (allocated(npartoftype)) then
    if (size(npartoftype(:,1)).lt.ntypes) then
       print "(a)",' ERROR: too many particle types for allocated array size in data read'
       ntypes = size(npartoftype(:,1))
    endif
    do i=1,ntypes
       do j=1,size(npartoftype(i,:))
          if (npartoftype(i,j).lt.0) then
             print "(a)",' ERROR: number of '//trim(labeltype(i))//' particles < 0 in data read'
             npartoftype(i,j) = 0
          endif
       enddo
    enddo
    ntoti = sum(npartoftype(:,1))
    if (ntoti > size(dat(:,1,1))) then
       print "(2(a,i10),a)",' ERROR: size of dat array (',size(dat(:,1,1)),&
             ') too small for number of particles (',ntoti,')'
       ivegotdata = .false.
    endif
    !if (debugmode) then
       !
       !--for mixed type storage, check that the number of particles
       !  of each type adds up to npartoftype
       !
       ntoti = sum(npartoftype(:,1))
       noftype(:) = 0
       nunknown   = 0
       if (size(iamtype(:,1)).ge.ntoti) then
          do i=1,ntoti
             itype = iamtype(i,1)
             if (itype.gt.0 .and. itype.le.ntypes) then
                noftype(itype) = noftype(itype) + 1
             else
                nunknown = nunknown + 1
             endif
          enddo
          do itype=1,ntypes
             if (npartoftype(itype,1).ne.noftype(itype)) then
                print "(a,i10,a,i10)",' ERROR in data read: got ',noftype(itype),' '//trim(labeltype(itype))// &
                                      ' particles from iamtype, but npartoftype = ',npartoftype(itype,1)
             endif
          enddo
          if (nunknown.gt.0) then
             print "(a,i10,a)",' ERROR in data read: got ',nunknown, &
                               ' particles of unknown type in iamtype array from data read'
          endif
       endif
    !endif
 endif

end subroutine check_data_read

!-------------------------------------
!
! simple utility to spit out native
! endian-ness
!
!-------------------------------------

subroutine endian_info
 implicit none
 logical :: bigendian

 bigendian = IACHAR(TRANSFER(1,"a")) == 0

 if (bigendian) then
    print "(a)",' native byte order on this machine is BIG endian'
 !--we no longer warn for little endian, as this is now most common
 !else
 !   print "(a)",' native byte order on this machine is LITTLE endian'
 endif

end subroutine endian_info

end module getdata
