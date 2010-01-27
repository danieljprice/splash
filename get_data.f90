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
 public           :: get_data, get_labels, set_coordlabels
 integer, private :: ncolumnsfirst
 integer, private :: icoordsprev = -1

 private

contains

subroutine get_data(ireadfile,gotfilenames,firsttime)
  use asciiutils,     only:ucase
  use exact,          only:read_exactparams
  use filenames,      only:rootname,nstepsinfile,nfiles,nsteps,maxfile,ifileopen
  use limits,         only:set_limits
  use settings_data,  only:ncolumns,iendatstep,ncalc,ivegotdata, &
                      DataisBuffered,iCalcQuantities,ndim,       &
                      iRescale,required,ipartialread,lowmemorymode
  use settings_part,  only:iexact
  use particle_data,  only:dat,time,npartoftype,maxcol
  use prompting,      only:prompt
  use labels,         only:labeltype
  use calcquantities, only:calc_quantities
  use settings_units, only:units
  implicit none
  integer, intent(in) :: ireadfile
  logical, intent(in) :: gotfilenames
  logical, intent(in), optional :: firsttime
  logical :: setlimits,isfirsttime
  logical, parameter  :: timing = .true.
  integer :: i,istart,ierr
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
     print "(/a)",' reading from all dumpfiles...'
     !call endian_info()
     
     do i=1,nfiles
        call read_data(rootname(i),istart,nstepsinfile(i))
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
     endif
     print "(a,i6,a,i3)",' >> Finished data read, nsteps = ',nsteps,' ncolumns = ',ncolumns

     !
     !--set labels (and units) for each column of data
     !
     print "(/a)",' setting plot labels...'
     if (ivegotdata .and. ncolumns.gt.0) then
        call get_labels
        call adjust_data_codeunits
     endif
     
     if (iRescale .and. any(abs(units(0:ncolumns)-1.0).gt.tiny(units))) then
        write(*,"(/a)") ' rescaling data...'
        do i=1,ncolumns
           if (abs(units(i)-1.0).gt.tiny(units) .and. abs(units(i)).gt.tiny(units)) then
              dat(:,i,1:nsteps) = dat(:,i,1:nsteps)*units(i)
           endif
        enddo
        time(1:nsteps) = time(1:nsteps)*units(0)
     endif     
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
     print "(/a)",' reading single dumpfile'
     !call endian_info()
     if (timing) call cpu_time(t1)
     call read_data(rootname(ireadfile),istart,nstepsinfile(ireadfile))
!--try different endian if failed the first time
     !if (nstepsinfile(ireadfile).eq.0) then
     !   print "(a)",' trying different endian'
     !   call read_data_otherendian(rootname(ireadfile),istart,nstepsinfile(ireadfile))    
     !endif
     if (timing) then
        call cpu_time(t2)
        if (t2-t1.gt.1.) then
           if (ipartialread) then
              print*,'time for (partial) data read = ',t2-t1,'s'
              print*
           else
              print*,'time for data read = ',t2-t1,'s'
              print*
           endif
        endif
!        do i=1,ncolumns+ncalc
!           print*,' required(',i,') = ',required(i)
!        enddo
     endif
     !!print*,'nsteps in file = ',nstepsinfile(ireadfile)
     if (ANY(nstepsinfile(1:ireadfile).gt.0)) ivegotdata = .true.
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
     endif
     
     if (iRescale .and. any(abs(units(0:ncolumns)-1.0).gt.tiny(units))) then
        write(*,"(/a)") ' rescaling data...'
        do i=1,min(ncolumns,maxcol)
           if (abs(units(i)-1.0).gt.tiny(units) .and. abs(units(i)).gt.tiny(units)) then
              dat(:,i,1:nstepsinfile(ireadfile)) = dat(:,i,1:nstepsinfile(ireadfile))*units(i)
           endif
        enddo
        time(1:nstepsinfile(ireadfile)) = time(1:nstepsinfile(ireadfile))*units(0)
     endif
     !
     !--calculate various additional quantities
     !
     if (nstepsinfile(ireadfile).gt.0 .and. iCalcQuantities) then
        if (ipartialread .and. .not.any(required(ncolumns+1:ncolumns+ncalc))) then
           !--for partial data reads do a "pretend" call to calc quantities 
           !  just to get ncalc and column labels right
           call calc_quantities(1,nstepsinfile(ireadfile),dontcalculate=.true.)
        else
           call calc_quantities(1,nstepsinfile(ireadfile))
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
        call set_limits(1,nstepsinfile(ireadfile),1,ncolumns+ncalc)
        !--also set iendatstep the first time around
        iendatstep = nsteps
     endif
  endif
  !
  !--check for errors in data read / print warnings
  !
  if (ndim.ne.0 .and. ncolumns.gt.0 .and. nsteps.gt.0) then
     if (sum(npartoftype(:,1)).gt.0 .and. npartoftype(1,1).eq.0) then
        print "(3(/,a),/)",' WARNING! DATA APPEARS TO CONTAIN NO '//trim(ucase(labeltype(1)))//' PARTICLES:', &
                           '  nothing will appear unless plotting of other particle ', &
                           '  types is turned on via the o)ptions menu'
     endif
  endif
  !
  !--reset coordinate and vector labels (depending on coordinate system)
  !
  call set_coordlabels(ncolumns+ncalc)

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
 use labels,         only:label
 use settings_data,  only:ncolumns,iRescale
 use settings_units, only:read_unitsfile,units,unitslabel
 use particle_data,  only:maxcol
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
    call read_asciifile(trim(fileprefix)//'.columns',nlabelsread,label(1:ncolumns))
    if (nlabelsread.lt.ncolumns) &
       print "(a,i3)",' WARNING: end of file in '//trim(fileprefix)//'.columns file: labels read to column ',nlabelsread
 endif
 !
 !--read units file and change units if necessary
 !
 call read_unitsfile(trim(unitsfile),ncolumns,ierr)
 !
 !--add units labels to labels
 !
 if (iRescale) then
    do i=1,min(ncolumns,maxcol)
       if (index(label(i),trim(unitslabel(i))).eq.0) label(i) = trim(label(i))//trim(unitslabel(i))
    enddo
 endif

end subroutine get_labels

!----------------------------------------------------------------
!
!  routine to set labels for vector quantities and spatial
!  coordinates depending on the coordinate system used.
!
!----------------------------------------------------------------
subroutine set_coordlabels(numplot)
 use geometry,       only:labelcoord
 use labels,         only:label,iamvec,labelvec,ix,labeldefault
 use settings_data,  only:icoords,icoordsnew,ndim,iRescale
 use settings_units, only:unitslabel
 implicit none
 integer, intent(in) :: numplot
 integer             :: i
!
!--sanity check on icoordsnew...
!  (should not be zero)
!
 if (icoordsnew.le.0) then
    if (icoords.gt.0) then
       icoordsnew = icoords
    else
       icoordsnew = 1
    endif
 endif
!
!--store the previous value of icoordsnew that was used
!  last time we adjusted the labels
!
 if (icoordsprev.lt.0) icoordsprev = icoordsnew

!
!--set coordinate and vector labels (depends on coordinate system)
!
 if (icoordsnew.ne.icoords) then
!
!--here we are using a coordinate system that differs from the original
!  one read from the code (must change labels appropriately)
!
    print*,' changing coordinate labels ...'
    do i=1,ndim
       label(ix(i)) = labelcoord(i,icoordsnew)
       if (iRescale .and. icoords.eq.icoordsnew) then
          label(ix(i)) = trim(label(ix(i)))//trim(unitslabel(ix(i)))
       endif
    enddo
 elseif (icoordsnew.ne.icoordsprev) then
!
!--here we are reverting back to the original coordinate system
!  so we have to re-read the original labels from the data read
!
    call get_labels
 endif
!
!--set vector labels if iamvec is set and the labels are the default
! 
 if (icoordsnew.gt.0) then
    do i=1,numplot
       if (iamvec(i).ne.0 .and. &
          (icoordsnew.ne.icoords .or. icoordsnew.ne.icoordsprev &
           .or. index(label(i),trim(labeldefault)).ne.0)) then
          label(i) = trim(labelvec(iamvec(i)))//'\d'//trim(labelcoord(i-iamvec(i)+1,icoordsnew))
          if (iRescale) then
             label(i) = trim(label(i))//'\u'//trim(unitslabel(i))
          endif
       endif
    enddo
 endif
 icoordsprev = icoordsnew
 
 return
end subroutine set_coordlabels

!----------------------------------------------------------------
!
!  utility to check that label settings are sensible
!
!----------------------------------------------------------------
subroutine check_labels
 use settings_data, only:ndim,ndimV,ncolumns
 use labels,        only:ix,irho,ih,ipmass
 use particle_data, only:masstype
 implicit none
 integer :: i

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
    if (irho.eq.0 .or. ih.eq.0) then
       print "(4(/,a))",' WARNING: Rendering capabilities cannot be enabled', &
                '  until positions of density, smoothing length and particle', &
                '  masses are known (specified using the integer variables ', &
                '  irho,ih and ipmass in the read_data routine)'
    elseif (irho.gt.0 .and. ih.gt.0 .and. ipmass.eq.0 .and. all(masstype(:,:).lt.tiny(0.))) then
       print "(2(/,a))",' WARNING: Particle masses not read as array but mass not set:', &
                        '          RENDERING WILL NOT WORK! '


    endif
 endif

end subroutine check_labels

!----------------------------------------------------
!
!  amend data after the data read based on
!  various environment variable settings
!
!  must be called AFTER the data has been read
!  but BEFORE rescaling to physical units is applied
!
!----------------------------------------------------
subroutine adjust_data_codeunits
 use system_utils,  only:renvironment
 use labels,        only:ih
 use settings_data, only:ncolumns
 use particle_data, only:dat
 implicit none
 real :: hmin
 
 if (ih.gt.0 .and. ih.le.ncolumns) then
    hmin = renvironment('SPLASH_HMIN_CODEUNITS',errval=-1.)
    if (hmin.gt.0.) then
       if (.not.allocated(dat)) then
          print*,' INTERNAL ERROR: dat not allocated in adjust_data_codeunits'
          return
       endif
       print*,' >> SETTING MINIMUM H TO ',hmin
       where (dat(:,ih,:) < hmin)
          dat(:,ih,:) = hmin
       end where
    endif
 endif
 
end subroutine adjust_data_codeunits

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
    print 10,'BIG-endian'
 else
    print 10,'LITTLE-endian'
 endif
10 format(' native byte order on this machine is ',a,/,&
          ' [byte order used to read data may be set by compiler flags/environment variables]')

end subroutine endian_info

end module getdata
