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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM
! THE NEXT GENERATION SPH CODE (sphNG)
!
! (also my Phantom SPH code which uses a similar format)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  COMMAND LINE FLAGS:
!
! --cm if set, then centre of mass is reset to origin
! --omega=3.142 if non-zero subtracts corotating velocities with omega as set
! --omegat=3.142 if non-zero subtracts corotating positions and velocities with omega as set
! --timeunits='yrs' sets default time units, either 's','min','hrs','yrs' or 'tfreefall'
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------
module sphNGread
 use params
 implicit none
 real(doub_prec) :: udist,umass,utime,umagfd
 real :: tfreefall,dtmax
 integer :: istartmhd,istartrt,nmhd,idivvcol,idivvxcol,icurlvxcol,icurlvycol,icurlvzcol,iHIIcol,iHeIIcol,iHeIIIcol
 integer :: nhydroreal4,istart_extra_real4
 integer :: itempcol = 0
 integer :: ncolstepfirst = 0
 integer :: nhydroarrays,nmhdarrays,ndustarrays,ndustlarge
 logical :: phantomdump,smalldump,mhddump,rtdump,usingvecp,igotmass,h2chem,rt_in_header
 logical :: usingeulr,cleaning
 logical :: batcode,tagged,debug
 integer, parameter :: maxarrsizes = 10
 integer, parameter :: maxinblock = 128 ! max allowed in each block
 integer, parameter :: lentag = 16
 character(len=lentag) :: tagarr(maxplot)
 integer, parameter :: itypemap_sink_phantom = 3
 integer, parameter :: itypemap_dust_phantom = 2
 integer, parameter :: itypemap_unknown_phantom = maxparttypes
 real, allocatable  :: grainsize(:)

 !------------------------------------------
 ! generic interface to utilities for tagged
 ! dump format
 !------------------------------------------
 interface extract
  module procedure extract_int, extract_real4, extract_real8, &
         extract_intarr, extract_real4arr, extract_real8arr
 end interface extract

contains

 !-------------------------------------------------------------------
 ! function mapping iphase setting in sphNG to splash particle types
 !-------------------------------------------------------------------
elemental integer function itypemap_sphNG(iphase)
 integer*1, intent(in) :: iphase

 select case(int(iphase))
 case(0)
    itypemap_sphNG = 1 ! gas
 case(11:)
    itypemap_sphNG = 2 ! dust
 case(1:9)
    itypemap_sphNG = 4 ! nptmass
 case(10)
    itypemap_sphNG = 5 ! star
 case default
    itypemap_sphNG = 6 ! unknown
 end select

end function itypemap_sphNG

 !---------------------------------------------------------------------
 ! function mapping iphase setting in Phantom to splash particle types
 !---------------------------------------------------------------------
elemental integer function itypemap_phantom(iphase)
 integer*1, intent(in) :: iphase

 select case(int(iphase))
 case(1:2)
    itypemap_phantom = iphase
 case(3:itypemap_unknown_phantom-1) ! put sinks as type 3, everything else shifted by one
    itypemap_phantom = iphase + 1
 case(-3) ! sink particles, either from external_binary or read from dump
    itypemap_phantom = itypemap_sink_phantom
 case default
    itypemap_phantom = itypemap_unknown_phantom
 end select

end function itypemap_phantom

 !------------------------------------------
 ! extraction of single integer variables
 !------------------------------------------
subroutine extract_int(tag,ival,intarr,tags,ntags,ierr,verbose)
 character(len=*),      intent(in)  :: tag
 integer,               intent(out) :: ival
 integer,               intent(in)  :: ntags,intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical, intent(in),   optional :: verbose
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ival = 0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i) then
          ival = intarr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(verbose)) &
    print "(a)",' WARNING: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_int

 !------------------------------------------
 ! extraction of single real*8 variables
 !------------------------------------------
subroutine extract_real8(tag,rval,r8arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real*8,                intent(out) :: rval
 real*8,                intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 rval = 0.d0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r8arr) >= i) then
          rval = r8arr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0) print "(a)",' WARNING: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real8

 !------------------------------------------
 ! extraction of single real*4 variables
 !------------------------------------------
subroutine extract_real4(tag,rval,r4arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real*4,                intent(out) :: rval
 real*4,                intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 rval = 0. ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r4arr) >= i) then
          rval = r4arr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0) print "(a)",' WARNING: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real4

 !------------------------------------------
 ! extraction of integer arrays
 !------------------------------------------
subroutine extract_intarr(tag,ival,intarr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 integer,               intent(out) :: ival(:)
 integer,               intent(in)  :: ntags,intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 ival(:) = 0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i .and. size(ival) > nmatched) then
          nmatched = nmatched + 1
          ival(nmatched) = intarr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(ival)) ierr = 0
 if (ierr /= 0) print "(a)",' WARNING: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_intarr

 !------------------------------------------
 ! extraction of real*8 arrays
 !------------------------------------------
subroutine extract_real8arr(tag,rval,r8arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real*8,                intent(out) :: rval(:)
 real*8,                intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 rval = 0.d0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r8arr) >= i .and. size(rval) > nmatched) then
          nmatched = nmatched + 1
          rval(nmatched) = r8arr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(rval)) ierr = 0
 if (ierr /= 0) print "(a)",' WARNING: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real8arr

 !------------------------------------------
 ! extraction of real*4 arrays
 !------------------------------------------
subroutine extract_real4arr(tag,rval,r4arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real*4,                intent(out) :: rval(:)
 real*4,                intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 rval = 0. ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r4arr) >= i .and. size(rval) > nmatched) then
          nmatched = nmatched + 1
          rval(nmatched) = r4arr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(rval)) ierr = 0
 if (ierr /= 0) print "(a)",' WARNING: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real4arr

 !----------------------------------------------------------------------
 ! Extract various options from the fileident string
 !----------------------------------------------------------------------
subroutine get_options_from_fileident(fileident,smalldump,tagged,phantomdump,&
                       usingvecp,usingeulr,cleaning,h2chem,rt_in_header,batcode)
 character(len=*), intent(in) :: fileident
 logical,          intent(out) :: smalldump,tagged,phantomdump,batcode
 logical,          intent(out) :: usingvecp,usingeulr,cleaning,h2chem,rt_in_header

 smalldump = .false.
 phantomdump = .false.
 usingvecp = .false.
 usingeulr = .false.
 cleaning = .false.
 h2chem = .false.
 rt_in_header = .false.
 batcode = .false.
 tagged = .false.
 if (fileident(1:1)=='S') then
    smalldump = .true.
 endif
 if (fileident(2:2)=='T') then
    tagged = .true.
 endif
 if (index(fileident,'Phantom') /= 0) then
    phantomdump = .true.
 else
    phantomdump = .false.
 endif
 if (index(fileident,'vecp') /= 0) then
    usingvecp = .true.
 endif
 if (index(fileident,'eulr') /= 0) then
    usingeulr = .true.
 endif
 if (index(fileident,'clean') /= 0) then
    cleaning = .true.
 endif
 if (index(fileident,'H2chem') /= 0) then
    h2chem = .true.
 endif
 if (index(fileident,'RT=on') /= 0) then
    rt_in_header = .true.
 endif
 if (index(fileident,'This is a test') /= 0) then
    batcode = .true.
 endif

end subroutine get_options_from_fileident

 !----------------------------------------------------------------------
 ! Set position of header items manually for older (untagged) formats
 !----------------------------------------------------------------------
subroutine fake_header_tags(nreals,realarr,tagsreal)
 integer, intent(in) :: nreals
 real,    intent(in) :: realarr(nreals)
 character(len=*), intent(out) :: tagsreal(:)
 integer, parameter :: ilocbinary = 24
 integer :: ipos
 !
 !--set the tags manually for older formats
 !  (but only the ones we care about)
 !
 if (phantomdump) then
    tagsreal(1) = 'time'
    tagsreal(3) = 'gamma'
    tagsreal(4) = 'rhozero'
    tagsreal(6) = 'hfact'
    tagsreal(7) = 'tolh'
    tagsreal(15:19) = 'massoftype'
    if (nreals >= ilocbinary + 14) then
       if (nreals >= ilocbinary + 15) then
          ipos = ilocbinary
       else
          print*,'*** WARNING: obsolete header format for external binary information ***'
          ipos = ilocbinary + 1
       endif
       if (debug) print*,'DEBUG: reading binary information from header ',ilocbinary
       if (any(realarr(ilocbinary:ilocbinary+14) /= 0.)) then
          tagsreal(ipos:ipos+2) = (/'x1','y1','z1'/)
          if (nreals >= ilocbinary+15) then
             tagsreal(ipos+3:ipos+9) = (/'m1','h1','x2','y2','z2','m2','h2'/)
             ipos = ipos + 10
          else
             tagsreal(ipos+3:ipos+7) = (/'h1','x2','y2','z2','h2'/)
             ipos = ipos + 8
          endif
          tagsreal(ipos:ipos+5) = (/'vx1','vy1','vz1','vx2','vy2','vz2'/)
       endif
    endif
 elseif (batcode) then
    tagsreal(1) = 'time'
    tagsreal(3) = 'gamma'
    tagsreal(4) = 'radL1'
    tagsreal(5) = 'PhiL1'
    tagsreal(15) = 'Er'
 else
    tagsreal(1) = 'gt'
    tagsreal(2) = 'dtmax'
    tagsreal(3) = 'gamma'
    tagsreal(4) = 'rhozero'
    tagsreal(5) = 'RK2'
    if (smalldump) then ! sphNG small dump
       if (nreals==15) then
          tagsreal(15) = 'pmassinitial'
       else
          tagsreal(23) = 'pmassinitial'
       endif
    endif
 endif

end subroutine fake_header_tags

subroutine set_grain_sizes(ntags,tags,vals,udist)
 integer, intent(in) :: ntags
 character(len=*), intent(in) :: tags(ntags)
 real, intent(inout) :: vals(ntags)
 real(doub_prec), intent(in) :: udist
 integer :: i,nd

 ! convert grain sizes to cm
 nd = 0
 do i=1,ntags
    if (index(tags(i),'grainsize') > 0) then
       nd = nd + 1
       vals(i) = vals(i)*udist
    endif
 enddo

 if (allocated(grainsize)) deallocate(grainsize)
 allocate(grainsize(nd))

 nd = 0
 do i=1,ntags
    if (index(tags(i),'grainsize') > 0) then
       nd = nd + 1
       grainsize(nd) = vals(i)
    endif
 enddo

end subroutine set_grain_sizes

!----------------------------------------------------------------------
! print information about dust grain sizes found in header
!----------------------------------------------------------------------
subroutine print_dustgrid_info(ntags,tags,vals,mgas)
 use asciiutils,     only:match_tag
 use labels,         only:get_label_grain_size
 integer, intent(in) :: ntags
 character(len=*), intent(in) :: tags(ntags)
 real, intent(in) :: vals(ntags),mgas
 integer :: i,nd

 nd = 0
 if (match_tag(tags,'grainsize1') > 0) then
    print "(/,a)",' Dust grid:'
    do i=1,ntags
       if (index(tags(i),'grainsize') > 0 .and. vals(i) > 0.) then
          nd = nd + 1
          print "(i3,a)",nd,': '//get_label_grain_size(vals(i))
       endif
    enddo
    if (nd > 0) print "(a)"
 endif

 ! nd = 0
 ! print *,' mgas = ',mgas/(2d33/umass)
 ! do i=1,ntags
 !    if (index(tags(i),'mdust_in') > 0) then
 !       nd = nd + 1
 !       if (vals(i) > 0.) print "(i3,a,1pg10.3,a,1pg10.3)",nd,': mass = ',vals(i)/(2d33/umass),&
 !                ' Msun; d/g = ',vals(i)/mgas
 !    endif
 ! enddo

end subroutine print_dustgrid_info

!----------------------------------------------------------------------
! Routine to read the header of sphNG dump files and extract relevant
! information
!----------------------------------------------------------------------
subroutine read_header(iunit,iverbose,debug,doubleprec,&
                        npart,npartoftypei,n1,ntypes,nblocks,&
                        narrsizes,realarr,tagsreal,nreals,ierr)
 use settings_data,  only:ndusttypes
 integer, intent(in)  :: iunit,iverbose
 logical, intent(in)  :: debug,doubleprec
 integer, intent(out) :: npart,npartoftypei(:),n1,ntypes,nblocks,narrsizes,nreals,ierr
 real,    intent(out) :: realarr(maxinblock)
 character(len=lentag), intent(out) :: tagsreal(maxinblock)
 character(len=lentag) :: tags(maxinblock)
 integer               :: intarr(maxinblock)
 real(doub_prec)       :: real8arr(maxinblock)
 real(sing_prec)       :: real4arr(maxinblock)
 integer               :: i,ierr1,ierr2,ierrs(4)
 integer               :: nints,ninttypes,nreal4s,nreal8s
 integer               :: n2,nreassign,naccrete,nkill

 ! initialise empty tag array
 tags(:) = ''
 intarr(:) = 0

 nblocks = 1 ! number of MPI blocks
 npartoftypei(:) = 0
 read(iunit,iostat=ierr) nints
 if (ierr /=0) then
    print "(a)",'error reading nints'
    return
 else
    if (tagged) then
       if (nints > maxinblock) then
          if (iverbose > 0) print*,'WARNING: number of ints in header exceeds splash array limit, ignoring some'
          nints = maxinblock
       endif
       read(iunit,iostat=ierr1) tags(1:nints)
       read(iunit,iostat=ierr2) intarr(1:nints)
       if (ierr1 /= 0 .or. ierr2 /= 0) then
          print "(a)",'error reading integer header'
          ierr = 1
          return
       endif
       if (debug) print*,'DEBUG: got tags = ',tags(1:nints)
       call extract('nblocks',nblocks,intarr,tags,nints,ierr)
       if (ierr /= 0) return
       call extract('nparttot',npart,intarr,tags,nints,ierr)
       if (ierr /= 0) return
       if (phantomdump) then
          call extract('ntypes',ntypes,intarr,tags,nints,ierr)
          if (ierr /= 0) return
          if (ntypes > maxparttypes) then
             if (iverbose > 0) &
                print "(a,i2)",' WARNING: number of particle types exceeds array limits: ignoring types > ',maxparttypes
             ntypes = maxparttypes
          endif
          call extract('npartoftype',npartoftypei(1:ntypes),intarr,tags,nints,ierr)
          if (ierr /= 0) return
       endif
       if (phantomdump .and. nints < 7) ntypes = nints - 1
       if (iverbose >= 2 .or. debug) print *,'npart = ',npart,' MPI blocks = ',nblocks
       if (phantomdump) then
          n1 = npartoftypei(1)
       else
          call extract('n1',n1,intarr,tags,nints,ierr)
       endif
    else
       if (nints < 3) then
          if (.not.phantomdump) print "(a)",'WARNING: npart,n1,n2 NOT IN HEADER??'
          read(iunit,iostat=ierr) npart
          npartoftypei(1) = npart
       elseif (phantomdump) then
          if (nints < 7) then
             ntypes = nints - 1
             read(iunit,iostat=ierr) npart,npartoftypei(1:ntypes)
          else
             ntypes = 5
             read(iunit,iostat=ierr) npart,npartoftypei(1:5),nblocks
          endif
          if (debug) then
             print*,'DEBUG: ntypes = ',ntypes,' npartoftype = ',npartoftypei(:)
          endif
          n1 = npartoftypei(1)
          n2 = 0
       elseif (nints >= 7) then
          read(iunit,iostat=ierr) npart,n1,n2,nreassign,naccrete,nkill,nblocks
       else
          print "(a)",'warning: nblocks not read from file (assuming non-MPI dump)'
          read(iunit,iostat=ierr) npart,n1,n2
       endif
       if (ierr /=0) then
          print "(a)",'error reading npart,n1,n2 and/or number of MPI blocks'
          return
       elseif (nblocks > 2000) then
          print *,'npart = ',npart,' MPI blocks = ',nblocks
          nblocks = 1
          print*,' corrupt number of MPI blocks, assuming 1 '
       else
          if (iverbose >= 1) print *,'npart = ',npart,' MPI blocks = ',nblocks
       endif
    endif
 endif

!--int*1, int*2, int*4, int*8
 ierr1 = 0
 ierr2 = 0
 do i=1,4
    read(iunit,iostat=ierr) ninttypes
    if (ninttypes > 0) then
       if (tagged) read(iunit,iostat=ierr1)
       read(iunit,iostat=ierr2)
    endif
    if (ierr /= 0 .or. ierr1 /= 0 .or. ierr2 /= 0) then
       print "(a)",'error skipping int types'
       return
    endif
 enddo
!--default reals
 read(iunit,iostat=ierr) nreals
 if (ierr /=0) then
    print "(a)",'error reading default reals'
    return
 else
    if (nreals > maxinblock) then
       print*,'WARNING: number of reals in header exceeds splash array limit, ignoring some'
       nreals = maxinblock
    endif
    if (tagged) read(iunit,iostat=ierr) tagsreal(1:nreals)
    if (doubleprec) then
       read(iunit,iostat=ierr) real8arr(1:nreals)
       realarr(1:nreals) = real(real8arr(1:nreals))
    else
       read(iunit,iostat=ierr) real4arr(1:nreals)
       realarr(1:nreals) = real(real4arr(1:nreals))
    endif
    if (.not.tagged) call fake_header_tags(nreals,realarr,tagsreal)
 endif
!
!--append integers to realarr so they can be used in
!  legends and calculated quantities
!
 if (nreals+nints <= maxinblock) then
    tagsreal(nreals+1:nreals+nints) = tags(1:nints)
    realarr(nreals+1:nreals+nints)  = real(intarr(1:nints))
    nreals = nreals + nints
 endif

!--real*4, real*8
 read(iunit,iostat=ierr) nreal4s
 if (nreal4s > 0) then
    if (tagged) read(iunit,iostat=ierr1)
    read(iunit,iostat=ierr2)
 endif
 if (ierr /= 0 .or. ierr1 /= 0 .or. ierr2 /= 0) then
    print "(a)",'error skipping real*4''s in header'
    return
 endif

 read(iunit,iostat=ierr) nreal8s
 if (ierr /= 0 .or. nreal8s < 0) then
    print "(a)",'error reading nreal8s'
    return
 endif
!   print "(a,i3)",' ndoubles = ',nreal8s
 if (iverbose >= 2 .or. debug) print "(4(a,i3),a)",' header contains ',nints,' ints, ',&
                           nreals,' reals,',nreal4s,' real4s, ',nreal8s,' doubles'
 if (tagged) then
    if (nreal8s > maxinblock) then
       print*,'WARNING: number of real8''s in header exceeds splash array limit, ignoring some'
       nreal8s = maxinblock
    endif
    read(iunit,iostat=ierr) tags(1:nreal8s)
    read(iunit,iostat=ierr) real8arr(1:nreal8s)
    call extract('udist',udist,real8arr,tags,nreal8s,ierrs(1))
    call extract('umass',umass,real8arr,tags,nreal8s,ierrs(2))
    call extract('utime',utime,real8arr,tags,nreal8s,ierrs(3))
    call extract('umagfd',umagfd,real8arr,tags,nreal8s,ierrs(4))

    ! extract the number of dustfrac arrays in the file
    ndusttypes = extract_ndusttypes(tags,tagsreal,intarr,nints)
!
!--append real*8s to realarr so they can be used in
!  legends and calculated quantities
!
    if (nreals+nreal8s <= maxinblock) then
       tagsreal(nreals+1:nreals+nreal8s) = tags(1:nreal8s)
       realarr(nreals+1:nreals+nreal8s)  = real(real8arr(1:nreal8s))
       nreals = nreals + nreal8s
    endif

 else
    if (nreal8s >= 4) then
       read(iunit,iostat=ierr) udist,umass,utime,umagfd
    elseif (nreal8s >= 3) then
       read(iunit,iostat=ierr) udist,umass,utime
       umagfd = 1.0
    else
       print "(a)",'*** WARNING: units not found in file'
       udist = 1.0
       umass = 1.0
       utime = 1.0
       umagfd = 1.0
    endif
 endif
 if (ierr /= 0) then
    print "(a)",'*** error reading units'
 endif
!
!--Total number of array blocks in the file
!
 read(iunit,iostat=ierr) narrsizes
 if (ierr /= 0) return
 if (debug) print*,' nblocks(total)=',narrsizes
 narrsizes = narrsizes/nblocks
 if (ierr /= 0) then
    print "(a)",'*** error reading number of array sizes ***'
    close(iunit)
    return
 elseif (narrsizes > maxarrsizes) then
    narrsizes = maxarrsizes
    print "(a,i2)",'WARNING: too many array sizes: reading only ',narrsizes
 endif
 if (narrsizes >= 4 .and. nreal8s < 4) then
    print "(a)",' WARNING: could not read magnetic units from dump file'
 endif
 if (debug) print*,' number of array sizes = ',narrsizes

end subroutine read_header

 !----------------------------------------------------------------------
 ! Read the header to each array block
 !----------------------------------------------------------------------
subroutine read_block_header(iunit,iblock,iarr,iverbose,debug,&
                              isize,nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8,&
                              ntotblock,npart,ntotal,nptmasstot,ncolstep,ierr)
 integer,   intent(in)    :: iunit,iblock,iarr,iverbose
 logical,   intent(in)    :: debug
 integer*8, intent(out)   :: isize(:)
 integer,   intent(out)   :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8,ierr
 integer,   intent(inout) :: ntotblock,npart,ntotal,nptmasstot,ncolstep

 read(iunit,iostat=ierr) isize(iarr),nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
 if (iarr==1) then
    ntotblock = isize(iarr)
    if (npart <= 0) npart = ntotblock
    ntotal = ntotal + ntotblock
 elseif (iarr==2) then
    nptmasstot = nptmasstot + isize(iarr)
 endif
 if (debug) print*,'DEBUG: array size ',iarr,' size = ',isize(iarr)
 if (isize(iarr) > 0 .and. iblock==1) then
    if (iverbose >= 2 .or. debug) print "(1x,a,i1,a,i12,a,5(i2,1x),a,3(i2,1x))", &
        'block ',iarr,' dim = ',isize(iarr),' nint =',nint,nint1,nint2,nint4,nint8,&
        'nreal =',nreal,nreal4,nreal8
 endif
!--we are going to read all real arrays but need to convert them all to default real
 if (iarr /= 2 .and. isize(iarr)==isize(1) .and. iblock==1) then
    ncolstep = ncolstep + nreal + nreal4 + nreal8
 endif

end subroutine read_block_header

 !----------------------------------------------------------------------
 ! Extract and print relevant variables from the header block
 !----------------------------------------------------------------------
subroutine extract_variables_from_header(tags,realarr,nreals,iverbose,debug,&
            gotbinary,nblocks,nptmasstot,npartoftypei,ntypes,&
            time,gamma,hfact,npart,ntotal,npartoftype,massoftype,dat,ix,ih,ipmass,ivx)
 character(len=lentag), intent(in) :: tags(maxinblock)
 real, intent(in) :: realarr(maxinblock)
 integer, intent(in) :: nreals,iverbose,nblocks,nptmasstot,npartoftypei(:),ntypes
 integer, intent(in) :: ix(3),ih,ipmass,ivx
 real, intent(out) :: time,gamma,hfact,massoftype(:)
 real, intent(inout) :: dat(:,:)
 integer, intent(out) :: npartoftype(:)
 integer, intent(inout) :: npart,ntotal
 logical, intent(in)  :: debug
 logical, intent(out) :: gotbinary
 real :: rhozero,tff,radL1,PhiL1,Er,RK2 !,dtmax
 real :: massoftypei(ntypes)
 integer :: i,ierrs(10)
 integer :: itype
 real, parameter :: pi=4.*atan(1.)

 if (phantomdump) then
    call extract('time',time,realarr,tags,nreals,ierrs(1))
 else
    call extract('gt',time,realarr,tags,nreals,ierrs(1))
 endif
 call extract('gamma',gamma,realarr,tags,nreals,ierrs(2))
 call extract('rhozero',rhozero,realarr,tags,nreals,ierrs(3))

!--extract required information from the first block header
 if (rhozero > 0.) then
    tfreefall = sqrt((3. * pi) / (32. * rhozero))
    tff = time/tfreefall
 else
    tfreefall = 0.
    tff = 0.
 endif
 if (phantomdump) then
    call extract('massoftype',massoftypei(1:ntypes),realarr,tags,nreals,ierrs(4))
    npartoftype(:) = 0
    do i=1,ntypes !--map from phantom types to splash types
       itype = itypemap_phantom(int(i,kind=1))
       if (debug) print*,'DEBUG: npart of type ',itype,' += ',npartoftypei(i)
       npartoftype(itype) = npartoftype(itype) + npartoftypei(i)
       massoftype(itype)  = massoftypei(i)
    enddo
    npartoftype(itypemap_sink_phantom) = nptmasstot  ! sink particles
    if (debug) print*,'DEBUG: npart of type sink ',itypemap_sink_phantom,' = ',nptmasstot
    !
    !--if Phantom calculation uses the binary potential
    !  then read this as two point mass particles
    !
    if (any(tags(1:nreals)=='x1')) then
       gotbinary = .true.
       npartoftype(itypemap_sink_phantom) = npartoftype(itypemap_sink_phantom) + 2
       ntotal = ntotal + 2
       call extract('x1',dat(npart+1,ix(1)),realarr,tags,nreals,ierrs(1))
       call extract('y1',dat(npart+1,ix(2)),realarr,tags,nreals,ierrs(2))
       call extract('z1',dat(npart+1,ix(3)),realarr,tags,nreals,ierrs(3))
       if (ipmass > 0) call extract('m1',dat(npart+1,ipmass),realarr,tags,nreals,ierrs(4))
       call extract('h1',dat(npart+1,ih),realarr,tags,nreals,ierrs(4))
       if (debug) print *,npart+1,npart+2
       if (iverbose >= 1) print *,'binary position:   primary: ',dat(npart+1,ix(1):ix(3))

       call extract('x2',dat(npart+2,ix(1)),realarr,tags,nreals,ierrs(1))
       call extract('y2',dat(npart+2,ix(2)),realarr,tags,nreals,ierrs(2))
       call extract('z2',dat(npart+2,ix(3)),realarr,tags,nreals,ierrs(3))
       if (ipmass > 0) call extract('m2',dat(npart+2,ipmass),realarr,tags,nreals,ierrs(4))
       call extract('h2',dat(npart+2,ih),realarr,tags,nreals,ierrs(5))
       if (iverbose >= 1 .and. ipmass > 0 .and. ih > 0) then
          print *,'                 secondary: ',dat(npart+2,ix(1):ix(3))
          print *,' m1: ',dat(npart+1,ipmass),' m2:',dat(npart+2,ipmass),&
                   ' h1: ',dat(npart+2,ipmass),' h2:',dat(npart+2,ih)
       endif
       if (ivx > 0) then
          call extract('vx1',dat(npart+1,ivx),realarr,tags,nreals,ierrs(1))
          call extract('vy1',dat(npart+1,ivx+1),realarr,tags,nreals,ierrs(2))
          call extract('vz1',dat(npart+1,ivx+2),realarr,tags,nreals,ierrs(3))
          call extract('vx2',dat(npart+2,ivx),realarr,tags,nreals,ierrs(4))
          call extract('vy2',dat(npart+2,ivx+1),realarr,tags,nreals,ierrs(5))
          call extract('vz2',dat(npart+2,ivx+2),realarr,tags,nreals,ierrs(6))
       endif
       npart  = npart  + 2
    endif
 else
    npartoftype(:) = 0
    npartoftype(1) = npart
    npartoftype(2) = max(ntotal - npart,0)
 endif
 hfact = 1.2
 if (phantomdump) then
    call extract('hfact',hfact,realarr,tags,nreals,ierrs(1))
    call extract('dtmax',dtmax,realarr,tags,nreals,ierrs(2))
    if (iverbose > 0) then
       print "(a,es12.4,a,f6.3,a,f5.2)", &
           ' time = ',time,' gamma = ',gamma,' hfact = ',hfact
    endif
 elseif (batcode) then
    call extract('radL1',radL1,realarr,tags,nreals,ierrs(1))
    call extract('PhiL1',PhiL1,realarr,tags,nreals,ierrs(2))
    call extract('Er',Er,realarr,tags,nreals,ierrs(3))
    if (iverbose > 0) then
       print "(a,es12.4,a,f9.5,a,f8.4,/,a,es12.4,a,es9.2,a,es10.2)", &
           '   time: ',time,  '   gamma: ',gamma, '   tsph: ',realarr(2), &
           '  radL1: ',radL1,'   PhiL1: ',PhiL1,'     Er: ',Er
    endif
 else
    call extract('RK2',RK2,realarr,tags,nreals,ierrs(1))
    call extract('dtmax',dtmax,realarr,tags,nreals,ierrs(2))
    if (iverbose > 0) then
       print "(a,es12.4,a,f9.5,a,f8.4,/,a,es12.4,a,es9.2,a,es10.2)", &
           '   time: ',time,  '   gamma: ',gamma, '   RK2: ',RK2, &
           ' t/t_ff: ',tff,' rhozero: ',rhozero,' dtmax: ',dtmax
    endif
 endif

end subroutine extract_variables_from_header

!---------------------------------------------------------------
! old subroutine for guessing labels in non-tagged sphNG format
!---------------------------------------------------------------
subroutine guess_labels(ncolumns,iamvec,label,labelvec,istartmhd, &
             istart_extra_real4,nmhd,nhydroreal4,ndimV,irho,iBfirst,ivx,&
             iutherm,idivB,iJfirst,iradenergy,icv,udist,utime,units,&
             unitslabel)
 use geometry, only:labelcoord
 use labels,   only:make_vector_label
 integer, intent(in) :: ncolumns,istartmhd,istart_extra_real4
 integer, intent(in) :: nmhd,nhydroreal4,ndimV,irho
 integer, intent(out) :: iBfirst,ivx,iutherm,idivB,iJfirst,iradenergy,icv
 integer, intent(inout) :: iamvec(:)
 character(len=*), intent(inout) :: label(:),labelvec(:),unitslabel(:)
 real(doub_prec),  intent(in)    :: udist,utime
 real,    intent(inout) :: units(:)
 real(doub_prec) :: uergg

!--the following only for mhd small dumps or full dumps
 if (ncolumns >= 7) then
    if (mhddump) then
       iBfirst = irho+1
       if (.not.smalldump) then
          ivx = iBfirst+ndimV
          iutherm = ivx+ndimV

          if (phantomdump) then
             !--phantom MHD full dumps
             if (nmhd >= 4) then
                call make_vector_label('A',istartmhd,ndimV,iamvec,labelvec,label,labelcoord(:,1))
                if (nmhd >= 7) then
                   label(istartmhd+3) = 'Euler beta_{x}'
                   label(istartmhd+4) = 'Euler beta_{x}'
                   label(istartmhd+5) = 'Euler beta_{y}'
                   idivB = istartmhd+2*ndimV
                else
                   idivB = istartmhd+ndimV
                endif
             elseif (nmhd >= 3) then
                label(istartmhd) = 'Euler alpha'
                label(istartmhd+1) = 'Euler beta'
                idivB = istartmhd + 2
             elseif (nmhd >= 2) then
                label(istartmhd) = 'Psi'
                idivB = istartmhd + 1
             elseif (nmhd >= 1) then
                idivB = istartmhd
             endif
             iJfirst = 0
             if (ncolumns >= idivB+1) then
                label(idivB+1) = 'alpha_{B}'
             endif

          else
             !--sphNG MHD full dumps
             label(iutherm+1) = 'grad h'
             label(iutherm+2) = 'grad soft'
             label(iutherm+3) = 'alpha'
             if (nmhd >= 7 .and. usingvecp) then
                call make_vector_label('A',istartmhd,ndimV,iamvec,labelvec,label,labelcoord(:,1))
                idivB = istartmhd+ndimV
             elseif (nmhd >= 6 .and. usingeulr) then
                label(istartmhd) = 'Euler alpha'
                label(istartmhd+1) = 'Euler beta'
                idivB = istartmhd + 2
             elseif (nmhd >= 6) then
                label(istartmhd) = 'psi'
                idivB = istartmhd + 1
                if (nmhd >= 8) then
                   label(istartmhd+2+ndimV+1) = '\eta_{real}'
                   label(istartmhd+2+ndimV+2) = '\eta_{art}'
                   units(istartmhd+2+ndimV+1:istartmhd+2+ndimV+2) = udist*udist/utime
                   unitslabel(istartmhd+2+ndimV+1:istartmhd+2+ndimV+2) = ' [cm^2/s]'
                endif
                if (nmhd >= 14) then
                   call make_vector_label('fsym',istartmhd+2+ndimV+3,ndimV,&
                         iamvec,labelvec,label,labelcoord(:,1))
                   call make_vector_label('faniso',istartmhd+2+ndimV+6,ndimV,&
                         iamvec,labelvec,label,labelcoord(:,1))
                endif
             elseif (nmhd >= 1) then
                idivB = istartmhd
             endif
             iJfirst = idivB + 1
             if (ncolumns >= iJfirst+ndimV) then
                label(iJfirst+ndimV) = 'alpha_{B}'
             endif
          endif
       else ! mhd small dump
          if (nhydroreal4 >= 3) iutherm = iBfirst+ndimV
       endif
    elseif (.not.smalldump) then
       ! pure hydro full dump
       ivx = irho+1
       iutherm = ivx + ndimV
       if (phantomdump) then
          if (istart_extra_real4 > 0 .and. istart_extra_real4 < 100) then
             label(istart_extra_real4) = 'alpha'
             label(istart_extra_real4+1) = 'alphau'
          endif
       else
          if (istart_extra_real4 > 0 .and. istart_extra_real4 < 100) then
             label(istart_extra_real4) = 'grad h'
             label(istart_extra_real4+1) = 'grad soft'
             label(istart_extra_real4+2) = 'alpha'
          endif
       endif
    endif

    if (phantomdump .and. h2chem) then
       if (smalldump) then
          label(nhydroarrays+nmhdarrays+1) = 'H_2 ratio'
       elseif (.not.smalldump .and. iutherm > 0) then
          label(iutherm+1) = 'H_2 ratio'
          label(iutherm+2) = 'HI abundance'
          label(iutherm+3) = 'proton abundance'
          label(iutherm+4) = 'e^- abundance'
          label(iutherm+5) = 'CO abundance'
       endif
    endif
    if (istartrt > 0 .and. istartrt <= ncolumns .and. rtdump) then ! radiative transfer dump
       iradenergy = istartrt
       label(iradenergy) = 'radiation energy'
       uergg = (udist/utime)**2
       units(iradenergy) = uergg
       if (smalldump) then
          icv = istartrt+1
       else
          label(istartrt+1) = 'opacity'
          units(istartrt+1) = udist**2/umass
          icv = istartrt+2
          label(istartrt+3) = 'lambda'
          units(istartrt+3) = 1.0

          label(istartrt+4) = 'eddington factor'
          units(istartrt+4) = 1.0
       endif
       if (icv > 0) then
          label(icv) = 'u/T'
          units(icv) = uergg
       endif
    else
       iradenergy = 0
       icv = 0
    endif
 endif
end subroutine guess_labels

integer function assign_column(tag,iarr,ipos,ikind,imaxcolumnread,idustarr,ncolstep) result(icolumn)
 use labels, only:ih,irho,ix,ipmass
 character(len=lentag), intent(in) :: tag
 integer,               intent(in) :: iarr,ipos,ikind,ncolstep
 integer,               intent(inout) :: imaxcolumnread
 integer,               intent(inout) :: idustarr

 if (tagged .and. len_trim(tag) > 0) then
    !
    ! use the tags to put certain arrays in an assigned place
    ! no matter what type is used for the variable in the file
    ! and no matter what order they appear in the dump file
    !
    select case(trim(tag))
    case('x')
       icolumn = ix(1)
    case('y')
       icolumn = ix(2)
    case('z')
       icolumn = ix(3)
    case('m')
       icolumn = ipmass
    case('h')
       icolumn = ih
    case('rho')
       icolumn = irho
    case('dustfracsum')
       idustarr = idustarr + 1
       icolumn = nhydroarrays + idustarr
    case('dustfrac')
       if (ndustarrays > 0) then
          idustarr = idustarr + 1
          icolumn = nhydroarrays + idustarr
       else
          icolumn = max(nhydroarrays + ndustarrays + nmhdarrays + 1,imaxcolumnread + 1)
       endif
    case('Bx')
       icolumn = nhydroarrays + ndustarrays + 1
    case('By')
       icolumn = nhydroarrays + ndustarrays + 2
    case('Bz')
       icolumn = nhydroarrays + ndustarrays + 3
    case default
       icolumn = max(nhydroarrays + ndustarrays + nmhdarrays + 1,imaxcolumnread + 1)
       if (icolumn > ncolstep) then
          ! check for dustfrac not being present
          if (idustarr == 0 .and. ndustarrays > 0) icolumn = nhydroarrays + 1
       endif
       if (iarr==1) then
          if (ikind==4) then  ! real*4 array
             istart_extra_real4 = min(istart_extra_real4,icolumn)
             if (debug) print*,' istart_extra_real4 = ',istart_extra_real4
          endif
       endif
    end select
 else
    !
    ! this is old code handling the non-tagged format where
    ! particular arrays are assumed to be in particular places
    !
    if (ikind==6) then   ! default reals
       if (iarr==1.and.((phantomdump.and.ipos==4) &
           .or.(.not.phantomdump.and.ipos==6))) then
          ! read x,y,z,m,h and then place arrays after always-present ones
          ! (for phantom read x,y,z only)
          icolumn = nhydroarrays+nmhdarrays + 1
       elseif (.not.phantomdump .and. (iarr==4 .and. ipos <= 3)) then
          icolumn = nhydroarrays + ipos
       else
          icolumn = imaxcolumnread + 1
       endif
    elseif (ikind==4) then  ! real*4s
       if (phantomdump) then
          if (iarr==1 .and. ipos==1) then
             icolumn = ih ! h is always first real4 in phantom dumps
             !!--density depends on h being read
             !required(ih) = .true.
          elseif (iarr==4 .and. ipos <= 3) then
             icolumn = nhydroarrays + ipos
          else
             icolumn = max(nhydroarrays+nmhdarrays + 1,imaxcolumnread + 1)
             if (iarr==1) then
                istart_extra_real4 = min(istart_extra_real4,icolumn)
                if (debug) print*,' istart_extra_real4 = ',istart_extra_real4
             endif
          endif
       else
          if (iarr==1 .and. ipos==1) then
             icolumn = irho ! density
          elseif (iarr==1 .and. smalldump .and. ipos==2) then
             icolumn = ih ! h which is real4 in small dumps
             !--this was a bug for sphNG files...
             !elseif (iarr==4 .and. i <= 3) then
             !   icolumn = nhydroarrays + ipos
          else
             icolumn = max(nhydroarrays+nmhdarrays + 1,imaxcolumnread + 1)
             if (iarr==1) then
                istart_extra_real4 = min(istart_extra_real4,icolumn)
                if (debug) print*,' istart_extra_real4 = ',istart_extra_real4
             endif
          endif
       endif
    else ! used for untagged format with real*8's
       icolumn = imaxcolumnread + 1
    endif
 endif
 if (icolumn < 1 .or. icolumn > ncolstep) then
    print*,' ERROR in column assignment for '//trim(tag)
    icolumn = ncolstep ! screw up last entry, but don't seg fault
 endif

 imaxcolumnread = max(imaxcolumnread,icolumn)

end function assign_column

!---------------------------------------------------------------
! function to extract the number of dust arrays
!---------------------------------------------------------------
integer function extract_ndusttypes(tags,tagsreal,intarr,nints) result(ndusttypes)
 character(len=lentag), intent(in) :: tags(maxinblock),tagsreal(maxinblock)
 integer, intent(in) :: intarr(:),nints
 integer :: i,idust,ierr,ndustsmall
 logical :: igotndusttypes

 ! Look for ndusttypes in the header
 igotndusttypes = .false.
 do i = 1,maxinblock
    if (trim(tags(i))=='ndusttypes') igotndusttypes = .true.
 enddo

 ! Retreive/guess the value of ndusttypes
 if (igotndusttypes) then
    call extract('ndusttypes',idust,intarr,tags,nints,ierr)
 else
    call extract('ndustsmall',ndustsmall,intarr,tags,nints,ierr,verbose=.false.)
    call extract('ndustlarge',ndustlarge,intarr,tags,nints,ierr,verbose=.false.)
    idust = ndustsmall+ndustlarge
    ! For older files where ndusttypes is not output to the header
    if (ierr /= 0) then
       idust = 0
       do i = 1,maxinblock
          if (tagsreal(i)=='grainsize') idust = idust + 1
       enddo
       if (idust > 0) then
          write(*,"(a)")    ' Warning! Could not find ndusttypes in header'
          write(*,"(a,I4)") '          ...counting grainsize arrays...ndusttypes =',idust
       endif
    endif
 endif
 ndusttypes = idust

end function extract_ndusttypes

subroutine get_rho_from_h(i1,i2,ih,ipmass,irho,required,npartoftype,massoftype,hfact,dat,iphase,nkilled)
 integer,            intent(in)    :: i1,i2,ih,ipmass,irho
 logical,            intent(in)    :: required(0:)
 integer,            intent(inout) :: npartoftype(:)
 real,               intent(in)    :: massoftype(:),hfact
 real,               intent(inout) :: dat(:,:)
 integer(kind=int1), intent(inout) :: iphase(:)
 integer,            intent(inout) :: nkilled
 integer :: itype,k
 real :: pmassi,hi,rhoi
 !
 !--dead particles have -ve smoothing lengths in phantom
 !  so use abs(h) for these particles and hide them
 !
 if (any(npartoftype(2:) > 0)) then
    if (.not.required(ih)) print*,'ERROR: need to read h, but required=F'
    !--need masses for each type if not all gas
    if (debug) print*,'DEBUG: phantom: setting h for multiple types ',i1,i2
    if (debug) print*,'DEBUG: massoftype = ',massoftype(:)
    do k=i1,i2
       itype = itypemap_phantom(iphase(k))
       pmassi = massoftype(itype)
       hi = dat(k,ih)
       if (ipmass > 0) pmassi = dat(k,ipmass)
       if (hi > 0.) then
          if (required(irho)) dat(k,irho) = pmassi*(hfact/hi)**3
       elseif (hi < 0.) then
          !print*,' accreted: ',k,' type was ',itype,iphase(k)
          npartoftype(itype) = npartoftype(itype) - 1
          npartoftype(itypemap_unknown_phantom) = npartoftype(itypemap_unknown_phantom) + 1
          if (required(irho)) dat(k,irho) = pmassi*(hfact/abs(hi))**3
          iphase(k) = -1
       else ! dead particles
          npartoftype(itype) = npartoftype(itype) - 1
          npartoftype(itypemap_unknown_phantom) = npartoftype(itypemap_unknown_phantom) + 1
          nkilled = nkilled + 1
          if (required(irho)) dat(k,irho) = 0.
          iphase(k) = -2
       endif
   enddo
else
   if (.not.required(ih)) print*,'ERROR: need to read h, but required=F'
   if (debug) print*,'debug: phantom: setting rho for all types'
   pmassi = massoftype(1)
   !--assume all particles are gas particles
   do k=i1,i2
      hi = dat(k,ih)
      if (ipmass > 0) pmassi = dat(k,ipmass)
      if (hi > 0.) then
         rhoi = pmassi*(hfact/hi)**3
      elseif (hi < 0.) then
         rhoi = pmassi*(hfact/abs(hi))**3
         npartoftype(1) = npartoftype(1) - 1
         npartoftype(itypemap_unknown_phantom) = npartoftype(itypemap_unknown_phantom) + 1
         iphase(k) = -1
      else ! if h = 0.
         rhoi = 0.
         npartoftype(1) = npartoftype(1) - 1
         npartoftype(itypemap_unknown_phantom) = npartoftype(itypemap_unknown_phantom) + 1
         iphase(k) = -2
      endif
      if (required(irho)) dat(k,irho) = rhoi
   enddo
endif

end subroutine get_rho_from_h

!----------------------------------------------------------------------
!  Set a negative smoothing length for merged sinks, so that
!  they can be ignored when plotting
!----------------------------------------------------------------------
subroutine set_sink_merged(i1,i2,ih,ipmass,dat)
 integer, intent(in) :: i1,i2,ih,ipmass
 real, intent(inout) :: dat(:,:)
 integer :: i

 if (ih > 0 .and. ipmass > 0) then
    do i=i1,i2
       if (dat(i,ipmass) < 0.) dat(i,ih) = -1.
    enddo
 endif

end subroutine set_sink_merged

!----------------------------------------------------------------------
!  Set density on sink particles based on the mass and radius
!  this is useful for opacity rendering, but also provides useful
!  information rather than just having zero density on sinks
!----------------------------------------------------------------------
subroutine set_sink_density(i1,i2,ih,ipmass,irho,dat)
 integer, intent(in) :: i1,i2,ih,ipmass,irho
 real, intent(inout) :: dat(:,:)
 integer :: i

 if (ih > 0 .and. ipmass > 0 .and. irho > 0) then
    do i=i1,i2
       if (dat(i,ih) > 0.) dat(i,irho) = dat(i,ipmass)/dat(i,ih)**3
    enddo
 endif

end subroutine set_sink_density

!----------------------------------------------------------------------
!  Map sink particle data to splash columns
!----------------------------------------------------------------------
integer function map_sink_property_to_column(k,ilocvx,ncolmax) result(iloc)
 use labels, only:ix,ipmass,ih,ivx
 integer, intent(in) :: k,ilocvx,ncolmax

 select case(k)
 case(1:3)
     iloc = ix(k)
 case(4)
     iloc = ipmass
 case(5)
     iloc = ih
 case default
     if (k >= ilocvx .and. k < ilocvx+3 .and. ivx > 0) then
        iloc = ivx + k-ilocvx ! put velocity into correct arrays
     else
        iloc = 0
     endif
 end select
 if (iloc > ncolmax) iloc = 0  ! error occurred

end function map_sink_property_to_column

!------------------------------------------------------------
! sanity check of the particle type accounting
!------------------------------------------------------------
subroutine check_iphase_matches_npartoftype(i1,i2,iphase,npartoftypei)
 use labels, only:labeltype
 use params, only:int1
 integer, intent(in) :: i1,i2
 integer(kind=int1), intent(in) :: iphase(i1:i2)
 integer, intent(inout) :: npartoftypei(:)
 integer :: npartoftype_new(size(npartoftypei))
 integer :: k,itype

 npartoftype_new(:) = 0
 do k=i1,i2
    itype = itypemap_phantom(iphase(k))
    npartoftype_new(itype) = npartoftype_new(itype) + 1
 enddo
 do k=1,size(npartoftypei)
    if (npartoftype_new(k) /= npartoftypei(k)) then
       print*,' WARNING: got ',npartoftype_new(k),&
              trim(labeltype(k))//' particles, expecting ',npartoftypei(k)
       npartoftypei(k) = npartoftype_new(k)
    endif
 enddo

end subroutine check_iphase_matches_npartoftype

!------------------------------------------------------------
! allocate and reallocate the iphase array as needed
!------------------------------------------------------------
subroutine allocate_iphase(iphase,nmax,phantomdump,gotbinary,nlocbinary)
 integer*1, allocatable, intent(inout) :: iphase(:)
 integer, intent(in) :: nmax,nlocbinary
 logical, intent(in) :: phantomdump,gotbinary
 integer*1, allocatable :: iphase_old(:)
 integer :: ncopy

 if (allocated(iphase)) then
    iphase_old = iphase ! allocates memory for iphase_old as well
    deallocate(iphase)
 endif
 allocate(iphase(nmax))

 if (phantomdump) then
    iphase(:) = 1
 else
    iphase(:) = 0
 endif

 if (gotbinary) then
    iphase(nlocbinary) = -3
    iphase(nlocbinary+1) = -3
 endif

 if (allocated(iphase_old)) then
    ncopy = min(size(iphase_old),size(iphase))
    iphase(1:ncopy) = iphase_old(1:ncopy)
    deallocate(iphase_old)
 endif

end subroutine allocate_iphase

end module sphNGread

!----------------------------------------------------------------------
!  Module for read_data_sphNG and set_labels_sphNG routines
!----------------------------------------------------------------------

module readdata_sphNG
 implicit none

 public :: read_data_sphNG, set_labels_sphNG, file_format_is_sphNG

 private
contains

!----------------------------------------------------------------------
!  Main read_data_sphNG routine for splash
!----------------------------------------------------------------------
subroutine read_data_sphNG(rootname,indexstart,iposn,nstepsread)
 use particle_data,  only:dat,gamma,time,headervals,&
                      iamtype,npartoftype,maxpart,maxstep,maxcol,masstype
 !use params,         only:int1,int8
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,&
                      lowmemorymode,ntypes,iverbose,ndusttypes
 use mem_allocation, only:alloc
 use system_utils,   only:lenvironment,renvironment
 use labels,         only:ipmass,irho,ih,ix,ivx,labeltype,print_types,headertags,&
                          iutherm,itemp,ikappa,irhorestframe,labelreq,nreq
 use calcquantities, only:calc_quantities
 use asciiutils,     only:make_tags_unique,match_tag
 use sphNGread
 use lightcurve_utils, only:get_temp_from_u,ionisation_fraction,get_opacity
 use read_kepler,      only:check_for_composition_file,read_kepler_composition
 integer, intent(in)  :: indexstart,iposn
 integer, intent(out) :: nstepsread
 character(len=*), intent(in) :: rootname
 integer :: i,j,k,ierr,iunit
 integer :: intg1,int2,int3,ilocvx,iversion
 integer :: i1,iarr,i2,iptmass1,iptmass2
 integer :: npart_max,nstep_max,ncolstep,icolumn,idustarr,nptmasstot
 integer :: narrsizes
 integer :: nskip,ntotal,npart,n1,ngas,nreals
 integer :: iblock,nblocks,ntotblock,ncolcopy
 integer :: ipos,nptmass,nptmassi,ndust,nstar,nunknown,ilastrequired
 integer :: imaxcolumnread,nhydroarraysinfile,nhdr,nkilled
 integer :: itype,iphaseminthistype,iphasemaxthistype,nthistype,iloc,idenscol
 integer :: icentre,icomp_col_start,ncomp
 integer, dimension(maxparttypes) :: npartoftypei
 real,    dimension(maxparttypes) :: massoftypei
 logical :: iexist, doubleprec,imadepmasscolumn,gotbinary,gotiphase

 character(len=len(rootname)+10) :: dumpfile,compfile
 character(len=100) :: fileident

 integer*8, dimension(maxarrsizes) :: isize
 integer, dimension(maxarrsizes) :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
 integer*1, dimension(:), allocatable :: iphase
 integer, dimension(:), allocatable :: listpm,level
 real(doub_prec), dimension(:), allocatable :: dattemp
 real*4, dimension(:), allocatable :: dattempsingle,massfac
 real(doub_prec) :: r8,unit_dens,unit_ergg
 real(sing_prec) :: r4
 real, dimension(:,:), allocatable :: dattemp2
 real, dimension(maxinblock) :: dummyreal
 real :: hfact,omega
 real(doub_prec) :: Xfrac,Yfrac
 real :: xHIi,xHIIi,xHeIi,xHeIIi,xHeIIIi,nei
 logical :: skip_corrupted_block_3,get_temperature,get_kappa,get_ionfrac,need_to_allocate_iphase
 character(len=lentag) :: tagsreal(maxinblock), tagtmp

 integer, parameter :: splash_max_iversion = 1
 real, parameter :: Xfrac_default=0.69843,Yfrac_default=0.28731

 nstepsread = 0
 nstep_max = 0
 npart_max = maxpart
 npart = 0
 iunit = 15
 ipmass = 4
 idivvcol = 0
 idivvxcol = 0
 icurlvxcol = 0
 icurlvycol = 0
 icurlvzcol = 0
 iHIIcol = 0
 iHeIIcol = 0
 iHeIIIcol = 0
 nhydroreal4 = 0
 umass = 1.d0
 utime = 1.d0
 udist = 1.d0
 umagfd = 1.d0
 istartmhd = 0
 istartrt  = 0
 istart_extra_real4 = 100
 nmhd      = 0
 igotmass    = .false.
 tfreefall   = 1.d0
 gotbinary   = .false.
 gotiphase   = .false.
 skip_corrupted_block_3 = .false.

 dumpfile = trim(rootname)
 !
 !--check if data file exists
 !
 inquire(file=dumpfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
    return
 endif
 !
 !--fix number of spatial dimensions
 !
 ndim = 3
 ndimV = 3

 j = indexstart
 nstepsread = 0
 doubleprec = .true.
 get_temperature = lenvironment("SPLASH_GET_TEMP")
 get_kappa = lenvironment("SPLASH_GET_KAPPA")
 get_ionfrac = lenvironment("SPLASH_GET_ION")
 if ((get_temperature .or. get_kappa) .and. itempcol > 0 .and. required(itempcol)) then
    required(irho) = .true.
    required(irhorestframe) = .true.
    required(iutherm) = .true.
 endif
 if (get_ionfrac .or. get_kappa) then
    required(irho) = .true.
    required(itemp) = .true.
    Xfrac = renvironment("SPLASH_XFRAC",Xfrac_default)
    Yfrac = renvironment("SPLASH_YFRAC",Yfrac_default)

    if ( Xfrac < 0. .or. Xfrac > 1.) then
       Xfrac = Xfrac_default
       print "(1x,a,f5.3)",'ERROR: Input Xfrac is not between 0 and 1, using default value of ',Xfrac
    endif
    if ( Yfrac < 0. .or. Yfrac > 1.) then
       Yfrac = Yfrac_default
       print "(1x,a,f5.3)",'ERROR: Input Yfrac is not between 0 and 1, using default value of ',Yfrac
    endif
    if ( Xfrac + Yfrac > 1.) then
      print "(1x,a,f5.3,a,f5.3,a,f5.3)",'ERROR: Xfrac + Yfrac = ',Xfrac+Yfrac,' exceeds 1. Using default values of Xfrac = ',&
         Xfrac_default,' and Yfrac = ',Yfrac_default
      Xfrac = Xfrac_default
      Yfrac = Yfrac_default
    endif
 endif
 ilastrequired = 0
 do i=1,size(required)-1
    if (required(i)) ilastrequired = i
 enddo

 if (iverbose >= 1) print "(1x,a)",'reading sphNG format'
 write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

 debug = lenvironment('SSPLASH_DEBUG')
 if (debug) iverbose = 1
!
!--open the (unformatted) binary file
!
 open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='unformatted')
 if (ierr /= 0) then
    print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
    return
 else
    !
    !--read header key to work out precision
    !
    doubleprec = .true.
    read(iunit,iostat=ierr) intg1,r8,int2,iversion,int3
    if (intg1 /= 690706 .and. intg1 /= 060769) then
       print "(a)",'*** ERROR READING HEADER: corrupt file/zero size/wrong endian?'
       close(iunit)
       return
    endif
    if (int2 /= 780806 .and. int2 /= 060878) then
       if (iverbose >= 2) print "(a)",' single precision dump'
       rewind(iunit)
       read(iunit,iostat=ierr) intg1,r4,int2,iversion,int3
       if (int2 /= 780806 .and. int2 /= 060878) then
          print "(a)",'ERROR determining single/double precision in file header'
       endif
       doubleprec = .false.
    elseif (int3 /= 690706) then
       print*,' got ',intg1,r4,int2,iversion,int3
       print "(a)",'*** WARNING: default int appears to be int*8: not implemented'
    else
       if (debug) print "(a)",' double precision dump' ! no need to print this
    endif
    if (iversion==690706) then ! handle old-format files (without version number) gracefully
       iversion = 0
    endif
 endif
 if (iversion > splash_max_iversion) then
    print "(/a,i2,/,a,i2)",&
        ' *** WARNING: this copy of splash can only read version ',splash_max_iversion, &
        '              but the file format version is ',iversion
    if (.not.lenvironment('SSPLASH_IGNORE_IVERSION')) then
       print "(2(/,a))",'   ** press any key to bravely proceed anyway ** ', &
                          '   (set SSPLASH_IGNORE_IVERSION=yes to silence this warning)'
       read*
    endif
 endif
!
!--read file ID
!
 read(iunit,iostat=ierr) fileident
 if (ierr /=0) then
    print "(a)",'*** ERROR READING FILE ID ***'
    close(iunit)
    return
 else
    if (iverbose >= 0) print "(1x,a)",trim(fileident)
 endif
 mhddump = .false.
 rtdump = .false.
 call get_options_from_fileident(fileident,smalldump,tagged,phantomdump,&
                                   usingvecp,usingeulr,cleaning,h2chem,rt_in_header,batcode)
 if (tagged .and. iversion < 1) print "(a)",'ERROR: got tagged format but iversion is ',iversion
!
!--read variables from header
!
 call read_header(iunit,iverbose,debug,doubleprec, &
                    npart,npartoftypei,n1,ntypes,nblocks,narrsizes,dummyreal,tagsreal,nreals,ierr)
 if (ierr /= 0) then
    print "(a)",' *** ERROR READING HEADER ***'
    close(iunit)
    return
 endif
!
!--Attempt to read all MPI blocks
!
 ntotal = 0
 ntotblock = 0
 nptmasstot = 0
 i2 = 0
 iptmass2 = 0
 igotmass = .true.
 imadepmasscolumn = .false.
 massoftypei(:) = 0.
 nkilled = 0

 over_MPIblocks: do iblock=1,nblocks
!
!--read array header from this block
!
    if (iblock==1) ncolstep = 0
    do iarr=1,narrsizes
       call read_block_header(iunit,iblock,iarr,iverbose,debug, &
           isize,nint(iarr),nint1(iarr),nint2(iarr),nint4(iarr),nint8(iarr),&
           nreal(iarr),nreal4(iarr),nreal8(iarr),&
           ntotblock,npart,ntotal,nptmasstot,ncolstep,ierr)
       if (ierr /= 0) then
          print "(a)",' *** ERROR READING ARRAY SIZES ***'
          close(iunit)
          return
       endif
    enddo
    if (debug) print*,'DEBUG: ncolstep=',ncolstep,' from file header, also nptmasstot = ',nptmasstot
!
!--this is a bug fix for a corrupt version of wdump outputting bad
!  small dump files
!
    if (smalldump .and. nreal(1)==5 .and. iblock==1 .and. lenvironment('SSPLASH_FIX_CORRUPT')) then
       print*,'FIXING CORRUPT HEADER ON SMALL DUMPS: assuming nreal=3 not 5'
       nreal(1) = 3
       ncolstep = ncolstep - 2
    endif

    npart_max = maxval(isize(1:narrsizes))
    npart_max = max(npart_max,npart+nptmasstot,ntotal)
!
!--work out from array header how many columns we are going to read
!  in order to allocate memory
!
    if (iblock==1) then
       igotmass = .true.
       if (smalldump .or. phantomdump) then
          if (phantomdump) then
             if (tagged) then
                call extract('massoftype',massoftypei(1:ntypes),dummyreal,tagsreal,nreals,ierr)
             else
                ! old phantom dumps had only 5 types
                call extract('massoftype',massoftypei(1:5),dummyreal,tagsreal,nreals,ierr)
             endif
          else
             call extract('pmassinitial',massoftypei(1),dummyreal,tagsreal,nreals,ierr)
             if (ierr /= 0) then
                print "(a)",' error extracting particle mass from small dump file'
                massoftypei(1) = 0.
                igotmass = .false.
             endif
          endif
          if (debug) print*,'DEBUG: got massoftype(gas) = ',massoftypei(1)
          if (any(massoftypei(1:ntypes) > tiny(0.)) .and. .not.lowmemorymode) then
             ncolstep = ncolstep + 1  ! make an extra column to contain particle mass
             imadepmasscolumn = .true.
          elseif (lowmemorymode) then
             igotmass = .false.
          else
             igotmass = .false.
          endif
          if (all(abs(massoftypei(1:ntypes)) < tiny(0.)) .and. nreal(1) < 4) then
             print "(a)",' error: particle masses not present in small dump file'
             igotmass = .false.
          endif
       endif
       if (debug) print*,'DEBUG: gotmass = ',igotmass, ' ncolstep = ',ncolstep
!
!--   to handle both small and full dumps, we need to place the quantities dumped
!     in both small and full dumps at the start of the dat array
!     quantities only in the full dump then come after
!     also means that hydro/MHD are "semi-compatible" in the sense that x,y,z,m,h
!     and rho are in the same place for both types of dump
!
       ix(1) = 1
       ix(2) = 2
       ix(3) = 3
       if (igotmass) then
          ipmass = 4
          ih = 5
          irho = 6
          nhydroarrays = 6 ! x,y,z,m,h,rho
       else
          ipmass = 0
          ih = 4
          irho = 5
          nhydroarrays = 5 ! x,y,z,h,rho
       endif
       nhydroarraysinfile = nreal(1) + nreal4(1) + nreal8(1)
       nhydroreal4 = nreal4(1)
       if (imadepmasscolumn) nhydroarraysinfile = nhydroarraysinfile + 1
       if (nhydroarraysinfile  < nhydroarrays .and. .not.phantomdump) then
          print "(a)",' ERROR: one of x,y,z,m,h or rho missing in small dump read'
          nhydroarrays = nreal(1)+nreal4(1)+nreal8(1)
       elseif (phantomdump .and. (nreal(1) < 3 .or. nreal4(1) < 1)) then
          print "(a)",' ERROR: x,y,z or h missing in phantom read'
       endif
       ndustarrays = ndusttypes
       if (debug) print*,' DEBUG: ndustarrays = ',ndustarrays
       if (narrsizes >= 4) then
          nmhdarrays = 3 ! Bx,By,Bz
          nmhd = nreal(4) + nreal4(4) + nreal8(4) - nmhdarrays ! how many "extra" mhd arrays
          if (debug) print*,'DEBUG: ',nmhd,' extra MHD arrays'
       else
          nmhdarrays = 0
       endif

       !--radiative transfer dump?
       if (narrsizes >= 3 .and. isize(3)==isize(1)) rtdump = .true.
       !--mhd dump?
       if (narrsizes >= 4) mhddump = .true.

       if (.not.(mhddump.or.smalldump)) then
          ivx = nhydroarrays+ndustarrays+1
       elseif (mhddump .and. .not.smalldump) then
          ivx = nhydroarrays+ndustarrays+nmhdarrays+1
       else
          ivx = 0
       endif
       !--need to force read of velocities e.g. for corotating frame subtraction
       if (any(required(ivx:ivx+ndimV-1))) required(ivx:ivx+ndimV-1) = .true.

       !--force read of h and rho if dustfrac is required
       if (ndustarrays > 0 .and. any(required(nhydroarrays+1:nhydroarrays+ndustarrays))) then
          if (debug) print*,' dustfrac in columns ',nhydroarrays+1,nhydroarrays+ndustarrays,' required = ',required(nhydroarrays+1)
          required(irho) = .true.
          required(ih) = .true.
       endif

       !--for phantom dumps, also make a column for density
       !  and divv, if a .divv file exists
       if (phantomdump) then
          ncolstep = ncolstep + 1
          ! make extra columns in the same place every time
          if (maxcol==0) then
             ncolstepfirst = ncolstep   !  save number of columns
          elseif (ncolstep < ncolstepfirst) then
             ncolstep = ncolstepfirst   !  used saved number of columns
          endif
          inquire(file=trim(dumpfile)//'.divv',exist=iexist)
          if (iexist) then
             idivvxcol   = ncolstep + 1
             icurlvxcol = ncolstep + 2
             icurlvycol = ncolstep + 3
             icurlvzcol = ncolstep + 4
             ncolstep   = ncolstep + 4
          endif
          if (get_temperature) then
            !add a column for the temperature
             ncolstep = ncolstep+1
             itempcol = ncolstep
          endif
          if (get_kappa) then
             ncolstep = ncolstep+1
             ikappa = ncolstep
          endif
          if (get_ionfrac) then
             iHIIcol   = ncolstep + 1
             iHeIIcol  = ncolstep + 2
             iHeIIIcol = ncolstep + 3
             ncolstep  = ncolstep + 3
          endif
          call check_for_composition_file(trim(dumpfile),&
               npart,ncolstep,icomp_col_start,ncomp,tagarr,compfile)
       endif
    endif
!
!--allocate memory now that we know the number of columns
!
    if (iblock==1) then
       ncolumns = ncolstep + ncalc
       if (ncolumns > maxplot) then
          print*,'ERROR with ncolumns = ',ncolumns,' in data read'
          return
       endif
       ilastrequired = 0
       do i=1,ncolumns
          if (required(i)) ilastrequired = i
       enddo
    endif

    need_to_allocate_iphase = (npart_max > maxpart) .or. .not.allocated(iphase)
    if (npart_max > maxpart .or. j > maxstep .or. ncolumns > maxcol) then
       if (lowmemorymode) then
          call alloc(max(npart_max+2,maxpart),j,ilastrequired)
       else
          call alloc(max(npart_max+2,maxpart),j,ncolumns,mixedtypes=.true.)
       endif
    endif
!
!--now that memory has been allocated, copy info from the header into
!  the relevant arrays
!
    if (iblock==1) then
       call extract_variables_from_header(tagsreal,dummyreal,nreals,iverbose,debug, &
           gotbinary,nblocks,nptmasstot,npartoftypei,ntypes,&
           time(j),gamma(j),hfact,npart,ntotal,npartoftype(:,j),masstype(:,j), &
           dat(:,:,j),ix,ih,ipmass,ivx)
       nhdr = min(nreals,maxhdr)
       headervals(1:nhdr,j) = dummyreal(1:nhdr)
       headertags(1:nhdr)   = tagsreal(1:nhdr)

       call set_grain_sizes(nhdr,headertags,headervals(:,j),udist) ! get grain sizes in cm
       call make_tags_unique(nhdr,headertags)
       if (iverbose > 0) call print_dustgrid_info(nhdr,headertags,headervals(:,j),masstype(1,j)*npartoftype(1,j))

       nstepsread = nstepsread + 1
       !
       !--stop reading file here if no columns required
       !
       if (ilastrequired==0) exit over_MPIblocks
    endif
!
!--allocate memory for iphase array now that gotbinary is known
!
    if (need_to_allocate_iphase) call allocate_iphase(iphase,max(npart_max+2,maxpart),phantomdump,gotbinary,npart_max+1)
!
!--Arrays
!
    imaxcolumnread = 0
    icolumn = 0
    idustarr   = 0
    istartmhd = 0
    istartrt = 0
    i1 = i2 + 1
    i2 = i1 + isize(1) - 1
    if (debug) then
       print "(1x,a10,i4,3(a,i12))",'MPI block ',iblock,':  particles: ',i1,' to ',i2,' of ',npart
    elseif (nblocks > 1) then
       if (iblock==1) write(*,"(a,i1,a)",ADVANCE="no") ' reading MPI blocks: .'
       write(*,"('.')",ADVANCE="no")
    endif
    iptmass1 = iptmass2 + 1
    iptmass2 = iptmass1 + isize(2) - 1
    nptmass = nptmasstot
    if (nptmass > 0 .and. debug) print "(15x,3(a,i12))",'  pt. masses: ',iptmass1,' to ',iptmass2,' of ',nptmass

    do iarr=1,narrsizes
       if (nreal(iarr) + nreal4(iarr) + nreal8(iarr) > 0) then
          if (iarr==4) then
             istartmhd = imaxcolumnread + 1
             if (debug) print*,' istartmhd = ',istartmhd
          elseif (iarr==3 .and. rtdump) then
             istartrt = max(nhydroarrays+nmhdarrays+1,imaxcolumnread + 1)
             if (debug) print*,' istartrt = ',istartrt
          endif
       endif
!--read iphase from array block 1
       if (iarr==1) then
          !--skip default int
          nskip = nint(iarr)
          do i=1,nskip
             if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
             read(iunit,end=33,iostat=ierr)
          enddo
          if (nint1(iarr) < 1) then
             if (.not.phantomdump .or. any(npartoftypei(2:) > 0)) then
                if (iverbose > 0) print "(a)",' WARNING: can''t locate iphase in dump'
             elseif (phantomdump) then
                if (iverbose > 0) print "(a)",' WARNING: can''t locate iphase in dump'
             endif
             gotiphase = .false.
             !--skip remaining integer arrays
             nskip = nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
          else
             gotiphase = .true.
             if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
             read(iunit,end=33,iostat=ierr) iphase(i1:i2)
             !--skip remaining integer arrays
             nskip = nint1(iarr) - 1 + nint2(iarr) + nint4(iarr) + nint8(iarr)
          endif
       elseif (smalldump .and. iarr==2 .and. isize(iarr) > 0 .and. .not.phantomdump) then
!--read listpm from array block 2 for small dumps (needed here to extract sink masses)
          if (allocated(listpm)) deallocate(listpm)
          allocate(listpm(isize(iarr)))
          if (nint(iarr) < 1) then
             print "(a)",'ERROR: can''t locate listpm in dump'
             nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
          else
             if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
             read(iunit,end=33,iostat=ierr) listpm(1:isize(iarr))
             nskip = nint(iarr) - 1 + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
          endif
       else
!--otherwise skip all integer arrays (not needed for plotting)
          nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
       endif

       if (iarr==3 .and. lenvironment('SSPLASH_BEN_HACKED')) then
          nskip = nskip - 1
          print*,' FIXING HACKED DUMP FILE'
       endif
       !print*,'skipping ',nskip
       do i=1,nskip
          if (tagged) read(iunit,end=33,iostat=ierr) tagtmp! skip tags
          !
          ! read the Adaptive Particle Refinement level array
          ! in order to correctly set particle masses, if present
          !
          select case(tagtmp)
          case('apr')
             allocate(level(isize(iarr)),massfac(isize(iarr)))
             read(iunit,end=33,iostat=ierr) level
             ! m = m/2**(level-1)
             massfac(1:isize(iarr)) = 1./2**(level(1:isize(iarr))-1)
             if (iblock==1) print "(a,i2)",' :: '//trim(tagtmp)//' max level = ',maxval(level)
             deallocate(level)
          case default
             read(iunit,end=33,iostat=ierr)
          end select
       enddo
!
!--real arrays
!
       if (iarr==2) then
!--read sink particles from phantom dumps
          if (phantomdump .and. iarr==2 .and. isize(iarr) > 0) then
             if (nreal(iarr) < 5) then
                print "(a)",'ERROR: not enough arrays written for sink particles in phantom dump'
                nskip = nreal(iarr)
             else
                if (debug) print*,'DEBUG: denoting ',npart,'->',npart+isize(iarr),' as sink particles'
                iphase(npart+1:npart+isize(iarr)) = -3
                ilocvx = nreal(iarr)-2 ! velocity is always last 3 numbers for phantom sinks
                if (doubleprec) then
                   !--convert default real to single precision where necessary
                   if (debug) print*,'DEBUG: reading sink data, converting from double precision ',isize(iarr)
                   if (allocated(dattemp)) deallocate(dattemp)
                   allocate(dattemp(isize(iarr)),stat=ierr)
                   if (ierr /= 0) then
                      print "(a)",'ERROR in memory allocation'
                      return
                   endif
                   tagtmp = ''
                   do k=1,nreal(iarr)
                      if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
                      if (debug) print*,'DEBUG: reading sink array ',k,isize(iarr),' tag = ',trim(tagtmp)
                      read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                      if (ierr /= 0) print*,' ERROR during read of sink particle data, array ',k

                      iloc = map_sink_property_to_column(k,ilocvx,size(dat(1,:,j)))
                      if (iloc > 0) then
                         do i=1,isize(iarr)
                            dat(npart+i,iloc,j) = real(dattemp(i))
                         enddo
                      elseif (trim(tagtmp)=='hsoft' .and. ih > 0) then
                         do i=1,isize(iarr)
                            if (abs(dat(npart+i,ih,j)) < tiny(0.)) then
                               dat(npart+i,ih,j) = real(dattemp(i))
                               if (i == 1) print*,'zero accretion radius: taking sink particle radius from softening length'
                            endif
                         enddo
                      else
                         if (debug) print*,'DEBUG: skipping sink particle array ',k
                      endif
                   enddo
                else
                   if (debug) print*,'DEBUG: reading sink data, converting from single precision ',isize(iarr)
                   if (allocated(dattempsingle)) deallocate(dattempsingle)
                   allocate(dattempsingle(isize(iarr)),stat=ierr)
                   if (ierr /= 0) then
                      print "(a)",'ERROR in memory allocation'
                      return
                   endif
                   do k=1,nreal(iarr)
                      if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
                      if (debug) print*,'DEBUG: reading sink array ',k,isize(iarr),' tag = ',trim(tagtmp)
                      read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                      if (ierr /= 0) print*,' ERROR during read of sink particle data, array ',k

                      iloc = map_sink_property_to_column(k,ilocvx,size(dat(1,:,j)))
                      if (iloc > 0) then
                         do i=1,isize(iarr)
                            dat(npart+i,iloc,j) = real(dattempsingle(i))
                         enddo
                      elseif (trim(tagtmp)=='hsoft' .and. ih > 0) then
                         do i=1,isize(iarr)
                            if (abs(dat(npart+i,ih,j)) < tiny(0.)) then
                               dat(npart+i,ih,j) = real(dattempsingle(i))
                               if (i == 1) print*,'zero accretion radius: taking sink particle radius from softening length'
                            endif
                         enddo
                      else
                         if (debug) print*,'DEBUG: skipping sink particle array ',k
                      endif
                   enddo
                endif
                ! PSEUDO-remove accreted sinks
                call set_sink_merged(npart+1,int(npart+isize(iarr)),ih,ipmass,dat(:,:,j))
                ! DEFINE density on sink particles (needed for opacity rendering)
                if (required(irho)) call set_sink_density(npart+1,int(npart+isize(iarr)),ih,ipmass,irho,dat(:,:,j))
                npart  = npart + isize(iarr)
             endif
          elseif (smalldump .and. iarr==2 .and. allocated(listpm)) then
!--for sphNG, read sink particle masses from block 2 for small dumps
             if (nreal(iarr) < 1) then
                if (isize(iarr) > 0) print "(a)",'ERROR: sink masses not present in small dump'
                nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
             else
                if (doubleprec) then
                   !--convert default real to single precision where necessary
                   if (allocated(dattemp)) deallocate(dattemp)
                   allocate(dattemp(isize(iarr)),stat=ierr)
                   if (ierr /=0) print "(a)",'ERROR in memory allocation'
                   if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
                   read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                   if (nptmass /= isize(iarr)) print "(a)",'ERROR: nptmass /= block size'
                   if (ipmass > 0) then
                      do i=1,isize(iarr)
                         dat(listpm(iptmass1+i-1),ipmass,j) = real(dattemp(i))
                      enddo
                   else
                      print*,'WARNING: sink particle masses not read because no mass array allocated'
                   endif
                else
                   !--convert default real to double precision where necessary
                   if (allocated(dattempsingle)) deallocate(dattempsingle)
                   allocate(dattempsingle(isize(iarr)),stat=ierr)
                   if (ierr /=0) print "(a)",'ERROR in memory allocation'
                   if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
                   read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                   if (nptmass /= isize(iarr)) print "(a)",'ERROR: nptmass /= block size'
                   if (ipmass > 0) then
                      do i=1,isize(iarr)
                         dat(listpm(iptmass1+i-1),ipmass,j) = real(dattempsingle(i))
                      enddo
                   else
                      print*,'WARNING: sink particle masses not read because no mass array allocated'
                   endif
                endif
                nskip = nreal(iarr) - 1 + nreal4(iarr) + nreal8(iarr)
             endif
          else
!--for other blocks, skip real arrays if size different
             nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
          endif
          do i=1,nskip
             if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
             read(iunit,end=33,iostat=ierr)
          enddo
          ! deallocate dattempsingle
          if (allocated(dattempsingle)) deallocate(dattempsingle)

       elseif (isize(iarr)==isize(1)) then
!
!--read all real arrays defined on all the particles (same size arrays as block 1)
!
          if ((doubleprec.and.nreal(iarr) > 0).or.nreal8(iarr) > 0) then
             if (allocated(dattemp)) deallocate(dattemp)
             allocate(dattemp(isize(iarr)),stat=ierr)
             if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattemp)'
          elseif (nreal(iarr) > 0 .or. nreal8(iarr) > 0) then
             if (allocated(dattempsingle)) deallocate(dattempsingle)
             allocate(dattempsingle(isize(iarr)),stat=ierr)
             if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattempsingle)'
          endif
!        default reals may need converting
          do i=1,nreal(iarr)
             tagtmp = ''
             if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
             icolumn = assign_column(tagtmp,iarr,i,6,imaxcolumnread,idustarr,ncolstep)
             if (tagged) tagarr(icolumn) = tagtmp
             ! force data read if label matches one of the required columns
             if (.not.required(icolumn)) then
                if (match_tag(labelreq(1:nreq),tagtmp) > 0) then
                   print "(1x,a,i2)",'-> found '//trim(tagtmp)//' in column ',icolumn
                   required(icolumn) = .true.
                endif
             endif

             if (debug)  print*,' reading real to col:',icolumn,' tag = ',trim(tagtmp), ' required = ',required(icolumn)
             if (required(icolumn)) then
                if (doubleprec) then
                   read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                   dat(i1:i2,icolumn,j) = real(dattemp(1:isize(iarr)))
                else
                   read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                   dat(i1:i2,icolumn,j) = real(dattempsingle(1:isize(iarr)))
                endif
             else
                read(iunit,end=33,iostat=ierr)
             endif
          enddo
!
!        set masses for equal mass particles (not dumped in small dump or in phantom)
!
          if (((smalldump.and.nreal(1) < ipmass).or.phantomdump).and. iarr==1) then
             if (any(abs(masstype(:,j)) > tiny(masstype))) then
                icolumn = ipmass
                if (required(ipmass) .and. ipmass > 0) then
                   if (phantomdump) then
                      dat(i1:i2,ipmass,j) = masstype(itypemap_phantom(iphase(i1:i2)),j)
                      if (allocated(massfac)) then
                         dat(i1:i2,ipmass,j) = dat(i1:i2,ipmass,j)*massfac(1:isize(iarr))
                         deallocate(massfac)
                      endif
                   else
                      where (iphase(i1:i2)==0) dat(i1:i2,icolumn,j) = masstype(1,j)
                   endif
                endif
                !--dust mass for phantom particles
                if (phantomdump .and. npartoftypei(itypemap_dust_phantom) > 0 .and. ipmass > 0) then
                   print*,'dust particle mass = ',masstype(itypemap_dust_phantom,j),&
                         ' ratio m_dust/m_gas = ',masstype(itypemap_dust_phantom,j)/masstype(1,j)
                endif
                if (debug) print*,'mass ',icolumn
             elseif (phantomdump .and. npartoftypei(1) > 0) then
                print*,' ERROR: particle mass zero in Phantom dump file!'
             endif
          endif
!
!        real4 arrays (may need converting if splash is compiled in double precision)
!
          if (nreal4(iarr) > 0 .and. kind(dat)==doub_prec) then
             if (allocated(dattempsingle)) deallocate(dattempsingle)
             allocate(dattempsingle(isize(iarr)),stat=ierr)
             if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattempsingle)'
          endif

          if (debug) print*,'DEBUG: SIZE of dattempsingle',size(dattempsingle)
!        real4s may need converting
          imaxcolumnread = max(imaxcolumnread,icolumn)
          if ((nreal(iarr)+nreal4(iarr)) > 6) imaxcolumnread = max(imaxcolumnread,6)

          do i=1,nreal4(iarr)
             tagtmp = ''
             if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
             icolumn = assign_column(tagtmp,iarr,i,4,imaxcolumnread,idustarr,ncolstep)
             if (debug) print*,'reading real4 to col:',icolumn,' tag = ',trim(tagtmp),' required = ',required(icolumn)
             if (tagged) tagarr(icolumn) = tagtmp

             if (phantomdump .and. icolumn==ih) required(ih) = .true. ! h always required for density

             if (required(icolumn)) then
                if (allocated(dattempsingle)) then
                   read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                   dat(i1:i2,icolumn,j) = real(dattempsingle(1:isize(iarr)))
                else
                   read(iunit,end=33,iostat=ierr) dat(i1:i2,icolumn,j)
                endif
             else
                read(iunit,end=33,iostat=ierr)
             endif
             !--construct density for phantom dumps based on h, hfact and particle mass
             if (phantomdump .and. icolumn==ih) then
                icolumn = irho ! density
                call get_rho_from_h(i1,i2,ih,ipmass,irho,required,npartoftype(:,j),&
                                    masstype(:,j),hfact,dat(:,:,j),iphase,nkilled)
             endif
          enddo
!        real 8's need converting
          do i=1,nreal8(iarr)
             tagtmp = ''
             if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
             icolumn = assign_column(tagtmp,iarr,i,8,imaxcolumnread,idustarr,ncolstep)
             if (debug) print*,'reading real8 to col:',icolumn,' tag = ',trim(tagtmp),' required = ',required(icolumn)
             if (tagged) tagarr(icolumn) = tagtmp
             if (required(icolumn)) then
                read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                dat(i1:i2,icolumn,j) = real(dattemp(1:isize(iarr)))
             else
                read(iunit,end=33,iostat=ierr)
             endif
          enddo
       endif
    enddo ! over array sizes
 enddo over_MPIblocks
!
!--reached end of file (during data read)
!
 goto 34
33 continue
 print "(/,1x,a,/)",'*** WARNING: END OF FILE DURING READ ***'
 print*,'Press any key to continue (but there is likely something wrong with the file...)'
 read*
34 continue

 call set_labels_sphNG
 !
 !--emit warning if we find particles with h = 0
 !
 if (nkilled > 0) print*,'WARNING: got ',nkilled,' dead (not accreted) particles, should not happen'
 !
 !--read .divv file for phantom dumps
 !
 if (phantomdump .and. idivvxcol /= 0 .and. any(required(idivvxcol:icurlvzcol))) then
    print "(a)",' reading divv from '//trim(dumpfile)//'.divv'
    open(unit=66,file=trim(dumpfile)//'.divv',form='unformatted',status='old',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(dumpfile)//'.divv'
    else
       read(66,iostat=ierr) dat(1:ntotal,idivvxcol,j)
       if (ierr /= 0) print "(a)",' WARNING: ERRORS reading divv from file'
       if (any(required(icurlvxcol:icurlvzcol))) then
          read(66,iostat=ierr) dat(1:ntotal,icurlvxcol,j)
          read(66,iostat=ierr) dat(1:ntotal,icurlvycol,j)
          read(66,iostat=ierr) dat(1:ntotal,icurlvzcol,j)
       endif
       if (ierr /= 0) print "(a)",' WARNING: ERRORS reading curlv from file'
       close(66)
    endif
 endif

 if (icomp_col_start > 0 .and. any(required(icomp_col_start:icomp_col_start+ncomp))) then
    call read_kepler_composition(compfile,ntotal,dat(:,:,j),icomp_col_start,ncomp)
 endif
 !
 !--calculate the temperature from density and internal energy (using physical units)
 !
 unit_dens = umass/(udist**3)
 !
 !--use primitive density for relativistic code
 !
 idenscol = irho
 if (irhorestframe > 0) idenscol = irhorestframe

 if (get_temperature .and. itempcol > 0 .and. required(itempcol)) then
    unit_ergg = (udist/utime)**2
    dat(1:ntotal,itempcol,j) = get_temp_from_u(dat(1:ntotal,idenscol,j)*unit_dens,dat(1:ntotal,iutherm,j)*unit_ergg) !irho = density
 endif
 if (get_kappa .and. ikappa > 0 .and. required(ikappa) .and. itemp > 0) then
    print*,'X,Y,Z = ',Xfrac,Yfrac,1.-Xfrac-Yfrac
    dat(1:ntotal,ikappa,j) = get_opacity(dat(1:ntotal,idenscol,j)*unit_dens,dat(1:ntotal,itemp,j)*1.d0,Xfrac,Yfrac)
 endif
 if (get_ionfrac .and. iHIIcol > 0 .and. iHeIIcol > 0 .and. iHeIIIcol > 0&
     .and. any(required(iHIIcol:iHeIIIcol))) then
    do i=1,ntotal
       call ionisation_fraction(real(dat(i,idenscol,j)*unit_dens),dat(i,itemp,j),&
                                real(Xfrac),real(Yfrac),xHIi,xHIIi,xHeIi,xHeIIi,xHeIIIi,nei)
       dat(i,iHIIcol,j)=xHIIi
       dat(i,iHeIIcol,j)=xHeIIi
       dat(i,iHeIIIcol,j)=xHeIIIi
    enddo
 endif

 !
 !--reset centre of mass to zero if environment variable "SSPLASH_RESET_CM" is set,
 !  or reset centre to densest clump if environment variable "SSPLASH_RESET_DENSE" is set
 !  the latter will override the former
 ! (updated from n1 to npart since order is not preserved when dumping data; JHW)
 icentre = 0
 if (lenvironment('SSPLASH_RESET_CM'))    icentre = 1
 if (lenvironment('SSPLASH_RESET_DENSE')) icentre = 2
 if (allocated(dat) .and. npart > 0 .and. npart <= size(dat(:,1,1)) .and. icentre > 0 .and. allocated(iphase)) then
    call reset_centre_of_mass(dat(1:npart,1:3,j),dat(1:npart,4,j),dat(1:npart,5,j),iphase(1:npart),npart,icentre)
 endif
 !
 !--reset corotating frame velocities if environment variable "SSPLASH_OMEGA" is set
 !
 if (allocated(dat) .and. n1 > 0 .and. all(required(1:2))) then
    omega = renvironment('SSPLASH_OMEGAT')
    if (abs(omega) > tiny(omega) .and. ndim >= 2) then
       call reset_corotating_positions(n1,dat(1:n1,1:2,j),omega,time(j))
    endif

    if (.not. smalldump) then
       if (abs(omega) < tiny(omega)) omega = renvironment('SSPLASH_OMEGA')
       if (abs(omega) > tiny(omega) .and. ivx > 0) then
          if (.not.all(required(1:2)) .or. .not.all(required(ivx:ivx+1))) then
             print*,' ERROR subtracting corotating frame with partial data read'
          else
             call reset_corotating_velocities(n1,dat(1:n1,1:2,j),dat(1:n1,ivx:ivx+1,j),omega)
          endif
       endif
    endif
 endif

 !--set flag to indicate that only part of this file has been read
 if (.not.all(required(1:ncolstep))) ipartialread = .true.

 nptmassi = 0
 nunknown = 0
 ngas  = 0
 ndust = 0
 nstar = 0
 !--can only do this loop if we have read the iphase array
 iphasealloc: if (allocated(iphase)) then
!
!--sanity check the iphase array
!
    if (gotiphase) call check_iphase_matches_npartoftype(1,npart,iphase,npartoftype(:,j))
!
!--translate iphase into particle types (mixed type storage)
!
    if (size(iamtype(:,j)) > 1) then
       if (phantomdump) then
          !
          !--phantom: translate iphase to splash types
          !
          do i=1,npart
             itype = itypemap_phantom(iphase(i))
             iamtype(i,j) = itype
             if (ih > 0 .and. required(ih)) then
                if (dat(i,ih,j) <= 0. .and. itype /= itypemap_sink_phantom) iamtype(i,j) = itypemap_unknown_phantom
             endif
          enddo
       else
          !
          !--sphNG: translate iphase to splash types
          !
          do i=1,npart
             itype = itypemap_sphNG(iphase(i))
             iamtype(i,j) = itype
             select case(itype)
             case(1)
                ngas = ngas + 1
             case(2)
                ndust = ndust + 1
             case(4)
                nptmassi = nptmassi + 1
             case(5)
                nstar = nstar + 1
             case default
                nunknown = nunknown + 1
             end select
          enddo
          do i=npart+1,ntotal
             iamtype(i,j) = 2
          enddo
       endif
       !print*,'mixed types: ngas = ',ngas,nptmassi,nunknown

    elseif (any(iphase(1:ntotal) /= 0)) then
       if (phantomdump) then
          print*,'ERROR: low memory mode will not work correctly with phantom + multiple types'
          print*,'press any key to ignore this and continue anyway (at your own risk...)'
          read*
       endif
!
!--place point masses after normal particles
!  if not storing the iamtype array
!
       print "(a)",' sorting particles by type...'
       nunknown = 0
       do i=1,npart
          if (iphase(i) /= 0) nunknown = nunknown + 1
       enddo
       ncolcopy = min(ncolstep,maxcol)
       allocate(dattemp2(nunknown,ncolcopy))

       iphaseminthistype = 0  ! to avoid compiler warnings
       iphasemaxthistype = 0
       do itype=1,3
          nthistype = 0
          ipos = 0
          select case(itype)
          case(1) ! ptmass
             iphaseminthistype = 1
             iphasemaxthistype = 9
          case(2) ! star
             iphaseminthistype = 10
             iphasemaxthistype = huge(iphasemaxthistype)
          case(3) ! unknown
             iphaseminthistype = -huge(iphaseminthistype)
             iphasemaxthistype = -1
          end select

          do i=1,ntotal
             ipos = ipos + 1
             if (iphase(i) >= iphaseminthistype .and. iphase(i) <= iphasemaxthistype) then
                nthistype = nthistype + 1
                !--save point mass information in temporary array
                if (nptmassi > size(dattemp2(:,1))) stop 'error: ptmass array bounds exceeded in data read'
                dattemp2(nthistype,1:ncolcopy) = dat(i,1:ncolcopy,j)
                !             print*,i,' removed', dat(i,1:3,j)
                ipos = ipos - 1
             endif
             !--shuffle dat array
             if (ipos /= i .and. i < ntotal) then
                !           print*,'copying ',i+1,'->',ipos+1
                dat(ipos+1,1:ncolcopy,j) = dat(i+1,1:ncolcopy,j)
                !--must also shuffle iphase (to be correct for other types)
                iphase(ipos+1) = iphase(i+1)
             endif
          enddo

          !--append this type to end of dat array
          do i=1,nthistype
             ipos = ipos + 1
             !          print*,ipos,' appended', dattemp2(i,1:3)
             dat(ipos,1:ncolcopy,j) = dattemp2(i,1:ncolcopy)
             !--we make iphase = 1 for point masses (could save iphase and copy across but no reason to)
             iphase(ipos) = iphaseminthistype
          enddo

          select case(itype)
          case(1)
             nptmassi = nthistype
             if (nptmassi /= nptmass) print *,'WARNING: nptmass from iphase =',nptmassi,'not equal to nptmass =',nptmass
          case(2)
             nstar = nthistype
          case(3)
             nunknown = nthistype
          end select
       enddo

    endif

 endif iphasealloc

 if (allocated(dattemp)) deallocate(dattemp)
 if (allocated(dattempsingle)) deallocate(dattempsingle)
 if (allocated(dattemp2)) deallocate(dattemp2)
 if (allocated(iphase)) deallocate(iphase)
 if (allocated(listpm)) deallocate(listpm)

 call set_labels_sphNG
 if (.not.phantomdump) then
    if (ngas /= npart - nptmassi - ndust - nstar - nunknown) &
           print*,'WARNING!!! ngas =',ngas,'but should be',npart-nptmassi-ndust-nstar-nunknown
    if (ndust /= npart - nptmassi - ngas - nstar - nunknown) &
           print*,'WARNING!!! ndust =',ndust,'but should be',npart-nptmassi-ngas-nstar-nunknown
    npartoftype(:,j) = 0
    npartoftype(1,j) = ngas  ! npart - nptmassi - ndust - nstar - nunknown
    npartoftype(2,j) = ndust ! npart - nptmassi - ngas  - nstar - nunknown
    npartoftype(3,j) = ntotal - npart
    npartoftype(4,j) = nptmassi
    npartoftype(5,j) = nstar
    npartoftype(6,j) = nunknown
 else
    if (debug) print*,' DEBUG: nunknown = ',nunknown
    npartoftype(1,j) = npartoftype(1,j) - nunknown
    npartoftype(itypemap_unknown_phantom,j) = npartoftype(itypemap_unknown_phantom,j) + nunknown
 endif

 if (iverbose > 0) call print_types(npartoftype(:,j),labeltype)

 close(15)
 if (debug) print*,' finished data read, npart = ',npart, ntotal, npartoftype(1:ntypes,j)

 return

contains

!
!--reset centre of mass to zero
!
subroutine reset_centre_of_mass(xyz,pmass,h,iphase,np,icentre)
 implicit none
 integer, intent(in) :: np,icentre
 real, dimension(np,3), intent(inout) :: xyz
 real, dimension(np), intent(in) :: h,pmass
 integer(kind=int1), dimension(np), intent(in) :: iphase
 real :: masstot,pmassi,minh
 real, dimension(3) :: xcm
 integer :: i,ctr
 !
 !--get centre of mass
 !
 xcm(:)  = 0.
 masstot = 0.
 minh    = huge(minh)
 do i=1,np
    if (iphase(i) >= 0) then
       pmassi  = pmass(i)
       masstot = masstot + pmass(i)
       where (required(1:3)) xcm(:) = xcm(:) + pmassi*xyz(i,:)
       minh = min(h(i),minh)
    endif
 enddo
 !
 !--if requested, find the location of the densest clump
 if (icentre==2) then
    xcm(:)  = 0.
    masstot = 0.
    ctr     = 0
    do i=1,np
       if (iphase(i) >= 0 .and. h(i) < minh*1.05) then
          ctr     = ctr + 1
          pmassi  = pmass(i)
          masstot = masstot + pmass(i)
          where (required(1:3)) xcm(:) = xcm(:) + pmassi*xyz(i,:)
       endif
    enddo
 endif

 xcm(:) = xcm(:)/masstot
 if (icentre==1) then
    print*,' RESETTING CENTRE OF MASS (',pack(xcm,required(1:3)),') TO ZERO '
 elseif (icentre==2) then
    print*,' RESETTING CENTRE OF DENSEST CLUMP (',pack(xcm,required(1:3)),') TO ZERO using ',ctr,' particles'
 endif
 if (required(1)) xyz(1:np,1) = xyz(1:np,1) - xcm(1)
 if (required(2)) xyz(1:np,2) = xyz(1:np,2) - xcm(2)
 if (required(3)) xyz(1:np,3) = xyz(1:np,3) - xcm(3)

 return
end subroutine reset_centre_of_mass

subroutine reset_corotating_velocities(np,xy,velxy,omeg)
 integer, intent(in) :: np
 real, dimension(np,2), intent(in) :: xy
 real, dimension(np,2), intent(inout) :: velxy
 real, intent(in) :: omeg
 integer :: ip

 print*,'SUBTRACTING COROTATING VELOCITIES, OMEGA = ',omeg
 do ip=1,np
    velxy(ip,1) = velxy(ip,1) + xy(ip,2)*omeg
 enddo
 do ip=1,np
    velxy(ip,2) = velxy(ip,2) - xy(ip,1)*omeg
 enddo

 return
end subroutine reset_corotating_velocities

subroutine reset_corotating_positions(np,xy,omeg,t)
 integer, intent(in) :: np
 real, dimension(np,2), intent(inout) :: xy
 real, intent(in) :: omeg,t
 real :: phii,phinew,r
 integer :: ip

 print*,'SUBTRACTING COROTATING POSITIONS, OMEGA = ',omeg,' t = ',t
!$omp parallel default(none) &
!$omp shared(xy,np) &
!$omp firstprivate(omeg,t) &
!$omp private(ip,r,phii,phinew)
!$omp do
 do ip=1,np
    r = sqrt(xy(ip,1)**2 + xy(ip,2)**2)
    phii = atan2(xy(ip,2),xy(ip,1))
    phinew = phii + omeg*t
    xy(ip,1) = r*cos(phinew)
    xy(ip,2) = r*sin(phinew)
 enddo
!$omp enddo
!$omp end parallel

 return
end subroutine reset_corotating_positions

end subroutine read_data_sphNG

!------------------------------------------------------------
! set labels for each column of data
!------------------------------------------------------------
subroutine set_labels_sphNG
 use labels, only:label,unitslabel=>unitslabel_default,&
              labelzintegration=>labelzintegration_default,labeltype,labelvec,iamvec, &
              ix,ipmass,irho,ih,iutherm,ipr,ivx,iBfirst,idivB,iJfirst,icv,iradenergy,&
              idustfrac,ideltav,idustfracsum,ideltavsum,igrainsize,igraindens, &
              ivrel,make_vector_label,get_label_grain_size,itemp,ikappa,ipmomx,irhorestframe
 use params
 use settings_data,   only:ndim,ndimV,ntypes,ncolumns,UseTypeInRenderings,debugmode
 use geometry,        only:labelcoord
 use settings_units,  only:units=>units_default,unitzintegration=>unitzintegration_default,&
                           get_nearest_length_unit,get_nearest_time_unit,&
                           get_nearest_mass_unit,get_nearest_velocity_unit
 use sphNGread
 use asciiutils,      only:lcase,make_tags_unique,match_tag
 use system_utils,    only:lenvironment,get_environment_or_flag
 integer :: i,j,idustlast
 real(doub_prec)   :: unitx,unitvel,unitmass
 character(len=20) :: string,unitlabelx,unitlabelv
 character(len=20) :: deltav_string

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_sphNG ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_sphNG ***'
    return
 endif
!--all formats read the following columns
 do i=1,ndim
    ix(i) = i
 enddo
 if (igotmass) then
    ipmass = 4   !  particle mass
    ih = 5       !  smoothing length
 else
    ipmass = 0
    ih = 4       !  smoothing length
 endif
 irho = ih + 1     !  density
 iutherm = 0
 idustfrac = 0

 !
 !--translate array tags into column labels, where necessary
 !
 if (tagged) then
    do i=1,ncolumns
       label(i) = tagarr(i)
       select case(trim(tagarr(i)))
       case('m')
          ipmass = i
       case('h')
          ih = i
       case('rho')
          irho = i
       case('vx')
          ivx = i
       case('pressure')
          ipr = i
       case('u')
          iutherm = i
       case('divv')
          idivvcol = i
          units(i) = 1./utime
          unitslabel(i) = ' [1/s]'
       case('divvx')
          idivvxcol = i
          units(i) = 1./utime
          unitslabel(i) = ' [1/s]'
       case('poten')
          units(i) = umass*(udist/utime)**2
          unitslabel(i) = ' [erg]'
       case('dt')
          units(i) = utime
          unitslabel(i) = ' [s]'
       case('curlvx')
          icurlvxcol = i
       case('curlvy')
          icurlvycol = i
       case('curlvz')
          icurlvzcol = i
       case('Bx')
          iBfirst = i
       case('divB')
          idivB = i
       case('curlBx')
          iJfirst = i
       case('psi')
          label(i) = '\psi'
          units(i) = umagfd*udist/utime
          unitslabel(i) = ' [G cm/s]'
       case('psi/c_h')
          units(i) = umagfd
          unitslabel(i) = ' [G]'
       case('dustfracsum')
          idustfracsum = i
       case('deltavsumx')
          ideltavsum = i
       case('deltavx')
          ideltav = i
       case('alpha')
          label(i) = '\alpha'
       case('alphaB')
          label(i) = '\alpha_B'
       case('EulerAlpha')
          label(i) = 'Euler \alpha'
       case('EulerBeta')
          label(i) = 'Euler \beta'
       case('EtaReal')
          label(i) = '\eta_{real}'
       case('EtaArtificial')
          label(i) = '\eta_{art}'
       case('Erad')
          iradenergy = i
          label(i) = 'radiation energy'
          unitslabel(i) = ' [erg/g]'
          units(i) = (udist/utime)**2
       case('opacity')
          label(i) = 'opacity'
          unitslabel(i) = ' [cm^2/g]'
          units(i) = udist**2/umass
       case('EddingtonFactor')
          label(i) = 'Eddington Factor'
       case('Cv')
          label(i) = 'u/T'
          icv = i
          units(i) = (udist/utime)**2
          unitslabel(i) = ' [erg/(g K)]'
       case('h2ratio')
          label(i) = 'H_2 ratio'
       case('abH1q','abHIq')
          label(i) = 'HI abundance'
       case('abhpq')
          label(i) = 'proton abundance'
       case('abeq')
          label(i) = 'e^- abundance'
       case('abco')
          label(i) = 'CO abundance'
       case('eta_{OR}')
          units(i:i+2) = udist**2/utime
          unitslabel(i:i+2) = ' [cm^2/s]'
       case('grainsize')
          igrainsize = i
       case('graindens')
          igraindens = i
       case('temperature')
          itemp = i
       case('vrel')
          ivrel = i
       case('px')
          ipmomx = i
       case('dens prim')
          irhorestframe = i
       case default
          if (debugmode) print "(a,i2)",' DEBUG: Unknown label '''//trim(tagarr(i))//''' in column ',i
          label(i) = tagarr(i)
       end select
    enddo
 else
    if (smalldump .and. nhydroreal4 >= 3) iutherm = irho+1
    call guess_labels(ncolumns,iamvec,label,labelvec,istartmhd,istart_extra_real4,&
          nmhd,nhydroreal4,ndimV,irho,iBfirst,ivx,iutherm,idivB,iJfirst,&
          iradenergy,icv,udist,utime,units,unitslabel)
 endif

 if (itemp==0 .or. itempcol > 0) itemp=itempcol ! if temperature not found in file, use computed one

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 if (irho > 0) label(irho) = 'density'
 if (iutherm > 0) label(iutherm) = 'u'
 if (ih > 0) label(ih) = 'h       '
 if (ipmass > 0) label(ipmass) = 'particle mass'
 if (idivB > 0) label(idivB) = 'div B'
 if (idivvxcol > 0) label(idivvxcol) = 'div vx'
 if (idivvcol > 0) label(idivvcol) = 'div v'
 if (itemp > 0) then
    label(itemp) = 'temperature'
    unitslabel(itemp) = ' [K]'
 endif
 if (itempcol > 0) then
    if (itempcol /= itemp) label(itempcol) = 'temperature (from u)'
    unitslabel(itempcol) = ' [K]'
 endif
 if (ikappa > 0) label(ikappa) = 'kappa'
 if (iHIIcol > 0) label(iHIIcol) = 'HII fraction'
 if (iHeIIcol > 0) label(iHeIIcol) = 'HeII fraction'
 if (iHeIIIcol > 0) label(iHeIIIcol) = 'HeIII fraction'
 if (icurlvxcol > 0 .and. icurlvycol > 0 .and. icurlvzcol > 0) then
    call make_vector_label('curl v',icurlvxcol,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 endif
 if (ideltavsum > 0) then
    ! Modify the deltavsum labels to have vector subscripts
    do j=1,ndimV
       label(ideltavsum+j-1) = 'deltavsum'//'_'//labelcoord(j,1)
    enddo
    ! Make N deltav labels with vector subscripts
    do i = ideltavsum+ndimV,ideltav,ndimV
       write(deltav_string,'(I10)') (i-ideltavsum)/ndimV
       write(deltav_string,'(A)') 'deltav'//trim(adjustl(deltav_string))
       do j=1,ndimV
          label(i+j-1) = trim(deltav_string)//'_'//labelcoord(j,1)
       enddo
    enddo
 endif
 !
 !--set labels for vector quantities
 !
 call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('B',iBfirst,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('J',iJfirst,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('p',ipmomx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 !
 !--ensure labels are unique by appending numbers where necessary
 !
 call make_tags_unique(ncolumns,label)
 !
 !--identify dust fraction in the case where there is only one species
 !
 idustfrac = match_tag(label,'dustfrac')
 !
 !--set units for plot data
 !
 call get_nearest_length_unit(udist,unitx,unitlabelx)
 call get_nearest_velocity_unit(udist/utime,unitvel,unitlabelv)
 if (ndim >= 3) then
    units(1:3) = unitx
    unitslabel(1:3) = unitlabelx
 endif
 if (ipmass > 0) then
    call get_nearest_mass_unit(umass,unitmass,unitslabel(ipmass))
    units(ipmass) = unitmass
 endif
 units(ih) = unitx
 unitslabel(ih) = unitlabelx
 if (ivx > 0) then
    units(ivx:ivx+ndimV-1) = unitvel
    unitslabel(ivx:ivx+ndimV-1) = unitlabelv
 endif
 if (ipmomx > 0) then
    units(ipmomx:ipmomx+ndimV-1) = unitvel
    unitslabel(ipmomx:ipmomx+ndimV-1) = unitlabelv
 endif
 if (ideltavsum > 0) then
    units(ideltavsum:ideltav+ndimV-1) = unitvel
    unitslabel(ideltavsum:ideltav+ndimV-1) = unitlabelv
 endif
 if (ideltav > 0) then
    units(ideltav:ideltav+ndimV-1) = unitvel
    unitslabel(ideltav:ideltav+ndimV-1) = unitlabelv
 endif
 if (iutherm > 0) then
    units(iutherm) = (udist/utime)**2
    unitslabel(iutherm) = ' [erg/g]'
 endif
 if (ipr > 0) then
    units(ipr) = umass/(udist*utime**2)
    unitslabel(ipr) = ' [g / (cm s^2)]'
 endif
 units(irho) = umass/udist**3
 unitslabel(irho) = ' [g/cm^3]'
 if (iBfirst > 0) then
    units(iBfirst:iBfirst+ndimV-1) = umagfd
    unitslabel(iBfirst:iBfirst+ndimV-1) = ' [G]'
 endif
 if (igrainsize > 0) then
    units(igrainsize) = udist
    unitslabel(igrainsize) = ' [cm]'
 endif
 if (igraindens > 0) then
    units(igraindens) = umass/udist**3
    unitslabel(igraindens) = ' [g/cm^3]'
 endif
 if (ivrel > 0) then
    units(ivrel) = udist/utime/100
    unitslabel(ivrel) = ' [m/s]'
 endif
 if (idivB > 0) then
    units(idivB:idivB+3) = umagfd/unitx
    unitslabel(idivB:idivB+3) = ' [G/'//trim(unitlabelx(3:))
 endif

 !--use the following two lines for time in years
 call get_environment_or_flag('SSPLASH_TIMEUNITS',string)
 select case(trim(lcase(adjustl(string))))
 case('s','seconds')
    units(0) = utime
    unitslabel(0) = trim(string)
 case('min','minutes','mins')
    units(0) = utime/60.d0
    unitslabel(0) = trim(string)
 case('h','hr','hrs','hours','hour')
    units(0) = utime/3600.d0
    unitslabel(0) = trim(string)
 case('y','yr','yrs','years','year')
    units(0) = utime/3.1536d7
    unitslabel(0) = trim(string)
 case('d','day','days')
    units(0) = utime/(3600.d0*24.d0)
    unitslabel(0) = trim(string)
 case('tff','freefall','tfreefall')
    !--or use these two lines for time in free-fall times
    units(0) = 1./tfreefall
    unitslabel(0) = ' '
 case default
    if (dtmax > 0. .and. (abs(utime-1d0) > 0.d0)) then ! use interval between dumps
       call get_nearest_time_unit(utime*dtmax,unitx,unitslabel(0))
       unitx = unitx/dtmax
    else
       call get_nearest_time_unit(utime,unitx,unitslabel(0))
    endif
    units(0) = unitx ! convert to real*4
 end select

 unitzintegration = udist
 labelzintegration = ' [cm]'
 !
 !--set labels for each particle type
 !
 if (phantomdump) then  ! phantom
    ntypes = itypemap_unknown_phantom
    labeltype(1) = 'gas'
    labeltype(2) = 'dust (old)'
    labeltype(3) = 'sink'
    labeltype(4) = 'ghost'
    labeltype(5) = 'star'
    labeltype(6) = 'dark matter'
    labeltype(7) = 'bulge'
    if (ndustlarge > 0) then
       ! try to label dust particles as 1cm dust, 10cm dust etc.
       idustlast = min(8+ndustlarge-1,size(labeltype)-1)
       if (allocated(grainsize)) then
          do i=1,ndustlarge
             j = 8 + i - 1
             if (j < size(labeltype)) then
                if (grainsize(i) > 0.) then
                   labeltype(j) = trim(get_label_grain_size(grainsize(i)))//' dust'
                else
                   labeltype(j) = 'dust'
                endif
             endif
          enddo
       else
          labeltype(8:idustlast) = 'dust'
       endif
    endif
    labeltype(itypemap_unknown_phantom) = 'unknown/dead'
    UseTypeInRenderings(:) = .true.
    UseTypeInRenderings(3) = .false.
    if (lenvironment('SSPLASH_PLOT_DUST')) then
       UseTypeInRenderings(2) = .false.
    endif
    if (lenvironment('SSPLASH_PLOT_STARS')) then
       UseTypeInRenderings(5) = .false.
    endif
    if (lenvironment('SSPLASH_PLOT_DM')) then
       UseTypeInRenderings(6) = .false.
    endif
 else
    ntypes = 6
    labeltype(1) = 'gas'
    labeltype(2) = 'dust'
    labeltype(3) = 'ghost'
    labeltype(4) = 'sink'
    labeltype(5) = 'star'
    labeltype(6) = 'unknown/dead'
    UseTypeInRenderings(1) = .true.
    UseTypeInRenderings(2) = .false.
    UseTypeInRenderings(3) = .true.
    UseTypeInRenderings(4) = .false.
    UseTypeInRenderings(5) = .true.
    UseTypeInRenderings(6) = .true.  ! only applies if turned on
 endif

 return
end subroutine set_labels_sphNG

!-----------------------------------------------------------
!
! check if a file is in phantom/sphNG format
!
!-----------------------------------------------------------
logical function file_format_is_sphNG(filename) result(is_sphNG)
 character(len=*), intent(in) :: filename
 integer :: iunit,intg1,ierr

 is_sphNG = .false.
 !
 ! open file and read the first line
 !
 open(newunit=iunit,iostat=ierr,file=filename,status='old',form='unformatted')
 if (ierr /= 0) return
 !
 ! check the magic integer to determine phantom/sphNG format
 !
 read(iunit,iostat=ierr) intg1
 if (intg1==690706 .or. intg1==060769) is_sphNG = .true.
 close(iunit)    ! close the file

end function file_format_is_sphNG

end module readdata_sphNG
