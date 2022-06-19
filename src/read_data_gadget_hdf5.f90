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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR HDF5 OUTPUT FROM THE GADGET CODE
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! GSPLASH_USE_Z if 'YES' uses redshift in the legend instead of time
! GSPLASH_USE_IDS if 'YES' resorts particles according to their ParticleIDs
! GSPLASH_DARKMATTER_HSOFT if given a value > 0.0 will assign a
!  smoothing length to dark matter particles which can then be
!  used in the rendering
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxpart,maxplot,maxstep) : main data array
!
! npartoftype(maxstep): number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------
!
!  The module below contains interface routines to c functions
!  that perform the actual calls to the HDF5 libs
!
!-------------------------------------------------------------------------
module gadgethdf5read
 use params, only:maxplot,doub_prec
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none
 real :: hsoft
 character(len=lenlabel), dimension(maxplot) :: blocklabelgas
 integer, dimension(maxplot) :: blocksize
 logical :: havewarned = .false.
 integer, parameter :: maxtypes = 6
 logical :: arepo = .false.

 interface
  subroutine read_gadget_hdf5_header(filename,maxtypes,npartoftypei,massoftypei,&
                                      timeh,zh,iFlagSfr,iFlagFeedback,Nall,iFlagCool, &
                                      igotids,ndim,ndimV,nfiles,ncol,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   integer(kind=c_int), intent(in), value :: maxtypes
   integer(kind=c_int), intent(out) :: iFlagSfr,iFlagFeedback,iFlagCool,igotids
   integer(kind=c_int), dimension(6), intent(out) :: npartoftypei,Nall
   real(kind=c_double), dimension(6), intent(out) :: massoftypei
   real(kind=c_double), intent(out) :: timeh,zh
   integer(kind=c_int), intent(out) :: ndim,ndimV,nfiles,ncol,ierr
  end subroutine read_gadget_hdf5_header

  subroutine read_gadget_hdf5_data(filename,maxtypes,npartoftypei,massoftypei,&
                                    ncol,isrequired,i0,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in)  :: filename
   integer(kind=c_int), intent(in), value :: maxtypes
   integer(kind=c_int), dimension(6), intent(in) :: npartoftypei
   real(kind=c_double), dimension(6), intent(in) :: massoftypei
   integer(kind=c_int), intent(in), value  :: ncol
   integer(kind=c_int), intent(out) :: ierr
   integer(kind=c_int), dimension(ncol), intent(in)  :: isrequired
   integer(kind=c_int), dimension(maxtypes), intent(in) :: i0
  end subroutine read_gadget_hdf5_data
 end interface

contains
 !---------------------------------------------------------------------------
 !
 ! function to safely convert a string from c format (ie. with a terminating
 ! ascii null character) back to a normal Fortran string
 !
 !---------------------------------------------------------------------------
function fstring(array)
 character(kind=c_char), dimension(:), intent(in) :: array
 character(len=size(array)-1) :: fstring
 integer :: i

 fstring = ''
 do i=1,size(array)
    if (array(i)==achar(0)) exit
    fstring(i:i) = array(i)
 enddo

end function fstring

 !---------------------------------------------------------------------------
 !
 ! function to reformat the HDF5 label into the splash column label
 ! by inserting a space whereever a capital letter occurs
 !
 !---------------------------------------------------------------------------
function reformatlabel(label)
 character(len=*), intent(in) :: label
 character(len=2*len(label)) :: reformatlabel
 integer :: is,ia,ib,ip

 reformatlabel = label
 ip = 1
 do is = 2, len_trim(label)
    ip = ip + 1
    ia = iachar(reformatlabel(ip:ip))
    ib = iachar(reformatlabel(ip-1:ip-1))
    if ((ia >= iachar('A').and.ia <= iachar('Z')) .and. .not. &
          (ib >= iachar('A').and.ib <= iachar('Z'))) then
       reformatlabel = reformatlabel(1:ip-1)//' '//reformatlabel(ip:)
       ip = ip + 1
    endif
 enddo

end function reformatlabel

end module gadgethdf5read

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------


module readdata_gadget_hdf5
 implicit none

 public :: read_data_gadget_hdf5, set_labels_gadget_hdf5

 private
contains

subroutine read_data_gadget_hdf5(rootname,istepstart,ipos,nstepsread)
 use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep
 use params,         only:doub_prec,maxparttypes,maxplot
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,required,ipartialread, &
                           ntypes,debugmode,iverbose
 use mem_allocation, only:alloc
 use labels,         only:ih,irho,ipmass
 use system_utils,   only:renvironment,lenvironment,ienvironment,envlist
 use asciiutils,     only:cstring
 use gadgethdf5read, only:hsoft,blocklabelgas,havewarned,read_gadget_hdf5_header, &
                           read_gadget_hdf5_data,maxtypes,arepo
 integer, intent(in)                :: istepstart,ipos
 integer, intent(out)               :: nstepsread
 character(len=*), intent(in)       :: rootname
 character(len=len(rootname)+10)    :: datfile,densfile,hfile
 character(len=20)                  :: string
 integer, dimension(maxparttypes)   :: npartoftypei,Nall
 integer               :: i,j,itype,ierr,ierrh,ierrrho,nhset,ifile
 integer               :: index1,index2
 integer               :: ncolstep,npart_max,nstep_max,ntoti,ntotall,idot
 integer               :: iFlagSfr,iFlagFeedback,iFlagCool,igotids,nfiles,nhfac
 integer, dimension(6) :: i0
 integer, parameter    :: iunit = 11, iunitd = 102, iunith = 103
 logical               :: iexist,reallocate,debug,goterrors
 real(doub_prec)                    :: timetemp,ztemp
 real(doub_prec), dimension(6)      :: massoftypei
 real :: hfact,hfactmean,pmassi
 real, parameter :: pi = 3.1415926536
 integer, dimension(maxplot) :: isrequired

 nstepsread = 0
 goterrors  = .false.
 if (maxparttypes < 6) then
    print*,' *** ERROR: not enough particle types for GADGET data read ***'
    print*,' *** you need to edit splash parameters and recompile ***'
    stop
 endif

 if (len_trim(rootname) > 0) then
    datfile = trim(rootname)
 else
    print*,' **** no data read **** '
    return
 endif
!
!--check if first data file exists
!
 print "(1x,a)",'reading GADGET HDF5 format'
 inquire(file=datfile,exist=iexist)
 if (.not.iexist) then
    !
    !--look for a file with .0 on the end for multiple-file reads
    !
    datfile=trim(rootname)//'.0.hdf5'
    inquire(file=datfile,exist=iexist)
    if (.not.iexist) then
       print "(a)",' *** error: '//trim(rootname)//': file not found ***'
       return
    endif
 endif
!
!--set parameters which do not vary between timesteps
!
 ndim  = 0
 ndimV = 0
!  idumpformat = ienvironment('GSPLASH_FORMAT')
!  checkids    = lenvironment('GSPLASH_CHECKIDS')
 debug       = lenvironment('GSPLASH_DEBUG') .or. debugmode
!
!--read data from snapshots
!
 i = istepstart
!
!--i0 is the offset used to read the data into the arrays
!  (non-zero for read from multiple files)
!  The offset is different for each particle type, somewhat
!  complicating the data read -- we shuffle the particles from
!  multiple files so that they are in type order.
!
 i0(:) = 0
!
!--loop over the number of files
!
 ifile = 0
 ntotall = 0
 over_files: do while(iexist)

    write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
    ifile = ifile + 1

    !
    !--open file and read header information
    !
    npartoftypei(:) = 0.
    Nall(:) = 0.
    massoftypei(:) = 0.
    if (debug) print*,'DEBUG: reading header...'
    call read_gadget_hdf5_header(cstring(datfile),maxtypes, &
       npartoftypei,massoftypei,timetemp,ztemp,iFlagSfr,iFlagFeedback,Nall,&
       iFlagCool,igotids,ndim,ndimV,nfiles,ncolstep,ierr)
    if (ierr /= 0) then
       print "(a)", '*** ERROR READING HEADER ***'
       return
    endif

    ! read(iunit,iostat=ierr) npartoftypei(1:6),massoftypei,timetemp,ztemp, &
    !     iFlagSfr,iFlagFeedback,Nall(1:6),iFlagCool,nfiles

    ntoti = int(sum(npartoftypei(1:6)))  ! int here is unnecessary, but avoids compiler warnings

    if (nfiles > 1) then
       ntotall = int(sum(Nall(1:6)))
    else
       ntotall = ntoti
    endif
    !
    !--if we are reading from multiple files,
    !  check that the sequence starts from the correct file
    !
    if (nfiles > 1) then
       idot = index(datfile,'.hdf5')
       idot = index(datfile(1:idot-1),'.',back=.true.)
       if (ifile==1 .and. datfile(idot:idot+1) /= '.0') then
          if (nfiles < 100) then
             string = "(/,a,i2,a,/,a,/)"
          else
             string = "(/,a,i7,a,/,a,/)"
          endif
          print string,' ERROR: read is from multiple files (nfiles = ',nfiles,')',&
                  '        but this is not the first file (does not end in .0.hdf5): skipping...'
          close(iunit)
          return
       endif
    endif

    if (ifile==1) then
       ncolumns = ncolstep
       !
       !--call set labels to get ih, ipmass, irho for use in the read routine
       !
       hsoft = 0. ! to avoid unset variable
       call set_labels_gadget_hdf5
    endif

    if (ifile==1) then
       print*,'time            : ',timetemp
       print "(1x,a,f8.2)",'z (redshift)    : ',ztemp
    endif
    print "(a,6(1x,i10))",' Npart (by type) : ',npartoftypei(1:6)
    if (any(massoftypei > 0.)) print "(a,6(1x,es10.3))",' Mass  (by type) : ',massoftypei
    print "(a,6(1x,i10))",' N_gas           : ',npartoftypei(1)
    print "(a,1x,i10)",' N_total         : ',ntoti
    if (ifile==1) print "(a,1x,i10)",' N data columns  : ',ncolstep
    if (nfiles > 1 .and. ifile==1) then
       print "(a,6(1x,i10))",' Nall            : ',Nall(1:6)
    endif

    if (nfiles > 1) then
       if (ifile==1) print "(a,i4,a)",' reading from ',nfiles,' files'
    elseif (nfiles <= 0) then
       print*,'*** ERROR: nfiles = ',nfiles,' in file header: aborting'
       return
    endif

    if (ifile==1) then
       !--Softening lengths for Dark Matter Particles...
       hsoft = renvironment('GSPLASH_DARKMATTER_HSOFT')
       !
       !--try to read dark matter and star particle smoothing lengths and/or density from a separate
       !  one column ascii file. If only density, use this to compute smoothing lengths.
       !
       densfile = trim(rootname)//'.dens'
       hfile = trim(rootname)//'.hsml'
       hfact = 1.2 ! related to the analytic neighbour number (hfact=1.2 gives 58 neighbours in 3D)
       open(unit=iunitd,file=densfile,iostat=ierrrho,status='old',form='formatted')
       open(unit=iunith,file=hfile,iostat=ierrh,status='old',form='formatted')
       if (ih==0 .and. (hsoft > tiny(hsoft) .or. ierrrho==0 .or. ierrh==0)) then
          ncolumns = ncolumns + 1
          blocklabelgas(ncolumns) = 'SmoothingLength'
          ih = ncolumns
          call set_labels_gadget_hdf5
       endif
       if (irho==0 .and. (hsoft > tiny(hsoft) .or. ierrrho==0 .or. ierrh==0)) then
          ncolumns = ncolumns + 1
          blocklabelgas(ncolumns) = 'Density'
          irho = ncolumns
          call set_labels_gadget_hdf5
       endif
       !
       !--if successfully read header, increment the nstepsread counter
       !
       nstepsread = nstepsread + 1
    endif
    !
    !--now read data
    !
    reallocate = .false.
    npart_max = maxpart
    nstep_max = max(maxstep,1)

    if (ntoti > maxpart) then
       reallocate = .true.
       if (maxpart > 0) then
          ! if we are reallocating, try not to do it again
          npart_max = int(1.1*ntotall)
       else
          ! if first time, save on memory
          npart_max = int(ntotall)
       endif
    endif
    if (i >= maxstep .and. i /= 1) then
       nstep_max = i + max(10,INT(0.1*nstep_max))
       reallocate = .true.
    endif
    !
    !--reallocate memory for main data array
    !
    if (reallocate .or. .not.(allocated(dat))) then
       if (igotids==1) then
          call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)
       else
          call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol))
       endif
    endif
    masstype(1:6,i) = massoftypei(1:6)
    !
    !--copy npartoftypei into allocated header arrays
    !  and set the offset position of particle types in the main data arrays
    !
    if (nfiles==1 .or. ifile==1) then
       i0(1) = 0
       do itype=2,ntypes
          if (nfiles==1) then
             i0(itype) = sum(npartoftypei(1:itype-1)) ! this is avoid depending on Nall at all for single file read
          else
             i0(itype) = sum(Nall(1:itype-1))
          endif
       enddo
       npartoftype(:,i) = npartoftypei
    else
       i0(1) = npartoftype(1,i)
       do itype=2,ntypes
          i0(itype) = sum(Nall(1:itype-1)) + npartoftype(itype,i)
       enddo
       npartoftype(:,i) = npartoftype(:,i) + npartoftypei
    endif
    if (debugmode) print*,'DEBUG: starting position for each type in data array: ',i0(:)
    !
    !--set time to be used in the legend
    !
    if (ifile==1) then
       !--use this line for code time
       time(i) = real(timetemp)
    else
       if (abs(real(timetemp)-time(i)) > tiny(0.)) print*,'ERROR: time different between files in multiple-file read '
       if (sum(Nall) /= ntotall) then
          print*,' ERROR: Nall differs between files'
          goterrors = .true.
       endif
    endif
    !
    !--read particle data
    !
    got_particles: if (ntoti > 0) then

       isrequired(:) = 0
       where (required(1:ncolumns)) isrequired(1:ncolumns) = 1

       call read_gadget_hdf5_data(cstring(datfile),maxtypes,npartoftypei,massoftypei,ncolumns,isrequired,i0,ierr)

    endif got_particles
!
!--now memory has been allocated, set arrays which are constant for all time
!
    gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read
!
    if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--for read from multiple files, work out the next file in the sequence
!
    iexist = .false.
    if (nfiles > 1 .and. ifile < nfiles) then
       !--see if the next file exists
       idot = index(datfile,'.hdf5')
       idot = index(datfile(1:idot-1),'.',back=.true.)
       if (idot <= 0) then
          print "(a)",' ERROR: read from multiple files but could not determine next file in sequence'
          goterrors = .true.
       else
          write(string,*) ifile
          if (ifile < 10) then
             write(datfile,"(a,i1)") trim(datfile(1:idot))//trim(adjustl(string))//'.hdf5'
          elseif (ifile < 100) then
             write(datfile,"(a,i2)") trim(datfile(1:idot))//trim(adjustl(string))//'.hdf5'
          else
             write(datfile,"(a,i3)") trim(datfile(1:idot))//trim(adjustl(string))//'.hdf5'
          endif
          iexist = .false.
          inquire(file=datfile,exist=iexist)
          if (.not.iexist) then
             print "(a)",' ERROR: read from multiple files '// &
           'but could not find '//trim(datfile)//': next in sequence'
             goterrors = .true.
          endif
       endif
    endif

 enddo over_files

 if (arepo) then
    dat(1:npartoftype(1,i),ih,i) = dat(1:npartoftype(1,i),ih,i)**(1./3.)
    print "(a)",' this is an AREPO snapshot, using Volume**(1./3.) as smoothing length'
 elseif (required(ih) .and. size(dat(1,:,:)) >= ih .and. npartoftype(1,i) > 0) then
    !
    !--for some reason the smoothing length output by GADGET is
    !  twice the usual SPH smoothing length
    !  (do this after we have read data from all of the files)
    !
    print "(a)",' converting GADGET h on gas particles to usual SPH definition (x 0.5)'
    dat(1:npartoftype(1,i),ih,i) = 0.5*dat(1:npartoftype(1,i),ih,i)
 endif

 if (nfiles > 1. .and. any(npartoftype(:,i) /= Nall(:))) then
    print*,'ERROR: sum of Npart across multiple files  /=  Nall in data read '
    print*,'Npart = ',npartoftype(:,i)
    print*,'Nall  = ',Nall(:)
    goterrors = .true.
 endif
 !
 !--look for dark matter smoothing length/density files
 !
 if (ierrh==0 .or. ierrrho==0) then
    if (ierrh==0) then
       print "(a)",' READING DARK MATTER SMOOTHING LENGTHS from '//trim(hfile)
       ierr = 0
       index1 = npartoftype(1,i)+1
       index2 = npartoftype(1,i)+sum(npartoftype(2:,i))
       read(iunith,*,iostat=ierr) (dat(j,ih,i),j=index1,index2)
       close(unit=iunith)
       if (ierr < 0) then
          nhset = 0
          do j=index1,index2
             if (dat(j,ih,i) > 0.) nhset = nhset + 1
          enddo
          print "(a,i10,a,/)",' *** END-OF-FILE: GOT ',nhset,' SMOOTHING LENGTHS ***'
       elseif (ierr > 0) then
          print "(a)", ' *** ERROR reading smoothing lengths from file'
          goterrors = .true.
       else
          print "(a,i10,a)",' SMOOTHING LENGTHS READ OK for ',index2-index1+1,' dark matter / star particles '
       endif
       hsoft = 1.0 ! just so dark matter rendering is allowed in set_labels_gadget_hdf5 routine
    endif

    if (ierrrho==0) then
       print "(a)",' READING DARK MATTER DENSITIES FROM '//trim(densfile)
       ierr = 0
       index1 = npartoftype(1,i)+1
       index2 = npartoftype(1,i)+sum(npartoftype(2:,i))
       read(iunitd,*,iostat=ierr) (dat(j,irho,i),j=index1,index2)
       close(iunitd)
       if (ierr < 0) then
          nhset = 0
          do j=index1,index2
             if (dat(j,irho,i) > 0.) nhset = nhset + 1
          enddo
          print "(a,i10,a,/)",' *** END-OF-FILE: GOT ',nhset,' DENSITIES ***'
       elseif (ierr > 0) then
          print "(a)", ' *** ERROR reading dark matter densities from file'
          goterrors = .true.
       else
          print "(a,i10,a)",' DENSITY READ OK for ',index2-index1+1,' dark matter / star particles '
       endif
       if (ierrh /= 0 .and. ipmass > 0) then
          where(dat(:,irho,i) > tiny(dat))
             dat(:,ih,i) = hfact*(dat(:,ipmass,i)/dat(:,irho,i))**(1./3.)
          elsewhere
             dat(:,ih,i) = 0.
          end where
          print "(a,i10,a,f5.2,a)", &
            ' SMOOTHING LENGTHS SET for ',index2-index1+1,' DM/star particles using h = ',hfact,'*(m/rho)**(1/3)'
       endif
       hsoft = 1.0 ! just so dark matter rendering is allowed in set_labels_gadget_hdf5 routine
    endif
 else
    !
    !--if a value for the dark matter smoothing length is set
    !  via the environment variable GSPLASH_DARKMATTER_HSOFT,
    !  give dark matter particles this smoothing length
    !  and a density of 1 (so column density plots work)
    !
    if (hsoft > tiny(hsoft)) then
       if (required(ih)) then
          print "(a,1pe10.3,a)",' ASSIGNING SMOOTHING LENGTH of h = ',hsoft, &
                                 ' to dark matter particles'
          !print*,'ih = ',ih,' npartoftype = ',npartoftype(1:2,i), shape(dat)
          if (ih > 0) then
             dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),ih,i) = hsoft
          else
             print*,' ERROR: smoothing length not found in data arrays'
             goterrors = .true.
          endif
       endif
       if (required(irho)) then
          if (irho > 0) then
             dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),irho,i) = 1.0
          else
             print*,' ERROR: place for density not found in data arrays'
             goterrors = .true.
          endif
       endif
    else
       if (npartoftype(1,i) <= 0 .and. sum(npartoftype(:,i)) > 0) then
          print "(66('*'),4(/,a),/)",'* NOTE!! For GADGET data using dark matter only, column density ',&
                              '* plots can be produced by setting the GSPLASH_DARKMATTER_HSOFT ',&
                              '* environment variable to give the dark matter smoothing length', &
                              '* (for a fixed smoothing length)'
          hsoft = (maxval(dat(:,1,i)) - minval(dat(:,1,i)))/sum(npartoftype(2:,i))**(1./3.)
          print*,' suggested value for GSPLASH_DARKMATTER_HSOFT = ',hsoft
          hsoft = 0.

          print "(7(/,a),/)",'* Alternatively, and for best results, calculate a number density', &
                                        '* on dark matter particles, set individual smoothing lengths from', &
                                        '* this using h = hfact*(n)**(-1/3), with hfact=1.2 and either ', &
                                        '* dump the results back into the HSML array in the original dump ', &
                                        '* file (if using the block-labelled format), or create an ascii ',&
                                        '* file called '//trim(hfile)//' containing the smoothing length ',&
                                        '* values for the dark matter particles.'
          print "(2(/,a),/,66('*'),/)",  '* Also make sure normalised interpolations are OFF when plotting ',&
                                        '* dark matter density '
       endif
    endif
 endif
!
!--pause with fatal errors
!
 if (goterrors .and. .not.lenvironment('GSPLASH_IGNORE_ERRORS')) then
    print "(/,a)",'*** ERRORS detected during data read: data will be corrupted'
    print "(a,/)",'    Please REPORT this and/or fix your file ***'
    print "(a)",'     (set GSPLASH_IGNORE_ERRORS=yes to skip this message)'
    if (iverbose >= 1) then
       print "(a)",'    > Press any key to bravely proceed anyway  <'
       read*
    endif
 endif
!
!--give a friendly warning about using too few or too many neighbours
!  (only works with equal mass particles because otherwise we need the number density estimate)
!
 if (ih > 0 .and. required(ih) .and. ipmass > 0 .and. required(ipmass) &
      .and. abs(massoftypei(1)) < tiny(0.) .and. ndim==3 .and. .not.havewarned) then
    nhfac = 100
    if (npartoftype(1,i) > nhfac) then
       hfactmean = 0.
       do j=1,nhfac
          pmassi = dat(j,ipmass,i)
          if (pmassi > 0.) then
             pmassi = 1./pmassi
          else
             pmassi = 0.
          endif
          hfact = dat(j,ih,i)*(dat(j,irho,i)*pmassi)**(1./ndim)
          hfactmean = hfactmean + hfact
       enddo
       hfact = hfactmean/real(nhfac)
       havewarned = .true.
       ! CHECK that h is not implausibly small due to the multiply-by-0.5 -- if so, reverse this
       if (hfact < 0.75) then
          dat(1:npartoftype(1,i),ih,i) = dat(1:npartoftype(1,i),ih,i)*2.0
          hfact = hfact*2.0
       endif
       if (hfact < 1.125 .or. hfact > 1.45) then
          print "(/,a)",'** FRIENDLY NEIGHBOUR WARNING! **'
          print "(3x,a,f5.1,a,/,3x,a,f4.2,a,i1,a)", &
                 'It looks like you are using around ',4./3.*pi*(2.*hfact)**3,' neighbours,', &
                 'corresponding to h = ',hfact,'*(m/rho)^(1/',ndim,') in 3D:'

          if (hfact < 1.15) then
             print "(4(/,3x,a))",'This is a quite a low number of neighbours for the cubic spline and ', &
                                  'may result in increased noise and inaccurate wave propagation speeds', &
                                  '(a cubic lattice is also an unstable initial configuration for the ',&
                                  ' particles in this regime -- see Morris 1996, Borve et al. 2004).'
          elseif (hfact > 1.45) then
             print "(4(/,3x,a))",'Using h >~ 1.5*(m/rho)^(1/3) with the cubic spline results in the', &
                                'particle pairing instability due to the first neighbour being placed under', &
                                'the hump in the kernel gradient. Whilst not fatal, it results in a', &
                                'loss of resolution so is a bit of a waste of cpu time.'
             print "(4(/,3x,a))",'If you are attempting to perform a "resolution study" by increasing the', &
                                'neighbour number, this is a *bad idea*, as you are also increasing h.',      &
                                '(a better way is to increase the smoothness of the integrals without changing h',   &
                                ' by adopting a smoother kernel such as the M6 Quintic that goes to 3h).'
          endif
          print "(/,3x,a,/,3x,a,/)", &
              'A good default is h = 1.2 (m/rho)^1/ndim ', &
              'corresponding to around 58 neighbours in 3D.'
       else
          print "(/,1x,a,f5.1,a,/,1x,a,f4.2,a,i1,a,/)", &
                'Simulations employ ',4./3.*pi*(2.*hfact)**3,' neighbours,', &
                'corresponding to h = ',hfact,'*(m/rho)^(1/',ndim,') in 3D'
       endif
    endif
 else
    !print*,'not true'
 endif
!
!--cover the special case where no particles have been read
!
 if (ntotall <= 0) then
    npartoftype(1,i) = 1
    dat(:,:,i) = 0.
 endif

 if (nstepsread > 0) then
    print "(a,i10,a)",' >> read ',sum(npartoftype(:,istepstart+nstepsread-1)),' particles'
 endif

end subroutine read_data_gadget_hdf5

subroutine read_gadgethdf5_data_fromc(icol,npartoftypei,temparr,id,itype,i0) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int,c_double
 use particle_data,  only:dat,iamtype
 use settings_data,  only:debugmode
 use labels,         only:label
 use system_utils,   only:lenvironment
 integer(kind=c_int), intent(in) :: icol,npartoftypei,itype,i0
 real(kind=c_double), dimension(npartoftypei), intent(in) :: temparr
 integer(kind=c_int), dimension(npartoftypei), intent(in) :: id
 integer(kind=c_int) :: i,icolput
 integer :: nmax,nerr,idi
 logical :: useids

 icolput = icol
 if (debugmode) print "(a,i2,a,i2,a,i8)",'DEBUG: reading column ',icol,' type ',itype,' -> '//trim(label(icolput))//', offset ',i0
 if (icolput > size(dat(1,:,1)) .or. icolput==0) then
    print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in receive_data_fromc'
    return
 endif
 nmax = size(dat(:,1,1))

 useids = lenvironment('GSPLASH_USE_IDS') .or. lenvironment('GSPLASH_CHECKIDS')
 if (all(id <= 0) .or. size(iamtype(:,1)) <= 1) useids = .false.
 if (debugmode) print*,'DEBUG: using particle IDs = ',useids,' max = ',nmax

 if (useids) then
    nerr = 0
    !print*,' id range is ',minval(id),' to ',maxval(id),' type ',itype+1,' column = ',trim(label(icolput))
    do i=1,npartoftypei
       if (id(i) < 1 .or. id(i) > nmax) then
          idi = id(i)
          !
          !--correct for particle IDs > 1e9 (used to represent recycled particles?)
          !
          if (idi > 1000000000) then
             idi = idi - 1000000000
             if (idi <= nmax .or. idi <= 0) then
                dat(idi,icolput,1) = real(temparr(i))
                iamtype(idi,1) = itype + 1
             else
                nerr = nerr + 1
                if (debugmode .and. nerr <= 10) print*,i,'fixed id = ',idi
             endif
          else
             nerr = nerr + 1
             if (debugmode .and. nerr <= 10) print*,i,' id = ',idi,idi-1000000000
          endif
       else
          dat(id(i),icolput,1) = real(temparr(i))
          iamtype(id(i),1) = itype + 1
       endif
    enddo
    if (nerr > 0) print*,'ERROR: got particle ids outside array dimensions ',nerr,' times'
 else
    if (i0 < 0) then
       print*,'ERROR: i0 = ',i0,' but should be positive: SOMETHING IS VERY WRONG...'
       return
    elseif (i0+npartoftypei > nmax) then
       print "(a,i8,a)",' ERROR: offset = ',i0,': read will exceed array dimensions in receive_data_fromc'
       nmax = nmax - i0
    else
       nmax = npartoftypei
    endif
    do i=1,nmax
       dat(i0+i,icolput,1) = real(temparr(i))
    enddo
    if (size(iamtype(:,1)) > 1) then
       do i=1,nmax
          iamtype(i0+i,1) = itype + 1
       enddo
    endif
 endif

end subroutine read_gadgethdf5_data_fromc

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels_gadget_hdf5
 use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass, &
                          ih,irho,ipr,iutherm,iBfirst,iax,make_vector_label
 use params
 use settings_data,  only:ndim,ndimV,ntypes,UseTypeInRenderings,debugmode
 use geometry,       only:labelcoord
 use system_utils,   only:envlist,ienvironment
 use gadgethdf5read, only:hsoft,blocklabelgas,blocksize,reformatlabel,arepo
 use asciiutils,     only:lcase
 integer :: i,icol,irank

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_gadget_hdf5 ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_gadget_hdf5 ***'
    return
 endif

 icol = 1
 ix = 0
 do i=1,size(blocklabelgas)
    irank = blocksize(i)
    if (irank > 0 .and. (len_trim(blocklabelgas(i)) > 0)) then
       if (debugmode) print*,i,trim(blocklabelgas(i))
       select case(blocklabelgas(i))
       case('Coordinates')
          ix(1) = icol
          ix(2) = icol + 1
          if (irank >= 3) ix(3) = icol + 2
       case('Velocities','Velocity')
          ivx = icol
       case('SmoothingLength','SmoothingLengths')
          ih = icol
       case('Volume')
          ih = icol
          arepo = .true.
       case('Masses','Mass')
          ipmass = icol
       case('InternalEnergy','InternalEnergies')
          iutherm = icol
       case('Density','Densities')
          irho = icol
       case('MagneticField')
          iBfirst = icol
       case('Pressures','Pressure')
          ipr = icol
       end select
       label(icol:icol+irank-1) = reformatlabel(blocklabelgas(i))

       if (irank==ndimV) then
          call make_vector_label(label(icol),icol,ndimV,iamvec,labelvec,label,labelcoord(:,1))
       endif
       icol = icol + irank
    endif
 enddo
 !
 !--set labels of the quantities read in
 !
 if (ix(1) > 0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
 if (irho > 0)    label(irho)       = 'density' ! needed to convert to "column density"
 if (iutherm > 0) label(iutherm)    = 'u'       ! otherwise \int u dz looks ugly
 !if (ipmass > 0)  label(ipmass)     = 'particle mass'
 !if (ih > 0)      label(ih)         = 'h'
 !
 !--set labels for vector quantities
 !
 call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('a',iax,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('B',iBfirst,ndimV,iamvec,labelvec,label,labelcoord(:,1))

 !--set labels for each particle type
 !
 ntypes = 6
 labeltype(1) = 'gas'
 labeltype(2) = 'dark matter'
 labeltype(3) = 'boundary 1'
 labeltype(4) = 'boundary 2'
 labeltype(5) = 'star'
 labeltype(6) = 'sink / black hole'
 UseTypeInRenderings(1) = .true.
 !
 !--dark matter particles are of non-SPH type (ie. cannot be used in renderings)
 !  unless they have had a smoothing length defined
 !
 if (hsoft > tiny(hsoft)) then
    UseTypeInRenderings(2) = .true.
 else
    UseTypeInRenderings(2) = .false.
 endif
 UseTypeInRenderings(3:6) = .false.

end subroutine set_labels_gadget_hdf5

subroutine set_blocklabel_gadget(icol,irank,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 use gadgethdf5read, only:blocklabelgas,blocksize,fstring
 integer(kind=c_int), intent(in) :: icol,irank
 character(kind=c_char), dimension(256), intent(in) :: name

 blocklabelgas(icol+1) = fstring(name)
 blocksize(icol+1) = irank
 !print*,icol+1,' name = ',trim(blocklabelgas(icol+1)),' x ',irank

end subroutine set_blocklabel_gadget
end module readdata_gadget_hdf5
