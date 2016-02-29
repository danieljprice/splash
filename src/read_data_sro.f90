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
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM STEPHAN ROSSWOG'S CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! RSPLASH_FORMAT can be 'MHD' or 'HYDRO'
! RSPLASH_RESET_COM if 'YES' then centre of mass is reset for n2=0 (ie single objects)
! RSPLASH_COROTATING if 'YES' then velocities are transformed to corotating frame
! RSPLASH_HFACT can be changed to give correct hfact value for particle masses
!  on minidumps: e.g. setenv RSPLASH_HFACT=1.2
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
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data, only:dat,time,npartoftype,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,iformat
  use mem_allocation, only:alloc
  use system_utils, only:lenvironment
  use system_commands, only:get_environment
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer, parameter :: max_spec = 7 ! number of species in abundance file
  real :: hfact, dhfact3,hfacttemp
  integer :: i,j,k,ierr,ierr1
  integer :: nprint,nptmass,npart_max,nstep_max
  integer :: n1,n2,idump,ncol
  logical :: iexist,magfield,minidump,doubleprec,iabunfileopen
  character(len=len(rootname)) :: dumpfile
  character(len=13) :: abunfile
  character(len=10) :: string
  real :: timei,tkin,tgrav,tterm,escap,rstar,mstar,Etot_burn_cgs
  real(doub_prec) :: timedb,tkindb,tgravdb,ttermdb
  real(doub_prec) :: escapdb,rstardb,mstardb,Etot_burn_cgsdb
  real(doub_prec), dimension(:,:), allocatable :: datdb

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  iabunfileopen = .false.
  hfact = 1.5
  dhfact3 = 1./hfact**3

  dumpfile = trim(rootname)
  !
  !--check if first data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
     return
  endif
  !
  !--use minidump format if minidump
  !
  minidump = .false.
  if (index(dumpfile,'minidump').ne.0) minidump = .true.
  !
  !--get hfact for minidumps
  !
  if (minidump) then
     call get_environment('RSPLASH_HFACT',string)
     read(string,*,iostat=ierr) hfacttemp
     if (hfacttemp.gt.0.5 .and. hfacttemp .lt.10.0 .and. ierr.eq.0) then
        hfact = hfacttemp
        print *,'setting hfact =',hfact,' from RSPLASH_HFACT environment variable'
     else
        print "(1x,a)",'error reading hfact from RSPLASH_HFACT environment variable'
     endif
  endif
  !
  !--try to guess full dump format from file names
  !
  magfield = .true.
!  if (.not.minidump) then
!     if (index(dumpfile,'SMBH').gt.0) magfield = .false.
!     if (index(dumpfile,'nsbh').gt.0) magfield = .false.
!     if (index(dumpfile,'NSBH').gt.0) magfield = .false.
!     if (index(dumpfile,'WD').gt.0) magfield = .false.
!  endif
  !
  !--override this with environment variable
  !
  call get_environment('RSPLASH_FORMAT',string)
  select case(trim(adjustl(string)))
  case('MHD','mhd','ns','NS')
     magfield = .true.
  case('WD','hydro','HYDRO','ns_bh_v2')
     magfield = .false.
  end select

  if (magfield) then
     print "(1x,a)",'reading MAGMA code format (set RSPLASH_FORMAT=hydro for hydro format)'
  else
     print "(1x,a)",'reading Stephan Rosswog (hydro) code format (set RSPLASH_FORMAT=MHD for MAGMA)'
  endif

  !
  !--fix number of spatial dimensions
  !
  ndim = 3
  ndimV = 3
  if (magfield) then
     if (minidump) then
        ncol = 11
        iformat = 1
     else
        ncol = 27
        iformat = 2
     endif
  else
     if (minidump) then
        ncol = 7  ! number of columns in file
        iformat = 3
     else
        ncol = 16
        iformat = 4
     endif
  endif
  n1 = 0
  n2 = 0
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,1)

  j = indexstart
  nstepsread = 0

  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
     return
  else
     !
     !--read the number of particles in the first step,
     !  allocate memory and rewind
     !
     doubleprec = .true.
     if (minidump) then
        !--try double precision first
        read(15,end=55,iostat=ierr) timedb,nprint,nptmass
        !--change to single precision if stupid answers
        if (nprint.le.0.or.nprint.gt.1e10 &
            .or.nptmass.lt.0.or.nptmass.gt.1e6) then
           doubleprec = .false.
           rewind(15)
           read(15,end=55,iostat=ierr) timei,nprint,nptmass
           if (magfield) then
              print "(a)",' single precision MHD minidump'
           else
              print "(a)",' single precision hydro minidump'
           endif
        else
           if (magfield) then
              print "(a)",' double precision MHD minidump'
           else
              print "(a)",' double precision hydro minidump'
           endif
           timei = real(timedb)
        endif
     else
        !--try double precision first
        read(15,end=55,iostat=ierr) nprint,rstardb,mstardb,n1,n2, &
                nptmass,timedb
        !--change to single precision if stupid answers
        if (n1.lt.0.or.n1.gt.1e10.or.n2.lt.0.or.n2.gt.1e10 &
           .or.nptmass.lt.0.or.nptmass.gt.1.e6) then
           doubleprec = .false.
           rewind(15)
           read(15,end=55,iostat=ierr) nprint,rstar,mstar,n1,n2, &
                nptmass,timei
           if (magfield) then
              print "(a)",' single precision full MHD dump'
           else
              print "(a)",' single precision full hydro dump'
           endif
        else
           if (magfield) then
              print "(a)",' double precision full MHD dump'
           else
              print "(a)",' double precision full hydro dump'
           endif
           timei = real(timedb)
        endif
     endif
     print "(a,f10.2,a,i9,a,i6)",' time: ',timei,' npart: ',nprint,' nptmass: ',nptmass
     !--barf if stupid answers in single and double precision
     if (nptmass.lt.0.or.nptmass.gt.1.e6 .or. nprint.lt.0 &
         .or. nprint.gt.1e10 .or. (nprint.eq.0 .and. nptmass.eq.0)) then
        print "(a)",' *** ERRORS IN TIMESTEP HEADER: NO DATA READ ***'
        close(15)
        return
     endif

     ncolumns = ncol
     !
     !--check if abundance files are present and read from them
     !
     if (.not.magfield) then
        !--extract dump number from filename (last 5 characters)
        read(dumpfile(len_trim(dumpfile)-4:len_trim(dumpfile)),*,iostat=ierr1) idump
        if (ierr1 /= 0) then
           print "(a)",' error extracting dump number from filename'
        else
           write(abunfile,"(a,'.',i5.5)") 'abun_a7',idump
           inquire(file=abunfile,exist=iexist)
           if (.not.iexist) then
              print "(a)",' abundance file '//trim(abunfile)//' NOT FOUND'
           else
              open(unit=41,file=abunfile,status='old',form='unformatted',iostat=ierr1)
              if (ierr1 /= 0) then
                 print "(a)",'*** ERROR OPENING '//trim(abunfile)
              else
                 ncolumns = ncol + max_spec + 2
                 iabunfileopen = .true.
              endif
           endif
        endif
     endif

     if (.not.allocated(dat) .or. (nprint+nptmass).gt.npart_max) then
        npart_max = max(npart_max,INT(1.1*(nprint+nptmass)))
        call alloc(npart_max,nstep_max,ncolumns)
     endif
     rewind(15)
  endif
  if (ierr /= 0) then
     print "(a)",'*** ERROR READING TIMESTEP HEADER ***'
  else
!
!--loop over the timesteps in this file
!
     npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
     if (j.gt.maxstep) then
        call alloc(maxpart,j+1,maxcol)
     endif
!
!--now read the timestep data in the dumpfile
!
     if (magfield) then
        if (minidump) then
           dat(:,:,j) = 0.

           if (doubleprec) then
              allocate(datdb(maxpart,ncol),stat=ierr)
              if (ierr /= 0) then
                 print*,"(a)",'*** error allocating memory for double conversion ***'
                 return
              else
                 datdb = 0.
              endif
              read(15,end=55,iostat=ierr) timedb,nprint,nptmass, &
                ((datdb(i,k),i=1,nprint),k=1,6),  &
                ((datdb(i,k),i=1,nprint),k=8,ncol), &
                (datdb(i,7), i=nprint+1, nprint+nptmass), &
                ((datdb(i,k), i=nprint+1, nprint+nptmass),k=1,3)

                if (ierr /= 0) print "(a)",'*** WARNING: ERRORS DURING READ ***'
                dat(:,1:ncol,j) = real(datdb(:,1:ncol))
                time(j) = real(timedb)

           else
              read(15,end=55,iostat=ierr) time(j),nprint,nptmass, &
                ((dat(i,k,j),i=1,nprint),k=1,6), &
                ((dat(i,k,j),i=1,nprint),k=8,ncol), &
                (dat(i,7,j), i=nprint+1, nprint+nptmass), &
                ((dat(i,k,j), i=nprint+1, nprint+nptmass),k=1,3)
           endif
           !
           !--because masses are not dumped, we need to reconstruct them
           !  from density and h (only strictly true for grad h code)
           !
           dat(1:nprint,7,j) = dat(1:nprint,4,j)**3*dat(1:nprint,5,j)*dhfact3
           print "(a,f5.2,a)", &
            ' WARNING: setting particle masses assuming h = ',hfact,'*(m/rho)^(1/3)'
        else

           dat(:,:,j) = 0. ! because ptmasses don't have all quantities
        !
        !--read full dump
        !
           if (doubleprec) then
              allocate(datdb(maxpart,27),stat=ierr)
              if (ierr /= 0) then
                 print*,"(a)",'*** error allocating memory for double conversion ***'
                 return
              else
                 datdb = 0.
              endif
              read(15,iostat=ierr) nprint,rstardb,mstardb,n1,n2, &
                nptmass,timedb,(datdb(i,7),i=1,nprint), &
                escapdb,tkindb,tgravdb,ttermdb, &
                ((datdb(i,k),i=1,nprint),k=1,6), &
                ((datdb(i,k),i=1,nprint),k=8,ncol), &
                (datdb(i,9), i=nprint+1, nprint+nptmass), &
                ((datdb(i,k), i=nprint+1, nprint+nptmass),k=1,6), &
                ((datdb(i,k), i=nprint+1, nprint+nptmass),k=21,23)

                if (ierr < 0) print "(a)",'*** WARNING: END OF FILE DURING READ ***'
                if (ierr > 0) print "(a)",'*** WARNING: ERRORS DURING READ ***'
                dat(:,1:ncol,j) = real(datdb(:,1:ncol))
                time(j) = real(timedb)

           else
              read(15,iostat=ierr) nprint,rstar,mstar,n1,n2, &
                nptmass,time(j),(dat(i,7,j),i=1,nprint), &
                escap,tkin,tgrav,tterm, &
                ((dat(i,k,j),i=1,nprint),k=1,6), &
                ((dat(i,k,j),i=1,nprint),k=8,ncol), &
                (dat(i,9,j), i=nprint+1, nprint+nptmass), &
                ((dat(i,k,j), i=nprint+1, nprint+nptmass),k=1,6), &
                ((dat(i,k,j), i=nprint+1, nprint+nptmass),k=21,23)

                if (ierr < 0) print "(a)",'*** WARNING: END OF FILE DURING READ ***'
                if (ierr > 0) print "(a)",'*** WARNING: ERRORS DURING READ ***'

           endif
        endif
     else
     !
     !--hydro minidumps
     !
        if (minidump) then
           if (doubleprec) then
              allocate(datdb(maxpart,7),stat=ierr)
              if (ierr /= 0) then
                 print*,"(a)",'*** error allocating memory for double conversion ***'
                 return
              else
                 datdb = 0.
              endif
              read(15,end=55,iostat=ierr) timedb,nprint,nptmass, &
              ((datdb(i,k), i=1, nprint),k=1,6), &
              (datdb(i,7), i=nprint+1, nprint+nptmass), &
              ((datdb(i,k), i=nprint+1, nprint+nptmass),k=1,3)

              if (ierr < 0) print "(a)",'*** WARNING: END OF FILE DURING READ ***'
              if (ierr > 0) print "(a)",'*** WARNING: ERRORS DURING READ ***'
              dat(:,1:ncol,j) = real(datdb(:,1:ncol))
              time(j) = real(timedb)
           else
              read(15,end=55,iostat=ierr) time(j),nprint,nptmass, &
              ((dat(i,k,j), i=1, nprint),k=1,6), &
              (dat(i,7,j), i=nprint+1, nprint+nptmass), &
              ((dat(i,k,j), i=nprint+1, nprint+nptmass),k=1,3)
           endif
           !
           !--because masses are not dumped, we need to reconstruct them
           !  from density and h (only strictly true for grad h code)
           !
           dat(1:nprint,7,j) = dat(1:nprint,4,j)**3*dat(1:nprint,5,j)*dhfact3
           print "(a,f5.2,a)", &
            ' WARNING: setting particle masses assuming h = ',hfact,'*(m/rho)^(1/3)'
        else
     !
     !--hydro full dumps
     !
           dat(:,:,j) = 0. ! because ptmasses don't have all quantities

           if (doubleprec) then
              allocate(datdb(maxpart,ncol),stat=ierr)
              if (ierr /= 0) then
                 print*,"(a)",'*** error allocating memory for double conversion ***'
                 return
              else
                 datdb = 0.
              endif
              read(15,iostat=ierr) nprint,rstardb,mstardb,n1,n2, &
                nptmass,timedb,(datdb(i,7),i=1,nprint), &
                escapdb,tkindb,tgravdb,ttermdb, &
                ((datdb(i,k),i=1,nprint),k=1,6), &
                ((datdb(i,k),i=1,nprint),k=8,ncol), &
                (datdb(i,9), i=nprint+1, nprint+nptmass), &
                ((datdb(i,k), i=nprint+1, nprint+nptmass),k=1,6), &
                ((datdb(i,k), i=nprint+1, nprint+nptmass),k=13,15)

                if (ierr < 0) print "(a)",'*** WARNING: END OF FILE DURING READ ***'
                if (ierr > 0) print "(a)",'*** WARNING: ERRORS DURING READ ***'
                dat(:,1:ncol,j) = real(datdb(:,1:ncol))
                time(j) = real(timedb)

           else
              read(15,iostat=ierr) nprint,rstar,mstar,n1,n2, &
                nptmass,time(j),(dat(i,7,j),i=1,nprint), &
                escap,tkin,tgrav,tterm, &
                ((dat(i,k,j),i=1,nprint),k=1,6), &
                ((dat(i,k,j),i=1,nprint),k=8,ncol), &
                (dat(i,9,j), i=nprint+1, nprint+nptmass), &
                ((dat(i,k,j), i=nprint+1, nprint+nptmass),k=1,6), &
                ((dat(i,k,j), i=nprint+1, nprint+nptmass),k=13,15)

                if (ierr < 0) print "(a)",'*** WARNING: END OF FILE DURING READ ***'
                if (ierr > 0) print "(a)",'*** WARNING: ERRORS DURING READ ***'

           endif
        endif
     endif

     if (ierr /= 0 ) then
        print "(a)",'|*** ERROR READING TIMESTEP ***'
!           return
!        else
!           nstepsread = nstepsread + 1
     endif
     nstepsread = nstepsread + 1

     npartoftype(1,j) = nprint
     npartoftype(2,j) = nptmass
!!        print*,j,' time = ',time(j)
     gamma(j) = 1.666666666667
!
!--read abundances from abundance file
!
     if (iabunfileopen) then
        print "(a)",' ... reading abundances from '//trim(abunfile)//' ...'
        read(41,iostat=ierr1) nprint
        if (ierr1 /= 0) then
           print "(a)",' *** ERROR READING ABUNDANCE FILE ***'
        elseif (nprint.ne.npartoftype(1,j)) then
           print "(a)",' *** ERROR: npart in abundance file differs from full dump ***'
        else
           rewind(41)
           if (doubleprec) then
              read(41,iostat=ierr1) nprint,((datdb(i,k),k=1,max_spec+2),i=1,nprint),Etot_burn_cgsdb
              dat(:,ncol+1:ncolumns,j) = real(datdb(:,1:max_spec+2))
              print*,' Etot_burn (cgs) = ',Etot_burn_cgsdb
           else
              read(41,iostat=ierr1) nprint,((dat(i,k,j),k=1,max_spec+2),i=1,nprint),Etot_burn_cgs
              print*,' Etot_burn (cgs) = ',Etot_burn_cgs
           endif
           if (ierr1 < 0) then
              print "(a)",' *** END OF FILE REACHED IN ABUNDANCE FILE ***'
           elseif (ierr1 > 0) then
              print "(a)",' *** ERRORS DURING ABUNDANCE FILE READ ***'
           elseif (nprint.ne.npartoftype(1,j)) then
              print "(a)",' *** ERROR: npart in abundance file differs from full dump ***'
           endif
        endif
        close(unit=41)
     endif

     j = j + 1

     if (allocated(datdb)) deallocate(datdb)

  endif

55 continue
  !
  !--reached end of file
  !
  close(15)

  !
  !--reset centre of mass to zero
  !
  if (allocated(dat) .and. n2.eq.0 .and. lenvironment('RSPLASH_RESET_CM')) then
     if (minidump) then
        call reset_centre_of_mass(dat(1:nprint,1:3,j-1),dat(1:nprint,7,j-1),nprint)
     else ! full dumps ipmass = 9
        call reset_centre_of_mass(dat(1:nprint,1:3,j-1),dat(1:nprint,9,j-1),nprint)
     endif
  endif

  !
  !--transform velocities to corotating frame
  !
  if (.not.minidump .and. allocated(dat) .and. lenvironment('RSPLASH_COROTATING')) then
     print*,'TRANSFORMING VELOCITIES TO CORORATING FRAME'
     call set_corotating_vels(dat(1:nprint,9,j-1),dat(1:nprint,4:5,j-1),n1,nprint)
  endif

  if (allocated(npartoftype)) then
     print*,'>> end of dump file: nsteps =',j-1,'ntot = ', &
        sum(npartoftype(:,j-1)),'nptmass=',npartoftype(2,j-1)
  endif
return

contains

!
!--reset centre of mass to zero
!
 subroutine reset_centre_of_mass(xyz,pmass,npart)
  implicit none
  integer, intent(in) :: npart
  real, dimension(npart,3), intent(inout) :: xyz
  real, dimension(npart) :: pmass
  real :: masstot,xcm,ycm,zcm

  !
  !--get centre of mass
  !
  masstot = SUM(pmass(1:npart))
  xcm = SUM(pmass(1:npart)*xyz(1:npart,1))/masstot
  ycm = SUM(pmass(1:npart)*xyz(1:npart,2))/masstot
  zcm = SUM(pmass(1:npart)*xyz(1:npart,3))/masstot

  print*,' centre of mass is at ',xcm,ycm,zcm
  print*,' resetting to zero...'
  xyz(1:npart,1) = xyz(1:npart,1) - xcm
  xyz(1:npart,2) = xyz(1:npart,2) - ycm
  xyz(1:npart,3) = xyz(1:npart,3) - zcm

  return
 end subroutine reset_centre_of_mass

!
!--adjust velocities to corotating frame
!
 subroutine set_corotating_vels(pmass,vxy,n1,npart)
  implicit none
  integer, intent(in) :: n1,npart
  real, dimension(npart,2), intent(inout) :: vxy !, xy
  real, dimension(npart) :: pmass
  real :: mass1,mass2 !,xcm1,ycm1,xcm2,ycm2
  real :: vxcm1,vycm1,vxcm2,vycm2

  !
  !--get centre of mass of star 1 and star 2
  !
  mass1 = SUM(pmass(1:n1))
!  xcm1 = SUM(pmass(1:n1)*xy(1:n1,1))/mass1
!  ycm1 = SUM(pmass(1:n1)*xy(1:n1,2))/mass1

  mass2 = SUM(pmass(n1+1:npart))
!  xcm2 = SUM(pmass(n1+1:npart)*xy(n1+1:npart,1))/mass2
!  ycm2 = SUM(pmass(n1+1:npart)*xy(n1+1:npart,2))/mass2
  !
  !--work out centre of mass velocities for each star
  !
  vxcm1 = SUM(pmass(1:n1)*vxy(1:n1,1))/mass1
  vycm1 = SUM(pmass(1:n1)*vxy(1:n1,2))/mass1

  vxcm2 = SUM(pmass(n1+1:npart)*vxy(n1+1:npart,1))/mass2
  vycm2 = SUM(pmass(n1+1:npart)*vxy(n1+1:npart,2))/mass2

  !
  !--subtract centre of mass velocities appropriately
  !
  vxy(1:n1,1) = vxy(1:n1,1) - vxcm1
  vxy(1:n1,2) = vxy(1:n1,2) - vycm1
  vxy(n1+1:npart,1) = vxy(n1+1:npart,1) - vxcm2
  vxy(n1+1:npart,2) = vxy(n1+1:npart,2) - vycm2

 end subroutine set_corotating_vels

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use filenames, only:rootname
  use labels, only:label,unitslabel,labelvec,labeltype,iamvec,&
              ix,ivx,ih,irho,iutherm,ipmass,iBfirst,idivB
  use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings,iformat
  use geometry, only:labelcoord
  use settings_units, only:units
  implicit none
  integer :: i
  logical :: minidump
  real :: udistcm,udistkm,utime,umass,uvelkms

  minidump = .false.
  if (index(rootname(1),'minidump').ne.0) minidump = .true.

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif

  do i=1,ndim
     ix(i) = i
  enddo
  if (minidump) then
     ivx = 0
     ih = 4        !  smoothing length
     irho = 5     ! location of rho in data array
     iutherm = 0  !  thermal energy
     ipmass = 7  !  particle mass
     label(6) = 'T'
     if (ncolumns.gt.7) then
        iBfirst = 8
        idivB = 11
     endif
     if (ncolumns.gt.11) then
        label(12) = 'grad h'
        label(13) = 'grad soft'
        label(14) = 'dsoft'
     endif
!
!--set transformation factors between code units/real units
!
     udistkm = 1.5  ! km
     udistcm = 1.5e5
     utime = 5.0415e-6
     umass = 1.99e33
     units(1:3) = udistkm
     unitslabel(1:3) = ' [km]'
     units(ih) = udistkm
     unitslabel(ih) = ' [km]'
     units(ipmass) = umass
     unitslabel(ipmass) = ' [g]'
     units(irho) = umass/udistcm**3
     unitslabel(irho) = ' [g/cm\u3\d]'
     if (iBfirst.gt.0) then
        units(iBfirst:iBfirst+ndimV-1) = 8.0988e14
        unitslabel(iBfirst:iBfirst+ndimV-1) = ' [G]'
        units(idivB) = units(iBfirst)/udistcm
        unitslabel(idivB) = ' [G/cm]'
     endif

  else !--full dump
     ivx = ndim+1
     ih = 7
     iutherm = 8
     ipmass = 9
     irho = 10
     label(11) = 'temperature [ MeV ]'
     label(12) = 'electron fraction (y\de\u)'

     if (iformat.eq.2) then ! MHD full dump
        iBfirst = 13
        label(16) = 'psi'
        idivB = 17
        iamvec(18:20) = 18
        labelvec(18:20) = 'force'
        do i=1,ndimV
           label(18+i-1) = trim(labelvec(18))//'\d'//labelcoord(i,1)
        enddo
        label(21) = 'Euler alpha'
        label(22) = 'Euler beta'
        label(23) = 'Bevol\dz'
        label(24) = 'grad h'
        label(25) = 'grad soft'
        label(26) = 'av  '
        label(27) = 'avB'
!
!--set transformation factors between code units/real units
!
        udistkm = 1.5  ! km
        udistcm = 1.5e5
        utime = 5.0415e-6
        umass = 1.99e33

        units(iBfirst:iBfirst+ndimV-1) = 8.0988e14
        unitslabel(iBfirst:iBfirst+ndimV-1) = ' [G]'
        units(idivB) = units(iBfirst)/udistcm
        unitslabel(idivB) = ' [G/cm]'

        units(1:3) = udistkm
        unitslabel(1:3) = ' [km]'
        units(4:6) = 1.0
        unitslabel(4:6) = '/c'
        units(7) = udistkm
        unitslabel(7) = ' [km]'
        units(8) = (udistcm/utime)**2
        unitslabel(8) = ' [erg/g]'
        units(9) = umass
        unitslabel(9) = ' [g]'
        units(10) = umass/udistcm**3
        unitslabel(10) = ' [g/cm\u3\d]'

     else
        iamvec(13:15) = 13
        labelvec(13:15) = 'force'
        do i=1,ndimV
           label(13+i-1) = trim(labelvec(13))//'\d'//labelcoord(i,1)
        enddo
        label(16) = 'dgrav'
        if (ncolumns.gt.16) then
           label(11) = 'temperature [ 10\u6\dK ]'
           do i=17,ncolumns
              write(label(i),"('species ',i2)") i-16
           enddo
           if (ncolumns.ge.25) then
              label(17) = 'He'
              label(18) = 'C'
              label(19) = 'O'
              label(20) = 'Ne'
              label(21) = 'Mg'
              label(22) = 'Si'
              label(23) = 'Fe'
              label(24) = 'mean A'
              label(25) = 'mean Z'
           endif
        endif

        udistcm = 1.0e9
        utime = 2.7443
        umass = 1.99e33
        uvelkms = (udistcm/utime)/1e5

        units(1:3) = 1.0 !!udistcm
        unitslabel(1:3) = ' [10\u9\d cm]'
        units(4:6) = uvelkms
        unitslabel(4:6) = ' [km/s]'
        units(ih) = units(1)
        unitslabel(ih) = unitslabel(1)
        units(8) = (udistcm/utime)**2
        unitslabel(8) = ' [erg/g]'
        units(9) = umass
        unitslabel(9) = ' [g]'
        units(10) = umass/udistcm**3
        unitslabel(10) = ' [g/cm\u3\d]'
     endif

  endif

  units(0) = utime*1000.
  unitslabel(0) = ' ms'

  if (ivx.ne.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
     enddo
  endif

  if (iBfirst.ne.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
     enddo
     label(idivB) = 'div B'
  endif

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = '\gr'
  if (iutherm.gt.0) label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'
  !
  !--set labels for each particle type
  !
  ntypes = 2 !!maxparttypes
  labeltype(1) = 'gas'
  labeltype(2) = 'point mass'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.

!-----------------------------------------------------------

  return
end subroutine set_labels
