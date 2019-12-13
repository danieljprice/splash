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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! THIS VERSION FOR SARAH MADDISON+MARK HUTCHISON'S DUSTY-SPH CODE
! -> Now automatically handles single/double precision
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
! npartoftype(maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,ipos,nstepsread)
 use particle_data,  only:npartoftype,time,gamma,dat,maxpart,maxstep,maxcol,iamtype
 use params
 use filenames,      only:nfiles
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,iverbose,debugmode
 use mem_allocation, only:alloc
 use system_utils,   only:lenvironment
 implicit none
 integer,          intent(in)  :: indexstart,ipos
 integer,          intent(out) :: nstepsread
 character(len=*), intent(in)  :: rootname
 character(len=len(rootname)+4) :: datfile
 integer :: i,icol,ierr,iunit,ilen,j,ilast,k
 integer :: ihead1,ihead2,ihead3,numdata ! # of items on header lines
 integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
 integer :: ncolstep,np,istep,nstepsinfile
 logical :: reallocate,finished

 integer :: norigin,ncr,istart,iout,nmlmax,index,ratio
 integer :: dimsw,fluidsw,itrace
 integer, dimension(:), allocatable :: head1
 integer, dimension(:), allocatable :: head2
 real(doub_prec), dimension(:), allocatable :: head3
 real(doub_prec) :: timei,dt,sopt0
 real(doub_prec) :: gammai
 real(doub_prec), dimension(:), allocatable :: dattemp
 integer, dimension(:), allocatable :: iam
 integer, dimension(:), allocatable :: iwas
 real :: dum

 iunit = 11 ! file unit number
 ndim_max = 1
 ndimV_max = 1
 nstepsread = 0
 if (rootname(1:1) /= ' ') then
    datfile = trim(rootname)
    !print*,'rootname = ',rootname
 else
    print*,' **** no data read **** '
    return
 endif

 if (iverbose >= 1) print "(1x,a)",'reading Maddison/Hutchison format'
 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
!
!--open data file and read data
!
 open(unit=iunit,iostat=ierr,file=datfile,status='old',form='unformatted')
 if (ierr /= 0) then
    print*,' *** Error opening '//trim(datfile)//' ***'
    return
 endif
!
!--read first header line
!
 read(iunit,iostat=ierr,end=80) ihead1,ihead2,ihead3,ncolstep,np

 allocate(head1(ihead1))
 allocate(head2(ihead2))
 allocate(head3(ihead3))

 read(iunit,iostat=ierr,end=80) head1(1:ihead1)
 read(iunit,iostat=ierr,end=80) head2(1:ihead2)
 read(iunit,iostat=ierr,end=80) head3(1:ihead3)

 ! Header one variables
 norigin = head1(1)
 ncr     = head1(2)
 istart  = head1(3)
 iout    = head1(4)
 nmlmax  = head1(5)
 index   = head1(6)
 ratio   = head1(7)
 itrace  = head1(8)

 ! Header two variables
 ndim       = abs(head2(1))
 fluidsw    = head2(2)
 !velsw      = head2(3)
 !dragsw     = head2(4)
 !Bkernsw    = head2(5)
 !Dkernsw    = head2(6)
 !freesw     = head2(7)
 !stopgas    = head2(8)
 !coolsw     = head2(9)
 !interpsw   = head2(10)
 !voidsw     = head2(11)
 !phasesw    = head2(12)
 !gravsw     = head2(13)
 !photosw    = head2(14)
 !eossw      = head2(15)
 !masssw     = head2(16)
 !usenumdens = head2(17)
 !varhsw     = head2(18)
 !tracesw    = head2(19)
 !timesw     = head2(20)
 !growsw     = head2(21)
 !shocksw    = head2(22)
 !wavesw     = head2(23)
 !plotsw     = head2(24)
 !disc1dsw   = head2(25)

 ! Header three variables
 dt        = head3(1)
 timei     = head3(2)
 sopt0     = head3(3)
 gammai    = head3(4)

 deallocate(head1)
 deallocate(head2)
 deallocate(head3)

 ncolstep  = ncolstep-2  ! minus 2 because iwas and iam read individually
 ndimV     = 3  ! always have 3 velocity components written to file

 print "(a,i2,a,f8.4)",' ncolumns: ',ncolstep,' gamma: ',gammai

 !
 !--check for basic errors in first line
 !

 if (ierr /= 0 .or. ncr < 0 .or. istart < 0 &
     .or. iout  < 0 .or. nmlmax < 0 .or. index < 0 .or. ratio < 0) then

    print "(a)",' *** Error reading header ***'
    print*,' norigin = ',norigin,' ncr = ',ncr,' istart =',istart,' iout = ',iout
    print*,' nmlmax = ',nmlmax,' index = ',index,' ratio =',ratio
    close(iunit)
    return
 endif
 !
 !--Check for errors
 !
 if (ierr /= 0 .or. np < 0 .or. np > 1.e9 .or. itrace > np) then
    print*,'n = ',np,' dt = ',dt,' time = ',timei,' i = ',itrace
    print*,'*** error reading timestep header ***'
    close(iunit)
    return
 endif

 !
 !--check for errors in 3rd line
 !
 if (ndim > 3 .or. ndimV > 3) then
    print*,'*** error in header: ndim or ndimV in file > 3'
    ndim  = 3
    ndimV = 3
    close(iunit)
    return
 endif

 nstepsinfile = 1 ! nmlmax/iout
 nstep_max = 1    ! max(nstepsinfile,maxstep)
 nstepsread = 0
 npart_max = maxpart
 ncol_max  = ncolstep
!
!--read first step
!
 over_steps: do i = indexstart,indexstart + nstepsinfile - 1
!
!--allocate memory for data arrays
!
    nstep_max = nstepsinfile !max(nstep_max,nfiles,maxstep,indexstart)
    npart_max = max(np,maxpart,norigin)

    if (.not.allocated(dat) .or. np > maxpart  &
          .or. nstep_max > maxstep .or. ncol_max > maxcol) then
       call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=.true.)
    endif
    !
    !--now that memory is allocated, put header quantities -> splash quantities
    !
    time(i) = timei
    gamma(i) = gammai
    npartoftype(1,i) = np
    if (iverbose >= 1) then
       print "(a,i5,a,f8.4,a,i8,a,f8.4)",' step:',i,' time:',time(i),' npart:',np,' dt:',dt
    else
       print "(a,i5,a,f8.4,a,i8,a,i8)",' step:',i,' time:',time(i),' npart:',np
    endif

    if (ncolstep /= ncol_max) then
       print*,'*** Warning number of columns not equal for timesteps'
       ncolumns = ncolstep
       if (iverbose >= 1) print*,'ncolumns = ',ncolumns,ncol_max
       if (ncolumns > ncol_max) ncol_max = ncolumns
    endif

    ncolumns = ncolstep
    nstepsread = nstepsread + 1

    !
    !--read data for this timestep
    !
    allocate(iwas(np))
    read(iunit,iostat=ierr,end=80)  iwas(1:np)

    if ( any(iwas < 1).or.any(iwas > norigin) ) then
       do j = 1,np
          iwas(j) = j
       enddo
    endif

    npartoftype(:,i) = 0
    allocate(iam(norigin))

    iam(:) = 3

    read(iunit,iostat=ierr,end=80) iam(iwas(:))

    do j=1,norigin
       select case(iam(j))
       case(1)       ! GAS
          npartoftype(1,i) = npartoftype(1,i) + 1
          iamtype(j,i) = 1_int1
       case(0)       ! DUST
          npartoftype(2,i) = npartoftype(2,i) + 1
          iamtype(j,i) = 2_int1
       case default  ! DELETED PARTICLES
          npartoftype(3,i) = npartoftype(3,i) + 1
          iamtype(j,i) = 3_int1
       end select
    enddo
    deallocate(iam)

    if (kind(dat) /= kind(dattemp)) then
       if (debugmode) print*,' converting kind from ',kind(dattemp),' to ',kind(dat)
       allocate(dattemp(norigin))

       !dattemp(:) = 0D0

       !--convert precision
       !do icol=1,ncolstep-1 ! all columns except h
       do icol=1,ncolstep ! all columns with h
          read(iunit,iostat=ierr,end=80) dattemp(iwas(:))
          dat(1:norigin,icol,i) = real(dattemp)
       enddo
       deallocate(dattemp)
    else
       !--read directly into dat array if data types are the same
       !do icol=1,ncolstep-1 ! all columns except h
       do icol=1,ncolstep ! all columns with h
          read(iunit,iostat=ierr,end=80) dat(iwas(:),icol,i)
       enddo
    endif

    deallocate(iwas)
 enddo over_steps

 close(unit=11)

 ilast = indexstart+nstepsinfile - 1
 ncolumns = ncol_max

 call set_labels

 if ( fluidsw < 0 .and. .not.lenvironment('NSPLASH_BARYCENTRIC') ) then
    call fake_twofluids
 endif

 if (npartoftype(2,ilast) > 0) then
    print*,' ngas = ',npartoftype(1,ilast),' ndust = ',npartoftype(2,ilast)
    if (npartoftype(3,ilast) > 0) print*,' nunknown = ',npartoftype(3,ilast)
 endif
 if (debugmode) print*,'DEBUG> Read steps ',indexstart,'->',indexstart + nstepsread - 1, &
       ' last step ntot = ',sum(npartoftype(:,indexstart+nstepsread-1))
 return

80 continue
 print*,' *** data file empty : no timesteps ***'
 return

contains

!--------------------------------------------
!--------------------------------------------
subroutine fake_twofluids
 use labels, only:idustfrac,irho,ix,ih,ipmass,ivx,ideltav
 implicit none
 integer :: ndust,jdust
 integer :: ntoti
 real    :: rhodust,rhogas,rhotot,dustfraci,pmassgas,pmassdust,pmassj
 real, dimension(ndimV) :: veli,vgas,vdust,deltav

 integer :: itemp
 integer :: ifx,ify,ifz
 integer :: ir_grain
 integer :: ics
 integer :: idhdt
 integer :: idepsdt

 !do i=1,ndim
 !   ix(i) = i             ! x,y,z positions
 !enddo
 !!--2D means x-z
 !if (ndim==2) ix(2) = 3 ! x,z positions if 2D
 !ivx       = 4            ! velocity (vector so it takes 3 rows)
 !irho      = 7            ! density
 !ipmass    = 8            ! mass
 !iutherm   = 9            ! thermal energy
 !idendt    = 10           ! time derivative of thermal energy
 !itemp     = 11           ! temperature
 !ifx       = 12           ! force in x direction
 !ify       = 13           ! force in y direction
 !ifz       = 14           ! force in z direction
 ir_grain  = 15            ! grain size
 !ics       = 16           ! sound speed
 !ih        = 17           ! smoothing length
 !idhdt     = 18           ! time derivative of h
 !idustfrac = 19           ! dust fraction
 !idepsdt   = 20           ! time derivative of dust fraction


 if (idustfrac > 0 .and. irho > 0) then
    do i=indexstart,indexstart+nstepsread-1
       ntoti = sum(npartoftype(:,i))
       if (.not.allocated(dat) .or. (ntoti + npartoftype(1,i)) > maxpart) then
          call alloc(ntoti + npartoftype(1,i),maxstep,maxcol,mixedtypes=.true.)
       endif
       ndust = 0
       !--zero the properties of newly created dust particles
       dat(ntoti+1:ntoti+npartoftype(1,i),:,i) = 0.
       do j=1,ntoti
          if (iamtype(j,i)==1) then
             ndust = ndust + 1 ! one dust particle for every gas particle
             rhotot  = dat(j,irho,i)
             dustfraci = dat(j,idustfrac,i)
             rhogas  = rhotot*(1. - dustfraci)
             rhodust = rhotot*dustfraci
             !--replace global properties with gas-only stuff
             dat(j,irho,i) = rhogas
             !--copy x, smoothing length onto dust particle
             jdust = ntoti + ndust
             !--fix dust fraction for viewing purposes
             !dat(jdust,idustfrac,i) = dustfraci
             !dat(j,idustfrac,i) = 1.-dustfraci
             !--fix grain size to be for dust only
             dat(jdust,ir_grain,i) = dat(j,ir_grain,i)
             dat(j,ir_grain,i) = 0.

             !--fill in dust properties
             if (ndim > 0) dat(jdust,ix(1:ndim),i) = dat(j,ix(1:ndim),i)
             if (ih > 0)   dat(jdust,ih,i)         = dat(j,ih,i)
             if (irho > 0) dat(jdust,irho,i)       = rhodust
             iamtype(ntoti + ndust,i) = 2

             !--particle masses
             if (ipmass > 0) then
                pmassj    = dat(j,ipmass,i)
                pmassgas  = pmassj*(1. - dustfraci)
                pmassdust = pmassj*dustfraci
                dat(j,ipmass,i)     = pmassgas
                dat(jdust,ipmass,i) = pmassdust
             endif

             !--velocities
             if (ideltav > 0 .and. ivx > 0 .and. ndimV > 0) then
                veli(:)   = dat(j,ivx:ivx+ndimV-1,i)
                deltav(:) = dat(j,ideltav:ideltav+ndimV-1,i)
                if ( rhodust < 1.e-30 ) then
                   vgas(:)   = veli(:)
                   !vdust(:)  = 0.
                   vdust(:)  = 1.e10  ! Dirty way to clean up axis
                else
                   vgas(:)   = veli(:) - rhodust/rhotot*deltav(:)
                   vdust(:)  = veli(:) + rhogas/rhotot*deltav(:)
                endif
                dat(j    ,ivx:ivx+ndimV-1,i) = vgas(:)
                dat(jdust,ivx:ivx+ndimV-1,i) = vdust(:)
             endif
          endif
       enddo
       if (iverbose >= 1) then
          print "(a,i10,a)",' Creating ',ndust,' fictional dust particles...'
          print "(a)",' (set NSPLASH_BARYCENTRIC=yes to plot barycentric values)'
       endif
       npartoftype(2,i) = npartoftype(2,i) + ndust
    enddo
 else
    print "(a)",' ERROR: could not locate dust-to-gas ratio and/or density'
 endif

end subroutine fake_twofluids
!--------------------------------------------
!--------------------------------------------
end subroutine read_data


!------------------------------------------------------------
! set labels for each column of data
!------------------------------------------------------------

subroutine set_labels
 use labels, only:ix,ivx,ih,irho,iutherm,ipmass,    &
                  iamvec,labelvec,label,labeltype,  &
                  idustfrac,ideltav
 use params
 use settings_data, only:ndim,ndimV,iformat,ntypes, &
                    UseTypeInRenderings
 use geometry, only:labelcoord
 implicit none
 integer :: i

 integer :: idendt,itemp
 integer :: ifx,ify,ifz
 integer :: ir_grain
 integer :: ics
 integer :: idhdt
 integer :: idepsdt
 integer :: iddeltavdt
 integer :: itcour
 integer :: itstop
 integer :: itdiff

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
    return
 endif

 do i=1,ndim
    ix(i) = i              ! x,y,z positions
 enddo
 !--2D means x-z
 if (ndim==2) ix(2) = 3  ! x,z positions if 2D
 ivx        = 4            ! velocity (vector so it takes 3 rows)
 irho       = 7            ! density
 ipmass     = 8            ! mass
 iutherm    = 9            ! thermal energy
 idendt     = 10           ! time derivative of thermal energy
 itemp      = 11           ! temperature
 ifx        = 12           ! force in x direction
 ify        = 13           ! force in y direction
 ifz        = 14           ! force in z direction
 ir_grain   = 15           ! grain size
 ics        = 16           ! sound speed
 ih         = 17           ! smoothing length
 idhdt      = 18           ! time derivative of h
 idustfrac  = 19           ! dust fraction
 idepsdt    = 20           ! time derivative of dust fraction

 ideltav    = 21
 iddeltavdt = 24

 itcour     = 27
 itstop     = 28
 itdiff     = 29


!!TESTING
! itemp     = 10           ! temperature
! ifx       = 11           ! force in x direction
! ify       = 12           ! force in y direction
! ifz       = 13           ! force in z direction
! ir_grain  = 14           ! grain size
! ics       = 15           ! sound speed
! ih        = 16           ! smoothing length
! idhdt     = 17           ! time derivative of h
! idustfrac = 18           ! dust fraction
! idepsdt   = 19           ! time derivative of dust fraction
!
! ideltav   = 20
!!TESTING


 iamvec(ifx:ifz) = ifx
 iamvec(ivx:ivx+ndimV-1)   = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'
 labelvec(ifx:ifz)         = 'f'

 label(1:3)                = labelcoord(1:3,1)
 label(ivx:ivx+ndimV-1)    = 'v_'//labelcoord(1:3,1)
 label(irho)               = '\rho'
 label(ipmass)             = 'm_{particle}'
 label(iutherm)            = 'u'
 label(idendt)             = 'du/dt'
 label(itemp)              = 'Temp'
 label(ifx:ifz)            = 'f_'//labelcoord(1:3,1)
 label(ir_grain)           = 's_{grain}'
 label(ics)                = 'cs'
 label(ih)                 = 'h'
 label(idhdt)              = 'dh/dt'
 label(idustfrac)          = '\epsilon'
 label(idepsdt)            = 'd\epsilon/dt'
 label(ideltav)            = '\Delta v_x'
 label(ideltav+1)          = '\Delta v_y'
 label(ideltav+2)          = '\Delta v_z'
 label(iddeltavdt)         = 'd\Delta v_x/dt'
 label(iddeltavdt+1)       = 'd\Delta v_y/dt'
 label(iddeltavdt+2)       = 'd\Delta v_z/dt'

 label(itcour)             = 't_{cour}'
 label(itstop)             = 't_{stop}'
 label(itdiff)             = 't_{diff}'

!
!--set labels for each type of particles
!
 ntypes = 3
 labeltype(1) = 'gas'
 labeltype(2) = 'dust'
 labeltype(3) = 'unknown'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .true.
 UseTypeInRenderings(3) = .false.

!-----------------------------------------------------------

 return
end subroutine set_labels
