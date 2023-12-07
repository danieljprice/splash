!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE STARSMASHER CODE
!
!-------------------------------------------------------------------------

module readdata_starsmasher
 implicit none

 public :: read_data_starsmasher, set_labels_starsmasher

 private
 integer :: frame = -1  ! default value of -1

contains


subroutine read_data_starsmasher(rootname,istepstart,ipos,nstepsread)
  use particle_data, only:dat,iamtype,npartoftype,time,gamma,headervals,&
                          maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread
  use mem_allocation, only:alloc
  use labels,       only:ih,irho,headertags
  use system_utils, only:renvironment,lenvironment
  use prompting,    only:prompt
  integer, intent(in) :: istepstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+10) :: datfile
  integer :: i,j,ierr
  integer :: ncolstep,npart_max,nstep_max,ntoti
  logical :: iexist,reallocate
  real(doub_prec) :: timetemp,dummy
  integer :: ntot, nnopt, nout, nit, nav, ngr, nrelax
  real(doub_prec) :: hmin, hmax, sep0, tf, dtout, alpha, beta, eta2, trelax, dt, omega2
  real(doub_prec) :: dx, dy, dz, dm, dh, drho, dvx, dvy, dvz, dudot
  real(doub_prec) :: duth, dmmu
  real(doub_prec) :: theta
  real(doub_prec) :: gram, sec, cm, kelvin, erg, boltz
  parameter(gram=1.d0,sec=1.d0,cm=1.d0,kelvin=1.d0)
  parameter(erg=gram*cm**2/sec**2)
  parameter(boltz=1.380658d-16*erg/kelvin)
  real(doub_prec) :: munit, runit, gravconst
  parameter(munit=1.9891d33,runit=6.9599d10)
  parameter(gravconst = 6.67390d-08)

  nstepsread = 0
  npart_max = maxpart

  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** '
     return
  endif
!
!--check if first data file exists
!
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(datfile)//': file not found ***'
     return
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim = 3
  ndimV = 3
!
!--read data from snapshots
!
  i = istepstart

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(11,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)", '*** ERROR OPENING FILE ***'
     return
  endif
  !
  !--read header for this timestep
  !

  read(11, iostat = ierr) &
       ntot, nnopt, hmin, hmax, sep0, tf, dtout, nout, nit, timetemp, &
       nav, alpha, beta, eta2, ngr, nrelax, trelax, dt, omega2

  if (omega2 > 0.d0 .and. frame == -1) then
     ! frame can equal -1 only the first time through, so this question
     ! will get asked (at most) only once
     print *, 'Period of rotating frame=',8*atan(1.d0)/omega2**0.5d0
     frame = 0
     call prompt('Do you want movie in the inertial frame? (1=yes)',frame,0,1)
  endif
  if (frame.eq.-1) frame=0

  if (ierr /= 0) then
     print "(a)", '*** ERROR READING TIMESTEP HEADER ***'
     return
  endif

  iformat = 0
  ncolstep = 13 ! 3 x pos, 3 x vel, pmass, rho, utherm, mean_mu, h, du/dt, temperature
  ncolumns = ncolstep

  irho = 8
  ih   = ncolstep

  ntoti = ntot    !int(sum(npartoftypei(1:6)))
  print*,'time             : ',timetemp
  print*,'N_total          : ',ntoti
  print*,'N data columns   : ',ncolstep

  !
  !--if successfully read header, increment the nstepsread counter
  !
  nstepsread = nstepsread + 1
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*ntoti)
     else
        ! if first time, save on memory
        npart_max = int(ntoti)
     endif
  endif
  if (i.ge.maxstep .and. i.ne.1) then
     nstep_max = i + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif

  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol),mixedtypes=.true.)
  endif
  !
  !--copy header into header arrays
  !
  npartoftype(:,i) = 0
  !
  !--set time to be used in the legend
  !
  time(i) = real(timetemp)
  !
  !--store other quantities from the file header
  !
  headertags(1:18) = (/'ntot  ','nnopt ','hmin  ','hmax  ',&
                       'sep0  ','tf    ','dtout ','nout  ',&
                       'nit   ','nav   ','alpha ','beta  ',&
                       'eta2  ','ngr   ','nrelax','trelax',&
                       'dt    ','omega2'/)
  headervals(1:18,i) = (/real(ntot),real(nnopt),real(hmin),real(hmax),&
                         real(sep0),real(tf),real(dtout),real(nout),&
                         real(nit),real(nav),real(alpha),real(beta), &
                         real(eta2),real(ngr),real(nrelax),real(trelax),&
                         real(dt),real(omega2)/)
  !
  !--read particle data
  !
  if (ntoti.gt.0) then
     !
     !--read particles' data
     !
     do j=1,ntot

        read(11, iostat = ierr) &
             dx, dy, dz,                           &   ! positions
             dm,                                   &   ! mass
             dh, drho,                             &   ! h & rho
             dvx, dvy, dvz,                        &   ! velocities
             dummy, dummy, dummy,                  &   ! vxdot, vydot, vzdot
             duth,                                 &   ! uthermal
             dudot,                                &   ! udot
             dummy, & ! dummy, dummy, dummy,           &   ! gx, gy, gz, gpot
             dmmu,                                 &   ! mean_mu
             dummy, dummy!, dummy, dummy                ! aa, bb, cc, divv

        if(frame.eq.1) then
           theta=omega2**0.5d0*timetemp
! rotation is counterclockwise:
           dat(j,1,i)=dx*cos(theta)-dy*sin(theta)
           dat(j,2,i)=dx*sin(theta)+dy*cos(theta)
           dat(j,4,i)=dvx*cos(theta)-dvy*sin(theta)
           dat(j,5,i)=dvx*sin(theta)+dvy*cos(theta)
        else
           dat(j,1,i)  = dx
           dat(j,2,i)  = dy
           dat(j,4,i)  = dvx
           dat(j,5,i)  = dvy
        endif

        dat(j,3,i)  = dz
        dat(j,6,i)  = dvz

        dat(j,7,i)  = dm

        dat(j,8,i)  = drho
        dat(j,9,i)  = duth
        dat(j,10,i) = dmmu
        dat(j,11,i) = dh
        dat(j,12,i) = dudot
        dat(j,13,i) = duth*gravconst*munit/runit/(1.5d0*boltz/dmmu)

        if (abs(duth) > tiny(duth)) then
           iamtype(j,i) = 1
           npartoftype(1,i) = npartoftype(1,i)+1
        else
           iamtype(j,i) = 2
           npartoftype(2,i) = npartoftype(2,i)+1
        endif

     end do
!    pos = 1,2,3
!    vel = 4,5,6
!    mass = 7
!    rho, u, mean_mu, h = 8, 9, 10, 11

  else
     ntoti = 1
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif

!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read
!
  if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--close data file and return
!
  close(unit=11)

  if (nstepsread.gt.0) then
     print*,'>> last step ntot =',sum(npartoftype(:,istepstart+nstepsread-1))
  endif
  return

end subroutine read_data_starsmasher

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels_starsmasher
  use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,ih,&
                          irho,ipr,iutherm,make_vector_label
  use params
  use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings
  use geometry,      only:labelcoord
  use system_utils,  only:renvironment
  integer :: i

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
  ivx = 4
  ipmass = 7
  irho = 8        ! location of rho in data array
  ipr = 0
  iutherm = 9     !  thermal energy
  ih  = 11

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = 'density'
  label(iutherm) = 'specific internal energy u'
  label(ih) = 'h'
  label(10) = 'mean_mu'
  label(ipmass) = 'particle mass'
  label(12) = 'du/dt'
  label(13) = 'temperature'
  !
  !--set labels for vector quantities
  !
  call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
  !
  !--set labels for each particle type
  !
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'compact object'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.

end subroutine set_labels_starsmasher

end module readdata_starsmasher
