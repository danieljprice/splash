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
  integer :: i,icol,ierr,iunit,ilen,j,ilast
  integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
  integer :: ncolstep,np,istep,nstepsinfile
  logical :: reallocate,finished

  integer :: norigin,ncr,istart,iout,nmlmax,index,ratio
  integer :: dimsw,nfluid,ipart
  real(sing_prec) :: h,xlen,yzlen,xfact,hrej,rejfac
  real(sing_prec) :: gammai,alpha,beta,pspred
  real(sing_prec) :: gconst,zeta,amfac,frac,tin,vinit
  real(sing_prec) :: Cd,psep,kdrag,tcoeff,Rd
  real(sing_prec) :: a1,a2,a3,vdlfrac,vgrfrac,dustden,gasden,D2G_ratio
  real(sing_prec) :: csmax2,omega
  real(sing_prec) :: timei,dt
  real(sing_prec), dimension(:), allocatable :: dattemp
  integer, dimension(:), allocatable :: iam
  real :: dum

  iunit = 11 ! file unit number
  ndim_max = 1
  ndimV_max = 1
  nstepsread = 0
  if (rootname(1:1).ne.' ') then
     datfile = trim(rootname)
     !print*,'rootname = ',rootname
  else
     print*,' **** no data read **** '
     return
  endif

  if (iverbose.ge.1) print "(1x,a)",'reading Maddison/Hutchison format'
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
  read(iunit,iostat=ierr,end=80) norigin,ncr,istart,iout,nmlmax,index,ratio
  read(iunit,iostat=ierr,end=80) h,xlen,yzlen,xfact,hrej,rejfac, &
      gammai,alpha,beta,pspred, &
      gconst,zeta,amfac,frac,tin,vinit, &
      Cd,psep,kdrag,tcoeff,Rd, &
      a1,a2,a3, &
      vdlfrac,vgrfrac, &
      dustden,gasden,D2G_ratio, &
      csmax2, omega
  read(iunit,iostat=ierr,end=80) ndim,nfluid
  print*,'dimsw = ',ndim,' nfluid = ',nfluid,' h/psep = ',h/psep

  ncolstep  = 14 ! number of columns in the file
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
  !--check for errors in 3rd line
  !
  if (ndim > 3 .or. ndimV > 3) then
     print*,'*** error in header: ndim or ndimV in file > 3'
     ndim  = 3
     ndimV = 3
     close(iunit)
     return
  endif

  nstepsinfile = nmlmax/iout
  nstep_max = max(nstepsinfile,maxstep)
  nstepsread = 0
  npart_max = maxpart
  ncol_max  = ncolstep
!
!--read first step
!
  over_steps: do i = indexstart,indexstart + nstepsinfile - 1
     !
     !--read header line for this timestep
     !
     read(iunit,iostat=ierr,end=80) np,dt,timei,ipart
     if (ierr /= 0 .or. np < 0 .or. np > 1.e9 .or. ipart > np) then
        print*,'n = ',np,' dt = ',dt,' time = ',timei,' i = ',ipart
        print*,'*** error reading timestep header ***'
        close(iunit)
        return
     endif
!
!--allocate memory for data arrays
!
     nstep_max = max(nstep_max,nfiles,maxstep,indexstart)
     npart_max = max(np,maxpart)
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
     if (iverbose.ge.1) then
        print "(a,i5,a,f8.4,a,i8,a,f8.4)",' step:',i,' time:',time(i),' npart:',np,' dt:',dt
     else
        print "(a,i5,a,f8.4,a,i8,a,i8)",' step:',i,' time:',time(i),' npart:',np
     endif

     if (ncolstep.ne.ncol_max) then
        print*,'*** Warning number of columns not equal for timesteps'
        ncolumns = ncolstep
        if (iverbose.ge.1) print*,'ncolumns = ',ncolumns,ncol_max
        if (ncolumns.gt.ncol_max) ncol_max = ncolumns
     endif

     ncolumns = ncolstep
     nstepsread = nstepsread + 1
     !
     !--read data for this timestep
     !
     if (kind(dat).ne.kind(dattemp)) then
        if (debugmode) print*,' converting kind from ',kind(dattemp),' to ',kind(dat)
        allocate(dattemp(np))
        !--convert precision
        do icol=1,ncolstep-1 ! all columns except h
           read(iunit,iostat=ierr,end=80) dattemp(1:np)
           dat(1:np,icol,i) = real(dattemp)
        enddo
        deallocate(dattemp)
     else
        !--read directly into dat array if data types are the same
        do icol=1,ncolstep-1 ! all columns except h
           read(iunit,iostat=ierr,end=80) dat(1:np,icol,i)
        enddo
     endif
     !--add column for the smoothing length
     dat(1:np,ncolstep,i) = h
     
     npartoftype(:,i) = 0
     allocate(iam(np))
     read(iunit,iostat=ierr,end=80) iam(1:np)
     do j=1,np
        select case(iam(j))
        case(1)
           npartoftype(1,i) = npartoftype(1,i) + 1
           iamtype(j,i) = 1_int1
        case(0)
           npartoftype(2,i) = npartoftype(2,i) + 1
           iamtype(j,i) = 2_int1
        case default  ! unknown
           npartoftype(3,i) = npartoftype(3,i) + 1
           iamtype(j,i) = 3_int1
        end select
     enddo
     deallocate(iam)
     read(iunit,iostat=ierr,end=80)  ! skip iwas

     if (np <= 0) then ! handle zero particle case just-in-case
        npartoftype(1,i) = 1
        dat(:,:,i) = 0.
     endif

   enddo over_steps

close(unit=11)

ilast = indexstart+nstepsinfile - 1
ncolumns = ncol_max

if (npartoftype(2,ilast).gt.0) then
   print*,' ngas = ',npartoftype(1,ilast),' ndust = ',npartoftype(2,ilast)
   if (npartoftype(3,ilast) > 0) print*,' nunknown = ',npartoftype(3,ilast)
endif
if (debugmode) print*,'DEBUG> Read steps ',indexstart,'->',indexstart + nstepsread - 1, &
       ' last step ntot = ',sum(npartoftype(:,indexstart+nstepsread-1))
return

80 continue
print*,' *** data file empty : no timesteps ***'
return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels, only:ix,ivx,ih,irho,iutherm,ipmass,&
                  iamvec,labelvec,label,labeltype
 use params
 use settings_data, only:ndim,ndimV,iformat,ntypes, &
                    UseTypeInRenderings
 use geometry, only:labelcoord
 implicit none
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
 !--2D means x-z
 if (ndim.eq.2) ix(2) = 3
 ivx = 4
 irho = 7                ! location of rho in data arra
 ipmass = 8
 iutherm = 9
 label(10) = 'Temp'
 iamvec(11:13) = 11
 labelvec(11:13) = 'f'  ! force (vector)
 ih = 14                ! smoothing length

 label(1:3) = labelcoord(1:3,1)
 !
 !--label vector quantities (e.g. velocity) appropriately
 !
 iamvec(ivx:ivx+ndimV-1) = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'

 label(irho) = 'density'
 label(iutherm) = 'u'
 label(ih) = 'h'
 label(ipmass) = 'particle mass'
!
!--set labels for each type of particles
!
 ntypes = 3
 labeltype(1) = 'gas'
 labeltype(2) = 'dust'
 labeltype(3) = 'unknown'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .false.
 UseTypeInRenderings(3) = .false.

!-----------------------------------------------------------

 return
end subroutine set_labels
