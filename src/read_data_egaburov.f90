!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE
!
! NOTE THAT THIS ONLY "OFFICIALLY" WORKS WITH THE PARALLEL CODE AS WE
! REQUIRE KNOWLEDGE OF THE PARTICLE SMOOTHING LENGTHS
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! GSPLASH_USE_Z if 'YES' uses redshift in the legend instead of time
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

module readdata_egaburov
 implicit none

 public :: read_data_egaburov, set_labels_egaburov

 private
contains

subroutine read_data_egaburov(rootname,istepstart,ipos,nstepsread)
 use particle_data, only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep
 use params
 use settings_data, only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread
 use mem_allocation, only:alloc
 use labels, only:ih,irho
 use system_utils, only:renvironment,lenvironment
 implicit none
 integer, intent(in) :: istepstart,ipos
 integer, intent(out) :: nstepsread
 character(len=*), intent(in) :: rootname
 character(len=len(rootname)+10) :: datfile
 integer :: i,j,ierr
 integer :: ncolstep,npart_max,nstep_max
 logical :: iexist,reallocate


 integer         :: proc, nread

 integer         :: myproc, nproc, npx, npy, npz;
 integer         :: global_n, local_n, ndim_data
 real(sing_prec) :: t_global, dt_global
 integer         :: iteration
 real(sing_prec) :: cfl_no, gamma_gas
 integer         :: periodic_flag
 real(sing_prec) :: xmin, ymin, zmin, xmax, ymax, zmax

 integer         :: idx
 real(sing_prec) :: posx, posy, posz, pvelx, pvely, pvelz
 real(sing_prec) :: dens, ethm, pres, pmag, vabs
 real(sing_prec) :: velx, vely, velz, Bx, By, Bz
 real(sing_prec) :: hsml, wght, Bpsi, divB, v1, v2, v3, v4



!!!!!!!!!!!!!!!!!!

 nstepsread = 0
 npart_max = maxpart

 if (len_trim(rootname) > 0) then
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
    stop
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
 open(11,iostat=ierr,file=datfile,status='old', form='unformatted')
 if (ierr /= 0) then
    print "(a)", '*** ERROR OPENING FILE ***'
    stop
 endif

 !
 !--read header for this timestep
 !
 read(11, iostat=ierr)  myproc, nproc, npx, npy, npz, &
       global_n, local_n, ndim_data, t_global, dt_global, iteration, &
       cfl_no, gamma_gas, periodic_flag, &
       xmin, ymin, zmin, xmax, ymax, zmax
 print *, global_n
 print *, local_n
 print *, nproc
 if (ierr /= 0) then
    print "(a)", '*** ERROR READING TIMESTEP HEADER ***'
    stop
 endif

!  t_global = t_global  - 0.13

 iformat = 0
 ncolstep = 30
 ncolumns = ncolstep

 print*,'nproc            : ',nproc
 print*,'npx, npy, npz    : ',npx, npy, npz
 print*,'time             : ',t_global
 print*,'gamma_gas        : ',gamma_gas
 print*,'N_total          : ',global_n
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

 if (global_n  >  maxpart) then
    reallocate = .true.
    if (maxpart > 0) then
       ! if we are reallocating, try not to do it again
       npart_max = int(1.1*global_n)
    else
       ! if first time, save on memory
       npart_max = int(global_n)
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
    call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol))
 endif

 npartoftype(:,i) = 0
 npartoftype(1,i) = global_n

 !--use this line for code time
 time(i) = real(t_global)

 !
 !--read particle data
 !
 nread = 0
 do proc = 1, nproc
    do j = 1, local_n
       read(11, iostat=ierr) &
             idx, &
             posx, posy, posz, pvelx, pvely, pvelz, &
             dens, ethm, pres, pmag,  vabs, &
             velx, vely, velz, Bx, By, Bz, &
             hsml, wght, Bpsi, divB, &
             v1, v2, v3, v4
!        if (pres > 10) then
!           print *, posx, posy, posz
!        endif
       if (ierr /= 0) then
          print *, '*** ERROR READING PARTICLE', i, 'PROC:', myproc
          stop
       endif
       dat(j + nread, 1, i) = posx
       dat(j + nread, 2, i) = posy
       dat(j + nread, 3, i) = posz
       dat(j + nread, 4, i) = velx
       dat(j + nread, 5, i) = vely
       dat(j + nread, 6, i) = velz

       dat(j + nread, 7,  i) = dens
       dat(j + nread, 8,  i) = pres
       dat(j + nread, 9,  i) = pmag
       dat(j + nread, 10, i) = vabs

       dat(j + nread, 11, i) = Bx
       dat(j + nread, 12, i) = By
       dat(j + nread, 13, i) = Bz
       if (pmag > 0) then
          dat(j + nread, 14, i) = abs(divB/sqrt(pmag*2.0))
       else
          dat(j + nread, 14, i) = 0.0;
       end  if

       dat(j + nread, 15, i) = Bpsi
       dat(j + nread, 16, i) = hsml*2
       dat(j + nread, 17, i) = wght

       dat(j + nread, 18, i) = pvelx
       dat(j + nread, 19, i) = pvely
       dat(j + nread, 20, i) = pvelz
       dat(j + nread, 21, i) = abs(divB)

       dat(j + nread, 22, i) = dens*wght
       dat(j + nread, 23, i) = pres/dens/(gamma_gas - 1.0)   ! uthermal
       dat(j + nread, 24, i) = sqrt(pmag*2.0)
       dat(j + nread, 25, i) = pres + pmag
       dat(j + nread, 26, i) = v1;
       dat(j + nread, 27, i) = v2;
       dat(j + nread, 28, i) = v3;
       dat(j + nread, 29, i) = v4;
       dat(j + nread, 30, i) = idx;

    enddo
    nread = nread + local_n;

    if (proc < nproc) then
       read(11, iostat=ierr) myproc, nproc, npx, npy, npz, &
             global_n, local_n, ndim_data, t_global, dt_global, iteration, &
             cfl_no, gamma_gas, periodic_flag, &
             xmin, ymin, zmin, xmax, ymax, zmax
       if (ierr /= 0) then
          print *, '*** ERROR READING TIMESTEP HEADER ***, proc=', proc
          stop
       endif
    endif
 enddo

 if (nread  /=  global_n) then
    print *, 'nread=     ', nread
    print *, 'global_n=  ', global_n
    print "(a)", ' *** SOMETHING WENT WRONG ***'
    stop
 endif

!!!!!!!!!!!!!!!!!!!!!
 gamma = gamma_gas


 irho = 7
 ih   = 16

!
!--set flag to indicate that only part of this file has been read
!
 if (.not.all(required(1:ncolstep))) ipartialread = .true.

!
!--close data file and return
!
 close(unit=11)

 if (nstepsread > 0) then
    print*,'>> last step ntot =',sum(npartoftype(:,istepstart+nstepsread-1))
 endif
 return

end subroutine read_data_egaburov

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels_egaburov
 use labels, only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,ih,irho,ipr,iutherm, idivb, iBfirst
 use params
 use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings
 use geometry, only:labelcoord
 use system_utils, only:renvironment
 implicit none
 integer :: i

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_egaburov ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_egaburov ***'
    return
 endif

 do i=1,ndim
    ix(i) = i
 enddo

 ivx = 4
 irho = 7
 iutherm = 23
 ipr =  8
 ih  = 16
 iBfirst = 11
 ipmass = 22
 idivb  = 14

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 label(irho) = 'density'
 label(iutherm) = 'cs'
 label(ih) = 'h'
 label(ipmass) = 'particle mass'
 label(ipr) = 'pressure'
 label(9)  = 'Pmag'
 label(10) = 'vabs'
 label(15) = "\gpsi"
 label(18) = 'pvelx'
 label(19) = 'pvely'
 label(20) = 'pvelz'
 label(11) = 'Bx'
 label(12) = 'By'
 label(13) = 'Bz'
 label(14) = 'divB'
 label(24) = '|B|'
 label(25) = 'Ptot'
 label(26) = 'scal0'
 label(27) = 'scal1'
 label(28) = 'scal2'
 label(29) = 'scal3'
 label(30) = 'idx'


 !
 !--set labels for vector quantities
 !
 iamvec(ivx:ivx+ndimV-1) = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'
 do i=1,ndimV
    label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
 enddo
 !--mag field
 if (iBfirst > 0) then
    iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
    labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
    do i=1,ndimV
       label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
    enddo
 endif
 !
 !--set labels for each particle type
 !
 ntypes = 1
 labeltype(1) = 'gas'
 UseTypeInRenderings(1) = .true.

 return
end subroutine set_labels_egaburov
end module readdata_egaburov
