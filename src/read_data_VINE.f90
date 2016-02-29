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
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM THE VINE CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! VINE_MHD or VSPLASH_MHD if 'YES' or 'TRUE', reads MHD dump files
! VINE_HFAC or VSPLASH_HFAC if 'YES' or 'TRUE', multiplies the
!  smoothing lengths read from the file by a factor of 2.8 for
!  compatibility with older VINE dump files.
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
module vineread
  implicit none
!
! These are the indices of various values saved in
! the header part of the dump file.
!
! SHOULD BE IDENTICAL TO THOSE DEFINED IN THE VINE io.F FILE...
!
  integer,parameter::id_iheadlen= 1, id_npart  = 2, id_npart_sph = 3
  integer,parameter::id_nstep   = 4, id_idump  = 5, id_makechkpnt= 6
  integer,parameter::id_idbg    = 7, id_itme   = 8, id_ibctype   = 9
  integer,parameter::id_indts   =10, id_iusebin=11, id_ipres     =12
  integer,parameter::id_ipdv    =13, id_ieos   =14, id_icool     =15
  integer,parameter::id_iheat   =16, id_ivisc  =17, id_ivtime    =18
  integer,parameter::id_ibals   =19, id_ishok  =20, id_igrav     =21
  integer,parameter::id_isoft   =22, id_imac   =23, id_ihts      =24
  integer,parameter::id_iexpand =25, id_ncdumps=26, id_maxclumpsize=27
  integer,parameter::id_intgrtr =28, id_ishift =29, id_ndim     =30
  integer,parameter::id_io_fmt  =31, id_iunits =32, id_maxbunchsize=33
  integer,parameter::id_npoim   =34, id_ndt_ana=35, id_maxbuild =36
  integer,parameter::id_fullextrap=37,id_revise=38, id_n_dtmax  =39
  integer,parameter::id_n_ana    =40
  integer,parameter::id_lastint1=41  !1 after last id for integers

  integer, parameter::id_t       = 1,id_dtmax   = 2,id_deltnew  = 3
  integer, parameter::id_gamma   = 4,id_ekin    = 5,id_egrav    = 6
  integer, parameter::id_etherm  = 7,id_ascale  = 8,id_cosbox   = 9
  integer, parameter::id_vlength =10,id_alfstar =11,id_betastar =12
  integer, parameter::id_tol     =13,id_cfl     =14,id_dhmax    =15
  integer, parameter::id_treeacc =16,id_sepmax  =17,id_hmin     =18
  integer, parameter::id_eps     =19,id_dtinit  =20,id_tstop    =21
  integer, parameter::id_dt_out  =22,id_uconst  =23,id_xmas     =24
  integer, parameter::id_ylen    =25,id_umin    =26,id_rcore    =27
  integer, parameter::id_h_0     =28,id_omegam  =29,id_gmasslim =30
  integer, parameter::id_hmax    =31,id_uexpon  =32,id_ftol     =33
  integer, parameter::id_clhfrac =34,id_xmin    =35,id_xmax     =36
  integer, parameter::id_ymin    =37,id_ymax    =38,id_zmin     =39
  integer, parameter::id_zmax    =40
  integer, parameter::id_eosK    =41,id_rhozero =42,id_ftolpm   =43
  integer, parameter::id_tolpm   =44,id_tnextout=45,id_alfmax   =46
  integer, parameter::id_vtol    =47,id_vtolpm  =48,id_dtmin    =49
  integer, parameter::id_tlastout=50,id_dt_ana  =51,id_bsepmax  =52
  integer, parameter::id_taucool =53
  integer, parameter::id_lastreal1=54 !1 after last id for real numbers

end module vineread

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data,  only:npartoftype,dat,time,gamma,maxcol,maxpart,maxstep
  use params,         only:doub_prec
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc
  use labels,         only:ivx, iBfirst,ih,ipmass
  use mem_allocation, only:alloc
  use system_utils,   only:lenvironment
  use vineread,       only:id_gamma,id_iheadlen,id_ndim,id_npart,id_npart_sph,&
                           id_npoim,id_t
  implicit none
  integer, intent(in)          :: indexstart,ipos
  integer, intent(out)         :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: iheadlength
  integer :: i,j,ierr,nparti,ntoti,i1,icol
  integer :: npart_max,nstep_max,ncolstep,nptmass
  logical :: iexist,mhdread,useipindx
  character(len=len(rootname)+10)    :: dumpfile
  integer, parameter                 :: maxheadlength = 1000
  integer, dimension(maxheadlength)  :: iheader
  integer, dimension(:), allocatable :: ipindx,itstepbin

  !--we are assuming dump is double precision
  real(doub_prec), dimension(maxheadlength)    :: dheader
  real(doub_prec), dimension(:,:), allocatable :: dattemp, dattempvec
  real :: dum

  real :: hfactor

  nstepsread = 0
  npart_max = maxpart
  !--this is the default header length
  iheadlength = maxheadlength

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
  !--fix number of spatial dimensions
  !
  ndim = 3
  ndimV = 3
  mhdread = .false.
  if (lenvironment('VSPLASH_MHD') .or. lenvironment('VINE_MHD')) then
     mhdread = .true.
  endif
  if (lenvironment('VSPLASH_HFAC') .or. lenvironment('VINE_HFAC')) then
     hfactor = 2.8
  else
     hfactor = 1.0
  endif
  nstep_max = max(indexstart,1)

  j = indexstart
  nstepsread = 0
  nparti = 0
  ncolstep = 0

  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  if (mhdread) then
     print "(a)",' reading VINE MHD format'
  else
     print "(a)",' reading default VINE format (set VINE_MHD=yes for MHD)'
  endif
  !
  !--open the (unformatted) binary file and read the number of particles
  !
  open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
  else
     !
     !--read timestep header (integers only)
     !
     read(15,iostat=ierr) (iheader(i),i=1,iheadlength)
     !
     !--get number of particles from header and allocate memory
     !
     iheadlength = iheader(id_iheadlen)
     if (iheadlength.gt.maxheadlength) print "(a)",' ERROR: header length too big!'
     ntoti   = iheader(id_npart    )
     nparti  = iheader(id_npart_sph)
     nptmass = iheader(id_npoim    )
     ndim    = iheader(id_ndim     )
     if (ntoti.lt.nparti) then
        print*,' *** WARNING: ntotal < npart_sph in header, setting n_total=n_sph'
        ntoti = nparti
     endif
     if (nptmass.lt.0) then
        print*,' *** WARNING: error in nptmass read from header, nptmass = ',nptmass,' setting to 0'
        nptmass = 0
     endif
     if (nparti.le.0) then
        print*,' *** WARNING: error in npart read from header, npart = ',nparti
        ierr = 2
     endif
     if (ndim.le.0 .or. ndim.gt.3) then
        print*,' *** WARNING: error in ndim read from header, ndim = ',ndim
        ierr = 1
     endif

     ndimV   = ndim
     if (ndim.ne.3) print "(a,i1)",' number of dimensions = ',ndim
     if (mhdread) then
        ncolstep = 2*ndim + 6 + ndim
     else
        ncolstep = 2*ndim + 6
     endif
     ncolumns = ncolstep
     if ((.not.allocated(dat) .or. ntoti+nptmass.gt.npart_max) .and. ierr.eq.0) then
        if (.not.allocated(dat)) then
           npart_max = ntoti + nptmass
        else
           npart_max = max(npart_max,INT(1.1*(ntoti+nptmass)))
        endif
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
     !
     !--rewind file
     !
     rewind(15)
  endif
  if (ierr /= 0) then
     print "(/,a)", '  *** ERROR READING TIMESTEP HEADER: wrong endian? ***'
     print "(/,a)", '   (see splash userguide for compiler-dependent'
     print "(a)", '    ways to change endianness on the command line)'
     print "(/,a)", '   (set environment variable VINE_MHD to yes or TRUE '
     print "(a,/)", '    if you are trying to read MHD format)'
  else

       npart_max = max(npart_max,ntoti)
!
!--allocate/reallocate memory if j > maxstep
!
       if (j.gt.maxstep) then
          call alloc(maxpart,j+2*nstepsread,maxcol)
       endif
!
!--allocate a temporary array for double precision variables
!
       if (allocated(dattemp)) deallocate(dattemp)
       allocate(dattemp(npart_max,ncolstep),stat=ierr)
       dattemp = 0.
       if (ierr /= 0) print*,'not enough memory in read_data (dattemp)'
!
!--allocate a temporary array for vectors
!
       if (allocated(dattempvec)) deallocate(dattempvec)
       if (mhdread) then
          allocate(dattempvec(3*ndim+1,npart_max),stat=ierr)
       else
          allocate(dattempvec(2*ndim+1,npart_max),stat=ierr)
       endif
       dattempvec = 0.
       if (ierr /= 0) print*,'not enough memory in read_data (dattempvec)'
!
!--allocate a temporary array for particle index
!
       if (allocated(ipindx)) deallocate(ipindx)
       allocate(ipindx(npart_max),stat=ierr)
       !ipindx = 0
       if (ierr /= 0) print*,'not enough memory in read_data (ipindx)'
!
!--allocate a temporary array for itstepbin (MHD or point masses only)
!
       if (mhdread .or. nptmass.gt.0) then
          if (allocated(itstepbin)) deallocate(itstepbin)
          allocate(itstepbin(npart_max),stat=ierr)
          !itstepbin = 0
          if (ierr /= 0) print*,'not enough memory in read_data (itstepbin)'
       endif
!
!--now read the timestep data in the dumpfile
!
       write(*,"(a,i5,a)",advance="no") '| step ',j,': '

       ivx = ndim + 2 ! location of vx in 'columns'
!      starting point for non position and velocity columns
       icol = ndim + 1 + ndimV + 1

       if (mhdread) then
          if (nptmass.gt.0) then
             print "(a)",' WARNING: MHD format but point masses are present'
             print "(a)",'          and reading of point masses is not implemented'
             print "(a)",' *** Please email a copy of io.F so I can fix this ***  '
          endif
          iBfirst = icol+5
          read(15,iostat=ierr) &
               (iheader(i),i=1,iheadlength), &
               (dheader(i),i=1,iheadlength), &
               (dattempvec(1:ndim+1,i),i=1,ntoti), &
               (dattempvec(ivx:ivx+ndimV-1,i),i=1,ntoti), &
               (dattemp(i,icol), i=1,ntoti), &
               (dattemp(i,icol+1), i=1,nparti), &
               (dattemp(i,icol+2), i=1,nparti), &
               (dattemp(i,icol+3), i=1,nparti), &
               (dattemp(i,icol+4), i=1,ntoti), &
               (ipindx(i), i=1,ntoti), &
               (itstepbin(i),i=1,ntoti), &
               (dattempvec(ivx+ndimV:ivx+2*ndimV-1,i),i=1,nparti)
       else
          if (nptmass.gt.0) then
             !
             !--read point mass information at the end of the dump file
             !
             read(15,iostat=ierr) &
                  (iheader(i),i=1,iheadlength), &
                  (dheader(i),i=1,iheadlength), &
                  (dattempvec(1:ndim+1,i),i=1,ntoti), &
                  (dattempvec(ivx:ivx+ndimV-1,i),i=1,ntoti), &
                  (dattemp(i,icol), i=1,ntoti), &
                  (dattemp(i,icol+1), i=1,nparti), &
                  (dattemp(i,icol+2), i=1,nparti), &
                  (dattemp(i,icol+3), i=1,nparti), &
                  (dattemp(i,icol+4), i=1,ntoti), &
                  (ipindx(i), i=1,ntoti), &
                  (dum, i=1,nparti), &
                  (itstepbin(i), i=1,ntoti), &
                  (dattempvec(1:ndim+1,i),i=ntoti+1,ntoti+nptmass), &
                  (dattempvec(ivx:ivx+ndimV-1,i),i=ntoti+1,ntoti+nptmass), &
                  ((dum, i1=1,3),i=1,nptmass), &
                  (dattemp(i,icol), i=ntoti+1,ntoti+nptmass)
          else
             !
             !--no point masses, so shorter read
             !
             read(15,iostat=ierr) &
                  (iheader(i),i=1,iheadlength), &
                  (dheader(i),i=1,iheadlength), &
                  (dattempvec(1:ndim+1,i),i=1,ntoti), &
                  (dattempvec(ivx:ivx+ndimV-1,i),i=1,ntoti), &
                  (dattemp(i,icol), i=1,ntoti), &
                  (dattemp(i,icol+1), i=1,nparti), &
                  (dattemp(i,icol+2), i=1,nparti), &
                  (dattemp(i,icol+3), i=1,nparti), &
                  (dattemp(i,icol+4), i=1,ntoti), &
                  (ipindx(i), i=1,ntoti)
          endif
       endif

       if (ierr < 0) then
          print "(a)",'*** END OF FILE IN READ DATA ***'
       elseif (ierr /= 0) then
          if (mhdread) then
             print "(a)",'*** ERROR READING DATA: MAYBE NOT AN MHD FILE?? ***'
          else
             print "(a)",'*** ERROR READING DATA ***'
             print "(/,a)", '   (set environment variable VINE_MHD to yes or TRUE '
             print "(a,/)", '    if you are trying to read MHD format)'
          endif
       endif
       nstepsread = nstepsread + 1
!
!--spit out time
!
       time(j)  = real(dheader(id_t    ))
       gamma(j) = real(dheader(id_gamma))
       print "(a,es10.3,3(a,i8))",'t = ',time(j),' n(SPH) = ',ntoti,' n(Nbody) = ',ntoti-nparti,' n(star) = ',nptmass
!
!--check sanity of ipindx array: do not sort particles if values not sensible
!
       useipindx = .true.
       if (any(ipindx(1:ntoti).le.0 .or. ipindx(1:ntoti).gt.ntoti)) then
          print*,'WARNING: ipindx array has values < 0 or > ntot: particles not sorted'
          useipindx = .false.
       endif
!
!--convert posm and velocity vectors to columns and double to single precision
!
       do i=1,2*ndim+1
          if (useipindx) then
             dat(ipindx(1:ntoti),i,j) = real(dattempvec(i,1:ntoti))
          else
             dat(1:ntoti,i,j) = real(dattempvec(i,1:ntoti))
          endif
       enddo
       if (nptmass.gt.0) then
          do i=1,2*ndim+1
             dat(ntoti+1:ntoti+nptmass,i,j) = real(dattempvec(i,ntoti+1:ntoti+nptmass))
          enddo
       endif
!
!--convert B vectors to columns and double to single precision
!
       if (mhdread) then
          i1 = iBfirst - 1
          do i=ivx+ndimV,ivx+2*ndimV-1
             i1 = i1 + 1
             if (useipindx) then
                dat(ipindx(1:nparti),i1,j) = real(dattempvec(i,1:nparti))
             else
                dat(1:nparti,i1,j) = real(dattempvec(i,1:nparti))
             endif
          enddo
       endif
!
!--now convert scalars
!
       if (useipindx) then
          dat(ipindx(1:ntoti),icol:ncolstep,j) = real(dattemp(1:ntoti,icol:ncolstep))
       else
          dat(1:ntoti,icol:ncolstep,j) = real(dattemp(1:ntoti,icol:ncolstep))
       endif
       if (nptmass.gt.0) then
          dat(ntoti+1:ntoti+nptmass,icol,j) = real(dattemp(ntoti+1:ntoti+nptmass,icol))
       endif

       call set_labels
       if (ih.gt.0 .and. hfactor.gt.1.0) then
          dat(1:ntoti+nptmass,ih,j) = hfactor*dat(1:ntoti+nptmass,ih,j)
       endif

       if (nptmass.lt.10) then
          do i=1,nptmass
             if (ndim.eq.2) then
                print "('| point mass ',i1,': pos = (',es10.2,',',es10.2,'), mass = ',es10.2)", &
                   i,dat(ntoti+i,1:ndim,j),dat(ntoti+i,ipmass,j)
             else
                print "('| point mass ',i1,': pos = (',2(es10.2,','),es10.2,'), mass = ',es10.2)", &
                   i,dat(ntoti+i,1:ndim,j),dat(ntoti+i,ipmass,j)
             endif
          enddo
       endif

!
!--set particle numbers
!
       npartoftype(1,j) = nparti
       npartoftype(2,j) = ntoti - nparti
       npartoftype(3,j) = nptmass
!
!--clean up
!
       if (allocated(dattemp)) deallocate(dattemp)
       if (allocated(dattempvec)) deallocate(dattempvec)
       if (allocated(ipindx)) deallocate(ipindx)
       if (allocated(itstepbin)) deallocate(itstepbin)

  endif

 !
 !--reached end of file
 !
 close(15)
 if (nstepsread .gt. 0) then
    print*,'>> end of dump file: ntotal = ',sum(npartoftype(:,j))
 endif


return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,ih,ipmass,ivx,iutherm,irho,ix,iBfirst, &
                   labelvec,iamvec,labeltype
  use params
  use settings_data, only:ndim,ndimV,UseTypeInRenderings,ntypes
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
  ipmass = ndim+1   !  particle mass
  ivx = ndim+2

  ih = ndim + 1 + ndimV + 1       !  smoothing length
  iutherm = ih+1  !  thermal energy
  irho = iutherm+1     ! location of rho in data array

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = 'density'
  label(iutherm) = 'u'
  label(ih) = 'h'
  label(ipmass) = 'particle mass'
  label(irho+1) = 'alpha'
  label(irho+2) = 'potential energy'
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 3
  labeltype(1) = 'gas'
  labeltype(2) = 'Nbody'
  labeltype(3) = 'point mass'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.
  UseTypeInRenderings(3) = .false.

!-----------------------------------------------------------

  return
end subroutine set_labels
