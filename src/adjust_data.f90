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
!  Copyright (C) 2005-2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module adjustdata
 implicit none

contains

!----------------------------------------------------
!
!  specify data dependencies for adjust data
!  procedures to ensure required quantities
!  are read from the data file
!
!----------------------------------------------------
subroutine get_adjust_data_dependencies(required)
 use labels,        only:idustfrac,irho,ivx,ideltav
 use settings_data, only:ndimV,UseFakeDustParticles
 logical, intent(inout) :: required(0:)
 integer :: i

 ! if making fake dust particles, read epsilon and deltav
 if (idustfrac > 0 .and. irho > 0 .and. UseFakeDustParticles) then
    if (required(irho)) required(idustfrac) = .true.
    if (ideltav > 0) then
       do i=ivx,ivx+ndimV-1
          if (required(i)) required(ideltav+i-ivx) = .true.
       enddo
       if (any(required(ivx:ivx+ndimV-1))) required(idustfrac) = .true.
    endif
 endif

end subroutine get_adjust_data_dependencies
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
 use system_utils,    only:renvironment,envlist,ienvironment,lenvironment,ienvlist,&
                           ienvstring,get_environment_or_flag,get_command_flag
 use labels,          only:ih,ix,ivx,get_sink_type,ipmass,idustfrac,irho,labeltype,label
 use settings_data,   only:ncolumns,ndimV,ndim,ntypes,iverbose,UseFakeDustParticles,&
                           UseFastRender,icoords,track_string
 use particle_data,   only:dat,npartoftype,iamtype
 use filenames,       only:ifileopen,nstepsinfile
 use geometry,        only:labelcoord
 use part_utils,      only:locate_first_two_of_type,locate_nth_particle_of_type,&
                           locate_particle_from_string,get_binary,got_particles_of_type
 real :: hmin,dphi,domega,period
 real, dimension(3) :: x0,v0
 integer :: i,j,ierr,isink,isinkpos,itype,nlist,nerr
 integer :: ntot,isink1,isink2,isinklist(2),iorigin
 logical :: centreonsink,dontCentreVelocity,got_sinks,no_dust_particles
 character(len=20), dimension(3) :: list
 character(len=30) :: string

 !
 !--environment variable setting to enforce a minimum h
 !
 if (ih > 0 .and. ih <= ncolumns) then
    hmin = max(renvironment('SPLASH_HMIN_CODEUNITS',errval=-1.),&  ! for backwards compatibility
               renvironment('SPLASH_BEAM',errval=-1.))
    if (hmin > 0.) then
       if (.not.allocated(dat)) then
          print*,' INTERNAL ERROR: dat not allocated in adjust_data_codeunits'
          return
       endif
       print "(/,a,es10.3)",' >> SETTING MINIMUM H TO ',hmin
       where (dat(:,ih,:) < hmin .and. dat(:,ih,:) > 0.)
          dat(:,ih,:) = hmin
       end where
       print "(a)",' >> Switching accelerated rendering ON'
       UseFastRender = .true.
    endif
 endif

 if (ndim > 0) then
    !
    !--environment variable to corotate with first two sink particles
    !  can be SPLASH_COROTATE=true (just picks first two sinks)
    !      or SPLASH_COROTATE=1,3
    !
    isinklist = ienvlist('SPLASH_COROTATE',2)
    got_sinks = all(isinklist > 0)
    if (lenvironment('SPLASH_COROTATE') .or. got_sinks) then
       itype = get_sink_type(ntypes)
       if (itype > 0) then
          if (all(npartoftype(itype,:) < 2)) then
             print "(a)",' ERROR: --corotate set but less than 2 sink particles'
          else
             if (iverbose >= 1) print*
             if (got_sinks) then
                print "(a,i3,a,i3,a)",' :: COROTATING FRAME WITH SINKS ',isinklist(1),&
                                      ', ',isinklist(2),' from --corotate flag'
             else
                print "(a)",' :: COROTATING FRAME WITH FIRST 2 SINKS from --corotate flag'
             endif
             do j=1,nstepsinfile(ifileopen)
                if (got_sinks) then
                   call locate_nth_particle_of_type(isinklist(1),isink1,itype,iamtype(:,j),npartoftype(:,j),ntot)
                   call locate_nth_particle_of_type(isinklist(2),isink2,itype,iamtype(:,j),npartoftype(:,j),ntot)
                else
                   !  find first two sink particles in the data
                   call locate_first_two_of_type(isink1,isink2,itype,iamtype(:,j),npartoftype(:,j),ntot)
                endif
                !  get properties of the binary
                call get_binary(isink1,isink2,dat(:,:,j),x0,v0,dphi,domega,ndim,ndimV,ncolumns,ix,ivx,ipmass,iverbose,ierr)
                !  rotate all the particles into this frame
                if (ierr==0) call rotate_particles(dat(:,:,j),ntot,dphi,domega,x0(1:ndim),ndim,ndimV,v0)
             enddo
          endif
       else
          print "(a,/,a)",' ERROR: --corotate set but could not determine type ', &
                          '        corresponding to sink particles'
       endif
    endif

    !
    !--environment variable setting to centre plots on a selected sink particle
    !
    !--can specify either just "true" for sink #1, or specify a number for a particular sink
    centreonsink = lenvironment('SPLASH_CENTRE_ON_SINK') .or. lenvironment('SPLASH_CENTER_ON_SINK')
    isink        = max(ienvironment('SPLASH_CENTRE_ON_SINK'),ienvironment('SPLASH_CENTER_ON_SINK'))
    dontCentreVelocity = lenvironment('SPLASH_DONTCENTREVEL') .or. lenvironment('SPLASH_DONTCENTERVEL')
    if (isink > 0 .or. centreonsink .and. all(ix(1:ndim) > 0)) then
       if (isink==0) isink = 1
       itype = get_sink_type(ntypes)
       if (itype > 0) then
          if (all(npartoftype(itype,:) < isink)) then
             print "(a,i10,a)",' ERROR: --sink = ',isink,' but not enough sink particles'
          else
             if (iverbose >= 1) print*
             if (isink < 10) then
                print "(a,i1,a)",' :: CENTREING ON SINK ',isink,' from --sink flag'
             else
                print "(a,i3,a)",' :: CENTREING ON SINK ',isink,' from --sink flag'
             endif
             if (dontCentreVelocity) print "(a)",' :: NOT CENTREING VELOCITY'
             do j=1,nstepsinfile(ifileopen)
                call locate_nth_particle_of_type(isink,isinkpos,itype,iamtype(:,j),npartoftype(:,j),ntot)
                call centre_on_particle(isinkpos,dat(:,:,j),ntot,ndim,ndimV,ncolumns,dontCentreVelocity,iverbose,label='sink')
             enddo
          endif
       else
          print "(a,/,a)",' ERROR: --sink set but could not determine type ', &
                          '        corresponding to sink particles'
       endif
    endif
    !
    !--center on a particular particle
    !
    call get_environment_or_flag('SPLASH_ORIGIN',string)
    iorigin = ienvstring(string) ! first try to read the string as an integer

    if (iorigin > 0 .or. (len_trim(string) > 0) .and. all(ix(1:ndim) > 0) .and. ndim > 0) then
       if (track_string(1:1) /= '0' .and. len_trim(track_string) > 0) then
          print "(/,a,/)",' ERROR: cannot use --track and --origin at the same time, disabling particle tracking'
          track_string = '0'
       endif
       if (centreonsink)  then
          print "(/,a,/)",' ERROR: cannot use --sink and --origin at the same time'
       else
          do j=1,nstepsinfile(ifileopen)
             !
             !--handle strings like --origin=maxdens to locate the particle for the origin
             !
             ntot = sum(npartoftype(:,j))
             if (len_trim(string) > 0 .and. iorigin==0) then
                iorigin = locate_particle_from_string(string,ntot,dat(:,:,j),irho)
                if (iorigin <= 0) exit ! quit loop over steps
             endif
             if (j==1 .or. ienvstring(string) == 0) then
                print "(a,i0,a)",' :: CENTREING ON PARTICLE ',iorigin,' from --origin='//trim(string)//' flag'
                if (dontCentreVelocity) print "(a)",' :: NOT CENTREING VELOCITY'
             endif
             !
             !--now centre on the chosen particle
             !
             call centre_on_particle(iorigin,dat(:,:,j),ntot,ndim,ndimV,ncolumns,dontCentreVelocity,iverbose,label='')
          enddo
       endif
    endif

    !
    !--environment variable setting to subtract a mean velocity
    !
    if (ivx > 0 .and. ivx+ndimV-1 <= ncolumns) then
       call envlist('SPLASH_VZERO',nlist,list)
       nerr = 0
       if (nlist > 0 .and. nlist < ndimV) then
          print "(/,2(a,i1))",' >> ERROR in --vzero setting: number of components = ',nlist,', needs to be ',ndimV
          nerr = 1
       elseif (nlist > 0) then
          if (nlist > ndimV) print "(a,i1,a,i1)",' >> WARNING! --vzero setting has ',nlist, &
                                                  ' components: using only first ',ndimV
          nerr = 0
          do i=1,ndimV
             read(list(i),*,iostat=ierr) v0(i)
             if (ierr /= 0) then
                print "(a)",' >> ERROR reading v'//trim(labelcoord(i,icoords))//&
                            ' component from --vzero setting'
                nerr = ierr
             endif
          enddo
          if (nerr==0) then
             print "(a)",' >> SUBTRACTING MEAN VELOCITY (from --vzero setting):'
             if (.not.allocated(dat) .or. size(dat(1,:,1)) < ivx+ndimV-1) then
                print*,' INTERNAL ERROR: dat not allocated in adjust_data_codeunits'
                return
             endif
             do i=1,ndimV
                print "(4x,a,es10.3)",trim(label(ivx+i-1))//' = '//trim(label(ivx+i-1))//' - ',v0(i)
                dat(:,ivx+i-1,:) = dat(:,ivx+i-1,:) - v0(i)
             enddo
          endif
       endif
       if (nerr /= 0) then
          print "(4x,a)",'WARNING: --vzero setting ignored'
       endif
    endif
 endif
 !
 !--fake a set of dust particles from the one-fluid dust formulation
 !
 no_dust_particles = .not.got_particles_of_type('dust',labeltype,npartoftype)
 if (idustfrac > 0 .and. irho > 0 .and. no_dust_particles) then
    if (UseFakeDustParticles) then
       call fake_twofluids(1,nstepsinfile(ifileopen),ndim,ndimV,dat,npartoftype,iamtype)
    elseif (iverbose >= 1) then
       print "(a)",' One fluid dust: set option in d) menu to make fake dust particles'
    endif
 endif

 !
 !--phase fold column 1 with a given period
 !
 period = renvironment('SPLASH_PHASEFOLD')
 if (period > 0.) then
    print "(a,es10.3)",' PHASE FOLDING column 1 with period of ',period
    dat(:,1,1) = dat(:,1,1) - period*int(dat(:,1,1)/period)
 endif


end subroutine adjust_data_codeunits

!-----------------------------------------------------------------
! routine to rotate particles with a given cylindrical angle dphi
!-----------------------------------------------------------------
pure subroutine rotate_particles(dat,np,dphi,domega,x0,ndim,ndimV,v0)
 use labels, only:ix,ivx
 integer, intent(in) :: np,ndim,ndimV
 real,    intent(in) :: dphi,domega
 real,    dimension(:,:), intent(inout) :: dat
 real, dimension(ndim), intent(in) :: x0
 real, dimension(ndimV), intent(in) :: v0
 real, dimension(ndim)  :: xi
 real, dimension(ndimV) :: vi
 real :: r,phi,xnew,ynew,cosp,sinp,vr,vphi
 integer :: i

 !--rotate positions
 do i=1,np
    xi  = dat(i,ix(1:ndim)) - x0(1:ndim)
    r   = sqrt(xi(1)**2 + xi(2)**2)
    phi = atan2(xi(2),xi(1))
    phi = phi + dphi
    cosp = cos(phi)
    sinp = sin(phi)
    xnew = r*cosp
    ynew = r*sinp
    dat(i,ix(1)) = xnew
    dat(i,ix(2)) = ynew
    !--rotate velocities, if present
    if (ivx > 0) then
       vi = dat(i,ivx:ivx+ndimV-1) - v0
       vr = vi(1)*xi(1)/r + vi(2)*xi(2)/r
       vphi = (vi(1)*(-xi(2)/r) + vi(2)*xi(1)/r) !- r*domega
       dat(i,ivx)   = vr*cosp - vphi*sinp
       dat(i,ivx+1) = vr*sinp + vphi*cosp
    endif
 enddo

end subroutine rotate_particles

!------------------------------------------------------
! routine to centre particle positions and velocities
! on a particle location
!------------------------------------------------------
subroutine centre_on_particle(ipart,dat,np,ndim,ndimV,ncol,dontCentreVelocity,iverbose,label)
 use labels, only:ix,ivx
 integer, intent(in) :: ipart,np,ndim,ndimV,ncol,iverbose
 logical, intent(in) :: dontCentreVelocity
 real, dimension(:,:), intent(inout) :: dat
 real, dimension(ndim)  :: x0
 real, dimension(ndimV) :: v0
 character(len=*), intent(in) :: label

 !--sanity check inputs
 if (ipart <= 0) then
    print "(a)",' ERROR: could not locate '//trim(label)//' particle in dat array'
    return
 endif
 if (ipart > np .or.  ipart > size(dat(:,1))) then
    print "(/,a,/)",' ERROR: '//trim(label)//' particle number exceeds number of particles'
    return ! skip if iorigin > npart
 endif
 !--make positions relative to chosen particle
 x0(1:ndim) = dat(ipart,ix(1:ndim))
 if (iverbose >= 1) print "(a,3(1x,es10.3))",' :: '//trim(label)//' position =',x0(1:ndim)
 !--make velocities relative to origin particle
 if ((ivx > 0 .and. ivx+ndimV-1 <= ncol) .and. (.not. dontCentreVelocity)) then
    v0(1:ndimV) = dat(ipart,ivx:ivx+ndimV-1)
    if (iverbose >= 1) print "(a,3(1x,es10.3))",' :: '//trim(label)//' velocity =',v0(1:ndimV)
 else
    v0 = 0.
 endif
 call shift_particles(dat(:,:),np,ndim,ndimV,ncol,x0,v0)

end subroutine centre_on_particle

!------------------------------------------------------
! routine to shift particle positions and velocities
! to new location
!------------------------------------------------------
pure subroutine shift_particles(dat,np,ndim,ndimV,ncol,x0,v0)
 integer, intent(in) :: np,ndim,ndimV,ncol
 real, dimension(:,:), intent(inout) :: dat
 real, dimension(ndim),  intent(in) :: x0
 real, dimension(ndimV), intent(in) :: v0

 call shift_positions(dat,np,ndim,x0)
 call shift_velocities(dat,np,ndimV,ncol,v0)

end subroutine shift_particles

!------------------------------------------------------
! routine to shift particle positions to new location
!------------------------------------------------------
pure subroutine shift_positions(dat,np,ndim,x0)
 use labels, only:ix
 integer, intent(in) :: np,ndim
 real,    dimension(:,:), intent(inout) :: dat
 real, dimension(ndim),  intent(in) :: x0
 integer :: icol

 !--shift positions
 do icol=1,ndim
    dat(1:np,ix(icol)) = dat(1:np,ix(icol)) - x0(icol)
 enddo

end subroutine shift_positions

!------------------------------------------------------
! routine to shift particle velocities by constant
!------------------------------------------------------
pure subroutine shift_velocities(dat,np,ndimV,ncol,v0)
 use labels, only:ivx
 integer, intent(in) :: np,ndimV,ncol
 real,    dimension(:,:), intent(inout) :: dat
 real, dimension(ndimV), intent(in) :: v0
 integer :: icol

 !--make velocities relative to sink particle
 if (ivx > 0 .and. ivx+ndimV-1 <= ncol) then
    do icol=1,ndimV
       dat(1:np,ivx+icol-1) = dat(1:np,ivx+icol-1) - v0(icol)
    enddo
 endif

end subroutine shift_velocities

!------------------------------------------------------
!
! routine to fake a set of dust particles from
! the one fluid dust method
!
!------------------------------------------------------
subroutine fake_twofluids(istart,iend,ndim,ndimV,dat,npartoftype,iamtype)
 use params,         only:int1
 use labels,         only:idustfrac,irho,ix,ih,ipmass,ivx,ideltav
 use mem_allocation, only:alloc
 use particle_data,  only:maxpart,maxstep,maxcol
 use settings_data,  only:iverbose,ndusttypes,idustfrac_plot,ideltav_plot
 integer,                       intent(in)    :: istart,iend,ndim,ndimV
 real,    dimension(:,:,:),     intent(inout), allocatable :: dat
 integer, dimension(:,:),       intent(inout), allocatable :: npartoftype
 integer(int1), dimension(:,:), intent(inout), allocatable :: iamtype
 integer :: ndust,jdust,ntoti,i,j
 integer :: idustfrac_temp,ideltav_temp
 real    :: rhodust,rhogas,rhotot,dustfraci,gasfraci,pmassgas,pmassdust,pmassj
 real, dimension(ndimV) :: veli,vgas,vdust,deltav
 logical :: use_vels

 if (idustfrac > 0 .and. irho > 0) then
    !
    !--determine which dust fraction is being used to create the fake dust particles
    !
    !if (iverbose >= 0) print*,' got dustfrac in column ',idustfrac
    if (ndusttypes>1) then
       if (idustfrac_plot == 0 ) then
          idustfrac_temp = idustfrac
          idustfrac_plot = idustfrac_plot
          ideltav_temp   = ideltav
          ideltav_plot   = ideltav
       else
          idustfrac_temp = idustfrac_plot
          ideltav_temp   = ideltav_plot
       endif
    else
       idustfrac_temp = idustfrac
       ideltav_temp   = ideltav
    endif
    !
    !--create explicit dust phase from gas and one-fluid dust properties
    !
    do i=istart,iend
       ntoti = sum(npartoftype(:,i))
       if (.not.allocated(dat) .or. (ntoti + npartoftype(1,i)) > maxpart) then
          call alloc(ntoti + npartoftype(1,i),maxstep,maxcol,mixedtypes=.true.)
       endif
       if (npartoftype(2,i) > 0) cycle
       ndust = 0
       !--zero the properties of newly created dust particles
       dat(ntoti+1:ntoti+npartoftype(1,i),:,i) = 0.
       if (idustfrac_temp > size(dat(1,:,1)) .or. idustfrac_temp <= 0) then
          print*,' ERROR: idustfrac out of range: cannot create fake dust particles'
          return
       endif
       use_vels = (ideltav_temp > 0 .and. ivx > 0 .and. ndimV > 0)
       do j=1,ntoti
          if (iamtype(j,i)==1) then
             ndust = ndust + 1 ! one dust particle for every gas particle
             rhotot  = dat(j,irho,i)
             dustfraci = dat(j,idustfrac_temp,i)
             gasfraci = 1. - dustfraci
             rhogas  = rhotot*gasfraci
             rhodust = rhotot*dustfraci
             !--replace global properties with gas-only stuff
             dat(j,irho,i) = rhogas
             dat(j,idustfrac,i) = 0.    ! dust fraction = 0 on gas particles
             !--copy x, smoothing length onto dust particle
             jdust = ntoti + ndust

             !--fill in dust properties
             if (ndim > 0) dat(jdust,ix(1:ndim),i) = dat(j,ix(1:ndim),i)
             if (ih > 0)   dat(jdust,ih,i)         = dat(j,ih,i)
             if (irho > 0) dat(jdust,irho,i)       = rhodust
             dat(jdust,idustfrac,i) = 1. ! dust fraction = 1 on dust particles
             iamtype(ntoti + ndust,i) = 2

             !--particle masses
             if (ipmass > 0) then
                pmassj    = dat(j,ipmass,i)
                pmassgas  = pmassj*gasfraci
                pmassdust = pmassj*dustfraci
                dat(j,ipmass,i)     = pmassgas
                dat(jdust,ipmass,i) = pmassdust
             endif

             !--velocities
             if (use_vels) then
                veli(:)      = dat(j,ivx:ivx+ndimV-1,i)
                deltav(:)    = dat(j,ideltav_temp:ideltav_temp+ndimV-1,i)
                vgas(:)      = veli(:) - (1. - gasfraci)*deltav(:)
                vdust(:)     = veli(:) + gasfraci*deltav(:)
                dat(j,ivx:ivx+ndimV-1,i)     = vgas(:)
                dat(jdust,ivx:ivx+ndimV-1,i) = vdust(:)
                if (ideltav /= 0) then
                    !--set deltav to zero on gas and dust particles
                    dat(j,ideltav_temp:ideltav_temp+ndimV-1,i) = deltav(:)
                    dat(jdust,ideltav:ideltav+ndimV-1,i) = deltav(:)
                 endif
                !if (abs(deltav(1)) > 0.3) print*,' particle ',j,'->',jdust,' vg = ',vgas(:),' vd = ',vdust(:),' v = ',veli,deltav
             endif
          endif
       enddo
       if (iverbose >= 1) then
          print "(a,i8,a)",' Creating ',ndust,' fake dust particles'
          if (.not.use_vels .and. ivx > 0) print "(a)",' WARNING: deltav not found in one fluid dust data: cannot get vels'
       endif
       npartoftype(2,i) = npartoftype(2,i) + ndust
    enddo
 else
    print "(a)",' ERROR: could not locate dust-to-gas ratio and/or density'
 endif

end subroutine fake_twofluids

end module adjustdata
