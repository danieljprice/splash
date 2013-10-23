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
module adjustdata
 implicit none

contains
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
 use system_utils,  only:renvironment,envlist,ienvironment,lenvironment
 use labels,        only:ih,ivx,label,ix,get_sink_type,ipmass
 use settings_data, only:ncolumns,ndimV,icoords,ndim,debugmode,ntypes,iverbose
 use particle_data, only:dat,npartoftype,iamtype
 use geometry,      only:labelcoord
 use filenames,     only:ifileopen,nstepsinfile
 use part_utils,    only:locate_first_two_of_type,locate_nth_particle_of_type
 implicit none
 real :: hmin,m1,m2,dmtot,dphi
 real, dimension(3) :: vsink,xyzsink,xyzsink1,xyzsink2,x0,dx,v0,v1,v2
 character(len=20), dimension(3) :: list
 integer :: i,j,nlist,nerr,ierr,isink,isinkpos,itype
 integer :: isink1,isink2,ntot
 logical :: centreonsink

 !
 !--environment variable setting to enforce a minimum h
 !
 if (ih.gt.0 .and. ih.le.ncolumns) then
    hmin = renvironment('SPLASH_HMIN_CODEUNITS',errval=-1.)
    if (hmin.gt.0.) then
       if (.not.allocated(dat)) then
          print*,' INTERNAL ERROR: dat not allocated in adjust_data_codeunits'
          return
       endif
       print "(/,a,es10.3)",' >> SETTING MINIMUM H TO ',hmin
       where (dat(:,ih,:) < hmin .and. dat(:,ih,:).gt.0.)
          dat(:,ih,:) = hmin
       end where
    endif
 endif

 !
 !--environment variable setting to subtract a mean velocity
 !
 if (ivx.gt.0 .and. ivx+ndimV-1.le.ncolumns) then
    call envlist('SPLASH_VZERO_CODEUNITS',nlist,list)
    nerr = 0
    if (nlist.gt.0 .and. nlist.lt.ndimV) then
       print "(/,2(a,i1))",' >> ERROR in SPLASH_VZERO_CODEUNITS setting: number of components = ',nlist,', needs to be ',ndimV
       nerr = 1
    elseif (nlist.gt.0) then
       if (nlist.gt.ndimV) print "(a,i1,a,i1)",' >> WARNING! SPLASH_VZERO_CODEUNITS setting has ',nlist, &
                                               ' components: using only first ',ndimV
       nerr = 0
       do i=1,ndimV
          read(list(i),*,iostat=ierr) v0(i)
          if (ierr.ne.0) then
             print "(a)",' >> ERROR reading v'//trim(labelcoord(i,icoords))//&
                         ' component from SPLASH_VZERO_CODEUNITS setting'
             nerr = ierr
          endif
       enddo
       if (nerr.eq.0) then
          print "(a)",' >> SUBTRACTING MEAN VELOCITY (from SPLASH_VZERO_CODEUNITS setting):'
          if (.not.allocated(dat) .or. size(dat(1,:,1)).lt.ivx+ndimV-1) then
             print*,' INTERNAL ERROR: dat not allocated in adjust_data_codeunits'
             return
          endif
          do i=1,ndimV
             print "(4x,a,es10.3)",trim(label(ivx+i-1))//' = '//trim(label(ivx+i-1))//' - ',v0(i)
             dat(:,ivx+i-1,:) = dat(:,ivx+i-1,:) - v0(i)
          enddo
       endif
    endif
    if (nerr.ne.0) then
       print "(4x,a)",'SPLASH_VZERO_CODEUNITS setting not used'
    endif
 endif
 if (ndim.gt.0) then
    !
    !--environment variable to corotate with first two sink particles
    !
    if (lenvironment('SPLASH_COROTATE')) then
       itype = get_sink_type(ntypes)
       if (itype.gt.0) then
          if (all(npartoftype(itype,:).lt.2)) then
             print "(a)",' ERROR: SPLASH_COROTATE set but less than 2 sink particles'
          else
             if (iverbose.ge.1) print*
             print "(a,i3,a)",' :: COROTATING FRAME WITH FIRST 2 SINKS from SPLASH_COROTATE setting'
             do j=1,nstepsinfile(ifileopen)
                !--find first two sink particles in the data
                call locate_first_two_of_type(isink1,isink2,itype,iamtype(:,j),npartoftype(:,j),ntot)

                if (isink1 <= 0 .or. isink2 <= 0) then
                   if (debugmode) print*,' sink1 = ',isink1,' sink2 = ',isink2
                   print "(a)",' ERROR locating sink particles in the data'
                else
                   xyzsink1  = 0.
                   xyzsink2 = 0.
                   xyzsink1(1:ndim) = dat(isink1,ix(1:ndim),j)
                   xyzsink2(1:ndim) = dat(isink2,ix(1:ndim),j)
                   !--get centre of mass
                   if (ipmass > 0 .and. ipmass <= ncolumns) then
                      m1 = dat(isink1,ipmass,j)
                      m2 = dat(isink2,ipmass,j)
                   else
                      m1 = 1.
                      m2 = 1.
                   endif
                   dmtot = 1./(m1 + m2)
                   x0 = (m1*xyzsink1 + m2*xyzsink2)*dmtot
                   if (iverbose >= 1) then
                      print "(a,3(1x,es10.3),a,es10.3)",' :: sink 1 pos =',xyzsink1(1:ndim),' m = ',m1
                      print "(a,3(1x,es10.3),a,es10.3)",' :: sink 2 pos =',xyzsink2(1:ndim),' m = ',m2
                      print "(a,3(1x,es10.3))",' :: c. of mass =',x0(1:ndim)
                   endif
                   !--work out angle needed to rotate into corotating frame
                   dx   = xyzsink1 - x0
                   dphi = atan2(dx(2),dx(1))

                   !--get velocities
                   if (ivx > 0 .and. ivx + ndimV <= ncolumns) then
                      v1 = dat(isink1,ivx:ivx+ndimV-1,j)
                      v2 = dat(isink2,ivx:ivx+ndimV-1,j)
                      v0 = (m1*v1 + m2*v2)*dmtot
                      if (iverbose >= 1) print "(a,3(1x,es10.3))",' :: vel c of m =',v0(1:ndimV)
                   else
                      v0 = 0.
                   endif
                   call rotate_particles(dat(:,:,j),ntot,dphi,x0(1:ndim),ndim,ndimV,v0)
                endif
             enddo
          endif
       else
          print "(a,/,a)",' ERROR: SPLASH_COROTATE set but could not determine type ', &
                          '        corresponding to sink particles'
       endif
    endif

    !
    !--environment variable setting to centre plots on a selected sink particle
    !
    !--can specify either just "true" for sink #1, or specify a number for a particular sink
    centreonsink = lenvironment('SPLASH_CENTRE_ON_SINK') .or. lenvironment('SPLASH_CENTER_ON_SINK')
    isink        = max(ienvironment('SPLASH_CENTRE_ON_SINK'),ienvironment('SPLASH_CENTER_ON_SINK'))
    if (isink.gt.0 .or. centreonsink) then
       if (isink.eq.0) isink = 1
       itype = get_sink_type(ntypes)
       if (itype.gt.0) then
          if (all(npartoftype(itype,:).lt.isink)) then
             print "(a,i10,a)",' ERROR: SPLASH_CENTRE_ON_SINK = ',isink,' but not enough sink particles'
          else
             if (iverbose.ge.1) print*
             if (isink.lt.10) then
                print "(a,i1,a)",' :: CENTREING ON SINK ',isink,' from SPLASH_CENTRE_ON_SINK setting'          
             else
                print "(a,i3,a)",' :: CENTREING ON SINK ',isink,' from SPLASH_CENTRE_ON_SINK setting'
             endif
             do j=1,nstepsinfile(ifileopen)
                call locate_nth_particle_of_type(isink,isinkpos,itype,iamtype(:,j),npartoftype(:,j),ntot)
                if (isinkpos.eq.0) then
                   print "(a)",' ERROR: could not locate sink particle in dat array'
                else
                   if (debugmode) print*,' SINK POSITION = ',isinkpos,npartoftype(1:itype,j)
                   !--make positions relative to sink particle
                   xyzsink(1:ndim) = dat(isinkpos,ix(1:ndim),j)
                   if (iverbose.ge.1) print "(a,3(1x,es10.3))",' :: sink position =',xyzsink(1:ndim)
                   !--make velocities relative to sink particle
                   if (ivx.gt.0 .and. ivx+ndimV-1.le.ncolumns) then
                      vsink(1:ndimV) = dat(isinkpos,ivx:ivx+ndimV-1,j)
                      if (iverbose.ge.1) print "(a,3(1x,es10.3))",' :: sink velocity =',vsink(1:ndimV)
                   else
                      vsink = 0.
                   endif
                   call shift_particles(dat(:,:,j),ntot,ndim,ndimV,ncolumns,xyzsink,vsink)
                endif
             enddo
          endif
       else
          print "(a,/,a)",' ERROR: SPLASH_CENTRE_ON_SINK set but could not determine type ', &
                          '        corresponding to sink particles'
       endif
    endif
 endif

end subroutine adjust_data_codeunits

!-----------------------------------------------------------------
! routine to rotate particles with a given cylindrical angle dphi
!-----------------------------------------------------------------
pure subroutine rotate_particles(dat,np,dphi,x0,ndim,ndimV,v0)
 use labels, only:ix,ivx
 integer, intent(in) :: np,ndim,ndimV
 real,    intent(in) :: dphi
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
    phi = phi - dphi
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
       vphi = vi(1)*(-xi(2)/r) + vi(2)*xi(1)/r !- r*omega
       dat(i,ivx)   = vr*cosp - vphi*sinp
       dat(i,ivx+1) = vr*sinp + vphi*cosp
    endif
 enddo

end subroutine rotate_particles

!------------------------------------------------------
! routine to shift particle positions and velocities
! to new location
!------------------------------------------------------
pure subroutine shift_particles(dat,np,ndim,ndimV,ncol,x0,v0)
 integer, intent(in) :: np,ndim,ndimV,ncol
 real,    dimension(:,:), intent(inout) :: dat
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
! routine to shift particle positions to new location
!------------------------------------------------------
pure subroutine shift_velocities(dat,np,ndimV,ncol,v0)
 use labels, only:ix,ivx
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

end module adjustdata
