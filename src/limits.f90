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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------
!
! subroutines to do with setting of plot limits from data
! and using only a subset of the particles according to a
! range in parameters
!
!----------------------------------------------------------
module limits
 use params
 implicit none
 real, dimension(maxplot,2) :: lim,range,lim2
 private :: warn_minmax
 
 public

contains

!----------------------------------------------------------
! set plot limits for all columns
! NB: does not differentiate between particle types at the moment
!----------------------------------------------------------
subroutine set_limits(ifromstep,itostep,ifromcol,itocol)
  use labels,        only:label,ix
  use geometry,      only:coord_transform_limits
  use particle_data, only:npartoftype,dat,maxcol
  use settings_data, only:ndim,icoords,icoordsnew
  integer, intent(in) :: ifromstep,itostep,ifromcol,itocol
  integer :: i,j,k,ntoti,itocoli

  print 100,ifromstep,itostep,ifromcol,itocol
100 format(/' setting plot limits: steps ',i5,'->',i5,' cols ',i2,'->',i3)
  if (ifromcol.gt.maxcol .or. maxcol.eq.0) then
     print "(a)",' *** error: set_limits: column > array size ***'
     return
  endif
  if (ifromcol.gt.itocol) then
     print "(a)",' *** error in call to set_limits: begin column > end column'
     return
  endif
  itocoli = itocol
  if (itocol.gt.maxcol) then
     print "(a,i3,a)",' *** warning: set_limits: only setting limits up to column ',maxcol,' ***'
     itocoli = maxcol
  endif
  !!--find limits of particle properties
  lim(ifromcol:itocol,1) = huge(lim)
  lim(ifromcol:itocol,2) = -huge(lim)
  do i=ifromstep,itostep
     ntoti = sum(npartoftype(:,i))
     do j=ifromcol,itocoli
        do k=1,ntoti
           lim(j,1) = min(lim(j,1),dat(k,j,i))
           lim(j,2) = max(lim(j,2),dat(k,j,i))
        enddo
     enddo
  enddo
  !
  !--warn if limits are the same
  !
  do j=ifromcol,itocol
     call warn_minmax(label(j),lim(j,1),lim(j,2))
  enddo
  !print "(a/)",' plot limits set'

  lim2(ifromcol:itocol,:) = 0.

  !
  !--transform coord limits into new coordinate system if coord transform is applied
  !
  if (icoordsnew.ne.icoords .and. ndim.gt.0) then
     if (ifromcol.le.ix(ndim)) then ! separate if is to avoid referencing ix(0) if ndim=0
        call coord_transform_limits(lim(ix(1):ix(ndim),1),lim(ix(1):ix(ndim),2), &
                                    icoords,icoordsnew,ndim)
     endif
  endif

  return
end subroutine set_limits

!----------------------------------------------------------
! save plot limits for all columns to a file
!----------------------------------------------------------
subroutine write_limits(limitsfile)
  use settings_data, only:numplot,ndataplots
  character(len=*), intent(in) :: limitsfile
  integer :: i

  print*,'saving plot limits to file ',trim(limitsfile)

  open(unit=55,file=limitsfile,status='replace',form='formatted',ERR=998)
  do i=1,numplot
     if (rangeset(i) .and. i.lt.ndataplots .and. lim2set(i)) then
        write(55,"(6(1x,1pe14.6))",err=999) lim(i,1),lim(i,2),range(i,1),range(i,2),lim2(i,1),lim2(i,2)
     elseif (lim2set(i) .and. i.lt.ndataplots) then
        write(55,"(6(1x,1pe14.6))",err=999) lim(i,1),lim(i,2),0.,0.,lim2(i,1),lim2(i,2)
     elseif (rangeset(i) .and. i.lt.ndataplots) then
        write(55,"(4(1x,1pe14.6))",err=999) lim(i,1),lim(i,2),range(i,1),range(i,2)
     else
        write(55,"(2(1x,1pe14.6))",err=999) lim(i,1),lim(i,2)
     endif
  enddo
  close(unit=55)

  return

998 continue
  print*,'*** error opening limits file: limits not saved'
  return
999 continue
  print*,'*** error saving limits'
  close(unit=55)
  return

end subroutine write_limits

!----------------------------------------------------------
! read plot limits for all columns from a file
!----------------------------------------------------------
subroutine read_limits(limitsfile,ierr)
  use labels,        only:label
  use settings_data, only:numplot,ncolumns,ncalc
  use asciiutils,    only:ncolumnsline
  character(len=*), intent(in) :: limitsfile
  integer,         intent(out) :: ierr
  integer                      :: i,ncolsline
  character(len=120)           :: line
  logical :: iexist

  ierr = 0

  inquire(file=limitsfile,exist=iexist)
  if (.not.iexist) then
     print "(1x,a)",trim(limitsfile)//' not found'
     ierr = 1
     return
  endif

  open(unit=54,file=limitsfile,status='old',form='formatted',err=997)
  print "(a)",' read '//trim(limitsfile)
  do i=1,numplot
     read(54,"(a)",err=998,end=999) line
     ncolsline = ncolumnsline(line)
     if (ncolsline.lt.2) then
        goto 998
     elseif (ncolsline.ge.6) then
        read(line,*,err=998,end=999) lim(i,1),lim(i,2),range(i,1),range(i,2),lim2(i,1),lim2(i,2)
     elseif (ncolsline.ge.4) then
        read(line,*,err=998,end=999) lim(i,1),lim(i,2),range(i,1),range(i,2)
     else
        read(line,*,err=998,end=999) lim(i,1),lim(i,2)
     endif
     call assert_sensible_limits(lim(i,1),lim(i,2))
     !
     !--warn if limits are the same
     !
     call warn_minmax(label(i),lim(i,1),lim(i,2))
  enddo
  close(unit=54)
  return

997 continue
  print*,trim(limitsfile),' not found'
  ierr = 1
  return
998 continue
  call print_rangeinfo()
  call print_lim2info()
  print*,'*** error reading limits from file'
  ierr = 2
  close(unit=54)
  return
999 continue
  !--only give error if we really do not have enough columns
  !  (on first call nextra is not set)
  if (i.le.ncolumns+ncalc) then
     print "(a,i3)",' end of file in '//trim(limitsfile)//': limits read to column ',i
     ierr = -1
  endif

  !--print info about range restrictions read from file
  call print_rangeinfo()
  call print_lim2info()
  close(unit=54)
  return

end subroutine read_limits

!----------------------------------------------------------
! get a subset of the particles by enforcing range restrictions
!----------------------------------------------------------
subroutine get_particle_subset(icolours,datstep,ncolumns)
 use labels, only:label
 integer, intent(inout) :: icolours(:)
 real,    intent(in)    :: datstep(:,:)
 integer, intent(in)    :: ncolumns
 integer :: icol

 if (anyrangeset()) then
    !--reset colours of all particles (to not hidden) if using range restriction
    where (icolours(:).eq.-1000)
       icolours(:) = 0
    elsewhere
       icolours(:) = abs(icolours(:))
    endwhere

    do icol=1,ncolumns
       if (rangeset(icol)) then
          print "(a,1pe10.3,a,1pe10.3,a)",' | using only particles in range ', &
                range(icol,1),' < '//trim(label(icol))//' < ',range(icol,2),' |'
       !
       !--loop over the particles and colour those outside the range
       !  NB: background colour (0) is set to -1000
       !
          where (datstep(:,icol).lt.range(icol,1) .or. &
                 datstep(:,icol).gt.range(icol,2))
             where (icolours.eq.0)
                icolours = -1000
             elsewhere
                icolours = -abs(icolours)
             end where
          end where
       endif
    enddo
 endif

 return
end subroutine get_particle_subset

!----------------------------------------------------------
! reset all range restrictions to zero
!----------------------------------------------------------
subroutine reset_all_ranges()
 use particle_data, only:icolourme

 print "(a)",' removing all range restrictions '
 where (icolourme(:).eq.-1000)
    icolourme(:) = 0
 elsewhere
    icolourme(:) = abs(icolourme(:))
 endwhere
 range(:,:) = 0.

 return
end subroutine reset_all_ranges

!----------------------------------------------------------
! function which returns whether or not a range
! has been set for a given column
!----------------------------------------------------------
logical function rangeset(icol)
 integer, intent(in) :: icol

 rangeset = .false.
 if (abs(range(icol,2)-range(icol,1)).gt.tiny(range)) rangeset = .true.

 return
end function rangeset

!----------------------------------------------------------
! function which returns whether or not lim2
! has been set for a given column
!----------------------------------------------------------
logical function lim2set(icol)
 integer, intent(in) :: icol

 lim2set = .false.
 if (abs(lim2(icol,2)).gt.tiny(lim2) .or. abs(lim2(icol,1)).gt.tiny(lim2)) lim2set = .true.

 return
end function lim2set

!----------------------------------------------------------
! reset all range restrictions to zero
!----------------------------------------------------------
subroutine reset_lim2(icol)
 integer, intent(in) :: icol

 print "(a)",' contour limits same as render limits'
 if (icol.gt.0 .and. icol.le.maxplot) lim2(icol,:) = 0

 return
end subroutine reset_lim2


!----------------------------------------------------------
! function which returns whether or not a range
! has been set for any column
!----------------------------------------------------------
logical function anyrangeset()
 use settings_data, only:ndataplots
 integer :: i

 anyrangeset = .false.
 do i=1,ndataplots
    if (rangeset(i)) anyrangeset = .true.
 enddo

 return
end function anyrangeset

!----------------------------------------------------------
! prints info about current range restriction settings
!----------------------------------------------------------
subroutine print_rangeinfo()
 use settings_data, only:ndataplots
 use labels,        only:label
 integer :: i

 if (anyrangeset()) then
    print "(/,a,/)",'>> current range restrictions set: '
    do i=1,ndataplots
       if (rangeset(i)) then
          print "(a,1pe10.3,a,1pe10.3,a)", &
          ' ( ',range(i,1),' < '//trim(label(i))//' < ',range(i,2),' )'
       endif
    enddo
    print "(/,2(a,/))",'>> only particles within this range will be plotted ', &
                     '   and/or used in interpolation routines'
 !else
    !print "(/,a,/)",'>> no current parameter range restrictions set '
 endif

end subroutine print_rangeinfo

!----------------------------------------------------------
! prints info about current range restriction settings
!----------------------------------------------------------
subroutine print_lim2info()
 use settings_data, only:ndataplots
 use labels,        only:label
 integer :: i

 do i=1,ndataplots
    if (lim2set(i)) then
       print "(a,1pe10.3,a,1pe10.3,a)", &
       ' ( contours use ',lim2(i,1),' < '//trim(label(i))//' < ',lim2(i,2),' )'
    endif
 enddo

end subroutine print_lim2info

!----------------------------------------------------------
! prints warning if min=max in limits setting
!----------------------------------------------------------
subroutine warn_minmax(labelx,xmin,xmax)
 character(len=*), intent(in) :: labelx
 real,             intent(in) :: xmin,xmax

 if (abs(xmin-xmax).lt.tiny(xmax)) then
    print "(a,a20,a,1pe9.2)",'  warning: ',labelx,' min = max = ',xmin
 endif

 return
end subroutine warn_minmax

!----------------------------------------------------------
! Makes sure that variable is within a given range
! If no range specified, ensures that it is within
! the allowed range for the variable type,
!  i.e. -0.5*huge(x)->0.5*huge(x)
!----------------------------------------------------------
subroutine assert_range(x,min,max)
 real, intent(inout) :: x
 real, intent(in), optional :: min,max
 real :: xmin,xmax
 
 xmin = -0.5*huge(xmin)  ! for limits need xmax - xmin to
 xmax = 0.5*huge(xmax)   ! be less than huge(x)
 if (present(min)) xmin = min
 if (present(max)) xmax = max
 if (x < xmin) x = xmin
 if (x > xmax) x = xmax
 if (x /= x) x = 0.

 return
end subroutine assert_range

!----------------------------------------------------------
! Interface to the above, but checks two numbers at once
! and checks that max > min
!----------------------------------------------------------
subroutine assert_sensible_limits(xmin,xmax)
 real, intent(inout) :: xmin,xmax
 real :: xtmp

 call assert_range(xmin)
 call assert_range(xmax)
 
 if (xmax < xmin) then
    xtmp = xmin
    xmin = xmax
    xmax = xtmp
 endif

 return
end subroutine assert_sensible_limits

end module limits
