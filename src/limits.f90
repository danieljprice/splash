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
 private :: warn_minmax,check_coord_column,set_coord_limits,&
            read_limits_file,override_limits_from_flags

 public

contains

!----------------------------------------------------------
! set plot limits for all columns
! NB: does not differentiate between particle types at the moment
!----------------------------------------------------------
subroutine set_limits(ifromstep,itostep,ifromcol,itocol,use_type)
 use labels,        only:label,ix
 use geometry,      only:coord_transform_limits
 use particle_data, only:npartoftype,dat,maxcol,iamtype
 use settings_data, only:ndim,icoords,icoordsnew,iverbose,debugmode
 integer, intent(in) :: ifromstep,itostep,ifromcol,itocol
 logical, intent(in), optional :: use_type(:)
 integer :: i,j,k,ntoti,itocoli
 logical :: first

 if (iverbose > 1 .or. debugmode) print 100,ifromstep,itostep,ifromcol,itocol
100 format(/' setting plot limits: steps ',i5,'->',i5,' cols ',i2,'->',i3)
 !--avoid tripping over array bounds (but these are valid do-nothing calls)
 if (ifromcol > maxcol .or. maxcol==0) return
 if (ifromcol > itocol) return
 itocoli = min(itocol,maxcol)

 !--find limits of particle properties
 lim(ifromcol:itocol,1) = huge(lim)
 lim(ifromcol:itocol,2) = -huge(lim)
 do i=ifromstep,itostep
    ntoti = sum(npartoftype(:,i))
    if (allocated(iamtype) .and. present(use_type) .and. size(iamtype(:,1)) > 1) then
       do j=ifromcol,itocoli
          do k=1,ntoti
             if (use_type(iamtype(k,i))) then
                lim(j,1) = min(lim(j,1),dat(k,j,i))
                lim(j,2) = max(lim(j,2),dat(k,j,i))
             endif
          enddo
       enddo
    else
       do j=ifromcol,itocoli
          do k=1,ntoti
             lim(j,1) = min(lim(j,1),dat(k,j,i))
             lim(j,2) = max(lim(j,2),dat(k,j,i))
          enddo
       enddo
    endif
 enddo
 !
 !--warn if limits are the same
 !
 first = .true.
 do j=ifromcol,itocol
    if (iverbose >= 0) call warn_minmax(label(j),lim(j,1),lim(j,2),first)
 enddo

 lim2(ifromcol:itocol,:) = 0.

 !
 !--transform coord limits into new coordinate system if coord transform is applied
 !
 if (icoordsnew /= icoords .and. ndim > 0) then
    if (ifromcol <= ix(ndim)) then ! separate if is to avoid referencing ix(0) if ndim=0
       call coord_transform_limits(lim(ix(1):ix(ndim),1),lim(ix(1):ix(ndim),2), &
                                    icoords,icoordsnew,ndim)
    endif
 endif

end subroutine set_limits

!----------------------------------------------------------
! rescale plot limits for all columns by factor
! e.g. if units have changed
!----------------------------------------------------------
subroutine rescale_limits(fac)
 real, intent(in) :: fac(:)
 integer :: j

 do j=1,min(size(fac),size(lim),size(lim2))
    lim(j,:) = lim(j,:)*fac(j)
    lim2(j,:) = lim2(j,:)*fac(j)
 enddo

end subroutine rescale_limits

!----------------------------------------------------------
! save plot limits for all columns to a file
!----------------------------------------------------------
subroutine write_limits(limitsfile,iverbose)
 use settings_data, only:numplot,ndataplots
 character(len=*), intent(in) :: limitsfile
 integer, intent(in) :: iverbose
 integer :: i

 if (iverbose >= 0) print*,'saving plot limits to file ',trim(limitsfile)

 open(unit=55,file=limitsfile,status='replace',form='formatted',ERR=998)
 do i=1,numplot
    if (rangeset(i) .and. i < ndataplots .and. lim2set(i)) then
       write(55,"(6(1x,1pe14.6))",err=999) lim(i,1),lim(i,2),range(i,1),range(i,2),lim2(i,1),lim2(i,2)
    elseif (lim2set(i) .and. i < ndataplots) then
       write(55,"(6(1x,1pe14.6))",err=999) lim(i,1),lim(i,2),0.,0.,lim2(i,1),lim2(i,2)
    elseif (rangeset(i) .and. i < ndataplots) then
       write(55,"(4(1x,1pe14.6))",err=999) lim(i,1),lim(i,2),range(i,1),range(i,2)
    else
       write(55,"(2(1x,1pe14.6))",err=999) lim(i,1),lim(i,2)
    endif
 enddo
 close(unit=55)

 return

998 continue
 if (iverbose >= 0) print*,'*** error opening limits file: limits not saved'
 return
999 continue
 if (iverbose >= 0) print*,'*** error saving limits'
 close(unit=55)

end subroutine write_limits

!----------------------------------------------------------
! read plot limits from file, then apply command-line overrides
!----------------------------------------------------------
subroutine read_limits(limitsfile,ierr,overridden)
 character(len=*), intent(in) :: limitsfile
 integer,         intent(out) :: ierr
 logical, intent(out), optional :: overridden

 call read_limits_file(limitsfile,ierr)
 call override_limits_from_flags(overridden)

end subroutine read_limits

!----------------------------------------------------------
! read plot limits for all columns from a file
!----------------------------------------------------------
subroutine read_limits_file(limitsfile,ierr)
 use labels,        only:label
 use settings_data, only:numplot,ncolumns,ncalc,iverbose
 use asciiutils,    only:ncolumnsline
 character(len=*), intent(in) :: limitsfile
 integer,         intent(out) :: ierr
 integer                      :: i,ncolsline
 character(len=120)           :: line
 logical :: iexist

 ierr = 0

 inquire(file=limitsfile,exist=iexist)
 if (.not.iexist) then
    !print "(1x,a)",trim(limitsfile)//' not found'
    ierr = 1
    return
 endif

 open(unit=54,file=limitsfile,status='old',form='formatted',err=997)
 if (iverbose >= 0) print "(a)",' read '//trim(limitsfile)
 do i=1,numplot
    read(54,"(a)",err=998,end=999) line
    ncolsline = ncolumnsline(line)
    if (ncolsline < 2) then
       goto 998
    elseif (ncolsline >= 6) then
       read(line,*,err=998,end=999) lim(i,1),lim(i,2),range(i,1),range(i,2),lim2(i,1),lim2(i,2)
    elseif (ncolsline >= 4) then
       read(line,*,err=998,end=999) lim(i,1),lim(i,2),range(i,1),range(i,2)
    else
       read(line,*,err=998,end=999) lim(i,1),lim(i,2)
    endif
    call assert_sensible_limits(lim(i,1),lim(i,2))
    !
    !--warn if limits are the same
    !
    if (iverbose >= 0) call warn_minmax(label(i),lim(i,1),lim(i,2))
 enddo
 close(unit=54)
 return

997 continue
 if (iverbose >= 0) print*,trim(limitsfile),' not found'
 ierr = 1
 return
998 continue
 if (iverbose >= 0) then
    call print_rangeinfo()
    call print_lim2info()
    print*,'*** error reading limits from file'
 endif
 ierr = 2
 close(unit=54)
 return
999 continue
 !--only give error if we really do not have enough columns
 !  (on first call nextra is not set)
 if (i < ncolumns+ncalc) then
    if (iverbose >= 0) print "(a,i3)",' end of file in '//trim(limitsfile)//': limits read to column ',i
    ierr = -1
 endif

 !--print info about range restrictions read from file
 if (iverbose >= 0) then
    call print_rangeinfo()
    call print_lim2info()
 endif
 close(unit=54)

end subroutine read_limits_file

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
    where (icolours(:)==-1000)
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
          where (datstep(:,icol) < range(icol,1) .or. &
                 datstep(:,icol) > range(icol,2))
             where (icolours==0)
                icolours = -1000
             elsewhere
                icolours = -abs(icolours)
             end where
          end where
       endif
    enddo
 endif

end subroutine get_particle_subset

!----------------------------------------------------------
! reset all range restrictions to zero
!----------------------------------------------------------
subroutine reset_all_ranges()
 use particle_data, only:icolourme

 print "(a)",' removing all range restrictions '
 where (icolourme(:)==-1000)
    icolourme(:) = 0
 elsewhere
    icolourme(:) = abs(icolourme(:))
 endwhere
 range(:,:) = 0.

end subroutine reset_all_ranges

!----------------------------------------------------------
! function which returns whether or not a range
! has been set for a given column
!----------------------------------------------------------
logical function rangeset(icol)
 integer, intent(in) :: icol

 rangeset = .false.
 if (abs(range(icol,2)-range(icol,1)) > tiny(range)) rangeset = .true.

end function rangeset

!----------------------------------------------------------
! function which returns whether or not lim2
! has been set for a given column
!----------------------------------------------------------
logical function lim2set(icol)
 integer, intent(in) :: icol

 lim2set = .false.
 if (abs(lim2(icol,2)) > tiny(lim2) .or. abs(lim2(icol,1)) > tiny(lim2)) lim2set = .true.

end function lim2set

!----------------------------------------------------------
! reset all range restrictions to zero
!----------------------------------------------------------
subroutine reset_lim2(icol)
 integer, intent(in) :: icol

 print "(a)",' contour limits same as render limits'
 if (icol > 0 .and. icol <= maxplot) lim2(icol,:) = 0

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
subroutine warn_minmax(labelx,xmin,xmax,first)
 character(len=*), intent(in) :: labelx
 real,             intent(in) :: xmin,xmax
 logical,          intent(inout), optional :: first

 if (abs(xmin-xmax) < tiny(xmax)) then
    if (present(first)) then
       if (first) print "(a)"
       first = .false.
    endif
    print "(a,a20,a,1pe9.2)",'  warning: ',labelx,' min = max = ',xmin
 endif

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
 if (isnan(x)) x = 0.

end subroutine assert_range

!----------------------------------------------------------
! Interface to the above, but checks two numbers at once
! and checks that max > min
!----------------------------------------------------------
subroutine assert_sensible_limits(xmin,xmax)
 real, intent(inout) :: xmin,xmax

 call assert_range(xmin)
 call assert_range(xmax)

end subroutine assert_sensible_limits

!----------------------------------------------------------
! Interface to the above, but checks two numbers at once
! and checks that max > min
!----------------------------------------------------------
subroutine fix_equal_limits(xmin,xmax)
 real, intent(inout) :: xmin,xmax

 call assert_sensible_limits(xmin,xmax)
 if (abs(xmax - xmin) < tiny(xmin)) then
    xmax = xmax + 1.0
    if (xmin > 0.) then
       xmin = max(xmin - 1.0,0.)
    else
       xmin = xmin - 1.0
    endif
 endif

end subroutine fix_equal_limits

!----------------------------------------------------------
! checks whether all the min and max limits in a
! range of columns are the same
!----------------------------------------------------------
logical function limits_are_equal(n,iplotx,iploty)
 integer, intent(in) :: n
 integer, intent(in), dimension(n) :: iplotx,iploty
 real :: xlim(2), ylim(2)
 integer :: j

 limits_are_equal = .true.
 if (n > 0) then
    xlim = lim(iplotx(1),:)
    ylim = lim(iploty(1),:)
 endif
 do j=1,2
    if (any((lim(iplotx(1:n),j) - xlim(j)) > epsilon(0.))) limits_are_equal = .false.
    if (any((lim(iploty(1:n),j) - ylim(j)) > epsilon(0.))) limits_are_equal = .false.
 enddo

end function limits_are_equal

!----------------------------------------------------------
! check that coordinate axis idim maps to a valid plot column
!----------------------------------------------------------
subroutine check_coord_column(idim,flag)
 use labels,        only:ix
 use settings_data, only:ndim,numplot
 integer, intent(in) :: idim
 character(len=*), intent(in) :: flag

 if (idim < 1 .or. idim > ndim .or. ix(idim) <= 0 .or. ix(idim) > numplot) then
    print "(a)",' ERROR: '//trim(flag)//' given but coordinate column is not set'
    stop
 endif

end subroutine check_coord_column

!----------------------------------------------------------
! set lim for coordinate axis idim after validating ix(idim)
!----------------------------------------------------------
subroutine set_coord_limits(idim,xlo,xhi,flag)
 use labels, only:ix
 integer, intent(in) :: idim
 real, intent(in) :: xlo,xhi
 character(len=*), intent(in) :: flag

 call check_coord_column(idim,flag)
 lim(ix(idim),1:2) = (/xlo,xhi/)

end subroutine set_coord_limits

!----------------------------------------------------------
! apply --xmin/--limits/--lim overrides (mutually exclusive)
!----------------------------------------------------------
subroutine override_limits_from_flags(overridden)
 use labels,        only:ix
 use geometry,      only:labelcoord
 use settings_data, only:ndim,numplot
 use system_utils,  only:get_command_flag,get_command_option,rflaglist
 integer, parameter :: maxvals = 2*maxplot
 real :: vals(maxvals),xmin,xmax,halfwidth
 logical, intent(out), optional :: overridden
 logical :: have_axis,have_limits,have_lim,have_min(3),have_max(3)
 integer :: i,nvals,nfam,icol

 do i=1,3
    have_min(i) = get_command_flag(trim(labelcoord(i,1))//'min')
    have_max(i) = get_command_flag(trim(labelcoord(i,1))//'max')
 enddo
 have_axis = any(have_min .or. have_max)
 have_limits = get_command_flag('limits')
 have_lim    = get_command_flag('lim')

 nfam = count((/have_axis,have_limits,have_lim/))
 if (present(overridden)) overridden = (nfam > 0)
 if (nfam > 1) then
    print "(a)",' ERROR: use only one of --xmin/--xmax/... or --limits=... or --lim=...'
    stop
 elseif (nfam == 0) then
    return
 endif

 if (have_lim) then
    !
    ! set x, y and z limits from --lim=halfwidth flag
    !
    halfwidth = get_command_option('lim',default=-1.)
    if (halfwidth < 0. .or. ndim < 1) then
       print "(a)",' ERROR: --lim=halfwidth requires identified coordinate columns'
       stop
    endif
    do i=1,ndim
       call set_coord_limits(i,-halfwidth,halfwidth,'--lim')
    enddo
    print "(2(a,1pg10.3),a)",' plot limits: centred box [',-halfwidth,',',halfwidth,']'
 elseif (have_limits) then
    !
    ! set individual column limits from --limits=[min,max,min,max,...]
    !
    vals = rflaglist('limits',maxvals,ngot=nvals)
    if (nvals < 2 .or. mod(nvals,2) /= 0) then
       print "(a)",' ERROR: --limits requires an even number of comma-separated values'
       stop
    endif
    if (nvals == 2 .and. ndim >= 1) then
       !
       ! if we have exactly 2 values then apply to all coordinate axes (like --lim)
       !
       do i=1,ndim
          call set_coord_limits(i,vals(1),vals(2),'--limits')
       enddo
       print "(a,2(1x,1pg10.3))",' plot limits: coordinate axes set to ',vals(1:2)
    elseif (nvals == 2*ndim .and. ndim >= 1) then
       ! 
       ! if we have exactly 2*ndim values then these
       ! are assumed to apply to coordinate columns
       !
       do i=1,ndim
          call set_coord_limits(i,vals(2*i-1),vals(2*i),'--limits')
       enddo
       print "(a,i2,a)",' plot limits: set from --limits for ',ndim,' axes'
    else
       !
       ! otherwise apply to all columns
       !
       do icol=1,nvals/2
          if (icol > numplot) then
             print "(a,i3)",' ERROR: --limits column index exceeds numplot: ',icol
             stop
          endif
          lim(icol,1:2) = vals(2*icol-1:2*icol)
       enddo
       print "(a,i3,a)",' plot limits: overrode ',nvals/2,' columns from --limits'
    endif
 else
    !
    ! set x,y and z limits from --xmin/--xmax flags
    !
    do i=1,ndim
       if (.not. have_min(i) .and. .not. have_max(i)) cycle
       call check_coord_column(i,'--'//trim(labelcoord(i,1))//'min/max')
       xmin = get_command_option(trim(labelcoord(i,1))//'min',default=lim(ix(i),1))
       xmax = get_command_option(trim(labelcoord(i,1))//'max',default=lim(ix(i),2))
       lim(ix(i),1:2) = (/xmin,xmax/)
       print "(a,2(1x,1pg10.3))",' plot limits: '//trim(labelcoord(i,1))//' set to',xmin,xmax
    enddo
 endif

end subroutine override_limits_from_flags

end module limits
