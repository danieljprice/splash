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
!  Copyright (C) 2005-2011 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------
!
! This module allows arbitrary functions to be computed from the
! particle data. Uses the fparser module to parse the functions.
!
!-----------------------------------------------------------------
module calcquantities
 use labels, only:lenlabel,lenunitslabel
 use params, only:maxplot
 implicit none
 public :: calc_quantities,setup_calculated_quantities
 public :: calc_quantities_use_x0

 integer, parameter, private :: maxcalc = 35
 character(len=lenlabel),      dimension(maxcalc) :: calcstring = ' '
 character(len=lenlabel),      dimension(maxcalc) :: calclabel = ' '
 character(len=lenunitslabel), dimension(maxcalc) :: calcunitslabel = ' '

 integer, parameter, private :: lenvars = 20 ! max length for any variable
 integer, parameter, private :: nextravars = 5
 character(len=lenvars), dimension(nextravars), parameter, private :: &
          extravars=(/'t    ','gamma','x0   ','y0   ','z0   '/)

 namelist /calcopts/ calcstring,calclabel,calcunitslabel
 
 public :: calcopts,calcstring,calclabel,calcunitslabel
 private
 
contains

!-----------------------------------------------------------------
!
!  Allow the user to edit the list of function strings used
!  to compute additional quantities from the particle data
!
!-----------------------------------------------------------------
subroutine setup_calculated_quantities(ncalc)
  use settings_data, only:ncolumns
 use prompting,     only:prompt
 implicit none
 integer, intent(out) :: ncalc
 integer              :: ncalctot,istart,iend, ipick, ninactive, i
 logical              :: done,first
 character(len=1)     :: charp
 integer, dimension(maxcalc) :: incolumn
 ipick = ncolumns + 1

 done = .false.
 first = .true.
 charp = 'a'
 calcmenu: do while (.not.done)
    call check_calculated_quantities(ncalc,ncalctot,incolumn)
    ninactive = ncalctot - ncalc   

    iend = maxcalc
    if (ncalctot.gt.0 .or. .not.first) then
       charp='a'
       if (.not.first) charp = 'q'
       print*
       call prompt(' a)dd to, e)dit, c)lear current list or q)uit/finish? ',&
                   charp,list=(/'a','e','c','q','s','S','Q'/),noblank=.true.)
       select case(charp)
       case('a')
          istart = ncalctot
          iend   = ncalctot+1
       case('e')
          if (ninactive.gt.0) then
             call prompt(' pick a function to edit ',ipick,-ninactive,-1,ncolumns+1,ncolumns+ncalc)
          else
             call prompt(' pick a function to edit ',ipick,ncolumns+1,ncolumns+ncalc)
          endif
          istart = 0
          do i=1,ncalctot
             if (incolumn(i).eq.ipick) istart = i - 1
          enddo
          iend = istart + 1
       case('c')
           calcstring(:) = ' '
           calclabel(:) = ' '
           calcunitslabel(:) = ' '
           cycle calcmenu
       case('q','Q','s','S')
          done = .true.
       case default
          istart = 0
          iend   = maxcalc
       end select
    else
       istart = 0
       iend   = 1
    endif
    
    if (.not.done) call add_calculated_quantities(istart,iend,ncalc,first,incolumn)
    first = .false.
 enddo calcmenu
 if (ncalc.lt.10) then
    print "(a,i1,a)",' setup ',ncalc,' additional quantities'
 else
    print "(a,i2,a)",' setup ',ncalc,' additional quantities'
 endif

end subroutine setup_calculated_quantities

!-----------------------------------------------------------------
!
!  utility (private) to add one or more calculated quantities
!  to the current list and/or edit previous settings
!
!-----------------------------------------------------------------
subroutine add_calculated_quantities(istart,iend,ncalc,printhelp,incolumn)
 use prompting,     only:prompt
 use fparser,       only:checkf
 use labels,        only:lenlabel
 use settings_data, only:ncolumns,iRescale
 implicit none
 integer, intent(in)  :: istart,iend
 integer, intent(out) :: ncalc
 logical, intent(in)  :: printhelp
 integer, dimension(maxcalc), intent(in) :: incolumn
 integer :: i,j,ntries,ierr,iequal,nvars
 logical :: iask
 character(len=120) :: string
 character(len=lenvars), dimension(maxplot+nextravars) :: vars

 i = istart + 1
 ntries = 0
 ncalc = istart
 if (i.gt.maxcalc) then
    print "(/,a,i2,a)",' *** Error, maximum number of calculated quantities (',maxcalc,') reached, cannot add any more.'
    print "(a)",       ' *** If you hit this limit, *please email me* so I can change the default limits!'
    print "(a)",       ' *** (and then edit calc_quantities.f90, changing the parameter "maxcalc" to something higher...)'
    return
 endif 

 if (printhelp) then
    print "(/,a)",' Specify a function to calculate from the data '
    print "(10(a))",' Valid variables are the column labels',(', '''//trim(extravars(j))//'''',j=1,nextravars-1),&
                ' and '''//trim(extravars(nextravars))//''' (origin setting) '
    print "(a)",' Spaces, escape sequences (\d) and units labels are removed from variable names'
    print "(a)",' Note that previously calculated quantities can be used in subsequent calculations'
 endif
 call print_example_quantities()
  
 overfuncs: do while(ntries.lt.3 .and. i.le.iend .and. i.le.maxcalc)
    if (len_trim(calcstring(i)).ne.0 .or. ncalc.gt.istart) then
       if (incolumn(i).gt.0) then
          write(*,"(a,i2,a)") '[Column ',incolumn(i),']'       
       else
          write(*,"(a)") '[Currently inactive]'
       endif
    endif

    if (len_trim(calclabel(i)).gt.0) then
       string = trim(calclabel(i))//' = '//trim(calcstring(i))    
    else
       string = trim(calcstring(i))
    endif
    call prompt('Enter function string to calculate (blank for none) ',string)
    if (len_trim(string).eq.0) then
       !
       !--if editing list and get blank string at prompt,
       !  remove the entry and shuffle the list appropriately
       !
       do j=i,maxcalc-1
          !print*,' shuffling ',j,' = '//trim(calcstring(j+1))
          calcstring(j) = trim(calcstring(j+1))
          if (j+1.lt.size(calclabel)) then
             calclabel(j) = calclabel(j+1)
          endif
          if (len_trim(calcstring(j)).eq.0) then
             if (j.eq.i) then
                exit overfuncs
             else
                cycle overfuncs
             endif
          endif
       enddo
    else
       !
       !--set label from what lies to the left of the equal sign
       !  and the calcstring from what lies to the right
       !  (if no equals sign, then prompt later for the label)
       !
       iequal = index(string,'=')
       if (iequal.ne.0) then
          calclabel(i)  = string(1:iequal-1)
          calcstring(i) = string(iequal+1:len_trim(string))
       else
          calcstring(i) = trim(string)
       endif

       !--remove preceding spaces
       calcstring(i) = trim(adjustl(calcstring(i)))
       !
       !--fill variable array with the list of valid variable
       !  names for this column
       !
       call get_variables(ncolumns+i-1,nvars,vars)
       !
       !--check for errors parsing function
       !
       ierr = checkf(shortlabel(calcstring(i)),vars(1:nvars))
       if (ierr.ne.0 ) then
          ntries = ntries + 1
          print "(a,i1,a)",' error parsing function string: try again (',ntries,' of 3)'
          if (ntries.eq.3) then
             iask = .false.
             call prompt(' Cannot parse function (after 3 attempts). Set as inactive function anyway?',iask)
             if (iask) then
                ierr = 0
             else
                calcstring(i) = ' '
             endif
          endif
       else
          write(*,"(a)",advance='no') 'Function parses OK: '
       endif
       if (ierr.eq.0) then
          !
          !--prompt for label if not set
          !
          if (iequal.eq.0) then
             call prompt(' Enter label for this quantity ',calclabel(i),noblank=.true.)
          endif
          if (iRescale) then
             call prompt(' Enter units label for this quantity ',calcunitslabel(i))
          endif
          !print "(a,a,i2,/)",'Setting '//trim(calclabel(i)), &
          !                   ' = '//trim(calcstring(i))//' in column ',ncolumns+i
          ncalc = i
          i = i + 1
       endif
    endif
 enddo overfuncs

end subroutine add_calculated_quantities

!---------------------------------------------------------------------
!
!  utility to give a nice list of examples to follow / cut and paste
!  [ this basically replaces what was hardwired into the
!    old calc_quantities routine ]
!
!---------------------------------------------------------------------
subroutine print_example_quantities
 use labels,        only:label,lenlabel,irho,iutherm,iBfirst,ix,icv,iradenergy,iamvec,labelvec
 use settings_data, only:ncolumns,ndim,icoordsnew,ndimV
 use settings_units,only:unitslabel
 use geometry,      only:labelcoord
 implicit none
 integer :: i,j,ivecstart

 print "(/,a)",' Examples based on current data: '
 !--radius
 if (ndim.gt.0 .and. icoordsnew.eq.1 .and. ncolumns.ge.ndim) then
    write(*,"(11x,a)",ADVANCE='NO') 'r = sqrt(('// &
          trim(shortlabel(label(ix(1)),unitslabel(ix(1))))//'-'//trim(labelcoord(1,1))//'0)**2'
    if (ndim.gt.1) then
       write(*,"(a,a,a)") (' + ('//trim(shortlabel(label(ix(i)),unitslabel(ix(i))))// &
                                                 '-'//trim(labelcoord(i,1))//'0)**2',i=2,ndim),')'
    else
       write(*,"(a)") ')'
    endif
 elseif (ncolumns.ge.2) then
 !--if ndim=0 give random example to give at least one
    print "(11x,a)",trim(shortlabel(label(1)))//'*'//trim(shortlabel(label(2)))
 endif
 !--pressure
 if (irho.gt.0 .and. iutherm.gt.0) then
    print "(11x,a)",'pressure = (gamma-1)*'//trim(shortlabel(label(irho),unitslabel(irho)))// &
                '*'//trim(shortlabel(label(iutherm),unitslabel(iutherm)))
 endif
 !
 !--magnitudes of all vector quantities (only if cartesian coords are set)
 !
 ivecstart = 0
 if (icoordsnew.eq.1 .and. ndim.gt.0 .and. ndimV.gt.0) then
    do i=1,ncolumns
       if (iamvec(i).gt.0 .and. iamvec(i).le.ncolumns .and. iamvec(i).ne.ivecstart) then
          ivecstart = iamvec(i)
          write(*,"(11x,a)",ADVANCE='NO') '|'//trim(labelvec(ivecstart))//'| '// &
            '= sqrt('//trim(shortlabel(label(ivecstart),unitslabel(ivecstart)))//'**2'
          if (ndimV.gt.1) then
             write(*,"(a,a,a)") (' + '//trim(shortlabel(label(j),unitslabel(j)))//'**2',&
                                 j=ivecstart+1,ivecstart+ndimV-1),')'
          else
             write(*,"(a)") ')'
          endif
       endif
    enddo
 endif
 !--magnetic pressure
 if (ndim.gt.0 .and. ndimV.gt.0 .and. iBfirst.gt.0 .and. icoordsnew.eq.1) then
    write(*,"(11x,a)",ADVANCE='NO') 'P\dmag = 0.5*('//trim(shortlabel(label(iBfirst),unitslabel(iBfirst)))//'**2'
    if (ndimV.gt.1) then
       write(*,"(a,a,a)") (' + '//trim(shortlabel(label(i),unitslabel(i)))//'**2',i=iBfirst+1,iBfirst+ndimV-1),')'
    else
       write(*,"(a)") ')'
    endif
 endif
 !--gas temperature if cv present
 if (ndim.gt.0 .and. iutherm.gt.0 .and. icv.gt.0) then
    print "(6x,a)",'T\dgas\u = '//trim(shortlabel(label(iutherm),unitslabel(iutherm)))//'/' &
                    //trim(shortlabel(label(icv),unitslabel(icv)))
 endif
 !--radiation temperature
 if (ndim.gt.0 .and. irho.gt.0 .and. iradenergy.gt.0) then
    print "(6x,a)",'T\drad\u = ('//trim(shortlabel(label(irho),unitslabel(irho)))//'*' &
                    //trim(shortlabel(label(iradenergy),unitslabel(iradenergy)))//'/7.5646e-15)**0.25)'
 endif
 print "(a)"

end subroutine print_example_quantities

!-----------------------------------------------------------------
!
!  utility (private) to print the current list of calculated
!  quantities, checking that they parse correctly
!
!-----------------------------------------------------------------
subroutine check_calculated_quantities(ncalcok,ncalctot,incolumn)
 use settings_data,  only:ncolumns,iRescale
 use fparser,        only:checkf
 use labels,         only:label
 use settings_units, only:unitslabel
 implicit none
 integer, intent(out) :: ncalcok,ncalctot
 integer, dimension(maxcalc), intent(out), optional :: incolumn
 integer :: i,ierr,nvars,indexinactive
 character(len=lenvars), dimension(maxplot+nextravars) :: vars

 ncalcok = 0
 ncalctot = 0
 indexinactive = 0
 i = 1
 if (present(incolumn)) incolumn(:) = 0
 print "(/,a)", ' Current list of calculated quantities:'
 do while(i.le.maxcalc .and. len_trim(calcstring(i)).ne.0)
    !
    !--get the list of valid variable names for this column
    !
    call get_variables(ncolumns+ncalcok,nvars,vars)
    !
    !--check that the function parses
    !
    ierr = checkf(shortlabel(calcstring(i)),vars(1:nvars),Verbose=.false.)

    if (ierr.eq.0) then
       ncalcok = ncalcok + 1
       print "(1x,i2,') ',a50,' [OK]')",ncolumns+ncalcok,trim(calclabel(i))//' = '//calcstring(i)
       if (present(incolumn)) incolumn(i) = ncolumns + ncalcok
       !
       !--set the label for the proposed column here
       !  so that subsequent calculations can use this variable
       !  (note that we don't need to set the units label here 
       !   as this is done in the actual calc_quantities call)
       !
       label(ncolumns+ncalcok) = trim(calclabel(i))
       unitslabel(ncolumns+ncalcok) = trim(calcunitslabel(i))
       if (iRescale) label(ncolumns+ncalcok) = trim(label(ncolumns+ncalcok))//trim(unitslabel(ncolumns+ncalcok))
    else
       indexinactive = indexinactive - 1
       print "(i3') ',a50,' [INACTIVE]')",indexinactive,trim(calclabel(i))//' = '//calcstring(i)
       if (present(incolumn)) incolumn(i) = indexinactive
    endif
    ncalctot = i
    i = i + 1
 enddo
 if (ncalcok.eq.0) print "(a)",' (none)'

end subroutine check_calculated_quantities

!-----------------------------------------------------------------
!
!  actually compute the extra quantities from the particle data
!
!-----------------------------------------------------------------
subroutine calc_quantities(ifromstep,itostep,dontcalculate)
  use labels,         only:label,labelvec,iamvec,ix
  use particle_data,  only:dat,npartoftype,gamma,time,maxpart,maxstep,maxcol
  use settings_data,  only:ncolumns,ncalc,iRescale,xorigin,debugmode,itrackpart,ndim
  use mem_allocation, only:alloc
  use settings_units, only:unitslabel,units
  use fparser,        only:checkf,parsef,evalf,EvalerrMsg,EvalErrType,rn,initf,endf
  use params,         only:maxplot
  implicit none
  integer, intent(in) :: ifromstep, itostep
  logical, intent(in), optional :: dontcalculate
  integer :: i,j,ncolsnew,ierr,icalc,ntoti,nvars,ncalctot
  logical :: skip
!  real, parameter :: mhonkb = 1.6733e-24/1.38e-16
!  real, parameter :: radconst = 7.5646e-15
!  real, parameter :: lightspeed = 3.e10   ! in cm/s (cgs)
  real(kind=rn), dimension(maxplot+nextravars)          :: vals
  character(len=lenvars), dimension(maxplot+nextravars) :: vars
  real, dimension(3) :: x0
  
  !
  !--allow dummy call to set labels without actually calculating stuff
  !
  if (present(dontcalculate)) then
     skip = dontcalculate
  else
     skip = .false.
  endif
  
  ierr = 0
  ncalc = 0
  call check_calculated_quantities(ncalc,ncalctot)
  
  if (.not.skip) print "(2(a,i2),a,/)",' Calculating ',ncalc,' of ',ncalctot,' additional quantities...'
  ncolsnew = ncolumns + ncalc
  if (ncolsnew.gt.maxcol) call alloc(maxpart,maxstep,ncolsnew) 

  !
  !--reset iamvec to zero for calculated columns
  !
  iamvec(ncolumns+1:ncolsnew) = 0
  labelvec(ncolumns+1:ncolsnew) = ' '

  !
  !--evaluate functions in turn
  !
  if (.not.skip .and. ncalc.gt.0) then
     call initf(ncalc)
     !
     !--compile each function into bytecode
     !
     icalc = 1
     do i=1,maxcalc
        if (icalc.le.ncalc) then

           !
           !--get the list of valid variable names for this column
           !
           call get_variables(ncolumns+icalc-1,nvars,vars)
           !
           !--now actually parse the function
           !
           call parsef(icalc,shortlabel(calcstring(i)),vars(1:nvars),err=ierr,Verbose=.false.)
           if (ierr.eq.0) then
              icalc = icalc + 1
           endif
        endif
     enddo
     !
     !--evaluate functions from particle data
     !
     do i=ifromstep,itostep
        ntoti = SUM(npartoftype(:,i))
        !
        !--set origin position
        !
        if (itrackpart.gt.0 .and. itrackpart.le.ntoti) then
           x0(:) = 0.
           if (ix(1).gt.0 .and. ix(1).le.ncolumns) then
              x0(1) = dat(itrackpart,ix(1),i)
           else
              print*,'** internal error: tracking particle set but cannot locate x coordinate'
           endif
           if (ix(2).gt.0 .and. ix(2).le.ncolumns) x0(2) = dat(itrackpart,ix(2),i)
           if (ix(3).gt.0 .and. ix(3).le.ncolumns) x0(3) = dat(itrackpart,ix(3),i)
           if (i.eq.ifromstep) then
              print "(a,i10)",' using position of tracked particle ',itrackpart
              print "(a,3(e11.3),/)",' (x0,y0,z0) = ',dat(itrackpart,ix(1:ndim),i)
           endif
        else
           x0(:) = xorigin(:)
        endif
        
        do icalc=1,ncalc
           if (debugmode) print*,'DEBUG: ',icalc,' calculating '//trim(label(ncolumns+icalc))
           !
           !--additional settings allowed in function evaluations
           !  i.e., time and gamma from dump file and current origin settings
           !  make sure the number here aligns with the "nextravars" setting
           !
           vals(ncolumns+icalc)   = time(i)
           vals(ncolumns+icalc+1) = gamma(i)
           vals(ncolumns+icalc+2) = x0(1)
           vals(ncolumns+icalc+3) = x0(2)
           vals(ncolumns+icalc+4) = x0(3)
           do j=1,ntoti
              vals(1:ncolumns+icalc-1) = dat(j,1:ncolumns+icalc-1,i)
              dat(j,ncolumns+icalc,i) = real(evalf(icalc,vals(1:ncolumns+icalc+nextravars-1)))
           enddo
           if (EvalErrType.ne.0) then
              print "(a)",' ERRORS evaluating '//trim(calcstring(icalc))//': ' &
                          //trim(EvalerrMsg())
           endif
        enddo
     enddo
     call endf
  endif
 
  !
  !--override units of calculated quantities if necessary
  !
  if (iRescale .and. any(abs(units(ncolumns+1:ncolumns+ncalc)-1.0).gt.tiny(0.)) &
      .and. .not.skip) then
     write(*,"(/a)") ' rescaling data...'
     do i=ncolumns+1,ncolumns+ncalc
        if (abs(units(i)-1.0).gt.tiny(0.) .and. abs(units(i)).gt.tiny(0.)) then
           dat(:,i,ifromstep:itostep) = dat(:,i,ifromstep:itostep)*units(i)
        endif
        if (index(label(i),trim(unitslabel(i))).eq.0) label(i) = trim(label(i))//trim(unitslabel(i))
     enddo
  elseif (iRescale) then
     do i=ncolumns+1,ncolumns+ncalc
        if (index(label(i),trim(unitslabel(i))).eq.0) label(i) = trim(label(i))//trim(unitslabel(i))
     enddo
  endif
  
  return
end subroutine calc_quantities

!-----------------------------------------------------------------
!
!  utility (private) to add a column (not currently used)
!
!-----------------------------------------------------------------
subroutine addcolumn(inewcolumn,labelin)
 use labels,        only:label
 use settings_data, only:ncolumns,ncalc 
 implicit none
 integer, intent(out) :: inewcolumn
 character(len=*), intent(in) :: labelin

 ncalc = ncalc + 1
 inewcolumn = ncolumns + ncalc
 if (inewcolumn.le.size(label)) then
    label(inewcolumn) = trim(labelin)
 else
    print*,' WARNING!!! too many columns for array dimensions'
    print*,' => change parameter ''maxplot'' in globaldata.f90 and recompile'
 endif

 return
end subroutine addcolumn

!-----------------------------------------------------------------
!
!  utility (private) to fill the array of variable names
!
!-----------------------------------------------------------------
subroutine get_variables(maxlabel,nvars,variables)
 use labels,         only:label
 use settings_units, only:unitslabel
 implicit none
 integer,                        intent(in)  :: maxlabel
 integer,                        intent(out) :: nvars
 character(len=*), dimension(:), intent(out) :: variables
 integer :: i
 !
 !--can use column labels up to the previous quantity calculated
 !
 variables(:) = ' '
 nvars = maxlabel + nextravars

 do i=1,maxlabel
    variables(i) = trim(adjustl(shortlabel(label(i),unitslabel(i))))
 enddo
 do i=1,nextravars
    variables(maxlabel+i) = trim(extravars(i))
 enddo

end subroutine get_variables

!-----------------------------------------------------------------
!
!  utility (private) to strip spaces, escape sequences and
!  units labels from variable names
!
!-----------------------------------------------------------------
elemental function shortlabel(string,unitslab)
 use labels, only:lenlabel
 implicit none
 character(len=lenlabel), intent(in)           :: string
 character(len=*),        intent(in), optional :: unitslab
 character(len=lenlabel)                       :: shortlabel
 integer :: ipos

 shortlabel = string
 !--strip off the units label
 if (present(unitslab)) then
    if (len_trim(unitslab).gt.0) then
    !--remove units label (only do this once)
       ipos = index(trim(shortlabel),trim(unitslab))
       if (ipos.ne.0) then
          shortlabel = shortlabel(1:ipos-1)//&
                       shortlabel(ipos+len_trim(unitslab)+1:len_trim(shortlabel))
       endif
    endif
 endif

 !--remove spaces
 ipos = index(trim(shortlabel),' ')
 do while (ipos.ne.0)
    shortlabel = shortlabel(1:ipos-1)//shortlabel(ipos+1:len_trim(shortlabel))
    ipos = index(trim(shortlabel),' ')
 enddo
 !--remove escape sequences (\d etc.)
 ipos = index(trim(shortlabel),'\')
 do while (ipos.ne.0)
    shortlabel = shortlabel(1:ipos-1)//shortlabel(ipos+2:len_trim(shortlabel))
    ipos = index(trim(shortlabel),'\')
 enddo

end function shortlabel

!-----------------------------------------------------------------
!
!  utility (public) to query whether or not the origin position
!  is actually used in the currently set quantities
!
!-----------------------------------------------------------------
logical function calc_quantities_use_x0()
 implicit none
 integer :: i

 calc_quantities_use_x0 = .false.
 do i=1,maxcalc
    if (index(calcstring(i),trim(extravars(3))).gt.0) calc_quantities_use_x0 = .true.
    if (index(calcstring(i),trim(extravars(4))).gt.0) calc_quantities_use_x0 = .true.
    if (index(calcstring(i),trim(extravars(5))).gt.0) calc_quantities_use_x0 = .true.
 enddo

end function calc_quantities_use_x0

end module calcquantities
