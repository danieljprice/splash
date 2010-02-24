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
!  Copyright (C) 2005-2010 Daniel Price. All rights reserved.
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
 implicit none
 public :: calc_quantities,setup_calculated_quantities

 integer, parameter, private :: maxcalc = 10
 character(len=60),            dimension(maxcalc) :: calcstring = ' '
 character(len=lenlabel),      dimension(maxcalc) :: calclabel = ' '
 character(len=lenunitslabel), dimension(maxcalc) :: calcunitslabel = ' '

 integer, parameter, private :: nextravars = 5
 character(len=lenlabel), dimension(nextravars), parameter, private :: &
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
 use prompting,     only:prompt
 implicit none
 integer, intent(out) :: ncalc
 integer              :: i,istart,iend
 logical              :: done,first
 character(len=1)     :: charp
 
 done = .false.
 first = .true.
 charp = 'a'
 calcmenu: do while (.not.done)
    call check_calculated_quantities(ncalc,i)

    iend = maxcalc
    if (i.gt.0 .or. .not.first) then
       charp='a'
       print*
       call prompt(' a)dd to, e)dit, c)lear current list or q)uit/finish? ',&
                   charp,list=(/'a','e','c','q','s','S','Q'/),noblank=.true.)
       select case(charp)
       case('a')
          istart = i
          iend   = i+1
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

    if (.not.done) call add_calculated_quantities(istart,iend,ncalc,first)
    first = .false.
 enddo calcmenu
 print*,' setup ',ncalc,' additional quantities'
 
end subroutine setup_calculated_quantities

!-----------------------------------------------------------------
!
!  utility (private) to add one or more calculated quantities
!  to the current list and/or edit previous settings
!
!-----------------------------------------------------------------
subroutine add_calculated_quantities(istart,iend,ncalc,printhelp)
 use prompting,     only:prompt
 use fparser,       only:checkf
 use labels,        only:label,lenlabel
 use settings_data, only:ncolumns,iRescale
 use settings_units,only:unitslabel
 implicit none
 integer, intent(in)  :: istart,iend
 integer, intent(out) :: ncalc
 logical, intent(in)  :: printhelp
 integer :: i,j,ntries,ierr,iequal
 logical :: iask

 if (printhelp) then
    print "(/,a)",' Specify a function to calculate from the data '
    print "(10(a))",' Valid variables are the column labels',(', '''//trim(extravars(i))//'''',i=1,nextravars-1),&
                ' and '''//trim(extravars(nextravars))//''' (origin setting) '
    print "(a)",' Spaces, escape sequences (\d) and units labels are removed from variable names'
    print "(a)",' Note that previously calculated quantities can be used in subsequent calculations'
 endif
 call print_example_quantities()
 
 i = istart + 1
 ntries = 0
 ncalc = istart
 overfuncs: do while(ntries.lt.3 .and. i.le.iend)
    if (len_trim(calcstring(i)).ne.0 .or. ncalc.gt.istart) write(*,"(a,i2,a)") '[Column ',ncolumns+i,']'
    call prompt('Enter function string to calculate (blank for none) ',calcstring(i))
    if (len_trim(calcstring(i)).eq.0) then
       !
       !--if editing list and get blank string at prompt,
       !  remove the entry and shuffle the list appropriately
       !
       do j=i,maxcalc-1
          !print*,' shuffling ',j,' = '//trim(calcstring(j+1))
          calcstring(j) = trim(calcstring(j+1))
          if (ncolumns+j+1.lt.size(label)) then
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
       !
       iequal = index(calcstring(i),'=')
       if (iequal.ne.0) then
          calclabel(i) = calcstring(i)(1:iequal-1)
          calcstring(i) = calcstring(i)(iequal+1:len_trim(calcstring(i)))
       endif
       !
       !--check for errors parsing function
       !
       ierr = checkf(shortlabel(calcstring(i)),&
              (/(shortlabel(label(i),unitslabel(i)),i=1,ncolumns+i-1),extravars/))
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
 use labels,        only:label,lenlabel,irho,iutherm,ivx,ix,icv,iradenergy
 use settings_data, only:ncolumns,ndim,icoordsnew,ndimV
 use settings_units,only:unitslabel
 use geometry,      only:labelcoord
 implicit none
 integer :: i

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
    print "(a)",'           pressure = (gamma-1)*'//trim(shortlabel(label(irho),unitslabel(irho)))// &
                '*'//trim(shortlabel(label(iutherm),unitslabel(iutherm)))
 endif
 !--magnitude of v
 if (ndim.gt.0 .and. ndimV.gt.0 .and. ivx.gt.0 .and. icoordsnew.eq.1) then
    write(*,"(11x,a)",ADVANCE='NO') '|v| = sqrt('//trim(shortlabel(label(ivx),unitslabel(ivx)))//'**2'
    if (ndimV.gt.1) then
       write(*,"(a,a,a)") (' + '//trim(shortlabel(label(i),unitslabel(i)))//'**2',i=ivx+1,ivx+ndimV-1),')'
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
subroutine check_calculated_quantities(ncalcok,ncalctot)
 use labels,         only:label
 use settings_data,  only:ncolumns
 use settings_units, only:unitslabel
 use fparser,        only:checkf
 implicit none
 integer, intent(out) :: ncalcok,ncalctot
 integer :: i,ierr

 ncalcok = 0
 ncalctot = 0
 i = 1
 print "(/,a)", ' Current list of calculated quantities:'
 do while(i.le.maxcalc .and. len_trim(calcstring(i)).ne.0)
    ierr = checkf(shortlabel(calcstring(i)),&
            (/(shortlabel(label(i),unitslabel(i)),i=1,ncolumns+ncalcok),extravars/),Verbose=.false.)

    if (ierr.eq.0) then
       ncalcok = ncalcok + 1
       print "(1x,i2,') ',a50,' [OK]')",ncolumns+ncalcok,trim(calclabel(i))//' = '//calcstring(i)
    else
       print "(1x,'XX) ',a50,' [INACTIVE]')",trim(calclabel(i))//' = '//calcstring(i)
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
  use labels,         only:label,labelvec,iamvec
  use particle_data,  only:dat,npartoftype,gamma,time,maxpart,maxstep,maxcol
  use settings_data,  only:ncolumns,ncalc,iRescale,xorigin,debugmode !,itrackpart
  use mem_allocation, only:alloc
  use settings_units, only:unitslabel,units
  use fparser,        only:checkf,parsef,evalf,EvalerrMsg,EvalErrType,rn,initf,endf
  use params,         only:maxplot
  implicit none
  integer, intent(in) :: ifromstep, itostep
  logical, intent(in), optional :: dontcalculate
  integer :: i,j,ncolsnew,ierr,icalc,ntoti
  logical :: skip
!  real, parameter :: mhonkb = 1.6733e-24/1.38e-16
!  real, parameter :: radconst = 7.5646e-15
!  real, parameter :: lightspeed = 3.e10   ! in cm/s (cgs)
  real(kind=rn), dimension(maxplot+nextravars) :: vals
  
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
  i = 1
  print*
  do while (i.le.maxcalc .and. len_trim(calcstring(i)).gt.0 .and.ncolumns+ncalc.lt.size(label))
     ierr = checkf(shortlabel(calcstring(i)),&
            (/(shortlabel(label(i),unitslabel(i)),i=1,ncolumns+ncalc),extravars/))
     if (ierr.eq.0 ) then
        print "(a,a10,' = ',a)",' calculating ',trim(calclabel(i)),trim(calcstring(i))
        ncalc = ncalc + 1
        !
        !--now actually assign the calculated quantity to a particular column
        !  and set required information for the column
        !
        label(ncolumns+ncalc) = trim(calclabel(i))
        if (iRescale) label(ncolumns+icalc) = trim(label(ncolumns+ncalc))//trim(calcunitslabel(i))
     else
        print "(a)",' error parsing function '//trim(calclabel(i))//' = '//trim(calcstring(i))
     endif
     i = i + 1
  enddo
  
  if (.not.skip) print*,'calculating ',ncalc,' additional quantities...'
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
           call parsef(icalc,shortlabel(calcstring(i)), &
                (/(shortlabel(label(i),unitslabel(i)),i=1,ncolumns+icalc-1),extravars/),err=ierr,Verbose=.false.)
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
        do icalc=1,ncalc
           if (debugmode) print*,'DEBUG: ',icalc,' calculating '//trim(label(ncolumns+icalc))
           !
           !--additional settings allowed in function evaluations
           !  i.e., time and gamma from dump file and current origin settings
           !  make sure the number here aligns with the "nextravars" setting
           !
           vals(ncolumns+icalc)   = time(i)
           vals(ncolumns+icalc+1) = gamma(i)
           vals(ncolumns+icalc+2) = xorigin(1)
           vals(ncolumns+icalc+3) = xorigin(2)
           vals(ncolumns+icalc+4) = xorigin(3)
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
!  utility (private) to strip spaces, escape sequences and
!  units labels from variable names
!
!-----------------------------------------------------------------
elemental function shortlabel(string,unitslab)
 use labels, only:lenlabel
 implicit none
 character(len=lenlabel), intent(in)  :: string
 character(len=lenlabel) :: shortlabel
 character(len=*), intent(in), optional :: unitslab
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

end module calcquantities
