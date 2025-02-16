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
!  Copyright (C) 2005-2018 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
!
! This module allows arbitrary functions to be computed from the
! particle data. Uses the fparser module to parse the functions.
!
!-----------------------------------------------------------------
module calcquantities
 use labels, only:lenlabel,lenunitslabel
 use params, only:maxplot,maxhdr
 implicit none
 public :: calc_quantities,setup_calculated_quantities
 public :: calc_quantities_use_x0,get_calc_data_dependencies
 public :: print_example_quantities

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
 integer, intent(inout) :: ncalc
 integer                :: ncalctot,istart,iend,ipick,ninactive, i
 logical                :: done,first
 character(len=1)       :: charp
 integer, dimension(maxcalc)   :: incolumn
 ipick = ncolumns + 1

 done = .false.
 first = .true.
!
!--on the first call to setup, prefill the list of calculated
!  quantities with ALL of the valid examples.
!
 if (ncalc==0) call print_example_quantities(.true.,ncalc)
 charp = 'a'
 calcmenu: do while (.not.done)
    call check_calculated_quantities(ncalc,ncalctot,incolumn,verbose=.true.)
    ninactive = ncalctot - ncalc

    iend = maxcalc
    if (ncalctot > 0 .or. .not.first) then
       charp='q' !'a'
       if (.not.first) charp = 'q'
       print*
       call prompt(' a)dd to, e)dit, c)lear current list, or q)uit/finish? ',&
                   charp,list=(/'a','e','c','q','s','S','Q'/),noblank=.true.)
       select case(charp)
       case('a')
          istart = ncalctot
          iend   = ncalctot+1
       case('e')
          if (ninactive > 0 .and. ncalc > 0) then
             call prompt(' pick a function to edit ',ipick,-ninactive,-1,ncolumns+1,ncolumns+ncalc)
          elseif (ncalc > 0) then
             call prompt(' pick a function to edit ',ipick,ncolumns+1,ncolumns+ncalc)
          endif
          istart = 0
          do i=1,ncalctot
             if (incolumn(i)==ipick) istart = i - 1
          enddo
          iend = istart + 1
       case('c')
          calcstring(:) = ' '
          calclabel(:) = ' '
          calcunitslabel(:) = ' '
          first = .false.
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

    if (.not.done) call add_calculated_quantities(istart,iend,ncalc,first,incolumn,verbose=.true.)
    first = .false.
 enddo calcmenu

 if (ncalc > 0) then
    if (ncalc < 10) then
       print "(a,i1,a)",' setup ',ncalc,' additional quantities'
    else
       print "(a,i2,a)",' setup ',ncalc,' additional quantities'
    endif
 endif

end subroutine setup_calculated_quantities

!-----------------------------------------------------------------
!
!  utility (private) to add one or more calculated quantities
!  to the current list and/or edit previous settings
!
!-----------------------------------------------------------------
subroutine add_calculated_quantities(istart,iend,ncalc,printhelp,incolumn,verbose)
 use prompting,     only:prompt
 use fparser,       only:checkf
 use settings_data, only:ncolumns,iRescale,required
 use labels,        only:shortstring,headertags,count_non_blank
 integer, intent(in)  :: istart,iend
 integer, intent(out) :: ncalc
 logical, intent(in)  :: printhelp
 integer, dimension(maxcalc), intent(in) :: incolumn
 logical, intent(in)  :: verbose
 integer :: i,j,ntries,ierr,iequal,nvars,nhdr
 logical :: iask
 character(len=120) :: string
 character(len=lenvars), dimension(maxplot+nextravars+maxhdr) :: vars

 i = istart + 1
 ntries = 0
 ncalc = istart

 if (i > maxcalc) then
    print "(/,a,i2,a)",' *** Error, maximum number of calculated quantities (',maxcalc,') reached, cannot add any more.'
    print "(a)",       ' *** If you hit this limit, *please email me* so I can change the default limits!'
    print "(a)",       ' *** (and then edit calc_quantities.f90, changing the parameter "maxcalc" to something higher...)'
    return
 endif

 if (printhelp) then
    print "(/,a)",' Special characters and spaces removed from labels.'
    print "(a)",' Previous quantities can be used in subsequent calculations.'
    print "(/,10(a))",' Valid variables: column labels',(', '''//trim(extravars(j))//'''',j=1,nextravars),&
                ' and header variables:'
    nhdr = count_non_blank(headertags)
    if (nhdr > 0) print "(43(2x,4(a32,1x),/))",headertags(1:nhdr)
 endif
 call print_example_quantities(verbose)

 overfuncs: do while(ntries < 3 .and. i <= iend .and. i <= maxcalc)
    if (len_trim(calcstring(i)) /= 0 .or. ncalc > istart) then
       if (incolumn(i) > 0) then
          write(*,"(a,i2,a)") '[Column ',incolumn(i),']'
       else
          write(*,"(a)") '[Currently inactive]'
       endif
    endif

    if (len_trim(calclabel(i)) > 0) then
       string = trim(calclabel(i))//' = '//trim(calcstring(i))
    else
       string = trim(calcstring(i))
    endif
    call prompt('Enter function string to calculate (blank for none) ',string)
    if (len_trim(string)==0) then
       !
       !--if editing list and get blank string at prompt,
       !  remove the entry and shuffle the list appropriately
       !
       do j=i,maxcalc-1
          !print*,' shuffling ',j,' = '//trim(calcstring(j+1))
          calcstring(j) = trim(calcstring(j+1))
          if (j+1 < size(calclabel)) then
             calclabel(j) = calclabel(j+1)
          endif
          if (len_trim(calcstring(j))==0) then
             if (j==i) then
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
       call splitstring(string,calclabel(i),calcstring(i))
       !
       !--fill variable array with the list of valid variable
       !  names for this column
       !
       call get_variables(ncolumns+i-1,nvars,vars)
       !
       !--check for errors parsing function
       !
       ierr = checkf(shortstring(calcstring(i)),vars(1:nvars))
       if (ierr /= 0 ) then
          ntries = ntries + 1
          print "(a,i1,a)",' error parsing function string: try again (',ntries,' of 3)'
          if (ntries==3) then
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
          required(ncolumns+i) = .true.
       endif
       if (ierr==0) then
          !
          !--prompt for label if not set
          !
          if (iequal==0) then
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
!  utility to split input string into label and function at the
!  equals sign
!
!---------------------------------------------------------------------
subroutine splitstring(string,calclabel,calcstring)
 character(len=*), intent(in)    :: string
 character(len=*), intent(inout) :: calclabel
 character(len=*), intent(out)   :: calcstring
 integer :: iequal

 iequal = index(string,'=')
 if (iequal /= 0) then
    calclabel  = string(1:iequal-1)
    calcstring = string(iequal+1:len_trim(string))
 else
    calcstring = trim(string)
 endif

 !--remove preceding spaces
 calcstring = trim(adjustl(calcstring))

end subroutine splitstring

!---------------------------------------------------------------------
!
!  utility to give a nice list of examples to follow / cut and paste
!  [ this basically replaces what was hardwired into the
!    old calc_quantities routine ]
!
!---------------------------------------------------------------------
subroutine print_example_quantities(verbose,ncalc)
 use labels,        only:label,unitslabel,shortlabel,lenlabel,irho,iutherm,ipr,iBfirst,&
                         ix,icv,idivB,ih,iradenergy,iamvec,labelvec,idustfrac,&
                         ideltav,ivx,headertags
 use settings_data, only:ncolumns,ndim,icoordsnew,ndimV
 use geometry,      only:labelcoord
 use asciiutils,    only:append_number,find_repeated_tags
 use particle_data, only:headervals
 use filenames,     only:ifileopen
 logical :: verbose
 integer, intent(inout), optional :: ncalc
 logical :: prefill
 character(len=lenlabel) :: string,ldfracsum
 integer :: i,j,ivecstart,ierr,ilen,idustfrac1,ndusttypes,nc
 logical :: gotpmag,gotpressure

 gotpmag = .false.
 gotpressure = .false.
 prefill = .false.
 if (present(ncalc)) then
    prefill = .true.
    nc = ncalc
 else
    nc = 0
 endif

 if (verbose) then
    if (prefill) then
       print "(/,a)",' Prefilling list with useful quantities from current data...'
    else
       print "(/,a)",' e.g. '
    endif
 endif

 !--radius
 string = ' '
 if (ndim > 0 .and. icoordsnew==1 .and. ncolumns >= ndim) then
    write(string,"(a)") 'r = sqrt(('// &
          trim(shortlabel(label(ix(1)),unitslabel(ix(1))))//'-'//trim(labelcoord(1,1))//'0)**2'
    ilen = len_trim(string)
    if (ndim > 1) then
       write(string(ilen+1:),"(a,a,a)",iostat=ierr) &
             (' + ('//trim(shortlabel(label(ix(i)),unitslabel(ix(i))))// &
              '-'//trim(labelcoord(i,1))//'0)**2',i=2,ndim),')'
    else
       write(string(ilen+1:),"(a)") ')'
    endif
    call print_or_prefill(prefill,string,nc,ulab=unitslabel(ix(1)))
 elseif (ncolumns >= 2 .and. .not.prefill) then
    !--if ndim=0 give random example to give at least one
    print "(11x,a)",trim(shortlabel(label(1),unitslabel(1)))//'*'&
                    //trim(shortlabel(label(2),unitslabel(2)))
 endif

 !
 !--total dust fraction if multiple dust phases
 !
 ldfracsum = ' '
 if (idustfrac == 0) then
    ! find number of dust species from how many times the label "dustfrac"
    ! is repeated. For gas-only calculations, this will give ndusttypes=0
    call find_repeated_tags('dustfrac',ncolumns,label,idustfrac1,ndusttypes)
    if (ndusttypes > 1) then
       string = 'dustfrac = '//trim(shortlabel(label(idustfrac1),unitslabel(idustfrac1)))
       do i=idustfrac1+1,idustfrac1+ndusttypes-1
          string = trim(string)//'+'//(trim(shortlabel(label(i),unitslabel(i))))
       enddo
       ldfracsum = 'dustfrac'
       call print_or_prefill(prefill,string,nc)
    endif
 elseif (idustfrac > 0) then
    ndusttypes = 1
    ldfracsum = shortlabel(label(idustfrac),unitslabel(idustfrac))
 endif

 !--pressure (requires total dust fraction)
 string = ' '
 if (irho > 0 .and. iutherm > 0) then
    gotpressure = .true.
    if (len_trim(ldfracsum) > 0) then
       write(string,"(a)",iostat=ierr) &
            'pressure = (gamma-1)*(1 - '//trim(ldfracsum)//')*' &
            //trim(shortlabel(label(irho),unitslabel(irho)))//  &
            '*'//trim(shortlabel(label(iutherm),unitslabel(iutherm)))
    else
       write(string,"(a)",iostat=ierr) &
            'pressure = (gamma-1)*'//trim(shortlabel(label(irho),unitslabel(irho)))// &
            '*'//trim(shortlabel(label(iutherm),unitslabel(iutherm)))
    endif
    if (index(unitslabel(irho),'g/cm^3') > 0 .and. &
        index(unitslabel(iutherm),'erg/g') > 0) then
       call print_or_prefill(prefill,string,nc,ulab=' [erg/cm^3]')
    else
       call print_or_prefill(prefill,string,nc)
    endif
 endif
 !
 !--one-fluid dust stuff
 !
 if (len_trim(ldfracsum) > 0 .and. irho > 0) then
    string = ' '
    !--gas density
    write(string,"(a)",iostat=ierr) '\rho_{g} = ' &
                    //trim(shortlabel(label(irho),unitslabel(irho))) &
                    //'*(1 - '//trim(ldfracsum)//')'
    call print_or_prefill(prefill,string,nc,ulab=unitslabel(irho))

    !--(total) dust density
    write(string,"(a)",iostat=ierr) '\rho_{d} = ' &
                    //trim(shortlabel(label(irho),unitslabel(irho))) &
                    //'*'//trim(ldfracsum)
    call print_or_prefill(prefill,string,nc,ulab=unitslabel(irho))

    !--densities of each individual dust phase
    if (ndusttypes > 1) then
       do i = idustfrac1,idustfrac1+ndusttypes-1
          string = '\rho_{d,'
          if (ifileopen > 0) call append_grain_size_label(string,i-idustfrac1+1,headertags,headervals(:,ifileopen),ierr)
          if (ierr /= 0) call append_number(string,i-idustfrac1+1)
          string = trim(string)//'} = ' &
               //trim(shortlabel(label(irho),unitslabel(irho)))//'*'//trim(label(i))
          call print_or_prefill(prefill,string,nc,ulab=unitslabel(irho))
       enddo
    endif
    !--dust-to-gas ratio
    write(string,"(a)",iostat=ierr) 'dust-to-gas ratio = ' &
                    //trim(ldfracsum)//'/(1. - '//trim(ldfracsum)//')'
    call print_or_prefill(prefill,string,nc)

    if (ideltav > 0 .and. ivx > 0 .and. ndimV > 0) then
       if (ndusttypes==1) then
          !--gas velocities
          do i=1,ndimV
             write(string,"(a)",iostat=ierr) trim(labelvec(ivx))//'_{gas,'//trim(labelcoord(i,icoordsnew))//'} = ' &
                             //trim(shortlabel(label(ivx + i-1),unitslabel(ivx + i-1))) &
                             //' - '//trim(shortlabel(label(idustfrac),unitslabel(idustfrac))) &
                             //'*'//trim(shortlabel(label(ideltav + i-1),unitslabel(ideltav + i-1)))
             call print_or_prefill(prefill,string,nc,ulab=unitslabel(ivx))
          enddo
          !--dust velocities
          do i=1,ndimV
             write(string,"(a)",iostat=ierr) trim(labelvec(ivx))//'_{dust,'//trim(labelcoord(i,icoordsnew))//'} = ' &
                             //trim(shortlabel(label(ivx + i-1),unitslabel(ivx + i-1))) &
                             //' + (1 - '//trim(shortlabel(label(idustfrac),unitslabel(idustfrac))) &
                             //')*'//trim(shortlabel(label(ideltav + i-1),unitslabel(ideltav + i-1)))
             call print_or_prefill(prefill,string,nc,ulab=unitslabel(ivx))
          enddo
       else
          ! Still needs to be implemented...
       endif
    endif
 endif

 !
 !--magnitudes of all vector quantities (only if cartesian coords are set)
 !
 ivecstart = 0
 if (icoordsnew==1 .and. ndim > 0 .and. ndimV > 0) then
    do i=1,ncolumns
       if (iamvec(i) > 0 .and. iamvec(i) <= ncolumns .and. iamvec(i) /= ivecstart) then
          ivecstart = iamvec(i)
          string = ' '
          write(string,"(a)",iostat=ierr) '|'//trim(labelvec(ivecstart))//'| '// &
            '= sqrt('//trim(shortlabel(label(ivecstart),unitslabel(ivecstart)))//'**2'
          ilen = len_trim(string)
          if (ndimV > 1) then
             write(string(ilen+1:),"(a,a,a)",iostat=ierr) &
                  (' + '//trim(shortlabel(label(j),unitslabel(j)))//'**2',&
                                 j=ivecstart+1,ivecstart+ndimV-1),')'
          else
             write(string(ilen+1:),"(a)",iostat=ierr) ')'
          endif
          call print_or_prefill(prefill,string,nc,ulab=unitslabel(ivecstart))
       endif
    enddo
 endif
 !--magnetic pressure
 string = ' '
 if (ndim > 0 .and. ndimV > 0 .and. iBfirst > 0 .and. icoordsnew==1) then
    gotpmag = .true.
    write(string,"(a)",iostat=ierr) &
        'P_{mag} = ('//trim(shortlabel(label(iBfirst),unitslabel(iBfirst)))//'**2'
    ilen = len_trim(string)
    if (ndimV > 1) then
       write(string(ilen+1:),"(a,a,a)",iostat=ierr) &
            (' + '//trim(shortlabel(label(i),unitslabel(i))) &
            //'**2',i=iBfirst+1,iBfirst+ndimV-1),')/(2*mu)'
    else
       write(string(ilen+1:),"(a)",iostat=ierr) ')'
    endif
    call print_or_prefill(prefill,string,nc)
 endif
 !--h*div B / B
 if (ndim > 0 .and. ndimV > 0 .and. ih > 0 .and. iBfirst > 0 .and. &
     icoordsnew==1 .and. idivB > 0) then
    write(string,"(a)",iostat=ierr) &
        'h|div B|/|B| = '//trim(shortlabel(label(ih),unitslabel(ih)))//'*abs(' &
                         //trim(shortlabel(label(idivB),unitslabel(idivB)))//')/' &
                         //'sqrt(('//trim(shortlabel(label(iBfirst),unitslabel(iBfirst)))//'^2'
    ilen = len_trim(string)
    if (ndimV > 1) then
       write(string(ilen+1:),"(a,a,a)",iostat=ierr) &
            (' + '//trim(shortlabel(label(i),unitslabel(i)))//'^2',i=iBfirst+1,iBfirst+ndimV-1),'))'
    else
       write(string(ilen+1:),"(a)",iostat=ierr) '))'
    endif
    call print_or_prefill(prefill,string,nc)
 endif
 !--Plasma beta
 if (ndim > 0 .and. ndimV > 0 .and. iBfirst > 0 .and. gotpmag .and. (gotpressure .or. ipr > 0)) then
    write(string,"(a)",iostat=ierr) 'plasma \beta = pressure/P_mag'
    call print_or_prefill(prefill,string,nc,&
         comment='[ assuming pressure and Pmag calculated ]')
 endif

 !--gas temperature if cv present
 if (ndim > 0 .and. iutherm > 0 .and. icv > 0) then
    string = ' '
    write(string,"(a)",iostat=ierr) 'T_{gas} = '//trim(shortlabel(label(iutherm),unitslabel(iutherm)))//'/' &
                    //trim(shortlabel(label(icv),unitslabel(icv)))
    call print_or_prefill(prefill,string,nc,ulab=' [K]')
 endif
 !--radiation temperature
 if (ndim > 0 .and. irho > 0 .and. iradenergy > 0) then
    string = ' '
    write(string,"(a)",iostat=ierr) 'T_{rad} = ('//trim(shortlabel(label(irho),unitslabel(irho)))//'*' &
                    //trim(shortlabel(label(iradenergy),unitslabel(iradenergy)))//'/7.5646e-15)**0.25'
    call print_or_prefill(prefill,string,nc,ulab=' [K]')
 endif

 if (present(ncalc)) ncalc = nc
 if (.not.prefill) print "(a)"

end subroutine print_example_quantities

!-----------------------------------------------------------------
!
!  utility (private) to append enumerated grainsize as a label
!
!-----------------------------------------------------------------
subroutine append_grain_size_label(string,idust,tags,vals,ierr)
 use labels, only:get_label_grain_size,count_non_blank
 character(len=*), intent(inout) :: string
 integer, intent(in)  :: idust
 character(len=*), intent(in) :: tags(:)
 real,    intent(in) :: vals(:)
 integer, intent(out) :: ierr
 integer :: ntags,nd,i

 ntags = count_non_blank(tags)
 nd = 0
 ierr = 1
 do i=1,ntags
    if (index(tags(i),'grainsize') > 0) then
       nd = nd + 1
       if (nd==idust) then
          ierr = 0
          string = trim(string)//get_label_grain_size(vals(i))
       endif
    endif
 enddo

end subroutine append_grain_size_label

!-----------------------------------------------------------------
!
!  utility (private) to get mass of a particular grain species
!  from the header tags
!
!-----------------------------------------------------------------
real function get_mass_of_species(string,idust,tags,vals,ierr)
 use labels, only:count_non_blank
 character(len=*), intent(inout) :: string
 integer, intent(in)  :: idust
 character(len=*), intent(in) :: tags(:)
 real,    intent(in) :: vals(:)
 integer, intent(out) :: ierr
 integer :: ntags,nd,i

 get_mass_of_species = 0.
 ntags = count_non_blank(tags)
 nd = 0
 ierr = 1
 do i=1,ntags
    if (index(tags(i),'mdust_in') > 0) then
       nd = nd + 1
       if (nd==idust) then
          ierr = 0
          get_mass_of_species = vals(i)
       endif
    endif
 enddo

end function get_mass_of_species

!-----------------------------------------------------------------
!
!  utility (private) to either print the example quantity or
!  add it to a prefilled list of calculated quantities
!
!-----------------------------------------------------------------
subroutine print_or_prefill(prefill,string,nc,comment,ulab)
 logical, intent(in) :: prefill
 character(len=*), intent(in) :: string
 integer, intent(inout) :: nc
 character(len=*), intent(in), optional :: comment,ulab
 logical :: already_used
 integer :: i

 if (prefill) then
    if (len_trim(string) > 0) then ! do not prefill blank strings
       nc = nc + 1
       call splitstring(string,calclabel(nc),calcstring(nc))
       if (present(ulab)) calcunitslabel(nc) = trim(ulab)
    endif
 else
    ! do not print strings already in the list
    already_used = .false.
    do i=1,size(calcstring)
       if (trim(calclabel(i))//' = '//trim(calcstring(i))==trim(string)) already_used = .true.
    enddo
    ! append comment if present
    if (present(comment)) then
       if (.not.already_used) print "(6x,a)",trim(string)//'    [ '//trim(comment)//' ]'
    else
       if (.not.already_used) print "(6x,a)",trim(string)
    endif
 endif

end subroutine print_or_prefill

!-----------------------------------------------------------------
!
!  utility (private) to print the current list of calculated
!  quantities, checking that they parse correctly
!
!-----------------------------------------------------------------
subroutine check_calculated_quantities(ncalcok,ncalctot,incolumn,verbose)
 use settings_data,  only:ncolumns,iRescale
 use fparser,        only:checkf
 use labels,         only:label,unitslabel,shortstring,irhodust_start,irhodust_end
 integer, intent(out) :: ncalcok,ncalctot
 integer, dimension(maxcalc), intent(out), optional :: incolumn
 logical, intent(in), optional :: verbose
 integer :: i,ierr,nvars,indexinactive
 character(len=lenvars), dimension(maxplot+nextravars+maxhdr) :: vars
 logical :: isverbose

 if (present(verbose)) then
    isverbose = verbose
 else
    isverbose = .true.
 endif

 ncalcok = 0
 ncalctot = 0
 indexinactive = 0
 i = 1
 if (present(incolumn)) incolumn(:) = 0

 if (all(len_trim(calcstring(:))==0)) return
 if (isverbose) print "(/,a)", ' Current list of calculated quantities:'
 do while(i <= maxcalc .and. len_trim(calcstring(i)) /= 0)
    !
    !--get the list of valid variable names for this column
    !
    call get_variables(ncolumns+ncalcok,nvars,vars)
    !
    !--check that the function parses
    !
    ierr = checkf(shortstring(calcstring(i)),vars(1:nvars),Verbose=.false.)
    if (ierr==0) then
       ncalcok = ncalcok + 1
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
       if (isverbose) then
          if (iRescale) then
             print "(1x,i2,') ',a50,' [OK]',1x,a)",ncolumns+ncalcok,trim(calclabel(i))//' = '//calcstring(i),trim(calcunitslabel(i))
          else
             print "(1x,i2,') ',a50,' [OK]')",ncolumns+ncalcok,trim(calclabel(i))//' = '//calcstring(i)
          endif
       endif
       !
       !--recognise the dust density in the list of calculated quantities
       !
       if (trim(calcstring(i))=='density*dustfrac1') irhodust_start = ncolumns+ncalcok
       if (calcstring(i)(1:16)=='density*dustfrac')    irhodust_end = ncolumns+ncalcok ! overwrite until last density*dustfrac
    else
       indexinactive = indexinactive - 1
       if (isverbose) then
          print "(i3,') ',a50,' [INACTIVE]')",indexinactive,trim(calclabel(i))//' = '//calcstring(i)
       endif
       if (present(incolumn)) incolumn(i) = indexinactive
    endif
    ncalctot = i
    i = i + 1
 enddo
 if (ncalcok==0 .and. isverbose) print "(a)",' (none)'

end subroutine check_calculated_quantities

!-----------------------------------------------------------------
!
!  utility (private) to check dependencies of calculated
!  quantities, so that we only read what is necessary
!  from the dump file
!
!-----------------------------------------------------------------
subroutine get_calc_data_dependencies(required)
 use params,         only:maxplot
 use settings_data,  only:debugmode
 use fparser,        only:checkf
 use labels,         only:label,shortlabel,shortstring
 logical, dimension(0:maxplot), intent(inout) :: required
 character(len=lenvars), dimension(maxplot+nextravars+maxhdr) :: vars
 integer, dimension(maxcalc) :: incolumn
 integer :: ncalcok,ncalctot,nvars,i,j

 call check_calculated_quantities(ncalcok,ncalctot,incolumn,verbose=.false.)

 do i=ncalctot,1,-1   ! go in REVERSE order to get recursive dependencies properly
    if (incolumn(i) > 0) then
       if (required(incolumn(i))) then
          if (debugmode) then
             print*,'DEBUG: computing dependencies for '//trim(label(incolumn(i)))//&
                    ' = '//trim(shortstring(calcstring(i)))
          endif
          !
          !--get the list of valid variable names for this column
          !
          call get_variables(incolumn(i),nvars,vars)
          !
          !--check if the string contains any preceding variables
          !
          do j=1,incolumn(i)-1
             !
             !--this could be smarter here (at the moment we just check for
             !  matching substrings, but we should check for the use of the
             !  string as an actual variable -- this is mainly an issue for
             !  single letter variables like x,y,z etc)
             !
             if (index(shortlabel(calcstring(i)),trim(vars(j))) /= 0) then
                if (debugmode) print*,'DEBUG: -> depends on '//trim(label(j))
                required(j) = .true.
             endif
          enddo
       endif
    endif
 enddo

end subroutine get_calc_data_dependencies

!-----------------------------------------------------------------
!
!  actually compute the extra quantities from the particle data
!
!-----------------------------------------------------------------
subroutine calc_quantities(ifromstep,itostep,dontcalculate)
 use labels,         only:label,unitslabel,labelvec,iamvec,ix,ivx,irho,shortstring, &
                           count_non_blank,headertags
 use particle_data,  only:dat,npartoftype,gamma,time,headervals,maxpart,maxstep,maxcol,iamtype
 use settings_data,  only:ncolumns,ncalc,iRescale,xorigin,debugmode,ndim,required,iverbose, &
                           icoords,icoordsnew,ipartialread,track_string
 use mem_allocation, only:alloc
 use settings_units, only:units,units_calc
 use fparser,        only:checkf,parsef,evalf,EvalerrMsg,EvalErrType,rn,initf,endf
 use params,         only:maxplot
 use timing,         only:wall_time,print_time
 use geomutils,      only:change_coords
 use part_utils,     only:get_tracked_particle
 integer, intent(in) :: ifromstep, itostep
 logical, intent(in), optional :: dontcalculate
 integer :: i,j,ncolsnew,ierr,icalc,ntoti,nvars,ncalctot,nused,itrackpart
 integer :: ndust,nhdr
 logical :: skip
!  real, parameter :: mhonkb = 1.6733e-24/1.38e-16
!  real, parameter :: radconst = 7.5646e-15
!  real, parameter :: lightspeed = 3.e10   ! in cm/s (cgs)
 real(kind=rn), dimension(maxplot+nextravars+maxhdr)          :: vals,unitvals
 character(len=lenvars), dimension(maxplot+nextravars+maxhdr) :: vars
 real, dimension(3) :: x0,v0
 real :: t1,t2

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
 call check_calculated_quantities(ncalc,ncalctot,verbose=(.not.skip .and. iverbose > 0))

 if (.not.skip .and. ncalc > 0) then
    nused = 0
    if (.not.ipartialread) then
       !
       !--need to be careful if data file has been read fully
       !  as in this case we also assume all calculated quantities
       !  have been done. So need to make sure that all quantities
       !  *are* actually calculated in this case.
       !
       required(:) = .true.
    endif
    do i=1,ncalc
       if (required(ncolumns+i)) nused = nused + 1
    enddo
    if (iverbose > 0) print "(2(a,i2),a,/)",' Calculating ',nused,' of ',ncalctot,' additional quantities...'
 endif
 ncolsnew = ncolumns + ncalc
 if (ncolsnew > maxcol) call alloc(maxpart,maxstep,ncolsnew)

 !
 !--reset iamvec to zero for calculated columns
 !
 iamvec(ncolumns+1:ncolsnew) = 0
 labelvec(ncolumns+1:ncolsnew) = ' '

 !
 !--evaluate functions in turn
 !
 if (.not.skip .and. ncalc > 0) then
    call initf(ncalc)
    !
    !--compile each function into bytecode
    !
    icalc = 1
    do i=1,maxcalc
       if (icalc <= ncalc) then

          !
          !--get the list of valid variable names for this column
          !
          call get_variables(ncolumns+icalc-1,nvars,vars)
          !
          !--now actually parse the function
          !
          call parsef(icalc,shortstring(calcstring(i)),vars(1:nvars),err=ierr,Verbose=.false.)
          if (ierr==0) then
             icalc = icalc + 1
          endif
       endif
    enddo
    !
    !--evaluate functions from particle data
    !
    call wall_time(t1)
    do i=ifromstep,itostep
       ntoti = SUM(npartoftype(:,i))
       ndust = npartoftype(2,i)
       !
       !--set origin position
       !
       v0(:) = 0.
       itrackpart = get_tracked_particle(track_string,npartoftype(:,i),iamtype(:,i),dat(:,:,i),irho)
       if (itrackpart > 0 .and. itrackpart <= ntoti) then
          x0(:) = 0.
          if (ix(1) > 0 .and. ix(1) <= ncolumns) then
             x0(1) = dat(itrackpart,ix(1),i)
          else
             print*,'** internal error: tracking particle set but cannot locate x coordinate'
          endif
          if (ix(2) > 0 .and. ix(2) <= ncolumns) x0(2) = dat(itrackpart,ix(2),i)
          if (ix(3) > 0 .and. ix(3) <= ncolumns) x0(3) = dat(itrackpart,ix(3),i)
          if (i==ifromstep) then
             print "(a,i10)",' using position of tracked particle ',itrackpart
             print "(a,3(e11.3),/)",' (x0,y0,z0) = ',dat(itrackpart,ix(1:ndim),i)
          endif
          if (ivx > 0 .and. ivx+ndim-1 <= ncolumns) then
             v0(1:ndim) = dat(itrackpart,ivx:ivx+ndim-1,i)
          endif
       else
          x0(:) = xorigin(:)
       endif

       do icalc=1,ncalc
          if (required(ncolumns+icalc)) then
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
             if (iRescale) then
                unitvals(ncolumns+icalc)   = units(0)
                unitvals(ncolumns+icalc+1) = 1.
                unitvals(ncolumns+icalc+2) = units(ix(1))
                unitvals(ncolumns+icalc+3) = units(ix(2))
                unitvals(ncolumns+icalc+4) = units(ix(3))
             endif
             nhdr = count_non_blank(headertags)
             do j=1,nhdr
                vals(ncolumns+icalc+4+j) = headervals(j,i)
                if (iRescale) unitvals(ncolumns+icalc+4+j) = 1.
             enddo
             if (icoordsnew /= icoords .and. ndim > 0 .and. all(ix(1:ndim) > 0)) then
                !
                !--if alternative coordinate system is in use, then we need to apply
                !  the coordinate transformations to the data BEFORE using it
                !  to calculate additional quantities
                !
                do j=1,ntoti
                   vals(1:ncolumns+icalc-1) = dat(j,1:ncolumns+icalc-1,i)
                   call change_coords(vals(1:ncolumns+icalc-1),ncolumns+icalc-1,&
                                       ndim,icoords,icoordsnew,x0(1:ndim),v0(1:ndim))
                   !--evaluate function with transformed values
                   dat(j,ncolumns+icalc,i) = real(evalf(icalc,vals(1:ncolumns+icalc+nextravars+nhdr-1)))
                enddo
                ! evaluate units of the new function
                if (iRescale) then
                   unitvals(1:ncolumns+icalc-1) = units(1:ncolumns+icalc-1)
                   call change_coords(unitvals(1:ncolumns+icalc-1),ncolumns+icalc-1,&
                                      ndim,icoords,icoordsnew,x0(1:ndim),v0(1:ndim))
                   ! DP: see comment below
                   units_calc(ncolumns+icalc) = units(ncolumns+icalc) !real(evalf(icalc,unitvals(1:ncolumns+icalc+nextravars+nhdr-1)))
                endif
             else
                !!$omp parallel do default(none) private(j,vals,icolumn) shared(dat,i,icalc,ncolumns)
                do j=1,ntoti
                   vals(1:ncolumns+icalc-1) = dat(j,1:ncolumns+icalc-1,i)
                   dat(j,ncolumns+icalc,i) = real(evalf(icalc,vals(1:ncolumns+icalc+nextravars+nhdr-1)))
                enddo
                if (iRescale) then
                   unitvals(1:ncolumns+icalc-1) = units(1:ncolumns+icalc-1)
                   !
                   ! DP: ideally we would fix the line below to compute the exact solution unit scaling directly
                   ! for an arbitrary function string, but we have to take out the additions and
                   ! subtractions, so would have to parse a reduced function, not the original function
                   !
                   ! For now, this is done manually by recognising the dimensionality of certain computed quantities
                   ! in identify_calculated_quantities
                   !
                   units_calc(ncolumns+icalc) = units(ncolumns+icalc) !real(evalf(icalc,unitvals(1:ncolumns+icalc+nextravars+nhdr-1)))
                endif
               !!$omp end parallel do
             endif
             if (EvalErrType /= 0) then
                print "(a)",' ERRORS evaluating '//trim(calcstring(icalc))//': ' &
                             //trim(EvalerrMsg())
             endif
             !
             !--identify calculated quantities based on the label
             !
             if (i==ifromstep) then
                call identify_calculated_quantity(label(ncolumns+icalc),ncolumns,ncolumns+icalc)
             endif
          else
             if (debugmode) print*,'DEBUG: ',icalc,' skipping '//trim(label(ncolumns+icalc))//' (not required)'
          endif
       enddo
    enddo
    call endf
    call wall_time(t2)
    if (t2-t1 > 1.) call print_time(t2-t1)
 endif
 !
 !--override units of calculated quantities if necessary
 !
 if (iRescale .and. any(abs(units(ncolumns+1:ncolumns+ncalc)-1.0) > tiny(0.)) &
      .and. .not.skip) then
    !write(*,"(/a)") ' rescaling data...'
    do i=ncolumns+1,ncolumns+ncalc
       if (abs(units(i)-1.0) > tiny(0.) .and. abs(units(i)) > tiny(0.)) then
          dat(:,i,ifromstep:itostep) = dat(:,i,ifromstep:itostep)*units(i)
       endif
       if (index(label(i),trim(unitslabel(i)))==0) label(i) = trim(label(i))//trim(unitslabel(i))
    enddo
 elseif (iRescale) then
    do i=ncolumns+1,ncolumns+ncalc
       if (index(label(i),trim(unitslabel(i)))==0) label(i) = trim(label(i))//trim(unitslabel(i))
       if (debugmode) print*,'DEBUG: column ',i,trim(label(i)),' units_calc = ',units_calc(i)
    enddo
 endif

 return
end subroutine calc_quantities

!-----------------------------------------------------------------
!
!  utility (private) to internally identify a calculated quantity
!  so that the relevant exact solutions can be plotted for that
!  quantity. This is mainly just the radius, but can include
!  other things also.
!
!-----------------------------------------------------------------
subroutine identify_calculated_quantity(labelcol,ncolumns,icolumn)
 use asciiutils,    only:lcase
 use labels,        only:irad,ike,ipr,ikappa,itemp,label_synonym,ix,ivx
 use settings_data, only:debugmode,iRescale
 use settings_units,only:units_calc,units
 character(len=*), intent(in) :: labelcol
 integer, intent(in) :: ncolumns,icolumn
 !
 !--identify quantities based on the label
 !  only do this if the location flags are not already set
 !  (e.g. in the data read) - but DO overwrite if they
 !  are calculated quantities as the locations can change
 !
 select case(label_synonym(labelcol))
 case('r','radius','rad')
    call assign_column(irad,icolumn,ncolumns,debugmode,'radius')
    if (iRescale .and. ix(1) > 0) units_calc(icolumn) = units(ix(1))
 case('kinetic energy','ke','1/2 v^2','v^2/2')
    call assign_column(ike,icolumn,ncolumns,debugmode,'kinetic energy')
    if (iRescale .and. ivx > 0) units_calc(icolumn) = units(ivx)**2
 case('pressure','pr','p')
    call assign_column(ipr,icolumn,ncolumns,debugmode,'pressure')
 case('kappa','opacity')
    call assign_column(ikappa,icolumn,ncolumns,debugmode,'opacity')
 case('temperature','temp')
    call assign_column(itemp,icolumn,ncolumns,debugmode,'temperature')
 end select

end subroutine identify_calculated_quantity

!-----------------------------------------------------------------
!
!  helper routine for above
!
!-----------------------------------------------------------------
subroutine assign_column(i,icolumn,ncolumns,debugmode,string)
 integer, intent(inout) :: i
 integer, intent(in)    :: icolumn,ncolumns
 logical, intent(in)    :: debugmode
 character(len=*), intent(in) :: string

 if (i <= 0 .or. i > ncolumns) then
    i = icolumn
    if (debugmode) print "(1x,a,i2,a)",'identifying column ',icolumn,' as the '//trim(string)
 endif

end subroutine assign_column

!-----------------------------------------------------------------
!
!  utility (private) to fill the array of variable names
!
!-----------------------------------------------------------------
subroutine get_variables(maxlabel,nvars,variables)
 use labels,         only:label,shortlabel,unitslabel,headertags,count_non_blank
 integer,                        intent(in)  :: maxlabel
 integer,                        intent(out) :: nvars
 character(len=*), dimension(:), intent(out) :: variables
 integer :: i,nheader
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
 nheader = count_non_blank(headertags)
 nvars = nvars + nheader
 do i=1,nheader
    variables(maxlabel+nextravars+i) = trim(headertags(i))
 enddo

end subroutine get_variables

!-----------------------------------------------------------------
!
!  utility (public) to query whether or not the origin position
!  is actually used in the currently set quantities
!
!-----------------------------------------------------------------
logical function calc_quantities_use_x0()
 integer :: i

 calc_quantities_use_x0 = .false.
 do i=1,maxcalc
    if (index(calcstring(i),trim(extravars(3))) > 0) calc_quantities_use_x0 = .true.
    if (index(calcstring(i),trim(extravars(4))) > 0) calc_quantities_use_x0 = .true.
    if (index(calcstring(i),trim(extravars(5))) > 0) calc_quantities_use_x0 = .true.
 enddo

end function calc_quantities_use_x0

end module calcquantities
