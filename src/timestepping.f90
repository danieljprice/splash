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

module timestepping
 implicit none
 public :: timestep_loop

 private

contains
!
! This subroutine drives the main plotting loop
!
subroutine timestep_loop(ipicky,ipickx,irender,icontourplot,ivecplot)
  use filenames,         only:nsteps,ifileopen
  use particle_data,     only:iamtype,npartoftype,masstype,time,gamma,dat,time_was_read
  use settings_data,     only:istartatstep,iendatstep,nfreq,DataIsBuffered, &
                              iUsesteplist,isteplist,ncolumns,ipartialread
  use settings_page,     only:interactive,nstepsperpage,iColourEachStep,iChangeStyles,nomenu
  use timestep_plotting, only:initialise_plotting,plotstep
  use plotlib,           only:plot_close
  implicit none
  integer, intent(in) :: ipicky,ipickx,irender,icontourplot,ivecplot
  integer :: ipos, istep, ilocindat, iadvance, istepsonpage, istepprev
  logical :: ipagechange

  call initialise_plotting(ipicky,ipickx,irender,icontourplot,ivecplot)

  !----------------------------------------------------------------------------
  ! loop over timesteps (flexible to allow going forwards/backwards in
  !                      interactive mode)
  !
  ! bookkeeping is as follows:
  !         ipos : current step number or position in steplist array
  !    ilocindat : current or requested location in the dat array
  !                (usually 1, but can be > 1 if more than one step per file
  !                 or data is buffered to memory)
  !                (requested location is sent to get_nextstep which skips files
  !                or steps appropriately and sets actual location)
  !       istep :  current step number
  !                (as if all steps were in memory in sequential order)
  ! istepsonpage : number of steps which have been plotted on current page
  !                this is used because steps are coloured/marked differently
  !                from this routine (call to colour_timestep)
  !
  !----------------------------------------------------------------------------
  ipos = istartatstep
  iadvance = nfreq   ! amount to increment timestep by (changed in interactive)
  istepsonpage = 0
  istep = istartatstep
  istepprev = 0
  !--if the current file has only been partially read,
  !  make sure we read the file again now that we may have different plotting options
  if (ipartialread) ifileopen = 0

  over_timesteps: do while (ipos.le.iendatstep)

     ipos = max(ipos,1) !--can't go further back than ipos=1

     if (iUseStepList) then
        istep = isteplist(ipos)
        if (istep.gt.nsteps) then
           print*,'ERROR: step > nsteps in step list, setting step = last'
           istep = nsteps
        elseif (istep.le.0) then
        !--this should never happen
           stop 'internal error: corrupted step list: please send bug report'
        endif
     else
        istep = ipos
        if (istep.ge.nsteps) then
           istep = nsteps
           iendatstep = istep
           ipos = min(ipos,iendatstep)
        endif
     endif

     !
     !--make sure we have data for this timestep
     !
     if (DataIsBuffered) then
        !--if data is in memory, we just go to the position in dat
        if (istep.gt.nsteps) then
           print*,'error: step # > nsteps, setting step = last'
           istep = nsteps
        endif
        ilocindat = istep
     else
        !--otherwise read file containing this step into memory and get position in dat array
        !  (note that nsteps can change in get_nextstep, so may need to re-evaluate
        !   whether we are on the last step or not, and adjust iendatstep to last step)
        !
        call get_nextstep(istep,ilocindat)

        if (.not.iUseStepList .and. istep.ge.nsteps) then
           !--reset step position to last useable timestep (ie. nsteps)
           istep = nsteps
           iendatstep = istep
           !--use interactive halt at last step (ie. set position = last position)
           if (interactive) then
              ipos = min(ipos,iendatstep)
              !--if istep has changed, may need to re-read step
              !  (get_nextstep does nothing if istep is the same)
              call get_nextstep(istep,ilocindat)
           endif
        endif
        !--this is a general "catch all" when step cannot be located
        if (ilocindat.le.0) then
           print*,'ERROR: could not locate timestep'
           exit over_timesteps
        endif
     endif

     !
     !--write timestepping log
     !
     if (.not.time_was_read(time(ilocindat))) then
        print 32, istep
     elseif (time(ilocindat).lt.1.e-2 .or. time(ilocindat).gt.1.e2) then
        print 33, time(ilocindat),istep
     else
        print 34, time(ilocindat),istep
     endif
32   format (5('-'),' t = (not read), dump #',i5,1x,18('-'))
33   format (5('-'),' t = ',es9.2,', dump #',i5,1x,18('-'))
34   format (5('-'),' t = ',f8.2,', dump #',i5,1x,18('-'))

     istepsonpage = istepsonpage + 1
!     if ((nstepsperpage.gt.1 .and. istepsonpage > 1 .and. istepsonpage.le.nstepsperpage).or.nstepsperpage.eq.0) then
     if ((nstepsperpage.gt.1 .and. istepsonpage.le.nstepsperpage).or.nstepsperpage.eq.0) then
        ipagechange = .false.
     else
        istepsonpage = 1
        ipagechange = .true.
     endif
     !--colour the timestep if appropriate
     if ((nstepsperpage.eq.0 .or. nstepsperpage.gt.1) .and. (iColourEachStep .or. iChangeStyles)) then
        call colour_timestep(istepsonpage,iColourEachStep,iChangeStyles)
     else
     !--otherwise set default colours for each particle type
     !  (do not call if repeating same step so interactive colours stick for same step)
        if (istep.ne.istepprev) call colourparts_default(npartoftype(:,ilocindat),iamtype(:,ilocindat))
        istepprev = istep
     endif

!     print*,'ipos = ',ipos,' istep = ',istep,' iposindat = ',ilocindat
     call plotstep(ipos,istep,istepsonpage,irender,icontourplot,ivecplot,iamtype(:,ilocindat), &
                   npartoftype(:,ilocindat),masstype(:,ilocindat),dat(:,:,ilocindat), &
                   time(ilocindat),gamma(ilocindat),ipagechange,iadvance)
!
!--increment timestep -- iadvance can be changed interactively
!
     if (iadvance.eq.-666) exit over_timesteps ! this is the interactive quit signal
     ipos = ipos + iadvance ! if ipos goes over iendatstep, this ends the loop

  enddo over_timesteps

  if (.not.interactive) then
     if (nomenu) then
        !--gracefully exit
        print "(/,a,/)",'Finished plotting: Many thankyous for your kind custom.'
     else ! prepare to return to main menu
        print*,'press return to finish'
        read*
        !--if somehow the data has become corrupted (e.g. last file full of rubbish)
        !  read in the first dump again
        if (ncolumns.le.0) then
           print*,'data is corrupted: re-reading first data file'
           call get_nextstep(1,ilocindat)
        endif
     endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call plot_close

  return

end subroutine timestep_loop

!-------------------------------------------------------------------
! works out whether or not we need to read another dump into memory
!-------------------------------------------------------------------
recursive subroutine get_nextstep(istep,ilocindat)
 use filenames,     only:nstepsinfile,nfiles,ifileopen,nsteps
 use getdata,       only:get_data
 use settings_data, only:buffer_steps_in_file
 implicit none
 integer, intent(in)  :: istep
 integer, intent(out) :: ilocindat
 integer :: ifile,nstepstotal,nstepsprev,iposinfile

 !
 !--request is for step istep
 !  need to determine which file step i is in,
 !  whether or not it is already in memory (and read it if not)
 !  and finally determine position of requested step in dat array
 !
 ifile = 0
 nstepstotal = 0
 nstepsprev = 0
 iposinfile = 1 ! this should always be overwritten anyway
 do while (nstepstotal.lt.istep .and. ifile.lt.nfiles)
    ifile = ifile + 1
    nstepsprev = nstepstotal
    nstepstotal = nstepstotal + nstepsinfile(ifile)
 enddo

 if (nstepstotal.ge.istep) then
 !--set position in dat array depending on how many steps are in the file
    iposinfile = istep - nstepsprev
 else
 !--this is where we cannot locate the timestep in the data (not enough steps)
 !  ie. ifile > nfiles
    print*,'reached last useable timestep'
    ilocindat = 0
    return
 endif

! print*,'step ',istep,' in file ',ifile,' nsteps = ',nsteps
! print*,'position in file = ',iposinfile
 if (istep.gt.nsteps) then
    ilocindat = 0
    return
 endif
 ilocindat = iposinfile

 !
 !--if data is not stored in memory, read next step from file
 !  At the moment assumes number of steps in each file are the same
 !

 !--neither or these two error conditions should occur
 if (ifile.gt.nfiles) then
    print*,'*** get_nextstep: error: ifile > nfiles'
 elseif (ifile.lt.1) then
    print*,'*** get_nextstep: error: request for file < 1'
 elseif (ifile.ne.ifileopen) then
 !
 !--read next data file and determine position in file
 !
    call get_data(ifile,.true.)
 !
 !--because nstepsinfile is predicted for files which have
 !  not been opened, we may have the situation where
 !  iposinfile does not point to a real timestep (ie. iposinfile > nstepsinfile).
 !  In this case we query the step again with our better knowledge of nstepsinfile.
 !
    if (iposinfile.gt.nstepsinfile(ifile)) then
       print*,'not enough steps in file... trying next file'
       call get_nextstep(istep,ilocindat)
    endif
 elseif (.not.buffer_steps_in_file) then
    if (iposinfile <= nstepsinfile(ifile)) then
       ! if all steps have not been read from the file, read just the next one
       call get_data(ifile,.true.,iposinfile=iposinfile)
       ilocindat = 1
    endif
 endif

 return
end subroutine get_nextstep

!-------------------------------------------------------------
! colours all the particles a given colour for this timestep
! and/or changes the marker type for type 1 particles
!-------------------------------------------------------------
subroutine colour_timestep(istep,iChangeColours,iChangeStyles)
  use particle_data, only:icolourme
  use settings_part, only:linecolourthisstep,linestylethisstep,imarktype
  use settings_page, only:maxlinestyle,modlinestyle,maxcolour,modcolour
  use plotlib,       only:plotlib_maxlinestyle,plotlib_maxlinecolour
  implicit none
  integer, intent(in) :: istep
  logical, intent(in) :: iChangeColours, iChangeStyles
  integer :: icolour,imarkernumber

  if (iChangeColours) then
     if (allocated(icolourme)) then
        icolour = istep + 1
        icolour = (icolour-2)/modcolour + 1
        if (icolour.gt.plotlib_maxlinecolour) then
           print "(a,i2,a)",'warning: step colour > ',plotlib_maxlinecolour,': re-using colours'
           icolour = mod(icolour-1,plotlib_maxlinecolour) + 1
        endif
        icolour = mod(icolour-1,min(maxcolour,plotlib_maxlinecolour)) + 1
        icolourme = icolour
        linecolourthisstep = icolour
     else
        print "(a)",'***error: array not allocated in colour_timestep***'
     endif
  endif
  if (iChangeStyles) then
     !--modlinestyle should not be greater than max line styles in plotting library
     linestylethisstep = mod((istep-1)/modlinestyle,min(maxlinestyle,plotlib_maxlinestyle)) + 1
     imarkernumber = istep
     select case(imarkernumber)
     case(1)
        imarktype(1) = 4
     case(2)
        imarktype(1) = 17
     case(3)
        imarktype(1) = 2
     case(4)
        imarktype(1) = 3
     case(5:16)
        imarktype(1) = imarkernumber
     case(17:)
        imarktype(1) = imarkernumber + 1
     end select
  endif

  return
end subroutine colour_timestep

!---------------------------------------------------------------------------------------
! colours all the particles using the default colour for their type
!---------------------------------------------------------------------------------------
subroutine colourparts_default(npartoftype,iamtype)
  use params, only:int1
  use settings_data, only:ntypes
  use particle_data, only:icolourme
  use settings_part, only:idefaultcolourtype
  implicit none
  integer, dimension(:), intent(in) :: npartoftype
  integer(kind=int1), dimension(:), intent(in) :: iamtype
  integer :: i,index1,index2,itype

  if (size(iamtype).gt.1) then
     do i=1,sum(npartoftype(1:ntypes))
        itype = iamtype(i)
        if (itype.gt.0 .and. itype.le.ntypes) then
           if (idefaultcolourtype(itype).ge.0) then
              icolourme(i) = idefaultcolourtype(itype)
           endif
        endif
     enddo
  else
     index1 = 1
     do itype=1,ntypes
        index2 = index1 + npartoftype(itype) - 1
        if (idefaultcolourtype(itype).ge.0) then
           icolourme(index1:index2) = idefaultcolourtype(itype)
        endif
        index1 = index2 + 1
     enddo
  endif

end subroutine colourparts_default

end module timestepping
