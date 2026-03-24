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
!  Copyright (C) 2005-2025 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------!-----------------------------------------------------------------
!
! This file is for reading the composition file for each particle
! it also performs a check to see if composition file exists or not
!
!-----------------------------------------------------------------
module readcomposition
 implicit none

 public :: check_for_composition_file,read_composition

 private

contains

! This function returns the prefix of the filename
function get_prefix(filename) result(prefix)
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: prefix

 integer :: iu

 iu = index(filename,'_',back=.true.)
 if (iu > 1) then
    prefix = filename(1:iu-1)
 else
    prefix = filename
 endif

end function get_prefix

subroutine check_for_composition_file(dumpfile,ntotal,ncolstep,icomp_col_start,ncomp,labels,filename)
 use asciiutils,    only:get_ncolumns,get_nrows,read_column_labels,basename
 use settings_data, only:debugmode
 character(len=*), intent(in) :: dumpfile
 integer, intent(in) :: ntotal
 integer, intent(inout) :: ncolstep
 character(len=*), intent(inout) :: labels(:)
 integer, intent(out) :: ncomp,icomp_col_start
 character(len=*), intent(out) :: filename

 integer :: iu,nrows,nheaderlines,nlabels,ierr
 character(len=len(dumpfile)) :: prefix

 logical :: iexist

 ncomp = 0
 icomp_col_start = 0

 ! first see if file_00000.cols exists
 filename = trim(dumpfile)//'.cols'
 inquire(file=filename,exist=iexist)
 if (debugmode) print*,' DEBUG: looking for '//trim(filename)//' exist = ',iexist

  ! then see if a local file_00000.comp file exists
 if (.not.iexist) then
    filename = trim(dumpfile)//'.comp'
    inquire(file=filename,exist=iexist)
    if (debugmode) print*,' DEBUG: looking for '//trim(filename)//' exist = ',iexist
 endif

 ! or see if a global file.comp file exists
 if (.not.iexist) then
    prefix = get_prefix(dumpfile)
    filename = trim(prefix)//'.comp'
    inquire(file=filename,exist=iexist)
    if (debugmode) print*,' DEBUG: looking for '//trim(filename)//' exist = ',iexist
 endif

 ! look for a file in the current directory, not the same directory as the data
 if (.not.iexist) then
    filename = basename(filename)
    inquire(file=filename,exist=iexist)
    if (debugmode) print*,' DEBUG: looking for '//trim(filename)//' exist = ',iexist
 endif

 ! look for an AV_00100 file to match the dumpfile number
 if (.not.iexist) then
    filename = get_av_filename(dumpfile)
    inquire(file=filename,exist=iexist)
    if (debugmode) print*,' DEBUG: looking for '//trim(filename)//' exist = ',iexist
    if (iexist) then
       icomp_col_start = ncolstep + 1
       ncolstep = ncolstep + 1
       ncomp = ncomp + 1
       labels(icomp_col_start) = 'A_V'
       print "(a,i0,a)", '> got ',ncomp,' extra columns from '//trim(filename)
       return
    endif
 endif

 !if (debugmode) print*,'DEBUG: looking for '//trim(filename)//' exist = ',iexist
 if (.not.iexist) return

 ! see if Kepler composition file exists
 open(newunit=iu,file=filename,iostat=ierr,status='old')
 if (ierr /= 0) then
    print "(1x,a)", 'ERROR opening '//trim(filename)
 endif

 if (ierr == 0) then
    call get_ncolumns(iu,ncomp,nheaderlines)
    call get_nrows(iu,nheaderlines,nrows)
    ! check nrows equals number of particles
    if (nrows /= ntotal) then
       print "(1x,a,i0,a,i0,a)",'ERROR: reading '//trim(filename)//' nrows (',nrows,') /= nparticles (',ntotal,')'
       return
    endif
    if (ncomp > 0) then
       icomp_col_start = ncolstep + 1
       ncolstep = ncolstep + ncomp
       call read_column_labels(iu,nheaderlines,ncomp,nlabels,&
            labels(icomp_col_start:icomp_col_start+ncomp-1))
       print "(a,i0,a)", '> got ',ncomp,' extra columns from '//trim(filename)
    endif
 endif
 close(iu)

end subroutine check_for_composition_file

subroutine read_composition(filename,ntotal,dat,icomp_col_start,ncomp,iorig)
 use asciiutils,      only:get_ncolumns
 use params,          only:int8
 real, intent(inout) :: dat(:,:)
 character(len=*), intent(in) :: filename
 integer, intent(in) :: ncomp,icomp_col_start
 integer, intent(inout) :: ntotal
 integer(kind=int8), allocatable, intent(in) :: iorig(:)
 integer :: iu,ierr,i,ncols,nhdr

 open(newunit=iu,file=trim(filename),form='formatted',status='old',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
 elseif (icomp_col_start+ncomp-1 <= size(dat(1,:))) then
    if (index(filename,'AV_') > 0) then
       call read_av_into_column(filename,ntotal,dat,icomp_col_start,iorig)
    else
       ! get number of columns
       call get_ncolumns(iu,ncols,nhdr)
       ! skip header lines
       do i=1,nhdr
          read(iu,*,iostat=ierr)
       enddo
       ! read data from file
       do i=1,ntotal
          read(iu,*,iostat=ierr) dat(i,icomp_col_start:icomp_col_start+ncomp-1)
          if (ierr /= 0) print*,' ERROR reading '//trim(filename)//' on line ',i+nhdr
       enddo
    endif
 else
    print "(a)",' ERROR: wrong number of columns in '//trim(filename)//'.comp'
 endif

end subroutine read_composition

!-----------------------------------------------------------------
! This function returns the filename of the AV file
! for the given dumpfile, e.g. dump_00700 -> AV_00700
!-----------------------------------------------------------------
function get_av_filename(dumpfile) result(avfile)
 character(len=*), intent(in) :: dumpfile
 character(len=len(dumpfile)+4) :: avfile
 integer :: i,islash

 i = index(trim(dumpfile), '_', back=.true.)
 avfile = ''
 if (i > 0 .and. i < len_trim(dumpfile)) then
    islash = index(trim(dumpfile), '/', back=.true.)
    if (islash > 0) then
       avfile = dumpfile(1:islash)//'AV_'//dumpfile(i+1:len_trim(dumpfile))
    else
       avfile = 'AV_'//dumpfile(i+1:len_trim(dumpfile))
    endif
 endif

end function get_av_filename

!-----------------------------------------------------------------
! Read the AV file (stream unformatted: npart_av, id(:), val(:)) into one column.
! id_tmp in the file are original particle IDs. If iorig is allocated, match by
! iorig(particle_index); otherwise id_tmp is treated as dump index 1..ntotal.
!-----------------------------------------------------------------
subroutine read_av_into_column(filename,ntotal,dat,icol,iorig)
 use params, only:doub_prec,int8
 real, intent(inout) :: dat(:,:)
 character(len=*), intent(in) :: filename
 integer, intent(in) :: icol
 integer, intent(in) :: ntotal
 integer(kind=int8), allocatable, intent(in) :: iorig(:)
 integer :: iu,ierr,npart_av,k,max_id,orig_id
 integer, allocatable :: id_tmp(:)
 real(doub_prec), allocatable :: val_tmp(:), lookup(:)

 dat(1:ntotal,icol) = 0.0

 ! open the AV file
 open(newunit=iu,file=trim(filename),access='stream',form='unformatted',status='old',iostat=ierr)
 if (ierr /= 0) then
    print*, "WARNING: Could not open A_V file "//trim(filename)//", using A_V = 0"
    return
 endif

 ! read number of particles
 read(iu,iostat=ierr) npart_av
 if (ierr /= 0 .or. npart_av <= 0) then
    print*, "WARNING: Invalid A_V header in "//trim(filename)//", using A_V = 0"
    close(iu)
    return
 endif

 ! read arrays
 allocate(id_tmp(npart_av), val_tmp(npart_av))
 read(iu,iostat=ierr) id_tmp
 read(iu,iostat=ierr) val_tmp
 close(iu)

 if (ierr /= 0) then
    print*, "WARNING: Error reading A_V data from "//trim(filename)//", using A_V = 0"
    deallocate(id_tmp, val_tmp)
    return
 endif

 max_id = maxval(id_tmp(1:npart_av))
 if (max_id < 1) then
    deallocate(id_tmp, val_tmp)
    return
 endif

 allocate(lookup(max_id))
 lookup = 0.0d0
 do k = 1, npart_av
    if (id_tmp(k) >= 1 .and. id_tmp(k) <= max_id) then
       lookup(id_tmp(k)) = val_tmp(k)
    endif
 enddo

 if (allocated(iorig) .and. size(iorig) >= ntotal) then
   print "(a,i0,a,i0,a)", ' Loaded A_V for ',npart_av,' particles from '//&
                          trim(filename)//' using ',max_id,' unique IDs'
   do k = 1, ntotal
       orig_id = int(iorig(k), kind(orig_id))
       if (orig_id >= 1 .and. orig_id <= max_id) then
          dat(k,icol) = real(lookup(orig_id))
       endif
    enddo
 else
    print "(a,i0,a)", ' Loaded A_V for ',npart_av,' particles from '// &
                      trim(filename)//' (no iorig found in dump)'
    do k = 1, ntotal
       if (k <= max_id) dat(k,icol) = real(lookup(k))
    enddo
 endif
 deallocate(id_tmp, val_tmp, lookup)

end subroutine read_av_into_column

end module readcomposition
