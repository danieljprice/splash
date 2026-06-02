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
!-----------------------------------------------------------------
!
! HDF5 extra columns for Phantom binary dumps (requires HDF5=yes build).
! Reads chemistry sidecar files (dump_XXXXX.h5) via read_phantom_hdf5_utils.c.
!
! C reads HDF5 and calls back into Fortran with no access to subroutine
! arguments. dat_slice_ptr and iorig_ptr bridge the caller's arrays; C passes
! icomp_col_start into read_extra_column_fromc_phantom for column indexing.
!
!-----------------------------------------------------------------
module readcomposition_hdf5
 use labelschem, only:format_chemistry_label
 use params, only:maxplot,int8
 use labels, only:lenlabel,ih,label
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none

 integer :: h5_npart_save = 0
 integer :: ntotal_save = 0
 integer, allocatable :: h5_particle_ids(:)
 integer(kind=int8), dimension(:), pointer :: iorig_ptr => null()
 real, dimension(:,:), pointer :: dat_slice_ptr => null()
 logical, save :: use_live_scatter = .false.
 logical, save :: scatter_mode_set = .false.
 character(len=lenlabel), dimension(maxplot) :: extra_labels_buf

 interface
  subroutine phantom_hdf5_extra_check(filename,ntotal,ncomp,npart_file,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   integer(kind=c_int), intent(inout) :: ntotal
   integer(kind=c_int), intent(out) :: ncomp,npart_file,ierr
  end subroutine phantom_hdf5_extra_check

  subroutine phantom_hdf5_read_particle_ids(filename,npart,ids,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   integer(kind=c_int), intent(in) :: npart
   integer(kind=c_int), dimension(npart), intent(out) :: ids
   integer(kind=c_int), intent(out) :: ierr
  end subroutine phantom_hdf5_read_particle_ids

  subroutine phantom_hdf5_extra_read(filename,ntotal,icomp_col_start,ncomp,isrequired,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   integer(kind=c_int), intent(in) :: ntotal,icomp_col_start,ncomp
   integer(kind=c_int), intent(out) :: ierr
   integer(kind=c_int), dimension(ncomp), intent(in) :: isrequired
  end subroutine phantom_hdf5_extra_read
 end interface

 public :: check_for_composition_hdf5, read_composition_hdf5

 private :: apply_extra_labels_from_buf

contains

!-----------------------------------------------------------------
! check for Phantom HDF5 sidecar; count species and set column labels
!-----------------------------------------------------------------
subroutine check_for_composition_hdf5(dumpfile,ntotal,ncolstep,icomp_col_start,ncomp,labels,filename,ierr)
 use asciiutils, only:cstring
 use params, only:maxplot
 character(len=*), intent(in) :: dumpfile
 integer, intent(in) :: ntotal
 integer, intent(inout) :: ncolstep
 character(len=*), intent(inout) :: labels(:)
 integer, intent(out) :: ncomp,icomp_col_start,ierr
 character(len=*), intent(out) :: filename
 integer(c_int) :: ntotal_c,ncomp_c,npart_c,ierr_c

 ncomp = 0
 icomp_col_start = 0
 ierr = 0
 h5_npart_save = 0
 filename = trim(dumpfile)//'.h5'

 ntotal_c = ntotal
 extra_labels_buf = ' '
 call phantom_hdf5_extra_check(cstring(filename),ntotal_c,ncomp_c,npart_c,ierr_c)
 if (ierr_c /= 0) then
    ierr = ierr_c
    return
 endif
 if (ncomp_c <= 0) return

 h5_npart_save = npart_c
 icomp_col_start = ncolstep + 1
 if (ncolstep + ncomp_c > maxplot) then
    print "(1x,a,i0,a,i0,a)",'ERROR: too many columns (',ncolstep+ncomp_c,' > maxplot=',maxplot,') from '//trim(filename)
    icomp_col_start = 0
    h5_npart_save = 0
    return
 endif

 ncomp = ncomp_c
 ncolstep = icomp_col_start + ncomp - 1
 call apply_extra_labels_from_buf(labels, icomp_col_start, ncomp)
 print "(a,i0,a)",'> got ',ncomp,' extra columns from '//trim(filename)

end subroutine check_for_composition_hdf5

!-----------------------------------------------------------------
! read extra composition columns from Phantom HDF5 sidecar into dat
!-----------------------------------------------------------------
subroutine read_composition_hdf5(filename,ntotal,dat,icomp_col_start,ncomp,required_mask,ierr,iorig)
 use asciiutils, only:cstring
 real, intent(inout), target :: dat(:,:)
 character(len=*), intent(in) :: filename
 integer, intent(in) :: ncomp,icomp_col_start,ntotal
 logical, intent(in) :: required_mask(:)
 integer, intent(out) :: ierr
 integer(kind=int8), target, intent(in), optional :: iorig(:)
 integer :: i,nreq
 integer, dimension(:), allocatable :: isrequired
 integer(c_int) :: ntotal_c,icomp_c,ncomp_c,ierr_c,ierr_ids

 ierr = 0
 if (ncomp <= 0) return

 nreq = min(ncomp,size(required_mask))
 if (nreq < 1 .or. .not.any(required_mask(1:nreq))) return

 allocate(isrequired(ncomp))
 isrequired = 0
 do i=1,nreq
    if (required_mask(i)) isrequired(i) = 1
 enddo

 ntotal_save = ntotal
 dat_slice_ptr => dat
 scatter_mode_set = .false.
 use_live_scatter = .false.

 if (associated(iorig_ptr)) nullify(iorig_ptr)
 if (present(iorig) .and. size(iorig) >= ntotal) iorig_ptr => iorig

 if (h5_npart_save > 0) then
    allocate(h5_particle_ids(h5_npart_save))
    call phantom_hdf5_read_particle_ids(cstring(filename),h5_npart_save,h5_particle_ids,ierr_ids)
    if (ierr_ids /= 0) then
       print "(a)",' WARNING: could not read particle ids from '//trim(filename)//'; using row order'
       deallocate(h5_particle_ids)
    endif
 endif

 do i=1,ncomp
    if (i > size(required_mask) .or. .not.required_mask(i)) cycle
    if (icomp_col_start+i-1 >= 1 .and. icomp_col_start+i-1 <= size(dat,2)) then
       dat(:,icomp_col_start+i-1) = 0.
    endif
 enddo

 ntotal_c = ntotal
 icomp_c = icomp_col_start
 ncomp_c = ncomp
 call phantom_hdf5_extra_read(cstring(filename),ntotal_c,icomp_c,ncomp_c,isrequired,ierr_c)
 ierr = ierr_c

 if (allocated(h5_particle_ids)) deallocate(h5_particle_ids)
 nullify(iorig_ptr)
 nullify(dat_slice_ptr)
 ntotal_save = 0

 deallocate(isrequired)

end subroutine read_composition_hdf5

!-----------------------------------------------------------------
! scatter one HDF5 sidecar column into dat by original particle id
!-----------------------------------------------------------------
subroutine scatter_extra_column_by_id(icolput, icomp_col_start, npart_h5, values)
 use params, only:sing_prec
 integer, intent(in) :: icolput,icomp_col_start,npart_h5
 real(kind=c_double), intent(in) :: values(npart_h5)
 integer :: k,p,orig_id,max_id,nmatched
 real(kind=sing_prec), allocatable :: lookup(:)
 logical, allocatable :: id_present(:)

 if (icolput < 1 .or. .not.associated(dat_slice_ptr)) return
 if (icolput > size(dat_slice_ptr,2)) return
 if (npart_h5 <= 0 .or. ntotal_save <= 0) return

 if (use_live_scatter .or. .not.allocated(h5_particle_ids) .or. .not.associated(iorig_ptr) &
     .or. npart_h5 > size(h5_particle_ids)) then
    call scatter_extra_by_live_particle_order(icolput, icomp_col_start, npart_h5, values)
    return
 endif

 max_id = maxval(h5_particle_ids(1:npart_h5))
 if (max_id < 1) return
 allocate(lookup(max_id),id_present(max_id))
 lookup = 0.
 id_present = .false.
 do k=1,npart_h5
    if (h5_particle_ids(k) >= 1 .and. h5_particle_ids(k) <= max_id) then
       lookup(h5_particle_ids(k)) = real(values(k), kind=sing_prec)
       id_present(h5_particle_ids(k)) = .true.
    endif
 enddo

 if (.not.scatter_mode_set) then
    nmatched = 0
    do p=1,ntotal_save
       orig_id = int(iorig_ptr(p))
       if (orig_id >= 1 .and. orig_id <= max_id) then
          if (id_present(orig_id)) nmatched = nmatched + 1
       endif
    enddo
    use_live_scatter = (nmatched < npart_h5/2)
    scatter_mode_set = .true.
    if (use_live_scatter) then
       print "(a,i0,a,i0,a)",' WARNING: HDF5 sidecar id match low (',nmatched,'/',npart_h5,&
          '); using live-particle order for chemistry columns'
       deallocate(lookup,id_present)
       call scatter_extra_by_live_particle_order(icolput, icomp_col_start, npart_h5, values)
       return
    endif
 endif

 do p=1,ntotal_save
    orig_id = int(iorig_ptr(p))
    if (orig_id >= 1 .and. orig_id <= max_id .and. id_present(orig_id)) then
       dat_slice_ptr(p,icolput) = real(lookup(orig_id), kind=kind(dat_slice_ptr(1,1)))
    endif
 enddo
 deallocate(lookup,id_present)

end subroutine scatter_extra_column_by_id

!-----------------------------------------------------------------
! scatter one HDF5 sidecar column using live-particle order in the dump
!-----------------------------------------------------------------
subroutine scatter_extra_by_live_particle_order(icolput, icomp_col_start, npart_h5, values)
 integer, intent(in) :: icolput,icomp_col_start,npart_h5
 real(kind=c_double), intent(in) :: values(npart_h5)
 integer :: p,nlive

 if (icolput < 1 .or. .not.associated(dat_slice_ptr)) return
 if (icolput > size(dat_slice_ptr,2)) return
 nlive = 0
 do p=1,ntotal_save
    if (ih > 0 .and. (icomp_col_start <= 0 .or. ih < icomp_col_start) &
        .and. ih <= size(dat_slice_ptr,2) .and. dat_slice_ptr(p,ih) <= 0.) cycle
    nlive = nlive + 1
    if (nlive > npart_h5) exit
    dat_slice_ptr(p,icolput) = real(values(nlive), kind=kind(dat_slice_ptr(1,1)))
 enddo

end subroutine scatter_extra_by_live_particle_order

!-----------------------------------------------------------------
! copy species names from extra_labels_buf into tagarr and label
!-----------------------------------------------------------------
subroutine apply_extra_labels_from_buf(tagarr, icomp_col_start, ncomp)
 character(len=*), intent(inout) :: tagarr(:)
 integer, intent(in) :: icomp_col_start,ncomp
 integer :: i,icol

 if (icomp_col_start < 1 .or. ncomp < 1) return
 if (icomp_col_start + ncomp - 1 > size(tagarr)) return
 do i=1,ncomp
    icol = icomp_col_start + i - 1
    if (icol > size(tagarr)) exit
    if (len_trim(extra_labels_buf(i)) == 0) cycle
    if (icol <= size(label)) label(icol) = trim(extra_labels_buf(i))
    tagarr(icol) = trim(extra_labels_buf(i))
 enddo

end subroutine apply_extra_labels_from_buf

!-----------------------------------------------------------------
! C callback: store formatted sidecar species name for plot labels
!-----------------------------------------------------------------
subroutine set_extra_column_label_phantom(icol,name) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int,c_char
 use asciiutils, only:fstring
 integer(kind=c_int), intent(in) :: icol
 character(kind=c_char), intent(in) :: name(256)

 if (icol >= 1 .and. icol <= size(extra_labels_buf)) then
    extra_labels_buf(icol) = format_chemistry_label(fstring(name))
 endif

end subroutine set_extra_column_label_phantom

!-----------------------------------------------------------------
! C callback: scatter one sidecar column read from HDF5 into dat
!-----------------------------------------------------------------
subroutine read_extra_column_fromc_phantom(icol,npart,temparr,icomp_col_start) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int,c_double
 integer(kind=c_int), intent(in) :: icol,npart,icomp_col_start
 real(kind=c_double), intent(in) :: temparr(npart)
 integer :: icolput

 icolput = icomp_col_start + icol - 1
 call scatter_extra_column_by_id(icolput, icomp_col_start, npart, temparr)

end subroutine read_extra_column_fromc_phantom

end module readcomposition_hdf5
