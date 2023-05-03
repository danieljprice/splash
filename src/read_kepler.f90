!-----------------------------------------------------------------
!
! This file is for reading the composition file for each particle
! that was read from KEPLER models.
! it also performs a check to see if composition file exists or not
!
!-----------------------------------------------------------------

module read_kepler
 implicit none

 public :: check_for_composition_file,read_kepler_composition

 private

contains

! This function returns the prefix of the filename
function get_prefix(filename) result(prefix)
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: prefix

 integer :: iu

 iu = index(filename,'_')
 if (iu > 1) then
    prefix = filename(1:iu-1)
 else
    prefix = filename
 endif

end function get_prefix

subroutine check_for_composition_file(dumpfile,ntotal,ncolstep,icomp_col_start,ncomp,labels)
 use asciiutils, only:get_ncolumns,get_nrows,read_column_labels
 character(len=*), intent(in) :: dumpfile
 integer, intent(in) :: ntotal
 integer, intent(inout) :: ncolstep
 character(len=*), intent(inout) :: labels(:)
 integer, intent(out) :: ncomp,icomp_col_start

 integer :: iu,nrows,nheaderlines,nlabels,ierr
 character(len=len(dumpfile)) :: prefix
 logical :: iexist

 ncomp = 0
 icomp_col_start = 0

 prefix = get_prefix(dumpfile)

 inquire(file=trim(prefix)//'.comp',exist=iexist)
 if (.not.iexist) return

 ! see if Kepler composition file exists
 open(newunit=iu,file=trim(prefix)//'.comp',iostat=ierr,status='old')
 if (ierr /= 0) then
     print "(a)", 'ERROR opening ', trim(prefix)//'.comp'
 endif

 if (ierr == 0) then
    call get_ncolumns(iu,ncomp,nheaderlines)
    call get_nrows(iu,nheaderlines,nrows)
    ! check nrows equals number of particles
    if (nrows /= ntotal) then
       print*,'ERROR number of rows should equal the number of particles'
       return
    endif

    if (ncomp > 0) then
       icomp_col_start = ncolstep + 1
       ncolstep = ncolstep + ncomp
       call read_column_labels(iu,nheaderlines,ncomp,nlabels,&
            labels(icomp_col_start:icomp_col_start+ncomp-1))
    endif
 endif
 close(iu)

end subroutine check_for_composition_file

subroutine read_kepler_composition(dumpfile,ntotal,dat,icomp_col_start,ncomp)
 use asciiutils, only:get_ncolumns
 real, intent(inout) :: dat(:,:)
 character(len=*), intent(in) :: dumpfile
 integer, intent(in) :: ncomp,icomp_col_start
 integer, intent(inout) :: ntotal
 character(len=len(dumpfile)) :: prefix
 integer :: iu,ierr,i,ncols,nhdr

 prefix = get_prefix(dumpfile)

 open(newunit=iu,file=trim(prefix)//'.comp',form='formatted',status='old',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(prefix)//'.comp'
 elseif (icomp_col_start+ncomp-1 <= size(dat(1,:))) then
    ! get number of columns
    call get_ncolumns(iu,ncols,nhdr)
    ! skip header lines
    do i=1,nhdr
       read(iu,*,iostat=ierr)
    enddo
    ! read data from file
    do i=1,ntotal
       read(iu,*,iostat=ierr) dat(i,icomp_col_start:icomp_col_start+ncomp-1)
       if (ierr /= 0) print*,' ERROR reading '//trim(prefix)//'.comp on line ',i+nhdr
    enddo
 else
    print "(a)",' ERROR: wrong number of columns in '//trim(prefix)//'.comp'
 endif

end subroutine read_kepler_composition

end module read_kepler
