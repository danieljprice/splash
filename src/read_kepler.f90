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

 ! first see if a global file.comp file exists
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

subroutine read_kepler_composition(filename,ntotal,dat,icomp_col_start,ncomp)
 use asciiutils, only:get_ncolumns
 real, intent(inout) :: dat(:,:)
 character(len=*), intent(in) :: filename
 integer, intent(in) :: ncomp,icomp_col_start
 integer, intent(inout) :: ntotal
 integer :: iu,ierr,i,ncols,nhdr

 open(newunit=iu,file=trim(filename),form='formatted',status='old',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
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
       if (ierr /= 0) print*,' ERROR reading '//trim(filename)//' on line ',i+nhdr
    enddo
 else
    print "(a)",' ERROR: wrong number of columns in '//trim(filename)//'.comp'
 endif

end subroutine read_kepler_composition

end module read_kepler
