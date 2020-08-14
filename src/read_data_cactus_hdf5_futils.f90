!-------------------------------------------------------------------------
!
!  Utility routines for the cactus hdf5 data read
!  Module contains interface routines to c functions
!  that perform the actual calls to the HDF5 libs
!
!-------------------------------------------------------------------------
module cactushdf5read
 use params, only:maxplot
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none
 character(len=lenlabel), dimension(maxplot) :: blocklabel
 character(len=130) :: datfileprev = ' '
 logical :: file_is_open = .false.
 integer :: ntoti_prev,ncol_prev,nstep_prev

 interface
  subroutine open_cactus_hdf5_file(filename,istep,npart,ncol,nstep_max,ndim,ndimV,time,ignoretl,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in) :: filename
   integer(kind=c_int), intent(in), value :: istep,ignoretl
   integer(kind=c_int), intent(out) :: npart,ncol,nstep_max,ndim,ndimV,ierr
   real(kind=c_double), intent(out) :: time
  end subroutine open_cactus_hdf5_file

  subroutine read_cactus_hdf5_data(filename,istep,npart,time,dx,ignoretl,ierr) bind(c)
   import
   character(kind=c_char), dimension(*), intent(in)  :: filename
   integer(kind=c_int), intent(in), value :: istep,ignoretl
   integer(kind=c_int), intent(out) :: npart,ierr
   real(kind=c_double), intent(out) :: time,dx
  end subroutine read_cactus_hdf5_data

  subroutine close_cactus_hdf5_file(ierr) bind(c)
   import
   integer(kind=c_int), intent(out) :: ierr
  end subroutine close_cactus_hdf5_file

 end interface

contains

!-------------------------------------------------------------------------
!
!  The following routines are callback routines called by the c
!  utilities
!
!-------------------------------------------------------------------------
subroutine set_blocklabel(icol,name,lenname) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 ! use cactushdf5read, only:blocklabel
 use asciiutils,     only:fstring
 implicit none
 integer(kind=c_int),    intent(in) :: icol,lenname
 character(kind=c_char), intent(in) :: name(lenname)
 character(len=24) :: temp
 integer :: ivar

 temp = fstring(name)
 ivar = index(temp,'::')
 if (ivar > 0) temp = temp(ivar+2:)
 if (icol <= size(blocklabel)) then
    blocklabel(icol) = trim(temp)
 else
    print*,'ERROR - too many columns in file'
 endif
 !print*,icol,' name = ',trim(blocklabel(icol))

end subroutine set_blocklabel

subroutine sort_cactus_data(n,iter,iorder) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int
 use sort, only:indexxi
 implicit none
 integer(kind=c_int), intent(in)  :: n
 integer(kind=c_int), intent(in)  :: iter(n)
 integer(kind=c_int), intent(out) :: iorder(n)

 call indexxi(n,iter,iorder)

end subroutine sort_cactus_data

!-------------------------------------------------------------------------
!
!  The following routines compute useful things, e.g. tr K
!
!-------------------------------------------------------------------------
subroutine calc_trK(gxxd,gxyd,gxzd,gyyd,gyzd,gzzd,kxxd,kxyd,kxzd,kyyd,kyzd,kzzd,trk)
 !
 ! Subroutine to calculate trace K ( trK = g^{ij} K_{ij} )
 ! Takes g_{ij} and K_{ij} at one position and time, returns trK
 !
 real, intent(in) :: gxxd,gxyd,gxzd,gyyd,gyzd,gzzd ! spatial down metric components
 real, intent(in) :: kxxd,kxyd,kxzd,kyyd,kyzd,kzzd ! down extrinsic curvature components
 real :: gxxu,gxyu,gxzu,gyyu,gyzu,gzzu             ! spatial up metric components
 real, intent(out) :: trk
 real, dimension(3,3) :: gijd, giju  ! 4x4 metric down and up respectively
 real :: det

 ! down components
 gijd(1,1) = gxxd
 gijd(1,2) = gxyd
 gijd(1,3) = gxzd
 gijd(2,1) = gxyd
 gijd(2,2) = gyyd
 gijd(2,3) = gyzd
 gijd(3,1) = gxzd
 gijd(3,2) = gyzd
 gijd(3,3) = gzzd

 call inv3x3(gijd,giju,det)

 ! up (inverse) components
 gxxu = giju(1,1)
 gxyu = giju(1,2)
 gxzu = giju(1,3)
 gyyu = giju(2,2)
 gyzu = giju(2,3)
 gzzu = giju(3,3)


 trk = (gxxu * kxxd) + (2. * gxyu * kxyd) + (2. * gxzu * kxzd) + &
 & (gyyu * kyyd) + (2. * gyzu * kyzd) + (gzzu * kzzd)

end subroutine calc_trK

pure subroutine inv3x3(A,B,det)
 real, intent(in), dimension(3,3) :: A
 real, intent(out), dimension(3,3) :: B ! inverse matrix
 real, intent(out) :: det

 det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - &
 & A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + &
 & A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

 B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
 B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
 B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
 B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
 B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
 B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
 B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
 B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
 B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

 B(:,:) = B(:,:)/det

end subroutine inv3x3

 !
 ! find location of metric and extrinsic curvature in columns
 !
subroutine find_metric(ncols,labelcol,igxx,igxy,igxz,igyy,igyz,igzz,ikxx,ikxy,ikxz,ikyy,ikyz,ikzz,&
      irho,ialp,ivel0,ivel1,ivel2,gotrho,gotalp)
 integer,          intent(in) :: ncols
 character(len=*), intent(in) :: labelcol(ncols)
 integer, intent(out) :: igxx,igxy,igxz,igyy,igyz,igzz
 integer, intent(out) :: ikxx,ikxy,ikxz,ikyy,ikyz,ikzz
 integer, intent(out) :: irho,ialp,ivel0,ivel1,ivel2
 integer :: i, idens
 logical, intent(out) :: gotrho, gotalp
 gotrho = .False.; gotalp = .False.

 igxx = 0; igxy = 0; igxz = 0; igyy = 0; igyz = 0; igzz = 0
 ikxx = 0; ikxy = 0; ikxz = 0; ikyy = 0; ikyz = 0; ikzz = 0
 irho = 0; ialp = 0; ivel0 = 0; ivel1 = 0; ivel2 = 0
 do i=1,ncols
    select case(labelcol(i))
    case('gxx')
       igxx = i
    case('gxy')
       igxy = i
    case('gxz')
       igxz = i
    case('gyy')
       igyy = i
    case('gyz')
       igyz = i
    case('gzz')
       igzz = i
    case('kxx')
       ikxx = i
    case('kxy')
       ikxy = i
    case('kxz')
       ikxz = i
    case('kyy')
       ikyy = i
    case('kyz')
       ikyz = i
    case('kzz')
       ikzz = i
    case('rho')
       irho = i
       gotrho = .True.
    case('dens')
       idens = i
    case('alp')
       ialp = i
       gotalp = .True.
    case('vel[0]')
       ivel0 = i
    case('vel[1]')
       ivel1 = i
    case('vel[2]')
       ivel2 = i
    end select
 enddo

 ! if we didn't find 'rho' in cols, use dens instead
 if (gotrho .eqv. .False.) irho = idens

end subroutine find_metric

subroutine compute_extra_columns(ncols,nextra,dat)
 integer, intent(in)  :: ncols
 integer, intent(out) :: nextra
 real, intent(inout), optional  :: dat(:,:)
 integer :: i,n,itrk
 integer :: igxx,igxy,igxz,igyy,igyz,igzz
 integer :: ikxx,ikxy,ikxz,ikyy,ikyz,ikzz
 integer :: irho,ialp,ivel0,ivel1,ivel2
 logical :: gotrho,gotalp

 nextra = 0
 call find_metric(ncols,blocklabel,igxx,igxy,igxz,igyy,igyz,igzz,ikxx,ikxy,ikxz,ikyy,ikyz,ikzz,&
       irho,ialp,ivel0,ivel1,ivel2,gotrho,gotalp)
 n = size(dat(:,1))
 if (igxx > 0 .and. igxy > 0 .and. igxz > 0 .and. igyy > 0 .and. igyz > 0 .and. igzz > 0 .and. &
      ikxx > 0 .and. ikxy > 0 .and. ikxz > 0 .and. ikyy > 0 .and. ikyz > 0 .and. ikzz > 0) then
    itrk = ncols + 1
    blocklabel(itrk) = 'tr K'
    nextra = 1
    if (present(dat)) then
       !print*,' getting trk ',igxx,igxy,igxz,igyy,igyz,igzz,ikxx,ikxy,ikxz,ikyy,ikyz,ikzz,itrk
       do i=1,n
          call calc_trK(dat(i,igxx),dat(i,igxy),dat(i,igxz),dat(i,igyy),dat(i,igyz),dat(i,igzz),&
                         dat(i,ikxx),dat(i,ikxy),dat(i,ikxy),dat(i,ikyy),dat(i,ikyz),dat(i,ikzz),dat(i,itrk))
       enddo
    endif
 endif

end subroutine compute_extra_columns

end module cactushdf5read
