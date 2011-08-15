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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!------------------------------------------------------------------------
!
! module containing subroutines for applying transformations to data
! prior to plotting.
!
! the transformations are:
!
!  1) log10(x)
!  2) |x|
!  3) 1/x
!  4) sqrt(x)
!  5) x^2
!  6) ln(x)
!
! * combinations of transformations are done when the
!   input number is > 10 (e.g. 321 means 1/x, then abs, then log10)
!
! subroutines contained within this module are the following:
!
!   transform        : applies transformation to a one dimensional array
!   transform2       : applies the transformation to a two dimensional array
!   transform_limits : transforms the plot limits appropriately
!   transform_label  : changes the plot label appropriately
!
! Written for use in supersphplot/splash by
! Daniel Price, Institute of Astronomy, Cambridge, UK 2002-2004
!               University of Exeter, UK              2004-
! dprice@astro.ex.ac.uk
!
!------------------------------------------------------------------------
module transforms
 integer, parameter, public :: ntrans = 6  ! this is the number of different transformations
 real, parameter, private :: zerolog = 1.e-12 ! this is minimum set if xmin = 0 and log
 public :: transform,transform_inverse,trans
 public :: transform_limits,transform_limits_inverse,transform_label
 public :: convert_to_ln_fac
 public :: islogged

 private

 interface transform
   module procedure transform,transform_limits,transform2,transforma
 end interface transform

 !--function interface
 interface trans
   module procedure transformarray
 end interface trans

 interface transform_inverse
   module procedure transform_inverse,transform_limits_inverse,transforma_inverse
 end interface transform_inverse

contains
!------------------------------------------------------------------------
!
!  subroutine returns log, 1/x of a given array
!
!  * can specify up to 9 individual operations to perform
!  * combinations of transformations are done when the
!    input number is > 10 (e.g. 321 means 1/x, then abs, then log10)
!
!------------------------------------------------------------------------
subroutine transform(array,itrans,errval)
  implicit none
  integer, intent(in) :: itrans
  real, dimension(:), intent(inout) :: array
  real, intent(in), optional :: errval
  real, dimension(size(array)) :: arraytemp
  real :: errvali
  character(len=20) :: string
  integer :: i
  !
  !--errval is the value to be set for errors
  !  (default is zero if not present)
  !
  if (present(errval)) then
     errvali = errval
  else
     errvali = 0.
  endif
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then

     write(string,*) itrans
     !
     !--do a transformation for each digit
     !
     arraytemp = array

     do i=1,len_trim(string)
        !
        !--perform transformation appropriate to this digit
        !
        select case(string(i:i))
        case('1')
           where (arraytemp > 0. .and. arraytemp.ne.errvali)
              arraytemp = log10(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        case('2')
           where (arraytemp.ne.errvali)
              arraytemp = abs(arraytemp)
           end where
        case('3')
           where (arraytemp .ne. 0. .and. arraytemp.ne.errvali)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = errvali
           end where
        case('4')
           where (arraytemp .gt. 0. .and. arraytemp.ne.errvali)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        case('5')
           where (arraytemp.ne.errvali)
              arraytemp = arraytemp**2
           end where
        case('6')
           where (arraytemp > 0. .and. arraytemp.ne.errvali)
              arraytemp = log(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        end select
     enddo

     array = arraytemp
  endif

end subroutine transform

!------------------------------------------------------------------------
!
!  interface to above for a single real number
!
!------------------------------------------------------------------------
subroutine transforma(aa,itrans,errval)
  implicit none
  integer, intent(in) :: itrans
  real, intent(inout) :: aa
  real, intent(in), optional :: errval
  real, dimension(1) :: array

  array(1) = aa
  if (present(errval)) then
     call transform(array,itrans,errval=errval)
  else
     call transform(array,itrans)
  endif
  aa = array(1)

  return
end subroutine transforma

!------------------------------------------------------------------------
!
!  function interface: returns array valued function
!
!------------------------------------------------------------------------
function transformarray(array,itrans,errval)
  implicit none
  integer, intent(in) :: itrans
  real, intent(in), dimension(:) :: array
  real, dimension(size(array)) :: transformarray
  real, intent(in), optional :: errval

  transformarray = array
  if (present(errval)) then
     call transform(transformarray,itrans,errval=errval)
  else
     call transform(transformarray,itrans)
  endif

  return
end function transformarray

!------------------------------------------------------------------------
!
!  inverse transform
!
!------------------------------------------------------------------------
subroutine transform_inverse(array,itrans,errval)
  implicit none
  integer, intent(in) :: itrans
  real, dimension(:), intent(inout) :: array
  real, intent(in), optional :: errval
  real, dimension(size(array)) :: arraytemp
  real :: errvali
  character(len=20) :: string
  integer :: i
  !
  !--errval is the value to be set for errors
  !  (default is zero if not present)
  !
  if (present(errval)) then
     errvali = errval
  else
     errvali = 0.
  endif
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then

     write(string,*) itrans
     !
     !--do a transformation for each digit
     !
     arraytemp = array

     do i=len_trim(string),1,-1
        !
        !--perform transformation appropriate to this digit
        !
        select case(string(i:i))
        case('1')
           where (arraytemp.ne.errvali)
              arraytemp = 10**arraytemp
           end where
        case('3')
           where (arraytemp .ne. 0. .and. arraytemp.ne.errvali)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = errvali
           end where
        case('4')
           where (arraytemp.ne.errvali)
              arraytemp = arraytemp**2
           end where
        case('5')
           where (arraytemp .gt. 0. .and. arraytemp.ne.errvali)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        case('6')
           where (arraytemp.ne.errvali)
              arraytemp = exp(arraytemp)
           end where
        end select
     enddo

     array = arraytemp
  endif

end subroutine transform_inverse


!------------------------------------------------------------------------
!
!  interface to above for a single real number
!
!------------------------------------------------------------------------
subroutine transforma_inverse(aa,itrans,errval)
  implicit none
  integer, intent(in) :: itrans
  real, intent(inout) :: aa
  real, intent(in), optional :: errval
  real, dimension(1) :: array

  array(1) = aa
  if (present(errval)) then
     call transform_inverse(array,itrans,errval=errval)
  else
     call transform_inverse(array,itrans)
  endif
  aa = array(1)

  return
end subroutine transforma_inverse

!------------------------------------------------------------------------
!
!  same as transform but for a two dimensional array
!  applies the transformation to the same array as was input
!
!------------------------------------------------------------------------
subroutine transform2(array,itrans,errval)
  implicit none
  integer, intent(in) :: itrans
  real, dimension(:,:), intent(inout) :: array
  real, intent(in), optional :: errval
  real, dimension(size(array(:,1)),size(array(1,:))) :: arraytemp
  real :: errvali
  character(len=20) :: string
  integer :: i
  !
  !--errval is the value to be set for errors
  !  (default is zero if not present)
  !
  if (present(errval)) then
     errvali = errval
  else
     errvali = 0.
  endif
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then

     write(string,*) itrans
     !
     !--do a transformation for each digit
     !
     arraytemp = array

     do i=1,len_trim(string)
        !
        !--perform transformation appropriate to this digit
        !
        select case(string(i:i))
        case('1')
           where (arraytemp > 0. .and. arraytemp.ne.errvali)
              arraytemp = log10(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        case('2')
           where (arraytemp.ne.errvali)
              arraytemp = abs(arraytemp)
           end where
        case('3')
           where (arraytemp .ne. 0. .and. arraytemp.ne.errvali)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = errvali
           end where
        case('4')
           where (arraytemp .gt. 0. .and. arraytemp.ne.errvali)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        case('5')
           where (arraytemp.ne.errvali)
              arraytemp = arraytemp**2
           end where
        case('6')
           where (arraytemp > 0. .and. arraytemp.ne.errvali)
              arraytemp = log(arraytemp)
           elsewhere
              arraytemp = errvali
           end where
        end select
     enddo

     array = arraytemp

  endif

end subroutine transform2

!------------------------------------------------------------------------
!
!  same as transform but for the plot limits
!  (min can become max and vice versa)
!
!------------------------------------------------------------------------
subroutine transform_limits(xmin,xmax,itrans)
  implicit none
  integer, intent(in) :: itrans
  real, intent(inout) :: xmin,xmax
  real :: xmintemp,xmaxtemp
  character(len=20) :: string
  integer :: i
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then

     write(string,*) itrans
     !
     !--do a transformation for each digit
     !
     xmintemp = xmin
     xmaxtemp = xmax

     do i=1,len_trim(string)
        !
        !--perform transformation appropriate to this digit
        !
        select case(string(i:i))
        case('1')
           if (xmintemp > 0) then
              xmintemp = log10(xmintemp)
           elseif (xmintemp.eq.0) then
              print*,' log10(xmin = 0): min set to ',zerolog
              xmintemp = log10(zerolog)
           endif
           if (xmaxtemp > 0) then
              xmaxtemp = log10(xmaxtemp)
           elseif (xmaxtemp.eq.0) then
              print*,' log10(xmax = 0): max set to ',zerolog
              xmaxtemp = log10(zerolog)
           endif
        case('2')
           if ((xmintemp.lt.0. .and. xmaxtemp.gt.0.) &
           .or.(xmaxtemp.lt.0. .and. xmintemp.gt.0.)) then
           !
           !--minimum is zero if limits have opposite signs
           !
              xmaxtemp = max(abs(xmintemp),abs(xmaxtemp))
              xmintemp = 0.
           else
           !
           !--or just take magnitude
           !
              xmintemp = abs(xmintemp)
              xmaxtemp = abs(xmaxtemp)
           endif
        case('3')
           if (xmintemp .ne. 0) then
              xmintemp = 1./xmintemp
           else
              xmintemp = 0.
           endif
           if (xmaxtemp .ne. 0) then
              xmaxtemp = 1./xmaxtemp
           else
              xmaxtemp = 0.
           endif
        case('4')
           if (xmintemp .ge. 0) then
              xmintemp = sqrt(xmintemp)
           else
              xmintemp = 0.
           endif
           if (xmaxtemp .ge. 0) then
              xmaxtemp = sqrt(xmaxtemp)
           else
              xmaxtemp = 0.
           endif
        case('5')
           xmintemp = xmintemp**2
           xmaxtemp = xmaxtemp**2
        case('6')
           if (xmintemp > 0) then
              xmintemp = log(xmintemp)
           elseif (xmintemp.eq.0) then
              print*,' ln(xmin = 0): min set to ',zerolog
              xmintemp = log(zerolog)
           endif
           if (xmaxtemp > 0) then
              xmaxtemp = log(xmaxtemp)
           elseif (xmaxtemp.eq.0) then
              print*,' ln(xmax = 0): max set to ',zerolog
              xmaxtemp = log(zerolog)
           endif
        end select
     enddo

     xmin = min(xmintemp,xmaxtemp)
     xmax = max(xmintemp,xmaxtemp)

  endif

end subroutine transform_limits

!------------------------------------------------------------------------
!
!  inverse transform for the plot limits
!  (so that we can change the transformed limits and set the
!   untransformed limits accordingly)
!
!------------------------------------------------------------------------
subroutine transform_limits_inverse(xmin,xmax,itrans)
  implicit none
  integer, intent(in) :: itrans
  real, intent(inout) :: xmin,xmax
  real :: xmintemp,xmaxtemp,xtemp
  character(len=20) :: string
  integer :: i
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then

     write(string,*) itrans
     !
     !--do a transformation for each digit
     !
     xmintemp = xmin
     xmaxtemp = xmax

     do i=len_trim(string),1,-1  ! do digits in reverse
        !
        !--perform transformation appropriate to this digit
        !
        select case(string(i:i))
        case('1')
           xmintemp = 10**xmintemp
           xmaxtemp = 10**xmaxtemp
        case('2')
           !
           !--if minimum is zero give limits opposite signs
           !  (but same magnitude), otherwise do nothing
           !
           if (xmintemp.eq.0.) then
              xtemp = max(abs(xmintemp),abs(xmaxtemp))
              xmintemp = -xtemp
              xmaxtemp = xtemp
           endif
        case('3')
           if (xmintemp .ne. 0) then
              xmintemp = 1./xmintemp
           else
              xmintemp = 0.
           endif
           if (xmaxtemp .ne. 0) then
              xmaxtemp = 1./xmaxtemp
           else
              xmaxtemp = 0.
           endif
        case('4')
           xmintemp = xmintemp**2
           xmaxtemp = xmaxtemp**2
        case('5')
           if (xmintemp.gt.0) then
              xmintemp = sqrt(xmintemp)
           else
              xmintemp = 0.
           endif
           if (xmaxtemp.gt.0) then
              xmaxtemp = sqrt(xmaxtemp)
           else
              xmaxtemp = 0.
           endif
        case('6')
           xmintemp = exp(xmintemp)
           xmaxtemp = exp(xmaxtemp)
        end select
     enddo

     xmin = min(xmintemp,xmaxtemp)
     xmax = max(xmintemp,xmaxtemp)

  endif

end subroutine transform_limits_inverse

!------------------------------------------------------------------------
!
!  function to adjust the label of a plot if log, 1/x etc
!
!  Note: *cannot* put print or write statements into this function
!        as it is used in the middle of write or print statements
!        this means that finding the digits is a bit trickier
!
!------------------------------------------------------------------------
function transform_label(label,itrans)
  implicit none
  integer, intent(in) :: itrans
  character(len=*), intent(in) :: label
  integer :: itransmulti,i,ndigits
  integer, dimension(5) :: digit
  character(len=len(label)+20) :: transform_label
  character(len=len(label)+20) :: temp_label
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then
     call get_digits(itrans,digit,ndigits)
     temp_label = label
     !
     !--do a transformation for each digit
     !
     do i=1,ndigits
        itransmulti = digit(i)
        !
        !--perform transformation appropriate to this digit
        !
        select case(itransmulti)
        case(1)
           temp_label = 'log '//trim(temp_label)
        case(2)
           temp_label = '|'//trim(temp_label)//'|'
        case(3)
           temp_label = '1/'//trim(temp_label)
        case(4)
           temp_label = 'sqrt('//trim(temp_label)//')'
        case(5)
           temp_label = trim(temp_label)//'\u2\d'
        case(6)
           temp_label = 'ln '//trim(temp_label)
        case default
           temp_label = trim(temp_label)
        end select
     enddo

     transform_label = temp_label
  else
     transform_label = label
  endif

end function transform_label

!------------------------------------------------------------------------
!     get_digits: for an integer i returns number of digits it contains
!     and a list of these *without* using write statements
!
!     i            : integer to split into digits
!     nmax           : dimensions of digits array
!     digits(nmax) : array of digits
!     ndigits      : number of digits in i
!------------------------------------------------------------------------

subroutine get_digits(i,digits,ndigits)
  implicit none
  integer, intent(in) :: i
  integer, intent(out) :: ndigits
  integer, intent(out), dimension(:) :: digits
  integer :: j,isubtract,idigit

  ndigits = 0

  isubtract = 0

  do j=size(digits),0,-1
     if (i.ge.10**j) then
        ndigits = ndigits + 1
        idigit = (i - isubtract)/10**j
        digits(ndigits) = idigit
        isubtract = isubtract + digits(ndigits)*10**j
     endif
  enddo

end subroutine get_digits

!------------------------------------------------------------------------
!     function returning the correction factor by which to multiply
!     to convert the log in the transform to a natural log
!------------------------------------------------------------------------
real function convert_to_ln_fac(itrans)
  implicit none
  integer, intent(in) :: itrans
  character(len=20) :: string
  integer :: i
  real :: xtemp

  !
  !--default conversion factor is unity
  !
  xtemp = 1.0
  !
  !--extract the digits from the input number
  !
  if (itrans.gt.0) then
     write(string,*) itrans
     do i=len_trim(string),1,-1  ! do digits in reverse
        !
        !--perform transformation appropriate to this digit
        !
        select case(string(i:i))
        !
        !--correction factor is ln(10.) but can be
        !  1./ln(10.), sqrt(1./ln(10.), etc. depending
        !  on other transformations
        !
        case('1')
           xtemp = log(10.)*xtemp
        case('2')
           xtemp = abs(xtemp)
        case('3')
           xtemp = 1./xtemp
        case('4')
           xtemp = sqrt(xtemp)
        case('5')
           xtemp = xtemp**2
        !case('6')
        !   xtemp = xtemp
        end select
     enddo
  endif

  convert_to_ln_fac = xtemp

end function convert_to_ln_fac

!---------------------------------------
! query function to return whether the
! transformation involves a log or not
!---------------------------------------
logical function islogged(itrans)
  implicit none
  integer, intent(in) :: itrans
  character(len=20) :: string

  write(string,"(i8)") itrans
  islogged = (index(string,'1').ne.0 .or. index(string,'6').ne.0)

end function islogged

end module transforms
