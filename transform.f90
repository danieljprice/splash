!------------------------------------------------------------------------
!
! module containing global parameter for transform subroutines
! this is the maximum number of transformations you can do in a row      
!
!------------------------------------------------------------------------
module transform_params
  integer, parameter :: nmax = 10
end module transform_params

!------------------------------------------------------------------------
!
!  subroutine returns log, 1/x of a given array
!
!  * can specify up to 9 individual operations to perform
!  * combinations of transformations are done when the 
!    input number is > 10 (e.g. 321 means 1/x, then abs, then log10)
!
!------------------------------------------------------------------------
subroutine transform(array,arrayout,itrans,isize)
  use transform_params
  implicit none
  integer :: i,ndigits,itransmulti
  integer, intent(IN) :: isize,itrans
  real, dimension(isize), intent(IN) :: array
  real, dimension(isize), intent(OUT) :: arrayout
  real, dimension(isize) :: arraytemp
  integer, dimension(nmax) :: digit     
  !
  !--extract the digits from the input number 
  !            
  if (itrans.gt.0) then      
     call get_digits(itrans,digit,ndigits,nmax)         
     !
     !--do a transformation for each digit     
     !
     arraytemp = array

     do i=1,ndigits
        itransmulti = digit(i)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(itransmulti)
        case(1)
           where (arraytemp > 0)
              arraytemp = log10(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case(2)
           arraytemp = abs(arraytemp)    
        case(3)
           where (arraytemp .ne. 0)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = 0.
           end where
        case(4) 
           where (arraytemp .gt. 0)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case(5)
           arraytemp = arraytemp**2	     
        end select
     enddo

     arrayout = arraytemp

  else
     arrayout = array
  endif

end subroutine transform

!------------------------------------------------------------------------
!
!  same as transform but for the plot limits
!  (min can become max and vice versa)
!
!------------------------------------------------------------------------
subroutine transform_limits(xminin,xmaxin,xminout,xmaxout,itrans)
  use transform_params
  implicit none
  integer :: i,ndigits,itransmulti
  integer, intent(in) :: itrans
  real, intent(in) :: xminin,xmaxin
  real, intent(out) :: xminout,xmaxout
  real :: xmintemp,xmaxtemp
  integer, dimension(nmax) :: digit     
  !
  !--extract the digits from the input number 
  !            
  if (itrans.gt.0) then      
     call get_digits(itrans,digit,ndigits,nmax)         
     !
     !--do a transformation for each digit     
     !
     xmintemp = xminin
     xmaxtemp = xmaxin

     do i=1,ndigits
        itransmulti = digit(i)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(itransmulti)
        case(1)
           if (xmintemp > 0) then
              xmintemp = log10(xmintemp)
           elseif (xmintemp.eq.0) then
              print*,' log10(xmin = 0): min set to 10-5'
              xmintemp = -5.
           endif
           if (xmaxtemp > 0) then
              xmaxtemp = log10(xmaxtemp)
           elseif (xmaxtemp.eq.0) then
              print*,' log10(xmax = 0): max set to 10-5'
              xmaxtemp = -5.
           endif
        case(2)
           xmintemp = abs(xmintemp)
           xmaxtemp = abs(xmaxtemp)
        case(3)
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
        case(4) 
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
        case(5)
           xmintemp = xmintemp**2
           xmaxtemp = xmaxtemp**2
        end select
     enddo

     xminout = min(xmintemp,xmaxtemp)
     xmaxout = max(xmintemp,xmaxtemp)

  else
     xminout = xminin
     xmaxout = xmaxin
  endif

end subroutine transform_limits

!------------------------------------------------------------------------
!
!  same as transform but for a two dimensional array
!  applies the transformation to the same array as was input
!
!------------------------------------------------------------------------
subroutine transform2(array,itrans,isizex,isizey)
  use transform_params
  implicit none
  integer :: i,ndigits,itransmulti
  integer, intent(IN) :: itrans,isizex,isizey
  real, dimension(isizex,isizey), intent(INOUT) :: array
!!  real, dimension(isizex,isizey), intent(OUT) :: arrayout
  real, dimension(isizex,isizey) :: arraytemp
  integer, dimension(nmax) :: digit     
  !
  !--extract the digits from the input number 
  !            
  if (itrans.gt.0) then      
     call get_digits(itrans,digit,ndigits,nmax)         
     !
     !--do a transformation for each digit     
     !
     arraytemp = array

     do i=1,ndigits
        itransmulti = digit(i)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(itransmulti)
        case(1)
           where (arraytemp > 0)
              arraytemp = log10(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case(2)
           arraytemp = abs(arraytemp)    
        case(3)
           where (arraytemp .ne. 0)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = 0.
           end where
        case(4) 
           where (arraytemp .gt. 0)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case(5)
           arraytemp = arraytemp**2	 
        end select
     enddo

     array = arraytemp

  else
!     do i = 1,isizex
!        do j = 1,isizey
!           print*,i,j
!           print*,array(i,j)
!        enddo
!     enddo
!     arrayout = array    
  endif

end subroutine transform2


!------------------------------------------------------------------------
!
!  function to adjust the label of a plot if log, 1/x etc
!
!  Note: *cannot* put print or write statements into this function
!        as it is used in the middle of write or print statements
!      
!------------------------------------------------------------------------
function transform_label(label,itrans)
  use transform_params
  implicit none
  integer :: itrans,itransmulti,i,itransprev,ndigits
  integer, dimension(nmax) :: digit
  character(LEN=*) :: label, transform_label
  character(LEN=120) :: temp_label      
  !
  !--extract the digits from the input number 
  !            
  if (itrans.gt.0) then      
     call get_digits(itrans,digit,ndigits,nmax)         
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
           temp_label = 'log\d10\u'//trim(temp_label)
        case(2)
           temp_label = '|'//trim(temp_label)//'|'
        case(3)
           temp_label = '1/'//trim(temp_label)
        case(4)
           temp_label = 'SQRT('//trim(temp_label)//')'
        case(5)
           temp_label = trim(temp_label)//'\u2\d'
        case DEFAULT
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
!     and a list of these
!
!     i            : integer to split into digits
!     nmax	   : dimensions of digits array
!     digits(nmax) : array of digits
!     ndigits      : number of digits in i
!------------------------------------------------------------------------

subroutine get_digits(i,digits,ndigits,nmax)
  implicit none
  integer, intent(IN) :: i,nmax
  integer, intent(OUT) :: ndigits
  integer, intent(OUT), dimension(nmax) :: digits
  integer :: j,isubtract,idigit

  ndigits = 0
  isubtract = 0      

  do j=nmax,0,-1
     if (i.ge.10**j) then 
        ndigits = ndigits + 1
        idigit = (i - isubtract)/10**j
        digits(ndigits) = idigit
        isubtract = isubtract + digits(ndigits)*10**j
     endif
  enddo

end subroutine get_digits
