!This is a small f90 module containing a generic subroutine for prompting
!integer, real, double and logical variables and strings.
!I use them  quite often and found the solution below very useful.
!
!The general syntax is:
!
!   call prompt(PROMPT,VAR,...)
!
!Action: writes the string PROMPT on the terminal plus the current
!        value of the variable VAR and reads VAR. If  <cr>
!        is pressed instead of a new value,
!        the variable VAR stays untouched.
!
!In addition there are a few optional parameter to the routine like
!setting defaults or limits etc... 
!         
!In principle three f90 features are used which are not available
!in f77:  recursion, non-advancing I/O and overloading.
!
!A detailed description of the syntax can be found in the header
!of the module below.

!---------------------------- please cut here -------------------------------------

!
! f90 Module 'prompting'                                     
!
! Definition of Generic Subroutine: prompt
!
! Syntax:    prompt(text, value, [min], [max]) 
!                text      character string
!                value     integer, real or double 
!                min, max  allowed range of same type as value (optional)
!
!            prompt(text, string, [length], [case])
!                text      character string
!                string    character string 
!                length    length of string (optional on return)
!                case      option 
!                          1 -> convert string to lower case
!                          2 -> convert string to upper case
!                          lower=1, upper=2 are defined public
!                          within this module
!
!            prompt(text, value, [default])
!                text      character string
!                value     logical 
!                default   logical (optional) 
!                          will always overwrite the current value
!
! Author: Th. S. Ullrich, University Heidelberg
! e-mail: ullrich@ceres1.physi.uni-heidelberg.de
! Last mod: 18 Aug 94
!

module prompting

   private                     

   !
   ! Options for string prompting routine
   !

   integer, parameter, public :: lower = 1, upper = 2
   
   ! 
   !  Create generic name 'prompt' 
   !
   
   interface prompt
      module procedure &
      integer_prompt, real_prompt, string_prompt, double_prompt, logical_prompt
   end interface
   public :: prompt 
    
contains

   ! 
   !  Integer prompting routine 
   !
   
   recursive subroutine integer_prompt(text, value, min, max)
      character(len=*), intent(in) :: text
      integer                      :: value, newvalue
      character(len=64)            :: string
      character(len=16)            :: chmin, chmax
      integer                      :: ios
      integer, optional            :: min, max
      logical                      :: error
      
      
      chmin = ''                        
      chmax = ''
      error = .false.
      
      !
      !  Pack arguments in strings for compact and nicer prompt
      !
      
      write(string,*) value                
      if (present(min)) write(chmin,"(g10.0)") min
      if (present(max)) write(chmax,"(g10.0)") max     
      
      !
      !  Write prompt string to terminal
      !
      
      if (present(min).or.present(max)) then
         write(*,"(a,1x,'([',a,':',a,'],',1x,'<cr>=',a,'):',1x)",advance='no') &
                 trim(adjustl(text)), trim(adjustl(chmin)), &
                 trim(adjustl(chmax)), trim(adjustl(string))
      else
         write(*,"(a,1x,'(<cr>=',a,'):',1x)",advance='no') &
                 trim(adjustl(text)), trim(adjustl(string))
      endif
      
      !
      !  Read new value, quit and keep old value if zero sized string
      !
      
      read(*,"(a)") string
      if (len(trim(adjustl(string))) == 0) return      
      read(string,"(g10.0)",iostat=ios) newvalue

      !
      !  Check if new string is of right type and within given range
      !
      
      if (ios /= 0) then
         print "(a)", "Error, not an integer number"
         error = .true.
      else
         if (present(min) .and. newvalue < min) then
            print "(a)", "Error, value out of range"
            error = .true.
         endif
         if (present(max) .and. newvalue > max) then
            print "(a)", "Error, value out of range"
            error = .true.
         endif
      endif
       
      !
      !  Assign new value if everything is ok, else prompt again
      !
      
      if (error) then
         call integer_prompt(text, value, min, max)
      else         
         value = newvalue
      endif
            
   end subroutine integer_prompt
   
   
   ! 
   !  Real prompting routine 
   !
   
   recursive subroutine real_prompt(text, value, min, max)
      character(len=*), intent(in) :: text
      real                         :: value, newvalue
      character(len=64)            :: string
      character(len=16)            :: chmin, chmax
      integer                      :: ios
      real, optional               :: min, max
      logical                      :: error
      
      
      chmin = ''                        
      chmax = ''
      error = .false.
      
      !
      !  Pack arguments in strings for compact and nicer prompt
      !
      
      write(string,"(g13.4)") value                
      if (present(min)) write(chmin,"(g13.4)") min
      if (present(max)) write(chmax,"(g13.4)") max     
      
      !
      !  Write prompt string to terminal
      !
      
      if (present(min).or.present(max)) then
         write(*,"(a,1x,'([',a,':',a,'],',1x,'<cr>=',a,'):',1x)",advance='no') &
                 trim(adjustl(text)), trim(adjustl(chmin)), &
                 trim(adjustl(chmax)), trim(adjustl(string))
      else
         write(*,"(a,1x,'(<cr>=',a,'):',1x)",advance='no') &
                 trim(adjustl(text)), trim(adjustl(string))
      endif
      
      !
      !  Read new value, quit and keep old value if zero sized string
      !
      
      read(*,"(a)") string
      if (len(trim(adjustl(string))) == 0) return      
      read(string,*,iostat=ios) newvalue

      !
      !  Check if new string is of right type and within given range
      !
      
      if (ios /= 0) then
         print "(a)", "Error, not a real number"
         error = .true.
      else
         if (present(min) .and. newvalue < min) then
            print "(a)", "Error, value out of range"
            error = .true.
         endif
         if (present(max) .and. newvalue > max) then
            print "(a)", "Error, value out of range"
            error = .true.
         endif
      endif
       
      !
      !  Assign new value if everything is ok, else prompt again
      !
      
      if (error) then
         call real_prompt(text, value, min, max)
      else         
         value = newvalue
      endif
            
   end subroutine real_prompt
         
   
   ! 
   !  Double precision prompting routine 
   !
   
   recursive subroutine double_prompt(text, value, min, max)
      integer, parameter           :: db = kind(0.d0)
      character(len=*), intent(in) :: text
      real(kind=db)                :: value, newvalue
      character(len=64)            :: string
      character(len=16)            :: chmin, chmax
      integer                      :: ios
      real(kind=db), optional      :: min, max
      logical                      :: error
      
      
      chmin = ''                        
      chmax = ''
      error = .false.
      
      !
      !  Pack arguments in strings for compact and nicer prompt
      !
      
      write(string,"(g13.4)") value                
      if (present(min)) write(chmin,"(g13.4)") min
      if (present(max)) write(chmax,"(g13.4)") max     
      
      !
      !  Write prompt string to terminal
      !
      
      if (present(min).or.present(max)) then
         write(*,"(a,1x,'([',a,':',a,'],',1x,'<cr>=',a,'):',1x)",advance='no') &
                 trim(adjustl(text)), trim(adjustl(chmin)), &
                 trim(adjustl(chmax)), trim(adjustl(string))
      else
         write(*,"(a,1x,'(<cr>=',a,'):',1x)",advance='no') &
                 trim(adjustl(text)), trim(adjustl(string))
      endif
      
      !
      !  Read new value, quit and keep old value if zero sized string
      !
      
      read(*,"(a)") string      
      if (len(trim(adjustl(string))) == 0) return
      read(string,*,iostat=ios) newvalue

      !
      !  Check if new string is of right type and within given range
      !
      
      if (ios /= 0) then
         print "(a)", "Error, not a real number"
         error = .true.
      else
         if (present(min) .and. newvalue < min) then
            print "(a)", "Error, value out of range"
            error = .true.
         endif
         if (present(max) .and. newvalue > max) then
            print "(a)", "Error, value out of range"
            error = .true.
         endif
      endif
       
      !
      !  Assign new value if everything is ok, else prompt again
      !
      
      if (error) then
         call double_prompt(text, value, min, max)
      else         
         value = newvalue
      endif
                  
   end subroutine double_prompt

   
   
   ! 
   !  Logical prompting routine 
   !
   
   recursive subroutine logical_prompt(text, lvalue, default)
      character(len=*), intent(in)  :: text
      logical                       :: lvalue
      logical, optional, intent(in) :: default
      character(len=32)             :: string
      
      !
      !  If present, set default 
      !
      
      if (present(default)) lvalue = default 
      
      !
      !  Default answer yes/no
      !
      
      if (lvalue) then
         string='yes'
      else
         string='no'
      endif
      
      
      !
      !  Write prompt string to terminal
      !
      
      write(*,"(a,1x,'(<cr>=',a,'):',1x)",advance='no') &
           trim(adjustl(text)), trim(adjustl(string))
           
      !
      !  Read new value, quit and keep old value if zero sized string
      !
      
      read(*,"(a)") string     
      if (len(trim(adjustl(string))) == 0) return
      
      !
      !  Translate answer in .true./.false., if invalid prompt again
      !
      
      select case (adjustl(string))
      case ('y')
           lvalue = .true.
      case ('yes')
           lvalue = .true.
      case ('on')
           lvalue = .true.
      case ('t')
           lvalue = .true.
      case ('true')
           lvalue = .true.
      case ('n')
           lvalue = .false.
      case ('no')
           lvalue = .false.
      case ('off')
           lvalue = .false.
      case ('f')
           lvalue = .false.
      case ('false')
           lvalue = .false.
      case default
         print "(a)", "Error, answer y(es)/t(rue)/on or n(o)/f(alse)/off"
         call logical_prompt(text, lvalue, default)
      end select
      
   end subroutine logical_prompt
      
   
   ! 
   !  String prompting routine 
   !
   
   subroutine string_prompt(text, string, length, case)
      character(len=*), intent(in)   :: text
      character(len=*)               :: string
      character(len=128)             :: newstring
      integer, optional, intent(out) :: length
      integer, optional              :: case
      integer                        :: is, ia
      integer, parameter             :: aoffset = 32
      
      
      !
      !  Write prompt string to terminal
      !
      
      write(*,"(a,1x,'(<cr>=',a,'):',1x)",advance='no') &
           trim(adjustl(text)), trim(adjustl(string))
           
      !
      !  Read new value, quit and keep old value if zero sized string
      !
      
      read(*,"(a)") newstring      
      if ( len(trim(adjustl(newstring))) /= 0 ) string = newstring
      if (present(length)) length = len(trim(string))
         
      !
      !  Convert string to upper/lower case if requested
      !
      
      if (present(case)) then
         if (case == upper) then
            do is = 1, len(string)
               ia = iachar(string(is:is))
               if (ia >= iachar('a').and.ia <= iachar('z')) &
                   string(is:is) = achar(ia-aoffset)
            enddo
         endif
         if (case == lower) then
            do is = 1, len(string)
               ia = iachar(string(is:is))
               if (ia >= iachar('A').and.ia <= iachar('Z')) &
                   string(is:is) = achar(ia+aoffset)
            enddo
         endif
      endif
      
   end subroutine string_prompt
   
end module prompting


