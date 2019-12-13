!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
module fparser
 !------- -------- --------- --------- --------- --------- --------- --------- -------
 ! Fortran 90 function parser v1.1
 !------- -------- --------- --------- --------- --------- --------- --------- -------
 !
 ! This function parser module is intended for applications where a set of mathematical
 ! fortran-style expressions is specified at runtime and is then evaluated for a large
 ! number of variable values. This is done by compiling the set of function strings
 ! into byte code, which is interpreted efficiently for the various variable values.
 !
 ! The source code is available from http://fparser.sourceforge.net
 !
 ! Please send comments, corrections or questions to the author:
 ! Roland Schmehl <roland.schmehl@alumni.uni-karlsruhe.de>
 !
 !------- -------- --------- --------- --------- --------- --------- --------- -------
 ! The function parser concept is based on a C++ class library written by  Juha
 ! Nieminen <warp@iki.fi> available from http://warp.povusers.org/FunctionParser/
 !------- -------- --------- --------- --------- --------- --------- --------- -------
 !
 ! Modifications by D, Price for integration in SPLASH:
 !  7th Aug 2009: added checkf interface routine to check syntax without compiling code
 !                added endf routine to stop memory leaks, also called from initf if needed
 !                bug fix with error message for sqrt(-ve)
 !
 !  9th Aug 2009: added Mathematical constant recognition (pi)
 !
 ! 27th Jan 2010: check for -ve numbers to fractional powers and zero to negative power added
 !
 ! 19th Oct 2016: added Fortran 2008 intrinsic functions
 !                added optional iErrType argument for error message printing
 !
 !  5th Mar 2019: added atan2 function (and parsing of functions with two arguments)

 IMPLICIT NONE
!--modification here by D.Price: define type parameters here rather than in a separate module
 integer, parameter, public  :: rn = KIND(0.0d0)          ! Precision of real numbers
!--modification (22/6/2018), adjusted integer kind to allow more than 128 variables
 integer, parameter, private :: is = SELECTED_INT_KIND(4) ! Data type of bytecode
!--end modification

 !------- -------- --------- --------- --------- --------- --------- --------- -------
 PUBLIC                     :: initf,    & ! Initialize function parser for n functions
                                parsef,   & ! Parse single function string
                                evalf,    & ! Evaluate single function
                                checkf,   & ! Check syntax in a function string
                                endf,     & ! Clean up memory once finished
                                EvalErrMsg  ! Error message (Use only when EvalErrType>0)
 integer, PUBLIC            :: EvalErrType ! =0: no error occured, >0: evaluation error
!--modification by D. Price: add parseErr parameter (used in checkf)
 integer, PRIVATE           :: ParseErrType ! =0: no error occured, >0: parse error
!--modification by D. Price: add verboseness internal variable (used in checkf)
 LOGICAL, PRIVATE           :: PrintErrors = .true. ! =0: no error occured, >0: parse error
!--modification by J. Wurster: added cgs and code versions of mu0
 real, PUBLIC               :: mu0=1.0_rn   ! =1 for code units, 4pi for cgs
 !------- -------- --------- --------- --------- --------- --------- --------- -------
 PRIVATE
 SAVE
 integer(is),                              PARAMETER :: cImmed   = 1,          &
                                                         cNeg     = 2,          &
                                                         cAdd     = 3,          &
                                                         cSub     = 4,          &
                                                         cMul     = 5,          &
                                                         cDiv     = 6,          &
                                                         cPow     = 7,          &
                                                         cAbs     = 8,          &
                                                         cExp     = 9,          &
                                                         cLog10   = 10,         &
                                                         cLog     = 11,         &
                                                         cSqrt    = 12,         &
                                                         cSinh    = 13,         &
                                                         cCosh    = 14,         &
                                                         cTanh    = 15,         &
                                                         cSin     = 16,         &
                                                         cCos     = 17,         &
                                                         cTan     = 18,         &
                                                         cAsin    = 19,         &
                                                         cAcos    = 20,         &
                                                         cAtan2   = 21,         &
                                                         cAtan    = 22,         &
                                                         cBesj0   = 23,         &
                                                         cBesj1   = 24,         &
                                                         cBesy0   = 25,         &
                                                         cBesy1   = 26,         &
                                                         cerfcs   = 27,         &
                                                         cerfc    = 28,         &
                                                         cerf     = 29,         &
                                                         cgamma   = 30,         &
                                                         VarBegin = 31
 CHARACTER (LEN=1), dimension(cAdd:cPow),  PARAMETER :: Ops      = (/ '+',     &
                                                                       '-',     &
                                                                       '*',     &
                                                                       '/',     &
                                                                       '^' /)
 CHARACTER (LEN=5), dimension(cAbs:cgamma), PARAMETER :: Funcs   = (/'abs  ', &
                                                                      'exp  ', &
                                                                      'log10', &
                                                                      'log  ', &
                                                                      'sqrt ', &
                                                                      'sinh ', &
                                                                      'cosh ', &
                                                                      'tanh ', &
                                                                      'sin  ', &
                                                                      'cos  ', &
                                                                      'tan  ', &
                                                                      'asin ', &
                                                                      'acos ', &
                                                                      'atan2', &
                                                                      'atan ', &
                                                                      'besj0', &
                                                                      'besj1', &
                                                                      'besy0', &
                                                                      'besy1', &
                                                                      'erfcs', &
                                                                      'erfc ', &
                                                                      'erf  ', &
                                                                      'gamf ' /)
 TYPE tComp
    integer(is), dimension(:), POINTER :: ByteCode
    integer                            :: ByteCodeSize
    real(rn),    dimension(:), POINTER :: Immed
    integer                            :: ImmedSize
    real(rn),    dimension(:), POINTER :: Stack
    integer                            :: StackSize, &
                                           StackPtr
 END TYPE tComp
 TYPE (tComp),  dimension(:),  POINTER :: Comp              ! Bytecode
 integer,   dimension(:),  ALLOCATABLE :: ipos              ! Associates function strings
 !
contains
 !
subroutine initf (n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Initialize function parser for n functions
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer, INTENT(in) :: n                                 ! Number of functions
 integer             :: i
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 if (ASSOCIATED(Comp)) then
    print "(a)",' fparser warning: initf called repeatedly without prior call to endf'
    CALL endf
 ENDIF
 ALLOCATE (Comp(n))
 DO i=1,n
    NULLIFY (Comp(i)%ByteCode,Comp(i)%Immed,Comp(i)%Stack)
 end do
end subroutine initf
 !
subroutine endf()
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Clean up memory at the end of the function parsing/evaluation calls (D. Price)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer             :: i
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 if (ASSOCIATED(Comp)) then
    DO i=1,size(Comp)
       if (ASSOCIATED(Comp(i)%ByteCode)) DEALLOCATE ( Comp(i)%ByteCode, &
                                                         Comp(i)%Immed,    &
                                                         Comp(i)%Stack     )
    enddo
    DEALLOCATE(Comp)
 ENDIF
end subroutine endf
 !
subroutine parsef (i, FuncStr, Var, err, Verbose)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Parse ith function string FuncStr and compile it into bytecode
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                         INTENT(in) :: i         ! Function identifier
 CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Function string
 CHARACTER (LEN=*), dimension(:), INTENT(in) :: Var       ! Array with variable names
 CHARACTER (LEN=LEN(FuncStr))                :: Func      ! Function string, local use
 integer, INTENT(OUT), OPTIONAL              :: err
 LOGICAL, INTENT(IN), OPTIONAL               :: Verbose   ! Turn error messages on/off
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 if (i < 1 .OR. i > SIZE(Comp)) then
    write(*,*) '*** Parser error: Function number ',i,' out of range'
    if (present(err)) err = 1
    return
 END IF
 EvalErrType  = 0   ! D. Price : to prevent accidental misuse
 ParseErrType = 0
 PrintErrors  = .true.
 if (present(Verbose)) PrintErrors = Verbose

 ALLOCATE (ipos(LEN(Func)))                       ! Char. positions in orig. string
 Func = FuncStr                                           ! Local copy of function string
 CALL Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
 CALL RemoveSpaces (Func)                                 ! Condense function string
 !CALL GetConstants (Func)
 CALL CheckSyntax (Func,FuncStr,Var)
 DEALLOCATE (ipos)
 if (present(err)) err = ParseErrType
 PrintErrors = .true.                                     ! reset this to true

 !--D. Price: return after ParseErr here instead of stop inside CheckSyntax
 if (ParseErrType /= 0) return
 CALL Compile (i,Func,Var)                                ! Compile into bytecode
end subroutine parsef
 !
integer FUNCTION checkf(FuncStr, Var, Verbose)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check syntax in a function string (added by D. Price) but do not compile it
 ! Returns an error code NOT related to ErrMsg
 ! Optional variable "verbose" determines whether or not error messages are printed
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),               INTENT(IN) :: FuncStr   ! Function string
 CHARACTER (LEN=*), dimension(:), INTENT(IN) :: Var       ! Array with variable names
 LOGICAL, INTENT(IN), OPTIONAL               :: Verbose   ! Turn error messages on/off
 CHARACTER (LEN=LEN(FuncStr))                :: Func      ! Function string, local use
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 EvalErrType  = 0   ! D. Price : to prevent accidental misuse
 ParseErrType = 0
 PrintErrors  = .true.
 if (present(Verbose)) PrintErrors = Verbose

 ALLOCATE (ipos(LEN(Func)))                               ! Char. positions in orig. string
 Func = FuncStr                                           ! Local copy of function string
 CALL Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
 CALL RemoveSpaces (Func)                                 ! Condense function string
 !CALL GetConstants (Func)
 CALL CheckSyntax (Func,FuncStr,Var)
 DEALLOCATE (ipos)

 PrintErrors = .true.                                     ! reset this to true
 checkf = ParseErrType
END FUNCTION checkf
 !
FUNCTION evalf (i, Val) RESULT (res)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Evaluate bytecode of ith function for the values passed in array Val(:)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                INTENT(in) :: i                  ! Function identifier
 real(rn), dimension(:), INTENT(in) :: Val                ! Variable values
 real(rn)                           :: res                ! Result
 integer                            :: IP,              & ! Instruction pointer
                                          DP,              & ! Data pointer
                                          SP                 ! Stack pointer
 real(rn),                PARAMETER :: zero = 0._rn
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 DP = 1
 SP = 0
 DO IP=1,Comp(i)%ByteCodeSize
    SELECT CASE (Comp(i)%ByteCode(IP))

    CASE (cImmed); SP=SP+1; Comp(i)%Stack(SP)=Comp(i)%Immed(DP); DP=DP+1
    CASE   (cNeg); Comp(i)%Stack(SP)=-Comp(i)%Stack(SP)
    CASE   (cAdd); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)+Comp(i)%Stack(SP); SP=SP-1
    CASE   (cSub); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)-Comp(i)%Stack(SP); SP=SP-1
    CASE   (cMul); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)*Comp(i)%Stack(SP); SP=SP-1
    CASE   (cDiv); if (Comp(i)%Stack(SP)==0._rn) then; EvalErrType=1; res=zero; return; ENDIF
       Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)/Comp(i)%Stack(SP); SP=SP-1
       ! D. Price: check for zero to negative powers and negative numbers to fractional powers
    CASE   (cPow); if (Comp(i)%Stack(SP-1)==0._rn .and.Comp(i)%Stack(SP)<0._rn) &
                      THEN; EvalErrType=1; res=zero; return; ENDIF
       if (Comp(i)%Stack(SP-1)<=0._rn .and.(Comp(i)%Stack(SP) /= nint(Comp(i)%Stack(SP)))) &
                      THEN; EvalErrType=5; res=zero; return; ENDIF
       Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)**Comp(i)%Stack(SP); SP=SP-1
    CASE   (cAbs); Comp(i)%Stack(SP)=ABS(Comp(i)%Stack(SP))
    CASE   (cExp); Comp(i)%Stack(SP)=exp(Comp(i)%Stack(SP))
    CASE (cLog10); if (Comp(i)%Stack(SP)<=0._rn) then; EvalErrType=3; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=LOG10(Comp(i)%Stack(SP))
    CASE   (cLog); if (Comp(i)%Stack(SP)<=0._rn) then; EvalErrType=3; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=log(Comp(i)%Stack(SP))
    CASE  (cSqrt); if (Comp(i)%Stack(SP)<0._rn) then; EvalErrType=2; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=sqrt(Comp(i)%Stack(SP))
    CASE  (cSinh); Comp(i)%Stack(SP)=SINH(Comp(i)%Stack(SP))
    CASE  (cCosh); Comp(i)%Stack(SP)=COSH(Comp(i)%Stack(SP))
    CASE  (cTanh); Comp(i)%Stack(SP)=TANH(Comp(i)%Stack(SP))
    CASE   (cSin); Comp(i)%Stack(SP)=sin(Comp(i)%Stack(SP))
    CASE   (cCos); Comp(i)%Stack(SP)=cos(Comp(i)%Stack(SP))
    CASE   (cTan); Comp(i)%Stack(SP)=TAN(Comp(i)%Stack(SP))
    CASE  (cAsin); if ((Comp(i)%Stack(SP)<-1._rn).OR.(Comp(i)%Stack(SP)>1._rn)) then
          EvalErrType=4; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=asin(Comp(i)%Stack(SP))
    CASE  (cAcos); if ((Comp(i)%Stack(SP)<-1._rn).OR.(Comp(i)%Stack(SP)>1._rn)) then
          EvalErrType=4; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=acos(Comp(i)%Stack(SP))
    CASE  (cAtan2); Comp(i)%Stack(SP-1)=atan2(Comp(i)%Stack(SP-1),Comp(i)%Stack(SP)); SP=SP-1
    CASE  (cAtan); Comp(i)%Stack(SP)=atan(Comp(i)%Stack(SP))
    CASE  (cBesj0); Comp(i)%Stack(SP)=bessel_j0(Comp(i)%Stack(SP))
    CASE  (cBesj1); Comp(i)%Stack(SP)=bessel_j1(Comp(i)%Stack(SP))
    CASE  (cBesy0); if (Comp(i)%Stack(SP)<=0._rn) then; EvalErrType=6; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=bessel_y0(Comp(i)%Stack(SP))
    CASE  (cBesy1); if (Comp(i)%Stack(SP)<=0._rn) then; EvalErrType=7; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=bessel_y1(Comp(i)%Stack(SP))
    CASE  (cerf);   Comp(i)%Stack(SP)=erf(Comp(i)%Stack(SP))
    CASE  (cerfc);  Comp(i)%Stack(SP)=erfc(Comp(i)%Stack(SP))
    CASE  (cerfcs); Comp(i)%Stack(SP)=erfc_scaled(Comp(i)%Stack(SP))
    CASE  (cgamma); if (Comp(i)%Stack(SP)==-abs(nint(Comp(i)%Stack(SP)))) then; EvalErrType=8; res=zero; return; ENDIF
       Comp(i)%Stack(SP)=gamma(Comp(i)%Stack(SP))
    CASE  DEFAULT; SP=SP+1; Comp(i)%Stack(SP)=Val(Comp(i)%ByteCode(IP)-VarBegin+1)
    END SELECT
 end do
 EvalErrType = 0
 res = Comp(i)%Stack(1)
END FUNCTION evalf
 !
subroutine CheckSyntax (Func,FuncStr,Var)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check syntax of function string,  returns 0 if syntax is ok
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),               INTENT(in) :: Func      ! Function string without spaces
 CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Original function string
 CHARACTER (LEN=*), dimension(:), INTENT(in) :: Var       ! Array with variable names
 integer(is)                                 :: n
 CHARACTER (LEN=1)                           :: c
 real(rn)                                    :: r
 LOGICAL                                     :: err
 integer                                     :: ParCnt, & ! Parenthesis counter
                                                   j,ib,in,lFunc,inold,ibold
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 j = 1
 ParCnt = 0
 lFunc = LEN_TRIM(Func)

 step: DO
    if (j > lFunc) then
       CALL ParseErrMsg (j, FuncStr)
       EXIT
    ENDIF
    c = Func(j:j)
    !-- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for valid operand (must appear)
    !-- -------- --------- --------- --------- --------- --------- --------- -------
    if (c == '-' .OR. c == '+') then                      ! Check for leading - or +
       j = j+1
       if (j > lFunc) then
          CALL ParseErrMsg (j, FuncStr, 'Missing operand')
          EXIT
       ENDIF
       c = Func(j:j)
       if (ANY(c == Ops) .or. c == ',') then
          CALL ParseErrMsg (j, FuncStr, 'Multiple operators')
          EXIT
       ENDIF
    END IF
    n = MathFunctionIndex (Func(j:))
    if (n > 0) then                                       ! Check for math function
       j = j+LEN_TRIM(Funcs(n))
       if (j > lFunc) then
          CALL ParseErrMsg (j, FuncStr, 'Missing function argument')
          EXIT
       ENDIF
       c = Func(j:j)
       if (c /= '(') then
          CALL ParseErrMsg (j, FuncStr, 'Missing opening parenthesis')
          EXIT
       ENDIF
       if (n == cAtan2) then                              ! Check #args for function with two arguments
          if (CountArgs(Func(j:)) /= 2) then
             CALL ParseErrMsg (j, FuncStr, 'Function expects two arguments')
             EXIT
          ENDIF
       ELSE
          if (CountArgs(Func(j:)) > 1) then
             CALL ParseErrMsg (j, FuncStr, 'Too many function arguments')
             EXIT
          ENDIF
       ENDIF
    END IF
    if (c == '(') then                                    ! Check for opening parenthesis
       ParCnt = ParCnt+1
       j = j+1
       CYCLE step
    END IF
    if (SCAN(c,'0123456789.') > 0) then                   ! Check for number
       r = RealNum (Func(j:),ib,in,err)
       if (err) then
          CALL ParseErrMsg (j, FuncStr, 'Invalid number format:  '//Func(j+ib-1:j+in-2))
          EXIT
       ENDIF
       j = j+in-1
       if (j > lFunc) EXIT
       c = Func(j:j)
    ELSE                                                  ! Check for variable
       n = VariableIndex (Func(j:),Var,ib,in)
       if (n == 0) then                                   ! DP: If not a variable, check for constants
          ibold = ib
          inold = in
          r = MathConst (Func(j:),ib,in,err)
          if (err) then                                   ! Return error if constants not found
             CALL ParseErrMsg (j, FuncStr, 'Invalid element: '//Func(j+ib-1:j+in-2))
             ib = ibold
             in = inold
             EXIT
          ENDIF
       ENDIF
       j = j+in-1
       if (j > lFunc) EXIT
       c = Func(j:j)
    END IF
    DO WHILE (c == ')')                                   ! Check for closing parenthesis
       ParCnt = ParCnt-1
       if (ParCnt < 0) CALL ParseErrMsg (j, FuncStr, 'Mismatched parenthesis')
       if (Func(j-1:j-1) == '(') CALL ParseErrMsg (j-1, FuncStr, 'Empty parentheses')
       j = j+1
       if (j > lFunc) EXIT
       c = Func(j:j)
    end do
    !-- -------- --------- --------- --------- --------- --------- --------- -------
    ! Now, we have a legal operand: A legal operator or end of string must follow
    !-- -------- --------- --------- --------- --------- --------- --------- -------
    if (j > lFunc) EXIT
    if (ANY(c == Ops) .or. c == ',') then                               ! Check for multiple operators
       if (j+1 > lFunc) CALL ParseErrMsg (j, FuncStr)
       if (ANY(Func(j+1:j+1) == Ops)) CALL ParseErrMsg (j+1, FuncStr, 'Multiple operators')
    ELSE                                                  ! Check for next operand
       CALL ParseErrMsg (j, FuncStr, 'Missing operator')
    END IF
    !-- -------- --------- --------- --------- --------- --------- --------- -------
    ! Now, we have an operand and an operator: the next loop will check for another
    ! operand (must appear)
    !-- -------- --------- --------- --------- --------- --------- --------- -------
    j = j+1
 end do step
 if (ParCnt > 0) CALL ParseErrMsg (j, FuncStr, 'Missing )')
end subroutine CheckSyntax
 !
FUNCTION EvalErrMsg ( ierrType ) RESULT (msg)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Return error message
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*), dimension(8), PARAMETER :: m = (/ 'Division by zero                   ', &
                                                         'Argument of SQRT negative          ', &
                                                         'Argument of LOG <= 0               ', &
                                                         'Argument of ASIN or ACOS illegal   ', &
                                                         '-ve number to fractional power     ', &
                                                         'Argument of Bessel_y0 <= 0         ', &
                                                         'Argument of Bessel_y1 <= 0         ', &
                                                         'Argument of Gamma function illegal '/)
 CHARACTER (LEN=LEN(m))                     :: msg
 integer,  INTENT(in), OPTIONAL :: ierrType
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 if (present(ierrType)) then
    if (iErrType < 1 .OR. iErrType > SIZE(m)) then
       msg = ''
    ELSE
       msg = m(iErrType)
    ENDIF
 ELSE
    if (EvalErrType < 1 .OR. EvalErrType > SIZE(m)) then
       msg = ''
    ELSE
       msg = m(EvalErrType)
    ENDIF
 ENDIF
END FUNCTION EvalErrMsg
 !
subroutine ParseErrMsg (j, FuncStr, Msg)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Print error message (modification by D.Price: do not terminate program,
 !                      also added option to not print error message)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                     INTENT(in) :: j
 CHARACTER (LEN=*),           INTENT(in) :: FuncStr       ! Original function string
 CHARACTER (LEN=*), OPTIONAL, INTENT(in) :: Msg
 integer                                 :: k
 !----- -------- --------- --------- --------- --------- --------- --------- -------

 if (PrintErrors) then
    if (PRESENT(Msg)) then
       write(*,*) '*** Error in syntax of function string: '//Msg
    ELSE
       write(*,*) '*** Error in syntax of function string:'
    ENDIF
    write(*,*)
    write(*,'(A)') ' '//FuncStr
    if (ALLOCATED(ipos)) then                               ! Avoid out-of-bounds-errors
       if (SIZE(ipos) >= j) then
          DO k=1,ipos(j)
             write(*,'(A)',ADVANCE='NO') ' '                ! Advance to the jth position
          end do
          write(*,'(A)') '?'
       ENDIF
    ENDIF
 ENDIF
 ParseErrType = 1

end subroutine ParseErrMsg
 !
FUNCTION OperatorIndex (c) RESULT (n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Return operator index
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=1), INTENT(in) :: c
 integer(is)                   :: n,j
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 n = 0
 DO j=cAdd,cPow
    if (c == Ops(j)) then
       n = j
       EXIT
    END IF
 end do
END FUNCTION OperatorIndex
 !
FUNCTION MathFunctionIndex (str) RESULT (n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Return index of math function beginnig at 1st position of string str
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*), INTENT(in) :: str
 integer(is)                   :: n,j
 integer                       :: k
 CHARACTER (LEN=LEN(Funcs))    :: fun
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 n = 0
 DO j=cAbs,cgamma                                         ! Check all math functions
    k = MIN(LEN_TRIM(Funcs(j)), LEN(str))
    CALL LowCase (str(1:k), fun)
    if (fun == Funcs(j)) then                             ! Compare lower case letters
       n = j                                              ! Found a matching function
       EXIT
    END IF
 end do
END FUNCTION MathFunctionIndex
 !
FUNCTION VariableIndex (str, Var, ibegin, inext) RESULT (n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Return index of variable at begin of string str (returns 0 if no variable found)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),               INTENT(in) :: str       ! String
 CHARACTER (LEN=*), dimension(:), INTENT(in) :: Var       ! Array with variable names
 integer(is)                                 :: n,j       ! Index of variable
 integer, OPTIONAL,              INTENT(out) :: ibegin, & ! Start position of variable name
                                                   inext     ! Position of character after name
 integer                                     :: ib,in,lstr
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 n = 0
 lstr = LEN_TRIM(str)
 if (lstr > 0) then
    DO ib=1,lstr                                          ! Search for first character in str
       if (str(ib:ib) /= ' ') EXIT                        ! When lstr>0 at least 1 char in str
    end do
    DO in=ib,lstr                                         ! Search for name terminators
       if (SCAN(str(in:in),'+-*/^,) ') > 0) EXIT
    end do
    DO j=1,SIZE(Var,kind=is)
       if (str(ib:in-1) == Var(j)) then
          n = j                                           ! Variable name found
          EXIT
       END IF
    end do
    !--else below added by D. Price - should never be required though
 ELSE                                                     ! blank string
    ib = 1  ! to avoid compiler warnings
    in = 2  ! and any possible seg fault
 END IF
 if (PRESENT(ibegin)) ibegin = ib
 if (PRESENT(inext))  inext  = in
END FUNCTION VariableIndex
 !
subroutine RemoveSpaces (str)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Remove Spaces from string, remember positions of characters in old string
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*), INTENT(inout) :: str
 integer                          :: k,lstr
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 lstr = LEN_TRIM(str)
 ipos = 0
 do k=1,lstr
    ipos(k) = k
 enddo
 k = 1
 DO WHILE (str(k:lstr) /= ' ')
    if (str(k:k) == ' ') then
       str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
       ipos(k:lstr) = (/ ipos(k+1:lstr), 0 /)             ! Move 1 element to left
       k = k-1
    END IF
    k = k+1
 end do
end subroutine RemoveSpaces
 !
subroutine Replace (ca,cb,str)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Replace ALL appearances of character set ca in string str by character set cb
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),       INTENT(in) :: ca
 CHARACTER (LEN=LEN(ca)), INTENT(in) :: cb                ! LEN(ca) must be LEN(cb)
 CHARACTER (LEN=*),    INTENT(inout) :: str
 integer                             :: j,lca
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 lca = LEN(ca)
 DO j=1,LEN_TRIM(str)-lca+1
    if (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
 end do
end subroutine Replace
 !
FUNCTION CountArgs(str) RESULT(n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Count number of arguments in a function string
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 CHARACTER (LEN=*), INTENT(in) :: str
 integer                       :: n,j
 n = 1
 DO j=1,len_trim(str)
    if (str(j:j) == ',') n = n + 1
 enddo
END FUNCTION CountArgs
 !
subroutine Compile (i, F, Var)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Compile i-th function string F into bytecode
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                         INTENT(in) :: i         ! Function identifier
 CHARACTER (LEN=*),               INTENT(in) :: F         ! Function string
 CHARACTER (LEN=*), dimension(:), INTENT(in) :: Var       ! Array with variable names
 integer                                     :: istat
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 if (ASSOCIATED(Comp(i)%ByteCode)) DEALLOCATE ( Comp(i)%ByteCode, &
                                                   Comp(i)%Immed,    &
                                                   Comp(i)%Stack     )
 Comp(i)%ByteCodeSize = 0
 Comp(i)%ImmedSize    = 0
 Comp(i)%StackSize    = 0
 Comp(i)%StackPtr     = 0
 CALL CompileSubstr (i,F,1,LEN_TRIM(F),Var)               ! Compile string to determine size
 ALLOCATE ( Comp(i)%ByteCode(Comp(i)%ByteCodeSize), &
               Comp(i)%Immed(Comp(i)%ImmedSize),       &
               Comp(i)%Stack(Comp(i)%StackSize),       &
               STAT = istat                            )
 if (istat /= 0) then
    write(*,*) '*** Parser error: Memory allocation for byte code failed'
    stop
 ELSE
    Comp(i)%ByteCodeSize = 0
    Comp(i)%ImmedSize    = 0
    Comp(i)%StackSize    = 0
    Comp(i)%StackPtr     = 0
    CALL CompileSubstr (i,F,1,LEN_TRIM(F),Var)            ! Compile string into bytecode
 END IF
 !
end subroutine Compile
 !
subroutine AddCompiledByte (i, b)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Add compiled byte to bytecode
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,     INTENT(in) :: i                             ! Function identifier
 integer(is), INTENT(in) :: b                             ! Value of byte to be added
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 Comp(i)%ByteCodeSize = Comp(i)%ByteCodeSize + 1
 if (ASSOCIATED(Comp(i)%ByteCode)) Comp(i)%ByteCode(Comp(i)%ByteCodeSize) = b
end subroutine AddCompiledByte
 !
FUNCTION MathItemIndex (i, F, Var) RESULT (n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Return math item index, if item is real number, enter it into Comp-structure
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                         INTENT(in) :: i         ! Function identifier
 CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
 CHARACTER (LEN=*), dimension(:), INTENT(in) :: Var       ! Array with variable names
 integer(is)                                 :: n         ! Byte value of math item
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 n = 0
 if (SCAN(F(1:1),'0123456789.') > 0) then                 ! Check for begin of a number
    Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
    if (ASSOCIATED(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = RealNum (F)
    n = cImmed
 ELSE                                                     ! Check for a variable
    n = VariableIndex (F, Var)
    if (n > 0) then
       n = VarBegin+n-1_is
    ELSE   ! Check for Mathematical constants
       n = MathConstIndex(i, F)
    ENDIF
 END IF
END FUNCTION MathItemIndex
 !
FUNCTION MathConstIndex (i, F, ibegin, inext) RESULT (n)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Routine added by D. Price
 ! Substitute values for Mathematical Constants (e.g. pi)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                         INTENT(in) :: i         ! Function identifier
 CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
 integer, OPTIONAL,              INTENT(out) :: ibegin, & ! Start position of real number
                                                   inext     ! 1st character after real number
 integer(is)                                 :: n         ! Byte value of math item
 real(rn)                                    :: res
 LOGICAL                                     :: err
 integer                                     :: ib,in
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 n = 0
 res = MathConst(F,ib,in,err)
 if (.not.err) then
    Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
    if (ASSOCIATED(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = res
    n = cImmed
 ELSE
    ib = 1
    in = 1
 END IF
 if (PRESENT(ibegin)) ibegin = ib
 if (PRESENT(inext))  inext  = in

END FUNCTION MathConstIndex
 !
FUNCTION CompletelyEnclosed (F, b, e) RESULT (res)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*), INTENT(in) :: F                       ! Function substring
 integer,           INTENT(in) :: b,e                     ! First and last pos. of substring
 LOGICAL                       :: res
 integer                       :: j,k
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 res=.false.
 if (F(b:b) == '(' .AND. F(e:e) == ')') then
    k = 0
    DO j=b+1,e-1
       IF     (F(j:j) == '(') then
          k = k+1
       ELSEif (F(j:j) == ')') then
          k = k-1
       END IF
       if (k < 0) EXIT
    end do
    if (k == 0) res=.true.                                ! All opened parenthesis closed
 END IF
END FUNCTION CompletelyEnclosed

FUNCTION TwoArgs (F, b, e, m) RESULT (res)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check if function substring F(b:e) has two arguments i.e. f(a,b)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*), INTENT(in) :: F                       ! Function substring
 integer,           INTENT(in) :: b,e                     ! First and last pos. of substring
 integer,           INTENT(out) :: m
 LOGICAL                       :: res

 res = .false.
 m = INDEX(F(b:e),',')
 if (m > 0) then
    res = .true.
    m = b + m - 1
 ENDIF

END FUNCTION TwoArgs
 !
RECURSIVE subroutine CompileSubstr (i, F, b, e, Var)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Compile i-th function string F into bytecode
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,                         INTENT(in) :: i         ! Function identifier
 CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
 integer,                         INTENT(in) :: b,e       ! Begin and end position substring
 CHARACTER (LEN=*), dimension(:), INTENT(in) :: Var       ! Array with variable names
 integer(is)                                 :: n
 integer                                     :: b2,j,k,io,m
 CHARACTER (LEN=*),                PARAMETER :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                                            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check for special cases of substring
 !----- -------- --------- --------- --------- --------- --------- --------- -------
!      write(*,*) 'PARSING F(b:e) = '//F(b:e)
 IF     (F(b:b) == '+') then                              ! Case 1: F(b:e) = '+...'
!      write(*,*)'1. F(b:e) = "+..."'
    CALL CompileSubstr (i, F, b+1, e, Var)
    return
 ELSEif (CompletelyEnclosed (F, b, e)) then               ! Case 2: F(b:e) = '(...)'
!      write(*,*)'2. F(b:e) = "(...)"'
    CALL CompileSubstr (i, F, b+1, e-1, Var)
    return
 ELSEif (SCAN(F(b:b),calpha) > 0) then
    n = MathFunctionIndex (F(b:e))
    if (n > 0) then
       b2 = b+INDEX(F(b:e),'(')-1
       if (CompletelyEnclosed(F, b2, e)) then             ! Case 3: F(b:e) = 'fcn(...)'
!            write(*,*)'3. F(b:e) = "fcn(...)"'
          if (n == cAtan2 .and. TwoArgs(F,b2+1,e-1,m)) then
!                print*,' SPLITTING ',b,m,e,' F(b:e)=',F(b2+1:e-1)
             CALL CompileSubstr (i, F, b2+1, m-1, Var)
             CALL CompileSubstr (i, F, m+1, e-1, Var)
          ELSE
             CALL CompileSubstr(i, F, b2+1, e-1, Var)
          ENDIF
          CALL AddCompiledByte (i, n)
          return
       END IF
    END IF
 ELSEif (F(b:b) == '-') then
    if (CompletelyEnclosed (F, b+1, e)) then              ! Case 4: F(b:e) = '-(...)'
!         write(*,*)'4. F(b:e) = "-(...)"'
       CALL CompileSubstr (i, F, b+2, e-1, Var)
       CALL AddCompiledByte (i, cNeg)
       return
    ELSEif (SCAN(F(b+1:b+1),calpha) > 0) then
       n = MathFunctionIndex (F(b+1:e))
       if (n > 0) then
          b2 = b+INDEX(F(b+1:e),'(')
          if (CompletelyEnclosed(F, b2, e)) then          ! Case 5: F(b:e) = '-fcn(...)'
!               write(*,*)'5. F(b:e) = "-fcn(...)"'
             if (n == cAtan2 .and. TwoArgs(F,b2+1,e-1,m)) then
                CALL CompileSubstr (i, F, b2+1, m-1, Var)
                CALL CompileSubstr (i, F, m+1, e-1, Var)
             ELSE
                CALL CompileSubstr(i, F, b2+1, e-1, Var)
             ENDIF
             CALL AddCompiledByte (i, n)
             CALL AddCompiledByte (i, cNeg)
             return
          END IF
       END IF
    ENDIF
 ENDIF
!   write(*,*) 'PROCEED TO OPERATORS',TwoArgs(F,b,e,m),F(b:e)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 DO io=cAdd,cPow                                          ! Increasing priority +-*/^
    k = 0
    DO j=e,b,-1
       IF     (F(j:j) == ')') then
          k = k+1
       ELSEif (F(j:j) == '(') then
          k = k-1
       END IF
       if (k == 0 .AND. F(j:j) == Ops(io) .AND. IsBinaryOp (j, F)) then
          if (ANY(F(j:j) == Ops(cMul:cPow)) .AND. F(b:b) == '-') then ! Case 6: F(b:e) = '-...Op...' with Op > -
!               write(*,*)'6. F(b:e) = "-...Op..." with Op > -'
             CALL CompileSubstr (i, F, b+1, e, Var)
             CALL AddCompiledByte (i, cNeg)
             return
          ELSE                                                        ! Case 7: F(b:e) = '...BinOp...'
!               write(*,*)'7. Binary operator',F(j:j)
             CALL CompileSubstr (i, F, b, j-1, Var)
             CALL CompileSubstr (i, F, j+1, e, Var)
             CALL AddCompiledByte (i, OperatorIndex(Ops(io)))
             Comp(i)%StackPtr = Comp(i)%StackPtr - 1
             return
          END IF
       END IF
    end do
 end do
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check for remaining items, i.e. variables or explicit numbers
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 b2 = b
 if (F(b:b) == '-') b2 = b2+1
 n = MathItemIndex(i, F(b2:e), Var)
!   write(*,*)'8. AddCompiledByte ',n
 CALL AddCompiledByte (i, n)
 Comp(i)%StackPtr = Comp(i)%StackPtr + 1
 if (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
 if (b2 > b) CALL AddCompiledByte (i, cNeg)
end subroutine CompileSubstr
 !
FUNCTION IsBinaryOp (j, F) RESULT (res)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Check if operator F(j:j) in string F is binary operator
 ! Special cases already covered elsewhere:              (that is corrected in v1.1)
 ! - operator character F(j:j) is first character of string (j=1)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 integer,           INTENT(in) :: j                       ! Position of Operator
 CHARACTER (LEN=*), INTENT(in) :: F                       ! String
 LOGICAL                       :: res                     ! Result
 integer                       :: k
 LOGICAL                       :: Dflag,Pflag
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 res=.true.
 if (F(j:j) == '+' .OR. F(j:j) == '-') then               ! Plus or minus sign:
    if (j == 1) then                                      ! - leading unary operator ?
       res = .false.
    ELSEif (SCAN(F(j-1:j-1),'+-*/^(') > 0) then           ! - other unary operator ?
       res = .false.
    ELSEif (SCAN(F(j+1:j+1),'0123456789') > 0 .AND. &     ! - in exponent of real number ?
               SCAN(F(j-1:j-1),'eEdD')       > 0) then
       Dflag=.false.; Pflag=.false.
       k = j-1
       DO WHILE (k > 1)                                   !   step to the left in mantissa
          k = k-1
          IF     (SCAN(F(k:k),'0123456789') > 0) then
             Dflag=.true.
          ELSEif (F(k:k) == '.') then
             if (Pflag) then
                EXIT                                      !   * EXIT: 2nd appearance of '.'
             ELSE
                Pflag=.true.                              !   * mark 1st appearance of '.'
             ENDIF
          ELSE
             EXIT                                         !   * all other characters
          END IF
       end do
       if (Dflag .AND. (k == 1 .OR. SCAN(F(k:k),'+-*/^(') > 0)) res = .false.
    END IF
 END IF
END FUNCTION IsBinaryOp
 !
FUNCTION RealNum (str, ibegin, inext, error) RESULT (res)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),  INTENT(in) :: str                    ! String
 real(rn)                       :: res                    ! Real number
 integer, OPTIONAL, INTENT(out) :: ibegin,              & ! Start position of real number
                                      inext                  ! 1st character after real number
 LOGICAL, OPTIONAL, INTENT(out) :: error                  ! Error flag
 integer                        :: ib,in,istat
 LOGICAL                        :: Bflag,               & ! .T. at begin of number in str
                                      InMan,               & ! .T. in mantissa of number
                                      Pflag,               & ! .T. after 1st '.' encountered
                                      Eflag,               & ! .T. at exponent identifier 'eEdD'
                                      InExp,               & ! .T. in exponent of number
                                      DInMan,              & ! .T. if at least 1 digit in mant.
                                      DInExp,              & ! .T. if at least 1 digit in exp.
                                      err                    ! Local error flag
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 Bflag=.true.; InMan=.false.; Pflag=.false.; Eflag=.false.; InExp=.false.
 DInMan=.false.; DInExp=.false.
 ib   = 1
 in   = 1
 DO WHILE (in <= LEN_TRIM(str))
    SELECT CASE (str(in:in))
    CASE (' ')                                            ! Only leading blanks permitted
       ib = ib+1
       if (InMan .OR. Eflag .OR. InExp) EXIT
    CASE ('+','-')                                        ! Permitted only
       IF     (Bflag) then
          InMan=.true.; Bflag=.false.                     ! - at beginning of mantissa
       ELSEif (Eflag) then
          InExp=.true.; Eflag=.false.                     ! - at beginning of exponent
       ELSE
          EXIT                                            ! - otherwise stop
       ENDIF
    CASE ('0':'9')                                        ! Mark
       IF     (Bflag) then
          InMan=.true.; Bflag=.false.                     ! - beginning of mantissa
       ELSEif (Eflag) then
          InExp=.true.; Eflag=.false.                     ! - beginning of exponent
       ENDIF
       if (InMan) DInMan=.true.                           ! Mantissa contains digit
       if (InExp) DInExp=.true.                           ! Exponent contains digit
    CASE ('.')
       IF     (Bflag) then
          Pflag=.true.                                    ! - mark 1st appearance of '.'
          InMan=.true.; Bflag=.false.                     !   mark beginning of mantissa
       ELSEif (InMan .AND..NOT.Pflag) then
          Pflag=.true.                                    ! - mark 1st appearance of '.'
       ELSE
          EXIT                                            ! - otherwise stop
       END IF
    CASE ('e','E','d','D')                                ! Permitted only
       if (InMan) then
          Eflag=.true.; InMan=.false.                     ! - following mantissa
       ELSE
          EXIT                                            ! - otherwise stop
       ENDIF
    CASE DEFAULT
       EXIT                                               ! stop at all other characters
    END SELECT
    in = in+1
 end do
 err = (ib > in-1) .OR. (.NOT.DInMan) .OR. ((Eflag.OR.InExp).AND..NOT.DInExp)
 if (err) then
    res = 0.0_rn
 ELSE
    read(str(ib:in-1),*,IOSTAT=istat) res
    err = istat /= 0
 END IF
 if (PRESENT(ibegin)) ibegin = ib
 if (PRESENT(inext))  inext  = in
 if (PRESENT(error))  error  = err
END FUNCTION RealNum

 !
FUNCTION MathConst (str, ibegin, inext, error) RESULT (res)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Return values of Mathematical constants in string
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),  INTENT(in) :: str                    ! String
 real(rn)                       :: res                    ! Real number
 integer, OPTIONAL, INTENT(out) :: ibegin,              & ! Start position of real number
                                      inext                  ! 1st character after real number
 LOGICAL, OPTIONAL, INTENT(out) :: error                  ! Error flag
 integer                        :: ib,in
 LOGICAL                        :: err                    ! Local error flag
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ib = 1
 in = 1
 err = .false.
 if (len(str) < 2) then
    res = 0.0_rn
    err = .true.
 ELSE
    if (str(1:2)=='pi') then
       res = 3.14159265358979323846_rn
       in = 3
    ELSEif (str(1:2)=='mu') then
       res = mu0
       in = 3
    ELSE
       res = 0.0_rn
       err = .true.
    ENDIF
 ENDIF

 if (PRESENT(ibegin)) ibegin = ib
 if (PRESENT(inext))  inext  = in
 if (PRESENT(error))  error  = err

END FUNCTION MathConst

 !
subroutine LowCase (str1, str2)
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 ! Transform upper case letters in str1 into lower case letters, result is str2
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 IMPLICIT NONE
 CHARACTER (LEN=*),  INTENT(in) :: str1
 CHARACTER (LEN=*), INTENT(out) :: str2
 integer                        :: j,k
 CHARACTER (LEN=*),   PARAMETER :: lc = 'abcdefghijklmnopqrstuvwxyz'
 CHARACTER (LEN=*),   PARAMETER :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 !----- -------- --------- --------- --------- --------- --------- --------- -------
 str2 = str1
 DO j=1,LEN_TRIM(str1)
    k = INDEX(uc,str1(j:j))
    if (k > 0) str2(j:j) = lc(k:k)
 end do
end subroutine LowCase
 !
end module fparser
