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
! Format Krome/Phantom chemistry species names for Giza plot labels
! (TeX subscripts/superscripts, e.g. C_{2}H_{5}OH^{+}).
!
!-----------------------------------------------------------------
module labelschem
 use labels, only:lenlabel
 use asciiutils, only:lcase,ucase,string_replace,isdigit
 implicit none

 public :: format_chemistry_label

 private :: chemistry_label_alias,strip_ion_suffix,parse_krome_species_name,prettify_chemistry_fallback
 private :: match_element_symbol,append_subscript_digits,capitalise_element
 private :: is_two_letter_element,is_one_letter_element,is_lowercase_letter
 private :: is_blocked_diatomic_pair,match_special_two_letter,allow_standard_two_letter

contains

!-----------------------------------------------------------------
! format raw Krome/Phantom species name for display in label()
!-----------------------------------------------------------------
function format_chemistry_label(raw) result(display)
 character(len=*), intent(in) :: raw
 character(len=lenlabel) :: display
 character(len=lenlabel) :: name,body,ion
 logical :: parsed

 name = lcase(trim(raw))
 if (len_trim(name) == 0) then
    display = ' '
    return
 endif
 if (index(name,'{') > 0) then
    display = trim(raw)
    return
 endif

 display = chemistry_label_alias(name)
 if (len_trim(display) > 0) return

 call strip_ion_suffix(name, body, ion)
 parsed = parse_krome_species_name(body, display)
 if (.not.parsed .or. len_trim(display) == 0) then
    display = prettify_chemistry_fallback(body)
 endif
 display = trim(display)//trim(ion)

end function format_chemistry_label

!-----------------------------------------------------------------
! known aliases that the parser should not treat literally
!-----------------------------------------------------------------
function chemistry_label_alias(name) result(alias)
 character(len=*), intent(in) :: name
 character(len=lenlabel) :: alias

 alias = ' '
 select case (trim(name))
 case ('h')
    alias = 'H'
 case ('av')
    alias = 'A_V'
 case ('co')
    alias = 'CO'
 case ('no')
    alias = 'NO'
 end select

end function chemistry_label_alias

!-----------------------------------------------------------------
! strip trailing _plus / _minus ion suffix into TeX charge string
!-----------------------------------------------------------------
subroutine strip_ion_suffix(name, body, ion)
 character(len=*), intent(inout) :: name
 character(len=*), intent(out) :: body,ion
 integer :: n

 body = trim(name)
 ion = ' '
 n = len_trim(body)
 if (n > 5) then
    if (body(n-4:n) == '_plus') then
       ion = '^{+}'
       body = body(1:n-5)
    endif
 endif
 if (len_trim(ion) == 0) then
    n = len_trim(body)
    if (n > 6) then
       if (body(n-5:n) == '_minus') then
          ion = '^{-}'
          body = body(1:n-6)
       endif
    endif
 endif
 name = trim(body)

end subroutine strip_ion_suffix

!-----------------------------------------------------------------
! parse concatenated element symbols and stoichiometric digits
!-----------------------------------------------------------------
logical function parse_krome_species_name(name, formatted)
 character(len=*), intent(in) :: name
 character(len=lenlabel), intent(out) :: formatted
 integer :: i,nlen,nsym,ndig
 character(len=3) :: sym
 character(len=12) :: digits
 logical :: anytok

 formatted = ' '
 nlen = len_trim(name)
 if (nlen < 1) then
    parse_krome_species_name = .false.
    return
 endif

 i = 1
 anytok = .false.
 do while (i <= nlen)
    if (match_element_symbol(name, i, nsym, sym)) then
       formatted = trim(formatted)//trim(sym)
       i = i + nsym
       ndig = append_subscript_digits(name, i, nlen, digits)
       if (ndig > 0) then
          formatted = trim(formatted)//'_{'//trim(digits)//'}'
          i = i + ndig
       endif
       anytok = .true.
    elseif (is_lowercase_letter(name, i, nlen)) then
       parse_krome_species_name = .false.
       return
    else
       i = i + 1
    endif
 enddo

 parse_krome_species_name = anytok

end function parse_krome_species_name

!-----------------------------------------------------------------
! match element symbol at position i (try two letters, then one)
!-----------------------------------------------------------------
logical function match_element_symbol(name, i, nsym, sym)
 character(len=*), intent(in) :: name
 integer, intent(in) :: i
 integer, intent(out) :: nsym
 character(len=*), intent(out) :: sym
 character(len=2) :: pair
 integer :: nlen

 match_element_symbol = .false.
 sym = ' '
 nsym = 0
 nlen = len_trim(name)
 if (i < 1 .or. i > nlen) return

 if (match_special_two_letter(name, i, nlen, nsym, sym)) then
    match_element_symbol = .true.
    return
 endif

 if (i + 1 <= nlen) then
    pair = lcase(name(i:i+1))
    if (is_two_letter_element(pair)) then
       if (.not.is_blocked_diatomic_pair(pair)) then
          if (allow_standard_two_letter(name, i, nlen, pair)) then
             sym = capitalise_element(pair)
             nsym = 2
             match_element_symbol = .true.
             return
          endif
       endif
    endif
 endif

 if (is_one_letter_element(name(i:i))) then
    sym = ucase(name(i:i))
    nsym = 1
    match_element_symbol = .true.
 endif

end function match_element_symbol

!-----------------------------------------------------------------
! context-specific two-letter matches (override generic IUPAC rules)
!-----------------------------------------------------------------
logical function match_special_two_letter(name, i, nlen, nsym, sym)
 character(len=*), intent(in) :: name
 integer, intent(in) :: i,nlen
 integer, intent(out) :: nsym
 character(len=*), intent(out) :: sym
 character(len=2) :: pair

 match_special_two_letter = .false.
 sym = ' '
 nsym = 0
 if (i + 1 > nlen) return

 pair = lcase(name(i:i+1))

 !--silicon: SiO, SiS, H_{2}SiO, ... (not S+I+...)
 if (pair == 'si') then
    if (i + 2 > nlen .or. is_lowercase_letter(name, i+2, nlen)) then
       sym = 'Si'
       nsym = 2
       match_special_two_letter = .true.
    endif
    return
 endif

 !--helium hydride: HeH (not mis-parsed He)
 if (pair == 'he') then
    if (i + 2 <= nlen) then
       if (lcase(name(i+2:i+2)) == 'h') then
          if (i + 3 > nlen .or. .not.is_lowercase_letter(name, i+3, nlen)) then
             sym = 'He'
             nsym = 2
             match_special_two_letter = .true.
          endif
       endif
    endif
    return
 endif

 !--halogen monoxides: ClO, BrO (not C+l+o)
 if (pair == 'cl' .or. pair == 'br') then
    if (i + 2 <= nlen) then
       if (lcase(name(i+2:i+2)) == 'o') then
          if (i + 3 > nlen .or. .not.is_lowercase_letter(name, i+3, nlen)) then
             if (pair == 'cl') then
                sym = 'Cl'
             else
                sym = 'Br'
             endif
             nsym = 2
             match_special_two_letter = .true.
          endif
       endif
    endif
 endif

end function match_special_two_letter

!-----------------------------------------------------------------
! standard IUPAC two-letter element at position i
!-----------------------------------------------------------------
logical function allow_standard_two_letter(name, i, nlen, pair)
 character(len=*), intent(in) :: name
 integer, intent(in) :: i,nlen
 character(len=2), intent(in) :: pair

 allow_standard_two_letter = .false.
 if (i < 1 .or. i + 1 > nlen) return

 !--h2cl: halogen at end after stoichiometric digit
 if (i > 1) then
    if (isdigit(name(i-1:i-1))) then
       if (pair == 'cl' .or. pair == 'br') then
          if (i + 2 > nlen) allow_standard_two_letter = .true.
       endif
       return
    endif
 endif

 !--do not take Co/No inside a formula (hcnoh, h2co)
 if (is_lowercase_letter(name, i+2, nlen)) return

 allow_standard_two_letter = .true.

end function allow_standard_two_letter

!-----------------------------------------------------------------
! true if position is a lowercase a-z letter in name
!-----------------------------------------------------------------
logical function is_lowercase_letter(name, pos, nlen)
 character(len=*), intent(in) :: name
 integer, intent(in) :: pos,nlen
 integer :: ia

 is_lowercase_letter = .false.
 if (pos < 1 .or. pos > nlen) return
 ia = iachar(name(pos:pos))
 if (ia >= iachar('a') .and. ia <= iachar('z')) is_lowercase_letter = .true.

end function is_lowercase_letter

!-----------------------------------------------------------------
! two-letter sequences that are atom pairs in Krome names, not elements
!-----------------------------------------------------------------
logical function is_blocked_diatomic_pair(pair)
 character(len=2), intent(in) :: pair
 character(len=2) :: s

 s = lcase(pair)
 is_blocked_diatomic_pair = .false.
 select case (s)
 case ('co','no','cn','nc','on','oc','oh','os','as','in','at','po','sc','cs','sn','ho')
    is_blocked_diatomic_pair = .true.
 end select

end function is_blocked_diatomic_pair

!-----------------------------------------------------------------
! return length of digit run starting at i
!-----------------------------------------------------------------
integer function append_subscript_digits(name, i, nlen, digits)
 character(len=*), intent(in) :: name
 integer, intent(in) :: i,nlen
 character(len=*), intent(out) :: digits
 integer :: j

 digits = ' '
 append_subscript_digits = 0
 if (i < 1 .or. i > nlen) return
 if (.not.isdigit(name(i:i))) return
 j = i
 do while (j <= nlen .and. isdigit(name(j:j)))
    j = j + 1
 enddo
 append_subscript_digits = j - i
 digits = name(i:j-1)

end function append_subscript_digits

!-----------------------------------------------------------------
! capitalise element symbol (e.g. fe -> Fe, c -> C)
!-----------------------------------------------------------------
function capitalise_element(sym) result(cap)
 character(len=*), intent(in) :: sym
 character(len=3) :: cap

 if (len_trim(sym) >= 2) then
    cap = ucase(sym(1:1))//lcase(sym(2:2))
 else
    cap = ucase(sym(1:1))
 endif

end function capitalise_element

!-----------------------------------------------------------------
! fallback label if formula parse fails
!-----------------------------------------------------------------
function prettify_chemistry_fallback(name) result(pretty)
 character(len=*), intent(in) :: name
 character(len=lenlabel) :: pretty

 pretty = trim(name)
 call string_replace(pretty,'_',' ')

end function prettify_chemistry_fallback

!-----------------------------------------------------------------
! two-letter element symbols (IUPAC)
!-----------------------------------------------------------------
logical function is_two_letter_element(pair)
 character(len=2), intent(in) :: pair
 character(len=2) :: s

 s = lcase(pair)
 is_two_letter_element = .false.
 select case (s)
 case ('he','li','be','ne','na','mg','al','si','cl','ar','ca','sc','ti','cr','mn','fe','co',&
       'ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','zr','nb','mo','ru','rh','pd',&
       'ag','cd','in','sn','sb','te','xe','cs','ba','la','ce','pr','nd','pm','sm','eu','gd',&
       'tb','dy','ho','er','tm','yb','lu','hf','ta','re','os','ir','pt','au','hg','tl','pb',&
       'bi','po','at','rn','fr','ra','ac','th','pa','np','pu','am','cm','bk','cf','es','fm',&
       'md','no','lr')
    is_two_letter_element = .true.
 end select

end function is_two_letter_element

!-----------------------------------------------------------------
! single-letter element symbols
!-----------------------------------------------------------------
logical function is_one_letter_element(letter)
 character(len=1), intent(in) :: letter
 character(len=1) :: c

 c = lcase(letter)
 is_one_letter_element = .false.
 select case (c)
 case ('h','b','c','n','o','f','p','s','k','v','y','i','w','u')
    is_one_letter_element = .true.
 end select

end function is_one_letter_element

end module labelschem
