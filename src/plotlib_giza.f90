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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
!  The plotlib module in SPLASH provides a consistent API so that SPLASH
!  can be compiled against different graphics libraries as the backend
!
!  This version provides an interface to giza, a plotting
!   library written by Daniel Price & James Wetter.
!
!  Giza implements basic 2D plotting functionality
!  on top of the cairo graphics library
!
! Interface written by James Wetter and Daniel Price (2010)
!---------------------------------------------------------------------------
module plotlib
  use giza, only: &
      plot_arro=>giza_arrow, &
      plot_annotate=>giza_annotate, &
      plot_band=>giza_band, &
      plot_bbuf=>giza_begin_buffer,&
      plot_box=>giza_box, &
      plot_circ=>giza_circle, &
      plot_close=>giza_close_device, &
      plot_curs=>giza_get_key_press, &
      plot_ebuf=>giza_end_buffer, &
      plot_end=>giza_close_device, &
      plot_env=>giza_set_environment, &
      plot_errb=>giza_error_bars, &
      plot_funx=>giza_function_x, &
      plot_label=>giza_label, &
      plot_line=>giza_line, &
      plot_lcur=>giza_mark_line, &
      plot_olin=>giza_mark_points, &
      plot_ncur=>giza_mark_points_ordered, &
      plot_page=>giza_change_page, &
      plot_poly=>giza_polygon, &
      plot_pt1=>giza_single_point, &
      plot_pt=>giza_points, &
      plot_ptxt=>giza_ptext, &
      plot_qch=>giza_get_character_height, &
      plot_qci=>giza_get_colour_index, &
      plot_qcir=>giza_get_colour_index_range, &
      plot_qcr=>giza_get_colour_representation, &
      plot_qfs=>giza_get_fill, &
      plot_qlw=>giza_get_line_width,&
      plot_qls=>giza_get_line_style,&
      plot_qlc=>giza_get_line_cap, &
      plot_qtxt=>giza_qtext, &
      plot_qwin=>giza_get_window, &
      plot_rect=>giza_rectangle, &
      plot_sah=>giza_set_arrow_style, &
      plot_scf=>giza_set_font, &
      plot_sch=>giza_set_character_height, &
      plot_sci=>giza_set_colour_index, &
      plot_scir=>giza_set_colour_index_range, &
      plot_scr=>giza_set_colour_representation, &
      plot_sfs=>giza_set_fill, &
      plot_slc=>giza_set_line_cap, &
      plot_sls=>giza_set_line_style,&
      plot_slw=>giza_set_line_width, &
      plot_stbg=>giza_set_text_background, &
      plot_svp=>giza_set_viewport, &
      plot_swin=>giza_set_window, &
      plot_text=>giza_text, &
      plot_wnad=>giza_set_window_equal_scale, &
      plot_qcur=>giza_device_has_cursor, &
      plot_rgb_from_table=>giza_rgb_from_table, &
      giza_contour,            &
      giza_get_character_size, &
      giza_get_surface_size,   &
      giza_get_viewport,       &
      giza_open_device,        &
      giza_open_device_size,   &
      giza_render,             &
      giza_render_gray,        &
      giza_render_transparent, &
      giza_set_colour_table,   &
      giza_stop_prompting,     &
      giza_start_prompting,    &
      plot_left_click=>giza_left_click_f,     &
      plot_right_click=>giza_right_click_f,   &
      plot_middle_click=>giza_middle_click_f, &
      plot_shift_click=>giza_shift_click_f,   &
      plot_scroll_up=>giza_scroll_up_f,       &
      plot_scroll_down=>giza_scroll_down_f,   &
      plot_scroll_left=>giza_scroll_left_f,   &
      plot_scroll_right=>giza_scroll_right_f, &
      giza_vector,             &
      giza_format_number,      &
      giza_query_device,       &
      giza_draw_pixels, &
      giza_colour_index_min,&
      giza_colour_index_max,&
      giza_extend_pad,&
      giza_extend_repeat, &
      giza_extend_reflect, &
      giza_extend_none
  implicit none
  logical, parameter :: plotlib_is_pgplot = .false.
  logical, parameter :: plotlib_supports_alpha = .true.
  integer, parameter :: plotlib_maxlinestyle = 6
  integer, parameter :: plotlib_maxfillstyle = 5
  integer, parameter :: plotlib_maxlinecolour = 16
  integer, parameter :: plotlib_extend_pad = giza_extend_pad
  integer, parameter :: plotlib_extend_repeat = giza_extend_repeat
  integer, parameter :: plotlib_extend_reflect = giza_extend_reflect
  integer, parameter :: plotlib_extend_none = giza_extend_none

public

contains

!---------------------------------------------
! initialise the plotting library
!---------------------------------------------
subroutine plot_init(devicein, ierr, papersizex, aspectratio, paperunits)
 use giza, only:giza_units_inches,giza_units_pixels,giza_units_mm
 implicit none

 character(len=*),intent(in)   :: devicein
 integer,intent(out)           :: ierr
 real, intent(in), optional    :: papersizex,aspectratio
 integer, intent(in), optional :: paperunits
 real                          :: width,height
 integer                       :: units, id

 if (present(papersizex)) then
    width = papersizex
    if (present(aspectratio)) then
       height = width*aspectratio
    else
       height = width/sqrt(2.)
    endif
    if (present(paperunits)) then
       select case(paperunits)
       case(0)
          units = giza_units_pixels
       case(1)
          units = giza_units_inches
       case(2)
          units = giza_units_mm
          width = 10.*width
          height = 10.*height
       end select
    else
       units = giza_units_inches
    endif
    id = giza_open_device_size(devicein, 'splash', width, height, units)
 else
    id = giza_open_device(devicein,'splash')
 endif
 ! id<0 should return an error, but +ve id is OK
 select case(id)
 case(:-1)
    ierr = id
 case default
    ierr = 0
 end select

 if(ierr.eq.0) then
    call giza_stop_prompting
 endif
end subroutine plot_init

subroutine plot_gray(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr, iextend)
  integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
  real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
  real               :: affine(6)
  integer, intent(in), optional :: iextend

  call convert_tr_to_affine(tr,affine)

  if (present(iextend)) then
     call giza_render_gray(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,a1,a2,iextend,affine)
  else
     call giza_render_gray(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,a1,a2,0,affine)
  endif

end subroutine plot_gray

subroutine plot_imag(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr, iextend)
  integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
  real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
  real               :: affine(6)
  integer, intent(in), optional :: iextend

  call convert_tr_to_affine(tr,affine)
  
  if (present(iextend)) then
     call giza_render(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,a1,a2,iextend,affine)
  else
     call giza_render(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,a1,a2,0,affine)
  endif

end subroutine plot_imag

subroutine plot_imag_alpha(dat, alpha, idim, jdim, i1, i2, j1, j2, a1, a2, tr, iextend)
  integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
  real,intent(in)    :: dat(IDIM,JDIM), alpha(IDIM,JDIM), A1, A2, TR(6)
  real               :: affine(6)
  integer, intent(in), optional :: iextend

  call convert_tr_to_affine(tr,affine)
  
  if (present(iextend)) then
     call giza_render(idim,jdim,dat,alpha,i1-1,i2-1,j1-1,j2-1,a1,a2,iextend,affine)
  else
     call giza_render(idim,jdim,dat,alpha,i1-1,i2-1,j1-1,j2-1,a1,a2,0,affine)
  endif

end subroutine plot_imag_alpha

subroutine plot_imag_transparent(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr)
  integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
  real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
  real               :: affine(6)

  call convert_tr_to_affine(tr,affine)
  call giza_render_transparent(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,a1,a2,0,affine)

end subroutine plot_imag_transparent

subroutine plot_ctab(l,r,g,b,nc,contra,bright)
  implicit none
  integer,intent(in) :: nc
  real,intent(in)    :: l(nc),r(nc),g(nc),b(nc),contra,bright

  call giza_set_colour_table(l,r,g,b,nc,contra,bright)

end subroutine plot_ctab

subroutine plot_qvsz(units,x1,x2,y1,y2)
  use giza, only:giza_get_paper_size,giza_units_device
  implicit none

  real, intent(out)   :: x1,x2,y1,y2
  integer, intent(in) :: units

  x1 = 0.
  y1 = 0.
  call giza_get_paper_size(units,x2,y2)

end subroutine plot_qvsz

subroutine plot_bins(nbin,x,data,centre)
  integer, intent(in)               :: nbin
  real, dimension(nbin), intent(in) :: x, data
  logical, intent(in)               :: centre

  print*,' WARNING: plot_bins not implemented in giza'

end subroutine plot_bins

subroutine plot_qvp(units, x1, x2, y1, y2)
 implicit none
 integer,intent(in) :: units
 real,intent(out)   :: x1, x2, y1, y2

 call giza_get_viewport(units_giza(units),x1,x2,y1,y2)

end subroutine plot_qvp

subroutine plot_qcs(units,xch,ych)
 implicit none
 integer,intent(in) :: units
 real,intent(out)   :: xch,ych

 call giza_get_character_size(units_giza(units),xch,ych)

end subroutine plot_qcs

subroutine plot_qcol(icolmin,icolmax)
  integer,intent(out) :: icolmin,icolmax

  icolmin = giza_colour_index_min
  icolmax = giza_colour_index_max

end subroutine plot_qcol

subroutine plot_scrn(ci,name,ier)
  implicit none
  integer,intent(in)          :: ci
  character(len=*),intent(in) :: name
  integer,intent(out)         :: ier

  print*,' WARNING: plot_scrn not implemented in giza'
  ier = 1

end subroutine plot_scrn

subroutine plot_qinf(item,value,length)
  implicit none
  character(len=*),intent(in)  :: item
  character(len=*),intent(out) :: value
  integer,intent(out)          :: length
  character(len=10) :: datestring,timestring

  select case(item)
  case('VERSION','version')
     value = 'giza-0.1'
  case('STATE','state')
     print*,' WARNING: query for STATE not yet implemented in giza'
  case('USER','user')
     print*,' WARNING: query for USER not yet implemented in giza'
  case('NOW','now')
     call date_and_time(datestring,timestring)
     value = datestring(7:8)//'-'//datestring(5:6)//'-'//datestring(1:4)// &
             ' '//timestring(1:2)//':'//timestring(3:4)
   case('DEVICE','device')
     print*,' WARNING: query for DEVICE not yet implemented in giza'
  case('FILE','file')
     print*,' WARNING: query for FILE not yet implemented in giza'
  case('TYPE','type')
     call giza_query_device('type',value)
  case('DEV/TYPE','dev/type')
     print*,' WARNING: query for DEV/TYPE not yet implemented in giza'
  case('HARDCOPY','hardcopy')
     call giza_query_device('hardcopy',value)
  case('TERMINAL','terminal')
     !--in giza the current device is never the terminal
     value = 'NO'
  case('CURSOR','cursor')
     call giza_query_device('cursor',value)
  case('SCROLL','scroll')
     !--no scroll capability in any current giza devices
     value = 'NO'
  case default
     value = ' '
  end select
  length = len_trim(value)

end subroutine plot_qinf

subroutine plot_numb(m,pp,form,string,nc)
  implicit none
  integer,intent(in)           :: m,pp,form
  character(len=*),intent(out) :: string
  integer,intent(out)          :: nc

  call giza_format_number(m,pp,form,string)
  nc = len_trim(string)

end subroutine plot_numb

subroutine plot_set_opacity(alpha)
  implicit none
  real, intent(in)           :: alpha
  integer :: ci
  real    :: red,green,blue

  call plot_qci(ci)
  call plot_qcr(ci,red,green,blue)
  call plot_scr(ci,red,green,blue,alpha)

end subroutine plot_set_opacity


subroutine plot_err1(dir,x,y,e,t)
  implicit none
  integer,intent(in) :: dir
  real,intent(in)    :: x,y,e
  real,intent(in)    :: t
  real, dimension(1) :: xi,yi,ei

  xi(1) = x
  yi(1) = y
  ei(1) = e
  call plot_errb(dir,1,xi,yi,ei,t)

end subroutine plot_err1

subroutine plot_conb(a,idim,jdim,i1,i2,j1,j2,c,nc,tr,blank)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
  real,intent(in)    :: a(idim,jdim),c(*),tr(6),blank
  real               :: affine(6)

  print*,' WARNING: blanking in contouring not implemented in giza'
  call convert_tr_to_affine(tr,affine)
  call giza_contour(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,nc,c,affine)

end subroutine plot_conb

subroutine plot_cons(a,idim,jdim,i1,i2,j1,j2,c,nc,tr)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
  real,intent(in)    :: a(idim,jdim),c(*),tr(6)
  real               :: affine(6)

  call convert_tr_to_affine(tr,affine)
  call giza_contour(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,nc,c,affine)

end subroutine plot_cons

subroutine plot_conl(a,idim,jdim,i1,i2,j1,j2,c,tr,label,intval,mininit)
  implicit none
  integer,intent(in)          :: idim,jdim,i1,i2,j1,j2,intval,mininit
  real,intent(in)             :: a(idim,jdim),c,tr(6)
  character(len=*),intent(in) :: label
  real               :: affine(6)
  integer, parameter :: nc = 1
  real, dimension(nc) :: clevel

  clevel(1) = c
  print*,' WARNING: labelled coutouring not implemented in giza'
  call convert_tr_to_affine(tr,affine)
  call giza_contour(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,nc,clevel,affine)
  print*,'nc = ',nc

end subroutine plot_conl

subroutine plot_vect(a,b,idim,jdim,i1,i2,j1,j2,c,nc,tr,blank)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
  real,intent(in)    :: a(idim,jdim),b(idim,jdim),tr(6),blank,c
  real               :: affine(6)

  call convert_tr_to_affine(tr,affine)
  call giza_vector(idim,jdim,a,b,i1-1,i2-1,j1-1,j2-1,c,nc,affine,blank)

end subroutine plot_vect

subroutine plot_pixl(ia,idim,jdim,i1,i2,j1,j2,x1,x2,y1,y2)
  use giza, only:giza_draw_pixels
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2
  integer,intent(in) :: ia(idim,jdim)
  real,intent(in)    :: x1,x2,y1,y2

  call giza_draw_pixels(IDIM, JDIM, IA, I1-1, I2-1, J1-1, J2-1, X1, X2, Y1, Y2, 0)

end subroutine plot_pixl

subroutine plot_pap(widthin,aspect,paperunits)
  use giza, only:giza_set_paper_size
  use giza, only:giza_units_inches,giza_units_pixels,giza_units_mm
  implicit none
  real,intent(in) :: widthin,aspect
  integer, intent(in), optional :: paperunits
  integer :: units
  real    :: width
  
  width = widthin
  units = giza_units_inches

  if (present(paperunits)) then
     select case(paperunits)
     case(0)
        units = giza_units_pixels
     case(1)
        units = giza_units_inches
     case(2)
        units = giza_units_mm
        width = 0.1*width
     end select
  endif
  
  call giza_set_paper_size(units,width,width*aspect)  

end subroutine plot_pap

!
!--this subroutine can be called  to
!  make sure that the viewport lies exactly on
!  pixel boundaries.
!
!  unnecessary for giza
!
subroutine plot_set_exactpixelboundaries()
 implicit none

end subroutine plot_set_exactpixelboundaries

!------------------------------------------------------------
! Function to convert PGPLOT units value to giza units value
!------------------------------------------------------------
 integer function units_giza(pgplotunits)
  use giza, only:giza_units_normalized,giza_units_inches, &
                 giza_units_mm,giza_units_pixels,giza_units_world
  implicit none
  integer, intent(in) :: pgplotunits

  select case(pgplotunits)
  case(0)
     units_giza = giza_units_normalized
  case(1)
     units_giza = giza_units_inches
  case(2)
     units_giza = giza_units_mm
  case(3)
     units_giza = giza_units_pixels
  case(4)
     units_giza = giza_units_world
  case default  ! giza will give an error
     units_giza = pgplotunits
  end select

 end function units_giza

subroutine convert_tr_to_affine(tr,affine)
 implicit none
 real, dimension(6), intent(in)  :: tr
 real, dimension(6), intent(out) :: affine

 affine(1) = TR(2)
 affine(2) = TR(3)
 affine(3) = TR(5)
 affine(4) = TR(6)
 affine(5) = TR(1) + 0.5 * TR(2)
 affine(6) = TR(4) + 0.5 * TR(6)

end subroutine convert_tr_to_affine

end module plotlib

