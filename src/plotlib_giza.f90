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
!  Copyright (C) 2005-2010 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
!  The plotlib module in SPLASH provides a consistent api so that SPLASH
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
      plot_svp=>giza_set_viewport, &
      plot_swin=>giza_set_window, &
      giza_get_surface_size, &
      plot_box=>giza_box, &
      plot_circ=>giza_circle, &
      plot_line=>giza_line, &
      giza_open_device, &
      plot_close=>giza_close_device, &
      plot_page=>giza_change_page, &
      plot_slw=>giza_set_line_width, &
      giza_open_device, &
      giza_open_device_size, &
      plot_sch=>giza_set_character_height, &
      plot_qch=>giza_get_character_height, &
      plot_scf=>giza_set_font, &
      plot_annotate=>giza_annotate, &
      giza_get_character_size, &
      plot_bbuf=>giza_begin_buffer,&
      plot_ebuf=>giza_end_buffer, &
      plot_qwin=>giza_get_window, &
      giza_get_viewport, &
      plot_sci=>giza_set_colour_index, &
      plot_scr=>giza_set_colour_representation, &
      plot_qcr=>giza_get_colour_representation, &
      plot_qlw=>giza_get_line_width,&
      plot_qls=>giza_get_line_style,&
      plot_sls=>giza_set_line_style,&
      plot_qci=>giza_get_colour_index, &
      plot_end=>giza_close_device, &
      plot_ptxt=>giza_ptext, &
      plot_slc=>giza_set_line_cap, &
      plot_qlc=>giza_get_line_cap, &
      giza_stop_prompting, &
      giza_start_prompting, &
      plot_curs=>giza_get_key_press, &
      plot_qcur=>giza_device_has_cursor, &
      giza_set_colour_table, &
      giza_render, &
      plot_wnad=>giza_set_window_equal_scale, &
      plot_pt1=>giza_single_point, &
      plot_pt=>giza_points, &
      plot_errb=>giza_error_bars, &
      plot_sfs=>giza_set_fill, &
      plot_qfs=>giza_get_fill, &
      giza_left_click_f, &
      plot_band=>giza_band, &
      giza_vector, &
      giza_arrow,  &
      giza_format_number, &
      giza_query_device
  implicit none
  
  character(len=1),parameter :: plot_left_click = giza_left_click_f
public

contains

!---------------------------------------------
! query whether or not the library is PGPLOT
!---------------------------------------------
logical function plot_lib_is_pgplot()
 implicit none
 
 plot_lib_is_pgplot = .false.

end function plot_lib_is_pgplot

!---------------------------------------------
! initialise the plotting library
!---------------------------------------------
subroutine plot_init(devicein, ierr, papersizex, aspectratio)
 implicit none

 character(len=*),intent(in) :: devicein
 integer,intent(out)         :: ierr
 real, intent(in), optional  :: papersizex,aspectratio
 real                        :: widthcm,heightcm
 real, parameter             :: inch_to_cm = 2.54
 
 if (present(papersizex)) then
    widthcm = papersizex*inch_to_cm
    if (present(aspectratio)) then
       heightcm = widthcm*aspectratio
    else
       heightcm = widthcm/sqrt(2.)
    endif
    ierr = giza_open_device_size(devicein, 'splash', widthcm, heightcm)
 else
    ierr = giza_open_device(devicein,'splash')
 endif

 if(ierr.eq.0) then
    call giza_stop_prompting
 endif
end subroutine plot_init

subroutine plot_imag(a, idim, jdim, i1, i2, j1, j2,&
                     a1, a2, tr)
  integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
  real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
  real               :: affine(6)

  call convert_tr_to_affine(tr,affine)
  call giza_render(idim,jdim,a,i1-1,i2-1,j1-1,j2-1,a1,a2,affine)
  
end subroutine plot_imag

subroutine plot_ctab(l,r,g,b,nc,contra,bright)
  implicit none
  integer,intent(in) :: nc
  real,intent(in)    :: l(nc),r(nc),g(nc),b(nc),contra,bright

  call giza_set_colour_table(l,r,g,b,nc)

end subroutine plot_ctab

subroutine plot_qvsz(units,x1,x2,y1,y2)
  implicit none

  real, intent(out)   :: x1,x2,y1,y2
  integer, intent(in) :: units

  call giza_get_surface_size(x1,x2,y1,y2)
end subroutine plot_qvsz

subroutine plot_bins(nbin,x,data,centre)
  integer, intent(in)               :: nbin 
  real, dimension(nbin), intent(in) :: x, data
  logical, intent(in)               :: centre

  print*,' WARNING: plot_bins not implemented in giza'

end subroutine plot_bins

subroutine plot_qvp(units, x1, x2, y1, y2)
  integer,intent(in) :: units
  real,intent(out)   :: x1, x2, y1, y2

  if(units.eq.0) then
     call giza_get_viewport(units,x1,x2,y1,y2)
  else if(units.eq.3) then
     call giza_get_viewport(1,x1,x2,y1,y2)
  else
     call giza_get_viewport(units,x1,x2,y1,y2)
  end if

end subroutine plot_qvp

subroutine plot_qcs(units,xch,ych)
  integer,intent(in) :: units
  real,intent(out)   :: xch,ych

  if(units.eq.0) then
     call giza_get_character_size(0,xch,ych)
  else if(units.eq.4) then
     call giza_get_character_size(1,xch,ych)
  else
     call giza_get_character_size(units,xch,ych)
  endif
end subroutine plot_qcs

subroutine plot_qcol(icolmin,icolmax)
  integer,intent(out) :: icolmin,icolmax

  icolmin = 0
  icolmax = 20
end subroutine plot_qcol

subroutine plot_qcir(icolmin,icolmax)
  integer,intent(out) :: icolmin,icolmax

  icolmin = 0
  icolmax = 19
end subroutine plot_qcir

subroutine plot_scir(icilo, icihi)
  integer,intent(in) :: icilo,icihi

  print*,' WARNING: plot_scir not implemented in giza'

end subroutine plot_scir

subroutine plot_rect(x1,x2,y1,y2)
  implicit none     
  real,intent(in) :: x1,x2,y1,y2

  print*,' WARNING: plot_rect not implemented in giza'

end subroutine plot_rect

subroutine plot_arro(x1,y1,x2,y2)
  implicit none
  real,intent(in) :: x1,y1,x2,y2

  call giza_arrow(x1,y1,x2,y2)

end subroutine plot_arro

subroutine plot_qtxt(x,y,angle,fjust,text,xbox,ybox)
  implicit none
  real,intent(in)               :: x, y, angle, fjust
  character(len=*),intent(in)   :: text
  real,intent(out),dimension(4) :: xbox,ybox
  
  print*,' WARNING: plot_qtxt not implemented in giza'
  xbox = 0.
  ybox = 0.
  
end subroutine plot_qtxt

subroutine plot_stbg(bg)
  implicit none
  integer,intent(in) :: bg

  print*,' WARNING: plot_stbg not implemented in giza'

end subroutine plot_stbg

subroutine plot_funx(fx,n,ymin,ymax,pgflags)
  implicit none      
  real,external            :: fx
  integer,intent(in)       :: n,pgflags
  real,intent(in)          :: ymin,ymax

  print*,' WARNING: plot_funx not implemented in giza'

end subroutine plot_funx

subroutine plot_label(xlbl,ylbl,toplbl)
  implicit none
  character(len=*),intent(in) :: xlbl,ylbl,toplbl

  print*,' WARNING: plot_label not implemented in giza'

end subroutine plot_label

subroutine plot_scrn(ci,name,ier)
  implicit none
  integer,intent(in)          :: ci
  character(len=*),intent(in) :: name
  integer,intent(out)         :: ier

  print*,' WARNING: plot_scrn not implemented in giza'
  ier = 1

end subroutine plot_scrn

subroutine plot_poly(n,xpts,ypts)
  implicit none
  integer,intent(in) :: n
  real,intent(in)    :: xpts(*),ypts(*)

  print*,' WARNING: plot_poly not implemented in giza'

end subroutine plot_poly

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
!     if (plot_qcur()) then
!        value = 'YES'
!     else
!        value = 'NO'
!     endif
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
  real :: x
  
  x = m*10.**pp
!  print*,'pp = ',pp,' m = ',m
  call giza_format_number(x,3,string)
  nc = len_trim(string)

end subroutine plot_numb

subroutine plot_text(x,y,text)
  implicit none
  real,intent(in)             :: x,y
  character(len=*),intent(in) :: text

  call plot_ptxt(x,y,0.,0.,text)

end subroutine plot_text

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

subroutine plot_conb(a,idim,jdim,i1,i2,j1,j2,c,nc,tr, &
     blank)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
  real,intent(in)    :: a(idim,jdim),c(*),tr(6),blank

  print*,' WARNING: plot_conb not implemented in giza'

end subroutine plot_conb

subroutine plot_cons(a,idim,jdim,i1,i2,j1,j2,c,nc,tr)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
  real,intent(in)    :: a(idim,jdim),c(*),tr(6)
  
  print*,' WARNING: plot_cons not implemented in giza'

end subroutine plot_cons

subroutine plot_conl(a,idim,jdim,i1,i2,j1,j2,c,tr, &
     label,intval,mininit)
  implicit none
  integer,intent(in)          :: idim,jdim,i1,i2,j1,j2,intval,mininit
  real,intent(in)             :: a(idim,jdim),c,tr(6)
  character(len=*),intent(in) :: label

  print*,' WARNING: plot_conl not implemented in giza'

end subroutine plot_conl

subroutine plot_sah(fs, angle, cutback)
  implicit none
  integer, intent(in) :: fs
  real, intent(in)    :: angle, cutback
  
  print*,' WARNING: plot_sah not implemented in giza'

end subroutine plot_sah

subroutine plot_vect(a,b,idim,jdim,i1,i2,j1,j2,c,nc,tr, &
     blank)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
  real,intent(in)    :: a(idim,jdim),b(idim,jdim),tr(6),blank,c
  real               :: affine(6)

  call convert_tr_to_affine(tr,affine)
  call giza_vector(idim,jdim,a,b,i1-1,i2-1,j1-1,j2-1,c,nc,affine,blank)

end subroutine plot_vect

subroutine plot_pixl(ia,idim,jdim,i1,i2,j1,j2, &
     x1,x2,y1,y2)
  implicit none
  integer,intent(in) :: idim,jdim,i1,i2,j1,j2
  integer,intent(in) :: ia(idim,jdim)
  real,intent(in)    :: x1,x2,y1,y2

  print*,' WARNING: plot_env not implemented in giza'

end subroutine plot_pixl

subroutine plot_env(xmin,xmax,ymin,ymax,just,axis)
  implicit none
  real,intent(in)    :: xmin,xmax,ymin,ymax
  integer,intent(in) :: just,axis

  print*,' WARNING: plot_env not implemented in giza'

end subroutine plot_env

subroutine plot_pap(width,aspect)
  implicit none
  real,intent(in) :: width,aspect

  print*,' WARNING: plot_pap not implemented in giza'

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

