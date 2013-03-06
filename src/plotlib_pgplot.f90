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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
!  The plotlib module in SPLASH provides a consistent API so that SPLASH
!  can be compiled against different graphics libraries as the backend
!
!  This version provides an interface to Tim Pearson's PGPLOT,
!   which was the original backend used in SPLASH v1.x. Thus,
!   the functions mostly translate directly to PGPLOT equivalents,
!   and basically this is a Fortran 90 interface for PGPLOT.
!
! Interface written by James Wetter and Daniel Price (2010)
!---------------------------------------------------------------------------
module plotlib
  implicit none
  logical, parameter :: plotlib_is_pgplot = .true.
  logical, parameter :: plotlib_supports_alpha = .false.
  integer, parameter :: plotlib_maxlinestyle = 5
  integer, parameter :: plotlib_maxfillstyle = 5
  integer, parameter :: plotlib_maxlinecolour = 16
  integer, parameter :: plotlib_extend_pad = 1     ! not implemented in PGPLOT
  integer, parameter :: plotlib_extend_repeat = 2  ! not implemented in PGPLOT
  integer, parameter :: plotlib_extend_reflect = 3 ! not implemented in PGPLOT
  integer, parameter :: plotlib_extend_none = 0    ! not implemented in PGPLOT

public

character(len=1),parameter :: plot_left_click = 'A'
character(len=1),parameter :: plot_right_click = 'X'
character(len=1),parameter :: plot_middle_click = 'D'
character(len=1),parameter :: plot_shift_click = achar(15)
character(len=1),parameter :: plot_scroll_up = achar(21)
character(len=1),parameter :: plot_scroll_down = achar(4)
character(len=1),parameter :: plot_scroll_left = achar(12)
character(len=1),parameter :: plot_scroll_right = achar(18)

interface plot_qci
   subroutine PGQCI (CI)
     integer,intent(out) ::  CI
   end subroutine PGQCI
end interface

interface plot_qlw
   subroutine PGQLW (LW)
     integer,intent(out) :: LW
   end subroutine PGQLW
end interface

interface plot_qcr
   subroutine PGQCR (CI, CR, CG, CB)
      integer,intent(in) :: CI
      real,intent(out)   :: CR, CG, CB
    end subroutine PGQCR
 end interface

interface plot_qvp
   subroutine PGQVP (UNITS, X1, X2, Y1, Y2)
     integer,intent(in) :: UNITS
     real,intent(out)   :: X1, X2, Y1, Y2
   end subroutine PGQVP
end interface

interface plot_qwin
   subroutine PGQWIN (X1, X2, Y1, Y2)
     real,intent(out) :: X1, X2, Y1, Y2
   end subroutine PGQWIN
end interface

interface plot_ebuf
   subroutine PGEBUF
   end subroutine PGEBUF
end interface

interface plot_bbuf
   subroutine PGBBUF
   end subroutine PGBBUF
end interface

interface plot_qcs
   subroutine PGQCS(UNITS, XCH, YCH)
     integer,intent(in) :: UNITS
     real,intent(out)   :: XCH, YCH
   end subroutine PGQCS
end interface

interface plot_annotate
   subroutine PGMTXT (SIDE, DISP, COORD, FJUST, TEXT)
     character(len=*),intent(in) :: SIDE, TEXT
     real,intent(in)             :: DISP, COORD, FJUST
   end subroutine PGMTXT
end interface

interface plot_sch
   subroutine PGSCH (SIZE)
     real,intent(in) :: SIZE
   end subroutine PGSCH
end interface

interface plot_sci
   subroutine PGSCI (CI)
     integer,intent(in) :: CI
   end subroutine PGSCI
end interface

interface plot_slw
   subroutine PGSLW (LW)
     integer,intent(in) :: lw
   end subroutine PGSLW
   
   module procedure plot_slw_float
end interface

interface plot_page
   subroutine pgpage
   end subroutine pgpage
end interface

interface plot_close
   subroutine PGEND
   end subroutine PGEND
end interface

interface plot_svp
   subroutine pgsvp(XLEFT, XRIGHT, YBOT, YTOP)
     real,intent(in) :: XLEFT, XRIGHT, YBOT, YTOP
   end subroutine pgsvp
end interface

interface plot_swin
   subroutine pgswin(X1, X2, Y1, Y2)
     real,intent(in) :: X1, X2, Y1, Y2
   end subroutine pgswin
end interface

interface plot_wnad
   subroutine pgwnad(X1, X2, Y1, Y2)
     real,intent(in) :: X1, X2, Y1, Y2
   end subroutine pgwnad
end interface

interface plot_qvsz
   subroutine PGQVSZ (UNITS, X1, X2, Y1, Y2)
     integer,intent(in) :: UNITS
     real,intent(out)   :: X1, X2, Y1, Y2
   end subroutine PGQVSZ
end interface

interface plot_line
   subroutine pgline(npts,xline,yline)
     integer, intent(in)               :: npts
     real, intent(in), dimension(npts) :: xline,yline
   end subroutine pgline
end interface

interface plot_box
   subroutine pgbox(XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
     character*(*),intent(in) :: XOPT, YOPT
     real,intent(in)          :: XTICK, YTICK
     integer,intent(in)       :: NXSUB, NYSUB
   end subroutine pgbox
end interface

interface plot_scr
   subroutine PGSCR (CI, CR, CG, CB)
     integer, intent(in) :: CI
     real, intent(in)    :: CR, CG, CB
   end subroutine pgscr
   
   module procedure pgscra
end interface

interface plot_bins
   subroutine PGBIN (NBIN, X, DATA, CENTER)
     integer,intent(in) :: NBIN
     real,intent(in)    :: X(*), DATA(*)
     LOGICAL,intent(in) :: CENTER
   end subroutine PGBIN
end interface

interface plot_imag
   module procedure plot_imag_transparent
end interface

interface plot_qcol
   subroutine pgqcol(icolmin,icolmax)
     integer,intent(out) :: icolmin,icolmax
   end subroutine pgqcol
end interface

interface plot_qcir
   subroutine pgqcir(icolmin,icolmax)
     integer,intent(out) :: icolmin,icolmax
   end subroutine pgqcir
end interface

interface plot_scir
   subroutine pgscir(icilo, icihi)
     integer,intent(in) :: icilo,icihi
   end subroutine pgscir
end interface

interface plot_ctab
   subroutine pgctab(l,r,g,b,nc,contra,bright)
     integer,intent(in) :: nc
     real,intent(in)    :: l(nc),r(nc),g(nc),b(nc),contra,bright
   end subroutine pgctab
end interface

interface plot_qls
   subroutine pgqls(ls)
     integer,intent(out) :: ls
   end subroutine pgqls
end interface

interface plot_qfs
   subroutine pgqfs(fs)
     integer,intent(out) :: fs
   end subroutine pgqfs
end interface

interface plot_sls
   subroutine pgsls(ls)
     integer,intent(in) :: ls
   end subroutine pgsls
end interface

interface plot_sfs
   subroutine pgsfs(fs)
     integer,intent(in) :: fs
   end subroutine pgsfs
end interface

interface plot_rect
   subroutine pgrect(x1,x2,y1,y2)
     real,intent(in) :: x1,x2,y1,y2
   end subroutine pgrect
   
   module procedure plot_rect_rounded
end interface

interface plot_arro
   subroutine pgarro(x1,y1,x2,y2)
     real,intent(in) :: x1,y1,x2,y2
   end subroutine pgarro
end interface

interface plot_circ
   subroutine pgcirc(xcent,ycent,radius)
     real,intent(in) :: xcent,ycent,radius
   end subroutine pgcirc
end interface

interface plot_lcur
  subroutine pglcur (maxpt, npt, x, y)
    implicit none
    integer, intent(in) :: maxpt
    integer, intent(inout) :: npt
    real,    intent(inout) :: x(*), y(*)
  end subroutine pglcur

  module procedure plot_clcur
end interface plot_lcur

interface plot_olin
   subroutine pgolin (maxpt, npt, x, y, symbol)
     implicit none
     integer, intent(in)    :: maxpt
     integer, intent(inout) :: npt
     real,    intent(inout) :: x(*), y(*)
     integer, intent(in)    :: symbol
   end subroutine pgolin
end interface plot_olin

interface plot_ncur
  subroutine pgncur(maxpt, npt, x, y, symbol)
    implicit none
    integer, intent(in)    :: maxpt
    integer, intent(inout) :: npt
    real,    intent(inout) :: x(*), y(*)
    integer, intent(in)    :: symbol
 end subroutine pgncur
end interface plot_ncur

interface plot_qtxt
   subroutine pgqtxt(x,y,angle,fjust,text,xbox,ybox)
     real,intent(in)               :: x, y, angle, fjust
     character(len=*),intent(in)   :: text
     real,intent(out),dimension(4) :: xbox,ybox
   end subroutine pgqtxt
end interface

interface plot_ptxt
   subroutine pgptxt(x,y,angle,fjust,text)
     real,intent(in)             :: x,y,angle,fjust
     character(len=*),intent(in) :: TEXT
   end subroutine pgptxt
end interface

interface plot_stbg
   subroutine pgstbg(bg)
     integer,intent(in) :: bg
   end subroutine pgstbg
end interface

interface plot_curs
 !  integer function pgcurs(x,y,ch)
 !    real,intent(inout)             :: x,y
 !    character*(*),intent(out) :: ch
 !  end function pgcurs
   module procedure pgcurs_sub
end interface

interface plot_pt
   subroutine pgpt(n,xpts,ypts,symbol)
     integer,intent(in) :: n
     real,intent(in)    :: xpts(*),ypts(*)
     integer,intent(in) :: symbol
   end subroutine pgpt
end interface

interface plot_funx
   subroutine pgfunx(fx,n,ymin,ymax,pgflags)
     real,external :: fx
     integer,intent(in) :: n,pgflags
     real,intent(in)    :: ymin,ymax
   end subroutine pgfunx
end interface

interface plot_label
   subroutine pglabel(xlbl,ylbl,toplbl)
     character(len=*),intent(in) :: xlbl,ylbl,toplbl
   end subroutine pglabel
end interface

interface plot_scrn
   subroutine pgscrn(ci,name,ier)
     integer,intent(in)          :: ci
     character(len=*),intent(in) :: name
     integer,intent(out)         :: ier
   end subroutine pgscrn
end interface

interface plot_poly
   subroutine pgpoly(n,xpts,ypts)
     integer,intent(in) :: n
     real,intent(in)    :: xpts(*),ypts(*)
   end subroutine pgpoly
end interface

interface plot_qinf
   subroutine pgqinf(item,value,length)
     character(len=*),intent(in)  :: item
     character(len=*),intent(out) :: value
     integer,intent(out)          :: length
   end subroutine pgqinf
end interface

interface plot_band
   module procedure pgband_sub
end interface

interface plot_pt1
   subroutine pgpt1(xpt,ypt,symbol)
     real,intent(in)     :: xpt,ypt
     integer,intent(in)  :: symbol
   end subroutine pgpt1
end interface

interface plot_numb
   subroutine pgnumb(m,pp,form,string,nc)
     integer,intent(in)           :: m,pp,form
     character(len=*),intent(out) :: string
     integer,intent(out)          :: nc
   end subroutine pgnumb
end interface

interface plot_qch
   subroutine pgqch(ch)
     real,intent(out) :: ch
   end subroutine pgqch
end interface

interface plot_text
   subroutine pgtext(x,y,text)
     real,intent(in)             :: x,y
     character(len=*),intent(in) :: text
   end subroutine pgtext
end interface

interface plot_err1
   subroutine pgerr1(dir,x,y,e,t)
     integer,intent(in) :: dir
     real,intent(in)    :: x,y,e
     real,intent(in)    :: t
   end subroutine pgerr1
end interface

interface plot_errb
   subroutine pgerrb(dir,n,x,y,e,t)
     integer,intent(in) :: dir,n
     real,intent(in)    :: x(n),y(n),e(n)
     real,intent(in)    :: t
   end subroutine pgerrb
end interface

interface plot_conb
   subroutine pgconb(a,idim,jdim,i1,i2,j1,j2,c,nc,tr,blank)
     integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
     real,intent(in)    :: a(idim,jdim),c(*),tr(6),blank
   end subroutine pgconb
end interface

interface plot_cons
   subroutine pgcons(a,idim,jdim,i1,i2,j1,j2,c,nc,tr)
     integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
     real,intent(in)    :: a(idim,jdim),c(*),tr(6)
   end subroutine pgcons
end interface

interface plot_conl
   subroutine pgconl(a,idim,jdim,i1,i2,j1,j2,c,tr,label,intval,mininit)
     integer,intent(in)          :: idim,jdim,i1,i2,j1,j2,intval,mininit
     real,intent(in)             :: a(idim,jdim),c,tr(6)
     character(len=*),intent(in) :: label
   end subroutine pgconl
end interface

interface plot_sah
   subroutine pgsah(fs, angle, cutback)
     integer, intent(in) :: fs
     real, intent(in)    :: angle, cutback
   end subroutine pgsah
end interface

interface plot_vect
   subroutine pgvect(a,b,idim,jdim,i1,i2,j1,j2,c,nc,tr,blank)
     integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
     real,intent(in)    :: a(idim,jdim),b(idim,jdim),tr(6),blank,c
   end subroutine pgvect
end interface

interface plot_pixl
   subroutine pgpixl(ia,idim,jdim,i1,i2,j1,j2,x1,x2,y1,y2)
     integer,intent(in) :: idim,jdim,i1,i2,j1,j2
     integer,intent(in) :: ia(idim,jdim)
     real,intent(in)    :: x1,x2,y1,y2
   end subroutine pgpixl
end interface

interface plot_env
   subroutine pgenv(xmin,xmax,ymin,ymax,just,axis)
     real,intent(in)    :: xmin,xmax,ymin,ymax
     integer,intent(in) :: just,axis
   end subroutine pgenv
end interface

contains

!---------------------------------------------
! initialise the plotting library
!---------------------------------------------
subroutine plot_init(devicein, ierr, papersizex, aspectratio, paperunits)
 implicit none

 character*(*), intent(in)    :: devicein
 integer, intent(out)         :: ierr
 real, intent(in), optional   :: papersizex,aspectratio
 integer, intent(in), optional :: paperunits
 integer                      :: pgopen
 real :: aspect

 if (devicein(1:1).eq.'?') then
    call pgbegin(0,'?',1,1)
    ierr = 1
 else
    ierr = pgopen(devicein)
 endif
 !--check if there is an error
 !  (be careful here: from PGPLOT zero or -ve indicates an error)
 if (ierr.le.0) then
    if (ierr.eq.0) ierr = -1   !--make sure we return an error
    return
 else
    ierr = 0
 endif
 !-- Turn off promting
 call pgask(.false.)

 !-- set paper size if given
 if (present(papersizex)) then
    if (present(aspectratio)) then
       aspect = aspectratio
    else
       aspect = sqrt(2.)
    endif
    if (present(paperunits)) then
    !--make sure that the units are in inches for PGPLOT
       if (paperunits.ne.1) return
    endif
    call plot_pap(papersizex,aspect)
 endif

end subroutine plot_init

subroutine plot_slc(lc)
  implicit none
  integer,intent(in) :: lc

  !--line cap has no effect in PGPLOT

end subroutine plot_slc

subroutine plot_pap(width,aspect,paperunits)
  real,intent(in) :: width,aspect
  integer, intent(in), optional :: paperunits

  if (present(paperunits)) then
     if (paperunits.ne.1) print "(a)",' WARNING: units not valid for PGPLOT'
  endif
  call pgpap(papersizex,aspect)

end subroutine plot_pap

subroutine plot_qlc(lc)
  implicit none
  integer,intent(out) :: lc

  lc = 0
end subroutine plot_qlc

subroutine plot_set_opacity(alpha)
  implicit none
  real, intent(in)           :: alpha

  !--opacity has no effect in PGPLOT

end subroutine plot_set_opacity

!--interface to set transparent colour
!  (not implemented in PGPLOT)
subroutine pgscra (CI, CR, CG, CB, CA)
  integer, intent(in) :: CI
  real, intent(in)    :: CR, CG, CB, CA
  
  !--just throw away the alpha value
  call PGSCR(CI,CR,CG,CB)

end subroutine pgscra

!--floating point line widths
!  (not implemented in PGPLOT)
subroutine plot_slw_float (LW)
  real,intent(in) :: lw
  
  call PGSLW(nint(lw))
  
end subroutine plot_slw_float

subroutine plot_rgb_from_table(frac,r,g,b)
  implicit none
  real, intent(in) :: frac
  real, intent(out) :: r,g,b

  !--rgb from table not implemented in PGPLOT

end subroutine plot_rgb_from_table

logical function plot_qcur()
  implicit none
  character(len=10) :: string
  integer           :: nc
  call pgqinf('CURSOR',string,nc)

  if(string(1:nc).eq.'YES') then
     plot_qcur = .true.
  else
     plot_qcur = .false.
  end if

end function plot_qcur

!
!--inverts the return value of pgcurs
!
function pgcurs_sub(x,y,ch)
  real,intent(inout)        :: x,y
  character*(*),intent(out) :: ch
  integer :: pgcurs_sub,ierr
  integer, external :: pgcurs

  ierr = pgcurs(x,y,ch)

  if (ierr.eq.0) then
    pgcurs_sub = 1
  else
    pgcurs_sub = 0
  endif

end function pgcurs_sub

!--transparent rendering does not work in PGPLOT, but
!  we give it an interface anyway
subroutine plot_imag_transparent(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr, iextend)
 implicit none
 integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
 real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
 integer, intent(in), optional :: iextend

 call pgimag(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr)

end subroutine plot_imag_transparent

subroutine plot_imag_alpha(dat, alpha, idim, jdim, i1, i2, j1, j2, a1, a2, tr, iextend)
  integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
  real,intent(in)    :: dat(IDIM,JDIM), alpha(IDIM,JDIM), A1, A2, TR(6)
  real               :: affine(6)
  integer, intent(in), optional :: iextend

  call pgimag(dat, idim, jdim, i1, i2, j1, j2, a1, a2, tr)

end subroutine plot_imag_alpha

!--giza version of plot_gray takes additional arguments
subroutine plot_gray(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr, iextend)
 implicit none
 integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
 real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
 integer, intent(in), optional :: iextend

 call pggray(a, idim, jdim, i1, i2, j1, j2, a1, a2, tr)

end subroutine plot_gray

!--version of lcur that returns last character pressed
subroutine plot_clcur(maxpt, npt, x, y, ch)
 implicit none
 integer, intent(in) :: maxpt
 integer, intent(inout) :: npt
 real,    intent(inout) :: x(*), y(*)
 character*(*),intent(out) :: ch

 call pglcur (maxpt, npt, x, y)
 ch = 'A'
 
end subroutine plot_clcur

!--rounded rectangle plotting
!  (not implemented -- just calls pgrect)
subroutine plot_rect_rounded(x1,x2,y1,y2,r)
 implicit none
 real,intent(in) :: x1,x2,y1,y2,r

 call pgrect(x1,x2,y1,y2)

end subroutine plot_rect_rounded
!
!--inverts the return value of pgband
!
function pgband_sub(mode, posn, xref, yref, x, y, ch)
  integer,intent(in) :: mode, posn
  real,intent(in)    :: xref, yref
  real,intent(inout) :: x, y
  character*(*),intent(out) :: ch
  integer :: ierr,pgband_sub
  integer,external :: pgband

  ierr = pgband(mode,posn,xref,yref,x,y,ch)

  if(ierr.eq.1) then
    pgband_sub = 0
  else
    pgband_sub = 1
  endif
end function pgband_sub

!
!--this subroutine can be called after PGSVP to
!  make sure that the viewport lies exactly on
!  pixel boundaries.
!
!  Queries PGPLOT routines directly so no need
!  for input/output
!

subroutine plot_set_exactpixelboundaries()
 implicit none
 real :: xminpix,xmaxpix,yminpix,ymaxpix
 real :: vptxmin,vptxmax,vptymin,vptymax
 real :: dv
 real, parameter :: tol = 1.e-6
 !
 ! setting axes adjusts the viewport, so query to get adjusted settings
 !
 call pgqvp(0,vptxmin,vptxmax,vptymin,vptymax)
 !print*,'got ',vptxmin,vptxmax,vptymin,vptymax
 !
 ! adjust viewport on pixel devices so that
 ! boundaries lie exactly on pixel boundaries
 !
 !  query viewport size in pixels
 call pgqvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
 !print*,' in pixels = ',xminpix,xmaxpix,yminpix,ymaxpix

 !  work out how many viewport coords/pixel
 dv = (vptymax - vptymin)/(ymaxpix-yminpix)

 !  adjust viewport min/max to lie on pixel boundaries
 vptymin = max((nint(yminpix)-tol)*dv,0.)
 vptymax = min((nint(ymaxpix)-tol)*dv,1.0-epsilon(1.0)) ! be careful of round-off errors

 !  same for x
 dv = (vptxmax - vptxmin)/(xmaxpix-xminpix)
 vptxmin = max((nint(xminpix)-tol)*dv,0.)
 vptxmax = min((nint(xmaxpix)-tol)*dv,1.0-epsilon(1.0)) ! be careful of round-off errors

 !  adjust viewport
 !print*,'adjusting ',vptxmin,vptxmax,vptymin,vptymax
 call pgsvp(vptxmin,vptxmax,vptymin,vptymax)

 !call pgqvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
 !print*,' in pixels = ',xminpix,xmaxpix,yminpix,ymaxpix

 return
end subroutine plot_set_exactpixelboundaries



end module plotlib
