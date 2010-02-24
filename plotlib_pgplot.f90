!---------------------------------------------------------------------------
! module containing application programming interfaces for basic
! plotting functions. The idea is to add more to this module to
! eventually use it to be able to change backends more easily.
!---------------------------------------------------------------------------
module plotlib 
  implicit none
  
public

interface plot_qci
   SUBROUTINE PGQCI (CI)
     INTEGER,intent(out) ::  CI
   end SUBROUTINE PGQCI
end interface

interface plot_qlw
   SUBROUTINE PGQLW (LW)
     INTEGER,intent(out) :: LW
   end SUBROUTINE PGQLW
end interface

interface plot_qcr
   SUBROUTINE PGQCR (CI, CR, CG, CB)
      INTEGER,intent(in) :: CI
      REAL,intent(out)   :: CR, CG, CB
    end SUBROUTINE PGQCR
 end interface

interface plot_qvp
   SUBROUTINE PGQVP (UNITS, X1, X2, Y1, Y2)
     INTEGER,intent(in) :: UNITS
     REAL,intent(out)   :: X1, X2, Y1, Y2
   end SUBROUTINE PGQVP
end interface

interface plot_qwin
   SUBROUTINE PGQWIN (X1, X2, Y1, Y2)
     REAL,intent(out) :: X1, X2, Y1, Y2
   end SUBROUTINE PGQWIN
end interface

interface plot_ebuf
   SUBROUTINE PGEBUF
   end SUBROUTINE PGEBUF
end interface

interface plot_bbuf
   SUBROUTINE PGBBUF
   end SUBROUTINE PGBBUF
end interface

interface plot_qcs
   SUBROUTINE PGQCS(UNITS, XCH, YCH)
     INTEGER,intent(in) :: UNITS
     REAL,intent(out)   :: XCH, YCH
   end SUBROUTINE PGQCS
end interface

interface plot_annotate
   SUBROUTINE PGMTXT (SIDE, DISP, COORD, FJUST, TEXT)
     CHARACTER(len=*),intent(in) :: SIDE, TEXT
     REAL,intent(in)             :: DISP, COORD, FJUST
   end SUBROUTINE PGMTXT
end interface

interface plot_sch
   SUBROUTINE PGSCH (SIZE)
     REAL,intent(in) :: SIZE
   end SUBROUTINE PGSCH
end interface

interface plot_sci
   SUBROUTINE PGSCI (CI)
     INTEGER,intent(in) :: CI
   end SUBROUTINE PGSCI
end interface

interface plot_slw
   SUBROUTINE PGSLW (LW)
     INTEGER,intent(in) :: lw
   end SUBROUTINE PGSLW
end interface

interface plot_page
   subroutine pgpage
   end subroutine pgpage
end interface

interface plot_close
   SUBROUTINE PGEND
   end SUBROUTINE PGEND
end interface

interface plot_svp
   subroutine pgsvp(XLEFT, XRIGHT, YBOT, YTOP)
     REAL,intent(in) :: XLEFT, XRIGHT, YBOT, YTOP
   end subroutine pgsvp
end interface

interface plot_swin
   subroutine pgswin(X1, X2, Y1, Y2)
     REAL,intent(in) :: X1, X2, Y1, Y2
   end subroutine pgswin
end interface

interface plot_qvsz
   subroutine PGQVSZ (UNITS, X1, X2, Y1, Y2)
     INTEGER,intent(in) :: UNITS
     REAL,intent(out)   :: X1, X2, Y1, Y2
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
     CHARACTER*(*),intent(in) :: XOPT, YOPT
     REAL,intent(in)          :: XTICK, YTICK
     INTEGER,intent(in)       :: NXSUB, NYSUB
   end subroutine pgbox
end interface

interface plot_scr
   SUBROUTINE PGSCR (CI, CR, CG, CB)
     INTEGER, intent(in) :: CI
     REAL, intent(in)    :: CR, CG, CB
   end subroutine pgscr
end interface

interface plot_bins
   SUBROUTINE PGBIN (NBIN, X, DATA, CENTER)
     INTEGER,intent(in) :: NBIN
     REAL,intent(in)    :: X(*), DATA(*)
     LOGICAL,intent(in) :: CENTER
   end SUBROUTINE PGBIN
end interface

interface plot_imag
   subroutine pgimag(a, idim, jdim, i1, i2, j1, j2,&
        a1, a2, tr)
     integer,intent(in) :: IDIM, JDIM, I1, I2, J1, J2
     real,intent(in)    :: A(IDIM,JDIM), A1, A2, TR(6)
   end subroutine pgimag
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
   subroutine pgcurs(x,y,ch)
     real,intent(out)             :: x,y 
     character(len=*),intent(out) :: ch
   end subroutine pgcurs
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
     real,external      :: fx
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
   subroutine pgband(mode,posn,xref,yref,x,y,ch)
     integer,intent(in) :: mode,posn
     real,intent(in)    :: xref,yref
     real,intent(inout) :: x,y
     character(len=*),intent(out) :: ch
   end subroutine pgband
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
   subroutine pgconb(a,idim,jdim,i1,i2,j1,j2,c,nc,tr, &
     blank)
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
   subroutine pgconl(a,idim,jdim,i1,i2,j1,j2,c,tr, &
        label,intval,mininit)
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
   subroutine pgvect(a,b,idim,jdim,i1,i2,j1,j2,c,nc,tr, &
        blank)
     integer,intent(in) :: idim,jdim,i1,i2,j1,j2,nc
     real,intent(in)    :: a(idim,jdim),b(idim,jdim),tr(6),blank,c
   end subroutine pgvect
end interface

interface plot_pixl
   subroutine pgpixl(ia,idim,jdim,i1,i2,j1,j2, &
        x1,x2,y1,y2)
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

interface plot_pap
   subroutine pgpap(width,aspect)
     real,intent(in) :: width,aspect
   end subroutine pgpap
end interface

contains

subroutine plot_init(devicein, ierr)
 implicit none

 character*(*), intent(in)    :: devicein
 integer, intent(out)         :: ierr
 integer                      :: pgopen

 if (devicein(1:1).eq.'?') then
    call pgbegin(0,'?',1,1)
    ierr = 1
 else
    ierr = pgopen(devicein)
 endif
 !--check if there is an error 
 if (ierr.le.0) return  ! zero or negative indicates an error

! !-- Check if the device is interactive  
! call pgqinf('CURSOR',string,ilen)

! if (string(1:ilen).eq.'YES') then
!  plot_deviceisinteractive = .true.
! else
!  plot_deviceisinteractive = .false.
! endif

 !-- Turn off promting
 call pgask(.false.)

 !-- Check if it is a vector device
!  call pgqinf('TYPE',string,ilen)
!  select case(string(1:ilen))
!  case('PS','CPS','VPS','VCPS','NULL','LATEX')
!     plot_deviceisvector = .true.
!  case default
!     plot_deviceisvector = .false.
!  end select

end subroutine plot_init

subroutine plot_slc(lw)
  implicit none
  integer,intent(in) :: lw
end subroutine plot_slc

subroutine plot_qlc(lw)
  implicit none
  integer,intent(out) :: lw

  lw = 0
end subroutine plot_qlc

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
