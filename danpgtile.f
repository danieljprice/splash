c
c  simple subroutine to tile graphs appropriately on a page in PGPLOT
c  divides up a single panel into subpanels, with a margin at the edge
c  should replace the call to PGENV and PGLABEL
c
c  the page setup looks like this:
c
c    |   |   |   | 
c  --+---+---+---+--
c    | 1 | 2 | 3 |
c  --+---+---+---+--
c    | 4 | 5 | 6 |
c  --+---+---+---+--
c    |   |   |   |
c
c  (ie. with margins in x and y)
c  Note that we divide up a single panel, so PGBEG should be called with nx=1,ny=1
c
c  Arguments:
c   iplot  : current plot number
c   nx     : number of panels in x direction
c   ny     : number of panels in y direction
c   xmin,  : xmax, ymin, ymax : plot limits (should be same for all plots)
c   labelx : x axis label (should be same for all plots)
c   labely : y axis label (should be same for all plots)
c   title  : current plot title (can differ between plots)
c
c  Daniel Price, Institute of Astronomy, Cambridge, 2004.
c
      SUBROUTINE DANPGTILE(iplot,nx,ny,xmin,xmax,ymin,ymax, 
     &                      labelx,labely,title,just)
      IMPLICIT NONE
      INTEGER iplot,nx,ny,ix,iy,just
      REAL xmin,xmax,ymin,ymax
      REAL vptsizeeffx,vptsizeeffy,panelsizex,panelsizey
      REAL vmargintop,vmarginbottom,vmarginleft,vmarginright
      REAL vptxmin,vptxmax,vptymin,vptymax
      REAL aspectratio,devaspectratio,x1,x2,y1,y2
      CHARACTER optsx*7, optsy*7
      CHARACTER*(*) labelx,labely,title
c
c new page if iplot > number of plots on page
c
      IF (iplot.GT.nx*ny) THEN
         IF (MOD(iplot,nx*ny).EQ.1) CALL PGPAGE
	 iplot = iplot - (nx*ny)*((iplot-1)/(nx*ny))
      ELSEIF (iplot.LE.0) THEN
         RETURN
      ENDIF
c
c check for errors in input
c      
      IF (nx.LE.0 .OR. ny.LE.0) RETURN
c
c set margin size in units of viewport dimensions
c
      vmarginleft = 0.06
      vmarginright = 0.001
      vmargintop = 0.001
      vmarginbottom = 0.08
c
c adjust effective viewport size if just=1 and graphs are not square
c      
      IF (just.eq.1) THEN
         IF (ymax.EQ.ymin) THEN
	    PRINT*,'DANPGTILE: Error: ymax=ymin'
	    RETURN
	 ENDIF
c
c query the current aspect ratio of the device and set aspect ratio appropriately
c
         CALL PGQVSZ(3,x1,x2,y1,y2)
	 devaspectratio = (x2-x1)/(y2-y1)
         aspectratio = ((xmax-xmin)*nx)/((ymax-ymin)*ny)/devaspectratio
      ELSE
	 aspectratio = 1.0
      ENDIF
c
c effective viewport size = size - margins
c
      vptsizeeffx = 1.0 - vmarginright - vmarginleft
      vptsizeeffy = 1.0 - vmargintop - vmarginbottom
      vptsizeeffx = aspectratio*vptsizeeffy   !!/devaspectratio
      panelsizex = vptsizeeffx/nx
      panelsizey = vptsizeeffy/ny 
      ix = iplot - ((iplot-1)/nx)*nx
      iy = (iplot-1)/nx + 1
!!      print*,i,ix,iy
!!      print*,panelsizex,panelsizey,vptsizeeffx,vptsizeeffy
      
      vptxmin = vmarginleft + (ix-1)*panelsizex
      vptxmax = vptxmin + panelsizex
      vptymin = vmarginbottom + (ny-iy)*panelsizey
      vptymax = vptymin + panelsizey
!!      print*,vptxmin,vptxmax,vptymin,vptymax
      CALL PGSVP(vptxmin,vptxmax,vptymin,vptymax)
c
c--set axes
c
      IF (just.EQ.1) THEN
         CALL PGWNAD(xmin,xmax,ymin,ymax)
      ELSE
         CALL PGSWIN(xmin,xmax,ymin,ymax)
      ENDIF
c
c set options for call to pgbox (draws axes) and label axes where appropriate
c
      IF (ix.EQ.1) THEN
         optsy = '1bvcnst'
	 call pgmtxt('L',2.5,0.5,0.5,labely)
      ELSE
         optsy = 'bcst'
      ENDIF  
      IF (iy.EQ.ny) THEN
         optsx = 'bcnsta' 
	 call pgmtxt('B',3.0,0.5,0.5,labelx)
      ELSE
         optsx = 'bcsta'
      ENDIF
      CALL PGBOX(optsx,0.0,0,optsy,0.0,0)
c
c plot the title inside the plot boundaries
c      
      call pgmtxt('T',-1.5,0.96,1.0,title)
      
      RETURN      
      END SUBROUTINE
