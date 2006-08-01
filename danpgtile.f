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
c   just   : just=1 gives equal aspect ratios (same as in PGENV)
c   axis   : axes options (same as in PGENV)
c   vmarginleft,right,bottom,top : initial margin sizes (% of page)
c            (default should be zero for these)
c
c  Daniel Price, Institute of Astronomy, Cambridge, 2004.
c
      SUBROUTINE DANPGTILE(iplotin,nx,ny,xmin,xmax,ymin,ymax, 
     &                     labelx,labely,title,just,axis,
     &                     vmarginleftin,vmarginrightin,
     &                     vmarginbottomin,vmargintopin) 
      IMPLICIT NONE
      INTEGER iplotin,nx,ny,just,axis
      INTEGER iplot,ix,iy
      REAL xmin,xmax,ymin,ymax,vmarginleftin,vmarginrightin
      REAL vmargintopin,vmarginbottomin
      REAL vptsizeeffx,vptsizeeffy,panelsizex,panelsizey
      REAL vmargintop,vmarginbottom,vmarginleft,vmarginright
      REAL vptxmin,vptxmax,vptymin,vptymax
      REAL aspectratio,devaspectratio,x1,x2,y1,y2
      REAL xch,ych,xlabeloffset,ylabeloffset
      CHARACTER xopts*10, yopts*10
      CHARACTER*(*) labelx,labely,title
c
c new page if iplot > number of plots on page
c
      IF (iplotin.GT.nx*ny) THEN
         IF (MOD(iplotin,nx*ny).EQ.1) CALL PGPAGE
         iplot = iplotin - (nx*ny)*((iplotin-1)/(nx*ny))
      ELSEIF (iplotin.LE.0) THEN
         RETURN
      ELSE
         iplot = iplotin
      ENDIF
c
c check for errors in input
c      
      IF (nx.LE.0 .OR. ny.LE.0) RETURN
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
c set positions of x and y labels in units of character height from edge
c
      xlabeloffset = 3.0
      ylabeloffset = 4.5
c
c query the character height as fraction of viewport
c
      CALL PGQCS(0,xch,ych)
c
c set margin size in units of viewport dimensions
c allow enough room for the plot labels if they are drawn
c NB: PGPLOT sets the character height as some fraction of the smallest
c     dimension
c
      vmargintop = vmargintopin
      vmarginright = vmarginrightin
      IF (axis.GE.0) THEN
         vmarginleft = vmarginleftin + (ylabeloffset+1.0)*xch
         vmarginbottom = vmarginbottomin + (xlabeloffset+1.0)*ych
      ELSE
         vmarginleft = vmarginleftin
         vmarginbottom = vmarginbottomin
      ENDIF
c
c effective viewport size = size - margins
c
      vptsizeeffx = 1.0 - vmarginright - vmarginleft
      vptsizeeffy = 1.0 - vmargintop - vmarginbottom
c     reduce x or y size if just=1 to get right aspect ratio
      IF (aspectratio.lt.1.0) THEN
         vptsizeeffx = aspectratio*vptsizeeffy
      ELSEIF (aspectratio.gt.1.0) THEN
         vptsizeeffy = vptsizeeffx/aspectratio
      ENDIF
c
c--set size of each panel
c      
      panelsizex = vptsizeeffx/nx
      panelsizey = vptsizeeffy/ny 
      ix = iplot - ((iplot-1)/nx)*nx
      iy = (iplot-1)/nx + 1
c      print*,i,ix,iy
c      print*,panelsizex,panelsizey,vptsizeeffx,vptsizeeffy
      
      vptxmin = vmarginleft + (ix-1)*panelsizex
      vptxmax = vptxmin + panelsizex
      vptymax = (1.0 - vmargintop) - (iy-1)*panelsizey
      vptymin = vptymax - panelsizey
c      print*,vptxmin,vptxmax,vptymin,vptymax
      CALL PGSVP(vptxmin,vptxmax,vptymin,vptymax)
c
c set axes
c
      IF (just.EQ.1) THEN
         CALL PGWNAD(xmin,xmax,ymin,ymax)
      ELSE
         CALL PGSWIN(xmin,xmax,ymin,ymax)
      ENDIF
c
c option to return before actually doing anything
c      
      IF (title.EQ.'NOPGBOX') RETURN
c
c set options for call to pgbox (draws axes) and label axes where appropriate
c (options are exactly as in PGENV apart from axis=-3 which I have added)
c
      yopts = '*'
      IF (axis.EQ.-3) THEN
         xopts = 'BCST'
      ELSEIF (axis.EQ.-2) THEN
         xopts = ' '
      ELSEIF (axis.EQ.-1) THEN
        xopts = 'BC'
      ELSEIF (axis.EQ.0) THEN
        xopts = 'BCST'
      ELSEIF (axis.EQ.1) THEN
        xopts = 'ABCST'
      ELSEIF (axis.EQ.2) THEN
        xopts = 'ABCGST'
      ELSEIF (axis.EQ.10) THEN
        xopts = 'BCSTL'
        yopts = 'BCST'
      ELSEIF (axis.EQ.20) THEN
        xopts = 'BCST'
        yopts = 'BCSTL'
      ELSEIF (axis.EQ.30) THEN
        xopts = 'BCSTL'
        yopts = 'BCSTL'
      ELSE
        CALL GRWARN('DANPGTILE: illegal AXIS argument.')
        xopts = 'BCNST'
      ENDIF
      IF (yopts.eq.'*') yopts = xopts
c
c decide whether to number and label the y axis
c      
      IF (ix.EQ.1 .AND. axis.GE.0) THEN
         yopts = '1VN'//yopts
         CALL PGMTXT('L',ylabeloffset,0.5,0.5,labely)
      ELSEIF (axis.GE.0) THEN
         yopts = yopts//'N'
      ENDIF  
c
c decide whether to number and label the x axis
c      
      IF (iy.EQ.ny .AND. axis.GE.0) THEN
         xopts = 'N'//xopts
         CALL PGMTXT('B',xlabeloffset,0.5,0.5,labelx)
      ENDIF
      
      CALL PGBOX(xopts,0.0,0,yopts,0.0,0)
c
c plot the title inside the plot boundaries
c
      CALL PGMTXT('T',-1.5,0.96,1.0,title)
      
      RETURN      
      END SUBROUTINE
