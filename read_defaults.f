!
! reads default options from file
!
      SUBROUTINE read_defaults
      USE multiplot
      USE settings
      USE exact_params
      IMPLICIT NONE
      
      INQUIRE (exist=iexist, file='defaults')
      IF (iexist) THEN
         OPEN(unit=1,file='defaults',status='old',form='formatted')
	    READ(1,*,END=7,ERR=8) animate,magfield,iadapt,xsec_nomulti,
     &                             flythru,plotcirc
	    READ(1,*,END=7,ERR=8) imark, imarkg, nacross, ndown,nyplotmulti     
	    READ(1,*,END=7,ERR=8) iplotline,iplotlinein,linestylein
	    READ(1,*,END=7,ERR=8) iexact, iplotav, nbins
	    READ(1,*,END=7,ERR=8) irender,ivecplot_nomulti,
     &	    iplotpartvec_nomulti,npix_nomulti,npixvec_nomulti
	    READ(1,*,END=7,ERR=8) iplotcont_nomulti,ncontours_nomulti,
     &	    icolours,iplotghost
            READ(1,*,END=7,ERR=8) ipapersize,papersizex,aspectratio
	    READ(1,*,END=7,ERR=8) delta,lambda
	    READ(1,*,END=7,ERR=8) Htstar,Atstar,Ctstar,sigma0,norder

            READ(1,*,END=7,ERR=8) itrans(:)
	    READ(1,*,END=7,ERR=8) multiplotx(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) multiploty(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) irendermulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) iplotcontmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) ncontoursmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) ivecplotmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) npixmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) npixvecmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) iplotpartvecmulti(1:nyplotmulti)
	    READ(1,*,END=7,ERR=8) x_secmulti(1:nyplotmulti)	    
	    READ(1,*,END=7,ERR=8) xsecposmulti(1:nyplotmulti)
	 CLOSE(unit=1)
	 PRINT*,'read default options from file '
      ENDIF
      GOTO 9
7     CONTINUE
      PRINT*,'**** Warning: end of file in defaults ****'
      CLOSE(unit=1)
      GOTO 9
8     CONTINUE
      PRINT*,'Error reading defaults from file'
      CLOSE(unit=1)     
9     CONTINUE

      RETURN
      END SUBROUTINE read_defaults
