!
!     writes default options to file (should match read_defaults)
!
      SUBROUTINE write_defaults
       USE exact_params
       USE settings
       USE multiplot
       IMPLICIT NONE
       
       OPEN(unit=1,file='defaults',status='replace',form='formatted')
          WRITE(1,*) animate,magfield,iadapt,xsec_nomulti,
     &	             flythru,plotcirc
          WRITE(1,*) imark,imarkg,nacross,ndown,nyplotmulti
          WRITE(1,*) iplotline,iplotlinein,linestylein
          WRITE(1,*) iexact,iplotav,nbins
          WRITE(1,*) irender,ivecplot_nomulti,iplotpartvec_nomulti,
     &	             npix_nomulti,npixvec_nomulti
          WRITE(1,*) iplotcont_nomulti,ncontours_nomulti,
     &	             icolours,iplotghost
          WRITE(1,*) ipapersize,papersizex,aspectratio
          WRITE(1,*) delta,lambda
          WRITE(1,*) Htstar,Atstar,Ctstar,sigma0,norder
          
          WRITE(1,*) itrans(:)
          WRITE(1,*) multiplotx(1:nyplotmulti)
          WRITE(1,*) multiploty(1:nyplotmulti)
          WRITE(1,*) irendermulti(1:nyplotmulti)
          WRITE(1,*) iplotcontmulti(1:nyplotmulti)
          WRITE(1,*) ncontoursmulti(1:nyplotmulti)
          WRITE(1,*) ivecplotmulti(1:nyplotmulti)
          WRITE(1,*) npixmulti(1:nyplotmulti)
          WRITE(1,*) npixvecmulti(1:nyplotmulti)
          WRITE(1,*) iplotpartvecmulti(1:nyplotmulti)    
	  WRITE(1,*) x_secmulti(1:nyplotmulti)
	  WRITE(1,*) xsecposmulti(1:nyplotmulti)
      CLOSE(unit=1)
      PRINT*,'default options saved to file'
    
      RETURN              
      END SUBROUTINE write_defaults
