      INTERFACE
      subroutine ainvsr2(msglvl,msgunit,n,ia,ja,a,ip,jp,ap,
     *  size_p,size_c,size_r,diagtol,
     *  drfl,mi,diag_one,droptyp,imodif,fill,fillmax,
     *  ifillmax,garrow,garcol,info)
      integer msglvl,msgunit,n,size_r,size_c,size_p
      integer ia(*),ja(*)
      double precision a(*)
      integer, pointer :: ip(:),jp(:)
      double precision, pointer :: ap(:)
      integer garcol,garrow,droptyp
      integer imodif,diag_one,fill,fillmax,ifillmax,info
      double precision mi,drfl,diagtol
      END
      END INTERFACE

      INTERFACE
      subroutine ainvini1(n,m,ia,ja,ic,jc,kc,ac,max_c,ir,jr,kr,max_r,
     *  iendru,iendlu,idistr,idistc,lsize,ierr)
      integer n,m,max_r,max_c,lsize,ierr
      integer idistr,idistc,len
      integer ia(*),ja(*)
      integer, pointer :: ir(:),jr(:),kr(:),ic(:),jc(:),kc(:)
      double precision, pointer :: ac(:)
      integer iendru,iendlu
      END
      END INTERFACE

      INTERFACE
      subroutine rank1n(scrlvl,n,m,ptr,cnr,lenr,
     *  max_r,ptc,cnc,lenc,h,max_c,
     *  indexk,wr01,wn03,wr03,ind3,wn02,wr02,ind2,
     *  iendru,iendlu,drfl,nit,if2,idist,garcol,garrow,droptyp,
     *  nreallocr,nreallocc)
      integer n,m,droptyp,nreallocr,nreallocc,scrlvl
      integer nit,if2,max_r,max_c
      integer, pointer :: ptr(:),cnr(:),lenr(:),ptc(:),cnc(:),lenc(:)
      double precision, pointer :: h(:)
      integer iendru,iendlu,nrow,garrow,garcol,indexk,ind2,ind3
      integer idist
      integer wn02(*),wn03(*)
      double precision wr01(*),wr02(*),wr03(*),temp
      double precision drfl
      END
      END INTERFACE

      INTERFACE
      subroutine rank1n2(scrlvl,n,m,ptr,cnr,lenr,
     *  max_r,ptc,cnc,lenc,h,max_c,
     *  indexk,wr01,wn03,wr03,ind3,wn02,wr02,ind2,
     *  iendru,iendlu,drfl,nit,if2,idist,garcol,garrow,droptyp,
     *  nreallocr,nreallocc)
      integer scrlvl,n,m,droptyp,nreallocr,nreallocc
      integer nit,if2,max_r,max_c
      integer, pointer :: ptr(:),cnr(:),lenr(:),ptc(:),cnc(:),lenc(:)
      double precision, pointer :: h(:)
      integer iendru,iendlu,nrow,garrow,garcol,indexk,ind2,ind3
      integer idist
      integer wn02(*),wn03(*)
      double precision wr01(*),wr02(*),wr03(*),temp
      double precision drfl
      END
      END INTERFACE

      INTERFACE
      subroutine cd2cs_f90(n,ia,da,ja,la,aa,ib,jb,ab,diag_one)
      integer n,diag_one
      integer, pointer :: ia(:),ja(:),la(:)
      integer ib(*),jb(*)
      double precision, pointer :: da(:),aa(:)
      double precision ab(*)
      END
      END INTERFACE

      INTERFACE
      subroutine cd2cs3_f90(n,ia,la,ja,aa,ib,jb,ab)
      integer n
      integer, pointer :: ia(:),ja(:),la(:)
      integer ib(*),jb(*)
      double precision, pointer :: aa(:)
      double precision ab(*)
      END
      END INTERFACE

      INTERFACE
      subroutine srtcs1_f90(n,ia,ja,a)
      integer n
      integer :: ia(*),ja(*)
      double precision :: a(*)
      END
      END INTERFACE

      INTERFACE
      subroutine realloci(ja,arrsize,arrincr,new_arrsize,ierr)
      integer arrsize,arrincr,new_arrsize,ierr
      integer, pointer, dimension(:) :: ja
      END
      END INTERFACE

      INTERFACE
      subroutine reallocr(aa,arrsize,arrincr,new_arrsize,ierr)
      integer arrsize,arrincr,new_arrsize,ierr
      double precision, pointer, dimension(:) :: aa
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE MVUCSR1_F90(M,IA,JA,A,X,Y)
      INTEGER M
      INTEGER, POINTER :: IA(:),JA(:)
      DOUBLE PRECISION, POINTER :: A(:)
      DOUBLE PRECISION X(*),Y(*)
      END
      END INTERFACE

      INTERFACE
      SUBROUTINE TVUCSR1_F90(M,N,IA,JA,A,X,Y)
      INTEGER M,N
      INTEGER, POINTER :: IA(:),JA(:)
      DOUBLE PRECISION, POINTER :: A(:)
      DOUBLE PRECISION X(*),Y(*)
      END
      END INTERFACE

