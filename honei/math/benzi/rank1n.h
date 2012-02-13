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
