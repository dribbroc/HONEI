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
