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
