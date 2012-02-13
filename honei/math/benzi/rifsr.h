      INTERFACE
      subroutine rifsr(msglvl,msgunit,n,ia,ja,a,size_c,
     *  size_r,ip,jp,ap,size_p,
     *  diagtol,drflic,drflai,mi,diag_one,droptyp,imodif,
     *  filldyn,fillmax,
     *  ifillmax,garrow,garcol,fillic,info)
      integer n,msglvl,msgunit,size_r,size_c,size_p
      integer garcol,garrow,droptyp,fillic
      integer imodif,diag_one,filldyn,fillmax,ifillmax
      integer info
      integer, pointer :: ip(:),jp(:)
      double precision, pointer :: ap(:)
      integer ia(*),ja(*)
      double precision a(*)
      double precision mi,drflai,drflic,diagtol
      END
      END INTERFACE
