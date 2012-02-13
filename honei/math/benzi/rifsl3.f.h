      INTERFACE
      subroutine rifsl3(msglvl,msgunit,n,ia,ja,aa,max_z,
     *  ip,jp,ap,max_p,diag,
     *  diagtol,drfl,drflai,imodif,diag_one,droptyp,fillai,fillic,info)
      integer msglvl,msgunit,n
      integer ia(*),ja(*)
      double precision aa(*)
      integer, pointer :: ip(:),jp(:)
      double precision, pointer :: ap(:),diag(:)
      integer max_z,max_p,droptyp,diag_one
      integer imodif,fillic,info,fillai
      double precision diagtol,drfl,drflai
      END
      END INTERFACE
