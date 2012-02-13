      INTERFACE
      subroutine cd2cs_f90(n,ia,da,ja,la,aa,ib,jb,ab,diag_one)
      integer n,diag_one
      integer, pointer :: ia(:),ja(:),la(:)
      integer ib(*),jb(*)
      double precision, pointer :: da(:),aa(:)
      double precision ab(*)
      END
      END INTERFACE
