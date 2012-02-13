      INTERFACE
      subroutine cd2cs3_f90(n,ia,la,ja,aa,ib,jb,ab)
      integer n
      integer, pointer :: ia(:),ja(:),la(:)
      integer ib(*),jb(*)
      double precision, pointer :: aa(:)
      double precision ab(*)
      END
      END INTERFACE
