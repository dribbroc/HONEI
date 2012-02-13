
* (c) sparslab module name=rifsr.main
*
* purpose:
*  a sample driver for computation of symmetric
*  right-looking rif based on
*  symmetric sainv preconditioner construction.
*  this driver does the following:
*    1) takes a simple matrix a stored in csr scheme
*    2) compute the robust incomplete factorization (rif)
*       m approximately equal to a
*    3) applies inverse of m to a given vector x of all ones
*
* history:
*   original version for sparslab - tu - 16/09/2000.
*
      integer n,i
      double precision zero
      parameter(zero=0.0d0)
      integer msglvl,msgunit,size_r,size_c,size_p
      integer, pointer :: ip(:),jp(:)
      double precision, pointer :: ap(:)
      integer garcol,garrow,droptyp
      integer imodif,diag_one,filldyn,fillmax,ifillmax,info
      integer fillic,fill
      double precision mi,drflai,diagtol
      integer, allocatable, dimension(:) :: ia,ja
      double precision, allocatable, dimension(:) :: aa,x,y,wr01
c
      include 'sparslab.i'
c
c start of rifsr.main
c
c  -- set up the initial matrix in csr format
c  -- it is the matrix poisson(2)
c
c
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  -- specify the matrix in csr format
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      n=4
      allocate(ia(n+1))
      ia(1)=1
      ia(2)=4
      ia(3)=7
      ia(4)=10
      ia(5)=13
      allocate(ja(ia(n+1)-1),aa(ia(n+1)-1))
      ja(1)=1
      ja(2)=2
      ja(3)=3
      ja(4)=1
      ja(5)=2
      ja(6)=4
      ja(7)=1
      ja(8)=3
      ja(9)=4
      ja(10)=2
      ja(11)=3
      ja(12)=4
      aa(1)=4
      aa(2)=-1
      aa(3)=-1
      aa(4)=-1
      aa(5)=4
      aa(6)=-1
      aa(7)=-1
      aa(8)=4
      aa(9)=-1
      aa(10)=-1
      aa(11)=-1
      aa(12)=4
c
c     %%%%%%%%%%%%%%%%%
c  -- to be set by user
c     %%%%%%%%%%%%%%%%%
c
c
c  -- set the drop tolerance
c  -- for the sainv processs
c
      drflai=0.1d0
c
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  -- end of the part to be set by user
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  -- start of defaults -- need not to be changed
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  -- sizes: an initial first guess
c  --        will be automatically adapted inside the code
c
      size_p=1000000
      size_r=1000000
      size_c=1000000
c
c  -- other parameters without the need to be changed
c
      diagtol=1.1d-16
      droptyp=0
      mi=0.1d0
      diag_one=1
c
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  -- end of defaults -- need not to be changed
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c  -- standard level of output
c
      msglvl=2
c
c  -- open unit for output
c
      msgunit=2
      open(unit=msgunit,file='sainvsr.main.out')
c
c  -- call sainv
c
      call ainvsr2(msglvl,msgunit,n,
     *  ia,ja,aa,ip,jp,ap,
     *  size_p,size_c,size_r,diagtol,
     *  drflai,mi,diag_one,droptyp,imodif,fill,fillmax,
     *  ifillmax,garrow,garcol,info)
c
c  -- output preconditioner info
c
      write(2,*)' preconditioner successfully computed'
      write(2,*)' preconditioner size',fillic
      write(2,*)' pointers'
      write(2,*) ip(1:n+1)
      write(2,*)' column indices'
      write(2,*) jp(1:ip(n+1)-1)
      write(2,*)' nonzero entries'
      write(2,*) ap(1:ip(n+1)-1)
c
c  -- get right-hand side computed by multiplying
c  -- vector x of all ones
c
      allocate(x(n),y(n),wr01(n))
      do i=1,n
        x(i)=1.0d0
      end do
      call mvucsr1(n,ia,ja,aa,x,y)
      deallocate(ia,ja,aa)
c
c  -- apply inverse of m (computed robust incomplete factorization) to
c  -- a vector y computed from a given vector x of all ones
c
c  -- note that we skip headers of f90 pointer arrays
c  -- using a dirty trick
c
      call mvucsr1(n,ip(1),jp(1),ap(1),y,wr01)
      call tvucsr1(n,n,ip(1),jp(1),ap(1),wr01,x)
c
c  -- we have the result in x
c  -- print solution x
c
      write(2,*)' solution x'
      write(2,*) x(1:n)
c
      deallocate(x,y,wr01)
c
      close(2)
c
c end of rifsr.main
c
      end


* (c) sparslab module name=ainvsr2
*
* purpose:
*  symmetric right-looking ainv stabilized preconditioner construction.
*  factorized approximate inverse is in the form: a^{-1}=z^t*d^{-1}*z,
*  where z is the lower triangular matrix stored by rows.
*
*       | *  *  .  .  .  |   | *              |   | *              |
*       |    *  .  .  ^  |   |    *           |   | *  *           |
*       |       *  .  |  | * |       *        | * | .  .  *        |
*       |          *  |  |   |          *     |   | .  .  .  *     |
*       |             *  |   |             *  |   | .  ---->    *  |
*
* a^{-1} =     z^t        *         d^{-1}     *         z
*
* history:
*   original version for sparslab - tu - 07/03/1996.
*
* references:
*   benzi/meyer/tuma, 1996
*
c
c make_explicit_interface
c
      subroutine ainvsr2(msglvl,msgunit,n,ia,ja,a,ip,jp,ap,
     *  size_p,size_c,size_r,diagtol,
     *  drfl,mi,diag_one,droptyp,imodif,fill,fillmax,
     *  ifillmax,garrow,garcol,info)
c
c parameters
c
      integer msglvl,msgunit,n,size_r,size_c,size_p
      integer ia(*),ja(*)
      double precision a(*)
      integer :: ip(*),jp(*)
      double precision ap(*)
      integer garcol,garrow,droptyp
      integer imodif,diag_one,fill,fillmax,ifillmax,info
      double precision mi,drfl,diagtol
c
c includes
c
      include 'ainvini1.h'
      include 'rank1n.h'
      include 'cd2cs_f90.h'
c
c internals
c
      integer, pointer :: ptr(:),cnr(:),lenr(:),ptc(:),cnc(:),lenc(:)
      double precision, pointer :: h(:),diag(:)
      integer, allocatable :: wn01(:),wn02(:),wn03(:),wn04(:)
      double precision, allocatable :: wr01(:),wr02(:),wr03(:)
c
      integer ierr,gillmur,ind2,ind3,indexr,indexc
      integer ifr,iendru,iendlu,eps,idistr,idistc,lsize
      integer nit,savefill,nreallocr,nreallocc,nja
      double precision size_mult,dv1,one,zero,delta
      parameter(size_mult=1.5d0,one=1.0d0,zero=0.0d0)
      parameter(eps=2.2d-16,delta=1.0d-60)
c
c start of ainvsr2
c
c  -- allocation
c
      nja=ia(n+1)-1
      size_r=max(3*nja,size_r)
      size_c=max(3*nja,size_c)
      allocate(wn01(n+1),wn02(n+1),wn03(n+1),wn04(n+1),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','ainvsr2',
     *  ' allocation error',8,2)
      allocate(wr01(n+1),wr02(n+1),wr03(n+1),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','ainvsr2',
     *  ' allocation error',8,2)
      allocate(ptr(n+1),lenr(n+1),cnr(size_r),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','ainvsr2',
     *  ' allocation error',8,2)
      allocate(ptc(n+1),lenc(n+1),cnc(size_c),h(size_c),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','ainvsr1',
     *  ' allocation error',8,2)
      size_r=size_r
c
c  -- allocate reallocatables
c
      allocate(diag(n),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','ainvsr2',
     *  ' allocation error',8,2)
c
c  -- initialization
c
      info=0
      idistr=1
      idistc=1
      lsize=5
      garcol=0
      garrow=0
      imodif=0
      nit=1
      diagtol=eps
      drfl=drfl+delta
      droptyp=droptyp
      call iwset(n,0,wn01,1)
      call iwset(n,0,wn04,1)
      call dwset(n,zero,wr01,1)
c
c  -- initialize the preconditioner using dynamic matrix formats
c  -- dsr and dsc.
c
      call ainvini1(n,n,ia,ja,ptc,cnc,lenc,h,size_c,
     *  ptr,cnr,lenr,size_r,
     *  iendru,iendlu,idistr,idistc,lsize,ierr)
c
c  -- set initial fill-in parameters
c
      fillmax=0
      ifr=n
      fill=n
c
c  -- echo
c
      if(msglvl.ge.2)
     *    write(*,*)
     *      ' compute symmetric right-looking ainv preconditioner'
c
c  -- the loop nit=1 ...
c
  333 continue
c
c    -- find a pivot
c
        indexr=nit
        indexc=nit
        call ainvrow2(n,ia,ja,a,ptc(1),lenc(1),cnc(1),h(1),
     *    ptr(1),lenr(1),cnr(1),
     *    wn01,wr01,wn02,wr02,wn03,wn04,ind2,indexr,dv1,drfl)
c
c    -- pivot modification and computation of updating vectors
c
        call aicol(ptc(1),lenc(1),cnc(1),h(1),wn03,wr03,ind3,
     *    wr02,ind2,indexc,dv1,diagtol,imodif,gillmur,mi,diag_one)
c
c    -- save diagonal factor
c    -- this block can be considered as an anachronism because in case
c    -- that diag_one.eq.1 this array is not used and in the opposite
c    -- case use of diag can be avoided by setting diagonal elements
c    -- in factor. but we keep that as an monument of the old style
c    -- routines until new classification of preconditioner formats
c    -- is in all the routines.
c
        if(diag_one.eq.1) then
          diag(nit)=one
        else
          diag(nit)=dv1
        end if
c
c    -- test numerical and structural singularity
c
        if(indexr.le.0.or.indexc.le.0) then
          call serr(0,2,'sparslab','ainvsr2_s',
     *      ' matrix is structurally singular',1041,2)
        end if
c
c    -- test whether to update
c
        if(nit.lt.n) then
c
c      -- compute size of fill-in
c
          ifr=ifr-lenc(indexc)
          savefill=fill
c
c      -- echo
c
          if(msglvl.ne.0) then
            if(msglvl.eq.2) then
              if(nit/100.eq.dble(nit)/100.0d0) then
                if(msgunit.eq.0) then
                  write(*,*)' nit =',nit, ' fill = ',fill
                elseif(msgunit.gt.0) then
                  write(msgunit,*)' nit =',nit, ' fill = ',fill
                end if
              end if
            elseif(msglvl.ge.3) then
              if(msgunit.eq.0) then
                write(*,*)' nit =',nit, ' fill = ',fill
              elseif(msgunit.gt.0) then
                write(msgunit,*)' nit =',nit, ' fill = ',fill
              end if
            end if
          end if
c
c    -- sparse rank-one update
c
          call rank1n(msglvl,n,n,ptr,cnr,lenr,size_r,
     *      ptc,cnc,lenc,h,size_c,
     *      indexc,wr01,wn03,wr03,ind3,wn02,wr02,ind2,
     *      iendru,iendlu,drfl,nit,fill,idistc,
     *      garcol,garrow,droptyp,nreallocr,nreallocc)
c
          ifr=ifr+fill-savefill
c
c      -- update fillmax - maximum size of reduced fill-in
c      -- which is obtained in the iteration ifillmax during the
c      -- preconditioner computation
c
          if(fillmax.lt.ifr) then
            ifillmax=nit
            fillmax=ifr
          endif
          nit=nit+1
          go to 333
        end if
c
c  -- format transformation
c
      size_p=fill
c      allocate(ip(n+1),stat=ierr)
c allocate(jp(size_p),stat=ierr)
c      allocate(ap(size_p),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','ainvsr2',
     *  ' allocation error',8,2)
      call cd2cs_f90(n,ptc,diag,cnc,lenc,h,ip,jp,ap,diag_one)
c
c  -- deallocation
c
      deallocate(wn01,wn02,wn03,wn04)
      deallocate(wr01,wr02,wr03,diag)
      deallocate(ptr,cnr,lenr,ptc,cnc,lenc,h)
c
c  -- return
c
      return
c
c end of ainvsr2
c
      end

* (c) sparslab module name=ainvini1
*
* purpose:
*   initialization of dynamic data structures for right-looking
*   ainv preconditioner.
*   input matrix is in csr format.
*
* history:
*   original version for sparslab - tu - 13/04/1996.
*
* parameters:
*  ii  n  number of input matrix columns.
*  ii  m  number of input matrix rows.
*  gi  ia(n+1)/ja(ia(nn+1)-1) input matrix graph.
*  mo  ic(n)/jc(ic(n)+kc(n)-1)/kc(n)/ac(ic(n)+kc(n)-1) output unit
*        matrix in dsc format.
*  ii  max_c  size of the vectors jc/ac.
*  mo  ir(n)/jr(ir(n)+kr(n)-1)/kr(n) output graph of unit matrix
*        in dsr format.
*  ii  max_r  size of the vector jr.
*  io  iendru  position of last entry in ir/jr/kr.
*  io  iendlu  position of last entry in ic/jc/kc/ac.
*  ii  idistr  type of the distribution pattern of the dsc format.
*  ii  idistc  type of the distribution pattern of the dsr format.
*  ii  lsize  size of the distance between entries in ds* formats.
*
c
c make_explicit_interface
c
      subroutine ainvini1(n,m,ia,ja,ic,jc,kc,ac,max_c,ir,jr,kr,max_r,
     *  iendru,iendlu,idistr,idistc,lsize,ierr)
c
c parameters
c
      integer n,m,max_r,max_c,lsize,ierr
      integer idistr,idistc,len
      integer ia(*),ja(*)
      integer, pointer :: ir(:),jr(:),kr(:),ic(:),jc(:),kc(:)
      double precision, pointer :: ac(:)
      integer iendru,iendlu
c
c include
c
c      include 'sparslab.i'
c
c internals
c
      integer i,j,k,l,jstrt,jstop,nja
      integer izero
      double precision one
      parameter(izero=0,one=1.0d0)
c
c start of ainvini1
c
c  -- test/set max_r and max_c
c
      nja=ia(n+1)-1
      if(idistc.eq.1) then
c
c  -- search the matrix rows
c  -- not necessary if appropriate pointers are passed from
c  -- the matrix stored in a triangle
c  -- the input matrix in ia/ja/a should be partially ordered
c
        if(ia(n+1)-1+n*lsize.ge.max_c-1) then
          max_c=2*nja+n*lsize
          deallocate(ac)
          allocate(ac(max_c),stat=ierr)
          if(ierr.ne.0) call serr(0,2,'sparslab','ainvini1',
     *      ' missing space to store the preconditioner',1603,2)
          deallocate(jc)
          allocate(jc(max_c),stat=ierr)
          if(ierr.ne.0) call serr(0,2,'sparslab','ainvini1',
     *      ' missing space to store the preconditioner',1603,2)
        end if
      end if
      if(idistr.eq.1) then
c
c  -- search the matrix rows
c  -- not necessary if appropriate pointers are passed from
c  -- the matrix stored in a triangle
c  -- the input matrix in ia/ja/a should be partially ordered
c
        if(ia(m+1)-1+n*lsize.ge.max_r-1) then
          max_r=2*nja+n*lsize
          deallocate(jr)
          allocate(jr(max_r),stat=ierr)
          if(ierr.ne.0) call serr(0,2,'sparslab','ainvini1',
     *      ' missing space to store the preconditioner',1603,2)
        end if
      end if
c
c  -- clear vectors for row and column indices
c
      call iwset(max_r-1,izero,jr(1),1)
      jr(max_r)=1
      call iwset(max_c-1,izero,jc(1),1)
      jc(max_c)=1
c
c  -- set column pointers
c  -- and insert an unit matrix stored by columns
c
      if(idistc.eq.0.or.idistc.eq.2) then
        if(idistc.eq.0) then
          len=1
          iendlu=len*(m-1)+1
        elseif(idistc.eq.2) then
          i=max_c/3
          len=max(1,i/m)
          iendlu=len*(m-1)+1
        end if
c
c    -- set initial values in the distance len
c
        k=1
        do i=1,n
          ic(i)=k
          kc(i)=1
          jc(k)=i
          ac(k)=one
          k=k+len
        end do
      elseif(idistc.eq.1) then
        k=1
        do i=1,n
          jstrt=ia(i)
          jstop=ia(i+1)-1
          len=jstop-jstrt+1
          do j=jstrt,jstop
            l=ja(j)
            if(l.eq.i) then
              len=j-jstrt+1
              go to 130
            end if
          end do
 130      continue
          ic(i)=k
          kc(i)=1
          jc(k)=i
          ac(k)=one
          k=k+len+lsize
        end do
        iendlu=k-len-lsize
      end if
      if(iendlu.ge.max_c-1) then
        stop 'internal error'
      end if
c
c  -- set row pointers
c  -- and insert row structure of the unit matrix
c
      if(idistr.eq.0.or.idistr.eq.2) then
        if(idistr.eq.0) then
          len=1
          iendru=len*(m-1)+1
        elseif(idistr.eq.2) then
          i=max_r/3
          len=max(1,i/m)
          iendru=len*(m-1)+1
        end if
c
c    -- set initial values
c
        k=1
        do i=1,m
          ir(i)=k
          kr(i)=1
          jr(k)=i
          k=k+len
        end do
      elseif(idistr.eq.1) then
        k=1
        do i=1,m
          jstrt=ia(i)
          jstop=ia(i+1)-1
          len=jstop-jstrt+1
          do j=jstrt,jstop
            l=ja(j)
            if(l.eq.i) then
              len=jstop-j+1
              go to 110
            end if
          end do
 110      continue
          ir(i)=k
          kr(i)=1
          jr(k)=i
          k=k+len+lsize
        end do
        iendru=k-len-lsize
      end if
      if(iendru.ge.max_r-1) then
        stop 'internal error'
      end if
c
c  -- return
c
      return
c
c end of ainvini1
c
      end

* (c) sparslab module name=ainvrow2
*
* purpose:
*   find a_{indexrc}^t*z^{(indexrc-1)}[i..n]
*   for right-looking stabilized ainv preconditioner.
*
* history:
*   original version for sparselab - tu - 29-03-96.
*   updated - tu - 18/5/1999.
*
* parameters:
*
      subroutine ainvrow2(n,ia,ja,a,ptc,lenc,cnc,h,ptr,lenr,cnr,
     *  wn01,wr01,wn02,wr02,wn03,wn04,ind2,indexrc,dv1,drfl)
c
      integer n
      integer ia(*),ja(*)
      double precision drfl
      double precision a(*)
      integer lenc(*),ptr(*),lenr(*),wn01(*),cnr(*)
      double precision h(*)
      double precision wr01(*),wr02(*)
      integer wn02(*),wn03(*),wn04(*),ptc(*),cnc(*)
c
c internals
c
      integer i,j,k,i1,ii,l
      integer ind3,ind
      integer ind2,indexrc
      integer indk,jstrt,jstop,istrt,istop
      double precision zero,dv1,sum,mi,mifl
      parameter(zero=0.0d0)
      parameter(mi=0.1d0)
c
c start of ainvrow2
c
c  -- find columns intersecting row a_{indexrc}.
c  -- store the indices in wn02(1..ind2).
c  -- indk denotes position of the pivot column indexrc in wn02.
c
c             p_{nit}^{(i-1)} = z_i^{(i-1)}^t*a*z_{nit}^{(i-1)}
c
      n=n
      ind3=0
      mifl=mi*drfl
      istrt=ptc(indexrc)
      istop=ptc(indexrc)+lenc(indexrc)-1
      do i=istrt,istop
        ii=cnc(i)
        jstrt=ia(ii)
        jstop=ia(ii+1)-1
        do j=jstrt,jstop
          k=ja(j)
          l=wn04(k)
          if(l.eq.0) then
            ind3=ind3+1
            wr01(k)=h(i)*a(j)
            wn04(k)=ind3
            wn03(ind3)=k
          else
            wr01(k)=wr01(k)+h(i)*a(j)
          end if
        end do
      end do
c
      ind2=0
      do i=1,ind3
        j=wn03(i)
        jstrt=ptr(j)
        jstop=jstrt+lenr(j)-1
        do j=jstrt,jstop
          k=cnr(j)
          if(wn01(k).eq.0) then
            wn01(k)=j
            ind2=ind2+1
            wn02(ind2)=k
            if(k.eq.indexrc) indk=ind2
          end if
        end do
      end do
c
c  -- loop of dot products between
c  -- columns intersecting row indexrc and
c  -- column z_{indexrc}
c
      do i=1,ind2
        sum=zero
        i1=wn02(i)
        wn01(i1)=0
        jstrt=ptc(i1)
        jstop=jstrt+lenc(i1)-1
        do j=jstrt,jstop
          k=cnc(j)
          if(wr01(k).ne.zero) then
            sum=sum+h(j)*wr01(k)
          end if
        end do
        wr02(i)=sum
      end do
c
      do i=1,ind3
        j=wn03(i)
        wn04(j)=0
        wr01(j)=zero
      end do
c
c  -- get the diagonal pivot dv1
c
      dv1=wr02(indk)
c
c  -- remove the pivot index/entry from the list
c  -- wr02(1..ind2)
c
      wn02(indk)=wn02(ind2)
      wr02(indk)=wr02(ind2)
      ind2=ind2-1
c
c  -- remove small elements
c
      ind=0
      do i=1,ind2
        if(abs(wr02(i)).ge.mifl) then
          ind=ind+1
          wr02(ind)=wr02(i)
          wn02(ind)=wn02(i)
        end if
      end do
      ind2=ind
c
c  -- return
c
      return
c
c end of ainvrow2
c
      end


* (c) sparslab module name=aicol
* purpose:
*   modify a pivot for sparse approximate inverse preconditioner ainv,
*   invert dv1
*   get a column vector wn03/wr03(1..ind3) and scale
*   a row vector wn02/wr02(1..ind2) by dv1.
*   these operations correspond to initialization
*   before sparse rank-1 update in each step
*   of the approximate inverse preconditioner procedure
*   modification strategies: plain
*                            gill-murray
*   the choice is based on the gillmur variable from ipar
*   if diag_one=1 then scale column z_i by sqrt(pivot_i)
*
* history:
*   original version for sparslab - tu - 12/04/1996.
*
* parameters:
*
      subroutine aicol(ptc,lenc,cnc,h,wn03,wr03,ind3,
     *  wr02,ind2,indexc,dv1,ptol,imodif,gillmur,mi,diag_one)
c
c parameters
c
c      integer, pointer :: ptc(:),lenc(:),cnc(:)
c      double precision, pointer :: h(:)
      integer :: ptc(*),lenc(*),cnc(*)
      double precision :: h(*)
      integer wn03(*)
      double precision wr02(*),wr03(*)
      double precision dv1,ptol,mi
      integer ind2,ind3,indexc,imodif,gillmur
      integer diag_one
c
c internals
c
      integer istrt,istop,i
      double precision zero,temp,temp1,one
      logical never
      parameter(zero=0.0d0,one=1.0d0,never=.false.)
c
c start of aicol
c
      istrt=ptc(indexc)
      istop=istrt+lenc(indexc)-1
c
c  -- test dv1 - if dv1 too small then modify
c  -- get vector wn03/wr03(1..ind3)
c
      ind3=0
      if(gillmur.eq.1.and.never) then
        temp=zero
        do i=istrt,istop
          ind3=ind3+1
          wn03(ind3)=cnc(i)
          wr03(ind3)=h(i)
          temp=max(temp,abs(h(i)))
        end do
        temp1=zero
        do i=1,ind2
          temp1=max(temp1,wr02(i))
        end do
        dv1=max(ptol,mi*temp*temp1)
      else
        if(dv1.lt.ptol) then
          imodif=imodif+1
          dv1=max(abs(dv1),ptol)
        end if
        do i=istrt,istop
          ind3=ind3+1
          wn03(ind3)=cnc(i)
          wr03(ind3)=h(i)
        end do
      end if
c
c  -- invert dv1
c
      dv1=one/dv1
c
c  -- if diag_one is set then scale pivot column
c
      if(diag_one.eq.1) then
        ind3=0
        temp=sqrt(dv1)
        dv1=temp
        do i=istrt,istop
          ind3=ind3+1
          h(i)=h(i)*temp
          wr03(ind3)=h(i)
        end do
      end if
c
c  -- scale the vector wn02/wr02(1..ind2)
c
      do i=1,ind2
        wr02(i)=wr02(i)*dv1
      end do
c
c  -- return
c
      return
c
c end of aicol
c
      end


* (c) sparslab module name=rank1n
*
* purpose:
*   sparse rank-one update in dynamic data structures.
*   it updates a submatrix given by row indices
*   wn03(1..ind3) and column indices wn02(1..ind2)
*   and removes pivot column indexk and row indexj
*   from the data structures.
*   nonsymmetric version.
*
*   vector wn03/wr03(1..ind3) corresponds to a pivot column.
*   input matrix is stored in the dynamic data structures
*   cnr,ptr,lenr/cnc,ptc,lenc
*
* history:
*   original version for sparslab - tu - 01/05/2001.
*   based on the codes sig.f, nig.f, aig.f, ainvupd1.f
*
c
c make_explicit_interface
c
      subroutine rank1n(scrlvl,n,m,ptr,cnr,lenr,
     *  max_r,ptc,cnc,lenc,h,max_c,
     *  indexk,wr01,wn03,wr03,ind3,wn02,wr02,ind2,
     *  iendru,iendlu,drfl,nit,if2,idist,garcol,garrow,droptyp,
     *  nreallocr,nreallocc)
c
c parameters
c
      integer n,m,droptyp,nreallocr,nreallocc,scrlvl
      integer nit,if2,max_r,max_c
      integer, pointer :: ptr(:),cnr(:),lenr(:),ptc(:),cnc(:),lenc(:)
      double precision, pointer :: h(:)
      integer iendru,iendlu,nrow,garrow,garcol,indexk,ind2,ind3
      integer idist
      integer wn02(*),wn03(*)
      double precision wr01(*),wr02(*),wr03(*),temp
      double precision drfl
c
c internals
c
      integer istrt,istop,ibeg,ncol,i1,ierr
      integer i,j,j1,j2,j3,iqcmpr,jstrt,jstop,incr
      integer ll1,k,len,newmax_r,newmax_c,oldmax_r,oldmax_c
      double precision zero,one,drfl_rel,avg_size
      character dropping*8
      parameter (zero=0.0d0,one=1.0d0)
c
c interfaces
c
      interface
      subroutine realloci(ja,arrsize,arrincr,new_arrsize,ierr)
      integer arrsize,arrincr,new_arrsize,ierr
      integer, pointer, dimension(:) :: ja
      end
      end interface

      interface
      subroutine reallocr(aa,arrsize,arrincr,new_arrsize,ierr)
      integer arrsize,arrincr,new_arrsize,ierr
      double precision, pointer, dimension(:) :: aa
      end
      end interface
c
c  -- start of rank1n
c
c  -- set dropping
c
      if(droptyp.eq.0) then
        dropping='absolute'
      elseif(droptyp.eq.1) then
        dropping='relative'
      end if
c
c  -- define new parameters
c
      ibeg=1
c
c  -- remove indexk from the rows
c  -- which intersect column indexk.
c
      istrt=ptc(indexk)
      istop=istrt+lenc(indexk)-1
      do i=istrt,istop
        k=cnc(i)
        jstrt=ptr(k)
        jstop=jstrt+lenr(k)-1
        do j=jstrt,jstop
          if(cnr(j).eq.indexk) then
            cnr(j)=cnr(jstop)
            cnr(jstop)=0
            lenr(k)=lenr(k)-1
            go to 25
          end if
        end do
 25     continue
      end do
c
c  -- main loop for columns with indices from wn02(1..ind2)
c
      do ll1=1,ind2
c
c    -- get index of the updated column and
c    -- corresponding column value of the rank-1 matrix
c
        ncol=wn02(ll1)
        temp=wr02(ll1)
c
c    -- scatter and scale the pivot column from wr03 into wr01
c
        do i=1,ind3
          j=wn03(i)
          wr01(j)=-wr03(i)*temp
        end do
c
c    -- define start and end of the column ncol
c
        istrt=ptc(ncol)
        istop=istrt+lenc(ncol)-1
c
c    -- compute average size of an original nonzero
c    -- for relative dropping
c
        if(dropping.eq.'relative') then
          avg_size=zero
          do j=istrt,istop
            temp=h(j)
            avg_size=avg_size+abs(temp)
          end do
          len=istop-istrt+1
          if(len.ne.0) then
            avg_size=avg_size/dble(len)
          end if
        end if
c
c    -- get drop tolerance
c
        if(dropping.eq.'relative') then
          drfl_rel=drfl*avg_size
        else
          drfl_rel=drfl
        end if
c
c    -- update existing entries of the column ncol
c
        do j=istrt,istop
          k=cnc(j)
          if(wr01(k).ne.zero) then
c
c        -- update
c
            h(j)=h(j)+wr01(k)
            wr01(k)=zero
            if(abs(h(j)).lt.drfl_rel.and.k.ne.ncol) then
c
c          -- mark the value which is to be dropped
c
              cnc(j)=0
c
c          -- define start and end of the corresponding row list
c
              jstrt=ptr(k)
              jstop=jstrt+lenr(k)-1
c
c          -- find index ncol in the row list corresponding
c          -- to an entry to be dropped
c
              do i1=jstrt,jstop
                if(cnr(i1).eq.ncol) then
c
c              -- index found
c
                  cnr(i1)=cnr(jstop)
                  cnr(jstop)=0
                  lenr(k)=lenr(k)-1
c
c              -- jump outside
c
                  go to 204
                end if
              end do
 204          continue
            end if
          end if
        end do
c
c -- remove marked candidates to be dropped
c -- shrink remaining entries in the column ncol between
c -- istrt and istop
c
        i=istrt
        do j=istrt,istop
          if(cnc(j).eq.0) then
            if2=if2-1
            go to 220
          else
            cnc(i)=cnc(j)
            h(i)=h(j)
            i=i+1
          end if
 220      continue
        end do
        lenc(ncol)=i-istrt
        do j=i,istop
          cnc(j)=0
        end do
        istop=i-1
c
c    -- insert fill-in into column ncol
c
        do 80 i=1,ind3
          nrow=wn03(i)
c
c      -- test fill-in size
c
          if(abs(wr01(nrow)).gt.drfl_rel) then
c
c        -- look for a free space at the beginning of a column list
c
            if(istrt.gt.1) then
              if(cnc(istrt-1).eq.0) then
c
c            -- add a the new entry at the beginning of the column list
c
                istrt=istrt-1
                cnc(istrt)=nrow
                h(istrt)=wr01(nrow)
                lenc(ncol)=lenc(ncol)+1
                ptc(ncol)=ptc(ncol)-1
                if2=if2+1
                go to 550
              end if
            endif
c
c        -- look for a free space at the end of a column list
c
 400        if(cnc(istop+1).eq.0) then
c
c          -- add a new element at the end of a column list
c
              iqcmpr=0
              istop=istop+1
              cnc(istop)=nrow
              h(istop)=wr01(nrow)
              if2=if2+1
              iendlu=max(iendlu,istop)
              lenc(ncol)=lenc(ncol)+1
            else
              j1=istop-istrt+1
              if((iendlu+j1).lt.max_c-1) then
c
c            -- copy a column at the current end
c            -- of the column working space
c
                j2=istrt-1-iendlu
                do j=iendlu+1,iendlu+j1
                  j3=j2+j
                  cnc(j)=cnc(j3)
                  h(j)=h(j3)
                  cnc(j3)=0
                end do
                iqcmpr=0
                istrt=iendlu+1
                ptc(ncol)=istrt
                istop=iendlu+j1
                go to 400
              else
c
c            -- garbage collection of the column working space
c
                if(iqcmpr.eq.1.or.
     *             (dble(iendlu)/dble(nit))*dble(n).gt.max_c) then
c
c              -- reallocation of cnc,h
c
                  oldmax_c=max_c
                  incr=min(max_c,int(0.5d0*(n*n-nit*(nit-1)/2)))
                  call realloci(cnc,max_c,incr,newmax_c,ierr)
                  call reallocr(h,max_c,incr,newmax_c,ierr)
                  max_c=newmax_c
                  nreallocc=nreallocc+1
                  call iwset(max_c-oldmax_c+1,0,cnc(oldmax_c),1)
                  cnc(max_c)=1
                end if
                iqcmpr=1
                garcol=garcol+1
                if(scrlvl.gt.0) then
                  write(*,*) ' column garbage collection no.',garcol
                end if
                call garclds(n,ptc(1),cnc(1),lenc(1),h(1),max_c,
     *            ibeg,iendlu,nit,idist)
                istrt=ptc(ncol)
                istop=istrt+lenc(ncol)-1
                go to 400
              end if
            end if
 550        continue
c
c        -- add a new index ncol into a row list corresponding to
c        -- the new fill-in position nrow
c
            jstrt=ptr(nrow)
            jstop=jstrt+lenr(nrow)-1
c
c        -- look for a free space at the beginning of a row list
c
            if(jstrt.gt.1) then
              if(cnr(jstrt-1).eq.0) then
c
c            -- add a new element at the beginning of a row list
c
                jstrt=jstrt-1
                cnr(jstrt)=ncol
                lenr(nrow)=lenr(nrow)+1
                ptr(nrow)=ptr(nrow)-1
                go to 450
              end if
            endif
c
c          -- look for a free space at the end of a row list
c
 440        if(cnr(jstop+1).eq.0) then
c
c          -- add a new element at the end of a row list
c
              iqcmpr=0
              jstop=jstop+1
              cnr(jstop)=ncol
              iendru=max(iendru,jstop)
              lenr(nrow)=lenr(nrow)+1
            else
c
c          -- copy a row at the current end of the row working space
c
              j1=jstop-jstrt+1
              i1=jstrt-1-iendru
              if((iendru+j1).lt.max_r-1) then
                do 540 j2=1+iendru,j1+iendru
                  j3=i1+j2
                  cnr(j2)=cnr(j3)
                  cnr(j3)=0
 540            continue
                iqcmpr=0
                jstrt=iendru+1
                ptr(nrow)=jstrt
                jstop=iendru+j1
                go to 440
              else
c
c            -- garbage collection of the row working space
c
                if(iqcmpr.eq.1.or.
     *             (dble(iendru)/dble(nit))*dble(n).gt.max_r) then
c
c              -- reallocation of cnr
c
                  oldmax_r=max_r
                  incr=min(max_r,int(0.5d0*(n*n-nit*(nit-1)/2)))
                  call realloci(cnr,max_r,incr,newmax_r,ierr)
                  max_r=newmax_r
                  nreallocr=nreallocr+1
                  call iwset(max_r-oldmax_r+1,0,cnr(oldmax_r),1)
                  cnr(max_r)=1
                end if
                iqcmpr=1
                garrow=garrow+1
                if(scrlvl.gt.0) then
                  write(*,*) ' row garbage collection no.',garrow
                end if
                call garclds2(m,ptr(1),cnr(1),lenr(1),max_r,ibeg,
     *            iendru,idist)
                jstrt=ptr(nrow)
                jstop=jstrt+lenr(nrow)-1
                go to 440
              end if
            end if
 450        continue
          end if
c
c      -- clear nrow position in wr01
c
          wr01(nrow)=zero
c
c      -- end of fill-in cycle
c
 80     continue
c
c   -- end of loop for columns
c
      end do
c
c  -- return
c
      return
c
c end of rank1n
c
      end


* (c) sparslab module name=rank1n2
*
* purpose:
*   sparse rank-one update in dynamic data structures.
*   it updates a submatrix given by row indices
*   wn03(1..ind3) and column indices wn02(1..ind2)
*   and removes pivot column indexk and row indexj
*   from the data structures.
*   nonsymmetric version.
*   removes processed columns from dynamic data structures.
*   advantageous for such algorithms which will not use these
*   columns later (like rif et al.)
*
*   vector wn03/wr03(1..ind3) corresponds to a pivot column.
*   input matrix is stored in the dynamic data structures
*   cnr,ptr,lenr/cnc,ptc,lenc
*
* history:
*   original version for sparslab - tu - 01/05/2001.
*   based on the codes sig.f, nig.f, aig.f, ainvupd1.f
*
c
c make_explicit_interface
c
      subroutine rank1n2(scrlvl,n,m,ptr,cnr,lenr,
     *  max_r,ptc,cnc,lenc,h,max_c,
     *  indexk,wr01,wn03,wr03,ind3,wn02,wr02,ind2,
     *  iendru,iendlu,drfl,nit,if2,idist,garcol,garrow,droptyp,
     *  nreallocr,nreallocc)
c
c parameters
c
      integer scrlvl,n,m,droptyp,nreallocr,nreallocc
      integer nit,if2,max_r,max_c
      integer, pointer :: ptr(:),cnr(:),lenr(:),ptc(:),cnc(:),lenc(:)
      double precision, pointer :: h(:)
      integer iendru,iendlu,nrow,garrow,garcol,indexk,ind2,ind3
      integer idist
      integer wn02(*),wn03(*)
      double precision wr01(*),wr02(*),wr03(*),temp
      double precision drfl
c
c interfaces
c
      interface
      subroutine realloci(ja,arrsize,arrincr,new_arrsize,ierr)
      integer arrsize,arrincr,new_arrsize,ierr
      integer, pointer, dimension(:) :: ja
      end
      end interface

      interface
      subroutine reallocr(aa,arrsize,arrincr,new_arrsize,ierr)
      integer arrsize,arrincr,new_arrsize,ierr
      double precision, pointer, dimension(:) :: aa
      end
      end interface
c
c internals
c
      integer istrt,istop,ibeg,ncol,i1,ierr
      integer i,j,j1,j2,j3,iqcmpr,jstrt,jstop,incr
      integer ll1,k,len,newmax_r,newmax_c,oldmax_r,oldmax_c
      double precision zero,one,drfl_rel,avg_size
      character dropping*8
      parameter (zero=0.0d0,one=1.0d0)
c
c  -- start of rank1n2
c
c  -- set dropping
c
      if(droptyp.eq.0) then
        dropping='absolute'
      elseif(droptyp.eq.1) then
        dropping='relative'
      end if
c
c  -- define new parameters
c
      ibeg=1
c
c  -- remove indexk from the rows
c  -- which intersect column indexk.
c
      istrt=ptc(indexk)
      istop=istrt+lenc(indexk)-1
      do i=istrt,istop
        k=cnc(i)
        jstrt=ptr(k)
        jstop=jstrt+lenr(k)-1
        do j=jstrt,jstop
          if(cnr(j).eq.indexk) then
            cnr(j)=cnr(jstop)
            cnr(jstop)=0
            lenr(k)=lenr(k)-1
            go to 25
          end if
        end do
 25     continue
      end do
c
c  -- main loop for columns with indices from wn02(1..ind2)
c
      do ll1=1,ind2
c
c    -- get index of the updated column and
c    -- corresponding column value of the rank-1 matrix
c
        ncol=wn02(ll1)
        temp=wr02(ll1)
c
c    -- scatter and scale the pivot column from wr03 into wr01
c
        do i=1,ind3
          j=wn03(i)
          wr01(j)=-wr03(i)*temp
        end do
c
c    -- define start and end of the column ncol
c
        istrt=ptc(ncol)
        istop=istrt+lenc(ncol)-1
c
c    -- compute average size of an original nonzero
c    -- for relative dropping
c
        if(dropping.eq.'relative') then
          avg_size=zero
          do j=istrt,istop
            temp=h(j)
            avg_size=avg_size+abs(temp)
          end do
          len=istop-istrt+1
          if(len.ne.0) then
            avg_size=avg_size/dble(len)
          end if
        end if
c
c    -- get drop tolerance
c
        if(dropping.eq.'relative') then
          drfl_rel=drfl*avg_size
        else
          drfl_rel=drfl
        end if
c
c    -- update existing entries of the column ncol
c
        do j=istrt,istop
          k=cnc(j)
          if(wr01(k).ne.zero) then
c
c        -- update
c
            h(j)=h(j)+wr01(k)
            wr01(k)=zero
            if(abs(h(j)).lt.drfl_rel) then
c
c          -- mark the value which is to be dropped
c
              cnc(j)=0
c
c          -- define start and end of the corresponding row list
c
              jstrt=ptr(k)
              jstop=jstrt+lenr(k)-1
c
c          -- find index ncol in the row list corresponding
c          -- to an entry to be dropped
c
              do i1=jstrt,jstop
                if(cnr(i1).eq.ncol) then
c
c              -- index found
c
                  cnr(i1)=cnr(jstop)
                  cnr(jstop)=0
                  lenr(k)=lenr(k)-1
c
c              -- jump outside
c
                  go to 204
                end if
              end do
 204          continue
            end if
          end if
        end do
c
c -- remove marked candidates to be dropped
c -- shrink remaining entries in the column ncol between
c -- istrt and istop
c
        i=istrt
        do j=istrt,istop
          if(cnc(j).eq.0) then
            if2=if2-1
            go to 220
          else
            cnc(i)=cnc(j)
            h(i)=h(j)
            i=i+1
          end if
 220      continue
        end do
        lenc(ncol)=i-istrt
        do j=i,istop
          cnc(j)=0
        end do
        istop=i-1
c
c    -- insert fill-in into column ncol
c
        do 80 i=1,ind3
          nrow=wn03(i)
c
c      -- test fill-in size
c
          if(abs(wr01(nrow)).gt.drfl_rel) then
c
c        -- look for a free space at the beginning of a column list
c
            if(istrt.gt.1) then
              if(cnc(istrt-1).eq.0) then
c
c            -- add a the new entry at the beginning of the column list
c
                istrt=istrt-1
                cnc(istrt)=nrow
                h(istrt)=wr01(nrow)
                lenc(ncol)=lenc(ncol)+1
                ptc(ncol)=ptc(ncol)-1
                if2=if2+1
                go to 550
              end if
            endif
c
c        -- look for a free space at the end of a column list
c
 400        if(cnc(istop+1).eq.0) then
c
c          -- add a new element at the end of a column list
c
              iqcmpr=0
              istop=istop+1
              cnc(istop)=nrow
              h(istop)=wr01(nrow)
              if2=if2+1
              iendlu=max(iendlu,istop)
              lenc(ncol)=lenc(ncol)+1
            else
              j1=istop-istrt+1
              if((iendlu+j1).lt.max_c-1) then
c
c            -- copy a column at the current end
c            -- of the column working space
c
                j2=istrt-1-iendlu
                do j=iendlu+1,iendlu+j1
                  j3=j2+j
                  cnc(j)=cnc(j3)
                  h(j)=h(j3)
                  cnc(j3)=0
                end do
                iqcmpr=0
                istrt=iendlu+1
                ptc(ncol)=istrt
                istop=iendlu+j1
                go to 400
              else
c
c            -- garbage collection of the column working space
c
                if(iqcmpr.eq.1.or.
     *             (dble(iendlu)/dble(nit))*dble(n).gt.max_c) then
c
c              -- reallocation of cnc,h
c
                  oldmax_c=max_c
                  incr=min(max_c,int(0.5d0*(n*n-nit*(nit-1)/2)))
                  call realloci(cnc,max_c,incr,newmax_c,ierr)
                  call reallocr(h,max_c,incr,newmax_c,ierr)
                  max_c=newmax_c
                  nreallocc=nreallocc+1
                  call iwset(max_c-oldmax_c+1,0,cnc(oldmax_c),1)
                  cnc(max_c)=1
                end if
                iqcmpr=1
                garcol=garcol+1
                if(scrlvl.gt.0) then
                  write(*,*) ' column garbage collection no.',garcol
                end if
                call garclds(n,ptc(1),cnc(1),lenc(1),h(1),max_c,
     *            ibeg,iendlu,nit,idist)
                istrt=ptc(ncol)
                istop=istrt+lenc(ncol)-1
                go to 400
              end if
            end if
 550        continue
c
c        -- add a new index ncol into a row list corresponding to
c        -- the new fill-in position nrow
c
            jstrt=ptr(nrow)
            jstop=jstrt+lenr(nrow)-1
c
c        -- look for a free space at the beginning of a row list
c
            if(jstrt.gt.1) then
              if(cnr(jstrt-1).eq.0) then
c
c            -- add a new element at the beginning of a row list
c
                jstrt=jstrt-1
                cnr(jstrt)=ncol
                lenr(nrow)=lenr(nrow)+1
                ptr(nrow)=ptr(nrow)-1
                go to 450
              end if
            endif
c
c          -- look for a free space at the end of a row list
c
 440        if(cnr(jstop+1).eq.0) then
c
c          -- add a new element at the end of a row list
c
              iqcmpr=0
              jstop=jstop+1
              cnr(jstop)=ncol
              iendru=max(iendru,jstop)
              lenr(nrow)=lenr(nrow)+1
            else
c
c          -- copy a row at the current end of the row working space
c
              j1=jstop-jstrt+1
              i1=jstrt-1-iendru
              if((iendru+j1).lt.max_r-1) then
                do 540 j2=1+iendru,j1+iendru
                  j3=i1+j2
                  cnr(j2)=cnr(j3)
                  cnr(j3)=0
 540            continue
                iqcmpr=0
                jstrt=iendru+1
                ptr(nrow)=jstrt
                jstop=iendru+j1
                go to 440
              else
c
c            -- garbage collection of the row working space
c
                if(iqcmpr.eq.1.or.
     *             (dble(iendru)/dble(nit))*dble(n).gt.max_r) then
c
c              -- reallocation of cnr
c
                  oldmax_r=max_r
                  incr=min(max_r,int(0.5d0*(n*n-nit*(nit-1)/2)))
                  call realloci(cnr,max_r,incr,newmax_r,ierr)
                  max_r=newmax_r
                  nreallocr=nreallocr+1
                  call iwset(max_r-oldmax_r+1,0,cnr(oldmax_r),1)
                  cnr(max_r)=1
                end if
                iqcmpr=1
                garrow=garrow+1
                if(scrlvl.gt.0) then
                  write(*,*) ' row garbage collection no.',garrow
                end if
                call garclds2(m,ptr(1),cnr(1),lenr(1),max_r,ibeg,
     *            iendru,idist)
                jstrt=ptr(nrow)
                jstop=jstrt+lenr(nrow)-1
                go to 440
              end if
            end if
 450        continue
          end if
c
c      -- clear nrow position in wr01
c
          wr01(nrow)=zero
c
c      -- end of fill-in cycle
c
 80     continue
c
c   -- end of loop for columns
c
      end do
c
c  -- remove column indexk
c
      istrt=ptc(indexk)
      istop=istrt+lenc(indexk)-1
      do i=istrt,istop
        cnc(i)=0
      end do
      if2=if2-lenc(indexk)
      lenc(indexk)=0
      ptc(indexk)=1
c
c  -- return
c
      return
c
c end of rank1n2
c
      end


* (c) sparslab module name=cd2cs_f90
*
* purpose:
*    matrix format transformation routine cd --> cs
*    (some variants of compressed dynamic --> compressed sparse)
*
* history:
*   original version for sparslab - tu - 25/07/2000.
*
c
c make_explicit_interface
c
      subroutine cd2cs_f90(n,ia,da,ja,la,aa,ib,jb,ab,diag_one)
c
c parameters
c
      integer n,diag_one
      integer, pointer :: ia(:),ja(:),la(:)
      double precision, pointer :: da(:),aa(:)
      integer ib(*), jb(*)
      double precision ab(*)
c
c internals
c
      include 'cd2cs3_f90.h'
      include 'srtcs1_f90.h'
      integer i
      double precision one
      parameter(one=1.0d0)
c
c start of cd2cs_f90
c
      call cd2cs3_f90(n,ia,la,ja,aa,ib,jb,ab)
      if(diag_one.eq.0) then
        call srtcs1_f90(n,ib,jb,ab)
        do i=1,n
          ab(ib(i+1)-1)=da(i)
        end do
      end if
c
c  -- return
c
      return
c
c end of cd2cs_f90
c
      end

* (c) sparslab module name=cd2cs3_f90
*
* purpose:
*  shrink dynamic data structures in ia/ja/la/aa into
*  ib/jb/ab such that i-th column is moved
*  from the positions ia(i), ..., ia(i)+la(i)-1
*  to the positions ib(i), ..., ib(i+1)-1
*
* history:
*   original version for sparslab - tu - 15/04/1996.
*
c
c make_explicit_interface
c
      subroutine cd2cs3_f90(n,ia,la,ja,aa,ib,jb,ab)
c
c parameters
c
      integer n
      integer, pointer :: ia(:),ja(:),la(:)
      integer ib(*), jb(*)
      double precision, pointer :: aa(:)
      double precision ab(*)
c
c internals
c
      integer ind,i,j,jstrt,jstop
c
c start of cd2cs3_f90
c
c  -- copy the data
c
      ind=1
      ib(1)=ind
      do i=1,n
        jstrt=ia(i)
        jstop=ia(i)+la(i)-1
        do j=jstrt,jstop
          jb(ind)=ja(j)
          ab(ind)=aa(j)
          ind=ind+1
        end do
        ib(i+1)=ind
      end do
c
c  -- return
c
      return
c
c end of cd2cs3_f90
c
      end


* (c) sparslab module name=srtcs1
*
* purpose:
*   sort a sparse matrix stored in a cs storage format.
*
* history:
*   original version for sparslab - tu - 01/06/1996.
*
* parameters:
*   ii  n  matrix dimension.
*   ou  ia(n+1)/ja(ia(n+1)-1/a(ia(n+1)-1  updated matrix cs object.
*
      subroutine srtcs1(n,ia,ja,a)
c
c globals
      integer n,ia(*),ja(*)
      double precision a(*)
c
c locals
      integer i,jlen,jstrt,jstop
c
c start of srtcs1
c
      jstop=ia(1)
      do i=1,n
        jstrt=jstop
        jstop=ia(i+1)
        jlen=jstop-jstrt
        call uxvsr4(jlen,ja(jstrt),a(jstrt))
      end do
c
      return
c
c end of srtcs1
c
      end

* (c) sparslab module name=srtcs1_f90
*
* purpose:
*   sort a sparse matrix stored in a cs storage format.
*
* history:
*   original version for sparslab - tu - 01/06/1996.
*
* parameters:
*   ii  n  matrix dimension.
*   ou  ia(n+1)/ja(ia(n+1)-1/a(ia(n+1)-1  updated matrix cs object.
c
c make_explicit_interface
c
      subroutine srtcs1_f90(n,ia,ja,a)
c
c parameters
c
      integer n
      integer ia(*),ja(*)
      double precision a(*)
c
c internals
c
      integer i,jlen,jstrt,jstop
c
c start of srtcs1_f90
c
      jstop=ia(1)
      do i=1,n
        jstrt=jstop
        jstop=ia(i+1)
        jlen=jstop-jstrt
        call uxvsr4(jlen,ja(jstrt),a(jstrt))
c        call uxvsr4_f90(jlen,ja(jstrt),a(jstrt))
      end do
c
c  -- return
c
      return
c
c end of srtcs1_f90
c
      end

* (c) sparslab module name=garclds
*
* purpose:
*   garbage collection of a matrix in ds* format.
*
* history:
*   original version for sparslab - tu - 13/4/1996.
*
* parameters:
*  ii  n  dimension of the matrix in ds* format.
*  mu  ia(n)/ja(ia(n)+ka(n)-1)/ka(n)/aa(ia(n)+ka(n)-1) updated
*        matrix in ds* format.
*  ii  max_a  size of the vectors ja/aa.
*  ii  ibeg  position of first entry in ia/ja/ka/aa.
*  io  iend  position of last entry in ia/ja/ka/aa.
*  ii  nit  step number (used for modifications of distance parameter).
*  ii  idist  type of the distribution pattern of the ds* format.
*
      subroutine garclds(n,ia,ja,ka,aa,max_a,
     *  ibeg,iend,nit,idist)
c
c parameters
c
      integer n,ibeg,iend,max_a,nit,idist
      integer ia(*),ka(*),ja(*)
      double precision aa(*)
c
c internals
c
      integer i,j,k,l,istrt,isize
      integer izero
      double precision temp
      parameter(izero=0)
c
c start of garclds
c
c  -- compute the minimum necessary size of the compressed output
c
      isize=0
      do i=1,n
        if(ka(i).ne.0) then
          k=ia(i)
          l=ja(k)
          ja(k)=-i
          ia(i)=l
          isize=isize+ka(i)
        else
          ia(i)=ibeg
        end if
      end do
c
c  -- compute distance parameter idist
c
      if(idist.eq.-1) then
c
c  -- -- original case ldist.eq.'a'
c
        idist=max(0,(max_a-1-isize)/(2*(n-nit+1)))
      elseif(idist.gt.0) then
c
c  -- -- original case ldist.eq.'p'
c
        idist=min(max(0,(max_a-1-isize)/(2*(n-nit+1))),idist)
      else
c
c  -- -- original case ldist.eq.'n'
c
        idist=0
      endif
c
c  -- the garbage collection
c
      istrt=ibeg
      i=ibeg
 300  continue
      k=ja(i)
c
c  -- -- find a start of some row/column
c
      if(k.lt.0) then
        temp=aa(i)
c
c  -- -- move the row/column
c
        do j=1,ka(-k)-1
          ja(istrt+j)=ja(i+j)
          aa(istrt+j)=aa(i+j)
        end do
c
c  -- -- move the head and redefine the starting pointer
c
        aa(istrt)=temp
        ja(istrt)=ia(-k)
        ia(-k)=istrt
c
c  -- -- redefine starting pointers of both the source and target
c
        istrt=istrt+ka(-k)
        i=i+ka(-k)
c
c  -- -- clear some space (typically not large)
c
        do j=1,min(i-istrt,idist)
          ja(istrt)=0
          istrt=istrt+1
        end do
      else
c
c  -- -- move further
c
        i=i+1
      end if
c
c  -- -- end of coded loop
c
      if(i.lt.iend) go to 300
c
c  -- clear up to the end
c
      call iwset(iend-istrt+1,izero,ja(istrt),1)
c
c  -- redefine end
c
      iend=istrt-1
c
c  -- return
c
      return
c
c end of garclds
c
      end

* (c) sparslab module name=garclds2
*
* purpose:
*   garbage collection of a graph in ds* format.
*
* history:
*   original version for sparslab - tu - 13/4/1996.
*
* parameters:
*  ii  n  dimension of the graph in ds* format.
*  mu  ia(n)/ja(ia(n)+ka(n)-1)/ka(n) updated graph in ds* format.
*  ii  max_a  size of the vectors ja/aa.
*  ii  ibeg  position of first entry in ia/ja/ka/aa.
*  io  iend  position of last entry in ia/ja/ka/aa.
*  ii  idist  type of the distribution pattern of the ds* format.
*
      subroutine garclds2(n,ia,ja,ka,max_a,ibeg,iend,idist)
c
c parameters
c
      integer n,ibeg,iend,max_a,idist
      integer ia(*),ka(*),ja(*)
c
c internals
c
      integer i,j,k,l,istrt,isize
      integer izero
      parameter(izero=0)
c
c start of garclds2
c
c  -- compute the minimum necessary size of the compressed output
c
      isize=0
      do i=1,n
        if(ka(i).ne.0) then
          k=ia(i)
          l=ja(k)
          ja(k)=-i
          ia(i)=l
          isize=isize+ka(i)
        else
          ia(i)=ibeg
        end if
      end do
c
c  -- compute distance parameter idist
c
      if(idist.eq.-1) then
c
c  -- -- original ldist.eq.'a'
c
        idist=max(0,(max_a-1-isize)/(2*(n+1)))
      elseif(idist.gt.0) then
c
c  -- -- original ldist.eq.'p'
c
        idist=min(max(0,(max_a-1-isize)/(2*(n+1))),idist)
      else
c
c  -- -- original case ldist.eq.'n'
c
        idist=0
      endif
c
c  -- the garbage collection
c
      istrt=ibeg
      i=ibeg
 300  continue
      k=ja(i)
c
c  -- -- find a start of some row/column
c
      if(k.lt.0) then
c
c  -- -- move the vertex
c
        do j=1,ka(-k)-1
          ja(istrt+j)=ja(i+j)
        end do
c
c  -- -- move the head and redefine the starting pointer
c
        ja(istrt)=ia(-k)
        ia(-k)=istrt
c
c  -- -- redefine starting pointers of both the source and target
c
        istrt=istrt+ka(-k)
        i=i+ka(-k)
c
c  -- -- clear some space (typically not large)
c
        do j=1,min(i-istrt,idist)
          ja(istrt)=0
          istrt=istrt+1
        end do
      else
c
c  -- -- move further
c
        i=i+1
      end if
c
c  -- -- end of coded loop
c
      if(i.lt.iend) go to 300
c
c  -- clear up to the end
c
      call iwset(iend-istrt+1,izero,ja(istrt),1)
c
c  -- redefine end
c
      iend=istrt-1
c
c  -- return
c
      return
c
c end of garclds2
c
      end

* (C) SPARSLAB MODULE NAME=DWSET
*
* PURPOSE:
*   SET ELEMENTS OF A DOUBLE PRECISION VECTOR TO A SCALAR.
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSLAB - TU - 18/5/1998.
*
* PARAMETERS:
*
      SUBROUTINE DWSET(N,ALPHA,DX,INCX)
C
C PARAMETERS
C
      INTEGER N,INCX
      DOUBLE PRECISION ALPHA,DX(*)
C
C INTERNALS
C
      INTEGER I,IX,M,MP1
C
C FUNCTIONS
C
      INTRINSIC MOD
C
C START OF DWSET
C
      IF(N.LE.0) RETURN
      IF(INCX.NE.1) THEN
C
C    -- CASE INCX.NE.1
C
        IX=1
        IF(INCX.LT.0) IX=(-N+1)*INCX+1
        DO I=1,N
          DX(IX)=ALPHA
          IX=IX+INCX
        END DO
        RETURN
      ELSE
C
C    -- CASE INCX.EQ.1
C
        M=MOD(N,7)
        IF(M.NE.0) THEN
          DO I=1,M
            DX(I)=ALPHA
          END DO
          IF(N.LT.7) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,7
          DX(I)=ALPHA
          DX(I+1)=ALPHA
          DX(I+2)=ALPHA
          DX(I+3)=ALPHA
          DX(I+4)=ALPHA
          DX(I+5)=ALPHA
          DX(I+6)=ALPHA
        END DO
        RETURN
      END IF
C
C END OF DWSET
C
      END

* (C) SPARSLAB MODULE NAME=IWSET
*
* PURPOSE:
*   SET ENTRIES OF AN INTEGER VECTOR TO IALPHA.
*   UNROLLED CODE.
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSELAB - TU - 14/1/1997
*
* PARAMETERS:
*
      SUBROUTINE IWSET(N,IALPHA,IX,INCX)
C
C GLOBALS
      INTEGER N,IALPHA,INCX
      INTEGER IX(*)
C
C LOCALS
      INTEGER I,IND,M,MP1
C
C FUNCTIONS
      INTRINSIC MOD
C
C START OF IWSET
C
      IF(N.LE.0) RETURN
      IF(INCX.NE.1) THEN
C
C INCX.NE.1
        IND=1
        IF(INCX.LT.0) IND=(-N+1)*INCX+1
        DO I=1,N
          IX(IND)=IALPHA
          IND=IND+INCX
        END DO
        RETURN
      ELSE
C
C INCX.EQ.1
        M=MOD(N,7)
        IF(M.NE.0) THEN
          DO I=1,M
            IX(I)=IALPHA
          END DO
          IF(N.LT.7) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,7
          IX(I)=IALPHA
          IX(I+1)=IALPHA
          IX(I+2)=IALPHA
          IX(I+3)=IALPHA
          IX(I+4)=IALPHA
          IX(I+5)=IALPHA
          IX(I+6)=IALPHA
        END DO
        RETURN
      END IF
C
C END OF IWSET
C
      RETURN
C
      END

* (C) SPARSLAB MODULE NAME=IWCOPY
*                DIRECTORY=SOURCE
*
* PURPOSE:
*   COPY AN INTEGER VECTOR.
*   UNROLLED CODE.
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSELAB - TU - 14/1/1998.
*   REWRITTEN FROM diniti.
*
* PARAMETERS:
*
      SUBROUTINE IWCOPY(N,IX,INCX,IY,INCY)
C
C GLOBALS
      INTEGER N,INCX,INCY,IX(*),IY(*)
C
C LOCALS
      INTEGER I,INTX,INTY,M,MP1
C
C FUNCTIONS
      INTRINSIC MOD
C
C START OF IWCOPY
C
      IF(N.LE.0) RETURN
      IF(.NOT.(INCX.EQ.1.AND.INCY.EQ.1)) THEN
C
C NONEQUAL INCREMENTS
        INTX=1
        INTY=1
        IF(INCX.LT.0) INTX=(-N+1)*INCX+1
        IF(INCY.LT.0) INTY=(-N+1)*INCY+1
        DO I=1,N
          IY(INTY)=IX(INTX)
          INTX=INTX+INCX
          INTY=INTY+INCY
        END DO
      ELSE
C
C INCREMENTS EQUAL TO 1
        M=MOD(N,7)
        IF(M.NE.0) THEN
          DO I=1,M
            IY(I)=IX(I)
          END DO
          IF(N.LT.7) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,7
          IY(I)=IX(I)
          IY(I+1)=IX(I+1)
          IY(I+2)=IX(I+2)
          IY(I+3)=IX(I+3)
          IY(I+4)=IX(I+4)
          IY(I+5)=IX(I+5)
          IY(I+6)=IX(I+6)
        END DO
      END IF
C
      RETURN
C
C END OF IWCOPY
C
      END

* (c) sparslab module name=realloci
*
* purpose:
*  reallocate an integer array.
*
* history:
*   original version for sparslab - tu - 16/04/2001.
*
c
c make_explicit_interface
c
      subroutine realloci(ja,arrsize,arrincr,new_arrsize,ierr)
c
c parameters
c
      integer arrsize,arrincr,new_arrsize,ierr
      integer, pointer, dimension(:) :: ja
c
c internals
c
      integer minincr
      integer, allocatable, dimension(:) :: ja_new
      double precision minincrfactor
      parameter(minincrfactor=0.1d0)
c
c start of realloci
c
      new_arrsize=arrsize+arrincr
      if(new_arrsize.le.arrsize) then
        minincr=minincrfactor*arrsize
        new_arrsize=arrsize+minincr
        if(new_arrsize.le.arrsize) new_arrsize=huge(arrsize)
      end if
      allocate(ja_new(arrsize),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','realloci',
     *  ' Allocation failure',1605,2)
      call iwcopy(arrsize,ja(1),1,ja_new(1),1)
      deallocate(ja)
      allocate(ja(new_arrsize),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','realloci',
     *  ' Allocation failure',1605,2)
      call iwcopy(arrsize,ja_new(1),1,ja(1),1)
      deallocate(ja_new)
c
c  -- return
c
      return
c
c end of realloci
c
      end

* (c) sparslab module name=reallocr
*
* purpose:
*  reallocate an integer array.
*
* history:
*   original version for sparslab - tu - 16/04/2001.
*
c
c make_explicit_interface
c
      subroutine reallocr(aa,arrsize,arrincr,new_arrsize,ierr)
c
c parameters
c
      integer arrsize,arrincr,new_arrsize,ierr
      double precision, pointer, dimension(:) :: aa
c
c internals
c
      integer minincr
      double precision, allocatable, dimension(:) :: aa_new
      double precision minincrfactor
      parameter(minincrfactor=0.1d0)
c
c start of reallocr
c
      new_arrsize=arrsize+arrincr
      if(new_arrsize.le.arrsize) then
        minincr=minincrfactor*arrsize
        new_arrsize=arrsize+minincr
        if(new_arrsize.le.arrsize) new_arrsize=huge(arrsize)
      end if
      allocate(aa_new(arrsize),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','reallocr',
     *  ' Allocation failure',1605,2)
      call dcopy(arrsize,aa(1),1,aa_new(1),1)
      deallocate(aa)
      allocate(aa(new_arrsize),stat=ierr)
      if(ierr.ne.0) call serr(0,2,'sparslab','reallocr',
     *  ' Allocation failure',1605,2)
      call dcopy(arrsize,aa_new(1),1,aa(1),1)
      deallocate(aa_new)
c
c  -- return
c
      return
c
c end of reallocr
c
      end

* SUBROUTINE UXVSRT                ALL SYSTEMS                88/12/01
C PORTABILITY : ALL SYSTEMS
C 88/12/01 TU : ORIGINAL VERSION
*
* PURPOSE :
* SHELLSORT
*
* PARAMETERS :
*  II  K LENGTH OF SORTED VECTOR.
*  IU  ARRAY(K) SORTED ARRAY.
*
      SUBROUTINE UXVSRT(K,ARRAY)
      INTEGER K
      INTEGER ARRAY(*)
      INTEGER IS,LA,LT,LS,LLS,I,J,JS,KHALF
C
C NOTHING TO BE SORTED
C
      IF(K.LE.1) GO TO 400
C
C SHELLSORT
C
      LS=131071
      KHALF=K/2
C      LS=2**(INT(LOG10(KHALF)/LOG10(2)))
      DO 300 LT=1,17
        IF(LS.GT.KHALF) THEN
          LS=LS/2
          GO TO 300
        END IF
        LLS=K-LS
        DO 200 I=1,LLS
          IS=I+LS
          LA=ARRAY(IS)
          J=I
          JS=IS
 100      IF(LA.GE.ARRAY(J)) THEN
            ARRAY(JS)=LA
            GO TO 200
          ELSE
            ARRAY(JS)=ARRAY(J)
            JS=J
            J=J-LS
          END IF
          IF(J.GE.1) GO TO 100
          ARRAY(JS)=LA
 200    CONTINUE
        LS=LS/2
 300  CONTINUE
 400  CONTINUE
      RETURN
      END

      SUBROUTINE UXVSR4(K,ARRAYI,ARRAYR)
      INTEGER ARRAYI(*)
      INTEGER IS,LA,LT,LS,LLS,I,J,K,JS,KHALF
      DOUBLE PRECISION LB
      DOUBLE PRECISION ARRAYR(*)
C
C NOTHING TO BE SORTED
C
      IF(K.LE.1) GO TO 400
C
C SHELLSORT
C
      LS=131071
      KHALF=K/2
      DO 300 LT=1,17
        IF(LS.GT.KHALF) THEN
          LS=LS/2
          GO TO 300
        END IF
        LLS=K-LS
        DO 200 I=1,LLS
          IS=I+LS
          LA=ARRAYI(IS)
          LB=ARRAYR(IS)
          J=I
          JS=IS
 100      IF(LA.GE.ARRAYI(J)) THEN
            ARRAYI(JS)=LA
            ARRAYR(JS)=LB
            GO TO 200
          ELSE
            ARRAYI(JS)=ARRAYI(J)
            ARRAYR(JS)=ARRAYR(J)
            JS=J
            J=J-LS
          END IF
          IF(J.GE.1) GO TO 100
          ARRAYI(JS)=LA
          ARRAYR(JS)=LB
 200    CONTINUE
        LS=LS/2
 300  CONTINUE
 400  CONTINUE
      RETURN
      END

* (C) SPARSLAB MODULE NAME=SERR
*
* PURPOSE:
*  ERROR OUTPUT ROUTINE.
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSLAB - TU - 05/10/1997.
*
* PARAMETERS :
*  II  EUNIT  ERROR UNIT NUMBER.
*  II  LUNIT  LOG UNIT NUMBER.
*  CI  LIBNAM  PACKAGE NAME.
*  CI  FILNAM  ROUTINE NAME.
*  CI  STRING  OUTPUT STRING.
*  II  ERR_NO  ERROR IDENTIFIER.
*  II  ERR_LEV  ERROR LEVEL:
*        =-1,0  FOR WARNINGS; =1 FOR RECOVERABLE ERRORS; =2 FOR FATAL
*        ERRORS.
*
      SUBROUTINE SERR(EUNIT,LUNIT,LIBNAM,FILNAM,STRING,ERR_NO,ERR_LEV)
C
C PARAMETERS
C
      INTEGER EUNIT,LUNIT,ERR_NO,ERR_LEV
      INTEGER LEN_OUTPUT_STRING_MAX
      PARAMETER(LEN_OUTPUT_STRING_MAX=100)
      CHARACTER*(LEN_OUTPUT_STRING_MAX) LIBNAM,FILNAM,STRING
C
C START OF SERR
C
      IF(EUNIT.EQ.0) THEN
        WRITE(*,1) '       Error No. ',ERR_NO
        WRITE(*,3) ' -- ',STRING
        WRITE(*,2) '         found in routine  ',FILNAM
        WRITE(*,2) '         called within package  ',LIBNAM
      ELSE
        WRITE(EUNIT,1) '       Error No. ',ERR_NO
        WRITE(EUNIT,3) ' -- ',STRING
        WRITE(EUNIT,2) '         found in routine  ',FILNAM
        WRITE(EUNIT,2) '         called within package  ',LIBNAM
      END IF
      IF(LUNIT.NE.EUNIT) THEN
        WRITE(LUNIT,1) '       Error No. ',ERR_NO
        WRITE(LUNIT,3) ' -- ',STRING
        WRITE(LUNIT,2) '         found in routine  ',FILNAM
        WRITE(LUNIT,2) '         called within package  ',LIBNAM
      END IF
      IF(ERR_LEV.GE.2) THEN
C
C    -- FATAL ERROR
C
        STOP
      END IF
C
C  -- FORMATS
C
 1    FORMAT((A),I7)
 2    FORMAT((A),2X,(A))
 3    FORMAT(9X,(A),(A))
C
C  -- RETURN
C
      RETURN
C
C END OF SERR
C
      END
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

      SUBROUTINE CS120(N,IA,JA)
C
C  PARAMETERS
C
      INTEGER N
      INTEGER IA(*),JA(*)
C
C INTERNALS
C
      INTEGER I
C
C START OF CS120
C
      DO I=1,IA(N+1)-1
        JA(I)=JA(I)-1
      END DO
      DO I=1,N+1
        IA(I)=IA(I)-1
      END DO
C
C  -- RETURN
C
      RETURN
C
C END OF CS120
C
      END

      SUBROUTINE CS021(N,IA,JA)
C
C  PARAMETERS
C
      INTEGER N
      INTEGER IA(*),JA(*)
C
C INTERNALS
C
      INTEGER I
C
C START OF CS021
C
      DO I=1,N+1
        IA(I)=IA(I)+1
      END DO
      DO I=1,IA(N+1)-1
        JA(I)=JA(I)+1
      END DO
C
C  -- RETURN
C
      RETURN
C
C END OF CS021
C
      END

* (c) sparslab module name=bkcch
*
* purpose:
*  back substitution for ic.
*  preconditioner in format cch (311).
*
* history:
*   original version for sparslab - tu - 06/07/2000.
*
* parameters :
*  ii  job  option.
*        job = 0 then  x = l^{-t} * l^{-1} * rhs
*        job > 0 then  x = l^{-1} * rhs
*        job < 0 then  x = l^{-t} * rhs
*
      subroutine bkcch(n,ia,ja,a,rhs,x,job)
c
c globals
c
      integer n
      integer ia(*),ja(*),job
      double precision a(*),x(*),rhs(*)
c
c internals
c
      integer i,j,k,jstrt,jstop
      double precision temp,zero
      parameter(zero=0.0d0)
c
c start of bkcch
c
c  -- x = rhs
c
      call dcopy(n,rhs,1,x,1)
c
c  -- x = l^{-1} * x
c
      if(job.ge.0) then
        do i=1,n
          jstrt=ia(i)
          x(i)=x(i)*a(jstrt)
          jstop=ia(i+1)-1
          temp=x(i)
          do j=jstrt+1,jstop
            k=ja(j)
            x(k)=x(k)-a(j)*temp
          end do
        end do
      end if
c
c  -- x = l^{-t} * x
c
      if(job.le.0) then
        do i=n,1,-1
          jstrt=ia(i)
          jstop=ia(i+1)-1
          temp=zero
          do j=jstrt+1,jstop
            k=ja(j)
            temp=temp-a(j)*x(k)
          end do
          x(i)=x(i)+temp
          x(i)=x(i)*a(jstrt)
        end do
      end if
c
c  -- return
c
      return
c
c end of bkcch
c
      end

* (C) SPARSLAB MODULE NAME=MVUCSR1_F90
*
* PURPOSE:
*   SPARSE MATRIX-VECTOR MULTIPLICATION Y=A*X.
*   CSR FORMAT.
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSLAB - TU - 01/05/1996.
*   MODIFIED - TU - 29/8/1997
*
* PARAMETERS:
*   II  M  NUMBER OF MATRIX ROWS.
*   MI  IA(M+1)/JA(IA(M+1)-1)/A(IA(M+1)-1)  INPUT MATRIX CSR OBJECT.
*   RI  X(N)  INPUT VECTOR.
*   RO  Y(M)  OUTPUT VECTOR Y=A*X.
*
C
C MAKE_EXPLICIT_INTERFACE
C
      SUBROUTINE MVUCSR1_F90(M,IA,JA,A,X,Y)
C
C PARAMETERS
C
      INTEGER M
      INTEGER, POINTER :: IA(:),JA(:)
      DOUBLE PRECISION, POINTER :: A(:)
      DOUBLE PRECISION X(*),Y(*)
C
C INTERNALS
C
      INTEGER I,J,K,JSTRT,JSTOP
      DOUBLE PRECISION TEMP,ZERO
      PARAMETER(ZERO=0.0D0)
C
C START OF MVUCSR1_F90
C
      JSTRT=1
      DO I=2,M+1
        JSTOP=IA(I)
        TEMP=ZERO
        DO J=JSTRT,JSTOP-1
          K=JA(J)
          TEMP=TEMP+A(J)*X(K)
        END DO
        Y(I-1)=TEMP
        JSTRT=JSTOP
      END DO
C
C  -- RETURN
C
      RETURN
C
C END OF MVUCSR1_F90
C
      END

* (C) SPARSLAB MODULE NAME=MVUCSR1
*
* PURPOSE:
*   SPARSE MATRIX-VECTOR MULTIPLICATION Y=A*X.
*   CSR FORMAT.
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSLAB - TU - 01/05/1996.
*   MODIFIED - TU - 29/8/1997
*
* PARAMETERS:
*   II  M  NUMBER OF MATRIX ROWS.
*   MI  IA(M+1)/JA(IA(M+1)-1)/A(IA(M+1)-1)  INPUT MATRIX CSR OBJECT.
*   RI  X(N)  INPUT VECTOR.
*   RO  Y(M)  OUTPUT VECTOR Y=A*X.
*
      SUBROUTINE MVUCSR1(M,IA,JA,A,X,Y)
C
C PARAMETERS
C
      INTEGER M
      INTEGER IA(*),JA(*)
      DOUBLE PRECISION A(*),X(*),Y(*)
C
C INTERNALS
C
      INTEGER I,J,K,JSTRT,JSTOP
      DOUBLE PRECISION TEMP,ZERO
      PARAMETER(ZERO=0.0D0)
C
C START OF MVUCSR1
C
      JSTRT=1
      DO I=2,M+1
        JSTOP=IA(I)
        TEMP=ZERO
        DO J=JSTRT,JSTOP-1
          K=JA(J)
          TEMP=TEMP+A(J)*X(K)
        END DO
        Y(I-1)=TEMP
        JSTRT=JSTOP
      END DO
C
C  -- RETURN
C
      RETURN
C
C END OF MVUCSR1
C
      END


* (C) SPARSLAB MODULE NAME=TVUCSR1_F90
*
*
* PURPOSE:
*   SPARSE MATRIX-VECTOR MULTIPLICATION Y=A^T*X; A IN CSR FORMAT.
*   (OR SPARSE MATRIX-VECTOR MULTIPLICATION Y=A*X; A IN CSC FORMAT.)
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSLAB - TU - 05-06-96
*   RENAMED FOR SPARSLAB - TU - 07-03-96
*   MODIFIED - TU - 29/8/1997
*
* PARAMETERS:
*   II  M  NUMBER OF MATRIX ROWS.
*   II  N  NUMBER OF MATRIX COLUMNS.
*   OI  IA(M+1)/JA(IA(M+1)-1)/A(IA(M+1)-1)  INPUT MATRIX CSR OBJECT.
*   RI  X(M)  INPUT VECTOR.
*   RO  Y(N)  OUTPUT VECTOR Y=A*X.
*
* SUBROUTINES AND FUNCTIONS
*   DWSET
C
C MAKE_EXPLICIT_INTERFACE
C
      SUBROUTINE TVUCSR1_F90(M,N,IA,JA,A,X,Y)
C
C PARAMETERS
C
      INTEGER M,N
      INTEGER, POINTER :: IA(:),JA(:)
      DOUBLE PRECISION, POINTER :: A(:)
      DOUBLE PRECISION X(*),Y(*)
C
C INTERNALS
C
      INTEGER I,J,K,JSTRT,JSTOP
      DOUBLE PRECISION TEMP,ZERO
      PARAMETER(ZERO=0.0D0)
C
C START OF TVUCSR1_F90
C
      CALL DWSET(N,ZERO,Y,1)
      JSTRT=1
      DO I=2,M+1
        JSTOP=IA(I)
        TEMP=X(I-1)
        DO J=JSTRT,JSTOP-1
          K=JA(J)
          Y(K)=Y(K)+A(J)*TEMP
        END DO
        JSTRT=JSTOP
      END DO
C
      RETURN
C
C END OF TVUCSR1_F90
C
      END

* (C) SPARSLAB MODULE NAME=TVUCSR1
*
*
* PURPOSE:
*   SPARSE MATRIX-VECTOR MULTIPLICATION Y=A^T*X; A IN CSR FORMAT.
*   (OR SPARSE MATRIX-VECTOR MULTIPLICATION Y=A*X; A IN CSC FORMAT.)
*
* HISTORY:
*   ORIGINAL VERSION FOR SPARSLAB - TU - 05-06-96
*   RENAMED FOR SPARSLAB - TU - 07-03-96
*   MODIFIED - TU - 29/8/1997
*
* PARAMETERS:
*   II  M  NUMBER OF MATRIX ROWS.
*   II  N  NUMBER OF MATRIX COLUMNS.
*   OI  IA(M+1)/JA(IA(M+1)-1)/A(IA(M+1)-1)  INPUT MATRIX CSR OBJECT.
*   RI  X(M)  INPUT VECTOR.
*   RO  Y(N)  OUTPUT VECTOR Y=A*X.
*
* SUBROUTINES AND FUNCTIONS
*   DWSET
C
      SUBROUTINE TVUCSR1(M,N,IA,JA,A,X,Y)
C
C PARAMETERS
C
      INTEGER M,N
      INTEGER IA(*),JA(*)
      DOUBLE PRECISION A(*),X(*),Y(*)
C
C LOCALS
C
      INTEGER I,J,K,JSTRT,JSTOP
      DOUBLE PRECISION TEMP,ZERO
      PARAMETER(ZERO=0.0D0)
C
C START OF TVUCSR1
C
      CALL DWSET(N,ZERO,Y,1)
      JSTRT=1
      DO I=2,M+1
        JSTOP=IA(I)
        TEMP=X(I-1)
        DO J=JSTRT,JSTOP-1
          K=JA(J)
          Y(K)=Y(K)+A(J)*TEMP
        END DO
        JSTRT=JSTOP
      END DO
C
      RETURN
C
C END OF TVUCSR1
C
      END
