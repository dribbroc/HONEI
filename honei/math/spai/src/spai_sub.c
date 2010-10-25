/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#define DECLAREEXTERNAL

#include "spai_sub.h"
#include "debug.h"

#include "local_lapack.h"

#include "config.h"

#define setbit(ptr,bit)  ((*ptr) |= (1 << (bit)))
#define btest(word,bit)  ((word) & (1 << (bit)))

/**********************************************************************/

int scalar_len
(index_set *I,
 matrix *A)
{
  int i,slen,global_index;

  for (i=0, slen=0; i<I->len; i++) {
    global_index = I->ptr[i];
    slen += A->block_sizes[global_index];
  }

  return(slen);
}

/**********************************************************************/

int precompute_column_square_inverses
(matrix *A)
{
  int i,j,b,b2,len;
  double *Aj_sq_inv;
  int *ipiv,n;
  double *work;
  int local_buf_size,buf_size,index,pe;
  double *local_si_buffer;
  int *buf_sizes,*buf_offsets,ierr;

  local_buf_size = 0;
  for (j=0; j<A->mnl; j++) {
    n = A->block_sizes[A->my_start_index+j];
    local_buf_size += (n*n);
  }

  local_si_buffer = (double *) mmalloc(local_buf_size*sizeof(double));
  ipiv = (int *) mmalloc(A->max_block_size*sizeof(int));
  work = (double *) mmalloc(A->max_block_size*sizeof(double));

  index = 0;
  for (j=0; j<A->mnl; j++) {
    len = A->lines->len[j];
    n = A->block_sizes[A->my_start_index+j];
    Aj_sq_inv = &local_si_buffer[index];

    if ((ierr = precompute_column_square_inverse1
      (A->lines->A[j],
       A->lines->ptrs[j],
       len,
       n,
       A,
       Aj_sq_inv,
       ipiv,
       work)) != 0)  return ierr;

    index += (n*n);

  }

  si_indices = (int *) mmalloc(A->n*sizeof(int));
  si_indices[0] = 0;

  for (i=1; i<A->n; i++) {
    n = A->block_sizes[i-1];
    si_indices[i] = si_indices[i-1] + (n*n);
  }

#ifdef MPI

  buf_sizes = (int *) mmalloc(A->numprocs*sizeof(int));
  buf_offsets = (int *) mmalloc(A->numprocs*sizeof(int));

  MPI_Barrier(A->comm);
  MPI_Allgather((void *) &local_buf_size,1,MPI_INT,
		(void *) buf_sizes,1,MPI_INT,
		A->comm);
  buf_size = buf_sizes[0];
  buf_offsets[0] = 0;
  for (pe=1; pe<A->numprocs; pe++) {
    buf_size += buf_sizes[pe];
    buf_offsets[pe] = buf_offsets[pe-1] + buf_sizes[pe-1];
  }
  si_buffer = (double *) mmalloc(buf_size*sizeof(double));
  MPI_Barrier(A->comm);

  MPI_Allgatherv((void *) local_si_buffer, local_buf_size, MPI_DOUBLE,
		 (void *) si_buffer, buf_sizes, buf_offsets,
		MPI_DOUBLE,A->comm);

  free(buf_sizes);
  free(buf_offsets);

#else

  /* For serial just just copy local buffer */
  si_buffer = (double *) mmalloc(local_buf_size*sizeof(double));

  for (i=0; i<local_buf_size; i++) {
    si_buffer[i] = local_si_buffer[i];
  }

#endif

  free(local_si_buffer);
  free(ipiv);
  free(work);
  return 0;
}

/**********************************************************************/

int precompute_column_square_inverse1
(double *Aj,
 int *ptrs,
 int len,
 int n,
 matrix *A,
 double *Aj_sq_inv,
 int *ipiv,
 double *work)
{
  int info;

  square_line(Aj,ptrs,len,
	      A,n,
	      Aj_sq_inv);
  if (n == 1) *Aj_sq_inv = 1/(*Aj_sq_inv); /* will do */
  else
   {
    F77_FUNC(dgetrf,DGETRF)(&n, &n, Aj_sq_inv, &n, ipiv, &info);

    if (info)
     {
      fprintf(stdout,"problem in DGETRF %d\n",info);
      return SPAI_BLAS_ERROR;
     }

    F77_FUNC(dgetri,DGETRI)(&n, Aj_sq_inv, &n, ipiv, work, &n, &info);

    if (info)
     {
      fprintf(stdout,"problem in DGETRI %d\n",info);
      fflush(stdout);
      return SPAI_BLAS_ERROR;
     }
   }
  return 0;
}

/**********************************************************************/
/*
   Compute {j | A(L,j) != 0} -> iset
*/

void getrows
(matrix *A,
 matrix *M,
 index_set *L,
 index_set *iset)
{
  int i,j,pos,bit,lenu,slenu,k,word;
  int row_pe,row_index;
  int len,rlen;
  int *buf;
  int *rbuf;
  double *vbuf;
  int ptr;

  /* Clear the bit vector */
  memset(bitvec,0,sizeof(unsigned int)*A->n);

  /* Look at every row of A indexed by L.
     The column index of every nonzero in those rows is or'd into bitvec. */
  lenu = slenu = 0;
  for (i=0; i<L->len; i++) {

    get_line(A, M, (L->ptr)[i], &len, &rlen, &buf, &rbuf, &vbuf);

    for (j=0; j<len; j++) {
      ptr = buf[j];
      pos = ptr >> log_nbits_of_int;
      bit = ptr % nbits_of_int;
      if (! btest( bitvec[pos], bit )) {
	setbit( &bitvec[pos], bit );
	iset->ptr[lenu] = ptr;
	lenu++;
	slenu += A->block_sizes[ptr];
      }
    }
  }

  iset->len  =  lenu;
  iset->slen = slenu;

}

/**********************************************************************/
/*
   Compute {i | A(i,L) != 0} -> iset
*/

void getcols
(matrix *A,
 matrix *M,
 index_set *L,
 index_set *iset)
{
  int i,j,pos,bit,lenu,slenu,k,word;
  int row_pe,row_index;
  int len,rlen;
  int *buf;
  int *rbuf;
  double *vbuf;
  int ptr;

  /* Clear the bit vectors */
  memset(bitvec,0,sizeof(unsigned int)*A->n);

  /* Look at every column of A indexed by L.
     The row index of every nonzero in those columns is or'd into bitvec. */
  lenu = slenu = 0;
  for (i=0; i<L->len; i++) {

    get_line(A, M, (L->ptr)[i], &len, &rlen, &buf, &rbuf, &vbuf);

    for (j=0; j<rlen; j++) {
      ptr = rbuf[j];
      pos = ptr >> log_nbits_of_int;
      bit = ptr % nbits_of_int;
      if (! btest( bitvec[pos], bit )) {
	setbit( &bitvec[pos], bit );
	iset->ptr[lenu] = ptr;
	lenu++;
	slenu += A->block_sizes[ptr];
      }
    }
  }

  iset->len  = lenu;
  iset->slen = slenu;

}

/**********************************************************************/
/* Construct the incremental part of the full matrix */

/* dimensions: Ahat(m,n), where
   m is the scalar length of I.
   n is the scalar length of J_tilde.
*/

void full_matrix
(matrix *A,
 matrix *M,
 int max_dim,
 double *Ahat)
{
  int nb,mb;
  int row_start,col_start,buf_start;
  int row,col,col2;
  int i,j,k,m,n;
  int len,rlen;
  int jb,ib;
  int line;
  int *buf;
  int *rbuf;
  double *vbuf;
  int ptr;
  int index,pe,global_index;

  for (m=0, k=0; k<I->len; k++)
   {
    global_index = I->ptr[k];
    m += A->block_sizes[global_index];
   }

  for (n=0, k=0; k<J_tilde->len; k++)
   {
    global_index = J_tilde->ptr[k];
    n += A->block_sizes[global_index];
   }

  if (n)
   {
    if ((n*max_dim + m) > Ahat_size)
     {
      printf("error in full_matrix: exceeded Ahat size\n");
      exit(1);
     }

    /* initialize */
    for (col=0, col2=0;
	 col<n;
	 col++, col2+=max_dim)
      for (row=0; row<m; row++)
        Ahat[col2+row]=0.0;

    for (col_start=0, j=0;
	 j<J_tilde->len;
	 j++, col_start += (nb*max_dim))
     {
      get_line(A, M, (J_tilde->ptr)[j], &len, &rlen, &buf, &rbuf, &vbuf);

      global_index = J_tilde->ptr[j];
      nb = A->block_sizes[global_index];

      for (k=0, buf_start=0; k<len; k++)
       {
	line = buf[k];

	/* Skip over unneeded entries */
	for (i=0, row_start=0;
	     (!(line == I->ptr[i]) && (i < I->len));
	     i++)
         {
	  global_index = I->ptr[i];
	  mb = A->block_sizes[global_index];
	  row_start += mb;
	 }

	global_index = I->ptr[i];
	mb = A->block_sizes[global_index];

	for (jb=0,
	       col=col_start;
	     jb<nb;
	     jb++, col+=max_dim)
         {
	  for (ib=0, row=row_start;
	       ib<mb;
	       ib++, row++)
	    Ahat[col + row]
	      = vbuf[buf_start + jb*mb + ib];
	 }

	global_index = buf[k];
	buf_start += (nb*A->block_sizes[global_index]);

       }

     }
   }

}

/**********************************************************************/
/* for debugging */

void write_full_matrix(FILE *fptr, double *Ahat, int m, int n, int max_dim)
{
  int i,j,jj,bmax;

  for (i=0; i<m; i++) {
    for (j=0, jj=0; j<n; j++, jj+=max_dim) {
      fprintf(fptr,"%12.4le ",Ahat[jj+i]);
    }
    fprintf(fptr,";\n");
  }
}

/**********************************************************************/

int augment_sparsity
(matrix *A,
 matrix *M,
 int col,
 int maxapi,
 double resnorm)
{
  double mean;

  /* copy I (residual indices) into isres */
  copyvv(I,isres);
  orderv(resb,
	 isres->ptr,
	 isres->len,
	 A,
	 A->block_sizes[col],
	 A->max_block_size);

  /* compute union of all nonzero entries in I -> is */

  getcols(A,M, I,is);

  deleter(J,is,A);

  if (! is->len) {
    if (message) {
      fprintf(message,"No new candidates found\n");
      return(0);
    }
  }

  /* compute 1D min to get rs */

  mean = calcrs(A,M,col,resnorm);
  if (debug) {
    fprintf(fptr_dbg,"  in augment_sparsity mean=%le\n",mean);
    fflush(fptr_dbg);
    }


  /* determine new nonzeros -> J_tilde */
  select_new_nonzeros(mean,A);

  return(1);
}

/**********************************************************************/
/*
  extracts from (is,rs) all elements smaller than val -> J_tilde
  if more than maxnew found, take smallest maxnew elements
*/

/* ???? work on this later */

void select_new_nonzeros_insertion
(double val,
 matrix *A)
{
  double dummy = 1.0e100;  /* Should be infinity */
  double curr;
  int i,k,icurr,new_slen;

  void insertion_sort(int, double *, int *);

  J_tilde->len = J_tilde->slen = 0;

  /* Extract all elements <= val => rz */
  new_slen = 0;
  for (i=0; i<is->len; i++)
    if (rs[i] <= val) {
      (J_tilde->ptr)[J_tilde->len] = (is->ptr)[i];
      rz[J_tilde->len++] = rs[i];
      new_slen += block_size(is->ptr[i],A);
    }


  /* ?????? */



}

void insertion_sort(int n, double *keys, int *indices)
{
  int j,i;
  double a,index;

  for (j=1; j<n; j++) {
    a = keys[j];
    index = indices[j];
    i = j-1;
    while ((i >= 0) && keys[i] > a) {
      keys[i+1] = keys[i];
      indices[i+1]=indices[i];
      i--;
    }
    keys[i+1] = a;
    indices[i+1] = index;
  }
}

/**********************************************************************/
/*
  extracts from (is,rs) all elements smaller than val -> J_tilde
  if more than maxnew found, take smallest maxnew elements
*/

void select_new_nonzeros
(double val,
 matrix *A)
{
  double dummy = 1.0e100;  /* Should be infinity */
  double curr;
  int i,k,icurr,new_slen;

  J_tilde->len = J_tilde->slen = 0;

  /* ????
  if (debug) {
    fprintf(fptr_dbg,"entered select_new_nonzeros\n");
    fprintf(fptr_dbg,"val=%le\n",val);
    for (i=0; i<is->len; i++) {
      if (rs[i] <= val) fprintf(fptr_dbg,"  *** ");
      else fprintf(fptr_dbg,"      ");
      fprintf(fptr_dbg,"%le %d\n",rs[i],is->ptr[i]);
    }
  }
  */

  /* Extract all elements <= val => rz */
  new_slen = 0;
  for (i=0; i<is->len; i++)
    if (rs[i] <= val) {
      (J_tilde->ptr)[J_tilde->len] = (is->ptr)[i];
      rz[J_tilde->len++] = rs[i];
      new_slen += block_size(is->ptr[i],A);
    }

  /* ????
  if (debug) {
    fprintf(fptr_dbg,"   all elements <= val => rz\n");
    for (i=0; i<new_slen; i++) {
      fprintf(fptr_dbg,"%le %d\n",rz[i],J_tilde->ptr[i]);
    }
  }
  */

  if (J_tilde->len > maxnew) {

    /* copy rz -> rs */
    for (i=0; i<J_tilde->len; i++) {
      (is->ptr)[i] = (J_tilde->ptr)[i];
      rs[i] = rz[i];
    }
    is->len = J_tilde->len;
    is->slen = J_tilde->slen;

    new_slen = 0;
    for (k=0; k<maxnew; k++) {
      curr = rs[0];
      icurr = 0;
      for (i=1; i<J_tilde->len; i++) {
        if (rs[i] < curr) {
	  /* ?????
	  if (debug) {
	    fprintf(fptr_dbg,"     i=%d rs[i]=%le curr=%le\n",i,rs[i],curr);
	  }
	  */
          curr = rs[i];
          icurr = i;
        }
      }
      /* ????
      if (debug) {
	fprintf(fptr_dbg,"  k=%d curr=%le icurr=%d\n",k,curr,icurr);
      }
      */
      /* smallest now in icurr */
      (J_tilde->ptr)[k] = (is->ptr)[icurr];
      new_slen += block_size(J_tilde->ptr[k],A);
      rs[icurr] = dummy;
    }
    J_tilde->len = maxnew;
  }

  J_tilde->slen = new_slen;

}

/**********************************************************************/
/* copies is_in into is_out */

void copyvv
(index_set *is_in,
 index_set *is_out)
{
  memcpy((void *) is_out->ptr,
         (void *) is_in->ptr,
         (is_in->len)*sizeof(int));
  is_out->len  = is_in->len;
  is_out->slen = is_in->slen;
}

/**********************************************************************/
/* appends isb at end of isa */
/* Returns 1 if success, 0 if failure */

int append
(index_set *isa,
 index_set *isb)
{
  if (isa->len + isb->len > max_dim) {
    if (message) {
      fprintf(message,
	      "Not enough memory in work arrays.  Increase max_dim.\n");
      return(0);
    }
  }

  memcpy((void *) &(isa->ptr)[isa->len],
         (void *) isb->ptr,
         (isb->len)*sizeof(int));
  isa->len += isb->len;
  isa->slen += isb->slen;

  return(1);
}

/**********************************************************************/
/* computes v^T * w */

double innprod(double *v, double *w, int n)
{
  double result;
  int i;

  for (result=0.0, i=0; i<n; i++)
    result += v[i]*w[i];

  return(result);
}

/**********************************************************************/
/* v is an mxn matrix */

double frobenius_norm
(double *v,
 int m,
 int n)
{
  return(sqrt(innprod(v,v,m*n)));
}

/**********************************************************************/
/*
  computes u \ r --> u
  no ordering is assumed
*/

void deleter
(index_set *r,
 index_set *u,
 matrix *A)
{
  int iu,ir,new_len,new_slen;
  int found;
  int iaddress;

  if (r->len) {
    new_len = new_slen = 0;
    for (iu=0; iu<u->len; iu++) {
      iaddress = *((int *) &(u->ptr)[iu]);
      for (ir=0, found=0; ir<r->len; ir++) {
	if (iaddress == (*((int *) &(r->ptr)[ir]))) {
	  found = 1;
	  break;
	}
      }
      if (! found) {
	u->ptr[new_len++] = u->ptr[iu];
	new_slen += block_size(u->ptr[iu],A);
      }
    }
    u->len = new_len;
    u->slen = new_slen;
  }
}

/**********************************************************************/

double calcrs
(matrix *A,
 matrix *M,
 int k,
 double resnrm)
{
  int i,ii,q,ia,ir,colj,colk;
  int rowtmp,restmp;
  double sumrs,resnrm_sq;
  int len,rlen;
  int *buf,*rbuf,isptr;
  double *vbuf;
  int bk,bj,bla,blr;
  double term,*Aj_sq_inv;
  int r_start_block,a_start_block;

  double *a,*b;
  int i1,ii1,j1,jj1,k1;
  int bs;

  bs = A->bs;

  resnrm_sq = resnrm*resnrm;
  sumrs = 0.0;
  bk = A->block_sizes[k];

  for (i=0; i<is->len; i++) {

    isptr = is->ptr[i];

    get_line(A, M, isptr, &len, &rlen, &buf, &rbuf, &vbuf);

    bj = block_size(isptr,A);
    zero_block(sum_block,bj,bk);
    Aj_sq_inv = &si_buffer[si_indices[isptr]];
    ia = ir = r_start_block = a_start_block = 0;
    /* Compute (Aj^T)r => sum_block (bj x bk) */
    do {

/*      rowtmp = (* ((int *) &buf[ia]));
	restmp = (* ((int *) &(isres->ptr)[ir])); */
      rowtmp = buf[ia];
      restmp = isres->ptr[ir];
      bla = block_size(buf[ia],A);
      blr = block_size(isres->ptr[ir],A);

      if (rowtmp == restmp) {

	if (bla != blr) {
	  printf("error in calcrs\n");
	  exit(1);
	}

	/* vbuf[*] is (bla x bj)
	   resb[*] is (blr x bk)
	   plus: bla == blr */
	/* Avoid function call for bs == 1 case */
	if (bs == 1)
	  sum_block[0] +=
	    (vbuf[a_start_block]*resb[r_start_block]);
	else
	  mult_blocks_TN
	    (&vbuf[a_start_block],&resb[r_start_block],
	     bj,bk,bla,sum_block);

	a_start_block += (bla*bj); ia++;
	r_start_block += (blr*bk); ir++;

      }

      else if (rowtmp > restmp) {
	r_start_block += (blr*bk);
	ir++;
      }
      else {
	a_start_block += (bla*bj);
	ia++;
      }
    }
    while ((ir < isres->len) && (ia < len));

    /* Avoid function calls for bs == 1 case */
    if (bs == 1) minus_Bj[0] = Aj_sq_inv[0]*sum_block[0];
    else {
      zero_block(minus_Bj,bj,bk);
      mult_blocks(Aj_sq_inv,sum_block,bj,bk,bj,minus_Bj);
    }
    term = 0.0;

    for (q=0, colj=0, colk=0;
	 q<bk;
	 q++, colj+=bj, colk+=bk) {
      for (ii=0; ii<bj; ii++) {
	term += minus_Bj[ii+colj]*sum_block[ii+colk];
      }
    }

    sumrs += (rs[i] = resnrm_sq - term);
  }

  /* return the mean */
  return( sumrs/(is->len) );

}

/**********************************************************************/
/* for debugging */

void write_line
(FILE *fptr,
 matrix *A,
 int ptr,
 int len,
 int rlen,
 int *buf,
 int *rbuf,
 double *vbuf)
{
  int index,pe,my_global_index,global_index,m,n,i,block_start;

  my_global_index = ptr;
  m = A->block_sizes[my_global_index];
  fprintf(fptr,"***********line = %d, block_width=%d\n",
	  my_global_index,m);

  fprintf(fptr,"buf and vbuf: len=%d\n",len);
  block_start = 0;
  for (i=0; i<len; i++) {
    global_index = buf[i];
    fprintf(fptr,"  i=%d global_index=%d\n",i,global_index);
    n = A->block_sizes[global_index];
    write_block(fptr,&vbuf[block_start],m,n);
    block_start += m*n;
  }

  fprintf(fptr,"rbuf: rlen=%d\n",rlen);
  for (i=0; i<len; i++) {
    global_index = rbuf[i];
    fprintf(fptr,"  i=%d global_index=%d\n",i,global_index);
  }

}

