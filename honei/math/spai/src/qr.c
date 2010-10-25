/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "basics.h"
#include "qr.h"
#include "debug.h"

#include "config.h"

#include "local_lapack.h"

char *Tchar = "T";
char *Nchar = "N";
char *Uchar = "U";
char *Lchar = "L";

#ifdef T3D
_fcd Tchar_fcd;
_fcd Nchar_fcd;
_fcd Uchar_fcd;
_fcd Lchar_fcd;
#endif

/**********************************************************************/

int qr
(matrix *A,
 int col,
 int nbq,
 int dimr)
{
  int i,j,jj,jjj,offset,info,indx;
  int nrhs,unblocked_len;
  double *Q;
  int nline_A,block_width,jj_start;
  int index,pe,global_index;

  nrhs = block_width = A->block_sizes[col];
  nline_A = scalar_len(J_tilde,A);

  /* Make sure R is big enough */
  if ((maxapi+1)*(dimr+nline_A) > R_size) {
    R_size = (maxapi+1)*(dimr+nline_A);
    R = new_double_array(R,R_size,"R");
  }

#ifdef T3D
  Tchar_fcd = _cptofcd(Tchar,1);
  Nchar_fcd = _cptofcd(Nchar,1);
  Uchar_fcd = _cptofcd(Uchar,1);
  Lchar_fcd = _cptofcd(Lchar,1);
#endif

  if (nbq >= 1) {

    multq(Tchar, nbq-1, Ahat,nline_A, block_width);

    for (j=0,
	   jj=maxapi*dimr,
	   jjj=0;
	 j<nline_A;
	 j++,
	   jj+=maxapi,
	   jjj+=max_dim)
      for (i=0; i<dimr; i++)
	R[i + jj] = Ahat[i+jjj];
  }

  /*C ------- compute QR of A (premultiplied by Q**T if necessary) */
  n1[nbq+1] = n1[nbq] + nline_A;
  offset = n1[nbq]-block_width;

  /* Make sure Z is big enough */
  if ((nline_A*max_dim + n2[nbq]) > Z_size) {
    Z_size = nline_A*max_dim + n2[nbq];
    Z = new_double_array(Z,Z_size,"Z");
  }

  for (j=0,
	 jj=0;
       j<nline_A;
       j++,
	 jj+=max_dim)
    for (i=0; i<n2[nbq]; i++)
      Z[i + jj] = Ahat[i+offset + jj];

  if (max_dim < n2[nbq]) {
    printf(".. illegal value in dgeqrf: \n");
    printf("   n2[nbq]=%d must be <= max_dim=%d\n",
	   n2[nbq],max_dim);
    printf(".. Try increasing the -mb parameter\n");
    return SPAI_BLAS_ERROR;
  }

  TAU_ptr[nbq+1] = TAU_ptr[nbq] + MIN(n2[nbq],nline_A);

  if (TAU_size <= TAU_ptr[nbq+1]) {
    TAU_size = 2*TAU_ptr[nbq+1];
    TAU    = new_double_array(TAU,TAU_size,"TAU");
  }

  F77_FUNC(dgeqrf,DGEQRF)(&n2[nbq],&nline_A,Z,&max_dim,&TAU[TAU_ptr[nbq]],rw,&max_dim,&info);

/* #if defined(T3D) */
/*   SGEQRF(&n2[nbq],&nline_A,Z,&max_dim,&TAU[TAU_ptr[nbq]],rw,&max_dim,&info); */
/* #elif defined(SP2) */
/*   dgeqrf(&n2[nbq],&nline_A,Z,&max_dim,&TAU[TAU_ptr[nbq]],rw,&max_dim,&info); */
/* #else */
/*   dgeqrf_(&n2[nbq],&nline_A,Z,&max_dim,&TAU[TAU_ptr[nbq]],rw,&max_dim,&info); */
/* #endif */

  if (info) {
    fprintf(stdout,"problem in DGEQRF %d",info);
    return SPAI_BLAS_ERROR;
  }

  /* C ------- copy Q and R into arrays...  */

  for (j=0,
	 jj=maxapi*dimr,
	 jjj=0;
       j<nline_A;
       j++,
	 jj+=maxapi,
	 jjj+=max_dim)
    for (i=0; i<=j; i++)
      R[dimr+i + jj] = Z[i + jjj];

  Qlist[nbq] = (double *) mmalloc(nline_A*n2[nbq]*sizeof(double));
  Q = Qlist[nbq];
  for (j=0,
	 jj=0,
	 jjj=0;
       j<nline_A;
       j++,
	 jj+=n2[nbq],
	 jjj+=max_dim)
    for (i=j+1; i<n2[nbq]; i++)
      Q[i + jj] = Z[i + jjj];

  dimr += nline_A;

  /* set X, like ek */
  /* !!!! do the following more efficiently */
  unblocked_len = scalar_len(I,A);
  fill_zeros(nrhs,max_dim,unblocked_len,x);
  indx = seek_ptr(A,I,col,&jj_start);
  if (indx < 0) {

    /* not found in I */
    (I->ptr)[(I->len)++] = col;
    fill_zeros(block_width,max_dim,jj_start,res);

    for (j=0,
	   jj=jj_start;
	 j<nrhs;
	 j++,
	   jj+=max_dim)
      for (i=0; i<block_width; i++) {
	if (i == j) res[jj + i] = 1.0;
	else res[jj + i] = 0.0;
      }
  }

  else {

    /* k is present in I, at indx */
    for (j=0,
	   jj=jj_start;
	 j<nrhs;
	 j++,
	   jj+=max_dim)
      for (i=0; i<block_width; i++) {
	if (i == j) x[jj + i] = 1.0;
	else x[jj + i] = 0.0;
      }

    multq(Tchar, nbq, x,nrhs, block_width);

    /* compute residual */
    fill_zeros(block_width,max_dim,unblocked_len,res);

    for (j=0,
	   jj=0;
	 j<nrhs;
	 j++,
	   jj+=max_dim) {
      for (i=dimr; i<unblocked_len; i++) {
	res[jj + i] = -x[jj + i];
      }
    }

    multq(Nchar, nbq,res,nrhs, block_width);

    /* solve upper triangular system, solution in x */

    F77_FUNC(dtrtrs,DTRTRS)(Uchar,Nchar,Nchar,
           &dimr,&nrhs,R,&maxapi,x,&max_dim,&info);

/* #if defined(T3D) */
/*     STRTRS(Uchar_fcd,Nchar_fcd,Nchar_fcd, */
/*            &dimr,&nrhs,R,&maxapi,x,&max_dim,&info); */
/* #elif defined(SP2) */
/*     dtrtrs(Uchar,Nchar,Nchar, */
/*            &dimr,&nrhs,R,&maxapi,x,&max_dim,&info); */
/* #else */
/*     dtrtrs_(Uchar,Nchar,Nchar, */
/*            &dimr,&nrhs,R,&maxapi,x,&max_dim,&info); */
/* #endif */

    if (info) {
      fprintf(stdout,"problems in DTRTRS %d",info);
      return SPAI_BLAS_ERROR;
    }
  }
  return 0;
}

/**********************************************************************/
/* for debugging */

void write_unblocked(double *v, int max_dim, int nrhs, int n)
{
  int i,j,jj;

  for (i=0; i<n; i++) {
    for (j=0,
	   jj=0;
	 j<nrhs;
	 j++,
	   jj+=max_dim) {
      printf("%le ",v[jj+i]);
    }
    printf("\n");
  }
}

/**********************************************************************/

void multq
(char *trans,
 int nbq,
 double *A,
 int nline_A,
 int block_width)
{
  int iq,offset,j,jj,i,dif,info;

  if (*trans == *Tchar)

    for (iq = 0; iq <= nbq; iq++) {
      offset = n1[iq]-block_width;

      /* Make sure Z is big enough */
      if ((nline_A*max_dim + n2[nbq]) > Z_size) {
	Z_size = nline_A*max_dim + n2[nbq];
	Z = new_double_array(Z,Z_size,"Z");
      }

      /* copy part of A that needs to be multiplied with Q(IQ) */
      for (j=0,
	     jj=0;
	   j<nline_A;
	   j++,
	     jj+=max_dim)
	for (i=0; i<n2[iq]; i++)
	  Z[i + jj] = A[i+offset + jj];

      /* perform multiplication */
      dif = n1[iq+1] - n1[iq];

/* #ifdef T3D */
/*       SORMQR(Lchar_fcd,Tchar_fcd,&n2[iq],&nline_A,&dif, */
/* 	     Qlist[iq], */
/* 	     &n2[iq], */
/* 	     &TAU[TAU_ptr[iq]], */
/* 	     Z,&max_dim,rw,&max_dim,&info); */

/* #elif ORIGIN */
/*       dormqr_(Lchar,Tchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	     Z,&max_dim,rw,&max_dim,&info); */

/* #elif INDY */
/*       dormqr_(Lchar,Tchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	      Z,&max_dim,rw,&max_dim,&info); */

/* #elif SP2 */
/*       dormqr(Lchar,Tchar,&n2[iq],&nline_A,&dif, */
/* 	     Qlist[iq], */
/* 	     &n2[iq], */
/* 	     &TAU[TAU_ptr[iq]], */
/* 	     Z,&max_dim,rw,&max_dim,&info); */

/* #elif Darwin */
/*       dormqr_(Lchar,Tchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	      Z,&max_dim,rw,&max_dim,&info); */
/* #else */
/*       dormqr_(Lchar,Tchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	      Z,&max_dim,rw,&max_dim,&info); */
/* #endif */

      F77_FUNC (dormqr,DORMQR)
	      (Lchar,Tchar,&n2[iq],&nline_A,&dif,
	       Qlist[iq],
	       &n2[iq],
	       &TAU[TAU_ptr[iq]],
	       Z,&max_dim,rw,&max_dim,&info);

      if (info) fprintf(stdout,"problem in DORMQR\n");

      for (j=0, jj=0;
	   j<nline_A;
	   j++, jj+=max_dim)
	for (i=0; i<n2[iq]; i++) {
	  A[i+offset + jj] = Z[i + jj];
	}

    }

  else {

    for (iq = nbq; iq >= 0; iq--) {
      offset = n1[iq]-block_width;

      /* Make sure Z is big enough */
      if ((nline_A*max_dim + n2[nbq]) > Z_size) {
	Z_size = nline_A*max_dim + n2[nbq];
	Z = new_double_array(Z,Z_size,"Z");
      }

      /* copy part of A that needs to be multiplied with Q(IQ) */
      for (j=0,
	     jj=0;
	   j<nline_A;
	   j++,
	     jj+=max_dim)
	for (i=0; i<n2[iq]; i++)
	  Z[i + jj] = A[i+offset + jj];

      /* perform multiplication */
      dif = n1[iq+1] - n1[iq];

      F77_FUNC(dormqr,DORMQR)(Lchar,Nchar,&n2[iq],&nline_A,&dif,
	     Qlist[iq],
	     &n2[iq],
	     &TAU[TAU_ptr[iq]],
	     Z,&max_dim,rw,&max_dim,&info);

/* #ifdef T3D */
/*       SORMQR(Lchar_fcd,Nchar_fcd,&n2[iq],&nline_A,&dif, */
/* 	     Qlist[iq], */
/* 	     &n2[iq], */
/* 	     &TAU[TAU_ptr[iq]], */
/* 	     Z,&max_dim,rw,&max_dim,&info); */

/* #elif ORIGIN */
/*       dormqr_(Lchar,Nchar,&n2[iq],&nline_A,&dif, */
/* 	     Qlist[iq], */
/* 	     &n2[iq], */
/* 	     &TAU[TAU_ptr[iq]], */
/* 	     Z,&max_dim,rw,&max_dim,&info); */

/* #elif INDY */
/*       dormqr_(Lchar,Nchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	      Z,&max_dim,rw,&max_dim,&info); */

/* #elif SP2 */
/*       dormqr(Lchar,Nchar,&n2[iq],&nline_A,&dif, */
/* 	     Qlist[iq], */
/* 	     &n2[iq], */
/* 	     &TAU[TAU_ptr[iq]], */
/* 	     Z,&max_dim,rw,&max_dim,&info); */

/* #elif Darwin */
/*       dormqr_(Lchar,Nchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	      Z,&max_dim,rw,&max_dim,&info); */

/* #else */
/*       dormqr_(Lchar,Nchar,&n2[iq],&nline_A,&dif, */
/* 	      Qlist[iq], */
/* 	      &n2[iq], */
/* 	      &TAU[TAU_ptr[iq]], */
/* 	      Z,&max_dim,rw,&max_dim,&info); */
/* #endif */

      if (info) fprintf(stdout,"problem in DORMQR\n");

      /* copy result back into A */
      for (j=0,
	     jj=0;
	   j<nline_A;
	   j++,
	     jj+=max_dim)
	for (i=0; i<n2[iq]; i++)
	  A[i+offset + jj] = Z[i + jj];

    }

  }

}

/**********************************************************************/

int seek_ptr
(matrix *A,
 struct index_set *s,
 int ptr,
 int *scalar_start_ptr)
{
  int i,global_index;

  *scalar_start_ptr = 0;
  for (i=0; i<s->len; i++) {
    if (s->ptr[i] == ptr)
      return(i);
    global_index = s->ptr[i];
    *scalar_start_ptr += A->block_sizes[global_index];
  }


  return(-1);
}

/**********************************************************************/

void fill_zeros(int nrhs, int max_dim, int len, double *v)
{
  int i,j,jj;

  for (j=0, jj=0; j<nrhs; j++, jj+=max_dim)
    for (i=0; i<len; i++)
      v[jj+i] = 0.0;
}

/**********************************************************************/
