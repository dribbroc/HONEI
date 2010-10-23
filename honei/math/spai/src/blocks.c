/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "blocks.h"
#include "debug.h"

/* (m x n) blocks are stored in column-major order */

/**********************************************************************/

void write_block
(FILE *fptr,
 double *a,
 int m,
 int n)
{
  int i,j;

  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      fprintf(fptr,"%le ",a[i+m*j]);
    }
    fprintf(fptr,";\n");
  }
}

/**********************************************************************/
/* given a vector of blocks v, computes c = v*v^T */

void square_line
(double *v,
 int *ptrs,
 int len,
 matrix *A,
 int n,            /* block width */
 double *c)
{
  int i,m,start_block,global_index;

  zero_block(c,n,n);

  start_block = 0;
  for (i=0; i<len; i++) {

    global_index = ptrs[i];
    m = A->block_sizes[global_index];

    mult_blocks_TN(&v[start_block],&v[start_block],n,n,m,c);

    start_block += (m*n);
  }

}

/**********************************************************************/

void zero_block
(double *a,
 int m,
 int n)
{
  int i;

  for (i=0; i<m*n; i++)
    a[i] = 0.0;
}

/**********************************************************************/
/* Given blocks a and b, compute c = c + a*b */

void mult_blocks
(double *a,  /* m x k */
 double *b,  /* k x n */
 int m,
 int n,
 int k,
 double *c)  /* m x n */
{
  double one=1.0;
  char *Tchar = "T";
  char *Nchar = "N";
  int i1,j1,jj1,jjj1,k1,kk1;

  /* The following loop is (almost) equivalent to:
  dgemm_(Nchar,Nchar, &m,&n,&k,
	 &one, a, &m, b, &k, &one, c, &m);
  It doesn't treat the array dimensions in exactly the same way.
  */

  for (i1=0; i1<m; i1++)
    for (j1=0, jj1=0, jjj1=0; j1<n; j1++, jj1+=k, jjj1+=m)
      for (k1=0, kk1=0; k1<k; k1++, kk1+=m) {
	c[i1+jjj1] += a[i1+kk1]*b[k1+jj1];
      }

}


/**********************************************************************/
/* Given blocks a and b, compute c = c + (a^T)*b */

void mult_blocks_TN
(double *a,  /* k x m */
 double *b,  /* k x n */
 int m,
 int n,
 int k,
 double *c)  /* m x n */
{
  double one=1.0;
  char *Tchar = "T";
  char *Nchar = "N";
  int i1,ii1,j1,jj1,jjj1,k1;

  /* The following loop is (almost) equivalent to:
  dgemm_(Tchar,Nchar, &m,&n,&k,
	 &one, a, &m, b, &k, &one, c, &m);
  It doesn't treat the array dimensions in exactly the same way.
  */

  for (i1=0, ii1=0; i1<m; i1++, ii1+=k)
    for (j1=0, jj1=0, jjj1=0; j1<n; j1++, jj1+=k, jjj1+=m)
      for (k1=0; k1<k; k1++) {
	c[i1+jjj1] += a[k1+ii1]*b[k1+jj1];
      }

}

/**********************************************************************/
/* Converts a scalar ordering to a block ordering. */

void convert_to_block
(double *x,
 double *xb,
 int col,
 int *ptrs,
 matrix *A,
 int max,
 int len)
{
  int start_x,start_xb,m,n,i,j,k,jb,js;
  int global_index;

  start_x = start_xb = 0;
  n = A->block_sizes[col];

  for (k=0; k<len; k++) {

    global_index = ptrs[k];
    m = A->block_sizes[global_index];

    for (j=0, jb=0, js=0;
	 j<n;
	 j++, jb+=m, js+=max) {

      for (i=0; i<m; i++) {

	xb[start_xb + jb + i] = x[start_x + js + i];

      }
    }
    start_x += m;
    start_xb += m*n;
  }
}
