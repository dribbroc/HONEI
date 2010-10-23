/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "matrix.h"
#include "basics.h"
#include "blocks.h"

/* Used for temporary block storage in orderv */
/* Assume no block is bigger that 1000 words */
double temp[1000];

/**********************************************************************/

void init_matrix
(matrix *A,
 SPAI_Comm comm)
{

  A->comm = comm;

#ifdef MPI
  MPI_Comm_size(comm,&A->numprocs);
  MPI_Comm_rank(comm,&A->myid);
  MPI_Barrier(comm);
#else
  A->numprocs = 1;
  A->myid = 0;
#endif

  A->mnl = 0;
  A->my_start_index = 0;

  A->transposed = 0;
  A->n = 0;
  A->bs = 0;
  A->max_block_size = 0;
  A->maxnz = 0;
  A->mnls = NULL;
  A->start_indices = NULL;
  A->my_nnz = 0;
  A->lines = NULL;
  A->pe = NULL;
  A->block_sizes = NULL;
}

/**********************************************************************/

matrix *new_matrix
(SPAI_Comm comm)
{
  matrix *A;

  A = (matrix *) mmalloc(sizeof(matrix));
  if (! A) {
    printf("\n failed to malloc matrix.");
    exit(1);
  }

  init_matrix(A,comm);
  return(A);
}

/**********************************************************************/

void init_compressed_lines(clines *lines)
{
  lines->A = NULL;
  lines->A_buf = NULL;
  lines->ptrs = NULL;
  lines->ptrs_buf = NULL;
  lines->rptrs = NULL;
  lines->rptrs_buf = NULL;
  lines->len = NULL;
  lines->rlen = NULL;
  lines->block_sizes = NULL;
}

/**********************************************************************/

clines *new_compressed_lines(int nl, int row_structure)
{
  clines *lines;
  int i;

  lines = (clines *)
    mmalloc(sizeof(clines));
  if (! lines) {
    printf("\n failed to malloc line.");
    exit(1);
  }

  init_compressed_lines(lines);

  /* ************** column structure ***************** */

  lines->A =
    (double **) mmalloc(sizeof(double *)*nl);
  for (i=0; i<nl; i++)
    lines->A[i] = NULL;

  lines->ptrs = (int **)
    mmalloc(sizeof(int *)*nl);
  for (i=0; i<nl; i++)
    lines->ptrs[i] = NULL;
  lines->ptrs_buf = NULL;

  lines->len = (int *)
    mmalloc(sizeof(int)*nl);

  lines->slen = (int *)
    mmalloc(sizeof(int)*nl);


  /* ************** row structure ***************** */

  if (row_structure) {
    lines->rptrs = (int **)
      mmalloc(sizeof(int *)*nl);
    for (i=0; i<nl; i++)
      lines->rptrs[i] = NULL;
    lines->rptrs_buf = NULL;

    lines->rlen = (int *)
      mmalloc(sizeof(int)*nl);
  }
  else {
    lines->rptrs = NULL;
    lines->rptrs_buf = NULL;
    lines->rlen = NULL;
  }

  if ((! lines->A)           ||
      (! lines->ptrs)        ||
      (! lines->len)           ) {
    printf("\n failed to malloc line buffers\n");
    exit(1);
  }

  /* initialize */
  for (i=0; i<nl; i++) {
    lines->len[i] = 0;
    if (row_structure) lines->rlen[i] = 0;
    lines->A[i] = NULL;
    lines->ptrs[i] = NULL;
  }

  return(lines);
}

/**********************************************************************/

void free_compressed_lines(matrix *M)
{
  clines *lines;
  int i;

  if (! M->lines) return;
  lines = M->lines;

  if (lines->A_buf) free(lines->A_buf);
  else {
    for (i=0; i<M->mnl; i++) {
      if (lines->A[i]) free(lines->A[i]);
    }
  }

  if (lines->ptrs_buf) free(lines->ptrs_buf);
  else
    for (i=0; i<M->mnl; i++) {
      if (lines->ptrs[i]) free(lines->ptrs[i]);
    }

  if (lines->rptrs) {
    if (lines->rptrs_buf) free(lines->rptrs_buf);
    else
      for (i=0; i<M->mnl; i++) {
	if (lines->rptrs[i]) free(lines->rptrs[i]);
      }
  }

  if (lines->A) free(lines->A);
  if (lines->ptrs) free(lines->ptrs);
  if (lines->len) free(lines->len);
  if (lines->slen) free(lines->slen);
  if (lines->rptrs) free(lines->rptrs);
  if (lines->rlen) free(lines->rlen);

  if (lines) free(lines);
}

/**********************************************************************/

matrix *clone_matrix
(matrix *A)
{
  matrix *M;
  int i,ii;
  SPAI_Comm comm;

  comm = A->comm;

  M = new_matrix(comm);

  M->transposed = A->transposed;

  M->n = A->n;
  M->bs = A->bs;
  M->max_block_size = A->max_block_size;

  /* maxnz is currently unknown */
  M->maxnz = -1;

  M->mnls = (int *) mmalloc(sizeof(int)*M->numprocs);
  M->start_indices = (int *) mmalloc(sizeof(int)*M->numprocs);

#ifdef MPI

  MPI_Barrier(comm);
  MPI_Allgather((void *) &A->mnl, 1, MPI_INT,
		(void *) M->mnls, 1, MPI_INT,
		comm);

#else

  M->mnls[0] = A->mnls[0];

#endif

  M->start_indices[0] = 0;
  for (i=1; i<M->numprocs; i++)
    M->start_indices[i] = M->start_indices[i-1] + M->mnls[i-1];


  M->mnl = M->mnls[M->myid];
  M->my_start_index = M->start_indices[M->myid];

  /* no row structure in a cloned matrix */
  M->lines = new_compressed_lines(M->mnl,0);

  M->pe = (int *) mmalloc(sizeof(int)*A->n);
  for (i=0; i<A->n; i++) {
    M->pe[i] = A->pe[i];
  }

  M->block_sizes = (int *) mmalloc(sizeof(int)*A->n);
  for (i=0; i<A->n; i++) M->block_sizes[i] = A->block_sizes[i];

#ifdef MPI
  MPI_Barrier(comm);
#endif

  return(M);

}

/**********************************************************************/
/* The "sp_" is used because "free_matrix" conflicts with a routine
   in MATLAB. */

void sp_free_matrix(matrix *M)
{

  if (M) {

    free_compressed_lines(M);
    if (M->mnls) free(M->mnls);
    if (M->start_indices) free(M->start_indices);
    if (M->pe) free(M->pe);
    if (M->block_sizes) free(M->block_sizes);

    free(M);

  }

}

/**********************************************************************/

void order_pointers(matrix *M)
{
  int k,i,start,len,pe;
  clines *lines;
  int col;

  lines = M->lines;

  for (k=0; k<M->mnl; k++) {

    col = M->my_start_index + k;

    orderv(lines->A[k],
	   lines->ptrs[k],
	   lines->len[k],
	   M,
	   M->block_sizes[col],
	   M->max_block_size);
  }

  /* row structure ? */
  if (M->lines->rptrs) {
    for (k=0; k<M->mnl; k++) {

      col = M->my_start_index + k;

      orderv(NULL,
	     lines->rptrs[k],
	     lines->rlen[k],
	     M,
	     M->block_sizes[col],
	     M->max_block_size);

    }
  }

}

/**********************************************************************/
/* orders V in increasing order. */

void orderv
(double *v,
 int *v_ind,
 int lenv,
 matrix *A,
 int this_block_size,
 int max_block_size)
{
  int done,i;
  int itemp;
  int nbytes;
  int start_im1,start_i,bs_im1,bs_i;

  if (lenv > 1) {
    do {
      done = 1;

      start_im1 = 0;
      bs_im1 = block_size(v_ind[0],A)*this_block_size;

      for (i=1; i<lenv; i++) {

	start_i = start_im1 + bs_im1;
	bs_i = block_size(v_ind[i],A)*this_block_size;

        if (v_ind[i] < v_ind[i-1]) {
          done = 0;
          itemp = v_ind[i];
          v_ind[i] = v_ind[i-1];
          v_ind[i-1] = itemp;

	  /* Swap the two blocks */
	  if (v) {
	    memcpy((void *) temp,
		   (void *) &v[start_i],
		   bs_i*sizeof(double));
	    memcpy((void *) &v[start_i],
		   (void *) &v[start_im1],
		   bs_im1*sizeof(double));
	    memcpy((void *) &v[start_im1],
		   (void *) temp,
		   bs_i*sizeof(double));
	  }

	  start_im1 += bs_i;

	}

	else {
	  start_im1 += bs_im1;
	  bs_im1 = bs_i;
	}


      }
    }
    while (! done);
  }

}

/**********************************************************************/
/* computes A*v -> w for compressed-column storage.
   First, w is set to zero.

   This is a simple version that does a global gather and does not
   overlap communication.

*/

void A_times_v_cc
(matrix *A,
 vector *v,
 vector *w)
{
  int k,i,j,len,index;
  clines *cols;
  double *wv,*vv,*vals;
  SPAI_Comm comm;

  comm = A->comm;

  /* gather the remote parts of v */

#ifdef MPI

  vv = (double *) mmalloc(A->n*sizeof(double));
  wv = (double *) mmalloc(A->n*sizeof(double));
  MPI_Barrier(comm);
  gather_vec(A,v,vv);

#else

  vv = v->v;
  wv = w->v;

#endif

  for (j=0; j<A->n; j++) wv[j] = 0.0;
  cols = A->lines;

  for (k=0, index=A->my_start_index;
       k<A->mnl;
	 k++,
	 index++) {

    vals = cols->A[k];
    len = cols->len[k];

    for (i=0; i<len; i++) {

      j = cols->ptrs[k][i];
      wv[j] += vals[i] * vv[index];

    }
  }

#ifdef MPI

  MPI_Barrier(comm);
  MPI_Allreduce(wv,vv,A->n, MPI_DOUBLE, MPI_SUM, comm);
  for (j=0,
	 index=A->my_start_index;
       j<A->mnl;
       j++,
	 index++) {
    w->v[j] = vv[index];
  }

  free(vv);
  free(wv);

#endif

}

/**********************************************************************/
/* computes A*v -> w for compressed-row storage.
   First, w is set to zero.

   This is a simple version that does a global gather and does not
   overlap communication.

*/

void A_times_v_rc
(matrix *A,
 vector *v,
 vector *w)
{
  int k,kb,i,j,jj,len;
  clines *rows;
  double *wv,*vv,*row;
  SPAI_Comm comm;

  comm = A->comm;

  /* gather the remote parts of v */

#ifdef MPI

  vv = (double *) mmalloc(A->n*sizeof(double));
  wv = (double *) mmalloc(A->n*sizeof(double));
  MPI_Barrier(comm);
  gather_vec(A,v,vv);

#else

  vv = v->v;
  wv = w->v;

#endif

  for (j=0; j<A->n; j++) wv[j] = 0.0;
  rows = A->lines;

  for (k=0, kb=A->my_start_index;
       k<A->mnl;
       k++, kb++) {

    row = rows->A[k];
    len = rows->len[k];

    for (i=0; i<len; i++) {
      j = rows->ptrs[k][i];
      wv[kb] += row[i] * vv[j];
    }
  }

#ifdef MPI

  MPI_Barrier(comm);
  MPI_Allreduce(wv,vv,A->n, MPI_DOUBLE, MPI_SUM, comm);
  for (j=0, jj=A->my_start_index;
       j<A->mnl;
       j++, jj++) {
    w->v[j] = vv[jj];
  }

  free(vv);
  free(wv);

#endif

}

/**********************************************************************/

#ifdef MPI
void gather_vec
(matrix *A,
 vector *v,
 double *gather_buf)
{

  int i;

  MPI_Allgatherv
    ((void *) v->v, v->mnl, MPI_DOUBLE,
     (void *) gather_buf, A->mnls, A->start_indices, MPI_DOUBLE,
     A->comm);

}
#endif

/**********************************************************************/

int count_nonzeros
(matrix *A)
{
  int i,j,mynnz=0,nnz;

  for (i=0; i<A->mnl; i++) mynnz +=A->lines->len[i];

#ifdef MPI

  MPI_Allreduce(&mynnz, &nnz, 1,
	     MPI_INT, MPI_SUM, A->comm);

#else

  nnz = mynnz;

#endif

  return(nnz);
}


/**********************************************************************/
/* writes in Matrix Market format */

#define MatrixMarketBanner "%%MatrixMarket"

void write_matrix_mm
(matrix *A,
 const char *filename,
 int transposed)
{
  int pe,i,j,len,k, row,col, nnz, next;
  int ptr;
  int *pbuf;
  double *abuf;
  int bs,bs2;
  int ib,jb;
  FILE *fptr;
  char fullfilename[1024];
  char cat_cmd[1024];
  char rm_cmd[1024];
  SPAI_Comm comm;

  comm = A->comm;

  bs = A->bs;
  bs2 = bs*bs;

  if (A->numprocs > 1) {
    sprintf(fullfilename,"%s_tmp%5.5d",filename,A->myid);
    sprintf(cat_cmd,"cat %s_tmp* > %s",filename,filename);
    sprintf(rm_cmd,"rm -f %s_tmp*",filename);
  }
  else
    sprintf(fullfilename,"%s",filename);

  fptr = fopen(fullfilename,"w");

  nnz = count_nonzeros(A);

  /* write Matrix-Market header */
  if (A->myid == 0) {
    fprintf(fptr, "%s ", MatrixMarketBanner);
    fprintf(fptr,"matrix coordinate real general\n");
    fprintf(fptr,"%d %d %d\n",A->n,A->n,nnz);
    fflush(fptr);
  }

  for (j=0; j<A->mnl; j++) {

    for (i=0, next=0; i<A->lines->len[j]; i++) {

      ptr = A->lines->ptrs[j][i];

      if (! transposed) {
	row = ptr;
	col = j+A->my_start_index;
      }
      else {
	col = ptr;
	row = j+A->my_start_index;
      }

      if (A->bs == 1) {
	fprintf(fptr,"%d %d %le\n",
		row+1, col+1,
		  A->lines->A[j][i]);
      }
      else {
	fprintf(fptr,"%d %d\n",
		row+1,col+1);
	write_block(fptr,&(A->lines->A[j][next]),
		    A->block_sizes[row],
		    A->block_sizes[col]);
	next += (A->block_sizes[row]*A->block_sizes[col]);
      }
    }
  }

  fflush(fptr);
  fclose(fptr);

#ifdef MPI

  MPI_Barrier(comm);

#endif

  if (A->numprocs > 1) {
    if (A->myid == 0) {
      system(cat_cmd);
    }
  }

#ifdef MPI

  MPI_Barrier(comm);

#endif

  if (A->numprocs > 1) {
    if (A->myid == 0) {
      system(rm_cmd);
    }
  }

}

/**********************************************************************/

void matrix_statistics
(matrix *A,
 char *s)
{
  int i,k,row;
  int nnz_local,mincol_local,minrow_local=-1,maxcol_local,maxrow_local=-1;
  int nnz,mincol,minrow=-1,maxcol,maxrow=-1;
  int *lencol,*lenrow=0;
  int index;
  int max_block_size,min_block_size;
  double avg_block_size;
  int row_structure;
  SPAI_Comm comm;

  comm = A->comm;

  if (A->lines->rptrs) row_structure = 1;
  else row_structure = 0;

  lencol = (int *) mmalloc(A->n*sizeof(int));
  if (row_structure) lenrow = (int *) mmalloc(A->n*sizeof(int));

  for (k=0; k<A->n; k++) lencol[k] = -1;
  if (row_structure)
    for (k=0; k<A->n; k++) lenrow[k] = -1;
  nnz_local = 0;
  for (k=0; k<A->mnl; k++) {
    index = k + A->my_start_index;
    lencol[index] = A->lines->len[k];
    if (row_structure) lenrow[index] = A->lines->rlen[k];
    nnz_local += lencol[index];
  }

  maxcol_local = lencol[0];
  if (row_structure) maxrow_local = lenrow[0];
  mincol_local = 0;
  if (row_structure) minrow_local = 0;
  for (k=0; k<A->n; k++) {
    if (lencol[k] != -1) {
      if (mincol_local) {
	if (mincol_local > lencol[k]) mincol_local = lencol[k];
      }
      else mincol_local = lencol[k];
      if (maxcol_local < lencol[k]) maxcol_local = lencol[k];
    }
    if (row_structure) {
      if (lenrow[k] != -1) {
	if (minrow_local) {
	  if (minrow_local > lenrow[k]) minrow_local = lenrow[k];
	}
	else minrow_local = lenrow[k];
	if (maxrow_local < lenrow[k]) maxrow_local = lenrow[k];
      }
    }
  }

#ifdef MPI

  MPI_Barrier(comm);
  MPI_Reduce(&nnz_local,&nnz,1,
	     MPI_INT,MPI_SUM,0,comm);
  MPI_Reduce(&mincol_local,&mincol,1,
	     MPI_INT,MPI_MIN,0,comm);
  if (row_structure) MPI_Reduce(&minrow_local,&minrow,1,
				MPI_INT,MPI_MIN,0,comm);
  MPI_Reduce(&maxcol_local,&maxcol,1,
	     MPI_INT,MPI_MAX,0,comm);
  if (row_structure) MPI_Reduce(&maxrow_local,&maxrow,1,
				MPI_INT,MPI_MAX,0,comm);


#else
  nnz = nnz_local;
  mincol = mincol_local;
  if (row_structure) minrow = minrow_local;
  maxcol = maxcol_local;
  if (row_structure) maxrow = maxrow_local;
#endif

  max_block_size = min_block_size = avg_block_size =
    A->block_sizes[0];
  for (k=1; k<A->n; k++) {
    if (A->block_sizes[k] < min_block_size)
      min_block_size = A->block_sizes[k];
    if (A->block_sizes[k] > max_block_size)
      max_block_size = A->block_sizes[k];
    avg_block_size += A->block_sizes[k];
  }
  A->max_block_size = max_block_size;
  avg_block_size /= A->n;
  A->tot_nnz = nnz;
  if (A->myid == 0) {
    printf("\n");
    printf("Matrix statistics for %s\n",s);
    printf("--------------------------------------\n");
    printf("Size of matrix              =  %d\n",A->n);
    printf("Total number of nnz entries =  %d\n",nnz);
    printf("Max. nb of nnz per column   =  %d\n",maxcol);
    if (row_structure)
      printf("Max. nb of nnz per row      =  %d\n",maxrow);
    else
      printf("Max. nb of nnz per row      =  ***\n");
    printf("Min. nb of nnz per column   =  %d\n",mincol);
    if (row_structure)
      printf("Min. nb of nnz per row      =  %d\n",minrow);
    else
      printf("Min. nb of nnz per row      =  ***\n");
    printf("Min. block size             =  %d\n",min_block_size);
    printf("Max. block size             =  %d\n",max_block_size);
    printf("Avg. block size             =  %le\n",avg_block_size);
    if (! row_structure)
      printf("  *** => this matrix has no explicit row structure\n");
    printf("\n");
  }

  free(lencol);
  if (row_structure) free(lenrow);

}

/**********************************************************************/

int block_size(int ptr, matrix *A)
{
  return(A->block_sizes[ptr]);
}

/**********************************************************************/

int calc_maxnz
(matrix *A)
{
  int maxnz,j,this_nz,this_bs,k,bs,tmp;
  int row,col,this_row,this_col;
  SPAI_Comm comm;

  comm = A->comm;

  maxnz = 0;

  for (j=0; j<A->mnl; j++) {

    /* check the column structure */
    this_nz = 0;
    this_col = j+A->my_start_index;
    this_bs = A->block_sizes[this_col];
    for (k=0; k<A->lines->len[j]; k++) {
      row = A->lines->ptrs[j][k];
      bs = A->block_sizes[row];
      this_nz += (this_bs*bs);
    }
    if (this_nz > maxnz) maxnz = this_nz;

    /* row structure ? */
    if (A->lines->rptrs) {
      this_nz = 0;
      this_row = j+A->my_start_index;
      this_bs = A->block_sizes[this_row];
      for (k=0; k<A->lines->rlen[j]; k++) {
	col = A->lines->rptrs[j][k];
	bs = A->block_sizes[col];
	this_nz += (this_bs*bs);
      }
      if (this_nz > maxnz) maxnz = this_nz;
    }

  }

#ifdef MPI
  MPI_Barrier(comm);
  MPI_Allreduce(&maxnz, &tmp, 1,
	     MPI_INT, MPI_MAX, comm);
  maxnz = tmp;
#endif

  return(maxnz);
}

