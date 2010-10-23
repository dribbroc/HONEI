/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "variable_blocks.h"

/**********************************************************************/
/* Convert a scalar matrix A to a block matrix */

matrix *block_matrix
(matrix *A,
 int bs,
 int upper_bs_limit,
 int verbose)
{
  matrix *B;

  if (verbose) start_timer(ident_block_matrix);

  if (bs) B = constant_block_matrix(A,bs);
  else    B = variable_block_matrix(A,upper_bs_limit);

  if (verbose) {
    stop_timer(ident_block_matrix);
    report_times(ident_block_matrix,"block_matrix",0,A->comm);

  }
  return(B);

}

/**********************************************************************/
/* Convert a scalar matrix A to a constant-block matrix
 (except for possibly the last few columns). */

matrix *constant_block_matrix
(matrix *A,
 int block_size)
{
  matrix *V;
  int nblocks,*block_sizes;

  find_constant_blocks(A,block_size,&block_sizes,&nblocks);
  V = convert_to_block_matrix(A,nblocks,block_sizes);

  V->bs = block_size;

  free(block_sizes);

  return(V);

}

/**********************************************************************/

void find_constant_blocks
(matrix *A,
 int block_size,
 int **block_sizes_ptr,
 int *nblocks_ptr)
{
  int *block_sizes;
  int nblocks,i,odd_block;

  nblocks = (A->mnl)/block_size;
  odd_block = (A->mnl) % block_size;
  if (odd_block) nblocks++;

  block_sizes =
    *block_sizes_ptr = (int *) mmalloc(nblocks*sizeof(int));
  for (i=0; i<nblocks; i++) block_sizes[i] = block_size;
  if (odd_block) block_sizes[nblocks-1] = odd_block;

  *nblocks_ptr = nblocks;
  *block_sizes_ptr = block_sizes;

}

/**********************************************************************/
/* Convert a scalar matrix A to a variable-block matrix */

matrix *variable_block_matrix
(matrix *A,
 int upper_bs_limit)
{
  matrix *V;
  int nblocks,*block_sizes;

  find_diagonal_blocks(A,upper_bs_limit,&block_sizes,&nblocks);
  V = convert_to_block_matrix(A,nblocks,block_sizes);
  V->bs = 0;

  free(block_sizes);

  return(V);

}

/**********************************************************************/

void find_diagonal_blocks
(matrix *A,
 int upper_bs_limit,
 int **block_sizes_ptr,
 int *nblocks_ptr)
{
  int *block_numbers,*block_sizes;
  int nblocks,i,iend,ii;

  block_numbers = (int *) mmalloc(A->mnl*sizeof(int));

  nblocks = 0;
  i = 0;

  do {
    iend = find_diagonal_block(A,i,upper_bs_limit);
    for (ii=i; ii<=iend; ii++)
      block_numbers[ii] = nblocks;
    nblocks++;
    i = iend+1;
  }
  while (i < A->mnl);

  block_sizes =
    *block_sizes_ptr = (int *) mmalloc(nblocks*sizeof(int));
  for (i=0; i<nblocks; i++) block_sizes[i] = 0;
  for (i=0; i<A->mnl; i++) {
    block_sizes[block_numbers[i]]++;
  }

  free(block_numbers);

  *nblocks_ptr = nblocks;

}

/**********************************************************************/

int find_diagonal_block(matrix *A, int i, int upper_bs_limit)
{
  int k,bs,max_bs;
  int initial_run_length(matrix *, int);
  int check_next_run(matrix *, int, int, int);

  bs = max_bs = initial_run_length(A,i);
  if (max_bs == 1)
    return(i);

  if (upper_bs_limit)
    if (max_bs > upper_bs_limit)
      bs = max_bs = upper_bs_limit;

  for (k=0; k<max_bs; k++) {

    if ((k > 0) && (i+k) == A->mnl) return(i+k-1);

    bs = check_next_run(A,i,i+k,bs);
    if (! bs) {
      return(i+k-1);
    }
  }

  return(i+max_bs-1);

}

/**********************************************************************/

int initial_run_length(matrix *A, int i)
{
  int k,index,bs;
  int ii;
  int ptr;

  /* Find the index of the diagonal entry, if any */
  for (k=0, index=-1; k<A->lines->len[i]; k++) {
    ptr = A->lines->ptrs[i][k];
    if ((A->pe[ptr] == A->myid) &&
	((ptr - A->my_start_index) == i)) {
      index = k;
      break;
    }
  }

  /* no diagonal entry */
  if (index < 0) return(1);

  bs = 1;

  /* Count consecutive entries */
  for (k=index+1,
	 ii=i+1;
       k<A->lines->len[i];
       k++) {
    ptr = A->lines->ptrs[i][k];
    if ((A->pe[ptr] == A->myid) &&
	((ptr - A->my_start_index) == ii)) {
      bs++;
      ii++;
    }
    else return(bs);
  }

  return(bs);

}


/**********************************************************************/

int check_next_run(matrix *A, int i, int i_next, int bs)
{
  int k,index,new_bs,ii;
  int ptr;

  /* Look for the index of entry (myid, i) */
  for (k=0, index=-1; k<A->lines->len[i_next]; k++) {
    ptr = A->lines->ptrs[i_next][k];
    if ((A->pe[ptr] == A->myid) &&
	((ptr - A->my_start_index) == i)) {
      index = k;
      break;
    }
  }

  /* return 0 if it doesn't exist */
  if (index < 0) {
    return(0);
  }

  /* Count consectutive entries */
  for (k=index+1,
	 ii=i+1,
	 new_bs=1;
       k<A->lines->len[i_next];
       k++) {
    ptr = A->lines->ptrs[i_next][k];
    if ((A->pe[ptr] == A->myid) &&
	((ptr - A->my_start_index) == ii)) {
      new_bs++;
      ii++;
    }
    else break;
  }

  /* Does the run include the diagonal entry of i_next? */
  if (new_bs > i_next-i) {

    /* Yes, it does */
    if (new_bs > bs) {
      return(bs);
    }
    else {
      return(new_bs);
    }
  }

  return(0);

}

/**********************************************************************/

matrix *convert_to_block_matrix
(matrix *A,
 int nblocks_local,
 int *block_sizes_local)
{
  matrix *B;
  double *fullcol;
  int *block_bitvec,*block_start,*block_index_map;
  int pe,row,col,i,j,ib,jb,starti,startj,jcol,k;
  int count,found_block,n,m,total_height;
  int max_n,max_mnl;

  int start,next;
  int start_indx;

  /* I'm assuming that no matrix element will be equal to 1.0E300 */
  double infinity = 1.0E300;

  int col_buf_size,row_buf_size,A_buf_size;
  int max_A_buf_size;
  int max_col_buf_size;
  int max_row_buf_size=1;
  int index,num;
  int *ptr_adr,*rptr_adr;
  double *A_adr;
  int next_ptr,next_rptr,next_Aptr;

  B = new_matrix(A->comm);

  B->transposed = A->transposed;

  B->mnls = (int *) mmalloc(sizeof(int)*B->numprocs);

#ifdef MPI

  MPI_Barrier(B->comm);
  MPI_Allgather((void *) &nblocks_local, 1, MPI_INT,
		(void *) B->mnls, 1, MPI_INT,
		B->comm);

#else

  B->mnls[0] = nblocks_local;

#endif

  max_mnl = 0;
  for (pe=0; pe<B->numprocs; pe++) {
    if (max_mnl < B->mnls[pe]) max_mnl = B->mnls[pe];
  }

  B->n = 0;
  for (i=0; i<B->numprocs; i++) B->n += B->mnls[i];

  B->start_indices = (int *) mmalloc(sizeof(int)*B->numprocs);
  B->start_indices[0] = 0;
  B->pe = (int *) mmalloc(sizeof(int)*B->n);
  for (pe=1; pe<B->numprocs; pe++)
    B->start_indices[pe] = B->start_indices[pe-1] + B->mnls[pe-1];


  B->mnl = B->mnls[B->myid];
  B->my_start_index = B->start_indices[B->myid];

  for (pe=0; pe<B->numprocs; pe++) {
    start_indx = B->start_indices[pe];
    for (i=0; i<B->mnls[pe]; i++)
      B->pe[start_indx+i] = pe;
  }

  B->block_sizes = (int *) mmalloc(B->n*sizeof(int));

#ifdef MPI

  MPI_Barrier(B->comm);
  MPI_Allgatherv
    ((void *) block_sizes_local, B->mnl, MPI_INT,
     (void *) B->block_sizes, B->mnls, B->start_indices, MPI_INT,
     B->comm);

#else

  for (i=0; i<B->n; i++) B->block_sizes[i] = block_sizes_local[i];

#endif

  B->bs = 0;

  max_n =0;
  for(ib=0; ib<B->n; ib++) {
    if (max_n < B->block_sizes[ib]) max_n = B->block_sizes[ib];
  }

  B->max_block_size = max_n;

  /* row structure ? */
  if (A->lines->rptrs)
    B->lines = new_compressed_lines(B->mnl,1);
  else
    B->lines = new_compressed_lines(B->mnl,0);


  fullcol = (double *) mmalloc(A->n*max_n*sizeof(double));
  block_bitvec = (int *) mmalloc(B->n*sizeof(int));
  for (k=0; k<A->n*max_n; k++) fullcol[k] = infinity;
  for (ib=0; ib<B->n; ib++) block_bitvec[ib] = 0;

  /* This maps scalar inidices to block indices. */
  block_index_map = (int *) mmalloc(A->n*sizeof(int));
  i = 0;
  for (ib=0;
       ib<B->n;
       ib++) {
    for (k=0; k<B->block_sizes[ib]; k++) {
      block_index_map[i+k] = ib;
    }
    i += B->block_sizes[ib];
  }

  /* This gives the scalar staring index for each block. */
  block_start = (int *) mmalloc(B->n*sizeof(int));
  block_start[0] = 0;
  for (ib=1; ib<B->n; ib++)
    block_start[ib] = block_start[ib-1] + B->block_sizes[ib-1];


  /* Install the columns */
  /* i and j are row and column indices into the scalar matrix.
     ib and jb are row and column indices into the block matrix */

  for (startj=0,
	 jb=0;
       jb<max_mnl;
       startj += block_sizes_local[jb],
	 jb++) {

    if (jb < B->mnl) {

      n = block_sizes_local[jb];
      count = 0; /* number of blocks in this column */
      total_height = 0;

      /* fill fullcol */
      for (j=startj; j<startj+n; j++) {
	for (k=0; k<A->lines->len[j]; k++) {

	  i = A->lines->ptrs[j][k];
	  jcol = j - startj;
	  fullcol[jcol*A->n + i] = A->lines->A[j][k];

	  ib = block_index_map[i];
	  if (! block_bitvec[ib]) {
	    block_bitvec[ib] = 1;
	    count++;
	    total_height += B->block_sizes[ib];
	  }

	}
      }

      /* Install block-column jb */
      B->lines->ptrs[jb] =
	(int *) mmalloc(count*sizeof(int));
      B->lines->A[jb] =
	(double *) mmalloc(total_height*n*sizeof(double));

      B->lines->len[jb] = count;
      B->lines->slen[jb] = total_height;
      for (ib=0, k=0, next=0;
	   ib<B->n;
	 ib++) {
	if (block_bitvec[ib]) {
	  B->lines->ptrs[jb][k] = ib;
	  start = block_start[ib];
	  m = B->block_sizes[ib];
	  for (col=0; col<n; col++) {
	    for (row=0; row<m; row++) {
	      if (fullcol[start + col*A->n + row] == infinity) {
		fullcol[start + col*A->n + row] = 0.0;
	      }
	      B->lines->A[jb][next + col*m + row] =
		fullcol[start + col*A->n + row];
	      fullcol[start + col*A->n + row] = infinity;
	    }
	  }

	  block_bitvec[ib] = 0;

	  k++;
	  next += m*n;
	}
      }
    }
  }

  /* i and j are column and row indices into the scalar matrix.
     ib and jb are column and row indices into the block matrix */
  /* row structure ? */
  if (B->lines->rptrs) {
    for (startj=0,
	   jb=0;
	 jb<max_mnl;
	 startj += block_sizes_local[jb],
	   jb++) {

      if (jb < B->mnl) {

	n = block_sizes_local[jb];
	count = 0; /* number of blocks in this row */

	/* fill fullcol (actually fullrow) */
	for (j=startj; j<startj+n; j++) {
	  for (k=0; k<A->lines->rlen[j]; k++) {

	    i = A->lines->rptrs[j][k];
	    jcol = j - startj;
	    fullcol[jcol*A->n + i] = 1.0;

	    ib = block_index_map[i];
	    if (! block_bitvec[ib]) {
	      block_bitvec[ib] = 1;
	      count++;
	      block_start[ib] = i;
	    }

	  }
	}

	/* Install block-row structure jb */
	B->lines->rptrs[jb] =
	  (int *) mmalloc(count*sizeof(int));
	B->lines->rlen[jb] = count;
	for (ib=0, k=0, next=0;
	     ib<B->n;
	     ib++) {
	  if (block_bitvec[ib]) {
	    B->lines->rptrs[jb][k] = ib;
	    k++;
	    next += m*n;
	    block_bitvec[ib] = 0;
	  }
	}
      }

    }
  }

  /* Convert the row and column data to consistent shmalloc buffers */
  col_buf_size = row_buf_size = A_buf_size = 0;
  for (j=0, index=B->my_start_index;
       j<B->mnl;
       j++, index++) {
    col_buf_size += B->lines->len[j];
    /* row structure ? */
    if (B->lines->rptrs)
      row_buf_size += B->lines->rlen[j];
    A_buf_size += (B->block_sizes[index]*B->lines->slen[j]);
  }

#ifdef MPI
  MPI_Barrier(B->comm);
  MPI_Allreduce(&col_buf_size,&max_col_buf_size,1,
		MPI_INT,MPI_MAX,B->comm);
  /* row structure ? */
  if (B->lines->rptrs)
    MPI_Allreduce(&row_buf_size,&max_row_buf_size,1,
		  MPI_INT,MPI_MAX,B->comm);
  MPI_Allreduce(&A_buf_size,&max_A_buf_size,1,
		MPI_INT,MPI_MAX,B->comm);
#else
  max_col_buf_size = col_buf_size;
  if (A->lines->rptrs)
    max_row_buf_size = row_buf_size;
  max_A_buf_size = A_buf_size;
#endif

#ifdef SHMEM
  B->lines->ptrs_buf  = (int *)    shmalloc(max_col_buf_size*sizeof(int));
  /* row structure ? */
  if (B->lines->rptrs)
    B->lines->rptrs_buf = (int *)    shmalloc(max_row_buf_size*sizeof(int));
  B->lines->A_buf     = (double *) shmalloc(max_A_buf_size*sizeof(double));
#else
  B->lines->ptrs_buf  = (int *)    mmalloc(max_col_buf_size*sizeof(int));
  /* row structure ? */
  if (B->lines->rptrs)
    B->lines->rptrs_buf = (int *)    mmalloc(max_row_buf_size*sizeof(int));
  B->lines->A_buf     = (double *) mmalloc(max_A_buf_size*sizeof(double));
#endif

  for (j=0,
	 index = B->my_start_index,
	 next_ptr = next_rptr = next_Aptr = 0,
	 ptr_adr = B->lines->ptrs_buf,
	 rptr_adr = B->lines->rptrs_buf,
	 A_adr = B->lines->A_buf;
       j<B->mnl;
       j++,
	 index++) {

    num = B->lines->len[j]*sizeof(int);
    memcpy(&(B->lines->ptrs_buf[next_ptr]),B->lines->ptrs[j],num);
    next_ptr += B->lines->len[j];
    free(B->lines->ptrs[j]);
    B->lines->ptrs[j] = ptr_adr;
    ptr_adr += B->lines->len[j];

    /* row structure ? */
    if (B->lines->rptrs) {
      num = B->lines->rlen[j]*sizeof(int);
      memcpy(&(B->lines->rptrs_buf[next_rptr]),B->lines->rptrs[j],num);
      next_rptr += B->lines->rlen[j];
      free(B->lines->rptrs[j]);
      B->lines->rptrs[j] = rptr_adr;
      rptr_adr += B->lines->rlen[j];
    }

    num = B->block_sizes[index]*B->lines->slen[j]*sizeof(double);
    memcpy(&(B->lines->A_buf[next_Aptr]),B->lines->A[j],num);
    next_Aptr += B->block_sizes[index]*B->lines->slen[j];

    free(B->lines->A[j]);
    B->lines->A[j] = A_adr;
    A_adr += B->block_sizes[index]*B->lines->slen[j];
  }

  B->maxnz = calc_maxnz(B);

  free(fullcol);
  free(block_bitvec);
  free(block_start);
  free(block_index_map);

  return(B);

}

/**********************************************************************/
/* Convert a variable-block matrix A to a scalar matrix */

matrix *scalar_matrix
(matrix *B,
 int verbose)
{
  matrix *S;
  int pe,n,m,mnl,i,j,ib,jb,ib_global,jb_global;
  int jstart,istart,i_global;
  int start_indx,height,kb,col,row,len,next;
  double val;
  int *block_start;
  double *Aval;

  if (verbose) start_timer(ident_scalar_matrix);

  S = new_matrix(B->comm);

  mnl = 0;
  for (j=0; j<B->mnl; j++) {
    mnl += B->block_sizes[j+B->my_start_index];
  }

  S->mnls = (int *) mmalloc(sizeof(int)*S->numprocs);

#ifdef MPI

  MPI_Barrier(S->comm);
  MPI_Allgather((void *) &mnl, 1, MPI_INT,
		(void *) S->mnls, 1, MPI_INT,
		S->comm);

#else

  S->mnls[0] = mnl;

#endif

  S->n = 0;
  for (i=0; i<S->numprocs; i++) S->n += S->mnls[i];

  S->start_indices = (int *) mmalloc(sizeof(int)*S->numprocs);
  S->start_indices[0] = 0;
  for (pe=1; pe<S->numprocs; pe++)
    S->start_indices[pe] = S->start_indices[pe-1] + S->mnls[pe-1];

  S->mnl = S->mnls[S->myid];
  S->my_start_index = S->start_indices[S->myid];

  S->pe = (int *) mmalloc(sizeof(int)*S->n);
  for (pe=0; pe<S->numprocs; pe++) {
    start_indx = S->start_indices[pe];
    for (i=0; i<S->mnls[pe]; i++)
      S->pe[start_indx+i] = pe;
  }

  S->block_sizes = (int *) mmalloc(S->n*sizeof(int));
  for (j=0; j<S->n; j++) S->block_sizes[j]=1;

  S->transposed = B->transposed;
  S->bs = 1;
  S->max_block_size = 1;

  /* No row structure */
  S->lines = new_compressed_lines(S->mnl,0);

  /* This gives the scalar staring index for each block. */
  block_start = (int *) mmalloc(B->n*sizeof(int));
  block_start[0] = 0;
  for (ib=1; ib<B->n; ib++)
    block_start[ib] = block_start[ib-1] + B->block_sizes[ib-1];

  /* i and j are local row and column indices into the scalar matrix.
     ib and jb are local row and column indices into the block matrix */
  for (jb=0,
	 jb_global=B->my_start_index,
	 jstart=0;
       jb<B->mnl;
       jb++,
	 jstart+=B->block_sizes[jb_global],
	 jb_global++) {

    n = B->block_sizes[jb_global];

    /* Determine scalar height of block column jb, */
    height = 0;
    for (kb=0; kb<B->lines->len[jb]; kb++) {
      ib_global = B->lines->ptrs[jb][kb];
      height += B->block_sizes[ib_global];
    }

    /* Allocate scalar columns in S associated with
       block column jb */
    for (j=jstart;
	 j<jstart+B->block_sizes[jb_global];
	 j++) {
      S->lines->ptrs[j] = (int *) mmalloc(height*sizeof(int));
      S->lines->A[j] = (double *) mmalloc(height*sizeof(double));
      S->lines->len[j] = 0;
    }

    /* Install scalar values associated with block column jb */
    for (kb=0,next=0; kb<B->lines->len[jb]; kb++) {
      Aval = &(B->lines->A[jb][next]);
      ib_global = B->lines->ptrs[jb][kb];
      istart = block_start[ib_global];
      m = B->block_sizes[ib_global];
      for (col=0; col<n; col++) {
	for (row=0; row<m; row++) {
	  val = Aval[col*m+row];
	  if (val != 0.0) {
	    j = jstart+col;
	    i_global = istart+row;
	    len = S->lines->len[j]++;
	    S->lines->ptrs[j][len] = i_global;
	    S->lines->A[j][len] = val;
	  }
	}
      }
      next += m*n;
    }
  }

  S->maxnz = calc_maxnz(S);

  free(block_start);

  if (verbose) {
    stop_timer(ident_scalar_matrix);
    report_times(ident_scalar_matrix,"scalar_matrix",0,S->comm);
  }

  return(S);
}
