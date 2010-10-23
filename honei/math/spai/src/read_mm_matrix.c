/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "read_mm_matrix.h"

/**********************************************************************/
/* Read a general Matrix Market file into a block representation. */

/* NOTE: For historical reasons, this routine can (in principle) read
   a matrix in "constant block" form; i.e., bs != 1. I recommend NOT
   doing this. Read all matrices as having scalar entries (bs == 1) and
   use the routines "constant_block_matrix" or "variable_block_matrix"
   to convert to block form.
*/


/* (transpose == 0) =>The matrix values are stored by column. */
/* (transpose == 1) =>The matrix values are stored by row.    */

matrix *read_mm_matrix
(char *matrix_file,
 int bs,
 int chunk_size,
 int symmetric_pattern,
 int transpose,
 int binary,
 int verbose,
 SPAI_Comm comm)
{
  FILE *f;
  int M,N,nnz,nrows,ncols;
  int nl,mnl,split_pe,split_indx,start_indx;
  mm_data *rows,*cols;
  int nnz_rows,nnz_cols;
  matrix *A;
  int actual_M,actual_N,actual_nnz;

  int numprocs,myid;
#ifdef MPI
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myid);
  MPI_Barrier(comm);
#else
  numprocs = 1;
  myid = 0;
#endif

  if (verbose) start_timer(ident_read_mm_matrix);

  f = fopen(matrix_file,"r");
  /* If the file doesn't exist (on the Origin2000) this generates a fatal error.
     I don't know why. */
  if (! f)
    if (myid == 0) {
      printf("Couldn't find matrix file: %s\n",matrix_file);
      exit(1);
    }

  if (binary)
    read_mm_matrix_header_binary(f, &actual_M,&actual_N,&actual_nnz);
  else
    read_mm_matrix_header_ascii(f, &actual_M,&actual_N,&actual_nnz);

  if (actual_M != actual_N) {
    if (myid == 0) {
      printf("Matrix must be square.\n");
      printf("  M=%d N=%d\n",actual_M,actual_N);
    }
    exit(1);
  }

  if (actual_N % bs) {
    if (myid == 0) {
      printf("WARNING: Number of cols not multiple of block size.\n");
      printf("         filling in with identity elements\n");
    }
    M = actual_M/bs + 1;
    M = N = M*bs;
    nnz = actual_nnz + (M-actual_M);
  }
  else {
    M = actual_M;
    N = actual_N;
    nnz = actual_nnz;
  }

  /* Determine the distribution of rows and columns across
     processors. */
  basic_distribution
    (comm,
     N/bs, chunk_size, &nl,&mnl,&split_pe,&split_indx,&start_indx);

  if (binary)
    read_mm_data_binary
      (f,
       transpose,
       N,nnz,
       actual_N,actual_nnz,
       mnl*bs,start_indx*bs,
       &rows,&nnz_rows,
       &cols,&nnz_cols);

  else
    read_mm_data_ascii
      (f,
       transpose,
       N,nnz,
       actual_N,actual_nnz,
       mnl*bs,start_indx*bs,
       &rows,&nnz_rows,
       &cols,&nnz_cols);

  fclose(f);

  A = mm_to_matrix
    (transpose,
     symmetric_pattern,
     N,
     bs,
     mnl,
     rows,nnz_rows,
     cols,nnz_cols,
     comm);

  free(rows);
  free(cols);

  order_pointers(A);

  if (verbose) {
    stop_timer(ident_read_mm_matrix);
    report_times(ident_read_mm_matrix,"read_mm_matrix",0,comm);
  }

  return(A);
}

/**********************************************************************/

void skip_mm_matrix_header_ascii(FILE *f)
{
  int M,N,nnz;

  read_mm_matrix_header_ascii(f,&M,&N,&nnz);
}

/**********************************************************************/

void skip_mm_matrix_header_binary(FILE *f)
{
  int M,N,nnz;

  fread(&M, sizeof(int), (size_t) 1, f);
  fread(&N, sizeof(int), (size_t) 1, f);
  fread(&nnz, sizeof(int), (size_t) 1, f);

}

/**********************************************************************/

void read_mm_matrix_header_ascii
(FILE *f,
 int *M_ptr,
 int *N_ptr,
 int *nnz_ptr)
{
    MM_typecode matcode;
    int M, N, nnz;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (!
	(mm_is_real(matcode)) &&
	(mm_is_matrix(matcode)) &&
	(mm_is_sparse(matcode)) &&
	(mm_is_general(matcode))) {
      printf("MM matrix must be real, sparse, and general\n");
      exit(1);
    }

    if (mm_read_mtx_crd_size(f, M_ptr, N_ptr, nnz_ptr)) {
      printf("error in  mm_read_mtx_crd_size\n");
      exit(1);
    }
}


/**********************************************************************/

void read_mm_matrix_header_binary
(FILE *f,
 int *M_ptr,
 int *N_ptr,
 int *nnz_ptr)
{

  fread(M_ptr, sizeof(int), (size_t) 1, f);
  fread(N_ptr, sizeof(int), (size_t) 1, f);
  fread(nnz_ptr, sizeof(int), (size_t) 1, f);

}

/**********************************************************************/

void write_mm_data
(FILE *fptr,
 mm_data *rows,
 int nnz_rows,
 mm_data *cols,
 int nnz_cols)
{
  int k;

  fprintf(fptr,"....... rows ......\n");
  for (k=0; k<nnz_rows; k++) {
    fprintf(fptr,"%4d: %4d %4d %4d %4d %lf\n",
	    k,rows[k].i,rows[k].j,rows[k].ib,rows[k].jb,rows[k].val);
  }

  fprintf(fptr,"....... cols ......\n");
  for (k=0; k<nnz_cols; k++) {
    fprintf(fptr,"%4d: %4d %4d %4d %4d %lf\n",
	    k,cols[k].i,cols[k].j,cols[k].ib,cols[k].jb,cols[k].val);
  }
}


/**********************************************************************/
/*
   Returns matrix data in distributed form.
   Both the matrix and the transpose are available.
   For the transpose just switch "rows" and "cols".

   In the external Matrix Market format indices are 1-based.
   All internal representations are 0-based.
*/

void read_mm_data_ascii
(FILE *f,
 int transpose,
 int N,
 int nnz,
 int actual_N,
 int actual_nnz,
 int mnl,
 int start_indx,
 mm_data **rows_ptr,
 int *nnz_rows_ptr,
 mm_data **cols_ptr,
 int *nnz_cols_ptr)
{
  int my_nnz_row,my_nnz_col,row,col,lenrow,lencol;
  int i,ii,k;
  matrix *A;
  int len,rlen;
  int pos,pos1,pos2;
  int rpos2;
  double val;
  mm_data *rows,*cols;
  char line[128];

  /* count row and column nonzeros on this processor */
  my_nnz_row = my_nnz_col = 0;
  for (i=0; i<actual_nnz; i++) {
    fgets(line,128,f);

    /* Change ',' to ' ' */
    for (ii=0; line[ii]; ii++) if (line[ii] == ',') line[ii] = ' ';

    if (transpose == 0)
      sscanf(line, "%d %d %le\n", &row, &col, &val);
    else
      sscanf(line, "%d %d %le\n", &col, &row, &val);

    row--;
    col--;
    if ((row >= start_indx) && (row < (start_indx+mnl)))
      my_nnz_row++;
    if ((col >= start_indx) && (col < (start_indx+mnl)))
      my_nnz_col++;
  }
  for (i=actual_N; i<N; i++) {
    row = i;
    col = i;
    if ((row >= start_indx) && (row < (start_indx+mnl)))
      my_nnz_row++;
    if ((col >= start_indx) && (col < (start_indx+mnl)))
      my_nnz_col++;
  }

  *nnz_rows_ptr = my_nnz_row;
  *nnz_cols_ptr = my_nnz_col;


  rewind(f);
  skip_mm_matrix_header_ascii(f);

  rows = (mm_data *) mmalloc(my_nnz_row*sizeof(mm_data));

  if (! rows) {
    printf("failed to malloc rows\n");
    fflush(stdout);
    exit(1);
  }

  cols = (mm_data *) mmalloc(my_nnz_col*sizeof(mm_data));

  if (! cols) {
    printf("failed to malloc cols\n");
    fflush(stdout);
    exit(1);
  }
  *rows_ptr = rows;
  *cols_ptr = cols;

  for (lenrow=0, lencol=0, i=0;
       i<actual_nnz;
       i++) {

    fgets(line,128,f);

    /* Change ',' to ' ' */
    for (ii=0; line[ii]; ii++) if (line[ii] == ',') line[ii] = ' ';

    if (transpose == 0)
      sscanf(line, "%d %d %le\n", &row, &col, &val);
    else
      sscanf(line, "%d %d %le\n", &col, &row, &val);

    row--;
    col--;

    if ((row >= start_indx) && (row < start_indx+mnl)) {
      rows[lenrow].i = row;
      rows[lenrow].j = col;
      rows[lenrow].val = val;
      lenrow++;
    }
    if ((col >= start_indx) && (col < start_indx+mnl)) {
      cols[lencol].i = row;
      cols[lencol].j = col;
      cols[lencol].val = val;
      lencol++;
    }
  }

  for (i=actual_N; i<N; i++) {
    printf("/// inserting identity element at %d,%d\n",i,i);
    row = col = i;
    val = 1.0;

    if ((row >= start_indx) && (row < start_indx+mnl)) {
      rows[lenrow].i = row;
      rows[lenrow].j = col;
      rows[lenrow].val = val;
      lenrow++;
    }
    if ((col >= start_indx) && (col < start_indx+mnl)) {
      cols[lencol].i = row;
      cols[lencol].j = col;
      cols[lencol].val = val;
      lencol++;
    }
  }
}

/**********************************************************************/
/*
   Returns matrix data in distributed form.
   Both the matrix and the transpose are available.
   For the transpose just switch "rows" and "cols".

   In the external Matrix Market format indices are 1-based.
   All internal representations are 0-based.
*/

void read_mm_data_binary
(FILE *f,
 int transpose,
 int N,
 int nnz,
 int actual_N,
 int actual_nnz,
 int mnl,
 int start_indx,
 mm_data **rows_ptr,
 int *nnz_rows_ptr,
 mm_data **cols_ptr,
 int *nnz_cols_ptr)
{
  int my_nnz_row,my_nnz_col,row,col,lenrow,lencol;
  int i,ii,k;
  matrix *A;
  int len,rlen;
  int pos,pos1,pos2;
  int rpos2;
  double val;
  mm_data *rows,*cols;
  int rowcol[2];

  /* count row and column nonzeros on this processor */
  my_nnz_row = my_nnz_col = 0;
  for (i=0; i<actual_nnz; i++) {
    fread(rowcol, sizeof(int), (size_t) 2, f);
    fread(&val, sizeof(double), (size_t) 1, f);

    if (transpose == 0) {
      row = rowcol[0];
      col = rowcol[1];
    }
    else {
      row = rowcol[1];
      col = rowcol[0];
    }

    row--;
    col--;
    if ((row >= start_indx) && (row < (start_indx+mnl)))
      my_nnz_row++;
    if ((col >= start_indx) && (col < (start_indx+mnl)))
      my_nnz_col++;
  }
  for (i=actual_N; i<N; i++) {
    row = i;
    col = i;
    if ((row >= start_indx) && (row < (start_indx+mnl)))
      my_nnz_row++;
    if ((col >= start_indx) && (col < (start_indx+mnl)))
      my_nnz_col++;
  }

  *nnz_rows_ptr = my_nnz_row;
  *nnz_cols_ptr = my_nnz_col;

  rewind(f);
  skip_mm_matrix_header_binary(f);

  rows = (mm_data *) mmalloc(my_nnz_row*sizeof(mm_data));

  if (! rows) {
    printf("failed to malloc rows\n");
    fflush(stdout);
    exit(1);
  }

  cols = (mm_data *) mmalloc(my_nnz_col*sizeof(mm_data));

  if (! cols) {
    printf("failed to malloc cols\n");
    fflush(stdout);
    exit(1);
  }
  *rows_ptr = rows;
  *cols_ptr = cols;

  for (lenrow=0, lencol=0, i=0;
       i<actual_nnz;
       i++) {

    fread(rowcol, sizeof(int), (size_t) 2, f);
    fread(&val, sizeof(double), (size_t) 1, f);

    if (transpose == 0) {
      row = rowcol[0];
      col = rowcol[1];
    }
    else {
      row = rowcol[1];
      col = rowcol[0];
    }

    row--;
    col--;

    if ((row >= start_indx) && (row < start_indx+mnl)) {
      rows[lenrow].i = row;
      rows[lenrow].j = col;
      rows[lenrow].val = val;
      lenrow++;
    }
    if ((col >= start_indx) && (col < start_indx+mnl)) {
      cols[lencol].i = row;
      cols[lencol].j = col;
      cols[lencol].val = val;
      lencol++;
    }
  }

  for (i=actual_N; i<N; i++) {
    printf("/// inserting identity element at %d,%d\n",i,i);
    row = col = i;
    val = 1.0;

    if ((row >= start_indx) && (row < start_indx+mnl)) {
      rows[lenrow].i = row;
      rows[lenrow].j = col;
      rows[lenrow].val = val;
      lenrow++;
    }
    if ((col >= start_indx) && (col < start_indx+mnl)) {
      cols[lencol].i = row;
      cols[lencol].j = col;
      cols[lencol].val = val;
      lencol++;
    }
  }
}

/**********************************************************************/
/* (transpose == 0) =>The matrix values are stored by column. */
/* (transpose == 1) =>The matrix values are stored by row.    */

matrix *mm_to_matrix
(int transpose,
 int symmetric_pattern,
 int N,
 int bs,
 int mnl,
 mm_data *rows,
 int nnz_rows,
 mm_data *cols,
 int nnz_cols,
 SPAI_Comm comm)
{
  matrix *A;
  int i,k,start_indx,pe,local_indx,ib,jb;
  int local_row_indx,col_indx;
  int local_col_indx,row_indx;
  int current_block,current_col_index,current_row_index;
  int row,col,local_row,local_col;
  int *nnz_per_col=NULL,*nnz_per_row=NULL;
  int *mapping=NULL;
  int pos,rpos,pos1,pos2,len,rlen;
  clines *lines;
  int current_col,current_row;
  int col_buf_size,row_buf_size,A_buf_size;
  int max_col_buf_size,max_A_buf_size;
  int max_row_buf_size=1;
  int index;

  A = new_matrix(comm);

  if (! transpose) A->transposed = 0;
  else             A->transposed = 1;

  A->n = N;
  A->bs = bs;
  A->max_block_size = bs;

  A->mnls = (int *) mmalloc(sizeof(int)*A->numprocs);

  A->start_indices = (int *) mmalloc(sizeof(int)*A->numprocs);

  A->pe = (int *) mmalloc(sizeof(int)*N);

  A->block_sizes = (int *) mmalloc(sizeof(int)*N);
  for (k=0; k<N; k++) A->block_sizes[k] = bs;

#ifdef MPI

  MPI_Barrier(comm);
  MPI_Allgather((void *) &mnl, 1, MPI_INT,
		(void *) A->mnls, 1, MPI_INT,
		comm);

#else

  A->mnls[0] = mnl;

#endif

  A->start_indices[0] = 0;
  for (pe=1; pe<A->numprocs; pe++)
    A->start_indices[pe] = A->start_indices[pe-1] + A->mnls[pe-1];

  for (pe=0; pe<A->numprocs; pe++) {
    start_indx = A->start_indices[pe];
    for (i=0; i<A->mnls[pe]; i++)
      A->pe[start_indx+i] = pe;
  }

  A->mnl = A->mnls[A->myid];
  A->my_start_index = A->start_indices[A->myid];

  nnz_per_col = (int *) mmalloc(N*sizeof(int));
  if (! nnz_per_col) {
    printf("failed to malloc nnz_per_col\n");
    fflush(stdout);
    exit(1);
  }

  nnz_per_row = (int *) mmalloc(N*sizeof(int));
  if (! nnz_per_row) {
    printf("failed to malloc nnz_per_row\n");
    fflush(stdout);
    exit(1);
  }

  for (i=0; i<N; i++) nnz_per_col[i] = 0;
  for (i=0; i<N; i++) nnz_per_row[i] = 0;

  start_indx = A->my_start_index;

  /* set block indices */
  for (k=0; k<nnz_rows; k++) {
    rows[k].ib = rows[k].i/bs;
    rows[k].jb = rows[k].j/bs;
  }
  for (k=0; k<nnz_cols; k++) {
    cols[k].ib = cols[k].i/bs;
    cols[k].jb = cols[k].j/bs;
  }

  /* Order by block indices */
  qsort((void *) rows, nnz_rows, sizeof(mm_data), &rank_row_major);
  qsort((void *) cols, nnz_cols, sizeof(mm_data), &rank_column_major);

  /* Count the number of nonzero blocks in each row */
  current_block = -1;
  current_row = rows[0].ib;
  for (k=0; k<nnz_rows; k++) {
    if (rows[k].ib != current_row) {
      current_row = rows[k].ib;
      current_block = -1;
    }
    current_col_index = rows[k].jb;
    if (current_col_index != current_block) {
      current_block = current_col_index;
      local_row_indx = rows[k].ib - start_indx;
      nnz_per_row[local_row_indx]++;
    }
  }

  /* Count the number of nonzeros blocks in each col */
  current_block = -1;
  current_col = cols[0].jb;
  for (k=0; k<nnz_cols; k++) {
    if (cols[k].jb != current_col) {
      current_col = cols[k].jb;
      current_block = -1;
    }
    current_row_index = cols[k].ib;
    if (current_row_index != current_block) {
      current_block = current_row_index;
      local_col_indx = cols[k].jb - start_indx;
      nnz_per_col[local_col_indx]++;
    }
  }

  /* Include a row structure ? */
  A->my_nnz = nnz_cols;
  if (symmetric_pattern)
    A->lines = new_compressed_lines(A->mnl,0);
  else
    A->lines = new_compressed_lines(A->mnl,1);

  /**************************************************/
  /* ************* Allocate space ****************** */
  /**************************************************/

  lines = A->lines;

  /* Determine lengths */
  for (i=0; i<A->mnl; i++) {
    len = nnz_per_col[i];
    lines->len[i] = len;
    lines->slen[i] = bs*len;
    if (! symmetric_pattern) {
      rlen = nnz_per_row[i];
      lines->rlen[i] = rlen;
    }
  }

  /* Convert the row and column data to consistent shmalloc buffers */
  col_buf_size = row_buf_size = A_buf_size = 0;
  for (i=0, index=A->my_start_index;
       i<A->mnl;
       i++, index++) {
    col_buf_size += A->lines->len[i];
    /* row structure ? */
    if (A->lines->rptrs)
      row_buf_size += A->lines->rlen[i];
    A_buf_size += (A->block_sizes[index]*A->lines->slen[i]);
  }

#ifdef MPI
  MPI_Barrier(comm);
  MPI_Allreduce(&col_buf_size,&max_col_buf_size,1,
		MPI_INT,MPI_MAX,comm);
  /* row structure ? */
  if (A->lines->rptrs)
    MPI_Allreduce(&row_buf_size,&max_row_buf_size,1,
		  MPI_INT,MPI_MAX,comm);
  MPI_Allreduce(&A_buf_size,&max_A_buf_size,1,
		MPI_INT,MPI_MAX,comm);
#else
  max_col_buf_size = col_buf_size;
  if (A->lines->rptrs)
    max_row_buf_size = row_buf_size;
  max_A_buf_size = A_buf_size;
#endif

#ifdef SHMEM
  A->lines->ptrs_buf  = (int *)    shmalloc(max_col_buf_size*sizeof(int));
  /* row structure ? */
  if (A->lines->rptrs)
    A->lines->rptrs_buf = (int *)    shmalloc(max_row_buf_size*sizeof(int));
  A->lines->A_buf     = (double *) shmalloc(max_A_buf_size*sizeof(double));
#else
  A->lines->ptrs_buf  = (int *)    mmalloc(max_col_buf_size*sizeof(int));
  /* row structure ? */
  if (A->lines->rptrs)
    A->lines->rptrs_buf = (int *)    mmalloc(max_row_buf_size*sizeof(int));
  A->lines->A_buf     = (double *) mmalloc(max_A_buf_size*sizeof(double));
#endif

  /* Set pointers. */
  for (i=0,
	 pos1 = pos2 = 0,
	 rpos = 0;
       i<A->mnl;
       i++) {

    len = nnz_per_col[i];
    lines->A[i] = &(lines->A_buf[pos1]);
    pos1 += bs*bs*len;
    lines->ptrs[i] = &(lines->ptrs_buf[pos2]);
    pos2 += len;

    if (! symmetric_pattern) {
      rlen = nnz_per_row[i];
      lines->rptrs[i] = &(lines->rptrs_buf[rpos]);
      rpos += rlen;
    }

  }

  /* Determine the mapping from global indices to pointers */
  mapping = (int *) mmalloc(A->n*sizeof(int));
  if (! mapping) {
    printf("failed to malloc mapping\n");
    fflush(stdout);
    exit(1);
  }
  pe = 0;
  local_indx = 0;
  for (i=0; i<A->n; i++) {
    if (local_indx >= A->mnls[pe]) {
      pe++;
      local_indx = 0;
    }
    mapping[i] = local_indx + A->start_indices[pe];
    local_indx++;
  }

  /* Install column structure and values. */
  current_col_index = cols[0].jb;
  current_block = -1;
  pos = -1;
  len = 0;
  for (k=0; k<nnz_cols; k++) {

    row = cols[k].ib;
    col = cols[k].jb;
    local_col = col - A->my_start_index;

    if (col != current_col_index) {
      current_block = -1;
      pos = -1;
      len = 0;
      current_col_index = col;
    }

    if (row != current_block) {
      current_block = row;
      pos++;
      lines->ptrs[local_col][pos] = mapping[row];
      len++;
      lines->len[local_col] = len;
    }

    ib = cols[k].i % bs;
    jb = cols[k].j % bs;

    lines->A[local_col][bs*bs*pos + bs*jb + ib] = cols[k].val;

  }

  /* Install row structure? */
  if (! symmetric_pattern) {
    current_row_index = rows[0].ib;
    current_block = -1;
    pos = -1;
    len = 0;
    for (k=0; k<nnz_rows; k++) {

      row = rows[k].ib;
      col = rows[k].jb;
      local_row = row - A->my_start_index;

      if (row != current_row_index) {
	current_block = -1;
	pos = -1;
	len = 0;
	current_row_index = row;
      }

      if (col != current_block) {
	current_block = col;
	pos++;
	lines->rptrs[local_row][pos] = mapping[col];
	len++;
	lines->rlen[local_row] = len;
      }

      ib = rows[k].i % bs;
      jb = rows[k].j % bs;

    }
  }

  A->maxnz = calc_maxnz(A);

  free(mapping);
  free(nnz_per_row);
  free(nnz_per_col);

  return(A);

}

/**********************************************************************/

int rank_column_major(const void *a, const void *b)
{
  mm_data *ma,*mb;

  ma = (mm_data *) a;
  mb = (mm_data *) b;

  if (ma->jb <  mb->jb) return(-1);
  if ((ma->jb == mb->jb) && (ma->ib < mb->ib)) return(-1);
  if ((ma->jb == mb->jb) && (ma->ib == mb->ib)) return(0);
  return(1);

}

/**********************************************************************/

int rank_row_major(const void *a, const void *b)
{
  mm_data *ma,*mb;

  ma = (mm_data *) a;
  mb = (mm_data *) b;

  if (ma->ib <  mb->ib) return(-1);
  if ((ma->ib == mb->ib) && (ma->jb < mb->jb)) return(-1);
  if ((ma->ib == mb->ib) && (ma->jb == mb->jb)) return(0);
  return(1);

}

/**********************************************************************/
/* Read a rhs to match matrix A.
   This must in in Matrix Market "array format"
*/

vector *read_rhs_for_matrix(char *rhs_file, matrix *A)
{
  vector *rhs;
  char line[128];
  int mnl,start_indx,i;
  double val;

  MM_typecode matcode;
  int M,N;

  FILE *f;

  f = fopen(rhs_file,"r");
  if (! f)
    if (A->myid == 0) {
      printf("Couldn't find rhs file: %s\n",rhs_file);
      exit(1);
    }


  if (mm_read_banner(f, &matcode) != 0)
    {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
    }

  if (!
	(mm_is_real(matcode)) &&
	(mm_is_matrix(matcode)) &&
	(mm_is_dense(matcode)) &&
	(mm_is_general(matcode))) {
      printf("MM rhs must be real, dense, and general\n");
      exit(1);
    }

  if (mm_read_mtx_array_size(f, &M, &N)) {
    printf("error in  mm_read_mtx_array_size\n");
    exit(1);
  }

  if (M != A->n) {
    printf("error: rhs does not match matrix size\n");
    printf("M=%d A->n=%d\n",M,A->n);\
    exit(1);
  }

  if (N != 1) {
    printf("error: only one rhs can be read\n");
    exit(1);
  }

  mnl = A->mnl;
  start_indx = A->my_start_index;

  rhs = new_vector(M, mnl);

  for (i=0; i<M; i++) {
    fgets(line,128,f);
    sscanf(line,"%le\n", &val);
    if ((i >= start_indx) && (i < (start_indx+mnl))) {
      rhs->v[i-start_indx] = val;
    }
  }

  fclose(f);

  return(rhs);
}

