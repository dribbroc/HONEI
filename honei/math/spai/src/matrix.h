/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __matrix_H
#define __matrix_H

#include "basics.h"
#include "vector.h"
#include "index_set.h"

#include "debug.h"

extern double *temp_block;

/* ********************************************************************

   This struct is used to represent both the rows and columns
   of a distributed nonsymmetric sparse matrix.

   The nonzeros are stored in compressed form.

******************************************************************** */

struct compressed_lines {

  /* Pointers to access the coefficients (by column) */
  double **A;
  /* Buffer for column-wise coefficients*/
  double *A_buf;

  /* Pointers to access the distributed indices */
  int **ptrs;
  /* Buffer for distributed indices */
  int *ptrs_buf;

  /* row indices */
  int **rptrs;
  int *rptrs_buf;

  /* The number of nonzeros in the columns */
  int *len;

  /* The scalar length of the columns */
  int *slen;

  /* The number of nonzeros in the rows */
  int *rlen;

  /* Adding stuff for variable-length blocks ???? */
  int *block_sizes;

};

typedef struct compressed_lines clines;

struct mv_buf {
  int len;
  int *inds;    /* local indices */
  double *vals;
};

typedef struct mv_buf mv_buf;


struct matrix {

  SPAI_Comm comm;
  int myid;
  int numprocs;
  int mnl;
  int my_start_index;

  int transposed;     /* 0 => distibution by columns => right preconditioenr
			 1 => distribution by rows => left preconditioner*/

  int n;              /* number of rows/cols */

  int bs;             /* block size (0 => variable) */
  int max_block_size; /* largest block size */
  int maxnz;          /* maximum number of scalar nonzeros in a column (or row) */
  int tot_nnz;        /* total number of non-zero elements */

  int *mnls;          /* local n in every pe */
  int *start_indices; /* starting index in every pe */

  int my_nnz;         /* nonzeros in this pe (column-wise) */

  clines *lines;      /* columns and rows */

  int *pe;            /* processor assignment for *every* row and column */

  int *block_sizes;   /* The size of every diagonal block */

};

typedef struct matrix matrix;

void init_matrix
(matrix *,
 SPAI_Comm);

matrix *new_matrix
(SPAI_Comm);

void init_compressed_lines
(clines *);

clines *new_compressed_lines
(int,
 int);

void free_compressed_lines
(matrix *);

matrix *clone_matrix
(matrix *);

void sp_free_matrix
(matrix *);

void order_pointers
(matrix *);

void orderv
(double *,
 int *,
 int,
 matrix *,
 int,
 int);

void A_times_v_cc
(matrix *,
 vector *,
 vector *);

void A_times_v_rc
(matrix *,
 vector *,
 vector *);

void gather_vec
(matrix *,
 vector *,
 double *);

int count_nonzeros
(matrix *);

void write_matrix_mm
(matrix *,
 char *,
 int);

void matrix_statistics
(matrix *,
 char *);

double *inverse_diagonal
(matrix *A);

void square_block_inverse
(int,
 double *,
 double *,
 int *,
 double *);

int block_size
(int,
 matrix *);

int calc_maxnz
(matrix *);

#endif
