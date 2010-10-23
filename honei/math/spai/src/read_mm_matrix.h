
/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __read_matrix_mm_H
#define __read_matrix_mm_H

#include "basics.h"
#include "matrix.h"
#include "mmio.h"
#include "timing.h"

struct mm_data {
  int i;
  int j;
  double val;
  int ib;
  int jb;
};

typedef struct mm_data mm_data;

matrix *read_mm_matrix
(const char *,
 int,
 int,
 int,
 int,
 int,
 int,
 SPAI_Comm comm);

void skip_mm_matrix_header_ascii
(FILE *);

void skip_mm_matrix_header_binary
(FILE *);

void read_mm_matrix_header_ascii
(FILE *,
 int *M,
 int *N,
 int *nnz);

void read_mm_matrix_header_binary
(FILE *,
 int *M,
 int *N,
 int *nnz);

void write_mm_data
(FILE *,
 mm_data *,
 int,
 mm_data *,
 int);

void read_mm_data_ascii
(FILE *,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 mm_data **,
 int *,
 mm_data **,
 int *);

void read_mm_data_binary
(FILE *,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 mm_data **,
 int *,
 mm_data **,
 int *);

matrix *mm_to_matrix
(int,
 int,
 int,
 int,
 int,
 mm_data *,
 int,
 mm_data *,
 int,
 SPAI_Comm comm);


int rank_column_major
(const void *a,
 const void *b);

int rank_row_major
(const void *a,
 const void *b);

vector *read_rhs_for_matrix
(char *rhs_file,
 matrix *A);

#endif

