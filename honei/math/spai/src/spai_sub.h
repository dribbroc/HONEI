/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __spai_sub_H
#define __spai_sub_H

#include "basics.h"
#include "vector.h"
#include "index_set.h"
#include "blocks.h"
#include "matrix.h"
#include "spai.h"
#include "com_server.h"

int scalar_len
(index_set *,
 matrix *);

int precompute_column_square_inverses
(matrix *);

int precompute_column_square_inverse1
(double *,
 int *,
 int,
 int,
 matrix *,
 double *,
 int *,
 double *);

void getrows
(matrix *,
 matrix *,
 index_set *,
 index_set *);

void getcols
(matrix *,
 matrix *,
 index_set *,
 index_set *);

void full_matrix
(matrix *,
 matrix *,
 int,
 double *);

void write_full_matrix
(FILE *,
 double *,
 int,
 int,
 int);

int augment_sparsity
(matrix *,
 matrix *,
 int,
 int,
 double);

void select_new_nonzeros
(double,
 matrix *);

void copyvv
(index_set *,
 index_set *);

int append
(index_set *,
 index_set *);

double innprod
(double *,
 double *,
 int);

double frobenius_norm
(double *,
 int,
 int);

void deleter
(index_set *,
 index_set *,
 matrix *);

double calcrs
(matrix *,
 matrix *,
 int,
 double);

int block_size
(int,
 matrix *);

void write_line
(FILE *,
 matrix *,
 int,
 int,
 int,
 int *,
 int *,
 double *);

#endif
