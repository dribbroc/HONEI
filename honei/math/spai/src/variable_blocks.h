#ifndef __variable_blocks_H
#define __variable_blocks_H

#include "basics.h"
#include "matrix.h"
#include "timing.h"

#include "debug.h"

matrix *block_matrix
(matrix *,
 int,
 int,
 int);

matrix *constant_block_matrix
(matrix *,
 int);

void find_constant_blocks
(matrix *,
 int,
 int **,
 int *);

matrix *variable_block_matrix
(matrix *,
 int);

void find_diagonal_blocks
(matrix *,
 int,
 int **,
 int *);

int find_diagonal_block
(matrix *,
 int,
 int);

int initial_run_length
(matrix *,
 int);

int check_next_run
(matrix *,
 int,
 int,
 int);

matrix *convert_to_block_matrix
(matrix *,
 int ,
 int *);

matrix *scalar_matrix
(matrix *,
 int);

#endif
