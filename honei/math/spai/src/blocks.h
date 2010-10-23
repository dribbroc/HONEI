/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __blocks_H
#define __blocks_H

#include "basics.h"
#include "index_set.h"
#include "matrix.h"

void write_block
(FILE *,
 double *,
 int,
 int);

void square_line
(double *,
 int *,
 int,
 matrix *,
 int,
 double *);

void zero_block
(double *,
 int,
 int);

void mult_blocks
(double *,
 double *,
 int,
 int,
 int,
 double *);

void mult_blocks_TN
(double *,
 double *,
 int,
 int,
 int,
 double *);

void convert_to_block
(double *,
 double *,
 int,
 int *,
 matrix *,
 int,
 int);

#endif
