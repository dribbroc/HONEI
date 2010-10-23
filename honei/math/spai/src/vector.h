/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __vector_H
#define __vector_H

#include "basics.h"
#include "mmio.h"
#include "debug.h"

struct vector {
  int n;          /* number of elements */
  int mnl;        /* my local n */
  double *v;      /* values */
};

typedef struct vector vector;

vector *new_vector
(int,
 int);

void free_vector
(vector *);

void rzeros
(vector *);

void v_plus_cw
(vector *,
 vector *,
 double,
 vector *);

void rcopy_vv
(vector *,
 vector *);

double norm
(vector *,
 SPAI_Comm comm);

double dot
(vector *,
 vector *,
 SPAI_Comm);

void write_vector
(FILE *,
 vector *);

vector *uniform_vector
(int,
 int,
 double);

void write_vector_mm
(vector *,
 char *,
 SPAI_Comm comm);

#endif
