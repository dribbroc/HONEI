/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __spai_H
#define __spai_H

#include "basics.h"
#include "vector.h"
#include "blocks.h"
#include "index_set.h"
#include "matrix.h"
#include "bicgstab.h"
#include "spai_sub.h"
#include "spai_error.h"
#include "qr.h"
#include "hash_table.h"
#include "com_server.h"
#include "load_balance.h"

/* parameters (external, so known to spai_sub) */
extern FILE *message;
extern double epsilon;
extern int nbsteps;
extern int max_dim;
extern int maxnew;

extern int maxapi;

extern index_set *J,*I,*J_tilde,*I_tilde,*is;
extern index_set *isres;
extern double *res,*resb;

extern int *remote_buf;
extern int *remote_rbuf;
extern double *remote_abuf;

extern int *n1,*n2;

extern int Ahat_size;
extern int TAU_size;
extern int R_size;
extern int Z_size;

extern int TAU_ub;
extern int R_ub;
extern int Z_ub;

#define init_val 1.0e-300

extern int *TAU_ptr;

extern double *Ahat;
extern double *R;
extern double *TAU;
extern double *Z;

extern double *rw,*rs,*rz;
extern double *x;
extern double *xb;

extern double **Qlist;

extern double *temp_block;
extern double *minus_Bj;
extern double *sum_block;

extern double *si_buffer;
extern int *si_indices;

extern int nbits_of_int;
extern int log_nbits_of_int;

extern unsigned int *bitvec;

extern int *len_all;
extern int *rlen_all;
extern int *slen_all;
extern int **ptr_addresses;
extern int **rptr_addresses;
extern double **A_addresses;
extern int *ptr_offsets;
extern int *rptr_offsets;
extern int *A_offsets;

#ifdef SHMEM
extern int **ptr_addresses;
extern int **rptr_addresses;
extern double **A_addresses;
#endif

int bspai
(matrix *, matrix **,
 FILE *,
 double,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 double);

int spai
(matrix *, matrix **,
 FILE *,
 double,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 int,
 double);

int spai_line
(matrix *,
 int,
 int,
 int,
 int,
 double,
 matrix *);

void allocate_globals
(matrix *);

void init_double_array
(double *,
 int,
 double);

int amount_touched
(double *,
 int,
 double);

void free_globals
();

int *new_int_array
(int *,
 int,
 char *);

double *new_double_array
(double *,
 int,
 char *);

#endif
