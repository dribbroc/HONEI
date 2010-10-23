/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __com_server_H
#define __com_server_H

#include "basics.h"
#include "vector.h"
#include "index_set.h"
#include "matrix.h"
#include "spai.h"
#include "load_balance.h"
#include "hash_table.h"

extern int ndone;
extern int Im_done;
extern int all_done;

extern hash_table *ht;

/* Messages with these tags are handled by com_server */
#define get_line_tag         1
#define request_Mline_tag    2
#define put_Mline_tag        3
#define Im_done_tag          4
#define done_signal_tag      5
#define exit_tag             6

/* Messages with these tags aren't */
#define sendlines_tag       10
#define sendrlines_tag      11
#define sendvals_tag        12
#define send_Mline_tag      13
#define put_Mline_inds_tag  14
#define put_Mline_vals_tag  15


void com_server
(matrix *,
 matrix *);

void handle_get_line
(matrix *,
 int);

void handle_request_Mline
(matrix *,
 int);

void handle_put_Mline
(matrix *,
 int);

void handle_request_alines
(matrix *,
 int);

void handle_Im_done
(matrix *,
 int);

void handle_done_signal
(matrix *,
 int);

void get_line
(matrix *,
 matrix *,
 int,
 int *,
 int *,
 int **,
 int **,
 double **);

int request_Mline
(matrix *,
 matrix *,
 int,
 int *);

void put_Mline
(matrix *,
 matrix *,
 int,
 int *,
 double *,
 int,
 int);

void say_Im_done
(matrix *,
 matrix *);

void check_done
(matrix *,
 matrix *);

#endif
