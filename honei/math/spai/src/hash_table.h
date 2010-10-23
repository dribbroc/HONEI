/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __hash_table_H
#define __hash_table_H

#include "basics.h"
#include "index_set.h"

#include "debug.h"

struct hash_table {

  int size;

  int    **index_table;
  int    **rindex_table;
  double      **vals_table;

};

typedef struct hash_table hash_table;


/* public */

hash_table *init_hash_table
(int);

void free_hash_table
(hash_table *ht);

int insert
(hash_table *,
 int,
 int *,
 int *,
 double *,
 int,
 int,
 int);

int lookup
(hash_table *,
 int,
 int **,
 int **,
 double **,
 int *,
 int *);

int findloc
(hash_table *,
 int index);

int getloc
(hash_table *,
 int index);

#endif
