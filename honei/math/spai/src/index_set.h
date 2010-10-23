/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __index_set_H
#define __index_set_H

#include "basics.h"

struct index_set {
  int len;  /* number of blocks */
  int slen; /* scalar length (sum of block heights) */
  int *ptr;
};

typedef struct index_set index_set;

index_set *new_index_set
(index_set *,
 int,
 char *);

void free_index_set
(index_set *);

void write_index_set
(FILE *,
 int *,
 index_set *);

#endif
