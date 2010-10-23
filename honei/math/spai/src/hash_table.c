/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "hash_table.h"

#define jump 5

/**********************************************************************/

hash_table *init_hash_table(int num)
{
  int i;
  hash_table *ht;

  /* Check that size and jump are relative primes */
  if (! (num % jump)) {
    fprintf(stderr,"Warning in hash_tables:init:\n");
    fprintf(stderr,"   num and jump are not relative primes\n");
  }

  ht = (hash_table *) mmalloc(sizeof(hash_table));

  ht->size = num;

  ht->index_table =
    (int **) mmalloc(sizeof(int *)*num);
  if (! ht->index_table) {
    fprintf(stderr,"failed to allocate index_table\n");
    exit(1);
  }

  for (i=0; i<ht->size; i++) ht->index_table[i] = NULL;

  ht->rindex_table =
    (int **) mmalloc(sizeof(int *)*num);
  if (! ht->index_table) {
    fprintf(stderr,"failed to allocate rindex_table\n");
    exit(1);
  }

  for (i=0; i<ht->size; i++) ht->rindex_table[i] = NULL;

  ht->vals_table = (double **) mmalloc(sizeof(double *)*num);
  if (! ht->vals_table) {
    fprintf(stderr,"failed to allocate vals_table\n");
    exit(1);
  }

  for (i=0; i<ht->size; i++) ht->vals_table[i] = NULL;

  return(ht);

}

/**********************************************************************/

void free_hash_table(hash_table *ht)
{
  int i;

  if (ht) {

    for (i=0; i<ht->size; i++) {
      if (ht->index_table[i]) {
	free(ht->index_table[i]);
      }
      if (ht->rindex_table[i]) {
	free(ht->rindex_table[i]);
      }
      if (ht->vals_table[i]) {
	free(ht->vals_table[i]);
      }
    }

    free(ht->index_table);
    free(ht->rindex_table);
    free(ht->vals_table);
    free(ht);

  }

}

/**********************************************************************/

int insert
(hash_table *ht,
 int global_index,
 int *index_buf,
 int *rindex_buf,
 double *vals_buf,
 int len,
 int rlen,
 int vcount)
{
 int *ibuf;
 int *ribuf;
 double *vbuf;
 int loc;

 if (! ht) return(0);

 /* Find the location */
 loc = findloc(ht,global_index);

 ibuf = (int *) mmalloc(sizeof(int)*(len+2));

 ibuf[0] = global_index;
 ibuf[1] = len;
 memcpy((void *) &ibuf[2], (void *) index_buf, len*sizeof(int));
 ht->index_table[loc] = ibuf;

 ribuf = (int *) mmalloc(sizeof(int)*(rlen+2));
 if (! ribuf) {
   fprintf(stderr,"failed to allocate buffer for rindex table.\n");
   exit(1);
 }

 ribuf[0] = global_index;
 ribuf[1] = rlen;
 memcpy((void *) &ribuf[2], (void *) rindex_buf, rlen*sizeof(int));
 ht->rindex_table[loc] = ribuf;

 vbuf = (double *) mmalloc(sizeof(double)*vcount);
 if (! vbuf) {
   fprintf(stderr,"failed to allocate buffer for vals table.\n");
   exit(1);
 }

 memcpy((void *) vbuf, (void *) vals_buf, vcount*sizeof(double));
 ht->vals_table[loc] = vbuf;

 return(1);

}

/**********************************************************************/

int lookup
(hash_table *ht,
 int global_index,
 int **index_buf,
 int **rindex_buf,
 double **vals_buf,
 int *len,
 int *rlen)
{
  int loc;

  if (! ht) return(0);

  if ((loc=getloc(ht,global_index)) > -1) {

    /* found them */
    *index_buf = ht->index_table[loc];
    *rindex_buf = ht->rindex_table[loc];
    *vals_buf = ht->vals_table[loc];

    *len = (*index_buf)[1];
    (*index_buf)++;
    (*index_buf)++;

    *rlen = (*rindex_buf)[1];
    (*rindex_buf)++;
    (*rindex_buf)++;

    return(1);
  }

  else {

    return(0);

  }
}

/**********************************************************************/
/* Finds the location in which to make an insertion.
   This will attempt at most 5 jumps of linear rehashing before
   returning a location.  In that case, whatever was already in the
   location will be discarded
*/

int findloc
(hash_table *ht,
 int global_index)
{
  unsigned long loc;
  int ntries;
  int *buf;
  int size;

  size = ht->size;

  loc = global_index % size;
  ntries = 1;
  buf = ht->index_table[loc];

  do {
    if (buf == NULL) {
      return loc;
    }
    else {
      /* linear rehash */
      loc += jump;
      loc %= size;
      ntries++;
      buf = ht->index_table[loc];
    }
  }
  while (ntries < 5);

  /* The 5 jumps are up.  Delete what's there and return loc. */
  if (! (buf == NULL)) {
    free(buf);
    free(ht->rindex_table[loc]);
    free(ht->vals_table[loc]);
  }

  return(loc);
}

/**********************************************************************/

int getloc
(hash_table *ht,
 int global_index)
{
  unsigned long loc;
  int ntries;
  int *buf;
  int size;

  size = ht->size;
  loc = global_index % size;
  ntries = 1;

  do {
    buf = ht->index_table[loc];
    if (buf == NULL) return -1;
    if (buf[0] == global_index) {
      return loc;
    }
    else {
      /* linear rehash */
      loc += jump;
      loc %= size;
      ntries++;
    }
  }
  while (ntries < size);

  return(-1);
}






