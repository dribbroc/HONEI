/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "index_set.h"
#include "spai_sub.h"

/**********************************************************************/

index_set *new_index_set(index_set *s, int n, char *str)
{

  if (debug) {
    fprintf(fptr_dbg,"new_index_array allocating %d for %s\n",
	    n,str);
  }
  if (! s) {
    s = (index_set *) mmalloc(sizeof(index_set));
    s->len     = 0;
    s->slen     = 0;
    s->ptr = NULL;
  }
  s->ptr = (int *) mmalloc(sizeof(int)*n);

  return(s);
}

/**********************************************************************/

void free_index_set(index_set *s)
{
  if (s->ptr) free(s->ptr);
  free(s);
}

/**********************************************************************/

void write_index_set(FILE *fptr, int *start_indices, index_set *s)
{
  int i;

  for (i=0; i<s->len; i++) {
    fprintf(fptr,"i=%d index=%d\n",
	    i, s->ptr[i]);
  }

  fflush(fptr);
}

