/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "basics.h"
#include "timing.h"
#include "debug.h"

/**********************************************************************/

void init_SPAI()
{

#ifdef SHMEM
  shmem_udcflush();
  shmem_set_cache_inv();
#endif

  init_timers();

}

/**********************************************************************/
/* Determines the distibution of n objects over all PEs.

input:
   n:          order of the matrix
   chunk_size: all pes (except possibly the last) have mnl as a
               multiple of chuck size (used when the matrix is to
	       be converted to block form)

output:
   nl:         maximum number of local objects (same for all pes)
   mnl:        my number of local objects
   split_pe:   the first pe for which mnl \= nl
               (0 if numprocs divides n evenly)
   split_indx: the first (global) index in split_pe
               (All indices are zero-based.)
   start_indx: the starting index in this pe.
*/

void basic_distribution
  (SPAI_Comm comm,
   int n,
   int chunk_size,
   int *nl,
   int *mnl,
   int *split_pe,
   int *split_indx,
   int *start_indx)
{
  int flr,cng,rem,the_mnl,the_partial_sum_mnl;
  int nc;

  int numprocs,myid;
#ifdef MPI
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myid);
  MPI_Barrier(comm);
#else
  numprocs = 1;
  myid = 0;
#endif

#ifdef MPI

  nc = ceil( (double) n/chunk_size);

  flr = floor( (double) nc/numprocs);
  cng = ceil( (double) nc/numprocs);
  rem = fmod( (double) nc, (double) numprocs);
  if (myid < rem) the_mnl = cng*chunk_size;
  else            the_mnl = flr*chunk_size;

  MPI_Barrier(comm);
  MPI_Scan
    (&the_mnl, &the_partial_sum_mnl, 1,
     MPI_INT, MPI_SUM, comm);

  *nl = cng*chunk_size;
  *mnl = the_mnl;
  *split_pe = rem;
  *split_indx = chunk_size*cng*rem;
  *start_indx = the_partial_sum_mnl - the_mnl;

  /* Adjust mnl in last processor
     if n is not a multiple of chuck size */
  if (myid == (numprocs-1))
    if (the_partial_sum_mnl != n)
      *mnl = the_mnl - (the_partial_sum_mnl - n);

#else

  *nl = n;
  *mnl = n;
  *split_pe = 0;
  *split_indx = 0;
  *start_indx = 0;

#endif

}

/**********************************************************************/

void* mmalloc(size_t size)
{
  void* p = malloc(size);
  if (p == NULL)
    {
     printf("fatal error: malloc failed to grant request of size %d\n", size);
     fflush(stdout);
    }
  return p;
}
