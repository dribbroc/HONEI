/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "load_balance.h"

int next_line;


/**********************************************************************/
/* Returns the index of the next column to process.
   or, returns -1 if unsuccessful in finding a column */

int grab_Mline
(matrix *A,
 matrix *M,
 SPAI_Comm comm)
{
  int i,pe,index,col;

  /* Do I still have local columns to process? */
  if (next_line < A->mnl) {

    pe = A->myid;
    index = next_line;
    col = index + A->start_indices[pe];
    next_line++;

    return(col); /* success */
  }

#ifdef MPI
  else {  /* look for a column on another pe */
    for (i=1; i<A->numprocs; i++) {
      pe = (A->myid+i) % A->numprocs;
      if (request_Mline(A, M, pe, &index)) {
	col = index + A->start_indices[pe];
	return(col); /* success */
      }
    }
  }
#endif

  /* No one has any more columns */
  return(-1);  /* failure */
}

/**********************************************************************/







