/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

/**********************************************************************

com_server is an MPI communications server that mimics the one-sided
communication needed by SPAI.

**********************************************************************/

#include "com_server.h"

#include "debug.h"
#include "timing.h"

int ndone,Im_done,all_done;

/* cache for remotely-referenced lines */
hash_table *ht;

/**********************************************************************/

#ifdef MPI

void com_server
(matrix *A,
 matrix *M)
{
  int tag,flag,requestor;
  MPI_Status status;
  MPI_Request request;
  int ierr;

  start_timer(ident_com_server);

  do {
    ierr = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG,
		      A->comm, &flag, &status);
    if (ierr) {
      printf("error in MPI_Iprobe\n");
      fflush(stdout);
      exit(1);
    }
    if (flag) {

      requestor = status.MPI_SOURCE;
      tag = status.MPI_TAG;

      switch (tag) {

      case get_line_tag:
	/* requesting a line of pointers and values */
	{if (A) handle_get_line(A,requestor); break;}

      case request_Mline_tag:
	/* requesting control of a line of M */
	{if (A) handle_request_Mline(A,requestor); break;}

      case put_Mline_tag:
	/* Someone's storing a line of M */
	{if (M) handle_put_Mline(M,requestor); break;}

      case Im_done_tag:
	/* A PE is done.  PE 0 counts it */
	{if (A) handle_Im_done(A,requestor); break;}

      case done_signal_tag:
	/* Everyone is done */
	{if (A) handle_done_signal(A,requestor); break;}

      case exit_tag:
	/* error exit */
	exit(1);

      }
    }
  }
  while (flag && (tag < 10));

  /* Is everyone done? */
  check_done(A,M);

  stop_timer(ident_com_server);\

}

#else

void com_server
(matrix *A,
 matrix *M)
{}

#endif

/**********************************************************************/

/* The following routines are handlers of requests. */

/**********************************************************************/

/* Someone is requesting a line of pointers and values of A */

#ifdef MPI

void handle_get_line
(matrix *A,
 int requestor)
{
  clines *lines;
  int index,len,rlen,count;
  int *buf,*rbuf;
  double *abuf;
  int block_width;
  MPI_Status status;
  MPI_Request request;
  int ierr;

  lines = A->lines;
  ierr = MPI_Recv((void *) &index, 1, MPI_INT, requestor, get_line_tag,
		  A->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }

  block_width = A->block_sizes[index+A->my_start_index];

  /* determine bufs and count */
  len = lines->len[index];
  count = block_width*lines->slen[index];
  rlen = lines->rlen[index];
  buf = lines->ptrs[index];
  rbuf = lines->rptrs[index];
  abuf = lines->A[index];

  /* Send then synchronously.  The receives have already been posted. */
  ierr = MPI_Send((int *) buf, len, MPI_INT, requestor, sendlines_tag,
		  A->comm);
  if (ierr) {
    printf("error in MPI_Send\n");
    fflush(stdout);
    exit(1);
  }

  /* row structure ? */
  if (A->lines->rptrs) {
    ierr = MPI_Send((int *) rbuf, rlen, MPI_INT, requestor, sendrlines_tag,
		    A->comm);
    if (ierr) {
      printf("error in MPI_Send\n");
      fflush(stdout);
      exit(1);
    }
  }

  ierr = MPI_Send((int *) abuf, count, MPI_DOUBLE,
		  requestor, sendvals_tag,
		  A->comm);
  if (ierr) {
    printf("error in MPI_Send\n");
    fflush(stdout);
    exit(1);
  }

}

#endif

/**********************************************************************/
/* Someone is requesting control of a line of M */
/* This happens because of the load balancing mechanism. */

#ifdef MPI

void handle_request_Mline
(matrix *A,
 int requestor)
{
  int index;
  MPI_Request request;
  MPI_Status status;
  int failure = -1;
  int ierr;

  ierr = MPI_Recv((void *) &index, 1, MPI_INT, requestor, request_Mline_tag,
		  A->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }

  if (next_line < A->mnl) {
    /* Send the next unfinished line. */
    ierr = MPI_Send((int *) &next_line, 1, MPI_INT, requestor, send_Mline_tag,
		    A->comm);
    if (ierr) {
      printf("error in MPI_SEND\n");
      fflush(stdout);
      exit(1);
    }

    next_line++;

  }
  else {

    /* No unfinished lines remain.  Return failure. */
    ierr = MPI_Send((int *) &failure, 1, MPI_INT, requestor, send_Mline_tag,
		    A->comm);
    if (ierr) {
      printf("error in MPI_SEND\n");
      fflush(stdout);
      exit(1);
    }
  }
}

#endif

/**********************************************************************/
/* Someone's storing a line of M */

#ifdef MPI

void handle_put_Mline
(matrix *M,
 int requestor)
{
  int index_and_len[3],index,len,slen;
  int *rptr;
  double *aptr;
  MPI_Status status;
  int b2;
  int ierr;
  int block_width;

  start_timer(ident_handle_put_Mline);

  /* Where does it go and how big is it? */
  ierr = MPI_Recv((void *) index_and_len, 3, MPI_INT,
		  requestor, put_Mline_tag,
		  M->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }
  index = index_and_len[0];
  len = index_and_len[1];
  slen = index_and_len[2];

  /* Allocate the line and receive it. */
  rptr =
    M->lines->ptrs[index] =
    (int *) mmalloc(len*sizeof(int));
  ierr = MPI_Recv((void *) rptr, len, MPI_INT,
		  requestor, put_Mline_inds_tag,
		  M->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }

  block_width = M->block_sizes[M->my_start_index+index];
  aptr =
    M->lines->A[index] =
    (double *) mmalloc(block_width*slen*sizeof(double));
  ierr = MPI_Recv((void *) aptr, block_width*slen, MPI_DOUBLE,
		  requestor, put_Mline_vals_tag,
		  M->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }

  M->lines->len[index] = len;

  stop_timer(ident_handle_put_Mline);

}

#endif

/**********************************************************************/
/* A PE is done.  PE 0 counts it */

#ifdef MPI

void handle_Im_done
(matrix *A,
 int requestor)
{
  int tmp;
  MPI_Status status;
  int ierr;

  start_timer(ident_handle_Im_done);

  ierr = MPI_Recv((void *) &tmp, 1, MPI_INT, requestor, Im_done_tag,
		  A->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }
  ndone++;

  stop_timer(ident_handle_Im_done);

}

#endif

/**********************************************************************/
/* PE 0 is telling me that everyone is done. */

#ifdef MPI

void handle_done_signal
(matrix *A,
 int requestor)
{
  int tmp;
  MPI_Status status;
  int ierr;

  start_timer(ident_handle_done_signal);

  ierr = MPI_Recv((void *) &tmp, 1, MPI_INT, 0, done_signal_tag,
		  A->comm, &status);
  if (ierr) {
    printf("error in MPI_Recv\n");
    fflush(stdout);
    exit(1);
  }
  all_done = 1;

  stop_timer(ident_handle_done_signal);

}

#endif

/**********************************************************************/

/* The following routines are initiators of requests. */

/**********************************************************************/
/* Get a line of A.
   It could be local to my myid or remote.
   If it's remote it may already be in the cache.
   */

void get_line
(matrix  *A,
 matrix  *M,
 int ind,
 int *len_ptr,
 int *rlen_ptr,
 int **buf_ptr,
 int **rbuf_ptr,
 double **vbuf_ptr)
{
 int global_index,send_flag,lines_flag,rlines_flag,vals_flag,index,pe;
 int k;
#ifdef MPI
 MPI_Status send_status,lines_status,rlines_status,vals_status;
 MPI_Request get_lines_request, get_rlines_request,
   get_vals_request, send_request;
#endif
 int *buf,*rbuf;
 double *vbuf;
 int count,rcount,vcount;
 int ierr;
 int offset;

 start_timer(ident_get_line);

 pe = A->pe[ind];
 index = ind - A->start_indices[pe];

 /* Is the line local? */
 if (pe == A->myid) {

   *len_ptr = A->lines->len[index];
   *buf_ptr = A->lines->ptrs[index];

   /* row structure ? */
   if (A->lines->rptrs) {
     *rlen_ptr = A->lines->rlen[index];
     *rbuf_ptr = A->lines->rptrs[index];
   }
   else {
     *rlen_ptr = A->lines->len[index];
     *rbuf_ptr = A->lines->ptrs[index];
   }

   *vbuf_ptr = A->lines->A[index];
 }

 else {

   /* look in the cache */
   global_index = ind;
   if (lookup(ht, global_index,
	      buf_ptr, rbuf_ptr, vbuf_ptr,
	      &count, &rcount)) {
     *len_ptr = count;
     *rlen_ptr = rcount;
   }

   else {

     /* It's not there.  Get it from a remote processor */

     count = len_all[ind];
     rcount = rlen_all[ind];
     vcount = A->block_sizes[ind]*slen_all[ind];

     #if (defined(MPI) && !defined(SHMEM))

     ierr = MPI_Irecv((void *) remote_buf,
		      A->maxnz, MPI_INT, pe, sendlines_tag,
		      A->comm, &get_lines_request);
     if (ierr) {
       printf("error in MPI_Irecv\n");
       fflush(stdout);
       exit(1);
     }

     /* row structure ? */
     if (A->lines->rptrs) {
       ierr = MPI_Irecv((void *) remote_rbuf,
			A->maxnz, MPI_INT, pe, sendrlines_tag,
			A->comm, &get_rlines_request);
       if (ierr) {
	 printf("error in MPI_Irecv\n");
	 fflush(stdout);
	 exit(1);
       }

     }

     ierr = MPI_Irecv((void *) remote_abuf,
		      A->maxnz, MPI_DOUBLE, pe, sendvals_tag,
		      A->comm, &get_vals_request);
     if (ierr) {
       printf("error in MPI_Irecv\n");
       fflush(stdout);
       exit(1);
     }

     /* Send the request. */
     ierr = MPI_Isend((void *) &index, 1, MPI_INT, pe, get_line_tag,
		      A->comm, &send_request);
     if (ierr) {
       printf("error in MPI_Isend\n");
       fflush(stdout);
       exit(1);
     }
     do {
       com_server(A,M);
       ierr = MPI_Test(&send_request, &send_flag, &send_status);
       if (ierr) {
	 printf("error in MPI_Test\n");
	 fflush(stdout);
	 exit(1);
       }
     }
     while (! send_flag);

     /* Service requests until the data comes back. */
     do {
       com_server(A,M);
       ierr = MPI_Test(&get_lines_request, &lines_flag, &lines_status);
       if (ierr) {
	 printf("error in MPI_Test\n");
	 fflush(stdout);
	 exit(1);
       }
     }
     while (! lines_flag);

     /* row structure ? */
     if (A->lines->rptrs) {
       do {
	 com_server(A,M);
	 ierr = MPI_Test(&get_rlines_request, &rlines_flag, &rlines_status);
	 if (ierr) {
	   printf("error in MPI_Test\n");
	   fflush(stdout);
	   exit(1);
	 }
       }
       while (! rlines_flag);
     }

     do {
       com_server(A,M);
       ierr = MPI_Test(&get_vals_request, &vals_flag, &vals_status);
       if (ierr) {
	 printf("error in MPI_Test\n");
	 fflush(stdout);
	 exit(1);
       }
     }
     while (! vals_flag);

     *len_ptr = count;
     *buf_ptr = remote_buf;

     /* row structure ? */
     if (A->lines->rptrs) {
       *rlen_ptr = rcount;
       *rbuf_ptr = remote_rbuf;
     }
     else {
       *rlen_ptr = count;
       *rbuf_ptr = remote_buf;
     }

     *vbuf_ptr = remote_abuf;

#endif

#if defined(SHMEM)

     pe = A->pe[global_index];

     offset = ptr_offsets[global_index];
     shmem_int_get(remote_buf,&(A->lines->ptrs_buf[offset]),
		   count,pe);

     /* row structure ? */
     if (A->lines->rptrs) {
       offset = rptr_offsets[global_index];
       shmem_int_get(remote_rbuf,&(A->lines->rptrs_buf[offset]),
		   rcount,pe);
     }

     offset = A_offsets[global_index];
     shmem_double_get(remote_abuf,&(A->lines->A_buf[offset]),
		      vcount,pe);

     *len_ptr = count;
     *buf_ptr = remote_buf;

     /* row structure ? */
     if (A->lines->rptrs) {
       *rlen_ptr = rcount;
       *rbuf_ptr = remote_rbuf;
     }
     else {
       *rlen_ptr = count;
       *rbuf_ptr = remote_buf;
     }

     *vbuf_ptr = remote_abuf;

#endif

     if (ht) {  /* Does the hash table exist? */

       /* Cache the line */

       /* row structure ? */
       if (A->lines->rptrs) {
	 insert(ht, global_index,
		remote_buf,remote_rbuf,remote_abuf,
		count,rcount,vcount);
       }
       else {
	 insert(ht, global_index,
		remote_buf,remote_buf,remote_abuf,
		count,count,vcount);
       }

     }

   }
 }

 stop_timer(ident_get_line);

}

/**********************************************************************/
/* Get a line of M. */

#define locked 1
#define not_locked 0

int lock = not_locked;


int request_Mline
(matrix *A,
 matrix *M,
 int pe,
 int *index_ptr)
{

#if (defined(MPI) && !defined(SHMEM))
  int one = 1;
  int index;
  MPI_Request recv_request, send_request;
  MPI_Status status;
  int flag;
  int ierr;
#endif

#if defined(SHMEM)
  int status,new_next_line;
  int index;
#endif

  int success = 0;

  start_timer(ident_request_Mline);

#if (defined(MPI) && !defined(SHMEM))

  ierr = MPI_Irecv((void *) &index, 1, MPI_INT, pe, send_Mline_tag,
		   A->comm, &recv_request);
  if (ierr) {
    printf("error in MPI_Irecv\n");
    fflush(stdout);
    exit(1);
  }

  /* Send the request. */
  ierr = MPI_Isend((void *) &one, 1, MPI_INT, pe, request_Mline_tag,
		   A->comm, &send_request);
  if (ierr) {
    printf("error in MPI_Isend\n");
    fflush(stdout);
    exit(1);
  }

  /* Service requests until the data comes back. */
  do {
    com_server(A,M);
    ierr = MPI_Test(&recv_request, &flag, &status);
    if (ierr) {
      printf("error in MPI_Test\n");
      fflush(stdout);
      exit(1);
    }
  }
  while (! flag);
  if (index >= 0) {
    *index_ptr = index;
    success = 1;
  }
  else {
    success = 0;
  }

#endif

#if defined(SHMEM)

  /* lock pe */
  for (;;) {
    status =
      shmem_int_cswap(&lock, not_locked, locked, pe);
    if (status == locked) break;
    com_server(A,M);
  }

  shmem_int_get(index_ptr, &next_line, 1, pe);

  if (*index_ptr < A->mnls[pe]) {

    /* update next_line */
    new_next_line = (*index_ptr) + 1;
    shmem_int_put(&next_line, &new_next_line, 1, pe);

    /* unlock pe */
    shmem_int_cswap(&lock, locked, not_locked, pe);

    /* return success */
    success = 1;

  }
  else {

    /* unlock pe */
    shmem_int_cswap(&lock, locked, not_locked, pe);

    /* return failure */
    success = 0;

  }

#endif

  stop_timer(ident_request_Mline);
  return(success);

}

/**********************************************************************/
/* Store a line of M. */

#ifdef MPI

void put_Mline
(matrix *A,
 matrix *M,
 int add,
 int *inds,
 double *vals,
 int len,
 int slen)
{
  int *rptr;
  double *aptr;
  int i,flag,index_and_len[3];
  MPI_Request requests[3];
  MPI_Status statuses[3];
  int ierr;
  int pe,index,block_width;

  start_timer(ident_put_Mline);

  pe = A->pe[add];
  index = add-A->start_indices[pe];
  block_width = block_size(add,A);

  if (pe == A->myid) {

    /* The line is local. Just allocate it and copy the data. */
    aptr =
      M->lines->A[index] =
      (double *) mmalloc(block_width*slen*sizeof(double));
    M->lines->A[index] = aptr;

    rptr =
      M->lines->ptrs[index] =
      (int *) mmalloc(len*sizeof(int));

    for (i=0; i<len; i++) {
      rptr[i] = inds[i];
    }
    for (i=0; i<block_width*slen; i++) {
      aptr[i] = vals[i];
    }
    M->lines->len[index] = len;

  }

  else {

    /* The line is remote. */
    index_and_len[0] = index;
    index_and_len[1] = len;
    index_and_len[2] = slen;

    /* Send the request. */
    ierr = MPI_Isend((void *) index_and_len, 3, MPI_INT,
		     pe, put_Mline_tag,
		     A->comm, &requests[0]);
    if (ierr) {
      printf("error in MPI_Isend\n");
      fflush(stdout);
      exit(1);
    }

    /* Send the data. */
    ierr = MPI_Isend((void *) inds, len, MPI_INT, pe,
		     put_Mline_inds_tag,
		     A->comm, &requests[1]);
    if (ierr) {
      printf("error in MPI_Isend\n");
      fflush(stdout);
      exit(1);
    }
    ierr = MPI_Isend((void *) vals, block_width*slen, MPI_DOUBLE, pe,
		     put_Mline_vals_tag,
		     A->comm, &requests[2]);
    if (ierr) {
      printf("error in MPI_Isend\n");
      fflush(stdout);
      exit(1);
    }
    do {
      com_server(A,M);
      ierr = MPI_Testall(3,requests,&flag,statuses);
      if (ierr) {
	printf("error in MPI_Testall\n");
	fflush(stdout);
	exit(1);
      }
    }
    while (! flag);
  }

  stop_timer(ident_put_Mline);

}

#else

void put_Mline
(matrix *A,
 matrix *M,
 int add,
 int *inds,
 double *vals,
 int len,
 int slen)
{
  int *rptr;
  double *aptr;
  int pe,index,block_width,i;

  pe = A->pe[add];
  index = add-A->start_indices[pe];
  block_width = block_size(add,A);

  aptr = (double *) mmalloc(block_width*slen*sizeof(double));
  M->lines->A[index] = aptr;

  rptr =
    M->lines->ptrs[index] =
    (int *) mmalloc(len*sizeof(int));

  for (i=0; i<len; i++) {
    rptr[i] = inds[i];
  }
  for (i=0; i<block_width*slen; i++) {
    aptr[i] = vals[i];
  }
  M->lines->len[index] = len;
}

#endif

/**********************************************************************/
/* I'm finished with my lines.  Tell PE 0 */

#ifdef MPI

void say_Im_done
(matrix *A,
 matrix *M)
{
  int one = 1;
  MPI_Request request;
  int ierr;
  int flag;
  MPI_Status status;

  start_timer(ident_say_Im_done);

  Im_done = 1;
  if (A->myid == 0) {
    ndone++;
  }
  else {
    ierr = MPI_Isend((void *) &one, 1, MPI_INT, 0, Im_done_tag,
		     A->comm, &request);
    if (ierr) {
      printf("error in MPI_Isend\n");
      fflush(stdout);
      exit(1);
    }
    do {
      com_server(A,M);
      ierr = MPI_Test(&request,&flag,&status);
      if (ierr) {
	printf("error in MPI_Test\n");
	fflush(stdout);
	exit(1);
      }
    }
    while (! flag);
  }

  stop_timer(ident_say_Im_done);

}

#else

void say_Im_done
(matrix *A,
 matrix *M)
{

  Im_done = 1;
  if (A->myid == 0) {
    ndone++;
  }
}

#endif


/**********************************************************************/
/* PE 0 checks whether everyone is done.  If so, tell everyone. */

#ifdef MPI

void check_done
(matrix *A,
 matrix *M)
{
  int i,one = 1;
  MPI_Request request;
  int ierr;
  int flag;
  MPI_Status status;

  start_timer(ident_check_done);

  if (A->myid==0)
    if ((! all_done) && (ndone == A->numprocs))
      {
	for (i=1; i<A->numprocs; i++) {
	  ierr = MPI_Isend((void *) &one, 1, MPI_INT, i, done_signal_tag,
			   A->comm, &request);
	  if (ierr) {
	    printf("error in MPI_Isend\n");
	    fflush(stdout);
	    exit(1);
	  }
	}

	all_done = 1;
      }

  stop_timer(ident_check_done);

}

#else

void check_done
(matrix *A,
 matrix *M)

{}

#endif


