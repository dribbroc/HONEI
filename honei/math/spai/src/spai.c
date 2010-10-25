/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/


#include "spai.h"
#include "variable_blocks.h"
#include "debug.h"
#include "timing.h"

/* parameters (external, so known to spai_sub and qr) */
FILE *message;        /* file for warning messages */
double epsilon;       /* tolerance */
int nbsteps;          /* max number of "improvement" steps per line */
int max_dim;              /* max dimensions of I, q, etc. */
int maxnew;           /* max number of new entries per step */

int maxapi;           /* maximum number of nonzeros for any column of M */

/* index sets */
index_set *J,*I,*J_tilde,*I_tilde,*is;
index_set *isres;

/* residual */
double *res;
double *resb;

int *remote_buf;
int *remote_rbuf;
double *remote_abuf;

int *n1,*n2;

int Ahat_size;
int TAU_size;
int R_size;
int Z_size;

int *TAU_ptr;

int TAU_ub;
int R_ub;
int Z_ub;

double *Ahat;
double *R;
double *TAU;
double *Z;

double *rw,*rs,*rz;
double *x;
double *xb;

double **Qlist;

double *temp_block;
double *minus_Bj;
double *sum_block;

double *si_buffer=NULL;
int *si_indices=NULL;

int nbits_of_int;
int log_nbits_of_int;

unsigned int *bitvec;

/*
   For recording the norms of res_k for the scalar and block columns
   (Just for some data analysis.)
*/
FILE *resplot_fptr=NULL;
int start_col;
double scalar_resnorm;
char *block_flag,*scalar_flag;
int num_bad_cols;

int *len_all=NULL;
int *rlen_all=NULL;
int *slen_all=NULL;

int *ptr_offsets=NULL;
int *rptr_offsets=NULL;
int *A_offsets=NULL;

/**********************************************************************/

int bspai
(matrix *A, matrix **bspai_matrix,
 FILE *messages_arg,   /* file for warning messages */
 double epsilon_arg,   /* tolerance */
 int nbsteps_arg,      /* max number of accepted elements per line
                          (sc = 0); is set to ld+lu+1 for spar = 2 */
 int max_arg,          /* max dimensions of I, q, etc. */
 int maxnew_arg,       /* max number of new entries per step */
 int bs,               /* block size */
 int cache_size_arg,   /* one of (1,2,3,4,5,6) indicting size of cache */
                       /* cache_size == 0 indicates no caching */
 int verbose_arg,
 int spar_arg,
 int lower_diag_arg,
 int upper_diag_arg,
 double tau_arg)
{
  matrix *B;
  matrix *M_B;
  matrix *M;
  extern int A_max_block_size;
  int total_bad_cols;
  int j,nnz,ierr;
  if (verbose_arg) matrix_statistics(A,"A");
  if (spar_arg && A->max_block_size > 1)	/* MH: ??? */
    {
	printf("Sorry, sc > 0 only supported for block size 1\n");
	exit(1);
    }
  /* Scalar case */
  if (bs == 1) {
    if (verbose_arg && (A->myid == 0))
      printf("\nConstructing scalar SPAI matrix M from A\n");

    if ((ierr = spai(A, &M,
	     messages_arg,
	     epsilon_arg,
	     nbsteps_arg,
	     max_arg,
	     maxnew_arg,
	     cache_size_arg,
	     verbose_arg,
             spar_arg,
             lower_diag_arg,
             upper_diag_arg,
		     tau_arg)) != 0)  return ierr;
#ifdef MPI
    MPI_Barrier(A->comm);
#endif

  }

  /* Block case */
  else {

    /* Convert A to either a constant or variable block matrix B. */
    /* The 3rd argument can be used to set an upper limit
       on the block size. For now it's unused (i.e., zero). */
    if (verbose_arg && (A->myid == 0))
      printf("\nConstructing block matrix B from A\n");
    B = block_matrix(A,bs,0,verbose_arg);

    if (verbose_arg) matrix_statistics(B,"B");

    if (verbose_arg && (A->myid == 0))
      printf("Constructing block SPAI matrix M_B from B\n");
    if ((ierr = spai(B, &M_B,
	       messages_arg,
	       epsilon_arg,
	       nbsteps_arg,
	       max_arg,
	       maxnew_arg,
	       cache_size_arg,
	       verbose_arg,
               spar_arg,
               lower_diag_arg,
               upper_diag_arg,
		     tau_arg)) != 0)  return ierr;

#ifdef MPI
    MPI_Barrier(A->comm);
#endif

    /* VE */
    if (num_bad_cols) printf("SPAI: # bad columns = %d\n",num_bad_cols);

    if (verbose_arg) matrix_statistics(M_B,"M_B");

    /*  Convert M_B to scalar matrix.  */
    if (verbose_arg && (A->myid == 0))
      printf("\nConstructing scalar SPAI matrix M from M_B\n");
    M = scalar_matrix(M_B,verbose_arg);

    sp_free_matrix(B);
    sp_free_matrix(M_B);

  }

  if (verbose_arg) matrix_statistics(M,"M");

  *bspai_matrix = M;
  return 0;

}

/**********************************************************************/

int spai
(matrix *A, matrix **spai_mat,
 FILE *messages_arg,   /* file for warning messages */
 double epsilon_arg,   /* tolerance */
 int nbsteps_arg,      /* max number of "improvement" steps per line */
 int max_arg,          /* max dimensions of I, q, etc. */
 int maxnew_arg,       /* max number of new entries per step */
 int cache_size_arg,   /* one of (1,2,3,4,5,6) indicting size of cache */
                       /* cache_size == 0 indicates no caching */
 int verbose_arg,
 int spar_arg,
 int lower_diag_arg,
 int upper_diag_arg,
 double tau_arg)
{
  matrix *M;
  int col,ierr;

  int cache_sizes[6];

  /* Only create resplot for the numprocs=1 case. */
  if (debug && (A->numprocs == 1)) {
    resplot_fptr = fopen("resplot","w");
    fprintf(resplot_fptr,
	    "ep=%5.5lf ns=%d mn=%d bs=%d\n",
	    epsilon_arg,nbsteps_arg,maxnew_arg,A->bs);
    fprintf(resplot_fptr,"\n");
    fprintf(resplot_fptr,"scol: scalar column number\n");
    fprintf(resplot_fptr,"srn:  scalar resnorm\n");
    fprintf(resplot_fptr,"bcol: block column number\n");
    fprintf(resplot_fptr,"brn:  block resnorm\n");
    fprintf(resplot_fptr,"* indicates epsilon not attained\n");
    fprintf(resplot_fptr,"\n");
    fprintf(resplot_fptr,"   scol   srn       bcol   brn\n");
  }


  start_col = 0;
  num_bad_cols = 0;

  cache_sizes[0] = 101;
  cache_sizes[1] = 503;
  cache_sizes[2] = 2503;
  cache_sizes[3] = 12503;
  cache_sizes[4] = 62501;
  cache_sizes[5] = 104743;



  if (verbose_arg && !A->myid) {
    if (spar_arg == 0)
       printf("\n\nComputing SPAI: epsilon = %f\n",epsilon_arg);
    else if (spar_arg == 1)
       printf("\n\nComputing SPAI: tau = %f\n",tau_arg);
    else if (spar_arg == 2)
       printf("\n\nComputing SPAI: # diagonals = %d\n",
       lower_diag_arg+upper_diag_arg+1);
    fflush(stdout);
  }

  epsilon = epsilon_arg;
  message = messages_arg;
  maxnew = maxnew_arg;
  max_dim = max_arg;

  /* Determine maximum number of scalar nonzeros
     for any column of M */
  if (spar_arg == 0)
    {
     nbsteps = nbsteps_arg;
     maxapi = A->max_block_size * (1 + maxnew*nbsteps);
    }
  else if(spar_arg == 1)
    {
     nbsteps = A->maxnz;
     maxapi = A->max_block_size * (1 + nbsteps);
    }
  else
    {
     nbsteps = lower_diag_arg+upper_diag_arg+1;
     maxapi = A->max_block_size * (1 + nbsteps);
    }
  allocate_globals(A);

#ifdef MPI
  MPI_Barrier(A->comm);
#endif

  if ((cache_size_arg < 0) || (cache_size_arg > 6)) {
    fprintf(stderr,"illegal cache size in spai\n");
    exit(1);
  }

  if (cache_size_arg > 0)
    ht = init_hash_table(cache_sizes[cache_size_arg-1]);

  M = clone_matrix(A);

  ndone = 0;
  Im_done = 0;
  all_done = 0;
  next_line = 0;

  /* Timing of SPAI starts here.
     In a "real production" code everything before this could be static.
  */
  if (verbose_arg) start_timer(ident_spai);

  if ((ierr = precompute_column_square_inverses(A)) != 0)  return ierr;

#ifdef MPI
  MPI_Barrier(A->comm);
#endif

  for (;;) {

    col = grab_Mline(A, M, A->comm);

    if (debug && col >= 0) {
      fprintf(fptr_dbg,"col=%d of %d\n",col,A->n);
      fflush(fptr_dbg);
    }

    if (col < 0 ) break;
    if ((ierr =
         spai_line(A,col,spar_arg,lower_diag_arg,upper_diag_arg,tau_arg,M)) != 0)  return ierr;

  }

#ifdef MPI

  say_Im_done(A,M);

  do {
    com_server(A,M);
  }
  while (! all_done);
  MPI_Barrier(A->comm);

#endif

#ifdef MPI
  MPI_Barrier(A->comm);
#endif

  if (verbose_arg) {
    stop_timer(ident_spai);
    report_times(ident_spai,"spai",0,A->comm);
  }

  free_globals(nbsteps);
  free_hash_table(ht);

  if (resplot_fptr) fclose(resplot_fptr);

  *spai_mat = M;
  return 0;

}

/**********************************************************************/

int spai_line
(matrix *A,
 int col,
 int spar,
 int lower_diag,
 int upper_diag,
 double tau,
 matrix *M)
{
  int s,nbq,nnz,dimr,block_width;
  double scalar_resnorm,block_resnorm,adjust_epsilon;

  int i,index,pe,len,ierr;
  int row_address;

  int *rptr;
  double *aptr;
  int j, k, ptr, low_c, up_c, ccol, row;
  int rlen;
  int *buf;
  int *rbuf;
  double *vbuf;
  double comp_max, tau_limit = 1 - tau;

  block_width = A->block_sizes[col];
  adjust_epsilon = epsilon*sqrt((double) block_width);

  if (spar == 1)   /* mark elements depending on tau parameter */
   {
    comp_max = 0;
/* find maximum in column resp. row if transposed */
    for (j=0; j<A->lines->len[col]; j++)
      {
       ptr = A->lines->ptrs[col][j];
       if (comp_max < fabs( A->lines->A[col][j]))
           comp_max = fabs( A->lines->A[col][j]);
      }
/* keep diagonal and elements about fraction of maximum */
    for (i=0, j=0; j<A->lines->len[col]; j++)
      {
       ptr = A->lines->ptrs[col][j];
       if (ptr == col + A->my_start_index
           || fabs(A->lines->A[col][j]/comp_max) > tau_limit)
       {
         n1[i] = A->block_sizes[j];
	 J->ptr[i++] = ptr;
	}
      }
     J->len = i;
     J->slen = i;
     dimr = nnz = 0;
    }
  else if (spar == 2)   /* set diagonals - mind switching cols and rows */
    {
     if ((low_c = col-upper_diag) < 0) low_c = 0;
     if ((up_c = col+lower_diag) > A->n-1) up_c = A->n-1;
     for (i=0, j=low_c; j<=up_c; j++,i++)
       {
        J->ptr[i] = j;
        n1[i] = A->block_sizes[j];
       }
     J->len = i;
     J->slen = i;
     dimr = nnz = 0;
    }
  else /* initial sparsity diagonal */
    {
     J->ptr[0] = col;
     J->len = 1;
     J->slen = block_width;
     n1[0] = block_width;
     dimr = nnz = 0;
    }
  /* compute I */
  getrows(A,M,J,I);

  copyvv(J,J_tilde);

  for (s=0,
	 nbq = 0,
	 TAU_ptr[0] = 0,
                            /* effectively infinity */
	 scalar_resnorm=block_resnorm=1000000*epsilon;
       (s < nbsteps);
       s++,
	 nbq++) {

    com_server(A,M);

    full_matrix(A,M,max_dim, Ahat);

    n2[s] = I->slen - dimr;

    /* compute solution -> x, residual, and update QR */
    if ((ierr = qr(A,col,nbq,dimr)) != 0)  return ierr;

    nnz = J->len;
    dimr = J->slen;

    /* is solution good enough? */
    /* Use Froebenius norm */
    convert_to_block
      (res,resb,col,I->ptr,A,max_dim,I->len);
    block_resnorm = frobenius_norm(resb,block_width,I->slen);

    if (debug) {
      fprintf(fptr_dbg,"  s=%d col=%d of %d block_resnorm=%12.4le\n",
	      s,col,A->n,block_resnorm);
      fflush(fptr_dbg);
    }
    if (spar == 1         /* row population with tau parameter */
     || spar == 2) break; /* fixed diagonals - no further ado */
    if (block_resnorm <= adjust_epsilon)  break;

    /* Don't bother with last augment_sparsity */
    if (s == (nbsteps-1)) break;

    if (! augment_sparsity(A,M,col,maxapi,block_resnorm)) break;

    getrows(A,M, J_tilde,I_tilde);

    deleter(I,I_tilde,A);
    if (! append(J,J_tilde)) break;   /* J <- J U J_tilde */
    if (! append(I,I_tilde)) break;   /* I <- I U I_tilde */

  }

  if (block_resnorm > adjust_epsilon && spar == 0) {
    num_bad_cols++;
    if (message) {
      fprintf(message,
	      "could not meet tol, col=%d resnorm = %le, adjust_epsilon = %le\n",
	      col+1,
	      block_resnorm/sqrt((double) block_width),
	      adjust_epsilon);
      fflush(message);
    }
  }

  if (resplot_fptr) {
    for (i=0; i<block_width; i++) {
      if (block_resnorm <= adjust_epsilon) block_flag = " ";
      else block_flag = "*";
      scalar_resnorm = frobenius_norm(&res[i*max_dim],1,I->slen);
      if (scalar_resnorm <= epsilon) scalar_flag = " ";
      else scalar_flag = "*";
      fprintf(resplot_fptr,"%6d   %5.3lf %s %6d   %5.3lf %s\n",
	      start_col+i,
	      scalar_resnorm,
	      scalar_flag,
	      col,
	      block_resnorm/sqrt((double) block_width),
	      block_flag);
    }
    start_col += block_width;
  }

  /* current solution in x, up to nnz, written to M(k,:) */
  /* convert x to block structure */
  convert_to_block
    (x,xb,col,J->ptr,A,max_dim,nnz);

  put_Mline(A,M, col, J->ptr, xb, nnz, J->slen);

  for (i=0; i<nbsteps; i++) {
    if (Qlist[i]) {
      free(Qlist[i]);
      Qlist[i] = NULL;
    }
    else break;
  }
  return 0;
}

/**********************************************************************/

void allocate_globals
(matrix *A)
{
  int n,i,bs,bs2;
  int *my_offsets;
  int num;
  int global_index;

  TAU_ub = R_ub = Z_ub = 0;

  n = A->n;
  bs = A->max_block_size;
  bs2 = bs*bs;

  J = new_index_set(NULL,max_dim,"J");
  I = new_index_set(NULL,max_dim,"I");
  J_tilde = new_index_set(NULL,max_dim,"J_tilde");
  I_tilde = new_index_set(NULL,max_dim,"I_tilde");

  remote_buf = new_int_array(NULL,A->maxnz,"remote_buf");
  remote_rbuf   = new_int_array(NULL,A->maxnz,"remote_rbuf");
  remote_abuf   = new_double_array(NULL,bs2*(A->maxnz+1),"remote_abuf");

  is = new_index_set(NULL,max_dim,"is");

  isres = new_index_set(NULL,max_dim,"isres");
  res = new_double_array(NULL,bs2*max_dim,"res");
  resb = new_double_array(NULL,bs2*max_dim,"resb");

  n1     = new_int_array(NULL,nbsteps+1,"n1");
  n2     = new_int_array(NULL,nbsteps,"n2");

  Ahat_size = max_dim*maxapi;
  Ahat   = new_double_array(NULL,Ahat_size,"Ahat");
  if (debug) init_double_array(Ahat,Ahat_size,init_val);

  TAU_ptr = (int *) new_int_array(NULL,nbsteps+1,"TAU_ptr");

  TAU_size = 1000;
  TAU    = new_double_array(NULL,TAU_size,"TAU");
  if (debug) init_double_array(TAU,TAU_size,init_val);

  R_size = 1000;
  R      = new_double_array(NULL,R_size,"R");
  if (debug) init_double_array(R,R_size,init_val);

  Z_size = 1000;
  Z      = new_double_array(NULL,Z_size,"Z");
  if (debug) init_double_array(Z,Z_size,init_val);

  rw = new_double_array(NULL,bs2*max_dim,"rw");
  rs     = new_double_array(NULL,bs2*max_dim,"rs");
  rz     = new_double_array(NULL,bs2*max_dim,"rz");
  x      = new_double_array(NULL,bs2*max_dim,"x");
  xb     = new_double_array(NULL,bs2*max_dim,"xb");

  Qlist = (double **) mmalloc(nbsteps*sizeof(double *));
  for (i=0; i<nbsteps; i++) Qlist[i] = NULL;

  temp_block = new_double_array(NULL,bs2,"temp_block");
  minus_Bj = new_double_array(NULL,bs2,"minus_Bj");
  sum_block = new_double_array(NULL,bs2,"sum_block");

  nbits_of_int = sizeof(int)*8;
  if (nbits_of_int == 32) {
    log_nbits_of_int = 5;}
  else if (nbits_of_int == 64) {
    log_nbits_of_int = 6;}
  else {
    if (A->myid == 0)
      fprintf(stderr,"unsupported word size: %d\n",nbits_of_int);
    exit(1);
  }

  bitvec = (unsigned int *) mmalloc(sizeof(unsigned int)*A->n);
  if (! bitvec) {
      fprintf(stderr,"failed to malloc bitvec\n");
      exit(1);
  }

  len_all = (int *) mmalloc(A->n*sizeof(int));
  rlen_all = (int *) mmalloc(A->n*sizeof(int));
  slen_all = (int *) mmalloc(A->n*sizeof(int));

#ifdef MPI

  ptr_offsets = (int *) mmalloc(A->n*sizeof(int *));
  rptr_offsets = (int *) mmalloc(A->n*sizeof(int *));
  A_offsets = (int *) mmalloc(A->n*sizeof(int *));

  MPI_Barrier(A->comm);

  MPI_Allgatherv
    ((void *) A->lines->len, A->mnl, MPI_INT,
     (void *) len_all, A->mnls, A->start_indices, MPI_INT,
     A->comm);

  MPI_Barrier(A->comm);

  MPI_Allgatherv
    ((void *) A->lines->rlen, A->mnl, MPI_INT,
     (void *) rlen_all, A->mnls, A->start_indices, MPI_INT,
     A->comm);

  MPI_Barrier(A->comm);

  MPI_Allgatherv
    ((void *) A->lines->slen, A->mnl, MPI_INT,
     (void *) slen_all, A->mnls, A->start_indices, MPI_INT,
     A->comm);

  /* offsets for ptrs */
  my_offsets = (int *) mmalloc(A->mnl*sizeof(int *));
  my_offsets[0] = 0;
  for (i=1; i<A->mnl; i++) {
    my_offsets[i] = my_offsets[i-1] + A->lines->len[i-1];
  }

  MPI_Barrier(A->comm);

  MPI_Allgatherv
    ((void *) my_offsets, A->mnl, MPI_INT,
     (void *) ptr_offsets, A->mnls, A->start_indices, MPI_INT,
     A->comm);

  /* offsets for rptrs */
  my_offsets[0] = 0;
  for (i=1; i<A->mnl; i++) {
    my_offsets[i] = my_offsets[i-1] + A->lines->rlen[i-1];
  }

  MPI_Barrier(A->comm);

  MPI_Allgatherv
    ((void *) my_offsets, A->mnl, MPI_INT,
     (void *) rptr_offsets, A->mnls, A->start_indices, MPI_INT,
     A->comm);

  /* offsets for coefficients */
  my_offsets[0] = 0;
  global_index = A->my_start_index;
  for (i=1; i<A->mnl; i++) {
    bs = A->block_sizes[global_index];
    num = bs*A->lines->slen[i-1];
    my_offsets[i] = my_offsets[i-1] + num;
    global_index++;
  }

  MPI_Barrier(A->comm);

  MPI_Allgatherv
    ((void *) my_offsets, A->mnl, MPI_INT,
     (void *) A_offsets, A->mnls, A->start_indices, MPI_INT,
     A->comm);

  free(my_offsets);

#endif

}

/**********************************************************************/

void init_double_array(double *a, int size, double val)
{
  int i;

  for (i=0; i<size; i++) a[i] = val;
}

/**********************************************************************/

int amount_touched(double *a, int size, double val)
{
  int i,last;

  last = size-1;
  for (i=last; i>=0; i--) {
    if (a[i] == val) last = i;
  }

  return(last);
}

/**********************************************************************/

void free_globals()
{
  int i;

  free_index_set(J);
  free_index_set(I);
  free_index_set(J_tilde);
  free_index_set(I_tilde);

  free(n1);
  free(n2);

  free(TAU_ptr);

  free(Ahat);
  free(R);
  free(TAU);
  free(Z);

  free(remote_buf);
  free(remote_rbuf);
  free(remote_abuf);

  free_index_set(is);
  free_index_set(isres);
  free(res);
  free(resb);

  free(rw);
  free(rs);
  free(rz);
  free(x);
  free(xb);

  for (i=0; i<nbsteps; i++)
    if (Qlist[i]) free(Qlist[i]);
  free(Qlist);

  free(temp_block);
  free(minus_Bj);
  free(sum_block);

  if (si_buffer) free(si_buffer);
  if (si_indices) free(si_indices);

  if (len_all) free(len_all);
  if (rlen_all) free(rlen_all);
  if (slen_all) free(slen_all);

  if (ptr_offsets) free(ptr_offsets);
  if (rptr_offsets) free(rptr_offsets);
  if (A_offsets) free(A_offsets);

  free(bitvec);

}

/**********************************************************************/

int *new_int_array(int *ptr, int n, char *s)
{

  if (debug) {
    fprintf(fptr_dbg,"new_int_array allocating %d for %s\n",
	    n,s);
  }
  ptr = (int *) realloc(ptr,sizeof(int)*n);
  if (! ptr) {
    printf("\n realloc failed in new_int_array\n");
    fflush(stdout);
    exit(1);
  }
  return(ptr);
}

/**********************************************************************/

double *new_double_array(double *ptr, int n,char *s)
{

  if (debug) {
    fprintf(fptr_dbg,"new_double_array allocating %d for %s\n",
	    n,s);
  }
  ptr = (double *) realloc(ptr,sizeof(double)*n);
  if (! ptr) {
    printf("\n realloc failed in new_double_array\n");
    fflush(stdout);
    exit(1);
  }
  return(ptr);
}
