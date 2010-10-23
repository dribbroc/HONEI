/*
   SPAI MPI/C version @Copyright 1996,  All Rights Reserved
   Stephen Barnard
*/

/* ********************************************************************
   ********************************************************************

      Main program illustrating the use of the SPAI preconditioner

      See README for instructions.

      This program will:
         - read a matrix A and right-hand-side b
         - compute the SPAI preconditioner M
         - solve the system Ax = b with BICGSTAB

   ********************************************************************
   ********************************************************************/

#include <string.h>
#include <stdio.h>
#include "basics.h"
#include "command_line.h"
#include "vector.h"
#include "index_set.h"
#include "matrix.h"
#include "bicgstab.h"
#include "spai.h"
#include "read_mm_matrix.h"
#include "variable_blocks.h"
#include "timing.h"
#include "debug.h"

int **alloc_list;
int alloc_size;
int alloc_curr;
int act_iter;
int precond_time = 0, solver_time = 0;

void testspai_start();
void time_stamp(char*);

time_t start_time, last_time;
struct tm* tm;

int main(int argc, char *argv[])
{
  /* SPAI parameters */
  double epsilon_param;
  int    nbsteps_param;
  int    max_param;
  int    maxnew_param;
  int    cache_size_param;
  int    block_size_param;
  int    symmetric_pattern_param;
  int    left_precon_param;
  int    verbose_param;
  int    spar_param;
  int    lower_diag_param;
  int    upper_diag_param;
  double tau_param;


  /* BICGSTAB parameters */
  int    max_iter_param;
  double tolerance_param;
  int    output  = 2;
  int j;

  /* Binary file ? */
  int    binary_param;

  matrix *A = NULL;
  matrix *M = NULL;
  vector *x = NULL;
  vector *rhs = NULL;

  char *matrix_file;
  char *rhs_file=NULL;
  char *message_file;
  char message_file_name[1024];
  FILE *fptr_message;
  FILE *fptr_rhs;
  FILE *fptr_x;
  FILE *fptr_summary;

  int numprocs,myid,ierr,title=0;
#ifdef MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Barrier(MPI_COMM_WORLD);
#else
#define MPI_COMM_WORLD NULL
  numprocs = 1;
  myid = 0;
#endif

  init_SPAI();
  testspai_start();

  /* The command-line parameters */
  parameters(argc,
             argv,
             &matrix_file,
             &rhs_file,
             &epsilon_param,
             &nbsteps_param,
             &max_param,
             &maxnew_param,
             &cache_size_param,
             &message_file,
             &max_iter_param,
             &tolerance_param,
             &block_size_param,
             &symmetric_pattern_param,
             &left_precon_param,
             &binary_param,
             &verbose_param,
             &debug,
             &spar_param,
             &lower_diag_param,
             &upper_diag_param,
             &tau_param,
             MPI_COMM_WORLD);

  if (debug) {
    sprintf(dbg_file_name,"dbg%d",myid);
    fptr_dbg = fopen(dbg_file_name,"w");
  }

  if (! strlen(message_file)) {
    fptr_message = NULL;
  }
  else {
    if (numprocs > 1) sprintf(message_file_name,"%s%d",message_file,myid);
    else sprintf(message_file_name,"%s",message_file);
    fptr_message = fopen(message_file_name,"w");
  }

  /********************************************************************/
  /********************************************************************/
  /*             Read the matrix A (scalar matrix)                    */
  /********************************************************************/
  /********************************************************************/

  A = read_mm_matrix(matrix_file,
                     1,
                     1,
                     symmetric_pattern_param,
                     left_precon_param,
                     binary_param,
                     verbose_param,
                     MPI_COMM_WORLD);

  /********************************************************************/
  /********************************************************************/
  /*           Read the right-hand side, if any                       */
  /********************************************************************/
  /********************************************************************/

  x = uniform_vector(A->n,A->mnls[myid],1.0);
  if (! rhs_file) {
    if (myid == 0 && verbose_param) {
      printf("No RHS file ... Using a rhs of A times all ones.\n");
    }
    rhs = new_vector(A->n,A->mnls[myid]);
    if (! left_precon_param)
      A_times_v_cc(A,x,rhs);
    else
      A_times_v_rc(A,x,rhs);
  }
  else {
    rhs = read_rhs_for_matrix(rhs_file, A);
  }

  if (verbose_param > 2) time_stamp("matrix read");

  if (verbose_param > 1) write_matrix_mm(A,"A.mm",left_precon_param);

  /********************************************************************/
  /********************************************************************/
  /*              Compute the SPAI preconditioner                     */
  /********************************************************************/
  /********************************************************************/

  if ((ierr = bspai
    (A, &M,
     fptr_message,
     epsilon_param,
     nbsteps_param,
     max_param,
     maxnew_param,
     block_size_param,
     cache_size_param,
     verbose_param,
     spar_param,
     lower_diag_param,
     upper_diag_param,
     tau_param)) != 0) exit(ierr);

  precond_time = last_time;
  if (verbose_param > 2) time_stamp("precon done");
  precond_time = (int)last_time - precond_time;

  if (verbose_param > 1) write_matrix_mm(M,"M.mm",left_precon_param);
  if (verbose_param > 2) time_stamp("matrix out");

  /********************************************************************/
  /********************************************************************/
  /*                         Run bi-cgstab.                           */
  /********************************************************************/
  /********************************************************************/

  /* Solution vector */
  rzeros(x);

  /* Solve for Ax=b (where b is the rhs) using M */

  if (! left_precon_param) {
    /* Right preconditioner */
    act_iter = bicgstab_R
      (A_times_v_cc,
       A,
       M,
       x,
       rhs,
       max_iter_param,
       tolerance_param,
       verbose_param);
  }

  else {
    /* Left preconditioner */
    act_iter = bicgstab_L
      (A_times_v_rc,
       A,
       M,
       x,
       rhs,
       max_iter_param,
       tolerance_param,
       verbose_param);
  }
  solver_time = last_time;
 if (verbose_param > 2)  time_stamp("after bicgstab");
  solver_time = (int)last_time - solver_time;

  write_vector_mm(x,"solution.mm",MPI_COMM_WORLD);
  if (verbose_param > 2) time_stamp("end");
 if (verbose_param > 3)
   {
       if ((fptr_summary = fopen("spai_summary.file", "r")))
        fclose(fptr_summary);
    else title = 1;
       if ((fptr_summary = fopen("spai_summary.file", "a")))
      {
        if (title)
          fprintf(fptr_summary,
                   "# dd/mm/yy hh.mm.ss size A->tot_nnz sp_param ld ud mn M->tot_nnz mi precond_time solver_time iter eps tau tol\n");
          fprintf(fptr_summary,
 "%02d/%02d/%02d %02d.%02d.%02d %d %d %d %d %d %d %d %d %d %d %d %12.4e %f %12.4e\n",
          tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
          tm->tm_hour, tm->tm_min, tm->tm_sec, A->n, A->tot_nnz, spar_param,
          lower_diag_param, upper_diag_param, maxnew_param, M->tot_nnz,
          max_iter_param, precond_time, solver_time,
          act_iter, epsilon_param, tau_param, tolerance_param);
      }
   }
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  sp_free_matrix(A);
  sp_free_matrix(M);

  free_vector(x);
  free_vector(rhs);

  if(fptr_message) fclose(fptr_message);

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
}


void testspai_start()
{

  time(&start_time); /* initialize timing */
  tm = localtime(&start_time); /* split system time */
  last_time = start_time;
  printf("\n  ++++++++++++++++++++++++++++++++\n");
  printf("  + %s %02d/%02d/%02d %02d.%02d.%02d +\n", "SPAI v.3.2",
         tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
         tm->tm_hour, tm->tm_min, tm->tm_sec);
  printf("  ++++++++++++++++++++++++++++++++\n\n");
}


void time_stamp(char* place)
{
  time_t now;
  int k, l;
  time(&now);    /* get system time */
  k = (int)now - (int)start_time;
  l = (int)now - (int)last_time;
  last_time = now;
  printf("place: %-16s sec.s since start: %d   since last call: %d\n",
         place, k, l);
}
