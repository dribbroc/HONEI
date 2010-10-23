/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "timing.h"
#include "debug.h"

double timer0[NUMTIMERS],timer1[NUMTIMERS];
double max_time[NUMTIMERS],sum_time[NUMTIMERS];

/**********************************************************************/

void init_timers()
{
  int i;

#ifdef MPI
  for (i=0; i<NUMTIMERS; i++) {
    max_time[i] = sum_time[i] = 0.0;
  }
#endif

}

/**********************************************************************/

void start_timer(int itime)
{

#ifdef MPI
  timer0[itime] = MPI_Wtime();
#endif


}

/**********************************************************************/

void stop_timer(int itime)
{
  double delta;

#ifdef MPI
  timer1[itime] = MPI_Wtime();
  delta = timer1[itime]-timer0[itime];
  if (delta > max_time[itime]) max_time[itime] = delta;
  sum_time[itime] += delta;
#endif

}

/**********************************************************************/

void report_times
(int itime,
 char *ident,
 int all_processors,
 SPAI_Comm comm)
{

  double max,*sum;
  int i;

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
  MPI_Barrier(comm);
  MPI_Reduce(&sum_time[itime],&max,1,
	     MPI_DOUBLE,MPI_MAX,0,comm);

  if (myid == 0) {
    printf("\n");
    printf("======= timing for %s: %14.7lf seconds =======\n",ident,max);
    printf("\n");
  }

  if (all_processors) {
    sum = (double *) mmalloc(numprocs*sizeof(double));
    MPI_Gather(&sum_time[itime],1,MPI_DOUBLE,
	       sum,1,MPI_DOUBLE,
	       0,comm);
    if (myid == 0) {
      for (i=0; i<numprocs; i++) {
	printf("  processor %d:  sum: %14.7lf\n",
	       i,sum[i]);
      }
      printf("\n");
    }
    free(sum);
  }

  if (myid == 0) {
    printf("\n");
  }

#endif

}
