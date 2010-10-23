/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "vector.h"

/**********************************************************************/
/* constructs a new vector of size n with same distribution as M */
/* If initial is null, initialize to zeros */

vector *new_vector
(int n,
 int mnl)
{
  vector *vec;
  int i,max_mnl;

  vec = (vector *) mmalloc(sizeof(vector));
  if (! vec) {
    printf("\n failed to allocate vector");
    exit(1);
  }

  vec->n = n;
  vec->mnl = mnl;

  vec->v = (double *) mmalloc(sizeof(double)*vec->mnl);
  if (! vec->v) {
    printf("\n failed to allocate vector buffer");
    exit(1);
  }

  rzeros(vec);

  return(vec);
  }

/**********************************************************************/

void free_vector(vector *vec)
{

  if (vec->v) free(vec->v);
  free(vec);
}

/**********************************************************************/

void rzeros(vector *v)
{
  int i;

  for (i=0; i<v->mnl; i++) {
    v->v[i]=0.0;
  }
}

/**********************************************************************/
/* computes v + c*w -> z
   c is a constant */

void v_plus_cw(vector *v,
	       vector *w,
	       double c,
	       vector *z)
{
  int i;

  for (i=0; i<v->mnl; i++)
    z->v[i] = v->v[i] + c*w->v[i];
}

/**********************************************************************/
/* copies v1 to v2 */

void rcopy_vv(vector *v1,
	      vector *v2)
{
  int i;

  for (i=0; i<v1->mnl; i++)
    v2->v[i] = v1->v[i];
}

/**********************************************************************/
/* computes 2 norm of v */

double norm(vector *v, SPAI_Comm comm)
{
  int i;
  double r,rsum;

  r = 0.0;
  for (i=0; i<v->mnl; i++)
    r += v->v[i]*v->v[i];

#ifdef MPI

  MPI_Allreduce
    (&r, &rsum, 1, MPI_DOUBLE, MPI_SUM, comm);

#else

  rsum = r;

#endif

  r = sqrt(rsum);
  return(r);
}

/**********************************************************************/
/* computes v^T * w */

double dot(vector *v, vector *w, SPAI_Comm comm)
{
  int i,num;
  double r,rsum;

  r=0.0;
  for (i=0; i<v->mnl; i++)
    r += v->v[i]*w->v[i];

#ifdef MPI

  MPI_Allreduce
    (&r, &rsum, 1, MPI_DOUBLE, MPI_SUM, comm);

#else

  rsum = r;

#endif

  return(rsum);
}

/**********************************************************************/
/* serial only -- for debugging */

void write_vector(FILE *fptr, vector *v)
{
  int i;

  for (i=0; i<v->mnl; i++) {
    fprintf(fptr,"i=%d v=%le\n",i,v->v[i]);
  }
}

/**********************************************************************/

vector *uniform_vector(int n, int mnl, double val)
{
  vector *v;
  int i;

  v = new_vector(n,mnl);
  for (i=0; i<mnl; i++) {
    v->v[i] = val;
  }

  return(v);
}


/**********************************************************************/
/* writes in Matrix Market array format */

#define MatrixMarketBanner "%%MatrixMarket"

void write_vector_mm(vector *v, char *filename, SPAI_Comm comm)
{
  int j;
  FILE *fptr;
  char fullfilename[1024];
  char cat_cmd[1024];
  char rm_cmd[1024];

  int numprocs,myid;
#ifdef MPI
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myid);
  MPI_Barrier(comm);
#else
  numprocs = 1;
  myid = 0;
#endif

  if (numprocs > 1) {
    sprintf(fullfilename,"%s_tmp%5.5d",filename,myid);
    sprintf(cat_cmd,"cat %s_tmp* > %s",filename,filename);
    sprintf(rm_cmd,"rm -f %s_tmp*",filename);
  }
  else
    sprintf(fullfilename,"%s",filename);

  fptr = fopen(fullfilename,"w");

  /* write Matrix-Market header */
  if (myid == 0) {
    fprintf(fptr, "%s ", MatrixMarketBanner);
    fprintf(fptr,"matrix array real general\n");
    fprintf(fptr,"%d %d\n",v->n,1);
    fflush(fptr);
  }

  for (j=0; j<v->mnl; j++)
    fprintf(fptr,"%le\n",v->v[j]);

  fflush(fptr);
  fclose(fptr);

#ifdef MPI

  MPI_Barrier(comm);

#endif

  if (numprocs > 1) {
    if (myid == 0) {
      system(cat_cmd);
    }
  }

#ifdef MPI

  MPI_Barrier(comm);

#endif

  if (numprocs > 1) {
    if (myid == 0) {
      system(rm_cmd);
    }
  }

}

