/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#include "bicgstab.h"

/**********************************************************************/
/* Solves A*x = b
   M is a right approximate inverse of A
   x is an initial guess
   right preconditioner
*/

int bicgstab_R
(void Av_proc(),
 matrix *A,
 matrix *M,
 vector *x,
 vector *b,
 int maxiter,
 double tol,
 int verbose)
{
  vector *r,*rhat,*p,*v,*w,*z;
  int k=0;
  double c=-1,rho,alpha,beta=0,omega=0,tmp;
  double res,res0;

  if (verbose) start_timer(ident_bicgstab);

  if ((A->myid == 0) && verbose) {
    printf("\n\n      bicgstab_R: tol=%e maxiter=%d\n",
	   tol,maxiter);
    printf("iteration     residual\n");
  }

  p = new_vector(A->n,A->mnl);
  v = new_vector(A->n,A->mnl);
  w = new_vector(A->n,A->mnl);
  z = new_vector(A->n,A->mnl);
  r = new_vector(A->n,A->mnl);
  rhat = new_vector(A->n,A->mnl);

  /* compute initial residual */
  Av_proc(M,x,w);
  Av_proc(A,w,z);
  v_plus_cw(b,z,c,r);
  rcopy_vv(r,rhat);
  res0 = norm(b,A->comm);

  do {

    k++;

    /* compute new p */
    v_plus_cw(p,v,-omega,z);
    v_plus_cw(r,z,beta,p);

    /* compute new v, r, and alpha */
    Av_proc(M,p,w);
    Av_proc(A,w,v);

    rho = dot(rhat,r,A->comm);
    alpha = rho/dot(rhat,v,A->comm);
    v_plus_cw(r,v,-alpha,w);

    /* compute new omega */
    Av_proc(M,w,r);
    Av_proc(A,r,z);

    if ((tmp = dot(z,z,A->comm)) < 1.e-20)  omega = 1;
    else omega = dot(z,w,A->comm)/tmp;

    /* compute new x and new r */
    v_plus_cw(x,p,alpha,r);
    v_plus_cw(r,w,omega,x);
    v_plus_cw(w,z,-omega,r);
    beta = (alpha/omega)*(dot(rhat,r,A->comm)/rho);

    /* compute exact residual -> w */
    Av_proc(M,x,w);
    Av_proc(A,w,z);
    v_plus_cw(b,z,c,w);
    res = norm(w,A->comm);

    if (!A->myid && verbose && (k%10 == 0 || verbose > 0)) {
      printf("%9d %.6e\n", k, res/res0);
      fflush(stdout);
    }
  }
  while ((res/res0  > tol) && (k < maxiter));

  if (!A->myid && verbose) printf("\n%9d %.6e\n", k,res/res0);

  /* compute solution M*x */
  Av_proc(M,x,w);
  rcopy_vv(w,x);

  free_vector(p);
  free_vector(v);
  free_vector(w);
  free_vector(z);
  free_vector(r);
  free_vector(rhat);

  if (verbose) {
    stop_timer(ident_bicgstab);
    report_times(ident_bicgstab,"bicgstab",0,A->comm);
  }
  return k;

}

/**********************************************************************/
/* Solves A*x = b
   M is a left preconditioner
   x is an initial guess
*/

int bicgstab_L
(void Av_proc(),
 matrix *A,
 matrix *M,
 vector *x,
 vector *b,
 int maxiter,
 double tol,
 int verbose)
{
  vector *bb,*r,*rhat,*p,*v,*w,*z;
  int k=0;
  double c=-1,rho,alpha,beta=0,omega=0,tmp;
  double res,res0;

  if (verbose) start_timer(ident_bicgstab);

  if ((A->myid == 0) && verbose) {
    printf("\n\n      bicgstab_R: tol=%e maxiter=%d\n",
	   tol,maxiter);
    printf("iteration     residual\n");
  }

  bb = new_vector(A->n,A->mnl);
  p = new_vector(A->n,A->mnl);
  v = new_vector(A->n,A->mnl);
  w = new_vector(A->n,A->mnl);
  z = new_vector(A->n,A->mnl);
  r = new_vector(A->n,A->mnl);
  rhat = new_vector(A->n,A->mnl);

  /* new rhs */
  Av_proc(M,b,bb);

  /* compute initial residual */
  Av_proc(A,x,w);
  Av_proc(M,w,z);
  v_plus_cw(bb,z,c,r);
  rcopy_vv(r,rhat);
  res0 = norm(bb,A->comm);

  do {

    k++;

    /* compute new p */
    v_plus_cw(p,v,-omega,z);
    v_plus_cw(r,z,beta,p);

    /* compute new v, r, and alpha */
    Av_proc(A,p,w);
    Av_proc(M,w,v);

    rho = dot(rhat,r,A->comm);
    alpha = rho/dot(rhat,v,A->comm);
    v_plus_cw(r,v,-alpha,w);

    /* compute new omega */
    Av_proc(A,w,r);
    Av_proc(M,r,z);

    if ((tmp = dot(z,z,A->comm)) < 1.e-20)  omega = 1;
    else omega = dot(z,w,A->comm)/tmp;

    /* compute new x and new r */
    v_plus_cw(x,p,alpha,r);
    v_plus_cw(r,w,omega,x);
    v_plus_cw(w,z,-omega,r);
    beta = (alpha/omega)*(dot(rhat,r,A->comm)/rho);

    /* compute exact residual -> w */
    Av_proc(A,x,w);
    Av_proc(M,w,z);
    v_plus_cw(bb,z,c,w);
    res = norm(w,A->comm);

    if (!A->myid && verbose && (k%10 == 0 || verbose > 0)) {
      printf("%9d %.4e\n", k, res / res0);
      fflush(stdout);
    }
  }
  while ((res/res0 > tol) && (k < maxiter));

  if (!A->myid && verbose) printf("\n%9d %.6e\n", k, res/res0);

  free_vector(bb);
  free_vector(p);
  free_vector(v);
  free_vector(w);
  free_vector(z);
  free_vector(r);
  free_vector(rhat);

  if (verbose) {
    stop_timer(ident_bicgstab);
    report_times(ident_bicgstab,"bicgstab",0,A->comm);
  }
  return k;

}
