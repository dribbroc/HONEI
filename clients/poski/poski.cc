#include <stdio.h>
#define HAVE_EPS_DOUBLE 0
#define  HAVE_C99_FPCLASSIFY 0
#include <poski/poski.h>

#define DIM 3
#define NUM_STORED_NZ 5

int main(int /*argc*/, char ** /*argv*/)
{
	printf("\n######################################################\n");
	printf("# A First Example:\n");
	printf("#   - This example computes SpMV (y=alpha*Ax + beta*y)\n");
	printf("######################################################\n");
	int Aptr[DIM+1]={0,1,3,5};
	int Aind[NUM_STORED_NZ]={0,0,1,0,2};
	double Aval[NUM_STORED_NZ]={1,-2,1,0.5,1};

	double x[DIM] = {.25, .45, .65};
	double y[DIM] = {1, 1, 1};
	double alpha=-1,beta=1;

	poski_mat_t A_tunable;
	poski_vec_t x_view, y_view;

	poski_Init();
	poski_threadarg_t *threadargs = poski_InitThreads();

	A_tunable = poski_CreateMatCSR(Aptr, Aind, Aval, DIM, DIM, NUM_STORED_NZ, SHARE_INPUTMAT, threadargs, NULL, 2, INDEX_ZERO_BASED, MAT_GENERAL);

	x_view = poski_CreateVecView(x, DIM, STRIDE_UNIT, NULL);
	y_view = poski_CreateVecView(y, DIM, STRIDE_UNIT, NULL);

	poski_MatMult(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view);

	poski_DestroyThreads(threadargs);
	poski_DestroyVec(x_view);
	poski_DestroyVec(y_view);
	poski_DestroyMat(A_tunable);

	poski_Close();

	printf(" Given matrix A:\n");
	//poski_report_sparse_CSR_va(Aptr, Aind, Aval, DIM, DIM, NUM_STORED_NZ);
	printf(" Given vectors:\n");
	printf("\t+ x = [ %.3lf; %.3lf; %.3lf]\n", x[0],x[1],x[2]);
	printf("\t+ y = [ 1.0; 1.0; 1.0]\n");
	printf(" Results:\n");
	printf("  Expect: y = [ 0.75; 1.05; 0.225 ]\n");
	printf("  Answer: y = [%.3lf; %.3lf; %.3lf ]\n",y[0],y[1],y[2]);
	printf("######################################################\n\n");

	return 0;
}

