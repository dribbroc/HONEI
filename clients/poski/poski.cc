/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <stdio.h>
#define HAVE_EPS_DOUBLE 0
#define HAVE_C99_FPCLASSIFY 0
#include <poski/poski.h>

#include <honei/util/time_stamp.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/sparse_matrix_csr.hh>
#include <honei/math/matrix_io.hh>
#include <honei/util/configuration.hh>

#include <iostream>

//#define DIM 3
//#define NUM_STORED_NZ 5

int main(int argc, char ** argv)
{
#ifdef HONEI_MPI
    MPI_Init(&argc, &argv);
#endif
    Configuration::instance()->set_value("csr::blocksize", 1);

  std::string filename(HONEI_SOURCEDIR);
  filename += "/honei/math/testdata/poisson_advanced4/sort_2/";
  filename += "prol_7.ell";
  SparseMatrixELL<double> aell(MatrixIO<io_formats::ELL>::read_matrix(filename, double(0)));
  SparseMatrixCSR<double> acsr(aell);

  int DIM = acsr.rows();
  int NUM_STORED_NZ = acsr.used_elements();

  int * Aptr = new int[DIM+1];
  for (int i(0) ; i < DIM+1 ; ++i)
      Aptr[i] = acsr.Ar()[i];

  int * Aind = new int[NUM_STORED_NZ];
  double * Aval = new double[NUM_STORED_NZ];

  for (int i(0) ; i < NUM_STORED_NZ ; ++i)
  {
      Aind[i] = (acsr.Aj())[i];
      Aval[i] = (acsr.Ax())[i];
  }

  double x[DIM];
  double y[DIM];;
  for (int i(0) ; i < DIM ; ++i)
  {
      x[i] = 1;
      y[i] = 2;
  }
  double alpha=1,beta=1;

  printf("Matrix read in finished!\n");

  poski_mat_t A_tunable;
  poski_vec_t x_view, y_view;

  poski_Init();
  poski_threadarg_t *threadargs = poski_InitThreads();

  A_tunable = poski_CreateMatCSR(Aptr, Aind, Aval, DIM, DIM, NUM_STORED_NZ, SHARE_INPUTMAT, threadargs, NULL, 2, INDEX_ZERO_BASED, MAT_GENERAL);

  x_view = poski_CreateVecView(x, DIM, STRIDE_UNIT, NULL);
  y_view = poski_CreateVecView(y, DIM, STRIDE_UNIT, NULL);

  TimeStamp at, bt;
  at.take();
  for (unsigned long i(0); i < 1000 ; ++i)
      poski_MatMult(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view);
  bt.take();
  std::cout<<"TOE: "<<bt.total()-at.total()<<std::endl;

  poski_DestroyThreads(threadargs);
  poski_DestroyVec(x_view);
  poski_DestroyVec(y_view);
  poski_DestroyMat(A_tunable);

  poski_Close();

  //poski_report_sparse_CSR_va(Aptr, Aind, Aval, DIM, DIM, NUM_STORED_NZ);

#ifdef HONEI_MPI
    MPI_Finalize();
#endif

    delete[] Aptr;
    delete[] Aind;
    delete[] Aval;
  return 0;
}

