#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mmio.h"

int main(int argc, char *argv[])
{
  FILE *ascii;
  FILE *binary;
  char *ascii_name;
  char *binary_name;
  char line[128];

  int M,N,nnz,k,ii,rowcol[2];
  double val;

  void read_mm_matrix_header
    (FILE *,
     int *,
     int *N,
     int *);

  if (argc < 3) {
    fprintf(stderr,"Not enough arguments.\n");
    exit(1);
  }

  ascii_name =   argv[1];
  binary_name =  argv[2];

  ascii = fopen(ascii_name,"r");
  if (! ascii) {
    fprintf(stderr,"error opening file: %s\n",ascii_name);
    exit(1);
  }

  /* Open the binary output file */
  binary = fopen(binary_name,"wb");
  if (! binary) {
    fprintf(stderr,"error opening file: %s\n",binary_name);
    exit(1);
  }

  read_mm_matrix_header(ascii, &M,&N,&nnz);

  /* write binary header */
  fwrite(&M,sizeof(int),1,binary);
  fwrite(&M,sizeof(int),1,binary);
  fwrite(&nnz,sizeof(int),1,binary);

  for (k=0; k<nnz; k++) {
    fgets(line,128,ascii);

    /* Change ',' to ' ' */
    for (ii=0; line[ii]; ii++) if (line[ii] == ',') line[ii] = ' ';

    sscanf(line,"%d %d %le\n", &rowcol[0],&rowcol[1],&val);
    fwrite(rowcol,sizeof(int),2,binary);
    fwrite(&val,sizeof(double),1,binary);
  }

  fclose(ascii);
  fclose(binary);

  exit(EXIT_SUCCESS);
}

/**********************************************************************/

void read_mm_matrix_header
(FILE *f,
 int *M_ptr,
 int *N_ptr,
 int *nnz_ptr)
{
    MM_typecode matcode;
    int M, N, nnz;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (!
	(mm_is_real(matcode)) &&
	(mm_is_matrix(matcode)) &&
	(mm_is_sparse(matcode)) &&
	(mm_is_general(matcode))) {
      printf("MM matrix must be real, sparse, and general\n");
      exit(1);
    }

    if (mm_read_mtx_crd_size(f, M_ptr, N_ptr, nnz_ptr)) {
      printf("error in  mm_read_mtx_crd_size\n");
      exit(1);
    }
}
