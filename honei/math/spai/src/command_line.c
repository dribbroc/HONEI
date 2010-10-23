/*
   SPAI MPI/C version @Copyright 1996,  All Rights Reserved
   Stephen Barnard
*/

#include "command_line.h"

/**********************************************************************/
/*
  parameter:     matrix file

The following data statement is used for help printout, and to define
all defaults != 0 (i.e. they must be defined there as shown) */
static char options[] =
"  options:   \n"
"          -bi:  read matrix as a binary file   \n"
"                default is 0                   \n"
"\n"
"          -bs:  block size                          \n"
"                default is bs = 1, variable is 0    \n"
"\n"
"          -cs:  cache size {0,1,2,3,4,5} in SPAI             \n"
"                default is cs = 5 (the biggest cache size)   \n"
"\n"
"          -db:  debugging level   \n"
"                default is 0      \n"
"\n"
"          -ep:  epsilon parameter for SPAI  \n"
"                default is ep = 0.6         \n"
"\n"
"          -ld:  no. of lower (below main) diagonals (for sc = 2)  \n"
"                default is 0                                      \n"
"\n"
"          -lp:  left preconditioning                    \n"
"                default is 0 => right preconditioning   \n"
"\n"
"          -mb:  size of various working buffers in SPAI  \n"
"                default is mb = 100000                   \n"
"\n"
"          -mf:  message file for warning messages   \n"
"                default is 0                        \n"
"\n"
"          -mi:  maximum number of iterations in BICGSTAB  \n"
"                default is mi = 500                       \n"
"\n"
"          -mn:  maximum number of new nonzero candidate per step  \n"
"                default is mn = 5                                 \n"
"\n"
"          -ns:  maximum number of improvement steps per row in SPAI  \n"
"                default is ns = 5                                    \n"
"\n"
"          -sc:  sparsity control: 0 = adaptive, 1 = specified tau, \n"
"                2 = fixed diagonals (-ud & -ld)                    \n"
"                (main diagonal always included)                    \n"
"                default is 0                                       \n"
"\n"
"          -sp:  symmetric pattern    \n"
"                default is 0         \n"
"\n"
"          -ta:  (tau) only those entries Mij included for which    \n"
"                |Aij| > (1-ta)*max_{j}|Aij|        (for -sc = 1)   \n"
"                i.e. ta=0: main diagonal, ta=1: full pattern of A  \n"
"                Main diagonal always included                      \n"
"                default is 0                                       \n"
"\n"
"          -to:  tolerance in BICGSTAB    \n"
"                default is to = 1.0e-8  \n"
"\n"
"          -ud:  no. of upper (above main) diagonals (for spar = 2)  \n"
"                default is 0                                        \n"
"\n"
"          -vb:  verbose: print parameters, timings and matrix statistics  \n"
"                for all values: write file solution.mm                    \n"
"                vb > 0: write matrix statistics, iteration process        \n"
"                vb > 1: write matrices A.mm and M.mm                      \n"
"                vb > 2: print timing after major steps                    \n"
"                vb > 3: add one line per execution to spai_summary.file   \n"
"                        with date, time, some input and run parameters    \n"
"                default is vb = 1                                         \n"
"\n\n"
;

int i_default(char* code)
{
  char *p;
  int val = 0;
  if ((p = strstr(options, code)))
    {
     p += strlen(code);
     sscanf(p, "%d", &val);
    }
  return val;
}

double d_default(char* code)
{
  char *p;
  double val = 0;
  if ((p = strstr(options, code)))
    {
     p += strlen(code);
     sscanf(p, "%le", &val);
    }
  return val;
}

void parameters
(int argc,
 char *argv[],
 char **matrix_file_ptr,
 char **rhs_file_ptr,
 double *ep_ptr,
 int *ns_ptr,
 int *mb_ptr,
 int *mn_ptr,
 int *cs_ptr,
 char **mf_ptr,
 int *mi_ptr,
 double *to_ptr,
 int *bs_ptr,
 int *sp_ptr,
 int *lp_ptr,
 int *bi_ptr,
 int *vb_ptr,
 int *db_ptr,
 int *sc_ptr,
 int *ld_ptr,
 int *ud_ptr,
 double *ta_ptr,
 MPI_Comm comm)

{
  int k;

  static char matrix_file[128];
  static char rhs_file[128];
  static double epsilon;
  static int nsteps;
  static int max_buf;
  static int max_new;
  static int cache_size;
  static char message_file[128] = "";
  static int max_iter;
  static double tolerance;
  static int symmetric_pattern = 0;
  static int block_size;
  static int left_precon = 0;
  static int binary = 0;
  static int verbose;
  static int debug_level = 0;
  static int spar = 0;
  static int lower_diag = 0;
  static int upper_diag = 0;
  static double tau = 0;

  static int found_matrix_file = 0;
  static int found_rhs_file = 0;

  int numprocs,myid;
#ifdef MPI
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myid);
  MPI_Barrier(comm);
#else
  numprocs = 1;
  myid = 0;
#endif
/* set defaults from definition above */
  epsilon = d_default("ep =");
  tolerance = d_default("to =");
  block_size = i_default("bs =");
  cache_size = i_default("cs =");
  max_buf = i_default("mb =");
  max_iter = i_default("mi =");
  max_new = i_default("mn =");
  nsteps = i_default("ns =");
  verbose = i_default("vb =");
  k=0;
/* Help required ? */
  if (argc < 2)
  {
   printf("\n usage: spai matrix.mm rhs.mm options\n");
   printf("    or: spai -h(elp) for help on input parameters\n\n");
   exit(0);
  }
  if (strncmp(argv[1], "-h", 2) == 0)
  {
   printf("\n input parameters (first is mandatory):\n\n");
   printf(" matrix:     matrix file in MatrixMarket format\n");
   printf(" rhs:        rhs file ditto (if missing rhs=A*1 used)\n");
   printf(" options (any order, always -code number, e.g. -ep 0.7):\n");
   printf("%s", options);

   exit(0);
  }
  /* First look for the matrix_file argument */
    if (argv[1][0] != '-') {
      found_matrix_file = 1;
      strcpy(matrix_file,argv[++k]);
      argc--;
    }

  /* Next look for the rhs_file argument */
  if (argc > 1) {
    if (argv[2][0] != '-') {
      found_rhs_file = 1;
      strcpy(rhs_file,argv[++k]);
      argc--;
    }
  }

  while (argc-- > 2) {

    switch (match_arg(argv[++k])) {

    case 0: /* -ep */
      epsilon = strtod(argv[++k],NULL);
      argc--;
      break;

    case 1: /* -ns */
      nsteps = strtod(argv[++k],NULL);
      argc--;
      break;

    case 2: /* -mb */
      max_buf = strtod(argv[++k],NULL);
      argc--;
      break;

    case 3: /* -mn */
      max_new = strtod(argv[++k],NULL);
      argc--;
      break;

    case 4: /* -cs */
      cache_size = strtod(argv[++k],NULL);
      argc--;
      break;

    case 5: /* -mf */
      strcpy(message_file,argv[++k]);
      argc--;
      break;

    case 6: /* -mi */
      max_iter = strtod(argv[++k],NULL);
      argc--;
      break;

    case 7: /* -to */
      tolerance = strtod(argv[++k],NULL);
      argc--;
      break;

    case 8: /* bs */
      block_size = strtod(argv[++k],NULL);
      argc--;
      break;

    case 9: /* sp */
      symmetric_pattern = strtod(argv[++k],NULL);
      argc--;
      break;

    case 10: /* lp */
      left_precon = strtod(argv[++k],NULL);
      argc--;
      break;

    case 11: /* bi */
      binary = strtod(argv[++k],NULL);
      argc--;
      break;

    case 12: /* vb */
      verbose = strtod(argv[++k],NULL);
      argc--;
      break;

    case 13: /* db */
      debug_level = strtod(argv[++k],NULL);
      argc--;
      break;

    case 14: /* sc */
      spar = strtod(argv[++k],NULL);
      argc--;
      break;

    case 15: /* ld */
      lower_diag = strtod(argv[++k],NULL);
      if (lower_diag < 0) lower_diag = 0;
      argc--;
      break;

    case 16: /* ud */
      upper_diag = strtod(argv[++k],NULL);
      if (upper_diag < 0) upper_diag = 0;
      argc--;
      break;

    case 17: /* ta */
      tau = strtod(argv[++k],NULL);
      if (tau < 0.) tau = 0;
      else if (tau > 1.) tau = 1;
      argc--;
      break;

    default:  /* illegal parameter */
      if (myid == 0) {
	fprintf(stderr,"illegal parameter: %s\n",argv[k]);
      }
      k++; argc--;
      break;

    }
  }

  if (found_matrix_file) *matrix_file_ptr = matrix_file;
  else *matrix_file_ptr = NULL;

  if (found_rhs_file) *rhs_file_ptr = rhs_file;
  else *rhs_file_ptr = NULL;

  *ep_ptr = epsilon;
  *ns_ptr = nsteps;
  *mb_ptr = max_buf;
  *mn_ptr = max_new;
  *cs_ptr = cache_size;
  *mf_ptr = message_file;
  *mi_ptr = max_iter;
  *to_ptr = tolerance;
  *bs_ptr = block_size;
  *sp_ptr = symmetric_pattern;
  *lp_ptr = left_precon;
  *bi_ptr = binary;
  *vb_ptr = verbose;
  *db_ptr = debug_level;
  *sc_ptr = spar;
  *ld_ptr = lower_diag;
  *ud_ptr = upper_diag;
  *ta_ptr = tau;

  if (verbose && (myid == 0)) {
    printf("numprocs                   %d\n",  numprocs);
    printf("matrix file                %s\n",  matrix_file);
    printf("rhs file                   %s\n",  rhs_file);
    printf("epsilon (ep)               %le\n", epsilon);
    printf("nbsteps (ns)               %d\n",  nsteps);
    printf("max (mb)                   %d\n",  max_buf);
    printf("maxnew (mn)                %d\n",  max_new);
    printf("cache size (cs)            %d\n",  cache_size);
    printf("message file (mf)          %s\n",  message_file);
    printf("max_iter (mi)              %d\n",  max_iter);
    printf("tolerance (to)             %le\n", tolerance);
    printf("block size (bs)            %d\n",  block_size);
    printf("symmetric pattern (sp)     %d\n",  symmetric_pattern);
    printf("left preconditioner (lp)   %d\n",  left_precon);
    printf("binary (bi)                %d\n",  binary);
    printf("verbose (vb)               %d\n",  verbose);
    printf("debug (db)                 %d\n",  debug_level);
    printf("sparsity parameter (sc)    %d\n",  spar);
    printf("lower_diag (ld)            %d\n",  lower_diag);
    printf("upper_diag (ud)            %d\n",  upper_diag);
    printf("tau (ta)                   %f\n",  tau);
  }

}

/**********************************************************************/

int match_arg(char *string)
{
  static char *param_table[] = {
    "-ep",
    "-ns",
    "-mb",
    "-mn",
    "-cs",
    "-mf",
    "-mi",
    "-to",
    "-bs",
    "-sp",
    "-lp",
    "-bi",
    "-vb",
    "-db",
    "-sc",
    "-ld",
    "-ud",
    "-ta",
    " "     /* must remain last - ends list */
  };
  int k;

  k = 0;
  do
    if (!strcmp(string,param_table[k]))
      return(k);
  while (*param_table[++k] != ' ');

  return(-1);
}

/**********************************************************************/

