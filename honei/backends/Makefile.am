CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2

if CELL
  CELLDIR = cell
endif

if CUDA
  CUDADIR = cuda
endif

if SSE
  SSEDIR = sse
endif

if MPI
  MPIDIR = mpi
endif

SUBDIRS = $(SSEDIR) $(CELLDIR) $(CUDADIR) $(MPIDIR)
