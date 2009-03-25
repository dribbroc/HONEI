CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2

  MULTICOREDIR = multicore

if CELL
  CELLDIR = cell
endif

if CUDA
  CUDADIR = cuda
endif

if SSE
  SSEDIR = sse
endif

if OPENCL
  OPENCLDIR = opencl
endif

if ITANIUM
  ITANIUMDIR = itanium
endif

if MPI
  MPIDIR = mpi
endif

SUBDIRS = $(SSEDIR) $(CELLDIR) $(CUDADIR) $(MPIDIR) $(MULTICOREDIR) $(OPENCLDIR) $(ITANIUMDIR)
