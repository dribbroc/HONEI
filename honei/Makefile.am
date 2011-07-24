CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = honeirc

if VISUAL_OLD
  VISUALDIR = visual
endif

if MPI
  MPIDIR = mpi
endif

SUBDIRS = util backends la math graph swe lbm fem $(MPIDIR) $(VISUALDIR)
