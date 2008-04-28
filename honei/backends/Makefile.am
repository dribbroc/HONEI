CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2

if CELL
  CELLDIR = cell
endif

if SSE
  SSEDIR = sse
endif

SUBDIRS = $(SSEDIR) $(CELLDIR)
