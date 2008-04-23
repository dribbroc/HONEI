CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = honeirc

if CELL
  CELLDIR = cell
endif

if VISUAL
  VISUALDIR = visual
endif

SUBDIRS = $(CELLDIR) util la math graph swe lbm $(VISUALDIR)
