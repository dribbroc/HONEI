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

SUBDIRS = $(CELLDIR) util libla math graph swe liblbm $(VISUALDIR)

bench:
	$(MAKE) $(AM_MAKEFLAGS) all
	$(MAKE) $(AM_MAKEFLAGS) benchm-recursive

bench-am:

benchm-am:

quickcheck-am:
quickcheck:
	$(MAKE) $(AM_MAKEFLAGS) quickcheck-recursive
quickcheck-sse-am:
quickcheck-sse:
	$(MAKE) $(AM_MAKEFLAGS) quickcheck-sse-recursive
quickcheck-cell-am:
quickcheck-cell:
	$(MAKE) $(AM_MAKEFLAGS) quickcheck-cell-recursive
quickcheck-mc-am:
quickcheck-mc:
	$(MAKE) $(AM_MAKEFLAGS) quickcheck-mc-recursive
check-sse-am:
check-sse:
	$(MAKE) $(AM_MAKEFLAGS) check-sse-recursive
check-cell-am:
check-cell:
	$(MAKE) $(AM_MAKEFLAGS) check-cell-recursive
check-mc-am:
check-mc:
	$(MAKE) $(AM_MAKEFLAGS) check-mc-recursive

