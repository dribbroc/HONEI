CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in configure config/* aclocal.m4 \
	config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash misc/do_m4.bash misc/generated-file.txt

SUBDIRS = benchmark libutil libla libmath unittest libgraph libswe

bench:
	$(MAKE) $(AM_MAKEFLAGS) all 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="bench-recursive" bench-recursive

bench-am:

quickcheck: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="quickcheck-recursive" quickcheck-recursive
quickcheck-sse: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="quickcheck-sse-recursive" quickcheck-sse-recursive
quickcheck-cell: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="quickcheck-cell-recursive" quickcheck-cell-recursive
quickcheck-mc: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="quickcheck-recursive" quickcheck-recursive
check-sse: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="check-sse-recursive" check-sse-recursive
check-cell: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="check-cell-recursive" check-sse-recursive
check-mc: 
	$(MAKE) $(AM_MAKEFLAGS) RECURSIVE_TARGETS+="check-mc-recursive" check-mc-recursive

doxygen:
	$(MAKE) -C doc $@
