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

doxygen:
	$(MAKE) -C doc $@
