CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in configure config/* aclocal.m4 \
	config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash misc/do_m4.bash misc/generated-file.txt

SUBDIRS = honei clients unittest

.PHONY: benchmark

benchmark:
	$(MAKE) -C benchmark all $@

doxygen:
	$(MAKE) -C doc $@
