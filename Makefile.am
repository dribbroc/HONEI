CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in configure config/* aclocal.m4 \
	config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash misc/do_m4.bash misc/generated-file.txt misc/make_sk.bash misc/autotroll.m4 misc/autotroll.mk

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = doc honei clients benchmark

.PHONY: benchmark

benchmark:
	$(MAKE) all
	$(MAKE) -C benchmark benchmark $@

doxygen:
	$(MAKE) -C doc $@

pdf:
	$(MAKE) -C doc $@
