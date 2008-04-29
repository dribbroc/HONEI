ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`sourceslist', `')dnl
define(`objlist', `')dnl
define(`add', `define(`sourceslist', sourceslist `kernels/$1-registrator.cc')dnl
define(`objlist', objlist `kernels/libcell-$1.lo')')dnl

include(`honei/backends/cell/spe/kernels/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~ *.body *.caps *.func *.cc
MAINTAINERCLEANFILES = Makefile.in Makefile.am
SUBDIRS = libutil libgraph libla libmath libswe kernels .
DISTCLEANFILES = cleanlist
EXTRA_DIST = \
	Makefile.am.m4
DEFS = \
	$(DEBUGDEF)

lib_LTLIBRARIES = libhoneibackendscellspe.la
libhoneibackendscellspe_la_SOURCES = sourceslist
libhoneibackendscellspe_la_LIBADD = objlist

Makefile.am : Makefile.am.m4 kernels/files.m4
	$(top_srcdir)/misc/do_m4.bash Makefile.am
