ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`benchmarklist', `')dnl
define(`addbench', `define(`benchmarklist', benchmarklist `$1_BENCHMARK')dnl
$1_BENCHMARK_SOURCES = $1_BENCHMARK.cc
$1_BENCHMARK_LDADD = \
	$(top_builddir)/benchmark/libbenchmark.a \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/math/libhoneimath.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/graph/libhoneigraph.la \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(DYNAMIC_LD_LIBS)
$1_BENCHMARK_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)
')dnl
define(`addthis', `dnl
ifelse(`$2', `bench', `addbench(`$1')', `')dnl
')dnl
define(`add', `addthis(`$1',`$2')addthis(`$1',`$3')addthis(`$1',`$4')addthis(`$1',`$5')addthis(`$1',`$6')addthis(`$1',`$7')')dnl

include(`benchmark/files.m4')

if CELL

CELLTESTLIBS = $(top_builddir)/honei/cell/libcell.la

endif

AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in Makefile.am
EXTRA_DIST = Makefile.am.m4 files.m4

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

noinst_LIBRARIES = libbenchmark.a
noinst_PROGRAMS = benchmarklist

libbenchmark_a_SOURCES = \
	benchmark.cc benchmark.hh

BENCHMARKS = benchmarklist

.PHONY: benchmark
benchmark: $(BENCHMARKS)
	@for b in $(BENCHMARKS) ; do \
	    echo ">>> $$b[BACKENDS=$(BACKENDS)]" ; \
	    ./$$b $(BACKENDS) ; \
	done

Makefile.am : Makefile.am.m4 files.m4
	$(top_srcdir)/misc/do_m4.bash Makefile.am
