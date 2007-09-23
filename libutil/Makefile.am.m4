ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`general_filelist', `')dnl
define(`general_testlist', `')dnl
define(`general_headerlist', `')dnl
define(`gpu_filelist', `')dnl
define(`gpu_testlist', `')dnl
define(`gpu_headerlist', `')dnl
define(`cell_filelist', `')dnl
define(`cell_testlist', `')dnl
define(`cell_headerlist', `')dnl
define(`addtest', `define(`$1_testlist', $1_testlist `$2_TEST')dnl
$2_TEST_SOURCES = $2_TEST.cc
$2_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	libutil.la \
	$(DYNAMIC_LD_LIBS)
$2_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)
')dnl
define(`addhh', `define(`$1_filelist', $1_filelist `$2.hh')define(`$1_headerlist', $1_headerlist `$2.hh')')dnl
define(`addcc', `define(`$1_filelist', $1_filelist `$2.cc')')dnl
define(`addthis', `dnl
ifelse(`$3', `hh', `addhh(`$1',`$2')', `')dnl
ifelse(`$3', `cc', `addcc(`$1',`$2')', `')dnl
ifelse(`$3', `test', `addtest(`$1',`$2')', `')')dnl
define(`add', `addthis(`$1',`$2',`$3')addthis(`$1',`$2',`$4')addthis(`$1',`$2',`$5')')dnl

include(`libutil/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in Makefile.am
EXTRA_DIST = Makefile.am.m4 files.m4
DEFS = \
	$(DEBUGDEF)

if GPU

GPUSOURCES = gpu_filelist
GPUTESTS = gpu_testlist
GPUHEADERS = gpu_headerlist
GPULIBS = -lX11 -lGL -lGLEW

else

GPUSOURCES =
GPUTESTS =
GPUHEADERS =
GPULIBS =

endif

if CELL

CELLSOURCES = cell_filelist
CELLTESTS = cell_testlist
CELLHEADERS = cell_headerlist
CELLLIBS = -lspe2

else

CELLSOURCES =
CELLTESTS =
CELLHEADERS =
CELLLIBS =

endif

lib_LTLIBRARIES = libutil.la

libutil_la_SOURCES = general_filelist $(GPUSOURCES) $(CELLSOURCES)
libutil_la_LIBADD = \
	-lpthread \
	$(GPULIBS) \
	$(CELLLIBS)

pg512_includedir = $(includedir)/pg512/
pg512_include_HEADERS = general_headerlist $(GPUHEADERS) $(CELLHEADERS)

TESTS = general_testlist $(GPUTESTS) $(CELLTESTS)
TESTS_ENVIRONMENT = bash $(top_builddir)/unittest/run.sh

check_PROGRAMS = $(TESTS)
	
bench:
	
quickcheck: $(TESTS)
	$(MAKE) $(AM_MAKEFLAGS) TESTS_ENVIRONMENT="bash $(top_builddir)/unittest/run_quick.sh" check
	
Makefile.am : Makefile.am.m4 files.m4
	$(top_srcdir)/misc/do_m4.bash Makefile.am
