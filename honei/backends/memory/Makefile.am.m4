ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`general_filelist', `')dnl
define(`general_testlist', `')dnl
define(`general_headerlist', `')dnl
define(`cuda_filelist', `')dnl
define(`cuda_testlist', `')dnl
define(`cuda_headerlist', `')dnl
define(`addtest', `define(`$1_testlist', $1_testlist `$2_TEST')dnl
$2_TEST_SOURCES = $2_TEST.cc
$2_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	libhoneibackendsmemory.la \
	$(DYNAMIC_LD_LIBS)
	$(top_builddir)/honei/util/libhoneiutil.la
$2_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)
')dnl
define(`addhh', `define(`$1_filelist', $1_filelist `$2.hh')define(`$1_headerlist', $1_headerlist `$2.hh')')dnl
define(`addimpl', `define(`$1_filelist', $1_filelist `$2-impl.hh')')dnl
define(`addcc', `define(`$1_filelist', $1_filelist `$2.cc')')dnl
define(`addthis', `dnl
ifelse(`$3', `hh', `addhh(`$1',`$2')', `')dnl
ifelse(`$3', `impl', `addimpl(`$1', `$2')', `')dnl
ifelse(`$3', `cc', `addcc(`$1',`$2')', `')dnl
ifelse(`$3', `test', `addtest(`$1',`$2')', `')')dnl
define(`add', `addthis(`$1',`$2',`$3')addthis(`$1',`$2',`$4')addthis(`$1',`$2',`$5')')dnl

include(`honei/backends/memory/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in Makefile.am
EXTRA_DIST = Makefile.am.m4 files.m4
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

if CUDA

CUDASOURCES = cuda_filelist
CUDATESTS = cuda_testlist
CUDAHEADERS = cuda_headerlist
CUDALIBS = $(top_builddir)/honei/backends/cuda/libhoneibackendscuda.la

else

CUDASOURCES =
CUDATESTS =
CUDAHEADERS =
CUDALIBS =

endif


lib_LTLIBRARIES = libhoneibackendsmemory.la

libhoneibackendsmemory_la_SOURCES = general_filelist $(CUDASOURCES)
libhoneibackendsmemory_la_LIBADD = \
	$(CUDALIBS) \
	$(top_builddir)/honei/util/libhoneiutil.la

libhoneibackendsmemory_includedir = $(includedir)/honei/backends/memory
libhoneibackendsmemory_include_HEADERS = general_headerlist $(CUDAHEADERS)

TESTS = general_testlist $(CUDATESTS)
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/unittest/run.sh

check_PROGRAMS = $(TESTS)

Makefile.am : Makefile.am.m4 files.m4
	cd $(top_srcdir) ; ./misc/do_m4.bash honei/backends/memory/Makefile.am
