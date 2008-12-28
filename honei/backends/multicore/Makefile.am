AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

lib_LTLIBRARIES = libhoneibackendsmulticore.la

libhoneibackendsmulticore_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la

libhoneibackendsmulticore_la_SOURCES = dispatch_policy.hh \
				 ticket.hh \
				 ticket.cc \
				 thread_pool.hh \
				 thread_pool.cc \
				 thread.hh \
				 thread.cc

thread_pool_TEST_SOURCES = thread_pool_TEST.cc

thread_pool_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	libhoneibackendsmulticore.la \
	$(top_builddir)/honei/util/libhoneiutil.la
	$(DYNAMIC_LD_LIBS)
thread_pool_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

libhoneibackendsmulticore_includedir = $(includedir)/honei/backends/multicore/
libhoneibackendsmulticore_include_HEADERS = dispatch_policy.hh ticket.hh thread_pool.hh thread.hh
TESTS = thread_pool_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/unittest/run.sh
check_PROGRAMS = $(TESTS)
