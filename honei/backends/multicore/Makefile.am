AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(DEBUGDEF) \
	$(BOOSTDEF) \
	$(PROFILERDEF)

lib_LTLIBRARIES = libhoneibackendsmulticore.la

libhoneibackendsmulticore_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la

libhoneibackendsmulticore_la_SOURCES = concurrent_list.hh \
				 concurrent_list-impl.hh \
				 dispatch_policy.hh \
				 numainfo.hh \
				 operation.hh \
				 ticket.hh \
				 ticket.cc \
				 thread_pool.hh \
				 thread_pool.cc \
				 thread_function.hh \
				 thread_function.cc \
				 thread_task.hh \
				 topology.hh \
				 topology.cc \
				 x86_spec.hh

numainfo_TEST_SOURCES = numainfo_TEST.cc

thread_pool_TEST_SOURCES = thread_pool_TEST.cc

topology_TEST_SOURCES = topology_TEST.cc

thread_pool_TEST_LDADD = \
	libhoneibackendsmulticore.la \
	$(top_builddir)/honei/util/libhoneiutil.la
	$(DYNAMIC_LD_LIBS)
thread_pool_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

topology_TEST_LDADD = \
libhoneibackendsmulticore.la \
$(top_builddir)/honei/util/libhoneiutil.la
$(DYNAMIC_LD_LIBS)
topology_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

numainfo_TEST_LDADD = \
libhoneibackendsmulticore.la \
$(top_builddir)/honei/util/libhoneiutil.la
$(DYNAMIC_LD_LIBS)
numainfo_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

libhoneibackendsmulticore_includedir = $(includedir)/honei/backends/multicore/
libhoneibackendsmulticore_include_HEADERS = concurrent_list.hh concurrent_list-impl.hh dispatch_policy.hh numainfo.hh operation.hh \
						 ticket.hh thread_pool.hh thread_function.hh thread_task.hh \
						 x86_spec.hh
TESTS = numainfo_TEST thread_pool_TEST topology_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/honei/util/run.sh
check_PROGRAMS = $(TESTS)
