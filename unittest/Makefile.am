unittest_TEST_SOURCES = unittest_TEST.cc
unittest_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(DYNAMIC_LD_LIBS)
unittest_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
EXTRA_DIST = run.sh
DEFS = \
	$(DEBUGDEF)

noinst_LIBRARIES = libunittest.a
libunittest_a_SOURCES = \
	unittest.cc unittest.hh

TESTS = unittest_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/unittest/run.sh
check_PROGRAMS = $(TESTS)


