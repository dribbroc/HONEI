dnl vim: set sw=8 sts=8 noet :

AC_INIT([pg512], [0.0])
AC_PREREQ(2.5)
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(1.9)

dnl check for required programs
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET

dnl check for cblas support
AC_CHECK_LIB(f77blas, dgemv_, [DYNAMIC_LD_LIBS="-lf77blas" have_dl=yes])

dnl output
AM_CONFIG_HEADER(config.h)
AC_OUTPUT(
	Makefile
	doc/Makefile
	libla/Makefile
	libutil/Makefile
)
