dnl vim: set sw=4 sts=4 ts=4 noet ft=config foldmethod=marker foldmarker={{{,}}} :

AC_INIT([pg512], [0.0])
AC_PREREQ(2.5)
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(1.9)
AC_LANG([C++])

dnl {{{ basic toolchain checks
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
dnl }}}

dnl {{{ check for libebt
AC_MSG_CHECKING([for libebt])
AC_COMPILE_IFELSE([
#include <libebt/libebt.hh>
#include <string>
struct Tag { };
int main(int, char **)
{
	libebt::BacktraceContext<Tag> context("pg512!");
}
],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])
	 AC_MSG_ERROR([libebt (http://libebt.berlios.de/) is required])])
dnl }}}

dnl {{{ check for a sensible C++ compiler
dnl {{{ check for partial specialisation of class templates
AC_MSG_CHECKING([for compiler support for partial specialisation of template classes])
AC_COMPILE_IFELSE([
template <typename T_>
struct S
{
	static const int s = 0;
};

template <typename T_>
struct S<const T_>
{
};

template <typename T_>
struct S<const T_ &>
{
	static const int t = 0;
};

int main(int, char **)
{
	return S<int &>::s + S<const int &>::t;
}
],
		[AC_MSG_RESULT([yes])],
		[
		 AC_MSG_RESULT([no])
		 AC_MSG_ERROR([Your C++ compiler doesn't compile C++ code])
		])

dnl {{{ check for tr1 shared pointers
AC_MSG_CHECKING([for std::tr1::shared_ptr<>])
AC_COMPILE_IFELSE([
#include <tr1/memory>

int main(int, char **)
{
	std::tr1::shared_ptr<int> x(new int(10));
	return *x;
}
],
		[AC_MSG_RESULT([yes])],
		[
		 AC_MSG_RESULT([no])
		 AC_MSG_ERROR([Your C++ Standard Library doesn't support std::tr1::shared_ptr<>.])
		])
dnl }}}

dnl {{{ check for numeric limits
AC_MSG_CHECKING([for std::numeric_limits<>])
AC_COMPILE_IFELSE([
#include <limits>
int main(int, char **)
{
	int result(std::numeric_limits<float>::digits >> 3);
	result += (std::numeric_limits<double>::digits >> 3);

	return result;
}
],
	[AC_MSG_RESULT([yes])],
	[
	 AC_MSG_RESULT([no])
	 AC_MSG_ERROR([Your C++ Standard Library doesn't support std::numeric_limits<>.])
	])
dnl }}}
dnl }}}

dnl {{{ check for build support conditionals
dnl {{{ check for GPU support
AC_MSG_CHECKING([whether to build GPU support])
AC_ARG_ENABLE([gpu],
	AS_HELP_STRING([--enable-gpu], [Build gpu support (needs packages: libX11, libGL, libGLEW)]),
		[ENABLE_GPU=$enableval
		  AC_MSG_RESULT([$enableval])],
		[ENABLE_GPU=no
		  AC_MSG_RESULT([no])])
AC_SUBST([ENABLE_GPU])
AM_CONDITIONAL([GPU], test "x$ENABLE_GPU" = "xyes")
if test "x$ENABLE_GPU" = "xyes" ; then
	dnl {{{ check for libX11
	AC_CHECK_LIB(X11, XOpenDisplay,
		[have_x11=yes])
	if test "x$have_x11" != "xyes" ; then
		AC_MSG_ERROR([libX11 is needed for GPU support.])
	fi
	dnl }}}
	dnl {{{ check for libGL
	AC_CHECK_LIB(GL, glXCreateContext,
		[have_gl=yes])
	if test "x$have_gl" != "xyes" ; then
		AC_MSG_ERROR([libGL is needed for GPU support.])
	fi
	dnl }}}
	dnl {{{ check for libGLEW
	AC_CHECK_LIB(GLEW, glewInit,
		[have_glew=yes])
	if test "x$have_glew" != "xyes" ; then
		AC_MSG_ERROR([libGLEW is needed for GPU support.])
	fi
	dnl }}}
fi
dnl }}}
dnl {{{ check for debug support
AC_MSG_CHECKING([whether to build debug support])
AC_ARG_ENABLE([debug],
	AS_HELP_STRING([--enable-debug], [Build debug support (needs packages: none)]),
		[ENABLE_DEBUG=$enableval
		  AC_MSG_RESULT([$enableval])],
		[ENABLE_DEBUG=yes
		  AC_MSG_RESULT([yes])])
AC_SUBST([ENABLE_DEBUG])
if test "x$ENABLE_DEBUG" = "xyes"; then
	DEBUGDEF="-DDEBUG"
fi
AC_SUBST([DEBUGDEF])
AM_CONDITIONAL([DEBUG], test "x$ENABLE_DEBUG" = "xyes")
dnl }}}
dnl }}}

dnl output
AM_CONFIG_HEADER(config.h)
AC_OUTPUT(
	Makefile
	doc/Makefile
	libgraph/Makefile
	libla/Makefile
	libswe/Makefile
	libutil/Makefile
	unittest/Makefile
)
