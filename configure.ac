dnl vim: set sw=4 sts=4 ts=4 noet ft=config foldmethod=marker foldmarker={{{,}}} :

AC_INIT([honei], [0.2])
AC_PREREQ(2.5)
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(1.9)
AC_LANG([C++])

dnl {{{ basic toolchain checks
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
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

dnl {{{ check for libpthread
AC_CHECK_LIB(pthread, pthread_mutex_init,
	[have_pthread=yes])
if test "x$have_pthread" != "xyes" ; then
	AC_MSG_ERROR([libpthread is required.])
fi
dnl }}}

dnl {{{ check for build support conditionals
dnl {{{ check for CELL support
AC_MSG_CHECKING([whether to build CELL BE support])
AC_ARG_ENABLE([cell],
	AS_HELP_STRING([--enable-cell], [Build Cell BE support (needs packages: libspe2)]),
		[ENABLE_CELL=$enableval
		  AC_MSG_RESULT([$enableval])],
		[ENABLE_CELL=no
		  AC_MSG_RESULT([no])])
AC_SUBST([ENABLE_CELL])
AM_CONDITIONAL([CELL], test "x$ENABLE_CELL" = "xyes")
if test "x$ENABLE_CELL" = "xyes" ; then
	dnl {{{ check for libspe2
	AC_CHECK_LIB(spe2, spe_context_create,
		[have_spe2=yes])
	if test "x$have_spe2" != "xyes" ; then
		AC_MSG_ERROR([libspe2 (http://sourceforge.net/projects/libspe) is required for CELL BE support.])
	fi
	dnl }}}
	dnl {{{ check for toolchain
	dnl {{{ check for spu-g++
	if test "x$SPU_CXX" = "x" ; then
		SPU_CXX=spu-g++
	fi
	AC_CHECK_PROG(have_spu_cxx, [$SPU_CXX], [yes], [no])
	if test "x$have_spu_cxx" != "xyes" ; then
		AC_MSG_ERROR([The SPU C++ compiler ($SPU_CXX) is required for CELL BE support.])
	fi
	AC_SUBST([SPU_CXX])
	dnl }}}
	dnl {{{ check for ppu-embedspu
	if test "x$PPU_EMBEDSPU" = "x" ; then
		PPU_EMBEDSPU=ppu-embedspu
	fi
	AC_CHECK_PROG(have_ppu_embedspu, [$PPU_EMBEDSPU], [yes], [no])
	if test "x$have_ppu_embedspu" != "xyes" ; then
		AC_MSG_ERROR([The SPU embedding program ($PPU_EMBEDSPU) is required for CELL BE support.])
	fi
	AC_SUBST([PPU_EMBEDSPU])
	dnl }}}
	dnl }}}
	dnl {{{ enable CELL support
	CELLDEF="-DHONEI_CELL"
	if test "x$SPU_CXXFLAGS" = "x" ; then
		SPU_CXXFLAGS="-O1"
	fi
	AC_SUBST([CELLDEF])
	AC_SUBST([SPU_CXXFLAGS])
	dnl }}}
fi
dnl }}}
dnl {{{ check for GPU support
AC_MSG_CHECKING([whether to build GPU support])
AC_ARG_ENABLE([gpu],
	AS_HELP_STRING([--enable-gpu], [Build GPU support (needs packages: libX11, libGL, libGLEW)]),
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
dnl {{{ check for SSE support
AC_MSG_CHECKING([whether to build SSE support])
AC_ARG_ENABLE([sse],
	AS_HELP_STRING([--enable-sse], [Build SSE support (needs packages: none)]),
		[ENABLE_SSE=$enableval
		  AC_MSG_RESULT([$enableval])],
		[ENABLE_SSE=no
		  AC_MSG_RESULT([no])])
AC_SUBST([ENABLE_SSE])
AM_CONDITIONAL([SSE], test "x$ENABLE_SSE" = "xyes")
if test "x$ENABLE_SSE" = "xyes"; then
	dnl {{{ enable sse support
	SSEDEF="-DHONEI_SSE"
	dnl }}}
fi
AC_SUBST([SSEDEF])
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
dnl {{{ check for profiler support
AC_MSG_CHECKING([whether to build profiler support])
AC_ARG_ENABLE([profiler],
	AS_HELP_STRING([--enable-profiler], [Build profiler support (needs packages: none)]),
		[ENABLE_PROFILER=$enableval
		  AC_MSG_RESULT([$enableval])],
		[ENABLE_PROFILER=yes
		  AC_MSG_RESULT([yes])])
AC_SUBST([ENABLE_PROFILER])
if test "x$ENABLE_PROFILER" = "xyes"; then
	PROFILERDEF="-DHONEI_PROFILER"
fi
AC_SUBST([PROFILERDEF])
AM_CONDITIONAL([PROFILER], test "x$ENABLE_PROFILER" = "xyes")
dnl }}}
dnl }}}

dnl {{{ check for optional packages
dnl {{{ check for HDF5 support
AC_MSG_CHECKING([whether to build HDF5 wrappers])
AC_ARG_WITH([hdf5],
	AS_HELP_STRING([--with-hdf5], [Build HDF5 wrappers (needs packages: libhdf5)]),
		[WITH_HDF5=$withval
		  AC_MSG_RESULT([$withval])],
		[WITH_HDF5=no
		  AC_MSG_RESULT([no])])
AC_SUBST([WITH_HDF5])
AM_CONDITIONAL([HDF5], test "x$WITH_HDF5" = "xyes")
if test "x$WITH_HDF5" = "xyes"; then
	dnl {{{ check for libhdf5
	AC_CHECK_LIB(hdf5, H5open,
		[have_hdf5=yes])
	if test "x$have_hdf5" != "xyes" ; then
		AC_MSG_ERROR([libhdf5 (hdf.ncsa.uiuc.edu/HDF5/) not found.])
	fi
	dnl }}}
fi
dnl }}}
dnl {{{
AC_MSG_CHECKING([whether to build the visualisation engine])
AC_ARG_WITH([visual],
	AS_HELP_STRING([--with--visual], [Build visualisation engine (needs packages: libGL, libGLU, libglut)]),
		[WITH_VISUAL=$withval
		  AC_MSG_RESULT([$withval])],
		[WITH_HDF5=no
		  AC_MSG_RESULT([no])])
AC_SUBST([WITH_VISUAL])
AM_CONDITIONAL([VISUAL], test "x$WITH_VISUAL" = "xyes")
if test "x$WITH_VISUAL" = "xyes"; then
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
		AC_MSG_ERROR([libGL (http://www.opengl.org) not found.])
	fi
	dnl }}}
	dnl {{{ check for libGLU
	AC_CHECK_LIB(GLU, gluPerspective,
		[have_glu=yes])
	if test "x$have_glu" != "xyes" ; then
		AC_MSG_ERROR([libGLU (http://www.opengl.org) not found.])
	fi
	dnl }}}
	dnl {{{ check for libglut
	AC_CHECK_LIB(glut, glutPostRedisplay,
		[have_glut=yes])
	if test "x$have_glut" != "xyes" ; then
		AC_MSG_ERROR([libglut (http://www.opengl.org) not found.])
	fi
	dnl }}}
fi
dnl }}}
dnl }}}

dnl {{{ global automake variables
BACKENDS="cpu mc"
if test "x$ENABLE_CELL" == "xyes" ; then
	BACKENDS="$BACKENDS cell"
fi
if test "x$ENABLE_SSE" = "xyes"; then
	BACKENDS="$BACKENDS sse mc-sse"
fi
TYPE="full"
AC_SUBST([BACKENDS])
AC_SUBST([TYPE])
dnl }}}

dnl output
AM_CONFIG_HEADER(config.h)
AC_OUTPUT(
	Makefile
	benchmark/Makefile
	clients/Makefile
	clients/config/Makefile
	clients/swe/Makefile
	doc/Makefile
	honei/Makefile
	honei/backends/Makefile
	honei/backends/cell/Makefile
	honei/backends/cell/ppe/Makefile
	honei/backends/cell/spe/Makefile
	honei/backends/cell/spe/kernels/Makefile
	honei/backends/cell/spe/libgraph/Makefile
	honei/backends/cell/spe/libla/Makefile
	honei/backends/cell/spe/libmath/Makefile
	honei/backends/cell/spe/libswe/Makefile
	honei/backends/cell/spe/libutil/Makefile
	honei/backends/sse/Makefile
	honei/graph/Makefile
	honei/la/Makefile
	honei/lbm/Makefile
	honei/math/Makefile
	honei/swe/Makefile
	honei/util/Makefile
	honei/visual/Makefile
	unittest/Makefile
)
