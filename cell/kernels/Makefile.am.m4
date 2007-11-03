ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`filelist', `')dnl
define(`sourceslist', `')dnl
define(`cleanlist', `')dnl
define(`ppeobjlist', `')dnl
define(`add', `define(`filelist', filelist `$1')dnl
define(`sourceslist', sourceslist `$1.cc')dnl
define(`cleanlist', cleanlist `$1.body' `$1.func' `$1.cc')dnl
define(`ppeobjlist', ppeobjlist `libcell_a-$1.o' `libcell_a-$1-env.o')dnl
$1.cc : $1.sk $(top_srcdir)/misc/make_sk.bash $2.cc.in
	if ! $(top_srcdir)/misc/make_sk.bash $1.sk $2.cc.in ; then rm -f $`'@ ; exit 1 ; fi

$1-env.cc : env.cc.in $1
	sed -e "s/@NAME@/$1/g" \
	    -e "s/@BEGIN@/$$(spu-readelf -s $1 | sed -ne "/_end/s/^[^:]*:[^0]*\([^ ]*\).*/0x\1/p")/" \
	    -e "s/@END@/0x35000/" \
	    -e "/@HEADER@/r $(top_srcdir)/misc/generated-file.txt" \
	    -e "/@HEADER@/d" \
	    $< > $`'@

$1_SOURCES = $1.cc
$1_CXXFLAGS = -O1 -Wall -msafe-dma -fno-exceptions -fno-rtti
$1_LDADD = \
	$(top_srcdir)/cell/libla/libla_spe.a \
	$(top_srcdir)/cell/libutil/libutil_spe.a

libcell_a-$1.o : $1
	ppu-embedspu $1 $< $`'@

libcell_a-$1-env.o : $1-env.cc
	ppu-g++ -o $`'@ -c $< $(CXXFLAGS) $(AM_CXXFLAGS)
')dnl

include(`cell/kernels/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CXX = spu-g++

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in Makefile.am
DISTCLEANFILES = cleanlist
EXTRA_DIST = \
	Makefile.am.m4 \
	files.m4 \
	kernel.cc.in
BUILT_SOURCES = sourceslist
DEFS = \
	$(DEBUGDEF)

noinst_PROGRAMS = filelist
noinst_LIBRARIES = libcell.a

libcell_a_SOURCES =
libcell_a_LIBADD = ppeobjlist

Makefile.am : Makefile.am.m4 files.m4
	$(top_srcdir)/misc/do_m4.bash Makefile.am

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
