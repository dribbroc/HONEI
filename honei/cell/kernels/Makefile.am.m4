ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') dnl include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`filelist', `')dnl
define(`sourceslist', `')dnl
define(`cleanlist', `')dnl
define(`objlist', `')dnl
define(`add', `define(`filelist', filelist `$1')dnl
define(`sourceslist', sourceslist `$1.cc' `$1-registrator.cc')dnl
define(`cleanlist', cleanlist `$1.body' `$1.functions' `$1.opcodes' `$1.cc')dnl
define(`objlist', objlist `libcell-$1.o')dnl
$1.cc : $1.sk $(top_srcdir)/misc/make_sk.bash $2-kernel.cc.in
	$(top_srcdir)/misc/make_sk.bash $1.sk
	sed -e "/@FUNCTIONS@/r $1.functions" \
	    -e "/@FUNCTIONS@/d" \
	    -e "/@BODY@/r $1.body" \
	    -e "/@BODY@/d" \
	    -e "/@HEADER@/r $(top_srcdir)/misc/generated-file.txt" \
	    -e "/@HEADER@/d" \
	    -e "/vim/s/set/set ro/" \
	    $2-kernel.cc.in > $`'@

$1-registrator.cc : registrator.cc.in $1
	sed -e "s/@BEGIN@/$$(spu-readelf -s $1 | sed -ne "/_end/s/^[^:]*:[^0]*\([^ ]*\).*/0x\1/p")/" \
	    -e "s/@END@/0x35000/" \
	    -e "s/@IDENTIFIER@/$1/g" \
	    -e "s/@NAME@/$$(echo $1 | sed -e "s/kernel_//")/" \
	    -e "/@OPCODES@/r $1.opcodes" \
	    -e "/@OPCODES@/d" \
	    -e "s/@OPCODECOUNT@/$$(wc -l $1.opcodes | cut -d " " -f 1)/" \
	    -e "s/@TYPE@/kt_$2/g" \
	    -e "/@HEADER@/r $(top_srcdir)/misc/generated-file.txt" \
	    -e "/@HEADER@/d" \
	    -e "/vim/s/set/set ro/" \
	    $< > $`'@

$1_SOURCES = $1.cc
$1_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti
$1_LDADD = \
	$(top_srcdir)/honei/cell/libutil/libutil_spe.a \
	$(top_srcdir)/honei/cell/libgraph/libgraph_spe.a \
	$(top_srcdir)/honei/cell/libla/libla_spe.a \
	$(top_srcdir)/honei/cell/libmath/libmath_spe.a \
	$(top_srcdir)/honei/cell/libswe/libswe_spe.a \
	$(top_srcdir)/honei/cell/libutil/libutil_spe.a

libcell-$1.o : $1
	$(PPU_EMBEDSPU) $1_handle $< $`'@
	sed -e "s/@NAME@/libcell-$1.o/"\
	    libtool-hack.in > libcell-$1.lo
')dnl

include(`honei/cell/kernels/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CXX = $(SPU_CXX)
CXXFLAGS = $(SPU_CXXFLAGS)

BUILT_SOURCES = sourceslist objlist
CLEANFILES = *~ *.body *.functions *.opcodes *.cc *.o
DISTCLEANFILES = cleanlist
MAINTAINERCLEANFILES = Makefile.in Makefile.am

EXTRA_DIST = \
	Makefile.am.m4 \
	files.m4 \
	registrator.cc.in \
	stand_alone-kernel.cc.in
DEFS = \
	$(DEBUGDEF)

noinst_PROGRAMS = filelist

Makefile.am : Makefile.am.m4 files.m4
	$(top_srcdir)/misc/do_m4.bash Makefile.am
