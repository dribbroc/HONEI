ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`filelist', `')dnl
define(`sourceslist', `')dnl
define(`cleanlist', `')dnl
define(`objlist', `')dnl
define(`add', `define(`filelist', filelist `$1')dnl
define(`sourceslist', sourceslist `$1.cc' `$1-registrator.cc')dnl
define(`cleanlist', cleanlist `$1.body' `$1.func' `$1.cc')dnl
define(`objlist', objlist `libcell-$1.o')dnl
$1.cc : $1.sk $(top_srcdir)/misc/make_sk.bash $2-kernel.cc.in
	if ! $(top_srcdir)/misc/make_sk.bash $1.sk $2-kernel.cc.in ; then rm -f $`'@ ; exit 1 ; fi

$1-registrator.cc : registrator.cc.in $1
	sed -e "s/@BEGIN@/$$(spu-readelf -s $1 | sed -ne "/_end/s/^[^:]*:[^0]*\([^ ]*\).*/0x\1/p")/" \
	    -e "s/@END@/0x35000/" \
	    -e "s/@IDENTIFIER@/$1/g" \
	    -e "s/@NAME@/$$(echo $1 | sed -e "s/kernel_//")/" \
	    -e "/@OPCODES@/r $1.caps" \
	    -e "/@OPCODES@/d" \
	    -e "s/@OPCODECOUNT@/$$(wc -l $1.caps | cut -d " " -f 1)/" \
	    -e "s/@TYPE@/kt_$2/g" \
	    -e "/@HEADER@/r $(top_srcdir)/misc/generated-file.txt" \
	    -e "/@HEADER@/d" \
	    $< > $`'@

$1_SOURCES = $1.cc
$1_CXXFLAGS = -O1 -Wall -msafe-dma -fno-exceptions -fno-rtti
$1_LDADD = \
	$(top_srcdir)/cell/libla/libla_spe.a \
	$(top_srcdir)/cell/libutil/libutil_spe.a

libcell-$1.o : $1
	ppu-embedspu $1_handle $< $`'@
	sed -e "s/@NAME@/libcell-$1.o/"\
	    libtool-hack.in > libcell-$1.lo
')dnl

include(`cell/kernels/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CXX = spu-g++

BUILT_SOURCES = sourceslist objlist
CLEANFILES = *~ *.body *.caps *.func *.cc *.o
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

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
