ifdef(`__gnu__',`',`errprint(`This is not GNU m4...
')m4exit(1)') include(`misc/generated-file.txt')

dnl vim: set ft=m4 noet :

define(`filelist', `')dnl
define(`sklist', `')dnl
define(`skcleanlist', `')dnl
define(`sksourceslist', `')dnl
define(`ppeobjlist', `')dnl
define(`addcc', `define(`filelist', filelist `$1.cc')')dnl
define(`addsk', `define(`sklist', sklist `$1')dnl
define(`sksourceslist', sksourceslist `$1.cc')dnl
define(`skcleanlist', skcleanlist `$1.body' `$1.cc')dnl
define(`ppeobjlist', ppeobjlist `libla_ppe_a-$1.o')dnl
$1.cc : $1.sk $(top_srcdir)/misc/make_sk.bash kernel.cc.in
	if ! $(top_srcdir)/misc/make_sk.bash $1.sk kernel.cc.in ; then rm -f $`'@ ; exit 1 ; fi
$1_SOURCES = $1.cc
$1_CXXFLAGS = -O1 -Wall -msafe-dma -fno-exceptions -fno-rtti
$1_LDADD = $(top_srcdir)/libla/cell/libla_spe.a

libla_ppe_a-$1.o : $1
	ppu-embedspu $1 $< $`'@
')dnl
define(`addthis', `dnl
ifelse(`$2', `cc', `addcc(`$1')', `')dnl
ifelse(`$2', `sk', `addsk(`$1')', `')')dnl
define(`add', `addthis(`$1',`$2')')dnl

include(`libla/cell/files.m4')

AM_CXXFLAGS = -I$(top_srcdir)

CXX = spu-g++

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in Makefile.am
DISTCLEANFILES = skcleanlist
EXTRA_DIST = \
	Makefile.am.m4 \
	files.m4 \
	kernel.cc.in
BUILT_SOURCES = sksourceslist
DEFS = \
	$(DEBUGDEF)

noinst_PROGRAMS = sklist
noinst_LIBRARIES = libla_spe.a libla_ppe.a

libla_spe_a_SOURCES = filelist
libla_spe_a_CXXFLAGS = -O3 -Wall -msafe-dma -fno-exceptions -fno-rtti

libla_ppe_a_SOURCES =
libla_ppe_a_LIBADD = ppeobjlist

Makefile.am : Makefile.am.m4 files.m4
	$(top_srcdir)/misc/do_m4.bash Makefile.am
