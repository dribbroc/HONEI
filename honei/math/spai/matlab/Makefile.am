MEXEXT = @MEXEXT@
SUFFIXES = .c .$(MEXEXT)

EXTRA_DIST = spai_full.c spai.m spaitau.m spaidiags.m m1.dat	\
		test_matlab_01.sh test_matlab_01.out

matlabdir   = $(prefix)/matlab/spai
matlab_DATA = spai_full.$(MEXEXT) spai.m spaitau.m spaidiags.m

CLEANFILES = spai_full.$(MEXEXT) test_matlab_01.check

.c.$(MEXEXT):
	@MEX@ -v -argcheck $< -I../src $(AM_CFLAGS) $(LDFLAGS) -L../src -lspai @LAPACK_LIBS@ @BLAS_LIBS@ @LIBS@ -lm

TESTS      = test_matlab_01.sh
