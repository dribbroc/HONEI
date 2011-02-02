
# --------------------------------------------------------------------- #
#  Flags for `aclocal'
# --------------------------------------------------------------------- #

#ACLOCAL_AMFLAGS = -I m4


# --------------------------------------------------------------------- #
#  Conditional automake subdirectories.
# --------------------------------------------------------------------- #

if BLAPACK
else
   OPT_BLAS=blas
endif

if BLAPACK
else
   OPT_LAPACK=lapack
endif

#if BUILD_MATLAB_IF
#   OPT_MATLAB=matlab
#endif


# --------------------------------------------------------------------- #
#  Subdirectory and distribution declarations.
# --------------------------------------------------------------------- #

SUBDIRS = $(OPT_BLAS) $(OPT_LAPACK) src
#$(OPT_MATLAB) doc data

#EXTRA_DIST = autogen.sh


# --------------------------------------------------------------------- #
#  Copy dist files to webserver directory.
# --------------------------------------------------------------------- #

#TAR_NAME = $(PACKAGE_NAME)-$(PACKAGE_VERSION)

#WWW_DIR = /mnt/Enterprise/WebServer/Computational/Documents/software/spai
#WWW_FILES = $(TAR_NAME).tar.gz $(TAR_NAME).tar.bz2

#www: dist
#	cp $(WWW_FILES) $(WWW_DIR)
#	for file in $(WWW_FILES); do  chmod a+r $(WWW_DIR)/$$file;  done


# --------------------------------------------------------------------- #

MAINTAINERCLEANFILES = *~ */*~
