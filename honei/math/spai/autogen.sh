#! /bin/sh

# ===================================================================== #
#									#
#  SPAI autogen.sh file.						#
#									#
#  Initializes the autoconf/automake environment.			#
#									#
# ===================================================================== #

# Options for various commands.
#

ACLOCAL_INCLUDE="-I m4"
AUTOMAKE_OPTIONS="-ac"


# --------------------------------------------------------------------- #

echo "Running aclocal ..."
aclocal ${ACLOCAL_INCLUDE}

echo "Running autoheader ..."
autoheader

echo "Running autoconf ..."
autoconf

# echo "Running libtoolize ..."
# libtoolize --force

echo "Running automake ${AUTOMAKE_OPTIONS} ..."
automake ${AUTOMAKE_OPTIONS}

# if test -x config.status -a "$#" = "0"; then
#     echo "Configuring with previous options ..."
#     ./config.status --recheck
# else
#     echo "Configuring ..."
#     ./configure $*
# fi
