dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`banded_matrix',             `hh')
add(`dense_matrix',              `hh', `test')
add(`dense_vector',              `hh', `test')
add(`element_iterator',          `hh')
add(`matrix',                    `hh')
add(`scalar_product',            `hh', `cc')
add(`tags',                      `hh')
add(`vector',                    `hh')
add(`vector_error',              `hh', `cc')
add(`vector_iterator',           `hh')
add(`vector_norm',               `hh', `test')

