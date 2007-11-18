dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`absolute',                      `hh', `cc', `test')
add(`algorithm',                     `hh')
add(`banded_matrix',                 `hh', `test')
add(`dense_matrix',                  `hh', `test')
add(`dense_matrix_tile',             `hh', `test')
add(`dense_vector',                  `hh', `cc', `test')
add(`dense_vector_range',            `hh', `cc', `test')
add(`dense_vector_slice',            `hh', `cc', `test')
add(`difference',                    `hh', `sse', `test')
add(`dot_product',                   `hh', `sse', `cell', `test', `benchmark')
add(`element_inverse',               `hh', `test', `benchmark')
add(`element_iterator',              `hh', `test')
add(`element_product',               `hh', `sse', `cell', `test', `benchmark')
add(`matrix',                        `hh')
add(`matrix_error',                  `hh', `cc')
add(`norm',                          `hh', `test')
add(`product',                       `hh', `sse', `cell', `test', `benchmark')
add(`reduction',                     `hh', `cell', `test', `benchmark')
add(`residual',                      `hh', `test')
add(`scale',                         `hh', `sse', `cell', `test', `benchmark')
add(`scaled_sum',                    `hh', `sse', `test', `benchmark')
add(`sparse_matrix',                 `hh', `test')
add(`sparse_vector',                 `hh', `cc', `test')
add(`sum',                           `hh', `sse', `cell', `test', `benchmark')
add(`vector',                        `hh')
add(`vector_error',                  `hh', `cc')
add(`vector_iterator',               `hh', `test')
