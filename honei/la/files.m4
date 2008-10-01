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
add(`algorithm',                     `hh', `test')
add(`band_iterator',                 `hh')
add(`banded_matrix',                 `fwd', `hh', `impl', `cc', `test')
add(`banded_matrix_q1',              `hh', `impl', `cc', `test')
add(`const_vector',                  `fwd', `hh', `cc', `impl', `test')
add(`dense_matrix',                  `fwd', `hh', `impl', `cc', `test')
add(`dense_matrix_tile',             `fwd', `hh', `test')
add(`dense_vector',                  `fwd', `hh', `impl', `cc', `test')
add(`dense_vector_base',             `hh')
add(`dense_vector_range',            `fwd', `hh', `impl', `cc', `test')
add(`dense_vector_slice',            `fwd', `hh', `impl', `cc', `test')
add(`difference',                    `hh', `sse', `cell', `cuda', `test')
add(`dot_product',                   `hh', `sse', `cell', `cuda', `test')
add(`element_inverse',               `hh', `sse', `cell', `cuda', `test')
add(`element_iterator',              `hh', `test')
add(`element_product',               `hh', `sse', `cell', `cuda', `test')
add(`matrix',                        `hh')
add(`matrix_error',                  `hh', `cc')
add(`norm',                          `hh', `fwd', `sse', `cell', `cuda', `test')
add(`product',                       `hh', `sse', `cell', `cuda', `test')
add(`reduction',                     `hh', `fwd', `sse', `cell', `test')
add(`residual',                      `hh', `test')
add(`scale',                         `hh', `sse', `cell', `cuda', `test')
add(`scaled_sum',                    `hh', `sse', `cell', `cuda', `test')
add(`sparse_matrix',                 `fwd', `hh', `test')
add(`sparse_vector',                 `hh', `impl', `cc', `test')
add(`sum',                           `hh', `sse', `cell', `cuda', `test')
add(`trace',                         `hh', `test')
add(`vector',                        `hh')
add(`vector_error',                  `hh', `cc')
add(`vector_iterator',               `hh', `test')
