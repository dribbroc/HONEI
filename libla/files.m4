dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`banded_matrix',             `hh', `test')
add(`dense_matrix',              `hh', `test')
add(`dense_vector',              `hh', `test')
add(`element_iterator',          `hh', `test')
add(`matrix',                    `hh')
add(`matrix_difference',         `hh', `test')
add(`matrix_error',              `hh', `cc')
add(`matrix_element_inverse',    `hh', `test')
add(`matrix_elementwise_product',`hh', `test')
add(`matrix_element_sum',        `hh', `test')
add(`matrix_product',            `hh', `test')
add(`matrix_row_sum_vector',     `hh', `test')
add(`matrix_sum',                `hh', `test')
add(`scalar_matrix_product',     `hh', `test')
add(`scalar_matrix_sum',         `hh', `test')
add(`scalar_product',            `hh', `cc', `test')
add(`scalar_vector_product',     `hh', `test')
add(`scalar_vector_sum',         `hh', `test')
add(`sparse_vector',             `hh', `test')
add(`vector',                    `hh')
add(`vector_element_sum',        `hh', `test')
add(`vector_absolute',           `hh', `cc', `test')
add(`vector_difference',         `hh', `test')
add(`vector_error',              `hh', `cc')
add(`vector_iterator',           `hh')
add(`vector_mask',               `hh', `test')
add(`vector_masked_max_index',   `hh', `test')
add(`vector_masked_min_index',   `hh', `test'
add(`vector_matrix_product',     `hh')
add(`vector_norm',               `hh', `test')
add(`vector_product',            `hh')
add(`vector_scaled_sum',         `hh', `test')
add(`vector_sum',                `hh', `test')
