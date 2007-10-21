dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `cc, `sk'.
dnl Note that there isn't much error checking done on this file at present...

add(`allocator',                               `cc')
add(`dense_dense_float_sum',                   `cc')
add(`dense_dense_float_dot_product',           `cc')
add(`dense_dense_float_element_product',       `cc')
add(`dense_dense_float_matrix_vector_product', `cc')
add(`kernel_reference',                        `sk')
add(`kernel_test',                             `sk')
