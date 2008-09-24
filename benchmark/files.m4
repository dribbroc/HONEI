dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`breadth_first_search',                      `bench')
add(`collide_stream',                            `bench')
add(`difference',                                `bench')
add(`dot_product',                               `bench')
add(`element_inverse',                           `bench')
add(`element_product',                           `bench')
add(`equilibrium_distribution',                  `bench')
add(`graph',                                     `bench')
add(`memory_arbiter',                            `bench')
add(`node_distance',                             `bench')
add(`poisson_cg_double_dense',                   `bench')
add(`poisson_cg_float_banded',                   `bench')
add(`poisson_iteref_cg_double_dense',            `bench')
add(`poisson_iteref_cg_float_banded',            `bench')
add(`poisson_iteref_pcg_double_dense',           `bench')
add(`poisson_iteref_pcg_float_banded',           `bench')
add(`poisson_jac_double_dense',                  `bench')
add(`poisson_jac_float_banded',                  `bench')
add(`poisson_jac_float_banded_q1',               `bench')
add(`poisson_jac_double_banded_q1',              `bench')
add(`poisson_jackernel_double_banded',           `bench')
add(`poisson_jackernel_double_banded_sse',       `bench')
add(`poisson_jackernel_float_banded',            `bench')
add(`poisson_jkc_float_banded',                  `bench')
add(`poisson_jkc_double_banded',                 `bench')
add(`poisson_jackernel_float_banded_sse',        `bench')
add(`poisson_mixedprec_cg',                      `bench')
add(`poisson_mg_float_banded',                   `bench')
add(`poisson_pcg_double_dense',                  `bench')
add(`poisson_pcg_float_banded',                  `bench')
add(`position',                                  `bench')
add(`product',                                   `bench')
add(`reduction',                                 `bench')
add(`relax_solver',                              `bench')
add(`relax_solver_mp1',                          `bench')
add(`relax_solver_mp2',                          `bench')
add(`scale',                                     `bench')
add(`scaled_sum',                                `bench')
add(`solver',                                    `bench')
add(`solver_labswe',                             `bench')
add(`source',                                    `bench')
add(`sum',                                       `bench')
