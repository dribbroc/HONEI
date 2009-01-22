dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`collide_stream',                  `hh', `test')
add(`collide_stream_grid',             `hh', `sse', `cuda', `test')
add(`dc_advanced',                           `test')
add(`equilibrium_distribution',        `hh', `test')
add(`equilibrium_distribution_grid',   `hh', `sse', `cuda', `test')
add(`extraction_grid',                 `hh', `sse', `cuda')
add(`force_grid',                      `hh', `cuda')
add(`grid',                            `hh')
add(`grid_packer',                     `hh', `test')
add(`grid_partitioner',                `hh', `test')
add(`partial_derivative',              `hh', `test')
add(`solver_labswe',                   `hh', `test')
add(`solver_labswe_grid',              `hh', `test')
add(`solver_labswe_grid_multi',              `test')
add(`solver_labswe_grid_multi_regression',   `test')
add(`solver_labswe_grid_regression',         `test')
add(`solver_labnavsto',                `hh', `test')
add(`source',                          `hh', `test')
add(`tags',                            `hh')
add(`update_velocity_directions_grid', `hh', `cuda')
