dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`boundary_init_fsi',               `hh', `test', `cuda')
add(`collide_stream',                  `hh', `test')
add(`collide_stream_grid',             `hh', `sse', `cuda', `cell', `test')
add(`collide_stream_fsi',              `hh', `test', `cuda')
add(`collide_stream_grid_regression',        `test')
add(`dc_advanced',                           `test')
add(`dc_advanced_grid',                      `test')
add(`dc_util',                         `hh')
add(`equilibrium_distribution',        `hh', `test')
add(`equilibrium_distribution_grid',   `hh', `sse', `cuda', `cell', `test')
add(`equilibrium_distribution_grid_regression',  `test')
add(`extraction_grid',                 `hh', `sse', `cuda', `cell', `test')
add(`extraction_grid_regression',            `test')
add(`force_grid',                      `hh', `test', `cuda')
add(`fluid_solid_interaction',               `test')
add(`grid',                            `hh')
add(`grid_packer',                     `hh', `test')
add(`grid_partitioner',                `hh', `test')
add(`lbm_limiter',                     `hh', `test')
add(`partial_derivative',              `hh', `test')
add(`scenario_collection',             `hh', `test')
add(`scan_conversion_fsi',             `hh', `test')
add(`solid',                           `hh', `test')
add(`solver_labswe',                   `hh', `test')
add(`solver_lbm_grid',                 `hh', `test')
add(`solver_lbm_fsi',                  `hh', `test')
add(`solver_lbm_grid_multi',                 `test')
add(`solver_lbm_grid_regression',            `test')
add(`solver_labnavsto',                `hh', `test')
add(`solver_labnavsto_grid',                 `test')
add(`source',                          `hh', `test')
add(`tags',                            `hh')
add(`update_velocity_directions_grid', `hh', `sse', `cuda', `cell')
add(`update_velocity_directions_grid_regression', `test')
