dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`engine',               `hh', `test')
add(`enginegraph',          `hh', `test')
add(`evolving',                   `test')
add(`enginesolve',          `hh', `test')
add(`enginesolveremote',    `hh', `test')
add(`engine_implicit',      `hh', `test')
add(`engine_client',        `hh', `test')
add(`engine_server',        `hh', `test')
add(`graphrandom',          `test')
add(`solver_server',        `hh', `test')
add(`solver_client',        `hh', `test')
