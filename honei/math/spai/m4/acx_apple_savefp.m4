AC_DEFUN([_STAR_RESTFP_FIX],
   [AC_CACHE_CHECK([whether we need any library fixups],
       [star_cv_restfp_fixup],
       [AC_REQUIRE([AC_CANONICAL_BUILD])
        AC_REQUIRE([AC_PROG_CC])
        AC_REQUIRE([AC_PROG_F77])
	# MH 2005-10-04 15:46:  Check for any Darwin (not only 7).
        if expr $build_os : "darwin" >/dev/null; then
dnl Only affects OSX/Darwin
            # Following uses undocumented (but probably fairly stable)
            # autoconf internal variable.
	    # MH 2005-10-04 18:53:  Use documented G77 shell variable.
            if test "$G77" = yes; then
dnl The problem only affects g77/gcc, so we know we're dealing with these below
                AC_LANG_PUSH(C)
                rm -f conftest*
                star_cv_restfp_fixup=unknown
                AC_LANG_CONFTEST(AC_LANG_PROGRAM([], restFP()))
                { AC_TRY_COMMAND($CC -o conftest.x -S conftest.c)
                  test $ac_status = 0
                } &&
                sed "s/_restFP/restFP/g" conftest.x>conftest.s &&
                { AC_TRY_COMMAND($CC -c -o conftest.$ac_objext conftest.s)
                  test $ac_status = 0
                } || star_cv_restfp_fixup=broken
                AC_LANG_POP(C)
                if test $star_cv_restfp_fixup = broken; then
                    AC_MSG_WARN([unable to assemble restFP test])
                else
                    # Link this with the C compiler
                    AC_TRY_COMMAND($CC -o conftest conftest.$ac_objext)
                    _s_cstatus=$ac_status
                    # Link this with the Fortran compiler
                    AC_TRY_COMMAND($F77 -o conftest conftest.$ac_objext)
                    if test $_s_cstatus = 0 -a $ac_status = 0; then
                        # both compilers can compile it
                        star_cv_restfp_fixup=no
                    elif test $_s_cstatus != 0 -a $ac_status != 0; then
                        # neither compiler can compile it
                        star_cv_restfp_fixup=no
                    elif test $_s_cstatus = 0; then
                        # The C compiler can, but the Fortran cannot
                        star_cv_restfp_fixup=yes
                    else
                        # The C compiler can't compile, but the Fortran can.
                        # Haven't heard of this case!  Don't know what to do.
                        star_cv_restfp_fixup=broken
                    fi
                fi
                # Link with -lcc_dynamic.
                # See http://www.astro.gla.ac.uk/users/norman/note/2004/restFP/
                if test $star_cv_restfp_fixup = yes; then
                    AC_TRY_COMMAND($F77 -o conftest conftest.$ac_objext -lcc_dynamic)
                    if test $ac_status = 0; then
                        star_cv_restfp_fixup=cc_dynamic
                    fi
                fi
                if test $star_cv_restfp_fixup = yes; then
                    # ooops
                    AC_MSG_WARN([unable to solve restFP problem])
                    star_cv_restfp_fixup=broken
                fi
                rm -f conftest*
            elif test -z "$F77"; then
                # not g77, and indeed no Fortran at all
                star_cv_restfp_fixup=nofortran
            else
                # There is a Fortran, but it's not g77, so either there's no
                # problem, or it's a mixed-compiler problem that's harder
                # than we know how to deal with.  But presumably the user
                # has worked this out.
                star_cv_restfp_fixup=no
            fi
        else # !Darwin
            star_cv_restfp_fixup=no
        fi
        ])
   case $star_cv_restfp_fixup in
     cc_dynamic)
       # Add the required libraries to C_F77_... variables, which are
       # generated in the required places by (our) automake.
       C_F77LINK_MAGIC="-lcc_dynamic"
       ;;
     nofortran)
       AC_MSG_NOTICE([No Fortran in path, so presumably no g77/gcc library problems])
       ;;
     *) ;;
   esac
   AC_SUBST(C_F77LINK_MAGIC)
])# _STAR_RESTFP_FIX
