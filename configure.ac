AC_INIT([ecpp], [1.0], [jean-pierre.flori@ssi.gouv.fr])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

dnl Automake
AM_INIT_AUTOMAKE()

dnl Libtool
LT_INIT([win32-dll])

dnl Magical macros
AC_CANONICAL_HOST
AC_PROG_CC
AC_PROG_CPP
AC_LANG([C])

dnl GMP
unset gmp_lib_path GMP_CFLAGS GMP_CC
AC_ARG_WITH(gmp_include,
   [  --with-gmp-include=DIR  GMP include directory ],
   CPPFLAGS="$CPPFLAGS -I$withval")
AC_ARG_WITH(gmp_lib,
   [  --with-gmp-lib=DIR      GMP lib directory ], [
   LDFLAGS="$LDFLAGS -L$withval"
  ])
AC_ARG_WITH(gmp,
   [  --with-gmp=DIR          GMP install directory ], [
   if test -z "$with_gmp_lib" && test -z "$with_gmp_include" ; then
      CPPFLAGS="$CPPFLAGS -I$withval/include"
      LDFLAGS="$LDFLAGS -L$withval/lib"
   else
      AC_MSG_FAILURE([Do not use --with-gmp and --with-gmp-include/--with-gmp-lib options simultaneously.])
   fi
  ])
AC_CHECK_LIB([gmp], [__gmpz_init])

dnl Windows shared library stuff
case $host in
  *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
   AC_MSG_CHECKING(for DLL/static GMP)
   if test "$enable_shared" = yes; then
     ECPP_LDFLAGS="$ECPP_LDFLAGS -no-undefined"
     LIBECPP_LDFLAGS="$LIBECPP_LDFLAGS"
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if !__GMP_LIBGMP_DLL
# error "Dead man"
error
#endif
     ]], [[]])],[AC_MSG_RESULT(DLL)],[
  AC_MSG_RESULT(static)
  AC_MSG_ERROR([gmp.h isn't a DLL: use --enable-static --disable-shared]) ])
   else
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if __GMP_LIBGMP_DLL
# error "Dead man"
error
#endif
     ]], [[]])],[AC_MSG_RESULT(static)],[
  AC_MSG_RESULT(DLL)
  AC_MSG_ERROR([gmp.h is a DLL: use --disable-static --enable-shared]) ])
  fi
  ;;
esac

AC_SUBST(ECPP_LDFLAGS)
AC_SUBST(LIBECPP_LDFLAGS)

dnl Output
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])
AC_OUTPUT
