dnl Process this file with autoconf to produce a configure script.

AC_INIT(shapelib-fortran, 1.3, [dcesari@arpa.emr.it])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([foreign])

AC_LANG([C])
AC_C_BIGENDIAN
AC_CHECK_SIZEOF(void *)

if test "x$ac_cv_sizeof_void_p" = "x1"; then
  sizeof_ptr_c="shp_int_b"
elif test "x$ac_cv_sizeof_void_p" = "x2"; then
  sizeof_ptr_c="shp_int_s"
elif test "x$ac_cv_sizeof_void_p" = "x4"; then
  sizeof_ptr_c="shp_int_l"
elif test "x$ac_cv_sizeof_void_p" = "x8"; then
  sizeof_ptr_c="shp_int_ll"
else
  sizeof_ptr_c="__error__"
fi
AC_DEFINE_UNQUOTED(SIZEOF_PTR_C, $sizeof_ptr_c,
 "Fortran KIND defining an INTEGER containing a C pointer")

dnl Locate Fortran 90 compiler among these if not specified
AC_PROG_FC([gfortran g95 pgf90 ifort xlf90 f90])
F77=$FC
FFLAGS=$FCFLAGS
AC_FC_WRAPPERS()dnl Determine how to call C from Fortran

AC_CHECK_LIB([shp], [SHPOpen], , AC_MSG_ERROR([shp library not found]), )

AC_PROG_LIBTOOL

pkgdocdir='${datadir}/doc/${PACKAGE}-${VERSION}'
AC_SUBST(pkgdocdir)

AC_CONFIG_FILES([
Makefile
])

AC_OUTPUT
