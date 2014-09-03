# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.00])
AC_INIT([kesitcik], [1.00], [onur@kreix.com])

#LT_INIT()
LT_INIT([disable-static])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
#AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([enable])
: ${CFLAGS=""}
: ${CXXFLAGS=""}

#AC_CONFIG_SRCDIR([include/prflag.h])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
dnl AC_PROG_CXXCPP
AC_PROG_CXX
dnl AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.

# Checks for header files.
AC_PATH_X
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_OUTPUT(Makefile)
#AC_CONFIG_FILES([Makefile
                 #main/makefile])
#AC_CONFIG_SUBDIRS([bak])
#AC_OUTPUT
