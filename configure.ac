# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ([2.64])
AC_INIT([LatSim],[0.1],[antonin.portelli@me.com],[LatSim])
AC_CONFIG_AUX_DIR([.buildutils])
AC_CONFIG_SRCDIR([lib/Global.cpp])
AC_CONFIG_SRCDIR([examples/exScalOpBench.cpp])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_SILENT_RULES([yes])
AC_CONFIG_HEADERS([config.h])

# Checks for programs
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S

# Options
AC_ARG_WITH([LatCore],
    [AS_HELP_STRING([--with-LatCore=prefix],
    [use this option for a non-standard install prefix of the LatCore library])],
    [AM_CXXFLAGS="$AM_CXXFLAGS -I$with_LatCore/include"]
    [AM_LDFLAGS="$AM_LDFLAGS -L$with_LatCore/lib"],
    []
)

# Get compilers informations
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([C_COMP_VENDOR],["$ax_cv_c_compiler_vendor"],
			[vendor of C compiler that will compile the code])
AM_CONDITIONAL([CC_GNU],[test $ax_cv_c_compiler_vendor = "gnu"])
AM_CONDITIONAL([CC_INTEL],[test $ax_cv_c_compiler_vendor = "intel"])
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CXX_COMP_VENDOR],["$ax_cv_cxx_compiler_vendor"],
			[vendor of C++ compiler that will compile the code])
AM_CONDITIONAL([CXX_GNU],[test $ax_cv_cxx_compiler_vendor = "gnu"])
AM_CONDITIONAL([CXX_INTEL],[test $ax_cv_cxx_compiler_vendor = "intel"])
AX_GCC_VERSION
AC_DEFINE_UNQUOTED([GCC_VERSION],["$GCC_VERSION"],
			[version of gcc that will compile the code])
AX_GXX_VERSION
AC_DEFINE_UNQUOTED([GXX_VERSION],["$GXX_VERSION"],
			[version of g++ that will compile the code])
LX_FIND_MPI

CXXFLAGS="$AM_CXXFLAGS $MPI_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $MPI_CXXLDFLAGS $LDFLAGS"
LIBS="$LIBS $MPI_CXXLDFLAGS"

# Checks for libraries.
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([libm library not found])])
AC_CHECK_LIB([LatCore],[_ZN7LatCore12testFunctionEv],[],
             [AC_MSG_ERROR([LatCore library not found])])

# Checks for header files.
AC_HEADER_STDC

AC_SUBST([LIBS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile lib/Makefile examples/Makefile])
AC_OUTPUT
