#################################################################################################
# Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
# Produced at the Lawrence Livermore National Laboratory  
# Written by Todd Gamblin, tgamblin@llnl.gov.
# LLNL-CODE-417602
# All rights reserved.  
# 
# This file is part of Libra. For details, see http://github.com/tgamblin/libra.
# Please also read the LICENSE file for further information.
# 
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
# 
#  * Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the disclaimer below.
#  * Redistributions in binary form must reproduce the above copyright notice, this list of
#    conditions and the disclaimer (as noted below) in the documentation and/or other materials
#    provided with the distribution.
#  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################################

#
# LX_FIND_MPI()
#  ------------------------------------------------------------------------
# This macro finds an MPI compiler and extracts includes and libraries from 
# it for use in automake projects.  The script exports the following variables:
#
# AC_DEFINE variables:
#     HAVE_MPI         AC_DEFINE'd to 1 if we found MPI
#
# AC_SUBST variables:
#     MPICC            Name of MPI compiler
#     MPI_CFLAGS       Includes and defines for MPI C compilation
#     MPI_CLDFLAGS     Libraries and library paths for linking MPI C programs
#
#     MPICXX           Name of MPI C++ compiler
#     MPI_CXXFLAGS     Includes and defines for MPI C++ compilation
#     MPI_CXXLDFLAGS   Libraries and library paths for linking MPI C++ programs
#
#     MPIF77           Name of MPI Fortran 77 compiler
#     MPI_F77FLAGS     Includes and defines for MPI Fortran 77 compilation
#     MPI_F77LDFLAGS   Libraries and library paths for linking MPI Fortran 77 programs
# 
#     MPIFC            Name of MPI Fortran compiler
#     MPI_FFLAGS       Includes and defines for MPI Fortran compilation
#     MPI_FLDFLAGS     Libraries and library paths for linking MPI Fortran programs
# 
# Shell variables output by this macro:
#     have_C_mpi       'yes' if we found MPI for C, 'no' otherwise
#     have_CXX_mpi     'yes' if we found MPI for C++, 'no' otherwise
#     have_F77_mpi     'yes' if we found MPI for F77, 'no' otherwise
#     have_F_mpi       'yes' if we found MPI for Fortran, 'no' otherwise
#
AC_DEFUN([LX_FIND_MPI], 
[
     AC_LANG_CASE(
     [C], [
         AC_REQUIRE([AC_PROG_CC])
         if [[ ! -z "$MPICC" ]]; then
             LX_QUERY_MPI_COMPILER(MPICC, [$MPICC], C)
         else
             LX_QUERY_MPI_COMPILER(MPICC, [openmpicc mpicc mpiicc mpixlc mpipgcc], C)
         fi
     ],
     [C++], [    
         AC_REQUIRE([AC_PROG_CXX])
         if [[ ! -z "$MPICXX" ]]; then
             LX_QUERY_MPI_COMPILER(MPICXX, [$MPICXX], CXX)
         else
             LX_QUERY_MPI_COMPILER(MPICXX, [openmpicxx mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC], CXX)
         fi
     ],
     [F77], [
         AC_REQUIRE([AC_PROG_F77])
         if [[ ! -z "$MPIF77" ]]; then
             LX_QUERY_MPI_COMPILER(MPIF77, [$MPIF77], F77)
         else
             LX_QUERY_MPI_COMPILER(MPIF77, [openmpif77 mpif77 mpiifort mpixlf77 mpixlf77_r], F77)
         fi
     ],
     [Fortran], [
         AC_REQUIRE([AC_PROG_FC])
         if [[ ! -z "$MPIFC" ]]; then
             LX_QUERY_MPI_COMPILER(MPIFC, [$MPIFC], F)
         else
             mpi_default_fc="mpif95 openmpif90 mpif90 mpigfortran mpif2003"
             mpi_intel_fc="mpiifort"
             mpi_xl_fc="mpixlf95 mpixlf95_r mpixlf90 mpixlf90_r mpixlf2003 mpixlf2003_r"
             mpi_pg_fc="mpipgf95 mpipgf90"
             LX_QUERY_MPI_COMPILER(MPIFC, [$mpi_default_fc $mpi_intel_fc $mpi_xl_fc $mpi_pg_fc], F)
         fi
     ])
])


#
# LX_QUERY_MPI_COMPILER([compiler-var-name], [compiler-names], [output-var-prefix])
#  ------------------------------------------------------------------------
# AC_SUBST variables:
#     MPI_<prefix>FLAGS       Includes and defines for MPI compilation
#     MPI_<prefix>LDFLAGS     Libraries and library paths for linking MPI C programs
# 
# Shell variables output by this macro:
#     found_mpi_flags         'yes' if we were able to get flags, 'no' otherwise
#
AC_DEFUN([LX_QUERY_MPI_COMPILER],
[
     # Try to find a working MPI compiler from the supplied names
     AC_PATH_PROGS($1, [$2], [not-found])
     
     # Figure out what the compiler responds to to get it to show us the compile
     # and link lines.  After this part of the macro, we'll have a valid 
     # lx_mpi_command_line
     printf "checking whether $$1 responds to '-showme:compile'... "
     lx_mpi_compile_line=`$$1 -showme:compile 2>/dev/null`
     if [[ "$?" -eq 0 ]]; then
         echo yes
         lx_mpi_link_line=`$$1 -showme:link 2>/dev/null`
     else
         echo no
         printf "checking whether $$1 responds to '-showme'... "
         lx_mpi_command_line=`$$1 -showme 2>/dev/null`
         if [[ "$?" -ne 0 ]]; then
             echo no
             printf "checking whether $$1 responds to '-compile-info'... "
             lx_mpi_compile_line=`$$1 -compile-info 2>/dev/null`
             if [[ "$?" -eq 0 ]]; then
                 echo yes
                 lx_mpi_link_line=`$$1 -link-info 2>/dev/null`
             else
                 echo no
                 printf "checking whether $$1 responds to '-show'... "
                 lx_mpi_command_line=`$$1 -show 2>/dev/null`
                 if [[ "$?" -eq 0 ]]; then
                     echo yes
                 else
                     echo no
                 fi
             fi
         else
             echo yes
         fi
     fi
          
     if [[ ! -z "$lx_mpi_compile_line" -a ! -z "$lx_mpi_link_line" ]]; then
         lx_mpi_command_line="$lx_mpi_compile_line $lx_mpi_link_line"
     fi

     if [[ ! -z "$lx_mpi_command_line" ]]; then
         # Now extract the different parts of the MPI command line.  Do these separately in case we need to 
         # parse them all out in future versions of this macro.
         lx_mpi_defines=`    echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-D\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_includes=`   echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-I\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_link_paths=` echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-L\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_libs=`       echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-l\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_link_args=`  echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-Wl,\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         
         # Create variables and clean up newlines and multiple spaces
         MPI_$3FLAGS="$lx_mpi_defines $lx_mpi_includes"
         MPI_$3LDFLAGS="$lx_mpi_link_paths $lx_mpi_libs $lx_mpi_link_args"
         MPI_$3FLAGS=`  echo "$MPI_$3FLAGS"   | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`
         MPI_$3LDFLAGS=`echo "$MPI_$3LDFLAGS" | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`

         OLD_CPPFLAGS=$CPPFLAGS
         OLD_LIBS=$LIBS
         CPPFLAGS=$MPI_$3FLAGS
         LIBS=$MPI_$3LDFLAGS

         AC_TRY_LINK([#include <mpi.h>],
                     [int rank, size;
                      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                      MPI_Comm_size(MPI_COMM_WORLD, &size);],
                     [# Add a define for testing at compile time.
                      AC_DEFINE([HAVE_MPI], [1], [Define to 1 if you have MPI libs and headers.])
                      have_$3_mpi='yes'],
                     [# zero out mpi flags so we don't link against the faulty library.
                      MPI_$3FLAGS=""
                      MPI_$3LDFLAGS=""
                      have_$3_mpi='no'])

         # AC_SUBST everything.
         AC_SUBST($1)
         AC_SUBST(MPI_$3FLAGS)
         AC_SUBST(MPI_$3LDFLAGS)

         LIBS=$OLD_LIBS
         CPPFLAGS=$OLD_CPPFLAGS
     else
         Echo Unable to find suitable MPI Compiler. Try setting $1.
         have_$3_mpi='no'         
     fi
])

AC_DEFUN([AX_GCC_OPTION], [
AC_REQUIRE([AC_PROG_CC])

AC_MSG_CHECKING([if gcc accepts $1 option])

AS_IF([ test "x$GCC" = "xyes" ],[
AS_IF([ test -z "$3" ],[
ax_gcc_option_test="int main()
{
return 0;
}"
],[
ax_gcc_option_test="$3"
])

# Dump the test program to file
cat <<EOF > conftest.c
$ax_gcc_option_test
EOF

# Dump back the file to the log, useful for debugging purposes
AC_TRY_COMMAND(cat conftest.c 1>&AS_MESSAGE_LOG_FD)

AS_IF([ AC_TRY_COMMAND($CC $2 $1 -c conftest.c 1>&AS_MESSAGE_LOG_FD) ],[
AC_MSG_RESULT([yes])
$4
],[
AC_MSG_RESULT([no])
$5
])
],[
AC_MSG_RESULT([no gcc available])
])
])

AC_DEFUN([AX_GCC_VERSION], [
GCC_VERSION=""
AX_GCC_OPTION([-dumpversion],[],[],[
ax_gcc_version_option=yes
],[
ax_gcc_version_option=no
])
AS_IF([test "x$GCC" = "xyes"],[
AS_IF([test "x$ax_gcc_version_option" != "xno"],[
AC_CACHE_CHECK([gcc version],[ax_cv_gcc_version],[
ax_cv_gcc_version="`$CC -dumpversion`"
AS_IF([test "x$ax_cv_gcc_version" = "x"],[
ax_cv_gcc_version=""
])
])
GCC_VERSION=$ax_cv_gcc_version
])
])
AC_SUBST([GCC_VERSION])
])

AC_DEFUN([AX_GXX_VERSION], [
GXX_VERSION=""
AX_GCC_OPTION([-dumpversion],[],[],[
ax_gcc_version_option=yes
],[
ax_gcc_version_option=no
])
AS_IF([test "x$GXX" = "xyes"],[
AS_IF([test "x$ax_gxx_version_option" != "no"],[
AC_CACHE_CHECK([gxx version],[ax_cv_gxx_version],[
ax_cv_gxx_version="`$CXX -dumpversion`"
AS_IF([test "x$ax_cv_gxx_version" = "x"],[
ax_cv_gxx_version=""
])
])
GXX_VERSION=$ax_cv_gxx_version
])
])
AC_SUBST([GXX_VERSION])
])

AC_DEFUN([AX_COMPILER_VENDOR],
[
AC_CACHE_CHECK([for _AC_LANG compiler vendor], ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor,
[ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor=unknown
# note: don't check for gcc first since some other compilers define __GNUC__
for ventest in intel:__ICC,__ECC,__INTEL_COMPILER ibm:__xlc__,__xlC__,__IBMC__,__IBMCPP__ pathscale:__PATHCC__,__PATHSCALE__ gnu:__GNUC__ sun:__SUNPRO_C,__SUNPRO_CC hp:__HP_cc,__HP_aCC dec:__DECC,__DECCXX,__DECC_VER,__DECCXX_VER borland:__BORLANDC__,__TURBOC__ comeau:__COMO__ cray:_CRAYC kai:__KCC lcc:__LCC__ metrowerks:__MWERKS__ sgi:__sgi,sgi microsoft:_MSC_VER watcom:__WATCOMC__ portland:__PGI; do
vencpp="defined("`echo $ventest | cut -d: -f2 | sed 's/,/) || defined(/g'`")"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
#if !($vencpp)
thisisanerror;
#endif
])], [ax_cv_]_AC_LANG_ABBREV[_compiler_vendor=`echo $ventest | cut -d: -f1`; break])
done
])
])

# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_11.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_11([ext|noext],[mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++11
#   standard; if necessary, add switches to CXXFLAGS to enable support.
#
#   The first argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for an extended mode.
#
#   The second argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline C++11 support is required and that the macro
#   should error out if no mode with that support is found.  If specified
#   'optional', then configuration proceeds regardless, after defining
#   HAVE_CXX11 if and only if a supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014 Alexey Sokolov <sokolov@google.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

m4_define([_AX_CXX_COMPILE_STDCXX_11_testbody], [[
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    struct Base {
    virtual void f() {}
    };
    struct Child : public Base {
    virtual void f() override {}
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);

    auto d = a;
    auto l = [](){};
]])

AC_DEFUN([AX_CXX_COMPILE_STDCXX_11], [dnl
  m4_if([$1], [], [],
        [$1], [ext], [],
        [$1], [noext], [],
        [m4_fatal([invalid argument `$1' to AX_CXX_COMPILE_STDCXX_11])])dnl
  m4_if([$2], [], [ax_cxx_compile_cxx11_required=true],
        [$2], [mandatory], [ax_cxx_compile_cxx11_required=true],
        [$2], [optional], [ax_cxx_compile_cxx11_required=false],
        [m4_fatal([invalid second argument `$2' to AX_CXX_COMPILE_STDCXX_11])])
  AC_LANG_PUSH([C++])dnl
  ac_success=no
  AC_CACHE_CHECK(whether $CXX supports C++11 features by default,
  ax_cv_cxx_compile_cxx11,
  [AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_11_testbody])],
    [ax_cv_cxx_compile_cxx11=yes],
    [ax_cv_cxx_compile_cxx11=no])])
  if test x$ax_cv_cxx_compile_cxx11 = xyes; then
    ac_success=yes
  fi

  m4_if([$1], [noext], [], [dnl
  if test x$ac_success = xno; then
    for switch in -std=gnu++11 -std=gnu++0x; do
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx11_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++11 features with $switch,
                     $cachevar,
        [ac_save_CXXFLAGS="$CXXFLAGS"
         CXXFLAGS="$CXXFLAGS $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_11_testbody])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXXFLAGS="$ac_save_CXXFLAGS"])
      if eval test x\$$cachevar = xyes; then
        CXXFLAGS="$CXXFLAGS $switch"
        ac_success=yes
        break
      fi
    done
  fi])

  m4_if([$1], [ext], [], [dnl
  if test x$ac_success = xno; then
    for switch in -std=c++11 -std=c++0x; do
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx11_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++11 features with $switch,
                     $cachevar,
        [ac_save_CXXFLAGS="$CXXFLAGS"
         CXXFLAGS="$CXXFLAGS $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_11_testbody])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXXFLAGS="$ac_save_CXXFLAGS"])
      if eval test x\$$cachevar = xyes; then
        CXXFLAGS="$CXXFLAGS $switch"
        ac_success=yes
        break
      fi
    done
  fi])
  AC_LANG_POP([C++])
  if test x$ax_cxx_compile_cxx11_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler with support for C++11 language features is required.])
    fi
  else
    if test x$ac_success = xno; then
      HAVE_CXX11=0
      AC_MSG_NOTICE([No compiler with C++11 support was found])
    else
      HAVE_CXX11=1
      AC_DEFINE(HAVE_CXX11,1,
                [define if the compiler supports basic C++11 syntax])
    fi

    AC_SUBST(HAVE_CXX11)
  fi
])

