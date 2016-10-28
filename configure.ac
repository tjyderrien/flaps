#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

## Copyright (C) 2016 N. Tancogne-Dejean
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>

AC_PREREQ([2.69])
AC_INIT([Flaps], [1.0], [derrien@fzu.cz])
AC_CONFIG_SRCDIR([src/maths.f90])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

#The computer
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_RANLIB

#Fortran stuff
AC_FC_LINE_LENGTH(132)
AC_PROG_F77

#Lets try to find the compiler
acx_save_fcflags="${FCFLAGS}"
AC_PROG_FC([], Fortran 90)
if test x"$FC" = x; then
  AC_MSG_ERROR([could not find Fortran 90 compiler])
fi
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
FCFLAGS="${acx_save_fcflags}"

AC_FC_LIBRARY_LDFLAGS

dnl check flag to include module path
AX_F90_MODULE_FLAG
F90_MODULE_FLAG=$ax_cv_f90_modflag
AC_SUBST([F90_MODULE_FLAG])

##OpenMP
#AC_ARG_ENABLE([openmp],
#[  --enable-openmp    Turn on OpenMP parallelization],
#[case "${enableval}" in
#  yes) openmp=true;;
#  no)  openmp=false;;
#  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
#esac],[openmp=false])
#if (test x$openmp = xtrue); then
#    AC_OPENMP
#    OMP_PREPRO="-DOMP"
#else
#    OMP_PREPRO=""
#fi
#AC_SUBST([OMP_PREPRO])

# Checks for libraries.

dnl check for BLAS
ACX_BLAS([], AC_MSG_ERROR([could not find required BLAS library]))


dnl check for LAPACK
ACX_LAPACK([], AC_MSG_ERROR([could not find required LAPACK library]))

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_PATH_PROG([PERL],perl,no)

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 external_libs/Makefile
                 external_libs/amos/Makefile
                 external_libs/gmsh/Makefile
                 src/Makefile])
AC_OUTPUT
