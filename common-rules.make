## Process this file with automake to produce Makefile.in

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

# ---------------------------------------------------------------
# Modules Path.
# ---------------------------------------------------------------

FCFLAGS_MODS =  \
      @F90_MODULE_FLAG@$(top_builddir)/external_libs      \
      @F90_MODULE_FLAG@$(top_builddir)/external_libs/amos \
      @F90_MODULE_FLAG@$(top_builddir)/external_libs/gmsh

# ---------------------------------------------------------------
# Define libraries here.
# ---------------------------------------------------------------

flaps_LIBS = \
       $(top_builddir)/external_libs/libflaps.a \
       $(top_builddir)/external_libs/amos/libamos.a \
       $(top_builddir)/external_libs/gmsh/libgmsh.a

core_LIBS = \
      @LIBS_LAPACK@ \
      @LIBS_BLAS@

all_LIBS = $(flaps_LIBS) $(core_LIBS)


# ---------------------------------------------------------------
# How to compile F90 files.
# ---------------------------------------------------------------
SUFFIXES = .f90 .f .o

.f90.o:
	@FC@ @FCFLAGS@ $(FCFLAGS_MODS) -c -o $@ $*.f90

.f.o:
	@FC@ @FCFLAGS@ $(FCFLAGS_MODS) -c -o $@ $*.f


# ---------------------------------------------------------------
# Miscellaneous.
# ---------------------------------------------------------------

# Cleaning.
CLEANFILES = *~ *.bak *.mod



