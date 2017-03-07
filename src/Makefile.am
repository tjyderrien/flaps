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

include $(srcdir)/../common-rules.make

bin_PROGRAMS = flaps

flaps_SOURCES =        \
          grid.f90     \
          input.f90    \
          material.f90 \
          maths.f90    \
          mesh.f90     \
          mie.f90      \
          output.f90   \
          poisson.f90  \
          source.f90   \
          solver.f90   \
          types.f90    \
          string.f90   \
          timer.f90    \
          utils.f90    \
          profiler.f90 \
          restart.f90  \
          main_explicit_Mie.f90

flaps_LDADD = $(all_LIBS)

CLEANFILES += flaps

depend:
	@makedepf90 (flaps_SOURCES) 

main_explicit_Mie.o: main_explicit_Mie.f90 input.o material.o mesh.o mie.o output.o source.o \
                     solver.o grid.o timer.o restart.o
material.o: material.f90 types.o maths.o source.o
mesh.o: mesh.f90 types.o output.o
mie.o: mie.f90 maths.o
output.o: output.f90 types.o
source.o: source.f90 maths.o types.o
solver.o: solver.f90 types.o profiler.o
grid.o: grid.f90 maths.o types.o
input.o: input.f90 types.o string.o utils.o
string.o: string.f90
timer.o: timer.f90
utils.o: utils.f90
profiler.o : profiler.f90 timer.o
restart.o : restart.f90 types.o
