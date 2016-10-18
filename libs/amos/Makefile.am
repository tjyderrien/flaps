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

include $(srcdir)/../../common-rules.make

noinst_LIBRARIES = libamos.a

libamos_a_SOURCES =              \
                  d1mach.f  zabs.f   zasyi.f  zbknu.f  zexp.f   zmlt.f  \
                  zshch.f  zuni1.f  zunk2.f  \
                  dgamln.f  zacai.f  zbesh.f  zbuni.f  zkscl.f  zrati.f \
                  zsqrt.f  zuni2.f  zuoik.f \
                  i1mach.f  zacon.f  zbesj.f  zbunk.f  zlog.f   zs1s2.f \
                  zuchk.f  zunik.f  zwrsk.f \
                  xerror.f  zairy.f  zbinu.f  zdiv.f \
                  zmlri.f  zseri.f  zunhj.f  zunk1.f

depend:
	@makedepf90 $(libamos_a_SOURCES)
