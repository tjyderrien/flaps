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

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = external_libs  \
          src \
          doc

DISTCLEANFILES = \
        ChangeLog config.log config.status

doc:
	@cd doc; $(MAKE)

