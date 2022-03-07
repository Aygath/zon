# Makefile for `zon' utility
# Copyright (C) 2021-2022  Michael Welter

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


bin_PROGRAMS = zon

AM_CFLAGS = -DSYSCONFDIR=\"$(sysconfdir)\"
zon_SOURCES = src/zon.c
EXTRA_DIST = src/astro.c src/zon_options.c zon.1 

exampledir = $(docdir)/examples
dist_example_DATA = src/doc/examples/boot-zon-op.service src/doc/examples/zon-op

man_MANS = zon.1
zon.1: src/zon_options.c README configure.ac
	@echo "W: !!!! WARNING !!!! Manpage zon.1 might be outdated - regenerate with 'make help2man' if you are NOT cross-compiling." 
help2man: 
	-$(HELP2MAN) --no-info -i README -o zon.1 -n "zon [Options...] -- scriptable time output about sun rise and set in UTC" ./zon 
