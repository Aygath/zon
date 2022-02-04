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
