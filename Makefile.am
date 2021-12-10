bin_PROGRAMS = zon

zon_SOURCES = src/zon.c
EXTRA_DIST = src/astro.c src/DEBIAN/control 
# EXTRA_DIST = src/astro.c src/doc/examples/boot-zon-op.service src/doc/examples/zon-op src/doc/changelog src/doc/copyright /src/DEBIAN/control

dist_doc_DATA = src/doc/copyright src/doc/changelog

exampledir = $(docdir)/examples
dist_example_DATA = src/doc/examples/boot-zon-op.service src/doc/examples/zon-op

man_MANS = zon.1
zon.1: README zon
	-$(HELP2MAN) --no-info -i README -o zon.1 -n "zon [Options...] -- scriptable time output about sun rise and set in UTC" ./zon 
