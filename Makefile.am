bin_PROGRAMS = zon
zon_SOURCES = src/zon.c
EXTRA_DIST = src/astro.c
man_MANS = zon.1
zon.1: README zon
	-$(HELP2MAN) --no-info -i README -o zon.1 -n "zon [Options...] -- scriptable time output about sun rise and set in UTC" ./zon ; rm ./zon.1.gz ; gzip --best -n ./zon.1
