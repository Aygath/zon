bin_PROGRAMS = zon
zon_SOURCES = src/zon.c
man_MANS = zon.1
zon.1: README zon
	-$(HELP2MAN) -i README -o zon.1 ./zon
