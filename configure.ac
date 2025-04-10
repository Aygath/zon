dnl Template file for GNU Autoconf
dnl Copyright 2021-2022  Michael Welter

dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.

dnl You should have received a copy of the GNU General Public License
dnl along with this program.  If not, see <http://www.gnu.org/licenses/>.

AC_INIT([zon], [1.2.2~alfa], [https://github.com/Aygath/zon])
dnl Allow us to name files in subdirs, so that we do not need to have Makelfiles there: 
AM_INIT_AUTOMAKE([subdir-objects])
dnl Change default prefix for all files to "/" (resulting in /etc, /usr ...), (default for autotools is /usr/local): 
AC_PREFIX_DEFAULT (/)
dnl Test if we have a C compiler:
AC_PROG_CC

dnl Haiku does not use libm for the math functions, they are part
dnl of the C library
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the cos() function])
])
dnl Besides a C compiler, we also need the "help2man" program, in order to generate man-page:
AC_PATH_PROG(HELP2MAN, help2man, false // No help2man //)
AC_CONFIG_FILES([Makefile])

PKG_PROG_PKG_CONFIG
AC_ARG_WITH([systemdsystemunitdir],
     [AS_HELP_STRING([--with-systemdsystemunitdir=DIR], [Directory for systemd service files])],,
     [with_systemdsystemunitdir=auto])
AS_IF([test "x$with_systemdsystemunitdir" = "xyes" -o "x$with_systemdsystemunitdir" = "xauto"], [
     def_systemdsystemunitdir=$($PKG_CONFIG --variable=systemdsystemunitdir systemd)

     AS_IF([test "x$def_systemdsystemunitdir" = "x"],
   [AS_IF([test "x$with_systemdsystemunitdir" = "xyes"],
    [AC_MSG_ERROR([systemd support requested but pkg-config unable to query systemd package])])
    with_systemdsystemunitdir=no],
   [with_systemdsystemunitdir="$def_systemdsystemunitdir"])])
AS_IF([test "x$with_systemdsystemunitdir" != "xno"],
      [AC_SUBST([systemdsystemunitdir], [$with_systemdsystemunitdir])
       AC_DEFINE(HAVE_SYSTEMD,, Define if you want to use systemd socket activation) ])
AM_CONDITIONAL([HAVE_SYSTEMD], [test "x$with_systemdsystemunitdir" != "xno"])

AC_OUTPUT
