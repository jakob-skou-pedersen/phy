dnl --------------------------------
dnl Initialization macros.
dnl --------------------------------

AC_INIT([phy], [0.1.0], [jsp@binf.ku.dk]) 
AM_INIT_AUTOMAKE([-Wall -Werror foreign]) 
AM_CONFIG_HEADER([config.h])


dnl --------------------------------
dnl Output 
dnl --------------------------------

AC_CONFIG_FILES([Makefile phy/Makefile phy/test/Makefile src/Makefile]) 

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

GENERIC_LIBRARY_NAME=phy
PACKAGE=$GENERIC_LIBRARY_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)

#release versioning
GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=1
GENERIC_MICRO_VERSION=0

VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(VERSION)
AC_SUBST(RELEASE)

API_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(API_VERSION)

LIB_VERSION=$GENERIC_MAJOR_VERSION:$GENERIC_MINOR_VERSION:$GENERIC_MICRO_VERSION
AC_SUBST(LIB_VERSION)

#shared library versioning
#GENERIC_LIBRARY_VERSION=0:0:1
#                        | | |
#                 +------+ | +---+
#                 |        |     |
#              current:revision:age
#                 |        |     |
#                 |        |     +- increment if interfaces have been added
#                 |        |        set to zero if interfaces have been removed
#                                   or changed
#                 |        +- increment if source code has changed
#                 |           set to zero if current is incremented
#                 +- increment if interfaces have been added, removed or changed


dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

AC_PROG_CC
AC_PROG_CXX
AM_PROG_AR
LT_INIT 
AM_SANITY_CHECK

AC_LANG_CPLUSPLUS

dnl -----------------------------------------------
dnl Generates Makefile's, configuration files and scripts
dnl -----------------------------------------------

AC_OUTPUT
