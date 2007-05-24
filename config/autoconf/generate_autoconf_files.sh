#!/bin/sh
##
## File:        generate_autoconf_files.m4
## Copyright:   (c) 2005-2006 Princeton University
## Authors:     Kevin T. Chu
## Revision:    $$
## Modified:    $$
## Description: shell script for generating configure script using autoconf
##
## NOTES:
## - This script is intended to be called from the config/autoconf directory.
##

# change directories to top level source directory
cd ../..

# generate aclocal.m4
aclocal --acdir=config/autoconf/ --output=config/autoconf/aclocal.m4 

# generate configure that searches config/autoconf for macro files
autoconf -Iconfig/autoconf configure.ac > configure
chmod a+x configure

# remove autom4te.cache
rm -rf autom4te.cache
