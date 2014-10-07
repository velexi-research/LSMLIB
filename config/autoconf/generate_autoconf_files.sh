#!/bin/sh
##
## File:        generate_autoconf_files.sh
## Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
##                  Regents of the University of Texas.  All rights reserved.
##              (c) 2009 Kevin T. Chu.  All rights reserved.
## Revision:    $Revision$
## Modified:    $Date$
## Description: shell script for generating configure script using autoconf
##
## NOTES:
## - This script is intended to be called from the config/autoconf directory.
##

# Find top-level directory
TOP_DIR=`pwd`
CHECK_FILE=pylsmlib
CHECK=`ls $TOP_DIR | grep $CHECK_FILE`
while [ -z "$CHECK" ]; do
    TOP_DIR=`dirname $TOP_DIR`
    CHECK=`ls $TOP_DIR | grep $CHECK_FILE`
done
cd $TOP_DIR

# generate aclocal.m4
aclocal -Iconfig/autoconf/ --output=config/autoconf/aclocal.m4

# generate configure that searches config/autoconf for macro files
autoconf -Iconfig/autoconf configure.ac > $TOP_DIR/configure
chmod a+x configure

# remove autom4te.cache
rm -rf autom4te.cache
