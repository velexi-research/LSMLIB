#!/bin/bash

#############################################################################
#
# This script builds the Google Test framework to support unit testing.
#
# Author: Kevin Chu
#
#############################################################################

# Set script name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Help function
function HELP {
    echo "Usage: ${SCRIPT} libdir tar_archive"
}

# Process command-line arguments
if [ $# -ne 2 ]; then
    echo "${SCRIPT}:ERROR: incorrect number of arguments"
    echo
    HELP
    echo
    exit 1
fi
LIBDIR=$1
TARBALL=$2

# Preparations
if [ -d "${LIBDIR}/gtest" ]; then
    rm -rf ${LIBDIR}/gtest
fi
mkdir -p ${LIBDIR}/gtest

# Build Google Test library
tar xfz ${TARBALL} -C ${LIBDIR}/gtest --strip-components=1
cd ${LIBDIR}/gtest
cmake . && make

# Return success code
exit 0
